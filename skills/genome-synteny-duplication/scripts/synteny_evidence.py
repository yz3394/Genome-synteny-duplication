#!/usr/bin/env python3
"""Generic, family-independent MCScanX and duplication evidence.

Coordinates are 1-based inclusive. Tandem uses annotation_rank difference 1;
block membership uses engine_rank between the first and last side anchors.
All homology coverage is the span of one best HSP, never alignment length.
The raw m8 remains the authoritative directional and multiple-HSP evidence.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path

BLOCK_FIELDS = 'comparison_id block_id genome_a seqid_a start_a end_a genome_b seqid_b start_b end_b orientation anchor_count score block_evalue'.split()
ANCHOR_FIELDS = 'comparison_id block_id genome_a gene_a engine_a seqid_a start_a end_a genome_b gene_b engine_b seqid_b start_b end_b hit_evalue'.split()
UNIQUE_FIELDS = 'comparison_id genome_a gene_a engine_a seqid_a start_a end_a genome_b gene_b engine_b seqid_b start_b end_b supporting_block_ids supporting_block_count best_anchor_evalue'.split()
MEMBER_FIELDS = 'comparison_id block_id side genome_id gene_id engine_id seqid start end engine_rank is_anchor anchor_partner_ids in_rank_interval overlaps_physical_span'.split()
HOMOLOGY_FIELDS = ('comparison_id genome_a gene_a protein_a engine_a seqid_a start_a end_a strand_a annotation_rank_a '
                   'genome_b gene_b protein_b engine_b seqid_b start_b end_b strand_b annotation_rank_b '
                   'pident alignment_length evalue bitscore best_hsp_query best_hsp_subject qstart qend sstart send '
                   'span_aa_a span_aa_b protein_length_a protein_length_b coverage_a coverage_b coverage_method observed_directions review_flags').split()
PAIR_FIELDS = HOMOLOGY_FIELDS + 'genome_id seqid relation_type intervening_genes rank_difference bp_gap genes_overlap quality_tier array_id'.split()
ARRAY_FIELDS = 'comparison_id array_id genome_id seqid start end gene_ids engine_ids gene_count edge_count relation_type'.split()
CLASS_FIELDS = 'comparison_id genome_id gene_id engine_id has_protein native_code native_label assessment_status has_independent_tandem has_collinear_anchor'.split()
HEADER_RE = re.compile(r'^##\s*Alignment\s+(\d+)\s*:\s*score=(\S+)\s+e_value=(\S+)\s+N=(\d+)\s+(\S+)&(\S+)\s+(plus|minus)\s*$')
ANCHOR_RE = re.compile(r'^\s*(\d+)\s*-\s*(\d+)\s*:\s*(\S+)\s+(\S+)\s+(\S+)\s*$')
NATIVE_LABELS = {0: 'singleton', 1: 'dispersed', 2: 'proximal', 3: 'tandem', 4: 'WGD_or_segmental'}


def _true(value):
    return value is True or str(value).lower() in {'1', 'true', 'yes', 'y'}


def _read(path):
    with Path(path).open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def _write(path, fields, rows):
    count = 0
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def _genes(genes):
    by_engine, by_locus, ranks = {}, {}, defaultdict(set)
    for original in genes:
        g = dict(original)
        for field in ('start', 'end', 'annotation_rank', 'engine_rank'):
            g[field] = int(g[field])
        g['protein_length'] = int(g.get('protein_length') or 0)
        g['has_protein'] = _true(g['has_protein'])
        if not g.get('engine_id') or g['start'] < 1 or g['end'] < g['start']:
            raise ValueError(f'Invalid gene coordinates/engine ID: {g}')
        if min(g['annotation_rank'], g['engine_rank']) < 1:
            raise ValueError(f'Ranks must be positive: {g["engine_id"]}')
        key = (g['genome_id'], g['gene_id'])
        if key in by_locus or g['engine_id'] in by_engine:
            raise ValueError(f'Duplicate locus or engine ID: {key}')
        by_engine[g['engine_id']] = by_locus[key] = g
        for rank in ('annotation_rank', 'engine_rank'):
            rk = (g['genome_id'], g['seqid'], rank)
            if g[rank] in ranks[rk]:
                raise ValueError(f'Duplicate {rank} on {rk[:2]}')
            ranks[rk].add(g[rank])
    for key, values in ranks.items():
        if values != set(range(1, len(values) + 1)):
            raise ValueError(f'Nonconsecutive full-annotation ranks: {key}')
    return by_engine, by_locus


def _order(a, b, comparison):
    """Canonical pair key; cross pairs follow configured a -> b, not file order."""
    ga, gb = comparison['a'], comparison['b']
    if ga == gb:
        if a['genome_id'] != ga or b['genome_id'] != ga:
            raise ValueError(f'Unexpected genome in self comparison {comparison["id"]}')
        return tuple(sorted((a, b), key=lambda g: (g['seqid'], g['annotation_rank'], g['gene_id'])))
    if a['genome_id'] == gb and b['genome_id'] == ga:
        return b, a
    if a['genome_id'] != ga or b['genome_id'] != gb:
        raise ValueError(f'Unexpected genome pair in comparison {comparison["id"]}')
    return a, b


def _side_fields(a, b):
    row = {}
    for letter, gene in (('a', a), ('b', b)):
        for field, out in [('genome_id','genome'), ('gene_id','gene'), ('protein_id','protein'), ('engine_id','engine'),
                           ('seqid','seqid'), ('start','start'), ('end','end'), ('strand','strand'),
                           ('annotation_rank','annotation_rank'), ('protein_length','protein_length')]:
            row[f'{out}_{letter}'] = gene.get(field, '')
    return row


def _blocks(path, by_engine, comparison):
    if path is None:
        return []
    blocks, current, ids = [], None, set()
    with Path(path).open() as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith('## Alignment'):
                match = HEADER_RE.fullmatch(line)
                if match is None:
                    raise ValueError(f'{path}:{lineno}: malformed MCScanX block header')
                bid, score, evalue, count, seq_a, seq_b, orientation = match.groups()
                if int(count) < 1 or not all(math.isfinite(float(v)) for v in (score,evalue)) or float(evalue) < 0:
                    raise ValueError(f'{path}:{lineno}: invalid block count/score/E value')
                if bid in ids:
                    raise ValueError(f'{path}:{lineno}: duplicate block {bid}')
                ids.add(bid)
                current = dict(block_id=f'{comparison["id"]}_B{int(bid):06d}', score=score, block_evalue=evalue,
                               expected=int(count), seq_a=seq_a, seq_b=seq_b, orientation=orientation, anchors=[])
                blocks.append(current)
            elif line.startswith('#'):
                continue
            else:
                match = ANCHOR_RE.fullmatch(line)
                if match is None or current is None:
                    raise ValueError(f'{path}:{lineno}: malformed/orphan anchor')
                raw_bid, _, aid, bid, evalue = match.groups()
                if current['block_id'] != f'{comparison["id"]}_B{int(raw_bid):06d}':
                    raise ValueError(f'{path}:{lineno}: anchor block ID mismatch')
                if aid not in by_engine or bid not in by_engine:
                    raise ValueError(f'{path}:{lineno}: anchor IDs absent from normalized gene map')
                a, b = by_engine[aid], by_engine[bid]
                _order(a, b, comparison)  # validates comparison membership
                if a['engine_seqid'] != current['seq_a'] or b['engine_seqid'] != current['seq_b']:
                    raise ValueError(f'{path}:{lineno}: anchor/header sequence mismatch')
                if aid == bid or not a['has_protein'] or not b['has_protein']:
                    raise ValueError(f'{path}:{lineno}: self or proteinless anchor')
                if not math.isfinite(float(evalue)) or float(evalue) < 0:
                    raise ValueError(f'{path}:{lineno}: invalid anchor E value')
                if comparison['a'] != comparison['b']:
                    a, b = _order(a, b, comparison)
                current['anchors'].append((a, b, evalue))
    for block in blocks:
        if len(block['anchors']) != block['expected'] or not block['anchors']:
            raise ValueError(f'{path}: anchor count mismatch/empty block {block["block_id"]}')
        keys = [(a['engine_id'], b['engine_id']) for a, b, _ in block['anchors']]
        if len(keys) != len(set(keys)):
            raise ValueError(f'{path}: repeated pair inside block {block["block_id"]}')
    return blocks


def _homology(raw_hits, by_engine, comparison, out_dir, max_evalue, max_intervening):
    """Disk-backed best unique pairs: preserves large search pools without huge RAM."""
    near = []
    stats = dict(raw_hsp_rows=0, passing_nonself_hsp_rows=0, self_hsp_rows=0, unique_homology_pairs=0)
    with tempfile.TemporaryDirectory(prefix='.homology-', dir=out_dir) as temp:
        conn = sqlite3.connect(str(Path(temp) / 'pairs.sqlite'))
        try:
            conn.execute('PRAGMA journal_mode=OFF')
            conn.execute('PRAGMA synchronous=OFF')
            conn.execute('CREATE TABLE pairs (a TEXT,b TEXT,bits REAL,ev REAL,alen INTEGER,directions INTEGER,data TEXT,PRIMARY KEY(a,b)) WITHOUT ROWID')
            sql = '''INSERT INTO pairs VALUES (?,?,?,?,?,?,?) ON CONFLICT(a,b) DO UPDATE SET
                directions=(pairs.directions | excluded.directions),
                data=CASE WHEN excluded.bits>pairs.bits OR (excluded.bits=pairs.bits AND excluded.ev<pairs.ev)
                    OR (excluded.bits=pairs.bits AND excluded.ev=pairs.ev AND excluded.alen>pairs.alen) THEN excluded.data ELSE pairs.data END,
                alen=CASE WHEN excluded.bits>pairs.bits OR (excluded.bits=pairs.bits AND excluded.ev<pairs.ev)
                    OR (excluded.bits=pairs.bits AND excluded.ev=pairs.ev AND excluded.alen>pairs.alen) THEN excluded.alen ELSE pairs.alen END,
                ev=CASE WHEN excluded.bits>pairs.bits OR (excluded.bits=pairs.bits AND excluded.ev<pairs.ev) THEN excluded.ev ELSE pairs.ev END,
                bits=MAX(pairs.bits,excluded.bits)'''
            for path in raw_hits:
                with Path(path).open() as handle:
                    for lineno, line in enumerate(handle, 1):
                        if not line.strip() or line.startswith('#'):
                            continue
                        fields = line.split()
                        if len(fields) != 12:
                            raise ValueError(f'{path}:{lineno}: expected exactly 12 m8 columns')
                        qid, sid = fields[:2]
                        if qid not in by_engine or sid not in by_engine:
                            raise ValueError(f'{path}:{lineno}: homology ID absent from normalized gene map')
                        q, s = by_engine[qid], by_engine[sid]
                        a, b = _order(q, s, comparison)
                        if not q['has_protein'] or not s['has_protein']:
                            raise ValueError(f'{path}:{lineno}: proteinless homology hit')
                        stats['raw_hsp_rows'] += 1
                        if qid == sid:
                            stats['self_hsp_rows'] += 1
                            continue
                        pident, evalue, bits = float(fields[2]), float(fields[10]), float(fields[11])
                        alen = int(fields[3])
                        qstart,qend,sstart,send = map(int, fields[6:10])
                        if not all(math.isfinite(v) for v in (pident,evalue,bits)) or not 0 <= pident <= 100 or min(evalue,bits) < 0 or alen < 1:
                            raise ValueError(f'{path}:{lineno}: invalid numeric homology fields')
                        if min(qstart,qend,sstart,send) < 1 or max(qstart,qend) > q['protein_length'] or max(sstart,send) > s['protein_length']:
                            raise ValueError(f'{path}:{lineno}: HSP coordinates outside protein sequence')
                        if evalue > max_evalue:
                            continue
                        stats['passing_nonself_hsp_rows'] += 1
                        direction = 1 if a['engine_id'] == qid else 2
                        payload = [pident,alen,fields[10],bits,qid,sid,qstart,qend,sstart,send]
                        conn.execute(sql, (a['engine_id'],b['engine_id'],bits,evalue,alen,direction,json.dumps(payload,separators=(',',':'))))
                conn.commit()
            def rows():
                for aid,bid,directions,payload in conn.execute('SELECT a,b,directions,data FROM pairs ORDER BY a,b'):
                    a,b = by_engine[aid],by_engine[bid]
                    pident,alen,evalue,bits,qid,sid,qs,qe,ss,se = json.loads(payload)
                    qspan,sspan = abs(qe-qs)+1,abs(se-ss)+1
                    span_a,span_b = (qspan,sspan) if aid == qid else (sspan,qspan)
                    ca,cb = span_a/a['protein_length'],span_b/b['protein_length']
                    flags = []
                    if ca < .5 or cb < .5:
                        flags.append('single_hsp_coverage_below_0.5')
                    for letter,gene in (('a',a),('b',b)):
                        if gene.get('model_flags'):
                            flags.append(f'model_{letter}:{gene["model_flags"]}')
                    row = dict(comparison_id=comparison['id'], **_side_fields(a,b),pident=pident,alignment_length=alen,
                               evalue=evalue,bitscore=bits,best_hsp_query=qid,best_hsp_subject=sid,qstart=qs,qend=qe,sstart=ss,send=se,
                               span_aa_a=span_a,span_aa_b=span_b,coverage_a=f'{ca:.8g}',coverage_b=f'{cb:.8g}',
                               coverage_method='single_best_HSP_span',observed_directions=int(bool(directions&1))+int(bool(directions&2)),review_flags=';'.join(flags))
                    if comparison['a'] == comparison['b'] and a['seqid'] == b['seqid']:
                        diff = abs(a['annotation_rank']-b['annotation_rank'])
                        if 1 <= diff <= max_intervening+1:
                            overlap = max(a['start'],b['start']) <= min(a['end'],b['end'])
                            if overlap:
                                flags.append('overlapping_gene_models')
                            candidate = dict(row,genome_id=a['genome_id'],seqid=a['seqid'],relation_type='strict_tandem' if diff == 1 else 'proximal',
                                             intervening_genes=diff-1,rank_difference=diff,bp_gap=b['start']-a['end']-1,
                                             genes_overlap=int(overlap),quality_tier='review' if flags else 'candidate_supported_by_homology',array_id='',review_flags=';'.join(flags))
                            near.append(candidate)
                    yield row
            stats['unique_homology_pairs'] = _write(out_dir/'homology_pairs.tsv',HOMOLOGY_FIELDS,rows())
        finally:
            conn.close()
    return near,stats


def process_comparison(genes, raw_hits, collinearity, out_dir, comparison, search_evalue=1e-5,
                       proximal_max_intervening=8, classifier_path=None):
    """Write normalized evidence tables. No family assignments influence inference.

    classifier_path is a native MCScanX *.gene_type DATA file, not an executable.
    The caller controls output versioning and verifies hashes before cache reuse.
    Return counts, explicit stage statuses, rank definitions, and output paths.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True,exist_ok=True)
    if not math.isfinite(search_evalue) or search_evalue < 0 or proximal_max_intervening < 0:
        raise ValueError('Invalid evidence E value or proximal threshold')
    by_engine,by_locus = _genes(genes)
    for name in ('id','a','b'):
        if not comparison.get(name):
            raise ValueError(f'Missing comparison {name}')
    self_analysis = comparison['a'] == comparison['b']
    selected = [g for g in by_engine.values() if g['genome_id'] in {comparison['a'],comparison['b']}]
    if {g['genome_id'] for g in selected} != {comparison['a'],comparison['b']}:
        raise ValueError('Comparison references a genome absent from normalized genes')
    # An empty search file cannot support a biological negative when a required
    # genome has no searchable representatives at all.
    no_proteins = any(not any(g['has_protein'] for g in selected if g['genome_id']==genome)
                      for genome in {comparison['a'],comparison['b']})
    blocks = _blocks(collinearity,by_engine,comparison)
    by_seq = defaultdict(list)
    for gene in selected:
        by_seq[(gene['genome_id'],gene['seqid'])].append(gene)
    for rows in by_seq.values():
        rows.sort(key=lambda g:g['engine_rank'])
    block_rows,anchor_rows,unique,side_windows = [],[],{},[]
    for block in blocks:
        sides = []
        for index,side in ((0,'A'),(1,'B')):
            partners = defaultdict(set)
            for pair in block['anchors']:
                partners[pair[index]['engine_id']].add(pair[1-index]['gene_id'])
            side_genes = [pair[index] for pair in block['anchors']]
            sequence_keys = {(g['genome_id'],g['seqid']) for g in side_genes}
            if len(sequence_keys) != 1:
                raise ValueError(f'Block {block["block_id"]} side spans multiple sequences')
            genome,seqid = next(iter(sequence_keys))
            window = dict(block_id=block['block_id'],side=side,genome_id=genome,seqid=seqid,
                          lo=min(g['engine_rank'] for g in side_genes),hi=max(g['engine_rank'] for g in side_genes),
                          start=min(g['start'] for g in side_genes),end=max(g['end'] for g in side_genes),partners=partners)
            side_windows.append(window)
            sides.append(window)
        a,b = sides
        block_rows.append(dict(comparison_id=comparison['id'],block_id=block['block_id'],genome_a=a['genome_id'],seqid_a=a['seqid'],start_a=a['start'],end_a=a['end'],
                               genome_b=b['genome_id'],seqid_b=b['seqid'],start_b=b['start'],end_b=b['end'],orientation=block['orientation'],anchor_count=len(block['anchors']),
                               score=block['score'],block_evalue=block['block_evalue']))
        for a,b,evalue in block['anchors']:
            anchor_rows.append(dict(comparison_id=comparison['id'],block_id=block['block_id'],**_side_fields(a,b),hit_evalue=evalue))
            ca,cb = _order(a,b,comparison)
            key = (ca['engine_id'],cb['engine_id'])
            if key not in unique:
                unique[key] = dict(comparison_id=comparison['id'],**_side_fields(ca,cb),blocks=set(),best_anchor_evalue=evalue)
            unique[key]['blocks'].add(block['block_id'])
            if float(evalue) < float(unique[key]['best_anchor_evalue']):
                unique[key]['best_anchor_evalue'] = evalue
    def member_rows(outside_only=False):
        for window in side_windows:
            seq_genes = by_seq[(window['genome_id'],window['seqid'])]
            # Every locus, including proteinless loci, has an engine rank.
            candidates = seq_genes if outside_only else seq_genes[window['lo']-1:window['hi']]
            for gene in candidates:
                in_rank = window['lo'] <= gene['engine_rank'] <= window['hi']
                physical = gene['start'] <= window['end'] and gene['end'] >= window['start']
                if outside_only and (in_rank or not physical):
                    continue
                yield dict(comparison_id=comparison['id'],block_id=window['block_id'],side=window['side'],genome_id=gene['genome_id'],gene_id=gene['gene_id'],
                           engine_id=gene['engine_id'],seqid=gene['seqid'],start=gene['start'],end=gene['end'],engine_rank=gene['engine_rank'],
                           is_anchor=int(gene['engine_id'] in window['partners']),anchor_partner_ids=';'.join(sorted(window['partners'].get(gene['engine_id'],()))),
                           in_rank_interval=int(in_rank),overlaps_physical_span=int(physical))
    counts = {}
    counts['blocks'] = _write(out_dir/'blocks.tsv',BLOCK_FIELDS,block_rows)
    counts['anchor_rows'] = _write(out_dir/'anchor_pairs.tsv',ANCHOR_FIELDS,anchor_rows)
    unique_rows = [dict(row,supporting_block_ids=';'.join(sorted(row['blocks'])),supporting_block_count=len(row['blocks'])) for _,row in sorted(unique.items())]
    counts['unique_anchor_pairs'] = _write(out_dir/'unique_anchor_pairs.tsv',UNIQUE_FIELDS,unique_rows)
    counts['block_gene_rows'] = _write(out_dir/'all_block_genes.tsv',MEMBER_FIELDS,member_rows())
    counts['outside_rank_physical_overlap_rows'] = _write(out_dir/'block_physical_overlaps.tsv',MEMBER_FIELDS,member_rows(True))
    near,hit_stats = _homology(raw_hits,by_engine,comparison,out_dir,search_evalue,proximal_max_intervening)
    counts.update(hit_stats)
    tandem = sorted((r for r in near if r['intervening_genes']==0),key=lambda r:(r['seqid'],r['annotation_rank_a'],r['gene_a']))
    proximal = [r for r in near if r['intervening_genes']>0]
    # Components of strict adjacent edges are arrays; no proximity-only edge joins them.
    adjacency = defaultdict(set)
    for row in tandem:
        adjacency[row['engine_a']].add(row['engine_b'])
        adjacency[row['engine_b']].add(row['engine_a'])
    arrays,seen,array_for = [],set(),{}
    for origin in sorted(adjacency,key=lambda eid:(by_engine[eid]['seqid'],by_engine[eid]['annotation_rank'])):
        if origin in seen:
            continue
        stack,component = [origin],set()
        while stack:
            eid = stack.pop()
            if eid in component:
                continue
            component.add(eid)
            stack.extend(adjacency[eid]-component)
        seen.update(component)
        group = sorted((by_engine[eid] for eid in component),key=lambda g:g['annotation_rank'])
        array_id = f'{comparison["id"]}_TA{len(arrays)+1:06d}'
        array_for.update({eid:array_id for eid in component})
        arrays.append(dict(comparison_id=comparison['id'],array_id=array_id,genome_id=group[0]['genome_id'],seqid=group[0]['seqid'],start=min(g['start'] for g in group),end=max(g['end'] for g in group),
                           gene_ids=';'.join(g['gene_id'] for g in group),engine_ids=';'.join(g['engine_id'] for g in group),gene_count=len(group),edge_count=sum(len(adjacency[eid]) for eid in component)//2,relation_type='strict_tandem'))
    for row in tandem:
        row['array_id'] = array_for[row['engine_a']]
    counts['tandem_pairs'] = _write(out_dir/'tandem_pairs.tsv',PAIR_FIELDS,tandem)
    counts['proximal_pairs'] = _write(out_dir/'proximal_pairs.tsv',PAIR_FIELDS,proximal)
    counts['tandem_arrays'] = _write(out_dir/'tandem_arrays.tsv',ARRAY_FIELDS,arrays)
    native = {}
    if classifier_path is not None:
        if not self_analysis:
            raise ValueError('Native duplication classifier is only accepted for self comparison')
        with Path(classifier_path).open() as handle:
            for lineno,line in enumerate(handle,1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields)!=2 or fields[0] not in by_engine or fields[0] in native:
                    raise ValueError(f'{classifier_path}:{lineno}: invalid/duplicate classifier gene')
                code = int(fields[1])
                if code not in NATIVE_LABELS or by_engine[fields[0]]['genome_id'] != comparison['a']:
                    raise ValueError(f'{classifier_path}:{lineno}: invalid classifier code/genome')
                native[fields[0]] = code
    missing_native = [g['engine_id'] for g in selected if g['has_protein'] and g['engine_id'] not in native] if classifier_path is not None else []
    counts['missing_native_protein_loci'] = len(missing_native)
    anchored = {eid for pair in unique for eid in pair}
    def classification():
        for g in sorted(selected,key=lambda x:(x['genome_id'],x['seqid'],x['annotation_rank'])):
            eid = g['engine_id']
            status = ('not_applicable_cross_comparison' if not self_analysis else
                      'not_assessable_no_protein' if not g['has_protein'] else
                      'not_run' if classifier_path is None else 'classified' if eid in native else 'missing_native_row')
            code = native.get(eid) if status=='classified' else None
            yield dict(comparison_id=comparison['id'],genome_id=g['genome_id'],gene_id=g['gene_id'],engine_id=eid,has_protein=int(g['has_protein']),
                       native_code='' if code is None else code,native_label='' if code is None else NATIVE_LABELS[code],assessment_status=status,
                       has_independent_tandem=int(eid in adjacency),has_collinear_anchor=int(eid in anchored))
    counts['classification_gene_rows'] = _write(out_dir/'gene_duplication_classification.tsv',CLASS_FIELDS,classification())
    summary = dict(comparison=dict(comparison),counts=counts,statuses=dict(collinearity='not_assessable_no_proteins' if no_proteins else 'not_run' if collinearity is None else 'complete',
                   homology='not_assessable_no_proteins' if no_proteins else 'complete' if raw_hits else 'not_run',tandem='not_applicable_cross_comparison' if not self_analysis else 'not_assessable_no_proteins' if no_proteins else 'complete' if raw_hits else 'not_run',
                   proximal='not_applicable_cross_comparison' if not self_analysis else 'not_assessable_no_proteins' if no_proteins else 'complete' if raw_hits else 'not_run',
                   classifier='not_applicable_cross_comparison' if not self_analysis else 'not_assessable_no_proteins' if no_proteins else 'not_run' if classifier_path is None else 'incomplete_missing_native_rows' if missing_native else 'complete'),
                   definitions=dict(tandem='same genome and seqid; annotation_rank difference = 1; nonself protein homology E <= search_evalue',
                                    proximal=f'same genome and seqid; 1 <= intervening annotation loci <= {proximal_max_intervening}; homology E <= search_evalue',
                                    block_membership='inclusive engine_rank interval between extreme anchors, including loci without protein; A/B canonicalized for cross comparisons',
                                    physical_overlap='separate supplementary rows outside rank interval in block_physical_overlaps.tsv',
                                    coverage='span of best single HSP divided by full protein length per side; <0.5 is a review flag, not exclusion',
                                    duplicate_class='native mutually exclusive classification does not replace independent tandem or anchor evidence'),
                   parameters=dict(search_evalue=search_evalue,proximal_max_intervening=proximal_max_intervening),
                   outputs={p.stem:str(p.resolve()) for p in sorted(out_dir.glob('*.tsv'))})
    (out_dir/'evidence_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    return summary


def overlay_families(result_dirs, genes, families, out_dir, isoforms=None):
    """Map arbitrary family labels to loci; never modify or recompute core evidence.

    protein_id may name a nonrepresentative isoform when isoforms provides its
    genome_id/gene_id/protein_id mapping. Both gene and protein, when supplied,
    must resolve to the same locus. Original IDs remain untouched.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True,exist_ok=True)
    by_engine,by_locus = _genes(genes)
    proteins = {}
    for row in list(genes)+list(isoforms or []):
        if not row.get('protein_id'):
            continue
        key = (row['genome_id'],row['protein_id'])
        locus = (row['genome_id'],row['gene_id'])
        if locus not in by_locus:
            raise ValueError(f'Isoform mapping references unknown locus {locus}')
        if key in proteins and proteins[key] != locus:
            raise ValueError(f'Ambiguous protein-to-locus mapping {key}')
        proteins[key] = locus
    labels = defaultdict(set)
    source_proteins = defaultdict(set)
    for row in families:
        genome = row.get('genome_id')
        gene_key = (genome,row['gene_id']) if row.get('gene_id') else None
        protein_key = (genome,row['protein_id']) if row.get('protein_id') else None
        if gene_key is not None and gene_key not in by_locus:
            raise ValueError(f'Family gene ID is absent: {gene_key}')
        if protein_key is not None and protein_key not in proteins:
            raise ValueError(f'Family protein ID is absent; provide explicit isoform mapping: {protein_key}')
        from_protein = proteins.get(protein_key)
        if gene_key and from_protein and gene_key != from_protein:
            raise ValueError(f'Family gene/protein ID disagreement: {row}')
        locus = gene_key or from_protein
        if locus is None:
            raise ValueError(f'Family row needs genome_id and gene_id or protein_id: {row}')
        family_labels = {x.strip() for x in str(row.get('family_label') or row.get('families') or '').split(';') if x.strip()}
        if not family_labels:
            raise ValueError(f'Empty family label: {row}')
        labels[locus].update(family_labels)
        if protein_key:
            source_proteins[locus].add(protein_key[1])
    membership = []
    for locus,names in sorted(labels.items()):
        g = by_locus[locus]
        membership.append(dict(genome_id=locus[0],gene_id=locus[1],protein_id=g.get('protein_id',''),engine_id=g['engine_id'],families=';'.join(sorted(names)),
                               family_source_protein_ids=';'.join(sorted(source_proteins[locus])),has_protein=int(g['has_protein']),seqid=g['seqid'],start=g['start'],end=g['end'],strand=g['strand']))
    _write(out_dir/'family_members.tsv','genome_id gene_id protein_id engine_id families family_source_protein_ids has_protein seqid start end strand'.split(),membership)
    counts,seen_comparisons = {},set()
    resolved_dirs = [str(Path(p).resolve()) for p in result_dirs]
    if len(resolved_dirs) != len(set(resolved_dirs)):
        raise ValueError('Duplicate result directories passed to overlay')
    pair_files = {'anchor_pairs.tsv':'family_anchor_pairs.tsv','unique_anchor_pairs.tsv':'family_unique_anchor_pairs.tsv',
                  'tandem_pairs.tsv':'family_tandem_pairs.tsv','proximal_pairs.tsv':'family_proximal_pairs.tsv','homology_pairs.tsv':'family_homology_pairs.tsv'}
    for source_name,dest_name in pair_files.items():
        output_rows = []
        fields = None
        for result_dir in result_dirs:
            path = Path(result_dir)/source_name
            if not path.is_file():
                continue
            with path.open(newline='') as handle:
                reader = csv.DictReader(handle,delimiter='\t')
                expected_fields = list(reader.fieldnames or [])
                if fields is not None and fields != expected_fields:
                    raise ValueError(f'Inconsistent result schemas: {path}')
                fields = expected_fields
                for row in reader:
                    la = labels.get((row['genome_a'],row['gene_a']),set())
                    lb = labels.get((row['genome_b'],row['gene_b']),set())
                    if not la and not lb:
                        continue
                    output_rows.append(dict(row,families_a=';'.join(sorted(la)),families_b=';'.join(sorted(lb)),shared_families=';'.join(sorted(la&lb)),
                                            both_in_focal_families=int(bool(la) and bool(lb)),same_family_pair=int(bool(la&lb))))
        fallback = {'anchor_pairs.tsv':ANCHOR_FIELDS,'unique_anchor_pairs.tsv':UNIQUE_FIELDS,'tandem_pairs.tsv':PAIR_FIELDS,'proximal_pairs.tsv':PAIR_FIELDS,'homology_pairs.tsv':HOMOLOGY_FIELDS}[source_name]
        counts[dest_name] = _write(out_dir/dest_name,(fields or fallback)+'families_a families_b shared_families both_in_focal_families same_family_pair'.split(),output_rows)
    members = []
    for result_dir in result_dirs:
        result_dir = Path(result_dir)
        summary = result_dir/'evidence_summary.json'
        if summary.is_file():
            comparison_id = json.loads(summary.read_text())['comparison']['id']
            if comparison_id in seen_comparisons:
                raise ValueError(f'Duplicate comparison results passed to overlay: {comparison_id}')
            seen_comparisons.add(comparison_id)
        path = result_dir/'all_block_genes.tsv'
        if not path.is_file():
            continue
        for row in _read(path):
            names = labels.get((row['genome_id'],row['gene_id']))
            if names:
                members.append(dict(row,families=';'.join(sorted(names))))
    counts['family_block_members.tsv'] = _write(out_dir/'family_block_members.tsv',MEMBER_FIELDS+['families'],members)
    report = dict(family_loci=len(labels),family_label_assignments=sum(map(len,labels.values())),family_labels=sorted({x for names in labels.values() for x in names}),counts=counts,
                  core_inference_modified=False,definitions=dict(family_pairs='at least one focal endpoint; both_in_focal_families and same_family_pair distinguish bilateral and shared-label evidence',
                     block_members='rank-interval membership only; is_anchor distinguishes direct anchors; no anchor does not imply absence of homology'),
                  outputs={p.stem:str(p.resolve()) for p in sorted(out_dir.glob('*.tsv'))})
    (out_dir/'family_overlay_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    return report
