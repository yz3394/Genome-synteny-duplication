#!/usr/bin/env python3
"""Prepare immutable, family-independent loci and proteins for MCScanX.

Coordinates remain 1-based inclusive. IDs are never inferred by trimming a
suffix. GFF3 Parent chains, explicit protein_id attributes, exact FASTA IDs, or
an explicit mapping table establish identity. One FASTA sequence is selected
per gene; separate loci with identical proteins are retained.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

SCHEMA_VERSION = "1.0"
MCSCANX_MAX_COORDINATE = 2_147_483_647  # Standard MCScanX stores positions as signed int.
GENE_FIELDS = [
    "genome_id", "species", "assembly_id", "annotation_id", "subgenome",
    "seqid", "engine_seqid", "gene_id", "protein_id", "engine_id", "start",
    "end", "strand", "annotation_rank", "engine_rank", "has_protein",
    "protein_length", "sequence_sha256", "representative_reason", "model_flags",
]
ISOFORM_FIELDS = [
    "genome_id", "gene_id", "transcript_id", "protein_id", "engine_id",
    "is_representative", "protein_length", "sequence_sha256", "mapping_source",
    "model_flags",
]


class PreparationError(ValueError):
    """Inputs cannot be interpreted without an explicit correction or mapping."""


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _identifier(value, label: str) -> str:
    if not isinstance(value, str) or not value or any(c in value for c in "\t\r\n"):
        raise PreparationError(f"Missing or invalid {label}: {value!r}")
    return value


def _table(path: Path, required: set[str]) -> list[dict]:
    with _open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise PreparationError(f"{path}: required columns {sorted(required)}")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise PreparationError(f"Malformed TSV row: {path}")
    return rows


def _fasta(path: Path, proteins: bool = True) -> tuple[dict[str, str], dict[str, list[str]]]:
    sequences, flags = {}, {}
    identifier, chunks = None, []

    def finish():
        if identifier is None:
            return
        if identifier in sequences:
            raise PreparationError(f"Duplicate FASTA ID {identifier!r}: {path}")
        sequence = "".join(chunks).upper()
        seqflags = []
        if proteins and sequence.endswith("*"):
            sequence = sequence[:-1]
            seqflags.append("terminal_stop_removed")
        if not sequence:
            raise PreparationError(f"Empty FASTA sequence {identifier!r}: {path}")
        if proteins:
            illegal = set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBZXJUO")
            if illegal:
                raise PreparationError(f"Protein {identifier}: internal stop or illegal residues {sorted(illegal)}")
            if set(sequence) & set("BZXJ"):
                seqflags.append("ambiguous_amino_acids")
        else:
            illegal = set(sequence) - set("ACGTURYSWKMBDHVN")
            if illegal:
                raise PreparationError(f"Assembly {identifier}: illegal nucleotides {sorted(illegal)}")
        sequences[identifier] = sequence
        flags[identifier] = seqflags

    with _open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            if line.startswith(">"):
                finish()
                parts = line[1:].split()
                if not parts:
                    raise PreparationError(f"Empty FASTA header: {path}:{line_number}")
                identifier, chunks = parts[0], []
            else:
                if identifier is None:
                    raise PreparationError(f"FASTA sequence before header: {path}:{line_number}")
                chunks.append("".join(line.split()))
    finish()
    if not sequences:
        raise PreparationError(f"Empty FASTA: {path}")
    return sequences, flags


def _attributes(text: str, path: Path, line: int) -> dict[str, list[str]]:
    attributes = {}
    if text == ".":
        return attributes
    for item in text.rstrip(";").split(";"):
        if not item:
            continue
        if "=" not in item:
            raise PreparationError(f"Expected GFF3 key=value attribute: {path}:{line}")
        key, value = item.split("=", 1)
        if key in attributes:
            raise PreparationError(f"Repeated GFF3 attribute {key}: {path}:{line}")
        attributes[key] = [unquote(v) for v in value.split(",")]
    return attributes


def _read_gff(path: Path, proteins: dict[str, str]):
    nodes, genes, checks, protein_features = {}, {}, [], []
    with _open(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise PreparationError(f"Expected nine GFF3 columns: {path}:{number}")
            seqid, _, kind, start, end, _, strand, _, attribute_text = fields
            seqid = _identifier(unquote(seqid), "GFF seqid")
            try:
                start, end = int(start), int(end)
            except ValueError as exc:
                raise PreparationError(f"Noninteger GFF coordinate: {path}:{number}") from exc
            if start < 1 or end < start or strand not in {"+", "-", ".", "?"}:
                raise PreparationError(f"Invalid GFF coordinates/strand: {path}:{number}")
            attributes = _attributes(attribute_text, path, number)
            ids, parents = attributes.get("ID", []), attributes.get("Parent", [])
            # A polypeptide may use Derives_from instead of Parent.
            if not parents and kind.lower() == "polypeptide":
                parents = attributes.get("Derives_from", [])
            if len(ids) > 1:
                raise PreparationError(f"Multiple GFF IDs on a feature: {path}:{number}")
            identifier = ids[0] if ids else ""
            is_gene = kind.lower() in {"gene", "pseudogene"}
            is_transcript = kind.lower() in {"mrna", "transcript"} or kind.lower().endswith("rna")
            if (is_gene or is_transcript) and not identifier:
                raise PreparationError(f"Gene/transcript lacks ID: {path}:{number}")
            if identifier:
                _identifier(identifier, "GFF feature ID")
                node = dict(seqid=seqid, start=start, end=end, strand=strand, kind=kind,
                            parents=tuple(parents), is_gene=is_gene, is_transcript=is_transcript)
                if identifier in nodes:
                    old = nodes[identifier]
                    if is_gene or is_transcript or any(old[k] != node[k] for k in
                            ("seqid", "strand", "kind", "parents")):
                        raise PreparationError(f"Conflicting/duplicate GFF ID: {identifier}")
                    # GFF3 permits a discontinuous CDS to repeat its ID.
                    old["start"], old["end"] = min(start, old["start"]), max(end, old["end"])
                else:
                    nodes[identifier] = node
                if is_gene:
                    biotype = next((attributes[k][0] for k in
                                    ("gene_biotype", "gene_type", "biotype") if k in attributes), "")
                    genes[identifier] = dict(gene_id=identifier, seqid=seqid, start=start,
                                             end=end, strand=strand, gene_type=kind, biotype=biotype)
            checks.append((identifier, tuple(parents), seqid, start, end, number, kind))
            if attributes.get("protein_id"):
                protein_features.append((attributes["protein_id"], (identifier,) if identifier else tuple(parents)))

    if not genes:
        raise PreparationError(f"No gene features found in {path}; supply genes_tsv for a normalized annotation")
    ancestors_cache = {}

    def ancestors(identifier: str, active: frozenset[str] = frozenset()) -> set[str]:
        if identifier in ancestors_cache:
            return ancestors_cache[identifier]
        if identifier not in nodes:
            raise PreparationError(f"GFF Parent/Derives_from does not exist: {identifier}")
        if identifier in active:
            raise PreparationError(f"Cycle in GFF Parent chain at {identifier}")
        node = nodes[identifier]
        result = {identifier} if node["is_gene"] else set()
        for parent in node["parents"]:
            result.update(ancestors(parent, active | {identifier}))
        ancestors_cache[identifier] = result
        return result

    for identifier in nodes:
        ancestors(identifier)
    coding_genes, transcripts = set(), defaultdict(set)
    for identifier, parents, seqid, start, end, number, kind in checks:
        roots = set()
        for parent in parents:
            roots.update(ancestors(parent))
            parent_node = nodes[parent]
            if seqid != parent_node["seqid"] or start < parent_node["start"] or end > parent_node["end"]:
                raise PreparationError(f"Child outside parent bounds: {path}:{number} / {parent}")
        if identifier and nodes[identifier]["is_transcript"]:
            if len(ancestors(identifier)) != 1:
                raise PreparationError(f"Transcript must resolve to one gene: {identifier}")
            transcripts[next(iter(ancestors(identifier)))].add(identifier)
        if kind.lower() in {"mrna", "cds"}:
            coding_genes.update(roots)
    for gene_id, gene in genes.items():
        if gene["biotype"].lower() in {"protein_coding", "protein-coding"}:
            coding_genes.add(gene_id)

    mapping = defaultdict(lambda: {"genes": set(), "transcripts": set(), "sources": set()})

    def bind(protein_id: str, identifiers: tuple[str, ...], source: str):
        if protein_id not in proteins:
            return
        roots, tx = set(), set()
        for identifier in identifiers:
            roots.update(ancestors(identifier))
            if nodes[identifier]["is_transcript"]:
                tx.add(identifier)
            for parent in nodes[identifier]["parents"]:
                if nodes[parent]["is_transcript"]:
                    tx.add(parent)
        if len(roots) != 1:
            raise PreparationError(f"Protein mapping is ambiguous or lacks gene: {protein_id} -> {sorted(roots)}")
        mapping[protein_id]["genes"].update(roots)
        mapping[protein_id]["transcripts"].update(tx)
        mapping[protein_id]["sources"].add(source)

    for identifier in nodes.keys() & proteins.keys():
        bind(identifier, (identifier,), "exact_GFF_feature_ID")
    for protein_ids, identifiers in protein_features:
        if not identifiers:
            raise PreparationError("protein_id feature needs ID or Parent")
        for protein_id in protein_ids:
            bind(protein_id, identifiers, "GFF_protein_id")
    return genes, mapping, coding_genes, transcripts, nodes, ancestors


def _read_normalized(path: Path, proteins: dict[str, str]):
    rows = _table(path, {"gene_id", "seqid", "start", "end", "strand", "protein_id"})
    genes, coding_genes, transcripts = {}, set(), defaultdict(set)
    mapping = defaultdict(lambda: {"genes": set(), "transcripts": set(), "sources": set()})
    seen = set()
    for row in rows:
        gene_id = _identifier(row["gene_id"], "gene_id")
        seqid = _identifier(row["seqid"], "seqid")
        try:
            start, end = int(row["start"]), int(row["end"])
        except ValueError as exc:
            raise PreparationError(f"Noninteger normalized coordinates for {gene_id}") from exc
        if start < 1 or end < start or row["strand"] not in {"+", "-", ".", "?"}:
            raise PreparationError(f"Invalid normalized coordinates/strand for {gene_id}")
        gene = dict(gene_id=gene_id, seqid=seqid, start=start, end=end, strand=row["strand"],
                    gene_type=row.get("gene_type", "gene"), biotype=row.get("biotype", ""))
        if gene_id in genes and genes[gene_id] != gene:
            raise PreparationError(f"Conflicting normalized gene rows: {gene_id}")
        protein_id, transcript_id = row["protein_id"], row.get("transcript_id", "")
        key = gene_id, transcript_id, protein_id
        if key in seen:
            raise PreparationError(f"Duplicate normalized gene/protein row: {key}")
        seen.add(key)
        genes[gene_id] = gene
        if transcript_id:
            transcripts[gene_id].add(transcript_id)
        if protein_id:
            if protein_id not in proteins:
                raise PreparationError(f"Normalized protein ID not in FASTA: {protein_id}")
            mapping[protein_id]["genes"].add(gene_id)
            if transcript_id:
                mapping[protein_id]["transcripts"].add(transcript_id)
            mapping[protein_id]["sources"].add("normalized_genes_tsv")
        if gene["biotype"].lower() in {"protein_coding", "protein-coding"} or protein_id:
            coding_genes.add(gene_id)
    if not genes:
        raise PreparationError(f"Empty normalized annotation: {path}")
    return genes, mapping, coding_genes, transcripts


def _lengths(spec: dict, genes: dict) -> tuple[dict[str, int], str]:
    selectors = [key for key in ("assembly_fasta", "lengths_tsv", "seq_lengths") if spec.get(key) is not None]
    if len(selectors) > 1:
        raise PreparationError("Use exactly one of assembly_fasta, lengths_tsv, or seq_lengths")
    if not selectors:
        lengths = defaultdict(int)
        for gene in genes.values():
            lengths[gene["seqid"]] = max(lengths[gene["seqid"]], gene["end"])
        return dict(lengths), "inferred_not_assembly"
    key = selectors[0]
    if key == "assembly_fasta":
        # Only lengths are required: avoid retaining a multi-gigabase assembly.
        lengths, current = {}, None
        with _open(Path(spec[key])) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                if line.startswith(">"):
                    parts = line[1:].split()
                    if not parts or parts[0] in lengths:
                        raise PreparationError(f"Empty/duplicate assembly FASTA ID at line {line_number}")
                    current = parts[0]
                    lengths[current] = 0
                else:
                    if current is None:
                        raise PreparationError("Assembly sequence before first FASTA header")
                    sequence = "".join(line.split()).upper()
                    if set(sequence) - set("ACGTURYSWKMBDHVN"):
                        raise PreparationError(f"Illegal assembly nucleotide at line {line_number}: {current}")
                    lengths[current] += len(sequence)
        if not lengths:
            raise PreparationError("Assembly FASTA is empty")
    elif key == "lengths_tsv":
        lengths = {}
        for row in _table(Path(spec[key]), {"seqid", "length"}):
            if row["seqid"] in lengths:
                raise PreparationError(f"Duplicate sequence length: {row['seqid']}")
            lengths[_identifier(row["seqid"], "lengths seqid")] = row["length"]
    else:
        if not isinstance(spec[key], dict):
            raise PreparationError("seq_lengths must map sequence IDs to lengths")
        lengths = dict(spec[key])
    for seqid, length in lengths.items():
        _identifier(seqid, "lengths seqid")
        try:
            integer = int(length)
        except (ValueError, TypeError) as exc:
            raise PreparationError(f"Invalid sequence length: {seqid}={length}") from exc
        if isinstance(length, bool) or str(integer) != str(length) or integer <= 0:
            raise PreparationError(f"Sequence lengths must be positive integers: {seqid}={length}")
        lengths[seqid] = integer
    return lengths, "assembly_fasta" if key == "assembly_fasta" else "supplied_lengths"


def _tsv(fields: list[str], rows: list[dict]) -> bytes:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


def _noncoding_label(label: str) -> bool:
    label = label.lower().replace("-", "_")
    return ("pseudogene" in label or
            (label.endswith("rna") and label != "mrna") or
            label in {"noncoding", "non_coding", "mirna_primary_transcript", "pre_mirna"})


def prepare_genome(spec: dict, out_dir: Path) -> dict:
    """Validate a genome and write reusable input tables without family state.

    Existing identical results are reusable. Existing changed files, changed
    configuration, or changed input hashes are rejected instead of overwritten.
    """
    allowed = {"id", "species", "assembly_id", "annotation_id", "subgenome", "gff",
               "genes_tsv", "proteins", "id_map_tsv", "representative_map_tsv", "rank_policy",
               "assembly_fasta", "lengths_tsv", "seq_lengths"}
    unknown = set(spec) - allowed
    if unknown:
        raise PreparationError(f"Unknown genome options (family inputs are not accepted): {sorted(unknown)}")
    for key in ("id", "species", "assembly_id", "annotation_id"):
        _identifier(spec.get(key), key)
    if bool(spec.get("gff")) == bool(spec.get("genes_tsv")):
        raise PreparationError("Provide exactly one of gff or genes_tsv")
    if not spec.get("proteins"):
        raise PreparationError("proteins FASTA is required")
    rank_policy = spec.get("rank_policy", "protein_coding")
    if rank_policy not in {"all_genes", "protein_coding"}:
        raise PreparationError("rank_policy must be all_genes or protein_coding")
    out_dir = Path(out_dir)
    canonical = dict(spec, rank_policy=rank_policy, subgenome=spec.get("subgenome", ""))
    sources = {}
    for key in ("gff", "genes_tsv", "proteins", "id_map_tsv", "representative_map_tsv", "assembly_fasta", "lengths_tsv"):
        if spec.get(key):
            path = Path(spec[key]).resolve()
            if not path.is_file():
                raise PreparationError(f"Missing source {key}: {path}")
            canonical[key] = str(path)
            sources[key] = {"path": str(path), "sha256": _sha(path), "size_bytes": path.stat().st_size}
    code_hash = _sha(Path(__file__))
    fingerprint_payload = {"schema_version": SCHEMA_VERSION, "config": canonical,
                           "sources": sources, "module_sha256": code_hash}
    fingerprint = hashlib.sha256(json.dumps(fingerprint_payload, sort_keys=True).encode()).hexdigest()
    manifest_path = out_dir / "input_manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old.get("input_fingerprint") != fingerprint:
            raise PreparationError(f"Existing inputs differ; choose a new output directory: {out_dir}")
        if set(old.get("outputs", {})) != {"genes.tsv", "all_isoforms.tsv", "proteins.faa", "positions.gff", "chromosomes.tsv"}:
            raise PreparationError(f"Incomplete prepared-output manifest: {manifest_path}")
        for name, info in old.get("outputs", {}).items():
            path = out_dir / name
            if not path.is_file() or _sha(path) != info["sha256"]:
                raise PreparationError(f"Existing prepared output is missing or changed: {path}")
        return old

    proteins, seqflags = _fasta(Path(spec["proteins"]))
    if spec.get("gff"):
        genes, mapping, coding, transcripts, nodes, ancestors = _read_gff(Path(spec["gff"]), proteins)
    else:
        genes, mapping, coding, transcripts = _read_normalized(Path(spec["genes_tsv"]), proteins)
        nodes, ancestors = None, None
    if spec.get("id_map_tsv"):
        seen = set()
        for row in _table(Path(spec["id_map_tsv"]), {"gene_id", "transcript_id", "protein_id"}):
            gene_id, transcript_id, protein_id = row["gene_id"], row["transcript_id"], row["protein_id"]
            if (gene_id, transcript_id, protein_id) in seen:
                raise PreparationError(f"Duplicate explicit ID mapping: {gene_id}/{protein_id}")
            seen.add((gene_id, transcript_id, protein_id))
            if gene_id not in genes or protein_id not in proteins:
                raise PreparationError(f"Explicit ID map refers to absent gene/protein: {gene_id}/{protein_id}")
            if transcript_id and nodes is not None:
                if transcript_id not in nodes or not nodes[transcript_id]["is_transcript"] or ancestors(transcript_id) != {gene_id}:
                    raise PreparationError(f"Explicit transcript does not belong to gene: {transcript_id}/{gene_id}")
            mapping[protein_id]["genes"].add(gene_id)
            mapping[protein_id]["transcripts"].update({transcript_id} if transcript_id else set())
            mapping[protein_id]["sources"].add("explicit_id_map")
            if transcript_id:
                transcripts[gene_id].add(transcript_id)
    unmapped = set(proteins) - set(mapping)
    if unmapped:
        raise PreparationError(f"{len(unmapped)} protein FASTA IDs cannot be mapped; supply id_map_tsv. Examples: {sorted(unmapped)[:8]}")
    isoforms = defaultdict(list)
    for protein_id, info in mapping.items():
        if len(info["genes"]) != 1:
            raise PreparationError(f"Protein maps to multiple loci: {protein_id} -> {sorted(info['genes'])}")
        gene_id = next(iter(info["genes"]))
        isoforms[gene_id].append(protein_id)
        coding.add(gene_id)
    explicitly_noncoding = set()
    for gene_id, gene in genes.items():
        labels = [gene["gene_type"], gene["biotype"]]
        if nodes is not None:
            labels += [nodes[tx]["kind"] for tx in transcripts[gene_id]]
        if any(_noncoding_label(label) for label in labels):
            explicitly_noncoding.add(gene_id)
    # Unknown gene-only loci stay in the rank universe so absent biotype or
    # translation evidence cannot create artificial adjacency of homologs.
    unresolved_biotype = set(genes) - coding - explicitly_noncoding
    coding_conflicts = coding & explicitly_noncoding
    all_gene_count = len(genes)
    lengths, length_source = _lengths(spec, genes)
    for gene_id, gene in genes.items():
        if gene["seqid"] not in lengths or gene["end"] > lengths[gene["seqid"]]:
            raise PreparationError(f"Gene outside supplied assembly lengths: {gene_id} ({gene['seqid']}:{gene['start']}-{gene['end']})")
    if rank_policy == "protein_coding":
        genes = {gene_id: gene for gene_id, gene in genes.items()
                 if gene_id in coding or gene_id in unresolved_biotype}
    if not genes:
        raise PreparationError("No genes eligible under the selected rank policy")
    for gene_id, gene in genes.items():
        if gene["start"] > MCSCANX_MAX_COORDINATE or gene["end"] > MCSCANX_MAX_COORDINATE:
            raise PreparationError(
                f"Gene coordinate exceeds MCScanX signed 32-bit limit ({MCSCANX_MAX_COORDINATE}): "
                f"{gene_id} ({gene['seqid']}:{gene['start']}-{gene['end']}). "
                "Use a validated alternative coordinate/engine workflow; coordinates were not truncated.")

    frozen = {}
    if spec.get("representative_map_tsv"):
        for row in _table(Path(spec["representative_map_tsv"]), {"gene_id", "protein_id"}):
            gene_id, protein_id = row["gene_id"], row["protein_id"]
            if gene_id in frozen or gene_id not in genes or protein_id not in isoforms[gene_id]:
                raise PreparationError(f"Invalid/duplicate frozen representative: {gene_id}/{protein_id}")
            frozen[gene_id] = protein_id
    namespace = "g" + hashlib.sha256(spec["id"].encode()).hexdigest()[:8]
    engine_ids = {gene_id: f"{namespace}_gene{i:08d}" for i, gene_id in enumerate(sorted(genes), 1)}
    engine_seqids = {seqid: f"{namespace}_chr{i:06d}" for i, seqid in enumerate(sorted(lengths), 1)}
    groups = defaultdict(list)
    gene_rows, isoform_rows, selected, shorter = [], [], {}, 0
    for gene_id, gene in genes.items():
        choices = sorted(isoforms[gene_id], key=lambda protein_id: (-len(proteins[protein_id]), protein_id))
        chosen = frozen.get(gene_id, choices[0] if choices else "")
        if chosen:
            selected[gene_id] = chosen
        shorter += bool(chosen and len(proteins[chosen]) < len(proteins[choices[0]]))
        flags = list(seqflags.get(chosen, []))
        if not chosen:
            flags.append("no_protein")
        if gene["strand"] in {".", "?"}:
            flags.append("unknown_strand")
        if "pseudo" in (gene["gene_type"] + gene["biotype"]).lower():
            flags.append("pseudogene_annotation")
        if gene_id in unresolved_biotype:
            flags.append("unknown_biotype_retained")
        if gene_id in explicitly_noncoding:
            flags.append("noncoding_annotation")
        if gene_id in coding_conflicts:
            flags.append("coding_evidence_with_noncoding_annotation")
        row = dict(genome_id=spec["id"], species=spec["species"], assembly_id=spec["assembly_id"],
                   annotation_id=spec["annotation_id"], subgenome=spec.get("subgenome", ""),
                   seqid=gene["seqid"], engine_seqid=engine_seqids[gene["seqid"]], gene_id=gene_id,
                   protein_id=chosen, engine_id=engine_ids[gene_id], start=gene["start"], end=gene["end"],
                   strand=gene["strand"], has_protein=int(bool(chosen)), protein_length=len(proteins.get(chosen, "")),
                   sequence_sha256=hashlib.sha256(proteins[chosen].encode()).hexdigest() if chosen else "",
                   representative_reason="explicit_frozen_map" if gene_id in frozen else
                       "longest_protein_lexical_tiebreak" if chosen else "no_available_protein",
                   model_flags=";".join(flags))
        groups[row["seqid"]].append(row)
        represented_tx = set()
        for protein_id in sorted(choices):
            info = mapping[protein_id]
            represented_tx.update(info["transcripts"])
            isoform_rows.append(dict(genome_id=spec["id"], gene_id=gene_id,
                transcript_id=";".join(sorted(info["transcripts"])), protein_id=protein_id,
                engine_id=engine_ids[gene_id], is_representative=int(protein_id == chosen),
                protein_length=len(proteins[protein_id]), sequence_sha256=hashlib.sha256(proteins[protein_id].encode()).hexdigest(),
                mapping_source=";".join(sorted(info["sources"])), model_flags=";".join(seqflags[protein_id])))
        for transcript_id in sorted(transcripts[gene_id] - represented_tx):
            isoform_rows.append(dict(genome_id=spec["id"], gene_id=gene_id, transcript_id=transcript_id,
                protein_id="", engine_id=engine_ids[gene_id], is_representative=0, protein_length=0,
                sequence_sha256="", mapping_source="annotation_without_protein", model_flags="no_protein"))
    for seqid in sorted(groups):
        rows = sorted(groups[seqid], key=lambda row: (row["start"], row["end"], row["gene_id"]))
        for rank, row in enumerate(rows, 1):
            row["annotation_rank"] = rank
        for rank, row in enumerate(sorted(rows, key=lambda row: (row["start"], row["engine_id"])), 1):
            row["engine_rank"] = rank
        gene_rows.extend(rows)
    chromosomes = [dict(seqid=seqid, engine_seqid=engine_seqids[seqid], length=lengths[seqid], genome_id=spec["id"])
                   for seqid in sorted(lengths)]
    fasta_lines = []
    for row in gene_rows:
        if row["has_protein"]:
            sequence = proteins[row["protein_id"]]
            fasta_lines.append(">" + row["engine_id"] + "\n" + "\n".join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + "\n")
    artifacts = {
        "genes.tsv": _tsv(GENE_FIELDS, gene_rows),
        "all_isoforms.tsv": _tsv(ISOFORM_FIELDS, sorted(isoform_rows, key=lambda row: (row["gene_id"], row["protein_id"], row["transcript_id"]))),
        "proteins.faa": "".join(fasta_lines).encode(),
        "positions.gff": "".join(f"{row['engine_seqid']}\t{row['engine_id']}\t{row['start']}\t{row['end']}\n" for row in gene_rows).encode(),
        "chromosomes.tsv": _tsv(["seqid", "engine_seqid", "length", "genome_id"], chromosomes),
    }
    warnings = []
    if length_source == "inferred_not_assembly":
        warnings.append("Sequence lengths inferred from maximum gene end: incomplete physical chromosome extents; provide assembly lengths for chromosome-scale figures.")
    no_protein = len(genes) - len(selected)
    if no_protein:
        warnings.append(f"{no_protein} eligible loci lack protein; kept in both annotation and MCScanX coordinate ranks.")
    if unresolved_biotype:
        warnings.append(f"{len(unresolved_biotype)} loci have unresolved coding biotype; retained and flagged to prevent artificial adjacency. Review annotation before interpreting local duplication.")
    if coding_conflicts:
        warnings.append(f"{len(coding_conflicts)} loci have coding evidence and a noncoding/pseudogene label; retained with a conflict flag for annotation review.")
    tied_starts = sum(len(rows) > 1 for rows in _group_starts(gene_rows).values())
    rank_differences = sum(row["annotation_rank"] != row["engine_rank"] for row in gene_rows)
    if tied_starts:
        warnings.append(f"{tied_starts} same-start groups; {rank_differences} loci have different annotation/engine rank. Both rank policies are exported.")
    if shorter:
        warnings.append(f"{shorter} explicit frozen representatives are shorter than the longest available protein.")
    manifest = dict(schema_version=SCHEMA_VERSION, input_fingerprint=fingerprint,
        module_sha256=code_hash, config=canonical, sources=sources, genome_id=spec["id"],
        namespace=namespace, coordinate_system="1-based inclusive", length_source=length_source,
        rank_policy=rank_policy, annotation_sort=["start", "end", "gene_id"], engine_sort=["start", "engine_id"],
        counts=dict(annotated_genes=all_gene_count, eligible_genes=len(genes), representative_proteins=len(selected),
                    missing_protein_genes=no_protein, input_proteins=len(proteins), sequences=len(lengths),
                    frozen_representatives=len(frozen), frozen_shorter_than_longest=shorter,
                    same_start_groups=tied_starts, rank_different_loci=rank_differences,
                    excluded_noncoding_loci=all_gene_count - len(genes),
                    unknown_biotype_retained_loci=len(unresolved_biotype),
                    coding_evidence_annotation_conflicts=len(coding_conflicts),
                    terminal_stops_removed=sum("terminal_stop_removed" in flags for flags in seqflags.values())),
        outputs={name: {"sha256": hashlib.sha256(content).hexdigest(), "size_bytes": len(content)} for name, content in artifacts.items()},
        warnings=warnings)
    artifacts["input_manifest.json"] = (json.dumps(manifest, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode()
    # Validate every existing file before writing anything, preventing partial
    # updates if a prior non-manifest output was changed or reused incorrectly.
    for name, content in artifacts.items():
        path = out_dir / name
        if path.exists() and path.read_bytes() != content:
            raise PreparationError(f"Refusing to replace differing prepared file: {path}")
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, content in artifacts.items():
        path = out_dir / name
        if not path.exists():
            with path.open("xb") as handle:
                handle.write(content)
    return manifest


def _group_starts(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["seqid"], row["start"]].append(row)
    return grouped
