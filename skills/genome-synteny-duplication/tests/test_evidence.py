"""Scientific invariants for the generic evidence module, without external tools."""
import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from synteny_evidence import process_comparison, overlay_families


def gene(genome,name,rank,start=None,end=None,protein=True,seqid='chr1'):
    start = start if start is not None else rank*100
    return dict(genome_id=genome,species=genome,assembly_id='v1',annotation_id='ann1',subgenome='',
                seqid=seqid,engine_seqid=f'{genome}|{seqid}',gene_id=name,protein_id=f'{name}.p' if protein else '',
                engine_id=f'{genome}|{name}',start=start,end=end if end is not None else start+49,strand='+',
                annotation_rank=rank,engine_rank=rank,has_protein=int(protein),protein_length=100 if protein else 0,
                sequence_sha256='same_sequence',representative_reason='longest',model_flags='')


def m8(q,s,bits=100,evalue='1e-20',qs=1,qe=100,ss=1,se=100,alen=100):
    return f'{q}\t{s}\t80\t{alen}\t0\t0\t{qs}\t{qe}\t{ss}\t{se}\t{evalue}\t{bits}\n'


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def hashes(folder):
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(folder).iterdir() if p.is_file()}


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.genes = [gene('A','g0',1,50,900,False),gene('A','g1',2,200,249),
                      gene('A','g2',3,400,449,False),gene('A','g3',4,600,649),gene('A','g4',5,800,849)]
        self.hit = self.root/'hits.m8'
        self.hit.write_text(m8('A|g1','A|g3',qs=1,qe=20,alen=90)+
                            m8('A|g3','A|g4',qs=1,qe=20,alen=90)+
                            m8('A|g4','A|g3',bits=200,qs=1,qe=60,ss=1,se=80,alen=90)+
                            m8('A|g4','A|g4'))
        self.col = self.root/'test.collinearity'
        self.col.write_text('## Alignment 0: score=100 e_value=1e-30 N=2 A|chr1&A|chr1 plus\n'
                            '0- 0: A|g1 A|g3 1e-20\n0- 1: A|g3 A|g4 1e-20\n')
        self.native = self.root/'test.gene_type'
        self.native.write_text('A|g1\t4\nA|g3\t4\nA|g4\t4\nA|g0\t0\nA|g2\t0\n')
        self.comp = dict(id='A_self',a='A',b='A',kind='self')

    def tearDown(self):
        self.temp.cleanup()

    def run_fixture(self):
        folder = self.root/'results'
        report = process_comparison(self.genes,[self.hit],self.col,folder,self.comp,classifier_path=self.native)
        return folder,report

    def test_missing_protein_separates_tandem_and_native_class_does_not_override(self):
        folder,report = self.run_fixture()
        tandem = read(folder/'tandem_pairs.tsv')
        self.assertEqual([(r['gene_a'],r['gene_b']) for r in tandem],[('g3','g4')])
        proximal = read(folder/'proximal_pairs.tsv')
        self.assertEqual([(r['gene_a'],r['gene_b'],r['intervening_genes']) for r in proximal],[('g1','g3','1')])
        classes = {r['gene_id']:r for r in read(folder/'gene_duplication_classification.tsv')}
        self.assertEqual(classes['g3']['native_code'],'4')
        self.assertEqual(classes['g3']['has_independent_tandem'],'1')
        self.assertEqual(classes['g3']['has_collinear_anchor'],'1')
        self.assertEqual(classes['g2']['assessment_status'],'not_assessable_no_protein')
        self.assertEqual(classes['g2']['native_code'],'')
        self.assertEqual(report['counts']['tandem_arrays'],1)

    def test_best_single_hsp_coordinates_and_reciprocal_evidence(self):
        folder,_ = self.run_fixture()
        homologs = {(r['gene_a'],r['gene_b']):r for r in read(folder/'homology_pairs.tsv')}
        best = homologs[('g3','g4')]
        self.assertEqual(best['best_hsp_query'],'A|g4')
        self.assertEqual(best['observed_directions'],'2')
        self.assertAlmostEqual(float(best['coverage_a']),.8)
        self.assertAlmostEqual(float(best['coverage_b']),.6)
        low = homologs[('g1','g3')]
        self.assertAlmostEqual(float(low['coverage_a']),.2)  # alen was 90, so alen/L is wrong
        self.assertIn('coverage_below_0.5',low['review_flags'])
        self.assertEqual(len(homologs),2)  # true self hit dropped; identical distinct-locus peptides retained

    def test_rank_membership_nonanchors_and_physical_overlap_separate(self):
        folder,_ = self.run_fixture()
        side_a = {r['gene_id']:r for r in read(folder/'all_block_genes.tsv') if r['side']=='A'}
        self.assertEqual(set(side_a),{'g1','g2','g3'})
        self.assertEqual(side_a['g2']['is_anchor'],'0')
        self.assertEqual(side_a['g1']['anchor_partner_ids'],'g3')
        overlaps = read(folder/'block_physical_overlaps.tsv')
        self.assertTrue(any(r['gene_id']=='g0' and r['side']=='A' and r['in_rank_interval']=='0' for r in overlaps))

    def test_cross_native_reverse_sides_canonicalized_and_tandem_not_applicable(self):
        genes = [gene('A','x',1),gene('A','y',2),gene('B','x',1),gene('B','y',2)]
        hits = self.root/'cross.m8'
        hits.write_text(m8('B|x','A|y')+m8('A|y','B|x')+m8('B|y','A|x'))
        col = self.root/'cross.collinearity'
        col.write_text('## Alignment 7: score=99 e_value=1e-20 N=2 B|chr1&A|chr1 minus\n'
                       '7- 0: B|x A|y 1e-20\n7- 1: B|y A|x 1e-20\n')
        folder = self.root/'cross'
        report = process_comparison(genes,[hits],col,folder,dict(id='AB',a='A',b='B',kind='cross'))
        block = read(folder/'blocks.tsv')[0]
        self.assertEqual((block['genome_a'],block['genome_b'],block['orientation']),('A','B','minus'))
        anchor = read(folder/'anchor_pairs.tsv')[0]
        self.assertEqual((anchor['gene_a'],anchor['gene_b']),('y','x'))
        self.assertEqual(report['statuses']['tandem'],'not_applicable_cross_comparison')
        self.assertEqual(read(folder/'tandem_pairs.tsv'),[])
        self.assertEqual(report['counts']['unique_homology_pairs'],2)

    def test_array_joins_only_adjacent_supported_edges(self):
        genes = [gene('A',f'g{i}',i) for i in range(1,5)]
        self.hit.write_text(m8('A|g1','A|g2')+m8('A|g2','A|g3')+m8('A|g2','A|g4'))
        folder = self.root/'arrays'
        process_comparison(genes,[self.hit],None,folder,self.comp)
        arrays = read(folder/'tandem_arrays.tsv')
        self.assertEqual(len(arrays),1)
        self.assertEqual(arrays[0]['gene_ids'],'g1;g2;g3')
        self.assertEqual((arrays[0]['gene_count'],arrays[0]['edge_count']),('3','2'))
        self.assertEqual(len(read(folder/'proximal_pairs.tsv')),1)

    def test_explicit_proximal_intervening_boundary(self):
        genes = [gene('A',f'g{i}',i) for i in range(1,6)]
        self.hit.write_text(m8('A|g1','A|g4')+m8('A|g1','A|g5'))
        folder = self.root/'prox'
        process_comparison(genes,[self.hit],None,folder,self.comp,proximal_max_intervening=2)
        self.assertEqual([r['gene_b'] for r in read(folder/'proximal_pairs.tsv')],['g4'])
        self.assertEqual(len(read(folder/'homology_pairs.tsv')),2)

    def test_family_overlay_multilabel_isoform_and_core_immutable(self):
        folder,_ = self.run_fixture()
        before = hashes(folder)
        families = [dict(genome_id='A',protein_id='g3.short',family_label='CAD'),
                    dict(genome_id='A',gene_id='g3',families='MDR;CAD'),
                    dict(genome_id='A',gene_id='g4',families='MDR'),
                    dict(genome_id='A',gene_id='g2',family_label='NEW')]
        out = self.root/'overlay'
        report = overlay_families([folder],self.genes,families,out,[dict(genome_id='A',gene_id='g3',protein_id='g3.short')])
        self.assertFalse(report['core_inference_modified'])
        self.assertEqual(before,hashes(folder))
        tandem = read(out/'family_tandem_pairs.tsv')[0]
        self.assertEqual((tandem['families_a'],tandem['shared_families'],tandem['same_family_pair']),('CAD;MDR','MDR','1'))
        self.assertTrue(any(r['gene_id']=='g2' and r['is_anchor']=='0' for r in read(out/'family_block_members.tsv')))
        self.assertEqual(report['family_loci'],3)

    def test_gene_protein_disagreement_and_unmapped_ids_rejected(self):
        folder,_ = self.run_fixture()
        with self.assertRaisesRegex(ValueError,'disagreement'):
            overlay_families([folder],self.genes,[dict(genome_id='A',gene_id='g1',protein_id='g3.p',family_label='F')],self.root/'bad')
        self.hit.write_text(m8('A|missing','A|g3'))
        with self.assertRaisesRegex(ValueError,'absent'):
            process_comparison(self.genes,[self.hit],None,self.root/'bad2',self.comp)

    def test_zero_results_headered_and_no_run_distinct_from_empty_run(self):
        folder = self.root/'zero'
        self.hit.write_text('')
        self.col.write_text('# no blocks\n')
        report = process_comparison(self.genes,[self.hit],self.col,folder,self.comp)
        self.assertEqual(report['statuses']['tandem'],'complete')
        self.assertEqual(report['counts']['tandem_pairs'],0)
        self.assertEqual(report['statuses']['collinearity'],'complete')
        self.assertTrue((folder/'tandem_pairs.tsv').read_text().startswith('comparison_id\t'))
        report = process_comparison(self.genes,[],None,self.root/'notrun',self.comp)
        self.assertEqual(report['statuses']['tandem'],'not_run')
        self.assertEqual(report['statuses']['collinearity'],'not_run')

    def test_no_searchable_proteins_is_not_assessable_not_a_negative_result(self):
        genes = [gene('A','g1',1,protein=False),gene('A','g2',2,protein=False)]
        self.hit.write_text('')
        self.col.write_text('# No selected hits; engine was not invoked.\n')
        folder = self.root/'proteinless'
        report = process_comparison(genes,[self.hit],self.col,folder,self.comp)
        for stage in ['homology','collinearity','tandem','proximal','classifier']:
            self.assertEqual(report['statuses'][stage],'not_assessable_no_proteins',stage)
        self.assertEqual(report['counts']['tandem_pairs'],0)
        self.assertEqual(read(folder/'homology_pairs.tsv'),[])
        self.assertTrue((folder/'homology_pairs.tsv').read_text().startswith('comparison_id\t'))
        self.assertEqual({r['assessment_status'] for r in read(folder/'gene_duplication_classification.tsv')},{'not_assessable_no_protein'})

    def test_cross_with_one_unsearchable_genome_is_not_assessable(self):
        genes = [gene('A','a',1),gene('B','b',1,protein=False)]
        self.hit.write_text('')
        report = process_comparison(genes,[self.hit],None,self.root/'half_proteinless',dict(id='AB',a='A',b='B',kind='cross'))
        self.assertEqual(report['statuses']['homology'],'not_assessable_no_proteins')
        self.assertEqual(report['statuses']['collinearity'],'not_assessable_no_proteins')
        self.assertEqual(report['statuses']['tandem'],'not_applicable_cross_comparison')

    def test_broken_rank_and_block_count_rejected(self):
        genes = [dict(g) for g in self.genes]
        genes[-1]['annotation_rank'] = 3
        with self.assertRaisesRegex(ValueError,'Duplicate annotation_rank'):
            process_comparison(genes,[],None,self.root/'bad_rank',self.comp)
        self.col.write_text('## Alignment 0: score=100 e_value=1e-30 N=3 A|chr1&A|chr1 plus\n0- 0: A|g1 A|g3 1e-20\n')
        with self.assertRaisesRegex(ValueError,'anchor count'):
            process_comparison(self.genes,[],self.col,self.root/'bad_count',self.comp)

    def test_same_start_engine_rank_is_separate_from_annotation_rank(self):
        genes = [gene('A','long',1,100,600),gene('A','short',2,100,200),gene('A','last',3,900,949)]
        genes[0]['engine_rank'],genes[1]['engine_rank'] = 2,1
        col = self.root/'ties.collinearity'
        col.write_text('## Alignment 0: score=100 e_value=1e-30 N=1 A|chr1&A|chr1 plus\n0- 0: A|long A|last 1e-20\n')
        folder = self.root/'ties'
        process_comparison(genes,[],col,folder,self.comp)
        side_a = [r['gene_id'] for r in read(folder/'all_block_genes.tsv') if r['side']=='A']
        self.assertEqual(side_a,['long'])
        self.assertTrue(any(r['gene_id']=='short' and r['side']=='A' for r in read(folder/'block_physical_overlaps.tsv')))

    def test_missing_native_rows_are_incomplete_not_singletons(self):
        self.native.write_text('A|g1\t4\n')
        folder,report = self.run_fixture()
        self.assertEqual(report['statuses']['classifier'],'incomplete_missing_native_rows')
        self.assertEqual(report['counts']['missing_native_protein_loci'],2)
        row = next(r for r in read(folder/'gene_duplication_classification.tsv') if r['gene_id']=='g3')
        self.assertEqual(row['assessment_status'],'missing_native_row')
        self.assertEqual(row['native_code'],'')

    def test_overlapping_adjacent_models_are_retained_for_review(self):
        genes = [gene('A','g1',1,100,600),gene('A','g2',2,200,400)]
        self.hit.write_text(m8('A|g1','A|g2'))
        folder = self.root/'overlap'
        process_comparison(genes,[self.hit],None,folder,self.comp)
        pair = read(folder/'tandem_pairs.tsv')[0]
        self.assertEqual(pair['quality_tier'],'review')
        self.assertEqual(pair['genes_overlap'],'1')
        self.assertIn('overlapping_gene_models',pair['review_flags'])


if __name__=='__main__':
    unittest.main()
