import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

P = Path(__file__).resolve().parents[1] / 'scripts' / 'synteny_workflow.py'
spec = importlib.util.spec_from_file_location('workflow', P)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)


class WorkflowTests(unittest.TestCase):
    def test_integer_gap_rejected(self):
        with self.assertRaises(ValueError):
            w.settings({'collinearity': {'gap_penalty': -0.5}})
        with self.assertRaises(ValueError):
            w.settings({'collinearity': {'min_anchors': True}})

    def test_cache_input_change_and_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            calls = []
            def builder(p):
                calls.append(p)
                (p/'evidence.tsv').write_text('a\tb\n')
                return {'status': 'no_detected_evidence'}
            p, reused = w.cache_step(d, {'input_sha': 'A', 'top': 5}, builder)
            self.assertFalse(reused)
            p2, reused = w.cache_step(d, {'input_sha': 'A', 'top': 5}, builder)
            self.assertTrue(reused)
            self.assertEqual(p, p2)
            p3, reused = w.cache_step(d, {'input_sha': 'A', 'top': 10}, builder)
            self.assertNotEqual(p, p3)
            self.assertEqual(len(calls), 2)
            (p/'evidence.tsv').write_text('tampered')
            with self.assertRaises(ValueError):
                w.cache_step(d, {'input_sha': 'A', 'top': 5}, builder)
            with self.assertRaises(ValueError):
                w.verify_complete(p)

    def test_failed_attempt_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            def fail(p):
                (p/'partial').write_text('evidence')
                raise RuntimeError('stopped')
            with self.assertRaises(RuntimeError):
                w.cache_step(d, {'a': 1}, fail)
            prior = list(Path(d).iterdir())[0]
            p, reused = w.cache_step(d, {'a': 1}, lambda p: (p/'done').write_text('ok') and {})
            self.assertNotEqual(prior, p)
            self.assertTrue((prior/'partial').exists())

    def test_distinct_targets_self_cap_and_best_hsp(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d)/'hits.m8'
            def row(q,s,score,e='1e-20'):
                return f'{q}\t{s}\t90\t100\t0\t0\t1\t100\t1\t100\t{e}\t{score}\n'
            raw.write_text(row('q','q',300)+row('q','identical_other_locus',300)+
                           row('q','s2',150)+row('q','s2',200)+row('q','s3',100))
            out = Path(d)/'selected'
            qc = w.select_hits([raw],out,2,1e-5,4)
            rows = [x.split('\t') for x in out.read_text().splitlines()]
            self.assertEqual([x[1] for x in rows],['identical_other_locus','s2'])
            self.assertEqual(rows[1][11],'200')
            self.assertEqual(qc['self_rows_removed'],1)
            self.assertEqual(qc['queries_at_search_cap'],1)

    def test_empty_selection_valid(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d)/'empty'
            raw.write_text('')
            out = Path(d)/'selected'
            self.assertEqual(w.select_hits([raw],out,5,1e-5,50)['selected_rows'],0)
            self.assertTrue(out.exists())

    def test_header_effective_gap(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'collinearity'
            col=w.settings({})[1]
            text='# MATCH_SCORE: 50\n# MATCH_SIZE: 5\n# GAP_PENALTY: -1\n# OVERLAP_WINDOW: 5\n# E_VALUE: 1e-05\n# MAX GAPS: 25\n'
            p.write_text(text)
            self.assertEqual(w.verify_header(p,col)['GAP_PENALTY'],-1)
            p.write_text(text.replace('GAP_PENALTY: -1','GAP_PENALTY: 0'))
            with self.assertRaises(ValueError):
                w.verify_header(p,col)


if __name__=='__main__':
    unittest.main()
