"""Focused figure contracts: source links, coordinate transforms, and no inventions."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "render_figures.py"
spec = importlib.util.spec_from_file_location("render_figures", SCRIPT)
figures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(figures)


def make_fixture(root):
    genes = []
    for genome in ("A", "B"):
        for name, start, end in (("g1", 101, 200), ("g2", 180, 250), ("g3", 401, 450)):
            genes.append({"genome_id": genome, "seqid": "chr1", "gene_id": name, "start": start, "end": end, "strand": "+", "species": genome, "assembly_id": "v1", "annotation_id": "v1"})
    figures.write_tsv(root / "genes.tsv", genes)
    figures.write_tsv(root / "chromosomes.tsv", [{"genome_id": g, "seqid": "chr1", "length": 1000} for g in ("A", "B")])
    figures.write_tsv(root / "families.tsv", [{"genome_id": g, "gene_id": name, "family_label": "F"} for g in ("A", "B") for name in ("g1", "g2")])
    blocks = [{"comparison_id": "A_B", "block_id": "B1", "genome_a": "A", "seqid_a": "chr1", "start_a": 101, "end_a": 450, "genome_b": "B", "seqid_b": "chr1", "start_b": 101, "end_b": 450, "anchor_count": 3, "orientation": "+"}]
    figures.write_tsv(root / "blocks.tsv", blocks)
    anchors = [{"comparison_id": "A_B", "block_id": "B1", "genome_a": "A", "gene_a": name, "genome_b": "B", "gene_b": name} for name in ("g1", "g3")]
    figures.write_tsv(root / "anchor_pairs.tsv", anchors)
    figures.write_tsv(root / "tandem_pairs.tsv", [{"genome_id": "A", "genome_a": "A", "gene_a": "g1", "genome_b": "A", "gene_b": "g2", "intervening_genes": 0, "quality_tier": "review", "review_flags": "overlap", "array_id": "A_T1"}])
    return SimpleNamespace(genes=[root / "genes.tsv"], chromosomes=[root / "chromosomes.tsv"], families=root / "families.tsv", results=[root], config=root / "config.json")


class CoordinateTests(unittest.TestCase):
    def test_reverse_preserves_physical_distance(self):
        forward = {"start": 100, "end": 1100}
        reverse = dict(forward, reverse=True)
        self.assertEqual(figures.position(forward, 100), 0)
        self.assertEqual(figures.position(reverse, 100), 1)
        self.assertAlmostEqual(figures.position(forward, 600), .5)
        for bp in (100, 200, 600, 1100):
            self.assertAlmostEqual(figures.position(forward, bp) + figures.position(reverse, bp), 1)
        self.assertEqual(figures.visible_strand(reverse, "+"), "-")
        self.assertEqual(figures.visible_strand(reverse, "."), ".")

    def test_overlap_lanes_are_inclusive(self):
        genes = [{"genome_id": "A", "gene_id": n, "start": s, "end": e} for n, s, e in (("long", 1, 100), ("inner", 10, 30), ("touch", 30, 70), ("after", 101, 110))]
        lanes = figures.assign_lanes(genes)
        self.assertNotEqual(lanes["A", "long"], lanes["A", "inner"])
        self.assertNotEqual(lanes["A", "inner"], lanes["A", "touch"])
        self.assertEqual(lanes["A", "after"], lanes["A", "long"])

    def test_native_inversion_keeps_engine_orientation(self):
        self.assertEqual(figures.block_orientation("minus"), "-")
        self.assertEqual(figures.block_orientation("plus"), "+")
        with self.assertRaises(ValueError):
            figures.block_orientation("unknown")


class SourceContractTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.args = make_fixture(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def test_namespaces_keep_identical_gene_ids_distinct(self):
        data = figures.FigureData(self.args, {})
        self.assertEqual(len(data.genes), 6)
        self.assertIn(("A", "g1"), data.genes)
        self.assertIn(("B", "g1"), data.genes)

    def test_repeated_support_deduplicates_only_edge(self):
        data = figures.FigureData(self.args, {})
        data.blocks.append(dict(data.blocks[0], block_id="B2"))
        data.block_map["A_B", "B2"] = data.blocks[-1]
        data.anchors.append(dict(data.anchors[0], block_id="B2"))
        pairs = data.unique_anchors()
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0]["block_id"], "B1;B2")
        self.assertEqual(data.omitted["duplicate_anchor_support_rows_collapsed"], 1)

    def test_background_threshold_does_not_hide_focal_anchor(self):
        data = figures.FigureData(self.args, {"min_block_anchors": 10})
        pairs = data.unique_anchors()
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["gene_a"], "g1")

    def test_unknown_family_gene_fails_instead_of_silent_drop(self):
        with self.args.families.open("a") as f:
            f.write("A\tunknown\tF\n")
        with self.assertRaisesRegex(ValueError, "Family gene missing"):
            figures.FigureData(self.args, {})

    def test_concurrent_source_edit_invalidates_figure(self):
        data = figures.FigureData(self.args, {})
        with self.args.families.open("a") as f:
            f.write("A\tg3\tF\n")
        with self.assertRaisesRegex(RuntimeError, "Source changed while plotting"):
            data.verify_sources_unchanged()

    def test_mixed_annotation_coordinates_fail(self):
        rows = figures.read_tsv(self.root / "anchor_pairs.tsv")
        rows[0]["start_a"] = "999"
        figures.write_tsv(self.root / "anchor_pairs.tsv", rows)
        with self.assertRaisesRegex(ValueError, "coordinate disagrees"):
            figures.FigureData(self.args, {})

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib runtime required")
    def test_local_no_invented_anchor_and_true_overlap_lane(self):
        config = {"kind": "local", "tracks": [{"genome_id": "A", "seqid": "chr1", "start": 1, "end": 1000, "reverse": True}, {"genome_id": "B", "seqid": "chr1", "start": 1, "end": 1000}], "dpi": 60}
        self.args.config.write_text(json.dumps(config))
        data = figures.FigureData(self.args, config)
        out = self.root / "out"
        out.mkdir()
        details = figures.render_local(data, out)
        figures.finalize(data, out, self.args, details)
        direct = [e for e in data.edges if e["relation_type"] == "direct_anchor"]
        tandem = [e for e in data.edges if e["relation_type"] == "strict_tandem"]
        # g2 occurs in both families and is close to g1; no g2-g2 link may appear.
        self.assertEqual({(e["gene_a"], e["gene_b"]) for e in direct}, {("g1", "g1"), ("g3", "g3")})
        self.assertEqual(len(tandem), 1)
        self.assertEqual(tandem[0]["quality_tier"], "review")
        self.assertNotEqual(data.drawn_genes["A", "g1"]["overlap_lane"], data.drawn_genes["A", "g2"]["overlap_lane"])
        self.assertEqual(data.drawn_genes["A", "g1"]["visible_strand"], "-")
        self.assertTrue(all((out / f"figure.{ext}").stat().st_size > 1000 for ext in ("svg", "png", "pdf")))
        self.assertEqual(len(figures.read_tsv(out / "figure_edges.tsv")), 3)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib runtime required")
    def test_window_without_anchor_keeps_zero_links(self):
        config = {"kind": "local", "tracks": [{"genome_id": g, "seqid": "chr1", "start": 205, "end": 290} for g in ("A", "B")], "dpi": 50}
        data = figures.FigureData(self.args, config)
        out = self.root / "zero"
        out.mkdir()
        figures.render_local(data, out)
        self.assertEqual(data.edges, [])
        self.assertEqual(len(data.drawn_genes), 2)
        self.assertTrue(all(g["clipped_to_window"] for g in data.drawn_genes.values()))


if __name__ == "__main__":
    unittest.main()
