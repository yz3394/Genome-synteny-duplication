"""Scientific input invariants: loci, isoforms, identity, coordinates and reuse."""
import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from synteny_inputs import PreparationError, prepare_genome


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, name, text):
        path = self.root / name
        if name.endswith(".gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return str(path)

    def spec(self, gff, fasta, **extra):
        spec = dict(id="A", species="Species alpha", assembly_id="v1", annotation_id="anno1",
                    gff=self.write("input.gff3", gff), proteins=self.write("proteins.faa", fasta),
                    seq_lengths={"chr1": 1000})
        spec.update(extra)
        return spec

    def rows(self, out, name="genes.tsv"):
        with (out / name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    @staticmethod
    def gene(identifier, start, end, tx=None, kind="mRNA", attrs=""):
        lines = f"chr1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID={identifier}{attrs}\n"
        if tx:
            lines += f"chr1\tx\t{kind}\t{start}\t{end}\t.\t+\t.\tID={tx};Parent={identifier}\n"
        return lines

    def test_isoforms_one_representative_per_locus_and_identical_loci_retained(self):
        gff = self.gene("g1", 1, 100, "t1")
        gff += "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1\n"
        gff += "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t3;Parent=g1\n"
        gff += self.gene("g2", 300, 400, "t4")
        spec = self.spec(gff, ">t1\nMAA\n>t2\nMAAAAA\n>t3\nMAAAAA\n>t4\nMAAAAA\n")
        out = self.root / "out"
        manifest = prepare_genome(spec, out)
        rows = self.rows(out)
        self.assertEqual([r["protein_id"] for r in rows], ["t2", "t4"])
        self.assertEqual(manifest["counts"]["representative_proteins"], 2)
        self.assertEqual(len(self.rows(out, "all_isoforms.tsv")), 4)
        self.assertEqual(rows[0]["sequence_sha256"], rows[1]["sequence_sha256"])
        self.assertNotEqual(rows[0]["engine_id"], rows[1]["engine_id"])

    def test_missing_protein_locus_preserves_tandem_spacing(self):
        gff = self.gene("left", 10, 50, "p1")
        gff += self.gene("gap", 100, 150, "no_translation")
        gff += self.gene("right", 300, 400, "p3")
        spec = self.spec(gff, ">p1\nMAAA\n>p3\nMAAA\n")
        out = self.root / "out"
        result = prepare_genome(spec, out)
        rows = self.rows(out)
        self.assertEqual([r["annotation_rank"] for r in rows], ["1", "2", "3"])
        self.assertEqual([r["engine_rank"] for r in rows], ["1", "2", "3"])
        self.assertEqual([r["has_protein"] for r in rows], ["1", "0", "1"])
        self.assertEqual(len((out / "positions.gff").read_text().splitlines()), 3)
        self.assertEqual(result["counts"]["missing_protein_genes"], 1)

    def test_same_start_has_explicit_distinct_sort_rules(self):
        gff = self.gene("a", 100, 400, "pa") + self.gene("b", 100, 200, "pb")
        out = self.root / "out"
        manifest = prepare_genome(self.spec(gff, ">pa\nMAAA\n>pb\nMAAA\n"), out)
        rows = self.rows(out)
        self.assertEqual([r["gene_id"] for r in rows], ["b", "a"])
        self.assertEqual([r["engine_rank"] for r in rows], ["2", "1"])
        self.assertEqual(manifest["counts"]["rank_different_loci"], 2)

    def test_namespaces_prevent_cross_genome_original_id_collision(self):
        spec = self.spec(self.gene("same", 1, 100, "same.p"), ">same.p\nMAAA\n")
        prepare_genome(spec, self.root / "a")
        prepare_genome(dict(spec, id="B"), self.root / "b")
        a, b = self.rows(self.root / "a")[0], self.rows(self.root / "b")[0]
        self.assertEqual(a["gene_id"], b["gene_id"])
        self.assertNotEqual(a["engine_id"], b["engine_id"])
        self.assertNotEqual(a["engine_seqid"], b["engine_seqid"])

    def test_ncbi_protein_id_and_discontinuous_cds(self):
        gff = self.gene("gene-1", 1, 100, "rna-XM1")
        gff += "chr1\tx\tCDS\t1\t30\t.\t+\t0\tID=cds-XP1;Parent=rna-XM1;protein_id=XP1\n"
        gff += "chr1\tx\tCDS\t70\t100\t.\t+\t0\tID=cds-XP1;Parent=rna-XM1;protein_id=XP1\n"
        out = self.root / "out"
        prepare_genome(self.spec(gff, ">XP1 protein description\nMAAA*\n"), out)
        row = self.rows(out)[0]
        self.assertEqual(row["protein_id"], "XP1")
        self.assertEqual(row["protein_length"], "4")
        self.assertIn("terminal_stop_removed", row["model_flags"])
        self.assertEqual(self.rows(out, "all_isoforms.tsv")[0]["transcript_id"], "rna-XM1")

    def test_no_suffix_guess_and_explicit_id_map_resolves(self):
        spec = self.spec(self.gene("g1", 1, 100, "g1.T1"), ">g1.1\nMAAA\n")
        with self.assertRaisesRegex(PreparationError, "cannot be mapped"):
            prepare_genome(spec, self.root / "out")
        spec["id_map_tsv"] = self.write("map.tsv", "gene_id\ttranscript_id\tprotein_id\ng1\tg1.T1\tg1.1\n")
        prepare_genome(spec, self.root / "out")
        self.assertEqual(self.rows(self.root / "out")[0]["protein_id"], "g1.1")

    def test_ambiguous_and_contradictory_gene_mapping_rejected(self):
        gff = self.gene("g1", 1, 100, "p1") + self.gene("g2", 300, 400, "p2")
        spec = self.spec(gff, ">p1\nMAAA\n>p2\nMAAA\n")
        spec["id_map_tsv"] = self.write("map.tsv", "gene_id\ttranscript_id\tprotein_id\ng2\tp2\tp1\n")
        with self.assertRaisesRegex(PreparationError, "multiple loci"):
            prepare_genome(spec, self.root / "out")

    def test_frozen_representative_is_explicit_and_independent_of_families(self):
        gff = self.gene("g1", 1, 100, "short") + "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=long;Parent=g1\n"
        spec = self.spec(gff, ">short\nMAA\n>long\nMAAAAA\n")
        spec["representative_map_tsv"] = self.write("frozen.tsv", "gene_id\tprotein_id\ng1\tshort\n")
        result = prepare_genome(spec, self.root / "out")
        self.assertEqual(result["counts"]["frozen_shorter_than_longest"], 1)
        self.assertEqual(self.rows(self.root / "out")[0]["protein_id"], "short")
        with self.assertRaisesRegex(PreparationError, "family"):
            prepare_genome(dict(spec, family_members="irrelevant.tsv"), self.root / "illegal")

    def test_rank_policy_retains_missing_coding_and_excludes_known_noncoding(self):
        gff = self.gene("g1", 1, 100, "p1")
        gff += self.gene("coding_missing", 200, 250, attrs=";gene_biotype=protein_coding")
        gff += self.gene("nc", 300, 400, "nc.tx", kind="ncRNA", attrs=";gene_biotype=lncRNA")
        spec = self.spec(gff, ">p1\nMAAA\n", rank_policy="protein_coding")
        result = prepare_genome(spec, self.root / "out")
        self.assertEqual(result["counts"]["eligible_genes"], 2)
        self.assertEqual([r["gene_id"] for r in self.rows(self.root / "out")], ["g1", "coding_missing"])

    def test_default_retains_unknown_locus_between_homologs(self):
        gff = self.gene("left", 10, 50, "p1")
        gff += self.gene("unknown", 100, 150)
        gff += self.gene("right", 300, 400, "p3")
        gff += self.gene("known_nc", 500, 550, attrs=";gene_biotype=lncRNA")
        spec = self.spec(gff, ">p1\nMAAA\n>p3\nMAAA\n")
        result = prepare_genome(spec, self.root / "out")
        rows = self.rows(self.root / "out")
        self.assertEqual(result["rank_policy"], "protein_coding")
        self.assertEqual([row["gene_id"] for row in rows], ["left", "unknown", "right"])
        self.assertEqual(int(rows[2]["annotation_rank"]) - int(rows[0]["annotation_rank"]), 2)
        self.assertIn("unknown_biotype_retained", rows[1]["model_flags"])
        self.assertEqual(result["counts"]["unknown_biotype_retained_loci"], 1)
        self.assertEqual(result["counts"]["excluded_noncoding_loci"], 1)
        legacy = prepare_genome(dict(spec, rank_policy="all_genes"), self.root / "all_genes")
        self.assertEqual(legacy["counts"]["eligible_genes"], 4)

    def test_lengths_and_parent_containment_checked(self):
        spec = self.spec(self.gene("g1", 1, 100, "p1"), ">p1\nMAAA\n", seq_lengths={"chr1": 99})
        with self.assertRaisesRegex(PreparationError, "outside supplied assembly"):
            prepare_genome(spec, self.root / "out")
        bad = self.gene("g1", 1, 100) + "chr1\tx\tmRNA\t1\t101\t.\t+\t.\tID=p1;Parent=g1\n"
        spec = self.spec(bad, ">p1\nMAAA\n")
        with self.assertRaisesRegex(PreparationError, "Child outside parent"):
            prepare_genome(spec, self.root / "out")

    def test_mcscanx_signed_coordinate_boundary_and_overflow(self):
        maximum = 2_147_483_647
        spec = self.spec(self.gene("g1", maximum - 20, maximum, "p1"), ">p1\nMAAA\n",
                         seq_lengths={"chr1": maximum + 1000})
        prepare_genome(spec, self.root / "boundary")
        self.assertEqual(self.rows(self.root / "boundary")[0]["end"], str(maximum))
        # A valid assembly length can exceed the engine coordinate range; a
        # locus beyond the engine limit must never wrap to a negative position.
        spec = self.spec(self.gene("g1", maximum + 1, maximum + 100, "p1"), ">p1\nMAAA\n",
                         seq_lengths={"chr1": maximum + 1000})
        with self.assertRaisesRegex(PreparationError, "signed 32-bit limit.*2147483647"):
            prepare_genome(spec, self.root / "overflow")
        self.assertFalse((self.root / "overflow").exists())

    def test_missing_parent_and_cycles_are_errors(self):
        spec = self.spec(self.gene("g1", 1, 100) + "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=p1;Parent=absent\n", ">p1\nMAAA\n")
        with self.assertRaisesRegex(PreparationError, "does not exist"):
            prepare_genome(spec, self.root / "out")
        gff = "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1;Parent=p1\nchr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=p1;Parent=g1\n"
        with self.assertRaisesRegex(PreparationError, "Cycle"):
            prepare_genome(self.spec(gff, ">p1\nMAAA\n"), self.root / "out")

    def test_identical_reuse_and_changed_config_input_or_output_rejected(self):
        spec = self.spec(self.gene("g1", 1, 100, "p1"), ">p1\nMAAA\n")
        out = self.root / "out"
        manifest = prepare_genome(spec, out)
        before = {path.name: path.read_bytes() for path in out.iterdir()}
        self.assertEqual(prepare_genome(spec, out), manifest)
        with self.assertRaisesRegex(PreparationError, "Existing inputs differ"):
            prepare_genome(dict(spec, annotation_id="anno2"), out)
        self.assertEqual(before, {path.name: path.read_bytes() for path in out.iterdir()})
        Path(spec["proteins"]).write_text(">p1\nMAAAA\n")
        with self.assertRaisesRegex(PreparationError, "Existing inputs differ"):
            prepare_genome(spec, out)
        Path(spec["proteins"]).write_text(">p1\nMAAA\n")
        (out / "genes.tsv").write_text("changed")
        with self.assertRaisesRegex(PreparationError, "output is missing or changed"):
            prepare_genome(spec, out)

    def test_gzip_and_inferred_lengths_are_explicit(self):
        spec = self.spec(self.gene("g1", 1, 100, "p1"), ">p1\nMAAA\n")
        spec["gff"] = self.write("input.gff3.gz", Path(spec["gff"]).read_text())
        spec["proteins"] = self.write("proteins.faa.gz", ">p1\nMAAA\n")
        del spec["seq_lengths"]
        result = prepare_genome(spec, self.root / "out")
        self.assertEqual(result["length_source"], "inferred_not_assembly")
        self.assertTrue(any("incomplete physical" in warning for warning in result["warnings"]))
        self.assertEqual(self.rows(self.root / "out", "chromosomes.tsv")[0]["length"], "100")

    def test_internal_stop_duplicate_id_and_missing_gene_id_rejected(self):
        gff = self.gene("g1", 1, 100, "p1")
        for fasta, expected in [(">p1\nMA*AA\n", "internal stop"), (">p1\nMAA\n>p1\nMAA\n", "Duplicate FASTA")]:
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(PreparationError, expected):
                    prepare_genome(self.spec(gff, fasta), self.root / "out")
        with self.assertRaisesRegex(PreparationError, "lacks ID"):
            prepare_genome(self.spec("chr1\tx\tgene\t1\t100\t.\t+\t.\tName=g1\n", ">p1\nMAA\n"), self.root / "out")

    def test_normalized_table_multisoform_and_missing_protein(self):
        spec = dict(id="A", species="Species alpha", assembly_id="v1", annotation_id="anno1",
                    proteins=self.write("proteins.faa", ">p1\nMAA\n>p2\nMAAAAA\n"),
                    genes_tsv=self.write("genes.tsv", "gene_id\tseqid\tstart\tend\tstrand\tprotein_id\ng1\tchr1\t1\t100\t+\tp1\ng1\tchr1\t1\t100\t+\tp2\ng2\tchr1\t200\t250\t-\t\n"),
                    lengths_tsv=self.write("lengths.tsv", "seqid\tlength\nchr1\t1000\n"))
        result = prepare_genome(spec, self.root / "out")
        self.assertEqual(result["counts"]["eligible_genes"], 2)
        self.assertEqual(self.rows(self.root / "out")[0]["protein_id"], "p2")
        self.assertEqual(result["length_source"], "supplied_lengths")
        self.assertIn("unknown_biotype_retained", self.rows(self.root / "out")[1]["model_flags"])


if __name__ == "__main__":
    unittest.main()
