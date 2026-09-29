"""Opt-in real MCScanX smoke: SYNTENY_MCSCANX=/path/to/MCScanX.

On some platforms getopt stops at the first positional argument. Non-default
parameters and pairwise-only output must be demonstrated on the actual binary,
not inferred from its zero exit status or the supplied command line.
"""
import hashlib
import json
import math
import re
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


def exercise_engine(binary, out_dir, prefix_first=False):
    """Retain a tiny engine run; legacy prefix-first behavior is observational."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    (out_dir / "x.gff").write_text(
        "aa_chr1\tg1\t1\t100\naa_chr1\tg2\t200\t300\n"
        "bb_chr1\th1\t1\t100\nbb_chr1\th2\t200\t300\n")
    (out_dir / "x.blast").write_text(
        "g1\th1\t90\t30\t0\t0\t1\t30\t1\t30\t1e-30\t100\n"
        "g2\th2\t90\t30\t0\t0\t1\t30\t1\t30\t1e-30\t100\n")
    options = ["-s", "2", "-m", "7", "-a"]
    command = [str(Path(binary).resolve()), *( ["x"] + options if prefix_first else options + ["x"] )]
    result = subprocess.run(command, cwd=out_dir, capture_output=True, text=True)
    (out_dir / "engine.log").write_text(result.stdout + result.stderr)
    header = {}
    collinearity = out_dir / "x.collinearity"
    if collinearity.is_file():
        for line in collinearity.read_text().splitlines():
            if line.startswith("## Alignment"):
                break
            if line.startswith("# ") and ":" in line:
                key, value = line[2:].split(":", 1)
                try:
                    header[key.strip()] = float(value.strip())
                except ValueError:
                    pass
    record = dict(command=command, returncode=result.returncode, requested={"MATCH_SIZE": 2, "MAX GAPS": 7, "pairwise_only": True},
                  effective_header=header, html_files=[str(p.relative_to(out_dir)) for p in out_dir.rglob("*.html") if p.is_file()],
                  prefix_first=prefix_first,
                  input_sha256={name: hashlib.sha256((out_dir / name).read_bytes()).hexdigest() for name in ("x.gff", "x.blast")})
    record["nondefault_parameters_effective"] = header.get("MATCH_SIZE") == 2 and header.get("MAX GAPS") == 7
    record["pairwise_only_effective"] = not record["html_files"]
    (out_dir / "run.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def exercise_engine_thresholds(binary, out_dir):
    """Measure score and E-value gates on the actual executable.

    Two synthetic anchors have score 100 and E=4: their coordinate distances
    and region spans are all one, N=m=2, giving 2*P(2,2)=4. Coordinates are a
    numerical engine fixture, not plausible biological gene models. Adjacent
    floating-point E thresholds distinguish equality from a rounded display.
    """
    binary = Path(binary).resolve()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    cases = {}
    settings = [("score_equal", 2, 10.0), ("score_below_cutoff", 3, 10.0),
                ("e_just_below", 2, math.nextafter(4.0, -math.inf)),
                ("e_equal", 2, 4.0), ("e_just_above", 2, math.nextafter(4.0, math.inf))]
    for name, size, evalue in settings:
        folder = out_dir / name
        folder.mkdir()
        (folder / "x.gff").write_text(
            "aa_chr1\tg1\t1\t1\naa_chr1\tg2\t2\t2\n"
            "bb_chr1\th1\t1\t1\nbb_chr1\th2\t2\t2\n")
        (folder / "x.blast").write_text(
            "g1\th1\t90\t30\t0\t0\t1\t30\t1\t30\t1e-30\t100\n"
            "g2\th2\t90\t30\t0\t0\t1\t30\t1\t30\t1e-30\t100\n")
        command = [str(binary), "-k", "50", "-s", str(size), "-g", "0", "-m", "7", "-w", "0", "-e", repr(evalue), "-a", "x"]
        proc = subprocess.run(command, cwd=folder, capture_output=True, text=True)
        (folder / "engine.log").write_text(proc.stdout + proc.stderr)
        output = folder / "x.collinearity"
        text = output.read_text() if output.is_file() else ""
        alignments = [line for line in text.splitlines() if line.startswith("## Alignment")]
        records = []
        for line in alignments:
            match = re.search(r"score=(\S+) e_value=(\S+) N=(\d+)", line)
            if match:
                score, block_e, anchors = match.groups()
                records.append(dict(score=float(score), displayed_block_evalue=float(block_e), anchor_count=int(anchors)))
        cases[name] = dict(command=command, returncode=proc.returncode, expected_cutoff_score=50*size,
                           requested_evalue=evalue, requested_evalue_hex=evalue.hex(), blocks=records,
                           input_sha256={item: hashlib.sha256((folder / item).read_bytes()).hexdigest() for item in ("x.gff", "x.blast")},
                           collinearity_sha256=hashlib.sha256(output.read_bytes()).hexdigest() if output.is_file() else None)
        (folder / "run.json").write_text(json.dumps(cases[name], indent=2) + "\n")
    accepted = lambda name: len(cases[name]["blocks"]) == 1 and cases[name]["returncode"] == 0
    rejected = lambda name: not cases[name]["blocks"] and cases[name]["returncode"] == 0
    score_ok = (accepted("score_equal") and cases["score_equal"]["blocks"][0]["score"] == 100.0
                and rejected("score_below_cutoff"))
    e_ok = rejected("e_just_below") and rejected("e_equal") and accepted("e_just_above")
    record = dict(binary=str(binary), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  fixture=dict(anchors=2, match_score=50, gap_penalty=0, redundancy_window=0,
                               theoretical_chain_score=100, theoretical_block_evalue=4.0,
                               evalue_formula="2 * P(2,2) * (1*1)/(1*1) = 4"), cases=cases,
                  checks=dict(score_equal_cutoff_is_accepted=score_ok, block_evalue_requires_strictly_less_than_cutoff=e_ok),
                  status="PASS" if score_ok and e_ok else "FAIL",
                  limitation="Synthetic numerical boundary fixture; does not validate biological sensitivity or genome-wide default parameters.")
    source = binary.parent / "dagchainer.cc"
    if source.is_file():
        source_text = source.read_text()
        record["adjacent_source_evidence"] = dict(path=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            score_guard_present="path_score[i] >= CUTOFF_SCORE" in source_text,
            evalue_guard_present="return sf->e_value < E_VALUE;" in source_text,
            note="Source beside the executable corroborates observed gates; binary behavior is tested independently.")
    (out_dir / "threshold_summary.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


@unittest.skipUnless(os.environ.get("SYNTENY_MCSCANX"), "Set SYNTENY_MCSCANX for actual-engine parameter smoke")
class EngineSmokeTests(unittest.TestCase):
    def test_options_before_prefix_apply_and_pairwise_flag_suppresses_html(self):
        with tempfile.TemporaryDirectory() as temp:
            result = exercise_engine(os.environ["SYNTENY_MCSCANX"], Path(temp) / "run")
        self.assertEqual(result["returncode"], 0)
        self.assertTrue(result["nondefault_parameters_effective"], result)
        self.assertTrue(result["pairwise_only_effective"], result)


    def test_exact_chain_score_cutoff_is_inclusive(self):
        with tempfile.TemporaryDirectory() as temp:
            result = exercise_engine_thresholds(os.environ["SYNTENY_MCSCANX"], Path(temp) / "thresholds")
        self.assertTrue(result["checks"]["score_equal_cutoff_is_accepted"], result)

    def test_block_evalue_cutoff_is_strict_at_adjacent_floating_values(self):
        with tempfile.TemporaryDirectory() as temp:
            result = exercise_engine_thresholds(os.environ["SYNTENY_MCSCANX"], Path(temp) / "thresholds")
        self.assertTrue(result["checks"]["block_evalue_requires_strictly_less_than_cutoff"], result)


if __name__ == "__main__":
    unittest.main()
