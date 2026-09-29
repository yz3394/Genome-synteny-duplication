#!/usr/bin/env python3
"""Build family-independent synteny evidence, then add versioned family overlays."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

VERSION = "1.0.0"
HERE = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False).encode()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def table(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_config(path):
    path = Path(path).resolve()
    if path.suffix.lower() in (".yaml", ".yml"):
        import yaml
        cfg = yaml.safe_load(path.read_text())
    else:
        cfg = json.loads(path.read_text())
    for g in cfg.get("genomes", []):
        for k in ("gff", "genes_tsv", "proteins", "assembly_fasta", "lengths_tsv",
                  "id_map_tsv", "representative_map_tsv"):
            if g.get(k):
                g[k] = str((path.parent / g[k]).resolve())
    if cfg.get("output_dir"):
        cfg["output_dir"] = str((path.parent / cfg["output_dir"]).resolve())
    for k, v in cfg.get("tools", {}).items():
        if k in ("diamond", "blastp", "mcscanx", "classifier") and "/" in v:
            cfg["tools"][k] = str((path.parent / v).resolve())
    return cfg


def integer(value, label, minimum=None, maximum=None):
    if type(value) is not int:
        raise ValueError(f"{label} must be an integer, got {value!r}; no silent truncation")
    if minimum is not None and value < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    if maximum is not None and value > maximum:
        raise ValueError(f"{label} must be <= {maximum}")
    return value


def positive(value, label):
    if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0:
        raise ValueError(f"{label} must be positive and finite")
    return float(value)


def settings(cfg):
    search = dict(engine="diamond", sensitivity="sensitive", evalue=1e-5,
                  max_target_seqs=50, threads=4)
    allowed = set(search)
    if set(cfg.get("search", {})) - allowed:
        raise ValueError("Unknown search settings: " + str(set(cfg["search"]) - allowed))
    search.update(cfg.get("search", {}))
    if search["engine"] not in ("diamond", "blastp"):
        raise ValueError("search.engine must be diamond or blastp")
    if search["sensitivity"] not in ("sensitive", "very-sensitive", "ultra-sensitive"):
        raise ValueError("Unsupported DIAMOND sensitivity")
    search["evalue"] = positive(search["evalue"], "search.evalue")
    integer(search["max_target_seqs"], "search.max_target_seqs", 0)
    integer(search["threads"], "search.threads", 1)
    if search["engine"] == "blastp" and search["max_target_seqs"] == 0:
        raise ValueError("BLASTP max_target_seqs must be positive")
    col = dict(top_k=5, min_anchors=5, max_gaps=25, overlap_window=5,
               block_evalue=1e-5, match_score=50, gap_penalty=-1)
    if set(cfg.get("collinearity", {})) - set(col):
        raise ValueError("Unknown collinearity parameters (only explicit per-comparison -b 0 supported)")
    col.update(cfg.get("collinearity", {}))
    for key, low in [("top_k", 0), ("min_anchors", 2), ("max_gaps", 0),
                     ("overlap_window", 0), ("match_score", 1)]:
        integer(col[key], "collinearity." + key, low)
    integer(col["gap_penalty"], "collinearity.gap_penalty", maximum=0)
    col["block_evalue"] = positive(col["block_evalue"], "block_evalue")
    dup = dict(proximal_max_intervening=8, classifier_n=10)
    if set(cfg.get("duplication", {})) - set(dup):
        raise ValueError("Unknown duplication settings")
    dup.update(cfg.get("duplication", {}))
    integer(dup["proximal_max_intervening"], "proximal_max_intervening", 0)
    integer(dup["classifier_n"], "classifier_n", 2)
    return search, col, dup


def executable(value):
    candidate = shutil.which(value)
    if not candidate:
        raise FileNotFoundError(f"Executable unavailable: {value}")
    p = Path(candidate).resolve()
    return {"path": str(p), "sha256": digest(p)}


def run_command(cmd, log, cwd=None, accept_codes=(0,)):
    started = datetime.now(timezone.utc).isoformat()
    with Path(log).open("w") as f:
        f.write(json.dumps({"argv": [str(x) for x in cmd], "started_utc": started}) + "\n")
        f.flush()
        proc = subprocess.run([str(x) for x in cmd], cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    if proc.returncode not in accept_codes:
        raise RuntimeError(f"Command exited {proc.returncode}; inspect {log}")


def cache_step(parent, key_data, builder):
    """Reuse only completed, fingerprint-matching, output-hash-verified artifacts."""
    parent = Path(parent)
    parent.mkdir(parents=True, exist_ok=True)
    key = fingerprint(key_data)
    candidates = sorted(parent.glob(key + "*"))
    for directory in candidates:
        if not directory.is_dir():
            continue
        marker = directory / "COMPLETE.json"
        if not marker.exists():
            continue
        record = json.loads(marker.read_text())
        if record.get("fingerprint") != key or record.get("status") != "complete":
            continue
        for name, expected in record["outputs_sha256"].items():
            p = directory / name
            if not p.is_file() or digest(p) != expected:
                raise ValueError(f"Cached artifact changed: {p}; preserve it and use a new output root")
        return directory, True
    directory = parent / (key if not candidates else key + f"-retry{len(candidates)}")
    directory.mkdir(exist_ok=False)
    write_json(directory / "REQUEST.json", key_data)
    try:
        details = builder(directory) or {}
        outputs = {str(p.relative_to(directory)): digest(p) for p in sorted(directory.rglob("*"))
                   if p.is_file() and p.name not in ("COMPLETE.json", "FAILED.json")}
        write_json(directory / "COMPLETE.json", dict(status="complete", fingerprint=key,
                   completed_utc=datetime.now(timezone.utc).isoformat(), details=details,
                   outputs_sha256=outputs))
    except Exception as exc:
        write_json(directory / "FAILED.json", {"type": type(exc).__name__, "error": str(exc)})
        raise
    return directory, False


def verify_complete(directory):
    directory = Path(directory)
    record = json.loads((directory / "COMPLETE.json").read_text())
    if record.get("status") != "complete":
        raise ValueError(f"Incomplete artifact: {directory}")
    for name, expected in record["outputs_sha256"].items():
        p = directory / name
        if not p.is_file() or digest(p) != expected:
            raise ValueError(f"Artifact integrity failure: {p}")
    return record


def input_key(spec):
    files = {k: digest(spec[k]) for k in ("gff", "genes_tsv", "proteins", "assembly_fasta",
             "lengths_tsv", "id_map_tsv", "representative_map_tsv") if spec.get(k)}
    return {"spec": spec, "input_sha256": files, "input_code": digest(HERE / "synteny_inputs.py")}


def prepare(cfg):
    from synteny_inputs import prepare_genome
    out = Path(cfg["output_dir"])
    genomes = cfg.get("genomes", [])
    ids = [g["id"] for g in genomes]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("genomes must have unique nonempty IDs")
    prepared = {}
    for g in genomes:
        prepared[g["id"]], _ = cache_step(out / "prepared", input_key(g),
                                           lambda d, g=g: prepare_genome(g, d))
    return prepared


def search_pair(qdir, tdir, search, tool, root):
    key = {"query_sha256": digest(qdir / "proteins.faa"),
           "target_sha256": digest(tdir / "proteins.faa"), "settings": search,
           "engine": tool, "runtime_sha256": digest(__file__)}

    def build(d):
        engine = search["engine"]
        db = d / "db"
        hits = d / "hits.m8"
        if (qdir / "proteins.faa").stat().st_size == 0 or (tdir / "proteins.faa").stat().st_size == 0:
            hits.write_text("")
            return {"status": "not_assessable_no_proteins", "rows": 0}
        if engine == "diamond":
            run_command([tool["path"], "version"], d / "version.log")
            run_command([tool["path"], "makedb", "--in", tdir / "proteins.faa", "--db", db], d / "makedb.log")
            cmd = [tool["path"], "blastp", "--query", qdir / "proteins.faa", "--db", db,
                   "--" + search["sensitivity"], "--evalue", search["evalue"],
                   "--max-target-seqs", search["max_target_seqs"], "--max-hsps", "1",
                   "--outfmt", "6", "--threads", search["threads"], "--out", hits]
        else:
            run_command([tool["path"], "-version"], d / "version.log")
            make = tool["makeblastdb"]
            run_command([make["path"], "-in", tdir / "proteins.faa", "-dbtype", "prot", "-out", db], d / "makedb.log")
            cmd = [tool["path"], "-query", qdir / "proteins.faa", "-db", db,
                   "-evalue", search["evalue"], "-max_target_seqs", search["max_target_seqs"],
                   "-max_hsps", "1", "-outfmt", "6", "-num_threads", search["threads"], "-out", hits]
        run_command(cmd, d / "search.log")
        if not hits.is_file():
            raise RuntimeError("Successful search did not create requested hits file")
        return {"status": "complete" if hits.stat().st_size else "no_detected_hits"}

    return cache_step(root / "searches", key, build)[0] / "hits.m8"


def select_hits(raws, output, top_k, evalue, cap):
    """Streaming grouped-query m8 selection; top-k counts distinct nonself subjects."""
    summary = {"input_rows": 0, "selected_rows": 0, "queries_with_raw_hits": 0,
               "queries_at_search_cap": 0, "self_rows_removed": 0, "top_k_nonself": top_k,
               "search_cap": cap, "saturated_queries": []}
    seen = set()
    with Path(output).open("w") as dest:
        for raw in raws:
            current = None
            best = {}
            raw_subjects = set()

            def flush():
                if current is None:
                    return
                summary["queries_with_raw_hits"] += 1
                if cap and len(raw_subjects) >= cap:
                    summary["queries_at_search_cap"] += 1
                    summary["saturated_queries"].append(current)
                chosen = sorted(best.values(), key=lambda f: (-float(f[11]), float(f[10]), f[1]))
                if top_k:
                    chosen = chosen[:top_k]
                for fields in chosen:
                    dest.write("\t".join(fields) + "\n")
                    summary["selected_rows"] += 1

            with Path(raw).open() as h:
                for line in h:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) != 12:
                        raise ValueError(f"Expected 12 m8 columns in {raw}")
                    summary["input_rows"] += 1
                    if f[0] != current:
                        flush()
                        if f[0] in seen:
                            raise ValueError("Query rows are noncontiguous or raw directions overlap; group each query before selection")
                        seen.add(f[0])
                        current, best, raw_subjects = f[0], {}, set()
                    raw_subjects.add(f[1])
                    if f[0] == f[1]:
                        summary["self_rows_removed"] += 1
                        continue
                    if float(f[10]) > evalue:
                        continue
                    old = best.get(f[1])
                    if old is None or (-float(f[11]), float(f[10]), tuple(f)) < (-float(old[11]), float(old[10]), tuple(old)):
                        best[f[1]] = f
            flush()
    return summary


def mc_parameters(col):
    return ["-s", str(col["min_anchors"]), "-m", str(col["max_gaps"]),
            "-w", str(col["overlap_window"]), "-e", str(col["block_evalue"]),
            "-k", str(col["match_score"]), "-g", str(col["gap_penalty"])]


def verify_header(path, col):
    found = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("## Alignment"):
            break
        if line.startswith("# ") and ":" in line:
            k, v = line[2:].split(":", 1)
            try:
                found[k.strip()] = float(v.strip())
            except ValueError:
                pass
    expected = {"MATCH_SCORE": col["match_score"], "MATCH_SIZE": col["min_anchors"],
                "GAP_PENALTY": col["gap_penalty"], "OVERLAP_WINDOW": col["overlap_window"],
                "E_VALUE": float(format(col["block_evalue"], ".6g")), "MAX GAPS": col["max_gaps"]}
    for k, v in expected.items():
        if k not in found or not math.isclose(found[k], v, rel_tol=1e-6):
            raise ValueError(f"MCScanX effective parameter mismatch: {k}: requested {v}, got {found.get(k)}")
    return found


def comparison_run(comp, prepared, search, col, dup, toolmeta, root):
    from synteny_evidence import process_comparison
    a, b = comp["a"], comp["b"]
    if a not in prepared or b not in prepared:
        raise ValueError(f"Unknown comparison genome in {comp}")
    if comp["kind"] not in ("intragenome", "inter_subgenome", "interspecies", "inter_assembly"):
        raise ValueError("Declare comparison.kind")
    if (a == b) != (comp["kind"] == "intragenome"):
        raise ValueError("Self comparisons require kind=intragenome; cross comparisons require explicit other kind")
    raws = [search_pair(prepared[a], prepared[b], search, toolmeta["search"], root)]
    if a != b:
        raws.append(search_pair(prepared[b], prepared[a], search, toolmeta["search"], root))
    genes = table(prepared[a] / "genes.tsv")
    if a != b:
        genes += table(prepared[b] / "genes.tsv")
    if len({g["engine_id"] for g in genes}) != len(genes):
        raise ValueError("Cross-genome engine ID collision; change genome IDs before analysis")
    seqmap = {}
    for g in genes:
        original = (g["genome_id"], g["seqid"])
        if g["engine_seqid"] in seqmap and seqmap[g["engine_seqid"]] != original:
            raise ValueError("Cross-genome chromosome namespace collision")
        seqmap[g["engine_seqid"]] = original
    key = {"comparison": comp, "prepared": {g: digest(prepared[g] / "genes.tsv") for g in {a, b}},
           "raw_sha256": [digest(p) for p in raws], "search": search, "collinearity": col,
           "duplication": dup, "tools": toolmeta, "runtime_sha256": digest(__file__),
           "evidence_sha256": digest(HERE / "synteny_evidence.py")}

    def build(d):
        positions = "".join((prepared[g] / "positions.gff").read_text() for g in dict.fromkeys([a, b]))
        (d / "analysis.gff").write_text(positions)
        hit_qc = select_hits(raws, d / "analysis.blast", col["top_k"], search["evalue"], search["max_target_seqs"])
        write_json(d / "hit_selection.json", hit_qc)
        raw_manifest = [{"path": str(p), "sha256": digest(p)} for p in raws]
        write_json(d / "raw_searches.json", raw_manifest)
        classfile = None
        effective = None
        coli = d / "analysis.collinearity"
        if hit_qc["selected_rows"]:
            # macOS getopt stops at the first positional argument: options must precede prefix.
            run_command([toolmeta["mcscanx"]["path"], *mc_parameters(col), "-b", "0", "-a", "analysis"], d / "mcscanx.log", d)
            effective = verify_header(coli, col)
        else:
            coli.write_text("# No selected non-self hits; MCScanX not invoked.\n")
        if a == b and hit_qc["selected_rows"] and "classifier" in toolmeta:
            cd = d / "classifier"
            cd.mkdir()
            for ext in ("gff", "blast"):
                (cd / ("dup." + ext)).symlink_to(Path("..") / ("analysis." + ext))
            run_command([toolmeta["classifier"]["path"], *mc_parameters(col),
                         "-n", dup["classifier_n"], "dup"], cd / "classifier.log", cd)
            classfile = cd / "dup.gene_type"
            if not classfile.is_file():
                raise RuntimeError("Classifier completed without gene_type output")
        summary = process_comparison(genes, raws, coli, d / "results", comp,
                    search_evalue=search["evalue"], proximal_max_intervening=dup["proximal_max_intervening"],
                    classifier_path=classfile)
        summary["native_tandem_status"] = "not_run_pairwise_only_mode"
        summary["classifier_status"] = summary["statuses"]["classifier"]
        summary["mcscanx_status"] = "complete" if effective else "not_run_no_selected_hits"
        summary["effective_parameters"] = effective
        write_json(d / "comparison_manifest.json", summary)
        return summary

    directory, reused = cache_step(root / "comparisons", key, build)
    return {"id": comp["id"], "a": a, "b": b, "kind": comp["kind"],
            "directory": str(directory), "results": str(directory / "results"), "reused": reused}


def build(cfg):
    search, col, dup = settings(cfg)
    root = Path(cfg["output_dir"])
    tools = cfg.get("tools", {})
    toolmeta = {"search": executable(tools.get(search["engine"], search["engine"])),
                "mcscanx": executable(tools.get("mcscanx", "MCScanX"))}
    if search["engine"] == "blastp":
        toolmeta["search"]["makeblastdb"] = executable(str(Path(toolmeta["search"]["path"]).parent / "makeblastdb"))
    if tools.get("classifier"):
        toolmeta["classifier"] = executable(tools["classifier"])
    comps = cfg.get("comparisons", [])
    if not comps or len({x["id"] for x in comps}) != len(comps):
        raise ValueError("comparisons must have unique IDs")
    prepared = prepare(cfg)
    results = []
    for comp in comps:
        print(f"Building {comp['id']} ({comp['kind']})", file=sys.stderr, flush=True)
        results.append(comparison_run(comp, prepared, search, col, dup, toolmeta, root))
    # No families or plot settings in this fingerprint or any prepared/search/chain input.
    key = {"prepared": {g: str(p) for g, p in prepared.items()},
           "comparisons": [{k: v for k, v in r.items() if k != "reused"} for r in results],
           "search": search, "collinearity": col, "duplication": dup, "tools": toolmeta,
           "workflow_sha256": digest(__file__)}

    def write(d):
        record = {"version": VERSION, "python": platform.python_version(), **key,
                  "output_root": str(root), "status": "complete"}
        write_json(d / "run_manifest.json", record)
        (d / "METHODS.md").write_text(
            "# Run-specific methods\n\n"
            f"Search engine: {search['engine']}; protein E <= {search['evalue']}; "
            f"reported target cap {search['max_target_seqs']}; one HSP per target. "
            f"DIAMOND mode (when applicable): {search['sensitivity']}. "
            f"Per-query non-self distinct targets retained for chaining: {col['top_k']} "
            "(0 means all reported candidates). Reciprocal directions are a union, not RBH.\n\n"
            f"MCScanX options: {' '.join(mc_parameters(col))} -b 0 -a. "
            f"Score cutoff: {col['match_score'] * col['min_anchors']}. "
            "Input selection defines biological comparison scope. Native tandem output was not run; "
            "independent tandem candidates use adjacency in the declared annotation rank and the wider search pool. "
            "Low coverage/model flags require review; candidate pairs are not counts of proven duplication events.\n\n"
            "Each prepared genome records the family-independent representative rule and rank universe. "
            "No-hit and not-assessable statuses must be interpreted before biological absence claims. "
            "WGD/segmental is a combined algorithmic evidence label. All paths, binary hashes, source hashes, "
            "commands and output hashes are preserved in this run and referenced step manifests.\n")
        return {"comparisons": len(results)}
    run_dir, reused = cache_step(root / "runs", key, write)
    return {"run_manifest": str(run_dir / "run_manifest.json"), "reused": reused,
            "comparisons": results, "prepared": {k: str(v) for k, v in prepared.items()}}


def overlay(run_manifest, family_file):
    from synteny_evidence import overlay_families
    verify_complete(Path(run_manifest).parent)
    manifest = json.loads(Path(run_manifest).read_text())
    prepared = [Path(p) for p in manifest["prepared"].values()]
    for p in prepared:
        verify_complete(p)
    for c in manifest["comparisons"]:
        verify_complete(c["directory"])
    genes = [g for p in prepared for g in table(p / "genes.tsv")]
    isoforms = [g for p in prepared for g in table(p / "all_isoforms.tsv")]
    families = table(family_file)
    key = {"run_sha256": digest(run_manifest), "family_sha256": digest(family_file),
           "evidence_sha256": digest(HERE / "synteny_evidence.py")}
    d, reused = cache_step(Path(manifest["output_root"]) / "families", key,
        lambda d: overlay_families([Path(c["results"]) for c in manifest["comparisons"]],
                                   genes, families, d, isoforms=isoforms))
    return {"overlay_directory": str(d), "reused": reused}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("prepare", "build"):
        s = sub.add_parser(name)
        s.add_argument("--config", required=True)
    o = sub.add_parser("overlay")
    o.add_argument("--run", required=True)
    o.add_argument("--families", required=True)
    args = p.parse_args()
    if args.command == "overlay":
        result = overlay(args.run, args.families)
    else:
        cfg = read_config(args.config)
        settings(cfg)
        if args.command == "prepare":
            result = {k: str(v) for k, v in prepare(cfg).items()}
        else:
            result = build(cfg)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
