#!/usr/bin/env python3
"""Render source-backed synteny figures without inferring new relationships.

Input coordinates are 1-based inclusive. Native BED/Circos inputs are converted
explicitly. This adapter never reruns inference or selects one best ortholog.
Run --help for the command line; plot contracts are in references/figures-and-qa.md.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#332288", "#44AA99", "#882255"]
GRAY = "#b6bdc4"
TANDEM = "#7b3294"
EDGE_FIELDS = ["figure_edge_id", "relation_type", "comparison_id", "block_ids", "genome_a", "gene_a", "seqid_a", "start_a", "end_a", "genome_b", "gene_b", "seqid_b", "start_b", "end_b", "shared_families", "display_family", "color", "quality_tier", "review_flags", "source_table", "source_row", "array_id"]


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def position(track, bp):
    """Physical coordinate in [0,1]; never invert the source coordinates."""
    span = track["end"] - track["start"]
    if span <= 0:
        raise ValueError("Plot intervals must span at least two bases")
    value = (bp - track["start"]) / span
    return 1 - value if track.get("reverse", False) else value


def visible_strand(track, strand):
    if strand not in ("+", "-"):
        return "."
    return ({"+": "-", "-": "+"}[strand] if track.get("reverse") else strand)


def block_orientation(value):
    if value in ("-", "minus", "reverse"):
        return "-"
    if value in ("+", "plus", "forward"):
        return "+"
    raise ValueError(f"Unknown block orientation: {value}")


def assign_lanes(genes):
    """Greedy interval coloring keeps overlapping gene models on separate lanes."""
    ends, lanes = [], {}
    for g in sorted(genes, key=lambda g: (g["start"], g["end"], g["gene_id"])):
        lane = next((i for i, end in enumerate(ends) if end < g["start"]), len(ends))
        if lane == len(ends):
            ends.append(g["end"])
        else:
            ends[lane] = g["end"]
        lanes[(g["genome_id"], g["gene_id"])] = lane
    return lanes


def load_inputs(args):
    genes, chromosomes, sources = {}, {}, []
    for path in args.genes:
        sources.append(Path(path))
        for row in read_tsv(path):
            key = row["genome_id"], row["gene_id"]
            for field in ("start", "end"):
                row[field] = int(row[field])
            if key in genes and genes[key] != row:
                raise ValueError(f"Conflicting gene record: {key}")
            genes[key] = row
    for path in args.chromosomes or []:
        sources.append(Path(path))
        for row in read_tsv(path):
            key = row["genome_id"], row["seqid"]
            row["length"] = int(row["length"])
            if row["length"] <= 0:
                raise ValueError(f"Nonpositive sequence length: {key}")
            if key in chromosomes and chromosomes[key] != row:
                raise ValueError(f"Conflicting chromosome length: {key}")
            chromosomes[key] = row
    families = defaultdict(set)
    if args.families:
        sources.append(Path(args.families))
        for row in read_tsv(args.families):
            key = row["genome_id"], row["gene_id"]
            if key not in genes:
                raise ValueError(f"Family gene missing from normalized input: {key}")
            label = row.get("family_label", row.get("family", row.get("families", "")))
            families[key].update(x for x in label.split(";") if x)
    blocks, anchors, tandems = [], [], []
    for directory in args.results:
        for filename, dest in (("blocks.tsv", blocks), ("anchor_pairs.tsv", anchors), ("tandem_pairs.tsv", tandems)):
            path = Path(directory) / filename
            if not path.exists():
                if filename == "tandem_pairs.tsv":
                    continue
                raise FileNotFoundError(path)
            sources.append(path)
            for number, row in enumerate(read_tsv(path), 2):
                row["source_table"], row["source_row"] = str(path.resolve()), number
                dest.append(row)
    seen_blocks = set()
    for row in blocks:
        key = row["comparison_id"], row["block_id"]
        if key in seen_blocks:
            raise ValueError(f"Duplicate block key: {key}")
        seen_blocks.add(key)
        row["anchor_count"] = int(row["anchor_count"])
    for row in anchors + tandems:
        for side in ("a", "b"):
            genome = row.get(f"genome_{side}") or row.get("genome_id")
            key = genome, row[f"gene_{side}"]
            if key not in genes:
                raise ValueError(f"Relationship endpoint missing: {key}")
            row[f"genome_{side}"] = genome
            for field in ("seqid", "start", "end"):
                value = row.get(f"{field}_{side}")
                if value not in (None, "") and str(value) != str(genes[key][field]):
                    raise ValueError(f"Relationship coordinate disagrees with gene input: {key}, {field}")
    for row in anchors:
        if (row["comparison_id"], row["block_id"]) not in seen_blocks:
            raise ValueError("Anchor references an unknown block")
    return genes, chromosomes, families, blocks, anchors, tandems, sources


class FigureData:
    def __init__(self, args, config):
        paths = list(args.genes) + list(args.chromosomes or []) + [Path(__file__)]
        if args.families:
            paths.append(args.families)
        if Path(args.config).exists():
            paths.append(args.config)
        paths += [Path(directory) / name for directory in args.results for name in ("blocks.tsv", "anchor_pairs.tsv", "tandem_pairs.tsv") if (Path(directory) / name).exists()]
        self.initial_hashes = {str(Path(p).resolve()): digest(p) for p in paths}
        self.genes, self.chromosomes, self.families, self.blocks, self.anchors, self.tandems, self.sources = load_inputs(args)
        self.verify_sources_unchanged()
        self.config = config
        self.edges, self.drawn_genes, self.omitted = [], {}, Counter()
        labels = sorted(set().union(*self.families.values())) if self.families else []
        self.focus = set(config.get("focus_families", labels))
        labels = [label for label in labels if label in self.focus]
        self.colors = {label: PALETTE[i % len(PALETTE)] for i, label in enumerate(labels)}
        self.colors.update({k: v for k, v in config.get("palette", config.get("families_palette", {})).items() if k in labels})
        self.priority = config.get("family_priority", labels)
        self.block_map = {(b["comparison_id"], b["block_id"]): b for b in self.blocks}
        self.block_ids = set(config.get("block_ids", []))
        if self.block_ids:
            available = {b["block_id"] for b in self.blocks}
            if self.block_ids - available:
                raise ValueError(f"Unknown requested blocks: {self.block_ids - available}")

    def family_style(self, a, b=None):
        labels = (self.families[a] if b is None else self.families[a] & self.families[b]) & self.focus
        ordered = [x for x in self.priority if x in labels] + sorted(labels - set(self.priority))
        family = ordered[0] if ordered else ""
        return ";".join(sorted(labels)), family, self.colors.get(family, GRAY)

    def verify_sources_unchanged(self):
        for path, expected in self.initial_hashes.items():
            if digest(path) != expected:
                raise RuntimeError(f"Source changed while plotting; use stable versioned inputs: {path}")

    def allowed_block(self, row, family=False):
        if self.block_ids and row["block_id"] not in self.block_ids:
            return False
        n = int(row.get("anchor_count", self.block_map[(row["comparison_id"], row["block_id"])]["anchor_count"]))
        minimum = int(self.config.get("min_family_block_anchors", 0) if family else self.config.get("min_block_anchors", 0))
        return n >= minimum

    def add_gene(self, key, role="gene", **extra):
        g = dict(self.genes[key])
        roles = set(filter(None, self.drawn_genes.get(key, {}).get("display_role", "").split(";")))
        roles.add(role)
        g.update(family_labels=";".join(sorted(self.families[key])), display_role=";".join(sorted(roles)), **extra)
        self.drawn_genes[key] = g

    def add_pair(self, row, kind="direct_anchor", **extra):
        a = row["genome_a"], row["gene_a"]
        b = row["genome_b"], row["gene_b"]
        shared, family, color = self.family_style(a, b)
        edge = dict(row)
        edge.update(relation_type=kind, block_ids=row.get("block_id", ""), shared_families=shared, display_family=family, color=color)
        edge.update(extra)
        for side, key in (("a", a), ("b", b)):
            for field in ("seqid", "start", "end"):
                edge[f"{field}_{side}"] = self.genes[key][field]
        self.edges.append(edge)
        return edge

    def unique_anchors(self):
        """Draw each gene pair once while retaining all supporting block IDs."""
        pairs = {}
        for row in self.anchors:
            a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
            shared = self.family_style(a, b)[0]
            if not self.allowed_block(row, bool(shared)):
                self.omitted["anchor_rows_excluded_by_block_filter"] += 1
                continue
            key = tuple(sorted((a, b)))
            if key in pairs:
                pairs[key]["block_id"] += ";" + row["block_id"]
                pairs[key]["source_row"] = str(pairs[key]["source_row"]) + ";" + str(row["source_row"])
                pairs[key]["source_table"] += ";" + row["source_table"]
                self.omitted["duplicate_anchor_support_rows_collapsed"] += 1
            else:
                pairs[key] = dict(row)
        return list(pairs.values())


def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none", "pdf.fonttype": 42, "text.usetex": False})
    return plt


def legend_handles(data, macro=False, tandem=True):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=GRAY, alpha=.35, label="Block span (not base-pair alignment)")] if macro else [Line2D([0], [0], color=GRAY, ls="--", lw=.8, label="Other direct anchor")]
    handles += [Line2D([0], [0], color=c, lw=2, label=f) for f, c in data.colors.items()]
    if tandem:
        handles.append(Line2D([0], [0], color=TANDEM, marker="|", label="Strict tandem candidate (same track)"))
    return handles


def save_figure(fig, out, config):
    for extension in ("png", "svg", "pdf"):
        fig.savefig(out / ("figure." + extension), dpi=int(config.get("dpi", 300)), facecolor="white", bbox_inches="tight", pad_inches=.15)


def render_local(data, out):
    plt = pyplot()
    from matplotlib.patches import FancyArrow, Rectangle
    tracks = data.config.get("tracks", [])
    if len(tracks) < 2:
        raise ValueError("Local plot requires at least two explicit tracks")
    tracks = [dict(t, start=int(t["start"]), end=int(t["end"])) for t in tracks]
    gene_track, lane_map, track_genes = {}, {}, []
    for i, t in enumerate(tracks):
        if t["end"] <= t["start"]:
            raise ValueError("Track end must exceed start")
        length = data.chromosomes.get((t["genome_id"], t["seqid"]), {}).get("length")
        if t["start"] < 1 or (length is not None and t["end"] > length):
            raise ValueError("Track outside assembly bounds")
        genes = [g for g in data.genes.values() if g["genome_id"] == t["genome_id"] and g["seqid"] == t["seqid"] and g["start"] <= t["end"] and g["end"] >= t["start"]]
        lanes = assign_lanes(genes)
        for g in genes:
            key = g["genome_id"], g["gene_id"]
            if key in gene_track:
                raise ValueError("A gene occurs in multiple local tracks; use separate figures for repeated regions")
            gene_track[key] = i
            lane_map[key] = lanes[key]
        track_genes.append(genes)
    gap = max(1.8, .15 * (max(lane_map.values(), default=0) + 1) + 1.4)
    ys = [(len(tracks) - 1 - i) * gap for i in range(len(tracks))]
    figsize = data.config.get("figsize", [12, 2 + 2.2 * len(tracks)])
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=.19, right=.96, bottom=.18, top=.87)
    ax.set_xlim(-.02, 1.02)
    ax.set_ylim(-.85, ys[0] + .95)
    ax.axis("off")
    endpoint = {}
    for key, i in gene_track.items():
        g, t = data.genes[key], tracks[i]
        midpoint = (g["start"] + g["end"]) / 2
        endpoint[key] = (position(t, midpoint), ys[i] + .13 * lane_map[key])
    # Cross-track links are drawn only when both midpoint endpoints are visible.
    for row in data.unique_anchors():
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        if a not in endpoint or b not in endpoint or gene_track[a] == gene_track[b] or not (0 <= endpoint[a][0] <= 1 and 0 <= endpoint[b][0] <= 1):
            data.omitted["anchors_outside_local_endpoints"] += 1
            continue
        e = data.add_pair(row)
        family = bool(e["shared_families"])
        ax.plot([endpoint[a][0], endpoint[b][0]], [endpoint[a][1], endpoint[b][1]], color=e["color"], lw=1.5 if family else .7, alpha=.9 if family else .6, ls="-" if family else (0, (3, 3)), zorder=1)
    drawn_tandem = set()
    tandem_scope = data.config.get("tandem_scope", "focal" if data.config.get("focus_families") else "all")
    if tandem_scope not in ("all", "focal"):
        raise ValueError("tandem_scope must be all or focal")
    for row in data.tandems:
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        key = tuple(sorted((a, b)))
        if int(row.get("intervening_genes", 0)) != 0:
            raise ValueError("Nonadjacent pair in strict tandem table")
        if tandem_scope == "focal" and not data.family_style(a, b)[0]:
            data.omitted["tandem_pairs_outside_focal_families"] += 1
            continue
        if a not in endpoint or b not in endpoint or gene_track[a] != gene_track[b] or not (0 <= endpoint[a][0] <= 1 and 0 <= endpoint[b][0] <= 1):
            data.omitted["tandem_pairs_outside_local_endpoints"] += 1
            continue
        if key in drawn_tandem:
            continue
        drawn_tandem.add(key)
        x1, x2 = sorted((endpoint[a][0], endpoint[b][0]))
        y = ys[gene_track[a]] - .19
        review = row.get("quality_tier") == "review"
        ax.plot([x1, x1, x2, x2], [y + .05, y, y, y + .05], color=TANDEM, lw=1.35, ls="--" if review else "-", zorder=4)
        if review:
            ax.text((x1 + x2) / 2, y - .07, "review", ha="center", fontsize=7, color=TANDEM)
        data.add_pair(row, "strict_tandem", color=TANDEM)
    for i, (t, genes) in enumerate(zip(tracks, track_genes)):
        y = ys[i]
        ax.plot([0, 1], [y, y], color="#71808a", lw=.7, zorder=2)
        label = t.get("label", t["genome_id"])
        metadata = sorted({" / ".join(filter(None, [g.get("species"), g.get("assembly_id"), g.get("annotation_id"), g.get("subgenome")])) for g in genes})
        small = metadata[0] if len(metadata) == 1 else "see figure_genes.tsv"
        ax.text(-.035, y, label + "\n" + t["seqid"], ha="right", va="center", fontsize=9, fontweight="bold")
        ax.text(0, y + .85, small, fontsize=7.5, color="#526371", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": 1})
        ticks = [t["start"] + q * (t["end"] - t["start"]) / 4 for q in range(5)]
        for bp in ticks:
            x = position(t, bp)
            ax.plot([x, x], [y - .30, y - .34], lw=.5, color="#60717b")
            ax.text(x, y - .40, f"{bp / 1e6:.3f}", ha="center", va="top", fontsize=7.5)
        ax.text(1, y - .57, "Mb; " + ("reverse genomic order" if t.get("reverse") else "forward genomic order"), ha="right", fontsize=7.5, color="#526371")
        # Use non-overlapping label slots with leader lines; source span is unchanged.
        focal = sorted((g for g in genes if data.family_style((g["genome_id"], g["gene_id"]))[0]), key=lambda g: position(t, (g["start"] + g["end"]) / 2))
        label_mode = data.config.get("label_mode", "focal")
        if label_mode not in ("focal", "none"):
            raise ValueError("label_mode must be focal or none")
        if label_mode == "none":
            focal = []
        labels = {}
        if focal:
            minimum = min(.14, .92 / max(1, len(focal)))
            xs = [min(.97, max(.03, position(t, (g["start"] + g["end"]) / 2))) for g in focal]
            for j in range(1, len(xs)):
                xs[j] = max(xs[j], xs[j - 1] + minimum)
            if xs[-1] > .98:
                shift = xs[-1] - .98
                xs = [x - shift for x in xs]
            if xs[0] < .02:
                xs = [.02 + .96 * j / max(1, len(xs) - 1) for j in range(len(xs))]
            labels = {(g["genome_id"], g["gene_id"]): (x, j % 2) for j, (g, x) in enumerate(zip(focal, xs))}
        for g in genes:
            key = g["genome_id"], g["gene_id"]
            left, right = max(t["start"], g["start"]), min(t["end"], g["end"])
            x1, x2 = sorted((position(t, left), position(t, right)))
            gy = y + lane_map[key] * .13
            shared, family, color = data.family_style(key)
            strand = visible_strand(t, g.get("strand", "."))
            height = .075 if family else .05
            width = x2 - x1
            if strand == "." or width == 0:
                ax.add_patch(Rectangle((x1, gy - height / 2), width, height, facecolor=color, edgecolor="none", zorder=4))
            else:
                ax.add_patch(FancyArrow(x1 if strand == "+" else x2, gy, width if strand == "+" else -width, 0, width=height, head_width=height, head_length=min(.009, width * .35), length_includes_head=True, color=color, linewidth=0, zorder=4))
            clipped = g["start"] < t["start"] or g["end"] > t["end"]
            if clipped:
                for bp in (g["start"], g["end"]):
                    if not t["start"] <= bp <= t["end"]:
                        xx = min(1, max(0, position(t, bp)))
                        ax.plot([xx, xx], [gy - .05, gy + .05], color="#202020", lw=.8, zorder=5)
            if key in labels:
                lx, label_lane = labels[key]
                ly = y + .28 + .15 * label_lane + .13 * max(lane_map.values(), default=0)
                ax.annotate(g["gene_id"], xy=(min(1, max(0, endpoint[key][0])), gy + .05), xytext=(lx, ly), ha="center", va="bottom", fontsize=7.5, color=color, arrowprops={"arrowstyle": "-", "lw": .5, "color": color})
            data.add_gene(key, track_index=i, overlap_lane=lane_map[key], clipped_to_window=int(clipped), visible_strand=strand)
    fig.suptitle(data.config.get("title", "Local gene order and direct collinear anchors"), fontsize=13, y=.98)
    direct_count = sum(e["relation_type"] == "direct_anchor" for e in data.edges)
    gray_count = sum(e["relation_type"] == "direct_anchor" and not e["shared_families"] for e in data.edges)
    td_count = sum(e["relation_type"] == "strict_tandem" for e in data.edges)
    fig.text(.5, .94, f"{len(data.drawn_genes)} genes | {direct_count} direct anchors ({gray_count} gray) | {td_count} {tandem_scope} tandem pairs", ha="center", fontsize=9)
    fig.legend(handles=legend_handles(data), loc="lower center", ncol=min(4, len(data.colors) + 2), frameon=False, fontsize=8, bbox_to_anchor=(.5, .045))
    fig.text(.5, .012, "Physical gene spans; arrowheads show strand. Tracks are independently scaled. Purple brackets: adjacent homolog candidates.\nDashed purple brackets need review. Gray gene marks do not imply a matched gene; only drawn links are anchors.", ha="center", fontsize=7.5)
    save_figure(fig, out, data.config)
    plt.close(fig)
    return {"backend": "Matplotlib", "physical_scale": "independent per track", "tracks": tracks, "tandem_scope": tandem_scope}


def selected_chromosomes(data, genome):
    available = {seq: row for (gid, seq), row in data.chromosomes.items() if gid == genome}
    declared = data.config.get("chromosomes", {}).get(genome)
    if declared is None:
        declared = sorted(available)
    if not declared or any(seq not in available for seq in declared):
        raise ValueError(f"Missing actual chromosome lengths for {genome}")
    return [available[x] for x in declared]


def render_dotplot(data, out):
    plt = pyplot()
    genomes = data.config.get("genomes") or [t["genome_id"] for t in data.config.get("tracks", [])]
    if len(genomes) != 2:
        raise ValueError("dotplot requires genomes: [x_genome, y_genome]")
    offsets, totals, chroms = [], [], []
    for genome in genomes:
        rows = selected_chromosomes(data, genome)
        off, total = {}, 0
        for row in rows:
            off[row["seqid"]] = total
            total += row["length"]
        offsets.append(off)
        totals.append(total)
        chroms.append(rows)
    fig, ax = plt.subplots(figsize=data.config.get("figsize", [9, 8]))
    points = defaultdict(list)
    for row in data.unique_anchors():
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        if a[0] == genomes[1] and b[0] == genomes[0] and genomes[0] != genomes[1]:
            a, b = b, a
        ga, gb = data.genes[a], data.genes[b]
        if a[0] != genomes[0] or b[0] != genomes[1] or ga["seqid"] not in offsets[0] or gb["seqid"] not in offsets[1]:
            data.omitted["anchors_outside_selected_chromosomes"] += 1
            continue
        color = data.family_style(a, b)[2]
        x = (offsets[0][ga["seqid"]] + (ga["start"] + ga["end"]) / 2) / 1e6
        y = (offsets[1][gb["seqid"]] + (gb["start"] + gb["end"]) / 2) / 1e6
        points[color].append((x, y))
        data.add_pair(row)
        data.add_gene(a, "anchor_endpoint")
        data.add_gene(b, "anchor_endpoint")
    for color, values in sorted(points.items(), key=lambda item: item[0] != GRAY):
        xx, yy = zip(*values)
        ax.scatter(xx, yy, c=color, s=5, edgecolors="none", rasterized=False)
    for side, (rows, off) in enumerate(zip(chroms, offsets)):
        centers = [(off[r["seqid"]] + r["length"] / 2) / 1e6 for r in rows]
        for value in off.values():
            (ax.axvline if side == 0 else ax.axhline)(value / 1e6, color="#e1e4e7", lw=.5, zorder=0)
        labels = [data.config.get("chromosome_labels", {}).get(genomes[side], {}).get(r["seqid"], r["seqid"]) for r in rows]
        (ax.set_xticks if side == 0 else ax.set_yticks)(centers, labels=labels, rotation=90 if side == 0 else 0, fontsize=7)
    ax.set_xlim(0, totals[0] / 1e6)
    ax.set_ylim(0, totals[1] / 1e6)
    ax.set_xlabel(genomes[0] + " — chromosome lengths in physical scale (Mb)")
    ax.set_ylabel(genomes[1] + " — chromosome lengths in physical scale (Mb)")
    ax.secondary_xaxis("top").set_xlabel("Concatenated assembly position (Mb)")
    ax.secondary_yaxis("right").set_ylabel("Concatenated assembly position (Mb)")
    ax.set_title(data.config.get("title", "Direct collinear anchors"))
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], ls="", marker=".", color=GRAY, label="Other direct anchor")]
    handles += [Line2D([0], [0], ls="", marker=".", color=c, label=f) for f, c in data.colors.items()]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .08), ncol=min(5, len(handles)), frameon=False, fontsize=8)
    fig.text(.5, .04, f"{len(data.edges):,} observed gene pairs; no mirrored or inferred matches added.", ha="center", fontsize=8)
    metadata = sorted({" / ".join(filter(None, [g.get("species"), g.get("assembly_id"), g.get("annotation_id"), g.get("subgenome")])) for g in data.genes.values() if g["genome_id"] in genomes})
    fig.text(.5, .005, " | ".join(metadata), ha="center", fontsize=7)
    fig.subplots_adjust(bottom=.29, top=.85)
    save_figure(fig, out, data.config)
    plt.close(fig)
    return {"backend": "Matplotlib", "genomes": genomes, "physical_scale": "concatenated true sequence lengths"}


def macro_objects(data, genomes):
    selected = {(g, r["seqid"]): r for g in genomes for r in selected_chromosomes(data, g)}
    blocks, anchors = [], []
    for block in data.blocks:
        keys = [(block["genome_a"], block["seqid_a"]), (block["genome_b"], block["seqid_b"])]
        if not all(k in selected for k in keys) or not data.allowed_block(block):
            data.omitted["blocks_outside_selection_or_display_threshold"] += 1
            continue
        blocks.append(block)
        data.edges.append(dict(block, relation_type="block_span", block_ids=block["block_id"], color=GRAY))
    for row in data.unique_anchors():
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        ga, gb = data.genes[a], data.genes[b]
        if (a[0], ga["seqid"]) not in selected or (b[0], gb["seqid"]) not in selected:
            data.omitted["anchors_outside_selected_chromosomes"] += 1
            continue
        if not data.family_style(a, b)[0]:
            data.omitted["individual_background_anchors_represented_by_block_spans"] += 1
            continue
        anchors.append(row)
        data.add_pair(row)
        data.add_gene(a, "anchor_endpoint")
        data.add_gene(b, "anchor_endpoint")
    return selected, blocks, anchors


def run_logged(command, cwd, log, env=None):
    process = subprocess.run(command, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log.write_text(process.stdout)
    if process.returncode:
        raise RuntimeError(f"Renderer failed ({process.returncode}); see {log}")
    return {"argv": command, "cwd": str(cwd), "returncode": process.returncode}


def render_parallel(data, out):
    genomes = data.config.get("genomes") or [t["genome_id"] for t in data.config.get("tracks", [])]
    if len(genomes) != 2 or genomes[0] == genomes[1]:
        raise ValueError("JCVI parallel requires two distinct genomes")
    selected, blocks, anchors = macro_objects(data, genomes)
    native = out / "native_jcvi"
    native.mkdir()
    aliases = {key: f"T{genomes.index(key[0])}_{i + 1}" for i, key in enumerate(selected)}
    ids = {key: f"G{i:08d}" for i, key in enumerate(sorted(data.genes))}
    beds = defaultdict(list)
    for key, gene in data.genes.items():
        chrom = gene["genome_id"], gene["seqid"]
        if chrom in selected:
            beds[key[0]].append([aliases[chrom], gene["start"] - 1, gene["end"], ids[key], 0, gene.get("strand", ".")])
    simple = []
    for i, block in enumerate(blocks):
        ends = {}
        for side in ("a", "b"):
            genome, seq = block[f"genome_{side}"], block[f"seqid_{side}"]
            ends[genome] = []
            for which in ("start", "end"):
                name = f"BLOCK_{i}_{side}_{which}"
                bp = int(block[f"{which}_{side}"])
                beds[genome].append([aliases[genome, seq], bp - 1, bp, name, 0, "+"])
                ends[genome].append(name)
        if set(ends) != set(genomes):
            raise ValueError("Parallel figure includes within-genome block; provide relevant comparison dirs")
        simple.append("\t".join(ends[genomes[0]] + ends[genomes[1]] + [str(block["anchor_count"]), block_orientation(block.get("orientation", "+"))]))
    for row in anchors:
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        if a[0] == genomes[1]:
            a, b = b, a
        if a[0] != genomes[0] or b[0] != genomes[1]:
            raise ValueError("Parallel figure includes within-genome anchor")
        color = data.family_style(a, b)[2]
        simple.append(f"{color}*{ids[a]}\t{ids[a]}\t{ids[b]}\t{ids[b]}\t1\t+")
    seqids, layout = [], []
    sizes, display_labels = {}, {}
    for i, genome in enumerate(genomes):
        seq = []
        for key, row in selected.items():
            if key[0] != genome:
                continue
            alias = aliases[key]
            seq.append(alias)
            sizes[alias] = row["length"]
            display_labels[alias] = data.config.get("chromosome_labels", {}).get(genome, {}).get(row["seqid"], row["seqid"])
            beds[genome].append([alias, row["length"] - 1, row["length"], f"END_{alias}", 0, "+"])
        with (native / f"track{i}.bed").open("w") as f:
            for fields in sorted(beds[genome], key=lambda r: (r[0], r[1], r[2], r[3])):
                f.write("\t".join(map(str, fields)) + "\n")
        seqids.append(",".join(seq))
        label = data.config.get("labels", {}).get(genome, genome).replace(",", ";")
        layout.append(f"{.77 - i * .50}, .10, .96, 0, #536878, {label}, {'top' if i == 0 else 'bottom'}, track{i}.bed, center")
    layout.append("e, 0, 1, blocks.simple")
    (native / "seqids").write_text("\n".join(seqids) + "\n")
    (native / "layout").write_text("\n".join(layout) + "\n")
    (native / "blocks.simple").write_text("\n".join(simple) + "\n")
    write_tsv(native / "id_mapping.tsv", [{"render_id": rid, "genome_id": key[0], "gene_id": key[1]} for key, rid in ids.items()])
    spec = dict(data.config, colors=data.colors, sizes=sizes, chromosome_labels=display_labels, block_count=len(blocks), family_anchor_count=len(anchors), source_metadata=sorted({" / ".join(filter(None, [g.get("species"), g.get("assembly_id"), g.get("annotation_id"), g.get("subgenome")])) for g in data.genes.values() if g["genome_id"] in genomes}))
    spec_path = native / "render_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n")
    python = data.config.get("tools", {}).get("jcvi_python", sys.executable)
    cmd = [str(python), str(Path(__file__).resolve()), "--native-jcvi", str(spec_path.resolve()), str(out.resolve())]
    command = run_logged(cmd, native, out / "renderer.log")
    return {"backend": "JCVI karyotype + Matplotlib annotations", "command": command, "physical_scale": "true sequence lengths; independently scaled genome tracks", "native_version": json.loads((out / "renderer_version.json").read_text())}


def native_jcvi(spec_path, out):
    plt = pyplot()
    import importlib.metadata
    import jcvi.graphics.karyotype as karyotype
    from jcvi.graphics.base import normalize_axes
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    spec = json.loads(Path(spec_path).read_text())
    out = Path(out)
    old_shade, old_name = karyotype.Shade, karyotype.make_circle_name
    def shade(*args, **kwargs):
        highlighted = bool(kwargs.get("highlight"))
        kwargs["alpha"] = .85 if highlighted else .25
        kwargs["lw"] = 1.2 if highlighted else 0
        return old_shade(*args, **kwargs)
    karyotype.Shade = shade
    karyotype.make_circle_name = lambda sid, rev: spec["chromosome_labels"][sid]
    try:
        fig = plt.figure(figsize=spec.get("figsize", [13, 7]))
        ax = fig.add_axes((0, 0, 1, 1))
        plot = karyotype.Karyotype(ax, "seqids", "layout", generank=False, sizes=spec["sizes"], shadestyle="line", chrstyle="rect", plot_circles=True)
        normalize_axes(ax)
        import math
        for i, track in enumerate(plot.tracks):
            target = track.total / 10
            magnitude = 10 ** math.floor(math.log10(target))
            bp = max(n * magnitude for n in (1, 2, 5) if n * magnitude <= target)
            x, y = .16 + .46 * i, .17
            length = bp * track.ratio
            ax.plot([x, x + length], [y, y], color="#435466", lw=1)
            ax.plot([x, x], [y - .004, y + .004], color="#435466", lw=1)
            ax.plot([x + length, x + length], [y - .004, y + .004], color="#435466", lw=1)
            ax.text(x + length / 2, y + .013, f"{track.label}: {bp / 1e6:g} Mb", ha="center", fontsize=8)
        fig.text(.5, .96, spec.get("title", "Collinear block spans and family anchors"), ha="center", fontsize=13)
        fig.text(.5, .91, f"{spec['block_count']} block spans; {spec['family_anchor_count']} distinct family anchor pairs", ha="center", fontsize=9)
        handles = [Patch(facecolor=GRAY, alpha=.4, label="Block span (not base-pair alignment)")]
        handles += [Line2D([0], [0], color=c, lw=2, label=f) for f, c in spec["colors"].items()]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .07), ncol=min(5, len(handles)), frameon=False, fontsize=8)
        fig.text(.5, .028, "True sequence lengths; genome tracks independently scaled. MCScanX inferred; JCVI rendered.\n" + " | ".join(spec["source_metadata"]), ha="center", fontsize=7)
        save_figure(fig, out, spec)
        plt.close(fig)
    finally:
        karyotype.Shade, karyotype.make_circle_name = old_shade, old_name
    (out / "renderer_version.json").write_text(json.dumps({"jcvi": importlib.metadata.version("jcvi"), "matplotlib": importlib.metadata.version("matplotlib")}, indent=2) + "\n")


def render_circle(data, out):
    from matplotlib.colors import to_rgb
    genomes = data.config.get("genomes") or [t["genome_id"] for t in data.config.get("tracks", [])]
    if len(genomes) != 1:
        raise ValueError("Circos internal plot requires genomes: [genome_id]")
    selected, blocks, anchors = macro_objects(data, genomes)
    native = out / "native_circos"
    native.mkdir()
    aliases = {key: f"chr{i + 1}" for i, key in enumerate(selected)}
    karyotype = []
    for key, row in selected.items():
        label = data.config.get("chromosome_labels", {}).get(key[0], {}).get(row["seqid"], row["seqid"]).replace(" ", "_")
        karyotype.append(f"chr - {aliases[key]} {label} 0 {row['length']} chr_fill")
    (native / "karyotype.txt").write_text("\n".join(karyotype) + "\n")
    lines = []
    for block in blocks:
        parts = []
        for side in ("a", "b"):
            parts += [aliases[block[f"genome_{side}"], block[f"seqid_{side}"]], str(int(block[f"start_{side}"]) - 1), str(block[f"end_{side}"])]
        twist = ",twist=yes" if block_orientation(block.get("orientation", "+")) == "-" else ""
        lines.append(" ".join(parts) + " color=background_a3" + twist)
    (native / "blocks.links").write_text("\n".join(lines) + "\n")
    color_names = {family: f"fam{i}" for i, family in enumerate(data.colors)}
    lines = []
    for row in anchors:
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        parts = []
        for key in (a, b):
            gene = data.genes[key]
            bp = (gene["start"] + gene["end"]) // 2
            parts += [aliases[key[0], gene["seqid"]], str(bp - 1), str(bp)]
        family = data.family_style(a, b)[1]
        lines.append(" ".join(parts) + f" color={color_names[family]}")
    (native / "family.links").write_text("\n".join(lines) + "\n")
    highlights, seen = [], set()
    for row in data.tandems:
        a, b = (row["genome_a"], row["gene_a"]), (row["genome_b"], row["gene_b"])
        ga, gb = data.genes[a], data.genes[b]
        key = tuple(sorted((a, b)))
        if key in seen or (a[0], ga["seqid"]) not in selected or (b[0], gb["seqid"]) not in selected:
            data.omitted["tandem_pairs_outside_selection_or_duplicate"] += 1
            continue
        if a[0] != b[0] or ga["seqid"] != gb["seqid"] or int(row.get("intervening_genes", 0)) != 0:
            raise ValueError("Invalid strict tandem pair for circle")
        seen.add(key)
        highlights.append(f"{aliases[a[0], ga['seqid']]} {min(ga['start'], gb['start']) - 1} {max(ga['end'], gb['end'])}")
        data.add_pair(row, "strict_tandem_span", color=TANDEM)
        data.add_gene(a, "tandem_highlight_endpoint")
        data.add_gene(b, "tandem_highlight_endpoint")
    (native / "tandem.highlights").write_text("\n".join(highlights) + "\n")
    rgb = lambda c: ",".join(str(round(v * 255)) for v in to_rgb(c))
    colors = "\n".join(f"{color_names[f]} = {rgb(c)}" for f, c in data.colors.items())
    link_sections = []
    if blocks:
        link_sections.append("<link>\nfile = blocks.links\nribbon = yes\ncolor = background_a3\nz = 1\n</link>")
    if anchors:
        link_sections.append("<link>\nfile = family.links\nribbon = no\nthickness = 2p\nz = 5\n</link>")
    plots = "<plots>\n<plot>\ntype = highlight\nfile = tandem.highlights\nr0 = 1.03r\nr1 = 1.06r\nfill_color = tandem\nstroke_thickness = 0\n</plot>\n</plots>" if highlights else ""
    links_config = "\n".join(link_sections)
    config = f"""karyotype = karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = yes
<ideogram>
<spacing>
default = 0.015r
</spacing>
radius = 0.84r
thickness = 20p
fill = yes
fill_color = chr_fill
stroke_color = white
show_label = yes
label_font = default
label_radius = 1.12r
label_size = 25p
label_parallel = yes
</ideogram>
show_ticks = yes
show_tick_labels = yes
<ticks>
radius = 0.90r
color = chr_fill
thickness = 1p
multiplier = 1e-6
format = %d
<tick>
spacing = {max(1, int(max(r['length'] for r in selected.values()) / 5e6))}u
size = 8p
show_label = yes
label_size = 20p
label_offset = 3p
</tick>
</ticks>
<links>
radius = 0.80r
bezier_radius = 0.15r
{links_config}
</links>
{plots}
<image>
dir = .
file = native.png
png = yes
svg = yes
radius = 900p
angle_offset = -90
auto_alpha_colors = yes
auto_alpha_steps = 5
background = white
</image>
<colors>
<<include etc/colors.conf>>
chr_fill = 65,84,100
background = 178,188,198
tandem = {rgb(TANDEM)}
{colors}
</colors>
<fonts>
<<include etc/fonts.conf>>
</fonts>
<patterns>
<<include etc/patterns.conf>>
</patterns>
<<include etc/housekeeping.conf>>
"""
    (native / "circos.conf").write_text(config)
    executable = str(data.config.get("tools", {}).get("circos", "circos"))
    env = dict(os.environ)
    if Path(executable).is_absolute():
        env["PATH"] = str(Path(executable).parent) + os.pathsep + env.get("PATH", "")
    cmd = [executable, "-conf", "circos.conf"]
    command = run_logged(cmd, native, out / "renderer.log", env)
    version = subprocess.run([executable, "-version"], env=env, capture_output=True, text=True).stdout.strip()
    decorate_circle(data, native, out, len(blocks), len(anchors), len(highlights))
    return {"backend": "Circos; vector text legend added by adapter", "command": command, "native_version": version, "physical_scale": "true assembly chromosome lengths; ticks in Mb", "formats": ["svg", "png"], "pdf_note": "SVG is the native vector deliverable; no raster PDF substitution is produced."}


def decorate_circle(data, native, out, blocks, anchors, tandems):
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib import font_manager
    import textwrap
    svgpath = native / "native.svg"
    tree = ET.parse(svgpath)
    root = tree.getroot()
    ns = "http://www.w3.org/2000/svg"
    ET.register_namespace("", ns)
    image = Image.open(native / "native.png").convert("RGB")
    w, h = image.size
    extra, pad = 650, 90
    root.set("width", str(w + extra))
    root.set("height", str(h + pad))
    root.set("viewBox", f"0 {-pad} {w + extra} {h + pad}")
    bg = ET.Element(f"{{{ns}}}rect", {"x": "0", "y": str(-pad), "width": str(w + extra), "height": str(h + pad), "fill": "white"})
    root.insert(0, bg)
    canvas = Image.new("RGB", (w + extra, h + pad), "white")
    canvas.paste(image, (0, pad))
    draw = ImageDraw.Draw(canvas)
    fontpath = font_manager.findfont("DejaVu Sans")
    font = ImageFont.truetype(fontpath, 29)
    def label(x, y, text, color="#243746", size=29):
        element = ET.SubElement(root, f"{{{ns}}}text", {"x": str(x), "y": str(y), "fill": color, "font-family": "DejaVu Sans", "font-size": str(size)})
        element.text = text
        draw.text((x, y + pad - size), text, fill=color, font=ImageFont.truetype(fontpath, size))
    title = data.config.get("title", "Within-genome collinearity")
    label(40, -25, title, size=36)
    x, y = w + 15, 220
    lines = [(f"{blocks} block spans", GRAY), (f"{anchors} family anchor pairs", "#243746"), (f"{tandems} strict tandem spans", TANDEM)]
    lines += [(family, color) for family, color in data.colors.items()]
    lines += [("", "#243746"), ("Gray ribbons: block envelopes", "#526371"), ("(not base-pair alignment).", "#526371"), ("Colored lines: direct anchors.", "#526371"), ("Purple outer marks: adjacent", TANDEM), ("homolog candidates, including", TANDEM), ("flagged models; see edge table.", TANDEM), ("Chromosome ticks: Mb.", "#526371"), ("MCScanX inferred; Circos rendered.", "#526371")]
    metadata = sorted({" / ".join(filter(None, [g.get("species"), g.get("assembly_id"), g.get("annotation_id"), g.get("subgenome")])) for g in data.genes.values() if g["genome_id"] in data.config.get("genomes", [])})
    for value in metadata:
        lines += [(s, "#526371") for s in textwrap.wrap(value, 36)]
    for value, color in lines:
        for part in textwrap.wrap(value, 37) or [""]:
            label(x, y, part, color)
            y += 46
    if y > h - 40:
        raise ValueError("Circle legend too long; reduce families or split figure")
    tree.write(out / "figure.svg", encoding="utf-8", xml_declaration=True)
    canvas.save(out / "figure.png", dpi=(300, 300))


def finalize(data, out, args, details):
    import importlib.metadata
    data.verify_sources_unchanged()
    data.omitted["input_genes_without_a_drawn_glyph_or_endpoint"] = len(data.genes) - len(data.drawn_genes)
    for i, edge in enumerate(data.edges, 1):
        edge["figure_edge_id"] = f"E{i:07d}"
    write_tsv(out / "figure_edges.tsv", data.edges, EDGE_FIELDS)
    genes = list(data.drawn_genes.values())
    fields = list(dict.fromkeys(k for g in genes for k in g)) or ["genome_id", "gene_id", "display_role", "family_labels"]
    write_tsv(out / "figure_genes.tsv", genes, fields)
    write_tsv(out / "omitted_counts.tsv", [{"reason": k, "count": v} for k, v in sorted(data.omitted.items())], ["reason", "count"])
    counts = dict(Counter(e["relation_type"] for e in data.edges))
    caption = (data.config.get("title", "Synteny visualization") + ".\n" + f"Rendered with {details['backend']}. " + f"Displayed relations: {counts}; distinct displayed gene endpoints or glyphs: {len(genes)}. " + "Family colors indicate supplied family membership, not demonstrated biochemical activity. " + "Direct anchors are observed MCScanX gene pairs; block spans summarize the corresponding anchor envelopes and do not assert nucleotide alignment. " + "Strict tandem marks require adjacency in the declared annotation ranking and protein homology. Tandem review flags remain in figure_edges.tsv; inspect homology_pairs.tsv for anchor sequence coverage/model evidence. " + "Assembly labels, source identifiers, exact coordinates, evidence rows, software details, and display omissions accompany this figure. " + "Display thresholds are not statistical confidence thresholds.\n")
    (out / "caption.txt").write_text(caption)
    sources = list(dict.fromkeys(data.sources + [Path(args.config)]))
    manifest = {"status": "rendered_requires_visual_review", "config": data.config, "renderer": details, "matplotlib_version": importlib.metadata.version("matplotlib"), "adapter_sha256": digest(__file__), "source_object_counts": {"genes": len(data.genes), "blocks": len(data.blocks), "anchor_support_rows": len(data.anchors), "tandem_pairs": len(data.tandems)}, "displayed_relation_counts": counts, "displayed_gene_count": len(genes), "input_sha256": {str(p.resolve()): digest(p) for p in sources}, "command": sys.argv, "all_family_labels_retained_in_gene_table": True, "color_priority": data.priority, "visual_review_required": ["Check labels at intended publication dimensions", "Inspect overlap lanes, window clipping and reverse strand", "Reconcile relation table with the visible figure", "Review weak homology and gene model flags before biological claims"], "output_sha256": {str(p.relative_to(out)): digest(p) for p in sorted(out.rglob("*")) if p.is_file()}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--native-jcvi":
        native_jcvi(sys.argv[2], sys.argv[3])
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genes", nargs="+", required=True, type=Path)
    p.add_argument("--chromosomes", nargs="+", type=Path)
    p.add_argument("--results", nargs="+", required=True, type=Path)
    p.add_argument("--families", type=Path)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite figure directory: {args.output}")
    config = json.loads(args.config.read_text())
    if config.get("kind") not in ("dotplot", "local", "circle", "parallel"):
        raise ValueError("kind must be dotplot, local, circle, or parallel")
    data = FigureData(args, config)
    args.output.mkdir(parents=True)
    try:
        details = {"local": render_local, "dotplot": render_dotplot, "circle": render_circle, "parallel": render_parallel}[config["kind"]](data, args.output)
        finalize(data, args.output, args, details)
    except Exception as error:
        (args.output / "FAILED.json").write_text(json.dumps({"error": str(error), "config": config}, indent=2) + "\n")
        raise
    print(json.dumps({"output": str(args.output), "displayed_relations": len(data.edges), "genes": len(data.drawn_genes)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
