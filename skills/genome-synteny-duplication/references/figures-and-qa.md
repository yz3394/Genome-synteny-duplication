# Figures and visual verification

Inference remains MCScanX. Rendering with JCVI does not turn the result into JCVI's independent MCscan inference. Use actual external tools when those backends are selected; missing tools are errors, not permission to silently substitute a different renderer.

## Choose a scale

| Question | Backend / `kind` | What is drawn |
|---|---|---|
| Whole-genome internal duplication | Native Circos / `circle` | Chromosomes at assembly lengths; pale block-span ribbons, focal direct anchors, an outer tandem track |
| Two-genome chromosome correspondence | Native JCVI karyotype / `parallel` | Two chromosome rows, block-span ribbons, focal anchors and individual genome scale bars |
| Rearrangement / many-to-many pattern | Matplotlib / `dotplot` | Actual direct anchor coordinates, with original chromosome boundaries |
| Local gene order and tandem candidates | Matplotlib / `local` | Physical gene spans, strand arrows, overlap lanes, actual cross-track anchors, independent within-track tandem brackets |

The local renderer is explicitly a custom evidence layer: it supports long overlapping models and simultaneous independent tandem evidence that the macro backends do not directly encode. Every backend exports the evidence actually drawn. Native Circos produces SVG + PNG; JCVI and the other backends produce PDF + SVG + PNG. Do not describe a bitmap embedded in a PDF as vector output.

```bash
python scripts/render_figures.py \
  --genes PREPARED_A/genes.tsv PREPARED_B/genes.tsv \
  --chromosomes PREPARED_A/chromosomes.tsv PREPARED_B/chromosomes.tsv \
  --results CROSS_AB/results \
  --families OVERLAY/family_members.tsv \
  --config plot.json --output NEW_FIGURE_DIRECTORY
```

Use a Python environment with Matplotlib. Native JCVI can use a different configured Python. `--families` is optional; supplied tables require `genome_id`, `gene_id`, and `families` (semicolon separated) or `family_label`. For protein-only family lists, first run the overlay mapping; the plotter consumes gene IDs. Output directories must be new.

In a restricted environment, point `MPLCONFIGDIR` to a writable project/temporary directory if the default font cache is unavailable.

Use cross-comparison results for the parallel view, and a self-comparison for a circle. For local views, additionally supply self-comparison result directories to draw tandem evidence. A whole-block overview can set `label_mode: "none"`; retain readable gene IDs in its paired close-up and all IDs in the exported table.

Minimal macro configuration:

```json
{
  "kind": "parallel",
  "genomes": ["A", "B"],
  "title": "Chromosome correspondence",
  "focus_families": ["Family1", "Family2"],
  "palette": {"Family1": "#D55E00", "Family2": "#0072B2"},
  "min_block_anchors": 20,
  "min_family_block_anchors": 0,
  "tools": {"jcvi_python": "/path/to/jcvi/env/bin/python"}
}
```

For `circle`, use one genome and `tools.circos`. Restrict/order chromosomes with `chromosomes: {"A": ["chr1", "chr2"]}`; optional `chromosome_labels` and `labels` change displayed text while source IDs remain in exports. Minimum anchor counts here are **display filters**, never a replacement for statistical detection thresholds. Focal anchors can remain visible in smaller blocks, independently of macro background filtering.

Minimal local configuration:

```json
{
  "kind": "local",
  "title": "Local gene order",
  "tracks": [
    {"genome_id": "A", "seqid": "chr1", "start": 100000, "end": 300000, "reverse": true, "label": "Species A"},
    {"genome_id": "B", "seqid": "chr3", "start": 400000, "end": 600000, "reverse": false, "label": "Species B"}
  ],
  "focus_families": ["Family1", "Family2"],
  "tandem_scope": "focal",
  "figsize": [14, 8]
}
```

Windows are 1-based closed coordinates. Tracks are independently scaled and labelled as such. Only genes physically overlapping the window are drawn; clipped spans are marked. Arrow direction reverses with the axis. A gene may overlap a window while its midpoint lies outside it: its span can be visible without drawing an off-window anchor endpoint. Use a wider overview panel to show block context outside a focal close-up.

## Visual grammar

- Gray gene marks mean other annotated loci. Gray dashed lines mean **observed direct anchor pairs**, not merely two nearby genes or block membership. If no gray anchors occur inside the chosen windows, draw none and state that explicitly.
- Pale macro ribbons summarize the envelope of a block's anchors. They do not assert aligned nucleotide sequence or that every intervening gene is paired.
- Family-colored links require a shared displayed family at their two anchor endpoints. Multiple family labels are retained in the tables; a declared `family_priority` chooses a single display color when needed.
- Purple same-track brackets denote independent strict tandem candidates. Dashed purple brackets mark review candidates. Physical proximity alone is never sufficient to add a bracket.
- Tandem review flags accompany the exported tandem edges. Anchor edge tables represent chain evidence; obtain their sequence coverage and model-review evidence by joining the corresponding `homology_pairs.tsv` on genome/gene pair. An anchor line does not certify full-length homology or model integrity.
- Family identity, homology, adjacency, collinearity and demonstrated enzyme activity remain separate claims. A missing anchor is not evidence of gene loss.

## Required QA and delivery

Inspect the generated PNG and vector file, including labels at the intended final publication width. Check reverse axes, endpoints, biological metadata, scale, overlap lanes, boundary clipping, foreground/background ordering, distinguishable colors and a complete legend. Dense whole-genome family overlays may need separate family panels; do not shrink labels until unreadable.

Reconcile `figure_edges.tsv` and `figure_genes.tsv` with source evidence and the visible counts. Save `omitted_counts.tsv`, `caption.txt`, native inputs/commands and `manifest.json`. The manifest begins with `rendered_requires_visual_review`; write a separate dated QA record of what was inspected rather than silently claiming rendering proves publication readiness. Preserve failed/earlier outputs and create a new output directory after a correction.
