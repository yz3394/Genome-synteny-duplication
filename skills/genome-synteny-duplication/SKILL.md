---
name: genome-synteny-duplication
description: 分析注释基因组内部串联与共线性重复、亚基因组或物种间共线性，并叠加关注家族及生成可追溯的科学图。Use for tandem duplication, WGD/segmental candidates, gene-level synteny and family-focused chromosome plots; consumes validated family lists rather than identifying enzyme function.
metadata:
  version: 1.0.0
---

# Genome synteny and duplication

Build a reusable genome-wide evidence set, then overlay arbitrary gene families without changing the search proteome. Preserve original gene/protein IDs and previous analysis versions.

## Choose the task

- **New genomes / new comparisons:** read [inputs and evidence](references/inputs-and-evidence.md) and [methods](references/methods-and-parameters.md). Prepare a config, audit its biological scope, then run `scripts/synteny_workflow.py build --config CONFIG.json`.
- **New family on an existing run:** map the supplied members to the existing gene/isoform map; run `scripts/synteny_workflow.py overlay --run RUN_MANIFEST.json --families FAMILIES.tsv`. Do not choose new analysis representatives from the family list.
- **New figure or locus review:** read [figures and QA](references/figures-and-qa.md). Use `scripts/render_figures.py` on existing normalized evidence. Recompute inference only if inputs/methods changed or a concrete sensitivity question requires it.
- **Legacy results:** adapt existing IDs and coordinates to the documented tables, preserve native outputs and parameters, and use `synteny_evidence.process_comparison` to import them. Mark unverified search completeness. Historical numbers are regression checks only, not targets for tuning a new analysis.

## Before computation

1. Identify assembly/annotation releases, species, ploidy, subgenomes and comparison scope. Distinguish two species, two subgenomes, alternative assemblies and phased haplotypes.
2. Use matched GFF3, proteins and sequence lengths/assembly. The helper supports GFF3 or an explicitly mapped gene table; do not silently treat arbitrary GTF or transcript IDs as gene IDs.
3. Freeze one representative per locus independently of family membership. A supplied canonical/curated representative map can override the deterministic longest-protein fallback. Keep the family evidence isoform separately.
4. Define the annotation rank universe and keep missing-protein loci in that universe. Inspect mapping failures, gene overlaps, ambiguous models and scaffold scope. Do not merge different loci because their proteins are identical.
5. Check available disk/CPU, executables and versions. The core helpers use Python standard libraries; plots require Matplotlib and the selected external backend. YAML is optional; JSON needs no extra parser. Resolve tool paths from the environment/config, not a previous project.

The example [configuration](assets/analysis_config.example.json) shows the schema; replace its paths and metadata. `prepare --config CONFIG.json` allows an input audit before the expensive searches.

## Run and interpret

- Search complete representative proteomes, not just focal families. Cross comparisons use both directions as a union, not an RBH or one-to-one orthology assertion.
- Retain a wide homology pool for tandem review and explicitly select per-query distinct nonself targets for MCScanX. Inspect target-cap saturation; run focused sensitivity checks when it could change the conclusion.
- Defaults are a documented starting profile, not universal biological thresholds. Keep sequence E-values separate from block E-values. Validate effective parameters in the output header.
- Put MCScanX options **before** the positional prefix: macOS `getopt` may ignore later options. Reject fractional integer parameters. The wrapper uses explicit comparison inputs and `-b 0`, avoiding two-character species inference.
- Keep independent tandem pairs/arrays alongside native mutually exclusive gene classes. Adjacency plus homology is candidate evidence; low coverage, short or overlapping models need review. Pair counts and array counts are not event counts.
- `block_member` is not `direct_anchor`. No anchor is not proof of gene loss; a shared family or a physical cluster is not proof of pathway activity.
- Report within-genome collinear duplicates as WGD/segmental candidates. Add Ks, synteny depth, outgroups or gene trees when the evolutionary question needs them; do not infer a specific WGD from a line plot alone.
- A successful empty result, unassessable input and an unexecuted step are different states. In pairwise-only mode native `.tandem` is not run; the independent tandem module has its own status.

## Deliver and verify

The workflow writes hashed `prepared/`, `searches/`, `comparisons/`, `runs/` and independent `families/` artifacts. Reuse requires matching inputs/parameters/code/tools and verified output hashes. Failed attempts are retained separately. Never accept an existing nonempty file as sufficient cache evidence.

Deliver the run-specific Methods, mapped IDs, blocks, anchors, all block genes, tandem pairs/arrays, review flags and family evidence. For figures, export their actual gene/edge tables, omitted counts, captions and vector/raster outputs. Inspect rendered figures at their intended publication size.

Write a project-level `ANALYSIS_RECORD.md` with input/run/figure paths, thresholds, software provenance, checks and unresolved limits. Update project state only under that project's existing authorization.

Run relevant invariant tests (`python -m unittest discover -s tests`) when modifying helpers or adapting unfamiliar input formats. The optional real-engine test uses `SYNTENY_MCSCANX=/absolute/path/MCScanX`. See [validated scope](references/coffee-regression.md) for measured tests and remaining limits. Do not claim a new dataset is validated merely because the packaged examples passed.
