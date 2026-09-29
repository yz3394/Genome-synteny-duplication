# Methods and parameter decisions

## Baseline and alternatives

The helper's initial search is DIAMOND sensitive, protein E≤1e−5, cap50 distinct targets/query and one HSP/target. MCScanX receives the best five distinct nonself targets/query/target genome, unioned across reciprocal directions. Tandem review retains the wider reported pool. Search caps can affect DIAMOND's heuristic search; neither cap50 nor top5 proves exhaustive homology retrieval.

MCScanX settings: `-s 5 -m 25 -w 5 -e 1e-5 -k 50 -g -1 -b 0 -a`, then the positional prefix. Match size is a nominal lower bound implemented through score filtering: score cutoff=50×5=250, plus gap and block-E conditions. `-w` collapses local redundant hits; it is not a tandem distance. `-k` here is match score, not a search hit cap.

Other explicit profiles:

- `top_k=0` retains all reported passing candidates, supporting legacy comparison.
- `top_k=10` (and 20 only if needed) tests chain-input sensitivity with the same raw search pool.
- Protein E=1e−10 provides a stricter reference comparison. Do not change representatives, hit filtering and algorithm together and then attribute the difference to one factor.
- `search.engine=blastp` selects BLASTP with a version-adapted `-max_target_seqs` cap and `-max_hsps 1`. The protocol's biological target is top5 nonself, not literal top6 in every new dataset.
- For more divergent comparisons or unresolved focal loci, use DIAMOND very-sensitive or targeted BLASTP and document actual changes. A change of engine is a method comparison, not a guaranteed identical reproduction.

Strict tandem means zero intervening loci in the declared rank universe, with homology support. `proximal_max_intervening=8` means exactly that count. Native `classifier_n=10` is a different parameter: this pinned MCScanX checks rank difference<10, so its non-tandem proximal range is rank differences2–9. Native classification can replace a tandem gene's label with code4; independent pairs are retained.

Low best-HSP coverage (<0.5) and overlapping gene models carry review flags, not universal truth filters. Protein lengths are exported; evaluate unusually short sequences against the expected architecture of their family rather than automatically excluding every protein below a universal length cutoff. Review key uncertain candidates with full-length/domain/model evidence; report a candidate set separately from experimentally supported genes. Identity or domain membership alone does not establish enzyme activity.

## Checks that affect conclusions

Inspect raw target saturation before self-hit removal, unique-target retention after filtering, all eligible/no-protein/scaffold loci and ambiguous annotation. If saturated focal queries or missing adjacent homologs matter, increase the search cap or do focused searches; preserve both settings. Compare anchor sets, block split/merge patterns, focal relationships and review flags, not the aesthetic density of a plot.

Fingerprint includes source content, representative/rank policy, code, parameters and executable hashes; BLASTP also records makeblastdb. Changed settings create a new artifact; corrupted cached outputs stop reuse. Successful empty output is valid, while an unexecuted/native pairwise-only step is `not_run`. No source deletion is necessary to change parameters.

Species names are metadata. The native MCScanX `-b 1/2` compares the first two seqid characters, so the helper instead constructs explicit comparison inputs and uses `-b 0`. Do not enable native inference on arbitrary namespaces.

## Literature and verified implementation notes

- [MCScanX, Wang et al., NAR 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336/): genome-wide homology plus gene order, top5 reference, and priority-based classification.
- [Nature Protocols 2024](https://doi.org/10.1038/s41596-024-00968-2): pp2211/2215 specify reciprocal searches and E1e−10 or1e−5; Step26D overlays gene families on existing blocks. The supplied 24-page full text was reviewed during development.
- That protocol's p2216 default gap/E sentence conflicts with its p2213 help. Use the installed executable/source and recorded output header. Its p2211 suggestion of `-g -0.5` cannot be represented in the pinned integer-parsing implementation (`atoi`); the wrapper rejects it. This is an implementation check, not an assertion of a formal author correction.
- [Nature Protocols 2026 Addendum](https://doi.org/10.1038/s41596-026-01380-8): unfiltered isoforms can inflate tandem calls; representative selection is a documented heuristic.
- [Zhao et al., Nature Communications 2021](https://doi.org/10.1038/s41467-021-23665-0): plant microsynteny using DIAMOND and MCScanX; its `b5s5m25` shorthand is not a CLI `-b 5` instruction.
- [DIAMOND official options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options): search sensitivity/caps and HSP behavior are version-specific.

The development snapshot used MCScanX commit `0956cb8f900c152e2be8ba3196829e50ab454b94`. New installs must retain their own version/source provenance; do not assume the same defaults or parameter types.

Specific WGD, duplication timing, transposed duplication and one-to-one orthology need additional evidence appropriate to the question. Optional Ks requires matching CDS, reviewed codon alignments and a verified estimator; do not blindly call the obsolete BioPerl DNAStatistics wrapper. Outgroup order and species phylogeny matter for transposition inference.
