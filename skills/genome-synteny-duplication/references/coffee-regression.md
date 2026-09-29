# Validated development scope

Development validation on 2026-09-28 used Arabica ET-39 Coffee Genome Hub CC/CE and Catharanthus roseus v3. This is a worked regression case, not hardcoded expected behavior for new organisms. Project validation artifacts remain outside this portable skill.

## Full inputs and real execution

- CC: 33,585 loci; 33,571 available representative proteins; 14 missing-protein loci retained in rank.
- CE: 35,394 loci; 35,390 representatives; 4 missing-protein loci retained.
- Catharanthus: 66,262 input proteins reduced to 26,347 locus representatives; four loci have differing annotation/engine ranks across three same-start groups.
- The final input helper regenerated identical gene, protein, position, isoform and chromosome tables for all three explicit `all_genes` inputs after coordinate/unknown-biotype safeguards were added.
- Fresh complete-proteome searches and five actual MCScanX comparisons completed. Top5 baseline blocks: CC internal134; CE internal149; CC–CE432; CC–CRO495; CE–CRO519. These are parameter-dependent detections, not counts of independent evolutionary events.
- Independent tandem evidence recovered 2,116 CC pairs/1,497 arrays and 1,901 CE pairs/1,430 arrays; 67 pairs share at least one of the eight supplied focal-family labels. At least one-endpoint family filtering returns a different count and must be labelled accordingly.

## Historical evidence and figures

Imports of all five historical comparisons matched native block extents/counts/orientations, unique anchor pairs, all rank-interval member/anchor flags, and internal tandem gene-pair sets. All 4,017 tandem pairs matched best-HSP bitscore and bilateral coverage within historical rounding precision. Of these, 672 CC and624 CE pairs have coverage below0.5 on at least one side and remain review candidates.

The historical Chr3 figure regression recovered191 overview anchors (187 gray,4 colored),44 local overlapping loci (29 gray,15 focal), no gray anchors inside the close-up, and3 focal tandem edges in2 arrays. Native Circos, native JCVI, dotplot and physical-span local/overview renderers were executed and visually reviewed. Figure data carry source rows, omission counts and hashes. Example figures use historical results and are explicitly not new-top5 figure products.

## Controlled sensitivity and independent use

Top5/10/all-reported-cap50 chaining was compared with identical representatives and raw homology pools in all five comparisons. For CC internal, blocks were134/231/2063; for CE149/199/2167. Cross-species block totals were1014/1038/1147. Some anchors were lost as well as gained when increasing the pool; the outputs are not nested confidence tiers. Tandem evidence uses the unchanged broader pool. Whole-proteome raw-query saturation is recorded; cap50 is not exhaustive.

A separate representative-map comparison holds genome namespace, coordinates/ranks, engine and chain policy fixed. Legacy vs canonical selection can change protein IDs without changing sequences; reports separate these cases from true sequence differences and distinguish focal-family priority from equal-length isoform tie-breaking.

That full-genome controlled comparison found3,048 representative-ID changes but only140 amino-acid sequence changes (61 shorter focal representatives and79 equal-length nonfocal sequence choices). With top5, the frozen legacy map lost2/added1 CC–CRO anchors and lost9 CE–CRO anchors relative to canonical longest. With all reported cap50 hits, both cross-species unique anchor sets were identical. The four Chr3 MATE/TDC direct anchors were retained under both representative maps and all tested top-k profiles. This measures parameter interaction, not universal insensitivity to isoform choice.

Six actual regional profiles used322 Coffee and236 Catharanthus loci with all original isoforms available: DIAMOND sensitive E1e−5/E1e−10, BLASTP E1e−5/E1e−10, and separate DIAMOND cap200 / very-sensitive comparisons. Cross-block counts were2/1/2/1/2/2, unique anchors108/98/112/99/108/108. Both focal MATE/TDC anchors and the CC TDC tandem pair were retained in all six profiles. Regional database size and scope matter; this does not establish full-proteome BLASTP equivalence or cap completeness.

An independent executor received only the skill and raw regional inputs, built self/cross evidence, rendered local figures and added a new protein-ID gene set. All118 core artifact files retained identical hashes and modification times. A nonrepresentative source isoform mapped to its unchanged analysis representative. Two reference/filename mistakes and overly broad caption wording found by this test were corrected.

## Boundary validation and remaining limits

The final combined suite passed52 tests, including3 opt-in real-engine tests; skill-creator structural validation and local Markdown-link checks passed. A final-source regional build reproduced the earlier regional counts and statuses after zero-protein status handling was repaired.

Tests cover isoforms, identical proteins at distinct loci, namespace collisions, missing proteins, unknown intervening loci, same-start ranks, overlap/clipping, zero results, incomplete native classification, strict vs proximal pairs, nested family labels, cache corruption/parameter invalidation, coordinate conversions and actual figure endpoint evidence.

Opt-in tests on the real macOS MCScanX prove that options after the prefix can be ignored, options before it apply, and `-a` suppresses the HTML path. A score100 two-anchor fixture is accepted at cutoff100 and rejected at150; block E4 is rejected at threshold4 and accepted immediately above4. Fractional gap penalties and coordinates outside signed32-bit range are rejected by the helpers.

This package does not verify gene-model biology, specific WGD events, duplication age, orthology exclusivity or pathway enzyme activity. Short/overlapping models and low-coverage pairs need family-specific review. The final intended journal size must be visually inspected. Software environments are external dependencies; this skill does not bundle MCScanX, DIAMOND, BLAST+, JCVI, Circos or a genome database.
