# Inputs, identity and evidence

## Configuration

`output_dir` is the project artifact root. Paths are resolved relative to the config file. `genomes` contains unique `id`, `species`, `assembly_id`, `annotation_id`, `proteins`, and exactly one of `gff` or `genes_tsv`. Optional fields: `subgenome`, `rank_policy`, `assembly_fasta`, `lengths_tsv`, `seq_lengths`, `id_map_tsv`, `representative_map_tsv`.

- `gff`: GFF3, plain or gzip. gene→mRNA/transcript→CDS Parent relations and exact protein identifiers are used. Repeated CDS segments are permitted. Unsupported naming requires an explicit map, not suffix guessing.
- `genes_tsv`: normalized loci with `gene_id,seqid,start,end,strand,protein_id`; multiple isoform rows may describe the same identical locus. It is an adapter input, not a second GFF interpretation.
- `id_map_tsv`: explicit `gene_id,transcript_id,protein_id`. Contradictory/ambiguous mappings fail.
- `representative_map_tsv`: frozen `gene_id,protein_id` overrides, independent of the family table. Other loci use longest protein, lexical smallest protein ID on equal lengths. Terminal stops are recorded/removed; internal stops require model review and explicit repaired input.
- Sequence lengths: use assembly FASTA, `lengths_tsv` with `seqid,length`, or explicit `seq_lengths`. With none, the maximum gene end is labeled inferred, and full chromosome physical lengths remain unverified.
- `rank_policy`: `protein_coding` is the helper default. It retains loci with protein/mRNA/CDS or explicit coding evidence, excludes explicitly noncoding/pseudogene loci without coding evidence, and conservatively retains unresolved generic loci with an `unknown_biotype_retained` flag, count and warning. Missing protein does not by itself remove a locus. `all_genes` retains all annotated loci and supports legacy comparisons. Record the actual universe, not just its label; unresolved intervening loci must not silently disappear and create adjacency.

`comparisons`: list of `id,a,b,kind`. Self comparisons use `intragenome`; cross comparisons use `inter_subgenome`, `interspecies` or `inter_assembly`. The latter describes sequence sets, not extra duplication events. Structural ancestry/homeology still needs biological context.

`tools`: executable names or paths for `diamond`/`blastp`, `mcscanx` and optional `classifier`. Native gene classes are explicitly not run when no classifier was configured. `search`, `collinearity` and `duplication` keys are in the example and methods reference. `families` and plot settings do not enter any core fingerprint.

## Fixed prepared artifacts

`genes.tsv` preserves original IDs, selected protein, coordinates, strand, annotation/engine ranks and quality flags. Engine IDs/seqids are generated within a genome namespace and never replace user-facing IDs. `all_isoforms.tsv` supports later mapping of family evidence proteins that were not selected for genome-wide search.

Annotation rank sorts by `(start,end,gene_id)`; engine rank sorts by `(start,engine_id)` per chromosome. Tandem adjacency uses the declared annotation rank. Block membership uses engine rank. Same-start and overlap flags must be inspected when these ranks differ. Positions include eligible missing-protein loci, while the search FASTA only includes available representatives.

`chromosomes.tsv` records lengths, namespace and length provenance. The four-column `positions.gff` is an MCScanX adapter output, not a replacement for the original GFF3. Original sources and hashes remain in `input_manifest.json`.

The tested MCScanX stores gene positions in signed 32-bit integers. Eligible coordinates above 2,147,483,647 are rejected explicitly; do not truncate or silently replace physical coordinates with ranks. Such genomes need a separately validated engine/coordinate strategy.

## Evidence tables

Coordinates in normalized tables are 1-based inclusive. Native BED/Circos inputs use an explicit conversion.

| Table | Grain and meaning |
|---|---|
| `homology_pairs.tsv` | One canonical best-HSP homology pair with directional support; raw m8 retains all reported HSP evidence. Coverage is coordinate span / each protein length, not aligned-length shortcut |
| `blocks.tsv` | One engine block with A/B extents, direction, anchor count, score and block E |
| `anchor_pairs.tsv` | One block×anchor relation; cross sides normalized to configured A/B |
| `unique_anchor_pairs.tsv` | A unique pair with all supporting block IDs |
| `all_block_genes.tsv` | block×side×gene inside the first/last anchor engine-rank interval; explicit anchor flag and real partners |
| `block_physical_overlaps.tsv` | Extra physical overlaps outside that rank interval; not silently added to rank membership |
| `tandem_pairs.tsv` | Same seqid, adjacent annotated loci and a passing homology hit, with model/coverage flags |
| `tandem_arrays.tsv` | Connected components of strict adjacent edges; gene_count and edge_count separate |
| `proximal_pairs.tsv` | Supported nonadjacent pairs within the explicitly configured intervening-gene limit |
| `gene_duplication_classification.tsv` | Native mutually exclusive category plus independent tandem/anchor evidence and assessment status |

Family TSV requires `genome_id`, at least one of `gene_id/protein_id`, and `family_label` or semicolon-separated `families`. If both IDs are supplied they must map to the same locus. A locus can have several labels. The overlay exports source protein IDs separately from the fixed analysis representative.

Family pair tables retain pairs with at least one focal endpoint; use `both_in_focal_families` and `same_family_pair/shared_families` to select bilateral or same-family relationships. Do not count all those rows as same-family duplications. Summaries must distinguish unique genes, pairs, arrays, block occurrences and family-label assignments.
