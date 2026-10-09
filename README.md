# Genome synteny and duplication

A Codex skill for genome synteny and duplication analysis. Build reusable synteny evidence from complete representative proteomes and gene coordinates, then overlay any focal gene family. Preserve original gene IDs, parameters, input/output hashes and figure source tables.

## Features and scope

- Strict tandem duplication candidates, proximal duplication candidates and WGD/segmental collinearity candidates within genomes.
- Syntenic blocks, direct anchor pairs and all genes within each block across subgenomes, species or assemblies.
- Reuse existing genome-wide analyses when adding gene families, without reselecting representative proteins based on family membership.
- Generate circular plots with Circos, two-genome chromosome plots with JCVI, and dotplots and local gene-order plots with Matplotlib. Export the genes and edges actually plotted, captions and omission counts.

The skill is in [`skills/genome-synteny-duplication`](skills/genome-synteny-duplication/SKILL.md). Keep synteny, adjacency and family membership as separate evidence types; these alone do not establish enzyme function, a specific WGD event or gene loss. Inspect figures at their final publication size.

## Install in Codex

Clone the repository and enter its root directory:

```bash
git clone https://github.com/yz3394/Genome-synteny-duplication.git
cd Genome-synteny-duplication
```

Copy the entire `skills/genome-synteny-duplication` directory into the Codex `skills` directory. **If a local skill with the same name already exists, preserve and compare local changes before replacing it.** The following commands skip copying when the destination already exists:

```bash
skill_target="${CODEX_HOME:-$HOME/.codex}/skills/genome-synteny-duplication"
if [ -e "$skill_target" ]; then
  echo "Destination already exists. Preserve and compare local changes first. No files were copied."
else
  mkdir -p "${CODEX_HOME:-$HOME/.codex}/skills"
  cp -R skills/genome-synteny-duplication "$skill_target"
fi
```

Invoke `$genome-synteny-duplication` in a Codex session that can read the skill. Provide matching GFF3, protein FASTA, sequence lengths or assembly FASTA, and validated family member lists.

## Requirements

| Purpose | Dependency |
|---|---|
| Input preparation, evidence processing, JSON configuration | Python ≥3.10, standard library |
| Whole-proteome searches | DIAMOND, or BLAST+ `blastp` and `makeblastdb` |
| Synteny inference | MCScanX; native gene classification also requires `duplicate_gene_classifier` |
| Local plots and dotplots | Matplotlib |
| Two-genome chromosome plots | JCVI and its runtime dependencies |
| Circular plots | Circos, Matplotlib, Pillow |
| YAML configuration (optional) | PyYAML |

This repository does not include these external tools, databases or genomes. Development validation was performed on macOS; other platforms and software versions require checks of effective output parameters and relevant tests. Development used MCScanX commit `0956cb8f900c152e2be8ba3196829e50ab454b94`; record the actual versions used in each new environment.

## Basic usage

Run from the repository root. First copy and edit [`analysis_config.example.json`](skills/genome-synteny-duplication/assets/analysis_config.example.json), replacing genome metadata, input locations and executable paths. Input paths are resolved relative to the configuration file.

```bash
python3 skills/genome-synteny-duplication/scripts/synteny_workflow.py prepare --config /path/to/analysis.json
python3 skills/genome-synteny-duplication/scripts/synteny_workflow.py build --config /path/to/analysis.json
python3 skills/genome-synteny-duplication/scripts/synteny_workflow.py overlay --run /path/to/run_manifest.json --families /path/to/families.tsv
```

The default starting profile uses DIAMOND sensitive, protein E≤1e−5 and a reporting cap of 50. MCScanX receives the top 5 distinct nonself targets per query. Tandem candidates retain a wider homology pool. These parameters are a starting point; perform sensitivity checks as needed for hit saturation, sequence divergence and focal loci.

- [Input formats and evidence tables](skills/genome-synteny-duplication/references/inputs-and-evidence.md)
- [Methods, parameters and primary literature](skills/genome-synteny-duplication/references/methods-and-parameters.md)
- [Plotting commands, configuration and checks](skills/genome-synteny-duplication/references/figures-and-qa.md)

## Validation scope

Development validation on 2026-09-28 included five whole-proteome comparisons involving Coffee CC/CE and Catharanthus roseus v3, regression against historical results, parameter and representative-transcript sensitivity checks, six regional search profiles, and actual rendering of Circos, JCVI and local plots. See [validated scope](skills/genome-synteny-duplication/references/coffee-regression.md) for details. Project raw data and complete run outputs are not distributed with this repository.

The fully enabled environment passed **52 tests** at that time. In a new environment, run:

```bash
python3 -m unittest discover -s skills/genome-synteny-duplication/tests
SYNTENY_MCSCANX=/absolute/path/to/MCScanX python3 -m unittest discover -s skills/genome-synteny-duplication/tests
```

Without `SYNTENY_MCSCANX`, 3 real-engine tests are skipped; without Matplotlib, 2 additional rendering tests are skipped. Historical validation does not validate new data, platforms or dependency versions, and the regional BLASTP checks do not establish whole-proteome search equivalence.

## Local updates and GitHub synchronization

Future updates of the installed local skill are validated and synchronized to this repository only when the repository owner explicitly requests a GitHub update. Hourly checks have been cancelled. See [MAINTENANCE.md](MAINTENANCE.md) for synchronization direction, validation requirements and conflict handling. File copying, Git commits and GitHub pushes are separate steps; the synchronization script does not push.
