# LTRquest

[![CI](https://github.com/cwb14/LTRquest/actions/workflows/ci.yml/badge.svg)](https://github.com/cwb14/LTRquest/actions/workflows/ci.yml)
[![container](https://github.com/cwb14/LTRquest/actions/workflows/docker.yml/badge.svg)](https://github.com/cwb14/LTRquest/actions/workflows/docker.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-brightgreen)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.9-blue)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue)](LICENSE)

**Find full-length LTR retrotransposons, including the ones hiding inside other ones.**

Most LTR-RT callers scan a genome once. An element that has another element
inserted into it no longer looks like an LTR-RT — the insertion breaks the
LTR–internal–LTR structure the scanner is matching — so nested elements are
missed, and so is the host they landed in.

LTRquest scans, masks what it found, and scans again. Each round exposes the
next layer out. The rounds are then reconciled into **depth buckets**:

```
depth0:  [===============]                     nothing nested inside
depth1:  [=====[depth0]=====]                  one element nested inside
depth2:  [===[depth1 [depth0] ]===]            two layers nested inside
```

On a simulated genome with 4,468 true intact elements (20 threads, scored at
≥90 % reciprocal overlap):

| Tool | Runtime | TP | FP | FN | F1 |
|---|---|---|---|---|---|
| **LTRquest** | **5:18** | 4,328 | 170 | 140 | **0.966** |
| EDTA raw LTR module (LTR_retriever) | 20:28 | 2,085 | 0 | 2,383 | 0.640 |

---

## Install

LTRquest orchestrates a stack of external tools (GenomeTools, LTR_finder,
MMseqs2, HMMER, miniprot, …), so the container is the path of least resistance:

```bash
docker pull ghcr.io/cwb14/ltrquest:1.0.0
```

<details>
<summary>conda / mamba</summary>

```bash
mamba env create -f environment.yml
mamba activate ltrquest
pip install .
```

Or, once the Bioconda recipe lands:

```bash
mamba create -n ltrquest -c conda-forge -c bioconda ltrquest
```
</details>

<details>
<summary>pip (Python parts only)</summary>

```bash
pip install ltrquest
```

This gives you the `ltrquest` driver and every stage, but you must supply the
external binaries yourself — see [`environment.yml`](environment.yml) for the
list. Kmer2LTR, TEsorter2, TRF-mod and sdust are fetched and built into
`--tools-dir` on first use; the container pre-builds them and never reaches the
network.
</details>

## Quickstart

The repo ships a real chromosome (*Arabidopsis thaliana* chr2) to check the
install. One round, a few minutes:

```bash
mkdir athal && cd athal

ltrquest \
  --genome   ../tests/data/Athal_tair10_chr2.fa.gz \
  --proteins ../tests/data/Athal.pep.gz \
  --threads  20 \
  --max-rounds 1
```

Drop `--max-rounds 1` for the full nested search. `--proteins` accepts a protein
FASTA from *any* related species — it does not have to match the genome.

## A worked example

Run the quickstart without `--max-rounds` and LTRquest iterates, masking between
rounds, until a round finds fewer than `--terminate_count` (default 100) new
elements:

```
Round 1 / 10 (r1) ......... detected <n1> LTR-RTs
  Masking original genome for next round: feature-character=N
Round 2 / 10 (r2) ......... detected <n2> LTR-RTs      # elements inside round-1 hits
  Masking original genome for next round: feature-character=R
Round 3 / 10 (r3) ......... detected <n3> LTR-RTs
  Terminate: <n3> < 100. No further rounds will be run.
Reconciling 3 round(s) into depth-bucketed outputs...
```

Each round gets its own IUPAC mask character, which is how the reconciler later
works out who is nested in whom. The elements land in one file per nesting
depth:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv     un-nested elements
Athal_tair10_chr2_LTRs_depth0_clean_ltr.fa
Athal_tair10_chr2_LTRs_depth1_clean_ltr.tsv     one element nested inside
Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa
Athal_tair10_chr2_LTRs_all_depth_LTR_cleaned.gff3
Athal_tair10_chr2_LTRs_plots/                   structure PDFs + TEGV browser
```

Each row of a `.tsv` is one element, keyed by coordinate and classification,
with its LTR–LTR divergence, an insertion age, the TSD, the protein domains and
its nesting relations. Two rows from a real run, abridged to the columns that
matter here:

```
#name                              LTR_len  …  K2P_d     K2P_T    …  domains                                          nest_status
chr1:907045-918753#LTR/Copia/…     340      …  0.104197  1736609  …  INT|Reina@907840-908727;RH|Reina@909001-909432   nest-outer:chr1:909422-915750
chr1:909422-915750#LTR/Copia/…     334      …  0.068221  1137014  …  RH|Bianca@910002-910364;RT|Bianca@910695-911483  nest-inner:chr1:907045-918753
```

Read: an element inserted ~1.74 Mya, with a younger one at ~1.14 Mya sitting
inside it — `nest-outer` on the host points at its tenant, `nest-inner` on the
tenant points back. In the `depth1` FASTA the host record carries that insertion
hard-masked as `N`, so the host's own sequence stays contiguous.

Full column and GFF3 reference: **[docs/outputs.md](docs/outputs.md)**.
Every flag: `ltrquest --help`.

## Nextflow

The round loop is also a Nextflow DSL2 pipeline, laid out to nf-core
conventions, for running many genomes in parallel on HPC or in the cloud:

```bash
nextflow run cwb14/LTRquest -profile test,docker --outdir results
nextflow run cwb14/LTRquest -profile docker --input samplesheet.csv --outdir results
nextflow run cwb14/LTRquest -profile awsbatch --input s3://bucket/samplesheet.csv --outdir s3://bucket/results
```

See **[docs/nextflow.md](docs/nextflow.md)** for the samplesheet format, the
profiles, and what the pipeline does and does not port from the CLI driver.

## How it fits together

```
ltrquest (driver)
  │
  ├── round 1..N ──┬─ ltrquest.detect     LTRharvest + LTR_finder → TEsorter2 → Kmer2LTR
  │                └─ ltrquest.mask       mask this round's hits → next round's genome
  │
  ├── ltrquest.reconcile                  pool rounds, resolve containment → depth buckets
  ├── Kmer2LTR + flag_fp_families         cluster into families, purge false-positive families
  ├── ltrquest.annotate                   add strand + family columns
  ├── ltrquest.gff3                       pooled GFF3
  └── ltrquest-plots                      structure PDFs, summary PDF, TEGV browser
```

Every stage is a module with its own `--help` and can be re-run on its own
(`python -m ltrquest.reconcile --help`).

## Citing

LTRquest was extracted from [synLTR](https://github.com/cwb14/synLTR) and shares
its false-positive family caller with
[Kmer2LTR](https://github.com/cwb14/Kmer2LTR). If you use it, please cite this
repository.

## License

[GPL-3.0-or-later](LICENSE).
