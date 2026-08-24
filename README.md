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

LTRquest drives a stack of external tools (GenomeTools, LTR_finder, MMseqs2,
HMMER, miniprot, …). You do not have to install any of them — the container has
them all, and needs neither `sudo` nor a package manager.

### Apptainer / Singularity — recommended

The right choice on a shared cluster, where you cannot install Docker and do not
have root. The whole tool is one file you own:

```bash
apptainer pull ltrquest.sif docker://ghcr.io/cwb14/ltrquest:1.0.1
```

That is the install. Put `apptainer exec ltrquest.sif` in front of any command:

```bash
apptainer exec ltrquest.sif ltrquest --help
```

Apptainer already mounts your working directory and already runs as you, so
there are no paths to bind or permissions to set — commands look exactly as they
would if the tool were installed natively.

> **On a cluster** you may need `module load apptainer` (or `module load
> singularity`) first. `singularity` works identically — just swap the command
> name. If your home directory is small, send the download cache to scratch
> before pulling:
> ```bash
> export APPTAINER_CACHEDIR=/your/scratch/.apptainer
> ```

### conda / mamba

If you would rather have the tools on your `PATH` than in a container:

```bash
mamba env create -f environment.yml
mamba activate ltrquest
pip install .
```

<details>
<summary><b>Docker</b> — for your own machine, where you have root</summary>

```bash
docker pull ghcr.io/cwb14/ltrquest:1.0.1
```

Docker has to be told which directory to work in and which user to be, so its
commands are longer than the Apptainer equivalents:

```bash
docker run --rm -v "$PWD:/data" -w /data -u "$(id -u):$(id -g)" \
  ghcr.io/cwb14/ltrquest:1.0.1 ltrquest --help
```
</details>

<details>
<summary><b>pip</b> — Python parts only, external tools not included</summary>

```bash
pip install ltrquest
```

You then have to supply the external binaries yourself; see
[`environment.yml`](environment.yml) for the list. Kmer2LTR, TEsorter2, TRF-mod
and sdust are fetched and built into `--tools-dir` on first use. The container
pre-builds all four and never reaches the network.
</details>

### One command, whichever route you took

Every example below is written as plain `ltrquest …`. If you installed with
conda or pip, that already works. If you are using a container, either put
`apptainer exec ltrquest.sif` in front, or set this once and forget it:

```bash
alias ltrquest='apptainer exec /full/path/to/ltrquest.sif ltrquest'
```

## Try it on real data

LTRquest ships a real chromosome — *Arabidopsis thaliana* chr2, 19.7 Mb — so you
can see it work on something real before pointing it at your own genome.

Copy-paste the whole block. It takes a few minutes on 20 cores:

```bash
mkdir ltrquest-demo && cd ltrquest-demo

# the tool
apptainer pull ltrquest.sif docker://ghcr.io/cwb14/ltrquest:1.0.1

# the data: one chromosome, and a protein set to guide the search
curl -sLO https://github.com/cwb14/LTRquest/raw/main/tests/data/Athal_tair10_chr2.fa.gz
curl -sLO https://github.com/cwb14/LTRquest/raw/main/tests/data/Athal.pep.gz

# go
apptainer exec ltrquest.sif ltrquest \
  --genome Athal_tair10_chr2.fa.gz --proteins Athal.pep.gz --threads 20 --max-rounds 1
```

<details>
<summary>The same run, if you installed another way</summary>

**conda / mamba, or pip:**

```bash
ltrquest --genome Athal_tair10_chr2.fa.gz --proteins Athal.pep.gz --threads 20 --max-rounds 1
```

**Docker:**

```bash
docker run --rm -v "$PWD:/data" -w /data -u "$(id -u):$(id -g)" \
  ghcr.io/cwb14/ltrquest:1.0.1 \
  ltrquest --genome Athal_tair10_chr2.fa.gz --proteins Athal.pep.gz --threads 20 --max-rounds 1
```

Already cloned the repo? The same files are in `tests/data/`.
</details>

When it finishes, open **`*_plots/*_TEGV.html`** in a web browser: a genome
browser of everything it found, in one self-contained file you can email to a
collaborator.

### Your own genome

Three flags cover almost every run:

```bash
ltrquest --genome your_genome.fa --proteins any_related_species.pep.fa --threads 20
```

| | |
|---|---|
| `--genome` | Your genome, `.fa` or `.fa.gz`. **The only required flag.** |
| `--proteins` | Optional but recommended. Can come from **any related species** — it does not have to match your genome. |
| `--threads` | How many cores to use. Default 20 — set it to what your machine actually has. |

Dropping `--max-rounds 1` is what turns on the nested search: LTRquest keeps
going until a round stops finding anything new. Everything else has a sensible
default — `ltrquest --help` lists them.

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

Have more than a couple of genomes? The round loop is also a Nextflow DSL2
pipeline, laid out to nf-core conventions, which runs them in parallel and
picks up where it left off if a node dies. Nextflow fetches the container for
you — there is nothing to install but Nextflow itself:

```bash
nextflow run cwb14/LTRquest -profile singularity --input samplesheet.csv --outdir results
```

The samplesheet is three columns:

```csv
sample,genome,proteins
athaliana,genomes/Athal.fa.gz,proteins/Athal.pep.gz
alyrata,genomes/Alyrata.fa.gz,proteins/Athal.pep.gz
```

Swap `-profile singularity` for `docker`, `conda`, or `awsbatch` to run the same
pipeline on your laptop or on AWS. See **[docs/nextflow.md](docs/nextflow.md)**
for the profiles, the samplesheet rules, and the one place the pipeline
deliberately differs from the CLI.

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

Every stage is also its own command, so any of them can be re-run alone without
redoing the rest — `ltrquest-reconcile --help`, `ltrquest-gff3 --help`, and so
on.

## Citing

LTRquest was extracted from [synLTR](https://github.com/cwb14/synLTR) and shares
its false-positive family caller with
[Kmer2LTR](https://github.com/cwb14/Kmer2LTR). If you use it, please cite this
repository.

## License

[GPL-3.0-or-later](LICENSE).
