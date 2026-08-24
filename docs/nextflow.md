# Running LTRquest under Nextflow

The CLI (`ltrquest`) handles one machine and a handful of genomes. The Nextflow
pipeline exists for the case the CLI is bad at: many genomes, a scheduler, and a
run long enough that you want `-resume` when a node dies.

Nextflow fetches the container itself, so the only thing you install is
Nextflow:

```bash
# 30 seconds, no bioinformatics tools needed -- checks the wiring
nextflow run cwb14/LTRquest -profile test,singularity --outdir results

# the real thing
nextflow run cwb14/LTRquest -profile singularity --input sheet.csv --outdir results
```

On a cluster, `module load apptainer` (or `singularity`) first. `-profile
singularity` drives Apptainer too -- Nextflow treats them as the same engine --
and there is a separate `-profile apptainer` if you would rather be explicit.

## Samplesheet

```csv
sample,genome,proteins
athaliana,genomes/Athal.fa.gz,proteins/Athal.pep.gz
alyrata,genomes/Alyrata.fa.gz,proteins/Athal.pep.gz
```

| Column | |
|---|---|
| `sample` | Names this genome's outputs. Becomes `<sample>_LTRs`. |
| `genome` | FASTA, plain or gzipped. |
| `proteins` | Optional but strongly recommended. May come from *any* related species — one file can serve every row. |

Relative paths resolve against the samplesheet's own directory, not the launch
directory, so a sheet and its data stay portable as a unit. Absolute paths and
`s3://` URIs are used as given.

> **Sequence IDs must be unique across genomes.** Family clustering is pooled
> and keys elements on `chrom:start-end`, so two genomes that both call a
> chromosome `Chr2` will cross-assign families and cross-purge real elements.
> The CLI refuses to start in that case; the Nextflow pipeline does not check,
> so rename first:
>
> ```bash
> awk '/^>/{sub(/^>/,">Aly_")}1' Alyrata.fa > Alyrata.renamed.fa
> ```

## Profiles

| Profile | |
|---|---|
| `singularity` / `apptainer` | The published image, no root needed. The usual choice on a shared cluster. Everything is baked in; runs offline. |
| `docker` | The same image, on a machine where you have root. |
| `conda` | Builds the per-module environments. Kmer2LTR and TEsorter2 are not on Bioconda, so this path still fetches them on first use. |
| `test` | Minimal wiring check. Pair with `-stub-run` and it needs no tools at all. |
| `test_full` | The real *Arabidopsis* chr2 that ships with the repo, all rounds. |
| `awsbatch` | AWS Batch against an S3 work directory — see below. |

Resource ceilings come from `process.resourceLimits`. Override them for your
site rather than editing the repo:

```groovy
// site.config
process {
    resourceLimits = [ cpus: 64, memory: '512.GB', time: '72.h' ]
    executor = 'slurm'
    queue    = 'genomics'
}
```

```bash
nextflow run cwb14/LTRquest -c site.config -profile singularity --input sheet.csv --outdir results
```

Nextflow caches the image once and reuses it for every task and every later run;
set `NXF_SINGULARITY_CACHEDIR` to somewhere with room if your home directory is
small.

## AWS Batch

```bash
nextflow run cwb14/LTRquest -profile awsbatch \
    --input      s3://my-bucket/samplesheet.csv \
    --outdir     s3://my-bucket/results \
    --aws_queue  my-batch-queue \
    --aws_region us-east-1 \
    -work-dir    s3://my-bucket/work
```

Credentials come from the usual chain (environment, `~/.aws/credentials`, or an
instance role). The Batch compute environment needs an AMI with Docker and
enough instance storage for a genome plus its suffix array — detection is the
memory- and disk-hungry step, and the `process_high` label asks for 16 CPUs and
64 GB. Spot instances are worth it: `conf/awsbatch.config` retries on Batch's
host-terminated exit code.

## What the pipeline does

```
per genome ─┬─ ROUND_01: detect ──► mask ─┐
            │                             ▼
            ├─ ROUND_02: detect ──► mask ─┐   ... up to --max_rounds
            │                             ▼
            └─ ROUND_0N: detect            
                   │
                   ▼
              reconcile            pool rounds, resolve containment → depth buckets
                   │
  ═════════════════▼══════════════ pooled across ALL genomes ═════════════════
              cluster              one Kmer2LTR pass → the shared family basis
                   │
              flag-fp              purge false-positive families
  ═════════════════▼══════════════ back to per genome ═══════════════════════
              annotate             strand + family columns
                   │
              gff3 ──► plots
```

The pooled middle is the point of running several genomes together: families
computed per species are not comparable between species, and a repeat that looks
convincing in one genome but wrong across several is only caught when they are
judged together. It is also a barrier — every genome's detection must finish
before clustering starts.

### The round chain

Round *N* scans the genome round *N−1* masked, so the rounds are a strict chain.
Nextflow gives no compact way to express that: a process may be invoked only
once per workflow, calling a workflow from inside a closure is rejected, and the
`recurse` operator is still behind `nextflow.preview.recursion` and accepts only
value channels — which would forfeit the per-genome parallelism that is the
reason to be here at all. So `subworkflows/local/ltrquest_detect_rounds` includes
the round subworkflow ten times and writes the chain down.

Nothing is wasted. Rounds past `--max_rounds` are starved at the gate, and a
round that finds fewer than `--terminate_count` elements emits no successor
state, so every later round sees an empty channel and is never scheduled.
Nextflow's laziness is the loop's `break`.

## Differences from the CLI

One difference, and it is deliberate.

When false positives are **pervasive** — above `--fp_mask_threshold` — the CLI
hard-masks those repeats in the genome and re-runs the entire pipeline on the
masked copy, up to ten times. The Nextflow pipeline does a single pass and
reports the fraction in `results/families/*_fpcheck.log` instead.

That outer loop is a chain of *whole-pipeline* iterations, each of whose
existence depends on the previous one's output. It is the same constraint as the
round chain but an order of magnitude more expensive to unroll, and it fires
rarely. If your `_fpcheck.log` reports a fraction above the threshold, either
raise it, or run that genome through the CLI, which handles the re-run properly:

```bash
grep 'FP fraction' results/families/*_fpcheck.log
```

Everything else — the rounds, the reconciliation, the pooled family basis, the
annotation, the GFF3, the plots — is the same code the CLI calls, invoked the
same way.

## Development

Every process carries a `stub:` block, so the whole DAG can be exercised with no
external tools at all:

```bash
nextflow run . -profile test -stub-run --outdir results
```

That is what CI runs on every push, against both the oldest supported Nextflow
(24.04) and the current release. The stubs are written to be faithful about the
things the DAG depends on — the detector's stub returns two elements in round 1,
one in round 2 and none after, which is enough to drive the termination gate;
the false-positive stub derives its output names from its inputs, because the
per-genome split downstream keys on exactly those names.
