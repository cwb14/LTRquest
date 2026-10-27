# Changelog

All notable changes to LTRquest are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Fixes the three ways a fresh install could fail on a machine that is not a
developer workstation. All three were found by installing from scratch — pulling
the published image and creating the conda environment — rather than by testing
the source tree.

Upgrades the bundled Kmer2LTR from its original k-mer/MAFFT implementation to
the rewritten package, and reorders the pipeline around what that rewrite
makes possible. Kmer2LTR now runs *ahead of* classification instead of after
it, so its boundary calls, K2P divergence figures and TSDs are the round's
single source of truth rather than figures LTRquest re-derived and could
disagree with. The new tool's CLI and output schema are not compatible with
the old one's, which forces the schema, flag and stage changes below.

### Added

- **FastLTR, a third LTR caller** (`--fastltr rescue|all|off`, default `rescue`,
  or `--no-fastltr`; `--fastltr-args`; `ltrquest-detect` also takes
  `--fastltr-min-ltr`; Nextflow `--fastltr`, `--ltr_timeout`).
  `--ltr-timeout` stops an LTRharvest or LTR_FINDER call that runs too long on a
  chunk, and neither tool writes anything before it finishes, so a capped 500 kb
  chunk lost every candidate of that caller, not only its satellite's: two CHM13
  chunks capped at 30 s gave 0 candidates from either tool, and a 30 s cap fired
  on 60 of 302 *A. thaliana* T2T chunks (its five centromeres and both NORs).
  LTRquest now ships its own adaptation of Gene Myers' FastLTR (FASTAN, the
  detector behind HiRepEat; `src/ltrquest/fastltr/FastLTR.c`, compiled against
  FASTAN 843d0ed in `--tools-dir` on first use, and in the image at build time;
  FASTAN's license ships beside it as `LICENSE.FASTAN`).
  It searches LTR spacings of 100-20,000 bp (FastLTR's own floor of 3 kb, 2 kb in
  HiRepEat's build, leaves out 63%, or 47%, of an earlier run's CHM13 calls; raw
  recall of those calls rose from 20% to 63%, 82% of the high-confidence ones); a
  tandem alignment marks the bases it covers so other spacings do not re-align
  through the same array (a CHM13 p-arm chunk went from 212 s to 0.97 s); seeds
  are counted over +-8 diagonals so indels between old LTRs do not hide them;
  masked bases are deleted before the call and hits lifted back, so FastLTR sees
  an outer's LTR pair once its inner is masked; and a pair whose copy recurs one
  spacing further out is dropped as a tandem-array unit (FastLTR reports every
  direct-repeat pair, where LTRharvest keeps the best of overlapping ones).
  `rescue` runs it only on chunks a cap fired on (with no cap it never runs and is
  not built); `all` runs it on every chunk. Its candidates join the pool before
  Kmer2LTR, and every later stage treats them like the others'. Re-running a
  directory whose genomes were detected with `ltr-timeout` before this change now
  stops, as for any changed detection setting, since their records lack FastLTR:
  pass `--no-fastltr` to reuse them as they were, or `--redetect`.

  Against uncapped runs, a 30 s cap with rescue recovers 68% (*A. thaliana*) and
  86% (rice T2T-NIP) of the curated-library-supported elements a 30 s cap alone
  loses entirely in fully capped chunks (0% and 1%), keeps 98.4% and 99.7% of
  those where no cap fired, and costs 0.37 h against 0.33 h for the cap alone and
  1.48 h uncapped (*A. thaliana*). Inside rice's CENH3 centromeres it restores 150
  elements (67 CRM; 120 curated-supported) where the cap alone leaves 97 and the
  uncapped run finds 155 (116 supported). Of 401 library LTR-RTs planted into real
  CEN180 and 45S rDNA arrays, cap plus rescue recovers 352 (88%) against 103 (26%).
  `all` is not the default: in chunks no cap fired on, the elements it adds to an
  uncapped rice run carry curated support 36% against that run's own 63% (rice:
  +392 supported, +770 unsupported elements). On a simulated genome with no arrays
  it is harmless (F1 0.951 against 0.946), so it stays available for maximum
  recall. Known gap: arrays whose unit is longer than 20 kb (human 45S rDNA,
  43 kb) escape both the array test and the tandem filter, and on CHM13 rescue's
  first attempt calls 375 `LTR/unknown` unit copies there (LTRquest 1.0.1 at a
  300 s cap: 86); after the FP correction 146 remain, 101 of them a provirus
  (MER39B LTRs, HERVIP10FH) that the rDNA unit carries.
  Benchmarks: `benchmarks/fastltr/`.

- A tandem-array filter, on by default (`--no-tandem-filter` turns it off).
  Complete assemblies include the rDNA arrays and centromeric satellites, and
  LTRharvest and LTR_FINDER pair any two copies of a repeat unit as "LTRs": the
  telomere-to-telomere *A. thaliana* Col-CC assembly gave 876 such calls from
  its two 45S rDNA arrays (~10 kb units), 43% of all its elements. The `_clean_`
  tables now drop a call, alongside the false-positive families, when its
  LTR-pair repeat carries on past both termini: at least half of the
  non-overlapping 15-mers in the 1 kb outside each LTR occur in the 2 kb of
  internal sequence next to the other LTR (`ltrquest-flagfp --tandem-genome`,
  `--tandem-min-frac`). A genuine element's flanks are unrelated target-site
  sequence and score ~0; the containment test, rather than a fixed-offset
  comparison, still matches rDNA units that differ in spacer-subrepeat copy
  number. Two elements sharing a middle LTR each continue on one side only and
  are kept; three or more head-to-tail copies count as an array. On Col-CC
  (same detection, filter off vs on) the clean set goes from 2,012 to 1,074
  calls: 859 of the 871 rDNA calls, 41 CEN180-satellite calls and 39 calls
  in other tandem duplications leave; the pooled stages run in 2:13 instead of 2:33.
  In three more plant genomes (*B. napus* A, *C. canephora*, *O. punctata*) it
  removes 4.8%, 2.0% and 0.2% of calls. The raw
  `depth<N>_ltr.tsv` tables keep these calls, `<prefix>_fpcheck.tandem.tsv`
  scores every call, and they never count toward the FP fraction that
  triggers masking. Detection is unchanged, so a genome reused from an earlier
  run gains the filter without `--redetect`. The Nextflow pipeline runs it too:
  `LTRQUEST_FLAGFP` now stages every sample's genome, and `--skip_tandem_filter`
  turns it off.

  Tried first and dropped: the same test as a gate on raw LTRharvest/LTR_FINDER
  candidates, before Kmer2LTR. Their unrefined LTR bounds let ~550 unit-length
  rDNA candidates through, and removing the rest let short spacer-subrepeat
  candidates that dedup had suppressed surface as calls: 814 rDNA-region calls
  survived.
- `--reboundary` / `ltrquest-reboundary`: move the ends of LTR-RT calls that have
  no TSD at their called ends, typically calls a detector stopped short at an indel
  or a mutation-dense patch. The templates are calls with an exact TSD, pooled over
  genomes, families, strands and nesting depths; one BLASTN finds each call's 5
  nearest (`--templates`). Where at least 2 agree (`--min-support`), an end moves
  out, or in by up to 10 bp (`--max-trim`). Calls with a TSD are never moved, and a
  second, staggered call of the same element is merged away. Kmer2LTR re-measures
  every changed element on its called LTR pair; bases with no partner in the other
  LTR count as gaps, so K2P and age are not inflated. Rewrites the `_clean_` tables
  and FASTAs, writes `<prefix>_reboundary.tsv`, and marks changed elements
  `boundary_source=templates` and `boundary_shift=<ext5>,<ext3>` in the GFF3.
  `--posthoc` updates a finished run in place with a one-time backup (`--restore`
  undoes it). Needs BLAST+. On rice (6,709 elements) it moved 1,705 calls and
  merged 21 in 35 s, and a curated rice LTR library agrees with 99.5% of the new
  ends. Benchmarks: `benchmarks/reboundary/`.

  Needs a Kmer2LTR with the homology-paired credit (`align.CREDIT_PAIRS_BY_HOMOLOGY`
  and `classify(..., spans=)`). An older one is refused: `ltrquest-reboundary` exits
  3 and the pipeline stops. The Docker image now pins Kmer2LTR `2c7e7f8`, the
  commit that adds it.

  Development builds only: the pre-release design's flags `--method`,
  `--references`, `--subfamily-jaccard`, `--min-copies`, `--max-ratio`,
  `--qc-min-n`, `--untested`, `--credit`, `--mafft` and `--cache` are gone, and its
  sidecars (decision `extended`) still load.
- Add a genome to an earlier multi-genome run without re-detecting the others:
  re-run in the same directory with the longer `--genome` list. Each detected
  genome leaves a `<prefix>.detect.json` (input checksums, detection settings);
  on any re-run, single-genome included, a genome whose record still matches is
  reused, and families are rebuilt over all genomes. A mismatch stops the run;
  `--redetect` ignores the records. Genomes are now pooled in name order, so
  families no longer depend on the order `--genome` lists them.
- `--strand-recovery conservative|balanced|sensitive`, off by default, which
  recovers strand for the LTR-RTs the strand cascade leaves at `.`. The
  unstranded element, in genome-forward orientation, is aligned against the
  elements that already carry a strand, held in coding sense, so the
  alignment's orientation is the query's own genomic strand. Four views run
  (dc-megablast and minimap2, each on the whole element and on the internal
  region); a view calls only when all of its own alignments agree, and
  disagreement between views vetoes the locus. Leave-one-out against
  LTRquest's own tesorter strand over 3616 elements: 58% recovered at 0.19%
  disagreement (conservative), 79% at 0.70% (balanced), 84% at 1.06%
  (sensitive). `--strand-recovery-ppt` adds a polypurine-tract fallback for
  loci homology cannot reach, measured at 93% and therefore off by default.

  The calls land in a new `<prefix>_strand_recovery.tsv`, which
  `ltrquest-annotate` reads as a fourth cascade tier via `--recovered-strands`
  — pinned by path, never globbed, so a sidecar left behind by an earlier run
  cannot re-strand a later one. Every primary output follows: `strand` in the
  depth tables, `strand_source=homology`/`ppt` in both GFF3s, and the depth
  FASTAs re-oriented so a recovered minus element is stored in coding sense
  exactly like a TEsorter2-called one, with `orientation` restated to match.
  Available on the Nextflow path as `--strand_recovery` and
  `--strand_recovery_ppt`, and standalone as `ltrquest-recover-strand`.

- `--mutation-rate`, forwarded to Kmer2LTR's `-u` for both the per-round
  table and the pooled family-clustering pass. Insertion age used to come
  from a fixed μ = 3×10⁻⁸ baked into the old Kmer2LTR script with no flag to
  change it; now `k2p_time` is `NA` unless you supply a rate that actually
  matches your organism.

### Changed

- **Overlapping calls follow one rule.** Detection used to compare only calls
  overlapping by >= 80% of the shorter (`--dedup-threshold`), and then kept every
  call that took part in a nesting, so 8.8-9.5% of the *Poa* clean calls crossed
  another call (98.8% of them from the same round), although real insertions are
  only ever disjoint or nested. Now any two calls that overlap without one lying
  strictly inside the other with its own LTR pair conflict, and only the stronger
  is kept: hallmarks at its own ends (TSD and TG..CA, then TSD, then TG..CA), then
  the LTR pair's `bitscore`, then its distinct protein domains (`ltrquest.overlaps`).
  The same rule replaces the cross-round shared-LTR merge in `ltrquest-reconcile`.
  Each round writes the calls it dropped to `<prefix>.conflicts.tsv` in its work
  directory. Scored against empty sites in related genomes (an independent test of
  a real insertion) on 1,087 crossing pairs with exactly one confirmed partner in
  three *Poa* genomes, the rule keeps the confirmed call 95% of the time (95.2,
  94.7, 94.8%), against 71-76% for the LTR alignment score alone. Re-detecting
  *P. supina* with and without the change (with the TSD fix below): crossing calls
  8.4% -> 0.09%, canonical calls 51.6% -> 63.2%, at a recall cost of 0.3 points on
  EDTA-confirmed insertions (95.6 -> 95.3%). The cost is partial calls that
  re-boundarying used to extend to the element's true ends, now dropped before it
  runs.
- **The element table gains a `confidence` column** (34 columns; before
  `domains`), and the GFF3 a `confidence` attribute: `high` when the reported
  ends are TG..CA with an exact TSD, `medium` with one of the two, `low` with
  neither. A young real insertion keeps both, so in a young element `low` marks a
  call that is very likely wrong; in *Poa*, 44% of calls were non-canonical and,
  among young ones, an empty site confirmed them 12-21% of the time against 73-78%
  for canonical calls of the same age.
- `ltrquest-plots`' `--min_ltr_aln` is renamed `--min_ltr_len`, and
  `flagged_false_positives.tsv`'s `aln_len_raw` column is renamed
  `ltr3_len_raw`. Both hold a 3' LTR length, not an alignment length; the
  names say what they measure.
- `ensure_trfmod` and `ensure_sdust` look in `--tools-dir`, then on `PATH`, then
  clone and build — the order `ensure_tools` has used since 1.0.1. `sdust` is
  packaged (`bioconda::sdust`) and is now in `environment.yml`, so `--run-sdust`
  no longer needs a compiler.
- Failures to find or build a helper name the conda package that provides it and
  the compiler toolchain that would build it, instead of only the tool.
- LTRquest's own `depth<N>_clean_ltr.tsv` is now Kmer2LTR's 33-column schema
  (see [docs/outputs.md](docs/outputs.md)), replacing the old 21-column one
  built around `LTR_len`, `raw_d`/`raw_T`, `JC69_d`/`JC69_T` and
  `left_trim`/`right_trim` — none of which the new boundary model or K2P
  estimator produces. `K2P_d`/`K2P_T` survive only as the GFF3 attribute
  names they always were; the table columns underneath them are now
  `k2p`/`k2p_time`.
- `flag_fp_families.py` moves out of the Kmer2LTR clone and into LTRquest as
  `ltrquest.flag_fp` / `ltrquest-flagfp`. Flagging false-positive families
  reads Kmer2LTR's clustering output; it was never a Kmer2LTR concern, only
  ever colocated with it because that happened to be where the script lived.

### Removed

- `ltrquest-detect --dedup-threshold` and the pre-classification purge that used
  the same 80% threshold; see *Overlapping calls follow one rule* above. Dropping
  the purge sends ~18% more candidates through TEsorter2 (about a minute on a
  1.2 Gb genome).
- Six flags that had nothing left to forward to once Kmer2LTR's internals
  changed underneath them: `--wfa-align` (the new tool always aligns with
  WFA, so there is no MAFFT default left to opt out of), `--kmer2ltr-max-win-
  overdisp` and `--kmer2ltr-min-retained-fraction` (tuning for the old k-mer
  window search; the new discovery stage is exact Smith-Waterman with no
  overdispersion or retained-fraction knobs), `--kmer2ltr-domains` (the old
  script's own domain calls; domains have only ever come from TEsorter2 in
  this pipeline), `--tsd-min-len` (the new Kmer2LTR accepts a TSD on an
  exact-match, two-distinct-base rule with no length floor to set), and
  `--tesorter-use-ret` (chose whether classification saw the ret span or the
  internal span of a candidate; structurally inert now that Kmer2LTR needs
  both LTRs to measure divergence, so classification must always receive the
  full-length element and the internals-only path can no longer feed it).
- `ltrquest-reconcile --scn`. It fed a shared-LTR lookup keyed on each round's
  pre-trim SCN locus, but Kmer2LTR renames every element to its post-trim
  locus before reconcile ever sees it, so the flag could no longer reach the
  elements it was meant to describe. Boundaries are read from each round's
  own `--tsv` table instead (`ltr_bounds_from_table`), which carries the
  right key and the right coordinates on the same row. The reconcile
  Nextflow module's fourth input path -- each round's SCN file, staged only
  to feed that flag -- goes with it, along with the workflow term that
  collected those files; the process no longer stages input it never reads.

### Fixed

- **A capped call that had written output would have failed its chunk.** On a
  timeout, Python's subprocess hands back what the killed call printed as bytes
  even in text mode, and the salvage wrote it with `write_text`, which raises on
  bytes. LTRharvest and LTR_FINDER print nothing before they finish, so the bug
  never fired -- and "salvaged partial SCN" always salvaged nothing. The partial
  output is now decoded; FastLTR, which prints each hit as it finds it, relies on
  that.

- **Ends 1-2 bp off a canonical boundary were reported where they were called.**
  Kmer2LTR records a TSD it finds after moving a boundary 1 bp (`tsd_offset`) but
  leaves the boundary, and detection never applied the shift, so 6.5-7.9% of *Poa*
  calls kept ends 1-2 bp off TG..CA-with-TSD in their ids, LTR coordinates, GFF3
  and FASTA (*P. supina* `Ps1B:542705-550914` for `Ps1B:542706-550914`). Its first
  hit was also not always the element's ends: a call 1 bp wide carries an exact
  4-bp TSD at the wrong ends, and a 5-mer matches by chance 1 bp off a correct call
  about one time in four. Each round now moves every element onto the nearest
  ends within 2 bp that are TG..CA with an exact 4-6 bp TSD (`ltrquest.ends`),
  re-measures it with Kmer2LTR, and removes a second call that reaches the same
  ends; elements already on such ends get that TSD in `tsd` with `tsd_offset`
  `0,0`, which also makes them re-boundarying templates. Re-boundarying and host
  rejoining read the TSD on canonical ends the same way.
- **`--reboundary` failed on a plain-gzip genome.** Its genome access used pyfaidx,
  which reads only BGZF, so `UnsupportedCompressionFormat` stopped the stage and the
  driver kept the calls as detected with a one-line warning. A gzipped genome is
  now expanded into the run directory for the stage and removed after it.
- **An element inserted into another element's LTR came out as two adjacent free
  elements.** The insertion splits the host's LTR into `[1..p] TENANT [p+1..L]`.
  LTRharvest and LTR_FINDER pair the two near-side fragments, so the host was
  called cut short at the tenant, with no TSD, and the tenant as a free element
  beside it. The full host was beyond them: its LTR alignment is broken, and with
  the tenant inside it is usually longer than their maximum LTR distance. Masking
  between rounds could not help, since host and tenant are mostly found in the
  same round. `ltrquest.reconcile` now rebuilds such hosts before it resolves
  containment (`ltrquest.rejoin`). For two calls that abut within 20 bp, it asks
  whether the bases beyond one continue the other's own LTR from the junction on:
  identity ≥ 0.70, E ≤ 1e-5, starting at both junctions. Where one side and only
  one does, the host is extended to the end of that continuation and the tenant
  becomes its `nest-inner`. Such a host's record masks the tenant inside its LTR;
  `ltr5_len`/`ltr3_len` count LTR bases only; and `motif` and the TSD columns are
  read at the new ends. Its divergence columns stay Kmer2LTR's measurement of the
  fragment pair it was called with: the same LTR pair, fewer sites. No rejoin is
  made where both sides continue, where the continuation would run through
  another call, where an element would be in two, or where a rebuilt LTR would be
  longer than the callers' 7 kb maximum. Past that length the "continuation" is
  the next unit of a tandem array, not an LTR. `<prefix>_rejoined.tsv` maps
  each old key to its new one; the strand loaders read it, so a rebuilt host
  keeps the TEsorter2 and pass-2 calls made on its fragment. Reconcile takes
  `--genome` and `--tools-dir` (Kmer2LTR's TSD rule); without `--genome` it
  rebuilds nothing and says so.

  On six *Poa* assemblies, with detection unchanged:
  - 692–1,631 hosts are rebuilt per genome, 10–12% of all nests.
  - The rebuilt ends carry a TSD in 32–40% of hosts, against 1.0–1.3% with the
    right end shifted 500 bp and 7–10% at the old, truncated ends.
  - Tenant excised, the rebuilt LTRs align end to end: identity 0.97–0.99,
    ≥ 90% coverage in 99.5% of hosts.
  - LTR length relative to the family median goes from 0.64–0.75 as called to
    1.03–1.05 rebuilt.
  - The step takes 45 s to 3 min per genome.

  On six PrinTE-simulated genomes (known truth, ≥ 90% reciprocal overlap):
  - intact-element recall goes from 0.92–0.94 to 0.94–0.97, and precision from
    0.89–0.91 to 0.91–0.93;
  - hosts with a tenant in an LTR are recovered in 61–67% of cases instead of
    15–21%; hosts with internal tenants only are unchanged (0.90–0.91);
  - 88–92% of rebuilt hosts are exactly the true host with that tenant in its
    LTR. At most 2 per genome replace a call that was already correct.

  On telomere-to-telomere *A. thaliana* Col-CC (2,172 calls, two 45S rDNA
  arrays): 6 hosts are rebuilt, none in either rDNA array. The LTR-length limit
  refuses 2 centromeric "hosts" whose continuation ran 7–8 kb into the array.

  Related changes:
  - Re-boundarying no longer targets a call whose LTR holds a nested element, as
    template selection already did.
  - `plot_summary`'s internal length subtracts only what is nested in the
    internal region.
  - Genomes detected before this change record no `rejoin` setting. A run that
    would pool them with new ones refuses them; re-detect them with
    `--redetect`.
- **A nested host could take its strand from its tenant.** Two paths copied it:
  - TEsorter2 classifies a round's elements (Step 9) before Step 9c masks the
    inners found in that same round, so a host still held its tenant's
    domains. Step 9b already rebuilt the host's clade from its own hits but
    kept TEsorter2's strand, which is the one strand all hits share. A host
    with no domain hits of its own therefore got its tenant's strand. Step 9b
    now sets that strand to `?` unless the host's own hits (outside its
    inners) all lie on it; it never invents one. Both of the round's cls
    tables (`*.cls.tsv` and `*.rexdb.cls.tsv`) are rewritten, because
    `ltrquest.annotate` reads both and a `+`/`-` in either beats a `?`. The
    library orientation, `ltrquest.annotate`, `ltrquest.gff3` and
    `ltrquest.recover_strand` all read the corrected strand.
  - `ltrquest.recover_strand` read each element from the plain genome, so a
    host's sequence still held its tenant, and a host and tenant aligned to
    each other at 100% identity. Each element, as query or reference, is now
    read with its `nest_status` inners masked to N.

  On *P. annua* (5 rounds), 1,408 hosts lose a strand. The polypurine tract
  agreed with that strand for only 20 of 45 PPT-called hosts (44%, binomial
  p = 0.55). For the 958 hosts that keep their strand it agrees 104/109
  (95%). In strand recovery on chromosome Pa1A (`--strand-recovery
  sensitive --strand-recovery-ppt`), homology-recovered tenants were 497/497
  co-oriented with their host before the fix and 180/247 after. Pairs of
  different clades come out at 54% (35/65, tenant with a TSD), and the stage
  now recovers 1,025 strands instead of 717. Unmasked tenants had been
  creating conflicting alignments.
- **FP masking crashed on a gzipped genome.** When a genome's FP fraction passed
  `--fp-mask-threshold`, `ltrquest.flag_fp` ran `makeblastdb` on
  `<prefix>.input_genome.fa`, a symlink to the user's genome; a `.fa.gz` input
  failed with "Input doesn't start with a defline", hours into the run. The
  genome is now streamed into `makeblastdb` on stdin, decompressed if need be,
  and the masked genome is written from plain or gzip input alike.
- Nextflow: a genome that completed a single detection round (`--max_rounds 1`,
  or a first round below `--terminate_count`) crashed `LTRQUEST_RECONCILE`. One
  table arrives as a Path, not a list, and iterating it walked the path's
  components, so the reconciler was called as `--tsv round_1 x_r1_ltr.tsv`. The
  three-round stub test never reached the case.
- **`--reboundary` could abort on one pair of calls and move nothing.** When
  detection leaves two calls of one element, one inside the other, one call can be
  extended and the other trimmed to the same ends. Neither move conflicted with the
  other call as it stood, so both were accepted, and the rewrite refused to give two
  elements one key: `ValueError: ... would both be written with the same key`, then
  `re-boundarying failed; the calls are kept as detected` for every genome in the
  run (6 of 9 allotetraploid runs, 4 genomes each). The first move in genome order
  now wins and the other call keeps its ends (reason `duplicates_move`); a trim onto
  the exact span of a call inside is refused (`duplicates_element`).
- **`./ltrquest.sif --genome x.fa` did not work.** The image declared
  `ENTRYPOINT []` with only a `CMD`, and Apptainer builds a SIF's runscript from
  those two: with no ENTRYPOINT it *replaces* the command with the user's
  arguments, so the first flag was run as a program —
  `FATAL: "--help": executable file not found in $PATH`. The image now has an
  entry point ([`bin/entrypoint.sh`](bin/entrypoint.sh)) that runs a first
  argument naming a command and treats anything else as `ltrquest` arguments.
  `apptainer exec IMG ltrquest …` and `docker run IMG ltrquest …` are unchanged
  — `exec` never reaches an entry point, and an explicit command still wins.
- **A working directory reached through a symlink left the run in `$HOME`.**
  Apptainer mounts your home and working directories and nothing else. Where
  home points elsewhere — `/home/you/data -> /scratch/you`, the usual cluster
  layout — it adds no mount for the target, cannot follow the link from inside,
  and falls back to `$HOME` with a warning that scrolls past. Every relative
  path in the command then resolved somewhere else, and the only symptom was
  `ERROR: Genome not found: your_genome.fa`. Two changes: a missing input inside
  a container now says so and names the fix, and
  [`bin/ltrquest-container`](bin/ltrquest-container) — shipped inside the image
  — works the mounts out from the command it is given.
- **The conda environment could not build the helpers it needs.**
  `environment.yml` had no `git`, `make` or compiler, but TRF-mod is compiled on
  first use and runs by default, so a clean environment died on `cc: not found`
  the first time it saw a genome. This was invisible on any machine with system
  build tools. The environment now carries them, and TRF-mod's build is retried
  with `CC=$CC` because its makefile hardcodes `CC=gcc`, which a conda toolchain
  does not provide.
- **GFF3 `long_terminal_repeat` sub-features landed on the wrong span for
  minus-strand elements.** This predates the rewrite: `gff3.py` applies
  `start + ltr5_end - 1` as a forward-genomic offset, but the pipeline used to
  measure those coordinates against a library where minus-strand records had
  already been reverse-complemented for classification, so the two
  sub-features were mis-sized whenever the LTRs differed in length or the
  flanks were asymmetric. It is fixed as a consequence of Kmer2LTR now
  running on the original forward-genomic sequence, ahead of classification,
  rather than of anything `gff3.py` itself had to change.
- **Same-round nest masking painted the wrong bases on reverse-complemented
  elements.** `nest_status` intervals are always forward-genomic, but on the
  ~45% of elements whose library record Kmer2LTR reverse-complemented for
  family clustering, the interval was never mirrored before being painted
  onto the flipped record. `mask_same_round_inners_in_fa` now mirrors it for
  every record `bounded_fasta_oriented` flipped.
- **Depth-bucketed nest masking had the same reverse-complement bug as the
  entry above.** `apply_depth_masking` is `mask_same_round_inners_in_fa`'s
  depth-level twin -- it paints the same kind of forward-genomic interval onto
  a library record -- but never carried the mirroring fix, so every
  minus-strand nest-outer got its depth mark at the mirror-image position in
  `{prefix}_depth<N>_ltr.fa`. It now reads the outer's own `orientation`
  column -- restated to match the library's flips earlier in the same round
  -- and mirrors for every row that reads `-`.
- **A truncated or extended element could be recorded as a nested insertion.**
  `_ltrs_shared` tells a truncation/extension variant from a genuine nested
  insertion by checking whether the two elements' LTRs overlap, and returns
  `False` ("distinct LTRs", i.e. a real nest) whenever either element's key is
  absent from the LTR-boundary map it is given. Reconcile's map came from
  `--scn`, keyed on each round's pre-trim locus, and was checked against
  elements Kmer2LTR had already renamed to their post-trim locus, so a missing
  key read as nesting instead of as a failed lookup -- disabling the
  truncation guard for 38% of elements (79/206 survivors on a full
  Arabidopsis chr2 run). `dedup_kmer2ltr_tsv` made the same call from a
  boundary map re-keyed through the rename rather than read from the
  post-trim table. Both now derive their map with `ltr_bounds_from_table`,
  which reads each element's `ltr5_start`/`ltr5_end`/`ltr3_start`/`ltr3_end`
  off its own row, so the key and the bounds it is checked against can no
  longer fall in different frames.

## [1.0.1] - 2026-08-24

Fixes two bugs that made the published 1.0.0 container unusable for a real run.
Both were found by pulling the released image and running it, rather than by
testing the source tree.

### Fixed

- **The container could not run a genome at all.** `ensure_tools` cloned and
  compiled minimap2 and miniprot from source *before* it would consider the
  copies already on `PATH`. The runtime image ships no `git` — deliberately, so
  a run never reaches the network — so detection died on
  `FileNotFoundError: 'git'` before it processed a single base. Resolution order
  is now: a prebuilt copy in `--tools-dir`, then `PATH`, then clone and build.
  Conda and pip users benefit too: they no longer pay a source build for tools
  their environment already provides.
- **`-profile singularity` could not pull its image.** Every Nextflow module
  pointed Singularity and Apptainer at `oras://…:1.0.0-singularity`, an artifact
  the release workflow never publishes, so the pull 404'd. Both engines convert
  an OCI image on the fly, so all engines now use the one published image.
- A missing `git` now raises a message naming the tool, the conda command that
  installs it, and the `--tools-dir` alternative, instead of a bare traceback
  from six frames inside `subprocess`.

### Added

- The container smoke test now resolves the helper binaries offline, which is
  the step `--help` never reached and the reason 1.0.0 shipped broken.
- Tests pinning the helper-resolution order, and guards against the package
  version drifting out of sync with `nextflow.config` or the module container
  tags, and against any module referencing an unpublished ORAS artifact.

### Changed

- README leads with Apptainer/Singularity rather than Docker, and gives
  copy-pasteable commands for the bundled test data on every install route.
  The audience is biologists on shared clusters, who cannot install Docker and
  do not have root.

## [1.0.0] - 2026-08-21

First release as a standalone tool. LTRquest was extracted from
[synLTR](https://github.com/cwb14/synLTR) `module2` with `git filter-repo`, so
the development history of every file it carries came with it.

### Added

- `ltrquest` console entry point, plus one per pipeline stage
  (`ltrquest-detect`, `-mask`, `-reconcile`, `-annotate`, `-gff3`, `-plots`).
- Installable Python package (`pip install .`), with the pipeline driver shipped
  as package data rather than a loose script.
- Nextflow DSL2 pipeline (`main.nf`) following nf-core conventions, with `stub:`
  blocks on every process so the whole DAG can be validated without the external
  toolchain.
- Container image with the helper repositories (Kmer2LTR, TEsorter2, TRF-mod,
  sdust) pre-built at pinned commits, so a run never reaches the network.
- `LTRQUEST_TOOLS_DIR` to point a run at a pre-built helper-tool directory.
- pytest suite covering the masking, reconciliation, annotation and GFF3 stages,
  and the installed package's wiring.
- GitHub Actions running the tests, lint, a Nextflow stub run, and a container
  build on every push.

### Changed

- Filenames collapsed to module names: `ltrharvest5.py` → `ltrquest.detect`,
  `mask_ltr.py` → `ltrquest.mask`, `reconcile_nests.py` → `ltrquest.reconcile`,
  `ltr_annotate.py` → `ltrquest.annotate`, `ltr_tsv_to_gff3.py` →
  `ltrquest.gff3`, `ltrharvest_wrapper2.sh` → the `ltrquest` driver.
- `--ltrharvest5-args` / `--ltrharvest5-args-from-round` renamed to
  `--detect-args` / `--detect-args-from-round`.
- `--script_path` removed. Stages are resolved as modules of the installed
  package, so there is nothing left to point at.
- GFF3 column 2 now reads `LTRquest` rather than `synLTR`.

### Removed

- The superseded detector and phasing revisions, the orphaned utilities, and the
  `solo_LTR` pipeline, none of which the entry point reached. They remain in the
  synLTR history.
