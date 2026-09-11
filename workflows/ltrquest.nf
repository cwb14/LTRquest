/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LTRquest: iterative nested LTR-RT detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LTRQUEST_DETECT_ROUNDS } from '../subworkflows/local/ltrquest_detect_rounds/main'
include { LTRQUEST_RECONCILE     } from '../modules/local/ltrquest/reconcile/main'
include { LTRQUEST_CLUSTER       } from '../modules/local/ltrquest/cluster/main'
include { LTRQUEST_FLAGFP        } from '../modules/local/ltrquest/flagfp/main'
include { LTRQUEST_RECOVERSTRAND_ALIGN } from '../modules/local/ltrquest/recoverstrand/main'
include { LTRQUEST_RECOVERSTRAND_APPLY } from '../modules/local/ltrquest/recoverstrand/main'
include { LTRQUEST_ANNOTATE      } from '../modules/local/ltrquest/annotate/main'
include { LTRQUEST_GFF3          } from '../modules/local/ltrquest/gff3/main'
include { LTRQUEST_PLOTS         } from '../modules/local/ltrquest/plots/main'

workflow LTRQUEST {

    take:
    ch_samples  // channel: [ meta, genome, proteins ]

    main:
    // nextflow_schema.json is documentation, not a gate: nothing calls
    // validateParameters(). Check here so a typo'd preset fails now rather than
    // at the recovery stage, which sits after every expensive one. The CLI does
    // the same at ltrquest.sh's argument validation.
    if (params.strand_recovery
            && !(params.strand_recovery in ['conservative', 'balanced', 'sensitive'])) {
        error "--strand_recovery must be conservative, balanced or sensitive " +
              "(got '${params.strand_recovery}')"
    }
    if (params.strand_recovery_ppt && !params.strand_recovery) {
        error "--strand_recovery_ppt has no effect without --strand_recovery"
    }

    ch_versions = Channel.empty()

    //
    // Phase 1 - per genome: detect, mask, detect again, until the well runs dry.
    //
    LTRQUEST_DETECT_ROUNDS(ch_samples)
    ch_versions = ch_versions.mix(LTRQUEST_DETECT_ROUNDS.out.versions)

    //
    // Phase 2 - per genome: pool the rounds and resolve containment.
    //
    // Round N cannot know whether what it just found sits inside something round
    // N-1 found, so containment is settled once, here, over the union. The two
    // lists are sorted by round so their indices line up.
    //
    ch_reconcile_in = LTRQUEST_DETECT_ROUNDS.out.results
        .map { meta, round, tsv, fasta, workdir -> [ meta, [ round, tsv, fasta, workdir ] ] }
        .groupTuple()
        .map { meta, per_round ->
            def ordered = per_round.sort { a, b -> a[0] <=> b[0] }
            [ meta,
              ordered.collect { it[1] },
              ordered.collect { it[2] } ]
        }

    LTRQUEST_RECONCILE(ch_reconcile_in)
    ch_versions = ch_versions.mix(LTRQUEST_RECONCILE.out.versions)

    //
    // Phase 3 - pooled across every genome: one clustering pass, one
    // false-positive call.
    //
    // Pooling is the point. Families computed per species are not comparable
    // between species, and a repeat that looks convincing in one genome but
    // wrong across several only gets caught when they are judged together.
    //
    def pool = [ id: params.family_prefix ]

    ch_pooled_fasta = LTRQUEST_RECONCILE.out.fasta.map { _meta, files -> files }.collect()
    ch_pooled_tsv   = LTRQUEST_RECONCILE.out.tsv.map   { _meta, files -> files }.collect()

    LTRQUEST_CLUSTER(ch_pooled_fasta.map { files -> [ pool, files ] })
    ch_versions = ch_versions.mix(LTRQUEST_CLUSTER.out.versions)

    // ltrquest.flag_fp needs a genome to write its (here unused) masked FASTA.
    // The CLI would re-run the whole pipeline on that mask when false positives
    // are pervasive; see docs/nextflow.md for why this pipeline reports instead.
    ch_first_genome = ch_samples.map { _meta, genome, _proteins -> genome }.first()

    LTRQUEST_FLAGFP(
        LTRQUEST_CLUSTER.out.consensus_cluster
            .join(LTRQUEST_CLUSTER.out.internal_cluster)
            .join(LTRQUEST_CLUSTER.out.consensus_fasta)
            .combine(ch_pooled_tsv.map { files -> [files] })
            .combine(ch_pooled_fasta.map { files -> [files] })
            .combine(ch_first_genome)
    )
    ch_versions = ch_versions.mix(LTRQUEST_FLAGFP.out.versions)

    //
    // Phase 4 - back to per genome, judged against the one pooled cluster table,
    // so `family` denotes the same family in every genome's output.
    //
    // The FP stage returns one flat pile of cleaned files for every genome at
    // once; they are re-attributed by the prefix in their own names, which is
    // the same prefix the per-genome tables were written with.
    //
    ch_clean_by_id = LTRQUEST_FLAGFP.out.tsv.map { _meta, files -> files }.flatten()
        .mix(LTRQUEST_FLAGFP.out.fasta.map { _meta, files -> files }.flatten())
        .map { f -> [ f.name.replaceFirst(/_depth\d+_clean_ltr\.(tsv|fa)$/, ''), f ] }
        .groupTuple()

    ch_raw_by_id = LTRQUEST_RECONCILE.out.tsv.map { meta, files -> [ meta.id, files ] }
        .join(LTRQUEST_RECONCILE.out.fasta.map { meta, files -> [ meta.id, files ] })
        .map { id, tsvs, fastas -> [ id, tsvs + fastas ] }

    ch_work_by_id = LTRQUEST_DETECT_ROUNDS.out.results
        .map { meta, _round, _tsv, _fasta, workdir -> [ meta.id, workdir ] }
        .groupTuple()

    // One bundle per genome: everything the annotation, GFF3 and plotting stages
    // read, keyed by sample so each stage can take just the slice it needs.
    ch_bundle = ch_samples
        .map { meta, genome, _proteins -> [ meta.id, meta, genome ] }
        .join(ch_raw_by_id)
        .join(ch_clean_by_id)
        .join(ch_work_by_id)
        .combine(LTRQUEST_CLUSTER.out.consensus_cluster.map { _meta, tsv -> tsv })
        .map { _id, meta, genome, raw, clean, workdirs, cluster ->
            def files = raw + clean
            [ meta,
              files.findAll { it.name.endsWith('.tsv') },
              files.findAll { it.name.endsWith('.fa')  },
              workdirs, genome, cluster ]
        }

    // Strand recovery, first half. Only when asked for. It calls strand for the
    // elements the annotator's cascade cannot reach and writes the sidecar the
    // annotator reads back as a fourth tier. With recovery off, an empty list
    // flows in its place: it stages nothing and the modules render no flag for
    // it, so the run is exactly what it was before.
    if (params.strand_recovery) {
        LTRQUEST_RECOVERSTRAND_ALIGN(
            ch_bundle.map { meta, tsvs, _fastas, workdirs, genome, _cluster ->
                [ meta, tsvs, workdirs, genome ]
            }
        )
        ch_versions = ch_versions.mix(LTRQUEST_RECOVERSTRAND_ALIGN.out.versions)
        ch_recovery = LTRQUEST_RECOVERSTRAND_ALIGN.out.recovery
    } else {
        ch_recovery = ch_bundle.map { meta, _t, _f, _w, _g, _c -> [ meta, [] ] }
    }

    LTRQUEST_ANNOTATE(
        ch_bundle.map { meta, tsvs, _fastas, workdirs, _genome, cluster ->
            [ meta, tsvs, workdirs, cluster ]
        }.join(ch_recovery)
    )
    ch_versions = ch_versions.mix(LTRQUEST_ANNOTATE.out.versions)

    // Second half: bring the depth FASTAs into line with the strand column the
    // annotator just wrote, so a recovered minus element is stored in coding
    // sense exactly like a TEsorter2-called one. The GFF3 carries no sequence,
    // so it is indifferent to this and could equally have run first.
    if (params.strand_recovery) {
        LTRQUEST_RECOVERSTRAND_APPLY(
            LTRQUEST_ANNOTATE.out.tsv
                .join(ch_bundle.map { meta, _t, fastas, _w, _g, _c -> [ meta, fastas ] })
                .join(ch_recovery)
        )
        ch_versions = ch_versions.mix(LTRQUEST_RECOVERSTRAND_APPLY.out.versions)
        ch_tables = LTRQUEST_RECOVERSTRAND_APPLY.out.tsv
        ch_fastas = LTRQUEST_RECOVERSTRAND_APPLY.out.fasta
    } else {
        ch_tables = LTRQUEST_ANNOTATE.out.tsv
        ch_fastas = ch_bundle.map { meta, _t, fastas, _w, _g, _c -> [ meta, fastas ] }
    }

    // From here on the tables are the ANNOTATED ones, not the reconciler's.
    ch_annotated = ch_tables
        .join(ch_fastas)
        .join(ch_bundle.map { meta, _tsvs, _fastas, workdirs, genome, cluster ->
            [ meta, workdirs, genome, cluster ]
        })

    LTRQUEST_GFF3(
        ch_annotated.map { meta, tables, _fastas, workdirs, genome, cluster ->
            [ meta, tables, workdirs, genome, cluster ]
        }.join(ch_recovery)
    )
    ch_versions = ch_versions.mix(LTRQUEST_GFF3.out.versions)

    //
    // Phase 5 - plots.
    //
    if (!params.skip_plots) {
        LTRQUEST_PLOTS(
            ch_annotated.map { meta, tables, fastas, workdirs, genome, cluster ->
                [ meta, tables, fastas, workdirs, genome, cluster ]
            }
        )
        ch_versions = ch_versions.mix(LTRQUEST_PLOTS.out.versions)
    }

    emit:
    depth_tables = ch_tables
    gff3         = LTRQUEST_GFF3.out.gff3
    versions     = ch_versions
}
