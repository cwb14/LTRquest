/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LTRquest: iterative nested LTR-RT detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LTRQUEST_DETECT_ROUNDS } from '../subworkflows/local/ltrquest_detect_rounds/main'
include { LTRQUEST_RECONCILE     } from '../modules/local/ltrquest/reconcile/main'
include { LTRQUEST_CLUSTER       } from '../modules/local/ltrquest/cluster/main'
include { LTRQUEST_FLAGFP        } from '../modules/local/ltrquest/flagfp/main'
include { LTRQUEST_ANNOTATE      } from '../modules/local/ltrquest/annotate/main'
include { LTRQUEST_GFF3          } from '../modules/local/ltrquest/gff3/main'
include { LTRQUEST_PLOTS         } from '../modules/local/ltrquest/plots/main'

workflow LTRQUEST {

    take:
    ch_samples  // channel: [ meta, genome, proteins ]

    main:
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
    // N-1 found, so containment is settled once, here, over the union. The three
    // lists are sorted by round so their indices line up.
    //
    ch_reconcile_in = LTRQUEST_DETECT_ROUNDS.out.results
        .map { meta, round, tsv, fasta, workdir -> [ meta, [ round, tsv, fasta, workdir ] ] }
        .groupTuple()
        .map { meta, per_round ->
            def ordered = per_round.sort { a, b -> a[0] <=> b[0] }
            [ meta,
              ordered.collect { it[1] },
              ordered.collect { it[2] },
              ordered.collect { it[3] }.collect { dir -> file("${dir}/*.ltrtools.stitched.scn") }.flatten() ]
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

    LTRQUEST_ANNOTATE(
        ch_bundle.map { meta, tsvs, _fastas, workdirs, _genome, cluster ->
            [ meta, tsvs, workdirs, cluster ]
        }
    )
    ch_versions = ch_versions.mix(LTRQUEST_ANNOTATE.out.versions)

    // From here on the tables are the ANNOTATED ones, not the reconciler's.
    ch_annotated = LTRQUEST_ANNOTATE.out.tsv
        .join(ch_bundle.map { meta, _tsvs, fastas, workdirs, genome, cluster ->
            [ meta, fastas, workdirs, genome, cluster ]
        })

    LTRQUEST_GFF3(
        ch_annotated.map { meta, tables, _fastas, workdirs, genome, cluster ->
            [ meta, tables, workdirs, genome, cluster ]
        }
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
    depth_tables = LTRQUEST_ANNOTATE.out.tsv
    gff3         = LTRQUEST_GFF3.out.gff3
    versions     = ch_versions
}
