//
// One detection round: find elements in the current (progressively masked)
// genome, then mask them to build the genome the next round will scan.
//
include { LTRQUEST_DETECT } from '../../../modules/local/ltrquest/detect/main'
include { LTRQUEST_MASK   } from '../../../modules/local/ltrquest/mask/main'

workflow LTRQUEST_ROUND {

    take:
    ch_state  // channel: [ meta, round_genome, original_genome, proteins, prior_libs ]
    round     // integer:  which round this is, 1-based

    main:
    ch_versions = Channel.empty()

    LTRQUEST_DETECT(
        ch_state.map { meta, round_genome, _orig, proteins, prior_libs ->
            [ meta, round_genome, proteins, prior_libs, round ]
        }
    )
    ch_versions = ch_versions.mix(LTRQUEST_DETECT.out.versions.first())

    // The three artifacts of a round travel together from here on: the
    // reconciler needs all of them, index-aligned by round.
    ch_results = LTRQUEST_DETECT.out.tsv
        .join(LTRQUEST_DETECT.out.fasta,   by: [0, 1])
        .join(LTRQUEST_DETECT.out.workdir, by: [0, 1])

    //
    // Termination gate. A round that turns up fewer than --terminate_count
    // elements is the last one worth running: it publishes its hits but no
    // successor state, so the next round's input channel stays empty and the
    // chain simply stops advancing.
    //
    // The last permitted round short-circuits the gate entirely. Masking a
    // genome is not free, and on the final round there is nothing left to mask
    // it for.
    //
    ch_continue = round >= (params.max_rounds as int)
        ? Channel.empty()
        : LTRQUEST_DETECT.out.fasta
            .filter { _meta, _round, fasta -> fasta.countFasta() >= params.terminate_count }
            .map { meta, _round, fasta -> [ meta, fasta ] }

    // Round N masks the ORIGINAL genome using its own hits plus every earlier
    // round's, so coverage accumulates while the unmasked flanks stay real
    // sequence.
    LTRQUEST_MASK(
        ch_state
            .map { meta, _round_genome, orig, _proteins, prior_libs -> [ meta, orig, prior_libs ] }
            .join(ch_continue)
            .map { meta, orig, prior_libs, lib -> [ meta, orig, lib, prior_libs, round ] }
    )
    ch_versions = ch_versions.mix(LTRQUEST_MASK.out.versions.first())

    // State for round N+1: the freshly masked genome, plus this round's library
    // appended to the running list used for pass-2 homology and for masking.
    ch_next = ch_state
        .map { meta, _round_genome, orig, proteins, prior_libs ->
            [ meta, orig, proteins, prior_libs ]
        }
        .join(LTRQUEST_MASK.out.genome)
        .join(ch_continue)
        .map { meta, orig, proteins, prior_libs, masked, lib ->
            [ meta, masked, orig, proteins, prior_libs + [lib] ]
        }

    emit:
    results  = ch_results                  // [ meta, round, tsv, fa, workdir ]
    genic    = LTRQUEST_DETECT.out.genic   // [ meta, gff ] - round 1 only, needs --proteins
    next     = ch_next                     // [ meta, masked_genome, orig, proteins, libs ]
    versions = ch_versions
}
