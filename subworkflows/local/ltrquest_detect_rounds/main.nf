//
// The iterative detection loop: scan, mask what was found, scan again.
//
include { LTRQUEST_ROUND as ROUND_01 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_02 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_03 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_04 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_05 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_06 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_07 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_08 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_09 } from '../ltrquest_round/main'
include { LTRQUEST_ROUND as ROUND_10 } from '../ltrquest_round/main'

// Starve every round past --max_rounds. Written as a function rather than a
// closure because the language spec only resolves call syntax against
// declarations.
def gate(ch, round) {
    return (params.max_rounds as int) >= round ? ch : Channel.empty()
}

workflow LTRQUEST_DETECT_ROUNDS {

    take:
    ch_samples  // channel: [ meta, genome, proteins ]

    main:
    //
    // Why the chain is written out rather than looped.
    //
    // Round N scans the genome round N-1 masked, so the rounds are a strict
    // chain, and Nextflow gives no way to express that compactly:
    //
    //   * a process or workflow may be invoked only once per workflow, so a
    //     Groovy loop over one ROUND is out;
    //   * calling a workflow from inside a closure is rejected outright, so a
    //     list of round-closures is out;
    //   * the `recurse` operator, which exists for exactly this, is still behind
    //     `nextflow.preview.recursion` and takes only value channels -- which
    //     would forfeit the per-sample parallelism that is the reason to be
    //     running this under Nextflow at all.
    //
    // So the round subworkflow is included once per possible round and the chain
    // is written down. Ten is the ceiling because the detector has ten IUPAC
    // characters to mark rounds with, and the reconciler reads a nested
    // element's round back off that character.
    //
    // The unroll costs nothing it does not use. `gate` starves any round past
    // --max_rounds, and a round that finds fewer than --terminate_count elements
    // emits no successor state, so every later round sees an empty channel and
    // is never scheduled. Nextflow's laziness is the loop's `break`.
    //
    def max_rounds = params.max_rounds as int

    ROUND_01( ch_samples.map { meta, genome, proteins -> [ meta, genome, genome, proteins, [] ] },  1 )
    ROUND_02( gate(ROUND_01.out.next,  2),  2 )
    ROUND_03( gate(ROUND_02.out.next,  3),  3 )
    ROUND_04( gate(ROUND_03.out.next,  4),  4 )
    ROUND_05( gate(ROUND_04.out.next,  5),  5 )
    ROUND_06( gate(ROUND_05.out.next,  6),  6 )
    ROUND_07( gate(ROUND_06.out.next,  7),  7 )
    ROUND_08( gate(ROUND_07.out.next,  8),  8 )
    ROUND_09( gate(ROUND_08.out.next,  9),  9 )
    ROUND_10( gate(ROUND_09.out.next, 10), 10 )

    if (max_rounds > 10) {
        log.warn "--max_rounds ${max_rounds} exceeds the 10 available IUPAC mask characters; capping at 10."
    }

    ch_results = ROUND_01.out.results.mix(
        ROUND_02.out.results, ROUND_03.out.results, ROUND_04.out.results,
        ROUND_05.out.results, ROUND_06.out.results, ROUND_07.out.results,
        ROUND_08.out.results, ROUND_09.out.results, ROUND_10.out.results)

    ch_genic = ROUND_01.out.genic.mix(
        ROUND_02.out.genic, ROUND_03.out.genic, ROUND_04.out.genic,
        ROUND_05.out.genic, ROUND_06.out.genic, ROUND_07.out.genic,
        ROUND_08.out.genic, ROUND_09.out.genic, ROUND_10.out.genic)

    ch_versions = ROUND_01.out.versions.mix(
        ROUND_02.out.versions, ROUND_03.out.versions, ROUND_04.out.versions,
        ROUND_05.out.versions, ROUND_06.out.versions, ROUND_07.out.versions,
        ROUND_08.out.versions, ROUND_09.out.versions, ROUND_10.out.versions)

    emit:
    results  = ch_results   // [ meta, round, tsv, fa, workdir ]
    genic    = ch_genic     // [ meta, gff ]
    versions = ch_versions
}
