#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cwb14/LTRquest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/cwb14/LTRquest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { LTRQUEST } from './workflows/ltrquest'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
    log.info """
    LTRquest ${workflow.manifest.version}

    Usage:
      nextflow run cwb14/LTRquest --input samplesheet.csv --outdir results -profile docker

    Required:
      --input                 CSV samplesheet: sample,genome,proteins
      --outdir                Where to publish results

    Detection:
      --max_rounds            Detection rounds, 1-10 (default: ${params.max_rounds})
                              1 finds only un-nested elements, 2 adds one layer, ...
      --terminate_count       Stop once a round finds fewer than this many
                              elements (default: ${params.terminate_count})
      --fp_mask_threshold     False-positive element fraction above which the run
                              is flagged (default: ${params.fp_mask_threshold})
      --family_prefix         Names the shared family namespace (default: ${params.family_prefix})
      --skip_plots            Skip the plotting stage

    Profiles:
      -profile docker|singularity|conda|test|awsbatch
    """.stripIndent()
}

// Resolve a samplesheet field: absolute paths and URLs as given, anything else
// relative to the samplesheet's own directory.
def resolveInput(String path, sheetDir) {
    def remote = path ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/
    return (remote || path.startsWith('/'))
        ? file(path, checkIfExists: true)
        : file(sheetDir.resolve(path), checkIfExists: true)
}

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    if (!params.input)  { error "--input is required: a CSV samplesheet with columns sample,genome,proteins" }
    if (!params.outdir) { error "--outdir is required" }

    def n_rounds = params.max_rounds as int
    if (n_rounds < 1) { error "--max_rounds must be at least 1 (got ${params.max_rounds})" }

    //
    // Sequence IDs must be unique across genomes: pooled clustering keys
    // elements on 'chrom:start-end', so a shared ID would cross-assign families
    // and cross-purge real elements between genomes. The CLI checks this up
    // front; here it is the user's responsibility, and it is called out in
    // docs/nextflow.md.
    //
    // Relative paths in a samplesheet resolve against the samplesheet, not the
    // launch directory, so a sheet and its data stay portable as a unit.
    def sheet = file(params.input, checkIfExists: true)

    ch_samples = Channel
        .fromPath(sheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) { error "samplesheet row is missing a 'sample' value: ${row}" }
            if (!row.genome) { error "samplesheet row '${row.sample}' is missing a 'genome' path" }
            [ [ id: "${row.sample}_LTRs" ],
              resolveInput(row.genome, sheet.parent),
              row.proteins ? resolveInput(row.proteins, sheet.parent) : [] ]
        }

    LTRQUEST(ch_samples)
}
