process LTRQUEST_MASK {
    tag "${meta.id}|r${round}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(genome), path(lib), path(prior_libs, stageAs: 'prior/*'), val(round)

    output:
    tuple val(meta), path("${prefix}_r${round}.fa"), emit: genome
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // The ORIGINAL genome is re-masked from scratch every round, never the
    // previous round's masked copy. Masking is cumulative in what it covers, but
    // the unmasked flanks have to come from real sequence or each round would
    // erode the context the next one needs to find a boundary.
    def round_char = ['N', 'R', 'D', 'Y', 'S', 'W', 'K', 'M', 'B', 'H'][round - 1]
    def extra      = round > 1 ? '--extra-features-fasta prior_libs.fa' : ''
    """
    ${round > 1 ? 'cat prior/* > prior_libs.fa' : ''}
    ltrquest-mask \\
        --features-fasta ${lib} \\
        --genome ${genome} \\
        --feature-character ${round_char} \\
        --far-character V \\
        --distance 15000 \\
        ${extra} \\
        ${args} \\
        > ${prefix}_r${round}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>chr1\\nACGTACGTACGTACGT\\n' > ${prefix}_r${round}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
