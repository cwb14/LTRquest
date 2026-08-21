process LTRQUEST_PLOTS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/cwb14/ltrquest:1.0.0-singularity' :
        'ghcr.io/cwb14/ltrquest:1.0.0' }"

    input:
    tuple val(meta), path(tables), path(fastas), path(workdirs), path(genome), path(consensus_cluster)

    output:
    tuple val(meta), path("${prefix}_plots"), emit: plots
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // Plots are a nice-to-have. `error_ignore` keeps a matplotlib hiccup from
    // discarding an otherwise finished annotation.
    """
    ltrquest-plots \\
        --prefix ${prefix} \\
        --genome ${genome} \\
        --indir . \\
        --out-dir ${prefix}_plots \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_plots/struct
    touch ${prefix}_plots/${prefix}_summary.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matplotlib: 3.8.0
    END_VERSIONS
    """
}
