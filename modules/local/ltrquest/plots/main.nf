process LTRQUEST_PLOTS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

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
