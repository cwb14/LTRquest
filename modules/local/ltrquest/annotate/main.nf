process LTRQUEST_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(tables, stageAs: 'in/*'), path(workdirs), path(consensus_cluster)

    output:
    tuple val(meta), path("*_depth*_ltr.tsv", arity: '1..*'), emit: tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // The annotator rewrites its tables in place, so they are copied out of the
    // staging directory first: Nextflow stages inputs as symlinks, and editing
    // through one would corrupt the upstream task's published outputs.
    //
    // Every genome is annotated against the SAME pooled cluster table, which is
    // what makes `family` mean the same thing across genomes.
    """
    cp in/*.tsv .

    ltrquest-annotate \\
        --prefix ${prefix} \\
        --indir . \\
        --consensus-cluster ${consensus_cluster} \\
        --family-prefix ${params.family_prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp in/*.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
