process LTRQUEST_GFF3 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(tables), path(workdirs), path(genome), path(consensus_cluster)

    output:
    tuple val(meta), path("${prefix}_all_depth_LTR_cleaned.gff3")        , emit: gff3
    tuple val(meta), path("${prefix}_all_depth_protein_LTR_cleaned.gff3"), emit: protein_gff3, optional: true
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // A second GFF3 carrying the miniprot alignments is written too, but only
    // when round 1 produced a genic GFF -- that is, only when the run had
    // --proteins. Hence `optional: true` above.
    """
    ltrquest-gff3 \\
        --prefix ${prefix} \\
        --indir . \\
        --genome ${genome} \\
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
    printf '##gff-version 3\\n' > ${prefix}_all_depth_LTR_cleaned.gff3
    printf 'chr1\\tLTRquest\\tLTR_retrotransposon\\t1000\\t6000\\t.\\t-\\t.\\tID=${prefix}_LTRRT_00001\\n' \\
      >> ${prefix}_all_depth_LTR_cleaned.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
