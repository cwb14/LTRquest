// Strand recovery, in the two halves it has to run in.
//
// ALIGN calls strand for the elements ltrquest.annotate's cascade cannot reach
// and publishes the sidecar the annotator reads as its fourth tier. APPLY then
// re-orients the depth FASTAs to the strand column the annotator wrote, so the
// stored orientation and the table agree by construction. The annotator runs
// between them, which is why this is two processes and not one.

process LTRQUEST_RECOVERSTRAND_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(tables), path(workdirs), path(genome)

    output:
    tuple val(meta), path("${prefix}_strand_recovery.tsv"), emit: recovery
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ppt  = params.strand_recovery_ppt ? '--ppt' : ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // Read-only: this phase writes the sidecar and nothing else, so the tables
    // do not have to be copied out of the staging directory the way the
    // annotator's do.
    """
    ltrquest-recover-strand \\
        --phase align \\
        --prefix ${prefix} \\
        --indir . \\
        --genome ${genome} \\
        --preset ${params.strand_recovery} \\
        --threads ${task.cpus} \\
        ${ppt} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
        blast: \$(blastn -version 2>&1 | sed -n 's/^blastn: //p')
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '#locus\\tstrand\\tsource\\tevidence\\n' > ${prefix}_strand_recovery.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}


process LTRQUEST_RECOVERSTRAND_APPLY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(tables, stageAs: 'in/*'), path(fastas, stageAs: 'in/*'), path(recovery)

    output:
    tuple val(meta), path("*_depth*_ltr.tsv", arity: '1..*'), emit: tsv
    tuple val(meta), path("*_depth*_ltr.fa" , arity: '1..*'), emit: fasta
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // Both the tables and the FASTAs are rewritten in place, so they are copied
    // out of the staging directory first: Nextflow stages inputs as symlinks,
    // and editing through one would corrupt the upstream task's published
    // outputs. LTRQUEST_ANNOTATE does the same for the same reason.
    """
    cp in/*.tsv in/*.fa .

    ltrquest-recover-strand \\
        --phase apply \\
        --prefix ${prefix} \\
        --indir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp in/*.tsv in/*.fa .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
