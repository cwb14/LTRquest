// Family-guided re-boundarying. One pooled task: a family's LTR model is built
// from its copies in every genome, so the samples cannot be split here the way
// the per-genome stages are. It rewrites the _clean_ tables and FASTAs and
// writes one <prefix>_reboundary.tsv per sample, which LTRQUEST_GFF3 reads to
// follow the elements it renamed.

process LTRQUEST_REBOUNDARY {
    tag "pooled"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    val(ids)
    path(tables, stageAs: 'in/*')
    path(fastas, stageAs: 'in/*')
    path(genomes, stageAs: 'genomes/g??/*')

    output:
    path("*_depth*_clean_ltr.tsv", arity: '1..*'), emit: tsv
    path("*_depth*_clean_ltr.fa" , arity: '1..*'), emit: fasta
    path("*_reboundary.tsv"      , arity: '1..*'), emit: sidecar
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_list = (genomes instanceof List ? genomes : [genomes]).join(' ')
    // Rewritten in place, so copied out of the staging directory first: editing
    // through Nextflow's input symlinks would corrupt the upstream task's outputs.
    """
    cp in/*_clean_ltr.tsv in/*_clean_ltr.fa .

    ltrquest-reboundary \\
        --indir . \\
        --prefix ${ids.join(' ')} \\
        --genome ${genome_list} \\
        --threads ${task.cpus} \\
        --mutation-rate ${params.mutation_rate} \\
        --tools-dir \${LTRQUEST_TOOLS_DIR:-/opt/ltrquest/tools} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
        mafft: \$(mafft --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    cp in/*_clean_ltr.tsv in/*_clean_ltr.fa .
    for id in ${ids.join(' ')}; do
      printf '#old_seq_id\\tnew_seq_id\\tdecision\\n' > \${id}_reboundary.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
