process LTRQUEST_CLUSTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(depth_fastas)

    output:
    tuple val(meta), path("${prefix}_all_ltr.fa")                        , emit: pooled
    tuple val(meta), path("${prefix}_all_ltr.consensus_id*_cluster.tsv") , emit: consensus_cluster
    tuple val(meta), path("${prefix}_all_ltr.internal_id*_cluster.tsv")  , emit: internal_cluster
    tuple val(meta), path("${prefix}_all_ltr.consensus.fa")              , emit: consensus_fasta
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: '--min-seq-id 0.75'
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    # Pool every depth, stripping the IUPAC characters the reconciler used to
    # mask nested inners, so each record is only the outermost element itself.
    cat ${depth_fastas} \\
      | awk '/^>/{printf "\\n%s\\n",\$0;next}{gsub(/[^ACGTacgt]/,"");printf "%s",\$0}END{print ""}' \\
      > ${prefix}_all_ltr.fa

    python \${LTRQUEST_TOOLS_DIR:-/opt/ltrquest/tools}/Kmer2LTR/Kmer2LTR.py \\
        -i ${prefix}_all_ltr.fa \\
        -o ${prefix}_all_ltr \\
        --ltr-cluster --internal-cluster \\
        -p ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>chr1:1000-6000#LTR/Gypsy/Tekay\\nACGT\\n' > ${prefix}_all_ltr.fa
    printf 'chr1:1000-6000\\tchr1:1000-6000\\n' > ${prefix}_all_ltr.consensus_id0.75_cluster.tsv
    printf 'chr1:1000-6000\\tchr1:1000-6000\\n' > ${prefix}_all_ltr.internal_id0.75_cluster.tsv
    printf '>chr1:1000-6000\\nACGT\\n' > ${prefix}_all_ltr.consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: 16.747c6
    END_VERSIONS
    """
}
