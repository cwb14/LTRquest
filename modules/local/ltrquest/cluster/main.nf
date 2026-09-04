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
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    # Pool every depth, stripping the IUPAC characters the reconciler used to
    # mask nested inners, so each record is only the outermost element itself.
    cat ${depth_fastas} \\
      | awk '/^>/{printf "\\n%s\\n",\$0;next}{gsub(/[^ACGTacgt]/,"");printf "%s",\$0}END{print ""}' \\
      > ${prefix}_all_ltr.fa

    # Routed through ltrquest.kmer2ltr rather than a hardcoded script path:
    # Kmer2LTR ships as an installable package, and resolve() finds it whether
    # that means a console script, a prior clone, or a fresh one.
    #
    # The .tsv is load-bearing: Kmer2LTR names the cluster tables after
    # `Path(-o).with_suffix("")`, which strips whatever follows the last dot.
    # A sample id carrying one -- GCF_000001735.4_TAIR10.1 and every other
    # assembly name -- would otherwise put them under a truncated basename that
    # the output patterns above never match.
    python -c '
import sys
from ltrquest.kmer2ltr import resolve, run
tools_dir, in_fa, out_tsv, threads, mutation_rate, min_seq_id = sys.argv[1:7]
run(resolve(tools_dir), in_fa, out_tsv,
    threads=int(threads), mutation_rate=float(mutation_rate),
    ltr_cluster=True, internal_cluster=True, min_seq_id=float(min_seq_id), verbose=True)
' \${LTRQUEST_TOOLS_DIR:-/opt/ltrquest/tools} ${prefix}_all_ltr.fa ${prefix}_all_ltr.tsv ${task.cpus} \\
  ${params.mutation_rate} ${params.min_seq_id}

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
