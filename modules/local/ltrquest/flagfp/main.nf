process LTRQUEST_FLAGFP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/cwb14/ltrquest:1.0.0-singularity' :
        'ghcr.io/cwb14/ltrquest:1.0.0' }"

    input:
    tuple val(meta), path(consensus_cluster), path(internal_cluster), path(consensus_fasta), path(depth_tsvs), path(depth_fastas), path(genome)

    output:
    tuple val(meta), path("*_depth*_clean_ltr.tsv", arity: '1..*'), emit: tsv
    tuple val(meta), path("*_depth*_clean_ltr.fa", arity: '1..*') , emit: fasta
    tuple val(meta), path("${prefix}_fpcheck.log") , emit: log
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # High-abundance non-LTR repeats can seed convincing but spurious LTR-RT
    # calls. Families whose members are mostly such repeats are purged here.
    python \${LTRQUEST_TOOLS_DIR:-/opt/ltrquest/tools}/Kmer2LTR/flag_fp_families.py \\
        --consensus-cluster ${consensus_cluster} \\
        --internal-cluster ${internal_cluster} \\
        --ltr-fasta ${consensus_fasta} \\
        --domains-tsv ${depth_tsvs} \\
        -o ${prefix}_fpcheck \\
        --genome ${genome} \\
        --masked-out ${prefix}_FP_masked.fa \\
        --threads ${task.cpus} \\
        --fp-mask-threshold ${params.fp_mask_threshold} \\
        ${args} \\
        2>&1 | tee ${prefix}_fpcheck.log >&2

    # Keep only the records whose id survived into the cleaned table.
    for tsv in ${depth_tsvs}; do
        clean="\${tsv%_ltr.tsv}_clean_ltr.tsv"
        fa="\${tsv%_ltr.tsv}_ltr.fa"
        [ -s "\$clean" ] && [ -s "\$fa" ] || continue
        awk -F'\\t' 'NR==FNR { if (\$0 !~ /^#/ && NF >= 2) keep[\$1]=1; next }
                     /^>/ { p = (substr(\$0, 2) in keep) } p' \\
            "\$clean" "\$fa" > "\${tsv%_ltr.tsv}_clean_ltr.fa"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Names are derived from the input tables, exactly as the real stage does:
    // downstream joins key on the per-genome prefix embedded in those names, so
    // a stub that invented its own would test the wrong wiring.
    """
    for tsv in ${depth_tsvs}; do
        cp "\$tsv" "\${tsv%_ltr.tsv}_clean_ltr.tsv"
    done
    for fa in ${depth_fastas}; do
        cp "\$fa" "\${fa%_ltr.fa}_clean_ltr.fa"
    done
    printf '[INFO] FP fraction: 0/1 = 0.0000 (threshold ${params.fp_mask_threshold})\\n' > ${prefix}_fpcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: 16.747c6
    END_VERSIONS
    """
}
