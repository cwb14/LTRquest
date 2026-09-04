process LTRQUEST_RECONCILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(tsvs, stageAs: 'round_*/*'), path(fastas, stageAs: 'round_*/*')

    output:
    tuple val(meta), path("${prefix}_depth*_ltr.tsv", arity: '1..*'), emit: tsv
    tuple val(meta), path("${prefix}_depth*_ltr.fa", arity: '1..*') , emit: fasta
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    // The two lists are round-ordered and index-aligned; the reconciler reads
    // the round index off the position, so they must be sorted identically.
    def tsv_list = tsvs.collect  { it.toString() }.sort().join(' ')
    def fa_list  = fastas.collect{ it.toString() }.sort().join(' ')
    """
    ltrquest-reconcile \\
        --out-prefix ${prefix} \\
        --tsv ${tsv_list} \\
        --fa ${fa_list} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '#seq_id\\tseq_len\\tstatus\\tltr5_start\\tltr5_end\\tltr3_start\\tltr3_end\\tltr5_len\\tltr3_len\\tflank5_len\\tflank3_len\\taln_len\\tn_sites\\tn_ts\\tn_tv\\tn_gapcols\\tidentity\\tp_dist\\tk2p\\tk2p_se\\tbitscore\\tflank_margin_bits\\tcigar\\tmotif\\tk2p_time\\torientation\\ttsd\\ttsd_offset\\ttsd_input\\tdomains\\tnest_status\\n' > ${prefix}_depth0_ltr.tsv
    # One field per header name: a stub row narrower than its own header is a
    # schema mismatch that -stub-run, the only automated exercise this DAG
    # gets, would otherwise read straight past.
    printf 'chr1:1000-6000#LTR/Gypsy/Tekay\\t5001\\tpass\\t1\\t500\\t4502\\t5001\\t500\\t500\\t0\\t0\\t500\\t500\\t8\\t4\\t0\\t0.9760\\t0.0240\\t0.0244\\t0.0070\\t900.0\\t12.5\\t500=\\ttg...ca\\t406667\\t+\\tTGCAA\\t0,0\\tTGCAA\\tGAG|Tekay@1600-2200;RT|Tekay@2500-3400\\t.\\n' >> ${prefix}_depth0_ltr.tsv
    printf '>chr1:1000-6000#LTR/Gypsy/Tekay\\nACGT\\n' > ${prefix}_depth0_ltr.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
