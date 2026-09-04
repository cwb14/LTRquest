process LTRQUEST_DETECT {
    tag "${meta.id}|r${round}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // One image for every engine. Singularity and Apptainer convert an OCI
    // image on the fly, so the nf-core habit of pointing them at a separate
    // `oras://…-singularity` artifact only helps if you actually publish one.
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    tuple val(meta), path(genome), path(proteins), path(prior_libs, stageAs: 'prior/*'), val(round)

    output:
    tuple val(meta), val(round), path("${prefix}_ltr.tsv") , emit: tsv
    tuple val(meta), val(round), path("${prefix}_ltr.fa")  , emit: fasta
    tuple val(meta), val(round), path("${prefix}.work")    , emit: workdir
    tuple val(meta), path("${prefix}.work/*.genic.gff")    , emit: genic, optional: true
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // The round index is copied into a local before `prefix` is assigned.
    // `prefix` must be an undeclared assignment so the output block can see it,
    // and in Nextflow 24.x an undeclared assignment that interpolates an INPUT
    // variable re-declares that input in the process scope -- after which every
    // later reference to it fails to compile with "already defined".
    def round_n = round as int
    prefix   = task.ext.prefix ?: "${meta.id}_r${round_n}"

    // Detection widens its search window each round, because a nested element's
    // host spans more sequence than the element itself: masking round N-1's hits
    // leaves a run of IUPAC characters that the round-N candidate has to cover.
    def scn_max_ret = 150000 + (round_n - 1) * 15000
    def scn_max_int = 140000 * round_n
    def maxdistltr  = 15000  + (round_n - 1) * 15000
    def ltrf_d      = 15000  + (round_n - 1) * 15000
    def overlap     = 25000  + (round_n - 1) * 15000

    // Round N paints its own same-round inners with IUPAC_SEQ[N-1]. From round 2
    // on, a candidate must contain at least one earlier round's character (it is
    // only interesting if something is nested in it) and must not contain 'V',
    // which marks sequence too far from any element to be worth revisiting.
    def iupac      = ['N', 'R', 'D', 'Y', 'S', 'W', 'K', 'M', 'B', 'H']
    def round_char = iupac[round_n - 1]
    def require    = round_n > 1 ? iupac[0..(round_n - 2)].join(',') : ''

    def protein_opt = proteins ? "--proteins ${proteins}" : ''
    def pass2_opt   = round_n > 1
        ? "--pass2-classified-fasta pass2_lib.fa --require-run-chars ${require} --exclude-run-char V"
        : ''
    """
    ${round_n > 1 ? """
    # Pass-2 reference library: every prior round's elements, reduced to A/C/G/T
    # so the IUPAC nest markers do not leak into the homology search.
    cat prior/* \\
      | awk '/^>/ {printf("\\n%s\\n",\$0);next;} {printf("%s",\$0);} END {printf("\\n");}' \\
      | sed '/^>/! s/[^ATCGatcg]//g' > pass2_lib.fa
    """ : ''}
    ltrquest-detect \\
        --genome ${genome} \\
        ${protein_opt} \\
        --threads ${task.cpus} \\
        --out-prefix ${prefix} \\
        --tools-dir \${LTRQUEST_TOOLS_DIR:-./tools} \\
        --scn-min-ltr-len 10 \\
        --scn-min-ret-len 80 \\
        --scn-max-ret-len ${scn_max_ret} \\
        --scn-min-int-len 0 \\
        --scn-max-int-len ${scn_max_int} \\
        --ltrharvest-args "-mindistltr 100 -minlenltr 100 -maxlenltr 7000 -mintsd 0 -maxtsd 0 -similar 70 -vic 60 -seed 15 -seqids yes -xdrop 10 -maxdistltr ${maxdistltr}" \\
        --ltrfinder-args "-w 2 -C -D ${ltrf_d} -d 100 -L 7000 -l 100 -p 20 -M 0.00 -S 0.0" \\
        --size 500000 \\
        --overlap ${overlap} \\
        --tesorter-rule 70-75-80 \\
        --tsd-pass2 \\
        --nested-flank-min 10 \\
        --nested-base-min 800 \\
        --same-round-inner-char ${round_char} \\
        --mutation-rate ${params.mutation_rate} \\
        ${pass2_opt} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
        genometools: \$(gt --version 2>&1 | head -n1 | sed 's/^gt (GenomeTools) //')
        ltr_finder: \$(ltr_finder -h 2>&1 | grep -oE 'v[0-9.]+' | head -n1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    def round_n = round as int
    prefix = task.ext.prefix ?: "${meta.id}_r${round_n}"
    // Two elements in round 1, one in round 2, none from round 3 on -- enough to
    // exercise the terminate-count gate and the cross-round nesting path.
    def n = round_n == 1 ? 2 : (round_n == 2 ? 1 : 0)
    """
    mkdir -p ${prefix}.work
    touch ${prefix}.work/${prefix}.ltrtools.stitched.scn
    printf '#seq_id\\tseq_len\\tstatus\\tltr5_start\\tltr5_end\\tltr3_start\\tltr3_end\\tltr5_len\\tltr3_len\\tflank5_len\\tflank3_len\\taln_len\\tn_sites\\tn_ts\\tn_tv\\tn_gapcols\\tidentity\\tp_dist\\tk2p\\tk2p_se\\tbitscore\\tflank_margin_bits\\tcigar\\tmotif\\tk2p_time\\torientation\\ttsd\\ttsd_offset\\ttsd_input\\tdomains\\tnest_status\\n' > ${prefix}_ltr.tsv
    : > ${prefix}_ltr.fa
    # One field per header name: a stub row narrower than its own header is a
    # schema mismatch that -stub-run, the only automated exercise this DAG
    # gets, would otherwise read straight past.
    for i in \$(seq 1 ${n}); do
        s=\$(( ${round_n} * 1000 + i * 100 ))
        e=\$(( s + 5000 ))
        g1=\$(( s + 600 ));  g2=\$(( s + 1200 ))
        r1=\$(( s + 1500 )); r2=\$(( s + 2400 ))
        printf 'chr1:%s-%s#LTR/Gypsy/Tekay\\t5001\\tpass\\t1\\t500\\t4502\\t5001\\t500\\t500\\t0\\t0\\t500\\t500\\t8\\t4\\t0\\t0.9760\\t0.0240\\t0.0244\\t0.0070\\t900.0\\t12.5\\t500=\\ttg...ca\\t406667\\t+\\tTGCAA\\t0,0\\tTGCAA\\tGAG|Tekay@%s-%s;RT|Tekay@%s-%s\\t.\\n' \\
            "\$s" "\$e" "\$g1" "\$g2" "\$r1" "\$r2" >> ${prefix}_ltr.tsv
        printf '>chr1:%s-%s#LTR/Gypsy/Tekay\\nACGTACGTAC\\n' "\$s" "\$e" >> ${prefix}_ltr.fa
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
        genometools: 1.6.6
        ltr_finder: 1.07
    END_VERSIONS
    """
}
