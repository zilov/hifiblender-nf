process MERFIN_HIST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merfin:1.0--h4ac6f70_2':
        'biocontainers/merfin:1.0--h4ac6f70_2' }"

    input:
    tuple val(meta), path(fasta_assembly)   // Required Input -sequence files can be FASTA or FASTQ; uncompressed, gz compressed.
    tuple val(meta1), path(meryl_db_reads)  // Required readmers (raw reads meryl db). As it comes from another tool, it might be relevant to mantain the meta.
    path(lookup_table)                      // Optional input vector of probabilities (obtained by genomescope2 with parameter --fitted_hist).
    path(seqmers)                           // Optional input for pre-built sequence meryl db (-seqmers).
    val(peak)                               // Required input to hard set copy 1 and infer multiplicity to copy number.

    output:
    tuple val(meta), path("*.hist")     , emit: hist
    tuple val(meta), path("*_qv_results.tsv"), emit: qv_results
    path("*.hist.stderr.log")           , emit: log_stderr
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def prefix                  = task.ext.prefix ?: "${meta.id}"
    def optional_lookup_table   = lookup_table ? "-prob ${lookup_table}" : ""
    def optional_seqmers        = seqmers ? "-seqmers ${seqmers}" : ""
    def mem_adjust              = "${task.memory}".replaceAll(' ', '')  // Memory parameter - essential to avoid crash
    """
    merfin -hist \\
        -threads $task.cpus \\
        -memory $mem_adjust \\
        $args \\
        -sequence $fasta_assembly \\
        -readmers $meryl_db_reads \\
        -peak $peak \\
        $optional_lookup_table \\
        $optional_seqmers \\
        -output ${prefix}.hist \\
        2> >( tee ${prefix}.hist.stderr.log )

    # Parse the stderr log and extract the desired parameters
    awk '
        /K-mers not found in reads \\(missing\\)/ {missing=\$NF}
        /K-mers overly represented in assembly/ {overrep=\$NF}
        /K-mers found in the assembly/ {found=\$NF}
        /Missing QV:/ {missing_qv=\$NF}
        /Merfin QV\\*:/ {merfin_qv=\$NF}
        END {
            print "Sample\t${prefix}" > "${prefix}_qv_results.tsv"
            print "Parameter\tValue" >> "${prefix}_qv_results.tsv"
            print "K-mers not found in reads (missing)\t" missing >> "${prefix}_qv_results.tsv"
            print "K-mers overly represented in assembly\t" overrep >> "${prefix}_qv_results.tsv"
            print "K-mers found in the assembly\t" found >> "${prefix}_qv_results.tsv"
            print "Missing QV\t" missing_qv >> "${prefix}_qv_results.tsv"
            print "Merfin QV*\t" merfin_qv >> "${prefix}_qv_results.tsv"
        }
    ' ${prefix}.hist.stderr.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merfin: \$( merfin --version |& sed 's/merfin //' )
    END_VERSIONS
    """

    stub:
    def args                    = task.ext.args ?: ''
    def prefix                  = task.ext.prefix ?: "${meta.id}"
    def optional_lookup_table   = lookup_table ? "-prob ${lookup_table}" : ""
    def optional_seqmers        = seqmers ? "-seqmers ${seqmers}" : ""
    """
    touch ${prefix}.hist
    touch ${prefix}.hist.stderr.log
    touch ${prefix}_qv_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merfin: \$( merfin --version |& sed 's/merfin //' )
    END_VERSIONS
    """
}