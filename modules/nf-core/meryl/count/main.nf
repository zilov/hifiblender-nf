process MERYL_COUNT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.4.1--h4ac6f70_0':
        'biocontainers/meryl:1.4.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(reads)
    val kvalue

    output:
    tuple val(meta), path("${meta.id}/*.meryl"), emit: meryl_db
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${meta.id}"
    """
    mkdir -p ${outdir}
    for READ in ${reads}; do
        meryl count \\
            k=${kvalue} \\
            threads=${task.cpus} \\
            memory=${task.memory.toGiga()} \\
            ${args} \\
            \$READ \\
            output ${outdir}/read.\${READ%.f*}.meryl
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${meta.id}"
    """
    mkdir -p ${outdir}
    for READ in ${reads}; do
        touch ${outdir}/read.\${READ%.f*}.meryl
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}