process CREATE_NEXTDENOVO_CONFIG {
    label 'process_assembly'
    
    input:
    val(meta)
    path(genomescope_summary)
    val(read_type)

    output:
    path "${meta.id}_nextdenovo.config", emit: config

    script:
    // Calculate resources based on task.cpus and task.memory
    def total_memory = task.memory.toGiga()
    def memory_per_job = (total_memory / task.cpus).intValue()
    def parallel_jobs = task.cpus

    """
    # Extract max haploid genome length from GenomeScope2 summary
    genome_size=\$(awk '/Genome Haploid Length/ {print \$(NF-1)}' ${genomescope_summary} | sed 's/,//g')

    cat <<EOF > "${meta.id}_nextdenovo.config"
    [General]
    job_type = local
    job_prefix = nextDenovo
    task = all
    rewrite = yes
    deltmp = yes
    parallel_jobs = ${parallel_jobs}
    input_type = raw
    read_type = ${read_type}
    input_fofn = input.fofn
    workdir = 01_rundir

    [correct_option]
    read_cutoff = 1k
    genome_size = \${genome_size}
    sort_options = -m ${memory_per_job}g -t ${task.cpus}
    pa_correction = ${(task.cpus / 2).intValue()}
    correction_options = -p ${task.cpus}

    [assemble_option]
    minimap2_options_cns = -t ${task.cpus}
    nextgraph_options = -a 1 -t ${task.cpus}

    EOF

    if [ "${read_type}" == "ont" ]; then
        echo "minimap2_options_raw = -t ${task.cpus}" >> "${meta.id}_nextdenovo.config"
    fi
    """
}


process NEXTDENOVO {
    tag "$meta.id"
    label 'process_high'
    label 'process_assembly'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextdenovo:2.5.0--py39h20169af_0' :
        'quay.io/biocontainers/nextdenovo:2.5.0--py39h20169af_0' }"

    input:
    tuple val(meta), path(reads)
    path config_file

    output:
    tuple val(meta), path("${prefix}/03.ctg_graph/nd.asm.fasta"), emit: assembly
    tuple val(meta), path("${prefix}/03.ctg_graph/nd.asm.fasta.stat"), emit: stats
    tuple val(meta), path("${prefix}/03.ctg_graph/nd.asm.fasta.gfa"), emit: graph
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = "input.fofn"
    """
    # Create input.fofn file
    echo "${reads}" > ${input_file}

    # Run NextDenovo
    nextDenovo ${config_file}

    # Convert GFA to FASTA
    awk '/^S/{print ">"\$2;print \$3}' ${prefix}/03.ctg_graph/nd.asm.fasta.gfa > ${prefix}/03.ctg_graph/nd.asm.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextdenovo: \$(nextDenovo --version 2>&1 | sed 's/^.*NextDenovo //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/03.ctg_graph
    touch ${prefix}/03.ctg_graph/nd.asm.fasta
    touch ${prefix}/03.ctg_graph/nd.asm.fasta.stat
    touch ${prefix}/03.ctg_graph/nd.asm.fasta.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextdenovo: \$(nextDenovo --version 2>&1 | sed 's/^.*NextDenovo //')
    END_VERSIONS
    """
}