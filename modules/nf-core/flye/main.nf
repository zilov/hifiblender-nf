process FLYE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("${meta.id}/*.fasta.gz"), emit: fasta, topic: assembly_fasta
    tuple val(meta), path("${meta.id}/*.gfa.gz")  , emit: gfa, topic: assembly_gfa
    tuple val(meta), path("${meta.id}/*.gv.gz")   , emit: gv
    tuple val(meta), path("${meta.id}/*.txt")     , emit: txt
    tuple val(meta), path("${meta.id}/*.log")     , emit: log
    tuple val(meta), path("${meta.id}/*.json")    , emit: json
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def outdir = "${meta.id}"
    def prefix = task.ext.prefix ?: "${meta.id}"

    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }
    """
    flye \\
        $mode \\
        $reads \\
        --out-dir $outdir \\
        --threads $task.cpus \\
        $args

    gzip -c $outdir/assembly.fasta > $outdir/${prefix}.assembly.fasta.gz
    gzip -c $outdir/assembly_graph.gfa > $outdir/${prefix}.assembly_graph.gfa.gz
    gzip -c $outdir/assembly_graph.gv > $outdir/${prefix}.assembly_graph.gv.gz
    mv $outdir/assembly_info.txt $outdir/${prefix}.assembly_info.txt
    mv $outdir/flye.log $outdir/${prefix}.flye.log
    mv $outdir/params.json $outdir/${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${meta.id}"
    """
    mkdir -p $outdir

    echo stub | gzip -c > $outdir/${prefix}.assembly.fasta.gz
    echo stub | gzip -c > $outdir/${prefix}.assembly_graph.gfa.gz
    echo stub | gzip -c > $outdir/${prefix}.assembly_graph.gv.gz
    echo contig_1 > $outdir/${prefix}.assembly_info.txt
    echo stub > $outdir/${prefix}.flye.log
    echo stub > $outdir/${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}