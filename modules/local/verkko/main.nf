process VERKKO {
    tag "$meta.id"
    label 'process_high'
    label 'process_assembly'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verkko%3A2.2--h45dadce_0' :
        'quay.io/biocontainers/verkko:2.2--h45dadce_0' }"

        

    input:
    val(meta)
    path(hifi_reads)
    path(ont_reads)
    path(hic_reads_1)
    path(hic_reads_2)
    path(porec_reads)
    path(pat_hapmer)
    path(mat_hapmer)
    val telomere_motif
    path reference

    output:
    tuple val(meta), path("${prefix}/assembly.fasta"), emit: assembly
    tuple val(meta), path("${prefix}/assembly.homopolymer-compressed.gfa"), emit: graph
    tuple val(meta), path("${prefix}/assembly.haplotype*.fasta"), optional: true, emit: haplotypes
    tuple val(meta), path("${prefix}/assembly.colors.csv"), optional: true, emit: colors
    tuple val(meta), path("${prefix}/assembly.paths.tsv"), optional: true, emit: paths
    tuple val(meta), path("${prefix}/assembly.*.fasta"), optional: true, emit: contaminants
    tuple val(meta), path("${prefix}/assembly.*.exemplar.fasta"), optional: true, emit: exemplars
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def snakemake_scratchdir = "$workDir/snakmake_scratch"
    def args = task.ext.args ?: params.verkko_args
    prefix = task.ext.prefix ?: "${meta.id}"
    def hifi_input = hifi_reads ? "--hifi $hifi_reads" : ""
    def ont_input = ont_reads ? "--nano $ont_reads" : ""
    def hic_input = hic_reads_1 && hic_reads_2 ? "--hic1 $hic_reads_1 --hic2 $hic_reads_2" : ""
    def porec_input = porec_reads ? "--porec $porec_reads" : ""
    def hapmer_input = pat_hapmer && mat_hapmer ? "--hap-kmers $pat_hapmer $mat_hapmer trio" : ""
    def telomere_arg = telomere_motif ? "--telomere-motif $telomere_motif" : ""
    def ref_arg = reference ? "--ref $reference" : ""

    """
    # set snakemake scratch dir in process dir see https://github.com/snakemake/snakemake/issues/1593
    mkdir \${PWD}/snakemake_scratch
    export XDG_CACHE_HOME=\${PWD}/snakemake_scratch

    verkko \\
        -d $prefix \\
        $hifi_input \\
        $ont_input \\
        $hic_input \\
        $porec_input \\
        $hapmer_input \\
        $telomere_arg \\
        $ref_arg \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verkko: \$(verkko --version 2>&1 | sed 's/^.*verkko //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/assembly.fasta
    touch ${prefix}/assembly.homopolymer-compressed.gfa
    touch ${prefix}/assembly.haplotype1.fasta
    touch ${prefix}/assembly.haplotype2.fasta
    touch ${prefix}/assembly.colors.csv
    touch ${prefix}/assembly.paths.tsv
    touch ${prefix}/assembly.contaminant.fasta
    touch ${prefix}/assembly.contaminant.exemplar.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verkko: \$(verkko --version 2>&1 | sed 's/^.*verkko //')
    END_VERSIONS
    """
}
