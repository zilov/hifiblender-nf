process PARSE_GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(model), path(summary)

    output:
    tuple val(meta), env(kmercov), env(estimated_genome_size)

    script:
    """
    # Extract kmercov from model file and format it as a whole number
    kmercov=\$(awk '/kmercov/ {printf "%.0f", \$2}' ${model} | awk '{printf "%d", \$1}')

    # Extract estimated genome size from summary file
    estimated_genome_size=\$(awk '/Genome Haploid Length/ {gsub(/,/, "", \$4); print \$4}' ${summary})

    # Output the results
    echo "kmercov=\$kmercov"
    echo "estimated_genome_size=\$estimated_genome_size"
    """

    stub:
    """
    kmercov=310.3
    estimated_genome_size=153603
    echo "kmercov: \$kmercov"
    echo "estimated_genome_size: \$estimated_genome_size"
    """
}