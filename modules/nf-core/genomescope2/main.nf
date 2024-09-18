process GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py311r42hdfd78af_6':
        'biocontainers/genomescope2:2.0--py311r42hdfd78af_6' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("${meta.id}/${prefix}_linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("${meta.id}/${prefix}_transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("${meta.id}/${prefix}_log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("${meta.id}/${prefix}_transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("${meta.id}/${prefix}_model.txt")                  , emit: model
    tuple val(meta), path("${meta.id}/${prefix}_summary.txt")                , emit: summary
    tuple val(meta), path("${meta.id}/${prefix}_lookup_table.txt")           , emit: lookup_table, optional: true
    tuple val(meta), path("${meta.id}/${prefix}_fitted_hist.png")            , emit: fitted_histogram_png, optional: true
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${meta.id}"
    """
    mkdir -p $outdir

    genomescope2 \\
        --input $histogram \\
        $args \\
        --output $outdir \\
        --name_prefix $prefix

    if [ -f "$outdir/fitted_hist.png" ]; then
        mv $outdir/fitted_hist.png ${outdir}/${prefix}_fitted_hist.png
    fi
    if [ -f "$outdir/lookup_table.txt" ]; then
        mv $outdir/lookup_table.txt ${outdir}/${prefix}_lookup_table.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "${meta.id}"
    """
    mkdir -p $outdir
    touch ${outdir}/${prefix}_linear_plot.png
    touch ${outdir}/${prefix}_transformed_linear_plot.png
    touch ${outdir}/${prefix}_log_plot.png
    touch ${outdir}/${prefix}_transformed_log_plot.png
    touch ${outdir}/${prefix}_model.txt
    touch ${outdir}/${prefix}_summary.txt
    touch ${outdir}/${prefix}_fitted_hist.png
    touch ${outdir}/${prefix}_lookup_table.txt

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """
}