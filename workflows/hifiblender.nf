/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hifiblender_pipeline'
include { MERYL_COUNT            } from '../modules/nf-core/meryl/count/main.nf'
include { PARSE_INPUT            } from "../subworkflows/local/parse_input.nf"
// include { GENOMESCOPE2          } from '../modules/genomescope2/main'
// include { VERKKO                } from '../modules/verkko/main'
// include { HIFIASM               } from '../modules/hifiasm/main'
// include { FLYE                  } from '../modules/flye/main'
// include { NEXTDENOVO            } from '../modules/nextdenovo/main'
// include { QUAST                 } from '../modules/quast/main'
// include { GFASTAT               } from '../modules/gfastat/main'
// include { BUSCO                 } from '../modules/busco/main'
// include { COMPLEASM             } from '../modules/compleasm/main'
// include { MERFIN                } from '../modules/merfin/main'
// include { SLIZER                } from '../modules/slizer/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HIFIBLENDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // THINK ABOUT: 
    // submodule for nanopore error correction
    // submodule for assembly polishing
    // submodule for hic-scaffolding
    // submodule for illumina assemmbly (prokaryotes)

    (ch_sample_map, ch_fastqc_input) = PARSE_INPUT(ch_samplesheet)        
    
    //  
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastqc_input
    )


    //
    // MODULE: Run Meryl db build 
    // takes as input illumina, ont or hifi reads, k value (default is 23)
    // outputs meryl kmer database and histogram
    // if several reads formats provided, combines databases with meryl union resulting in hybrid database
    // for parental phasing creates meryl database for each parent - used in verkko
    //
    // MERYL_COUNT (
    //     ch_samplesheet
    // )

    //
    // MODULE: Run genomescope2
    // takes as input kmers histogram from meryl and k value (default is 23)
    // outputs kmers stats - coverage value for QV estimation and others
    //

    // GENOMESCOPE2 (
    //     ch_samplesheet
    // )

    //
    // MODULE: Run hifiasm 
    // input for several modes - parental and maternal, ont-only, hifi-only, ont-hifi hybrid, hic
    // outputs assembly in several formats
    //

    // HIFIASM (
    //     ch_samplesheet
    // )

    //
    // MODULE: Run verkko
    // input for several modes - ont-only, hifi-only, ont-hifi hybrid, hic
    // outputs assembly in several formats
    //

    // VERKKO (
    //     ch_samplesheet
    // )

    //
    // MODULE: Run nextdenovo
    //

    // NEXTDENOVO (
    //     ch_samplesheet
    // )

    //
    // MODULE: Run flye
    //

    // FLYE (
    //     ch_samplesheet
    // )

    //
    // MODULE: Run QUAST
    //

    // QUAST (
    //     ch_assembly
    // )

    //
    // MODULE: Run BUSCO
    // on local lineage db if provided, or will download db (could define it with auto mode)
    //

    // BUSCO (
    //     ch_assembly
    // )

    //
    // MODULE: Run COMPLEASM
    // on local lineage db if provided, or will download db (could define it with auto mode)
    //

    // COMPLEASM (
    //     ch_assembly
    // )

    //
    // MODULE: Run SLIZER
    // assembly and reads are input if ONT and HIFI reads, ONT is used
    // outputs large structural variations in VCF
    //

    // SLIZER (
    //     ch_samplesheet
    //     ch_assembly
    // )

    //
    // MODULE: Run MERFIN
    // input assembly, coverage value and meryl database (hybrid if several reads types)
    // outputs QV, QV*, K completeness values
    //
    
    // MERFIN (
    //     ch_samplesheet
    //     ch_assembly
    //     ch_meryl
    //     ch_genomescope
    // )


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    //
    // MULTIQC: combine stats
    //

    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/