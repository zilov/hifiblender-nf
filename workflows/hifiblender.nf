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
include { MERYL_HISTOGRAM        } from "../modules/nf-core/meryl/histogram/main.nf"
include { MERYL_UNIONSUM         } from "../modules/nf-core/meryl/unionsum/main.nf"
include { GENOMESCOPE2           } from "../modules/nf-core/genomescope2/main.nf"
include { PARSE_INPUT            } from "../subworkflows/local/parse_input.nf"
// include { GENOMESCOPE2          } from '../modules/genomescope2/main'
// include { VERKKO                } from '../modules/verkko/main'
include { HIFIASM               } from '../modules/nf-core/hifiasm/main.nf'
include { FLYE                  } from '../modules/nf-core/flye/main.nf'
// include { NEXTDENOVO            } from '../modules/nextdenovo/main'
 include { QUAST                 } from '../modules/nf-core/quast/main.nf'
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

    // TODO: 
    // submodule for nanopore error correction
    // submodule for assembly polishing
    // submodule for hic-scaffolding
    // submodule for illumina assemmbly (prokaryotes)

    (ch_sample_map, ch_sample_reads, ch_fastqc_input) = PARSE_INPUT(ch_samplesheet)        

    // Parse tools list

    // Define valid tools and assemblers
    def valid_tools = ['qv', 'hifiasm', 'flye', 'nextdenovo', 'verkko', 'busco', 'compleasm', 'sniffles', 'slizer']
    def assemblers = ['hifiasm', 'flye', 'nextdenovo', 'verkko']

    // Parse the tools parameter
    def tools = params.tools instanceof List ? params.tools : params.tools.split(',').collect{ it.trim().toLowerCase() }

    // Validate tools
    def invalid_tools = tools - valid_tools
    if (!invalid_tools.isEmpty()) {
        error "Invalid tool(s) specified: ${invalid_tools.join(', ')}. Valid tools are: ${valid_tools.join(', ')}"
    }

    // Check if any assembler is specified, if not, add hifiasm
    if (tools.intersect(assemblers).isEmpty()) {
        tools.add('hifiasm')
        log.info "No assembler specified. Adding default assembler: hifiasm"
    }

    log.info "Tools to be run: ${tools.join(', ')}"

    
    //  
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastqc_input
    )

    //
    // SUBMODULE: QV
    // builds k-mer database with meryl, calculates coverage with genomescope and assess qv with merfin after assembly of genome 
    //

    if ('qv' in tools) {

        //
        // MODULE: Run Meryl db build 
        // takes as input illumina, ont or hifi reads, k value (default is 23)
        // outputs meryl kmer database and histogram
        // if several reads formats provided, combines databases with meryl union resulting in hybrid database
        // for parental phasing creates meryl database for each parent - used in verkko
        //

        MERYL_COUNT (
            ch_sample_reads, params.k
        )

        ch_meryl_db = MERYL_COUNT.out.meryl_db.view()

        // sum pair end meryl databases
        MERYL_UNIONSUM(
            ch_meryl_db, params.k
        )
        
        ch_meryl_sample_db = MERYL_UNIONSUM.out.meryl_db.view((sample, dbs) -> "$sample dbs $dbs")

        MERYL_HISTOGRAM (
            ch_meryl_sample_db, params.k
        )

        ch_meryl_sample_histogram = MERYL_HISTOGRAM.out.hist.view((sample, hist) -> "$sample histogram $hist")

        //
        // MODULE: Run genomescope2
        // takes as input kmers histogram from meryl and k value (default is 23)
        // outputs kmers stats - coverage value for QV estimation and others
        //

        GENOMESCOPE2 (
            ch_meryl_sample_histogram
        )
        
        
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

    }


    if ('hifiasm' in tools) {
        //
        // MODULE: Run hifiasm 
        // input for several modes - parental and maternal, ont-only, hifi-only, ont-hifi hybrid, hic
        // outputs assembly in several formats
        //

        hifiasm_input_ch = ch_sample_map
            .filter { sample -> sample.hifi || sample.ont }
            .map { sample ->
                def reads = sample.hifi ?: sample.ont
                [
                    [sample.sample, reads],
                    sample.parental ?: [],
                    sample.maternal ?: [],
                    sample.hic_1 ?: [],
                    sample.hic_2 ?: []
                ]
            }
        
        HIFIASM (
            hifiasm_input_ch.map { it[0] },  // meta, reads
            hifiasm_input_ch.map { it[1] },  // paternal_kmer_dump
            hifiasm_input_ch.map { it[2] },  // maternal_kmer_dump
            hifiasm_input_ch.map { it[3] },  // hic_read1
            hifiasm_input_ch.map { it[4] }   // hic_read2
        )
    }

    if ('verkko' in tools) {

        //
        // MODULE: Run verkko
        // input for several modes - ont-only, hifi-only, ont-hifi hybrid, hic
        // outputs assembly in several formats
        //

        // VERKKO (
        //     ch_samplesheet
        // )

    }
    
    if ('nextdenovo' in tools) {
        
        //
        // MODULE: Run nextdenovo
        //

        // NEXTDENOVO (
        //     ch_samplesheet
        // )

    }

    if ('flye' in tools) {

        //
        // MODULE: Run flye
        //

        flye_input_ch = ch_sample_map
        .filter { sample -> (sample.hifi || sample.ont) }
        .map { sample ->
            def reads = sample.hifi ?: sample.ont
            def mode = sample.hifi ? '--pacbio-hifi' : '--nano-raw'
            [
                sample.sample,  // meta
                reads,                  // reads
                mode                    // mode
            ]
        }

        FLYE (
            flye_input_ch.map { it[0..1] },  // meta, reads
            flye_input_ch.map { it[2] }      // mode
        )

    }



    //
    // MODULE: Run QUAST
    //

    quast_input_ch = ch_sample_map
        .map { sample ->
            def meta = sample.sample
            def consensus = sample.assembly ?: []  // Path to assembled genome
            def reference = sample.reference ?: []  // Path to reference genome, if available
            def gff = sample.gff ?: []  // Path to GFF file, if available
            
            [
                [ meta, consensus ],
                [ meta, reference ],
                [ meta, gff ]
            ]
        }

    QUAST (
        quast_input_ch.map { it[0] },  // meta, consensus
        quast_input_ch.map { it[1] },  // meta, reference (fasta)
        quast_input_ch.map { it[2] }   // meta, gff
    )

    if ('busco' in tools) {

        //
        // MODULE: Run BUSCO
        // on local lineage db if provided, or will download db (could define it with auto mode)
        //

        // BUSCO (
        //     ch_assembly
        // )

    }

    if ('compleasm' in tools) {

        //
        // MODULE: Run COMPLEASM
        // on local lineage db if provided, or will download db (could define it with auto mode)
        //

        // COMPLEASM (
        //     ch_assembly
        // )

    }


    if ('sniffles' in tools) {

        //
        // MODULE: Run SNIFFLES
        // assembly and reads are input if ONT and HIFI reads, ONT is used
        // outputs large structural variations in VCF
        //

        // SNIFFLES (
        //     ch_samplesheet
        //     ch_assembly
        // )

    }

    if ('slizer' in tools) {

        //
        // MODULE: Run SLIZER
        // assembly and reads are input if ONT and HIFI reads, HIFI is used
        // outputs low and high coverage regions based on z-score
        //

        // SLIZER (
        //     ch_samplesheet
        //     ch_assembly
        // )

    }


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