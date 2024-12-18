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
include { PARSE_GENOMESCOPE2           } from "../modules/local/parse_genomescope2/main.nf"
include { GFASTATS               } from "../modules/nf-core/gfastats/main.nf"
include { PARSE_INPUT            } from "../subworkflows/local/parse_input.nf"
include { MINIMAP2 } from "../subworkflows/local/minimap2.nf"
include { VERKKO                 } from '../modules/local/verkko/main.nf'
include { CREATE_NEXTDENOVO_CONFIG } from '../modules/local/nextdenovo/main.nf'
include { NEXTDENOVO } from '../modules/local/nextdenovo/main.nf'
include { HIFIASM                } from '../modules/nf-core/hifiasm/main.nf'
include { FLYE                   } from '../modules/nf-core/flye/main.nf'
include { QUAST                 } from '../modules/nf-core/quast/main.nf'
include { BUSCO_BUSCO                 } from '../modules/nf-core/busco/busco/main'
// include { COMPLEASM             } from '../modules/compleasm/main'
include { MERFIN_HIST                } from '../modules/nf-core/merfin/hist/main.nf'
// include { SLIZER                } from '../modules/slizer/main'


// features

nextflow.preview.topic = true

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

    ch_meryl_db = MERYL_COUNT.out.meryl_db

    // sum pair end meryl databases
    MERYL_UNIONSUM(
        ch_meryl_db, params.k
    )
    
    ch_meryl_sample_db = MERYL_UNIONSUM.out.meryl_db


    MERYL_HISTOGRAM (
        ch_meryl_sample_db, params.k
    )

    ch_meryl_sample_histogram = MERYL_HISTOGRAM.out.hist

    //
    // MODULE: Run genomescope2
    // takes as input kmers histogram from meryl and k value (default is 23)
    // outputs kmers stats - coverage value for QV estimation and others
    //

    ch_genomescope_out = GENOMESCOPE2 (
        ch_meryl_sample_histogram
    )
    
    // prepare input for PARSE_GENOMESCOPE2
    ch_genomescope_parse_input = ch_genomescope_out.model
                                    .combine(ch_genomescope_out.summary, by: 0)

    // extracts kmercov peak, estimated length (required for nextdenovo run) and lookup table
    ch_parsed_genomescope = PARSE_GENOMESCOPE2(ch_genomescope_parse_input)

    // Process samples that need assembly
    
    def assembler_channels_fa = []
    def assembler_channels_gfa = []

    ch_to_assemble = ch_sample_map.filter { sample -> !sample.assembly_fasta and !sample.assembly_gfa }
    
    input_assembly_ch_fa = ch_sample_map.filter {sample -> sample.assembly_fasta}
        .map {sample -> [[id: "${sample.sample.id}_user_provided", sample_id: sample.sample.id, assembler: "user_provided"], sample.assembly_fasta]}
    
    input_assembly_ch_gfa = ch_sample_map.filter {sample -> sample.assembly_gfa}
        .map {sample -> [[id: "${sample.sample.id}_user_provided", sample_id: sample.sample.id, assembler: "user_provided"], sample.assembly_gfa]}

    assembler_channels_fa << input_assembly_ch_fa
    assembler_channels_gfa << input_assembly_ch_gfa


    if ('hifiasm' in tools) {
        //
        // MODULE: Run hifiasm 
        // input for several modes - parental and maternal, ont-only, hifi-only, ont-hifi hybrid, hic
        // outputs assembly in several formats
        //

        hifiasm_input_ch = ch_to_assemble
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
        
        hifiasm_output_ch = HIFIASM (
            hifiasm_input_ch.map { it[0] },  // meta, reads
            hifiasm_input_ch.map { it[1] },  // paternal_kmer_dump
            hifiasm_input_ch.map { it[2] },  // maternal_kmer_dump
            hifiasm_input_ch.map { it[3] },  // hic_read1
            hifiasm_input_ch.map { it[4] }   // hic_read2
        )

        hifiasm_assembly_fasta = hifiasm_output_ch.primary_contigs_fa
            .map { meta, fasta -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_hifiasm"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'hifiasm'
                [ newMeta, fasta ]
            }
        
        assembler_channels_fa << hifiasm_assembly_fasta

        hifiasm_assembly_gfa = hifiasm_output_ch.primary_contigs
            .map { meta, gfa -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_hifiasm"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'hifiasm'
                [ newMeta, gfa ]
            }

        assembler_channels_gfa << hifiasm_assembly_gfa

    }

    if ('verkko' in tools) {

        //
        // MODULE: Run verkko
        // input for several modes - ont-only, hifi-only, ont-hifi hybrid, hic
        // outputs assembly in several formats
        //

        verkko_input_ch = ch_to_assemble
            .filter { sample -> (sample.hifi) }
            .map { sample ->
                [       
                        sample.sample,
                        sample.hifi ? file(sample.hifi) : [], // hifi
                        sample.ont ? file(sample.ont) : [], //ont
                        sample.hic_1 ? file(sample.hic_1) : [], // hic_1
                        sample.hic_2 ? file(sample.hic_2) : [], // hic_2
                        [], // porec
                        sample.parental_kmers ? file(sample.parental_kmers) : [], // paternal hapmers
                        sample.maternal_kmers ? file(sample.maternal_kmers) : [], // maternal hapmers
                        [],
                        []
                ]
            }

        verkko_output_ch = VERKKO(
            verkko_input_ch.map { it[0] },  // meta
            verkko_input_ch.map { it[1] },  // hifi
            verkko_input_ch.map { it[2] },  // ont
            verkko_input_ch.map { it[3] },  // hic_1
            verkko_input_ch.map { it[4] },  // hic_2
            verkko_input_ch.map { it[5] },  // porec
            verkko_input_ch.map { it[6] },  // paternal hapmers
            verkko_input_ch.map { it[7] },  // maternal hapmers
            verkko_input_ch.map { it[8] },  // telomere_motif
            verkko_input_ch.map { it[9] },  // reference
        )

        verkko_assembly_fasta = verkko_output_ch.assembly
            .map { meta, fasta -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_verkko"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'verkko'
                [ newMeta, fasta ]
            }

        assembler_channels_fa << verkko_assembly_fasta

        verkko_assembly_gfa = verkko_output_ch.graph
            .map { meta, gfa -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_verkko"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'verkko'
                [ newMeta, gfa ]
            }
        
        assembler_channels_gfa << verkko_assembly_gfa

    }
    
    if ('nextdenovo' in tools) {
        
        //
        // MODULE: Run NextDenovo
        //

        nextdenovo_input_ch = ch_to_assemble
            .filter { sample -> sample.hifi || sample.ont }
            .map { sample ->
                def reads = sample.hifi ?: sample.ont
                def read_type = sample.hifi ? 'hifi' : 'ont'
                [
                    sample.sample,  // meta
                    reads,          // reads
                    read_type       // read_type
                ]
            }

        // Combine the nextdenovo input channel with the GenomeScope results
        nextdenovo_combined_ch = nextdenovo_input_ch
            .join(ch_genomescope_out.summary, by: 0)  // Join based on the sample ID (assuming it's in the first position of both channels)

        // Create NextDenovo config
        nextdenovo_config_ch = CREATE_NEXTDENOVO_CONFIG(
            nextdenovo_combined_ch.map { it[0] },  // meta
            nextdenovo_combined_ch.map { it[3] },  // GenomeScope summary
            nextdenovo_combined_ch.map { it[2] }   // read_type
        )

        // Run NextDenovo
        nextdenovo_output_ch = NEXTDENOVO(
            nextdenovo_input_ch.map { it[0..1] },  // [meta, reads]
            nextdenovo_config_ch.config
        )

        nextdenovo_assembly_fasta = nextdenovo_output_ch.assembly
            .map { meta, fasta -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_nextdenovo"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'nextdenovo'
                [ newMeta, fasta ]
            }

        assembler_channels_fa << nextdenovo_assembly_fasta

        nextdenovo_assembly_gfa = nextdenovo_output_ch.graph
            .map { meta, gfa -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_nextdenovo"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'nextdenovo'
                [ newMeta, gfa ]
            }
        
        assembler_channels_gfa << nextdenovo_assembly_gfa
    }

    if ('flye' in tools) {

        //
        // MODULE: Run flye
        //

        flye_input_ch = ch_to_assemble
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

        flye_output_ch = FLYE (
            flye_input_ch.map { it[0..1] },  // meta, reads
            flye_input_ch.map { it[2] }      // mode
        )

        flye_assembly_fasta = flye_output_ch.fasta
            .map { meta, fasta -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_flye"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'flye'
                [ newMeta, fasta ]
            }

        assembler_channels_fa << flye_assembly_fasta

        flye_assembly_gfa = flye_output_ch.gfa
            .map { meta, gfa -> 
                def newMeta = meta.clone()
                newMeta.id = "${meta.id}_flye"
                newMeta.sample_id = meta.id
                newMeta.assembler = 'flye'
                [ newMeta, gfa ]
            }
        
        assembler_channels_gfa << flye_assembly_gfa
   }


    // Combine fasta assembly channels
    combined_fa_channel = Channel.empty()
    assembler_channels_fa.each { channel ->
        combined_fa_channel = combined_fa_channel.mix(channel)
    }

    // Combine gfa assembly channels
    combined_gfa_channel = Channel.empty()
    assembler_channels_gfa.each { channel ->
        combined_gfa_channel = combined_gfa_channel.mix(channel)
    }
    

    //
    // MODULE: Run QUAST
    //

    quast_out_ch = QUAST (
        combined_fa_channel,  // meta, consensus
        [[:], []],  // meta2, reference (fasta)
        [[:], []]   // meta3, gff
    )

    labeled_reads = ch_to_assemble.map { sampleMap -> [sampleMap.sample + [sample_id: sampleMap.sample["id"]], sampleMap] }
    result = MINIMAP2(combined_fa_channel, labeled_reads)


    


    // MINIMAP2(assembler_channels_fa, flat_reads)

    //
    // MODULE: Run GFASTATS
    //

    gfastats_out_ch = GFASTATS(combined_gfa_channel, 'gfa', '', '', [], [], [], [])

    if ('busco' in tools) {

        
        // MODULE: Run BUSCO
        // on local lineage db if provided, or will download db (could define it with auto mode)

        // Extract BUSCO lineage information from ch_sample_map
        busco_lineage_info = ch_sample_map
            .map { sample ->
                def lineage = sample.busco_lineage ?: 'auto'
                def lineage_local = sample.busco_lineage_local ?: []
                [sample.sample.id, lineage, lineage_local]
            }

        // Combine combined_fa_channel with BUSCO lineage information
        busco_input_channel = combined_fa_channel
            .map { meta, fasta -> [meta.sample_id, meta, fasta] }
            .combine(busco_lineage_info, by: 0)
            .map { sample_id, meta, fasta, lineage, lineage_local ->
                [meta, fasta, lineage, lineage_local]
            }
        

        BUSCO_BUSCO (
            busco_input_channel.map { it[0..1] },  // meta, fasta
            'genome', // mode
            busco_input_channel.map { it[2] }, // lineage
            busco_input_channel.map { it[3] }, // lineage_dir
            [], // busco config
        )

    }

    if ('qv' in tools) {

        //
        // MODULE: Run MERFIN
        // input assembly, coverage value and meryl database (hybrid if several reads types)
        // outputs QV, QV*, K completeness values
        //

        meryl_db_for_merfin = ch_meryl_sample_db.map { meta, db -> [meta.id, db]}
        parse_genomescope2_for_merfin = ch_parsed_genomescope.map{ meta, cov, length -> [meta.id, cov]}
        lookup_table = ch_genomescope_out.lookup_table ? ch_sample_map.map{sample -> [sample.sample.id, []]} : ch_genomescope_out.lookup_table

        // First, let's prepare the main input channel
        merfin_input_channel = combined_fa_channel.map { meta, fasta -> 
                [meta.sample_id, meta, fasta]
            }
            .combine(meryl_db_for_merfin, by:0)
            .combine(parse_genomescope2_for_merfin, by:0)
            .combine(lookup_table, by:0)
            .map { sample_id, meta, fasta, meryl_db, kmercov, lookup ->
                 [meta, fasta, meryl_db, kmercov, lookup]
             }.view()

        merfin_output = MERFIN_HIST (
            merfin_input_channel
                .map{meta, fasta, meryl_db, kmercov, lookup -> [meta, fasta]},
            merfin_input_channel
                .map{meta, fasta, meryl_db, kmercov, lookup -> [meta, meryl_db]},
            merfin_input_channel
                .map{meta, fasta, meryl_db, kmercov, lookup -> lookup},
            [], //seqmers
            merfin_input_channel
                .map{meta, fasta, meryl_db, kmercov, lookup -> kmercov},
        )

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