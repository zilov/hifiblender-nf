/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERYL                 } from '../modules/nf-core/meryl/main'
include { GENOMESCOPE2          } from '../modules/nf-core/genomescope2/main'
include { VERKKO                } from '../modules/nf-core/verkko/main'
include { HIFIASM               } from '../modules/nf-core/hifiasm/main'
include { FLYE                  } from '../modules/nf-core/flye/main'
include { NEXTDENOVO            } from '../modules/nf-core/nextdenovo/main'
include { QUAST                 } from '../modules/nf-core/quast/main'
include { GFASTAT               } from '../modules/nf-core/gfastat/main'
include { BUSCO                 } from '../modules/nf-core/busco/main'
include { COMPLEASM             } from '../modules/nf-core/compleasm/main'
include { MERFIN                } from '../modules/nf-core/merfin/main'
include { SLIZER                } from '../modules/nf-core/slizer/main'

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

    //
    // MODULE: Run Meryl db build 
    // takes as input illumina, ont or hifi reads, k value (default is 23)
    // outputs meryl kmer database and histogram
    // if several reads formats provided, combines databases with meryl union resulting in hybrid database
    // for parental phasing creates meryl database for each parent - used in verkko
    //
    MERYL (
        ch_samplesheet
    )

    //
    // MODULE: Run genomescope2
    // takes as input kmers histogram from meryl and k value (default is 23)
    // outputs kmers stats - coverage value for QV estimation and others
    //

    GENOMESCOPE2 (
        ch_samplesheet
    )

    //
    // MODULE: Run hifiasm 
    // input for several modes - parental and maternal, ont-only, hifi-only, ont-hifi hybrid, hic
    // outputs assembly in several formats
    //
    HIFIASM (
        ch_samplesheet
    )

    //
    // MODULE: Run verkko
    // input for several modes - ont-only, hifi-only, ont-hifi hybrid, hic
    // outputs assembly in several formats
    //
    VERKKO (
        ch_samplesheet
    )

    //
    // MODULE: Run nextdenovo
    //
    NEXTDENOVO (
        ch_samplesheet
    )

    //
    // MODULE: Run flye
    //
    FLYE (
        ch_samplesheet
    )

    //
    // MODULE: Run QUAST
    //
    QUAST (
        ch_assembly
    )

    //
    // MODULE: Run BUSCO
    // on local lineage db if provided, or will download db (could define it with auto mode)
    //
    BUSCO (
        ch_assembly
    )

    //
    // MODULE: Run COMPLEASM
    // on local lineage db if provided, or will download db (could define it with auto mode)
    //
    COMPLEASM (
        ch_assembly
    )

    //
    // MODULE: Run SLIZER
    // assembly and reads are input if ONT and HIFI reads, ONT is used
    // outputs large structural variations in VCF
    //
    SLIZER (
        ch_samplesheet
        ch_assembly
    )

    //
    // MODULE: Run MERFIN
    // input assembly, coverage value and meryl database (hybrid if several reads types)
    // outputs QV, QV*, K completeness values
    //
    MERFIN (
        ch_samplesheet
        ch_assembly
        ch_meryl
        ch_genomescope
    )


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    //
    // MULTIQC: combine stats
    //

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
