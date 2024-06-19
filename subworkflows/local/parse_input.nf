workflow PARSE_INPUT {
    take:
    ch_samplesheet

    main:
        // Process the sample sheet and create a map for each sample
        ch_sample_maps = ch_samplesheet
            .map { sample, assembler, hifi, ont, illumina_1, illumina_2, hic_1, hic_2, busco_lineage, busco_lineage_local, meryl_db_local, assembly, k ->
                def sampleMap = [
                    sample: sample,
                    assembler: assembler,
                    hifi: hifi ?: null,
                    ont: ont ?: null,
                    illumina_1: illumina_1 ?: null,
                    illumina_2: illumina_2 ?: null,
                    hic_1: hic_1 ?: null,
                    hic_2: hic_2 ?: null,
                    busco_lineage: busco_lineage,
                    busco_lineage_local: busco_lineage_local ?: null,
                    meryl_db_local: meryl_db_local ?: null,
                    assembly: assembly ?: null,
                    k: k
                ]
            }

        
        // Extract reads for FastQC
        ch_hifi_fastqc = ch_sample_maps
            .filter { it.hifi }
            .map { sampleMap ->
                [ sampleMap.sample, [ sampleMap.hifi ] ]
            }
        
        ch_ont_fastqc = ch_sample_maps
            .filter { it.ont }
            .map { sampleMap ->
                [ sampleMap.sample, [ sampleMap.ont ] ]
            }

        ch_illumina_fastqc = ch_sample_maps
            .filter { it.illumina_1 }
            .map { sampleMap ->
                def reads = sampleMap.illumina_2 != null ? [ sampleMap.illumina_1, sampleMap.illumina_2 ] : [ sampleMap.illumina_1 ]
                [ sampleMap.sample, reads ]
            }

        ch_hic_fastqc = ch_sample_maps
            .filter { it.hic_1 && it.hic_2 }
            .map { sampleMap ->
                [ sampleMap.sample, [ sampleMap.hic_1, sampleMap.hic_2 ] ]
            }

        // Combine all FastQC channels
        ch_fastqc_input = ch_hifi_fastqc.concat(ch_ont_fastqc, ch_illumina_fastqc, ch_hic_fastqc).view()


    emit:
        sampleMap = ch_sample_maps
        fastqc_reads = ch_fastqc_input
}
