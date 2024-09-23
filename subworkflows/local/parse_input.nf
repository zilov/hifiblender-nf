workflow PARSE_INPUT {
    take:
    ch_samplesheet

    main:
        ch_sample_maps = ch_samplesheet
            .map { sample, hifi, ont, illumina_1, illumina_2, hic_1, hic_2, busco_lineage, busco_lineage_local, meryl_db_local, assembly_fasta, assembly_gfa, k ->
                def sampleMap = [
                    sample: sample,
                    hifi: hifi ? file(hifi) : null,
                    ont: ont ? file(ont) : null,
                    illumina_1: illumina_1 ? file(illumina_1) : null,
                    illumina_2: illumina_2 ? file(illumina_2) : null,
                    hic_1: hic_1 ? file(hic_1) : null,
                    hic_2: hic_2 ? file(hic_2) : null,
                    busco_lineage: busco_lineage,
                    busco_lineage_local: busco_lineage_local ? file(busco_lineage_local) : null,
                    meryl_db_local: meryl_db_local ? file(meryl_db_local) : null,
                    assembly_fasta: assembly_fasta ? file(assembly_fasta) : null,
                    assembly_gfa: assembly_gfa ? file(assembly_gfa) : null,
                    k: k
                ]
                return sampleMap
            }

        ch_sample_reads = ch_sample_maps
            .map { sampleMap ->
                def reads = [sampleMap.hifi, sampleMap.ont, sampleMap.illumina_1, sampleMap.illumina_2] - null
                return [ sampleMap.sample, reads.flatten() ]
            }
            .filter { meta, reads -> !reads.isEmpty() }

        // Combine all FastQC channels
        ch_fastqc_input = ch_sample_maps
            .flatMap { sampleMap ->
                def reads = []
                if (sampleMap.illumina_1) {
                    reads << [sampleMap.sample, sampleMap.illumina_2 ? [file(sampleMap.illumina_1), file(sampleMap.illumina_2)] : [file(sampleMap.illumina_1)] ]
                }
                if (sampleMap.hic_1 && sampleMap.hic_2) {
                    reads << [sampleMap.sample, [file(sampleMap.hic_1), file(sampleMap.hic_2)] ]
                }
                if (sampleMap.ont && sampleMap.ont.toString().endsWith('.fastq.gz') || sampleMap.ont.toString().endsWith('.fastq')) {
                    reads << [sampleMap.sample, [file(sampleMap.ont)] ]
                }
                if (sampleMap.hifi && sampleMap.hifi.toString().endsWith('.fastq.gz') || sampleMap.hifi.toString().endsWith('.fastq.gz') )  {
                    reads << [sampleMap.sample, [file(sampleMap.hifi)] ]
                }
                return reads
            }

    emit:
        sampleMap = ch_sample_maps
        ch_sample_reads = ch_sample_reads
        ch_fastqc_input = ch_fastqc_input
}
