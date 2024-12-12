include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main.nf'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main.nf'

workflow MINIMAP2 {
    take:
        reference_ch
        reads_ch

    main:
        flat_reads = reads_ch.flatMap { meta, sampleMap ->
                def hifi_reads = [meta + [read_type: "hifi"], [sampleMap.hifi] - null]
                def ont_reads = [meta + [read_type: "ont"], [sampleMap.ont] - null]
                def illumina_reads = [meta + [read_type: "illumina-reads"] , [sampleMap.illumina_1, sampleMap.illumina_2] - null]
                def hic_reads = [meta + [read_type: "hic_reads"], [sampleMap.hic_1, sampleMap.hic_2] - null]
                return [ hifi_reads, ont_reads, illumina_reads, hic_reads ]
            }
            .filter { _meta, read_files ->
                !read_files.isEmpty()
            }

        product = reference_ch.cross(flat_reads) { v -> v[0]["sample_id"] }

        references = product.map { it[0] }
        reads = product.map { it[1] }

        (_paf, bam, _index, _versions_ch) = MINIMAP2_ALIGN(reads, references, true, 'bai', false, false)
    emit:
        bam_ch = bam
}
