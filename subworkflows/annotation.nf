include { fastq_to_fasta } from '../modules/fastq_to_fasta'
include { igblast } from '../modules/igblast'

process merge_chunked_tsvs {
    tag "$sequence_id"
    label 'process_low'
    publishDir "${params.out_dir}/original_igblast_tsv", mode: 'copy', pattern: "*_merged.tsv"

    input:
    tuple val(sequence_id), path(tsvs)

    output:
    tuple val(sequence_id), path("*_merged.tsv")

    script:
    """
    cat *.tsv > ${sequence_id}_merged.tsv
    """
}

workflow annotation {
    take:
        chunked_merged_reads
        igblast_databases

    main:
        // convert fastq to fasta (needed to run igblast)
        chunked_merged_fastas = fastq_to_fasta(chunked_merged_reads)

        // set up environment variables
        igdata_dir="${igblast_databases}/igdata/"
        igblastdb_dir="${igblast_databases}/databases/"

        // annotate reads using igblast
        igblast_tsv = igblast(
            chunked_merged_fastas,
            igblast_databases,
            igdata_dir,
            igblastdb_dir).airr_table

        grouped_tsvs = igblast_tsv.groupTuple(by: 0) // group by first element (the sample ID)
        merged_tsvs = merge_chunked_tsvs(grouped_tsvs) // cat all of them

    emit:
        merged_tsvs
 }