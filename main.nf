#!/usr/bin/env nextflow

/*
 * A nextflow pipeline for the pre-processing of 
 * 2x300 illumina sequencing data from nanobodies
 */

version = "v0.1.0"

if(params.help == true){
log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩ 
ξ ･×･ ξ 
ξ ~  ξ   		★ Nanobody pre-processing ★
ξ　  ξ              			   v0.1.0 
ξ　  “~~~~~~_		
ξ　          ξ 
 ξ ξ ξ~~~ξ ξ       		
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙

Usage: nextflow run ./nanobody_preprocessing/main.nf --fastq_dir [input path] --sample_sheet [sample sheet]
--help                : prints this help message
Required arguments:
--out_dir             : where the output files will be written to (default: "$projectDir/results)
--fastq_dir           : where the input fastq files are located
--sample_sheet        : location of the .csv sample sheet (format: sample_num,library,antigen,round,replicate)
Optional (only needed for advanced users)
--igblast_databases   : location of the igblast databases (default: "$projectDir/igblast_refs/")
--adapter_r1          : pattern for trimgalore to trim off R1 (default: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
--adapter_r2          : pattern for trimgalore to trim off R2 (default: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
--maximum_overlap     : maximum overlap (in bp) of the two reads (default: 200)
--chunk_size          : size of chunks to use for IgBLAST processing (default: 100000)
"""
System.exit(0)
}

// TODO: parameter validation
// validate that these files/directories exist
fastq_dir = params.fastq_dir
sample_sheet = params.sample_sheet
igblast_databases = params.igblast_databases

// validate that these are strings 
adapter_r1 = params.adapter_r1
adapter_r2 = params.adapter_r2

// validate these are numbers with appropriate limits
maximum_overlap = params.maximum_overlap
chunk_size = params.chunk_size

log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩ 
ξ ･×･ ξ 
ξ ~  ξ   		★ Nanobody pre-processing ★
ξ　  ξ              			   v0.1.0 
ξ　  “~~~~~~_		
ξ　          ξ 
 ξ ξ ξ~~~ξ ξ       		
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
★ read directory           : ${params.fastq_dir}
★ sample sheet             : ${params.sample_sheet}
★ output directory         : ${params.out_dir}
"""

/*
 * Bring in modules
 */
include { parse_sample_sheet } from './subworkflows/file_import'
include { trim_merge } from './subworkflows/trim_merge' 
include { quality_control } from './subworkflows/quality_control'
include { annotation } from './subworkflows/annotation'
include { r_processing } from './subworkflows/r_processing'

/*
 * Run the workflow
 */

workflow{
    // read the sample sheet to associate sample names and fastq files
    // output is a tuple with [sample_num, [R1, R2]]
    parse_sample_sheet(fastq_dir, sample_sheet)
    sample_read_pairs = parse_sample_sheet.out.samples_R1_R2

    // trim and merge the data
    trim_merge(sample_read_pairs, adapter_r1, adapter_r2, maximum_overlap)
    trimmed_and_merged_reads = trim_merge.out.trimmed_and_merged_reads
    fastqc_reports = trim_merge.out.fastqc_reports.collect()
    
    // perform quality control (multiQC on the fastQC output),
    // as well as simple checking of how many reads made it 
    // through the trim/merge process
    quality_control(sample_read_pairs, trimmed_and_merged_reads, fastqc_reports)
    
    // annotation using IgBLAST
    // first let's split each file into n chunks (default: 100,000)
    trimmed_and_merged_reads
        .splitFastq(by: chunk_size, file: true)
        .set{chunked_merged_reads}

    igblast_tsvs = annotation(chunked_merged_reads, igblast_databases)
    //igblast_tsvs.view()
    // R processing of the IgBLAST output
    r_processing(igblast_tsvs)
}
/* 
log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩          
ξ ･×･ ξ         ~~ THANK-YOU FOR RUNNING ~~
ξ ~  ξ   		★ Nanobody pre-processing ★
ξ　  ξ              			   v0.1.0 
ξ　  “~~~~~~_		
ξ　          ξ 
 ξ ξ ξ~~~ξ ξ       		
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
You might like to look at the following:
★ output directory                             : ${params.out_dir}
★ multi QC report                              : ${params.out_dir}/multiqc_report.html
★ how many reads passed the trimming & merging : ${params.out_dir}/percentage_passing_trim_merge.tsv
""" */

