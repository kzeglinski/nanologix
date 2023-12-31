// module for running IgBLAST
// general layout is based on the nf-core modules
process igblast {
    tag "$sequence_id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::igblast=1.19.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igblast%3A1.19.0--pl5321h3928612_0' :
        'quay.io/biocontainers/igblast:1.19.0--pl5321h3928612_0' }"

    input:
    tuple val(sequence_id), path(reads)
    val igblast_databases
    env IGDATA
    env IGBLASTDB

    output:
    tuple val(sequence_id), path('*_igblast.tsv'), emit: airr_table

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read_base_name=\$(basename "$reads" .fastq)
    # run igblast
    # outfmt 19 = AIRR format (tsv, easy to use in downstream steps)
    # num_alignments_* = only report the best hit
    igblastn -germline_db_V ${igblast_databases}/databases/imgt_alpaca_ighv \
        -germline_db_J ${igblast_databases}/databases/imgt_alpaca_ighj \
        -germline_db_D ${igblast_databases}/databases/imgt_alpaca_ighd \
        -organism alpaca \
        -query $reads \
        -num_threads $task.cpus \
        -auxiliary_data ${igblast_databases}/igdata/optional_file/alpaca_gl.aux \
        -show_translation \
        -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 \
        -outfmt 19 > \${read_base_name}_igblast.tsv
    """
}