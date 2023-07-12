// convert fastq to fasta (needed for IgBLAST)
process fastq_to_fasta {
    tag "$sequence_id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(sequence_id), path(reads)

    output:
    tuple val(sequence_id), path("*.fasta")

    script:
    """
    read_base_name=\$(basename "$reads" .fastq)
    seqkit fq2fa --threads $task.cpus $reads -o "\${read_base_name}.fasta"
    """
}