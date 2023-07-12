process run_r_script {
    tag "$sequence_id"
    label 'process_high'
    publishDir "${params.out_dir}/processed_tsv", mode: 'copy', pattern: "*.tsv"
    module 'R/4.2.3'

    input:
    tuple val(sequence_id), path(merged_tsv)

    output:
    tuple val(sequence_id), path('*.tsv'), emit: processed_tsv

    script:
    """
    #!/usr/bin/env Rscript
    #install.packages('stringr', repos='http://cran.us.r-project.org')
    #install.packages('dplyr', repos='http://cran.us.r-project.org')
    #install.packages('vroom', repos='http://cran.us.r-project.org')
    #install.packages('tidyr', repos='http://cran.us.r-project.org')
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(vroom)

    working_dir <- getwd()

    # CURRENT TSV FILE
    sample_id <- '$sequence_id'

    # reading in the data
    entire_data <- vroom('$merged_tsv',
                     col_select = c("stop_codon", "vj_in_frame", "v_frameshift",
                                    "productive", "v_call", "d_call", "j_call", "cdr3_aa",
                                    "sequence_alignment_aa"))

    # categories of productivity
    entire_data %>%
        select(c("productive", "stop_codon", "vj_in_frame", "v_frameshift")) %>%
        group_by(productive, stop_codon, vj_in_frame, v_frameshift) %>%
        drop_na() %>%
        summarise(count = n()) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(percentage = prop.table(count) * 100) %>%
        mutate(sample_id = sample_id) -> productivity_categories

    # germline V genes
    entire_data %>%
        filter(productive == TRUE) %>%
        select("v_call") %>%
        table() %>%
        as_tibble() %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        mutate(sample_id = sample_id) -> v_calls

    # germline D genes
    entire_data %>%
        select("d_call") %>%
        table() %>%
        as_tibble() %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        mutate(sample_id = sample_id) -> d_calls

    # germline J genes
    entire_data %>%
        filter(productive == TRUE) %>%
        select("j_call") %>%
        table() %>%
        as_tibble() %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        mutate(sample_id = sample_id) -> j_calls

    # CDR3 counts
    entire_data %>%
        filter(productive == TRUE) %>%
        select(cdr3_aa) %>%
        drop_na() %>%
        group_by(cdr3_aa) %>%
        summarise(count = n()) %>%
        mutate(proportion = count / sum(count)) %>%
        mutate(sample_id = sample_id) %>%
        arrange(desc(count)) -> cdr3_counts

    # whole nanobody aa sequence
    entire_data %>%
        filter(productive == TRUE) %>%
        select(c(sequence_alignment_aa, cdr3_aa)) %>%
        drop_na() %>%
        summarise(count = n(), cdr3_aa = cdr3_aa[1], .by = sequence_alignment_aa) %>%
        mutate(cpm = (count / sum(count)) * 1000000) %>%
        group_by(cdr3_aa) %>% # choose the most common aa seq for each cdr3
        summarise(
            sequence_alignment_aa = sequence_alignment_aa[which.max(cpm)],
            cdr3_count = n()) %>%
        filter(cdr3_count > 1) %>%
        mutate(cdr3_cpm = (cdr3_count / sum(cdr3_count)) * 1000000) %>%
        select(-cdr3_count) %>%
        mutate(sample_id = sample_id) %>%
        arrange(desc(cdr3_cpm)) -> whole_aa_seq

    # write out
    vroom_write(productivity_categories, "${sequence_id}_productivity_categories.tsv")
    vroom_write(v_calls, "${sequence_id}_v_calls.tsv")
    vroom_write(d_calls, "${sequence_id}_d_calls.tsv")
    vroom_write(j_calls, "${sequence_id}_j_calls.tsv")
    vroom_write(cdr3_counts, "${sequence_id}_cdr3_counts.tsv")
    vroom_write(whole_aa_seq, "${sequence_id}_whole_aa_seq.tsv")
    """
}

workflow r_processing {
    take:
        igblast_tsvs

    main:
        // convert fastq to fasta (needed to run igblast)
        processed_tsv = run_r_script(igblast_tsvs)

    emit:
        processed_tsv
 }