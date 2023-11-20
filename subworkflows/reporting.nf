process render_report {
    tag "creating report"
    label 'process_high'
    publishDir "${params.out_dir}/report/", mode: 'copy', pattern: "*.html"
   // stageInMode 'copy'
    container "library://kzeglinski/nanologix/nanologix-report:v0.3.0"

    input:
    path(processed_tsv)
    path(sample_sheet)
    path(multiqc_plots)
    path(percentage_passing_trim_merge)
    path(template_dir)
    path(extensions_dir)
    path(qmd_templates)
    val(analysis_name)

    output:
    path('*.html'), emit: report

    script:
    """
    #!/usr/bin/env bash

    export DENO_DIR="\$PWD"
    export XDG_CACHE_HOME="/tmp/quarto_cache_home"
    export XDG_DATA_HOME="/tmp/quarto_data_home"

    tar -xvf _extensions.tar
    tar -xvf _template.tar

    quarto render qc_report.qmd --log qc_report.log
    quarto render wnp_report.qmd --log wnp_report.log

    """
}

process prepare_report_templates {
    tag "preparing report templates"
    label 'process_low'
    stageInMode 'copy'
    container "library://kzeglinski/nanologix/nanologix-report:v0.3.0"

    input:
    path(sample_sheet)
    path(qmd_templates)
    val(analysis_name)

    output:
    path('*.qmd'), emit: report_templates

    script:
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(vroom)

    # read in the sample sheet
    metadata <- vroom(fs::dir_ls(glob = "*.csv"),
        col_select = c("sample_num", "library", "antigen", "round", "replicate"))

    # create replicate ID if required
    if(!all(is.na(metadata[["replicate"]]))){
        metadata <- metadata %>%
            filter(round != 0) %>% # temporarily remove R0 samples
            group_by(library, antigen, round) %>%
            mutate(replicate_id = cur_group_id()) %>%
            # and give our replicates an informative name, not just a number
            mutate(replicate_id_informative =
                paste0(library, "_", antigen, "_round_", round)) %>%
            ungroup() %>%
            bind_rows(filter(metadata, round == 0))
    }

    # create panning ID if required
    if(!all(metadata[["round"]] == 0)){
        metadata <- metadata %>%
            filter(round != 0) %>% # temporarily remove R0 samples
            group_by(library, antigen) %>%
            mutate(panning_id = as.character(cur_group_id())) %>%
            ungroup() %>%
            bind_rows(filter(metadata, round == 0))

        metadata <- metadata %>%
            group_by(library) %>%
            summarise(panning_id = paste0(unique(na.omit(panning_id)), collapse = "_")) %>%
            mutate(round = 0, replicate = NA) %>%
            right_join(metadata, by = c("library", "round", "replicate")) %>%
            mutate(panning_id = coalesce(panning_id.x, panning_id.y)) %>%
            select(-panning_id.x, -panning_id.y)
    }

    # qc_report
    readLines("template_qc_report.qmd") %>%
        stringr::str_replace(pattern = "param_analysis_name", replace = "${analysis_name}") %>%
        writeLines(con = "qc_report.qmd")

    # wnp_report
    readLines("template_wnp_report.qmd") %>%
        stringr::str_replace(pattern = "param_analysis_name", replace = "${analysis_name}") %>%
        writeLines(con = "wnp_report.qmd")

    """
}

