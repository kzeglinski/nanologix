// adapted from https://github.com/stevekm/nextflow-demos/blob/master/parse-samplesheet/main.nf
workflow parse_sample_sheet {
    take:
        fastq_dir
        sample_sheet

    main:
        Channel.fromPath( file(sample_sheet) )
        .splitCsv(header: true, sep: ',')
        .set {split_sample_sheet}

        split_sample_sheet
        .map{row ->
            def sample_ID = row['sample_num']
            def reads1 = "${fastq_dir}/" + row['r1_file_name']
            def reads2 = "${fastq_dir}/" + row['r2_file_name']
            return [ sample_ID, [reads1, reads2]]
        }
        .set { samples_R1_R2 } // set of all fastq R1 R2 per sample

        split_sample_sheet
        .map{row ->
            def report_ID = row['report_id']
            def sample_ID = row['sample_num']
            return [ sample_ID, report_ID ]
        }
        .set { report_sample }

    emit:
        samples_R1_R2
        report_sample

}