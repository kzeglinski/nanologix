// adapted from https://github.com/stevekm/nextflow-demos/blob/master/parse-samplesheet/main.nf
workflow parse_sample_sheet {
    take:
        fastq_dir
        sample_sheet

    main:
        Channel.fromPath( file(sample_sheet) )
        .splitCsv(header: true, sep: ',')
        .map{row ->
            def sample_ID = row['sample_num']
            def reads1 = "${fastq_dir}/" + row['r1_file_name']
            def reads2 = "${fastq_dir}/" + row['r2_file_name']
            return [ sample_ID, [reads1, reads2]]
        }
        .set { samples_R1_R2 } // set of all fastq R1 R2 per sample

    emit:
        samples_R1_R2

}