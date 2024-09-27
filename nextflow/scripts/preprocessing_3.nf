#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Defining the parameters
params.samplefile='/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/test_data3/samples_3.csv'


process PREPROCESSING {
    publishDir "output", mode: 'copy', pattern: '*', overwrite: 'true'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    path "*_merged.fastq.gz", emit: merge_fastq 

    script:
    """
    bash /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/preprocessing_3.sh \
    -s ${sample_id} \
    -r ${fastq1} \
    -R ${fastq2}
    """
}

workflow {
    def samples = Channel
                    .fromPath("${params.samplefile}", type: 'file')
                    .splitCsv( header: true, sep: ',' )
                    .map {row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2))}

    samples.view()
    PREPROCESSING(samples)

    PREPROCESSING.out.merge_fastq.view()
}
