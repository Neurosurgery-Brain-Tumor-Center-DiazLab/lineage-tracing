#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Defining the parameters
params.samplefile='/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/test_data3/samples_3.csv'
params.intermediate_dir='/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate_5/'
params.output_dir="/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results_4/"


process PREPROCESSING {
    // publishDir "${params.intermediate_dir}", mode: 'copy', pattern: '*', overwrite: 'true'
    publishDir params.output_dir, mode: 'copy', pattern: '*_merged.fastq.gz'
    
    // Terminate if any error comes
    errorStrategy 'terminate'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    path intermediate
    path output

    output:
    // path "SB28_DAISY_DAY0_sub1_A1_S120_L007_I1_001.fastq.gz_merged.fastq.gz"
    path "*"

    script:
    """
    bash /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/preprocessing_2.sh \
    -s ${sample_id} \
    -r ${fastq1} \
    -R ${fastq2} \
    -m ${intermediate} \
    -o ${output}
    """
}

workflow {
    def samples = Channel
                    .fromPath("${params.samplefile}", type: 'file')
                    .splitCsv( header: true, sep: ',' )
                    .map {row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2))}
    def interm = Channel.fromPath("${params.intermediate_dir}", type: 'dir')
    def out = Channel.fromPath("${params.output_dir}", type: 'dir')

    samples.view()
    interm.view()
    out.view()

    PREPROCESSING(samples, interm, out)
}
