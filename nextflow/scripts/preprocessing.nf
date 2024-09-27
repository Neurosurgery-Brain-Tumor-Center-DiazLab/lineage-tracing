#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Defining the parameters
// Please make these directories before running the pipeline
// params.input_dir = '/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/test_data4/'
// params.intermediate_dir='/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate_4/' 
// params.output_dir="/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results_4/"
params.samplefile = '/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/test_data3/samples_4.csv'

// Define the channels here it might work
process preprocessing {
    // publishDir "${params.intermediate_dir}", mode: 'copy'
    // publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(input), path(intermediate), path(output)
    // path input
    // path intermediate
    // path output

    script:
    """
    bash /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/preprocessing.sh -i ${input} -m ${intermediate} -o ${output}
    """
}

workflow {
    def samples = Channel
                    .fromPath("${params.samplefile}", type: 'file')
                    .splitCsv( header: true, sep: ',' )
                    .map {row -> tuple(row.input, file(row.intermediate), file(row.output))}
    // def interm = Channel.fromPath("${params.intermediate_dir}", type: 'dir')
    // def out = Channel.fromPath("${params.output_dir}", type: 'dir')

    samples.view()
    // interm.view()
    // out.view()

    preprocessing(samples)
}
