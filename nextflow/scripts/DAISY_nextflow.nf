#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplefile = '/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/test_data3/samples_3.csv'
params.cellranger = '/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/test_data3/cellrangerpath.csv'

process PREPROCESSING {
    publishDir "output_preprocessing", mode: 'copy', pattern: '*', overwrite: true

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), path(cellranger)

    output:
    tuple val(sample_id), path(cellranger), path("*_merged.fastq.gz")

    script:
    """
    bash /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/Nextflow/preprocessing_3.sh \
    -s ${sample_id} \
    -r ${fastq1} \
    -R ${fastq2}
    """
}

process READLENGTH {
    publishDir "output_readlength", mode: 'copy', pattern: '*.png'

    input:
    tuple val(sample_id),  path(cellranger), path(merge_fastq)

    output:
    path "*.png"

    script:
    """
    seqkit watch --fields ReadLen ${merge_fastq} -O ${sample_id}_readlength.png
    """
}

process HOMOMATCH {
    publishDir "output_homodomain", mode: 'copy', pattern: '*_homodomain.fastq.gz'

    input:
    tuple val(sample_id), path(cellranger), path(merge_fastq)

    output:
    tuple val(sample_id), path(cellranger), path("*_homodomain.fastq.gz")

    script:
    """
    seqkit grep -s -m 2 -j 12 -p "TAGTTACGCCAAGCTTGAATTC" ${merge_fastq} -o ${sample_id}_homodomain.fastq.gz
    """
}

process FQ2FA {
    publishDir "output_homodomain", mode: 'copy', pattern: '*_homodomain.fasta'

    input:
    tuple val(sample_id), path(cellranger), path(homofq)

    output:
    tuple val(sample_id), path(cellranger), path("*_homodomain.fasta") 

    script:
    """
    seqkit fq2fa ${homofq} -o ${sample_id}_homodomain.fasta
    """
}

process BARCODEMATCH {
    publishDir "output_barcodechunk", mode: 'copy',  overwrite: true
    
    input:
    tuple val(sample_id), path(cellranger), path(homofasta)

    output:
    tuple val(sample_id), path(cellranger), path "*_10x_barcode_combined_nolinebreak.fasta"

    script:
    """
    zcat ${cellranger} | sed 's/-1//g' | awk '{print "^"\$0}'  > barcode.txt
    split -l 500 barcode.txt barcode_chunk_
    ls barcode_chunk_* | xargs -I {} seqkit grep -s -r -m 0 -j 4 -f {} ${homofasta} -o "{}_${sample_id}"
    cat barcode_chunk_*_${sample_id} > ${sample_id}_10x_barcode_combined.fasta
    seqtk seq -l 0 ${sample_id}_10x_barcode_combined.fasta > ${sample_id}_10x_barcode_combined_nolinebreak.fasta
    """
}

process ALIGNMENT {
    publishDir "output_chunkprocess", mode: 'copy',  overwrite: true
    
    input:
    tuple val(sample_id), path(cellranger), path(barcode_fasta)
    
    output:
    tuple val(sample_id), path(cellranger), path("*_alignment.csv"), path("${sample_id}")
    
    script:
    """
    python3 /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/sequence_extract.py ${barcode_fasta} ${sample_id}_10x_UMI_static_mut.csv
    awk -F, 'NF==4 && \$1!="" && \$2!="" && \$3!="" && \$4!=""' ${sample_id}_10x_UMI_static_mut.csv > ${sample_id}_10x_UMI_static_mut_clean.csv
    sort ${sample_id}_10x_UMI_static_mut_clean.csv | uniq -c | awk -F " " '{if (\$(NF-1) > 2) print \$0}' > ${sample_id}_10x_UMI_static_mut_clean_remove_PCR.csv
    cut -d "," -f1,3,4 ${sample_id}_10x_UMI_static_mut_clean_remove_PCR.csv | awk '{print \$NF}' | sort | uniq -c | awk -F " " '{print \$(NF-1)","\$NF}' > ${sample_id}_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv
    python3 /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/max_selection.py ${sample_id}_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv
    python3 /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/alignment_ATGCindel.py ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_alignment.csv
    python3 /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/alignment0_ATGC.py ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_0_alignment.csv
    python3 /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/alignment_0_1_ATGC.py ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_0_1_alignment.csv
    python3 /diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/DAISY/AJ/scripts/BCsort.py ${sample_id}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_alignment.csv ${sample_id}
    """
}

workflow {
    samples = Channel
        .fromPath(params.samplefile, type: 'file')
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_id, file(row.fastq1), file(row.fastq2), file(row.cellranger)) }

    processed = PREPROCESSING(samples)
    readlength = READLENGTH(processed)
    homomatch = HOMOMATCH(processed)
    homofasta = FQ2FA(homomatch)
    barcode_fasta = BARCODEMATCH(homofasta)
    alignment = ALIGNMENT(barcode_fasta)
}
