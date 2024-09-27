#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -m <intermediate_directory> -c <cellranger_directory>"
    exit 1
}

# Parse command-line arguments
while getopts ":m:c:" opt; do
  case ${opt} in
    m )
      intermediate_dir=$OPTARG
      ;;
    c )
      cellranger_dir=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      usage
      ;;
  esac
done
shift $((OPTIND -1))

# To check the required arguments are provided
if [ -z "${intermediate_dir}" ] || [ -z "${cellranger_dir}" ]; then
    echo "Missing required arguments"
    usage
fi

# Extracting the 10x barcode from the cellranger output
# Remove '-1' from barcodes and add '^' to grep from the start
zcat "${cellranger_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" | sed 's/-1//g' | awk '{print "^"$0}' > "${intermediate_dir}/barcode.txt"
barcode_10x="${intermediate_dir}/barcode.txt"

split -l 500 "${barcode_10x}" "${intermediate_dir}/barcode_chunk_" # Splitting barcode file into chunks

# Get the list of samples
samples=$(ls "${intermediate_dir}"/*sub*_merged_homdomain.fasta)

for sample in $samples; do
    base_output_file="${intermediate_dir}/$(basename "${sample%.fasta}")"
    parallel -j 12 seqkit grep -s -r -m 0 -j 4 -f {} "$sample" -o "{}_$(basename "$sample")" ::: "${intermediate_dir}"/barcode_chunk_* # Running seqkit grep in chunk to parallelize
    cat "${intermediate_dir}"/barcode_chunk_*_"$(basename "${sample}")" > "${base_output_file}_10x_barcode_combined.fasta" # Combining output files
    rm "${intermediate_dir}"/barcode_chunk_*_"$(basename "${sample}")" ## cleaning the intermediate files

    # converting 10x barcode combined.fasta seqtk without a line break. Since in the next step we have extract out 10x BC, UMI Static, Mutable barcode
    seqtk seq -l 0 "${base_output_file}_10x_barcode_combined.fasta" > "${base_output_file}_10x_barcode_combined_nolinebreak.fasta"
done

echo "All seqkit grep commands have completed."