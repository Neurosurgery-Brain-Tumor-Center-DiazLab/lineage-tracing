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

# Check if the required arguments are provided
if [ -z "${intermediate_dir}" ] || [ -z "${cellranger_dir}" ]; then
    echo "Missing required arguments"
    usage
fi

# Extracting the 10x barcode from the cellranger output
# Remove '-1' from barcodes and add '^' to grep from the start
zcat "${cellranger_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" | sed 's/-1//g' | awk '{print "^"$0}' > "${intermediate_dir}/barcode.txt"
barcode_10x="${intermediate_dir}/barcode.txt"

# Split the barcode file into chunks
split -l 500 "${barcode_10x}" "${intermediate_dir}/barcode_chunk_"

# Get the list of samples
samples=$(ls "${intermediate_dir}"/*sub*_merged_homdomain.fasta)

# Run seqkit grep for each sample
for sample in $samples; do
    # Define the base output filename
    base_output_file="$(basename "${sample%.fasta}")"

    # Run parallel seqkit grep
    parallel -j 12 seqkit grep -s -r -m 0 -j 4 -f {} "$sample" -o "{}_$base_output_file.fasta" ::: "${intermediate_dir}"/barcode_chunk_*

    # Combine output files
    cat "${intermediate_dir}"/barcode_chunk_*_"${base_output_file}".fasta > "${intermediate_dir}"/"${base_output_file"}_10x_barcode_chunks_combined.fasta
    
    # Optionally, clean up the intermediate chunk files
    rm "${intermediate_dir}"/barcode_chunk_*_"${base_output_file}".fasta
done

echo "All seqkit grep commands have completed."
