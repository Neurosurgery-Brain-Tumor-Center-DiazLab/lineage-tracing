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

# Check if the required argument is provided
if [ -z "${intermediate_dir}" ] || [ -z "${cellranger_dir}" ]; then
    echo "Missing required arguments"
    usage
fi

### Extracting the 10x barcode from the cellranger output
## since cellranger out barcode.tsv file has -1 at end of the cellname so to match it we have to remove it.
## We are also adding ^ to grep it from starting. We would be able to capture the sequence for R2
zcat "${cellranger_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" | sed 's/-1//g' | awk '{print "^"$0}' > "${intermediate_dir}/barcode.txt"
barcode_10x="${intermediate_dir}/barcode.txt"

# Get the list of samples
samples=$(ls "${intermediate_dir}"/*sub*_homdomain.fastq.gz)

# Run seqkit grep with parameter s for sequence, r regex, m mismatch, j threads, f file, o output
for sample in $samples; do
    output_file="${sample%.fastq.gz}_10x_barcode.fastq.gz"
    seqkit grep -s -r -m 0 -j 12 -f "${barcode_10x}" "${sample}" -o "${output_file}" &
done

# Wait for all background jobs to complete
wait

echo "All seqkit grep commands have completed."
