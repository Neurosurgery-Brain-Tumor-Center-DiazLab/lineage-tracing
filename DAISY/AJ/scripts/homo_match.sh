#!/bin/bash
## We have find the reads with match "TAGTTACGCCAAGCTTGAATTC" 

usage() {
    echo "Usage: $0 -i <intermediate_directory>"
    exit 1
}

# Parse command-line arguments
while getopts ":m:" opt; do
  case ${opt} in
    m )
      intermediate_dir=$OPTARG
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
if [ -z "${intermediate_dir}" ]; then
    echo "Missing required argument"
    usage
fi

# Getting the sample names
samples=$(ls ${intermediate_dir}/*sub*ed*.fastq.gz)

# Parameters:
#   -d, --degenerate                  pattern/motif contains degenerate base
#       --delete-matched              delete a pattern right after being matched, this keeps the firstly
#                                     matched data and speedups when using regular expressions
#   -m, --max-mismatch int            max mismatch when matching by seq. For large genomes like human
#                                     genome, using mapping/alignment tools would be faster
#   -p, --pattern strings             search pattern (multiple values supported. Attention: use double
#                                     quotation marks for patterns containing comma, e.g., -p '"A{2,}"')
for sample in $samples; do
    output_file="${sample%.fastq.gz}_homdomain.fastq.gz"
    seqkit grep -s -m 2 -j 12 -p "TAGTTACGCCAAGCTTGAATTC" ${sample} -o ${output_file} &
done

# Wait for all background jobs to complete
wait

echo "All seqkit grep commands have completed."