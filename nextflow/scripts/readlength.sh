#!/bin/bash
# Function to display usage information

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

# Get the list of samples
samples=$(ls ${intermediate_dir}/*sub*_R*.fastq.gz)

# Run seqkit watch in parallel
for sample in $samples; do
  output_file="${sample%.fastq.gz}_readlength.png"
  seqkit watch --fields ReadLen ${sample} -O ${output_file} &
done

# Wait for all background jobs to complete
wait

echo "All seqkit watch commands have completed."