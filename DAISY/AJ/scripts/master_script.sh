#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_directory> -m <intermediate_directory> -o <output_directory> -c <cellranger_directory>"
    exit 1
}

# Parse command-line arguments
while getopts ":i:m:o:c:" opt; do
  case ${opt} in
    i )
      input_dir=$OPTARG
      ;;
    m )
      intermediate_dir=$OPTARG
      ;;
    o )
      output_dir=$OPTARG
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

# Check if all required arguments are provided
if [ -z "${input_dir}" ] || [ -z "${intermediate_dir}" ] || [ -z "${output_dir}" ] || [ -z "${cellranger_dir}" ]; then
    echo "Missing required arguments"
    usage
fi

# Determine the script directory
script_dir=$(dirname "$0")

# Create directories if they don't exist
mkdir -p "$intermediate_dir"
mkdir -p "$output_dir"

# Run step 1: Preprocessing Step
# "$script_dir/preprocessing.sh" -i "$input_dir" -m "$intermediate_dir" -o "$output_dir"

# # Check if preprocessing step was successful
# if [ $? -ne 0 ]; then
#     echo "Preprocessing step failed. Exiting."
#     exit 1
# fi

# Run step 2: Read Length Step
# "$script_dir/readlength.sh" -m "$intermediate_dir"

# # Check if read length step was successful
# if [ $? -ne 0 ]; then
#     echo "Read length step failed. Exiting."
#     exit 1
# fi
# echo "All steps completed successfully."

# Run step 3: grepping the TAGTTACGCCAAGCTTGAATTC in the merged, paired, unpaired files
# "$script_dir/homo_match.sh" -m "$intermediate_dir"

# # # Check if homo domain step was successful
# if [ $? -ne 0 ]; then
#     echo "Homodomain match step failed. Exiting."
#     exit 1
# fi
# echo "All steps completed successfully."

# Run step 4: Selecting the 10x barcode sequence
"$script_dir/barcode_matching.sh" -m "$intermediate_dir" -c "${cellranger_dir}"

# # Check if homo domain step was successful
if [ $? -ne 0 ]; then
    echo "barcode matching step failed. Exiting."
    exit 1
fi
echo "All steps completed successfully."
