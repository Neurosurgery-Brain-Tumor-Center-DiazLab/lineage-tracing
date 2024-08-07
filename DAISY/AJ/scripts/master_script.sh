#!/bin/bash

# Exit if any command fails
set -e

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
"$script_dir/preprocessing.sh" -i "$input_dir" -m "$intermediate_dir" -o "$output_dir"

# # Check if preprocessing step was successful
if [ $? -ne 0 ]; then
    echo "Preprocessing step failed. Exiting."
    exit 1
fi

# Run step 2: Read Length Step
"$script_dir/readlength.sh" -m "$intermediate_dir"

# # Check if read length step was successful
if [ $? -ne 0 ]; then
    echo "Read length step failed. Exiting."
    exit 1
fi
echo "Read length step completed successfully."

# Run step 3: grepping the TAGTTACGCCAAGCTTGAATTC in the merged, paired, unpaired files
# AD: the homodomain sequence should be a parameter
"$script_dir/homo_match.sh" -m "$intermediate_dir"

# # # Check if homo domain step was successful
if [ $? -ne 0 ]; then
    echo "Homodomain match step failed. Exiting."
    exit 1
fi
echo "Homo domain step completed successfully."

# Run Step 4: Fastq to Fasta
# Before barcode matching that will increase the speed of the seqkit grep using seqtk
"$script_dir/homo_fq2fa.sh" -m "$intermediate_dir"

# # # Check if homo domain step was successful
if [ $? -ne 0 ]; then
    echo "Homodomain fastq 2 fasta step failed. Exiting."
    exit 1
fi
echo "All steps completed successfully."


# Run step 5: Selecting the 10x barcode sequence
"$script_dir/barcode_matching_chunks2.sh" -m "$intermediate_dir" -c "${cellranger_dir}"

# # # Check if homo domain step was successful
if [ $? -ne 0 ]; then
    echo "barcode matching step failed. Exiting."
    exit 1
fi

# Run the Python script for each file ending with "merged_homdomain_10x_barcode_combined.fasta"
for fasta_file in "${intermediate_dir}"/*merged_homdomain_10x_barcode_combined.fasta; do
    python3 "$script_dir/sequence_extract.py" "$fasta_file" "${fasta_file%.fasta}_10x_UMI_static_mut.csv"
    
    # Cleaning the file
    awk -F, 'NF==4 && $1!="" && $2!="" && $3!="" && $4!=""' "${fasta_file%.fasta}_10x_UMI_static_mut.csv" > "${fasta_file%.fasta}_10x_UMI_static_mut_clean.csv"

    Removing the PCR Bias
    sort "${fasta_file%.fasta}_10x_UMI_static_mut_clean.csv" | uniq -c | awk -F " " '{if ($(NF-1) > 2) print $0}' > "${fasta_file%.fasta}_10x_UMI_static_mut_clean_remove_PCR.csv"

    Selecting the 10x barcode, Static and Daisy column
    cut -d "," -f1,3,4 "${fasta_file%.fasta}_10x_UMI_static_mut_clean_remove_PCR.csv" | awk '{print $NF}' | sort | uniq -c | awk -F " " '{print $(NF-1)","$NF}' > "${fasta_file%.fasta}_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv"

    Selecting the maximum 10x barcode
    python3 "$script_dir/max_selection.py" "${fasta_file%.fasta}_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv"

    Global Alignment nucleotide 
    python3 "$script_dir/alignment_ATGCindel.py" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_alignment.csv"

    Global Alignment 0 and deletion
    python3 "$script_dir/alignment0_ATGC.py" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_0_alignment.csv"

    Global Alignment 0, 1 and deletion
    python3 "$script_dir/alignment_0_1_ATGC.py" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_0_1_alignment.csv"

    BC sorting
    python3 "$script_dir/BCsort.py" "${fasta_file%.fasta}_max_10x_UMI_static_mut_clean_remove_PCR_umi_col_ATGC_alignment.csv" "${fasta_file%_merged_homdomain_10x_barcode_combined.fasta}"
done

# Function to process a BC file
process_bc_file() {
  BC=$1
  script_dir=$2
  echo "Processing BC file: $BC"
  python3 "$script_dir/njmodified.py" "$BC" --output "${BC%.txt}.newick" || echo "Error running njmodified.py on $BC"
  python3 "$script_dir/txt2nexus.py" "$BC" || echo "Error running txt2nexus.py on $BC"
}

export -f process_bc_file

# Loop through fasta files and generate commands for parallel execution
for fasta_file in "${intermediate_dir}"/*_merged_homdomain_10x_barcode_combined.fasta; do
  sample_dir="${fasta_file%_merged_homdomain_10x_barcode_combined.fasta}"
  echo "Processing phylogeny for file: $sample_dir"
  find "$sample_dir" -name "BC*.txt" | parallel -j 20 process_bc_file {} "$script_dir"
done

echo "All steps completed successfully."
