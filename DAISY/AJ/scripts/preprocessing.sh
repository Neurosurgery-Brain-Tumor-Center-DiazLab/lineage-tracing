#!/bin/bash
# Function to display usage

usage() {
    echo "Usage: $0 -i <input_directory> -m <intermediate_directory> -o <output_directory>"
    exit 1
}

# Parse command-line arguments
while getopts ":i:m:o:" opt; do
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
if [ -z "${input_dir}" ] || [ -z "${intermediate_dir}" ] || [ -z "${output_dir}" ]; then
    echo "Missing required arguments"
    usage
fi

# Activate the conda environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda init
# conda activate daisy

# Check if fastp is available
if ! command -v fastp &> /dev/null; then
    echo "fastp could not be found. Please make sure it is installed in the conda environment."
    exit 1
fi

# Creating the output directories
echo "Making the result and intermediate directory"
mkdir -p "$intermediate_dir"
mkdir -p "$output_dir"

# Getting the sample names
samples=$(ls ${input_dir}/*sub*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//g' | sed "s:${input_dir}/::g" | sort | uniq)

# Fastp processing
for sample in $samples; do
  echo "Processing $sample"

  # Input
  R1="${sample}_R1_001.fastq.gz"
  R2="${sample}_R2_001.fastq.gz"

  # Unmerged output
  QC_R1="${sample}_QC_R1.fastq.gz"
  QC_R2="${sample}_QC_R2.fastq.gz"
  
  # Merged output
  merged="${sample}_merged.fastq.gz"
  unmerged_R1="${sample}_unmerged_R1.fastq.gz"
  unmerged_R2="${sample}_unmerged_R2.fastq.gz"

  # Unpaired output
  unpaired_R1="${sample}_unpaired_R1.fastq.gz"
  unpaired_R2="${sample}_unpaired_R2.fastq.gz"

  # Failed output
  failed_read="${sample}_fail.fastq.gz"

  # Parameters:
  # -q: Phred Quality > 20 is qualified
  # -u: how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%
  # -n: remove NNNNN
  # -e: average quality > 20
  # Default it perform adapter trimming if you like to disable it by -A or --disable_adapter_trimming
  # --merged_out shouuld be given to specify the file to store merged reads, otherwise you should enable --stdout to stream the merged reads to STDOUT. The merged reads are also filtered.
  # --out1 and --out2 will be the reads that cannot be merged successfully, but both pass all the filters.
  # --unpaired1 will be the reads that cannot be merged, read1 passes filters but read2 doesn't.
  # --unpaired2 will be the reads that cannot be merged, read2 passes filters but read1 doesn't.
  # --include_unmerged can be enabled to make reads of --out1, --out2, --unpaired1 and --unpaired2 redirected to --merged_out. So you will get a single output file. This option is disabled by default.

  # Perform fastp
  fastp -i "${input_dir}/$R1" \
    -I "${input_dir}/$R2" \
    -o "${intermediate_dir}/$QC_R1" \
    -O "${intermediate_dir}/$QC_R2" \
    --unpaired1 "${intermediate_dir}/$unpaired_R1" \
    --unpaired2 "${intermediate_dir}/$unpaired_R2" \
    --failed_out "${intermediate_dir}/$failed_read" \
    --merge --merged_out "${intermediate_dir}/$merged" \
    --out1 "${intermediate_dir}/$unmerged_R1" \
    --out2 "${intermediate_dir}/$unmerged_R2" \
    -q 15 -u 40 -n 5 -e 20

  echo "Quality filtering and merging of ${sample} R1 and R2 completed."
done
