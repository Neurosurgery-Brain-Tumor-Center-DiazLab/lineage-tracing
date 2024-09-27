#!/bin/bash
# Function to display usage

usage() {
    echo "Usage: $0 -s <sampleID> -r <Read1_fastq> -R <Read2_fastq> -m <intermediate_directory> -o <output_directory>"
    exit 1
}

# Parse command-line arguments
while getopts ":s:r:R:m:o:" opt; do
  case ${opt} in
    s )
      samplename=$OPTARG
      ;;
    r )
      R1_fastq=$OPTARG
      ;;
    R )
      R2_fastq=$OPTARG
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

# Debugging output to check variable assignments
echo "Performing Quality for"
echo "Sample Name: ${samplename}"
echo "Read1 Fastq: ${R1_fastq}"
echo "Read2 Fastq: ${R2_fastq}"
echo "Intermediate Dir: ${intermediate_dir}"
echo "Output Dir: ${output_dir}"

# Check if all required arguments are provided
if [ -z "${intermediate_dir}" ] || [ -z "${output_dir}" ] || [ -z "${samplename}" ] || [ -z "${R1_fastq}" ] || [ -z "${R2_fastq}" ]; then
    echo "Missing required arguments"
    usage
fi

# Create necessary directories
mkdir -p "${intermediate_dir}/${samplename}"
mkdir -p "${output_dir}/${samplename}"

# Fastp processing
echo "Processing ${samplename}"

# Check if fastp is available
if ! command -v fastp &> /dev/null; then
    echo "fastp could not be found. Please make sure it is installed in the conda environment."
    exit 1
fi

# Getting the sample names
sample=$(ls $R1_fastq | sed 's:/:\t:g' | awk '{print $NF}' | sed 's/_R1_001.fastq.gz//g')

# Fastp processing
echo "Processing $samplename"

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

# Perform fastp
# fastp -i "${R1_fastq}" \
#     -I "${R2_fastq}" \
#     -o "${intermediate_dir}/${samplename}/${QC_R1}" \
#     -O "${intermediate_dir}/${samplename}/${QC_R2}" \
#     --unpaired1 "${intermediate_dir}/${samplename}/${unpaired_R1}" \
#     --unpaired2 "${intermediate_dir}/${samplename}/${unpaired_R2}" \
#     --failed_out "${intermediate_dir}/${samplename}/${failed_read}" \
#     --merge --merged_out "${intermediate_dir}/${samplename}/${merged}" \
#     --out1 "${intermediate_dir}/${samplename}/${unmerged_R1}" \
#     --out2 "${intermediate_dir}/${samplename}/${unmerged_R2}" \
#     -q 15 -u 40 -n 5 -e 20

fastp -i "${R1_fastq}" \
    -I "${R2_fastq}" \
    -o "${QC_R1}" \
    -O "${QC_R2}" \
    --unpaired1 "${unpaired_R1}" \
    --unpaired2 "${unpaired_R2}" \
    --failed_out "${failed_read}" \
    --merge --merged_out "${merged}" \
    --out1 "${unmerged_R1}" \
    --out2 "${unmerged_R2}" \
    -q 15 -u 40 -n 5 -e 20


echo "Quality filtering and merging of ${sample} R1 and R2 completed."
