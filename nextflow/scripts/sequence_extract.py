import argparse
from Bio import SeqIO
import csv
import regex as re

def extract_sequences(seq):
    # Extract the first 16 bases
    ten_x_barcode = seq[:16]
    
    # Extract the next 12 bases
    umi = seq[16:28]
    
    # Regex pattern with up to one mismatch for "AGCTTGAATTC"
    pattern_agctttgaattc = "(AGCTTGAATTC){s<=1}"
    
    # Regex pattern with up to one mismatch for "TTTGCT"
    pattern_tttgct = "(TTTGCT){s<=1}"
    
    # Find the match for "AGCTTGAATTC" with up to one mismatch
    match_agctttgaattc = re.search(pattern_agctttgaattc, seq)
    if match_agctttgaattc:
        start_static = match_agctttgaattc.end()
    else:
        start_static = -1
    
    # Find the match for "TTTGCT" with up to one mismatch, starting after the first pattern
    match_tttgct = re.search(pattern_tttgct, seq[start_static:])
    if match_tttgct:
        end_static = match_tttgct.start() + start_static
    else:
        end_static = -1
    
    # Extract the bases between "AGCTTGAATTC" and "TTTGCT" with up to one mismatch, allowing 9-11 bases
    if start_static != -1 and end_static != -1:
        static_length = end_static - start_static
        if 9 <= static_length <= 11:
            static = seq[start_static:end_static]
        else:
            static = ""
    else:
        static = ""
    
    # Extract all bases after "TTTGCT", including "TTTGCT"
    if end_static != -1:
        mutable_barcode = seq[end_static:]  # Including "TTTGCT"
    else:
        mutable_barcode = ""
    
    return ten_x_barcode, umi, static, mutable_barcode

# Input FASTA file
# input_file = "/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007_merged_homdomain_10x_barcode_combined_nolinebreak.fasta"

# Output CSV file
# output_file = "/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007_merged_homdomain_10x_barcode_combined_nolinebreak_static_mutable.csv"

# Open the output CSV file for writing
def main(input_file, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['10x_barcode', 'UMI', 'static', 'mutable_barcode']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for record in SeqIO.parse(input_file, "fasta"):
            seq = str(record.seq)
            ten_x_barcode, umi, static, mutable_barcode = extract_sequences(seq)
            
            writer.writerow({
                '10x_barcode': ten_x_barcode,
                'UMI': umi,
                'static': static,
                'mutable_barcode': mutable_barcode
                })

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences from FASTA file")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_file", help="Path to the output CSV file")
    
    args = parser.parse_args()
    main(args.input_file, args.output_file)
