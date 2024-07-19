from Bio import SeqIO

def extract_sequences(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        # Extract first 16 bases
        first_16 = seq[:16]
        # Extract next 12 bases
        next_12 = seq[16:28]
        # Extract sequence between patterns
        start_index = seq.find("AGCTTGAATTC") + len("AGCTTGAATTC")
        end_index = seq.find("TTTGCT")
        between_patterns = seq[start_index:end_index]
        # Extract sequence including and after "TTTGCT"
        after_pattern = seq[seq.find("TTTGCT"):]

        # Output or save the results as needed
        print(f">First_16\n{first_16}")
        print(f">Next_12\n{next_12}")
        print(f">Between_Patterns\n{between_patterns}")
        print(f">After_Pattern\n{after_pattern}")

# Replace 'input.fasta' with your actual file
extract_sequences("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007_unpaired_R1_homdomain_10x_barcode_2.fasta")

from Bio import SeqIO
import re
import csv

def extract_sequences(seq):
    # Extract the first 16 bases
    ten_x_barcode = seq[:16]
    
    # Extract the next 12 bases
    umi = seq[16:28]
    
    # Extract the bases between "AGCTTGAATTC" and "TTTGCT"
    agctttgaattc = "GAATTC"
    tttgct = "TTTGCT"
    
    start_static = seq.find(agctttgaattc)
    end_static = seq.find(tttgct, start_static + len(agctttgaattc))
    
    if start_static != -1 and end_static != -1:
        start_static += len(agctttgaattc)
        static = seq[start_static:end_static]
    else:
        static = ""
    
    # Extract all bases after "TTTGCT"
    if end_static != -1:
        mutable_barcode = seq[end_static:]
    else:
        mutable_barcode = ""
    
    return ten_x_barcode, umi, static, mutable_barcode

# Input FASTA file
input_file = "/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007_unpaired_R1_homdomain_10x_barcode.fasta"

# Output CSV file
output_file = "/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/output.csv"

# Open the output CSV file for writing
with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['10x_barcode', 'UMI', 'static', 'mutable_barcode']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    # Parse the FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq)
        ten_x_barcode, umi, static, mutable_barcode = extract_sequences(seq)
        
        # Write to CSV
        writer.writerow({
            '10x_barcode': ten_x_barcode,
            'UMI': umi,
            'static': static,
            'mutable_barcode': mutable_barcode
        })

