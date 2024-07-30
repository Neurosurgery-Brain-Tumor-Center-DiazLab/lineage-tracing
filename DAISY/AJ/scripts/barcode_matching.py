from Bio import SeqIO

# Load barcodes
with open('/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/barcode.txt', 'r') as barcode_file:
    barcodes = [line.strip() for line in barcode_file]

# Open FASTA files
input_fasta = '/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007_merged_homdomain.fasta'
output_fasta = '/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007_merged_homdomain_10x_barcode_biopython.fasta'

# Create a set of barcodes for quick lookup
barcode_set = set(barcodes)

# Open output file
with open(output_fasta, 'w') as out_file:
    # Iterate over sequences in the FASTA file
    for record in SeqIO.parse(input_fasta, 'fasta'):
        # Check if the sequence starts with any of the barcodes
        if any(record.seq.startswith(barcode) for barcode in barcode_set):
            # Write the record to the output file
            SeqIO.write(record, out_file, 'fasta')

print(f"Filtering complete. Results saved to {output_fasta}.")
