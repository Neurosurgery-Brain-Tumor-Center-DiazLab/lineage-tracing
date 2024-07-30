
import argparse
import sys
import pandas as pd

def main(input_file, output_file):
    # Read the file
    df = pd.read_csv(input_file, delimiter=',', header=None, names=["Total Number", "10x Barcode", "Static Barcode", "Mutable Barcode"])

    # Identify the maxium UMI for the 10x BC
    filtered_df = df.loc[df.groupby("10x Barcode")["Total Number"].idxmax()]
    
    # Save the file
    filtered_df.to_csv(output_file, sep=',', index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select the maximum 10x barcode")
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_file", help="Path to the output CSV file")
    
    args = parser.parse_args()
    main(args.input_file, args.output_file)


