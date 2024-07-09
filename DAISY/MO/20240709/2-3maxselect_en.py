import sys
import pandas as pd

def main(input_file, output_file):

    df = pd.read_csv(input_file, sep='\t')

    result = df.loc[df.groupby('tenx')['overlapmutate'].idxmax()]

    result.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python maxselect.py input.txt output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)