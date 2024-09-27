import os
import shutil
from collections import Counter

def main(input_file, output_folder):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(input_file, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().split(',')
    data = [line.strip().split(',') for line in lines[1:]]

    static_index = header.index('static')

    static_count = Counter([line[static_index] for line in data])

    sorted_statics = [static for static, _ in static_count.most_common()]

    header.append('BCnumber')
    for i in range(len(data)):

        static = data[i][static_index]
        data[i].append(f"BC{sorted_statics.index(static) + 1}")


    bc1_data = [line for line in data if line[-1] == 'BC1']
    bc1_file = os.path.join(output_folder, "BC1.txt")
    with open(bc1_file, 'w') as f:

        f.write(','.join(header) + '\n')

        for line in bc1_data:
            f.write(','.join(line) + '\n')

    bcall_file = os.path.join(output_folder, "BCall.txt")
    with open(bcall_file, 'w') as f:

        f.write(','.join(header) + '\n')

        for line in data:
            f.write(','.join(line) + '\n')


    for i, static in enumerate(sorted_statics, start=1):

        static_data = [line for line in data if line[static_index] == static]

        output_file = os.path.join(output_folder, f"BC{i}.txt")
        with open(output_file, 'w') as f:

            f.write(','.join(header) + '\n')

            for line in static_data:
                f.write(','.join(line) + '\n')

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python BCsort.py input_file output_folder")
        sys.exit(1)

    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    main(input_file, output_folder)
