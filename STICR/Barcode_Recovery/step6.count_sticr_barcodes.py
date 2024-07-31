##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import path_programs_and_files as path

inputdir = path.COLLAPSED_DIR
outputdir = path.COUNT_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

umi_tools = path.UMI_TOOLS

flat_files = dict()
infiles = os.listdir(inputdir)
infiles = list(filter(lambda x: (".Bit_hits.flat.txt" in x), infiles))
for infile in infiles:
	sample = infile.split('.')[0]
	flat_files[sample] = inputdir + infile
	print(sample, infile)

# Count STICR barcodes by CBC and UMI_counts
for sample in flat_files.keys():
        flat_file = flat_files[sample]
        command = umi_tools + " count_tab --per-cell --edit-distance-threshold=1 -I " + flat_file + " -S " + outputdir + sample + "_counts.tsv -L " + outputdir + sample + "_counts.log"
        print(command)
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode == 0: print("Command succeeded:", result.stdout)
        else: print("Command failed:", result.stderr)

