##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import path_programs_and_files as path

thread = path.THREAD
memory = path.MEMORY

inputdir = path.RAW_SEQUENCING_DIR
outputdir = path.WHITELIST_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

umi_tools = path.UMI_TOOLS
expected_cell_number = path.EXPECTED_CELL_NUMBER
pattern = path.PATTERN

read1_fastqs = dict()
fastqs = os.listdir(inputdir)
fastqs = list(filter(lambda x: ("_R1.001.fastq.gz" in x), fastqs))
for fastq in fastqs:
	sample = fastq.split('.')[0]
	read1_fastq = inputdir + fastq
	read1_fastqs[sample] = read1_fastq
	#print (sample, read1_fastq)

sbatch_dir = path.SBATCH_DIR + "step1.make_cell_barcode_whitelist/"
if not os.path.exists(sbatch_dir): os.makedirs(sbatch_dir)

# Find cell barcode whitelist with expected cell number
for sample in read1_fastqs.keys():
	read1_fastq = read1_fastqs[sample]
	#print (sample, read1_fastq)
	sbatch = open(sbatch_dir + sample + '.sh', 'w')
	sbatch.write(path.SBATCH(thread, memory))
	command = umi_tools + ' whitelist --stdin ' + read1_fastq + ' --bc-pattern=' + pattern + ' --set-cell-number=' + str(path.EXPECTED_CELL_NUMBER[sample]) + ' --log2stderr > ' + outputdir + sample + '.cell_barcode_whitelist.txt\n'
	sbatch.write(command)
	sbatch.close()

