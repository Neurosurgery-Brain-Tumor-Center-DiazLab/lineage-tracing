##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import path_programs_and_files as path

thread = path.THREAD
memory = path.MEMORY

inputdir = path.RAW_SEQUENCING_DIR
whitelist_dir = path.WHITELIST_DIR
outputdir = path.EXTRACTED_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

umi_tools = path.UMI_TOOLS
pattern = path.PATTERN

fastq_files = dict()
fastqs = os.listdir(inputdir)
fastqs = list(filter(lambda x: (".fastq.gz" in x), fastqs))
fastqs.sort()
for i in range(0, len(fastqs), 2):
	sample = fastqs[i].split('.')[0]
	read1_fq = inputdir + fastqs[i]; read2_fq = inputdir + fastqs[i+1]
	fastq_files[sample] = [read1_fq, read2_fq]
	#print (sample, read1_fq, read2_fq)

sbatch_dir = path.SBATCH_DIR + "step2.add_cell_barcodes_UMIs/"
if not os.path.exists(sbatch_dir): os.makedirs(sbatch_dir)

# Extract real cell barcodes and UMIs from READ1 and add it to READ2 read head
for sample in fastq_files.keys():
	read1_fq = fastq_files[sample][0]
	read2_fq = fastq_files[sample][1]
	#print (sample, read1_fq, read2_fq)
	whitelist = whitelist_dir + sample + ".cell_barcode_whitelist.txt"
	if os.path.exists(whitelist):
		sbatch = open(sbatch_dir + sample + '.sh', 'w')
		sbatch.write(path.SBATCH(thread, memory))
		command = umi_tools + ' extract --bc-pattern=' + pattern + ' --stdin ' + read1_fq + ' --stdout ' + outputdir + sample + '.R1.extracted.fastq.gz --read2-in ' + read2_fq + ' --read2-out=' + outputdir + sample + '.R2.extracted.fastq.gz --whitelist=' + whitelist + '\n'
		sbatch.write(command)
		sbatch.close()
	else: print ("Please check cell barcode whitelist: ", sample)

