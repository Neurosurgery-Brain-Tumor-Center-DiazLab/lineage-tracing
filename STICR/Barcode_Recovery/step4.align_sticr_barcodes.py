##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import path_programs_and_files as path

thread = path.THREAD
memory = path.MEMORY

inputdir = path.TRIMMED_DIR
outputdir = path.ALIGNED_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

star = path.STAR
star_index = path.STAR_INDEX

read2_fastqs = dict()
fastqs = os.listdir(inputdir)
fastqs = list(filter(lambda x: (".cutadapt.trimmed.fastq" in x), fastqs))
for fastq in fastqs:
	sample = fastq.split('.')[0]
	read2_fastq = inputdir + fastq
	read2_fastqs[sample] = read2_fastq
	#print (sample, read2_fastq)

sbatch_dir = path.SBATCH_DIR + "step4.align_sticr_barcodes/"
if not os.path.exists(sbatch_dir): os.makedirs(sbatch_dir)

# Extract BIT sequences from READ2 and align it to pre-built references
for sample in read2_fastqs.keys():
	read2_fastq = read2_fastqs[sample]
	for bit in ["Bit1", "Bit2", "Bit3"]:
		sbatch = open(sbatch_dir + sample + '.' + bit + '.sh', 'w')
		sbatch.write(path.SBATCH(thread, memory))
		sbatch.write(star + " --runThreadN " + str(thread) + " --genomeDir " + star_index + "STICR_" + bit + "/ --readFilesIn " + read2_fastq + " --outSAMattributes NH HI AS nM NM XS --outFileNamePrefix " + outputdir + sample + "_" + bit + "/ --outSAMtype SAM --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 1 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterMultimapNmax 1 --outReadsUnmapped Fastx\n")
		sbatch.close()

