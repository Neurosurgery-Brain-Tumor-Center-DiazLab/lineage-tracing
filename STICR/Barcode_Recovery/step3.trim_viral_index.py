##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import subprocess
import path_programs_and_files as path

thread = path.thread
memory = path.memory

inputdir = path.EXTRACTED_DIR
outputdir = path.TRIMMED_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

cutadapt = path.CUTADAPT

read2_fastqs = dict()
fastqs = os.listdir(inputdir)
fastqs = list(filter(lambda x: (".R2.extracted.fastq.gz" in x), fastqs))
for fastq in fastqs:
	sample = fastq.split('.')[0]
	read2_fastq = inputdir + fastq
	read2_fastqs[sample] = read2_fastq
	#print (sample, read2_fastq)

samples = list(read2_fastqs.keys())
samples.sort()

# Split READ2 by STICR viral index
for sample in samples:
        read2_fastq = read2_fastqs[sample]
        command = cutadapt + " -j " + str(thread) + " -g " + path.VIRAL_INDEX + " -O 8 --action none -o " + outputdir + sample + ".cutadapt.trimmed.fastq --untrimmed-output " + outputdir + sample + ".cutadapt.untrimmed.fastq " + read2_fastq
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print("Command succeeded:", result.stdout)
        else:
            print("Command failed:", result.stderr)

