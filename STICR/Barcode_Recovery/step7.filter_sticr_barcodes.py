##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import path_programs_and_files as path

import matplotlib.pyplot as plt
import numpy as np

inputdir = path.COUNT_DIR
outputdir = path.FILTER_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

count_files = dict()
infiles = os.listdir(inputdir)
infiles = list(filter(lambda x: ("_counts.tsv" in x), infiles))
for infile in infiles:
	sample = infile.split('_counts.tsv')[0]
	count_files[sample] = inputdir + infile
	#print(sample, infile)

# Filter STICR barcodes with UMI threshold
samples = list(count_files.keys())
samples.sort()
for sample in samples:
	print(sample)
	outfile = open(outputdir + sample + ".Final_Barcodes.tsv", 'w')
	outfile.write("CBC\tbarcode\tUMI_Count\n")
	count_file = count_files[sample]
	infile = open(count_file, 'r')
	readlines = infile.readlines(); infile.close()
	x = []; y = []
	for readline in readlines[1:]:
		readline = readline.strip().split('\t')
		cell_barcode = readline[0].split("'")[1]
		sticr_barcode = readline[1]; umi_counts = readline[2]
		x.append(int(umi_counts))
		if int(umi_counts) <= 20: y.append(int(umi_counts))
		if int(umi_counts) >= path.UMI_THRESHOLD:
			outfile.write(cell_barcode+'\t'+sticr_barcode+'\t'+umi_counts+'\n')
	outfile.close()

	hist_dir = outputdir + "histogram_umi/"
	if not os.path.exists(hist_dir): os.makedirs(hist_dir)

	par = max(x)
	plt.style.use('ggplot')
	plt.hist(x, bins=par, edgecolor="black")
	plt.xticks(np.arange(0,par,10))
	plt.savefig(hist_dir + sample + ".hist_max.pdf", format = "pdf", bbox_inches = "tight")
	plt.show()
	plt.clf()

	par = max(y)
	plt.style.use('ggplot')
	plt.hist(y, bins=par, edgecolor="black")
	plt.xticks(np.arange(0,21,1))
	plt.axvline(path.UMI_THRESHOLD, color='k', linestyle='dotted', linewidth=2)
	plt.savefig(hist_dir + sample + ".hist_max20.pdf", format = "pdf", bbox_inches = "tight")
	plt.show()
	plt.clf()

