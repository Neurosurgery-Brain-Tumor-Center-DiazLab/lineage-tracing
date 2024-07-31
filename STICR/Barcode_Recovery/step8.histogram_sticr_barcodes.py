##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os
import path_programs_and_files as path
import matplotlib.pyplot as plt
import numpy as np

inputdir = path.FILTER_DIR
outputdir = path.HISTOGRAM_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

sticr_files = dict()
infiles = os.listdir(inputdir)
infiles = list(filter(lambda x: (".Final_Barcodes.tsv" in x), infiles))
for infile in infiles:
	sample = infile.split('.Final_Barcodes.tsv')[0]
	sticr_files[sample] = inputdir + infile
	print(sample, infile)

# Valid STICR barcodes with QC PASS
samples = list(sticr_files.keys())
samples.sort()
for sample in samples:
	barcode_file = path.BARCODE_DIR + sample + ".Cell_Barcodes.tsv"
	infile = open(barcode_file, 'r')
	readlines = infile.readlines(); infile.close()
	cell_barcodes = []
	for readline in readlines: cell_barcodes.append(readline.strip())
	print(sample, len(cell_barcodes))
	sticr_file = sticr_files[sample]
	infile = open(sticr_file, 'r')
	readlines = infile.readlines(); infile.close()
	sticr_barcodes = dict(); unique_barcodes = []; valid_barcodes = []
	for readline in readlines[1:]:
		readline = readline.strip().split('\t')
		cell_barcode = readline[0]; sticr_barcode = readline[1]
		if not cell_barcode in unique_barcodes:
			unique_barcodes.append(cell_barcode)
		if cell_barcode in cell_barcodes:
			if not cell_barcode in valid_barcodes:
				valid_barcodes.append(cell_barcode)
			if not sticr_barcode in sticr_barcodes.keys():
				sticr_barcodes[sticr_barcode] = []
			sticr_barcodes[sticr_barcode].append(cell_barcode)
	keys = sorted(sticr_barcodes, key = lambda key: len(sticr_barcodes[key]), reverse=True)
	#for key in keys: print (key, len(sticr_barcodes[key]))
	x = []
	for key in sticr_barcodes.keys(): x.append(len(sticr_barcodes[key]))
	par = max(x)
	plt.style.use('ggplot')
	if par <= 10:
		plt.hist(x, bins=10, edgecolor="black")
		plt.xticks(np.arange(0,10,1))
	else:
		plt.hist(x, bins=par, edgecolor="black")
		if par <= 100: plt.xticks(np.arange(0,par,10))
		elif par <= 200: plt.xticks(np.arange(0,par,20))
		elif par <= 300: plt.xticks(np.arange(0,par,30))
		elif par <= 400: plt.xticks(np.arange(0,par,40))
		else: plt.xticks(np.arange(0,par,50))
	plt.savefig(outputdir + sample + ".hist_clone.pdf", format = "pdf", bbox_inches = "tight")
	plt.show()
	plt.clf()

