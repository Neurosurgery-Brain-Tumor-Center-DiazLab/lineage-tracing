##conda activate base
#!/c4/home/bhyu0217/anaconda3/bin/python

import os, sys
import path_programs_and_files as path

inputdir = path.ALIGNED_DIR
outputdir = path.COLLAPSED_DIR
if not os.path.exists(outputdir): os.makedirs(outputdir)

class bit:
	def __init__(self, read_id, flag, bit_id, pos):
		self.__read_id = read_id
		self.__flag = flag
		self.__bit_id = bit_id
		self.__pos = pos
	def read_id(self): return (self.__read_id)
	def flag(self): return (self.__flag)
	def bit_id(self): return (self.__bit_id)
	def pos(self): return (self.__pos)

def make_dict(samfile):
	bit_dict = dict()
	samfile = open(samfile, 'r')
	readlines = samfile.readlines(); samfile.close()
	for readline in readlines:
		if not readline.startswith('@'):
			readline = readline.strip().split('\t')
			read_id = '_'.join([readline[0].split('_')[0], readline[0].split('_')[2], readline[0].split('_')[1]])
			flag = readline[1]
			bit_id = readline[2]
			pos = readline[5].split('S')[0]
			#print (read_id, cell_barcode, umi, flag, bit_id, pos)
			bit_dict[readline[0]] = bit(read_id, flag, bit_id, pos)
	return (bit_dict)

sample_name = sys.argv[1]
bit1_dict = make_dict(inputdir + sample_name + "_Bit1/Aligned.out.sam")
bit2_dict = make_dict(inputdir + sample_name + "_Bit2/Aligned.out.sam")
bit3_dict = make_dict(inputdir + sample_name + "_Bit3/Aligned.out.sam")
print(sample_name, len(bit1_dict.keys()), len(bit2_dict.keys()), len(bit3_dict.keys()))

bit_dict = dict()
for bit_id in bit1_dict.keys():
	if bit_id in bit2_dict.keys() and bit_id in bit3_dict.keys():
		bit_dict[bit_id] = [bit1_dict[bit_id], bit2_dict[bit_id], bit3_dict[bit_id]]
print(sample_name, len(bit_dict.keys()))

new_bit_dict = dict()
for bit_id in bit_dict.keys():
	bit1 = bit_dict[bit_id][0]; bit2 = bit_dict[bit_id][1]; bit3 = bit_dict[bit_id][2]
	if bit1.flag() == "0" and bit2.flag() == "0" and bit3.flag() == "0":
		bit_hit = '+'.join([bit1.bit_id(), bit2.bit_id(), bit3.bit_id()])
		if not bit_hit in new_bit_dict.keys(): new_bit_dict[bit_hit] = []
		new_bit_dict[bit_hit].append(bit1.read_id())
	else:
		#print (bit_id, bit1.bit_id(), bit1.flag(), bit2.bit_id(), bit2.flag(), bit3.bit_id(), bit3.flag())
		pass

outfile = open(outputdir + sample_name + ".Bit_hits.flat.txt", 'w')
bit_keys = list(new_bit_dict.keys())
bit_keys.sort()
for bit_key in bit_keys:
	reads = new_bit_dict[bit_key]; reads.sort()
	for read in reads: outfile.write(read+'\t'+bit_key+'\n')

outfile.close()

