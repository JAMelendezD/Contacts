import numpy as np
import argparse

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Input text file')
parser.add_argument('fmap', type=str, help='Map with value and residue')
parser.add_argument('chain', type=str, help='Chain identifier to match residues')
parser.add_argument('out', type=str, help='Output file without extension')
args = parser.parse_args()

def splitm(line):
	'''
	Correctly split pdb file
	'''
	return([line[0:6].strip(),line[6:11].strip(),line[12:16].strip(),line[16:17].strip(),line[17:20].strip(),
			line[21:22].strip(),line[22:26].strip(),line[26:27].strip(),line[30:38].strip(),line[38:46].strip(),
			line[46:54].strip(),line[54:60].strip(),line[60:66].strip(),line[76:78].strip(),line[78:80].strip()])


def create_dict(fname):
	map = {}
	with open(fname,'r') as f:
		for line in f:
			data = line.split()
			map[data[0][3:]] = float(data[1])*100
	return(map)

def run(inp,out,fmap,chain):
	print('Suggested spectrum')
	print('spectrum b, 0xF4F3F3 0xD28288 0xF6A377 0xFBDF66 0xCFF66A 0x77FB74')
	map = create_dict(fmap)
	keys = map.keys() 
	pdb_format = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
	with open(inp,'r') as f:
		with open(out, 'w') as fw:
			for line in f:
				data = splitm(line)
				if data[0] == 'ATOM':
					data[1] = int(data[1])
					data[6] = int(data[6])
					data[8] = float(data[8])
					data[9] = float(data[9])
					data[10] = float(data[10])
					data[11] = float(data[11])
					if data[5] == chain:
						if str(data[6]) in keys:
							data[12] = map[str(data[6])]
							fw.write(pdb_format.format(*data))
						else:
							if len(data[12]) == 0:
								data[12] = 0.0
								fw.write(pdb_format.format(*data))
							else:
								fw.write(line)
					else:
						if len(data[12]) == 0:
							data[12] = 0.0
							fw.write(pdb_format.format(*data))
						else:
							fw.write(line)
				else:
					fw.write(line)

if __name__ == '__main__':
	run(args.input,args.out,args.fmap,args.chain)