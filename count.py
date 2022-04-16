import numpy as np
import argparse

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Input text file')
parser.add_argument('fix_sele1', type=int, help='Number to add to the default resnum for the first selection')
parser.add_argument('fix_sele2', type=int, help='Number to add to the default resnum for the second selection')
parser.add_argument('total_frames', type=int, help='Total number of frames')
parser.add_argument('out', type=str, help='Output file without extension')
parser.add_argument('--mode', default=0,required=False, type=int, help='0 counts individual. 1 counts pairs.')
args = parser.parse_args()


def count_pair(fname,fix1,fix2):
	'''
	Count to get probability of a pair interacting 
	'''
	dict_count = {}
	previous_frame = 1e10
	with open(fname,"r") as f:
		for line in f:
			data = line.split()
			frame = data[0]
			if data[4] == 'A':
				pair = data[1]+str(int(data[2])+(fix1-1))+'-'+data[5]+str(int(data[6])+(fix2-1))
			elif data[4] == 'B':
				pair = data[5]+str(int(data[6])+(fix1-1))+'-'+data[1]+str(int(data[2])+(fix2-1))

			if frame != previous_frame:
				previous_frame = frame
				pair_in_frame = []
			if pair not in pair_in_frame:
				try:
					dict_count[pair] += 1
					pair_in_frame.append(pair)
				except:
					dict_count[pair] = 1
					pair_in_frame.append(pair)
	return(dict_count)


def count_ind(fname,fix1,fix2):
	'''
	Count to get probability of individual at least interacting once
	'''
	dict_count1 = {}
	dict_count2 = {}
	previous_frame = 1e10
	with open(fname,"r") as f:
		for line in f:
			data = line.split()
			frame = data[0]
			if data[4] == 'A':
				name1 = data[1]+str(int(data[2])+(fix1-1))
				name2 = data[5]+str(int(data[6])+(fix2-1))
			elif data[4] == 'B':
				name1 = data[5]+str(int(data[6])+(fix1-1))
				name2 = data[1]+str(int(data[2])+(fix2-1))
			if frame != previous_frame:
				previous_frame = frame
				ind1_in_frame = []
				ind2_in_frame = []
			if name1 not in ind1_in_frame:
				try:
					dict_count1[name1] += 1
					ind1_in_frame.append(name1)
				except:
					dict_count1[name1] = 1
					ind1_in_frame.append(name1)
			if name2 not in ind2_in_frame:
				try:
					dict_count2[name2] += 1
					ind2_in_frame.append(name2)
				except:
					dict_count2[name2] = 1
					ind2_in_frame.append(name2)
	return(dict_count1,dict_count2)

def write(dict_,frames,fout):
	with open(fout,'w') as f:
		for item in dict_.items():
			f.write(f'{item[0]:>16s}{item[1]/frames:8.4f}\n')

def run():
	if args.mode == 0:
		count = count_pair(args.input,args.fix_sele1,args.fix_sele2)
		count = {k: v for k, v in sorted(count.items(), key=lambda item: item[1],reverse=True)}
		write(count,args.total_frames,args.out+'.dat')
	elif args.mode == 1:
		count1, count2 = count_ind(args.input,args.fix_sele1,args.fix_sele2)
		count1 = {k: v for k, v in sorted(count1.items(), key=lambda item: item[1],reverse=True)}
		count2 = {k: v for k, v in sorted(count2.items(), key=lambda item: item[1],reverse=True)}
		write(count1,args.total_frames,args.out+'_1.dat')
		write(count2,args.total_frames,args.out+'_2.dat')

if __name__ == '__main__':
	run()