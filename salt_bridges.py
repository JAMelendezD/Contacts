import numpy as np
import MDAnalysis as mda
import argparse
from numba import jit
from tqdm import tqdm
from os.path import exists

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('top', type=str, help='tpr file')
parser.add_argument('traj', type=str, help='trajectory file')
parser.add_argument('first', type=int, help='First frame starts at 0')
parser.add_argument('last', type=int, help='Last frame inclusive')
parser.add_argument('d_cutoff', type=float, help='D-A distance cutoff to define the contacts')
parser.add_argument('out', type=str, help='output file')
args = parser.parse_args()

def names(selection,num_atoms):
	'''
	Function to get the residue names and store them in lists for a given selection
	'''
	sel_names = []
	for i in range(num_atoms):
		sel_names.append(selection.atoms[i].resname)
	return(sel_names)

def log(sel1_resnums,sel2_resnums,sel1_names,sel2_names,sel1_atoms,sel2_atoms,sele1,sele2,output,indeces,frame):
	'''
	Appends to a log file evertime a contact is found
	'''
	with open(output,'a') as f:
		for ele in indeces:
			i = int(ele[0])
			j = int(ele[1])
			res_num1 = sel1_resnums[i]
			res_num2 = sel2_resnums[j]
			res_name1 = sel1_names[i]
			res_name2 = sel2_names[j]
			atom_name1 = sel1_atoms[i]
			atom_name2 = sel2_atoms[j]
			f.write(f'{frame:8d}{res_name1:>5s}{res_num1:5d}{atom_name1:>5s}{sele1:>3s}{res_name2:>5s}{res_num2:5d}{atom_name2:>5s}{sele2:>3s}{ele[2]:6.2f}\n')

@jit(nopython=True)
def distance(atom1, atom2):
	'''
	Computes the euclidian distance between 2 vectors
	'''
	dx = atom2[0] - atom1[0]
	dy = atom2[1] - atom1[1]
	dz = atom2[2] - atom1[2]
	r = (dx * dx + dy * dy + dz * dz) ** 0.5
	return r

@jit(nopython=True)
def contacts(positions1,positions2,cutoff):
	'''
	Computes all contacts defined by a cutoff for a given frame
	'''
	out = []
	for i in range(len(positions1)):
		for j in range(len(positions2)):
			r = distance(positions1[i], positions2[j])
			if r <= cutoff:
				out.append([i,j,r])
	return out

def run():
	if exists(args.out+'.txt'):
		raise FileExistsError(f'File {args.out+".txt"} exists in current directory')
	if args.top.endswith('.tpr'):
		pass
	else:
		raise ValueError("Extension for topology must be .tpr")
	
	print(f'MDA version: {mda.__version__}')

	u = mda.Universe(args.top,args.traj) # Works only with .gro files since mda renumerates the .tpr
	len_traj = len(u.trajectory)
	
	print(f'The number of frames are:\t\t{len_traj:8d}')

	sele1 = 'bynum 1:6557'
	sele2 = 'bynum 13115:14918'

	basic1 = u.select_atoms(f'{sele1} and (resname ARG LYS) and (name NH* NZ)')
	acid2 = u.select_atoms(f'{sele2} and (resname ASP GLU) and (name OE* OD*)')
	basic2 = u.select_atoms(f'{sele2} and (resname ARG LYS) and (name NH* NZ)')
	acid1 = u.select_atoms(f'{sele1} and (resname ASP GLU) and (name OE* OD*)')

	num_atoms_basic1 = len(basic1)
	num_atoms_basic2 = len(basic2)
	num_atoms_acid1 = len(acid1)
	num_atoms_acid2 = len(acid2)

	print(f'The number of atoms basic 1:\t\t{num_atoms_basic1:8d}')
	print(f'The number of atoms basic 2:\t\t{num_atoms_basic2:8d}')
	print(f'The number of atoms acid 1:\t\t{num_atoms_acid1:8d}')
	print(f'The number of atoms acid 2:\t\t{num_atoms_acid2:8d}')

	basic1_resnums = list(basic1.atoms.resnums)
	basic2_resnums = list(basic2.atoms.resnums)
	acid1_resnums = list(acid1.atoms.resnums)
	acid2_resnums = list(acid2.atoms.resnums)

	print(f'The first and last resnums for basic 1:\t{basic1_resnums[0]:5d}{basic1_resnums[-1]:5d}')
	print(f'The first and last resnums for basic 2:\t{basic2_resnums[0]:5d}{basic2_resnums[-1]:5d}')
	print(f'The first and last resnums for acid 1:\t{acid1_resnums[0]:5d}{acid1_resnums[-1]:5d}')
	print(f'The first and last resnums for acid 2:\t{acid2_resnums[0]:5d}{acid2_resnums[-1]:5d}')

	basic1_names = names(basic1,num_atoms_basic1)
	basic2_names = names(basic2,num_atoms_basic2)
	acid1_names = names(acid1,num_atoms_acid1)
	acid2_names = names(acid2,num_atoms_acid2)

	basic1_atoms = list(basic1.atoms.names)
	basic2_atoms = list(basic2.atoms.names)
	acid1_atoms = list(acid1.atoms.names)
	acid2_atoms = list(acid2.atoms.names)

	for ts in tqdm(u.trajectory[args.first:args.last+1],colour='green',desc='Frames'):
		temp = contacts(basic1.positions,acid2.positions,args.d_cutoff)
		log(basic1_resnums,acid2_resnums,basic1_names,acid2_names,basic1_atoms,acid2_atoms,'A','B',args.out+'.txt',temp,int(ts.frame))
		temp = contacts(basic2.positions,acid1.positions,args.d_cutoff)
		log(basic2_resnums,acid1_resnums,basic2_names,acid1_names,basic2_atoms,acid1_atoms,'B','A',args.out+'.txt',temp,int(ts.frame))


if __name__ == '__main__':
	run()