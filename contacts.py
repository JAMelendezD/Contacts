import numpy as np
import MDAnalysis as mda
import argparse
from numba import jit
from tqdm import tqdm
from os.path import exists

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('top', type=str, help='gro or tpr file')
parser.add_argument('traj', type=str, help='trajectory file')
parser.add_argument('first', type=int, help='First frame starts at 0')
parser.add_argument('last', type=int, help='Last frame inclusive')
parser.add_argument('cutoff', type=float, help='cutoff to define the contacts')
parser.add_argument('out', type=str, help='output file')
args = parser.parse_args()

def res_and_names(selection,num_atoms):
	'''
	Function to get the residue ids and the residue names and store them in lists for a given selection
	'''
	sel_resids = np.zeros(num_atoms,dtype=int)
	sel_names = []
	for i in range(num_atoms):
		sel_resids[i] = selection.atoms[i].resindex
		sel_names.append(selection.atoms[i].resname)
	return(sel_resids,sel_names)

def log(output,indeces,frame):
	'''
	Appends to a log file evertime a contact is found
	'''
	with open(output,'a') as f:
		for ele in indeces:
			i = int(ele[0])
			j = int(ele[1])
			res_num1 = sel1_resnums[i]
			res_num2 = sel2_resnums[j]
			res_id1 = sel1_resids[i]
			res_id2 = sel2_resids[j]
			res_name1 = sel1_names[i]
			res_name2 = sel2_names[j]
			atom_name1 = sel1_atoms[i]
			atom_name2 = sel2_atoms[j]
			f.write(f'{frame:8d}{res_name1:>5s}{res_num1:5d}{res_id1:5d}{atom_name1:>5s}{res_name2:>5s}{res_num2:5d}{res_id2:5d}{atom_name2:>5s}{ele[2]:5.2f}\n')

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

if __name__ == '__main__':
	if exists(args.out):
		raise FileExistsError(f'File {args.out} exists in current directory')
	
	print(f'MDA version: {mda.__version__}')

	u = mda.Universe(args.top,args.traj) # Works only with .gro files since mda renumerates the .tpr
	len_traj = len(u.trajectory)
	
	print(f'The number of frames are:\t\t\t{len_traj:8d}')

	sel1 = u.select_atoms('bynum 1:6557 and not name H*')
	sel2 = u.select_atoms('bynum 6558:13114 and not name H*')

	num_atoms_sel1 = len(sel1)
	num_atoms_sel2 = len(sel2)

	print(f'The number of atoms in selection 1:\t\t{num_atoms_sel1:8d}')
	print(f'The number of atoms in selection 2:\t\t{num_atoms_sel2:8d}')

	sel1_resnums = list(sel1.atoms.resnums)
	sel2_resnums = list(sel2.atoms.resnums)

	print(f'The first and last resnums for selection 1:\t{sel1_resnums[0]:5d}{sel1_resnums[-1]:5d}')
	print(f'The first and last resnums for selection 2:\t{sel2_resnums[0]:5d}{sel2_resnums[-1]:5d}')

	sel1_resids,sel1_names = res_and_names(sel1,num_atoms_sel1)
	sel2_resids,sel2_names = res_and_names(sel2,num_atoms_sel2)
	
	print(f'The first and last resids for selection 1:\t{sel1_resids[0]:5d}{sel1_resids[-1]:5d}')
	print(f'The first and last resids for selection 2:\t{sel2_resids[0]:5d}{sel2_resids[-1]:5d}')

	sel1_atoms = list(sel1.atoms.names)
	sel2_atoms = list(sel2.atoms.names)

	print(f'The first and last atoms for selection 1:\t{sel1_atoms[0]:>5s}{sel1_atoms[-1]:>5s}')
	print(f'The first and last atoms for selection 2:\t{sel2_atoms[0]:>5s}{sel2_atoms[-1]:>5s}')

	for ts in tqdm(u.trajectory[args.first:args.last+1],colour='green',desc='Frames'):
		temp = contacts(sel1.positions,sel2.positions,args.cutoff)
		log(args.out,temp,int(ts.frame))