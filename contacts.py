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
parser.add_argument('start_sele1', type=int, help='Value to shift residues of first selection to match gro')
parser.add_argument('start_sele2', type=int, help='Value to shift residues of second selection to match gro')
parser.add_argument('out', type=str, help='output file')
args = parser.parse_args()

def res_and_names(selection,num_atoms):
	'''
	Function to get the residue ids and the residue names and store them in lists for a given selection
	'''
	sel_resids = np.zeros(num_atoms,dtype=int)
	sel_names = []
	for i in range(num_atoms):
		sel_resids[i] = selection.atoms[i].resindex+1
		sel_names.append(selection.atoms[i].resname)
	return(sel_resids,sel_names)

def log(resnames1,resnames2,resids1,resids2,names1,names2,atoms1,atoms2,output,indeces,frame):
	'''
	Appends to a log file evertime a contact is found
	'''
	with open(output,'a') as f:
		for ele in indeces:
			i = int(ele[0])
			j = int(ele[1])
			res_num1 = resnames1[i]
			res_num2 = resnames2[j]
			res_id1 = resids1[i]
			res_id2 = resids2[j]
			res_name1 = names1[i]
			res_name2 = names2[j]
			atom_name1 = atoms1[i]
			atom_name2 = atoms2[j]
			f.write(f'{frame:8d}{res_name1:>5s}{res_num1:5d}{res_id1:5d}{atom_name1:>5s}{res_name2:>5s}{res_num2:5d}{res_id2:5d}{atom_name2:>5s}{ele[2]:6.2f}\n')

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

	sel1_resids,sel1_names = res_and_names(sel1,num_atoms_sel1)
	sel2_resids,sel2_names = res_and_names(sel2,num_atoms_sel2)
	
	sel1_resnums = sel1.resnums
	sel2_resnums = sel2.resnums

	print(f'The first and last resnums for selection 1:\t{sel1_resnums[0]:5d}{sel1_resnums[-1]:5d}')
	print(f'The first and last resnums for selection 2:\t{sel2_resnums[0]:5d}{sel2_resnums[-1]:5d}')

	print(f'The first and last resids for selection 1:\t{sel1_resids[0]:5d}{sel1_resids[-1]:5d}')
	print(f'The first and last resids for selection 2:\t{sel2_resids[0]:5d}{sel2_resids[-1]:5d}')

	sel1_atoms = list(sel1.atoms.names)
	sel2_atoms = list(sel2.atoms.names)

	for ts in tqdm(u.trajectory[args.first:args.last+1],colour='green',desc='Frames'):
		temp = contacts(sel1.positions,sel2.positions,args.cutoff)
		log(sel1_resnums,sel2_resnums,sel1_resids,sel2_resids,
		    sel1_names,sel2_names,sel1_atoms,sel2_atoms,
			args.out,temp,int(ts.frame))

if __name__ == '__main__':
	run()