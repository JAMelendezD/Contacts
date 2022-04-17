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
parser.add_argument('d_cutoff', type=float, help='distance cutoff to define the contacts')
parser.add_argument('sele1', type=str, help='main selection 1')
parser.add_argument('sele2', type=str, help='main selection 2')
parser.add_argument('--mode', default=0,required=False, type=int, help='0 hydrophobic contacts. 1 all contacts')
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
	
	print(f'The number of frames are:\t\t\t{len_traj:8d}')

	sele1 = args.sele1
	sele2 = args.sele2

	if args.mode == 0:
		sel1 = u.select_atoms(f'{sele1} and name C* and (resname ALA or resname VAL or resname PHE or resname LEU or resname ILE or resname TRP or resname PRO or resname MET)')
		sel2 = u.select_atoms(f'{sele2} and name C* and (resname ALA or resname VAL or resname PHE or resname LEU or resname ILE or resname TRP or resname PRO or resname MET)')
	elif args.mode == 1:
		sel1 = u.select_atoms(f'{sele1} and not name H*')
		sel2 = u.select_atoms(f'{sele2}  and not name H*')
	else:
		raise ValueError('Mode does not exists')

	num_atoms_sel1 = len(sel1)
	num_atoms_sel2 = len(sel2)

	print(f'The number of atoms in selection 1:\t\t{num_atoms_sel1:8d}')
	print(f'The number of atoms in selection 2:\t\t{num_atoms_sel2:8d}')

	sel1_names = names(sel1,num_atoms_sel1)
	sel2_names = names(sel2,num_atoms_sel2)
	
	sel1_resnums = list(sel1.atoms.resnums)
	sel2_resnums = list(sel2.atoms.resnums)

	print(f'The first and last resnums for selection 1:\t{sel1_resnums[0]:5d}{sel1_resnums[-1]:5d}')
	print(f'The first and last resnums for selection 2:\t{sel2_resnums[0]:5d}{sel2_resnums[-1]:5d}')

	sel1_atoms = list(sel1.atoms.names)
	sel2_atoms = list(sel2.atoms.names)

	for ts in tqdm(u.trajectory[args.first:args.last+1],colour='green',desc='Frames'):
		temp = contacts(sel1.positions,sel2.positions,args.d_cutoff)
		log(sel1_resnums,sel2_resnums,sel1_names,sel2_names,sel1_atoms,sel2_atoms,'A','B',args.out+'.txt',temp,int(ts.frame))

if __name__ == '__main__':
	run()