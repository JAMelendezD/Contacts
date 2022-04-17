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
parser.add_argument('ang_cutoff', type=float, help='H-D-A angle cutoff to define the contacts')
parser.add_argument('sele1', type=str, help='main selection 1')
parser.add_argument('sele2', type=str, help='main selection 2')
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
			f.write(f'{frame:8d}{res_name1:>5s}{res_num1:5d}{atom_name1:>5s}{sele1:>3s}{res_name2:>5s}{res_num2:5d}{atom_name2:>5s}{sele2:>3s}{ele[2]:6.2f}{ele[3]:6.2f}\n')

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

def angle(hbonds,donors,acceptors_pos,universe):
	result = []
	for ele in hbonds:
		donor = donors[int(ele[0])]
		donor_pos = donor.position
		acceptor_pos = acceptors_pos[int(ele[1])]
		donor = donors[int(ele[0])]
		bonds = donor.bonds.indices
		possible_hs = universe.atoms[bonds[bonds != donor.index]]
		for possible_h in possible_hs:
			if possible_h.name[0] == 'H': 
				H_pos = possible_h.position
				dif1 = acceptor_pos-donor_pos
				dif2 = H_pos-donor_pos
				angle = np.rad2deg(np.arccos(np.dot(dif1,dif2)/(np.linalg.norm(dif1)*np.linalg.norm(dif2))))
				if angle <= args.ang_cutoff:
					ele.append(angle)
					result.append(ele)
	return(result)

def run():
	if exists(args.out+'.txt'):
		raise FileExistsError(f'File {args.out+".txt"} exists in current directory')
	if args.top.endswith('.tpr'):
		pass
	else:
		raise ValueError("Extension for topology must be .tpr")
	
	print(f'MDA version: {mda.__version__}')

	u = mda.Universe(args.top,args.traj)
	len_traj = len(u.trajectory)
	
	print(f'The number of frames are:\t\t\t{len_traj:8d}')

	sele1 = args.sele1
	sele2 = args.sele2

	donor1 = u.select_atoms(f'{sele1} and name O* and bonded name H*', f'{sele1} and name N* and bonded name H*',f'{sele1} and name S* and bonded name H*')
	acceptor2 = u.select_atoms(f'{sele2} and name O*', f'{sele2} and name N*')
	donor2 = u.select_atoms(f'{sele2} and name O* and bonded name H*', f'{sele2} and name N* and bonded name H*',f'{sele2} and name S* and bonded name H*')
	acceptor1 = u.select_atoms(f'{sele1} and name O*', f'{sele1} and name N*')

	num_atoms_donor1 = len(donor1)
	num_atoms_donor2 = len(donor2)
	num_atoms_acceptor1 = len(acceptor1)
	num_atoms_acceptor2 = len(acceptor2)

	print(f'The number of atoms donors 1:\t\t\t{num_atoms_donor1:8d}')
	print(f'The number of atoms donors 2:\t\t\t{num_atoms_donor2:8d}')
	print(f'The number of atoms acceptors 1:\t\t{num_atoms_acceptor1:8d}')
	print(f'The number of atoms acceptors 2:\t\t{num_atoms_acceptor2:8d}')

	donor1_resnums = list(donor1.atoms.resnums)
	donor2_resnums = list(donor2.atoms.resnums)
	acceptor1_resnums = list(acceptor1.atoms.resnums)
	acceptor2_resnums = list(acceptor2.atoms.resnums)

	print(f'The first and last resnums for donors 1:\t{donor1_resnums[0]:5d}{donor1_resnums[-1]:5d}')
	print(f'The first and last resnums for donors 2:\t{donor2_resnums[0]:5d}{donor2_resnums[-1]:5d}')
	print(f'The first and last resnums for acceptors 1:\t{acceptor1_resnums[0]:5d}{acceptor1_resnums[-1]:5d}')
	print(f'The first and last resnums for acceptors 2:\t{acceptor2_resnums[0]:5d}{acceptor2_resnums[-1]:5d}')

	donor1_names = names(donor1,num_atoms_donor1)
	donor2_names = names(donor2,num_atoms_donor2)
	acceptor1_names = names(acceptor1,num_atoms_acceptor1)
	acceptor2_names = names(acceptor2,num_atoms_acceptor2)

	donor1_atoms = list(donor1.atoms.names)
	donor2_atoms = list(donor2.atoms.names)
	acceptor1_atoms = list(acceptor1.atoms.names)
	acceptor2_atoms = list(acceptor2.atoms.names)

	for ts in tqdm(u.trajectory[args.first:args.last+1],colour='green',desc='Frames'):
		temp = contacts(donor1.positions,acceptor2.positions,args.d_cutoff)
		temp = angle(temp,donor1,acceptor2.positions,u)
		log(donor1_resnums,acceptor2_resnums,donor1_names,acceptor2_names,donor1_atoms,acceptor2_atoms,'A','B',args.out+'.txt',temp,int(ts.frame))
		temp = contacts(donor2.positions,acceptor1.positions,args.d_cutoff)
		temp = angle(temp,donor2,acceptor1.positions,u)
		log(donor2_resnums,acceptor1_resnums,donor2_names,acceptor1_names,donor2_atoms,acceptor1_atoms,'B','A',args.out+'.txt',temp,int(ts.frame))


if __name__ == '__main__':
	run()