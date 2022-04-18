import MDAnalysis as mda
import argparse

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('top', type=str, help='tpr file')
parser.add_argument('traj', type=str, help='trajectory file')
args = parser.parse_args()

def run():

	u = mda.Universe(args.top,args.traj)
	len_traj = len(u.trajectory)
	print(len_traj)

if __name__ == '__main__':
	run()