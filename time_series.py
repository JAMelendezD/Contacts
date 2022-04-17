import numpy as np
import argparse
import matplotlib.pyplot as plt

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Input text file with raw data')
parser.add_argument('input_map', type=str, help='Input dat file with probability map')
parser.add_argument('fix_sele1', type=int, help='Number to add to the default resnum for the first selection')
parser.add_argument('fix_sele2', type=int, help='Number to add to the default resnum for the second selection')
parser.add_argument('total_frames', type=int, help='Total number of frames')
parser.add_argument('dt', type=float, help='Time interval for each frame in ns')
parser.add_argument('cutoff', type=float, help='Numeric cutoff to only make connections with values above abs(cutoff).')
parser.add_argument('--mode', default=0,required=False, type=int, help='0 to count without repeating pairs. 1 to count all.')
parser.add_argument('out', type=str, help='Name for output files.')
args = parser.parse_args()

def count_total_pair(frames,dt,fname):
	'''
	Count series considering residues only once
	'''
	time = np.arange(0,frames)*dt
	count = np.zeros(frames)
	previous_frame = 1e10
	with open(fname, 'r') as f:
		for line in f:
			data = line.split()
			frame = int(data[0])
			if frame != previous_frame:
				previous_frame = frame
				pair_in_frame = []
			pair = data[2]+'-'+data[6]
			if pair not in pair_in_frame: 
				count[frame] += 1
				pair_in_frame.append(pair)
	return(np.array(list(zip(time,count)))) 

def count_total_all(frames,dt,fname):
	'''
	Count series includes everything same result as GROMACS
	'''
	time = np.arange(0,frames)*dt
	count = np.zeros(frames)
	with open(fname, 'r') as f:
		for line in f:
			data = line.split()
			frame = int(data[0])
			count[frame] += 1
	return(np.array(list(zip(time,count)))) 


def init_matrix(fname,cutoff,frames):
	'''
	Initializes the matrix with all the pairs and the count at 0
	'''
	with open(fname,'r') as f:
		pairs = np.loadtxt(fname,usecols=0,dtype=str)
		probs = np.loadtxt(fname,usecols=1,dtype=float)

		indices = np.where(np.abs(probs)>=cutoff)
		pair_cut = pairs[indices]

	size = len(pair_cut)
	matrix = np.zeros((size,frames))

	possible_dic = {pair_cut[i]: i for i in range(size)}
	return(matrix,possible_dic)

def fill_matrix(fname,matrix,possible,fix1,fix2,mode):
	'''
	Fills the matrix by counting based on two possible modes
	'''
	if mode == 0:
		with open(fname, 'r') as f:
			for line in f:
				data = line.split()
				frame = int(data[0])
				if data[4] == 'A':
					pair = data[1]+str(int(data[2])+(fix1-1))+'-'+data[5]+str(int(data[6])+(fix2-1))
				elif data[4] == 'B':
					pair = data[5]+str(int(data[6])+(fix1-1))+'-'+data[1]+str(int(data[2])+(fix2-1))
				
				keys = possible.keys()
				if pair in keys:
					matrix[possible[pair]][frame] += 1
	elif mode == 1:
		previous_frame = 1e10
		with open(fname, 'r') as f:
			for line in f:
				data = line.split()
				frame = int(data[0])
				if frame != previous_frame:
					previous_frame = frame
					pair_in_frame = []
				if data[4] == 'A':
					pair = data[1]+str(int(data[2])+(fix1-1))+'-'+data[5]+str(int(data[6])+(fix2-1))
				elif data[4] == 'B':
					pair = data[5]+str(int(data[6])+(fix1-1))+'-'+data[1]+str(int(data[2])+(fix2-1))

				keys = possible.keys()
				if pair in keys:
					if pair not in pair_in_frame: 
						matrix[possible[pair]][frame] += 1
						pair_in_frame.append(pair)
	return(matrix,keys) 

def save_plot_matrix(matrix,cmap,names,frames,dt,output):
	upper_lim = np.max(matrix)
	cmap = plt.get_cmap(cmap, upper_lim+1)
	fig = plt.figure(figsize = (15,8))
	fig = plt.imshow(matrix,interpolation='none',cmap=cmap, aspect = 'auto')
	ax = plt.gca()
	cbar = plt.colorbar(fig, fraction=0.05, pad=0.05, orientation = 'horizontal', ticks=np.arange(0.45,upper_lim+1,1.01))
	plt.clim(0,upper_lim+1)
	cbar.ax.set_xticklabels(np.arange(upper_lim+1))
	ax.set_yticks(np.arange(len(names)))
	ax.set_yticklabels(names)
	ax.set_xticks(np.arange(0,frames,5000))
	x_labels = np.array(ax.get_xticks().tolist())*dt
	ax.set_xticklabels(np.array(x_labels,dtype=int))
	ax.xaxis.set_label_position('top')
	ax.set_xlabel('Time (ns)')
	ax.set_ylim(len(names)-0.5, -0.5)
	ax.xaxis.tick_top()
	plt.tight_layout()
	plt.savefig(output+'_mat.png', dpi = 300)
	return

def save_plot_total(series,out):
	fig = plt.figure(figsize=(8,4))
	x = series[:,0]
	y = series[:,1]
	plt.plot(x,y,color='coral')
	plt.axhline(np.mean(y),0,1,color='k',ls='--')
	plt.xlim(np.min(x),np.max(x))
	plt.text(np.max(x)+np.max(x)/100,np.mean(y),str(round(np.mean(y),1)))
	plt.savefig(out+'.png',dpi=300,bbox_inches='tight')
	return

def run():
	if args.mode ==  0:
		series = count_total_all(args.total_frames,args.dt,args.input)
	elif args.mode == 1:
		series = count_total_pair(args.total_frames,args.dt,args.input)
	else:
		raise ValueError('Mode does not exists')
	np.savetxt(args.out+'.dat',series,fmt=['%10.2f','%4d'])
	save_plot_total(series,args.out)

	matrix, possible_dic = init_matrix(args.input_map,args.cutoff,args.total_frames)
	matrix, names = fill_matrix(args.input,matrix,possible_dic,args.fix_sele1,args.fix_sele2,args.mode)
	save_plot_matrix(matrix,'Blues',names,args.total_frames,args.dt,args.out)

if __name__ == '__main__':
	run()