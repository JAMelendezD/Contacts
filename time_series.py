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
	'''
	Makes a image plot of each pair as a row as a function of time
	'''
	upper_lim = np.max(matrix)
	cmap = plt.get_cmap(cmap, upper_lim+1)
	fig = plt.figure(figsize = (16,8))
	fig = plt.imshow(matrix,interpolation='none',cmap=cmap, aspect = 'auto')
	ax = plt.gca()
	cbar = plt.colorbar(fig, fraction=0.05, pad=0.05, orientation = 'horizontal', ticks=np.arange(0.45,upper_lim+1,1.01))
	plt.clim(0,upper_lim+1)
	cbar.ax.set_xticklabels(np.arange(upper_lim+1))
	ax.set_yticks(np.arange(len(names)))
	ax.set_yticklabels(names)
	ax.set_xticks(np.arange(0,frames,frames//10))
	x_labels = np.array(ax.get_xticks().tolist())*dt
	ax.set_xticklabels(np.array(x_labels,dtype=int))
	ax.xaxis.set_label_position('top')
	ax.set_xlabel('Time (ns)')
	ax.set_ylim(len(names)-0.5, -0.5)
	ax.xaxis.tick_top()
	plt.tight_layout()
	plt.savefig(output+'_mat.png', dpi = 300)
	return

def save_plot_total(series,out,sem):
	'''
	Plot the time series
	'''
	fig = plt.figure(figsize=(8,4))
	x = series[:,0]
	y = series[:,1]
	plt.plot(x,y,color='coral')
	plt.axhline(np.mean(y),0,1,color='k',ls='--')
	plt.xlim(np.min(x),np.max(x))
	plt.xlabel('Time (ns)')
	label = str(round(np.mean(y),1))+r'$\pm$'+str(round(sem,1))
	plt.text(np.max(x)-np.max(x)/9,np.mean(y)*1.02,label)
	plt.savefig(out+'.png',dpi=300,bbox_inches='tight')
	return

def save_plot_block(blocks,error,block_mean,out):
	'''
	Plots the standard error of the mean as a function of block size
	'''
	fig = plt.figure(figsize = (8,4))
	plt.plot(blocks, error, color ='k')
	plt.scatter(blocks, error, edgecolor = 'k', color = 'green', zorder=10)
	plt.xlabel('Block Size')
	plt.ylabel('SEM')
	plt.savefig(out+'block_average.png',dpi=300,bbox_inches='tight')

	fig = plt.figure(figsize = (8,4))
	plt.errorbar(blocks, block_mean, yerr=error, capsize = 2, color ='coral')
	plt.scatter(blocks, block_mean, edgecolor = 'k', color ='coral', zorder=10)
	plt.xlabel('Block Size')
	plt.ylabel('Average <X>')
	plt.savefig(out+'block_average_error.png',dpi=300,bbox_inches='tight')


def make_block_average(data):
	'''
	Block averaging procedure for a time series
	'''
	sdata = len(data)
	# If size data is odd remove the first data point. Even data points divide better in equal size blocks
	if sdata%2 == 1:
		data = data[1:] 
		sdata = sdata-1

	minblocksize = 1
	maxblocksize = int(sdata/4)

	average = np.mean(data)
	error = []
	block_mean = []
	blocks = []

	for blocksize in range(minblocksize,maxblocksize+1,1):
		s = []
		if sdata%blocksize == 0:
			nblocks = int(sdata/blocksize)
			for i in range(nblocks):
				s.append(np.mean(data[i*blocksize:(i+1)*blocksize]))
			block_mean.append(np.mean(s))
			summation = 0
			for val in s:
				summation += (val-average)**2
			sigma_sq = summation/(nblocks-1) 
			error.append(np.sqrt(sigma_sq/nblocks))
			blocks.append(blocksize)

	return(blocks,error,block_mean)

def run():
	if args.mode ==  0:
		series = count_total_all(args.total_frames,args.dt,args.input)
	elif args.mode == 1:
		series = count_total_pair(args.total_frames,args.dt,args.input)
	else:
		raise ValueError('Mode does not exists')
	np.savetxt(args.out+'.dat',series,fmt=['%10.2f','%4d'])

	blocks,error,block_mean = make_block_average(series[:,1])
	save_plot_block(blocks,error,block_mean,args.out)
	sem = error[-1]
	save_plot_total(series,args.out,sem)

	matrix, possible_dic = init_matrix(args.input_map,args.cutoff,args.total_frames)
	matrix, names = fill_matrix(args.input,matrix,possible_dic,args.fix_sele1,args.fix_sele2,args.mode)
	save_plot_matrix(matrix,'Blues',names,args.total_frames,args.dt,args.out)

if __name__ == '__main__':
	run()