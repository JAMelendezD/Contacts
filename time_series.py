import numpy as np
import argparse
import matplotlib.pyplot as plt

### Arguments ###

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Input text file')
parser.add_argument('total_frames', type=int, help='Total number of frames')
parser.add_argument('dt', type=float, help='Time interval for each frame in ns')
parser.add_argument('cutoff', type=float, help='Numeric cutoff to only make connections with values above abs(cutoff).')
parser.add_argument('--mode', default=0,required=False, type=int, help='0 to count without repeating pairs. 1 to count all. 2 to count unique based on left of file')
parser.add_argument('out', type=str, help='Name for output files.')
args = parser.parse_args()


def count_total_first(frames,dt,fname):
	'''
	Count considering only unique left
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
			pair = data[2]
			if pair not in pair_in_frame: 
				count[frame] += 1
				pair_in_frame.append(pair)
	return(np.array(list(zip(time,count)))) 

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
	Count series includes everything
	'''
	time = np.arange(0,frames)*dt
	count = np.zeros(frames)
	with open(fname, 'r') as f:
		for line in f:
			data = line.split()
			frame = int(data[0])
			count[frame] += 1
	return(np.array(list(zip(time,count)))) 

def save_plot(series,out):
	fig = plt.figure(figsize=(8,4))
	x = series[:,0]
	y = series[:,1]
	plt.plot(x,y,color='green')
	plt.axhline(np.mean(y),0,1,color='k',ls='--')
	plt.xlim(np.min(x),np.max(x))
	plt.text(np.max(x)+np.max(x)/100,np.mean(y),str(round(np.mean(y),1)))
	plt.savefig(out+'.png',dpi=300,bbox_inches='tight')

def run():
	if args.mode ==  0:
		series = count_total_all(args.total_frames,args.dt,args.input)
	elif args.mode == 1:
		series = count_total_pair(args.total_frames,args.dt,args.input)
	elif args.mode == 2:
		series = count_total_first(args.total_frames,args.dt,args.input)
	else:
		raise ValueError('Mode does not exists')
	np.savetxt(args.out+'.dat',series,fmt=['%10.2f','%4d'])
	save_plot(series,args.out)


if __name__ == '__main__':
	run()