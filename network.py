import numpy as np
import matplotlib.pyplot as plt
import pygraphviz as pgv
import matplotlib as cm
import os
import matplotlib as mpl
import subprocess
import argparse
from matplotlib.ticker import FuncFormatter

plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'font.family': 'serif',
					"font.serif" : "Times New Roman", 
					"text.usetex": True})

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str,
					help='Name and path of the input file.')
parser.add_argument('output', type=str, 
					help='Name for the output files.')
parser.add_argument('cutoff', type=float, 
					help='Numeric cutoff to only make connections with values above abs(cutoff).')
parser.add_argument('--sort', default=0,required=False, type=int, 
					help='0 minimizes crossings. 1 sorts group1. 2 sorts both groups.')
parser.add_argument('--rev', default=0,required=False, type=int, 
					help='0 nothing. 1 inverts the top nodes when already sorted.')
parser.add_argument('--auxmark', default=0,required=False, type=int, 
					help='0 default to one marker. 1 adds a different marker.')
parser.add_argument('--auxfile', default=0,required=False, type=str, 
					help='File for second marker')
parser.add_argument('--prob', default=0,required=False, type=int, 
					help='0 used for strictly probabilities values [0,1]. 1 for values between [-1,1].')
parser.add_argument('--rowsep', default=1.5,required=False, type=float, 
					help='Spacing between the rows of nodes.')
parser.add_argument('--nodesep', default=0.05,required=False, type=float, 
					help='Spacing between the nodes.')
args = parser.parse_args()


if args.auxmark and not args.auxfile:
	raise FileNotFoundError("Provided the auxiliary marker but no file for it")

inp = args.input
cutoff = args.cutoff
out = args.output
sort_both = args.sort
reverse = args.rev
dif = args.prob
row_sep = str(args.rowsep)
node_sep = str(args.nodesep)
aux = args.auxmark
auxfile = args.auxfile

translate = {'ACE': ' ','CYS': 'C','ASP': 'D','SER': 'S','GLN': 'Q', 
     		 'LYS': 'K','ILE': 'I','PRO': 'P','THR': 'T','PHE': 'F', 
     		 'ASN': 'N','GLY': 'G','HIS': 'H','LEU': 'L','ARG': 'R', 
     		 'TRP': 'W','ALA': 'A','VAL': 'V','GLU': 'E','TYR': 'Y', 
     		 'MET': 'M','NME': ' '}

if dif == 1:
	cmap = plt.cm.bwr_r
else:
	cmap = plt.cm.plasma_r
colors = cmap(np.linspace(0, 1, 101))

def bash_command(cmd):
	'''
	Run bash commands
	'''
	p = subprocess.Popen(['/bin/bash', '-c', cmd])
	if p.wait() != 0:
		print("There was an error")
	else:
		return

def scale_lin(num,min_num,max_num,lower_bound):
	'''
	Gives a linear interpolation 
	'''
	norm = (1-lower_bound)*((num-min_num)/(max_num-min_num))+lower_bound
	return(norm)

def make_map(pairs,values,aux):
	'''
	Creates a dictionary of the appropiate connections and probabilities from 2 lists
	'''
	network = {}
	group2 = []
	n = len(pairs)
	for pair in pairs:
		sep = pair.split('-')
		res1 = translate[sep[0][0:3]]+sep[0][3:]
		res2 = translate[sep[1][0:3]]+sep[1][3:]
		if res1 not in network:
			network[res1] = []
		if res2 not in group2:
			group2.append(res2)

	group1 = list(network.keys())
	for res in group1:
		for i in range(n):
			sep = pairs[i].split('-')
			res1 = translate[sep[0][0:3]]+sep[0][3:]
			if res == res1:
				res2 = translate[sep[1][0:3]]+sep[1][3:]
				if aux == True:
					network[res].append((res2,values[i],['both','odot','odot']))
				else:
					network[res].append((res2,values[i],['None','None','None']))
	return(network,group1,group2)

def mod_base_on_aux(map,ordered1,size1,mapaux,aux1,aux2):
	'''
	Replaces elements of current map with equal elements found in second map
	'''
	for i in range(size1):
		res = ordered1[i]
		if res in aux1:
			contacts = map[res]
			for j in range(len(contacts)):
				contact = map[res][j]
				if contact[0] in aux2:
					for k in range(j,-1,-1):
						try:
							map[res][j] = mapaux[res][k]
						except:
							pass
	return(map)

if aux == 1:
	# TOFIX For some reason needs more than one element to read as array
	pairs = np.loadtxt(auxfile,usecols=0,dtype=str)
	probs = np.loadtxt(auxfile,usecols=1,dtype=float)
	indices = np.where(np.abs(probs)>=cutoff)
	prob_cut_aux = probs[indices]
	pair_cut_aux = pairs[indices]
	min_prob_aux = np.min(prob_cut_aux)
	max_prob_aux = np.max(prob_cut_aux)
	aux_net, aux1, aux2 = make_map(pair_cut_aux,prob_cut_aux,aux=True)

pairs = np.loadtxt(inp,usecols=0,dtype=str)
probs = np.loadtxt(inp,usecols=1,dtype=float)

indices = np.where(np.abs(probs)>=cutoff)
prob_cut = probs[indices]
pair_cut = pairs[indices]
n = len(prob_cut)
min_prob = np.min(prob_cut)
max_prob = np.max(prob_cut)
lim = np.max([abs(min_prob),abs(max_prob)])

if aux == 1:
	min_prob = np.min([abs(min_prob),abs(min_prob_aux)])
	max_prob = np.max([abs(max_prob),abs(max_prob_aux)])
	lim = np.max([abs(min_prob),abs(max_prob)])

print('Number of pairs under the cutoff ({}): {}'.format(cutoff,n))

network, group1, group2 = make_map(pair_cut, prob_cut,aux=False)

len_group1 = len(group1)
len_group2 = len(group2)

ordered1 = sorted(group1, key=lambda x: int(x[1:]))
ordered2 = sorted(group2, key=lambda x: int(x[1:]))

'''
Assumes probability of first map is always smaller if found in second map.
Therefore for a given cutoff for first map will always find one in second map.
If it doesnt appear in first map it will not appear at all even if present in second map.
'''
if aux == 1:
	network = mod_base_on_aux(network,ordered1,len_group1,aux_net,aux1,aux2)

C = pgv.AGraph() 
C.node_attr['style']='filled'
C.node_attr['shape']='circle'
C.node_attr['height'] =0.80
C.node_attr['fixedsize']='true'
C.node_attr['fontcolor']='#000000'
C.node_attr['fillcolor']='#FFE4B5'

if aux == 1:
	for i in range(len_group1):
		res = ordered1[i]
		contacts = network[res]
		for j in range(len(contacts)):
			contact = network[res][j]

for i in range(len_group1):
	res = ordered1[i]
	contacts = network[res]
	if sort_both != 0:
		if i < len_group1-1:
			C.add_edge(' '+ordered1[i]+' ',' '+ordered1[i+1]+' ',style='invis')
	for j in range(len(contacts)):
		contact = network[res][j]
		if dif == 0:
			C.add_edge(contact[0],' '+res+' ',
				color=cm.colors.to_hex(colors[int(scale_lin(contact[1],min_prob,max_prob,0)*100)]),
				penwidth=2.0, dir=contact[2][0],arrowtail=contact[2][1],arrowhead=contact[2][2])
		elif dif == 1:
			C.add_edge(contact[0],' '+res+' ',
				color=cm.colors.to_hex(colors[int(scale_lin(contact[1],-lim,lim,0)*100)]),
				penwidth=2.0, dir=contact[2][0],arrowtail=contact[2][1],arrowhead=contact[2][2])
	n=C.get_node(' '+res+' ')
	n.attr['fillcolor']='#B0C4DE'

if sort_both == 2:
	if reverse ==1:
		for i in range(len_group2):
			if i < len_group2-1:
				C.add_edge(ordered2[i+1],ordered2[i],style='invis')
	else:
		for i in range(len_group2):
			if i < len_group2-1:
				C.add_edge(ordered2[i],ordered2[i+1],style='invis')

if sort_both == 0:
	C.write('{}.dot'.format(out))
else:
	C.write('tmp.dot')
	com = "len=$(cat tmp.dot | wc -l);let len=$len-1;head -n $len tmp.dot > {}.dot".format(out)
	bash_command(com)
	com = "rm tmp.dot"
	bash_command(com)


if sort_both == 2:
	with open('{}.dot'.format(out),'a') as f:
		f.write("	{\n")
		f.write("		rank=same\n")
		for i in range(len_group2):
			f.write('"{}"\n'.format(ordered2[i]))
		f.write("	}\n")
		f.write("	{\n")
		f.write("		rank=same\n")
		for i in range(len_group1):
			f.write('" {} "\n'.format(ordered1[i]))
		f.write("	}\n")
		f.write("	}\n")
	f.close()
elif sort_both == 1:
	with open('{}.dot'.format(out),'a') as f:
		f.write("	{\n")
		f.write("		rank=same\n")
		for i in range(len_group1):
			f.write('" {} "\n'.format(ordered1[i]))
		f.write("	}\n")
		f.write("	}\n")
	f.close()
elif sort_both == 0:
	pass

G=pgv.AGraph('{}.dot'.format(out),ranksep = row_sep, nodesep=node_sep)
G.draw('{}.pdf'.format(out), prog = 'dot')
fig = plt.figure(figsize=(10,2))
ax = fig.add_axes([0.05, 0.80, 0.9, 0.1])
fmt = lambda x, pos: '{:.2f}'.format(x)
if dif == 1:
	cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal',
	cmap=cmap,norm=mpl.colors.Normalize(-lim,lim),format=FuncFormatter(fmt))
else:
	cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal',
	cmap=cmap,norm=mpl.colors.Normalize(cutoff, 1.0),format=FuncFormatter(fmt))
plt.savefig("{}_colorbar.pdf".format(out),bbox_inches='tight')