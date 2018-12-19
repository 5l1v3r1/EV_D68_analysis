import numpy as np
import json
from Bio import Phylo
from collections import defaultdict
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
with open('auspice/enterovirus_d68_genome_tree.json') as fh:
	T = json.load(fh)


current_state = defaultdict(lambda :'ancestral')
ancestral_states = {}
state_counts = defaultdict(lambda :defaultdict(int))
mutation_count = defaultdict(list)

pre_order_node_list = []
def pre_order(node):
	pre_order_node_list.append(node)
	if 'children'  in node and node["children"]:
		for c in node["children"]:
			c["parent"] = node
			pre_order(c)


post_order_node_list = []
def post_order(node):
	if 'children'  in node and node["children"]:
		for c in node["children"]:
			post_order(c)
	post_order_node_list.append(node)


post_order(T)
pre_order(T)
L=7500
fs=16
all_pos = np.arange(7500)
mut_count = []
terminal_mut_count = []
pre_order_node_list[0]["seq"]=['Y']*L
for n in pre_order_node_list[1:]:
	n["seq"] = list(n["parent"]["seq"])
	for m in n["muts"]:
		a,pos,d = m[0], int(m[1:-1]), m[-1]
		n["seq"][pos]=d
		if d!='N':
			mutation_count[pos].append(m)
	if not ('children'  in n and n["children"]):
		for pos, x in enumerate(n["seq"]):
			state_counts[pos][n["seq"][pos]]+=1

nmuts = np.array([len(mutation_count[i]) for i in range(L)])
state_freqs = [np.array([n for n in state_counts[pos].values()]) for pos in range(L)]
state_freqs = [x/x.sum() for x in state_freqs]
entropy = np.maximum(0,-np.array([np.sum(x*np.log(x+1e-10)) for x in state_freqs]))

proteins = {'VP4':[733,940], 'VP3':[940,1684], 'VP2':[1684,2390], 'VP1':[2390,3316],
            '2A':[3316,3757], '2B':[3757, 4054], '2C':[4054,5044], '3A':[5044,5311],
            #'3B':[5311,5377],
            '3C':[5377, 5926], '3D':[5926,7297]}

loops = {'BC':[89*3,103*3], 'DE':[140*3, 153*3]}
for p,pos in proteins.items():
    proteins[p]=(np.array(pos)-733)/3
for p,pos in loops.items():
    loops[p]=(np.array(pos))/3

fig, axs = plt.subplots(3,2, sharex='col', figsize=(12,6))
cp_labels = {0:'1st', 1:'2nd', 2:'3rd'}
for ci, annos in [[0, proteins], [1,loops]]:
    for i in range(3):
        if ci==0:
            rf = (all_pos>733)&(all_pos<7297)&((all_pos-733)%3==i)
        else:
            rf = (all_pos>2390)&(all_pos<3316)&((all_pos-733)%3==i)

        axs[i, ci].plot(nmuts[rf], '-o', c='C%d'%(i+1),
                        label=cp_labels[i]+' position')
        axs[i, ci].set_ylim([0,30])
        if ci:
            axs[i, ci].text(180, 25, cp_labels[i]+' codon position', fontsize=fs)
        axs[i, ci].tick_params(labelsize=fs*0.8)
        if i==1:
            if ci==0:
                axs[i,ci].set_ylabel('# of changes', fontsize=fs)
            pi = 0
            y1,y2= (0, 30) if ci else (25,30)
            for p, pos in sorted(annos.items(), key=lambda x:x[1][0]):
                r = Rectangle((pos[0], y1),
                              pos[1]-pos[0],
                              y2-y1,
                              facecolor=[0.6+0.2*(pi%2)] * 3,
                              edgecolor='k',
                              label=p)
                pi+=1
                xt = pos[0] + 0.5 * (pos[1]-pos[0])
                yt = y2-4

                axs[i,ci].add_patch(r)
                axs[i,ci].text(xt, yt,
                    p,
                    color='k',
                    fontsize=fs*0.7,
                    ha='center')



axs[-1,0].set_xlabel('polyprotein', fontsize=fs)
axs[-1,1].set_xlabel('VP1', fontsize=fs)
plt.tight_layout()
plt.savefig('/home/richard/scicore/GROUP/data/2017_Karolinska_EV-D68/figs/between_host_diversity.pdf')
