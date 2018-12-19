import numpy as np
import json
from Bio import Phylo


with open('auspice/enterovirus_d68_genome_tree.json') as fh:
	T = json.load(fh)

pre_order_node_list = []
def pre_order(node):
	pre_order_node_list.append(node)
	if 'children'  in node and node["children"]:
		for c in node["children"]:
			pre_order(c)


post_order_node_list = []
def post_order(node):
	if 'children'  in node and node["children"]:
		for c in node["children"]:
			post_order(c)
	post_order_node_list.append(node)


post_order(T)
pre_order(T)

mut_count = []
terminal_mut_count = []
for n in pre_order_node_list:
	if "aa_muts" in n:
		n_aa = sum([len(x) for x in n["aa_muts"].values()])
	else:
		n_aa = 0
	n_total = len(n["muts"])
	pos = np.array([int(x[1:-1]) for x in n["muts"]])
	n_coding = np.sum((pos>699)&(pos<7263))
	mut_count.append([n_total, n_coding, n_aa] + [np.sum(pos%3==i) for i in range(3)])
	if not ('children'  in n and n["children"]):
		terminal_mut_count.append(mut_count[-1])

mut_count = np.array(mut_count)
terminal_mut_count = np.array(terminal_mut_count)

print("dN/dS",mut_count[:,2].sum()/(mut_count[:,1].sum()-mut_count[:,2].sum()))
print("terminal dN/dS",terminal_mut_count[:,2].sum()/(terminal_mut_count[:,1].sum()-terminal_mut_count[:,2].sum()))