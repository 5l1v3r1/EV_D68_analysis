from __future__ import division
import sys, glob
sys.path.append('/home/richard/scicore/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
sys.path.append('/scicore/home/neher/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
#import seaborn as sns

from create_allele_counts import load_allele_counts
from coverage_consensus_diversity import coverage, consensus
from minor_variant import trim_ac

def running_average(a, ws=10):
    return np.convolve(a, np.ones(ws, dtype=float)/ws, mode='same')


plt.ion()
if __name__ == '__main__':

    freqs = {}
    diversity = {}
    cov = {}
    samples = glob.glob('mapped_data/1*')

    for sample in samples:
        ac,ins = load_allele_counts(sample)
        sname = sample.rstrip('/').split('/')[-1]
        cov[sname] = coverage(ac[0][1] )
        freqs[sname] = trim_ac(ac, n_states=5)
        diversity[sname] = {ref:1 - np.sum(x**2, axis=0) for ref, x in freqs[sname].items()}

    min_cov = 1000
    snames = sorted(list(diversity.keys()))
    diversity_smooth = {}
    ref='KX675261.1'
    for sname in snames:
        good_ind = cov[sname]>min_cov
        all_pos = np.arange(good_ind.shape[0])
        tmp = []
        for i in range(3):
            rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==i)
            tmp.append(running_average(diversity[sname][ref][rf]*good_ind[rf], ws=20))
        diversity_smooth[sname] = np.array(tmp)


    coinf = ##### REDACTED SAMPLE NAMES
    average_diversity = np.ma.array([np.ma.array(x[ref], mask=cov[s]<min_cov)
                                     for s,x in diversity.items() if s not in coinf]).mean(axis=0)
    plt.figure()
    for i in range(3):
        rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==i)
        plt.plot(running_average(average_diversity[rf]), label='pos '+str(i+1))
    plt.yscale('log')
    plt.ylim([0.0003, 0.1])


    for i in range(3):
        plt.figure()
        for sname in diversity_smooth:
            if sname not in coinf:
                plt.plot(diversity_smooth[sname][i])
        plt.yscale('log')
        plt.ylim([0.001, 0.1])

