import sys
sys.path.append('/home/richard/scicore/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns

from create_allele_counts import load_allele_counts
from coverage_consensus_diversity import coverage, consensus
from minor_variant import trim_ac
from helpers import add_panel_label

plt.ion()
if __name__ == '__main__':

    freqs = {}
    major_freqs = {}
    cov = {}
    for sample in ['JA-A', 'JA-B', 'QC_JA-A', 'QC_JA-B']:
        ac,ins = load_allele_counts('mapped_data/'+sample)
        cov[sample] = coverage(ac[0][1] )
        freqs[sample] = trim_ac(ac, n_states=5)
        major_freqs[sample] = {ref:np.max(x, axis=0) for ref, x in freqs[sample].items()}


    fs=24
    fig, axs = plt.subplots(1,2, figsize=(12,6))
    axs[0].plot([0.001, 1.0], [0.001, 1.0], c='k')
    min_cov=2000
    for s in 'AB':
        for ref in major_freqs['JA-'+s]:
            good_ind = (cov['JA-'+s]>min_cov)&(cov['QC_JA-'+s]>min_cov)
            axs[0].scatter(1.0-major_freqs['JA-'+s][ref][good_ind],
                        1.0-major_freqs['QC_JA-'+s][ref][good_ind], label='sample '+s)
    axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    axs[0].set_ylim(0.001, .5)
    axs[0].set_xlim(0.001, .5)
    axs[0].set_ylabel('minor freqencies in run 2', fontsize=fs)
    axs[0].set_xlabel('minor freqencies in run 1', fontsize=fs)
    axs[0].tick_params(labelsize=0.8*fs)
    axs[0].legend(fontsize=fs, frameon=True, edgecolor='k')
    add_panel_label(axs[0], 'A', fs, x_offset=-0.18)

    all_pos = np.arange(good_ind.shape[0])

    bins = np.logspace(-4,-2,11)
    bc=0.5*(bins[:-1]+bins[1:])
    lc = 0
    for s in 'AB':
        for ref in major_freqs['JA-'+s]:
            for ri, r in enumerate(['JA', 'QC_JA']):
                lc+=1
                good_ind = cov[r+'-'+s]>min_cov
                counts = np.bincount(np.searchsorted(bins, 1.0-major_freqs[r+'-'+s][ref][good_ind]))[1:-1]
                axs[1].plot(bc[:len(counts)], 1.0*counts/counts.sum(),
                            label=s+", run %s"%(ri+1), lw=3, c="C%d"%lc)

                # there is essentially no between codon positions in the bulk
                # rf2 = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==2)
                # counts = np.bincount(np.searchsorted(bins, 1.0-major_freqs[r+'-'+s][ref][good_ind&rf2]))[1:-1]
                # axs[1].plot(bc[:len(counts)], 1.0*counts/counts.sum(),
                #             ls='--', lw=2, c="C%d"%lc)

    axs[1].set_xscale('log')
    axs[1].set_xlabel('minor freqency',fontsize=fs)
    axs[1].set_ylabel('fraction of sites',fontsize=fs)
    axs[1].tick_params(labelsize=0.8*fs)
    axs[1].legend(fontsize=fs*0.8, frameon=True, edgecolor='k')
    add_panel_label(axs[1], 'B', fs, x_offset=-0.18)

    plt.tight_layout()
    plt.savefig('figs/variant_consistency.pdf')
