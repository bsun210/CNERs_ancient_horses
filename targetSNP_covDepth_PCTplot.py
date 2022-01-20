#Balaji Sundararaman
#11-24-21
#Usage
#python3 targetSNP_covDepth_percent_plot.py perSNP_bedCov.tsv suffix

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import sys
import pandas as pd

dfx = pd.read_csv(sys.argv[1],header=None,sep="\t")
#print(dfx.head())

bins = list(np.arange(0,31,1))
bins.append(1000)
#print(len(bins))

tot_targets = dfx.shape[0]
#print(sum(dfx[6]))
#print(dfx[6]/sum(dfx[6]))
cov_hist = np.histogram(dfx[6], bins=bins)
cov_hist_pct = (cov_hist[0]/tot_targets)*100.0
#print(cov_hist_pct)
#print(sum(cov_hist_pct))
#print(len(cov_hist_pct))
cum_cov = [100-sum(cov_hist_pct[:xc]) for xc in range(len(cov_hist_pct))]
#print(len(cum_cov))

for xk in bins[:-1]:
    print (xk, cum_cov[bins[:-1].index(xk)])
xlabels = ['0','5','10','15','20','25','30+'] #,'35','40','45','50+']

fig, ax1 = plt.subplots()
ax1.set_title(sys.argv[2])
ax1.set_xlabel('Coverage', fontweight='bold', fontsize=15, family='Sans')
ax1.set_xlim(-1, 31)
ax1.set_xticks([xx for xx in range(0,31,5)])
ax1.set_xticklabels(xlabels) #, fontsize=12, family='Sans')
#ax1.set_xticklabels(fontsize=12, family='Sans')
ax1.plot(bins[:-1],cum_cov, '.-', color='k')
ax1.plot(range(-1, 31), [90 for i in range(32)], ls='--',color='k')
ax1.set_ylabel('% Target SNPs with > X coverage \n[% Target SNPs]', color='k', fontweight='bold', fontsize=15, family='Sans')
ax1.set_ylim(-5, 105)

ax2 = ax1.twinx()
ax2.bar(bins[:-1],cov_hist_pct, color='#bbbbbb')
#ax2.set_ylabel('% Target SNPs', color='#bbbbbb', fontweight='bold', fontsize=15, family='Sans')
ax2.set_ylim(-5, 105)
ax2.set_yticklabels([])

#for axis in ['bottom','left']:
for axk in [ax1,ax2]:
    axk.spines['bottom'].set_linewidth(2)
    axk.spines['left'].set_linewidth(2)
    axk.spines['top'].set_visible(False)
    axk.spines['right'].set_visible(False)

plt.tick_params(labelsize=12, length=6, width=2, colors='k', pad=10)
plt.tick_params(axis='both', which='both', top=False, right=False)

fn = str(sys.argv[2])+"_targetSNP_covDepth.png"
fig.tight_layout()
plt.savefig(fn, dpi=600)
