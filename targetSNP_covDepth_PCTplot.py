#Balaji Sundararaman
#03-28-22

#Usage
#python3 targetSNP_covDepth_percent_plot.py perSNP_bedCov.tsv suffix

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
#import matplotlib.image as mpimg
import sys
import pandas as pd
#import seaborn as sns
#sns.set(style="ticks")

#print(sys.argv[1])
#dfx = pd.read_csv(sys.argv[1],header=False,index_col=3,sep="\t")
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
#print(bins[:-1],"\n" cum_cov)
#print(len(cum_cov))
#print('90xy', np.interp(90,cum_cov,bins[:-1]))
#print('yx', np.interp(90,bins[:-1],cum_cov))

for xk in bins[:-1]:
    print (xk, cum_cov[bins[:-1].index(xk)])
xlabels = ['0','5','10','15','20','25','30+'] #,'35','40','45','50+']
#print(xlabels)
#print(len(xlabels))

plt.figure(figsize=(5,4))
ax1=plt.axes([0.22,0.15,0.75,0.75])

#fig, ax1 = plt.subplots(figsize=(5,4))
#pos1= ax1.get_position()
#pos2=[pos1.x0 - 0.1, pos1.y0,  pos1.width, pos1.height]

ax1.set_title(sys.argv[2], fontsize=14)
ax1.set_xlabel('Coverage', fontsize=14) # fontweight='bold', fontsize=15, family='Sans')
ax1.set_xlim(-1, 31)
ax1.set_xticks([xx for xx in range(0,31,5)])
ax1.set_xticklabels(xlabels) #, fontsize=12, family='Sans')
ax1.tick_params(axis='both', which='major', labelsize=12, length=6, width=1)
#ax1.set_xticklabels(fontsize=12, family='Sans')
ax1.plot(bins[:-1],cum_cov, '-', color='grey')
ax1.scatter(bins[:-1],cum_cov, s=10, color='k',zorder=2.5)
ax1.plot(range(-1, 31), [90 for i in range(32)], ls='--',color='grey')
#ax1.set_ylabel('% Target SNPs with >= X coverage \n(% Target SNPs at X coverage)', fontsize=16) #color='k', fontweight='bold', fontsize=15, family='Sans')
ax1.set_ylim(-5, 105)
#ax1.set_position(pos2)

ax2 = ax1.twinx()
ax2.bar(bins[:-1],cov_hist_pct, color='#bbbbbb')
#ax2.set_ylabel('% Target SNPs', color='#bbbbbb', fontweight='bold', fontsize=15, family='Sans')
ax2.set_ylim(-5, 105)
#ax2.set_yticklabels([])
ax2.axes.yaxis.set_visible(False)
#ax2.set_position(pos2)

plt.text(-0.22, 0.5, '% Target SNPs with ≥X coverage', fontsize=14, color='k', transform=ax1.transAxes, rotation='vertical',va='center',ha='center')
plt.text(-0.15, 0.5, '% Target SNPs at X coverage', fontsize=14, color='grey', transform=ax1.transAxes, rotation='vertical',va='center',ha='center')

fn = str(sys.argv[2])+"_targetSNP_covDepth.png"
plt.savefig(fn, size=(5,4), dpi=600) #, bbox_inches="tight")
print("Done plotting!")
