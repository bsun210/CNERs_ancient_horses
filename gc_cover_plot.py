#Balaji Sundararaman
#01-15-2022
#usage python3 <perCENR_coverage_from_picard.tsv> <file_name_prefix>

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import matplotlib.image as mpimg
import sys
import pandas as pd

print(sys.argv[1])
dfx = pd.read_csv(sys.argv[1],header=0,index_col=4,sep="\t")

#print(dfx.head())

bins = list(b/100.0 for b in range(0,101,1))
#print(len(bins))
gc_hist = np.histogram(dfx['%gc'], bins=bins)

gc_cov_weight = np.histogram(dfx['%gc'], weights=dfx['normalized_coverage'], bins=bins)
normCov_gc = np.nan_to_num(gc_cov_weight[0]/gc_hist[0],copy=False)
for i in range(len(bins[1:])):
    print(bins[i], normCov_gc[i], gc_hist[0][i])

fig1, ax1 = plt.subplots()

#ax1.set_title(sys.argv[2], fontsize=14)
ax1.set_xlabel('GC Content', fontsize=16) # fontweight='semibold', fontsize=16, family='Arial', labelpad=15)
new_y = np.ma.masked_where(normCov_gc==0, normCov_gc)
ax1.scatter(bins[1:], new_y, s=3, color='black')
ax1.set_ylabel('Normalized Coverage', fontsize=16) # fontweight='semibold', fontsize=16, family='Arial', labelpad=15)
ax1.axhline(y=1, color='black', alpha = 0.7, lw = 1, linestyle='-.')
ax1.tick_params(axis='both', which='major', labelsize=13, length=6, width=1)

ax2 = ax1.twinx()
ax2.hist(dfx['%gc'], bins=bins, color='black', alpha=0.35)
ax2.set_ylim(0, 4000)
ax2.set_yticks([yy for yy in range(0,4001,1000)])
ax2.set_yticklabels(['0','1000','2000','3000','4000'])
ax2.set_ylabel('Target GC Histogram', fontsize=16) #fontweight='bold', fontsize=16, family='Sans', labelpad=10)
ax2.tick_params(axis='both', which='major', labelsize=13, length=6, width=1)

fn = str(sys.argv[2])+"_normCov_GC.png"
fig1.set_size_inches(5,4)
plt.savefig(fn, dpi=600, bbox_inches='tight')
