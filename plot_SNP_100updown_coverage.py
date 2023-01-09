#Balaji
#03/28/22
#use python3
#python3 plot_SNP_100uodown_coverage.py <SNP_100updown_coverage.tsv> <prefix>

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import statistics, sys

print ('Loading SNP coverage file...')
with open(sys.argv[1], 'r') as covFile:
    covFile = covFile.readlines()

cov_pos = [[] for x in range(201)]
for sc in range(0,len(covFile),201):
    snp_cov = covFile[sc:sc+201]
    for si in range(201):
        covi = snp_cov[si].split('\t')[4].strip()
        cov_pos[si].append(int(covi))

print ('Plotting SNP coverage...')
x_pos = [xi for xi in range(201)]
y_pos2 = [statistics.median(cv) for cv in cov_pos]
y_pos1 = [np.percentile(cv, 25) for cv in cov_pos]
y_pos3 = [np.percentile(cv, 75) for cv in cov_pos]

plt.figure(figsize=(5,4))
ax=plt.axes([0.15,0.15,0.8,0.75])

ax.plot(x_pos, y_pos1, c='k', linestyle=':',linewidth=1, label='First quartile')
ax.plot(x_pos, y_pos2, c='k', linestyle='-',linewidth=1.5, label='Median')
ax.plot(x_pos, y_pos3, c='k', linestyle='--',linewidth=1, label='Third quartile')
ax.fill_between(x_pos, y_pos1, y_pos3, facecolor='grey', alpha=0.1)
ax.legend(frameon=False, loc='upper right', fontsize=10)
plt.axvline(x=100, ls='-.', lw=1, c='k')

ax.set_xlim(-5,205)
ymx = np.around(np.max(y_pos3)*1.2,0)
ax.set_ylim(-1,ymx)
ax.set_xticks(range(0,205,25))
#ax.set_yticks(range(0,21,2))
lx = [clx for clx in range(-100,101,25)]
#ly = [cly for cly in range(0,21,2)]
ax.set_xticklabels(lx)
#ax.set_yticklabels(ly)
ax.set_title(str(sys.argv[2]),fontsize=14)
ax.set_xlabel('Position relative to SNP', fontsize=14)
ax.set_ylabel('Coverage Depth', fontsize=14)
ax.tick_params(axis = 'both', which = 'major', labelsize=12, length=6, width=1, top=False, right=False, direction='out')

fn = sys.argv[2]+'_SNP_positionCov.png'
plt.savefig(fn, size=(5,4), dpi=600)#bbox_inches='tight')
print("Done plotting!")
