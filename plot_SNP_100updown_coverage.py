#Balaji Sundararaman
#06/23/21
#use python3
#python3 plot_SNP_100updown_coverage.py <SNP_100updown_coverage.tsv> <prefix>

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import statistics, sys

#print ('Loading SNP coverage file...')
with open(sys.argv[1], 'r') as covFile:
    covFile = covFile.readlines()

cov_pos = [[] for x in range(201)]
for sc in range(0,len(covFile),201):
    snp_cov = covFile[sc:sc+201]
    for si in range(201):
        covi = snp_cov[si].split('\t')[4].strip()
        cov_pos[si].append(int(covi))

#print ('Plotting SNP coverage...')
x_pos = [xi for xi in range(201)]
y_pos2 = [statistics.median(cv) for cv in cov_pos]
y_pos1 = [np.percentile(cv, 25) for cv in cov_pos]
y_pos3 = [np.percentile(cv, 75) for cv in cov_pos]

cv_fig, cv_ax = plt.subplots()
cv_ax.plot(x_pos, y_pos1, c='k', linestyle=':',linewidth=1, label='First quartile')
cv_ax.plot(x_pos, y_pos2, c='k', linestyle='-',linewidth=1.5, label='Median')
cv_ax.plot(x_pos, y_pos3, c='k', linestyle='--',linewidth=1, label='Third quartile')
cv_ax.fill_between(x_pos, y_pos1, y_pos3, facecolor='grey', alpha=0.1)
cv_ax.legend(frameon=False)
plt.axvline(x=100, ls='-.', lw=1, c='k')

cv_ax.set_xlim(-5,205)
#cv_ax.set_ylim(-1,21)
cv_ax.set_xticks(range(0,205,25))
#cv_ax.set_yticks(range(0,21,2))
cv_lx = [clx for clx in range(-100,101,25)]
#cv_ly = [cly for cly in range(0,21,2)]
cv_ax.set_xticklabels(cv_lx)
#cv_ax.set_yticklabels(cv_ly)
cv_ax.set_title(str(sys.argv[2]),fontsize=15)
cv_ax.set_xlabel('Position relative to SNP', fontsize=15)
cv_ax.set_ylabel('Coverage Depth', fontsize=15)
cv_ax.tick_params(axis = 'both', which = 'major', labelsize = 12, top='off', right='off', direction='out')

fn = sys.argv[2]+'_SNP_positionCov.png'
plt.savefig(fn,figsize=(6,4), dpi=300, bbox_inches='tight')
