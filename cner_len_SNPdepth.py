#Balaji Sundararaman
#01-15-2022
#usage python3 <all23771SNPs_depth.tsv> <file_name_prefix>

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys, csv, time

cner_len = ['50bp','80bp','100bp']
cov_len = {l:[] for l in cner_len}

with open(sys.argv[1], 'r') as f:
	snps = csv.reader(f, delimiter='\t')
		for snp in snps:
			c_len = snp[3].split('_')[-1]
			cov = int(snp[6])
            if cov > 0:
				cov_len[c_len].append(cov)

#for k,v in cov_len.items():
#    print(k, len(v), min(v), max(v), np.median(v))

fig, ax = plt.subplots()

if max(cov_len['80bp']) >100:
        bin_max = 100
else: bin_max = 50 #max(cov_len['80bp']) + 1

bins = list(np.arange(0,bin_max,1))
#bins.append(1000)

plt.hist([cov_len['50bp'],cov_len['80bp'],cov_len['100bp']], bins=bins, label=cner_len, color = ["#999999","#0072B2","#E69F00"], normed=True)
plt.axvline(x=np.median(cov_len['50bp']), color="#999999", lw = 1, linestyle=':')
plt.axvline(x=np.median(cov_len['80bp']), color="#0072B2", lw = 1, linestyle=':')
plt.axvline(x=np.median(cov_len['100bp']), color="#E69F00", lw = 1, linestyle=':')

s,p58 = stats.mannwhitneyu(cov_len['50bp'],cov_len['80bp'])
s,p51 = stats.mannwhitneyu(cov_len['50bp'],cov_len['100bp'])
s,p81 = stats.mannwhitneyu(cov_len['80bp'],cov_len['100bp'])

#print('{:.4e}'.format(p58))
#print('{:.4e}'.format(p51))
#print('{:.4e}'.format(p81))

plt.text(.97,.87,"80bp vs 50bp   MWW p={:.2e}".format(p58), ha="right", va="top", transform=plt.gca().transAxes, fontsize='medium')
plt.text(.97,.80,"80bp vs 100bp MWW p={:.2e}".format(p81), ha="right", va="top", transform=plt.gca().transAxes, fontsize='medium')
plt.text(.97,.73,"50bp vs 100bp MWW p={:.2e}".format(p51), ha="right", va="top", transform=plt.gca().transAxes, fontsize='medium')

plt.legend(loc='upper right', ncol=3, frameon=False, fontsize=10)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel('SNP Coverage Depth',fontsize=16)
plt.ylabel('Coverage Density',fontsize=16)
#plt.title(str(sys.argv[2]),fontsize=12)

fn = str(sys.argv[2])+"_cner_len_cov_hist.png"
fig.set_size_inches(6,4)
plt.savefig(fn, dpi=600, bbox_inches='tight')
