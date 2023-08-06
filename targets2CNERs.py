#python targets2CNERs_v3.py <targets.fasta> <kmer.tsv>

import sys, csv, math, timeit
import ahocorasick

start = timeit.default_timer()

CNERs = 0
#print 'Loading SNPs fasta file...'
with open(sys.argv[1], 'r') as FA:
    FA = FA.readlines()
print('Complete making target stacks')

with open(sys.argv[2], 'r') as KMC:
        KMC = KMC.read().splitlines()
print('Complete making kmer stacks')

A = ahocorasick.Automaton()
for key in KMC:
    A.add_word(key, key)
A.make_automaton()

homo_poly = ['NNNNNN','AAAAAA','GGGGGG','CCCCCC','TTTTTT']
with open('CNERs_filter_metrics.csv', 'w+') as baits_file:
    for l in range(0, len(FA), 2):
        rsID = FA[l:l+2][0].strip()[1:]
        seq = FA[l:l+2][1].strip()
        slen = len(seq)
        gc = round(((seq.count('G')+seq.count('C'))*1.0/slen),4)
        homo = str(any(homo in seq for homo in homo_poly))
        qcs = [rsID, seq, slen, gc, homo]
        for na in [1,0.75,0.5,0.3,0.15,0.03]:
            qcs.append(round(81.5+(16.6*math.log10(na))+(0.41*gc)-(600/slen),2))

        n_found = 0
        for end_ind, found in A.iter(seq):
            n_found += 1
            if n_found == 16: break
        qcs.append(n_found)

        qcwrite = csv.writer(baits_file, delimiter='\t')
        qcwrite.writerow(qcs)
        CNERs += 1

print('Finished processing ' + str(CNERs) + ' CNERs.')
stop = timeit.default_timer()
print('Time: ' + str((stop - start)/3600) + ' hrs')
