#!/usr/bin/env python

import pandas as pd
import os



macpath = "/Users/Amber/Dropbox (ASU)/"
droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
fullpath = macpath+droppath
'''
AAAfile = fullpath + 'S_7_AAA_2_norm_plus.tsv.cluster'
ACCfile = fullpath + 'S_7_ACC_2_norm_plus.tsv.cluster'
AGGfile = fullpath + 'S_7_AAA_2_norm_plus.tsv.cluster'


AAA = pd.read_csv(AAAfile)
ACC = pd.read_csv(ACCfile)
AGG = pd.read_csv(AGGfile)
'''


clustids = {0: 'AAA', 1: 'ACC', 2: 'AGG', 3: 'ATT',
            4: 'CAC', 5: 'CCG', 6: 'CTA',
            7: 'GAG', 8: 'GCT', 9: 'GTC',
            10: 'TAT', 11: 'TCA', 12: 'TGC', 13: 'TTG'}


def seqpercluster(*argv):
    commclust={}
    argnum=len(argv)
    shared={}
    excluded = []
    for arg in argv:
        clustid = sorted(arg['cluster_id'].unique())
        for x in clustid:
            if x not in commclust:
                commclust[x] = 1
            else:
                commclust[x] += 1
    for x, y in commclust.items():
        if y == argnum:
            shared[x] = 0
        else:
            excluded.append(x)
    for arg in argv:
        for x,y in shared.items():
            subclust = arg[arg['cluster_id'] == x]
            shared[x] += len(subclust)
    print('Excluded Clusters : ' + str(excluded))
    return shared


def outfilegen(clustid, types, dataset=0):
    if types == 'i':
        clustname = clustids[dataset]
        outfile = fullpath + str(clustname) + '/' + str(clustname) + str(clustid) + '.fasta'
    if types == 'f':
        outfile = fullpath + 'full' + str(clustid) + '.fasta'
    return outfile


def writefasta(seqs, file, clustid):
    y=open(file, 'w+')
    count = 1
    for x in seqs:
        print('>'+'seq' + str(count) + '-' + str(clustid), file=y)
        print(x, file=y)
        count += 1
    y.close()


def extractseq(clustlist, *argv):
    for x in clustlist:
        fullclustout = outfilegen(x, 'f')
        fullclust = []
        for cid, arg in enumerate(argv):
            seqonly = []
            subclust = arg[arg['cluster_id'] == x]
            for index, row in subclust.iterrows():
                seqonly.append(row['aminoAcid'])
                fullclust.append(row['aminoAcid'])
            subclustout = outfilegen(x, 'i', cid)
            writefasta(seqonly, subclustout, x)
        writefasta(fullclust, fullclustout, x)
    print("All fasta files written to " + fullpath)

def common_seq_checker(*argv):
    # Take in Lists of Seqs and compares all lists for common sequences #
    for arg1id, arg1 in enumerate(argv):
        for arg2id, arg2 in enumerate(argv):
            if arg1id == arg2id or arg2id < arg1id:
                continue
            full = arg1 + arg2
            uniq = set(full)
            print('Set ' + str(arg1id) + ' and Set ' + str(arg2id) + ' TotalSeqs: ' + str(len(full)))
            print('CommonSeq = ' + str(2*(len(full)-len(uniq))))


'''
comm = seqpercluster(AAA, ACC, AGG)
print(comm)
loi = [5, 8, 24, 25, 38, 46, 49,60]

# >10000 24,46 -- Highest amount of sequences
# ~4000  38,8
# ~800   5,49
# ~300   60,25

extractseq(loi, AAA, ACC, AGG)
'''