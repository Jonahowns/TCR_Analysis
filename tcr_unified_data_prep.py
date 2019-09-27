#!/usr/bin/env python

import pandas as pd
import numpy as np
import statistics as stats
import os
from collections import OrderedDict

# DATA IMPORT
macpath = "/Users/Amber/Dropbox (ASU)/"
ubuntpath = "/home/jonah/Dropbox (ASU)/"
droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
# fullpath = macpath+droppath
# fullpath = ubuntpath + droppath
fullpath = ubuntpath + droppath
family_files = {'AAA': fullpath + 'S_7_AAA_2_norm_plus.tsv.cluster',
                'ACC': fullpath + 'S_7_ACC_2_norm_plus.tsv.cluster',
                'AGG': fullpath + 'S_7_AGG_2_norm_plus.tsv.cluster',
                'ATT': fullpath + 'S_7_ATT_2_norm_plus.tsv.cluster',
                'CAC': fullpath + 'S_7_CAC_2_norm_plus.tsv.cluster',
                'CCG': fullpath + 'S_7_CCG_2_norm_plus.tsv.cluster',
                'CTA': fullpath + 'S_7_CTA_2_norm_plus.tsv.cluster',
                'GAG': fullpath + 'S_7_GAG_2_norm_plus.tsv.cluster',
                'GCT': fullpath + 'S_7_GCT_2_norm_plus.tsv.cluster',
                'GTC': fullpath + 'S_7_GTC_2_norm_plus.tsv.cluster',
                'TAT': fullpath + 'S_7_TAT_2_norm_plus.tsv.cluster',
                'TCA': fullpath + 'S_7_TCA_2_norm_plus.tsv.cluster',
                'TGC': fullpath + 'S_7_TGC_2_norm_plus.tsv.cluster',
                'TTG': fullpath + 'S_7_TTG_2_norm_plus.tsv.cluster',
}
clustids = {1: 'AAA', 2: 'ACC', 3: 'AGG', 4: 'ATT',
            5: 'CAC', 6: 'CCG', 7: 'CTA',
            8: 'GAG', 9: 'GCT', 10: 'GTC',
            11: 'TAT', 12: 'TCA', 13: 'TGC', 14: 'TTG'}
# A Family
AAA = pd.read_csv(family_files['AAA']).dropna()
ACC = pd.read_csv(family_files['ACC']).dropna()
AGG = pd.read_csv(family_files['AGG']).dropna()
ATT = pd.read_csv(family_files['ATT']).dropna()
# C Family
CAC = pd.read_csv(family_files['CAC']).dropna()
CCG = pd.read_csv(family_files['CCG']).dropna()
CTA = pd.read_csv(family_files['CTA']).dropna()
# G Family
GAG = pd.read_csv(family_files['GAG']).dropna()
GCT = pd.read_csv(family_files['GCT']).dropna()
GTC = pd.read_csv(family_files['GTC']).dropna()
# T Family
TAT = pd.read_csv(family_files['TAT']).dropna()
TCA = pd.read_csv(family_files['TCA']).dropna()
TGC = pd.read_csv(family_files['TGC']).dropna()
TTG = pd.read_csv(family_files['TTG']).dropna()

# Cluster Selection and Seperation
Totalcomm = di.seqpercluster(AAA, ACC, AGG, ATT, CAC, CCG, CTA, GAG, GCT, GTC, TAT, TCA, TGC, TTG)
clustersofInterest = list(Totalcomm.keys())

print(clustersofInterest)


def outfilegen(clustid, types, dataset=0):
    if types == 'i':
        clustname = clustids[dataset]
        outfile = fullpath + str(clustname) + '/' + str(clustname) + str(clustid) + '.fasta'
    if types == 'f':
        outfile = fullpath + 'full' + str(clustid) + '.fasta'
    return outfile


def writefasta(seqs, file, clustid):
    y=open(file, 'w+')
    for seq, aff in seqs:
        print('>'+'clust' + str(clustid) + '-' + str(aff), file=y)
        print(seq, file=y)
    y.close()


'''
def extract_seq_full(cluster, *argv):
    # Extract Sequences for all families
    for cid, arg in enumerate(argv):
        seqonly = []
        subclust = arg[arg['cluster_id'] == int(cluster)]
        for index, row in subclust.iterrows():
            seqonly.append(row['aminoAcid'])

        writefasta(seqonly, subclustout, x)
    writefasta(fullclust, fullclustout, x)

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
'''


def length_checker(seqs, mode='check'):
    # Checks that Lengths of all Sequences are uniform
    seql = []
    for seq in seqs:
        seql.append(len(list(seq)))
    nseql = np.array(seql)
    mostcommonlen = int(stats.mode(nseql))
    lencatcher = 0
    print('Total Sequences = ' + str(len(seqs)) + '\n')
    for xid, x in enumerate(nseql):
        if x != mostcommonlen:
            if mode == 'fix':
                del seqs[xid-lencatcher]
            lencatcher += 1
    print('Not exactly ' + str(mostcommonlen) + ' bases, percentage = ' + str(lencatcher / len(nseql) * 100) + '\n')
    percentremoved = str(lencatcher / len(nseql) * 100)
    if mode == 'fix':
        print('Sequences Removed\n')
    return seqs, percentremoved


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


def clust_to_lists(clust):
    cd = clust[['aminoAcid', 'count (templates/reads)']].to_dict('list')
    seq = []
    for xid, x in enumerate(cd['aminoAcid']):
        seq.append((x, cd['count (templates/reads)'][xid]))
    return seq


class SeqHandler:
    def __init__(self, aaa, acc, agg, att, cac, ccg, cta, gag, gct, gtc, tat, tca, tgc, ttg, clusters):
        #File Handling
        self.logfile = open(fullpath+'/log.txt', 'w')
        #Datasets
        self.logger("CLUSTERS:")
        for clust in clusters:
            self.logger(str(int(clust)))
        self.logger("END OF CLUSTERS")
        self.clusters = clusters
        self.aaa = aaa
        self.acc = acc
        self.agg = agg
        self.att = att
        self.cac = cac
        self.ccg = ccg
        self.cta = cta
        self.gag = gag
        self.gct = gct
        self.gtc = gtc
        self.tat = tat
        self.tca = tca
        self.tgc = tgc
        self.ttg = ttg
        self.clust_dict = {}
        self.clust_by_fam = np.zeros((5, len(clusters)), dtype=object)
        self.clust_all = np.zeros((15, len(clusters)), dtype=object)
        self.full_clusts = np.zeros((2, len(clusters)), dtype=object)
        # Store Unique Sequences of the family
        self.uniq_seq_fam = np.zeros((4, len(clusters)), dtype=object)
        # Stores Statistics of Sequences (N of uniq seq, N of shared seq) per family per cluster
        self.uniq_seq_fam_stats = np.zeros((4, len(clusters)), dtype=object)

        for xid, clustid in enumerate(clusters):
            if xid == len(clusters)/2:
                print('Seqs Halfway Processed')
            if xid == len(clusters)-1:
                print('Seqs Processed')
            aaaSubClust = self.aaa[self.aaa['cluster_id'] == clustid]
            accSubClust = self.acc[self.acc['cluster_id'] == clustid]
            aggSubClust = self.agg[self.agg['cluster_id'] == clustid]
            attSubClust = self.att[self.att['cluster_id'] == clustid]
            cacSubClust = self.cac[self.cac['cluster_id'] == clustid]
            ccgSubClust = self.ccg[self.ccg['cluster_id'] == clustid]
            ctaSubClust = self.cta[self.cta['cluster_id'] == clustid]
            gagSubClust = self.gag[self.gag['cluster_id'] == clustid]
            gctSubClust = self.gct[self.gct['cluster_id'] == clustid]
            gtcSubClust = self.gtc[self.gtc['cluster_id'] == clustid]
            tatSubClust = self.tat[self.tat['cluster_id'] == clustid]
            tcaSubClust = self.tca[self.tca['cluster_id'] == clustid]
            tgcSubClust = self.tgc[self.tgc['cluster_id'] == clustid]
            ttgSubClust = self.ttg[self.ttg['cluster_id'] == clustid]

            aaa_seq = self.clust_to_lists(aaaSubClust)
            acc_seq = self.clust_to_lists(accSubClust)
            agg_seq = self.clust_to_lists(aggSubClust)
            att_seq = self.clust_to_lists(attSubClust)
            cac_seq = self.clust_to_lists(cacSubClust)
            ccg_seq = self.clust_to_lists(ccgSubClust)
            cta_seq = self.clust_to_lists(ctaSubClust)
            gag_seq = self.clust_to_lists(gagSubClust)
            gct_seq = self.clust_to_lists(gctSubClust)
            gtc_seq = self.clust_to_lists(gtcSubClust)
            tat_seq = self.clust_to_lists(tatSubClust)
            tca_seq = self.clust_to_lists(tcaSubClust)
            tgc_seq = self.clust_to_lists(tgcSubClust)
            ttg_seq = self.clust_to_lists(ttgSubClust)

            all_afam_seq = aaa_seq + acc_seq + agg_seq + att_seq
            tmp = OrderedDict(all_afam_seq).items()
            uniq_afam_seq = list(tmp)
            shared_a = str(len(all_afam_seq) - len(uniq_afam_seq))
            uniq_num_a = str(len(uniq_afam_seq))

            all_cfam_seq = cac_seq + ccg_seq + cta_seq
            tmp = OrderedDict(all_cfam_seq).items()
            uniq_cfam_seq = list(tmp)
            shared_c = str(len(all_cfam_seq) - len(uniq_cfam_seq))
            uniq_num_c = str(len(uniq_cfam_seq))

            all_gfam_seq = gag_seq + gct_seq + gtc_seq
            tmp = OrderedDict(all_gfam_seq).items()
            uniq_gfam_seq = list(tmp)
            shared_g = str(len(all_gfam_seq) - len(uniq_gfam_seq))
            uniq_num_g = str(len(uniq_gfam_seq))

            all_tfam_seq = tat_seq + tca_seq + tgc_seq + ttg_seq
            tmp = OrderedDict(all_tfam_seq).items()
            uniq_tfam_seq = list(tmp)
            shared_t = str(len(all_tfam_seq) - len(uniq_tfam_seq))
            uniq_num_t = str(len(uniq_tfam_seq))

            full_seq = all_afam_seq + all_cfam_seq + all_gfam_seq + all_tfam_seq
            tmp = OrderedDict(full_seq).items()
            uniq_full_seq = list(tmp)
            shared_full = str(len(full_seq) - len(uniq_full_seq))
            uniq_num_full = str(len(uniq_full_seq))

            self.logger('Cluster' + str(int(clustid)) + ' A FAM SHARED: ' + shared_a + ' UNIQ: ' + uniq_num_a,
                        'Cluster' + str(int(clustid)) + ' C FAM SHARED: ' + shared_c + ' UNIQ: ' + uniq_num_c,
                        'Cluster' + str(int(clustid)) + ' G FAM SHARED: ' + shared_g + ' UNIQ: ' + uniq_num_g,
                        'Cluster' + str(int(clustid)) + ' T FAM SHARED: ' + shared_t + ' UNIQ: ' + uniq_num_t,
                        'Cluster' + str(int(clustid)) + ' FULL SHARED: ' + shared_full + ' UNIQ: ' + uniq_num_full
                        )

            self.clust_by_fam[0, xid] = clustid
            self.clust_by_fam[1, xid] = uniq_afam_seq
            self.clust_by_fam[2, xid] = uniq_cfam_seq
            self.clust_by_fam[3, xid] = uniq_gfam_seq
            self.clust_by_fam[4, xid] = uniq_tfam_seq

            self.full_clusts[0, xid] = clustid
            self.full_clusts[1, xid] = uniq_full_seq

            self.clust_all[0, xid] = clustid
            self.clust_all[1, xid] = aaa_seq
            self.clust_all[2, xid] = acc_seq
            self.clust_all[3, xid] = agg_seq
            self.clust_all[4, xid] = att_seq
            self.clust_all[5, xid] = cac_seq
            self.clust_all[6, xid] = ccg_seq
            self.clust_all[7, xid] = cta_seq
            self.clust_all[8, xid] = gag_seq
            self.clust_all[9, xid] = gct_seq
            self.clust_all[10, xid] = gtc_seq
            self.clust_all[11, xid] = tat_seq
            self.clust_all[12, xid] = tca_seq
            self.clust_all[13, xid] = tgc_seq
            self.clust_all[14, xid] = ttg_seq

    # def clust_to_seq(self, cluster, mode='m'):
    #     # mode i -> individual cluster, mode m -> multiple clusters
    #     seqs = []
    #     if mode == 'i':
    #         for index, row in cluster.iterrows():
    #             seqs.append(row['aminoAcid'])
    #     else:
    #         for subclust in cluster:
    #             for index, row in subclust.iterrows():
    #                 seqs.append(row['aminoAcid'])
    #     return seqs

    def clust_to_lists(self, clust):
        cd = clust[['aminoAcid', 'count (templates/reads)']].to_dict('list')
        seq = []
        for xid, x in enumerate(cd['aminoAcid']):
            seq.append((x, cd['count (templates/reads)'][xid]))
        return seq

    def length_checker_by_ind(self):
        # Checks that Lengths of all Sequences are uniform
        self.logger('LENGTH INFO:')
        for xid, x in enumerate(self.clusters):
            for y in range(1, 14):
                seq, rem = length_checker(self.clust_all[y, xid])
                self.logger((str(int(x))+'-'+clustids[y+1]), rem)
        self.logger('END OF LENGTH INFO')

    def logger(self, *strings):
        for string in strings:
            print(string, file = self.logfile)

    def write_seqs(self):
        # os.mkdir(fullpath + 'SeqwAff')
        seqdir = fullpath + 'Ind/'
        for xid, x in enumerate(self.clusters):
            os.mkdir(seqdir+'Clust'+str(int(x)))
            afamfile = seqdir + 'Clust' + str(int(x)) + '/' + 'afam.fasta'
            cfamfile = seqdir + 'Clust' + str(int(x)) + '/' + 'cfam.fasta'
            gfamfile = seqdir + 'Clust' + str(int(x)) + '/' + 'gfam.fasta'
            tfamfile = seqdir + 'Clust' + str(int(x)) + '/' + 'tfam.fasta'
            fullfile = seqdir + 'Clust' + str(int(x)) + '/' + 'full.fasta'
            aaafile = seqdir + 'Clust' + str(int(x)) + '/' + 'aaa.fasta'
            accfile = seqdir + 'Clust' + str(int(x)) + '/' + 'acc.fasta'
            aggfile = seqdir + 'Clust' + str(int(x)) + '/' + 'agg.fasta'
            attfile = seqdir + 'Clust' + str(int(x)) + '/' + 'att.fasta'
            cacfile = seqdir + 'Clust' + str(int(x)) + '/' + 'cac.fasta'
            ccgfile = seqdir + 'Clust' + str(int(x)) + '/' + 'ccg.fasta'
            ctafile = seqdir + 'Clust' + str(int(x)) + '/' + 'cta.fasta'
            gagfile = seqdir + 'Clust' + str(int(x)) + '/' + 'gag.fasta'
            gctfile = seqdir + 'Clust' + str(int(x)) + '/' + 'gct.fasta'
            gtcfile = seqdir + 'Clust' + str(int(x)) + '/' + 'gtc.fasta'
            tatfile = seqdir + 'Clust' + str(int(x)) + '/' + 'tat.fasta'
            tcafile = seqdir + 'Clust' + str(int(x)) + '/' + 'tca.fasta'
            tgcfile = seqdir + 'Clust' + str(int(x)) + '/' + 'tgc.fasta'
            ttgfile = seqdir + 'Clust' + str(int(x)) + '/' + 'ttg.fasta'

            clustid = str(int(x))
            # w = 'WROTE: UNIQ SEQ CLUSTER ' + clustid
            # wa = w + '-A TO ' + afamfile
            # wc = w + '-C TO ' + cfamfile
            # wg = w + '-G TO ' + gfamfile
            # wt = w + '-T TO ' + tfamfile
            # wf = w + '-FULL TO ' + fullfile
            #
            writefasta(self.clust_by_fam[1, xid], afamfile, 'A')
            # self.logger(wa)
            writefasta(self.clust_by_fam[2, xid], cfamfile, 'C')
            # self.logger(wc)
            writefasta(self.clust_by_fam[3, xid], gfamfile, 'G')
            # self.logger(wg)
            writefasta(self.clust_by_fam[4, xid], tfamfile, 'T')
            # self.logger(wt)
            writefasta(self.full_clusts[1, xid], fullfile, 'full')
            # self.logger(wf)

            writefasta(self.clust_all[1, xid], aaafile, 'AAA')
            writefasta(self.clust_all[2, xid], accfile, 'ACC')
            writefasta(self.clust_all[3, xid], aggfile, 'AGG')
            writefasta(self.clust_all[4, xid], attfile, 'ATT')
            writefasta(self.clust_all[5, xid], cacfile, 'CAC')
            writefasta(self.clust_all[6, xid], ccgfile, 'CCG')
            writefasta(self.clust_all[7, xid], ctafile, 'CTA')
            writefasta(self.clust_all[8, xid], gagfile, 'GAG')
            writefasta(self.clust_all[9, xid], gctfile, 'GCT')
            writefasta(self.clust_all[10, xid], gtcfile, 'GTC')
            writefasta(self.clust_all[11, xid], tatfile, 'TAT')
            writefasta(self.clust_all[12, xid], tcafile, 'TCA')
            writefasta(self.clust_all[13, xid], tgcfile, 'TGC')
            writefasta(self.clust_all[14, xid], ttgfile, 'TTG')



    def __del__(self):
         print('Logfile Written to ' + fullpath)


handle = SeqHandler(AAA, ACC, AGG, ATT, CAC, CCG, CTA, GAG, GCT, GTC, TAT, TCA, TGC, TTG, clustersofInterest)
print(clustersofInterest)
# handle.length_checker_by_ind()
# handle.write_seqs()

