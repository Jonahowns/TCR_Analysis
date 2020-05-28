#!/usr/bin/env python

import pandas as pd
import numpy as np
import statistics as stats
import os
from collections import OrderedDict
import dcamethods as dca
import matplotlib.pyplot as plt
import math

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

npath = ubuntpath + 'Projects/TCR/'
nfam_files = {'c1': npath + 'Spleen_1Cluster.csv',
              'c2': npath + 'Spleen_2Cluster.csv',
              'c3': npath + 'Spleen_3Cluster.csv'}




# clustids = {1: 'AAA', 2: 'ACC', 3: 'AGG', 4: 'ATT',
#             5: 'CAC', 6: 'CCG', 7: 'CTA',
#             8: 'GAG', 9: 'GCT', 10: 'GTC',
#             11: 'TAT', 12: 'TCA', 13: 'TGC', 14: 'TTG'}
# A Family
# AAA = pd.read_csv(family_files['AAA']).dropna()
# ACC = pd.read_csv(family_files['ACC']).dropna()
# AGG = pd.read_csv(family_files['AGG']).dropna()
# ATT = pd.read_csv(family_files['ATT']).dropna()
# # C Family
# CAC = pd.read_csv(family_files['CAC']).dropna()
# CCG = pd.read_csv(family_files['CCG']).dropna()
# CTA = pd.read_csv(family_files['CTA']).dropna()
# # G Family
# GAG = pd.read_csv(family_files['GAG']).dropna()
# GCT = pd.read_csv(family_files['GCT']).dropna()
# GTC = pd.read_csv(family_files['GTC']).dropna()
# # T Family
# TAT = pd.read_csv(family_files['TAT']).dropna()
# TCA = pd.read_csv(family_files['TCA']).dropna()
# TGC = pd.read_csv(family_files['TGC']).dropna()
# TTG = pd.read_csv(family_files['TTG']).dropna()

def writefasta(seqs, file, clustid):
    y=open(file, 'w+')
    sid = 0
    for seq, aff in seqs:
        sid += 1
        print('>'+'clust' + str(clustid) + 'seq' + str(sid) + '-' + str(aff), file=y)
        print(seq, file=y)
    y.close()

def nfam_fixer(jseqs):
    nj = []
    stars_replaced = 0
    maxj = len(max(jseqs, key=len))
    for x in jseqs:
        if '*' in x:
            tmp1 = x.replace('*', '-')
            stars_replaced += 1
        else:
            tmp1 = x
        # Gaps needed to add to end
        toadd = maxj - len(tmp1)
        if toadd > 0:
            tmp2 = tmp1 + ''.join(['-' for x in range(toadd)])
        else:
            tmp2 = tmp1
        nj.append(tmp2)
    logmessages = []
    logmessages.append('Stars_Replaced = ' + str(stars_replaced))
    return nj, logmessages

def length_checker(seqs, logfile, mode='check'):
    # Checks that Lengths of all Sequences are uniform
    seql = []
    for seq in seqs:
        seql.append(len(list(seq)))
    nseql = np.array(seql)
    mostcommonlen = int(stats.mode(nseql))
    lencatcher = 0
    print('Total Sequences = ' + str(len(seqs)), file=logfile)
    for xid, x in enumerate(nseql):
        if x != mostcommonlen:
            if mode == 'fix':
                del seqs[xid-lencatcher]
            lencatcher += 1
    print('Most Common Length: ' + str(mostcommonlen) + ' Residues, percentage not most common length = ' + str(lencatcher / len(nseql) * 100), file=logfile)
    percentremoved = str(lencatcher / len(nseql) * 100)
    if mode == 'fix':
        print('Sequences Removed\n')
    return seqs, percentremoved

def import_dataset(file_dictionary):
    for xid, x in enumerate(file_dictionary.items()):
        if xid == 0:
            main = pd.read_csv(x[1]).dropna()
        else:
            main.append(pd.read_csv(x[1]).dropna())
    clusters = list(main["cluster_number"].unique())
    return main, clusters

def import_single_dataset(file):
    main = pd.read_csv(file).dropna()
    clusters = list(main["cluster_number"].unique())
    return main, clusters


def prep_data(dataset, clusters, logfile='', write_data=False, prefix=''):
    if logfile:
        l=open(logfile, 'w')
    for x in clusters:
        cdata = dataset[dataset["cluster_number"] == x]
        #Preps data as list of list tuples
        fasta_ready = cdata.values.tolist()
        if logfile:
            print('Cluster', x, 'Total Seqs:', len(fasta_ready), file=l)
        #Separate sequences, cp number, and cluster into lists
        jseqs, cpnum, clust = zip(*fasta_ready)
        # Check the length of all sequences in cluster
        length_checker(jseqs, l)
        # Check that all sequences are unique
        nonuniq = len(jseqs) - len(set(jseqs))
        print('Cluster', x, 'Repeated Sequences:', nonuniq, file=l)
        #Unique Fixes for dataset in this case for the nfam
        njseqs, logmessages = nfam_fixer(jseqs)
        for y in logmessages:
            print(y, file=l)
        # Writes the data into our favorite fasta form
        if write_data:
            outfile = npath + prefix + 'Cluster' + str(x) + '.fasta'
            writable = list(zip(njseqs, cpnum))
            writefasta(writable, outfile, x)



def npath_gen(dlabel, clabel):
    np = npath + dlabel + '/' + 'c' + clabel + '/'
    return np


def load_seqs_path(clustlist, datalabel, filelabel):
    S = []
    for i in clustlist:
        np = npath_gen(datalabel, str(i))
        fnp = np + filelabel + 'Cluster' + str(i) + '.fasta'
        seqs = dca.Fasta_Read_SeqOnly(fnp)
        S.append(seqs)
    return S


def load_seqs(clustlist, prefix='', postfix=''):
    S = []
    for i in clustlist:
        seqp = prefix + nclust_targets[i] + postfix + 'Cluster' + str(i) + '.fasta'
        seqs = dca.Fasta_Read_SeqOnly(seqp)
        S.append(seqs)
    return S

def load_HJ_from_Clust(clustid, **kwargs):
    params = False
    prefix = ''
    for key, value in kwargs.items():
        if key == 'parameters':
            params = value
        if key == 'prefix':
            prefix = value
    x = clustid
    np = npath_gen(prefix, str(clustid))
    jp = np + prefix + 'Cluster' + str(x) + '.j'
    hp = np + prefix + 'Cluster' + str(x) + '.h'
    H, n, q = dca.sorthmat_plmDCA_autoNandq(hp, gaps='yes')
    J = dca.sortjmat_plmDCA(jp, n, q, gaps='yes')
    if params:
        return H, J, n, q
    else:
        return H, J

def Clust_ALLJs(clusterlist):
    rn = len(clusterlist) % 2
    if rn:
        rownum = math.floor(len(clusterlist) / 2) + 1
    else:
        rownum = math.floor(len(clusterlist) / 2)
    fig, ax = plt.subplots(rownum, 2)
    fig.set_figheight(20)
    fig.set_figwidth(8)
    if rn:
        fig.delaxes(ax[rownum - 1, 1])

    for xid, i in enumerate(clusterlist):
        H, J, n, q = load_HJ_from_Clust(i, parameters=True)
        rp, cp = math.floor(xid / 2), xid % 2
        # dca.Fig_FullJ(ax[rp, cp], i, J, n, q, title=('Cluster ' + str(i) + ' Pairwise Parameters'), cbar=True)
        fp = nclust_paths[0] + 'clust' + str(i) + '_wl.png'
        dca.Fig_SeqLogo(fp, ax[rp, cp], i, title="Cluster " + str(i) + "SeqLogo")
        # plt.colorbar(cax=ax[rp, cp], fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(nclust_paths[0] + "ALL_SLs.png", dpi=800)


def CalcE_AllCLuster_wCoord(clusterlist, dprefix, postfix='', returnn=False, noJ=False):
    S = load_seqs_path(clusterlist, dprefix, dprefix)
    MasterCoord = []
    ns = []
    hj = []
    for xid, i in enumerate(clusterlist):
        H, J, n, q = load_HJ_from_Clust(i, parameters=True, prefix=dprefix)
        hj.append((H, J))
        ns.append(n)
        tmpclustlist = [x for x in clusterlist if x != i]
        tmpS = [x for x in S if S.index(x) != xid]
        print(tmpclustlist)
        print(tmpS)
        COI = [x for x in S if S.index(x) == xid]
        fullcoord = []
        for cid, clust in enumerate(tmpS):
            for x in clust:
                fullcoord.append((tmpclustlist[cid], dca.Calc_Energy_TCR(x, J, H, n, noJ=noJ)))
        for clust in COI:
            for x in clust:
                fullcoord.append((i, dca.Calc_Energy_TCR(x, J, H, n, noJ=noJ)))
        MasterCoord.append((i, fullcoord))
    #Master Coord holds all x and y values for all comparitive figure
    if returnn:
        return MasterCoord, ns, hj
    else:
        return MasterCoord


def CalcE_AllCLuster_SecondaryCoord(primaryclusterlist, nlist, myclusterlist, Hs, Js, dprefix, postfix='', noJ=False):
    S = load_seqs_path(myclusterlist, dprefix, dprefix)
    SecondaryCoord = []
    for xid, i in enumerate(primaryclusterlist):
        fullcoord = []
        for cid, clust in enumerate(S):
            for x in clust:
                fullcoord.append((myclusterlist[cid], dca.Calc_Energy_TCR(x, Js[xid], Hs[xid], nlist[xid], noJ=noJ)))
        SecondaryCoord.append((i, fullcoord))
    #Master Coord holds all x and y values for all comparitive figures
    return SecondaryCoord


def extract_data(tmp):
    means, stdevs = {}, {}
    xdata = list(set([x for x, y in tmp]))
    for i in xdata:
        ydata = [y for x, y in tmp if x==i]
        mean = stats.mean(ydata)
        err = stats.stdev(ydata)
        means[i] = mean
        stdevs[i] = err
    return means, stdevs


def EnergyCoord_COI_ScorePlot_AllClusters_Comparison(clusterlist1, clusterlist2, p1, p2, d1, d2, noJ=False):
    # This is the one who's cluster number and HJ will be dominant
    DOI, ns, hj = CalcE_AllCLuster_wCoord(clusterlist1, d1, postfix=p1, returnn=True, noJ=noJ)
    H, J = zip(*hj)
    # This one will compare
    MC2 = CalcE_AllCLuster_SecondaryCoord(clusterlist1, ns, clusterlist2, H, J, d2, postfix=p2, noJ=noJ)
    mclustlist = list(set(clusterlist1 + clusterlist2))
    mclustlist.sort()
    print(mclustlist)
    print('All Energies Calculated')
    rn = len(clusterlist1) % 2
    if rn:
        rownum = math.floor(len(clusterlist1)/2) + 1
    else:
        rownum = math.floor(len(clusterlist1) / 2)
    print(rownum)
    fig, ax = plt.subplots(rownum, 2)
    fig.set_figheight(20)
    fig.set_figwidth(8)
    if rn:
        fig.delaxes(ax[rownum - 1, 1])
    # Generate Figure for each Cluster
    clusterlist1.sort()
    for iid, i in enumerate(clusterlist1):
        tmp1 = [y for x, y in DOI if x == i][0]
        tmp2 = [y for x, y in MC2 if x == i][0]
        m1, s1 = extract_data(tmp1)
        m2, s2 = extract_data(tmp2)
        rp, cp = math.floor(iid / 2), iid % 2
        AvgPlot_sub_dict(ax[rp, cp], i, mclustlist, clusterlist1, clusterlist2, m1, m2, s1, s2, noJ=noJ)

    plt.tight_layout()
    if noJ:
        plt.savefig(nclust_paths[1] + nclust_targets[0] + p1 + 'v' + p2 + "ALL_COI_noJ.png", dpi=800)
    else:
        plt.savefig(nclust_paths[1] + nclust_targets[0] + p1 + 'v' + p2 + "ALL_COI.png", dpi=800)



# ScatPoints
#     scatpoints = [x for sub in specificcoord for x in sub]
#     ScatPlot(scatpoints, i)


def ScatPlot(coord, i):
    fig, ax = plt.subplots()
    ax.scatter(*zip(*coord), c=list(zip(*coord))[0], s=2)
    ax.title.set_text('Cluster ' + str(i) + ' HJ Comparison')
    plt.savefig(nclust_paths[0] + 'Cluster' + str(i) + '.png', dpi=400)
    plt.close()

def AvgPlot(xrange, means, stdevs, i):
    for xid, x in enumerate(xrange):
        plt.errorbar(x, means[xid], yerr=stdevs[xid], marker="o")
        plt.annotate("c" + str(x), (x+0.2, means[xid]-0.1))
    plt.xlabel("Cluster")
    plt.ylabel("Energy")
    plt.suptitle('Cluster ' + str(i) + ' plmDCA Comparison')
    plt.savefig(nclust_paths[0] + 'Cluster' + str(i) + '_Avg.png', dpi=400)
    plt.close()

def AvgPlot_sub(ax, xrange, means, stdevs, i):
    for xid, x in enumerate(xrange):
        ax.errorbar(x, means[xid], yerr=stdevs[xid], marker="o")
        ax.annotate("c" + str(x), (x+0.2, means[xid]-0.1), fontsize=6)
    xticks = xrange
    xticks.append(max(xrange)+1)
    ax.set_xticks(xticks)
    ax.set_yticks(list(np.arange(0, 70, 10)))
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Energy")
    ax.title.set_text('Cluster ' + str(i) + ' plmDCA Comparison')

def AvgPlot_sub_dict(ax, i, mcl, cl1, cl2, m1, m2, s1, s2, prefix='', noJ=False):

    m1s, s1s = [], []
    for x in cl1:
        d1m, s1m = m1[x], s1[x]
        m1s.append(d1m)
        s1s.append(s1m)

    m2s, s2s = [], []
    for x in cl2:
        d2m, s2m = m2[x], s2[x]
        m2s.append(d2m)
        s2s.append(s2m)

    c1 = 'aquamarine'
    c2 = 'firebrick'
    c3 = 'darkviolet'
    ax.errorbar(cl1, m1s, yerr=s1s, marker="^", color=c1, linestyle='None')
    # ax.annotate("c" + str(x), (x + 0.2, d1m + 0.5), fontsize=4)

    adjx = [x + 0.4 for x in cl2]
    ax.errorbar(adjx, m2s, yerr=s2s, marker="s", color=c2, linestyle='None')

    xticks = mcl
    xticks.append(max(mcl))
    xticks.append(0)
    offset = [x-0.3 for x in xticks]
    ax.set_xticks(xticks)
    ax.set_xticks(offset, minor=True)
    if noJ:
        ax.set_yticks(list(np.arange(0, 10, 2)))
    else:
        ax.set_yticks(list(np.arange(0, 70, 10)))
    ax.legend(['Dataset1', 'Dataset2'])
    ax.grid(axis='x', ls='--', lw=1, which='minor')
    ax.grid(axis='y', ls='--', lw=1)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Energy")
    ax.title.set_text(prefix + 'Cluster ' + str(i) + ' plmDCA Comparison')



d1, cd1 = import_single_dataset(nfam_files['c1'])
prep_data(d1, cd1, logfile=npath+'d1.log', write_data=False, prefix='/d1/d1')
print(cd1)
d2, cd2 = import_single_dataset(nfam_files['c2'])
prep_data(d2, cd2, logfile=npath+'d2.log', write_data=False, prefix='/d2/d2')
print(cd2)
d3, cd3 = import_single_dataset(nfam_files['c3'])
prep_data(d3, cd3, logfile=npath+'d3.log', write_data=False, prefix='/d3/d3')
print(cd3)




nclust_paths = {1: npath+'d1/',
                2: npath+'d2/',
                3: npath+'d3/',
                4: npath+'merged/'}


nclust_targets = {0: 'analysis/',
                1: 'c1/',
                2: 'c2/',
                3: 'c3/',
                4: 'c4/',
                5: 'c5/',
                6: 'c6/',
                7: 'c7/',
                8: 'c8/',
                9: 'c9/',
                10: 'c10/',
                11: 'c11/',
                12: 'c12/',
                13: 'c13/',
                -1: 'c-1'}




nclusts = list(np.arange(1, 14, 1))
nclusts.append(-1)




EnergyCoord_COI_ScorePlot_AllClusters_Comparison(cd1, cd2, 'd1', 'd2', 'd1', 'd2', noJ=True)
# EnergyCoord_COI_ScorePlot_AllClusters_Comparison(cd3, cd2, 'd3', 'd2', 'd3', 'd2', noJ=True)
# EnergyCoord_COI_ScorePlot_AllClusters_Comparison(cd3, cd1, 'd3', 'd1', 'd3', 'd1', noJ=True)
# EnergyCoord_COI_ScorePlot_AllClusters_Comparison(cd3, cd2, 'd3', 'd2', 'd3', 'd2', noJ=True)
# EnergyCoord_COI_ScorePlot_AllClusters_Comparison(cd, cd3, 'd2', 'd3', 'd2', 'd3', noJ=True)
# EnergyCoord_COI_ScorePlot_AllClusters_Comparison(cd2, cd1, 'd2', 'd1', 'd2', 'd1', noJ=True)



# Clust_ALLJs(nclusts)


#Visualize H and J's of each Family
# for x in nclusts:
#     jp = nclust_paths[x] + 'Cluster' + str(x) + '.j'
#     hp = nclust_paths[x] + 'Cluster' + str(x) + '.h'
#     H, n, q = dca.sorthmat_plmDCA_autoNandq(hp, gaps='yes')
#     J = dca.sortjmat_plmDCA(jp, n, q, gaps='yes')
#     fig, ax = plt.subplots(2)
#     jvis = dca.FullJ_disp(J, n, q)
#     dca.Fig_FullJ(ax[0], 'Cluster: ' + str(x), J, n, q)
#     dca.Fig_FullH(ax[1], 'Cluster: ' + str(x), H, n, q)
#     plt.savefig(nclust_paths[x] + 'Cluster'+str(x)+'plmparams.png', dpi=600)
#     plt.close()



# lfile =  + 'log.dat'
# df, clustids = import_dataset(nfam_files)
# print(clustids)
# prep_data(df, clustids, logfile=lfile, write_data=True)








# Cluster Selection and Seperation
# Totalcomm = di.seqpercluster(AAA, ACC, AGG, ATT, CAC, CCG, CTA, GAG, GCT, GTC, TAT, TCA, TGC, TTG)
# clustersofInterest = list(Totalcomm.keys())

# print(clustersofInterest)


def outfilegen(clustid, types, dataset=0):
    if types == 'i':
        clustname = clustids[dataset]
        outfile = fullpath + str(clustname) + '/' + str(clustname) + str(clustid) + '.fasta'
    if types == 'f':
        outfile = fullpath + 'full' + str(clustid) + '.fasta'
    return outfile


def writefasta(seqs, file, clustid):
    y=open(file, 'w+')
    sid = 0
    for seq, aff in seqs:
        sid += 1
        print('>'+'clust' + str(clustid) + 'seq' + str(sid) + '-' + str(aff), file=y)
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


# handle = SeqHandler(AAA, ACC, AGG, ATT, CAC, CCG, CTA, GAG, GCT, GTC, TAT, TCA, TGC, TTG, clustersofInterest)
# print(clustersofInterest)
# handle.length_checker_by_ind()
# handle.write_seqs()

