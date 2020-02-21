import gmethods as gm
import numpy as np
import dcamethods as dca
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import sklearn
import math
from scipy import stats


def read_likeli_new(filep):
    o = open(filep)
    ls, afs = [], []
    for line in o:
        data = line.split()
        ls.append(float(data[1]))
        afs.append(int(data[0]))
    o.close()
    return afs, ls


def adjust_plot_affinites(affs, offset=0):
    a_s = max(list(set(affs)))
    if a_s > 10000:
        rangeEnd = math.ceil(a_s / 1000) * 1000
        lr = 1000
    elif a_s > 1000:
        rangeEnd = math.ceil(a_s / 100) * 100
        lr = 100
    elif a_s > 100:
        rangeEnd = math.ceil(a_s / 10) * 10
        lr = 10
    else:
        rangeEnd = math.ceil(a_s)
        lr = 1
    plpoints = np.arange(0, rangeEnd, lr)
    affadj = []
    for x in affs:
        adj = False
        for pid, p in enumerate(plpoints):
            if adj:
                break
            if x < p and x > plpoints[pid - 1]:
                na = plpoints[pid - 1]
                affadj.append(na + offset)
                adj = True
    return affadj, lr


def Diff_Avg_RBM(outpath, labels, *infile, **kwargs):
    colors = ['g', 'm', 'c', 'b', 'y']
    title = 'Affinity vs Energy'
    for key, value in kwargs.items():
        if key == 'title':
            title = value
    adata = []
    apis = []
    data = []
    lps = []
    for iid, i in enumerate(infile):
        print(infile)
        afs, ls = read_likeli_new(i)
        linreg = stats.linregress(afs, ls)
        lps.append(linreg)
        if iid == 0:
            nafs, a_width = adjust_plot_affinites(afs)
        else:
            nafs, tmp = adjust_plot_affinites(afs, a_width/10)
        #nafs = afs
        print(iid, nafs)
        adata.append(list(set(nafs)))
        # oAff.update(affs)
        # print(oAff)
        api = list(zip(nafs, ls))
        apis += api
        data.append(api)
    for aid, ai in enumerate(data):
        avg = []
        err = []
        thaff = 0
        for aff in adata[aid]:
            if aff > thaff: thaff = aff
            yvals = np.array([y for (x, y) in ai if x == aff])
            yavg = yvals.mean()
            yerr = np.std(yvals)
            avg.append(yavg)
            err.append(yerr)
        print(len(avg))
        print(len(err))
        plt.errorbar(adata[aid], avg, err, linestyle='None', marker='^', label=labels[aid], c=colors[aid])
        xl = np.linspace(0, thaff, 100)
        plt.plot(xl, xl * lps[aid][0] + lps[aid][1], c=colors[aid], linestyle='-.')
        plt.text(xl[round(len(xl)/2)], xl[round(len(xl)/2)]*lps[aid][0] + lps[aid][1], 'RSCORE: ' + str(round(lps[aid][2],3)), c=colors[aid])
    plt.xlabel('Affinity')
    plt.ylabel('Likelihood')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0.)
    plt.suptitle(title)
    plt.savefig(outpath, dpi=600)
    plt.close()

wdir = '/home/jonah/Dropbox (ASU)/Projects/DCA/gam/'

# infiles = ['PAL_Anna_R' + str(x) + '_counts.txt' for x in np.arange(7, 16, 1)]
# mseqs, maffs = [], []
# for x in infiles:
#     aa, ss = gm.read_gfile_alldata(wdir + x)
#     ea, es = gm.prep_data(aa, ss, 0, 100, 'higher')
#     mseqs += es
#     maffs += ea
# print(len(mseqs))
# print(len(maffs))
# mas = zip(maffs, mseqs)
# fig, ax = plt.subplots(1, 1)
# dca.Fig_Distribution_w_Cutoff(ax, 'GAM Distribution', maffs, 1000)
# plt.savefig(wdir + 'datadist.png', dpi=600)

# dca.write_fasta_aff(mseqs, maffs, wdir + 'gam_motif_control.txt')


# xt, yt, xv, yv = gm.stratify(mas, wdir+'test_t.txt', wdir+'test_v.txt', True, 25000)
# fig, ax = plt.subplots(1,2)
# dca.Fig_Distribution_w_Cutoff(ax[0], 'train', yt, 0)
# dca.Fig_Distribution_w_Cutoff(ax[1], 'test', yv, 0)
# plt.savefig(wdir+'testvstrain.png')
# print(len(mseqs))

# maffs, mseqs = dca.Fasta_Read_Aff(wdir + 'gam_master.txt')
# nods = list(set(mseqs))
# noas = []
# for x in nods:
#     i = mseqs.index(x)
#     noas.append(maffs[i])
#
# a = dca.Motif_Aligner(wdir + 'm_mat.txt', 4, 10)
# a.load_seqs(nods, noas)
# a.align_seqs()
#
#
# start = a.start_positions
# sp = [[x] for x in start]
#
# db = DBSCAN(eps=2, metric='euclidean', min_samples=1000).fit(sp)
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True
# labels = db.labels_
# #
# dis_mat = sklearn.metrics.pairwise_distances(sp, metric='euclidean')
# print(dis_mat)
# print(len(mseqs), len(labels))
# #
# # # Number of clusters in labels, ignoring noise if present.
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# n_noise_ = list(labels).count(-1)
# for clust in set(labels):
#     print('Clust', clust, 'Length', list(labels).count(clust))
# #
# print('Estimated number of clusters: %d' % n_clusters_)
# print('Estimated number of noise points: %d' % n_noise_)
#
# usqs = a.unaligned_seqs
# sqs = a.aligned_seqs
# afs = a.affinities
# c0, c1, noise = [], [], []
# c0u, c1u = [], []
# c0sl, c1sl, noise = [], [], []
# for xid, x in enumerate(db.labels_):
#     if x == 0:
#         c0sl.append(start[xid])
#         c0.append((sqs[xid], afs[xid]))
#         c0u.append((usqs[xid], afs[xid]))
#     elif x == 1:
#         c1.append((sqs[xid], afs[xid]))
#         c1sl.append(start[xid])
#         c1u.append((usqs[xid], afs[xid]))
#     elif x == -1:
#         noise.append((sqs[xid], afs[xid]))
#     else:
#         continue
#
# c0s, c0a = zip(*c0)
# c1s, c1a = zip(*c1)
# c0adj = gm.adj_clust_len(c0sl, c0s, 38, 41)
# c1adj = gm.adj_clust_len(c1sl, c1s, 38, 41)
#
# c0_f = zip(c0a, c0adj)
# c1_f = zip(c1a, c1adj)
#
# c0_us, c0_ua = zip(*c0u)
# c1_us, c1_ua = zip(*c1u)
#
# dca.write_fasta_aff(usqs, afs, wdir + 'ua_t.txt')
# dca.write_fasta_aff(c0_us, c0_ua, wdir +'ua_c0.txt')
# dca.write_fasta_aff(c1_us, c1_ua, wdir + 'ua_c1.txt')

# gm.stratify(c0_f, wdir + 'c0_t.txt', wdir +'c0_v.txt', False, 10000, True, outw=wdir +'c0_w.txt')
# gm.stratify(c1_f, wdir + 'c1_t.txt', wdir +'c1_v.txt', False, 10000, True, outw=wdir +'c1_w.txt')



ydir = '/home/jonah/Dropbox (ASU)/Projects/DCA/thrombin_rbm/'
gdir = '/home/jonah/Dropbox (ASU)/Projects/DCA/gam/rbm_f/'
gdata = '/home/jonah/Dropbox (ASU)/Projects/DCA/gam/'

yrbmin_tmp = [str(x) + '_h.dat' for x in np.arange(2, 14, 4)]
yrbmwsq_tmp = [str(x.split('.')[0]) + '_w.dat' for x in yrbmin_tmp]
yrbmins = [None]*(2*len(yrbmin_tmp))
yrbmins[::2] = yrbmin_tmp
yrbmins[1::2] = yrbmwsq_tmp

grbmin_tmp = [str(x) + '_h.dat' for x in np.arange(10, 110, 50)]
grbmwsq_tmp = [str(x.split('.')[0]) + '_w.dat' for x in grbmin_tmp]
grbmins = [None]*(2*len(grbmin_tmp))
grbmins[::2] = grbmin_tmp
grbmins[1::2] = grbmwsq_tmp

gc0rbmin = ['c0_' + str(x) for x in grbmins]
# gc1rbmin = ['c1_' + str(x) for x in grbmins]

trainfiles = [[ydir+'s90.txt', ydir+'s10.txt'],[gdata + 'c0_t.txt', gdata + 'c0_v.txt']]
rbms = [yrbmins, gc0rbmin]
rdest = [ydir, gdir]

rbmins = gc0rbmin
baserbmout = [str(x.split('.')[0]) + '_li.txt' for x in rbmins]
trainrbmout = [str(x.split('.')[0]) + '_t.txt' for x in baserbmout]
testrbmout = [str(x.split('.')[0]) + '_v.txt' for x in baserbmout]
for i in range(len(gc0rbmin)):
    trainf = rdest[1] + trainrbmout[i]
    print(trainf)
    testf = rdest[1] + testrbmout[i]
    Diff_Avg_RBM(rdest[1] + str(baserbmout[i].split('.')[0]) + '.png', ['train', 'test'], trainf, testf)



# yinr = [str(x) + '_li.txt' for x in np.arange(2, 14, 2)]
# yinw = [str(x) + '_h_w.dat' for x in np.arange(2, 14, 2)]