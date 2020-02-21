import dcamethods as dca
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy import stats
import math


upath = "/home/jonah/Dropbox (ASU)/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
datap = "Projects/DCA/GunterAptamers/Selex_Data/"
datathc = 'Projects/THC/'
datat = "Projects/DCA/ThrombinAptamers/"
datao = "Projects/DCA/ThrombinAptamers/v4/split/"
datarbm = "Projects/DCA/rbm_rna_v1/"
analysisrbm = upath + "LabFolders/Jonah_projects/RBM/"

def read_likeli_new(filep):
    o = open(filep)
    ls, afs = [], []
    for line in o:
        data = line.split()
        ls.append(float(data[1]))
        afs.append(int(data[0]))
    o.close()
    return afs, ls

def likelihood_plot_rmb_wRscore(affs, likeli, title, outpath, cutoff='no'):
    a_s = list(set(affs))
    api = list(zip(affs, likeli))
    highestaff = 1
    datax, datae = [], []
    for x in a_s:
        if x > highestaff: highestaff = x
        prospects = [l for (aff, l) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    xl = np.linspace(0, highestaff, 100)
    plt.plot(xl, xl * linreg[0] + linreg[1], ':r')
    if cutoff == 'yes':
        cutoff = max([y for x, y in api if x == 1])
        plt.plot(xl, [cutoff for i in xl], ':b')
    plt.scatter(affs, likeli, color='r', s=0.5)
    plt.title(title)
    plt.ylabel('Likelihood')
    plt.xlabel('Affinity, Calc R Score: ' + str(linreg[2]))
    plt.savefig(outpath, dpi=600)
    plt.close()


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


plmA = 'plmAnalysis/'

computer_path = upath


plmp = computer_path + datap
v2p = plmp + 'v2_aligned/'
v3p = plmp + 'v3_fullalign/rbm_d/c0/'
v3dest = plmp + 'v3_fullalign/rbm_f/'
thcp = computer_path + datathc + 'v1_r15/rbm_s/'
thcd = computer_path + datathc + 'v1_r15/rbm_f/'
plma = computer_path + datap + plmA
rbmp = computer_path + datarbm
trainp = computer_path + datap + 'r15_train.txt'
testp = computer_path + datap + 'r15_test.txt'

rbmin_tmp = [str(x) + '_h_li' for x in np.arange(10, 70, 10)]
rbmin_sqtmp = [str(x) + '_h_w_li' for x in np.arange(10, 70, 10)]
nosq_all = [(x+'_t.txt', x+'_v.txt') for x in rbmin_tmp]
sq_all = [(x+'_t.txt', x+'_v.txt') for x in rbmin_sqtmp]

all_likelis = sq_all + nosq_all
all_titles = rbmin_sqtmp + rbmin_tmp
outs = [thcd + x + '.png' for x in all_titles]
print(all_titles)






# '''
cid = 0
for train, test in all_likelis:
    # tra, trl = read_likeli_new(v2p+train)
    # tea, tel = read_likeli_new(v2p+test)
    Diff_Avg_RBM(outs[cid], ['train', 'test'], thcp+train, thcp+test, title=all_titles[cid])
    # likelihood_plot_rmb_wRscore(tra, trl, all_titles[cid], v2p + train.split('.')[0] + '_plot.png', cutoff='no')
    cid += 1
    # likelihood_plot_rmb_wRscore(tea, tel, all_titles[cid], v2p + test.split('.')[0] + '_plot.png', cutoff='no')
    # cid += 1
# '''

# r15hp = plmp + 'r15_train.h'
# r15jp = plmp + 'r15_train.j'
#
# cutoffs = ['c100', 'c300', 'c500', 'c1k']
# jps = [plmp + x + '_r15_train_s.j' for x in cutoffs]
# hps = [plmp + x + '_r15_train_s.h' for x in cutoffs]
# trainps = [plmp + 'r15_train_' + x +'.txt' for x in cutoffs]
# testps = [plmp + 'r15_test_' + x +'.txt' for x in cutoffs]
# trainws = [plmp + 'r15_weights2_' + x +'.txt' for x in cutoffs]
# outs = ['c100_s.png', 'c300_nw2.png', 'c500_nw2.png', 'c1k_nw2.png']

'''
affs, seqs = dca.Fasta_Read_Aff(trainps[0])
weights = [(x/1000.)**2. for x in affs]
alldata = list(zip(affs, weights, seqs))
alldata.sort(key=lambda tup: tup[0])
alldata.reverse()
data_sep = list(zip(*alldata))
saffs, sweights, sseqs = list(data_sep[0]), list(data_sep[1]), list(data_sep[2])
wfile = plmp + 'r15c100_sorted_weights.txt'
sfile = plmp + 'r15c100_sorted_train.txt'
dca.write_fasta_aff(sseqs, saffs, sfile)
w = open(wfile, 'w')
for w8 in sweights:
    print(round(w8, 2), file=w)
w.close()
'''

# N, q = 40, 5
# for i in range(1):
#     r15j = dca.sortjmat_plmDCA(jps[i], N, q)
#     r15h = dca.sorthmat_plmDCA(hps[i], N, q)
#     jdisp = dca.FullJ_disp(r15j, N, q)
#     fig, ax = plt.subplots(1, 2)
#     dca.Fig_FullJ(ax[0], 'c100', jdisp, N, q)
#     dca.Fig_FullH(ax[1], 'c100', r15h, N, q)
#     # dca.Raw_wRscore_subplot(ax[0], r15j, r15h, trainps[i])
#     # dca.Raw_wRscore_subplot(ax[1], r15j, r15h, testps[i])
#     plt.savefig(plmp + 'c100_plmparameters.png', dpi=600)
#     plt.close()

'''
for i in range(4):
    affs, seqs = dca.Fasta_Read_Aff(trainps[i])
    weights = [(x/1000.)**2. for x in affs]
    o = open(trainws[i], 'w')
    for x in weights:
        print(round(x, 2), file=o)
    o.close()
'''

# r15jE, gvals, gdist = dca.TopJNorms_Jmatrix(r15j, N, q, 150)
# fig, ax = plt.subplots(1, 2)
# dca.Raw_wRscore_subplot(ax[0], r15j, r15h, trainp)
# dca.Raw_wRscore_subplot(ax[1], r15j, r15h, testp)
# plt.savefig(plmp + 'r15_train_vs_test.png', dpi=600)




# a_s, s_s = read_gfile_alldata(r15p)
# app, spp = prep_data(a_s, s_s, 2, 100)
# dca.write_fasta_aff(spp, app, v2p + 'r15_g100_all_lengths.txt')

# afs, seqs = read_gfile_alldata(thcp)
# app, spp = prep_data(afs, seqs, 2, 10, cutofftype='lower')
# api = list(zip(app, spp))
# sel = random.sample(api, 6000)
# af, sf = zip(*sel)
# dca.write_fasta_aff(sf, af, thcd + 'bcontrol.txt')


# dca.write_fasta_aff(spp, app, thcd + 'thc_gb.txt')

affs, seqs = dca.Fasta_Read_Aff(thcd+'thc_gb.txt')



# affs, seqs = dca.Fasta_Read_Aff(v2p +'r15_g100_all_lengths.txt')
# w1 = dca.Motif_Aligner(thcd + 'motif_finder/m1_mat.txt', 4, 11, gaps=False)
# w1.load_seqs(seqs, affs)
# w1.align_seqs()
# #
# start = w1.start_positions
# sp = [[x] for x in start]
# # print(sp)
# #
# #
# db = DBSCAN(eps=2, metric='euclidean', min_samples=1000).fit(sp)
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True
# labels = db.labels_
# #
# dis_mat = sklearn.metrics.pairwise_distances(sp, metric='euclidean')
# print(dis_mat)
# print(len(seqs), len(labels))
# #
# # # Number of clusters in labels, ignoring noise if present.
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# n_noise_ = list(labels).count(-1)
# for clust in set(labels):
#     print('Clust', clust, 'Length', list(labels).count(clust))
# #
# print('Estimated number of clusters: %d' % n_clusters_)
# print('Estimated number of noise points: %d' % n_noise_)
# #
# sqs = w1.aligned_seqs
# afs = w1.affinities
# c0, c1, c2 = [], [], []
# c0sl, c1sl, c2sl = [], [], []
# for xid, x in enumerate(db.labels_):
#     if x == 0:
#         c0sl.append(start[xid])
#         c0.append((sqs[xid], afs[xid]))
# #     elif x == 1:
# #         c1.append((sqs[xid], afs[xid]))
# #         c1sl.append(start[xid])
# #     elif x == 2:
# #         c2sl.append(start[xid])
# #         c2.append((sqs[xid], afs[xid]))
#     else:
#         continue
# # #
def adj_clust_len(sl, seqs, minl, maxl):
    ci, fi = 0, 0
    ci = min(sl)
    furthest = max(sl) - 1
    align_len = len(seqs[0])
    print(ci, furthest)
    for x in np.arange(minl, maxl+2, 1):
        trialindx = furthest + x
        print(trialindx)
        for xid, x in enumerate(seqs):
            try:
                if trialindx == align_len - 1:
                    fi = trialindx
                    break
                elif list(x)[trialindx] != '-' and xid + 1 != len(seqs) and trialindx != align_len - 1:
                    continue
                elif list(x)[trialindx] == '-' and xid+1 == len(seqs):
                    print('does this ever happen')
                    fi = trialindx
                    break
                elif list(x)[trialindx] != '-' and xid + 1 == len(seqs):
                    fi = trialindx
                    break
                else:
                    fi = align_len - 1
            except IndexError:
                fi = align_len - 1
                break
    trim_seqs = []
    print(ci, fi)
    for s in seqs:
        trim_seqs.append(''.join(list(s)[ci:fi+1]))
    return trim_seqs
# #
# # c1s, c1a = zip(*c1)
# c0s, c0a = zip(*c0)
# print(c0s[0])
# # c2s, c2a = zip(*c2)
# #
# # # dca.write_fasta_aff(c1s, c1a, v3p + 'v3_c1_all.txt')
# # # dca.write_fasta_aff(c0s, c0a, v3p + 'v3_c0_all.txt')
# # # dca.write_fasta_aff(c2s, c2a, v3p + 'v3_c2_all.txt')
# #
# nsc0 = adj_clust_len(c0sl, c0s, 38, 43)
# print(nsc0[0])
# fin = zip(c0a, nsc0)
# stratify(fin, thcd + 'thc_c0_t.txt', thcd + 'thc_c0_v.txt', weights=True, outw=thcd + 'thc_c0_w.txt')




# nsc1 = adj_clust_len(c1sl, c1s, 39, 40)
# nsc2 = adj_clust_len(c2sl, c2s, 39, 41)
#
# # dca.write_fasta_aff(nsc1, c1a, v3p + 'v3_c1_all.txt')
# dca.write_fasta_aff(nsc0, c0a, thcd + 'c0_trial.txt')
# c0e = zip(nsc0, c0a)
# c1e = zip(nsc1, c1a)
# c2e = zip(nsc2, c2a)
# cs = [c0e, c1e, c2e]
# out_t = [v3p+'c' + str(x) + '_train.txt' for x in np.arange(0,3,1)]
# out_v = [v3p+'c' + str(x) + '_test.txt' for x in np.arange(0,3,1)]
#
# stratify(c0e, out_t[0], out_v[0], weights=True, outw=v3p+'c0_train_weights.txt')
# stratify(c1e, out_t[1], out_v[1], weights=True, outw=v3p+'c1_train_weights.txt')
# stratify(c2e, out_t[2], out_v[2], weights=True, outw=v3p+'c2_train_weights.txt')
    # X_train, X_test, Y_train, Y_test = train_test_split(s, a, test_size=0.1, stratify=binned, random_state=0)
    # dca.write_fasta_aff(X_train, Y_train, v3p + out_t[xid])
    # dca.write_fasta_aff(X_test, Y_test, v3p + out_v[xid])

# dca.write_fasta_aff(nsc2, c2a, v3p + 'v3_c2_all.txt')


#
#
#
#
#

#
#
#
#
#


'''
afs, seqs = dca.Fasta_Read_Aff(v3p + 'v3_c0_all.txt')
a_s = max(list(set(afs)))
rangeEnd = math.ceil(a_s / 1000) * 1000
poss = np.arange(0, rangeEnd + 1000, 1000)
affadj = []
for x in afs:
    adj = False
    for pid, p in enumerate(poss):
        if adj:
            break
        if x < p and x > poss[pid - 1]:
            affadj.append(pid)
print(affadj)

X_train, X_test, Y_train, Y_test = train_test_split(seqs, afs, test_size=0.1, stratify=affadj, random_state=0)
dca.write_fasta_aff(X_train, Y_train, v3p+'v3_c1_train.txt')
'''



# a_pro, s_pro = prep_data(a_s, s_s, 40, 1000, cutofftype='higher')
# a_all, s_all = prep_data(a_s, s_s, 40, 0, cutofftype='higher')
# b100_a, b100_s = prep_data(a_all, s_all, 40, 10, cutofftype='lower')
# rbad_a, rbad_s = [], []
# past = []
# for i in range(3000):
#     r = random.randint(0, len(b100_s))
#     if r in past:
#         continue
#     else:
#         past.append(r)
#         if 'N' in b100_s[r]:
#             continue
        # else:
        #     rbad_a.append(b100_a[r])
        #     rbad_s.append(b100_s[r])

# X_train, X_sep, Y_train, Y_sep = train_test_split(s_pro, a_pro, test_size=0.2, random_state=0)
# X_test, X_ver, Y_test, Y_ver = train_test_split(X_sep, Y_sep, test_size=0.5, random_state=0)
# X_test += rbad_s
# Y_test += rbad_a
# xcheck = set(X_test)
# if len(xcheck) == len(X_test):
#     dca.write_fasta_aff(X_test, Y_test, upath+datap+'r15_test_c1k.txt')
# dca.write_fasta_aff(b100_s, b100_a, upath+datap+'r15_control.txt')
# dca.write_fasta_aff(X_ver, Y_ver, upath+datap+'r15_ver_c1k.txt')

# a_s, s_s = dca.Fasta_Read_Aff(r15p)

# nas, nss = prep_data(a_s, s_s, 40, 100)
# dca.write_fasta_aff(ls, a_s, upath+datao+'fam7_lh.txt')
# dca.write_fasta_aff(rs, a_s, upath+datao+'fam7_rh.txt')


# gouts = ['r15_g_' + str(i) + 'hidden.dat' for i in np.arange(10, 50, 10)]
# li = [str(x.split('.')[0]) + '_testlikeli.txt' for x in gouts]
# lic = [str(x.split('.')[0]) + '_trainlikeli.txt' for x in gouts]
# rbmout = ['Lplot_' + str(x.split('.')[0]) + '_test_wR.png' for x in gouts]
# lbmout = ['Lplot_' + str(x.split('.')[0]) + '_train_wR.png' for x in gouts]
# titlestest = [str(x.split('.')[0]) + '_test' for x in gouts]
# titlestrain = [str(x.split('.')[0]) + '_train' for x in gouts]
# train_affs, tseqs = dca.Fasta_Read_Aff(upath + datarbm + 'r15_training.txt')
# test_affs, t2seqs = dca.Fasta_Read_Aff(upath + datarbm + 'r15_test.txt')


# print(t2seqs)
# for i in range(len(li)):
#     l_test, ltseqs = read_likeli(upath + datarbm+lic[i], seqs='yes')
#     l_train = read_likeli(upath + datarbm+li[i])
#     print(len(set(t2seqs)))
#     rltseqs = [dca.rna2dna(x) for x in ltseqs]
#     ras = list(set(t2seqs) - set(rltseqs))
#     print(ras)
#     print(len(l_test), len(l_train))
#     print(len(tseqs), len(t2seqs))
#     likelihood_plot_rmb_wRscore(train_affs, l_train, titlestrain, analysisrbm+lbmout[i])
#     likelihood_plot_rmb_wRscore(test_affs, l_test, titlestest, analysisrbm+rbmout[i])

