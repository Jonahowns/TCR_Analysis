import dcamethods as dca
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy import stats
import math
from sklearn.model_selection import train_test_split
from sklearn.cluster import DBSCAN
import sklearn


upath = "/home/jonah/Dropbox (ASU)/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
datap = "Projects/DCA/GunterAptamers/Selex_Data/"
v2p = upath + datap + 'v2_aligned/'
v3p = upath + datap + 'v3_fullalign/'
datat = "Projects/DCA/ThrombinAptamers/"
datao = "Projects/DCA/ThrombinAptamers/v4/split/"
datarbm = "Projects/DCA/rbm_rna_v1/"
analysisrbm = upath + "LabFolders/Jonah_projects/RBM/"

r15p = upath + datap + 'Orig_data/' + 'PAL_Anna_R15_counts.txt'
r14p = upath + datap + 'PAL_Anna_R14_counts.txt'
r13p = upath + datap + 'PAL_Anna_R13_counts.txt'
r7dnap = upath + datao + '7gb.txt'
tenp = upath + datarbm + '10hidden_likelihoods.txt'

thcd = upath + 'Projects/THC/v1_r15/'
thcp = thcd + 'THCA_R15_counts.txt'



def read_gfile_alldata(filep):
    o = open(filep, 'r')
    affs, seqs = [], []
    for line in o:
        d = line.split()
        aff, seq = int(d[0]), str(d[1])
        affs.append(aff)
        seqs.append(seq)
    o.close()
    return affs, seqs


def data_prop(affs, seqs):
    a2 = set(affs)
    a2s = sorted(a2)
    for a in a2s:
        c = 0
        for aas in affs:
            if aas == a:
                c += 1
        print('Copy Number:', a, 'Count:',  c)
    c = 0
    ltot = []
    for aid, a in enumerate(affs):
        ls = []
        if a > 100:
            c += 1
            l = len(seqs[aid])
            ls.append(l)
            ltot += ls
    print('Higher than 100:', c)
    lst = sorted(set(ltot))
    for a in lst:
        c = 0
        for sid, s in enumerate(seqs):
            if affs[sid] > 100:
                if len(s) == a:
                    c += 1
        print('Affinty > 100', 'Length:', a, 'Number of Sequences', c)
    ltotal = []
    for s in seqs:
        l = len(s)
        ltotal.append(l)
    lp = set(ltotal)
    lps = sorted(lp)
    for x in lps:
        c = 0
        for aas in seqs:
            if len(aas) == x:
                c+=1
        print('Length:', x, 'Number of Sequences', c)


def prep_data(affs, seqs, loi, afcut, cutofftype='higher'):
    faffs, fseqs = [], []
    for sid, s in enumerate(seqs):
        if affs[sid] < afcut and cutofftype == 'higher':
            continue
        if affs[sid] > afcut and cutofftype == 'lower':
            continue
        # if len(s) != loi:
        #     continue
        # else:
        faffs.append(affs[sid])
        fseqs.append(s)
    return faffs, fseqs


def split_data(seqs, split_n):
    lhs, rhs = [], []
    for aid, a in enumerate(seqs):
        n = list(a)
        lhs.append(''.join(n[:split_n]))
        rhs.append(''.join(n[split_n:]))
    return lhs, rhs


def read_likeli(filep, seqs='no'):
    o = open(filep)
    ls, seq = [], []
    for line in o:
        data = line.split()
        interest = float(data[1])
        ls.append(interest)
        if seqs == 'yes':
            seq.append(str(data[0]))
    o.close()
    if seqs == 'yes':
        return ls, seq
    else:
        return ls


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


def bin_affs(binwidth, afs):
    a_s = max(list(set(afs)))
    rangeEnd = math.ceil(a_s / binwidth) * binwidth
    poss = np.arange(0, rangeEnd + binwidth, binwidth)
    affadj = []
    for x in afs:
        adj = False
        for pid, p in enumerate(poss):
            if adj:
                break
            if x < p and x > poss[pid - 1]:
                affadj.append(pid)
    exceptions = []
    for x in set(affadj):
        if affadj.count(x) == 1:
            ex = affadj.index(x)
            exceptions.append(ex)
    return affadj, exceptions


def stratify(z_sanda, outtrain, outtest, weights=False, outw='null'):
    a, s = zip(*z_sanda)
    sl, al = list(s), list(a)
    binned, exes = bin_affs(1000, a)
    exes.sort()
    exes.reverse()
    print(exes)
    special, atmp, stmp = [], [], []
    for x in exes:
        atmp.append(al[x])
        stmp.append(sl[x])
        print(sl[x])
        del sl[x]
        del al[x]
        del binned[x]
    X_train, X_test, Y_train, Y_test = train_test_split(sl, al, test_size=0.1, stratify=binned, random_state=0)
    X_train += stmp
    Y_train += atmp
    dca.write_fasta_aff(X_train, Y_train, outtrain)
    dca.write_fasta_aff(X_test, Y_test, outtest)
    if weights:
        ws = [x/1000 for x in Y_train]
        o = open(outw, 'w')
        for x in ws:
            print(x, file=o)
        o.close()




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
w1 = dca.Motif_Aligner(thcd + 'motif_finder/m1_mat.txt', 4, 11, gaps=False)
w1.load_seqs(seqs, affs)
w1.align_seqs()
#
start = w1.start_positions
sp = [[x] for x in start]
# print(sp)
#
#
db = DBSCAN(eps=2, metric='euclidean', min_samples=1000).fit(sp)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
#
dis_mat = sklearn.metrics.pairwise_distances(sp, metric='euclidean')
print(dis_mat)
print(len(seqs), len(labels))
#
# # Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
for clust in set(labels):
    print('Clust', clust, 'Length', list(labels).count(clust))
#
print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of noise points: %d' % n_noise_)
#
sqs = w1.aligned_seqs
afs = w1.affinities
c0, c1, c2 = [], [], []
c0sl, c1sl, c2sl = [], [], []
for xid, x in enumerate(db.labels_):
    if x == 0:
        c0sl.append(start[xid])
        c0.append((sqs[xid], afs[xid]))
#     elif x == 1:
#         c1.append((sqs[xid], afs[xid]))
#         c1sl.append(start[xid])
#     elif x == 2:
#         c2sl.append(start[xid])
#         c2.append((sqs[xid], afs[xid]))
    else:
        continue
# #
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
#
# c1s, c1a = zip(*c1)
c0s, c0a = zip(*c0)
print(c0s[0])
# c2s, c2a = zip(*c2)
#
# # dca.write_fasta_aff(c1s, c1a, v3p + 'v3_c1_all.txt')
# # dca.write_fasta_aff(c0s, c0a, v3p + 'v3_c0_all.txt')
# # dca.write_fasta_aff(c2s, c2a, v3p + 'v3_c2_all.txt')
#
nsc0 = adj_clust_len(c0sl, c0s, 38, 43)
print(nsc0[0])
fin = zip(c0a, nsc0)
stratify(fin, thcd + 'thc_c0_t.txt', thcd + 'thc_c0_v.txt', weights=True, outw=thcd + 'thc_c0_w.txt')




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
