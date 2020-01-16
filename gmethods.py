import dcamethods as dca
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split


upath = "/home/jonah/Dropbox (ASU)/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
datap = "Projects/DCA/GunterAptamers/Selex_Data/"
datat = "Projects/DCA/ThrombinAptamers/"
datao = "Projects/DCA/ThrombinAptamers/v4/split/"
datarbm = "Projects/DCA/rbm_rna_v1/"
analysisrbm = upath + "LabFolders/Jonah_projects/RBM/"

r15p = upath + datap + 'PAL_Anna_R15_counts.txt'
r14p = upath + datap + 'PAL_Anna_R14_counts.txt'
r13p = upath + datap + 'PAL_Anna_R13_counts.txt'
r7dnap = upath + datao + '7gb.txt'
tenp = upath + datarbm + '10hidden_likelihoods.txt'

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
        if len(s) != loi:
            continue
        else:
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

a_s, s_s = read_gfile_alldata(r15p)
a_pro, s_pro = prep_data(a_s, s_s, 40, 1000, cutofftype='higher')
a_all, s_all = prep_data(a_s, s_s, 40, 0, cutofftype='higher')
b100_a, b100_s = prep_data(a_all, s_all, 40, 10, cutofftype='lower')
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
dca.write_fasta_aff(b100_s, b100_a, upath+datap+'r15_control.txt')
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
