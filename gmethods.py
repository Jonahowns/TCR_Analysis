import dcamethods as dca
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split



upath = "/home/jonah/Dropbox (ASU)/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
datap = "Projects/DCA/GuntherAptamers/Selex_Data/"
datat = "Projects/DCA/ThrombinAptamers/"
datao = "Projects/DCA/ThrombinAptamers/v4/split/"
datarbm = "Projects/DCA/rbm_rna_v1/"
analysisrbm = upath + "LabFolders/Jonah_projects/RBM/"

r15p = wpath + datap + 'PAL_Anna_R15_counts.txt'
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




a_s, s_s = read_gfile_alldata(r15p)
a_pro, s_pro = prep_data(a_s, s_s, 40, 100, cutofftype='higher')
a_all, s_all = prep_data(a_s, s_s, 40, 0, cutofftype='higher')
b100_a, b100_s = prep_data(a_all, s_all, 40, 99, cutofftype='lower')
rbad_a, rbad_s = [], []
for i in range(3000):
    r = random.randint(0, len(b100_s))
    rbad_a.append(b100_a[r])
    rbad_s.append(b100_s[r])

X_train, X_test, Y_train, Y_test = train_test_split(s_pro, a_pro, test_size=0.1, random_state=0)
X_test += rbad_s
Y_test += rbad_a
dca.write_fasta_aff(X_test, Y_test, wpath+datap+'r15_test.txt')

# nas, nss = prep_data(a_s, s_s, 40, 100)
# dca.write_fasta_aff(ls, a_s, upath+datao+'fam7_lh.txt')
# dca.write_fasta_aff(rs, a_s, upath+datao+'fam7_rh.txt')



