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



def read_gfile(filep):
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
                adj = True
            elif x < p and pid==len(poss)-1:
                affadj.append(pid-1)
                adj = True
    exceptions = []
    for x in set(affadj):
        if affadj.count(x) == 1:
            ex = affadj.index(x)
            exceptions.append(ex)
    return affadj, exceptions


def stratify(z_sanda, outtrain, outtest, returnsets=False, binwidth=1000, weights=False, outw='null'):
    a, s = zip(*z_sanda)
    sl, al = list(s), list(a)
    binned, exes = bin_affs(binwidth, al)
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
    print(len(sl), len(al), len(binned))
    X_train, X_test, Y_train, Y_test = train_test_split(sl, al, test_size=0.1, stratify=binned, random_state=0)
    X_train += stmp
    Y_train += atmp
    if weights:
        ws = [x/1000 for x in Y_train]
        o = open(outw, 'w')
        for x in ws:
            print(x, file=o)
        o.close()
    if returnsets:
        return X_train, Y_train, X_test, Y_test
    else:
        dca.write_fasta_aff(X_train, Y_train, outtrain)
        dca.write_fasta_aff(X_test, Y_test, outtest)


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