import numpy as np
import math
import copy
import matplotlib.pyplot as plt
import sys
from scipy import stats
import random
import multiprocessing as mp
from statistics import mean
from scipy.stats import gaussian_kde

################################################
## Universal Methods for Analysis of DCA Data ##
################################################
aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aad = {'-': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12,
       'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20}
rna = ['-', 'A', 'C', 'G', 'U']
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U': 4, 'T': 4}
rnan = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'U'}
dna = ['-', 'A', 'C', 'G', 'T']
nucs = ['A', 'C', 'G', 'U']
nucd = ['A', 'C', 'G', 'T']
nuc_to_id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
base_flip_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


########################################################################################################################
########################################################################################################################
# Seq Prep Methods


def sim_score(masterseq, testseq):
    ms = list(rna2dna(masterseq))
    ts = list(rna2dna(testseq))
    sc = 0
    for xid, x in enumerate(ms):
        if x == ts[xid]:
            sc += 1
    return sc/len(ms)


def prune_alignment(seqs, simt=0.99, names=[]):
    final_choice_names = []
    final_choice_seqs = []
    for sid, seq in enumerate(seqs):
        print(sid)
        append = True
        seqoi = list(seq)
        for existseq in final_choice_seqs:
            es = list(existseq)
            seq_similarity = 0.
            for i in range(len(es)):
                if seqoi[i] == es[i]:
                    seq_similarity += 1.
            seq_similarity /= len(seq)
            if seq_similarity >= simt:
                append = False
                break
        if append:
            if names:
                final_choice_names.append(names[sid])
            final_choice_seqs.append(seq.upper())
    print('INFO: reduced length of alignment from %d to %d due to sequence similarity' % (
    len(seqs), len(final_choice_seqs)), file=sys.stderr),
    if names:
        return final_choice_names, final_choice_seqs
    else:
        return final_choice_seqs


def remove_diff_len(fullseqaff):
    seql = []
    for key, value in fullseqaff.items():
        seql.append(len(key))
    nseql = np.array(seql)
    mostcommonlen = int(stats.mode(nseql)[0])
    rm = []
    for xid, dict in enumerate(fullseqaff.items()):
        key, value = dict
        if len(key) != mostcommonlen:
            rm.append(key)
    for x in rm:
        fullseqaff.pop(x)
    return fullseqaff


def prep_full_fam_seqs(famid, sim):
    csvfile = '/home/jonah/Downloads/2HX_' + str(famid) + 'th_new.csv'
    bb1 = '/home/jonah/Downloads/' + str(famid) + 'bb1.txt'
    bb2 = '/home/jonah/Downloads/' + str(famid) + 'bb2.txt'
    bb3 = '/home/jonah/Downloads/' + str(famid) + 'bb3.txt'
    bb4 = '/home/jonah/Downloads/' + str(famid) + 'bb4.txt'
    bb5 = '/home/jonah/Downloads/' + str(famid) + 'bb5.txt'
    gb = '/home/jonah/Downloads/' + str(famid) + 'gb.txt'

    cutoff = 10.0
    o = open(csvfile, 'r')
    next(o)
    fullseqaff = {}
    print('better get here at least')
    enough = False
    count = 0
    while enough is False:
        count += 1
        line = next(o)
        seq, aff = line.split(';')
        rs = []
        for sid, x in enumerate(list(seq)):
            if x == 'T' or x == 't':
                rs.append('U')
            else:
                rs.append(x.upper())
        fullseqaff[''.join(rs)] = aff.rstrip()
        if count >= 40000:
            enough = True
    o.close()
    print('file read')

    fullseqaff = remove_diff_len(fullseqaff)
    print('removed diff lens')
    check = False
    while check is False:
        # abovecutoff = 0
        abovecutoff = sum(1 for x in fullseqaff.values() if float(x) >= cutoff)
        if abovecutoff > 300:
            check = True
        else:
            cutoff -= 1
            check = False

    print('Sequences w/ affinity >= ' + str(cutoff) + ' -- ' + str(abovecutoff))
    gbtitles = []
    gbseqs = []
    bbtitles = []
    bbseqs = []
    for xid, x in enumerate(fullseqaff.items()):
        key, value = x
        if float(value) >= cutoff:
            gbtitles.append(value)
            gbseqs.append(key)
        else:
            bbtitles.append(value)
            bbseqs.append(key)

    g = open(gb, 'w')
    for xid, x in enumerate(gbseqs):
        print('>seq' + str(xid) + '-' + str(gbtitles[xid]), file=g)
        print(x, file=g)
    g.close()

    print(len(bbtitles))
    print(len(bbseqs))

    chosenbbtitles = []
    chosenbbseqs = []

    for x in range(25000):
        z = random.randint(0, len(bbseqs) - 1)
        print(z)
        chosenbbtitles.append(bbtitles[z])
        chosenbbseqs.append(bbseqs[z])

    fbbt, fbbs = prune_alignment(chosenbbtitles, chosenbbseqs, sim)

    totalseqs = len(fbbs)
    seqsperfile = math.floor(totalseqs / 5)
    f1 = open(bb1, 'w')
    f2 = open(bb2, 'w')
    f3 = open(bb3, 'w')
    f4 = open(bb4, 'w')
    f5 = open(bb5, 'w')
    for x in range(len(fbbs)):
        if x < seqsperfile:
            print('>seq' + str(x) + '-' + str(fbbt[x]), file=f1)
            print(fbbs[x], file=f1)
        elif x > seqsperfile and x < 2 * seqsperfile:
            print('>seq' + str(x) + '-' + str(fbbt[x]), file=f2)
            print(fbbs[x], file=f2)
        elif x > 2 * seqsperfile and x < 3 * seqsperfile:
            print('>seq' + str(x) + '-' + str(fbbt[x]), file=f3)
            print(fbbs[x], file=f3)
        elif x > 3 * seqsperfile and x < 4 * seqsperfile:
            print('>seq' + str(x) + '-' + str(fbbt[x]), file=f4)
            print(fbbs[x], file=f4)
        elif x > 4 * seqsperfile:
            print('>seq' + str(x) + '-' + str(fbbt[x]), file=f5)
            print(fbbs[x], file=f5)
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()

    print('Total Sequences --' + str(totalseqs) + ' w/ ' + str(seqsperfile) + ' seqs per file')


# Given the Percentage of highest negative and positive values you want
# Function removes all other values from the J Matrix
def Rm_Vals_Percentage_J(J, pct, N, q, **kwargs):
    type = 'both'
    for key, value in kwargs.items():
        if key == 'type':
            type = value
    tmpJ = copy.deepcopy(J)
    pct1valsJ = math.ceil(N * (N - 1) * (q - 1) * 2 * (pct/100))
    vals = tmpJ.flatten()
    if type == '+' or type == 'both':
        pos = [i for i in vals if i > 0]
        t1pctJpos = 100 - pct1valsJ / len(pos) * 100
        pc = np.percentile(pos, t1pctJpos)
    if type == '-' or type == 'both':
        neg = [i for i in vals if i < 0]
        t1pctJneg = pct1valsJ / len(neg) * 100
        nc = np.percentile(neg, t1pctJneg)
    for i in range(N - 1):
        for j in range(N - 1):
            for k in range(q):
                for l in range(q):
                    if i > j:
                        tmpJ[i, j, k, l] = 0.0
                    if type == 'both':
                        if pc > J[i, j, k, l] > nc:
                            tmpJ[i, j, k, l] = 0.0
                    elif type == '+':
                        if pc > J[i, j, k, l]:
                            tmpJ[i, j, k, l] = 0.0
                    elif type == '-':
                        if J[i, j, k, l] > nc:
                            tmpJ[i, j, k, l] = 0.0
    return tmpJ


def Rm_Vals_Percentage_H(H, pct, N, q):
    tmpH = copy.deepcopy(H)
    pct1valsH = math.ceil(N * (q - 1) * (pct/100))
    hvals = tmpH.flatten()
    hpos = [i for i in hvals if i > 0]
    hneg = [i for i in hvals if i < 0]
    t1pctHpos = 100 - pct1valsH / len(hpos) * 100
    t1pctHneg = pct1valsH / len(hneg) * 100
    hpc = np.percentile(hpos, t1pctHpos)
    hnc = np.percentile(hneg, t1pctHneg)
    for i in range(N - 1):
        for j in range(q):
            if hpc > tmpH[i, j] > hnc:
                tmpH[i, j] = 0.0
    return tmpH


# Method to split up files the way i want them too # lol


########################################################################################################################
########################################################################################################################
# General Methods

# IN PROGRESS
def sortjmat_blDCA(file, N, q):
    o = open(file, 'r')
    jmatu = np.full((N, N, q, q), 0.0)
    for line in o:
        data = line.split(',')
        jmatu[int(data[0]), int(data[1]), int(data[2]), int(data[3])] = float(data[4].rstrip())
    o.close()
    # Fix up the Matrix to be exactly the same as the plmDCA Matrix
    jmate = np.full((N - 1, N - 1, q, q), 0.0)
    for i in range(N - 1):
        for j in range(N - 1):
            if i > j:
                continue
            jmate[i, j, :, :] = jmatu[i, j + 1, :, :]
    normed = Normalize_JMatrix(np.negative(jmate), N, q)
    return normed


# Takes plmDCA J Matrix File and inputs the values into a N-1, N-1, q, q matrix
def sortjmat_plmDCA(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N - 1, N - 1, q, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]) - 1, int(data[1]) - 2, int(data[2]) - 1, int(data[3]) - 1] = float(data[4].rstrip())
    o.close()
    fullmatrix[:, :, 0, :] = 0.0
    fullmatrix[:, :, :, 0] = 0.0
    #normed = Normalize_JMatrix(fullmatrix, N, q)
    return fullmatrix

def sort_2d_gen_matrix(file, x, y):
    o = open(file, 'r')
    filled = 0
    fullmatrix = np.full((x, y), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]), int(data[1])] = float(data[2].rstrip())
        if float(data[2].rstrip()) > 0:
            filled += 1
    o.close()
    print(filled)
    return fullmatrix

def sort_4d_gen_matrix(file, x, y, w, z):
    o = open(file, 'r')
    filled = 0
    fullmatrix = np.full((x, y, w, z), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]), int(data[1]), int(data[2]), int(data[3])] = float(data[4].rstrip())
        if float(data[4].rstrip()) > 0:
            filled += 1
    o.close()
    print(filled)
    return fullmatrix

def sorthmat_blDCA(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]), int(data[1])] = float(data[2].rstrip())
    o.close()
    normed = Normalize_HMatrix(np.negative(fullmatrix), N, q)
    return normed


def sorthmat_plmDCA_autoNandq(file):
    o = open(file, 'r')
    linelist = o.readlines()
    o.close()
    ind = linelist[-1].split(',')
    N = int(ind[0])
    q = int(ind[1])
    fullmatrix = np.full((N, q), 0.0)
    for line in linelist:
        data = line.split(',')
        fullmatrix[int(data[0]) - 1, int(data[1]) - 1] = float(data[2].rstrip())
    return fullmatrix, N, q


def Get_Nq_TCR(clustid):
    upath = "/home/jonah/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/FullSeq/"
    fpath = upath + 'Clust' + str(clustid) + '/' + str(clustid) + 'full.h'
    o = open(fpath, 'r')
    linelist = o.readlines()
    o.close()
    ind = linelist[-1].split(',')
    N = int(ind[0])
    q = int(ind[1])
    return N, q


# Takes plmDCA H Matrix File and inputs the values into a N-1, q matrix
def sorthmat_plmDCA(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]) - 1, int(data[1]) - 1] = float(data[2].rstrip())
    o.close()
    fullmatrix[:, 0] = 0.0
    return fullmatrix


# Returns H Matrix with only the top values specified by a percentile
# by default returns values greater than the 80th percentile
# Keyword argument percentile sets the percentile
def TopH_values_disp(H, N, q, **kwargs):
    percentile = 80
    for key, value in kwargs.items():
        if key == 'percentile':
            percentile = value
        else:
            print('No keyword argument ' + key + ' found')
    Hdisp = np.full((N, q), 0.0)
    val = np.percentile(H, percentile)
    for i in range(0, N):
        for j in range(0, q):
            if H[i, j] > val:
                Hdisp[i, j] = H[i, j]
    return Hdisp


# Returns J Matrix with only the top NORMS of individual Jij Matrices specified by a percentile
# by default returns values greater than the 80th percentile
# Keyword argument percentile sets the percentile
# Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# norm values returned as list of vals and cutoff returned as single value
def TopJ_Norms_disp(J, N, **kwargs):
    percentile = 80
    dist = False
    for key, value in kwargs.items():
        if key == 'percentile':
            percentile = value
        elif key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    jnorm = np.full((N - 1, N - 1), 0.0)
    jdisp = np.full((N - 1, N - 1), 0.0)
    vals = []
    for i in range(N - 1):
        for j in range(N - 1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
    tval = np.percentile(jnorm, percentile)
    for i in range(N - 1):
        for j in range(N - 1):
            if jnorm[i, j] >= tval:
                jdisp[i, j] = jnorm[i, j]
            if jnorm[i, j] != 0.0:
                vals.append(jnorm[i, j])
    if dist:
        return jdisp, tval, vals
    else:
        return jdisp


# Return Full J Matrix in 2D so it can be used by matplotlib and shown easily
def FullJ_disp(J, N, q):
    Jdisp = np.full(((N - 1) * q, (N - 1) * q), 0.0)
    for i in range(N - 1):
        for j in range(N - 1):
            for k in range(q):
                for l in range(q):
                    if J[i, j, k, l] != 0.0:
                        Jdisp[i * q + k, j * q + l] = J[i, j, k, l]
                    else:
                        Jdisp[i * q + k, j * q + l] = 0.0
    return Jdisp

def PortionJ_disp(J, N, q, si, ei, sj, ej):
    idim, jdim = ei - si, ej - sj
    Jdisp = np.full(((idim) * q, (jdim) * q), 0.0)
    for i in range(si, ei):
        for j in range(sj, ej):
            for k in range(q):
                for l in range(q):
                    if J[i, j, k, l] != 0.0:
                        Jdisp[(i-si) * q + k, (j-sj) * q + l] = J[i, j, k, l]
                    else:
                        Jdisp[(i-si) * q + k, (j-sj) * q + l] = 0.0
    return Jdisp
# Returns J Matrix with H values added
def HJ_Mutant(J, H, N, q):
    mutt = copy.deepcopy(J)
    for x in range(N - 1):
        for k in range(q):  # J Indices
            mutt[x, x:N, k, :] += H[x, k]
    for y in range(N - 1):
        for l in range(q):  # y states
            mutt[0:y + 1, y, :, l] += H[y + 1, l]
    return mutt


# Takes in Good Binders J and H and Bad Binders J and H
# if htype = 'ind' compares top individual values of GB H and BB H and if they are in common removes them
# if htype = 'norm' compares top norms of GB H and BB H and if they are in common removes them
# Set pcts with hnormpct and jnormpct
def Binder_Comp_JH(J, bJ, H, bH, N, q, **kwargs):
    filledJnorms = (N - 1) * (N - 2) / 2 + N - 1
    filledHind = (N * q)
    nxj = 10  # Number of J Norms used in Comparison
    nxh = 10  # Number of H Norms used n Comparison
    htype = 'good'
    hdist = False
    jdist = False
    for key, value in kwargs.items():
        if key == 'jnormpct':
            nxj = math.ceil(value / 100 * filledJnorms)
        elif key == 'htype':
            htype = value
        elif key == 'hnormpct':
            if htype == 'norm':
                nxh = math.ceil(value / 100 * N)
            if htype == 'ind':
                nxh = math.ceil(value / 100 * filledHind)
        elif key == 'hdist':
            hdist = value
        elif key == 'jdist':
            jdist = value
        else:
            print('No keyword argument ' + key + ' found')
    diff = np.subtract(J, bJ)
    diffnorms = []
    for x in range(N - 1):
        for y in range(N - 1):
            if x > y:
                continue
            diffnorms.append(np.linalg.norm(diff[x, y, :, :]))
    cut = np.percentile(diffnorms, 50)
    topJ, valsJ, tvalJ = TopX_JNorms(J, N, nxj, pct=20, dist=True)
    topBJ, valsBJ, tvalBJ = TopX_JNorms(bJ, N, nxj, pct=20, dist=True)
    for xj, yj, val in topJ:
        for xb, yb, valb in topBJ:
            if xj == xb and yj == yb:
                if np.linalg.norm(diff[xj, yj, :, :]) < cut:
                    J[xj, yj, :, :] = 0.0
    if htype != 'good':
        topH, valsH, tvalH = TopX_HVals(H, N, nxh, pct=20, htype=htype, dist=True)
        topBH, valsBH, tvalBH = TopX_HVals(bH, N, nxh, pct=20, htype=htype, dist=True)
        if htype == 'norm':
            for xi, val in topH:
                for yb, valb in topBH:
                    if xi == yb:
                        H[xi, :] = 0.0
        elif htype == 'ind':
            for xi, yi, val in topH:
                for xb, yb, valb in topBH:
                    if xi == xb and yi == yb:
                        H[xi, yi] = 0.0
        if hdist:
            Hdist = [valsH, tvalH, valsBH, tvalBH]
    else:
        Hdist = [0]
    if jdist:
        Jdist = [valsJ, tvalJ, valsBJ, tvalBJ]
    if hdist and jdist:
        return H, J, Hdist, Jdist
    elif hdist and not jdist:
        return H, J, hdist
    elif not hdist and jdist:
        return H, J, jdist
    elif not hdist and not jdist:
        return H, J


# Reads Fasta files and specifically takes affinity from how its outputted in seqprep
def Fasta_Read_Aff(fastafile):
    o = open(fastafile)
    titles = []
    seqs = []
    for line in o:
        if line.startswith('>'):
            titles.append(float(line.rstrip().split('-')[1]))
        else:
            seqs.append(line.rstrip())
    o.close()
    return titles, seqs

def Fasta_Read_Aff_wC(fastafile, affc):
    baffs, bseqs = [], []
    affs, seqs = Fasta_Read_Aff(fastafile)
    for xid, aff in enumerate(affs):
        if aff >= affc and len(seqs[xid]) == 40:
            baffs.append(aff)
            bseqs.append(seqs[xid])
    return baffs, bseqs

def classify_seq(seql):
    if ''.join(seql[0:4]) == 'AGGG' and ''.join(seql[5:10]) == 'TGATG' and ''.join(seql[11:20]) == 'GTGGTAGGC':
        sid = 0
    elif ''.join(seql[0:8]) == 'AGGGTAGG' and ''.join(seql[9:20]) == 'GTGGATGATGC':
        sid = 1
    elif ''.join(seql[0:6]) == 'TAGGTT' and ''.join(seql[7:20]) == 'TGGGTAGGGTGGT':
        sid = 2
    elif ''.join(seql[1:5]) == 'TAGG' and ''.join(seql[6:19]) == 'TGGGTAGGGTGGT':
        sid = 3
    else:
        sid = 4
    return sid

def classify_fullseq(seql):
    lseq = seql[0:20]
    rseq = seql[20:40]
    lid, rid = classify_seq(lseq), classify_seq(rseq)
    return lid, rid


def Fasta_Read_Aff_fams(fastafile):
    o = open(fastafile)
    titles = []
    lseqs, rseqs = [[], [], [], [], []], [[], [], [], [], []]
    seqs = []
    for line in o:
        if line.startswith('>'):
            aff = float(line.rstrip().split('-')[1])
            titles.append(aff)
        else:
            bseq = list(line.rstrip())
            if len(bseq) == 40 and aff > 1:
                lid, rid = classify_fullseq(bseq)
                # lseq = bseq[0:20]
                # rseq = bseq[20:40]
                # lid, rid = classify_seq(lseq), classify_seq(rseq)
                if lid == 3 and rid == 3:
                    seqs.append((''.join(bseq), aff))
                # lseqs[lid].append((''.join(lseq), aff))
                # rseqs[rid].append((''.join(rseq), aff))
    o.close()
    # return lseqs, rseqs
    return seqs

def write_fams(outp, fastafile):
    lseqs, rseqs = Fasta_Read_Aff_fams(fastafile)
    outlprefix = ['la.txt', 'lb.txt', 'lc.txt', 'ld.txt', 'le.txt']
    outrprefix = ['ra.txt', 'rb.txt', 'rc.txt', 'rd.txt', 're.txt']
    outl = [open(outp + x, 'w') for x in outlprefix]
    outr = [open(outp + x, 'w') for x in outrprefix]
    for i in range(5):
        count = 0
        for s, a in lseqs[i]:
            count += 1
            print('>seq', count, '-', a, sep='', file=outl[i])
            print(s, file=outl[i])
        count = 0
        for s, a in rseqs[i]:
            count += 1
            print('>seq', count, '-', a, sep='', file=outr[i])
            print(s, file=outr[i])
        outl[i].close()
        outr[i].close()





def Fasta_Read_GB(fastafile):
    baffs, bseqs = [], []
    affs, seqs = Fasta_Read_Aff(fastafile)
    for xid, aff in enumerate(affs):
        if aff != 1:
            baffs.append(aff)
            bseqs.append(seqs[xid])
    return baffs, bseqs

def Fasta_Read_BB(fastafile):
    baffs, bseqs = [], []
    affs, seqs = Fasta_Read_Aff(fastafile)
    for xid, aff in enumerate(affs):
        if aff == 1:
            baffs.append(aff)
            bseqs.append(seqs[xid])
    return baffs, bseqs


def Fasta_Read_Aff_wE_above80(fastafile, cutoff, J, H):
    o = open(fastafile)
    titles = []
    seqs = []
    for line in o:
        if line.startswith('>'):
            titles.append(float(line.rstrip().split('-')[1]))
        else:
            e = Calc_Energy(line.rstrip(), J, H)
            if e > cutoff:
                seqs.append(line.rstrip())
            else:
                del titles[-1]
    o.close()
    return titles, seqs


def Fasta_Read_MC(fastafile):
    o = open(fastafile)
    energies = []
    seqs = []
    for line in o:
        if line.startswith('s'):
            energies.append(float(line.rstrip().split('-')[1]))
        else:
            seqs.append(line.rstrip())
    o.close()
    return energies, seqs




def Fasta_Read_SeqOnly(fastafile):
    o = open(fastafile)
    seqs = []
    for line in o:
        if line.startswith('>'):
            continue
        else:
            seqs.append(line.rstrip())
    o.close()
    return seqs


# Calculates 'energy' of give sequence according to provided J and H Matrices
def Calc_Energy(seq, J, H):
    full = list(seq)
    dist = len(full)
    Jenergy = 0
    Henergy = 0
    for x in range(0, dist):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        for y in range(x + 1, dist):
            jbase = rnad[seq[y]]
            Jenergy += J[x, y - 1, ibase, jbase]
    energy = Jenergy + Henergy
    # print(Jenergy)
    # print(Henergy)
    return energy


def Point_mutation_better_binder_checker(seq, J, H):
    startingE = Calc_Energy(seq, J, H)
    tseq = list(copy.deepcopy(seq))
    for iid, i in enumerate(tseq):
        alternatives = [x for x in nucs if x != i]
        for y in alternatives:
            tseq[iid] = y
            tenergy = Calc_Energy(tseq, J, H)
            print(iid+1, i, y, tenergy)
            if tenergy > startingE:
                print(iid+1, i, 'to', y)
        tseq[iid] = i

def TwoMutation_Bad_binder_checker(seq, J, H, allseqpath):
    startingE = Calc_Energy(seq, J, H)
    tseq = list(copy.deepcopy(seq))
    results = []
    for iid, i in enumerate(tseq):
        alternatives = [x for x in nucs if x != i]
        altseqs = []
        for an in alternatives:
            tseq[iid] = an
            altseqs.append(copy.deepcopy(tseq))
        for altseq in altseqs:
            for yid, y in enumerate(altseq):
                alt2 = [x for x in nucs if x != y]
                for j in alt2:
                    altseq[yid] = j
                    tenergy = Calc_Energy(altseq, J, H)
                    # tmpa, simscore, tmpb = ensemble_checker(allseqpath, altseq)[0]
                    # if simscore != 0.95:
                    results.append((''.join(altseq), tenergy))
                altseq[yid] = y
        tseq[iid] = i
        results.sort(key=lambda tup: tup[1])
        highest = results[1:5]
        lowest = results[-5:-1]
        return highest, lowest


def Weighting_Positions_highEseq(seq, J, H):
    tseq = list(copy.deepcopy(seq))
    avgs, score = [], []
    for iid, i in enumerate(tseq):
        alternatives = [x for x in nucs if x != i]
        altE = []
        for y in alternatives:
            tseq[iid] = y
            tenergy = Calc_Energy(tseq, J, H)
            altE.append(tenergy)
        avgs.append(mean(altE))
        tseq[iid] = i
    tavg = mean(avgs)
    savgs = copy.deepcopy(avgs)
    savgs.sort()
    results = np.full((len(tseq), 2), 0.0)
    print(avgs)
    print(savgs)
    for x in savgs[:5]:
        xid = avgs.index(x)
        score = -x/tavg
        results[xid, 0], results[xid, 1] = xid+1, score
    for x in savgs[-5:]:
        xid = avgs.index(x)
        score = x/tavg
        results[xid, 0], results[xid, 1] = xid+1, score
    return results




def Calc_Energy_Breakdown(seq, J, H):
    full = list(seq)
    dist = len(full)
    Jenergy = 0
    Henergy = 0
    iE = np.full((dist, 4), 0.0)
    hpe = 0.0
    for x in range(0, dist):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        iE[x, 0] += H[x, ibase]
        iE[x, 3] = x+1
        for y in range(x + 1, dist):
            jbase = rnad[seq[y]]
            Jenergy += J[x, y - 1, ibase, jbase]
            if J[x, y-1, ibase, jbase] > hpe:
                hpe = Jenergy
                highestpair = [x+1, y+1]
            iE[x, 1] += 0.5*J[x, y - 1, ibase, jbase]
            iE[y, 1] += 0.5*J[x, y - 1, ibase, jbase]
        iE[x, 2] = iE[x, 0] + iE[x, 1]
    energy = Jenergy + Henergy
    return energy, iE, highestpair


def Calc_Energy_Breakdown_Ind_Position(seq, J, H, pos):
    full = list(seq)
    dist = len(full)
    Jenergy = 0
    Henergy = 0
    iE = np.full((dist, 4), 0.0)
    for x in range(0, dist):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        iE[x, 0] += H[x, ibase]
        iE[x, 3] = x+1
        for y in range(x + 1, dist):
            jbase = rnad[seq[y]]
            Jenergy += J[x, y - 1, ibase, jbase]
            iE[x, 1] += 0.5*J[x, y - 1, ibase, jbase]
            iE[y, 1] += 0.5*J[x, y - 1, ibase, jbase]
        iE[x, 2] = iE[x, 0] + iE[x, 1]
    energy = Jenergy + Henergy
    return energy, iE

# N and q correspond to the matrix the
def Calc_Energy_TCR(seq, J, H, N):
    full = list(seq)
    Jenergy = 0
    Henergy = 0
    for x in range(N):
        try:
            ibase = aad[full[x]]
        except IndexError:
            break
        Henergy += H[x, ibase]
        for y in range(x + 1, N):
            try:
                jbase = aad[full[y]]
            except IndexError:
                break
            Jenergy += J[x, y - 1, ibase, jbase]
    energy = Jenergy + Henergy
    # print(Jenergy)
    # print(Henergy)
    return energy


def Calc_Energy3(seq, J, H, K):
    full = list(seq)
    dist = len(full)
    Jenergy = 0
    Henergy = 0
    Kenergy = 0
    for x in range(dist):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        for y in range(x + 1, dist):
            jbase = rnad[seq[y]]
            Jenergy += J[x, y - 1, ibase, jbase]
            if y <= dist-2:
                for z in range(x+2, dist):
                    kbase = rnad[seq[z]]
                    Kenergy += K[x, y-1, z-2, ibase, jbase, kbase]
    energy = Jenergy + Henergy + Kenergy
    print(Jenergy)
    print(Henergy)
    print(Kenergy)
    return energy


# Returns the highest x amount of J Norms
# Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# Highest norms returned in descending list of tuples (i, j, normval) by magnitude of normval
# norm values returned as list of vals and cutoff returned as single value
# Keyword argument percentile determines the percentile used as a cutoff
def TopX_JNorms(J, N, x, **kwargs):
    dist = False
    for key, value in kwargs.items():
        if key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    jnorm = np.full((N - 1, N - 1), 0.0)
    vals = []
    jvals = []
    for i in range(N - 1):
        for j in range(N - 1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
            if jnorm[i, j] != 0.0:
                vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
                jvals.append(jnorm[i, j])
    vals.sort(key=lambda tup: tup[2])
    ind = int(-x)
    top10 = vals[ind:]
    top10.reverse()
    tval = top10[-1][2]
    if dist:
        return top10, jvals, tval
    else:
        return top10


def Strip_nearest_neighbors(jmat, N, q):
    for x in range(N-1):
        for y in range(N-1):
            for k in range(q):
                for l in range(q):
                    if abs(x - y) > 2:
                        continue
                    else:
                        jmat[x, y, k, l] = 0.0


def TopX_Pos_JNorms(J, N, q, x):
    pct = 80
    jnorm = np.full((N - 1, N - 1), 0.0)
    vals = []
    Jpos = Sign_Seperator(J, N, q, mattype='j', sign='+')
    for i in range(N - 1):
        for j in range(N - 1):
            jnorm[i, j] = np.linalg.norm(Jpos[i, j, :, :])
            if jnorm[i, j] != 0.0:
                vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
    vals.sort(key=lambda tup: tup[2])
    ind = int(-x)
    top10 = vals[ind:-1]
    return top10


def TopX_Hnorms(H, N, x, **kwargs):
    dist = False
    for key, value in kwargs.items():
        if key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    hnorm = np.full((N), 0.0)
    vals = []
    hvals = []
    for i in range(N):
        hnorm[i] = np.linalg.norm(H[i, :])
        if hnorm[i] != 0.0:
            vals.append((i, hnorm[i]))  # 0, 0 -> 1, 2
            hvals.append(hnorm[i])
    tot = N
    tval = np.percentile(hnorm, (tot - x)/tot*100)
    vals.sort(key=lambda tup: tup[1])
    ind = int(-x)
    top10 = vals[ind:-1]
    if dist:
        return top10, hvals, tval
    else:
        return top10


# Returns the highest x amount of H Norms
# Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# norm values returned as list of vals and cutoff returned as single value
# Keyword argument percentile determines the percentile used as a cutoff, default is 80
# Keyword htype has options 'ind', 'norm'
# 'ind' returns the top individual values in H
# 'norm' returns the top norms out of N norms in H
def TopX_HVals(H, N, x, **kwargs):
    pct = 80
    htype = 'ind'
    dist = False
    for key, value in kwargs.items():
        if key == 'percentile':
            pct = value
        elif key == 'htype':
            htype = value
        elif key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    if htype == 'norm':
        hnorm = np.full((N - 1), 0.0)
        vals = []
        hvals = []
        for i in range(N - 1):
            hnorm[i] = np.linalg.norm(H[i, :])
            if hnorm[i] != 0.0:
                vals.append((i, hnorm[i]))
                hvals.append(hnorm[i])
        tval = np.percentile(hvals, pct)
        vals.sort(key=lambda tup: tup[1])
    elif htype == 'ind':
        vals = []
        hvals = []
        for i in range(N - 1):
            for j in range(1, 5):
                vals.append((i, j, abs(H[i, j])))
                hvals.append(H[i, j])
        tval = np.percentile(hvals, pct)
        vals.sort(key=lambda tup: tup[2])
    ind = int(-x)
    top10 = vals[ind:-1]
    if dist:
        return top10, hvals, tval
    else:
        return top10


# Quick Method to get N by reading just the length of the first seq.. will only work if all seqs are the same length
def getn(fastafile):
    o = open(fastafile, 'r')
    o.readline()
    seq = o.readline().rstrip()
    n = len(list(seq))
    o.close()
    return n


########################################################################################################################
########################################################################################################################
# Subplot Methods


# Shows Entire J Matrix
# Check Function for Keyword Arguments, available are lw, fontsize, xlabel, ylabel, vmin, vmax, title, and cmap
def Fig_FullJ(subplot, id, J, n, q, **kwargs):
    cmap = 'seismic'
    vml = -1
    vmg = 1
    lw = 0.1
    xlabel = 'j'
    ylabel = 'i'
    fontsize = 6
    title = 'JMat ID: ' + str(id)
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        elif key == 'lw':
            lw = value
        elif key == 'xlabel':
            xlabel = value
        elif key == 'ylabel':
            ylabel = value
        elif key == 'ticksize':
            fontsize = value
        elif key == 'vmin':
            vml = value
        elif key == 'vmax':
            vmg = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=6)
    subplot.imshow(J, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
    subplot.set_xticks(np.arange(-.5, (n - 2) * q, q))
    subplot.set_yticks(np.arange(-.5, (n - 2) * q, q))
    subplot.set_xticklabels(np.arange(2, n, 1))
    subplot.set_yticklabels(np.arange(1, n - 1, 1))
    subplot.grid(True, color='g', lw=lw)
    subplot.set_ylabel(ylabel)
    subplot.set_xlabel(xlabel)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)


# Produces figure on given subplot of J norms Matrix which is a parameter
# Check Function for Keyword Arguments, available are lw, fontsize, xlabel, ylabel, vmin, vmax, title, and cmap
def Fig_Jnorm(subplot, id, J, n, **kwargs):
    cmap = 'seismic'
    vml = 0
    vmg = 4
    lw = 1.0
    xlabel = 'j'
    ylabel = 'i'
    fontsize = 6
    title = 'Jmat Norms ID: ' + str(id)
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        elif key == 'lw':
            lw = value
        elif key == 'xlabel':
            xlabel = value
        elif key == 'ylabel':
            ylabel = value
        elif key == 'ticksize':
            fontsize = value
        elif key == 'vmin':
            vml = value
        elif key == 'vmax':
            vmg = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=6)
    subplot.imshow(J, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
    subplot.set_xticks(np.arange(-.5, (n - 1), 1))
    subplot.set_yticks(np.arange(-.5, (n - 1), 1))
    subplot.set_xticklabels(np.arange(2, n + 1, 1))
    subplot.set_yticklabels(np.arange(1, n, 1))
    subplot.grid(True, color='g', lw=lw)
    subplot.set_ylabel(ylabel)
    subplot.set_xlabel(xlabel)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)


# Produces Figure of Full H Matrix on a given subplot
# Check Function for Keyword Arguments, available are fontsize, xlabel, ylabel, vmin, vmax, title, and cmap
def Fig_FullH(subplot, id, H, n, q, **kwargs):
    cmap = 'seismic'
    vml = -1
    vmg = 1
    xl = False
    xlabel = 'hello'
    ylabel = 'i'
    fontsize = 6
    title = 'Hmat ID: ' + str(id)
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        elif key == 'xlabel':
            xl = True
            xlabel = value
        elif key == 'ylabel':
            ylabel = value
        elif key == 'fontsize':
            fontsize = value
        elif key == 'vmax':
            vmg = value
        elif key == 'vmin':
            vml = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.imshow(H, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=fontsize)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)
    subplot.set_xticks(np.arange(0, q + 1, 1))
    subplot.set_yticks(np.arange(0, n + 1, 1))
    if q == 21 and xl is False:
        subplot.set_xticklabels(aa)
        subplot.set_xlabel('Amino Acid')
    elif q == 5 and xl is False:
        subplot.set_xticklabels(rna)
        subplot.set_xlabel('Base')
    else:
        subplot.set_xlabel(xlabel)
    subplot.set_yticklabels(np.arange(1, n + 1, 1))
    subplot.set_ylabel(ylabel)


# Produces Figure of Distribution of Values on a given subplot
# A Vertical Line is produced at the Cutoff Value.. typically used in showing which values in a distribution were used
# Check Function for Keyword Arguments, available are fontsize, xlabel, xmin, xmax, title, and pcolor
def Fig_Distribution_w_Cutoff(subplot, id, Values, Cutoff, **kwargs):
    title = 'Distribution ID: ' + str(id)
    plotcolor = 'r'
    xml = 0
    xmg = 2
    xlabel = 'Value'
    fontsize = 6
    for key, value in kwargs.items():
        if key == 'xlabel':
            xlabel = value
        elif key == 'pcolor':
            plotcolor = value
        elif key == 'fontsize':
            fontsize = value
        elif key == 'xmin':
            xml = value
        elif key == 'xmax':
            xmg = value
        elif key == 'title':
            title = value
        else:
            print('No keyword argument ' + key + ' found')
    deN = gaussian_kde(Values)
    xd1 = np.linspace(xml, xmg, 100)
    subplot.plot(xd1, deN(xd1), color=plotcolor)
    subplot.plot(Values, [0.01] * len(Values), '|', color='k')
    subplot.set_xlabel(xlabel)
    subplot.grid(True)
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=fontsize)
    subplot.axvline(x=Cutoff)


# Shows Premade SeqLogo on a given Subplot
# Keyword Arguments are title and fontsize
def Fig_SeqLogo(Filepath, Subplot, id):
    title = 'SeqLogo ID: ' + str(id)
    fontsize = 6
    for key, value in kwargs.items():
        if key == 'title':
            title = value
        elif key == 'fontsize':
            fontsize = value
        else:
            print('No keyword argument ' + key + ' found')
    fsl1 = mpimg.imread(Filepath)
    Subplot.imshow(fsl1)
    Subplot.axis('off')
    Subplot.title.set_text(title)
    Subplot.title.set_size(fontsize=fontsize)


# On Specified Subplot shows Individual Jij at specified x and y **NOTE J12 is equivalent to J[0, 0, :, :]
# Check Function for Keyword Arguments, available are cbar, lw, fontsize, vmin, vmax, title, cmap, and type
# Three Types are available: 'dna', 'rna' and 'pep'
def Fig_IndJij(subplot, J, x, y, id, **kwargs):
    vml = -0.5
    vmg = 0.5
    cmap = 'seismic'
    fontsize = 4
    type = 'rna'
    cbar = False
    lw = 0.1
    title = 'p: ' + str(x + 1) + ' and p: ' + str(y + 2)
    for key, value in kwargs.items():
        if key == 'vmax':
            vmg = value
        elif key == 'vmin':
            vml = value
        elif key == 'cmap':
            cmap = value
        elif key == 'fontsize':
            fontsize = value
        elif key == 'type':
            type = value
        elif key == 'lw':
            lw = value
        elif key == 'title':
            title = value
        elif key == 'cbar':
            cbar = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.imshow(J[x, y, :, :], cmap=cmap, vmin=vml, vmax=vmg)
    if cbar:
        plt.colorbar(pos, ax=subplot, fraction=0.046, pad=0.04)
    if type == 'rna' or type == 'dna':
        subplot.set_xticks(np.arange(-.5, 4.5, 1))
        subplot.set_yticks(np.arange(-.5, 4.5, 1))
        if type == 'rna':
            subplot.set_xticklabels(rna)
            subplot.set_yticklabels(rna)
        else:
            subplot.set_xticklabels(dna)
            subplot.set_yticklabels(dna)
    elif type == 'pep':
        subplot.set_xticks(np.arange(-.5, 20.5, 1))
        subplot.set_yticks(np.arange(-.5, 20.5, 1))
        subplot.set_xticklabels(aa)
        subplot.set_yticklabels(aa)
    subplot.tick_params(axis='both', which='major', labelsize=fontsize)
    subplot.tick_params(axis='both', which='minor', labelsize=fontsize)
    subplot.grid(True, color='r', lw=lw)
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=(fontsize + 2))


def TopJNorms_Jmatrix(J, N, q, numberofnorms):
    TJ = np.full((N-1, N-1, q, q), 0.0)
    t10, jvals, dist = TopX_JNorms(J, N, numberofnorms, dist='True')
    for x, y, val in t10:
        TJ[x, y, :, :] = copy.deepcopy(J[x, y, :, :])
    return TJ, jvals, dist

def TopHNorms_Hmatrix(H, N, q, numberofnorms):
    TH = np.full((N, q), 0.0)
    t10, jvals, dist = TopX_Hnorms(H, N, numberofnorms, dist='True')
    for x, val in t10:
        TH[x, :] = copy.deepcopy(H[x, :])
    return TH, jvals, dist

########################################################################################################################
########################################################################################################################
# Full Figure Methods


# Returns a figure showing the Individual Jijs of the top 10 Norms of the provided J Matrix
# OutPath is the directory the figure is being saved to
def Top10norms_figure_RNA(id, J, N, OutPath):
    # Get Indices of top 10 norms
    jx = TopX_JNorms(J, N, 30)
    fig, ax = plt.subplots(2, 5, constrained_layout=True)
    for i in range(10):
        x, y, z = jx[i]
        j = i % 5
        k = 0
        if i == 0:
            Fig_IndJij(ax[k, j], J, x, y, id, vmin=-0.5, vmax=.5, type='rna', cbar=False)
        else:
            if i > 4:
                k = 1
            Fig_IndJij(ax[k, j], J, x, y, id, vmin=-0.5, vmax=.5, type='rna')

    fig.suptitle('Highest Jij Norms ID: ' + str(id))
    plt.savefig(OutPath + str(id) + 'JNormt10.png', dpi=600)

def Top10norms_figure_DNA(id, J, N, OutPath):
    # Get Indices of top 10 norms
    jx = TopX_JNorms(J, N, 20)
    fig, ax = plt.subplots(2, 5, constrained_layout=True)
    for i in range(10):
        x, y, z = jx[i]
        j = i % 5
        k = 0
        if i == 0:
            Fig_IndJij(ax[k, j], J, x, y, id, vmin=-0.5, vmax=.5, type='dna', cbar=False)
        else:
            if i > 4:
                k = 1
            Fig_IndJij(ax[k, j], J, x, y, id, vmin=-0.5, vmax=.5, type='dna')

    fig.suptitle('Highest Jij Norms ID: ' + str(id))
    plt.savefig(OutPath + str(id) + 'JNormt10.png', dpi=600)

# Can input multiple fasta files and it will combine all seqs, and plot their energies according to provided J and H
# keyword argument is title
def Plot_Seq_Aff_v_E(J, H, outpath, *argv, **kwargs):
    title = 'Affinity vs Energy'
    for key, value in kwargs.items():
        if key == 'title':
            title = value
    titles = []
    seqs = []
    for arg in argv:
        tmpt, tmps = Fasta_Read_Aff(arg)
        titles.extend(tmpt)
        seqs.extend(tmps)
    energies = []
    for x in seqs:
        nrg = Calc_Energy(x, J, H)
        energies.append(nrg)
    api = list(zip(titles, energies))
    x = list(set([x for (x, y) in api]))
    x.sort()
    avg = []
    err = []
    highestaff = 1
    for aff in x:
        if aff > highestaff:
            highestaff = aff
        yvals = np.array([y for (x, y) in api if x == aff])
        yavg = yvals.mean()
        yerr = np.std(yvals)
        avg.append(yavg)
        err.append(yerr)
    linreg = stats.linregress(x, avg)
    xl = np.linspace(0, highestaff, 100)
    plt.plot(xl, xl*linreg[0]+linreg[1], ':r')
    plt.errorbar(x, avg, err, linestyle='None', marker='^')
    plt.xlabel('R-Score: ' + str(linreg[2]))
    plt.ylabel('Energy')
    plt.suptitle(title)
    plt.savefig(outpath, dpi=600)


def seq_breakdown_by_aff(seqfile, J, H, aff):
    affs, seqs = Fasta_Read_Aff(seqfile)
    energies, seqoi = [], []
    for sid, seq in enumerate(seqs):
        if affs[sid] == aff:
            energies.append(Calc_Energy(seq, J, H))
            seqoi.append(seq)
    lv = min(energies)
    hv = max(energies)
    rv = (hv - lv)/50
    results = np.full((50, 2), 0.0)
    for i in range(1,51):
        cutoff = lv+i*rv
        pct = len([x for x in energies if x >= cutoff])/len(energies)
        results[i-1, 0], results[i-1, 1] = .02*i*hv, pct*len(energies)
    return results


def ensemble_checker(seqfile, *seqs):
    eaffs, eseqs = Fasta_Read_Aff(seqfile)
    results = []
    for seq in seqs:
        tseq = rna2dna(seq)
        seqoi = list(tseq)
        highest_sim_score = 0.
        msseq = seq
        for eseq in eseqs:
            es = list(eseq)
            seq_similarity = 0.
            for i in range(len(es)):
                if seqoi[i] == es[i]:
                    seq_similarity += 1.
            seq_similarity /= len(seq)
            if seq_similarity > highest_sim_score:
                highest_sim_score = seq_similarity
                msseq = eseq
        results.append((tseq, highest_sim_score, msseq))
    return results

def get_affinity(seqfile, seq):
    eaffs, eseqs = Fasta_Read_Aff(seqfile)
    results = []
    tseq = rna2dna(seq)
    id = eseqs.index(tseq)
    return eaffs[id]



def Jnorm_finder(J, N, q, sn, en, step, seqfile, z):
    norml = np.arange(sn, en, step)
    affs, seqs = Fasta_Read_Aff(seqfile)
    emptyH = np.full((N, q), 0.0)
    cR = -1.0
    results = []
    for norm_number in norml:
        tJ, jvals, dist = TopJNorms_Jmatrix(J, N, q, norm_number)
        rscore = R_SCORE(affs, seqs, emptyH, tJ)
        if rscore > cR:
            cR = rscore
            results.append((norm_number, rscore))
    z.put(results)


def Hnorm_finder(H, N, q, sn, en, step, seqfile, z):
    norml = np.arange(sn, en, step)
    affs, seqs = Fasta_Read_Aff(seqfile)
    emptyJ = np.full((N, N, q, q), 0.0)
    cR = -1.0
    results = []
    for norm_number in norml:
        tH, hvals, dist = TopHNorms_Hmatrix(H, N, q, norm_number)
        rscore = R_SCORE(affs, seqs, tH, emptyJ)
        if rscore > cR:
            cR = rscore
            results.append((norm_number, rscore))
    z.put(results)


def Jnorm_finder_bJ(bJ, N, q, sn, en, step, seqfile, z):
    norml = np.arange(sn, en, step)
    affs, seqs = Fasta_Read_Aff(seqfile)
    emptyH = np.full((N, q), 0.0)
    cR = 1
    results = []
    for norm_number in norml:
        tJ, jvals, dist = TopJNorms_Jmatrix(bJ, N, q, norm_number)
        rscore = R_SCORE(affs, seqs, emptyH, tJ)
        if rscore < cR:
            cR = rscore
            results.append((norm_number, rscore))
    z.put(results)


def rna2dna(seq):
    tseq = ''
    for iid, i in enumerate(list(seq)):
        if i == 'U':
            tseq += 'T'
        else:
            tseq += i
    return tseq

########################################################################################################################
########################################################################################################################
# Generate Seq Methods

def check_vals(x, y, pvals):
    if len(pvals) == 0:
        return 0, 0, 'none'
    else:
        for idd, pack in enumerate(pvals):
            xp, yp, rx, ry, pval = pack
            if x == xp:
                return 1, idd, 'xs'
            elif y == yp:
                return 1, idd, 'ys'
            elif x == yp:
                return 1, idd, 'xnyp'
            elif y == xp:
                return 1, idd, 'ynxp'
        return 0, 0, 'none'


# Donp.full(())6666esn't Work for the End Scoring Methods
def past_entry_comp_goodseq(J, pvals, xn, yn):
    tmppvals = copy.deepcopy(pvals)
    ind, xid, stype = check_vals(xn, yn, pvals)
    if ind == 0:
        pos = [x for x in J[xn, yn, :, :].flatten() if x > 0]
        if len(pos) > 0:
            # tmpxn, tmpyn = list(np.where(J[xn, yn, :, :] == np.amax(J[xn, yn, :, :])))
            tmp = np.argmax(J[xn, yn, :, :].flatten())
            print(tmp)
            tmpxn = int(math.floor(tmp/5))
            tmpyn = int(tmp % 5)
            print(tmpxn)
            print(tmpyn)
            rxn = tmpxn
            ryn = tmpyn
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
            return [xn, yn, rxn, ryn], tmppvals
        else:
            tmppvals.append((100, 100, 0, 0, 0))
            return [100, 100, 0, 0, 0], tmppvals
    elif ind == 1:
        xp, yp, rxp, ryp, val = tmppvals[xid]
        pos = [x for x in J[xn, yn, :, :].flatten() if x > 0]
        tmpxn, tmpyn = list(np.where(J[xn, yn, :, :] == np.amax(J[xn, yn, :, :])))  # Indices of highest in N
        rxn = int(tmpxn)
        ryn = int(tmpyn)
        if stype == 'xs':
            pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, rxp, :])
            tchoice = np.amax(J[xp, yp, rxn, :]) + J[xn, yn, rxn, ryn]

        if stype == 'ys':
            pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, :, ryp])
            tchoice = np.amax(J[xp, yp, :, ryn]) + J[xn, yn, rxn, ryn]

        if stype == 'xnyp':
            pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, ryp, :])
            tchoice = np.amax(J[xp, yp, ryn, :]) + J[xn, yn, rxn, ryn]

        if stype == 'ynxp':
            pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, :, rxp])
            tchoice = np.amax(J[xp, yp, :, rxn]) + J[xn, yn, rxn, ryn]

        if pchoice > tchoice:
            if stype == 'xs':
                rxn = rxp
                ryn = int(np.where(J[xn, yn, rxn, :] == np.amax(J[xn, yn, rxn, :]))[0])
            if stype == 'ys':
                ryn = ryp
                rxn = int(np.where(J[xn, yn, :, ryn] == np.amax(J[xn, yn, :, ryn]))[0])
            if stype == 'xnyp':
                rxn = ryp
                ryn = int(np.where(J[xn, yn, rxn, :] == np.amax(J[xn, yn, rxn, :]))[0])
            if stype == 'ynxp':
                ryn = rxp
                rxn = int(np.where(J[xn, yn, :, ryn] == np.amax(J[xn, yn, :, ryn]))[0])
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice > pchoice:
            if stype == 'xs':
                rxp = rxn
                ryp = int(np.where(J[xp, yp, rxp, :] == np.amax(J[xp, yp, rxp, :]))[0])
            if stype == 'ys':
                ryp = ryn
                rxp = int(np.where(J[xp, yp, :, ryp] == np.amax(J[xp, yp, :, ryp]))[0])
            if stype == 'xnyp':
                ryp = rxn
                rxp = int(np.where(J[xp, yp, :, ryp] == np.amax(J[xp, yp, :, ryp]))[0])
            if stype == 'ynxp':
                rxp = ryn
                ryp = int(np.where(J[xp, yp, rxp, :] == np.amax(J[xp, yp, rxp, :]))[0])
            tmppvals[xid] = (xp, yp, rxp, ryp, J[xp, yp, rxp, ryp])
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice == pchoice:
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        vals = [xp, yp, rxp, ryp, xn, yn, rxn, ryn]
        return vals, tmppvals


def past_entry_comp_badseq(J, pvals, xn, yn):
    tmppvals = copy.deepcopy(pvals)
    ind, xid, stype = check_vals(xn, yn, pvals)
    if ind == 0:
        tmpxn, tmpyn = list(np.where(J[xn, yn, 1:5, 1:5] == np.amin(J[xn, yn, 1:5, 1:5])))
        rxn = int(tmpxn) + 1
        ryn = int(tmpyn) + 1
        tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
        return [xn, yn, rxn, ryn], tmppvals
    elif ind == 1:
        xp, yp, rxp, ryp, val = tmppvals[xid]
        tmpxn, tmpyn = list(np.where(J[xn, yn, 1:5, 1:5] == np.amin(J[xn, yn, 1:5, 1:5])))  # Indices of highest in N
        rxn = int(tmpxn) + 1
        ryn = int(tmpyn) + 1
        if stype == 'xs':
            pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, rxp, 1:5])
            tchoice = np.amin(J[xp, yp, rxn, 1:5]) + J[xn, yn, rxn, ryn]

        if stype == 'ys':
            pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, 1:5, ryp])
            tchoice = np.amin(J[xp, yp, 1:5, ryn]) + J[xn, yn, rxn, ryn]

        if stype == 'xnyp':
            pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, ryp, 1:5])
            tchoice = np.amin(J[xp, yp, ryn, 1:5]) + J[xn, yn, rxn, ryn]

        if stype == 'ynxp':
            pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, 1:5, rxp])
            tchoice = np.amin(J[xp, yp, 1:5, rxn]) + J[xn, yn, rxn, ryn]

        if pchoice < tchoice:
            if stype == 'xs':
                rxn = rxp
                ryn = int(np.where(J[xn, yn, rxn, 1:5] == np.amin(J[xn, yn, rxn, 1:5]))[0]) + 1
            if stype == 'ys':
                ryn = ryp
                rxn = int(np.where(J[xn, yn, 1:5, ryn] == np.amin(J[xn, yn, 1:5, ryn]))[0]) + 1
            if stype == 'xnyp':
                rxn = ryp
                ryn = int(np.where(J[xn, yn, rxn, 1:5] == np.amin(J[xn, yn, rxn, 1:5]))[0]) + 1
            if stype == 'ynxp':
                ryn = rxp
                rxn = int(np.where(J[xn, yn, 1:5, ryn] == np.amin(J[xn, yn, 1:5, ryn]))[0]) + 1
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice < pchoice:
            if stype == 'xs':
                rxp = rxn
                ryp = int(np.where(J[xp, yp, rxp, 1:5] == np.amin(J[xp, yp, rxp, 1:5]))[0]) + 1
            if stype == 'ys':
                ryp = ryn
                rxp = int(np.where(J[xp, yp, 1:5, ryp] == np.amin(J[xp, yp, 1:5, ryp]))[0]) + 1
            if stype == 'xnyp':
                ryp = rxn
                rxp = int(np.where(J[xp, yp, 1:5, ryp] == np.amin(J[xp, yp, 1:5, ryp]))[0]) + 1
            if stype == 'ynxp':
                rxp = ryn
                ryp = int(np.where(J[xp, yp, rxp, 1:5] == np.amin(J[xp, yp, rxp, 1:5]))[0]) + 1
            tmppvals[xid] = (xp, yp, rxp, ryp, J[xp, yp, rxp, ryp])
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice == pchoice:
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        vals = [xp, yp, rxp, ryp, xn, yn, rxn, ryn]
        return vals, tmppvals


def Seq_edit_past_entry_comp(array, gseq):
    if len(array) >= 4:
        gseq[array[0] + 0] = rna[int(array[2])]
        gseq[array[1] + 1] = rna[int(array[3])]
        if len(array) == 8:
            gseq[array[4] + 0] = rna[int(array[6])]
            gseq[array[5] + 1] = rna[int(array[7])]
    return gseq


def gen_goodseq(J, H, N, norms):
    # Get Indices of top 10 norms
    gseq = np.full(40, ['X'], dtype=str)
    tval = TopX_Pos_JNorms(J, N, 5, norms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_goodseq(J, pvals, x, y)
        xv = vals[0]
        yv = vals[1]
        if xv < N and yv < N:
            gseq = Seq_edit_past_entry_comp(vals, gseq)
    for xid, x in enumerate(gseq):
        if x == 'X':
            gpx = np.argmax(H[xid, 1:5])
            gseq[xid] = rna[int(gpx) + 1]
    print(''.join(gseq))
    print(Calc_Energy(gseq, J, H))
    return ''.join(gseq)


def gen_badseq(J, H, N, norms):
    # Get Indices of top 10 norms
    bseq = np.full(40, ['X'], dtype=str)
    tval = TopX_JNorms(J, N, norms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_badseq(J, pvals, x, y)
        bseq = Seq_edit_past_entry_comp(vals, bseq)
    for xid, x in enumerate(bseq):
        if x == 'X':
            gpx = np.argmin(H[xid, 1:5])
            bseq[xid] = rna[int(gpx) + 1]
    print(pvals)
    print(Calc_Energy(bseq, J, H))
    print(''.join(bseq))
    return ''.join(bseq)


def Average_J(N, q, *argv):
    hybrid = np.full((N - 1, N - 1, q, q), 0.0)
    for i in range(N - 1):
        for j in range(N - 1):
            if i > j:
                continue
            for k in range(1, q):
                for l in range(1, q):
                    val = []
                    for arg in argv:
                        val.append(arg[i, j, k, l])
                    hybrid[i, j, k, l] = np.average(val)
    return hybrid


def Average_H(N, q, *argv):
    hybrid = np.full((N,  q), 0.0)
    for i in range(N):
        for j in range(1, q):
            val = []
            for arg in argv:
                val.append(arg[i, j])
                hybrid[i, j] = np.average(val)
    return hybrid


def Normalize_HMatrix(H, N, q):
    for i in range(N - 1):
        for j in range(1, q):
            H[i, 0] = 0.0
    d = 2. * (H - np.min(H)) / np.ptp(H) - 1
    return d


# Normalizes Values of J matrix to be between -1 and 1
def Normalize_JMatrix(J, N, q):
    for i in range(N - 1):
        for j in range(N - 1):
            for k in range(1, q):
                if i > j:
                    continue
                J[i, j, 0, :] = 0.0
                J[i, j, :, 0] = 0.0
    d = 2. * (J - np.min(J)) / np.ptp(J) - 1
    return d


def Sign_Seperator(mat, N, q, **kwargs):
    matrix = copy.deepcopy(mat)
    sign = 1
    type = 'j'
    for key, value in kwargs.items():
        if key == 'sign':
            if value == '+':
                sign = 1
            elif value == '-':
                sign = 0
            else:
                print('value error')
        elif key == 'mattype':
            type = value
        else:
            print('No key ' + key + ' found')
    if type == 'j':
        for i in range(N - 1):
            for j in range(N - 1):
                for k in range(q):
                    for l in range(q):
                        if matrix[i, j, k, l] < 0 and sign == 1:
                            matrix[i, j, k, l] = 0.0
                        elif matrix[i, j, k, l] > 0 and sign == 0:
                            matrix[i, j, k, l] = 0.0
    elif type == 'h':
        for i in range(N):
            for j in range(q):
                if matrix[i, j] < 0 and sign == 1:
                    matrix[i, j] = 0.0
                elif matrix[i, j] > 0 and sign == 0:
                    matrix[i, j] = 0.0
    return matrix


# Finds the ideal percentage of J values that make up the final J that results in the highest R score
# Uses Jonah's Scoring Method ;)
def Pct_Finder_JS(gJ, bJ, gH, bH, fasta, N, q):
    pctresults = []
    titles, seqs = Fasta_Read_Aff(fasta)
    for x in range(1, 500):
        pct = x/100
        H = 2*gH - bH
        J = Rm_Vals_Percentage_J(2*gJ - bJ, pct, N, q)
        energies = []
        for x in seqs:
            nrg = Calc_Energy(x, J, H)
            energies.append(nrg)
        api = list(zip(titles, energies))
        x = list(set([x for (x, y) in api]))
        x.sort()
        avg = []
        for aff in x:
            yvals = np.array([y for (x, y) in api if x == aff])
            yavg = yvals.mean()
            avg.append(yavg)
        linreg = stats.linregress(x, avg)
        pctresults.append((linreg[2], pct))
    pctresults.sort(key=lambda tup: tup[0])
    print(pctresults[0:2])
    print(pctresults[-3:-1])
    return pctresults[-1][1]


def best_seperation(J, H, N, q, fastafile):
    titles, seqs = Fasta_Read_Aff(fastafile)
    tmpH = H/4
    results = []
    for i in range(2, 300, 2):
        pct = i/100
        print('pct = ' + str(pct))
        tmpJ = Rm_Vals_Percentage_J(J, pct, N, q)
        energies = []
        for x in seqs:
            energies.append(Calc_Energy(x, tmpJ, tmpH))
        api = list(zip(titles, energies))
        affs = list(set(titles))
        datax = []
        datae = []
        for x in affs:
            prospects = [nrg for (aff, nrg) in api if aff == x]
            datax.append(x)
            datae.append(max(prospects))
        linreg = stats.linregress(datax, datae)
        results.append((pct, linreg[2]))
    results.sort(key=lambda tup: tup[1])
    print(results[0:2])
    print(results[-3:-1])
    bpct = results[-1][0]
    tmpJ = Rm_Vals_Percentage_J(J, bpct, N, q)
    Raw_wRscore(tmpJ, tmpH, '/home/jonah/Desktop/bestfit.png', fastafile)
    return bpct


def designed_GBJmatrix(gJ, N, q, fastafile):
    H = np.full((N, q), 0.0)
    finalJ = np.full((N-1, N-1, q, q), 0.0)
    titles, seqs = Fasta_Read_Aff(fastafile)
    currentR = -1
    for i in range(N - 1):
        for j in range(N - 1):
            for k in range(q):
                for l in range(q):
                    if i > j:
                        continue
                    finalJ[i, j, k, l] = gJ[i, j, k, l]
                    rscore = R_SCORE(titles, seqs, H, finalJ)
                    if rscore > currentR:
                        currentR = rscore
                    else:
                        finalJ[i, j, k, l] = 0.0
    print('designed Matrix R score: ' + str(currentR))
    return finalJ


def designed_GBHmatrix(gH, N, q, fastafile):
    finalH = np.full((N, q), 0.0)
    J = np.full((N - 1, N - 1, q, q), 0.0)
    titles, seqs = Fasta_Read_Aff(fastafile)
    currentR = -1
    for i in range(N):
        for j in range(1, q):
            finalH[i, j] = gH[i, j]
            rscore = R_SCORE(titles, seqs, finalH, J)
            if rscore > currentR:
                currentR = rscore
            else:
                finalH[i, j] = 0.0
    print('designed Matrix R score: ' + str(currentR))
    return finalH


def designerJ(N, q, fastafile, **kwargs):
    type = 'dna'
    for key, value in kwargs.items():
        if key == 'type':
            type = value
    H = np.full((N, q), 0.0)
    finalJ = np.full((N - 1, N - 1, q, q), 0.0)
    titles, seqs = Fasta_Read_Aff(fastafile)
    currentR = -1
    if type == 'pep':
        for i in range(N - 1):
            for j in range(N - 1):
                for k in range(q):
                    for l in range(q):
                        if i > j:
                            continue
                        finalJ[i, j, k, l] = 1
                        rscore = R_SCORE(titles, seqs, H, finalJ, type='pep')
                        if rscore > currentR:
                            currentR = rscore
                        else:
                            finalJ[i, j, k, l] = -1
                            rscore = R_SCORE(titles, seqs, H, finalJ, type='pep')
                            if rscore > currentR:
                                currentR = rscore
                            else:
                                finalJ[i, j, k, l] = 0.0
    else:
        for i in range(N - 1):
            for j in range(N - 1):
                for k in range(q):
                    for l in range(q):
                        if i > j:
                            continue
                        finalJ[i, j, k, l] = 1
                        rscore = R_SCORE(titles, seqs, H, finalJ)
                        if rscore > currentR:
                            currentR = rscore
                        else:
                            finalJ[i, j, k, l] = -1
                            rscore = R_SCORE(titles, seqs, H, finalJ)
                            if rscore > currentR:
                                currentR = rscore
                            else:
                                finalJ[i, j, k, l] = 0.0
    print('designed Matrix R score: ' + str(currentR))
    return finalJ


def designerJ_Seperation_TCR(N, q, SEQHANDLER, CoreNum, outfile):
    comm = round((N - 1) / CoreNum)
    ind = [x * comm for x in range(CoreNum)]
    ind.append(N - 1)
    instructions = []
    for i in range(CoreNum):
        instructions.append([i, N, q, ind[i], ind[i + 1], SEQHANDLER])
    p = mp.Pool(CoreNum)
    z = p.starmap(dJ_Worker, instructions)
    r = Write_Output_dJ(z, N, q, outfile)
    print(r)


def Write_Output_dJ(z, N, q, outfile):
    aPos = []
    aNeg = []
    rs = []
    J = np.full((N-1, N-1, q, q), 0)
    for i in z:
        for pos, neg, cR in i:
            for p in pos:
                aPos.append(p)
            for n in neg:
                aNeg.append(n)
            rs.append(cR)
    for d in aPos:
        x, y, i, j = d
        J[x, y, i, j] = 1
    for d in aNeg:
        x, y, i, j = d
        J[x, y, i, j] = -1
    export_J(J, N, q, outfile)
    return rs

class Sequence:
    def __init__(self, N, q, aff, seq='null', cv='no', **kwargs):
        self.type = 'dna'
        self.mat = 'J'
        for key, value in kwargs.items():
            if key == 'type':
                self.type = value
            if key == 'mat':
                self.mat = value
        self.affinity = aff
        if seq == 'null':
            self.seq = np.full(40, ['X'], dtype=str)
        else:
            self.seq = np.array(list(seq))

        self.energyp = None
        self.energyn = None

        if self.mat == 'K':
            # Automatiically Fill in the K Matrix
            self.K = []
            dist = len(self.seq)
            for x in range(dist):
                if self.type == 'dna':
                    ibase = rnad[seq[x]]
                elif self.type == 'pep':
                    ibase = aad[seq[x]]
                for y in range(x + 1, dist):
                    if self.type == 'dna':
                        jbase = rnad[seq[y]]
                    elif self.type == 'pep':
                        jbase = aad[seq[y]]
                    if y <= dist - 2:
                        for z in range(x + 2, dist):
                            if self.type == 'dna':
                                kbase = rnad[seq[z]]
                            elif self.type == 'pep':
                                kbase = aad[seq[z]]
                            self.K.append((x, y-1, z-2, ibase, jbase, kbase))
                            if cv == 'yes' and self.type == 'dna':
                                cibase = rnad[base_flip_rna[seq[x]]]
                                cjbase = rnad[base_flip_rna[seq[y]]]
                                ckbase = rnad[base_flip_rna[seq[z]]]
                                self.K.append((x, y-1, z-2, cibase, cjbase, ckbase))

        elif self.mat == 'J':
            #Automatically Fill in the J Matrix
            self.J = []
            dist = len(self.seq)
            for x in range(dist):
                ibase = aad[seq[x]]
                for y in range(x + 1, dist):
                    jbase = aad[seq[y]]
                    self.J.append((x, y - 1, ibase, jbase))

        elif self.mat == 'H':
            #Automatically Fill in the H Matrix
            self.H = []
            dist = len(self.seq)
            for x in range(dist):
                ibase = aad[seq[x]]
                self.H.append((x, ibase))


    def score(self, ScoringMatrixPos, ScoringMatrixNeg):
        if self.mat == 'K':
            Pos = [x for x in self.K if x in ScoringMatrixPos]
            Neg = [x for x in self.K if x in ScoringMatrixNeg]
        elif self.mat == 'J':
            Pos = [x for x in self.J if x in ScoringMatrixPos]
            Neg = [x for x in self.J if x in ScoringMatrixNeg]
        elif self.mat == 'H':
            Pos = [x for x in self.H if x in ScoringMatrixPos]
            Neg = [x for x in self.H if x in ScoringMatrixNeg]
        return len(Pos)*1, len(Neg)*-1


def dJ_Worker(id, N, q, bi, ei, SEQHANDLER):
    ScoringMatrixPos = []
    ScoringMatrixNeg = []
    cS = -1
    for i in range(bi, ei):
        for j in range(N-1):
            for k in range(1, q):
                for l in range(1, q):
                    if i > j:
                        continue
                    if i == math.floor((ei-bi)/2 + bi) and j == i+1:
                        print('50% on Worker ' + str(id))
                    ScoringMatrixPos.append((i, j, k, l))
                    ScoringMatrixNeg.append((i, j, k, l))
                    Sscorepos, Sscoreneg = SEP_OPTIMIZED(SEQHANDLER, ScoringMatrixPos, ScoringMatrixNeg)
                    if Sscoreneg < cS and Sscorepos < cS:
                        ScoringMatrixNeg.remove((i, j, k, l))
                        ScoringMatrixPos.remove((i, j, k, l))
                    if Sscorepos > Sscoreneg and Sscorepos > cS:
                        cS = Sscorepos
                        ScoringMatrixNeg.remove((i, j, k, l))
                    if Sscoreneg > Sscorepos and Sscoreneg > cS:
                        cS = Sscoreneg
                        ScoringMatrixPos.remove((i, j, k,l))
    print('Worker ' + str(id))
    print('Optimal R Score Obtained: ' + str(cS))
    return ScoringMatrixPos, ScoringMatrixNeg, cS


def SEP_OPTIMIZED(SEQHANDLER, ScoringMatrixPos, ScoringMatrixNeg):
    pEg, pEb = [], []
    nEg, nEb = [], []
    for x in SEQHANDLER:
        if x.affinity == 1:
            pos, neg = x.score(ScoringMatrixPos, ScoringMatrixNeg)
            pEb.append(pos)
            nEb.append(neg)
        else:
            pos, neg = x.score(ScoringMatrixPos, ScoringMatrixNeg)
            pEg.append(pos)
            nEg.append(neg)
    ntb, nbb = max(nEb), min(nEb)
    ptb, pbb = max(pEb), min(pEb)
    ntg, nbg = max(nEg), min(nEg)
    ptg, pbg = max(pEg), min(pEg)
    sepposscore = (ptg - ptb) + 2*(pbg - ptb)
    sepnegscore = (ntg - ntb) + 2*(nbg - ntb)
    return sepposscore, sepnegscore


def designerH(N, q, fastafile, **kwargs):
    type = 'dna'
    for key, value in kwargs.items():
        if key == 'type':
            type = value
    finalH = np.full((N, q), 0.0)
    J = np.full((N-1, N-1, q, q), 0.0)
    titles, seqs = Fasta_Read_Aff(fastafile)
    currentR = -1
    if type == 'pep':
        for i in range(N):
            for j in range(1, q):
                finalH[i, j] = 1
                rscore = R_SCORE(titles, seqs, finalH, J, type='pep')
                if rscore > currentR:
                    currentR = rscore
                else:
                    finalH[i, j] = -1
                    rscore = R_SCORE(titles, seqs, finalH, J, type='pep')
                    if rscore > currentR:
                        currentR = rscore
                    else:
                        finalH[i, j] = 0.0
    else:
        for i in range(N):
            for j in range(1, q):
                finalH[i, j] = 1
                rscore = R_SCORE(titles, seqs, finalH, J)
                if rscore > currentR:
                    currentR = rscore
                else:
                    finalH[i, j] = -1
                    rscore = R_SCORE(titles, seqs, finalH, J)
                    if rscore > currentR:
                        currentR = rscore
                    else:
                       finalH[i, j] = 0.0
    print('designed H Matrix R score: ' + str(currentR))
    return finalH


def designerK(N, q, fastafile):
    finalK = np.full((N-2, N-2, N-2, q, q, q), 0.0)
    J = np.full((N - 1, N - 1, q, q), 0.0)
    H = np.full((N, q), 0.0)
    titles, seqs = Fasta_Read_Aff(fastafile)
    currentR = -1
    for i in range(N-2):
        for j in range(N - 2):
            for k in range(N - 2):
                for x in range(1, q):
                    for y in range(1, q):
                        for z in range(1, q):
                            if i > j or j > k:
                                continue
                            finalK[i, j, k, x, y, z] = 1
                            rscore = R_SCORE3(titles, seqs, H, J, finalK)
                            if rscore > currentR:
                                currentR = rscore
                            else:
                                finalK[i, j, k, x, y, z] = -1
                                rscore = R_SCORE3(titles, seqs, H, J, finalK)
                                if rscore > currentR:
                                    currentR = rscore
                                else:
                                    finalK[i, j, k, x, y, z] = 0.0
    print('designed K Matrix R score: ' + str(currentR))
    return finalK


def export_K(K, N, q, outfile):
    o = open(outfile, 'w')
    for i in range(N-2):
        for j in range(N - 2):
            for k in range(N - 2):
                for x in range(1, q):
                    for y in range(1, q):
                        for z in range(1, q):
                            if i > j or j > k:
                                continue
                            if K[i, j, k, x, y, z] != 0.0:
                                content = str(i)+','+str(j+1)+','+str(k+2)+','+str(x)+','+str(y)\
                                          +','+str(z)+','+str(K[i, j, k, x, y, z])
                                print(content, file=o)
    print('K written')
    o.close()

def export_J(J, N, q, outfile):
    o = open(outfile, 'w')
    c = ','
    for i in range(N-1):
        for j in range(N-1):
            for k in range(1, q):
                for l in range(1, q):
                    if J[i, j, k, l] != 0.0:
                        content = str(i) + c + str(j) + c + str(k) + c + str(l) + c + str(J[i, j, k, l])
                        print(content, file=o)
    o.close()


def export_H(H, N, q, outfile):
    o = open(outfile, 'w')
    c = ','
    for i in range(N):
        for j in range(1, q):
            if H[i, j] != 0.0:
                content = str(i) + c + str(j) + c + str(H[i, j])
                print(content, file=o)
    o.close()



def R_SCORE(titles, seqs, H, J, **kwargs):
    type = 'dna'
    for key, value in kwargs.items():
        if key == 'type':
            type = value
    energies = []
    for xid, x in enumerate(seqs):
        if xid == 0:
            N = len(list(x))
        if type == 'pep':
            energies.append(Calc_Energy_TCR(x, J, H, N))
        else:
            energies.append(Calc_Energy(x, J, H))
    api = list(zip(titles, energies))
    affs = list(set(titles))
    datax = []
    datae = []
    for x in affs:
        prospects = [nrg for (aff, nrg) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    rscore = linreg[2]
    return rscore


def R_SCORE3(titles, seqs, H, J, K):
    energies = []
    for x in seqs:
        energies.append(Calc_Energy3(x, J, H, K))
    api = list(zip(titles, energies))
    affs = list(set(titles))
    datax = []
    datae = []
    for x in affs:
        prospects = [nrg for (aff, nrg) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    rscore = linreg[2]
    return rscore


def R_SCORE_w_SEQHANDLER(SEQHANDLER, ScoringMatrix, titles):
    for x in SEQHANDLER:
        x.score(ScoringMatrix)
    affs = list(set(titles))
    datax = []
    datae = []
    for a in affs:
        maxprop = max([x.energy for x in SEQHANDLER if x.affinity == a])
        datax.append(a)
        datae.append(maxprop)
    linreg = stats.linregress(datax, datae)
    rscore = linreg[2]
    return rscore



def Raw_Aff_v_E(J, H, title, outpath, infile):
    titles, seqs = Fasta_Read_Aff(infile)
    energies = []
    for x in seqs:
        energies.append(Calc_Energy(x, J, H))
    plt.scatter(titles, energies, color='r', s=0.5)
    plt.title(title)
    plt.ylabel('Energy')
    plt.xlabel('Affinity')
    plt.savefig(outpath, dpi=600)


def Stats_Energy_Seqs(infile, J, H, N):
    seqs = Fasta_Read_SeqOnly(infile)
    energies = []
    for x in seqs:
        energies.append(Calc_Energy_TCR(x, J, H, N))
    return mean(energies), np.std(energies), min(energies), max(energies)


def write_stat_file_TCR(filename, clusters, nameL, means, stds, mins, maxs):
    o = open(filename, 'w')
    for x in range(len(clusters)):
        print('Cluster:', clusters[x],file=o)
        for y in range(len(nameL)):
            print('Fam:',nameL[y],'Mean:',means[x,y],'Std:',stds[x,y],'Min',mins[x,y],'Max',maxs[x,y],file=o)
    o.close()


def write_fasta_aff(seqs, affs, out):
    o = open(out, 'w')
    for xid, x in enumerate(affs):
        print('>seq' + str(xid) + '-' + str(x), file=o)
        print(seqs[xid], file=o)
    o.close()

def write_fasta_seqonly(seqs, out):
    o = open(out, 'w')
    for xid, x in enumerate(seqs):
        print('>seq' + str(xid), file=o)
        print(x, file=o)
    o.close()

def Raw_wRscore(J, H, outpath, infile):
    titles, seqs = Fasta_Read_Aff(infile)
    energies = []
    for x in seqs:
        energies.append(Calc_Energy(x, J, H))
    datax = []
    datae = []
    affs = list(set(titles))
    api = list(zip(titles, energies))
    highestaff = 1
    for x in affs:
        if x > highestaff: highestaff = x
        prospects = [nrg for (aff, nrg) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    xl = np.linspace(0, highestaff, 100)
    plt.plot(xl, xl * linreg[0] + linreg[1], ':r')
    cutoff = max([y for x, y in api if x == 1])
    plt.plot(xl, [cutoff for i in xl], ':b')
    plt.scatter(titles, energies, color='r', s=0.5)
    plt.ylabel('Energy')
    plt.xlabel('R Score: ' + str(linreg[2]))
    plt.savefig(outpath, dpi=600)
    print("Cutoff",cutoff)
    plt.close()


def Diff_Avg(J, H, outpath, labels, *infile, **kwargs):
    colors = ['g', 'm', 'c', 'b', 'y']
    title = 'Affinity vs Energy'
    for key, value in kwargs.items():
        if key == 'title':
            title = value
    oHa = 0
    oAff = set()
    adata = []
    apis = []
    data = []
    for i in infile:
        titles, seqs = Fasta_Read_Aff(i)
        energies = []
        for x in seqs:
            energies.append(Calc_Energy(x, J, H))
        affs = list(set(titles))
        adata.append(affs)
        oAff.update(affs)
        print(oAff)
        api = list(zip(titles, energies))
        apis += api
        data.append(api)
    for aid, ai in enumerate(data):
        avg = []
        err = []
        thaff = 0
        for aff in adata[aid]:
            if aff > oHa: oHa = aff
            if aff > thaff: thaff = aff
            yvals = np.array([y for (x, y) in ai if x == aff])
            yavg = yvals.mean()
            yerr = np.std(yvals)
            avg.append(yavg)
            err.append(yerr)
        xd = [x + .2*aid for x in adata[aid]]
        plt.errorbar(xd, avg, err, linestyle='None', marker='^', label=labels[aid])
        # x = list(set([x for (x, y) in data[aid]]))
        # x.sort()
        # linreg = stats.linregress(x, avg)
        # xl = np.linspace(0, thaff, 100)
        # plt.plot(xl, xl * linreg[0] + linreg[1], c=colors[aid])
    plt.xlabel('Avg Diff Plot')
    plt.ylabel('Energy')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0.)
    plt.suptitle(title)
    plt.savefig(outpath, dpi=600)


def Diff_Raw_wRscore(J, H, outpath, labels, *infile):
    colors = ['g', 'm', 'c', 'b', 'y']
    oHa = 0
    oAff = set()
    apis = []
    data = []
    for i in infile:
        titles, seqs = Fasta_Read_Aff(i)
        energies = []
        for x in seqs:
            energies.append(Calc_Energy(x, J, H))
        affs = list(set(titles))
        oAff.update(affs)
        print(oAff)
        api = list(zip(titles, energies))
        apis += api
        data.append(api)
    datax = []
    datae = []
    for x in oAff:
        if x > oHa: oHa = x
        prospects = [nrg for (aff, nrg) in apis if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    xl = np.linspace(0, oHa, 100)
    plt.plot(xl, xl * linreg[0] + linreg[1], ':r')
    cutoff = max([y for x, y in apis if x == 1])
    plt.plot(xl, [cutoff for i in xl], ':b')
    for xid, i in enumerate(data):
        titles, energies = zip(*i)
        plt.scatter(titles, energies, c=colors[xid], label=labels[xid])
    plt.ylabel('Energy')
    plt.xlabel('R Score: ' + str(linreg[2]))
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
    plt.savefig(outpath, dpi=600)
    print("Cutoff",cutoff)
    plt.close()


def Raw_wRscore_subplot(subplot, J, H, infile):
    titles, seqs = Fasta_Read_Aff(infile)
    energies = []
    for x in seqs:
        energies.append(Calc_Energy(x, J, H))
    datax = []
    datae = []
    affs = list(set(titles))
    api = list(zip(titles, energies))
    highestaff = 1
    for x in affs:
        if x > highestaff: highestaff = x
        prospects = [nrg for (aff, nrg) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    xl = np.linspace(0, highestaff, 100)
    subplot.plot(xl, xl * linreg[0] + linreg[1], ':r')
    #cutoff = max([y for x, y in api if x == 1])
    #subplot.plot(xl, [cutoff for i in xl], ':b')
    subplot.scatter(titles, energies)
    subplot.set_ylabel('Energy')
    subplot.set_xlabel('R Score: ' + str(linreg[2]))
    #print("Cutoff", cutoff)



def Raw_wRscore_TCR(J, H, outpath, infile):
    titles, seqs = Fasta_Read_Aff(infile)
    energies = []
    for xid, x in enumerate(seqs):
        if xid == 0:
            N = len(list(x))
        energies.append(Calc_Energy_TCR(x, J, H, N))
    datax = []
    datae = []
    affs = list(set(titles))
    api = list(zip(titles, energies))
    highestaff = 1
    for x in affs:
        if x > highestaff: highestaff = x
        prospects = [nrg for (aff, nrg) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    xl = np.linspace(0, highestaff, 100)
    plt.plot(xl, xl * linreg[0] + linreg[1], ':r')
    plt.scatter(titles, energies)
    plt.ylabel('Energy')
    plt.xlabel('R Score: ' + str(linreg[2]))
    plt.savefig(outpath, dpi=600)


def gen_badseq_mutt(J, H, JMutt, N, numberofnorms, **kwargs):
    ns = 'J'
    for key, value in kwargs.items():
        if key == 'normsource':
            ns = value
        else:
            print('No keyword argument ' + key + ' found')
    # Get Indices of top 10 norms
    bseq = np.full(40, ['X'], dtype=str)
    if ns == 'J':
        tval = TopX_JNorms(J, N, numberofnorms)
    if ns == 'mutt':
        tval = TopX_JNorms(JMutt, N, numberofnorms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_badseq(JMutt, pvals, x, y)
        bseq = Seq_edit_past_entry_comp(vals, bseq)
    for xid, x in enumerate(bseq):
        if x == 'X':
            gpx = np.argmin(H[xid, 1:5])
            bseq[xid] = rna[int(gpx) + 1]
    print(pvals)
    print(''.join(gseq))
    print(Calc_Energy(gseq, J, H))
    return ''.join(gseq)


# Input J and number of norms
def gen_goodseq_mutt(J, H, JMutt, N, numberofnorms, **kwargs):
    ns = 'J'
    for key, value in kwargs.items():
        if key == 'normsource':
            ns = value
        else:
            print('No keyword argument ' + key + ' found')
    # Get Indices of top 10 norms
    gseq = np.full(40, ['X'], dtype=str)
    if ns == 'J':
        tval = TopX_JNorms(J, N, numberofnorms)
    if ns == 'mutt':
        tval = TopX_JNorms(JMutt, N, numberofnorms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_goodseq(JMutt, pvals, x, y)
        gseq = Seq_edit_past_entry_comp(vals, gseq)
    for xid, x in enumerate(gseq):
        if x == 'X':
            gpx = int(np.where(H[xid, :] == np.amax(H[xid, :]))[0])
            gseq[xid] = rna[int(gpx)]
    print(pvals)
    print(''.join(gseq))
    print(Calc_Energy(gseq, J, H))
    return ''.join(gseq)


# Monte Carlo sampling for better binder
class GenerSeq:
    def __init__(self, N, T, mut_steps=5, out_after=100, steps=10000):
        self._history = []
        self._T = T
        self._beta = 1. / T
        self._mut_steps = mut_steps
        self._out_after = out_after
        self._steps = steps
        self._seq = np.random.choice(nucs, N)

    def calculate_energy(self, J, h):
        J_energy = 0.
        h_energy = 0.
        for i in range(0, len(self._seq)):
            t1 = nuc_to_id[self._seq[i]]
            h_energy += h[i, t1]
            for j in range(i + 1, len(self._seq)):
                t2 = nuc_to_id[self._seq[j]]
                J_energy += J[i, j - 1, t1, t2]
        return J_energy + h_energy

    def mutate_seq(self):
        for m in range(self._mut_steps):
            pos = np.random.choice(range(len(self._seq)))
            self._seq[pos] = np.random.choice(nucs)
        return self._seq

    def run_adaptive_sampling(self, J, h, outpath):
        total, acc = 0, 0
        out = open(outpath, 'w')
        oldene = self.calculate_energy(J, h)
        oldseq = copy.deepcopy(self._seq)
        for i in range(self._steps):
            total += 1
            self.mutate_seq()
            newene = self.calculate_energy(J, h)
            p = math.exp(self._beta * (newene - oldene))
            if .3 < acc/total < .8:
                continue
            elif acc/total < 0.3:
                self._T *= 1.05
            elif acc/total > .8:
                self._T *= 0.95
            if random.random() < p:
                acc+=1
                # accept move
                oldene = newene
                oldseq = copy.deepcopy(self._seq)
                if i % self._out_after:
                    if ''.join(self._seq) not in self._history:
                        self._history.append(''.join(self._seq))
                        print('>' + str(i) + '-' + str(newene), file=out)
                        print(''.join(self._seq), file=out)
                    print(str((i / self._steps) * 100) + ' Percent Done')
            else:
                self._seq = copy.deepcopy(oldseq)
        out.close()

    def run_adaptive_sampling_wT(self, J, h, outpath, outTempPlot):
        total, acc = 0, 0
        xdata, tdata = [], []
        out = open(outpath, 'w')
        oldene = self.calculate_energy(J, h)
        oldseq = copy.deepcopy(self._seq)
        for i in range(self._steps):
            total += 1
            self.mutate_seq()
            newene = self.calculate_energy(J, h)
            p = math.exp(self._beta * (newene - oldene))
            if .3 < acc / total < .8:
                continue
            elif acc / total < 0.3:
                self._T *= 1.05
            elif acc / total > .8:
                self._T *= 0.95
            if random.random() < p:
                acc += 1
                # accept move
                oldene = newene
                oldseq = copy.deepcopy(self._seq)
                if i % self._out_after:
                    if ''.join(self._seq) not in self._history:
                        xdata.append(i)
                        tdata.append(self._T)
                        self._history.append(''.join(self._seq))
                        print('>' + str(i) + '-' + str(newene), file=out)
                        print(''.join(self._seq), file=out)
                    print(str((i / self._steps) * 100) + ' Percent Done')
            else:
                self._seq = copy.deepcopy(oldseq)
        plt.plot(xdata, tdata, c='r')
        plt.savefig(outTempPlot, dpi=400)
        out.close()

    def run_sampling(self, J, h, outpath):
        total, acc = 0, 0
        out = open(outpath, 'w')
        oldene = self.calculate_energy(J, h)
        oldseq = copy.deepcopy(self._seq)
        for i in range(self._steps):
            total += 1
            self.mutate_seq()
            newene = self.calculate_energy(J, h)
            p = math.exp(self._beta * (newene - oldene))
            # if .25 < acc/total < .35:
            #     continue
            # elif acc/total < 0.25:
            #     self._T *= 1.05
            # elif acc/total > .35:
            #     self._T *= 0.95
            if random.random() < p:
                acc+=1
                # accept move
                oldene = newene
                oldseq = copy.deepcopy(self._seq)
                if i % self._out_after:
                    if ''.join(self._seq) not in self._history:
                        self._history.append(''.join(self._seq))
                        print('>' + str(i) + '-' + str(newene), file=out)
                        print(''.join(self._seq), file=out)
                    print(str((i / self._steps) * 100) + ' Percent Done')
            else:
                self._seq = copy.deepcopy(oldseq)
        out.close()

    def run_bad_sampling(self, J, h, outpath):
        out = open(outpath, 'w')
        oldene = self.calculate_energy(J, h)
        oldseq = copy.deepcopy(self._seq)
        for i in range(self._steps):
            self.mutate_seq()
            newene = self.calculate_energy(J, h)
            p = math.exp(self._beta * (newene - oldene))
            if random.random() > p:
                # accept move
                oldene = newene
                oldseq = copy.deepcopy(self._seq)
                if i % self._out_after:
                    if ''.join(self._seq) not in self._history:
                        self._history.append(''.join(self._seq))
                        print('>' + str(i) + '-' + str(newene), file=out)
                        print(''.join(self._seq), file=out)
                    print((str((i / self._steps) * 100)) + ' Percent Done')
            else:
                self._seq = copy.deepcopy(oldseq)
        out.close()


# Monte Carlo sampling for Secondary Structure Prediction from plmDCA output


def import_correlations(file, N):
    o = open(file, 'r')
    corr = np.full((N-1, N-1), 0.0)
    for line in o:
        data = line.split(',')
        corr[int(data[0])-1, int(data[1])-2] = float(data[2].rstrip())
    o.close()
    return corr

# RULES each base can only bind to one other base
# Not all bases need by bound
def sample_corr(corr, N):
    tcorr = copy.deepcopy(corr).flatten()
    xs, ys = [], []
    for i in range(50):
        tmp = np.argmax(tcorr)
        x = math.floor(tmp/(N-1))
        y = tmp%(N-1)
        tcorr[tmp]=0.0
        if abs(x-y) > 3:
            xs.append(x)
            ys.append(y)
    for i in range(len(xs)):
        print(xs[i]+1,ys[i]+2)


def mc_analysis_plot(mcfile, ensemblefile, N, J, H, out):
    seqs = Fasta_Read_SeqOnly(mcfile)
    energy, mutations = [], []
    for seq in seqs:
        escore = ensemble_checker(ensemblefile, seq)[0][1]
        mut_num = N-escore*N
        E = Calc_Energy(seq, J, H)
        energy.append(E)
        mutations.append(mut_num)
    fig, ax = plt.subplots(1,2)
    x = np.arange(0, len(seqs), 1)
    ax[0].plot(x, energy, color='g')
    ax[0].set_ylabel('Energy')
    ax[0].set_xlabel('Seq')
    ax[1].plot(x, mutations, color='r')
    ax[1].set_ylabel('Mutation #')
    ax[1].set_xlabel('Seq')
    plt.savefig(out, dpi=400)
    plt.close()


def find_closest_seqs(seq, ensemblefile, N):
    affs, seqs = Fasta_Read_Aff(ensemblefile)
    sl = list(rna2dna(seq))
    bs, gs = 0, 0
    for sid, es in enumerate(seqs):
        el = list(es)
        tmp = 0
        if affs[sid] == 1:
            for i in range(N):
                if el[i] == sl[i]:
                    tmp += 1
            if tmp > bs:
                cbs = es
                bs = tmp
        elif affs[sid] > 1:
            for i in range(N):
                if el[i] == sl[i]:
                    tmp += 1
            if tmp > gs:
                ga = affs[sid]
                cgs = es
                gs = tmp
        else:
            print(affs[sid])
    gd, bd = dis_seqs(sl, cgs), dis_seqs(sl, cbs)
    return cbs, cgs, bd, gd, ga


def avg_mut_energy(seq, mutnum, seqnum, N, J, H):
    sseq = list(seq)
    es = []
    for x in range(seqnum):
        tseq = copy.deepcopy(sseq)
        for i in range(mutnum):
            pos = np.random.choice(range(len(tseq)))
            tseq[pos] = np.random.choice(nucd)
        es.append(Calc_Energy(tseq, J, H))
    return mean(es), np.std(es)





def analyze_MC_seqs(cutoff, ensemblefile, N, J, H, epath, mutpath, *outfiles):
    eaffs, eseqs = Fasta_Read_Aff_wE_above80(ensemblefile, 80.0, J, H)
    mseqs, mes, info = [], [], []
    fseqs, fes = [], []
    ef = open(epath, 'w+')
    mf = open(mutpath, 'w+')
    for oid, out in enumerate(outfiles):
        tes, tseqs = Fasta_Read_MC(out)
        for xid, x in enumerate(tes):
            if x > cutoff:
                mes.append(tes[xid])
                mseqs.append(tseqs[xid])
        if oid > 0:
            mes, mseqs = prune_alignment(mseqs, 1.0, mes)
    removed = 0
    for sid, seq in enumerate(mseqs):
        if sid % 10 == 0:
            print('seq', sid)
        present = [x for x in eseqs if x == seq]
        if not present:
            fseqs.append(seq)
            fes.append(mes[sid])
        else:
            removed += 1
    print('Removed', removed, 'Sequences')
    for sid, seq in enumerate(fseqs):
        sl = list(seq)
        rm = False
        bs, gs = 0, 0
        cbs, gbs = "A", "A"
        for esid, es in enumerate(eseqs):
            el = list(es)
            gtmp, btmp = 0, 0
            if eaffs[esid] == 1:
                for i in range(N):
                    if el[i] == sl[i]:
                        btmp += 1
                if btmp == N:
                    rm = True
                if btmp > bs:
                    cbs = es
                    bs = btmp
            elif eaffs[esid] > 1:
                for i in range(N):
                    if el[i] == sl[i]:
                        gtmp += 1
                if gtmp == N:
                    rm = True
                if gtmp > gs:
                    cgs = es
                    gs = gtmp
        print('seq', sid, flush=True)
        if rm == True:
            continue
        else:
            if gs > bs:
                muts = N - gs
                info.append((seq, fes[sid], muts, '*' + cgs, cbs))
            elif bs > gs:
                muts = N - bs
                info.append((seq, fes[sid], muts, cgs, '*' + cbs))
            else:
                muts = N - gs
                info.append((seq, fes[sid], muts, '*' + cgs, cbs))
    eorder = copy.deepcopy(info)
    eorder.sort(key=lambda tup: tup[1])
    info.sort(key=lambda tup: tup[2])
    for s, e, m, gs, bs in eorder:
        print(s, e, m, gs, bs, file=ef)
    for s, e, m, gs, bs in info:
        print(s, e, m, gs, bs, file=mf)
    ef.close()
    mf.close()


def dis_seqs(seq1, seq2):
    d = 0
    u = zip(seq1, seq2)
    for i, j in u:
        if i != j:
            d += 1
    return d

def diffs_seqs(seq1, seq2, affinity=0):
    diff = []
    u = zip(seq1, seq2)
    for pos, val in enumerate(u):
        i, j = val
        if i != j:
            if affinity == 0:
                diff.append((pos+1, j))
            else:
                char = str(affinity) + j
                diff.append((pos + 1, char))
    return diff

def sort_diffs(totaldiffs, N):
    adiff = ['' for x in range(N)]
    td = []
    for diff in totaldiffs:
        for x, m in diff:
            adiff[x-1] += m
    for aid, a in enumerate(adiff):
        td.append((aid+1, a))
    return td


def local_env(seq, ensemblefile, basesaway, N):
    eaffs, eseqs = Fasta_Read_Aff(ensemblefile)
    blseqs, glseqs, tgd, tbd = [], [], [], []
    g1, b1, g2, b2, g3, b3 = [], [], [], [], [], []
    bad_muts, good_muts, neut = [], [], []
    for eid, es in enumerate(eseqs):
        d = dis_seqs(seq, es)
        aff = eaffs[eid]
        if d > basesaway:
            continue
        else:
            if d ==1:
                if eaffs[eid] == 1:
                    b1.append(es)
                else:
                    g1.append((aff, es))
            if d == 2:
                if eaffs[eid] == 1:
                    b2.append(es)
                else:
                    g2.append((aff, es))
            if d == 3:
                if eaffs[eid] == 1:
                    b3.append(es)
                else:
                    g3.append(es)
    for b in b1:
        bad_muts.append(diffs_seqs(seq, b))
    for aff, g in g1:
        good_muts.append(diffs_seqs(seq, g, aff))
    for b in b2:
        d = diffs_seqs(seq, b)
        comm = [x for x in d if x in bad_muts]
        if comm:
            d.remove(comm[0])
            neut.append(d)
        else:
            tbd.append(d)
    for aff, g in g2:
        d = diffs_seqs(seq, g, aff)
        comm = [x for x in d if x in good_muts]
        if comm:
            d.remove(comm[0])
            neut.append(d)
        else:
            tgd.append(d)
    gd1 = sort_diffs(good_muts, N)
    bd1 = sort_diffs(bad_muts, N)
    gd2 = sort_diffs(tgd, N)
    bd2 = sort_diffs(tbd, N)
    return gd1, gd2, bd1, bd2


def avgdis_ensemble(seq, ensemblefile):
    eaffs, eseqs = Fasta_Read_Aff(ensemblefile)
    td = 0
    for es in eseqs:
        td += dis_seqs(seq, es)
    return td / len(eseqs)


def avgdis_gb(seq, ensemblefile):
    gaffs, gseqs = Fasta_Read_GB(ensemblefile)
    td = 0
    for gs in gseqs:
        td += dis_seqs(seq, gs)
    return td / len(gseqs)


def avgdis_bb(seq, ensemblefile):
    baffs, bseqs = Fasta_Read_BB(ensemblefile)
    td = 0
    for bs in bseqs:
        td += dis_seqs(seq, bs)
    return td/len(bseqs)

def rand_mut(seq, numofmuts, J, H):
    tseq = list(seq)
    for x in range(numofmuts):
        pos = np.random.choice(range(len(tseq)))
        tseq[pos] = np.random.choice(nucd)
    E = Calc_Energy(''.join(tseq), J, H)
    return ''.join(tseq), E

def all_mut_possibilities(seql):
    allpos = []
    for iid, i in enumerate(seql):
        alternatives = [x for x in nucd if x != i]
        for an in alternatives:
            seql[iid] = an
            allpos.append(copy.deepcopy(seql))
        seql[iid] = i
    return allpos

def ThreeMutation_checker(seq, J, H):
    tseq = list(copy.deepcopy(seq))
    results = []
    alts = all_mut_possibilities(tseq)
    alt2 = []
    for x in alts:
        alt2 += all_mut_possibilities(x)
    alt3 = []
    for x in alt2:
        alt3 += all_mut_possibilities(x)
    for x in alt3:
        if dis_seqs(x, seq) == 3:
            results.append((''.join(x), Calc_Energy(x, J, H)))
    results.sort(key=lambda tup: tup[1])
    highest = results[:5]
    lowest = results[-5:]
    return highest, lowest

def score_dist(ensemblefile, J, H, path):
    eaffs, eseqs = Fasta_Read_Aff(ensemblefile)
    affs = set(eaffs)
    pas, pes, alle = [], [], []
    for a in affs:
        es = []
        for xid, x in enumerate(eseqs):
            if eaffs[xid] == a:
                es.append(Calc_Energy(x, J, H))
        alle += es
        pas.append(a)
        pes.append(mean(es))
    fig, ax = plt.subplots()
    deN = gaussian_kde(alle)
    xd1 = np.linspace(min(alle), max(alle), 100)
    ax.plot(xd1, deN(xd1), color='r')
    ax.plot(alle, [0.00] * len(alle), '|', color='k')
    cs = ['b', 'g', 'c', 'm', 'y']
    for i in range(len(pas)):
        if i == 0:
            ax.axvline(x=pes[i], c='k', lw=0.5)
        else:
            ax.axvline(x=pes[i], c=cs[i % 5], lw=0.5)
        plt.text(pes[i] + 0.1, 0.002 + (i % 10 * 0.005), str(pas[i]), rotation=90, c=cs[i % 5])
        plt.text(20, .054 - i*0.002, str(pas[i]) + ': ' + str(round(pes[i], 2)), c=cs[i % 5])
    ax.set_xlabel('Energy')
    ax.grid(True)
    ax.title.set_text('Score Distribution')
    ax.title.set_size(fontsize=10)
    plt.savefig(path, dpi=600)


def mut_loc(mseq, tseq):
    u = zip(mseq, tseq)
    d = []
    for xid, x in enumerate(u):
        i, j = x
        if i != j:
            d.append((xid+1, str(j)))
    return d

#def mutation_checker(seq, J, H, mutnum, seqfile):
    # tseq = list(copy.deepcopy(seq))
    # tmp, results = [], []
    # for x in range(mutnum):
    #     if x == 0:
    #         alts = all_mut_possibilities(tseq)
    #     else:
    #         for x in alts:
    #             tmp += all_mut_possibilities(x)
    #             if x != mutnum-1:
    #                 alts = tmp
    #             else:
    #                 for x in tmp:
    #                     if dis_seqs(x, seq) == 3:
    #                         results.append((''.join(x), Calc_Energy(x, J, H)))
    # results.sort(key=lambda tup: tup[1])
    # lowest, le = results[0][0], results[0][1]
    # highest, he = results[-1][0], results[-1][1]
    # la = ensemble_checker(seqfile, lowest)
    # ha = ensemble_checker(seqfile, highest)
    # return la, le, ha, he

