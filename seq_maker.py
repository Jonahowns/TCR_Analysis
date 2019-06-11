#!/usr/bin/env python
import numpy as np
import copy


droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
fullpath = upath + droppath
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U':4}
rna = ['-', 'A', 'C', 'G', 'U']


def sortjmat(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N-1, N-1, q, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0])-1, int(data[1])-2, int(data[2])-1, int(data[3])-1] = float(data[4].rstrip())
    o.close()
    return fullmatrix


def sorthmat(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]) - 1, int(data[1]) - 1] = float(data[2].rstrip())
    o.close()
    return fullmatrix


def topxjnorms(J, N, x):
    jnorm = np.full((N-1, N-1), 0.0)
    vals = []
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
            if jnorm[i, j] != 0.0:
                vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
    vals.sort(key=lambda tup: tup[2])
    ind = int(-x)
    top10 = vals[ind:-1]
    print(ind, -1)
    print(vals)
    print(vals[ind:-1])
    return top10


def HJ_mutant_RNA(J, H, N):
    mutt = copy.deepcopy(J)
    for x in range(N - 1):
        for k in range(5):  # J Indices
            mutt[x, x:N, k, :] += H[x, k]
    for y in range(N - 1):
        for l in range(5):  # y states
            mutt[0:y + 1, y, :, l] += H[y + 1, l]
    return mutt


def HJ_mutant_Pep(J,H,N):
    mutt = copy.deepcopy(J)
    for x in range(N - 1):
        for k in range(21):  # J Indices
            mutt[x, x:N, k, :] += H[x, k]
    for y in range(N - 1):
        for l in range(21):  # y states
            mutt[0:y + 1, y, :, l] += H[y + 1, l]
    return mutt


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


def past_entry_comp_goodseq(J, pvals, xn, yn):
    tmppvals = copy.deepcopy(pvals)
    ind, xid, stype = check_vals(xn, yn, pvals)
    if ind == 0:
        tmpxn, tmpyn = list(np.where(J[xn, yn, :, :] == np.amax(J[xn, yn, :, :])))
        rxn = int(tmpxn)
        ryn = int(tmpyn)
        tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
        return [xn, yn, rxn, ryn], tmppvals
    elif ind == 1:
        xp, yp, rxp, ryp, val = tmppvals[xid]
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
        rxn = int(tmpxn)+1
        ryn = int(tmpyn)+1
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
                ryn = int(np.where(J[xn, yn, rxn, 1:5] == np.amin(J[xn, yn, rxn, 1:5]))[0])+1
            if stype == 'ys':
                ryn = ryp
                rxn = int(np.where(J[xn, yn, 1:5, ryn] == np.amin(J[xn, yn, 1:5, ryn]))[0])+1
            if stype == 'xnyp':
                rxn = ryp
                ryn = int(np.where(J[xn, yn, rxn, 1:5] == np.amin(J[xn, yn, rxn, 1:5]))[0])+1
            if stype == 'ynxp':
                ryn = rxp
                rxn = int(np.where(J[xn, yn, 1:5, ryn] == np.amin(J[xn, yn, 1:5, ryn]))[0])+1
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice < pchoice:
            if stype == 'xs':
                rxp = rxn
                ryp = int(np.where(J[xp, yp, rxp, 1:5] == np.amin(J[xp, yp, rxp, 1:5]))[0])+1
            if stype == 'ys':
                ryp = ryn
                rxp = int(np.where(J[xp, yp, 1:5, ryp] == np.amin(J[xp, yp, 1:5, ryp]))[0])+1
            if stype == 'xnyp':
                ryp = rxn
                rxp = int(np.where(J[xp, yp, 1:5, ryp] == np.amin(J[xp, yp, 1:5, ryp]))[0])+1
            if stype == 'ynxp':
                rxp = ryn
                ryp = int(np.where(J[xp, yp, rxp, 1:5] == np.amin(J[xp, yp, rxp, 1:5]))[0])+1
            tmppvals[xid] = (xp, yp, rxp, ryp, J[xp, yp, rxp, ryp])
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice == pchoice:
            tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        vals = [xp, yp, rxp, ryp, xn, yn, rxn, ryn]
        return vals, tmppvals


def Seq_edit_past_entry_comp(array, gseq):
    if len(array) >= 4:
        gseq[array[0]+0] = rna[int(array[2])]
        gseq[array[1]+1] = rna[int(array[3])]
        if len(array) == 8:
            gseq[array[4]+0] = rna[int(array[6])]
            gseq[array[5]+1] = rna[int(array[7])]
    return gseq


def gen_goodseq(famid, norms):
    analysispath = fullpath
    # Matrix Paths
    Jp = fullpath + str(famid) + 'j'
    Hp = fullpath + str(famid) + 'h'
    # N
    N = 40
    # Get Matrix Ready
    J = sortjmat(Jp, N, 5)
    H = sorthmat(Hp, N, 5)
    # Get Indices of top 10 norms
    gseq = np.full(40, ['X'], dtype=str)
    # bseq = np.full(40, ['X'], dtype=str)
    tval = topxjnorms(J, N, norms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_goodseq(J, pvals, x, y)
        gseq = Seq_edit_past_entry_comp(vals, gseq)
    for xid, x in enumerate(gseq):
        if x == 'X':
            gpx = int(np.where(H[xid, :] == np.amax(H[xid, :]))[0])
            gseq[xid] = rna[int(gpx)]
    print(''.join(gseq))
    return ''.join(gseq)


def gen_badseq(famid, norms):
    analysispath = fullpath
    # Matrix Paths
    Jp = fullpath + str(famid) + 'j'
    Hp = fullpath + str(famid) + 'h'
    # N
    N = 40
    # Get Matrix Ready
    J = sortjmat(Jp, N, 5)
    H = sorthmat(Hp, N, 5)
    # Get Indices of top 10 norms
    bseq = np.full(40, ['X'], dtype=str)
    tval = topxjnorms(J, N, norms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_badseq(J, pvals, x, y)
        bseq = Seq_edit_past_entry_comp(vals, bseq)
    for xid, x in enumerate(bseq):
        if x == 'X':
            gpx = int(np.where(H[xid, 1:5] == np.amin(H[xid, 1:5]))[0])
            bseq[xid] = rna[int(gpx)+1]
    print(pvals)
    print(''.join(bseq))
    return ''.join(bseq)


def gen_badseq_mutt(famid, norms):
    analysispath = fullpath
    # Matrix Paths
    Jp = fullpath + str(famid) + 'j'
    Hp = fullpath + str(famid) + 'h'
    # N
    N = 40
    # Get Matrix Ready
    J = sortjmat(Jp, N, 5)
    H = sorthmat(Hp, N, 5)
    jhmutt = HJ_mutant_RNA(J, H, N)
    # Get Indices of top 10 norms
    bseq = np.full(40, ['X'], dtype=str)
    tval = topxjnorms(J, N, norms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_badseq(jhmutt, pvals, x, y)
        bseq = Seq_edit_past_entry_comp(vals, bseq)
    for xid, x in enumerate(bseq):
        if x == 'X':
            gpx = int(np.where(H[xid, 1:5] == np.amin(H[xid, 1:5]))[0])
            bseq[xid] = rna[int(gpx)+1]
    print(pvals)
    print(''.join(bseq))
    return ''.join(bseq)


def gen_goodseq_mutt(famid, norms):
    analysispath = fullpath
    # Matrix Paths
    Jp = fullpath + str(famid) + 'j'
    Hp = fullpath + str(famid) + 'h'
    N = 40
    # Get Matrix Ready
    J = sortjmat(Jp, N, 5)
    H = sorthmat(Hp, N, 5)
    # Get Matrix Ready
    jhmutt = HJ_mutant_RNA(J, H, N)
    # Get Indices of top 10 norms
    gseq = np.full(40, ['X'], dtype=str)
    # bseq = np.full(40, ['X'], dtype=str)
    tval = topxjnorms(J, N, norms)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp_goodseq(jhmutt, pvals, x, y)
        gseq = Seq_edit_past_entry_comp(vals, gseq)
    for xid, x in enumerate(gseq):
        if x == 'X':
            gpx = int(np.where(H[xid, :] == np.amax(H[xid, :]))[0])
            gseq[xid] = rna[int(gpx)]
    print(pvals)
    print(''.join(gseq))
    return ''.join(gseq)


def Calc_Energy(seq, J, H):
    full = list(seq)
    dist = len(full)
    Jenergy = 0
    Henergy = 0
    for x in range(1, dist):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        for y in range(x+1, dist):
            jbase = rnad[seq[y]]
            Jenergy += J[x-1, y-2, ibase, jbase]
    energy = Jenergy + Henergy
    return energy


famid = 5
Jp = fullpath + str(famid) + 'j'
Hp = fullpath + str(famid) + 'h'
# N
N = 40

# Get Matrix Ready
J = sortjmat(Jp, N, 5)
H = sorthmat(Hp, N, 5)
'''
tseq ='AGGGGUUGGUGGGGUUGGAAAGGUGCUGGUUGGGACGGGG'
print(tseq)
best5 = gen_goodseq_mutt(famid, J, H, N, 10)
b5en = Calc_Energy(best5, J, H)
ten = Calc_Energy(tseq, J, H)
print(b5en)
print(ten)
'''
# bseq5 = gen_badseq(5, 250)
# bsen = Calc_Energy(bseq5, J, H)
tbseq = 'ACCAAUAUCACUCCCCGUUUUAAAUUUUAUACUAUAACAU'
tbben = Calc_Energy(tbseq, J, H)
bs5 = gen_badseq_mutt(5, 50)
b5bad = Calc_Energy(bs5, J, H)
print(b5bad)
print(tbben)

