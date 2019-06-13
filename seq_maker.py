#!/usr/bin/env python
import numpy as np
import copy
import numpy.linalg
import matplotlib.pyplot as plt

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


def gen_goodseq_custom_HJ(famid, norms, H, J):
    # N
    N = 40
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
            gpx = np.argmin(H[xid, 1:5])
            print(gpx)
            gseq[xid] = rna[int(gpx)+1]
    print(Calc_Energy(gseq, J, H))
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


def pot_energy(J, H, N): # This is Pretty Useless
    pEH = np.full((N, 5), 0.0)
    pEJ = np.full((N, 5), 0.0)
    for i in range(N):  # Fill H Values
        '''
        pEH[i, 0] += abs(H[i, 1])  # A
        pEH[i, 1] += abs(H[i, 2])  # C
        pEH[i, 2] += abs(H[i, 3])  # G
        pEH[i, 3] += abs(H[i, 4])  # T
        pEH[i, 4] += np.sum(pEH[i, 0:4])
        '''
        pEH[i, 0] += H[i, 1]  # A
        pEH[i, 1] += H[i, 2]  # C
        pEH[i, 2] += H[i, 3]  # G
        pEH[i, 3] += H[i, 4]  # T
        pEH[i, 4] += np.sum(pEH[i, 0:4])
    for i in range(N-1):  # Fill J Values
        for j in range(N-1):
            '''
            # X Contrib
            pEJ[i, 0] += np.sum(np.absolute(J[i, j, 1, :]))
            pEJ[i, 1] += np.sum(np.absolute(J[i, j, 2, :]))
            pEJ[i, 2] += np.sum(np.absolute(J[i, j, 3, :]))
            pEJ[i, 3] += np.sum(np.absolute(J[i, j, 4, :]))
            # Y Contrib
            pEJ[i+1, 0] += np.sum(np.absolute(J[j, i, :, 1]))
            pEJ[i+1, 1] += np.sum(np.absolute(J[j, i, :, 2]))
            pEJ[i+1, 2] += np.sum(np.absolute(J[j, i, :, 3]))
            pEJ[i+1, 3] += np.sum(np.absolute(J[j, i, :, 4]))
            '''
            pEJ[i, 0] += np.sum(J[i, j, 1, :])
            pEJ[i, 1] += np.sum(J[i, j, 2, :])
            pEJ[i, 2] += np.sum(J[i, j, 3, :])
            pEJ[i, 3] += np.sum(J[i, j, 4, :])
            # Y Contrib
            pEJ[i + 1, 0] += np.sum(J[j, i, :, 1])
            pEJ[i + 1, 1] += np.sum(J[j, i, :, 2])
            pEJ[i + 1, 2] += np.sum(J[j, i, :, 3])
            pEJ[i + 1, 3] += np.sum(J[j, i, :, 4])

        pEJ[i, 4] += np.sum(pEJ[i, 0:4])
    pEJ[39, 4] += np.sum(pEJ[39, 0:4])
    TPE = np.add(pEH, pEJ)
    t80 = np.percentile(TPE, 60, axis=0)
    seq = np.full(40, 'x', dtype=str)
    print(t80)
    for xid, row in enumerate(TPE):
        a, c, g, t, tot = row
        if a > t80[0] and c < t80[1] and g < t80[2] and t < t80[3]:
            print(xid+1, 'A')
        if a > t80[0] and c > t80[1] and g < t80[2] and t < t80[3]:
            print(xid+1, 'C')
        if a < t80[0] and c < t80[1] and g > t80[2] and t < t80[3]:
            print(xid+1, 'G')
        if a < t80[0] and c < t80[1] and g < t80[2] and t > t80[3]:
            print(xid+1, 'T')
        if a > t80[0] and c > t80[1] and g < t80[2] and t < t80[3]:
            print(xid+1, 'A or C')
        if a > t80[0] and c < t80[1] and g > t80[2] and t < t80[3]:
            print(xid+1, 'A or G')
        if a > t80[0] and c < t80[1] and g < t80[2] and t > t80[3]:
            print(xid+1, 'A or T')
        if a < t80[0] and c > t80[1] and g < t80[2] and t < t80[3]:
            print(xid+1, 'C or G')
        if a < t80[0] and c > t80[1] and g < t80[2] and t > t80[3]:
            print(xid+1, 'C or T')
        if a < t80[0] and c < t80[1] and g > t80[2] and t > t80[3]:
            print(xid+1, 'G or T')
        if a < t80[0] and c > t80[1] and g > t80[2] and t > t80[3]:
            print(xid+1, 'C G or T')
        if a > t80[0] and c < t80[1] and g > t80[2] and t > t80[3]:
            print(xid+1, 'A G or T')
        if a > t80[0] and c > t80[1] and g > t80[2] and t < t80[3]:
            print(xid+1, 'A C or G')
        if a > t80[0] and c > t80[1] and g < t80[2] and t > t80[3]:
            print(xid+1, 'A C or T')

    '''
    TPEn = np.divide(TPE[:, 4], nc)
    seqxvals = np.arange(1, 41, 1)
    viz = np.concatenate((seqxvals, TPEn))
    vizuL = np.reshape(viz, (2, 40))
    print(vizuL)
    '''   # T    #


def Rank_Test_seq(famid, outfile, J, H):
    o=open(fullpath + str(famid) + 'thgs.txt')
    titles = []
    seqs = []
    for line in o:
        if line.startswith('>'):
            titles.append(float(line.rstrip().split('-')[1]))
        else:
            seqs.append(line.rstrip())
    o.close()
    energies = []
    for x in seqs:
        nrg = Calc_Energy(x, J, H)
        energies.append(nrg)
    fnp = list(zip(energies, titles, seqs))
    op = open(fullpath + outfile, 'w')
    for line in fnp:
        print(line, file=op)
    op.close()
    return fnp


def Plot_seq_aff_v_E(famid, J, H):
    o=open(fullpath + str(famid) + 'thfull.txt')
    titles = []
    seqs = []
    for line in o:
        if line.startswith('>'):
            titles.append(float(line.rstrip().split('-')[1]))
        else:
            seqs.append(line.rstrip())
    o.close()
    energies = []
    for x in seqs:
        nrg = Calc_Energy(x, J, H)
        energies.append(nrg)
    api = list(zip(titles, energies))
    x = list(set([x for (x,y) in api]))
    x.sort()
    avg = []
    err = []
    for aff in x:
        yvals = np.array([y for (x, y) in api if x==aff])
        yavg = yvals.mean()
        yerr = np.std(yvals)
        avg.append(yavg)
        err.append(yerr)
    plt.errorbar(x, avg, err, linestyle='None', marker='^')
    plt.xlabel('affinity')
    plt.ylabel('Energy')
    plt.suptitle('Family: ' + str(famid) + ' Affinity vs Energy')
    plt.savefig(fullpath + str(famid) + 'affvsEnormscore.png', dpi=600)


def mixed_HJ(famid):
    analysispath = fullpath
    # Matrix Paths
    Jp = fullpath + str(famid) + 'j'
    Hp = fullpath + str(famid) + 'h'
    bJp = fullpath + str(famid) + 'bj'
    bHp = fullpath + str(famid) + 'bh'
    # N
    N = 40
    # Get Matrix Ready
    J = sortjmat(Jp, N, 5)
    H = sorthmat(Hp, N, 5)
    bJ = sortjmat(bJp, N, 5)
    bH = sorthmat(bHp, N, 5)
    topJ = topxjnorms(J, N, 80)
    topBJ = topxjnorms(bJ, N, 80)
    for xj, yj, val in topJ:
        for xb, yb, valb in topBJ:
            if xj == xb and yj == yb:
                J[xj, yj, :, :] = 0.0
    # J += 1
    # bJ += 1
    # J *= np.divide(1, np.sum(J))
    # bJ *= np.divide(1, np.sum(bJ))
    # J /= np.divide(J, bJ)
    # bJ /= np.divide(bJ, J)
    # print(np.sum(J))
    # print(np.sum(bJ))
    H += 1
    bH += 1
    H /= np.divide(H, bH)
    H *= np.divide(1, np.sum(H))
    bH *= np.divide(1, np.sum(bH))
    H /= np.divide(H, bH)
    bH /= np.divide(bH, H)
    print(np.sum(H))
    print(np.sum(bH))
    truH = H - bH
    return truH, J


famid = 8
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

# bseq5 = gen_badseq(5, 250)
# bsen = Calc_Energy(bseq5, J, H)
tbseq = 'ACCAAUAUCACUCCCCGUUUUAAAUUUUAUACUAUAACAU'
tbben = Calc_Energy(tbseq, J, H)
bs5 = gen_badseq_mutt(5, 50)
b5bad = Calc_Energy(bs5, J, H)
print(b5bad)
print(tbben)
'''




