

import numpy as np



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


def past_entry_comp(J, pvals, xn, yn):
    ind, xid, stype = check_vals(xn, yn, pvals)
    if ind == 0:
        rxn, ryn = list(np.where(J[xn, yn, :, :] == np.amax(J[xn, yn, :, :])))
        pvals.append(xn, yn, rxn, ryn, J[xn, yn, rxn, ryn])
        return [0], pvals
    if ind == 1:
        xp, yp, rxp, ryp, val = pvals[xid]
        val = np.amax(J[xn, yn, :, :])  # Highest in N w/ no constraints
        rxn, ryn = list(np.where(J[xn, yn, :, :] == val))  # Indices of highest in N
        if stype == 'xs':
            pchoice = J[xp, yp, rxp, ryp] + np.max(J[xn, yn, rxp, :])
            tchoice = np.max(J[xp, yp, rxn, :]) + J[xn, yn, rxn, ryn]
        elif stype == 'ys':
            pchoice = J[xp, yp, rxp, ryp] + np.max(J[xn, yn, :, ryp])
            tchoice = np.max(J[xp, yp, :, ryn]) + J[xn, yn, rxn, ryn]
        elif stype == 'xnyp':
            pchoice = J[xp, yp, rxp, ryp] + np.max(J[xn, yn, ryp, :])
            tchoice = np.max(J[xp, yp, ryn, :]) + J[xn, yn, rxn, ryn]
        elif stype == 'ynxp':
            pchoice = J[xp, yp, rxp, ryp] + np.max(J[xn, yn, :, rxp])
            tchoice = np.max(J[xp, yp, :, rxn]) + J[xn, yn, rxn, ryn]

        if pchoice > tchoice:
            if stype == 'xs':
                rxn = rxp
                ryn = int(np.where(J[xn, yn, rxn, :] == np.amax(J[xn, yn, rxn, :])))
            if stype == 'ys':
                ryn = ryp
                rxn = int(np.where(J[xn, yn, :, ryn] == np.amax(J[xn, yn, :, ryn])))
            if stype == 'xnyp':
                rxn = ryp
                ryn = int(np.where(J[xn, yn, rxn, :] == np.amax(J[xn, yn, rxn, :])))
            if stype == 'ynxp':
                ryn = rxp
                rxn = int(np.where(J[xn, yn, :, ryn] == np.amax(J[xn, yn, :, ryn])))
            pvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))

        if tchoice > pchoice:
            if stype == 'xs':
                rxp = rxn
                ryp = int(np.where(J[xp, yp, rxp, :] == np.amax(J[xp, yp, rxp, :])))
            if stype == 'ys':
                ryp = ryn
                rxp = int(np.where(J[xp, yp, :, ryp] == np.amax(J[xp, yp, :, ryp])))
            if stype == 'xnyp':
                ryp = rxn
                rxp = int(np.where(J[xp, yp, :, ryp] == np.amax(J[xp, yp, :, ryp])))
            if stype == 'ynxp':
                rxp = ryn
                ryp = int(np.where(J[xp, yp, rxp, :] == np.amax(J[xp, yp, rxp, :])))
            pvals[xid] = (xp, yp, rxp, ryp, J[xp, yp, rxp, ryp])
            pvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
        vals = [xp, yp, rxp, ryp, xn, yn, rxn, ryn]
        return vals, pvals


def Seq_edit_past_entry_comp(array, gseq):
    gseq[array[0]] = gseq[array[2]]
    gseq[array[1]] = gseq[array[3]]
    gseq[array[4]] = gseq[array[6]]
    gseq[array[5]] = gseq[array[7]]
    return gseq


