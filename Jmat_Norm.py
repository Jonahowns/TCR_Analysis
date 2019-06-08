import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import gaussian_kde
import seaborn as sns
import sys

#macpath = "/Users/Amber/Dropbox (ASU)/"
# droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
fullpath = upath + droppath
clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
rna = ['-', 'A', 'C', 'G', 'U']
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U':4}
rnan = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'U'}


def htopval(H, N, q):
    Hdisp = np.full((N, q), 0.0)
    val = np.percentile(H, 90)
    for i in range(0, N):
        for j in range(0, q):
            if H[i, j] > val:
                Hdisp[i, j] = H[i, j]
    return Hdisp


def getn(fastafile):
    o =open(fastafile,'r')
    o.readline()
    seq = o.readline().rstrip()
    n = len(list(seq))
    o.close()
    return n


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


def jnorm(J, N):
    jnorm = np.full((N-1, N-1), 0.0)
    jdisp = np.full((N-1, N-1), 0.0)
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
    tval = np.percentile(jnorm, 80)
    for i in range(N-1):
        for j in range(N-1):
            if jnorm[i, j] >= tval:
                jdisp[i, j] = jnorm[i, j]
    return jdisp


def topxjnorms(J, N, x):
    jnorm = np.full((N-1, N-1), 0.0)
    vals=[]
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
            vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
    vals.sort(key=lambda tup: tup[2])
    ind = -2 - x
    top10 = vals[ind:-1]
    return top10


def jnormtvalwdist(J, N, q):
    jnorm = np.full((N-1, N-1), 0.0)
    jdisp = np.full((N-1, N-1), 0.0)
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
    tval = np.percentile(jnorm, 80)
    vals=[]
    for i in range(N-1):
        for j in range(N-1):
            if i < j:
                vals.append(jnorm[i, j])
            if jnorm[i, j] >= tval:
                jdisp[i, j] = jnorm[i, j]
    return jdisp, vals, tval


def fig_fullJnorm(subplot, clustid, mat, n, cmap):
    subplot.title.set_text('Jmat Top 80% Norms Cluster: ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.imshow(mat, cmap=cmap, aspect='equal', vmin=0, vmax=4)
    subplot.set_xticks(np.arange(-.5, (n - 1), 1))
    subplot.set_yticks(np.arange(-.5, (n - 1), 1))
    subplot.set_xticklabels(np.arange(2, n+1, 1))
    subplot.set_yticklabels(np.arange(1, n, 1))
    subplot.grid(True, color='g', lw=1.0)
    subplot.set_ylabel('i')
    subplot.set_xlabel('j')
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)


def fig_fullJnorm_RNA(subplot, clustid, mat, n, cmap):
    subplot.title.set_text('Jmat Top 80% Family: ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.imshow(mat, cmap=cmap, aspect='equal', vmin=0, vmax=1)
    subplot.set_xticks(np.arange(-.5, (n - 1), 1))
    subplot.set_yticks(np.arange(-.5, (n - 1), 1))
    subplot.set_xticklabels(np.arange(2, n+1, 1))
    subplot.set_yticklabels(np.arange(1, n, 1))
    subplot.grid(True, color='g', lw=0.5)
    subplot.set_ylabel('i')
    subplot.set_xlabel('j')
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)


def fig_fullH(subplot, clustid, mat, n, cmap):
    # H1
    subplot.imshow(mat, cmap=cmap, aspect='equal')
    subplot.title.set_text('Hmat Cluster: ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)
    subplot.set_xticks(np.arange(0, 21, 1))
    subplot.set_yticks(np.arange(0, n, 1))
    subplot.set_xticklabels(aa)
    subplot.set_yticklabels(np.arange(1, n+1, 1))
    subplot.set_xlabel('Amino Acid')
    subplot.set_ylabel('i')


def fig_fullH_RNA(subplot, clustid, mat, n, cmap):
    # H1
    subplot.imshow(mat.T, cmap=cmap, aspect='equal')
    subplot.title.set_text('Hmat Family: ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)
    subplot.set_yticks(np.arange(0, 5, 1))
    subplot.set_xticks(np.arange(0, n, 1))
    subplot.set_yticklabels(rna)
    subplot.set_xticklabels(np.arange(1, n+1, 1))
    subplot.set_ylabel('Base ID')
    subplot.set_xlabel('i')


def distofnorms(subplot, clustid, vals, tval):
    deN = gaussian_kde(vals)
    xd1 = np.linspace(0, 2, 100)
    subplot.plot(xd1, deN(xd1), color='r')
    subplot.plot(vals, [0.01] * len(vals), '|', color='k')
    subplot.set_xlabel('Norm Value')
    subplot.grid(True)
    subplot.title.set_text('Distribution of Norms Clust ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.axvline(x=tval)


def distofnorms_RNA(subplot, clustid, vals, tval):
    deN = gaussian_kde(vals)
    xd1 = np.linspace(0, 2, 100)
    subplot.plot(xd1, deN(xd1), color='r')
    subplot.plot(vals, [0.01] * len(vals), '|', color='k')
    subplot.set_xlabel('Norm Value')
    subplot.grid(True)
    subplot.title.set_text('Distribution of Norms Family ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.axvline(x=tval)


def seqlogoplot(filepath, subplot, clustid):
    fsl1 = mpimg.imread(filepath)
    subplot.imshow(fsl1)
    subplot.axis('off')
    subplot.title.set_text('SeqLogo Clust ' + str(clustid))
    subplot.title.set_size(fontsize=6)


def seqlogoplot_RNA(filepath, subplot, clustid):
    fsl1 = mpimg.imread(filepath)
    subplot.imshow(fsl1)
    subplot.axis('off')
    subplot.title.set_text('SeqLogo Family: ' + str(clustid))
    subplot.title.set_size(fontsize=6)


def jmatshow_clust(clust1, clust2, clust3, clust4):
    # File Paths
    clustpath1 = fullpath + 'FullSeq/Clust' + str(clust1) + '/'
    clustpath2 = fullpath + 'FullSeq/Clust' + str(clust2) + '/'
    clustpath3 = fullpath + 'FullSeq/Clust' + str(clust3) + '/'
    clustpath4 = fullpath + 'FullSeq/Clust' + str(clust4) + '/'
    analysispath = fullpath + 'FullSeq/Analysis/'
    afasta1 = clustpath1 + 'afam.fasta'
    afasta2 = clustpath2 + 'afam.fasta'
    afasta3 = clustpath3 + 'afam.fasta'
    afasta4 = clustpath4 + 'afam.fasta'
    # SeqLogos
    logo1 = clustpath1 + 'fsl.png'
    logo2 = clustpath2 + 'fsl.png'
    logo3 = clustpath3 + 'fsl.png'
    logo4 = clustpath4 + 'fsl.png'
    # Matrix Paths
    J1p = clustpath1 + str(clust1) + 'full.j'
    H1p = clustpath1 + str(clust1) + 'full.h'
    J2p = clustpath2 + str(clust2) + 'full.j'
    H2p = clustpath2 + str(clust2) + 'full.h'
    J3p = clustpath3 + str(clust3) + 'full.j'
    H3p = clustpath3 + str(clust3) + 'full.h'
    J4p = clustpath4 + str(clust4) + 'full.j'
    H4p = clustpath4 + str(clust4) + 'full.h'
    # Get N
    N = getn(afasta1)
    N = getn(afasta2)
    N = getn(afasta3)
    N = getn(afasta4)
    # Get Matrices Ready
    J1 = sortjmat(J1p, N, 21)
    H1 = sorthmat(H1p, N, 21)
    J1d, vals1, rc1 = jnormtvalwdist(J1, N, 21)
    H1d = htopval(H1, N, 21)

    J2 = sortjmat(J2p, N, 21)
    H2 = sorthmat(H2p, N, 21)
    J2d, vals2, rc2 = jnormtvalwdist(J2, N, 21)
    H2d = htopval(H2, N, 21)

    J3 = sortjmat(J3p, N, 21)
    H3 = sorthmat(H3p, N, 21)
    J3d, vals3, rc3 = jnormtvalwdist(J3, N, 21)
    H3d = htopval(H3, N, 21)

    J4 = sortjmat(J4p, N, 21)
    H4 = sorthmat(H4p, N, 21)
    J4d, vals4, rc4 = jnormtvalwdist(J4, N, 21)
    H4d = htopval(H4, N, 21)

    #Create Figure
    fig, ax = plt.subplots(4, 4, figsize=(10, 8), constrained_layout=True)
    cmap = 'YlOrRd'
    # J Matrices
    fig_fullJnorm(ax[0,0], clust1, J1d, N, cmap)
    fig_fullJnorm(ax[0,1], clust2, J2d, N, cmap)
    fig_fullJnorm(ax[0,2], clust3, J3d, N, cmap)
    fig_fullJnorm(ax[0,3], clust4, J4d, N, cmap)

    # H Matrices
    fig_fullH(ax[1,0], clust1, H1d, N, cmap)
    fig_fullH(ax[1,1], clust2, H2d, N, cmap)
    fig_fullH(ax[1,2], clust3, H3d, N, cmap)
    fig_fullH(ax[1,3], clust4, H4d, N, cmap)

    # Seq Logos
    seqlogoplot(logo1, ax[2, 0], clust1)
    seqlogoplot(logo2, ax[2, 1], clust2)
    seqlogoplot(logo3, ax[2, 2], clust3)
    seqlogoplot(logo4, ax[2, 3], clust4)

    # Distribution of Norms
    distofnorms(ax[3, 0], clust1, vals1, rc1)
    distofnorms(ax[3, 1], clust2, vals2, rc2)
    distofnorms(ax[3, 2], clust3, vals3, rc3)
    distofnorms(ax[3, 3], clust4, vals4, rc4)

    figname = str(clust1) + '-' + str(clust2) + '-' + str(clust3) + '-' + str(clust4) + 'top80norms.png'
    plt.savefig(analysispath + figname, dpi=600)


def jmatshow_genseqs():
    analysispath = fullpath
    # SeqLogos
    logo1 = fullpath + '5sl.png'
    logo2 = fullpath + '6sl.png'
    logo3 = fullpath + '7sl.png'
    logo4 = fullpath + '8sl.png'
    # Matrix Paths
    J1p = fullpath + '5j'
    H1p = fullpath + '5h'
    J2p = fullpath + '6j'
    H2p = fullpath + '6h'
    J3p = fullpath + '7j'
    H3p = fullpath + '7h'
    J4p = fullpath + '8j'
    H4p = fullpath + '8h'
    # N
    N = 40
    # Get Matrices Ready
    J1 = sortjmat(J1p, N, 5)
    H1 = sorthmat(H1p, N, 5)
    J1d, vals1, rc1 = jnormtvalwdist(J1, N, 5)
    H1d = htopval(H1, N, 5)

    J2 = sortjmat(J2p, N, 5)
    H2 = sorthmat(H2p, N, 5)
    J2d, vals2, rc2 = jnormtvalwdist(J2, N, 5)
    H2d = htopval(H2, N, 5)

    J3 = sortjmat(J3p, N, 5)
    H3 = sorthmat(H3p, N, 5)
    J3d, vals3, rc3 = jnormtvalwdist(J3, N, 5)
    H3d = htopval(H3, N, 5)

    J4 = sortjmat(J4p, N, 5)
    H4 = sorthmat(H4p, N, 5)
    J4d, vals4, rc4 = jnormtvalwdist(J4, N, 5)
    H4d = htopval(H4, N, 5)

    #Create Figure
    fig, ax = plt.subplots(4, 4, figsize=(10, 8), constrained_layout=True, gridspec_kw={'height_ratios': [1, 0.5, 0.5, 0.5]})
    cmap = 'YlOrRd'
    # J Matrices
    fig_fullJnorm_RNA(ax[0,0], 5, J1d, N, cmap)
    fig_fullJnorm_RNA(ax[0,1], 6, J2d, N, cmap)
    fig_fullJnorm_RNA(ax[0,2], 7, J3d, N, cmap)
    fig_fullJnorm_RNA(ax[0,3], 8, J4d, N, cmap)

    # H Matrices
    fig_fullH_RNA(ax[1,0], 5, H1d, N, cmap)
    fig_fullH_RNA(ax[1,1], 6, H2d, N, cmap)
    fig_fullH_RNA(ax[1,2], 7, H3d, N, cmap)
    fig_fullH_RNA(ax[1,3], 8, H4d, N, cmap)

    # Seq Logos
    seqlogoplot_RNA(logo1, ax[2, 0], 5)
    seqlogoplot_RNA(logo2, ax[2, 1], 6)
    seqlogoplot_RNA(logo3, ax[2, 2], 7)
    seqlogoplot_RNA(logo4, ax[2, 3], 8)

    # Distribution of Norms
    distofnorms_RNA(ax[3, 0], 5, vals1, rc1)
    distofnorms_RNA(ax[3, 1], 6, vals2, rc2)
    distofnorms_RNA(ax[3, 2], 7, vals3, rc3)
    distofnorms_RNA(ax[3, 3], 8, vals4, rc4)

    figname = 'RNAFAMs.png'
    plt.savefig(analysispath + figname, dpi=600)


def IndJij(subplot, J, x, y, famid):
    subplot.imshow(J[x, y, :, :], cmap='seismic', vmin=-0.5, vmax=0.5)
    subplot.set_xticks(np.arange(-.5, 4.5, 1))
    subplot.set_yticks(np.arange(-.5, 4.5, 1))
    subplot.set_xticklabels(['-', 'A', 'C', 'G', 'U'])
    subplot.set_yticklabels(['-', 'A', 'C', 'G', 'U'])
    # subplot.tick_params(axis='both', which='major', labelsize=4)
    # subplot.tick_params(axis='both', which='minor', labelsfiize=4)
    subplot.grid(True, color='r', lw=0.1)
    subplot.title.set_text('Fam ' + str(famid) + ' ' + 'Pair: ' + str(x + 1) + ' and ' + str(y + 2))
    subplot.title.set_size(fontsize=6)


def IndJij_wColorBar(subplot, J, x, y, famid):
    pos = subplot.imshow(J[x, y, :, :], cmap='seismic', vmin=-0.5, vmax=0.5)
    plt.colorbar(pos, ax=subplot, fraction=0.046, pad=0.04)
    subplot.set_xticks(np.arange(-.5, 4.5, 1))
    subplot.set_yticks(np.arange(-.5, 4.5, 1))
    subplot.set_xticklabels(['-', 'A', 'C', 'G', 'U'])
    subplot.set_yticklabels(['-', 'A', 'C', 'G', 'U'])
    # subplot.tick_params(axis='both', which='major', labelsize=4)
    # subplot.tick_params(axis='both', which='minor', labelsfiize=4)
    subplot.grid(True, color='r', lw=0.1)
    subplot.title.set_text('Fam ' + str(famid) + ' ' + 'Pair: ' + str(x + 1) + ' and ' + str(y + 2))
    subplot.title.set_size(fontsize=6)


def top10norms_figure(famid):
    analysispath = fullpath
    # Matrix Paths
    Jp = fullpath + str(famid) + 'j'
    # N
    N = 40
    # Get Matrix Ready
    J = sortjmat(Jp, N, 5)
    # Get Indices of top 10 norms
    jx = topxjnorms(J, N, 10)

    fig, ax = plt.subplots(2, 5, constrained_layout=True)
    for i in range(10):
        x, y, z = jx[i]
        j = i % 5
        k = 0
        if i == 0:
            IndJij_wColorBar(ax[k, j], J, x, y, famid)
        else:
            if i > 4:
                k = 1
            IndJij(ax[k, j], J, x, y, famid)

    fig.suptitle('Highest Jij Norms')
    plt.savefig(analysispath + str(famid) + 'famtop10.png', dpi=600)


# def checkprevVals_goodseq(J, pvals, x, y, gseq):
#     for idd, pack in enumerate(pvals):
#         xp, yp, rx, ry, pval = pack
#         val = np.amax(J[x, y, :, :])  # Highest in N w/ no constraints
#         gpx, gpy = list(np.where(J[x, y, :, :] == val))  # Indices of highest in N
#         if xp == x: # If base x is the same as base x in another Top 10 Jij Interaction
#             pchoice = J[xp, yp, rx, ry] + np.max(J[x, y, rx, :]) # Highest of P + Highest in N w x in state from P
#             tchoice = np.max(J[xp, yp, gpx, :]) + J[x, y, gpx, gpy] # Highest of P w x in state from N + highest in N
#             if pchoice > tchoice:
#                 xN = rx    # Chose Previous x state as x state for both N and P
#                 gseq[x] = rna[xN]
#                 gseq[yp] = rna[ry]
#                 val = np.amax(J[x, y, xN, :])
#                 gpy = int(np.where(J[x, y, xN, :] == val)[0])
#                 yN = int(gpy)
#                 gseq[y] = rna[int(gpy)]
#                 # Already have indices and value for P just need to Add N
#                 pvals.append((x, y, xN, yN, val))
#             else:
#                 gseq[y] = rna[int(gpy)]
#                 gseq[x] = rna[int(gpx)]
#                 xP, xN = int(gpx), int(gpx)
#                 yN = int(gpy)
#                 val = np.amax(J[x, yp, gpx, :])
#                 gpy = int(np.where(J[x, yp, gpx, :] == val)[0])
#                 gseq[yp] = rna[int(gpy)]
#                 yP = int(gpy)
#                 # Remove P pval arguments
#                 del pvals[idd]
#                 # Append new P arguments
#                 pvals.append((xp, yp, xP, yP, val))
#                 # Append new N arguments
#                 pvals.append((x, y, xN, yN, J[x, y, xN, yN]))
#         if yp == y: # If base y is the same as base y in another Top 10 Jij Interaction
#             pchoice = J[xp, yp, rx, ry] + np.max(J[x, y, :, ry]) # Highest of P + Highest in N w x in state from P
#             tchoice = np.max(J[xp, yp, :, gpy]) + J[x, y, gpx, gpy] # Highest of P w x in state from N + highest in N
#             if pchoice > tchoice:
#                 gseq[xp] = rna[rx]
#                 gseq[yp] = rna[ry]
#                 xP, yP = rx, ry
#                 yN = yP
#                 val = np.amax(J[x, y, :, yN])
#                 gpx = int(np.where(J[x, y, :, yN] == val)[0])
#                 xN = int(gpx)
#                 gseq[x] = rna[int(gpx)]
#                 # Already have indices and value for P just need to Add N
#                 pvals.append((x, y, xN, yN, val))
#             else:
#                 gseq[y] = rna[int(gpy)]
#                 gseq[x] = rna[int(gpx)]
#                 xN, yN = int(gpx), int(gpy)
#                 yP = yN
#                 val = np.amax(J[xp, yp, :, yP])
#                 gpx = int(np.where(J[xp, yp, :, yP] == val)[0])
#                 xP = int(gpx)
#                 gseq[xp] = rna[xP]
#                 # Remove P pval arguments
#                 del pvals[idd]
#                 # Append new P arguments
#                 pvals.append((xp, yp, xP, yP, val))
#                 # Append new N arguments
#                 pvals.append((x, y, xN, yN, J[x, y, xN, yN]))
#         if xp==y or yp==x:
#             if xp == y:
#                 pchoice = J[xp, yp, rx, ry] + np.max(J[x, y, :, rx])  # Highest of P + Highest in N w y in highest x state from P
#                 tchoice = np.max(J[xp, yp, gpy, :]) + J[x, y, gpx, gpy]  # Highest of P w x in highest y state + highest in N
#                 if pchoice > tchoice:
#                     gseq[xp] = rna(int(rx))
#                     gseq[yp] = rna(int(ry))
#                     xP, yP = rx, ry
#                     xP = yN
#                     val = np.amax(J[x, y, :, yN])
#                     gpx = int(np.where(J[x, y, :, yN] == val)[0])
#                     xN = int(gpx)
#                     gseq[x] = rna[int(gpy)]
#                     # Already have indices and value for P just need to Add N
#                     pvals.append((x, y, xN, yN, val))
#                 else:
#                     gseq[y] = rna[int(gpy)]
#                     gseq[x] = rna[int(gpx)]
#                     xN, yN = int(gpx), int(gpy)
#                     xP = yN    # xp == y
#                     val = np.amax(J[xp, yp, xP, :])
#                     gpy = int(np.where(J[x, y, :, gpx] == val)[0])
#                     yP = int(gpy)
#                     gseq[yp] = rna[yP]
#                     # Remove P pval arguments
#                     del pvals[idd]
#                     # Append new P arguments
#                     pvals.append((xp, yp, xP, yP, val))
#                     # Append new N arguments
#                     pvals.append((x, y, xN, yN, J[x, y, xN, yN]))
#             if yp == x:
#                 pchoice = J[xp, yp, rx, ry] + np.max(J[x, y, ry, :])  # Highest of P + Highest in N w x in state from P
#                 tchoice = np.max(J[xp, yp, :, gpx]) + J[x, y, gpx, gpy]  # Highest of P w x in state from N + highest in N
#                 if pchoice > tchoice:
#                     gseq[xp] = rna[int(rx)]
#                     gseq[yp] = rna[int(ry)]
#                     xP, yP = rx, ry
#                     xN = yP
#                     val = np.amax(J[x, y, :, yP])
#                     gpy = int(np.where(J[x, y, :, yP] == val)[0])
#                     yN = int(gpy)
#                     gseq[x] = rna[yN]
#                     # Already have indices and value for P just need to Add N
#                     pvals.append((x, y, xN, yN, val))
#                 else:
#                     gseq[y] = rna[int(gpy)]
#                     gseq[x] = rna[int(gpx)]
#                     xN, yN = int(gpx), int(gpy)
#                     yP = xN
#                     val = np.amax(J[xp, yp, :, yP])
#                     gpx = int(np.where(J[xp, yp, :, yP] == val)[0])
#                     xP = int(gpx)
#                     gseq[yp] = rna[xP]
#                     # Remove P pval arguments
#                     del pvals[idd]
#                     # Append new P arguments
#                     pvals.append((xp, yp, xP, yP, val))
#                     # Append new N arguments
#                     pvals.append((x, y, xN, yN, J[x, y, xN, yN]))
#     return pvals


def check_vals(x, y, pvals):
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
        else:
            return 0, idd, 'none'


def past_entry_comp(J, pvals, xn, yn):
    ind, xid, stype = check_vals(xn, yn, pvals)
    if ind == 0:
        return 0
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
        return [xp, yp, rxp, ryp, xn, yn, rxn, ryn], pvals


def Seq_edit_past_entry_comp(array, gseq):
    gseq[array[0]] = gseq[array[2]]
    gseq[array[1]] = gseq[array[3]]
    gseq[array[4]] = gseq[array[6]]
    gseq[array[5]] = gseq[array[7]]
    return gseq


def gen_goodseq(famid):
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
    tval = topxjnorms(J, N, 5)
    pvals = []
    for i in range(len(tval)):
        x, y, z = tval[i]
        vals, pvals = past_entry_comp(J, pvals, x, y)
        if len(vals) == 8:
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
    for x in range(dist-1):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        for y in range(dist-1):
            jbase = rnad[seq[y]]
            Jenergy += J[x, y, ibase, jbase]
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

best5 = gen_goodseq(5)
b5en = Calc_Energy(best5, J, H)
print(b5en)