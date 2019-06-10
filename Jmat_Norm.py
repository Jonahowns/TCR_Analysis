#!/usr/bin python


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import gaussian_kde
import sys
import copy

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


def full_jdisplay(J, N, q):
    Jdisp = np.full(((N-1)*q, (N-1)*q), 0.0)
    for i in range(N-1):
        for j in range(N-1):
            for k in range(q):
                for l in range(q):
                    if J[i, j, k, l] != 0.0:
                        Jdisp[i*q+k, j*q+l] = J[i, j, k, l]
                    else:
                        Jdisp[i*q+k, j*q+l] = 0.0
    return Jdisp


def topxjnorms(J, N, x):
    jnorm = np.full((N-1, N-1), 0.0)
    vals = []
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
            if jnorm[i, j] != 0.0:
                vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
    vals.sort(key=lambda tup: tup[2])
    ind = int(-2 - x)
    top10 = vals[ind:-1]
    print(ind, -1)
    print(vals)
    print(vals[ind:-1])
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


def fig_fullJ(subplot, clustid, mat, n, cmap):
    subplot.title.set_text('Jmat Top 90% Values Cluster: ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.imshow(mat, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
    subplot.set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    subplot.set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    subplot.set_xticklabels(np.arange(2, n, 1))
    subplot.set_yticklabels(np.arange(1, n - 1, 1))
    subplot.grid(True, color='g', lw=0.1)
    subplot.set_ylabel('i')
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)


def fig_fullJ_RNA(subplot, famid, mat, n, cmap):
    subplot.title.set_text('Jmat Full RNA Fam: ' + str(famid))
    subplot.title.set_size(fontsize=6)
    subplot.imshow(mat, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
    subplot.set_xticks(np.arange(-.5, (n - 2) * 5, 5))
    subplot.set_yticks(np.arange(-.5, (n - 2) * 5, 5))
    subplot.set_xticklabels(np.arange(2, n, 1))
    subplot.set_yticklabels(np.arange(1, n - 1, 1))
    subplot.grid(True, color='g', lw=0.1)
    subplot.set_ylabel('i')
    subplot.set_xlabel('j')
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)


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


def fig_distofnorms(subplot, clustid, vals, tval):
    deN = gaussian_kde(vals)
    xd1 = np.linspace(0, 2, 100)
    subplot.plot(xd1, deN(xd1), color='r')
    subplot.plot(vals, [0.01] * len(vals), '|', color='k')
    subplot.set_xlabel('Norm Value')
    subplot.grid(True)
    subplot.title.set_text('Distribution of Norms Clust ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.axvline(x=tval)


def fig_distofnorms_RNA(subplot, clustid, vals, tval):
    deN = gaussian_kde(vals)
    xd1 = np.linspace(0, 2, 100)
    subplot.plot(xd1, deN(xd1), color='r')
    subplot.plot(vals, [0.01] * len(vals), '|', color='k')
    subplot.set_xlabel('Norm Value')
    subplot.grid(True)
    subplot.title.set_text('Distribution of Norms Family ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.axvline(x=tval)


def fig_seqlogoplot(filepath, subplot, clustid):
    fsl1 = mpimg.imread(filepath)
    subplot.imshow(fsl1)
    subplot.axis('off')
    subplot.title.set_text('SeqLogo Clust ' + str(clustid))
    subplot.title.set_size(fontsize=6)


def fig_seqlogoplot_RNA(filepath, subplot, clustid):
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


def IndJij_RNA(subplot, J, x, y, famid):
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


def IndJij_RNA_wColorBar(subplot, J, x, y, famid):
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


def top10norms_figure_RNA(famid):
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


