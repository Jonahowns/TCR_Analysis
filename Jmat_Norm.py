import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import gaussian_kde
import seaborn as sns
import sys

#macpath = "/Users/Amber/Dropbox (ASU)/"
droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
fullpath = upath + droppath
clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


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


def distofnorms(subplot, clustid, vals, tval):
    den1 = gaussian_kde(vals)
    xd1 = np.linspace(0, 2, 100)
    subplot.plot(xd1, den1(xd1), color='r')
    subplot.plot(vals, [0.01] * len(vals), '|', color='k')
    subplot.set_xlabel('Norm Value')
    subplot.grid(True)
    subplot.title.set_text('Distribution of Norms Clust ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.axvline(x=tval)


def seqlogoplot(filepath, subplot, clustid):
    fsl1 = mpimg.imread(filepath)
    subplot.imshow(fsl1)
    subplot.axis('off')
    subplot.title.set_text('SeqLogo Clust ' + str(clustid))
    subplot.title.set_size(fontsize=6)


def jmatshow(clust1, clust2, clust3, clust4):
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
    n1 = getn(afasta1)
    n2 = getn(afasta2)
    n3 = getn(afasta3)
    n4 = getn(afasta4)
    # Get Matrices Ready
    J1 = sortjmat(J1p, n1, 21)
    H1 = sorthmat(H1p, n1, 21)
    J1d, vals1, rc1 = jnormtvalwdist(J1, n1, 21)
    H1d = htopval(H1, n1, 21)

    J2 = sortjmat(J2p, n2, 21)
    H2 = sorthmat(H2p, n2, 21)
    J2d, vals2, rc2 = jnormtvalwdist(J2, n2, 21)
    H2d = htopval(H2, n2, 21)

    J3 = sortjmat(J3p, n3, 21)
    H3 = sorthmat(H3p, n3, 21)
    J3d, vals3, rc3 = jnormtvalwdist(J3, n3, 21)
    H3d = htopval(H3, n3, 21)

    J4 = sortjmat(J4p, n4, 21)
    H4 = sorthmat(H4p, n4, 21)
    J4d, vals4, rc4 = jnormtvalwdist(J4, n4, 21)
    H4d = htopval(H4, n4, 21)

    #Create Figure
    fig, ax = plt.subplots(4, 4, figsize=(10, 8), constrained_layout=True)
    cmap = 'YlOrRd'
    # J Matrices
    fig_fullJnorm(ax[0,0], clust1, J1d, n1, cmap)
    fig_fullJnorm(ax[0,1], clust2, J2d, n2, cmap)
    fig_fullJnorm(ax[0,2], clust3, J3d, n3, cmap)
    fig_fullJnorm(ax[0,3], clust4, J4d, n4, cmap)

    # H Matrices
    fig_fullH(ax[1,0], clust1, H1d, n1, cmap)
    fig_fullH(ax[1,1], clust2, H2d, n2, cmap)
    fig_fullH(ax[1,2], clust3, H3d, n3, cmap)
    fig_fullH(ax[1,3], clust4, H4d, n4, cmap)

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

clusters.append(1)
it = iter(clusters)
show = list(zip(it, it, it, it))
for x in show:
    jmatshow(x[0], x[1], x[2], x[3])

