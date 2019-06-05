import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys

#macpath = "/Users/Amber/Dropbox (ASU)/"
droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
fullpath = upath + droppath


aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


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


def jdisplay(J, N, q):
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


def jnormtval(J, N, q):
    jnorm = np.full((N-2, N-2), 0.0)
    jdisp = np.full((N - 2, N - 2), 0.0)
    for i in range(0, N-2):
        for j in range(0, N-2):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
    tval = np.percentile(jnorm, 80)
    for i in range(0, N-2):
        for j in range(0, N-2):
            if jnorm[i, j] >= tval:
                jdisp[i, j] = jnorm[i, j]
    return jdisp


def jtopval(J, N, q):
    Jdisp = np.full(((N - 2) * q, (N - 2) * q), 0.0)
    val_80 = np.percentile(J, 99.9)
    for i in range(0, N - 2):
        for j in range(0, N - 2):
            for k in range(q):
                for l in range(q):
                    if J[i, j, k, l] > val_80:
                        Jdisp[i * q + k, j * q + l] = J[i, j, k, l]
                    else:
                        Jdisp[i * q + k, j * q + l] = 0.0
    return Jdisp


def htopval(H, N, q):
    Hdisp = np.full((N, q), 0.0)
    val = np.percentile(H, 90)
    for i in range(0, N):
        for j in range(0, q):
            if H[i, j] > val:
                Hdisp[i, j] = H[i, j]
    return Hdisp


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


def fig_fullH(subplot, clustid, mat, n, cmap):
    # H1
    subplot.imshow(mat, cmap=cmap, aspect='equal')
    subplot.title.set_text('Hmat Cluster: ' + str(clustid))
    subplot.title.set_size(fontsize=6)
    subplot.set_xticks(np.arange(0, 21, 1))
    subplot.set_yticks(np.arange(0, n, 1))
    subplot.set_xticklabels(aa)
    subplot.set_yticklabels(np.arange(1, n+1, 1))
    subplot.set_xlabel('Amino Acid')
    subplot.set_ylabel('i')


def comp_mats_full_J(clustid1, clustid2):
    clustpath1 = fullpath + 'FullSeq/Clust' + str(clustid1) + '/'
    analysispath = fullpath + 'FullSeq/Analysis/'
    clustpath2 = fullpath + 'FullSeq/Clust' + str(clustid2) + '/'
    afasta1 = clustpath1 + 'afam.fasta'
    afasta2 = clustpath2 + 'afam.fasta'
    # SeqLogos
    logo1 = clustpath1 + 'fsl.png'
    logo2 = clustpath2 + 'fsl.png'
    # Import Photos
    fsl1 = mpimg.imread(logo1)
    fsl2 = mpimg.imread(logo2)
    # Matrix Paths
    J1p = clustpath1 + str(clustid1) + 'full.j'
    H1p = clustpath1 + str(clustid1) + 'full.h'
    J2p = clustpath2 + str(clustid2) + 'full.j'
    H2p = clustpath2 + str(clustid2) + 'full.h'

    n1 = getn(afasta1)
    n2 = getn(afasta2)

    if n1 != n2:
        print('These Sequences are not of the same length')
        sys.exit()

    J1 = sortjmat(J1p, n1, 21)
    H1 = sorthmat(H1p, n1, 21)
    J1d = jnormtval(J1, n1, 21)
    H1d = htopval(H1, n1, 21)

    J2 = sortjmat(J2p, n2, 21)
    H2 = sorthmat(H2p, n2, 21)
    J2d = jnormtval(J2, n2, 21)
    H2d = htopval(H2, n2, 21)

    JDIFF = np.subtract(J1d, J2d)
    HDIFF = np.subtract(H1d, H2d)

    # Plotting
    fig, ax = plt.subplots(3, 3, figsize=(18, 12), gridspec_kw={'width_ratios': [1, 1, 2], 'height_ratios': [2, 0.75, 1]}
                           , constrained_layout=True)
    # fig.subplots_adjust(hspace=1.0)
    # fig.tight_layout()
    # J Params
    plt.setp(ax[0, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 0].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[2, 0].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[0, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    # H params
    plt.setp(ax[0, 1].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[2, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[0, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 2].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[1, 2].get_yticklabels(), rotation='horizontal', fontsize=6)

    plt.suptitle('J and H Matrix Differences ' + str(clustid1) + ' v ' + str(clustid2))
    cmap = 'RdGy'
    # J1
    ax[0, 0].title.set_text('Jmat Top 90% Values Cluster: ' + str(clustid1))
    ax[0, 0].title.set_size(fontsize=6)
    pos = ax[0, 0].imshow(J1d, cmap='YlOrRd', aspect='equal', vmin=0, vmax=4.0)
    j1bar = fig.colorbar(pos, ax=ax[0, 0], fraction=0.046, pad=0.04)
    ax[0, 0].set_xticks(np.arange(-.5, (n1 - 2), 1))
    ax[0, 0].set_yticks(np.arange(-.5, (n1 - 2), 1))
    ax[0, 0].set_xticklabels(np.arange(2, n1, 1))
    ax[0, 0].set_yticklabels(np.arange(1, n1 - 1, 1))
    ax[0, 0].grid(True, color='g', lw=1.0)
    ax[0, 0].set_ylabel('i')
    # H1
    ax[1, 0].imshow(H1d, cmap=cmap, aspect='equal')
    ax[1, 0].title.set_text('Hmat Cluster: ' + str(clustid1))
    ax[1, 0].title.set_size(fontsize=6)
    ax[1, 0].set_xticks(np.arange(0, 21, 1))
    ax[1, 0].set_yticks(np.arange(0, n1, 1))
    ax[1, 0].set_xticklabels(aa)
    ax[1, 0].set_yticklabels(np.arange(1, n1+1, 1))
    ax[1, 0].set_xlabel('Amino Acid')
    ax[1, 0].set_ylabel('i')
    # J2
    ax[0, 1].title.set_text('Jmat Top 90% Values Cluster: ' + str(clustid2))
    ax[0, 1].title.set_size(fontsize=6)
    pos = ax[0, 1].imshow(J2d, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
    #j2bar = fig.colorbar(pos, ax=ax[0, 1], fraction=0.046, pad=0.04)
    ax[0, 1].set_xticks(np.arange(-.5, (n2 - 2), 1))
    ax[0, 1].set_yticks(np.arange(-.5, (n2 - 2), 1))
    ax[0, 1].set_xticklabels(np.arange(2, n2, 1))
    ax[0, 1].set_yticklabels(np.arange(1, n2 - 1, 1))
    ax[0, 1].grid(True, color='g', lw=0.5)
    ax[0, 1].set_ylabel('i')
    # H2
    ax[1, 1].imshow(H2d, cmap=cmap, aspect='equal')
    ax[1, 1].title.set_text('Hmat Cluster: ' + str(clustid2))
    ax[1, 1].title.set_size(fontsize=6)
    ax[1, 1].set_xticks(np.arange(0, 21, 1))
    ax[1, 1].set_yticks(np.arange(0, n2, 1))
    ax[1, 1].set_xticklabels(aa)
    ax[1, 1].set_yticklabels(np.arange(1, n2 + 1, 1))
    ax[1, 1].set_xlabel('Amino Acid')
    ax[1, 1].set_ylabel('i')
    # JDiff
    ax[0, 2].title.set_text('Jmat Top 90% Values Difference')
    ax[0, 2].title.set_size(fontsize=6)
    pos = ax[0, 2].imshow(JDIFF, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
    jdbar = fig.colorbar(pos, ax=ax[0, 2], fraction=0.046, pad=0.04)
    jdbar.set_ticks([-1, 0, 1])
    jdbar.set_ticklabels([str(clustid1), 'Neither', str(clustid2)])
    ax[0, 2].set_xticks(np.arange(-.5,  n2-2, 1))
    ax[0, 2].set_yticks(np.arange(-.5, (n2 - 2), 1))
    ax[0, 2].set_xticklabels(np.arange(2, n2, 1))
    ax[0, 2].set_yticklabels(np.arange(1, n2 - 1, 1))
    ax[0, 2].grid(True, color='g', lw=0.5)
    ax[0, 2].set_ylabel('i')
    # HDiff
    pos = ax[1, 2].imshow(HDIFF, cmap=cmap, aspect='equal')
    hdbar = fig.colorbar(pos, ax=ax[1, 2], fraction=0.046, pad=0.04)
    hdbar.set_ticks([-1, 0, 1])
    hdbar.set_ticklabels([str(clustid1), 'Neither', str(clustid2)])
    ax[1, 2].title.set_text('Hmat Cluster: ' + str(clustid2))
    ax[1, 2].title.set_size(fontsize=6)
    ax[1, 2].set_xticks(np.arange(0, 21, 1))
    ax[1, 2].set_yticks(np.arange(0, n2, 1))
    ax[1, 2].set_xticklabels(aa)
    ax[1, 2].set_yticklabels(np.arange(1, n2 + 1, 1))
    ax[1, 2].set_xlabel('Amino Acid')
    ax[1, 2].set_ylabel('i')
    # ColorBar Params
    j1bar.ax.tick_params(labelsize=8)
    #j2bar.ax.tick_params(labelsize=2)
    jdbar.ax.tick_params(labelsize=10)
    # IMAGE SEQUENCE LOGOS
    ax[2, 0].imshow(fsl1)
    ax[2, 1].imshow(fsl2)
    fig.delaxes(ax[2, 2])
    # ax[2, 2].imshow(gsl)
    plt.savefig(analysispath + 'difftrial.png', dpi=600)
    # plt.savefig(analysispath + 'Clust' + str(clustid) + 'analysis.png', dpi=800)

clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]



comp_mats_full_J(1,4)