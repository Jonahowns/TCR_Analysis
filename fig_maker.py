import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

macpath = "/Users/Amber/Dropbox (ASU)/"
droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
fullpath = macpath+droppath


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


def create_figure(clustid):
    clustpath = macpath + droppath + 'FullSeq/Clust' + str(clustid) + '/'
    analysispath = macpath + droppath + 'FullSeq/Analysis/'
    afasta = clustpath + 'afam.fasta'
    # SeqLogos
    alogo = clustpath + 'asl.png'
    clogo = clustpath + 'csl.png'
    glogo = clustpath + 'gsl.png'
    tlogo = clustpath + 'tsl.png'
    fulllogo = clustpath + 'fsl.png'
    # Import Photos
    asl = mpimg.imread(alogo)
    csl = mpimg.imread(clogo)
    gsl = mpimg.imread(glogo)
    tsl = mpimg.imread(tlogo)
    fsl = mpimg.imread(fulllogo)
    # Matrix Paths
    aji = clustpath + str(clustid) + 'a.j'
    ahi = clustpath + str(clustid) + 'a.h'
    cji = clustpath + str(clustid) + 'c.j'
    chi = clustpath + str(clustid) + 'c.h'
    gji = clustpath + str(clustid) + 'g.j'
    ghi = clustpath + str(clustid) + 'g.h'
    tji = clustpath + str(clustid) + 't.j'
    thi = clustpath + str(clustid) + 't.h'
    fullji = clustpath + str(clustid) + 'full.j'
    fullhi = clustpath + str(clustid) + 'full.h'
    n = getn(afasta)
    aj = sortjmat(aji, n, 21)
    ah = sorthmat(ahi, n, 21)
    ajdisp = jdisplay(aj, n, 21)
    cj = sortjmat(cji, n, 21)
    ch = sorthmat(chi, n, 21)
    cjdisp = jdisplay(cj, n, 21)
    gj = sortjmat(gji, n, 21)
    gh = sorthmat(ghi, n, 21)
    gjdisp = jdisplay(gj, n, 21)
    tj = sortjmat(tji, n, 21)
    th = sorthmat(thi, n, 21)
    tjdisp = jdisplay(tj, n, 21)
    fullj = sortjmat(fullji, n, 21)
    fullh = sorthmat(fullhi, n, 21)
    fulljdisp = jdisplay(fullj, n, 21)

    # Plotting
    fig, ax = plt.subplots(5, 3, figsize=(12, 22), gridspec_kw={'width_ratios': [1, 0.5, 0.5]},constrained_layout=True)
    #fig.subplots_adjust(hspace=1.0)
    #fig.tight_layout()
    plt.setp(ax[0, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[3, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[4, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[0, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[3, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[4, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[0, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[1, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[2, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[3, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[4, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[0, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[3, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[4, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.suptitle('J and H Matrices Cluster ' + str(clustid))
    cmap = 'jet'
    # A J
    ax[0, 0].title.set_text('Jmat Afam Cluster: ' + str(clustid))
    ax[0, 0].title.set_size(fontsize=6)
    pos = ax[0, 0].imshow(ajdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    abar = fig.colorbar(pos, ax=ax[0, 0], fraction=0.046, pad=0.04)
    ax[0, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[0, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[0, 0].set_xticklabels(np.arange(2, n, 1))
    ax[0, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[0, 0].grid(True, color='g', lw=0.1)
    ax[0, 0].set_ylabel('i')
    # A H
    ax[0, 1].imshow(ah, cmap=cmap, aspect='equal')
    ax[0, 1].title.set_text('H mat Afam Cluster: ' + str(clustid))
    ax[0, 1].title.set_size(fontsize=6)
    ax[0, 1].set_xticks(np.arange(0, 21, 1))
    ax[0, 1].set_yticks(np.arange(0, n, 1))
    ax[0, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[0, 1].set_yticklabels(np.arange(1, n+1, 1))
    ax[0, 1].set_xlabel('Amino Acid')
    ax[0, 1].set_ylabel('i')
    # ax[0, 1].grid(True, color='b', lw=0.1)
    # C J
    ax[1, 0].title.set_text('Jmat Cfam Cluster: ' + str(clustid))
    ax[1, 0].title.set_size(fontsize=6)
    pos = ax[1, 0].imshow(cjdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    cbar = fig.colorbar(pos, ax=ax[1, 0], fraction=0.046, pad=0.04)
    ax[1, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[1, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[1, 0].set_xticklabels(np.arange(2, n, 1))
    ax[1, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[1, 0].grid(True, color='g', lw=0.1)
    ax[1, 0].set_ylabel('i')
    # C H
    ax[1, 1].imshow(ch, cmap=cmap, aspect='equal')
    ax[1, 1].title.set_text('H mat Cfam Cluster: ' + str(clustid))
    ax[1, 1].title.set_size(fontsize=6)
    ax[1, 1].set_xticks(np.arange(0, 21, 1))
    ax[1, 1].set_yticks(np.arange(0, n, 1))
    ax[1, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[1, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[1, 1].set_xlabel('Amino Acid')
    ax[1, 1].set_ylabel('i')
    # G J
    ax[2, 0].title.set_text('Jmat Gfam Cluster: ' + str(clustid))
    ax[2, 0].title.set_size(fontsize=6)
    pos = ax[2, 0].imshow(gjdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    gbar = fig.colorbar(pos, ax=ax[2, 0], fraction=0.046, pad=0.04)
    ax[2, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[2, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[2, 0].set_xticklabels(np.arange(2, n, 1))
    ax[2, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[2, 0].grid(True, color='g', lw=0.1)
    ax[2, 0].set_ylabel('i')
    # G H
    ax[2, 1].imshow(gh, cmap=cmap, aspect='equal')
    ax[2, 1].title.set_text('H mat Gfam Cluster: ' + str(clustid))
    ax[2, 1].title.set_size(fontsize=6)
    ax[2, 1].set_xticks(np.arange(0, 21, 1))
    ax[2, 1].set_yticks(np.arange(0, n, 1))
    ax[2, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[2, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[2, 1].set_xlabel('Amino Acid')
    ax[2, 1].set_ylabel('i')
    # T J
    ax[3, 0].title.set_text('Jmat Tfam Cluster: ' + str(clustid))
    ax[3, 0].title.set_size(fontsize=6)
    pos = ax[3, 0].imshow(tjdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    tbar = fig.colorbar(pos, ax=ax[3, 0], fraction=0.046, pad=0.04)
    ax[3, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[3, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[3, 0].set_xticklabels(np.arange(2, n, 1))
    ax[3, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[3, 0].grid(True, color='g', lw=0.1)
    ax[3, 0].set_ylabel('i')
    # T H
    ax[3, 1].imshow(th, cmap=cmap, aspect='equal')
    ax[3, 1].title.set_text('H mat Afam Cluster: ' + str(clustid))
    ax[3, 1].title.set_size(fontsize=6)
    ax[3, 1].set_xticks(np.arange(0, 21, 1))
    ax[3, 1].set_yticks(np.arange(0, n, 1))
    ax[3, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[3, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[3, 1].set_xlabel('Amino Acid')
    ax[3, 1].set_ylabel('i')
    # FULL J
    ax[4, 0].title.set_text('Jmat FULL Cluster: ' + str(clustid))
    ax[4, 0].title.set_size(fontsize=6)
    pos = ax[4, 0].imshow(fulljdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    fullbar = fig.colorbar(pos, ax=ax[4, 0], fraction=0.046, pad=0.04)
    ax[4, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[4, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[4, 0].set_xticklabels(np.arange(2, n, 1))
    ax[4, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[4, 0].set_xlabel('j')
    ax[4, 0].grid(True, color='g', lw=0.1)
    ax[4, 0].set_ylabel('i')
    # FULL H
    ax[4, 1].imshow(fullh, cmap=cmap, aspect='equal')
    ax[4, 1].title.set_text('H mat FULL Cluster: ' + str(clustid))
    ax[4, 1].title.set_size(fontsize=6)
    ax[4, 1].set_xticks(np.arange(0, 21, 1))
    ax[4, 1].set_yticks(np.arange(0, n, 1))
    ax[4, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[4, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[4, 1].set_xlabel('Amino Acid')
    ax[4, 1].set_ylabel('i')
    abar.ax.tick_params(labelsize=2)
    cbar.ax.tick_params(labelsize=2)
    gbar.ax.tick_params(labelsize=2)
    tbar.ax.tick_params(labelsize=2)
    fullbar.ax.tick_params(labelsize=2)
    #IMAGE SEQUENCE LOGOS
    ax[0, 2].imshow(asl)
    ax[1, 2].imshow(csl)
    ax[2, 2].imshow(gsl)
    ax[3, 2].imshow(tsl)
    ax[4, 2].imshow(fsl)
    plt.savefig(clustpath + 'Clust' + str(clustid) + 'analysis.png', dpi=800)
    plt.savefig(analysispath + 'Clust' + str(clustid) + 'analysis.png', dpi=800)


clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]


for clust in clusters:
    create_figure(clust)
