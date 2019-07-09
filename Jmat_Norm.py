#!/usr/bin python


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import gaussian_kde
import sys
import copy
import math
import dcamethods as dca


#macpath = "/Users/Amber/Dropbox (ASU)/"
# droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/"
droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
testpath = upath + droppath
blpath = "/home/jonah/bl-DCA/"
bldrop = "bl78/"
fullpath = upath + droppath + bldrop
clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
rna = ['-', 'A', 'C', 'G', 'U']
dna = ['-', 'A', 'C', 'G', 'T']
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U':4}
rnan = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'U'}



# def htopval(H, N, q):
#     Hdisp = np.full((N, q), 0.0)
#     val = np.percentile(H, 90)
#     for i in range(0, N):
#         for j in range(0, q):
#             if H[i, j] > val:
#                 Hdisp[i, j] = H[i, j]
#     return Hdisp


# def getn(fastafile):
#     o =open(fastafile,'r')
#     o.readline()
#     seq = o.readline().rstrip()
#     n = len(list(seq))
#     o.close()
#     return n


# def sortjmat(file, N, q):
#     o = open(file, 'r')
#     fullmatrix = np.full((N-1, N-1, q, q), 0.0)
#     for line in o:
#         data = line.split(',')
#         fullmatrix[int(data[0])-1, int(data[1])-2, int(data[2])-1, int(data[3])-1] = float(data[4].rstrip())
#     o.close()
#     return fullmatrix
#
#
# def sorthmat(file, N, q):
#     o = open(file, 'r')
#     fullmatrix = np.full((N, q), 0.0)
#     for line in o:
#         data = line.split(',')
#         fullmatrix[int(data[0]) - 1, int(data[1]) - 1] = float(data[2].rstrip())
#     o.close()
#     return fullmatrix


# def jnorm(J, N):
#     jnorm = np.full((N-1, N-1), 0.0)
#     jdisp = np.full((N-1, N-1), 0.0)
#     for i in range(N-1):
#         for j in range(N-1):
#             jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
#     tval = np.percentile(jnorm, 80)
#     for i in range(N-1):
#         for j in range(N-1):
#             if jnorm[i, j] >= tval:
#                 jdisp[i, j] = jnorm[i, j]
#     return jdisp


# def full_jdisplay(J, N, q):
#     Jdisp = np.full(((N-1)*q, (N-1)*q), 0.0)
#     for i in range(N-1):
#         for j in range(N-1):
#             for k in range(q):
#                 for l in range(q):
#                     if J[i, j, k, l] != 0.0:
#                         Jdisp[i*q+k, j*q+l] = J[i, j, k, l]
#                     else:
#                         Jdisp[i*q+k, j*q+l] = 0.0
#     return Jdisp


# def topxjnorms(J, N, x):
#     jnorm = np.full((N-1, N-1), 0.0)
#     vals = []
#     for i in range(N-1):
#         for j in range(N-1):
#             jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
#             if jnorm[i, j] != 0.0:
#                 vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
#     vals.sort(key=lambda tup: tup[2])
#     ind = int(-x)
#     top10 = vals[ind:-1]
#     print(ind, -1)
#     print(vals)
#     print(vals[ind:-1])
#     return top10


# # Returns the highest x amount of J Norms
# # Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# # norm values returned as tuple list (i, j, val) and cutoff returned as single value
# # Keyword argument percentile determines the percentile used as a cutoff
# def topxjnorms_w_dist(J, N, x, **kwargs):
#     dist=False
#     pct = 80
#     for key, value in kwargs.items():
#         if key == 'percentile':
#             pct = value
#         if key == 'dist':
#             dist = True
#     jnorm = np.full((N-1, N-1), 0.0)
#     vals = []
#     jvals =[]
#     for i in range(N-1):
#         for j in range(N-1):
#             jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
#             if jnorm[i, j] != 0.0:
#                 vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
#                 jvals.append(jnorm[i, j])
#     tval = np.percentile(jnorm, pct)
#     vals.sort(key=lambda tup: tup[2])
#     ind = int(-x)
#     top10 = vals[ind:-1]
#     if dist:
#         return top10, jvals, tval
#     else:
#         return top10


# def topxhnorms_w_dist(H, N, x, **kwargs):
#     pct = 75
#     htype = 'ind'
#     for key, value in kwargs.items():
#         if key == 'pct':
#             pct = value
#         if key == 'htype':
#             htype = value
#     if htype == 'norm':
#         hnorm = np.full((N - 1), 0.0)
#         vals = []
#         hvals = []
#         for i in range(N-1):
#             hnorm[i] = np.linalg.norm(H[i, :])
#             if hnorm[i] != 0.0:
#                 vals.append((i, hnorm[i]))
#                 hvals.append(hnorm[i])
#         tval = np.percentile(hvals, pct)
#         vals.sort(key=lambda tup: tup[1])
#     elif htype == 'ind':
#         vals = []
#         hvals = []
#         for i in range(N - 1):
#             for j in range(1, 5):
#                 vals.append((i, j, abs(H[i, j])))
#                 hvals.append(H[i, j])
#         tval = np.percentile(hvals, pct)
#         vals.sort(key=lambda tup: tup[2])
#     ind = int(-x)
#     top10 = vals[ind:-1]
#     return top10, hvals, tval


# def jnormtval(J, N, q):
#     jnorm = np.full((N-1, N-1), 0.0)
#     jdisp = np.full((N-1, N-1), 0.0)
#     for i in range(N-1):
#         for j in range(N-1):
#             jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
#     tval = np.percentile(jnorm, 80)
#     vals=[]
#     for i in range(N-1):
#         for j in range(N-1):
#             if i < j:
#                 vals.append(jnorm[i, j])
#             if jnorm[i, j] >= tval:
#                 jdisp[i, j] = jnorm[i, j]
#     return jdisp, vals, tval


# def Fig_FullJ(subplot, id, J, n, q, **kwargs):
#     cmap = 'seismic'
#     vml = -1
#     vmg = 1
#     lw = 0.1
#     xlabel = 'j'
#     ylabel = 'i'
#     fontsize = 6
#     title = 'Jmat Top 90% Values ID: ' + str(id)
#     for key, value in kwargs.items():
#         if key == 'cmap':
#             cmap = value
#         elif key == 'lw':
#             lw = value
#         elif key == 'xlabel':
#             xlabel = value
#         elif key == 'ylabel':
#             ylabel = value
#         elif key == 'ticksize':
#             fontsize = value
#         elif key == 'vml':
#             vml = value
#         elif key == 'vmg':
#             vmg = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=6)
#     subplot.imshow(J, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
#     subplot.set_xticks(np.arange(-.5, (n - 2) * q, q))
#     subplot.set_yticks(np.arange(-.5, (n - 2) * q, q))
#     subplot.set_xticklabels(np.arange(2, n, 1))
#     subplot.set_yticklabels(np.arange(1, n - 1, 1))
#     subplot.grid(True, color='g', lw=lw)
#     subplot.set_ylabel(ylabel)
#     supplot.set_xlabel(xlabel)
#     plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
#     plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)


# def fig_fullJ_RNA(subplot, famid, mat, n, cmap):
#     subplot.title.set_text('Jmat Full RNA Fam: ' + str(famid))
#     subplot.title.set_size(fontsize=6)
#     subplot.imshow(mat, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
#     subplot.set_xticks(np.arange(-.5, (n - 2) * 5, 5))
#     subplot.set_yticks(np.arange(-.5, (n - 2) * 5, 5))
#     subplot.set_xticklabels(np.arange(2, n, 1))
#     subplot.set_yticklabels(np.arange(1, n - 1, 1))
#     subplot.grid(True, color='g', lw=0.1)
#     subplot.set_ylabel('i')
#     subplot.set_xlabel('j')
#     plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
#     plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)

#
# def Fig_Jnorm(subplot, id, J, n, **kwargs):
#     cmap = 'seismic'
#     vml = 0
#     vmg = 4
#     lw = 1.0
#     xlabel = 'j'
#     ylabel = 'i'
#     fontsize = 6
#     title = 'Jmat Norms ID: ' + str(id)
#     for key, value in kwargs.items():
#         if key == 'cmap':
#             cmap = value
#         elif key == 'lw':
#             lw = value
#         elif key == 'xlabel':
#             xlabel = value
#         elif key == 'ylabel':
#             ylabel = value
#         elif key == 'ticksize':
#             fontsize = value
#         elif key == 'vml':
#             vml = value
#         elif key == 'vmg':
#             vmg = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=6)
#     subplot.imshow(J, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
#     subplot.set_xticks(np.arange(-.5, (n - 1), 1))
#     subplot.set_yticks(np.arange(-.5, (n - 1), 1))
#     subplot.set_xticklabels(np.arange(2, n+1, 1))
#     subplot.set_yticklabels(np.arange(1, n, 1))
#     subplot.grid(True, color='g', lw=lw)
#     subplot.set_ylabel(ylabel)
#     subplot.set_xlabel(xlabel)
#     plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
#     plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)


# def fig_fullJnorm_RNA(subplot, famid, mat, n, **kwargs):
#     title = 'Jmat Top 80% Family: ' + str(famid)
#     cmap = 'seismic'
#     vml = 0
#     vmg = 2
#     lw = 0.5
#     for key, value in kwargs.items():
#         if key == 'title':
#             title = value
#         if key == 'cmap':
#             cmap = value
#         if key == 'vmin':
#             vml = value
#         if key == 'vmax':
#             vmg = value
#         if key == 'lw':
#             lw = value
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=6)
#     subplot.imshow(mat, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
#     subplot.set_xticks(np.arange(-.5, (n - 1), 1))
#     subplot.set_yticks(np.arange(-.5, (n - 1), 1))
#     subplot.set_xticklabels(np.arange(2, n+1, 1))
#     subplot.set_yticklabels(np.arange(1, n, 1))
#     subplot.grid(True, color='g', lw=lw)
#     subplot.set_ylabel('i')
#     subplot.set_xlabel('j')
#     plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
#     plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)


# def fig_fullH(subplot, id, H, n, q,  **kwargs):
#     cmap = 'seismic'
#     vml = 0
#     vmg = 4
#     xl = False
#     xlabel = 'hello'
#     ylabel = 'i'
#     fontsize = 6
#     title = 'Hmat ID: ' + str(id)
#     for key, value in kwargs.items():
#         if key == 'cmap':
#             cmap = value
#         elif key == 'xlabel':
#             xl = True
#             xlabel = value
#         elif key == 'ylabel':
#             ylabel = value
#         elif key == 'ticksize':
#             fontsize = value
#         elif key == 'vml':
#             vml = value
#         elif key == 'vmg':
#             vmg = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     subplot.imshow(H, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=fontsize)
#     plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
#     plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)
#     subplot.set_xticks(np.arange(0, q+1, 1))
#     subplot.set_yticks(np.arange(0, n+1, 1))
#     if q == 21 and xl is False:
#         subplot.set_xticklabels(aa)
#         subplot.set_xlabel('Amino Acid')
#     elif q == 5 and xl is False:
#         subplot.set_xticklabels(rna)
#         subplot.set_xlabel('Base')
#     else:
#         subplot.set_xlabel(xlabel)
#     subplot.set_yticklabels(np.arange(1, n+1, 1))
#     subplot.set_ylabel(ylabel)


# def fig_fullH_RNA(subplot, famid, mat, n, **kwargs):
#     title = 'Hmat Family: ' + str(famid)
#     cmap = 'seismic'
#     vml = 0
#     vmg = 2
#     lw = 0.5
#     for key, value in kwargs.items():
#         if key == 'title':
#             title = value
#         if key == 'cmap':
#             cmap = value
#         if key == 'vmin':
#             vml = value
#         if key == 'vmax':
#             vmg = value
#         if key == 'lw':
#             lw = value
#     subplot.imshow(mat.T, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=6)
#     plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=6)
#     plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=6)
#     subplot.set_yticks(np.arange(0, 5, 1))
#     subplot.set_xticks(np.arange(0, n, 1))
#     subplot.set_yticklabels(rna)
#     subplot.set_xticklabels(np.arange(1, n+1, 1))
#     subplot.set_ylabel('Base ID')
#     subplot.set_xlabel('i')

#
# def Fig_Distribution_w_Cutoff(subplot, id, vals, tval):
#     deN = gaussian_kde(vals)
#     xd1 = np.linspace(0, 2, 100)
#     subplot.plot(xd1, deN(xd1), color='r')
#     subplot.plot(vals, [0.01] * len(vals), '|', color='k')
#     subplot.set_xlabel('Norm Value')
#     subplot.grid(True)
#     subplot.title.set_text('Distribution of Norms Clust ' + str(id))
#     subplot.title.set_size(fontsize=6)
#     subplot.axvline(x=tval)


# def Fig_Distribution_w_Cutoff(subplot, id, Values, Cutoff, **kwargs):
#     title = 'Distribution ID: ' + str(id)
#     plotcolor = 'r'
#     xml = 0
#     xmg = 2
#     xlabel = 'Value'
#     fontsize = 6
#     for key, value in kwargs.items():
#         if key == 'xlabel':
#             xlabel = value
#         elif key == 'pcolor':
#             plotcolor = value
#         elif key == 'fontsize':
#             fontsize = value
#         elif key == 'xmin':
#             xml = value
#         elif key == 'xmax':
#             xmg = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     deN = gaussian_kde(Values)
#     xd1 = np.linspace(xml, xmg, 100)
#     subplot.plot(xd1, deN(xd1), color=plotcolor)
#     subplot.plot(Values, [0.01] * len(Values), '|', color='k')
#     subplot.set_xlabel(xlabel)
#     subplot.grid(True)
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=fontsize)
#     subplot.axvline(x=Cutoff)


# def Fig_SeqLogo(Filepath, Subplot, id):
#     title = 'SeqLogo ID: ' + str(id)
#     fontsize = 6
#     for key, value in kwargs.items():
#         if key == 'title':
#             title = value
#         elif key == 'fontsize':
#             fontsize = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     fsl1 = mpimg.imread(Filepath)
#     Subplot.imshow(fsl1)
#     Subplot.axis('off')
#     Subplot.title.set_text(title)
#     Subplot.title.set_size(fontsize=fontsize)

# Specific To Cluster Analysis
def JNorm_Clust(clust1, clust2, clust3, clust4):
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




# Specific to RNA
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

# Three Types are available: 'dna', 'rna' and 'pep'
# def IndJij_RNA(subplot, J, x, y, id, **kwargs):
#     vml = -0.5
#     vmg = 0.5
#     cmap = 'seismic'
#     fontsize = 4
#     type = 'rna'
#     lw = 0.1
#     title = 'Jij ID: ' + str(id) + ' ' + 'Pair: ' + str(x + 1) + ' and ' + str(y + 2)
#     for key, value in kwargs.items():
#         if key == 'vmax':
#             vmg = value
#         elif key == 'vmin':
#             vml = value
#         elif key == 'cmap':
#             cmap = value
#         elif key == 'fontsize':
#             fontsize = value
#         elif key == 'type':
#             type = value
#         elif key == 'lw':
#             lw = value
#         elif key == 'title':
#             title = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     subplot.imshow(J[x, y, :, :], cmap=cmap, vmin=vml, vmax=vmg)
#     if type == 'rna' or type == 'dna'
#         subplot.set_xticks(np.arange(-.5, 4.5, 1))
#         subplot.set_yticks(np.arange(-.5, 4.5, 1))
#         if type == 'rna':
#             subplot.set_xticklabels(rna)
#             subplot.set_yticklabels(rna)
#         else:
#             subplot.set_xticklabels(dna)
#             subplot.set_yticklabels(dna)
#     elif type == 'pep':
#         subplot.set_xticks(np.arange(-.5, 20.5, 1))
#         subplot.set_yticks(np.arange(-.5, 20.5, 1))
#         subplot.set_xticklabels(aa)
#         subplot.set_yticklabels(aa)
#     subplot.tick_params(axis='both', which='major', labelsize=fontsize)
#     subplot.tick_params(axis='both', which='minor', labelsize=fontsize)
#     subplot.grid(True, color='r', lw=lw)
#     subplot.title.set_text(title)
#     subplot.title.set_size(fontsize=(fontsize+2))



# def IndJij_DNA_mutt_wColorBar(subplot, Jmutt, x, y, famid):
#     pos = subplot.imshow(Jmutt[x, y, :, :], cmap='seismic', vmin=-2, vmax=2)
#     plt.colorbar(pos, ax=subplot, fraction=0.046, pad=0.04)
#     subplot.set_xticks(np.arange(-.5, 4.5, 1))
#     subplot.set_yticks(np.arange(-.5, 4.5, 1))
#     subplot.set_xticklabels(['-', 'A', 'C', 'G', 'T'])
#     subplot.set_yticklabels(['-', 'A', 'C', 'G', 'T'])
#     # subplot.tick_params(axis='both', which='major', labelsize=4)
#     # subplot.tick_params(axis='both', which='minor', labelsfiize=4)
#     subplot.grid(True, color='r', lw=0.1)
#     subplot.title.set_text('Fam ' + str(famid) + ' ' + 'Pair: ' + str(x + 1) + ' and ' + str(y + 2))
#     subplot.title.set_size(fontsize=6)


# def IndJij_RNA_wColorBar(subplot, J, x, y, famid):
#     pos = subplot.imshow(J[x, y, :, :], cmap='seismic', vmin=-0.5, vmax=0.5)
#     plt.colorbar(pos, ax=subplot, fraction=0.046, pad=0.04)
#     subplot.set_xticks(np.arange(-.5, 4.5, 1))
#     subplot.set_yticks(np.arange(-.5, 4.5, 1))
#     subplot.set_xticklabels(['-', 'A', 'C', 'G', 'U'])
#     subplot.set_yticklabels(['-', 'A', 'C', 'G', 'U'])
#     # subplot.tick_params(axis='both', which='major', labelsize=4)
#     # subplot.tick_params(axis='both', which='minor', labelsfiize=4)
#     subplot.grid(True, color='r', lw=0.1)
#     subplot.title.set_text('Fam ' + str(famid) + ' ' + 'Pair: ' + str(x + 1) + ' and ' + str(y + 2))
#     subplot.title.set_size(fontsize=6)


# def top10norms_figure_RNA(id, J, N):
#     # Get Indices of top 10 norms
#     jx = topxjnorms(J, N, 10)
#
#     fig, ax = plt.subplots(2, 5, constrained_layout=True)
#     for i in range(10):
#         x, y, z = jx[i]
#         j = i % 5
#         k = 0
#         if i == 0:
#             IndJij_wColorBar(ax[k, j], J, x, y, famid)
#         else:
#             if i > 4:
#                 k = 1
#             IndJij(ax[k, j], J, x, y, famid)
#
#     fig.suptitle('Highest Jij Norms')
#     plt.savefig(analysispath + str(famid) + 'famtop10.png', dpi=600)


# def HJ_Mutant(J, H, N, q):
#     mutt = copy.deepcopy(J)
#     for x in range(N - 1):
#         for k in range(q):  # J Indices
#             mutt[x, x:N, k, :] += H[x, k]
#     for y in range(N - 1):
#         for l in range(q):  # y states
#             mutt[0:y + 1, y, :, l] += H[y + 1, l]
#     return mutt


# def top10norms_figure_DNA_mutt(famid):
#     analysispath = fullpath
#     # Matrix Paths
#     Jp = fullpath + str(famid) + 'j'
#     Hp = fullpath + str(famid) + 'h'
#     # N
#     N = 40
#     # Get Matrix Ready
#     J = sortjmat(Jp, N, 5)
#     H = sorthmat(Hp, N, 5)
#     # Get Indices of top 10 norms
#     jmutt = HJ_mutant_RNA(J, H, N)
#     jx = topxjnorms(J, N, 10)
#
#     fig, ax = plt.subplots(2, 5, constrained_layout=True)
#     for i in range(10):
#         x, y, z = jx[i]
#         j = i % 5
#         k = 0
#         if i == 0:
#             IndJij_DNA_mutt_wColorBar(ax[k, j], jmutt, x, y, famid)
#         else:
#             if i > 4:
#                 k = 1
#             IndJij_DNA_mutt(ax[k, j], jmutt, x, y, famid)
#
#     fig.suptitle('Highest Jij Mutt Norms')
#     plt.savefig(analysispath + str(famid) + 'famTop10muttdisp.png', dpi=600)


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

    H *= np.divide(1, np.sum(H))
    bH *= np.divide(1, np.sum(bH))

    print(np.sum(H))
    print(np.sum(bH))
    truH = H - bH
    return truH, J


# def Calc_Energy(seq, J, H):
#     full = list(seq)
#     dist = len(full)
#     Jenergy = 0
#     Henergy = 0
#     for x in range(1, dist):
#         ibase = rnad[seq[x]]
#         Henergy += H[x, ibase]
#         for y in range(x+1, dist):
#             jbase = rnad[seq[y]]
#             Jenergy += J[x-1, y-2, ibase, jbase]
#     energy = Jenergy + Henergy
#     return energy


# # Takes in Good Binders J and H and Bad Binders J and H
# # if htype = 'ind' compares top individual values of GB H and BB H and if they are in common removes them
# # if htype = 'norm' compares top norms of GB H and BB H and if they are in common removes them
# def mixed_HJ_w_dist(J, bJ, H, bH, N, q, **kwargs):
#     filledJnorms = (N-1)*(N-2)/2 + N-1
#     filledHind = (N * q)
#     nxj = 10 # Number of J Norms used in Comparison
#     nxh = 10 # Number of H Norms used n Comparison
#     htype = 'ind'
#     hdist = False
#     jdist = False
#     hnorms = 10
#     for key, value in kwargs.items():
#         if key == 'jnormpct':
#             nxj = math.ceil(value/100 * filledJnorms)
#         elif key == 'htype':
#             htype = value
#         elif key == 'hnormpct':
#             if htype == 'norm':
#                 nxh = math.ceil(value/100 * N)
#             if htype == 'ind':
#                 nxh = math.ceil(value/100 * filledHind)
#         elif key == 'hdist':
#             hdist = value
#         elif key == 'jdist':
#             jdist = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     topJ, valsJ, tvalJ = topxjnorms_w_dist(J, N, nxj, pct=20)
#     topBJ, valsBJ, tvalBJ = topxjnorms_w_dist(bJ, N, nxj, pct=20)
#     for xj, yj, val in topJ:
#         for xb, yb, valb in topBJ:
#             if xj == xb and yj == yb:
#                 J[xj, yj, :, :] = 0.0
#     if htype != 'good':
#         topH, valsH, tvalH = topxhnorms_w_dist(H, N, nxh, pct=20, htype=htype)
#         topBH, valsBH, tvalBH = topxhnorms_w_dist(bH, N, nxh, pct=20, htype=htype)
#         if htype == 'norm':
#             for xi, val in topH:
#                 for yb, valb in topBH:
#                     if xi == yb:
#                         H[xi, :] = 0.0
#         elif htype == 'ind':
#             for xi, yi, val in topH:
#                 for xb, yb, valb in topBH:
#                     if xi == xb and yi == yb:
#                         H[xi, yi] = 0.0
#         if hdist:
#             Hdist = [valsH, tvalH, valsBH, tvalBH]
#     else:
#         Hdist = [0]
#     if jdist:
#         Jdist = [valsJ, tvalJ, valsBJ, tvalBJ]
#     if hdist and jdist:
#         return H, J, Hdist, Jdist
#     elif hdist and not jdist:
#         return H, J, hdist
#     elif not hdist and jdist:
#         return H, J, jdist
#     elif not hdist and not jdist:
#         return H, J


def subplot_seq_aff_v_E(subplot, famid, J, H, **kwargs):
    title = 'Family: ' + str(famid) + ' Affinity vs Energy'
    for key, value in kwargs.items():
        if key == 'title':
            title = value
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
    subplot.errorbar(x, avg, err, linestyle='None', marker='^')
    subplot.set_xlabel('affinity')
    subplot.set_ylabel('Energy')
    subplot.set_title(title)


def pct_comp_fig(famid):
    N = 40
    H10, J10 = mixed_HJ_w_dist(famid, pctnorms=10, hnorms=4, htype='norm')
    H20, J20 = mixed_HJ_w_dist(famid, pctnorms=20, hnorms=8, htype='norm')
    H30, J30 = mixed_HJ_w_dist(famid, pctnorms=30, hnorms=12, htype='norm')
    H40, J40 = mixed_HJ_w_dist(famid, pctnorms=40, hnorms=16, htype='norm')
    Harr = np.array([H10, H20, H30, H40])
    Jarr = np.array([J10, J20, J30, J40])
    pctarr = np.arange(10, 50, 10)
    fig, ax = plt.subplots(7, 4, figsize=(16, 24))
    for x in range(4):
        subplot_seq_aff_v_E(ax[0, x], famid, Jarr[x, 0], Harr[x, 0], title=('Mixed Scoring pct ' + str(pctarr[x])))
        fig_fullJnorm_RNA(ax[1, x], famid, jnorm(Jarr[x, 0], N), N, lw=0.1, vmin=0, vmax=1.5, title=('J Norm pct ' + str(pctarr[x])), htype='norm')
        fig_fullH_RNA(ax[2, x], famid, Harr[x, 0], N, vmin=-0.5, vmax=0.5, title=('H Matt pct ' + str(pctarr[x])))
        fig_distofnorms_RNA(ax[3, x], famid, Jarr[x, 1], Jarr[x, 2], title=('Good Binders pct ' + str(pctarr[x])))
        fig_distofnorms_RNA(ax[4, x], famid, Jarr[x, 3], Jarr[x, 4], title=('Bad Binders pct ' + str(pctarr[x])))
        fig_distofnorms_RNA(ax[5, x], famid, Harr[x, 1], Harr[x, 2], title=('Good Binders pct ' + str(pctarr[x])))
        fig_distofnorms_RNA(ax[6, x], famid, Harr[x, 3], Harr[x, 4], title=('Bad Binders pct ' + str(pctarr[x])))
    plt.suptitle('Family ' + str(famid) + 'using Editied H Norms')
    plt.savefig("/home/jonah/Downloads/pctcompNORMH" + str(famid) + ".png", dpi=600)


def mix_score_dist(famid):
    N=40
    H, J, valsJ, valsBJ, tvalJ, tvalBJ = mixed_HJ_w_dist(famid, pctnorms=20)
    fig, ax = plt.subplots(3, 2)
    subplot_seq_aff_v_E(ax[0, 0], famid, J, H, title=('Mixed Scoring Family ' + str(famid)))
    # Jg = sortjmat(J, N, 5)
    jgnorm = jnorm(J, N)
    fig_fullJnorm_RNA(ax[0, 1], famid, jgnorm, N, title='Mixed J Norm', lw=0.1, vmin=0, vmax=2)
    fig_fullH_RNA(ax[1, 0], famid, H, N, vmin=-0.005, vmax=0.005, title='Good Binders H - Bad Binders H')
    slpath = fullpath + str(famid) + 'fullsl.png'
    fig_seqlogoplot_RNA(slpath, ax[1, 1], famid)
    fig_distofnorms_RNA(ax[2, 0], famid, valsJ, tvalJ, title='Good Binders J Norm Dist')
    fig_distofnorms_RNA(ax[2, 1], famid, valsBJ, tvalBJ, title='Bad Binders J Norm Dist')
    plt.savefig("/home/jonah/Downloads/scoring.png")


# def get_energy_seqs(famid):
#     o=open(fullpath + str(famid) + 'thfull.txt')
#     titles = []
#     seqs = []
#     for line in o:
#         if line.startswith('>'):
#             titles.append(float(line.rstrip().split('-')[1]))
#         else:
#             seqs.append(line.rstrip())
#     o.close()
#     return titles, seqs


def subplot_seq_aff_v_E_w_OtherFamSeqs(subplot, famid, J, H, **kwargs):
    title = 'Family: ' + str(famid) + ' Affinity vs Energy including Other Family Seqs'
    for key, value in kwargs.items():
        if key == 'title':
            title = value
    pfams = [5, 7, 8]
    titles = []
    seqs = []
    for x in pfams:
        titlesub, seqsub = get_energy_seqs(x)
        titles.append(titlesub)
        seqs.append(seqsub)
    energies = []
    for x in seqs:
        tmp = []
        for s in x:
            nrg = Calc_Energy(s, J, H)
            tmp.append(nrg)
        energies.append(tmp)
    api5 = list(zip(titles[0], energies[0]))
    api7 = list(zip(titles[1], energies[1]))
    api8 = list(zip(titles[2], energies[2]))
    x5 = list(set([x for (x, y) in api5]))
    x7 = list(set([x for (x, y) in api7]))
    x8 = list(set([x for (x, y) in api8]))
    x5.sort()
    x7.sort()
    x8.sort()
    y5A, y5E, y7A, y7E, y8A, y8E = ([] for i in range(6))
    for aff in x5:
        yvals5 = np.array([y for (x, y) in api5 if x==aff])
        yavg5 = yvals5.mean()
        yerr5 = np.std(yvals5)
        y5A.append(yavg5)
        y5E.append(yerr5)
    for aff in x7:
        yvals7 = np.array([y for (x, y) in api7 if x==aff])
        yavg7 = yvals7.mean()
        yerr7 = np.std(yvals7)
        y7A.append(yavg7)
        y7E.append(yerr7)
    for aff in x8:
        yvals8 = np.array([y for (x, y) in api8 if x==aff])
        yavg8 = yvals8.mean()
        yerr8 = np.std(yvals8)
        y8A.append(yavg8)
        y8E.append(yerr8)
    subplot.errorbar(x5, y5A, y5E, linestyle='None', marker='o', ecolor='r', alpha=0.5)
    subplot.errorbar(x7, y7A, y7E, linestyle='None', marker='P', ecolor='b', alpha=0.5)
    subplot.errorbar(x8, y8A, y8E, linestyle='None', marker='^', ecolor='y', alpha=0.5)
    subplot.set_xlabel('affinity')
    subplot.set_ylabel('Energy')
    subplot.set_title(title)


N = 40
q = 5
bJp = fullpath + 'jb7_vals'
bHp = fullpath + 'hb7_vals'
gJp = fullpath + 'j7_vals'
gHp = fullpath + 'h7_vals'


bJ = dca.sortjmat_blDCA(bJp, N, q)
bH = dca.sorthmat_blDCA(bHp, N, q)
gH = dca.sorthmat_blDCA(gHp, N, q)
gJ = dca.sortjmat_blDCA(gJp, N, q)

print(np.average(gH))
print(np.average(bH))
print(np.average(gJ))
print(np.average(bJ))


bHpos = dca.Sign_Seperator(bH, N, q, mattype='h', sign='+')
bHneg = dca.Sign_Seperator(bH, N, q, mattype='h', sign='-')
gHpos = dca.Sign_Seperator(gH, N, q, mattype='h', sign='+')
gHneg = dca.Sign_Seperator(gH, N, q, mattype='h', sign='-')
bJpos = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='+')
bJneg = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='-')
gJpos = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='+')
gJneg = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='-')

testseqpath = testpath + '7thfull.txt'

# H = 2*gH - bH # S1
# J = 2*gJ - bJ # S1
# H = gH # S2
# J = 2*gJ - bJpos # S2
# H = np.add(np.subtract(2*gHpos, bHpos, bHneg), 2*gHneg) #S3
# J = np.add(np.subtract(3*gJpos, 2*bJpos, bJneg), 2*gJneg) #S3

jout = '/home/jonah/Desktop/Jdesigner.txt'
hout = '/home/jonah/Desktop/Hdesigner.txt'


# H = gH
# Jdes = dca.designerJ(N, q, testseqpath)
# dca.export_J(Jdes, N, q, jout)
Hdes = dca.designerH(N, q, testseqpath)
dca.export_H(Hdes, N, q, hout)



# dca.best_seperation(J, H, N, q, testseqpath)
# J = dca.Rm_Vals_Percentage_J(2*gJ - bJ, 1.0, N, q)
# H = gH
# J = gJ
# For Family 7 it looks like 1 percent of the top values give the best R score
# For Family 8 it looks like 0.1 percent of the top values give the best R score

# H = gH
# J = dca.Rm_Vals_Percentage_J(gJ, 1, N, q)
# H = np.add(dca.Rm_Vals_Percentage_H(np.subtract(gH*Hx, bHpos), 10, N, q), bHneg)
# H = dca.Rm_Vals_Percentage_H(np.add(np.subtract(gHpos*Hx, bHpos), bHneg), 10, N, q)
# J = np.add(dca.Rm_Vals_Percentage_J(np.subtract(gJ*Jx, bJpos), 1, N, q), bJneg)
# J = dca.Rm_Vals_Percentage_J(np.add(np.subtract(gJpos*Jx, bJpos), bJneg), 1, N, q)

# (J*1.5 - bJ + J)
# J = dca.Rm_Vals_Percentage_J(np.add(np.subtract(3*gJpos, 2*bJpos, bJneg), 2*gJneg), 0.5, N, q)
# J = 2*dca.Rm_Vals_Percentage_J(np.add(np.subtract(gJpos, bJpos, bJneg), gJneg), 2, N, q)
# H = dca.Rm_Vals_Percentage_H(np.add(np.subtract(2*gHpos, bHpos, bHneg), 2*gHneg), 2, N, q)

# Honestly Not too Shabby
# J = np.add(np.subtract(gJpos, bJpos, bJneg), gJneg)
# H = np.add(np.subtract(3*gHpos, 2*bHpos, bHneg), 2*gHneg)
# H = 2*np.add(np.subtract(gHpos, bHpos, bHneg), gHneg)
# H = np.add(np.subtract(2*gHpos, bHpos, bHneg), gHneg)
# EHHHHH
# J = gJ - bJ
# H = gH - bH






# (J*1.5 - bJ + J)

# J = dca.Rm_Vals_Percentage_J(J, 2, N, q)


# J *= 2
jdisp = dca.FullJ_disp(Jdes, N, q)
fig, ax = plt.subplots(1, 2)
dca.Fig_FullJ(ax[0], 8, jdisp, N, q, vmin=-1, vmax=1, cmap='seismic')
dca.Fig_FullH(ax[1], 8, Hdes, N, q, vmin=-1, vmax=1, cmap='seismic')
plt.savefig('/home/jonah/Desktop/designerJ.png', dpi=600)
# dca.Fig_FullH(ax[1], 8, H, N, q, vmin=-1, vmax=1)
# plt.savefig('/home/jonah/Desktop/8gbpos.png', dpi=600)
# rdseq806 = 'GGGCUAAGGGCGUAGUCGGCGUAUGUUGGGUAGUUAAGUC'  #2.76366
# rdseq808 = 'GGGCUAAGGGCGUAGUCGGCGUAUGUUGGGUAGUUAAGUG'  #2.76366
# rdseq809 = 'GGGCUAAGGGCGUAGUCGGCGUAUGUUGGGUAGUUUGGUC'  #2.76366
# rdhonly  = 'AGGGUAUGGGUGUGGUGGGCUUUCGGUGGUUUGGUUGGUC'  #14.865
# rdseq801 = 'GGGCUAAGGGCGUAGUCGGCGUAUGUUGGGUAGUUAAGUG'  #2.5
# rdseq807 = 'GGGCUAAGGGCGUAGUCGGCGUAUGUUGGGUAGUUAUGUG'  #2.5
# rdseq804 = 'GGGCUAAGGGCGCGCUCGGCGUAUGUUGGGUAGUUAAGUG'  #2.07
# rdseq802 = 'GGCCUGAGGGCGUAGUCGGCGUAUGUUGGGUAGUUAAGUG'  #0.4
# rdseq803 = 'GGCCUGGGGGCGCGCUCGGCGUAUGUUGGGUAGUUAAGUG'  #0.8
#
#
# rdseq805 = 'GGCCUGGGGGCGCGCUGGGUGAUGUGUGGGCACCUAUGUC' # -6 Yikes
# E = [dca.Calc_Energy(rdseq806, J, H), dca.Calc_Energy(rdseq808, J, H), dca.Calc_Energy(rdseq809, J, H), dca.Calc_Energy(rdhonly, J, H), dca.Calc_Energy(rdseq805, J, H)]
# print(E)



# testseqpath = testpath + '7thfull.txt'
# dca.Plot_Seq_Aff_v_E(J, H, ('/home/jonah/Desktop/testenergy.png'), testseqpath)
# dca.Raw_Aff_v_E(J, H, ('/home/jonah/Desktop/raw7.png'), testseqpath)


# mix_score_dist(7)
# pct_comp_fig(7)