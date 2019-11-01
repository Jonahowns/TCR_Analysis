import dcamethods as dca
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

upath = "/home/jonah/Dropbox (ASU)/"
droppathv3 = "Projects/DCA/v3/"
v3path = upath + droppathv3
cpath = "/home/jonah/Desktop/Current/"
tenpath = cpath + 'v4/10percent/'

# corrfile = '/home/jonah/plmDCA/plmDCA_asymmetric_v2/8gs2.out'
# N = 20
# corr = dca.import_correlations(corrfile,N)
# dca.sample_corr(corr, N)
def mix_2comp_J(J1, J2, N, q):
    J = np.full((N-1, N-1, q, q), 0.0)
    J[0:19,0:19,:,:] = J1
    J[20:39, 20:39,:,:] =J2
    return J

def mix_2comp_H(H1, H2, N, q):
    H = np.full((N, q), 0.0)
    H[0:20,:] = H1
    H[20:40,:] = H2
    return H


N = 40
q = 5
fam = 8
j901p = tenpath + 's901.j'
h901p = tenpath + 's901.h'
j901 = dca.sortjmat_plmDCA(j901p, N, q)
h901 = dca.sorthmat_plmDCA(h901p, N, q)
#OG
allseqpath = cpath + str(fam) + 'all.txt'

J, vals, co = dca.TopJNorms_Jmatrix(j901, N, q, 329)
H = 2*h901
eH = np.full((N, q), 0.0)
jdisp = dca.FullJ_disp(J, N, q)
# dca.Raw_wRscore(J, H, '/home/jonah/Desktop/Current/aptamerfinal/figs/v5allseqsepwr.png', allseqpath)
# plt.close()
# dca.Raw_Aff_v_E(J, H, 'Family 8 Sequences', '/home/jonah/Desktop/Current/aptamerfinal/figs/v5allseqsep.png', allseqpath)
# plt.close()
# fig, ax = plt.subplots(1)
# dca.Fig_FullJ(ax, 'Norm J', jdisp, N, q)
# plt.savefig('/home/jonah/Desktop/Current/aptamerfinal/figs/v5JNorm.png', dpi=600)
# plt.close()
fig1, ax1 = plt.subplots(1)
dca.Fig_Distribution_w_Cutoff(ax1, 'J Norm Distribution', vals, co)
plt.savefig('/home/jonah/Desktop/Current/aptamerfinal/figs/v5JNormdist.png')
plt.close()
# dca.Top10norms_figure_DNA('Top 10 Norms Family 8', J, N, '/home/jonah/Desktop/Current/aptamerfinal/figs/v5top10.png')
# plt.close()


# bJsel, val, cut2 = dca.TopJNorms_Jmatrix(Jmix, N, q, 245)




# dca.mc_analysis_plot('/home/jonah/Desktop/Current/12gba.txt', allseqpath, N, Jsel, H, '/home/jonah/Desktop/Current/mcav2.png')
# bs = 'AGGGATGATGTGTGGTAGGCCTAGGGTGGGGAGGGTGGTG'
# mutt = 'AGGGATGATGTGTGGTAGGCGTAGGGTGGTGGGTGGGGTG'
# mut3 = 'AGGGATGATGTGTGGTAGGCGTAGGGTGGGGAGTGTGGTG'

# print(dca.Calc_Energy(bs, Jsel, H), dca.Calc_Energy(mutt, Jsel, H), dca.Calc_Energy(mut3, Jsel, H))

# results = dca.seq_breakdown_by_aff(allseqpath, Jsel, H, 1)
# print(results)

# z = mp.Queue()
# p1 = mp.Process(target=(dca.Jnorm_finder), args=(gJ, N, q, 300, 305, 1, allseqpath, z))
# p2 = mp.Process(target=(dca.Jnorm_finder), args=(gJ, N, q, 305, 310, 1, allseqpath, z))
# p3 = mp.Process(target=(dca.Jnorm_finder_bJ), args=(Jmix, N, q, 230, 250, 5, allseqpath, z))
# p4 = mp.Process(target=(dca.Jnorm_finder_bJ), args=(Jmix, N, q, 200, 250, 10, allseqpath, z))

# p1.start()
# p2.start()
# p3.start()
# p4.start()
# results = []
# r1, r2, r3 = z.get(), z.get(), z.get()
# p1.join()
# p2.join()
# p3.join()
# print(r1, r2, r3)


