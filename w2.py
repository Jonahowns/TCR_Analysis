import dcamethods as dca
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

upath = "/home/jonah/Dropbox (ASU)/"
droppathv3 = "Projects/DCA/v3/"
v3path = upath + droppathv3


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
j1p = '/home/jonah/plmDCA/plmDCA_asymmetric_v2/8gs1.j'
j2p = '/home/jonah/plmDCA/plmDCA_asymmetric_v2/8gs2.j'
h1p = '/home/jonah/plmDCA/plmDCA_asymmetric_v2/8gs1.h'
h2p = '/home/jonah/plmDCA/plmDCA_asymmetric_v2/8gs2.h'
J1 = dca.sortjmat_plmDCA(j1p, 20, 5)
J2 = dca.sortjmat_plmDCA(j2p, 20, 5)
H1 = dca.sorthmat_plmDCA(h1p, 20, 5)
H2 = dca.sorthmat_plmDCA(h2p, 20, 5)


Jmix = mix_2comp_J(J1, J2, N, q)
Hmix = mix_2comp_H(H1, H2, N , q)

#OG
fam = 8
bJp = v3path + str(fam) + 'ba.j'
bHp = v3path + str(fam) + 'ba.h'
gHp = v3path + str(fam) + 'go.h'
gJp = v3path + str(fam) + 'go.j'
bJ = dca.sortjmat_plmDCA(bJp, N, q)
bH = dca.sorthmat_plmDCA(bHp, N, q)
gH = dca.sorthmat_plmDCA(gHp, N, q)
gJ = dca.sortjmat_plmDCA(gJp, N, q)

cpath = "/home/jonah/Desktop/Current/"
allseqpath = cpath + str(fam) + 'all.txt'

Jsel, vals, cut = dca.TopJNorms_Jmatrix(gJ, N, q, 305)
bJsel, val, cut2 = dca.TopJNorms_Jmatrix(Jmix, N, q, 245)

JF = Jsel

H = 2*gH

bs = 'AGGGATGATGTGTGGTAGGCCTAGGGTGGGGAGGGTGGTG'
mutt = 'AGGGATGATGTGTGGTAGGCGTAGGGTGGTGGGTGGGGTG'
mut3 = 'AGGGATGATGTGTGGTAGGCGTAGGGTGGGGAGTGTGGTG'

print(dca.Calc_Energy(bs, Jsel, H), dca.Calc_Energy(mutt, Jsel, H), dca.Calc_Energy(mut3, Jsel, H))

results = dca.seq_breakdown_by_aff(allseqpath, Jsel, H, 1)
print(results)

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


