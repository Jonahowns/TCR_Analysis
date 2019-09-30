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
#Paths

# Goodbinders
j901p = tenpath + 's901.j'
h901p = tenpath + 's901.h'
j902p = tenpath + 's902.j'
h902p = tenpath + 's902.h'
j900p = tenpath + 's900.j'
h900p = tenpath + 's900.h'
j901 = dca.sortjmat_plmDCA(j901p, N, q)
j902 = dca.sortjmat_plmDCA(j902p, N, q)
h901 = dca.sorthmat_plmDCA(h901p, N, q)
h902 = dca.sorthmat_plmDCA(h902p, N, q)
j900 = dca.sortjmat_plmDCA(j900p, N, q)
h900 = dca.sorthmat_plmDCA(h900p, N, q)
# gHp = v3path + str(fam) + 'go.h'
# gJp = v3path + str(fam) + 'go.j'
# gH = dca.sorthmat_plmDCA(gHp, N, q)
# gJ = dca.sortjmat_plmDCA(gJp, N, q)

J = dca.TopJNorms_Jmatrix(j901, N, q, 332)[0]
H = 2*h901

allseqpath = cpath + str(fam) + 'all.txt'
s10p = tenpath + 's10.txt'
s90p = tenpath + 's90.txt'


sseqpath = '/home/jonah/Desktop/Current/v4/8v4poss.txt'

cutoff = 83
actualhighestscoring = 'AGGGATGATGTGTGGTAGGCCTATGATGGGTAGGGTGGTG'
ssels = []
seqs = dca.Fasta_Read_SeqOnly(sseqpath)
print(seqs)
for xid, x in enumerate(seqs):
    if dca.Calc_Energy(seqs[xid], J, H) > cutoff:
        ssels.append(seqs[xid])

# Get rid of any duplicate Seqs
fseqs = dca.prune_alignment(ssels, 1.0)
afseqs = []
nrgs, simm = [], []
for x in fseqs:
    tmpa, simscore, tmpb = dca.ensemble_checker(allseqpath, x)[0]
    if simscore != 1.0:
        afseqs.append(x)
        simm.append(simscore)

for x in afseqs:
    nrgs.append(dca.Calc_Energy(x, J, H))
    # simm.append(dca.sim_score(actualhighestscoring, x))

out = '/home/jonah/Desktop/Current/v4/n332gbrelaxed.txt'
o = open(out, 'w')
for xid, x in enumerate(afseqs):
    print(x, simm[xid], nrgs[xid], file=o)
o.close()


