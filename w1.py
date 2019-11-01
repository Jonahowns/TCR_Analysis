import dcamethods as dca
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

upath = "/home/jonah/Dropbox (ASU)/"
droppathv3 = "Projects/DCA/v3/"
v3path = upath + droppathv3
cpath = "/home/jonah/Desktop/Current/"
tenpath = cpath + 'v4/10percent/'
spath = cpath + 'v4/byfam/'

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
H = 2*h901
J = dca.TopJNorms_Jmatrix(j901, N, q, 329)[0]
Jnfp = spath + 'hbnf.j'
Hnfp = spath + 'hbnf.h'
adjp = spath + 'ad.j'
adhp = spath + 'ad.h'
ddjp = spath + 'dd.j'
ddhp = spath + 'dd.h'
bdjp = spath + 'bd.j'
bdhp = spath + 'bd.h'
adj = dca.sortjmat_plmDCA(adjp, N, q)
adh = dca.sorthmat_plmDCA(adhp, N, q)
bdj = dca.sortjmat_plmDCA(bdjp, N, q)
bdh = dca.sorthmat_plmDCA(bdhp, N, q)
ddj = dca.sortjmat_plmDCA(ddjp, N, q)
ddh = dca.sorthmat_plmDCA(ddhp, N, q)
jnf = dca.sortjmat_plmDCA(Jnfp, N, q)
hnf = dca.sorthmat_plmDCA(Hnfp, N, q)

# fig, ax = plt.subplots(1, 3)
# addisp = dca.FullJ_disp(adj, N, q)
# dddisp = dca.FullJ_disp(ddj, N, q)
# bddisp = dca.FullJ_disp(bdj, N, q)
# dca.Fig_FullJ(ax[0], 'D-D', dddisp, N, q)
# dca.Fig_FullJ(ax[1], 'A-D', addisp, N, q)
# dca.Fig_FullJ(ax[2], 'B-D', bddisp, N, q)
# plt.savefig(spath + 'crosstalkin.png', dpi=600)
hb =  'AGGGATGATGTGTGGTAGGCCTAGGATGGGTAGGGTGGTG'
hbm = 'AGGGATGATGTGTGGTGGGCCTAGGATGGGTAGGGTTGTG'
#      1   5    10   15   20   25   30   35
f1 = 'AGGGATGATGTGGATGATGCCTAGGATGGGTAGGGTGGTG'
f2 = 'AGGGTAGGTGTGGGGTATGCCTAGGATGGGTAGGGTGGTG'
f3 = 'AGGGTAGGTGTGGAGTAGGCCTAGGATGGGTAGGGTGGTG'
f4 = 'AGGGTAGATGTGTAGGATGCCTAGGATGGGTAGGGTGGTG'

f5 = 'AGGGATGATGTGGATTAGGCCTAGGATGGGTAGGGTGGTG'
f6 = 'AGGGATGATGTGTTGGAGGCCTAGGATGGGTAGGGTGGTG'
f7 = 'AGGGAAGATGTGGGGTAGGCCTAGGATGGGTAGGGTGGTG'

f8 = 'AGGGTAGGTGTGGATTATGCCTAGGATGGGTAGGGTGGTG'
f9 = 'AGGGTAGGTGTGGATTATGCCTAGGATGGGTAGGGTGGTT'
f10 = 'CGGGTAGGTGTGGATTATGCCTAGCATGGGTAGGGTGGTG'

f11 = 'AGGGATGATGTGTGGTAGGGCTAGGATGGGTAGGGTGGTG'
f12 = 'AGGGATGATGGGTGGTAGGGCTAGGATGGGTAGGGTGGTG'
f13 = 'TGGGATGATGTGTGGTAGGCCTAGCATGGGTAGGGTGGTG'

# print(dca.Calc_Energy(hbm, J, H))


# 17G - 87.53



allseqpath = cpath + str(fam) + 'all.txt'

# r= dca.ensemble_checker(allseqpath, 'AGGGATGATGTGTGGTAGGCTAGGGTTGGTGTGGGTGGCG')
# print(r)
results = dca.find_closest_seqs('CGGGATGGGGTGAGGTAGGCATAGGATGGGTAGGGTGGTT', allseqpath, N)
print(results)
# print(dca.avg_mut_energy('AGGGTAGGTGTGGATGATGCCTAGGATGGGTAGGGTGGTG', 4, 100000, N, J, H))
# print(dca.avg_mut_energy('AGGGTAGGTGTGGATGATGCCTAGGATGGGTAGGGTGGTG', 5, 100000, N, J, H))
# print(dca.avg_mut_energy('AGGGTAGGTGTGGATGATGCCTAGGATGGGTAGGGTGGTG', 6, 100000, N, J, H))

# h, l = dca.TwoMutation_Bad_binder_checker('AGGGTAGGTGTGGATTATGCCTAGGATGGGTAGGGTGGTG', J, H, allseqpath)
# print(h, l)
# dca.Raw_wRscore(jnf, hnf, '/home/jonah/Desktop/nftrial.png', allseqpath)
# dca.Diff_Raw_wRscore(jnf, hnf, '/home/jonah/Desktop/nftrial.png', ['nofam', 'other'], allseqpath, spath + 'hbnofam.txt')
# cxxnorms = dca.sort_4d_gen_matrix('/home/jonah/Desktop/cxxnorms.txt', N-1, N-1, q, q)
# plt.imshow(cxxnorms)
# plt.savefig('/home/jonah/Desktop/cxxnormvis.png', dpi=400)
# gHp = v3path + str(fam) + 'go.h'
# gJp = v3path + str(fam) + 'go.j'
# gH = dca.sorthmat_plmDCA(gHp, N, q)
# gJ = dca.sortjmat_plmDCA(gJp, N, q)



# J, vals, cutoff = dca.TopJNorms_Jmatrix(j901, N, q, 331)
# diff = J - cxxnorms
# indices =[]
# for i in range(N-1):
#     for j in range(N-1):
#         for k in range(q):
#             for l in range(q):
#                 if abs(diff[i,j,k,l]) < 0.0002:
#                     diff[i,j,k,l] = 0.0
#                 if diff[i, j, k, l] != 0.0:
#                     print(i,j)
#                     found = False
#                     if indices:
#                         for w in indices:
#                             if w == (i, j):
#                                 found = True
#                                 break
#                     if found == False:
#                         indices.append((i, j))
# print(indices)
# allseqpath = cpath + str(fam) + 'all.txt'
# outp = spath

# dca.write_fams(outp, allseqpath)
# seqs = dca.Fasta_Read_Aff_fams(allseqpath)
# outp = spath + 'ddcombo.txt'
# o = open(outp, 'w')
# count = 0
# for s, a in seqs:
#     count += 1
#     print('>seq', count, '-', a, sep='', file=o)
#     print(s, file=o)







# print(diff[0, 0, 1, 1])
# print(diff.sum())
# disp = dca.FullJ_disp(j901, N, q)
# plt.imshow(disp, cmap='seismic')
# plt.savefig('/home/jonah/Desktop/nnstripped.png', dpi=400)
# dca.analyze_MC_seqs(86.0, allseqpath, N, J, H, '/home/jonah/Desktop/Current/aptamerfinal/esortv3.txt', '/home/jonah/Desktop/Current/aptamerfinal/msortv3.txt',
#                     '/home/jonah/Desktop/Current/aptamerfinal/outrun1v3.txt', '/home/jonah/Desktop/Current/aptamerfinal/outrun2v3.txt', '/home/jonah/Desktop/Current/aptamerfinal/outrun3v2.txt',
#                     '/home/jonah/Desktop/Current/aptamerfinal/outrun4v3.txt', '/home/jonah/Desktop/Current/aptamerfinal/outrun5v3.txt')
# H = 2*h901
'''
allseqpath = cpath + str(fam) + 'all.txt'
s10p = tenpath + 's10.txt'
s90p = tenpath + 's90.txt'

jdisp = dca.FullJ_disp(J, N, q)

# fig, ax = plt.subplots(1)
# dca.Plot_Seq_Aff_v_E(J, H, '/home/jonah/Desktop/aptfigs/avgHw.png', allseqpath)
# dca.Fig_Distribution_w_Cutoff(ax, 'Distribution of Norms J Matrix', vals, cutoff)
# dca.Fig_FullJ(ax, 'Initial J', jdisp, N, q)
# dca.Fig_FullH(ax, 'Iniitial H', h901, N, q)
# plt.savefig('/home/jonah/Desktop/aptfigs/distJ.png', dpi=400)
actualhighestscoring = 'AGGGTAGGTGTGGATGATGCCTAGGATGGGTAGGTTGGTG'

high, low = dca.TwoMutation_Bad_binder_checker(actualhighestscoring, J, H, allseqpath)
print(high)
# print(low)
s1 = 'CGGGATGATGTGTGGTAGGCCTATCATGGGTAGGGTGGTG'
s2 = 'UGGGATGATGTGTGGTAGCCCTATGATGGGTAGGGTGGTG'
s3 = 'UGGGATCATGTGTGGTAGGCCTATGATGGGTAGGGTGGTG'
s4 = 'UGGGATGATGTGTGGTAGGCCTATAATGGGTAGGGTGGTG'
'''
# results = dca.ensemble_checker(allseqpath, s1, s2, s3, s4)
# print(results)

# b, g = dca.find_closest_seqs('AGGGAUGAUGUGUAUGAUGCCUAGGUUGGGUAGGGUGGUG', allseqpath, N)
# print(dca.Calc_Energy(b, J, H), dca.Calc_Energy(g, J, H))
# sseqpath = '/home/jonah/Desktop/Current/v4/8v4poss.txt'

# cutoff = 83
# actualhighestscoring = 'AGGGATGATGTGTGGTAGGCCTATGATGGGTAGGGTGGTG'
# ssels = []
# seqs = dca.Fasta_Read_SeqOnly(sseqpath)
# print(seqs)
# for xid, x in enumerate(seqs):
#     if dca.Calc_Energy(seqs[xid], J, H) > cutoff:
#         ssels.append(seqs[xid])
#
# Get rid of any duplicate Seqs
# fseqs = dca.prune_alignment(ssels, 1.0)
# afseqs = []
# nrgs, simm = [], []
# for x in fseqs:
#     tmpa, simscore, tmpb = dca.ensemble_checker(allseqpath, x)[0]
#     if simscore != 1.0:
#         afseqs.append(x)
#         simm.append(simscore)
#
# for x in afseqs:
#     nrgs.append(dca.Calc_Energy(x, J, H))
#     simm.append(dca.sim_score(actualhighestscoring, x))
#
# out = '/home/jonah/Desktop/Current/v4/n332gbrelaxed.txt'
# o = open(out, 'w')
# for xid, x in enumerate(afseqs):
#     print(x, simm[xid], nrgs[xid], file=o)
# o.close()


