import dcamethods as dca
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

upath = "/home/jonah/Dropbox (ASU)/"
droppathv3 = "Projects/DCA/ThrombinAptamers/v4/10percent/"
apath = upath + "Projects/DCA/ThrombinAptamers/v3/"
vpath = upath + droppathv3
cpath = "/home/jonah/Desktop/Current/"
# tenpath = cpath + 'v4/10percent/'
# spath = cpath + 'v4/byfam/'

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
j901p = vpath + 's901.j'
h901p = vpath + 's901.h'
j902p = vpath + 's902.j'
h902p = vpath + 's902.h'
j900p = vpath + 's900.j'
h900p = vpath + 's900.h'
j901 = dca.sortjmat_plmDCA(j901p, N, q)
j902 = dca.sortjmat_plmDCA(j902p, N, q)
h901 = dca.sorthmat_plmDCA(h901p, N, q)
h902 = dca.sorthmat_plmDCA(h902p, N, q)
j900 = dca.sortjmat_plmDCA(j900p, N, q)
h900 = dca.sorthmat_plmDCA(h900p, N, q)
H = 2*h901
J = dca.TopJNorms_Jmatrix(j901, N, q, 329)[0]


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


r1 = 'GGGGATGGGGGGGTGGAGGACTAGGTTGGGTAGGGTGGTG'
r2 = 'AGGGTAGGAGTGGATGATGCGTAGGTTGGGTAGGGTGGTC'
r3 = 'GTAGGATGGGTAGGGTGGTCAGGGATGATGTGTGGTAGGC'
r4 = 'AGGGATGATGTGTGGTAGGCGCGGGCGGTACGTGGGTGTG'
h1 = 'AGGGATGGGGTGGTGGAGGCCTAGGTTGGGTAGGGTGGTG'
h2 = 'AGGGTAGGTGTGGATGATGCCTAGGTTGGGTAGGGTGGTG'
h3 = 'CTAGGATGGGTAGGGTGGTGAGGGTTGATGTGTGGTAGGC'
h4 = 'AGGGATGATGTGTGGTAGGCGGGGTCGGTACGTGGGTGGG'
l1 = 'GGCUATGGGGGGGTGGAGGACTAGGTTGGGTAGGGTCGTG'
l2 = 'AGTTTAGGAGTGGATGATGCGTAGCTTGGGTAGGGTGGTC'
l3 = 'GTCGGACGGGTAGGGCGGTCAGGGATGATGTGTGGTAGGC'
l4 = 'AGGGATCATGTGTGGCCGGCGCGGGCGGTACGTGGGTGTG'

allseqpath = apath + str(fam) + 'all.txt'
# print(dca.Calc_Energy(h2, J, H))
baf, beq = dca.Fasta_Read_GB(allseqpath)
# results = []
selseqs = []
sela = []
# for i in range(10):
#     sels = np.random.choice(range(len(beq)))
#     selseqs.append(beq[sels])
#     sela.append(baf[sels])
# for xid, x in enumerate(selseqs):
#     print(x, 'E=', round(dca.Calc_Energy(x, J, H), 2), 'A =', sela[xid])
    # a, s = dca.avg_mut_energy(x, 3, 10000, N, J, H)
    # print('3 Mutations: avg=', a, 'std=', s)
    # l, h = dca.ThreeMutation_checker(x, J, H)
    # print('Highest:', h[-1][0], 'E=', h[-1][1])
    # print('Lowest:', l[0][0], 'E=', l[0][1])
x = 'GTAGGATGGGTAGGGTGGTCGTAGGATGGGTAGGGTGGTC'
print(x, 'E=', round(dca.Calc_Energy(x, J, H), 2), 'A =', 10)
# l, h = dca.ThreeMutation_checker(x, J, H)
# print('Highest:', h[-1][0], 'E=', h[-1][1])
# print('Lowest:', l[0][0], 'E=', l[0][1])

hi = 'GTAGGATGGGTAGGGTGGTCCTAGGTTGGGTAGGGTGGTG'
lo = 'GTAGGACGGGTAGGGCGGTCGTAGCATGGGTAGGGTGGTC'
print(dca.Calc_Energy(lo, J, H))

# r1, r2 = dca.ensemble_checker(allseqpath, hi), dca.ensemble_checker(allseqpath, lo)
# print(r1, r2)
# print(dca.mut_loc(x, lo))
# all7 = apath + '7all.txt'
# gseqs = [f1, f2, f3, f4, f5]
# hiseqs = [f1, f2]
# rseqs = [r1, r2, r3, r4]
# hseqs = [h1, h2, h3, h4]
# lseqs = [l1, l2, l3, l4]
# for x in gseqs:
#     print(round(dca.avgdis_bb(x, allseqpath), 2))
# for x in rseqs:
#     h, l = dca.ThreeMutation_checker(x, J, H)
#     print(h, l)
# for x in hiseqs:
#     print(dca.ThreeMutation_checker(x, J, H))

# dca.Raw_wRscore(J, H, '/home/jonah/Desktop/fam7scoring.png', all7)
# dca.Plot_Seq_Aff_v_E(J, H, '/home/jonah/Desktop/fam7avg.png', all7)
f8 = 'AGGGTAGGTGTGGATTATGCCTAGGATGGGTAGGGTGGTG'
f9 = 'AGGGTAGGTGTGGATTATGCCTAGGATGGGTAGGGTGGTT'
f10 = 'CGGGTAGGTGTGGATTATGCCTAGCATGGGTAGGGTGGTG'

f11 = 'AGGGATGATGTGTGGTAGGGCTAGGATGGGTAGGGTGGTG'
f12 = 'AGGGATGATGGGTGGTAGGGCTAGGATGGGTAGGGTGGTG'
f13 = 'TGGGATGATGTGTGGTAGGCCTAGCATGGGTAGGGTGGTG'

# print(dca.Calc_Energy(hbm, J, H))


# 17G - 87.53

b1 = 'GGGGGTTGGGCGGGATGGGCTTGGGTGGTGTAGGTTGGCG'
b2 = 'GCGGGTTGGGCAGGATCAGCTTGGGTGGTGCAGGTTCGCG'
# bs = [b1, b2]
# for x in bs:
#     print(dca.ensemble_checker(allseqpath, x))
t = 'AGGGATGATGTGGATGACGCCTAGGATGGGTAGGGTGGTG'
m = 'AGGGATGATGTGTGGTAGGCGTAGGATGGGTAGGGTGGTC'
m2 = 'AGGGTAGGTGTGGATGATGCCTAGGATGGGTAGGGTGGTG'

p2 = 'AGGGTAGGTGTGGATGATGCCTAGGATGGGTGGTGGGGTG'
p1 =   'AGGGATGATGTGTGGTAGGCGTAGGATGGGTGGGGTGGGA'
pmix = 'AGGGATGATGGTTGGTAGGCGTAGGATGGGTAGGGTGGTA'
#       1   5    10   15   20   25   30   35
# pos = dca.all_mut_possibilities(list(p2))
# results = []
# for x in pos:
#     if dca.ensemble_checker(allseqpath, x)[0][1] == 0.925:
rg = 'GGGGATGGGGGGGTGGAGGACTAGGTTGGGTAGGGTGGTG'
lg = 'GGCUATGGGGGGGTGGAGGACTAGGTTGGGTAGGGTCGTG'
rb = 'GCGGGTTGGGCAGGATGGGCTTGGGTGGTGCAGGTTGGCG'
lb = 'GCGGGTTGGGCAGGATCAGCTTGGGTGGTGCAGGTTCGCG'
ub = 'GGGGGTTGGGCGGGATGGGCTTGGGTGGTGTAGGTTGGCG'

cl ='AGGGTGGGAGTGGATGATGCCTAGGATGGGTAGGGTGGTG'
# r = dca.find_closest_seqs('AGGGTGGGAGCGGGGGACGCCTAGGTTGGGTAGGGTGGTG', allseqpath, N)
# print(r)
# print(dca.mut_loc(cl, 'AGGGTGGGAGCGGGGGACGCCTAGGTTGGGTAGGGTGGTG'))
# print(dca.Calc_Energy(pmix, J, H))
# print(dca.mut_loc(rg, lg))
# print(dca.mut_loc(rb, lb))
# print(dca.mut_loc(rb, ub))
#         results.append((x, e))
# results.sort(key=lambda tup: tup[1])
# interest = results[-10:]
# for x in interest:
#     print(x)
#     print(dca.mut_loc(p2, x[0]))

# r= dca.ensemble_checker(allseqpath, 'AGGGATGATGTGTGGTAGGCTAGGGTTGGTGTGGGTGGCG')
# print(r)

# s, e = dca.rand_mut('GTAGGATGGGTAGGGTGGTCAGGGATGATGTGTGGTAGGC', 2, J, H)
# print(s,'E=', round(e,2))


# print(dca.Calc_Energy('AGGGATGATGTGTGGTAGGCGCGGGCGGTACGTGGGTGTG', J, H))
# for i in range(6, 7):
#     e, s = dca.avg_mut_energy('AGGGATGATGTGTGGTAGGCGCGGGCGGTACGTGGGTGT', i, 100000, N, J, H)
#     print('E=', round(e, 2), 'std=', round(s, 2))

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

g1 = 'AGGGTAGGTGTGGGGTATGCCTAGGATGGGTAGGGTGGTG'
g2 = 'AGGGATGATGTGTGGTAGGCGTAGGATGGGTGGGGTGGGA'
g3 = 'AGGGTAGATGTGTAGGATGCCTAGGATGGGTAGGGTGGTG'
g4 = 'AGGGATGATGGTTGGTAGGCGTAGGATGGGTAGGGTGGTA'
g5 = 'AGGGATGATGTGGATTAGGCCTAGGATGGGTAGGGTGGTG'
g6 = 'AGGGTGGGAGCGGGGGACGCCTAGGTTGGGTAGGGTGGTG'
g7 = 'CGGGTAGGTGTGGATTATGCCTAGCATGGGTAGGGTGGTG'
g9 = 'GTAGGACGGGTAGGGCGGTCGTAGCATGGGTAGGGTGGTC'
g10 = 'GGGGGTTGGGCGGGATGGGCTTGGGTGGTGTAGGTTGGCG'
g11 = 'GCGGGTTGGGCAGGATCAGCTTGGGTGGTGCAGGTTCGCG'
gs = [g1, g2, g3, g4, g5, g6, g7, g9, g10, g11]
for x in gs:
    print('avg_dis_gb=', round(dca.avgdis_gb(x, allseqpath), 2), 'avg_dis_bb=', round(dca.avgdis_bb(x, allseqpath), 2))
