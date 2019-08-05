#!/usr/bin/env python
import numpy as np
import copy
import numpy.linalg
import matplotlib.pyplot as plt
import dcamethods as dca
import math
import random

droppath = "Projects/DCA/v2/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
fullpath = upath + droppath
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U':4}
rna = ['-', 'A', 'C', 'G', 'U']
nucs = ['A', 'C', 'G', 'U']
nuc_to_id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
g2path = upath + droppath + 'FamHJ/'
g3path = upath + droppath + '3GHJ/'


# def check_vals(x, y, pvals):
#     if len(pvals) == 0:
#         return 0, 0, 'none'
#     else:
#         for idd, pack in enumerate(pvals):
#             xp, yp, rx, ry, pval = pack
#             if x == xp:
#                 return 1, idd, 'xs'
#             elif y == yp:
#                 return 1, idd, 'ys'
#             elif x == yp:
#                 return 1, idd, 'xnyp'
#             elif y == xp:
#                 return 1, idd, 'ynxp'
#         return 0, 0, 'none'
#
#
# def past_entry_comp_goodseq(J, pvals, xn, yn):
#     tmppvals = copy.deepcopy(pvals)
#     ind, xid, stype = check_vals(xn, yn, pvals)
#     if ind == 0:
#         tmpxn, tmpyn = list(np.where(J[xn, yn, :, :] == np.amax(J[xn, yn, :, :])))
#         rxn = int(tmpxn)
#         ryn = int(tmpyn)
#         tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#         return [xn, yn, rxn, ryn], tmppvals
#     elif ind == 1:
#         xp, yp, rxp, ryp, val = tmppvals[xid]
#         tmpxn, tmpyn = list(np.where(J[xn, yn, :, :] == np.amax(J[xn, yn, :, :])))  # Indices of highest in N
#         rxn = int(tmpxn)
#         ryn = int(tmpyn)
#         if stype == 'xs':
#             pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, rxp, :])
#             tchoice = np.amax(J[xp, yp, rxn, :]) + J[xn, yn, rxn, ryn]
#
#         if stype == 'ys':
#             pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, :, ryp])
#             tchoice = np.amax(J[xp, yp, :, ryn]) + J[xn, yn, rxn, ryn]
#
#         if stype == 'xnyp':
#             pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, ryp, :])
#             tchoice = np.amax(J[xp, yp, ryn, :]) + J[xn, yn, rxn, ryn]
#
#         if stype == 'ynxp':
#             pchoice = J[xp, yp, rxp, ryp] + np.amax(J[xn, yn, :, rxp])
#             tchoice = np.amax(J[xp, yp, :, rxn]) + J[xn, yn, rxn, ryn]
#
#         if pchoice > tchoice:
#             if stype == 'xs':
#                 rxn = rxp
#                 ryn = int(np.where(J[xn, yn, rxn, :] == np.amax(J[xn, yn, rxn, :]))[0])
#             if stype == 'ys':
#                 ryn = ryp
#                 rxn = int(np.where(J[xn, yn, :, ryn] == np.amax(J[xn, yn, :, ryn]))[0])
#             if stype == 'xnyp':
#                 rxn = ryp
#                 ryn = int(np.where(J[xn, yn, rxn, :] == np.amax(J[xn, yn, rxn, :]))[0])
#             if stype == 'ynxp':
#                 ryn = rxp
#                 rxn = int(np.where(J[xn, yn, :, ryn] == np.amax(J[xn, yn, :, ryn]))[0])
#             tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#
#         if tchoice > pchoice:
#             if stype == 'xs':
#                 rxp = rxn
#                 ryp = int(np.where(J[xp, yp, rxp, :] == np.amax(J[xp, yp, rxp, :]))[0])
#             if stype == 'ys':
#                 ryp = ryn
#                 rxp = int(np.where(J[xp, yp, :, ryp] == np.amax(J[xp, yp, :, ryp]))[0])
#             if stype == 'xnyp':
#                 ryp = rxn
#                 rxp = int(np.where(J[xp, yp, :, ryp] == np.amax(J[xp, yp, :, ryp]))[0])
#             if stype == 'ynxp':
#                 rxp = ryn
#                 ryp = int(np.where(J[xp, yp, rxp, :] == np.amax(J[xp, yp, rxp, :]))[0])
#             tmppvals[xid] = (xp, yp, rxp, ryp, J[xp, yp, rxp, ryp])
#             tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#
#         if tchoice == pchoice:
#             tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#
#         vals = [xp, yp, rxp, ryp, xn, yn, rxn, ryn]
#         return vals, tmppvals
#
#
# def past_entry_comp_badseq(J, pvals, xn, yn):
#     tmppvals = copy.deepcopy(pvals)
#     ind, xid, stype = check_vals(xn, yn, pvals)
#     if ind == 0:
#         tmpxn, tmpyn = list(np.where(J[xn, yn, 1:5, 1:5] == np.amin(J[xn, yn, 1:5, 1:5])))
#         rxn = int(tmpxn) + 1
#         ryn = int(tmpyn) + 1
#         tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#         return [xn, yn, rxn, ryn], tmppvals
#     elif ind == 1:
#         xp, yp, rxp, ryp, val = tmppvals[xid]
#         tmpxn, tmpyn = list(np.where(J[xn, yn, 1:5, 1:5] == np.amin(J[xn, yn, 1:5, 1:5])))  # Indices of highest in N
#         rxn = int(tmpxn)+1
#         ryn = int(tmpyn)+1
#         if stype == 'xs':
#             pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, rxp, 1:5])
#             tchoice = np.amin(J[xp, yp, rxn, 1:5]) + J[xn, yn, rxn, ryn]
#
#         if stype == 'ys':
#             pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, 1:5, ryp])
#             tchoice = np.amin(J[xp, yp, 1:5, ryn]) + J[xn, yn, rxn, ryn]
#
#         if stype == 'xnyp':
#             pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, ryp, 1:5])
#             tchoice = np.amin(J[xp, yp, ryn, 1:5]) + J[xn, yn, rxn, ryn]
#
#         if stype == 'ynxp':
#             pchoice = J[xp, yp, rxp, ryp] + np.amin(J[xn, yn, 1:5, rxp])
#             tchoice = np.amin(J[xp, yp, 1:5, rxn]) + J[xn, yn, rxn, ryn]
#
#         if pchoice < tchoice:
#             if stype == 'xs':
#                 rxn = rxp
#                 ryn = int(np.where(J[xn, yn, rxn, 1:5] == np.amin(J[xn, yn, rxn, 1:5]))[0])+1
#             if stype == 'ys':
#                 ryn = ryp
#                 rxn = int(np.where(J[xn, yn, 1:5, ryn] == np.amin(J[xn, yn, 1:5, ryn]))[0])+1
#             if stype == 'xnyp':
#                 rxn = ryp
#                 ryn = int(np.where(J[xn, yn, rxn, 1:5] == np.amin(J[xn, yn, rxn, 1:5]))[0])+1
#             if stype == 'ynxp':
#                 ryn = rxp
#                 rxn = int(np.where(J[xn, yn, 1:5, ryn] == np.amin(J[xn, yn, 1:5, ryn]))[0])+1
#             tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#
#         if tchoice < pchoice:
#             if stype == 'xs':
#                 rxp = rxn
#                 ryp = int(np.where(J[xp, yp, rxp, 1:5] == np.amin(J[xp, yp, rxp, 1:5]))[0])+1
#             if stype == 'ys':
#                 ryp = ryn6
#                 rxp = int(np.where(J[xp, yp, 1:5, ryp] == np.amin(J[xp, yp, 1:5, ryp]))[0])+1
#             if stype == 'xnyp':
#                 ryp = rxn
#                 rxp = int(np.where(J[xp, yp, 1:5, ryp] == np.amin(J[xp, yp, 1:5, ryp]))[0])+1
#             if stype == 'ynxp':
#                 rxp = ryn
#                 ryp = int(np.where(J[xp, yp, rxp, 1:5] == np.amin(J[xp, yp, rxp, 1:5]))[0])+1
#             tmppvals[xid] = (xp, yp, rxp, ryp, J[xp, yp, rxp, ryp])
#             tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#
#         if tchoice == pchoice:
#             tmppvals.append((xn, yn, rxn, ryn, J[xn, yn, rxn, ryn]))
#
#         vals = [xp, yp, rxp, ryp, xn, yn, rxn, ryn]
#         return vals, tmppvals
#
#
# def Seq_edit_past_entry_comp(array, gseq):
#     if len(array) >= 4:
#         gseq[array[0]+0] = rna[int(array[2])]
#         gseq[array[1]+1] = rna[int(array[3])]
#         if len(array) == 8:
#             gseq[array[4]+0] = rna[int(array[6])]
#             gseq[array[5]+1] = rna[int(array[7])]
#     return gseq
#
#
# def gen_goodseq(J, H, N, norms):
#     # Get Indices of top 10 norms
#     gseq = np.full(40, ['X'], dtype=str)
#     tval = topxjnorms(J, N, norms)
#     pvals = []
#     for i in range(len(tval)):
#         x, y, z = tval[i]
#         vals, pvals = past_entry_comp_goodseq(J, pvals, x, y)
#         gseq = Seq_edit_past_entry_comp(vals, gseq)
#     for xid, x in enumerate(gseq):
#         if x == 'X':
#             gpx = np.argmax(H[xid, 1:5])
#             print(gpx)
#             gseq[xid] = rna[int(gpx) + 1]
#     print(''.join(gseq))
#     print(dca.Calc_Energy(gseq, J, H))
#     return ''.join(gseq)
#
#
# def gen_badseq(J, H, N, norms):
#     # Get Indices of top 10 norms
#     bseq = np.full(40, ['X'], dtype=str)
#     tval = topxjnorms(J, N, norms)
#     pvals = []
#     for i in range(len(tval)):
#         x, y, z = tval[i]
#         vals, pvals = past_entry_comp_badseq(J, pvals, x, y)
#         bseq = Seq_edit_past_entry_comp(vals, bseq)
#     for xid, x in enumerate(bseq):
#         if x == 'X':
#             gpx = np.argmin(H[xid, 1:5])
#             bseq[xid] = rna[int(gpx)+1]
#     print(pvals)
#     print(dca.Calc_Energy(gseq, J, H))
#     print(''.join(bseq))
#     return ''.join(bseq)
#
#
# def gen_badseq_mutt(J, H, JMutt, N, numberofnorms, **kwargs):
#     ns = 'J'
#     for key, value in kwargs.items():
#         if key == 'normsource':
#             ns = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     # Get Indices of top 10 norms
#     bseq = np.full(40, ['X'], dtype=str)
#     if ns == 'J':
#         tval = dca.TopX_JNorms(J, N, numberofnorms)
#     if ns == 'mutt':
#         tval = dca.TopX_JNorms(JMutt, N, numberofnorms)
#     pvals = []
#     for i in range(len(tval)):
#         x, y, z = tval[i]
#         vals, pvals = past_entry_comp_badseq(JMutt, pvals, x, y)
#         bseq = Seq_edit_past_entry_comp(vals, bseq)
#     for xid, x in enumerate(bseq):
#         if x == 'X':
#             gpx = np.argmin(H[xid, 1:5])
#             bseq[xid] = rna[int(gpx) + 1]
#     print(pvals)
#     print(''.join(gseq))
#     print(dca.Calc_Energy(gseq, J, H))
#     return ''.join(gseq)
#
# # Input J and number of norms
# def gen_goodseq_mutt(J, H, JMutt, N, numberofnorms, **kwargs):
#     ns = 'J'
#     for key, value in kwargs.items():
#         if key == 'normsource':
#             ns = value
#         else:
#             print('No keyword argument ' + key + ' found')
#     # Get Indices of top 10 norms
#     gseq = np.full(40, ['X'], dtype=str)
#     if ns == 'J':
#         tval = dca.TopX_JNorms(J, N, numberofnorms)
#     if ns == 'mutt':
#         tval = dca.TopX_JNorms(JMutt, N, numberofnorms)
#     pvals = []
#     for i in range(len(tval)):
#         x, y, z = tval[i]
#         vals, pvals = past_entry_comp_goodseq(JMutt, pvals, x, y)
#         gseq = Seq_edit_past_entry_comp(vals, gseq)
#     for xid, x in enumerate(gseq):
#         if x == 'X':
#             gpx = int(np.where(H[xid, :] == np.amax(H[xid, :]))[0])
#             gseq[xid] = rna[int(gpx)]
#     print(pvals)
#     print(''.join(gseq))
#     print(dca.Calc_Energy(gseq, J, H))
#     return ''.join(gseq)
#
#
# class GenerSeq:
#     def __init__(self, N, T, mut_steps=5, out_after=100, steps=10000):
#         self._history = []
#         self._T = T
#         self._beta = 1. / T
#         self._mut_steps = mut_steps
#         self._out_after = out_after
#         self._steps = steps
#         self._seq = np.random.choice(nucs, N)
#
#     def calculate_energy(self, J, h):
#         J_energy = 0.
#         h_energy = 0.
#         for i in range(1, len(self._seq)):
#             t1 = nuc_to_id[self._seq[i]]
#             h_energy += h[i, t1]
#             for j in range(i + 1, len(self._seq)):
#                 t2 = nuc_to_id[self._seq[j]]
#                 J_energy += J[i - 1, j - 2, t1, t2]
#         return J_energy + h_energy
#
#     def mutate_seq(self):
#         for m in range(self._mut_steps):
#             pos = np.random.choice(range(len(self._seq)))
#             self._seq[pos] = np.random.choice(nucs)
#         return self._seq
#
#     def run_sampling(self, J, h, outpath):
#         out = open(outpath, 'w')
#         oldene = self.calculate_energy(J, h)
#         oldseq = copy.deepcopy(self._seq)
#         for i in range(self._steps):
#             self.mutate_seq()
#             newene = self.calculate_energy(J, h)
#             p = math.exp(self._beta * (newene - oldene))
#             if random.random() < p:
#                 # accept move
#                 oldene = newene
#                 oldseq = copy.deepcopy(self._seq)
#                 if i % self._out_after:
#                     if ''.join(self._seq) not in self._history:
#                         self._history.append(''.join(self._seq))
#                         print('>' + str(i) + '-' + str(newene), file = out)
#                         print(''.join(self._seq), file = out)
#                     print(str((i / self._steps) * 100) + ' Percent Done')
#             else:
#                 self._seq = copy.deepcopy(oldseq)
#         out.close()
#
#     def run_bad_sampling(self, J, h, outpath):
#         out = open(outpath, 'w')
#         oldene = self.calculate_energy(J, h)
#         oldseq = copy.deepcopy(self._seq)
#         for i in range(self._steps):
#             self.mutate_seq()
#             newene = self.calculate_energy(J, h)
#             p = math.exp(self._beta * (newene - oldene))
#             if random.random() > p:
#                 # accept move
#                 oldene = newene
#                 oldseq = copy.deepcopy(self._seq)
#                 if i % self._out_after:
#                     if ''.join(self._seq) not in self._history:
#                         self._history.append(''.join(self._seq))
#                         print('>' + str(i) + '-' + str(newene), file=out)
#                         print(''.join(self._seq), file=out)
#                     print((str((i / self._steps) * 100)) + ' Percent Done')
#             else:
#                 self._seq = copy.deepcopy(oldseq)
#         out.close()


N = 40
q = 5
# 2 Group Method


fam = 7
#Paths
#BadBinders Import
bJp = g2path + str(fam) + 'BP.j'
bHp = g2path + str(fam) + 'BP.h'
# Goodbinders
gHp = g3path + str(fam) + 'gg.h'
gJp = g3path + str(fam) + 'gg.j'
# Best Binders
vJp = g3path + str(fam) + 'vg.j'
vHp = g3path + str(fam) + 'vg.h'
# Import Matrices
bJ = dca.sortjmat_plmDCA(bJp, N, q)
bH = dca.sorthmat_plmDCA(bHp, N, q)

gH = dca.sorthmat_plmDCA(gHp, N, q)
gJ = dca.sortjmat_plmDCA(gJp, N, q)

vJ = dca.sortjmat_plmDCA(vJp, N, q)
vH = dca.sorthmat_plmDCA(vHp, N, q)


bJpos = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='+')
bJneg = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='-')

gJpos = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='+')
gJneg = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='-')

vJneg = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='-')
vJpos = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='+')

# testseqpath = testpath + '7thfull.txt'
testseqpath = g2path + '7test.txt'
testseqpath8 = g2path + '8test.txt'

print(np.sum(vJ))
print(np.sum(gJ))
print(np.sum(bJ))


H = (vH + gH - bH) # S1
J = (vJ + gJ - bJ) # S1

Jmutt = dca.HJ_Mutant(J, H, N, q)
gmuttseq = dca.gen_goodseq(Jmutt, H, N, 200)
dca.Calc_Energy(gmuttseq, J, H)
dca.gen_badseq(J, H, N, 10)

# hybridH, hybridJ = dca.Binder_Comp_JH(J, bJ, H, bH, N, 5, htype='good', jnormpct=40)
# hybridmuttJ = dca.HJ_Mutant(hybridJ, hybridH, N, q)
# dca.gen_goodseq_mutt(hybridJ, hybridH, hybridmuttJ, 40, 500, normsource='J')


# T = 0.1
# gen = GenerSeq(N, T, mut_steps=5, out_after=10000, steps=100000000)
# gen.run_sampling(hybridJ, hybridH, out)





