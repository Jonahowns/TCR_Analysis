#!/usr/bin python


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import gaussian_kde
import sys
import copy
import math
import dcamethods as dca
from multiprocessing import Process

droppathw = "Projects/DCA/working/"
droppathv = "Projects/DCA/v2/"
droppathv3 = "Projects/DCA/v3/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
testpath = upath + 'Projects/DCA/GenSeqs/'
blpath = "/home/jonah/bl-DCA/"
bldrop = "bl78/"
cpath = "/home/jonah/Desktop/Current/"
g22path = upath + droppathv +"FamHJ/"
g2path = upath + droppathw
g3path = upath + droppathv + "3GHJ/"
v3path = upath + droppathv3

clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
rna = ['-', 'A', 'C', 'G', 'U']
dna = ['-', 'A', 'C', 'G', 'T']
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U': 4, 'T': 4}
rnan = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'U'}


N = 40
q = 5
# 2 Group Method


fam = 7
#Paths
#BadBinders Import
bJp = v3path + str(fam) + 'ba.j'
bHp = v3path + str(fam) + 'ba.h'
# Goodbinders
gHp = v3path + str(fam) + 'go.h'
gJp = v3path + str(fam) + 'go.j'
# Allseqs
# oHp = v3path + str(fam) + 'aa.h'
# oJp = v3path + str(fam) + 'aa.j'
#High scoring bad binders
# bSJp = v3path + str(fam) + 'bh.j'
# bSHp = v3path + str(fam) + 'bh.h'

# tJp = v3path + str(fam) + 'tb.j'
# tHp = v3path + str(fam) + 'tb.h'

#3 Group
#BadBinders Import
# bJp = g2path + str(fam) + 'ba.j'
# bHp = g2path + str(fam) + 'ba.h'
# #GoodBinders
# gHp = g3path + str(fam) + 'gg.h'
# gJp = g3path + str(fam) + 'gg.j'
# # Best Binders
# vJp = g3path + str(fam) + 'vg.j'
# vHp = g3path + str(fam) + 'vg.h'
# Import Matrices
bJ = dca.sortjmat_plmDCA(bJp, N, q)
bH = dca.sorthmat_plmDCA(bHp, N, q)

# bSJ = dca.sortjmat_plmDCA(bSJp, N, q)
# bSH = dca.sorthmat_plmDCA(bSHp, N, q)

# tJ = dca.sortjmat_plmDCA(tJp, N, q)
# tH = dca.sorthmat_plmDCA(tHp, N, q)

gH = dca.sorthmat_plmDCA(gHp, N, q)
gJ = dca.sortjmat_plmDCA(gJp, N, q)

# oH = dca.sorthmat_plmDCA(gHp, N, q)
# oJ = dca.sortjmat_plmDCA(gJp, N, q)




#
# vJ = dca.sortjmat_plmDCA(vJp, N, q)
# vH = dca.sorthmat_plmDCA(vHp, N, q)

# testseqpath = testpath + '7thfull.txt'
allseqpath = cpath + str(fam) + 'all.txt'

bJpos = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='+')
bJneg = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='-')
bHpos = dca.Sign_Seperator(bH, N, q, mattype='h', sign='+')
bHneg = dca.Sign_Seperator(bH, N, q, mattype='h', sign='-')

# H = np.full((N,q), 0.0)
H = 2*(2*(gH) - bH) # S1
J = 0.5*(2*(gJ) - bJ) # S1

# H = (oH + 0.5*(gH) - 0.5*bH) # S1
# J = (oJ + 0.5*(gJ) - 0.5*bJ) # S1

# H = (oH - bH + 2*gH)*5
# J = oJ - bJ + 2*gJ

# H = (2*oH - bSH)*5
# J = 2*oJ - bSJ

# H -= 5*bSH
# J -= 2*bSJ

# H = gH - bSH
# J = gJ - bSJ

# H = (2*(vH + gH) - bH) # S1
# J = (2*(vJ + gJ) - bJ) # S1
# norms = dca.TopX_JNorms(J, N, 10)
# print(norms)

tseq = 'AGGGAUGAUGUGUGGUAGGCCUAGGUUGGGGAGGGUGGUG'
Orig = 'AGGGUAGGUGUGGAUGAUGCCUAGGAUGGGUAGGGUGGUG' #259
M1 =   'AGGGUUGGUGGGGAUGAUGAGUAGGAUGGGUAGGGUGGUA' #226 4 modifications
M2 =   'AGAGUAGGUGUGGAUGAUGCCUAGGAUGGGCAAGGUGGUG' #211 3 modifications attacked more highly weighted, 31,32,33





tfseq =  'AGGGATGATGTGTGGTAGGCCTAGGGTGGGGAGGGTGGTG'
tf2seq = 'AGGGATGATTTGTGGTAGGCCTAGGGTGGGTAGGGTGGTG'

#         1        10   15   20   25   30
# dca.Raw_A ff_v_E(J, H, ('A vs E Family ' + str(fam)), ('/home/jonah/Desktop/' + str(fam) + 'v3trial.png'), allseqpath)
# dca.Raw_wRscore(J, H, ('/home/jonah/Desktop/' + str(fam) + 'v3Otrial.png'), allseqpath)
# dca.Plot_Seq_Aff_v_E(J, H, '/home/jonah/Desktop/7avgv3.png', allseqpath)
results = dca.seq_breakdown_by_aff(allseqpath, J, H, 1)
# results = dca.ensemble_checker(allseqpath, M2)
#results = dca.ensemble_checker(allseqpath, tf2seq)
# tot, results = dca.Calc_Energy_Breakdown(tfseq, J, H)
# dca.Point_mutation_better_binder_checker(tfseq, J, H)
np.set_printoptions(suppress=True)

# results = dca.Weighting_Positions_highEseq(Orig, J, H)
print(results)



# e = dca.Calc_Energy(M2, J, H)
# print(e)

# saffs, ssels = [], []
# affs, seqs = dca.Fasta_Read_Aff(allseqpath)
# for xid, x in enumerate(affs):
#     if x == 1 and dca.Calc_Energy(seqs[xid], J, H) > 55:
#         ssels.append(seqs[xid])
# dca.write_fasta_seqonly(ssels, '/home/jonah/Desktop/Current/targetedbad8.txt')