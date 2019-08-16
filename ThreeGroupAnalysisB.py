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

droppath = "Projects/DCA/v2/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
testpath = upath + 'Projects/DCA/GenSeqs/'
cpath = "/home/jonah/Current/"
bldrop = "bl78/"
g2path = upath + droppath + 'FamHJ/'
g3path = upath + droppath + '3GHJ/'

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


fam = 8
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
allseqpath = cpath + str(fam) + 'all.txt'


print(np.sum(vJ))
print(np.sum(gJ))
print(np.sum(bJ))


H = (vH + gH - bH) # S1
J = (vJ + gJ - bJ) # S1

# norms = dca.TopX_JNorms(J, N, 10)
# print(norms)

dca.Raw_Aff_v_E(J, H, ('A vs E Family ' + str(fam)), ('/home/jonah/Desktop/' + str(fam) + '3Gtrial.png'), allseqpath)
# dca.Raw_wRscore(J, H, ('/home/jonah/Desktop/' + str(fam) + 'normfam7.png'), testseqpath)
