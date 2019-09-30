#!/usr/bin python


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import gaussian_kde
import sys
import copy
import math
import dcamethods as dca
from multiprocessing import Process, Queue

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
tenpath = cpath + 'v4/10percent/'

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

Jsel = dca.TopJNorms_Jmatrix(j901, N, q, 332)[0]
H = h901

allseqpath = cpath + str(fam) + 'all.txt'
s10p = tenpath + 's10.txt'
s90p = tenpath + 's90.txt'



eH = np.full((N, q), 0.0)
# fig, ax = plt.subplots(1,2)
# dca. Raw_wRscore(Jsel, 2*H, '/home/jonah/Desktop/Current/v4/n332.png', allseqpath)
# dca.Diff_Raw_wRscore(Jsel, 2*H, '/home/jonah/Desktop/Current/v4/n3329010.png', ['all', 'training', '10%'], allseqpath, s90p, s10p)
dca.Diff_Avg(Jsel, 2*H, '/home/jonah/Desktop/Current/v4/n3329010avg.png', ['all',  'training', '10%'], allseqpath, s90p, s10p)
# dca.Raw_wRscore_subplot(ax[0], Jsel, 2*H, allseqpath)
# dca.Raw_wRscore_subplot(ax[1], Jsel, 2*H, s10p)
# plt.savefig('/home/jonah/Desktop/normtrial.png', dpi=400)


# H = (2*(gH) - bH) # S1 x2
# J = (2*(gJ) - bJ) # S1 x0.5
# z = Queue()
# p1 = Process(target=dca.Jnorm_finder, args=(j900, N, q, 50, 300, 50, allseqpath, z))
# p2 = Process(target=dca.Jnorm_finder, args=(j900, N, q, 300, 550, 50, allseqpath, z))
# p1.start()
# p2.start()
# results = []
# results += z.get()
# results += z.get()
# print(results)

