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

# Goodbinders
gHp = v3path + str(fam) + 'go.h'
gJp = v3path + str(fam) + 'go.j'

gH = dca.sorthmat_plmDCA(gHp, N, q)
gJ = dca.sortjmat_plmDCA(gJp, N, q)


allseqpath = cpath + str(fam) + 'all.txt'


eH = np.full((N, q), 0.0)
# H = (2*(gH) - bH) # S1 x2
# J = (2*(gJ) - bJ) # S1 x0.5
z = Queue()
p1 = Process(target=dca.Jnorm_finder, args=(gJ, N, q, 50, 300, 50, allseqpath, z))
p2 = Process(target=dca.Jnorm_finder, args=(gJ, N, q, 300, 550, 50, allseqpath, z))
p1.start()
p2.start()
results = []
results += z.get()
results += z.get()
print(results)

