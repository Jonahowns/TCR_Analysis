#!/usr/bin python

# SCRIPT TO SAMPLE FINAL ENERGY LANDSCAPE
import dcamethods as dca
import numpy as np
import multiprocessing as mp

cpath = "/home/jprocyk/mc_sampling/"



fam = 8
N = 40
q = 5
gHp = cpath + str(fam) + 'go.h'
gJp = cpath + str(fam) + 'go.j'
gH = dca.sorthmat_plmDCA(gHp, N, q)
gJ = dca.sortjmat_plmDCA(gJp, N, q)

allseqpath = cpath + str(fam) + 'all.txt'

Jsel, vals, cut = dca.TopJNorms_Jmatrix(gJ, N, q, 305)
Hsel, hvals, hcut = dca.TopHNorms_Hmatrix(gH, N, q, 5)


Ts = np.arange(0.2, 0.55, 0.05)
muts = np.arange(6, 20, 3)
out = []
outdir = '/home/jprocyk/mc_sampling/8v2seqs/'
for i in range(22):
    if i < 10:
        out.append(outdir + str(i) + 'gb.txt')
    elif i < 20:
        out.append(outdir + str(i) + 'gba.txt')
    else:
        out.append(outdir + str(i) + 'bb.txt')

jobs = []
for i in range(22):
    jobs.append(dca.GenerSeq(N, Ts[i%6], mut_steps=muts[i%3], out_after=10000, steps=10000000))

procs = []
for xid, x in enumerate(jobs):
    if xid < 10:
        procs.append(mp.Process(target=jobs[xid].run_sampling, args=(Jsel, 2*gH, out[xid])))
    elif xid < 20:
        procs.append(mp.Process(target=jobs[xid].run_adaptive_sampling, args=(Jsel, 2*gH, out[xid])))
    else:
        procs.append(mp.Process(target=jobs[xid].run_bad_sampling, args=(Jsel, 2*gH, out[xid])))


for p in procs:
    p.start()
for p in procs:
    p.join()




