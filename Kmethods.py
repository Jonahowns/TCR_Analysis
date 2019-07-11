import numpy as np
import dcamethods as dca
from scipy import stats
import copy
import multiprocessing as mp
import math
from sys import getsizeof


nuc_to_id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
base_flip_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U': 4}

class Sequence:
    def __init__(self, N, q, aff, seq='null', cv='no'):
        self.affinity = aff
        self.K = []
        if seq == 'null':
            self.seq = np.full(40, ['X'], dtype=str)
        else:
            self.seq = np.array(list(seq))

        self.energyp = None
        self.energyn = None
        # Automatiically Fill in the K Matrix
        dist = len(self.seq)
        for x in range(dist):
            ibase = rnad[seq[x]]
            for y in range(x + 1, dist):
                jbase = rnad[seq[y]]
                if y <= dist - 2:
                    for z in range(x + 2, dist):
                        kbase = rnad[seq[z]]
                        self.K.append((x, y-1, z-2, ibase, jbase, kbase))
                        if cv == 'yes':
                            cibase = rnad[base_flip_rna[seq[x]]]
                            cjbase = rnad[base_flip_rna[seq[y]]]
                            ckbase = rnad[base_flip_rna[seq[z]]]
                            self.K.append((x, y-1, z-2, cibase, cjbase, ckbase))

    def score(self, ScoringMatrixPos, ScoringMatrixNeg):
        Pos = [x for x in self.K if x in ScoringMatrixPos]
        Neg = [x for x in self.K if x in ScoringMatrixNeg]
        return len(Pos)*1, len(Neg)*-1


    def tmpK(self):
        tmpK = np.zeros((N-2, N-2, N-2, q, q, q), dtype='int')
        for entry in self.K:
            x, y, z, i, j, k, val = entry
            tmpK[x, y, z, i, j, k] = val
        return tmpK


droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
testpath = upath + droppath
testseqpath = testpath + '8thfull.txt'
titles, seqs = dca.Fasta_Read_Aff(testseqpath)
N = 40
q = 5

SEQHANDLER = []
for i in range(len(seqs)):
    x = Sequence(N, q, titles[i], seqs[i])
    SEQHANDLER.append(x)


def R_SCORE_OPTIMIZED(SEQHANDLER, ScoringMatrixPos, ScoringMatrixNeg):
    pE = []
    nE = []
    As = []
    for x in SEQHANDLER:
        pos, neg = x.score(ScoringMatrixPos, ScoringMatrixNeg)
        pE.append((x.affinity, pos))
        nE.append((x.affinity, neg))
        As.append(x.affinity)
    affs = list(set(As))
    datax = []
    datap = []
    datan = []
    for a in affs:
        maxpos = max([y for (x, y) in pE if x == a])
        maxneg = min([y for (x, y) in nE if x == a])
        datax.append(a)
        datap.append(maxpos)
        datan.append(maxneg)
    linregpos = stats.linregress(datax, datap)
    linregneg = stats.linregress(datax, datan)
    rscorepos = linregpos[2]
    rscoreneg = linregneg[2]
    return rscorepos, rscoreneg


def Worker(id, N, q, bi, ei, SEQHANDLER):
    ScoringMatrixPos = []
    ScoringMatrixNeg = []
    cR = -1
    for i in range(bi, ei):
        for j in range(N-2):
            for k in range(N-2):
                for x in range(1, q):
                    for y in range(1, q):
                        for z in range(1, q):
                            if i > j or j > k:
                                continue
                            if i == math.floor((ei-bi)/2 + bi):
                                print('50% on Worker ' + str(id))
                            ScoringMatrixPos.append((i, j, k, x, y, z))
                            ScoringMatrixNeg .append((i, j, k, x, y, z))
                            rscorepos, rscoreneg = R_SCORE_OPTIMIZED(SEQHANDLER, ScoringMatrixPos, ScoringMatrixNeg)
                            if rscoreneg < cR and rscorepos < cR:
                                ScoringMatrixNeg.remove((i, j, k, x, y, z))
                                ScoringMatrixPos.remove((i, j, k, x, y, z))
                            if rscorepos > rscoreneg and rscorepos > cR:
                                cR = rscorepos
                                ScoringMatrixNeg.remove((i, j, k, x, y, z))
                            if rscoreneg > rscorepos and rscoreneg > cR:
                                cR = rscoreneg
                                ScoringMatrixPos.remove((i, j, k, x, y, z))
    print('Worker ' + str(id))
    print('Optimal R Score Obtained: ' + str(cR))
    return ScoringMatrixPos, ScoringMatrixNeg, cR


def Write_Output(z, N, q, outfile):
    aPos = []
    aNeg = []
    rs = []
    K = np.full((N-2, N-2, N-2, q, q, q), 0)
    for i in z:
        for pos, neg, cR in i:
            for p in pos:
                aPos.append(p)
            for n in neg:
                aNeg.append(n)
            rs.append(cR)
    for d in aPos:
        x, y, z, i, j, k  = d
        K[x, y, z, i, j, k] = 1
    for d in aNeg:
        x, y, z, i, j, k  = d
        K[x, y, z, i, j, k] = -1
    dca.export_K(K, N, q, outfile)
    return rs

def Parallelized_ScoringK(N, q, SEQHANDLER, CoreNum, outfile):
    comm = round((N - 2) / CoreNum)
    ind = [x * comm for x in range(CoreNum)]
    ind.append(N - 2)
    instructions = []
    for i in range(CoreNum):
        instructions.append([i, N, q, ind[i], ind[i+1], SEQHANDLER])
    p = mp.Pool(CoreNum)
    z = p.starmap(Worker, instructions)
    r = Write_Output(z, N, q, outfile)
    print(r)




Parallelized_ScoringK(N, q, SEQHANDLER, 4,  '/home/jonah/Desktop/DesignerK.txt')










