import numpy as np
import dcamethods as dca
from scipy import stats
import copy
import multiprocessing as mp
import math
from sys import getsizeof

aad = {'-': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12,
       'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 20}
nuc_to_id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
base_flip_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U': 4}

class Sequence:
    def __init__(self, N, q, aff, seq='null', cv='no', **kwargs):
        self.type = 'dna'
        self.mat = 'J'
        for key, value in kwargs.items():
            if key == 'type':
                self.type = value
            if key == 'mat':
                self.mat = value
        self.affinity = aff
        if seq == 'null':
            self.seq = np.full(40, ['X'], dtype=str)
        else:
            self.seq = np.array(list(seq))

        self.energyp = None
        self.energyn = None

        if self.mat == 'K':
            # Automatiically Fill in the K Matrix
            self.K = []
            dist = len(self.seq)
            for x in range(dist):
                if self.type == 'dna':
                    ibase = rnad[seq[x]]
                elif self.type == 'pep':
                    ibase = aad[seq[x]]
                for y in range(x + 1, dist):
                    if self.type == 'dna':
                        jbase = rnad[seq[y]]
                    elif self.type == 'pep':
                        jbase = aad[seq[y]]
                    if y <= dist - 2:
                        for z in range(x + 2, dist):
                            if self.type == 'dna':
                                kbase = rnad[seq[z]]
                            elif self.type == 'pep':
                                kbase = aad[seq[z]]
                            self.K.append((x, y-1, z-2, ibase, jbase, kbase))
                            if cv == 'yes' and self.type == 'dna':
                                cibase = rnad[base_flip_rna[seq[x]]]
                                cjbase = rnad[base_flip_rna[seq[y]]]
                                ckbase = rnad[base_flip_rna[seq[z]]]
                                self.K.append((x, y-1, z-2, cibase, cjbase, ckbase))

        elif self.mat == 'J':
            #Automatically Fill in the J Matrix
            self.J = []
            dist = len(self.seq)
            for x in range(dist):
                ibase = aad[seq[x]]
                for y in range(x + 1, dist):
                    jbase = aad[seq[y]]
                    self.J.append((x, y - 1, ibase, jbase))

        elif self.mat == 'H':
            #Automatically Fill in the H Matrix
            self.H = []
            dist = len(self.seq)
            for x in range(dist):
                ibase = aad[seq[x]]
                self.H.append((x, ibase))


    def score(self, ScoringMatrixPos, ScoringMatrixNeg):
        if self.mat == 'K':
            Pos = [x for x in self.K if x in ScoringMatrixPos]
            Neg = [x for x in self.K if x in ScoringMatrixNeg]
        elif self.mat == 'J':
            Pos = [x for x in self.J if x in ScoringMatrixPos]
            Neg = [x for x in self.J if x in ScoringMatrixNeg]
        elif self.mat == 'H':
            Pos = [x for x in self.H if x in ScoringMatrixPos]
            Neg = [x for x in self.H if x in ScoringMatrixNeg]
        return len(Pos)*1, len(Neg)*-1


clustid = 1
droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/SeqwAff/"
cpath = upath + 'Clust' + str(clustid) + '/'
seqp = cpath + 'full.fasta'
titles, seqs = dca.Fasta_Read_Aff(seqp)
N = 40
q = 5

SEQHANDLER = []
for i in range(len(seqs)):
    x = Sequence(N, q, titles[i], seqs[i], type='dna', mat='K')
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


def Node_Worker(id, N, q, i, bj, ej, SEQHANDLER):
    ScoringMatrixPos = []
    ScoringMatrixNeg = []
    cR = -1
    for j in range(bj, ej):
        for k in range(N-2):
            for x in range(1, q):
                for y in range(1, q):
                    for z in range(1, q):
                        if i > j or j > k:
                            continue
                        if j == math.floor((ej-bj)/2 + bj) and k == N-2 and z == 1:
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


def Write_Output_Node(i, z):
    nPos = []
    nNeg = []
    rs = []
    for x in z:
        for pos, neg, cR in x:
            for p in pos:
                nPos.append(p)
            for n in neg:
                nNeg.append(n)
            rs.append((i, cR))
    return nPos, nNeg, rs

def Write_Overall_Output(Pos, Neg, rs, outk, outrs):
    K = np.full((N - 2, N - 2, N - 2, q, q, q), 0)
    for d in Pos:
        x, y, z, i, j, k = d
        K[x, y, z, i, j, k] = 1
    for d in Neg:
        x, y, z, i, j, k = d
        K[x, y, z, i, j, k] = -1
    o = open(outrs, 'w')
    for d in rs:
        i, rs = d
        print(i, file=o)
        print(rs, file=o)
    o.close()
    dca.export_K(K, N, q, outk)


def gen_indices(bPni):
    ind = []
    for x in bPni:
        sep = math.floor((N - 2) / x)
        for i in range(x):
            if i == 0:
                ind.append(0)
                ind.append(sep)
            elif i == x-1:
                ind.append(ind[-1])
                ind.append(N-2)
            else:
                ind.append(ind[-1])
                ind.append(ind[-1]+sep)
    it = iter(ind)
    tup = list(zip(it, it))
    return tup



def Parallelized_ScoringK(N, q, SEQHANDLER, CoreNum, outk, outr):
    bpn = math.floor(CoreNum / (N - 2))
    lo = CoreNum - (N - 2) * bpn
    bPni = [(bpn + 1) if x < lo else bpn for x in range(0, N - 2)]
    ind = gen_indices(bPni)
    Pos, Neg, r = [], [], []
    for i in range(N-2):
        instructions = []
        for j in range(CoreNum):
            instructions.append([i, N, q, i, ind[j][0], ind[j][1], SEQHANDLER])
        p = mp.Pool(CoreNum)
        z = p.starmap(Node_Worker, instructions)
        nPos, nNeg, rs = Write_Output_Node(i, z)
        Pos += nPos
        Neg += nNeg
        r += rs
    Write_Overall_Output(Pos, Neg, r, outk, outr)



Parallelized_ScoringK(N, q, SEQHANDLER, 4,  '/home/jonah/Desktop/DesignerK.txt')






