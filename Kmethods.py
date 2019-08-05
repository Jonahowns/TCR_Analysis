import numpy as np
import dcamethods as dca
from scipy import stats
import math
from mpi4py import MPI

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


def Node_Worker(rank, N, q, i, j, bk, ek, SEQHANDLER):
    ScoringMatrixPos = []
    ScoringMatrixNeg = []
    cR = -1
    for k in range(bk, ek):
        if k == (ek - 1):
            print('Worker', rank, 'Node', i, j, 'k = ', k)
        for x in range(1, q):
            for y in range(1, q):
                for z in range(1, q):
                    print('Worker ', rank, 'z= ', z)
                    if i > j or j > k:
                        continue
                    if j == math.floor((ek-bk)/2 + bk) and x == 1 and y == 1 and z == 1:
                        print("50% on Worker " + str(rank))
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
    print('Worker ' + str(rank))
    print('Optimal R Score Obtained: ' + str(cR))
    return ScoringMatrixPos, ScoringMatrixNeg, cR


def Write_Overall_Output(Pos, Neg, rs, outk, outrs, N, q):
    K = np.full((N - 2, N - 2, N - 2, q, q, q), 0)
    for npos in Pos:
        for d in npos:
            x, y, z, i, j, k = d
            K[x, y, z, i, j, k] = 1
    for nneg in Neg:
        for d in nneg:
            x, y, z, i, j, k = d
            K[x, y, z, i, j, k] = -1
    o = open(outrs, 'w')
    for nrs in rs:
        for d in nrs:
            r = d
            print(r, file=o)
    o.close()
    dca.export_K(K, N, q, outk)


def gen_indices(bPni, i):
    ind = []
    for xid, x in enumerate(bPni):
        sep = math.floor((N - 2) / x)
        for y in range(x):
            if y == 0:
                ind.append(xid+i)
                ind.append(0)
                ind.append(sep)
            elif y == x-1:
                ind.append(ind[-1])
                ind.insert(-1, xid+i)
                ind.append(N-3)
            else:
                ind.append(ind[-1])
                ind.insert(-1, xid+i)
                ind.append(ind[-1]+sep)
    it = iter(ind)
    tup = list(zip(it, it, it))
    return tup


# clustid = 1
# droppath = "Projects/DCA/GenSeqs/"
# fullpath = macpath+droppath
# upath = "/home/jonah/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/SeqwAff/"
# cpath = upath + 'Clust' + str(clustid) + '/'
# seqp = cpath + 'full.fasta'


if __name__ == "__main__":
    N = 40
    q = 5

    comm = MPI.COMM_WORLD
    rank = comm.rank
    CoreNum = comm.size

    print('Corenum: ', CoreNum)
    print('My rank is ', rank)

    seqfile = "/home/jprocyk/Kmat/8thfull.txt"
    titles, seqs = dca.Fasta_Read_Aff(seqfile)
    SEQHANDLER = []
    for i in range(len(seqs)):
        x = Sequence(N, q, titles[i], seqs[i], type='dna', mat='K')
        SEQHANDLER.append(x)


    PosContrib = []
    NegContrib = []
    rscores = []
    instructions = []
    for i in range(N-2):
        if rank == 0:
            bpn = math.floor(CoreNum/(N-2-i))
            if bpn > N-2:
                bpn = 38
                bPni = [bpn for i in range(N-2-i)]
                ind = gen_indices(bPni, i)
                used = len(bPni)*bpn
            else:
                lo = CoreNum - (N-2-i)*bpn
                bPni = [(bpn+1) if x < lo else bpn for x in range(0, N-2-i)]
                ind = gen_indices(bPni, i)
                used = CoreNum
            instructions = []
            for j in range(used):
                instructions.append([j+1, N, q, i, ind[j][0], ind[j][1], ind[j][2]])
        instructions = comm.bcast(instructions, root=0)
        Pos, Neg, r = [], [], []
        if rank != 0:
            try:
                si  = [(rk, N, q, i, j, bk, ek) for (rk, N, q, i, j, bk, ek) in instructions if rk == rank][0]
                Pos, Neg, r = Node_Worker(si[0], si[1], si[2], si[3], si[4], si[5], si[6], SEQHANDLER)
            except IndexError:
                print('Worker ', rank, 'not used i =',i)
        allpos = comm.gather(Pos, root = 0)
        allneg = comm.gather(Neg, root = 0)
        allrs = comm.gather(r, root = 0)
        if rank == 0:
            PosContrib += allpos
            NegContrib += allneg
            rscores += r

    if rank == 0:
        Write_Overall_Output(PosContrib, NegContrib, rscores, '/home/jprocyk/Kmat/DesignerK.txt', '/home/jprocyk/Kmat/rscores.txt', N, q)
        print('done')