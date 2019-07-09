import numpy as np
import dcamethods as dca
from scipy import stats
from sys import getsizeof


nuc_to_id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
base_flip_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U': 4}

class Sequence:
    def __init__(self, N, q, aff, seq='null', cv='no'):
        self.affinity = aff
        self.K = np.full((N-2, N-2, N-2, q, q, q), 0.0)
        if seq == 'null':
            self.seq = np.full(40, ['X'], dtype=str)
        else:
            self.seq = np.array(list(seq))

        self.score = None
        # Automatiically Fill in the K Matrix
        dist = len(self.seq)
        for x in range(dist):
            ibase = rnad[seq[x]]
            for y in range(x + 1, dist):
                jbase = rnad[seq[y]]
                if y <= dist - 2:
                    for z in range(x + 2, dist):
                        kbase = rnad[seq[z]]
                        self.K[x, y - 1, z - 2, ibase, jbase, kbase] = 1
                        if cv == 'yes':
                            cibase = rnad[base_flip_rna[seq[x]]]
                            cjbase = rnad[base_flip_rna[seq[y]]]
                            ckbase = rnad[base_flip_rna[seq[z]]]
                            self.K[x, y - 1, z - 2, cibase, cjbase, ckbase] = 1

    def score(self, ScoringMatrix):
        self.score = np.sum(np.multiply(self.K, ScoringMatrix))


droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
testpath = upath + droppath
testseqpath = testpath + '8thfull.txt'
titles, seqs = dca.Fasta_Read_Aff(testseqpath)
N = 40
q = 5



# SEQHANDLER = []
# x = Sequence(N, q, titles[0], seqs[0])
# SEQHANDLER.append(x)
# getsizeof(x)
# x = Sequence(N, q, titles[1], seqs[1])
# SEQHANDLER.append(x)

SEQHANDLER = []
for i in range(len(seqs)):
    x = Sequence(N, q, titles[i], seqs[i])
    SEQHANDLER.append(x)

currentR = -1
aff = list(set(titles))
ScoringMatrix = np.full((N-2, N-2, N-2, q, q, q), 0.0)
for i in range(N-2):
    for j in range(N - 2):
        for k in range(N - 2):
            for x in range(1, q):
                for y in range(1, q):
                    for z in range(1, q):
                        if i > j or j > k:
                            continue
                        ScoringMatrix[i, j, k, x, y, z] = 1
                        rscore = dca.R_SCORE_w_SEQHANDLER(SEQHANDLER, ScoringMatrix, titles)
                        if rscore > currentR:
                            currentR = rscore
                        else:
                            ScoringMatrix[i, j, k, x, y, z] = -1
                            rscore = dca.R_SCORE_w_SEQHANDLER(SEQHANDLER, ScoringMatrix, titles)
                            if rscore > currentR:
                                currentR = rscore
                            else:
                                ScoringMatrix[i, j, k, x, y, z] = 0.0
print('designed K Matrix R score: ' + str(currentR))
dca.export_K(ScoringMatrix, N, q, '/home/jonah/Desktop/DesignerK.txt')