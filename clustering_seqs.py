import dcamethods as dca
import matplotlib.pyplot as plt
from leven import levenshtein
import numpy as np
import itertools as it
from sklearn.cluster import DBSCAN
from sklearn import preprocessing
from sklearn import metrics
import sklearn

upath = "/home/jonah/Dropbox (ASU)/"
datap = "Projects/DCA/GunterAptamers/Selex_Data/"
plmA = 'plmAnalysis/'
computer_path = upath

plmp = computer_path + datap

fasta = plmp + 'r15_train_c100.txt'
affs, seqs = dca.Fasta_Read_Aff(fasta)

# One Hot Encoding Sequences
# enc = preprocessing.OneHotEncoder()
# X_all = []
# for x in seqs:
#     X_all.append(list(x))
# enc.fit(X_all)
# onehotlabels = enc.transform(X_all).toarray()


'''
def lev_metric(x, y):
    i, j = int(x[0]), int(y[0])     # extract indices
    return levenshtein(data[i], data[j])
'''
# X = np.arange(len(data)).reshape(-1, 1)
# X = onehotlabels
# print(X[0])
# db = DBSCAN(eps=0.4, metric='dice', min_samples=15).fit(X)
# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
# core_samples_mask[db.core_sample_indices_] = True
# labels = db.labels_
#
# dis_mat = sklearn.metrics.pairwise_distances(X, metric='dice')
# print(dis_mat)
# print(len(seqs), len(labels))
#
# # Number of clusters in labels, ignoring noise if present.
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# n_noise_ = list(labels).count(-1)
# for clust in set(labels):
#     print('Clust', clust, 'Length', list(labels).count(clust))
#
# print('Estimated number of clusters: %d' % n_clusters_)
# print('Estimated number of noise points: %d' % n_noise_)
# if(len(set(labels)) > 2):
#     print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels))
#



m_len, m_states = 10, 4

class Motif_Aligner():
    def __init__(self, memefile, motif_states, motif_len, alphabet='dna', gaps=False):
        dna = ['A', 'C', 'G', 'T']
        dnag = ['-', 'A', 'C', 'G', 'T']
        if alphabet == 'dna' and not gaps:
            self.alpha = dna
        elif alphabet == 'dna' and gaps:
            self.alpha = dnag
        self.mfile = memefile
        self.m_states = motif_states
        self.m_len = motif_len
        self.m_mat = self.import_meme_mat()
        self.results = self.alt_generator()
        self.all_motifs = self.decode_motifs()
        self.unaligned_seqs = []
        self.unaligned_indxs = []
        self.affinities = []
        self.aligned_seqs = None
        self.start_positions = []
    def import_meme_mat(self):
        o = open(self.mfile, 'r')
        meme_mat = np.full((self.m_len, self.m_states), 0.0)
        lid = 0
        for line in o:
            data = line.split()
            for i in range(self.m_states):
                meme_mat[lid, i] = float(data[i])
            lid += 1
        o.close()
        return meme_mat
    def alt_generator(self):
        poss = list(it.product(np.arange(self.m_states), repeat=self.m_len))
        results = []
        for p in poss:
            score = 0
            nz = True
            for xid, x in enumerate(p):
                if self.m_mat[xid, int(x)] != 0:
                    score += self.m_mat[xid, int(x)]
                else:
                    nz = False
                    break
            if nz == True:
                results.append((p, score))
        results.sort(key=lambda tup: tup[1])
        results.reverse()
        return results
    def decode_motifs(self):
        ams = []
        for c, w in self.results:
            l = []
            for x in c:
                l.append(self.alpha[int(x)])
            ams.append(''.join(l))
        return ams
    def load_seqs(self, seqs, affs):
        ds = []
        for sid, s in enumerate(seqs):
            found = False
            for x in self.all_motifs:
                indx = s.find(x)
                if indx != -1:
                    found = True
                    ds.append(indx)
                    break
            if not found:
                ds.append(-1)
        for iid, idx in enumerate(ds):
            if idx > -1:
                self.affinities.append(affs[iid])
                self.unaligned_seqs.append(seqs[iid])
                self.unaligned_indxs.append(idx)
        print('No Motif Found in', len(seqs) - len(self.unaligned_seqs), 'of', len(seqs), 'sequences')
    def align_seqs(self):
        g_indx = max(self.unaligned_indxs)
        self.aligned_seqs = ['x']*len(self.unaligned_seqs)
        print(len(self.aligned_seqs))
        ls = []
        for uid, ua in enumerate(self.unaligned_indxs):
            gaps_to_add = g_indx - ua
            self.start_positions.append(gaps_to_add)
            if gaps_to_add > 0:
                ns = ''.join(['-' for x in range(gaps_to_add)]) + self.unaligned_seqs[uid]
                ls.append(len(ns))
                self.aligned_seqs[uid] = ns
            else:
                self.aligned_seqs[uid] = self.unaligned_seqs[uid]
        print(self.aligned_seqs)
        maxl = max(ls)
        for aid, al in enumerate(self.aligned_seqs):
            print(al)
            end_gaps = maxl - len(al)
            if end_gaps > 0:
                ns = self.aligned_seqs[aid] + ''.join(['-' for x in range(end_gaps)])
                self.aligned_seqs[aid] = ns
    def write_msa(self, outfile):
        o = open(outfile, 'w')
        for iid, i in enumerate(self.aligned_seqs):
            print('>Seq' + str(iid) + '-' + str(self.affinities[iid]), file=o)
            print(i, file=o)
        o.close()






w1 = Motif_Aligner(plmp + 'meme_mat.txt', 4, 10)
print(w1.all_motifs)
w1.load_seqs(seqs, affs)
print(len(w1.unaligned_indxs))
print(len(w1.unaligned_seqs))
w1.align_seqs()
w1.write_msa('/home/jonah/Desktop/tmp/msatrial.txt')
print(w1.aligned_seqs)
print(w1.start_positions)

start = w1.start_positions
sp = [[x] for x in start]
print(sp)


db = DBSCAN(eps=2, metric='euclidean', min_samples=100).fit(sp)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

dis_mat = sklearn.metrics.pairwise_distances(sp, metric='euclidean')
print(dis_mat)
print(len(seqs), len(labels))

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
for clust in set(labels):
    print('Clust', clust, 'Length', list(labels).count(clust))

print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of noise points: %d' % n_noise_)
sqs = w1.aligned_seqs
afs = w1.affinities
c0, c1 = [], []
for xid, x in enumerate(db.labels_):
    if x == 0:
        c0.append((sqs[xid], afs[xid]))
    elif x == 1:
        c1.append((sqs[xid], afs[xid]))
    else:
        continue

c1s, c1a = zip(*c0)
c0s, c0a = zip(*c1)


dca.write_fasta_aff(c1s, c1a, plmp + 'c100_train_c1.txt')
dca.write_fasta_aff(c0s, c0a, plmp + 'c100_train_c0.txt')
#print(c1s[0], c1a[0])



'''
def motif_align(seqs, affs, motif, meme_position_specific_probability_matrix, type='dna'):
    if type == 'dna':
        alphabet = dna
    if type == 'rna':
        alphabet = rna
    if type == 'protein':
        alphabet = protein
'''



