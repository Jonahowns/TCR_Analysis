import dcamethods as dca
import matplotlib.pyplot as plt
from leven import levenshtein
import numpy as np
from itertools import combinations_with_replacement
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

def import_meme_mat(memefile, motif_states, motif_len):
    o = open(memefile, 'r')
    meme_mat = np.full((motif_len, motif_states), 0.0)
    lid = 0
    for line in o:
        data = line.split()
        for i in range(motif_states):
            meme_mat[lid, i] = float(data[i])
        lid += 1
    o.close()
    return meme_mat

mmat = import_meme_mat(plmp + 'meme_mat.txt', 4, 10)
print(mmat)
nucd = ['A', 'C', 'G', 'T']

def generate_motif_alternatives(meme_mat, motif_states, motif_len):
    diffs, consensus, max_w = [], [], 0
    for i in range(motif_len):
        consensus.append((i, np.argmax(meme_mat[i, :])))
        hw = np.max(meme_mat[i, :])
        max_w += float(hw)
        for j in range(motif_states):
            if (meme_mat[i, j] != 0.0 and (i, j) not in consensus):
                diffs.append((i, j, hw-meme_mat[i, j]))
    motif = ''.join([nucd[j] for (i, j) in consensus])
    return motif

def alt_generator(motif, meme_mat, motif_states, motif_len):
    # poss = permutations(list(np.arange(motif_states)), motif_len)
    poss = list(combinations_with_replacement(''.join([str(x) for x in np.arange(motif_states)]), motif_len))
    print(poss)
    results = []
    for p in poss:
        score = 0
        nz = True
        for xid, x in enumerate(p):
            print(xid, int(x), meme_mat[xid, int(x)])
            if meme_mat[xid, x] != 0:
                score += meme_mat[xid, int(x)]
            # else:
            #     nz = False
        if nz == True:
            results.append((p, score))
    results.sort(key=lambda tup: tup[1])
    print(results)


motif = generate_motif_alternatives(mmat, m_states, m_len)

alt_generator(motif, mmat, m_states, m_len)








'''
def motif_align(seqs, affs, motif, meme_position_specific_probability_matrix, type='dna'):
    if type == 'dna':
        alphabet = dna
    if type == 'rna':
        alphabet = rna
    if type == 'protein':
        alphabet = protein
'''



