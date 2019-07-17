import dcamethods as dca
import matplotlib.pyplot as plt
import numpy as np

AnalysisPath = "C:/Users/Amber/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/FullSeq/Analysis/Comparisons/"
wpath = "C:/Users/Amber/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/FullSeq"
upath = "/home/jonah/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/SeqwAff/"
clusterlist = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

def load_HJ_from_Clust(clustid, **kwargs):
    params = False
    for key, value in kwargs.items():
        if key == 'parameters':
            params = value
    Clustp = wpath + '/Clust' + str(clustid) + '/'
    hp = Clustp + str(clustid) + 'full.h'
    jp = Clustp + str(clustid) + 'full.j'
    H, N, q = dca.sorthmat_plmDCA_autoNandq(hp)
    J = dca.sortjmat_plmDCA(jp, N, q)
    if params:
        return H, J, N, q
    else:
        return H, J

def load_BindersHJ_from_Clust(clustid, **kwargs):
    params = False
    for key, value in kwargs.items():
        if key == 'parameters':
            params = value
    Clustp = upath + '/Clust' + str(clustid) + '/fullbinders/'
    ghp = Clustp + str(clustid) + 'gb.h'
    gjp = Clustp + str(clustid) + 'gb.j'
    bhp = Clustp + str(clustid) + 'bb.h'
    bjp = Clustp + str(clustid) + 'bb.j'
    gH, N, q = dca.sorthmat_plmDCA_autoNandq(ghp)
    gJ = dca.sortjmat_plmDCA(gjp, N, q)
    bH = dca.sorthmat_plmDCA(bhp, N, q)
    bJ = dca.sortjmat_plmDCA(bjp, N, q)
    if params:
        return gH, gJ, bH, bJ, N, q
    else:
        return gH, gJ, bH, bJ


# Returns List of 64 Sublists which contain all seqs energy calculated
def EnergyCoord_COI_ScorePlot_AllClusters(clusterlist):
    MC = CalcE_AllCLuster_wCoord(clusterlist)
    print('All Energies Calculated')
    for i in clusterlist:
        print("Making Figure for Cluster " + str(i))
        specificcoord = [y for x, y in MC if x == i]
        scatpoints = [x for sub in specificcoord for x in sub]
        ScatPlot(scatpoints, i, clusterlist)


def ScatPlot(coord, i, clusterlist):
    fig, ax = plt.subplots()
    ax.scatter(*zip(*coord), c=list(zip(*coord))[0], s=2)
    ax.title.set_text('Cluster ' + str(i) + ' HJ Comparison')
    plt.savefig(AnalysisPath + 'Cluster' + str(i) + '.png', dpi=400)


# Returns Energy Coordinates of all seqs in each cluster based calculated off of one cluster of interests H and J Matrices
# There is a set returned for each clusater taking a turn being the cluster of interest
def CalcE_AllCLuster_wCoord(clusterlist):
    S = load_seqs(clusterlist)
    MasterCoord = []
    for xid, i in enumerate(clusterlist):
        H, J = load_HJ_from_Clust(i)
        tmpclustlist = [x for x in clusterlist if x != i]
        tmpS = [x for x in S if S.index(x) != xid]
        print(tmpclustlist)
        print(tmpS)
        COI = [x for x in S if S.index(x) == xid]
        fullcoord = []
        for cid, clust in enumerate(tmpS):
            for x in clust:
                fullcoord.append((tmpclustlist[cid], dca.Calc_Energy_TCR(x, J, H)))
        for clust in COI:
            for x in clust:
                fullcoord.append((0, dca.Calc_Energy_TCR(x, J, H)))
        MasterCoord.append((i, fullcoord))
    return MasterCoord


def load_seqs(clustlist):
    S = []
    for i in clustlist:
        Clustp = wpath + '/Clust' + str(i) + '/'
        seqp = Clustp + 'full.fasta'
        seqs = dca.Fasta_Read_SeqOnly(seqp)
        S.append(seqs)
    print("All Seqs Loaded")
    return S


#EnergyCoord_COI_ScorePlot_AllClusters(clusterlist)



clustid = 1
gH, gJ, bH, bJ, N, q = load_BindersHJ_from_Clust(clustid, parameters=True)
print(N)
print(q)
bHpos = dca.Sign_Seperator(bH, N, q, mattype='h', sign='+')
bHneg = dca.Sign_Seperator(bH, N, q, mattype='h', sign='-')
gHpos = dca.Sign_Seperator(gH, N, q, mattype='h', sign='+')
gHneg = dca.Sign_Seperator(gH, N, q, mattype='h', sign='-')
bJpos = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='+')
bJneg = dca.Sign_Seperator(bJ, N, q, mattype='j', sign='-')
gJpos = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='+')
gJneg = dca.Sign_Seperator(gJ, N, q, mattype='j', sign='-')




# Rscore 0.18746
J = 2*gJ - bJpos
H = gH + bHneg
# J = gJ + gJneg


# H = np.add(np.subtract(5*gHpos, 3*bHpos), 7*gHneg) #S3
# J = np.add(np.subtract(5*gJpos, 3*bJpos), 7*gJneg) #S3




Clustp = upath + '/Clust' + str(clustid) + '/'
seqfile = Clustp + 'full.fasta'
dca.Raw_wRscore_TCR(J, H, '/home/jonah/Desktop/clust1goodbinder.png', seqfile)

# dca.designerJ(N, q, seqfile, type='pep')