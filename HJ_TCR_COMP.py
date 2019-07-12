import dcamethods as dca
import matplotlib.pyplot as plt
import numpy as np

AnalysisPath = "C:/Users/Amber/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/FullSeq/Analysis/Comparisons/"
wpath = "C:/Users/Amber/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/FullSeq"
clusterlist = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

def load_HJ_from_Clust(clustid):
    Clustp = wpath + '/Clust' + str(clustid) + '/'
    hp = Clustp + str(clustid) + 'full.h'
    jp = Clustp + str(clustid) + 'full.j'
    H, N, q = dca.sorthmat_plmDCA_autoNandq(hp)
    J = dca.sortjmat_plmDCA(jp, N, q)
    return H, J


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


EnergyCoord_COI_ScorePlot_AllClusters(clusterlist)

