import dcamethods as dca
import matplotlib.pyplot as plt
import numpy as np

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


def ScatPlot(coord, i):
    fig, ax = plt.subplots()
    ax.scatter(*zip(*coord), c=list(zip(*coord))[0], s=2)
    ax.title.set_text('Cluster ' + str(i) + ' HJ Comparison')
    plt.savefig(AnalysisPath + 'Cluster' + str(i) + '.png', dpi=400)
    plt.close()


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

##################################################
#Data preprocessing for TCR DATA
def prune_alignment_peptides(seqs, simt=0.99):
    final_choice_seqs = []
    for sid, seq in enumerate(seqs):
        append = True
        seqoi = list(seq)
        for existseq in final_choice_seqs:
            es = list(existseq)
            seq_similarity = 0.
            for i in range(len(es)):
                if seqoi[i] == es[i]:
                    seq_similarity += 1.
            seq_similarity /= len(seq)
            if seq_similarity >= simt:
                append = False
                break
        if append:
            final_choice_seqs.append(seq.upper())
    print('INFO: reduced length of alignment from %d to %d due to sequence similarity' % (
    len(seqs), len(final_choice_seqs)), file=sys.stderr),
    return final_choice_seqs


def common_seq_checker(*argv):
    # Take in Lists of Seqs and compares all lists for common sequences #
    for arg1id, arg1 in enumerate(argv):
        for arg2id, arg2 in enumerate(argv):
            if arg1id == arg2id or arg2id < arg1id:
                continue
            full = arg1 + arg2
            uniq = set(full)
            print('Set ' + str(arg1id) + ' and Set ' + str(arg2id) + ' TotalSeqs: ' + str(len(full)))
            print('CommonSeq = ' + str(2*(len(full)-len(uniq))))


def lengthchecker(inputfile, mode='check'):
    o = open(inputfile, 'r')
    seqs = []
    seql = []
    for line in o:
        if line.startswith('>'):
            continue
        else:
            seqs.append(line)
            seql.append(len(list(line)))
    nseql = np.array(seql)
    mostcommonlen = int(stats.mode(nseql))
    lencatcher = 0
    print('Total Sequences = ' + str(len(seqs)) + '\n')
    for xid, x in enumerate(nseql):
        if x != mostcommonlen:
            if mode == 'fix':
                del seqs[xid-lencatcher]
            lencatcher += 1
    o.close()
    print('Not exactly ' + str(mostcommonlen) + ' bases, percentage = ' + str(lencatcher / len(nseql) * 100) + '\n')
    if mode == 'fix':
        print('Sequences Removed\n')
    return seqs


def writefasta(seqs, file, clustid):
    y = open(file, 'w+')
    count = 1
    for x in seqs:
        print('>' + 'seq' + str(count) + '-' + str(clustid), file=y)
        print(x, file=y)
        count += 1
    y.close()


if __name__ == '__main__':

    # Variables
    clusterid = 5
    macpath = "/Users/Amber/Dropbox (ASU)/"
    droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/Trial/"
    clusterpath = 'Cluster' + str(clusterid) + '/'
    fullpath = macpath + droppath + clusterpath
    AAA = fullpath + 'AAA' + str(clusterid) + '.fasta'
    ACC = fullpath + 'ACC' + str(clusterid) + '.fasta'
    AGG = fullpath + 'AGG' + str(clusterid) + '.fasta'
    comb = fullpath + 'full' + str(clusterid) + '.fasta'

    aaa = lengthchecker(AAA)
    acc = lengthchecker(ACC)
    agg = lengthchecker(AGG)
    common_seq_checker(aaa, acc, agg)

    #all = lengthchecker(comb)
    #aggp = prune_alignment_peptides(all, 0.99)



####################################################
#Data Import
clustids = {0: 'AAA', 1: 'ACC', 2: 'AGG', 3: 'ATT',
            4: 'CAC', 5: 'CCG', 6: 'CTA',
            7: 'GAG', 8: 'GCT', 9: 'GTC',
            10: 'TAT', 11: 'TCA', 12: 'TGC', 13: 'TTG'}


def seqpercluster(*argv):
    commclust={}
    argnum=len(argv)
    shared={}
    excluded = []
    for arg in argv:
        clustid = sorted(arg['cluster_id'].unique())
        for x in clustid:
            if x not in commclust:
                commclust[x] = 1
            else:
                commclust[x] += 1
    for x, y in commclust.items():
        if y == argnum:
            shared[x] = 0
        else:
            excluded.append(x)
    for arg in argv:
        for x,y in shared.items():
            subclust = arg[arg['cluster_id'] == x]
            shared[x] += len(subclust)
    print('Excluded Clusters : ' + str(excluded))
    return shared


def outfilegen(clustid, types, dataset=0):
    if types == 'i':
        clustname = clustids[dataset]
        outfile = fullpath + str(clustname) + '/' + str(clustname) + str(clustid) + '.fasta'
    if types == 'f':
        outfile = fullpath + 'full' + str(clustid) + '.fasta'
    return outfile


def writefasta(seqs, file, clustid):
    y=open(file, 'w+')
    count = 1
    for x in seqs:
        print('>'+'seq' + str(count) + '-' + str(clustid), file=y)
        print(x, file=y)
        count += 1
    y.close()


def extractseq(clustlist, *argv):
    for x in clustlist:
        fullclustout = outfilegen(x, 'f')
        fullclust = []
        for cid, arg in enumerate(argv):
            seqonly = []
            subclust = arg[arg['cluster_id'] == x]
            for index, row in subclust.iterrows():
                seqonly.append(row['aminoAcid'])
                fullclust.append(row['aminoAcid'])
            subclustout = outfilegen(x, 'i', cid)
            writefasta(seqonly, subclustout, x)
        writefasta(fullclust, fullclustout, x)
    print("All fasta files written to " + fullpath)

def common_seq_checker(*argv):
    # Take in Lists of Seqs and compares all lists for common sequences #
    for arg1id, arg1 in enumerate(argv):
        for arg2id, arg2 in enumerate(argv):
            if arg1id == arg2id or arg2id < arg1id:
                continue
            full = arg1 + arg2
            uniq = set(full)
            print('Set ' + str(arg1id) + ' and Set ' + str(arg2id) + ' TotalSeqs: ' + str(len(full)))
            print('CommonSeq = ' + str(2*(len(full)-len(uniq))))


########################################################
#TCR CLUSTER STUFF
def create_figure(clustid):
    clustpath = macpath + droppath + 'FullSeq/Clust' + str(clustid) + '/'
    analysispath = macpath + droppath + 'FullSeq/Analysis/'
    afasta = clustpath + 'afam.fasta'
    # SeqLogos
    alogo = clustpath + 'asl.png'
    clogo = clustpath + 'csl.png'
    glogo = clustpath + 'gsl.png'
    tlogo = clustpath + 'tsl.png'
    fulllogo = clustpath + 'fsl.png'
    # Import Photos
    asl = mpimg.imread(alogo)
    csl = mpimg.imread(clogo)
    gsl = mpimg.imread(glogo)
    tsl = mpimg.imread(tlogo)
    fsl = mpimg.imread(fulllogo)
    # Matrix Paths
    aji = clustpath + str(clustid) + 'a.j'
    ahi = clustpath + str(clustid) + 'a.h'
    cji = clustpath + str(clustid) + 'c.j'
    chi = clustpath + str(clustid) + 'c.h'
    gji = clustpath + str(clustid) + 'g.j'
    ghi = clustpath + str(clustid) + 'g.h'
    tji = clustpath + str(clustid) + 't.j'
    thi = clustpath + str(clustid) + 't.h'
    fullji = clustpath + str(clustid) + 'full.j'
    fullhi = clustpath + str(clustid) + 'full.h'
    n = dca.getn(afasta)
    aj = dca.sortjmat(aji, n, 21)
    ah = sorthmat(ahi, n, 21)
    ajdisp = jdisplay(aj, n, 21)
    cj = sortjmat(cji, n, 21)
    ch = sorthmat(chi, n, 21)
    cjdisp = jdisplay(cj, n, 21)
    gj = sortjmat(gji, n, 21)
    gh = sorthmat(ghi, n, 21)
    gjdisp = jdisplay(gj, n, 21)
    tj = sortjmat(tji, n, 21)
    th = sorthmat(thi, n, 21)
    tjdisp = jdisplay(tj, n, 21)
    fullj = sortjmat(fullji, n, 21)
    fullh = sorthmat(fullhi, n, 21)
    fulljdisp = jdisplay(fullj, n, 21)

    # Plotting
    fig, ax = plt.subplots(5, 3, figsize=(12, 22), gridspec_kw={'width_ratios': [1, 0.5, 0.5]},constrained_layout=True)
    #fig.subplots_adjust(hspace=1.0)
    #fig.tight_layout()
    plt.setp(ax[0, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[3, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[4, 0].get_xticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[0, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[3, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[4, 0].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[0, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[1, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[2, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[3, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[4, 1].get_xticklabels(), rotation='vertical', fontsize=6)
    plt.setp(ax[0, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[1, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[2, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[3, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.setp(ax[4, 1].get_yticklabels(), rotation='horizontal', fontsize=6)
    plt.suptitle('J and H Matrices Cluster ' + str(clustid))
    cmap = 'jet'
    # A J
    ax[0, 0].title.set_text('Jmat Afam Cluster: ' + str(clustid))
    ax[0, 0].title.set_size(fontsize=6)
    pos = ax[0, 0].imshow(ajdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    abar = fig.colorbar(pos, ax=ax[0, 0], fraction=0.046, pad=0.04)
    ax[0, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[0, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[0, 0].set_xticklabels(np.arange(2, n, 1))
    ax[0, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[0, 0].grid(True, color='g', lw=0.1)
    ax[0, 0].set_ylabel('i')
    # A H
    ax[0, 1].imshow(ah, cmap=cmap, aspect='equal')
    ax[0, 1].title.set_text('H mat Afam Cluster: ' + str(clustid))
    ax[0, 1].title.set_size(fontsize=6)
    ax[0, 1].set_xticks(np.arange(0, 21, 1))
    ax[0, 1].set_yticks(np.arange(0, n, 1))
    ax[0, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[0, 1].set_yticklabels(np.arange(1, n+1, 1))
    ax[0, 1].set_xlabel('Amino Acid')
    ax[0, 1].set_ylabel('i')
    # ax[0, 1].grid(True, color='b', lw=0.1)
    # C J
    ax[1, 0].title.set_text('Jmat Cfam Cluster: ' + str(clustid))
    ax[1, 0].title.set_size(fontsize=6)
    pos = ax[1, 0].imshow(cjdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    cbar = fig.colorbar(pos, ax=ax[1, 0], fraction=0.046, pad=0.04)
    ax[1, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[1, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[1, 0].set_xticklabels(np.arange(2, n, 1))
    ax[1, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[1, 0].grid(True, color='g', lw=0.1)
    ax[1, 0].set_ylabel('i')
    # C H
    ax[1, 1].imshow(ch, cmap=cmap, aspect='equal')
    ax[1, 1].title.set_text('H mat Cfam Cluster: ' + str(clustid))
    ax[1, 1].title.set_size(fontsize=6)
    ax[1, 1].set_xticks(np.arange(0, 21, 1))
    ax[1, 1].set_yticks(np.arange(0, n, 1))
    ax[1, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[1, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[1, 1].set_xlabel('Amino Acid')
    ax[1, 1].set_ylabel('i')
    # G J
    ax[2, 0].title.set_text('Jmat Gfam Cluster: ' + str(clustid))
    ax[2, 0].title.set_size(fontsize=6)
    pos = ax[2, 0].imshow(gjdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    gbar = fig.colorbar(pos, ax=ax[2, 0], fraction=0.046, pad=0.04)
    ax[2, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[2, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[2, 0].set_xticklabels(np.arange(2, n, 1))
    ax[2, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[2, 0].grid(True, color='g', lw=0.1)
    ax[2, 0].set_ylabel('i')
    # G H
    ax[2, 1].imshow(gh, cmap=cmap, aspect='equal')
    ax[2, 1].title.set_text('H mat Gfam Cluster: ' + str(clustid))
    ax[2, 1].title.set_size(fontsize=6)
    ax[2, 1].set_xticks(np.arange(0, 21, 1))
    ax[2, 1].set_yticks(np.arange(0, n, 1))
    ax[2, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[2, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[2, 1].set_xlabel('Amino Acid')
    ax[2, 1].set_ylabel('i')
    # T J
    ax[3, 0].title.set_text('Jmat Tfam Cluster: ' + str(clustid))
    ax[3, 0].title.set_size(fontsize=6)
    pos = ax[3, 0].imshow(tjdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    tbar = fig.colorbar(pos, ax=ax[3, 0], fraction=0.046, pad=0.04)
    ax[3, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[3, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[3, 0].set_xticklabels(np.arange(2, n, 1))
    ax[3, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[3, 0].grid(True, color='g', lw=0.1)
    ax[3, 0].set_ylabel('i')
    # T H
    ax[3, 1].imshow(th, cmap=cmap, aspect='equal')
    ax[3, 1].title.set_text('H mat Afam Cluster: ' + str(clustid))
    ax[3, 1].title.set_size(fontsize=6)
    ax[3, 1].set_xticks(np.arange(0, 21, 1))
    ax[3, 1].set_yticks(np.arange(0, n, 1))
    ax[3, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[3, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[3, 1].set_xlabel('Amino Acid')
    ax[3, 1].set_ylabel('i')
    # FULL J
    ax[4, 0].title.set_text('Jmat FULL Cluster: ' + str(clustid))
    ax[4, 0].title.set_size(fontsize=6)
    pos = ax[4, 0].imshow(fulljdisp, cmap='seismic', aspect='equal', vmin=-1, vmax=1)
    fullbar = fig.colorbar(pos, ax=ax[4, 0], fraction=0.046, pad=0.04)
    ax[4, 0].set_xticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[4, 0].set_yticks(np.arange(-.5, (n - 2) * 21, 21))
    ax[4, 0].set_xticklabels(np.arange(2, n, 1))
    ax[4, 0].set_yticklabels(np.arange(1, n - 1, 1))
    ax[4, 0].set_xlabel('j')
    ax[4, 0].grid(True, color='g', lw=0.1)
    ax[4, 0].set_ylabel('i')
    # FULL H
    ax[4, 1].imshow(fullh, cmap=cmap, aspect='equal')
    ax[4, 1].title.set_text('H mat FULL Cluster: ' + str(clustid))
    ax[4, 1].title.set_size(fontsize=6)
    ax[4, 1].set_xticks(np.arange(0, 21, 1))
    ax[4, 1].set_yticks(np.arange(0, n, 1))
    ax[4, 1].set_xticklabels(np.arange(1, 22, 1))
    ax[4, 1].set_yticklabels(np.arange(1, n + 1, 1))
    ax[4, 1].set_xlabel('Amino Acid')
    ax[4, 1].set_ylabel('i')
    abar.ax.tick_params(labelsize=2)
    cbar.ax.tick_params(labelsize=2)
    gbar.ax.tick_params(labelsize=2)
    tbar.ax.tick_params(labelsize=2)
    fullbar.ax.tick_params(labelsize=2)
    #IMAGE SEQUENCE LOGOS
    ax[0, 2].imshow(asl)
    ax[1, 2].imshow(csl)
    ax[2, 2].imshow(gsl)
    ax[3, 2].imshow(tsl)
    ax[4, 2].imshow(fsl)
    plt.savefig(clustpath + 'Clust' + str(clustid) + 'analysis.png', dpi=800)
    plt.savefig(analysispath + 'Clust' + str(clustid) + 'analysis.png', dpi=800)


clusters = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]


for clust in clusters:
    create_figure(clust)



# TCR BINDERS LOL
mdir = "C:/Users/Amber/Dropbox (ASU)/LabFolders/fernando_tcr_cluster/Data_with_cluster_id/SeqwAff/"

def writefasta(seqs, file, clustid):
    y=open(file, 'w+')
    for seq, aff in seqs:
        print('>'+'clust' + str(clustid) + '-' + str(aff), file=y)
        print(seq, file=y)
    y.close()


def Sep_G_AND_B_BINDERS(clustid):
    os.mkdir(mdir+'Clust'+str(clustid)+'/fullbinders/')
    fbpath = mdir + 'Clust' + str(clustid) + '/fullbinders/'
    gbpath = fbpath + 'gb.fasta'
    bbpath = fbpath + 'bb.fasta'
    logp = fbpath + 'clust' + str(clustid) + 'log.txt'
    log = open(logp, 'w')
    cpath = mdir + 'Clust' + str(clustid)+ '/'
    titles, seqs = dca.Fasta_Read_Aff(cpath + 'full.fasta')
    bseqs, gseqs, affs = [], [], []
    for xid, i in enumerate(titles):
        if int(i) == 1:
            bseqs.append((seqs[xid], titles[xid]))
        else:
            gseqs.append((seqs[xid], titles[xid]))
            affs.append(i)
    affmaster = set(affs)
    pctg = len(gseqs)/len(seqs)*100
    pctb = len(bseqs)/len(seqs)*100
    print('Clust ID: ' + str(clustid), file=log)
    print('PCT GB: ' + str(pctg) + ' # ' + str(len(gseqs)), file=log)
    print('PCT BB: ' + str(pctb) + ' # ' + str(len(bseqs)), file=log)
    print('Affinity Breakdown Good Binders: ', file=log)
    for i in affmaster:
        num = len([x for x in affs if x == i])
        print('A: ' + str(i) + ' # ' + str(num), file=log)
    writefasta(gseqs, gbpath, clustid)
    writefasta(bseqs, bbpath, clustid)
    print('All Files Written Cluster: ' + str(clustid))
    log.close()


clusterlist = [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37, 38,
            42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

for i in clusterlist:
    Sep_G_AND_B_BINDERS(i)


#MORE TCR BINDERS LOOOOLLLL
import dcamethods as dca
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from statistics import mean


upath = "/home/jonah/Dropbox (ASU)/"
dpath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/Ind/"
cid = 1
apath = upath + dpath + 'Analysis/'


def import_allseqpaths_cluster(cid):
    cpath = upath + dpath + 'Clust' + str(cid) + '/'
    # Seqs
    aaas = cpath + 'aaa.fasta'
    accs = cpath + 'acc.fasta'
    aggs = cpath + 'agg.fasta'
    atts = cpath + 'att.fasta'
    cacs = cpath + 'cac.fasta'
    ccgs = cpath + 'ccg.fasta'
    ctas = cpath + 'cta.fasta'
    gags = cpath + 'gag.fasta'
    gcts = cpath + 'gct.fasta'
    gtcs = cpath + 'gtc.fasta'
    tats = cpath + 'tat.fasta'
    tcas = cpath + 'tca.fasta'
    tgcs = cpath + 'tgc.fasta'
    ttgs = cpath + 'ttg.fasta'
    seqL = [aaas, accs, aggs, atts, cacs, ccgs, ctas, gags, gcts, gtcs, tats, tcas, tgcs, ttgs]
    nL = dca.getn(seqL[0])
    return seqL, nL

def import_allHJ_from_cluster(cid):
    cpath = upath + dpath + 'Clust' + str(cid) +'/'
    #H and J Matrices
    aaahp = cpath + str(cid) + 'aaa.h'
    aaajp = cpath + str(cid) + 'aaa.j'
    accjp = cpath + str(cid) + 'acc.j'
    acchp = cpath + str(cid) + 'acc.h'
    aggjp = cpath + str(cid) + 'agg.j'
    agghp = cpath + str(cid) + 'agg.h'
    attjp = cpath + str(cid) + 'att.j'
    atthp = cpath + str(cid) + 'att.h'

    cacjp = cpath + str(cid) + 'cac.j'
    cachp = cpath + str(cid) + 'cac.h'
    ccgjp = cpath + str(cid) + 'ccg.j'
    ccghp = cpath + str(cid) + 'ccg.h'
    ctajp = cpath + str(cid) + 'cta.j'
    ctahp = cpath + str(cid) + 'cta.h'

    gagjp = cpath + str(cid) + 'gag.j'
    gaghp = cpath + str(cid) + 'gag.h'
    gctjp = cpath + str(cid) + 'gct.j'
    gcthp = cpath + str(cid) + 'gct.h'
    gtcjp = cpath + str(cid) + 'gtc.j'
    gtchp = cpath + str(cid) + 'gtc.h'

    tathp = cpath + str(cid) + 'tat.h'
    tatjp = cpath + str(cid) + 'tat.j'
    tcajp = cpath + str(cid) + 'tca.j'
    tcahp = cpath + str(cid) + 'tca.h'
    tgcjp = cpath + str(cid) + 'tgc.j'
    tgchp = cpath + str(cid) + 'tgc.h'
    ttgjp = cpath + str(cid) + 'ttg.j'
    ttghp = cpath + str(cid) + 'ttg.h'

    aaah, N, q = dca.sorthmat_plmDCA_autoNandq(aaahp)
    aaaj = dca.sortjmat_plmDCA(aaajp, N, q)
    acch = dca.sorthmat_plmDCA(acchp, N, q)
    accj = dca.sortjmat_plmDCA(accjp, N, q)
    aggh = dca.sorthmat_plmDCA(agghp, N, q)
    aggj = dca.sortjmat_plmDCA(aggjp, N, q)
    atth = dca.sorthmat_plmDCA(atthp, N, q)
    attj = dca.sortjmat_plmDCA(attjp, N, q)
    cach = dca.sorthmat_plmDCA(cachp, N, q)
    cacj = dca.sortjmat_plmDCA(cacjp, N, q)
    ccgh = dca.sorthmat_plmDCA(ccghp, N, q)
    ccgj = dca.sortjmat_plmDCA(ccgjp, N, q)
    ctah = dca.sorthmat_plmDCA(ctahp, N, q)
    ctaj = dca.sortjmat_plmDCA(ctajp, N, q)
    gagh = dca.sorthmat_plmDCA(gaghp, N, q)
    gagj = dca.sortjmat_plmDCA(gagjp, N, q)
    gcth = dca.sorthmat_plmDCA(gcthp, N, q)
    gctj = dca.sortjmat_plmDCA(gctjp, N, q)
    gtch = dca.sorthmat_plmDCA(gtchp, N, q)
    gtcj = dca.sortjmat_plmDCA(gtcjp, N, q)
    tath = dca.sorthmat_plmDCA(tathp, N, q)
    tatj = dca.sortjmat_plmDCA(tatjp, N, q)
    tcah = dca.sorthmat_plmDCA(tcahp, N, q)
    tcaj = dca.sortjmat_plmDCA(tcajp, N, q)
    tgch = dca.sorthmat_plmDCA(tgchp, N, q)
    tgcj = dca.sortjmat_plmDCA(tgcjp, N, q)
    ttgh = dca.sorthmat_plmDCA(ttghp, N, q)
    ttgj = dca.sortjmat_plmDCA(ttgjp, N, q)

    JmatL = [aaaj, accj, aggj, attj, cacj, ccgj, ctaj, gagj, gctj, gtcj, tatj, tcaj, tgcj, ttgj]
    HmatL = [aaah, acch, aggh, atth, cach, ccgh, ctah, gagh, gcth, gtch, tath, tcah, tgch, ttgh]

    return JmatL, HmatL

nameL = ['aaa', 'acc', 'agg', 'att', 'cac', 'ccg', 'cta', 'gag', 'gct', 'gtc', 'tat', 'tca', 'tgc', 'ttg']
nameD = {'aaa':1, 'acc':2, 'agg':3, 'att':4, 'cac':5, 'ccg':6, 'cta':7, 'gag':8, 'gct':9, 'gtc':10, 'tat':11, 'tca':12, 'tgc':13, 'ttg':14}

JmatMaster = []

seqMaster = []
Nlist = []

clusters= [1, 3, 4, 5, 7, 8, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 29, 30, 31, 32, 34, 37,
           38, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]

for clust in clusters:
    seqL, nL = (import_allseqpaths_cluster(clust))
    Nlist.append(nL)
    seqMaster.append(seqL)
# for xid, x in enumerate(nameL):
#     J = JmatL[xid]
#     H = HmatL[xid]
#     datax = []
#     dataE = []
#     for yid, y in enumerate(seqL):
#         seqs = dca.Fasta_Read_SeqOnly(y)
#         nrgs = []
#         for seq in seqs:
#             dataE.append(dca.Calc_Energy_TCR(seq, J, H, N))
#             datax.append(nameD[nameL[yid]])
#     fig, ax = plt.subplots()
#     ax.scatter(datax, dataE)
#     ax.title.set_text(x + ' HJ Comparison')
#     ax.set_xticks(np.arange(1,15,1))
#     ax.set_xticklabels(nameL)
#     ax.set_ylabel('Energy')
#     ax.set_xlabel('Families')
#     plt.savefig(cpath+str(cid)+str(x)+'indcomp.png', dpi=600)
#     plt.close()

viridis = cm.get_cmap('viridis',len(clusters))
for cid, clust in enumerate(clusters):
    if clust == 50: ## All HJs we have currently
        break
    # Open Folder for Cluster Analysis
    jmatlist, hmatlist = import_allHJ_from_cluster(clust)
    os.mkdir(apath + 'Clust' + str(clust))
    clpath = apath + 'Clust' + str(clust)  + '/'
    filelist = [clpath + str(clust) + nameL[x] + '.png' for x in range(len(nameL))]
    for fid, jmat in enumerate(jmatlist):
        J = jmat
        H = hmatlist[fid]
        # Calculate All sequences/ All Clusters mean energies
        plotdata = np.full((len(clusters), len(nameL)), 0.0)
        stddata = np.full((len(clusters), len(nameL)), 0.0)
        mindata = np.full((len(clusters), len(nameL)), 0.0)
        maxdata = np.full((len(clusters), len(nameL)), 0.0)
        for sfid, seqL in enumerate(seqMaster): #iterates through clusters
            for sid, seqF in enumerate(seqL): #iterates through families
                plotdata[sfid, sid], stddata[sfid, sid], mindata[sfid, sid], maxdata[sfid, sid] = dca.Stats_Energy_Seqs(seqF, J, H, Nlist[cid])
        allx, ally, allz = [], [], []
        statfile = clpath + nameL[fid] + '.stats'
        dca.write_stat_file_TCR(statfile, clusters, nameL, plotdata, stddata, mindata, maxdata)
        fig = plt.figure(figsize=(5, 15))
        ax = fig.add_subplot(111)
        for i in range(14): #iterate through families
            cdata = list(plotdata[:, i])
            allz += clusters
            if i == 0:
                for aid, cluste in enumerate(clusters):
                    ax.annotate(cluste, (1, cdata[aid]))
            elif i == 13:
                for aid, cluste in enumerate(clusters):
                    ax.annotate(cluste, (14, cdata[aid]))
            pdata = [i+1 for x in cdata]
            allx += pdata
            ally += cdata
        ax.scatter(allx, ally, c=allz, s=10)
        for ccid, cclust in enumerate(clusters):
            x = np.arange(1, 15, 1)
            y = plotdata[ccid, :]
            col = (cclust - min(clusters))/(max(clusters) - min(clusters))
            ax.plot(x,y,c=viridis(col))
        ax.title.set_text(nameL[fid] + ' HJ Comparison')
        ax.set_xticks(np.arange(1,15,1))
        ax.set_xticklabels(nameL)
        ax.set_ylabel('Energy')
        ax.set_xlabel('Families')
        plt.savefig(filelist[fid], dpi=600)
        plt.close()


# WEBLOGOS
ubuntpath = "/home/jonah/Dropbox (ASU)/"
# droppath = "LabFolders/fernando_tcr_cluster/Data_with_cluster_id/FullSeq/"
droppath = "Projects/DCA/GenSeqs/"
fullpath = ubuntpath + droppath
blpath = "/home/jonah/bl-dca/"
plmpath = ubuntpath + droppath
'''
for clust in clusters:
    clustid = str(int(clust))
    clustpath = fullpath + 'Clust' + clustid

    # input files
    afasta = clustpath + '/afam.fasta'
    cfasta = clustpath + '/cfam.fasta'
    gfasta = clustpath + '/gfam.fasta'
    tfasta = clustpath + '/tfam.fasta'
    fullfasta = clustpath + '/full.fasta'

    # output files
    aout = clustpath + '/asl.png'
    cout = clustpath + '/csl.png'
    gout = clustpath + '/gsl.png'
    tout = clustpath + '/tsl.png'
    fullout = clustpath + '/fsl.png'



    sp.check_call(['weblogo --errorbars NO --resolution 600 -F png < ' + cfasta + ' > ' + cout], stdout=sp.DEVNULL, shell=True)
    sp.check_call(['weblogo --errorbars NO --resolution 600 -F png < ' + gfasta + ' > ' + gout], stdout=sp.DEVNULL, shell=True)
    sp.check_call(['weblogo --errorbars NO --resolution 600 -F png < ' + tfasta + ' > ' + tout], stdout=sp.DEVNULL, shell=True)
    sp.check_call(['weblogo --errorbars NO --resolution 600 -F png < ' + fullfasta + ' > ' + fullout], stdout=sp.DEVNULL, shell=True)
'''
'''
# input files
in5 = fullpath + '/5thgs.txt'
in6 = fullpath + '/6thgs.txt'
in7 = fullpath + '/7thgs.txt'
in8 = fullpath + '/8thgs.txt'
out5 = fullpath + '/5sl.png'
out6 = fullpath + '/6sl.png'
out7 = fullpath + '/7sl.png'
out8 = fullpath + '/8sl.png'
'''

# SEL SEQS
sseqpath = '/home/jonah/Desktop/Current/v4/8selseqs.txt'

cutoff = 89
actualhighestscoring = 'AGGGATGATGTGTGGTAGGCCTATGATGGGTAGGGTGGTG'
ssels = []
seqs = dca.Fasta_Read_SeqOnly(sseqpath)
print(seqs)
for xid, x in enumerate(seqs):
    if dca.Calc_Energy(seqs[xid], J, H) > cutoff:
        ssels.append(seqs[xid])

# Get rid of any duplicate Seqs
fseqs = dca.prune_alignment(ssels, 1.0)
afseqs = []
nrgs, simm = [], []
for x in fseqs:
    tmpa, simscore, tmpb = dca.ensemble_checker(allseqpath, x)[0]
    if simscore != 1.0:
        afseqs.append(x)
        simm.append(simscore)

for x in afseqs:
    nrgs.append(dca.Calc_Energy(x, J, H))
    # simm.append(dca.sim_score(actualhighestscoring, x))

out = '/home/jonah/Desktop/Current/v4/8gbposs.txt'
o = open(out, 'w')
for xid, x in enumerate(afseqs):
    print(x, simm[xid], nrgs[xid], file = o)
o.close()



import dcamethods as dca

fam = 8
cpath = "/home/jonah/Desktop/Current/"
allseqpath = cpath + str(fam) + 'all.txt'

infile = '/home/jonah/Desktop/Current/aptamerfinal/figs/t3k.txt'
out = '/home/jonah/Desktop/Current/aptamerfinal/figs/hb1poss3k.txt'


hb = 'AGGGTAGGTGTGGATGATGCCTAGGATGGGTAGGGTGGTG'
i = open(infile, 'r')
o = open(out, 'w')
count=1
for line in i:
    data = line.split()
    seq, energy, mut_num = data[0], data[1], data[2]
    if data[3] == '*'+hb:
        count+=1
        de, db, dg = dca.avgdis_ensemble(seq, allseqpath), dca.avgdis_bb(seq, allseqpath), dca.avgdis_gb(seq, allseqpath)
        mut_loc = dca.diffs_seqs(hb, seq)
        print('Seq', count, ':', seq, file=o)
        print('        1   5    10   15   20   25   30   35    ', file=o)
        print('Energy', energy, 'Mutation#', mut_num, 'Mut_Loc:',mut_loc, file=o)
        print("Avg_Ensemble:", de, 'Avg_dist_bb', db, "Avg_dist_gb", dg, file=o)
        print('\n', file=o)
i.close()
o.close()


