import os
import dcamethods as dca

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







