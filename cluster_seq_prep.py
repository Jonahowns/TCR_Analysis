#!/usr/bin/env python
#Data preprocessing for TCR DATA

import numpy as np
import statistics as stats
import sys


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



'''
    tmpfile = '/home/jonah/fastatmp.txt'


    sim = 0.75

    affs = []
    seqs = []

    # Select for only high affinity sequences
    seql = []
    lencatcher = 0





    if count > 25000:
        print("That is probably too many sequences for mafft... exiting now")
        sys.exit()
        
    print(str(count) + ' sequences are being written to '+tmpfile+'\n')

    #Write Fastafile for mafft
    o=open(tmpfile,'w+')
    for sid in range(len(finalnames)):
        print('>'+str(finalnames[sid]),file=o)
        print(finalrecords[sid],file=o)
    o.close()

    #Align Sequences
    #sp.check_call(['mafft --auto --leavegappyregion ' + tmpfile + ' > /home/jonah/Downloads/aligned.txt'],stdout=sp.DEVNULL,shell=True)
    sp.check_call(['clustalw '+tmpfile > ' > /home/jonah/Downloads/aligned.txt'],stdout=sp.DEVNULL,shell=True)


    #Read Alignment File and prune it
    alignedseqs=[]
    alignednames=[]
    full=open('/home/jonah/Downloads/aligned.txt','r').read().splitlines()
    for line in full:
        if line.startswith('>'):
            alignednames.append(line)
            full[full.index(line)]='>'
    seqflat=''.join(full)
    alignedseqs=seqflat.split('>')
    del(alignedseqs[0])  
    print(len(alignedseqs))
    

    # Prune Similar aligned Sequences
    newnames, newseqs = prune_alignment(finalnames, finalrecords, sim)

    print('Percent of Sequences left after pruning --' + str(len(newseqs) / len(finalrecords) * 100) + '%')

    # Final Sequences Aligned
    f = open(outfile, 'w')
    for x in range(len(newnames)):
        print('>' + newnames[x], file=f)
        print(newseqs[x], file=f)
    f.close()
    '''