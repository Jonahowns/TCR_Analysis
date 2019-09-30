#!/usr/bin/env python


import sys
import random
import dcamethods as dca
import math
import numpy as np
import random



# CSV to fasta converter as well as data preprocessing for RNA aptamer selection

def prune_alignment(names, seqs, simt=0.99, prevseq=[], prevnames=[]):
    final_choice_names = []
    final_choice_seqs = []
    if prevseq:
        final_choice_seqs += prevseq
        final_choice_names += prevnames
    for sid, seq in enumerate(seqs):
        if sid < len(prevseq):
            continue
        if sid % 1000 == 0:
            print(sid)
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
            final_choice_names.append(names[sid])
            final_choice_seqs.append(seq.upper())
    print('INFO: reduced length of alignment from %d to %d due to sequence similarity' % (
    len(seqs), len(final_choice_seqs)), file=sys.stderr),
    return final_choice_names, final_choice_seqs

def import_allseqs_lowaff(csvfile, sim):
    o = open(csvfile)
    allseqs = []
    allaffs = []
    c = 1
    prevseqs = []
    prevaffs = []
    while True:
        if c % 4 == 0:
            print('Large Pruning-')
            prunedseqs, prunedaffs = prune_alignment(allseqs, allaffs, sim, prevseqs, prevaffs)
            prevseqs = allseqs
            prevaffs = allaffs
            allseqs = prunedseqs
            allaffs = prunedaffs
        print('Processing Chunk  ' + str(c))
        c += 1
        seqs, affs = [], []
        chunk = []
        try:
            for x in range(10000):
                chunk.append(next(o))
        except StopIteration:
            break
        if not chunk:
            break
        for line in chunk:
            seq, aff = line.split(',')
            aff.rstrip()
            if int(aff) > 1:
                continue
            if len(list(seq)) != 40:
                continue
            seqs.append(seq)
            affs.append(aff)
        selseq, selaffs = prune_alignment(affs, seqs, sim)
        allseqs += selseq
        allaffs += selaffs
        print('Current Seq Number :' + str(len(allseqs)))
    return allseqs, allaffs

def import_allseqs_highaff(csvfile):
    o = open(csvfile, 'r')
    chunk = []
    for x in range(10000):
        chunk.append(next(o))
    seqs, affs = [], []
    for line in chunk:
        seq, aff = line.rstrip().split(';')
        if int(aff) == 1:
            break
        if len(list(seq)) != 40:
            continue
        seqs.append(seq)
        affs.append(aff)
    o.close()
    return seqs, affs

def import_allseqs_highaff_twogroups(csvfile):
    o = open(csvfile, 'r')
    chunk = []
    for x in range(10000):
        chunk.append(next(o))
    vseqs, vaffs, gseqs, gaffs = [], [], [], []
    for line in chunk:
        seq, aff = line.rstrip().split(',')
        if int(aff) == 1:
            break
        if len(list(seq)) != 40:
            continue
        if int(aff) > 3:
            vseqs.append(seq)
            vaffs.append(aff)
        else:
            gseqs.append(seq)
            gaffs.append(aff)
    o.close()
    return vseqs, vaffs, gseqs, gaffs


def final_pruning_lowaff(fastafile, sim):
    o = open(fastafile)
    seqs = []
    for line in o:
        if line.startswith('>'):
            continue
        else:
            seqs.append(line.rstrip())
    o.close()
    affs = [1 for x in seqs]
    selaffs, selseqs = prune_alignment(affs, seqs, sim)
    return selseqs, selaffs

def create_test_seqs(csvfile):
    o = open(csvfile, 'r')
    seqs, affs = [], []
    for line in o:
        seq, aff = line.rstrip().split(',')
        if len(list(seq)) != 40:
            continue
        seqs.append(seq)
        affs.append(aff)
    print(len(seqs))
    o.close()
    sseq, saff = [], []
    check = False
    while check is False:
        for x in range(20000):
            r = random.randint(0, len(seqs)-1)
            print(r)
            sseq.append(seqs[r])
            saff.append(affs[r])
        paffs, pseqs = prune_alignment(saff, sseq)
        sseq = pseqs
        saff = paffs
        if len(sseq) < 100000:
            continue
        else:
            check = True
    return sseq, saff

def all_seqs(csvfile):
    o = open(csvfile, 'r')
    seqs, affs = [], []
    for line in o:
        seq, aff = line.rstrip().split(';')
        if len(list(seq)) != 40:
            continue
        seqs.append(seq)
        affs.append(aff)
    return seqs, affs


def ensemble_checker(seqfile, *seqs):
    eaffs, eseqs = dca.Fasta_Read_Aff(seqfile)
    results = []
    for seq in seqs:
        seqoi = list(seq)
        highest_sim_score = 0.
        msseq = seq
        for eseq in eseqs:
            es = list(eseq)
            seq_similarity = 0.
            for i in range(len(es)):
                if seqoi[i] == es[i]:
                    seq_similarity += 1.
            seq_similarity /= len(seq)
            if seq_similarity > highest_sim_score:
                highest_sim_score = seq_similarity
                msseq = eseq
        results.append((seq, highest_sim_score, msseq))
    return results

def seq_splitter(seqs):
    l1, l2 = [], []
    pos = 20
    for x in seqs:
        seql = list(x)
        seq1 = seql[:pos]
        seq2 = seql[pos:]
        l1.append(''.join(seq1))
        l2.append(''.join(seq2))
    return l1, l2






if __name__ == '__main__':
    # Variables
    csvfile8 = '/home/jonah/Downloads/2HX_8th_new.csv'
    csvfile7 = '/home/jonah/Downloads/2HX_7th_new.csv'
    csvfile6 = '/home/jonah/Downloads/2HX_6th_new.csv'
    csvfile5 = '/home/jonah/Downloads/2HX_5th_new.csv'
    fastafile7 = '/home/jonah/7allbadbinders.txt'
    fastafile8 = '/home/jonah/8allbadbinders.txt'
    outfile8 = '/home/jonah/8test.txt'
    outfile7 = '/home/jonah/7test.txt'
    outall7 = '/home/jonah/7all.txt'
    outall8 = '/home/jonah/8all.txt'
    all6 = '/home/jonah/6all.txt'
    all5 = '/home/jonah/5all.txt'
    all7 = '/home/jonah/7all.txt'
    all8 = '/home/jonah/8all.txt'
    outgood8 = '/home/jonah/Desktop/8gv2.txt'
    outg1 = '/home/jonah/Desktop/split/8gsplit1.txt'
    outg2 = '/home/jonah/Desktop/split/8gsplit2.txt'

    droppath = "Projects/DCA/v2/"
    upath = "/home/jonah/Dropbox (ASU)/"
    g2path = upath + droppath + 'FamHJ/'
    g3path = upath + droppath + '3GHJ/'
    out10 = '/home/jonah/Desktop/Current/v4/10percent/s10.txt'
    out90 = '/home/jonah/Desktop/Current/v4/10percent/s90.txt'


    seqs, affs = import_allseqs_highaff(csvfile8)
    seq10, aff10, pc = [], [], []
    for i in range(math.floor(len(seqs)/10)):
        r = random.randint(0, len(seqs)-1)
        if r not in pc:
            pc.append(r)
            seq10.append(seqs[r])
            aff10.append(affs[r])
        else:
            continue

    seq90 = [x for x in seqs if x not in seq10]
    aff90 = []
    for seq in seq90:
        ai = seqs.index(seq)
        aff90.append(affs[ai])

    # o = open(out10, 'w')
    # for iid, i in enumerate(aff10):
    #     print('>' + 'seq' + str(iid) + '-' + str(i), file=o)
    #     print(seq10[iid], file=o)
    # o.close()

    x = open(out90, 'w')
    for iid, i in enumerate(aff90):
        print('>' + 'seq' + str(iid) + '-' + str(i), file=x)
        print(seq90[iid], file=x)
    x.close()






