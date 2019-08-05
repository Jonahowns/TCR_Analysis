#!/usr/bin/env python


import sys
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
        seq, aff = line.rstrip().split(',')
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



if __name__ == '__main__':
    # Variables
    csvfile8 = '/home/jonah/Downloads/2HX_8th_new.csv'
    csvfile7 = '/home/jonah/Downloads/2HX_7th_new.csv'
    fastafile7 = '/home/jonah/7allbadbinders.txt'
    fastafile8 = '/home/jonah/8allbadbinders.txt'
    outfile8 = '/home/jonah/8test.txt'
    outfile7 = '/home/jonah/7test.txt'


    sim = 0.85
    seqt, afft = create_test_seqs(csvfile8)

    f = open(outfile8, 'w')
    for x in range(len(afft)):
        print('>seq' + str(x) + '-' + str(afft[x]), file=f)
        print(seqt[x], file=f)
    f.close()


    # Final Sequences Aligned
    # f = open(outfile7, 'w')
    # for x in range(len(affs7)):
    #     print('>seq' + str(x)+'-'+affs7[x], file=f)
    #     print(seqs7[x], file=f)
    # f.close()


