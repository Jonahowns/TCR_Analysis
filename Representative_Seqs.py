#!/usr/bin/env python


import sys



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




if __name__ == '__main__':
    # Variables
    csvfile8 = '/home/jonah/Downloads/2HX_8th_new.csv'
    csvfile7 = '/home/jonah/Downloads/2HX_7th_new.csv'
    fastafile7 = '/home/jonah/7allbadbinders.txt'
    fastafile8 = '/home/jonah/8allbadbinders.txt'
    outfile8 = '/home/jonah/8allBbindersP.txt'
    outfile7 = '/home/jonah/7allBbindersP.txt'
    x = open(outfile7, 'w')
    y = open(outfile8, 'w')

    sim = 0.85
    #seqs7, affs7 = final_pruning_lowaff(fastafile7, sim)
    seqs8, affs8 = final_pruning_lowaff(fastafile8, sim)

    # Final Sequences Aligned
    # f = open(outfile7, 'w')
    # for x in range(len(affs7)):
    #     print('>seq' + str(x)+'-'+affs7[x], file=f)
    #     print(seqs7[x], file=f)
    # f.close()

    f = open(outfile8, 'w')
    for x in range(len(affs8)):
        print('>seq' + str(x) + '-' + affs8[x], file=f)
        print(seqs8[x], file=f)
    f.close()
