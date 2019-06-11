import numpy as np
import pandas as pd
import sys

droppath = "Projects/DCA/GenSeqs/"
#fullpath = macpath+droppath
upath = "/home/jonah/Dropbox (ASU)/"
fullpath = upath + droppath
fivegb = fullpath + '5thgoodbinders'
fivebb = fullpath + '5thbadbinders'
sixgb = fullpath + '6thgoodbinders'
sixbb = fullpath + '6thbadbinders'
sevgb = fullpath + '7thgoodbinders'
sevbb = fullpath + '7thbadbinders'
eitgb = fullpath + '8thgoodbinders'
eitbb = fullpath + '8thbadbinders'


def prune_alignment(names, seqs, simt=0.99):
    final_choice_names = []
    final_choice_seqs = []
    for sid, seq in enumerate(seqs):
        # print(sid)
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

    print('INFO: reduced length of alignment from %d to %d due to sequence similarity' % (len(seqs), len(final_choice_seqs)),
          file=sys.stderr)
    return final_choice_names, final_choice_seqs

def val_sort(inputfile, gorb):
    o = open(inputfile, 'r')
    titles =[]
    seqs = []
    for line in o:
        if line.startswith('>'):
            if gorb == 'g':
                titles.append(float(line.rstrip().split('-')[1]))
            elif gorb == 'b':
                titles.append(float(line.rstrip().split('--')[1]))
        else:
            seqs.append(line.rstrip())
    o.close()
    seqs.reverse()
    titles.reverse()
    ftitles, fseqs = prune_alignment(titles, seqs, 0.85) # 0.85 is a decent number for this!
    full = list(zip(ftitles, fseqs))
    full.sort(key=lambda tup: tup[0])
    for xid, x in enumerate(full):
        num, seq = x
        ls = list(seq)
        for xidx, i in enumerate(ls):
            if i == 'U':
                ls[xidx] = 'T'
        full[xid] = (num, ''.join(ls))
    fnp = np.array(full)
    return fnp[-50:-1]

fivebad = val_sort(fivebb, 'b')[-4:-1]
fivegood = val_sort(fivegb, 'g')[-4:-1]
sixbad = val_sort(sixbb, 'b')[-4:-1]
sixgood = val_sort(sixgb, 'g')[-4:-1]
sevbad = val_sort(sevbb, 'b')[-4:-1]
sevgood = val_sort(sevgb, 'g')[-4:-1]
eitbad = val_sort(eitbb, 'b')[-4:-1]
eitgood = val_sort(eitgb, 'g')[-4:-1]
header5 = ['5g','5g','5b','5b']
header6 = ['6g','6g','6b','6b']
header7 = ['7g','7g','7b','7b']
header8 = ['8g','8g','8b','8b']
table=np.full((16,4), 'x',dtype=str)
table[0]=header5
table[1:4,0:2]= fivegood
table[1:4,2:5]= fivebad
table[4]=header6
table[5:8,0:2]= sixgood
table[5:8,2:5]= sixbad
table[8]=header7
table[9:12,0:2]= sevgood
table[9:12,2:5]= sevbad
table[12]=header8
table[13:17,0:2]= eitgood
table[13:17,2:5]= eitbad
print(table)