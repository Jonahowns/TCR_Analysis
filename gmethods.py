import dcamethods as dca
from sklearn.model_selection import train_test_split


upath = "/home/jonah/Dropbox (ASU)/"
datap = "Projects/DCA/GuntherAptamers/Selex_Data/"

r15p = upath + datap + 'PAL_Anna_R15_counts.txt'

def read_gfile_alldata(filep):
    o = open(filep, 'r')
    affs, seqs = [], []
    for line in o:
        d = line.split()
        aff, seq = int(d[0]), str(d[1].rsplit())
        affs.append(aff)
        seqs.append(seq)
    o.close()
    return affs, seqs

def data_prop(affs, seqs):
    a2 = set(affs)
    a2s = sorted(a2)
    for a in a2s:
        c = 0
        for aas in affs:
            if aas == a:
                c += 1
        print('Copy Number:', a, 'Count:',  c)
    c = 0
    ltot = []
    for aid, a in enumerate(affs):
        ls = []
        if a > 100:
            c += 1
            l = len(seqs[aid])
            ls.append(l)
            ltot += ls
    print('Higher than 100:', c)
    lst = sorted(set(ltot))
    for a in lst:
        c = 0
        for sid, s in enumerate(seqs):
            if affs[sid] > 100:
                if len(s) == a:
                    c += 1
        print('Affinty > 100', 'Length:', a, 'Number of Sequences', c)
    ltotal = []
    for s in seqs:
        l = len(s)
        ltotal.append(l)
    lp = set(ltotal)
    lps = sorted(lp)
    for x in lps:
        c = 0
        for aas in seqs:
            if len(aas) == x:
                c+=1
        print('Length:', x, 'Number of Sequences', c)


a_s, s_s = read_gfile_alldata(r15p)
data_prop(a_s, s_s)
X_train, X_test, Y_train, Y_test = train_test_split(s_s, a_s, test_size=0.2, random_state=0)


