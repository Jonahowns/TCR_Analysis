import dcamethods as dca
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy import stats


upath = "/home/jonah/Dropbox (ASU)/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
datap = "Projects/DCA/GunterAptamers/Selex_Data/"
datat = "Projects/DCA/ThrombinAptamers/"
datao = "Projects/DCA/ThrombinAptamers/v4/split/"
datarbm = "Projects/DCA/rbm_rna_v1/"
analysisrbm = upath + "LabFolders/Jonah_projects/RBM/"

def read_likeli_new(filep):
    o = open(filep)
    ls, afs = [], []
    for line in o:
        data = line.split()
        ls.append(float(data[1]))
        afs.append(int(data[0]))
    o.close()
    return afs, ls

def likelihood_plot_rmb_wRscore(affs, likeli, title, outpath, cutoff='no'):
    a_s = list(set(affs))
    api = list(zip(affs, likeli))
    highestaff = 1
    datax, datae = [], []
    for x in a_s:
        if x > highestaff: highestaff = x
        prospects = [l for (aff, l) in api if aff == x]
        datax.append(x)
        datae.append(max(prospects))
    linreg = stats.linregress(datax, datae)
    xl = np.linspace(0, highestaff, 100)
    plt.plot(xl, xl * linreg[0] + linreg[1], ':r')
    if cutoff == 'yes':
        cutoff = max([y for x, y in api if x == 1])
        plt.plot(xl, [cutoff for i in xl], ':b')
    plt.scatter(affs, likeli, color='r', s=0.5)
    plt.title(title)
    plt.ylabel('Likelihood')
    plt.xlabel('Affinity, Calc R Score: ' + str(linreg[2]))
    plt.savefig(outpath, dpi=600)
    plt.close()



plmA = 'plmAnalysis/'

computer_path = upath

plmp = computer_path + datap + plmA
rbmp = computer_path + datarbm
trainp = computer_path + datap + 'r15_train.txt'
testp = computer_path + datap + 'r15_test.txt'

rbmin_tmp = ['r15_g_' + str(x) + 'hidden_likelihood' for x in np.arange(10, 60, 10)]
rbmin_sqtmp = ['r15_g_' + str(x) + 'hidden_Wsq_likelihood' for x in np.arange(10, 60, 10)]
nosq_all = [(x+'_train.txt', x+'_test.txt') for x in rbmin_tmp]
sq_all = [(x+'_train.txt', x+'_test.txt') for x in rbmin_sqtmp]

all_likelis = sq_all + nosq_all
all_titles = rbmin_sqtmp + rbmin_tmp
'''
cid= 0
for train, test in all_likelis:
    tra, trl = read_likeli_new(rbmp+train)
    tea, tel = read_likeli_new(rbmp+test)
    likelihood_plot_rmb_wRscore(tra, trl, all_titles[cid], analysisrbm + train.split('.')[0] + '_plot.png', cutoff='no')
    cid += 1
    likelihood_plot_rmb_wRscore(tea, tel, all_titles[cid], analysisrbm + test.split('.')[0] + '_plot.png', cutoff='no')
    cid += 1
'''

r15hp = plmp + 'r15_train.h'
r15jp = plmp + 'r15_train.j'

cutoffs = ['c100', 'c300', 'c500', 'c1k']
jps = [plmp + x + '_r15_train_s.j' for x in cutoffs]
hps = [plmp + x + '_r15_train_s.h' for x in cutoffs]
trainps = [plmp + 'r15_train_' + x +'.txt' for x in cutoffs]
testps = [plmp + 'r15_test_' + x +'.txt' for x in cutoffs]
trainws = [plmp + 'r15_weights2_' + x +'.txt' for x in cutoffs]
outs = ['c100_s.png', 'c300_nw2.png', 'c500_nw2.png', 'c1k_nw2.png']

'''
affs, seqs = dca.Fasta_Read_Aff(trainps[0])
weights = [(x/1000.)**2. for x in affs]
alldata = list(zip(affs, weights, seqs))
alldata.sort(key=lambda tup: tup[0])
alldata.reverse()
data_sep = list(zip(*alldata))
saffs, sweights, sseqs = list(data_sep[0]), list(data_sep[1]), list(data_sep[2])
wfile = plmp + 'r15c100_sorted_weights.txt'
sfile = plmp + 'r15c100_sorted_train.txt'
dca.write_fasta_aff(sseqs, saffs, sfile)
w = open(wfile, 'w')
for w8 in sweights:
    print(round(w8, 2), file=w)
w.close()
'''

N, q = 40, 5
for i in range(1):
    r15j = dca.sortjmat_plmDCA(jps[i], N, q)
    r15h = dca.sorthmat_plmDCA(hps[i], N, q)
    jdisp = dca.FullJ_disp(r15j, N, q)
    fig, ax = plt.subplots(1, 2)
    dca.Fig_FullJ(ax[0], 'c100', jdisp, N, q)
    dca.Fig_FullH(ax[1], 'c100', r15h, N, q)
    # dca.Raw_wRscore_subplot(ax[0], r15j, r15h, trainps[i])
    # dca.Raw_wRscore_subplot(ax[1], r15j, r15h, testps[i])
    plt.savefig(plmp + 'c100_plmparameters.png', dpi=600)
    plt.close()

'''
for i in range(4):
    affs, seqs = dca.Fasta_Read_Aff(trainps[i])
    weights = [(x/1000.)**2. for x in affs]
    o = open(trainws[i], 'w')
    for x in weights:
        print(round(x, 2), file=o)
    o.close()
'''

# r15jE, gvals, gdist = dca.TopJNorms_Jmatrix(r15j, N, q, 150)
# fig, ax = plt.subplots(1, 2)
# dca.Raw_wRscore_subplot(ax[0], r15j, r15h, trainp)
# dca.Raw_wRscore_subplot(ax[1], r15j, r15h, testp)
# plt.savefig(plmp + 'r15_train_vs_test.png', dpi=600)


