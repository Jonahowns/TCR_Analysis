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
N, q = 40, 5
r15j = dca.sortjmat_plmDCA(r15jp, N, q)
r15h = dca.sorthmat_plmDCA(r15hp, N, q)

# r15jE, gvals, gdist = dca.TopJNorms_Jmatrix(r15j, N, q, 150)
fig, ax = plt.subplots(1, 2)
dca.Raw_wRscore_subplot(ax[0], r15j, r15h, trainp)
dca.Raw_wRscore_subplot(ax[1], r15j, r15h, testp)
plt.savefig(plmp + 'r15_train_vs_test.png', dpi=600)


