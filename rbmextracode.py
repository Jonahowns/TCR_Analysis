import dcamethods as dca
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
# from sklearn.model_selection import train_test_split


upath = "/home/jonah/Dropbox (ASU)/"
wpath = "C:/Users/Amber/Dropbox (ASU)/"
datap = "Projects/DCA/GuntherAptamers/Selex_Data/"
datat = "Projects/DCA/ThrombinAptamers/"
datao = "Projects/DCA/ThrombinAptamers/v4/split/"
datarbm = "Projects/DCA/rbm_rna_v1/"
analysisrbm = upath + "LabFolders/Jonah_projects/RBM/"

def likelihood_plot_rmb(affs, like, title, out):
    plt.scatter(affs, like, c='r', s=0.2)
    plt.title('Affinity vs. Likelihood :' + title)
    plt.xlabel('Affinity')
    plt.ylabel('Likelihood')
    plt.savefig(out)
    plt.close()


def likelihood_plot_rmb_wRscore(affs, likeli, title, outpath):
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
    cutoff = max([y for x, y in api if x == 1])
    #plt.plot(xl, [cutoff for i in xl], ':b')
    plt.scatter(affs, likeli, color='r', s=0.5)
    plt.title(title)
    plt.ylabel('Likelihood')
    plt.xlabel('Affinity, Calc R Score: ' + str(linreg[2]))
    plt.savefig(outpath, dpi=600)
    plt.close()


def read_likeli(filep):
    o = open(filep)
    ls = []
    for line in o:
        data = line.split()
        interest = float(data[1])
        ls.append(interest)
    o.close()
    return ls

afs, seqs = dca.Fasta_Read_Aff(upath + datarbm + '8all.txt')

rawout = 'rbm_10hidden_raw_weights.dat'
normout = 'rbm_10hidden_norm_weights.dat'
invout = 'rbm_10hidden_inv_weights.dat'
invnor = 'rbm_10hidden_norm_inv_weights.dat'
rbmins = [rawout, normout, invout, invnor]
rbmins2 = [upath + datarbm + str(x.split('.')[0]) + '_likelihoods.txt' for x in rbmins]
rbmout = ['Lplot_' + str(x.split('.')[0]) + '_wR.png' for x in rbmins]
# lbmtitles = [str(x) + ' Hidden Nodes' for x in np.arange(10, 80, 10)]
# lbmouts = ['Lplot_' + str(x) + 'hiddenwR.png' for x in np.arange(10, 80, 10)]
# lbmins = [wpath + datarbm + str(x) + 'hidden_likelihoods.txt' for x in np.arange(10, 80, 10)]
for i in range(len(rbmins2)):
    l = read_likeli(rbmins2[i])
    likelihood_plot_rmb_wRscore(afs, l, rbmins2[i], analysisrbm+rbmout[i])