import numpy as np
import math
import copy
import matplotlib.pyplot

################################################
## Universal Methods for Analysis of DCA Data ##
################################################
aa = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
rna = ['-', 'A', 'C', 'G', 'U']
rnad = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'U':4}
rnan = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'U'}
dna = ['-', 'A', 'C', 'G', 'T']
nucs = ['A', 'C', 'G', 'U']
nuc_to_id = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}

# Takes plmDCA J Matrix File and inputs the values into a N-1, N-1, q, q matrix
def sortjmat_plmDCA(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N-1, N-1, q, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0])-1, int(data[1])-2, int(data[2])-1, int(data[3])-1] = float(data[4].rstrip())
    o.close()
    return fullmatrix

# Takes plmDCA H Matrix File and inputs the values into a N-1, q matrix
def sorthmat_plmDCA(file, N, q):
    o = open(file, 'r')
    fullmatrix = np.full((N, q), 0.0)
    for line in o:
        data = line.split(',')
        fullmatrix[int(data[0]) - 1, int(data[1]) - 1] = float(data[2].rstrip())
    o.close()
    return fullmatrix

# Returns H Matrix with only the top values specified by a percentile
# by default returns values greater than the 80th percentile
# Keyword argument percentile sets the percentile
def TopH_values_disp(H, N, q, **kwargs):
    percentile = 80
    for key, value in kwargs.items():
        if key == 'percentile':
            percentile = value
        else:
            print('No keyword argument ' + key + ' found')
    Hdisp = np.full((N, q), 0.0)
    val = np.percentile(H, percentile)
    for i in range(0, N):
        for j in range(0, q):
            if H[i, j] > val:
                Hdisp[i, j] = H[i, j]
    return Hdisp

# Returns J Matrix with only the top NORMS of individual Jij Matrices specified by a percentile
# by default returns values greater than the 80th percentile
# Keyword argument percentile sets the percentile
# Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# norm values returned as list of vals and cutoff returned as single value
def TopJ_Norms_disp(J, N, **kwargs):
    percentile = 80
    dist = False
    for key, value in kwargs.items():
        if key == 'percentile':
            percentile = value
        elif key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    jnorm = np.full((N-1, N-1), 0.0)
    jdisp = np.full((N-1, N-1), 0.0)
    vals = []
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
    tval = np.percentile(jnorm, percentile)
    for i in range(N-1):
        for j in range(N-1):
            if jnorm[i, j] >= tval:
                jdisp[i, j] = jnorm[i, j]
            if jnorm[i, j] != 0.0:
                vals.append(jnorm[i, j])
    if dist:
        return jdisp, tval, vals
    else:
        return jdisp


# Return Full J Matrix in 2D so it can be used by matplotlib and shown easily
def FullJ_disp(J, N, q):
    Jdisp = np.full(((N-1)*q, (N-1)*q), 0.0)
    for i in range(N-1):
        for j in range(N-1):
            for k in range(q):
                for l in range(q):
                    if J[i, j, k, l] != 0.0:
                        Jdisp[i*q+k, j*q+l] = J[i, j, k, l]
                    else:
                        Jdisp[i*q+k, j*q+l] = 0.0
    return Jdisp


# Returns J Matrix with H values added
def HJ_Mutant(J, H, N, q):
    mutt = copy.deepcopy(J)
    for x in range(N - 1):
        for k in range(q):  # J Indices
            mutt[x, x:N, k, :] += H[x, k]
    for y in range(N - 1):
        for l in range(q):  # y states
            mutt[0:y + 1, y, :, l] += H[y + 1, l]
    return mutt


# Takes in Good Binders J and H and Bad Binders J and H
# if htype = 'ind' compares top individual values of GB H and BB H and if they are in common removes them
# if htype = 'norm' compares top norms of GB H and BB H and if they are in common removes them
# Set pcts with hnormpct and jnormpct
def Binder_Comp_JH(J, bJ, H, bH, N, q, **kwargs):
    filledJnorms = (N-1)*(N-2)/2 + N-1
    filledHind = (N * q)
    nxj = 10  # Number of J Norms used in Comparison
    nxh = 10  # Number of H Norms used n Comparison
    htype = 'good'
    hdist = False
    jdist = False
    for key, value in kwargs.items():
        if key == 'jnormpct':
            nxj = math.ceil(value/100 * filledJnorms)
        elif key == 'htype':
            htype = value
        elif key == 'hnormpct':
            if htype == 'norm':
                nxh = math.ceil(value/100 * N)
            if htype == 'ind':
                nxh = math.ceil(value/100 * filledHind)
        elif key == 'hdist':
            hdist = value
        elif key == 'jdist':
            jdist = value
        else:
            print('No keyword argument ' + key + ' found')
    topJ, valsJ, tvalJ = TopX_JNorms(J, N, nxj, pct=20, dist=True)
    topBJ, valsBJ, tvalBJ = TopX_JNorms(bJ, N, nxj, pct=20, dist=True)
    for xj, yj, val in topJ:
        for xb, yb, valb in topBJ:
            if xj == xb and yj == yb:
                J[xj, yj, :, :] = 0.0
    if htype != 'good':
        topH, valsH, tvalH = TopX_HVals(H, N, nxh, pct=20, htype=htype, dist=True)
        topBH, valsBH, tvalBH = TopX_HVals(bH, N, nxh, pct=20, htype=htype, dist=True)
        if htype == 'norm':
            for xi, val in topH:
                for yb, valb in topBH:
                    if xi == yb:
                        H[xi, :] = 0.0
        elif htype == 'ind':
            for xi, yi, val in topH:
                for xb, yb, valb in topBH:
                    if xi == xb and yi == yb:
                        H[xi, yi] = 0.0
        if hdist:
            Hdist = [valsH, tvalH, valsBH, tvalBH]
    else:
        Hdist = [0]
    if jdist:
        Jdist = [valsJ, tvalJ, valsBJ, tvalBJ]
    if hdist and jdist:
        return H, J, Hdist, Jdist
    elif hdist and not jdist:
        return H, J, hdist
    elif not hdist and jdist:
        return H, J, jdist
    elif not hdist and not jdist:
        return H, J



# Returns the highest x amount of J Norms
# Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# Highest norms returned in descending list of tuples (i, j, normval) by magnitude of normval
# norm values returned as list of vals and cutoff returned as single value
# Keyword argument percentile determines the percentile used as a cutoff
def TopX_JNorms(J, N, x, **kwargs):
    dist = False
    pct = 80
    for key, value in kwargs.items():
        if key == 'percentile':
            pct = value
        elif key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    jnorm = np.full((N-1, N-1), 0.0)
    vals = []
    jvals =[]
    for i in range(N-1):
        for j in range(N-1):
            jnorm[i, j] = np.linalg.norm(J[i, j, :, :])
            if jnorm[i, j] != 0.0:
                vals.append((i, j, jnorm[i, j]))  # 0, 0 -> 1, 2
                jvals.append(jnorm[i, j])
    tval = np.percentile(jnorm, pct)
    vals.sort(key=lambda tup: tup[2])
    ind = int(-x)
    top10 = vals[ind:-1]
    if dist:
        return top10, jvals, tval
    else:
        return top10


# Returns the highest x amount of H Norms
# Keyword argument dist determines wheter the cutoff and values are returned for use in distribution figures
# norm values returned as list of vals and cutoff returned as single value
# Keyword argument percentile determines the percentile used as a cutoff, default is 80
# Keyword htype has options 'ind', 'norm'
# 'ind' returns the top individual values in H
# 'norm' returns the top norms out of N norms in H
def TopX_HVals(H, N, x, **kwargs):
    pct = 80
    htype = 'ind'
    dist = False
    for key, value in kwargs.items():
        if key == 'percentile':
            pct = value
        elif key == 'htype':
            htype = value
        elif key == 'dist':
            dist = True
        else:
            print('No keyword argument ' + key + ' found')
    if htype == 'norm':
        hnorm = np.full((N - 1), 0.0)
        vals = []
        hvals = []
        for i in range(N-1):
            hnorm[i] = np.linalg.norm(H[i, :])
            if hnorm[i] != 0.0:
                vals.append((i, hnorm[i]))
                hvals.append(hnorm[i])
        tval = np.percentile(hvals, pct)
        vals.sort(key=lambda tup: tup[1])
    elif htype == 'ind':
        vals = []
        hvals = []
        for i in range(N - 1):
            for j in range(1, 5):
                vals.append((i, j, abs(H[i, j])))
                hvals.append(H[i, j])
        tval = np.percentile(hvals, pct)
        vals.sort(key=lambda tup: tup[2])
    ind = int(-x)
    top10 = vals[ind:-1]
    if dist:
        return top10, hvals, tval
    else:
        return top10


# Quick Method to get N by reading just the length of the first seq.. will only work if all seqs are the same length
def getn(fastafile):
    o =open(fastafile,'r')
    o.readline()
    seq = o.readline().rstrip()
    n = len(list(seq))
    o.close()
    return n


########################################################################################################################
########################################################################################################################
# Subplot Methods


# Shows Entire J Matrix
# Check Function for Keyword Arguments, available are lw, fontsize, xlabel, ylabel, vmin, vmax, title, and cmap
def Fig_FullJ(subplot, id, J, n, q, **kwargs):
    cmap = 'seismic'
    vml = -1
    vmg = 1
    lw = 0.1
    xlabel = 'j'
    ylabel = 'i'
    fontsize = 6
    title = 'JMat ID: ' + str(id)
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        elif key == 'lw':
            lw = value
        elif key == 'xlabel':
            xlabel = value
        elif key == 'ylabel':
            ylabel = value
        elif key == 'ticksize':
            fontsize = value
        elif key == 'vmin':
            vml = value
        elif key == 'vmax':
            vmg = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=6)
    subplot.imshow(J, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
    subplot.set_xticks(np.arange(-.5, (n - 2) * q, q))
    subplot.set_yticks(np.arange(-.5, (n - 2) * q, q))
    subplot.set_xticklabels(np.arange(2, n, 1))
    subplot.set_yticklabels(np.arange(1, n - 1, 1))
    subplot.grid(True, color='g', lw=lw)
    subplot.set_ylabel(ylabel)
    supplot.set_xlabel(xlabel)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)


# Produces figure on given subplot of J norms Matrix which is a parameter
# Check Function for Keyword Arguments, available are lw, fontsize, xlabel, ylabel, vmin, vmax, title, and cmap
def Fig_Jnorm(subplot, id, J, n, **kwargs):
    cmap = 'seismic'
    vml = 0
    vmg = 4
    lw = 1.0
    xlabel = 'j'
    ylabel = 'i'
    fontsize = 6
    title = 'Jmat Norms ID: ' + str(id)
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        elif key == 'lw':
            lw = value
        elif key == 'xlabel':
            xlabel = value
        elif key == 'ylabel':
            ylabel = value
        elif key == 'ticksize':
            fontsize = value
        elif key == 'vmin':
            vml = value
        elif key == 'vmax':
            vmg = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=6)
    subplot.imshow(J, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
    subplot.set_xticks(np.arange(-.5, (n - 1), 1))
    subplot.set_yticks(np.arange(-.5, (n - 1), 1))
    subplot.set_xticklabels(np.arange(2, n+1, 1))
    subplot.set_yticklabels(np.arange(1, n, 1))
    subplot.grid(True, color='g', lw=lw)
    subplot.set_ylabel(ylabel)
    subplot.set_xlabel(xlabel)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)


# Produces Figure of Full H Matrix on a given subplot
# Check Function for Keyword Arguments, available are fontsize, xlabel, ylabel, vmin, vmax, title, and cmap
def Fig_FullH(subplot, id, H, n, q,  **kwargs):
    cmap = 'seismic'
    vml = 0
    vmg = 4
    xl = False
    xlabel = 'hello'
    ylabel = 'i'
    fontsize = 6
    title = 'Hmat ID: ' + str(id)
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        elif key == 'xlabel':
            xl = True
            xlabel = value
        elif key == 'ylabel':
            ylabel = value
        elif key == 'fontsize':
            fontsize = value
        elif key == 'vmax':
            vml = value
        elif key == 'vmin':
            vmg = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.imshow(H, cmap=cmap, aspect='equal', vmin=vml, vmax=vmg)
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=fontsize)
    plt.setp(subplot.get_xticklabels(), rotation='vertical', fontsize=fontsize)
    plt.setp(subplot.get_yticklabels(), rotation='horizontal', fontsize=fontsize)
    subplot.set_xticks(np.arange(0, q+1, 1))
    subplot.set_yticks(np.arange(0, n+1, 1))
    if q == 21 and xl is False:
        subplot.set_xticklabels(aa)
        subplot.set_xlabel('Amino Acid')
    elif q == 5 and xl is False:
        subplot.set_xticklabels(rna)
        subplot.set_xlabel('Base')
    else:
        subplot.set_xlabel(xlabel)
    subplot.set_yticklabels(np.arange(1, n+1, 1))
    subplot.set_ylabel(ylabel)


# Produces Figure of Distribution of Values on a given subplot
# A Vertical Line is produced at the Cutoff Value.. typically used in showing which values in a distribution were used
# Check Function for Keyword Arguments, available are fontsize, xlabel, xmin, xmax, title, and pcolor
def Fig_Distribution_w_Cutoff(subplot, id, Values, Cutoff, **kwargs):
    title = 'Distribution ID: ' + str(id)
    plotcolor = 'r'
    xml = 0
    xmg = 2
    xlabel = 'Value'
    fontsize = 6
    for key, value in kwargs.items():
        if key == 'xlabel':
            xlabel = value
        elif key == 'pcolor':
            plotcolor = value
        elif key == 'fontsize':
            fontsize = value
        elif key == 'xmin':
            xml = value
        elif key == 'xmax':
            xmg = value
        elif key == 'title':
            title = value
        else:
            print('No keyword argument ' + key + ' found')
    deN = gaussian_kde(Values)
    xd1 = np.linspace(xml, xmg, 100)
    subplot.plot(xd1, deN(xd1), color=plotcolor)
    subplot.plot(Values, [0.01] * len(Values), '|', color='k')
    subplot.set_xlabel(xlabel)
    subplot.grid(True)
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=fontsize)
    subplot.axvline(x=Cutoff)


# Shows Premade SeqLogo on a given Subplot
# Keyword Arguments are title and fontsize
def Fig_SeqLogo(Filepath, Subplot, id):
    title = 'SeqLogo ID: ' + str(id)
    fontsize = 6
    for key, value in kwargs.items():
        if key == 'title':
            title = value
        elif key == 'fontsize':
            fontsize = value
        else:
            print('No keyword argument ' + key + ' found')
    fsl1 = mpimg.imread(Filepath)
    Subplot.imshow(fsl1)
    Subplot.axis('off')
    Subplot.title.set_text(title)
    Subplot.title.set_size(fontsize=fontsize)


# On Specified Subplot shows Individual Jij at specified x and y **NOTE J12 is equivalent to J[0, 0, :, :]
# Check Function for Keyword Arguments, available are cbar, lw, fontsize, vmin, vmax, title, cmap, and type
# Three Types are available: 'dna', 'rna' and 'pep'
def IndJij(subplot, J, x, y, id, **kwargs):
    vml = -0.5
    vmg = 0.5
    cmap = 'seismic'
    fontsize = 4
    type = 'rna'
    cbar = False
    lw = 0.1
    title = 'Jij ID: ' + str(id) + ' ' + 'Pair: ' + str(x + 1) + ' and ' + str(y + 2)
    for key, value in kwargs.items():
        if key == 'vmax':
            vmg = value
        elif key == 'vmin':
            vml = value
        elif key == 'cmap':
            cmap = value
        elif key == 'fontsize':
            fontsize = value
        elif key == 'type':
            type = value
        elif key == 'lw':
            lw = value
        elif key == 'title':
            title = value
        elif key == 'cbar':
            cbar = value
        else:
            print('No keyword argument ' + key + ' found')
    subplot.imshow(J[x, y, :, :], cmap=cmap, vmin=vml, vmax=vmg)
    if cbar:
        plt.colorbar(pos, ax=subplot, fraction=0.046, pad=0.04)
    if type == 'rna' or type == 'dna':
        subplot.set_xticks(np.arange(-.5, 4.5, 1))
        subplot.set_yticks(np.arange(-.5, 4.5, 1))
        if type == 'rna':
            subplot.set_xticklabels(rna)
            subplot.set_yticklabels(rna)
        else:
            subplot.set_xticklabels(dna)
            subplot.set_yticklabels(dna)
    elif type == 'pep':
        subplot.set_xticks(np.arange(-.5, 20.5, 1))
        subplot.set_yticks(np.arange(-.5, 20.5, 1))
        subplot.set_xticklabels(aa)
        subplot.set_yticklabels(aa)
    subplot.tick_params(axis='both', which='major', labelsize=fontsize)
    subplot.tick_params(axis='both', which='minor', labelsize=fontsize)
    subplot.grid(True, color='r', lw=lw)
    subplot.title.set_text(title)
    subplot.title.set_size(fontsize=(fontsize+2))


########################################################################################################################
########################################################################################################################
# Full Figure Methods


# Returns a figure showing the Individual Jijs of the top 10 Norms of the provided J Matrix
# OutPath is the directory the figure is being saved to
def top10norms_figure_RNA(id, J, N, OutPath):
    # Get Indices of top 10 norms
    jx = TopX_JNorms(J, N, 10)
    fig, ax = plt.subplots(2, 5, constrained_layout=True)
    for i in range(10):
        x, y, z = jx[i]
        j = i % 5
        k = 0
        if i == 0:
            IndJij(ax[k, j], J, x, y, id, vmin=0, vmax=2, type='rna', cbar=True)
        else:
            if i > 4:
                k = 1
            IndJij(ax[k, j], J, x, y, id, vmin=0, vmax=2, type='rna')

    fig.suptitle('Highest Jij Norms ID: ' + str(id))
    plt.savefig(OutPath + str(id) + 'JNormt10.png', dpi=600)



########################################################################################################################
########################################################################################################################
# Generate Seq Methods

# Calculates 'energy' of a seq based on provided J and H
# Higher Energy indicates a better binder
def Calc_Energy(seq, J, H):
    full = list(seq)
    dist = len(full)
    Jenergy = 0
    Henergy = 0
    for x in range(1, dist):
        ibase = rnad[seq[x]]
        Henergy += H[x, ibase]
        for y in range(x+1, dist):
            jbase = rnad[seq[y]]
            Jenergy += J[x-1, y-2, ibase, jbase]
    energy = Jenergy + Henergy
    return energy