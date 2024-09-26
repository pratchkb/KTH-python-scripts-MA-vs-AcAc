import numpy as np
import os
import os.path
import glob
# load proper modules
import get_optim
import geom_param
from itertools import islice
from itertools import chain
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt


def calculate_pyr_atoms(filename, molecule, natoms, hh1, hh2, pyrc):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

    pyr = []
    pyrmeth = []

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    before = []
    timstp = []

    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if " t= " in line: # Initial condition geometry
            t = line.split()
            timstp.append(float(t[2]))

    data = get_optim.read_xyzs_as_np(filename)
    for i in range(len(data[:])):
        p = geom_param.compute_pyramidalization(data[i],*pyrc)
        pyr.append(abs(p))

        h1 = geom_param.compute_bond(data[i],hh1)
        h2 = geom_param.compute_bond(data[i],hh2)

        if molecule == "MA":
            if h1 < h2:
                order = [4,3,1,7]
            else:
                order = [2,0,3,5]

        if molecule == "AcAc":
            if h1 < h2:
                order = [3,2,4,9]
            else:
                order = [0,1,2,5]

        pm = geom_param.compute_pyramidalization(data[i],*order)
        pyrmeth.append(abs(pm))

    return pyr, pyrmeth, timstp


#############################################################################################

if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    #Ntraj = [10,10]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/PYRAMIDALIZATION-both-compare/'

    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        if molecule == "MA":

            hh1 = [1,8]
            hh2 = [0,8]

            pyrs = [3,2,4,6]
            #pyr = calculate_pyr_atoms(filename, molecule, numatoms, hh1, hh2)
           

            #MAdata = collect_geom_param(run_path + "crossing-s1-s0.xyz", numatoms, h1, h2, pyrs, pyr1, pyr2)
            
            MAdata = calculate_pyr_atoms(run_path + "crossing-s1-s0.xyz", molecule, numatoms, hh1, hh2, pyrs)

        if molecule == "AcAc":
            
            hh1 = [4,13]
            hh2 = [1,13]
            
            pyrs = [2,0,3,14]
            #pyr1 = [0,1,2,5]
            #pyr2 = [3,2,4,9]

            #AcAcdata = collect_geom_param(run_path + "crossing-s1-s0.xyz", numatoms, h1, h2, pyrs, pyr1, pyr2)

            AcAcdata = calculate_pyr_atoms(run_path + "crossing-s1-s0.xyz", molecule, numatoms, hh1, hh2, pyrs)

    os.chdir(path_script)

    MAbeforet = []
    MAbeforepyrc = []
    MAbeforepyrM = []
    MAaftert = []
    MAafterpyrc = []
    MAafterpyrM = []

    for i,j,k in zip(MAdata[0], MAdata[1], MAdata[2]):
        if k < 75:
            MAbeforet.append(k)
            MAbeforepyrc.append(i)
            MAbeforepyrM.append(j)
        else:
            MAaftert.append(k)
            MAafterpyrc.append(i)
            MAafterpyrM.append(j)

    AcAcbeforet = []
    AcAcbeforepyrc = []
    AcAcbeforepyrM = []
    AcAcaftert = []
    AcAcafterpyrc = []
    AcAcafterpyrM = []

    for i,j,k in zip(AcAcdata[0], AcAcdata[1], AcAcdata[2]):
        if k < 75:
            AcAcbeforet.append(k)
            AcAcbeforepyrc.append(i)
            AcAcbeforepyrM.append(j)
        else:
            AcAcaftert.append(k)
            AcAcafterpyrc.append(i)
            AcAcafterpyrM.append(j)



    fig, ax = plt.subplots(2,1,figsize=(4,6), sharex=True)

    ax[0].scatter(MAdata[0], MAdata[1], s=2.0, c='green')
    ax[0].scatter(45.0,23.4,color='gold', edgecolor='k', marker='s', label=r'S$_1$/S$_0$-CIpy')
    ax[0].scatter(0.8,14.5,color='gold', edgecolor='k', marker='D', label=r'S$_1$/S$_0$-Tw-MECI')
    ax[0].set_ylabel(r'PyrM ($^{\circ}$)')
    ax[0].set_ylim(0,90)
    ax[0].set_xlim(0,90)

    ax[1].scatter(AcAcdata[0], AcAcdata[1], c='green', s=2.0)
    ax[1].scatter(43.4,5.5, color='gold', edgecolor='k', marker='s', label=r'S$_1$/S$_0$-CIpy')
    ax[1].scatter(0.8,0.3, color='gold', edgecolor='k', marker='D', label=r'S$_1$/S$_0$-Tw-MECI')
    ax[1].set_ylabel(r'PyrM ($^{\circ}$)')
    ax[1].set_xlabel(r'PyrC ($^{\circ}$)')
    ax[1].set_ylim(0,90)
    ax[1].set_xlim(0,90)

    fig.legend(loc='upper center', frameon=False)
    plt.tight_layout()
    plt.savefig(path_script + 'pyrC-vs-Pyrother.png', dpi=500)


    fig, ax = plt.subplots(2,1,figsize=(4,6), sharex=True)

    ax[0].scatter(MAbeforepyrc, MAbeforepyrM, s=2.0, c='green')
    ax[0].scatter(MAafterpyrc, MAafterpyrM, s=2.0, c='limegreen')
    ax[0].scatter(45.0,23.4,color='gold', edgecolor='k', marker='s', label=r'S$_1$/S$_0$-CIpy')
    ax[0].scatter(0.8,14.5,color='gold', edgecolor='k', marker='D', label=r'S$_1$/S$_0$-Tw-MECI')
    ax[0].set_ylabel(r'PyrM ($^{\circ}$)')
    ax[0].set_ylim(0,90)
    ax[0].set_xlim(0,90)

    ax[1].scatter(AcAcbeforepyrc, AcAcbeforepyrM, c='green', s=2.0)
    ax[1].scatter(AcAcafterpyrc, AcAcafterpyrM, c='limegreen', s=2.0)
    ax[1].scatter(43.4,5.5, color='gold', edgecolor='k', marker='s', label=r'S$_1$/S$_0$-CIpy')
    ax[1].scatter(0.8,0.3, color='gold', edgecolor='k', marker='D', label=r'S$_1$/S$_0$-Tw-MECI')
    ax[1].set_ylabel(r'PyrM ($^{\circ}$)')
    ax[1].set_xlabel(r'PyrC ($^{\circ}$)')
    ax[1].set_ylim(0,90)
    ax[1].set_xlim(0,90)

    fig.legend(loc='upper center', frameon=False)
    plt.tight_layout()
    plt.savefig(path_script + 'pyrC-vs-Pyrother-sep.png', dpi=500)

