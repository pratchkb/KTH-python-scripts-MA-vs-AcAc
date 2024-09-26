import numpy as np
import os
import glob
# load proper modules
import get_optim
import geom_param
from itertools import islice
import matplotlib.pyplot as plt

def calc_KE(filename):
#    """
#    """

    import numpy as np
    amber_to_autime = (1.E3 / 20.455) * 41.341374575751
    au_to_bohr = 1.88973

    amber_to_au = au_to_bohr/amber_to_autime

    data = np.genfromtxt(filename, skip_header=2, dtype='unicode')
    data[:,1:] = data[:,1:].astype(np.float64)

    # masses in atomic units
    mH = 1822.88696289
    mC = 21874.64355469
    mO = 29166.19140625

    KE = 0
    # velocities in atomic units
    for line in data:
        if "C" in  line:
            KE += 0.5 * mC * (np.square(float(line[1])*amber_to_au) + np.square(float(line[2])*amber_to_au) + np.square(float(line[3])*amber_to_au))
        if "O" in line:
            KE += 0.5 * mO * (np.square(float(line[1])*amber_to_au) + np.square(float(line[2])*amber_to_au) + np.square(float(line[3])*amber_to_au))
        if "H" in line:
            KE += 0.5 * mH * (np.square(float(line[1])*amber_to_au) + np.square(float(line[2])*amber_to_au) + np.square(float(line[3])*amber_to_au))

    return KE


if __name__ == "__main__":

    au_to_fs = 0.02418884254
    au_to_ev = 27.211399

    timstp = []
    KE_list = []
    tim = []
    KE_list_12 = []
    avgCH = []

    CH_pairs = [[2,5],[3,6],[4,7]]
    num = 9 # number of atoms for this particular system # please change for your system

    path = "/data/projects/Pratip_MA_vs_AcAc/MA-data/S0-min-b3lyp-d3/WIGNER/wigner_298.15K_IC5000/"
    pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S6/"

    for i in range(0,5000,1):
        i = str(i).zfill(4)
        timstp.append(i)

        filename = path + 'v' + i + '.xyz'

        KE_list.append(calc_KE(filename)*au_to_ev)


        data = get_optim.read_xyz_as_np(path + 'x' + str(i) + '.xyz')        # coordinates loaded
        # average of all CH bonds
        CH = []
        for CH_pair in CH_pairs:
            CH.append(geom_param.compute_bond(data, CH_pair))

        avgCH.append(np.average(CH))

        if np.average(CH) < 1.2:

            tim.append(i)

            filename = path + 'v' + i + '.xyz'
            KE_list_12.append(calc_KE(filename)*au_to_ev)

    #KE_dist = open('final-KE-dist.dat','w')
    #for i,j in zip(timstp, KE_list):
    #    KE_dist.write('{:.2f}\t{:.10f}\n'.format(float(i),j))

    #KE_dist.close()

    from scipy.stats import norm
    # normal distribution
    mu, std = norm.fit(KE_list)
    mu2, std2 = norm.fit(KE_list_12)

    plt.figure()
    plt.hist(KE_list, bins=30, range=[0,3.5], density=True, weights=None, color = 'orange', alpha = 0.9)
    #plt.hist(KE_list_12, bins=30, range=[0,4.5], density=True, weights=None, color = 'b', histtype='step', stacked=True, linewidth = 3.0, label='Avg. r(CH) < 1.2 $\AA$', fill=False)
    plt.axvline(x=0.9235, color = 'black', linestyle='--', label='0.5*ZPE')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    p2 = norm.pdf(x, mu2, std2)

    #plt.plot(x, p, 'k', linewidth=2)

    #plt.plot(x, p2, 'g', linewidth=2)
    #title = "means and stdevs: {:.2f}, {:.2f} and {:.2f}, {:.2f}".format(mu, std, mu2, std2)
    #plt.title(title)

    print(mu, std, mu2, std2)
    plt.xlabel('Kinetic Energy (eV)',fontsize=18)
    plt.ylabel('Normalized Distribution', fontsize=18)
    plt.xlim(0.0,3.5)
    plt.ylim(0,1.8)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)


    plt.legend(fontsize="18", loc ="upper right", frameon=False)

    plt.tight_layout()
    plt.savefig(pathscript + 'test-KE-wigner-298.15K.pdf',dpi=300)

#data = np.genfromtxt('final-torsion-classical.dat', dtype='float64')
#data = data[10335:,5]   # from 5 ps onwards for equlibration

#plt.figure()
#plt.hist(tmin_new, bins=45, range=[-90,90], density=True, weights=None)
#plt.axvline(x=30, color = 'black')
#plt.axvline(x=-30, color = 'black')
#plt.xticks([-90,-70,-50,-30,-10,10,30,50,70,90])
#plt.xlabel(r'Min torsion (H$_7$C$_6$C$_{10}$[H$_{11}$,H$_{12}$,H$_{13}$])$\degree$')
#plt.ylabel('Norm. dist.')
#plt.savefig('test-quantum-all.png',dpi=300)

