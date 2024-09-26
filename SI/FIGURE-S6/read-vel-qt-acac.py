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
            KE += 0.5 * mC * (np.square(float(line[1])) + np.square(float(line[2])) + np.square(float(line[3])))
        if "O" in line:
            KE += 0.5 * mO * (np.square(float(line[1])) + np.square(float(line[2])) + np.square(float(line[3])))
        if "H" in line:
            KE += 0.5 * mH * (np.square(float(line[1])) + np.square(float(line[2])) + np.square(float(line[3])))

    return KE


if __name__ == "__main__":

    au_to_fs = 0.02418884254
    au_to_ev = 27.211399

    path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/test-script-geom/code/velocity-code/Kinetic-energy-final/"
    pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S6/"
    os.chdir(path)
    
    timstp = []

    list_file = ['pulse-total-0.24-vel.xyz', 'complete-IC-final-QT-velocities.xyz', 'complete-QT-traj-combined-velocities.xyz']

    for file in list_file:
        # break up the trajectory in individual xyz files
        fdata = open(str(file), 'r')   # use your data file

        num = 15 # number of atoms for this particular system # please change for your system

        before = []
        KE_list = []

        os.system('mkdir tmp')
        for line in fdata:
            before.append(line)
            if len(before) > 2:
                before.pop(0)
            if "Time step:" in line:
                t = line.split()
                l = t[2]                    # 7th column gives the actual timestep 
        #l = float(l) * au_to_fs    # timestep converted to fs
                timstp.append(float(l)*20*au_to_fs)           # timestep appended
                dat = 'tmp/' + str(l) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                KE_list.append(calc_KE(dat)*au_to_ev)

        os.system('rm -rf tmp')
        fdata.close()


        #KE_dist = open('final-KE-dist.dat','w')
        #for i,j in zip(timstp, KE_list):
        #    KE_dist.write('{:.2f}\t{:.10f}\n'.format(float(i),j))

        #KE_dist.close()
        print(np.mean(KE_list))

        from scipy.stats import norm
        # normal distribution
        mu, std = norm.fit(KE_list)

        plt.figure()
        plt.hist(KE_list, bins=30, range=[0.0,4.5], density=True, weights=None, color = 'orange', alpha = 0.9)
        plt.axvline(x=1.675, color = 'black', linestyle='--', label='0.5*ZPE')

        # Plot the PDF.
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)

        #plt.plot(x, p, 'k', linewidth=2)
        #title = "mean and stdev: {:.2f} and {:.2f}".format(mu, std)
        #plt.title(title)
        print(mu, std)
        plt.xlabel('Kinetic Energy (eV)', fontsize=18)
        plt.ylabel('Normalized Distribution', fontsize=18)
        plt.xlim(0.0,4.5)
        plt.ylim(0,1.2)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)

        
        plt.legend(fontsize="18", loc ="upper right", frameon=False)
        plt.tight_layout()
        plt.savefig(pathscript + str(file) + '-acac.pdf',dpi=300)

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

    os.chdir(pathscript)
