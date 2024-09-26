import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.colors as colors
#from matplotlib.ticker import MultipleLocator
#from scipy import ndimage
import get_optim
import geom_param
from itertools import islice
import matplotlib.pyplot as plt

path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/2D-plot-data/BAGEL-lint-ecup-symCI1412-ecupreversed/rotate-methyl/"
pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/MAIN-PAPER-FIGURE/FIGURE-4/ecup-stag/"

def read_XMSCAST2_results(nliic,totrot,outname='bagel.out'):
	
    liic = []
    deg = []
    result_S0 = []
    result_S1 = []
    result_S2 = []

    before = []
    timstp = []

    Hbond1 = []
    Hbond2 = []

    htran = []

    #Iterate through the LIIC of H-transfer
    for iliic in range(0,nliic,1):

        #Iterate through methyl rotations
        for ideg in range(0,totrot,10):

            # Grab excited state energies and oscillator strengths
            results_file = os.getcwd() + '/' + str(iliic) + '/' + str(ideg) + '/' + str(outname)  
            with open(results_file, 'r') as f_in:
                for line in f_in.readlines():

                    if "MS-CASPT2 energy : state  0" in line:
                        state1, energy1 = int(line.split()[5]),float(line.split()[6]) 
                        result_S0.append(energy1)
                        liic.append(iliic)
                        deg.append(ideg)
 
                    if "MS-CASPT2 energy : state  1" in line:
                        state2, energy2 = int(line.split()[5]),float(line.split()[6])
                        result_S1.append(energy2)

                    if "MS-CASPT2 energy : state  2" in line:
                        state3, energy3 = int(line.split()[5]),float(line.split()[6])
                        result_S2.append(energy3)

            # break up the trajectory in individual xyz files
            filename = os.getcwd() + '/' + str(iliic) + '/' + str(ideg) + '/' + str(ideg) + '.xyz'
            fdata = open(filename, 'r')   # use your data file

            num = 15 # number of atoms for this particular system # please change for your system

            os.system('mkdir tmp')
            for line in fdata:
                #print(line)
                before.append(line)
                if len(before) > 2:
                    before.pop(0)
                if "Angle " in line:
                    t = line.split()
                    t = t[1]
            #t = float(t)     # timestep converted to fs
        #t = t[1].replace(',','')                    # 7th column gives the actual timestep
                    timstp.append(t)           # timestep appended
                    dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
                    f=open(dat,'w')
                    f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                    f.close()
                    data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')
       #O-H distances
                    Hbond1.append(geom_param.compute_bond(data,[1,13]))
                    Hbond2.append(geom_param.compute_bond(data,[13,4]))

                    htran.append(geom_param.compute_bond(data,[1,13])-geom_param.compute_bond(data,[13,4]))

            fdata.close()
            os.system('rm -r tmp')


    print(len(htran))
    f=open(pathscript + 'XMS-CASPT2-Htransfer-methrotation-ecup-updated.dat','w')
    for i,j,k,l,m,n in zip(liic, deg, result_S0, result_S1, result_S2, htran):
        f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(i,j,k,l,m,n))
    f.close()

os.chdir(path)
#call function to hrt the data file 
read_XMSCAST2_results(21,61)
os.chdir(pathscript)

#plot
datatot = np.genfromtxt(pathscript + 'XMS-CASPT2-Htransfer-methrotation-ecup-updated.dat', dtype='float64')

au_per_ev = 3.6749330610942925E-02
ev_per_au = 1 / au_per_ev

#assign reference value to E(S0_min)
min_e = -344.8015591313 # using the staggered energy  #datatot[0,2]

E_S0 = datatot[:,2] 
E_S1 = datatot[:,3]
E_S2 = datatot[:,4]


s0 = 0
s1 = 1
s2 = 2

def get_e(data):
    new = []

    for e in data:
        new.append((e - min_e) * ev_per_au)

    return new    


def get_de():

    de = []
    
    for e1,e2 in zip(datatot[:,3],datatot[:,4]):
        de.append((e2 - e1) * ev_per_au)

    return de
 
#    if state1 > state2: state1, state2 = state2, state1
#    return np.array([e[state2]-e[state1] for e  in d['energies']]) * ev_per_au

def get_levels(state):
    if state == s2:
        return np.arange(4.3, 4.9, 0.02) 
    elif state == s1:
        return np.arange(4.1, 4.6, 0.02)
    else:
        return np.arange(0, 0.5, 0.02)


# each state surface
x = datatot[:,0]
y = datatot[:,1]
for state in [s0, s1, s2]:
     #z = get_e(data, state)
     if state == s2:
     	z = get_e(E_S2)
     elif state == s1:
        z = get_e(E_S1)
     else:
        z = get_e(E_S0)
     fig = plt.figure(figsize=(4,3))
     #plt.suptitle(f'S{state}')
     plt.tricontourf(x, y, z, cmap='gist_heat_r', levels=get_levels(state))
     cbar = plt.colorbar()
     cbar.ax.tick_params(labelsize=12)
     cbar.set_label("Energy (eV)", fontsize=12)
     plt.gca().xaxis.set_tick_params(labelsize=12)
     plt.gca().yaxis.set_tick_params(labelsize=12)
     plt.xlim([0, 20])
     plt.ylim([0, 60])
     plt.xticks([0,4,8,12,16,20])
     plt.yticks([0,10,20,30,40,50,60])

     plt.xlabel(r"H-transfer", fontsize=12)
     plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
     plt.subplots_adjust(left=0.24, bottom=0.16, hspace=0)
     plt.savefig(f"energy_surf_S{state}.pdf")


def geom_param_all(filename):

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = 15 # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    BLA = []
    angles = []

    Hbond1 = []
    Hbond2 = []

    torsion1 = []
    torsion2 = []

    htran = []
    pyr = []

    os.system('mkdir tmp')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "Frame " in line:
            t = line.split()
            t = t[1]
            #t = float(t)     # timestep converted to fs
        #t = t[1].replace(',','')                    # 7th column gives the actual timestep
            timstp.append(t)           # timestep appended
            dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
            f=open(dat,'w')
            f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
            f.close()
            data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

        #compute pyramidalization
            #pyr.append(geom_param.compute_pyramidalization(data,2,0,3,14))

        #compute sum of angles
            angles.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))

        #O-H distances
            Hbond1.append(geom_param.compute_bond(data,[13,4]))
            Hbond2.append(geom_param.compute_bond(data,[1,13]))

            htran.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))

        #torsion
            torsion1.append(geom_param.compute_torsion(data,1,0,2,3))
            torsion2.append(geom_param.compute_torsion(data,0,2,3,4))


    fdata.close()
    os.system('rm -r tmp')

    return htran



plt.cla()
# S2 - S1 energy diff surface
x = datatot[:,0]
y = datatot[:,1]
#print(x)
#print(y)
z = get_de()
#plt.suptitle(f'S{state}S{state-1}')
plt.tricontourf(x, y, z, levels=np.arange(0,0.51,0.02),  cmap='gist_heat_r' )
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=12)
cbar.set_label("Energy difference (eV)", fontsize=12)
plt.gca().xaxis.set_tick_params(labelsize=12)
plt.gca().yaxis.set_tick_params(labelsize=12)
plt.xlim([0, 20])
plt.ylim([0, 60])
plt.xticks([0,4,8,12,16,20])
plt.yticks([0,10,20,30,40,50,60])
plt.xlabel(r"H-transfer", fontsize=12)
plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
plt.subplots_adjust(left=0.24, bottom=0.16, hspace=0)
plt.savefig(pathscript + "energydiff_surf_S2S1_ecup_stag.pdf", dpi=1000)

  
plt.cla()
# S2 - S1 energy diff surface
x = datatot[:,5]
y = datatot[:,1]

maxx = np.max(datatot[:,5])
minx = np.min(datatot[:,5])
#print(x)
#print(y)
z = get_de()
#print(z)
fig = plt.figure(figsize=(4,3))
#plt.suptitle(f'S{state}S{state-1}')
plt.tricontourf(x, y, z, levels=np.arange(0,0.51,0.02),  cmap='gist_heat_r' )
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=12)
cbar.set_label("Energy difference (eV)", fontsize=12)
plt.gca().xaxis.set_tick_params(labelsize=12)
plt.gca().yaxis.set_tick_params(labelsize=12)
plt.xlim([minx, maxx])
plt.ylim([0, 60])
plt.minorticks_on()
#plt.xticks([0,4,8,12,16,20])
plt.yticks([0,10,20,30,40,50,60])
plt.xlabel(r"H-transfer ($\AA$)", fontsize=12)
plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
plt.subplots_adjust(left=0.24, bottom=0.16, hspace=0)
plt.savefig(pathscript + "energydiff_surf_S2S1_HT_ecup_stag.pdf", dpi=1000)

