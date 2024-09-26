import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.colors as colors
#from matplotlib.ticker import MultipleLocator
#from scipy import ndimage

path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/2D-plot-data/BAGEL-lint-ecup-symCI1412-ecupreversed/rotate-methyl-Edown/"
pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/MAIN-PAPER-FIGURE/FIGURE-4/stag_edown/"


def read_XMSCAST2_results(nliic,totrot,outname='bagel.out'):
	
    liic = []
    deg = []
    result_S0 = []
    result_S1 = []
    result_S2 = []

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

    f=open(pathscript + 'XMS-CASPT2-Htransfer-methrotation-estag-updated.dat','w')
    for i,j,k,l,m in zip(liic, deg, result_S0, result_S1, result_S2):
        f.write('{}\t{}\t{}\t{}\t{}\n'.format(i,j,k,l,m))
    f.close()

os.chdir(path)
#call function to hrt the data file 
read_XMSCAST2_results(21,61)
os.chdir(pathscript)

#plot
datatot = np.genfromtxt(pathscript + 'XMS-CASPT2-Htransfer-methrotation-estag-updated.dat', dtype='float64')

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

plt.cla()
# S2 - S1 energy diff surface
x = datatot[:,0]
y = datatot[:,1]
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
plt.xlim([0, 20])
plt.ylim([0, 60])
plt.xticks([0,4,8,12,16,20])
plt.yticks([0,10,20,30,40,50,60])
plt.xlabel(r"H-transfer", fontsize=12)
plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
plt.subplots_adjust(left=0.24, bottom=0.16, hspace=0)
plt.savefig(f"energydiff_surf_S2S1.pdf", dpi=1000)

  

# plot 3 states together 
fig, ax = plt.subplots(1,1, subplot_kw={'projection': '3d'})
z2 = get_e(E_S2)
z1 = get_e(E_S1)
z0 = get_e(E_S0)

# ax.tricontourf(x, y, z1, levels=14,  cmap='gist_heat' )
surf = ax.plot_trisurf(x, y, z2, cmap='gist_heat', linewidth=0.1, vmin=0, vmax=5)
ax.plot_trisurf(x, y, z1, cmap='gist_heat', linewidth=0.1, vmin=0, vmax=5)
ax.plot_trisurf(x, y, z0, cmap='gist_heat', linewidth=0.1, vmin=0, vmax=5)
cbar = fig.colorbar(surf, ax=ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label("Energy (eV)", fontsize=12)
ax.view_init(elev=20, azim=-60)
plt.xlim([0, 20])
plt.ylim([0, 60])
plt.xticks([0,4,8,12,16,20])
plt.yticks([0,10,20,30,40,50,60])
plt.xlabel(r"H-transfer", fontsize=12)
plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
plt.savefig("energy_allsurf.pdf")

# plot s1 and s2 together
fig, ax = plt.subplots(1,1, subplot_kw={'projection': '3d'})
surf = ax.plot_trisurf(x, y, z2, cmap='gist_heat', linewidth=0.1, vmin=0, vmax=5)
ax.plot_trisurf(x, y, z1, cmap='gist_heat', linewidth=0.1, vmin=0, vmax=5)
cbar = fig.colorbar(surf, ax=ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label("Energy (eV)", fontsize=12)
ax.view_init(elev=20, azim=-60)
plt.xlim([0, 20])
plt.ylim([0, 60])
plt.xticks([0,4,8,12,16,20])
plt.yticks([0,10,20,30,40,50,60])
plt.xlabel(r"H-transfer", fontsize=12)
plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
plt.savefig("S1_S2_energy_surf.pdf")

# plot s1 and s2 together diff color
fig, ax = plt.subplots(1,1, subplot_kw={'projection': '3d'})
ax.plot_trisurf(x, y, z2, cmap='gist_heat', linewidth=0.1)
ax.plot_trisurf(x, y, z1, cmap='gist_gray', linewidth=0.1)
ax.view_init(elev=20, azim=-60)
plt.xlim([0, 20])
plt.ylim([0, 60])
plt.xticks([0,4,8,12,16,20])
plt.yticks([0,10,20,30,40,50,60])
plt.xlabel(r"H-transfer", fontsize=12)
plt.ylabel(r"Methyl rotation ($\degree$)", fontsize=12)
plt.savefig("S1_S2_energy_surf_color.pdf")

