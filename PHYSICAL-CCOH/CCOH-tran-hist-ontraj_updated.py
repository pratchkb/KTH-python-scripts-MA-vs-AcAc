import numpy as np
import matplotlib.pyplot as plt
import os

import get_optim
import geom_param

from itertools import islice

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as colors
from scipy import ndimage


def calculate_distance(atom1_coord, atom2_coord):
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12


def calculate_bla4(data,atom_order):
    """     
    Calculates all 3 types of BLA from the xyz data file and the list of atom order
    """     
            
    BLA4 = calculate_distance(data[atom_order[0]],data[atom_order[1]]) - calculate_distance(data[atom_order[1]],data[atom_order[2]]) + calculate_distance(data[atom_order[2]],data[atom_order[3]]) - calculate_distance(data[atom_order[3]],data[atom_order[4]])
    BLACO = calculate_distance(data[atom_order[0]],data[atom_order[1]]) - calculate_distance(data[atom_order[3]],data[atom_order[4]])
    BLACC = - calculate_distance(data[atom_order[1]],data[atom_order[2]]) + calculate_distance(data[atom_order[2]],data[atom_order[3]])
            
    return BLA4, BLACO, BLACC


def collect_geom_param(filename, atom_order, traj, natoms):

    num = natoms # number of atoms for this particular system # please change for your system

    data = get_optim.read_xyzs_as_np(filename)
    
    for i in range(len(data[:])):
        s = calculate_bla4(data[0],atom_order)
        bla = s[0]

    return bla


def calculate_ccoh_torsion(filename, molecule, natoms, hh1, hh2):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

    data = get_optim.read_xyzs_as_np(filename)
    tim = []
    tor = []    

    for i in range(len(data[:])):

        h1 = geom_param.compute_bond(data[i],hh2)
        h2 = geom_param.compute_bond(data[i],hh1)

        if molecule == "MA":
            if h1 < h2:
                tim.append(i*0.5)
                tor.append(abs(geom_param.compute_torsion(data[i],3,4,1,8)))
            else:
                tim.append(i*0.5)
                tor.append(abs(geom_param.compute_torsion(data[i],3,2,0,8)))

        if molecule == "AcAc":
            if h1 < h2:
                tim.append(i*0.5)
                tor.append(abs(geom_param.compute_torsion(data[i],2,3,4,13)))
            else:
                tim.append(i*0.5)
                tor.append(abs(geom_param.compute_torsion(data[i],2,0,1,13)))

    return tim, tor


def calculate_param_statewise(filename, atom_order, natoms):

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstps2 = []
    timstps1 = []

    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
            t = line.split()
            ti = t[1]
            st = t[2]
            #print(ti)
            if int(st) == 3:         # S2 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps2.append(l)           # timestep appended

            if int(st) == 2 or int(st) == 1:         # S1 or S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps1.append(l)           # timestep appended


    return timstps2, timstps1


def calculate_param_nostate(filename, atom_order, traj, natoms):

    import get_optim
    import geom_param
    
    num = natoms # number of atoms for this particular system # please change for your system

    time = []
    bla = []

    data = get_optim.read_xyzs_as_np(filename)
    time = []
    bla = []
    for i in range(len(data[:])):
        time.append(i*0.5)
        s = calculate_bla4(data[i],atom_order)
        bla.append(s[0])

    return time, bla, data

def collect_geom_param_hops(filename, atom_order, traj, natoms):

    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    tim = []

    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
            t = line.split()
            if traj == t[0][-10:]:        # checking traj number
                t = t[2]
                tim.append(float(t)) #* au_to_fs    # timestep converted to fs

    fdata.close()

    return tim

if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/PHYSICAL-CCOH/'
    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        #plots
        if molecule == "MA":
    	    hh1 = [0,8]; hh2 = [1,8]

        if molecule == "AcAc":
            hh1 = [1,13]; hh2 = [4,13]

        #print(hh1, hh2)
        #print(ntraj)
        #print(run_path)
        #print(numatoms)
     
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=(6,3))
        gs = GridSpec(1,1)
        clr = ['#0072BD', 'orange', '#66BF2E', 'k', 'green', 'darkgreen', 'darkblue', 'orangered', '#D45500', 'crimson']
        ax1 = fig.add_subplot(gs[0])
        #ax2 = fig.add_subplot(gs[1])


        os.chdir(run_path)
        #for21 = collect_geom_param(run_path + 'crossing-s2-s1.xyz', numatoms)
        #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms)
        #back21 = collect_geom_param(run_path + 'crossing-s1-s2.xyz', numatoms)

        HT1 = []
        HT2 = []
        FC = []

        
        t21 = []
        s21 = []
        t12 = []
        s12 = []

        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir('TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):
                os.chdir(run_path)
            else:
            	#densities = calculate_param_statewise_color("output.xyz", atom_order, 15)
                densities = calculate_ccoh_torsion("output.xyz", molecule, numatoms, hh1, hh2)
                pos = np.where(np.abs(np.diff(densities[0])) > 0.5)[0] + 1
                time = np.insert(densities[0], pos, np.nan)
                torr = np.insert(densities[1], pos, np.nan)
                #print(time,ht)
                ax1.plot(time, torr, color = clr[1], alpha=0.9, linewidth=0.4, zorder=0)
 
                #times21 = collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz',atom_order,'TRAJ_'+str(i),numatoms)
                #print(for21)
                #times10 = collect_geom_param_hops(run_path + 'crossing-s1-s0.xyz',atom_order,'TRAJ_'+str(i),numatoms)
                #print(for10)
                #times12 = collect_geom_param_hops(run_path + 'crossing-s1-s2.xyz',atom_order,'TRAJ_'+str(i),numatoms)

#                tim21 = []
#                bla21 = []
#                for time in times21:
#                    frame = int(time*2)
#                    tim21.append(time)
#                    bla21.append(calculate_bla4(densities[2][frame],atom_order)[0])
#                    t21.append(time)
#                    s21.append(calculate_bla4(densities[2][frame],atom_order)[0])
#
#
#                tim12 = []
#                bla12 = []
#                for time in times12:
#                    frame = int(time*2)
#                    tim12.append(time)
#                    bla12.append(calculate_bla4(densities[2][frame],atom_order)[0])
#                    t12.append(time)
#                    s12.append(calculate_bla4(densities[2][frame],atom_order)[0])
#
#
#
#                densit = calculate_param_statewise("output.xyz", atom_order, numatoms)
#                tim2 = []
#                bla2 = []
#                for time in densit[0]:
#                    tim2.append(time)
#                    bla2.append(calculate_bla4(densities[2][int(time*2.0)],atom_order)[0])
#                for elem in bla2:
#                    HT1.append(elem) 
#
#                #densit = calculate_param_statewise("output.xyz", atom_order, numatoms)
#                tim1 = []
#                bla1 = []
#                for time in densit[1]:
#                    tim1.append(time)
#                    bla1.append(calculate_bla4(densities[2][int(time*2.0)],atom_order)[0])                                
#                for elem in bla1[:10]:      # first 5 fs
#                    HT2.append(elem) 
#
#                IC = collect_geom_param("output.xyz", atom_order,'TRAJ_'+str(i), numatoms)                
#                FC.append(IC)
                
                os.chdir(run_path)    


#        ax2.hist(FC, bins=60, range=[-1,1], density=True, weights=None, orientation='horizontal', color = 'orange', alpha = 0.9, zorder=10)
#        ax2.hist(HT1, bins=60, range=[-1,1], density=True, histtype='step', orientation='horizontal', weights=None, color = 'blue', alpha = 0.9, zorder=11)
#        ax2.hist(HT2, bins=60, range=[-1,1], density=True, histtype='step', orientation='horizontal', weights=None, color = 'limegreen', alpha = 0.9, zorder=12)

        os.chdir(path_script)
        #ax1.scatter(tim21, bla21, color = clr[8], s=0.6, label=r'S$_2$/S$_1$ hops', zorder=60)
    	#ax1.scatter(for10[0], for10[1], s=0.6, color = clr[5], label=r'S$_1$/S$_0$ hops', zorder=61)
        #ax1.scatter(tim12, bla12, s=0.6, color = clr[6], label=r'S$_1$/S$_2$ hops', zorder=62)
#        if molecule == "MA":
#            ax1.axhline(y=-0.151, color='k', linestyle='dashed', lw = 0.7, label='FC', zorder=48)
        #    ax1.axhline(y=0.661, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=49)
#            ax2.axhline(y=-0.151, color='k', linestyle='dashed', lw = 0.7, label='FC',zorder=50)
        #    ax2.axhline(y=0.661, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=51)
#        if molecule == "AcAc":
#            ax1.axhline(y=-0.150, color='k', linestyle='dashed', lw = 0.7, label='FC', zorder=48)
        #    ax1.axhline(y=0.591, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=49)
#            ax2.axhline(y=-0.150, color='k', linestyle='dashed', lw = 0.7, label='FC',zorder=50)
        #    ax2.axhline(y=0.591, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=51)

        ax1.set_ylabel(r"Torsion CCOH ($\AA$)")
        ax1.set_xlabel(r"Time (fs)")
        ax1.set_xlim([0,200])
        #ax1.set_ylim([-0.5,0.5])
        #ax2.set_ylim([-0.5,0.5])
        #ax1.set_xticks([0,25,50,75,100,125,150,175])
        #ax2.set_xlim([0.0,6.0])
        #ax2.set_xticks([0.0,2.0,4.0,6.0])
        #ax2.set_yticks([])
        ax1.minorticks_on()
        #ax2.minorticks_on()
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.05)
        plt.savefig(path_script + f'{molecule}-torccoh-density-colored-updated.png', dpi=500)
        plt.savefig(path_script + f'{molecule}-torccoh-density-colored-updated.pdf', dpi=500)

