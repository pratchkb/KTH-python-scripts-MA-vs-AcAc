import numpy as np
import matplotlib.pyplot as plt
import os
import get_optim
import geom_param
from itertools import islice

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


#def collect_geom_param(filename, atom_order, traj, natoms):
#
#    num = natoms # number of atoms for this particular system # please change for your system
#
#    data = get_optim.read_xyzs_as_np(filename)
#
#    for i in range(len(data[:])):
#        s = calculate_bla4(data[0],atom_order)
#        bla = s[0]
#
#    return bla

def calculate_order_atoms(filename, molecule, natoms, hh1, hh2):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

    data = get_optim.read_xyzs_as_np(filename)
    data = data[0]
    h1 = geom_param.compute_bond(data,hh2)
    h2 = geom_param.compute_bond(data,hh1)

    if molecule == "MA":
        if h1 < h2:
            order = [0,2,3,4,1]
        else:
            order = [1,4,3,2,0]

    if molecule == "AcAc":
        if h1 < h2:
            order = [1,0,2,3,4]
        else:
            order = [4,3,2,0,1]

    return order


def calculate_sort(filename, natoms):

    from itertools import islice

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstps2 = []
    timstps1 = []
    timstps0 = []

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
            
            if int(st) == 2:         # S1 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps1.append(l)           # timestep appended

            if int(st) == 1:         # S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps0.append(l)           # timestep appended
    
    #os.chdir(path)
    return timstps2, timstps1, timstps0

def collect_geom_param_hops(filename, atom_orders, traj, natoms):

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

def calculate_param_nostate(filename, atom_order, traj, natoms, hh1, hh2, tim):

    num = natoms # number of atoms for this particular system # please change for your system
    htran = []
    bla = []
    time = []

    data = get_optim.read_xyzs_as_np(filename)
    for i in tim:
        t = int(i*2)
        time.append(t)
        h1 = geom_param.compute_bond(data[t],hh2)
        h2 = geom_param.compute_bond(data[t],hh1)
        htran.append(h1-h2)
        bla.append(calculate_bla4(data[t],atom_order)[0])

    return time, htran, bla


def collect_geom_param_hops(filename, trajname, atom_order, traj, natoms, hh1, hh2):

    data = get_optim.read_xyzs_as_np(trajname)
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    tim = []
    ht = []
    bla = []

    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
            t = line.split()
            if traj == t[0][-10:]:        # checking traj number
                t = float(t[2])
                tim.append(t) #* au_to_fs    # timestep converted to fs 
                frame = int(float(t)*2)

                bla.append(calculate_bla4(data[frame],atom_order)[0])
                ht.append(geom_param.compute_bond(data[frame],hh2)-geom_param.compute_bond(data[frame],hh1))



    return tim, ht, bla



if __name__ == "__main__":


    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    #Ntraj = [20,20]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    fig, ax = plt.subplots(2,1, figsize=(5,8), sharex=True)

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S17/'

    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        if molecule == "MA":
            hh1 = [0,8]; hh2 = [1,8]

        if molecule == "AcAc":
            hh1 = [1,13]; hh2 = [4,13]

        os.chdir(run_path)
        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir('TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):

                #print(str(i))
                os.chdir(run_path)

            else:

                # calling functions
                #timeall = calculate_sort("output.xyz", numatoms)
                #s2 = timeall[0][-1]
                #s2 = f"{s2}{'0'*4}"
                #s1 = timeall[1][-1]
                #print(s1) 
                #s1 = f"{s1}{'0'*4}"
                orders = calculate_order_atoms("output.xyz", molecule, numatoms, hh1, hh2) 

                for21 = collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz',run_path + 'TRAJ_'+str(i) + '/output.xyz', orders,'TRAJ_'+str(i),numatoms,hh1,hh2)
                #print(for21)
                for10 = collect_geom_param_hops(run_path + 'crossing-s1-s0.xyz',run_path + 'TRAJ_'+str(i) + '/output.xyz',orders,'TRAJ_'+str(i),numatoms,hh1,hh2)
                #print(for10)
                back21 = collect_geom_param_hops(run_path + 'crossing-s1-s2.xyz', run_path + 'TRAJ_'+str(i) + '/output.xyz', orders,'TRAJ_'+str(i),numatoms,hh1,hh2)

				#### tor is hydrogen-transfer, pyr is BLA in the next bit

                tim21 = []
                ht21 = []
                bla21 = []
                for elem,torelem,pyrelem in zip(for21[0],for21[1],for21[2]):
                    tim21.append(elem)
                    ht21.append(abs(torelem))
                    bla21.append(abs(pyrelem))
            
                timless = []
                htless = []
                blaless = []
                timbig = []
                htbig = []
                blabig = []
             
                threshold = 75
                #thfor10 = threshold_hop(for10[0],for10[1],for10[2],75)
                for elem,torelem,pyrelem in zip(for10[0],for10[1],for10[2]):
                    if elem < threshold:
                        timless.append(elem)
                        htless.append(abs(torelem))
                        blaless.append(abs(pyrelem))
                    else:
                        timbig.append(elem)
                        htbig.append(abs(torelem))
                        blabig.append(abs(pyrelem))

                if molecule == "MA":
                    if for21[0] != 0:
                        ax[0].scatter(ht21, bla21, s=15.0, c='#D45500', zorder=5)
                    if for10[0] != 0:
                        ax[0].scatter(htless,blaless, s=15.0, c='darkgreen', zorder=15)
                        ax[0].scatter(htbig,blabig, s=15.0, c='limegreen', zorder=16)

                if molecule == "AcAc":
                    if for21[0] != 0:
                        ax[1].scatter(ht21, bla21, s=15.0, c='#D45500', zorder=5)
                    if for10[0] != 0:
                        ax[1].scatter(htless,blaless, s=15.0, c='darkgreen', zorder=15)
                        ax[1].scatter(htbig,blabig, s=15.0, c='limegreen', zorder=16)
                 
                os.chdir(run_path)


    ax[0].set_ylabel(r"|BLA| ($\AA$)")
    ax[1].set_ylabel(r"|BLA| ($\AA$)")
    #ax[0].set_xlabel(r"Torsion ($^{\circ}$)")
    ax[1].set_xlabel(r"|HT| ($\AA$)") 
    ax[0].scatter(2.456,0.131, marker='D', label=r'S$_1$/S$_0$-Tw', color='gold', edgecolors='k', zorder=20)
    ax[0].scatter(0.898,0.102, marker='s', label=r'S$_1$/S$_0$-CI$_\text{pyrC}$', color='gold',edgecolors='k', zorder=21)       
    ax[0].scatter(0.001,0.000, marker='*', s=100.0, label=r'S$_2$/S$_1$-sym-HTI', color='gold',edgecolors='k', zorder=22)       
    ax[0].scatter(2.956,0.020, marker='P', s=60.0, label=r'S$_2$/S$_1$-Tw', color='gold',edgecolors='k', zorder=22)       
    ax[1].scatter(2.344,0.143, marker='D', label=r'S$_1$/S$_0$-Tw', color='gold', edgecolors='k', zorder=23)
    ax[1].scatter(0.711,0.049, marker='s', label=r'S$_1$/S$_0$-CI$_\text{pyrC}$', color='gold',edgecolors='k',zorder=24)
    ax[1].scatter(0.001,0.000, marker='*', s=100.0, label=r'S$_2$/S$_1$-sym-HTI', color='gold',edgecolors='k', zorder=25)       
    ax[1].scatter(2.134,0.017, marker='P', s=60.0, label=r'S$_2$/S$_1$-Tw', color='gold',edgecolors='k', zorder=22)       
    ax[0].set_xlim(-0.1,3.5)
    ax[0].set_ylim(-0.01,0.5)
    ax[1].set_xlim(-0.1,3.5)
    ax[1].set_ylim(-0.01,0.5)
    #ax[0].set_yticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0])
    #ax[0].set_xticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0,105.0,120.0])
    #ax[1].set_yticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0])
    #ax[1].set_xticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0,105.0,120.0])
    ax[0].minorticks_on()
    ax[1].minorticks_on()

    ax[0].legend(frameon=False, loc = 'upper left', ncols=2)
    ax[1].legend(frameon=False, loc = 'upper left', ncols=2)
    plt.tight_layout()
    fig.savefig(path_script + 'HT-vs-BLA.png', dpi=500)
