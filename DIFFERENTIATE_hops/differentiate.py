import numpy as np
import matplotlib.pyplot as plt
import os
#import sys
#sys.path.append('/Users/pratip/')


#from aimsrelated.tcutil.code.utils import units
#from aimsrelated.tcutil.code.atom_data import atom_data
import get_optim
import geom_param

from itertools import islice
from itertools import chain


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


def calculate_order_atoms(filename, hh1, hh2, molecule, natoms):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

    import get_optim
    import geom_param
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    order = []

    os.system('mkdir tmp')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=         0.00000" in line: # Initial condition geometry
             t = line.split()
             t = t[1]
             l = float(t) #* au_to_fs    # timestep converted to fs
             timstp.append(l)           # timestep appended
             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

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

    fdata.close()

    os.system('rm -r tmp')

    return order

def calculate_IC_BLA(filename, atom_order, natoms):
    """
    calculate all the BLAs along the trajectories depending on the numbering of atoms decided 
    """

    import numpy as np
    # load proper modules
    import get_optim
    import geom_param
    from itertools import islice

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    bla4 = []

    os.system('mkdir tmp2')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=         0.00000" in line:
             t = line.split()
             t = t[1]
             l = float(t) #* au_to_fs    # timestep converted to fs
             #timstp.append(l)           # timestep appended
             dat = 'tmp2/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp2/' + str(t) + '.xyz')        # coordinates loaded

             bla4 = calculate_bla4(data,atom_order)[0]

    fdata.close()

    os.system('rm -r tmp2')

    return bla4

def calculate_param(filename, molecule, natoms, hh1, hh2, A1, A2, A3, pyrc):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """
    ht = []
    SOA = []
    pyr = []

    data = get_optim.read_xyzs_as_np(filename)
    for i in range(len(data[:])):
        h1 = geom_param.compute_bond(data[i],hh2)
        h2 = geom_param.compute_bond(data[i],hh1)
        a1 = geom_param.compute_angle(data[i],*A1)
        a2 = geom_param.compute_angle(data[i],*A2)
        a3 = geom_param.compute_angle(data[i],*A3) 
        p = geom_param.compute_pyramidalization(data[i],*pyrc)
        
        ht.append(h1-h2)
        pyr.append(p)
        SOA.append(a1+a2+a3)

    return ht, SOA, pyr


def calculate_pyr_atoms(filename, molecule, natoms, pyrc):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """
    pyr = []
   
    data = get_optim.read_xyzs_as_np(filename)
    for i in range(len(data[:])):
        p = geom_param.compute_pyramidalization(data[i],*pyrc)
        pyr.append(p)

    return pyr


def calculate_main_torsion(filename, molecule, natoms, hh1, hh2):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

    data = get_optim.read_xyzs_as_np(filename)
    tor = []

    for i in range(len(data[:])):

        h1 = geom_param.compute_bond(data[i],hh2)
        h2 = geom_param.compute_bond(data[i],hh1)

        if molecule == "MA":
            if h1 < h2:
                tor.append(geom_param.compute_torsion(data[i],2,3,4,1))
            else:
                tor.append(geom_param.compute_torsion(data[i],4,3,2,0))

        if molecule == "AcAc":
            if h1 < h2:
                tor.append(geom_param.compute_torsion(data[i],0,2,3,4))
            else:
                tor.append(geom_param.compute_torsion(data[i],3,2,0,1))

    return tor


if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    #plots
    fig, ax = plt.subplots(3,2, figsize=(8,10))
   
    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/DIFFERENTIATE_hops/'

    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        
        #plots
        if molecule == "MA":
            hh1 = [0,8]; hh2 = [1,8]
            A1 = [0,2,3]
            A2 = [2,3,4]
            A3 = [3,4,1]
            pyrs = [3,2,4,6]
            os.chdir(run_path)
            #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms)
            S21 = calculate_param(run_path + 'crossing-s2-s1.xyz', molecule, numatoms, hh1, hh2, A1, A2, A3, pyrs)
            S12 = calculate_param(run_path + 'crossing-s1-s2.xyz', molecule, numatoms, hh1, hh2, A1, A2, A3, pyrs)
            S10 = calculate_param(run_path + 'crossing-s1-s0.xyz', molecule, numatoms, hh1, hh2, A1, A2, A3, pyrs)

            ax[0,0].hist(S21[0], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-3.0,3.0], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[0,0].hist(S12[0], bins=40, histtype='step', color='limegreen', alpha=0.9, linewidth=1.5, density=True, range=[-3.0,3.0], label=r'S$_1$$\rightarrow$S$_2$ hops')
            ax[0,0].hist(S10[0], bins=40, histtype='step', color='blue', alpha=0.9, density=True, range=[-3.0,3.0], label=r'S$_1$$\rightarrow$S$_0$ hops')

            ax[1,0].hist(S21[1], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[320.0,400.0], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[1,0].hist(S12[1], bins=40, histtype='step', color='limegreen', alpha=0.9, linewidth=1.5, density=True, range=[320.0,400.0], label=r'S$_1$$\rightarrow$S$_2$ hops')
            ax[1,0].hist(S10[1], bins=40, histtype='step', color='blue', alpha=0.9, density=True, range=[320.0,400.0], label=r'S$_1$$\rightarrow$S$_0$ hops')

            ax[2,0].hist(S21[2], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-90.0,90.0], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[2,0].hist(S12[2], bins=40, histtype='step', color='limegreen', alpha=0.9, linewidth=1.5, density=True, range=[-90.0,90.0], label=r'S$_1$$\rightarrow$S$_2$ hops')
            ax[2,0].hist(S10[2], bins=40, histtype='step', color='blue', alpha=0.9, density=True, range=[-90.0,90.0], label=r'S$_1$$\rightarrow$S$_0$ hops')

            os.chdir(run_path)


        if molecule == "AcAc":
            hh1 = [1,13]; hh2 = [4,13]
            A1 = [1,0,2]
            A2 = [0,2,3]
            A3 = [2,3,4]
            pyrs = [2,0,3,14]
            os.chdir(run_path)
            
            S21 = calculate_param(run_path + 'crossing-s2-s1.xyz', molecule, numatoms, hh1, hh2, A1, A2, A3, pyrs)
            S12 = calculate_param(run_path + 'crossing-s1-s2.xyz', molecule, numatoms, hh1, hh2, A1, A2, A3, pyrs)
            S10 = calculate_param(run_path + 'crossing-s1-s0.xyz', molecule, numatoms, hh1, hh2, A1, A2, A3, pyrs)

            ax[0,1].hist(S21[0], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-3.0,3.0], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[0,1].hist(S12[0], bins=40, histtype='step', color='limegreen', alpha=0.9, linewidth=1.5, density=True, range=[-3.0,3.0], label=r'S$_1$$\rightarrow$S$_2$ hops')
            ax[0,1].hist(S10[0], bins=40, histtype='step', color='blue', alpha=0.9, density=True, range=[-3.0,3.0], label=r'S$_1$$\rightarrow$S$_0$ hops')

            ax[1,1].hist(S21[1], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[320.0,400.0], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[1,1].hist(S12[1], bins=40, histtype='step', color='limegreen', alpha=0.9, linewidth=1.5, density=True, range=[320.0,400.0], label=r'S$_1$$\rightarrow$S$_2$ hops')
            ax[1,1].hist(S10[1], bins=40, histtype='step', color='blue', alpha=0.9, density=True, range=[320.0,400.0], label=r'S$_1$$\rightarrow$S$_0$ hops')

            ax[2,1].hist(S21[2], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-90.0,90.0], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[2,1].hist(S12[2], bins=40, histtype='step', color='limegreen', alpha=0.9, linewidth=1.5, density=True, range=[-90.0,90.0], label=r'S$_1$$\rightarrow$S$_2$ hops')
            ax[2,1].hist(S10[2], bins=40, histtype='step', color='blue', alpha=0.9, density=True, range=[-90.0,90.0], label=r'S$_1$$\rightarrow$S$_0$ hops')
 
            os.chdir(run_path)

        
        #plt.legend(frameon=False, loc='upper right')
        ax[0,0].set_ylim(0,2)
        ax[1,0].set_ylim(0,0.15)
        ax[2,0].set_ylim(0,0.08)
        ax[0,1].set_ylim(0,2.2)
        ax[1,1].set_ylim(0,0.15)
        ax[2,1].set_ylim(0,0.08)
        ax[0,0].set_xlim(-3.0,3.0)
        ax[1,0].set_xlim(320,380)
        ax[2,0].set_xlim(-90.0,90.0)
        ax[0,1].set_xlim(-3.0,3.0)
        ax[1,1].set_xlim(320,380)
        ax[2,1].set_xlim(-90.0,90.0)
        ax[0,0].set_xticks([-3,-2,-1,0,1,2,3])
        ax[0,1].set_xticks([-3,-2,-1,0,1,2,3])
        ax[2,0].set_xticks([-90,-60,-30,0,30,60,90])
        ax[2,1].set_xticks([-90,-60,-30,0,30,60,90])
        #ax[1,1].set_xticks([-1.2,-0.6,0.0,0.6,1.2])
        #ax[1,1].set_yticks([0.0,0.4,0.8,1.2,1.6])
        #ax[0,2].set_xticks([340,350,360,370,380])
        #ax[0,2].set_yticks([0.0,0.05,0.1,0.15])
        #ax[1,2].set_xticks([340,350,360,370,380])
        #ax[1,2].set_yticks([0.0,0.05,0.1,0.15])
        #ax[0,0].set_xticks([-0.4,-0.2,0.0,0.2,0.4])
        #ax[1,0].set_xticks([-0.4,-0.2,0.0,0.2,0.4])
        #ax[0,0].set_yticks([0,1,2,3,4,5,6])
        #ax[1,0].set_yticks([0,1,2,3,4,5,6])
        #plt.xlim(-1.5,1.5)
        #plt.xticks(fontsize=14)
        #plt.yticks(fontsize=14)
        ax[0,0].set_ylabel('Norm. dist.', fontsize=14)
        ax[1,0].set_ylabel('Norm. dist.', fontsize=14)
        ax[2,0].set_ylabel('Norm. dist.', fontsize=14)
        ax[0,0].set_xlabel(r'HT ($\AA$)', fontsize=14)
        ax[1,0].set_xlabel(r'SOA ($^{\circ}$)', fontsize=14)
        ax[2,0].set_xlabel(r'PyrC ($^{\circ}$)', fontsize=14)
        ax[0,1].set_xlabel(r'HT ($\AA$)', fontsize=14)
        ax[1,1].set_xlabel(r'SOA ($^{\circ}$)', fontsize=14)
        ax[2,1].set_xlabel(r'PyrC ($^{\circ}$)', fontsize=14)

        #ax[1,0].set_xlabel(r'PyrC ($^{\circ}$)', fontsize=14)
        #plt.xlabel(r"H-transfer Coordinate ($\AA$)", fontsize=14)
        ax[0,0].minorticks_on()
        ax[1,0].minorticks_on()
        ax[2,0].minorticks_on()
        ax[0,1].minorticks_on()
        ax[1,1].minorticks_on()
        ax[2,1].minorticks_on()
        ax[0,1].tick_params(axis='x', labelsize=14)
        ax[0,1].tick_params(axis='y', labelsize=14)
        ax[1,1].tick_params(axis='x', labelsize=14)
        ax[1,1].tick_params(axis='y', labelsize=14)
        ax[0,0].tick_params(axis='x', labelsize=14)
        ax[0,0].tick_params(axis='y', labelsize=14)
        ax[1,0].tick_params(axis='x', labelsize=14)
        ax[1,0].tick_params(axis='y', labelsize=14)
        ax[2,0].tick_params(axis='x', labelsize=14)
        ax[2,0].tick_params(axis='y', labelsize=14)
        ax[2,1].tick_params(axis='x', labelsize=14)
        ax[2,1].tick_params(axis='y', labelsize=14)

        plt.subplots_adjust(wspace=0.15,hspace=0.0)
        ax[2,1].legend(fontsize=14, frameon=False, loc='upper right')
        fig.tight_layout()
        fig.savefig(path_script + 'differentiate-hops-updated.pdf',dpi=500)
