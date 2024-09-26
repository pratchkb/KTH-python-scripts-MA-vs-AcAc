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

#def calculate_s1s0hops(filename):
#
#    import numpy as np
#    # load proper modules
#    #import get_optim
#    #import geom_param
#    #from itertools import islice
#    import matplotlib.pyplot as plt
#    #import os
#    import sys
#    sys.path.append('/Users/pratip/')
#    from aimsrelated.tcutil.code.utils import units
#    from aimsrelated.tcutil.code.atom_data import atom_data
#    from aimsrelated.tcutil.code.get_optim import get_optim
#    from aimsrelated.tcutil.code.geom_param import geom_param
#    from itertools import islice
#
#    # break up the trajectory in individual xyz files
#    fdata = open(filename, 'r')   # use your data file
#
#    num = 15 # number of atoms for this particular system # please change for your system
#
#    before = []
#    traj = []
#
#    #h1 = []
#    #h2 = []
#
#    for line in fdata:
#        before.append(line)
#        if len(before) > 2:
#            before.pop(0)
#        if "t=  " in line:
#             t = line.split()
#             t = t[0][-5:]
#             traj.append(t)           # timestep appended
#
#    return traj

def collect_geom_param_hops(filename, atom_order, traj, natoms):

    import numpy as np
    import os
    import os.path
    import glob
    # load proper modules
    import get_optim
    import geom_param

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    bla4 = 0

    #print(traj)
    os.system('mkdir tmp')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
            t = line.split()
            if traj == t[0][-10:]:        # checking traj number
                t = t[2]
                #tim.append(float(t)) #* au_to_fs    # timestep converted to fs 
                dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded
                bla4 = calculate_bla4(data,atom_order)[0]


    fdata.close()
    os.system('rm -r tmp')

    return bla4


def collect_geom_param(filename, natoms, h1, h2, central_angles):

    import get_optim
    import geom_param
    #order1 = calculate_order_atoms_torsion(filename, natoms)
    #print(order1)
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    #timstp = []
    Htran = []
    angle = []

    os.system('mkdir tmp')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
             t = line.split()
             t = t[2]
             l = float(t) #* au_to_fs    # timestep converted to fs
             #timstp.append(l)           # timestep appended
             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

             Htran.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))
             
             ang = 0
             for idx in central_angles:
                 ang += geom_param.compute_angle(data,*idx)
             angle.append(ang)

        if "Time step:" in line:
             t = line.split()
             t = t[1]
             #l = float(t) #* au_to_fs    # timestep converted to fs
             #timstp.append(t)           # timestep appended
             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

             ang = 0
             #compute sum of angles
             for idx in central_angles:
                 ang += geom_param.compute_angle(data,*idx)
             angle.append(ang)

             Htran.append((geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1)))
             #if h1 < h2:
                 #bla4.append(calculate_distance(data[1],data[0]) - calculate_distance(data[0],data[2]) + calculate_distance(data[2],data[3]) - calculate_distance(data[3],data[4]))
             #else:
                 #bla4.append(calculate_distance(data[4],data[3]) - calculate_distance(data[3],data[2]) + calculate_distance(data[2],data[0]) - calculate_distance(data[1],data[0]))



    fdata.close()

    os.system('rm -r tmp')

    return Htran, angle


if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    #plots
    fig, ax = plt.subplots(2,3, figsize=(11,7))
   
    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/TRANFORMATION-FC-S21HTI/'

    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        
        #plots
        if molecule == "MA":
            hh1 = [0,8]; hh2 = [1,8]
            central_angles = [[0,2,3], [2,3,4], [3,4,1]]
            os.chdir(run_path)
            FC = collect_geom_param(run_path + 'total-FC.xyz', numatoms, hh1, hh2, central_angles)
            #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms)
            S21 = collect_geom_param(run_path + 'crossing-s2-s1.xyz', numatoms, hh1, hh2, central_angles)


            ax[0,1].hist(FC[0], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-1.5,1.5], label='ICs')
            ax[0,1].hist(S21[0], bins=40, histtype='step', color='blue', alpha=0.9, linewidth=1.5, density=True, range=[-1.5,1.5], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[0,2].hist(FC[1], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[340,380], label='ICs')
            ax[0,2].hist(S21[1], bins=40, histtype='step', color='blue', alpha=0.9, linewidth=1.5, density=True, range=[340,380], label=r'S$_2$$\rightarrow$S$_1$ hops')

            bla4IC = []
            blas21 = []
            for i in range(1,ntraj):
                i = "{:05d}".format(i)
                os.chdir('TRAJ_'+str(i))
                if os.path.exists('DONT_ANALYZE'):
                    os.chdir(run_path)
                else:

                    atom_order = calculate_order_atoms("output.xyz", hh1, hh2, molecule, numatoms)
                    bla4IC.append(calculate_IC_BLA("output.xyz", atom_order, numatoms))
                    #print(collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz', atom_order,'TRAJ_'+str(i), numatoms))                   
                    blas21.append(collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz', atom_order,'TRAJ_'+str(i), numatoms))                   

                os.chdir(run_path)


            ax[0,0].hist(bla4IC, bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-0.5,0.5], label='ICs')
            ax[0,0].hist(blas21, bins=40, histtype='step', color='blue', alpha=0.9, linewidth=1.5, density=True, range=[-0.5,0.5], label='S$_2$$\rightarrow$S$_1$ hops')
            

        if molecule == "AcAc":
            hh1 = [1,13]; hh2 = [4,13]
            central_angles = [[1,0,2], [0,2,3], [2,3,4]]
            os.chdir(run_path)
            FC = collect_geom_param(run_path + 'total-FC.xyz', numatoms, hh1, hh2, central_angles)
            #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms)
            S21 = collect_geom_param(run_path + 'crossing-s2-s1.xyz', numatoms, hh1, hh2, central_angles)
            ax[1,1].hist(FC[0], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-1.5,1.5], label='ICs')
            ax[1,1].hist(S21[0], bins=40, histtype='step', color='blue', alpha=0.9, linewidth=1.5, density=True, range=[-1.5,1.5], label=r'S$_2$$\rightarrow$S$_1$ hops')
            ax[1,2].hist(FC[1], bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[340,380], label='ICs')
            ax[1,2].hist(S21[1], bins=40, histtype='step', color='blue', alpha=0.9, linewidth=1.5, density=True, range=[340,380], label=r'S$_2$$\rightarrow$S$_1$ hops')

            bla4IC = []
            blas21 = []
            for i in range(1,ntraj):
                i = "{:05d}".format(i)
                os.chdir('TRAJ_'+str(i))
                if os.path.exists('DONT_ANALYZE'):
                    os.chdir(run_path)
                else:
                    atom_order = calculate_order_atoms("output.xyz", hh1, hh2, molecule, numatoms)
                    bla4IC.append(calculate_IC_BLA("output.xyz", atom_order, numatoms))
                    blas21.append(collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz', atom_order,'TRAJ_'+str(i), numatoms))                   
 
                os.chdir(run_path)


            ax[1,0].hist(bla4IC, bins=40, histtype='bar', color='orange', alpha=0.9, density=True, range=[-0.5,0.5], label='ICs')
            ax[1,0].hist(blas21, bins=40, histtype='step', color='blue', alpha=0.9, linewidth=1.5, density=True, range=[-0.5,0.5], label='S$_2$$\rightarrow$S$_1$ hops')
       

        ax[0,1].axvline(x=-0.661, color='black', linestyle='dashed', label='FC')
        ax[1,1].axvline(x=-0.591, color='black', linestyle='dashed', label='FC')
        ax[0,1].axvline(x=0.661, color='black', linestyle='dashed', label='FC')
        ax[1,1].axvline(x=0.591, color='black', linestyle='dashed', label='FC')
        ax[0,2].axvline(x=366.5, color='black', linestyle='dashed', label='FC')
        ax[1,2].axvline(x=363.9, color='black', linestyle='dashed', label='FC')
        ax[0,0].axvline(x=-0.151, color='black', linestyle='dashed', label='FC')
        ax[1,0].axvline(x=-0.150, color='black', linestyle='dashed', label='FC')
        #plt.axvline(x=-0.416, color='gray', linestyle=':', label='Avg. Asym-HTI-MECIs')
        #plt.axvline(x=0.203, color='green', linestyle=':', label='Asym Ecdown')
        #plt.axvline(x=0.226, color='darkred', linestyle=':', label='Asym Ecup')
        #plt.axvline(x=0.001, color='magenta', linestyle='dashdot', label='Avg. Sym-HTI-MECIs')
        #plt.axvline(x=357.1, color='gray', linestyle='dashdot', label='Asym HTI down')
        
        #plt.legend(frameon=False, loc='upper right')
        ax[0,1].set_ylim(0,1.8)
        ax[1,1].set_ylim(0,1.8)
        ax[0,2].set_ylim(0,0.16)
        ax[1,2].set_ylim(0,0.16)
        ax[0,0].set_ylim(0,6.5)
        ax[1,0].set_ylim(0,6.5)
        ax[0,1].set_xlim(-1.5,1.5)
        ax[1,1].set_xlim(-1.5,1.5)
        ax[0,2].set_xlim(340,380)
        ax[1,2].set_xlim(340,380)
        ax[0,0].set_xlim(-0.5,0.5)
        ax[1,0].set_xlim(-0.5,0.5)
        ax[0,1].set_xticks([-1.2,-0.6,0.0,0.6,1.2])
        ax[0,1].set_yticks([0.0,0.4,0.8,1.2,1.6])
        ax[1,1].set_xticks([-1.2,-0.6,0.0,0.6,1.2])
        ax[1,1].set_yticks([0.0,0.4,0.8,1.2,1.6])
        ax[0,2].set_xticks([340,350,360,370,380])
        ax[0,2].set_yticks([0.0,0.05,0.1,0.15])
        ax[1,2].set_xticks([340,350,360,370,380])
        ax[1,2].set_yticks([0.0,0.05,0.1,0.15])
        ax[0,0].set_xticks([-0.4,-0.2,0.0,0.2,0.4])
        ax[1,0].set_xticks([-0.4,-0.2,0.0,0.2,0.4])
        ax[0,0].set_yticks([0,1,2,3,4,5,6])
        ax[1,0].set_yticks([0,1,2,3,4,5,6])
        #plt.xlim(-1.5,1.5)
        #plt.xticks(fontsize=14)
        #plt.yticks(fontsize=14)
        ax[0,0].set_ylabel('Norm. dist.', fontsize=16)
        ax[1,0].set_ylabel('Norm. dist.', fontsize=16)
        ax[1,1].set_xlabel(r'HT ($\AA$)', fontsize=16)
        ax[1,2].set_xlabel(r'SOA ($^{\circ}$)', fontsize=16)
        ax[1,0].set_xlabel(r'BLA ($\AA$)', fontsize=16)
        #plt.xlabel(r"H-transfer Coordinate ($\AA$)", fontsize=14)
        ax[0,1].minorticks_on()
        ax[1,1].minorticks_on()
        ax[0,2].minorticks_on()
        ax[1,2].minorticks_on()
        ax[0,0].minorticks_on()
        ax[1,0].minorticks_on()
        ax[0,1].tick_params(axis='x', labelsize=14)
        ax[0,1].tick_params(axis='y', labelsize=14)
        ax[1,1].tick_params(axis='x', labelsize=14)
        ax[1,1].tick_params(axis='y', labelsize=14)
        ax[0,2].tick_params(axis='x', labelsize=14)
        ax[0,2].tick_params(axis='y', labelsize=14)
        ax[1,2].tick_params(axis='x', labelsize=14)
        ax[1,2].tick_params(axis='y', labelsize=14)
        ax[0,0].tick_params(axis='x', labelsize=14)
        ax[0,0].tick_params(axis='y', labelsize=14)
        ax[1,0].tick_params(axis='x', labelsize=14)
        ax[1,0].tick_params(axis='y', labelsize=14)

        plt.subplots_adjust(wspace=0.15,hspace=0.0)
        fig.tight_layout()
        fig.savefig(path_script + 'crossing-Htransfer-updated.pdf',dpi=500)
