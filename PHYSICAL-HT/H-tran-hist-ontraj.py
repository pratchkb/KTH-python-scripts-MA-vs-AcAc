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

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes


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


def collect_geom_param(filename, natoms):


    import get_optim
    import geom_param
    #order1 = calculate_order_atoms_torsion(filename, natoms)
    #print(order1)
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system


    before = []
    timstp = []
    Htran = []

    os.system('mkdir tmp')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
             t = line.split()
             t = t[2]
             l = float(t) #* au_to_fs    # timestep converted to fs
             timstp.append(l)           # timestep appended
             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded
             
             Htran.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

    fdata.close() 

    os.system('rm -r tmp')

    return timstp, Htran 

#def calculate_order_atoms(filename, molecule, natoms):
#    """
#    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
#    """
#
#    import get_optim
#    import geom_param
#    # break up the trajectory in individual xyz files
#    fdata = open(filename, 'r')   # use your data file
#
#    num = natoms # number of atoms for this particular system # please change for your system
#
#    before = []
#    timstp = []
#
#    order = []
#
#    os.system('mkdir tmp')
#    for line in fdata:
#        before.append(line)
#        if len(before) > 2:
#            before.pop(0)
#        if "t=         0.00000" in line: # Initial condition geometry
#             t = line.split()
#             t = t[1]
#             l = float(t) #* au_to_fs    # timestep converted to fs
#             timstp.append(l)           # timestep appended
#             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
#             f=open(dat,'w')
#             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
#             f.close()
#             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded
#
#             h1 = geom_param.compute_bond(data,hh1)
#             h2 = geom_param.compute_bond(data,hh2)
#
#             if h1 > h2:
#                 if molecule == "AcAc":
#                     order = [1,0,2,3,4]
#                 else:
#                     order = [0,2,3,4,1]
#             else:
#                 if molecule == "AcAc":
#                     order = [4,3,2,0,1]
#                 else:
#                     order = [1,4,3,2,0]
#
#    fdata.close()
#
#    os.system('rm -r tmp')
#
#    return order

def calculate_s1s0hops(filename, natoms):

    import numpy as np
    # load proper modules
    import get_optim
    import geom_param
    #from itertools import islice
    from itertools import islice

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    traj = []

    #h1 = []
    #h2 = []

    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
             t = line.split()
             t = t[0][-5:]
             traj.append(t)           # timestep appended

    return traj

def calculate_param_statewise(filename, natoms):

    import get_optim
    import geom_param
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstps2 = []
    timstps1 = []
    timstps0 = []

    Htrans2 = []
    Htrans1 = []
    Htrans0 = []

    os.system('mkdir tmp3')
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
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

                Htrans2.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

            if int(st) == 2 or int(st) == 1:         # S1 or S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps1.append(l)           # timestep appended
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

                Htrans1.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

    os.system('rm -r tmp3')
    #os.chdir(path)

    return timstps2, Htrans2, timstps1, Htrans1


def calculate_param_nostate(filename, natoms):

    import get_optim
    import geom_param
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstp = []
    Htran = []

    os.system('mkdir tmp3')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
            t = line.split()
            ti = t[1]
            st = t[2]
            l = float(ti) #* au_to_fs    # timestep converted to fs
            timstp.append(l)           # timestep appended
            dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
            f=open(dat,'w')
            f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
            f.close()
            data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

            Htran.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

    os.system('rm -r tmp3')
    #os.chdir(path)

    return timstp, Htran


def calculate_param_statewise_color(filename, atom_order, natoms):


    import get_optim
    import geom_param
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstps2 = []
    timstps1 = []
    timstps0 = []
    Htrans2 = []
    Htrans1 = []
    Htrans0 = []

    os.system('mkdir tmp3')
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
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

                Htrans2.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

            if int(st) == 2:         # S1 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps1.append(l)           # timestep appended
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

                Htrans1.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

            if int(st) == 1:         # S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps0.append(l)           # timestep appended
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w') 
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded
                 
                Htrans0.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

    os.system('rm -r tmp3')
    #os.chdir(path)

    return timstps2, Htrans2, timstps1, Htrans1, timstps0, Htrans0


if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/PHYSICAL-HT/'
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
        gs = GridSpec(1,2, width_ratios=[3, 1])
        clr = ['#0072BD', 'orange', '#66BF2E', 'k', 'green', 'darkgreen', 'darkblue', 'orangered', '#D45500']
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])


        os.chdir(run_path)
        for21 = collect_geom_param(run_path + 'crossing-s2-s1.xyz', numatoms)
        #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms)
        back21 = collect_geom_param(run_path + 'crossing-s1-s2.xyz', numatoms)

        HT1 = []
        HT2 = []

        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir('TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):
                os.chdir(run_path)
            else:
            
                #atom_order = calculate_order_atoms("output.xyz", molecule, numatoms)
            	#densities = calculate_param_statewise_color("output.xyz", atom_order, 15)
                densities = calculate_param_nostate("output.xyz", numatoms)
                pos = np.where(np.abs(np.diff(densities[0])) > 0.5)[0] + 1
                time = np.insert(densities[0], pos, np.nan)
                ht = np.insert(densities[1], pos, np.nan)
                ax1.plot(time, ht, color = clr[1], alpha=0.9, linewidth=0.4, zorder=0)
  
                densit = calculate_param_statewise("output.xyz", numatoms)
                posht1 = np.where(np.abs(np.diff(densit[0])) > 0.5)[0] + 1
                timeht1 = np.insert(densit[0], posht1, np.nan)
                htht1 = np.insert(densit[1], posht1, np.nan)
                for elem in htht1:
                    HT1.append(elem) 

           	    #poss1 = np.where(np.abs(np.diff(densities[2])) > 0.5)[0] + 1
            	#times1 = np.insert(densities[2], poss1, np.nan)
            	#hts1 = np.insert(densities[3], poss1, np.nan)
            	#poss0 = np.where(np.abs(np.diff(densities[4])) > 0.5)[0] + 1
            	#times0 = np.insert(densities[4], poss0, np.nan)
            	#hts0 = np.insert(densities[5], poss0, np.nan) 

                posht2 = np.where(np.abs(np.diff(densit[2])) > 0.5)[0] + 1
                timeht2 = np.insert(densit[2], posht2, np.nan)
                htht2 = np.insert(densit[3], posht2, np.nan)
                for elem in htht2[:10]:        # 5 fs # first 10 elements
                    HT2.append(elem)

                #ax1.plot(time, ht, color = clr[1], alpha=0.9, linewidth=0.4, zorder=2)
            	#ax1.plot(times2, hts2, color = clr[0], alpha=0.9, linewidth=0.4, zorder=2)
            	#ax1.plot(times1, hts1, color = clr[1], alpha=0.9, linewidth=0.4, zorder=0)
            	#ax1.plot(times0, hts0, color = clr[4], linewidth=0.4, zorder=1)

                os.chdir(run_path)    


        ax2.hist(HT1, bins=40, range=[-3,3], density=True, weights=None, orientation='horizontal', color = 'orange', alpha = 0.9, zorder=10)
        ax2.hist(HT2, bins=40, range=[-3,3], density=True, histtype='step', orientation='horizontal', weights=None, color = 'blue', alpha = 0.9, zorder=10)

        ax1.scatter(for21[0], for21[1], color = clr[8], s=0.6, label=r'S$_2$/S$_1$ hops', zorder=60)
    	#ax1.scatter(for10[0], for10[1], s=0.6, color = clr[5], label=r'S$_1$/S$_0$ hops', zorder=61)
        ax1.scatter(back21[0], back21[1], s=0.6, color = clr[6], label=r'S$_1$/S$_2$ hops', zorder=62)
        if molecule == "MA":
            ax1.axhline(y=-0.661, color='k', linestyle='dashed', lw = 0.7, label='FC', zorder=48)
            ax1.axhline(y=0.661, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=49)
            ax2.axhline(y=-0.661, color='k', linestyle='dashed', lw = 0.7, label='FC',zorder=50)
            ax2.axhline(y=0.661, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=51)
        if molecule == "AcAc":
            ax1.axhline(y=-0.591, color='k', linestyle='dashed', lw = 0.7, label='FC', zorder=48)
            ax1.axhline(y=0.591, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=49)
            ax2.axhline(y=-0.591, color='k', linestyle='dashed', lw = 0.7, label='FC',zorder=50)
            ax2.axhline(y=0.591, color='k', linestyle='dashed', lw= 0.7, label='FC',zorder=51)

        ax1.set_ylabel(r"HT ($\AA$)")
        ax1.set_xlabel(r"Time (fs)")
        ax2.set_xlabel(r"Norm. dist.")
        ax1.set_xlim([0,200])
        ax1.set_ylim([-3.0,3.0])
        ax2.set_ylim([-3.0,3.0])
        ax2.set_xlim([0.0,1.5])
        ax2.set_xticks([0.0,0.5,1.0,1.5])
        ax1.minorticks_on()
        ax2.minorticks_on()
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.15)
        plt.savefig(path_script + f'{molecule}-H-transfer-density-colored.png', dpi=500)
        plt.savefig(path_script + f'{molecule}-H-transfer-density-colored.pdf', dpi=500)


