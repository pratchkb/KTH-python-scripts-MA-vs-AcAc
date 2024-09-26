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

def spectrum(HT,weight,sigma,y):
    gE=[]
    for Ei in y:
        tot=0
        for ht,wt in zip(HT,weight):
            tot+=wt*np.exp(-((((ht-Ei)/sigma)**2)))
        gE.append(tot)
    return gE

def calculate_distance(atom1_coord, atom2_coord):
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12


def collect_geom_param(filename, natoms, hh1, hh2):

    num = natoms # number of atoms for this particular system # please change for your system

    time = []
    htran = []

    data = get_optim.read_xyzs_as_np(filename)
    time = []
    htran = []
    for i in range(len(data[:])):
        time.append(i*0.5)
        h1 = geom_param.compute_bond(data[i],hh2)
        h2 = geom_param.compute_bond(data[i],hh1)
        htran.append(h1-h2)

    return time, htran, data


def calculate_param_statewise(filename, natoms):

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


def collect_geom_param_hops(filename, natoms, hh1, hh2):

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstps2 = []
    Htrans2 = []

    os.system('mkdir tmp3')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
            t = line.split()
            ti = t[2]
            l = float(ti) #* au_to_fs    # timestep converted to fs
            timstps2.append(l)           # timestep appended
            dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
            f=open(dat,'w')
            f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
            f.close()
            data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

            Htrans2.append(geom_param.compute_bond(data,hh2)-geom_param.compute_bond(data,hh1))

    return timstps2, Htrans2



def calculate_param_statewise_color(filename, natoms, hh1, hh2):


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
    Nused=[200,203]
    #Ntraj = [10,10]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/HT-with-backhop/'
    for molecule, ntraj, nused, run_path, numatoms  in zip(molecules, Ntraj, Nused, runpaths, Natoms):

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

        fig,ax = plt.subplots(figsize=(4.5,3))
        clr = ['#0072BD', 'orange', '#66BF2E', 'k', 'limegreen', 'darkgreen', 'darkblue', 'orangered', '#D45500', 'crimson', 'lime', 'cyan']


        os.chdir(run_path)
        #IC = collect_geom_param(run_path + 'total-FC.xyz', numatoms, hh1, hh2)  #ICs
        #for21 = collect_geom_param(run_path + 'crossing-s2-s1.xyz', numatoms, hh1, hh2)
        #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms, hh1, hh2)
        #back21 = collect_geom_param(run_path + 'crossing-s1-s2.xyz', numatoms, hh1, hh2)


        t21 = []
        s21 = []
        
        t12 = []
        s12 = []


        f=open(path_script + f'{molecule}-ICs-long.dat','w')
        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir(run_path + 'TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):
                os.chdir(run_path)
            else:
            
                #atom_order = calculate_order_atoms("output.xyz", molecule, numatoms)
            	#densities = calculate_param_statewise_color("output.xyz", atom_order, 15)
                densities = calculate_param_statewise_color("output.xyz", numatoms, hh1, hh2)
                pos = np.where(np.abs(np.diff(densities[0])) > 0.5)[0] + 1
                time = np.insert(densities[0], pos, np.nan)
                ht = np.insert(densities[1], pos, np.nan)
                ax.plot(time, ht, color = clr[1], alpha=0.9, linewidth=0.4, zorder=0)
               
                if len(time) == 0:
                    continue 
                if time[-1] > 25:               
                    f.write('{}\t{}\n'.format(molecule,i))
                    print(molecule,i)
                    
                os.chdir(run_path)
        
        f.close() 

