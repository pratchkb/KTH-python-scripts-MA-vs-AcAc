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

def calculate_order_atoms_torsion(filename, s1last, molecule, natoms, hh1, hh2):

    num = natoms # number of atoms for this particular system # please change for your system

    data = get_optim.read_xyzs_as_np(filename)
    frame = int(float(s1last)*2)
        
    h1 = geom_param.compute_bond(data[frame],hh2)
    h2 = geom_param.compute_bond(data[frame],hh1)
    if molecule == "MA":
        if h1 < h2:
            order = [2,3,4,1]
        else:
            order = [4,3,2,0]
    if molecule == "AcAc":
        if h1 < h2:
            order = [0,2,3,4]
        else:
            order = [3,2,0,1]

    return order




#def calculate_order_atoms_torsion(filename,s1last,molecule,natoms,hh1,hh2):
#    """
#    calculate the numbering of atoms depending on which O the H is attahced to, in the final geometry to back track 
#    """
#    #import get_optim
#    #import geom_param
#    #from itertools import islice
#    
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
#    #print(s1last)
#    os.system('mkdir tmp')
#    for line in fdata:
#        before.append(line)
#        if len(before) > 2:
#            before.pop(0)
#        if " t=       " in line: # Initial condition geometry
#            t = line.split()
#            if t[1] == s1last:
#                #print(t[1])
#                l = float(t[1]) #* au_to_fs    # timestep converted to fs
#                timstp.append(l)           # timestep appended
#                dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
#                f=open(dat,'w')
#                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
#                f.close()
#                data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded
#            
#                h1 = geom_param.compute_bond(data,hh2)
#                h2 = geom_param.compute_bond(data,hh1)
#                if molecule == "MA":
#                    if h1 < h2:
#                        order = [2,3,4,1]
#                    else:
#                        order = [4,3,2,0]
#                if molecule == "AcAc":
#                    if h1 < h2:
#                        order = [0,2,3,4]
#                    else:
#                        order = [3,2,0,1]
#
#
#    fdata.close()
#
#    os.system('rm -r tmp')
#
#    return order


def calculate_sort(filename, natoms):

    import numpy as np
    # load proper modules
    #import get_optim
    #import geom_param
    #from itertools import islice
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

def collect_geom_param_hops(filename, trajname, atom_order, traj, natoms, pyrs):

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    data = get_optim.read_xyzs_as_np(trajname)

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    tim = []
    tor = []
    pyr = []

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
 
                tor.append(geom_param.compute_torsion(data[frame],orders[0],orders[1],orders[2],orders[3]))
                pyr.append(geom_param.compute_pyramidalization(data[frame],*pyrs))


    fdata.close()

    return tim, tor, pyr


#def collect_geom_param_hops(filename, atom_order, traj, natoms, pyrs):
#
#    # load proper modules
#    #import get_optim
#    #import geom_param
#
#    # break up the trajectory in individual xyz files
#    fdata = open(filename, 'r')   # use your data file
#
#    num = natoms # number of atoms for this particular system # please change for your system
#
#    before = []
#    tim = []
#    tor = []
#    pyr = []
#
#    #print(traj)
#    os.system('mkdir tmp')
#    for line in fdata:
#        before.append(line)
#        if len(before) > 2:
#            before.pop(0)
#        if "t=  " in line:
#            t = line.split()
#            if traj == t[0][-10:]:        # checking traj number
#                t = t[2]
#                tim.append(float(t)) #* au_to_fs    # timestep converted to fs 
#                dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
#                f=open(dat,'w')
#                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
#                f.close()
#                data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded
#                tor.append(geom_param.compute_torsion(data,orders[0],orders[1],orders[2],orders[3]))
#                pyr.append(geom_param.compute_pyramidalization(data,*pyrs))
#
#
#    fdata.close()
#    os.system('rm -r tmp')
#
#    return tim, tor, pyr
#


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

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/EARLY-LATE-PLOTS-SCATTER/'

    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        if molecule == "MA":
            hh1 = [0,8]; hh2 = [1,8]
            pyrs = [3,2,4,6]

        if molecule == "AcAc":
            hh1 = [1,13]; hh2 = [4,13]
            pyrs = [2,0,3,14]

        os.chdir(run_path)
        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir('TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):

                #print(str(i))
                os.chdir(run_path)

            else:

                # calling functions
                timeall = calculate_sort("output.xyz", numatoms)
                #s2 = timeall[0][-1]
                #s2 = f"{s2}{'0'*4}"
                s1 = timeall[1][-1]
                #print(s1) 
                s1 = f"{s1}{'0'*4}"
                orders = calculate_order_atoms_torsion("output.xyz",str(s1),molecule,numatoms,hh1,hh2)
                #print(orders)
                for21 = collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz',run_path + 'TRAJ_'+str(i) + '/output.xyz', orders,'TRAJ_'+str(i),numatoms,pyrs)
                #print(for21)
                for10 = collect_geom_param_hops(run_path + 'crossing-s1-s0.xyz',run_path + 'TRAJ_'+str(i) + '/output.xyz',orders,'TRAJ_'+str(i),numatoms,pyrs)
                #print(for10)
                back21 = collect_geom_param_hops(run_path + 'crossing-s1-s2.xyz', run_path + 'TRAJ_'+str(i) + '/output.xyz', orders,'TRAJ_'+str(i),numatoms,pyrs)

                tim21 = []
                tor21 = []
                pyr21 = []
                for elem,torelem,pyrelem in zip(for21[0],for21[1],for21[2]):
                    tim21.append(elem)
                    tor21.append(abs(torelem))
                    pyr21.append(abs(pyrelem))
            
                timless = []
                torless = []
                pyrless = []
                timbig = []
                torbig = []
                pyrbig = []
             
                threshold = 75
                #thfor10 = threshold_hop(for10[0],for10[1],for10[2],75)
                for elem,torelem,pyrelem in zip(for10[0],for10[1],for10[2]):
                    if elem < threshold:
                        timless.append(elem)
                        torless.append(abs(torelem))
                        pyrless.append(abs(pyrelem))
                    else:
                        timbig.append(elem)
                        torbig.append(abs(torelem))
                        pyrbig.append(abs(pyrelem))

                if molecule == "MA":
                    if for21[0] != 0:
                        ax[0].scatter(tor21, pyr21, s=15.0, c='#D45500', zorder=5)
                    if for10[0] != 0:
                        ax[0].scatter(torless,pyrless, s=15.0, c='darkgreen', zorder=15)
                        ax[0].scatter(torbig,pyrbig, s=15.0, c='limegreen', zorder=16)

                if molecule == "AcAc":
                    if for21[0] != 0:
                        ax[1].scatter(tor21, pyr21, s=15.0, c='#D45500', zorder=5)
                    if for10[0] != 0:
                        ax[1].scatter(torless,pyrless, s=15.0, c='darkgreen', zorder=15)
                        ax[1].scatter(torbig,pyrbig, s=15.0, c='limegreen', zorder=16)
                 
                os.chdir(run_path)


    ax[0].set_ylabel(r"|PyrC| ($^{\circ}$)")
    ax[1].set_ylabel(r"|PyrC| ($^{\circ}$)")
    #ax[0].set_xlabel(r"Torsion ($^{\circ}$)")
    ax[1].set_xlabel(r"|Torsion| ($^{\circ}$)") 
    ax[0].scatter(96.1,0.8, marker='D', label=r'S$_1$/S$_0$-Tw', color='gold', edgecolors='k', zorder=20)
    ax[0].scatter(60.5,43.3, marker='P', s=60.0, label=r'S$_1$/S$_0$-CI$_\text{pyrC}$', color='gold',edgecolors='k', zorder=21)      
    ax[0].scatter(0,0, marker='*', s=60.0, label=r'S$_2$/S$_1$-HTI', color='gold',edgecolors='k', zorder=22)       
    ax[0].scatter(111.8,0.6, marker='s', label=r'S$_2$/S$_1$-Tw', color='gold',edgecolors='k', zorder=22)       
    ax[1].scatter(99.0,0.7, marker='D', label=r'S$_1$/S$_0$-Tw', color='gold', edgecolors='k', zorder=23)
    ax[1].scatter(55.6,45.8, marker='P', s=60.0, label=r'S$_1$/S$_0$-CI$_\text{pyrC}$', color='gold',edgecolors='k',zorder=24)
    ax[1].scatter(0,0, marker='*', s=60.0, label=r'S$_2$/S$_1$-HTI', color='gold',edgecolors='k', zorder=25)       
    ax[1].scatter(78.6,2.1, marker='s', label=r'S$_2$/S$_1$-Tw', color='gold',edgecolors='k', zorder=22)       
    ax[0].set_xlim(-5,120)
    ax[0].set_ylim(-5,90)
    ax[1].set_xlim(-5,120)
    ax[1].set_ylim(-5,90)
    ax[0].set_yticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0])
    ax[0].set_xticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0,105.0,120.0])
    ax[1].set_yticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0])
    ax[1].set_xticks([0.0,15.0,30.0,45.0,60.0,75.0,90.0,105.0,120.0])
    ax[0].minorticks_on()
    ax[1].minorticks_on()

    ax[0].legend(frameon=False, loc = 'upper left', ncols=2)
    ax[1].legend(frameon=False, loc = 'upper left', ncols=2)
    plt.tight_layout()
    fig.savefig(path_script + 'Torsionback-vs-pyr.png', dpi=500)
