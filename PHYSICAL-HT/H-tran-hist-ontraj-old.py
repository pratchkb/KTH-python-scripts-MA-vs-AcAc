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


#CACHE_FILENAME = "pop.pickle"
def calculate_distance(atom1_coord, atom2_coord):
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12

def calculate_order_atoms_torsion(filename, natoms):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the final geometry to back track 
    """


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
        if " t=       200.00000" in line: # Initial condition geometry
             t = line.split()
             t = t[1]
             #print(t)
             l = float(t) #* au_to_fs    # timestep converted to fs
             timstp.append(l)           # timestep appended
             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded
            

             h1 = geom_param.compute_bond(data,[4,13])
             h2 = geom_param.compute_bond(data,[1,13])

             if h1 < h2:
                 order = [0,2,3,4]
             else:
                 order = [3,2,0,1]

    fdata.close()

    os.system('rm -r tmp')

    return order


def calculate_bla4(data,atom_order):
    """     
    Calculates all 3 types of BLA from the xyz data file and the list of atom order
    """     
            
    BLA4 = calculate_distance(data[atom_order[0]],data[atom_order[1]]) - calculate_distance(data[atom_order[1]],data[atom_order[2]]) + calculate_distance(data[atom_order[2]],data[atom_order[3]]) - calculate_distance(data[atom_order[3]],data[atom_order[4]])
    BLACO = calculate_distance(data[atom_order[0]],data[atom_order[1]]) - calculate_distance(data[atom_order[3]],data[atom_order[4]])
    BLACC = - calculate_distance(data[atom_order[1]],data[atom_order[2]]) + calculate_distance(data[atom_order[2]],data[atom_order[3]])
            
    return BLA4, BLACO, BLACC


def collect_geom_param(filename, natoms):


    order1 = calculate_order_atoms_torsion(filename, natoms)
    #print(order1)
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system


    before = []
    timstp = []

    Hbond1 = []
    Hbond2 = []
    Htran = []
    pyr = []
    angles = []
    bla4 = []
    torsion = []
    torsion1 = []
    torsion2 = []
    
    tbelow20 = []
    htranbelow20 = []
    pyrbelow20 = []
    anglesbelow20 = []

    tabove20 = []
    htranabove20 = []
    pyrabove20 = []
    anglesabove20 = []

    tbelow70 = []
    htranbelow70 = []
    pyrbelow70 = []
    anglesbelow70 = []

    tabove70 = []
    htranabove70 = []
    pyrabove70 = []
    anglesabove70 = []

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

             #O-H distances
             Hbond1.append(geom_param.compute_bond(data,[13,4]))
             Hbond2.append(geom_param.compute_bond(data,[1,13]))
             #H-coord
             Htran.append((geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13])))
             #compute pyramidalization
             pyr.append(geom_param.compute_pyramidalization(data,2,0,3,14))
             #compute sum of angles
             angles.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))
             bla4.append(calculate_distance(data[1],data[0]) - calculate_distance(data[0],data[2]) + calculate_distance(data[2],data[3]) - calculate_distance(data[3],data[4]))
             #O-H distances
             h1 = geom_param.compute_bond(data,[4,13])
             h2 = geom_param.compute_bond(data,[1,13])

             #Htran.append((geom_param.compute_bond(data,[1,13])-geom_param.compute_bond(data,[13,4]))/2.0)
             if h1 < h2:
                 torsion.append(geom_param.compute_torsion(data,0,2,3,4))
             else:
                 torsion.append(geom_param.compute_torsion(data,3,2,0,1))

             torsion1.append(geom_param.compute_torsion(data,0,2,3,4))
             torsion2.append(geom_param.compute_torsion(data,3,2,0,1))

             #torback.append(geom_param.compute_torsion(data,order1))


             if l < 20:       #20 fs
                 tbelow20.append(l)
                 htranbelow20.append((geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13])))
                 pyrbelow20.append(geom_param.compute_pyramidalization(data,2,0,3,14))
                 anglesbelow20.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))
             else:            #>= 20 fs
                 tabove20.append(l)
                 htranabove20.append((geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13])))
                 pyrabove20.append(geom_param.compute_pyramidalization(data,2,0,3,14))
                 anglesabove20.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))

             if l < 65:       #65 fs
                 tbelow70.append(l)
                 htranbelow70.append((geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13])))
                 pyrbelow70.append(geom_param.compute_pyramidalization(data,2,0,3,14))
                 anglesbelow70.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))
             else:            #>= 65 fs
                 tabove70.append(l)
                 htranabove70.append((geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13])))
                 pyrabove70.append(geom_param.compute_pyramidalization(data,2,0,3,14))
                 anglesabove70.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))


    fdata.close() 

    os.system('rm -r tmp')

    return timstp, Hbond1, Hbond2, Htran, pyr, angles, tbelow20, htranbelow20, pyrbelow20, anglesbelow20, tabove20, htranabove20, pyrabove20, anglesabove20, tbelow70, htranbelow70, pyrbelow70, anglesbelow70, tabove70, htranabove70, pyrabove70, anglesabove70, bla4, torsion, torsion1, torsion2 

def calculate_order_atoms(filename, natoms):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

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

             h1 = geom_param.compute_bond(data,[1,13])
             h2 = geom_param.compute_bond(data,[4,13])

             if h1 > h2:
                 order = [1,0,2,3,4]
             else:
                 order = [4,3,2,0,1]

    fdata.close()

    os.system('rm -r tmp')

    return order

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

def calculate_param_statewise(filename, atom_order, natoms):


    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = natoms # number of atoms for this particular system # please change for your system

    before = []
    timstps2 = []
    timstps1 = []
    timstps0 = []



    bla4 = []
    blaCO = []
    blaCC = []
    bla4s1 = []
    blaCOs1 = []
    blaCCs1 = []
    bla4s0 = []
    blaCOs0 = []
    blaCCs0 = []
    Htrans2 = []
    Htrans1 = []
    Htrans0 = []
    pyrs0 = []
    pyrs1 = []
    pyrs2 = []
    angless0 = []
    angless1 = []
    angless2 = []

    pyrc3s0 = []
    pyrc3s1 = []
    pyrc3s2 = []
    pyrc5s0 = []
    pyrc5s1 = []
    pyrc5s2 = []

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

                s2 = calculate_bla4(data,atom_order)
                bla4.append(s2[0])
                blaCO.append(s2[1])
                blaCC.append(s2[2])

                Htrans2.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))
                #compute pyramidalization
                pyrs2.append(geom_param.compute_pyramidalization(data,2,0,3,14))
                pyrc3s2.append(geom_param.compute_pyramidalization(data,0,1,2,5))
                pyrc5s2.append(geom_param.compute_pyramidalization(data,3,2,4,9))
                #compute sum of angles
                angless2.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))

            if int(st) == 2 or int(st) == 1:         # S1 or S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps1.append(l)           # timestep appended
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

                s1 = calculate_bla4(data,atom_order)
                #bla4.append(s2[0])
                #blaCO.append(s2[1])
                #blaCC.append(s2[2])

                Htrans1.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))
                #compute pyramidalization
                pyrs1.append(geom_param.compute_pyramidalization(data,2,0,3,14))
                pyrc3s1.append(geom_param.compute_pyramidalization(data,0,1,2,5))
                pyrc5s1.append(geom_param.compute_pyramidalization(data,3,2,4,9))
                #compute sum of angles
                angless1.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))


    os.system('rm -r tmp3')
    #os.chdir(path)

    return timstps2, bla4, blaCO, blaCC, Htrans2, pyrs2, angless2, pyrc3s2, pyrc5s2, timstps1, Htrans1, pyrs1, angless1, pyrc3s1, pyrc5s1


def calculate_param_nostate(filename, atom_order, natoms):

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

            Htran.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))

    os.system('rm -r tmp3')
    #os.chdir(path)

    return timstp, Htran


def calculate_param_statewise_color(filename, atom_order, natoms):


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

                Htrans2.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))

            if int(st) == 2:         # S1 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps1.append(l)           # timestep appended
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded

                Htrans1.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))

            if int(st) == 1:         # S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps0.append(l)           # timestep appended
                dat = 'tmp3/' + str(ti) + '.xyz'       # coordinates saved
                f=open(dat,'w') 
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = get_optim.read_xyz_as_np('tmp3/' + str(ti) + '.xyz')        # coordinates loaded
                 
                Htrans0.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))

    os.system('rm -r tmp3')
    #os.chdir(path)

    return timstps2, Htrans2, timstps1, Htrans1, timstps0, Htrans0


if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    runpath = [
        "/data/projects/MAvsAcAc-project/AcAc-data/PROD_Final_Feb27/",
        "",
    ]

    #plots
    path_data1 = '/data/projects/MAvsAcAc-project/AcAc-data/PROD_Final_Feb27/'
    path_script = '/data/projects/MAvsAcAc-project/PLOTS/PHYSICAL-HT/'

    for21 = collect_geom_param(path_data1 + 'crossing-s2-s1.xyz', 15)
    for10 = collect_geom_param(path_data1 + 'crossing-s1-s0.xyz', 15)
    back21 = collect_geom_param(path_data1 + 'crossing-s1-s2.xyz', 15)

    #plot
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(6,3))
    gs = GridSpec(1,2, width_ratios=[3, 1])
    clr = ['#0072BD', 'orange', '#66BF2E', 'k', 'green', 'darkgreen', 'darkblue', 'orangered', '#D45500']
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])


    os.chdir(path_data1)

    N=5
    HT1 = []
    HT2 = []

    for i in range(1,N):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            os.chdir(path_data1)
        else:
            
            atom_order = calculate_order_atoms("output.xyz", 15)
            #densities = calculate_param_statewise_color("output.xyz", atom_order, 15)
            densities = calculate_param_nostate("output.xyz", atom_order, 15)
            pos = np.where(np.abs(np.diff(densities[0])) > 0.5)[0] + 1
            time = np.insert(densities[0], pos, np.nan)
            ht = np.insert(densities[1], pos, np.nan)

            densit = calculate_param_statewise("output.xyz", atom_order,15)
            posht1 = np.where(np.abs(np.diff(densit[0])) > 0.5)[0] + 1
            timeht1 = np.insert(densit[0], posht1, np.nan)
            htht1 = np.insert(densit[4], posht1, np.nan)
            for elem in htht1:
                HT1.append(elem) 

            #poss1 = np.where(np.abs(np.diff(densities[2])) > 0.5)[0] + 1
            #times1 = np.insert(densities[2], poss1, np.nan)
            #hts1 = np.insert(densities[3], poss1, np.nan)
            #poss0 = np.where(np.abs(np.diff(densities[4])) > 0.5)[0] + 1
            #times0 = np.insert(densities[4], poss0, np.nan)
            #hts0 = np.insert(densities[5], poss0, np.nan) 

            posht2 = np.where(np.abs(np.diff(densit[9])) > 0.5)[0] + 1
            timeht2 = np.insert(densit[9], posht2, np.nan)
            htht2 = np.insert(densit[10], posht2, np.nan)
            for elem in htht2[:10]:        # 5 fs # first 10 elements
                HT2.append(elem)

            ax1.plot(time, ht, color = clr[1], alpha=0.9, linewidth=0.4, zorder=2)
            #ax1.plot(times2, hts2, color = clr[0], alpha=0.9, linewidth=0.4, zorder=2)
            #ax1.plot(times1, hts1, color = clr[1], alpha=0.9, linewidth=0.4, zorder=0)
            #ax1.plot(times0, hts0, color = clr[4], linewidth=0.4, zorder=1)

            os.chdir(path_data1)    


    ax2.hist(HT1, bins=40, range=[-3,3], density=True, weights=None, orientation='horizontal', color = 'orange', alpha = 0.9, zorder=10)
    ax2.hist(HT2, bins=40, range=[-3,3], density=True, histtype='step', orientation='horizontal', weights=None, color = 'blue', alpha = 0.9, zorder=10)

    ax1.scatter(for21[0], for21[3], color = clr[8], s=0.6, label=r'S$_2$/S$_1$ hops', zorder=5)
    #ax1.scatter(for10[0], for10[3], s=0.6, color = clr[5], label=r'S$_1$/S$_0$ hops', zorder=6)
    ax1.scatter(back21[0], back21[3], s=0.6, color = clr[6], label=r'S$_1$/S$_2$ hops', zorder=7)
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


