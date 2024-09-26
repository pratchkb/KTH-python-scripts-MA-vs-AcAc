import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.pyplot as plt
#import os
import sys
sys.path.append('/Users/pratip/')
from aimsrelated.tcutil.code.utils import units
from aimsrelated.tcutil.code.atom_data import atom_data
from aimsrelated.tcutil.code.get_optim import get_optim
from aimsrelated.tcutil.code.geom_param import geom_param
from itertools import islice



#CACHE_FILENAME = "pop.pickle"
def calculate_distance(atom1_coord, atom2_coord):
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12

def calculate_order_atoms_torsion(filename):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the final geometry to back track 
    """
    #import get_optim
    #import geom_param
    #from itertools import islice
    import matplotlib.pyplot as plt
    #import os
    import sys
    sys.path.append('/Users/pratip/')
    from aimsrelated.tcutil.code.utils import units
    from aimsrelated.tcutil.code.atom_data import atom_data
    from aimsrelated.tcutil.code.get_optim import get_optim
    from aimsrelated.tcutil.code.geom_param import geom_param
    from itertools import islice


    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = 15 # number of atoms for this particular system # please change for your system

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




def collect_geom_param(filename):

    import numpy as np
    import os
    import os.path
    import glob
    # load proper modules
    #import get_optim
    #import geom_param
    #from itertools import islice
    import matplotlib.pyplot as plt
    #import os
    import sys
    sys.path.append('/Users/pratip/')
    from aimsrelated.tcutil.code.utils import units
    from aimsrelated.tcutil.code.atom_data import atom_data
    from aimsrelated.tcutil.code.get_optim import get_optim
    from aimsrelated.tcutil.code.geom_param import geom_param
    from itertools import islice


    order1 = calculate_order_atoms_torsion(filename)
    #print(order1)
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = 15 # number of atoms for this particular system # please change for your system


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

#def load_cache(filename):
#    try:
#        with open(filename, 'rb') as f:
#            return pickle.load(f)
#    except:
#        return {}

#def save_cache(filename: str, data):
#    """ Will save the data into a pickle file
#
#    Args:
#        filename : string : the filename in which to write
#        data : 
#    """
#    try:
#        with open(filename, 'wb') as f:
#            pickle.dump(data, f)
#    except Exception as e:
#        print(f"damn - could not save, {e}")

#def load_or_process_population_data(path):
#    raw_data = load_cache(CACHE_FILENAME)
    # Read all the raw data
#    raw_data = read_population_data(data_path, raw_data)
#    save_cache(CACHE_FILENAME, raw_data)
#    return raw_data


if __name__ == "__main__":

#    data_path = '/cfs/klemming/projects/snic/x_nanho/users/pratip/ACETYLACETONE/PROD_Final_Feb27/pop_data/'
#    raw_data = load_or_process_population_data(data_path)

    for21 = collect_geom_param('crossing-s2-s1.xyz')
    for10 = collect_geom_param('crossing-s1-s0.xyz')
    back21 = collect_geom_param('crossing-s1-s2.xyz')

    #plots
    path = '/Users/pratip/PROD_Final_Feb27/'
    

    fig, ax = plt.subplots(1,1, figsize=(4,3))#, sharex=True)

    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            #print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], (fdata[:,4]-fdata[:,3]), color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)


    ax.scatter(for21[0], for21[3], s=2.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[3], s=2.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    #ax.scatter(back21[0], back21[3], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"H-transfer coordinate ($\AA$)")
    ax.set_xlabel(r"Time (fs)")
    ax.set_ylim(-4.0,4.0)
    #ax.set_yticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('H-transfer-coord-again.png', dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(5,3.5), sharex=True)

    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], fdata[:,7], color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    ax.scatter(for21[0], for21[4], s=2.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[4], s=2.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    #ax.scatter(back21[0], back21[4], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"PyrC ($^{\circ}$)")
    ax.set_xlabel(r"Time (fs)")
    ax.set_ylim(-100.0,100.0)
    ax.set_yticks([-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('pyrC-coord.png', dpi=500)


    # dpyr/dt
    fig, ax = plt.subplots(1,1, figsize=(5,3.5), sharex=True)

    for i in range(3,4):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            dx = 0.5
            deriv = np.gradient(fdata[:,7], dx) 
            ax.plot(fdata[:,0], deriv, color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    #ax.scatter(for21[0], for21[4], s=2.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    #ax.scatter(for10[0], for10[4], s=2.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    #ax.scatter(back21[0], back21[4], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"deriv PyrC ($^{\circ}$)")
    ax.set_xlabel(r"Time (fs)")
    #ax.set_ylim(-100.0,100.0)
    #ax.set_yticks([-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('deriv-pyrC.png', dpi=500)


    fig, ax = plt.subplots(1,1, figsize=(4,3))#, sharex=True)

    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], fdata[:,2], color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    ax.scatter(for21[0], for21[5], s=2.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[5], s=2.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    #ax.scatter(back21[0], back21[5], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"Sum of Angles ($^{\circ}$)")
    ax.set_xlabel(r"Time (fs)")
    #ax.set_ylim(-2.5,2.5)
    #ax.set_yticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('SumAngle-coord.png', dpi=500)


    fig, ax = plt.subplots(1,1, figsize=(4,3))#, sharex=True)

    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], fdata[:,11], color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    ax.scatter(for21[0], for21[23], s=1.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[23], s=1.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    ax.scatter(back21[0], back21[23], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"Torsion ($^{\circ}$)")
    ax.set_xlabel(r"Time (fs)")
    #ax.set_ylim(-2.5,2.5)
    #ax.set_yticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('Torsion-coord-AcAc.png', dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))#, sharex=True)

    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], fdata[:,5], color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    ax.scatter(for21[0], for21[24], s=1.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[24], s=1.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    ax.scatter(back21[0], back21[24], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"Torsion1 ($^{\circ}$)")
    ax.set_xlabel(r"Time (fs)")
    #ax.set_ylim(-2.5,2.5)
    #ax.set_yticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('Torsion1-coord-AcAc.png', dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))#, sharex=True)

    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], fdata[:,6], color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    ax.scatter(for21[0], for21[25], s=1.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[25], s=1.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    ax.scatter(back21[0], back21[25], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"Torsion2 ($^{\circ}$)")
    ax.set_xlabel(r"Time (fs)")
    #ax.set_ylim(-2.5,2.5)
    #ax.set_yticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('Torsion2-coord-AcAc.png', dpi=500)


    for i in range(1,272):
        i = "{:05d}".format(i)
        os.chdir('TRAJ_'+str(i))
        if os.path.exists('DONT_ANALYZE'):
            print(str(i))
            os.chdir(path)
        else:
            fdata=np.genfromtxt('geometric_param_pos1.dat',dtype='float64')
            ax.plot(fdata[:,0], fdata[:,8], color = 'orange', linewidth=0.4, zorder=0)
            os.chdir(path)

    #print(for21[22])
    ax.scatter(for21[0], for21[22], s=1.0, color = 'darkred', label=r'S$_2$/S$_1$ hops', zorder=1)
    ax.scatter(for10[0], for10[22], s=1.0, color = 'blue', label=r'S$_1$/S$_0$ hops', zorder=3)
    ax.scatter(back21[0], back21[22], s=1.0, color = 'limegreen', label=r'S$_1$/S$_2$ hops', zorder=2)
    ax.set_ylabel(r"BLA4 ($\AA$)")
    ax.set_xlabel(r"Time (fs)")
    ax.set_ylim(-1,1)
    #ax.set_yticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax.set_xlim(0,200)
    ax.minorticks_on()
    #plt.legend(frameon=False, loc = 'upper left')
    plt.tight_layout()
    plt.savefig('BLA4-coord.png', dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))
    plt.hist(for21[3], bins=30, density=True, range=[-2.5,2.5], weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops')
    plt.hist(for10[3], bins=30, density=True, range=[-2.5,2.5], histtype='step', weights=None, color = 'blue', alpha = 0.9, label=r'S$_1$/S$_0$ hops')
    plt.hist(back21[3], bins=30, density=True, range=[-2.5,2.5], histtype='step', weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
    plt.xlabel(r'H-transfer coordinate ($\AA$)')
    plt.ylabel('Normalized distribution')
    plt.xlim(-2.0,2.0)
    plt.ylim(0,2.2)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    plt.legend(loc ="upper right", frameon=False)
    plt.tight_layout()
    ax.minorticks_on()
    #plt.savefig('H-transfer-histogram-again.png',dpi=500)
    plt.savefig('H-transfer-histogram-count.png',dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))
    plt.hist(for21[4], bins=50, range=[-90.0,90.0], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops')
    plt.hist(for10[4], bins=50, range=[-90.0,90.0], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_1$/S$_0$ hops')
    plt.hist(back21[4], bins=50, range=[-90.0,90.0], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')

    plt.xlabel(r'PyrC ($^{\circ}$)')
    plt.ylabel('Normalized Distribution')
    plt.xlim(-90.0,90.0)
    plt.ylim(0,0.07)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    plt.xticks([-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0])
    plt.legend(loc ="upper right", frameon=False)
    ax.minorticks_on()
    plt.tight_layout()
    plt.savefig('pyrC-histogram-again.png',dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))
    plt.hist(for21[5], bins=40, range=[320,400], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops')
    plt.hist(for10[5], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_1$/S$_0$ hops')
    plt.hist(back21[5], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
    plt.xlabel(r'Sum of Angles ($^{\circ}$)')
    plt.ylabel('Normalized Distribution')
    plt.xlim(320,400)
    ax.minorticks_on()
    plt.ylim(0,0.12)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    plt.legend(loc ="upper right", frameon=False)
    plt.tight_layout()
    plt.savefig('SumAngle-histogram-again.png',dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))
    plt.hist(for21[22], bins=40, range=[-1,1], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops')
    plt.hist(for10[22], bins=40, range=[-1,1], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_1$/S$_0$ hops')
    plt.hist(back21[22], bins=40, range=[-1,1], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
    plt.xlabel(r'BLA4 ($\AA$)')
    plt.ylabel('Normalized Distribution')
    plt.xlim(-1,1.0)
    ax.minorticks_on()
    plt.ylim(0,10.0)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    plt.legend(loc ="upper right", frameon=False)
    plt.tight_layout()
    plt.savefig('BLA4-histogram-again.png',dpi=500)

    # Sum of Angles seperated
    fig, ax = plt.subplots(2,1, figsize=(4.5,5.5), sharex=True)
    ax[0].hist(for21[9], bins=40, range=[320,400], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops < 20 fs')
    ax[0].hist(for21[13], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_2$/S$_1$ hops >= 20 fs')
    ax[1].hist(for10[17], bins=40, range=[320,400], density=True, weights=None, color = 'green', alpha = 0.9, label=r'S$_1$/S$_0$ hops < 65 fs')
    ax[1].hist(for10[21], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'red', alpha = 0.9, label=r'S$_1$/S$_0$ hops >= 65 fs')
    #plt.hist(back21[5], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
    ax[1].set_xlabel(r'Sum of Angles ($^{\circ}$)')
    ax[0].set_ylabel('Normalized Distribution')
    ax[1].set_ylabel('Normalized Distribution')
    ax[0].set_xlim(320,400)
    ax[1].set_xlim(320,400)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[0].set_ylim(0,0.16)
    ax[1].set_ylim(0,0.16)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    ax[0].legend(loc ="upper right", frameon=False, fontsize=10)
    ax[1].legend(loc ="upper right", frameon=False, fontsize=10)
    plt.tight_layout()
    plt.savefig('SumAngle-histogram-seperated.png',dpi=500)

    # H coordinate seperated
    fig, ax = plt.subplots(2,1, figsize=(4.5,5.5), sharex=True)
    ax[0].hist(for21[7], bins=30, range=[-2.5,2.5], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops < 20 fs')
    ax[0].hist(for21[11], bins=30, range=[-2.5,2.5], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_2$/S$_1$ hops >= 20 fs')
    ax[1].hist(for10[15], bins=30, range=[-2.5,2.5], density=True, weights=None, color = 'green', alpha = 0.9, label=r'S$_1$/S$_0$ hops < 65 fs')
    ax[1].hist(for10[19], bins=30, range=[-2.5,2.5], histtype='step', density=True, weights=None, color = 'red', alpha = 0.9, label=r'S$_1$/S$_0$ hops >= 65 fs')
    #plt.hist(back21[5], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
    ax[1].set_xlabel(r'H-transfer coordinate ($\AA$)')
    ax[0].set_ylabel('Normalized Distribution')
    ax[1].set_ylabel('Normalized Distribution')
    ax[0].set_xlim(-2.5,2.5)
    ax[1].set_xlim(-2.5,2.5)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[0].set_ylim(0,3.5)
    ax[1].set_ylim(0,3.5)

    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    ax[0].set_xticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax[1].set_xticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
    ax[0].legend(loc ="upper right", frameon=False, fontsize=10)
    ax[1].legend(loc ="upper right", frameon=False, fontsize=10)
    plt.tight_layout()
    plt.savefig('H-coord-histogram-seperated.png',dpi=500)


    # PyrC seperated
    fig, ax = plt.subplots(2,1, figsize=(4.5,5.5), sharex=True)
    ax[0].hist(for21[8], bins=50, range=[-90,90], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops < 20 fs')
    ax[0].hist(for21[12], bins=50, range=[-90,90], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_2$/S$_1$ hops >= 20 fs')
    ax[1].hist(for10[16], bins=50, range=[-90,90], density=True, weights=None, color = 'green', alpha = 0.9, label=r'S$_1$/S$_0$ hops < 65 fs')
    ax[1].hist(for10[20], bins=50, range=[-90,90], histtype='step', density=True, weights=None, color = 'red', alpha = 0.9, label=r'S$_1$/S$_0$ hops >= 65 fs')
    #plt.hist(back21[5], bins=40, range=[320,400], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
    ax[1].set_xlabel(r'PyrC ($^{\circ}$)')
    ax[0].set_ylabel('Normalized Distribution')
    ax[1].set_ylabel('Normalized Distribution')
    ax[0].set_xlim(-90.0,90.0)
    ax[1].set_xlim(-90.0,90.0)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[0].set_ylim(0,0.08)
    ax[1].set_ylim(0,0.08)

    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    ax[0].set_xticks([-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0])
    ax[1].set_xticks([-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0])
    ax[0].legend(loc ="upper right", frameon=False, fontsize=10)
    ax[1].legend(loc ="upper right", frameon=False, fontsize=10)
    plt.tight_layout()
    plt.savefig('pyrC-histogram-seperated.png',dpi=500)



    fig, ax = plt.subplots(1,1, figsize=(4,3))
    plt.hist(for21[0], bins=50, range=[0,200], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops')
    plt.hist(for10[0], bins=50, range=[0,200], histtype='step', density=True, weights=None, color = 'blue', alpha = 0.9, label=r'S$_1$/S$_0$ hops')
    plt.hist(back21[0], bins=50, range=[0,200], histtype='step', density=True, weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')

    plt.xlabel(r'Time (fs)')
    plt.ylabel('Normalized Distribution')
    plt.xlim(0,200)
    plt.ylim(0,0.1)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    plt.legend(loc ="upper right", frameon=False)
    plt.tight_layout()
    ax.minorticks_on()
    plt.savefig('H-transfer-histogram-time.png',dpi=500)

    fig, ax = plt.subplots(1,1, figsize=(4,3))
    plt.hist(for21[0], bins=50, range=[0,200], weights=None, color = 'orange', alpha = 0.9, label=r'S$_2$/S$_1$ hops')
    plt.hist(for10[0], bins=50, range=[0,200], histtype='step', weights=None, color = 'blue', alpha = 0.9, label=r'S$_1$/S$_0$ hops')
    plt.hist(back21[0], bins=50, range=[0,200], histtype='step', weights=None, color = 'limegreen', alpha = 0.9, label=r'S$_1$/S$_2$ hops')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')

    plt.xlabel(r'Time (fs)')
    plt.ylabel('Count')
    plt.xlim(0,200)
    plt.ylim(0,70)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
    plt.legend(loc ="upper right", frameon=False)
    plt.tight_layout()
    ax.minorticks_on()
    plt.savefig('H-transfer-count-histogram-time.png',dpi=500)



