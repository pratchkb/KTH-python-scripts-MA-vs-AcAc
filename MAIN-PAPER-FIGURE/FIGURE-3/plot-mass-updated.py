import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import os.path
from itertools import islice
# load proper modules
import get_optim
import geom_param

# activate latex text rendering
#rc('text', usetex=True)
def calculate_distance(atom1_coord, atom2_coord):
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12

def load_coord(dat):
    """
    loads coordinates from an xyz file named dat
    """

    import get_optim

    data = get_optim.read_xyz_as_np(dat)        # coordinates loaded
    return data

def calculate_bla4(data,atom_order):
    """
    Calculates all 3 types of BLA from the xyz data file and the list of atom order
    """

    BLA4 = calculate_distance(data[atom_order[0]],data[atom_order[1]]) - calculate_distance(data[atom_order[1]],data[atom_order[2]]) + calculate_distance(data[atom_order[2]],data[atom_order[3]]) - calculate_distance(data[atom_order[3]],data[atom_order[4]])
    BLACO = calculate_distance(data[atom_order[0]],data[atom_order[1]]) - calculate_distance(data[atom_order[3]],data[atom_order[4]])
    BLACC = - calculate_distance(data[atom_order[1]],data[atom_order[2]]) + calculate_distance(data[atom_order[2]],data[atom_order[3]])

    return BLA4, BLACO, BLACC

def calculate_order_atoms(filename):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry
    """

    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = 9 # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    order = []

    os.system('mkdir tmp')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "Frame " in line: # Initial condition geometry
             t = line.split()
             t = t[1]
             l = float(t) #* au_to_fs    # timestep converted to fs
             timstp.append(l)           # timestep appended
             dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

             h1 = geom_param.compute_bond(data,[0,8])
             h2 = geom_param.compute_bond(data,[1,8])

             if h1 > h2:
                 order = [0,2,3,4,1]
             else:
                 order = [1,4,3,2,0]

    fdata.close()

    os.system('rm -r tmp')

    return order


def calculate_param(filename, atom_order):
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

    num = 9 # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    #h1 = []
    #h2 = []
    bla4 = []
    blaCO = []
    blaCC = []
    Htran = []
    pyr = []
    angles = []
    tor = []

    os.system('mkdir tmp2')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "Frame " in line:
             t = line.split()
             t = t[1]
             l = float(t) #* au_to_fs    # timestep converted to fs
             timstp.append(l)           # timestep appended
             dat = 'tmp2/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp2/' + str(t) + '.xyz')        # coordinates loaded

             h1 = geom_param.compute_bond(data,[0,8])
             h2 = geom_param.compute_bond(data,[1,8])

             if h1 > h2: 
                 tor.append(geom_param.compute_torsion(data,2,3,4,1))
             else:
                 tor.append(geom_param.compute_torsion(data,4,3,2,0))

             bla4.append(calculate_bla4(data,atom_order)[0])
             blaCO.append(calculate_bla4(data,atom_order)[1])
             blaCC.append(calculate_bla4(data,atom_order)[2])
             Htran.append((geom_param.compute_bond(data,[1,8])-geom_param.compute_bond(data,[0,8])))
             pyr.append(geom_param.compute_pyramidalization(data,3,2,4,6))
             angles.append(geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4)+geom_param.compute_angle(data,3,4,1))


    fdata.close()

    os.system('rm -r tmp2')

    return timstp, bla4, blaCO, blaCC, Htran, pyr, angles, tor

if __name__ == "__main__":
    #singlets

    path = "/data/projects/Pratip_MA_vs_AcAc/MA-data/2d-plot-data/MALONALDEHYDE/"
    pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/MAIN-PAPER-FIGURE/FIGURE-3/"

    time = []
    bondla = []
    ht = []
    pyra = []
    ang = []
    torr = []

    # calling functions
    for i in range(0,21):
        atom_order = calculate_order_atoms(path + str(i) + "/" + str(i) + ".xyz")
        densities = calculate_param(path + str(i) + "/" + str(i) + ".xyz", atom_order)
    #values=calculate_param("output.xyz",atom_order
        time.append(densities[0])
        bondla.append(densities[1])
        ht.append(densities[4])
        pyra.append(densities[5])
        ang.append(densities[6])
        torr.append(densities[7])

    print(ht)
    data = np.genfromtxt(path + 'final-result-au-nm.dat', dtype='float64')

    #s0min energy staggered
    num = -266.4219620937 

    fig, ax = plt.subplots(3,1,figsize=(5,5.5), sharex=True, gridspec_kw={'height_ratios': [3,1.5,1]})
    ax[0].plot(ht,(data[:,1]-num)*27.2114, color='orange', label=r'S$_0$')
    ax[0].plot(ht,(data[:,2]-num)*27.2114, color='blue', label='S$_1$') ##eb6d5cff
    ax[0].plot(ht,(data[:,3]-num)*27.2114, color='limegreen', label='S$_2$')
    ax[0].tick_params(axis='both', labelsize=12)
    ax[0].minorticks_on()
    ax[0].set_ylim(0,6)
    ax[0].set_ylabel("Relative energy (eV)", fontsize=12)
    #ax[0].set_xlabel("Linear Interpolation")
#ax[1].plot(charges[:,0],charges[:,1], color='purple', linestyle="--", label='O')
#ax[1].plot(charges[:,0],charges[:,2], color='blue', linestyle="--", label='Central')
#ax[1].plot(charges[:,0],charges[:,3], color='pink', linestyle="--", label='Carbonyl')
#ax[0].plot(charges[:,0],charges[:,4]-charges[:,2], color='#B34C57', linestyle="--",zorder=1)
#ax[0].axhline(y=0.0, color='black', alpha=0.4, linewidth=0.3, zorder=0)
#ax[0].set_ylabel(r'S$_2$ $\Delta$q ($\it{e}$)', fontsize=12)
#ax[0].set_ylim(-1.1,1.1)
#ax[0].set_xlim([0.0,0.6])
#ax[0].set_xticks([0.0,0.5,1.0,1.5,2.0])
#ax[0].set_yticks([-1.0,-0.5,0.0,0.5,1.0])
#ax[0].tick_params(labelsize=12)

    ax[0].legend(loc="upper center", ncol=3, fontsize=12, frameon=False)

    ax[1].tick_params(axis='both', labelsize=12)
    ax[1].minorticks_on()
    ax[1].set_ylim(350,370)
    ax[1].set_ylabel(r'SOA ($^{\circ}$)', fontsize=12)
    ax[1].plot(ht, ang, color='k')

    ax[2].plot(ht, bondla, color='k')
    ax[2].tick_params(axis='both', labelsize=12)
    ax[2].minorticks_on()
    ax[2].set_ylim(-0.2,0.1)
    ax[2].set_yticks([-0.2,-0.1,0.0,0.1])
    ax[2].set_ylabel(r'BLA ($\AA$)', fontsize=12)

    
    #ax[3].plot(ht, pyra, color= 'gray', linestyle='--', label='PyrC', zorder=1)
    #ax[3].plot(ht, torr, color= 'gray', linestyle=':', label='Torsion', zorder=2)
    #ax[3].set_ylim(-20,20)
    #ax[3].minorticks_on()
    #ax[3].legend(loc='upper center', ncol=2, fontsize=10, frameon=False)
    #ax[3].tick_params(axis='both', labelsize=10)
    #ax[3].set_ylabel(r'Angles ($^{\circ}$)', fontsize=10)
    ax[2].set_xlabel(r'HT ($\AA$)', fontsize=12)

#ax[0].legend(loc="upper left", ncol=2, fontsize=10, frameon=False)
    plt.tight_layout()
    #plt.xlabel(r'H-transfer coordinate ($\AA$)')
    plt.xlim(-0.6664973812370789,0.6664973812370789)
    plt.subplots_adjust(wspace=0, hspace=0.10)
    plt.savefig(pathscript + 'Trial_malonaldehyde_Htransfer_linearint.pdf', dpi=1000)


#fig, ax = plt.subplots(figsize=(5.7,2.5))
#ax1=ax.twinx()
#ax.plot(geom[:,0], geom[:,4], label='Torsion', linestyle='-', color = 'black', zorder=1)
#ax.plot(geom[:,0], geom[:,1], label='PyrT', linestyle='--',color='gray',zorder=2)
#ax.plot(geom[:,0], geom[:,2], label='PyrC', linestyle=':',color = 'silver',zorder=3)
#ax.set_ylim(-100,100)
#ax.set_yticks([-90,-45,0,45,90])
#ax.set_ylabel(r'Angle ($\degree$)', fontsize=12)
#ax.set_xlim(0.0,0.6)
#ax.set_xticks([0.0,0.5,1.0,1.5,2.0])
#ax.minorticks_on()
#ax.yaxis.set_tick_params(which='minor', bottom=False)
#ax.set_xlabel(r"Mass-weighted distance (amu$^{1/2}\AA$)", fontsize=12)
#ax.tick_params(labelsize=12)

#ax1.plot(geom[:,0], geom[:,3], label='BLA', linestyle='-.', color = '#B34C57',zorder=4)
#ax.legend(loc="upper left", ncol=3, fontsize=10, frameon=False)
#ax1.set_ylim(1.0,1.5)
#ax1.set_yticks([1.0,1.1,1.2,1.3,1.4,1.5])
#ax1.tick_params(labelsize=12)
#ax1.tick_params(axis='y', colors='#B34C57')
#ax1.set_ylabel(r'BLA ($\AA$)', fontsize=12)
#ax1.yaxis.label.set_color('#B34C57')
#ax1.legend(loc="lower left", ncol=2, fontsize=10, frameon=False)
#plt.tight_layout()
#plt.savefig('AC_S0min_to_S1min_interpolated_geom_updated.png', dpi=300)
#plt.savefig('AC_S0min_to_S1min_interpolated_geom_updated.pdf', dpi=300)

