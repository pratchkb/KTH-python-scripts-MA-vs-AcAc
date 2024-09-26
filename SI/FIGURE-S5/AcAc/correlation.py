import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.0
import os.path
import glob
# load proper modules
import get_optim
import geom_param
from itertools import islice


# Speed of light in atomic units
h_evs = 4.13566733e-15 #eV * s
c_m_per_s = 299792458.
c_au = 137.035999679
au_per_fs = 41.341374575751
ev_per_au = 27.211386245988


def load_coord(dat):
    
    import get_optim

    data = get_optim.read_xyz_as_np(dat)        # coordinates loaded
    return data


def calculate_distance(atom1_coord, atom2_coord):
    """
    Calculates distance between two bonds
    """

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

def calculate_order_atoms(filename):
    """
    calculate the numbering of atoms depending on which O the H is attahced to, in the initial geometry 
    """

    import get_optim
    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = 15 # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    order = []

    data = get_optim.read_xyz_as_np(filename)        # coordinates loaded

    h1 = geom_param.compute_bond(data,[1,13])
    h2 = geom_param.compute_bond(data,[4,13])

    if h1 > h2:
        order = [1,0,2,3,4]
    else:
        order = [4,3,2,0,1]

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

    num = 15 # number of atoms for this particular system # please change for your system


    #h1 = []
    #h2 = []
    bla4 = []
    blaCO = []
    blaCC = []
             
    data = get_optim.read_xyz_as_np(filename)        # coordinates loaded

    bla4.append(calculate_bla4(data,atom_order)[0])
    blaCO.append(calculate_bla4(data,atom_order)[1])
    blaCC.append(calculate_bla4(data,atom_order)[2])

    fdata.close()

    return bla4, blaCO, blaCC


def read_bagel_xmscaspt2_results(nstates):
  
    files = []
    results_dict = {}
    results_dict_npi = {}
    results_dict_pipi = {}
    results_dict_diff = {}
    results_dict_geom = {}
    #count = 0

    # Iterate through all frameids sampled for absorption within current BOMD frame
    for isample in range(1,1785,1):          # 1785 very specific right now needs to change
   
        #x = str(isample).zfill(4) 
    # Grab excited state energies and oscillator strengths
        results_file = '/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ABSORPTION-updated/CALC-SSSR/' + str(isample) + '/input' + str(isample) + '.out'  #os.getcwd() + '/x{:04d}/{}'.format(isample,outname)
        xyz = '/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ABSORPTION-updated/CALC-SSSR/' + str(isample) + '/' + str(isample) + '.xyz'
        #print(xyz)

        results = []
        results_npi = []
        results_pipi = []
        results_diff = []
        results_geom = []
        energies = []
        energy_npi = []
        energy_pipi = []
        oscs = []
        oscs_npi = []
        oscs_pipi = []
        
        ener_diff = []
        osc_diff = []

        result_S0 = []
        result_S1 = []
        result_S2 = []

        BLA = []
        BLACO = []
        BLACC = []
        angles = []
        pyr = []
        Hbond1 = []
        Hbond2 = []
        torsion1 = []
        torsion2 = []

        dat = load_coord(xyz)
        #print(data)

        #append bla
        #BLA.append(geom_param.compute_bla(data,1,0,2,3))      # need to be changed for other molecules
        #print(BLA)

        #compute sum of angles
        #angles.append(geom_param.compute_angle(data,1,0,2)+geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4))

        #compute pyramidalization
        #pyr.append(geom_param.compute_pyramidalization(data,2,0,3,14))

        #O-H distances
        #Hbond1.append(geom_param.compute_bond(data,[1,13]))
        #Hbond2.append(geom_param.compute_bond(data,[13,4]))

        #torsion
        #torsion1.append(geom_param.compute_torsion(data,1,0,2,3))
        #torsion2.append(geom_param.compute_torsion(data,0,2,3,4))
       
        with open(results_file, 'r') as f_in:
            for line in f_in:

                if "MS-CASPT2 energy : state  0" in line:
                    state1, energy1 = int(line.split()[5]),float(line.split()[6])
                    result_S0.append(energy1)

                if "MS-CASPT2 energy : state  1" in line:
                    state2, energy2 = int(line.split()[5]),float(line.split()[6])
                    result_S1.append(energy2)

                if "MS-CASPT2 energy : state  2" in line:
                    state3, energy3 = int(line.split()[5]),float(line.split()[6])
                    result_S2.append(energy3)

                if ' * Transition    1 - 0 :' in line:
                    osc1 = float(next(f_in, '').strip().split()[4])

                if ' * Transition    2 - 0 :' in line:
                    osc2 = float(next(f_in, '').strip().split()[4])
                    
                    energy_S1 = (energy2 - energy1)*ev_per_au
                    energy_S2 = (energy3 - energy1)*ev_per_au
                  
                    energy_npi.append(energy_S1)
                    oscs_npi.append(osc1)
                    energy_pipi.append(energy_S2)
                    oscs_pipi.append(osc2)
                    energies.append(energy_S1)
                    oscs.append(osc1)
                    energies.append(energy_S2)
                    oscs.append(osc2)
                    ener_diff.append(energy_S2 - energy_S1)
                    #print(ener_diff)
                    osc_diff.append(osc2-osc1)
                    
                    #append bla
                    #BLA.append(geom_param.compute_bla(data,1,0,2,3))      # need to be changed for other molecules
                    #print(BLA)


                    #print(geom_param.compute_angle(data,0,2,3)+geom_param.compute_angle(data,2,3,4)+geom_param.compute_angle(data,3,4,1))
                    #compute sum of angles
                    angles.append(geom_param.compute_angle(dat,1,0,2)+geom_param.compute_angle(dat,0,2,3)+geom_param.compute_angle(dat,2,3,4))
                     
                   
                    #compute pyramidalization
                    pyr.append(geom_param.compute_pyramidalization(dat,2,0,3,14))

                    #O-H distances
                    Hbond1.append(geom_param.compute_bond(dat,[1,13]))
                    Hbond2.append(geom_param.compute_bond(dat,[4,13]))

                    #torsion
                    torsion1.append(geom_param.compute_torsion(dat,0,2,3,4))
                    torsion2.append(geom_param.compute_torsion(dat,3,2,0,1))

                    atom_order = calculate_order_atoms(xyz)
                    #print(atom_order)
                    BLAdensities = calculate_param(xyz, atom_order)
                    #print(BLAdensities)
                    #timestep.append(BLAdensities[0])
                    BLA.append(*BLAdensities[0])
                    BLACO.append(*BLAdensities[1])
                    BLACC.append(*BLAdensities[2])


        for energy, osc in zip(energies, oscs):
            results.append([energy, osc])
            results_dict[isample] = results 

        for energy, osc in zip(energy_npi, oscs_npi):
            results_npi.append([energy, osc])
            results_dict_npi[isample] = results_npi

        for energy, osc in zip(energy_pipi, oscs_pipi):
            results_pipi.append([energy, osc])
            results_dict_pipi[isample] = results_pipi

        for energy, osc, s2 in zip(ener_diff, osc_diff, energy_pipi):
            results_diff.append([energy, osc, s2])
            results_dict_diff[isample] = results_diff

        for energy, ang, pr, h1, h2, t1, t2, bla, blaco, blacc in zip(energy_pipi, angles, pyr, Hbond1, Hbond2, torsion1, torsion2, BLA, BLACO, BLACC):
            results_geom.append([energy, ang, pr, h1, h2, t1, t2, bla, blaco, blacc])
            results_dict_geom[isample] = results_geom

    #print(results_dict_geom)
    return results_dict, results_dict_npi, results_dict_pipi, results_dict_diff, results_dict_geom

#timestep = []
# Iterate through all frameids sampled for absorption within current BOMD frame
#for isample in range(0,5000,1):          # 1785 very specific right now needs to change

#    x = str(isample).zfill(4)
    # Grab excited state energies and oscillator strengths
#    xyz = '../' + str(x) + '/x' + str(x) + '.xyz'
    #print(xyz)

    #data = loadcoord(xyz)

#    atom_order = calculate_order_atoms(xyz)
#    BLAdensities = calculate_param_statewise(xyz, atom_order)
#    timestep.append(BLAdensities[0])
#    BLA.append(BLAdensities[1])
#    BLACO.append(BLAdensities[2])
#    BLACC.append(BLAdensities[3])

# Based on fitting to experiment 
FWHM = 0.24

# General variables
eV_shift = 0.096 ####0.095 #0.101
nsamples = 1785
nstates = 2

# Plot setup
#fig, ax1 = plt.subplots(figsize=(6,4))
#ax2 = ax1.twinx()

data = []
origdir = os.getcwd()

# Read energies and oscillator strengths
data = read_bagel_xmscaspt2_results(nstates)


#save S1 energy and osc into list #data[1] contains dictionary for S1 states
S1en = []
S1osc = []
for idx, states_data in data[1].items():
        for state_data in states_data:
            S1en.append(state_data[0]+eV_shift)
            S1osc.append(state_data[1])

#save S2 energy and osc into list #data[2] contains dictionary for S2 states
S2en = []
S2osc = []
for idx, states_data in data[2].items():
        for state_data in states_data:
            S2en.append(state_data[0]+eV_shift)
            S2osc.append(state_data[1])

#save S2-S1 energy and osc2-osc1 into list #data[3] contains dictionary for S2-S1 states
endiff = []
oscdiff = []
ends2 = []
for idx, states_data in data[3].items():
        for state_data in states_data:
            endiff.append(state_data[0])
            oscdiff.append(state_data[1])
            ends2.append(state_data[2])

#save S2-S1 energy and bla into list #data[3] contains dictionary for S2-S1 states
endiff2 = []

#blas = []
hh1 = []
hh2 = []
tt1 = []
tt2 = []
angs = []
pyrs = []
ht = []
bl = []
blco = []
blcc = []
 
for idx, states_data in data[4].items():
        for state_data in states_data:
            endiff2.append(state_data[0])
            #blas.append(state_data[1])
            angs.append(state_data[1])
            pyrs.append(abs(state_data[2]))
            hh1.append(state_data[3])
            hh2.append(state_data[4])
            tt1.append(abs(state_data[5]))
            tt2.append(abs(state_data[6]))
            ht.append((state_data[3]-state_data[4]))     
            bl.append(state_data[7])
            blco.append(state_data[8])
            blcc.append(state_data[9])


plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('Osc diff', fontsize=16)
plt.ylim(0.0,0.4)
plt.savefig('diff.png', dpi=300)

plt.figure()
plt.hist2d(ends2, oscdiff, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('Osc diff', fontsize=16)
plt.ylim(0.0,0.4)
plt.savefig('s2-oscdiff.png', dpi=300)


#plt.figure()
#plt.scatter(endiff2, blas, color='k')
#plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
#plt.ylabel('BLA', fontsize=16)
#plt.ylim(1.0,1.6)
#plt.savefig('diff-bla.png', dpi=300)

plt.figure()
#plt.scatter(endiff2, angs, color='k')
plt.hist2d(endiff2, angs, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel(r'Sum of Angles ($^{\circ}$)', fontsize=16)
plt.ylim(350,380)
plt.xlim(3.5,6.5)
plt.minorticks_on()
plt.savefig('diff-angs.png', dpi=300)

plt.figure()
plt.hist2d(endiff2, pyrs, bins=100, cmap='gist_heat_r')
plt.xlabel('Energy diff (eV)', fontsize=16)
plt.ylabel('PyrC', fontsize=16)
plt.ylim(-40,40)
plt.savefig('diff-pyr.png', dpi=300)

plt.figure()
plt.hist2d(endiff2, hh1, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('OH1', fontsize=16)
plt.ylim(0.5,2.0)
plt.savefig('diff-h1.png', dpi=300)

plt.figure()
plt.hist2d(endiff2, hh2, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('OH2', fontsize=16)
plt.ylim(0.5,2.0)
plt.savefig('diff-h2.png', dpi=300)

plt.figure()
plt.hist2d(endiff2, tt1, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('tor1', fontsize=16)
plt.ylim(-40,40)
plt.savefig('diff-tor1.png', dpi=300)

plt.figure()
plt.hist2d(endiff2, tt2, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('tor2', fontsize=16)
plt.ylim(-40,40)
plt.savefig('diff-tor2.png', dpi=300)


plt.figure()
#plt.scatter(endiff2, angs, color='k')
plt.hist2d(endiff2, ht, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel('H-transfer coordinate', fontsize=16)
plt.ylim(-1.5,1.5)
plt.xlim(3.5,6.5)
plt.minorticks_on()
plt.savefig('diff-ht.png', dpi=300)


plt.figure()
plt.hist2d(endiff2, bl, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel(r'BLA ($\AA$)', fontsize=16)
plt.ylim(-1.0,1.0)
plt.xlim(3.5,6.5)
plt.minorticks_on()
plt.savefig('diff-BLA.png', dpi=300)


plt.figure()
plt.hist2d(endiff2, blco, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel(r'BLA-CO ($\AA$)', fontsize=16)
plt.ylim(-1.0,1.0)
plt.xlim(3.5,6.5)
plt.minorticks_on()
plt.savefig('diff-BLACO.png', dpi=300)

plt.figure()
plt.hist2d(endiff2, blcc, bins=100, cmap='gist_heat_r')
plt.xlabel(r'S$_2$ energy (eV)', fontsize=16)
plt.ylabel(r'BLA-CC ($\AA$)', fontsize=16)
plt.ylim(-1.0,1.0)
plt.xlim(3.5,6.5)
plt.minorticks_on()
plt.savefig('diff-BLACC.png', dpi=300)



###############################################3
import matplotlib.pyplot as plt
#import matplotlib.colors
#import matplotlib.colors.PowerNorm
#import matplotlib.axes.Axes.hist2d
#import matplotlib.pyplot.hist2d
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator
from scipy import ndimage

mb = endiff2
mt = angs

xedges = np.linspace(3.5,6.5,80)
yedges = np.linspace(350,380,80)
X, Y = np.meshgrid(xedges, yedges)
sigma_x = 0.10
sigma_y = 1.
gau2d = np.exp(-0.5 * (X-5.0)**2/sigma_x**2 - 0.5 * (Y-365)**2/sigma_y**2 )
H, xedges, yedges = np.histogram2d(endiff2, angs, bins=80, density=True, range=[[3.5,6.5],[350,380]])
#
#
#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.
H = H.T
N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)   
# 

fig, ax = plt.subplots(figsize=(5,4))
cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#    cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
#    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
#    fig.colorbar(cbar, ax=ax)
#    #ax.set_aspect('equal')
plt.minorticks_on()
plt.xlim(3.5,6.5)
plt.xticks([3.5,4.0,4.5,5.0,5.5,6.0,6.5])
plt.ylim(350,380)
plt.xlabel(r"S$_2$ energy (eV)", fontsize=18)
plt.ylabel(r"SOA ($^{\circ}$)", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('cont-S2en-contourf-gau-ang.pdf', dpi=300)


mb = endiff2
mt = bl

xedges = np.linspace(3.5,6.5,80)
yedges = np.linspace(-1.0,1.0,80)
X, Y = np.meshgrid(xedges, yedges)
sigma_x = 0.10
sigma_y = 0.067
gau2d = np.exp(-0.5 * (X-5.0)**2/sigma_x**2 - 0.5 * (Y-0.0)**2/sigma_y**2 )
H, xedges, yedges = np.histogram2d(endiff2, bl, bins=80, density=True, range=[[3.5,6.5],[-1.0,1.0]])
#
#
#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.
H = H.T
N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)
# 

fig, ax = plt.subplots(figsize=(5,4))
cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#    cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
#    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
#    fig.colorbar(cbar, ax=ax)
#    #ax.set_aspect('equal')
plt.minorticks_on()
plt.xlim(3.5,6.5)
plt.xticks([3.5,4.0,4.5,5.0,5.5,6.0,6.5])
plt.ylim(-0.6,0.4)
plt.yticks([-0.6,-0.4,-0.2,0.0,0.2,0.4])
plt.xlabel(r"S$_2$ energy (eV)", fontsize=18)
plt.ylabel(r"BLA ($\AA$)", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('cont-S2en-contourf-gau-BLA.pdf', dpi=300)


mb = endiff2
mt = ht

xedges = np.linspace(3.5,6.5,80)
yedges = np.linspace(-1.5,1.5,80)
X, Y = np.meshgrid(xedges, yedges)
sigma_x = 0.10
sigma_y = 0.10
gau2d = np.exp(-0.5 * (X-5.0)**2/sigma_x**2 - 0.5 * (Y-0.0)**2/sigma_y**2 )
H, xedges, yedges = np.histogram2d(endiff2, ht, bins=80, density=True, range=[[3.5,6.5],[-1.5,1.5]])
#
#
#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.
H = H.T
N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)
# 

fig, ax = plt.subplots(figsize=(5,4))
cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#    cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
#    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
#fig.colorbar(cbar, ax=ax)
#    #ax.set_aspect('equal')
plt.minorticks_on()
plt.xlim(3.5,6.5)
plt.xticks([3.5,4.0,4.5,5.0,5.5,6.0,6.5])
plt.ylim(-1.5,1.5)
plt.yticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
plt.xlabel(r"S$_2$ energy (eV)", fontsize=18)
plt.ylabel(r"HT ($\AA$)", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('cont-S2en-contourf-gau-ht.pdf', dpi=300)


#############################################
mb = endiff
mt = oscdiff

xedges = np.linspace(0,4.4,80)
yedges = np.linspace(0.0,0.4,80)
X, Y = np.meshgrid(xedges, yedges)
sigma_x = 0.10
sigma_y = 0.01
gau2d = np.exp(-0.5 * (X-2.2)**2/sigma_x**2 - 0.5 * (Y-0.2)**2/sigma_y**2 )
H, xedges, yedges = np.histogram2d(endiff, oscdiff, bins=80, density=True, range=[[0,4.4],[0.0,0.4]])
#
#
#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.
H = H.T
N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)
# 

fig, ax = plt.subplots(figsize=(5,4))
cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#    cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
#    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
fig.colorbar(cbar, ax=ax)
#    #ax.set_aspect('equal')
plt.minorticks_on()
plt.xlim(0,2.5)
plt.xticks([0.0,0.5,1.0,1.5,2.0,2.5])
plt.ylim(0.0,0.4)
plt.yticks([0.0,0.1,0.2,0.3,0.4])
plt.xlabel(r"S$_2$-S$_1$ energy difference (eV)", fontsize=18)
plt.ylabel(r"Oscillator strength", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('cont-S2en-contourf-gau-osc.pdf', dpi=300)

#############################################
mb = ends2
mt = oscdiff

xedges = np.linspace(3.5,6.5,80)
yedges = np.linspace(0.0,0.4,80)
X, Y = np.meshgrid(xedges, yedges)
sigma_x = 0.10
sigma_y = 0.01
gau2d = np.exp(-0.5 * (X-5.0)**2/sigma_x**2 - 0.5 * (Y-0.2)**2/sigma_y**2 )
H, xedges, yedges = np.histogram2d(ends2, oscdiff, bins=80, density=True, range=[[3.5,6.5],[0.0,0.4]])
#
#
#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.
H = H.T
N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)
# 

fig, ax = plt.subplots(figsize=(5,4))
cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#    cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
#    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
#fig.colorbar(cbar, ax=ax)
#    #ax.set_aspect('equal')
plt.minorticks_on()
plt.xlim(3.5,6.5)
plt.xticks([3.5,4.0,4.5,5.0,5.5,6.0,6.5])
#plt.xticks([0.0,1.0,2.0,3.0,4.0])
plt.ylim(0.0,0.4)
plt.yticks([0.0,0.1,0.2,0.3,0.4])
#plt.yticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
plt.xlabel(r"S$_2$ energy (eV)", fontsize=18)
plt.ylabel(r"Oscillator strength", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.tight_layout()
plt.savefig('cont-S2-contourf-gau-osc.pdf', dpi=300)

