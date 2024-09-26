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


def calculate_param_nostate(filename, natoms, hh1, hh2):

    num = natoms # number of atoms for this particular system # please change for your system

    data = get_optim.read_xyzs_as_np(filename)
    time = []
    htran = []
    for i in range(len(data[:])):
        time.append(i*0.5)
        h1 = geom_param.compute_bond(data[i],hh2)
        h2 = geom_param.compute_bond(data[i],hh1)
        htran.append(h1-h2)

    return time, htran, data

def collect_geom_param_hops(filename, traj, natoms):

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
        clr = ['#0072BD', 'orange', '#66BF2E', 'k', 'green', 'darkgreen', 'darkblue', 'orangered', '#D45500', 'crimson', 'lime', 'cyan']
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])


        os.chdir(run_path)
        IC = collect_geom_param(run_path + 'total-FC.xyz', numatoms, hh1, hh2)  #ICs
        #for21 = collect_geom_param(run_path + 'crossing-s2-s1.xyz', numatoms, hh1, hh2)
        #for10 = collect_geom_param(runpath + 'crossing-s1-s0.xyz', numatoms, hh1, hh2)
        #back21 = collect_geom_param(run_path + 'crossing-s1-s2.xyz', numatoms, hh1, hh2)

        HT1 = []
        HT2 = []

        time_total = []
        density_total = []

        t21 = []
        s21 = []
        
        t12 = []
        s12 = []



        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir('TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):
                os.chdir(run_path)
            else:
            
                #atom_order = calculate_order_atoms("output.xyz", molecule, numatoms)
            	#densities = calculate_param_statewise_color("output.xyz", atom_order, 15)
                densities = calculate_param_nostate("output.xyz", numatoms, hh1, hh2)
                pos = np.where(np.abs(np.diff(densities[0])) > 0.5)[0] + 1
                time = np.insert(densities[0], pos, np.nan)
                ht = np.insert(densities[1], pos, np.nan)
                ax1.plot(time, ht, color = clr[1], alpha=0.9, linewidth=0.4, zorder=0)

                for time,density in zip(densities[0], densities[1]):
                    time_total.append(time)
                    density_total.append(density)
              
                times21 = collect_geom_param_hops(run_path + 'crossing-s2-s1.xyz','TRAJ_'+str(i),numatoms)
                #print(for21)
                times10 = collect_geom_param_hops(run_path + 'crossing-s1-s0.xyz','TRAJ_'+str(i),numatoms)
                #print(for10)
                times12 = collect_geom_param_hops(run_path + 'crossing-s1-s2.xyz','TRAJ_'+str(i),numatoms)

                tim21 = []
                bla21 = []
                for time in times21:
                    frame = int(time*2)
                    tim21.append(time)
                    bla21.append(geom_param.compute_bond(densities[2][frame],hh2)-geom_param.compute_bond(densities[2][frame],hh1))
                    t21.append(time)
                    s21.append(geom_param.compute_bond(densities[2][frame],hh2)-geom_param.compute_bond(densities[2][frame],hh1))


                tim12 = []
                bla12 = []
                for time in times12:
                    frame = int(time*2)
                    tim12.append(time)
                    bla12.append(geom_param.compute_bond(densities[2][frame],hh2)-geom_param.compute_bond(densities[2][frame],hh1))   
                    t12.append(time)
                    s12.append(geom_param.compute_bond(densities[2][frame],hh2)-geom_param.compute_bond(densities[2][frame],hh1))



                ax1.scatter(tim21, bla21, s=0.6, c=clr[9], zorder=5)
                ax1.scatter(tim12, bla12, s=0.6, c=clr[6], zorder=15)

 
                densit = calculate_param_statewise("output.xyz", numatoms)
                tim2 = []
                ht2 = []
                for time in densit[0]:
                    tim2.append(time)
                    ht2.append(geom_param.compute_bond(densities[2][int(time*2.0)],hh2)-geom_param.compute_bond(densities[2][int(time*2.0)],hh1))                           
                for elem in ht2:
                    HT1.append(elem)

                tim1 = []
                ht1 = []
                for time in densit[1]:
                    tim1.append(time)
                    ht1.append(geom_param.compute_bond(densities[2][int(time*2.0)],hh2)-geom_param.compute_bond(densities[2][int(time*2.0)],hh1))
                for elem in ht1[:10]:      # first 5 fs
                    HT2.append(elem)


                #IC = collect_geom_param(run_path + 'total-FC.xyz', numatoms, hh1, hh2)  #ICs

                os.chdir(run_path)

        ax2.hist(IC[1], bins=50, range=[-3,3], density=True, weights=None, orientation='horizontal', color = 'orange', alpha = 0.9, zorder=10)
        ax2.hist(HT1, bins=50, range=[-3,3], density=True, histtype='step', orientation='horizontal', weights=None, color = 'blue', alpha = 0.9, zorder=11)
        ax2.hist(HT2, bins=50, range=[-3,3], density=True, histtype='step', orientation='horizontal', weights=None, color = 'limegreen', alpha = 0.9, zorder=12)


        #ax1.scatter(for21[0], for21[1], color = clr[8], s=0.6, label=r'S$_2$/S$_1$ hops', zorder=60)
    	#ax1.scatter(for10[0], for10[1], s=0.6, color = clr[5], label=r'S$_1$/S$_0$ hops', zorder=61)
        #ax1.scatter(back21[0], back21[1], s=0.6, color = clr[6], label=r'S$_1$/S$_2$ hops', zorder=62)
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
        ax1.set_xticks([0,25,50,75,100,125,150,175])
        ax1.set_ylim([-3.0,3.0])
        ax2.set_ylim([-3.0,3.0])
        #ax2.set_ticklabels([])
        ax2.set_yticks([])
        ax2.set_xlim([0.0,1.5])
        ax2.set_xticks([0.0,0.5,1.0,1.5])
        ax1.minorticks_on()
        ax2.minorticks_on()
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.05)
        plt.savefig(path_script + f'{molecule}-H-transfer-density-colored-updated.png', dpi=500)
        plt.savefig(path_script + f'{molecule}-H-transfer-density-colored-updated.pdf', dpi=500)

       
        ####################################################
        mb = np.mean(time_total)
        mt = np.mean(density_total)

        f=open(path_script + f'{molecule}-data.dat','w')
        for i,j in zip(time_total, density_total):
            f.write('{}\t{}\n'.format(i,j))
        f.close()

        data = np.genfromtxt(path_script + f'{molecule}-data.dat',dtype='float64')

        xedges = np.linspace(0,200,100)
        yedges = np.linspace(-3.0,3.0,100)
        X, Y = np.meshgrid(xedges, yedges)

        H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=100, density=True, range=[[0,200],[-3.0,3.0]])

        # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
        #    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
        #    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
        #    # therefore transpose H for visualization purposes.
        H = H.T

        #print(len(H))
        #print(len(N_gau_conv))
#        fig, ax = plt.subplots(figsize=(7,5))
#        cbar = ax.contourf(X,Y,H,cmap='gist_heat_r',levels=np.linspace(H.min(),H.max(),51))
#        fig.colorbar(cbar, ax=ax)
        #cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')  # levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#        ax.set_xlim(0,200)
#        ax.set_ylim(-3.0,3.0)
#        ax.set_xlabel(r"Time (fs)", fontsize=18)
#        ax.set_ylabel(r"HT ($\AA$)", fontsize=18)
        #ax.set_xticks(fontsize=18)
        #ax.set_yticks(fontsize=18)
#        fig.tight_layout()
#        fig.savefig(path_script + f'{molecule}-contourf-density.pdf', dpi=300)

        
        sigma_x = 1.0
        sigma_y = 0.01
        gau2d = np.exp(-0.5 * (X-100)**2/sigma_x**2 - 0.5 * (Y)**2/sigma_y**2 )
        N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)

# 
        fig, ax = plt.subplots(figsize=(7,5))
        cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(0.00,0.0001716,51), zorder=10) #levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
        fig.colorbar(cbar, ax=ax)
        plt.xlim(0,200)
        plt.ylim(-3.0,3.0)
        plt.xlabel("Time (fs)", fontsize=18)
        plt.ylabel(r"HT ($\AA$)", fontsize=18)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.minorticks_on()
        plt.tight_layout()
        plt.savefig(path_script + f'{molecule}-contourf-density-broaden.pdf', dpi=300) 
      

        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=(6,3))
        gs = GridSpec(1,2, width_ratios=[3, 1])
        clr = ['#0072BD', 'orange', '#66BF2E', 'k', 'green', 'darkgreen', 'darkblue', 'orangered', '#D45500', 'crimson', 'dodgerblue', 'limegreen']
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        cbar = ax1.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(0.00,0.0001716,51), zorder=1)

        ax2.hist(IC[1], bins=50, range=[-3,3], density=True, weights=None, orientation='horizontal', color = 'orange', alpha = 0.9, zorder=10)
        ax2.hist(HT1, bins=50, range=[-3,3], density=True, histtype='step', orientation='horizontal', weights=None, color = 'blue', alpha = 0.9, zorder=11)
        ax2.hist(HT2, bins=50, range=[-3,3], density=True, histtype='step', orientation='horizontal', weights=None, color = 'limegreen', alpha = 0.9, zorder=12)

        ax1.scatter(t21, s21, s=0.6, c=clr[10], zorder=60)
        ax1.scatter(t12, s12, s=0.6, c=clr[11], zorder=61)
        
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
        ax1.set_xticks([0,25,50,75,100,125,150,175,200])
        ax1.set_ylim([-3.0,3.0])
        ax2.set_ylim([-3.0,3.0])
        #ax2.set_ticklabels([])
        ax2.set_yticks([-3,-2,-1,0,1,2,3])
        ax2.set_yticklabels(['' for _ in ax.get_yticks()])
        ax2.set_xlim([0.0,1.5])
        ax2.set_xticks([0.5,1.0,1.5])
        ax1.minorticks_on()
        ax2.minorticks_on()
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.04)
        fig.savefig(path_script + f'{molecule}-H-transfer-truedensity-colored-updated.png', dpi=500)
        fig.savefig(path_script + f'{molecule}-H-transfer-truedensity-colored-updated.pdf', dpi=500)

