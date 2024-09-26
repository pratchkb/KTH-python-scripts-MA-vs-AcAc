import numpy as np
import os
import os.path
import glob
# load proper modules
import get_optim
import geom_param
from itertools import islice
from itertools import chain
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

def collect_geom_param(filename, natoms, pyrc):

    num = natoms # number of atoms for this particular system # please change for your system

    time = []
    pyr = []

    data = get_optim.read_xyzs_as_np(filename)
    time = []
    htran = []
    for i in range(len(data[:])):
        time.append(i*0.5)
        pyr.append(geom_param.compute_pyramidalization(data[i],*pyrc))

    return time, pyr, data

def collect_geom_param_time(filename, natoms, pyrc, time_chunk):

    num = natoms # number of atoms for this particular system # please change for your system

    time = []
    pyr = []
     
    data = get_optim.read_xyzs_as_np(filename)

    addedtime = 200 - time_chunk[-1]
    for i in np.arange(time_chunk[0], time_chunk[-1], 0.5):
        time.append(i + addedtime)
        pyr.append(geom_param.compute_pyramidalization(data[int(i*2)],*pyrc))

    return time, pyr, data

def calculate_s1s0hops(filename, natoms):

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

            if int(st) == 1:         # S1 or S0 state
                l = float(ti) #* au_to_fs    # timestep converted to fs
                timstps0.append(l)           # timestep appended

    return timstps2, timstps1, timstps0


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



def calculate_param(filename, atom_order):
    """
    calculate all the BLAs along the trajectories depending on the numbering of atoms decided 
    """


    # break up the trajectory in individual xyz files
    fdata = open(filename, 'r')   # use your data file

    num = 15 # number of atoms for this particular system # please change for your system

    before = []
    timstp = []

    #h1 = []
    #h2 = []
    bla4 = []
    blaCO = []
    blaCC = []
    Htran = []

    os.system('mkdir tmp2')
    for line in fdata:
        before.append(line)
        if len(before) > 2:
            before.pop(0)
        if "t=  " in line:
             t = line.split()
             t = t[1]
             l = float(t) #* au_to_fs    # timestep converted to fs
             timstp.append(l)           # timestep appended
             dat = 'tmp2/' + str(t) + '.xyz'       # coordinates saved
             f=open(dat,'w')
             f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
             f.close()
             data = get_optim.read_xyz_as_np('tmp2/' + str(t) + '.xyz')        # coordinates loaded

             bla4.append(calculate_bla4(data,atom_order)[0])
             blaCO.append(calculate_bla4(data,atom_order)[1])
             blaCC.append(calculate_bla4(data,atom_order)[2])
             Htran.append(geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))

    fdata.close()
    
    os.system('rm -r tmp2')

    return timstp, bla4, blaCO, blaCC

#############################################################################################

if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    #Ntraj = [10,10]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    fig, ax = plt.subplots(2,5,figsize=(17,7))

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/PYRAMIDALIZATION-Central-derivative/'

    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        if molecule == "MA":
            pyrs = [3,2,4,6]

        if molecule == "AcAc":
            pyrs = [2,0,3,14]

        traj_hop = calculate_s1s0hops(run_path + "crossing-s1-s0.xyz", numatoms)

        os.chdir(run_path)
        bes10 = []
        besno = []

        avgbes10 = []
        avgbesno = []
       
        collects21 = []
        collectnos21 = []
        avgcollects21 = []
        avgcollectnos21 = []


        for i in range(1,ntraj):
            i = "{:05d}".format(i)
            os.chdir(run_path + 'TRAJ_'+str(i))
            if os.path.exists('DONT_ANALYZE'):

                os.chdir(run_path)

            else:
        
                if i in traj_hop:  
                    os.chdir(run_path + 'TRAJ_'+str(i))
                    #densities = collect_geom_param("output.xyz", numatoms, pyrs)
                    time_divide = calculate_param_statewise("output.xyz", numatoms)
                    s1 = collect_geom_param_time("output.xyz", numatoms, pyrs, time_divide[1])
                    #print(s1[0],s1[1])

                    pos = np.where(np.abs(np.diff(s1[0])) > 0.5)[0] + 1
                    time = np.insert(s1[0], pos, np.nan)
                    pyr = np.insert(s1[1], pos, np.nan)

                    ax[0,0].plot(time, pyr, color = 'orange', linewidth=0.4, zorder=0)
                    ax[0,0].set_ylabel(r"PyrC ($^{\circ}$)")
                    ax[0,0].set_ylim(-90.0,90.0)
                    ax[0,0].set_xlabel(r"Time (fs)")
                    ax[0,0].set_xlim(0,200)
                    ax[0,0].minorticks_on()

                    dx = 0.5   # timestep 
   
                    deriv = np.gradient(pyr, dx)
                    
                    ax[0,1].plot(time, deriv, color = 'orange', linewidth=0.4, zorder=0)
                    ax[0,1].set_ylabel(r"$\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
                    ax[0,1].set_ylim(-10.0,10.0)
                    ax[0,1].set_xlabel(r"Time (fs)")
                    ax[0,1].set_xlim(0,200)
                    ax[0,1].minorticks_on()

                    for elem in deriv[-10:]:
                        bes10.append(elem)    #5 fs before reaching s1-s0 hop 
                        
                    avgbes10.append(np.average(deriv[-10:]))
                
                    for elem in deriv[:10]:
                        collects21.append(elem) #5 fs after s2-s1 hop 
                    
                    avgcollects21.append(np.average(deriv[:10]))
           
                    os.chdir(run_path)

                else:
                    os.chdir(run_path + 'TRAJ_'+str(i))
                    #densities = collect_geom_param("output.xyz", numatoms, pyrs)
                    time_divide = calculate_param_statewise("output.xyz", numatoms)
                    s1 = collect_geom_param_time("output.xyz", numatoms, pyrs, time_divide[1])

                    pos = np.where(np.abs(np.diff(s1[0])) > 0.5)[0] + 1
                    time = np.insert(s1[0], pos, np.nan)
                    pyr = np.insert(s1[1], pos, np.nan)

                    ax[1,0].plot(time, pyr, color = 'orange', linewidth=0.4, zorder=0)
                    ax[1,0].set_ylabel(r"PyrC ($^{\circ}$)")
                    ax[1,0].set_ylim(-90.0,90.0)
                    ax[1,0].set_xlabel(r"Time (fs)")
                    ax[1,0].set_xlim(0,200)
                    ax[1,0].minorticks_on()


                    dx = 0.5
     
                    deriv2 = np.gradient(pyr, dx)

                    ax[1,1].plot(time, deriv2, color = 'orange', linewidth=0.4, zorder=0)
                    ax[1,1].set_ylabel(r"$\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
                    ax[1,1].set_ylim(-10.0,10.0)
                    ax[1,1].set_xlabel(r"Time (fs)")
                    ax[1,1].set_xlim(0,200)
                    ax[1,1].minorticks_on()
                    
                    #besno.append(deriv2[-10:])
                    for elem in deriv2[-10:]:
                        besno.append(elem)
                   
                    avgbesno.append(np.average(deriv2[-10:]))

                    for elem in deriv2[:10]:
                        collectnos21.append(elem)

                    avgcollectnos21.append(np.average(deriv2[:10]))

                    os.chdir(run_path)

        #print(collects21)
        #print(collectnos21)
        
        ax[0,2].hist(bes10, color='orange', density=True, range=[-10.0,10.0], bins = 40)
        ax[0,2].set_ylabel("Norm. dist.")
        #ax[0,1].set_ylim(-10.0,10.0)
        ax[0,2].set_xlabel(r"$\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[0,2].set_xlim(-10.0,10.0)
        ax[0,2].minorticks_on()

        ax[0,3].hist(collects21, color='orange', density=True, range=[-10.0,10.0], bins = 40)
        ax[0,3].set_ylabel("Norm. dist.")
        #ax[0,1].set_ylim(-10.0,10.0)
        ax[0,3].set_xlabel(r"$\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[0,3].set_xlim(-10.0,10.0)
        ax[0,3].minorticks_on()

        ax[0,4].scatter(avgcollects21, avgbes10, color='orange')
        ax[0,4].set_xlabel(r"After S$_2$/S$_1$ $\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[0,4].set_ylim(-10.0,10.0)
        ax[0,4].set_ylabel(r"Before S$_1$/S$_0$ $\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[0,4].set_xlim(-10.0,10.0)
        ax[0,4].minorticks_on()

        ax[1,2].hist(besno, color='orange', density=True, range=[-10.0,10.0], bins = 40)
        ax[1,2].set_ylabel("Norm. dist.")
        #ax[0,1].set_ylim(-10.0,10.0)
        ax[1,2].set_xlabel(r"$\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[1,2].set_xlim(-10.0,10.0)
        ax[1,2].minorticks_on()

        ax[1,3].hist(collectnos21, color='orange', density=True, range=[-10.0,10.0], bins = 40)
        ax[1,3].set_ylabel("Norm. dist.")
        #ax[0,1].set_ylim(-10.0,10.0)
        ax[1,3].set_xlabel(r"$\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[1,3].set_xlim(-10.0,10.0)
        ax[1,3].minorticks_on()

        ax[1,4].scatter(avgcollectnos21, avgbesno, color='orange')
        ax[1,4].set_xlabel(r"After S$_2$/S$_1$ $\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[1,4].set_ylim(-10.0,10.0)
        ax[1,4].set_ylabel(r"S$_1$ last 5 fs $\frac{d}{dt}$PyrC ($^{\circ}$/fs)")
        ax[1,4].set_xlim(-10.0,10.0)
        ax[1,4].minorticks_on()

        plt.tight_layout()
        plt.savefig(path_script + f'{molecule}-hopped_traj_pyr_s10sync.png', dpi=500)

