import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from itertools import islice
import code_orthogonal_fit as ortho
import geom_param
import get_optim
from get_optim import read_xyz_as_np
#%matplotlib widget

def min_abs(t1, t2, t3):

    list_tor = [t1, t2, t3]
     # list of the 3 torsions
    tmin1 = min(abs(t1), abs(t2), abs(t3)) # min of absolute values of the 3 torsions
    tmin2 = 0
    if tmin1 in list_tor:
        tmin2 = tmin1
    else:
        tmin2 = -tmin1

    return tmin1, tmin2


def ave_tor(data):
#     """
#     
#     """
    
    torsion = []
    Htran = 0
    
    OH1 = geom_param.compute_bond(data, [4,13])
    OH2 = geom_param.compute_bond(data, [1,13])
        
    Htran = OH1-OH2

    tor1 = geom_param.compute_torsion(data,6,5,9,10)
    tor2 = geom_param.compute_torsion(data,6,5,9,11)
    tor3 = geom_param.compute_torsion(data,6,5,9,12)
                    
    return (Htran, min_abs(tor1, tor2, tor3)[0]) 


if __name__ == "__main__":

    au_to_fs = 0.02418884254

    path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ALL-QT-level/"
    pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S8/"


    os.chdir(path)
    list_data = ['complete-QT-traj-combined','movie-cl-updated', 'wigner_298.15K', 'complete-final-QT-geom', 'pulse-total-0.24']  #'complete-QT-traj-combined' already ran
    # break up the trajectory in individual xyz files
    for file in list_data:
        print(file)
        fdata = open(str(file) + ".xyz", 'r')   # use your data file

        num = 15 # number of atoms for this particular system # please change for your system
        
        before = []
        timstp = []
        #tim = []
        OH_avg = []
        tormin = []
        tormin2 = []

        os.system('mkdir tmp')
        for line in fdata:
            before.append(line)
            if len(before) > 2:
                before.pop(0)
            if "Time step: " in line:
                t = line.split()
                #print(t)
                t = t[2]                    # 7th column gives the actual timestep
                #print(t)
                #l = float(t) * au_to_fs    # timestep converted to fs
                timstp.append(t)           # timestep appended
                dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

                OH, tor = ave_tor(data)
                OH_avg.append(OH)
                tormin.append(tor)
                #tormin2.append(tor2)
                    
        fdata.close()
        os.system('rm -rf tmp')

        #print(CH_avg)
        #print(timstp)
    

        torf = open(pathscript + 'final-OH-newcoord-torsion-' + str(file) + '.dat','w')   
        for a,b,c in zip(timstp, OH_avg, tormin):        # timestep in au
            torf.write('{}\t{:.3f}\t{:.1f}\n'.format(a,b,c))

        torf.close()    

    os.chdir(pathscript)
