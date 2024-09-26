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


# # calculate average of 3 CH bonds depending on O-H distance and also calculate minimum of the corresponding torsion O-C-C-H
def ave_min(data):
#     """
#     Takes xyz coordinates and returns the average of CH bonds of lowest frequency methyl group
#     """
    CH = []
    CH_pair1 = [[5,6],[5,7],[5,8]]
    torsion = []
    #torsion_pair1 = [[1,0,5,6],[1,0,5,7],[1,0,5,8]]
    OH1 = geom_param.compute_bond(data, [4,13])
    CH_pair2 = [[9,10],[9,11],[9,12]]
    OH2 = geom_param.compute_bond(data, [1,13])
    #torsion_pair2 = [[4,3,9,10],[4,3,9,11],[4,3,9,12]]

    if OH1 < OH2: 
        for CH_pair in CH_pair1:
            CH.append(geom_param.compute_bond(data, CH_pair))
      
        #torsion.append(abs(geom_param.compute_torsion(data,1,0,5,6)))
        #torsion.append(abs(geom_param.compute_torsion(data,1,0,5,7)))
        #torsion.append(abs(geom_param.compute_torsion(data,1,0,5,8)))
        tor1 = geom_param.compute_torsion(data,1,0,5,6)
        tor2 = geom_param.compute_torsion(data,1,0,5,7)
        tor3 = geom_param.compute_torsion(data,1,0,5,8)
                    
        return (np.average(CH), min_abs(tor1, tor2, tor3)[0]) 
    else:

        for CH_pair in CH_pair2:
            CH.append(geom_param.compute_bond(data, CH_pair))
        
        #torsion.append(abs(geom_param.compute_torsion(data,4,3,9,10)))
        #torsion.append(abs(geom_param.compute_torsion(data,4,3,9,11)))
        #torsion.append(abs(geom_param.compute_torsion(data,4,3,9,12)))
                    
        tor1 = geom_param.compute_torsion(data,4,3,9,10)
        tor2 = geom_param.compute_torsion(data,4,3,9,11)
        tor3 = geom_param.compute_torsion(data,4,3,9,12)
                    
        return (np.average(CH), min_abs(tor1, tor2, tor3)[0]) 



if __name__ == "__main__":

    au_to_fs = 0.02418884254

    path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ALL-QT-level/"
    pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S9/"

    os.chdir(path)
    list_data = ['movie-cl-updated', 'wigner_298.15K', 'complete-final-QT-geom', 'pulse-total-0.24', 'complete-QT-traj-combined'] #ran it already no need to run again
    # break up the trajectory in individual xyz files
    for file in list_data:
        fdata = open(str(file) + ".xyz", 'r')   # use your data file

        num = 15 # number of atoms for this particular system # please change for your system
        
        before = []
        timstp = []
        #tim = []
        CH_avg = []
        tormin = []

        os.system('mkdir tmp2')
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
                dat = 'tmp2/' + str(t) + '.xyz'       # coordinates saved
                f=open(dat,'w')
                f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
                f.close()
                data = read_xyz_as_np('tmp2/' + str(t) + '.xyz')        # coordinates loaded

                #print(geom_param.compute_bond(data, [4,13]))
                #print(geom_param.compute_bond(data, [1,13]))
                CHavg, tor = ave_min(data)
                CH_avg.append(CHavg)
                tormin.append(tor)
                    
               
        fdata.close()
        os.system('rm -rf tmp2')

        #print(CH_avg)
        #print(timstp)
    

        torf = open(pathscript + 'final-CH-torsion-' + str(file) + '.dat','w')   
        for a,b,c in zip(timstp, CH_avg, tormin):        # timestep in au
            torf.write('{}\t{:.3f}\t{:.1f}\n'.format(a,b,c))

        torf.close()    

    os.chdir(pathscript)
