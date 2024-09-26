import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from itertools import islice
import code_orthogonal_fit as ortho
import geom_param
import get_optim
from get_optim import read_xyzs_as_np
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

#data = read_xyz_as_np('x3966.xyz')

# calculate average of all CH bonds
def average_CH(data, CH_pairs):
    """
    Takes a list of all CH bond pairs and returns the average of all of them
    """
    CH = []
    for CH_pair in CH_pairs:
        CH.append(geom_param.compute_bond(data, CH_pair))

    return np.average(CH)



if __name__ == "__main__":

    au_to_fs = 0.02418884254

    path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ALL-QT-level/"
    pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S7/"

    os.chdir(path)
    list_data = ['wigner_298.15K'] #, 'complete-QT-traj-combined'] #ran it already no need to run again
    #list_data = ['complete-final-QT-geom'] #, 'complete-QT-traj-combined'] #ran it already no need to run again

    for file in list_data:
        #fdata = open(str(file) + ".xyz", 'r')   # use your data file
        data = get_optim.read_xyzs_as_np(str(file) + ".xyz")
        num = 15 # number of atoms for this particular system # please change for your system
    
        CH_pairs = [[5,6],[5,7],[5,8],[9,10],[9,11],[9,12],[2,14]]
        
        before = []
        timstp = []
        tim = []
        t1 = 0
        t2 = 0
        t3 = 0
        t4 = 0
        t5 = 0
        t6 = 0
        tor1 = []
        tor2 = []
        tor3 = []
        tor4 = []
        tor5 = []
        tor6 = []
        tmin = []
        tmin2 = []
        tmin_new = []
        tmin2_new = []



        for i in range(len(data[:])):    
                # normal
                #geometry = ortho.coord(data,[1,0,3,4])
                #normal_to_ortho_fit_plane = np.squeeze(np.asarray(ortho.orthogonal_plane_fit(geometry)[1]))
            if average_CH(data[i], CH_pairs) < 1.2:
                normal_to_ortho_fit_plane = ortho.compute_torsion_orthogonal_plane(data[i], [5,6,0], [1,0,3,4])[2]

                timstp.append(i)  # frame no

                t1 = ortho.compute_torsion_orthogonal_plane(data[i], [5,6,0], [1,0,3,4])  # torsion of one methyl group to the chelate ring
                #print(t1[0])
                #print(t1[1])
                tor1.append(t1[0])
                
                
                t2 = ortho.compute_torsion_orthogonal_plane(data[i], [5,7,0], [1,0,3,4])  # torsion of one methyl group to the chelate ring
                tor2.append(t2[0])

                t3 = ortho.compute_torsion_orthogonal_plane(data[i], [5,8,0], [1,0,3,4])  # torsion of one methyl group to the chelate ring
                tor3.append(t3[0])

                tmin.append(min_abs(t1[0], t2[0], t3[0])[0])       # list of absolute values
                tmin_new.append(min_abs(t1[0], t2[0], t3[0])[1])   # list of actual values
 
                t4 = ortho.compute_torsion_orthogonal_plane(data[i], [9,3,10], [1,0,3,4])  # torsion of one methyl group to the chelate ring
                tor4.append(t4[0])

                t5 = ortho.compute_torsion_orthogonal_plane(data[i], [9,3,11], [1,0,3,4])  # torsion of one methyl group to the chelate ring
                tor5.append(t5[0])

                t6 = ortho.compute_torsion_orthogonal_plane(data[i], [9,3,12], [1,0,3,4])  # torsion of one methyl group to the chelate ring
                tor6.append(t6[0])
    
                tmin2.append(min_abs(t4[0], t5[0], t6[0])[0])       # list of absolute values
                tmin2_new.append(min_abs(t4[0], t5[0], t6[0])[1])   # list of actual values



        torf = open(pathscript + 'final-torsion-ortho-plane-QT-all-final-' + str(file) + '1.2.dat','w')
        for a,b,c,d,e,f,g,h,i,j,k in zip(timstp, tor1, tor2, tor3, tor4, tor5, tor6, tmin, tmin_new, tmin2, tmin2_new):
            torf.write('{:.2f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\n'.format(a,b,c,d,e,f,g,h,i,j,k))

        torf.close()    
