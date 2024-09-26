import numpy as np
import os
import glob
# load proper modules
import get_optim
import geom_param


#final_dir = "/proj/nhlist/users/x_prcha/CRITICAL-POINTS-ENONES/S0min-S1min-BLA-CENTRAL-ANGLE"
def calculate_distance(atom1_coord, atom2_coord):
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12


#O-index based
def calc_bla4m(filename):

    import get_optim
    import geom_param

    data = get_optim.read_xyz_as_np(filename)
    h1 = geom_param.compute_bond(data,[1,8])
    h2 = geom_param.compute_bond(data,[0,8])
    if h1 < h2:
        bla = calculate_distance(data[0],data[2]) - calculate_distance(data[2],data[3]) + calculate_distance(data[3],data[4]) - calculate_distance(data[4],data[1])
    else:
        bla = calculate_distance(data[4],data[1]) - calculate_distance(data[0],data[2]) + calculate_distance(data[2],data[3]) + calculate_distance(data[3],data[4])

    return bla

#O-index based
def calc_bla4ac(filename):

    import get_optim
    import geom_param

    data = get_optim.read_xyz_as_np(filename)
    h1 = geom_param.compute_bond(data,[4,13])
    h2 = geom_param.compute_bond(data,[1,13])
    if h1 < h2:
        bla = calculate_distance(data[1],data[0]) - calculate_distance(data[0],data[2]) + calculate_distance(data[2],data[3]) - calculate_distance(data[3],data[4])
    else:
        bla = calculate_distance(data[4],data[3]) - calculate_distance(data[3],data[2]) + calculate_distance(data[2],data[0]) - calculate_distance(data[1],data[0])

    return bla





#molecules = ['acrolein', 'cro', 'MVK', 'methacro']

crit_pts = ['S0-min', 'S1-min', 'T1-min', 'MECI-HTI-asym', 'MECI-HTI-sym', 'MECI-S2S1-tw', 'MECI-S1S0-tw', 'CIpy']

#O-index based
data_path = "/data/projects/Pratip_MA_vs_AcAc/CONVERGED-GEOMS/MA/CAS_108/"
current = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/geometrical-parameters/"

fac=open(current + 'MA.dat','w')
os.chdir(data_path)
for crit_pt in crit_pts:
    data = get_optim.read_xyz_as_np(str(crit_pt) + ".xyz")
    fac.write(crit_pt + '\n')
    fac.write('{}\t{:0.3f}\n'.format("O1-O2", geom_param.compute_bond(data,[0,1])))
    fac.write('{}\t{:0.3f}\n'.format("O1-H9", geom_param.compute_bond(data,[0,8])))
    fac.write('{}\t{:0.3f}\n'.format("O2-H9", geom_param.compute_bond(data,[1,8])))
    fac.write('{}\t{:0.3f}\n'.format("O1-C3", geom_param.compute_bond(data,[0,2])))
    fac.write('{}\t{:0.3f}\n'.format("C3-C4", geom_param.compute_bond(data,[2,3])))
    fac.write('{}\t{:0.3f}\n'.format("C4-C5", geom_param.compute_bond(data,[3,4])))
    fac.write('{}\t{:0.3f}\n'.format("C5-O2", geom_param.compute_bond(data,[4,1])))
    fac.write('{}\t{:0.1f}\n'.format("O1C3C4", geom_param.compute_angle(data,0,2,3)))
    fac.write('{}\t{:0.1f}\n'.format("C3C4C5", geom_param.compute_angle(data,2,3,4)))
    fac.write('{}\t{:0.1f}\n'.format("C4C5O2", geom_param.compute_angle(data,3,4,1)))
    fac.write('{}\t{:0.1f}\n'.format("C4pyr", geom_param.compute_pyramidalization(data,3,2,4,6)))
    fac.write('{}\t{:0.1f}\n'.format("O1C3C4C5", geom_param.compute_torsion(data,0,2,3,4)))
    fac.write('{}\t{:0.1f}\n'.format("C3C4C5O2", geom_param.compute_torsion(data,2,3,4,1)))         #2341
    fac.write('{}\t{:0.1f}\n'.format("SOA", geom_param.compute_angle(data,0,2,3) + geom_param.compute_angle(data,2,3,4) + geom_param.compute_angle(data,3,4,1)))
    fac.write('{}\t{:0.3f}\n'.format("H-coord", (geom_param.compute_bond(data,[1,8])-geom_param.compute_bond(data,[0,8]))))
    fac.write('{}\t{:0.3f}\n'.format("BLA", calc_bla4m(str(crit_pt) + ".xyz")))
    fac.write('\n')

os.chdir(current)
fac.close()


crit_pts = ['S0-min', 'S1-min', 'T1-min', 'MECI-HTI-asym-ecdown', 'MECI-HTI-asym-stag', 'MECI-HTI-asym-ecup', 'MECI-HTI-sym-ecdown', 'MECI-HTI-sym-ecup', 'MECI-S2S1-tw', 'MECI-S1S0-tw', 'CIpy']

data_path = "/data/projects/Pratip_MA_vs_AcAc/CONVERGED-GEOMS/AcAc/CAS_108/"

fac=open(current + 'ac-ac.dat','w')
os.chdir(data_path)
for crit_pt in crit_pts:
    data = get_optim.read_xyz_as_np(str(crit_pt) + ".xyz")
    fac.write(crit_pt + '\n')
    fac.write('{}\t{:0.3f}\n'.format("O2-O5", geom_param.compute_bond(data,[1,4])))
    fac.write('{}\t{:0.3f}\n'.format("O5-H14", geom_param.compute_bond(data,[4,13])))
    fac.write('{}\t{:0.3f}\n'.format("O2-H14", geom_param.compute_bond(data,[1,13])))
    fac.write('{}\t{:0.3f}\n'.format("O2-C1", geom_param.compute_bond(data,[1,0])))
    fac.write('{}\t{:0.3f}\n'.format("C1-C3", geom_param.compute_bond(data,[0,2])))
    fac.write('{}\t{:0.3f}\n'.format("C3-C4", geom_param.compute_bond(data,[2,3])))
    fac.write('{}\t{:0.3f}\n'.format("C4-O5", geom_param.compute_bond(data,[3,4])))
    fac.write('{}\t{:0.3f}\n'.format("C1-C6", geom_param.compute_bond(data,[0,5])))
    fac.write('{}\t{:0.3f}\n'.format("C4-C10", geom_param.compute_bond(data,[3,9])))
    fac.write('{}\t{:0.1f}\n'.format("O2C1C3", geom_param.compute_angle(data,1,0,2)))
    fac.write('{}\t{:0.1f}\n'.format("C1C3C4", geom_param.compute_angle(data,0,2,3)))
    fac.write('{}\t{:0.1f}\n'.format("C3C4O5", geom_param.compute_angle(data,2,3,4)))
    fac.write('{}\t{:0.1f}\n'.format("C3pyr", geom_param.compute_pyramidalization(data,2,0,3,14)))
    fac.write('{}\t{:0.1f}\n'.format("O2C1C3C4", geom_param.compute_torsion(data,1,0,2,3)))
    fac.write('{}\t{:0.1f}\n'.format("C1C3C4O5", geom_param.compute_torsion(data,0,2,3,4)))
    fac.write('{}\t{:0.1f}\n'.format("SOA", geom_param.compute_angle(data,1,0,2) + geom_param.compute_angle(data,0,2,3) + geom_param.compute_angle(data,2,3,4)))
    fac.write('{}\t{:0.3f}\n'.format("H-coord", (geom_param.compute_bond(data,[13,4])-geom_param.compute_bond(data,[1,13]))))
    fac.write('{}\t{:0.3f}\n'.format("BLA", calc_bla4ac(str(crit_pt) + ".xyz")))
    fac.write('{}\t{:0.3f}\n'.format("pyrM", geom_param.compute_pyramidalization(data,3,2,4,9)))
    fac.write('\n')

os.chdir(current)
fac.close()


