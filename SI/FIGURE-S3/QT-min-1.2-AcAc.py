import numpy as np
import os
import glob
# load proper modules
import get_optim
import geom_param
from itertools import islice
import matplotlib.pyplot as plt

#bond_pairs = [[0,2],[2,3],[3,4],[4,1],[1,8],[8,0]]

#CH_pairs = [[2,5],[3,6],[4,7]]
bond_pairs = [[0,2],[2,3],[3,4],[4,13],[13,1],[0,5],[3,9]]

CH_pairs = [[5,6],[5,7],[5,8],[9,10],[9,11],[9,12]]

#pair = [1,2]
#triple = [3,7,6]

#angle
#angles_trio = [[7,3,1],[3,1,2],[1,2,5]]
#pyra_four = [[3,1,6,7],[1,2,4,3],[2,0,1,5]]

au_to_fs = 0.02418884254

timstp = []
bond1 = []
bond2 = []
bond3 = []
bond4 = []
bond5 = []
bond6 = []
bond7 = []
Htran = []

avgCH = []
before = []

# break up the trajectory in individual xyz files
#fdata = open("movie.xyz", 'r')   # use your data file

num = 15 # number of atoms for this particular system # please change for your system

path = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/test-script-geom/code/distribution/"
pathscript = "/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S3/"

#path = "/storage/pratip/MA/S0-min-b3lyp-d3/WIGNER/wigner_298.15K_IC5000/"

# break up the trajectory in individual xyz files
fdata = open(path + 'total-movie-710fs-all-correct.xyz', 'r')   # use your data file

os.system('mkdir tmp')
for line in fdata:
    before.append(line)
    if len(before) > 2:
        before.pop(0)
    if "Time step: " in line:
         t = line.split()
         t = t[2]
         l = float(t) #* au_to_fs    # timestep converted to fs
         timstp.append(l)           # timestep appended
         dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
         f=open(dat,'w')
         f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
         f.close()
         data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

         bond = []
         for bond_pair in bond_pairs:
             bond.append(geom_param.compute_bond(data, bond_pair))        # internal coordinates calculated

         bond1.append(bond[0]) 
         bond2.append(bond[1])
         bond3.append(bond[2])
         bond4.append(bond[3])
         bond5.append(bond[4])
         bond6.append(bond[5])

         Htran.append(bond[3]-bond[4])
        
             # average of all CH bonds
         CH = []
         for CH_pair in CH_pairs:
             CH.append(geom_param.compute_bond(data, CH_pair))

         avgCH.append(np.average(CH))

    #if np.average(CH) < 1.2:

#        tim.append(i)
os.system('rm -rf tmp')

#bondf = open('final-bondpairs-wigner.dat','w')
#for i,j,k,l,m,n,o,p in zip(timstp, bond1, bond2, bond3, bond4, bond5, bond6, bond7):
#    bondf.write('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(int(i),j,k,l,m,n,o,p))

#bondf.close()

#avgf = open('final-CHavg-wigner.dat','w')
#for i,j in zip(timstp, avgCH):
#    avgf.write('{}\t{:.3f}\n'.format(int(i),j))

torf = open(pathscript + 'final-QT-120fs-298.15K-AcAc.dat','w')
for i,j,k in zip(timstp, Htran, avgCH):
    torf.write('{}\t{:.3f}\t{:.3f}\n'.format(i,j,k))

torf.close()


timstpk = []
bond1k = []
bond2k = []
bond3k = []
bond4k = []
bond5k = []
bond6k = []
bond7k = []
Htrank = []

avgCHk = []
beforek = []

# break up the trajectory in individual xyz files
fdata = open(path + 'complete-final-QT-geom.xyz', 'r')   # use your data file

os.system('mkdir tmp')
for line in fdata:
    beforek.append(line)
    if len(beforek) > 2:
        beforek.pop(0)
    if "Time step: " in line:
         t = line.split()
         t = t[2]
         l = float(t) #* au_to_fs    # timestep converted to fs
         timstpk.append(l)           # timestep appended
         dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
         f=open(dat,'w')
         f.write("".join(beforek) + "".join(islice(fdata, num)))  # writing the file
         f.close()
         data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

         bond = []
         for bond_pair in bond_pairs:
             bond.append(geom_param.compute_bond(data, bond_pair))        # internal coordinates calculated

         bond1k.append(bond[0])
         bond2k.append(bond[1])
         bond3k.append(bond[2])
         bond4k.append(bond[3])
         bond5k.append(bond[4])
         bond6k.append(bond[5])

         Htrank.append(bond[3]-bond[4])

             # average of all CH bonds
         CH = []
         for CH_pair in CH_pairs:
             CH.append(geom_param.compute_bond(data, CH_pair))

         avgCHk.append(np.average(CH))

    #if np.average(CH) < 1.2:

#        tim.append(i)
os.system('rm -rf tmp')

#bondf = open('final-bondpairs-wigner.dat','w')
#for i,j,k,l,m,n,o,p in zip(timstp, bond1, bond2, bond3, bond4, bond5, bond6, bond7):
#    bondf.write('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(int(i),j,k,l,m,n,o,p))

#bondf.close()

#avgf = open('final-CHavg-wigner.dat','w')
#for i,j in zip(timstp, avgCH):
#    avgf.write('{}\t{:.3f}\n'.format(int(i),j))

torf = open(pathscript + 'final-QT-abs-298.15K-AcAc.dat','w')
for i,j,k in zip(timstpk, Htrank, avgCHk):
    torf.write('{}\t{:.3f}\t{:.3f}\n'.format(i,j,k))

torf.close()




timstps = []
bond1s = []
bond2s = []
bond3s = []
bond4s = []
bond5s = []
bond6s = []
bond7s = []
Htrans = []

avgCHs = []
befores = []

# break up the trajectory in individual xyz files
#fdata = open("movie.xyz", 'r')   # use your data file

num = 15 # number of atoms for t

# break up the trajectory in individual xyz files
fdata = open(path + 'pulse-total-0.24.xyz', 'r')   # use your data file

os.system('mkdir tmp')
for line in fdata:
    befores.append(line)
    if len(befores) > 2:
        befores.pop(0)
    if "Time step: " in line:
         t = line.split()
         t = t[2]
         l = float(t) #* au_to_fs    # timestep converted to fs
         timstps.append(l)           # timestep appended
         dat = 'tmp/' + str(t) + '.xyz'       # coordinates saved
         f=open(dat,'w')
         f.write("".join(before) + "".join(islice(fdata, num)))  # writing the file
         f.close()
         data = get_optim.read_xyz_as_np('tmp/' + str(t) + '.xyz')        # coordinates loaded

         bond = []
         for bond_pair in bond_pairs:
             bond.append(geom_param.compute_bond(data, bond_pair))        # internal coordinates calculated

         bond1s.append(bond[0]) 
         bond2s.append(bond[1])
         bond3s.append(bond[2])
         bond4s.append(bond[3])
         bond5s.append(bond[4])
         bond6s.append(bond[5])

         Htrans.append(bond[3]-bond[4])
        
             # average of all CH bonds
         CH = []
         for CH_pair in CH_pairs:
             CH.append(geom_param.compute_bond(data, CH_pair))

         avgCHs.append(np.average(CH))

    #if np.average(CH) < 1.2:

#        tim.append(i)
os.system('rm -rf tmp')

#bondf = open('final-bondpairs-wigner.dat','w')
#for i,j,k,l,m,n,o,p in zip(timstp, bond1, bond2, bond3, bond4, bond5, bond6, bond7):
#    bondf.write('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(int(i),j,k,l,m,n,o,p))

#bondf.close()

#avgf = open('final-CHavg-wigner.dat','w')
#for i,j in zip(timstp, avgCH):
#    avgf.write('{}\t{:.3f}\n'.format(int(i),j))

torfs = open(pathscript + 'final-QT-pulse-298.15K-AcAc.dat','w')
for i,j,k in zip(timstps, Htrans, avgCHs):
    torfs.write('{}\t{:.3f}\t{:.3f}\n'.format(i,j,k))

torfs.close()

#fig, ax = plt.subplots(1,1, figsize=(4,3))
#plt.axvline(x=0.625, linestyle = '--', color = 'black')
#plt.hist(Htran, bins=30, range=[-2.5,2.5], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'Wigner 298.15K')
    #plt.axvline(x=1.675, color = 'gray', linestyle='--', label='0.5*ZPE')
#plt.xlabel(r'H-transfer coordinate ($\AA$)')
#plt.ylabel('Normalized Distribution')
#plt.xlim(-2.0,2.0)
#plt.ylim(0,2.5)
  #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
#plt.xticks([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0])
#plt.legend(loc ="upper right", frameon=False)
#plt.tight_layout()
#ax.minorticks_on()
#plt.savefig('MA-Htran-298.15K-QTv1.png',dpi=500)

#fig, ax = plt.subplots(1,1, figsize=(4,3))
#plt.hist(avgCH, bins=40, range=[0.0,2.0], density=True, weights=None, color = 'orange', alpha = 0.9, label=r'Wigner 298.15K')
#plt.axvline(x=1.092, linestyle = '--', color = 'black')
#plt.axvline(x=30, color = 'black')
#plt.axvline(x=-30, color = 'black')
#plt.xticks([-90,-70,-50,-30,-10,10,30,50,70,90])
#plt.xlim(0,2.0)
#plt.ylim(0,10.0)
#plt.xlabel(r'Average CH distance ($\AA$)')
#plt.ylabel('Normalized distribution')
#plt.legend(loc ="upper right", frameon=False)
#plt.tight_layout()
#ax.minorticks_on()
#plt.savefig('MA-avdCH-298.15K-QTv1.png',dpi=500)


timstpp = []
bond11 = []
bond22 = []
bond33 = []
bond44 = []
bond55 = []
bond66 = []
bond77 = []
Htran1 = []

avgCH1 = []
before1 = []

num = 15 # number of atoms for this particular system # please change for your system

path2 = "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ACETYLACETONE-ground-state-tests/staggered/B3LYP-D3/HESSIAN-bin/Wigner_298.15K_5000/wigner_298.15K_IC5000/"

for i in range(0,5000,1):
    i = str(i).zfill(4)
    timstpp.append(i)
#
    data = get_optim.read_xyz_as_np(path2 + 'x' + str(i) + '.xyz')        # coordinates loaded
    bond = []
    for bond_pair in bond_pairs:
        bond.append(geom_param.compute_bond(data, bond_pair))        # internal coordinates calculated

        #print(bond)
    bond11.append(bond[0]) 
    bond22.append(bond[1])
    bond33.append(bond[2])
    bond44.append(bond[3])
    bond55.append(bond[4])
    bond66.append(bond[5])

    Htran1.append(bond[3]-bond[4])
        
    # average of all CH bonds
    CH = []
    for CH_pair in CH_pairs:
        CH.append(geom_param.compute_bond(data, CH_pair))

    avgCH1.append(np.average(CH))

   #if np.average(CH) < 1.2:



os.system('rm -rf tmp')

#bondf = open('final-bondpairs-wigner.dat','w')
#for i,j,k,l,m,n,o,p in zip(timstp, bond1, bond2, bond3, bond4, bond5, bond6, bond7):
#    bondf.write('{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(int(i),j,k,l,m,n,o,p))

#bondf.close()

#avgf = open('final-CHavg-wigner.dat','w')
#for i,j in zip(timstp, avgCH):
#    avgf.write('{}\t{:.3f}\n'.format(int(i),j))

torff = open(pathscript + 'final-wigner-298.15K-AcAc.dat','w')
for i,j,k in zip(timstpp, Htran1, avgCH1):
    torff.write('{}\t{:.3f}\t{:.3f}\n'.format(i,j,k))

torff.close()

#data = np.genfromtxt('../complete-QT-version2-298.15K.dat', dtype='float')

fig, ax = plt.subplots(1,1, figsize=(4,3))
#plt.axvline(x=30, color = 'black')
plt.hist(Htran, bins=30, range=[-2.5,2.5], density=True, weights=None, color = 'orange', alpha = 0.9, label='QT(all)')
plt.hist(Htrank, bins=30, range=[-2.5,2.5], density=True, weights=None, histtype = 'step', color = 'blue', alpha = 0.9, label='QT(abs)')
plt.hist(Htrans, bins=30, range=[-2.5,2.5], density=True, weights=None, histtype = 'step', color = 'limegreen', alpha = 0.9, label='QT(pulse)')
plt.hist(Htran1, bins=30, range=[-2.5,2.5], density=True, weights=None, histtype='step', color = 'darkred', alpha = 0.9, label='Wigner')
plt.axvline(x=-0.624, linestyle = '--', color = 'gray', label='FC')
plt.xlabel(r'HT ($\AA$)', fontsize=14)
plt.ylabel('Norm. dist.', fontsize=14)
plt.xlim(-2.0,2.0)
plt.ylim(0,2.5)
plt.yticks(fontsize=14)
plt.xticks([-2.0,-1.0,0.0,1.0,2.0],fontsize=14)
plt.legend(loc ="upper right", frameon=False, fontsize=8)
plt.tight_layout()
ax.minorticks_on()
plt.savefig('AcAc-Htran-298.15K_QTWIG_1890.png',dpi=500)
plt.savefig('AcAc-Htran-298.15K_QTWIG_1890.pdf',dpi=500)

fig, ax = plt.subplots(1,1, figsize=(4,3))
plt.hist(avgCH, bins=40, range=[0.0,2.0], density=True, weights=None, color = 'orange', alpha = 0.9, label='QT(all)')
plt.hist(avgCHk, bins=40, range=[0.0,2.0], density=True, weights=None, histtype = 'step', color = 'blue', alpha = 0.9, label='QT(abs)')
plt.hist(avgCHs, bins=40, range=[0.0,2.0], density=True, weights=None, histtype = 'step', color = 'limegreen', alpha = 0.9, label='QT(pulse)')
plt.hist(avgCH1, bins=40, range=[0.0,2.0], density=True, weights=None, histtype='step', color = 'darkred', alpha = 0.9, label='Wigner')
#plt.axvline(x=1.092, linestyle = '--', label='FC', color='gray')
plt.axvline(x=30, color = 'black')
plt.axvline(x=-30, color = 'black')
#plt.xticks([-90,-70,-50,-30,-10,10,30,50,70,90])
plt.yticks(fontsize=14)
plt.xlim(0.5,2.5)
plt.xticks([0.5,1.0,1.5,2.0,2.5],fontsize=14)
plt.ylim(0,15.0)
plt.xlabel(r'Average CH distance ($\AA$)', fontsize=14)
plt.ylabel('Norm. dist.', fontsize=14)
plt.legend(loc ="upper right", frameon=False, fontsize=14)
plt.tight_layout()
ax.minorticks_on()
plt.savefig('AcAc-avdCH-298.15K_QTWIG_1890.png',dpi=500)
plt.savefig('AcAc-avdCH-298.15K_QTWIG_1890.pdf',dpi=500)
