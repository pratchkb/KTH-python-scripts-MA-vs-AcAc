import os

pathMA = '/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/'
pathAcAc = '/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/'
path = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/BLA_PLATEAU/'

os.chdir(pathMA)
for i in range(1,257):
    i = "{:05d}".format(i)
    os.chdir('TRAJ_'+str(i))
    if os.path.exists('DONT_ANALYZE'):
        #print(str(i))
        os.chdir(pathMA)
    else:
        if os.path.exists('BACK_HOP'):
        with open('output.lis', 'r') as f:
            for line in f:
                if "Surface Hop/Pointer State Switch: new state=    3 old state=    2" in line:
                    os.system('touch BACKHOP')

        os.chdir(pathMA)

#os.chdir(pathAcAc)

#for i in range(1,271):
#    i = "{:05d}".format(i)
#    os.chdir('TRAJ_'+str(i))
#    if os.path.exists('DONT_ANALYZE'):
#        #print(str(i))
#        os.chdir(pathAcAc)
#    else:
#        with open('output.lis', 'r') as f:
#            for line in f:
#                if "Surface Hop/Pointer State Switch: new state=    3 old state=    2" in line:
#                    os.system('touch BACKHOP')
#
#        os.chdir(pathAcAc)
#
os.chdir(path)

