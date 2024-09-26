import os
import numpy as np

if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    Ntraj = [257,272]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
              ]
    Natoms = [9, 15]

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/Population/'
    for molecule, ntraj, run_path, numatoms  in zip(molecules, Ntraj, runpaths, Natoms):

        os.chdir(run_path)

        for i in range(1,ntraj):
            i = "{:05d}".format(i)

            os.mkdir(path_script + f'{molecule}_pop_data/' + str(i))
            os.system('echo $(pwd)')	
            os.chdir('TRAJ_' + str(i))
            if os.path.exists('DONT_ANALYZE'):
                os.chdir(run_path)
            else:
                os.chdir('output_data')

                #os.system('echo $(pwd)')
                s0 = [] 
                s1 = []
                s2 = []
                tot = []
                time = []

                # reading classical population
                fpop=open('coeff_class_MCH.out', 'r')
                for iline, line in enumerate(fpop):
                    if iline > 2:
                        s0.append(float(line.split()[2]))
                        s1.append(float(line.split()[3]))
                        s2.append(float(line.split()[4]))
                        tot.append(float(line.split()[2])+float(line.split()[3])+float(line.split()[4]))
                        time.append(float(line.split()[0]))

                fpop.close()

                os.chdir(path_script + f'{molecule}_pop_data/' + str(i))
                fwpop=open('N.dat', 'w+')
                for i,j,k,l,m in zip(time, s0, s1, s2, tot): 
                    fwpop.write('{:.2f}\t{:.10f}\t{:.10f}\t{:.10f}\t{:.10f}\n'.format(i,j,k,l,m))    

                fwpop.close()
    
                os.chdir(run_path)

