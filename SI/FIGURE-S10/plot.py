import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    clr = ['black', 'limegreen', 'orange','blue','cyan']

    fig, ax = plt.subplots(2,3, figsize=(10,6), sharex=True, sharey=True)

    # CAS(14,12)
    for nstate in range(3,6):
        data = np.genfromtxt('/data/projects/Pratip_MA_vs_AcAc/AcAc-data/state_averaging_test_IC110/XMS-14-12-sssr-imag/SA' + str(nstate) + '/final-result-au-nm.dat', dtype='float')
        for i in range(1,nstate+1):
            ax[0,nstate-3].plot(data[:,0], (data[:,i]-data[1,1])*27.2114, c=clr[i-1])
            ax[0,nstate-3].set_ylim([-0.2,7.5])
            ax[0,nstate-3].set_xlim([33,50])
            ax[0,nstate-3].set_xticks([35,40,45,50])
            ax[0,nstate-3].set_yticks([0,1,2,3,4,5,6,7])
            ax[0,0].set_ylabel('Energy(eV)',fontsize=14)
            ax[0,nstate-3].minorticks_on()
            ax[0,nstate-3].tick_params(axis='both', which='major', labelsize=14)
            ax[0,nstate-3].tick_params(axis='both', which='minor', labelsize=14)

    # CAS(10,8)
    for nstate in range(3,6):
        data = np.genfromtxt('/data/projects/Pratip_MA_vs_AcAc/AcAc-data/state_averaging_test_IC110/SA' + str(nstate) + '-10-8/final-result-au-nm.dat', dtype='float')
        for i in range(1,nstate+1):
            ax[1,nstate-3].plot(data[:,0], (data[:,i]-data[1,1])*27.2114, c=clr[i-1], label=f"S$_{i-1}$")          # 33 fs is the reference
            ax[1,nstate-3].set_ylim([-0.2,7.5])
            ax[1,nstate-3].set_xlim([33,50])
            ax[1,nstate-3].set_xticks([35,40,45,50])
            ax[1,nstate-3].set_yticks([0,1,2,3,4,5,6,7])
            ax[1,nstate-3].set_xlabel('Time (fs)',fontsize=14)
            ax[1,0].set_ylabel('Energy(eV)',fontsize=14)
            ax[1,nstate-3].minorticks_on()
            ax[1,nstate-3].tick_params(axis='both', which='major', labelsize=14)
            ax[1,nstate-3].tick_params(axis='both', which='minor', labelsize=14)
    
    ax[1,2].legend(frameon=False)
    ax[1,0].axvline(x=36.5, linestyle='--', color='gray')
    fig.subplots_adjust(hspace=0.05,wspace=0.05)
    #fig.legend()
    fig.tight_layout()
    
    fig.savefig('TRAJ110-testing-updated.pdf', dpi=500)

