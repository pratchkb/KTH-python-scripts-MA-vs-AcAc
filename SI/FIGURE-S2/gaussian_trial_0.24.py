import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.0

# Speed of light in atomic units
h_evs = 4.13566733e-15 #eV * s
c_m_per_s = 299792458.
c_au = 137.035999679
au_per_fs = 41.341374575751
ev_per_au = 27.211386245988 

def read_bagel_xmscaspt2_results(nsamples,nstates):
  
    files = []
    results_dict = {}
    results_dict_npi = {}
    results_dict_pipi = {}
    count = 0

    # Iterate through all frameids sampled for absorption within current BOMD frame
    for isample in range(1,nsamples,1):          # very specific right now needs to change

    # Grab excited state energies and oscillator strengths
        results_file = str(isample) + '/input' + str(isample) + '.out'  #os.getcwd() + '/x{:04d}/{}'.format(isample,outname)
        #print(results_file)

        results = []
        results_npi = []
        results_pipi = []
        energies = []
        energy_npi = []
        energy_pipi = []
        oscs = []
        oscs_npi = []
        oscs_pipi = []

        result_S0 = []
        result_S1 = []
        result_S2 = []
       
        with open(results_file, 'r') as f_in:
            for line in f_in:

                if "MS-CASPT2 energy : state  0" in line:
                    state1, energy1 = int(line.split()[5]),float(line.split()[6])
                    result_S0.append(energy1)

                if "MS-CASPT2 energy : state  1" in line:
                    state2, energy2 = int(line.split()[5]),float(line.split()[6])
                    result_S1.append(energy2)

                if "MS-CASPT2 energy : state  2" in line:
                    state3, energy3 = int(line.split()[5]),float(line.split()[6])
                    result_S2.append(energy3)

                if ' * Transition    1 - 0 :' in line:
                    osc1 = float(next(f_in, '').strip().split()[4])

                if ' * Transition    2 - 0 :' in line:
                    osc2 = float(next(f_in, '').strip().split()[4])
                    
                    energy_S1 = (energy2 - energy1)*ev_per_au
                    energy_S2 = (energy3 - energy1)*ev_per_au
                  
                    energy_npi.append(energy_S1)
                    oscs_npi.append(osc1)
                    energy_pipi.append(energy_S2)
                    oscs_pipi.append(osc2)
                    energies.append(energy_S1)
                    oscs.append(osc1)
                    energies.append(energy_S2)
                    oscs.append(osc2)
        

        for energy, osc in zip(energies, oscs):
            results.append([energy, osc])
            results_dict[isample] = results 

        for energy, osc in zip(energy_npi, oscs_npi):
            results_npi.append([energy, osc])
            results_dict_npi[isample] = results_npi

        for energy, osc in zip(energy_pipi, oscs_pipi):
            results_pipi.append([energy, osc])
            results_dict_pipi[isample] = results_pipi

    return results_dict, results_dict_npi, results_dict_pipi

def plot_stem(energies_eV, oscs, clr):

    markerline, stemlines, baseline = plt.stem(energies_eV, 
                                               oscs,
                                               color = clr,
                                               markerfmt= ' ')
                                               #label=label)
    plt.setp(baseline, color='k', linewidth=1, alpha=0.)
    plt.setp(stemlines, linewidth=1, color=clr,alpha=0.5)


def plot_spectrum(x, data, FWHM, nsamples, norm=True,norm_const=None):
    """
    x = Frequency grid
    data = [frameid, energy, osc]
    FWHM = experimental value
    """

    totspec = np.zeros_like(x)
    prefactor = 2. * np.pi**2./c_au
    sigma = FWHM / (2 * (2 * np.log(2))**(0.5))

    for idx, states_data in data.items():
        for state_data in states_data:
            totspec += gaussian_broadening(x, state_data, sigma)
    totspec *= prefactor
    normalized_totspec = totspec/nsamples

    norm_constant = 1.
    if norm == True:
        if norm_const == None:
            norm_constant = max(normalized_totspec)
            normalized_totspec /= max(normalized_totspec)
        else:
            normalized_totspec /= norm_const
            norm_constant = norm_const
    return max(totspec), norm_constant, normalized_totspec.tolist()

def gaussian_broadening(x, val, sigma):
    """Normalized Gaussian"""
    return 1./np.sqrt(2.*np.pi*sigma**2.)*np.array( val[1] * np.exp(-0.5*(((x - val[0])) / (sigma))**2))



if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    nsamples = [1891,1785]
    runpaths = [
             "/data/projects/Pratip_MA_vs_AcAc/MA-data/ABSORPTION/CALC-SSSR/",
             "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/ABSORPTION-updated/CALC-SSSR/"
              ]
    Natoms = [9, 15]


    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/SI/FIGURE-S2/'

    # Helpful conversions
    eV_per_au = 27.2114
    eV_per_nm = 1239.84193

    greens = ['#4C9A2A','#68BB59','#1E5631','#4C9A2A', '#ACDF87', '#76BA1B', '#A4DE02']

    # Color list to cycle through
    color_list = ['b','orange', 'g','r','k','purple','yellow']
    color_ctr = 0

    # Based on fitting to experiment 
    FWHM = 0.24

    # General variables
    eV_shift = 0.096 ####0.095 #0.101
    #nsamples = 1784
    nstates = 2

    os.chdir(path_script)
    for molecule, npt, run_path  in zip(molecules, nsamples, runpaths):

        data = []
        # Read energies and oscillator strengths
        os.chdir(run_path)
        print(run_path)
        data = read_bagel_xmscaspt2_results(npt,nstates)
        os.chdir(path_script)
        print(run_path)
        print(data)

        #MAIN SCRIPT #data[0] contains the dictionary of all states
        freq_grid = np.linspace(2.5, 8.5, 8000)
        wavel = h_evs*c_m_per_s / (freq_grid+eV_shift) * 1.e9
        totmax, norm_constant, broadened = plot_spectrum(freq_grid, data[0], FWHM, npt, norm=True)


        #save S1 energy and osc into list #data[1] contains dictionary for S1 states
        S1en = []
        S1osc = []
        for idx, states_data in data[1].items():
            for state_data in states_data:
                S1en.append(state_data[0]+eV_shift)
                S1osc.append(state_data[1])

        #print out spectra
        f=open(f'osc_S1_{molecule}.dat','w')
        for i,j in zip(S1en,S1osc):
            f.write('{}\t{}\n'.format(i,j))
        f.close()

        #save S2 energy and osc into list #data[2] contains dictionary for S2 states
        S2en = []
        S2osc = []
        for idx, states_data in data[2].items():
            for state_data in states_data:
                S2en.append(state_data[0]+eV_shift)
                S2osc.append(state_data[1])

        f=open(f'osc_S2_{molecule}.dat','w')
        for i,j in zip(S2en,S2osc):
            f.write('{}\t{}\n'.format(i,j))
        f.close()
            
        #print out spectra
        f_nm=open(f'spec_nm_{molecule}_298.15K_0.24.dat','w')
        for i,j in zip(wavel,broadened):
            f_nm.write('{}\t{}\n'.format(i,j))
        f_nm.close()

        f_ev=open(f'spec_ev_{molecule}_298.15K_0.24.dat','w')
        for i,j in zip(freq_grid+eV_shift,broadened):
            f_ev.write('{}\t{}\n'.format(i,j))
        f_ev.close()

        print('totmax', totmax)
        #os.chdir(run_path)

#        if molecule == "MA":
#            nsamples = 1890
#
#            dataa = []
#           # Read energies and oscillator strengths
#            os.chdir(run_path)
#            dataa = read_bagel_xmscaspt2_results(nsamples,nstates)
#            os.chdir(path_script)
#
#            #MAIN SCRIPT #data[0] contains the dictionary of all states
#            freq_grid = np.linspace(2.5, 8.5, 8000)
#            wavel = h_evs*c_m_per_s / (freq_grid+eV_shift) * 1.e9
#            totmax, norm_constant, broadened = plot_spectrum(freq_grid, dataa[0], FWHM, nsamples, norm=True)
#
#
#            #save S1 energy and osc into list #data[1] contains dictionary for S1 states
#            S1ena = []
#            S1osca = []
#            for idx, states_data in dataa[1].items():
#                for state_data in states_data:
#                    S1ena.append(state_data[0]+eV_shift)
#                    S1osca.append(state_data[1])
#
#            f=open(f'osc_S1_{molecule}.dat','w')
#            for i,j in zip(S1ena,S1osca):
#                f.write('{}\t{}\n'.format(i,j))
#            f.close()
#            
#            #save S2 energy and osc into list #data[2] contains dictionary for S2 states
#            S2ena = []
#            S2osca = []
#            for idx, states_data in dataa[2].items():
#                for state_data in states_data:
#                    S2ena.append(state_data[0]+eV_shift)
#                    S2osca.append(state_data[1])
#
#            f=open(f'osc_S2_{molecule}.dat','w')
#            for i,j in zip(S2ena,S2osca):
#                f.write('{}\t{}\n'.format(i,j))
#            f.close()
#            
#            #print out spectra
#            f_nm=open(f'spec_nm_{molecule}_298.15K_0.24.dat','w')
#            for i,j in zip(wavel,broadened):
#                f_nm.write('{}\t{}\n'.format(i,j))
#            f_nm.close()
#
#            f_ev=open(f'spec_ev_{molecule}_298.15K_0.24.dat','w')
#            for i,j in zip(freq_grid+eV_shift,broadened):
#                f_ev.write('{}\t{}\n'.format(i,j))
#
#            f_ev.close()
#
#            print('totmax', totmax)
#            os.chdir(run_path)
#
#        #print(S1en)
#        #print(S1osc)
#        #print(S1ac_en)
#        #print(S1ac_osc)

    os.chdir(path_script)
        # Plotting all 3 axes
    fig, (ax0,ax1,ax2) = plt.subplots(3, figsize=(5,9), sharex=True)

    #normalized experimental spectra
    exp_eV, exp_norm = np.loadtxt('data_UV_normalizedintensityVSeV_NAGAKURA_1977_Bull_Chme_Soc_Japan.txt', unpack=True)
    ma_eV, ma_norm = np.loadtxt('spec_ev_MA_298.15K_0.24.dat', unpack=True)
    acac_eV, acac_norm = np.loadtxt('spec_ev_AcAc_298.15K_0.24.dat', unpack=True)

    s1enma, s1oscma = np.loadtxt('osc_S1_MA.dat', unpack=True)
    s2enma, s2oscma = np.loadtxt('osc_S2_MA.dat', unpack=True)
      
    s1enacac, s1oscacac = np.loadtxt('osc_S1_AcAc.dat', unpack=True)
    s2enacac, s2oscacac = np.loadtxt('osc_S2_AcAc.dat', unpack=True)

        #area under experimental cutoff at around 5.82 
        #exper = np.trapz(exp_norm[0:93], exp_eV[0:93])

    # Plot decor
    plt.subplots_adjust(right=0.85, top=0.95, bottom=0.15)

    y = np.linspace(0.0, 1.2, 8000)
    # Plot broadened (and shifted) spectrum, filling distinct excitation regions
    ax0.fill_betweenx(y, x1=4.61, x2=4.71, facecolor='gold', alpha=1.0, zorder=1)
    ax0.plot(exp_eV, exp_norm, color='k', linewidth=3.0, label='Expt.', zorder=4)   #(Expt.)
    ax0.plot(acac_eV, acac_norm, color = 'crimson', lw=3.0, zorder=5, linestyle='--')  
    ax0.plot(ma_eV, ma_norm, color='darkgray', linewidth=3.0, zorder=6, linestyle='-.')
    ax0.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
    ax0.set_xticks([3.0,4.0,5.0,6.0,7.0])
    ax0.legend(loc='upper right', frameon=False, framealpha=0.4, fontsize=12)
    ax0.set_ylabel('Norm. Intensity', fontsize=16)
    #ax0.set_xlabel('Energy (eV)', fontsize=16)
    ax0.set_ylim(0.0,1.10)
    ax0.set_xlim(3.0,7.0)
    ax0.tick_params(axis='both', which='major', labelsize=12)
    ax0.tick_params(axis='both', which='minor', labelsize=12)
    ax0.minorticks_on()

    #every 2nd point
        #print(S2ena)
    #S2ena = S2ena[::2]
    #S2osca = S2osca[::2]

    ax3 = ax1.twinx()
    y = np.linspace(0.0, 2.5, 20000)
    ax1.plot(ma_eV, ma_norm, color = 'darkgray', label = 'MA Calc.\n (+0.096 eV)', lw=3.0, zorder=5, linestyle='-.')
    ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
    ax1.set_xticks([3.0,4.0,5.0,6.0,7.0])
    ax1.legend(loc='upper right', frameon=False, framealpha=0.4, fontsize=12)
    ax1.set_ylabel('Norm. Intensity', fontsize=16)
    #ax1.set_xlabel('Energy (eV)', fontsize=16)
    ax1.set_ylim(0.0,1.10)
    ax1.set_xlim(3.0,7.0)
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.tick_params(axis='both', which='minor', labelsize=12)
    ax1.set_zorder(ax3.get_zorder()+1)
    ax1.set_frame_on(False)
    ax3.vlines(x=s2enma, ymin=0, ymax=s2oscma, colors='limegreen', alpha = 0.2, ls='-', lw=0.5, label=r'S$_2$', zorder=2)
    ax3.fill_betweenx(y, x1=4.61, x2=4.71, facecolor='gold', alpha=1.0, zorder=1)
    ax3.vlines(x=s1enma, ymin=0, ymax=s1oscma, colors='blue', alpha = 0.2, ls='-', lw=0.5, label=r'S$_1$', zorder=3)  ##eb6d5cff
    ax3.minorticks_on()
    ax3.set_ylabel('Osc. Strength', fontsize=16)
    ax3.set_ylim(0.0,1.40)
    ax3.set_xlim(2.0,9.0)
    ax3.set_xlim(3.0,7.0)
    ax3.tick_params(axis='both', which='major', labelsize=12)
    ax3.tick_params(axis='both', which='minor', labelsize=12)
    ax3.legend(loc='upper left', fontsize=14, frameon=False)


        #every 2nd point
        #S2en = S2en[::2]
        #S2osc = S2osc[::2]

    ax4 = ax2.twinx()
    y = np.linspace(0.0, 2.5, 20000)
    ax2.plot(acac_eV, acac_norm, color = 'crimson', label = 'AcAc Calc.\n (+0.096 eV)', lw=3.0, zorder=5, linestyle='--')
    ax2.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
    ax2.set_xticks([3.0,4.0,5.0,6.0,7.0])
    ax2.legend(loc='upper right', frameon=False, framealpha=0.4, fontsize=12)
    ax2.set_ylabel('Norm. Intensity', fontsize=16)
    ax2.set_xlabel('Energy (eV)', fontsize=16)
    ax2.set_ylim(0.0,1.10)
    ax2.set_xlim(3.0,7.0)
    ax2.minorticks_on()
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='minor', labelsize=12)
    ax2.set_zorder(ax4.get_zorder()+1)
    ax2.set_frame_on(False)

    ax4.vlines(x=s2enacac, ymin=0, ymax=s2oscacac, colors='limegreen', alpha = 0.2, ls='-', lw=0.5, label=r'S$_2$', zorder=2)
    ax4.fill_betweenx(y, x1=4.61, x2=4.71, facecolor='gold', alpha=1.0, zorder=1)
    ax4.vlines(x=s1enacac, ymin=0, ymax=s1oscacac, colors='blue', alpha = 0.2, ls='-', lw=0.5, label=r'S$_1$', zorder=3)  ##eb6d5cff
    ax4.minorticks_on()
    ax4.set_ylabel('Osc. Strength', fontsize=16)
    ax4.set_ylim(0.0,1.40)
    ax4.set_xlim(2.0,9.0)
    ax4.set_xlim(3.0,7.0)
    ax4.legend(loc='upper left', fontsize=14, frameon=False)
    ax4.tick_params(axis='both', which='major', labelsize=12)
    ax4.tick_params(axis='both', which='minor', labelsize=12)

    fig.subplots_adjust(hspace=0.02)
    fig.tight_layout()
    fig.savefig('abs_spec_xmscaspt2_sssr_prelim_ev_298.15K_1784IC_0.24_pub2.png', dpi=400)
    fig.savefig('abs_spec_xmscaspt2_sssr_prelim_ev_298.15K_1784IC_0.24_pub2.pdf', dpi=400)

#clean axis for new plot
#ax.cla()
#ax2.cla()


