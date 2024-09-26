import numpy as np
import sys
import os
import logging
import pickle
import time
import glob
from multiprocessing import Pool
import matplotlib.pyplot as plt

# special external packages goes here

# Our own packages
import bootstrap
import plot


#CACHE_FILENAME = "cacheacac.pickle"

def read_ndat(filename):
    try:
        time.sleep(1)
        print(f"opening {filename}")
        with open(filename, 'r') as f:
            logging.info(f"Opens {filename}")
            return np.loadtxt(f)
    except ValueError as e:
        logging.error(f"Failed to read file {filename} due to error {e}")
        return None

def read_population_data(data_path, pops):
    files_to_read = []
    for filename in glob.glob(f"{data_path}/?????/N.dat"):
        if filename in pops:
            continue
        files_to_read.append(filename)
    
    with Pool(14) as p:
        populations = p.map(read_ndat, files_to_read)
    for filename, p in zip(files_to_read, populations):
        if p is not None:
            pops[filename] = p

    return pops

def load_cache(filename):
    try: 
        with open(filename, 'rb') as f:
            return pickle.load(f)
    except:
        return {}

def save_cache(filename: str, data):
    """ Will save the data into a pickle file

    Args:
        filename : string : the filename wo which to write
        data : 
    """
    try:
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
    except Exception as e:
        logging.error(f"damn - could not save, {e}")

def load_or_process_population_data(path):
    raw_data = load_cache(CACHE_FILENAME)
    # Read all the raw data
    raw_data = read_population_data(data_path, raw_data)
    save_cache(CACHE_FILENAME, raw_data)
    return raw_data


def find_last_time(data):

    #print(data.values())
    return max([max(ndat[:,0]) for ndat in data.values()]) 


def interpolate_data(raw_data, times, states):
    interpolated_data = {
            s : [] for s in states
        }
    for ndat in raw_data.values():
        t = ndat[:,0]
        for s in states:
            s_data = ndat[:,1+s]
            s_interp = np.interp(times, t, s_data)
            interpolated_data[s].append(s_interp)

    return interpolated_data


if __name__ == "__main__":

    molecules = ['MA', 'AcAc']
    #Ntraj = [257,272]
    #runpaths = [
    #         "/data/projects/Pratip_MA_vs_AcAc/MA-data/production/PROD/ANALYSIS-MAY8/",
    #         "/data/projects/Pratip_MA_vs_AcAc/AcAc-data/PROD_Final_Feb27/"
    #          ]
    #Natoms = [9, 15]
    axis = [0,1]

    path_script = '/data/projects/Pratip_MA_vs_AcAc/PLOTS/Population/'

    fig,ax = plt.subplots(figsize=(6,5))

    for molecule in molecules:
        CACHE_FILENAME = f"cache_{molecule}.pickle"    

        logging.basicConfig(
            filename='tut.log',
            format="[%(levelname)s] %(asctime)s : %(message)s",
            encoding='utf-8',
            level=logging.INFO,
            )
        logging.info(f'Started running test.py: args: {sys.argv}, running git commit...')
        data_path = path_script + f'{molecule}_pop_data/'
        dt = 0.5

        print(data_path)
        raw_data1 = {}
        raw_data1 = load_or_process_population_data(data_path)
        last_time1 = find_last_time(raw_data1)
        times1 = np.arange(0, last_time1, dt)
        interpolated_data1 = interpolate_data(raw_data1, times1, [0, 1, 2])

        # Calculate the boostrapped population
        bootstrapped_pop1 = bootstrap.population(interpolated_data1, 2000)
        ax.clear()
        ax = plot.bootstrapped_dict(bootstrapped_pop1, times1, [0, 1, 2])

        # Plot it
        ax.set_xlabel('Time (fs)', fontsize=16)
        ax.set_ylabel('Adiabatic Population', fontsize=16)
        ax.set_ylabel('Adiabatic Population', fontsize=16)
        ax.set_xlim(0,200)
        ax.set_ylim(0,1)
        ax.set_ylim(0,1)
        ax.tick_params(labelsize=16)
        ax.minorticks_on()
        ax.minorticks_on()
        plt.legend(loc = 'center right', frameon=False, fontsize=16)
        plt.tight_layout()
        plt.savefig(path_script + f'population_{molecule}.pdf', dpi=500)
