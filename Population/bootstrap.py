import copy
import numpy as np


def population(
    pop_dict: list[dict],
    nsamples: int,
    column_wise=False,
    ) -> dict[int: np.ndarray]:
    """Bootstraps the given simulation for for 'nsamples'

    Args:
        pop_dict (list of dicts): The merged and interpolated populations
        nsamples (int): number of resamples for the bootstrap
        column_wise (bool): Bootstrap individual values between ICs
                            instead of whole population decays

    Returns:
        Dict(int: np.ndarray): Dict containing 'nsamples' bootstrapped popoulations per state
    """

    # bootstrap the popoulation
    if column_wise:
        return bootstrap_ic_population_column(pop_dict, nsamples)
            
    
    return bootstrap_ic_population_row(pop_dict, nsamples)


def bootstrap_ic_population_column(pop_dict, nsamples):
    ret_dict = {}
    
    # Loop over states
    for I, values in pop_dict.items():
        mat = np.array(values)
        
        vals = np.zeros((nsamples, mat.shape[1]))
        for i in range(mat.shape[1]):
            
            # Each column is sampled 'nsamples' times
            # So e.g. 50 ICs and 1000 samples
            sample = np.random.choice(mat[:,i], (mat.shape[0], nsamples), replace=True)

            # generate a 50x1000 per column
            # each column then takes the average of the 50 vals
            # so giving 1000 rows per column
            vals[:,i] = np.average(sample, axis=0)
            
        ret_dict[I] = vals
        # Finally ending up with an array of 
        # Column x Samples array per state
    return ret_dict



def bootstrap_ic_population_row(pop_dict, nsamples, rest_pop_array=None):
    ret_dict = {}
    rest_pop = {}
    # Loop over states
    for I, values in pop_dict.items():
        mat = np.array(values)
        
        rest_pop_samples = []
        
        vals = np.zeros((nsamples, mat.shape[1]))
        for i in range(nsamples):
            
            # For nsamples a random batch of IC samples are averaged
            # So e.g. 50 ICs and 1000 samples
            rand_ics = np.random.choice(range(mat.shape[0]), (mat.shape[0]))
            sample = mat[rand_ics]
            
            # If present, handle rest population
            # Create a 1:1 mapping of the bootstrapped population to the leaked/stuck population
            if rest_pop_array is not None:
                rest_pop_samples.append(np.average([rest_pop_array[ic][I-1] for ic in rand_ics]))

            # generate a 50x1000 per column
            # each column then takes the average of the 50 vals
            # so giving 1000 rows per column
            vals[i,:] = np.average(sample, axis=0)
            
        ret_dict[I] = vals
        
        # If present, handle rest population
        if rest_pop_array is not None:
            rest_pop[I] = rest_pop_samples
        
    # Finally ending up with an array of 
    # Column x Samples array per state
    if rest_pop_array is not None:
        return ret_dict, rest_pop
    
    return ret_dict


