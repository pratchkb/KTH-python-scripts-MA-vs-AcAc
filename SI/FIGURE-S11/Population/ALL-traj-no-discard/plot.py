from typing import List, Optional

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

def bootstrapped_dict(
    bootstrapped_pops, 
    times, 
    states, 
    state_clrs={0:'black', 1:'limegreen', 2:'orange'},
    plot_accum=False,
    ):
    accum = None
    for I in states:
        vals = np.array(bootstrapped_pops[I])
        avg = np.average(vals, axis=0)
        plt.plot(times, avg, color=state_clrs[I], lw=3.0, label=f"S$_{I}$")
        if accum is None:
            accum = np.zeros_like(avg)
        accum += avg
        std = np.std(vals, axis=0)
        plt.fill_between(times, avg-std, y2=avg+std, alpha=0.4, color=state_clrs[I], linewidth=0)
    if plot_accum:
        plt.plot(times, accum, 'k')
    plt.xlim([times[0], times[-1]])
    plt.ylim([0, 1.01])
    return plt.gca()
