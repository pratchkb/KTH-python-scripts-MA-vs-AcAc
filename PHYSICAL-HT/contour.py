import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator
from scipy import ndimage

list_file = ['MA-data.dat', 'AcAc-data.dat']

#
for file in list_file: 
    data = np.genfromtxt(file, dtype='float64')

    xedges = np.linspace(0,200,100)
    yedges = np.linspace(-3.0,3.0,100)
    X, Y = np.meshgrid(xedges, yedges)
    sigma_x = 1.0
    sigma_y = 0.06
    gau2d = np.exp(-0.5 * (X-100)**2/sigma_x**2 - 0.5 * (Y)**2/sigma_y**2 )

    H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=100, density=True, range=[[0,200],[-3.0,3.0]])

#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.

    H = H.T
    N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)   

# 
    fig, ax = plt.subplots(figsize=(7,5))
    cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
    plt.xlim(0,200)
    plt.ylim(-3.0,3.0)
    plt.xlabel("Time (fs)", fontsize=18)
    plt.ylabel(r"HT ($\AA$)", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(str(file) + '-contourf-gau.pdf', dpi=300)
