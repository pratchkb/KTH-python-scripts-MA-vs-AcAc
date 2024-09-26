import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#import matplotlib.colors
#import matplotlib.colors.PowerNorm
#import matplotlib.axes.Axes.hist2d
#import matplotlib.pyplot.hist2d
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator
from scipy import ndimage

list_file = ['final-torsion-ortho-plane-QT-all-final-pulse-total-0.24.dat','final-torsion-ortho-plane-QT-all-final-wigner_298.15K1.2.dat','final-torsion-ortho-plane-QT-all-final-complete-final-QT-geom.dat']

for file in list_file: 
    data = np.genfromtxt(file, dtype='float64')

    xedges = np.linspace(-60,120,180)
    yedges = np.linspace(-60,120,180)
    X, Y = np.meshgrid(xedges, yedges)
    sigma_x = 2.
    sigma_y = 2.
    gau2d = np.exp(-0.5 * (X-30)**2/sigma_x**2 - 0.5 * (Y-30)**2/sigma_y**2 )
    H, xedges, yedges = np.histogram2d(data[:,7], data[:,9], bins=180, density=True, range=[[-60,120],[-60,120]])


    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
    # therefore transpose H for visualization purposes.
    H = H.T

    N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)    
 
    fig, ax = plt.subplots(figsize=(5,5))
    plt.rcParams['scatter.marker'] = 'x'
    plt.rcParams['lines.markersize'] = '10'
    
    #cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),21), zorder=10)
    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
    #fig.colorbar(cbar, ax=ax)
    ax.set_aspect('equal')
    ax.axvline(x=30, color='gray', linestyle='--', linewidth=0.5, zorder=10)
    ax.axvline(x=60, color='gray', linestyle='--', linewidth=0.5, zorder=6)
    #ax.axvline(x=-120, color='gray', linestyle='--', zorder=2, linewidth=0.5)
    #ax.axvline(x=120, color='gray', linestyle='--', zorder=3, linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, zorder=7)
    ax.axhline(y=30, color='gray', linestyle='--', zorder=11, linewidth=0.5)
    ax.axhline(y=60, color='gray', linestyle='--', linewidth=0.5, zorder=8)
    #ax.axhline(y=-120, color='gray', linestyle='--', zorder=7, linewidth=0.5)
    #ax.axhline(y=120, color='gray', linestyle='--', zorder=8, linewidth=0.5)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, zorder=9)
    fig.gca()
    cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),21))
    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
    #cbar2 = ax.contour(X, Y, N_gau_conv, colors='g', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),21)
#    fig.colorbar(cbar, ax=ax)
    plt.scatter(1.1,58.4, marker = 'x', color = 'k', linewidth=2, zorder=2)
    plt.scatter(58.2,58.5, marker = 'x',  color = 'k',  linewidth=2, zorder=3)
    plt.scatter(1.1,0.6, marker = 'x', color = 'k', linewidth=2, zorder=4)
    plt.scatter(57.2,0.6, marker = 'x', color ='k', linewidth=2, zorder=5)
    plt.xlim([-30.0,90.0])
    plt.ylim([-30.0,90.0])
    plt.xticks([-30,0,30,60,90])
    plt.yticks([-30,0,30,60,90])
    plt.xlabel(r"min$_{x}\angle$A ($\degree$)", fontsize=18)
    plt.ylabel(r"min$_{y}\angle$B ($\degree$)", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.savefig(f'{file}-contourf-gau.pdf', dpi=300)






