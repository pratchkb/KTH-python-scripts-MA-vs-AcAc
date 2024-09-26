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

list_file = ['final-OH-newcoord-torsion-wigner_298.15K.dat', 'final-OH-newcoord-torsion-complete-QT-traj-combined.dat', 'final-OH-newcoord-torsion-complete-final-QT-geom.dat', 'final-OH-newcoord-torsion-pulse-total-0.24.dat']

#'final-OH-newcoord-torsion-complete-QT-traj-combined.dat'



for file in list_file: 
    data = np.genfromtxt(file, dtype='float64')

#    print(data)

    mb = np.mean(data[:,1])
    mt = np.mean(data[:,2])

    xedges = np.linspace(-2.0,2.0,80)
    yedges = np.linspace(0,80,80)
    X, Y = np.meshgrid(xedges, yedges)
    sigma_x = 0.15
    sigma_y = 1.
    gau2d = np.exp(-0.5 * (X-0.0)**2/sigma_x**2 - 0.5 * (Y-40)**2/sigma_y**2 )

    H, xedges, yedges = np.histogram2d(data[:,1], data[:,2], bins=80, density=True, range=[[-2.0,2.0],[0,80]])
#
#
#    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
#    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
#    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
#    # therefore transpose H for visualization purposes.
    H = H.T
    N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)   
# 
    fig, ax = plt.subplots(figsize=(5,5))
    cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#    cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#    #pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
#    #fig.colorbar(cbar, ax=np.linspace(0,1,21))
#    fig.colorbar(cbar, ax=ax)
#    #ax.set_aspect('equal')
    plt.xlim(-1.8,1.8)
    plt.xticks([-1.8,-1.2,-0.6,0.0,0.6,1.2,1.8])
    plt.ylim(0,80)
    plt.xlabel(r"HT ($\AA$)", fontsize=18)
    plt.ylabel(r"TorsionC ($^{\circ}$)", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.savefig(str(file) + '-gau-' + '-contourf-gau.pdf', dpi=300)


bond = []
torsions = []
datap = np.genfromtxt('wigner-CH-OH-newcoord-torsion-corr-HT.dat', dtype='float64')

for i in range(len(datap[:,1])):
    if datap[i,1] < 1.2:
        bond.append(datap[i,4])
        torsions.append(datap[i,5])

mb = np.mean(bond)
mt = np.mean(torsions)
 
xedges = np.linspace(-2.0,2.0,80)
yedges = np.linspace(0,80,80)
X, Y = np.meshgrid(xedges, yedges)
sigma_x = 0.2
sigma_y = 1.
gau2d = np.exp(-0.5 * (X-0.0)**2/sigma_x**2 - 0.5 * (Y-40.0)**2/sigma_y**2 )

#    sigma_x = 1.
#    sigma_y = 5.
#    gau2d = np.exp(-0.5 * (X)**2/sigma_x**2 - 0.5 * (Y)**2/sigma_y**2 )

H, xedges, yedges = np.histogram2d(bond, torsions, bins=80, density=True, range=[[-2.0,2.0],[0,80]])


    # Histogram2D in matplotlib does not follow Cartesian convention (see Notes),
    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa and y values on the ordinate axis. 
    # Rather, x is histogrammed along the first dimension of the array (vertical), and y along the second dimension of the array (horizontal).
    # therefore transpose H for visualization purposes.
H = H.T

N_gau_conv = ndimage.convolve(H, gau2d, mode='constant', cval=0.0)

fig, ax = plt.subplots(figsize=(5,5))
#cbar = ax.contourf(X, Y, gau2d, cmap='gist_heat_r')

print(N_gau_conv.min())
print(N_gau_conv.max())
cbar = ax.contourf(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r', levels=np.linspace(N_gau_conv.min(),N_gau_conv.max(),51), zorder=10)
#cbar = ax.contourf(X, Y, H, norm=colors.Normalize(vmin=H.min(), vmax=H.max()), cmap='gist_heat_r', levels=np.linspace(H.min(),H.max(),21), zorder=10)
#fig.colorbar(cbar, ax=ax)
#pcm = ax.pcolor(X, Y, H)
#pcm = ax.pcolor(X, Y, N_gau_conv, norm=colors.Normalize(vmin=N_gau_conv.min(), vmax=N_gau_conv.max()), cmap='gist_heat_r')
plt.xlim([-1.8,1.8])
plt.xticks([-1.8,-1.2,-0.6,0.0,0.6,1.2,1.8])
plt.xlabel(r"HT ($\AA$)", fontsize=18)
plt.ylabel(r"TorsionC ($^{\circ}$)", fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
    #plt.ylim([-60.0,120.0])
    #plt.xticks([-60,0,60,120])
    #plt.yticks([-60,0,60,120])
    #plt.xlabel(r"torsion-plane-min ($\degree$)")
    #plt.ylabel(r"torsion-plane-min ($\degree$)")
plt.tight_layout()
plt.savefig(f'wigner_less_1.2-OH-newcoord-contourf-gau.pdf', dpi=300)




