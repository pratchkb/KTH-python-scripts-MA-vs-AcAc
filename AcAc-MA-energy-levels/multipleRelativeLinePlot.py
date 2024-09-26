import sys
import matplotlib.pyplot as plt
import matplotlib 
from matplotlib import rc, rcParams
import math
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 2.0
#matplotlib.rc('text', usetex = True)

#rcParams['text.usetex'] = True
#rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
#
## Attempt to override usetex font BS!
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Helvetica']

au2ev = 27.211386245988
GREEN = [0, 0.69, 0.31] 

class Plotdata(object):
    f = None
    def __init__(self, filename):
        self.title          = ""
        self.states         = []
        self.structures     = []
        self.structure_data = {}

        self.parsefile(filename)

    def parsefile(self, filename):
        try:
            f = open(filename)

            i = 0
            for line in f:
                if i == 0:
                    self.title = line.strip()
                elif i == 1:
                    self.states = line.split()
                    self.data = [[] for _ in range(len(self.states))]
                    for state in self.states:
                        self.structure_data[state] = {}
                else:
                    if line[0] == '#':
                        continue
                    fields = line.split()
                    values = fields[-len(self.states)::]
                    for i, val in enumerate(values):
                        self.data[i].append(float(val))
                    structure = " ".join(fields[:-len(self.states)])
                    self.structures.append(structure)
                    for state,val in zip(self.states, values):
                        self.structure_data[state][structure] = float(val)
                i += 1
            print(self.title)
            print(self.states)
            print(self.structures)
            print(self.structure_data)

        except Exception as error:
            print('Error while reading file {} : {}'.format(filename, error))

        finally:
            if not f.closed:
                f.close()


def connectedLinePlot(yValues, clr, lStyle,base,lbl="", skip=0, alpha=1.):
    global au2ev
    yValues = [(x - base)*au2ev for x in yValues]
    offset = 0+(skip*2)
    l = len(yValues)
    handle = None
    print('lbl', lbl)
    for i in range(l):
        # Plot horisontal values
        handle, = plt.plot([offset+2*i, offset+2*i+1], 
			  [yValues[i], yValues[i]], 
			  color=clr, 
			  linestyle=lStyle,
			  linewidth=2.8,
			  label = lbl, 
			  alpha=alpha,zorder=99999)

        # Connects with lines in between
        if (i+1) < l:
            plt.plot([offset+2*i+1, offset+2*(i+1)], [yValues[i], yValues[i+1]], color=clr, linestyle='--',linewidth=0.5, alpha=alpha)

    return handle 

if __name__ == '__main__' and len(sys.argv) > 0:
    structure_order = ['MECI-P','S1-P','S1planar','S0','S1planar','S1-I','MECI-I','S0-trans','S1planar-trans']
    labels = [r'MECI-P', r'S$_1$-P', r'S$_1$-planar', r'S$_0$-min',r'S$_1$-planar',r'S$_1$-I', r'MECI-I',r'S$_0$-min',r'S$_1$-planar'] 
    fig = plt.figure(figsize=(12,6))
    names = [name.split('.')[0] for name in sys.argv[1:]]
    legend_names = []
    savename = "_".join(names)
    files = sys.argv[1:]
    # Create a list of alpha channels values based on the number of files
    #alpha_channels = [float(x)/len(files) for x in range(1,len(files)+1)]
    alpha_channels = [1.0 for x in range(1,len(files)+1)]
    plthandle = []
    all_colors = [['#000000ff','#9e9e9eff'],['#9f1722ff','#eb6d5cff'], ['#1470b0ff','#53bde7ff']]
    for colors, a_ch, arg in zip(all_colors, alpha_channels, files):
        plotfile = Plotdata(arg) 
        ax = plt.gca()
        plt.subplots_adjust(bottom=0.20,top=0.95)
        states_to_plot = ['S0','S1']
        for clr,state in zip(colors,states_to_plot):
#            legend_names.append("{} {}".format(arg.split('.')[0], state)) 
            data = []
            S0 = plotfile.structure_data['S0']['S0']
            for structure in structure_order:
                data.append(plotfile.structure_data[state][structure] - S0)  
                print('{0:15s} {1:4s} {2:6.3f}'.format(structure, state, (plotfile.structure_data[state][structure] - S0)*au2ev))
            plthandle.append(connectedLinePlot(data, clr,'-',0, "{} {}".format(arg.split('.')[0], state), skip=0, alpha=a_ch))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.axhline(y=0, color='k',linestyle='-',zorder=1)
        ax.axhline(y=(plotfile.structure_data['S1']['S0']-S0)*au2ev, color='grey',linestyle='--',zorder=1,lw=0.6)

        ## x/y-limits
        plt.xlim([-0.3, 2*len(labels)-0.5])
        ## Y-label
        plt.ylabel('Relative energy (eV)',fontsize = 16,labelpad=10)
        ## Moves the ticks to center under each horisontal line
        ax.set_xticks([0.5+2*i for i in range(len(labels))])
        ax.tick_params(axis=u'x', which=u'both',length=0)
        #ax.tick_params(axis=u'y', which=u'both',left=False,length=0)
        ax.tick_params(axis=u'y',right=False)
        ax.tick_params(axis=u'y',left=True)
        ax.set_xticklabels(labels)
        #ax.spines['left'].set_visible(False)
        #
        ## And changes tick labels from numbers to labels
        #ax.set_yticklabels(["{:.1f}".format(x-1) for x in range(int(ymax+2))])
        ax.tick_params(labelsize=16,length=8,width=2.0)

    plt.legend(plthandle,legend_names,frameon=False,loc=2, ncol=len(legend_names), bbox_to_anchor=(0.08, 1.06))
    plt.savefig('relative_energy_plot.pdf', dpi=300,transparent=True)
