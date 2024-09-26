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
    yValues = [(x - base) for x in yValues]
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
    structure_order = ['FC','Asym HTI','Sym HTI', 'S2/S1 tw', 'T1 min', 'CIpy', 'S1/S0 tw']
    labels = [r'S$_0$-min', r'S$_2$/S$_1$-asym-HTI', r'S$_2$/S$_1$-sym-HTI', r'S$_2$/S$_1$-tw', r'T$_1$-min', r'CIpy', r'S$_1$/S$_0$-tw'] 

    fig = plt.figure(figsize=(6,8))
    names = [name.split('.')[0] for name in sys.argv[1:]]
    legend_names = [r'S$_0$', r'S$_1$', r'S$_2$', r'T$_1$', r'T$_2$']
    savename = "_".join(names)
    files = sys.argv[1:]
    # Create a list of alpha channels values based on the number of files
    #alpha_channels = [float(x)/len(files) for x in range(1,len(files)+1)]
    alpha_channels = [1.0 for x in range(1,len(files)+1)]
    plthandle = []
    #all_colors = [['#800000','#800000','#800000'], ['#FF0000','#FF0000','#FF0000'], ['#FF6600','#FF6600','#FF6600'], ['#FFCC00','#FFCC00','#FFCC00']] #[['#000000ff','#9e9e9eff','orange'],['#9f1722ff','#eb6d5cff','blue'], ['#1470b0ff','#53bde7ff','green']]
    #all_colors = [['black','black','black'],['#eb6d5cff','#eb6d5cff','#eb6d5cff'], ['orange','orange','orange'], ['#2A7FFF','#2A7FFF','#2A7FFF']]
    all_colors = [['darkred','orangered','orange','limegreen','darkgreen']]#,['#BA232D','#BA232D','#BA232D'], ['orange','orange','orange'], ['#0076D2','#0076D2','#0076D2']]
    for colors, a_ch, arg in zip(all_colors, alpha_channels, files):
        plotfile = Plotdata(arg) 
        ax = plt.gca()
        plt.subplots_adjust(bottom=0.20,top=0.95)
        states_to_plot = ['S0', 'S1', 'S2', 'T1', 'T2']
        for clr,state in zip(colors,states_to_plot):
#            legend_names.append("{} {}".format(arg.split('.')[0], state)) 
            data = []
            S0 = plotfile.structure_data['S0']['FC']
            for structure in structure_order:
                data.append(plotfile.structure_data[state][structure])  
                print('{0:15s} {1:4s} {2:6.3f}'.format(structure, state, (plotfile.structure_data[state][structure])))
            plthandle.append(connectedLinePlot(data, clr,'-',0, "{} {}".format(arg.split('.')[0], state), skip=0, alpha=a_ch))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.axhline(y=0, color='k',linestyle='-',zorder=1)
        
        ax.axhline(y=4.717,color='gray', linestyle='--', zorder=1, lw=0.6)
        ax.axhline(y=3.643, color='lightgray',linestyle='--',zorder=2,lw=0.6)
        #ax.axhline(y=(plotfile.structure_data['S1']['S1min']-S0)*au2ev, color='gray',linestyle='--',zorder=1,lw=0.6)
        #ax.axhline(y=(plotfile.structure_data['S2']['S0']-S0)*au2ev, color='lightgray',linestyle='--',zorder=2,lw=0.6)

        ## x/y-limits
        plt.xlim([-0.3, 2*len(labels)-0.5])
        #plt.xticks([0,1,2,3,4,5,6,7])
        ## Y-label
        plt.ylim([1,6.4])
        plt.yticks([0,1,2,3,4,5,6])
        plt.ylabel('Relative energy (eV)',fontsize = 20,labelpad=10)
        ## Moves the ticks to center under each horisontal line
        ax.set_xticks([0.5+2*i for i in range(len(labels))])
        ax.tick_params(axis=u'x', which=u'both',length=0)
        #ax.tick_params(axis=u'y', which=u'both',left=False,length=0)
        #ax.tick_params(axis=u'y',right=False)
        #ax.tick_params(axis=u'y',left=True)
        ax.set_xticklabels(labels)
        ax.tick_params(axis='x', labelrotation=45)
        #ax.spines['left'].set_visible(False)
        #
        ## And changes tick labels from numbers to labels
        #ax.set_yticklabels(["{:.1f}".format(x-1) for x in range(int(ymax+2))])
        ax.tick_params(labelsize=20,length=8,width=2.0)

    #plt.legend(plthandle,legend_names,frameon=False,fontsize=20, loc ='lower right')# bbox_to_anchor=(0.08, 1.06))    #ncol=len(legend_names)
    plt.tight_layout()
    plt.savefig('MA_relative_energy_plot_test_all_version2.pdf', dpi=1000,transparent=True)
