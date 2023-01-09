from __future__ import division
import sys
# path_name = 'C:/Users/mkbli/Dropbox/Deep Dendrite/mods/'
# sys.path.append(path_name)
# path_name = 'C:/Users/mkbli/Dropbox/Deep Dendrite/'
# sys.path.append(path_name)
# sys.path.append("C://nrn")
# sys.path.append("C://nrn//bin")
# sys.path.append("C://nrn/lib/python/")
#from neuron import gui
from neuron import h
#h.nrn_load_dll('C:/Users/mkbli/Dropbox/Neuron/NEURON_tutorial-master/Python_scripts/nrnmech.dll')
import cell_parent as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import utils
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcyk')
from scipy.optimize import curve_fit
#h.nrn_load_dll('C:/Users/mkbli/Dropbox/Deep Dendrite/nrnmech.dll')

class stick_and_ball_model(cp.cell_parent):
    '''
        stick and ball model parent class
        
        parameters:
        nseg: number of dendritic segments
        rec_locations: List of locations on dendrite to record
        
        next few lines are same for dend
        soma_L: soma length [um]
        soma_diam: soma diameter [um]
        soma_Ra: soma resistance [ohm*cm]
        soma_g_pas: specific membrane resistance [ohm*cm^2]
        soma_cm: capacitance 

        pas: passive properties of the neuron, default = True 
        rec_soma_v: record soma voltage
        rec_stim = record stimulation current
    '''
    def __init__(self, soma_L = 23, soma_diam = 718, soma_Ra = 100, dend_L = 200, dend_diam = 0.75, 
                 dend_Ra = 150, soma_cm = 1, dend_cm = 2, pas = True, soma_g_pas = 0.0000338 , dend_g_pas = 0.0000467, soma_e_pas = -77, 
                 dend_e_pas = -77,):#Hay model - soma_g_pas=0.0000338, dend_g_pas = 0.0000467
        super().__init__()
        self.model_name = 'ball_and_stick'
        'creating the model'
        'soma'
        self.soma = h.Section(name="soma")
        self.soma.L = soma_L
        self.soma.diam = soma_diam
        self.soma.Ra = soma_Ra
        self.soma.Ra = 0.1
#         print(self.soma.Ra)
        self.soma.cm = soma_cm
        'dendrite'
        self.dendrites = []
        self.dend = h.Section(name="dend")
        self.dend.L = dend_L
        self.dend.diam = dend_diam
        self.dend.Ra = dend_Ra
        self.dend.cm = dend_cm
        self.dend.nseg = dend_L//4
        self.dendrites.append(self.dend)
        
        'spines'
        self.heads = []
        self.necks = []
        
        'adding the passive properties'
        if pas: 
            self.soma.insert('pas')
            self.soma.g_pas = soma_g_pas
            self.soma.e_pas = soma_e_pas
            self.dend.insert('pas')
            self.dend.g_pas = dend_g_pas
            self.dend.e_pas = dend_e_pas
        
        'connecting the soma to the dendrite'
        self.dend.connect(self.soma, 1, 0)
        
        'work this out'
        self.dend_TRs, self.dend_ratios = self.get_dend_TRs() #Matrix of pairwise transfer resistances
        self.soma_TRs, self.soma_ratios = self.get_soma_TRs() #Matrix of soma
#         print(self.dend.Ra, self.dend.diam)

# def absexp(x, a, c, b):
#     return (a+b*x)*np.exp(np.abs(c*x))
# print('aa')
# # if __name__ == "__main__":
# print('a')
# plt.ion()
# cell = stick_and_ball_model(dend_L = 500)
# cell.plot_TRs()
# plt.show()
# ('end')


        


       