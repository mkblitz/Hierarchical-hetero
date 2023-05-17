# from __future__ import division
import sys
import cell_parent as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import utils
from neuron import h, gui

class many_sticks_model(cp.cell_parent):
    def __init__(self, num_of_sections, soma_L = 23, soma_diam = 718, soma_Ra = 100, dend_L = 50, dend_diam = 0.75, 
                 dend_Ra = 150, soma_cm = 1, dend_cm = 1*2, pas = True, soma_g_pas=0.0000338 , dend_g_pas = 0.0000467, soma_e_pas = -77, #0.00005*2
                 dend_e_pas = -77,):
        super().__init__()
#         print(soma_g_pas, dend_g_pas )
        self.model_name = 'many_sticks'
        self.num_of_sections = num_of_sections
        'creating the model'
        'soma'
        self.soma = h.Section(name="soma")
        self.soma.L = soma_L
        self.soma.diam = soma_diam
        self.soma.Ra = soma_Ra
        self.soma.cm = soma_cm
        
        'dendrites'
        self.dendrites = []
        for sec_num in range(num_of_sections):
            dend = h.Section(name="dend_sec_"+str(sec_num))
            dend.L = dend_L#should maybe be random?
            dend.diam = dend_diam
            dend.Ra = dend_Ra
            dend.cm = dend_cm
            dend.nseg = dend_L//5
            if pas:#adding the passive properties
                dend.insert('pas')
                dend.g_pas = dend_g_pas
                dend.e_pas = dend_e_pas
            self.dendrites.append(dend)
        
        'spines'
        self.heads = []
        self.necks = []
        
        'adding the passive properties'
        if pas: 
            self.soma.insert('pas')
            self.soma.g_pas = soma_g_pas
            self.soma.e_pas = soma_e_pas
        
        'connecting the soma to the dendrite'
        for sec_num in range(len(self.dendrites)):
            if sec_num == 0:
                self.dendrites[0].connect(self.soma, 1, 0)
            else:
                connect_to = utils.find_dend_sec_path(sec_num)[-2]
                self.dendrites[sec_num].connect(self.dendrites[connect_to], 1, 0)
        
        self.dend_TRs, self.dend_ratios = self.get_dend_TRs() #Matrix of pairwise transfer resistances
        self.soma_TRs, self.soma_ratios = self.get_soma_TRs() #Matrix of soma
        
#         print(self.dendrites)
# cell = many_sticks_model(31)

        
        
        
        
        
        