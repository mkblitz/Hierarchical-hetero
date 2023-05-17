import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)
import model_run as mr
import matplotlib.pyplot as plt

sim_time = 100 
vs = []

num_of_sections = 15#how many sections the dendrite will have. for a symmetrical tree (1 turns to 2), symmetrical options are 1,3,7,15,31,63,127...
branch_L = 50#length of each section of the dendrite
num_of_stims = 1#number of stimulations at every active synapse
freq = 200#frequency of stimulations
syn_step = 1000/freq
start_1 = 10#time of first stimulation after start of simulation
syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)

num_of_active_prox_syns = 18#number of active proximal (branch) synapses
num_of_active_dist_syns = 0#number of active distal (branch) synapses

stims_prox = {1+i:syn_times for i in range(num_of_active_prox_syns)}
stims_dist = {1+i:syn_times for i in range(num_of_active_dist_syns)}

'set what branchs get active synapses. branches are counted from soma = #0, and then left to right in ascending order.' 
'in 2D on branch #4 there are 18 active synapses'
synapses = {1:{0.5:[]}, 4:{0.5:stims_prox}, 9:{0.5:stims_dist}}

'for voltage trace, set True'
voltage_trace = True
branch_num= 4 #1,3,4,7 

print('synapses = ',synapses)#will print: {branch_number: {location_on_branch: {active_synapses_number: array_of_stimulation_times_at_that_synapse}}}
    
cell, voltage_matrices = mr.multiple_dendrite_heatmaps(num_of_sections, synapses, num_of_active_prox_syns, num_of_active_dist_syns, 
                                                       sim_time = sim_time, dend_L = branch_L, voltage_trace = voltage_trace, branch_num = branch_num)


