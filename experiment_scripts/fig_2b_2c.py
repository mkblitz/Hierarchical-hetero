import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)
import model_run as mr

'choose model, ball_and_stick (bas) or Hay'
model = 'bas'#'bas','Hay'

'setting the parameters'
num_of_stims = 1#number of stimulations at every active synapse
num_of_active_syns = 18#number of active synapses
freq = 200#frequency of stimulations

syn_step = 1000/freq
start_1 = 10#time of first stimulation after start of simulation
syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)

stims = {1+i:syn_times for i in range(num_of_active_syns)}

'adding the synapses - see structure on line 28'
active_loc = 0.5#location of active synapses, as a fraction of branch length
none = {1:[]}
synapses = {0:{}}
synapses[0][0.3] = none#adds an inactive spine at 60nm
synapses[0][active_loc] = stims#adds the active synapses at 100nm
synapses[0][0.99] = none#adds an inactive spine at 200nm

'if you want to scale Ih, set True. to control scaling go to model_run -> dendrite_heatmap, otherwise: seg.gIhbar_Ih = 1e-4*(-2+4.28*np.exp(distance/(0.135*323)))'
scale_ih = False

print('synapses = ',synapses)#will print: {branch_number: {location_on_branch: {active_synapses_number: array_of_stimulation_times_at_that_synapse}}}

'running the program'
cell, voltage_matrix = mr.dendrite_heatmap(synapses, cell_model = model, scale_ih = scale_ih, active_loc = active_loc)

'to confirm gIh density, comment in'
# for seg in cell.dend:
#     print(seg.gIhbar_Ih)


