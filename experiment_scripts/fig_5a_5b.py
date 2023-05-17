import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)
import plot_functions as pf
import model_run as mr

num_of_stims = 1#number of stimulations at every active synapse
freq = 200
syn_step = 1000/freq
start_1 = 10#time of first stimulation after start of simulation    
syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)
sim_time = 200

voltage_matrices, prox_syn_voltage_traces, mid_syn_voltage_traces, dist_syn_voltage_traces, inactive_syn_indices = [],[],[],[],[]

count = 0
prox_loc, mid_loc, dist_loc = 0.3,0.65,1
locations = [prox_loc, mid_loc, dist_loc]

'to control number of proximal and distal synapses'
num_of_prox_syns_i = 20
num_of_dist_syns_i = 20
num_of_prox_syns_f = num_of_prox_syns_i + num_of_dist_syns_i
num_of_dist_syns_f = num_of_prox_syns_i + num_of_dist_syns_i

inactive_syn_indices = [num_of_prox_syns_i,num_of_prox_syns_f,1]

for subplot_num in range(3):
    synapses = {0:{}}
    if count == 0:
        stims = {1+i:syn_times for i in range(num_of_prox_syns_f)}
        synapses[0][locations[0]] = stims
        synapses[0][locations[1]] = {1:[]}
        synapses[0][locations[2]] = {1:[]}
    if count == 1:
        stims_a = {1+i:syn_times for i in range(num_of_prox_syns_i)}
        stims_b = {1+i:syn_times for i in range(num_of_dist_syns_i)}
        synapses[0][locations[0]] = stims_a
        synapses[0][locations[1]] = {1:[]}
        synapses[0][locations[2]] = stims_b
    if count == 2:
        stims = {1+i:syn_times for i in range(num_of_dist_syns_f)}
        synapses[0][locations[0]] = {1:[]}
        synapses[0][locations[1]] = {1:[]}
        synapses[0][locations[2]] = stims        
      
    print('experiment {}, synapses ='.format(count), synapses)#will print: {branch_number: {location_on_branch: {active_synapses_number: array_of_stimulation_times_at_that_synapse}}}
      
    cell, vm = mr.dendrite_heatmap(synapses, plot = False, sim_time = sim_time)
    voltage_matrices.append(vm)
              
    prox_syn_voltage_traces.append(cell.vsyn[0])
    mid_syn_voltage_traces.append(cell.vsyn[inactive_syn_indices[subplot_num]])
    dist_syn_voltage_traces.append(cell.vsyn[-1])

    count += 1

syn_voltage_traces = [prox_syn_voltage_traces, mid_syn_voltage_traces, dist_syn_voltage_traces]    
   
pf.plot_dendrite_heatmaps_fig_5(cell, voltage_matrices, syn_voltage_traces, locations, num_of_prox_syns_i, num_of_dist_syns_i, num_of_prox_syns_f,
                                 num_of_dist_syns_f, )
