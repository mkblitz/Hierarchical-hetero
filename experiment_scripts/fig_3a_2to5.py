import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)
import Calcium_plasticity_model_run as cpmr

sim_time = 1000
dt = 0.5

num_of_stims = 1#number of stimulations at every active synapse
dend_L = 200#length of branch
num_of_active_syns = [27]#number of active synapses
highlighted_trace = num_of_active_syns[0]

act_vec = [0]
num_of_sections = 1

theta_d_GB = 0.5#depression threshold
theta_p_GB = 1#potentiation threshold

"add 'rho_GB' to see the plasticity. other options can be found in Calcium_plasticity_model_run -> calcium_plasticity"
rec_vecs = ['vsyn','ica_NMDA','ica_VDCC','effcai_GB',]

freq = 200
step = 1000/freq
start_1 = 10#time of first stimulation after start of simulation    
syn_locations = [0.3,0.45,0.5,0.7]#where along the dendrite to place the synapses, as a fraction of the branch length
active_syn_loc = 0.5#location of active synapses, as a fraction of branch length
active_loc = (0, active_syn_loc)
fig = None
for i in range(len(num_of_active_syns)):
    num_of_traces = (len(num_of_active_syns), highlighted_trace, i)#number of traces, trace we want to highlight, current trace number to plot
    active_synapses = num_of_active_syns[i]
    syn_times = np.arange(start_1,start_1+num_of_stims*step,step)
    stims = {1+i:syn_times for i in range(active_synapses)}#only works because weights irrelevant. can change to weight+i with i = 0.00001...
    for j in range(len(act_vec)):
        synapses = {0:{}} # for many synapse locations on single branch, next two lines as well
        for loc in syn_locations:
            synapses[0][loc] = {1:[]}           
        synapses[act_vec[j]][active_loc[1]] = stims
        print('synapses = ',synapses)#will print: {branch_number: {location_on_branch: {active_synapses_number: array_of_stimulation_times_at_that_synapse}}}
        cell,fig = cpmr.calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt, traces = True, 
                                       active_loc = active_loc, fig = fig, num_of_traces = num_of_traces, theta_d_GB = theta_d_GB, theta_p_GB = theta_p_GB,
                                       rec_vecs = rec_vecs)
        
        