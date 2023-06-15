import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)
import plot_functions as pf
import matplotlib.pyplot as plt
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

freq = 200
step = 1000/freq
start_1 = 10#time of first stimulation after start of simulation    
'control where along the branch to place the synapses, as a fraction of the branch length'
syn_locations = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
active_syn_loc = 0.5#location of active synapses, as a fraction of branch length
active_loc = (0, active_syn_loc)

fig, axes = plt.subplots()#figsize = (15,10))
plt.xticks([])
plt.yticks([])
for i in ['top','bottom','right','left']:
    axes.spines[i].set_visible(False)
rows,cols = len(act_vec), len(num_of_active_syns),
subplot_mat = np.matrix([list(range(1 + cols * i, 1 + cols * (i + 1))) for i in range(rows)])
plot_order_vec = []
[plot_order_vec.extend(np.array(i)[0]) for i in np.transpose(subplot_mat)]
weight = 1
none = {weight:[]}
count = 0 
for i in range(len(num_of_active_syns)):
    num_of_traces = (len(num_of_active_syns), highlighted_trace, i)#number of traces, trace we want to highlight, current trace number to plot
    active_synapses = num_of_active_syns[i]
    syn_times = np.arange(start_1,start_1+num_of_stims*step,step)
    stims = {weight+i:syn_times for i in range(active_synapses)}#only works because weights irrelevant. can change to weight+i with i = 0.00001...
    for j in range(len(act_vec)):
        synapses = {0:{}} # for many synapse locations on single branch, next two lines as well
        for loc in syn_locations:
            synapses[0][loc] = {weight:[]}           
        synapses[act_vec[j]][active_loc[1]] = stims
        print('synapses = ',synapses)#will print: {branch_number: {location_on_branch: {active_synapses_number: array_of_stimulation_times_at_that_synapse}}}
        cell,fig = cpmr.calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt, traces = False, 
                                       active_loc = active_loc, fig = fig, num_of_traces = num_of_traces, theta_d_GB = theta_d_GB, theta_p_GB = theta_p_GB)
        fig_num = j + count
        fig_num = plot_order_vec[fig_num]
        
        if fig_num <= cols:
            col_title = num_of_active_syns[fig_num-1]
        if fig_num == 1:
            plt.suptitle('dt: '+str(dt)+', # of stims: '+str(num_of_stims)+', freq: '
                         +str(1000/step)+'Hz, theta_d: '+ str(cell.e_synapses[0].theta_d_GB) + ', theta_p: ' + str(cell.e_synapses[0].theta_p_GB))
        fig, axes = pf.plot_plasticity_tree(cell, fig, axes, fig_num, rows, cols, col_title = col_title, syn_size = 350, soma_size = 600)
        col_title = None
    count += rows



















