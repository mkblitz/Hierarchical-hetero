import matplotlib.pyplot as plt
import numpy as np
# import sys
# from neuron import h
import ball_and_stick as bas
import ball_and_many_sticks as bams
# import pickle as pkl
import plot_functions as pf
import utils
import L5PC
# plt.ion()
# import time
# import winsound
# t_initial = time.time()
from Calcium_plasticity_model_run import calcium_plasticity as cp 

print('start')
def colored_branches():
    num_of_stims = 1
#     num_of_active_syns = 15
    freq = 200
    syn_step = 1000/freq
    start_1 = 10
    syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)
#     stims = {1+i:syn_times for i in range(num_of_active_syns)}
    dt = 0.5
    sim_time = 1000
    spines = True
    
    dend_L = 50
    num_of_active_syns = [1,20,40,60]
    act_vec = [0,1,3,7,]#on which segments the activation should happen 
    num_of_sections = 15#1,3,7,15,31,63,127
    weight = 1
    none = {weight:[]}
    sec_fracs = np.round(np.arange(0.1,1.01,0.2),3)
#     sec_fracs = np.round(np.arange(0.25,1.01,0.5),3)
#     synapses = {sec_num:{sec_frac:none for sec_frac in sec_fracs} for sec_num in range(num_of_sections)}
    active_synapse = (0,0.5)#branch number, location fraction
    
    fig, axes = plt.subplots(figsize = (15,10), dpi = 250)
    plt.xticks([])
    plt.yticks([])
    for i in ['top','bottom','right','left']:
        axes.spines[i].set_visible(False)
    rows,cols = len(act_vec), len(num_of_active_syns),
    subplot_mat = np.matrix([list(range(1 + cols * i, 1 + cols * (i + 1))) for i in range(rows)])
    plot_order_vec = []
    [plot_order_vec.extend(np.array(i)[0]) for i in np.transpose(subplot_mat)]
    count = 0
#     print(synapses)

    for i in range(len(num_of_active_syns)):
        active_synapses = num_of_active_syns[i]
        stims = {weight+i:syn_times for i in range(active_synapses)}
        for j in range(len(act_vec)):
#             active_synapse_index = len(sec_fracs)*(act_vec[j] + 1)
#             print(active_synapse_index)
            synapses = {sec_num:{sec_frac:none for sec_frac in sec_fracs} for sec_num in range(num_of_sections)}
            synapses[act_vec[j]][active_synapse[1]] = stims
#             print('synapses ', synapses)

            cell,_ = cp(num_of_sections, synapses, sim_time = sim_time, cell_model = 'bams', dt = dt, dend_L = dend_L)

            fig_num = count + j
            fig_num = plot_order_vec[fig_num]
    
            if fig_num <= cols:
                col_title = num_of_active_syns[fig_num-1]
#                 col_title = num_of_active_syns
            if fig_num == 1:
                plt.suptitle('dt: '+str(dt)+', # of stims: '+str(num_of_stims)+', freq: '
                             +str(1000/syn_step)+'Hz, theta_d: '+ str(cell.e_synapses[0].theta_d_GB) + ', theta_p: ' + str(cell.e_synapses[0].theta_p_GB))
            fig, axes = pf.plot_colored_branches(cell, fig, axes, fig_num, rows, cols, col_title = col_title,)# active_syn_index = active_synapse_index)
            col_title = None
        count += rows
#             if plot_many_dendrite_subplot:
#                 count += rows
#     fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Calcium in Dendrites - Paper/Figures/Figure 4/fig.png')
    return cell

cell = colored_branches()
# print(len(cell.e_synapses))
print('end')
# for v in cell.rho_GB:
#     print(v[-1])
#     
# plt.plot(cell.t,cell.vsyn[200])
# plt.plot(cell.t,cell.effcai_GB[5])
# plt.plot(cell.t,cell.effcai_GB[200])
    
    
    