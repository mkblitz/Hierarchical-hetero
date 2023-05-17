import matplotlib.pyplot as plt
import numpy as np
import ball_and_stick as bas
import ball_and_many_sticks as bams
import plot_functions as pf
import utils
import L5PC
from Calcium_plasticity_model_run import calcium_plasticity as cp 

def colored_branches(num_of_active_syns, act_vec, syn_times, num_of_branches, branch_L, theta_d_GB = 0.5, theta_p_GB = 1,active_synapse = (0,0.5), 
                     sim_time = 1000, dt = 0.5, num_of_stims = 1, freq = 200):
    
    syn_step = 1000/freq
    fig, axes = plt.subplots()#figsize = (15,10), dpi = 300)
    plt.xticks([])
    plt.yticks([])
    for i in ['top','bottom','right','left']:
        axes.spines[i].set_visible(False)
    rows,cols = len(act_vec), len(num_of_active_syns),
    subplot_mat = np.matrix([list(range(1 + cols * i, 1 + cols * (i + 1))) for i in range(rows)])
    plot_order_vec = []
    [plot_order_vec.extend(np.array(i)[0]) for i in np.transpose(subplot_mat)]
    count = 0
    
    for i in range(len(num_of_active_syns)):
        active_synapses = num_of_active_syns[i]
        stims = {1+i:syn_times for i in range(active_synapses)}
        for j in range(len(act_vec)):

            none = {1:[]}
            sec_fracs = np.round(np.arange(0.1,1.01,0.2),3)
            synapses = {sec_num:{sec_frac:none for sec_frac in sec_fracs} for sec_num in range(num_of_branches)}
            
            synapses[act_vec[j]][active_synapse[1]] = stims

            cell,_ = cp(num_of_branches, synapses, sim_time = sim_time, cell_model = 'bams', dt = dt, dend_L = branch_L, theta_d_GB = theta_d_GB,
                        theta_p_GB = theta_p_GB)

            fig_num = count + j
            fig_num = plot_order_vec[fig_num]
    
            if fig_num <= cols:
                col_title = num_of_active_syns[fig_num-1]
            if fig_num == 1:
                plt.suptitle('dt: '+str(dt)+', # of stims: '+str(num_of_stims)+', freq: '
                             +str(1000/syn_step)+'Hz, theta_d: '+ str(cell.e_synapses[0].theta_d_GB) + ', theta_p: ' + str(cell.e_synapses[0].theta_p_GB))
            fig, axes = pf.plot_colored_branches(cell, fig, axes, fig_num, rows, cols, col_title = col_title,)# active_syn_index = active_synapse_index)
            col_title = None
        count += rows
#             if plot_many_dendrite_subplot:
#                 count += rows
    return cell

if __name__ == "__main__":
    
    start_1 = 10
    syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)
    num_of_active_syns = [20,30,40,50]
    act_vec = [0,1,3,7,]#on which segments the activation should happen 
    num_of_branches = 15#1,3,7,15,31,63,127
    branch_L = 50
    cell = colored_branches(num_of_active_syns, act_vec, syn_times, num_of_branches, branch_L)
    print('end')
    