'''
proximal and distal active synapse clusters, and central inactive synapses. checking peak calcium and plasticity at all three locations for different 
strength activations.
'''
import ball_and_many_sticks as bams
import numpy as np
import plot_functions as pf
import matplotlib.pyplot as plt
import utils
from matplotlib import cm 
import matplotlib.gridspec as gridspec
import time
import matplotlib.colors as colors
import math 
import winsound
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

sim_time = 200
dt = 0.1
weight = 1

'{syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}'
start_1,start_2,start_3,start_4,start_5 = 10,4500,1500,2250,3000

def calcium_plasticity(num_of_sections, synapses, sim_time = 100, cell_model = 'bams', nmda_ratio = 2.38, excitatory = True, num_of_rec_locations = 30,
                       dend_L = 200, dt = 0.025, synapse_type = 'glu syn', traces = False, menorah = False, active_loc = (0, 0.5),
                                fig = None, num_of_traces = None, theta_d_GB = 0.2, theta_p_GB = 0.4):
        '''
        synapses: {syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}
        syn_locations: [syn1, syn2,...]
        syn_times, weights: [[syn1_1,syn1_2,...],[syn2_1,syn2_2,...]]
        '''
        dend_sec_L = dend_L
        if cell_model == 'bams':
            cell = bams.many_sticks_model(num_of_sections, dend_L = dend_sec_L,)
        
        cell.synapses = synapses
        cell.active_synapses = {}
        for branch in cell.synapses:
            for loc in cell.synapses[branch]:
                if len(cell.synapses[branch][loc][1])>0:
                    cell.active_synapses[branch] = loc
        rec_locations = np.linspace(0, 1, num_of_rec_locations)
        for syn_sec_num in list(synapses):
            for syn_loc in list(synapses[syn_sec_num]):
                for syn_weight in list(synapses[syn_sec_num][syn_loc]):
                    if excitatory:
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_sec_num][syn_loc][syn_weight]], dend_sec_num = syn_sec_num,
                                                        NMDA_ratio = nmda_ratio, synapse_type = synapse_type, theta_d_GB = theta_d_GB, theta_p_GB = theta_p_GB)
    #     rec_vecs = ['i','cai_CR','g_NMDA','g_AMPA','gmax_d_AMPA','ica_NMDA','ica_VDCC']
        rec_vecs = ['vsyn','ica_NMDA','ica_VDCC','effcai_GB','rho_GB']#,'m_VDCC', 'h_VDCC','gca_bar_VDCC','effcai_GB','gca_VDCC','gmax_AMPA'
        cell.run_sim(rec_locations, sim_time, dt, recording_vecs = rec_vecs, v_init =  None)
    #     pf.plot_multiple_dendrite_heatmaps(cell, annotate = True, dt = dt, sim_time = sim_time)
        save_fig = False
        if traces:
            fig = pf.plot_conductances(cell,rec_vecs, save_fig, protocol, active_loc = active_loc, fig = fig, num_of_traces = num_of_traces)
        if menorah:
            pf.plot_plasticity_tree(cell)    
        return cell, fig
    
def run_experiment(max_synapses, syn_branches, num_of_sections = 15, branch_L = 50, theta_d_GB = 0.2, theta_p_GB = 0.4):  
    num_of_stims = 1
    active_synapses_vec = np.arange(0, max_synapses+1, 5)#y axis 

    syn_locs = [0.5,0.5,0.5]
    
    num_of_experiments = 1
    num_of_syn_locs = len(syn_locs)#number of columns
    
    freq = 200
    syn_step = 1000//freq
    syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
    mean_max = 'max'
    fig_title = mean_max +' eff_cai + rho , num_of_stims: '+str(num_of_stims) +', freq: ' +str(freq) +'hz'
    fig, axes = plt.subplots(2,3, figsize = (14,8))
    cmap_rho = cm.get_cmap('bwr', 3)
    num_of_ticklabels = 2# number of jumps between xticklabels. The higher the number, the fewer tick labes
    eff_cai_matrices_1, eff_cai_matrices_inact, eff_cai_matrices_2 = [],[],[]
    rho_matrices_1, rho_matrices_inact, rho_matrices_2 = [],[],[]
    
    for experiment_num in range(num_of_experiments):
        syn1_location = syn_locs[0]
        inactive_syn_loc = syn_locs[1]
        syn2_location = syn_locs[2]
        
        eff_cai_mat_1, eff_cai_mat_inact, eff_cai_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
        rho_mat_1, rho_mat_inact, rho_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
    
        for active_syn_num1 in range(len(active_synapses_vec)):
            col = active_syn_num1
            active_synapses1 = active_synapses_vec[active_syn_num1]
            print('row: '+str(active_syn_num1))
            for active_syn_num2 in range(len(active_synapses_vec)):
                
                row = len(active_synapses_vec) - active_syn_num2 - 1
                active_synapses2 = active_synapses_vec[active_syn_num2]
                stims_1 = {weight+i:syn_times for i in range(active_synapses1)}
                stims_2 = {weight+i:syn_times for i in range(active_synapses2)}
                if active_synapses1 == 0: stims_1 = {weight:[]}
                if active_synapses2 == 0: stims_2 = {weight:[]}
    
                synapses = {syn_branches[0]:{syn1_location:stims_1}, syn_branches[1]:{inactive_syn_loc:{1:[]}}, syn_branches[2]:{syn2_location:stims_2}}
    
                print(synapses)
                cell,_ = calcium_plasticity(num_of_sections, synapses, dend_L = branch_L, sim_time = sim_time, dt = dt, traces = False,
                                           active_loc = (0, syn1_location), theta_d_GB = theta_d_GB, theta_p_GB = theta_p_GB)
                '''ONLY WORKS FOR 2 SYNAPSE LOCATIONS!!!! For more, need to work out mapping to relevant rho_GB vecs'''
    #             if mean_max == 'mean':
    #                 eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.mean(utils.npa(cell.effcai_GB[0])), np.mean(utils.npa(cell.effcai_GB[-1]))
                print(len(list(synapses[list(synapses)[0]][syn1_location])))
    
                eff_cai_mat_1[row, col], eff_cai_mat_inact[row,col], eff_cai_mat_2[row, col] = np.max(utils.npa(cell.effcai_GB[0])), np.max(utils.npa(cell.effcai_GB[len(list(synapses[list(synapses)[0]][syn1_location]))])), np.max(utils.npa(cell.effcai_GB[-1]))
    
                '''add lines for extra matrices from above'''
                rho_mat_1[row, col], rho_mat_inact[row, col], rho_mat_2[row, col] = np.heaviside(cell.rho_GB[0][-1]-0.5,0.5), np.heaviside(cell.rho_GB[len(list(synapses[list(synapses)[0]][syn1_location]))][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
    
        eff_cai_matrices_1.append(eff_cai_mat_1), eff_cai_matrices_inact.append(eff_cai_mat_inact), eff_cai_matrices_2.append(eff_cai_mat_2)
        rho_matrices_1.append(rho_mat_1), rho_matrices_inact.append(rho_mat_inact), rho_matrices_2.append(rho_mat_2)
    
    min_val = np.min((np.min(eff_cai_matrices_1), np.min(eff_cai_matrices_inact), np.min(eff_cai_matrices_2)))
    max_val = np.max((np.max(eff_cai_matrices_1), np.max(eff_cai_matrices_inact), np.max(eff_cai_matrices_2)))
    vmin = np.round(min_val,3)
    vmax = np.round(max_val,3)
     
    gamma = 1
    tick1,tick2,tick3,tick4,tick5 = 0,vmax/4,2*vmax/4,3*vmax/4,vmax
    
    cmap = colors.LinearSegmentedColormap.from_list('some name', [(tick1/vmax, 'white'),(tick2/vmax,'yellow'),(tick3/vmax, 'orange'),
                                                                  (tick4/vmax, 'red'),(tick5/vmax,'black')],
                                                    gamma = gamma)
    cmap = 'hot_r'
    print(tick1,tick2,tick3,tick4,tick5,)
    for experiment_num in range(num_of_experiments):
        syn1_location = syn_locs[0]
        inactive_syn_loc = syn_locs[1]
        syn2_location = syn_locs[2]
    
        'synapse 1'
        ax1 = axes[0,0]
        ax2 = axes[1,0]
        'inactive synapse'
        ax3 = axes[0,1]
        ax4 = axes[1,1]
        'synapse 2'
        ax5 = axes[0,2]
        ax6 = axes[1,2]
    
        ax1.set_title('Proximal synapse')# - loc: {:.0f}'.format(syn1_location*dend_L))
    
        fig1 = ax1.imshow(eff_cai_matrices_1[experiment_num], cmap = cmap, )#vmin = 0, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
        fig2 = ax2.imshow(rho_matrices_1[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
        ax2.set_xlabel('# active syns (1)')
        plt.colorbar(fig1, ax = ax1, fraction=0.046, pad=0.02, orientation = 'horizontal')
    
        ax3.set_title('Inactive synapse')# - loc: {:.0f}'.format(syn2_location*dend_L))
        fig3 = ax3.imshow(eff_cai_matrices_inact[experiment_num], cmap = cmap, )#vmin = 0, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
        fig4 = ax4.imshow(rho_matrices_inact[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
        ax4.set_xlabel('# active syns (1)')
        plt.colorbar(fig3, ax = ax3, fraction=0.046, pad=0.02, orientation = 'horizontal')
        ax5.set_title('Distal synapse')
        fig5 = ax5.imshow(eff_cai_matrices_2[experiment_num], cmap = cmap, )#vmin = 0, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
        fig6 = ax6.imshow(rho_matrices_2[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
        ax6.set_xlabel('# active syns (1)')
        plt.colorbar(fig5, ax = ax5, fraction=0.046, pad=0.02, orientation = 'horizontal')

        ax1.set_yticks(np.arange(len(active_synapses_vec)))
        ax1.set_yticklabels(active_synapses_vec[::-1])
        ax1.set_ylabel('eff_ca\n# active syns (2)')
    
        ax2.set_yticks(np.arange(len(active_synapses_vec)))
        ax2.set_yticklabels(active_synapses_vec[::-1])
        ax2.set_ylabel('rho\n# active syns (2)')
        
        fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB) + '\n branches = '+str(syn_branches))
        
        ax3.set_yticks(np.arange(len(active_synapses_vec)))
        ax3.set_yticklabels([])
        ax4.set_yticks(np.arange(len(active_synapses_vec)))
        ax4.set_yticklabels([])
        ax5.set_yticks(np.arange(len(active_synapses_vec)))
        ax5.set_yticklabels([])
        ax6.set_yticks(np.arange(len(active_synapses_vec)))
        ax6.set_yticklabels([])

        ax1.set_xticks(np.arange(len(active_synapses_vec)))
        ax1.set_xticklabels([])
        
        ax2.set_xticks(np.arange(len(active_synapses_vec)))
        ax2.set_xticklabels(active_synapses_vec)
        
        ax3.set_xticks(np.arange(len(active_synapses_vec)))
        ax3.set_xticklabels([])
    
        ax4.set_xticks(np.arange(len(active_synapses_vec)))
        ax4.set_xticklabels(active_synapses_vec)
        
        ax5.set_xticks(np.arange(len(active_synapses_vec)))
        ax5.set_xticklabels([])
    
        ax6.set_xticks(np.arange(len(active_synapses_vec)))
        ax6.set_xticklabels(active_synapses_vec)
        print(np.shape(np.squeeze(axes)))
        kws = dict(ha='center', va='center', color='deepskyblue', fontsize=15)
        for axi in axes:
            for ax in axi:
                for i in range(len(active_synapses_vec)):
                    for j in range(len(active_synapses_vec)):
                            ax.text(j, i, active_synapses_vec[-1-i]+active_synapses_vec[j], **kws)
    
    winsound.Beep(1000,200)
    print('end')
        
    