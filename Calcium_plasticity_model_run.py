import ball_and_many_sticks as bams
import numpy as np
import plot_functions as pf
import matplotlib.pyplot as plt
import utils
from matplotlib import cm 
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors

sim_time = 1000
dt = 0.5
weight = 1
# step = 10
'{syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}'
start_1,start_2,start_3,start_4,start_5 = 210,4500,1500,2250,3000

def calcium_plasticity(num_of_sections, synapses, sim_time = 100, cell_model = 'bams', nmda_ratio = 2.38, excitatory = True, num_of_rec_locations = 30,
                       dend_L = 200, dt = 0.025, synapse_type = 'glu syn', traces = False, menorah = False, active_loc = (0, 0.5),
                                fig = None, num_of_traces = None, theta_d_GB = 0.5, theta_p_GB = 1,
                                rec_vecs = ['vsyn','ica_NMDA','ica_VDCC','effcai_GB','rho_GB']):
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
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_sec_num][syn_loc][syn_weight]], 
                                                        dend_sec_num = syn_sec_num, NMDA_ratio = nmda_ratio, synapse_type = synapse_type, 
                                                        theta_d_GB = theta_d_GB, theta_p_GB = theta_p_GB)
#         rec_vecs options: 'vsyn','ica_NMDA','ica_VDCC','effcai_GB','rho_GB','m_VDCC', 'h_VDCC','gca_bar_VDCC','effcai_GB','gca_VDCC','gmax_AMPA',
#        'i','cai_CR','g_NMDA','g_AMPA','gmax_d_AMPA'
        cell.run_sim(rec_locations, sim_time, dt, recording_vecs = rec_vecs, v_init =  None)
    #     pf.plot_multiple_dendrite_heatmaps(cell, annotate = True, dt = dt, sim_time = sim_time)
        save_fig = False
        if traces:
            fig = pf.plot_conductances(cell,rec_vecs, save_fig, protocol, active_loc = active_loc, fig = fig, num_of_traces = num_of_traces)
        if menorah:
            pf.plot_plasticity_tree(cell)    
        return cell, fig

protocol = 'multiple_branch_multiple_syn'#'multiple_branch_multiple_syn','multiple_branch_single_syn','2_syn_heatmap','eff_cai_heatmap',
#eff_cai_heatmap_activeSyns_vs_dt,eff_cai_heatmap_activeSyns_vs_activeSyns
if __name__ == "__main__":
    print('start')
    if protocol in ['multiple_branch_multiple_syn','multiple_branch_single_syn']:
        num_of_stims = 1
        if protocol == 'multiple_branch_single_syn': 
            dend_L = 50
            num_of_active_syns = [30,50,70,90,]
            act_vec = [0,1,3,7,]#on which segments the activation should happen 
            num_of_sections = 1#1,3,7,15,31,63,127  
        elif protocol == 'multiple_branch_multiple_syn': 
            dend_L = 200
            num_of_active_syns = [27]
            highlighted_trace = num_of_active_syns[0]
            
            act_vec = [0]
            num_of_sections = 1
        freq = 200
        step = 1000/freq
        start_1 = 10    
        plot_many_dendrite_subplot = True #turn this on to get NxM subplots, each one with different active syn location and different num of active synapses
        traces = False
        if plot_many_dendrite_subplot:
            fig, axes = plt.subplots(figsize = (15,10))
            plt.xticks([])
            plt.yticks([])
            for i in ['top','bottom','right','left']:
                axes.spines[i].set_visible(False)
            rows,cols = len(act_vec), len(num_of_active_syns),
            subplot_mat = np.matrix([list(range(1 + cols * i, 1 + cols * (i + 1))) for i in range(rows)])
            plot_order_vec = []
            [plot_order_vec.extend(np.array(i)[0]) for i in np.transpose(subplot_mat)]
        syn_locations = [0.3,0.45,0.5,0.7]
        none = {weight:[]}
        active_syn_loc = 0.5
        active_loc = (0, active_syn_loc)
        count = 0 
        if traces:
            fig = None
        
        for i in range(len(num_of_active_syns)):
            num_of_traces = (len(num_of_active_syns), highlighted_trace, i)#number of traces, trace we want to highlight, current trace number to plot
            active_synapses = num_of_active_syns[i]
            syn_times = np.arange(start_1,start_1+num_of_stims*step,step)
            stims = {weight+i:syn_times for i in range(active_synapses)}#only works because weights irrelevant. can change to weight+i with i = 0.00001...
            stims_2 = {weight+i:syn_times+20 for i in range(active_synapses)}
            for j in range(len(act_vec)):
                if protocol == 'multiple_branch_single_syn':
                    synapses = {k:{active_loc[1]:none} for k in range(num_of_sections)} #for single syn loc on each branch
                elif protocol == 'multiple_branch_multiple_syn':
                    synapses = {0:{}} # for many synapse locations on single branch, next two lines as well
                    for loc in syn_locations:#[0.3,0.4,0.5,0.6,0.7]:
                        synapses[0][loc] = {weight:[]}           
                synapses[act_vec[j]][active_loc[1]] = stims
                cell,fig = calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt, traces = traces, 
                                               active_loc = active_loc, fig = fig, num_of_traces = num_of_traces, )# menorah = menorah)
                if plot_many_dendrite_subplot:
                    fig_num = j + count
                    fig_num = plot_order_vec[fig_num]
                    
                    if fig_num <= cols:
                        col_title = num_of_active_syns[fig_num-1]
                    if fig_num == 1:
                        plt.suptitle('dt: '+str(dt)+', # of stims: '+str(num_of_stims)+', freq: '
                                     +str(1000/step)+'Hz, theta_d: '+ str(cell.e_synapses[0].theta_d_GB) + ', theta_p: ' + str(cell.e_synapses[0].theta_p_GB))
                    fig, axes = pf.plot_plasticity_tree(cell, fig, axes, fig_num, rows, cols, col_title = col_title, syn_size = 350, soma_size = 600)
                    col_title = None
            if plot_many_dendrite_subplot:
                count += rows
    
    elif protocol == '2_syn_heatmap':
        num_of_sections = 1 
        dend_L = 200 
        dx_vec = np.float16(np.arange(-0.15,0.16,0.05))+0.0001#x axis
        dt_vec = np.linspace(-300,300,len(dx_vec))#y axis
        syn_1_locs = [0.4,0.5,0.6]
        num_of_syn1_locs = len(syn_1_locs)#number of columns
        stim_protocol = '1_1'#'1_0','1_1'
        variable = 'frequency'#'frequency','num_of_active_syns'
        if variable == 'num_of_active_syns':
            active_synapse_options = [6,8,10]#10]#at active location, how many active syns should be firing
            num_of_rows = len(active_synapse_options)
            num_of_stims = 1
            syn_step = 1#number of ms between consecutive stimulations
            syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
            fig_title = 'stim_protocol: ' + stim_protocol 
    
        elif variable == 'frequency':
            active_synapse_options = [1]
            active_synapses = active_synapse_options[0]
            num_of_stims = 2
            syn_steps = [500,333,250]
            num_of_rows = len(syn_steps)
            fig_title = 'stim_protocol: ' + stim_protocol +', num_active_syns: ' + str(num_of_stims)
            
        cmap = cm.get_cmap('coolwarm',)# 3)
        fig = plt.figure()
        outer = gridspec.GridSpec(num_of_rows, num_of_syn1_locs, )
        num_of_ticklabels = 2# number of jumps between xticklabels. The higher the number, the fewer tick labes
        
        for experiment_num in range(num_of_syn1_locs*num_of_rows):
            syn_1_loc = syn_1_locs[experiment_num%num_of_syn1_locs]
            if variable == 'num_of_active_syns':
                active_synapses = active_synapse_options[np.int0(np.floor(experiment_num/num_of_syn1_locs))]
            elif variable == 'frequency':
                syn_step = syn_steps[np.int0(np.floor(experiment_num/num_of_syn1_locs))]
                
                frequency = 1000//syn_step
            syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
            stims_1 = {weight+i:syn_times for i in range(active_synapses)}
            rho_mat_1, rho_mat_2 = np.zeros((len(dt_vec),len(dx_vec))), np.zeros((len(dt_vec),len(dx_vec)))
            
            for loc_num in range(len(dx_vec)):
                syn_2_loc = syn_1_loc + dx_vec[loc_num] 
                if syn_2_loc < 0: syn_2_loc = 0 #to make sure 0 <= syn 2 loc <=1
                elif syn_2_loc > 1: syn_2_loc = 1
                col = loc_num
                for delta_t_num in range(len(dt_vec)):
                    row = delta_t_num
                    stims_2 = {weight+i:syn_times + dt_vec[delta_t_num] for i in range(active_synapses)}
                    if stim_protocol == '1_0': synapses = {0:{syn_1_loc:stims_1, syn_2_loc:{weight:[]}}}#for 1 0
                    elif stim_protocol == '1_1': synapses= {0:{syn_1_loc:stims_1, syn_2_loc:stims_2}}# for 1 1
                    cell = calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt)
                    '''ONLY WORKS FOR 2 SYNAPSE LOCATIONS!!!! For more, need to work out mapping to relevant rho_GB vecs'''
                    rho_mat_1[row, col], rho_mat_2[row, col] = np.heaviside(cell.rho_GB[0][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
    #                 rho_mat_1[row, col], rho_mat_2[row, col] = cell.rho_GB[0][-1], cell.rho_GB[-1][-1]
            inner = gridspec.GridSpecFromSubplotSpec(1, 2,#creates inner axes 1x2
                            subplot_spec=outer[experiment_num],)
            
            ax = plt.Subplot(fig, inner[0])
            if experiment_num < num_of_syn1_locs:
                ax.set_title('Synapse 1 - loc: {}'.format(syn_1_loc*dend_L))
            ax.imshow(rho_mat_1, vmin = 0, vmax = 1, cmap = cmap)
            ax.set_xticks(np.arange(len(dx_vec))[::num_of_ticklabels])
            ax.set_xticklabels((dend_L*np.round(dx_vec, 2))[::num_of_ticklabels])
            ax.set_yticks(np.arange(len(dt_vec)))
            ax.set_yticklabels(dt_vec)
            if experiment_num >= num_of_syn1_locs*num_of_rows - num_of_syn1_locs:
                ax.set_xlabel('loc of syn 2 relative to syn 1')
            if experiment_num % num_of_syn1_locs == 0:
                if variable == 'num_of_active_syns':
                    ax.set_ylabel('# of active synapses: {}\ndt'.format(active_synapses))
                elif variable == 'frequency':
                    ax.set_ylabel('Frequency: {}hz\ndt'.format(frequency))
            else:
                ax.set_ylabel('dt')
            fig.add_subplot(ax)
            
            ax = plt.Subplot(fig, inner[1])
            if experiment_num < num_of_syn1_locs:
                ax.set_title('Synapse 2')
            ax.imshow(rho_mat_2, vmin = 0, vmax = 1, cmap = cmap)
            ax.set_xticks(np.arange(len(dx_vec))[::num_of_ticklabels])
            ax.set_xticklabels((dend_L*np.round(dx_vec, 2))[::num_of_ticklabels])
            ax.set_yticks(np.arange(len(dt_vec)))
            ax.set_yticklabels(dt_vec)
            if experiment_num >= num_of_syn1_locs*num_of_rows - num_of_syn1_locs:
                ax.set_xlabel('loc of syn 2 relative to syn 1')        
            fig.add_subplot(ax)
        
        fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB))
        cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), cax=cbar_ax)
        cbar.set_ticks([1/6, 3/6, 5/6])  
        cbar.set_ticklabels(['depression','no change','potentiation'])
    
    elif protocol == 'eff_cai_heatmap':
        num_of_sections = 1 
        dend_L = 200 
        num_of_stims = 3
        freq_vec = [1,10,100,1000]#np.arange(1,100,10)#x axis
        active_synapses_vec = [0,10,20,30]
        dt_vec = np.linspace(-1000,1000,5)#y axis
        syn_locs = [0.25,0.5,0.75]
        prox_loc_options = [0,0,0,1,1,2]
        dist_loc_options = [0,1,2,1,2,2]
        num_of_syn_locs = len(syn_locs)#number of columns
        stim_protocol = '1_1'#'1_0','1_1'
    
        eff_cai_fig_title = '<eff_cai>, stim_protocol: ' + stim_protocol +', num_of_stims: '+str(num_of_stims) # + ', num_active_syns: ' + str(active_synapses)
        rho_fig_title = '<rho>, stim_protocol: ' + stim_protocol +', num_of_stims: '+str(num_of_stims)
        
        trace_fig, trace_axis = plt.subplots(len(active_synapses_vec), len(active_synapses_vec), sharey = True)
        trace_fig.suptitle('Synapse 1, num_active_syns: ' + str(active_synapses_vec))
        trace_axis_right = []#list of second y axis for trace_axis
        
        trace_fig2, trace_axis2 = plt.subplots(len(active_synapses_vec), len(active_synapses_vec), sharey = True)
        trace_axis2_insets1, trace_axis2_insets2 = [], []
        trace_fig2.suptitle('Synapse 2, num_active_syns: ' + str(active_synapses_vec))
        
        cmap = cm.get_cmap('coolwarm',)# 3)
        eff_cai_fig = plt.figure()
        rho_fig = plt.figure()
        outer = gridspec.GridSpec(1, 6,)
        left1, bottom1, width1, height1 = [0.6, 0.63, 0.35, 0.35]#inset 1 coordinates 
        left2, bottom2, width2, height2 = [0.6, 0.23, 0.35, 0.35]#inset 2 coordinates 
        for experiment_num in range(6):
            
            trace_col = 0
            syn1_location = syn_locs[prox_loc_options[experiment_num]]
            syn2_location = syn_locs[dist_loc_options[experiment_num]]+0.0001

            rho_mat_1, rho_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
            eff_cai_mat_1, eff_cai_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
            subplot_num = 0
            for active_syn_num1 in range(len(active_synapses_vec)):
                col = active_syn_num1
                active_synapses1 = active_synapses_vec[active_syn_num1]
                trace_row = 0 
 
                for active_syn_num2 in range(len(active_synapses_vec)):
                    row = len(active_synapses_vec) - active_syn_num2 - 1
                    active_synapses2 = active_synapses_vec[active_syn_num2]
                    freq = 100
                    syn_step = 1000//freq
                    syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
                    stims_1 = {weight+i:syn_times for i in range(active_synapses1)}
                    if active_synapses1 == 0:
                        stims_1 = {weight:[]}
                    stims_2 = {weight+i:syn_times for i in range(active_synapses2)}
                    if active_synapses2 == 0:
                        stims_2 = {weight:[]}
                    if stim_protocol == '1_0': synapses = {0:{syn1_location:stims_1, syn2_location:{weight:[]}}}#for 1 0
                    elif stim_protocol == '1_1': synapses= {0:{syn1_location:stims_1, syn2_location:stims_2}}# for 1 1
                    cell = calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt)
                    '''ONLY WORKS FOR 2 SYNAPSE LOCATIONS!!!! For more, need to work out mapping to relevant rho_GB vecs'''
                    eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.mean(utils.npa(cell.effcai_GB[0])), np.mean(utils.npa(cell.effcai_GB[-1]))
                    rho_mat_1[row, col], rho_mat_2[row, col] = np.heaviside(cell.rho_GB[0][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
                    trace_row += 1
                    subplot_num += 1
                trace_col += 1
            inner = gridspec.GridSpecFromSubplotSpec(2, 1,#creates inner axes 2x1
                            subplot_spec=outer[experiment_num],)
            
            ax_eff_cai = plt.Subplot(eff_cai_fig, inner[0])
            ax_eff_cai.set_title('Synapse 1 - loc: {}'.format(syn1_location*dend_L))
            ax_eff_cai.imshow(eff_cai_mat_1, cmap = cmap, vmin = 0, vmax = 0.4,)
            ax_eff_cai.set_xlabel('syn 1')
            
            ax_rho = plt.Subplot(rho_fig, inner[0])
            ax_rho.set_title('Synapse 1 - loc: {}'.format(syn1_location*dend_L))
            ax_rho.imshow(rho_mat_1, cmap = cmap, vmin = 0, vmax = 1,)
            ax_rho.set_xlabel('syn 1')
            
            if experiment_num == 0:
                ax_eff_cai.set_yticks(np.arange(len(active_synapses_vec)))
                ax_eff_cai.set_yticklabels(active_synapses_vec[::-1])
                ax_eff_cai.set_ylabel('syn 2')
                eff_cai_fig.suptitle(eff_cai_fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB))
    
                ax_rho.set_yticks(np.arange(len(active_synapses_vec)))
                ax_rho.set_yticklabels(active_synapses_vec[::-1])
                ax_rho.set_ylabel('syn 2')
                rho_fig.suptitle(rho_fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB))
                
                cbar_ax_eff_cai = eff_cai_fig.add_axes([0.92, 0.15, 0.01, 0.7])
                cbar_eff_cai = eff_cai_fig.colorbar(cm.ScalarMappable(cmap=cmap), cax=cbar_ax_eff_cai)
                cbar_eff_cai.set_ticks([1/6, 3/6, 5/6]) 
                cbar_eff_cai.set_ticklabels([])
                
                cbar_ax_rho = rho_fig.add_axes([0.92, 0.15, 0.01, 0.7])
                cbar_rho = rho_fig.colorbar(cm.ScalarMappable(cmap=cmap), cax=cbar_ax_rho)
                cbar_rho.set_ticks([1/6, 3/6, 5/6])  
                cbar_rho.set_ticklabels([])
            #     cbar_rho.set_ticklabels(['depression','no change','potentiation'])
            else:
                ax_eff_cai.set_yticklabels([])
                ax_rho.set_yticklabels([])
            ax_eff_cai.set_xticks(np.arange(len(active_synapses_vec)))
            ax_eff_cai.set_xticklabels(active_synapses_vec)
            
            ax_rho.set_xticks(np.arange(len(active_synapses_vec)))
            ax_rho.set_xticklabels(active_synapses_vec)

            eff_cai_fig.add_subplot(ax_eff_cai)
            rho_fig.add_subplot(ax_rho)
            
            ax_eff_cai = plt.Subplot(eff_cai_fig, inner[1])
            ax_eff_cai.set_title('Synapse 2 - loc: {:.0f}um'.format(syn2_location*dend_L))
            ax_eff_cai.imshow(eff_cai_mat_2, cmap = cmap, vmin = 0, vmax = 0.4)
            ax_eff_cai.set_xlabel('syn 1')
            
            ax_rho = plt.Subplot(rho_fig, inner[1])
            ax_rho.set_title('Synapse 2 - loc: {:.0f}um'.format(syn2_location*dend_L))
            ax_rho.imshow(rho_mat_2, cmap = cmap, vmin = 0, vmax = 1)
            ax_rho.set_xlabel('syn 1')
            
            ax_eff_cai.set_xticks(np.arange(len(active_synapses_vec)))
            ax_eff_cai.set_xticklabels(active_synapses_vec)
            
            ax_rho.set_xticks(np.arange(len(active_synapses_vec)))
            ax_rho.set_xticklabels(active_synapses_vec)
            
            if experiment_num == 0:
                ax_eff_cai.set_yticks(np.arange(len(active_synapses_vec)))
                ax_eff_cai.set_yticklabels(active_synapses_vec[::-1])
                ax_eff_cai.set_ylabel('syn 2')
                
                ax_rho.set_yticks(np.arange(len(active_synapses_vec)))
                ax_rho.set_yticklabels(active_synapses_vec[::-1])
                ax_rho.set_ylabel('syn 2')
            else:
                ax_eff_cai.set_yticklabels([])
                
                ax_rho.set_yticklabels([])
    
            eff_cai_fig.add_subplot(ax_eff_cai)
            
            rho_fig.add_subplot(ax_rho)
             
    elif protocol == 'eff_cai_heatmap_activeSyns_vs_dt':
        num_of_sections = 1 
        dend_L = 200 
        num_of_stims = 1
        active_synapses_vec = np.arange(1, 23, 3)#y axis
        start_1 = 210
        dt_vec = np.int0(np.linspace(-100,100,9))#x axis
        syn_locs = [0.3,0.9]
        prox_loc_options = [0,0,1]
        dist_loc_options = [0,1,1]
        num_of_syn_locs = len(syn_locs)#number of columns
        num_of_experiments = 3
        stim_protocol = '1_1'#'1_0','1_1'
        freq = 100
        syn_step = 1000//freq
        syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
        mean_max = 'max'
        fig_title = mean_max +' eff_cai + rho, stim_protocol: ' + stim_protocol +', num_of_stims: '+str(num_of_stims) +', freq: ' +str(freq) +'hz'
        fig = plt.figure()
        outer = gridspec.GridSpec(1, num_of_experiments,)
        cmap_rho = cm.get_cmap('bwr', 3)
        num_of_ticklabels = 2# number of jumps between xticklabels. The higher the number, the fewer tick labes
        if mean_max == 'max': vmax = 3
        elif mean_max == 'mean': vmax = 1.2
        eff_cai_matrices_1, eff_cai_matrices_2 = [],[]
        rho_matrices_1, rho_matrices_2 = [],[]
        
        plot_fig5_6_supp = 0#this in order to create new heatmaps of 3 passive syns at 3 locations (prox, center, dist) for each experiment
        if plot_fig5_6_supp:
            extra_eff_cai_matrices_1, extra_eff_cai_matrices_2, extra_eff_cai_matrices_3 = [],[],[]
            extra_rho_matrices_1, extra_rho_matrices_2, extra_rho_matrices_3 = [],[],[]
            term = 1#used this to get the right vectors after I add the 3 passive synapses from the line above
        else: term = 0
        
        for experiment_num in range(num_of_experiments):
            print('exp num: ',experiment_num)
            syn1_location = syn_locs[prox_loc_options[experiment_num]]
            syn2_location = syn_locs[dist_loc_options[experiment_num]]+0.0001
            
            eff_cai_mat_1, eff_cai_mat_2 = np.zeros((len(active_synapses_vec),len(dt_vec))), np.zeros((len(active_synapses_vec),len(dt_vec)))
            rho_mat_1, rho_mat_2 = np.zeros((len(active_synapses_vec),len(dt_vec))), np.zeros((len(active_synapses_vec),len(dt_vec)))
            if plot_fig5_6_supp:
                extra_eff_cai_mat_1, extra_eff_cai_mat_2, extra_eff_cai_mat_3 = np.zeros((len(active_synapses_vec),len(dt_vec))), np.zeros((len(active_synapses_vec),len(dt_vec))), np.zeros((len(active_synapses_vec),len(dt_vec)))
                extra_rho_mat_1, extra_rho_mat_2, extra_rho_mat_3 = np.zeros((len(active_synapses_vec),len(dt_vec))), np.zeros((len(active_synapses_vec),len(dt_vec))), np.zeros((len(active_synapses_vec),len(dt_vec)))
            for active_syn_num in range(len(active_synapses_vec)):
                row = len(active_synapses_vec) - active_syn_num - 1
                active_synapses = active_synapses_vec[active_syn_num]
                print('row: '+str(active_syn_num))
                for delta_t_num in range(len(dt_vec)):
                    col = delta_t_num
                    stims_1 = {weight+i:syn_times + dt_vec[delta_t_num] for i in range(active_synapses)}
                    stims_2 = {weight+i:syn_times for i in range(active_synapses)}
                    if active_synapses == 0:
                        stims_1 = {weight:[]}
                        stims_2 = {weight:[]}
                    if stim_protocol == '1_0': synapses = {0:{syn1_location:stims_1, syn2_location:{weight:[]}}}#for 1 0
                    elif stim_protocol == '1_1': synapses= {0:{syn1_location:stims_1, syn2_location:stims_2}}# for 1 1
                    if plot_fig5_6_supp:
                        term = -3
                        for loc in [0.2502,0.502,0.7502]:
                            synapses[0][loc] = {1:[]}
                    cell,_ = calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt, traces = False,
                                               active_loc = (0, syn1_location))
                    '''ONLY WORKS FOR 2 SYNAPSE LOCATIONS!!!! For more, need to work out mapping to relevant rho_GB vecs'''
                    if mean_max == 'mean':
                        eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.mean(utils.npa(cell.effcai_GB[0])), np.mean(utils.npa(cell.effcai_GB[-1]))
                    elif mean_max == 'max':
                        eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.max(utils.npa(cell.effcai_GB[0])), np.max(utils.npa(cell.effcai_GB[-1+term]))
#                         extra_eff_cai_mat_1[row, col], extra_eff_cai_mat_2[row, col], extra_eff_cai_mat_3[row, col] = np.max(utils.npa(cell.effcai_GB[-3])), np.max(utils.npa(cell.effcai_GB[-2])), np.max(utils.npa(cell.effcai_GB[-1]))
                    rho_mat_1[row, col], rho_mat_2[row, col] = np.heaviside(cell.rho_GB[0][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1+term][-1]-0.5,0.5)
#                     extra_rho_mat_1[row, col], extra_rho_mat_2[row, col], extra_rho_mat_3[row, col] = np.heaviside(cell.rho_GB[-3][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-2][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
            eff_cai_matrices_1.append(eff_cai_mat_1), eff_cai_matrices_2.append(eff_cai_mat_2)
            rho_matrices_1.append(rho_mat_1), rho_matrices_2.append(rho_mat_2)
#             extra_eff_cai_matrices_1.append(extra_eff_cai_mat_1), extra_eff_cai_matrices_2.append(extra_eff_cai_mat_2), extra_eff_cai_matrices_3.append(extra_eff_cai_mat_3)
#             extra_rho_matrices_1.append(extra_rho_mat_1), extra_rho_matrices_2.append(extra_rho_mat_2), extra_rho_matrices_3.append(extra_rho_mat_3),
    #                 subplot_num += 1
        
        min_val = np.min((np.min(eff_cai_matrices_1),np.min(eff_cai_matrices_2)))
        max_val = np.max((np.max(eff_cai_matrices_1),np.max(eff_cai_matrices_2)))
        vmin = np.round(min_val,3)
        vmax = np.round(max_val,3)
        
        gamma = 0.6

        cmap = colors.LinearSegmentedColormap.from_list('some name', [(0, 'white'),(0.25,'yellow'),(0.5, 'orange'),(0.75, 'red'),(1,'black')],)#gamma = 1000)#(zero_max*0.5,'yellow'),(zero_max,'orange')
        for experiment_num in range(num_of_experiments):
            syn1_location = syn_locs[prox_loc_options[experiment_num]]
            syn2_location = syn_locs[dist_loc_options[experiment_num]]
            inner = gridspec.GridSpecFromSubplotSpec(4, 1,#creates inner axes 4x1
                            subplot_spec=outer[experiment_num],)
            
            'synapse 1'
            ax1 = plt.Subplot(fig, inner[0])
            ax2 = plt.Subplot(fig, inner[1])
            'synapse 2'
            ax3 = plt.Subplot(fig, inner[2])
            ax4 = plt.Subplot(fig, inner[3])
            
            ax1.set_title('Synapse 1 - loc: {:.0f}'.format(syn1_location*dend_L))
            fig1 = ax1.imshow(eff_cai_matrices_1[experiment_num], cmap = cmap, vmin = vmin, vmax = vmax)#norm=colors.PowerNorm(gamma=gamma, ,))
            fig2 = ax2.imshow(rho_matrices_1[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
    
            ax3.set_title('Synapse 2 - loc: {:.0f}'.format(syn2_location*dend_L))
            fig3 = ax3.imshow(eff_cai_matrices_2[experiment_num], cmap = cmap, vmin = vmin, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, ))
            fig4 = ax4.imshow(rho_matrices_2[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
            ax4.set_xlabel('dt')
            
            if experiment_num == 0:
                ax1.set_yticks(np.arange(len(active_synapses_vec)))
                ax1.set_yticklabels(2*active_synapses_vec[::-1])
                ax1.set_ylabel('# active syns')
    
                ax2.set_yticks(np.arange(len(active_synapses_vec)))
                ax2.set_yticklabels(2*active_synapses_vec[::-1])
                ax2.set_ylabel('# active syns')
                
                ax3.set_yticks(np.arange(len(active_synapses_vec)))
                ax3.set_yticklabels(2*active_synapses_vec[::-1])
                ax3.set_ylabel('# active syns')
                
                ax4.set_yticks(np.arange(len(active_synapses_vec)))
                ax4.set_yticklabels(2*active_synapses_vec[::-1])
                ax4.set_ylabel('# active syns')
                
                fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB))
                
                cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
                cbar = fig.colorbar(fig3, cax=cbar_ax, )#cm.ScalarMappable(cmap=cmap, norm=colors.PowerNorm(gamma=gamma,))
                ctick1, ctick2, ctick3, ctick4 = vmin, np.round(vmin+(vmax-vmin)/3, 3), np.round(vmin+2*(vmax-vmin)/3, 3), vmax
                cbar.set_ticks([ctick1, ctick2, ctick3, ctick4]) 
                cbar.set_ticklabels([ctick1, ctick2, ctick3, ctick4])
            else:
                ax1.set_yticks(np.arange(len(active_synapses_vec)))
                ax1.set_yticklabels([])
                ax2.set_yticks(np.arange(len(active_synapses_vec)))
                ax2.set_yticklabels([])
                ax3.set_yticks(np.arange(len(active_synapses_vec)))
                ax3.set_yticklabels([])
                ax4.set_yticks(np.arange(len(active_synapses_vec)))
                ax4.set_yticklabels([])
            
            ax1.set_xticks(np.arange(len(dt_vec)))
            ax1.set_xticklabels([])
            
            ax2.set_xticks(np.arange(len(dt_vec)))
            ax2.set_xticklabels([])
            
            ax3.set_xticks(np.arange(len(dt_vec)))
            ax3.set_xticklabels([])
            
            ax4.set_xticks(np.arange(len(dt_vec)))
            ax4.set_xticklabels(dt_vec)
            
            for label in ax4.xaxis.get_ticklabels()[1::2]:
                    label.set_visible(False)
            
            fig.add_subplot(ax1)
            fig.add_subplot(ax2)
            fig.add_subplot(ax3)
            fig.add_subplot(ax4)
    
        if plot_fig5_6_supp:
            extra_eff_cai_matrices_list = [extra_eff_cai_matrices_1, extra_eff_cai_matrices_2, extra_eff_cai_matrices_3]
            extra_rho_matrices_list = [extra_rho_matrices_1, extra_rho_matrices_2, extra_rho_matrices_3]
            for experiment_num in range(num_of_experiments):
                syn1_location = syn_locs[prox_loc_options[experiment_num]]
                syn2_location = syn_locs[dist_loc_options[experiment_num]]
                fig, axes = plt.subplots(2,3, figsize = (15,10))
                fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB) +
                             '\nActive locations: ' + str(syn1_location*dend_L)+'um, ' + str(syn2_location*dend_L) + 'um')
                for col_num in range(len(axes[0])):
                    fig1 = axes[0,col_num].imshow(extra_eff_cai_matrices_list[col_num][experiment_num], cmap = cmap, norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
                    axes[1,col_num].imshow(extra_rho_matrices_list[col_num][experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
                    axes[0,col_num].set_title('Passive synapse loc: {:.0f}um'.format(syn_locs[col_num]*dend_L))
                    if col_num == 0:
                        axes[0,col_num].set_yticks(np.arange(len(active_synapses_vec)))
                        axes[0,col_num].set_yticklabels(active_synapses_vec[::-1])
                        axes[0,col_num].set_ylabel('# active syns')
            
                        axes[1,col_num].set_yticks(np.arange(len(active_synapses_vec)))
                        axes[1,col_num].set_yticklabels(active_synapses_vec[::-1])
                        axes[1,col_num].set_ylabel('# active syns')
                    else:
                        axes[0,col_num].set_yticklabels([])
                        axes[1,col_num].set_yticklabels([])
            
                    axes[0,col_num].set_xticks(np.arange(len(dt_vec)))
                    axes[0,col_num].set_xticklabels([])
                    axes[1,col_num].set_xticks(np.arange(len(dt_vec)))
                    axes[1,col_num].set_xticklabels(dt_vec)
                    axes[1,col_num].set_xlabel('dt')
                
                cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
                cbar = fig.colorbar(fig1, cax=cbar_ax, )#cm.ScalarMappable(cmap=cmap, norm=colors.PowerNorm(gamma=gamma,))
    
    elif protocol == 'eff_cai_heatmap_activeSyns_vs_activeSyns':
        num_of_sections = 1
        dend_L = 200 
        num_of_stims = 1
        active_synapses_vec = np.arange(0, 22, 3)#y axis
        syn_locs = [0.3,0.9]
        prox_loc_options = [0,0,1]
        dist_loc_options = [0,1,1]
        num_of_experiments = 3
        num_of_syn_locs = len(syn_locs)#number of columns
        stim_protocol = '1_1'#'1_0','1_1'
        freq = 200
        syn_step = 1000//freq
        syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
        mean_max = 'max'
        fig_title = mean_max +' eff_cai + rho , stim_protocol: ' + stim_protocol +', num_of_stims: '+str(num_of_stims) +', freq: ' +str(freq) +'hz'
        fig = plt.figure()
        outer = gridspec.GridSpec(1, num_of_experiments,)
        cmap_rho = cm.get_cmap('bwr', 3)
        num_of_ticklabels = 2# number of jumps between xticklabels. The higher the number, the fewer tick labes
        eff_cai_matrices_1, eff_cai_matrices_2 = [],[]
        rho_matrices_1, rho_matrices_2 = [],[]
        
        plot_fig5_6_supp = 0#this in order to create new heatmaps of 3 passive syns at 3 locations (prox, center, dist) for each experiment
        if plot_fig5_6_supp:
            extra_eff_cai_matrices_1, extra_eff_cai_matrices_2, extra_eff_cai_matrices_3 = [],[],[]
            extra_rho_matrices_1, extra_rho_matrices_2, extra_rho_matrices_3 = [],[],[]
            term = 1#used this to get the right vectors after I add the 3 passive synapses from the line above
        else: term = 0
        
        for experiment_num in range(num_of_experiments):
            print('exp num: ',experiment_num)
            syn1_location = syn_locs[prox_loc_options[experiment_num]]
            syn2_location = syn_locs[dist_loc_options[experiment_num]]+0.0001
            
            eff_cai_mat_1, eff_cai_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
            rho_mat_1, rho_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
            if plot_fig5_6_supp:
                extra_eff_cai_mat_1, extra_eff_cai_mat_2, extra_eff_cai_mat_3 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
                extra_rho_mat_1, extra_rho_mat_2, extra_rho_mat_3 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
            for active_syn_num1 in range(len(active_synapses_vec)):
#                 row = len(active_synapses_vec) - active_syn_num1 - 1
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
                    if stim_protocol == '1_0': synapses = {0:{syn1_location:stims_1, syn2_location:{weight:[]}}}#for 1 0
                    elif stim_protocol == '1_1': synapses= {0:{syn1_location:stims_1, syn2_location:stims_2}}# for 1 1
                    if plot_fig5_6_supp:
                        term = -3
                        for loc in [0.2502,0.502,0.7502]:
                            synapses[0][loc] = {1:[]}
                    cell,_ = calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt, traces = False,
                                               active_loc = (0, syn1_location))
                    '''ONLY WORKS FOR 2 SYNAPSE LOCATIONS!!!! For more, need to work out mapping to relevant rho_GB vecs'''
                    if mean_max == 'mean':
                        eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.mean(utils.npa(cell.effcai_GB[0])), np.mean(utils.npa(cell.effcai_GB[-1]))
                    elif mean_max == 'max':
                        eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.max(utils.npa(cell.effcai_GB[0])), np.max(utils.npa(cell.effcai_GB[-1+term]))
#                         extra_eff_cai_mat_1[row, col], extra_eff_cai_mat_2[row, col], extra_eff_cai_mat_3[row, col] = np.max(utils.npa(cell.effcai_GB[-3])), np.max(utils.npa(cell.effcai_GB[-2])), np.max(utils.npa(cell.effcai_GB[-1]))
                        '''add lines for extra matrices from above'''
                    rho_mat_1[row, col], rho_mat_2[row, col] = np.heaviside(cell.rho_GB[0][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1+term][-1]-0.5,0.5)
#                     extra_rho_mat_1[row, col], extra_rho_mat_2[row, col], extra_rho_mat_3[row, col] = np.heaviside(cell.rho_GB[-3][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-2][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
    #                 subplot_num += 1
            eff_cai_matrices_1.append(eff_cai_mat_1), eff_cai_matrices_2.append(eff_cai_mat_2)
            rho_matrices_1.append(rho_mat_1), rho_matrices_2.append(rho_mat_2)
#             extra_eff_cai_matrices_1.append(extra_eff_cai_mat_1), extra_eff_cai_matrices_2.append(extra_eff_cai_mat_2), extra_eff_cai_matrices_3.append(extra_eff_cai_mat_3)
#             extra_rho_matrices_1.append(extra_rho_mat_1), extra_rho_matrices_2.append(extra_rho_mat_2), extra_rho_matrices_3.append(extra_rho_mat_3),

        min_val = np.min((np.min(eff_cai_matrices_1),np.min(eff_cai_matrices_2)))
        max_val = np.max((np.max(eff_cai_matrices_1),np.max(eff_cai_matrices_2)))
        vmin = np.round(min_val,3)
        vmax = np.round(max_val,3)
        
        zero_max_1 = np.max([np.max(mat[:,0]) for mat in eff_cai_matrices_1])
        zero_max_2 = np.max([np.max(mat[-1,:]) for mat in eff_cai_matrices_2])
        zero_max = np.max((zero_max_1,zero_max_2))#max value for col/row with 0 active synapses
        
        non_zero_min_1 = np.min([np.min(mat[:,1:]) for mat in eff_cai_matrices_1])
        non_zero_min_2 = np.min([np.min(mat[:-1,:]) for mat in eff_cai_matrices_2])
        non_zero_min = np.min((non_zero_min_1,non_zero_min_2))
        gamma = 1
        tick1,tick2,tick3,tick4,tick5 = 0,vmax/4,2*vmax/4,3*vmax/4,vmax

        cmap = colors.LinearSegmentedColormap.from_list('some name', [(tick1/vmax, 'white'),(tick2/vmax,'yellow'),(tick3/vmax, 'orange'),
                                                                      (tick4/vmax, 'red'),(tick5/vmax,'black')],
                                                        gamma = gamma)
        cmap = 'hot_r'
        print(tick1,tick2,tick3,tick4,tick5,)
        for experiment_num in range(num_of_experiments):
            syn1_location = syn_locs[prox_loc_options[experiment_num]]
            syn2_location = syn_locs[dist_loc_options[experiment_num]]
            inner = gridspec.GridSpecFromSubplotSpec(4, 1,#creates inner axes 4x1
                            subplot_spec=outer[experiment_num],)
            
            'synapse 1'
            ax1 = plt.Subplot(fig, inner[0])
            ax2 = plt.Subplot(fig, inner[1])
            'synapse 2'
            ax3 = plt.Subplot(fig, inner[2])
            ax4 = plt.Subplot(fig, inner[3])
            ax1.set_title('Synapse 1 - loc: {:.0f}'.format(syn1_location*dend_L))
            fig1 = ax1.imshow(eff_cai_matrices_1[experiment_num], cmap = cmap,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
            fig2 = ax2.imshow(rho_matrices_1[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)

            ax3.set_title('Synapse 2 - loc: {:.0f}'.format(syn2_location*dend_L))
            fig3 = ax3.imshow(eff_cai_matrices_2[experiment_num], cmap = cmap,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
            fig4 = ax4.imshow(rho_matrices_2[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
            ax4.set_xlabel('# active syns (1)')
            
            if experiment_num == 0:
                ax1.set_yticks(np.arange(len(active_synapses_vec)))
                ax1.set_yticklabels(active_synapses_vec[::-1])
                ax1.set_ylabel('# active syns (2)')
    
                ax2.set_yticks(np.arange(len(active_synapses_vec)))
                ax2.set_yticklabels(active_synapses_vec[::-1])
                ax2.set_ylabel('# active syns (2)')
                
                ax3.set_yticks(np.arange(len(active_synapses_vec)))
                ax3.set_yticklabels(active_synapses_vec[::-1])
                ax3.set_ylabel('# active syns (2)')
                
                ax4.set_yticks(np.arange(len(active_synapses_vec)))
                ax4.set_yticklabels(active_synapses_vec[::-1])
                ax4.set_ylabel('# active syns (2)')
                
                fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB))
                
                cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
                cbar = fig.colorbar(fig1, cax=cbar_ax, )#cm.ScalarMappable(cmap=cmap, norm=colors.PowerNorm(gamma=gamma,))
                cbar.set_ticks([tick1,tick2,tick3,tick4,tick5,],3)#tick2,tick3,tick4,tick5,tick6]) 
                cbar.set_ticklabels(np.round([tick1,tick2,tick3,tick4,tick5,],3))
    
            else:
                ax1.set_yticks(np.arange(len(active_synapses_vec)))
                ax1.set_yticklabels([])
                ax2.set_yticks(np.arange(len(active_synapses_vec)))
                ax2.set_yticklabels([])
                ax3.set_yticks(np.arange(len(active_synapses_vec)))
                ax3.set_yticklabels([])
                ax4.set_yticks(np.arange(len(active_synapses_vec)))
                ax4.set_yticklabels([])
            ax1.set_xticks(np.arange(len(active_synapses_vec)))
            ax1.set_xticklabels([])
            
            ax2.set_xticks(np.arange(len(active_synapses_vec)))
            ax2.set_xticklabels([])
            
            ax3.set_xticks(np.arange(len(active_synapses_vec)))
            ax3.set_xticklabels([])
    
            ax4.set_xticks(np.arange(len(active_synapses_vec)))
            ax4.set_xticklabels(active_synapses_vec)

            fig.add_subplot(ax1)
            fig.add_subplot(ax2)
            fig.add_subplot(ax3)
            fig.add_subplot(ax4)

        if plot_fig5_6_supp:
            
            extra_eff_cai_matrices_list = [extra_eff_cai_matrices_1, extra_eff_cai_matrices_2, extra_eff_cai_matrices_3]
            extra_rho_matrices_list = [extra_rho_matrices_1, extra_rho_matrices_2, extra_rho_matrices_3]
            for experiment_num in range(num_of_experiments):
                syn1_location = syn_locs[prox_loc_options[experiment_num]]
                syn2_location = syn_locs[dist_loc_options[experiment_num]]
                fig, axes = plt.subplots(2,3, figsize = (15,10))
                fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB) +
                             '\nActive locations: ' + str(syn1_location*dend_L)+'um, ' + str(syn2_location*dend_L) + 'um')

                for col_num in range(len(axes[0])):
                    fig1 = axes[0,col_num].imshow(extra_eff_cai_matrices_list[col_num][experiment_num], cmap = cmap, norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
                    axes[1,col_num].imshow(extra_rho_matrices_list[col_num][experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
                    axes[0,col_num].set_title('Passive synapse loc: {:.0f}um'.format(syn_locs[col_num]*dend_L))
                    if col_num == 0:
                        axes[0,col_num].set_yticks(np.arange(len(active_synapses_vec)))
                        axes[0,col_num].set_yticklabels(active_synapses_vec[::-1])
                        axes[0,col_num].set_ylabel('# active syns (1)')
            
                        axes[1,col_num].set_yticks(np.arange(len(active_synapses_vec)))
                        axes[1,col_num].set_yticklabels(active_synapses_vec[::-1])
                        axes[1,col_num].set_ylabel('# active syns (1)')
                    else:
                        axes[0,col_num].set_yticklabels([])
                        axes[1,col_num].set_yticklabels([])
            
                    axes[0,col_num].set_xticks(np.arange(len(active_synapses_vec)))
                    axes[0,col_num].set_xticklabels([])
                    axes[1,col_num].set_xticks(np.arange(len(active_synapses_vec)))
                    axes[1,col_num].set_xticklabels(active_synapses_vec)
                    axes[1,col_num].set_xlabel('# active syns (2)')
                
                cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
                cbar = fig.colorbar(fig1, cax=cbar_ax, )#cm.ScalarMappable(cmap=cmap, norm=colors.PowerNorm(gamma=gamma,))
                cbar.set_ticks([tick1,tick6],3)#tick2,tick3,tick4,tick5,tick6]) 
                cbar.set_ticklabels(np.round([tick1,tick6],3))
    print('end')
