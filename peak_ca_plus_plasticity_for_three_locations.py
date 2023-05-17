'''
proximal and distal active synapse clusters, and central inactive synapses. checking peak calcium and plasticity at all three locations for different 
strength activations. made for new figure 5
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

sim_time = 1000
dt = 0.25
weight = 1
# step = 10
'{syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}'
start_1,start_2,start_3,start_4,start_5 = 10,4500,1500,2250,3000

def calcium_plasticity(num_of_sections, synapses, sim_time = 100, cell_model = 'bams', nmda_ratio = 2.38, excitatory = True, num_of_rec_locations = 30,
                       dend_L = 200, dt = 0.025, synapse_type = 'glu syn', traces = False, menorah = False, active_loc = (0, 0.5),
                                fig = None, num_of_traces = None):
        '''
        synapses: {syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}
        syn_locations: [syn1, syn2,...]
        syn_times, weights: [[syn1_1,syn1_2,...],[syn2_1,syn2_2,...]]
        '''
#         synapse_type = 'Behabadi'
#         synapse_type = 'glu syn'
#         print('Synapse type:',synapse_type)
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
#         cell.plot_TRs()
    #     pf.plot_ratios(cell)
#         print('Seg length = ',cell.dendrites[0].L)
    #     print(list(np.squeeze(synapses)))
#         print('active syns ', cell.active_synapses)
        for syn_sec_num in list(synapses):
    #         print(syn_sec_num)
    #         print(list(synapses[syn_sec_num]))
            for syn_loc in list(synapses[syn_sec_num]):
                for syn_weight in list(synapses[syn_sec_num][syn_loc]):
                    if excitatory:
#                         print([syn_loc])
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_sec_num][syn_loc][syn_weight]], dend_sec_num = syn_sec_num,
                                                        NMDA_ratio = nmda_ratio, synapse_type = synapse_type)
    #     rec_vecs = ['i','cai_CR','g_NMDA','g_AMPA','gmax_d_AMPA','ica_NMDA','ica_VDCC']
        rec_vecs = ['vsyn','ica_NMDA','ica_VDCC','effcai_GB','rho_GB']#,'m_VDCC', 'h_VDCC','gca_bar_VDCC','effcai_GB','gca_VDCC','gmax_AMPA'
#         rec_vecs = []
        cell.run_sim(rec_locations, sim_time, dt, recording_vecs = rec_vecs, v_init =  None)
    #     pf.plot_multiple_dendrite_heatmaps(cell, annotate = True, dt = dt, sim_time = sim_time)
        save_fig = False
        if traces:
            fig = pf.plot_conductances(cell,rec_vecs, save_fig, protocol, active_loc = active_loc, fig = fig, num_of_traces = num_of_traces)
#             print(cell.e_synapses[0].NMDA_ratio)
#             print(cell.dend_vs)
#             plt.figure()
#             plt.plot(cell.t, cell.dend_vs[2])
            
        if menorah:
            pf.plot_plasticity_tree(cell)    
        return cell, fig
    
def run_experiment(max_synapses, syn_locs = [0.3,0.65,1]):
    num_of_sections = 1
    dend_L = 200 
    num_of_stims = 1
    active_synapses_vec = np.arange(0, max_synapses+1, 5)#y axis 
    # active_synapses_vec = np.arange(3)
    # active_synapses_vec = np.arange(0, 10, 3)#y axis
#     syn_locs = [0.3,0.65,1]
    # prox_loc_options = [0,0,1]
    # dist_loc_options = [0,1,1]
    num_of_experiments = 1
    num_of_syn_locs = len(syn_locs)#number of columns
    # stim_protocol = '1_1'#'1_0','1_1'
    freq = 200
    syn_step = 1000//freq
    syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
    mean_max = 'max'
    fig_title = mean_max +' eff_cai + rho , num_of_stims: '+str(num_of_stims) +', freq: ' +str(freq) +'hz'
    fig, axes = plt.subplots(2,3, figsize = (14,8),)
    # outer = gridspec.GridSpec(2, 3,)
    cmap_rho = cm.get_cmap('bwr', 3)
    #         cmap_rho = cm.get_cmap('coolwarm', 3)
    num_of_ticklabels = 2# number of jumps between xticklabels. The higher the number, the fewer tick labes
    eff_cai_matrices_1, eff_cai_matrices_inact, eff_cai_matrices_2 = [],[],[]
    rho_matrices_1, rho_matrices_inact, rho_matrices_2 = [],[],[]
    #     if mean_max == 'max': vmax = 0.694
    #     elif mean_max == 'mean': vmax = 1.9
    
    # plot_fig5_6_supp = 0#this in order to create new heatmaps of 3 passive syns at 3 locations (prox, center, dist) for each experiment
    # if plot_fig5_6_supp:
    #     extra_eff_cai_matrices_1, extra_eff_cai_matrices_2, extra_eff_cai_matrices_3 = [],[],[]
    #     extra_rho_matrices_1, extra_rho_matrices_2, extra_rho_matrices_3 = [],[],[]
    #     term = 1#used this to get the right vectors after I add the 3 passive synapses from the line above
    # else: term = 0
    
    for experiment_num in range(num_of_experiments):
    #     print('exp num: ',experiment_num)
        syn1_location = syn_locs[0]
        inactive_syn_loc = syn_locs[1]
        syn2_location = syn_locs[2]
        
        eff_cai_mat_1, eff_cai_mat_inact, eff_cai_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
        rho_mat_1, rho_mat_inact, rho_mat_2 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
    #     if plot_fig5_6_supp:
    #         extra_eff_cai_mat_1, extra_eff_cai_mat_2, extra_eff_cai_mat_3 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
    #         extra_rho_mat_1, extra_rho_mat_2, extra_rho_mat_3 = np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec))), np.zeros((len(active_synapses_vec),len(active_synapses_vec)))
        for active_syn_num1 in range(len(active_synapses_vec)):
    #                 row = len(active_synapses_vec) - active_syn_num1 - 1
            col = active_syn_num1
            active_synapses1 = active_synapses_vec[active_syn_num1]
            print('row: '+str(active_syn_num1))
            for active_syn_num2 in range(len(active_synapses_vec)):
                
    #                     col = active_syn_num2
                row = len(active_synapses_vec) - active_syn_num2 - 1
    #                     print('r,c ', row, col)
                active_synapses2 = active_synapses_vec[active_syn_num2]
                stims_1 = {weight+i:syn_times for i in range(active_synapses1)}
                stims_2 = {weight+i:syn_times for i in range(active_synapses2)}
                if active_synapses1 == 0: stims_1 = {weight:[]}
                if active_synapses2 == 0: stims_2 = {weight:[]}
    #                 print(stims_1)
    #                 print(stims_2)
    #             if stim_protocol == '1_0': synapses = {0:{syn1_location:stims_1, syn2_location:{weight:[]}}}#for 1 0
    #             elif stim_protocol == '1_1': synapses= {0:{syn1_location:stims_1, syn2_location:stims_2}}# for 1 1
                synapses= {0:{syn1_location:stims_1, inactive_syn_loc:{1:[]}, syn2_location:stims_2}}
    #             if plot_fig5_6_supp:
    #                 term = -3
    #                 for loc in [0.2502,0.502,0.7502]:
    #                     synapses[0][loc] = {1:[]}
    # #                     print(synapses[0])
                print(synapses)
                cell,_ = calcium_plasticity(num_of_sections, synapses, dend_L = dend_L, sim_time = sim_time, dt = dt, traces = False,
                                           active_loc = (0, syn1_location))
                '''ONLY WORKS FOR 2 SYNAPSE LOCATIONS!!!! For more, need to work out mapping to relevant rho_GB vecs'''
    #             if mean_max == 'mean':
    #                 eff_cai_mat_1[row, col], eff_cai_mat_2[row, col] = np.mean(utils.npa(cell.effcai_GB[0])), np.mean(utils.npa(cell.effcai_GB[-1]))
    #             elif mean_max == 'max':
    #             print(len(cell.effcai_GB))
#                 print(len(synapses[0][syn1_location]),len(cell.effcai_GB))
                eff_cai_mat_1[row, col], eff_cai_mat_inact[row,col], eff_cai_mat_2[row, col] = np.max(utils.npa(cell.effcai_GB[0])), np.max(utils.npa(cell.effcai_GB[len(synapses[0][syn1_location])])), np.max(utils.npa(cell.effcai_GB[-1]))
    #                         extra_eff_cai_mat_1[row, col], extra_eff_cai_mat_2[row, col], extra_eff_cai_mat_3[row, col] = np.max(utils.npa(cell.effcai_GB[-3])), np.max(utils.npa(cell.effcai_GB[-2])), np.max(utils.npa(cell.effcai_GB[-1]))
                '''add lines for extra matrices from above'''
                rho_mat_1[row, col], rho_mat_inact[row, col], rho_mat_2[row, col] = np.heaviside(cell.rho_GB[0][-1]-0.5,0.5), np.heaviside(cell.rho_GB[len(synapses[0][syn1_location])][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
    #                     extra_rho_mat_1[row, col], extra_rho_mat_2[row, col], extra_rho_mat_3[row, col] = np.heaviside(cell.rho_GB[-3][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-2][-1]-0.5,0.5), np.heaviside(cell.rho_GB[-1][-1]-0.5,0.5)
    #                 subplot_num += 1
        eff_cai_matrices_1.append(eff_cai_mat_1), eff_cai_matrices_inact.append(eff_cai_mat_inact), eff_cai_matrices_2.append(eff_cai_mat_2)
        rho_matrices_1.append(rho_mat_1), rho_matrices_inact.append(rho_mat_inact), rho_matrices_2.append(rho_mat_2)
    #             extra_eff_cai_matrices_1.append(extra_eff_cai_mat_1), extra_eff_cai_matrices_2.append(extra_eff_cai_mat_2), extra_eff_cai_matrices_3.append(extra_eff_cai_mat_3)
    #             extra_rho_matrices_1.append(extra_rho_mat_1), extra_rho_matrices_2.append(extra_rho_mat_2), extra_rho_matrices_3.append(extra_rho_mat_3),
    #         plt.figure()
    #         for i in range(len(cell.rho_GB)):
    #             plt.plot(cell.rho_GB[i])
    # print(np.shape(eff_cai_matrices_inact))
    min_val = np.min((np.min(eff_cai_matrices_1), np.min(eff_cai_matrices_inact), np.min(eff_cai_matrices_2)))
    max_val = np.max((np.max(eff_cai_matrices_1), np.max(eff_cai_matrices_inact), np.max(eff_cai_matrices_2)))
    vmin = np.round(min_val,3)
    vmax = np.round(max_val,3)
    #         print(eff_cai_matrices_1, eff_cai_matrices_2)
    #         print(np.max(eff_cai_matrices_1))
    #         print(np.max(eff_cai_matrices_2))
    zero_max_1 = np.max([np.max(mat[:,0]) for mat in eff_cai_matrices_1])
    zero_max_inact = np.max([np.max(mat[:,0]) for mat in eff_cai_matrices_inact])
    zero_max_2 = np.max([np.max(mat[-1,:]) for mat in eff_cai_matrices_2])
    zero_max = np.max((zero_max_1,zero_max_inact,zero_max_2))#max value for col/row with 0 active synapses
    
    non_zero_min_1 = np.min([np.min(mat[:,1:]) for mat in eff_cai_matrices_1])
    non_zero_min_inact = np.min([np.min(mat[:,1:]) for mat in eff_cai_matrices_inact])
    #         print(non_zero_min_1, eff_cai_matrices_1)
    non_zero_min_2 = np.min([np.min(mat[:-1,:]) for mat in eff_cai_matrices_2])
    non_zero_min = np.min((non_zero_min_1,non_zero_min_inact,non_zero_min_2))
    #     print(zero_max,vmax)
    gamma = 1
    #         tick1,tick2,tick3,tick4,tick5 = 0,zero_max,zero_max+non_zero_min,zero_max + (vmax + non_zero_min)/2, vmax
    tick1,tick2,tick3,tick4,tick5 = 0,vmax/4,2*vmax/4,3*vmax/4,vmax
    #     print(vmax)
    #         tick1,tick2,tick3,tick4,tick5,tick6 = 0,(vmax-zero_max)/(2*vmax),(vmax-zero_max)/(1*vmax),2*(vmax-zero_max)/(1*vmax),vmax/2,vmax
    
    #     cmap = colors.LinearSegmentedColormap.from_list('some name', [(tick1, 'white'),(tick2*4,'yellow'),(tick3*4, 'orange'),(tick4, 'red'),(tick5,'black')],)#gamma = 1000)#(zero_max*0.5,'yellow'),(zero_max,'orange')
    cmap = colors.LinearSegmentedColormap.from_list('some name', [(tick1/vmax, 'white'),(tick2/vmax,'yellow'),(tick3/vmax, 'orange'),
                                                                  (tick4/vmax, 'red'),(tick5/vmax,'black')],
                                                    gamma = gamma)
    cmap = 'hot_r'
    print(tick1,tick2,tick3,tick4,tick5,)
    for experiment_num in range(num_of_experiments):
        syn1_location = syn_locs[0]
        inactive_syn_loc = syn_locs[1]
        syn2_location = syn_locs[2]
    #     inner = gridspec.GridSpecFromSubplotSpec(4, 1,#creates inner axes 4x1
    #                     subplot_spec=outer[experiment_num],)
        
        'synapse 1'
        ax1 = axes[0,0]
        ax2 = axes[1,0]
        'inactive synapse'
        ax3 = axes[0,1]
        ax4 = axes[1,1]
        'synapse 2'
        ax5 = axes[0,2]
        ax6 = axes[1,2]
    #             print('a')
    #         eff_cai_mat_1 = np.log((eff_cai_mat_1-np.min(eff_cai_mat_1))/(np.max(eff_cai_mat_1)-np.min(eff_cai_mat_1))+1)
        ax1.set_title('Proximal synapse')# - loc: {:.0f}'.format(syn1_location*dend_L))
    #         a = np.vstack((eff_cai_matrices_1[experiment_num][2,:],eff_cai_matrices_1[experiment_num][2,:],eff_cai_matrices_1[experiment_num][2,:]))
        fig1 = ax1.imshow(eff_cai_matrices_1[experiment_num], cmap = cmap, )#vmin = 0, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
        fig2 = ax2.imshow(rho_matrices_1[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
        ax2.set_xlabel('# active syns (1)')
        plt.colorbar(fig1, ax = ax1, fraction=0.046, pad=0.02, orientation = 'horizontal')
    #         ax1.imshow(a, cmap = cmap)
    #         eff_cai_mat_2 = np.log((eff_cai_mat_2-np.min(eff_cai_mat_2))/(np.max(eff_cai_mat_2)-np.min(eff_cai_mat_2))+1)
    #         b = np.vstack((eff_cai_matrices_2[experiment_num][:,0],eff_cai_matrices_2[experiment_num][:,0],eff_cai_matrices_2[experiment_num][:,0]))
    #         print(np.shape(a))
    #             print('b')
        ax3.set_title('Inactive synapse')# - loc: {:.0f}'.format(syn2_location*dend_L))
        fig3 = ax3.imshow(eff_cai_matrices_inact[experiment_num], cmap = cmap, )#vmin = 0, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
        fig4 = ax4.imshow(rho_matrices_inact[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
        ax4.set_xlabel('# active syns (1)')
        plt.colorbar(fig3, ax = ax3, fraction=0.046, pad=0.02, orientation = 'horizontal')
    #         ax2.imshow(b, cmap = cmap)
        ax5.set_title('Distal synapse')
        fig5 = ax5.imshow(eff_cai_matrices_2[experiment_num], cmap = cmap, )#vmin = 0, vmax = vmax,)# norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
        fig6 = ax6.imshow(rho_matrices_2[experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
        ax6.set_xlabel('# active syns (1)')
        plt.colorbar(fig5, ax = ax5, fraction=0.046, pad=0.02, orientation = 'horizontal')
    #     if experiment_num == 0:
    #                 print('c')
        ax1.set_yticks(np.arange(len(active_synapses_vec)))
        ax1.set_yticklabels(active_synapses_vec[::-1])
        ax1.set_ylabel('eff_ca\n# active syns (2)')
    
        ax2.set_yticks(np.arange(len(active_synapses_vec)))
        ax2.set_yticklabels(active_synapses_vec[::-1])
        ax2.set_ylabel('rho\n# active syns (2)')
        
    #     ax3.set_yticks(np.arange(len(active_synapses_vec)))
    #     ax3.set_yticklabels(active_synapses_vec[::-1])
    #     ax3.set_ylabel('# active syns (2)')
    #     
    #     ax4.set_yticks(np.arange(len(active_synapses_vec)))
    #     ax4.set_yticklabels(active_synapses_vec[::-1])
    #     ax4.set_ylabel('# active syns (2)')
        
        fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB) + '\n locations = '+str(syn_locs))
        
    #     cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
    #     cbar = fig.colorbar(fig1, cax=cbar_ax, )#cm.ScalarMappable(cmap=cmap, norm=colors.PowerNorm(gamma=gamma,))
    #             cbar = fig.colorbar(fig2, cax=cbar_ax, )
    #     cbar.set_ticks([tick1,tick2,tick3,tick4,tick5,],3)#tick2,tick3,tick4,tick5,tick6]) 
    #     cbar.set_ticklabels(np.round([tick1,tick2,tick3,tick4,tick5,],3))
    #             print(np.round([tick1,tick2,tick3,tick4,tick5,tick6],3))
        #     cbar.set_ticklabels(['depression','no change','potentiation'])
    #                 cbar_ax_rho = fig.add_axes([0.72, 0.15, 0.01, 0.7])
    #                 cbar_rho = fig.colorbar(fig2, cax=cbar_ax_rho)
    #                 cbar_rho.set_ticks([1/6, 3/6, 5/6])  
    #                 cbar_rho.set_ticklabels([])
    
    #     else:
        ax3.set_yticks(np.arange(len(active_synapses_vec)))
        ax3.set_yticklabels([])
        ax4.set_yticks(np.arange(len(active_synapses_vec)))
        ax4.set_yticklabels([])
        ax5.set_yticks(np.arange(len(active_synapses_vec)))
        ax5.set_yticklabels([])
        ax6.set_yticks(np.arange(len(active_synapses_vec)))
        ax6.set_yticklabels([])
    #             print('d')
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
        kws = dict(ha='center', va='center', color='deepskyblue', fontsize=15)#weight='semibold'
        for axi in axes:
            for ax in axi:
                for i in range(len(active_synapses_vec)):
                    for j in range(len(active_synapses_vec)):
                            ax.text(j, i, active_synapses_vec[-1-i]+active_synapses_vec[j], **kws)
    #         for label in ax2.xaxis.get_ticklabels()[1::2]:
    #                 label.set_visible(False)
    #             print('e')
    #     fig.add_subplot(ax1)
    #     fig.add_subplot(ax2)
    #     fig.add_subplot(ax3)
    #     fig.add_subplot(ax4)
    #         fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Calcium in Dendrites - Paper/Figures/Figure 5/'+str(time.time())[-3:]+', max ca =  ' +str(vmax) + ', gamma = ' + str(gamma) +'.png')
    #         print('f')
    # if plot_fig5_6_supp:
    #     
    #     extra_eff_cai_matrices_list = [extra_eff_cai_matrices_1, extra_eff_cai_matrices_2, extra_eff_cai_matrices_3]
    #     extra_rho_matrices_list = [extra_rho_matrices_1, extra_rho_matrices_2, extra_rho_matrices_3]
    #     for experiment_num in range(num_of_experiments):
    #         syn1_location = syn_locs[prox_loc_options[experiment_num]]
    #         syn2_location = syn_locs[dist_loc_options[experiment_num]]
    #         fig, axes = plt.subplots(2,3, figsize = (15,10))
    #         fig.suptitle(fig_title + ', theta_d = ' + str(cell.e_synapses[0].theta_d_GB) + ', theta_p = ' + str(cell.e_synapses[0].theta_p_GB) +
    #                      '\nActive locations: ' + str(syn1_location*dend_L)+'um, ' + str(syn2_location*dend_L) + 'um')
    # #             for row_num in range(len(axes)):
    # #             print(np.shape(extra_eff_cai_matrices_list))
    # #             print(np.shape(extra_eff_cai_matrices_list[0]))
    #         for col_num in range(len(axes[0])):
    #             fig1 = axes[0,col_num].imshow(extra_eff_cai_matrices_list[col_num][experiment_num], cmap = cmap, norm=colors.PowerNorm(gamma=gamma, vmin = 0, vmax = vmax,))
    #             axes[1,col_num].imshow(extra_rho_matrices_list[col_num][experiment_num], cmap = cmap_rho, vmin = 0, vmax = 1)
    #             axes[0,col_num].set_title('Passive synapse loc: {:.0f}um'.format(syn_locs[col_num]*dend_L))
    #             if col_num == 0:
    #                 axes[0,col_num].set_yticks(np.arange(len(active_synapses_vec)))
    #                 axes[0,col_num].set_yticklabels(active_synapses_vec[::-1])
    #                 axes[0,col_num].set_ylabel('# active syns (1)')
    #     
    #                 axes[1,col_num].set_yticks(np.arange(len(active_synapses_vec)))
    #                 axes[1,col_num].set_yticklabels(active_synapses_vec[::-1])
    #                 axes[1,col_num].set_ylabel('# active syns (1)')
    #             else:
    #                 axes[0,col_num].set_yticklabels([])
    #                 axes[1,col_num].set_yticklabels([])
    #     
    #             axes[0,col_num].set_xticks(np.arange(len(active_synapses_vec)))
    #             axes[0,col_num].set_xticklabels([])
    #             axes[1,col_num].set_xticks(np.arange(len(active_synapses_vec)))
    #             axes[1,col_num].set_xticklabels(active_synapses_vec)
    #             axes[1,col_num].set_xlabel('# active syns (2)')
    #         
    #         cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
    #         cbar = fig.colorbar(fig1, cax=cbar_ax, )#cm.ScalarMappable(cmap=cmap, norm=colors.PowerNorm(gamma=gamma,))
    #         cbar.set_ticks([tick1,tick6],3)#tick2,tick3,tick4,tick5,tick6]) 
    #         cbar.set_ticklabels(np.round([tick1,tick6],3))
    # #                 fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Calcium in Dendrites - Paper/Figures/Figure 5/fig '+str(experiment_num)+'.png')
    
#     fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Paper/Figures/Figure 5 new/calcium_plasticity_40_deepbluesky.png') 
    winsound.Beep(1000,200)
    print('end')









    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    