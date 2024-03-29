import matplotlib.pyplot as plt
import numpy as np
import utils
plt.ion()
import string
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from matplotlib.cm import ScalarMappable
import matplotlib.cm 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
from pylab import *

def npa(vec):
    """
    takes a neuron vector and returns a numpy array
    """
    return vec.as_numpy()


def plot_voltages(time_vec, syn_currents, dend_voltages, soma_voltages, run_type, weights = None, syn_loc = None, dend_L = 500):
    '''
    plot voltages at different synapses and the soma, and plots the synapse currents as well
    parameters - 
    
    '''
    if run_type == 'single':
        fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8,5))
    elif run_type == 'multiple':
        fig, ax = plt.subplots(3, 1, sharex=True, figsize=(9,9))
        
    ax[1].set_title('Currents')
    
    if run_type == 'single':
        ax[0].set_title('Voltages')
        ax[0].plot(time_vec, soma_voltages,color='k',label='soma(0.5)')
        for i,v in list(enumerate(dend_voltages)):
            ax[0].plot(time_vec, v, label = 'dend({:.2})'.format(rec_locations[i]))
        for i,v in list(enumerate(syn_currents)):
            ax[1].plot(time_vec, syn_currents[v][weights[i]], label = 'dend({:.2})'.format(rec_locations[i]))
    
    elif run_type == 'multiple':
        fig.suptitle('Synapse location = {:.0f}um'.format(syn_loc*dend_L), fontsize=16)
        ax[0].set_title('Dendrite voltages')
        for i,v in list(enumerate(dend_voltages[syn_loc])):
            ax[0].plot(time_vec, dend_voltages[syn_loc][v],
                     label = 'weight({:.2f})'.format(v))
        for i,v in list(enumerate(syn_currents[syn_loc])):
            ax[1].plot(time_vec, syn_currents[syn_loc][v], 
                       label = 'weight({:.1f})'.format(v))
        ax[2].set_title('Soma voltages')
        for i,v in list(enumerate(soma_voltages)):
            ax[2].plot(time_vec, soma_voltages[v], 
                     label = 'weight({:.2f})'.format(v))
        
    for i in range(len(ax)):
        
        ax[i].spines["top"].set_visible(False)
        ax[i].spines["right"].set_visible(False)
    ax[0].legend()
    plt.show()


def plot_voltages_multiple_distances(cell, time_vec, syn_currents, dend_voltages, soma_voltages, syn_loc, save_fig = False,
                                     fig_name = None):
    
    fig, ax = plt.subplots(len(syn_loc), 3, sharex=True, sharey = 'col', figsize=(9,9))
    ax[0,0].set_title('Dendrite voltages')
    ax[0,2].set_title('Soma voltages')
    
    for row in range(len(syn_loc)):
        syn_dist = syn_loc[row][0]*cell.dend.L
        ax[row,1].set_title('Synapse location = {:.0f}um\n Currents'.format(syn_dist),)
        for i,v in list(enumerate(dend_voltages[syn_loc[row][0]])):
            voltage = dend_voltages[syn_loc[row][0]][v]
            right_index = utils.get_index_at_half_value(voltage)
            ax[row,0].plot(voltage,
                     label = 'weight({:.2f})'.format(v))
            ax[row,0].scatter(right_index,voltage[right_index])
        for i,v in list(enumerate(syn_currents[syn_loc[row][0]])):
            voltage = syn_currents[syn_loc[row][0]][v]
            right_index = utils.get_index_at_half_value(voltage, eps = 0.000005)
            ax[row,1].plot(voltage, 
                       label = 'weight({:.1f})'.format(v))
            ax[row,1].scatter(right_index,voltage[right_index])
        for i,v in list(enumerate(soma_voltages[syn_loc[row][0]])):
            voltage = soma_voltages[syn_loc[row][0]][v]
            right_index = utils.get_index_at_half_value(voltage)
            ax[row,2].plot(voltage, 
                     label = 'weight({:.2f})'.format(v))
            ax[row,2].scatter(right_index,voltage[right_index])

    ax[0,0].legend()
    plt.tight_layout()
    if save_fig:
        fig_folder = []
        plt.savefig(fig_folder+fig_name+', traces.png')
    plt.show()
    
    
def plot_proximal_vs_distal(dend_voltages, soma_v, time_vec, dend_L = 200):
    locations = list(dend_voltages)
    weights = list(dend_voltages[locations[0]])
    fig1, axes1 = plt.subplots(1,3, sharey = True, figsize = (8,3), )
    fig2, axes2 = plt.subplots(1,3, sharey = True, figsize = (8,3,))
    fig1.suptitle('Dendrite voltages')
    fig2.suptitle('Soma voltages')
    for loc_num in range(len(locations)):
        loc = locations[loc_num]
        ax1 = axes1[loc_num]
        ax2 = axes2[loc_num]
        ax1.set_title('Synapse location = '+str(round(loc*dend_L))+'um')
        ax2.set_title('Synapse location = '+str(round(loc*dend_L))+'um')
        for weight_num in range(len(weights)):
            weight = weights[weight_num]
            ax1.plot(time_vec, dend_voltages[loc][weight], label = round(weight))
            ax2.plot(time_vec, soma_v[loc][weight], label = round(weight))
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
    ax1.legend(title = 'Weights')
    ax2.legend(title = 'Weights')
    fig1.tight_layout()
    fig2.tight_layout()
    plt.show()
   
        
def plot_peak_responses(cell, dend_traces, soma_traces, weights, dend_trace_type, peak_type = 'mean', save_fig = False, 
                        fig_name = None,): 
    '''
    parameters - 
    dend_traces: dict of the form {syn_location: {weight:voltage_vec}}, 
    soma_traces: 'dendrites':{location: {weight:voltage_vec}}}
     
    weights: array of weights
    '''
    
    peaks_soma = utils.get_trace_avgs(soma_traces, peak_type, 'voltage')
    peaks_dendrite = utils.get_trace_avgs(dend_traces, peak_type, dend_trace_type)
    syn_locations = list(peaks_soma)
    rest = soma_traces[syn_locations[0]][weights[0]][0]
    TRInds = [cell.loc_2_seg(loc) for loc in syn_locations]
    TRs = cell.soma_TRs
    v_ratios = cell.soma_ratios
    if dend_trace_type == 'current':
        scaling_factors = {cell.seg_2_loc(ind):TRs[ind] for ind in TRInds}
    else:
        scaling_factors = {cell.seg_2_loc(ind):v_ratios[ind] for ind in TRInds}    
    fig, axes = plt.subplots(1,3, sharey = False, figsize = (12,3))
    ax0 = axes[0]
    ax1 = axes[1]
    ax2 = axes[2]

    cmap = plt.get_cmap('jet')
    colors = cmap(np.linspace(0, 1, len(syn_locations)))
    for i in range(len(syn_locations)):
        syn_dist = cell.loc_2_dist(cell.dend(syn_locations[i]))
        loc = syn_locations[i]        
        color = colors[i]
        ax0.plot(weights,peaks_soma[loc], label = str(round(syn_dist)) + 'um', c = color)
        ax0.set_xlabel('Synaptic conductance (nS)')
        if dend_trace_type == 'current':
            if i == 0:
                ax0twin = ax0.twinx()
            ax0.set_ylabel('Soma voltage (mV)')
            ax0twin.plot(weights,peaks_dendrite[loc],  linestyle = '--', 
                    label = 'dendrite ' + str(round(syn_dist)) + 'um', c = color )
            ax0twin.set_ylabel('Dendritic current (nA)')
        else:
            ax0.set_ylabel('Voltage (mV)')
            ax0.plot(weights,peaks_dendrite[loc], linestyle = '--', c = color)            
 
        locs = np.linspace(cell.loc_2_dist(cell.dend(0)),cell.loc_2_dist(cell.dend(1)),len(cell.soma_TRs))
        if dend_trace_type == 'current':
            ax1.plot(locs, cell.soma_TRs)
            ax1.set_ylabel('Transfer resistance (Mega Ohms)')
            ax1.set_title('Transfer resistance from soma')
        else: 
            ax1.plot(locs, cell.soma_ratios)
            ax1.set_title('Transfer voltage ratio (Mega Ohms)')
            ax1.set_ylabel('Ratio')
        ax1.set_xlabel('Distance from soma (microns)')

        ax2.plot(weights,np.array(peaks_soma[loc])-rest,  label = 'soma ' + str(int(syn_locations[i]*cell.dend.L)) + 'um', linewidth = 4, c = color, alpha = 0.5)
        if dend_trace_type == 'current':
            ax2.plot(weights,(np.array(peaks_dendrite[loc]))*scaling_factors[loc], linestyle = '--', 
                    label = 'dendrite ' + str(round(syn_dist)) + 'um', alpha = 1, c = color)
        else:
            ax2.plot(weights,((np.array(peaks_dendrite[loc]))-rest)*scaling_factors[loc],  linestyle = '--', 
            label = 'dendrite ' + str(round(syn_dist)) + 'um', alpha = 1, c = color)
        ax2.set_xlabel('Synaptic conductance (nS)')
        ax2.set_ylabel('Scaled voltage (mV)')
         
    ax0.legend()
    plt.tight_layout()
    if save_fig:
        fig_folder = []
        plt.savefig(fig_folder+fig_name+', '+dend_trace_type+'.png')
    plt.show()
    return scaling_factors


def plot_fig_3(voltages, syn_currents, weights, mean_max = 'mean', dend_L = 200, num_of_synapse_locations = 3, 
               save_fig = True, NMDA_ratio = 2.38):
    '''
    voltages should be = {'dendrite:{(l1,l2):{(w1,w2):vec},,,}, {(l2,l1):{(w2,w1):vec,,,}},'soma':{...}}
    syn_currents should be = {(l1,l2):{(w1,w2):vec,(w1,w2):vec...},(l1,l2):{...}}
    ''' 
    print('starting to plot')
    soma_locations = list(voltages['soma'])
    dend_locations = list(voltages['dendrite'])
    weights = list(voltages['soma'][soma_locations[0]])
    num_of_weights = int(np.sqrt(len(weights)))
    num_of_soma_matrices = int(num_of_synapse_locations*(num_of_synapse_locations+1)*0.5)#num of plots needed for soma. Should be n(n+1)/2
    num_of_dend_matrices = int(num_of_synapse_locations*num_of_synapse_locations)

    soma_peak_matrices = [np.zeros((num_of_weights,num_of_weights)) for i in range(num_of_soma_matrices)]
    dend_peak_matrices = [np.zeros((num_of_weights,num_of_weights)) for i in range(num_of_dend_matrices)]
    syn_current_peak_matrices = [np.zeros((num_of_weights,num_of_weights)) for i in range(num_of_dend_matrices)]
    duration_matrices = [np.zeros((num_of_weights,num_of_weights)) for i in range(num_of_dend_matrices)]
    
    for num in range(num_of_soma_matrices):
        mat = soma_peak_matrices[num]
        l1,l2 = soma_locations[num]

        weight_index = 0
        for num1 in range(num_of_weights):
            for num2 in range(num_of_weights):
                v = voltages['soma'][l1,l2][weights[weight_index]]
                right_index = utils.get_index_at_half_value(v)

                mat[num1,num2] = np.mean(v[:right_index])
                weight_index += 1

    for num in range(num_of_dend_matrices):
        mat = dend_peak_matrices[num]
        l1,l2 = dend_locations[num]
        weight_index = 0
        for num1 in range(num_of_weights):
            for num2 in range(num_of_weights):
                v = voltages['dendrite'][l1,l2][weights[weight_index]]
                right_index = utils.get_index_at_half_value(v)
                mat[num1,num2] = np.mean(v[:right_index])
                weight_index += 1
     
    for num in range(num_of_dend_matrices):
        mat = syn_current_peak_matrices[num]
        l1,l2 = dend_locations[num]
        weight_index = 0
        for num1 in range(num_of_weights):
            for num2 in range(num_of_weights):
                v = -1*syn_currents[l1,l2][weights[weight_index]]
                right_index = utils.get_index_at_half_value(v)

                mat[num1,num2] = np.mean(v[:np.squeeze(right_index)])
                weight_index += 1
    
    for num in range(num_of_dend_matrices):
        mat = duration_matrices[num]
        l1,l2 = dend_locations[num]
        weight_index = 0
        for num1 in range(num_of_weights):
            for num2 in range(num_of_weights):
                v = -1*syn_currents[l1,l2][weights[weight_index]]
                right_index = utils.get_index_at_half_value(v)

                mat[num1,num2] = right_index
                weight_index += 1
                
    if num_of_synapse_locations == 5:
        dend_layout = """
        uvwxy
        pqrst
        klmno
        fghij
        abcde
        """
        soma_layout = """
        XXXXo
        XXXmn
        XXjkl
        Xfghi
        abcde
        """
    elif num_of_synapse_locations == 3:
         dend_layout = """
         ghi
         def
         abc
         """
         soma_layout = """
         XXf
         Xde
         abc
         """
    elif num_of_synapse_locations == 1:
         dend_layout = """
         a
         """
         soma_layout = """
         a
         """
    fig, ax = plt.subplots()
    loc = list(syn_currents.keys())[-1]
    n = len(weights)//4
    eps = 0.0000001
    
    ax.plot(syn_currents[loc][weights[0]], label = 'weights ='+str(weights[0][0])+','+str(weights[0][1]))
    right_index = utils.get_index_at_half_value(syn_currents[loc][weights[0]], eps = eps)
    ax.scatter(right_index, syn_currents[loc][weights[0]][right_index])
    
    ax.plot(syn_currents[loc][weights[n]], label = 'weights ='+str(weights[n][0])+','+str(weights[n][1]))
    right_index = utils.get_index_at_half_value(syn_currents[loc][weights[n]], eps = eps)
    ax.scatter(right_index,syn_currents[loc][weights[n]][right_index])
    
    ax.plot(syn_currents[loc][weights[2*n]], label = 'weights ='+str(weights[2*n][0])+','+str(weights[2*n][1]))
    right_index = utils.get_index_at_half_value(syn_currents[loc][weights[2*n]], eps = eps)
    ax.scatter(right_index,syn_currents[loc][weights[2*n]][right_index])
    
    ax.plot(syn_currents[loc][weights[-1]], label = 'weights ='+str(weights[-1][0])+','+str(weights[-1][1]))
    right_index = utils.get_index_at_half_value(syn_currents[loc][weights[-1]], eps = eps)
    ax.scatter(right_index,syn_currents[loc][weights[-1]][right_index])
    
    ax.plot(syn_currents[loc][weights[-15]], label = 'weights ='+str(weights[-15][0])+','+str(weights[-15][1]))
    right_index = utils.get_index_at_half_value(syn_currents[loc][weights[-15]], eps = eps)
    ax.scatter(right_index,syn_currents[loc][weights[-15]][right_index])    
    
    ax.set_title('Driver loc = {}, mod/clamp location = {}'.format(loc[0]*dend_L,loc[1]*dend_L))
    plt.legend()
    
    num_of_tick_labels = 4
    x_tick_labels = np.repeat(None,num_of_weights)
    x_tick_labels[0],x_tick_labels[-1] = weights[0][0],weights[-1][0]
    y_tick_labels = np.repeat(None,num_of_weights)
    y_tick_labels[0],y_tick_labels[-1] = weights[0][1],weights[-1][1]
    
    fig = plt.figure(constrained_layout=True, figsize = (10,10))
    axes = fig.subplot_mosaic(dend_layout)
    fig.suptitle(f'Synapse currents, NMDA_ratio = {NMDA_ratio}', size = 20)
    vmin, vmax = np.min(syn_current_peak_matrices),np.max(syn_current_peak_matrices)

    plot = True
    if plot:
        for i in range(num_of_dend_matrices):
            l1,l2 = list(syn_currents.keys())[i]
            ax = axes[string.ascii_lowercase[i]]
            s = ax.imshow(np.rot90(syn_current_peak_matrices[i]),cmap=plt.get_cmap('jet'), vmin = vmin, vmax = vmax)
            ax.set_xticks(np.arange(num_of_weights))
            ax.set_xticklabels(x_tick_labels)
            ax.set_yticks(np.arange(num_of_weights))
            ax.set_yticklabels(np.flip(y_tick_labels))
            ax.set_xlabel('driver_w', size = 10)
            ax.set_ylabel('mod_value', size = 10)
            ax.set_title('Driver = {:.0f}um, modulator =  {:.0f}um'.format(l1*dend_L, l2*dend_L), size = 12)
        if save_fig:
            fig_folder = ''
            plt.savefig(fig_folder+'nmda_ratio = '+str(NMDA_ratio)+'.png')
         
        fig = plt.figure(constrained_layout=True)
        axes = fig.subplot_mosaic(dend_layout)
        fig.suptitle('Dendrite voltages', size = 20)
        vmin, vmax = np.min(dend_peak_matrices),np.max(dend_peak_matrices)
        for i in range(num_of_dend_matrices):
            l1,l2 = dend_locations[i]
            ax = axes[string.ascii_lowercase[i]]
            s = ax.imshow(np.rot90(dend_peak_matrices[i]),cmap=plt.get_cmap('jet'), vmin = vmin, vmax = vmax)
            ax.set_xticks(np.arange(num_of_weights))
            ax.set_xticklabels(x_tick_labels)
            ax.set_yticks(np.arange(num_of_weights))
            ax.set_yticklabels(np.flip(y_tick_labels))
            ax.set_xlabel('driver_w', size = 10)
            ax.set_ylabel('mod_value', size = 10)
            ax.set_title('Driver = {:.0f}um, modulator =  {:.0f}um'.format(l1*dend_L, l2*dend_L), size = 12)

        if save_fig:
            fig_folder = ''
            plt.savefig(fig_folder+'dend, fig_3.png')
    
        fig = plt.figure(constrained_layout=True)
        axes = fig.subplot_mosaic(soma_layout, empty_sentinel = "X")
        fig.suptitle('Soma voltages', size = 20)
        vmin, vmax = np.min(soma_peak_matrices),np.max(soma_peak_matrices)
        for i in range(num_of_soma_matrices):
            l1,l2 = soma_locations[i]
            ax = axes[string.ascii_lowercase[i]]
            s = ax.imshow(np.rot90(soma_peak_matrices[i]),cmap=plt.get_cmap('jet'), vmin = vmin, vmax = vmax)
            ax.set_xticks(np.arange(num_of_weights))
            ax.set_xticklabels(x_tick_labels)
            ax.set_yticks(np.arange(num_of_weights))
            ax.set_yticklabels(np.flip(y_tick_labels))
            ax.set_xlabel('driver_w', size = 10)
            ax.set_ylabel('mod_value', size = 10)
            ax.set_title('Driver = {:.0f}um, modulator =  {:.0f}um'.format(l1*dend_L, l2*dend_L), size = 12)
        if save_fig:
            fig_folder = ''
            plt.savefig(fig_folder+'soma, fig_3.png')
        
        fig = plt.figure(constrained_layout=True)
        axes = fig.subplot_mosaic(dend_layout)
        fig.suptitle('Durations', size = 20)
        vmin, vmax = np.min(duration_matrices),np.max(duration_matrices)
        for i in range(num_of_dend_matrices):
            l1,l2 = dend_locations[i]
            ax = axes[string.ascii_lowercase[i]]
            s = ax.imshow(np.rot90(duration_matrices[i]),cmap=plt.get_cmap('jet'), vmin = vmin, vmax = vmax)
            ax.set_xticks(np.arange(num_of_weights))
            ax.set_xticklabels(x_tick_labels)
            ax.set_yticks(np.arange(num_of_weights))
            ax.set_yticklabels(np.flip(y_tick_labels))
            ax.set_xlabel('driver_w', size = 10)
            ax.set_ylabel('mod_value', size = 10)
            ax.set_title('Driver = {:.0f}um, modulator =  {:.0f}um'.format(l1*dend_L, l2*dend_L), size = 12)
            
        fig = plt.figure(constrained_layout=True)
        axes = fig.subplot_mosaic(dend_layout)
        fig.suptitle('Voltage duration scatter', size = 20)
        for i in range(num_of_dend_matrices):
            ax = axes[string.ascii_lowercase[i]]
            rows, cols = np.shape(duration_matrices[i])

            ax.scatter(np.array(duration_matrices[i]).ravel(),np.array(syn_current_peak_matrices[i]).ravel(), color = 'k', s = 4)
            ax.set_xlabel('Duration', size = 10)
            ax.set_ylabel('Mean syn voltage', size = 10)
            ax.set_title('Driver = {:.0f}um, modulator =  {:.0f}um'.format(l1*dend_L, l2*dend_L), size = 12)
            ax.get_shared_x_axes().join(ax, axes['a'])
            ax.get_shared_y_axes().join(ax, axes['a'])
    return syn_current_peak_matrices


def plot_currents(syn_currents, num_of_locations = 5, dend_L = 200, v_clmp = False):
    '''
    plot of nXn subplots showing current at driver synapse (with constant weight) while modulator changes weights. Each subplot is a different 
    loc tuple
    syn_c format: {(driver_loc1, modulator_loc1):{(driver_weight,modulator_weight1:vec,(driver_weight,modulator_weight2}:vec,...,
    (driver_loc1, modulator_loc2):{(driver_weight,modulator_weight1):vec}}
    '''
    if num_of_locations == 3:
         plot_layout = """
         ghi
         def
         abc
         """
    elif num_of_locations == 5:
        plot_layout = """
        uvwxy
        pqrst
        klmno
        fghij
        abcde
        """
    if num_of_locations == 6:
        plot_layout = """
        456789
        yz0123
        stuvwx
        mnopqr
        ghijkl
        abcdef
        """
    locations_list = list(syn_currents.keys())
    weights_list = list(syn_currents[locations_list[0]].keys())
    num_of_weights = len(syn_currents[locations_list[0]].keys())
    
    fig_main = plt.figure(constrained_layout=True, figsize = (10,10),)
    axes_main = fig_main.subplot_mosaic(plot_layout, )#subplot_kw = dict(sharex = axes[0])
    fig_main.suptitle('Synapse currents, Driver weight = {}'.format(weights_list[0][0]), size = 20)
    fig_1 = plt.figure(constrained_layout=True, figsize = (10,10),)
    fig_2 = plt.figure(constrained_layout=True, figsize = (10,10),)
    axes_1 = fig_1.subplot_mosaic(plot_layout, )
    axes_2 = fig_2.subplot_mosaic(plot_layout, )
    fig_1.suptitle('Synapse current 1st derivative', size = 20)
    fig_2.suptitle('Synapse current 2nd derivative', size = 20)
    fig_3 = plt.figure(constrained_layout=True, figsize = (10,10),)
    axes_3 = fig_3.subplot_mosaic(plot_layout, )
    fig_3.suptitle('Synapse current 3rd derivative', size = 20)
    for i in range(num_of_locations**2):
        driver_loc, mod_loc = locations_list[i]
        letter_digits = string.ascii_lowercase + string.digits
        ax_main = axes_main[letter_digits[i]]
        ax_1 = axes_1[letter_digits[i]]
        ax_2 = axes_2[letter_digits[i]]
        ax_3 = axes_3[letter_digits[i]]
        for weights in weights_list:
            current = gaussian_filter1d(-syn_currents[driver_loc,mod_loc][weights],5)

            right_idx = utils.get_index_at_half_value(current, eps = 0.0001)
            peaks = find_peaks(current)

            first_derivative, second_derivative, third_derivative = utils.get_derivative(-syn_currents[driver_loc,mod_loc][weights],)
            ax_1.plot(first_derivative, label = 'weight({:.2f})'.format(weights[1]))
            ax_1.get_shared_x_axes().join(ax_1, axes_1['a'])
            ax_1.get_shared_y_axes().join(ax_1, axes_1['a'])
            ax_2.plot(second_derivative, label = 'weight({:.2f})'.format(weights[1]))
            ax_2.get_shared_x_axes().join(ax_2, axes_2['a'])
            ax_2.get_shared_y_axes().join(ax_2, axes_2['a'])
            ax_3.plot(third_derivative, label = 'weight({:.2f})'.format(weights[1]))
            ax_3.get_shared_x_axes().join(ax_3, axes_3['a'])
            ax_3.get_shared_y_axes().join(ax_3, axes_3['a'])
            ax_main.plot(current, label = 'weight({:.2f})'.format(weights[1]))
            ax_main.get_shared_x_axes().join(ax_main, axes_main['a'])
            ax_main.get_shared_y_axes().join(ax_main, axes_main['a'])
            
            first_der_signchange = np.where(((np.roll(np.sign(first_derivative), 1) - np.sign(first_derivative)) != 0).astype(int))
            second_der_signchange = np.where(((np.roll(np.sign(second_derivative), 1) - np.sign(second_derivative)) != 0).astype(int))
            third_der_signchange = np.where(((np.roll(np.sign(third_derivative), 1) - np.sign(third_derivative)) != 0).astype(int))

            ax_main.scatter(right_idx, current[right_idx], marker = '*')

        ax_1.set_title('locations: {:.0f}um, {:.0f}um'.format(driver_loc*dend_L, mod_loc*dend_L), size = 12)
        ax_2.set_title('locations: {:.0f}um, {:.0f}um'.format(driver_loc*dend_L, mod_loc*dend_L), size = 12)
        ax_3.set_title('locations: {:.0f}um, {:.0f}um'.format(driver_loc*dend_L, mod_loc*dend_L), size = 12)
        ax_main.set_title('locations: {:.0f}um, {:.0f}um'.format(driver_loc*dend_L, mod_loc*dend_L), size = 12)
    ax_main.legend()
    ax_1.legend()
    ax_2.legend()
    ax_3.legend()
    close = True
    if close:
        plt.close(fig_1)      
        plt.close(fig_2)
        plt.close(fig_3)
    save = False
    if save:
        fig_folder = ''
        plt.savefig(fig_folder+'currents, 2 syn.png')

    
def plot_mean_vs_duration_scatters(syn_currents, voltages, syn_times, varying, columns, currents_to_plot = ['total'], dend_L = 200, stim_start = 110, ):
    '''
    plots a whole bunch (mean vs duration for current and v, syn currents, soma_v, dend_v,...) as subplots depending on locations and current clamp amp,
    as a function of weights or frequencies
    '''
    locations_list = list(syn_currents[currents_to_plot[0]])
    num_of_locations = len(locations_list)
    clamp_amps_list = list(syn_currents[currents_to_plot[0]][locations_list[0]])
    num_of_clamp_amps = len(clamp_amps_list)
    if columns == 'locations':
        num_of_columns = num_of_locations
    elif columns == 'clamp amps':
        num_of_columns = num_of_clamp_amps
    for current_to_plot in currents_to_plot:
        cutoff_type = 'width_at_half_peak'
        eps = 0.01
        if varying == 'weights':
            var_values_list = list(syn_currents[current_to_plot][locations_list[0]][clamp_amps_list[0]])
            num_of_var_values = len(var_values_list)
        elif varying == 'frequency':
            var_values_list = list(syn_currents[current_to_plot][locations_list[0]][clamp_amps_list[0]])
            num_of_var_values = len(var_values_list)
        if current_to_plot == 'AMPA':
            cutoff_type = 'trace_value'
        syn_c_means_durations, soma_v_means_durations = {}, {}
        for loc in locations_list:
            syn_c_means_durations[loc],soma_v_means_durations[loc] = {}, {}
            for amp in clamp_amps_list:
                syn_c_means_durations[loc][amp],soma_v_means_durations[loc][amp] = {}, {}
                for val in var_values_list:
                    syn_c_means_durations[loc][amp][val],soma_v_means_durations[loc][amp][val] = {}, {}
                    syn_current = -1*syn_currents[current_to_plot][loc][amp][val]
                    soma_v = voltages['soma'][loc][amp][val]
                    syn_c_right_index = utils.get_index_at_half_value(syn_current, cutoff_type = cutoff_type, eps = eps)
                    soma_v_right_index = utils.get_index_at_half_value(soma_v, cutoff_type = cutoff_type, eps = eps)
                    syn_c_duration, soma_v_duration = syn_c_right_index, soma_v_right_index
                    syn_c_mean, soma_v_mean = np.mean(syn_current[stim_start:syn_c_right_index]), np.mean(soma_v[stim_start:soma_v_right_index])
                    syn_c_means_durations[loc][amp][val][syn_c_duration], soma_v_means_durations[loc][amp][val][soma_v_duration] = syn_c_mean, soma_v_mean
        if num_of_columns == 3:
             plot_layout = """
             adg
             beh
             cfi
             """
        if num_of_columns == 4:
             plot_layout = """
             aeim
             bfjn
             cgko
             dhlp
             """
    
        syn_current_mean_duration_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        syn_current_mean_duration_axes = syn_current_mean_duration_fig.subplot_mosaic(plot_layout, )
        syn_current_mean_duration_fig.suptitle('Synapse current duration vs mean scatter - {}'.format(current_to_plot), size = 20)
        
        soma_v_mean_duration_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        soma_v_mean_duration_axes = soma_v_mean_duration_fig.subplot_mosaic(plot_layout, )
        soma_v_mean_duration_fig.suptitle('Soma_V duration vs mean scatter - {}'.format(current_to_plot), size = 20)
        
        currents_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        currents_axes = currents_fig.subplot_mosaic(plot_layout, )
        currents_fig.suptitle('Synapse currents - {}'.format(current_to_plot), size = 20)
        
        dend_v_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        dend_v_axes = dend_v_fig.subplot_mosaic(plot_layout, )
        dend_v_fig.suptitle('Dendrite voltages - {}'.format(current_to_plot), size = 20)
        
        soma_v_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        soma_v_axes = soma_v_fig.subplot_mosaic(plot_layout, )
        soma_v_fig.suptitle('Soma voltages - {}'.format(current_to_plot), size = 20)
        
        axes = [syn_current_mean_duration_axes,soma_v_mean_duration_axes,currents_axes,dend_v_axes,soma_v_axes]
        y_labels = ['mean','mean','current','dend_v','soma_v']
        x_labels = ['Duration','Duration','Time','Time','Time']
        
        letter_digits = string.ascii_lowercase + string.digits
        colors = np.linspace(0.3, 1, num_of_var_values)
        count = 0
        x = np.arange(0,100.1,0.1)
        spike_train_x = np.arange(0,100,1)
        spike_train = np.zeros(len(spike_train_x))
        if varying == 'weights':
            for i in syn_times[0]:
                spike_train[i] = 0.1
        cmap = matplotlib.cm.get_cmap('YlOrRd')
        for locs in range(num_of_locations):
            syn_loc, clamp_loc = locations_list[locs]
            [axi[letter_digits[count]].set_title('Synapse/clamp location = {:.0f}um'.format(syn_loc*dend_L), fontsize=16) for axi in axes]
            for amp in range(num_of_clamp_amps):
                clamp_amp = clamp_amps_list[amp]
                if count < num_of_clamp_amps:
                    [axes[i][letter_digits[count]].set_ylabel('Clamp amp = {:.1f}\n{}'.format(clamp_amp,y_labels[i]), fontsize=12) for i in range(len(axes))]
                else:
                    [axi[letter_digits[count]].set_yticklabels([]) for axi in axes]
                current_ax = currents_axes[letter_digits[count]]
                current_ax.get_shared_x_axes().join(current_ax, currents_axes['a'])
                current_ax.get_shared_y_axes().join(current_ax, currents_axes['a'])
                
                syn_current_mean_duration_ax = syn_current_mean_duration_axes[letter_digits[count]]
                syn_current_mean_duration_ax.get_shared_x_axes().join(syn_current_mean_duration_ax, syn_current_mean_duration_axes['a'])
                syn_current_mean_duration_ax.get_shared_y_axes().join(syn_current_mean_duration_ax, syn_current_mean_duration_axes['a'])
                
                soma_v_mean_duration_ax = soma_v_mean_duration_axes[letter_digits[count]]
                soma_v_mean_duration_ax.get_shared_x_axes().join(soma_v_mean_duration_ax, soma_v_mean_duration_axes['a'])
                soma_v_mean_duration_ax.get_shared_y_axes().join(soma_v_mean_duration_ax, soma_v_mean_duration_axes['a'])
                
                dend_v_ax = dend_v_axes[letter_digits[count]]
                dend_v_ax.get_shared_x_axes().join(dend_v_ax, dend_v_axes['a'])
                dend_v_ax.get_shared_y_axes().join(dend_v_ax, dend_v_axes['a'])
                
                soma_v_ax = soma_v_axes[letter_digits[count]]
                soma_v_ax.get_shared_x_axes().join(soma_v_ax, soma_v_axes['a'])
                soma_v_ax.get_shared_y_axes().join(soma_v_ax, soma_v_axes['a'])
                
                single_axes = [syn_current_mean_duration_ax,soma_v_mean_duration_ax,dend_v_ax,soma_v_ax,current_ax]
                syn_cur_durations, soma_v_durations = [], []
                syn_cur_means, soma_v_means = [], []
                if count == 6 and varying == 'weights':#set insert with spike train
                    for ax in [current_ax,dend_v_ax,soma_v_ax]:
                        ins_ax = inset_axes(ax,
                                    width="30%", # width = 30% of parent_bbox
                                    height=1., # height : 1 inch
                                    loc=1)
                        ins_ax.plot(spike_train_x, spike_train)
                        ins_ax.spines['right'].set_visible(False)
                        ins_ax.spines['top'].set_visible(False)
                        ins_ax.spines['bottom'].set_visible(False)
                        ins_ax.spines['left'].set_visible(False)
                        ins_ax.set_yticklabels([])
                        ins_ax.set_xticklabels(syn_times[0])
                        ins_ax.set_yticks([])
                        ins_ax.set_xticks(syn_times[0])
                for val in var_values_list:
                    if varying == 'weights':
                        label = '({:.2f})'.format(val)
                    elif varying == 'frequency':
                        label = '{}'.format(len(val))
                    color = [cmap(colors[var_values_list.index(val)])]
                    syn_c_duration = list(syn_c_means_durations[(syn_loc, clamp_loc)][clamp_amp][val])[0]#synapse current durations
                    syn_cur_durations.append(syn_c_duration)
                    syn_cur_means.append(syn_c_means_durations[(syn_loc, clamp_loc)][clamp_amp][val][syn_c_duration])
                    
                    soma_v_duration = list(soma_v_means_durations[(syn_loc, clamp_loc)][clamp_amp][val])[0]#synapse current durations
                    soma_v_durations.append(soma_v_duration)
                    soma_v_means.append(soma_v_means_durations[(syn_loc, clamp_loc)][clamp_amp][val][soma_v_duration])
 
                    dend_v = voltages['dend'][syn_loc, clamp_loc][clamp_amp][val]
                    dend_v_ax.plot(x, dend_v,  label = label, c = cmap(colors[var_values_list.index(val)]))
                    soma_v = voltages['soma'][syn_loc, clamp_loc][clamp_amp][val]
                    soma_v_ax.plot(x, soma_v,  label = label, c =cmap(colors[var_values_list.index(val)]))
                    soma_v_ax.scatter(soma_v_duration*0.1, soma_v[soma_v_duration], marker = '8', c = [cmap(colors[var_values_list.index(val)])])#*0.1 because of dt = 0.1
                    soma_v_ax.scatter(stim_start*0.1, soma_v[stim_start], marker = '8',c = [cmap(colors[var_values_list.index(val)])])
                    
                    current = -syn_currents[current_to_plot][syn_loc, clamp_loc][clamp_amp][val]
                    current_ax.plot(x, current,  label = label, c =cmap(colors[var_values_list.index(val)]))
                    current_ax.scatter(syn_c_duration*0.1, current[syn_c_duration], marker = '8',c = [cmap(colors[var_values_list.index(val)])])#*0.1 because of dt = 0.1
                    current_ax.scatter(stim_start*0.1, current[stim_start], marker = '8',c = [cmap(colors[var_values_list.index(val)])])
                syn_current_mean_duration_ax.scatter(np.array(syn_cur_durations)*0.1, syn_cur_means, c =colors, edgecolor = 'k', cmap = cmap,)
                soma_v_mean_duration_ax.scatter(np.array(soma_v_durations)*0.1, soma_v_means,c = colors, edgecolor = 'k', cmap = cmap,)
                if (count+1)%3 != 0:
                    [axi.set_xticklabels([]) for axi in single_axes]
                else:
                    [axes[i][letter_digits[count]].set_xlabel('{} [ms]'.format(x_labels[i]), fontsize=12) for i in range(len(axes))]  
                count += 1
        if varying == 'frequency':
            [axes[i+2]['a'].legend(title = 'spikes/100 ms', ncol=2) for i in range(len(axes)-2)]
        elif varying == 'weights':
            [axes[i+2]['a'].legend(title = 'weight', ncol=2) for i in range(len(axes)-2)]
            norm = plt.Normalize(var_values_list[0], var_values_list[-1])
            sm =  ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            syn_c_cbar = syn_current_mean_duration_fig.colorbar(sm,)  
            soma_v_cbar = soma_v_mean_duration_fig.colorbar(sm,)  
        plt.tight_layout()
        plt.show()

      
def plot_mean_vs_weights_scatters(syn_currents, voltages, syn_times, varying, columns, currents_to_plot = ['total'], dend_L = 200, stim_start = 110,
                                  ):
    '''
    like plot_mean_vs_durations_scatters from above, but now each column is either loc or clamp amp, rows are mean and duration, x_axis is weights, 
    curves are diff locations or clamp amps 
    '''
    cutoff_type = 'width_at_half_peak'
    eps = 0.01
    locations_list = list(voltages['soma'])
    num_of_locations = len(locations_list)
    clamp_amps_list = list(voltages['soma'][locations_list[0]])
    num_of_clamp_amps = len(clamp_amps_list)
    values = [locations_list,clamp_amps_list]
    if columns == 'locations':
        num_of_columns = num_of_locations
        curve_variable = values[1]#the difference between each curve. Can be diff loc or clamp amp values
        col_title = 'Syn/clamp loc = '
        legend_title = 'Clamp amp [nS]'
        col_title_values = [val[0]*dend_L for val in locations_list]
        col_title_units = 'um'
    elif columns == 'clamp amps':
        num_of_columns = num_of_clamp_amps
        curve_variable = values[0]
        col_title = 'Clamp amp = '
        legend_title = 'Syn/clamp loc [um] '
        title_values = clamp_amps_list
        col_title_values = clamp_amps_list
        col_title_units = 'nS'

    for current_to_plot in currents_to_plot:
        if varying == 'weights':
            var_values_list = list(syn_currents[current_to_plot][locations_list[0]][clamp_amps_list[0]])
            num_of_var_values = len(var_values_list)
        elif varying == 'frequency':
            var_values_list = list(syn_currents[current_to_plot][locations_list[0]][clamp_amps_list[0]])
            num_of_var_values = len(var_values_list)
        if current_to_plot == 'AMPA':
            cutoff_type = 'trace_value'
        if columns == 'locations':
            gain_threshold_x_label = 'Clamp amp [nS]'
            means_durations = {}#{loc1:{amp1:[[mean],[durations]]}}
            for loc in locations_list:
                means_durations[loc] = {}
                for amp in clamp_amps_list:
                    syn_c_means_vec, syn_c_durations_vec = [], []
                    for val in var_values_list:
                        syn_current = -1*syn_currents[current_to_plot][loc][amp][val]
                        soma_v = voltages['soma'][loc][amp][val]
                        syn_c_right_index = utils.get_index_at_half_value(syn_current, cutoff_type = cutoff_type, eps = eps)
                        soma_v_right_index = utils.get_index_at_half_value(soma_v, cutoff_type = cutoff_type, eps = eps)
                        syn_c_duration, soma_v_duration = syn_c_right_index, soma_v_right_index
                        syn_c_mean, soma_v_mean = np.mean(syn_current[stim_start:syn_c_right_index]), np.mean(soma_v[stim_start:soma_v_right_index])
                        syn_c_means_vec.append(syn_c_mean)
                        syn_c_durations_vec.append(syn_c_duration)
                    means_durations[loc][amp] = [syn_c_means_vec,syn_c_durations_vec]
        elif columns == 'clamp amps':
            gain_threshold_x_label = 'Location [um]'
            means_durations = {}#{amp1:{loc1:[[mean],[durations]]}}
            for amp in clamp_amps_list:
                means_durations[amp] = {}
                for loc in locations_list:
                    syn_c_means_vec, syn_c_durations_vec = [], []
                    for val in var_values_list:
                        syn_current = -1*syn_currents[current_to_plot][loc][amp][val]
                        soma_v = voltages['soma'][loc][amp][val]
                        syn_c_right_index = utils.get_index_at_half_value(syn_current, cutoff_type = cutoff_type, eps = eps)
                        soma_v_right_index = utils.get_index_at_half_value(soma_v, cutoff_type = cutoff_type, eps = eps)
                        syn_c_duration, soma_v_duration = syn_c_right_index, soma_v_right_index
                        syn_c_mean, soma_v_mean = np.mean(syn_current[stim_start:syn_c_right_index]), np.mean(soma_v[stim_start:soma_v_right_index])
                        syn_c_means_vec.append(syn_c_mean)
                        syn_c_durations_vec.append(syn_c_duration)
                    means_durations[amp][loc] = [syn_c_means_vec,syn_c_durations_vec]
        if num_of_columns == 3:
             plot_layout = """
             adg
             beh
             cfi
             """
        elif num_of_columns == 4:
             plot_layout = """
             adgj
             behk
             cfil
             """
        num_of_lines = len(curve_variable)
        cmap = matplotlib.cm.get_cmap('YlOrRd')
        colors = np.squeeze(np.linspace(0.3, 1, num_of_lines).reshape(1,-1))
    
        fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        axes = fig.subplot_mosaic(plot_layout, )
        fig.suptitle('Synapse current vs weights - {}'.format(current_to_plot), size = 20)
        
        gain_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        gain_axes = gain_fig.subplot_mosaic(plot_layout, )
        gain_fig.suptitle('Gain - {}'.format(current_to_plot), size = 20)
        
        threshold_fig = plt.figure(constrained_layout=True, figsize = (10,10),)
        threshold_axes = threshold_fig.subplot_mosaic(plot_layout, )
        threshold_fig.suptitle('Threshold - {}'.format(current_to_plot), size = 20)
        
        letter_digits = string.ascii_lowercase + string.digits
        count = 0
        for col in range(num_of_columns):
            syn_loc, clamp_loc = locations_list[col]
            clamp_amp = clamp_amps_list[col]
            
            mean_ax = axes[letter_digits[count]]
            duration_ax = axes[letter_digits[count+1]]
            integral_ax = axes[letter_digits[count+2]]
            gain_mean_ax = gain_axes[letter_digits[count]]
            gain_duration_ax = gain_axes[letter_digits[count+1]]
            gain_integral_ax = gain_axes[letter_digits[count+2]]
            threshold_mean_ax = threshold_axes[letter_digits[count]]
            threshold_duration_ax = threshold_axes[letter_digits[count+1]]
            threshold_integral_ax = threshold_axes[letter_digits[count+2]]
            
            mean_thresholds = []
            duration_thresholds = []
            integral_thresholds = []
            
            'create lists of axes'
            top_axes = [mean_ax, gain_mean_ax, threshold_mean_ax]
            bottom_axes = [integral_ax, gain_integral_ax, threshold_integral_ax]
            main_row_axes = [mean_ax, duration_ax, integral_ax]
            gain_row_axes = [gain_mean_ax, gain_duration_ax, gain_integral_ax]
            threshold_row_axes = [threshold_mean_ax, threshold_duration_ax, threshold_integral_ax]
            
            'set title on columns'
            [ax.set_title(col_title+'{:.1f}'.format(col_title_values[col])+col_title_units, fontsize=16) for ax in top_axes]
            
            'share all axes with first column'
            [main_row_axes[i].get_shared_y_axes().join(main_row_axes[i], axes[letter_digits[i]]) for i in range(len(main_row_axes))]
            [gain_row_axes[i].get_shared_y_axes().join(gain_row_axes[i], gain_axes[letter_digits[i]]) for i in range(len(gain_row_axes))]
            [threshold_row_axes[i].get_shared_y_axes().join(threshold_row_axes[i], threshold_axes[letter_digits[i]]) for i in range(len(threshold_row_axes))]
            [main_row_axes[i].get_shared_x_axes().join(main_row_axes[i], axes[letter_digits[i]]) for i in range(len(main_row_axes))]
            [gain_row_axes[i].get_shared_x_axes().join(gain_row_axes[i], gain_axes[letter_digits[i]]) for i in range(len(gain_row_axes))]
            [threshold_row_axes[i].get_shared_x_axes().join(threshold_row_axes[i], threshold_axes[letter_digits[i]]) for i in range(len(threshold_row_axes))]

            y_labels = ['Mean', 'Duration [ms]', 'Integral']
             
            if count == 0:
                [main_row_axes[i].set_ylabel(y_labels[i]) for i in range(len(y_labels))]
                [gain_row_axes[i].set_ylabel(y_labels[i]) for i in range(len(y_labels))]
                [threshold_row_axes[i].set_ylabel(y_labels[i]) for i in range(len(y_labels))]
            else:
                [ax.set_yticklabels([]) for ax in main_row_axes]
                [ax.set_yticklabels([]) for ax in gain_row_axes]
                [ax.set_yticklabels([]) for ax in threshold_row_axes]
            [ax.set_xticklabels([]) for ax in main_row_axes[:-1]]
            [ax.set_xticklabels([]) for ax in gain_row_axes[:-1]]
            [ax.set_xticklabels([]) for ax in threshold_row_axes[:-1]]
            
            integral_ax.set_xlabel('Weights')
            gain_integral_ax.set_xlabel(gain_threshold_x_label)
            threshold_integral_ax.set_xlabel(gain_threshold_x_label)
            gain_x_values = []
            for val in curve_variable:#loc or clamp amp
                color_count = curve_variable.index(val)
                if columns == 'locations':
                    means = means_durations[(syn_loc, clamp_loc)][val][0]
                    durations = 0.1*np.array(means_durations[(syn_loc, clamp_loc)][val][1])
                    gain_x_value = clamp_amp
                elif columns == 'clamp amps':
                    means = means_durations[clamp_amp][val][0]
                    durations = 0.1*np.array(means_durations[clamp_amp][val][1])
                    val = val[0]*dend_L
                    gain_x_value = val
                integral = means*durations
                gain_x_values.append(gain_x_value)
                mean_derivs,_,_ = utils.get_derivative(means[1:],)
                duration_derivs,_,_ = utils.get_derivative(durations[1:],)
                integral_derivs,_,_ = utils.get_derivative(integral[1:])
                
                mean_threshold_index = utils.find_threshold(means[1:])+1
                duration_threshold_index = utils.find_threshold(durations[1:])+1
                integral_threshold_index = utils.find_threshold(integral[1:])+1
                
                mean_thresholds.append(var_values_list[mean_threshold_index])
                duration_thresholds.append(var_values_list[duration_threshold_index])
                integral_thresholds.append(var_values_list[integral_threshold_index])
                
                plot_type = 'lines'
                if plot_type == 'lines':
                    point_size = 2.5
                    color = cmap(colors[color_count])#for plot
                    mean_ax.plot(var_values_list, means, '--bo', label = '{:.2f}'.format(val), c = color, markersize = point_size)
                    mean_ax.scatter(var_values_list[mean_threshold_index], means[mean_threshold_index])
                    duration_ax.plot(var_values_list, durations, '--bo', label = '{:.2f}'.format(val), c = color, markersize = point_size)
                    duration_ax.scatter(var_values_list[duration_threshold_index], durations[duration_threshold_index])
                    integral_ax.plot(var_values_list, integral, '--bo', label = '{:.2f}'.format(val), c = color, markersize = point_size)
                    integral_ax.scatter(var_values_list[integral_threshold_index], integral[integral_threshold_index])                    

                elif plot_type == 'scatter':
                    point_size = 5
                    color = np.array([cmap(colors[color_count])])#for scatter
                    mean_ax.scatter(var_values_list, means, label = '{:.2f}'.format(val), c = color, s= point_size)
                    duration_ax.scatter(var_values_list, durations, label = '{:.2f}'.format(val), c = color, s= point_size)
                    integral_ax.scatter(var_values_list, means*durations, label = '{:.2f}'.format(val), c = color, s= point_size)
                gain_mean_ax.scatter(gain_x_value, means[-1], c = 'k')
                gain_duration_ax.scatter(gain_x_value, durations[-1], c = 'k')
                gain_integral_ax.scatter(gain_x_value, integral[-1], c = 'k')
            
            threshold_mean_ax.scatter(gain_x_values, mean_thresholds, c= 'k')
            threshold_duration_ax.scatter(gain_x_values, duration_thresholds, c= 'k')
            threshold_integral_ax.scatter(gain_x_values, integral_thresholds, c= 'k')
            count += 3
        
        mean_ax.legend(title = legend_title, ncol=2)
        
    plt.tight_layout()
    plt.show()
   
    
def plot_ratios(cell):
    locs = np.linspace(cell.loc_2_dist(cell.dend(0)),cell.loc_2_dist(cell.dend(1)),len(cell.soma_TRs))
    fig, ax = plt.subplots()
    ax.scatter(locs, cell.soma_ratios)
    ax.set_title('Transfer voltage ratio (Mega Ohms)')
    ax.set_ylabel('Ratio')
    ax.set_xlabel('Distance from soma (microns)')
    plt.show()
    
    
def plot_dendrite_heatmap(cell, synapses = None, annotate = False, dt = 0.1, section_to_plot = 0, dendrite_or_spines = 'dendrite', dend_loc_to_plot = 0.5,
                          plot_spine_v_traces = True, active_loc = 0.5, plot = True):
    if dendrite_or_spines == 'dendrite':
        voltage_matrix = np.vstack((cell.soma_v, cell.dend_vs[section_to_plot]))
    
    if plot:
        fig = plt.figure()#figsize = (4,4), dpi = 300)
        im = plt.imshow(voltage_matrix, origin='lower', cmap=plt.get_cmap('inferno'), vmin = -77, vmax = 0)
        ax = plt.gca()
        x_positions = np.arange(0,len(cell.soma_v),100)
        x_labels = np.int0(np.linspace(0,dt*len(cell.t),len(x_positions)))
        plt.xticks(x_positions[::3], x_labels[::3])
        ax.set_aspect(25)
        r, c = np.shape(voltage_matrix)
        ax.set_yticks(np.arange(0, r, 1)[::2])
        ylabels = list(np.linspace(0,cell.dend.L,r))
        ylabels = [round(i) for i in ylabels]
        ylabels[0] = 'soma'
        spikes = synapses[list(synapses)[0]][list(synapses[list(synapses)[0]])[0]][list(synapses[list(synapses)[0]][list(synapses[list(synapses)[0]])[0]])[0]]
        freq = '1 spike'
        neck_resistance = utils.calculate_resistance(cell.necks[0].Ra, cell.necks[0].diam, cell.necks[0].L)
        if len(spikes) > 1:
            freq = np.round(1000/(spikes[1]-spikes[0]))
        
        title = ('\nSyn branches: {}, syn loc: {}, # of active syns: {}, \n# of stims: {}, theta_d: {}, theta_p: {}, freq: {}hz, \nNeck resistance = {}MOhm, neck_diam = {}um'
                 .format(list(synapses),
                                                                                                      list(synapses[section_to_plot]),
                                                                                                      len(list(synapses[section_to_plot][active_loc])),
                                                                                                      len(list(synapses[section_to_plot][active_loc][1])),
                                                                                                      cell.e_synapses[0].theta_d_GB,cell.e_synapses[0].theta_p_GB,
                                                                                                      freq, neck_resistance, cell.necks[0].diam))
        ax.set_title('Voltage heatmap '+ title)
        ax.set_yticklabels(ylabels[::2])
        ax.set_ylabel('Distance from soma [um]')
        ax.set_xlabel('Time [ms]')
        plt.colorbar() 
        
        if annotate:
            for loc in list(synapses):
                for weight in list(synapses[loc]):
                    time = synapses[loc][weight][0]/dt
                    text1 = ax.text(time, loc*r, weight,
                                       ha="center", va="center",
                                       color = 'y')   
        'comment in block for traces'
        fig, ax = plt.subplots()#figsize = (5,4), dpi = 300)
        ax.set_title('Voltages' + title )

        plt.plot(cell.t, cell.soma_v, label = 'Soma', c = 'k')
        seg_num = round(dend_loc_to_plot*len(cell.dend_vs[section_to_plot]))
        v = cell.dend_vs[section_to_plot][seg_num]
        plt.plot(cell.t, v, label = 'dend '+str(int(dend_loc_to_plot*cell.dend.L))+'um', ls ='--', zorder = 2,)

        if plot_spine_v_traces:
            plt.plot(cell.t, cell.vsyn[0], label = 'syn 60um', c = 'g')
            plt.plot(cell.t, cell.vsyn[1], label = 'syn '+str(int(dend_loc_to_plot*cell.dend.L))+'um', c = 'r',)
            plt.plot(cell.t, cell.vsyn[-1], label = 'syn 200um', c = 'b', alpha = 0.4, zorder = 1, )#linewidth = 4)
      
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('Time [ms]')
        ax.set_ylabel('Voltage [mV]')
        ax.set_ylim([-79,0])
          
        'inset axis'
    #     left, bottom, width, height = [0.7, 0.35, 0.3, 0.3]
    #     inset_ax = ax.inset_axes([left, bottom, width, height],)
    #     inset_ax.spines['top'].set_visible(False)
    #     inset_ax.spines['right'].set_visible(False)
    #     inset_ax.set_xticks([])
    #     inset_ax.set_xticklabels([])
    #     inset_ax.set_yticks([])
    #     inset_ax.set_yticklabels([]) 
    #     inset_ax.plot(cell.t, v, ls ='--', c = 'r')
    #     inset_ax.plot(cell.t, cell.vsyn[1], c = 'r')
    # #     inset_ax.plot(cell.t, cell.dend_vs[section_to_plot][-1], c = 'b')    
    #     inset_ax.plot(cell.t, cell.vsyn[-1], c = 'b')
    #     inset_ax.set_xlim([11.5, 16])  
    #     inset_ax.set_ylim([-29, -13])
          
        plt.legend()
    return voltage_matrix


def plot_dendrite_heatmap_hay(cell, synapses = None, annotate = False, dt = 0.1, section_to_plot = 0, dendrite_or_spines = 'dendrite', dend_loc_to_plot = 0.5,
                          plot_spine_v_traces = True, active_loc = 0.5, plot = True):
    
    if dendrite_or_spines == 'dendrite':
        print(np.shape(cell.dend_vs[0]))
        voltage_matrix = np.vstack((cell.soma_v, cell.dend_vs[section_to_plot]))
    
    if plot:
        fig = plt.figure(figsize = (4,4), dpi = 300)
        im = plt.imshow(voltage_matrix, origin='lower', cmap=plt.get_cmap('inferno'), vmin = -77, vmax = 0)
        ax = plt.gca()
        x_positions = np.arange(0,len(cell.soma_v),100)
        x_labels = np.int0(np.linspace(0,dt*len(cell.t),len(x_positions)))
        plt.xticks(x_positions[::3], x_labels[::3])
        ax.set_aspect(25)
        r, c = np.shape(voltage_matrix)
        ax.set_yticks(np.arange(0, r, 1)[::2])
        ylabels = list(np.linspace(0,cell.dend.L,r))
        ylabels = [round(i) for i in ylabels]
        ylabels[0] = 'soma'
        spikes = synapses[list(synapses)[0]][list(synapses[list(synapses)[0]])[0]][list(synapses[list(synapses)[0]][list(synapses[list(synapses)[0]])[0]])[0]]
        freq = '1 spike'
        neck_resistance = utils.calculate_resistance(cell.necks[0].Ra, cell.necks[0].diam, cell.necks[0].L)
        if len(spikes) > 1:
            freq = np.round(1000/(spikes[1]-spikes[0]))
        
        title = ('\nSyn branches: {}, syn loc: {}, # of active syns: {}, \n# of stims: {}, theta_d: {}, theta_p: {}, freq: {}hz, \nNeck resistance = {}MOhm, neck_diam = {}um'
                 .format(list(synapses),
                                                                                                      list(synapses[section_to_plot]),
                                                                                                      len(list(synapses[section_to_plot][active_loc])),
                                                                                                      len(list(synapses[section_to_plot][active_loc][1])),
                                                                                                      cell.e_synapses[0].theta_d_GB,cell.e_synapses[0].theta_p_GB,
                                                                                                      freq, neck_resistance, cell.necks[0].diam))
        ax.set_title('Voltage heatmap '+ title)
        ax.set_yticklabels(ylabels[::2])
        ax.set_ylabel('Distance from soma [um]')
        ax.set_xlabel('Time [ms]')
        plt.colorbar() 
        
        if annotate:
            for loc in list(synapses):
                for weight in list(synapses[loc]):
                    time = synapses[loc][weight][0]/dt
                    text1 = ax.text(time, loc*r, weight,
                                       ha="center", va="center",
                                       color = 'y')   
        'comment in block for traces'
        fig, ax = plt.subplots(figsize = (5,4), dpi = 300)
        ax.set_title('Voltages' + title )

        plt.plot(cell.t, cell.soma_v, label = 'Soma', c = 'k')
        seg_num = round(dend_loc_to_plot*len(cell.dend_vs[section_to_plot]))
        v = cell.dend_vs[section_to_plot][seg_num]
        plt.plot(cell.t, v, label = 'dend '+str(int(np.round(dend_loc_to_plot*cell.dend.L)))+'um', ls ='--', zorder = 2,)

        if plot_spine_v_traces:
            plt.plot(cell.t, cell.vsyn[0], label = 'syn 60um', c = 'g')
            plt.plot(cell.t, cell.vsyn[1], label = 'syn '+str(int(np.round(dend_loc_to_plot*cell.dend.L)))+'um', c = 'r',)
            plt.plot(cell.t, cell.vsyn[-1], label = 'syn 200um', c = 'b', alpha = 0.4, zorder = 1, )#linewidth = 4)
      
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('Time [ms]')
        ax.set_ylabel('Voltage [mV]')
        ax.set_ylim([-79,0])
          
        'inset axis'
    #     left, bottom, width, height = [0.7, 0.35, 0.3, 0.3]
    #     inset_ax = ax.inset_axes([left, bottom, width, height],)
    #     inset_ax.spines['top'].set_visible(False)
    #     inset_ax.spines['right'].set_visible(False)
    #     inset_ax.set_xticks([])
    #     inset_ax.set_xticklabels([])
    #     inset_ax.set_yticks([])
    #     inset_ax.set_yticklabels([]) 
    #     inset_ax.plot(cell.t, v, ls ='--', c = 'r')
    #     inset_ax.plot(cell.t, cell.vsyn[1], c = 'r')
    # #     inset_ax.plot(cell.t, cell.dend_vs[section_to_plot][-1], c = 'b')    
    #     inset_ax.plot(cell.t, cell.vsyn[-1], c = 'b')
    #     inset_ax.set_xlim([11.5, 16])  
    #     inset_ax.set_ylim([-29, -13])
          
        plt.legend()
    return voltage_matrix


def plot_multiple_dendrite_heatmaps(cell, annotate = False, dt = 0.1, sec_numbers = [0], sim_time = 100, active_syn_protocol = 'NA', voltage_trace = False,
                                    branch_num = 4):
    num_of_rows = np.int0(np.log2(len(cell.dendrites)))+1
    num_of_columns = 2**num_of_rows - 1
    voltage_matrices = []
    
    fig, axes = plt.subplots(num_of_rows, num_of_columns)#, dpi = 400, figsize = (24,12),)
    synapses = cell.synapses
    active_branches = list(synapses)
    active_locations = [(i,list(synapses[i])) for i in synapses]
    num_of_active_syns = len(list(synapses[active_branches[0]][active_locations[0][1][0]]))
    title = ('\nActive locations: {}, # of active syns - prox/dist: {}, \ntheta_d: {}, theta_p: {}'.format(active_locations,
                                                                                                  active_syn_protocol,
                                                                                                  cell.e_synapses[0].theta_d_GB,cell.e_synapses[0].theta_p_GB))
    fig.suptitle(title)
    plotted_axes = []
    x_positions = np.linspace(0,len(cell.soma_v),3)    
    x_labels = np.int0(np.linspace(0,sim_time,len(x_positions)))
    fig.text(0.5,0.01, "Time [ms]", ha="center", va="center", size = 15)

    for sec_num in range(len(cell.dend_vs)):
        if sec_num == 0:
            voltage_matrix = np.vstack((cell.soma_v, cell.dend_vs[sec_num]))
        else:
            voltage_matrix = np.vstack(cell.dend_vs[sec_num])
        voltage_matrices.append(voltage_matrix)
        
        bottom_up_row = np.int0(np.floor(np.log2(sec_num+1)))#starting to count from the bottom
        
        col = np.int0(np.floor(num_of_columns/2**(bottom_up_row+1))+(sec_num-2**bottom_up_row+1)*2**(num_of_rows-bottom_up_row))
        row = num_of_rows - 1 - bottom_up_row
        if num_of_rows == 1:
            ax = axes
        else:
            ax = axes[row,col]
        im = ax.imshow(voltage_matrix, origin='lower', cmap=plt.get_cmap('inferno'), vmin = -77, vmax = 0)
        plotted_axes.append((row,col))
        ax.set_aspect('auto')
        r, c = np.shape(voltage_matrix)
        ylabels = list(np.linspace(0,cell.dendrites[sec_num].L,len(np.arange(0, r, 1))))
        ylabels = [round(i) for i in ylabels]
 
        if sec_num == 0:
            
            ylabels[0] = 'soma'
            ax.set_xticks(x_positions)
            ax.set_xticklabels(x_labels)
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])
        ax.yaxis.set_label_position("right")
        
        ax.set_yticks(np.arange(0, r, 1))
        if sec_num == 2**bottom_up_row - 1: 
            ax.set_yticklabels([])
        else:    
            ax.set_yticklabels([])
          
    for r in range(num_of_rows):
        for c in range(num_of_columns):
             if (r,c) not in plotted_axes:
                 axes[r,c].axis('off')

    if voltage_trace:
        if branch_num == 1:
            left, bottom, width, height = [0.35, 0.34, 0.1, 0.1]
            loc = (1,0.5)
        elif branch_num == 3:
            left, bottom, width, height = [0.245, 0.54, 0.1, 0.1]
            loc = (3,0.5)
        elif branch_num == 4:
            left, bottom, width, height = [0.445, 0.54, 0.1, 0.1]
            loc = (4,0.5)
        elif branch_num == 7:
            left, bottom, width, height = [0.02, 0.745, 0.1, 0.1]#[0.45, 0.53, 0.1, 0.1]
            loc = (7,0.5)
            
        inset_ax = fig.add_axes([left, bottom, width, height],)
        inset_ax.spines['top'].set_visible(False)
        inset_ax.spines['left'].set_linewidth(4)
        inset_ax.spines['right'].set_visible(False)
        inset_ax.spines['bottom'].set_linewidth(4)
        inset_ax.set_xticks([])
        inset_ax.set_xticklabels([])
        inset_ax.set_yticks([])
        inset_ax.set_yticklabels([]) 
        seg_num = round(loc[1]*len(cell.dend_vs[loc[0]]))
        v = cell.dend_vs[loc[0]][seg_num]
        inset_ax.plot(cell.t, v, c = 'k', linewidth = 4)
        inset_ax.set_ylim([-80,0])

    return voltage_matrices

     
def plot_conductances(cell, rec_vecs, save_fig = False, protocol = None, active_loc = (0, 0.5), fig = None, num_of_traces = (1, 20, 0)):
    
    num_of_synapses = len(cell.e_synapses)
    branches = list(cell.synapses)
    num_of_figs = len(branches)
    synapses_per_branch = np.int0(np.zeros(len(branches)))
    unique_synapses_per_branch = np.zeros(len(branches))
    start_idx = 0
    'inset coordinates'
    left, bottom, width, height = [0.625, 0.3, 0.35, 0.6]
    colors = ['k','b','b','b','g']

    colors = plt.cm.autumn(np.linspace(0, 1, num_of_traces[0])[::-1])
    len_of_x_display = 400
    cell_t_short = npa(cell.t)[:len_of_x_display]
    inset_ax_max = np.max(cell.effcai_GB)
    for branch_num in range(len(branches)):
        branch = branches[branch_num]
        unique_synapses_per_branch[branch_num] += len(list(cell.synapses[branch]))
        for loc in list(cell.synapses[branch]):
            synapses_per_branch[branch_num] += len(list(cell.synapses[branch][loc]))
    for branch_num in range(len(branches)):
        branch = branches[branch_num]
        if fig:
            pass
        else:
            fig = plt.figure()#figsize = (15,10),)
            
        if len(cell.synapses[active_loc[0]][active_loc[1]]) == num_of_traces[1]:
            ls = '-'
        else:
            ls = '--'
        
        num_of_synapses_to_plot = np.int0(unique_synapses_per_branch[branch_num])
        
        unique_synapse_indices = np.unique(cell.e_locations[start_idx:start_idx + synapses_per_branch[branch_num]], return_index=True)[1] + start_idx

        start_idx += synapses_per_branch[branch_num]
        if len(cell.synapses[active_loc[0]][active_loc[1]][1]) < 2:
            fig.suptitle('branch #: '+str(branch)+', dt: '+str(cell.dt)+', # of stims: '+str(len(cell.synapses[active_loc[0]][active_loc[1]][1]))+', freq: '
                     +str(1000/(1))+'Hz, # of active synapses: '+str(len(cell.synapses[active_loc[0]][active_loc[1]]))+', active synapse loc: '+str(active_loc)
                     + ', theta_d: '+ str(cell.e_synapses[0].theta_d_GB) + ', theta_p: ' + str(cell.e_synapses[0].theta_p_GB) 
                     + ', sim length: ' + str(np.int0(len(cell.t)*cell.dt)))
        else:
            fig.suptitle('branch #: '+str(branch)+', dt: '+str(cell.dt)+', # of stims: '+str(len(cell.synapses[active_loc[0]][active_loc[1]][1]))+', freq: '
                     +str(1000/(cell.synapses[active_loc[0]][active_loc[1]][1][1]-cell.synapses[active_loc[0]][active_loc[1]][1][0]))
                     +'Hz, # of active synapses: '+str(len(cell.synapses[active_loc[0]][active_loc[1]]))+', active synapse loc: '+str(active_loc) 
                     + ', theta_d: '+ str(cell.e_synapses[0].theta_d_GB) + ', theta_p: ' + str(cell.e_synapses[0].theta_p_GB) 
                     + ', sim length: ' + str(np.int0(len(cell.t)*cell.dt)))
        axes_list = []
        plot_num = 1
        for row_num in range(1,len(rec_vecs)+1):
            trace_name = rec_vecs[row_num-1]
            syn_locations_frac = [0.3,0.45,0.5,0.7]
            for syn_num in range(num_of_synapses_to_plot):
                if syn_num == 0: 
                    ax = plt.subplot(len(rec_vecs),num_of_synapses_to_plot,plot_num,)
                    ax.set_ylabel(trace_name, size = 15)
                else:
                    if row_num == 4:
                        ax = plt.subplot(len(rec_vecs),num_of_synapses_to_plot,plot_num, sharex = axes_list[-1],)# sharey = axes_list[-1])
                    else:
                        ax = plt.subplot(len(rec_vecs),num_of_synapses_to_plot,plot_num, sharex = axes_list[-1], sharey = axes_list[-1])
                        plt.setp(ax.get_yticklabels(), visible=False)             
                axes_list.append(ax)
                syn_loc = list(cell.synapses[branch])[syn_num]
                if row_num in [1,2,3]: 
                    if row_num == 1:
                        ax.title.set_text('Syn loc = '+str(int(syn_loc*cell.dendrites[0].L))+'um')#\nActivations (weight:stim_times) = '+str(cell.synapses[0][syn_loc]))

                factor = 1
                if trace_name[0] == 'i': factor = -1
                if row_num == 5:
                    vec1 = eval('factor*cell.'+trace_name+'[unique_synapse_indices[syn_num]]')
                    vec2 = eval('factor*cell.'+trace_name+'[unique_synapse_indices[syn_num]]')
                    trace_a = eval('factor*cell.'+trace_name+'[unique_synapse_indices[syn_num]]')
                    long_trace = vec1.resample(vec1,0.01)
                    trace = npa(vec2.c(0,len_of_x_display))[:-1]

                else:
                    long_trace = eval('factor*cell.'+trace_name+'[unique_synapse_indices[syn_num]]')
                    trace = npa(long_trace.c(0,len_of_x_display))[:-1]

                if row_num in [1,2,3,]:
                    ax.plot(cell_t_short, trace, ls = ls, c = 'k')# c = colors[num_of_traces[2]])#'c = colors[row_num-1])# label = trace_name)

                if row_num == 5: 
                    ax.plot(cell.t, trace_a, ls = ls)

                    ax.set_xlabel('dt [ms]'), 
                    ax.axhline(0.5, alpha = 0.25, ls = '--')
                    inset_ax = ax.inset_axes([left, bottom+0.4, width, height/2],)
                    inset_ax.spines['top'].set_visible(False)
                    inset_ax.spines['right'].set_visible(False)
                    inset_ax.plot(npa(cell.t)[::100][:-1], long_trace, )#c = colors[row_num-1])
                    inset_ax.set_xticklabels([])
                    inset_ax.set_yticks([])
                    inset_ax.set_yticklabels([])
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                if trace_name == 'effcai_GB':
                    if len(cell.synapses[active_loc[0]][active_loc[1]]) == num_of_traces[1]:
                        plt.axhline(cell.e_synapses[0].theta_p_GB, alpha = 1, ls = '--', c = 'r', label = r'$\theta_p$', zorder = 2)
                        plt.axhline(cell.e_synapses[0].theta_d_GB, alpha = 1, ls = '--', c = 'b', label = r'$\theta_d$', zorder = 1)

                    ax.plot(cell.t, long_trace, ls = ls, c = 'k')#'c = colors[row_num-1]),#inset_ax.set_xlim([0,1000])

                    if syn_num == 0:
                        plt.legend(loc = 'best')#loc = 'upper left')
                plot_num += 1

    plt.tight_layout()    
    return fig
  
        
def plot_heterosynaptic_plasticity(): 
    '''
    bar chart
    '''
    fig = plt.figure(dpi = 300)
    bar_values = [[2.5,2.5],[2.5,1.5],[2.5,0.5],[1.5,2.5],[1.5,1.5],[1.5,0.5],[0.5,2.5],[0.5,1.5],[0.5,0.5]]
    for plot_num in range(9):
        ax = plt.subplot(3,3,plot_num+1,)

        ax.bar([1,2],bar_values[plot_num], color = ['dimgray','lightgray'])
        ax.set_xticklabels([])
        ax.set_yticks([1,2])
        if plot_num in [0,3,6]: ax.set_yticklabels([r'$\theta_D$',r'$\theta_P$'])
        else: ax.set_yticklabels([])
        ax.set_ylim([0,3])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.axhline(1, alpha = 0.5, ls = '--', c = 'b', )
        ax.axhline(2, alpha = 0.5, ls = '--', c = 'r', )
        if plot_num in [3,6,7,]:
            ax.text(0.7,0,'x', size = 140, alpha = 0.3)
    return 

    
def plot_plasticity_tree(cell, fig, axes, fig_num, rows = 2, cols = 2, col_title = None, x = 5, y = 1, syn_size = 50, soma_size = 60):
    num_of_branches = cell.num_of_sections
    num_of_branch_levels = int(np.log2(num_of_branches) + 1)
    points_for_tree = [[[x,y]],]
    for i in range(1, num_of_branch_levels):
        points_for_tree.append([])
        for syn in range(2**(i-1)):
            factor = 2**i
            points_for_tree[i].append([points_for_tree[i-1][syn][0]-2/factor,points_for_tree[i-1][syn][1]+1])
            points_for_tree[i].append([points_for_tree[i-1][syn][0]+2/factor,points_for_tree[i-1][syn][1]+1])
    points_for_tree.insert(0, [[x,y-1]])

    'building the tree'
    ax = fig.add_subplot(rows,cols,fig_num)
    for branch_lvl in range(1,len(points_for_tree)):
        for point_num in range(len(points_for_tree[branch_lvl])):
            point_0 = points_for_tree[branch_lvl-1][int(np.floor(point_num/2))]
            point_1 = points_for_tree[branch_lvl][point_num]
            plt.plot([point_0[0],point_1[0]],[point_0[1],point_1[1]], 'k', zorder = 1)
    
    ax.set_xticks([])
    
    ylabels = list(range(0,int(num_of_branch_levels*cell.dendrites[0].L) +1, int(cell.dendrites[0].L)))
    ylabels[0] = 'soma'
    
    plt.scatter(x, 0, marker = 's', zorder = 2, s = soma_size, c = 'k')
    if (fig_num - 1)%cols == 0:
        ax.set_yticks(np.arange(num_of_branch_levels + 1))
        ax.set_yticklabels(ylabels)
    else:
        ax.set_yticks([])
        ax.set_yticklabels([])

    if col_title:
        ax.set_title('# of active synapses: '+str(col_title))
        
    point_list = []
    for branch in points_for_tree:
        [point_list.append(point) for point in branch]
    
    'scattering the_synapses'
    x_locs, y_locs = [],[]
    branches = list(cell.synapses)
    synapses_per_branch = np.int0(np.zeros(len(branches)))
    start_idx = 0
    plot_syn_indices = []
    for branch_num in range(len(branches)):
        branch = branches[branch_num]
        
        x_start, y_start = point_list[np.int0(np.ceil(branch/2))]
        x_end, y_end = point_list[branch + 1]
         
        for syn_num in range(len(cell.synapses[branch])):
            fraction_loc = list(cell.synapses[branch])[syn_num]
            x_loc = x_start + fraction_loc*(x_end-x_start)
            y_loc = y_start + fraction_loc*(y_end-y_start)
            x_locs.append(x_loc)
            y_locs.append(y_loc)
            
            synapses_per_branch[branch_num] += len(list(cell.synapses[branch][fraction_loc]))
        unique_synapse_indices = np.unique(cell.e_locations[start_idx:start_idx + synapses_per_branch[branch_num]], return_index=True)[1] + start_idx
        start_idx += synapses_per_branch[branch_num]
        plot_syn_indices.extend(unique_synapse_indices)
    
    color_vals = [np.heaviside(cell.rho_GB[i][-1]-0.5,0.5) for i in plot_syn_indices]
    shapes = ['o' for i in range(len(plot_syn_indices))]
    count = 0
    
    cmap = cm.get_cmap('coolwarm', 3)
    colors = cmap(color_vals)
    for i in range(len(x_locs)):
        scatter(x_locs[i], y_locs[i], s = syn_size, color=colors[i], marker = shapes[i], edgecolor = 'k', cmap=cmap, zorder = 2)
    if fig_num == 1:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
        cbar = plt.colorbar(cm.ScalarMappable(cmap=cmap), cax=cbar_ax)
        plt.clim(0,1)
        cbar.set_ticks([1/6, 3/6, 5/6])   
        cbar.set_ticklabels(['depression','no change','potentiation'])
    return fig, axes


def make_cho_plot():
    '''
    function to make specific plot from Cho et. al. (2001) for Toviah presentation
    "An experimental test of the role of postsynaptic calcium levels in determining synaptic strength using perirhinal cortex of rat"
    '''
    fig_num = 2
    if fig_num == 2:
        x_vec = [0.2, 0.5, 5, 10, 20]
        y_vec = [43, 40, -7, -40, 6]
    
    fig, ax = plt.subplots()
    plt.plot(x_vec, y_vec, '-ko')
    ax.invert_xaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('[EGTA](mM)')
    ax.set_ylabel('Peak EPSC \n(baseline % change)')
    plt.axhline(ls = '--', alpha = 0.4)
    plt.tight_layout()
    
    
def plot_colored_branches(cell, fig, axes, fig_num, rows = 2, cols = 2, col_title = None, x = 5, y = 1,):# active_syn_index = -1):
    num_of_branches = cell.num_of_sections
    num_of_branch_levels = int(np.log2(num_of_branches) + 1)
    points_for_tree = [[[x,y]],]
    for i in range(1, num_of_branch_levels):
        points_for_tree.append([])
        for syn in range(2**(i-1)):
            factor = 2**i
            points_for_tree[i].append([points_for_tree[i-1][syn][0]-2/factor,points_for_tree[i-1][syn][1]+1])
            points_for_tree[i].append([points_for_tree[i-1][syn][0]+2/factor,points_for_tree[i-1][syn][1]+1])
    points_for_tree.insert(0, [[x,y-1]])

    'building the tree'
    ax = fig.add_subplot(rows,cols,fig_num)
    
    for branch_lvl in range(1,len(points_for_tree)):
        for point_num in range(len(points_for_tree[branch_lvl])):
            point_0 = points_for_tree[branch_lvl-1][int(np.floor(point_num/2))]
            point_1 = points_for_tree[branch_lvl][point_num]
            plt.plot([point_0[0],point_1[0]],[point_0[1],point_1[1]], 'k', zorder = 1)
    
    ax.set_xticks([])
    
    ylabels = list(range(0,int(num_of_branch_levels*cell.dendrites[0].L) +1, int(cell.dendrites[0].L)))
    ylabels[0] = 'soma'
    
    if fig_num == 1:
        ax.set_yticks(np.arange(num_of_branch_levels + 1))
        ax.set_yticklabels(ylabels)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    else:
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
    if col_title:
        ax.set_title('# of active synapses: '+str(col_title))
        
    point_list = []
    for branch in points_for_tree:
        [point_list.append(point) for point in branch]
    
    'scattering the_synapses'
    x_locs, y_locs = [],[]
    branches = list(cell.synapses)
    synapses_per_branch = np.int0(np.zeros(len(branches)))
    start_idx = 0
    plot_syn_indices = []
    for branch_num in range(len(branches)):
        branch = branches[branch_num]
        
        x_start, y_start = point_list[np.int0(np.ceil(branch/2))]
        x_end, y_end = point_list[branch + 1]
         
        for syn_num in range(len(cell.synapses[branch])):
            fraction_loc = list(cell.synapses[branch])[syn_num]
            x_loc = x_start + fraction_loc*(x_end-x_start)
            y_loc = y_start + fraction_loc*(y_end-y_start)
            x_locs.append(x_loc)
            y_locs.append(y_loc)
            
            synapses_per_branch[branch_num] += len(list(cell.synapses[branch][fraction_loc]))
        unique_synapse_indices = np.unique(cell.e_locations[start_idx:start_idx + synapses_per_branch[branch_num]], return_index=True)[1] + start_idx
        start_idx += synapses_per_branch[branch_num]
        plot_syn_indices.extend(unique_synapse_indices)

    plot_syn_indices.sort()

    color_vals = [np.heaviside(cell.rho_GB[i][-1]-0.5,0.5) for i in plot_syn_indices]

    shapes = ['o' for i in range(len(plot_syn_indices))]
    count = 0
    for branch in cell.synapses:
        for loc in cell.synapses[branch]:
            for weight in cell.synapses[branch][loc]:
                if len(cell.synapses[branch][loc][weight])==0:
                     count+=1
                else:
                    if count in plot_syn_indices:
                        active_syn_index = count
                    count += 1

    cmap = cm.get_cmap('bwr', 3)
    colors = cmap(color_vals)
    for i in range(len(x_locs)):
        if colors[i].all() == 1:
            lw = 0.5
        else:
            lw = 0
        if i == active_syn_index:
            scatter(x_locs[i], y_locs[i], s = 20, color=colors[i], marker = shapes[i], cmap=cmap, zorder = 2, edgecolor = 'k', linewidth = lw)
        else:
            scatter(x_locs[i], y_locs[i], s = 20, color=colors[i], marker = shapes[i], cmap=cmap, zorder = 2, edgecolor = 'k', linewidth = lw)

    
    'inset axis'
    left, bottom, width, height = [0.7, 0.03, 0.3, 0.3]
    inset_ax = ax.inset_axes([left, bottom, width, height],)
    inset_ax.spines['top'].set_visible(False)
    inset_ax.spines['right'].set_visible(False)
    inset_ax.set_xticks([])
    inset_ax.set_xticklabels([])
    inset_ax.set_yticks([])
    inset_ax.set_yticklabels([]) 
    inset_ax.plot(npa(cell.t)[:500], npa(cell.vsyn[active_syn_index])[:500], 'k')
    inset_ax.set_ylim([-77,0])
    
    if fig_num == 1:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.01, 0.7])
        cbar = plt.colorbar(cm.ScalarMappable(cmap=cmap), cax=cbar_ax)
        plt.clim(0,1)
        cbar.set_ticks([1/6, 3/6, 5/6])   
        cbar.set_ticklabels(['depression','no change','potentiation'])
    return fig, axes


def multiple_resistance_traces_subplot(cells, resistances, syn_options, section_to_plot = 0):
    '''
    function to check voltage traces with/wo NMDA spike, and at different spine neck resistances
    '''
    num_of_resistances = len(resistances)
    num_of_syn_options = len(syn_options)
    active_loc = list(cells[0].synapses[section_to_plot])
    spikes = cells[0].synapses[list(cells[0].synapses)[0]][list(cells[0].synapses[list(cells[0].synapses)[0]])[0]][list(cells[0].synapses[list(cells[0].synapses)[0]][list(cells[0].synapses[list(cells[0].synapses)[0]])[0]])[0]]
    freq = '1 spike'
    if len(spikes) > 1:
        freq = np.round(1000/(spikes[1]-spikes[0]))
    count = 0
    fig, axes = plt.subplots(num_of_syn_options, num_of_resistances, sharey = True)#figsize = (5,4), dpi = 300)
    left, bottom, width, height = [0.7, 0.7, 0.25, 0.25]
    for resistance_num in range(num_of_resistances):
        for syn_option_num in range(num_of_syn_options):
            cell = cells[count]
            ax = axes[syn_option_num,resistance_num]
            if resistance_num == 0:
                ax.set_ylabel('# of active syns = ' + str(syn_options[syn_option_num]))
            if syn_option_num == 0:
                ax.set_title('resistance = ' + str(resistances[resistance_num]) + 'M*Ohm')
  
            ax.plot(cell.t, cell.soma_v, label = 'Soma', c = 'k')
            ax.plot(cell.t, cell.dend_vs[section_to_plot][round(0.3*len(cell.dend_vs[section_to_plot]))], label = 'dend 60um', c = 'g')
            ax.plot(cell.t, cell.vsyn[0], label = 'syn 60um', ls ='--', c = 'g')
            loc = list(cell.synapses[section_to_plot])[0]
            seg_num = round(loc*len(cell.dend_vs[section_to_plot]))
            v = cell.dend_vs[section_to_plot][seg_num]
            ax.plot(cell.t, v, label = 'dend '+str(int(loc*cell.dend.L))+'um', c = 'r')
            ax.plot(cell.t, cell.vsyn[1], label = 'syn '+str(int(loc*cell.dend.L))+'um', ls ='--', c = 'r')
            ax.plot(cell.t, cell.dend_vs[section_to_plot][-1], label = 'dend 200um', c = 'b')    
            ax.plot(cell.t, cell.vsyn[-1], label = 'syn 200um', ls ='--', c = 'b')
            
            'inset axis'
            
            inset_ax = ax.inset_axes([left, bottom, width, height],)
            inset_ax.spines['top'].set_visible(False)
            inset_ax.spines['right'].set_visible(False)

            argmax = np.argmax(cell.vsyn[1])
            inset_ax.plot(npa(cell.t)[argmax-250:argmax+1000], npa(v)[argmax-250:argmax+1000], c = 'r')
            inset_ax.plot(npa(cell.t)[argmax-250:argmax+1000], npa(cell.vsyn[1])[argmax-250:argmax+1000], ls ='--', c = 'r')
            syn_max = cell.vsyn[1][argmax]
            dend_max = v[argmax]
            inset_ax.set_ylim([dend_max-0.5*(syn_max-dend_max),syn_max+0.5*(syn_max-dend_max)])
            count += 1


def plot_dendrite_heatmaps_fig_5(cell, voltage_matrices, syn_voltage_traces, locations, num_of_prox_syns_i, num_of_dist_syns_i, num_of_prox_syns_f,
                                 num_of_dist_syns_f, dt = 0.1):
    
    fig, axes = plt.subplots(1,3, figsize = (20,6))
    fig.suptitle(('prox,center,dist = ',locations[0], locations[1], locations[2]))
    count = 0
    x_positions = np.arange(0,len(cell.soma_v),100)
    x_labels = np.int0(np.linspace(0,dt*len(cell.t),len(x_positions)))
    r, c = np.shape(voltage_matrices[0])
    ylabels = list(np.linspace(0,cell.dend.L,r))
    ylabels = [round(i) for i in ylabels]
    ylabels[0] = 'soma'
    exp_labels = [str(num_of_prox_syns_f)+' proximal',str(num_of_prox_syns_i)+' proximal, '+str(num_of_dist_syns_i)+' distal',str(num_of_dist_syns_f)+' distal']
    
    for col in range(3):
    
        ax = axes[col]
        im = ax.imshow(voltage_matrices[count], origin='lower', cmap=plt.get_cmap('inferno'), vmin = -77, vmax = 0)
        ax.set_xticks(x_positions[::3])
        ax.set_xticklabels(x_labels[::3])
    
        ax.set_aspect(50)
         
        ax.set_yticks(np.arange(0, r, 1)[::2])
         
        ax.set_yticklabels(ylabels[::2])
        ax.set_ylabel('Distance from soma [um]')
        ax.set_xlabel('Time [ms]')
        ax.set_title(exp_labels[count])
        count += 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    plt.colorbar(im, cax=cbar_ax)

    fig, axes = plt.subplots(1,3, figsize = (20,6), sharey = True)
    trace_labels = ['prox_syn_traces','inactive syn trace','dist_syn_traces']
    
    for col in range(3):
        ax = axes[col]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        for trace in range(3):
            ax.plot(cell.t, syn_voltage_traces[col][trace], label = exp_labels[trace])
        ax.set_title(trace_labels[col])
        ax.set_xticks([])
        ax.set_xticklabels([])
    
    ax.legend()
    
    fig.suptitle(('prox,center,dist = ',locations[0], locations[1], locations[2]))
