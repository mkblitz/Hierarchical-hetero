import matplotlib.pyplot as plt
import numpy as np
import sys
from neuron import h
import ball_and_stick as bas
import ball_and_many_sticks as bams
import pickle as pkl
import plot_functions as pf
import utils
import L5PC
plt.ion()
import time
import winsound
from matplotlib.pyplot import cm

'''
single: for a single run, but you can control number of synapses, and nuber of spikes and spike times per synapse
multiple: check one synapse, but for many different weights
multiple_distances: check multiple synapses each with multiple weights
'''
t_initial = time.time()

dt = 0.1#0.1

clamp_type = None
e_syn_locations = None
e_syn_weights = None

run_type = 'multiple_dendrite_heatmaps'
spines = True

if __name__ == "__main__":
    if run_type == 'single':
        parameters = 'Jadi'
        if parameters == 'Jadi':
            plot_v = False
            e_syn_locations = [0.2, 0.4,]
            e_syn_times = [[10,100],[15],]
            rec_locations = e_syn_locations
            e_syn_weights = [0.5, 0.6,]
    elif run_type == 'multiple':
        parameters = 'multiple'
        e_syn_locations = [0.5]
        e_syn_times = [[15]]
        rec_locations = e_syn_locations
        e_syn_weights = np.arange(0.0,5,1)
    elif run_type == 'multiple_distances':
        parameters = 'multiple_distances'
        e_syn_locations = [[0.4],[0.6],[0.8]]
        e_syn_times = [[10,]]
        e_syn_weights = np.arange(20,70,10)
        cell_model = 'bas' #'l5pc','bas', 'Hay'
        nmda_ratio =2.38
        save_fig = False
        two_synapses = False
        clamp = False
        clamp_params = {'clamp_type': 'i', 'cell_part': 'dend', 'amp': 0.3, 'loc': 0.8}
        excitatory = True
    elif run_type == 'multiple_currents':
        parameters = 'multiple_currents'
        e_syn_locations = [[0.2],[0.5],[0.8]]
        e_syn_times = [[100,]]
        e_syn_weights = np.arange(0.1,61,7)
        e_syn_weights = np.insert(e_syn_weights,0,0)
        cell_model = 'bas' #'l5pc','bas', 'Hay'
        nmda_ratio = 2.38
        save_fig = False
        two_synapses = True
        clamp = True
        v_heatmap = False
        excitatory = True
        only_NMDA = True
        if clamp:
            clamp_type = 'i'#'v','i'
            if clamp_type == 'v':
                clamp_voltages = np.linspace(-75, 0, len(e_syn_weights))
            elif clamp_type == 'i':
                clamp_voltages = np.linspace(0, 0.4, len(e_syn_weights))
            modulator_values = clamp_voltages
        else:
            modulator_values = e_syn_weights
    elif run_type == 'fig_3':
        parameters = 'fig_3'
        num_of_loc = 3
        if num_of_loc == 3:
            e_syn_locations = [0.4,0.6,0.8]
            e_syn_times = [[10],[10],[10]]
        elif num_of_loc == 5:
            e_syn_locations = [0.1,0.3,0.5,0.7,0.9]
            e_syn_times = [[10],[10],[10],[10],[10]]
        e_syn_weights = np.arange(0,60,1)
        dend_L = 200
        save_fig = False
    elif run_type == 'mean_vs_duration_scatters':
        columns = 'locations'#'clamp amps','locations'
        if columns == 'locations':
            syn_locations = [0.2,0.5,0.8,]
            clamp_amps = np.linspace(0, 0.4, 3)
        elif columns == 'clamp amps':
            syn_locations = [0.3,0.4,0.5,0.6,0.7,0.8,0.9]
            clamp_amps = np.linspace(0, 0.4, 3)
        current_to_plot = ['total']#,'NMDA',]
        nmda_ratio = 2.38
        varying = 'frequency'
        if varying == 'weights':
            syn_times = [[10]]
            syn_weights = np.arange(0,60,7)
        elif varying == 'frequency':
            syn_weights = [10]
            dur = 5
            syn_times = [[np.arange(10,i,dur) for i in range(10+dur,75,dur)]]
    elif run_type == 'dendrite_heatmap':
        num_of_stims = 1
        freq = 200
        syn_step = 1000/freq
        start_1 = 10
        syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)
        
        active_loc = 0.5
        num_of_active_syns = 18
        stims = {1+i:syn_times for i in range(num_of_active_syns)}
        sec_fracs = np.round(np.arange(0.05,1.01,0.05),3)
        none = {1:[]}
        synapses = {0:{}}
        synapses[0][active_loc] = stims
        print(synapses)
        spines = True
    elif run_type == 'multiple_dendrite_heatmaps':
        num_of_sections = 15#1,3,7,15,31,63,127
        weight = 45
        dend_L = 50
        '{syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}'
        num_of_stims = 1
        freq = 200
        syn_step = 1000/freq
        start_1 = 10
        syn_times = np.arange(start_1,start_1+num_of_stims*syn_step,syn_step)
        
        num_of_active_prox_syns = 0
        num_of_active_dist_syns = 40
    
        stims_prox = {1+i:syn_times for i in range(num_of_active_prox_syns)}
        stims_dist = {1+i:syn_times for i in range(num_of_active_dist_syns)}
        synapses = {4:{0.5:stims_dist}}
        spines = True
        print(synapses)
    
    elif run_type == 'multiple_resistances':
        num_of_stims = 1
        freq = 200
        syn_step = 1000/freq
        start_1 = 10
        syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)
        
        active_loc = 0.5
        
        spines = True
        neck_resistances = [50,500,1000]
        num_of_active_syns = [7,15]
    
def single_run(e_syn_locations = e_syn_locations, e_syn_weights = e_syn_weights, rec_locations = e_syn_locations, 
               plot_voltages = False):
    '''
    runs simulation a single time, with parameters from above
    dend_voltages = {(location, weight):[voltage_vec],(location, weight):[voltage_vec],,,}
    soma_voltages = {((locations),(weights)):[soma_voltage_vec]}
    '''
    stick_ball = bas.stick_and_ball_model(dend_L = 200,)
    stick_ball.create_excitatory_synapses(e_syn_locations, e_syn_weights,e_syn_times)
    stick_ball.run_sim(rec_locations, sim_time, dt,)
    dend_voltages = {}
    soma_voltages = {}
    syn_currents = {}
    for i in range(len(e_syn_locations)):
        dend_voltages[(e_syn_locations[i],e_syn_weights[i])] = stick_ball.dend_vs[i]
    soma_voltages[tuple(e_syn_locations),tuple(e_syn_weights)] = stick_ball.soma_v
    for i in range(len(rec_locations)):
        syn_currents[rec_locations[i]] = {e_syn_weights[i]:stick_ball.syn_currents[i]}
    if plot_voltages:
        pf.plot_voltages(stick_ball.t, syn_currents, dend_voltages, soma_voltages[list(soma_voltages.keys())[0]], 
                         run_type, weights = e_syn_weights,)
    return dend_voltages, soma_voltages, syn_currents , stick_ball


def multiple_run():
    '''
    runs simulation multiple times, getting the voltage at a single synapse with different weights
    voltages = {'soma':{syn_location:{weight:voltage_vec}}, 'dendrites': {location: {weight:voltage_vec}}}
    syn_currents = {loc: {weight: current_vec}}
    '''
    
    voltages = {'soma':{rec_locations[0]:{}}, 'dendrites':{rec_locations[0]:{}}}
    syn_currents = {rec_locations[0]:{}}
    
    for weight in e_syn_weights:
        stick_ball = bas.stick_and_ball_model()
        stick_ball.create_excitatory_synapses(e_syn_locations, [weight], stimTimes = e_syn_times)
        stick_ball.run_sim(rec_locations,sim_time, dt,)
        voltages['soma'][rec_locations[0]][weight] = stick_ball.soma_v
        voltages['dendrites'][rec_locations[0]][weight] = stick_ball.dend_vs[0]
        syn_currents[rec_locations[0]][weight] = stick_ball.syn_currents[0]
        stick_ball.clear_synapses()
    #meta_parameters = stick_ball.get_parameters()
    pf.plot_voltages(stick_ball.t, syn_currents, voltages['dendrites'],voltages['soma'][rec_locations[0]], run_type, 
                     syn_loc = e_syn_locations[0], dend_L =  stick_ball.dend.L)
#     pf.plot_peak_responses(voltages, e_syn_weights, what_to_plot = 'max')
    return voltages, syn_currents ,stick_ball.t #stick_ball


def multiple_distances_run(e_syn_times, cell_model = 'l5pc', nmda_ratio = 2.38, save_fig = False, dt = 0.025, duration = 0.2, two_synapses = False,
                           excitatory = True, clamp = False, clamp_params = {'cell_part': 'dend', 'voltage':0, 'loc':0.5}):
    '''
    runs simulation multiple times, getting the voltage at a single synapse in different locations and with different weights
    voltages = {'soma':{'syn_location':{weight:voltage_vec}}, 'dendrites': {location: {weight:voltage_vec}}}
    syn_currents = {loc: {weight: current_vec}}
    '''
    if two_synapses:
        e_syn_times = [[10,],[10]]
    fig_name = 'cell_type = '+cell_model+', NMDA_ratio = '+str(nmda_ratio)
    soma_voltages = {}
    dend_voltages = {}
    syn_currents = {}
    TR_to_soma = {}
    e_pas = -75
    v_init = -50
    if cell_model == 'l5pc':
        cell = L5PC.L5PC_model()
    elif cell_model == 'Hay':
        cell = L5PC.L5PC_model(Hay = True)
    elif cell_model == 'bas':
        cell = bas.stick_and_ball_model(dend_L = 200, dend_Ra = 100, soma_g_pas=0.03, soma_e_pas = e_pas, 
                                        dend_e_pas = e_pas, dend_diam = 1)
    
    cell.plot_TRs()

    print('Cell length = ',cell.dend.L)
    for syn_loc in e_syn_locations:
        print('syn loc = ', syn_loc)
        e_syn_location = syn_loc
        soma_voltages[syn_loc[0]] = {}
        dend_voltages[syn_loc[0]] = {}
        syn_currents[syn_loc[0]] = {}
        for weight in e_syn_weights:            
            if excitatory:
                cell.create_excitatory_synapses(e_syn_location, [weight], stimTimes = e_syn_times, NMDA_ratio = nmda_ratio)
            else:
                cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times,)
            if clamp:
                if clamp_params['clamp_type'] == 'v':
                    cell.set_voltage_clamp(loc = clamp_params['loc'], amp = clamp_params['amp'], dur = sim_time, cell_part = clamp_params['cell_part'])
                elif clamp_params['clamp_type'] == 'i':
                    cell.set_current_clamp(loc = clamp_params['loc'], amp = clamp_params['amp'], dur = sim_time, cell_part = clamp_params['cell_part'])                        
            #cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times)
            cell.run_sim([syn_loc[0],], duration, dt, None)
            soma_voltages[syn_loc[0]][weight] = utils.npa(cell.soma_v)
            dend_voltages[syn_loc[0]][weight] = utils.npa(cell.dend_vs[0])
            syn_currents[syn_loc[0]][weight] = utils.npa(cell.syn_currents[0])

            cell.clear_synapses()
    #meta_parameters = cell.get_parameters()
    plot = False
    if plot:
        pf.plot_voltages_multiple_distances(cell, cell.t, syn_currents, dend_voltages, soma_voltages, syn_loc = e_syn_locations,
                                            save_fig = save_fig, fig_name = fig_name)
        scaling_factors_V = pf.plot_peak_responses(cell, dend_voltages, soma_voltages, e_syn_weights, 'voltage', peak_type = 'mean',
                                                    save_fig = save_fig, fig_name = fig_name)
        scaling_factors_I = pf.plot_peak_responses(cell, syn_currents, soma_voltages, e_syn_weights, 'current',  peak_type = 'mean',
                                                    save_fig = save_fig, fig_name = fig_name)
#     pf.plot_proximal_vs_distal(dend_voltages, soma_voltages, cell.t)
    
    return soma_voltages, dend_voltages, syn_currents ,cell.t, syn_loc, cell


def mean_vs_duration_scatters(syn_locs, clamp_amps, weight_values, syn_times, varying, excitatory = True, cell_model = 'bas', nmda_ratio = 2.38, 
                              current_to_plot = 'total', columns = 'locations'):
    '''
    creates NxM subplots (N=num of current clamp amps, M=num of syn locations), where each subplot has syn_loc_m = clamp_loc_m, clamp_amp_val_n,
    scatter of spike duration as a function of mean voltage over duration, with colors representing weight values
    
    votages = {'dend':{(syn_loc, clamp_loc): {clamp_amp: {weight_1: {vec1},...}}}}
    if with frequency: votages = {'dend':{(syn_loc, clamp_loc): {clamp_amp: {{frequncy1: {vec1},...}}}}}
    '''
    clamp_locs = syn_locs
    e_pas = -75
    voltages = {'dend':{},'soma':{}}
    syn_currents = {'total':{},'AMPA':{}, 'NMDA':{}}
    if cell_model == 'l5pc':
        cell = L5PC.L5PC_model()
    elif cell_model == 'Hay':
        cell = L5PC.L5PC_model(Hay = True)
    elif cell_model == 'bas':
        cell = bas.stick_and_ball_model(dend_L = 200, dend_Ra = 100, soma_g_pas=0.03, soma_e_pas = e_pas, 
                                        dend_e_pas = e_pas, dend_diam = 1)
    elif cell_model == 'bams':
        cell = bams.many_sticks_model(num_of_sections = 1, dend_L = 200, dend_Ra = 100, soma_g_pas=0.03, soma_e_pas = e_pas, 
                                        dend_e_pas = e_pas, dend_diam = 1,)
    print('Cell length = ',cell.dend.L)
    for syn_loc in syn_locations:
        clamp_loc = syn_loc
        for v in list(voltages):
            voltages[v][syn_loc,clamp_loc] = {}
        for current in list(syn_currents):
            syn_currents[current][syn_loc,clamp_loc] = {}
        for clamp_amp in clamp_amps:
            for current in list(syn_currents):
                syn_currents[current][syn_loc,clamp_loc][clamp_amp] = {}
            for v in list(voltages):
                voltages[v][syn_loc,clamp_loc][clamp_amp] = {}
            if varying == 'weights':
                for syn_weight in weight_values:
                    if excitatory:
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = syn_times, NMDA_ratio = nmda_ratio)
                    else:
                        cell.create_inhibitory_synapses([syn_loc], [syn_weight], stimTimes = syn_times)
                    cell.set_current_clamp(loc = clamp_loc, amp = clamp_amp, dur = sim_time)
                    cell.run_sim([syn_loc,clamp_loc], sim_time, dt,)
                    voltages['dend'][(syn_loc,clamp_loc)][clamp_amp][syn_weight] = utils.npa(cell.dend_vs[0])
                    voltages['soma'][(syn_loc,clamp_loc)][clamp_amp][syn_weight] = utils.npa(cell.soma_v)
                    syn_currents['total'][(syn_loc,clamp_loc)][clamp_amp][syn_weight] = utils.npa(cell.syn_currents[0])
                    syn_currents['AMPA'][(syn_loc,clamp_loc)][clamp_amp][syn_weight] = utils.npa(cell.syn_currents_AMPA[0])
                    syn_currents['NMDA'][(syn_loc,clamp_loc)][clamp_amp][syn_weight] = utils.npa(cell.syn_currents_NMDA[0])
                    cell.clear_synapses()
            elif varying == 'frequency':
                syn_weight = weight_values[0]
                for freq in syn_times[0]:
                    if excitatory:
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [freq], NMDA_ratio = nmda_ratio)
                    else:
                        cell.create_inhibitory_synapses([syn_loc], [syn_weight], stimTimes = freq)
                    
                    cell.set_current_clamp(loc = clamp_loc, amp = clamp_amp, dur = sim_time)
                    cell.run_sim([syn_loc,clamp_loc], sim_time, dt, recording_vecs = ['i','i_NMDA','i_AMPA'])
                    freq_tup = tuple(freq)
                    voltages['dend'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.dend_vs[0][0])
                    voltages['soma'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.soma_v)
                    
                    cell.clear_synapses()
    print(voltages)
    print(syn_currents)
#     pf.plot_mean_vs_duration_scatters(syn_currents, voltages, syn_times, varying, currents_to_plot = current_to_plot, stim_start = 111, columns = columns)
#     pf.plot_mean_vs_weights_scatters(syn_currents, voltages, syn_times, varying, currents_to_plot = current_to_plot, columns = columns)
    
    return cell 


def multiple_currents_run(e_syn_times, modulator_values, cell_model = 'l5pc', nmda_ratio = 2.38, save_fig = False, dt = 0.025, duration = 0.2, 
                          two_synapses = False, clamp = False, clamp_type = None, excitatory = True, only_NMDA = True):
    '''
    runs simulation multiple times, getting the voltage at a single synapse in different locations and with different weights
    voltages = {'soma':{'syn_location':{weight:voltage_vec}}, 'dendrites': {location: {weight:voltage_vec}}}
    syn_currents = {loc: {weight: current_vec}}
    '''
    if two_synapses:
        e_syn_times = [[10,],[10]]
    fig_name = 'cell_type = '+cell_model+', NMDA_ratio = '+str(nmda_ratio)
    soma_voltages = {}
    dend_voltages = {}
    syn_currents = {}
    TR_to_soma = {}
    e_pas = -75
    if cell_model == 'l5pc':
        cell = L5PC.L5PC_model()
    elif cell_model == 'Hay':
        cell = L5PC.L5PC_model(Hay = True)
    elif cell_model == 'bas':
        cell = bas.stick_and_ball_model(dend_L = 200, dend_Ra = 100, soma_g_pas=0.03, soma_e_pas = e_pas, 
                                        dend_e_pas = e_pas, dend_diam = 1)
    print('Cell length = ',cell.dend.L)
    for nmda_ratio in [1]:
        print('nmda_ratio = ',nmda_ratio)
        for driver_loc in e_syn_locations:
            print('driver loc = ',driver_loc)
            driver_loc = driver_loc[0]
            for modulator_loc in e_syn_locations:
                modulator_loc = modulator_loc[0]
                dend_voltages[driver_loc,modulator_loc] = {}
                if driver_loc <= modulator_loc:
                    soma_voltages[driver_loc,modulator_loc] = {}
                syn_currents[driver_loc,modulator_loc] = {}
                for num in range(len(modulator_values)):
                    if two_synapses:
                        if clamp:
                            if v_heatmap:
                                driver_weight = e_syn_weights[num]
                                for modulator_weight in modulator_values:
                                    if excitatory:
                                        cell.create_excitatory_synapses([driver_loc], [driver_weight], stimTimes = e_syn_times, NMDA_ratio = nmda_ratio)
                                    else:
                                        cell.create_inhibitory_synapses([driver_loc], [driver_weight], stimTimes = e_syn_times)
                                    if clamp_type == 'v':
                                        cell.set_voltage_clamp(loc = modulator_loc, amp = modulator_weight, dur = sim_time)
                                    elif clamp_type == 'i':
                                        cell.set_current_clamp(loc = modulator_loc, amp = modulator_weight, dur = sim_time)
                                    cell.run_sim([driver_loc,modulator_loc], sim_time, dt, only_NMDA)
                                    dend_voltages[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.dend_vs[0])
                                    syn_currents[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.syn_currents[0])
                                    if driver_loc <= modulator_loc:
                                        soma_voltages[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.soma_v)
                                    cell.clear_synapses()
                            else:#not v_heatmap
                                driver_weight = 50
                                modulator_value = modulator_values[num]
                                cell.create_excitatory_synapses([driver_loc], [driver_weight], stimTimes = e_syn_times, 
                                                            NMDA_ratio = nmda_ratio)
                                if clamp_type == 'v':
                                    cell.set_voltage_clamp(loc = modulator_loc, amp = modulator_value, dur = sim_time)
                                elif clamp_type == 'i':
                                    cell.set_current_clamp(loc = modulator_loc, amp = modulator_value, dur = sim_time)
                        else:#not clamp
                            driver_weight = e_syn_weights[num]
                            for modulator_weight in modulator_values:
                                    if excitatory:
                                        cell.create_excitatory_synapses([driver_loc,modulator_loc], [driver_weight,modulator_weight], stimTimes = e_syn_times,
                                                                        NMDA_ratio = nmda_ratio)
                                    else:
                                        cell.create_inhibitory_synapses([driver_loc,modulator_loc], [driver_weight,modulator_weight], stimTimes = e_syn_times)
                                    cell.run_sim([driver_loc,modulator_loc], sim_time, dt, only_NMDA)
                                    dend_voltages[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.dend_vs[0])
                                    syn_currents[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.syn_currents[0])
                                    if driver_loc <= modulator_loc:
                                        soma_voltages[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.soma_v)
                                    cell.clear_synapses()
                    else:#not two_synapses
                        weight = e_syn_weights[num]
                        cell.create_excitatory_synapses([driver_loc], [weight], stimTimes = e_syn_times, NMDA_ratio = nmda_ratio)
                    if not v_heatmap:
                        cell.run_sim([driver_loc,modulator_loc], duration, dt, only_NMDA)
                        if two_synapses:
                            syn_currents[driver_loc,modulator_loc][driver_weight,modulator_value] = utils.npa(cell.syn_currents[0])
                        else:
                            syn_currents[driver_loc][weight] = utils.npa(cell.syn_currents[0])
                        cell.clear_synapses()
        #meta_parameters = cell.get_parameters()
        plot = False
        if plot:
            print('a')
            pf.plot_voltages_multiple_distances(cell, cell.t, syn_currents, dend_voltages, soma_voltages, syn_loc = e_syn_locations,
                                                save_fig = save_fig, fig_name = fig_name)
            scaling_factors_V = pf.plot_peak_responses(cell, dend_voltages, soma_voltages, e_syn_weights, 'voltage', peak_type = 'mean',
                                                        save_fig = save_fig, fig_name = fig_name)
            scaling_factors_I = pf.plot_peak_responses(cell, syn_currents, soma_voltages, e_syn_weights, 'current',  peak_type = 'mean',
                                                        save_fig = save_fig, fig_name = fig_name)
        voltages = {'dendrite':dend_voltages, 'soma':soma_voltages}
        if v_heatmap:
            _ = pf.plot_fig_3(voltages, syn_currents, modulator_values, num_of_synapse_locations = len(e_syn_locations), NMDA_ratio = nmda_ratio)
        else:
            pf.plot_currents(syn_currents, num_of_locations = len(e_syn_locations))
    return soma_voltages, dend_voltages, syn_currents ,cell.t, cell, voltages
    
    
def data_for_fig_3(weights, dend_L = 200, save_fig = False):

    '''
    gets data for figure 3 from excitation paper, for synapses
    dend_voltages_l1 = {(loc1, loc2):{(w1, w2):[voltage_vec], (w1, w2):[,],,,}}
    dend_voltages_l2 = {(loc2, loc1):{(w2, w1):[voltage_vec], (w2, w1):[,],,,}}
    soma_voltages = {(loc1,loc2):{(w1,w2):[soma_voltage_vec], (w1, w2):[,],,,}}
    data for file = ({(loc1,loc2):{(w1,w2): voltage_vec, (w1,w2):...},(loc2,loc1):{(w2,w1):voltage_vec,...}, soma})
    '''
    dend_voltages = {}
    soma_voltages = {}
    syn_currents = {}
    
    cell = bas.stick_and_ball_model(dend_L = dend_L,)
    
    for l1 in e_syn_locations:
        for l2 in e_syn_locations:
            dend_voltages[l1,l2] = {}
            syn_currents[l1,l2] = {}
            if l1 <= l2:
                soma_voltages[l1,l2] = {}
            for w1 in weights:
                for w2 in weights:
                    
                    cell.create_excitatory_synapses([l1,l2], [w1,w2], stimTimes = e_syn_times)
                    cell.run_sim([l1,l2], sim_time, dt,)
                    #dend_v, soma_v,_,_ = single_run(e_syn_locations = e_syn_locations, e_syn_weights = [w1, w2], rec_locations = [l1,l2],)
                    dend_voltages[(l1,l2)][(w1, w2)] = utils.npa(cell.dend_vs[0])
                    syn_currents[(l1,l2)][(w1, w2)] = utils.npa(cell.syn_currents[0])
                    if l1 <= l2:
                        soma_voltages[(l1,l2)][(w1, w2)] = utils.npa(cell.soma_v)
                    cell.clear_synapses()

    voltages = {'dendrite':dend_voltages, 'soma':soma_voltages}
#     pf.plot_voltages_multiple_distances(stick_ball, stick_ball.t, syn_currents, dend_voltages, soma_voltages, 
#                                         syn_loc = e_syn_locations,)
    syn_c_mat = pf.plot_fig_3(voltages, syn_currents, weights, dend_L = dend_L, num_of_synapse_locations = len(e_syn_locations),
                               save_fig = save_fig)
    return dend_voltages, soma_voltages, syn_currents, cell, syn_c_mat


def dendrite_heatmap(synapses, cell_model = 'bas', nmda_ratio = 2.38, excitatory = True, num_of_rec_locations = 25, with_spines = spines, plot = True, 
                     active_loc = 0.5, scale_ih = False, sim_time = 100):
    '''
    creates a heatmap of voltage as a function of time and location on the dendrite. 
    syn_locations: [syn1, syn2,...]
    syn_times, weights: [[syn1_1,syn1_2,...],[syn2_1,syn2_2,...]]
    '''
    dend_L = 200
    
    if cell_model == 'l5pc':
        cell = L5PC.L5PC_model()
    elif cell_model == 'Hay':
        cell = L5PC.L5PC_model(Hay = True)
    elif cell_model == 'bas':
        if with_spines:#l_list = [cell.hoc_cell.apic[i].L for i in range(len(cell.hoc_cell.apic))]
            dend_cm = 2
        cell = bas.stick_and_ball_model(dend_L = dend_L,
                                         dend_cm = dend_cm)
        if scale_ih:
            cell.dend.insert('Ih')
            constant_ih = False
            for seg in cell.dend:
#                 seg.insert('Ih')
                if constant_ih:
                    seg.gIhbar_Ih = 0.0002
                else:
                    distance = h.distance(seg, cell.soma(1))
#                     seg.gIhbar_Ih = 1e-7*(-0.8696+2.087*np.exp(distance/323))
                    seg.gIhbar_Ih = 1e-4*(-2+4.28*np.exp(distance/(0.135*323)))
                
    rec_locations = np.linspace(0, 1, num_of_rec_locations)
    'creates an inactive synapse at 0.3*L um'
#     cell.create_excitatory_synapses([0.3], [1], stimTimes = [[]], with_spine = spines, glu syn = 'glu syn',)
    for syn_branch in list(synapses):
        for syn_loc in list(synapses[syn_branch]):
            for syn_weight in list(synapses[syn_branch][syn_loc]):
                if excitatory:
                    cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_branch][syn_loc][syn_weight]],
                                                    with_spine = spines, synapse_type = 'glu syn', )
                else:
                    cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times,)
    'creates an inactive synapse at distal end'
#     cell.create_excitatory_synapses([1], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn', )
    cell.run_sim(rec_locations, sim_time, dt, recording_vecs = ['vsyn','ica_NMDA','ica_VDCC'])
    if cell_model == 'bas':
        voltage_matrix = pf.plot_dendrite_heatmap(cell, synapses, annotate = False, dt = dt, dendrite_or_spines = 'dendrite', plot = plot)
    elif cell_model == 'Hay':
        voltage_matrix = pf.plot_dendrite_heatmap_hay(cell, synapses, annotate = False, dt = dt, dendrite_or_spines = 'dendrite', plot = plot, 
                                                      active_loc = active_loc)
    return cell, voltage_matrix


def multiple_dendrite_heatmaps(num_of_sections, synapses, num_of_active_prox_syns, num_of_active_dist_syns,
                               sim_time = 100, cell_model = 'bams', excitatory = True, num_of_rec_locations = 10, 
                               dend_L = 50, spines = spines, voltage_trace = False, branch_num = 4):
    '''
    creates a heatmap of voltage as a function of time and location on the dendrite. 
    synapses: {syn1_sec_num:{syn1_location:{weight:time}},syn2_sec_num:{syn2_location:{weight:time}},...}
    syn_locations: [syn1, syn2,...]
    syn_times, weights: [[syn1_1,syn1_2,...],[syn2_1,syn2_2,...]]
    '''
    dend_sec_L = dend_L
    if cell_model == 'bas':
        cell = bas.stick_and_ball_model(dend_L = dend_L, dend_Ra = 100,
                                         dend_diam = 1)
    elif cell_model == 'bams':
        cell = bams.many_sticks_model(num_of_sections, dend_L = dend_sec_L,)
    cell.synapses = synapses
    rec_locations = np.linspace(0, 1, num_of_rec_locations)
    print('Branch length = ',cell.dendrites[0].L)
    for syn_sec_num in list(synapses):
        for syn_loc in list(synapses[syn_sec_num]):
            for syn_weight in list(synapses[syn_sec_num][syn_loc]):
                if excitatory:
                    cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_sec_num][syn_loc][syn_weight]], dend_sec_num = syn_sec_num,
                                                    with_spine = spines)
    cell.run_sim(rec_locations, sim_time, dt, recording_vecs = ['vsyn'])
    active_syn_protocol = str(num_of_active_prox_syns) +'_'+ str(num_of_active_dist_syns)
    voltage_matrices = pf.plot_multiple_dendrite_heatmaps(cell, annotate = True, dt = dt, sim_time = sim_time, active_syn_protocol = active_syn_protocol,
                                                          voltage_trace = voltage_trace, branch_num = branch_num)
    return cell, voltage_matrices


def multiple_resistance_voltage_traces(neck_resistances, num_of_active_syns, cell_model = 'bas', num_of_rec_locations = 25, with_spines = spines, ):
    dend_L = 200
    print('Cell length = ',dend_L)
    soma_diam = 10
    cells = []
    for resistance in neck_resistances:
        for num_of_syns in num_of_active_syns:
            stims = {1+i:syn_times for i in range(num_of_syns)}
            synapses = {0:{active_loc:stims}}
            if cell_model == 'l5pc':
                cell = L5PC.L5PC_model()
            elif cell_model == 'Hay':
                cell = L5PC.L5PC_model(Hay = True)
            elif cell_model == 'bas':
                if with_spines:#l_list = [cell.hoc_cell.apic[i].L for i in range(len(cell.hoc_cell.apic))]
                    dend_cm = 2
                cell = bas.stick_and_ball_model(dend_L = dend_L,
                                                 dend_cm = dend_cm)
            cell.synapses = synapses
            rec_locations = np.linspace(0, 1, num_of_rec_locations)
            'creates an inactive synapse at 0.3*L um'
            cell.create_excitatory_synapses([0.3], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn', neck_Rn = resistance)
            for syn_branch in list(synapses):
                for syn_loc in list(synapses[syn_branch]):
                    for syn_weight in list(synapses[syn_branch][syn_loc]):
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_branch][syn_loc][syn_weight]],
                                                        with_spine = spines, synapse_type = 'glu syn', neck_Rn = resistance)
        
            'creates an inactive synapse at distal end'
            cell.create_excitatory_synapses([1], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn', neck_Rn = resistance)
            
            cell.run_sim(rec_locations, sim_time, dt, recording_vecs = ['vsyn'])
            cells.append(cell)
    pf.multiple_resistance_traces_subplot(cells, neck_resistances, num_of_active_syns)
    return cells

if __name__ == "__main__":
    if run_type == 'single':
        dend_voltages, soma_voltages, syn_c, sb = single_run(plot_voltages = plot_v)
    elif run_type == 'multiple':
        voltages, syn_c,t = multiple_run()
    elif run_type == 'multiple_distances':
        soma_voltages, dend_voltages, syn_c,t,sl, cell = multiple_distances_run(e_syn_times, cell_model, nmda_ratio, save_fig, dt, sim_time,
                                                                                two_synapses, excitatory,clamp, clamp_params)
    elif run_type == 'multiple_currents':
        soma_voltages, dend_voltages, syn_c,t, cell, voltages = multiple_currents_run(e_syn_times, modulator_values, cell_model, nmda_ratio, save_fig, dt, sim_time,
                                                                                two_synapses, clamp, clamp_type, excitatory, only_NMDA)
    elif run_type == 'fig_3':
        dend_v, soma_v, syn_c, cell, syn_c_mat = data_for_fig_3(e_syn_weights, dend_L, save_fig)
    elif run_type == 'mean_vs_duration_scatters':
        cell = mean_vs_duration_scatters(syn_locations, clamp_amps, syn_weights, syn_times, varying, current_to_plot = current_to_plot, nmda_ratio = nmda_ratio,
                                         columns = columns)
    elif run_type == 'dendrite_heatmap':
        cell, vm = dendrite_heatmap(synapses,)
            
    elif run_type == 'multiple_dendrite_heatmaps':
        cell, voltage_matrices = multiple_dendrite_heatmaps(num_of_sections, synapses, num_of_active_prox_syns, num_of_active_dist_syns,
                                                            sim_time = sim_time, dend_L = dend_L,)
    elif run_type == 'multiple_resistances':
        cells = multiple_resistance_voltage_traces(neck_resistances, num_of_active_syns)
