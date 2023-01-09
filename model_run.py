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
t_initial = time.time()
from matplotlib.pyplot import cm

# exec('nrnivmodl mods')

'''
single: for a single run, but you can control number of synapses, and nuber of spikes and spike times per synapse
multiple: check one synapse, but for many different weights
multiple_distances: check multiple synapses each with multiple weights
'''

sim_time = 200 
dt = 0.1#0.1

clamp_type = None
e_syn_locations = None
e_syn_weights = None

run_type = 'dendrite_heatmap'
# cell = bas.stick_and_ball_model(dend_L = 200, dend_cm = 2,)
spines = True
# cell.create_excitatory_synapses([0.5], [1], stimTimes = [[10]],
#                                                     with_spine = spines, synapse_type = 'glu syn from Andras')
# cell.run_sim([0.5], sim_time, dt, recording_vecs = ['vsyn'])
# print('area ', cell.e_synapses[0].area_mk)
# print('V ', cell.e_synapses[0].volume_CR_mk)
# print('surface area = ', cell.e_synapses[0].get_segment().sec.L*np.pi*cell.e_synapses[0].get_segment().sec.diam)
# print('Volume = ', cell.e_synapses[0].get_segment().sec.L*np.pi* ((1/2)*(cell.e_synapses[0].get_segment().sec.diam))**2)


# cell = bams.many_sticks_model(num_of_sections = 15, soma_diam = 1000)
# N = cell.dendrites[0].nseg
# locs = np.arange(0,cell.dendrites[0].L, cell.dendrites[0].L/N)
# colors = cm.Oranges(np.linspace(0, 1, 15))
# fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()
# cell.plot_TRs()
# for i in range(15):
#     TRs, ratios, dend_to_soma_TRs, dend_to_soma_ratios = cell.two_way_soma_TRs(i)
# #     print(i, TR)
# #     ax1.plot(locs, TRs, label = 'branch ' + str(i), color = colors[i])
# #     ax2.plot(locs, dend_to_soma_TRs, label = 'branch ' + str(i), color = colors[i])
#     ax1.plot(locs, ratios, label = 'branch ' + str(i), color = colors[i])
#     ax2.plot(locs, dend_to_soma_ratios, label = 'branch ' + str(i), color = colors[i])
# ax1.legend(), ax2.legend()


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
#     e_syn_locations = [[0.1]]
    e_syn_times = [[10,]]
    e_syn_weights = np.arange(20,70,10)
#     e_syn_weights = np.insert(e_syn_weights,0,0)
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
#     e_syn_locations = [[0.1]]
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
#     syn_locations = [0.2,0.5,0.8,]
#     clamp_amps = np.linspace(0, 0.4, 3)
    current_to_plot = ['total']#,'NMDA',]
    nmda_ratio = 2.38
    varying = 'frequency'
    if varying == 'weights':
        syn_times = [[10]]
        syn_weights = np.arange(0,60,7)
    elif varying == 'frequency':
        syn_weights = [10]
#         freq_list = [1,5,10,15,20,25,30,35,40,45,50]#in spikes/100ms
#         syn_times = [np.linspace(10,100,i+1) for i in freq_list]
#         syn_times = [[syn_time[:-1] for syn_time in syn_times]]
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
#     synapses = {0:{sec_frac:none for sec_frac in sec_fracs}}
    synapses[0][active_loc] = stims
    print(synapses)
#     synapses = {0:{}}
#     for loc in active_loc:
#         synapses[0][loc] = stims
#     print(synapses)
    spines = True
#     print(len(synapses[0]))
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
#     print(stims)
    synapses = {1:{0.5:stims_prox}, 3:{0.5:[]}, 7:{0.5:stims_dist}}#1:{0.5:{weight:[30],}},}
#     synapses = {7:{0.5:stims_dist}}
#     synapses = {3:{0.5:stims_prox,},4:{0.5:stims_dist}}
#     synapses = {3:{0.5:stims_prox,},4:{0.5:stims_dist}}
#     synapses = {4:{0.5:stims_prox}}
    spines = True
    print(synapses)

elif run_type == 'multiple_resistances':
    num_of_stims = 1
#     num_of_active_syns = 14
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
    #voltages = {'soma':stick_ball.soma_v, 'dendrites':stick_ball.dend_vs}
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
        #print(weight)
        stick_ball = bas.stick_and_ball_model()
        stick_ball.create_excitatory_synapses(e_syn_locations, [weight], stimTimes = e_syn_times)
        stick_ball.run_sim(rec_locations,sim_time, dt,)
        voltages['soma'][rec_locations[0]][weight] = stick_ball.soma_v
        #print(npa(voltages['soma'][weight])[650])
        voltages['dendrites'][rec_locations[0]][weight] = stick_ball.dend_vs[0]
        syn_currents[rec_locations[0]][weight] = stick_ball.syn_currents[0]
        #print(type(syn_currents[rec_locations[0]][weight]))
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
    
#     for sec in cell.hoc_cell.all:
#         try:
#             sec.gIhbar_Ih = 0
#         except:
#             print(sec.name())
    print('Cell length = ',cell.dend.L)
    for syn_loc in e_syn_locations:
        print('syn loc = ', syn_loc)
        e_syn_location = syn_loc
        soma_voltages[syn_loc[0]] = {}
        dend_voltages[syn_loc[0]] = {}
        syn_currents[syn_loc[0]] = {}
        for weight in e_syn_weights:
            #print(weight)
            #print('weight = ',weight)
#             if two_synapses:
#                 cell.create_excitatory_synapses([e_syn_location[0],0.9], [weight,59], stimTimes = e_syn_times, NMDA_ratio = nmda_ratio)
#             else:
            
            if excitatory:
                cell.create_excitatory_synapses(e_syn_location, [weight], stimTimes = e_syn_times, NMDA_ratio = nmda_ratio)
            else:
                cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times,)
            if clamp:
                if clamp_params['clamp_type'] == 'v':
                    cell.set_voltage_clamp(loc = clamp_params['loc'], amp = clamp_params['amp'], dur = sim_time, cell_part = clamp_params['cell_part'])
                elif clamp_params['clamp_type'] == 'i':
                    cell.set_current_clamp(loc = clamp_params['loc'], amp = clamp_params['amp'], dur = sim_time, cell_part = clamp_params['cell_part'])                        
#                 cell.set_voltage_clamp(loc = e_syn_location[0], amp = clamp_params['voltage'], dur = sim_time, cell_part = clamp_params['cell_part'])
            #cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times)
            cell.run_sim([syn_loc[0],], duration, dt, None)
            soma_voltages[syn_loc[0]][weight] = utils.npa(cell.soma_v)
            dend_voltages[syn_loc[0]][weight] = utils.npa(cell.dend_vs[0])
            syn_currents[syn_loc[0]][weight] = utils.npa(cell.syn_currents[0])
                #syn_currents[(driver_loc,modulator_loc)][(driver_weight, modulator_weight)] = utils.npa(cell.syn_currents[0])
            #print(type(syn_currents[rec_locations[0]][weight]))
#             print(cell.e_synapses)
#             print('syn_loc = ',cell.e_locations)
#             print('rec_loc = ',cell.rec_locations)
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
#     cell.plot_TRs()
#     pf.plot_ratios(cell)
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
#                         print(freq)
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [freq], NMDA_ratio = nmda_ratio)
                    else:
                        cell.create_inhibitory_synapses([syn_loc], [syn_weight], stimTimes = freq)
                    
                    cell.set_current_clamp(loc = clamp_loc, amp = clamp_amp, dur = sim_time)
                    cell.run_sim([syn_loc,clamp_loc], sim_time, dt, recording_vecs = ['i','i_NMDA','i_AMPA'])
#                     print(weight_values)
                    freq_tup = tuple(freq)
#                     print(cell.e_synapses[0])
                    voltages['dend'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.dend_vs[0][0])
                    voltages['soma'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.soma_v)
#                     syn_currents['total'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.syn_currents[0])
#                     syn_currents['AMPA'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.syn_currents_AMPA[0])
#                     syn_currents['NMDA'][(syn_loc,clamp_loc)][clamp_amp][freq_tup] = utils.npa(cell.syn_currents_NMDA[0])
                    
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
    
#     cell = L5PC.L5PC_model()
#     voltages['soma'][rec_locations[0]][weight] = stick_ball.soma_v
#     voltages['dendrites'][rec_locations[0]][weight] = stick_ball.dend_vs[0]
#     syn_currents[rec_locations[0]][weight] = stick_ball.syn_currents[0]
    for l1 in e_syn_locations:
        for l2 in e_syn_locations:
#             print(l1,l2)
            dend_voltages[l1,l2] = {}
            syn_currents[l1,l2] = {}
            if l1 <= l2:
                        #print(l1,l2)
                soma_voltages[l1,l2] = {}
#     print(soma_voltages)
            for w1 in weights:
                for w2 in weights:
                    
                    cell.create_excitatory_synapses([l1,l2], [w1,w2], stimTimes = e_syn_times)
                    cell.run_sim([l1,l2], sim_time, dt,)
                    #dend_v, soma_v,_,_ = single_run(e_syn_locations = e_syn_locations, e_syn_weights = [w1, w2], rec_locations = [l1,l2],)
                    dend_voltages[(l1,l2)][(w1, w2)] = utils.npa(cell.dend_vs[0])
                    #print(cell.syn_currents)
                    syn_currents[(l1,l2)][(w1, w2)] = utils.npa(cell.syn_currents[0])
                    if l1 <= l2:
                        #print(l1,l2)
                        soma_voltages[(l1,l2)][(w1, w2)] = utils.npa(cell.soma_v)
                    cell.clear_synapses()
    #print(soma_voltages)
#     file_name = 'location 1 = '+str(l1)+', location 2 = '+str(l2)+f', weights range = {e_syn_weights[0],e_syn_weights[-1]}'.format()
    voltages = {'dendrite':dend_voltages, 'soma':soma_voltages}
#     print(syn_currents)
#     outfile = open(path_name+'/pkl files/Figure 3 from Excitation paper/'+file_name+'.pkl','wb')
#     pkl.dump(data_for_file,outfile)
#     outfile.close()
#     print(data_for_file)
#     pf.plot_voltages_multiple_distances(stick_ball, stick_ball.t, syn_currents, dend_voltages, soma_voltages, 
#                                         syn_loc = e_syn_locations,)
    syn_c_mat = pf.plot_fig_3(voltages, syn_currents, weights, dend_L = dend_L, num_of_synapse_locations = len(e_syn_locations),
                               save_fig = save_fig)
    return dend_voltages, soma_voltages, syn_currents, cell, syn_c_mat


def dendrite_heatmap(synapses, cell_model = 'bas', nmda_ratio = 2.38, excitatory = True, num_of_rec_locations = 25, with_spines = spines):
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
    rec_locations = np.linspace(0, 1, num_of_rec_locations)
#     print(cell.dendrites[0].diam, cell.dendrites[0].Ra,)
#     cell.plot_TRs()
#     pf.plot_ratios(cell)
#     print('Cell length = ',cell.dend.L)
    'creates an inactive synapse at 0.3*L um'
#     cell.create_excitatory_synapses([0.3], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn from Andras',)
    for syn_branch in list(synapses):
        for syn_loc in list(synapses[syn_branch]):
#             print(syn_loc)
            for syn_weight in list(synapses[syn_branch][syn_loc]):
#                 print(syn_weight)
    #             print(synapses[syn_loc][syn_weight])
    #             print(synapses[syn_loc][syn_weight])
                if excitatory:
    #                 print([syn_loc])
                    cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_branch][syn_loc][syn_weight]],
                                                    with_spine = spines, synapse_type = 'glu syn from Andras', )
                else:
                    cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times,)
    'creates an inactive synapse at distal end'
#     cell.create_excitatory_synapses([1], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn from Andras', )
#     print(sim_time, dt)
    cell.run_sim(rec_locations, sim_time, dt, recording_vecs = ['vsyn','ica_NMDA','ica_VDCC'])
    voltage_matrix = pf.plot_dendrite_heatmap(cell, synapses, annotate = False, dt = dt, dendrite_or_spines = 'dendrite')
#     print(len(cell.vsyn))
    return cell, voltage_matrix


def multiple_dendrite_heatmaps(num_of_sections, synapses, sim_time = 100, cell_model = 'bams', excitatory = True, num_of_rec_locations = 10, 
                               dend_L = 50, spines = spines):
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
#     cell.plot_TRs()
#     pf.plot_ratios(cell)
    print('Branch length = ',cell.dendrites[0].L)
    for syn_sec_num in list(synapses):
        for syn_loc in list(synapses[syn_sec_num]):
            for syn_weight in list(synapses[syn_sec_num][syn_loc]):
                if excitatory:
                    cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_sec_num][syn_loc][syn_weight]], dend_sec_num = syn_sec_num,
                                                    with_spine = spines)
#             else:
#                 cell.create_inhibitory_synapses(e_syn_location, [weight], stimTimes = e_syn_times,)
    cell.run_sim(rec_locations, sim_time, dt, recording_vecs = ['vsyn'])
    active_syn_protocol = str(num_of_active_prox_syns) +'_'+ str(num_of_active_dist_syns)
    voltage_matrices = pf.plot_multiple_dendrite_heatmaps(cell, annotate = True, dt = dt, sim_time = sim_time, active_syn_protocol = active_syn_protocol)
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
#             print(synapses)
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
        #     cell.dend.L = dend_L
        #     cell.plot_TRs()
        #     pf.plot_ratios(cell)
            
            'creates an inactive synapse at 0.3*L um'
            cell.create_excitatory_synapses([0.3], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn from Andras', neck_Rn = resistance)
            for syn_branch in list(synapses):
                for syn_loc in list(synapses[syn_branch]):
        #             print(syn_loc)
                    for syn_weight in list(synapses[syn_branch][syn_loc]):
        #                 print(syn_weight)
            #             print(synapses[syn_loc][syn_weight])
            #             print(synapses[syn_loc][syn_weight])
        
                        cell.create_excitatory_synapses([syn_loc], [syn_weight], stimTimes = [synapses[syn_branch][syn_loc][syn_weight]],
                                                        with_spine = spines, synapse_type = 'glu syn from Andras', neck_Rn = resistance)
        
            'creates an inactive synapse at distal end'
            cell.create_excitatory_synapses([1], [1], stimTimes = [[]], with_spine = spines, synapse_type = 'glu syn from Andras', neck_Rn = resistance)
            
            cell.run_sim(rec_locations, sim_time, dt, recording_vecs = ['vsyn'])
            cells.append(cell)
#     print(cells)
    pf.multiple_resistance_traces_subplot(cells, neck_resistances, num_of_active_syns)
#     print(cell)
    return cells

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
    'code for new figure 5'
#     voltage_matrices = []
#      
#     prox_syn_voltage_traces = []
#     mid_syn_voltage_traces = []
#     dist_syn_voltage_traces = []
#      
#     prox_syn_ica_NMDA_traces = []
#     mid_syn_ica_NMDA_traces = []
#     dist_syn_ica_NMDA_traces = []
#      
#     prox_syn_ica_VDCC_traces = []
#     mid_syn_ica_VDCC_traces = []
#     dist_syn_ica_VDCC_traces = []
#      
#     inactive_syn_indices = []
#     count = 0
#     prox_loc, mid_loc, dist_loc = 0.3,0.65,1
#     loc = [prox_loc, mid_loc, dist_loc]
#     num_of_prox_syns_i = 20
#     num_of_dist_syns_i = 20
#     num_of_prox_syns_f = num_of_prox_syns_i + num_of_dist_syns_i
#     num_of_dist_syns_f = num_of_prox_syns_i + num_of_dist_syns_i
#     inactive_syn_indices = [num_of_prox_syns_i,1,num_of_prox_syns_i,num_of_prox_syns_f,1]
#  
#     for subplot_num in range(5):
#         synapses = {0:{}}
#         if count == 0:
#             stims = {1+i:syn_times for i in range(num_of_prox_syns_i)}
#             synapses[0][prox_loc] = stims
#             synapses[0][mid_loc] = {1:[]}
#             synapses[0][dist_loc] = {1:[]}
#         if count == 1:
#             stims = {1+i:syn_times for i in range(num_of_dist_syns_i)}
#             synapses[0][prox_loc] = {1:[]}
#             synapses[0][mid_loc] = {1:[]}
#             synapses[0][dist_loc] = stims
#         if count == 3:
#             stims_a = {1+i:syn_times for i in range(num_of_prox_syns_i)}
#             stims_b = {1+i:syn_times for i in range(num_of_dist_syns_i)}
#             synapses[0][prox_loc] = stims_a
#             synapses[0][mid_loc] = {1:[]}
#             synapses[0][dist_loc] = stims_b
#         if count == 2:
#             stims = {1+i:syn_times for i in range(num_of_prox_syns_f)}
#             synapses[0][prox_loc] = stims
#             synapses[0][mid_loc] = {1:[]}
#             synapses[0][dist_loc] = {1:[]}
#         if count == 4:
#             stims = {1+i:syn_times for i in range(num_of_dist_syns_f)}
#             synapses[0][prox_loc] = {1:[]}
#             synapses[0][mid_loc] = {1:[]}
#             synapses[0][dist_loc] = stims        
#          
#         print(synapses)
#          
#         cell, vm = dendrite_heatmap(synapses,)
#         voltage_matrices.append(vm)
#                  
#         prox_syn_voltage_traces.append(cell.vsyn[0])
#         mid_syn_voltage_traces.append(cell.vsyn[inactive_syn_indices[subplot_num]])
#         dist_syn_voltage_traces.append(cell.vsyn[-1])
#          
#         prox_syn_ica_NMDA_traces.append(cell.ica_NMDA[0])
#         mid_syn_ica_NMDA_traces.append(cell.ica_NMDA[inactive_syn_indices[subplot_num]])
#         dist_syn_ica_NMDA_traces.append(cell.ica_NMDA[-1])
#          
#         prox_syn_ica_VDCC_traces.append(cell.ica_VDCC[0])
#         mid_syn_ica_VDCC_traces.append(cell.ica_VDCC[inactive_syn_indices[subplot_num]])
#         dist_syn_ica_VDCC_traces.append(cell.ica_VDCC[-1])
#          
#         count += 1
#      
#     print(inactive_syn_indices)
#      
#     fig, axes = plt.subplots(1,3, figsize = (20,6))
#     fig.suptitle(('prox,center,dist = ',prox_loc, mid_loc, dist_loc))
# #     ax[0,0] = plt.subplot(231)
#     count = 0
#     count = 2
#     x_positions = np.arange(0,len(cell.soma_v),100)
#     x_labels = np.int0(np.linspace(0,dt*len(cell.t),len(x_positions)))
#     r, c = np.shape(voltage_matrices[0])
#     ylabels = list(np.linspace(0,cell.dend.L,r))
#     ylabels = [round(i) for i in ylabels]
#     ylabels[0] = 'soma'
#     exp_labels = [str(num_of_prox_syns_i)+' proximal,',str(num_of_dist_syns_i)+' distal',str(num_of_prox_syns_f)+' proximal',
#                   str(num_of_prox_syns_i)+' proximal, '+str(num_of_dist_syns_i)+' distal',str(num_of_dist_syns_f)+' distal']
#     for row in range(2):
#         for col in range(3):
#             if count < 5:
# #                 ax = axes[row,col]
#                 ax = axes[col]
#                 im = ax.imshow(voltage_matrices[count], origin='lower', cmap=plt.get_cmap('inferno'), vmin = -77, vmax = 0)
#                 ax.set_xticks(x_positions[::3])
#                 ax.set_xticklabels(x_labels[::3])
# #                 ax.set_yticks([])
# #                 ax.set_yticklabels([])
#                 ax.set_aspect(50)
# #     
#                 ax.set_yticks(np.arange(0, r, 1)[::2])
# #     
#                 ax.set_yticklabels(ylabels[::2])
#                 ax.set_ylabel('Distance from soma [um]')
#                 ax.set_xlabel('Time [ms]')
#                 ax.set_title(exp_labels[count])
#                 count += 1
#     fig.subplots_adjust(right=0.8)
#     cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#     plt.colorbar(im, cax=cbar_ax)#fig, ax = ax[1,1],)
# #     plt.tight_layout()
# #     fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Calcium in Dendrites - Paper/Figures/Figure 5 new/heatmaps_1.png') 
#      
#     fig, axes = plt.subplots(1,3, figsize = (20,6), sharey = True)
#     trace_labels = ['prox_syn_traces','inactive syn trace','dist_syn_traces']
#     traces = [prox_syn_voltage_traces, mid_syn_voltage_traces, dist_syn_voltage_traces]
# #     for col in range(2,3):
# #     ax = axes[1,col]
#      
#     for col in range(3):
# #         ax = axes[1,col]
#         ax = axes[col]
#         ax.spines['right'].set_visible(False)
#         ax.spines['top'].set_visible(False)
#         for trace in range(2,5):
#             ax.plot(cell.t, traces[col][trace], label = exp_labels[trace])
#         ax.set_title(trace_labels[col])
#         ax.set_xticks([])
#         ax.set_xticklabels([])
# #         ax.set_yticklabels([])
#     ax.legend()
#      
#     traces_nmda = [prox_syn_ica_NMDA_traces, mid_syn_ica_NMDA_traces, dist_syn_ica_NMDA_traces]
#     traces_vdcc = [prox_syn_ica_VDCC_traces, mid_syn_ica_VDCC_traces, dist_syn_ica_VDCC_traces]
# #     plt.tight_layout()
# #     fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Calcium in Dendrites - Paper/Figures/Figure 5 new/traces_1.png') 
#      
# #     fig, axes = plt.subplots(3,3, sharey = 'row')
#     fig.suptitle(('prox,center,dist = ',prox_loc, mid_loc, dist_loc))
#     axes[1,0].set_ylabel('ica_NMDA')
#     axes[2,0].set_ylabel('ice_VDCC')
#     for col in range(3):
#         ax_v = axes[0,col]
#         ax_nmda = axes[1,col]
#         ax_vdcc = axes[2,col]
#         for trace in range(5):
#             ax_v.plot(cell.t, traces[col][trace], label = exp_labels[trace])
#             ax_nmda.plot(cell.t, -traces_nmda[col][trace], label = exp_labels[trace])
#             ax_vdcc.plot(cell.t, -traces_vdcc[col][trace], label = exp_labels[trace])
#         ax_v.set_title(trace_labels[col])
#         ax_v.legend()
#         ax_nmda.legend()
#         ax_vdcc.legend()
        
elif run_type == 'multiple_dendrite_heatmaps':
    cell, voltage_matrices = multiple_dendrite_heatmaps(num_of_sections, synapses, sim_time = sim_time, dend_L = dend_L,)
    'code for figures 6/7'
#     prox_syns = [40,20,0]
#     dist_syns = [0,20,40]
#     
#     branches_to_plot = [4,1,3]
# #     branches_to_plot = [7,3,1]
#      
#     fig, axes = plt.subplots(3,4, sharey = 'col', figsize = (18,12))
# 
#     for exp_num in range(len(prox_syns)):
#         num_of_active_prox_syns = prox_syns[exp_num]
#         num_of_active_dist_syns = dist_syns[exp_num]
#         
#         stims_prox = {1+i:syn_times for i in range(num_of_active_prox_syns)}
#         if num_of_active_prox_syns == 0:
#             stims_prox = {1:[]}
#         stims_dist = {1+i:syn_times for i in range(num_of_active_dist_syns)}
#         if num_of_active_dist_syns == 0:
#             stims_dist = {1:[]}
#     #     print(stims)
#         if branches_to_plot == [7,3,1]:
#             synapses = {1:{0.5:stims_prox}, 3:{0.5:{1:[]}}, 7:{0.5:stims_dist}}
#             trace_numbers = [-1,np.max((num_of_active_prox_syns,1)),0]
#         elif branches_to_plot == [4,1,3]:
#             synapses = {1:{0.5:{1:[]}}, 3:{0.5:stims_prox}, 4:{0.5:stims_dist}, }
#             trace_numbers = [-1,0,1]
#             
#         print(synapses)
#         cell, voltage_matrices = multiple_dendrite_heatmaps(num_of_sections, synapses, sim_time = sim_time, dend_L = dend_L,)
#         x_positions = np.linspace(0,len(cell.soma_v),3)    
#         x_labels = np.int0(np.linspace(0,sim_time,len(x_positions)))
#         for row in range(3):
#             ax = axes[row, exp_num]
#             if row == 0:
#                 ax.set_title('prox/dist # of syns: '+str(num_of_active_prox_syns)+'/'+str(num_of_active_dist_syns))
#             im = ax.imshow(voltage_matrices[branches_to_plot[row]], origin='lower', cmap=plt.get_cmap('inferno'), vmin = -77, vmax = 0)
#             ax.set_xticks(x_positions)
#             ax.set_xticklabels([])
#             if row == 2:
#                 ax.set_xticklabels(x_labels)
#             ax.set_yticks([])
#             ax.set_yticklabels([])
#             trace_ax = axes[row, 3]
#             if exp_num == 0:
#                 ax.set_ylabel('branch #'+str(branches_to_plot[row]))
#                 
#                 trace_ax.spines['top'].set_visible(False)
#                 trace_ax.spines['right'].set_visible(False)
#             
#     #         ax.set_aspect('auto')
#             ax.set_aspect(250)
#             
#             seg_num = round(0.5*len(cell.dend_vs[3]))    
#             'for spine base voltage' 
# #             v = cell.dend_vs[branches_to_plot[row]][seg_num]
#             'for spine head voltage'
#             print(trace_numbers[row])
#             v = cell.vsyn[trace_numbers[row]]
#             
#             trace_ax.plot(cell.t, v, label = str(num_of_active_prox_syns)+'_'+str(num_of_active_dist_syns))
#             if exp_num == 2:
#                 trace_ax.legend()
# #     fig.savefig('C:/Users/mkbli/Dropbox/deep_dendrite/Calcium in Dendrites - Paper/Figures/Figure 7/heatmap_plus_traces_'+str(branches_to_plot)+'.png')
#         print(cell.necks[0].diam)
# #     prox_syn_traces = []
#     inactive_syn_traces = []
#     dist_syn_traces = []
#     
#     syn_nums = [(40,0),(20,20),(0,40)]
#     for syn_num in range(len(syn_nums)):
#         num_of_active_prox_syns = syn_nums[syn_num][0]
#         num_of_active_dist_syns = syn_nums[syn_num][1]
#     
#         stims_prox = {1+i:syn_times for i in range(num_of_active_prox_syns)}
#         stims_dist = {1+i:syn_times for i in range(num_of_active_dist_syns)}
#         synapses = {1:{0.5:stims_prox}, 3:{0.5:[]}, 7:{0.5:stims_dist}}#1:{0.5:{weight:[30],}},}
# #         print(synapses)
#         cell = multiple_dendrite_heatmaps(num_of_sections, synapses, sim_time = sim_time, dend_L = dend_L,)
#         seg_num = round(0.5*len(cell.dend_vs[3]))
#         
#         v_prox = cell.dend_vs[1][seg_num]
#         prox_syn_traces.append(v_prox)
#         
#         v = cell.dend_vs[3][seg_num]
#         inactive_syn_traces.append(v)
#         
#         v_dist = cell.dend_vs[7][seg_num]
#         dist_syn_traces.append(v_dist)
#         
#     fig, axes = plt.subplots(1,3, sharey = True)# figsize = (20,6))
#     trace_labels = ['prox_syn_traces','inactive syn trace','dist_syn_traces']
#     
#     for ax_num in range(len(axes)):
#         axes[ax_num].set_title(trace_labels[ax_num])
#         
#     for trace_num in range(len(inactive_syn_traces)):
#         axes[0].plot(cell.t, prox_syn_traces[trace_num], label = str(syn_nums[trace_num]))
#         axes[1].plot(cell.t, inactive_syn_traces[trace_num], label = str(syn_nums[trace_num]))
#         axes[2].plot(cell.t, dist_syn_traces[trace_num], label = str(syn_nums[trace_num]))
#     axes[2].legend()
    
elif run_type == 'multiple_resistances':
    cells = multiple_resistance_voltage_traces(neck_resistances, num_of_active_syns)

#utils.get_fourier_transform(syn_c[0.0][9], dt)
# file_name = parameters
# parameters_for_file = {'run_type':run_type, 'nseg':nseg, 'sim_time':sim_time, 'dt':dt, 'e_syn_locations':e_syn_locations, 
#                        'rec_locations':rec_locations, 'e_syn_weights':e_syn_weights}
# data_for_file = {'voltages':voltages, 'currents':syn_c, 'parameters':parameters_for_file}
# outfile = open(path_name+'/pkl files/'+file_name+'.pkl','wb')
# pkl.dump(data_for_file,outfile)
# outfile.close()
#print('lambda_c = ',cell.get_length_constant(),'cm')
# plt.show()
t_final = time.time()
 
print('Time elapsed = '+str(t_final-t_initial)+' seconds')
# print(syn_times)
#winsound.Beep(750, 800)


# cell = bams.many_sticks_model(15, dend_L = 50,)
# cell.plot_TRs()
# cell_bas = bas.stick_and_ball_model(dend_L = 200, dend_Ra = 100, soma_diam = 1700, soma_L = 10,
#                                          dend_diam = 1, soma_g_pas = 0.0000338, dend_g_pas = 0.0000467, dend_cm = 2,
#                                          soma_e_pas = -77, dend_e_pas = -77)
# cell_hay.dend_TRs-cell_bas.dend_TRs