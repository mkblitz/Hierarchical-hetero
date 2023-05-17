from __future__ import division
import sys
from neuron import h
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import utils
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcyk')
from scipy.optimize import curve_fit

class cell_parent():
    def __init__(self):
        self.clear_synapses()

    def create_excitatory_synapses(self, locations, weights, stimTimes, dend_sec_num = 0, NMDA_ratio = 2.38, synapse_type = 'glu syn',
                                   with_spine = True, neck_Rn = 226, theta_d_GB = 0.5, theta_p_GB = 1):
        '''
        parameters - 
        locations/weights: of each synapse 
        stimTimes: list of lists [[syn1_t1,syn1_t2,...],[syn2_t1,...],[syn3_t1,syn3_t2,...]] 
        time_between_spikes: mean time between spikes [default: 50ms = 20Hz]
        tau1/2: rise/decay time
        stim_timing: 'continuous' or 'controlled' depending on the timing of stimulus you want
            time_between_spikes: for 'continuous' spikes, with spaces in [ms]
            stim_start/end: for 'controlled' spikes, start and end time
            num_stims: number of stimulations in 'controlled' time
        '''
        self.synapse_type = synapse_type
        dend = self.dendrites[dend_sec_num]
        for i in range(len(locations)):
            if with_spine:
                head = h.Section()
                neck = h.Section()
                 
                head.insert('pas')
                neck.insert('pas')
     
                head.connect(neck(1))
                neck.connect(dend(locations[i]))
                 
                
                head.diam = 1*0.4#konur et al 2003
                head.L = 1*0.4
                
                head.g_pas = (2*0.0000467)#1000GigaOhm, p. 286 BioPhysics of Computation
                head.Ra = 100
                neck.Ra = 1*150
                neck.L = 1*0.66#arellano et al 2007a 
                neck.diam = utils.get_diam_from_MegaOhms(neck_Rn, neck_length = neck.L, Ra = neck.Ra) #1*0.2, arellano et al 2007a, multiplied by 0.5 with Idan to increase resistance 
                neck.g_pas = (2*0.0000467)#500MegaOhm, Harnett et al 2012
                self.heads.append(head)
                self.necks.append(neck)

            if synapse_type == 'Behabadi':
                self.e_synapses.append(h.double_exp_AMPA_NMDA_Behabadi(head(1)))
                self.e_synapses[-1].NMDA_ratio = NMDA_ratio
            elif synapse_type == 'glu syn':
                if with_spine:
                    syn = h.GluSynapse_TM_MK(head(1))
                    syn.spine_head_L = head.L
                    syn.spine_head_diam = head.diam
                    self.e_synapses.append(syn)
                else:
                    self.e_synapses.append(h.GluSynapse_TM_MK(self.dendrites[-1](locations[i])))

                self.e_synapses[-1].gmax0_AMPA = 1.5
                self.e_synapses[-1].gmax_d_AMPA = 1
                self.e_synapses[-1].gmax_p_AMPA = 2.0
                self.e_synapses[-1].gca_bar_VDCC = 0.4# 20 um-2 * 3.72 pS = 0.0744nS/um^2. 20 um-2 * 20 pS = 0.4nS/um^2
            
                self.e_synapses[-1].slope_NMDA = 1*0.072#0.062#
                self.e_synapses[-1].mgo_NMDA = 1
                self.e_synapses[-1].scale_NMDA = 2.552#3.57#
                self.e_synapses[-1].tau_r_NMDA = 0.29     
                self.e_synapses[-1].tau_d_NMDA = 70 
                self.e_synapses[-1].gmax_NMDA = 1.5#1*0.55, 3.56 from Bhabadi 2012
                self.e_synapses[-1].cao_CR = 1*2
                
                self.e_synapses[-1].theta_d_GB = theta_d_GB#
                self.e_synapses[-1].theta_p_GB = theta_p_GB#

                self.e_synapses[-1].rho0_GB = (self.e_synapses[-1].gmax0_AMPA - self.e_synapses[-1].gmax_d_AMPA)/(self.e_synapses[-1].gmax_p_AMPA - self.e_synapses[-1].gmax_d_AMPA)

            for time in stimTimes[i]:
                stim = h.NetStim()                   
                stim.number = 1
                stim.start = time
                self.e_netstims.append(stim)
                self.e_netcons.append(h.NetCon(self.e_netstims[-1], self.e_synapses[-1])) 
                                   
                self.e_netcons[-1].weight[0] = weights[i]
            self.e_locations.append(locations[i])
    
    def create_inhibitory_synapses(self, locations, weights, stimTimes,):# time_between_spikes = 50, tau1 = 0.01, tau2 = 9e9):
        '''
        parameters - 
        locations/weights: of each synapse 
        time_between_spikes: mean time between spikes [default: 50ms = 20Hz]
        tau1/2: rise/decay time
        '''
        for i in range(len(locations)):
            self.i_synapses.append(h.ProbGABAA(self.dend(locations[i])))
            self.i_synapses[i].Dep = 0
            self.i_synapses[i].Fac = 0 
          
            self.i_locations.append(locations[i])
            for time in stimTimes[i]:
                stim = h.NetStim()                   
                stim.number = 1
                stim.start = time
                self.i_netstims.append(stim)
                self.i_netcons.append(h.NetCon(self.i_netstims[-1], self.i_synapses[i]))                    
                self.i_netcons[-1].weight[0] = weights[i]
            
            self.i_locations.append(locations[i])
        
    def set_voltage_clamp(self, loc, amp, dur, cell_part = 'dend'):
        if cell_part == 'soma':
            clamp = h.SEClamp(self.soma(loc))
        elif cell_part == 'dend':
            clamp = h.SEClamp(self.dend(loc))
        clamp.dur1 = dur
        clamp.amp1 = amp
        clamp.rs = 0.00000001
        self.SE_clamps.append(clamp)
        
    def set_current_clamp(self, loc, amp, dur, cell_part = 'dend'):
        '''
        set current clamp at location in soma
        parameters -
        soma_stim_loc: location on soma to inject current
        stim_delay: start of the current injection [ms]
        stim_duration: duration of stimulation [ms]
        stim_amp: amplitude of stimulation [nA]
        '''
        if cell_part == 'soma':
            clamp = h.IClamp(self.soma(loc))
        elif cell_part == 'dend':
            clamp = h.IClamp(self.dendrites[0](loc))
        clamp.delay = 2
        clamp.dur = dur
        clamp.amp = amp
        self.I_clamps.append(clamp)
            
    def clear_synapses(self,syn_type = 'all'):
        '''
        syn_type default = 'all' to clear all synapses, or choose 'e'/'i'
        '''
        if syn_type in ['all','e']: 
            self.e_locations = [] 
            self.e_synapses = []
            self.e_netstims = []
            self.e_netcons = [] 
            self.SE_clamps = []#voltage clamp
            self.I_clamps = []#current clamp
        if syn_type in ['all','i']:
            self.i_locations = []
            self.i_synapses = []
            self.i_netstims = [] 
            self.i_netcons = []
        
    def create_recording_vectors(self, recording_vecs = ['i'], rec_soma_v = True, rec_dend_v = True, rec_syn_current = True,
                                 rec_stim = False, rec_conductances = True,):
        '''
        recording vectors
        parameters - 
        soma_stim_loc: location on soma to inject current
        rec_soma_v: record voltage at soma
        rec_dend_v: record dendritic voltage in dendrite
        rec_stim: record stimulation current
        '''
        self.t = h.Vector()#time vector
        self.t.record(h._ref_t)
        if rec_soma_v:
            self.soma_v = h.Vector()
            self.soma_v.record(self.soma(0.5)._ref_v)
        if rec_dend_v:
            if self.synapse_type == 'glu syn':
                self.dend_vs = {}
            elif self.synapse_type == 'Behabadi':
                self.dend_vs = []
            for dend_num in range(len(self.dendrites)):
                self.dend_vs[dend_num] = []
                dend = self.dendrites[dend_num]
                for loc in self.rec_locations:
                    self.dend_vs[dend_num].append(h.Vector())
                    self.dend_vs[dend_num][-1].record(dend(loc)._ref_v)
        for val in recording_vecs:
            exec('self.'+val+' = []')
        for syn in self.e_synapses:
            if self.synapse_type == 'glu syn':
                for val in recording_vecs:
                    exec('self.'+val+'.append(h.Vector())')
                    exec('self.'+val+'[-1].record(syn._ref_'+val+')')

            elif self.synapse_type == 'Behabadi':
                self.dend_vs.append(h.Vector())
                self.dend_vs[-1].record(syn._ref_vsyn)
        for syn in self.i_synapses:
            self.syn_currents.append(h.Vector())
            self.syn_currents[-1].record(syn._ref_i)
        if rec_stim:
            self.stim_current = h.Vector()
            self.stim_current.record(self.stim._ref_i)
        
    def run_sim(self, rec_locations, sim_time, dt, recording_vecs = ['i'], v_init = None, ):
        '''
        running the simulation
        parameteres -
        sim_time: length of simulation [ms]
        
        dt: timestep
        v_init: initial voltage
        soma_stim_loc: location on soma to inject current
        stim_duration: duration of stimulation [ms]

        '''
        self.dt = dt
        self.rec_locations = rec_locations
        self.create_recording_vectors(recording_vecs)
        h.tstop = sim_time
        h.steps_per_ms = dt #sets the dt, both lines are necessary
        h.dt = dt
        if (self.model_name == 'ball_and_stick' or self.model_name == 'many_sticks') and v_init == None:
            h.v_init = self.soma.e_pas
        else:
            h.v_init = v_init
        h.run()        
        
    def loc_2_seg(self, loc):
        #we have 1 transfer resistance per segment. input: location between 0 and 1. output: segment number (0 to nseg)
        nSeg = self.dend.nseg
        return int(loc*nSeg)
    
    def seg_2_loc(self, seg, dend_num = 0):
        #input: segnumber between 0 and nSeg, Output: loc between 0 and 1
        return int(seg)*1/self.dendrites[dend_num].nseg
    
    def loc_2_dist(self, loc, ref_loc = None):
        #return loc * self.dend.L
        if ref_loc is None:
            ref_loc = self.soma(0.5)
        return h.distance(loc, ref_loc)
    
    def seg_2_dist(self, seg):
        return self.loc_2_dist(self.seg_2_loc(seg))
        
    def get_dend_TRs(self, dend_num = 0, freq = 0,):
        N = self.dendrites[dend_num].nseg
        TRs = np.zeros((N,N))
        ratios = np.zeros((N,N))
        pairwise_dists = np.zeros(N*N)
        TRs_by_dist = np.zeros(N*N)
        midpoints = np.zeros(N*N)
        L1s = np.zeros(N*N)
        for ind1 in range(N):
            loc1 = ind1/float(N)
            imp = self.getTR(loc1, sec = self.dendrites[dend_num], freq = freq)
            for ind2 in range(N):
                loc2 = ind2/float(N)
                TRs[ind1][ind2] = imp.transfer(loc2, sec = self.dendrites[dend_num])
                ratios[ind1][ind2] = imp.ratio(loc2, sec = self.dendrites[dend_num])                

        self.dend_TRs, self.dend_ratios = TRs, ratios #Matrix of pairwise transfer resistances
        return TRs, ratios
    
    def getTR(self, loc, sec, freq = 0):
        imp = h.Impedance()
        imp.loc(loc, sec=sec)
        imp.compute(freq,1)
        return imp
    
    def get_input_resistances(self, dend_num = 0, freq = 0):
        N = self.dendrites[dend_num].nseg
        IRs = np.zeros(N)
        ratios = np.zeros(N)
        
        for ind1 in range(N):
            loc1 = ind1/float(N)
            imp = self.getTR(loc1, sec = self.dendrites[dend_num],)
            IRs[ind1] = imp.transfer(loc1, sec = self.dendrites[dend_num])

        return IRs,
    
    def get_soma_TRs(self, dend_num = 0, freq = 0):
        N = self.dendrites[dend_num].nseg
        TRs = np.zeros(N)
        ratios = np.zeros(N)
        imp = self.getTR(0.5, sec = self.soma, freq = freq)
        for ind1 in range(N):
            loc1 = ind1/float(N)
            TRs[ind1] = imp.transfer(loc1, sec = self.dendrites[dend_num])
            ratios[ind1] = imp.ratio(loc1, sec  = self.dendrites[dend_num])
        self.soma_TRs, self.soma_ratios = TRs, ratios #Matrix of soma
        return TRs, ratios 
    
    def plot_TRs(self, dend_num = 0):
        N = self.dendrites[dend_num].nseg
        locs = np.arange(0,self.dendrites[dend_num].L, self.dendrites[dend_num].L/N)
        plt.figure('soma')
        plt.plot(locs, self.soma_TRs)
        plt.xlabel('Distance from soma [um]')
        plt.ylabel('Transfer resistance [Mega Ohms]')
        plt.title('Transfer resistance from soma')
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        fig, axes = plt.subplots(1,3,figsize = (11,5))
        extent = [0, self.dendrites[dend_num].L, 0, self.dendrites[dend_num].L]
        im0 = axes[0].imshow(self.dend_TRs, extent = extent, origin='lower', aspect = 'auto')
        axes[0].set_xlabel('L1 Distance from soma')
        axes[0].set_ylabel('L2 Distance from soma')
        axes[0].set_title('Pairwise transfer resistances')
        axes[0].set_aspect('auto')

        cbar = plt.colorbar(im0, ax = axes[0])
        cbar.ax.set_ylabel('Transfer resistance (Mega ohms)')
        
        ncurves = 5
        for i in range(ncurves):
            sliceInd = int(i*self.dendrites[dend_num].nseg/ncurves)
            slice = self.dend_TRs[sliceInd][:]

            axes[1].plot(locs, slice, label = str(int(self.seg_2_loc(sliceInd)*self.dendrites[dend_num].L)) + 'um')
            axes[2].plot(locs-int(self.seg_2_loc(sliceInd)*self.dendrites[dend_num].L), slice)

        axes[1].set_xlabel('L2 Distance from soma')
        axes[1].set_ylabel('Transfer resistance (Mega ohms)')
        axes[1].legend(title = 'L1 Distance \nfrom soma')
        axes[2].set_xlabel('L2')
        axes[2].set_ylabel('Transfer resistance (Mega ohms)')
        plt.tight_layout(h_pad = 4) 
        plt.show()
            
    def get_length_constant(self):
        '''returns in cm'''
        lambda_c = np.sqrt((1/(self.dend.g_pas*self.dend.Ra))*((self.dend.diam*0.0001)/(4)))
        return lambda_c 
               