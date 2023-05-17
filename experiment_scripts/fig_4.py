import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)
import colored_branches as cb

num_of_stims = 1#number of stimulations at every active synapse
freq = 200
syn_step = 1000/freq
start_1 = 10#time of first stimulation after start of simulation    
syn_times = np.arange(start_1,start_1+num_of_stims*syn_step-1,syn_step)

theta_d_GB = 0.5#depression threshold
theta_p_GB = 1#potentiation threshold

branch_L = 50#length of each dendritic segment 
num_of_active_syns = [20,30,40,50]#number of active spines in each column, going left to right 
act_vec = [0,1,3,7,]#on which segments the activation should happen per row, going top to bottom 
num_of_branches = 15#how many sections the dendrite will have. for a symmetrical tree (1 parent 2 daughters), options are 1,3,7,15,31,63,127...

cell = cb.colored_branches(num_of_active_syns, act_vec, syn_times, num_of_branches, branch_L, theta_d_GB, theta_p_GB)