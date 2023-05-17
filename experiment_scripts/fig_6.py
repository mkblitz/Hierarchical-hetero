import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)

import bams_peak_ca_plus_plasticity as bpc

max_synapses = 40

'make sure to put inactive synapse in middle of syn_branches. in 6A we used 1,3,7, and in 6D we used 7,3,8'
# syn_branches = [1,3,7]#vertical. this is branch 1,3,7 starting from the soma, and left to right. So Soma branch is #0, then left to right #1,#2, next layer L to R #3,4,5,6, etc.
# syn_branches = [3,1,4]#horizontal example using other branches 
syn_branches = [7,3,8]#horizontal

theta_d_GB = 0.5 #for 6D-F set: 0.2
theta_p_GB = 1 #for 6D-F set: 0.4

print('start')
bpc.run_experiment(max_synapses, syn_branches, theta_d_GB = theta_d_GB, theta_p_GB = theta_p_GB)
print('end')















