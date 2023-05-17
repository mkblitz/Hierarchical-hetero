import sys
import os
import numpy as np
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
sys.path.append(path_parent)

import peak_ca_plus_plasticity_for_three_locations as pcp

max_synapses = 40#maximum number of synapses per cluster, so in total there will be 2*max_synapses
syn_locs = [0.3,0.65,1]

print('start')
pcp.run_experiment(max_synapses, syn_locs)
print('end')



















