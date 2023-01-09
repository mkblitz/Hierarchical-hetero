import numpy as np
from neuron import h
from scipy import fft
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

def get_dist(frac1, frac2 = 0, dend_length = 500):
	#returns the distance between two dendritic locations, 
	#specified by fraction (between 0 and 1) along the dendrite
	return frac1*dend_length - frac2*dend_length
	
	
def npa(vec):
    """
    takes a neuron vector and returns a numpy array
    """
    return vec.as_numpy().copy()


def loc_2_TR_ind(cell, loc):
	nSeg = cell.dend.nseg
	return int(loc*nSeg)
   
   
def get_trace_avgs(trace, get_what = 'mean', trace_type = 'voltage'):
	#gets the mean, min max of a dict of traces
# 	#trace is {syn_location: {weight:voltage_vec}}}
	flip = 1
	if trace_type == 'current':
		flip = -1
	locs = list(trace)
	if get_what == 'mean':
		peaks = {loc:[np.mean(flip*trace[loc][weight]) for weight in list(trace[loc])] for loc in locs}
	elif get_what == 'max': 
		peaks = {loc:[np.max(flip*trace[loc][weight]) for weight in list(trace[loc])] for loc in locs}
# 	elif get_what == 'min': 
# 		peaks = {loc:[np.min(trace[loc][weight]) for weight in list(trace[loc])] for loc in locs}
	return peaks
	
# def computeTR(cell, refLoc):
#     imp = h.Impedance()
#     imp.loc(0.5, sec=cell.trees[refLoc[0]][refLoc[1]])
#     imp.compute(0, 1)
#     return imp
#     
# def getTR(cell, key, refLoc):
#     imp = computeTR(cell, refLoc)
#     tree = cell.trees[key[0]]        
#     section = tree[key[1]]
#     segLoc = cell.getSegLoc(key)       
#     return imp.transfer(segLoc, sec=section)
   
# def getTRs(cell):
# 	N = cell.nsegs
# 	TRs = np.zeros(N,N)
# 	for loc1 in np.arange(0,1, float(1/N))
# 		imp = getTR(cell, loc1)
# 		for loc2 in np.arange(0,1, float(1/N))
# 			TR[loc1, loc2] = imp.transfer(loc2, cell.dend)
# 	return TRs
   
   
def get_scaling_factors(peaks_dend, peaks_soma, dend_trace_type = 'current', rest = -75):
	locs = list(peaks_dend)
	if dend_trace_type == 'current':
 		dend_soma_ratio = {loc:(np.array(peaks_soma[loc])-rest)/np.array(peaks_dend[loc]) for loc in locs}
	else:
		dend_soma_ratio = {loc:(np.array(peaks_soma[loc])-rest)/(np.array(peaks_dend[loc])-rest) for loc in locs}
	scaling_factors = {loc:np.nanmean(dend_soma_ratio[loc]) for loc in locs}
	return dend_soma_ratio, scaling_factors


def get_longest_section_connected_to_soma():
	for sec_num in range(len(cell.trees[1])):
	    sec = cell.trees[1][sec_num]
	    if cell.loc_2_dist(sec(0)) == 0.0:
	        print(sec_num,sec.L)

	        
def get_fourier_transform(vec, dt = 0.025):
	
	#sampling_rate = round(1/dt)
	#t = np.arange(0,duration,dt) #0 to 100 milliseconds
	#x = (t*sampling_rate)/duration
	dt /= 1000 #dt = dt/1000
	ft = fft(vec,)# norm='forward')           #fourier transform
	plt.figure()
	plt.subplot(2,1,1)
	plt.plot(vec)
	
	#plt.xlabel('time')
	plt.ylabel('current')
	
	plt.subplot(2,1,2)
	#plt.plot(x[:len(x)//2], abs(ft[:len(vec)//2]))
	n = vec.size
	fft_freq = fft.fftfreq(n, dt)
	plt.plot(fft_freq[:len(fft_freq)//2], abs(ft[:len(vec)//2]))
	plt.xlabel('freq (Hz)')
	plt.ylabel('|Y(freq)|')
	
	plt.tight_layout()
	plt.show()	
	return ft


def get_derivative(f_t, dt = 1):
# 	f_t_plus_dt = f_t[1:]
# 	f_t_plus_dt = np.append(f_t_plus_dt,0)
# 	print(f_t)
# 	print(f_t_plus_dt)
# 	first_derivative = (f_t_plus_dt-f_t)/dt
	f_t = gaussian_filter1d(f_t, 5)
	first_derivative = np.gradient(f_t)
	second_derivative = np.gradient(first_derivative)
	third_derivative = np.gradient(second_derivative)
# 	first_derivative = [(f_t[i+1] - f_t[i])/dt for i in range(len(f_t)-1)]
# 	first_derivative = np.append(first_derivative,first_derivative[-1])
# 	second_derivative = [(first_derivative[i+1] - first_derivative[i])/dt for i in range(len(first_derivative)-1)]
# 	second_derivative = np.append(second_derivative,second_derivative[-1])
	return first_derivative, second_derivative, third_derivative


def get_index_at_half_value(vec, cutoff_type = 'width_at_half_peak', eps = 2, stim_start = 110):
# 	half_max = max(vec) / 2.
# 	d = np.sign(half_max - vec[0:-1]) - np.sign(half_max - vec[1:])
# 	left_indx = np.where(d > 0)[0]
# 	right_indx_1 = np.where(d < 0)[-1]
 
# 	print(len(peaks[0]))
	#print(np.shape(vec))
	if cutoff_type == '1st_derivative':
		abs_derivative,_,_ = np.abs(get_derivative(vec))
		pos_abs_deriv = np.where(abs_derivative > eps)
		right_indx = pos_abs_deriv[-1]
		if len(right_indx)>0:
			right_indx = right_indx[-1]
		else:
			right_indx = 1

	elif cutoff_type == 'width_at_half_peak':
		vec_normed = vec - vec[stim_start]#vec[0] is v_rest
		#abs_vec = np.abs(vec)
		abs_vec_normed = np.abs(vec_normed)
		peaks = find_peaks(abs_vec_normed)
		if len(peaks[0]) == 0:
			peaks = [[stim_start]]
		peak_index = peaks[0][-1]
		right_indx = np.argmin(np.abs(abs_vec_normed[peak_index:]-abs_vec_normed[peak_index]/2)) + peak_index

	elif cutoff_type == 'trace_value':
# 		eps = 0.04
		abs_trace = np.abs(vec)
		pos_abs_trace = np.where(abs_trace > eps)
		
		right_indx = pos_abs_trace[-1]
		if len(right_indx)>0:
			right_indx = right_indx[-1]
		else:
			right_indx = stim_start
# 		print(right_indx)
	return right_indx, peak_index


def find_threshold(vec, eps = 1):
# 	print(vec)
	derivative,_,_ = get_derivative(vec)
# 	print(derivative)
	peaks = find_peaks(derivative)
	index = peaks[0]
# 	print(peaks)
# 	print(index)
# 	pos_deriv = np.where(derivative > eps)
# 	print(pos_abs_deriv)
# 	index = pos_deriv[0]
# # 	print(index)
	if len(index)>0:
		index = index[0]
	else:
		index = 0
	return index


def plot_single_current(current, weights, n, ax, eps):
	ax.plot(current, label = 'weights ='+str(weights[n][0])+','+str(weights[n][1]))
	right_index = get_index_at_half_value(current, eps)
	ax.scatter(right_index, current[right_index])
	return
	
	
def find_dend_sec_path(sec):
	'''
	function created to find the connecting path between a dendrite section number and the soma for the ball_and_many_sticks model
	'''
	path = [sec]
	while sec !=0:
# 		print(sec)
		sec = np.int0(np.ceil(sec/2)-1)
		path.insert(0,sec)
	return path


def find_forward_brach_sec_numbers(sec):
	'''
	unnecessary: can use dend_sec.children()
	function created to find the branch numbers connected (forward, not back) to a dendrite section for the ball_and_many_sticks model
	'''
	return (2*sec+1, 2*sec+2)


def transpose_value_transform(value, rows, cols):
	'''
	takes a value of a matrix and returns the transpose value
	e.g. on a 4x4, 2 becomes 5, 13 becomes 4, etc.
	'''
	orig_matrix = np.matrix([list(range(1 + cols * i, 1 + cols * (i + 1))) for i in range(rows)])
	new_val = np.array(orig_matrix[np.where(np.transpose(orig_matrix) == value)])[0][0]
	return new_val


def calculate_resistance(Ra, neck_diam, neck_length):
	'''
	function to calculate spine resistance
	Rn = (4 * axial resistivity * neck length)/(pi * neck diameter)
	input units: [ohm*cm, um, um]
	output units: M*ohm
	'''
	Rn = (4 * Ra * 1e-2 * neck_length * 1e-6)/(np.pi * neck_diam**2 * 1e-12)
	return Rn*1e-6

# calculate_resistance(150, 0.05, 1.5)
# 1145.9155902616462

# calculate_resistance(150, 0.075, 1.5)
# 509.2958178940651

# calculate_resistance(150, 0.5, 1.5)
# 11.459155902616464

# calculate_resistance(150, 0.17, 1.5)
# 99.12764621640537

def get_diam_from_MegaOhms(Rn, neck_length = 1.5, Ra = 150):
#Rn *1e-6  = (4 * Ra * 1e-2 * neck_length * 1e-6)/(np.pi * (neck_diam*1e-6)**2)
	return np.sqrt((4 * Ra * 1e-2 * neck_length * 1e-6)/(1e6*Rn*np.pi*1e-12))


