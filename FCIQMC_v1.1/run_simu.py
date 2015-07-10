def run(orb_num, e_num, max_iter_num, crit_walker_num, \
	avg_shift_iter_num = 0, change_shift_step = 5, cutoff = 0.01, \
	damping = 0.005, h = None, init_crit_w_num = 10, init_shift = 0, tau = 0.0005):
	"""Main FCIQMC function."""
	
	import file_io
	import math
	import vec_ops
	import w_dist_ops
	import walker
	import matplotlib
	matplotlib.rcParams['backend'] = 'Qt4Agg'
	import matplotlib.pyplot as plt
		
	### Initialization.
	# Read Hamiltonian from file.
	if h == None:
		h = file_io.read_hamiltonian(orb_num, e_num)
	exact_gnd_state = h.exact_gnd_state() # The exact eigenvector.
	exact_corr = -h.ref_energy # The exact correlation energy.
	
	# Intializes walker distribution as D0 and puts 10 walkers on it.
	# Picks the sign so that it is more likely to converge to the right-signed solution.
	init_walker_num = 50
	w_dist_p = {} # Walker distribution.
	try:
		init_sign = (exact_gnd_state[h.ref_key] > 0)
	except KeyError:
		init_sign = True
	if init_sign:
		w_dist_p[h.ref_key] = walker.Walker(init_walker_num, True)
	else:
		w_dist_p[h.ref_key] = walker.Walker(-init_walker_num, True)
	
	# Some parameters.
	if avg_shift_iter_num == 0:
		avg_shift_iter_num = int(max_iter_num/10) # Number of iterations for an averaged shift.
	shift = init_shift
	avg_shift = 0 # Averaged shift.
	curr_walker_num = init_walker_num # Current unsigned number of walkers.
	prev_walker_num = crit_walker_num # Previous unsigned number of walkers.
	
	curr_walker_num_list = [curr_walker_num] # List of walker numbers.
	shift_list = [0] # List of shift values.
	corr_proj_list = [0] # List of projected correlation energies.
	
	change_shift_flag = False # The flag of starting to change the shift.
	
	### Figures.
	# Figure 0: eigenvector.
	f0 = plt.figure(0)
	plt.axis([0, h.para_list[4], -1, 1])
	plt.ion()
	plt.xlabel('Dimension')
	plt.ylabel('Value')
	plt.show()
	
	# Figure 1: correlation energy.
	f1 = plt.figure(1)
	min_shift = -h.ref_energy # It is used to decide the y-axis range of figure 1.
	plt.axis([0, max_iter_num, min_shift*1.5, 0.5])
	plt.ion()
	plt.xlabel('Iteration')
	plt.ylabel('Corr. energy')
	plt.show()
	
	# Figure 2: walker number.
	f2 = plt.figure(2)
	max_curr_walker_num = curr_walker_num # It is used to decide the y-axis range of figure 2.
	plt.axis([0, max_iter_num, 0, max_curr_walker_num*1.2])
	plt.ion()
	plt.xlabel('Iteration')
	plt.ylabel('Walker number')
	plt.show()
	
	print 'Press enter to start simulation.'
	raw_input() # Just a pulse.
	
	### Iterations.
	for count in range(1, max_iter_num+1):
		w_dist_ops.w_dist_single_step(h, w_dist_p, tau, shift, init_crit_w_num) # Spawning, death and annihilation.
		
		# Walker number.
		curr_walker_num = w_dist_ops.w_dist_count_u_num(w_dist_p)
		curr_walker_num_list.append(curr_walker_num)
		
		# Projected energy.
		proj_corr = w_dist_ops.w_dist_corr_by_proj(h, w_dist_p)
		corr_proj_list.append(proj_corr)
		
		# Shift.
		shift_list.append(shift)
		
		# Outputs some information.
		# disp_info(count, curr_walker_num, shift, proj_corr)
		
		if curr_walker_num > max_curr_walker_num:
			max_curr_walker_num = curr_walker_num
		
		# Updates the plots every $max_iter_num/200$ iterations.
		if count % (math.ceil(max_iter_num/200.0)) == 0:
			# Calculates normalized ground state from walker distribution.
			cal_gnd_state = w_dist_ops.w_dist_2_vec(w_dist_p)
			# Disgards entries below some cutoff value.
			vec_ops.vec_dom_entries(cal_gnd_state, cutoff)
			
			# Figure 0: eigenvector.
			plt.figure(0)
			plt.clf()
			plt.axis([0, h.para_list[4], -1, 1])
			# Changes sparse distributions into dense vectors.
			exact_dense_vec = vec_ops.vec_sparse_2_dense(h.key_list, exact_gnd_state, h.para_list[4])
			cal_dense_vec = vec_ops.vec_sparse_2_dense(h.key_list, cal_gnd_state, h.para_list[4])
			# Draws exact and calculated eigenvectors.
			draw_curve(range(h.para_list[4]), exact_dense_vec, 'r')
			draw_curve(range(h.para_list[4]), cal_dense_vec, 'b')
			plt.xlabel('Dimension')
			plt.ylabel('Value')
			plt.pause(0.01)
			
			# Figure 1: correlation energy.
			plt.figure(1)
			plt.clf()
			plt.axis([0, max_iter_num, min_shift*1.5, 0.5])
			# Draws exact, shift and projected correlation energies.
			draw_curve(range(0, max_iter_num+1), [exact_corr] * (max_iter_num+1), 'r')
			draw_curve(range(0, count+1), shift_list, 'b')
			draw_curve(range(0, count+1), corr_proj_list, 'g')	
			plt.xlabel('Iteration')
			plt.ylabel('Corr. energy')
			plt.pause(0.01)
			
			# Figure 2: walker number.
			plt.figure(2)
			plt.clf()
			plt.axis([0, max_iter_num, 0, max_curr_walker_num*1.2])
			# Draws walker number.
			draw_curve(range(0, count+1), curr_walker_num_list, 'b')
			plt.xlabel('Iteration')
			plt.ylabel('Walker number')
			plt.pause(0.01)
			
		if curr_walker_num > crit_walker_num:
			# Starts to change the shift.
			change_shift_flag = True
		
		# Changes the shift every $change_shift_step$ iterations.
		if change_shift_flag:
			if count % change_shift_step == 0:
				correction = -damping / change_shift_step / tau\
					* math.log(curr_walker_num / float(prev_walker_num))
				shift += correction
				if shift < min_shift:
					min_shift = shift
				prev_walker_num = curr_walker_num
		
		# Averages the shift.
		if count > max_iter_num - avg_shift_iter_num:
			avg_shift += shift
	
	### End of the loop.
	# Averages the shift.
	avg_shift /= avg_shift_iter_num
	
	return h, cal_gnd_state, avg_shift
      
def make_h(orb_num, e_num, bit_num, sparsity):
	"""Generates a random sparse Hamiltonian."""
	
	import file_io
	
	file_io.write_hamiltonian(orb_num, e_num, bit_num, sparsity)
	h = file_io.read_hamiltonian(orb_num, e_num)
	return h

def disp_info(count, walker_num, shift, proj_corr):
	"""Displays some information."""
	
	print 'Iteration number',
	print count
	print 'Walker number',
	print walker_num
	print 'Shift',
	print shift
	print 'Projected correlation energy',
	print proj_corr

def draw_curve(x_list, y_list, color):
	"""Draws a curve."""
	
	import matplotlib
	matplotlib.rcParams['backend'] = 'Qt4Agg'
	import matplotlib.pyplot as plt
	
	plt.plot(x_list, y_list, color)