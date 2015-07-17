import det
import det_ops
import integrals
import math
import key_ops
import ctrl_panel
import test
import visual

import matplotlib
matplotlib.rcParams['backend'] = 'Qt4Agg'
import matplotlib.pyplot as plt

def run():
	"""Main FCIQMC function."""
	
	# Initialization.
	
	bit_num = integrals.para_list[0]
	chunk_num = integrals.para_list[1]
	dim = len(ctrl_panel.key_list)
	
	change_shift_crit_num = ctrl_panel.change_shift_crit_num
	w_num = ctrl_panel.init_walker_num
	init_w_num = w_num
	max_iter_num = ctrl_panel.max_iter_num
	
	damping = ctrl_panel.damping
	ref_key = ctrl_panel.ref_key
	shift = ctrl_panel.init_shift
	tau = ctrl_panel.tau
	
	change_shift_step = ctrl_panel.change_shift_step
	update_plots_step = ctrl_panel.update_plots_step
	
	exact_corr = ctrl_panel.exact_corr
	
	aver_flag = False
	change_shift_flag = False
	
	aver_numer = 0
	aver_denom = 1
	aver_shift = 0
	aver_times = 0
	
	old_w_num = change_shift_crit_num
	old_shift = shift
	
	dets_p = dict({ref_key : det.Det(w_num, True)})
	dets_p[ref_key].diag_entry = 0
	dets_p_old = {}
	for key in dets_p:
		dets_p_old[key] = det.Det(dets_p[key].value, True)
	
	aver_iter_num_list = []
	aver_proj_list = []
	aver_shift_list = []
	shift_list = []
	init_w_num_list = []
	iter_num_list = []
	log_init_w_num_list = []
	log_w_num_list = []
	proj_list = []
	w_num_list = []
	
	spawn_map = visual.new_spawn_map(dim)
	
	for iter_num in range(0, max_iter_num+1):	
		
		if iter_num % update_plots_step == 0:
			
			iter_num_list.append(iter_num)
			
			# Figure 4: 2D spawning map.
			visual.plot_update_2D(dim, spawn_map)
			
			# Figure 0: log of walker number.
			log_w_num_list.append(math.log(w_num))
			log_init_w_num_list.append(math.log(init_w_num))
			data = (iter_num_list, log_w_num_list, log_init_w_num_list)
			color = ('b', 'r')
			ax_range = [0, iter_num] + ctrl_panel.y_axis_log_w_num_plot
			label = ('Iteration', 'Log of walker number')
			visual.plot_update_1D(0, data, color, ax_range, label)
			
			# Figure 1: walker number.
			w_num_list.append(w_num)
			init_w_num_list.append(init_w_num)
			data = (iter_num_list, w_num_list, init_w_num_list)
			color = ('b', 'r')
			ax_range = [0, iter_num, 0, w_num*1.2]
			label = ('Iteration', 'Walker number')
			visual.plot_update_1D(1, data, color, ax_range, label)
			
			# Figure 2: correlation energy.
			shift_list.append(shift)
			numer, denom = det_ops.corr_by_proj(dets_p, ref_key)
			proj_energy = numer / denom
			proj_list.append(proj_energy)
			if aver_flag:
				aver_iter_num_list.append(iter_num)
				aver_shift_list.append(aver_shift)
				aver_proj_list.append(aver_proj)
			data = (iter_num_list, shift_list, proj_list, [exact_corr]*len(iter_num_list))
			color = ('b', 'g', 'r')
			label = ('Iteration', 'Energy')
			ax_range = [0, iter_num] + ctrl_panel.y_axis_energy_plot
			visual.plot_update_1D(2, data, color, ax_range, label)
			
			# Figure 3: |Psi(t)>.
			vec = det_ops.dets_2_vec(dets_p)
			data = (range(len(vec)), vec)
			color = ('b',)
			label = ('Dimension', 'Value')
			ax_range = [0, dim] + ctrl_panel.y_axis_eig_vec_plot
			visual.plot_update_1D(3, data, color, ax_range, label)
			
			# ref_key_cand = det_ops.find_ref(dets_p);
			
			print iter_num, w_num, dets_p[ref_key].value, shift, proj_energy
			
		if (not change_shift_flag) and w_num > change_shift_crit_num:
			change_shift_flag = True
			shift = proj_energy
			crit_iter_num = iter_num
		if change_shift_flag:
			if iter_num % change_shift_step == 0:
				correction = -damping / change_shift_step / tau \
					* math.log(w_num / float(old_w_num))
				shift += correction
				old_w_num = w_num
				if (not aver_flag) and (iter_num - crit_iter_num) > ctrl_panel.wait_for_aver_num:
					if abs(shift - old_shift) < 0.02:
						aver_flag = True
				if aver_flag:
					aver_shift = (aver_shift*aver_times+shift) / (aver_times+1)
					aver_numer = (aver_numer*aver_times+numer) / (aver_times+1)
					aver_denom = (aver_denom*aver_times+denom) / (aver_times+1)
					aver_proj = aver_numer / aver_denom
					aver_times += 1
				old_shift = shift
		
		det_ops.single_step(dets_p, spawn_map, tau, shift, True)
		w_num = det_ops.count_u_num(dets_p)
		init_w_num = abs(dets_p[ref_key].value)
		
	# return aver_shift, aver_proj, dets_p, vec

def smooth(x_list, half_width):
	"""Smoothens x_list by averaging over a range of 2*half_width+1."""
	
	y_list = []
	for i in range(half_width, len(x_list)-half_width):
		average = 0
		for j in range(i-half_width, i+half_width+1):
			average += x_list[j]
		average /= (2*half_width + 1.0)
		y_list.append(average)
	return y_list

def derivative(x_list, half_width):
	"""Takes derivative of x_list over a range of 2*half_width+1."""
	
	y_list = []
	for i in range(half_width, len(x_list)-half_width):
		y_list.append(x_list[i+half_width]-x_list[i-half_width])
	return y_list

def draw_curve(x_list, y_list, color):
	"""Draws a curve."""
	
	plt.plot(x_list, y_list, color)