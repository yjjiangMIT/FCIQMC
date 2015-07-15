import det
import det_ops
import integrals
import math
import key_ops
import ctrl_panel
import test

import matplotlib
matplotlib.rcParams['backend'] = 'Qt4Agg'
import matplotlib.pyplot as plt

def run():
	"""Main FCIQMC function."""
	
	# Initialization.
	
	bit_num = integrals.para_list[0]
	chunk_num = integrals.para_list[1]
	
	change_shift_crit_num = ctrl_panel.change_shift_crit_num
	init_walker_num = ctrl_panel.init_walker_num
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
	w_num = init_walker_num
	
	dets_p = dict({ref_key : det.Det(w_num, True)})
	dets_p[ref_key].diag_entry = 0
	dets_p_old = {}
	for key in dets_p:
		dets_p_old[key] = det.Det(dets_p[key].value, True)
	key_list = test.gen_key_list()
	
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
	
	# Figure 0: walker number in log.
	init_figure(0, [0, 20] + ctrl_panel.y_axis_log_w_num_plot, 'Iteration', 'Log of walker number')
	
	# Figure 1: walker number.
	init_figure(1, [0, 20, 0, w_num*1.2], 'Iteration', 'Walker number')
	
	# Figure 2: energy.
	init_figure(2, [0, 20] + ctrl_panel.y_axis_energy_plot, 'Iteration', 'Energy')
	
	# Figure 3: |Psi(t)>.
	init_figure(3, [0, len(key_list)] + ctrl_panel.y_axis_eig_vec_plot, 'Dimension', 'Value')
	
	for iter_num in range(1, max_iter_num+1):
		det_ops.single_step(dets_p, tau, shift)
		w_num = det_ops.count_u_num(dets_p)
		init_w_num = abs(dets_p[ref_key].value)
		
		if iter_num % update_plots_step == 0:
			
			iter_num_list.append(iter_num)
			
			log_w_num_list.append(math.log(w_num))
			log_init_w_num_list.append(math.log(init_w_num))
			# Figure 0: walker number.
			plt.figure(0)
			plt.clf()
			plt.axis([0, iter_num] + ctrl_panel.y_axis_log_w_num_plot)
			# Draws log of walker number.
			draw_curve(iter_num_list, log_w_num_list, 'b')
			draw_curve(iter_num_list, log_init_w_num_list, 'r')
			plt.xlabel('Iteration')
			plt.ylabel('Log of walker number')
			plt.pause(0.01)
			
			w_num_list.append(w_num)
			init_w_num_list.append(init_w_num)
			# Figure 1: walker number.
			plt.figure(1)
			plt.clf()
			plt.axis([0, iter_num, 0, w_num*1.2])
			# Draws walker number.
			draw_curve(iter_num_list, w_num_list, 'b')
			draw_curve(iter_num_list, init_w_num_list, 'r')
			plt.xlabel('Iteration')
			plt.ylabel('Walker number')
			plt.pause(0.01)
			
			shift_list.append(shift)
			numer, denom = det_ops.corr_by_proj(dets_p, ref_key)
			proj_energy = numer / denom
			proj_list.append(proj_energy)
			if aver_flag:
				aver_iter_num_list.append(iter_num)
				aver_shift_list.append(aver_shift)
				aver_proj_list.append(aver_proj)
			# Figure 2: correlation energy.
			plt.figure(2)
			plt.clf()
			plt.axis([0, iter_num] + ctrl_panel.y_axis_energy_plot)
			# Draws shift.
			draw_curve([0, max_iter_num], [exact_corr, exact_corr], 'k')
			draw_curve(iter_num_list, shift_list, 'b')
			draw_curve(iter_num_list, proj_list, 'r')
			draw_curve(aver_iter_num_list, aver_shift_list, 'c')
			draw_curve(aver_iter_num_list, aver_proj_list, 'm')
			plt.xlabel('Iteration')
			plt.ylabel('Energy')
			plt.pause(0.01)
			
			vec = det_ops.dets_2_vec(key_list, dets_p)
			# Figure 3: |Psi(t)>.
			plt.figure(3)
			plt.clf()
			plt.axis([0, len(key_list)] + ctrl_panel.y_axis_eig_vec_plot)
			# Draws eigenvector.
			draw_curve(range(1, len(vec)+1), vec, 'b')
			plt.xlabel('Dimension')
			plt.ylabel('Value')
			plt.pause(0.01)
			
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
	
	"""
	# Figure 3: Derivative of log(w_num).
	s_log_w_num_list = smooth(log_w_num_list, 2)
	d_log_w_num_list = derivative(s_log_w_num_list, 2)
	f3 = plt.figure(3)
	plt.axis([0, max_iter_num, -0.1, 0.2])
	plt.ion()
	plt.clf()
	draw_curve(iter_num_list[:len(s_log_w_num_list)], s_log_w_num_list, 'b')
	draw_curve(iter_num_list[:len(d_log_w_num_list)], d_log_w_num_list, 'r')
	plt.xlabel('Iteration')
	plt.ylabel('Derivative of log(w_num)')
	plt.show()
	"""
	
	return aver_shift, aver_proj, dets_p, vec

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

def init_figure(index, axis_range, x_label, y_label):
	"""Initializes a figure."""
	
	f = plt.figure(index)
	plt.axis(axis_range)
	plt.ion()
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.show()

def draw_curve(x_list, y_list, color):
	"""Draws a curve."""
	
	plt.plot(x_list, y_list, color)