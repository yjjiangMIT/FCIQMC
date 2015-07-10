import det
import det_ops
import black_box
import math
import key_ops
import gl_consts

import matplotlib
matplotlib.rcParams['backend'] = 'Qt4Agg'
import matplotlib.pyplot as plt

def run(old_vec = []):
	"""Main FCIQMC function."""
	
	### Initialization.
	
	# Intializes walker distribution as D0 and puts 10 walkers on it.
	ref_key = (8, 10, 6, 4)
	
	tau = gl_consts.tau
	shift = gl_consts.shift
	max_iter_num = gl_consts.max_iter_num
	bit_num = gl_consts.para_list[0]
	chunk_num = gl_consts.para_list[1]
	
	shift = -6
	old_shift = shift
	gl_consts.init_crit_w_num = 0
	change_shift_step = 5
	update_plots_step = 20
	damping = 0.05
	change_shift_flag = False
	aver_flag = False
	change_shift_crit_num = 50000
	old_w_num = change_shift_crit_num
	tau = 0.002
	exact_gnd = -8.892162249
	
	aver_times = 0
	aver_shift = 0
	aver_numer = 0
	aver_denom = 1
	
	max_iter_num = 100000
	
	init_walker_num = 50
	init_walker = det.Det(init_walker_num, True)
	w_num = init_walker_num
	dets_p = {}
	det_ops.merge(dets_p, dict({ref_key : init_walker}))
	ref_energy = dets_p[ref_key].diag_entry
	
	w_num_list = []
	log_w_num_list = []
	
	iter_num_list = []
	aver_iter_num_list = []
	shift_list = []
	aver_shift_list = []
	proj_list = []
	aver_proj_list = []
	
	dets_p_old = {}
	for key in dets_p:
		dets_p_old[key] = det.Det(dets_p[key].value, True)
	
	# Figure 0: walker number.
	f0 = plt.figure(0)
	# plt.axis([0, max_iter_num, 0, 20000])
	plt.axis([0, max_iter_num, 3, 11])
	plt.ion()
	plt.xlabel('Iteration')
	# plt.ylabel('Walker number')
	plt.ylabel('Log of walker number')
	plt.show()
	
	# Figure 1: shift.
	f2 = plt.figure(1)
	plt.axis([0, max_iter_num, -12, 0])
	plt.ion()
	plt.xlabel('Iteration')
	plt.ylabel('Energy')
	plt.show()
	
	""""""
	# Figure 2: |Psi(t)>.
	f2 = plt.figure(2)
	plt.axis([0, (2**bit_num)**chunk_num - 1, -1, 1])
	plt.ion()
	plt.xlabel('Dimension')
	plt.ylabel('Value')
	plt.show()
	""""""
	
	for iter_num in range(1, max_iter_num+1):
		det_ops.single_step(dets_p, tau, shift)
		w_num = det_ops.count_u_num(dets_p)
		
		if iter_num % update_plots_step == 0:
			
			iter_num_list.append(iter_num)
			
			log_w_num_list.append(math.log(w_num))
			# Figure 0: walker number.
			plt.figure(0)
			plt.clf()
			plt.axis([0, max_iter_num, 3, 11])
			# Draws walker number.
			draw_curve(iter_num_list, log_w_num_list, 'b')
			plt.xlabel('Iteration')
			plt.ylabel('Log of walker number')
			plt.pause(0.01)
			
			shift_list.append(shift)
			numer, denom = det_ops.corr_by_proj(dets_p, ref_key)
			proj_energy = numer / denom + ref_energy
			proj_list.append(proj_energy)
			if aver_flag:
				aver_iter_num_list.append(iter_num)
				aver_shift_list.append(aver_shift)
				aver_proj_list.append(aver_proj)
			# Figure 1: shift.
			plt.figure(1)
			plt.clf()
			plt.axis([0, max_iter_num, -12, 0])
			# Draws shift.
			draw_curve([0, max_iter_num], [exact_gnd, exact_gnd], 'k')
			draw_curve(iter_num_list, shift_list, 'b')
			draw_curve(iter_num_list, proj_list, 'r')
			draw_curve(aver_iter_num_list, aver_shift_list, 'c')
			draw_curve(aver_iter_num_list, aver_proj_list, 'm')
			plt.xlabel('Iteration')
			plt.ylabel('Energy')
			plt.pause(0.01)
			
			""""""
			vec = det_ops.dets_2_vec(dets_p)
			# Figure 2: |Psi(t)>.
			plt.figure(2)
			plt.clf()
			plt.axis([0, len(vec)-1, -1, 1])
			# Draws eigenvector.
			if len(old_vec) == len(vec):
				draw_curve(range(len(old_vec)), old_vec, 'r')
			draw_curve(range(len(vec)), vec, 'b')
			plt.xlabel('Dimension')
			plt.ylabel('Value')
			plt.pause(0.01)
			""""""
			
			ref_key_cand = det_ops.find_ref(dets_p);
			
			print iter_num, w_num, dets_p[ref_key].value, shift, proj_energy, ref_key_cand
			
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
				if (not aver_flag) and (iter_num - crit_iter_num) > 500:
					if abs(shift - old_shift) < 0.02:
						aver_flag = True
				if aver_flag:
					aver_shift = (aver_shift*aver_times+shift) / (aver_times+1)
					aver_numer = (aver_numer*aver_times+numer) / (aver_times+1)
					aver_denom = (aver_denom*aver_times+denom) / (aver_times+1)
					aver_proj = aver_numer / aver_denom + ref_energy
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
	
	return aver_shift, dets_p, vec

def smooth(x_list, half_width):
	y_list = []
	for i in range(half_width, len(x_list)-half_width):
		average = 0
		for j in range(i-half_width, i+half_width+1):
			average += x_list[j]
		average /= (2*half_width + 1.0)
		y_list.append(average)
	return y_list

def derivative(x_list, half_width):
	y_list = []
	for i in range(half_width, len(x_list)-half_width):
		y_list.append(x_list[i+half_width]-x_list[i-half_width])
	return y_list

def draw_curve(x_list, y_list, color):
	"""Draws a curve."""
	
	plt.plot(x_list, y_list, color)

def gen_matrix():
	"""Creates a MATLAB file of the full Hamiltonian as a sparse matrix."""
	
	dets_p = {}
	dets_c = {}
	
	spatial_orbs_list = []
	key_list = []
	
	for i in range(0, 8):
		for j in range(i+1, 8):
			for k in range(j+1, 8):
				spatial_orbs_list.append([k, j, i])
	for sp_orbs_i in spatial_orbs_list:
		for sp_orbs_j in spatial_orbs_list:
			sp_orbs_temp = sp_orbs_j[:]
			for k in range(3):
				sp_orbs_temp[k] += 8
			orbs = tuple(sp_orbs_temp + sp_orbs_i)
			key = key_ops.orbs_2_key(orbs)
			key_list.append(key)
	for key in key_list:
		dets_c[key] = det.Det(1, True)
	
	det_ops.merge(dets_p, dets_c)
	
	f = open('matrix_m8n6.m', 'w')
	f.write('data = [')
	for i in range(len(key_list)):
		key_i = key_list[i]
		for j in range(len(key_list)):
			key_j = key_list[j]
			orbs_i, sign, orbs_diff = key_ops.difference(key_i, key_j)
			if orbs_i != None:
				entry = black_box.sandwich(orbs_i, orbs_diff)
				if entry == 0:
					continue
				elif not sign:
					entry = -entry
				f.write(str(i+1))
				f.write(' ')
				f.write(str(j+1))
				f.write(' ')
				f.write('{0:.4f}'.format(entry))
				f.write(' ')
	f.write('];\r\ndata = reshape(data, 3, length(data)/3);\r\ni = data(1,:);\r\nj = data(2,:);\r\n')
	f.write('s = data(3,:);\r\nH = sparse(i,j,s);\r\neig_vals = eigs(H);\r\ngnd = min(eig_vals);')
	f.close()
	
	return dets_p