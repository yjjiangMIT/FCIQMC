import random_h
import random_vec
import walker
import det
import math
import global_var
import matplotlib
matplotlib.rcParams['backend'] = 'Qt4Agg'
import matplotlib.pyplot as pyplot

h = random_h.RandomH(global_var.rand_hamiltonian_flag)
recursion_num_list = []
unsigned_num_list = []

def run_fciqmc(h, max_recursion, recursion_num_list, unsigned_num_list):
	global_var.change_shift_flag = False
	# h.disp_mat()
	
	# Creates a single walker on D0.
	Walker = walker.Walker
	parent_list = []
	w = Walker(h, det.Det(h.basis_list[0]))
	w.signed_num = 2
	parent_list.append(w)
	count = 0
	j = 1
	unsigned_num_1 = global_var.nc
	correction = 0
	osc_count = 0
	corr_energy = 0
	aver_sample_size = max_recursion / 5
	ground_vec = [0] * global_var.dim
	while count < max_recursion:
		count += 1
		do_recursion(h, parent_list, recursion_num_list, unsigned_num_list)
		
		if Walker.unsigned_num > global_var.nc:
			global_var.change_shift_flag = True
		if global_var.change_shift_flag:
			if j == global_var.change_shift_step:
				unsigned_num_0 = Walker.unsigned_num
				old_correction = correction
				correction = -global_var.damping / global_var.change_shift_step / global_var.tau\
					* math.log(unsigned_num_0 / float(unsigned_num_1))
				global_var.shift += correction
				if correction * old_correction < 0 and global_var.damping > 5e-3:
					osc_count += 1
				if osc_count == 10:
					# Oscillation.
					global_var.damping /= 1.2
					osc_count = 0
				unsigned_num_1 = unsigned_num_0
				
				# For test only.
				print 'Shift  ', global_var.shift
				print 'Damping', global_var.damping
				
				j = 0
			j += 1
		if count > max_recursion - aver_sample_size - 1:
			corr_energy += global_var.shift
			vec = Walker.cal_ground_eig_vec(parent_list)
			if len(vec) == global_var.dim:
				for k in range(global_var.dim):
					ground_vec[k] += vec[k]
	corr_energy /= aver_sample_size
	total = 0
	for entry in ground_vec:
		total += entry ** 2
	norm = math.sqrt(total)
	for i in range(global_var.dim):
		ground_vec[i] /= norm
	print ''
	print corr_energy
	random_vec.RandomVec.disp_vec(ground_vec)
	
	draw_unsigned_num(recursion_num_list, unsigned_num_list)
	return parent_list

def do_recursion(h, parent_list, recursion_num_list, unsigned_num_list):
	Walker = walker.Walker
	child_list = Walker.spawn_and_die_for_all(h, parent_list)
	Walker.merge_walkers(parent_list, child_list)
	Walker.count_unsigned_num(parent_list)
	recursion_num_list.append(len(recursion_num_list)+1)
	unsigned_num_list.append(Walker.unsigned_num)
	
	# For test only.
	print 'Number of recursions', recursion_num_list[-1]
	print 'Amount of walkers', Walker.unsigned_num

def draw_unsigned_num(recursion_num_list, unsigned_num_list):
	pyplot.plot(recursion_num_list, unsigned_num_list)
	pyplot.show()