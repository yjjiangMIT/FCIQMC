import det
import ctrl_panel
import key_ops
import random
import integrals
import math
import test

# Provides all operations on a dictionary whose values are Det objects.
# The keys are tuples whose bit representations correspond to orbital occupations.

def merge(dets_p, dets_c):
	"""Merges dets_c into dets_p."""
	
	for key in dets_c:
		if key in dets_p:
			# Parent list includes this determinant.
			dets_p[key].value += dets_c[key].value
			if dets_p[key].value == 0:
				del dets_p[key]
			else:
				dets_p[key].flag = (abs(dets_p[key].value) >= ctrl_panel.init_crit_w_num)
		else:
			# Parent list does not include this determinant.
			if dets_c[key].if_survive_as_child():
				# The merged walker survives if the child meets survival criterions.
				dets_p[key] = dets_c[key]
				dets_p[key].set_diag_entry(key)
				dets_p[key].flag = (abs(dets_p[key].value) >= ctrl_panel.init_crit_w_num)

def spawn(dets_p, dets_c, key_p, spawn_map, tau, update_spawn_map_flag):
	"""Spawning of all parents on the determinant indexed by key_p."""
	
	# Spawning.
	p_sign = (dets_p[key_p].value > 0)
	p_u_num = abs(dets_p[key_p].value)
	count = 0
	
	for count in range(p_u_num):
		# Single or double excitation.
		rand_val = random.random()
		if rand_val < ctrl_panel.single_prob:
			# Single excitation.
			orbs_p, key_c, sign, p_gen, orbs_diff = key_ops.single_excite(key_p)
			p_sing_or_doub = ctrl_panel.single_prob
		else:
			# Double excitation.
			orbs_p, key_c, sign, p_gen, orbs_diff = key_ops.double_excite(key_p)
			p_sing_or_doub = 1 - ctrl_panel.single_prob
		mat_element = integrals.sandwich(orbs_p, orbs_diff)
		if not sign:
			# Accounts for a possible negative sign from permuting the spawned determinant.
			mat_element = -mat_element
		prob = tau * mat_element / p_gen / p_sing_or_doub
		
		# Gets the sign of children.
		if prob < 0:
			c_sign = p_sign
			prob = -prob
		else:
			c_sign = not p_sign
		
		c_num = int(prob)
		# At least we have c_num children.
		
		prob_frac = prob - c_num
		rand_val = random.random()
		if rand_val < prob_frac:
			# One more child.
			c_num += 1
		
		if c_num != 0:
			# Add c_num children to the list.
			if not c_sign:
				c_num = -c_num
			if key_c in dets_c:
				dets_c[key_c].value += c_num
				if dets_c[key_c].value == 0:
					del dets_c[key_c]
				else:
					dets_c[key_c].flag = dets_p[key_p].flag
			else:
				dets_c[key_c] =  det.Det(c_num, dets_p[key_p].flag)
			if update_spawn_map_flag:
				key_index_p = ctrl_panel.key_list.index(key_p)
				spawn_map[key_index_p][ctrl_panel.key_list.index(key_c)] += abs(c_num)
	
def die(dets_p, key_p, tau, shift):
	"""Dying/cloning of all parents on the determinant indexed by key."""
	
	prob = tau * (dets_p[key_p].diag_entry - shift)
	p_sign = (dets_p[key_p].value > 0)
	p_u_num = abs(dets_p[key_p].value)
	count = 0
	
	if prob < 0:
		# Cloning.
		prob = -prob
		prob_frac = prob - int(prob)
		clone_num = int(prob) * p_u_num
		# At least we have clone_num clones.
		
		for count in range(p_u_num):
			rand_val = random.random()
			if rand_val < prob_frac:
				# One more cloning event.
				clone_num += 1
		if clone_num != 0:
			if not p_sign:
				clone_num = -clone_num
			dets_p[key_p].value += clone_num
	else:
		# Dying.
		for count in range(p_u_num):
			rand_val = random.random()
			if rand_val < prob:
				# Dead.
				if p_sign:
					dets_p[key_p].value -= 1
				else:
					dets_p[key_p].value += 1
	if dets_p[key_p].value == 0:
		del dets_p[key_p]

def annih(dets_p, dets_c):
	"""Overall annihilation step."""
	
	merge(dets_p, dets_c)

def single_step(dets_p, spawn_map, tau, shift, update_spawn_map_flag = False):
	"""The whole process of spawning, dying and annihilation."""
	
	dets_c = {}
	
	# Spawning and dying/cloning.
	p_keys = dets_p.keys()
	for key_p in p_keys:
		spawn(dets_p, dets_c, key_p, spawn_map, tau, update_spawn_map_flag)
		die(dets_p, key_p, tau, shift)
	
	# Annihilation.
	annih(dets_p, dets_c)

def count_u_num(dets):
	"""Count the unsigned number of walkers."""
	
	total_u_num = 0
	for key in dets:
		total_u_num += int(abs(dets[key].value))
	return total_u_num

def dets_2_vec(dets):
	"""From a determinant list to a normalized eigenvector."""
	
	vec = [0] * len(ctrl_panel.key_list)
	norm_sq = 0
	for i in range(len(ctrl_panel.key_list)):
		key = ctrl_panel.key_list[i]
		if key in dets:
			norm_sq += dets[key].value ** 2
	for i in range(len(ctrl_panel.key_list)):
		key = ctrl_panel.key_list[i]
		if key in dets:
			vec[i] = dets[key].value / math.sqrt(norm_sq)
	return vec

def dot_prod(dets_1, dets_2, norm_flag):
	"""Calculate the dot product between two determinant lists."""
	
	prod = 0
	for key in dets_1:
		if key in dets_2:
			prod += dets_1[key].value * dets_2[key].value
	if norm_flag:
		norm_sq_1 = dot_prod(dets_1, dets_1, False)
		norm_sq_2 = dot_prod(dets_2, dets_2, False)
		prod = prod / math.sqrt(norm_sq_1) / math.sqrt(norm_sq_2)
	return prod

def corr_by_proj(dets_p, ref_key):
	"""Calculates correlation energy by projection."""
	
	# E_corr = sum_j!=o(<D_j|H|D_0>*(<N_j>/<N_0>))
	# Returns numerator and denominator separately, such that they can both be averaged over iterations.
	ref_orbs = key_ops.key_2_orbs(ref_key)
	numer = 0 # Initializes numerator.
	
	for key in dets_p:
		if key != ref_key:
			orbs_gnd, sign_exc, orbs_diff, count = key_ops.difference(ref_key, key)
			if orbs_gnd != None:
				term = integrals.sandwich(orbs_gnd, orbs_diff) * dets_p[key].value
				if not sign_exc:
					term = -term
				numer += term
	denom = float(dets_p[ref_key].value) # Denominator.
	return numer, denom

def find_ref(dets_p):
	"""Finds the determinant that has the most population."""
	
	# This is used when a true D0 is not applicable, such that an arbitrary determinant is used as reference.
	# By running this method after sufficient iterations, the true D0 can be obtained.
	max_num = 0
	for key in dets_p:
		if max_num < abs(dets_p[key].value):
			max_num = abs(dets_p[key].value)
			ref_key = key
	return ref_key