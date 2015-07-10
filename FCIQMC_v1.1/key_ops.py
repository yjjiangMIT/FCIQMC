# Provides all operations on a key.
# Keys are strings like '4, 8, 2'. The example means a configuration [0100,1000,0010].
# Keys are used to label vector (Hash table) entries.

def occup_2_indices(occup, para_list):
	"""Converts from occup to indices. Inverse of indices_2_occup."""
	
	# e.g. from [4,8,2] = [(0100)_dec,(1000)_dec,(0010)_dec] to 0100,1000,0010 = [10,7,1].
	
	bit_num = para_list[0]
	chunk_num = para_list[1]
	indices = []
	for i in range(len(occup)):
		for j in range(bit_num-1, -1, -1):
			if (occup[i]>>j) & 1:
				indices.append((chunk_num-i-1) * bit_num + j)
	return indices
	
def indices_2_occup(indices, para_list):
	"""Converts from indices to occup. Inverse of occup_2_indices."""
	
	# e.g. from [10,7,1] = 0100,1000,0010 to [(0100)_dec,(1000)_dec,(0010)_dec] = [4,8,2].
	
	bit_num = para_list[0]
	chunk_num = para_list[1]
	occup = [0] * chunk_num
	for i in indices:
		chunk_index = chunk_num - i/bit_num - 1
		power = i % bit_num
		occup[chunk_index] |= 2**power
	return occup
	
def next_config_indices(curr_indices, para_list, first_run_flag = True):
	"""Computes the next available configuration, in the form of indices."""
	
	# Configurations are arranged in dictionary order.
	# e.g. [1,0]->[2,0]->[2,1]->[3,0]->[3,1]->[3,2]->...->[5,4]
	
	spin_orb_num = para_list[2] * 2
	e_num = para_list[3]
	if first_run_flag and curr_indices[-1] == spin_orb_num - e_num:
		return None
	else:
		next_indices = curr_indices[:]
		if len(next_indices) == 1:
			return [next_indices[0] + 1]
		elif next_indices[-1] == next_indices[-2] - 1:
			next_indices = next_config_indices(next_indices[:-1], para_list, False)\
				+ [e_num - len(next_indices)]
		else:
			next_indices[-1] += 1
	return next_indices
	
def next_config_occup(curr_occup, para_list):
	"""Computes the next available configuration, in the form of occup."""
	
	curr_indices = occup_2_indices(curr_occup, para_list)
	next_indices = next_config_indices(curr_indices, para_list)
	if next_indices == None:
		return None
	else:
		next_occup = indices_2_occup(next_indices, para_list)
		return next_occup

def create_keys(key_list, para_list):
	"""Create an ordered list of keys."""
	
	e_num = para_list[3]
	dim = para_list[4]
	indices = range(e_num-1, -1, -1)
	occup = indices_2_occup(indices, para_list)
	for i in range(dim):   
		key = str(occup)
		key = key[1:-1]
		key_list.append(key)
		occup = next_config_occup(occup, para_list)
	key_list.sort()

def key_2_index(key_list, key):
	"""Convert from key to index in the matrix."""
	
	# Can use bisearch, but for now just use build-in search.
	
	try:
		index = key_list.index(key)
		return index
	except ValueError:
		return None

def excited_key_and_prob(h, key_gnd):
	"""From a ground state, gets a possible excited state, and the probability rate."""
	
	import key_ops
	import random
	
	index = key_ops.key_2_index(h.key_list, key_gnd)
	vec_excited_list = h.matrix[index]
	choice_num = len(vec_excited_list) - 1
	if choice_num > 0:
		# Connected to others.
		key_exc = key_gnd
		while key_exc == key_gnd:
			key_exc = random.sample(vec_excited_list, 1)[0]
		entry = vec_excited_list[key_exc]
		prob_rate = choice_num * entry
		return key_exc, prob_rate
	else:
		return None, 0

def dead_prob(h, key):
	"""Gets the probability rate of death."""
	
	import key_ops
	
	index = key_ops.key_2_index(h.key_list, key)
	vec_excited_list = h.matrix[index]
	try:
		diag_entry = vec_excited_list[key]
	except KeyError:
		diag_entry = 0
	prob_rate = diag_entry - h.ref_energy
	return prob_rate