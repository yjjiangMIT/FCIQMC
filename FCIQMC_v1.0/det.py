import random
import load_data
import global_var

class Det:
	# Class of determinants.
	# occup: an integer list recording the occupation of orbitals.
	# e.g. if global_var.bit_num = 4, then the occupation 0100,1000,0010 corresponds to occup = [2,8,4].	
	chunk_num = global_var.spin_orb_num/global_var.bit_num + 1 # Number of chunks required to store a configuration.
	
	def __init__(self, occup):
		# Initiates a Det object.
		self.occup = occup[:]
	
	def get_indices(self):
		# Returns two integer lists:
		# (1) occup_indices: an integer list recording the indices of occupied orbitals.
		# e.g. if global_var.bit_num = 4, then the occupation 0100,1000,0010 corresponds to occup_indices = [1,7,10].
		# (2) unoccup_indices: an integer list recording the indices of unoccupied orbitals.
		# e.g. if global_var.bit_num = 4, then the occupation 0100,1000,0010 corresponds to unoccup_indices = [0,2,3,4,5,6,8,9,11].
		occup_indices = Det.occup_2_indices(self.occup)
		unoccup_indices = [i for i in range(global_var.spin_orb_num) if i not in occup_indices]
		return occup_indices, unoccup_indices
	
	def double_excite(self):
		# Generates a double-excited configuration.
		excited = Det(self.occup)
		ground_occup_indices, ground_unoccup_indices = self.get_indices()
		destructed = random.sample(ground_occup_indices, 2)
		created = random.sample(ground_unoccup_indices, 2)
		
		# For test only.
		# print 'Destruct', destructed
		# print 'Create  ', created
		
		for i in [0, 1]:
			dest_chunk_index = destructed[i] / global_var.bit_num
			excited.occup[Det.chunk_num-dest_chunk_index-1] &= ~2**(destructed[i] - dest_chunk_index * global_var.bit_num)
			crea_chunk_index = created[i] / global_var.bit_num
			excited.occup[Det.chunk_num-crea_chunk_index-1] |= 2**(created[i] - crea_chunk_index * global_var.bit_num)
		return excited
	
	def random_excite(self, basis_list):
		ground_index = load_data.get_index(basis_list, self.occup)
		excited_index = ground_index
		while excited_index == ground_index:
			excited_index = random.sample(range(len(basis_list)), 1)[0]
		excited = Det(basis_list[excited_index])
		return excited
	
	def excite(self, basis_list = [], excite_type = 'd'):
		if excite_type == 'd':
			excited = self.double_excite()
			pq_prob = 1.0 / (2*global_var.e_num*(global_var.e_num-1));
			rs_prob = 1.0 / (2*(global_var.spin_orb_num-global_var.e_num)*(global_var.spin_orb_num-global_var.e_num-1));
			select_prob = pq_prob * rs_prob
		elif excite_type == 'r':
			excited = self.random_excite(basis_list)
			select_prob = 1.0 / (len(basis_list)-1)
		return select_prob, excited
	
	@classmethod
	def occup_2_indices(self, occup):
		# Converts from occup to occup_indices. Inverse of indices_2_occup.
		# e.g. from [4,8,2] = [(0100)_dec,(1000)_dec,(0010)_dec] to 0100,1000,0010 = [10,7,1].
		occup_indices = []
		for i in range(len(occup)):
			for j in range(global_var.bit_num-1, -1, -1):
				if (occup[i]>>j) & 1:
					occup_indices.append((Det.chunk_num-i-1) * global_var.bit_num + j)
		return occup_indices
	
	@classmethod
	def indices_2_occup(self, occup_indices):
		# Converts from occup_indices to occup. Inverse of occup_2_indices.
		# e.g. from [10,7,1] = 0100,1000,0010 to [(0100)_dec,(1000)_dec,(0010)_dec] = [4,8,2].
		occup = [0] * Det.chunk_num
		for i in occup_indices:
			chunk_index = Det.chunk_num - i / global_var.bit_num - 1
			power = i % global_var.bit_num
			occup[chunk_index] |= 2**power
		return occup
	
	@classmethod
	def next_config_indices(self, curr_occup_indices, first_run_flag):
		# Computes the next available configuration, in the form of occup_indices.
		# Configurations are arranged in dictionary order.
		# e.g. [1,0]->[2,0]->[2,1]->[3,0]->[3,1]->[3,2]->...->[5,4]
		if first_run_flag and curr_occup_indices[-1] == global_var.spin_orb_num - global_var.e_num:
			return None
		else:
			next_occup_indices = curr_occup_indices[:]
			if len(next_occup_indices) == 1:
				return [next_occup_indices[0] + 1]
			elif next_occup_indices[-1] == next_occup_indices[-2] - 1:
				next_occup_indices = Det.next_config_indices(next_occup_indices[:-1], False) +\
					[global_var.e_num - len(next_occup_indices)]
			else:
				next_occup_indices[-1] += 1
		return next_occup_indices
	
	@classmethod
	def next_config_occup(self, curr_occup):
		# Computes the next available configuration, in the form of occup.
		curr_occup_indices = Det.occup_2_indices(curr_occup)
		next_occup_indices = Det.next_config_indices(curr_occup_indices, True)
		if next_occup_indices == None:
			return None
		else:
			next_occup = Det.indices_2_occup(next_occup_indices)
			return next_occup