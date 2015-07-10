# Provides all operations required to write a Hamiltonian into a file.

def write_into_file(file_name, para_list, sparsity):
	"""Writes a Hamiltonian into a file."""
	
	h_file = open(file_name, 'w')
	
	# Writes some parameters.
	for i in para_list:
		h_file.write(str(i) + '  ')
	h_file.write('\r\n')
	eig_vecs = create_eig_vecs(para_list, sparsity)
	
	# Writes the eigenvectors.
	for vec in eig_vecs:
		for key in vec: 
			h_file.write(key + '  ')
			h_file.write(str(vec[key]) + '   ')
		h_file.write('\r\n')
	eig_vals = create_eig_vals(para_list)
	
	# Writes the eigenvalues.
	for val in eig_vals:
		h_file.write(str(val) + '  ')
	h_file.write('\r\n')
	h_file.close()

def create_eig_vals(para_list):
	"""Creates random eigenvalues."""
	
	import random
	
	eig_vals = [0] # The smallest eigenvalue vanishes.
	dim = para_list[4]
	for i in range(dim - 1):
		eig_vals.append(random.uniform(1, 2))
	eig_vals.sort()
	return eig_vals

def create_eig_vecs(para_list, sparsity):
	"""Creates random sparse eigenvectors that are orthonormal."""
	
	import key_ops
	import vec_ops
	
	e_num = para_list[3]
	dim = para_list[4]
	indices = range(e_num-1, -1, -1) 
	eig_vecs = []
	for i in range(dim):
		while True:
			occup = key_ops.indices_2_occup(indices, para_list)
			vec = {}
			while occup != None:
				key = str(occup)
				key = key[1:-1]
				value = get_value(sparsity)
				if value != 0:
					vec[key] = value
				occup = key_ops.next_config_occup(occup, para_list)
			create_suc_flag = vec_ops.vec_schmidt(vec, eig_vecs) # Orthonormalize.
			if create_suc_flag:
				break
	return eig_vecs

def get_value(sparsity):
	"""Gets a random value."""
	
	import random
	
	if random.random() > sparsity:
		# Zero entry.
		return 0
	else:
		return random.uniform(-1, 1)