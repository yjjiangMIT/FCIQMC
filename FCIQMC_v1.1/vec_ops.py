# Provides all operations on a sparse vector whose entries are float-point numbers or integers.
# Used to construct a Hamiltonian.
# It is implemented by Python $dict$ structure, i.e. a Hash table.
# The keys are Slater determinants, stored in strings.

def vec_normalize(vec):
	"""Normalizes a sparse vector."""
	
	# Returns True if it does not vanish, otherwise returns False.
	
	import math
	
	norm_sq = 0
	for key in vec:
		norm_sq += vec[key] ** 2
	norm = math.sqrt(norm_sq)
	if norm > 1e-15:
		vec_scalar_prod(vec, 1/norm, True)
		return True
	else:
		return False

def vec_dot_prod(vec1, vec2):
	"""Dot product of $vec1$ and $vec2$."""
	
	prod = 0
	for key in vec1:
		try:
			entry2 = vec2[key]
			entry1 = vec1[key]
			prod += entry1 * entry2
		except KeyError:
			pass
	return prod

def vec_scalar_prod(vec, multiplier, change_flag):
	"""Scalar product of $vec$ with a scalar multiplier."""
	
	# Modifies $vec$ if $change_flag$, else returns the result as a new vector.
	
	if change_flag:
		if abs(multiplier) < 1e-15:
			vec = {}
		else:
			for key in vec:
				vec[key] *= multiplier
		return None
	else:
		prod = {}
		if abs(multiplier) > 1e-15:
			for key in vec:
				prod[key] = vec[key] * multiplier
		return prod

def vec_add(vec1, vec2, change_flag):
	"""Entry-wise summation of $vec1$ and $vec2$."""
	
	# Modifies $vec1$ if $change_flag$, else returns the result as a new vector.
	
	if change_flag:
		for key in vec1:
			try:
				vec1[key] += vec2[key]
			except KeyError:
				pass
		for key in vec2:
			if not vec1.has_key(key):
				vec1[key] = vec2[key]
		vec1_keys = vec1.keys()
		for key in vec1_keys:
			if vec1[key] == 0:
				del vec1[key]
	else:
		vec_sum = {}
		for key in vec1:
			vec_sum[key] = vec1[key]
		vec_add(vec_sum, vec2, True)
		return vec_sum

def vec_subtract(vec1, vec2, change_flag):
	"""Entry-wise subtraction of $vec1$ and $vec2$."""
	
	# Modifies $vec1$ if $change_flag$, else returns the result as a new vector.
	
	if change_flag:
		for key in vec1:
			try:
				vec1[key] -= vec2[key]
			except KeyError:
				pass
		for key in vec2:
			if not vec1.has_key(key):
				vec1[key] = -vec2[key]
		vec1_keys = vec1.keys()
		for key in vec1_keys:
			if vec1[key] == 0:
				del vec1[key]
	else:
		vec_diff = {}
		for key in vec1:
			vec_diff[key] = vec1[key]
		vec_subtract(vec_diff, vec2, True)
		return vec_diff

def vec_schmidt(vec, orth_vec_list):
	"""Schmidt orthogonalization."""
	
	# Makes $vec$ orthonormal to $orth_vec_list$, a list of vectors, and appends $vec$ to it.
	
	for other_vec in orth_vec_list:
		overlap = vec_dot_prod(vec, other_vec)
		projection = vec_scalar_prod(other_vec, overlap, False)
		vec_subtract(vec, projection, True)
	nonvanish_flag = vec_normalize(vec)
	if nonvanish_flag:
                orth_vec_list.append(vec)
        return nonvanish_flag

def vec_sparse_2_dense(key_list, vec_dist, dim):
	"""Convertes a sparse vector to a dense one."""
	
	import key_ops
	
	dense_vec = [0] * dim
	for key in vec_dist:
		index = key_ops.key_2_index(key_list, key)
		dense_vec[index] = vec_dist[key]
	return dense_vec

def vec_dom_entries(vec, cutoff, normalize_flag = True):
	"""Picks out dominant entries."""
	
	key_list = vec.keys()
	for key in key_list:
		if abs(vec[key]) < cutoff:
			del vec[key]
	if normalize_flag:
		vec_normalize(vec)