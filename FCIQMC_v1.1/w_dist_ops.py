# Provides all operations on a sparse vector whose entries are $Walker$ objects.
# Used to deal with a walker distribution.
# It is implemented by Python $dict$ structure, i.e. a Hash table.
# The keys are Slater determinants, stored in strings.

def w_dist_normalize(w_dist):
	"""Normalizes a sparse vector."""
	
	# Returns True if it does not vanish, otherwise returns False.
	
	import math
	
	norm_sq = 0
	for key in w_dist:
		norm_sq += w_dist[key].value ** 2
	norm = math.sqrt(norm_sq)
	if norm > 1e-15:
		w_dist_scalar_prod(w_dist, 1/norm, True)
		return True
	else:
		return False

def w_dist_dot_prod(w_dist1, w_dist2):
	"""Dot product of $w_dist1$ and $w_dist2$."""
	
	prod = 0
	for key in w_dist1:
		try:
			entry2 = w_dist2[key].value
			entry1 = w_dist1[key].value
			prod += entry1 * entry2
		except KeyError:
			pass
	return prod

def w_dist_scalar_prod(w_dist, multiplier, change_flag):
	"""Scalar product of $w_dist$ with a scalar multiplier."""
	
	import walker
	
	# Modifies $w_dist$ if $change_flag$, else returns the result as a new vector.

	if change_flag:
		if abs(multiplier) < 1e-15:
			w_dist = {}
		else:
			for key in w_dist:
				w_dist[key].value *= multiplier
		return None
	else:
		prod = {}
		if abs(multiplier) > 1e-15:
			for key in w_dist:
				prod[key] = walker.Walker(w_dist[key].value * multiplier)
		return prod

def w_dist_add(w_dist1, w_dist2, change_flag):
	"""Entry-wise summation of $w_dist1$ and $w_dist2$."""
	
	import walker
	
	# Modifies $w_dist1$ if $change_flag$, else returns the result as a new vector.
	
	if change_flag:
		for key in w_dist1:
			try:
				w_dist1[key].value += w_dist2[key].value
			except KeyError:
				pass
		for key in w_dist2:
			if not w_dist1.has_key(key):
				w_dist1[key] = walker.Walker(w_dist2[key].value, w_dist2[key].flag)
		w_dist1_keys = w_dist1.keys()
		for key in w_dist1_keys:
			if w_dist1[key].value == 0:
				del w_dist1[key]
	else:
		w_dist_sum = {}
		for key in w_dist1:
			w_dist_sum[key] = walker.Walker(w_dist1[key].value, w_dist1[key].flag)
		w_dist_add(w_dist_sum, w_dist2, True)
		return w_dist_sum

def w_dist_subtract(w_dist1, w_dist2, change_flag):
	"""Entry-wise subtraction of $w_dist1$ and $w_dist2$."""
	
	import walker
	
	# Modifies $w_dist1$ if $change_flag$, else returns the result as a new vector.
	
	if change_flag:
		for key in w_dist1:
			try:
				w_dist1[key].value -= w_dist2[key].value
			except KeyError:
				pass
		for key in w_dist2:
			if not w_dist1.has_key(key):
				w_dist1[key] = walker.Walker(-w_dist2[key].value)
		w_dist1_keys = w_dist1.keys()
		for key in w_dist1_keys:
			if w_dist1[key].value == 0:
				del w_dist1[key]
	else:
		w_dist_diff = {}
		for key in w_dist1:
			w_dist_diff[key] = walker.Walker(w_dist1[key].value)
		w_dist_subtract(w_dist_diff, w_dist2, True)
		return w_dist_diff

def w_dist_2_vec(w_dist, normalize_flag = True):
	"""Converts from $w_dist$ to $vec$."""
	
	import vec_ops
	
	vec = {}
	# Only extracts $value$ from each $Walker$ object.
	for key in w_dist:
		vec[key] = w_dist[key].value
	# Normalization.
	if normalize_flag:
		vec_ops.vec_normalize(vec)
	return vec

def w_dist_corr_by_proj(h, w_dist_p):
	"""Calculates correlation energy by projection."""
	
	# E_corr = sum_j(<D_j|H|D_0>*(<N_j>/<N_0>)) - E_0
	corr = -h.ref_energy
	w_dist = h.matrix[h.ref_index]
	for key in w_dist:
		proj = w_dist[key]
		try:
			N0 = w_dist_p[h.ref_key].value
			Nj = w_dist_p[key].value
			corr += Nj/float(N0) * proj
		except KeyError:
			pass
	return corr

def w_dist_spawn(h, w_dist_p, w_dist_c, vec_n_init_c_counts, key_p, tau, is_n_init_flag):
	"""Spawning of all ps on the determinant indexed by $key$."""
	
	import key_ops
	import random
	import walker
	
	# Excitation.
	key_c, prob_rate = key_ops.excited_key_and_prob(h, key_p)
	
	# Spawning.
	if key_c != None:
		p_sign = (w_dist_p[key_p].value > 0)
		p_u_num = abs(w_dist_p[key_p].value)
		count = 0
		for count in range(p_u_num):
			prob = tau * prob_rate
			if prob < 0:
				c_sign = p_sign
				prob = -prob
			else:
				c_sign = not p_sign
			c_num = int(prob)
			# At least we have $c_num$ cren.
			prob -= c_num
			rand_value = random.random()
			if rand_value < prob:
				# One more c.
				c_num += 1
			if c_num != 0:
				if is_n_init_flag:
					try:
						vec_n_init_c_counts[key_c] += c_num
					except KeyError:
						vec_n_init_c_counts[key_c] = c_num
				if not c_sign:
					c_num = -c_num
				w_dist_this_c = {key_c: walker.Walker(c_num)}
				w_dist_add(w_dist_c, w_dist_this_c, True)

def w_dist_tell_init(w_dist):
	w_dist_n_init = {}
	w_dist_init = {}
	for key in w_dist:
		if w_dist[key].flag:
			w_dist_init[key] = w_dist[key]
		else:
			w_dist_n_init[key] = w_dist[key]
	return w_dist_n_init, w_dist_init

def w_dist_c_survive(w_dist_p, w_dist_n_init_c, vec_n_init_c_counts):
	n_init_c_keys = w_dist_n_init_c.keys()
	for key_c in n_init_c_keys:
		if vec_n_init_c_counts[key_c] == 1:
			if not w_dist_p.has_key(key_c):
				# Child cannot survive because:
				# (1) It is spawned by a non-initiator.
				# (2) Its determinant is not occupied by a parent.
				# (3) It is alone on its determinant.
				del w_dist_n_init_c[key_c]

def w_dist_new_init(w_dist, crit_num):
	for key in w_dist:
		if w_dist[key].flag and abs(w_dist[key].value) < crit_num:
			w_dist[key].flag = False
		elif not w_dist[key].flag and abs(w_dist[key].value) >= crit_num:
			w_dist[key].flag = True

def w_dist_die(h, w_dist_p, key_p, tau, shift):
	"""Dying/cloning of all ps on the determinant indexed by $key$."""
	
	import key_ops
	import random
	
	prob_rate = key_ops.dead_prob(h, key_p)
	prob = tau * (prob_rate - shift)
	p_sign = (w_dist_p[key_p].value > 0)
	p_u_num = abs(w_dist_p[key_p].value)
	count = 0
	if prob < 0:
		# Cloning.
		for count in range(p_u_num):
			prob = -tau * (prob_rate - shift)
			clone_num = int(prob)
			# At least we have $clone_num$ clones.
			prob -= clone_num
			rand_value = random.random()
			if rand_value < prob:
				# One more cloning event.
				clone_num += 1
			if clone_num != 0:
				if not p_sign:
					clone_num = -clone_num
				w_dist_p[key_p].value += clone_num
	else:
		# Dying.
		for count in range(p_u_num):
			rand_value = random.random()
			if rand_value < prob:
				# Dead.
				if p_sign:
					w_dist_p[key_p].value -= 1
				else:
					w_dist_p[key_p].value += 1
	if w_dist_p[key_p].value == 0:
		del w_dist_p[key_p]

def w_dist_annih(w_dist_p, w_dist_c):
	"""Overall annihilation step."""
	
	w_dist_add(w_dist_p, w_dist_c, True)

def w_dist_count_u_num(w_dist):
	"""Count the unsigned number of walkers."""
	
	total_u_num = 0
	for key in w_dist:
		total_u_num += int(abs(w_dist[key].value))
	return total_u_num

def w_dist_single_step(h, w_dist_p, tau, shift, init_crit_w_num):
	"""The whole process of spawning, dying and annihilation."""
	
	vec_n_init_c_counts = {}
	w_dist_c = {}
	
	w_dist_n_init_p, w_dist_init_p = w_dist_tell_init(w_dist_p)
	
	# Spawning of non-initiators.
	n_init_p_keys = w_dist_n_init_p.keys()
	for key_p in n_init_p_keys:
		w_dist_spawn(h, w_dist_n_init_p, w_dist_c, vec_n_init_c_counts, key_p, tau, True)
	# Children of non_initiators need to pass a survival test.
	w_dist_c_survive(w_dist_p, w_dist_c, vec_n_init_c_counts)
	
	# Spawning of initiators.
	init_p_keys = w_dist_init_p.keys()
	for key_p in init_p_keys:
		w_dist_spawn(h, w_dist_init_p, w_dist_c, {}, key_p, tau, False)
	
	p_keys = w_dist_p.keys()
	for key_p in p_keys:
		w_dist_die(h, w_dist_p, key_p, tau, shift)
	w_dist_annih(w_dist_p, w_dist_c)
	
	w_dist_new_init(w_dist_p, init_crit_w_num)