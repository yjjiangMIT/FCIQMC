import read_fcidump
para_list, sym_tuple, black_box_one_body, black_box_two_body = read_fcidump.fcidump()

import gl_consts

def spatial_one_body_int(orb_i, orb_j):
	"""Calculates spatial one-body integral (i|h|j)."""
	
	if orb_i <= orb_j:
		key = (orb_i, orb_j)
	else:
		key = (orb_j, orb_i)
	if key in black_box_one_body:
		return black_box_one_body[key]
	else:
		return 0

def spatial_two_body_int(orb_ij, orb_kl):
	"""Calculates spatial two-body integral (ij|kl)."""
	
	orb_ij = sort_2_tuple(orb_ij)
	orb_kl = sort_2_tuple(orb_kl)
	if cmp(orb_ij, orb_kl) <= 0:
		key = orb_ij + orb_kl
	else:
		key = orb_kl + orb_ij
	if key in black_box_two_body:
		return black_box_two_body[key]
	else:
		return 0

def sort_2_tuple(orb_ij):
	"""Sort $orb_ij$."""
	
	orb_i, orb_j = orb_ij
	if orb_i > orb_j:
		orb_ij = (orb_j, orb_i)
	return orb_ij

def det_one_body_int(orbs_gnd, orbs_diff):
	"""Calculates determinantal one-body integral <Psi|O1|Psi> or <Psi|O1|Psi_m^p>."""
	
	# |Psi_m^p> may differ from the basis by a sign, which must be added outside this method.
	orb_num = gl_consts.para_list[2]
	e_num = gl_consts.para_list[3]
	
	if len(orbs_diff) == 2:
		# Single excitation.
		# |K> = |...mn...>.
		# |L> = |...pn...>.
		# <K|O1|L> = (m|h|p).
		orb_m, orb_p = orbs_diff
		if orb_m >= orb_num:
			# m and p are alpha spins.
			orb_m -= orb_num
			orb_p -= orb_num
		return spatial_one_body_int(orb_m, orb_p)
		
	elif len(orbs_diff) == 0:
		# Diagonal entry.
		# <K|O1|K> = sum_(m_a)^(N/2)(m_a|h|m_a) + sum_(m_b)^(N/2)(m_b|h|m_b).
		result = 0
		for orb in orbs_gnd[ : e_num/2]:
			# m and p are alpha spins.
			orb_m = orb - orb_num
			result += spatial_one_body_int(orb_m, orb_m)
		for orb_m in orbs_gnd[e_num/2 : ]:
			# m and p are beta spins.
			result += spatial_one_body_int(orb_m, orb_m)
		return result
		
	else:
		# If more than single excitation, the one-body integral vanishes.
		return 0

def det_two_body_int(orbs_gnd, orbs_diff):
	"""Calculates determinantal two-body integral <Psi|O2|Psi> or <Psi|O2|Psi_m^p> or <Psi|O2|Psi_mn^pq>."""
	
	# |Psi_m^p> or |Psi_mn^pq> may differ from the basis by a sign, which must be added outside this method.
	orb_num = gl_consts.para_list[2]
	e_num = gl_consts.para_list[3]
	
	if len(orbs_diff) == 4:
		# Double excitation.
		# |K> = |...mn...>.
		# |L> = |...pq...>.
		# <K|O2|L> = (mp|nq) - (mq|np) if m and n have the same spin.
		# <K|O2|L> = (mp|nq)           if m and n have opposite spins.
		orb_m, orb_n, orb_p, orb_q = orbs_diff
		if orb_n >= orb_num:
			# m and n are both alpha spins.
			orb_m -= orb_num
			orb_n -= orb_num
			orb_p -= orb_num
			orb_q -= orb_num
			return spatial_two_body_int((orb_m,orb_p), (orb_n,orb_q)) - \
				spatial_two_body_int((orb_m,orb_q), (orb_n,orb_p))
		elif orb_m < orb_num:
			# m and n are both beta spins.
			return spatial_two_body_int((orb_m,orb_p), (orb_n,orb_q)) - \
				spatial_two_body_int((orb_m,orb_q), (orb_n,orb_p))
		else:
			# m is an alpha spin and n is a beta spin.
			orb_m -= orb_num
			orb_p -= orb_num
			return spatial_two_body_int((orb_m,orb_p), (orb_n,orb_q))
		
	elif len(orbs_diff) == 2:
		# Single excitation.
		# |K> = |...mn...>.
		# |L> = |...pn...>.
		# <K|O2|L> = sum_(n_a)^(N/2)(mp|n_a n_a)-(m n_a|n_a p) + sum_(n_b)^(N/2)(mp|n_b n_b) if m and p are alpha.
		# <K|O2|L> = sum_(n_a)^(N/2)(mp|n_a n_a) + sum_(n_b)^(N/2)(mp|n_b n_b)-(m n_b|n_b p) if m and p are beta.
		result = 0
		orb_m, orb_p = orbs_diff
		if orb_m >= orb_num:
			# Differ in an alpha spin.
			orb_m -= orb_num
			orb_p -= orb_num
			for orb in orbs_gnd[ : e_num/2]:
				# n is an alpha spin.
				orb_n = orb - orb_num
				result += spatial_two_body_int((orb_m,orb_p), (orb_n,orb_n))
				result -= spatial_two_body_int((orb_m,orb_n), (orb_n,orb_p))
			for orb_n in orbs_gnd[e_num/2 : ]:
				# n is a beta spin.
				result += spatial_two_body_int((orb_m,orb_p), (orb_n,orb_n))
		else:
			# Differ in a beta spin.
			for orb in orbs_gnd[ : e_num/2]:
				# n is an alpha spin.
				orb_n = orb - orb_num
				result += spatial_two_body_int((orb_m,orb_p), (orb_n,orb_n))
			for orb_n in orbs_gnd[e_num/2 : ]:
				# n is a beta spin.
				result += spatial_two_body_int((orb_m,orb_p), (orb_n,orb_n))
				result -= spatial_two_body_int((orb_m,orb_n), (orb_n,orb_p))
		return result
		
	elif len(orbs_diff) == 0:
		# Diagonal entry.
		# <K|O2|K> = 1/2*sum_(m_a)^(N/2)sum_(n_a)^(N/2)(m_a m_a|n_a n_a)-(m_a n_a|n_a m_a) + 
		#            1/2*sum_(m_b)^(N/2)sum_(n_b)^(N/2)(m_b m_b|n_b n_b)-(m_b n_b|n_b m_b) + 
		#                sum_(m_a)^(N/2)sum_(n_b)^(N/2)(m_a m_a|n_b n_b)
		result = 0
		for orb_1 in orbs_gnd[ : e_num/2]:
			# m is an alpha spin.
			orb_m = orb_1 - orb_num
			for orb_2 in orbs_gnd[ : e_num/2]:
				# n is an alpha spin.
				orb_n = orb_2 - orb_num
				result += 0.5 * spatial_two_body_int((orb_m,orb_m), (orb_n,orb_n))
				result -= 0.5 * spatial_two_body_int((orb_m,orb_n), (orb_n,orb_m))
			for orb_n in orbs_gnd[e_num/2 : ]:
				# n is a beta spin.
				result += spatial_two_body_int((orb_m,orb_m), (orb_n,orb_n))
		for orb_m in orbs_gnd[e_num/2 : ]:
			# m is a beta spin.
			for orb_n in orbs_gnd[e_num/2 : ]:
				# n is a beta spin.
				result += 0.5 * spatial_two_body_int((orb_m,orb_m), (orb_n,orb_n))
				result -= 0.5 * spatial_two_body_int((orb_m,orb_n), (orb_n,orb_m))
		return result
		
	else:
		# If more than double excitation, the two-body integral vanishes.
		return 0

def sandwich(orbs_gnd, orbs_diff):
	return det_one_body_int(orbs_gnd, orbs_diff) + \
		det_two_body_int(orbs_gnd, orbs_diff)