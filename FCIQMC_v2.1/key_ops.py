import integrals
import random
import gen_integrals
# Provides all operations on a key.
# Keys are tuples like (4,8,2), i.e. decimal version of bit representations of orbital occupation.
# The example means a configuration [0100,1000,0010].
# Keys are used to refer to determinants in a dictionary.
# Orbs are tuples that carry the same information as keys. They tell which orbitals are taken.
# The standard ordering of an orbs tuple is from largest/leftmost to smallest/rightmost.

def key_2_orbs(key):
	"""Converts from key to orbs. Inverse of orbs_2_key."""
	
	# e.g. from (4,8,2) = ((0100)_dec,(1000)_dec,(0010)_dec) to 0100,1000,0010 = (10,7,1).
	
	bit_num = integrals.para_list[0]
	chunk_num = integrals.para_list[1]
	orbs = ()
	for i in range(len(key)):
		for j in range(bit_num-1, -1, -1):
			if (key[i]>>j) & 1:
				orbs += (((chunk_num-i-1)*bit_num + j),)
	return orbs
	
def orbs_2_key(orbs):
	"""Converts from orbs to key. Inverse of key_2_orbs."""
	
	# e.g. from (10,7,1) = 0100,1000,0010 to ((0100)_dec,(1000)_dec,(0010)_dec) = (4,8,2).
	
	bit_num = integrals.para_list[0]
	chunk_num = integrals.para_list[1]
	key_list = [0] * chunk_num
	for i in orbs:
		chunk_index = chunk_num - i/bit_num - 1
		power = i % bit_num
		key_list[chunk_index] |= 2**power
	return tuple(key_list)

def double_excite(key_gnd):
	"""Generate a double excitation of the determinant labled by key_gnd."""
	
	chunk_num = integrals.para_list[1]
	orb_num = integrals.para_list[2]
	e_num = integrals.para_list[3]
	
	orbs_gnd = key_2_orbs(key_gnd)

	alpha_range = range(orb_num, 2*orb_num)
	beta_range = range(0, orb_num)
	
	# Randomly choose two electrons (m,n), such that m>n.
	pos_m, pos_n = random.sample(range(e_num), 2)
	if pos_m < pos_n:
		pos_mn = (pos_m, pos_n) # This ordering is correct, for a smaller tuple index corresponds to a larger number.
	else:
		pos_mn = (pos_n, pos_m)
	
	orb_mn = (orbs_gnd[pos_mn[0]], orbs_gnd[pos_mn[1]])
	if orb_mn[1] >= orb_num:
		# Both alpha spins.
		p_gen_double = 4.0 / (e_num*(e_num-1)*(orb_num-0.5*e_num)*(orb_num-0.5*e_num-1))
		p_range = alpha_range
		q_range = alpha_range
	elif orb_mn[0] < orb_num:
		# Both beta spins.
		p_gen_double = 4.0 / (e_num*(e_num-1)*(orb_num-0.5*e_num)*(orb_num-0.5*e_num-1))
		p_range = beta_range
		q_range = beta_range
	else:
		# One alpha spin and one beta spin.
		p_gen_double = 2.0 / (e_num*(e_num-1)*(orb_num-0.5*e_num)*(orb_num-0.5*e_num))
		p_range = alpha_range
		q_range = beta_range
	
	# Generate double excitation (p,q), such that p>q.
	while True:
		orb_p = random.sample(p_range, 1)[0]
		if orb_p not in orbs_gnd:
			break
	while True:
		orb_q = random.sample(q_range, 1)[0]
		if (orb_q not in orbs_gnd) and (orb_q != orb_p):
			if orb_q > orb_p:
				orb_pq = (orb_q, orb_p)
			else:
				orb_pq = (orb_p, orb_q)
			break
	
	# m is excited to p, while n is excited to q.
	orbs_exc_list = list(orbs_gnd)
	orbs_exc_list[pos_mn[0]] = orb_pq[0]
	orbs_exc_list[pos_mn[1]] = orb_pq[1]
	orbs_exc = tuple(orbs_exc_list) # The excited orbs tuple.
	
	sign_exc = cal_sign(orbs_exc, pos_mn)
	key_exc = orbs_2_key(orbs_exc)
	orbs_diff = orb_mn + orb_pq # A 4-tuple (m,n,p,q), such that m>n and p>q.
	
	return orbs_gnd, key_exc, sign_exc, p_gen_double, orbs_diff

def single_excite(key_gnd):
	"""Generate a single excitation of the determinant labled by key_gnd."""
	
	chunk_num = integrals.para_list[1]
	orb_num = integrals.para_list[2]
	e_num = integrals.para_list[3]
	
	orbs_gnd = key_2_orbs(key_gnd)
	
	alpha_range = range(orb_num, 2*orb_num)
	beta_range = range(0, orb_num)
	
	pos_m = random.randint(0, e_num-1)
	orb_m = orbs_gnd[pos_m]
	if orb_m >= orb_num:
		# Alpha spin.
		p_range = alpha_range
	else:
		# Beta spin.
		p_range = beta_range
	p_gen_single = 1.0 / (e_num*(orb_num-0.5*e_num))
	
	# Generate single excitation.
	while True:
		orb_p = random.sample(p_range, 1)[0]
		if orb_p not in orbs_gnd:
			break
	
	# m is excited to p.
	orbs_exc_list = list(orbs_gnd)
	orbs_exc_list[pos_m] = orb_p
	orbs_exc = tuple(orbs_exc_list)
	
	sign_exc = cal_sign(orbs_exc, (pos_m,))
	key_exc = orbs_2_key(orbs_exc)
	orbs_diff = (orb_m, orb_p)
	
	return orbs_gnd, key_exc, sign_exc, p_gen_single, orbs_diff

def cal_sign(orbs, pos_wrong):
	"""Calculate the sign of an orbs tuple."""
	
	# pos_wrong indicates which indices in orbs tuple are not ordered in the standard (largest-to-smallest) way.
	# This method returns the parity of the orbs tuple. In other words, how many swaps needed to order it.
	# For example, orbs = (5,8,4,2,3,0) and pos_wrong = (1,4) means the 8 and 3 are wrongly posed.
	# The method returns True, which means an even parity.
	
	# Note: how do I get the parity of a tuple?
	# A standard ordered tuple is sorted in reverse order, e.g. Dg = (7,5,3,1).
	# (1) Single excitation:
	# An unsorted single excitation of the tuple looks like De' = (7,5,3,4), where the excitation is 1 --> 4.
	# To get the parity of De', count (a) the number of indices before 4 that are less than 4,
	#                                 (b) the number of indices after 4 that are greater than 4.
	# There is one, namely 3<4.
	# This means one swap will make De' ordered. So the parity is (-1)^1 = -1.
	# In other words, De = (7,5,4,3) = -De'.
	# (2) Double excitation:
	# It is not much different. Let De' = (6,5,3,4) be the excited state. Note that (6,4) is reversely ordered.
	# Always sort the two excitation spots in reverse order before the next step.
	# Count (a) the number of indices before 6 that are less than 6,
	#       (b) the number of indices after 6 that are greater than 6,
	#       (c) the number of indices before 4 that are less than 4,
	#       (d) the number of indices after 4 that are greater than 4.
	# There is one, namely 3<4.
	# This means one swap will make De' ordered. So the parity is (-1)^1 = -1.
	# In other words, De = (6,5,4,3) = -De'.
	
	rev_num = 0 # Initiates the reverse number.
	
	if len(pos_wrong) == 1:
		# Single excitation.
		pos_p = pos_wrong[0]
		orb_p = orbs[pos_p]
		for pos in range(pos_p-1, -1, -1):
			if orbs[pos] < orb_p:
				rev_num += 1
			else:
				break
		for i in orbs[pos_p+1 : ]:
			if i > orb_p:
				rev_num += 1
			else:
				break
	elif len(pos_wrong) == 2:
		# Double excitation.
		pos_p, pos_q = pos_wrong
		orb_p = orbs[pos_p]
		orb_q = orbs[pos_q]
		for pos in range(pos_p-1, -1, -1):
			if orbs[pos] < orb_p:
				rev_num += 1
			else:
				break
		for i in orbs[pos_p+1 : ]:
			if i > orb_p:
				rev_num += 1
			elif i != orb_q:
				break
		for pos in range(pos_q-1, -1, -1):
			if orbs[pos] < orb_q:
				rev_num += 1
			elif orbs[pos] != orb_p:
				break
		for i in orbs[pos_q+1 : ]:
			if i > orb_q:
				rev_num += 1
			else:
				break
	return (rev_num % 2 == 0)

def difference(key_gnd, key_exc):
	"""Calculate the difference between key_gnd and key_exc."""
	
	# Returns the sign of excitation as well as the difference for matrix element calculation.
	# What is the tuple diff:
	# (1) for an input of a ground state |Psi> and an excited state |Psi_m^p>, diff = (m,p).
	# (2) for an input of a ground state |Psi> and an excited state |Psi_mn^pq>, diff = (m,n,p,q), such that m>n and p>q.
	bit_num = integrals.para_list[0]
	chunk_num = integrals.para_list[1]
	
	orb_mn = ()
	orb_pq = ()
	pos_wrong = ()
	orbs_gnd = key_2_orbs(key_gnd)
	diff = [0] * chunk_num
	count = 0
	
	for i in range(chunk_num):
		diff[i] = key_gnd[i] ^ key_exc[i]
		for j in range(bit_num-1, -1, -1):
			if (diff[i]>>j) & 1:
				orb = (chunk_num-i-1) * bit_num + j
				count += 0.5
				if count <= 2:
					if (key_gnd[i]>>j) & 1:
						pos_wrong = pos_wrong + (orbs_gnd.index(orb),)
						orb_mn = orb_mn + (orb,)
					else:
						orb_pq = orb_pq + (orb,)
				else:
					# Differ by more than double excitation.
					return None, None, None
	orbs_diff = orb_mn + orb_pq
	
	orbs_exc_list = list(orbs_gnd)
	for i in range(int(count)):
		orbs_exc_list[pos_wrong[i]] = orb_pq[i]
	orbs_exc = tuple(orbs_exc_list)
	sign_exc = cal_sign(orbs_exc, pos_wrong)
	
	return orbs_gnd, sign_exc, orbs_diff