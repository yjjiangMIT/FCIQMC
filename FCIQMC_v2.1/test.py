import det
import key_ops
import det_ops
import integrals

def gen_key_list():
	"""Creates keys of all the determinants."""
	
	orb_num = integrals.para_list[2]
	e_num = integrals.para_list[3]
	
	dets_p = {}
	dets_c = {}
	spatial_orbs_list = []
	key_list = []
	
	# Creates all spatial orbital occupations.
	add_digit(0, orb_num, e_num/2, [], spatial_orbs_list)
	
	for sp_orbs_i in spatial_orbs_list:
		for sp_orbs_j in spatial_orbs_list:
			sym = 0
			for orb_i in sp_orbs_i:
				sym ^= integrals.sym_tuple[orb_i]
			for orb_j in sp_orbs_j:
				sym ^= integrals.sym_tuple[orb_j]
			if sym == 0:
				# Only picks out those determinants with the correct symmetry.
				sp_orbs_temp = sp_orbs_j[:]
				for k in range(e_num/2):
					sp_orbs_temp[k] += orb_num
				orbs = tuple(sp_orbs_temp + sp_orbs_i) # Creates all spin orbital occupations.
				key = key_ops.orbs_2_key(orbs)
				key_list.append(key)
	
	return key_list

def write_mat_2_file(file_name = 'Hamiltonian_N2.m'):
	"""Writes the full Hamiltonian as well as some MATLAB codes into a .m file."""
	
	key_list = gen_key_list()
	
	f = open(file_name, 'w')
	f.write('data = [')
	for i in range(len(key_list)):
		key_i = key_list[i]
		for j in range(len(key_list)):
			key_j = key_list[j]
			orbs_i, sign, orbs_diff = key_ops.difference(key_i, key_j)
			if orbs_i != None:
				entry = integrals.sandwich(orbs_i, orbs_diff)
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

def add_digit(range_start, range_end, max_len, list_last_level, list_ensemble):
	"""Recursively build up all the lists of spatial orbital occupations."""
	
	if len(list_last_level) == max_len:
		list_ensemble.append(list_last_level)
	else:
		for i in range(range_start, range_end):
			list_this_level = [i] + list_last_level
			add_digit(i+1, range_end, max_len, list_this_level, list_ensemble)
