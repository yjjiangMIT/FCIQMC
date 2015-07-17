import key_ops
import integrals
import ctrl_panel

def gen_key_list(ref_key):
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
	key_list = classify_key_list(key_list, ref_key)
	return key_list

def classify_key_list(key_list, ref_key):
	key_list = sort_key_list(key_list)
	keys_S = []
	keys_D = []
	keys_T = []
	keys_Q = []
	keys_P = []
	keys_H = []
	for key in key_list:
		if key != ref_key:
			orbs_i, sign, orbs_diff, count = key_ops.difference(ref_key, key)
			if count == 1:
				keys_S.append(key)
			elif count == 2:
				keys_D.append(key)
			elif count == 3:
				keys_T.append(key)
			elif count == 4:
				keys_Q.append(key)
			elif count == 5:
				keys_P.append(key)
			elif count == 6:
				keys_H.append(key)
	key_list = [ref_key] + keys_S + keys_D + keys_T + keys_Q + keys_P + keys_H
	
	return key_list

def sort_key_list(key_list):
	ke_pair_list = []
	for key in key_list:
		energy = integrals.mat_element(key, key)
		ke_pair_list.append((energy, key))
	ke_pair_list.sort()
	sorted_key_list = []
	for ke in ke_pair_list:
		sorted_key_list.append(ke[1])
	return sorted_key_list

def write_mat_2_file(file_name = 'Hamiltonian_new.m'):
	"""Writes the full Hamiltonian as well as some MATLAB codes into a .m file."""
	
	f = open(file_name, 'w')
	f.write('data = [')
	for i in range(len(key_list)):
		key_i = ctrl_panel.key_list[i]
		for j in range(len(key_list)):
			key_j = ctrl_panel.key_list[j]
			orbs_i, sign, orbs_diff, count = key_ops.difference(key_i, key_j)
			if orbs_i != None:
				entry = integrals.sandwich(orbs_i, orbs_diff)
				if abs(entry) < 1e-8:
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
