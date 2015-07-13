import det
import key_ops
import det_ops
import black_box

def gen_matrix():
	"""Creates a MATLAB file of the full Hamiltonian as a sparse matrix."""
	
	dets_p = {}
	dets_c = {}
	
	spatial_orbs_list = []
	key_list = []
	
	for i in range(0, 8):
		for j in range(i+1, 8):
			for k in range(j+1, 8):
				for l in range(k+1, 8):
					for m in range(l+1, 8):
						spatial_orbs_list.append([m, l, k, j, i])
	for sp_orbs_i in spatial_orbs_list:
		for sp_orbs_j in spatial_orbs_list:
			sp_orbs_temp = sp_orbs_j[:]
			for k in range(5):
				sp_orbs_temp[k] += 8
			orbs = tuple(sp_orbs_temp + sp_orbs_i)
			key = key_ops.orbs_2_key(orbs)
			key_list.append(key)
	for key in key_list:
		dets_c[key] = det.Det(1, True)
	
	det_ops.merge(dets_p, dets_c)
	
	f = open('Hamiltonian_N2.m', 'w')
	f.write('data = [')
	for i in range(len(key_list)):
		key_i = key_list[i]
		for j in range(len(key_list)):
			key_j = key_list[j]
			orbs_i, sign, orbs_diff = key_ops.difference(key_i, key_j)
			if orbs_i != None:
				entry = black_box.sandwich(orbs_i, orbs_diff)
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
	
	return dets_p