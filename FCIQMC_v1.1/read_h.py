# Provides all operations required to read a Hamiltonian from a file.

def read_outof_file(file_name, para_list, eig_vecs, eig_vals):
	"""Reads a Hamiltonian from a file."""
	
	h_file = open(file_name, 'r')
	
	# Reads the parameters.
	line = h_file.readline()
	str_para_list = line.split('  ')
	del str_para_list[-1]
	for i in range(len(str_para_list)):
		para_list.append(int(str_para_list[i]))
	
	# Reads the eigenvectors.
	for i in range(para_list[4]):
		vec = {}
		line = h_file.readline()
		str_list = line.split('   ')
		del str_list[-1]
		for string in str_list:
			str_sublist = string.split('  ')
			key = str_sublist[0]
			value = float(str_sublist[1])
			vec[key] = value
		eig_vecs.append(vec)
	
	# Reads the eigenvalues
	line = h_file.readline()
	str_list = line.split('  ')
	for i in range(para_list[4]):
		eig_vals.append(float(str_list[i]))
	
	h_file.close()