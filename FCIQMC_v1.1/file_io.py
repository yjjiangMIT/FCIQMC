# Summarizes reading and writing of files.

def write_hamiltonian(orb_num, e_num, bit_num, sparsity):
	"""Writes a Hamiltonian into a file."""
	
	import write_h
	import math
	
	spin_orb_num = 2 * orb_num
	# Number of chunks required to store a configuration.
	chunk_num = int(math.ceil(spin_orb_num/float(bit_num)))
	# Dimension of the FCI space.
	dim = math.factorial(spin_orb_num) / math.factorial(e_num) / math.factorial(spin_orb_num - e_num)
	
	para_list = [bit_num, chunk_num, orb_num, e_num, dim]
	file_name = 'Documents/FCIQMC_sparse/h_file_' + str(orb_num) + '_' + str(e_num) + '.txt'
	
	write_h.write_into_file(file_name, para_list, sparsity)
	
def read_hamiltonian(orb_num, e_num):
	"""Reads a Hamiltonian from a file."""
	
	import hamiltonian
	
	file_name = 'Documents/FCIQMC_sparse/h_file_' + str(orb_num) + '_' + str(e_num) + '.txt'
	h = hamiltonian.Hamiltonian(file_name)
	return h