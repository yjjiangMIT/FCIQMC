import os
import math

def fcidump(file_name = 'FCIDUMP_00031788'):
	'''Load one- and two-particle integrals from FCIDUMP file.'''
	
	# Opens FCIDUMP file.
	dir_name = '/home/jiang/Documents/FCIQMC/FCIQMC_v2.1/'
	current_dir = os.getcwd()
	dir_name = dir_name[len(current_dir)+1:]
	f = open(dir_name + file_name, 'r')
	
	# Gets number of orbitals, electrons and twice the spin value.
	string = f.readline()
	string = string.split(',')
	orb_num = int(string[0].split(' ')[-1])
	e_num = int(string[1].split(' ')[-1])
	twice_spin = int(string[2].split(' ')[-1])
	
	bit_num = 4
	chunk_num = int(math.ceil(2.0 * orb_num / bit_num))
	dim = (math.factorial(orb_num) / math.factorial(e_num/2) / math.factorial(orb_num-e_num/2)) ** 2
	
	para_list = [bit_num, chunk_num, orb_num, e_num, dim, twice_spin]
	
	# Gets symmetry of the orbitals.
	string = f.readline()
	string = string.split('=')[-1]
	string = string.split(',')
	sym_tuple = ()
	for orb in range(orb_num):
		sym_tuple += ((int(string[orb]) - 1),)
	
	# Reads two useless lines.
	string = f.readline()
	string = f.readline()
	
	# Reads the integrals.
	black_box_one_body = {}
	black_box_two_body = {}
	
	string = f.readline()
	string = string.split('   ')
	string[-1] = string[-1][:-1]
	
	while int(string[3]) != 0:
		# Two-particle integrals.
		key = (int(string[1])-1, int(string[2])-1, int(string[3])-1, int(string[4])-1)
		value = float(string[0])
		black_box_two_body[key] = value
		
		string = f.readline()
		string = string.split('   ')
		string[-1] = string[-1][:-1]
	while int(string[1]) != 0:
		# One-particle integrals.
		key = (int(string[1])-1, int(string[2])-1)
		value = float(string[0])
		black_box_one_body[key] = value
		
		string = f.readline()
		string = string.split('   ')
		string[-1] = string[-1][:-1]
	
	f.close()
	
	return para_list, sym_tuple, black_box_one_body, black_box_two_body

def disp_ints(black_box):
	"""Print $black_box$ to screen."""
	
	for key in black_box:
		print key,
		print black_box[key]