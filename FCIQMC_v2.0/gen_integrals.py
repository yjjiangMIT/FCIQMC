import random
import math

def gen_spatial_one_body_ints(para_list):
	"""Generate random numbers for all one-body integrals."""
	
	# $black_box_one_body$ is a dictionary, whose keys are tuples of the two spatial orbitals 
	# involved in the integrals, and whose values are randomly generated matrix elements.
	# e.g. key = (i,j), value = (i|h|j).
	
	orb_num = para_list[2]
	black_box_one_body = {}
	
	for orb_i in range(orb_num):
		for orb_j in range(orb_i, orb_num):
			key = (orb_i, orb_j)
			rand_val = random.random()
			if orb_i == orb_j:
				if rand_val < 0.02:
					black_box_one_body[key] = random.random() * 10
				else:
					black_box_one_body[key] = random.random()
			else:
				if rand_val < 0.01:
					black_box_one_body[key] = random.random()*2 - 1
	return black_box_one_body

def gen_spatial_two_body_ints(para_list):
	"""Generate random numbers for all two-body integrals."""
	
	# $black_box_two_body$ is a dictionary, whose keys are tuples of the four spatial orbitals 
	# involved in the integrals, and whose values are randomly generated matrix elements.
	# e.g. key = (i,j,k,l), value = (ij|kl).
	
	orb_num = para_list[2]
	black_box_two_body = {}
	
	for orb_i in range(orb_num):
		for orb_j in range(orb_i, orb_num):
			orb_ij = (orb_i, orb_j)
			for orb_k in range(orb_num):
				for orb_l in range(orb_k, orb_num):
					orb_kl = (orb_k, orb_l)
					rand_val = random.random()
					key = orb_ij + orb_kl
					if cmp(orb_ij, orb_kl) < 0:
						if rand_val < 0.01:
							black_box_two_body[key] = random.random()*2 - 1
					elif cmp(orb_ij, orb_kl) == 0:
						if rand_val < 0.02:
							black_box_two_body[key] = random.random() * 10
						else:
							black_box_two_body[key] = random.random()
	return black_box_two_body

def disp_ints(black_box):
	"""Print $black_box$ to screen."""
	
	for key in black_box:
		print key,
		print black_box[key]

def write_into_file(bit_num, orb_num, e_num, file_name = 'integrals'):
	"""Writes integrals into a file."""
	
	chunk_num = int(math.ceil(2.0 * orb_num / bit_num))
	dim = (math.factorial(orb_num) / math.factorial(e_num/2) / math.factorial(orb_num-e_num/2)) ** 2
	
	para_list = [bit_num, chunk_num, orb_num, e_num, dim]
	black_box_one_body = gen_spatial_one_body_ints(para_list)
	black_box_two_body = gen_spatial_two_body_ints(para_list)
	
	data_file = open('Documents/FCIQMC/FCIQMC_v2.0/' + file_name + '.txt', 'w')
	
	# Writes the parameters.
	for i in para_list:
		data_file.write(str(i) + ' ')
	data_file.write('\r\n')
	
	# Writes the data.
	for key in black_box_one_body:
		data_file.write(str(key)[1:-1] + ': ')
		data_file.write(str(black_box_one_body[key]) + '; ')
	data_file.write('\r\n')
	for key in black_box_two_body:
		data_file.write(str(key)[1:-1] + ': ')
		data_file.write(str(black_box_two_body[key]) + '; ')
	data_file.write('\r\n')
	
	data_file.close()

def read_outof_file(file_name = 'integrals'):
	"""Reads integrals from a file."""
	
	data_file = open('Documents/FCIQMC/FCIQMC_v2.0/' + file_name + '.txt', 'r')
	
	# Reads the parameters.
	line = data_file.readline()
	str_para_list = line.split(' ')
	del str_para_list[-1]
	para_list = []
	for i in str_para_list:
		para_list.append(int(i))
	
	# Reads one-body integrals.
	black_box_one_body = {}
	
	line = data_file.readline()
	str_one_body_ints = line.split('; ')
	del str_one_body_ints[-1]
	for i in str_one_body_ints:
		temp = i.split(': ')
		key = (int(temp[0][0]), int(temp[0][3]))
		value = float(temp[1])
		black_box_one_body[key] = value
	
	# Reads two-body integrals.
	black_box_two_body = {}
	
	line = data_file.readline()
	str_two_body_ints = line.split('; ')
	del str_two_body_ints[-1]
	for i in str_two_body_ints:
		temp = i.split(': ')
		key = (int(temp[0][0]), int(temp[0][3]), int(temp[0][6]), int(temp[0][9]))
		value = float(temp[1])
		black_box_two_body[key] = value
	
	data_file.close()
	return para_list, black_box_one_body, black_box_two_body