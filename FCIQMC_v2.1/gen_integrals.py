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