from random_h import *
from det import *

def get_index(basis, occup):
	if len(basis) == 0:
		return None
	else:
		middle = len(basis) / 2
		comparison = cmp(occup, basis[middle])
		if comparison < 0:
			index = get_index(basis[:middle], occup)
		elif comparison > 0:
			index = get_index(basis[middle+1:], occup)
			if index != None:
				index += middle + 1
		else:
			index = middle
		return index

def get_entry(h, occup1, occup2):
	index1 = get_index(h.basis_list, occup1)
	index2 = get_index(h.basis_list, occup2)
	if index1 == None or index2 == None:
		return None
	else:
		return h.matrix[index1].vec[index2]