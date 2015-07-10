import random_vec
import global_var

# Write/read a Hamiltonian into/from a file.
def write_vec(h_file, vec):
	for i in vec:
		h_file.write(str(i) + ' ')
	h_file.write('\r\n')

def line_2_vec(line):
	vec = []
	strvec = line.split(' ')
	strvec.pop()
	for string in strvec:
		vec.append(float(string))
	return vec

def write_h_file(h, h_file):
	parameters = [global_var.spin_orb_num/2, global_var.e_num, global_var.dim]
	write_vec(h_file, parameters)
	h_file.write('\r\n')
	for vec in h.eig_vecs:
		write_vec(h_file, vec.vec)
	h_file.write('\r\n')
	write_vec(h_file, h.eig_vals)
      
def read_h_file(h_file):
	line = h_file.readline()
	vec = line_2_vec(line)
	global_var.spin_orb_num = 2 * int(vec[0])
	global_var.e_num = int(vec[1])
	global_var.dim = int(vec[2])
	line = h_file.readline()
	
	eig_vecs = []
	line = h_file.readline()
	while len(line) > 2:
		vec = line_2_vec(line)
		eig_vecs.append(random_vec.RandomVec(True, False, vec))
		line = h_file.readline()
	line = h_file.readline()
	eig_vals = line_2_vec(line)
	return eig_vecs, eig_vals