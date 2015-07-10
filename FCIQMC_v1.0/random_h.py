import random_vec
import random
import det
import math
import global_var
import stored_h

class RandomH:
	# Class of random Hamiltonian: symmetric, positive semi-definite and whose smallest eigenvalue equals zero.
	# Its eigenvectors and eigenvalues are stored as eig_vecs and eig_vals.
	# The Hamiltonian itself is stored as matrix.
	
	def __init__(self, rand_flag = True):
		# Initiates a RandomH object.
		self.eig_vecs = []
		self.eig_vals = [0]
		self.matrix = []
		if rand_flag:
			# Generates random eigenvectors and eigenvalues.
			# The following two lines do not make dominant entries in eigenvectors.
			# for i in range(global_var.dim):
			# 	random_vec.RandomVec().schmidt(self.eig_vecs)
			
			# The first eigenvector have n dominant entries.
			vec = random_vec.RandomVec()
			n = 10
			indices = random.sample(range(global_var.dim), n)
			for i in indices:
				vec.vec[i] = random.sample([-1,1], 1)[0] * random.uniform(2, 3)
			vec.schmidt(self.eig_vecs)
			for i in range(global_var.dim - 1):
				random_vec.RandomVec().schmidt(self.eig_vecs)
			
			# The first eigenvalue is zero.
			for i in range(global_var.dim - 1):
				self.eig_vals.append(random.uniform(1, 2))
			self.eig_vals.sort()
		else:
			# Load data from a file.
			h_file_read = open(global_var.file_name, 'r')
			self.eig_vecs, self.eig_vals = stored_h.read_h_file(h_file_read)
			h_file_read.close()
		for i in range(global_var.dim):
			self.matrix.append(random_vec.RandomVec(False))
		# Computes H=sum(Ei*|vi><vi|)
		for i in range(global_var.dim):
			for j in range(global_var.dim):
				vec = self.eig_vecs[i]
				self.matrix[j].vec_add(vec.scalar_prod(vec.vec[j]*self.eig_vals[i], False), True)
		
		self.ref_energy = self.matrix[0].vec[0]
		# Labels the basis by occup.
		self.basis_occup()
	
	def disp_mat(self):
		# Prints the matrix.
		for i in range(global_var.dim):
			for vec in self.matrix:
				print '{0:6.3f}'.format(vec.vec[i]),
			print '   ',
			print self.basis_list[i]
	
	def disp_eig_vals(self):
		# Prints the eigenvalues.
		for val in self.eig_vals:
			print '{0:6.3f}'.format(val)
	
	def disp_eig_vecs(self):
		# Prints the eigenvectors.
		for i in range(global_var.dim):
			for vec in self.eig_vecs:
				print '{0:6.3f}'.format(vec.vec[i]),
			print ''
	
	def acc_ground_state(self):
		return self.eig_vecs[0].vec
		
	def basis_occup(self):
		# Gets orbital occupation in each basis determinant.
		occup_indices = range(global_var.e_num-1, -1, -1)
		self.basis_list = []
		while True:
			occup = det.Det.indices_2_occup(occup_indices)
			self.basis_list.append(occup)
			occup_indices = det.Det.next_config_indices(occup_indices, True)
			if occup_indices == None:
				break
	
	def disp_ground_eig_vec(self):
		vec = self.eig_vecs[0].vec
		print ''
		print 'Exact eigenvectors'
		random_vec.RandomVec.disp_vec(vec)