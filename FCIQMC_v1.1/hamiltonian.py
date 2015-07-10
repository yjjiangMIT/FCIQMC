class Hamiltonian:
	# Random Hamiltonian class: symmetric, positive semi-definite and whose smallest eigenvalue vanishes.
	# Its eigenvectors and eigenvalues are stored as $eig_vecs$ and $eig_vals$.
	# The Hamiltonian itself is stored as a matrix.
	# A sparse matrix is a full list of sparse vectors.
	
	def __init__(self, file_name):
		"""Initializes a $Hamiltonian$ object."""
		
		import key_ops
		import read_h
		import vec_ops
		
		self.para_list = []
		self.eig_vecs = []
		self.eig_vals = []
		self.key_list = []
		
		# Loads data from a file.
		read_h.read_outof_file(file_name, self.para_list, self.eig_vecs, self.eig_vals)
		# Makes an ordered list of keys.
		key_ops.create_keys(self.key_list, self.para_list)
		
		# Computes H=sum(Ei*|vi><vi|)
		self.matrix = [{}] * self.para_list[4]
		for i in range(self.para_list[4]):
			vec = self.eig_vecs[i]
			for key in vec:
				index = key_ops.key_2_index(self.key_list, key)
				column = vec_ops.vec_scalar_prod(vec, vec[key]*self.eig_vals[i] , False)
				self.matrix[index] = vec_ops.vec_add(self.matrix[index], column, False)
		
		# Reference determinant and reference energy.
		# Selects the lowest-energy determinant.
		self.ref_index = 0
		self.ref_key = self.key_list[0]
		self.ref_energy = self.matrix[0][self.ref_key]
		for i in range(1, self.para_list[4]):
			try:
				temp_ref_energy = self.matrix[i][self.key_list[i]]
			except KeyError:
				temp_ref_energy = 0
			if temp_ref_energy < self.ref_energy:
				self.ref_index = i
				self.ref_key = self.key_list[i]
				self.ref_energy = temp_ref_energy
	
	def exact_gnd_state(self):
		"""Gets the exact ground state."""
		return self.eig_vecs[0]
	     
	def disp_mat(self):
		"""Displays the Hamiltonian as a dense matrix."""
		
		import key_ops
		
		for vec in self.matrix:
			array = [0] * self.para_list[4]
			for key in vec:
				index = key_ops.key_2_index(self.key_list, key)
				array[index] = vec[key]
			for i in array:
				print '{0:6.3f}'.format(i),
			print ''