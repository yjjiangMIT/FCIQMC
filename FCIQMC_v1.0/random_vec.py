import random
import math
import global_var

class RandomVec:
	# Class of random vectors.
	lower = 0
	upper = 1
	
	def __init__(self, init_flag = True, rand_flag = True, vec = []):
		# Initiation of a random vector of length global_var.dim if init_flag == True,
		# else creates an all-zero vector of length global_var.dim.
		self.vec = []
		if init_flag:
			if rand_flag:
				for i in range(global_var.dim):
					self.vec.append(random.uniform(RandomVec.lower, RandomVec.upper))
			else:
				self.vec = vec
			self.normalize()
		else:
			for i in range(global_var.dim):
				self.vec.append(0)
	
	def normalize(self):
		# Normalization. Modifies self.
		norm = math.sqrt(self.dot_prod(self))
		self.scalar_prod(1/norm, True)
	
	def dot_prod(self, vec):
		# DOt product of self with a RandomVec object vec.
		result = 0
		for i in range(global_var.dim):
			result += self.vec[i] * vec.vec[i]
		return result
	
	def scalar_prod(self, multiplier, change_flag):
		# Scalar product of self with a number multiplier.
		# Modifies self if change_flag == True, else returns the result as a new RandomVec object.
		if change_flag:
			for i in range(global_var.dim):
				self.vec[i] *= multiplier
		else:
			result = RandomVec(False)
			for i in range(global_var.dim):
				result.vec[i] = (self.vec[i] * multiplier)
			return result
	
	def vec_add(self, vec, change_flag):
		# Entry-wise summation of self with a RandomVec object vec.
		# Modifies self if change_flag == True, else returns the result as a new RandomVec object.
		if change_flag:
			for i in range(global_var.dim):
				self.vec[i] += vec.vec[i]
		else:
			result = RandomVec(False)
			for i in range(global_var.dim):
				result.vec[i] = self.vec[i] + vec.vec[i]
			return result
	
	def vec_subtract(self, vec, change_flag):
		# Entry-wise subtraction of self with a RandomVec object vec.
		# Modifies self if change_flag == True, else returns the result as a new RandomVec object.
		if change_flag:
			for i in range(global_var.dim):
				self.vec[i] -= vec.vec[i]
		else:
			result = RandomVec(False)
			for i in range(global_var.dim):
				result.vec[i] = self.vec[i] - vec.vec[i]
			return result
	
	def schmidt(self, orth_vec_list):
		# Schmidt orthogonalization.
		# Makes self orthonormal to orth_vec_list, a list of RandomVec objects and appends self to it.
		for vec in orth_vec_list:
			self.vec_subtract(vec.scalar_prod(self.dot_prod(vec), False), True)
		self.normalize()
		orth_vec_list.append(self)
	
	@classmethod
	def disp_vec(self, vec):
		for i in range(len(vec)):
			print '{0:6.3f}'.format(vec[i]),
		print ''
