import key_ops
import black_box

class Det:
	# Class of a walker site.
	# $value$ is the unsigned number of walkers on this site.
	# $flag$ is used to indicate initiators.
	
	def __init__(self, value, flag):
		"""Initiates a $Walker$ object."""
		
		self.value = value
		self.flag = flag
		self.diag_entry = 0
	
	def if_survive_as_child(self):
		return self.flag or abs(self.value) > 1
	
	def set_diag_entry(self, key):
		"""Sets diagonal matrix element of this determinant."""
		orbs = key_ops.key_2_orbs(key)
		self.diag_entry = black_box.sandwich(orbs, ())