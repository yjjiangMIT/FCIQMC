import key_ops
import integrals

class Det:
	# Class of a determinant.
	# value is the unsigned number of walkers on this determinant.
	# flag is used to indicate (1) if a parent is an initiator (2) if a child is spawned by an initiator.
	
	def __init__(self, value, flag):
		"""Initiates a Det object."""
		
		self.value = value
		self.flag = flag
		self.diag_entry = 0
	
	def if_survive_as_child(self):
		"""Judges if the spawned walkers on this determinant can survive."""
		
		# Survival criterion: (1) child spawned by an initiator (self.flag == True)
		#                  or (2) more than one children
		return self.flag or abs(self.value) > 1
	
	def set_diag_entry(self, key):
		"""Sets diagonal matrix element of this determinant."""
		
		# Excludes the reference energy of D0.
		orbs = key_ops.key_2_orbs(key)
		self.diag_entry = integrals.sandwich(orbs, ()) - integrals.ref_energy