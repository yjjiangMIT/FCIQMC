class Walker:
	# Class of a walker site.
	# $value$ is the unsigned number of walkers on this site.
	# $flag$ is used to indicate initiators.
	
	def __init__(self, value, flag = False):
		"""Initiates a $Walker$ object."""
		
		self.value = value
		self.flag = flag