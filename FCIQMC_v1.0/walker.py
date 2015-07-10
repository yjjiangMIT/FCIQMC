import load_data
import random
import global_var
import math

class Walker:
        # Class of walker sites. Each site can have multiple walkers.
        # The controlling parameters here are tau and shift.
	unsigned_num = 1
	
	def __init__(self, h, loc, sign = True):
                # Initiates a Walker object.
		self.loc = loc
		if sign:
			self.signed_num = 1
		else:
			self.signed_num = -1
		self.diag = load_data.get_entry(h, self.loc.occup, self.loc.occup)
		
		# For test only.
		if global_var.test_print_flag:
			print 'Created walker at  ',
			print self.loc.occup
		
	def spawn(self, h):
                # Spawning of one parent on this walker site.
                # Returns a Walker object named child.
		select_prob, spawn_to = self.loc.excite(h.basis_list, 'r')
		child = None
		child_num = 0
		
		# For test only.
		if global_var.test_print_flag:
			print 'Start of spawning step'
			print 'Attempt to spawn at',
			print spawn_to.occup

		# Gets the probability of spawning.
		# It may be larger than unity, which may imply multiple spawning.
		prob = global_var.tau * load_data.get_entry(h, self.loc.occup, spawn_to.occup)\
			/ select_prob
		if prob < 0:
                        # Children have the same sign as their parent.
			child_sign = (self.signed_num > 0)
			prob = -prob
		else:
                        # Children have the opposite sign as their parent.
			child_sign = (self.signed_num < 0)
		child_num = int(prob)
		if child_num > 0:
			# Spawns child_num children successfully.
			child = Walker(h, spawn_to)
			child.signed_num = child_num
			if not child_sign:
				child.signed_num = -child.signed_num
		rand_val = random.random()
		
		# For test only.
		if global_var.test_print_flag:
			for i in range(child_num):
				print 'Probability is  ',
				print prob - i
				print 'Spawn successfully'
			print 'Probability is  ',
			print prob - child_num
			print 'Random number is',
			print rand_val
		
		if rand_val < prob - child_num:
			# Another spawning is successful.
			
			# For test only.
			if global_var.test_print_flag:
				print 'Spawn successfully'
			
			if child_num == 0:
				child = Walker(h, spawn_to, child_sign)
			elif child_sign:
				child.signed_num += 1
			else:
				child.signed_num -= 1
		else:
			# Spawning is unsuccessful.
			
			# For test only.
			if global_var.test_print_flag:
				print 'Spawn unsuccessfully'
		
		# For test only.
		if global_var.test_print_flag:
			if child != None:
				print 'Number of living walkers on this site:',
				print child.signed_num
			print 'End of spawning step'
			print ''
		
		return child
	 
	def die(self, h):
                # The death/cloning process of one parent on this walker site.
		prob = global_var.tau * (load_data.get_entry(h, self.loc.occup, self.loc.occup)\
			- h.ref_energy - global_var.shift)
		rand_val = random.random()
		
		# For test only.
		if global_var.test_print_flag:
			print 'Start of dying step'
			print 'Probability is  ',
			print prob
			print 'Random number is',
			print rand_val
		
		if prob < 0:
			# Cloning.
			
			# For test only.
			if global_var.test_print_flag:
				print 'In fact it is cloning'
			
			prob = -prob
			clone_num = int(prob)
			if rand_val < prob - clone_num:
				# Another cloning is successful.
				
				# For test only.
				if global_var.test_print_flag:
					print 'One walker is cloned.'
				
				clone_num += 1
			if self.signed_num > 0:
				self.signed_num += clone_num
			else:
				self.signed_num -= clone_num
		else:
			# Death.
			if rand_val < prob:
				# Dead.
				
				# For test only.
				if global_var.test_print_flag:
					print 'One walker dies.'
				
				if self.signed_num > 0:
					self.signed_num -= 1
				else:
					self.signed_num += 1
			else:
				
				# For test only.
				if global_var.test_print_flag:
					print 'No one dies.'
		
		# For test only.
		if global_var.test_print_flag:
			print 'Number of living walkers on this site:',
			print self.signed_num
			print 'End of dying step'
			print ''
	
	def is_dead(self):
                # Decides if a walker site contains no survivers.
		if self.signed_num == 0:
			# No walker survives.
			return True
		else:
			# Still some walkers survive.
			return False
	      
	def annihilate(self, w):
		self.signed_num += w.signed_num
		
		# For test only.
		if global_var.test_print_flag:
			print 'Start of annihilation step'
			print 'Site',
			print self.loc.occup
			print 'Number of living walkers on this site:',
			print self.signed_num
			print 'End of annihilation step'
			print ''
	
	def compare(self, w):
		for i in range(len(self.loc.occup)):
			diff = self.loc.occup[i] - w.loc.occup[i]
			if diff < 0:
				# self is smaller.
				return -1
			elif diff > 0:
				# self is larger.
				return 1
		# They are the same for each digit.
		return 0
	
	@classmethod
	def spawn_and_die_for_all(self, h, parent_list):
		child_list = []
		i = 0
		while i < len(parent_list):
			# Loops over all the parent sites.
			parent_site = parent_list[i]
			number = parent_site.signed_num
			if number < 0:
				number = -number
			for j in range(number):
				# Loops over all the parents on one site.
				# Parents first spawn, then die.
				child = parent_site.spawn(h)
				if child != None:
					child_list.append(child)
				parent_site.die(h)
			if parent_site.is_dead():
				parent_list.pop(i)
			else:
				i += 1
		return child_list
	
	@classmethod
	def bisec_search(self, w_list, w):
                # Bisection search.
                # Given a walker w, search for its appearence in a sorted list.
                # If w exists in the list, returns True and its position.
                # Else returns False and its insert position.
                # This method is called recursively.
                length = len(w_list)
                if length == 0:
                        # Base case 1.
                        return False, 0
                if length == 1:
                        # Base case 2.
                        comparison = w.compare(w_list[0])
                        if comparison == 0:
                                # Finds a match.
                                return True, 0
                        elif comparison > 0:
                                # No match, and w is behind current position.
                                return False, 1
                        else:
                                # No match, and w is at current position.
                                return False, 0
                else:
                        # Recursion.
                        middle = length / 2
                        comparison = w.compare(w_list[middle])
                        if comparison == 0:
                                # Finds a match.
                                return True, middle
                        elif comparison > 0:
                                # Search the second half of the list.
                                is_found, pos_in_sub = Walker.bisec_search(w_list[middle+1:], w)
                                return is_found, pos_in_sub + middle + 1
                        else:
                                # Search the first half of the list.
                                return Walker.bisec_search(w_list[:middle], w)

	@classmethod
	def merge_walkers(self, parent_list, child_list):
                # Merge the parent list with the newly born child list.
                # Deal with annihilations in the process.
                for child in child_list:
                        is_found, position = Walker.bisec_search(parent_list, child)
                        if is_found:
                                # Annihilation.
                                parent = parent_list[position]
                                parent.annihilate(child)
                                if parent.is_dead():
                                        # Annihilated.
                                        parent_list.pop(position)
                        else:
                                # Create a new site.
                                parent_list.insert(position, child)

	@classmethod
        def disp_w_list(self, w_list):
                for w in w_list:
                        print w.loc.occup, w.signed_num
                        
        @classmethod
        def count_unsigned_num(self, w_list):
		total = 0
		for w in w_list:
			total += int(math.fabs(w.signed_num))
		Walker.unsigned_num = total
		
	@classmethod
	def cal_ground_eig_vec(self, w_list):
		norm_sq = 0
		vec = []
		for w in w_list:
			norm_sq += w.signed_num ** 2
		norm = math.sqrt(norm_sq)
		for w in w_list:
			vec.append(w.signed_num / norm)
		return vec
