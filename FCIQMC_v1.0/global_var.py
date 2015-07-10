import math
# All the knobs to turn.

test_print_flag = False
rand_hamiltonian_flag = False
change_shift_flag = False
bit_num = 4 # The number of bits in a storage unit. e.g. 64.
spin_orb_num = 10 # The number of spin-orbitals. i.e. 2M.
e_num = 3 # The number of electrons, spin not taken into account at this stage. i.e. N = N_alpha + N_beta.
dim = math.factorial(spin_orb_num) / math.factorial(e_num) / math.factorial(spin_orb_num - e_num)

tau = 1e-2
shift = 0
nc = 1000
change_shift_step = 5
damping = 0.1

file_name = 'Documents/FCIQMC/h_file_5_3.txt'