import integrals
import key_ops
import test

# All the knobs to tune.
change_shift_crit_num = 20000
init_crit_w_num = 5
init_walker_num = 50
max_iter_num = 10000
wait_for_aver_num = 2000

damping = 0.05
init_shift = 0.05
# ref_key = (1, 15, 1, 15)
ref_key = (4, 5, 1)
single_prob = 0.5
tau = 0.005

change_shift_step = 5
update_plots_step = 20

# The exact correlation energy is found by directly diagonalizing the Hamiltonian on MATLAB.
# exact_corr = -0.1543290876
exact_corr = -0.3118599099

# The ranges of various plots.
orb_num = integrals.para_list[2]
e_num = integrals.para_list[3]
y_axis_eig_vec_plot = [-1, 1]
y_axis_energy_plot = [-0.5, 0.1]
y_axis_log_w_num_plot = [3, 11]

# Sets the reference energy.
integrals.ref_energy = integrals.sandwich(key_ops.key_2_orbs(ref_key), ())

key_list = test.gen_key_list(ref_key)

# Begins the simulation.
import run_simu
aver_shift, aver_proj, dets_p, vec = run_simu.run()