import black_box
import key_ops

change_shift_crit_num = 20000
init_crit_w_num = 5
init_walker_num = 50
max_iter_num = 50000

damping = 0.05
init_shift = 2
ref_key = (1, 15, 1, 15)
single_prob = 0.5
tau = 0.0025

change_shift_step = 5
update_plots_step = 20

exact_gnd = -0.1543290876

axis_eig_vec_plot = [0, (2**black_box.para_list[0])**black_box.para_list[1] - 1, -1, 1]
axis_energy_plot = [0, max_iter_num, -0.3, 0.1]
axis_walker_num_plot = [0, max_iter_num, 4, 11]

black_box.ref_energy = black_box.sandwich(key_ops.key_2_orbs(ref_key), ())

import run_simu
aver_shift, dets_p, vec = run_simu.run()