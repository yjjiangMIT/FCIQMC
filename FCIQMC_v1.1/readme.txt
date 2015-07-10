This is a toy code for FCIQMC algorithm, version 1.1.
Description: FCIQMC with initiators. Reads in a randomly generated Hamiltonian from file and stores it as a sparse matrix in a Python dictionary.

Author: Yijun Jiang
Starting date: 06-17-2015
Finishing date: 06-26-2015

Here is an example of running the code.
>>> import run_simu
>>> h, cal_gnd_state, avg_shift = run_simu.run(orb_num, e_num, max_iter_num, crit_walker_num)