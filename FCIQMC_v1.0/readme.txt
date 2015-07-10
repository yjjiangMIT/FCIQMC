This is a toy code for FCIQMC algorithm, version 1.0.
Description: FCIQMC without initiators. Reads in a randomly generated Hamiltonian from file and stores it as a dense matrix.

Author: Yijun Jiang
Starting date: 06-01-2015
Finishing date: 06-16-2015

Here is an example of running the code.
>>> import fciqmc
>>> parent_list = fciqmc.run_fciqmc(h, max_recursion, recursion_num_list, unsigned_num_list)