import global_var
global_var.rand_hamiltonian_flag = True
import fciqmc
import stored_h

h_file_write = open(global_var.file_name, 'w')
stored_h.write_h_file(fciqmc.h, h_file_write)
h_file_write.close()