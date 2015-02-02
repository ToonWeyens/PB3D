cp ../PB3D . && mpirun -np $1 ./PB3D  cdxu cdxu -st_pc_factor_shift_type NONZERO -st_pc_type lu -st_pc_factor_mat_solver_package mumps $2
