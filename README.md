PB3D
====

Peeling Ballooning in 3D

Doctoral work by Toon Weyens
Universidad Carlos III de Madrid
2012-2016

To set up:
1. git clone https://ToonWeyens@bitbucket.org/ToonWeyens/pb3d.git
2. change PB3D_DIR in the makefile
3. copy or add symbolic links to the output files (e.g. wout_cdxu or cbm18a)

CHANGELOG

0.71: - Added MPI support for yes_no
      - First steps to do away with module test. From now on test will be done within the program at runtime.
      - Improved greatly the run script. Now -d and -s can be used to run Valgrind debugging and optionally also error source tracking.
      - Extended calc_eqd_grid to handle 3D cases and implemented this throughout.
      - Fixed a bug in calc_ang_grid when using toroidal flux: zeta should be the 2nd index, not the first.
      - Moved extend_grid to grid_ops.
      - Implemented trim_grid, which trims the ghost regions of a grid.
      - Resorted to using trim_grid and ext_grid in many of the routines where plots are made.

0.72: - Implemented function to return relative or absolute difference between two inputs in utilities.
      - Implemented routine that plots 2 variables and their rel. and abs. differences
      - T_VF, jac_F is checked now for VMEC and it is found to be correct.
      - D1 p and D3 p are checked for both VMEC and HELENA and it is NOT correct.

0.73: - Fixed the confusion about the normal variables in F and E coords.
      - Improved plotting with GNUPlot, mostly in color handling.
      - Renamed yes_no to get_log and extended it with get_real and get_int.
      - Renamed broadcast_l to broadcast_log and extended it with broadcast_int, broadcast_real
      - Reimplemented generic tests of calc_deriv and conv_FHM to module test.
      - Fixed bug in calc_ang_grid: the parallel angle is always the first angle.
      - Tests now use global variables "debug_x" where "x" is the name of the routine tested.
      - Implemented successful testing of g_V.
      - Moved checking of HELENA to testing routine.
      - Implemented checking of Jac_F also for Helena and it is found to be correct.
