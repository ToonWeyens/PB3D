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
4. Libraries to install:
    - HDF5: libhdf5-openmpi-dev
    - FFTW: libfftw3-dev

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

0.74: - Split off MPI_utilities from MPI_ops, containing the numerical utilities that have to do with MPI.
      - Put diff from utilities in output_ops, so now output_ops is indeed below utilities (but not MPI_utilities)
      - Con2dis and dis2con are now an error-reporting function.
      - Restructured the code to calculate HELENA on HELENA specified grid and afterwards interpolating on X grid.
      - Results are similar to previous ones when poloidal flux is used, but discrepancies for toroidal flux.
      - Simplified normalization: now input is normalized directly.
      - Fixed bug in DU_1: factor i/n or i/m was duplicated.
      - Fixed confusion about rho: Its profile is indeed free to be chosen and has no influence on the marginal stability.
      - Faulty solutions are removed by default; Can be overriden using retain_all_sol.

0.75: - Simplified HDF5 plotting: Group dimensions now read from variable size, individual version calls array version.
      - Corrected bug in HDF5 plotting: Symmetry is checked explicitely when one of the dimensions is 1.
      - Plotting of q, iota and other flux quantities now has correct normalization in both normal axis and function values.
      - Run scripts now create a new folder based on the date and time.
      - Added check on jacobian and magnetic field: Persistent error of a few percent, for various VMEC input parameters; Probably due to Fourier transform. Much better for HELENA.
      - Corrected a bug: the normal component of the curvature had the wrong sign in PV.
      - Improved normalization (and fixed bug for HELENA): Now mu_0 is also normalized so the equations don't change form.
      - As HELENA already provides normalized outputs, normalization is not necessary any more.
      - Improved running by putting output name in run script and not in PB3D itself.

0.76: - Implemented a test for the derivatives of g_E. There are some serious issues with the quality of the higher order numerical derivatives.
      - Fixed a bug in the calculation of the jacobian for Helena: The derivatives of the determinant of the transformation matrix are zero. Now the pressure balance for HELENA is good.
      - Reorganized the code by creating a new module for post processing, which is called after the main driver.
      - Moved the calculation of the extra quantities, such as shear, curvature, ... to the equilibrium module and improved the calculation of sigma.
      - Fixed a bug in the calculation of sigma.
      - HDF5 routines now use a type corresponding to REAL64 instead of deprecated H5T_NATIVE_DOUBLE.
      - Reorganized basic structure of programme.
      - Now there is a global variable prog_style that indicates whether some common routines are used for PB3D or PB3D_PP.
      - Now there is a global variable group_output that indicates whether all group mastes can output. If true, the output on screen indicates the outputting group and the file output is directed to the correct one.
      - MPI_wait can now also wait on all the groups.
      - Outputs revised: One text output per group, one HDF5 output and EV output per alpha job and level of Richardson's extrapolation.
      - Deallocation revised.
      - ERROR in richardson extrapolation: You can only do this for constant Eigenvalues! So need an algorithm to select these!
