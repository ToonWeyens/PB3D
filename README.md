PB3D
====

Peeling Ballooning in 3D

Doctoral work by Toon Weyens
Universidad Carlos III de Madrid
ITER Organization
Technische Universiteit Eindhoven
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

0.77: - Created postprocessing programme PB3D_PP.
      - Removed the postprocessing module.
      - print_output_eq and print_output_X ready.
      - Implemented reading of variables from PB3D.
      - Prog_version is now a real variable.
      - The main PP driver is under process: It has three PB3D variables: one from output, one field-aligned and one plot-extended.
      - Removed usage of grp_r_X_loc in sol_ops, as X grid has no ghost regions and also because the routines should be called using a trimmed grid.
      - Implemented the post-processing program until decompose_energy.
      - However: Untested routines!

0.78: - Calc_real_XUQ now becomes calc_XUQ because the output is complex.
      - Moved calculation of rho to calc_flux_q, since it is a flux quantity.
      - Fixed a bug in plot_HDF5 when the variable name 'var' was used.
      - Fixed two bugs in calc_XUQ considering fac_0 and fac_1.
      - Implemented decompose_energy, but some output still missing.
      - Moved calc_int_magn to grid_ops and implemented calc_int_vol.

0.79: - When using Richardson extrapolation, the EV's for different levels are not displayed in log.
      - Implemented higher order expression (only) for d/dx in EV problem, to be used with norm_disc_style 2.
      - Changed the name from PB3D_PP to PB3D_POST.

0.80: - Fixed a bug in divide_X_grid for a high number of processes which would cause X_limits(2) to rise above n_r_X.
      - Added vacuum module vac, which will be used to calculate the vacuum response. Temporarily, this is set to zero.
      - Reordered the SLEPC routines to accomodate iteration for inverse.
      - Rewrote the SLEPC routines with the updated theory taking into account the plasma edge.
      - Fixed a bug in trim_grid when non-overlapping grids were used.
      - Now no more ghost points are needed for the calculations. However, for plotting and writing output a ghost region of one is kept.
      - HDF5 is now the only possibility, but some outputs are given with GNUPlot.
      - draw_GP and draw_GP_animated now have the possibility to draw decoupled 3D plots using draw_dim which replaces is_2D.
      - Introduced parameter GP_max_size, above which GNUPlot is not used.
      - Changed plot_jq to plot_resonance. It also outputs the resonating flux surfaces in 2D HDF5.
      - Harmonic plot now also plots the resonatic flux surfaces in 2D.

0.81: - Fixed inconsistencies in plot_X_vec and decompose_energy: the numerical derivatives need two-sided ghost region.
      - Introduced symmetric ghost regions in PB3D_POST with width ghost_width_POST, as this is needed for the numerical derivatives of X.
      - Accompanying this, trim_grid is extended with shift_grid.
      - Fixed bug in calc_XYZ_real where flux_p_H was not divided by 2 pi.

0.82: - Fixed an error in fill_mat for order 2: factor 12 in stead of 2.
      - Fixed errors in calculating d_nz and o_nz.
      - Changed norm_disc_style to norm_disc_ord, only referring to the order.
      - Made fill_mat and set_BC generic for any order. The latter now employs set_left_BC and set_right_BC.

0.83: - BC_style chooses which style for left and right boundary.
      - Eigenvectors are now normalized with respect to their maximum modulus.

0.84: - The output name of draw_GP(_animated) is now consistent with the system used for HDF5: a new variable draw_name is introduced.
      - Energy decomposition output is now written in a file PB3D_out_EV.txt
      - HDF5 output now shows less messages.
      - PB3D_POST run script now takes PB3D_out.h5 from a directory specified and by default runs in that directory.

0.85: - Reduced (eliminated?)  memory leaks by properly nullifying and deallocating pointers.
      - Use_normalization cannot any more be passed to PB3D_POST, but is instead read from PB3D output.
      - Corrected bug where mu_0 was forgotten in kappa_n.

0.86: - Added new option 'debug_X_grid', which sets the perturbation grid equal to the equilibrium grid.
      - Added new test 'test_diff' which looks for the effects of introducing artificial numerical diffusion.
      - Implemented the 'debug_store_results' which checks whether the solution stored corresponds to the relative error returned.

0.87: - Changed the routine 'divide_X_grid', which now returns the whole r_F instead of just grp_r_F.
      - Changed the way the routine 'solve_EV_system_SLEPC' calculates the step size: It now makes use of grid_X%r_F.
      - Corrected a VERY IMPORTANT bug by setting V_2 to - V_2.
      - extremely uncorrect unstable spectrum largely corrected.

0.88: - Introduced variables 'tol_slepc' and 'max_n_it_slepc', which are displayed before starting the solver.
      - Fixed bugs considering BC's.
      - Energies are not scaled by E_kin any more.

0.89: - Improved makefile, adding support for quadrivium compilation.
      - Changed things in 'output_ops' and 'sol_ops' that provoked a compilation bug on quadrivium: an array could not be used as a building block of a bigger array as in y = [x,5] where x = [1,2].
      - Split the perturbation output in perturbation output and solution output, which are written to 2 different HDF5 groups.
      - 'exp_ang_par_F' is now 'J_exp_ang_par_F' and contains the Jacobian as well.
      - For PB3D, the metric variables are no longer interpolated on the field-aligned grid, as they are no longer necessary with the introduction of 'J_exp_ang_par_F'.
      - Introduced a new run script, for qsub on quadrivium.

0.90: - LAST VERSION OF OLD PB3D WITH PARALLELIZATION IN THE NORMAL COORDINATE.
      - Small bug fixes.
      - Started implementing memory check.

0.91: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS.
      - Replaced the routines that calculate the possible derivatives of a certain order by a generic one.
      - Implemented calc_derivs_1D_id which calculates the 1D index of derivatives. A table of these is stored in an array d for quick access.
      - 'driver_rich.f90' has been renamed 'driver.f90' and the original 'driver.f90' was deleted.
      - Introduced 'norm_disc_style' for eq, X and sol, so derivatives are done coherently.
      - 'trim_grid' now optionally returns 'norm_id', the normal indices of the trimmed variables.
      - The derivatives of VMEC outputs is now done in 'prepare_RZL' because the precision is not yet known before reading the input.
      - 'calc_eq_r_range', 'divide_X_grid' and 'split_MPI' are now bundled in 'calc_norm_range'.
      - 'reconstuct_PB3D' is also called in the pertubation phase, but only equilibrium variables are reconstructed.
      - X variables now are all allocatable.
      - The X jobs can either be found in a straightforward way, requiring the use of the mpirun option '--mca osc pt2pt', or in an intelligent way, but this requires an up-to-date version of openmpi. The flag 'lold_MPI' can be used to switch between the two versions.

0.92: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
      - Removed the whole system using MPI to transmit which process does which X job. Now, it works using an external file as well as an external lock file.
      - Fixed a bug when the wrong grid variables were written out.
      - Split of the transformation of derivatives in F.
      - Moved the plotting of flux quantities as these require knowledge of the derivative in F.
      - Removed plot rank, so now rank is used everywhere.
      - Renamed 'check_modes' to 'check_X_modes'.
      - 'check_X_modes' and 'resonance_plot' do not use an X_type anymore.
      - 'divide_X_jobs' now contains 'calc_memory'. Also it can be used in vector or tensor mode and in theory higher orders as well, if information about the number of variables at these orders is provided.
      - The perturbation phase is to be split in two: a vector phase and a tensor phase. This is not yet implemented.

0.93: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
      - 'prepare_X' has been replace by 'calc_X' which calculates the perturbation variables for a certain order.
      - There are two orders of perturbation variables. These are calculated separately and the results of order 1 are read from disk when calculating those of order 2.
      - 'print_HDF5_arrs' now detects whether it is an individual or collective call. Also, it checks whether a group already exists.
      - PB3D files are not opened anymore like normal input nor are their 1D variables broadcasted: Every process opens the file on its own.
      - Restored 'prog_style' to only 2 possibilities: 1 (PB3D) or 2 (POST)
      - Split 'X_type' into 'X_1_type' (vectorial perturbation), 'X_2_type' (tensorial perturbation) and 'sol_type' (solution).
      - Code working up to calcalation of tensorial perturbation quantities.

0.94: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
      - HDF5 variables are now bundled in 'HDF5_vars'.
      - 'print_output_sol' is now situated in 'sol_ops'.
      - 'fourier_ops' is now called 'fourier'.
      - 'read_PB3D' takes optional mode information which is passed to 'read_HDF5_arrs' which takes optional string arguments.
      - there are now 5 types of variables--misc, eq, X_1, X_2 and sol--that can be printed and reconstructed independently.
      - 'open_output' now also writes miscelleaneous variables to HDF5 file.
      - Fixed bug in 'read_HDF5_arrs' where head group was not closed which causd problems for consequent reads.
      - 'c' now does its job also for submatrices, when the limits are provided.
      - Tensorial perturbation variables are now correctly calculated up to KV, PV.

0.95: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
      - Matching modes now lead to less reading of variables.
      - File and variable names have changed to consistent 'eq', 'X', 'sol' and 'POST'.
      - Solution driver template introduced.
      - In the perturbation driver, the only thing missing is the integration along the magnetic field lines.
      - In the solution driver evetything needs to be pieced together.

0.96: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
      - 'interp_HEL_on_grid' now overwrites input quantities if no outputs are provided.
      - This way the doubling of the memory requirement for HELENA is avoided.
      - The calculation of magnetic integrals now happens completely in 'X_ops'.
      - 'PB3D_type' is not used anymore; Instead, individual variables are used.
      - Perturbation driver has been completed.

0.97: - NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
      - Updated 'X_ops' to not use 'PB3D_type'.
      - Implemented the X_2 versions of the routines that read and reconstruct PB3D variables.
      - Changed 'perturbation grid' to 'solution grid'.
      - 'setup_and_calc_grid_B' has been renamed to 'setup_and_calc_grid_eq_B' and 'setup_and_calc_grid_sol' has been implemented.
      - Completed the solution driver routine.
      - Updated SLEPC routines to handle the solution type.
      - Fixed a bug in 'trim_grid' where the limits of the trimmed grid were not updated.
      - Introduced 'wait_file' which makes use of a lock file and used this in the writing of HDF5 files.

0.98: - FIRST NEW USABLE VERSION: NEW PARALLELIZATION STRATEGY.
      - The miscellaneous variables are now read in 'open_output' for POST, just as they are written there for PB3D.
      - The PB3D operations have to be initialized before using them. This sets the equilibrium style.
      - Modified 'broadcast_input_vars' to also broadcast miscellaneous variables in POST.
      - The maximum fluxes are now written in the equilibrium variables.

0.99: - Implemented a plot of the comparison between the Eigenvalue and the energy fraction in 'plot_X_val_comp'.
      - Bug fixes in calculation of ranges and 'no_plots' for POST.
      - Bug fix in plotting of resonant surfaces.
      - 'norm_disc_prec_sol' is now used in PB3D to indicate the discretization precision of the matrices A and B.
      - SLEPC routines now need a trimmed solution grid, and an untrimmed solution grid is nowhere necessary in PB3D.

1.00: - FIRST WORKING VERSION
      - Bug fixed where Hermitian conjugate was stored in stead of the real matrices in SLEPC.
      - Checked consistency of results when "intuitive" version of plasma potential energy is used.

1.01: - UNUSABLE VERSION: See last comment.
      - Plots can now have an optional different size.
      - Corrected bug with wrong sign in calculation of jacobian in HELENA coordinates.
      - Merged contents fourier with VMEC modules.
      - Fixed bugs in the HELENA interpolation.
      - Changed 'read_HEL' for asymmetric cases to duplicate information for theta equal to 0 into theta equal to 2pi.
      - End-points are now included for the parallel angle, to simplify working with trapezoidal integration. This eliminates the discrepancies between results for different parallel intervals. N.B.: NOT VALID for some other integration rules!
      - Added a test folder to be used for FORTRAN tests.
      - Started implementing the transformation of the metric factors. This is NOT finished and will probably be abandoned

1.02: - Metric variables are interpolated still for HELENA, but in a different, more clever way. Not yet implemented for toroidal flux.
      - Bounds are now checked.
      - Fixed and finished implementing the comparison between U and U_inf.
      - Extended 'debug_store_results' so that it displays X*AX and X*BX separately.
      - Extended debugging of calculation of integral. An alternative calculation for orthogonal computational grids shows identical results.

1.03: - Fixed a bug in the calculation of the volume integrals: ghost region was needed.
      - Consistency of volume integral with both repeated trapezoidal integration and repeated simple integration shown.
      - Introduced normalization style 'norm_style' that now also allows for normalization of only the normal component of the perturbation.
      - Various bugfixes when no resonant surfaces are found or when solutions are removed.
      - Solution vectors are now normalized in store loop, to ensure consistency of X*VX step_size and the integrals of POST.
      - Fixed a bug in 'calc_memory' where overflow occured.

1.04: - UNUSABLE VERSION: CHECKING INDIVIDUAL TERMS WITH ENERGY DECOMPOSITION.
      - A more detailed grid for X quantities is needed!
      - Cleaned up the code when 'ldebug' is not used.
      - Implemented debugging of the perturbation driver so that the vectorial instead of the tensorial perturbation variables are interpolated.
      - Implemented the option to choose order of geodesic perturbation.
      - Implemented a test to plot |X|^2 after solving the Eigenvalue problem and compare this with the energy reconstruction.
      - Fixed a bug where the transpose of a matrix should be used when feeding into SLEPC.

1.05: - ONLY PB3D USUABLE.
      - POST has to be adapted still.
      - Reintroduced all terms in the energy equations.
      - New module for Richardson extrapolation.
      - Extrapolation is now done in main program, as X variables have to recalculated on a different grid for each Richardson level.
      - Involved routines "calc_U", "calc_KV" and "calc_PV" have been rewritten. "calc_U" is still not tested for VMEC.
      - Updated the user output considering Richardson extrapolation in general.
      - Failure to read input file now results in error.
      - 'test_max_mem' now indicates whether memory is to be tested or not.

1.06: - FIRST USUABLE VERSION WITH CORRECT ENERGY DECOMPOSITION: PERTURBATION VARIABLES NOW CORRECTLY TABULATED IN OWN GRID.
      - Adapted POST to new X storage convention as well.
      - Fixed some bugs considering Richardson extrapolation.
      - POST now by default loads the highest found Richardson level, but this can be overriden using "PB3D_rich_lvl".
      - Bugfixes considering the new tabulation of perturbation variables.
      - "debug_store_results" is now part of the standard programme.
      - HELENA always uses normalization, so "use_normalization" is always set to true.

1.07: - FIRST VERSION WITH FAST OPTION. It seems to work for cbm18a if "max_r_sol" is not larger than 0.9.
      - Introduced variable "X_style", which allows the user to choose between prescribed modes (earlier default) and a fast version.
      - "check_X_modes" now has alternative actions for X style 2.
      - The original versions of "get_suffix", "set_nn_mod" and "is_necessary_X" do not exist any more: Everything works with X limits now.
      - Split some routines of PB3D_vars into new module PB3D_utilities.
      - Split some routines of X_vars into new module X_utilities.
      - Split some routines of SLEPC_ops into new module SLEPC_utilities.
      - Split some routines of input_ops into new module input_utilities.
      - Split some routines of grid_ops into new module grid_utilities.
      - Split some routines of sol_ops into new module sol_utilities, which can also translate between local and total solution vectors.
      - Derivatives of the solution vectors have to be done on total variables.
      - Suffixes to names are now integer variables, not characters.
      - The limits of the modes and their values now have to be calculated from the equilibrium and have to be recalculated using "setup_nm_X".
      - Fixed a bug in "interp_fun".
      - Improved debugging of "calc_zero_NR".
      - Fixed a bug in "calc_res_surf" for dfun.
      - Introduced local versus total mode numbers, which may differ for X style 2.
      - "plot_X_vec" has been modified as well as "resonance_plot" to correctly implement the difference between local and total mode numbers.
      - Fixed inconsistency with naming "X_vec" and "X_val", etc. "sol" Is now used.

1.08: - NOT USABLE FOR VMEC: ERRORS PRESENT!
      - Implemented routine to get info about file, though unused for now.
      - New flag "no_execute_command_line" that disables the execution of command line, as this sometimes causes a freeze, for unknown reasons.
      - Identified bug when reading VMEC with "lrfp" flag. Now, netcdf output is used instead, as it does not have this problem.
      - Fixed some bugs in "reconstruct_PB3D" considering the presence of field-aligned grids. These are now automatically set for both equilibrium styles.
      - Fixed a bug in the first derivative of "flux_p_E" for HELENA: factor 2pi is now taken into account.
      - Changed the storage convention for "PV_int" and "KV_int".
      - Generalized and automatized the "calc_deriv" routines and put them in "grid_utilities"; They can take any order and precision now.
      - Cleaned up "calc_flux_q".
      - Split HELENA in HELENA_ops and HELENA_vars.
      - Added interpolation capability similar to the derivation capability.
      - Both interpolation and derivation use "apply_disc" to apply the discretization operator.
      - "interp_fun" is still available to be used only when the interpolations are not repetitive and linear interpolation suffices.

1.09: - FIRST VERSION THAT WORKS FOR AXISYMMETRIC VMEC. However, post-processing is not working well, and comparision with HELENA not yet perfect.
      - Normalization factors "R_0", "B_0", "pres_0", "psi_0" can now be user-provided, just like "rho_0", but NO checking for consistency!
      - Fixed a bug in "calc_interp_data" where the bounds of "x" were not respected.
      - Fixed bug in Post-processing for VMEC input where the plot output is calculated.
      - Removed "test_Dg_E" as it does not make sense: In general the equilibrium grid is not straight.
      - Removed the "repack" routine and started using the natural VMEC structure.
      - Fixed an error in the testing of B where only components in equilibrium coordinates were not properly translated to flux coordinates.
      - The structure of "SLEPC_ops" has been slightly rewritten for the future accomodation of shell matrices.
      - Fixed bug in "calc_U" where the part ~ n_frac was not used in the calculation of DU_1 and DU_2, as well as some erroneous indices.
      - Moved some of the routines from "calc_eq" to "calc_met". In fact, these routine structures should be rethought, possibly merged.
