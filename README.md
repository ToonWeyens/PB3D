# PB3D

*Peeling Ballooning in 3D*

Doctoral work by Toon Weyens
Universidad Carlos III de Madrid
Technische Universiteit Eindhoven
ITER Organization
2012-2017

## To set up:
1. git clone https://ToonWeyens@bitbucket.org/ToonWeyens/pb3d.git
2. change PB3D_DIR in the makefile
3. copy or add symbolic links to the output files (e.g. wout_cdxu or cbm18a)
4. Install
    * petsc
    * slepc
    * hdf5
    * netcdf
    * libstell
    * StrumPack
5. make all

## To run
1. Go to Run/
1. run 'setup_Run.sh' and point to the folder.

## Changelog

## 0.71:
* Added MPI support for yes_no
* First steps to do away with module test. From now on test will be done within the program at runtime.
* Improved greatly the run script. Now -d and -s can be used to run Valgrind debugging and optionally also error source tracking.
* Extended calc_eqd_grid to handle 3D cases and implemented this throughout.
* Fixed a bug in calc_ang_grid when using toroidal flux: zeta should be the 2nd index, not the first.
* Moved extend_grid to grid_ops.
* Implemented trim_grid, which trims the ghost regions of a grid.
* Resorted to using trim_grid and ext_grid in many of the routines where plots are made.

## 0.72:
* Implemented function to return relative or absolute difference between two inputs in utilities.
* Implemented routine that plots 2 variables and their rel. and abs. differences
* T_VF, jac_F is checked now for VMEC and it is found to be correct.
* D1 p and D3 p are checked for both VMEC and HELENA and it is NOT correct.

## 0.73:
* Fixed the confusion about the normal variables in F and E coords.
* Improved plotting with GNUPlot, mostly in color handling.
* Renamed yes_no to get_log and extended it with get_real and get_int.
* Renamed broadcast_l to broadcast_log and extended it with broadcast_int, broadcast_real
* Reimplemented generic tests of calc_deriv and conv_FHM to module test.
* Fixed bug in calc_ang_grid: the parallel angle is always the first angle.
* Tests now use global variables 'debug_x' where 'x' is the name of the routine tested.
* Implemented successful testing of g_V.
* Moved checking of HELENA to testing routine.
* Implemented checking of Jac_F also for Helena and it is found to be correct.

## 0.74:
* Split off MPI_utilities from MPI_ops, containing the numerical utilities that have to do with MPI.
* Put diff from utilities in output_ops, so now output_ops is indeed below utilities (but not MPI_utilities)
* Con2dis and dis2con are now an error-reporting function.
* Restructured the code to calculate HELENA on HELENA specified grid and afterwards interpolating on X grid.
* Results are similar to previous ones when poloidal flux is used, but discrepancies for toroidal flux.
* Simplified normalization: now input is normalized directly.
* Fixed bug in DU_1: factor i/n or i/m was duplicated.
* Fixed confusion about rho: Its profile is indeed free to be chosen and has no influence on the marginal stability.
* Faulty solutions are removed by default; Can be overriden using retain_all_sol.

## 0.75:
* Simplified HDF5 plotting: Group dimensions now read from variable size, individual version calls array version.
* Corrected bug in HDF5 plotting: Symmetry is checked explicitely when one of the dimensions is 1.
* Plotting of q, iota and other flux quantities now has correct normalization in both normal axis and function values.
* Run scripts now create a new folder based on the date and time.
* Added check on jacobian and magnetic field: Persistent error of a few percent, for various VMEC input parameters; Probably due to Fourier transform. Much better for HELENA.
* Corrected a bug: the normal component of the curvature had the wrong sign in PV.
* Improved normalization (and fixed bug for HELENA): Now mu_0 is also normalized so the equations don't change form.
* As HELENA already provides normalized outputs, normalization is not necessary any more.
* Improved running by putting output name in run script and not in PB3D itself.

## 0.76:
* Implemented a test for the derivatives of g_E. There are some serious issues with the quality of the higher order numerical derivatives.
* Fixed a bug in the calculation of the jacobian for Helena: The derivatives of the determinant of the transformation matrix are zero. Now the pressure balance for HELENA is good.
* Reorganized the code by creating a new module for post processing, which is called after the main driver.
* Moved the calculation of the extra quantities, such as shear, curvature, ... to the equilibrium module and improved the calculation of sigma.
* Fixed a bug in the calculation of sigma.
* HDF5 routines now use a type corresponding to REAL64 instead of deprecated H5T_NATIVE_DOUBLE.
* Reorganized basic structure of programme.
* Now there is a global variable prog_style that indicates whether some common routines are used for PB3D or PB3D_PP.
* Now there is a global variable group_output that indicates whether all group mastes can output. If true, the output on screen indicates the outputting group and the file output is directed to the correct one.
* MPI_wait can now also wait on all the groups.
* Outputs revised: One text output per group, one HDF5 output and EV output per alpha job and level of Richardson's extrapolation.
* Deallocation revised.
* ERROR in richardson extrapolation: You can only do this for constant Eigenvalues! So need an algorithm to select these!

## 0.77:
* Created postprocessing programme PB3D_PP.
* Removed the postprocessing module.
* print_output_eq and print_output_X ready.
* Implemented reading of variables from PB3D.
* Prog_version is now a real variable.
* The main PP driver is under process: It has three PB3D variables: one from output, one field-aligned and one plot-extended.
* Removed usage of grp_r_X_loc in sol_ops, as X grid has no ghost regions and also because the routines should be called using a trimmed grid.
* Implemented the post-processing program until decompose_energy.
* However: Untested routines!

## 0.78:
* Calc_real_XUQ now becomes calc_XUQ because the output is complex.
* Moved calculation of rho to calc_flux_q, since it is a flux quantity.
* Fixed a bug in plot_HDF5 when the variable name 'var' was used.
* Fixed two bugs in calc_XUQ considering fac_0 and fac_1.
* Implemented decompose_energy, but some output still missing.
* Moved calc_int_magn to grid_ops and implemented calc_int_vol.

## 0.79:
* When using Richardson extrapolation, the EV's for different levels are not displayed in log.
* Implemented higher order expression (only) for d/dx in EV problem, to be used with norm_disc_style 2.
* Changed the name from PB3D_PP to PB3D_POST.

## 0.80:
* Fixed a bug in divide_X_grid for a high number of processes which would cause X_limits(2) to rise above n_r_X.
* Added vacuum module vac, which will be used to calculate the vacuum response. Temporarily, this is set to zero.
* Reordered the SLEPC routines to accomodate iteration for inverse.
* Rewrote the SLEPC routines with the updated theory taking into account the plasma edge.
* Fixed a bug in trim_grid when non-overlapping grids were used.
* Now no more ghost points are needed for the calculations. However, for plotting and writing output a ghost region of one is kept.
* HDF5 is now the only possibility, but some outputs are given with GNUPlot.
* draw_GP and draw_GP_animated now have the possibility to draw decoupled 3D plots using draw_dim which replaces is_2D.
* Introduced parameter GP_max_size, above which GNUPlot is not used.
* Changed plot_jq to plot_resonance. It also outputs the resonating flux surfaces in 2D HDF5.
* Harmonic plot now also plots the resonatic flux surfaces in 2D.

## 0.81:
* Fixed inconsistencies in plot_X_vec and decompose_energy: the numerical derivatives need two-sided ghost region.
* Introduced symmetric ghost regions in PB3D_POST with width ghost_width_POST, as this is needed for the numerical derivatives of X.
* Accompanying this, trim_grid is extended with shift_grid.
* Fixed bug in calc_XYZ_real where flux_p_H was not divided by 2 pi.

## 0.82:
* Fixed an error in fill_mat for order 2: factor 12 in stead of 2.
* Fixed errors in calculating d_nz and o_nz.
* Changed norm_disc_style to norm_disc_ord, only referring to the order.
* Made fill_mat and set_BC generic for any order. The latter now employs set_left_BC and set_right_BC.

## 0.83:
* BC_style chooses which style for left and right boundary.
* Eigenvectors are now normalized with respect to their maximum modulus.

## 0.84:
* The output name of draw_GP(_animated) is now consistent with the system used for HDF5: a new variable draw_name is introduced.
* Energy decomposition output is now written in a file PB3D_out_EV.txt
* HDF5 output now shows less messages.
* PB3D_POST run script now takes PB3D_out.h5 from a directory specified and by default runs in that directory.

## 0.85:
* Reduced (eliminated?)  memory leaks by properly nullifying and deallocating pointers.
* Use_normalization cannot any more be passed to PB3D_POST, but is instead read from PB3D output.
* Corrected bug where mu_0 was forgotten in kappa_n.

## 0.86:
* Added new option 'debug_X_grid', which sets the perturbation grid equal to the equilibrium grid.
* Added new test 'test_diff' which looks for the effects of introducing artificial numerical diffusion.
* Implemented the 'debug_store_results' which checks whether the solution stored corresponds to the relative error returned.

## 0.87:
* Changed the routine 'divide_X_grid', which now returns the whole r_F instead of just grp_r_F.
* Changed the way the routine 'solve_EV_system_SLEPC' calculates the step size: It now makes use of grid_X%r_F.
* Corrected a VERY IMPORTANT bug by setting V_2 to - V_2.
* extremely uncorrect unstable spectrum largely corrected.

## 0.88:
* Introduced variables 'tol_slepc' and 'max_n_it_slepc', which are displayed before starting the solver.
* Fixed bugs considering BC's.
* Energies are not scaled by E_kin any more.

## 0.89:
* Improved makefile, adding support for quadrivium compilation.
* Changed things in 'output_ops' and 'sol_ops' that provoked a compilation bug on quadrivium: an array could not be used as a building block of a bigger array as in y = [x,5] where x = [1,2].
* Split the perturbation output in perturbation output and solution output, which are written to 2 different HDF5 groups.
* 'exp_ang_par_F' is now 'J_exp_ang_par_F' and contains the Jacobian as well.
* For PB3D, the metric variables are no longer interpolated on the field-aligned grid, as they are no longer necessary with the introduction of 'J_exp_ang_par_F'.
* Introduced a new run script, for qsub on quadrivium.

## 0.90:
* LAST VERSION OF OLD PB3D WITH PARALLELIZATION IN THE NORMAL COORDINATE.
* Small bug fixes.
* Started implementing memory check.

## 0.91:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS.
* Replaced the routines that calculate the possible derivatives of a certain order by a generic one.
* Implemented calc_derivs_1D_id which calculates the 1D index of derivatives. A table of these is stored in an array d for quick access.
* 'driver_rich.f90' has been renamed 'driver.f90' and the original 'driver.f90' was deleted.
* Introduced 'norm_disc_style' for eq, X and sol, so derivatives are done coherently.
* 'trim_grid' now optionally returns 'norm_id', the normal indices of the trimmed variables.
* The derivatives of VMEC outputs is now done in 'prepare_RZL' because the precision is not yet known before reading the input.
* 'calc_eq_r_range', 'divide_X_grid' and 'split_MPI' are now bundled in 'calc_norm_range'.
* 'reconstuct_PB3D' is also called in the pertubation phase, but only equilibrium variables are reconstructed.
* X variables now are all allocatable.
* The X jobs can either be found in a straightforward way, requiring the use of the mpirun option '--mca osc pt2pt', or in an intelligent way, but this requires an up-to-date version of openmpi. The flag 'lold_MPI' can be used to switch between the two versions.

## 0.92:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
* Removed the whole system using MPI to transmit which process does which X job. Now, it works using an external file as well as an external lock file.
* Fixed a bug when the wrong grid variables were written out.
* Split of the transformation of derivatives in F.
* Moved the plotting of flux quantities as these require knowledge of the derivative in F.
* Removed plot rank, so now rank is used everywhere.
* Renamed 'check_modes' to 'check_X_modes'.
* 'check_X_modes' and 'resonance_plot' do not use an X_type anymore.
* 'divide_X_jobs' now contains 'calc_memory'. Also it can be used in vector or tensor mode and in theory higher orders as well, if information about the number of variables at these orders is provided.
* The perturbation phase is to be split in two: a vector phase and a tensor phase. This is not yet implemented.

## 0.93:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
* 'prepare_X' has been replace by 'calc_X' which calculates the perturbation variables for a certain order.
* There are two orders of perturbation variables. These are calculated separately and the results of order 1 are read from disk when calculating those of order 2.
* 'print_HDF5_arrs' now detects whether it is an individual or collective call. Also, it checks whether a group already exists.
* PB3D files are not opened anymore like normal input nor are their 1D variables broadcasted: Every process opens the file on its own.
* Restored 'prog_style' to only 2 possibilities: 1 (PB3D) or 2 (POST)
* Split 'X_type' into 'X_1_type' (vectorial perturbation), 'X_2_type' (tensorial perturbation) and 'sol_type' (solution).
* Code working up to calcalation of tensorial perturbation quantities.

## 0.94:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
* HDF5 variables are now bundled in 'HDF5_vars'.
* 'print_output_sol' is now situated in 'sol_ops'.
* 'fourier_ops' is now called 'fourier'.
* 'read_PB3D' takes optional mode information which is passed to 'read_HDF5_arrs' which takes optional string arguments.
* there are now 5 types of variables--misc, eq, X_1, X_2 and sol--that can be printed and reconstructed independently.
* 'open_output' now also writes miscelleaneous variables to HDF5 file.
* Fixed bug in 'read_HDF5_arrs' where head group was not closed which causd problems for consequent reads.
* 'c' now does its job also for submatrices, when the limits are provided.
* Tensorial perturbation variables are now correctly calculated up to KV, PV.

## 0.95:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
* Matching modes now lead to less reading of variables.
* File and variable names have changed to consistent 'eq', 'X', 'sol' and 'POST'.
* Solution driver template introduced.
* In the perturbation driver, the only thing missing is the integration along the magnetic field lines.
* In the solution driver evetything needs to be pieced together.

## 0.96:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
* 'interp_HEL_on_grid' now overwrites input quantities if no outputs are provided.
* This way the doubling of the memory requirement for HELENA is avoided.
* The calculation of magnetic integrals now happens completely in 'X_ops'.
* 'PB3D_type' is not used anymore; Instead, individual variables are used.
* Perturbation driver has been completed.

## 0.97:
* NON-USABLE VERSION: SWITCHING BETWEEN PARALLELIZATIONS. 
* Updated 'X_ops' to not use 'PB3D_type'.
* Implemented the X_2 versions of the routines that read and reconstruct PB3D variables.
* Changed 'perturbation grid' to 'solution grid'.
* 'setup_and_calc_grid_B' has been renamed to 'setup_and_calc_grid_eq_B' and 'setup_and_calc_grid_sol' has been implemented.
* Completed the solution driver routine.
* Updated SLEPC routines to handle the solution type.
* Fixed a bug in 'trim_grid' where the limits of the trimmed grid were not updated.
* Introduced 'wait_file' which makes use of a lock file and used this in the writing of HDF5 files.

## 0.98:
* FIRST NEW USABLE VERSION: NEW PARALLELIZATION STRATEGY.
* The miscellaneous variables are now read in 'open_output' for POST, just as they are written there for PB3D.
* The PB3D operations have to be initialized before using them. This sets the equilibrium style.
* Modified 'broadcast_input_vars' to also broadcast miscellaneous variables in POST.
* The maximum fluxes are now written in the equilibrium variables.

## 0.99:
* Implemented a plot of the comparison between the Eigenvalue and the energy fraction in 'plot_X_val_comp'.
* Bug fixes in calculation of ranges and 'no_plots' for POST.
* Bug fix in plotting of resonant surfaces.
* 'norm_disc_prec_sol' is now used in PB3D to indicate the discretization precision of the matrices A and B.
* SLEPC routines now need a trimmed solution grid, and an untrimmed solution grid is nowhere necessary in PB3D.

## 1.00:
* FIRST WORKING VERSION
* Bug fixed where Hermitian conjugate was stored in stead of the real matrices in SLEPC.
* Checked consistency of results when 'intuitive' version of plasma potential energy is used.

## 1.01:
* UNUSABLE VERSION: See last comment.
* Plots can now have an optional different size.
* Corrected bug with wrong sign in calculation of jacobian in HELENA coordinates.
* Merged contents fourier with VMEC modules.
* Fixed bugs in the HELENA interpolation.
* Changed 'read_HEL' for asymmetric cases to duplicate information for theta equal to 0 into theta equal to 2pi.
* End-points are now included for the parallel angle, to simplify working with trapezoidal integration. This eliminates the discrepancies between results for different parallel intervals.
* N.B.: NOT VALID for some other integration rules!
* Added a test folder to be used for FORTRAN tests.
* Started implementing the transformation of the metric factors. This is NOT finished and will probably be abandoned

## 1.02:
* Metric variables are interpolated still for HELENA, but in a different, more clever way. Not yet implemented for toroidal flux.
* Bounds are now checked.
* Fixed and finished implementing the comparison between U and U_inf.
* Extended 'debug_store_results' so that it displays X*AX and X*BX separately.
* Extended debugging of calculation of integral. An alternative calculation for orthogonal computational grids shows identical results.

## 1.03:
* Fixed a bug in the calculation of the volume integrals: ghost region was needed.
* Consistency of volume integral with both repeated trapezoidal integration and repeated simple integration shown.
* Introduced normalization style 'norm_style' that now also allows for normalization of only the normal component of the perturbation.
* Various bugfixes when no resonant surfaces are found or when solutions are removed.
* Solution vectors are now normalized in store loop, to ensure consistency of X*VX step_size and the integrals of POST.
* Fixed a bug in 'calc_memory' where overflow occured.

## 1.04:
* UNUSABLE VERSION: CHECKING INDIVIDUAL TERMS WITH ENERGY DECOMPOSITION.
* A more detailed grid for X quantities is needed!
* Cleaned up the code when 'ldebug' is not used.
* Implemented debugging of the perturbation driver so that the vectorial instead of the tensorial perturbation variables are interpolated.
* Implemented the option to choose order of geodesic perturbation.
* Implemented a test to plot |X|^2 after solving the Eigenvalue problem and compare this with the energy reconstruction.
* Fixed a bug where the transpose of a matrix should be used when feeding into SLEPC.

## 1.05:
* ONLY PB3D USUABLE.
* POST has to be adapted still.
* Reintroduced all terms in the energy equations.
* New module for Richardson extrapolation.
* Extrapolation is now done in main program, as X variables have to recalculated on a different grid for each Richardson level.
* Involved routines 'calc_U', 'calc_KV' and 'calc_PV' have been rewritten. 'calc_U' is still not tested for VMEC.
* Updated the user output considering Richardson extrapolation in general.
* Failure to read input file now results in error.
* 'test_max_mem' now indicates whether memory is to be tested or not.

## 1.06:
* FIRST USUABLE VERSION WITH CORRECT ENERGY DECOMPOSITION: PERTURBATION VARIABLES NOW CORRECTLY TABULATED IN OWN GRID.
* Adapted POST to new X storage convention as well.
* Fixed some bugs considering Richardson extrapolation.
* POST now by default loads the highest found Richardson level, but this can be overriden using 'PB3D_rich_lvl'.
* Bugfixes considering the new tabulation of perturbation variables.
* 'debug_store_results' is now part of the standard programme.
* HELENA always uses normalization, so 'use_normalization' is always set to true.

## 1.07:
* FIRST VERSION WITH FAST OPTION. It seems to work for cbm18a if 'max_r_sol' is not larger than 0.9.
* Introduced variable 'X_style', which allows the user to choose between prescribed modes (earlier default) and a fast version.
* 'check_X_modes' now has alternative actions for X style 2.
* The original versions of 'get_suffix', 'set_nn_mod' and 'is_necessary_X' do not exist any more: Everything works with X limits now.
* Split some routines of PB3D_vars into new module PB3D_utilities.
* Split some routines of X_vars into new module X_utilities.
* Split some routines of SLEPC_ops into new module SLEPC_utilities.
* Split some routines of input_ops into new module input_utilities.
* Split some routines of grid_ops into new module grid_utilities.
* Split some routines of sol_ops into new module sol_utilities, which can also translate between local and total solution vectors.
* Derivatives of the solution vectors have to be done on total variables.
* Suffixes to names are now integer variables, not characters.
* The limits of the modes and their values now have to be calculated from the equilibrium and have to be recalculated using 'setup_nm_X'.
* Fixed a bug in 'interp_fun'.
* Improved debugging of 'calc_zero_NR'.
* Fixed a bug in 'calc_res_surf' for dfun.
* Introduced local versus total mode numbers, which may differ for X style 2.
* 'plot_X_vec' has been modified as well as 'resonance_plot' to correctly implement the difference between local and total mode numbers.
* Fixed inconsistency with naming 'X_vec' and 'X_val', etc. 'sol' Is now used.

## 1.08:
* NOT USABLE FOR VMEC: ERRORS PRESENT!
* Implemented routine to get info about file, though unused for now.
* New flag 'no_execute_command_line' that disables the execution of command line, as this sometimes causes a freeze, for unknown reasons.
* Identified bug when reading VMEC with 'lrfp' flag. Now, netcdf output is used instead, as it does not have this problem.
* Fixed some bugs in 'reconstruct_PB3D' considering the presence of field-aligned grids. These are now automatically set for both equilibrium styles.
* Fixed a bug in the first derivative of 'flux_p_E' for HELENA: factor 2pi is now taken into account.
* Changed the storage convention for 'PV_int' and 'KV_int'.
* Generalized and automatized the 'calc_deriv' routines and put them in 'grid_utilities'; They can take any order and precision now.
* Cleaned up 'calc_flux_q'.
* Split HELENA in HELENA_ops and HELENA_vars.
* Added interpolation capability similar to the derivation capability.
* Both interpolation and derivation use 'apply_disc' to apply the discretization operator.
* 'interp_fun' is still available to be used only when the interpolations are not repetitive and linear interpolation suffices.

## 1.09:
* FIRST VERSION THAT WORKS FOR AXISYMMETRIC VMEC. However, post-processing is not working well, and comparision with HELENA not yet perfect.
* Normalization factors 'R_0', 'B_0', 'pres_0', 'psi_0' can now be user-provided, just like 'rho_0', but NO checking for consistency!
* Fixed a bug in 'calc_interp_data' where the bounds of 'x' were not respected.
* Fixed bug in Post-processing for VMEC input where the plot output is calculated.
* Removed 'test_Dg_E' as it does not make sense: In general the equilibrium grid is not straight.
* Removed the 'repack' routine and started using the natural VMEC structure.
* Fixed an error in the testing of B where only components in equilibrium coordinates were not properly translated to flux coordinates.
* The structure of 'SLEPC_ops' has been slightly rewritten for the future accomodation of shell matrices.
* Fixed bug in 'calc_U' where the part ~ n_frac was not used in the calculation of DU_1 and DU_2, as well as some erroneous indices.
* Moved some of the routines from 'calc_eq' to 'calc_met'. In fact, these routine structures should be rethought, possibly merged.

## 1.10:
* VMEC IS NOT WORKING WELL YET.
* Added restart functionality within PB3D using 'rich_restart_lvl' to enable for Richardson extrapolation to continue (2...max_it_rich) or to skip the pre-perturbation phase (1). If it is 0, no restart is done.
* The broadcast routines now only pass user options. All the rest is handled by HDF5 output and input.
* Fixed confusion about the normalization: Now HELENA uses MISHKA normalization by default.
* The miscellaneous output variables are merged with the input output variables, i.e. output variables due to the input of equilibrium code.
* The system of reading and then reconstructing variables has been replaced by a direct one. 'PB3D_vars' ceased to exist.
* Split 'rich' module in two parts: 'rich_vars' and 'rich_ops'.
* The deallocation of equilibrium input (VMEC or HELENA) now happens through 'dealloc_in'.
* Standardized the outline of driver routines, by making them start by reconstructing the necessary PB3D output and finish by deallocating.

## 1.11:
* Bug fixed: The input phase now passes on a trimmed range. The 'calc_norm_range' routine has been extended accordingly.
* 'misc_eq', 'misc_eq_V' and 'misc_eq_H' have been renamed to 'misc_in', 'misc_in_V' and 'misc_in_H'. 'n_r_in' has been added.
* 'max_flux' now only refers to the maximum flux in the coordinate system used (i.e. poloidal or toroidal).
* Fixed a bug in 'calc_res_surf' where instead of the rotational transform the safety factor was used.
* 'flux_p_V' and 'Dflux_p_V' are now calculated inside 'read_VMEC' as the full input grid is needed.
* 'calc_XYZ' now requires the equilibrium grid for the routine to know how the input variables are tabulated.
* 'coord_E2F' and 'coord_F2E' do not need the equilibrium variables any more. Also, 'calc_loc_r' has been deleted.
* VMEC now also uses magnetic field on axis as normalization, saved in 'B_0_V' in VMEC.

## 1.12:
* VMEC IS WORKING AND HAS BEEN TESTED FOR AXISYMMETRIC CBM18A.
* HOWEVER, THE PROCESS OF CHANGING RICHARDSON EXTRAPOLATION IS NOT COMPLETED YET.
* 'eq_type' and 'met_type' have been merged and split differently again in flux equilibrium variables and metric equilibrium variables.
* 'eq_ops' and 'met_ops' have been merged and some routines have been split off in 'eq_utilities'.
* 'interp_HEL_on_grid' interpolates metric equilibrium variables only, though flux equilibrium variables are needed to do it.
* 'calc_flux_q' has been absorbed in 'calc_eq_1'.
* Moved 'print_output_in' and 'read_eq' to input_ops and 'dealloc_in' to input_utilities.
* 'rich_restart_lvl' now can take on values 1..max_lvl_rich, not 0.
* Streamlined output from reconstruct_PB3D routines.

## 1.13:
* X_2 now holds only 6 variables, and is reused for the field-averaged variables as well.
* Got rid of unnecessary and cluttering default case selection for variables that are checked in initialization.
* Print_output routines need data set names provided, similar for reconstruction.
* Now shell commands are saved in a log file that can be run afterwards, which is useful in combination with '--no_execute_command_line'.
* 'rich_info_short' now always returns the Richardson level in the output.
* Fixed bug in 'find_max_lvl_rich' where groups were not properly closed.

## 1.14:
* BOTH VMEC AND HELENA WORK FINE WITH RICHARDSON EXTRAPOLATION
* Passing the Richardson level to be used in the data name for HDF5 output is now done using an optional integer.
* 'rich_info_short' has been removed.
* For VMEC, variables saved in HDF5 from different levels are now combined correctly in the reconstruct_PB3D routines through 'tot_rich'.
* 'tol_SLEPC' is now an array to be passed for every Richardson level, or is set intelligently by default.
* The guess for different Richardson levels is set correctly again, keeping in mind that every level has the same nr. of normal points.
* Migrated away completely from 'interp_fun'. The functionality is coverd by 'setup_interp_data' and 'apply_disc'.
* The user can now provide the default relaxation factor for Newton-Rhapson.
* Tweaked the calculation of magnetic field lines to make it faster
* HDF5 datasets can be overwritten when Richardson restart is used.
* Richardson variables are not written to HDF5 any more and reconstructed for Richardson restart but are set up from solutions.

## 1.15:
* IMPROVED MEMORY USAGE. CURRENTLY BEHAVIOR IS NOT COMPLETELY UNDERSTOOD BUT MEMORY USE SEEMS TO BE LIMITED.
* Removed splines.
* Fixed memory leaks concerning 'disc', 'grid', 'eq_1', 'eq_2', 'X_1', 'X_2' and 'sol' variables and added a check.
* Fixed memory leaks concerning the printing of variables to HDF5: previously there was no deallocation of var_1D.
* Added color to the output, using the FOUL module. There is now the option 'warning' to subroutine 'writo'.
* Added option '--mem_usage' to print information about memory usage at the end of every message.
* Simplified 'create_grid'.

## 1.16:
* Fixed a bug in POST where not the total Richardson variables were taken.
* Fixed a bug in 'reconstruct_PB3D_X_1' where the variable names were read incorrectly.
* 'calc_XUQ' does not need 'grid_sol' any more, as stated in the header.
* Calculation of extended plot grids in POST is now much more economical.
* Fixed a bug for negative normal coordinates, but the code has not been debugged properly for this!
* min and max of theta and zeta_plot is now an input variable.
* Moved 'grid_plot_real' to the magnetic integral phase of driver_X because it needs the full field-aligned equilibrium grid.
* Fixed bug concerning 'test_p': The test should be done after F derivatives are calculated from E derivatives.
* Command line is not any more executed by default.
* Introduced multiple tries for Newton-Rhapson, with different relaxation factors.

## 1.17:
* Fixed a bug in the calculation of the normal range for high discretization orders.
* Fixed a bug in the interpolation routines for high discretization orders.
* The routines that adapt input variables now have correct error handling.
* New input flag 'POST_style' that allows choice between POST output on extended grid (1) or field-aligned grid (2).
* Reorganized the POST driver structure.
* User can now use input variable 'slab_plots' to optionally generate slab plots.
* Fixed bug in 'plot_HDF5' where plot was wrongly identified as having poloidal symmetry while it was just a slab plot.
* For debugging, the angular coordinates in POST plotting can be swapped for style 2 (field-aligned).

## 1.18:
* FIRST WORKING 3-D VERSION, though there might be bugs.
* The swapping of angular coordinates in POST plotting is now runtime option '--swap_angles', not only for debug.
* Fixed bug where alpha was twice multiplied by pi in output.
* Fixed bug where SLEPC tolerance was wrongly set.
* The phase of the Eigenvectors is also plotted now.
* Fixed bug in the calculation of trigonometric factors, where nfp was forgotten.
* Improved the routine that calculates numerical interpolation by allowing it to switch to lower orders.

## 1.19:
* Both the real and the imaginary part, as well as the phase, of the output vectors are now plotted in HDF5.
* Fixed a bug with the radial coordinate of the plots for POST when multiple processes are used.
* Replaced the workings of 'setup_interp_data' to use Barycentric Lagrangian polynomials.
* Renamed 'utilities' to 'num_utilities'.
* Fine-tuned usage of 'calc_zero_NR' when calculating the magnetic field lines by using a better guess.
* Removed the restriction of the 3D decoupled GNUPlot plot of the modes as by default nothing gets plotted.

## 1.20:
* Updated the manner in which Richardson extrapolation is done by reusing the previous magnetic integrals, and not the values.
* The values of 'X_2' are not any more written to HDF5, only the magnetic integrals.
* For HELENA, the grid is also halved for higher Richardson levels.
* Apart from the trapezoidal rule, also Simpson's 3/8 rule can be used for magnetic integrals.
* Renamed 'plot_grid' and 'plot_grid_real' to 'plot_magn_grid' and 'magn_grid_plot' for consistency.
* Moved 'magn_grid_plot' back to the equilibrium driver. The full grid is reconstructed specially.
* Reorganized perturbation driver: It consistts now of 3 subdrivers.

## 1.21:
* Fixed a bug in the calculation of n and m for toroidal flux.
* Fixed a bug in the calculation of the F derivatives of the flux quantities, which was ugly but did not have many consequences.
* Removed unnecessary reconstruction of equilibrium variables in solution driver.
* Replaced pointers in eq_types by allocatables.
* Rewrote the initialization and deallocation routines for eq, X and sol variables.

## 1.22:
* NOW THE MEMORY LIMITS ARE CORRECTLY HANDLED BY INTRODUCING EQUILIBRIUM JOBS.
* Equilibrium jobs now handle a subset of the parallel grid.
* Added option to deallocate unused variables on the fly in 'calc_eq_2'. Using this, the memory spikes for the equilibrium phase are lower.
* Fixed bug where memory info was printed before the file existed.
* Moved 'divide_X_jobs' to 'X_utilities'.
* Created an equilibrium version of 'divide_X_jobs'.
* Removed the routine 'test_max_mem'.
* Input variable 'max_mem_per_proc' has been renamed to 'max_tot_mem_per_proc', to distinguish from internal variable 'max_X_mem_per_proc'.
* 'read_HDF5_arrs' now has an array version as well, which uses the individual version.
* The same for 'retrieve_var_1D_id'.
* Using 'minim_output', the output file size can be minimized, by not saving eq_2 and X_1 variables between Richardson levels.
* 'tol_SLEPC' for next Richardson levels is adapted to the maximum relative Richardson error.
* Memory information is extended with the limits.

## 1.23:
* FIRST OPTIMIZED VERSION TO FIND MAGNETIC FIELD LINES, BUT RESULTS ARE NOT COMPLETELY EQUAL TO BEFORE, AND POST DOES NOT WORK.
* Removed shell matrix things, as they are not yet implemented.
* Fixed a bug in 'insert_block_mat' where 'block_loc' was used erroneously.
* Clean up a bit the SLEPC routines, all using n_mod_X now, fixed some minor memory leaks.
* Time information in memory info now comes from MPI_Wtime.
* 'fourier2real' now also has a version that does not make use of trigonometric factors, but of theta and zeta directly.
* Improved debug of 'calc_ang_grid_eq_B' by checking whether F variables are recovered.
* There is now also a 3D equivalent of 'calc_zero_NR'. It is used in 'coord_F2E', so now the magnetic field lines are calculated faster.
* Implemented Zhang's method for root-finding, which is now used to calculate the resonant flux surfaces.

## 1.24:
* FIRST COMPLETE PB3D VERSION. SUSPECTED PROBLEMS WITH 1.23 WERE DUE TO LOW 'tol_SLEPC'.
* Minimal outputs are also handled by POST now.
* Fixed a bug where incompatibility between HELENA and minimal output was ignored.
* Improved the handling of 'tol_SLEPC'.

## 1.25:
* Added files from Alpha study.
* The solution driver now does not need the equilbrium grid any more.
* Split off some of the procedures in 'HDF5_ops' into 'HDF5_utilities'.
* Implemented 'set_1D_vars' which sets the hyperslab of the 1D equivalent of a variable in multiple dimensions and/or chunk variables.
* Fixed a bug in the conversion of half to full mesh of VMEC variables.
* 'norm_style' is now called 'K_style' and 'norm_style' is used to set the style of normalization (e.g. MISHKA, COBRA, ...).
* Updated the run scripts to use Dr. Memory in stead of Valgrind.
* 'get_suffix' is renamed to 'get_sec_ind_tot' as suffixes are not necessary any more, but total secondary index are.
* Implemented 'get_sec_X_range' that seeks a contiguous range for tensorial perturbation variables.
* 'read_HFD5_arrs' is now 'read_HDF5_arr' and returns just one 1-D variable.
* 'retrieve_var_1D' is not necessary any more and has been removed.
* Improved the MISHKA normalization by using the true value of B on axis, extrapolating half-mesh in VMEC. This changes T_0 slightly.
* In the PB3D reconstruction routines, the normal limits are not passed any more, as this information is encoded in the grid already.

## 1.26:
* Added a command-line variable 'jump_to_sol' that can be used to jump straight to solution driver for first Richardson level.
* The integrated tensorial perturbation are deallocated in the SLEPC routines, instead of in the solution driver, to save memory.
* The same is true for the previous solution variables.
* Fixed bugs in the PB3D reconstruction routines, where the equilibrium jobs were not appropriately calculated for every Richardson level.
* Fixed bug where individual write was neglected for multiple processes.
* Changed R_H and Z_H at first normal position to raxis, respectively 0 for HELENA.

## 1.27:
* Fixed some confusion about parallel acces and transfer property lists in HDF5.
* Individual writes of PB3D variables now use standard I/O driver so that more than 2GB can be written independently.
* Changed structure of the the X_2 part of the perturbation driver for HELENA: The X_1 variables are written to HDF5, not the X_2.
* Slightly improved user output for perturbation driver.

## 1.28:
* The HDF5 routines do not work if only one process doing HDF5 output. So now there is a duplication in the printing of grids and solutions.
* Fixed a bug concerning the reconstruction of full variables for Richardson levels greater than 1, which affects POST only.
* The Energy Reconstruction output now carries the Richardson level '_E'.
* Fixed a bug in the calculation of the normal ranges for POST, where the sign of the tolerance was taken wrongly.

## 1.29:
* Fixed a bug where the reading of X variables where the local parallel limits were erroneously set. Now they are set to -1.
* In the reading of HDF5 variables, if there is a negative upper local limit, the total limit is taken.
* Reading HDF5 files now also checks for the lock file that indicates the file might be being written at the time.
* Fixed a bug in the writing and reading of BC_style, which is an array of size 2, not a scalar.
* Cleared confusion about COBRA normalization and fixed the erroneous situation that called for modification of the equations.

## 1.30:
* 'minim_output' is now also possible for POST.
* Tweaked setting up of guess using now a tolerance.
* Fixed a bug considering plotting the grid from HELENA equilibria.
* Fixed a bug determining the normal range in POST.
* Changed the output of the modes at midplane to use lines, not points.
* Updated 'write_flux_q_in_file_for_VMEC' to also give the Fourier coefficients of the boundary shape.
* read_HEL now calculates the toroidal flux on the full grid.
* Moved read_HEL to HELENA_ops, as it needs grid utilities.

## 1.31:
* Created new script to run PB3D, making use of cases, that bundles all the previous scripts.
* run_PB3D.sh can also run on Uranus, making use of qsub. Selection of machine is automatical but new machines have to be added manually.
* Fixed some memory leaks.
* 'lold_MPI' is no longer available.

## 1.32:
* Makefile now uses shared libraries.
* Fixed bug in the run scripts that used kb instead of mb.
* Lock files are now explicitely removed at initialization in the code itself.
* Improved the run scripts further by making use of common functions for PB3D and POST in run_aux.sh.
* The run scripts now have the option to run for multiple variations of an input file, making use of an additional optional flag and file.
* Fixed a bug when running HELENA in release mode considering the metric variables h_H.
* Removed the need for lockfiles, as everything is now performed internally using a mutex system made with MPI.
* The mutex system can be debugged using a debug flag in MPI_utilities.
* For large numbers of procs, the code is slow as every HDF5 read and write is currently isolated using mutex. This is necessary as it is not possible to read from a file that is being written to, even though multiple procs can indeed access the same file at the same time. This needs to be resolved; either using SRMW from version 1.10, or another mechanism.

## 1.33:
* Renamed 'str_ops' to 'str_utilities' for consistency.
* Split num_utilities into num_utilities and num_ops, where the latter requires output_ops.
* Bubble sort implemented in num_utilities.
* Renamed 'mutex' to lock.
* Implemented tests for the lock system.
* Completely redesigned the lock system.
* The lock system now supports blocking as well as non-blocking locks, which is good for HDF5 performance with many reads and few writes.
* Sometimes opening HDF5 files still fails. Multiple attempts will therefore be performed.
* The timing functions now use system_clock with integer of kind 8, which is precise.

## 1.34:
* Fixed a bug where no NB file lock was requested for probe_HDF5_group.
* Fixed a bug where jump_to_sol no longer worked because n_X and m_X were not set. A part of the X driver is now also run.

## 1.35:
* Reorganized the external library structure (foul and fftpack), using better compilation flags as well.
* Fixed a bug in con2dis_reg, where no correct range was found for negative values.
* Added the possibility to perform trigonometric interpolation, but only for an even order.
* Added tests for interpolation routines.
* Fixed the calculation of the fourier modes for exportation of HELENA to VMEC. Simulations have to be rerun.
* Export of HELENA to VMEC now happens through the flag '--export_HEL'.
* Procedure merge_GP was deleted.
* GNUPlot 'GP' procedures are renamed to external 'ex' procedures.
* Array versions of draw_ex and draw_aminated_ex, corresponding to multi-file versions, are deleted.
* The external plotting now takes multiple draw options.
* External animation is only possible for 2D.
* Added external plotting possibility: Bokeh for 2D and Mayavi2 for 3D.
* There is a problem when Mayavi plots are generated through system calls in the program, but they do work by calling them afterwards.

## 1.36:
* Fixed bugs with the drawing options.
* Multiple 3D plots are again supported, for GNUPlot and Mayavi.
* Unified external program drawing into one procedure.
* A HELENA equilibrium can now be perturbed before exporting to VMEC format.
* The modification of the equilibrium can be provided in a file as well. Examples are given in 'pert_*.dat'.

## 1.37:
* The plotting variables of HELENA equilibria for VMEC export can now be adapted.
* HELENA export now outputs a VMEC input file that can be further adapted.

## 1.38:
* VERSION NOT INTENDED FOR USE.
* Many changes to improve the resilience for input / output errors in order to work on the ITER cluster.
* New precompile include that allows to run any command (typically I/O) to be repeated a number of times, using the CPP flag 'lrIO' for 'resilient I/O'. This is going to be removed in the next commit, this commit serves solely to store the idea.
* Bug fixes.
* Compatibility fixes with older compilers and INTEL compilers.
* Support for INTEL compiler (checked with 12.0.1), using a wrapper include file and a compile flag 'lwith_intel', as opposed to 'lwith_gnu'.
* Workaround for recursive function 'derivs' that did not work with INTEL (https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/270127).
* 'search_file' has been deleted, as its functionality can be done using open.
* Improved run scripts, with better support for TORQUE and SLURM clusters. This solves the I/O problems.
* New script for extraction of results from array jobs.
* New script for inspection of SLURM jobs.
* Cleaned stellinstall and added it to this repository.

## 1.39:
* Removed the 'lrIO' precompiler option.
* '--minim_output' is set by the run script by default. Remove it in 'run_PB3D.sh' if you don't want this.
* stellinstall is now configured for 1. XPS-L501X, 2. ITER.
* Updated run scripts and bundled everything into 'run.sh'.

## 1.40:
* Fixed bug in output routines where no spaces were left in output data files.
* Fixed bug in run script variable substitution where min_r_sol was replaced instead of n_r_sol.
* Add global multiplication factor to HELENA modification and export.
* Changed the run script slightly.
* NUFFT now happens in 'grid_utilities'.
* The routine HELENA export now also plots the approximate proportionality constant between the toroidal field ripple and the plasma displacement, as well as the corresponding Fourier series.
* HELENA export now also is able to automatically use this for the edge modification shape.

## 1.41:
* Introduced variable 'pert_mult_factor_POST', which is used to deterimine by how much to perturb the X, Y and Z values of the plot grid in POST. It works only for VMEC on an extended grid (POST style 1).
* In 'calc_pert_cart_comp' the Cartesian components of the perturbation are calculated. If 'pert_mult_factor_POST' is not zero, these are then used to perturb the position.

## 1.42:
* Fixed a bug in the calculation of the Cartesian components of perturbation where the trimmed grid was erroneously used.
* Fixed a bug when checking the input normal range, erroneously making use of the number of points in the equilibrium, which is not yet known.
* The solution normal range can now be chosen for POST. Note that this obviously affects the energy reconstruction.
* In the POST driver, the nm variables are now setup for the entire ranges.

## 1.43:
* Fixed a bug when 'calc_tot_sol_vec' did not work for partial normal grids.
* Added script that can concatenate HDF5 files.

## 1.44:
* Changed institution to ITER Organization.
* Fixed a bug where uninitialized pointers were used in the POST driver for non-full outputs. Now, the harmonics are plot separately.
* 'slab_plots' is now an integer variable 'slab_plots_style' that can take value 2 as well, which stands for slab plots with the angular coordinates wrapped to the fundamental intervals.
* Improved the run script for the ITER cluster, making use now of two auxiliary scripts 'gen_node_file.sh' and 'get_disk_space.sh'.

## 1.45:
* The run script does not copy the HDF5 output file for POST on the ITER cluster, but instead creates a symbolic link.
* Updated some comments in the code, for example concerning the wrong statement that HDF5 data read has to be contiguous. This is not any more the case.
* The routines that call 'divide_eq_jobs' now need to include derivatives the array size.
* There is a new function 'copy_grid' that allows the user to copy a grid to a new, unallocated grid, taking a possible subset of variables.
* Part of the tasks of the POST driver is now done in an initalisation routine 'init_POST'.
* The other parts are done for multiple equilibrium jobs.
* These jobs are divided according to memory requirements.
* 'setup_par_id' has been generalized and can also output the indices in HDF5 memory.
* In HDF5 array reading, the indices of the dummy variables now start at 1 always. It is the job of the initialization of the various variables (e.g. eq_1, eq_2, ...) to set the range limits correct (e.g. starting at 0 for derivatives, ...).
* This is an untested version.

## 1.46:
* In 'print_ex_2D' and 'print_ex_3D', x (and y) variables can now be provided for 1 plot, and are then copied for the possible others.
* Implemented a test procedure 'test_read_HDF5_subset' for the subset reading capacities.
* This routine can test for the correct reading of subsets, with or without normally divided grids.
* Some bugs have been fixed.
* Tests have been concluded.
* There appears to be a problem when reading the variables from HELENA equilibria for testing: Performance drops very much. This needs to be investigated.

## 1.47:
* The problem was due to too small chunk sizes for the variables that were written by multiple processes, not just for HELENA.
* Chunk sizes are now calculated different, by multiplying the total sizes for increasing dimensions until reaching more than the minimum requested (10kb set for now).
* Chunk caching is not done any more, so 'a_plist_id' from 'set_1D_vars' is not available any more.
* Added support for scalasca in the compilation, gprof is deprecated as it does not work properly, even for only 1 process.
* Added a file 'scorep.fil' with common scorep filters for the user-input routines that take up time. Use with -f scorep.fill when analyzing with scalasca.

## 1.48:
* Fixed a bug in the reading of subsets for X_1.
* Fixed some bugs for driver_POST.
* For POST, the magnetic integration style is now set to 1, no matter what PB3D did, as the integration is done by volume, using the trapezoidal rule always.
* Changed the structure of driver_POST a little bit, as the integrations of energy have to be added up for different parallel jobs.
* The energy reconstruction now takes into account multiple equilibrium parallel jobs.
* There is still one HDF5 output per equilibrium parallel job. This could be changed in the future. In Paraview, 'group datasets' can be used.
* A note considering the implementation: The division for POST is always done in the first coordinate, which commonly corresponds to theta when poloidal flux is used as normal coordinate. This makes the individual plots a bit weird when not viewed in their entirety. This should not be an issue if all the jobs are plotted together.

## 1.49:
* Symmetry type can now be forced in HDF5 output and is also written to the file xdmf in an information element.
* There is now an option 'cont_plot' that indicates that the plot in HDF5 is a continuation of a previously started plot. This can be used to (over-)write HDF5 data, for example when dividing into jobs in POST.
* POST now makes use of continued plots so that the plot outputs are single files.

## 1.50:
* Version 1.49 suffered from an important flaw: It only worked when there were no multiple equilibrium jobs in PB3D.
* The system of equilibrium jobs is now overhauled: There is only one group per Richardson level.
* Many routines have been simplified.
* 'eq_job' has been removed in many routines.
* 'read_HDF5_arr' now only has one version, and 'conv_1D2ND' only has individual versions.

## 1.51:
* The full output grids that were used in POST to later extract subsets for the different equilibrium jobs are not used any more. Everything is done for the subset directly.
* The option 'minim_output' does not exist any more, as it is the only option.
* To get field-aligned output, the variables are calculated again in POST, without much penalty of time.
* 'slab_plots_style' is now 'plot_grid_style' and admits a new style 3 that corresponds to a straight cylinder by unwrapping the torus.

## 1.52:
* Metric equilibrium variables are not written to file any more for VMEC, as this is slow for multiple equilibrium parallel jobs.
* Both flux and metric equilibrium variables are passed on through the drivers for PB3D, and not deallocated and read.
* The flux equilibrium variables are deallocated at the end of the perturbation driver of every parallel job.
* This happens only to the metric equilibrium variables for VMEC. For HELENA, they are kept until the end.
* Vectorial perturbation variables are still written to HDF5 output, but for HELENA this is the full HDF5 output file, and for VMEC, this is the temporary HDF5 output file, so no parallel subset is needed. For HELENA, there are no equilibrium jobs for the first Richardson level, so nothing was changed there.
* These changes aleviate the problem that existed with big simulations with a lot of equilibrium jobs. In these cases, the whole 1-D variable system does not work satisfactorily, as the storage order in these 1D variables is unfavorable. Note that if these systems are to be used at some point, this should be looked at.
* Simplified perturbation driver somewhat.

## 1.53:
* The division of in parallel jobs now takes into account the fact that we have to be able to calculate the perturbation variables, which takes some memory as well.
* A bit less output for perturbation driver.
* Fixed a bug in the reconstruction of PB3D input variables.
* INTEL has a strange bug in 'calc_E', which is solved by setting an allocated array to zero.
* Fixed a bug that had not yet appeared where loc_n_r was used in the calculation of the X jobs size. This must be something that is not dependent on the process.
* The first steps are taken for the plotting of the magnetic field.

## 1.54:
* 'calc_XYZ_grid' now also optionally outputs R.
* Fixed a bug in 'calc_XYZ_grid' where the inverse toroidal variable was taken wrongly.
* Previously in XDMF output, the dimensions and number of elements were only written if they were larger than 1. This was removed as it seemed not to work with ParaView.
* 'plot_HDF5' can now also handle vector plots with col=4. The implementation could be generalized and made more elegant and robust.

## 1.55:
* Fixed a bug in the calculation of required memory.
* Fixed bugs in the plotting of B for HELENA.
* The routines in 'calc_B' have been generalized for any vector in 'calc_vec_comp' in grid_utilities and are called by 'B_plot'.
* 'calc_pert_cart_comp' should be replaced by 'calc_vec_comp' but for this, the calc_vec_comp needs to interpolate transformation matrices.

## 1.56:
* Extended 'calc_vec_comp' to include arbitrary grids.
* 'calc_pert_cart_comp' is not used any more.
* 'plot_HDF5' needs to be rewritten for general vectors, not using collection type 4, so that it can also take into account the possible projection of a vector on the symmetry plane. This is for the future.
* 'adapt_inoutput_POST' is removed, as also the solution grid can be perturbed for field-aligned grids as well.

## 1.57:
* 'plot_sol' is now split in two parts, 'plot_sol_xi' and 'plot_sol_Q' that can be used to indicate that the plasma perturbation and the magnetic perturbation have to be plotted.
* 'plot_kappa' can now be used to plot the curvature.
* Fixed bug in 'plot_sol_vec', where the normalization constants were not taken into account.
* The output routines now should produce output that has been transformed back to unnormalized values and so does 'calc_XYZ_grid'. 
* In 'plot_kappa' there is the debug option to show the center of gravity, as a curiosity.

## 1.58:
* Fixed some small bugs concerning the non-debug version.
* Implemented 'J_plot' that serves to plot the current, similar to 'B_plot' that plotted the magnetic field.

## 1.59:
* Fixed a bug for the tests on input variables. They were called too late.
* Added a normalization to be used with 'pert_mult_factor_POST' to provide X_0.

## 1.60:
* Intel makefile now specifies that the heap should be used for arrays larger than 100kB, to avoid overflowing the stack.
* 'dealloc_vars' is back in use for the metric equilibrium routines, and has been extended to 'broadcast_output_eq_2', where it is most critical.
* The initialization of equilibrium variables is done more carefully now, avoiding unnecessary E variables where possible.
* The calculation of memory in 'divide_eq_jobs' was wrong in assuming that the normal perturbation range would be divided. Also, total memory usage is not available any more.

## 1.61:
* Extensive rewrite of the whole overarching system with drivers.
* Variables come from the main program and are passed through the drivers. They are deallocated finally in 'stop_MPI'.
* Input variables are now saved throughout the whole program, after reading them in 'init_rich'.
* X jobs do not longer exist. They are not done in a large batch as before, but by iterating over the individual modes, possibly in blocks.
* The normal grids in the perturbation phase are now also divided and equal to the trimmed equilibrium grids.
* Implemented procedures to redistribute grids and equilibrium variables over a new normal range.
* The division in equilibrium jobs is now different: For PB3D HELENA only nchi parallel values are calculated and their interpolation does not happen in batch.
* 'tol_norm' is not used any more to calculate the normal extent of the input variables, so that they match with the grids used later.
* The maximum memory in total is now passed as input variable 'max_tot_mem'.
* Started usin Valgrind again, in stead of Dr. Memory.
* Fixed a bug where 'pert_mult_factor_POST' was wrongly used when it was zero.

## 1.62:
* Well-functioning version.
* Fixed bug: X_2 output is only written to HDF5 for last equilibrium job.
* Updated and cleaned up the job distribution procedures.
* Threw away many unnecessairy 'wait_MPI' commands.
* Updated the runtime parameters for SLEPC.

## 1.63:
* In this version, the shell matrices were started to be implemented, but this is incomplete. Don't use this functionality.
* The solution grid does have a ghost region now. The older versions crash for too small solution grids.
* New file 'SLEPC_vars.f90' that defines contexts for shell matrices as well as explicit interfaces.
* New routine 'get_ghost_vec' in SLEPC_utilities that allows one to get the ghost regions to the left and right.
* The solution variables are now stored in trimmed grid.

## 1.64:
* Removed shell functionality again. In the future it might be restored and finalized.
* Removed SLEPC_vars.f90 and its routines, as well as 'get_ghost_vec'.
* Added new input variable 'sol_n_procs' to control how many MPI processes are used for SLEPC. Often, 1 is better han multiple. A negative number sets it equal to all available.
* Started using block matrices in SLEPC, to improve legibility. No marked change in performance.
* By default, the generalized Davidson method is used. If there is convergence to a positive (stable) mode, maybe you should choose ncv higher.
* The run script has been changed somewhat, using only ncv=16 now.
* Checking whether the resulting eigenpairs are valid is now also done for the release version.
* Fixed a bug where the vacuum contribution was not copied to the integrated perturbation quantities.
* Changed the 'calc_norm_range' procedures somewhat, to take into account that the solution grid might have a different number of processes.
* Added hard-coded option to use Hermiticity, but this does not work yet, as the SLEPC tolerance is much too low (see http://lists.mcs.anl.gov/pipermail/petsc-users/2016-October/030781.html). Forcing the matrix to be Hermitian artificially also does not work.

## 1.65:
* Fixed a bug in the plotting of the magnetic grid.
* Even if no solution found, POST can still output equilibrium quantities.
* MUMPS warning has been removed, as it is not correct. The true source of slow SLEPC convergence was not using an optimzed version.
* Now a target value is used through 'EV_guess' in combination with shift-invert. By default, it is -1E-1. This speeds up convergence greatly.
* Introduced a new input variable 'solver_SLEPC_style' that can choose between a default setting for Krylov-Schur (1) or GD (2).
* Some more user messages concerning '--jump_to_sol' as well as while setting up the SLEPC solver.
* Changed message formatting somewhat.
* B_aligned was not intialized in POST; it has been removed.
* Extended grids can have monotomously decreasing coordinates now.
* 'calc_vec_comp' now takes into account possible bad points at r=0.
* 'calc_vec_comp' can now calculate fluxes.

## 1.66:
* The numerical derivatives in the normal direction are now performed before saving the HELENA and VMEC variables, in the full input grid.
* The procedure 'prepare_RZL' has been deleted and its functionality is now done in 'read_VMEC'.
* The procedure 'calc_eq_1' is now much lighter, and consists mainly of copying.
* Fixed a bug in 'B_plot' and 'J_plot' where 'plot_fluxes' was used wrongly while optional.
* Fixed some small bugs.
* The module VMEC has been split in three modules, the usual VMEC_ops, VMEC_utilities and VMEC_vars.

## 1.67:
* Fixed some output bugs.
* Bokeh now outputs checkboxes.
* The VMEC files modules were wrongly absent in the previous version.

## 1.68:
* Fixed bug when solution was reconstructed even when it did not exist.
* The 'calc_vec_comp' procedure now also calculates the 'Magnetic' components, which are in the (psi,theta,zeta) direction. This streamlines the integrated flux calculation as well.
* Fixed some bugs where normalization factors were not properly taken into account when calculating fluxes.
* 'extend_grid_E' is now called 'extend_grid_F' and works on the Flux variables. This is very important for when integration happens in these grids, so that the non-varying coordinate should be constant.
* Fixed some small bug in the division of equilibrium jobs for POST, where no maximum was provided.
* There is still an issue with the pressure balance. This is due to an inaccuracy in the calculation of the current from the magnetic field. This probably also causes a deviation from the correct results for the integrated current fluxes.

## 1.69:
* Extended 'debug_calc_derived_q' somewhat.
* Splines have been implemented for derivatives. Only the third order (cubic) has been implemented. This is now fixed.
* 'conv_FHM' has been removed, replaced by splines.
* Fixed issue in HDF5 plotting, where the symmetry type was not recoverable for vector plots.

## 1.70:
* Fixed an important bug in the 'fourier2real' where the incorrect sign was taken for the derivatives due to an incorrect application of integer division.
* The pressure balance is more correct, as is the test on D3sigma.
* Changed the tests in 'calc_derived_q' somewhat.
* Added a new test on the calculation of R, Z and L for VMEC, more specifically for the correct treatment of the derivatives.
* 'plot_HDF5' now also takes single X, Y and Z for a collection.

## 1.71:
* The discrepancy between the fluxes is finally fixed, probably. Might have to test with bigger 3-D effects.
* Fixed a bug in the calculation of fluxes. 
* Using '--plot_VMEC_modes' the decay of the VMEC Fourier modes can be investigated.

## 1.72:
* Adapted to XPS 9360.
* Adapted markdown file.
* New script to setup Run.

## 1.73: 
* Removed stellinstall, as it now has its own repository.
* Fixed small bugs considering release version.
* Fixed small bugs for new Petsc version.
* Fixed small bugs for mayavi. On Ubuntu, it now works *without* usin pip at all, and *only* with packages.

## 1.74:
* Fixed a bug in the export of HELENA, where modes with nonzero m were counted double.
* As 'prop_B_tor' has to be multiplied by the position ripple to yield the magnetic field ripple, its inverse needs to be used. This has been fixed.
* Currently, the implementation expects implicitly that the toroidal field ripple is constant for all poloidal positions. This is going to be changed.

## 1.75: 
* The translation between magnetic and position perturbation is not correct.
* 'nufft' now correctly takes into account functions that are defined in any interval.
* 'export_HEL' now also works for non top-bottom symmetric configurations.

## 1.76:
* New option 'compare_tor_pos' to compare quantities B, J and kappa at different toroidal positions.
* 'compare_tor_pos' also gives a resulting factor 'prop_B_tor' that has to be multiplied by the B ripple to get the position ripple.
* The old 'prop_B_tor' that was used in 'write_flux_q_in_file_for_VMEC' is not working correctly, as shown by direct comparison.
* This will be changed in the near future.
* 'write_flux_q_in_file_for_VMEC' now has some more output.
* Changed the default EV guess.
* In 'write_flux_q_in_file_for_VMEC' the cases with full normal output are now interpolated.

## 1.77:
* UNFINISHED VERSION: The 'write_flux_q_in_file_for_VMEC' is not functioning.
* Bug fixes for POST.
* If 'compare_tor_pos' is not called properly, an error now results.
* 'compare_tor_pos' now outputs a ripple map file.
* 'write_flux_q_in_file_for_VMEC' is now renamed to 'create_VMEC_input' and has been overhauled.
* 'create_VMEC_input' can now also read perturbations defined on R, Z and phi, corresponding to a custom origin R, Z. This is the output of DESCUR, for example.

## 1.78:
* The ability to use a proportionality file to translate a B_tor perturbation to a position perturbation is now implemented but needs testing.
* This file can be created using 'compare_tor_pos'.
* The splines routines have been replaced by the 'bspline_module' package by Jacob Williams (https://github.com/jacobwilliams/bspline-fortran.git), which is multi-dimensional.
* Tweaked external output a bit.

## 1.79:
* Bug fixes.
* Further tweaks to the system of perturbations of HELENA equilibria.
* Perturbation map can be shifted vertically, which is important for eqdisk files processed by HELENA, which shift the geometric axis to Z=0.
* There is now a check for physical consistency for the sign of the Alfven time.
* 'pause_prog' now takes a secret message 'stop' that can be used to stop the program, which is useful when controlling PB3D in an automated way.
* Improved the guess for theta_E as a function of theta_F in the procedure 'coord_F2E' by using lambda. This greatly improves convergence sometimes.

## 1.80:
* Bug fixes.
* The proportionality factor between magnetic and position ripples now use the geometric angle.
* 'compare_tor_pos' now needs 3 toroidal points, one in the middle. This makes toroidal comparisons easier.

## 1.81:
* Fixed a bug in the interpolation of non-symmetric prop_B_tor.
* Introduced 'order_per_fun' to order periodic functions and possible add an extra overlap left and right.

## 1.82:
* There is a large error with the calculation of the modes for 'create_VMEC_input'. It can most easily be seen when plotting the cross-section inverted.
* Fixed some bugs and added some checks in the calculation of the proportionality factors prop_B_tor.
* 'delta_r_plot' now returns half the differences of r and delta_B/B, which corresponds to the ripple definition in the litterature.
* Improved the whole 'create_VMEC_input' procedure, fixed the plotting.

## 1.83:
* Rewrote and greatly simplified 'create_VMEC_input'. It now uses the geometrical poloidal angle, as this angle needs to be kept constant for any reasonable comparison.
* Did the same thing for the 'delta_r_plot'.
* Changed 'calc_vec_comp' so that it also uses the geometric angle and radius.
* Splines is now merged into the more logical num_utilities module.
* The origin of the geometrical poloidal angle now has to be provided on the command line for compare_tor_pos.
* For compare_tor_pos, the fundamental theta interval is now needed always. Otherwise the interpolation in geometrical poloidal angle becomes quite messy.

## 1.84:
* Bug fixes for unperturbed HELENA equilibrium export.
* Removed the object oriented bspline module as it is unused and more difficult to compile in certain clusters.
* Bug fixes in SLEPC_ops.
* The maximum perturbation on axis is now specified absolutely for position perturbation.
* Fixed a bug in the creation of VMEC input files where modes that resulted in the same (m,n) were not added up but overwritten.

## 1.85: 
* There seems to be an error in energy reconstruction.
* Time is not stopped and started any more in output_ops.
* POST driver was doing some work twice. This has been removed.
* Slightly changed 'calc_E' to avoid segmentation faults in intel compilers caused by large temporary arrays.
* Fixed a bug for 'compare_tor_pos' with multiple processes: RZ_0 was not broadcasted.
* If no solution is found for 'compare_tor_pos', the difference is set to 0 and a warning is displayed'

## 1.86:
* Fixed a bug where no solution variables were found if the first level failed.
* Tweaked the calculation of the Jacobian in HELENA coordinates, as well as h_H_12.
* Now the Jacobian in VMEC coordinates is taken directly from VMEC, as at psi=0 there is a singularity using the formulas. The testing routine is now different as well.
* The plotting of fluxes now no longer fails for entire normal range.
* Implemented a routine 'solve_vand' that solves a Vandermonde matrix system.
* Using 'solve_vand' now in 'calc_deriv_data', instead of doing it manually with LAPACK.

## 1.87: 
* Inverse iteration for the solution of the generalized eigenvalue problem is not done any more.
* The minimization of surface energy now uses new theory, but it is still not correct.
* There is a new possibility for the boundary condition to be used, by employing asymmetric finite differences, but it does not work properly either.
* Looks like the best option is to use BC_style 3 that just neglects exterior points in the finite differences.
* There is a new possibility for the boundary condition to be used, by employing asymmetric finite differences, but it does not work properly either.
* Looks like the best option is to use BC_style 3 that just neglects exterior points in the finite differences.
* 'calc_coeff_fin_diff' now uses the Vandermonde method, and it can optionally produce asymmetric formula's.
* Fixed the bug where the vacuum contribution was not saved and therefore could not be used when jumping to the solution.

## 1.88:
* 1 to set to zero.
* 2 to use asymmetric finite differences close to the edge and delta_vac on the edge.
* 3 to extend the normal grid to accomodate finite differences on the edge, and delta_vac on the edge.
* 4 to explicitely impose the boundary condition on the edge.
* Intoruced 'norm_disc_style_sol' to enable for left finite differences as well.
* Boundary style 4 does not work well yet, which might be due to the absence of vacuum.
* The best method is currently the default: left finite differences with boundary condition style 2 (3 is identical for left differences).

## 1.89:
* NON-USABLE VERSION: IMPLEMENTING VACUUM.
* eq_ops does not deallocate R_E and Z_E anymore, as they are needed for vacuum.
* New module 'dtorh' that calculates toroidal functions.
* Strumpack is going to be used to solve structured matrices resulting from the vacuum.
* 'vac' is now split into 'vac_ops' and 'vac_vars'.

## 1.90:
* Slightly extended 'export_HEL' to also give output for current coil matching with Javier Artola's routines in Jorek.

## 1.91:
* First steps taken to use Doxygen.
* Bug fixes in 'export_HEL'.
* Fixed bug in splines with low order.
* coord_F2E_rtz now has an optional option to choose the order.
* Tuned the parameters in 'coord_F2E_rtz', because the current ones were too aggressive for pointy plasmas.
* Fixed bug in 'store_vac' where 'calc_vec_comp' was used incorrectly and very redundantly.
* 'calc_zero_HH' now shows output.
* eq_ops does not deallocate T_VC and jac_E anymore, as they are needed for vacuum.

## 1.92: 
* Changed to Doxygen formatting for comments.
* An accompanying document is still missing.

## 1.93:
* Automated the makfile, including 'finalize_version' to get a commit ready, which calls 'tag' and 'doc' as well as 'clean', 'PB3D' and 'POST'.
* 'make doc' creates the documentation, with the current version number.
* 'make tag' creates the git tag, with the current version number.
* Added the tags to the respective commits.
* Created documentation for installation and for inputs.

## 1.94:
* The Household methods have been modified by introducing backtracking instead of multiple tries.
* The poloidal flux is used as default, as the toroidal flux is experimental and untested.
* Fixed up some of the formatting in Doxygen for interfaces.
* References to interfaces and custom types are now correctly displayed, prepended by the module name.
* Added documentation for outputs and  general code structure and updated the main page.
* Changed the Doxygen documentation style and looks.

## 1.95:
* Fixed some formatting for code fragments.
* Changed 'examples' to 'tutorial' and moved the automated page of examples to files.
* Changed colorscheme and other design tweaks.
* Removed latex mistakes so that it compiles again.
* Polished up the latex files so that we now have a manual.
* Introduced 'clean_html_and_latex.sh' that does this.
* 'create_VMEC_input' is now a public routine, which makes more sense.

## 1.96:
* The tutorial is not yet ready yet, as there were numerical problems.
* Fixed some errors in manual pdf generation.
* Fixed small bugs in how many solutions are plot: By default all, and if none are chosen, it does not crash any more.
* Rewrote procedure 'set_nonzeros', fixing many bugs.
* Fixed bugs in setting options of SLEPC.
* Fixed bug in 'calc_zero_HH_3D' where zero_guess was not initialized at some points because the correction was zero.
* Fised bug in Bokeh drawing where wrong slash was used and the comma was forgotten.

## 1.97:
* Fixed doxygen bugs and added files in examples.

## 1.98:
* '.javifile' has been renamed to '.jorek'.
* HELENA now should also return B0 in the mapping file.
* Cleaned up, fixed, and improved the creation of VMEC file in 'create_VMEC_input'.
* 'create_VMEC_input' now can do either ncurr = 0 or ncurr = 1.
* Fixed bug in the JOREK export. There is a good match now.

## 1.99:
* New module 'vac_utilities' to house the utilities that are split off from 'vac_ops'.
* 'vac_vars' now include the angles along the magnetic field lines.
* Implemented 'calc_GH_2'.
* Implemented test on G and H, making use of test potential (R e^zeta)^n. Results are positive.
* Implemented second test, making use of a spherical potential of dubious validity, as also indicated by negative results.

## 2.00:
* Fixed some bugs in vacuum.
* For now, can only use poloidal flux for HELENA, because otherwise G and H of the vacuum would have to be calculated for every n, which is to be avoided for efficiency.
* Implemented vacuum for HELENA, but untested.
* It looks like H is almost equal to unity matrix times -2pi. Is this physical?

## 2.01:
* Implemented vacuum for axisymmetry, but not sure about the verification.
* There are rather large inconsistencies concerning the energy reconstruction, which have to be looked in to.
* Slightly changed the routine that calculates the optimal fast modes.
* Fixed an important bug in 'calc_GH_2' that introduced a strong asymmetry.
* 'calc_GH_int_2' now uses factors 1/2 instead of 1/3 and 1/6 to make the matrix G more symmetric.
* Fixed some bugs in the reading of vacuum variables.
* Added test on 'calc_vac' that looks at whether the vacuum response is positive definite through the eigenvalues.
* Improved output on the EV on the midplane, taking correct mode number now and using less HDF5 files.
* 'calc_res_surf' now outputs the real mode number, not the index.

## 2.02:
* Fixed an important bug in energy reconstruction where the pressure gradient was wrongly interpolated.
* The resonance plot now also outputs at least the safety factor or rotational transform when there is no resonance found.
* Increased the maximum size for an external plot.
* Fixed an bug in 'Weyl_fac', where a factor pi was left out.
* Fixed an inconsistency in 'calc_int_vol' where 2pi was left out.
* The energy reconstruction of the normal terms as well as the vacuum terms are now more correct for HELENA.
* Relaxed the test on admissible modes, displaying just warnings now.

## 2.03:
* Fixed a bug where the plot size was not correctly set in the external output routines.
* Implemented 'vac_pot_plot', which can be called through 'plot_vac_pot' to plot the vacuum potential.
* This procedure only works for 2-D vacua and is still experimental and untested.

## 2.04:
* Reorganized some vacuum things.
* Fixed bugs concerning the vacuum.
* Improved makefile for ITER.
* vacuum is now skipped for VMEC, without causing an error.

## 2.05:
* Implemented multiple field lines.
* Implemented new input variable 'alpha_style' that indicates whether to use a single field line with many turns, or multiple field lines with a single turn.
* In the latter case, 'min/max_par_X' is set to be a period of 2[pi], and 'n_alpha' is used to indicate the number of field lines.
* 'alpha' is now an array and situated in grid_vars.
* Updated 'magn_grid_plot' to handle multiple field lines.
* For VMEC, the current profiles can now be checked using 'J_plot' in POST.

## 2.06:
* UNUSABLE VERSION FOR VMEC.
* The old H and G results are now copied to the new vacuum in a next Richardson extrapolation level.
* Changed the structure of 'store_vac_VMEC'.
* Implemented new procedure 'interlaced_vac_copy' that performs a copy of vacuum variables G and H from a previous Richardson level in an interlaced way.

## 2.07:
* This version should run without error but is not debugged properly yet.
* 'eq_ops' does not deallocate T_EF anymore, as it is needed for vacuum.
* 'ang' is only copied for vacuum when using HELENA. For VMEC, the grid needs to be equidistant.
* The plan to only calculate the nonzero elements of G and H in 3-D vacuums has been abandoned. It is not general enough.
* Fixed some errors in the 3-D vacuum case.
* Updated the 'interlaced_vac_copy' to include all variables necessary and tested it as well.
* It performs its task for G and H as well, but they are later overwritten, as stated in bullet 3.
* Added trap for NaN to debug version on laptop.
* Changed definition of relative error due to insistence of Daan.

## 2.08:
* Vacuum not yet working for VMEC and numerical problems for top-bottom asymmetric equilibria. It is not sure whether these occur for HELENA as well as VMEC.
* Vacuum potential plot is now an input variable.
* 'calc_vac_res' does no need the grid any more.
* Added the forgotten step size in alpha to the vacuum response in 3-D.
* 'min_alpha' and 'max_alpha' are now saved without the factor pi for consistency with the parallel angle.
* Added tolerance to 'read_HEL' to prevent infinity at axis.
* There is an uncertainty about what to do with 'n_alpha' when there is only 1 field line. This should be accurate for axisymmetric cases but it is probably impossible.
* Fixed small bug where proportionality factor for ripple was set to zero for negative values.
* Added possibility to plot toroidal dependency of ripple, but quite rudimentary.
* Added output of ripple map in delta_r.
* Improved 'run.sh' so that it also works fine with absolute paths and gives more reasonable job names.
* Temporarily set vacuum calculation to zero for VMEC equilibria.
* Slighly changed 'test_metrics_H'.
* 'debug_calc_derived_q' now plots curvatures and shear in equilibrium grid. For VMEC this might not make sense.
* 'setup_deriv_data' now also can take an optional flag that indicates whether the signal is periodic. The test has been adapted appropriately. This option is now used in HELENA.
* The equidistant version of 'setup_deriv_data' now just calls the regular version with a dummy x.

## 2.09:
* The pressure derivative calculation for HELENA now checks whether it deviates much from HELENA output on first and last point.
* HELENA equilibrium calculations now calculate g_E before h_E.
* The HELENA equilibrium profiles are now much more accurate, especially at the plasma edge.

## 2.10:
* Alpha variables are now not stored in HDF5 but broadcasted with MPI.
* Implemented 'test_harm_cont_H' to plot Harmonic content in R_H and Z_H.
* 'debug_run_driver_X_2' now plots the full integrated coefficients PV and KV.

## 2.11:
* THIS VERSION DOES NOT WORK. FOR PB3D, NO INSTABILITIES ARE FOUND AND FOR POST, MORE DEVELOPMENT IS NEEDED.
* Implemented flag 'X_grid_style' that lets the user choose how to determine the normal component of the perturbation grid.
* For X_grid_style 1, the equilibrium normal grid is used, which corresponds to what, at the time of writing, seems like the best solution to avoid numerical problems when choosing a finer solution grid.
* For X_grid_style 2, the solution normal grid is used, which corresponds to the situation before version 2.11.
* 'calc_norm_range' now does not have a version for X any more, as this is taken by either the eq or the sol version, depending on 'X_grid_style'.
* Grids, equilibria and perturbations now have a copy procedure for deep copy.
* Unified 'setup_grid_X' and 'setup_grid_sol', which have complementary functions, depending on 'X_grid_style'.
* 'debug_X_norm' has been removed from solution driver.
* Simplified the SLEPC solution routines so that they only use the solution grid.
* The solution grid now is not just purely a normal grid, but has 1 parallel point and a number of geodesic points equal to the number of field lines, n_alpha.
* 'setup_nm_X' now also sets up a variable 'r_X' that indicates at which normal coordinate the n_X, m_X and sec_X_ind variables are set up.
* These variables can be interpolated on the new grid in 'interpolate_nm_x', and later restored in 'restore_nm_X'.
* 'debug_run_driver_X_2' now multiplies (hard-coded) the Z-axis by 100 to have easier plots.
* Fixed bugs in 'debug_run_driver_X_2' where multiple X jobs gave a wrong result. It is now all done after all the X jobs are done.

## 2.12:
* THIS VERSION HAS PB3D WORKING BUT NOT YET POST.
* Fixed bug when using X_grid_style 2, where the normal interpolation in the solution driver was not done correctly because the last normal point was always chosen.
* Fixed bug when using X_grid_style 2, where the solution grid was not correctly initialized.
* 'n_X', 'm_X', 'sec_X_ind' and 'r_X' are now saved in a custom type 'modes_type'.
* They are set in 'setup_modes' (formerly 'setup_nm_X'), of which 'init_modes' has been split off.
* 'init_nm_X' now only sets up the minimum and maximum mode numbers, in the equilibrium grid.
* 'interpolate_nm_X' and 'restore_nm_X' have been removed as they are no longer needed.
* Solution driver now also needs equilibrium grid to set up the solution modes.
* 'calc_XUQ' still outputs XUQ in the solution grid, but now also possibly interpolates the perturbation quantities if X_grid_style is 1.
* For POST_style 1, 'setup_out_grids' does not need full grids, as it will extend them.
* plot_sol_vec now uses less memory in the interpolated metric coefficients by selecting only the necessary ones.

## 2.13:
* UNUSABLE VERSION: BUG IN INTERPOLATION OF INTEGRATED X_2 QUANTITIES FOR FAST VERSION.
* Only quantities with same mode number (combination) can be interpolated between, which is currently not done.
* Alpha variables are once again stored in HDF5, because they are needed in POST.
* 'alpha' is now allocated in init_rich.
* For Richardson extrapolation, a bug is fixed when solution variables were set up to be used as guess: the grid has to be trimmed.
* Output concerning which Richardson level used for POST is now more helpful.
* The solution grid now again just a normal grid with 0 angular points. A local hybrid grid was implemented in driver_sol to interpolate for X_grid_style 2.
* If energy reconstruction not requested, it is not plotted any more.
* Implemented new option through 'V_interp_style' that allows the user to switch between 1 (finite difference, new default) and 2 (spline, previous default).

## 2.14:
* UNUSABLE VERSION: BUG IN INTERPOLATION OF INTEGRATED X_2 QUANTITIES FOR FAST VERSION FIXED, BUT STILL UNABLE TO REPRODUCE PREVIOUS RESULTS.
* Changed way in which secondary modes are stored: The index is kept constant now for a certain mode.
* 'setup_modes' has been adapted and now calculates the variable 'sec' which is part of the class 'modes' and which indicates the normal limits and table index of every mode. It also works for non-monotomous safety factors with possibly mulpiple ranges of same total mode.
* Also, the interpolation that is uses is now hard coded of precision 1 (i.e. linear). This is necessary to ensure that the mode ranges are consistent between the perturbation and solution grid for the case of X_grid_style 1.
* Improved Bokeh external output plotting so that it is more easily legible and does not throw an error for more than 255 plots.
* Debug information for X_1 and X_2 drivers are now procedures in X_ops, so they can be called externally as well as is done in solution driver.
* The old debug information for X_2 is not available any more, as it has been superseded by the real X_2 debug information that is also valid for X_style fast.
* 'setup_interp_data' now accepts extrapolation, and this is used for solution driver.
* 'setup_interp_data' now also accepts precision 0, which means using the constant single value.

## 2.15:
* USUABLE VERSION.
* There is still a problem when the equilibrium grid is not fine enough to allow for real interpolation in X_grid style 1 (eq) with X_style 2 (fast).
* A solution would be to define X_grid style 3 where the X grid is intermediary between the eq and the sol grid, in order to guarantee that the secondary mode range does not vary too fast in the normal direction.
* Another possible solution would be to limit the actual number of modes of the solution to to a number lesser than n_mod_X and to throw away those modes that cannot be treated appropriately with interpolation.
* Fixed bug in 'insert_block_mat' where k was used as index instead of m.
* Rewrote 'print_debug_X_1' in the style of 'print_debug_X_2'.
* 'setup_modes' now also plots the limits and the normal extent in debug mode.
* Minimal normal extent for modes is now monitored in 'interp_V'.
* For X_grid_style 1, 'init_modes' now sets the limits for the modes for the last normal point equal to the ones for the previous point, so that interpolation works better.
* Added option to remove possible previously present arrs in 'print_HDF5_arrs', which is currently used only to print the solution grid and variables in case of X_grid_style 2.

## 2.16:
* PB3D WORKS WELL NOW, BUT NEEDS SOME DEBUGGING TO MAKE SURE. THE NEW X_grid_style 3 WORKS VERY FAST AND ACCURATELY.
* POST NEEDS SOME MORE PROFOUND DEBUGGING, AND ALSO calc_E DOES NOT WORK YET.
* Implemented new X_grid_style 3 (optimized) where the perturbation variable are tabulated in a modified copy of the equilibrium grid that has been enriched in places where the safety factor changes too fast. It only makes sense to use this with X_style 2 (fast).
* Implemented new variable 'max_jq_change' that can be used with X_grid_style 3 (optimized) and X_style 2 (fast) to determine the maximum allowable change for the safety factor (pol. flux) or rotational transform (tor. flux) per normal grid point..
* Fixed bug in 'redistribute_output_X' where complex variables were not correctly redistributed.
* Fixed bug in external Bokeh drawing where interactivity of plots was not correct.
* Changed 'calc_norm_range' to make it more flexible and added PB3D_X.

## 2.17:
* INTERPOLATION ROUTINES ARE BAD: CONTINUITY IS NOT GUARANTEED. SPLINES ARE NECESSARY.
* Fixed bug in 'extend_grid_F' where total r_E was not set.
* Fixed bug in 'setup_grid_sol' where the perturbation grid was used wrongly as the equilibrium grid.
* Fixed bug where in the solution grid, the Equilibrium coordinates were not set and thus not written to output and therefore not correctly used in the postprocessing.
* In 'calc_norm_range' for POST, the perturbation limits are set explicitely to the solution limits for X_grid_style 2 (perturbation).
* Fixed bug in 'calc_XUQ' where normal derivative was done on the whole normal grid, even for X_style 2 (fast). The appropriate subset is now taken.
* Fixed bug in initialization of X_grid_style.
