!------------------------------------------------------------------------------!
!> Numerical variables used by most other modules.
!------------------------------------------------------------------------------!
module num_vars
    use ISO_FORTRAN_ENV
    
    implicit none
    private
    public dp, dpi, max_str_ln, max_name_ln, max_deriv, prog_name, &
        &output_name, prog_version, prog_style, min_PB3D_version, &
        &debug_version, &
        &shell_commands_name, mem_usage_name, &
        &mem_usage_count, weight_dp, rank, n_procs, sol_n_procs, time_start, &
        &max_tot_mem, max_X_mem, X_jobs_lims, X_job_nr, &
        &eq_jobs_lims, eq_job_nr, mem_scale_fac, pi, mu_0_original, iu, &
        &EV_style, eq_style, rho_style, U_style, norm_style, BC_style, &
        &X_style, X_grid_style, matrix_SLEPC_style, plot_resonance, &
        &plot_magn_grid, plot_B, plot_J, plot_flux_q, plot_kappa, plot_sol_xi, &
        &plot_sol_Q, plot_E_rec, plot_vac_pot, ltest, use_pol_flux_E, &
        &use_pol_flux_F, use_normalization, EV_BC, tol_SLEPC, max_it_slepc, &
        &norm_disc_prec_eq, K_style, EV_guess, norm_disc_prec_X, &
        &norm_disc_prec_sol, norm_disc_style_sol, &
        &POST_style, alpha_style, magn_int_style, solver_SLEPC_style, &
        &max_njq_change, &
        &max_it_rich, tol_rich, max_it_zero, max_nr_backtracks_HH, &
        &tol_zero, tol_norm, &
        &ex_max_size, eq_name, &
        &no_plots, no_output, plot_dir, script_dir, data_dir, n_theta_plot, &
        &n_zeta_plot, min_theta_plot, max_theta_plot, min_zeta_plot, &
        &max_zeta_plot, min_r_plot, max_r_plot, n_sol_requested, &
        &n_sol_plotted, retain_all_sol, do_execute_command_line, &
        &print_mem_usage, input_name, plot_grid_style, swap_angles, &
        &rich_restart_lvl, plot_size, jump_to_sol, compare_tor_pos, &
        &export_HEL, plot_VMEC_modes, ex_plot_style, pert_mult_factor_POST, &
        &min_Rvac_plot, max_Rvac_plot, min_Zvac_plot, max_Zvac_plot, &
        &n_vac_plot, invert_top_bottom_H, &
        &RZ_0, &
        &POST_output_full, POST_output_sol, &
        &shell_commands_i, mem_usage_i, output_EV_i, decomp_i, &
        &HEL_pert_i, HEL_export_i, input_i, PB3D_i, PB3D_name, eq_i, output_i, &
        &prop_B_tor_i

    ! technical variables
    !integer, parameter :: dp = kind(1.d0)                                       !< double precision
    !integer, parameter :: qp = selected_real_kind (32)                          !< quadruple precision
    integer, parameter :: dp = REAL64                                           !< double precision
    !integer, parameter :: qp = REAL128                                          !< quadruple precision
    integer, parameter :: dpi = INT64                                           !< double precision
    real(dp), parameter :: weight_dp = 0.008                                    !< size of double precision in kB
    integer, parameter :: max_str_ln = 120                                      !< maximum length of strings
    integer, parameter :: max_name_ln = 30                                      !< maximum length of filenames
    integer, parameter :: max_deriv = 2                                         !< highest derivatives for metric factors in Flux coords.
    integer :: prog_style                                                       !< program style (1: PB3D, 2: PB3D_POST)
    character(len=4) :: prog_name                                               !< name of program, used for info
    character(len=3), parameter :: output_name = 'out'                          !< name of output file
    character(len=14), parameter :: shell_commands_name = 'shell_commands'      !< name of shell commands file
    character(len=9), parameter :: mem_usage_name = 'mem_usage'                 !< name of memory usage file
    integer :: mem_usage_count                                                  !< counter for memory usage output
    real(dp), parameter :: prog_version = 2.26_dp                               !< version number
    real(dp), parameter :: min_PB3D_version = 2.24_dp                           !< minimum PB3D version for POST
#if ldebug
    logical :: debug_version = .true.                                           !< debug version used
#else
    logical :: debug_version = .false.                                          !< not debug version
#endif

    ! MPI variables
    integer :: rank                                                             !< MPI rank
    integer :: n_procs                                                          !< nr. of MPI processes
    integer :: sol_n_procs                                                      !< nr. of MPI processes for solution with SLEPC
    integer(kind=8) :: time_start                                               !< start time of simulation
    
    ! job variables
    real(dp) :: max_tot_mem                                                     !< maximum total memory for all processes [MB]
    real(dp) :: max_X_mem                                                       !< maximum memory for perturbation calculations for all processes [MB]
    integer, allocatable :: X_jobs_lims(:,:)                                    !< data about X jobs: [\f$\min_k\f$, \f$\max_k\f$, \f$\min_m\f$, \f$\max_m\f$] for all jobs
    integer, allocatable :: eq_jobs_lims(:,:)                                   !< data about eq jobs: [\f$\min_\theta\f$, \f$\max_\theta\f$] for all jobs
    integer :: X_job_nr                                                         !< nr. of X job
    integer :: eq_job_nr                                                        !< nr. of eq job
    real(dp), parameter :: mem_scale_fac = 6.0                                  !< empirical scale factor of memory to calculate eq compared to just storing it 

    ! physical and mathematical variables
    real(dp), parameter :: pi=4_dp*datan(1.0_dp)                                !< \f$\pi\f$
    real(dp), parameter :: mu_0_original = 4E-7_dp*pi                           !< permeability of free space
    complex(dp), parameter :: iu = (0,1)                                        !< complex unit

    ! concerning runtime
    integer :: EV_style                                                         !< determines the method used for solving an EV problem
    integer :: eq_style                                                         !< either 1 (VMEC) or 2 (HELENA)
    integer :: rho_style                                                        !< style for equilibrium density profile
    integer :: U_style                                                          !< style for calculation of U (1: ord.2, 2: ord.1, 1: ord.0)
    integer :: norm_style                                                       !< style for normalization
    integer :: K_style                                                          !< style for kinetic energy
    integer :: BC_style(2)                                                      !< style for BC left and right
    integer :: X_style                                                          !< style for secondary mode numbers (1: prescribed, 2: fast)
    integer :: matrix_SLEPC_style                                               !< style for matrix storage (1: sparse, 2: shell)
    integer :: solver_SLEPC_style                                               !< style for solver (1: Krylov-Schur, 2: GD)
    integer :: POST_style                                                       !< style for POST (1: extended grid, 2: B-aligned grid)
    integer :: X_grid_style                                                     !< style for normal component of X grid (1: eq, 2: sol, 3: enriched)
    integer :: alpha_style                                                      !< style for alpha (1: one field line, many turns, 2: many field lines, one turn)
    integer :: max_it_slepc                                                     !< maximum nr. of iterations for SLEPC
    logical :: plot_resonance                                                   !< whether to plot the q-profile or iota-profile with resonances
    logical :: plot_magn_grid                                                   !< whether to plot the grid in real coordinates
    logical :: plot_B                                                           !< whether to plot the magnetic field in real coordinates
    logical :: plot_J                                                           !< whether to plot the current in real coordinates
    logical :: plot_flux_q                                                      !< whether to plot flux quantities in real coordinates
    logical :: plot_kappa                                                       !< whether to plot curvature
    logical :: plot_sol_xi                                                      !< whether to plot plasma perturbation of solution in POST
    logical :: plot_sol_Q                                                       !< whether to plot magnetic perturbation of solution in POST
    logical :: plot_vac_pot                                                     !< whether to plot vacuum potential in POST
    logical :: plot_E_rec                                                       !< whether to plot energy reconstruction in POST
    logical :: ltest                                                            !< whether or not to call the testing routines
    logical :: use_pol_flux_E                                                   !< whether poloidal flux is used in E coords.
    logical :: use_pol_flux_F                                                   !< whether poloidal flux is used in F coords.
    logical :: use_normalization                                                !< whether to use normalization or not
    real(dp) :: EV_BC                                                           !< value of artificial Eigenvalue for boundary condition
    real(dp) :: EV_guess                                                        !< first guess for eigenvalue
    real(dp), allocatable :: tol_SLEPC(:)                                       !< tolerance for SLEPC for different Richardson levels
    real(dp) :: max_njq_change                                                  !< maximum change of prim. mode number times saf.  fac. / rot. transf. when using X_style 2 (fast)
    integer :: norm_disc_prec_eq                                                !< precision for normal discretization for equilibrium
    integer :: norm_disc_prec_X                                                 !< precision for normal discretization for perturbation
    integer :: norm_disc_prec_sol                                               !< precision for normal discretization for solution
    integer :: norm_disc_style_sol                                              !< style for normal discretization for solution (1: central fin. diff., 2: left fin. diff.)
    integer :: magn_int_style                                                   !< style for magnetic integrals (1: trapezoidal, 2: Simpson 3/8)
    
    ! concerning Richardson extrapolation
    integer :: max_it_rich                                                      !< number of levels for Richardson extrapolation
    real(dp) :: tol_rich                                                        !< tolerance for Richardson extrapolation
    
    ! concerning finding the magnetic field lines
    integer :: max_it_zero                                                      !< maximum number of iterations to find zeros
    integer :: max_nr_backtracks_HH                                             !< maximum number of backtracks for Householder, relax. factors
    real(dp) :: tol_zero                                                        !< tolerance for zeros
    real(dp) :: tol_norm                                                        !< tolerance for normal range (normalized to 0..1)

    ! concerning input / output
    integer, parameter :: ex_max_size = 10000                                   !< maximum size of matrices for external plot
    character(len=max_str_ln) :: eq_name                                        !< name of equilibrium file from VMEC or HELENA
    character(len=max_str_ln) :: PB3D_name                                      !< name of PB3D output file
    logical :: no_plots = .false.                                               !< no plots made
    logical :: jump_to_sol = .false.                                            !< jump to solution
    logical :: export_HEL = .false.                                             !< export HELENA
    logical :: plot_VMEC_modes = .false.                                        !< plot VMEC modes
    logical :: invert_top_bottom_H = .false.                                    !< invert top and bottom for HELENA equilibria
    logical :: no_output = .false.                                              !< no output shown
    logical :: POST_output_full = .false.                                       !< POST has output on full grids
    logical :: POST_output_sol = .false.                                        !< POST has outputs of solution
    logical :: do_execute_command_line = .false.                                !< call "execute_command_line" inside program
    logical :: print_mem_usage = .false.                                        !< print memory usage is printed
    logical :: swap_angles = .false.                                            !< swap angles theta and zeta in plots (only for POST)
    logical :: compare_tor_pos = .false.                                        !< compare quantities at toroidal positions (only for POST)
    logical :: retain_all_sol                                                   !< retain also faulty solutions
    character(len=5) :: plot_dir = 'Plots'                                      !< directory where to save plots
    character(len=7) :: script_dir = 'Scripts'                                  !< directory where to save scripts for plots
    character(len=4) :: data_dir = 'Data'                                       !< directory where to save data for plots
    integer :: n_theta_plot                                                     !< nr. of poloidal points in plot
    integer :: n_zeta_plot                                                      !< nr. of toroidal points in plot
    integer :: n_vac_plot(2)                                                    !< nr. of points  in R and Z in vacuum
    real(dp) :: min_theta_plot                                                  !< min. of \c theta_plot
    real(dp) :: max_theta_plot                                                  !< max. of \c theta_plot
    real(dp) :: min_zeta_plot                                                   !< min. of \c zeta_plot
    real(dp) :: max_zeta_plot                                                   !< max. of \c zeta_plot
    real(dp) :: min_r_plot                                                      !< min. of \c r_plot
    real(dp) :: max_r_plot                                                      !< max. of \c r_plot
    real(dp) :: min_Rvac_plot                                                   !< min. of \c R in which to plot vacuum
    real(dp) :: max_Rvac_plot                                                   !< max. of \c R in which to plot vacuum
    real(dp) :: min_Zvac_plot                                                   !< min. of \c Z in which to plot vacuum
    real(dp) :: max_Zvac_plot                                                   !< max. of \c Z in which to plot vacuum
    integer :: plot_grid_style                                                  !< style for POST plot grid (0: 3-D plots, 1: slab plots, 2: slab plots with folding, 3: straight cylinder))
    integer :: n_sol_requested                                                  !< how many solutions requested
    integer :: n_sol_plotted(4)                                                 !< how many solutions to be plot (first unstable, last unstable, first stable, last stable)
    integer :: plot_size(2)                                                     !< size of plot in inches
    integer :: rich_restart_lvl                                                 !< starting Richardson level (0: none [default])
    character(len=max_str_ln) :: input_name                                     !< will hold the full name of the input file
    integer :: ex_plot_style                                                    !< external plot style (1: GNUPlot, 2: Bokeh for 2D, Mayavi for 3D)
    real(dp) :: pert_mult_factor_POST                                           !< factor with which to multiply perturbation strength for POST
    
    ! Concerning comparing toroidal positions
    real(dp) :: RZ_0(2)                                                         !< origin of geometrical poloidal coordinate
    
    ! Concerning file numbers
    integer, parameter :: input_i = 50                                          !< file number of input file
    integer, parameter :: eq_i = 51                                             !< file number of equilibrium file from VMEC or HELENA
    integer, parameter :: PB3D_i = 52                                           !< file number of PB3D output file
    integer, parameter :: output_i = 53                                         !< file number of output file
    integer, parameter :: shell_commands_i = 54                                 !< file number of shell commands file
    integer, parameter :: mem_usage_i = 55                                      !< file number of memory usage file
    integer, parameter :: output_EV_i = 56                                      !< file number of output of EV
    integer, parameter :: decomp_i = 57                                         !< file number of output of EV decomposition
    integer, parameter :: HEL_export_i = 59                                     !< file number of output of HELENA equilibrium export file
    integer, parameter :: HEL_pert_i = 58                                       !< file number of HELENA equilibrium perturbation file
    integer, parameter :: prop_B_tor_i = 60                                     !< file number of \f$B_\phi\f$ proportionality factor file
end module num_vars
