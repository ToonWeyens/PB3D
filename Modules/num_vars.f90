!------------------------------------------------------------------------------!
!   Numerical variables used by most other modules                             !
!------------------------------------------------------------------------------!
module num_vars
    use ISO_FORTRAN_ENV
    
    implicit none
    private
    public dp, dpi, max_str_ln, max_name_ln, max_deriv, prog_name, &
        &output_name, prog_version, prog_style, min_PB3D_version, &
        &shell_commands_name, mem_usage_name, mem_usage_i, mem_usage_count, &
        &weight_dp, &
        &rank, n_procs, time_start, &
        &max_tot_mem_per_proc, max_X_mem_per_proc, X_jobs_lims, X_jobs_taken, &
        &X_job_nr, X_jobs_file_name, eq_jobs_lims, eq_job_nr, mem_scale_fac, &
        &pi, mu_0_original, iu, &
        &EV_style, eq_style, rho_style, U_style, norm_style, BC_style, &
        &X_style, matrix_SLEPC_style, plot_resonance, plot_magn_grid, &
        &plot_flux_q, ltest, use_pol_flux_E, use_pol_flux_F, &
        &use_normalization, EV_BC, tol_SLEPC, max_it_slepc, &
        &norm_disc_prec_eq, K_style, norm_disc_prec_X, norm_disc_prec_sol, &
        &POST_style, magn_int_style, &
        &max_it_rich, tol_rich, &
        &max_it_inv, &
        &max_it_zero, max_nr_tries_HH, relax_fac_HH, tol_zero, tol_norm, &
        &def_relax_fac_HH, &
        &ex_max_size, input_i, PB3D_i, PB3D_name, eq_i, eq_name, output_i, &
        &no_plots, no_output, plot_dir, script_dir, data_dir, n_theta_plot, &
        &n_zeta_plot, min_theta_plot, max_theta_plot, min_zeta_plot, &
        &max_zeta_plot, n_sol_requested, n_sol_plotted, retain_all_sol, &
        &do_execute_command_line, print_mem_usage, input_name, slab_plots, &
        &swap_angles, rich_restart_lvl, plot_size, minim_output, PB3D_name_eq, &
        &jump_to_sol, export_HEL, ex_plot_style

    ! technical variables
    !integer, parameter :: dp = kind(1.d0)                                       ! double precision
    !integer, parameter :: qp = selected_real_kind (32)                          ! quadruple precision
    integer, parameter :: dp = REAL64                                           ! double precision
    !integer, parameter :: qp = REAL128                                          ! quadruple precision
    integer, parameter :: dpi = INT64                                           ! double precision
    real(dp), parameter :: weight_dp = 0.008                                    ! size of double precision in kB
    integer, parameter :: max_str_ln = 120                                      ! maximum length of strings
    integer, parameter :: max_name_ln = 30                                      ! maximum length of filenames
    integer, parameter :: max_deriv = 2                                         ! highest derivatives for metric factors in Flux coords.
    integer :: prog_style                                                       ! program style (1: PB3D, 2: PB3D_POST)
    character(len=4) :: prog_name                                               ! name of program, used for info
    character(len=3), parameter :: output_name = 'out'                          ! name of output file
    character(len=14), parameter :: shell_commands_name = 'shell_commands'      ! name of shell commands file
    character(len=9), parameter :: mem_usage_name = 'mem_usage'                 ! name of memory usage file
    integer :: mem_usage_count                                                  ! counter for memory usage output
    integer, parameter :: mem_usage_i = 100                                     ! has to be fixed, so should be chosen high enough
    real(dp), parameter :: prog_version = 1.35_dp                               ! version number
    real(dp), parameter :: min_PB3D_version = 1.32_dp                           ! minimum PB3D version for POST

    ! MPI variables
    integer :: rank                                                             ! MPI rank
    integer :: n_procs                                                          ! nr. of MPI processes
    integer(kind=8) :: time_start                                               ! start time of simulation
    
    ! job variables
    real(dp) :: max_tot_mem_per_proc                                            ! maximum total memory per process [MB]
    real(dp) :: max_X_mem_per_proc                                              ! maximum memory for perturbation calculations per process [MB]
    integer, allocatable :: X_jobs_lims(:,:)                                    ! data about X jobs: [min_k, max_k, min_m, max_m] for all jobs
    logical, allocatable :: X_jobs_taken(:)                                     ! X jobs taken
    integer, allocatable :: eq_jobs_lims(:,:)                                   ! data about eq jobs: [min_theta,max_theta] for all jobs
    integer :: X_job_nr                                                         ! nr. of X job
    integer :: eq_job_nr                                                        ! nr. of eq job
    character(len=10) :: X_jobs_file_name = 'X_jobs.txt'                        ! name of X jobs file
    real(dp), parameter :: mem_scale_fac = 1.5                                  ! empirical scale factor of memory (because operations are done)

    ! physical and mathematical variables
    real(dp), parameter :: pi=4_dp*datan(1.0_dp)                                ! pi
    real(dp), parameter :: mu_0_original = 4E-7_dp*pi                           ! permeability of free space
    complex(dp), parameter :: iu = (0,1)                                        ! complex unit

    ! concerning runtime
    integer :: EV_style                                                         ! determines the method used for solving an EV problem
    integer :: eq_style                                                         ! either 1 (VMEC) or 2 (HELENA)
    integer :: rho_style                                                        ! style for equilibrium density profile
    integer :: U_style                                                          ! style for calculation of U (1: ord.2, 2: ord.1, 1: ord.0)
    integer :: norm_style                                                       ! style for normalization
    integer :: K_style                                                          ! style for kinetic energy
    integer :: BC_style(2)                                                      ! style for BC left and right
    integer :: X_style                                                          ! style for secondary mode numbers (1: prescribed, 2: fast)
    integer :: matrix_SLEPC_style                                               ! style for matrix storage (1: sparse, 2: shell)
    integer :: POST_style                                                       ! style for POST (1: extended grid, 2: B-aligned grid)
    integer :: max_it_slepc                                                     ! maximum nr. of iterations for SLEPC
    logical :: plot_resonance                                                   ! whether to plot the q-profile or iota-profile with resonances
    logical :: plot_magn_grid                                                   ! whether to plot the grid in real coordinates
    logical :: plot_flux_q                                                      ! whether to plot flux quantities in real coordinates
    logical :: ltest                                                            ! whether or not to call the testing routines
    logical :: use_pol_flux_E                                                   ! whether poloidal flux is used in E coords.
    logical :: use_pol_flux_F                                                   ! whether poloidal flux is used in F coords.
    logical :: use_normalization                                                ! whether to use normalization or not
    real(dp) :: EV_BC                                                           ! value of artificial Eigenvalue for boundary condition
    real(dp), allocatable :: tol_SLEPC(:)                                       ! tolerance for SLEPC for different Richardson levels
    integer :: norm_disc_prec_eq                                                ! precision for normal discretization for equilibrium
    integer :: norm_disc_prec_X                                                 ! precision for normal discretization for perturbation
    integer :: norm_disc_prec_sol                                               ! precision for normal discretization for solution
    integer :: magn_int_style                                                   ! style for magnetic integrals (1: trapezoidal, 2: Simpson 3/8)
    
    ! concerning Richardson extrapolation
    integer :: max_it_rich                                                      ! number of levels for Richardson extrapolation
    real(dp) :: tol_rich                                                        ! tolerance for Richardson extrapolation
    
    ! concerning finding the inverse
    integer :: max_it_inv                                                       ! maximum number of iterations to find the inverse

    ! concerning finding the magnetic field lines
    integer :: max_it_zero                                                      ! maximum number of iterations to find zeros
    integer :: max_nr_tries_HH                                                  ! maximum number of tries for Householder, relax. factors
    real(dp), parameter :: def_relax_fac_HH = 0.5                               ! default relax_fac_HH (can be chosen higher for higher orders)
    real(dp) :: relax_fac_HH                                                    ! standard relaxation factor for Householder iterations
    real(dp) :: tol_zero                                                        ! tolerance for zeros
    real(dp) :: tol_norm                                                        ! tolerance for normal range (normalized to 0..1)

    ! concerning input / output
    integer, parameter :: ex_max_size = 300                                     ! maximum size of matrices for external plot
    integer :: input_i                                                          ! file number of input file
    integer :: eq_i                                                             ! file number of equilibrium file from VMEC or HELENA
    character(len=max_str_ln) :: eq_name                                        ! name of equilibrium file from VMEC or HELENA
    integer :: PB3D_i                                                           ! file number of PB3D output file
    character(len=max_str_ln) :: PB3D_name                                      ! name of PB3D output file
    character(len=max_str_ln) :: PB3D_name_eq                                   ! name of PB3D output file for vars on eq grid (see minim_output)
    integer :: output_i                                                         ! file number of output file
    logical :: no_plots = .false.                                               ! no plots made
    logical :: jump_to_sol = .false.                                            ! jump to solution
    logical :: export_HEL = .false.                                             ! export HELENA
    logical :: no_output = .false.                                              ! no output shown
    logical :: do_execute_command_line = .false.                                ! call "execute_command_line" inside program
    logical :: print_mem_usage = .false.                                        ! print memory usage is printed
    logical :: swap_angles = .false.                                            ! swap angles theta and zeta in plots (only for POST)
    logical :: minim_output = .false.                                           ! minimize output file size
    logical :: retain_all_sol                                                   ! retain also faulty solutions
    logical :: slab_plots                                                       ! slab plots (only for POST)
    character(len=5) :: plot_dir = 'Plots'                                      ! directory where to save plots
    character(len=7) :: script_dir = 'Scripts'                                  ! directory where to save scripts for plots
    character(len=4) :: data_dir = 'Data'                                       ! directory where to save data for plots
    integer :: n_theta_plot                                                     ! nr. of poloidal points in plot
    integer :: n_zeta_plot                                                      ! nr. of toroidal points in plot
    real(dp) :: min_theta_plot, max_theta_plot                                  ! min. and max. of theta_plot
    real(dp) :: min_zeta_plot, max_zeta_plot                                    ! min. and max. of zeta_plot
    integer :: n_sol_requested                                                  ! how many solutions requested
    integer :: n_sol_plotted(4)                                                 ! how many solutions to be plot (first unstable, last unstable, first stable, last stable)
    integer :: plot_size(2)                                                     ! size of plot in inches
    integer :: rich_restart_lvl                                                 ! starting Richardson level (0: none [default])
    character(len=max_str_ln) :: input_name                                     ! will hold the full name of the input file
    integer :: ex_plot_style                                                    ! external plot style (1: GNUPlot, 2: Bokeh for 2D, Mayavi for 3D)
end module num_vars
