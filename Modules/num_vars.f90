!------------------------------------------------------------------------------!
!   Numerical variables used by most other modules                             !
!------------------------------------------------------------------------------!
module num_vars
    use ISO_FORTRAN_ENV
    
    implicit none
    private
    public dp, qp, max_str_ln, max_name_ln, max_deriv, prog_name, output_name, &
        &prog_version, prog_style, &
        &max_mem_per_proc, n_procs, rank, X_jobs_lims, X_jobs_taken, X_job_nr, &
        &X_jobs_file_name, X_jobs_lock_file_name, HDF5_lock_file_name, &
        &pi, mu_0_original, iu, &
        &EV_style, eq_style, rho_style, BC_style, plot_resonance, plot_grid, &
        &plot_flux_q, ltest, use_pol_flux_E, use_pol_flux_F, &
        &use_normalization, EV_BC, tol_SLEPC, max_it_slepc, rich_lvl_nr, &
        &norm_disc_prec_eq, norm_disc_prec_X, norm_disc_prec_sol, &
        &max_it_r, tol_r, no_guess, &
        &max_it_inv, &
        &max_it_NR, tol_NR, nyq_fac, tol_norm_r, &
        &GP_max_size, input_i, PB3D_i, PB3D_name, eq_i, eq_name, output_i, &
        &no_plots, no_messages, plot_dir, script_dir, data_dir, n_theta_plot, &
        &n_zeta_plot, n_sol_requested, n_sol_plotted, retain_all_sol, &
        &spline_type

    ! technical variables
    !integer, parameter :: dp = kind(1.d0)                                       ! double precision
    !integer, parameter :: qp = selected_real_kind (32)                          ! quadruple precision
    integer, parameter :: dp = REAL64                                           ! double precision
    integer, parameter :: qp = REAL128                                          ! quadruple precision
    integer, parameter :: max_str_ln = 120                                      ! maximum length of strings
    integer, parameter :: max_name_ln = 30                                      ! maximum length of filenames
    integer, parameter :: max_deriv = 2                                         ! highest derivatives calculated for metric factors
    integer :: prog_style                                                       ! program style (1: PB3D, 2: PB3D_POST)
    character(len=4) :: prog_name                                               ! name of program, used for info
    character(len=3), parameter :: output_name = 'out'                          ! name of output file
    real(dp), parameter :: prog_version = 0.99_dp                               ! version number

    ! MPI variables
    real(dp) :: max_mem_per_proc                                                ! maximum memory per process [MB]
    integer :: rank                                                             ! MPI rank
    integer :: n_procs                                                          ! nr. of MPI processes
    integer, allocatable :: X_jobs_lims(:,:)                                    ! data about X jobs: [min_k, max_k, min_m, max_m] for all jobs
    logical, allocatable :: X_jobs_taken(:)                                     ! X jobs taken
    integer :: X_job_nr                                                         ! nr. of X job
    character(len=10) :: X_jobs_file_name = 'X_jobs.txt'                        ! name of X jobs file
    character(len=12) :: X_jobs_lock_file_name = '.lock_file_X'                 ! name of lock file for X jobs
    character(len=15) :: HDF5_lock_file_name = '.lock_file_HDF5'                ! name of lock file for HDF5 operations

    ! physical and mathematical variables
    real(dp), parameter :: pi=4_dp*datan(1.0_dp)                                ! pi
    real(dp), parameter :: mu_0_original = 4E-7_dp*pi                           ! permeability of free space
    complex(dp), parameter :: iu = (0,1)                                        ! complex unit

    ! concerning runtime
    integer :: EV_style                                                         ! determines the method used for solving an EV problem
    integer :: eq_style                                                         ! either 1 (VMEC) or 2 (HELENA)
    integer :: rho_style                                                        ! style for equilibrium density profile
    integer :: BC_style(2)                                                      ! style for BC left and right
    integer :: max_it_slepc                                                     ! maximum nr. of iterations for SLEPC
    logical :: plot_resonance                                                   ! whether to plot the q-profile or iota-profile with resonances
    logical :: plot_grid                                                        ! whether to plot the grid in real coordinates
    logical :: plot_flux_q                                                      ! whether to plot flux quantities in real coordinates
    logical :: ltest                                                            ! whether or not to call the testing routines
    logical :: use_pol_flux_E                                                   ! whether poloidal flux is used in E coords.
    logical :: use_pol_flux_F                                                   ! whether poloidal flux is used in F coords.
    logical :: use_normalization                                                ! whether to use normalization or not
    real(dp) :: EV_BC
    real(dp) :: tol_SLEPC                                                       ! tolerance for SLEPC
    integer :: rich_lvl_nr                                                      ! which Richardson's level is being calculated
    integer :: norm_disc_prec_eq                                                ! precision for normal discretization for equilibrium
    integer :: norm_disc_prec_X                                                 ! precision for normal discretization for perturbation
    integer :: norm_disc_prec_sol                                               ! precision for normal discretization for solution
    
    ! concerning Richardson extrapolation
    integer :: max_it_r                                                         ! number of levels for Richardson extrapolation
    real(dp) :: tol_r                                                           ! tolerance for Richardson extrapolation
    logical :: no_guess = .false.                                               ! disable guessing Eigenfunction from previous level of Richardson
    
    ! concerning finding the inverse
    integer :: max_it_inv                                                       ! maximum number of iterations to find the inverse

    ! concerning finding the magnetic field lines
    integer :: max_it_NR                                                        ! maximum number of Newton-Rhapson iterations
    real(dp) :: tol_NR                                                          ! tolerance for Newton-Rhapson
    integer :: nyq_fac                                                          ! Nyquist factor to avoid aliasing in perturbation integrals
    real(dp) :: tol_norm_r                                                      ! tolerance for normal range (normalized to 0..1)

    ! concerning input / output
    integer, parameter :: GP_max_size = 300                                     ! maximum size of matrices for GNUPlot
    integer :: input_i                                                          ! file number of input file
    integer :: eq_i                                                             ! file number of equilibrium file from VMEC or HELENA
    character(len=max_str_ln) :: eq_name                                        ! name of equilibrium file from VMEC or HELENA
    integer :: PB3D_i                                                           ! file number of PB3D output file
    character(len=max_str_ln) :: PB3D_name                                      ! name of PB3D output file
    integer :: output_i                                                         ! file number of output file
    logical :: no_plots = .false.                                               ! true if no plots should be made
    logical :: no_messages = .false.                                            ! true if no messages should be shown
    character(len=5) :: plot_dir = 'Plots'                                      ! directory where to save plots
    character(len=7) :: script_dir = 'Scripts'                                  ! directory where to save scripts for plots
    character(len=4) :: data_dir = 'Data'                                       ! directory where to save data for plots
    integer :: n_theta_plot                                                     ! nr. of poloidal points in plot
    integer :: n_zeta_plot                                                      ! nr. of toroidal points in plot
    integer :: n_sol_requested                                                  ! how many solutions requested
    integer :: n_sol_plotted(4)                                                 ! how many solutions to be plot (first unstable, last unstable, first stable, last stable)
    logical :: retain_all_sol                                                   ! retain also faulty solutions
    
    ! concerning  spline interpolation
    ! The type of the spline is determined by "spline_type":
    !   1:  cubic
    ! The coefficients that are stored here are:
    !   x:  the abscissa at which the spline was calculated
    !   y:  the ordinates at which the spline was calculated
    !   z:  the calculated second order derivatives at the points x
    type :: spline_type
        integer :: spline_type                                                  ! determines type of spline
        real(dp), allocatable :: x(:)                                           ! abscissa
        real(dp), allocatable :: y(:)                                           ! ordinate
        real(dp), allocatable :: z(:)                                           ! second order derivatives of spline
    end type spline_type
end module num_vars
