!------------------------------------------------------------------------------!
!   numerical variables used by most other modules                             !
!------------------------------------------------------------------------------!
module num_vars
    implicit none
    private
    public dp, qp, minim_style, max_str_ln, n_seq_0, max_args, max_opts, &
        &prog_name, max_it_r, tol_r, ltest, pi, max_it_NR, tol_NR, no_guess, &
        &no_plots, input_i, output_i, eq_i, min_alpha, max_alpha, n_alpha, &
        &max_deriv, mu_0, calc_mesh_style, iu, EV_style, n_procs_per_alpha, &
        &n_procs, MPI_Comm_groups, MPI_Comm_masters, glb_rank, glb_n_procs, &
        &grp_rank, grp_n_procs, grp_nr, n_groups, output_name, next_job, &
        &next_job_win, plot_jq, n_sol_requested, min_n_r_X, min_r_X, max_r_X, &
        &nyq_fac, max_n_plots, alpha_job_nr, use_pol_flux, plot_grid, &
        &output_style, eq_style, eq_name, plot_dir, data_dir, script_dir, &
        &spline_type, plot_flux_q

    ! technical variables
    integer, parameter :: dp=kind(1.d0)                                         ! double precision
    integer, parameter :: qp = selected_real_kind (32)                          ! quadruple precision
    integer, parameter :: max_str_ln = 100                                      ! maximum length of filenames
    integer, parameter :: n_seq_0 = 10                                          ! start of index of file numbers for opening
    integer, parameter :: max_args = 10                                         ! maximum number of input arguments
    integer, parameter :: max_opts = 8                                          ! maximum number of options in input arguments
    integer, parameter, dimension(3) :: max_deriv = [2,2,2]                     ! highest derivatives that are tabulated for VMEC amplitudes R, Z, L in r,theta,zeta)
    character(len=max_str_ln) :: prog_name = 'PB3D'                             ! name of program, used for info
    character(len=max_str_ln) :: output_name                                    ! will hold name of output file
    integer :: n_sol_requested                                                  ! how many solutions requested

    ! MPI variables
    integer :: n_procs_per_alpha                                                ! how many processors are used per field line alpha
    integer, allocatable :: n_procs(:)                                          ! hwo many processors per group of alpha
    integer :: MPI_Comm_groups                                                  ! communicator for the groups of alpha
    integer :: MPI_Comm_masters                                                 ! communicator for the masters of the groups of alpha
    integer :: glb_rank                                                         ! global MPI rank
    integer :: glb_n_procs                                                      ! global nr. MPI processes
    integer :: grp_rank                                                         ! alpha group MPI rank
    integer :: grp_n_procs                                                      ! alpha gropu nr. MPI processes
    integer :: grp_nr                                                           ! group nr.
    integer :: n_groups                                                         ! nr. of groups
    integer :: next_job                                                         ! next job to be done
    integer :: next_job_win                                                     ! window to next_job

    ! physical and mathematical variables
    real(dp), parameter :: pi=4_dp*datan(1.0_dp)                                ! pi
    real(dp), parameter :: mu_0 = 4E-7_dp*pi                                    ! permeability of free space
    complex(dp), parameter :: iu = (0,1)                                        ! complex unit

    ! concerning runtime
    integer :: minim_style                                                      ! determines the method used for minimization
        ! 1 [def] : Euler-Lagrange min., finite diff and Richardson's method
    logical :: ltest                                                            ! whether or not to call the testing routines
    integer :: calc_mesh_style                                                  ! how equilibrium mesh is calculated
    integer :: EV_style                                                         ! determines the method used for solving an EV problem
    integer :: eq_style                                                         ! either 1 (VMEC) or 2 (HELENA)
    integer :: alpha_job_nr                                                     ! which alpha job is being calculated
    logical :: use_pol_flux                                                     ! whether or not the poloidal flux is used as radial variable
    logical :: plot_jq                                                          ! whether to plot the q-profile with nq-m = 0 or iota-profile with n-iotam = 0 (only global master)
    logical :: plot_grid                                                        ! whether to plot the grid in real coordinates (only group masters)
    logical :: plot_flux_q                                                      ! whether to plot flux quantities in real coordinates (only global master)
    
    ! concerning Richardson extrapolation
    integer :: max_it_r                                                         ! number of levels for Richardson extrapolation
    real(dp) :: tol_r                                                           ! tolerance for Richardson extrapolation
    logical :: no_guess = .false.                                               ! disable guessing Eigenfunction from previous level of Richardson

    ! concerning finding the magnetic field lines
    integer :: max_it_NR                                                        ! maximum number of Newton-Rhapson iterations
    real(dp) :: tol_NR                                                          ! tolerance for Newton-Rhapson
    integer :: nyq_fac                                                          ! Nyquist factor to avoid aliasing in perturbation integrals

    ! input / output
    integer :: input_i                                                          ! file number of input file
    integer :: eq_i                                                             ! file number of equilibrium file from VMEC or HELENA
    character(len=max_str_ln) :: eq_name                                        ! name of equilibrium file
    integer :: output_i                                                         ! file number of output file
    integer :: max_n_plots                                                      ! max. nr. of modes for which to output a plot
    logical :: no_plots = .false.                                               ! true if no plots should be made
    integer :: output_style                                                     ! style of output (GNUPlot, HDF5, ...)
    character(len=5) :: plot_dir = 'Plots'                                      ! directory where to save plots
    character(len=7) :: script_dir = 'Scripts'                                  ! directory where to save scripts for plots
    character(len=4) :: data_dir = 'Data'                                       ! directory where to save data for plots
    
    ! concerning the various field lines for which to do the calculations
    integer :: n_alpha                                                          ! how many field lines
    real(dp) :: min_alpha, max_alpha                                            ! min. and max. value for alpha, should be (0...2pi)
    
    ! concerning the perturbation grid
    integer :: min_n_r_X                                                        ! min. of n_r_X (e.g. first value in Richardson loop)
    real(dp) :: min_r_X, max_r_X                                                ! min. and max. normal coord. for pert. (either pol. or tor., depending on use_pol_flux from VMEC)
    
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
