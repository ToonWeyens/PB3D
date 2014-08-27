!------------------------------------------------------------------------------!
!   numerical variables used by most other modules                             !
!------------------------------------------------------------------------------!
module num_vars
    implicit none
    private
    public dp, qp, style, max_str_ln, n_seq_0, max_args, max_opts, prog_name, &
        &max_it_r, tol_r, ltest, pi, max_it_NR, tol_NR, &
        &input_i, output_i, VMEC_i, min_alpha, max_alpha, n_alpha, &
        &theta_var_along_B, max_deriv, mu_0, calc_mesh_style, iu, EV_style, &
        &n_procs_per_alpha, n_procs, MPI_Comm_groups, MPI_Comm_masters, &
        &glob_rank, glob_n_procs, group_rank, group_n_procs, group_nr, &
        &n_groups,  output_name, next_job, next_job_win, plot_q, &
        &n_sol_requested, min_n_r_X, min_r_X, max_r_X, nyq_fac, reuse_r

    ! technical variables
    integer, parameter :: dp=kind(1.d0)                                         ! double precision
    integer, parameter :: qp = selected_real_kind (32)                          ! quadruple precision
    integer, parameter :: max_str_ln = 100                                      ! maximum length of filenames
    integer, parameter :: n_seq_0 = 10                                          ! start of index of file numbers for opening
    integer, parameter :: max_args = 10                                         ! maximum number of input arguments
    integer, parameter :: max_opts = 8                                          ! maximum number of options in input arguments
    integer, parameter, dimension(3) :: max_deriv = [2,2,2]                     ! highest derivatives that are tabulated for VMEC amplitudes R, Z, L in theta,zeta,r)
    character(len=max_str_ln) :: prog_name = 'PB3D'                             ! name of program, used for info
    character(len=max_str_ln) :: output_name                                    ! will hold name of output file
    logical :: plot_q                                                           ! whether to plot the q-profile with nq-m = 0
    integer :: n_sol_requested                                                  ! how many solutions requested

    ! MPI variables
    integer :: n_procs_per_alpha                                                ! how many processors are used per field line alpha
    integer, allocatable :: n_procs(:)                                          ! hwo many processors per group of alpha
    integer :: MPI_Comm_groups                                                  ! communicator for the groups of alpha
    integer :: MPI_Comm_masters                                                 ! communicator for the masters of the groups of alpha
    integer :: glob_rank                                                        ! global MPI rank
    integer :: glob_n_procs                                                     ! global nr. MPI processes
    integer :: group_rank                                                       ! alpha group MPI rank
    integer :: group_n_procs                                                    ! alpha gropu nr. MPI processes
    integer :: group_nr                                                         ! group nr.
    integer :: n_groups                                                         ! nr. of groups
    integer :: next_job                                                         ! next job to be done
    integer :: next_job_win                                                     ! window to next_job

    ! physical and mathematical variables
    real(dp), parameter :: pi=4_dp*datan(1.0_dp)                                ! pi
    real(dp), parameter :: mu_0 = 4E-7_dp*pi                                    ! permeability of free space
    complex(dp), parameter :: iu = (0,1)                                        ! complex unit

    ! considering runtime
    integer :: style                                                            ! determines the method used for minimization
        ! 1 [def] : Euler-Lagrange min., finite diff and Richardson's method
    logical :: ltest                                                            ! whether or not to call the testing routines
    integer :: calc_mesh_style                                                  ! how equilibrium mesh is calculated
    integer :: EV_style                                                         ! determines the method used for solving an EV problem
    
    ! considering Richardson extrapolation
    integer :: max_it_r                                                         ! number of levels for Richardson extrapolation
    real(dp) :: tol_r                                                           ! tolerance for Richardson extrapolation
    logical :: reuse_r                                                          ! whether to reuse the matrices A and B from previous Richardson level

    ! considering finding the magnetic field lines
    integer :: max_it_NR                                                        ! maximum number of Newton-Rhapson iterations
    real(dp) :: tol_NR                                                          ! tolerance for Newton-Rhapson
    logical :: theta_var_along_B                                                ! true if theta is used as the parallel variable
    integer, parameter :: nyq_fac = 5                                           ! Nyquist factor to avoid aliasing in perturbation integrals

    ! input / output
    integer :: input_i                                                          ! file number of input file
    integer :: VMEC_i                                                           ! file number of VMEC file
    integer :: output_i                                                         ! file number of output file
    
    ! considering the various field lines for which to do the calculations
    integer :: n_alpha                                                          ! how many field lines
    real(dp) :: min_alpha, max_alpha                                            ! min. and max. value for alpha, should be (0...2pi)
    
    ! considering the perturbation grid
    integer :: min_n_r_X                                                        ! min. of n_r_X (e.g. first value in Richardson loop)
    real(dp) :: min_r_X, max_r_X                                                ! mimimum and maximum radius for perturbations
end module num_vars
