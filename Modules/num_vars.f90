!------------------------------------------------------------------------------!
!   numerical variables used by most other modules                             !
!------------------------------------------------------------------------------!
module num_vars
    implicit none
    private
    public max_it, dp, qp, style, max_str_ln, n_seq_0, max_args, &
        &max_opts, prog_name, max_it_r, ltest, pi, max_it_NR, tol_NR, &
        &input_i, output_i, VMEC_i, min_alpha, max_alpha, n_alpha, &
        &theta_var_along_B, max_deriv, mu_0, calc_mesh_style, iu, EV_style, &
        &n_procs_per_alpha, n_procs, MPI_comm_groups, glob_rank, glob_n_procs

    ! technical variables
    integer, parameter :: dp=kind(1.d0)                                         ! double precision
    integer, parameter :: qp = selected_real_kind (32)                          ! quadruple precision
    integer, parameter :: max_str_ln = 100                                      ! maximum length of filenames
    integer, parameter :: n_seq_0 = 10                                          ! start of index of file numbers for opening
    integer, parameter :: max_args = 10                                         ! maximum number of input arguments
    integer, parameter :: max_opts = 8                                          ! maximum number of options in input arguments
    integer, parameter, dimension(3) :: max_deriv = [4,4,4]                     ! highest derivatives that are tabulated for VMEC amplitudes R, Z, L in theta,zeta,r)
    character(len=max_str_ln) :: prog_name = 'PB3D'                             ! name of program, used for info
    integer :: n_procs_per_alpha                                                ! how many processors are used per field line alpha
    integer, allocatable :: n_procs(:)                                          ! hwo many processors per group of alpha
    integer :: MPI_comm_groups                                                  ! communicator for the groups of alpha
    integer :: glob_rank                                                        ! global MPI rank
    integer :: glob_n_procs                                                     ! global nr. MPI processes

    ! considering runtime
    integer :: max_it_r                                                         ! number of levels for Richardson's extrapolation
    integer :: max_it                                                           ! max. nr. iterations
    integer :: style                                                            ! determines the method used for minimization
        ! 1 [def] : Euler-Lagrange min., finite diff and Richardson's method
    logical :: ltest                                                            ! whether or not to call the testing routines
    integer :: calc_mesh_style
    integer :: EV_style                                                         ! determines the method used for solving an EV problem

    ! physical and mathematical variables
    real(dp), parameter :: pi=4_dp*datan(1.0_dp)                                ! pi
    real(dp), parameter :: mu_0 = 4E-7_dp*pi                                    ! permeability of free space
    complex(dp), parameter :: iu = (0,1)                                        ! complex unit

    ! considering finding the magnetic field lines
    integer :: max_it_NR                                                        ! maximum number of Newton-Rhapson iterations
    real(dp) :: tol_NR                                                          ! tolerance for Newton-Rhapson
    logical :: theta_var_along_B                                                ! true if theta is used as the parallel variable

    ! input / output
    integer :: input_i                                                          ! will hold the file number of input file
    integer :: VMEC_i                                                           ! will hold the file number of VMEC file
    integer :: output_i                                                         ! will hold the file number of output file
    
    ! considering the various field lines for which to do the calculations
    integer :: n_alpha                                                          ! how many field lines
    real(dp) :: min_alpha, max_alpha                                            ! min. and max. value for alpha, should be (0...2pi)
end module num_vars
