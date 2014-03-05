module num_vars
    implicit none
    private
    public max_it, dp, qp, style, max_str_ln, n_seq_0, max_args, &
        &max_opts, prog_name, max_r, ltest, pi, min_theta, max_theta, &
        &min_zeta, max_zeta, n_theta, n_zeta

    ! technical variables
    integer, parameter :: dp=kind(1.d0)                                         ! double precision
    integer, parameter :: qp = selected_real_kind (32)                        ! quadruple precision
    integer, parameter :: max_str_ln = 80                                       ! maximum length of filenames
    integer, parameter :: n_seq_0 = 10                                          ! start of index of file numbers for opening
    integer, parameter :: max_args = 10                                         ! maximum number of input arguments
    integer, parameter :: max_opts = 8                                          ! maximum number of options in input arguments
    character(len=max_str_ln) :: prog_name = 'PB3D'                             ! name of program, used for info

    ! considering runtime
    integer :: max_r                                                            ! number of levels for Richardson's extrapolation
    integer :: max_it                                                           ! max. nr. iterations
    integer :: style                                                            ! determines the method used for minimization
        ! 1 [def] : Euler-Lagrange min., finite diff and Richardson's method
    logical :: ltest                                                            ! whether or not to call the testing routines

    ! global variables
    real(dp) :: pi=4_dp*datan(1.0_dp)                                           ! pi
    real(dp) :: min_theta, max_theta, min_zeta, max_zeta                        ! angular values at which table of metrics is calculated
    integer :: n_theta, n_zeta                                                  ! number of points in table

!contains
end module num_vars
