!------------------------------------------------------------------------------!
!> Driver for the PB3D full-stack tests, built on the test-drive framework.
!!
!! In contrast to the unit tests (tests/unit), these tests link against the
!! complete pb3d_modules library, and therefore require the full dependency
!! stack (MPI, ScaLAPACK, STRUMPACK, ...). They are intended to be run with a
!! single MPI process; the executable initializes MPI itself, so it can simply
!! be run directly (MPI singleton mode).
!!
!! Usage:
!!  - no arguments:         run all test suites
!!  - one argument:         run only that suite (e.g. "vac_kernels")
!!  - two arguments:        run one test of one suite
!------------------------------------------------------------------------------!
program pb3d_fullstack_tester
    use, intrinsic :: iso_fortran_env, only: error_unit
    use MPI
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
        &select_suite, run_selected, get_argument
    use test_vac_kernels, only: collect_vac_kernels
    use test_vac_greens, only: collect_vac_greens
    use test_vac_3d, only: collect_vac_3d
    use num_vars, only: dp, rank, n_procs, prog_name, max_it_zero, tol_zero, &
        &max_nr_backtracks_HH, rich_restart_lvl
    use rich_vars, only: rich_lvl
    use X_vars, only: n_mod_X
    use messages, only: init_output

    implicit none

    ! local variables
    integer :: stat                                                             ! cumulative failure count
    integer :: is                                                               ! suite index
    integer :: ierr                                                             ! error variable
    character(len=:), allocatable :: suite_name, test_name                      ! command-line selection
    type(testsuite_type), allocatable :: testsuites(:)                          ! all test suites
    character(len=*), parameter :: fmt = '("#", *(1x, a))'                      ! output format

    ! minimal PB3D global state so that library routines (messages, BLACS
    ! grids, ...) work outside the full PB3D/POST startup sequence
    call MPI_init(ierr)
    if (ierr.ne.0) error stop 'MPI_init failed'
    call MPI_Comm_rank(MPI_Comm_world,rank,ierr)
    call MPI_Comm_size(MPI_Comm_world,n_procs,ierr)
    prog_name = 'TEST'
    n_mod_X = 1                                                                 ! vac%res is allocated with this size; overridden by the response test
    max_it_zero = 100                                                           ! zero-finder settings, normally set during input
    tol_zero = 1.e-10_dp                                                        ! processing (see input_ops); calc_GH_2 needs them
    max_nr_backtracks_HH = 20                                                   ! for its singularity tolerance
    rich_lvl = 1                                                                ! no Richardson extrapolation in the tests
    rich_restart_lvl = 1
    call init_output()

    stat = 0

    testsuites = [ &
        &new_testsuite("vac_kernels", collect_vac_kernels), &
        &new_testsuite("vac_greens", collect_vac_greens), &
        &new_testsuite("vac_3d", collect_vac_3d) &
        &]

    call get_argument(1, suite_name)
    call get_argument(2, test_name)

    if (allocated(suite_name)) then
        is = select_suite(testsuites, suite_name)
        if (is.gt.0 .and. is.le.size(testsuites)) then
            if (allocated(test_name)) then
                write(error_unit, fmt) "Suite:", testsuites(is)%name
                call run_selected(testsuites(is)%collect, test_name, &
                    &error_unit, stat)
                if (stat.lt.0) error stop 1
            else
                write(error_unit, fmt) "Testing:", testsuites(is)%name
                call run_testsuite(testsuites(is)%collect, error_unit, stat)
            end if
        else
            write(error_unit, fmt) "Available testsuites"
            do is = 1, size(testsuites)
                write(error_unit, fmt) "-", testsuites(is)%name
            end do
            error stop 1
        end if
    else
        do is = 1, size(testsuites)
            write(error_unit, fmt) "Testing:", testsuites(is)%name
            call run_testsuite(testsuites(is)%collect, error_unit, stat)
        end do
    end if

    call MPI_finalize(ierr)

    if (stat.gt.0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop 1
    end if
end program pb3d_fullstack_tester
