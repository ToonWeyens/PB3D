!------------------------------------------------------------------------------!
!> Driver for the PB3D unit tests, built on the test-drive framework.
!!
!! Usage:
!!  - no arguments:         run all test suites
!!  - one argument:         run only that suite (e.g. "dtorh")
!!  - two arguments:        run one test of one suite (e.g. "dtorh reference")
!!
!! Each suite is also registered individually with CTest, so the normal entry
!! point is simply \c ctest from the build directory.
!------------------------------------------------------------------------------!
program pb3d_tester
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
        &select_suite, run_selected, get_argument
    use test_str_utilities, only: collect_str_utilities
    use test_files_utilities, only: collect_files_utilities
    use test_dtorh, only: collect_dtorh
#ifdef PB3D_TEST_NUM
    use test_num_utilities, only: collect_num_utilities
    use test_X_utilities, only: collect_X_utilities
#endif
#ifdef PB3D_TEST_NUM_OPS
    use test_num_ops, only: collect_num_ops
#endif
    use test_PB3D_utilities, only: collect_PB3D_utilities
    use num_vars, only: dp, rank, n_procs, prog_name, max_it_zero, tol_zero, &
        &max_nr_backtracks_HH
    use messages, only: init_output

    implicit none

    ! local variables
    integer :: stat                                                             ! cumulative failure count
    integer :: is                                                               ! suite index
    character(len=:), allocatable :: suite_name, test_name                      ! command-line selection
    type(testsuite_type), allocatable :: testsuites(:)                          ! all test suites
    character(len=*), parameter :: fmt = '("#", *(1x, a))'                      ! output format

    ! minimal PB3D global state so that library routines (messages, ...) work
    ! outside the full PB3D/POST startup sequence
    rank = 0
    n_procs = 1
    prog_name = 'TEST'
    max_it_zero = 100                                                           ! zero-finder settings, normally set during input
    tol_zero = 1.e-10_dp                                                        ! processing (see input_ops)
    max_nr_backtracks_HH = 20
    call init_output()

    stat = 0

    testsuites = [ &
        &new_testsuite("str_utilities", collect_str_utilities), &
        &new_testsuite("files_utilities", collect_files_utilities), &
#ifdef PB3D_TEST_NUM
        &new_testsuite("num_utilities", collect_num_utilities), &
        &new_testsuite("X_utilities", collect_X_utilities), &
#endif
#ifdef PB3D_TEST_NUM_OPS
        &new_testsuite("num_ops", collect_num_ops), &
#endif
        &new_testsuite("PB3D_utilities", collect_PB3D_utilities), &
        &new_testsuite("dtorh", collect_dtorh) &
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

    if (stat.gt.0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop 1
    end if
end program pb3d_tester
