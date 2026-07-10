!------------------------------------------------------------------------------!
!> Unit tests of the zero finders in num_ops (Householder with backtracking,
!! Zhang bracketing) on functions with known roots.
!------------------------------------------------------------------------------!
module test_num_ops
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str
    use num_ops, only: calc_zero_HH, calc_zero_Zhang

    implicit none
    private
    public collect_num_ops

contains
    !> Collect the tests.
    subroutine collect_num_ops(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("zero_HH", test_zero_HH), &
            &new_unittest("zero_Zhang", test_zero_Zhang) &
            &]
    end subroutine collect_num_ops

    !> Householder zero finder on f = x^2 - 2 (root sqrt(2)), orders 1-3.
    subroutine test_zero_HH(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ord                                                          ! Householder order
        real(dp) :: zero                                                        ! found zero
        character(len=:), allocatable :: msg                                    ! error message

        do ord = 1,3
            zero = 0._dp
            msg = trim(calc_zero_HH(zero,f_sq2,ord,1.0_dp))
            call check(error, len_trim(msg), 0, 'calc_zero_HH returned: '//msg)
            if (allocated(error)) return
            call check(error, zero, sqrt(2._dp), thr=1.e-8_dp, rel=.true., &
                &message='HH order '//trim(r2str(real(ord,dp)))//' zero: '//&
                &trim(r2str(zero)))
            if (allocated(error)) return
        end do
    end subroutine test_zero_HH

    !> Zhang bracketing zero finder on f = cos (root pi/2).
    subroutine test_zero_Zhang(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp) :: zero                                                        ! found zero
        character(len=:), allocatable :: msg                                    ! error message

        zero = 0._dp
        msg = trim(calc_zero_Zhang(zero,f_cos,[1.0_dp,2.0_dp]))
        call check(error, len_trim(msg), 0, 'calc_zero_Zhang returned: '//msg)
        if (allocated(error)) return
        call check(error, zero, pi/2._dp, thr=1.e-8_dp, rel=.true., &
            &message='Zhang zero: '//trim(r2str(zero)))
        if (allocated(error)) return
    end subroutine test_zero_Zhang

    !> f = x^2 - 2 and derivatives, in the calc_zero_HH callback form.
    function f_sq2(x,ord)
        real(dp), intent(in) :: x
        integer, intent(in) :: ord
        real(dp) :: f_sq2

        select case (ord)
            case (0)
                f_sq2 = x**2 - 2._dp
            case (1)
                f_sq2 = 2._dp*x
            case (2)
                f_sq2 = 2._dp
            case default
                f_sq2 = 0._dp
        end select
    end function f_sq2

    !> f = cos, in the calc_zero_Zhang callback form.
    function f_cos(x)
        real(dp), intent(in) :: x
        real(dp) :: f_cos

        f_cos = cos(x)
    end function f_cos
end module test_num_ops
