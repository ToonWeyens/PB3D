!------------------------------------------------------------------------------!
!> Full-stack tests of the spline interpolation wrapper
!! (spline_utilities.spline, a PB3D-flavored interface to PSPLINE/EZspline),
!! converted from the interactive legacy check (Modules/test.f90,
!! test_splines) into automated assertions with implementation-independent
!! references:
!!
!!  - reproduction of polynomials that lie exactly in the interpolation
!!    space: linear data for order 1, and a cubic polynomial (with exact
!!    prescribed endpoint derivatives) for order 3, including all
!!    derivatives; the extrapolation outside the domain is a *quadratic*
!!    Taylor extension from the boundary (see calc_extrap in
!!    spline_utilities), so it is pinned with a quadratic polynomial;
!!  - convergence at the expected rate on sin(2 pi x) for all three orders
!!    (1: linear, 2: Akima Hermite, 3: cubic);
!!  - periodic boundary conditions on a full sin period;
!!  - refusal to extrapolate when extrap is not set (this check lives
!!    outside ldebug, so it is active in all builds).
!------------------------------------------------------------------------------!
module test_splines
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str, i2str
    use spline_utilities, only: spline

    implicit none
    private
    public collect_splines

contains
    !> Collect the tests.
    subroutine collect_splines(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("linear_exact_on_linear", test_linear_exact), &
            &new_unittest("cubic_exact_on_cubic", test_cubic_exact), &
            &new_unittest("cubic_extrapolation", test_cubic_extrap), &
            &new_unittest("sin_convergence", test_sin_convergence), &
            &new_unittest("periodic_bc_sin", test_periodic_sin), &
            &new_unittest("extrapolation_refused", test_extrap_refused) &
            &]
    end subroutine collect_splines

    !--------------------------------------------------------------------------
    ! helpers
    !--------------------------------------------------------------------------

    !> The cubic test polynomial and its derivatives.
    pure function cubic_poly(x,deriv) result(res)
        real(dp), intent(in) :: x
        integer, intent(in) :: deriv
        real(dp) :: res

        select case (deriv)
            case (0)
                res = x**3 - 2._dp*x**2 + 3._dp*x - 1._dp
            case (1)
                res = 3._dp*x**2 - 4._dp*x + 3._dp
            case (2)
                res = 6._dp*x - 4._dp
            case (3)
                res = 6._dp
            case default
                res = 0._dp
        end select
    end function cubic_poly

    !> Maximum error of the spline of order ord and derivative deriv against
    !! sin(2 pi x) sampled on nx points, evaluated on a fine interior grid,
    !! normalized by the amplitude (2 pi)^deriv of the derivative.
    integer function sin_err(nx,ord,deriv,err) result(ierr)
        integer, intent(in) :: nx, ord, deriv                                   ! knots, order, derivative
        real(dp), intent(out) :: err                                            ! normalized maximum error

        integer, parameter :: nx_int = 301                                      ! evaluation points
        integer :: kd                                                           ! counter
        real(dp) :: x(nx), y(nx)                                                ! knots
        real(dp) :: x_int(nx_int), y_int(nx_int), y_ref(nx_int)                 ! evaluation

        do kd = 1,nx
            x(kd) = (kd-1._dp)/(nx-1._dp)
            y(kd) = sin(2._dp*pi*x(kd))
        end do
        do kd = 1,nx_int
            x_int(kd) = (kd-1._dp)/(nx_int-1._dp)
            select case (deriv)
                case (0)
                    y_ref(kd) = sin(2._dp*pi*x_int(kd))
                case (1)
                    y_ref(kd) = 2._dp*pi*cos(2._dp*pi*x_int(kd))
            end select
        end do

        ierr = spline(x,y,x_int,y_int,ord=ord,deriv=deriv)
        if (ierr.ne.0) return

        err = maxval(abs(y_int-y_ref))/(2._dp*pi)**deriv
    end function sin_err

    !--------------------------------------------------------------------------
    ! tests
    !--------------------------------------------------------------------------

    !> Order 1 reproduces linear data exactly, including the derivative.
    subroutine test_linear_exact(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 7, nx_int = 50                               ! knots, evaluation points
        integer :: ierr, kd                                                     ! error status, counter
        real(dp) :: x(nx), y(nx)                                                ! knots
        real(dp) :: x_int(nx_int), y_int(nx_int)                                ! evaluation

        do kd = 1,nx
            x(kd) = 2._dp*(kd-1._dp)/(nx-1._dp)                                 ! 0..2
            y(kd) = 2._dp*x(kd) + 1._dp
        end do
        do kd = 1,nx_int
            x_int(kd) = 2._dp*(kd-1._dp)/(nx_int-1._dp)
        end do

        ierr = spline(x,y,x_int,y_int,ord=1,deriv=0)
        call check(error, ierr, 0, 'spline failed')
        if (allocated(error)) return
        call check(error, maxval(abs(y_int-(2._dp*x_int+1._dp))).lt.1.e-12_dp,&
            &'order 1 not exact on linear data')
        if (allocated(error)) return

        ierr = spline(x,y,x_int,y_int,ord=1,deriv=1)
        call check(error, ierr, 0, 'spline failed')
        if (allocated(error)) return
        call check(error, maxval(abs(y_int-2._dp)).lt.1.e-11_dp, &
            &'order 1 derivative not exact on linear data')
        if (allocated(error)) return
    end subroutine test_linear_exact

    !> Order 3 with exact prescribed endpoint first derivatives reproduces a
    !! cubic polynomial exactly, for all derivatives 0..3 (the interpolant
    !! is unique in the space the data lie in).
    subroutine test_cubic_exact(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 12, nx_int = 100                             ! knots, evaluation points
        integer :: ierr, kd, deriv                                              ! error status, counters
        real(dp) :: x(nx), y(nx)                                                ! knots
        real(dp) :: x_int(nx_int), y_int(nx_int), y_ref(nx_int)                 ! evaluation
        real(dp) :: err                                                         ! maximum error

        do kd = 1,nx
            x(kd) = 2._dp*(kd-1._dp)/(nx-1._dp)                                 ! 0..2
            y(kd) = cubic_poly(x(kd),0)
        end do
        do kd = 1,nx_int
            x_int(kd) = 2._dp*(kd-1._dp)/(nx_int-1._dp)
        end do

        do deriv = 0,3
            ierr = spline(x,y,x_int,y_int,ord=3,deriv=deriv,bcs=[1,1],&
                &bcs_val=[cubic_poly(x(1),1),cubic_poly(x(nx),1)])
            call check(error, ierr, 0, 'spline failed for derivative '//&
                &trim(i2str(deriv)))
            if (allocated(error)) return
            do kd = 1,nx_int
                y_ref(kd) = cubic_poly(x_int(kd),deriv)
            end do
            err = maxval(abs(y_int-y_ref))/maxval(abs(y_ref))
            call check(error, err.lt.1.e-9_dp, &
                &'cubic spline not exact on cubic polynomial, derivative '//&
                &trim(i2str(deriv))//': '//trim(r2str(err)))
            if (allocated(error)) return
        end do
    end subroutine test_cubic_exact

    !> The extrapolation is the quadratic Taylor extension from the boundary
    !! (calc_extrap in spline_utilities), so it stays exact for a quadratic
    !! polynomial with exact prescribed endpoint first derivatives.
    subroutine test_cubic_extrap(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 12, nx_int = 60                              ! knots, evaluation points
        integer :: ierr, kd                                                     ! error status, counter
        real(dp) :: x(nx), y(nx)                                                ! knots
        real(dp) :: x_int(nx_int), y_int(nx_int), y_ref(nx_int)                 ! evaluation
        real(dp) :: err                                                         ! maximum error

        do kd = 1,nx
            x(kd) = 2._dp*(kd-1._dp)/(nx-1._dp)                                 ! 0..2
            y(kd) = x(kd)**2 - 3._dp*x(kd) + 2._dp
        end do
        do kd = 1,nx_int
            x_int(kd) = -0.5_dp + 3._dp*(kd-1._dp)/(nx_int-1._dp)               ! -0.5..2.5: 0.5 extrapolated on both sides
            y_ref(kd) = x_int(kd)**2 - 3._dp*x_int(kd) + 2._dp
        end do

        ierr = spline(x,y,x_int,y_int,ord=3,deriv=0,bcs=[1,1],&
            &bcs_val=[2._dp*x(1)-3._dp,2._dp*x(nx)-3._dp],extrap=.true.)
        call check(error, ierr, 0, 'spline failed')
        if (allocated(error)) return

        err = maxval(abs(y_int-y_ref))/maxval(abs(y_ref))
        call check(error, err.lt.1.e-9_dp, &
            &'quadratic extrapolation not exact on quadratic polynomial: '//&
            &trim(r2str(err)))
        if (allocated(error)) return
    end subroutine test_cubic_extrap

    !> Convergence on sin(2 pi x) when doubling the number of knots, for all
    !! three orders: the error has to drop at least by the conservative
    !! factors below (theoretical rates: 4 for linear, 8 for Akima, 16 for
    !! cubic), and be small in absolute terms at the finer resolution.
    subroutine test_sin_convergence(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp), parameter :: ratio_min(3) = [3._dp,5._dp,10._dp]              ! minimum error-reduction factors
        real(dp), parameter :: err_max(3) = [4.e-3_dp,4.e-4_dp,3.e-5_dp]        ! maximum normalized errors at nx = 41 (ord 1 theory: h^2 |f''| / 8 = 3.1e-3)

        integer :: ierr, ord, deriv                                             ! error status, counters
        real(dp) :: err_lo, err_hi                                              ! errors at nx = 21 and 41

        do ord = 1,3
            do deriv = 0,1
                ierr = sin_err(21,ord,deriv,err_lo)
                call check(error, ierr, 0, 'spline failed')
                if (allocated(error)) return
                ierr = sin_err(41,ord,deriv,err_hi)
                call check(error, ierr, 0, 'spline failed')
                if (allocated(error)) return
                write(error_unit,'(A,I2,A,I2,A,ES10.3,A,ES10.3)') &
                    &'   [spline convergence] ord ',ord,', deriv ',deriv,&
                    &': err 21 = ',err_lo,', err 41 = ',err_hi
                if (deriv.eq.0) then
                    call check(error, err_hi.lt.err_max(ord), &
                        &'error too large for order '//trim(i2str(ord))//&
                        &': '//trim(r2str(err_hi)))
                    if (allocated(error)) return
                    call check(error, err_hi.lt.err_lo/ratio_min(ord), &
                        &'no convergence for order '//trim(i2str(ord))//&
                        &': '//trim(r2str(err_lo))//' -> '//trim(r2str(err_hi)))
                    if (allocated(error)) return
                else
                    call check(error, err_hi.lt.err_lo, &                       ! derivative errors converge more slowly; just require decrease
                        &'derivative error grows for order '//&
                        &trim(i2str(ord))//': '//trim(r2str(err_lo))//' -> '//&
                        &trim(r2str(err_hi)))
                    if (allocated(error)) return
                end if
            end do
        end do
    end subroutine test_sin_convergence

    !> Periodic boundary conditions on a full sin period: the periodic cubic
    !! spline has to be accurate up to the second derivative and exactly
    !! periodic at the seam.
    subroutine test_periodic_sin(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 21, nx_int = 301                             ! knots, evaluation points
        real(dp), parameter :: tol(0:2) = [1.e-4_dp,1.e-3_dp,1.e-2_dp]          ! normalized tolerances per derivative

        integer :: ierr, kd, deriv                                              ! error status, counters
        real(dp) :: x(nx), y(nx)                                                ! knots
        real(dp) :: x_int(nx_int), y_int(nx_int), y_ref(nx_int)                 ! evaluation
        real(dp) :: err                                                         ! maximum error
        real(dp) :: y_seam(2)                                                   ! values at both ends

        do kd = 1,nx
            x(kd) = (kd-1._dp)/(nx-1._dp)
            y(kd) = sin(2._dp*pi*x(kd))
        end do
        do kd = 1,nx_int
            x_int(kd) = (kd-1._dp)/(nx_int-1._dp)
        end do

        do deriv = 0,2
            ierr = spline(x,y,x_int,y_int,ord=3,deriv=deriv,bcs=[-1,-1])
            call check(error, ierr, 0, 'spline failed for derivative '//&
                &trim(i2str(deriv)))
            if (allocated(error)) return
            do kd = 1,nx_int
                y_ref(kd) = (2._dp*pi)**deriv*&
                    &sin(2._dp*pi*x_int(kd)+deriv*pi/2._dp)
            end do
            err = maxval(abs(y_int-y_ref))/(2._dp*pi)**deriv
            write(error_unit,'(A,I2,A,ES10.3)') &
                &'   [periodic spline] deriv ',deriv,': err = ',err
            call check(error, err.lt.tol(deriv), &
                &'periodic spline error too large for derivative '//&
                &trim(i2str(deriv))//': '//trim(r2str(err)))
            if (allocated(error)) return

            ! seam periodicity
            y_seam(1) = y_int(1)
            y_seam(2) = y_int(nx_int)
            call check(error, abs(y_seam(2)-y_seam(1)).lt.1.e-10_dp*&
                &(2._dp*pi)**deriv, &
                &'periodic spline not periodic at the seam for derivative '//&
                &trim(i2str(deriv)))
            if (allocated(error)) return
        end do
    end subroutine test_periodic_sin

    !> Evaluation outside the domain without extrap has to fail (this check
    !! is active in all builds, not only ldebug).
    subroutine test_extrap_refused(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nx = 7                                            ! knots
        integer :: ierr, kd                                                     ! error status, counter
        real(dp) :: x(nx), y(nx)                                                ! knots
        real(dp) :: x_int(3), y_int(3)                                          ! evaluation with an exterior point

        do kd = 1,nx
            x(kd) = (kd-1._dp)/(nx-1._dp)
            y(kd) = x(kd)**2
        end do
        x_int = [0.2_dp,0.8_dp,1.5_dp]                                          ! last point outside

        ierr = spline(x,y,x_int,y_int,ord=3,deriv=0)
        call check(error, ierr.ne.0, &
            &'extrapolation without extrap=.true. was not refused')
        if (allocated(error)) return
    end subroutine test_extrap_refused
end module test_splines
