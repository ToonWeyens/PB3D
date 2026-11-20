!------------------------------------------------------------------------------!
!> Full-stack tests of the volume integral (grid_utilities.calc_int_vol),
!! converted from the interactive legacy check (Modules/test.f90,
!! test_calc_int_vol) into automated assertions.
!!
!! The reference is the analytical integral of
!!    f = 1 - r^2 + i cos(theta)
!! over a torus with geometric axis R_0 and minor radius 1, with the flux
!! Jacobian J = r (R_0 + r cos(theta)):
!!    int f J dtheta dzeta dr = R_0 pi^2 + i 2 pi^2 / 3 .
!! The grid follows PB3D's convention (theta varies in dimension 1, zeta in
!! dimension 2, the normal variable in dimension 3, angles 0..2 pi inclusive,
!! r in 0..1).
!!
!! Covered: absolute agreement at a moderate resolution, second-order
!! convergence (the cell-averaged rule is trapezoidal for independent
!! coordinates), and the behavior for a single point in dimension 2: the
!! implementation substitutes the full turn 2 pi for the missing angular
!! extent (transf_J is set to 2 pi for a singleton dimension), so an
!! axisymmetric integrand yields the complete toroidal integral.
!------------------------------------------------------------------------------!
module test_calc_int_vol
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi, iu
    use str_utilities, only: r2str, i2str
    use grid_utilities, only: calc_eqd_grid, calc_int_vol

    implicit none
    private
    public collect_calc_int_vol

    real(dp), parameter :: R_0 = 2._dp                                          ! geometric axis of torus

contains
    !> Collect the tests.
    subroutine collect_calc_int_vol(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("torus_analytic", test_torus_analytic), &
            &new_unittest("torus_convergence", test_torus_convergence), &
            &new_unittest("decoupled_dim_2", test_decoupled_dim_2) &
            &]
    end subroutine collect_calc_int_vol

    !--------------------------------------------------------------------------
    ! helpers
    !--------------------------------------------------------------------------

    !> Numerical integral of f = 1 - r^2 + i cos(theta) with J = r (R_0 +
    !! r cos(theta)) on a (theta, zeta, r) grid of the given dimensions.
    integer function torus_int(dims,f_int) result(ierr)
        integer, intent(in) :: dims(3)                                          ! grid size
        complex(dp), intent(out) :: f_int(1)                                    ! integral

        integer :: id, kd                                                       ! counters
        real(dp) :: r                                                           ! minor radius
        real(dp), allocatable :: theta(:,:,:), zeta(:,:,:)                      ! angles
        real(dp), allocatable :: norm(:)                                        ! normal variable
        real(dp), allocatable :: J(:,:,:)                                       ! Jacobian
        complex(dp), allocatable :: fun(:,:,:,:)                                ! integrand

        allocate(theta(dims(1),dims(2),dims(3)))
        allocate(zeta(dims(1),dims(2),dims(3)))
        allocate(J(dims(1),dims(2),dims(3)))
        allocate(fun(dims(1),dims(2),dims(3),1))
        allocate(norm(dims(3)))

        ierr = calc_eqd_grid(theta,0._dp,2._dp*pi,1)
        if (ierr.ne.0) return
        if (dims(2).gt.1) then
            ierr = calc_eqd_grid(zeta,0._dp,2._dp*pi,2)
            if (ierr.ne.0) return
        else
            zeta = 0._dp
        end if
        ierr = calc_eqd_grid(norm,0._dp,1._dp)
        if (ierr.ne.0) return

        do kd = 1,dims(3)
            r = norm(kd)
            fun(:,:,kd,1) = 1._dp - r**2 + iu*cos(theta(:,:,kd))
            do id = 1,dims(1)
                J(id,:,kd) = r*(R_0+r*cos(theta(id,:,kd)))
            end do
        end do

        f_int = 0._dp
        ierr = calc_int_vol(theta,zeta,norm,J,fun,f_int)
    end function torus_int

    !--------------------------------------------------------------------------
    ! tests
    !--------------------------------------------------------------------------

    !> Agreement with the analytical value at a moderate resolution.
    subroutine test_torus_analytic(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp), parameter :: tol = 5.e-3_dp                                   ! relative tolerance

        integer :: ierr                                                         ! error status
        complex(dp) :: f_int(1)                                                 ! numerical integral
        complex(dp) :: f_ana                                                    ! analytical integral
        real(dp) :: err                                                         ! relative error

        ierr = torus_int([33,34,21],f_int)
        call check(error, ierr, 0, 'calc_int_vol failed')
        if (allocated(error)) return

        f_ana = R_0*pi**2 + iu*2._dp*pi**2/3._dp
        err = abs(f_int(1)-f_ana)/abs(f_ana)
        write(error_unit,'(A,2ES12.4,A,2ES12.4,A,ES10.3)') &
            &'   [int_vol] numerical = ',f_int(1),', analytical = ',f_ana,&
            &', rel err = ',err
        call check(error, err.lt.tol, &
            &'torus integral deviates from analytical: '//trim(r2str(err)))
        if (allocated(error)) return
    end subroutine test_torus_analytic

    !> Doubling the resolution has to reduce the error by ~4 (second order);
    !! require at least a factor 3.
    subroutine test_torus_convergence(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ierr                                                         ! error status
        complex(dp) :: f_int(1)                                                 ! numerical integral
        complex(dp) :: f_ana                                                    ! analytical integral
        real(dp) :: err_lo, err_hi                                              ! errors at the two resolutions

        f_ana = R_0*pi**2 + iu*2._dp*pi**2/3._dp

        ierr = torus_int([17,18,11],f_int)
        call check(error, ierr, 0, 'calc_int_vol failed')
        if (allocated(error)) return
        err_lo = abs(f_int(1)-f_ana)/abs(f_ana)

        ierr = torus_int([33,34,21],f_int)
        call check(error, ierr, 0, 'calc_int_vol failed')
        if (allocated(error)) return
        err_hi = abs(f_int(1)-f_ana)/abs(f_ana)

        write(error_unit,'(A,ES10.3,A,ES10.3)') &
            &'   [int_vol convergence] err coarse = ',err_lo,&
            &', err fine = ',err_hi
        call check(error, err_hi.lt.err_lo/3._dp, &
            &'no second-order convergence: '//trim(r2str(err_lo))//' -> '//&
            &trim(r2str(err_hi)))
        if (allocated(error)) return
    end subroutine test_torus_convergence

    !> A single point in dimension 2: the implementation substitutes the
    !! full turn 2 pi for the missing angular extent, so the axisymmetric
    !! integrand (zeta does not appear in f or J) has to give the complete
    !! toroidal integral again.
    subroutine test_decoupled_dim_2(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp), parameter :: tol = 5.e-3_dp                                   ! relative tolerance

        integer :: ierr                                                         ! error status
        complex(dp) :: f_int(1)                                                 ! numerical integral
        complex(dp) :: f_ana                                                    ! analytical integral
        real(dp) :: err                                                         ! relative error

        ierr = torus_int([65,1,41],f_int)
        call check(error, ierr, 0, 'calc_int_vol failed')
        if (allocated(error)) return

        f_ana = R_0*pi**2 + iu*2._dp*pi**2/3._dp
        err = abs(f_int(1)-f_ana)/abs(f_ana)
        write(error_unit,'(A,2ES12.4,A,2ES12.4,A,ES10.3)') &
            &'   [int_vol dim-2 = 1] numerical = ',f_int(1),&
            &', analytical = ',f_ana,', rel err = ',err
        call check(error, err.lt.tol, &
            &'singleton-dimension integral deviates from analytical: '//&
            &trim(r2str(err)))
        if (allocated(error)) return
    end subroutine test_decoupled_dim_2
end module test_calc_int_vol
