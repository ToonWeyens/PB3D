!------------------------------------------------------------------------------!
!> Unit tests of the pure numerical utilities in num_utilities, enabled by the
!! split of the PSPLINE-dependent spline wrappers into spline_utilities.
!!
!! Reference values are analytical or hand-computed; tolerances reflect the
!! order of the methods, not platform noise.
!------------------------------------------------------------------------------!
module test_num_utilities
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str, i2str
    use num_utilities, only: GCD, LCM, fac, calc_int, calc_coeff_fin_diff, &
        &solve_vand, c, is_sym, bubble_sort, calc_det, calc_inv, calc_ext_var

    implicit none
    private
    public collect_num_utilities

contains
    !> Collect the tests.
    subroutine collect_num_utilities(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("integer_helpers", test_integer_helpers), &
            &new_unittest("calc_int", test_calc_int), &
            &new_unittest("fin_diff_coeffs", test_fin_diff_coeffs), &
            &new_unittest("solve_vand", test_solve_vand), &
            &new_unittest("sym_indexing", test_sym_indexing), &
            &new_unittest("bubble_sort", test_bubble_sort), &
            &new_unittest("det_inv_0D", test_det_inv_0D), &
            &new_unittest("ext_var", test_ext_var) &
            &]
    end subroutine collect_num_utilities

    !> GCD, LCM and factorial.
    subroutine test_integer_helpers(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, GCD(12,18), 6, 'GCD(12,18)')
        if (allocated(error)) return
        call check(error, GCD(17,5), 1, 'GCD(17,5)')
        if (allocated(error)) return
        call check(error, LCM(4,6), 12, 'LCM(4,6)')
        if (allocated(error)) return
        call check(error, LCM(7,13), 91, 'LCM(7,13)')
        if (allocated(error)) return
        call check(error, fac(0), 1, 'fac(0)')
        if (allocated(error)) return
        call check(error, fac(5), 120, 'fac(5)')
        if (allocated(error)) return
    end subroutine test_integer_helpers

    !> Trapezoidal integration, equidistant and general grids, against
    !! int_0^pi sin = 2 (second-order rule: tolerance ~ h^2).
    subroutine test_calc_int(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 201                                           ! grid points
        integer :: id, ierr                                                     ! counter, error status
        real(dp) :: x(n), f(n), f_int(n)                                        ! abscissa, integrand, integral
        real(dp) :: step                                                        ! step size

        step = pi/(n-1)
        do id = 1,n
            x(id) = (id-1)*step
            f(id) = sin(x(id))
        end do

        ! equidistant version
        ierr = calc_int(f,step,f_int)
        call check(error, ierr, 0, 'calc_int (eqd) failed')
        if (allocated(error)) return
        call check(error, f_int(n), 2._dp, thr=2._dp*step**2, &
            &message='eqd integral of sin over [0,pi]: '//trim(r2str(f_int(n))))
        if (allocated(error)) return

        ! general-grid version on the same grid
        ierr = calc_int(f,x,f_int)
        call check(error, ierr, 0, 'calc_int (reg) failed')
        if (allocated(error)) return
        call check(error, f_int(n), 2._dp, thr=2._dp*step**2, &
            &message='reg integral of sin over [0,pi]: '//trim(r2str(f_int(n))))
        if (allocated(error)) return
    end subroutine test_calc_int

    !> Finite-difference weights against the classical values (for unit step
    !! size; the callers scale by the actual step).
    subroutine test_fin_diff_coeffs(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ierr, id                                                     ! error status, counter
        real(dp), allocatable :: coeff(:)                                       ! coefficients
        real(dp), parameter :: ref_c1(3) = [-0.5_dp,0._dp,0.5_dp]               ! central first derivative
        real(dp), parameter :: ref_c2(3) = [1._dp,-2._dp,1._dp]                 ! central second derivative
        real(dp), parameter :: ref_l1(2) = [-1._dp,1._dp]                       ! left first derivative

        ierr = calc_coeff_fin_diff(1,3,2,coeff)                                 ! d/dx, 3 points, centered
        call check(error, ierr, 0, 'calc_coeff_fin_diff failed')
        if (allocated(error)) return
        do id = 1,3
            call check(error, coeff(id), ref_c1(id), thr=1.e-12_dp, &
                &message='central 1st-derivative weight '//trim(i2str(id)))
            if (allocated(error)) return
        end do

        ierr = calc_coeff_fin_diff(2,3,2,coeff)                                 ! d2/dx2, 3 points, centered
        call check(error, ierr, 0, 'calc_coeff_fin_diff failed')
        if (allocated(error)) return
        do id = 1,3
            call check(error, coeff(id), ref_c2(id), thr=1.e-12_dp, &
                &message='central 2nd-derivative weight '//trim(i2str(id)))
            if (allocated(error)) return
        end do

        ierr = calc_coeff_fin_diff(1,2,2,coeff)                                 ! d/dx, 2 points, at right point
        call check(error, ierr, 0, 'calc_coeff_fin_diff failed')
        if (allocated(error)) return
        do id = 1,2
            call check(error, coeff(id), ref_l1(id), thr=1.e-12_dp, &
                &message='left 1st-derivative weight '//trim(i2str(id)))
            if (allocated(error)) return
        end do
    end subroutine test_fin_diff_coeffs

    !> Björck-Pereyra Vandermonde solver against a directly-constructed
    !! system: sum_j a_i^(j-1) x_j = b_i.
    subroutine test_solve_vand(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4                                             ! system size
        real(dp), parameter :: a(n) = [0.5_dp,1.0_dp,2.0_dp,3.0_dp]             ! nodes
        real(dp), parameter :: x_ref(n) = [1._dp,-2._dp,0.5_dp,3._dp]           ! chosen solution
        integer :: id, jd                                                       ! counters
        real(dp) :: b(n), x(n)                                                  ! right-hand side, solution

        ! construct b_i = sum_j a_i^(j-1) x_j (rows are powers of a_i, the
        ! layout shown in the solve_vand documentation)
        b = 0._dp
        do id = 1,n
            do jd = 1,n
                b(id) = b(id) + a(id)**(jd-1)*x_ref(jd)
            end do
        end do
        call solve_vand(n,a,b,x)
        do id = 1,n
            call check(error, x(id), x_ref(id), thr=1.e-10_dp, rel=.true., &
                &message='Vandermonde solution component '//trim(i2str(id)))
            if (allocated(error)) return
        end do
    end subroutine test_solve_vand

    !> The symmetric-storage index helper c and is_sym: for a symmetric 2x2
    !! matrix the storage is [(1,1),(2,1),(2,2)] and (1,2) maps onto (2,1).
    subroutine test_sym_indexing(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ierr                                                         ! error status
        logical :: sym                                                          ! whether symmetric

        call check(error, c([1,1],.true.,2), 1, 'c(1,1) sym')
        if (allocated(error)) return
        call check(error, c([2,1],.true.,2), 2, 'c(2,1) sym')
        if (allocated(error)) return
        call check(error, c([1,2],.true.,2), 2, 'c(1,2) sym = c(2,1)')
        if (allocated(error)) return
        call check(error, c([2,2],.true.,2), 3, 'c(2,2) sym')
        if (allocated(error)) return
        call check(error, c([2,1],.false.,2), 2, 'c(2,1) full (column-major)')
        if (allocated(error)) return
        call check(error, c([1,2],.false.,2), 3, 'c(1,2) full (column-major)')
        if (allocated(error)) return

        ierr = is_sym(2,3,sym)                                                  ! 2x2 with 3 elements: symmetric
        call check(error, ierr, 0, 'is_sym failed')
        if (allocated(error)) return
        call check(error, sym, 'is_sym(2,3)')
        if (allocated(error)) return
        ierr = is_sym(2,4,sym)                                                  ! 2x2 with 4 elements: full
        call check(error, ierr, 0, 'is_sym failed')
        if (allocated(error)) return
        call check(error, .not.sym, 'is_sym(2,4)')
        if (allocated(error)) return
    end subroutine test_sym_indexing

    !> Sorting with pivot tracking.
    subroutine test_bubble_sort(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: id                                                           ! counter
        real(dp) :: a(5)                                                        ! array to sort
        integer :: piv(5)                                                       ! pivots

        a = [3._dp,1._dp,-2._dp,5._dp,0._dp]
        call bubble_sort(a,piv)
        do id = 2,5
            call check(error, a(id).ge.a(id-1), 'array not sorted at '//&
                &trim(i2str(id)))
            if (allocated(error)) return
        end do
        ! pivots reconstruct the sorted array from the original
        call check(error, all(piv.eq.[3,5,2,1,4]), 'pivots incorrect')
        if (allocated(error)) return
    end subroutine test_bubble_sort

    !> LAPACK-backed 0-D determinant and inverse on a known 3x3 matrix.
    subroutine test_det_inv_0D(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ierr, id, jd                                                 ! error status, counters
        real(dp) :: A(3,3), A_inv(3,3), A_id(3,3)                               ! matrix, inverse, product
        real(dp) :: det_A                                                       ! determinant

        A = reshape([2._dp,0._dp,1._dp,&
                    &1._dp,3._dp,0._dp,&
                    &0._dp,1._dp,4._dp],[3,3])
        ierr = calc_det(det_A,A)
        call check(error, ierr, 0, 'calc_det failed')
        if (allocated(error)) return
        call check(error, det_A, 25._dp, thr=1.e-12_dp, rel=.true., &
            &message='det: '//trim(r2str(det_A)))
        if (allocated(error)) return

        A_inv = A
        ierr = calc_inv(A_inv,A)
        call check(error, ierr, 0, 'calc_inv failed')
        if (allocated(error)) return
        A_id = matmul(A,A_inv)
        do id = 1,3
            do jd = 1,3
                call check(error, A_id(id,jd), &
                    &merge(1._dp,0._dp,id.eq.jd), thr=1.e-12_dp, &
                    &message='A A^-1 not identity at ('//trim(i2str(id))//&
                    &','//trim(i2str(jd))//')')
                if (allocated(error)) return
            end do
        end do
    end subroutine test_det_inv_0D

    !> Polynomial extrapolation: exact for polynomials up to the number of
    !! points minus one, including derivatives.
    subroutine test_ext_var(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ierr, id                                                     ! error status, counter
        real(dp) :: x(4), y(4)                                                  ! sample points of y = x^3 - 2 x
        real(dp) :: y_ext                                                       ! extrapolated value

        do id = 1,4
            x(id) = real(id-1,dp)*0.5_dp
            y(id) = x(id)**3 - 2._dp*x(id)
        end do
        ierr = calc_ext_var(y_ext,y,x,2.0_dp)                                   ! value at x = 2
        call check(error, ierr, 0, 'calc_ext_var failed')
        if (allocated(error)) return
        call check(error, y_ext, 4._dp, thr=1.e-10_dp, rel=.true., &
            &message='cubic extrapolation: '//trim(r2str(y_ext)))
        if (allocated(error)) return
        ierr = calc_ext_var(y_ext,y,x,2.0_dp,1)                                 ! derivative at x = 2
        call check(error, ierr, 0, 'calc_ext_var failed')
        if (allocated(error)) return
        call check(error, y_ext, 10._dp, thr=1.e-9_dp, rel=.true., &
            &message='cubic derivative extrapolation: '//trim(r2str(y_ext)))
        if (allocated(error)) return
    end subroutine test_ext_var
end module test_num_utilities
