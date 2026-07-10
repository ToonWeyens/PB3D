!------------------------------------------------------------------------------!
!> Full-stack tests of the assembled axisymmetric vacuum matrices G and H
!! (vac_ops.calc_GH, style 2) through Green's identities, and of the boundary
!! potential solve (vac_ops.solve_Phi_BEM).
!!
!! For potentials phi e^{i n zeta}, with dphi := - norm . grad phi, the
!! boundary element discretization satisfies the jump relations
!!    H phi            = G dphi   for phi harmonic inside the boundary,
!!    (H + 4 pi) phi   = G dphi   for phi harmonic outside and decaying,
!! the first of which is the check that the ldebug branch of calc_GH plots
!! for manual inspection. Here both are automated with assertions on a
!! synthetic circular boundary
!!    R(t) = R_0 + a cos t,   Z(t) = a sin t,  t in [0,2 pi],
!! built the same way store_vac_HEL builds its boundary (the last point
!! duplicates the first).
!!
!! Test potentials:
!!  - phi = R^n        (regular inside the boundary),
!!  - phi = R^n Z      (regular inside the boundary),
!!  - phi = Q_{n-1/2}(gamma(x,x_src))/sqrt(R R_src) with x_src inside the
!!    plasma: the n-mode Green's function, regular and decaying outside --
!!    the class the vacuum response solve works on.
!!
!! The residuals are normalized by max |H phi| (or max |4 pi phi|).
!------------------------------------------------------------------------------!
module test_vac_greens
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str, i2str
    use dtorh, only: dtorh1
    use vac_vars, only: vac_type
    use vac_ops, only: calc_GH, solve_Phi_BEM

    implicit none
    private
    public collect_vac_greens

    ! circular boundary parameters
    real(dp), parameter :: R_0 = 3.0_dp                                         ! major radius
    real(dp), parameter :: a_min = 1.0_dp                                       ! minor radius
    integer, parameter :: prim_X_test = 2                                       ! toroidal mode number

contains
    !> Collect the tests.
    subroutine collect_vac_greens(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("greens_identity", test_greens_identity), &
            &new_unittest("greens_identity_convergence", &
            &   test_greens_convergence), &
            &new_unittest("solve_Phi_BEM_roundtrip", test_solve_roundtrip), &
            &new_unittest("vac_response_cylinder", test_response_cylinder) &
            &]
    end subroutine collect_vac_greens

    !--------------------------------------------------------------------------
    ! helpers
    !--------------------------------------------------------------------------

    !> Set up the synthetic circular vacuum boundary (single process).
    subroutine setup_circle_vac(vac,n_bnd,ierr)
        type(vac_type), intent(inout) :: vac                                    ! vacuum variables
        integer, intent(in) :: n_bnd                                            ! number of boundary points (last = first)
        integer, intent(out) :: ierr                                            ! error status

        integer :: id                                                           ! counter
        real(dp) :: t, R                                                        ! angle, major radius

        ierr = vac%init(2,n_bnd,prim_X_test,[n_bnd,1],1._dp)
        if (ierr.ne.0) return

        do id = 1,n_bnd
            t = 2._dp*pi*(id-1)/(n_bnd-1)
            R = R_0 + a_min*cos(t)
            vac%ang(id,1) = t
            vac%x_vec(id,:) = [R, a_min*sin(t)]
            vac%norm(id,:) = [-R*a_min*cos(t), -R*a_min*sin(t)]
            vac%dnorm(id,:) = [a_min*sin(t)*(a_min*cos(t)+R), &
                &a_min*(a_min*sin(t)**2-R*cos(t))]
        end do
    end subroutine setup_circle_vac

    !> Test potential and its dphi = -norm.grad phi on the boundary.
    !!
    !! Kinds: 1: R^n, 2: R^n Z, 3: n-mode Green's function with an interior
    !! source (exterior harmonic).
    subroutine eval_potential(vac,kind,phi,dphi,ierr)
        type(vac_type), intent(in) :: vac                                       ! vacuum variables
        integer, intent(in) :: kind                                             ! which potential
        real(dp), intent(out) :: phi(:), dphi(:)                                ! potential and -norm.grad
        integer, intent(out) :: ierr                                            ! error status

        integer :: id                                                           ! counter
        integer :: n                                                            ! mode number
        real(dp) :: R, Z                                                        ! coordinates
        real(dp) :: x_pert(2)                                                   ! perturbed point
        real(dp) :: fp, fm                                                      ! function values
        real(dp), parameter :: eps = 1.e-7_dp                                   ! step for central difference
        real(dp), parameter :: x_src(2) = [R_0-0.3_dp*a_min,0.1_dp*a_min]       ! interior source for kind 3

        ierr = 0
        n = vac%prim_X
        do id = 1,vac%n_bnd
            R = vac%x_vec(id,1)
            Z = vac%x_vec(id,2)
            select case (kind)
                case (1)                                                        ! R^n
                    phi(id) = R**n
                    dphi(id) = -vac%norm(id,1)*n*R**(n-1)
                case (2)                                                        ! R^n Z
                    phi(id) = R**n*Z
                    dphi(id) = -vac%norm(id,1)*n*R**(n-1)*Z - &
                        &vac%norm(id,2)*R**n
                case (3)                                                        ! Green's function source
                    ierr = green_pot(vac%x_vec(id,:),x_src,n,phi(id))
                    if (ierr.ne.0) return
                    ! dphi by central difference along norm
                    x_pert = vac%x_vec(id,:) + eps*vac%norm(id,:)
                    ierr = green_pot(x_pert,x_src,n,fp)
                    if (ierr.ne.0) return
                    x_pert = vac%x_vec(id,:) - eps*vac%norm(id,:)
                    ierr = green_pot(x_pert,x_src,n,fm)
                    if (ierr.ne.0) return
                    dphi(id) = -(fp-fm)/(2._dp*eps)
            end select
        end do
    end subroutine eval_potential

    !> The n-mode Green's function potential
    !! Q_{n-1/2}(gamma(x,x_src))/sqrt(R R_src).
    integer function green_pot(x,x_src,n,res) result(ierr)
        real(dp), intent(in) :: x(2), x_src(2)                                  ! evaluation and source points
        integer, intent(in) :: n                                                ! mode number
        real(dp), intent(out) :: res                                            ! potential

        real(dp) :: gam                                                         ! argument
        real(dp) :: pl_loc(0:n), ql_loc(0:n)                                    ! toroidal harmonics
        integer :: newn                                                         ! maximum reached degree

        gam = 1._dp + sum((x-x_src)**2)/(2._dp*x(1)*x_src(1))
        ierr = dtorh1(gam,0,n,pl_loc,ql_loc,newn)
        if (ierr.eq.0 .and. newn.lt.n) ierr = 1
        res = ql_loc(n)/sqrt(x(1)*x_src(1))
    end function green_pot

    !> Assemble G and H for a circular boundary and return the Green's
    !! identity residuals max|H phi - G dphi| / max|H phi| for the three test
    !! potentials.
    subroutine greens_residuals(n_bnd,res,ierr)
        integer, intent(in) :: n_bnd                                            ! number of boundary points
        real(dp), intent(out) :: res(3)                                         ! residuals for the potentials
        integer, intent(out) :: ierr                                            ! error status

        type(vac_type) :: vac                                                   ! vacuum variables
        integer :: kind                                                         ! potential kind
        real(dp), allocatable :: phi(:), dphi(:)                                ! potential and normal derivative
        real(dp), allocatable :: lhs(:), rhs(:)                                 ! H phi and G dphi

        call setup_circle_vac(vac,n_bnd,ierr)
        if (ierr.ne.0) return

        ierr = calc_GH(vac)
        if (ierr.ne.0) return

        allocate(phi(n_bnd),dphi(n_bnd),lhs(n_bnd),rhs(n_bnd))
        do kind = 1,3
            call eval_potential(vac,kind,phi,dphi,ierr)
            if (ierr.ne.0) return
            lhs = matmul(vac%H,phi)                                             ! single process: local = global
            rhs = matmul(vac%G,dphi)
            select case (kind)
                case (1,2)                                                      ! interior: H phi = G dphi
                    res(kind) = maxval(abs(lhs-rhs))/maxval(abs(lhs))
                case (3)                                                        ! exterior: (H + 4 pi) phi = G dphi
                    res(kind) = maxval(abs(lhs+4._dp*pi*phi-rhs))/&
                        &maxval(abs(4._dp*pi*phi))
            end select
            write(error_unit,'(A,I6,A,I2,A,ES10.3)') '   [greens residual] &
                &n_bnd = ',n_bnd,', potential ',kind,': ',res(kind)
        end do

        call vac%dealloc()
    end subroutine greens_residuals

    !--------------------------------------------------------------------------
    ! tests
    !--------------------------------------------------------------------------

    !> Green's identities for the three test potentials.
    !!
    !! The discretization converges at first order (the trapezoidal rule
    !! meets the logarithmic kernel singularity right outside the tiny
    !! analytically-integrated region), so the residuals at n_bnd = 101 are
    !! at the percent level; convergence is checked separately.
    subroutine test_greens_identity(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_bnd = 101                                       ! boundary points
        real(dp), parameter :: tol = 5.e-2_dp                                   ! relative tolerance

        integer :: ierr, kind                                                   ! error status, counter
        real(dp) :: res(3)                                                      ! residuals

        call greens_residuals(n_bnd,res,ierr)
        call check(error, ierr, 0, 'assembly failed')
        if (allocated(error)) return

        do kind = 1,3
            call check(error, res(kind).lt.tol, &
                &'Green identity residual for potential '//&
                &trim(i2str(kind))//' too large: '//trim(r2str(res(kind))))
            if (allocated(error)) return
        end do
    end subroutine test_greens_identity

    !> The identity residual has to decrease with resolution (first order:
    !! it halves when the number of points doubles).
    subroutine test_greens_convergence(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: ierr, kind                                                   ! error status, counter
        real(dp) :: res_lo(3), res_hi(3)                                        ! residuals at two resolutions

        call greens_residuals(51,res_lo,ierr)
        call check(error, ierr, 0, 'assembly failed')
        if (allocated(error)) return
        call greens_residuals(101,res_hi,ierr)
        call check(error, ierr, 0, 'assembly failed')
        if (allocated(error)) return

        do kind = 1,3
            call check(error, res_hi(kind).lt.0.7_dp*res_lo(kind), &
                &'no convergence for potential '//trim(i2str(kind))//': '//&
                &trim(r2str(res_lo(kind)))//' -> '//trim(r2str(res_hi(kind))))
            if (allocated(error)) return
        end do
    end subroutine test_greens_convergence

    !> Solving (H + 4 pi) Phi = G R with R = dphi of a known vacuum-side
    !! harmonic potential has to return its boundary trace Phi = phi.
    !!
    !! Note that this also covers the duplicated seam point of the closed
    !! boundary: H alone is exactly singular there (two identical rows), but
    !! the 4 pi diagonal term of the exterior operator lifts the degeneracy
    !! and enforces the equality of the two coinciding solution values.
    subroutine test_solve_roundtrip(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_bnd = 101                                       ! boundary points
        real(dp), parameter :: tol = 5.e-2_dp                                   ! relative tolerance

        type(vac_type) :: vac                                                   ! vacuum variables
        integer :: ierr                                                         ! error status
        integer :: desc_RPhi(9)                                                 ! descriptor for R and Phi
        real(dp), allocatable :: phi(:), dphi(:)                                ! potential and normal derivative
        real(dp), allocatable :: R_mat(:,:), Phi_mat(:,:)                       ! right-hand side and solution
        real(dp) :: err_sol                                                     ! solution error

        call setup_circle_vac(vac,n_bnd,ierr)
        call check(error, ierr, 0, 'setup failed')
        if (allocated(error)) return

        ierr = calc_GH(vac)
        call check(error, ierr, 0, 'assembly failed')
        if (allocated(error)) return

        allocate(phi(n_bnd),dphi(n_bnd))
        call eval_potential(vac,3,phi,dphi,ierr)                                ! exterior (Green's function) potential
        call check(error, ierr, 0, 'potential evaluation failed')
        if (allocated(error)) return

        allocate(R_mat(n_bnd,1),Phi_mat(n_bnd,1))
        R_mat(:,1) = dphi
        Phi_mat = 0._dp
        call descinit(desc_RPhi,n_bnd,1,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
            &max(1,vac%n_loc(1)),ierr)
        call check(error, ierr, 0, 'descinit failed')
        if (allocated(error)) return

        ierr = solve_Phi_BEM(vac,R_mat,Phi_mat,[n_bnd,1],[vac%n_loc(1),1],&
            &reshape([1,1],[2,1]),desc_RPhi)
        call check(error, ierr, 0, 'solve_Phi_BEM failed')
        if (allocated(error)) return

        err_sol = maxval(abs(Phi_mat(:,1)-phi))/maxval(abs(phi))
        call check(error, err_sol.lt.tol, &
            &'solve round-trip error too large: '//trim(r2str(err_sol)))
        if (allocated(error)) return

        call vac%dealloc()
    end subroutine test_solve_roundtrip

    !> The vacuum response matrix on the circular boundary against the
    !! analytical large-aspect-ratio (cylinder) limit.
    !!
    !! In the cylinder limit, the decaying exterior solution for poloidal
    !! Fourier boundary data e^{i m theta} is proportional to r^-|m|, so the
    !! full chain of calc_vac_res (Neumann data (n q - m) e^{-i m theta},
    !! exterior solve, projection with the integration rule) reduces to
    !!    res_{m m'} = -2 pi (n q - m)^2 / (R_0 |m| mu_0) delta_{m m'},
    !! with toroidal corrections of order a/R_0. This pins sign, scaling and
    !! mode structure of the response that enters the SLEPc boundary
    !! condition (set_BC_4).
    subroutine test_response_cylinder(error)
        use X_vars, only: n_mod_X, modes_type
        use num_vars, only: use_pol_flux_F, eq_style
        use eq_vars, only: vac_perm
        use vac_ops, only: calc_vac_res

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_bnd = 201                                       ! boundary points
        integer, parameter :: n_mod = 5                                         ! number of poloidal modes, m = 1..n_mod
        real(dp), parameter :: R_big = 20.0_dp                                  ! major radius (large aspect ratio)
        real(dp), parameter :: jq = 1.35_dp                                     ! safety factor at edge (nonresonant: n jq - m /= 0)
        real(dp), parameter :: tol_diag = 0.25_dp                               ! relative tolerance on the diagonal (aspect-ratio corrections)
        real(dp), parameter :: tol_offdiag = 0.10_dp                            ! off-diagonal tolerance, relative to largest diagonal

        type(vac_type) :: vac                                                   ! vacuum variables
        type(modes_type) :: mds                                                 ! minimal modes variables
        integer :: ierr, id, jd                                                 ! error status, counters
        integer :: n_mod_X_old                                                  ! original n_mod_X
        real(dp) :: t, R                                                        ! angle, major radius
        real(dp) :: res_ana                                                     ! analytical response
        real(dp) :: rel_diff                                                    ! relative difference

        ! the response needs several modes; restore module state afterwards
        n_mod_X_old = n_mod_X
        n_mod_X = n_mod
        use_pol_flux_F = .true.
        eq_style = 2                                                            ! HELENA (axisymmetric)

        ! minimal modes tables: only the last row of m is used
        allocate(mds%m(1,n_mod),mds%n(1,n_mod))
        mds%m(1,:) = [(id, id=1,n_mod)]
        mds%n(1,:) = prim_X_test

        ! circular boundary at large aspect ratio
        ierr = vac%init(2,n_bnd,prim_X_test,[n_bnd,1],jq)
        call check(error, ierr, 0, 'init failed')
        if (allocated(error)) return
        do id = 1,n_bnd
            t = 2._dp*pi*(id-1)/(n_bnd-1)
            R = R_big + a_min*cos(t)
            vac%ang(id,1) = t
            vac%x_vec(id,:) = [R, a_min*sin(t)]
            vac%norm(id,:) = [-R*a_min*cos(t), -R*a_min*sin(t)]
            vac%dnorm(id,:) = [a_min*sin(t)*(a_min*cos(t)+R), &
                &a_min*(a_min*sin(t)**2-R*cos(t))]
        end do

        ierr = calc_vac_res(mds,vac)
        call check(error, ierr, 0, 'calc_vac_res failed')
        if (allocated(error)) return

        do id = 1,n_mod                                                         ! single process: last rank = rank 0 has vac%res
            do jd = 1,n_mod
                if (id.eq.jd) then
                    res_ana = -2._dp*pi*(prim_X_test*jq-mds%m(1,id))**2/&
                        &(R_big*abs(mds%m(1,id))*vac_perm)
                    rel_diff = abs(real(vac%res(id,id))-res_ana)/&
                        &abs(res_ana)
                    write(error_unit,'(A,I2,A,ES12.5,A,ES12.5,A,F6.3)') &
                        &'   [vac response] m = ',mds%m(1,id),': res = ',&
                        &real(vac%res(id,id)),', cylinder = ',res_ana,&
                        &', rel diff = ',rel_diff
                    call check(error, rel_diff.lt.tol_diag, &
                        &'diagonal response for m = '//trim(i2str(id))//&
                        &' deviates from cylinder limit by '//&
                        &trim(r2str(rel_diff)))
                else
                    call check(error, abs(vac%res(id,jd)).lt.tol_offdiag*&
                        &2._dp*pi*(prim_X_test*jq-1._dp)**2/&
                        &(R_big*vac_perm), &
                        &'off-diagonal response ('//trim(i2str(id))//','//&
                        &trim(i2str(jd))//') too large: '//&
                        &trim(r2str(abs(vac%res(id,jd)))))
                end if
                if (allocated(error)) exit
            end do
            if (allocated(error)) exit
        end do

        ! clean up and restore module state
        call vac%dealloc()
        call mds%dealloc()
        n_mod_X = n_mod_X_old
        if (allocated(error)) return
    end subroutine test_response_cylinder
end module test_vac_greens
