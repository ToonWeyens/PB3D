!------------------------------------------------------------------------------!
!> Full-stack tests of the 3-D (field-line aligned, style 1) vacuum: the
!! singular interval kernel, the assembled G and H matrices through Green's
!! identities, and the vacuum response against the (independently verified)
!! axisymmetric style-2 calculation of the same boundary.
!!
!! The geometry is an analytical circular torus
!!    x = (R cos zeta, R sin zeta, a sin theta),   R = R_0 + a cos theta,
!! covered by field lines zeta = alpha + q theta on a product grid
!! (theta, alpha) in [0,2 pi] x [0,2 pi] (both ends included; the map is a
!! shear, so this tiles the torus once, with the standard trapezoidal
!! half-weights on the duplicated edges). In PB3D's flux coordinates
!! (alpha, psi, theta) the exact geometric quantities are
!!    e_theta = dx/dtheta|_alpha,   e_alpha = dx/dalpha,
!!    norm    = J nabla psi = e_theta x e_alpha,
!!    h_fac   = (g_aa, g_at, g_tt, ...) = (|e_alpha|^2, e_alpha.e_theta,
!!              |e_theta|^2, ...),
!! all available in closed form.
!!
!! For the assembled matrices, with dphi := norm . grad phi, the jump
!! relations analogous to the style-2 ones are measured empirically; note
!! that calc_GH_1 constructs the H diagonal through the constant-potential
!! row-sum identity (H 1 = -4 pi 1 by construction).
!------------------------------------------------------------------------------!
module test_vac_3d
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str, i2str
    use vac_vars, only: vac_type
    use vac_utilities, only: calc_GH_int_1
    use vac_ops, only: calc_GH

    implicit none
    private
    public collect_vac_3d

    ! circular torus parameters
    real(dp), parameter :: R_0 = 3.0_dp                                         ! major radius
    real(dp), parameter :: a_min = 1.0_dp                                       ! minor radius
    real(dp), parameter :: q_saf = 1.35_dp                                      ! safety factor at the boundary

contains
    !> Collect the tests.
    subroutine collect_vac_3d(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("GH_int_1_singular", test_int_1_singular), &
            &new_unittest("greens_identity_3d", test_greens_identity_3d), &
            &new_unittest("response_1_vs_2", test_response_1_vs_2) &
            &]
    end subroutine collect_vac_3d

    !--------------------------------------------------------------------------
    ! helpers
    !--------------------------------------------------------------------------

    !> Position and flux-coordinate tangents on the field line.
    pure subroutine torus(theta,alpha,x,e_theta,e_alpha)
        real(dp), intent(in) :: theta, alpha                                    ! angles
        real(dp), intent(out) :: x(3)                                           ! position
        real(dp), intent(out), optional :: e_theta(3), e_alpha(3)               ! tangents

        real(dp) :: R, zeta                                                     ! major radius, toroidal angle

        R = R_0 + a_min*cos(theta)
        zeta = alpha + q_saf*theta
        x = [R*cos(zeta), R*sin(zeta), a_min*sin(theta)]
        if (present(e_theta)) e_theta = &
            &[-a_min*sin(theta)*cos(zeta) - R*sin(zeta)*q_saf, &
            & -a_min*sin(theta)*sin(zeta) + R*cos(zeta)*q_saf, &
            &  a_min*cos(theta)]
        if (present(e_alpha)) e_alpha = [-R*sin(zeta), R*cos(zeta), 0._dp]
    end subroutine torus

    !> Cross product.
    pure function cross(u,v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)

        w = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
    end function cross

    !> Set up the synthetic style-1 vacuum on the torus (the geometry arrays
    !! are global, so every process fills them completely; only G, H and res
    !! are distributed), and the grid module variables that calc_GH_1 reads.
    subroutine setup_torus_vac(vac,n_par,n_alpha_loc,ierr)
        use grid_vars, only: min_par_X, max_par_X, min_alpha, max_alpha, n_alpha
        use rich_vars, only: n_par_X

        type(vac_type), intent(inout) :: vac                                    ! vacuum variables
        integer, intent(in) :: n_par, n_alpha_loc                               ! points per field line, number of field lines
        integer, intent(out) :: ierr                                            ! error status

        integer :: id, jd, kd                                                   ! counters
        real(dp) :: theta, alpha                                                ! angles
        real(dp) :: x(3), e_t(3), e_a(3)                                        ! position and tangents

        ! grid module variables (in units of pi, as in the input file)
        min_par_X = 0._dp
        max_par_X = 2._dp
        min_alpha = 0._dp
        max_alpha = 2._dp
        n_par_X = n_par
        n_alpha = n_alpha_loc

        ierr = vac%init(1,n_par*n_alpha_loc,2,[n_par,n_alpha_loc],q_saf)
        if (ierr.ne.0) return

        do jd = 1,n_alpha_loc
            alpha = 2._dp*pi*(jd-1)/(n_alpha_loc-1)
            do id = 1,n_par
                theta = 2._dp*pi*(id-1)/(n_par-1)
                kd = id + (jd-1)*n_par
                call torus(theta,alpha,x,e_theta=e_t,e_alpha=e_a)
                vac%x_vec(kd,:) = x
                vac%norm(kd,:) = cross(e_t,e_a)                                 ! J nabla psi
                vac%h_fac(kd,1) = sum(e_a**2)                                   ! g_alpha,alpha
                vac%h_fac(kd,2) = sum(e_a*e_t)                                  ! g_alpha,theta
                vac%h_fac(kd,3) = sum(e_t**2)                                   ! g_theta,theta
                vac%h_fac(kd,4) = 0._dp                                         ! only enters the (discarded) singular H diagonal
            end do
        end do
    end subroutine setup_torus_vac

    !--------------------------------------------------------------------------
    ! tests
    !--------------------------------------------------------------------------

    !> The singular branch of the 3-D interval kernel against brute-force
    !! quadrature of -1/distance over the half-cell, with the metric-form
    !! distance. The inner (alpha) integral is done analytically with the
    !! elementary asinh-type primitive; the outer integral by panel-refined
    !! Gauss quadrature toward the singular corner.
    subroutine test_int_1_singular(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp), parameter :: h_fac(4) = [4.0_dp,0.9_dp,2.5_dp,7.7_dp]         ! g_aa, g_at, g_tt, H factor
        real(dp), parameter :: steps(2) = [0.02_dp,0.03_dp]                     ! dpar, dalpha
        real(dp), parameter :: gl_x(8) = [&
            &-0.9602898564975363_dp, -0.7966664774136267_dp, &
            &-0.5255324099163290_dp, -0.1834346424956498_dp, &
            & 0.1834346424956498_dp,  0.5255324099163290_dp, &
            & 0.7966664774136267_dp,  0.9602898564975363_dp]
        real(dp), parameter :: gl_w(8) = [&
            &0.1012285362903763_dp, 0.2223810344533745_dp, &
            &0.3137066458778873_dp, 0.3626837833783620_dp, &
            &0.3626837833783620_dp, 0.3137066458778873_dp, &
            &0.2223810344533745_dp, 0.1012285362903763_dp]

        integer :: id, jd                                                       ! counters
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(2,3), x_in(3), norm_s(2,3)                              ! geometry (only distances matter)
        real(dp) :: s_no, sq_a, sq_t                                            ! nonorthogonality and metric roots
        real(dp) :: I_ref                                                       ! brute-force integral
        real(dp) :: p_lo, p_hi, u_q, w_q                                        ! quadrature panels and nodes
        real(dp) :: V_hc                                                        ! alpha half-cell in metric units
        real(dp) :: G_ref                                                       ! reference value
        real(dp) :: r_far                                                       ! distance to far end

        ! a singular interval: influence point = left source point; the right
        ! point sits at par distance steps(1) in the field-line direction
        ! (only its separation from x_in matters for the regular part)
        x_in = [1.0_dp,2.0_dp,3.0_dp]
        x_s(1,:) = x_in
        x_s(2,:) = x_in + sqrt(h_fac(3))*steps(1)*[1._dp,0._dp,0._dp]
        norm_s(1,:) = [0.4_dp,0.5_dp,0.6_dp]
        norm_s(2,:) = [0.4_dp,0.5_dp,0.7_dp]

        call calc_GH_int_1(G,H,x_s,x_in,norm_s,h_fac,steps,1.e-30_dp)

        ! brute force: I = int_0^{U} du int_{-V}^{V} dv / sqrt(u^2+2suv+v^2),
        ! outer u integral with geometrically refined panels toward u = 0,
        ! inner integral analytical:
        !   int dv/sqrt(...) = ln(v + su + sqrt(u^2+2suv+v^2))
        sq_a = sqrt(h_fac(1))
        sq_t = sqrt(h_fac(3))
        s_no = h_fac(2)/(sq_a*sq_t)
        V_hc = sq_a*steps(2)/2._dp
        I_ref = 0._dp
        p_hi = sq_t*steps(1)/2._dp                                              ! U: half cell in par direction
        do id = 1,60
            p_lo = p_hi*0.6_dp
            if (id.eq.60) p_lo = 0._dp
            do jd = 1,8
                u_q = 0.5_dp*(p_lo+p_hi) + 0.5_dp*(p_hi-p_lo)*gl_x(jd)
                w_q = 0.5_dp*(p_hi-p_lo)*gl_w(jd)
                if (u_q.le.0._dp) cycle
                I_ref = I_ref + w_q*log(&
                    &(V_hc+s_no*u_q+sqrt(u_q**2+2*s_no*u_q*V_hc+V_hc**2))/&
                    &(-V_hc+s_no*u_q+sqrt(u_q**2-2*s_no*u_q*V_hc+V_hc**2)))
            end do
            p_hi = p_lo
            if (p_hi.le.0._dp) exit
        end do
        G_ref = -2._dp*I_ref/(steps(1)*steps(2)*sq_a*sq_t)

        call check(error, G(1), G_ref, thr=1.e-8_dp*abs(G_ref), &
            &message='singular G: '//trim(r2str(G(1)))//' vs brute force '//&
            &trim(r2str(G_ref)))
        if (allocated(error)) return

        ! the H value of the singular endpoint is h_fac(4) G by construction
        ! (it never enters the final system, see calc_GH_1)
        call check(error, H(1), h_fac(4)*G(1), thr=1.e-14_dp*abs(H(1)), &
            &message='singular H not h_fac(4) G')
        if (allocated(error)) return

        ! the far endpoint keeps the regular kernel
        r_far = sqrt(sum((x_s(2,:)-x_in)**2))
        call check(error, G(2), -1._dp/r_far, thr=1.e-14_dp, rel=.true., &
            &message='far-end G not regular kernel')
        if (allocated(error)) return
        call check(error, H(2), sum(norm_s(2,:)*(x_s(2,:)-x_in))/r_far**3, &
            &thr=1.e-14_dp, rel=.true., &
            &message='far-end H not regular dipole kernel')
        if (allocated(error)) return
    end subroutine test_int_1_singular

    !> Green's identities for the assembled 3-D matrices on the torus.
    !!
    !! With dphi := norm . grad phi, the offsets c in
    !!    H phi - G dphi = c phi
    !! are measured for an interior-harmonic and an exterior-harmonic
    !! potential. By the constant-potential construction of the H diagonal,
    !! c = -4 pi exactly for phi = 1; the tests establish the offsets and
    !! the discretization error for nontrivial potentials.
    subroutine test_greens_identity_3d(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_pars(2) = [41,81]                               ! points per field line, two resolutions
        integer, parameter :: n_alphas(2) = [40,80]                             ! field lines

        integer :: ierr, id                                                     ! error status, counter
        real(dp) :: res_lo(3), res_hi(3)                                        ! identity residuals at the two resolutions

        call greens_residuals_3d(n_pars(1),n_alphas(1),res_lo,ierr)
        call check(error, ierr, 0, 'assembly failed')
        if (allocated(error)) return
        call greens_residuals_3d(n_pars(2),n_alphas(2),res_hi,ierr)
        call check(error, ierr, 0, 'assembly failed')
        if (allocated(error)) return

        ! constant potential: exact by construction of the H diagonal
        call check(error, res_hi(1).lt.1.e-10_dp, &
            &'H 1 /= -4 pi 1: '//trim(r2str(res_hi(1))))
        if (allocated(error)) return

        ! interior harmonic: c = -4 pi; exterior harmonic: c = 0
        call check(error, res_hi(2).lt.5.e-2_dp, &
            &'interior identity residual too large: '//&
            &trim(r2str(res_hi(2))))
        if (allocated(error)) return
        call check(error, res_hi(3).lt.5.e-2_dp, &
            &'exterior identity residual too large: '//&
            &trim(r2str(res_hi(3))))
        if (allocated(error)) return

        ! convergence with resolution for the nontrivial potentials
        do id = 2,3
            call check(error, res_hi(id).lt.0.85_dp*res_lo(id), &
                &'no convergence for potential '//trim(i2str(id))//': '//&
                &trim(r2str(res_lo(id)))//' -> '//trim(r2str(res_hi(id))))
            if (allocated(error)) return
        end do
    end subroutine test_greens_identity_3d

    !> Assemble the style-1 matrices on the torus at a given resolution and
    !! return the identity residuals for the three test potentials, each
    !! measured against its established offset (constant and interior
    !! harmonic: c = -4 pi; exterior harmonic: c = 0). The full offset scan
    !! is printed for diagnosis.
    !!
    !! The products with the distributed matrices are evaluated globally on
    !! all processes (see fullstack_utils), so the assertions are identical
    !! everywhere; processes outside the BLACS context report zero residuals.
    subroutine greens_residuals_3d(n_par,n_alpha_loc,res,ierr)
        use num_vars, only: rank
        use vac_vars, only: in_context
        use fullstack_utils, only: dis_matvec

        integer, intent(in) :: n_par, n_alpha_loc                               ! resolution
        real(dp), intent(out) :: res(3)                                         ! residuals for the potentials
        integer, intent(out) :: ierr                                            ! error status

        real(dp), parameter :: x_src(3) = [R_0,0._dp,0._dp]                     ! source on magnetic axis for the exterior potential

        type(vac_type) :: vac                                                   ! vacuum variables
        integer :: id, kind                                                     ! counters
        integer :: n_bnd                                                        ! number of boundary points
        real(dp) :: r_src                                                       ! distance to source
        real(dp), allocatable :: phi(:), dphi(:)                                ! potential and norm . grad phi
        real(dp), allocatable :: lhs(:), rhs(:)                                 ! H phi and G dphi
        real(dp) :: res_off(-2:2)                                               ! residuals vs offsets c = k 2 pi
        real(dp) :: scale_phi                                                   ! norm of 4 pi phi

        res = 0._dp
        call setup_torus_vac(vac,n_par,n_alpha_loc,ierr)
        if (ierr.ne.0) return
        n_bnd = vac%n_bnd

        ierr = calc_GH(vac)
        if (ierr.ne.0) return

        if (in_context(vac%ctxt_HG)) then
            allocate(phi(n_bnd),dphi(n_bnd),lhs(n_bnd),rhs(n_bnd))

            do kind = 1,3
                do id = 1,n_bnd
                    select case (kind)
                        case (1)                                                ! constant (interior harmonic)
                            phi(id) = 1._dp
                            dphi(id) = 0._dp
                        case (2)                                                ! x (interior harmonic)
                            phi(id) = vac%x_vec(id,1)
                            dphi(id) = vac%norm(id,1)
                        case (3)                                                ! 1/|x - x_src| (exterior harmonic, source inside)
                            r_src = sqrt(sum((vac%x_vec(id,:)-x_src)**2))
                            phi(id) = 1._dp/r_src
                            dphi(id) = -sum(vac%norm(id,:)*&
                                &(vac%x_vec(id,:)-x_src))/r_src**3
                    end select
                end do
                ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
                if (ierr.ne.0) return
                ierr = dis_matvec(vac,vac%G,vac%desc_G,dphi,rhs)
                if (ierr.ne.0) return
                lhs = lhs - rhs
                scale_phi = maxval(abs(4._dp*pi*phi))
                do id = -2,2
                    res_off(id) = maxval(abs(lhs - id*2._dp*pi*phi))/scale_phi
                end do
                if (rank.eq.0) write(error_unit,'(A,I6,A,I2,A,5ES10.2)') &
                    &'   [3d identity] n = ',n_bnd,', potential ',kind,&
                    &': residuals vs c = (-4,-2,0,2,4)pi: ',res_off
                select case (kind)
                    case (1,2)                                                  ! interior: c = -4 pi
                        res(kind) = res_off(-2)
                    case (3)                                                    ! exterior: c = 0
                        res(kind) = res_off(0)
                end select
            end do
        end if

        call vac%dealloc()
    end subroutine greens_residuals_3d

    !> The vacuum response of the axisymmetric circular boundary, calculated
    !! with the field-line 3-D machinery (style 1), against the same response
    !! from the axisymmetric machinery (style 2), which is itself verified
    !! against the analytical cylinder limit elsewhere.
    !!
    !! The two styles use entirely different kernels (free-space Green's
    !! function on field lines vs analytically toroidally integrated
    !! Q_{n-1/2}), different singular-integral treatments and different
    !! H-diagonal constructions, so agreement here validates the whole
    !! style-1 chain end to end. The Fourier data are identical:
    !! exp(i n alpha + i (n q - m) theta) = exp(i n zeta - i m theta).
    subroutine test_response_1_vs_2(error)
        use X_vars, only: n_mod_X, modes_type
        use num_vars, only: use_pol_flux_F, eq_style, rank
        use vac_ops, only: calc_vac_res
        use fullstack_utils, only: gather_res

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_par = 61                                        ! points per field line
        integer, parameter :: n_alpha_loc = 60                                  ! field lines
        integer, parameter :: n_bnd_2 = 201                                     ! boundary points for style 2
        integer, parameter :: n_mod = 5                                         ! poloidal modes m = 1..n_mod
        real(dp), parameter :: tol_diag = 0.15_dp                               ! relative tolerance on the diagonal
        real(dp), parameter :: tol_offdiag = 0.10_dp                            ! off-diagonal tolerance, rel. to largest diagonal

        type(vac_type) :: vac1, vac2                                            ! style-1 and style-2 vacuum variables
        type(modes_type) :: mds                                                 ! minimal modes variables
        integer :: ierr, id, jd                                                 ! error status, counters
        integer :: n_mod_X_old                                                  ! original n_mod_X
        real(dp) :: t, R                                                        ! angle, major radius
        real(dp) :: rel_diff                                                    ! relative difference
        real(dp) :: diag_max                                                    ! largest diagonal magnitude
        complex(dp), allocatable :: res1(:,:), res2(:,:)                        ! responses

        n_mod_X_old = n_mod_X
        n_mod_X = n_mod
        use_pol_flux_F = .true.
        eq_style = 1                                                            ! VMEC-like (style 1); only affects output naming

        allocate(mds%m(1,n_mod),mds%n(1,n_mod))
        mds%m(1,:) = [(id, id=1,n_mod)]
        mds%n(1,:) = 2

        ! style 1: field-line covered torus
        call setup_torus_vac(vac1,n_par,n_alpha_loc,ierr)
        call check(error, ierr, 0, 'style-1 setup failed')
        if (allocated(error)) return
        ierr = calc_vac_res(mds,vac1)
        call check(error, ierr, 0, 'style-1 response failed')
        if (allocated(error)) return
        allocate(res1(n_mod,n_mod))
        ierr = gather_res(vac1,res1)                                            ! response lives on the last process only
        call check(error, ierr, 0, 'gathering style-1 response failed')
        if (allocated(error)) return

        ! style 2: same boundary, axisymmetric machinery, same q
        ierr = vac2%init(2,n_bnd_2,2,[n_bnd_2,1],q_saf)
        call check(error, ierr, 0, 'style-2 init failed')
        if (allocated(error)) return
        do id = 1,n_bnd_2
            t = 2._dp*pi*(id-1)/(n_bnd_2-1)
            R = R_0 + a_min*cos(t)
            vac2%ang(id,1) = t
            vac2%x_vec(id,:) = [R, a_min*sin(t)]
            vac2%norm(id,:) = [-R*a_min*cos(t), -R*a_min*sin(t)]
            vac2%dnorm(id,:) = [a_min*sin(t)*(a_min*cos(t)+R), &
                &a_min*(a_min*sin(t)**2-R*cos(t))]
        end do
        eq_style = 2
        ierr = calc_vac_res(mds,vac2)
        call check(error, ierr, 0, 'style-2 response failed')
        if (allocated(error)) return
        allocate(res2(n_mod,n_mod))
        ierr = gather_res(vac2,res2)
        call check(error, ierr, 0, 'gathering style-2 response failed')
        if (allocated(error)) return

        do id = 1,n_mod
            if (rank.ne.0) exit
            write(error_unit,'(A,I2,A,ES12.5,A,ES12.5,A,F7.3)') &
                &'   [response 1 vs 2] m = ',mds%m(1,id),&
                &': style 1 = ',real(res1(id,id)),&
                &', style 2 = ',real(res2(id,id)),&
                &', ratio = ',real(res1(id,id))/real(res2(id,id))
        end do

        diag_max = 0._dp
        do id = 1,n_mod
            diag_max = max(diag_max,abs(res2(id,id)))
        end do

        do id = 1,n_mod
            do jd = 1,n_mod
                if (id.eq.jd) then
                    rel_diff = abs(res1(id,id)-res2(id,id))/abs(res2(id,id))
                    call check(error, rel_diff.lt.tol_diag, &
                        &'diagonal response mismatch for m = '//&
                        &trim(i2str(id))//': '//trim(r2str(rel_diff)))
                else
                    call check(error, abs(res1(id,jd)).lt.tol_offdiag*&
                        &diag_max, &
                        &'style-1 off-diagonal too large at ('//&
                        &trim(i2str(id))//','//trim(i2str(jd))//'): '//&
                        &trim(r2str(abs(res1(id,jd)))))
                end if
                if (allocated(error)) exit
            end do
            if (allocated(error)) exit
        end do

        ! clean up and restore module state
        call vac1%dealloc()
        call vac2%dealloc()
        call mds%dealloc()
        n_mod_X = n_mod_X_old
        if (allocated(error)) return
    end subroutine test_response_1_vs_2
end module test_vac_3d
