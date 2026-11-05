!------------------------------------------------------------------------------!
!> Full-stack tests of the 3-D (field-line aligned, style 1) vacuum: the
!! singular interval kernel, the assembled G and H matrices through Green's
!! identities, and the vacuum response against the (independently verified)
!! axisymmetric style-2 calculation of the same boundary.
!!
!! The geometry is an analytical circular torus
!!    x = (R cos zeta, R sin zeta, a sin theta),   R = R_0 + a cos theta,
!! covered by field lines zeta = alpha + q theta on a product grid
!! (theta, alpha) in [0,2 pi] x [0,2 pi), theta with both ends included
!! (trapezoidal half-weights on the duplicated edge) and alpha without its
!! end point, exactly as PB3D lays out its field lines (rich_ops:
!! calc_eqd_grid with excl_last; the map is a shear, so this tiles the torus
!! once). In PB3D's flux coordinates
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
    ! 8-point Gauss-Legendre rule on [-1,1], for brute-force references
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

contains
    !> Collect the tests.
    subroutine collect_vac_3d(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("GH_int_1_singular", test_int_1_singular), &
            &new_unittest("GH_int_1_near", test_int_1_near), &
            &new_unittest("greens_identity_3d", test_greens_identity_3d), &
            &new_unittest("response_1_vs_2", test_response_1_vs_2), &
            &new_unittest("point_source_3d", test_point_source_3d), &
            &new_unittest("exact_exterior_2d", test_exact_exterior_2d), &
            &new_unittest("solid_angle_node", test_solid_angle_node), &
            &new_unittest("H_sing_node", test_H_sing_node), &
            &new_unittest("diag_direct_3d", test_diag_direct_3d), &
            &new_unittest("point_source_direct_3d", test_point_source_direct) &
            &]
    end subroutine collect_vac_3d

    !--------------------------------------------------------------------------
    ! helpers
    !--------------------------------------------------------------------------

    !> Position and flux-coordinate tangents on the field line, for the
    !! circular torus or, with a ripple amplitude delta, for the rippled torus
    !! with minor radius a(zeta) = a_min (1 + delta cos(N_rip zeta)).
    pure subroutine torus(theta,alpha,x,e_theta,e_alpha,delta)
        real(dp), intent(in) :: theta, alpha                                    ! angles
        real(dp), intent(out) :: x(3)                                           ! position
        real(dp), intent(out), optional :: e_theta(3), e_alpha(3)               ! tangents
        real(dp), intent(in), optional :: delta                                 ! ripple amplitude

        integer, parameter :: N_rip = 5                                         ! ripple periods
        real(dp) :: R, zeta                                                     ! major radius, toroidal angle
        real(dp) :: a, da                                                       ! minor radius and its zeta derivative
        real(dp) :: dR_t, dz_t, dR_a, dz_a                                      ! derivatives of R and z along theta and alpha

        zeta = alpha + q_saf*theta
        a = a_min
        da = 0._dp
        if (present(delta)) then
            a = a_min*(1._dp + delta*cos(N_rip*zeta))
            da = -a_min*delta*N_rip*sin(N_rip*zeta)
        end if
        R = R_0 + a*cos(theta)
        x = [R*cos(zeta), R*sin(zeta), a*sin(theta)]
        if (present(e_theta)) then
            dR_t = da*q_saf*cos(theta) - a*sin(theta)
            dz_t = da*q_saf*sin(theta) + a*cos(theta)
            e_theta = [dR_t*cos(zeta) - R*sin(zeta)*q_saf, &
                &dR_t*sin(zeta) + R*cos(zeta)*q_saf, dz_t]
        end if
        if (present(e_alpha)) then
            dR_a = da*cos(theta)
            dz_a = da*sin(theta)
            e_alpha = [dR_a*cos(zeta) - R*sin(zeta), &
                &dR_a*sin(zeta) + R*cos(zeta), dz_a]
        end if
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
    subroutine setup_torus_vac(vac,n_par,n_alpha_loc,ierr,delta)
        use grid_vars, only: min_par_X, max_par_X, min_alpha, max_alpha, n_alpha
        use rich_vars, only: n_par_X

        type(vac_type), intent(inout) :: vac                                    ! vacuum variables
        integer, intent(in) :: n_par, n_alpha_loc                               ! points per field line, number of field lines
        integer, intent(out) :: ierr                                            ! error status
        real(dp), intent(in), optional :: delta                                 ! ripple amplitude

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
            alpha = 2._dp*pi*(jd-1)/n_alpha_loc                                 ! end point excluded, as in PB3D
            do id = 1,n_par
                theta = 2._dp*pi*(id-1)/(n_par-1)
                kd = id + (jd-1)*n_par
                call torus(theta,alpha,x,e_theta=e_t,e_alpha=e_a,delta=delta)
                vac%x_vec(kd,:) = x
                vac%norm(kd,:) = cross(e_t,e_a)                                 ! J nabla psi
                vac%h_fac(kd,1) = sum(e_a**2)                                   ! g_alpha,alpha
                vac%h_fac(kd,2) = sum(e_a*e_t)                                  ! g_alpha,theta
                vac%h_fac(kd,3) = sum(e_t**2)                                   ! g_theta,theta
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

        real(dp), parameter :: h_fac(3) = [4.0_dp,0.9_dp,2.5_dp]                ! g_aa, g_at, g_tt
        real(dp), parameter :: steps(2) = [0.02_dp,0.03_dp]                     ! dpar, dalpha

        integer :: id, jd                                                       ! counters
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(4,3), x_in(3), norm_s(4,3)                              ! geometry (only distances matter)
        real(dp) :: d(3), dn(3)                                                 ! interval vector and normal increment
        real(dp) :: s_no, sq_a, sq_t                                            ! nonorthogonality and metric roots
        real(dp) :: I_ref                                                       ! brute-force integral
        real(dp) :: p_lo, p_hi, u_q, w_q                                        ! quadrature panels and nodes
        real(dp) :: V_hc                                                        ! alpha half-cell in metric units
        real(dp) :: G_ref                                                       ! reference value
        real(dp) :: r_far                                                       ! length of the interval

        ! a singular interval on a straight stencil: influence point = stencil
        ! point 2 = left point of the interval; the right point sits at par
        ! distance steps(1) in the field-line direction, the normal varies
        ! linearly (the cubic interpolant reproduces both exactly)
        x_in = [1.0_dp,2.0_dp,3.0_dp]
        d = sqrt(h_fac(3))*steps(1)*[1._dp,0._dp,0._dp]
        dn = [0._dp,0._dp,0.1_dp]
        do id = 1,4
            x_s(id,:) = x_in + (id-2)*d
            norm_s(id,:) = [0.4_dp,0.5_dp,0.6_dp] + (id-2)*dn
        end do

        call calc_GH_int_1(G,H,x_s,x_in,norm_s,h_fac,steps,1.e-30_dp,2)

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

        ! the singular endpoint carries no H contribution: calc_GH_1 sets the
        ! diagonal of H afterwards
        call check(error, H(1), 0._dp, thr=1.e-300_dp, &
            &message='singular H not zero')
        if (allocated(error)) return

        ! the regular endpoint is near (one interval length) and receives the
        ! kernel integrated over its half of the interval: on the straight
        ! line with r = t |d|, 2 int_{1/2}^1 (-1/r) dt = -2 ln 2 / |d| and
        ! 2 int_{1/2}^1 (n_0 + t dn).(t d) / r^3 dt
        !    = 2 (n_0.d + ln 2 dn.d) / |d|^3
        r_far = sqrt(sum(d**2))
        call check(error, G(2), -2._dp*log(2._dp)/r_far, thr=1.e-9_dp, &
            &rel=.true., message='near-end G: '//trim(r2str(G(2)))//' vs '//&
            &trim(r2str(-2._dp*log(2._dp)/r_far)))
        if (allocated(error)) return
        call check(error, H(2), 2._dp*(sum(norm_s(2,:)*d)+&
            &log(2._dp)*sum(dn*d))/r_far**3, thr=1.e-9_dp, rel=.true., &
            &message='near-end H: '//trim(r2str(H(2))))
        if (allocated(error)) return
    end subroutine test_int_1_singular

    !> The near-field branch of the 3-D interval kernel against brute-force
    !! quadrature over the exact curve: four stencil points on a circle with
    !! the radius of curvature of the W7-X bean tip (0.27 R_0) at the spacing
    !! of the 121-point field lines (0.06 R_0), the influence point a quarter
    !! of an interval length off the arc, where the kernel varies by a factor
    !! three across a half-interval. The reference integrates the kernels over
    !! the exact arc with panel-refined Gauss quadrature; the routine sees only
    !! the four points (cubic interpolant, 8 Gauss points per half-interval).
    subroutine test_int_1_near(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp), parameter :: R_c = 0.27_dp                                    ! radius of the arc
        real(dp), parameter :: ds = 0.06_dp                                     ! arc length between stencil points
        real(dp), parameter :: w = 0.015_dp                                     ! distance of the influence point from the arc
        real(dp), parameter :: tol = 3.e-3_dp                                   ! relative tolerance (measured 6e-4 for G, 1.2e-3 for H: the cubic interpolant of the arc is 1.3e-5 off in position, against the distance 0.015)
        integer, parameter :: n_pan = 32                                        ! panels per half-interval

        integer :: id, jd, kd, ip                                               ! counters
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(4,3), x_in(3), norm_s(4,3)                              ! geometry
        real(dp) :: G_ref(2), H_ref(2)                                          ! brute-force values
        real(dp) :: t_lo, t_hi, t_q, w_q                                        ! panel limits, node and weight
        real(dp) :: x_q(3), n_q(3), r_q                                         ! position, normal and distance at a node

        ! stencil on the arc, interval between points 2 and 3; the normal is
        ! the (outward) radial direction, scaled to 1.5
        do id = 1,4
            call arc((id-2)*ds,x_s(id,:),norm_s(id,:))
        end do
        ! influence point: foot at 0.3 of the interval, w off the arc inward
        call arc(0.3_dp*ds,x_in,n_q)
        x_in = x_in - w*n_q/1.5_dp

        call calc_GH_int_1(G,H,x_s,x_in,norm_s,[1._dp,0._dp,1._dp,0._dp],&
            &[ds,ds],1.e-30_dp,2)

        ! brute force over the exact arc: 2 int_half K ds/ds, Gauss on panels
        do kd = 1,2
            G_ref(kd) = 0._dp
            H_ref(kd) = 0._dp
            do ip = 1,n_pan
                t_lo = 0.5_dp*(kd-1) + 0.5_dp*(ip-1)/n_pan
                t_hi = 0.5_dp*(kd-1) + 0.5_dp*ip/n_pan
                do jd = 1,8
                    t_q = 0.5_dp*(t_lo+t_hi) + 0.5_dp*(t_hi-t_lo)*gl_x(jd)
                    w_q = 0.5_dp*(t_hi-t_lo)*gl_w(jd)
                    call arc(t_q*ds,x_q,n_q)
                    r_q = sqrt(sum((x_q-x_in)**2))
                    G_ref(kd) = G_ref(kd) - 2._dp*w_q/r_q
                    H_ref(kd) = H_ref(kd) + 2._dp*w_q*sum(n_q*(x_q-x_in))/r_q**3
                end do
            end do
            call check(error, G(kd), G_ref(kd), thr=tol*abs(G_ref(kd)), &
                &message='near G('//trim(i2str(kd))//'): '//trim(r2str(G(kd)))&
                &//' vs brute force '//trim(r2str(G_ref(kd))))
            if (allocated(error)) return
            call check(error, H(kd), H_ref(kd), thr=tol*abs(H_ref(kd)), &
                &message='near H('//trim(i2str(kd))//'): '//trim(r2str(H(kd)))&
                &//' vs brute force '//trim(r2str(H_ref(kd))))
            if (allocated(error)) return
        end do
    contains
        ! point and scaled outward normal on the arc at arc length s
        subroutine arc(s,x,n)
            real(dp), intent(in) :: s
            real(dp), intent(out) :: x(3), n(3)
            n = [cos(s/R_c),sin(s/R_c),0._dp]
            x = [2._dp,0.5_dp,-1._dp] + R_c*n
            n = 1.5_dp*n
        end subroutine arc
    end subroutine test_int_1_near

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

    !> Manufactured exterior potential through the whole style-1 chain.
    !!
    !! A point source on the magnetic axis, phi = 1/|x - x_src|, is harmonic
    !! outside the surface and decays, and its normal data norm . grad phi
    !! are exact. solve_Phi_BEM (the row-sum diagonal of calc_GH_1 is the
    !! vacuum-side corner term, so H phi = G dphi holds) must then return
    !! phi on the surface to discretization accuracy, and the error must fall
    !! with resolution. Run on the circular torus and on a rippled torus
    !! a(zeta) = a (1 + 0.1 cos 5 zeta), both with PB3D's field-line grid
    !! (alpha without its end point). This pins the orientation of the
    !! normal, the alpha weights and the exterior solve at once (issue #4,
    !! step 6b): with the old alpha spacing (n_alpha - 1) and a grid with the
    !! end point duplicated, the system had two identical rows.
    subroutine test_point_source_3d(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_pars(2) = [41,81]                               ! points per field line, two resolutions
        integer, parameter :: n_alphas(2) = [40,80]                             ! field lines
        real(dp), parameter :: deltas(2) = [0._dp,0.1_dp]                       ! ripple amplitudes
        real(dp), parameter :: tol = 5.e-2_dp                                   ! relative tolerance at the fine resolution

        integer :: ierr, kind, ires                                             ! error status, counters
        real(dp) :: err(2)                                                      ! solution errors at the two resolutions

        do kind = 1,2
            do ires = 1,2
                call point_source_error(n_pars(ires),n_alphas(ires),&
                    &deltas(kind),err(ires),ierr)
                call check(error, ierr, 0, 'point-source solve failed')
                if (allocated(error)) return
            end do
            call check(error, err(2).lt.tol, &
                &'point-source solution error too large for ripple '//&
                &trim(r2str(deltas(kind)))//': '//trim(r2str(err(2))))
            if (allocated(error)) return
            call check(error, err(2).lt.0.85_dp*err(1), &
                &'no convergence for ripple '//trim(r2str(deltas(kind)))//&
                &': '//trim(r2str(err(1)))//' -> '//trim(r2str(err(2))))
            if (allocated(error)) return
        end do
    end subroutine test_point_source_3d

    !> Solve the style-1 exterior problem for the point-source potential on
    !! the (rippled) torus at a given resolution and return the relative
    !! error of the boundary trace, max |Phi - phi| / max |phi|. Processes
    !! outside the BLACS context report zero.
    subroutine point_source_error(n_par,n_alpha_loc,delta,err,ierr)
        use num_vars, only: rank
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: vec_dis2loc
        use vac_ops, only: solve_Phi_BEM
        use vac_utilities, only: n_col1_loc, vec_glob2dis

        integer, intent(in) :: n_par, n_alpha_loc                               ! resolution
        real(dp), intent(in) :: delta                                           ! ripple amplitude
        real(dp), intent(out) :: err                                            ! relative solution error
        integer, intent(out) :: ierr                                            ! error status

        real(dp), parameter :: x_src(3) = [R_0,0._dp,0._dp]                     ! source on the magnetic axis

        type(vac_type) :: vac                                                   ! vacuum variables
        integer :: id                                                           ! counter
        integer :: n_bnd                                                        ! number of boundary points
        integer :: n_col1                                                       ! local number of columns of R and Phi
        integer :: desc_RPhi(9)                                                 ! descriptor for R and Phi
        integer, allocatable :: lims_c_R(:,:)                                   ! local column limits of R and Phi
        real(dp) :: r_src                                                       ! distance to source
        real(dp), allocatable :: phi(:), dphi(:)                                ! potential and norm . grad phi
        real(dp), allocatable :: R_mat(:,:), Phi_mat(:,:)                       ! distributed right-hand side and solution
        real(dp), allocatable :: phi_num(:)                                     ! gathered numerical solution

        err = 0._dp
        call setup_torus_vac(vac,n_par,n_alpha_loc,ierr,delta=delta)
        if (ierr.ne.0) return
        n_bnd = vac%n_bnd

        ierr = calc_GH(vac)
        if (ierr.ne.0) return

        if (in_context(vac%ctxt_HG)) then
            allocate(phi(n_bnd),dphi(n_bnd))
            do id = 1,n_bnd
                r_src = sqrt(sum((vac%x_vec(id,:)-x_src)**2))
                phi(id) = 1._dp/r_src
                dphi(id) = -sum(vac%norm(id,:)*(vac%x_vec(id,:)-x_src))/&
                    &r_src**3
            end do

            n_col1 = n_col1_loc(vac)
            call set_loc_lims(n_col1,vac%bs,vac%ind_p(2),vac%n_p(2),lims_c_R)
            allocate(R_mat(max(1,vac%n_loc(1)),n_col1))
            allocate(Phi_mat(max(1,vac%n_loc(1)),n_col1))
            if (n_col1.gt.0) call vec_glob2dis(vac,dphi,R_mat(:,1))
            Phi_mat = 0._dp
            call descinit(desc_RPhi,n_bnd,1,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
                &max(1,vac%n_loc(1)),ierr)
            if (ierr.ne.0) return

            ierr = solve_Phi_BEM(vac,R_mat,Phi_mat,[n_bnd,1],&
                &[vac%n_loc(1),n_col1],lims_c_R,desc_RPhi)
            if (ierr.ne.0) return

            allocate(phi_num(n_bnd))
            ierr = vec_dis2loc(vac%ctxt_HG,&
                &reshape(Phi_mat(1:vac%n_loc(1),:),[vac%n_loc(1)*n_col1]),&
                &vac%lims_r,phi_num)
            if (ierr.ne.0) return

            err = maxval(abs(phi_num-phi))/maxval(abs(phi))
            if (rank.eq.0) write(error_unit,'(A,I6,A,F5.2,A,ES10.2)') &
                &'   [point source] n = ',n_bnd,', ripple ',delta,&
                &': max |Phi - phi| / max |phi| = ',err
        end if

        call vac%dealloc()
    end subroutine point_source_error

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
        use vac_utilities, only: dis_matvec

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

    !> Exact harmonics of the circular torus through the axisymmetric (style
    !! 2) chain: catalogue item T4 of issue #4.
    !!
    !! In toroidal coordinates (eta, theta_t, phi) with focal radius
    !! c = sqrt(R_0^2 - a^2) the circular torus is the surface cosh eta_0 =
    !! R_0 / a, and
    !!    sqrt(cosh eta - cos theta_t) F^n_{k-1/2}(cosh eta) cos(k theta_t) cos(n phi)
    !! is harmonic for the toroidal functions F = P (first kind: regular on
    !! the axis and at infinity, so an exterior potential) and F = Q (second
    !! kind: regular at the focal ring, so an interior one). Their boundary
    !! data are exact (the gradient by central differences in (R, Z) of the
    !! closed form, the toroidal functions from dtorh1 at order n), so the
    !! assembled operators must satisfy H phi - G dphi = -4 pi phi for the Q
    !! family and H phi = G dphi for the P family (one vacuum-side convention
    !! for both styles since November 2026), and solve_Phi_BEM must return
    !! the P-family boundary traces, with errors that fall with the number of
    !! boundary points. This pins the axisymmetric vacuum at finite aspect
    !! ratio, poloidal mode coupling included, where the cylinder test (the
    !! limit a/R_0 -> 0) and MISHKA (1 %) cannot. Measured on 31 October 2026
    !! (issue #4): identity residuals 1e-3 to 1.2e-2 and solution errors 0.9 to
    !! 2.4 % at 301 points for n = 2, k = 1 to 5 and n = 10, k = 6 to 15,
    !! halving at 601 points. With the plasma-side corner angle in calc_GH_2
    !! (2 November 2026), at 301 points: solution errors 0.19 to 1.07 %, P
    !! (exterior) identity residuals 1.0e-3 to 5.4e-3 (4.5e-3 to 1.2e-2
    !! before), Q (interior) ones 4.2e-3 to 1.2e-2 (1.0e-3 to 5.3e-3 before);
    !! all still halving at 601 points.
    subroutine test_exact_exterior_2d(error)
        use num_vars, only: eq_style, rank
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: vec_dis2loc
        use vac_ops, only: solve_Phi_BEM
        use vac_utilities, only: n_col1_loc, vec_glob2dis, dis_matvec

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_bnds(2) = [301,601]                             ! boundary points, two resolutions
        integer, parameter :: n_tors(2) = [2,10]                                ! toroidal mode numbers
        integer, parameter :: k_min(2) = [1,6]                                  ! lowest poloidal mode number per n
        integer, parameter :: k_max(2) = [5,15]                                 ! highest poloidal mode number per n
        real(dp), parameter :: tol_fine = 3.e-2_dp                              ! residuals and solution error at the fine resolution (measured: up to 1.5e-2, 7.7e-3 since 2 November 2026)
        real(dp), parameter :: conv = 0.8_dp                                    ! largest allowed error ratio fine / coarse (measured 0.50 to 0.72: first order)

        type(vac_type) :: vac                                                   ! style-2 vacuum variables
        integer :: ierr, id, kd, ires, itor                                     ! error status, counters
        integer :: n_bnd, n_tor                                                 ! current resolution and toroidal mode number
        integer :: n_col1                                                       ! local columns of the right-hand side
        integer :: desc_RPhi(9)                                                 ! descriptor of right-hand side and solution
        integer, allocatable :: lims_c_R(:,:)                                   ! local column limits
        integer :: eq_style_old                                                 ! saved module state
        real(dp) :: t, R                                                        ! angle and major radius
        real(dp) :: errs(3,2)                                                   ! interior residual, exterior residual, solution error at the two resolutions
        real(dp), allocatable :: phi(:), dphi(:), lhs(:), rhs(:), phi_num(:)    ! boundary data, products, solution
        real(dp), allocatable :: R_mat(:,:), Phi_mat(:,:)                       ! distributed right-hand side and solution

        eq_style_old = eq_style
        eq_style = 2                                                            ! axisymmetric machinery

        tor: do itor = 1,2
            n_tor = n_tors(itor)
            do kd = k_min(itor),k_max(itor)
                do ires = 1,2
                    n_bnd = n_bnds(ires)
                    errs(:,ires) = 0._dp

                    ierr = vac%init(2,n_bnd,n_tor,[n_bnd,1],q_saf)
                    call check(error, ierr, 0, 'style-2 init failed')
                    if (allocated(error)) exit tor
                    do id = 1,n_bnd
                        t = 2._dp*pi*(id-1)/(n_bnd-1)
                        R = R_0 + a_min*cos(t)
                        vac%ang(id,1) = t
                        vac%x_vec(id,:) = [R, a_min*sin(t)]
                        vac%norm(id,:) = [-R*a_min*cos(t), -R*a_min*sin(t)]
                        vac%dnorm(id,:) = [a_min*sin(t)*(a_min*cos(t)+R), &
                            &a_min*(a_min*sin(t)**2-R*cos(t))]
                    end do
                    ierr = calc_GH(vac)
                    call check(error, ierr, 0, 'style-2 matrices failed')
                    if (allocated(error)) exit tor

                    if (in_context(vac%ctxt_HG)) then
                        allocate(phi(n_bnd),dphi(n_bnd),lhs(n_bnd),rhs(n_bnd),&
                            &phi_num(n_bnd))

                        ! interior family: H phi - G dphi = -4 pi phi
                        do id = 1,n_bnd
                            call tor_harm(vac%x_vec(id,1),vac%x_vec(id,2),&
                                &vac%norm(id,:),n_tor,kd,.false.,phi(id),dphi(id))
                        end do
                        ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
                        if (ierr.eq.0) ierr = dis_matvec(vac,vac%G,vac%desc_G,&
                            &dphi,rhs)
                        call check(error, ierr, 0, 'matrix-vector product failed')
                        if (allocated(error)) exit tor
                        errs(1,ires) = maxval(abs(lhs-rhs+4._dp*pi*phi))/&
                            &maxval(abs(4._dp*pi*phi))

                        ! exterior family: H phi = G dphi, and the solve
                        do id = 1,n_bnd
                            call tor_harm(vac%x_vec(id,1),vac%x_vec(id,2),&
                                &vac%norm(id,:),n_tor,kd,.true.,phi(id),dphi(id))
                        end do
                        ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
                        if (ierr.eq.0) ierr = dis_matvec(vac,vac%G,vac%desc_G,&
                            &dphi,rhs)
                        call check(error, ierr, 0, 'matrix-vector product failed')
                        if (allocated(error)) exit tor
                        errs(2,ires) = maxval(abs(lhs-rhs))/&
                            &maxval(abs(4._dp*pi*phi))
                        n_col1 = n_col1_loc(vac)
                        call set_loc_lims(n_col1,vac%bs,vac%ind_p(2),vac%n_p(2),&
                            &lims_c_R)
                        allocate(R_mat(max(1,vac%n_loc(1)),n_col1))
                        allocate(Phi_mat(max(1,vac%n_loc(1)),n_col1))
                        if (n_col1.gt.0) call vec_glob2dis(vac,dphi,R_mat(:,1))
                        Phi_mat = 0._dp
                        call descinit(desc_RPhi,n_bnd,1,vac%bs,vac%bs,0,0,&
                            &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
                        if (ierr.eq.0) ierr = solve_Phi_BEM(vac,R_mat,Phi_mat,&
                            &[n_bnd,1],[vac%n_loc(1),n_col1],lims_c_R,desc_RPhi)
                        if (ierr.eq.0) ierr = vec_dis2loc(vac%ctxt_HG,&
                            &reshape(Phi_mat(1:vac%n_loc(1),:),&
                            &[vac%n_loc(1)*n_col1]),vac%lims_r,phi_num)
                        call check(error, ierr, 0, 'exterior solve failed')
                        if (allocated(error)) exit tor
                        errs(3,ires) = maxval(abs(phi_num-phi))/maxval(abs(phi))

                        if (rank.eq.0) write(error_unit,&
                            &'(A,I3,A,I3,A,I5,A,3ES10.2)') &
                            &'   [exact 2d] n = ',n_tor,', k = ',kd,&
                            &', n_bnd = ',n_bnd,&
                            &': residual interior, exterior, solve error ',&
                            &errs(:,ires)

                        deallocate(phi,dphi,lhs,rhs,phi_num,R_mat,Phi_mat,&
                            &lims_c_R)
                    end if
                    call vac%dealloc()
                end do

                ! processes outside the context report zeros and pass
                do id = 1,3
                    call check(error, errs(id,2).lt.tol_fine, &
                        &'n = '//trim(i2str(n_tor))//', k = '//trim(i2str(kd))//&
                        &': error '//trim(i2str(id))//' too large at '//&
                        &trim(i2str(n_bnds(2)))//' points: '//&
                        &trim(r2str(errs(id,2))))
                    if (allocated(error)) exit tor
                    call check(error, errs(id,2).le.conv*errs(id,1), &
                        &'n = '//trim(i2str(n_tor))//', k = '//trim(i2str(kd))//&
                        &': error '//trim(i2str(id))//' does not fall: '//&
                        &trim(r2str(errs(id,1)))//' -> '//trim(r2str(errs(id,2))))
                    if (allocated(error)) exit tor
                end do
            end do
        end do tor

        eq_style = eq_style_old
    contains
        !> Boundary value and normal data of the toroidal harmonic of order
        !! n and degree k - 1/2 at (R, Z): the gradient by central differences
        !! of the closed form (step 1e-5: truncation 1e-10, dtorh1 noise 1e-6).
        subroutine tor_harm(R,Z,norm,n,k,exterior,phi,dphi)
            real(dp), intent(in) :: R, Z                                        ! position
            real(dp), intent(in) :: norm(2)                                     ! normal vector (R, Z components)
            integer, intent(in) :: n, k                                         ! toroidal and poloidal mode numbers
            logical, intent(in) :: exterior                                     ! first kind (exterior) or second (interior)
            real(dp), intent(out) :: phi, dphi                                  ! value and normal derivative

            real(dp), parameter :: h = 1.e-5_dp                                 ! difference step

            phi = tor_val(R,Z,n,k,exterior)
            ! both styles take the normal data along the stored normal, the
            ! (inward) -J nabla psi of the fixture: dphi = norm . grad phi
            dphi = norm(1)*(tor_val(R+h,Z,n,k,exterior)-&
                &tor_val(R-h,Z,n,k,exterior))/(2._dp*h) + &
                &norm(2)*(tor_val(R,Z+h,n,k,exterior)-&
                &tor_val(R,Z-h,n,k,exterior))/(2._dp*h)
        end subroutine tor_harm

        !> The toroidal harmonic at (R, Z).
        real(dp) function tor_val(R,Z,n,k,exterior) result(val)
            use dtorh, only: dtorh1

            real(dp), intent(in) :: R, Z                                        ! position
            integer, intent(in) :: n, k                                         ! toroidal and poloidal mode numbers
            logical, intent(in) :: exterior                                     ! first kind (exterior) or second (interior)

            integer :: newn, ierr_loc                                           ! dtorh1 outputs
            real(dp) :: c, d1, d2, eta, s, cth, sth, th                         ! toroidal coordinates
            real(dp) :: pl(0:k), ql(0:k)                                        ! toroidal functions

            c = sqrt(R_0**2-a_min**2)
            d1 = sqrt((R+c)**2+Z**2)
            d2 = sqrt((R-c)**2+Z**2)
            eta = log(d1/d2)
            s = sqrt((R**2+Z**2-c**2)**2 + 4._dp*c**2*Z**2)
            cth = (R**2+Z**2-c**2)/s
            sth = 2._dp*c*Z/s
            th = atan2(sth,cth)
            ierr_loc = dtorh1(cosh(eta),n,k,pl,ql,newn)
            if (ierr_loc.ne.0 .or. newn.lt.k) then
                val = huge(1._dp)
                return
            end if
            if (exterior) then
                val = sqrt(cosh(eta)-cth)*pl(k)*cos(k*th)
            else
                val = sqrt(cosh(eta)-cth)*ql(k)*cos(k*th)
            end if
        end function tor_val
    end subroutine test_exact_exterior_2d

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

        ! the style-1 response converges at first order towards the style-2
        ! one (ratios 0.952 to 0.849 for m = 1 to 5 at 61 x 60, 0.972 to
        ! 0.908 at 101 x 100, with PB3D's field-line grid; issue #4, step 6b)
        integer, parameter :: n_par = 81                                        ! points per field line
        integer, parameter :: n_alpha_loc = 80                                  ! field lines
        integer, parameter :: n_bnd_2 = 301                                     ! boundary points for style 2
        integer, parameter :: n_mod = 5                                         ! poloidal modes m = 1..n_mod
        real(dp), parameter :: tol_diag = 0.15_dp                               ! relative tolerance on the diagonal
        real(dp), parameter :: tol_offdiag = 0.10_dp                            ! off-diagonal tolerance, rel. to largest diagonal
        ! the ratios at 81 x 80 as recorded on 31 October 2026 with the
        ! near-field quadrature along the field lines and the G diagonal from
        ! the linear-potential identity (before them: 0.965, 0.946, 0.926,
        ! 0.906, 0.886; gfortran 14 and 16, ScaLAPACK and STRUMPACK agree to
        ! 1e-6); asserted to 0.5 % so that a change of the surface weights
        ! shows: the field-line spacing 2 pi / (n_alpha - 1) instead of
        ! 2 pi / n_alpha moved the old ones by 1.25 %, inside tol_diag but not
        ! inside this (issue #4, step 7 mutation check); recorded again on 2
        ! November 2026, when the corner term of calc_GH_2 got the plasma-side
        ! angle: the style-2 response at these 300 boundary points dropped by
        ! 1.2 to 1.4 % (the estimate 4/N for N points around is 1.3 %), so
        ! the ratios rose from 0.994, 0.992, 0.988, 0.985, 0.982
        real(dp), parameter :: ratio_rec(5) = &
            &[1.006_dp,1.005_dp,1.002_dp,0.999_dp,0.996_dp]
        real(dp), parameter :: tol_ratio = 0.005_dp                             ! tolerance on the recorded ratios

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
                    if (allocated(error)) exit
                    rel_diff = abs(real(res1(id,id))/real(res2(id,id))-&
                        &ratio_rec(id))/ratio_rec(id)
                    call check(error, rel_diff.lt.tol_ratio, &
                        &'style 1 over style 2 ratio for m = '//&
                        &trim(i2str(id))//' moved from its recorded value '//&
                        &trim(r2str(ratio_rec(id)))//' by '//trim(r2str(rel_diff)))
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
    
    !> The exact solid angle at a faceted node (vac_utilities.solid_angle_node,
    !! eq. 45 of the 2019 manuscript) against three independent values: a flat
    !! node (2 pi), the axisymmetric corner of its eq. (47) with the tilts
    !! g_2 and g_4 of the two poloidal edges, whose plasma-side solid angle is
    !! 2 pi - 2 (g_2 + g_4) (erratum E7), and a node of the circular torus
    !! against Girard's theorem for the spherical quadrilateral bounded by the
    !! four facet planes (2 pi minus the sum of the signed angles between
    !! consecutive facet normals).
    subroutine test_solid_angle_node(error)
        use vac_utilities, only: solid_angle_node
        type(error_type), allocatable, intent(out) :: error
        
        real(dp), parameter :: g2 = 0.05_dp, g4 = 0.03_dp                       ! tilts of the axisymmetric corner
        real(dp), parameter :: th0 = 0.9_dp, al0 = 0.3_dp                       ! node on the torus
        real(dp), parameter :: dth = 0.2_dp, dal = 0.15_dp                      ! mesh steps
        integer :: id                                                           ! counter
        real(dp) :: u(4,3), z(3)                                                ! edges and axis
        real(dp) :: omega                                                       ! solid angle
        real(dp) :: x0(3), e_t(3), e_a(3)                                       ! node geometry
        real(dp) :: a(4,3), tau                                                 ! facet normals, turning angle sum
        
        ! flat node: the four edges in the tangent plane
        u(1,:) = [-1._dp,0.2_dp,0._dp]
        u(2,:) = [0.3_dp,-1._dp,0._dp]
        u(3,:) = [1._dp,0.1_dp,0._dp]
        u(4,:) = [-0.2_dp,1._dp,0._dp]
        z = [0._dp,0._dp,1._dp]
        call solid_angle_node(u,z,omega)
        call check(error, omega, 2._dp*pi, thr=1.e-12_dp, &
            &message='flat node: '//trim(r2str(omega)))
        if (allocated(error)) return
        
        ! axisymmetric corner: u_1 = -u_3 toroidal, u_2 and u_4 poloidal,
        ! tilted towards z by g2 and g4
        u(1,:) = [0._dp,1._dp,0._dp]
        u(2,:) = [-cos(g2),0._dp,sin(g2)]
        u(3,:) = [0._dp,-1._dp,0._dp]
        u(4,:) = [cos(g4),0._dp,sin(g4)]
        call solid_angle_node(u,z,omega)
        call check(error, omega, 2._dp*pi-2._dp*(g2+g4), thr=1.e-12_dp, &
            &message='axisymmetric corner: '//trim(r2str(omega))//' vs '//&
            &trim(r2str(2._dp*pi-2._dp*(g2+g4))))
        if (allocated(error)) return
        
        ! torus node against Girard: edges in counterclockwise order about
        ! the inward normal, facets between consecutive edges
        call torus(th0,al0,x0,e_theta=e_t,e_alpha=e_a)
        z = cross(e_t,e_a)                                                      ! inward for this fixture
        z = z/norm2(z)
        call torus(th0,al0-dal,u(1,:)); u(1,:) = u(1,:) - x0
        call torus(th0-dth,al0,u(2,:)); u(2,:) = u(2,:) - x0
        call torus(th0,al0+dal,u(3,:)); u(3,:) = u(3,:) - x0
        call torus(th0+dth,al0,u(4,:)); u(4,:) = u(4,:) - x0
        if (sum(cross(u(1,:),u(2,:))*z).lt.0._dp) then                          ! make the order counterclockwise about z
            u([2,4],:) = u([4,2],:)
        end if
        call solid_angle_node(u,z,omega)
        tau = 0._dp
        do id = 1,4
            a(id,:) = cross(u(id,:),u(mod(id,4)+1,:))
            a(id,:) = a(id,:)/norm2(a(id,:))
        end do
        do id = 1,4
            ! signed angle from a_{id-1} to a_id about the shared edge u_id
            tau = tau + atan2(sum(cross(a(mod(id+2,4)+1,:),a(id,:))*u(id,:))/&
                &norm2(u(id,:)),sum(a(mod(id+2,4)+1,:)*a(id,:)))
        end do
        call check(error, omega, 2._dp*pi-tau, thr=1.e-10_dp, &
            &message='torus node: (45) '//trim(r2str(omega))//' vs Girard '//&
            &trim(r2str(2._dp*pi-tau)))
        if (allocated(error)) return
    end subroutine test_solid_angle_node
    
    !> The cell integral of the H kernel at a node (vac_utilities.H_sing_node,
    !! the corrected eq. 69 of the 2019 manuscript, erratum E22) against
    !! brute-force quadrature of the true kernel over the cell of a torus node
    !! at two cell sizes (64-point Gauss-Legendre in the polar coordinates of
    !! the cell; 16 points leave 2.5 %): the ratio is 1 to 1e-5.
    subroutine test_H_sing_node(error)
        use vac_utilities, only: H_sing_node
        type(error_type), allocatable, intent(out) :: error
        
        real(dp), parameter :: th0 = 0.9_dp, al0 = 0.3_dp                       ! node on the torus
        real(dp), parameter :: steps(2,2) = reshape([0.04_dp,0.03_dp,&
            &0.02_dp,0.015_dp],[2,2])                                           ! (dpar, dalpha) at two sizes
        real(dp), parameter :: h = 1.e-4_dp                                     ! step of the second differences
        integer :: id, jd, kd, ld                                               ! counters
        real(dp) :: x0(3), e_t(3), e_a(3), n0(3)                                ! node geometry
        real(dp) :: xp(3), xm(3), xpp(3), xpm(3), xmp(3), xmm(3)                ! displaced points
        real(dp) :: g(3), dd(3,3)                                               ! metric factors and second derivatives
        real(dp) :: I_form(2), I_bf(2)                                          ! formula and brute force
        real(dp) :: dpar, dalpha                                                ! steps
        real(dp) :: t_c(4), t_lims(6)                                           ! corner angles and the pieces of the angular range
        real(dp) :: t, tt, w_t, r_max, r, w_r                                   ! quadrature variables
        real(dp) :: x(3), et(3), ea(3), nr(3), d(3)                             ! source geometry
        real(dp) :: ratio(2)                                                    ! brute force over formula
        integer, parameter :: n_q = 64                                          ! quadrature points per direction of the brute force
        real(dp) :: q_x(n_q), q_w(n_q)                                          ! Gauss-Legendre rule
        
        call gauss_legendre(n_q,q_x,q_w)
        call torus(th0,al0,x0,e_theta=e_t,e_alpha=e_a)
        n0 = cross(e_t,e_a)                                                     ! stored normal of the fixture (inward)
        g = [sum(e_a**2),sum(e_a*e_t),sum(e_t**2)]
        call torus(th0+h,al0,xp); call torus(th0-h,al0,xm)
        dd(1,:) = (xp-2*x0+xm)/h**2
        call torus(th0,al0+h,xp); call torus(th0,al0-h,xm)
        dd(3,:) = (xp-2*x0+xm)/h**2
        call torus(th0+h,al0+h,xpp); call torus(th0+h,al0-h,xpm)
        call torus(th0-h,al0+h,xmp); call torus(th0-h,al0-h,xmm)
        dd(2,:) = (xpp-xpm-xmp+xmm)/(4*h**2)
        
        do kd = 1,2
            dpar = steps(1,kd)
            dalpha = steps(2,kd)
            call H_sing_node(n0,g,dd,dpar,dalpha,I_form(kd))
            ! brute force in polar coordinates of the (theta, alpha) plane
            ! over the cell |theta-th0| <= dpar/2, |alpha-al0| <= dalpha/2:
            ! kernel times r is bounded; the angular range is split at the
            ! corners of the cell
            t_c = [atan2(dalpha,dpar), atan2(dalpha,-dpar), &
                &atan2(-dalpha,-dpar)+2*pi, atan2(-dalpha,dpar)+2*pi]
            t_lims = [0._dp, t_c, 2*pi]
            I_bf(kd) = 0._dp
            do id = 1,5
                do jd = 1,n_q
                    t = 0.5_dp*(t_lims(id)+t_lims(id+1)) + &
                        &0.5_dp*(t_lims(id+1)-t_lims(id))*q_x(jd)
                    w_t = 0.5_dp*(t_lims(id+1)-t_lims(id))*q_w(jd)
                    r_max = min(0.5_dp*dpar/max(abs(cos(t)),1.e-300_dp), &
                        &0.5_dp*dalpha/max(abs(sin(t)),1.e-300_dp))
                    do ld = 1,n_q
                        r = 0.5_dp*r_max*(1._dp+q_x(ld))
                        w_r = 0.5_dp*r_max*q_w(ld)
                        call torus(th0+r*cos(t),al0+r*sin(t),x,e_theta=et,&
                            &e_alpha=ea)
                        nr = cross(et,ea)
                        d = x - x0
                        tt = sum(nr*d)/norm2(d)**3
                        I_bf(kd) = I_bf(kd) + w_t*w_r*r*tt
                    end do
                end do
            end do
            ratio(kd) = I_bf(kd)/I_form(kd)
            if (rank_out()) write(error_unit,'(A,2F8.4,A,ES12.4,A,ES12.4,A,F8.5)') &
                &'   [H_sing_node] steps ',dpar,dalpha,': formula ',I_form(kd),&
                &' brute force ',I_bf(kd),' ratio ',ratio(kd)
        end do
        ! the cell is symmetric about the node, so the first correction to
        ! the leading-order kernel integrates to zero: the agreement is second
        ! order and already 1e-5 at these cell sizes
        do kd = 1,2
            call check(error, abs(ratio(kd)-1._dp).lt.1.e-3_dp, &
                &'cell integral off by more than 0.1 %: ratio '//&
                &trim(r2str(ratio(kd))))
            if (allocated(error)) return
        end do
    end subroutine test_H_sing_node
    
    !> The direct diagonal of H (solid angle plus cell integral,
    !! vac_utilities.calc_diag_direct) against the row-sum diagonal that
    !! calc_GH_1 uses, on the circular torus at two resolutions: the
    !! difference, which measures the quadrature error of the off-diagonal
    !! row, falls at first order.
    subroutine test_diag_direct_3d(error)
        use vac_vars, only: in_context
        use vac_utilities, only: calc_diag_direct, dis_diag, vac_steps
        type(error_type), allocatable, intent(out) :: error
        
        integer, parameter :: n_pars(2) = [41,81]                               ! points per field line
        integer, parameter :: n_alphas(2) = [40,80]                             ! field lines
        integer :: ierr, ires                                                   ! error status, counter
        type(vac_type) :: vac                                                   ! vacuum variables
        real(dp) :: dpar, dalpha                                                ! steps
        real(dp) :: d_max(2), d_mean(2), f_max(2), s_max(2), o_max(2)           ! max and mean smooth difference, max faceted difference, max cell integral, max angle deviation
        real(dp), allocatable :: diag_rs(:), diag_dir(:), diag_fac(:), I_sing(:), omega_v(:)  ! diagonals and their parts
        logical, allocatable :: avail(:)                                        ! direct value available
        
        do ires = 1,2
            call setup_torus_vac(vac,n_pars(ires),n_alphas(ires),ierr)
            call check(error, ierr, 0, 'setup failed')
            if (allocated(error)) return
            ierr = calc_GH(vac)
            call check(error, ierr, 0, 'assembly failed')
            if (allocated(error)) return
            allocate(diag_rs(vac%n_bnd),diag_dir(vac%n_bnd),diag_fac(vac%n_bnd),&
                &I_sing(vac%n_bnd),omega_v(vac%n_bnd),avail(vac%n_bnd))
            ierr = dis_diag(vac,vac%H,diag_rs)
            call check(error, ierr, 0, 'gathering the diagonal failed')
            if (allocated(error)) return
            call vac_steps(dpar,dalpha)
            ierr = calc_diag_direct(vac,dpar,dalpha,diag_dir,diag_fac,I_sing,&
                &omega_v,avail)
            call check(error, ierr, 0, 'direct diagonal failed')
            if (allocated(error)) return
            d_max(ires) = 0._dp
            d_mean(ires) = 0._dp
            f_max(ires) = 0._dp
            s_max(ires) = 0._dp
            o_max(ires) = 0._dp
            if (in_context(vac%ctxt_HG)) then
                d_max(ires) = maxval(abs(diag_rs-diag_dir)/(4*pi),mask=avail)
                d_mean(ires) = sum(abs(diag_rs-diag_dir)/(4*pi),mask=avail)/&
                    &count(avail)
                f_max(ires) = maxval(abs(diag_rs-diag_fac)/(4*pi),mask=avail)
                s_max(ires) = maxval(abs(I_sing)/(4*pi),mask=avail)
                o_max(ires) = maxval(abs(omega_v/(4*pi)-0.5_dp),mask=avail)
            end if
            if (rank_out()) write(error_unit,'(A,I6,A,5ES10.2)') &
                &'   [diag direct] n = ',vac%n_bnd,': max, mean |rowsum - &
                &smooth|/4pi, max |rowsum - faceted|/4pi, max |I_sing|/4pi, &
                &max |omega_v/4pi - 1/2| = ',d_max(ires),d_mean(ires),&
                &f_max(ires),s_max(ires),o_max(ires)
            deallocate(diag_rs,diag_dir,diag_fac,I_sing,omega_v,avail)
            call vac%dealloc()
        end do
        if (in_context(vac%ctxt_HG)) then
            call check(error, d_max(2).lt.0.1_dp, &
                &'row-sum and direct diagonals differ by more than 10 % of &
                &4 pi: '//trim(r2str(d_max(2))))
            if (allocated(error)) return
            call check(error, d_max(2).lt.0.75_dp*d_max(1), &
                &'no first-order convergence of the diagonal difference: '//&
                &trim(r2str(d_max(1)))//' -> '//trim(r2str(d_max(2))))
            if (allocated(error)) return
        end if
    end subroutine test_diag_direct_3d
    
    !> The point-source solve of the field-line style with the direct
    !! diagonals (vac_diag_style 2: -2 pi plus the cell integral for H, the
    !! analytical cell integral for G) against the default indirect ones, at
    !! two resolutions: both must converge, and their errors are printed for
    !! the comparison of the two options.
    subroutine test_point_source_direct(error)
        use num_vars, only: vac_diag_style
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n_pars(2) = [41,81]                               ! points per field line
        integer, parameter :: n_alphas(2) = [40,80]                             ! field lines
        real(dp), parameter :: tol = 5.e-2_dp                                   ! relative tolerance at the fine resolution
        integer :: ierr, ires, istyle                                           ! error status, counters
        integer :: style_old                                                    ! saved module state
        real(dp) :: err(2,2)                                                    ! solution errors (resolution, option)
        
        style_old = vac_diag_style
        do istyle = 1,2
            vac_diag_style = istyle
            do ires = 1,2
                call point_source_error(n_pars(ires),n_alphas(ires),0._dp,&
                    &err(ires,istyle),ierr)
                if (ierr.ne.0) exit
            end do
            if (ierr.ne.0) exit
        end do
        vac_diag_style = style_old
        call check(error, ierr, 0, 'point-source solve failed')
        if (allocated(error)) return
        if (rank_out()) write(error_unit,'(A,2ES10.2,A,2ES10.2)') &
            &'   [point source direct] errors indirect ',err(:,1),&
            &', direct ',err(:,2)
        call check(error, err(2,2).lt.tol, &
            &'point-source solution error with the direct diagonals too &
            &large: '//trim(r2str(err(2,2))))
        if (allocated(error)) return
        call check(error, err(2,2).lt.0.85_dp*err(1,2), &
            &'no convergence with the direct diagonals: '//&
            &trim(r2str(err(1,2)))//' -> '//trim(r2str(err(2,2))))
        if (allocated(error)) return
    end subroutine test_point_source_direct
    

    !> Gauss-Legendre rule with n points on [-1, 1], by Newton iteration on
    !! the Legendre polynomial (Numerical Recipes, gauleg).
    subroutine gauss_legendre(n,x,w)
        integer, intent(in) :: n                                                ! number of points
        real(dp), intent(out) :: x(n), w(n)                                     ! abscissas and weights
        
        integer :: id, jd                                                       ! counters
        real(dp) :: z, z1, p1, p2, p3, pp                                       ! Newton iterate and polynomial values
        
        do id = 1,(n+1)/2
            z = cos(pi*(id-0.25_dp)/(n+0.5_dp))
            do
                p1 = 1._dp
                p2 = 0._dp
                do jd = 1,n
                    p3 = p2
                    p2 = p1
                    p1 = ((2*jd-1)*z*p2 - (jd-1)*p3)/jd
                end do
                pp = n*(z*p1-p2)/(z**2-1._dp)
                z1 = z
                z = z1 - p1/pp
                if (abs(z-z1).lt.1.e-15_dp) exit
            end do
            x(id) = -z
            x(n+1-id) = z
            w(id) = 2._dp/((1._dp-z**2)*pp**2)
            w(n+1-id) = w(id)
        end do
    end subroutine gauss_legendre
    

    !> Whether this process prints test diagnostics.
    logical function rank_out()
        use num_vars, only: rank
        rank_out = rank.eq.0
    end function rank_out
end module test_vac_3d
