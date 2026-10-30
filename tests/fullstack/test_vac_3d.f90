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

contains
    !> Collect the tests.
    subroutine collect_vac_3d(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("GH_int_1_singular", test_int_1_singular), &
            &new_unittest("greens_identity_3d", test_greens_identity_3d), &
            &new_unittest("response_1_vs_2", test_response_1_vs_2), &
            &new_unittest("point_source_3d", test_point_source_3d), &
            &new_unittest("exact_exterior_2d", test_exact_exterior_2d) &
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

    !> Manufactured exterior potential through the whole style-1 chain.
    !!
    !! A point source on the magnetic axis, phi = 1/|x - x_src|, is harmonic
    !! outside the surface and decays, and its normal data norm . grad phi
    !! are exact. solve_Phi_BEM (style 1: no diagonal shift, the row-sum
    !! diagonal of calc_GH_1 already selects the exterior) must then return
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
        use fullstack_utils, only: n_col1_loc, vec_glob2dis

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
    !! assembled operators must satisfy H phi = G dphi for the Q family and
    !! (H + 4 pi) phi = G dphi for the P family, and solve_Phi_BEM must return
    !! the P-family boundary traces, with errors that fall with the number of
    !! boundary points. This pins the axisymmetric vacuum at finite aspect
    !! ratio, poloidal mode coupling included, where the cylinder test (the
    !! limit a/R_0 -> 0) and MISHKA (1 %) cannot. Measured on 31 October 2026
    !! (issue #4): identity residuals 1e-3 to 1.2e-2 and solution errors 0.9 to
    !! 2.4 % at 301 points for n = 2, k = 1 to 5 and n = 10, k = 6 to 15,
    !! halving at 601 points.
    subroutine test_exact_exterior_2d(error)
        use num_vars, only: eq_style, rank
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: vec_dis2loc
        use vac_ops, only: solve_Phi_BEM
        use fullstack_utils, only: n_col1_loc, vec_glob2dis, dis_matvec

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_bnds(2) = [301,601]                             ! boundary points, two resolutions
        integer, parameter :: n_tors(2) = [2,10]                                ! toroidal mode numbers
        integer, parameter :: k_min(2) = [1,6]                                  ! lowest poloidal mode number per n
        integer, parameter :: k_max(2) = [5,15]                                 ! highest poloidal mode number per n
        real(dp), parameter :: tol_fine = 3.e-2_dp                              ! residuals and solution error at the fine resolution (measured: up to 1.5e-2)
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

                        ! interior family: H phi = G dphi
                        do id = 1,n_bnd
                            call tor_harm(vac%x_vec(id,1),vac%x_vec(id,2),&
                                &vac%norm(id,:),n_tor,kd,.false.,phi(id),dphi(id))
                        end do
                        ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
                        if (ierr.eq.0) ierr = dis_matvec(vac,vac%G,vac%desc_G,&
                            &dphi,rhs)
                        call check(error, ierr, 0, 'matrix-vector product failed')
                        if (allocated(error)) exit tor
                        errs(1,ires) = maxval(abs(lhs-rhs))/maxval(abs(4._dp*pi*phi))

                        ! exterior family: (H + 4 pi) phi = G dphi, and the solve
                        do id = 1,n_bnd
                            call tor_harm(vac%x_vec(id,1),vac%x_vec(id,2),&
                                &vac%norm(id,:),n_tor,kd,.true.,phi(id),dphi(id))
                        end do
                        ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
                        if (ierr.eq.0) ierr = dis_matvec(vac,vac%G,vac%desc_G,&
                            &dphi,rhs)
                        call check(error, ierr, 0, 'matrix-vector product failed')
                        if (allocated(error)) exit tor
                        errs(2,ires) = maxval(abs(lhs+4._dp*pi*phi-rhs))/&
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
            ! the axisymmetric kernels take the normal data with the opposite
            ! sign of the field-line ones: -norm . grad phi, with norm the
            ! (inward) J nabla psi of the fixture; both identities below hold
            ! to 1e-3 with this sign and fail at order one with the other
            dphi = -(norm(1)*(tor_val(R+h,Z,n,k,exterior)-&
                &tor_val(R-h,Z,n,k,exterior))/(2._dp*h) + &
                &norm(2)*(tor_val(R,Z+h,n,k,exterior)-&
                &tor_val(R,Z-h,n,k,exterior))/(2._dp*h))
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
        ! the ratios at 81 x 80 as recorded on 30 October 2026 (gfortran 14 and
        ! 16, ScaLAPACK and STRUMPACK agree to 1e-6); asserted to 0.5 % so that
        ! a change of the surface weights shows: the field-line spacing
        ! 2 pi / (n_alpha - 1) instead of 2 pi / n_alpha moves them by 1.25 %,
        ! inside tol_diag but not inside this (issue #4, step 7 mutation check)
        real(dp), parameter :: ratio_rec(5) = &
            &[0.965_dp,0.946_dp,0.926_dp,0.906_dp,0.886_dp]
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
end module test_vac_3d
