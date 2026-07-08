!------------------------------------------------------------------------------!
!> Full-stack tests of the vacuum Green's function interval kernels
!! (vac_utilities.calc_GH_int_1 / calc_GH_int_2) against brute-force
!! references.
!!
!! The references are constructed independently of the implementation:
!!  - the regular G kernel is checked against a direct evaluation of the
!!    toroidal harmonics (dtorh),
!!  - the regular H kernel is checked against a numerical directional
!!    derivative of the G kernel along the source normal, which validates the
!!    analytical derivative algebra (the Aij helper variable) without
!!    re-deriving it,
!!  - the (near-)singular analytical integrals are checked against adaptive
!!    numerical quadrature of the true kernel, using the independently
!!    verified asymptote Q_{n-1/2}(1+x) ~ -1/2 ln(x/32) - b_n (see
!!    test_asymptote) to subtract the logarithmic singularity exactly.
!!
!! All geometry is an analytical circular toroidal boundary
!!    R(t) = R_0 + a cos t,   Z(t) = a sin t,
!! for which the (unnormalized) normal vector in PB3D convention is
!!    norm = -R (Z_t, -R_t) = (-R a cos t, -R a sin t).
!------------------------------------------------------------------------------!
module test_vac_kernels
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi
    use str_utilities, only: r2str, i2str
    use dtorh, only: dtorh1
    use vac_utilities, only: calc_GH_int_1, calc_GH_int_2

    implicit none
    private
    public collect_vac_kernels

    ! circular boundary parameters
    real(dp), parameter :: R_0 = 3.0_dp                                         ! major radius
    real(dp), parameter :: a_min = 1.0_dp                                       ! minor radius

    ! 16-point Gauss-Legendre rule on [-1,1]
    ! (Abramowitz & Stegun table 25.4, symmetric)
    real(dp), parameter :: gl_x(16) = [&
        &-0.9894009349916499_dp, -0.9445750230732326_dp, &
        &-0.8656312023878318_dp, -0.7554044083550030_dp, &
        &-0.6178762444026438_dp, -0.4580167776572274_dp, &
        &-0.2816035507792589_dp, -0.0950125098376374_dp, &
        & 0.0950125098376374_dp,  0.2816035507792589_dp, &
        & 0.4580167776572274_dp,  0.6178762444026438_dp, &
        & 0.7554044083550030_dp,  0.8656312023878318_dp, &
        & 0.9445750230732326_dp,  0.9894009349916499_dp]
    real(dp), parameter :: gl_w(16) = [&
        &0.0271524594117541_dp, 0.0622535239386479_dp, &
        &0.0951585116824928_dp, 0.1246289712555339_dp, &
        &0.1495959888165767_dp, 0.1691565193950025_dp, &
        &0.1826034150449236_dp, 0.1894506104550685_dp, &
        &0.1894506104550685_dp, 0.1826034150449236_dp, &
        &0.1691565193950025_dp, 0.1495959888165767_dp, &
        &0.1246289712555339_dp, 0.0951585116824928_dp, &
        &0.0622535239386479_dp, 0.0271524594117541_dp]

contains
    !> Collect the tests.
    subroutine collect_vac_kernels(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("asymptote", test_asymptote), &
            &new_unittest("G_regular", test_G_regular), &
            &new_unittest("H_regular_numderiv", test_H_regular_numderiv), &
            &new_unittest("GH_singular_quadrature", &
            &   test_GH_singular_quadrature), &
            &new_unittest("GH_edge_correction", test_GH_edge_correction), &
            &new_unittest("GH_int_1_regular", test_GH_int_1_regular) &
            &]
    end subroutine collect_vac_kernels

    !--------------------------------------------------------------------------
    ! geometry and kernel helpers
    !--------------------------------------------------------------------------

    !> Position, normal and its poloidal derivative on the circular boundary.
    pure subroutine circ(t,x,norm,dnorm)
        real(dp), intent(in) :: t                                               ! poloidal angle
        real(dp), intent(out) :: x(2)                                           ! (R,Z)
        real(dp), intent(out), optional :: norm(2)                              ! -R (Z_t, -R_t)
        real(dp), intent(out), optional :: dnorm(2)                             ! d norm / dt

        real(dp) :: R                                                           ! major radius of point

        R = R_0 + a_min*cos(t)
        x = [R, a_min*sin(t)]
        if (present(norm)) norm = [-R*a_min*cos(t), -R*a_min*sin(t)]
        if (present(dnorm)) dnorm = &
            &[a_min*sin(t)*(a_min*cos(t)+R), a_min*(a_min*sin(t)**2-R*cos(t))]
    end subroutine circ

    !> Argument gamma = 1 + rho^2/(2 R_s R_in) of the toroidal harmonics.
    pure function gam_arg(x_s,x_in) result(gam)
        real(dp), intent(in) :: x_s(2), x_in(2)
        real(dp) :: gam

        gam = 1._dp + sum((x_s-x_in)**2)/(2._dp*x_s(1)*x_in(1))
    end function gam_arg

    !> Toroidal harmonics Q_{n-3/2}(gam) and Q_{n-1/2}(gam) through dtorh.
    subroutine q_pair(gam,n,q,ierr)
        real(dp), intent(in) :: gam                                             ! argument
        integer, intent(in) :: n                                                ! toroidal mode number (>0)
        real(dp), intent(out) :: q(2)                                           ! Q_{n-3/2}, Q_{n-1/2}
        integer, intent(out) :: ierr                                            ! error status

        real(dp) :: pl_loc(0:n), ql_loc(0:n)                                    ! toroidal harmonics
        integer :: newn                                                         ! maximum reached degree

        ierr = dtorh1(gam,0,n,pl_loc,ql_loc,newn)
        if (ierr.eq.0 .and. newn.lt.n) ierr = 1
        q = ql_loc(n-1:n)
    end subroutine q_pair

    !> Regular-point helper variable Aij.
    !!
    !! This is the expression set up by the caller (vac_ops.calc_GH_2) at
    !! regular points; it is validated jointly with the H kernel formula by
    !! test_H_regular_numderiv.
    pure function aij_reg(x_s,norm_s,x_in,n) result(Aij)
        real(dp), intent(in) :: x_s(2), norm_s(2), x_in(2)
        integer, intent(in) :: n
        real(dp) :: Aij

        real(dp) :: rho2                                                        ! squared distance

        rho2 = sum((x_s-x_in)**2)
        Aij = 2._dp*x_in(1)*x_s(1)/(4._dp*x_in(1)*x_s(1)+rho2)*(n-0.5_dp)*&
            &(-2._dp/rho2*sum(norm_s*(x_s-x_in)) + norm_s(1)/x_s(1))
    end function aij_reg

    !> Singular-point helper variable Aij (limit of aij_reg for source ->
    !! influence point), as set up by the caller at (near-)singular points.
    pure function aij_sing(x_in,norm_in,dnorm_in,n) result(Aij)
        real(dp), intent(in) :: x_in(2), norm_in(2), dnorm_in(2)
        integer, intent(in) :: n
        real(dp) :: Aij

        Aij = 0.5_dp*(n-0.5_dp)*(x_in(1)*&
            &(norm_in(1)*dnorm_in(2)-norm_in(2)*dnorm_in(1))/&
            &sum(norm_in**2) + norm_in(1)/x_in(1))
    end function aij_sing

    !> b_n = sum_{k=1}^n 2/(2k-1).
    pure function b_of(n) result(b)
        integer, intent(in) :: n
        real(dp) :: b

        integer :: kd

        b = 0._dp
        do kd = 1,n
            b = b + 2._dp/(2*kd-1)
        end do
    end function b_of

    !> The G kernel  k_G(t) = -2 Q_{n-1/2}(gam(t)) / sqrt(R_in R_s(t)),  such
    !! that  G_ij = int hat_i(t) k_G(t) dt  over the source subintervals.
    subroutine k_G(t,x_in,n,res,ierr)
        real(dp), intent(in) :: t                                               ! source angle
        real(dp), intent(in) :: x_in(2)                                         ! influence point
        integer, intent(in) :: n                                                ! toroidal mode number
        real(dp), intent(out) :: res                                            ! kernel value
        integer, intent(out) :: ierr                                            ! error status

        real(dp) :: x_s(2)                                                      ! source point
        real(dp) :: q(2)                                                        ! toroidal harmonics

        call circ(t,x_s)
        call q_pair(gam_arg(x_s,x_in),n,q,ierr)
        res = -2._dp*q(2)/sqrt(x_in(1)*x_s(1))
    end subroutine k_G

    !> The H kernel  k_H(t) = 2/sqrt(R_in R_s) *
    !!  [ -norm_R/(2 R_s) Q_{n-1/2} - Aij (gam Q_{n-1/2} - Q_{n-3/2}) ],
    !! i.e. the directional derivative of the G kernel along the source
    !! normal (validated independently in test_H_regular_numderiv).
    subroutine k_H(t,x_in,n,res,ierr)
        real(dp), intent(in) :: t                                               ! source angle
        real(dp), intent(in) :: x_in(2)                                         ! influence point
        integer, intent(in) :: n                                                ! toroidal mode number
        real(dp), intent(out) :: res                                            ! kernel value
        integer, intent(out) :: ierr                                            ! error status

        real(dp) :: x_s(2), norm_s(2)                                           ! source point and normal
        real(dp) :: q(2)                                                        ! toroidal harmonics
        real(dp) :: gam                                                         ! argument
        real(dp) :: Aij                                                         ! helper variable

        call circ(t,x_s,norm=norm_s)
        gam = gam_arg(x_s,x_in)
        call q_pair(gam,n,q,ierr)
        Aij = aij_reg(x_s,norm_s,x_in,n)
        res = 2._dp/sqrt(x_in(1)*x_s(1)) * &
            &(-norm_s(1)/(2._dp*x_s(1))*q(2) - Aij*(gam*q(2)-q(1)))
    end subroutine k_H

    !--------------------------------------------------------------------------
    ! tests
    !--------------------------------------------------------------------------

    !> Verify the asymptote  Q_{n-1/2}(1+x) = -1/2 ln(x/32) - b_n + O(x ln x)
    !! against dtorh. The same asymptote (in an equivalent chord-based form)
    !! underlies the analytical integration of the singular subintervals in
    !! calc_GH_int_2, and it is used below to subtract the log singularity
    !! from the brute-force quadratures.
    subroutine test_asymptote(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: kd, n, ierr                                                  ! counters and error status
        real(dp) :: x, q(2), q_mod                                              ! argument, harmonics, model
        real(dp) :: max_diff                                                    ! worst absolute deviation

        do n = 1,3
            max_diff = 0._dp
            do kd = 4,8                                                         ! x = 1e-4 .. 1e-8
                x = 10._dp**(-kd)
                call q_pair(1._dp+x,n,q,ierr)
                call check(error, ierr, 0, 'dtorh1 failed')
                if (allocated(error)) return
                q_mod = -0.5_dp*log(x/32._dp) - b_of(n)
                max_diff = max(max_diff,abs(q(2)-q_mod))
            end do
            ! theory: error O(x ln x) ~ 1e-3 at x = 1e-4
            call check(error, max_diff.lt.2.e-3_dp, &
                &'asymptote deviates by '//trim(r2str(max_diff))//&
                &' for n = '//trim(i2str(n)))
            if (allocated(error)) return
        end do
    end subroutine test_asymptote

    !> Regular subinterval: G against direct evaluation of the kernel at the
    !! interval ends (the implementation is the trapezoidal rule, so the
    !! comparison is exact up to roundoff).
    subroutine test_G_regular(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2                                             ! toroidal mode number
        real(dp), parameter :: t_s_in(2) = [0.30_dp,0.35_dp]                    ! source interval
        real(dp), parameter :: t_i = 2.0_dp                                     ! influence angle

        integer :: kd, ierr                                                     ! counter, error status
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(2,2), norm_s(2,2)                                       ! source geometry
        real(dp) :: x_in(2), norm_in(2)                                         ! influence geometry
        real(dp) :: Aij(2), ql(2,2), q(2)                                       ! helpers
        real(dp) :: G_ref                                                       ! reference

        do kd = 1,2
            call circ(t_s_in(kd),x_s(kd,:),norm=norm_s(kd,:))
        end do
        call circ(t_i,x_in,norm=norm_in)
        do kd = 1,2
            call q_pair(gam_arg(x_s(kd,:),x_in),n,ql(kd,:),ierr)
            call check(error, ierr, 0, 'dtorh1 failed')
            if (allocated(error)) return
            Aij(kd) = aij_reg(x_s(kd,:),norm_s(kd,:),x_in,n)
        end do

        G = 0._dp
        H = 0._dp
        call calc_GH_int_2(G,H,t_s_in,t_i,x_s,x_in,norm_s,norm_in,Aij,ql,&
            &1.e-30_dp,[b_of(n-1),b_of(n)],n)

        do kd = 1,2
            call q_pair(gam_arg(x_s(kd,:),x_in),n,q,ierr)
            G_ref = -(t_s_in(2)-t_s_in(1))*q(2)/sqrt(x_in(1)*x_s(kd,1))
            call check(error, G(kd), G_ref, thr=1.e-13_dp, rel=.true., &
                &message='G('//trim(i2str(kd))//') mismatch')
            if (allocated(error)) return
        end do
    end subroutine test_G_regular

    !> Regular subinterval: H against a numerical directional derivative of
    !! the G kernel along the (unnormalized) source normal:
    !!   H(kd) = dt * d/deps [ Q_{n-1/2}(gam)/sqrt(R_in R_s) ]
    !!                       (x_s -> x_s + eps norm_s),
    !! which independently validates the H kernel algebra including Aij.
    subroutine test_H_regular_numderiv(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2                                             ! toroidal mode number
        real(dp), parameter :: t_s_in(2) = [0.30_dp,0.35_dp]                    ! source interval
        real(dp), parameter :: t_i = 2.0_dp                                     ! influence angle
        real(dp), parameter :: eps = 1.e-7_dp                                   ! step for central difference

        integer :: kd, ierr                                                     ! counter, error status
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(2,2), norm_s(2,2)                                       ! source geometry
        real(dp) :: x_in(2), norm_in(2)                                         ! influence geometry
        real(dp) :: Aij(2), ql(2,2)                                             ! helpers
        real(dp) :: x_pert(2), q(2)                                             ! perturbed point, harmonics
        real(dp) :: fp, fm, H_ref                                               ! function values, reference

        do kd = 1,2
            call circ(t_s_in(kd),x_s(kd,:),norm=norm_s(kd,:))
        end do
        call circ(t_i,x_in,norm=norm_in)
        do kd = 1,2
            call q_pair(gam_arg(x_s(kd,:),x_in),n,ql(kd,:),ierr)
            call check(error, ierr, 0, 'dtorh1 failed')
            if (allocated(error)) return
            Aij(kd) = aij_reg(x_s(kd,:),norm_s(kd,:),x_in,n)
        end do

        G = 0._dp
        H = 0._dp
        call calc_GH_int_2(G,H,t_s_in,t_i,x_s,x_in,norm_s,norm_in,Aij,ql,&
            &1.e-30_dp,[b_of(n-1),b_of(n)],n)

        do kd = 1,2
            x_pert = x_s(kd,:) + eps*norm_s(kd,:)
            call q_pair(gam_arg(x_pert,x_in),n,q,ierr)
            call check(error, ierr, 0, 'dtorh1 failed')
            if (allocated(error)) return
            fp = q(2)/sqrt(x_in(1)*x_pert(1))
            x_pert = x_s(kd,:) - eps*norm_s(kd,:)
            call q_pair(gam_arg(x_pert,x_in),n,q,ierr)
            fm = q(2)/sqrt(x_in(1)*x_pert(1))
            H_ref = (t_s_in(2)-t_s_in(1))*(fp-fm)/(2._dp*eps)
            call check(error, H(kd), H_ref, thr=1.e-6_dp, rel=.true., &
                &message='H('//trim(i2str(kd))//') vs numerical derivative: '&
                &//trim(r2str(H(kd)))//' vs '//trim(r2str(H_ref)))
            if (allocated(error)) return
        end do
    end subroutine test_H_regular_numderiv

    !> Brute-force reference for the integrals over a subinterval with the
    !! influence point at the left end (t_in = t_s(1)):
    !!   I(kd) = int hat_kd(t) k(t) dt,
    !! by subtracting the exactly integrable logarithmic model
    !!   s_kd(t) = c hat_kd(t) (ln(a_r (t-t_1)) + b-type constants)
    !! and integrating the smooth remainder with panel-refined Gauss rules.
    !!
    !! Shared by the singular and edge-correction tests.
    subroutine brute_force_GH_sing(t_s_in,n,I_G,I_H,ierr)
        real(dp), intent(in) :: t_s_in(2)                                       ! source interval
        integer, intent(in) :: n                                                ! toroidal mode number
        real(dp), intent(out) :: I_G(2), I_H(2)                                 ! integrals
        integer, intent(out) :: ierr                                            ! error status

        integer :: kd, id, jd                                                   ! counters
        real(dp) :: x_in(2), norm_in(2), dnorm_in(2)                            ! influence geometry (= left end)
        real(dp) :: x_s2(2)                                                     ! right-end position
        real(dp) :: delta                                                       ! interval length
        real(dp) :: a_r                                                         ! chord coefficient in the log model
        real(dp) :: c_G, c_H_log, c_H_const                                     ! model coefficients
        real(dp) :: I_ln(2)                                                     ! int hat_kd ln(u) du
        real(dp) :: I_hat                                                       ! int hat_kd du
        real(dp) :: u_lo                                                        ! truncation of the remainder integral
        real(dp) :: p_lo, p_hi                                                  ! panel ends
        real(dp) :: t_q, u_q, w_q                                               ! quadrature node in t and u, weight
        real(dp) :: kv, sv(2)                                                   ! kernel and model values
        real(dp) :: aij0                                                        ! Aij at the singular point
        integer, parameter :: n_pan = 60                                        ! quadrature panels
        real(dp), parameter :: pan_fac = 0.7_dp                                 ! geometric panel refinement

        ierr = 0
        delta = t_s_in(2)-t_s_in(1)
        call circ(t_s_in(1),x_in,norm=norm_in,dnorm=dnorm_in)
        call circ(t_s_in(2),x_s2)
        aij0 = aij_sing(x_in,norm_in,dnorm_in,n)

        ! log model: Q_{n-1/2}(gam(t)) ~ -ln(a_r (t-t_1)) - b_n with
        ! a_r = |x'(t_1)|/(8 R_in) (from the asymptote Q = -1/2 ln(x/32) - b_n
        ! with x = gam-1 ~ (|x'| u)^2/(2 R_in^2))
        a_r = a_min/(8._dp*x_in(1))

        ! exact log moments: int_0^delta hat_kd(u) ln(u) du
        I_ln(1) = delta*(0.5_dp*log(delta)-0.75_dp)
        I_ln(2) = delta*(0.5_dp*log(delta)-0.25_dp)
        I_hat = 0.5_dp*delta

        ! G model: k_G ~ c_G (-ln(a_r u) - b_n), c_G = -2/sqrt(R_in R_1)
        c_G = -2._dp/x_in(1)                                                    ! R_1 = R_in
        ! H model: k_H ~ c_H_log (-ln(a_r u) - b_n) + c_H_const with
        ! c_H_log = -norm_R(t_1)/R_1 / R_in and the finite part of the
        ! -Aij (gam Q_{n-1/2} - Q_{n-3/2}) term: (gam Q - Q') -> -2/(2n-1)
        c_H_log = -norm_in(1)/(x_in(1)**2)
        c_H_const = 2._dp/x_in(1)*aij0*2._dp/(2._dp*n-1._dp)

        do kd = 1,2
            I_G(kd) = c_G*(-I_ln(kd)-(log(a_r)+b_of(n))*I_hat)
            I_H(kd) = c_H_log*(-I_ln(kd)-(log(a_r)+b_of(n))*I_hat) + &
                &c_H_const*I_hat
        end do

        ! truncation: gamma-1 >= ~3e-10 for dtorh reliability
        u_lo = sqrt(3.e-10_dp*2._dp*x_in(1)*x_s2(1))/a_min

        ! panel-refined quadrature of the smooth remainders on [u_lo,delta]
        p_hi = delta
        do id = 1,n_pan
            p_lo = max(u_lo,p_hi*pan_fac)
            if (id.eq.n_pan) p_lo = u_lo
            do jd = 1,16
                u_q = 0.5_dp*(p_lo+p_hi) + 0.5_dp*(p_hi-p_lo)*gl_x(jd)
                w_q = 0.5_dp*(p_hi-p_lo)*gl_w(jd)
                t_q = t_s_in(1) + u_q
                ! model values
                sv(1) = c_G*(-log(a_r*u_q)-b_of(n))
                sv(2) = c_H_log*(-log(a_r*u_q)-b_of(n)) + c_H_const
                ! true kernels minus models, weighted by the hat functions
                call k_G(t_q,x_in,n,kv,ierr)
                if (ierr.ne.0) return
                I_G(1) = I_G(1) + w_q*(1._dp-u_q/delta)*(kv-sv(1))
                I_G(2) = I_G(2) + w_q*(u_q/delta)*(kv-sv(1))
                call k_H(t_q,x_in,n,kv,ierr)
                if (ierr.ne.0) return
                I_H(1) = I_H(1) + w_q*(1._dp-u_q/delta)*(kv-sv(2))
                I_H(2) = I_H(2) + w_q*(u_q/delta)*(kv-sv(2))
            end do
            p_hi = p_lo
            if (p_hi.le.u_lo) exit
        end do
    end subroutine brute_force_GH_sing

    !> Singular subinterval (influence point = left end, both ends within the
    !! analytical-approximation region): the analytical integrals of
    !! calc_GH_int_2 against brute-force quadrature, at two interval sizes to
    !! also confirm convergence of the approximation.
    subroutine test_GH_singular_quadrature(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2                                             ! toroidal mode number
        real(dp), parameter :: t_1 = 0.4_dp                                     ! left end = influence point
        real(dp), parameter :: deltas(2) = [5.e-3_dp,5.e-4_dp]                  ! interval sizes
        real(dp), parameter :: tols(2) = [2.e-2_dp,2.e-3_dp]                    ! relative tolerances

        integer :: id, kd, ierr                                                 ! counters, error status
        real(dp) :: t_s_in(2)                                                   ! source interval
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(2,2), norm_s(2,2)                                       ! source geometry
        real(dp) :: x_in(2), norm_in(2), dnorm_in(2)                            ! influence geometry
        real(dp) :: Aij(2), ql(2,2)                                             ! helpers
        real(dp) :: rho2_2                                                      ! squared distance to right end
        real(dp) :: I_G(2), I_H(2)                                              ! brute-force references

        do id = 1,size(deltas)
            t_s_in = [t_1,t_1+deltas(id)]
            do kd = 1,2
                call circ(t_s_in(kd),x_s(kd,:),norm=norm_s(kd,:))
            end do
            call circ(t_1,x_in,norm=norm_in,dnorm=dnorm_in)
            rho2_2 = sum((x_s(2,:)-x_in)**2)

            ! at singular points the caller supplies the limit form of Aij
            ! and does not evaluate the toroidal harmonics
            Aij = aij_sing(x_in,norm_in,dnorm_in,n)
            ql = 0._dp

            G = 0._dp
            H = 0._dp
            call calc_GH_int_2(G,H,t_s_in,t_1,x_s,x_in,norm_s,norm_in,Aij,ql,&
                &4._dp*rho2_2,[b_of(n-1),b_of(n)],n)                            ! tol > rho2_2: both ends in analytical region

            call brute_force_GH_sing(t_s_in,n,I_G,I_H,ierr)
            call check(error, ierr, 0, 'brute-force quadrature failed')
            if (allocated(error)) return

            do kd = 1,2
                call check(error, G(kd), I_G(kd), &
                    &thr=tols(id)*abs(I_G(kd)), &
                    &message='singular G('//trim(i2str(kd))//'), delta = '//&
                    &trim(r2str(deltas(id)))//': '//trim(r2str(G(kd)))//&
                    &' vs '//trim(r2str(I_G(kd))))
                if (allocated(error)) return
                call check(error, H(kd), I_H(kd), &
                    &thr=tols(id)*abs(I_H(kd)), &
                    &message='singular H('//trim(i2str(kd))//'), delta = '//&
                    &trim(r2str(deltas(id)))//': '//trim(r2str(H(kd)))//&
                    &' vs '//trim(r2str(I_H(kd))))
                if (allocated(error)) return
            end do
        end do
    end subroutine test_GH_singular_quadrature

    !> Near-singular subinterval with only the left end inside the
    !! analytical-approximation region (tol between the two rho^2): the
    !! implementation adds a trapezoidal correction for the difference
    !! between the true toroidal function and its approximation at the far
    !! end, which needs the toroidal harmonics there.
    subroutine test_GH_edge_correction(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2                                             ! toroidal mode number
        real(dp), parameter :: t_1 = 0.4_dp                                     ! left end = influence point
        real(dp), parameter :: delta = 5.e-3_dp                                 ! interval size
        real(dp), parameter :: tol_rel = 2.e-2_dp                               ! relative tolerance

        integer :: kd, ierr                                                     ! counter, error status
        real(dp) :: t_s_in(2)                                                   ! source interval
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(2,2), norm_s(2,2)                                       ! source geometry
        real(dp) :: x_in(2), norm_in(2), dnorm_in(2)                            ! influence geometry
        real(dp) :: Aij(2), ql(2,2)                                             ! helpers
        real(dp) :: rho2_2                                                      ! squared distance to right end
        real(dp) :: I_G(2), I_H(2)                                              ! brute-force references

        t_s_in = [t_1,t_1+delta]
        do kd = 1,2
            call circ(t_s_in(kd),x_s(kd,:),norm=norm_s(kd,:))
        end do
        call circ(t_1,x_in,norm=norm_in,dnorm=dnorm_in)
        rho2_2 = sum((x_s(2,:)-x_in)**2)

        ! left end singular: limit Aij; right end regular: caller supplies
        ! the toroidal harmonics and the regular Aij there
        Aij(1) = aij_sing(x_in,norm_in,dnorm_in,n)
        Aij(2) = aij_reg(x_s(2,:),norm_s(2,:),x_in,n)
        ql = 0._dp
        call q_pair(gam_arg(x_s(2,:),x_in),n,ql(2,:),ierr)
        call check(error, ierr, 0, 'dtorh1 failed')
        if (allocated(error)) return

        G = 0._dp
        H = 0._dp
        call calc_GH_int_2(G,H,t_s_in,t_1,x_s,x_in,norm_s,norm_in,Aij,ql,&
            &0.25_dp*rho2_2,[b_of(n-1),b_of(n)],n)                              ! tol < rho2_2: right end outside analytical region

        call brute_force_GH_sing(t_s_in,n,I_G,I_H,ierr)
        call check(error, ierr, 0, 'brute-force quadrature failed')
        if (allocated(error)) return

        do kd = 1,2
            call check(error, G(kd), I_G(kd), thr=tol_rel*abs(I_G(kd)), &
                &message='edge-corrected G('//trim(i2str(kd))//'): '//&
                &trim(r2str(G(kd)))//' vs '//trim(r2str(I_G(kd))))
            if (allocated(error)) return
            call check(error, H(kd), I_H(kd), thr=tol_rel*abs(I_H(kd)), &
                &message='edge-corrected H('//trim(i2str(kd))//'): '//&
                &trim(r2str(H(kd)))//' vs '//trim(r2str(I_H(kd))))
            if (allocated(error)) return
        end do
    end subroutine test_GH_edge_correction

    !> Regular subinterval of the 3-D (field-line aligned) kernels: G is the
    !! free-space Green's function -1/|x_s - x_in| and H its directional
    !! derivative along the source normal, checked against a central
    !! difference.
    subroutine test_GH_int_1_regular(error)
        type(error_type), allocatable, intent(out) :: error

        real(dp), parameter :: eps = 1.e-7_dp                                   ! step for central difference

        integer :: kd                                                           ! counter
        real(dp) :: G(2), H(2)                                                  ! kernel results
        real(dp) :: x_s(2,3), norm_s(2,3)                                       ! source geometry
        real(dp) :: x_in(3)                                                     ! influence point
        real(dp) :: x_pert(3)                                                   ! perturbed point
        real(dp) :: fp, fm, r, H_ref                                            ! function values, distance, reference

        x_s(1,:) = [2.9_dp,0.1_dp,0.2_dp]
        x_s(2,:) = [2.85_dp,0.15_dp,0.25_dp]
        norm_s(1,:) = [0.9_dp,0.1_dp,0.3_dp]
        norm_s(2,:) = [0.8_dp,0.2_dp,0.4_dp]
        x_in = [1.5_dp,-1.0_dp,0.5_dp]

        G = 0._dp
        H = 0._dp
        call calc_GH_int_1(G,H,x_s,x_in,norm_s,[1._dp,0.1_dp,1._dp,0.5_dp],&
            &[0.1_dp,0.1_dp],1.e-30_dp)

        do kd = 1,2
            r = sqrt(sum((x_s(kd,:)-x_in)**2))
            call check(error, G(kd), -1._dp/r, thr=1.e-14_dp, rel=.true., &
                &message='G('//trim(i2str(kd))//') mismatch')
            if (allocated(error)) return

            ! H = norm . grad_s (-1/r), via central difference
            x_pert = x_s(kd,:) + eps*norm_s(kd,:)
            fp = -1._dp/sqrt(sum((x_pert-x_in)**2))
            x_pert = x_s(kd,:) - eps*norm_s(kd,:)
            fm = -1._dp/sqrt(sum((x_pert-x_in)**2))
            H_ref = (fp-fm)/(2._dp*eps)
            call check(error, H(kd), H_ref, thr=1.e-6_dp, rel=.true., &
                &message='H('//trim(i2str(kd))//') vs numerical derivative: '&
                &//trim(r2str(H(kd)))//' vs '//trim(r2str(H_ref)))
            if (allocated(error)) return
        end do
    end subroutine test_GH_int_1_regular
end module test_vac_kernels
