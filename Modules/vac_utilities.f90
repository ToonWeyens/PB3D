!------------------------------------------------------------------------------!
!> Numerical utilities related to the vacuum response.
!!
!! \see vac_ops for information.
!------------------------------------------------------------------------------!
module vac_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use vac_vars, only: vac_type, BLACSCTXTSIZE, loc_idx

    implicit none
    
    ! 16-point Gauss-Legendre rule on [-1, 1], for the cell integral of H_sing_node
    integer, parameter :: n_gl16 = 16                                           !< number of points
    real(dp), parameter :: gl16_x(n_gl16) = [&
        &-0.9894009349916499_dp, -0.9445750230732326_dp, &
        &-0.8656312023878316_dp, -0.7554044083550031_dp, &
        &-0.6178762444026438_dp, -0.4580167776572274_dp, &
        &-0.2816035507792589_dp, -0.0950125098376374_dp, &
        & 0.0950125098376374_dp,  0.2816035507792589_dp, &
        & 0.4580167776572274_dp,  0.6178762444026438_dp, &
        & 0.7554044083550031_dp,  0.8656312023878316_dp, &
        & 0.9445750230732326_dp,  0.9894009349916499_dp]                        !< abscissas
    real(dp), parameter :: gl16_w(n_gl16) = [&
        &0.0271524594117541_dp, 0.0622535239386476_dp, &
        &0.0951585116824930_dp, 0.1246289712555339_dp, &
        &0.1495959888165767_dp, 0.1691565193950026_dp, &
        &0.1826034150449236_dp, 0.1894506104550686_dp, &
        &0.1894506104550686_dp, 0.1826034150449236_dp, &
        &0.1691565193950026_dp, 0.1495959888165767_dp, &
        &0.1246289712555339_dp, 0.0951585116824930_dp, &
        &0.0622535239386476_dp, 0.0271524594117541_dp]                          !< weights
    private
    public calc_GH_int_1, calc_GH_int_2, vec_dis2loc, mat_dis2loc, &
        &n_col1_loc, vec_glob2dis, dis_matvec, &
        &interlaced_vac_copy, vac_steps, solid_angle_node, H_sing_node, &
        &calc_diag_direct, dis_diag, cross3
    
contains
    !> Calculate  G_ij  and   H_ij  on  an  interval  for   field  line  aligned
    !! configurations.
    !!
    !! The source variables  are given on a stencil of  four consecutive points
    !! of the field line, <tt>[1:4,dim]</tt>  with \c dim the Cartesian dimension;
    !! the interval  is the one between  stencil points \c i_int  and \c i_int+1.
    !! Each endpoint  receives the kernel  integrated over  its half  of the
    !! interval,  normalized to  the  half-cell weight  \f$\Delta p  \Delta a /
    !! 2\f$ that  calc_GH_1 applies, so that  a constant kernel  returns its
    !! value.
    !!
    !! Far from the influence point, more  than \c n_near interval lengths away,
    !! the  kernel  value  at  the  endpoint  is used,  which  amounts  to  the
    !! trapezoidal  rule along  the field  line. Closer,  the kernel  is
    !! integrated  over  the  half-interval  with 8-point  Gauss-Legendre  on
    !! the cubic  interpolant of position  and normal through  the stencil.
    !! This resolves  the variation of  the kernel along  a neighbouring field
    !! line: where the field lines run almost toroidally, as at the bean tips
    !! of W7-X,  they are packed  far closer  (0.007 \f$R_0\f$) than  the points
    !! along them (0.06 \f$R_0\f$), the  kernel seen along the neighbouring line
    !! is a peak much narrower than the point spacing, and sampling it at the
    !! points  gives  errors  of order  one  whose  sign depends  on how  the
    !! points of neighbouring lines  happen to line up (issue #4,  31 October
    !! 2026: the fractional part of \f$q \Delta\theta / \Delta\alpha\f$ sorted
    !! every W7-X grid into far too large and plausible responses).
    !!
    !! For  subintegrals adjacent  to the  singularity (i.e.  with the  influence
    !! point as  one of  their endpoints), the  contribution of  the singular
    !! endpoint is  replaced by the analytical  integral of the kernel  over its
    !! grid half-cell.  In field-aligned coordinates,  with \f$\delta p\f$ along
    !! the field  line (metric \f$g_{\theta\theta}\f$  = <tt>h_fac(3)</tt>)  and
    !! \f$\delta a\f$  the field  line label  (metric \f$g_{\alpha\alpha}\f$  =
    !! <tt>h_fac(1)</tt>, coupling \f$g_{\alpha\theta}\f$ = <tt>h_fac(2)</tt>),
    !! the squared distance near the singular point is the quadratic form
    !! \f$ g_{\theta\theta} \delta p^2 + 2 g_{\alpha\theta} \delta p \delta a
    !!  + g_{\alpha\alpha} \delta a^2 \f$,
    !! and with  the metric-normalized  half-cell sizes  \f$U = \sqrt{g_{\theta
    !! \theta}} \Delta p/2\f$, \f$V = \sqrt{g_{\alpha\alpha}} \Delta a/2\f$ and
    !! nonorthogonality  \f$s  = g_{\alpha\theta}  /  \sqrt{g_{\alpha\alpha}
    !! g_{\theta\theta}}\f$, the exact primitive is
    !! \f[ \int_0^U \int_{-V}^{V} \frac{\text{d}u \text{d}v}
    !!  {\sqrt{u^2+2suv+v^2}} =
    !!  U \ln\left(\frac{V+sU+D_+}{-V+sU+D_-}\right) +
    !!  V \ln\left(\frac{U+sV+D_+}{-U+sV+D_-}\right), \quad
    !!  D_\pm = \sqrt{U^2 \pm 2sUV + V^2}. \f]
    !! In the  same units  as the  regular kernel  values, the singular endpoint
    !! contribution is then
    !! \f$-2 F / (\Delta p \Delta a \sqrt{g_{\alpha\alpha} g_{\theta\theta}})
    !! \f$.
    !!
    !! (The previous  expression had  the same  logarithm twice,  diverged for
    !! orthogonal  coordinates  and  had  the  wrong  sign;  it  is  verified
    !! against brute-force quadrature in the full-stack tests.)
    !!
    !! \note This routine does not calculate the contribution \f$2\beta\f$ nor
    !! the  cell integral  of the  H kernel  at the  singular point:  calc_GH_1
    !! sets the  diagonal of H  afterwards, from  the row-sum identity  or from
    !! the direct value of calc_diag_direct.
    !!
    !! \note With \c use_near false the near-field rule is not used. For
    !! points off the boundary (vac_ops.eval_pot_ext()) it pays only within
    !! about an interval length of the boundary: farther out the trapezoidal
    !! rule along the lines and across them is more accurate, because the
    !! errors of the trapezoidal rule at the ends of the lines cancel between
    !! the lines, whose ends interleave, while the edge of the near-field zone
    !! leaves an error of second order in the step along the lines (on the
    !! circular VMEC twin at 64 lines of 100 points: 1e-3 of the potential at
    !! half a minor radius from the boundary, against 1e-7 with the
    !! trapezoidal rule alone; issue #4, 6 November 2026).
    subroutine calc_GH_int_1(G,H,x_s,x_in,norm_s,h_fac_in,step_size,tol,i_int,&
        &use_near)
        ! input / output
        real(dp), intent(inout) :: G(2)                                         !< G
        real(dp), intent(inout) :: H(2)                                         !< H
        real(dp), intent(in) :: x_s(4,3)                                        !< position vector at the four stencil points
        real(dp), intent(in) :: x_in(3)                                         !< position vector at which to calculate influence
        real(dp), intent(in) :: norm_s(4,3)                                     !< normal vector at the four stencil points
        real(dp), intent(in) :: h_fac_in(3)                                     !< metric factors g_aa, g_at, g_tt at the influence point
        real(dp), intent(in) :: step_size(2)                                    !< step sizes in parallel direction and alpha
        real(dp), intent(in) :: tol                                             !< tolerance on distance between points
        integer, intent(in) :: i_int                                            !< stencil index of the left point of the interval (1, 2 or 3)
        logical, intent(in), optional :: use_near                               !< whether the near-field rule may be used (default true)

        ! local variables
        integer, parameter :: n_gl = 8                                          ! Gauss-Legendre points per half-interval
        real(dp), parameter :: n_near = 4._dp                                   ! near-field range in interval lengths
        real(dp), parameter :: gl_x(n_gl) = [&
            &-0.9602898564975363_dp, -0.7966664774136267_dp, &
            &-0.5255324099163290_dp, -0.1834346424956498_dp, &
            & 0.1834346424956498_dp,  0.5255324099163290_dp, &
            & 0.7966664774136267_dp,  0.9602898564975363_dp]                    ! Gauss-Legendre abscissas on [-1,1]
        real(dp), parameter :: gl_w(n_gl) = [&
            &0.1012285362903763_dp, 0.2223810344533745_dp, &
            &0.3137066458778873_dp, 0.3626837833783620_dp, &
            &0.3626837833783620_dp, 0.3137066458778873_dp, &
            &0.2223810344533745_dp, 0.1012285362903763_dp]                      ! Gauss-Legendre weights on [-1,1]
        integer :: kd, jd                                                       ! counters
        real(dp) :: r2(2)                                                       ! distance between source and influence points
        real(dp) :: s_no                                                        ! nonorthogonality g_at/sqrt(g_aa g_tt)
        real(dp) :: U_hc, V_hc                                                  ! metric-normalized half-cell sizes
        real(dp) :: D_pm(2)                                                     ! D_+ and D_-
        real(dp) :: F_int                                                       ! primitive of the half-cell integral
        real(dp) :: s, l(4)                                                     ! stencil parameter and cubic Lagrange weights
        real(dp) :: d_q(3), n_q(3), r_q                                         ! relative position, normal and distance at a quadrature point
        logical :: near                                                         ! interval within n_near lengths of the influence point

        ! initialize
        G = 0._dp
        H = 0._dp

        ! calculate distance between source points
        do kd = 1,2
            r2(kd) = sum((x_s(i_int+kd-1,:)-x_in)**2)
        end do
        near = minval(r2).lt.n_near**2*sum((x_s(i_int+1,:)-x_s(i_int,:))**2)
        if (present(use_near)) near = near .and. use_near

        do kd = 1,2
            if (r2(kd).le.tol) then                                             ! singular endpoint: analytical half-cell integral
                s_no = h_fac_in(2)/sqrt(h_fac_in(1)*h_fac_in(3))
                U_hc = sqrt(h_fac_in(3))*step_size(1)/2._dp
                V_hc = sqrt(h_fac_in(1))*step_size(2)/2._dp
                D_pm(1) = sqrt(U_hc**2+2*s_no*U_hc*V_hc+V_hc**2)
                D_pm(2) = sqrt(U_hc**2-2*s_no*U_hc*V_hc+V_hc**2)
                F_int = U_hc*log(abs((V_hc+s_no*U_hc+D_pm(1))/&
                    &(-V_hc+s_no*U_hc+D_pm(2)))) + &
                    &V_hc*log(abs((U_hc+s_no*V_hc+D_pm(1))/&
                    &(-U_hc+s_no*V_hc+D_pm(2))))
                G(kd) = -2._dp*F_int/(step_size(1)*step_size(2)*&
                    &sqrt(h_fac_in(1)*h_fac_in(3)))
                H(kd) = 0._dp                                                   ! the diagonal of H is set by calc_GH_1
            else if (near) then                                                 ! near endpoint: kernel integrated over its half-interval
                do jd = 1,n_gl
                    s = i_int-1 + 0.25_dp*(2*kd-1+gl_x(jd))                     ! stencil parameter, nodes at 0, 1, 2, 3
                    l = [-(s-1)*(s-2)*(s-3)/6._dp, s*(s-2)*(s-3)/2._dp, &
                        &-s*(s-1)*(s-3)/2._dp, s*(s-1)*(s-2)/6._dp]
                    d_q = matmul(l,x_s) - x_in
                    n_q = matmul(l,norm_s)
                    r_q = sqrt(sum(d_q**2))
                    G(kd) = G(kd) - 0.5_dp*gl_w(jd)/r_q
                    H(kd) = H(kd) + 0.5_dp*gl_w(jd)*sum(n_q*d_q)/r_q**3
                end do
            else                                                                ! far endpoint: kernel value
                G(kd) = -1._dp/sqrt(r2(kd))
                ! dipole kernel norm.(x_s-x_in)/|x_s-x_in|^3, i.e. the
                ! directional derivative of G along the source normal
                H(kd) = sum(norm_s(i_int+kd-1,:)*(x_s(i_int+kd-1,:)-x_in))*&
                    &(-G(kd))**3
            end if
        end do
    end subroutine calc_GH_int_1
    
    !> Calculate G_ij and H_ij on an interval for axisymmetric configurations.
    !!
    !! The indices for the  source variables are <tt>[left:right,dim]</tt> where
    !! \c dim is the Cartesian dimension. The same holds for \c ql.
    !!
    !! For subintegrals close to the singularity (i.e. where the tolerance is no
    !! being  met), the  routine approximates  the toroidal  functions with  the
    !! analytical approximation.
    !! If  one of  the boundaries  of  the subintegral  is far  enough from  the
    !! singularity, there  is an additional  contribution due to  the difference
    !! between  the  actual  toroidal   function  and  the  approximation.  This
    !! contribution is zero  if both boundaries are close, because  in this case
    !! it  is assumed  that  the  entire integral  is  given  by the  analytical
    !! approximation.
    !!
    !! \note The kernels are those of eqs. (130)-(131) of the 2019 manuscript
    !! with  the sign  of H  reversed (the leading minus signs  of the  three H
    !! formulas): the  normal derivative is  taken along the  stored normal
    !! \f$-\mathcal{J}\nabla\psi\f$, out of the vacuum, for both vacuum styles
    !! (erratum E1 of that manuscript). This routine does not
    !! calculate the corner term: calc_GH_2 adds \f$-2\beta\f$ to the diagonal
    !! of H, with \f$\beta\f$ the vacuum's internal angle
    !! on the vacuum side.
    subroutine calc_GH_int_2(G,H,t_s,t_in,x_s,x_in,norm_s,norm_in,Aij,ql,tol,&
        &b_coeff,n_tor)
        ! input / output
        real(dp), intent(inout) :: G(2)                                         !< G
        real(dp), intent(inout) :: H(2)                                         !< H
        real(dp), intent(in) :: t_s(2)                                          !< pol. angle at left and right limit of source interval
        real(dp), intent(in) :: t_in                                            !< pol. angle at which to calculate influence
        real(dp), intent(in) :: x_s(2,2)                                        !< position vector at left and right limit of source interval
        real(dp), intent(in) :: x_in(2)                                         !< position vector at which to calculate influence
        real(dp), intent(in) :: norm_s(2,2)                                     !< normal vector at left and right limit of source interval
        real(dp), intent(in) :: norm_in(2)                                      !< normal vector at which to calculate influence
        real(dp), intent(in) :: Aij(2)                                          !< Aij helper variable for H
        real(dp), intent(in) :: ql(2,2)                                         !< ql for n-1 and n at left and right limit of source interval
        real(dp), intent(in) :: tol                                             !< tolerance on rho2
        real(dp), intent(in) :: b_coeff(2)                                      !< \f$b = \sum_{k=1}^n \frac{2}{2k-1}\f$ for n and n-1
        integer, intent(in) :: n_tor                                            !< toroidal mode number
        
        ! local variables
        integer :: kd, jd                                                       ! counters
        integer :: kd_bound                                                     ! counter for boundary
        real(dp) :: a_coeff                                                     ! coefficient a = |norm|/8R^2 at influence point
        real(dp) :: t_in_loc                                                    ! local t_in, possibly shifted by 2pi
        real(dp) :: delta_t(2)                                                  ! t_in - t_s
        real(dp) :: rho2(2)                                                     ! square of distance in projected poloidal plane
        real(dp) :: dlogd(2)                                                    ! delta log delta
        real(dp) :: n_loc(2)                                                    ! norm_s / |norm_s|
        real(dp) :: dn_loc(2)                                                   ! d n_loc / dtheta
        
        ! initialize
        G = 0._dp
        H = 0._dp
        do kd = 1,2
            rho2(kd) = sum((x_s(kd,:)-x_in)**2)
        end do
        if (pi .lt. sum(t_s)/2 - t_in) then                                     ! t_s closer to t_in + 2pi
            t_in_loc = t_in + 2*pi
        else if (-pi .gt. sum(t_s)/2 - t_in) then                               ! t_s closer to t_in - 2pi
            t_in_loc = t_in - 2*pi
        else                                                                    ! t_s closer to t_in
            t_in_loc = t_in
        end if
        n_loc = (norm_s(2,:)/sqrt(sum(norm_s(2,:)**2))+&
            &norm_s(1,:)/sqrt(sum(norm_s(1,:)**2))) / 2._dp
        dn_loc = (norm_s(2,:)/sqrt(sum(norm_s(2,:)**2))-&
            &norm_s(1,:)/sqrt(sum(norm_s(1,:)**2))) / (t_s(2)-t_s(1))
        
        ! check for (near-)singular point
        if (minval(rho2).le.tol) then
           ! There is at least one singular point: use analytical approach
            a_coeff = sqrt(sum(norm_in**2))/(8._dp*x_in(1)**2)
            do kd = 1,2
                if (kd.eq.1) then
                    delta_t(1) = t_s(1) - t_in_loc
                    delta_t(2) = t_s(2) - t_in_loc
                else
                    delta_t(1) = -(t_s(2) - t_in_loc)
                    delta_t(2) = -(t_s(1) - t_in_loc)
                end if
                do jd = 1,2
                    if (abs(delta_t(jd)).le.0._dp) then
                        dlogd(jd) = 0._dp
                    else
                        dlogd(jd) = log(a_coeff*abs(delta_t(jd)))
                    end if
                end do
                G(kd) = 1._dp/sqrt(x_in(1)*x_s(kd,1)) * (&
                    &b_coeff(2)*(delta_t(2)-delta_t(1)) + &
                    &0.5_dp*delta_t(1) - 1.5_dp*delta_t(2) + &
                    &1._dp/(delta_t(2)-delta_t(1)) * &
                    &(delta_t(1)*(delta_t(1)-2*delta_t(2))*dlogd(1) + &
                    &delta_t(2)**2 * dlogd(2)) &
                    &)
                H(kd) = -0.5*norm_s(kd,1)/x_s(kd,1) * G(kd) - &
                    &(delta_t(2)-delta_t(1)) * (n_tor-0.5_dp)**(-1) * &
                    &Aij(kd)/sqrt(x_in(1)*x_s(kd,1))
            end do
            
            ! Check if  we are at  the boundary  of the analytical  approach. In
            ! this  case there  is  an  extra component  due  to the  difference
            ! between the toroidal function and the analytical approximation.
            ! This  will be  implemented through  a trapezoidal  contribution of
            ! which one side is zero.
            if (maxval(rho2).gt.tol) then
                ! fill G's and H's
                if (rho2(2).gt.rho2(1)) then                                    ! right boundary
                    kd_bound = 2
                else                                                            ! left boundary
                    kd_bound = 1
                end if
                G(kd_bound) = G(kd_bound) - 0.5_dp*(t_s(2)-t_s(1))*2._dp/&
                    &sqrt(x_in(1)*x_s(kd_bound,1)) * (ql(kd_bound,2)+&
                    &log(a_coeff*abs(t_s(kd_bound)-t_in_loc))+b_coeff(2))
                H(kd_bound) = H(kd_bound) - 0.5_dp*(t_s(2)-t_s(1))*2._dp/&
                    &sqrt(x_in(1)*x_s(kd_bound,1)) * &
                    &( (-0.5_dp*norm_s(kd_bound,1)/x_s(kd_bound,1) - &
                    &(1._dp+rho2(kd_bound)/(2*x_in(1)*x_s(kd_bound,1)))*&
                    &Aij(kd_bound)) * (ql(kd_bound,2)+&
                    &log(a_coeff*abs(t_s(kd_bound)-t_in_loc))+b_coeff(2)) + &
                    &Aij(kd_bound) * (ql(kd_bound,1)+&
                    &log(a_coeff*abs(t_s(kd_bound)-t_in_loc))+b_coeff(1)))
            end if
        else
            do kd = 1,2
                ! fill G's and H's
                G(kd) = - (t_s(2)-t_s(1))/sqrt(x_in(1)*x_s(kd,1)) * &
                    &ql(kd,2)
                H(kd) = -(t_s(2)-t_s(1))/sqrt(x_in(1)*x_s(kd,1)) * &
                    &( (-0.5_dp*norm_s(kd,1)/x_s(kd,1) - &
                    &(1._dp+rho2(kd)/(2*x_in(1)*x_s(kd,1)))*&
                    &Aij(kd))*ql(kd,2) + Aij(kd)*ql(kd,1) )
            end do
        end if
    end subroutine calc_GH_int_2
    
    !> Gathers a distributed vector on a single process.
    !!
    !! For the distributed  vector, the limits of the different  subrows need to
    !! be provided, as discussed in init_vac().
    !!
    !! By  default, the  result lands  on the  master process,  but this  can be
    !! changed  using \c  proc. If  it is  negative, all  processes receive  the
    !! result. In any case, \c vec_loc will be utilized.
    !!
    !! \note \c vec_loc needs to be allocated to the total dimensions, even when
    !! the results is not received by the process.
    integer function vec_dis2loc(ctxt,vec_dis,lims,vec_loc,proc) result(ierr)
        use num_vars, only: n_procs
        use vac_vars, only: in_context

        ! input / output
        integer, intent(in) :: ctxt                                             !< context for vector
        real(dp), intent(in) :: vec_dis(:)                                      !< distributed vector
        integer, intent(in) :: lims(:,:)                                        !< limits for different subrows
        real(dp), intent(inout) :: vec_loc(:)                                   !< local vector
        integer, intent(in), optional :: proc                                   !< which process receives result
        
        ! local variables
        integer :: id                                                           ! index of subrow
        integer :: limsl(2)                                                     ! local limits
        integer :: proc_loc                                                     ! local proc
        integer :: proc_dest(2)                                                 ! destination process index
        
        ! initialize ierr
        ierr = 0
        
        ! select only processes that are in the context
        if (.not.in_context(ctxt)) return
        
        ! set up destination process index
        proc_loc = 0
        if (present(proc)) proc_loc = proc
        if (proc_loc.gt.0 .and. proc_loc.lt.n_procs) then
            call blacs_pcoord(ctxt,proc_loc,proc_dest(1),proc_dest(2)) 
        else
            proc_dest = [-1,-1]
        end if
        
        ! set up total copy of vector
        vec_loc = 0._dp
        do id = 1,size(lims,2)
            ! set local limits
            if (size(vec_dis) .gt. 0) then
                limsl = sum(lims(2,1:id-1)-lims(1,1:id-1)+1) + &
                    &lims(:,id)-lims(1,id)+1
                vec_loc(lims(1,id):lims(2,id)) = vec_dis(limsl(1):limsl(2))
            end if
        end do
        
        ! add them together
        call dgsum2d(ctxt,'all',' ',size(vec_loc),1,vec_loc,size(vec_loc),&
            &proc_dest(1),proc_dest(2))
    end function vec_dis2loc
    
    !> Gathers a distributed vector on a single process.
    !!
    !! \see See vec_dis2loc() for exaplanation.
    integer function mat_dis2loc(ctxt,mat_dis,lims_r,lims_c,mat_loc,proc) &
        &result(ierr)
        use num_vars, only: n_procs
        use vac_vars, only: in_context

        ! input / output
        integer, intent(in) :: ctxt                                             !< context for matrix
        real(dp), intent(in) :: mat_dis(:,:)                                    !< distributed matrix
        integer, intent(in) :: lims_r(:,:)                                      !< limits for different subrows
        integer, intent(in) :: lims_c(:,:)                                      !< limits for different subcolumns
        real(dp), intent(inout) :: mat_loc(:,:)                                 !< local matrix
        integer, intent(in), optional :: proc                                   !< which process receives result
        
        ! local variables
        integer :: i_rd, i_cd                                                   ! index of subrow and subcolumn
        integer :: limsl_r(2)                                                   ! local row limits
        integer :: limsl_c(2)                                                   ! local column limits
        integer :: proc_loc                                                     ! local proc
        integer :: proc_dest(2)                                                 ! destination process index
        
        ! initialize ierr
        ierr = 0
        
        ! select only processes that are in the context
        if (.not.in_context(ctxt)) return
        
        ! set up destination process index
        proc_loc = 0
        if (present(proc)) proc_loc = proc
        if (proc_loc.gt.0 .and. proc_loc.lt.n_procs) then
            call blacs_pcoord(ctxt,proc_loc,proc_dest(1),proc_dest(2)) 
        else
            proc_dest = [-1,-1]
        end if
        
        ! set up total copy of matrix
        mat_loc = 0._dp
        do i_rd = 1,size(lims_r,2)
            if (size(mat_dis,1) .gt. 0) then
                limsl_r = sum(&
                    &lims_r(2,1:i_rd-1)-lims_r(1,1:i_rd-1)+1) + &
                    &lims_r(:,i_rd)-lims_r(1,i_rd)+1
                do i_cd = 1,size(lims_c,2)
                    if (size(mat_dis,2) .gt. 0) then
                        limsl_c = sum(&
                            &lims_c(2,1:i_cd-1)-lims_c(1,1:i_cd-1)+1) + &
                            &lims_c(:,i_cd)-lims_c(1,i_cd)+1
                        mat_loc(lims_r(1,i_rd):lims_r(2,i_rd),&
                            &lims_c(1,i_cd):lims_c(2,i_cd)) = mat_dis(&
                            &limsl_r(1):limsl_r(2),limsl_c(1):limsl_c(2))
                    end if
                end do
            end if
        end do
        
        ! add them together
        call dgsum2d(ctxt,'all',' ',size(mat_loc,1),size(mat_loc,2),mat_loc,&
            &size(mat_loc,1),proc_dest(1),proc_dest(2))
    end function mat_dis2loc
    
    !> Copies  vacuum  variables  of  a previous  level  into  the  appropriate,
    !! interlaced place of a next level.
    !! 
    !! This is done  by multiplying left by  \f$\overline{\text{T}}\f$ and right
    !! by \f$\overline{\text{T}}^T\f$ with
    !! \f[
    !!  a_{ij} =
    !!      \left\{\begin{aligned}
    !!          1 \quad &\text{for}  \left(i,j\right) =
    !!              \left(2j-1+(i-1)n_\text{old},j+(i-1)n_\text{old}\right) \\
    !!          0 \quad &\text{otherwise} 
    !!      \end{aligned}\right.
    !! \f]
    !! where \f$i\f$ ranges from \f$1\f$  to \f$n_\text{old}\f$ and \f$j\f$ from
    !! \f$1\f$ to \f$m_\text{old}\f$, with the size of \f$\overline{\text{T}}\f$
    !! equal to \f$n_\text{new} \times n_\text{old}  \f$ where \f$n\f$ refers to
    !! \c n_bnd(1) and \f$m\f$ to \c n_bnd(2).
    integer function interlaced_vac_copy(vac_old,vac) result(ierr)
#if ldebug
        use num_vars, only: ltest, n_procs, rank
        use input_utilities, only: get_log
#endif
        
        character(*), parameter :: rout_name = 'interlaced_vac_copy'
        
        ! input / output
        type(vac_type), intent(in) :: vac_old                                   !< old vacuum
        type(vac_type), intent(inout) :: vac                                    !< vacuum
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: i_cd                                                         ! index of subcol
        integer :: i_rd                                                         ! index of subrow
        integer :: rd, cd                                                       ! global counter for row and col
        integer :: rdl, cdl                                                     ! local counter for row and col
        integer :: alpha_id                                                     ! field line on which an index is situated
        integer :: par_id                                                       ! index point on field line
        integer :: n_ang(2)                                                     ! number of points in angular directions
        integer :: n_ang_old(2)                                                 ! number of points in angular directions in old vacuum
        real(dp), allocatable :: loc_A(:,:)                                     ! local transformation matrix
        real(dp), allocatable :: loc_B(:,:)                                     ! local dummy matrix
        integer :: desc_loc(BLACSCTXTSIZE)                                      ! descriptor for loc_A and loc_B
        integer :: desc_loc_old(BLACSCTXTSIZE)                                  ! descriptor for loc_A and loc_B in old vacuum context
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        logical :: test_redist                                                  ! whether to test the redistribution of H and G
        real(dp), allocatable :: HG_ser(:,:)                                    ! serial versions of H and G
        real(dp), allocatable :: HG_ser_old(:,:)                                ! serial versions of H and G in old vacuum
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! initialize
        n_ang = vac%n_ang
        n_ang_old = vac_old%n_ang
        
        ! test
        if (n_ang(2).ne.n_ang_old(2)) then
            ierr = 1
            err_msg = 'Old and new vacuum need to have an equal number of &
                &field lines'
            CHCKERR(err_msg)
        end if
        if (n_ang(1).ne.2*n_ang_old(1)-1) then
            ierr = 1
            err_msg = 'Old and new vacuum need to have a compatible number of &
                &points along the field lines'
            CHCKERR(err_msg)
        end if
        
        ! loop over all field lines
        do id = 1,n_ang(2)
            ! copy normal and position vector
            vac%norm((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = &
                &vac_old%norm((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
            vac%x_vec((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = &
                &vac_old%x_vec((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
            
            vac%h_fac((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = vac_old%&
                &h_fac((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
            vac%dnorm((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = &
                &vac_old%dnorm((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
        end do
        ! (vac%ang, points along a line times lines, is set for the whole new
        ! level by store_vac_VMEC; the flattened indices above would write
        ! past it from the second field line on)
        
        ! set up transformation matrix A and dummy matrix B
        allocate(loc_A(vac%n_loc(1),vac_old%n_loc(2)))
        allocate(loc_B(vac%n_loc(1),vac_old%n_loc(2)))
        loc_A = 0._dp
        loc_B = 0._dp
        subcols: do i_cd = 1,size(vac_old%lims_c,2)
            col: do cd = vac_old%lims_c(1,i_cd),vac_old%lims_c(2,i_cd)
                ! set field line index and point index and calculate row
                alpha_id = (cd-1)/n_ang_old(1)+1
                par_id = mod(cd-1,n_ang_old(1))+1
                rd = 2*par_id-1 + (alpha_id-1)*n_ang(1)
                
                ! set local column index
                cdl = loc_idx(vac_old%lims_c,i_cd,cd)
                
                subrows: do i_rd = 1,size(vac%lims_r,2)
                    ! set local row index if in this subrow
                    if (rd.ge.vac%lims_r(1,i_rd) .and. &
                        &rd.le.vac%lims_r(2,i_rd)) then
                        rdl = loc_idx(vac%lims_r,i_rd,rd)
                        
                        loc_A(rdl,cdl) = 1._dp
                    end if
                end do subrows
            end do col
        end do subcols
        call descinit(desc_loc,vac%n_bnd,vac_old%n_bnd,vac%bs,&
            &vac%bs,0,0,vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
        CHCKERR('descinit failed for loc')
        call descinit(desc_loc_old,vac%n_bnd,vac_old%n_bnd,vac%bs,&
            &vac%bs,0,0,vac_old%ctxt_HG,max(1,vac%n_loc(1)),ierr)
        CHCKERR('descinit failed for loc_old')
        
        ! treat H
        CHCKERR('')
        call pdgemm('N','N',vac%n_bnd,vac_old%n_bnd,vac_old%n_bnd,&
            &1._dp,loc_A,1,1,desc_loc_old,vac_old%H,1,1,&
            &vac_old%desc_H,0._dp,loc_B,1,1,desc_loc_old)
        call pdgemm('N','T',vac%n_bnd,vac%n_bnd,vac_old%n_bnd,&
            &1._dp,loc_B,1,1,desc_loc,loc_A,1,1,&
            &desc_loc,0._dp,vac%H,1,1,vac%desc_H)
        
        ! treat G
        call pdgemm('N','N',vac%n_bnd,vac_old%n_bnd,vac_old%n_bnd,&
            &1._dp,loc_A,1,1,desc_loc_old,vac_old%G,1,1,&
            &vac_old%desc_G,0._dp,loc_B,1,1,desc_loc_old)
        call pdgemm('N','T',vac%n_bnd,vac%n_bnd,vac_old%n_bnd,&
            &1._dp,loc_B,1,1,desc_loc,loc_A,1,1,&
            &desc_loc,0._dp,vac%G,1,1,vac%desc_G)
        
#if ldebug
        ! test
        if (ltest) then
            call writo('test resdistribution of H?')
            test_redist = get_log(.false.)
        else
            test_redist = .false.
        end if
        if (test_redist) then
            allocate(HG_ser(vac%n_bnd,vac%n_bnd))
            allocate(HG_ser_old(vac_old%n_bnd,vac_old%n_bnd))
            ierr = mat_dis2loc(vac_old%ctxt_HG,vac_old%H,&
                &vac_old%lims_r,vac_old%lims_c,HG_ser_old,&
                &proc=n_procs-1)
            CHCKERR('')
            ierr = mat_dis2loc(vac%ctxt_HG,vac%H,&
                &vac%lims_r,vac%lims_c,HG_ser,&
                &proc=n_procs-1)
            CHCKERR('')
            if (rank.eq.n_procs-1) then
                call writo('old H:',persistent=rank.eq.n_procs-1)
                call print_ar_2(HG_ser_old)
                call writo('new H:',persistent=rank.eq.n_procs-1)
                call print_ar_2(HG_ser)
                call plot_HDF5('HG_ser_old','HG_ser_old',&
                    &reshape(HG_ser_old,[vac_old%n_bnd,vac_old%n_bnd,1]))
                call plot_HDF5('HG_ser','HG_ser',&
                    &reshape(HG_ser,[vac%n_bnd,vac%n_bnd,1]))
            end if
        end if
#endif
        
        ! clean up
        deallocate(loc_A,loc_B)
    end function interlaced_vac_copy
    !> \public Number of local columns this process owns of a global column vector:
    !! the single column sits in the first block, owned by process column 0.
    integer function n_col1_loc(vac) result(res)
        use vac_vars, only: in_context
        
        type(vac_type), intent(in) :: vac                                       !< vacuum variables

        if (in_context(vac%ctxt_HG) .and. vac%ind_p(2).eq.0) then
            res = 1
        else
            res = 0
        end if
    end function n_col1_loc

    !> \public Scatter a global column vector into its local block-cyclic rows.
    !!
    !! \c x_dis is the flattened local array of size n_loc(1) * n_col1_loc();
    !! processes that do not own the column receive nothing.
    subroutine vec_glob2dis(vac,x_glob,x_dis)
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: x_glob(:)                                       !< global vector
        real(dp), intent(out) :: x_dis(:)                                       !< local part (flattened)

        integer :: i_rd, rd, rdl                                                ! subrow index, global and local row

        x_dis = 0._dp
        if (n_col1_loc(vac).eq.0) return
        do i_rd = 1,size(vac%lims_r,2)
            do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                rdl = loc_idx(vac%lims_r,i_rd,rd)
                x_dis(rdl) = x_glob(rd)
            end do
        end do
    end subroutine vec_glob2dis

    !> \public Distributed matrix-vector product y = A x.
    !!
    !! \c A is a local part of a distributed n_bnd x n_bnd matrix with
    !! descriptor \c desc_A (e.g. vac%H with vac%desc_H), \c x_glob is
    !! globally replicated, and \c y_glob is returned complete on all
    !! processes in the context (zero outside it). Used by the exact-potential
    !! checks of the vacuum operators (vac_ops.check_vac_exact() and the
    !! full-stack tests).
    integer function dis_matvec(vac,A,desc_A,x_glob,y_glob) result(ierr)
        use vac_vars, only: in_context

        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: A(:,:)                                          !< local part of distributed matrix
        integer, intent(in) :: desc_A(BLACSCTXTSIZE)                            !< descriptor of A
        real(dp), intent(in) :: x_glob(:)                                       !< global input vector
        real(dp), intent(out) :: y_glob(:)                                      !< global result vector

        integer :: n_col1                                                       ! local number of columns of the vector
        integer :: n_flat                                                       ! true local size of the vector
        integer :: desc_x(BLACSCTXTSIZE)                                        ! descriptor for x and y
        real(dp), allocatable :: x_dis(:), y_dis(:)                             ! local vector parts

        ierr = 0
        y_glob = 0._dp
        if (.not.in_context(vac%ctxt_HG)) return

        n_col1 = n_col1_loc(vac)
        n_flat = vac%n_loc(1)*n_col1
        allocate(x_dis(max(1,n_flat)),y_dis(max(1,n_flat)))                     ! at least 1 element for ScaLAPACK dummies
        call vec_glob2dis(vac,x_glob,x_dis)
        y_dis = 0._dp

        call descinit(desc_x,vac%n_bnd,1,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
            &max(1,vac%n_loc(1)),ierr)
        if (ierr.ne.0) return

        call pdgemv('N',vac%n_bnd,vac%n_bnd,1._dp,A,1,1,desc_A,&
            &x_dis,1,1,desc_x,1,0._dp,y_dis,1,1,desc_x,1)

        ierr = vec_dis2loc(vac%ctxt_HG,y_dis(1:n_flat),vac%lims_r,y_glob)       ! all processes receive
    end function dis_matvec
    
    !> \public Step sizes of the field-line grid in the parallel coordinate
    !! and in the field-line label, as used by the field-line vacuum: the
    !! parallel points include both ends, the field lines exclude the end
    !! point in alpha (rich_ops).
    subroutine vac_steps(dpar,dalpha)
        use grid_vars, only: min_par_X, max_par_X, min_alpha, max_alpha, n_alpha
        use rich_vars, only: n_par_X
        
        ! input / output
        real(dp), intent(out) :: dpar                                           !< step in the parallel coordinate
        real(dp), intent(out) :: dalpha                                         !< step in the field-line label
        
        dpar = (max_par_X-min_par_X)/(n_par_X-1)*pi
        dalpha = (max_alpha-min_alpha)/max(n_alpha,1)*pi
    end subroutine vac_steps
    
    !> \public Cross product of two 3-vectors.
    pure function cross3(u,v) result(w)
        ! input / output
        real(dp), intent(in) :: u(3), v(3)                                      !< vectors
        real(dp) :: w(3)                                                        !< their cross product
        
        w = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
    end function cross3
    
    !> \public Solid angle at a node of the mesh, on the side of \c n_hat.
    !!
    !! The surface of the  mesh at a node is the cone of  four flat facets, one
    !! between  each pair  of edges  that are  consecutive in  azimuth about  \c
    !! n_hat. The  solid angle of the  region on the \c  n_hat side of all  four
    !! facets is  eq. (45) of the  2019 manuscript on the  vacuum energy (exact;
    !! not its small-tilt form (46), which fails on skewed meshes, erratum E6):
    !! with  the facet  normals  \f$\mathbf{a}_i\f$ tilted  by \f$\beta_i\f$  from
    !! \c n_hat towards  the azimuth \f$\bar\varphi_i\f$, and  \f$\varphi_i\f$ the
    !! azimuths of the edges,
    !! \f[ \Omega = 2\pi + \sum_{i=1}^4 \left[\sin^{-1}\left(\sin\beta_i
    !!  \sin(\varphi-\bar\varphi_i)\right)\right]_{\varphi_i}^{\varphi_{i+1}} . \f]
    !! The edges may be given in any order. With \c n_hat the unit normal into
    !! the plasma, the result  is the plasma's solid angle  \f$\Omega_p\f$ and
    !! the vacuum's is \f$4\pi - \Omega_p\f$.
    pure subroutine solid_angle_node(u,n_hat,omega)
        ! input / output
        real(dp), intent(in) :: u(4,3)                                          !< the four edge vectors from the node to its neighbours
        real(dp), intent(in) :: n_hat(3)                                        !< unit normal of the side whose solid angle is returned
        real(dp), intent(out) :: omega                                          !< solid angle on that side
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: ord(4)                                                       ! edges sorted by azimuth
        integer :: i1, i2                                                       ! edges of a facet
        real(dp) :: z(3), x(3), y(3)                                            ! local frame
        real(dp) :: uh(4,3)                                                     ! unit edges
        real(dp) :: ph(4)                                                       ! azimuths of the edges
        real(dp) :: a(3)                                                        ! facet normal
        real(dp) :: c, beta, phb                                                ! cos of the tilt, tilt and azimuth of the facet normal
        real(dp) :: p1, p2                                                      ! azimuths of the sector
        
        ! local frame: z along n_hat, x along the tangential part of the first edge
        z = n_hat/norm2(n_hat)
        x = u(1,:) - sum(u(1,:)*z)*z
        x = x/norm2(x)
        y = cross3(z,x)
        
        ! unit edges and their azimuths
        do id = 1,4
            uh(id,:) = u(id,:)/norm2(u(id,:))
            ph(id) = atan2(sum(uh(id,:)*y),sum(uh(id,:)*x))
        end do
        
        ! sort the edges by azimuth (insertion sort)
        ord = [1,2,3,4]
        do id = 2,4
            kd = ord(id)
            jd = id-1
            do while (jd.ge.1)
                if (ph(ord(jd)).le.ph(kd)) exit
                ord(jd+1) = ord(jd)
                jd = jd-1
            end do
            ord(jd+1) = kd
        end do
        
        ! sum over the four sectors
        omega = 2*pi
        do id = 1,4
            i1 = ord(id)
            i2 = ord(mod(id,4)+1)
            a = cross3(uh(i1,:),uh(i2,:))
            a = a/norm2(a)
            c = sum(a*z)
            if (c.lt.0._dp) then                                                ! normal of the facet on the n_hat side
                a = -a
                c = -c
            end if
            beta = acos(min(c,1._dp))
            phb = atan2(sum(a*y),sum(a*x))
            p1 = ph(i1)
            p2 = ph(i2)
            if (p2.le.p1) p2 = p2 + 2*pi                                        ! the sector that closes the circle
            omega = omega + asin(sin(beta)*sin(p2-phb)) - &
                &asin(sin(beta)*sin(p1-phb))
        end do
    end subroutine solid_angle_node
    
    !> \public Integral of the H kernel over the cell of a node of the
    !! field-line mesh.
    !!
    !! The kernel  \f$\mathbf{n}(\mathbf{r})\cdot(\mathbf{r}-\mathbf{r}_0)/
    !! |\mathbf{r}-\mathbf{r}_0|^3\f$ (with \f$\mathbf{n}\f$ the stored normal
    !! \f$-\mathcal{J}\nabla\psi\f$ at the source) is  bounded on the surface
    !! near its  own node, because  the surface leaves its  tangent plane at
    !! second  order: to  leading order  \f$\mathbf{n}(\mathbf{r})\cdot
    !! (\mathbf{r}-\mathbf{r}_0) = -\tfrac12\, \mathbf{n}_0\cdot\partial_j
    !! \mathbf{e}_k\, \Delta u^j \Delta u^k\f$ (\f$j,k\f$ the parallel and
    !! field-line coordinates). Over the cell of the node, \f$|\theta-\theta_0|
    !! \le \Delta\theta/2\f$, \f$|\alpha-\alpha_0| \le \Delta\alpha/2\f$,
    !! in the local polar  coordinates of eqs. (49)-(72) of  the 2019 manuscript
    !! (with \f$A = g_{\alpha\theta}/|\mathbf{n}|\f$, \f$B = g_{\alpha\alpha}/
    !! |\mathbf{n}|\f$ and \f$\rho(\beta)\f$ the  radius of the cell), the
    !! integral is
    !! \f[ I = -\frac{\sqrt{g_{\alpha\alpha}}}{2|\mathbf{n}|^2}
    !!  \int_{-\pi/2}^{\pi/2} q(\beta)\,\rho(\beta)\,\mathrm{d}\beta, \quad
    !!  q = \mathbf{n}_0\cdot\partial_j\mathbf{e}_k\, t^j t^k, \quad
    !!  t = \left(\cos\beta, \frac{\sin\beta - A\cos\beta}{B}\right), \f]
    !! (with \f$A\f$ and \f$B\f$ taken positive: the mirror cell, the frame of the
    !! opposite orientation, is an isometric image and gives the same integral),
    !! which is  the corrected eq.  (69) of the manuscript  (erratum E22: the
    !! printed one has a  wrong factor and only the umbilic part  of the second
    !! fundamental form). The half-cell of the manuscript equals the full cell
    !! of the node for this kernel.  Verified against direct integration on
    !! sphere and torus nodes (first-order convergence), and by the full-stack
    !! test \c H_sing_node.
    pure subroutine H_sing_node(norm0,g,dd,dpar,dalpha,I_sing)
        ! input / output
        real(dp), intent(in) :: norm0(3)                                        !< stored normal at the node
        real(dp), intent(in) :: g(3)                                            !< metric factors g_aa, g_at, g_tt at the node
        real(dp), intent(in) :: dd(3,3)                                         !< second derivatives of the position: d_t e_t, d_t e_a, d_a e_a (rows)
        real(dp), intent(in) :: dpar                                            !< parallel step
        real(dp), intent(in) :: dalpha                                          !< step in the field-line label
        real(dp), intent(out) :: I_sing                                         !< integral of the kernel over the cell
        
        ! local variables
        integer :: ip, jd                                                       ! counters
        real(dp) :: nn                                                          ! |norm0|
        real(dp) :: A, B                                                        ! coefficients of the local frame
        real(dp) :: q_tt, q_ta, q_aa                                            ! second fundamental form along the stored normal
        real(dp) :: lims(4)                                                     ! limits of the three pieces of the cell boundary
        real(dp) :: p_lo, p_hi                                                  ! limits of a piece
        real(dp) :: beta, c, s                                                  ! polar angle, cos, sin-like coordinate
        real(dp) :: q, rho                                                      ! quadratic form and radius
        real(dp) :: tot                                                         ! integral
        
        nn = norm2(norm0)
        A = g(2)/nn
        B = g(1)/nn
        q_tt = sum(norm0*dd(1,:))
        q_ta = sum(norm0*dd(2,:))
        q_aa = sum(norm0*dd(3,:))
        lims = [-pi/2, atan(-B*dalpha/dpar+A), atan(B*dalpha/dpar+A), pi/2]
        tot = 0._dp
        do ip = 1,3
            p_lo = lims(ip)
            p_hi = lims(ip+1)
            do jd = 1,n_gl16
                beta = 0.5_dp*(p_lo+p_hi) + 0.5_dp*(p_hi-p_lo)*gl16_x(jd)
                c = cos(beta)
                s = (sin(beta)-A*c)/B                                           ! field-line coordinate of the direction
                q = q_tt*c**2 + 2*q_ta*c*s + q_aa*s**2
                select case (ip)
                    case (1)                                                    ! -pi/2 .. beta_1: lower edge alpha = alpha_0 - dalpha
                        rho = B*dalpha/(A*c-sin(beta))
                    case (2)                                                    ! beta_1 .. beta_2: far edge theta = theta_0 + dpar
                        rho = dpar/c
                    case default                                                ! beta_2 .. pi/2: upper edge alpha = alpha_0 + dalpha
                        rho = B*dalpha/(sin(beta)-A*c)
                end select
                tot = tot + 0.5_dp*(p_hi-p_lo)*gl16_w(jd)*q*rho
            end do
        end do
        I_sing = -sqrt(g(1))/(2*nn**2)*tot
    end subroutine H_sing_node
    
    !> \public Direct diagonals of H for the field-line vacuum.
    !!
    !! Two variants, both with the cell integral of the kernel (H_sing_node):
    !!  - smooth: \f$-2\pi\f$ for the corner term, the solid angle of the
    !!    vacuum at a point of the smooth surface, consistent with the kernels
    !!    of calc_GH_int_1, which are evaluated on the smooth surface;
    !!  - faceted: minus the vacuum's solid angle at the node of the faceted
    !!    surface, eq. (45) of the 2019 manuscript (solid_angle_node), which
    !!    belongs to a discretization with flat elements; on a skewed mesh,
    !!    such as the field-aligned one, the facets tilt strongly (erratum E6)
    !!    and this value is far from \f$2\pi\f$.
    !! The four edges run to the neighbours along the field line and on the
    !! neighbouring lines (periodic in the field-line label). The second
    !! derivatives of the position come from second differences of the stored
    !! points. The end points of the field lines have only three edges and
    !! half a cell, and are left out (\c avail false).
    integer function calc_diag_direct(vac,dpar,dalpha,diag_direct,diag_facet,&
        &I_sing,omega_v,avail) result(ierr)
        
        character(*), parameter :: rout_name = 'calc_diag_direct'
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: dpar                                            !< parallel step
        real(dp), intent(in) :: dalpha                                          !< step in the field-line label
        real(dp), intent(out) :: diag_direct(:)                                 !< direct diagonal of H for the smooth surface: -2 pi plus the cell integral
        real(dp), intent(out) :: diag_facet(:)                                  !< the same for the faceted surface: minus the vacuum's solid angle at the node plus the cell integral
        real(dp), intent(out) :: I_sing(:)                                      !< its cell-integral part
        real(dp), intent(out) :: omega_v(:)                                     !< the vacuum's solid angle at the node
        logical, intent(out) :: avail(:)                                        !< whether the point has a direct value
        
        ! local variables
        integer :: n_par, n_lin                                                 ! points per line, number of lines
        integer :: ip, il, id                                                   ! point, line and global indices
        integer :: ilm, ilp                                                     ! previous and next line (periodic)
        integer :: idm, idp                                                     ! previous and next point on the line
        real(dp) :: x0(3)                                                       ! position of the node
        real(dp) :: u(4,3)                                                      ! edges
        real(dp) :: n_hat(3)                                                    ! unit normal into the plasma
        real(dp) :: dd(3,3)                                                     ! second derivatives
        real(dp) :: omega_p                                                     ! the plasma's solid angle
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize
        ierr = 0
        diag_direct = 0._dp
        diag_facet = 0._dp
        I_sing = 0._dp
        omega_v = 0._dp
        avail = .false.
        n_par = vac%n_ang(1)
        n_lin = vac%n_ang(2)
        if (vac%style.ne.1) then
            ierr = 1
            err_msg = 'the direct diagonal of calc_diag_direct is for the &
                &field-line vacuum'
            CHCKERR(err_msg)
        end if
        if (n_lin.lt.3 .or. n_par.lt.3) then
            ierr = 1
            err_msg = 'the direct diagonal needs at least 3 field lines and &
                &3 points per line'
            CHCKERR(err_msg)
        end if
        
        do il = 1,n_lin
            ilm = modulo(il-2,n_lin)+1
            ilp = modulo(il,n_lin)+1
            do ip = 2,n_par-1
                id = ip + (il-1)*n_par
                idm = id-1
                idp = id+1
                x0 = vac%x_vec(id,:)
                ! edges: -alpha, -theta, +alpha, +theta
                u(1,:) = vac%x_vec(ip+(ilm-1)*n_par,:) - x0
                u(2,:) = vac%x_vec(idm,:) - x0
                u(3,:) = vac%x_vec(ip+(ilp-1)*n_par,:) - x0
                u(4,:) = vac%x_vec(idp,:) - x0
                ! solid angle of the plasma, then of the vacuum
                n_hat = vac%norm(id,:)/norm2(vac%norm(id,:))
                call solid_angle_node(u,n_hat,omega_p)
                omega_v(id) = 4*pi - omega_p
                ! second derivatives from second differences
                dd(1,:) = (vac%x_vec(idp,:) - 2*x0 + vac%x_vec(idm,:))/dpar**2
                dd(3,:) = (u(3,:) + u(1,:))/dalpha**2
                dd(2,:) = (vac%x_vec(ip+1+(ilp-1)*n_par,:) - &
                    &vac%x_vec(ip-1+(ilp-1)*n_par,:) - &
                    &vac%x_vec(ip+1+(ilm-1)*n_par,:) + &
                    &vac%x_vec(ip-1+(ilm-1)*n_par,:))/(4*dpar*dalpha)
                call H_sing_node(vac%norm(id,:),vac%h_fac(id,1:3),dd,dpar,&
                    &dalpha,I_sing(id))
                diag_direct(id) = I_sing(id) - 2*pi
                diag_facet(id) = I_sing(id) - omega_v(id)
                avail(id) = .true.
            end do
        end do
    end function calc_diag_direct
    
    !> \public Gathers the diagonal of a distributed matrix on all processes.
    integer function dis_diag(vac,A,diag_glob) result(ierr)
        use vac_vars, only: in_context
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: A(:,:)                                          !< local part of distributed matrix
        real(dp), intent(out) :: diag_glob(:)                                   !< its diagonal, global
        
        ! local variables
        integer :: i_rd, i_cd                                                   ! index of subrow and subcolumn
        integer :: rd                                                           ! global index on the diagonal
        
        ierr = 0
        diag_glob = 0._dp
        if (.not.in_context(vac%ctxt_HG)) return
        do i_rd = 1,size(vac%lims_r,2)
            do i_cd = 1,size(vac%lims_c,2)
                do rd = max(vac%lims_r(1,i_rd),vac%lims_c(1,i_cd)), &
                    &min(vac%lims_r(2,i_rd),vac%lims_c(2,i_cd))
                    diag_glob(rd) = A(loc_idx(vac%lims_r,i_rd,rd),&
                        &loc_idx(vac%lims_c,i_cd,rd))
                end do
            end do
        end do
        call dgsum2d(vac%ctxt_HG,'all',' ',size(diag_glob),1,diag_glob,&
            &size(diag_glob),-1,-1)
    end function dis_diag
end module vac_utilities
