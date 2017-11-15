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
    use vac_vars, only: vac_type

    implicit none
    private
    public calc_GH_int_2, vec_dis2loc, mat_dis2loc, in_context
    
contains
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
    !! \note This routine does not calculate the contribution \f$2\beta\f$.
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
                H(kd) = 0.5*norm_s(kd,1)/x_s(kd,1) * G(kd) + &
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
                H(kd_bound) = H(kd_bound) + 0.5_dp*(t_s(2)-t_s(1))*2._dp/&
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
                H(kd) = (t_s(2)-t_s(1))/sqrt(x_in(1)*x_s(kd,1)) * &
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
        
        character(*), parameter :: rout_name = 'vec_dis2loc'
        
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
        
        character(*), parameter :: rout_name = 'mat_dis2loc'
        
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
    
    !> Indicates whether current process is in the context.
    logical function in_context(ctxt) result(res)
        ! input / output
        integer, intent(in) :: ctxt                                             !< context for vector
        
        ! local variables
        integer :: n_p(2)                                                       ! nr. of processes in grid
        integer :: ind_p(2)                                                     ! index of local process in grid
        
        call blacs_gridinfo(ctxt,n_p(1),n_p(2),ind_p(1),ind_p(2)) 
        res = ind_p(1).ge.0
    end function in_context
end module vac_utilities
