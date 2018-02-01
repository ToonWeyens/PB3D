!------------------------------------------------------------------------------!
!> Numerical utilities related to the vacuum response.
!!
!! \see vac_ops for information.
!------------------------------------------------------------------------------!
module vac_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
    use StrumpackDensePackage
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use vac_vars, only: vac_type

    implicit none
    private
    public calc_GH_int_1, calc_GH_int_2, vec_dis2loc, mat_dis2loc, &
        &interlaced_vac_copy
    
contains
    !> Calculate  G_ij  and   H_ij  on  an  interval  for   field  line  aligned
    !! configurations.
    !!
    !! The indices for the  source variables are <tt>[left:right,dim]</tt> where
    !! \c dim is the Cartesian dimension. The same holds for \c ql.
    !!
    !! For  subintegrals  close  to  the singularity,  the  procedure  uses  the
    !! analytical approximation.
    !!
    !! \note This routine does not calculate the contribution \f$2\beta\f$.
    subroutine calc_GH_int_1(G,H,x_s,x_in,norm_s,h_fac_in,step_size,tol)
        ! input / output
        real(dp), intent(inout) :: G(2)                                         !< G
        real(dp), intent(inout) :: H(2)                                         !< H
        real(dp), intent(in) :: x_s(2,3)                                        !< position vector at left and right limit of source interval
        real(dp), intent(in) :: x_in(3)                                         !< position vector at which to calculate influence
        real(dp), intent(in) :: norm_s(2,3)                                     !< normal vector at left and right limit of source interval
        real(dp), intent(in) :: h_fac_in(4)                                     !< metric factors at which to calculate influence
        real(dp), intent(in) :: step_size(2)                                    !< step sizes in parallel direction and alpha
        real(dp), intent(in) :: tol                                             !< tolerance on distance between points
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp) :: r2(2)                                                       ! distance between source and influence points
        real(dp) :: E_fac                                                       ! factor E
        real(dp) :: sineps                                                      ! sin(eps)
        real(dp) :: dum_fac(2)                                                  ! dummy factors
        
        ! initialize
        G = 0._dp
        H = 0._dp
        
        ! calculate distance between source points
        do kd = 1,2
            r2(kd) = sum((x_s(kd,:)-x_in)**2)
        end do
        
        ! check for (near-)singular point
        if (minval(r2).le.tol) then
            E_fac = sqrt(h_fac_in(1)/h_fac_in(3))*step_size(2)/step_size(1)     ! |e_alpha|/|e_theta| d_par_X / d_par_X
            sineps = h_fac_in(2)/sqrt(h_fac_in(1)*h_fac_in(3))
            dum_fac(1) = sqrt(1+2*E_fac*sineps+E_fac**2)
            dum_fac(2) = sqrt(1-2*E_fac*sineps+E_fac**2)
            G = -0.5_dp/sqrt(h_fac_in(1)) * step_size(1) * (&
                &log(abs((dum_fac(2) + E_fac + sineps)/&
                &(dum_fac(1) - E_fac + sineps))) &
                &- E_fac/sineps**2* &
                &log(abs((dum_fac(2) + E_fac + sineps)/&
                &(dum_fac(1) - E_fac + sineps))) &
                &)
            H = h_fac_in(4)*G
        else
            G = -1._dp/sqrt(r2)
            do kd = 1,2
                H(kd) = sum(norm_s(kd,:)*(x_s(kd,:)-x_in))*(-G(kd))**(-3)
            end do
        end if
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
        use vac_vars, only: in_context
        
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
        use num_vars, only: n_procs
        use vac_vars, only: in_context
        
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
            
            ! copy variables specific to style
            select case (vac%style)
                case (1)                                                        ! field-line 3-D
                    vac%h_fac((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = vac_old%&
                        &h_fac((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
                case (2)                                                        ! axisymmetric
                    vac%dnorm((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = &
                        &vac_old%dnorm((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
                    vac%ang((id-1)*n_ang(1)+1:id*n_ang(1):2,:) = &
                        &vac_old%ang((id-1)*n_ang_old(1)+1:id*n_ang_old(1):1,:)
            end select
        end do
        
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
                cdl = sum(vac_old%lims_c(2,1:i_cd-1)-&
                    &vac_old%lims_c(1,1:i_cd-1)+1) + &
                    &cd-vac_old%lims_c(1,i_cd)+1
                
                subrows: do i_rd = 1,size(vac%lims_r,2)
                    ! set local row index if in this subrow
                    if (rd.ge.vac%lims_r(1,i_rd) .and. &
                        &rd.le.vac%lims_r(2,i_rd)) then
                        rdl = sum(vac%lims_r(2,1:i_rd-1)-&
                            &vac%lims_r(1,1:i_rd-1)+1) + &
                            &rd-vac%lims_r(1,i_rd)+1
                        
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
end module vac_utilities
