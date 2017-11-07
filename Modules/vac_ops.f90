!------------------------------------------------------------------------------!
!> Operations and variables pertaining to the vacuum response.
!!
!! The vacuum  response is calculated  using the Boundary Element  Method. There
!! are two possibilities, indicated by the variable \c style in \c vac_vars.
!!  -# 3-D field aligned:
!!      - Makes  use of  the collocation  method for the  grid points  along the
!!      magnetic field lines.
!!      - The 3-D integrals are therefore integrated using Weyl's theorem
!!      \cite Helander2014, p. 5.
!!  -# axisymmetric:
!!      - Makes use of an analytical  expression for the toroidal integration of
!!      the Green's functions, as shown in \cite Cohl1999.
!!      - The  integration in the poloidal  angle is done using  the collocation
!!      method.
!!
!! The final matrix equation is solved through Strumpack \cite Meiser2016.
!!
!! \see See \cite Weyens3D.
!!
!! \todo The vacuum part of PB3D is still under construction and not usable yet.
!------------------------------------------------------------------------------!
module vac_ops
#include <PB3D_macros.h>
#include <wrappers.h>
    use StrumpackDensePackage
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_2_type
    use vac_vars, only: copy_vac, &
        &vac_type

    implicit none
    private
    public store_vac, calc_vac, print_output_vac
#if ldebug
    public debug_calc_GH, debug_calc_vac
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_GH = .false.                                          !< plot debug information for calc_GH() \ldebug
    logical :: debug_calc_vac = .false.                                          !< plot debug information for calc_vac() \ldebug
#endif

contains
    !> Stores  calculation  of the  vacuum  response  by storing  the  necessary
    !! variables.
    !!
    !! This is  done by  the last process,  which is the  one that  contains the
    !! variables at the  plasma edge, and then this is  broadcasted to the other
    !! processes.
    !!
    !! For the first  equilibrium job, this routine copies the  results from the
    !! previous  Richardson level  into  the appropriate  subranges  of the  new
    !! vacuum variables.
    !!
    !! If the previous \c G and \c  H variables are empty, they are regenerated.
    !! This indicates that  either a Richardson restart was performed,  or a new
    !! Richardson  level was  started after  a previous  level in  which it  was
    !! jumped straight to the solution.
    !!
    !! \note For HELENA, there is only 1  equilibrium job, and the vacuum has to
    !! be calculated only  once. If there is  a jump to solution or  if there is
    !! Richardson restart, it needs to be reconstructed only.
    !!
    !! \return ierr
    integer function store_vac(grid,eq_1,eq_2,vac) result(ierr)
        use PB3D_ops, only: reconstruct_PB3D_vac
        use num_vars, only: rich_restart_lvl, eq_job_nr, eq_jobs_lims, &
            &eq_style, rank, n_procs, use_pol_flux_F
        use rich_vars, only: rich_lvl
        use MPI_utilities, only: broadcast_var
        use num_utilities, only: c
        use grid_vars, only: n_r_eq
        use X_vars, only: prim_X
        
        character(*), parameter :: rout_name = 'store_vac'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        integer :: id, jd                                                       ! counter
        integer :: r_id                                                         ! normal index
        real(dp) :: jq                                                          ! iota or q
        
        ! initialize ierr
        ierr = 0
        
        ! set jq
        if (use_pol_flux_F) then
            jq = eq_1%q_saf_FD(grid%loc_n_r,0)
        else
            jq = eq_1%rot_t_FD(grid%loc_n_r,0)
        end if
        ierr = broadcast_var(jq,n_procs-1)
        CHCKERR('')
        
        ! call the approriate procedure
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = store_vac_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = store_vac_HEL()
                CHCKERR('')
        end select
    contains
        ! VMEC version
        !> \private
        integer function store_vac_VMEC() result(ierr)
            use rich_vars, only: n_par_X
            use eq_utilities, only: calc_inv_met
            
            character(*), parameter :: rout_name = 'store_vac_VMEC'
            
            ! local variables
            integer :: eq_job_lims(2)                                           ! local eq_jobs_lims
            real(dp), allocatable :: norm_com_C(:,:)                            ! Cylindrical components of norm
            real(dp), allocatable :: T_CV_loc(:,:,:,:,:,:,:)                    ! local T_CV
            type(vac_type) :: vac_old                                           ! old vacuum variables
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Start storing vacuum quantities')
            
            call lvl_ud(1)
            
            ! if start of Richardson level  that is not the first, copy previous
            ! results
            if (eq_job_nr.eq.1) then
                if (rich_lvl.gt.1) then
                    if (rich_lvl.eq.rich_restart_lvl) then
                        ! reconstruct old vacuum
                        ierr = reconstruct_PB3D_vac(vac_old,'vac',&
                            &rich_lvl=rich_lvl-1)
                        CHCKERR('')
                    else
                        ! copy old vacuum
                        ierr = copy_vac(vac,vac_old)
                        CHCKERR('')
                        
                        ! deallocate
                        call vac%dealloc()
                    end if
                    
                    if (.not.allocated(vac_old%H)) then
                        ! recalculate old G and H
                        ierr = -1
                        CHCKERR('SET UP OLD G AND H')
                    end if
                end if
                
                ! allocate
                ierr = vac%init(eq_style,n_par_X,prim_X,[n_par_X,1],jq)
                CHCKERR('')
                
                !! copy previous results back in appropriate, interlaced place
                !if (rich_lvl.gt.1) then
                    !vac%norm(1:vac%n_bnd:2,:) = vac_old%norm
                    !vac%x_vec(1:vac%n_bnd:2,:) = vac_old%x_vec
                    !ierr = 1
                    !CHCKERR('DO THIS WITH PDGEMR2D')
                !end if
            end if
            
            ! add results from current equilibrium job to the variables
            eq_job_lims = eq_jobs_lims(:,eq_job_nr)
            if (rank.eq.n_procs-1) then
                ! grid point to save is last local
                r_id = grid%loc_n_r
                
                ! save X, Y and Z
                vac%x_vec(eq_job_lims(1):eq_job_lims(2),1) = &
                    &eq_2%R_E(:,1,r_id,0,0,0) * cos(grid%zeta_E(:,1,r_id))
                vac%x_vec(eq_job_lims(1):eq_job_lims(2),2) = &
                    &eq_2%R_E(:,1,r_id,0,0,0) * sin(grid%zeta_E(:,1,r_id))
                vac%x_vec(eq_job_lims(1):eq_job_lims(2),3) = &
                    &eq_2%Z_E(:,1,r_id,0,0,0)
                
                ! calculate Cartesian components of J nabla psi
                !   = -J_V / (1+Lambda_theta) (T_C^V_^T) (e_V)
                allocate(norm_com_C(grid%n(1),3))
                norm_com_C = 0._dp
                allocate(T_CV_loc(size(eq_2%T_VC,1),1:1,1:1,&
                    &size(eq_2%T_VC,4),0:0,0:0,0:0))
                ierr = calc_inv_met(T_CV_loc,eq_2%T_VC(:,1:1,1:1,:,:,:,:),&
                    &[0,0,0])
                CHCKERR('')
                do id = 1,3
                    norm_com_C(:,id) = -eq_2%jac_E(:,1,r_id,0,0,0) * &
                        &T_CV_loc(:,1,1,c([id,1],.false.),0,0,0)                ! Cylindrical components
                end do
                deallocate(T_CV_loc)
                vac%norm(eq_job_lims(1):eq_job_lims(2),1) = &
                    &norm_com_C(:,1)*cos(grid%zeta_E(:,1,r_id)) - &
                    &norm_com_C(:,2)*sin(grid%zeta_E(:,1,r_id))                 ! ~ e_X
                vac%norm(eq_job_lims(1):eq_job_lims(2),2) = &
                    &norm_com_C(:,1)*sin(grid%zeta_E(:,1,r_id)) + &
                    &norm_com_C(:,2)*cos(grid%zeta_E(:,1,r_id))                 ! ~ e_Y
                vac%norm(eq_job_lims(1):eq_job_lims(2),3) = norm_com_C(:,3)     ! ~ e_Z
            end if
            
            ! broadcast to other processes
            do id = 1,size(vac%x_vec,2)
                ierr = broadcast_var(&
                    &vac%x_vec(eq_job_lims(1):eq_job_lims(2),id),n_procs-1)
                CHCKERR('')
                ierr = broadcast_var(&
                    &vac%norm(eq_job_lims(1):eq_job_lims(2),id),n_procs-1)
                CHCKERR('')
            end do
            
            call lvl_ud(-1)
            
            call writo('Done storing vacuum quantities')
        end function store_vac_VMEC
        ! HELENA version
        !> \private
        integer function store_vac_HEL() result(ierr)
            use HELENA_vars, only: R_H, Z_H, nchi, chi_H, ias
            use grid_utilities, only: nufft
            
            character(*), parameter :: rout_name = 'store_vac_HEL'
            
            ! local variables
            integer :: n_theta                                                  ! total number of points, including overlap
            real(dp), allocatable :: x_vec_F(:,:)                               ! Fourier components of x_vec
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Start storing vacuum quantities')
            
            call lvl_ud(1)
            
            ! if  restarting  a  Richardson  level, reconstruct.  If  not,  only
            ! calculate for the first level
            ! HELENA should only have 1 equilibrium job for Richardson level 1
            if (eq_job_nr.eq.1 .and. rich_lvl.eq.rich_restart_lvl) then
                if (rich_restart_lvl.eq.1) then                                 ! starting from zero
                    if (ias.eq.0) then                                          ! top-bottom symmmetric
                        n_theta = 2*(nchi-1)+1
                    else
                        n_theta = nchi
                    end if
                    
                    ! allocate
                    ierr = vac%init(eq_style,n_theta,prim_X,[n_theta,1],jq)
                    CHCKERR('')
                    ! save points
                    if (rank.eq.n_procs-1) then
                        ! grid point to save is last total
                        r_id = n_r_eq
                        
                        ! save only R and Z
                        if (ias.eq.0) then
                            vac%x_vec(1:nchi,1) = R_H(1:nchi,r_id)
                            vac%x_vec(1:nchi,2) = Z_H(1:nchi,r_id)
                            vac%ang(1:nchi,1) = chi_H(1:nchi)
                            vac%x_vec(nchi+1:n_theta,1) = R_H(nchi-1:1:-1,r_id)
                            vac%x_vec(nchi+1:n_theta,2) = -Z_H(nchi-1:1:-1,r_id)
                            vac%ang(nchi+1:n_theta,1) = chi_H(2:nchi) + pi
                        else
                            vac%x_vec(:,1) = R_H(:,r_id)
                            vac%x_vec(:,2) = Z_H(:,r_id)
                            vac%ang(:,1) = chi_H
                        end if
                        !call print_ex_2D(['R','Z'],'',vac%x_vec,&
                            !&x=vac%ang(:,1:1),persistent=.true.)
                        !call print_ex_2D('cs','',vac%x_vec(:,2),&
                            !&x=vac%x_vec(:,1),persistent=.true.)
                        
                        ! calculate Cylindrical components of - J nabla psi
                        !   = - R (Z_theta nabla R - R_theta nabla Z)
                        ! and its poloidal derivative
                        vac%norm = 0._dp
                        vac%dnorm = 0._dp
                        ! nufft of Z
                        ierr = nufft(vac%ang(1:n_theta-1,1),&
                            &vac%x_vec(1:n_theta-1,2),x_vec_F)
                        CHCKERR('')
                        ! dZ/dtheta
                        do id = 1,size(x_vec_F,1)
                            vac%norm(1:n_theta-1,1) = &
                                &vac%norm(1:n_theta-1,1) + (id-1) * &
                                &(-x_vec_F(id,1)*sin((id-1)*&
                                &vac%ang(1:n_theta-1,1))&
                                &+x_vec_F(id,2)*cos((id-1)*&
                                &vac%ang(1:n_theta-1,1)))
                        end do
                        ! nufft of -R
                        ierr = nufft(vac%ang(1:n_theta-1,1),&
                            &-vac%x_vec(1:n_theta-1,1),x_vec_F)
                        CHCKERR('')
                        ! d(-R)/dtheta
                        do id = 1,size(x_vec_F,1)
                            vac%norm(1:n_theta-1,2) = &
                                &vac%norm(1:n_theta-1,2) + (id-1) * &
                                &(-x_vec_F(id,1)*sin((id-1)*&
                                &vac%ang(1:n_theta-1,1))&
                                &+x_vec_F(id,2)*cos((id-1)*&
                                &vac%ang(1:n_theta-1,1)))
                        end do
                        vac%norm(n_theta,:) = vac%norm(1,:)
                        ! multiply by -R
                        do id = 1,2
                            vac%norm(:,id) = - vac%x_vec(:,1)*vac%norm(:,id)
                        end do
                        ! nufft of norm
                        do jd = 1,2
                            ierr = nufft(vac%ang(1:n_theta-1,1),&
                                &vac%norm(1:n_theta-1,jd),x_vec_F)
                            CHCKERR('')
                            do id = 1,size(x_vec_F,1)
                                vac%dnorm(1:n_theta-1,jd) = &
                                    &vac%dnorm(1:n_theta-1,jd) + (id-1) * &
                                    &(-x_vec_F(id,1)*sin((id-1)*&
                                    &vac%ang(1:n_theta-1,1))&
                                    &+x_vec_F(id,2)*cos((id-1)*&
                                    &vac%ang(1:n_theta-1,1)))
                            end do
                        end do
                        vac%dnorm(n_theta,:) = vac%dnorm(1,:)
                    end if
                else                                                            ! restarting
                    ! reconstruct old vacuum
                    ierr = reconstruct_PB3D_vac(vac,'vac')
                    CHCKERR('')
                end if
            end if
            
            ! broadcast to other processes
            do id = 1,2
                ierr = broadcast_var(vac%x_vec(:,id),n_procs-1)
                CHCKERR('')
                ierr = broadcast_var(vac%norm(:,id),n_procs-1)
                CHCKERR('')
                ierr = broadcast_var(vac%dnorm(:,id),n_procs-1)
                CHCKERR('')
            end do
            ierr = broadcast_var(vac%ang(:,1),n_procs-1)
            CHCKERR('')
            
            call lvl_ud(-1)
            
            call writo('Done storing vacuum quantities')
        end function store_vac_HEL
    end function store_vac
    
    !> Calculate the matrices \c G and \c H.
    !!
    !! This  is  done  by  iterating  over  the  subintervals,  calculating  the
    !! contributions to the  potential values on the left and  right boundary of
    !! the interval, and adding these to the appropriate values in G and H.
    !!
    !! As each process in  the blacs grid only has access to  its own values, an
    !! extra interval  is is calculated  to the left  and right of  the internal
    !! process grid edges.
    !!
    !! \return ierr
    integer function calc_GH(vac) result(ierr)
        use num_vars, only: rank
        use vac_utilities, only: vec_dis2loc, in_context
        use MPI_utilities, only: wait_MPI
#if ldebug
        use num_vars, only: n_procs
        use vac_utilities, only: mat_dis2loc
#endif
        
        character(*), parameter :: rout_name = 'calc_GH'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        integer :: id, jd                                                       ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        integer :: rd                                                           ! global counter for row
        integer :: rdl                                                          ! local counter for row
        integer :: i_rd                                                         ! index of subrow
        integer :: desc_phi(BLACSCTXTSIZE)                                      ! descriptor for phi
        integer :: desc_dphi(BLACSCTXTSIZE)                                     ! descriptor for dphi
        integer :: desc_res(BLACSCTXTSIZE)                                      ! descriptor for the result
        real(dp), allocatable :: loc_ser(:,:)                                   ! dummy serial variable
        real(dp), allocatable :: r_sph(:)                                       ! spherical r with r^2 = R^2+Z^2
        real(dp), allocatable :: phi(:,:)                                       ! test phi
        real(dp), allocatable :: dphi(:,:)                                      ! test dphi/dn
        real(dp), allocatable :: res(:,:,:)                                     ! H phi, G dphi
        character(len=max_str_ln) :: plot_title(3)                              ! plot title
        character(len=max_str_ln) :: plot_name                                  ! plot name
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate the G and H matrices')
        call lvl_ud(1)
        
        ! call specific procedure for vacuum style
        select case (vac%style)
            case (1)                                                            ! field-line 3-D
                ierr = calc_GH_1(vac)
                CHCKERR('')
            case (2)                                                            ! axisymmetric
                ierr = calc_GH_2(vac)
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No vacuum style '//trim(i2str(vac%style))//&
                    &' possible'
                CHCKERR(err_msg)
        end select
        
#if ldebug
        if (debug_calc_GH .and. in_context(vac%ctxt_HG)) then
            ! initialize
            call descinit(desc_phi,vac%n_bnd,1,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for phi')
            call descinit(desc_dphi,vac%n_bnd,1,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for dphi')
            call descinit(desc_res,vac%n_bnd,1,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for res')
            
            ! plot H and G in HDF5
            allocate(loc_ser(vac%n_bnd,vac%n_bnd))
            ierr = mat_dis2loc(vac%ctxt_HG,vac%H,vac%lims_r,vac%lims_c,loc_ser,&
                &proc=n_procs-1)
            CHCKERR('')
            if (rank.eq.n_procs-1) then
                call plot_HDF5('H','H',reshape(loc_ser,[vac%n_bnd,vac%n_bnd,1]))
                call plot_diff_HDF5(reshape(loc_ser,[vac%n_bnd,vac%n_bnd,1]),&
                    &reshape(transpose(loc_ser),[vac%n_bnd,vac%n_bnd,1]),&
                    &'H_transp')
            end if
            ierr = mat_dis2loc(vac%ctxt_HG,vac%G,vac%lims_r,vac%lims_c,loc_ser,&
                &proc=n_procs-1)
            CHCKERR('')
            if (rank.eq.n_procs-1) then
                call plot_HDF5('G','G',reshape(loc_ser,[vac%n_bnd,vac%n_bnd,1]))
                call plot_diff_HDF5(reshape(loc_ser,[vac%n_bnd,vac%n_bnd,1]),&
                    &reshape(transpose(loc_ser),[vac%n_bnd,vac%n_bnd,1]),&
                    &'G_transp')
            end if
            deallocate(loc_ser)
            
            ! test whether G and H hold for test potentials phi
            if (vac%lims_c(1,1).eq.1) then                                      ! this process owns (part of) first column
                ! set up variables
                allocate(res(vac%n_loc(1),3,2))
                allocate(r_sph(vac%n_loc(1)))
                allocate(phi(vac%n_loc(1),2))
                allocate(dphi(vac%n_loc(1),2))
                
                ! set up rhs for H (phi) and G (dphi)
                subrows: do i_rd = 1,size(vac%lims_r,2)
                    row: do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                        ! set local row index
                        rdl = sum(vac%lims_r(2,1:i_rd-1)-&
                            &vac%lims_r(1,1:i_rd-1)+1) + &
                            &rd-vac%lims_r(1,i_rd)+1
                        
                        ! spherical radius
                        r_sph(rdl) = sqrt(sum(vac%x_vec(rd,:)**2))
                        
                        ! test 1
                        phi(rdl,1) = vac%x_vec(rd,1)**vac%prim_X
                        dphi(rdl,1) = - vac%prim_X*vac%norm(rd,1)/&
                            &vac%x_vec(rd,1) * phi(rdl,1)
                        
                        ! test 2
                        phi(rdl,2) = (vac%x_vec(rd,1)/&
                            &(r_sph(rdl)+abs(vac%x_vec(rd,2))))**vac%prim_X
                        dphi(rdl,2) = vac%prim_X*(vac%norm(rd,2)-&
                            &vac%norm(rd,1)*vac%x_vec(rd,2)/vac%x_vec(rd,1)) / &
                            &r_sph(rdl)*sign(1._dp,vac%x_vec(rd,2)) * phi(rdl,2)
                    end do row
                end do subrows
            else
                allocate(res(0,3,2))
                allocate(r_sph(0))
                allocate(phi(0,2))
                allocate(dphi(0,2))
            end if
            
            ! output
            allocate(loc_ser(vac%n_bnd,3))
            ierr = vec_dis2loc(vac%ctxt_HG,r_sph,vac%lims_r,loc_ser(:,1),&
                &proc=n_procs-1)
            CHCKERR('')
            if (rank.eq.n_procs-1) then
                plot_title(1) = 'spherical r'
                plot_name = 'r_sph'
                call print_ex_2D(plot_title(1),plot_name,loc_ser(:,1),&
                    &x=vac%ang(:,1),draw=.false.)
                call draw_ex([plot_title],plot_name,1,1,.false.)
            end if
            do jd = 1,size(phi,2)
                ierr = vec_dis2loc(vac%ctxt_HG,phi(:,jd),vac%lims_r,&
                    &loc_ser(:,jd),proc=n_procs-1)
                CHCKERR('')
            end do
            if (rank.eq.n_procs-1) then
                plot_title(1) = 'test potential ϕ'
                plot_name = 'phi'
                call print_ex_2D([plot_title(1)],plot_name,&
                    &loc_ser(:,1:size(phi,2)),x=vac%ang(:,1:1),draw=.false.)
                call draw_ex([plot_title(1)],plot_name,size(phi,2),1,.false.)
            end if
            do jd = 1,size(dphi,2)
                ierr = vec_dis2loc(vac%ctxt_HG,dphi(:,jd),vac%lims_r,&
                    &loc_ser(:,jd),proc=n_procs-1)
                CHCKERR('')
            end do
            if (rank.eq.n_procs-1) then
                plot_title(1) = 'normal derivative of test potential dϕ/dn'
                plot_name = 'dphi'
                call print_ex_2D([plot_title(1)],plot_name,&
                    &loc_ser(:,1:size(dphi,2)),x=vac%ang(:,1:1),draw=.false.)
                call draw_ex([plot_title(1)],plot_name,size(dphi,2),1,.false.)
            end if
            deallocate(loc_ser)
            
            ! loop over test potentials
            do jd = 1,size(phi,2)
                ! calculate H phi
                call pdgemv('N',vac%n_bnd,vac%n_bnd,1._dp,vac%H,1,1,&
                    &vac%desc_H,phi(:,jd),1,1,desc_phi,1,0._dp,res(:,1,jd),&
                    &1,1,desc_res,1)
                
                ! calculate - G dphi
                call pdgemv('N',vac%n_bnd,vac%n_bnd,-1._dp,vac%G,1,1,&
                    &vac%desc_G,dphi(:,jd),1,1,desc_dphi,1,0._dp,res(:,2,jd),&
                    &1,1,desc_res,1)
                
                ! add them together
                res(:,3,jd) = sum(res(:,1:2,jd),2)
            end do
            
            ! output of tests
            allocate(loc_ser(vac%n_bnd,3))
            do id = 1,size(phi,2)
                do jd = 1,3
                    ierr = vec_dis2loc(vac%ctxt_HG,res(:,jd,id),vac%lims_r,&
                        &loc_ser(:,jd),proc=n_procs-1)
                    CHCKERR('')
                end do
                if (rank.eq.n_procs-1) then
                    select case (id)
                        case (1)
                            plot_title(1) = 'H R^n'
                            plot_title(2) = '-G n Z_θ R^n'
                            plot_title(3) = '(H - G n Z_θ) R^n'
                            plot_name = 'test_vac_1'
                        case (2)
                            plot_title(1) = 'H (R/(r+|Z|))^n'
                            plot_title(2) = '- G |Z|/Z n dr/dθ (R/(r+|Z|))^n'
                            plot_title(3) = '(H - G |Z|/Z n dr/dθ) &
                                &(R/(r+|Z|))^n'
                            plot_name = 'test_vac_2'
                    end select
                    call print_ex_2D(plot_title,trim(plot_name)//'_all',&
                        &loc_ser(:,:),x=vac%ang(:,1:1),&
                        &draw=.false.)
                    call draw_ex(plot_title,trim(plot_name)//'_all',3,1,.false.)
                    call print_ex_2D(plot_title(3),plot_name,&
                        &loc_ser(:,3),x=vac%ang(:,1),&
                        &draw=.false.)
                    call draw_ex(plot_title(3:3),plot_name,1,1,.false.)
                end if
            end do
            deallocate(loc_ser)
        end if
#endif
        
        call lvl_ud(-1)
    end function calc_GH

    !> \public Calculates matrices \c G and \c H in 3-D configuration.
    !!
    !! \see calc_gh.
    integer function calc_GH_1(vac) result(ierr)
        character(*), parameter :: rout_name = 'calc_GH_1'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! initialize ierr
        ierr = 0
        
        call writo('Using field-line 3-D Boundary Element Method')
        
    end function calc_GH_1
    
    !> \public Calculates matrices \c G and \c H in axisymmetric configurations.
    !!
    !! It makes use of toroidal functions.
    !!
    !! \see calc_gh.
    integer function calc_GH_2(vac) result(ierr)
        use dtorh, only: dtorh1
        use vac_utilities, only: calc_GH_int_2
        use MPI_utilities, only: get_ser_var
        use num_ops, only: calc_zero_HH, calc_zero_Zhang
        use num_vars, only: rank
        
        character(*), parameter :: rout_name = 'calc_GH_2'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: i_rd, i_cd                                                   ! index of subrow and column
        integer :: rd, cd                                                       ! global counters for row and column
        integer :: rdl, cdl                                                     ! local counters for row and column, used for G and H
        integer :: max_n                                                        ! maximum degree
        integer :: lims_c_loc(2)                                                ! local lims_c
        real(dp) :: b_coeff(2)                                                  ! sum_k=1^n 2/2k-1 for n and n-1
        real(dp) :: tol_loc                                                     ! local tolerance on rho for singular points
        real(dp) :: del_r(2,2)                                                  ! unit vectors of next and previous subinterval
        real(dp) :: G_loc(2)                                                    ! local G
        real(dp) :: H_loc(2)                                                    ! local H
        real(dp) :: beta_comp                                                   ! complementary beta
        real(dp), parameter :: tol_Q = 1.E-6_dp                                 ! tolerance on Q approximation
        real(dp), allocatable :: beta(:)                                        ! beta
        real(dp), allocatable :: loc_ser(:)                                     ! dummy serial variable
        real(dp), allocatable :: rho2(:)                                        ! square of distance in projected poloidal plane for one column
        real(dp), allocatable :: arg(:)                                         ! argument 1 + rho^2/(RiRj)(
        real(dp), allocatable :: Aij(:)                                         ! Aij helper variable for H
        real(dp), allocatable :: ql(:,:,:)                                      ! ql for all rho2's
        real(dp), allocatable :: pl_loc(:), ql_loc(:)                           ! local toroidal functions
        logical, allocatable :: was_calc(:,:)                                   ! was already calculated
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        real(dp) :: rho2_lims(2)                                                ! limits on rho2
        real(dp) :: gamma_lims(2)                                               ! limits on gamma
        integer :: beta_lims(2)                                                 ! limits to plot beta
        !integer :: lims_cl(2)                                                   ! vac%lims_c in local coordinates
        !character(len=max_str_ln) :: plot_title                                 ! plot title
        !character(len=max_str_ln) :: plot_name                                  ! plot name
#endif
        
        ! initialize ierr
        ierr = 0
        
        call writo('Using axisymmetric Boundary Element Method')
        
        ! initialize
        vac%G = 0._dp
        vac%H = 0._dp
        
        ! set up b_coeff and local tolerance
        b_coeff = 0._dp
        do kd = 1,vac%prim_X
            b_coeff(2) = b_coeff(2) + 2._dp/(2*kd-1)
        end do
        b_coeff(1) = b_coeff(2) - 2._dp/(2*vac%prim_X-1)
        err_msg = calc_zero_Zhang(tol_loc,fun_tol,&
            &[1.E-10_dp,&                                                       ! very low lower limit
            &exp(-2*b_coeff(2))])                                               ! upper limit for -1/2ln(x) = b
        if (trim(err_msg).ne.'') then
            ierr = 1
            CHCKERR(trim(err_msg))
        end if
        tol_loc = 8*minval(vac%x_vec(:,1))*tol_loc                              ! tolerance on rho
        ierr = get_ser_var([tol_loc],loc_ser,scatter=.true.)
        CHCKERR('')
        tol_loc = maxval(loc_ser)
        deallocate(loc_ser)
        call writo('Tolerance on rho for analytical integral: '//&
            &trim(r2str(tol_loc)))
        
        ! allocate helper variables
        allocate(was_calc(1:vac%n_bnd,1:vac%n_bnd))
        was_calc = .false.
        allocate(ql(1:vac%n_bnd,1:vac%n_bnd,2))
        allocate(pl_loc(0:max(vac%prim_X,1)))
        allocate(ql_loc(0:max(vac%prim_X,1)))
#if ldebug
        rho2_lims = [huge(1._dp),-huge(1._dp)]
        gamma_lims = [huge(1._dp),-huge(1._dp)]
#endif
        
        ! iterate over all subrows: position of singular point
        subrows: do i_rd = 1,size(vac%lims_r,2)
            ! initialize helper variables
            allocate(beta(vac%lims_r(1,i_rd):vac%lims_r(2,i_rd)))               ! global index
#if ldebug
            beta_lims = [vac%lims_r(2,i_rd)+1,vac%lims_r(1,i_rd)-1]
#endif
            
            ! iterate over all rows of this subrow
            row: do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                ! set local row index
                rdl = sum(vac%lims_r(2,1:i_rd-1)-vac%lims_r(1,1:i_rd-1)+1) + &
                    &rd-vac%lims_r(1,i_rd)+1
                
                ! iterate over all subcolumns: subintegrals
                subcols: do i_cd = 1,size(vac%lims_c,2)
                    ! set up local limits, including possible overlap
                    lims_c_loc(1) = max(1,vac%lims_c(1,i_cd)-1)
                    lims_c_loc(2) = min(vac%n_bnd,vac%lims_c(2,i_cd)+1)
                    
                    ! initialize helper variables
                    allocate(rho2(lims_c_loc(1):lims_c_loc(2)))                 ! global index
                    allocate(arg(lims_c_loc(1):lims_c_loc(2)))                  ! global index
                    allocate(Aij(lims_c_loc(1):lims_c_loc(2)))                  ! global index
                    
                    ! iterate over all columns  of this subcolumn, possibly with
                    ! an overlap left and/or right to set up the Aij
                    col: do cd = lims_c_loc(1),lims_c_loc(2)
                        ! check rho2 and set arg to calculate tor. functions
                        rho2(cd) = sum((vac%x_vec(cd,:)-vac%x_vec(rd,:))**2)
                        arg(cd) = 1._dp + &
                            &rho2(cd)/(2*vac%x_vec(cd,1)*vac%x_vec(rd,1))
                        if (rho2(cd).le.tol_loc**2) then
                            ! don't set it
                        else if (was_calc(cd,rd)) then                          ! check symmetric element
                            ql(rd,cd,:) = ql(cd,rd,:)                           ! set it
                        else                                                    ! not yet calculated by this process
#if ldebug
                            ! update extrema
                            if (rho2(cd).lt.rho2_lims(1)) &
                                &rho2_lims(1) = rho2(cd)
                            if (rho2(cd).gt.rho2_lims(2)) &
                                &rho2_lims(2) = rho2(cd)
                            if (arg(cd).lt.gamma_lims(1)) &
                                &gamma_lims(1) = arg(cd)
                            if (arg(cd).gt.gamma_lims(2)) &
                                &gamma_lims(2) = arg(cd)
#endif
                            err_msg = 'The solutions did not converge for &
                                &rho^2 = '//trim(r2str(rho2(cd)))//&
                                &' and n = '//trim(i2str(vac%prim_X))
                            if (vac%prim_X.gt.0) then                           ! positive n
                                ierr = dtorh1(arg(cd),0,vac%prim_X,pl_loc,&
                                    &ql_loc,max_n)
                                CHCKERR('')
                                if (max_n.lt.vac%prim_X) then
                                    ierr = 1
                                    CHCKERR(err_msg)
                                end if
                                ql(rd,cd,:) = ql_loc(vac%prim_X-1:vac%prim_X)
                            else if (vac%prim_X.eq.0) then
                                ierr = dtorh1(arg(cd),0,1,pl_loc,ql_loc,max_n)
                                CHCKERR('')
                                if (max_n.lt.vac%prim_X) then
                                    ierr = 1
                                    CHCKERR(err_msg)
                                end if
                                ql(rd,cd,1) = ql_loc(1)                         ! Q_-3/2 = Q_1/2
                                ql(rd,cd,2) = ql_loc(0)                         ! Q_-1/2
                            else                                                ! negative n
                                ierr = 1
                                err_msg = 'negative n not yet implemented'      ! SHOULD JUST BE Q_{n-1/2} = Q_{-n-1/2}
                                CHCKERR('')
                            end if
                            was_calc(rd,cd) = .true.
                        end if
                        
                        ! set up Aij
                        if (rho2(cd).le.tol_loc**2) then
                            Aij(cd) = 0.5_dp*(vac%prim_X - 0.5_dp) * &
                                &(vac%x_vec(rd,1) * &
                                &(vac%norm(rd,1)*vac%dnorm(rd,2)-&
                                &vac%norm(rd,2)*vac%dnorm(rd,1))/&
                                &sum(vac%norm(rd,:)**2) &
                                &+vac%norm(rd,1)/vac%x_vec(rd,1))
                        else
                            Aij(cd) = 2._dp*vac%x_vec(rd,1)*vac%x_vec(cd,1) / &
                                &(4._dp*vac%x_vec(rd,1)*vac%x_vec(cd,1)+&
                                &rho2(cd)) *(vac%prim_X - 0.5_dp) * &
                                &(- 2._dp/rho2(cd) * sum(vac%norm(cd,:)*&
                                &(vac%x_vec(cd,:)-vac%x_vec(rd,:))) + &
                                &vac%norm(cd,1)/vac%x_vec(cd,1))
                        end if
                    end do col
                    
#if ldebug
                    if (debug_calc_GH) then
                        !! output
                        !plot_title = 'for row '//trim(i2str(rd))//&
                            !&' and starting column '//trim(i2str(lims_c_loc(1)))
                        !plot_name = '_'//trim(i2str(rd))//'_'//&
                            !&trim(i2str(lims_c_loc(1)))
                        !call print_ex_2D('rho^2 '//trim(plot_title),&
                            !&'rho2'//trim(plot_name),rho2,&
                            !&x=vac%ang(lims_c_loc(1):lims_c_loc(2),1),&
                            !&draw=.false.)
                        !call draw_ex(['rho^2 '//trim(plot_title)],&
                            !&'rho2'//trim(plot_name),1,1,.false.)
                        !call print_ex_2D('γ '//trim(plot_title),&
                            !&'gamma'//trim(plot_name),arg,&
                            !&x=vac%ang(lims_c_loc(1):lims_c_loc(2),1),&
                            !&draw=.false.)
                        !call draw_ex(['γ '//trim(plot_title)],&
                            !&'gamma'//trim(plot_name),1,1,.false.)
                        !call print_ex_2D('A_{ij} '//trim(plot_title),&
                            !&'Aij'//trim(plot_name),arg,&
                            !&x=vac%ang(lims_c_loc(1):lims_c_loc(2),1),&
                            !&draw=.false.)
                        !call draw_ex(['A_{ij} '//trim(plot_title)],&
                            !&'Aij'//trim(plot_name),1,1,.false.)
                    end if
#endif
                    
                    ! iterate over all pairs of columns of this subcolumn
                    ! (upper limit is one lower as previous)
                    col2: do cd = lims_c_loc(1),lims_c_loc(2)-1
                        ! set local column index
                        cdl = sum(vac%lims_c(2,1:i_cd-1)-&
                            &vac%lims_c(1,1:i_cd-1)+1) + &
                            &cd-vac%lims_c(1,i_cd)+1
                        
                        ! calculate contributions to H and G
                        call calc_GH_int_2(G_loc,H_loc,&
                            &vac%ang(cd:cd+1,1),vac%ang(rd,1),&
                            &vac%x_vec(cd:cd+1,:),vac%x_vec(rd,:),&
                            &vac%norm(cd:cd+1,:),vac%norm(rd,:),&
                            &Aij(cd:cd+1),ql(rd,cd:cd+1,:),&
                            &tol_loc**2,b_coeff,vac%prim_X)
                        
                        do kd = 0,1
                            if (cd+kd.ge.vac%lims_c(1,i_cd) .and. &
                                &cd+kd.le.vac%lims_c(2,i_cd)) then
                                vac%G(rdl,cdl+kd) = vac%G(rdl,cdl+kd) + G_loc(1+kd)
                                vac%H(rdl,cdl+kd) = vac%H(rdl,cdl+kd) + H_loc(1+kd)
                            end if
                        end do
                    end do col2
                    
#if ldebug
                    if (debug_calc_GH) then
                        !! local limits
                        !lims_cl = sum(vac%lims_c(2,1:i_cd-1)-&
                            !&vac%lims_c(1,1:i_cd-1)+1) + &
                            !&vac%lims_c(:,i_cd)-vac%lims_c(1,i_cd)+1
                        !! output
                        !call print_ex_2D('G '//trim(plot_title),&
                            !&'G'//trim(plot_name),&
                            !&vac%G(rdl,lims_cl(1):lims_cl(2)),x&
                            !&=vac%ang(vac%lims_c(1,i_cd):vac%lims_c(2,i_cd),1),&
                            !&draw=.false.)
                        !call draw_ex(['G '//trim(plot_title)],&
                            !&'G'//trim(plot_name),1,1,.false.)
                        !call print_ex_2D('H '//trim(plot_title),&
                            !&'H'//trim(plot_name),&
                            !&vac%H(rdl,lims_cl(1):lims_cl(2)),x&
                            !&=vac%ang(vac%lims_c(1,i_cd):vac%lims_c(2,i_cd),1),&
                            !&draw=.false.)
                        !call draw_ex(['H '//trim(plot_title)],&
                            !&'H'//trim(plot_name),1,1,.false.)
                    end if
#endif
                    
                    ! calculate contribution due to fundamental solution
                    if (vac%lims_c(1,i_cd).le.rd .and. &
                        &rd.le.vac%lims_c(2,i_cd)) then                         ! this subcolumn includes the diagonal
                        ! set local column index
                        cdl = sum(vac%lims_c(2,1:i_cd-1)-&
                            &vac%lims_c(1,1:i_cd-1)+1) + &
                            &rd-vac%lims_c(1,i_cd)+1
                        
                        if (rd.eq.1 .or. rd.eq.vac%n_bnd) then                  ! first or last global point
                            ! loop around to get previous direction vector
                            del_r(1,:) = vac%x_vec(vac%n_bnd,:)-&
                                &vac%x_vec(vac%n_bnd-1,:)
                            del_r(2,:) = vac%x_vec(2,:)-vac%x_vec(1,:)
                        else                                                    ! all points in between
                            del_r(1,:) = vac%x_vec(rd,:)-vac%x_vec(rd-1,:)
                            del_r(2,:) = vac%x_vec(rd+1,:)-vac%x_vec(rd,:)
                        end if
                        do kd = 1,2
                            del_r(kd,:) = del_r(kd,:)/sqrt(sum(del_r(kd,:)**2))
                        end do
                        beta_comp = acos(sum(product(del_r,1)))
                        if (del_r(1,1)*del_r(2,2) .lt. &
                            &del_r(1,2)*del_r(2,1)) beta_comp = - beta_comp     ! concave
                        beta(rd) = beta_comp + pi
#if ldebug
                        if (debug_calc_GH) then
                            if (rd.lt.beta_lims(1)) beta_lims(1) = rd
                            if (rd.gt.beta_lims(2)) beta_lims(2) = rd
                        end if
#endif
                        
                        vac%H(rdl,cdl) = vac%H(rdl,cdl) - 2*beta(rd)
                    end if
                    
                    ! clean up
                    deallocate(arg,rho2,Aij)
                end do subcols
            end do row
            
#if ldebug
            if (debug_calc_GH) then
                !! output
                !if (beta_lims(2).ge.beta_lims(1)) then
                    !plot_title = 'beta/π  for starting row '//&
                        !&trim(i2str(beta_lims(1)))
                    !plot_name = 'beta_'//trim(i2str(beta_lims(1)))
                    !call print_ex_2D(plot_title,plot_name,&
                        !&beta(beta_lims(1):beta_lims(2))/pi,&
                        !&x=vac%ang(beta_lims(1):beta_lims(2),1),draw=.false.)
                    !call draw_ex([plot_title],plot_name,1,1,.false.)
                !end if
            end if
#endif
            
            ! clean up
            deallocate(beta)
        end do subrows
        
#if ldebug
        if (debug_calc_GH) then
            ! output
            ierr = get_ser_var([rho2_lims(1)],loc_ser)
            CHCKERR('')
            if (rank.eq.0) rho2_lims(1) = minval(loc_ser)
            deallocate(loc_ser)
            ierr = get_ser_var([rho2_lims(2)],loc_ser)
            CHCKERR('')
            if (rank.eq.0) rho2_lims(2) = maxval(loc_ser)
            deallocate(loc_ser)
            ierr = get_ser_var([gamma_lims(1)],loc_ser)
            CHCKERR('')
            if (rank.eq.0) gamma_lims(1) = minval(loc_ser)
            deallocate(loc_ser)
            ierr = get_ser_var([gamma_lims(2)],loc_ser)
            CHCKERR('')
            if (rank.eq.0) gamma_lims(2) = maxval(loc_ser)
            deallocate(loc_ser)
            if (rank.eq.0) then
                call writo('limits in rho: '//trim(r2str(sqrt(rho2_lims(1))))//&
                    &'..'//trim(r2str(sqrt(rho2_lims(2)))))
                call writo('limits in gamma: '//trim(r2str(gamma_lims(1)))//&
                    &'..'//trim(r2str(gamma_lims(2))))
            end if
        end if
#endif
    contains
        ! function that returns f = -1/2 ln(x)-b - x/T.
        !
        ! It uses b_coeff (b) and tol_Q (T)  from the parent function where T is
        ! a hard-coded positive value << 1.
        ! 
        !> \private
        function fun_tol(val)
            ! input / output
            real(dp), intent(in) :: val
            real(dp) :: fun_tol
            
            fun_tol = -0.5_dp*log(val)-b_coeff(2)-val/tol_Q
        end function fun_tol
    end function calc_GH_2
    
    !> Calculates the vacuum response.
    !!
    !! First  \c G and \c  H are completed if  this is not the  first Richardson
    !! level.
    !! 
    !! Then the vacuum contribution is calculated using the two-fold strategy of
    !! first solving 
    !!  \f$\overline{\text{H}}\overline{\text{C}} =
    !!  \overline{\text{G}}\overline{\text{E}} \overline{\text{P}}\f$,
    !! followed by left-multiplication of \f$\overline{\text{C}}\f$ by
    !!  \f$\overline{\text{P}}\overline{\text{E}}^\dagger\f$.
    !!
    !! The  non-square matrix  \f$\overline{\text{E}}\f$ contains  the exponents
    !! for different  mode numbers and  different poloidal grid  points, whereas
    !! \f$\overline{\text{P}}\f$ are diagonal matrices  that contain the factors
    !! \f$(nq-m)\f$.
    !!
    !! In practice, the complex matrix \f$\overline{\text{E}}\f$ is split in the
    !! two components of a real matrix twice the width.
    !!
    !! If \c  jump_to_sol is used for  the current Richardson level,  the vacuum
    !! quantities are not calculated, but just restored.
    !!
    !! \return ierr
    integer function calc_vac(vac) result(ierr)
        use MPI
        use rich_vars, only: rich_lvl
        use num_vars, only: jump_to_sol, rich_restart_lvl, eq_style, &
            &use_pol_flux_F, rank, n_procs
        use PB3D_ops, only: reconstruct_PB3D_vac
        use X_vars, only: n_mod_X, n_X, m_X
        use vac_vars, only: set_loc_lims
        use vac_utilities, only: in_context, mat_dis2loc
        use MPI_utilities, only: wait_MPI
        use eq_vars, only: vac_perm
        use num_utilities, only: LCM
        
        character(*), parameter :: rout_name = 'calc_vac'
        
        ! input / output
        type(vac_type), intent(inout), target :: vac                            !< vacuum variables
        
        ! local variables
        integer :: jd                                                           ! counter
        integer :: rd, cd                                                       ! global counter for row and col
        integer :: rdl, cdl                                                     ! local counter for row and col
        integer :: i_rd, i_cd                                                   ! index of subrow and subcol
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        integer :: sec_X_loc                                                    ! local sec_X
        integer :: n_loc(2)                                                     ! n_loc for n_mod_X instead of n_bnd
        integer :: lims_cl(2)                                                   ! lims_c_loc in local coordinates
        integer, allocatable :: lims_r_loc(:,:)                                 ! local row limits for n_loc
        integer, allocatable :: lims_c_loc(:,:)                                 ! local column limits for n_loc
        integer, target :: desc_EP(BLACSCTXTSIZE)                               ! descriptor for EP
        integer, target :: desc_GEP(BLACSCTXTSIZE)                              ! descriptor for GEP
        integer, target :: desc_Phi(BLACSCTXTSIZE)                              ! descriptor for Phi
        integer, target :: desc_res2(BLACSCTXTSIZE)                             ! descriptor for res2
        real(dp) :: arg_loc                                                     ! local argument
        real(dp), allocatable, target :: EP(:,:)                                ! EP
        real(dp), allocatable, target :: GEP(:,:)                               ! GEP
        real(dp), allocatable, target :: Phi(:,:)                               ! intermediate Phi
        real(dp), allocatable, target :: res2(:,:)                              ! vacuum response, real version
        real(dp), allocatable, target :: res2_loc(:,:)                          ! local version of res2
        type(StrumpackDensePackage_F90_double) :: SDP_loc                       ! Strumpack object for local solve H Phi = GEP
        type(C_PTR) :: CdescH, CH                                               ! C pointers to descriptor of H and H itself
        type(C_PTR) :: CdescGEP, CGEP                                           ! C pointers to descriptor of GEP and GEP itself
        type(C_PTR) :: CdescPhi, CPhi                                           ! C pointers to descriptor of Phi and Phi itself
        real(dp) :: SDP_err                                                     ! error
        logical :: rank_o                                                       ! output rank
#if ldebug
        integer :: lwork                                                        ! size of work array
        real(dp), allocatable :: res2_ev(:,:)                                   ! eigenvalues of res2
        real(dp), allocatable :: tau(:)                                         ! array tau in upper Hessenberg form
        real(dp), allocatable :: work(:)                                        ! work array
        complex(dp), allocatable :: res_loc(:,:)                                ! local vac%res
        complex(dp), allocatable :: res_ev(:)                                   ! eigenvalues of res2
        complex(dp), allocatable :: work_c(:)                                   ! work array
        complex(dp), allocatable :: tau_c(:)                                    ! array tau in upper Hessenberg form
        !integer :: lims_rl(2)                                                   ! vac%lims_r in local coordinates
        !character(len=max_str_ln) :: plot_title                                 ! plot title
        !character(len=max_str_ln) :: plot_name                                  ! plot name
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! decide whether to append Richardson level to name
        select case (eq_style)
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        if (rich_lvl.eq.rich_restart_lvl .and. jump_to_sol) then
            ! reconstruct old vacuum
            ierr = reconstruct_PB3D_vac(vac,'vac',rich_lvl=rich_lvl_name)
            CHCKERR('')
        end if
        
        call writo('Start calculation of vacuum')
        
        call lvl_ud(1)
        
        ! set ouptut rank persistency
        if (rank.eq.n_procs-1) then
            rank_o = .true.
        else
            rank_o = .false.
        end if
        
        ! calculate matrices G and H
        ierr = calc_GH(vac)
        CHCKERR('')
        
        ! Only for processes that are in the blacs context
        if (in_context(vac%ctxt_HG)) then
            ! user output
            call writo('Launch STRUMPack',persistent=rank_o)
            call lvl_ud(1)
            
            ! set secondary mode numbers
            if (use_pol_flux_F) then
                vac%lim_sec_X(1) = m_X(size(m_X,1),1)
                vac%lim_sec_X(2) = m_X(size(m_X,1),size(m_X,2))
            else
                vac%lim_sec_X(1) = n_X(size(n_X,1),1)
                vac%lim_sec_X(2) = n_X(size(n_X,1),size(n_X,2))
            end if
            
            ! set C pointers
            CH = C_LOC(vac%H)
            CdescH = C_LOC(vac%desc_H)
            
            ! initialize the solver and set parameters
            call SDP_F90_double_init(SDP_loc,vac%MPI_Comm)
            SDP_loc%use_HSS=1
            SDP_loc%levels_HSS=4
            SDP_loc%min_rand_HSS=10
            SDP_loc%lim_rand_HSS=5
            SDP_loc%inc_rand_HSS=10
            SDP_loc%max_rand_HSS=100
            SDP_loc%tol_HSS=1D-12
            SDP_loc%steps_IR=10
            SDP_loc%tol_IR=1D-10
            
            ! compression
            call SDP_F90_double_compress_A(SDP_loc,CH,CdescH)
            
#if ldebug
            ! accuracy checking
            SDP_err = SDP_F90_double_check_compression(SDP_loc,CH,CdescH)
#endif
            
            ! factorization
            call SDP_F90_double_factor(SDP_loc,CH,CdescH)
            
            ! set sizes
            do jd = 1,2
                n_loc(jd) = numroc(2*n_mod_X,vac%bs,vac%ind_p(jd),0,vac%n_p(jd))
            end do
            call set_loc_lims(n_loc(1),vac%bs,vac%ind_p(1),vac%n_p(1),&
                &lims_r_loc)
            call set_loc_lims(n_loc(2),vac%bs,vac%ind_p(2),vac%n_p(2),&
                &lims_c_loc)
            
            ! initialize descriptors
            call descinit(desc_EP,vac%n_bnd,2*n_mod_X,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for EP')
            call descinit(desc_GEP,vac%n_bnd,2*n_mod_X,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for GEP')
            call descinit(desc_Phi,vac%n_bnd,2*n_mod_X,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for Phi')
            call descinit(desc_res2,2*n_mod_X,2*n_mod_X,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,n_loc(1)),ierr)
            CHCKERR('descinit failed for res2')
            
            ! allocate auxililary matrices
            allocate(EP(vac%n_loc(1),n_loc(2)))
            allocate(GEP(vac%n_loc(1),n_loc(2)))
            allocate(Phi(vac%n_loc(1),n_loc(2)))
            allocate(res2(n_loc(1),n_loc(2)))
            allocate(res2_loc(2*n_mod_X,2*n_mod_X))
            
            ! set EP
            ! iterate over all subcolumns of matrix with 2n_mod_X columns
            subcols: do i_cd = 1,size(lims_c_loc,2)
                col: do cd = lims_c_loc(1,i_cd),lims_c_loc(2,i_cd)
                    ! set local column index
                    cdl = sum(lims_c_loc(2,1:i_cd-1)-&
                        &lims_c_loc(1,1:i_cd-1)+1) + &
                        &cd-lims_c_loc(1,i_cd)+1
                    
                    ! set local secondary mode number and angle
                    sec_X_loc = vac%lim_sec_X(1)-1+(mod(cd-1,n_mod_X)+1)
                    arg_loc = arg(vac%prim_X,sec_X_loc)
                    
                    ! iterate over all subrows of matrix with n_bnd rows
                    subrows: do i_rd = 1,size(vac%lims_r,2)
                        row: do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                            ! set local row index
                            rdl = sum(vac%lims_r(2,1:i_rd-1)-&
                                &vac%lims_r(1,1:i_rd-1)+1) + &
                                &rd-vac%lims_r(1,i_rd)+1
                                
                            ! save EP
                            if (cd.le.n_mod_X) then                             ! real part of rhs
                                EP(rdl,cdl) = cos(vac%ang(rd,1)*arg_loc)        ! rp(i exp)
                            else
                                EP(rdl,cdl) = sin(vac%ang(rd,1)*arg_loc)        ! ip(i exp)
                            end if
                            EP(rdl,cdl) = (vac%prim_X*vac%jq-sec_X_loc)*&
                                &EP(rdl,cdl)
                        end do row
                    end do subrows
                end do col
            end do subcols
            
            ! calculate GEP
            call pdgemm('N','N',vac%n_bnd,2*n_mod_X,vac%n_bnd,1._dp,vac%G,1,1,&
                &vac%desc_G,EP,1,1,desc_EP,0._dp,GEP,1,1,desc_GEP)
            
            ! set C pointers
            CGEP = C_LOC(GEP)
            CdescGEP = C_LOC(desc_GEP)
            CPhi = C_LOC(Phi)
            CdescPhi = C_LOC(desc_Phi)
            
            ! solve H c = G EP
            call SDP_F90_double_solve(SDP_loc,CPhi,CdescPhi,CGEP,CdescGEP)
            
#if ldebug
            ! accuracy checking
            SDP_err = SDP_F90_double_check_solution(SDP_loc,CH,CdescH,CPhi,&
                &CdescPhi,CGEP,CdescGEP)
#endif
            
            ! iterative refinement
            SDP_err = SDP_F90_double_refine(SDP_loc,CH,CdescH,CPhi,CdescPhi,&
                &CGEP,CdescGEP)
            
#if ldebug
            ! statistics
            call SDP_F90_double_print_statistics(SDP_loc)
            
            !! user output
            !subcols2: do i_cd = 1,size(lims_c_loc,2)
                !col2: do cd = lims_c_loc(1,i_cd),lims_c_loc(2,i_cd)
                    !! set local column index
                    !cdl = sum(lims_c_loc(2,1:i_cd-1)-&
                        !&lims_c_loc(1,1:i_cd-1)+1) + &
                        !&cd-lims_c_loc(1,i_cd)+1
                    
                    !subrows2: do i_rd = 1,size(vac%lims_r,2)
                        !! set local row limits
                        !lims_rl = sum(vac%lims_r(2,1:i_rd-1)-&
                            !&vac%lims_r(1,1:i_rd-1)+1) + &
                            !&vac%lims_r(:,i_rd)-vac%lims_r(1,i_rd)+1
                        
                        !! output
                        !plot_title = 'for column '//trim(i2str(cd))//&
                            !&' and starting row '//&
                            !&trim(i2str(vac%lims_r(1,i_rd)))
                        !plot_name = '_'//trim(i2str(cd))//'_'//&
                            !&trim(i2str(vac%lims_r(1,i_rd)))
                        !call print_ex_2D('Phi '//trim(plot_title),&
                            !&'Phi'//trim(plot_name),&
                            !&Phi(lims_rl(1):lims_rl(2),cdl),&
                            !&x=vac%ang(vac%lims_r(1,i_rd):&
                            !&vac%lims_r(2,i_rd),1),draw=.false.)
                        !call draw_ex(['Phi '//trim(plot_title)],&
                            !&'Phi'//trim(plot_name),1,1,.false.)
                    !end do subrows2
                !end do col2
            !end do subcols2
#endif
            
            ! convert EP into IEP where I contains the integration rule:
            !         (θ_{i+1}-θ_{i-1})/2   for i = 2..n-1
            !   I_i = (θ_2-θ_1}/2           for i = 1
            !         (θ_{n}-θ_{n-1}}/2     for i = n
            ! with n = vac%n_bnd.
            subrows3: do i_rd = 1,size(vac%lims_r,2)
                row3: do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                    ! set local row index
                    rdl = sum(vac%lims_r(2,1:i_rd-1)-&
                        &vac%lims_r(1,1:i_rd-1)+1) + &
                        &rd-vac%lims_r(1,i_rd)+1
                    
                    subcols3: do i_cd = 1,size(lims_c_loc,2)
                        ! set local column limits
                        lims_cl = sum(lims_c_loc(2,1:i_cd-1)-&
                            &lims_c_loc(1,1:i_cd-1)+1) + &
                            &lims_c_loc(:,i_cd)-lims_c_loc(1,i_cd)+1
                        
                        EP(rdl,lims_cl(1):lims_cl(2)) = &
                            &EP(rdl,lims_cl(1):lims_cl(2)) * 0.5_dp*&
                            &(vac%ang(min(rd+1,vac%n_bnd),1)-&
                            &vac%ang(max(rd-1,1),1))
                    end do subcols3
                end do row3
            end do subrows3
            
            ! calculate (EP)^T Phi
            ! Note that this means that the second  half of the rows in EP^T are
            ! still missing a  minus sign if they are to  represent EP^†. This
            ! will be taken into account when res2_loc is converted to vac%res.
            call pdgemm('T','N',2*n_mod_X,2*n_mod_X,vac%n_bnd,1._dp,EP,1,1,&
                &desc_EP,Phi,1,1,desc_Phi,0._dp,res2,1,1,desc_res2)
            
            ! get serial version on last process
            ierr = mat_dis2loc(vac%ctxt_HG,res2,lims_r_loc,lims_c_loc,res2_loc,&
                &proc=n_procs-1)
            CHCKERR('')
            
#if ldebug
            ! user output
            if (debug_calc_vac) then
                call writo('Calculate the eigenvalues of res_loc to &
                    &investigate whether it is positive definite',&
                    &persistent=rank_o)
                call lvl_ud(1)
            end if
#endif
            
            ! convert to complex numbers and save in vac%res
            if (rank.eq.n_procs-1) then
                vac%res = res2_loc(1:n_mod_X,1:n_mod_X)-&
                    &(-res2_loc(n_mod_X+1:2*n_mod_X,n_mod_X+1:2*n_mod_X))       ! minus sign for lower rows
                vac%res = vac%res + iu * &
                    &(res2_loc(1:n_mod_X,n_mod_X+1:2*n_mod_X)-&
                    &res2_loc(n_mod_X+1:2*n_mod_X,1:n_mod_X))                   ! minus sign for lower rows
                vac%res = vac%res/vac_perm
#if ldebug
                call print_ex_2D(['real vacuum response'],'vac_res_Re',&
                    &rp(vac%res),draw=.false.)
                call draw_ex(['real vacuum response'],'vac_res_Re',n_mod_X,1,&
                    &.false.)
                call print_ex_2D(['imag vacuum response'],'vac_res_Im',&
                    &ip(vac%res),draw=.false.)
                call draw_ex(['imag vacuum response'],'vac_res_Im',n_mod_X,1,&
                    &.false.)
                call plot_HDF5(['real part','imag part','imag diff'],'vac_res',&
                    &reshape([rp(vac%res),ip(vac%res),&
                    &ip(vac%res+transpose(vac%res))],[n_mod_X,n_mod_X,1,3]))
                
                if (debug_calc_vac) then
                    ! copy vacuum response
                    allocate(res_loc(n_mod_X,n_mod_X))
                    res_loc = vac%res
                    
                    ! get size of work matrix required for zgehrd
                    allocate(work_c(1))
                    allocate(tau_c(n_mod_X-1))
                    call zgehrd(n_mod_X,1,n_mod_X,res_loc,n_mod_X,tau_c,work_c,&
                        &-1,ierr)
                    CHCKERR('zgehrd failed')
                    lwork = nint(rp(work_c(1)))
                    deallocate(work_c)
                    allocate(work_c(lwork))
                    
                    ! convert to upper Heisenberg form
                    call zgehrd(n_mod_X,1,n_mod_X,res_loc,n_mod_X,tau_c,work_c,&
                        &lwork,ierr)
                    CHCKERR('zgehrd failed')
                    deallocate(work_c)
                    
                    ! get size of work matrix required for zhseqr
                    allocate(work_c(1))
                    allocate(res_ev(n_mod_X))
                    call zhseqr('E','N',n_mod_X,1,n_mod_X,res_loc,n_mod_X,&
                        &res_ev,0._dp,n_mod_X,work_c,-1,ierr)
                    CHCKERR('zhseqr failed')
                    lwork = nint(rp(work_c(1)))
                    deallocate(work_c)
                    allocate(work_c(lwork))
                    
                    ! calculate eigenvalues
                    call zhseqr('E','N',n_mod_X,1,n_mod_X,res_loc,n_mod_X,&
                        &res_ev,0._dp,n_mod_X,work_c,lwork,ierr)
                    CHCKERR('zhseqr failed')
                    call print_ex_2D(['real part','imag part'],'res_ev',&
                        &reshape([rp(res_ev),ip(res_ev)],[n_mod_X,2]),&
                        &draw=.false.)
                    call draw_ex(['real part','imag part'],'res_ev',2,1,&
                        &.false.)
                    deallocate(work_c)
                    deallocate(res_ev)
                end if
#endif
            end if
#if ldebug
            
            ! calculate the Eigenvalues of the original real res2
            if (debug_calc_vac) then
                ! get size of work matrix required for pdgehrd
                allocate(work(1))
                allocate(tau(2*n_mod_X-1))
                call pdgehrd(2*n_mod_X,1,2*n_mod_X,res2,1,1,desc_res2,tau,work,&
                    &-1,ierr)
                CHCKERR('pdgehrd failed')
                lwork = nint(work(1))
                deallocate(work)
                allocate(work(lwork))
                
                ! convert to upper Heisenberg form
                call pdgehrd(2*n_mod_X,1,2*n_mod_X,res2,1,1,desc_res2,tau,&
                    &work,lwork,ierr)
                CHCKERR('pdgehrd failed')
                deallocate(work)
                deallocate(tau)
                
                ! calculate eigenvalues
                lwork = 2*n_mod_X / vac%bs
                if (lwork*vac%bs.lt.2*n_mod_X) lwork = lwork + 1
                lwork = 7*lwork / LCM(vac%n_p(1),vac%n_p(2))
                lwork = 3*2*n_mod_X + max(2*2*n_mod_X+2*n_loc(2) , lwork)
                allocate(work(lwork))
                allocate(res2_ev(2*n_mod_X,2))
                call pdlahqr(.false.,.false.,2*n_mod_X,1,2*n_mod_X,res2,&
                    &desc_res2,res2_ev(:,1),res2_ev(:,2),1,2*n_mod_X,[0._dp],&
                    &desc_res2,work,lwork,[0],1,ierr)
                CHCKERR('pdlahqr failed')
                call print_ex_2D(['real part','imag part'],'res2_ev',&
                    &res2_ev,draw=.false.)
                call draw_ex(['real part','imag part'],'res2_ev',2,1,.false.)
                deallocate(work)
                deallocate(res2_ev)
                
                ! user output
                call lvl_ud(-1)
                call writo('They should all be positive',persistent=rank_o)
            end if
#endif
        end if
        
        call lvl_ud(-1)
        
        call writo('Done calculating vacuum quantities')
    contains
        ! Returns the argument to be used in EP
        !>  \private
        real(dp) function arg(n,m)
            
            ! input / output
            integer, intent(in) :: n, m                                         ! toroidal and poloidal mode number
            
            ! local variables
            real(dp) :: fac_n                                                   ! factor to multiply with n
            real(dp) :: fac_m                                                   ! factor to multiply with m
            
            ! set factor for n and m
            select case (vac%style)
                case (1)                                                        ! field-line 3-D
                    if (use_pol_flux_F) then
                        fac_n = vac%jq
                        fac_m = 1.0_dp
                    else
                        fac_n = 1.0_dp
                        fac_m = vac%jq
                    end if
                case (2)                                                        ! axisymmetric
                    fac_n = 0._dp
                    fac_m = 1._dp
            end select
            
            ! set argument
            arg = n*fac_n - m*fac_m
        end function arg
    end function calc_vac
    
    !> Print vacuum quantities of a certain order to an output file
    !!
    !!  - vac:
    !!      - \c misc_vac
    !!      - \c norm
    !!      - \c x_vec
    !!      - \c vac_res
    !!
    !! If \c  rich_lvl is provided, \c  "_R_[rich_lvl]" is appended to  the data
    !! name if it is \c >0.
    !!
    !! \note
    !!  -#  This routine  should  be called  by a  single  process, in  contrast
    !!  to  the   other  output   routines  such   as  eq_ops.print_output_eq(),
    !!  print_output_sol(), ...
    !!  -# This  process should  be the last  one, as it  will set  the boundary
    !!  contribution.
    !!  -# This routine  does not save the  \c H and \c G  matrices because they
    !!  are large and foreseen to be  saved directly as Hierarchical matrices in
    !!  the future. They therefore need  to be regenerated if Richardson restart
    !!  is done, or  when after a jump to solution,  another Richardson level is
    !!  attempted. In calc_vac, it is  checked when copying the vacuum variables
    !!  from previous  to current Richardson  level, whether the  \c G and  \c H
    !!  matrices are allocated and if not, they are calculated.
    !!
    !! \return ierr
    integer function print_output_vac(vac,data_name,rich_lvl) result(ierr)
        use X_vars, only: n_mod_X
        use num_vars, only: PB3D_name
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        
        character(*), parameter :: rout_name = 'print_output_vac'
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables 
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to print
        
        ! local variables
        type(var_1D_type), allocatable, target :: vac_1D(:)                     ! 1D equivalent of vac variables
        type(var_1D_type), pointer :: vac_1D_loc => null()                      ! local element in vac_1D
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write vacuum variables to output file')
        call lvl_ud(1)
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(vac_1D(max_dim_var_1D))
        
        ! set up variables vac_1D
        id = 1
        
        ! misc_vac
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'misc_vac'
        allocate(vac_1D_loc%tot_i_min(1),vac_1D_loc%tot_i_max(1))
        allocate(vac_1D_loc%loc_i_min(1),vac_1D_loc%loc_i_max(1))
        vac_1D_loc%loc_i_min = [1]
        vac_1D_loc%loc_i_max = [8]
        vac_1D_loc%tot_i_min = vac_1D_loc%loc_i_min
        vac_1D_loc%tot_i_max = vac_1D_loc%loc_i_max
        allocate(vac_1D_loc%p(8))
        vac_1D_loc%p = [vac%style*1._dp,vac%n_bnd*1._dp,vac%prim_X*1._dp,&
            &vac%lim_sec_X(1)*1._dp,vac%lim_sec_X(2)*1._dp,&
            &shape(vac%ang)*1._dp,vac%jq]
        
        ! ang
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'ang'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = shape(vac%ang)
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%ang)))
        vac_1D_loc%p = reshape(vac%ang,[size(vac%ang)])
        
        ! norm
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'norm'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = shape(vac%norm)
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%norm)))
        vac_1D_loc%p = reshape(vac%norm,[size(vac%norm)])
        
        ! dnorm
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'dnorm'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = shape(vac%dnorm)
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%dnorm)))
        vac_1D_loc%p = reshape(vac%dnorm,[size(vac%dnorm)])
        
        ! x_vec
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'x_vec'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = shape(vac%x_vec)
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%x_vec)))
        vac_1D_loc%p = reshape(vac%x_vec,[size(vac%x_vec)])
        
        ! RE_res
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'RE_res'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = [n_mod_X,n_mod_X]
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%res)))
        vac_1D_loc%p = reshape(rp(vac%res),[size(vac%res)])
        
        ! IM_res
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'IM_res'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = [n_mod_X,n_mod_X]
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%res)))
        vac_1D_loc%p = reshape(ip(vac%res),[size(vac%res)])
        
        ! write
        ierr = print_HDF5_arrs(vac_1D(1:id-1),PB3D_name,trim(data_name),&
            &rich_lvl=rich_lvl,ind_print=.true.)
        CHCKERR('')
        
        ! clean up
        call dealloc_var_1D(vac_1D)
        nullify(vac_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_vac
end module vac_ops
