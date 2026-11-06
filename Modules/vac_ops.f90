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
!! The final  matrix equation is  solved through Strumpack  \cite Meiser2016,
!! or through ScaLAPACK when PB3D was built without Strumpack.
!!
!! \see See \cite Weyens3D.
!!
!! Status: verified by  the full-stack test suite  in tests/fullstack (see
!! Documentation/testing.md): the axisymmetric building blocks (kernels,
!! singular integrals, assembled G and H, boundary potential solve,
!! response against the analytical cylinder limit), the free-boundary
!! chain end-to-end (regression_cbm18a_free_bnd), and the 3-D field-line
!! machinery (singular half-cell kernel, jump relations, response of an
!! axisymmetric boundary cross-checked between the two styles). The 3-D
!! style requires at least 2 field lines (n_alpha > 1, i.e. alpha_style 2,
!! which is the recommended mode for 3-D runs anyway): single-field-line
!! Weyl coverage would need per-point transverse spacings (three-distance
!! theorem) and is deliberately not implemented.
!------------------------------------------------------------------------------!
module vac_ops
#include <PB3D_macros.h>
#include <wrappers.h>
#ifdef PB3D_WITH_STRUMPACK
    use StrumpackDensePackage
#endif
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_2_type, modes_type
    use vac_vars, only: copy_vac, BLACSCTXTSIZE, &
        &vac_type

    implicit none
    private
    public store_vac, calc_vac_res, print_output_vac, vac_pot_plot, &
        &calc_GH, solve_Phi_BEM, eval_pot_ext                                   ! calc_GH, solve_Phi_BEM and eval_pot_ext are exported for the full-stack test suite
#if ldebug
    public debug_calc_GH, debug_calc_vac_res, debug_vac_pot_plot
#endif

#ifndef PB3D_WITH_STRUMPACK
    ! (comes from StrumpackDensePackage when PB3D is built with STRUMPACK)
    integer, external :: numroc                                                 !< ScaLAPACK index computation function
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_GH = .false.                                          !< plot debug information for calc_GH() \ldebug
    logical :: debug_calc_vac_res = .false.                                     !< plot debug information for calc_vac_res() \ldebug
    logical :: debug_vac_pot_plot = .false.                                     !< plot debug information for vac_pot_plot() \ldebug
#endif

contains
    !> Stores  calculation  of the  vacuum  response  by storing  the  necessary
    !! variables.
    !!
    !! This is  done by  the last process,  which is the  one that  contains the
    !! variables at the  plasma edge, and then this is  broadcasted to the other
    !! processes.
    !!
    !! The workings of this routine  are summarized in the following diagram for
    !! VMEC:
    !! @dot
    !! digraph G {
    !!  eq_job -> restart [label="yes"]
    !!  eq_job -> store [label="no"]
    !!  restart -> jump [label="yes"]
    !!  restart -> copy [label="no"]
    !!  jump -> reconstructcur [label="yes"]
    !!  jump -> restartone [label="no"]
    !!  reconstructcur -> return
    !!  restartone -> reconstruct [label="no"]
    !!  copy -> store
    !!  reconstruct -> store
    !!  restartone -> store [label="yes"]
    !!  store -> return
    !!  
    !!  eq_job [shape=diamond,label="eq_job = 1"]
    !!  store [label="store new vec, res"]
    !!  restart [shape=diamond,label="restart"]
    !!  jump [shape=diamond,label="jump"]
    !!  copy [label="copy previous level"]
    !!  reconstructcur [label="reconstruct current level"]
    !!  return [label="return"]
    !!  restartone [shape=diamond,label="restart = 1"]
    !!  reconstruct [label="reconstruct previous level"]
    !! }
    !! @enddot
    !!
    !! Before  storing the  new variables,  the procedure  checks the  following
    !! things:
    !!  - For  the first equilibrium job,  this routine copies the  results from
    !!  the previous Richardson level into  the appropriate subranges of the new
    !!  vacuum variables if no restart of Richardson levels was done.
    !!  - If a  restart was done and the  level is greater than one,  there is a
    !!  reconstruction of the  previous level's variables. But  this will happen
    !!  in calc_vac_res(), not in this procedure.
    !!
    !! If  there is  a jump  straight to  the solution,  however, the  procedure
    !! returns after reconstructing the current level's variables.
    !!
    !! If the previous \c  G and \c H variables are  empty, they are regenerated
    !! later,  in  calc_gh(). This  indicates  that  a reconstruction  happened,
    !! either  because a  Richardson restart  was  performed, or  because a  new
    !! Richardson  level was  started after  a previous  level in  which it  was
    !! jumped straight to the solution.
    !!
    !! For  HELENA, there  is only  1 equilibrium  job at  the first  Richardson
    !! level, and the vacuum has to be  calculated only once then. If there is a
    !! jump to solution or  if there is Richardson restart, it  only needs to be
    !! reconstructed.
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
        use num_vars, only: jump_to_sol
        use grid_vars, only: n_alpha
        use rich_vars, only: n_par_X
        use eq_utilities, only: calc_inv_met
        use vac_utilities, only: interlaced_vac_copy, vac_steps
        use grid_vars, only: min_par_X
        
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
        real(dp) :: dpar_X, dalpha                                              ! steps of the field-line grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: par_id(3)                                                ! total parallel index
        real(dp), allocatable :: norm_com_C(:,:,:)                          ! Cylindrical components of norm
        real(dp), allocatable :: T_CV_loc(:,:,:,:,:,:,:)                    ! local T_CV
        type(vac_type) :: vac_old                                           ! old vacuum variables
        logical :: interlaced_vac_copy_needed                               ! interlaced copy vacuum needed
        
        ! initialize ierr
        ierr = 0
        
        ! set jq
        if (use_pol_flux_F) then
            jq = eq_1%q_saf_FD(grid%loc_n_r,0)
        else
            jq = eq_1%rot_t_FD(grid%loc_n_r,0)
            ierr = 1
            err_msg = 'TOR. FLUX HAS NOT BEEN TESTED AND IS MOST PROBABLY NOT &
                &CORRECT YET'
            CHCKERR(err_msg)
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
            
            character(*), parameter :: rout_name = 'store_vac_VMEC'
            
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Start storing vacuum quantities')
            
            call lvl_ud(1)
            
            ! the vacuum of the previous Richardson level, kept or interlaced
            ierr = restore_previous_level()
            CHCKERR('')
            
            ! the points, normals, angles and metric factors of this equilibrium
            ! job, on the last process
            ierr = store_geometry()
            CHCKERR('')
            
            ! the same on the other processes
            ierr = broadcast_geometry()
            CHCKERR('')
            
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
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0

            ! user output
            call writo('Start storing vacuum quantities')
            
            call lvl_ud(1)
            
            ! if  restarting  a  Richardson  level, reconstruct.  If  not,  only
            ! calculate for the first level
            if (rich_lvl.eq.rich_restart_lvl) then
                ! check for start from zero versus restart
                if (rich_restart_lvl.eq.1) then                                 ! starting from zero
                    ! test
                    if (eq_job_nr.ne.1) then
                        ierr = 1
                        err_msg = 'HELENA should only have 1 equilibrium job &
                            &for Richardson level 1'
                        CHCKERR(err_msg)
                    end if
                    
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
                else                                                            ! restarting
                    ! reconstruct old vacuum
                    ierr = reconstruct_PB3D_vac(vac,'vac')
                    CHCKERR('')
                end if
            end if
            
            call lvl_ud(-1)
            
            call writo('Done storing vacuum quantities')
        end function store_vac_HEL
        ! Keeps or interlaces the vacuum of the previous Richardson level, or reconstructs it on a restart.
        !> \private
        integer function restore_previous_level() result(ierr)
            ierr = 0
            ! if start of Richardson level  that is not the first, copy previous
            ! results
            if (eq_job_nr.eq.1) then
                if (rich_lvl.eq.rich_restart_lvl) then                          ! restarting
                    if (jump_to_sol) then                                       ! jumping straight to solution
                        ! nothing needs to be done here
                        return
                    else                                                        ! not jumping to solution
                        if (rich_restart_lvl.eq.1) then                         ! not actually restarting, just starting
                            ! nothing extra needed
                            interlaced_vac_copy_needed = .false.
                        else                                                    ! restarting
                            ! reconstruct old vacuum
                            ierr = reconstruct_PB3D_vac(vac_old,'vac',&
                                &rich_lvl=rich_lvl-1)
                            CHCKERR('')
                        
                            ! interlaced copy necessary
                            interlaced_vac_copy_needed = .true.
                        end if
                    end if
                else                                                            ! not restarting
                    ! copy old vacuum temporarily
                    ierr = copy_vac(vac,vac_old)
                    CHCKERR('')
                
                    ! deallocate
                    call vac%dealloc()
                
                    ! interlaced copy necessary
                    interlaced_vac_copy_needed = .true.
                end if
            
                ! allocate
                ierr = vac%init(1,n_par_X*n_alpha,prim_X,[n_par_X,n_alpha],jq)
                CHCKERR('')
            
                ! interlaced copy if needed
                if (interlaced_vac_copy_needed) then
                    ierr = interlaced_vac_copy(vac_old,vac)
                    CHCKERR('')
                end if

                ! the parallel angle of every point, as calc_vac_res sets it
                ! up for the normal data (written to the output through the
                ! limits in misc_vac, from which POST rebuilds it)
                call vac_steps(dpar_X,dalpha)
                do id = 1,n_par_X
                    vac%ang(id,:) = pi*min_par_X + (id-1)*dpar_X
                end do
            end if

        end function restore_previous_level
        
        ! Stores the Cartesian points, the normals -J nabla psi, the parallel angles and the metric factors of the current equilibrium job (last process).
        !> \private
        integer function store_geometry() result(ierr)
            ierr = 0
            ! add results from current equilibrium job to the variables
            if (rank.eq.n_procs-1) then
                ! grid point to save is last local
                r_id = grid%loc_n_r
            
                ! calculate Cartesian components of - J nabla psi
                !   = - J_F (Dpsi_pol/Dr_V) / 2pi (T_C^V)^T (e^C)
                !   = - J_F T_F^V(1,2) (T_C^V)^T (e^C)
                ! The  minus sign  is  the convention  of  the kernels:  calc_GH
                ! takes  the inward  normal, as  store_vac_HEL provides  it for
                ! the  axisymmetric style  and as  the analytic  torus of  the
                ! full-stack test  does (its normal  is e_theta x  e_alpha, which
                ! points inward  for a  left-handed (alpha, psi,  theta) frame).
                ! With +J nabla psi the off-diagonal H has the wrong sign, the
                ! row-sum identity then sets a  wrong diagonal, and the response
                ! of  the circular  VMEC twin  came out  3 to  4 times  too small
                ! against the  HELENA one;  with -J nabla  psi the  two agree
                ! within  the first-order  discretization error,  as  in  the
                ! full-stack test response_1_vs_2 (issue #4, step 6e).
                allocate(norm_com_C(grid%n(1),grid%n(2),3))
                norm_com_C = 0._dp
                allocate(T_CV_loc(size(eq_2%T_VC,1),size(eq_2%T_VC,2),1:1,&
                    &size(eq_2%T_VC,4),0:0,0:0,0:0))
                ierr = calc_inv_met(T_CV_loc,eq_2%T_VC(:,:,r_id:r_id,:,:,:,:),&
                    &[0,0,0])
                CHCKERR('')
                do id = 1,3
                    norm_com_C(:,:,id) = - eq_2%jac_FD(:,:,r_id,0,0,0) * &
                        &eq_2%T_EF(:,:,r_id,c([1,2],.false.),0,0,0)*&
                        &T_CV_loc(:,:,1,c([id,1],.false.),0,0,0)                ! Cylindrical contravariant components
                end do
                deallocate(T_CV_loc)
            
                ! iterate over all field lines
                do jd = 1,n_alpha
                    ! set total parallel id
                    if (rich_lvl.eq.1) then
                        par_id = [eq_jobs_lims(:,eq_job_nr),1]
                    else
                        par_id = [2*eq_jobs_lims(:,eq_job_nr),2]
                    end if
                    par_id(1:2) = par_id(1:2) + (jd-1)*n_par_X
                
                    ! save X, Y and Z
                    vac%x_vec(par_id(1):par_id(2):par_id(3),1) = &
                        &eq_2%R_E(:,jd,r_id,0,0,0) * cos(grid%zeta_E(:,jd,r_id))
                    vac%x_vec(par_id(1):par_id(2):par_id(3),2) = &
                        &eq_2%R_E(:,jd,r_id,0,0,0) * sin(grid%zeta_E(:,jd,r_id))
                    vac%x_vec(par_id(1):par_id(2):par_id(3),3) = &
                        &eq_2%Z_E(:,jd,r_id,0,0,0)
                
                    ! save norm
                    vac%norm(par_id(1):par_id(2):par_id(3),1) = &
                        &norm_com_C(:,jd,1)*cos(grid%zeta_E(:,jd,r_id)) - &
                        &norm_com_C(:,jd,2)*sin(grid%zeta_E(:,jd,r_id)) / &
                        &eq_2%R_E(:,jd,r_id,0,0,0)                              ! ~ e^X
                    vac%norm(par_id(1):par_id(2):par_id(3),2) = &
                        &norm_com_C(:,jd,1)*sin(grid%zeta_E(:,jd,r_id)) + &
                        &norm_com_C(:,jd,2)*cos(grid%zeta_E(:,jd,r_id)) / &
                        &eq_2%R_E(:,jd,r_id,0,0,0)                              ! ~ e^X
                    vac%norm(par_id(1):par_id(2):par_id(3),3) = &
                        &norm_com_C(:,jd,3)                                     ! ~ e^Z

                    ! save metric factors
                    vac%h_fac(par_id(1):par_id(2):par_id(3),1) = &
                        &eq_2%g_FD(:,jd,r_id,c([1,1],.true.),0,0,0)             ! g_11
                    vac%h_fac(par_id(1):par_id(2):par_id(3),2) = &
                        &eq_2%g_FD(:,jd,r_id,c([1,3],.true.),0,0,0)             ! g_13
                    vac%h_fac(par_id(1):par_id(2):par_id(3),3) = &
                        &eq_2%g_FD(:,jd,r_id,c([3,3],.true.),0,0,0)             ! g_33
                end do
            end if
        
        end function store_geometry
        
        ! Broadcasts the stored geometry to the other processes.
        !> \private
        integer function broadcast_geometry() result(ierr)
            ierr = 0
            ! broadcast to other processes
            do id = 1,size(vac%x_vec,2)
                ierr = broadcast_var(vac%x_vec(:,id),n_procs-1)
                CHCKERR('')
                ierr = broadcast_var(vac%norm(:,id),n_procs-1)
                CHCKERR('')
            end do
            do id = 1,size(vac%h_fac,2)
                ierr = broadcast_var(vac%h_fac(:,id),n_procs-1)
                CHCKERR('')
            end do
        
        end function broadcast_geometry
        
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
    !! Optionally, this  procedure can be  used to calculate the  coefficients G
    !! and H  with respect  to other  points, outside of  the boundary.  This is
    !! useful when the potential is to be calculated at external points.
    !!
    !! In this case, however, no (near-)  singular points are calculated, and if
    !! they appear any way, perhaps by accident, they will not be accurate.
    !!
    !! \Note Multiple  field lines are  stored sequentially, which  implies that
    !! the integral between the  last point on a field line  and the first point
    !! on the next field lines needs to be left out.
    !!
    !! \note For style 1 the kernel  is integrated along the field line within
    !! four interval  lengths of the  influence point (vac_utilities.calc_gh_int_1),
    !! because  neighbouring field  lines can  lie much  closer together  than
    !! the points along them.
    !!
    !! \todo For  3-D vacuums,  step sizes  are constant.  Subsequent Richardson
    !! levels should therefore make use of the fact that the contribution to the
    !! points inherited from  the previous levels can be just  scaled by 1/2 and
    !! do not need to be recalculated.  Currently, the copy is done correctly in
    !! interlaced_vac_copy(), but they are afterwards overwritten.
    !!
    !! \return ierr
    integer function calc_GH(vac,n_r_in,lims_r_in,x_vec_in,G_in,H_in) &
        &result(ierr)
        
        use vac_vars, only: in_context
        use vac_utilities, only: vec_dis2loc
#if ldebug
        use num_vars, only: rank, n_procs
        use vac_utilities, only: mat_dis2loc
#endif

        character(*), parameter :: rout_name = 'calc_GH'
        
        ! input / output
        type(vac_type), intent(inout), target :: vac                            !< vacuum variables
        integer, intent(in), optional :: n_r_in                                 !< total number of rows of external poins of influence
        integer, intent(in), optional, target :: lims_r_in(:,:)                 !< row limits of external points of influence
        real(dp), intent(in), optional, target :: x_vec_in(:,:)                 !< external position of influence 
        real(dp), intent(in), optional, target :: G_in(:,:)                     !< external G
        real(dp), intent(in), optional, target :: H_in(:,:)                     !< external H
        
        ! local variables
        integer :: n_r                                                          ! number of rows
        integer, pointer :: lims_r(:,:)                                         ! row limits
        real(dp), pointer :: G(:,:)                                             ! G
        real(dp), pointer :: H(:,:)                                             ! H
        real(dp), pointer :: x_vec(:,:)                                         ! x_vec used in row
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: ext_in                                                       ! external influence point
#if ldebug
        real(dp), allocatable :: loc_ser(:,:)                                   ! serial copy of a matrix
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate the G and H matrices')
        call lvl_ud(1)
        
        ! test
        if (present(n_r_in) .and. .not.present(x_vec_in) .or. &
            &present(n_r_in) .and. .not.present(lims_r_in) .or. &
            &present(n_r_in) .and. .not.present(G_in) .or. &
            &present(n_r_in) .and. .not.present(H_in)) then
            ierr = 1
            CHCKERR('Need all optional variables')
        end if
        
        ! set limits and variables
        if (present(lims_r_in)) then
            n_r = n_r_in
            lims_r => lims_r_in
            x_vec => x_vec_in
            G => G_in
            H => H_in
            ext_in = .true.
        else
            n_r = vac%n_bnd
            allocate(lims_r(2,size(vac%lims_r,2)))
            allocate(x_vec(vac%n_bnd,size(vac%x_vec,2)))
            lims_r = vac%lims_r
            x_vec = vac%x_vec
            G => vac%G
            H => vac%H
            ext_in = .false.
            ! both styles accumulate into G and H; after a Richardson step
            ! the field-line arrays still hold the interlaced copy of the
            ! previous level (store_vac_VMEC), which would be counted twice
            G = 0._dp
            H = 0._dp
        end if
        
        ! call specific procedure for vacuum style
        select case (vac%style)
            case (1)                                                            ! field-line 3-D
                ierr = calc_GH_1(vac,n_r,lims_r,x_vec,G,H,ext_in)
                CHCKERR('')
            case (2)                                                            ! axisymmetric
                ierr = calc_GH_2(vac,n_r,lims_r,x_vec,G,H,ext_in)
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No vacuum style '//trim(i2str(vac%style))//&
                    &' possible'
                CHCKERR(err_msg)
        end select
        
        ! clean up
        if (.not.present(lims_r_in)) then
            deallocate(lims_r)
            deallocate(x_vec)
        end if
        nullify(lims_r)
        nullify(x_vec)
        nullify(G)
        nullify(H)
        
#if ldebug
        if (debug_calc_GH .and. in_context(vac%ctxt_HG)) then
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
            
        end if
#endif
        
        call lvl_ud(-1)
    end function calc_GH

    !> \public Calculates matrices \c G and \c H in 3-D configuration.
    !!
    !! \see calc_gh.
    integer function calc_GH_1(vac,n_r_in,lims_r_in,x_vec_in,G_in,H_in,ext_in) &
        &result(ierr)
        
        use vac_utilities, only: calc_GH_int_1, dis_matvec, vac_steps, &
            &calc_diag_direct
        use grid_vars, only: n_alpha
        use num_vars, only: tol_zero, vac_diag_style
        use vac_vars, only: in_context, loc_idx
        
        character(*), parameter :: rout_name = 'calc_GH_1'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        integer, intent(in) :: n_r_in                                           !< total number of rows of external poins of influence
        integer, intent(in) :: lims_r_in(:,:)                                   !< row limits of external points of influence
        real(dp), intent(in) :: x_vec_in(:,:)                                   !< position of influence 
        real(dp), intent(in), pointer :: G_in(:,:)                              !< G at position of influence
        real(dp), intent(in), pointer :: H_in(:,:)                              !< H at position of influence
        logical, intent(in) :: ext_in                                           !< position of influence is external
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: i_rd, i_cd                                                   ! index of subrow and column
        integer :: rd, cd                                                       ! global counters for row and column
        integer :: rdl, cdl                                                     ! local counters for row and column, used for G and H
        integer :: lims_c_loc(2)                                                ! local column limits
        integer :: st                                                           ! first point of the four-point stencil along the field line
        real(dp) :: dpar_X                                                      ! step size in par_X 
        real(dp) :: dalpha                                                      ! step size in alpha
        real(dp) :: tol_loc                                                     ! local tolerance on rho for singular points
        real(dp) :: G_loc(2)                                                    ! local G
        real(dp) :: H_loc(2)                                                    ! local H
        real(dp) :: h_fac_in(3)                                                 ! metric factors for the singular half-cell
        integer :: j_min, j_nb                                                  ! boundary point nearest to an external point and its neighbour on the line
        logical :: near_row                                                     ! near-field rule allowed for this row
        real(dp), allocatable :: loc_res(:)                                     ! global vector: diagonal values or data
        real(dp), allocatable :: lin_res(:)                                     ! global vector: matrix-vector products
        real(dp), allocatable :: diag_dir(:), diag_fac(:), I_sing(:), omega_v(:) ! direct diagonals of H and their parts
        logical, allocatable :: avail(:)                                        ! direct value available
        character(len=max_str_ln) :: err_msg                                    ! error message

        ! initialize ierr
        ierr = 0
        
        call writo('Using field-line 3-D Boundary Element Method')

        ! initialize variables
        call vac_steps(dpar_X,dalpha)
        if (n_alpha.le.1) then
            ! a single field line makes the transverse cell size (and with it
            ! the singular self-integrals and the integration rule weights of
            ! calc_vac_res) degenerate; covering the surface with one long
            ! field line through Weyl's theorem needs an effective transverse
            ! spacing that is not implemented
            ierr = 1
            err_msg = '3-D vacuum needs at least 2 field lines (n_alpha > 1)'
            CHCKERR(err_msg)
        end if
        if (vac%n_ang(1).lt.4) then
            ierr = 1
            err_msg = '3-D vacuum needs at least 4 points per field line'
            CHCKERR(err_msg)
        end if
        
        ! set up local tolerance
        tol_loc = tol_zero*10                                                   ! tolerance for singular interval
        
        ! iterate over all subrows: position of singular point
        subrows: do i_rd = 1,size(lims_r_in,2)
            ! iterate over all rows of this subrow
            row: do rd = lims_r_in(1,i_rd),lims_r_in(2,i_rd)
                ! set local row index
                rdl = loc_idx(lims_r_in,i_rd,rd)
                
                ! the near-field rule along the lines: always for a boundary
                ! point; for an external point only if it lies within an
                ! interval length of the nearest boundary point (see
                ! calc_GH_int_1; every process holds all boundary points)
                near_row = .true.
                if (ext_in) then
                    j_min = minloc(sum((vac%x_vec-spread(x_vec_in(rd,:),1,&
                        &vac%n_bnd))**2,2),1)
                    j_nb = j_min+1                                              ! neighbour on the same line
                    if (mod(j_min,vac%n_ang(1)).eq.0) j_nb = j_min-1
                    near_row = sum((vac%x_vec(j_min,:)-x_vec_in(rd,:))**2).lt.&
                        &sum((vac%x_vec(j_nb,:)-vac%x_vec(j_min,:))**2)
                end if
                
                ! iterate over all subcolumns: subintegrals
                subcols: do i_cd = 1,size(vac%lims_c,2)
                    ! set up local limits, including possible overlap
                    lims_c_loc(1) = max(1,vac%lims_c(1,i_cd)-1)
                    lims_c_loc(2) = min(vac%n_bnd,vac%lims_c(2,i_cd)+1)
                    
                    ! iterate over all pairs of columns of this subcolumn
                    ! (upper limit is one lower as lims_c_loc
                    col: do cd = lims_c_loc(1),lims_c_loc(2)-1
                        ! skip the interval that bridges two field lines: its
                        ! geometry connects unrelated points
                        if (mod(cd,vac%n_ang(1)).eq.0) cycle col
                        
                        ! set local column index
                        cdl = loc_idx(vac%lims_c,i_cd,cd)
                        
                        ! four-point stencil along the field line, one-sided
                        ! at its ends
                        st = (cd-1)/vac%n_ang(1)*vac%n_ang(1)                   ! last point of the previous field line
                        st = max(st+1,min(cd-1,st+vac%n_ang(1)-3))

                        ! metric factors at the influence point, used only
                        ! for the singular half-cell; an external point has
                        ! none (its row index does not number a boundary
                        ! point), and the source point's are used should it
                        ! coincide with one
                        if (ext_in) then
                            h_fac_in = vac%h_fac(cd,:)
                        else
                            h_fac_in = vac%h_fac(rd,:)
                        end if

                        ! calculate contributions to H and G
                        call calc_GH_int_1(G_loc,H_loc,&
                            &vac%x_vec(st:st+3,:),x_vec_in(rd,:),&
                            &vac%norm(st:st+3,:),h_fac_in,&
                            &[dpar_X,dalpha],tol_loc**2,cd-st+1,&
                            &use_near=near_row)

                        ! loop over left and right side of interval
                        ! (the contributions are weighted with the half-cell
                        ! surface element dpar dalpha / 2, so that G and H are
                        ! true discretized surface operators, consistently
                        ! with the absolute -4 pi row-sum identity that sets
                        ! the H diagonal below, and with the axisymmetric
                        ! style, which carries its weights internally as well)
                        do kd = 0,1
                            if (cd+kd.ge.vac%lims_c(1,i_cd) .and. &             ! lower local bound
                                &cd+kd.le.vac%lims_c(2,i_cd)) then              ! upper local bound
                                G_in(rdl,cdl+kd) = G_in(rdl,cdl+kd) + &
                                    &G_loc(1+kd)*dpar_X*dalpha/2._dp
                                H_in(rdl,cdl+kd) = H_in(rdl,cdl+kd) + &
                                    &H_loc(1+kd)*dpar_X*dalpha/2._dp
                            end if
                        end do
                    end do col
                end do subcols
            end do row
        end do subrows
        
        ! if in context and not external, set the diagonals.
        !  - indirect (vac_diag_style 1, the default): H through
        !    the row-sum identity for the constant potential and G through the
        !    identities for the linear interior harmonics, which absorb the
        !    first-order quadrature error of the off-diagonal entries into the
        !    diagonal;
        !  - direct (vac_diag_style 2): H from the corner term of the smooth
        !    surface, -2 pi, plus the cell integral of the kernel
        !    (calc_diag_direct; the row sum at the end points of the lines,
        !    which have no full cell), G from the analytical cell integral of
        !    calc_GH_int_1, already in place.
        if (in_context(vac%ctxt_HG) .and. .not.ext_in) then
            allocate(loc_res(vac%n_bnd),lin_res(vac%n_bnd))
            
            ! H: the row-sum identity H 1 = -4 pi 1 (interior harmonic 1)
            loc_res = 0._dp
            call set_diag(H_in,loc_res,.false.)                                 ! diagonal of H to zero
            loc_res = 1._dp
            ierr = dis_matvec(vac,vac%H,vac%desc_H,loc_res,lin_res)             ! sum of rows of H
            CHCKERR('')
            loc_res = -lin_res-4*pi
            if (vac_diag_style.eq.2) then                                       ! direct where available
                allocate(diag_dir(vac%n_bnd),diag_fac(vac%n_bnd),&
                    &I_sing(vac%n_bnd),omega_v(vac%n_bnd),avail(vac%n_bnd))
                ierr = calc_diag_direct(vac,dpar_X,dalpha,diag_dir,diag_fac,&
                    &I_sing,omega_v,avail)
                CHCKERR('')
                where (avail) loc_res = diag_dir
                deallocate(diag_dir,diag_fac,I_sing,omega_v,avail)
            end if
            call set_diag(H_in,loc_res,.false.)
            
            ! G: the identity H x - G norm_x = -4 pi x for the linear interior
            ! harmonics x, y, z, projected on the normal: the analytical
            ! half-cell integral of calc_GH_int_1 is exact for the singular
            ! cell, but the one-point rule across the neighbouring field lines
            ! leaves a first-order error that depends on the lattice of
            ! points; the identity absorbs it into the diagonal, as the row
            ! sum does for H (issue #4, 31 October 2026)
            if (vac_diag_style.ne.2) then
                lin_res = 0._dp
                do kd = 1,3
                    ierr = dis_matvec(vac,vac%H,vac%desc_H,vac%x_vec(:,kd),&
                        &loc_res)
                    CHCKERR('')
                    lin_res = lin_res + vac%norm(:,kd)*&
                        &(loc_res+4*pi*vac%x_vec(:,kd))
                    ierr = dis_matvec(vac,vac%G,vac%desc_G,vac%norm(:,kd),&
                        &loc_res)
                    CHCKERR('')
                    lin_res = lin_res - vac%norm(:,kd)*loc_res
                end do
                lin_res = lin_res/sum(vac%norm**2,2)
                call set_diag(G_in,lin_res,.true.)
            end if
            
            ! clean up
            deallocate(loc_res,lin_res)
        end if
    contains
        ! Sets the diagonal of the local matrix A_in to vec, or adds vec to it.
        subroutine set_diag(A_in,vec,add)
        use vac_vars, only: loc_idx
            ! input / output
            real(dp), intent(in), pointer :: A_in(:,:)                          ! local matrix
            real(dp), intent(in) :: vec(:)                                      ! global vector of diagonal values
            logical, intent(in) :: add                                          ! add instead of set
            
            ! local variables
            integer :: i_rd, i_cd                                               ! index of subrow and column
            integer :: rd                                                       ! global row
            integer :: rdl, cdl                                                 ! local row and column
            
            do i_rd = 1,size(lims_r_in,2)
                do rd = lims_r_in(1,i_rd),lims_r_in(2,i_rd)
                    rdl = loc_idx(lims_r_in,i_rd,rd)
                    do i_cd = 1,size(vac%lims_c,2)
                        if (vac%lims_c(1,i_cd).le.rd .and. &
                            &rd.le.vac%lims_c(2,i_cd)) then                     ! this subcolumn includes the diagonal
                            cdl = loc_idx(vac%lims_c,i_cd,rd)
                            if (add) then
                                A_in(rdl,cdl) = A_in(rdl,cdl) + vec(rd)
                            else
                                A_in(rdl,cdl) = vec(rd)
                            end if
                        end if
                    end do
                end do
            end do
        end subroutine set_diag
    end function calc_GH_1
    
    !> \public Calculates matrices \c G and \c H in axisymmetric configurations.
    !!
    !! It makes use of toroidal functions.
    !!
    !! \see calc_gh.
    integer function calc_GH_2(vac,n_r_in,lims_r_in,x_vec_in,G_in,H_in,ext_in) &
        &result(ierr)
        
        use dtorh, only: dtorh1
        use vac_utilities, only: calc_GH_int_2
        use num_ops, only: calc_zero_HH, calc_zero_Zhang
        use num_vars, only: vac_diag_style
        use vac_vars, only: in_context, loc_idx
#if ldebug
        use num_vars, only: rank
        use MPI_utilities, only: get_ser_var
#endif

        character(*), parameter :: rout_name = 'calc_GH_2'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        integer, intent(in) :: n_r_in                                           !< total number of rows of external poins of influence
        integer, intent(in) :: lims_r_in(:,:)                                   !< row limits of external points of influence
        real(dp), intent(in) :: x_vec_in(:,:)                                   !< position of influence 
        real(dp), intent(in), pointer :: G_in(:,:)                              !< G at position of influence
        real(dp), intent(in), pointer :: H_in(:,:)                              !< H at position of influence
        logical, intent(in) :: ext_in                                           !< position of influence is external
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: i_rd, i_cd                                                   ! index of subrow and column
        integer :: rd, cd                                                       ! global counters for row and column
        integer :: rdl, cdl                                                     ! local counters for row and column, used for G and H
        integer :: max_n                                                        ! maximum degree
        integer :: n_abs                                                        ! absolute value of toroidal mode number
        integer :: lims_c_loc(2)                                                ! local column limits
        real(dp) :: b_coeff(2)                                                  ! sum_k=1^n 2/2k-1 for n and n-1
        real(dp) :: tol_loc                                                     ! local tolerance on rho for singular points
        real(dp) :: del_r(2,2)                                                  ! unit vectors of next and previous subinterval
        real(dp) :: G_loc(2)                                                    ! local G
        real(dp) :: H_loc(2)                                                    ! local H
        real(dp) :: beta_comp                                                   ! complementary beta
        real(dp) :: ang                                                         ! local angle
        real(dp), parameter :: tol_Q = 1.E-6_dp                                 ! tolerance on Q approximation
        real(dp), allocatable :: beta(:)                                        ! beta
        real(dp), allocatable :: rho2(:)                                        ! square of distance in projected poloidal plane for one column
        real(dp), allocatable :: arg(:)                                         ! argument 1 + rho^2/(RiRj)(
        real(dp), allocatable :: Aij(:)                                         ! Aij helper variable for H
        real(dp), allocatable :: ql(:,:,:)                                      ! ql for all rho2's
        real(dp), allocatable :: pl_loc(:), ql_loc(:)                           ! local toroidal functions
        logical, allocatable :: was_calc(:,:)                                   ! was already calculated
        integer :: j_sing                                                       ! boundary point closest to an external point
        real(dp) :: norm_in(2)                                                  ! normal at the influence point
        logical :: was_calc_sym                                                 ! symmetrical was_calc
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: dH(:), dG(:)                                   ! indirect diagonals
#if ldebug
        real(dp), allocatable :: loc_ser(:)                                     ! dummy serial variable
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

        ! set up b_coeff and local tolerance
        ! (the kernels only depend on |n|: the toroidal decomposition of the
        ! Green's function contains Q_{|m|-1/2}, equivalently
        ! Q_{-nu-1} = Q_{nu} for the half-integer degrees involved)
        n_abs = abs(vac%prim_X)
        b_coeff = 0._dp
        do kd = 1,n_abs
            b_coeff(2) = b_coeff(2) + 2._dp/(2*kd-1)
        end do
        b_coeff(1) = b_coeff(2) - 2._dp/(2*n_abs-1)
        err_msg = calc_zero_Zhang(tol_loc,fun_tol,&
            &[1.E-10_dp,&                                                       ! very low lower limit
            &exp(-2*b_coeff(2))])                                               ! upper limit for -1/2ln(x) = b
        if (trim(err_msg).ne.'') then
            ierr = 1
            CHCKERR(trim(err_msg))
        end if
        tol_loc = 8*minval(vac%x_vec(:,1))*tol_loc                              ! tolerance on rho
        ! (the same on every process, which all hold every boundary point:
        ! no gather, so that the processes of the BLACS grid alone can call
        ! this routine with external points, as check_vac_exact_2 does)
        call writo('Tolerance on rho for analytical integral: '//&
            &trim(r2str(tol_loc)))
        
        ! allocate helper variables
        allocate(was_calc(1:n_r_in,1:vac%n_bnd))
        was_calc = .false.                                                      ! will remain false if external influence
        allocate(ql(1:n_r_in,1:vac%n_bnd,2))
        allocate(pl_loc(0:max(n_abs,1)))
        allocate(ql_loc(0:max(n_abs,1)))
#if ldebug
        rho2_lims = [huge(1._dp),-huge(1._dp)]
        gamma_lims = [huge(1._dp),-huge(1._dp)]
#endif
        
        ! iterate over all subrows: position of singular point
        subrows: do i_rd = 1,size(lims_r_in,2)
            ! initialize helper variables
            allocate(beta(lims_r_in(1,i_rd):lims_r_in(2,i_rd)))                 ! global index
#if ldebug
            beta_lims = [lims_r_in(2,i_rd)+1,lims_r_in(1,i_rd)-1]
#endif
            
            ! iterate over all rows of this subrow
            row: do rd = lims_r_in(1,i_rd),lims_r_in(2,i_rd)
                ! set local row index
                rdl = loc_idx(lims_r_in,i_rd,rd)

                ! An external point within the tolerance of a boundary point
                ! (checked against all of them, which every process holds) is
                ! that point: its row gives the boundary value, 4 pi phi_j
                ! with 4 pi in the column of the point and nothing else. The
                ! other external points are farther than the tolerance from
                ! every boundary point, so calc_GH_int_2 never takes its
                ! analytical branch, which needs the angle and normal of an
                ! influence point on the boundary.
                if (ext_in) then
                    j_sing = minloc(sum((vac%x_vec-spread(x_vec_in(rd,:),1,&
                        &vac%n_bnd))**2,2),1)
                    if (sum((vac%x_vec(j_sing,:)-x_vec_in(rd,:))**2).le.&
                        &tol_loc**2) then
                        do i_cd = 1,size(vac%lims_c,2)
                            do cd = vac%lims_c(1,i_cd),vac%lims_c(2,i_cd)
                                cdl = loc_idx(vac%lims_c,i_cd,cd)
                                G_in(rdl,cdl) = 0._dp
                                H_in(rdl,cdl) = 0._dp
                                if (cd.eq.j_sing) H_in(rdl,cdl) = 4*pi
                            end do
                        end do
                        cycle row
                    end if
                end if

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
                        rho2(cd) = sum((vac%x_vec(cd,:)-x_vec_in(rd,:))**2)
                        arg(cd) = 1._dp + &
                            &rho2(cd)/(2*vac%x_vec(cd,1)*x_vec_in(rd,1))
                        was_calc_sym = .false.
                        if (.not.ext_in) was_calc_sym = was_calc(cd,rd)
                        if (rho2(cd).le.tol_loc**2) then
                            ! don't set it
                        else if (was_calc_sym) then                             ! check symmetric element
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
                            if (n_abs.gt.0) then                                ! nonzero n: kernels for |n|
                                ierr = dtorh1(arg(cd),0,n_abs,pl_loc,&
                                    &ql_loc,max_n)
                                CHCKERR('')
                                if (max_n.lt.n_abs) then
                                    ierr = 1
                                    CHCKERR(err_msg)
                                end if
                                ql(rd,cd,:) = ql_loc(n_abs-1:n_abs)
                            else                                                ! n = 0
                                ierr = dtorh1(arg(cd),0,1,pl_loc,ql_loc,max_n)
                                CHCKERR('')
                                if (max_n.lt.n_abs) then
                                    ierr = 1
                                    CHCKERR(err_msg)
                                end if
                                ql(rd,cd,1) = ql_loc(1)                         ! Q_-3/2 = Q_1/2
                                ql(rd,cd,2) = ql_loc(0)                         ! Q_-1/2
                            end if
                            if (.not.ext_in) was_calc(rd,cd) = .true.
                        end if
                        
                        ! set up Aij (an external row never gets here with a
                        ! singular column, see above)
                        if (rho2(cd).le.tol_loc**2) then
                            Aij(cd) = 0.5_dp*(n_abs - 0.5_dp) * &
                                &(x_vec_in(rd,1) * &
                                &(vac%norm(rd,1)*vac%dnorm(rd,2)-&
                                &vac%norm(rd,2)*vac%dnorm(rd,1))/&
                                &sum(vac%norm(rd,:)**2) &
                                &+vac%norm(rd,1)/x_vec_in(rd,1))
                        else
                            Aij(cd) = 2._dp*x_vec_in(rd,1)*vac%x_vec(cd,1) / &
                                &(4._dp*x_vec_in(rd,1)*vac%x_vec(cd,1)+&
                                &rho2(cd)) *(n_abs - 0.5_dp) * &
                                &(- 2._dp/rho2(cd) * sum(vac%norm(cd,:)*&
                                &(vac%x_vec(cd,:)-x_vec_in(rd,:))) + &
                                &vac%norm(cd,1)/vac%x_vec(cd,1))
                        end if
                    end do col
                    
#if ldebug
                    if (debug_calc_GH) then
                    end if
#endif
                    
                    ! iterate over all pairs of columns of this subcolumn
                    ! (upper limit is one lower as previous)
                    col2: do cd = lims_c_loc(1),lims_c_loc(2)-1
                        ! set local column index
                        cdl = loc_idx(vac%lims_c,i_cd,cd)

                        ! angle and normal of the influence point, used only
                        ! in the analytical branch of calc_GH_int_2, which an
                        ! external point never takes (its row index does not
                        ! number a boundary point)
                        if (ext_in) then
                            ang = 0._dp
                            norm_in = vac%norm(cd,:)
                        else
                            ang = vac%ang(rd,1)
                            norm_in = vac%norm(rd,:)
                        end if

                        ! calculate contributions to H and G
                        call calc_GH_int_2(G_loc,H_loc,&
                            &vac%ang(cd:cd+1,1),ang,&
                            &vac%x_vec(cd:cd+1,:),x_vec_in(rd,:),&
                            &vac%norm(cd:cd+1,:),norm_in,&
                            &Aij(cd:cd+1),ql(rd,cd:cd+1,:),&
                            &tol_loc**2,b_coeff,n_abs)
                        
                        do kd = 0,1
                            if (cd+kd.ge.vac%lims_c(1,i_cd) .and. &
                                &cd+kd.le.vac%lims_c(2,i_cd)) then
                                G_in(rdl,cdl+kd) = G_in(rdl,cdl+kd) + &
                                    &G_loc(1+kd)
                                H_in(rdl,cdl+kd) = H_in(rdl,cdl+kd) + &
                                    &H_loc(1+kd)
                            end if
                        end do
                    end do col2
                    
#if ldebug
                    if (debug_calc_GH) then
                        !! local limits
                    end if
#endif
                    
                    ! calculate contribution due to fundamental solution
                    ! (only if not external)
                    if (.not.ext_in .and. vac%lims_c(1,i_cd).le.rd .and. &
                        &rd.le.vac%lims_c(2,i_cd)) then                         ! this subcolumn includes the diagonal
                        ! set local column index
                        cdl = loc_idx(vac%lims_c,i_cd,rd)
                        
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
                        
                        ! corner term of the left-hand side of the Green
                        ! identity, -2 beta with beta the vacuum's internal
                        ! angle (> pi where the plasma is convex, eq. 146 of
                        ! the 2019 manuscript): with the kernels taken along
                        ! the stored normal -J nabla psi (calc_GH_int_2), the
                        ! assembled H is the vacuum-side operator, as for the
                        ! field-line style, and solve_Phi_BEM needs no shift
                        H_in(rdl,cdl) = H_in(rdl,cdl) - 2*beta(rd)
                    end if
                    
                    ! clean up
                    deallocate(arg,rho2,Aij)
                end do subcols
            end do row
            
#if ldebug
            if (debug_calc_GH) then
            end if
#endif
            
            ! clean up
            deallocate(beta)
        end do subrows
        
        ! indirect diagonals (vac_diag_style 1): H_ii and G_ii from the
        ! identities for the test potentials of set_pot_2 (diag_indirect_2),
        ! the counterpart of the identities of the field-line style; the
        ! direct diagonals (the corner angle and the analytical near-singular
        ! integrals, vac_diag_style 2) are first order in the poloidal points,
        ! the indirect ones third order (full-stack test exact_exterior_2d_options)
        if (in_context(vac%ctxt_HG) .and. .not.ext_in .and. &
            &vac_diag_style.eq.1) then
            allocate(dH(vac%n_bnd),dG(vac%n_bnd))
            ierr = diag_indirect_2(vac,dH,dG)
            CHCKERR('')
            call set_diag_vac(vac,H_in,dH)
            call set_diag_vac(vac,G_in,dG)
            deallocate(dH,dG)
        end if
        
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
    
    !> Test potentials of the axisymmetric vacuum, varying as
    !! \f$e^{\imath n\zeta}\f$, with their normal data along the stored normal
    !! \f$-\mathcal{J}\nabla\psi\f$ and the free term \c c_free of the identity
    !! \f$\overline{\text{H}}\phi - \overline{\text{G}}\,\text{d}\phi = c\,\phi\f$
    !! that they satisfy:
    !!  1. \f$R^n\f$, harmonic inside the boundary, \f$c = -4\pi\f$ (the shell
    !!     at infinity contributes \f$-4\pi\phi(\mathbf{r}_i)\f$, as for the
    !!     constant potential; eq. 149 of the 2019 manuscript lacks this term,
    !!     erratum E16);
    !!  2. \f$R^n Z\f$, harmonic inside the boundary, \f$c = -4\pi\f$;
    !!  3. \f$Q_{|n|-1/2}(\gamma(\mathbf{x},\mathbf{x}_\text{src}))/\sqrt{R
    !!     R_\text{src}}\f$, the n-mode Green's function of a ring source at
    !!     \c x_src inside the plasma, harmonic outside the boundary and
    !!     decaying, \f$c = 0\f$; its normal data by a central difference.
    !! (The spherical potential of eq. 152 of the manuscript, which the debug
    !! self-test of calc_GH used until November 2026, is not harmonic across
    !! \f$Z = 0\f$, erratum E17.)
    integer function set_pot_2(vac,kd,x_src,phi,dphi,c_free,x_ext,phi_ext) &
        &result(ierr)
        use dtorh, only: dtorh1
        
        character(*), parameter :: rout_name = 'set_pot_2'
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        integer, intent(in) :: kd                                               !< which potential
        real(dp), intent(in) :: x_src(2)                                        !< ring source (R, Z)
        real(dp), intent(out) :: phi(:)                                         !< potential at the boundary points
        real(dp), intent(out) :: dphi(:)                                        !< its normal data
        real(dp), intent(out) :: c_free                                         !< free term of the identity
        real(dp), intent(in), optional :: x_ext(:,:)                            !< points off the boundary, (R, Z)
        real(dp), intent(out), optional :: phi_ext(:)                           !< potential at \c x_ext
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_abs                                                        ! |n|
        integer :: max_n                                                        ! maximum degree reached by dtorh1
        real(dp) :: fp, fm                                                      ! potential values for the central difference
        real(dp) :: x_pert(2)                                                   ! point displaced along the normal
        real(dp), allocatable :: pl_loc(:), ql_loc(:)                           ! toroidal functions
        real(dp), parameter :: eps_src = 1.E-7_dp                               ! step of the central difference
        
        ierr = 0
        n_abs = abs(vac%prim_X)
        c_free = -4._dp*pi
        select case (kd)
            case (1)                                                            ! R^n
                phi = vac%x_vec(:,1)**n_abs
                dphi = n_abs*vac%norm(:,1)*vac%x_vec(:,1)**(n_abs-1)
                if (present(x_ext)) phi_ext = x_ext(:,1)**n_abs
            case (2)                                                            ! R^n Z
                phi = vac%x_vec(:,1)**n_abs*vac%x_vec(:,2)
                dphi = n_abs*vac%norm(:,1)*vac%x_vec(:,1)**(n_abs-1)*&
                    &vac%x_vec(:,2) + vac%norm(:,2)*vac%x_vec(:,1)**n_abs
                if (present(x_ext)) phi_ext = x_ext(:,1)**n_abs*x_ext(:,2)
            case (3)                                                            ! ring source
                c_free = 0._dp
                allocate(pl_loc(0:max(n_abs,1)),ql_loc(0:max(n_abs,1)))
                do id = 1,vac%n_bnd
                    ierr = ring_pot(vac%x_vec(id,:),phi(id))
                    CHCKERR('')
                    x_pert = vac%x_vec(id,:) + eps_src*vac%norm(id,:)
                    ierr = ring_pot(x_pert,fp)
                    CHCKERR('')
                    x_pert = vac%x_vec(id,:) - eps_src*vac%norm(id,:)
                    ierr = ring_pot(x_pert,fm)
                    CHCKERR('')
                    dphi(id) = (fp-fm)/(2*eps_src)
                end do
                if (present(x_ext)) then
                    do id = 1,size(x_ext,1)
                        ierr = ring_pot(x_ext(id,:),phi_ext(id))
                        CHCKERR('')
                    end do
                end if
            case default
                ierr = 1
                CHCKERR('no such test potential')
        end select
    contains
        ! The n-mode Green's function potential of the ring source,
        ! Q_{|n|-1/2}(gamma)/sqrt(R R_src), gamma = 1 + |x-x_src|^2/(2 R R_src)
        ! (eq. 112 of the 2019 manuscript).
        !> \private
        integer function ring_pot(x,res) result(ierr)
            ! input / output
            real(dp), intent(in) :: x(2)                                        ! point (R, Z)
            real(dp), intent(out) :: res                                        ! potential
            
            ! local variables
            real(dp) :: gam                                                     ! argument of the toroidal function
            
            gam = 1._dp + sum((x-x_src)**2)/(2._dp*x(1)*x_src(1))
            ierr = dtorh1(gam,0,max(n_abs,1),pl_loc,ql_loc,max_n)
            if (ierr.eq.0 .and. max_n.lt.n_abs) ierr = 1
            res = ql_loc(n_abs)/sqrt(x(1)*x_src(1))
        end function ring_pot
    end function set_pot_2
    
    !> Indirect diagonals of the axisymmetric vacuum operators.
    !!
    !! \f$H_{ii}\f$ and \f$G_{ii}\f$ from the least-squares solution, row by
    !! row, of the identities \f$\overline{\text{H}}\phi -
    !! \overline{\text{G}}\,\text{d}\phi = c\,\phi\f$ for the three test
    !! potentials of set_pot_2, with the off-diagonal parts of the products
    !! moved to the right-hand side (each potential scaled to unit maximum).
    !! The counterpart of the identities that the field-line style uses by
    !! default (calc_GH_1). The matrices keep their current diagonals during
    !! the computation.
    integer function diag_indirect_2(vac,dH,dG) result(ierr)
        use vac_vars, only: in_context
        use vac_utilities, only: dis_matvec, dis_diag
        
        character(*), parameter :: rout_name = 'diag_indirect_2'
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(out) :: dH(:)                                          !< indirect diagonal of H (global)
        real(dp), intent(out) :: dG(:)                                          !< indirect diagonal of G (global)
        
        ! local variables
        integer :: kd, rd                                                       ! counters
        real(dp) :: x_src(2)                                                    ! ring source
        real(dp) :: c_free                                                      ! free term of the identity
        real(dp) :: rr                                                          ! residual of an identity at a point
        real(dp) :: det_loc                                                     ! determinant of the normal equations
        real(dp), allocatable :: phi(:), dphi(:)                                ! test potential and its normal data
        real(dp), allocatable :: yH(:), yG(:)                                   ! products with the full H and G
        real(dp), allocatable :: dH0(:), dG0(:)                                 ! current diagonals
        real(dp), allocatable :: a11(:), a12(:), a22(:), b1(:), b2(:)           ! normal equations per row
        
        ierr = 0
        dH = 0._dp
        dG = 0._dp
        if (.not.in_context(vac%ctxt_HG)) return
        allocate(phi(vac%n_bnd),dphi(vac%n_bnd),yH(vac%n_bnd),yG(vac%n_bnd))
        allocate(dH0(vac%n_bnd),dG0(vac%n_bnd))
        allocate(a11(vac%n_bnd),a12(vac%n_bnd),a22(vac%n_bnd),b1(vac%n_bnd),&
            &b2(vac%n_bnd))
        x_src = [sum(vac%x_vec(:,1)),sum(vac%x_vec(:,2))]/vac%n_bnd             ! centroid of the boundary points
        ierr = dis_diag(vac,vac%H,dH0)
        CHCKERR('')
        ierr = dis_diag(vac,vac%G,dG0)
        CHCKERR('')
        a11 = 0._dp
        a12 = 0._dp
        a22 = 0._dp
        b1 = 0._dp
        b2 = 0._dp
        do kd = 1,3
            ierr = set_pot_2(vac,kd,x_src,phi,dphi,c_free)
            CHCKERR('')
            dphi = dphi/maxval(abs(phi))
            phi = phi/maxval(abs(phi))
            ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,yH)
            CHCKERR('')
            ierr = dis_matvec(vac,vac%G,vac%desc_G,dphi,yG)
            CHCKERR('')
            do rd = 1,vac%n_bnd
                ! h phi - g dphi = rr for the unknown diagonals h and g
                rr = c_free*phi(rd) - (yH(rd)-dH0(rd)*phi(rd)) + &
                    &(yG(rd)-dG0(rd)*dphi(rd))
                a11(rd) = a11(rd) + phi(rd)**2
                a12(rd) = a12(rd) - phi(rd)*dphi(rd)
                a22(rd) = a22(rd) + dphi(rd)**2
                b1(rd) = b1(rd) + phi(rd)*rr
                b2(rd) = b2(rd) - dphi(rd)*rr
            end do
        end do
        do rd = 1,vac%n_bnd
            det_loc = a11(rd)*a22(rd) - a12(rd)**2
            dH(rd) = (a22(rd)*b1(rd) - a12(rd)*b2(rd))/det_loc
            dG(rd) = (a11(rd)*b2(rd) - a12(rd)*b1(rd))/det_loc
        end do
    end function diag_indirect_2
    
    !> Sets the diagonal of a distributed matrix of the vacuum to a global
    !! vector.
    subroutine set_diag_vac(vac,A_in,vec)
        use vac_vars, only: loc_idx
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in), pointer :: A_in(:,:)                              !< local matrix
        real(dp), intent(in) :: vec(:)                                          !< global vector of diagonal values
        
        ! local variables
        integer :: i_rd, i_cd                                                   ! index of subrow and column
        integer :: rd                                                           ! global index on the diagonal
        
        do i_rd = 1,size(vac%lims_r,2)
            do i_cd = 1,size(vac%lims_c,2)
                do rd = max(vac%lims_r(1,i_rd),vac%lims_c(1,i_cd)), &
                    &min(vac%lims_r(2,i_rd),vac%lims_c(2,i_cd))
                    A_in(loc_idx(vac%lims_r,i_rd,rd),&
                        &loc_idx(vac%lims_c,i_cd,rd)) = vec(rd)
                end do
            end do
        end do
    end subroutine set_diag_vac
    
    !> Writes a per-point quantity of the boundary with the point's parallel
    !! and field-line indices and its position, for the localization of errors.
    subroutine write_points_vac(vac,file_name,val)
        use num_vars, only: rank, n_procs
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        character(*), intent(in) :: file_name                                   !< file name
        real(dp), intent(in) :: val(:)                                          !< value per boundary point
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: file_i                                                       ! file unit
        integer :: kd                                                           ! counter
        
        if (rank.ne.n_procs-1) return
        open(newunit=file_i,file=file_name,status='replace')
        write(file_i,'(A)') '# id i_par i_line x y z value'
        do id = 1,vac%n_bnd
            write(file_i,'(3I7,4ES15.7)') id, mod(id-1,vac%n_ang(1))+1, &
                &(id-1)/vac%n_ang(1)+1, vac%x_vec(id,:), &
                &[(0._dp, kd=size(vac%x_vec,2)+1,3)], val(id)
        end do
        close(file_i)
    end subroutine write_points_vac
    

    !> Calculates the vacuum response.
    !!
    !! First  \c G and \c  H are completed if  this is not the  first Richardson
    !! level.
    !! 
    !! Then the vacuum contribution is calculated using the two-fold strategy of
    !! first solving 
    !!  \f$\overline{\text{H}}\overline{\Phi} =
    !!  \overline{\text{G}}\overline{\text{E}} \overline{\text{P}}\f$,
    !! followed by left-multiplication of \f$\overline{\Phi}\f$ by
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
    !! For vacuum style  1, the poloidal grid points correspond  to the paralell
    !! grid points and have to be provied by a \c grid variable.
    !!
    !! If \c  jump_to_sol is used for  the current Richardson level,  the vacuum
    !! quantities are not calculated, but just restored.
    !!
    !! Currently, this procedure only works for vacuum style 2 (axisymmetric).
    !!
    !! \return ierr
    integer function calc_vac_res(mds,vac) result(ierr)
        use grid_vars, only: min_par_X, min_alpha, n_alpha
        use MPI
        use rich_vars, only: rich_lvl
        use num_vars, only: jump_to_sol, rich_restart_lvl, eq_style, &
            &use_pol_flux_F, rank, n_procs, check_vac
        use PB3D_ops, only: reconstruct_PB3D_vac
        use X_vars, only: n_mod_X
        use vac_vars, only: set_loc_lims, in_context, loc_idx
        use vac_utilities, only: mat_dis2loc, vac_steps
        use eq_vars, only: vac_perm
        use num_utilities, only: LCM
        use MPI_utilities, only: broadcast_var
        
        character(*), parameter :: rout_name = 'calc_vac_res'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        integer :: jd                                                           ! counter
        integer :: rd, cd                                                       ! global counter for row and col
        integer :: rdl, cdl                                                     ! local counter for row and col
        integer :: i_rd, i_cd                                                   ! index of subrow and subcol
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        integer :: sec_X_loc                                                    ! local sec_X
        integer :: n_loc(2)                                                     ! n_loc for n_mod_X instead of n_bnd
        integer :: lims_cl(2)                                                   ! lims_c_loc in local coordinates
        integer :: par_id                                                       ! parallel index
        integer :: alpha_id                                                     ! field line label
        integer, allocatable :: lims_r_loc(:,:)                                 ! local row limits for n_loc
        integer, allocatable :: lims_c_loc(:,:)                                 ! local column limits for n_loc
        integer, target :: desc_PhiEP(BLACSCTXTSIZE)                            ! descriptor for Phi and EP
        integer, target :: desc_GEP(BLACSCTXTSIZE)                              ! descriptor for GEP
        integer, target :: desc_res2(BLACSCTXTSIZE)                             ! descriptor for res2
        real(dp) :: arg_loc                                                     ! local argument
        real(dp) :: dpar_X                                                      ! step size in par_X 
        real(dp) :: dalpha                                                      ! step size in alpha
        real(dp) :: I_loc                                                       ! local integration rule
        real(dp) :: par_loc                                                     ! local parallel angle
        real(dp) :: alpha_loc                                                   ! local field line label
        real(dp), allocatable, target :: EP(:,:)                                ! EP
        real(dp), allocatable, target :: GEP(:,:)                               ! GEP
        real(dp), allocatable, target :: Phi(:,:)                               ! intermediate Phi
        real(dp), allocatable, target :: res2(:,:)                              ! vacuum response, real version
        real(dp), allocatable, target :: res2_loc(:,:)                          ! local version of res2
        logical :: rank_o                                                       ! output rank
#if ldebug
        integer :: lwork                                                        ! size of work array
        real(dp), allocatable :: res2_ev(:,:)                                   ! eigenvalues of res2
        real(dp), allocatable :: tau(:)                                         ! array tau in upper Hessenberg form
        real(dp), allocatable :: work(:)                                        ! work array
        real(dp), allocatable :: X_plot(:,:)                                    ! X of plot
        real(dp), allocatable :: Y_plot(:,:)                                    ! Y of plot
        complex(dp), allocatable :: res_loc(:,:)                                ! local vac%res
        complex(dp), allocatable :: res_ev(:)                                   ! eigenvalues of res2
        complex(dp), allocatable :: work_c(:)                                   ! work array
        complex(dp), allocatable :: tau_c(:)                                    ! array tau in upper Hessenberg form
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
            return
        end if
        
        call writo('Start calculation of vacuum quantities')
        
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
        
        ! check the operators with exact potentials if requested
        if (check_vac) then
            ierr = check_vac_exact(vac)
            CHCKERR('')
        end if
        
        ! Only for processes that are in the blacs context
        if (in_context(vac%ctxt_HG)) then
            ! user output
#ifdef PB3D_WITH_STRUMPACK
            call writo('Use STRUMPack to solve system',persistent=rank_o)
#else
            call writo('Use ScaLAPACK to solve system',persistent=rank_o)
#endif
            call lvl_ud(1)
            
            ! set step sizes if 3-D vacuum
            if (vac%style.eq.1) call vac_steps(dpar_X,dalpha)
            
            ! the normal data of the modes, EP, and the distributed layout
            call set_up_data()
            CHCKERR('')
#if ldebug
            if (debug_calc_vac_res) then
                call plot_HDF5('EP','EP',reshape(EP,[vac%n_bnd,n_loc(2),1]))
            end if
#endif

            ! the boundary potential of every mode, Phi, and G EP
            call solve_potential()
            CHCKERR('')
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Combine results into vacuum response',persistent=rank_o)
            call lvl_ud(1)
            
            ! the projection of Phi on the data with the integration rule
            call project_response()
            CHCKERR('')
            
#if ldebug
            ! user output
            if (debug_calc_vac_res) then
                call writo('Calculate the eigenvalues of res_loc to &
                    &investigate whether it is positive definite',&
                    &persistent=rank_o)
                call lvl_ud(1)
            end if
#endif
            
            ! the complex response, its factor and its Hermitian part
            if (rank.eq.n_procs-1) then
                call finish_response()
#if ldebug
                allocate(X_plot(n_mod_X,n_mod_X))
                allocate(Y_plot(n_mod_X,n_mod_X))
                do jd = 1,n_mod_X
                    X_plot(jd,:) = vac%sec_X(jd)
                    Y_plot(:, jd) = vac%sec_X(jd)
                end do
                call print_ex_2D(['real vacuum response'],'vac_res_Re',&
                    &rp(vac%res),x=X_plot,draw=.false.)
                call draw_ex(['real vacuum response'],'vac_res_Re',n_mod_X,1,&
                    &.false.)
                call print_ex_2D(['imag vacuum response'],'vac_res_Im',&
                    &ip(vac%res),x=X_plot,draw=.false.)
                call draw_ex(['imag vacuum response'],'vac_res_Im',n_mod_X,1,&
                    &.false.)
                call plot_HDF5(['real part','imag part','imag diff'],'vac_res',&
                    &reshape([rp(vac%res),ip(vac%res),&
                    &ip(vac%res+transpose(vac%res))],[n_mod_X,n_mod_X,1,3]),&
                    &X=reshape(X_plot,[n_mod_X,n_mod_X,1,1]),&
                    &Y=reshape(Y_plot,[n_mod_X,n_mod_X,1,1]))
                deallocate(X_plot,Y_plot)
                
                if (debug_calc_vac_res) then
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
            if (debug_calc_vac_res) then
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
        
            call lvl_ud(-1)
        end if
        
        call lvl_ud(-1)
        
        call writo('Done calculating vacuum quantities')
    contains
        ! The normal data of the modes along the stored normal, -dX/dpar, as the real and imaginary columns of EP, with the distributed layout and descriptors of EP, G EP, Phi and the response.
        !> \private
        subroutine set_up_data()
            ! set secondary mode numbers
            ! (can be already allocated if the vacuum was reconstructed, e.g.
            ! on the HELENA Richardson-restart path)
            if (allocated(vac%sec_X)) deallocate(vac%sec_X)
            allocate(vac%sec_X(n_mod_X))
            if (use_pol_flux_F) then
                vac%sec_X = mds%m(size(mds%m,1),:)
            else
                vac%sec_X = mds%n(size(mds%n,1),:)
            end if
            
            ! set sizes
            do jd = 1,2
                n_loc(jd) = numroc(2*n_mod_X,vac%bs,vac%ind_p(jd),0,vac%n_p(jd))
            end do
            call set_loc_lims(n_loc(1),vac%bs,vac%ind_p(1),vac%n_p(1),&
                &lims_r_loc)
            call set_loc_lims(n_loc(2),vac%bs,vac%ind_p(2),vac%n_p(2),&
                &lims_c_loc)
            
            ! initialize descriptors
            call descinit(desc_PhiEP,vac%n_bnd,2*n_mod_X,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for EP')
            call descinit(desc_GEP,vac%n_bnd,2*n_mod_X,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for GEP')
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
                    cdl = loc_idx(lims_c_loc,i_cd,cd)
                    
                    ! set local secondary mode number
                    sec_X_loc = vac%sec_X(mod(cd-1,n_mod_X)+1)
                    
                    ! iterate over all subrows of matrix with n_bnd rows
                    subrows: do i_rd = 1,size(vac%lims_r,2)
                        row: do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                            ! set local row index
                            rdl = loc_idx(vac%lims_r,i_rd,rd)
                            
                            ! set parallel index and field line label
                            par_id = mod(rd-1,vac%n_ang(1))+1
                            alpha_id = (rd-1)/vac%n_ang(1) + 1
                            
                            ! set up local argument
                            select case (vac%style)
                                case (1)                                        ! field-line 3-D
                                    par_loc = pi*min_par_X + &
                                        &(par_id-1)*dpar_X
                                    alpha_loc = pi*min_alpha + &
                                        &(alpha_id-1)*dalpha
                                    if (use_pol_flux_F) then
                                        arg_loc = vac%prim_X*alpha_loc + &
                                            &(vac%prim_X*vac%jq-sec_X_loc)*&
                                            &par_loc
                                    else
                                        arg_loc = vac%prim_X*alpha_loc + &
                                            &(sec_X_loc-vac%prim_X*vac%jq)*&
                                            &par_loc
                                    end if
                                case (2)                                        ! axisymmetric
                                    arg_loc = -sec_X_loc*&
                                        &vac%ang(par_id,alpha_id)
                            end select
                            
                            ! set up argument
                            
                            ! save EP: the normal data of the vacuum
                            ! potential along the stored normal -J nabla psi,
                            ! out of the vacuum, dphi = -J B.grad X = -dX/dpar
                            ! = -i (nq-m) X for a mode (nq-m) (the imaginary
                            ! unit is carried by the real and imaginary
                            ! columns)
                            if (cd.le.n_mod_X) then                             ! real part of rhs
                                EP(rdl,cdl) = cos(arg_loc)                      ! rp(i exp)
                            else
                                EP(rdl,cdl) = sin(arg_loc)                      ! ip(i exp)
                            end if
                            if (use_pol_flux_F) then
                                EP(rdl,cdl) = -(vac%prim_X*vac%jq-sec_X_loc)*&
                                    &EP(rdl,cdl)
                            else
                                EP(rdl,cdl) = -(sec_X_loc-vac%prim_X*vac%jq)*&
                                    &EP(rdl,cdl)
                            end if
                        end do row
                    end do subrows
                end do col
            end do subcols
        end subroutine set_up_data
        
        ! The boundary potential of every mode from the boundary element system, and the product G EP.
        !> \private
        subroutine solve_potential()
            ! solve for Phi
            ierr = solve_Phi_BEM(vac,EP,Phi,[vac%n_bnd,2*n_mod_X],&
                &[vac%n_loc(1),n_loc(2)],lims_c_loc,desc_PhiEP)
            CHCKERR('')
            
            ! calculate GEP
            call pdgemm('N','N',vac%n_bnd,2*n_mod_X,vac%n_bnd,1._dp,vac%G,1,1,&
                &vac%desc_G,EP,1,1,desc_PhiEP,0._dp,GEP,1,1,desc_GEP)
        end subroutine solve_potential
        
        ! The projection of Phi on the data with the integration rule (the trapezoidal weights along the lines; the average over the lines for the field-line style), gathered on the last process.
        !> \private
        subroutine project_response()
            ! Convert EP into IEP where I contains the integration rule:
            !         (θ_{i+1}-θ_{i-1})/2   for i = 2..n-1
            !   I_i = (θ_2-θ_1}/2           for i = 1
            !         (θ_{n}-θ_{n-1}}/2     for i = n
            ! with n = vac%n_bnd.
            ! For vacua of type 1, the step  size is constant. Also, there is an
            ! additional step size in the alpha direction.
            subrows3: do i_rd = 1,size(vac%lims_r,2)
                row3: do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                    ! set local row index
                    rdl = loc_idx(vac%lims_r,i_rd,rd)
                    
                    subcols3: do i_cd = 1,size(lims_c_loc,2)
                        ! set local column limits
                        lims_cl = [loc_idx(lims_c_loc,i_cd,lims_c_loc(1,i_cd)), &
                        &loc_idx(lims_c_loc,i_cd,lims_c_loc(2,i_cd))]
                        
                        select case (vac%style)
                            case (1)                                            ! field-line 3-D
                                ! the projection integrates over the field
                                ! lines as well: dalpha divided by the alpha
                                ! span gives the average over the lines, the
                                ! per-toroidal-mode normalization of the
                                ! axisymmetric style
                                I_loc = dpar_X/n_alpha
                                if (mod(rd-1,vac%n_ang(1))+1.eq.1 .or. &        ! first point on field line
                                    &mod(rd-1,vac%n_ang(1))+1.eq.vac%n_ang(1)) &! last point on field line
                                    &I_loc = I_loc*0.5_dp
                            case (2)                                            ! axisymmetric
                                I_loc = 0.5_dp*&
                                    &(vac%ang(min(rd+1,vac%n_bnd),1)-&
                                    &vac%ang(max(rd-1,1),1))
                        end select
                        
                        EP(rdl,lims_cl(1):lims_cl(2)) = &
                            &EP(rdl,lims_cl(1):lims_cl(2)) * I_loc
                    end do subcols3
                end do row3
            end do subrows3
            
            ! calculate (EP)^T Phi
            ! Note that this means that the second  half of the rows in EP^T are
            ! still missing a  minus sign if they are to  represent EP^†. This
            ! will be taken into account when res2_loc is converted to vac%res.
            call pdgemm('T','N',2*n_mod_X,2*n_mod_X,vac%n_bnd,1._dp,EP,1,1,&
                &desc_PhiEP,Phi,1,1,desc_PhiEP,0._dp,res2,1,1,desc_res2)
            
            ! get serial version on last process
            ierr = mat_dis2loc(vac%ctxt_HG,res2,lims_r_loc,lims_c_loc,res2_loc,&
                &proc=n_procs-1)
            CHCKERR('')
        end subroutine project_response
        
        ! The response as a complex matrix, divided by the permeability, multiplied by 2 pi (the toroidal integral) and made Hermitian.
        !> \private
        subroutine finish_response()
            ! convert to complex numbers and save in vac%res
                vac%res = res2_loc(1:n_mod_X,1:n_mod_X)-&
                    &(-res2_loc(n_mod_X+1:2*n_mod_X,n_mod_X+1:2*n_mod_X))       ! minus sign for lower rows
                vac%res = vac%res + iu * &
                    &(res2_loc(1:n_mod_X,n_mod_X+1:2*n_mod_X)-&
                    &res2_loc(n_mod_X+1:2*n_mod_X,1:n_mod_X))                   ! minus sign for lower rows
                vac%res = vac%res/vac_perm
                ! Bring  res to  the surface-energy  Hessian as  it enters  the
                ! potential-energy  matrix  A,  so  that  the  free  boundary
                ! condition adds  X_e^H res X_e  (set_BC_2, set_BC_4)  and POST
                ! adds the same term to  E_pot. The matrix A corresponds to 2 dW
                ! with the  toroidal integral included  (alpha_fac = 2  pi for a
                ! single field line). With the normal data along the stored
                ! normal -J nabla psi, out of the vacuum, the vacuum energy is
                ! dW_vac = 1/(2 mu_0) int dphi^* phi dS >= 0, so res up to here
                ! is the positive Neumann-to-Dirichlet map per unit toroidal
                ! angle and 2 dW_vac = 2 pi X_e^H res X_e. In the
                ! large-aspect-ratio cylinder  the vacuum  energy of  an edge
                ! displacement equals the field-line bending energy of a plasma
                ! layer of  thickness a/|m|  at the edge, 4 pi^2 (nq-m)^2 |X_e|^2
                ! / (mu_0 R_0 |m|) in A (full-stack test vac_response_cylinder).
                ! Checked against MISHKA on the circular external kink (issue
                ! #4, step 3).
                vac%res = 2._dp*pi*vac%res
                
                ! The exact  operator is Hermitian;  the discretization leaves
                ! an asymmetry  of about 1E-5  that makes SLEPC reject  tightly
                ! converged eigenvalues as complex. Keep the Hermitian part.
                vac%res = 0.5_dp*(vac%res+conjg(transpose(vac%res)))
        end subroutine finish_response
        
    end function calc_vac_res
    
    !> Checks the style-1 (3-D) vacuum operators on the equilibrium surface
    !! itself with potentials whose boundary data are exact (issue #4: the
    !! checks T5 and T6 of the validation catalogue on the real surface, not
    !! on the analytic torus of the full-stack tests).
    !!
    !! For a constant, the three Cartesian coordinates and \f$x^2 - y^2\f$,
    !! harmonic inside the surface, and for two potentials harmonic outside
    !! it and decaying, a point source on the magnetic axis at \f$\zeta = 0\f$
    !! and a ring of point sources along the whole magnetic axis, with unit
    !! weights and with weights \f$\cos(n \zeta)\f$ for the primary mode number
    !! \f$n\f$ (the toroidal structure of the modes), the residual
    !! of the Green identity \f$\overline{\text{H}} \phi -
    !! \overline{\text{G}} \partial_n \phi - c \phi\f$ is printed for
    !! \f$c = -4\pi, -2\pi, 0, 2\pi, 4\pi\f$, scaled by \f$\max |4 \pi \phi|\f$
    !! (interior potentials satisfy the identity with \f$c = -4\pi\f$,
    !! exterior ones with \f$c = 0\f$, as calc_GH_1 builds the diagonal; see
    !! the full-stack test greens_identity_3d). Then the exterior problem is
    !! solved for the two exterior potentials with solve_Phi_BEM() and the
    !! boundary trace is compared with the exact one. The point source is
    !! the sharper test when the axis comes close to the surface, the ring
    !! varies on the scale of the minor radius everywhere.
    !!
    !! Every number is printed on a line starting with CHECKVAC so that it can
    !! be collected from the log, in lines short enough for writo, which cuts
    !! beyond about 100 characters. Switched on with the input \c check_vac.
    !! For the field-line style (check_vac_exact_1): the Cartesian harmonics and
    !! the axis sources below; for the axisymmetric style (check_vac_exact_2):
    !! the potentials of set_pot_2, the ring-source solve, and the assembled
    !! diagonals against the indirect ones of diag_indirect_2.
    integer function check_vac_exact(vac) result(ierr)
        character(*), parameter :: rout_name = 'check_vac_exact'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        select case (vac%style)
            case (1)
                ierr = check_vac_exact_1(vac)
                CHCKERR('')
            case (2)
                ierr = check_vac_exact_2(vac)
                CHCKERR('')
        end select
    end function check_vac_exact
    
    !> Field-line style of check_vac_exact.
    integer function check_vac_exact_1(vac) result(ierr)
        use num_vars, only: rank, n_procs, vac_diag_style
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: vec_dis2loc, n_col1_loc, vec_glob2dis, &
            &dis_matvec, calc_diag_direct, dis_diag, vac_steps
        use VMEC_vars, only: R_V_c, R_V_s, Z_V_c, Z_V_s, mn_V, mnmax_V

        character(*), parameter :: rout_name = 'check_vac_exact_1'

        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables

        ! local variables
        integer, parameter :: n_pot = 8                                         ! number of test potentials
        integer, parameter :: n_ring = 256                                      ! sources along the magnetic axis
        integer, parameter :: n_dir = 8                                         ! boundary points of the exterior check
        real(dp), parameter :: f_ext(6) = [-0.2_dp,0.05_dp,0.1_dp,0.2_dp,&
            &0.5_dp,1._dp]                                                      ! their distances, in units of half the larger width
        integer :: id, jd, kd, ld                                               ! counters
        integer :: n_src                                                        ! number of sources of a potential
        integer :: n_bnd                                                        ! number of boundary points
        integer :: n_col1                                                       ! local number of columns of a global column vector
        integer :: desc_RPhi(BLACSCTXTSIZE)                                     ! descriptor for right-hand side and solution
        integer, allocatable :: lims_c_R(:,:)                                   ! local column limits of right-hand side and solution
        logical :: rank_o                                                       ! output rank
        real(dp) :: x_ring(n_ring,3)                                            ! sources along the magnetic axis
        real(dp) :: w_ring(n_ring)                                              ! their weights
        real(dp) :: zeta_k                                                      ! toroidal angle of a source
        real(dp) :: d_min                                                       ! smallest distance of the surface to the sources
        real(dp) :: res_off(-2:2)                                               ! identity residuals for offsets c = k 2 pi
        real(dp) :: err_max, err_2                                              ! errors of the exterior solution
        real(dp), allocatable :: phi(:), dphi(:)                                ! potential and norm . grad phi
        real(dp), allocatable :: lhs(:), rhs(:)                                 ! H phi and G dphi
        real(dp), allocatable :: R_mat(:,:), Phi_mat(:,:)                       ! distributed right-hand side and solution
        real(dp), allocatable :: phi_num(:)                                     ! numerical boundary trace
        real(dp), allocatable :: phi_num_ext(:,:)                               ! numerical traces of the three exterior potentials
        real(dp), allocatable :: x_ext(:,:), f_d(:)                             ! points of the exterior check and their distances
        real(dp), allocatable :: R_b(:)                                         ! major radius of the boundary points
        real(dp), allocatable :: Phi_glob(:,:), D_glob(:,:), val_ex(:,:)        ! boundary columns and exact values of the exterior check
        real(dp) :: scale_ext(10)                                               ! scales of its errors
        real(dp) :: a_w                                                         ! half the larger width of the boundary
        character(len=16) :: col_names(10)                                      ! names of its columns
        character(len=max_str_ln) :: pot_name                                   ! name of test potential
        real(dp) :: c_exp                                                       ! expected offset of the identity
        real(dp) :: dpar, dalpha                                                ! steps of the field-line grid
        real(dp), allocatable :: diag_rs(:), diag_dir(:), diag_fac(:)           ! row-sum, smooth direct and faceted direct diagonal of H
        real(dp), allocatable :: I_sing(:), omega_v(:)                          ! cell integral and vacuum's solid angle
        real(dp), allocatable :: d_diag(:), d_fac(:)                            ! differences per point
        logical, allocatable :: avail(:)                                        ! direct value available

        ! initialize ierr
        ierr = 0

        if (.not.in_context(vac%ctxt_HG)) return

        rank_o = rank.eq.n_procs-1
        n_bnd = vac%n_bnd

        ! sources on the magnetic axis: the m = 0 modes of R and Z at the
        ! first radial point, cos(m theta - n zeta) convention, placed with
        ! the same toroidal angle as the boundary points (R_V_c is normalized
        ! like them); source 1 is at zeta = 0
        do kd = 1,n_ring
            zeta_k = 2._dp*pi*(kd-1)/n_ring
            x_ring(kd,:) = 0._dp
            do id = 1,mnmax_V
                if (mn_V(id,1).eq.0) then
                    x_ring(kd,1) = x_ring(kd,1) + &
                        &R_V_c(id,1,0)*cos(mn_V(id,2)*zeta_k) - &
                        &R_V_s(id,1,0)*sin(mn_V(id,2)*zeta_k)
                    x_ring(kd,3) = x_ring(kd,3) + &
                        &Z_V_c(id,1,0)*cos(mn_V(id,2)*zeta_k) - &
                        &Z_V_s(id,1,0)*sin(mn_V(id,2)*zeta_k)
                end if
            end do
            x_ring(kd,2) = x_ring(kd,1)*sin(zeta_k)
            x_ring(kd,1) = x_ring(kd,1)*cos(zeta_k)
        end do

        call writo('Check of the vacuum operators with exact potentials',&
            &persistent=rank_o)
        call lvl_ud(1)
        call writo('CHECKVAC source n_bnd='//trim(i2str(n_bnd))//' x_src='//&
            &trim(r2strt(x_ring(1,1)))//' '//trim(r2strt(x_ring(1,2)))//' '//&
            &trim(r2strt(x_ring(1,3))),persistent=rank_o)

        allocate(phi(n_bnd),dphi(n_bnd),lhs(n_bnd),rhs(n_bnd),phi_num(n_bnd))

        ! Green identities
        do kd = 1,n_pot
            call set_pot(kd)
            ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
            CHCKERR('')
            ierr = dis_matvec(vac,vac%G,vac%desc_G,dphi,rhs)
            CHCKERR('')
            lhs = lhs - rhs
            do id = -2,2
                res_off(id) = maxval(abs(lhs - id*2._dp*pi*phi))/&
                    &maxval(abs(4._dp*pi*phi))
            end do
            ! (two lines: writo cuts lines beyond about 100 characters)
            call writo('CHECKVAC identity '//trim(pot_name)//' n_bnd='//&
                &trim(i2str(n_bnd))//' res(c=-4pi,-2pi,0)='//&
                &trim(r2strt(res_off(-2)))//' '//trim(r2strt(res_off(-1)))//&
                &' '//trim(r2strt(res_off(0))),persistent=rank_o)
            call writo('CHECKVAC identity '//trim(pot_name)//' n_bnd='//&
                &trim(i2str(n_bnd))//' res(c=2pi,4pi)='//&
                &trim(r2strt(res_off(1)))//' '//trim(r2strt(res_off(2))),&
                &persistent=rank_o)
            c_exp = 0._dp
            if (kd.le.5) c_exp = -4._dp*pi                                      ! interior harmonics
            call write_points_vac(vac,'check_vac_identity_'//trim(pot_name)//'.dat',&
                &abs(lhs-c_exp*phi)/maxval(abs(4._dp*pi*phi)))
        end do

        ! direct diagonals of H next to the row-sum one: the cell integral of
        ! the kernel (the corrected eq. 69 of the 2019 manuscript) with the
        ! corner term of the smooth surface, -2 pi, or with minus the vacuum's
        ! exact solid angle at the faceted node (its eq. 45); the difference
        ! from the row sum measures the quadrature error of the off-diagonal
        ! row that the row sum absorbs into the diagonal
        if (vac%style.eq.1) then
            allocate(diag_rs(n_bnd),diag_dir(n_bnd),diag_fac(n_bnd),&
                &I_sing(n_bnd),omega_v(n_bnd),d_diag(n_bnd),d_fac(n_bnd),&
                &avail(n_bnd))
            ierr = dis_diag(vac,vac%H,diag_rs)
            CHCKERR('')
            call vac_steps(dpar,dalpha)
            ierr = calc_diag_direct(vac,dpar,dalpha,diag_dir,diag_fac,I_sing,&
                &omega_v,avail)
            CHCKERR('')
            d_diag = 0._dp
            where (avail) d_diag = abs(diag_rs-diag_dir)/(4._dp*pi)
            d_fac = 0._dp
            where (avail) d_fac = abs(diag_rs-diag_fac)/(4._dp*pi)
            ! (four lines: writo cuts lines beyond about 100 characters)
            call writo('CHECKVAC diagonal n_bnd='//trim(i2str(n_bnd))//&
                &' n_points='//trim(i2str(count(avail)))//&
                &' vac_diag_style='//trim(i2str(vac_diag_style)),&
                &persistent=rank_o)
            call writo('CHECKVAC diagonal max|assembled-direct|/4pi='//&
                &trim(r2strt(maxval(d_diag)))//' mean='//&
                &trim(r2strt(sum(d_diag)/max(1,count(avail)))),&
                &persistent=rank_o)
            call writo('CHECKVAC diagonal_parts max|assembled-faceted|/4pi='//&
                &trim(r2strt(maxval(d_fac))),persistent=rank_o)
            call writo('CHECKVAC diagonal_parts max|I_sing|/4pi='//&
                &trim(r2strt(maxval(abs(I_sing)/(4._dp*pi))))//&
                &' max|omega_v/4pi-1/2|='//&
                &trim(r2strt(maxval(abs(omega_v/(4._dp*pi)-0.5_dp),&
                &mask=avail))),persistent=rank_o)
            call write_points_vac(vac,'check_vac_diag_diff.dat',d_diag)
            call write_points_vac(vac,'check_vac_diag_diff_facet.dat',d_fac)
            call write_points_vac(vac,'check_vac_diag_Ising.dat',abs(I_sing)/(4._dp*pi))
            call write_points_vac(vac,'check_vac_diag_omega.dat',omega_v/(4._dp*pi))
            deallocate(diag_rs,diag_dir,diag_fac,I_sing,omega_v,d_diag,d_fac,&
                &avail)
        end if
        

        ! exterior solves
        n_col1 = n_col1_loc(vac)
        call set_loc_lims(n_col1,vac%bs,vac%ind_p(2),vac%n_p(2),lims_c_R)
        allocate(R_mat(max(1,vac%n_loc(1)),n_col1))
        allocate(Phi_mat(max(1,vac%n_loc(1)),n_col1))
        call descinit(desc_RPhi,n_bnd,1,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
            &max(1,vac%n_loc(1)),ierr)
        CHCKERR('descinit failed')
        allocate(phi_num_ext(n_bnd,3))
        do kd = n_pot-2,n_pot
            call set_pot(kd)
            if (n_col1.gt.0) call vec_glob2dis(vac,dphi,R_mat(:,1))
            Phi_mat = 0._dp
            ierr = solve_Phi_BEM(vac,R_mat,Phi_mat,[n_bnd,1],&
                &[vac%n_loc(1),n_col1],lims_c_R,desc_RPhi)
            CHCKERR('')
            ierr = vec_dis2loc(vac%ctxt_HG,&
                &reshape(Phi_mat(1:vac%n_loc(1),:),[vac%n_loc(1)*n_col1]),&
                &vac%lims_r,phi_num)
            CHCKERR('')
            phi_num_ext(:,kd-n_pot+3) = phi_num
            err_max = maxval(abs(phi_num-phi))/maxval(abs(phi))
            err_2 = sqrt(sum((phi_num-phi)**2)/sum(phi**2))
            call writo('CHECKVAC solve '//trim(pot_name)//' n_bnd='//&
                &trim(i2str(n_bnd))//' min_dist='//trim(r2strt(d_min))//&
                &' err_max='//trim(r2strt(err_max))//' err_2='//&
                &trim(r2strt(err_2)),persistent=rank_o)
            call write_points_vac(vac,'check_vac_solve_'//trim(pot_name)//'.dat',&
                &abs(phi_num-phi)/maxval(abs(phi)))
        end do
        
        ! the potential off the boundary (check_vac_ext), from the exact
        ! boundary values and from the solved ones, on the normals of eight
        ! boundary points spread over the field lines and along them; inward
        ! (negative distance) the representation of an exterior potential is
        ! zero and that of the constant -1
        allocate(x_ext(n_dir*size(f_ext),3),f_d(n_dir*size(f_ext)))
        allocate(R_b(n_bnd))
        R_b = sqrt(vac%x_vec(:,1)**2+vac%x_vec(:,2)**2)
        a_w = max(maxval(R_b)-minval(R_b),&
            &maxval(vac%x_vec(:,3))-minval(vac%x_vec(:,3)))/2
        do kd = 1,n_dir
            id = ((kd-1)*vac%n_ang(2))/n_dir*vac%n_ang(1) + 1 + &
                &((2*kd-1)*(vac%n_ang(1)-1))/(2*n_dir)
            do jd = 1,size(f_ext)
                ld = (kd-1)*size(f_ext)+jd
                f_d(ld) = f_ext(jd)
                x_ext(ld,:) = vac%x_vec(id,:) - f_ext(jd)*a_w*&
                    &vac%norm(id,:)/sqrt(sum(vac%norm(id,:)**2))
            end do
        end do
        allocate(Phi_glob(n_bnd,10),D_glob(n_bnd,10),&
            &val_ex(n_dir*size(f_ext),10))
        do kd = n_pot-2,n_pot
            call set_pot(kd)
            jd = 2*(kd-n_pot+2)+1
            Phi_glob(:,jd) = phi
            Phi_glob(:,jd+1) = phi_num_ext(:,kd-n_pot+3)
            D_glob(:,jd) = dphi
            D_glob(:,jd+1) = dphi
            col_names(jd) = trim(pot_name)//'_ex'
            col_names(jd+1) = trim(pot_name)//'_sol'
            scale_ext(jd:jd+1) = maxval(abs(phi))
            n_src = n_ring
            if (kd.eq.n_pot-2) n_src = 1
            do ld = 1,size(f_d)
                val_ex(ld,jd) = 0._dp
                if (f_d(ld).gt.0._dp) val_ex(ld,jd) = sum(w_ring(1:n_src)/&
                    &sqrt(sum((spread(x_ext(ld,:),1,n_src)-&
                    &x_ring(1:n_src,:))**2,2)))
            end do
            val_ex(:,jd+1) = val_ex(:,jd)
        end do
        Phi_glob(:,7) = 1._dp
        D_glob(:,7) = 0._dp
        val_ex(:,7) = 0._dp
        where (f_d.lt.0._dp) val_ex(:,7) = -1._dp
        col_names(7) = 'constant'
        scale_ext(7) = 1._dp
        do kd = 1,3                                                             ! x, y, z: interior harmonics, -x inside
            Phi_glob(:,7+kd) = vac%x_vec(:,kd)
            D_glob(:,7+kd) = vac%norm(:,kd)
            val_ex(:,7+kd) = 0._dp
            where (f_d.lt.0._dp) val_ex(:,7+kd) = -x_ext(:,kd)
            scale_ext(7+kd) = max(maxval(abs(vac%x_vec(:,kd))),1.E-12_dp)
        end do
        col_names(8:10) = ['x','y','z']
        ierr = check_vac_ext(vac,f_d,x_ext,Phi_glob,D_glob,val_ex,col_names,&
            &scale_ext,7,cols_x=[8,9,10])
        CHCKERR('')

        call lvl_ud(-1)
    contains
        
        ! Sets phi and dphi = norm . grad phi for test potential kd.
        !> \private
        subroutine set_pot(kd)
            integer, intent(in) :: kd                                           ! which potential

            integer :: id, ld                                                   ! counters
            integer :: n_src                                                    ! number of sources
            real(dp) :: r_src                                                   ! distance to a source

            d_min = huge(1._dp)
            do id = 1,n_bnd
                select case (kd)
                    case (1)                                                    ! constant (interior harmonic)
                        pot_name = 'constant'
                        phi(id) = 1._dp
                        dphi(id) = 0._dp
                    case (2)                                                    ! x (interior harmonic)
                        pot_name = 'x'
                        phi(id) = vac%x_vec(id,1)
                        dphi(id) = vac%norm(id,1)
                    case (3)                                                    ! y (interior harmonic)
                        pot_name = 'y'
                        phi(id) = vac%x_vec(id,2)
                        dphi(id) = vac%norm(id,2)
                    case (4)                                                    ! z (interior harmonic)
                        pot_name = 'z'
                        phi(id) = vac%x_vec(id,3)
                        dphi(id) = vac%norm(id,3)
                    case (5)                                                    ! x^2 - y^2 (interior harmonic)
                        pot_name = 'x2-y2'
                        phi(id) = vac%x_vec(id,1)**2 - vac%x_vec(id,2)**2
                        dphi(id) = 2._dp*(vac%x_vec(id,1)*vac%norm(id,1) - &
                            &vac%x_vec(id,2)*vac%norm(id,2))
                    case (6,7,8)                                                ! point source(s) inside: exterior harmonic
                        select case (kd)
                            case (6)
                                pot_name = 'point_source'
                                n_src = 1
                                w_ring(1) = 1._dp
                            case (7)
                                pot_name = 'axis_ring'
                                n_src = n_ring
                                w_ring = 1._dp
                            case (8)
                                pot_name = 'axis_ring_n'
                                n_src = n_ring
                                w_ring = [(cos(vac%prim_X*2._dp*pi*(ld-1)/&
                                    &n_ring),ld=1,n_ring)]
                        end select
                        phi(id) = 0._dp
                        dphi(id) = 0._dp
                        do ld = 1,n_src
                            r_src = sqrt(sum((vac%x_vec(id,:)-x_ring(ld,:))**2))
                            d_min = min(d_min,r_src)
                            phi(id) = phi(id) + w_ring(ld)/r_src
                            dphi(id) = dphi(id) - w_ring(ld)*sum(vac%norm(id,:)*&
                                &(vac%x_vec(id,:)-x_ring(ld,:)))/r_src**3
                        end do
                end select
            end do
        end subroutine set_pot
    end function check_vac_exact_1

    
    !> Axisymmetric style of check_vac_exact.
    integer function check_vac_exact_2(vac) result(ierr)
        use num_vars, only: rank, n_procs, vac_diag_style
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: vec_dis2loc, n_col1_loc, vec_glob2dis, &
            &dis_matvec, dis_diag
        
        character(*), parameter :: rout_name = 'check_vac_exact_2'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: n_bnd                                                        ! number of boundary points
        integer :: n_col1                                                       ! local number of columns of a global column vector
        integer :: desc_RPhi(BLACSCTXTSIZE)                                     ! descriptor for right-hand side and solution
        integer, allocatable :: lims_c_R(:,:)                                   ! local column limits of right-hand side and solution
        logical :: rank_o                                                       ! output rank
        real(dp) :: x_src(2)                                                    ! ring source
        real(dp) :: c_free                                                      ! free term of the identity
        real(dp) :: res_off(-2:2)                                               ! identity residuals for offsets c = k 2 pi
        real(dp) :: err_max, err_2                                              ! errors of the exterior solution
        real(dp), allocatable :: phi(:), dphi(:)                                ! potential and norm . grad phi
        real(dp), allocatable :: lhs(:), rhs(:)                                 ! H phi and G dphi
        real(dp), allocatable :: R_mat(:,:), Phi_mat(:,:)                       ! distributed right-hand side and solution
        real(dp), allocatable :: phi_num(:)                                     ! numerical boundary trace
        real(dp), allocatable :: dH0(:), dG0(:), dH(:), dG(:)                   ! assembled and indirect diagonals
        integer, parameter :: n_dir = 8                                         ! boundary points of the exterior check
        real(dp), parameter :: f_ext(6) = [-0.2_dp,0.05_dp,0.1_dp,0.2_dp,&
            &0.5_dp,1._dp]                                                      ! their distances, in units of half the larger width
        integer :: jd, ld                                                       ! counters
        real(dp), allocatable :: x_ext(:,:), f_d(:), phi_ext(:)                 ! points of the exterior check, their distances and a potential there
        real(dp), allocatable :: Phi_glob(:,:), D_glob(:,:), val_ex(:,:)        ! boundary columns and exact values of the exterior check
        real(dp) :: scale_ext(3)                                                ! scales of its errors
        real(dp) :: a_w                                                         ! half the larger width of the boundary
        character(len=max_str_ln) :: pot_name                                   ! name of test potential
        
        ! initialize ierr
        ierr = 0
        
        if (.not.in_context(vac%ctxt_HG)) return
        rank_o = rank.eq.n_procs-1
        n_bnd = vac%n_bnd
        x_src = [sum(vac%x_vec(:,1)),sum(vac%x_vec(:,2))]/n_bnd                 ! centroid of the boundary points
        
        call writo('Check of the axisymmetric vacuum operators with exact &
            &potentials',persistent=rank_o)
        call lvl_ud(1)
        call writo('CHECKVAC source n_bnd='//trim(i2str(n_bnd))//' x_src='//&
            &trim(r2strt(x_src(1)))//' '//trim(r2strt(x_src(2))),&
            &persistent=rank_o)
        allocate(phi(n_bnd),dphi(n_bnd),lhs(n_bnd),rhs(n_bnd),phi_num(n_bnd))
        
        ! Green identities
        do kd = 1,3
            ierr = set_pot_2(vac,kd,x_src,phi,dphi,c_free)
            CHCKERR('')
            select case (kd)
                case (1)
                    pot_name = 'Rn'
                case (2)
                    pot_name = 'RnZ'
                case (3)
                    pot_name = 'ring_source'
            end select
            ierr = dis_matvec(vac,vac%H,vac%desc_H,phi,lhs)
            CHCKERR('')
            ierr = dis_matvec(vac,vac%G,vac%desc_G,dphi,rhs)
            CHCKERR('')
            lhs = lhs - rhs
            do id = -2,2
                res_off(id) = maxval(abs(lhs - id*2._dp*pi*phi))/&
                    &maxval(abs(4._dp*pi*phi))
            end do
            ! (two lines: writo cuts lines beyond about 100 characters)
            call writo('CHECKVAC identity '//trim(pot_name)//' n_bnd='//&
                &trim(i2str(n_bnd))//' res(c=-4pi,-2pi,0)='//&
                &trim(r2strt(res_off(-2)))//' '//trim(r2strt(res_off(-1)))//&
                &' '//trim(r2strt(res_off(0))),persistent=rank_o)
            call writo('CHECKVAC identity '//trim(pot_name)//' n_bnd='//&
                &trim(i2str(n_bnd))//' res(c=2pi,4pi)='//&
                &trim(r2strt(res_off(1)))//' '//trim(r2strt(res_off(2))),&
                &persistent=rank_o)
            call write_points_vac(vac,'check_vac_identity_'//trim(pot_name)//&
                &'.dat',abs(lhs-c_free*phi)/maxval(abs(4._dp*pi*phi)))
        end do
        
        ! exterior solve for the ring source (phi and dphi are its data)
        n_col1 = n_col1_loc(vac)
        call set_loc_lims(n_col1,vac%bs,vac%ind_p(2),vac%n_p(2),lims_c_R)
        allocate(R_mat(max(1,vac%n_loc(1)),n_col1))
        allocate(Phi_mat(max(1,vac%n_loc(1)),n_col1))
        call descinit(desc_RPhi,n_bnd,1,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
            &max(1,vac%n_loc(1)),ierr)
        CHCKERR('descinit failed')
        if (n_col1.gt.0) call vec_glob2dis(vac,dphi,R_mat(:,1))
        Phi_mat = 0._dp
        ierr = solve_Phi_BEM(vac,R_mat,Phi_mat,[n_bnd,1],&
            &[vac%n_loc(1),n_col1],lims_c_R,desc_RPhi)
        CHCKERR('')
        ierr = vec_dis2loc(vac%ctxt_HG,&
            &reshape(Phi_mat(1:vac%n_loc(1),:),[vac%n_loc(1)*n_col1]),&
            &vac%lims_r,phi_num)
        CHCKERR('')
        err_max = maxval(abs(phi_num-phi))/maxval(abs(phi))
        err_2 = sqrt(sum((phi_num-phi)**2)/sum(phi**2))
        call writo('CHECKVAC solve ring_source n_bnd='//trim(i2str(n_bnd))//&
            &' err_max='//trim(r2strt(err_max))//' err_2='//&
            &trim(r2strt(err_2)),persistent=rank_o)
        call write_points_vac(vac,'check_vac_solve_ring_source.dat',&
            &abs(phi_num-phi)/maxval(abs(phi)))
        
        ! the assembled diagonals against the indirect ones from the
        ! identities (equal by construction when vac_diag_style is 1; the
        ! difference measures the quadrature error of the off-diagonal rows
        ! when the assembled diagonals are the direct ones)
        allocate(dH0(n_bnd),dG0(n_bnd),dH(n_bnd),dG(n_bnd))
        ierr = dis_diag(vac,vac%H,dH0)
        CHCKERR('')
        ierr = dis_diag(vac,vac%G,dG0)
        CHCKERR('')
        ierr = diag_indirect_2(vac,dH,dG)
        CHCKERR('')
        ! (three lines: writo cuts lines beyond about 100 characters)
        call writo('CHECKVAC diagonal n_bnd='//trim(i2str(n_bnd))//&
            &' vac_diag_style='//trim(i2str(vac_diag_style)),&
            &persistent=rank_o)
        call writo('CHECKVAC diagonal max|H_assembled-H_indirect|/4pi='//&
            &trim(r2strt(maxval(abs(dH0-dH))/(4._dp*pi)))//' mean='//&
            &trim(r2strt(sum(abs(dH0-dH))/(4._dp*pi*n_bnd))),&
            &persistent=rank_o)
        call writo('CHECKVAC diagonal max|G_assembled-G_indirect|/max|G|='//&
            &trim(r2strt(maxval(abs(dG0-dG))/maxval(abs(dG0)))),&
            &persistent=rank_o)
        call write_points_vac(vac,'check_vac_diag_diff.dat',&
            &abs(dH0-dH)/(4._dp*pi))
        call write_points_vac(vac,'check_vac_diag_diff_G.dat',&
            &abs(dG0-dG)/maxval(abs(dG0)))
        
        ! the potential off the boundary (check_vac_ext): the ring source from
        ! its exact boundary values and from the solved ones, and R^n, on the
        ! normals of eight boundary points; inward (negative distance) the
        ! representation of the ring source is zero and that of R^n is -R^n
        allocate(x_ext(n_dir*size(f_ext),2),f_d(n_dir*size(f_ext)),&
            &phi_ext(n_dir*size(f_ext)))
        a_w = max(maxval(vac%x_vec(:,1))-minval(vac%x_vec(:,1)),&
            &maxval(vac%x_vec(:,2))-minval(vac%x_vec(:,2)))/2
        do kd = 1,n_dir
            id = 1 + ((kd-1)*(n_bnd-1))/n_dir
            do jd = 1,size(f_ext)
                ld = (kd-1)*size(f_ext)+jd
                f_d(ld) = f_ext(jd)
                x_ext(ld,:) = vac%x_vec(id,:) - f_ext(jd)*a_w*&
                    &vac%norm(id,:)/sqrt(sum(vac%norm(id,:)**2))
            end do
        end do
        allocate(Phi_glob(n_bnd,3),D_glob(n_bnd,3),val_ex(n_dir*size(f_ext),3))
        ierr = set_pot_2(vac,3,x_src,phi,dphi,c_free,x_ext,phi_ext)
        CHCKERR('')
        Phi_glob(:,1) = phi
        Phi_glob(:,2) = phi_num
        D_glob(:,1) = dphi
        D_glob(:,2) = dphi
        val_ex(:,1) = 0._dp
        where (f_d.gt.0._dp) val_ex(:,1) = phi_ext
        val_ex(:,2) = val_ex(:,1)
        scale_ext(1:2) = maxval(abs(phi))
        ierr = set_pot_2(vac,1,x_src,phi,dphi,c_free,x_ext,phi_ext)
        CHCKERR('')
        Phi_glob(:,3) = phi
        D_glob(:,3) = dphi
        val_ex(:,3) = 0._dp
        where (f_d.lt.0._dp) val_ex(:,3) = -phi_ext
        scale_ext(3) = maxval(abs(phi))
        ierr = check_vac_ext(vac,f_d,x_ext,Phi_glob,D_glob,val_ex,&
            &['ring_ex ','ring_sol','Rn      '],scale_ext,3)
        CHCKERR('')
        
        call lvl_ud(-1)
    end function check_vac_exact_2

    !> \public Potential at points off the boundary, from its values and its
    !! normal data on the boundary.
    !!
    !! The kernels are taken along  the stored normal \f$-\mathcal{J}\nabla\psi\f$,
    !! which points out  of the vacuum. Green's second identity  then gives, at a
    !! point \f$\vec{x}\f$ off the boundary,
    !!  \f[ \sum_j \left( H(\vec{x},\vec{x}_j)\phi_j -
    !!  G(\vec{x},\vec{x}_j)\,\text{d}\phi_j \right) = \left\{
    !!  \begin{array}{ll}
    !!  4\pi\phi(\vec{x}) & \text{in the vacuum,} \\
    !!  0 & \text{in the plasma}
    !!  \end{array} \right. \f]
    !! for a potential \f$\phi\f$ that is  harmonic in the vacuum and decays at
    !! infinity, and \f$0\f$ in the  vacuum and \f$-4\pi\phi(\vec{x})\f$ in the
    !! plasma for a potential  harmonic inside the boundary (the constant for the
    !! field-line style, \f$R^{|n|}\f$ for the axisymmetric one).  On the boundary
    !! the same sum, with the diagonal  that calc_GH() adds there, is the identity
    !! \f$\overline{\text{H}}\phi - \overline{\text{G}}\,\text{d}\phi = c\phi\f$
    !! of solve_Phi_BEM(), \f$c = 0\f$ and \f$-4\pi\f$ for the two kinds.
    !!
    !! This routine returns the sum divided by \f$4\pi\f$ for every column of
    !! \c Phi_s and \c D_s. The rows  of G and H are assembled by calc_GH() for
    !! the external points  (no diagonal; the trapezoidal rule with the near-field
    !! rule along the field lines for the field-line style, the toroidal-function
    !! kernels for the axisymmetric style), in chunks that bound the memory.
    !! Within a few grid spacings of the boundary the quadrature does not
    !! resolve the kernels.
    !!
    !! Every process of the BLACS grid of the vacuum has to call it.
    !!
    !! \return ierr
    integer function eval_pot_ext(vac,x_ext,Phi_s,D_s,desc_s,pot) result(ierr)
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: mat_dis2loc

        character(*), parameter :: rout_name = 'eval_pot_ext'

        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        real(dp), intent(in) :: x_ext(:,:)                                      !< points: (R,Z) for the axisymmetric style, (x,y,z) for the field-line style; the same on every process
        real(dp), intent(in) :: Phi_s(:,:)                                      !< local part of the boundary potential, n_bnd by n_col
        real(dp), intent(in) :: D_s(:,:)                                        !< local part of its normal data along the stored normal
        integer, intent(in) :: desc_s(BLACSCTXTSIZE)                            !< descriptor of \c Phi_s and \c D_s
        real(dp), intent(inout) :: pot(:,:)                                     !< potential at the points, n_ext by n_col, on every process of the BLACS grid

        ! local variables
        integer, parameter :: max_entries = 20000000                            ! largest block of external rows of G or H
        integer :: n_ext                                                        ! number of points
        integer :: n_col                                                        ! number of columns
        integer :: n_ch                                                         ! points per chunk
        integer :: i_1, i_2, n_r                                                ! first and last point of a chunk and their number
        integer :: n_row_loc, n_col_loc                                         ! local numbers of rows and columns
        integer :: desc_GH(BLACSCTXTSIZE)                                       ! descriptor of the external rows of G and H
        integer :: desc_P(BLACSCTXTSIZE)                                        ! descriptor of the potential of a chunk
        integer, allocatable :: lims_r_loc(:,:)                                 ! local row limits of a chunk
        integer, allocatable :: lims_c_loc(:,:)                                 ! local column limits
        real(dp), allocatable :: G_loc(:,:), H_loc(:,:)                         ! external rows of G and H
        real(dp), allocatable :: P_loc(:,:)                                     ! local part of the potential of a chunk
        real(dp), allocatable :: P_ser(:,:)                                     ! potential of a chunk, gathered

        ! initialize ierr
        ierr = 0

        n_ext = size(x_ext,1)
        n_col = size(pot,2)
        pot = 0._dp
        n_ch = max(1,min(n_ext,max_entries/max(1,vac%n_bnd)))

        ! local columns of the result
        n_col_loc = 0
        if (in_context(vac%ctxt_HG)) n_col_loc = &
            &numroc(n_col,vac%bs,vac%ind_p(2),0,vac%n_p(2))
        call set_loc_lims(n_col_loc,vac%bs,vac%ind_p(2),vac%n_p(2),lims_c_loc)

        ! chunks of points
        i_1 = 1
        do while (i_1.le.n_ext)
            i_2 = min(n_ext,i_1+n_ch-1)
            n_r = i_2-i_1+1
            n_row_loc = 0
            if (in_context(vac%ctxt_HG)) n_row_loc = &
                &numroc(n_r,vac%bs,vac%ind_p(1),0,vac%n_p(1))
            call set_loc_lims(n_row_loc,vac%bs,vac%ind_p(1),vac%n_p(1),&
                &lims_r_loc)

            ! external rows of G and H
            allocate(G_loc(max(1,n_row_loc),max(1,vac%n_loc(2))))
            allocate(H_loc(max(1,n_row_loc),max(1,vac%n_loc(2))))
            G_loc = 0._dp
            H_loc = 0._dp
            ierr = calc_GH(vac,n_r_in=n_r,lims_r_in=lims_r_loc,&
                &x_vec_in=x_ext(i_1:i_2,:),G_in=G_loc,H_in=H_loc)
            CHCKERR('')

            ! (H Phi - G D) / 4 pi, gathered
            if (in_context(vac%ctxt_HG)) then
                call descinit(desc_GH,n_r,vac%n_bnd,vac%bs,vac%bs,0,0,&
                    &vac%ctxt_HG,max(1,n_row_loc),ierr)
                CHCKERR('descinit failed for GH')
                call descinit(desc_P,n_r,n_col,vac%bs,vac%bs,0,0,&
                    &vac%ctxt_HG,max(1,n_row_loc),ierr)
                CHCKERR('descinit failed for P')
                allocate(P_loc(max(1,n_row_loc),max(1,n_col_loc)))
                P_loc = 0._dp
                call pdgemm('N','N',n_r,n_col,vac%n_bnd,1._dp/(4*pi),H_loc,&
                    &1,1,desc_GH,Phi_s,1,1,desc_s,0._dp,P_loc,1,1,desc_P)
                call pdgemm('N','N',n_r,n_col,vac%n_bnd,-1._dp/(4*pi),G_loc,&
                    &1,1,desc_GH,D_s,1,1,desc_s,1._dp,P_loc,1,1,desc_P)
                allocate(P_ser(n_r,n_col))
                ierr = mat_dis2loc(vac%ctxt_HG,P_loc(1:n_row_loc,1:n_col_loc),&
                    &lims_r_loc,lims_c_loc,P_ser)
                CHCKERR('')
                pot(i_1:i_2,:) = P_ser
                deallocate(P_loc,P_ser)
            end if

            deallocate(G_loc,H_loc,lims_r_loc)
            i_1 = i_2+1
        end do
    end function eval_pot_ext

    !> Copies columns  with the boundary points  as rows, the same  on every
    !! process, into their local block-cyclic part.
    subroutine glob2dis_cols(vac,A_glob,lims_c_loc,A_dis)
        use vac_vars, only: in_context, loc_idx

        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: A_glob(:,:)                                     !< global columns, n_bnd by n_col
        integer, intent(in) :: lims_c_loc(:,:)                                  !< local column limits
        real(dp), intent(inout) :: A_dis(:,:)                                   !< local part

        ! local variables
        integer :: i_rd, i_cd                                                   ! index of subrow and subcolumn
        integer :: rd, cd                                                       ! global row and column

        if (.not.in_context(vac%ctxt_HG)) return
        do i_cd = 1,size(lims_c_loc,2)
            do cd = lims_c_loc(1,i_cd),lims_c_loc(2,i_cd)
                do i_rd = 1,size(vac%lims_r,2)
                    do rd = vac%lims_r(1,i_rd),vac%lims_r(2,i_rd)
                        A_dis(loc_idx(vac%lims_r,i_rd,rd),&
                            &loc_idx(lims_c_loc,i_cd,cd)) = A_glob(rd,cd)
                    end do
                end do
            end do
        end do
    end subroutine glob2dis_cols
    
    !> Subtracts from the potential off the boundary the representation of
    !! the interior harmonic that equals it at the nearest boundary point.
    !!
    !! In the vacuum the representation of a potential that is harmonic inside
    !! the boundary is zero (eval_pot_ext()), so the subtraction changes
    !! nothing exactly. In the quadrature it removes the error that the double
    !! layer makes with the potential at the nearest boundary point, the part
    !! that grows without bound as the point approaches the boundary; what
    !! remains comes from the variation of the potential along the boundary and
    !! from the single layer. With \c psi_bnd the boundary values of the
    !! interior harmonic (the constant for the field-line style, \f$R^{|n|}\f$
    !! for the axisymmetric one) and \c pot_psi its representation, every
    !! column \f$c\f$ becomes \f$\phi(\vec{x},c) -
    !! \Phi_j^c/\psi_j\,\phi_\psi(\vec{x})\f$, \f$j\f$ the boundary point
    !! nearest to \f$\vec{x}\f$. For points in the vacuum only: inside the
    !! plasma the result is the potential of the nearest boundary point
    !! instead of zero.
    !!
    !! For the field-line style, with the gradient \c grad of the columns at
    !! the boundary points (grad_bnd_1()) and the representations \c pot_x of
    !! the coordinates x, y and z, the subtracted harmonic is the linear
    !! function \f$\Phi_j + \vec{g}_j\cdot(\vec{x}-\vec{x}_j)\f$, which
    !! matches the potential to first order at the nearest boundary point and
    !! removes the error of the potential's variation along the boundary as
    !! well. The constant alone leaves that error where the nearest boundary
    !! point lies on a neighbouring field line, a large part of a wavelength
    !! away in a mode of high toroidal mode number.
    subroutine sub_nearest(vac,x_ext,Phi_glob,psi_bnd,pot_psi,pot,grad,pot_x)
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: x_ext(:,:)                                      !< points
        real(dp), intent(in) :: Phi_glob(:,:)                                   !< boundary values of the columns, n_bnd by n_c
        real(dp), intent(in) :: psi_bnd(:)                                      !< boundary values of the interior harmonic
        real(dp), intent(in) :: pot_psi(:)                                      !< its representation at the points
        real(dp), intent(inout) :: pot(:,:)                                     !< representation of the columns at the points, n_ext by n_c
        real(dp), intent(in), optional :: grad(:,:,:)                           !< gradient of the columns at the boundary points, n_bnd by 3 by n_c (field-line style)
        real(dp), intent(in), optional :: pot_x(:,:)                            !< representation of x, y and z at the points, n_ext by 3
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: j_min                                                        ! nearest boundary point
        
        do id = 1,size(x_ext,1)
            j_min = minloc(sum((vac%x_vec-spread(x_ext(id,:),1,vac%n_bnd))**2,&
                &2),1)
            pot(id,:) = pot(id,:) - Phi_glob(j_min,:)/psi_bnd(j_min)*pot_psi(id)
            if (present(grad) .and. present(pot_x)) then
                do kd = 1,3
                    pot(id,:) = pot(id,:) - grad(j_min,kd,:)*&
                        &(pot_x(id,kd)-vac%x_vec(j_min,kd)*pot_psi(id))
                end do
            end if
        end do
    end subroutine sub_nearest
    
    !> Gradient of columns of boundary values at the boundary points of the
    !! field-line style.
    !!
    !! At every point the gradient \f$\vec{g}\f$ of a potential follows from
    !! \f$\vec{g}\cdot\Delta\vec{x}_p = \Delta\Phi_p\f$, \f$\vec{g}\cdot
    !! \Delta\vec{x}_\alpha = \Delta\Phi_\alpha\f$ and \f$\vec{g}\cdot
    !! \vec{n} = \text{d}\phi\f$: central differences of the position and of
    !! the boundary values along the field line (one-sided at its ends) and
    !! across the field lines at the same parallel point (periodic in the
    !! field-line label when the lines cover a full turn, \c alpha_periodic,
    !! one-sided otherwise), and the normal data along the stored normal.
    subroutine grad_bnd_1(vac,Phi_glob,D_glob,alpha_periodic,grad)
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        real(dp), intent(in) :: Phi_glob(:,:)                                   !< boundary values, n_bnd by n_c
        real(dp), intent(in) :: D_glob(:,:)                                     !< normal data, n_bnd by n_c
        logical, intent(in) :: alpha_periodic                                   !< the field lines cover a full turn in alpha
        real(dp), intent(out) :: grad(:,:,:)                                    !< gradient, n_bnd by 3 by n_c
        
        ! local variables
        integer :: rd                                                           ! boundary point
        integer :: par_id, alpha_id                                             ! its parallel index and field line
        integer :: n_par, n_al                                                  ! points per line and lines
        integer :: jp(2), ja(2)                                                 ! neighbours along and across the lines
        integer :: ipiv(3), info                                                ! pivots and status of the solve
        real(dp) :: A(3,3)                                                      ! system
        real(dp), allocatable :: b(:,:)                                         ! right-hand sides and solution
        
        n_par = vac%n_ang(1)
        n_al = vac%n_ang(2)
        allocate(b(3,size(Phi_glob,2)))
        do rd = 1,vac%n_bnd
            par_id = mod(rd-1,n_par)+1
            alpha_id = (rd-1)/n_par+1
            jp = [rd-1,rd+1]
            if (par_id.eq.1) jp(1) = rd
            if (par_id.eq.n_par) jp(2) = rd
            ja = [rd-n_par,rd+n_par]
            if (alpha_id.eq.1) then
                ja(1) = rd
                if (alpha_periodic) ja(1) = rd+(n_al-1)*n_par
            end if
            if (alpha_id.eq.n_al) then
                ja(2) = rd
                if (alpha_periodic) ja(2) = rd-(n_al-1)*n_par
            end if
            if (ja(1).eq.ja(2)) ja(1) = rd                                      ! two lines
            A(1,:) = vac%x_vec(jp(2),:)-vac%x_vec(jp(1),:)
            A(2,:) = vac%x_vec(ja(2),:)-vac%x_vec(ja(1),:)
            A(3,:) = vac%norm(rd,:)
            b(1,:) = Phi_glob(jp(2),:)-Phi_glob(jp(1),:)
            b(2,:) = Phi_glob(ja(2),:)-Phi_glob(ja(1),:)
            b(3,:) = D_glob(rd,:)
            call dgesv(3,size(b,2),A,3,ipiv,b,3,info)
            if (info.ne.0) b = 0._dp                                            ! degenerate: no first-order term
            grad(rd,:,:) = b
        end do
    end subroutine grad_bnd_1

    !> Checks the potential that eval_pot_ext() gives off the boundary against
    !! exact values, for check_vac_exact().
    !!
    !! Every column holds a potential's values \c Phi_glob and normal data \c
    !! D_glob on the boundary and its exact values \c val_ex at the points \c
    !! x_ext. The points lie on the normals of some boundary points, at the
    !! distance \c f_d times half the larger width of the boundary in R and Z,
    !! outward for positive \c f_d and inward for negative. For each column and
    !! distance the largest error over the points, divided by \c scale, is
    !! printed on a line starting with CHECKVAC exterior; for the points in the
    !! vacuum also the error after sub_nearest() with the interior harmonic
    !! of column \c col_psi (sub=). All values go to check_vac_exterior.dat.
    !!
    !! \return ierr
    integer function check_vac_ext(vac,f_d,x_ext,Phi_glob,D_glob,val_ex,&
        &col_names,scale,col_psi,cols_x) result(ierr)
        use grid_vars, only: min_alpha, max_alpha
        use num_vars, only: rank, n_procs
        use vac_vars, only: in_context, set_loc_lims

        character(*), parameter :: rout_name = 'check_vac_ext'

        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        real(dp), intent(in) :: f_d(:)                                          !< distance of each point from the boundary, in units of half the larger width
        real(dp), intent(in) :: x_ext(:,:)                                      !< points
        real(dp), intent(in) :: Phi_glob(:,:)                                   !< boundary values, n_bnd by n_col
        real(dp), intent(in) :: D_glob(:,:)                                     !< normal data, n_bnd by n_col
        real(dp), intent(in) :: val_ex(:,:)                                     !< exact values at the points, n_ext by n_col
        character(len=*), intent(in) :: col_names(:)                            !< names of the columns
        real(dp), intent(in) :: scale(:)                                        !< scale of the errors, per column
        integer, intent(in) :: col_psi                                          !< column of the interior harmonic
        integer, intent(in), optional :: cols_x(3)                              !< columns of x, y and z (field-line style): first-order subtraction as well

        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: n_col, n_col_loc                                             ! number of columns and local number
        integer :: file_i                                                       ! file unit
        integer :: desc_s(BLACSCTXTSIZE)                                        ! descriptor of the boundary columns
        integer, allocatable :: lims_c_loc(:,:)                                 ! local column limits
        real(dp), allocatable :: Phi_s(:,:), D_s(:,:)                           ! local boundary columns
        real(dp), allocatable :: pot(:,:)                                       ! potential at the points
        real(dp), allocatable :: pot_sub(:,:)                                   ! the same after sub_nearest
        real(dp), allocatable :: pot_sub1(:,:)                                  ! the same after the first-order sub_nearest
        real(dp), allocatable :: grad(:,:,:)                                    ! gradient of the columns on the boundary
character(len=max_str_ln) :: sub_str                                    ! error after sub_nearest as string
        logical, allocatable :: done(:)                                         ! distance already printed
        logical, allocatable :: same_d(:)                                       ! points at the distance being printed
        logical :: rank_o                                                       ! output rank

        ! initialize ierr
        ierr = 0

        rank_o = rank.eq.n_procs-1
        n_col = size(Phi_glob,2)

        ! boundary columns
        n_col_loc = 0
        if (in_context(vac%ctxt_HG)) n_col_loc = &
            &numroc(n_col,vac%bs,vac%ind_p(2),0,vac%n_p(2))
        call set_loc_lims(n_col_loc,vac%bs,vac%ind_p(2),vac%n_p(2),lims_c_loc)
        allocate(Phi_s(max(1,vac%n_loc(1)),max(1,n_col_loc)))
        allocate(D_s(max(1,vac%n_loc(1)),max(1,n_col_loc)))
        Phi_s = 0._dp
        D_s = 0._dp
        desc_s = 0
        call glob2dis_cols(vac,Phi_glob,lims_c_loc,Phi_s)
        call glob2dis_cols(vac,D_glob,lims_c_loc,D_s)
        if (in_context(vac%ctxt_HG)) then
            call descinit(desc_s,vac%n_bnd,n_col,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed')
        end if

        ! potential at the points, and after subtraction of the interior
        ! harmonic matching the nearest boundary value
        allocate(pot(size(x_ext,1),n_col))
        ierr = eval_pot_ext(vac,x_ext,Phi_s,D_s,desc_s,pot)
        CHCKERR('')
        pot_sub = pot
        call sub_nearest(vac,x_ext,Phi_glob,Phi_glob(:,col_psi),&
            &pot(:,col_psi),pot_sub)
        if (present(cols_x)) then
            allocate(grad(vac%n_bnd,3,n_col))
            call grad_bnd_1(vac,Phi_glob,D_glob,&
                &abs(max_alpha-min_alpha-2._dp).lt.1.E-8_dp,grad)
            pot_sub1 = pot
            call sub_nearest(vac,x_ext,Phi_glob,Phi_glob(:,col_psi),&
                &pot(:,col_psi),pot_sub1,grad=grad,pot_x=pot(:,cols_x))
        end if

        ! largest error per column and distance (one line each: writo cuts
        ! long lines)
        allocate(done(size(f_d)),same_d(size(f_d)))
        done = .false.
        do id = 1,size(f_d)
            if (done(id)) cycle
            same_d = abs(f_d-f_d(id)).lt.1.E-12_dp
            done = done .or. same_d
            do jd = 1,n_col
                sub_str = ''
                if (f_d(id).gt.0._dp .and. jd.ne.col_psi) sub_str = ' sub='//&
                    &trim(r2strt(maxval(abs(pot_sub(:,jd)-val_ex(:,jd)),&
                    &mask=same_d)/scale(jd)))
                call writo('CHECKVAC exterior '//trim(col_names(jd))//&
                    &' n_bnd='//trim(i2str(vac%n_bnd))//' d/a='//&
                    &trim(r2strt(f_d(id)))//' err='//trim(r2strt(&
                    &maxval(abs(pot(:,jd)-val_ex(:,jd)),mask=same_d)/&
                    &scale(jd)))//trim(sub_str),persistent=rank_o)
                if (present(cols_x)) then                                       ! (no short-circuit .and. in Fortran)
                    if (f_d(id).gt.0._dp .and. jd.ne.col_psi .and. &
                        &all(cols_x.ne.jd)) &
                        &call writo('CHECKVAC exterior1 '//&
                        &trim(col_names(jd))//' n_bnd='//&
                        &trim(i2str(vac%n_bnd))//' d/a='//&
                        &trim(r2strt(f_d(id)))//' sub1='//trim(r2strt(&
                        &maxval(abs(pot_sub1(:,jd)-val_ex(:,jd)),mask=same_d)/&
                        &scale(jd))),persistent=rank_o)
                end if
            end do
        end do

        ! all values
        if (rank_o) then
            open(newunit=file_i,file='check_vac_exterior.dat',&
                &status='replace')
            write(file_i,'(A)') '# d/a, point, then per column: numerical, &
                &after sub_nearest and exact potential'
            do id = 1,size(f_d)
                write(file_i,'(*(ES16.8))') f_d(id), x_ext(id,:), &
                    &(pot(id,jd), pot_sub(id,jd), val_ex(id,jd), jd = 1,n_col)
            end do
            close(file_i)
        end if
    end function check_vac_ext

    !> \public Plots  the vacuum  potential of  computed eigenmodes  outside the
    !! plasma.
    !!
    !! For every  eigenmode of \c X_ids  the normal data of  the vacuum
    !! potential on the boundary are those of calc_vac_res(): along the stored
    !! normal \f$-\mathcal{J}\nabla\psi\f$, out of the vacuum, \f$\text{d}\phi =
    !! -\mathcal{J}\vec{B}\cdot\nabla X = -i(nq-m)X\f$, summed over the modes of
    !! the solution at the plasma edge with the basis \f$e^{i(n\alpha +
    !! (nq-m)\theta)} = e^{i(n\zeta_\text{F} - m\theta)}\f$ of the field-line
    !! style, or \f$e^{-im\theta}\f$ for the axisymmetric style, which splits
    !! off \f$e^{in\zeta_\text{F}}\f$. These data are \f$-\mathcal{J}Q_n\f$, with
    !! \f$Q_n\f$ the normal component of the magnetic perturbation that POST
    !! plots: the normal component of the perturbed field is continuous across
    !! the boundary. The boundary potential follows from the exterior system of
    !! solve_Phi_BEM(), the potential in the vacuum from Green's representation
    !! (eval_pot_ext()): for the field-line style the complex potential at the
    !! 3-D points, for the axisymmetric style \f$\hat{\phi}(R,Z)\f$ times
    !! \f$e^{in\zeta_\text{F}}\f$.
    !!
    !! The points lie in the planes of constant \f$\zeta_\text{F}\f$ of the
    !! plasma plots, \c n_zeta_plot values from \c min_zeta_plot to \c
    !! max_zeta_plot in units of \f$\pi\f$ (the geometrical angle is
    !! \f$-\zeta_\text{F}\f$ for VMEC and HELENA), on a grid of \c n_vac_plot
    !! points in \f$R\f$ and \f$Z\f$ from \c min_Rvac_plot to \c max_Rvac_plot
    !! and from \c min_Zvac_plot to \c max_Zvac_plot, in metres. Where a maximum
    !! does not exceed its minimum (the default), the range is that of the
    !! boundary points widened on both sides by half the larger of their widths
    !! in \f$R\f$ and \f$Z\f$.
    !!
    !! Written for every eigenmode: the real and imaginary parts of the
    !! potential at time zero, <tt>[id]_vac_pot_RE</tt> and
    !! <tt>[id]_vac_pot_IM</tt> (HDF5), and the boundary points with their
    !! normal data and potential, <tt>Data/[id]_vac_bnd.dat</tt>. Written once:
    !! <tt>vac_in_plasma</tt>, 1 inside the plasma and 0 in the vacuum (for the
    !! field-line style the side of the normal of the nearest boundary point
    !! close to it and minus the representation of the constant farther out,
    !! for the axisymmetric style the boundary polygon); the potential is set
    !! to zero there. Near the boundary the quadrature does not resolve the
    !! kernels; sub_nearest() removes the largest part of that error
    !! (check_vac_exact() prints the error against the distance, with and
    !! without it).
    !!
    !! \return ierr
    integer function vac_pot_plot(grid_sol,vac,sol,X_ids) result(ierr)
        use num_vars, only: min_zeta_plot, max_zeta_plot, min_Rvac_plot, &
            &max_Rvac_plot, min_Zvac_plot, max_Zvac_plot, n_vac_plot, &
            &n_zeta_plot, n_procs, rank, use_pol_flux_F, use_normalization, &
            &data_dir
        use eq_vars, only: R_0
        use sol_vars, only: sol_type
        use MPI_utilities, only: broadcast_var
        use vac_vars, only: in_context, set_loc_lims
        use vac_utilities, only: mat_dis2loc, vac_steps
        use X_vars, only: n_mod_X
        use grid_vars, only: min_alpha, max_alpha, min_par_X, max_par_X
        use rich_vars, only: n_par_X
        
        character(*), parameter :: rout_name = 'vac_pot_plot'

        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        integer, intent(in) :: X_ids(:)                                         !< indices of the eigenmodes

        ! local variables
        integer :: id, jd, kd, ld, pd                                           ! counters
        integer :: n_X                                                          ! number of eigenmodes
        integer :: n_col                                                        ! boundary columns: real and imaginary part per eigenmode (and the constant)
        integer :: n_col_loc                                                    ! local number of columns
        integer :: n_RZ                                                         ! points per plane
        integer :: n_ext                                                        ! points evaluated
        integer :: rd                                                           ! boundary point
        integer :: par_id, alpha_id                                             ! parallel index and field line of a boundary point
        integer :: file_i                                                       ! file unit
        integer :: desc_s(BLACSCTXTSIZE)                                        ! descriptor of the boundary columns
        integer, allocatable :: lims_c_loc(:,:)                                 ! local column limits
        real(dp) :: L_fac                                                       ! unit of length of the plot: R_0 when normalized
        real(dp) :: lims_R(2), lims_Z(2)                                        ! ranges in R and Z, normalized
        real(dp) :: w_box                                                       ! widening of the default ranges
        real(dp) :: dpar_X, dalpha                                              ! steps of the field-line grid
        real(dp) :: arg_loc                                                     ! phase of a mode at a boundary point
        real(dp), allocatable :: R_b(:), Z_b(:)                                 ! R and Z of the boundary points
        real(dp), allocatable :: psi_bnd(:)                                     ! interior harmonic on the boundary
        real(dp), allocatable :: grad(:,:,:)                                    ! gradient of the eigenmode columns on the boundary (field-line style)
logical, allocatable :: in_ext(:)                                       ! points in the plasma
        integer :: j_min, j_nb                                                  ! boundary point nearest to a point and its neighbour on the line
        real(dp), allocatable :: R_pl(:), Z_pl(:), zeta_pl(:)                   ! R, Z and zeta_F of the plot grid
        real(dp), allocatable :: x_ext(:,:)                                     ! points
        real(dp), allocatable :: D_s(:,:), Phi_s(:,:)                           ! local normal data and boundary potential
        real(dp), allocatable :: D_ser(:,:), Phi_ser(:,:)                       ! the same, gathered
        real(dp), allocatable :: pot(:,:)                                       ! potential at the points
        real(dp), allocatable :: f_plot(:,:,:,:)                                ! real and imaginary part on the plot grid
        real(dp), allocatable :: in_pl(:,:,:)                                   ! 1 inside the plasma, 0 outside
        real(dp), allocatable :: XYZ(:,:,:,:)                                   ! Cartesian coordinates of the plot grid
        complex(dp), allocatable :: sol_vec(:,:)                                ! edge solution, modes by eigenmodes
        complex(dp), allocatable :: dphi(:,:)                                   ! normal data, boundary points by eigenmodes
        complex(dp) :: phi_c                                                    ! complex potential
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: rank_o                                                       ! output rank

        ! initialize ierr
        ierr = 0

        n_X = size(X_ids)
        if (n_X.eq.0) return
        rank_o = rank.eq.n_procs-1

        call writo('Vacuum potential of '//trim(i2str(n_X))//' eigenmodes')
        call lvl_ud(1)

        if (.not.use_pol_flux_F) then
            ierr = 1
            CHCKERR('The vacuum potential plot needs the poloidal flux')
        end if

        ! field-line style: the grid variables that vac_steps reads, here and
        ! in calc_GH_1, from the parallel angle of the vacuum output
        if (vac%style.eq.1) then
            if (maxval(abs(vac%ang)).le.0._dp) then
                ierr = 1
                err_msg = 'No parallel angle in the vacuum output: rerun PB3D'
                CHCKERR(err_msg)
            end if
            min_par_X = vac%ang(1,1)/pi
            max_par_X = vac%ang(vac%n_ang(1),1)/pi
            n_par_X = vac%n_ang(1)
            call vac_steps(dpar_X,dalpha)
        end if

        ! edge solution of the eigenmodes, on every process
        allocate(sol_vec(n_mod_X,n_X))
        sol_vec = 0._dp
        do pd = 1,n_X
            if (rank_o) sol_vec(:,pd) = sol%vec(:,grid_sol%loc_n_r,X_ids(pd))
            ierr = broadcast_var(sol_vec(:,pd),n_procs-1)
            CHCKERR('')
        end do

        ! normal data at the boundary points along the stored normal, as in
        ! calc_vac_res: -i (nq-m) X for every mode
        allocate(dphi(vac%n_bnd,n_X))
        dphi = 0._dp
        do rd = 1,vac%n_bnd
            par_id = mod(rd-1,vac%n_ang(1))+1
            alpha_id = (rd-1)/vac%n_ang(1)+1
            do jd = 1,n_mod_X
                if (vac%style.eq.1) then                                        ! field-line 3-D
                    arg_loc = vac%prim_X*(pi*min_alpha+(alpha_id-1)*dalpha) + &
                        &(vac%prim_X*vac%jq-vac%sec_X(jd))*&
                        &vac%ang(par_id,alpha_id)
                else                                                            ! axisymmetric
                    arg_loc = -vac%sec_X(jd)*vac%ang(par_id,alpha_id)
                end if
                dphi(rd,:) = dphi(rd,:) - iu*(vac%prim_X*vac%jq-vac%sec_X(jd))*&
                    &exp(iu*arg_loc)*sol_vec(jd,:)
            end do
        end do

        ! the boundary operators
        ierr = calc_GH(vac)
        CHCKERR('')

        ! boundary columns: real and imaginary part of every eigenmode, and an
        ! interior harmonic (the constant for the field-line style, R^|n| for
        ! the axisymmetric one), whose representation corrects the potential
        ! near the boundary (sub_nearest) and, for the field-line style, marks
        ! the plasma
        ! (and for the field-line style the coordinates x, y and z, whose
        ! representations make the correction first order)
        n_col = 2*n_X+1
        if (vac%style.eq.1) n_col = n_col+3
        allocate(D_ser(vac%n_bnd,n_col),Phi_ser(vac%n_bnd,n_col))
        allocate(psi_bnd(vac%n_bnd))
        D_ser = 0._dp
        Phi_ser = 0._dp
        do pd = 1,n_X
            D_ser(:,2*pd-1) = rp(dphi(:,pd))
            D_ser(:,2*pd) = ip(dphi(:,pd))
        end do
        if (vac%style.eq.1) then                                                ! field-line 3-D
            psi_bnd = 1._dp
            D_ser(:,2*n_X+2:2*n_X+4) = vac%norm
        else                                                                    ! axisymmetric
            psi_bnd = vac%x_vec(:,1)**abs(vac%prim_X)
            D_ser(:,2*n_X+1) = abs(vac%prim_X)*vac%norm(:,1)*&
                &vac%x_vec(:,1)**(abs(vac%prim_X)-1)
        end if
        n_col_loc = 0
        if (in_context(vac%ctxt_HG)) n_col_loc = &
            &numroc(n_col,vac%bs,vac%ind_p(2),0,vac%n_p(2))
        call set_loc_lims(n_col_loc,vac%bs,vac%ind_p(2),vac%n_p(2),lims_c_loc)
        allocate(D_s(max(1,vac%n_loc(1)),max(1,n_col_loc)))
        allocate(Phi_s(max(1,vac%n_loc(1)),max(1,n_col_loc)))
        D_s = 0._dp
        Phi_s = 0._dp
        desc_s = 0
        call glob2dis_cols(vac,D_ser,lims_c_loc,D_s)
        if (in_context(vac%ctxt_HG)) then
            call descinit(desc_s,vac%n_bnd,n_col,vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed')

            ! boundary potential from the exterior system
            ierr = solve_Phi_BEM(vac,D_s,Phi_s,[vac%n_bnd,n_col],&
                &[vac%n_loc(1),n_col_loc],lims_c_loc,desc_s)
            CHCKERR('')

            ! gathered on every process of the grid
            ierr = mat_dis2loc(vac%ctxt_HG,&
                &Phi_s(1:vac%n_loc(1),1:n_col_loc),vac%lims_r,lims_c_loc,&
                &Phi_ser)
            CHCKERR('')
        end if
        Phi_ser(:,2*n_X+1) = psi_bnd                                            ! the interior harmonics: set, not solved
        if (vac%style.eq.1) Phi_ser(:,2*n_X+2:2*n_X+4) = vac%x_vec
        call glob2dis_cols(vac,Phi_ser,lims_c_loc,Phi_s)

        ! plot grid, normalized like the boundary points
        L_fac = 1._dp
        if (use_normalization) L_fac = R_0
        if (vac%style.eq.1) then                                                ! field-line 3-D
            R_b = sqrt(vac%x_vec(:,1)**2+vac%x_vec(:,2)**2)
            Z_b = vac%x_vec(:,3)
        else                                                                    ! axisymmetric
            R_b = vac%x_vec(:,1)
            Z_b = vac%x_vec(:,2)
        end if
        w_box = max(maxval(R_b)-minval(R_b),maxval(Z_b)-minval(Z_b))/2
        if (max_Rvac_plot.gt.min_Rvac_plot) then
            lims_R = [min_Rvac_plot,max_Rvac_plot]/L_fac
        else
            lims_R = [max(minval(R_b)-w_box,0.1_dp*minval(R_b)),&
                &maxval(R_b)+w_box]
        end if
        if (max_Zvac_plot.gt.min_Zvac_plot) then
            lims_Z = [min_Zvac_plot,max_Zvac_plot]/L_fac
        else
            lims_Z = [minval(Z_b)-w_box,maxval(Z_b)+w_box]
        end if
        allocate(R_pl(n_vac_plot(1)),Z_pl(n_vac_plot(2)),zeta_pl(n_zeta_plot))
        do id = 1,n_vac_plot(1)
            R_pl(id) = lims_R(1) + (lims_R(2)-lims_R(1))*&
                &(id-1._dp)/max(1,n_vac_plot(1)-1)
        end do
        do jd = 1,n_vac_plot(2)
            Z_pl(jd) = lims_Z(1) + (lims_Z(2)-lims_Z(1))*&
                &(jd-1._dp)/max(1,n_vac_plot(2)-1)
        end do
        do kd = 1,n_zeta_plot
            zeta_pl(kd) = pi*(min_zeta_plot + (max_zeta_plot-min_zeta_plot)*&
                &(kd-1._dp)/max(1,n_zeta_plot-1))
        end do
        call writo('R from '//trim(r2strt(L_fac*lims_R(1)))//' to '//&
            &trim(r2strt(L_fac*lims_R(2)))//', Z from '//&
            &trim(r2strt(L_fac*lims_Z(1)))//' to '//&
            &trim(r2strt(L_fac*lims_Z(2)))//', '//trim(i2str(n_zeta_plot))//&
            &' planes')

        ! points: 3-D points in every plane for the field-line style, one (R,Z)
        ! grid for the axisymmetric style
        n_RZ = product(n_vac_plot)
        if (vac%style.eq.1) then                                                ! field-line 3-D
            n_ext = n_RZ*n_zeta_plot
            allocate(x_ext(n_ext,3))
            do kd = 1,n_zeta_plot
                do jd = 1,n_vac_plot(2)
                    do id = 1,n_vac_plot(1)
                        ld = id + (jd-1)*n_vac_plot(1) + (kd-1)*n_RZ
                        x_ext(ld,:) = [R_pl(id)*cos(-zeta_pl(kd)),&
                            &R_pl(id)*sin(-zeta_pl(kd)),Z_pl(jd)]
                    end do
                end do
            end do
        else                                                                    ! axisymmetric
            n_ext = n_RZ
            allocate(x_ext(n_ext,2))
            do jd = 1,n_vac_plot(2)
                do id = 1,n_vac_plot(1)
                    x_ext(id+(jd-1)*n_vac_plot(1),:) = [R_pl(id),Z_pl(jd)]
                end do
            end do
        end if

        ! potential at the points
        allocate(pot(n_ext,n_col))
        ierr = eval_pot_ext(vac,x_ext,Phi_s,D_s,desc_s,pot)
        CHCKERR('')

        ! output
        if (rank_o) then
            ! The points in the plasma: for the axisymmetric style inside the
            ! boundary polygon; for the field-line style, within an interval
            ! length of the nearest boundary point, where the representation
            ! of the constant is not resolved, on the plasma side of its
            ! normal, and farther out where that representation is below -1/2.
            allocate(in_ext(n_ext))
            do ld = 1,n_ext
                if (vac%style.eq.1) then                                        ! field-line 3-D
                    j_min = minloc(sum((vac%x_vec-spread(x_ext(ld,:),1,&
                        &vac%n_bnd))**2,2),1)
                    j_nb = j_min+1                                              ! neighbour on the same line
                    if (mod(j_min,vac%n_ang(1)).eq.0) j_nb = j_min-1
                    if (sum((vac%x_vec(j_min,:)-x_ext(ld,:))**2).lt.&
                        &sum((vac%x_vec(j_nb,:)-vac%x_vec(j_min,:))**2)) then
                        in_ext(ld) = sum((x_ext(ld,:)-vac%x_vec(j_min,:))*&
                            &vac%norm(j_min,:)).gt.0._dp                        ! the stored normal points into the plasma
                    else
                        in_ext(ld) = pot(ld,2*n_X+1).lt.-0.5_dp
                    end if
                else                                                            ! axisymmetric
                    in_ext(ld) = in_polygon(x_ext(ld,1),x_ext(ld,2))
                end if
            end do
            
            ! the potential corrected near the boundary (to first order for
            ! the field-line style), and zero in the plasma
            if (vac%style.eq.1) then                                            ! field-line 3-D
                allocate(grad(vac%n_bnd,3,2*n_X))
                call grad_bnd_1(vac,Phi_ser(:,1:2*n_X),D_ser(:,1:2*n_X),&
                    &abs(max_alpha-min_alpha-2._dp).lt.1.E-8_dp,grad)
                call sub_nearest(vac,x_ext,Phi_ser(:,1:2*n_X),psi_bnd,&
                    &pot(:,2*n_X+1),pot(:,1:2*n_X),grad=grad,&
                    &pot_x=pot(:,2*n_X+2:2*n_X+4))
            else                                                                ! axisymmetric
                call sub_nearest(vac,x_ext,Phi_ser(:,1:2*n_X),psi_bnd,&
                    &pot(:,2*n_X+1),pot(:,1:2*n_X))
            end if
            do pd = 1,2*n_X
                where (in_ext) pot(:,pd) = 0._dp
            end do
            
            ! Cartesian coordinates and the plasma
            allocate(XYZ(n_vac_plot(1),n_vac_plot(2),n_zeta_plot,3))
            allocate(in_pl(n_vac_plot(1),n_vac_plot(2),n_zeta_plot))
            allocate(f_plot(n_vac_plot(1),n_vac_plot(2),n_zeta_plot,2))
            do kd = 1,n_zeta_plot
                do jd = 1,n_vac_plot(2)
                    do id = 1,n_vac_plot(1)
                        XYZ(id,jd,kd,:) = L_fac*[R_pl(id)*cos(-zeta_pl(kd)),&
                            &R_pl(id)*sin(-zeta_pl(kd)),Z_pl(jd)]
                        ld = id + (jd-1)*n_vac_plot(1)
                        if (vac%style.eq.1) ld = ld + (kd-1)*n_RZ
                        in_pl(id,jd,kd) = 0._dp
                        if (in_ext(ld)) in_pl(id,jd,kd) = 1._dp
                    end do
                end do
            end do
            call plot_HDF5('in plasma','vac_in_plasma',in_pl,&
                &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                &descr='1 inside the plasma, 0 in the vacuum')

            ! potential of every eigenmode
            do pd = 1,n_X
                do kd = 1,n_zeta_plot
                    do jd = 1,n_vac_plot(2)
                        do id = 1,n_vac_plot(1)
                            if (vac%style.eq.1) then                            ! field-line 3-D
                                ld = id + (jd-1)*n_vac_plot(1) + (kd-1)*n_RZ
                                phi_c = pot(ld,2*pd-1) + iu*pot(ld,2*pd)
                            else                                                ! axisymmetric
                                ld = id + (jd-1)*n_vac_plot(1)
                                phi_c = (pot(ld,2*pd-1) + iu*pot(ld,2*pd))*&
                                    &exp(iu*vac%prim_X*zeta_pl(kd))
                            end if
                            f_plot(id,jd,kd,:) = [rp(phi_c),ip(phi_c)]
                        end do
                    end do
                end do
                call plot_HDF5('RE Phi',trim(i2str(X_ids(pd)))//&
                    &'_vac_pot_RE',f_plot(:,:,:,1),X=XYZ(:,:,:,1),&
                    &Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                    &descr='real part of the vacuum potential')
                call plot_HDF5('IM Phi',trim(i2str(X_ids(pd)))//&
                    &'_vac_pot_IM',f_plot(:,:,:,2),X=XYZ(:,:,:,1),&
                    &Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                    &descr='imaginary part of the vacuum potential')

                ! boundary points with their normal data and potential
                open(newunit=file_i,file=data_dir//'/'//&
                    &trim(i2str(X_ids(pd)))//'_vac_bnd.dat',status='replace')
                write(file_i,'(A)') '# i_par i_line x(1:dim) '//&
                    &'norm(1:dim) Re(dphi) Im(dphi) Re(phi) Im(phi); '//&
                    &'norm = -J nabla psi, dphi = norm.grad phi, '//&
                    &'normalized units'
                do rd = 1,vac%n_bnd
                    write(file_i,'(2I7,*(ES16.8))') &
                        &mod(rd-1,vac%n_ang(1))+1, (rd-1)/vac%n_ang(1)+1, &
                        &vac%x_vec(rd,:), vac%norm(rd,:), rp(dphi(rd,pd)), &
                        &ip(dphi(rd,pd)), Phi_ser(rd,2*pd-1), Phi_ser(rd,2*pd)
                end do
                close(file_i)

                call writo('Eigenmode '//trim(i2str(X_ids(pd)))//&
                    &': max |phi| on the boundary '//trim(r2strt(&
                    &maxval(abs(Phi_ser(:,2*pd-1:2*pd)))))//', on the grid '//&
                    &trim(r2strt(maxval(abs(f_plot)))))
            end do
        end if

        call lvl_ud(-1)
    contains
        ! Whether (R,Z) lies inside the polygon of the boundary points of the
        ! axisymmetric style (crossing number).
        !> \private
        logical function in_polygon(R,Z) result(res)
            ! input / output
            real(dp), intent(in) :: R, Z                                        ! point

            ! local variables
            integer :: kd, kd_n                                                 ! counter and next point

            res = .false.
            do kd = 1,vac%n_bnd
                kd_n = mod(kd,vac%n_bnd)+1
                if ((vac%x_vec(kd,2).gt.Z) .neqv. (vac%x_vec(kd_n,2).gt.Z)) then
                    if (R.lt.vac%x_vec(kd,1)+(Z-vac%x_vec(kd,2))*&
                        &(vac%x_vec(kd_n,1)-vac%x_vec(kd,1))/&
                        &(vac%x_vec(kd_n,2)-vac%x_vec(kd,2))) res = .not.res
                end if
            end do
        end function in_polygon
    end function vac_pot_plot
    
    !> \public Calculates potential  \c Phi on boundary of BEM  in terms of some
    !! right-hand side .
    !!
    !! This is done by solving the system of equations
    !!  \f$\overline{\text{H}}\overline{\Phi} =
    !!  \overline{\text{G}}\overline{\text{R}}\f$,
    !! where \f$\overline{\text{R}}\f$ is the right-hand side, for example equal
    !! to \f$\overline{\text{E}}\overline{\text{P}}\f$  to solve in  function of
    !! the individual Fourier modes.
    !!
    !! Both vacuum styles assemble  \f$\overline{\text{H}}\f$ in the vacuum-side
    !! convention: the  kernels take  the normal derivative  along the  stored
    !! normal  \f$-\mathcal{J}\nabla\psi\f$, out  of  the vacuum,  and the  diagonal
    !! carries  \f$-2\beta\f$ with  \f$\beta\f$ the  vacuum's internal  angle (the
    !! row-sum identity for the field-line style). Potentials that are harmonic
    !! in the  vacuum and  decay at infinity  then satisfy
    !! \f$\overline{\text{H}}\overline{\phi} =
    !! \overline{\text{G}}\overline{\text{D}\phi}\f$, and potentials that are
    !! harmonic inside the boundary satisfy
    !! \f$\overline{\text{H}}\overline{\phi} - \overline{\text{G}}\overline{\text{D}\phi}
    !! = -4\pi\overline{\phi}\f$ (both verified in the full-stack tests). The
    !! corner  term  also  makes  the  system  nonsingular  when  the  first  and
    !! last points of a closed boundary coincide, as they do for HELENA: the two
    !! rows differ only in the column of that term.
    !!
    !! The linear  system is solved  with STRUMPack when PB3D  was configured
    !! with it, and with ScaLAPACK LU (\c pdgesv) otherwise.
    integer function solve_Phi_BEM(vac,R,Phi,n_RPhi,n_loc_RPhi,lims_c_RPhi,&
        &desc_RPhi) result(ierr)

        use vac_vars, only: in_context, loc_idx
        use num_vars, only: check_vac, rank, n_procs

        character(*), parameter :: rout_name = 'solve_Phi_BEM'

        ! input / output
        type(vac_type), intent(in), target :: vac                               !< vacuum variables
        real(dp), intent(inout), target :: Phi(:,:)                             !< \f$\overline{\Phi}\f$
        real(dp), intent(in) :: R(:,:)                                          !< \f$\overline{\text{R}}\f$
        integer, intent(in) :: n_RPhi(2)                                        !< \c n for \f$\overline{\text{R}}\f$ and \f$\overline{\Phi}\f$
        integer, intent(in) :: n_loc_RPhi(2)                                    !< local \c n for \f$\overline{\text{R}}\f$ and \f$\overline{\Phi}\f$
        integer, intent(in) :: lims_c_RPhi(:,:)                                 !< local limits for row and columns of \f$\overline{\text{R}}\f$ and \f$\overline{\Phi}\f$
        integer, intent(in), target :: desc_RPhi(BLACSCTXTSIZE)                 !< descriptor for \f$\overline{\text{R}}\f$ and \f$\overline{\Phi}\f$

        ! local variables
        integer, target :: desc_GR(BLACSCTXTSIZE)                               ! descriptor for GR
        integer :: info_loc                                                     ! info variable for ScaLAPACK
        integer, allocatable :: ipiv(:)                                         ! pivot variable for the LU
        integer :: iwork_q(1)                                                   ! integer workspace query
        integer, allocatable :: iwork(:)                                        ! integer workspace of the condition estimate
        real(dp) :: anorm, rcond                                                ! 1-norm of the operator and its reciprocal condition number
        real(dp) :: work_q(1)                                                   ! workspace query
        real(dp), allocatable :: work(:)                                        ! workspace of the condition estimate
        real(dp), external :: pdlange                                           ! ScaLAPACK matrix norm
        logical :: solved_ok                                                    ! whether a verified solution was obtained
        real(dp), allocatable, target :: GR(:,:)                                ! GR
        real(dp), allocatable, target :: H_ext(:,:)                             ! copy of H, overwritten by the factorization
        character(len=max_str_ln) :: err_msg                                    ! error message
#ifdef PB3D_WITH_STRUMPACK
        real(dp) :: SDP_err                                                     ! error
        real(dp) :: norms_loc(2)                                                ! squared norms of residual and right-hand side
        real(dp), allocatable :: res_check(:,:)                                 ! residual of the compressed solve
        real(dp), allocatable, target :: H_sdp(:,:)                             ! copy of H_ext for the compressed solver
        type(C_PTR) :: CdescH, CH                                               ! C pointers to descriptor of H_sdp and H_sdp itself
        type(C_PTR) :: CdescGR, CGR                                             ! C pointers to descriptor of GR and GR itself
        type(C_PTR) :: CdescPhi, CPhi                                           ! C pointers to descriptor of Phi and Phi itself
        type(StrumpackDensePackage_F90_double) :: SDP_loc                       ! Strumpack object for local solve H_sdp Phi = GR
#endif
#if ldebug
        integer :: lims_rl(2)                                                   ! vac%lims_r in local coordinates
        integer :: cd                                                           ! global counter for row and col
        integer :: cdl                                                          ! local counter for col
        integer :: i_rd, i_cd                                                   ! index of subrow and subcol
        character(len=max_str_ln) :: plot_title                                 ! plot title
        character(len=max_str_ln) :: plot_name                                  ! plot name
        real(dp), allocatable :: ang(:)                                         ! angle
#endif

        ! initialize ierr
        ierr = 0

        ! Only for processes that are in the blacs context
        if (in_context(vac%ctxt_HG)) then
            ! allocate auxililary matrices
            allocate(GR(vac%n_loc(1),n_loc_RPhi(2)))

            ! initialize descriptors
            call descinit(desc_GR,vac%n_bnd,n_RPhi(2),vac%bs,vac%bs,0,0,&
                &vac%ctxt_HG,max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed for GR')

            ! calculate GR
            call pdgemm('N','N',vac%n_bnd,n_RPhi(2),vac%n_bnd,1._dp,vac%G,1,1,&
                &vac%desc_G,R,1,1,desc_RPhi,0._dp,GR,1,1,desc_GR)

            ! copy the operator: the factorization overwrites it, and the
            ! original stays available for the residual check below
            allocate(H_ext(size(vac%H,1),size(vac%H,2)))
            H_ext = vac%H

            solved_ok = .false.

#ifdef PB3D_WITH_STRUMPACK
            ! try the compressed (HSS) solver first, on a copy of the
            ! operator so that the original stays available for the residual
            ! check and the possible fallback below
            allocate(H_sdp(size(H_ext,1),size(H_ext,2)))
            H_sdp = H_ext

            ! set C pointers
            CH = C_LOC(H_sdp)
            CdescH = C_LOC(vac%desc_H)
            CGR = C_LOC(GR)
            CdescGR = C_LOC(desc_GR)
            CPhi = C_LOC(Phi)
            CdescPhi = C_LOC(desc_RPhi)

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

            ! solve (H + shift I) Phi = G R
            call SDP_F90_double_solve(SDP_loc,CPhi,CdescPhi,CGR,CdescGR)

#if ldebug
            ! accuracy checking
            SDP_err = SDP_F90_double_check_solution(SDP_loc,CH,CdescH,CPhi,&
                &CdescPhi,CGR,CdescGR)
#endif

            ! iterative refinement
            SDP_err = SDP_F90_double_refine(SDP_loc,CH,CdescH,CPhi,CdescPhi,&
                &CGR,CdescGR)
            deallocate(H_sdp)

            ! Verify the solution: the randomized HSS compression can fail
            ! silently when the boundary operator has higher off-diagonal
            ! ranks than the hard-coded sampling parameters can capture
            ! (observed for the field-line 3-D style at moderate sizes). The
            ! residual is checked against the uncompressed operator, and on
            ! failure the solve falls back to ScaLAPACK LU below.
            allocate(res_check(vac%n_loc(1),n_loc_RPhi(2)))
            res_check = GR
            call pdgemm('N','N',vac%n_bnd,n_RPhi(2),vac%n_bnd,1._dp,H_ext,1,1,&
                &vac%desc_H,Phi,1,1,desc_RPhi,-1._dp,res_check,1,1,desc_GR)
            norms_loc(1) = sum(res_check**2)
            norms_loc(2) = sum(GR**2)
            call dgsum2d(vac%ctxt_HG,'all',' ',2,1,norms_loc,2,-1,-1)
            deallocate(res_check)
            if (norms_loc(1).le.1.E-12_dp*norms_loc(2)) then                    ! relative residual <= 1e-6
                solved_ok = .true.
            else
                call writo('The STRUMPACK solution has a relative residual '//&
                    &trim(r2strt(sqrt(norms_loc(1)/max(norms_loc(2),&
                    &tiny(1._dp)))))//' - falling back to ScaLAPACK LU',&
                    &warning=.true.)
            end if
#endif

            if (.not.solved_ok) then
                ! solve (H + shift I) Phi = G R with ScaLAPACK LU (pdgetrf and
                ! pdgetrs, what pdgesv does): the right-hand side is
                ! overwritten with the solution and the matrix destroyed
                ! (which is why H_ext is a copy). With check_vac the
                ! reciprocal condition number in the 1-norm is estimated
                ! between the two (pdgecon) and printed as a CHECKVAC line.
                Phi(1:vac%n_loc(1),1:n_loc_RPhi(2)) = GR                        ! same distribution
                allocate(ipiv(vac%n_loc(1)+vac%bs))
                ipiv = 0
                if (check_vac) then
                    allocate(work(vac%n_loc(2)+vac%bs))
                    anorm = pdlange('1',vac%n_bnd,vac%n_bnd,H_ext,1,1,&
                        &vac%desc_H,work)
                    deallocate(work)
                end if
                call pdgetrf(vac%n_bnd,vac%n_bnd,H_ext,1,1,vac%desc_H,ipiv,&
                    &info_loc)
                if (info_loc.ne.0) then
                    ierr = 1
                    err_msg = 'pdgetrf failed with info '//trim(i2str(info_loc))
                    CHCKERR(err_msg)
                end if
                if (check_vac) then
                    call pdgecon('1',vac%n_bnd,H_ext,1,1,vac%desc_H,anorm,&
                        &rcond,work_q,-1,iwork_q,-1,info_loc)
                    allocate(work(int(work_q(1))),iwork(iwork_q(1)))
                    call pdgecon('1',vac%n_bnd,H_ext,1,1,vac%desc_H,anorm,&
                        &rcond,work,size(work),iwork,size(iwork),info_loc)
                    deallocate(work,iwork)
                    call writo('CHECKVAC condition n_bnd='//&
                        &trim(i2str(vac%n_bnd))//' norm1='//trim(r2strt(anorm))//&
                        &' cond1='//trim(r2strt(1._dp/max(rcond,tiny(1._dp)))),&
                        &persistent=rank.eq.n_procs-1)
                end if
                call pdgetrs('N',vac%n_bnd,n_RPhi(2),H_ext,1,1,vac%desc_H,ipiv,&
                    &Phi,1,1,desc_RPhi,info_loc)
                if (info_loc.ne.0) then
                    ierr = 1
                    err_msg = 'pdgetrs failed with info '//trim(i2str(info_loc))
                    CHCKERR(err_msg)
                end if
                deallocate(ipiv)
            end if

#if ldebug
#ifdef PB3D_WITH_STRUMPACK
            ! statistics
            call SDP_F90_double_print_statistics(SDP_loc)
#endif

            ! user output
            ! set global angle (the same for every subrow and column; the
            ! subrow slices below use global indices, so ang has to span
            ! the whole boundary, also when the rows are distributed over
            ! multiple block-cyclic subrows)
            allocate(ang(vac%n_bnd))
            ang = reshape(vac%ang,[vac%n_bnd])

            subcols2: do i_cd = 1,size(lims_c_RPhi,2)
                col: do cd = lims_c_RPhi(1,i_cd),lims_c_RPhi(2,i_cd)
                    ! set local column index
                    cdl = loc_idx(lims_c_RPhi,i_cd,cd)

                    subrows2: do i_rd = 1,size(vac%lims_r,2)
                        ! set local row limits
                        lims_rl = [loc_idx(vac%lims_r,i_rd,vac%lims_r(1,i_rd)), &
                        &loc_idx(vac%lims_r,i_rd,vac%lims_r(2,i_rd))]

                        ! output
                        plot_title = 'for column '//trim(i2str(cd))//&
                            &' and starting row '//&
                            &trim(i2str(vac%lims_r(1,i_rd)))
                        plot_name = '_'//trim(i2str(cd))//'_'//&
                            &trim(i2str(vac%lims_r(1,i_rd)))
                        call print_ex_2D('Phi '//trim(plot_title),&
                            &'Phi'//trim(plot_name),&
                            &Phi(lims_rl(1):lims_rl(2),cdl),&
                            &x=ang(vac%lims_r(1,i_rd):vac%lims_r(2,i_rd)),&
                            &draw=.false.)
                        call draw_ex(['Phi '//trim(plot_title)],&
                            &'Phi'//trim(plot_name),1,1,.false.)
                    end do subrows2
                end do col
            end do subcols2

            ! clean up
            deallocate(ang)
#endif

#ifdef PB3D_WITH_STRUMPACK
            ! destroy Strumpack object
            call SDP_F90_double_destroy(SDP_loc)
#endif
        end if
    end function solve_Phi_BEM
    
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
        use grid_vars, only: min_par_X, max_par_X
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
        
        ! misc_vac (the last two, the limits of the parallel angle in units of
        ! pi, let POST rebuild the parallel angle of the field-line style,
        ! which the vacuum potential plot needs; reconstruct_PB3D_vac reads
        ! them only when present)
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
            &vac%n_ang*1._dp,vac%jq,min_par_X,max_par_X]
        
        ! sec_X
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'sec_X'
        allocate(vac_1D_loc%tot_i_min(1),vac_1D_loc%tot_i_max(1))
        allocate(vac_1D_loc%loc_i_min(1),vac_1D_loc%loc_i_max(1))
        vac_1D_loc%tot_i_min = [1]
        vac_1D_loc%tot_i_max = n_mod_X
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(n_mod_X))
        vac_1D_loc%p = vac%sec_X*1._dp
        
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
        
        ! copy variables specific to style
        select case (vac%style)
            case (1)                                                            ! field-line 3-D
                ! h_fac
                vac_1D_loc => vac_1D(id); id = id+1
                vac_1D_loc%var_name = 'h_fac'
                allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
                allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
                vac_1D_loc%tot_i_min = [1,1]
                vac_1D_loc%tot_i_max = shape(vac%h_fac)
                vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
                vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
                allocate(vac_1D_loc%p(size(vac%h_fac)))
                vac_1D_loc%p = reshape(vac%h_fac,[size(vac%h_fac)])
            case (2)                                                            ! axisymmetric
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
        end select
        
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
