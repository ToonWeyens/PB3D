!------------------------------------------------------------------------------!
!   Calculates the equilibrium quantities, making use of the metric_ops,       !
!   eq_vars, etc                                                               !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use num_vars, only: pi, dp, max_str_ln
    use message_ops, only: print_ar_2, lvl_ud, writo
    use output_ops, only: print_GP_3D, draw_GP_animated, draw_GP, &
        &print_HDF5_3D, print_GP_2D
    use str_ops, only: i2str
    
    implicit none
    private
    public calc_eq, read_eq

contains
    ! calculate  the equilibrium  quantities on  a grid  determined by  straight
    ! field lines. This  grid has the dimensions  (n_par,grp_n_r_eq) where n_par
    ! is  the  number  of  points  taken along  the  magnetic  field  lines  and
    ! grp_n_r_eq .le.  n_r_eq is the  normal extent  in the equilibrium  grid of
    ! this rank. It is determined so  that the perturbation quantities that will
    ! be needed in this rank are fully covered, so no communication is necessary
    integer function calc_eq(alpha) result(ierr)
        use eq_vars, only: calc_mesh, calc_flux_q, dealloc_eq, calc_RZL, &
            &normalize_eq_vars, prepare_RZL, check_and_limit_mesh, init_eq, &
            &adapt_HEL_to_eq, &
            &q_saf_V, q_saf_FD, rot_t_V, rot_t_FD, flux_p_V, theta_V, &
            &flux_p_FD,flux_t_V, flux_t_FD, pres_V, pres_FD, n_par, zeta_V, &
            &theta_H, zeta_H, q_saf_H, rot_t_H, pres_H, flux_p_H, flux_t_H
        use metric_ops, only: calc_g_C, calc_g_C, calc_T_VC, calc_g_V, &
            &init_metric, calc_T_VF, calc_inv_met, calc_g_F, calc_jac_C, &
            &calc_jac_V, calc_jac_F, calc_f_deriv, dealloc_metric, &
            &normalize_metric_vars, calc_jac_H, calc_T_HF, calc_h_H, &
            &T_VF, T_FV, g_F, h_F, det_T_VF, det_T_FV, jac_F, g_FD, h_FD, &
            &jac_FD, T_HF, T_FH, det_T_HF, det_T_FH, g_H, h_H
        use utilities, only: derivs
        use num_vars, only: max_deriv, ltest, use_pol_flux, plot_grid, &
            &eq_style, grp_rank, use_normalization
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), pointer :: T_FE(:,:,:,:,:,:,:)                                ! E(quilibrium) (VMEC or HELENA) to F(lux) (lower)
        real(dp), pointer :: q_saf_E(:,:)                                       ! safety factor in E(equilibrium) coordinates
        real(dp), pointer :: rot_t_E(:,:)                                       ! rot. transform in E(equilibrium) coordinates
        real(dp), pointer :: flux_p_E(:,:), flux_t_E(:,:), pres_E(:,:)          ! pol. flux, tor. flux and pressure, and norm. Deriv. in E(quilibrium) coords.
        
        ! initialize ierr
        ierr = 0
        
        call writo('Start setting up equilibrium quantities in '//&
            &trim(i2str(n_par))//' discrete parallel points')
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(1)
            call writo('Start initializing variables')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! initialize equilibrium quantities
            call writo('Initialize equilibrium quantities...')
            ierr = init_eq()
            CHCKERR('')
            
            ! initialize metric quantities
            call writo('Initialize metric quantities...')
            ierr = init_metric()
            CHCKERR('')
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Variables initialized')
            
            call writo('Start determining the equilibrium grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! calculate  angular mesh  points (theta_V,zeta_V)  that follow  the
            ! magnetic field line determined by alpha
            ! Note: The normal mesh is determined by VMEC, it can either use the
            ! toroidal flux or the  poloidal flux (VMEC_use_pol_flux). This does
            ! NOT have to coincide with  use_pol_flux used by PB3D. However, the
            ! mesh is tabulated in the VMEC normal coordinate
            ierr = calc_mesh(alpha)
            CHCKERR('')
            
            ! adapt the HELENA variables to the equilibrium mesh
            if (eq_style.eq.2) then
                ierr = adapt_HEL_to_eq()
                CHCKERR('')
            end if
            
            ! plot grid if requested
            if (plot_grid .and. grp_rank.eq.0) then                             ! only group masters
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ierr = plot_grid_real(theta_V,zeta_V,0._dp,1._dp)
                        CHCKERR('')
                    case (2)                                                    ! HELENA
                        ierr = plot_grid_real(theta_H,zeta_H,0._dp,1._dp)
                        CHCKERR('')
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            else
                call writo('Grid plot not requested')
            end if
            
            ! check whether the mesh has been calculated correctly
            ! if so, limit  the normal range of theta_V and  zeta_V to the local
            ! range
            ierr = check_and_limit_mesh(alpha)
            CHCKERR('')
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium grid determined')
            
            call writo('Calculating equilibrium quantities on equilibrium grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! calculate the  cylindrical variables  R, Z and  lambda and
                    ! derivatives
                    call writo('Calculate R,Z,L...')
                    ierr = prepare_RZL()
                    CHCKERR('')
                    do id = 0,2
                        ierr = calc_RZL(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate flux quantities
                    call writo('Calculate flux quantities...')
                    ierr = calc_flux_q()
                    CHCKERR('')
                    
                    ! calculate the metrics in the cylindrical coordinate system
                    call writo('Calculate g_C...')                              ! h_C is not necessary
                    do id = 0,1
                        ierr = calc_g_C(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate  the  jacobian  in  the  cylindrical  coordinate
                    ! system
                    call writo('Calculate jac_C...')
                    do id = 0,1
                        ierr = calc_jac_C(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate   the  transformation  matrix  C(ylindrical)  ->
                    ! V(MEC)
                    call writo('Calculate T_VC...')
                    do id = 0,1
                        ierr = calc_T_VC(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the metric factors in the VMEC coordinate system
                    call writo('Calculate g_V...')
                    do id = 0,1
                        ierr = calc_g_V(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the jacobian in the VMEC coordinate system
                    call writo('Calculate jac_V...')
                    do id = 0,1
                        ierr = calc_jac_V(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the transformation matrix V(MEC) -> F(lux)
                    call writo('Calculate T_VF...')
                    do id = 0,1
                        ierr = calc_T_VF(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the inverse of the transformation matrix T_VF
                    call writo('Calculate T_FV...')
                    do id = 0,1
                        ierr = calc_inv_met(T_FV,T_VF,derivs(id))
                        CHCKERR('')
                        ierr = calc_inv_met(det_T_FV,det_T_VF,derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! point T_FE, q_saf_E, rot_t_E, pres_E, flux_p_E, flux_t_E
                    T_FE => T_FV
                    pres_E => pres_V
                    flux_p_E => flux_p_V
                    flux_t_E => flux_t_V
                    if (use_pol_flux) then
                        q_saf_E => q_saf_V
                    else
                        rot_t_E => rot_t_V
                    end if
                case (2)                                                        ! HELENA
                    ! calculate flux quantities
                    ierr = calc_flux_q()
                    CHCKERR('')
                    
                    ! calculate the jacobian in the HELENA coordinate system
                    call writo('Calculate jac_H...')
                    do id = 0,1
                        ierr = calc_jac_H(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate  the  metric  factors in  the HELENA  coordinate
                    ! system
                    call writo('Calculate h_H...')
                    do id = 0,1
                        ierr = calc_h_H(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the inverse g_H of the metric factors h_H
                    call writo('Calculate g_H...')
                    do id = 0,1
                        ierr = calc_inv_met(g_H,h_H,derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the transformation matrix H(ELENA) -> F(lux)
                    call writo('Calculate T_HF...')
                    do id = 0,1
                        ierr = calc_T_HF(derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! calculate the inverse of the transformation matrix T_VH
                    call writo('Calculate T_FH...')
                    do id = 0,1
                        ierr = calc_inv_met(T_FH,T_HF,derivs(id))
                        CHCKERR('')
                        ierr = calc_inv_met(det_T_FH,det_T_HF,derivs(id))
                        CHCKERR('')
                    end do
                    
                    ! point T_FE, q_saf_E, rot_t_E, pres_E, flux_p_E, flux_t_E
                    T_FE => T_FH
                    pres_E => pres_H
                    flux_p_E => flux_p_H
                    flux_t_E => flux_t_H
                    if (use_pol_flux) then
                        q_saf_E => q_saf_H
                    else
                        rot_t_E => rot_t_H
                    end if
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! calculate the metric factors in the Flux coordinate system
            call writo('Calculate g_F...')
            do id = 0,1
                ierr = calc_g_F(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the inverse h_F of the metric factors g_F
            call writo('Calculate h_F...')
            do id = 0,1
                ierr = calc_inv_met(h_F,g_F,derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the jacobian in the Flux coordinate system
            call writo('Calculate jac_F...')
            do id = 0,1
                ierr = calc_jac_F(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the derivatives in Flux coordinates from the derivatives
            ! in VMEC coordinates
            call writo('Calculate derivatives in flux coordinates...')
            do id = 0,1
                ! g_FD
                ierr = calc_f_deriv(g_F,T_FE,g_FD,max_deriv-[1,1,1],&
                    &derivs(id))
                CHCKERR('')
                
                ! h_FD
                ierr = calc_f_deriv(h_F,T_FE,h_FD,max_deriv-[1,1,1],&
                    &derivs(id))
                CHCKERR('')
                
                ! jac_FD
                ierr = calc_f_deriv(jac_F,T_FE,jac_FD,max_deriv-[1,1,1],&
                    &derivs(id))
                CHCKERR('')
                
                ! pres_FD
                ierr = calc_f_deriv(pres_E,T_FE(1,:,2,1,:,0,0),pres_FD(:,id),&
                    &max_deriv(1)-1,id)
                CHCKERR('')
                    
                ! flux_p_FD
                ierr = calc_f_deriv(flux_p_E,T_FE(1,:,2,1,:,0,0),&
                    &flux_p_FD(:,id),max_deriv(1)-1,id)
                CHCKERR('')
                    
                ! flux_t_FD
                ierr = calc_f_deriv(flux_t_E,T_FE(1,:,2,1,:,0,0),&
                    &flux_t_FD(:,id),max_deriv(1)-1,id)
                CHCKERR('')
                if (eq_style.eq.1) flux_t_FD = - flux_t_FD                      ! conversion VMEC LH -> RH coord. system
                
                if (use_pol_flux) then
                    ! q_saf_FD
                    ierr = calc_f_deriv(q_saf_E,T_FE(1,:,2,1,:,0,0),&
                        &q_saf_FD(:,id),max_deriv(1)-1,id)
                    CHCKERR('')
                    if (eq_style.eq.1) q_saf_FD = - q_saf_FD                    ! conversion VMEC LH -> RH coord. system
                else
                    ! rot_t_FD
                    ierr = calc_f_deriv(rot_t_E,T_FE(1,:,2,1,:,0,0),&
                        &rot_t_FD(:,id),max_deriv(1)-1,id)
                    CHCKERR('')
                    if (eq_style.eq.1) rot_t_FD = - rot_t_FD                    ! conversion VMEC LH -> RH coord. system
                end if
            end do
            
            ! normalize the quantities
            if (use_normalization) then
                call writo('Normalize the equilibrium and metric quantities...')
                call normalize_eq_vars
                call normalize_metric_vars
            end if
            
            ! deallocate unused equilibrium quantities
            if (.not.ltest) then
                call writo('Deallocate unused equilibrium and metric &
                    &quantities...')
                ierr = dealloc_eq()
                CHCKERR('')
                ierr = dealloc_metric()
                CHCKERR('')
            end if
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium quantities calculated on equilibrium grid')
            
        call lvl_ud(-1)
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call writo('Done setting up equilibrium quantities')
    end function calc_eq
    
    ! reads the equilibrium input file
    integer function read_eq() result(ierr)
        use num_vars, only: eq_style, glb_rank
        use VMEC_vars, only: read_VMEC
        use HEL_vars, only: read_HEL
        use eq_vars, only: n_r_eq, eq_use_pol_flux
        
        character(*), parameter :: rout_name = 'read_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only do this for the group master
        if (glb_rank.eq.0) then
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = read_VMEC(n_r_eq,eq_use_pol_flux)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ierr = read_HEL(n_r_eq,eq_use_pol_flux)
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
    end function read_eq
    
    ! plots the grid in real space
    ! The  variables theta  and  zeta follow  the magnetic  field  lines in  the
    ! equilibrium grid  and are tabulated along  the magnetic field lines  for a
    ! series  of flux  surfaces. The  start  and end  index of  the normal  flux
    ! surfaces in the equilibrium grid has  to be provided with min_r and max_r.
    ! Ideally, this routine should be run  with the full equilibrium grid normal
    ! extent, which implies for example, min_r = 0, max_r = 1.
    ! [MPI] Collective call
    integer function plot_grid_real(theta,zeta,min_r,max_r) &
        &result(ierr)
        use eq_vars, only: n_par, calc_XYZ_grid
        use num_vars, only: output_style, alpha_job_nr
        
        character(*), parameter :: rout_name = 'plot_grid_real'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:)                                      ! theta in VMEC or HELENA coordinates
        real(dp), intent(in) :: zeta(:,:)                                       ! zeta in VMEC or HELENA coordinates
        real(dp), intent(in) :: min_r, max_r                                    ! minimum and maximum normal range
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=*), parameter :: anim_name = &
            &'Magnetic field in flux surfaces'                                  ! name of animation
        real(dp), allocatable :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)             ! X, Y and Z of surface in Axisymmetric coordinates
        real(dp), allocatable :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)             ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
        real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)            ! theta and zeta for plots
        real(dp), allocatable :: r_plot(:)                                      ! r for plots
        integer :: id                                                           ! counter
        integer :: loc_n_par, loc_n_r                                           ! local nr. of parallel and normal points
        integer :: n_theta_plot = 40                                            ! nr. of poloidal points in plot
        integer :: n_zeta_plot = 160                                            ! nr. of toroidal points in plot
        
        ! initialize ierr
        ierr = 0
        
        call writo('Plotting magnetic field and flux surfaces')
        
        call lvl_ud(1)
        
        ! some tests
        if (size(theta,1).ne.size(zeta,1) .or. &
            &size(theta,1).ne.n_par .or. &
            &size(theta,2).ne.size(zeta,2)) then
            err_msg = 'theta and zeta have to have the correct sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! 1. plot flux surfaces
        call writo('writing flux surfaces')
        
        ! initialize loc_n_par and loc_n_r
        loc_n_par = size(theta,1)
        loc_n_r = size(theta,2)
        
        ! initialize theta_plot, zeta_plot and r_plot
        allocate(theta_plot(n_theta_plot,n_zeta_plot,loc_n_r))
        do id = 1,n_theta_plot
            theta_plot(id,:,:) = pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)      ! better to start from pi for the plot
        end do
        allocate(zeta_plot(n_theta_plot,n_zeta_plot,loc_n_r))
        do id = 1,n_zeta_plot
            zeta_plot(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
        end do
        allocate(r_plot(loc_n_r))
        do id = 1,loc_n_r
            r_plot(id) = min_r + (id-1._dp)/(loc_n_r-1)*(max_r-min_r)
        end do
        
        ! calculate X,Y and Z
        ierr = calc_XYZ_grid(theta_plot,zeta_plot,r_plot,X_1,Y_1,Z_1)
        CHCKERR('')
        
        ! deallocate
        deallocate(theta_plot,zeta_plot)
        
        ! 2. plot field lines
        call writo('writing field lines')
        
        ! initialize theta_plot and zeta_plot
        allocate(theta_plot(1,loc_n_par,loc_n_r),zeta_plot(1,loc_n_par,loc_n_r))
        do id = 1,loc_n_par
            theta_plot(1,id,:) = theta(id,:)
            zeta_plot(1,id,:) = zeta(id,:)
        end do
        
        ! calculate X,Y and Z
        ierr = calc_XYZ_grid(theta_plot,zeta_plot,r_plot,X_2,Y_2,Z_2)
        CHCKERR('')
        
        ! deallocate
        deallocate(theta_plot,zeta_plot)
        
        ! delegate to child routines
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                call plot_grid_real_GP(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)
            case(2)                                                             ! HDF5 output
                ierr = plot_grid_real_HDF5(x_1,x_2,y_1,y_2,z_1,z_2,anim_name)
                CHCKERR('')
            case default
                err_msg = 'No style associated with '//trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! deallocate
        deallocate(X_1,Y_1,Z_1)
        deallocate(X_2,Y_2,Z_2)
        
        call lvl_ud(-1)
        
        call writo('Done plotting magnetic field and flux surfaces')
    contains
        subroutine plot_grid_real_GP(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)         ! GNUPlot version
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: file_name(2), plot_title(2)            ! file name and title
            character(len=max_str_ln) :: draw_ops(2)                            ! individual plot command
            integer :: n_r                                                      ! total number of normal points
           
            ! set up n_r
            n_r = size(X_1,3)                                                   ! should be same for all other X_i, Y_i and Z_i
            
            ! user output
            call writo('Drawing animation with GNUPlot')
            call lvl_ud(1)
            
            ! set names
            plot_title(1) = 'Magnetic Flux Surface for alpha job '//&
                &trim(i2str(alpha_job_nr))
            file_name(1) = 'Flux_surfaces_'//trim(i2str(alpha_job_nr))//'.dat'
            plot_title(2) = 'Magnetic Field Line for alpha job '//&
                &trim(i2str(alpha_job_nr))
            file_name(2) = 'B_field_'//trim(i2str(alpha_job_nr))//'.dat'
            
            ! write flux surfaces
            call print_GP_3D(trim(plot_title(1)),trim(file_name(1)),&
                &Z_1(:,:,1:n_r-1),x=X_1(:,:,1:n_r-1),&
                &y=Y_1(:,:,1:n_r-1),draw=.false.)                               ! plot all but last points
            
            ! initialize n_r
            n_r = size(X_2,3)
            
            ! write magnetic field lines
            call print_GP_3D(trim(plot_title(2)),trim(file_name(2)),&
                &Z_2(:,:,2:n_r),x=X_2(:,:,2:n_r),y=Y_2(:,:,2:n_r),&
                &draw=.false.)                                                  ! plot all bit first points
            
            ! plot both output files
            call writo('creating GNUPlot animation')
            draw_ops(1) = 'linecolor rgb ''#d3d3d3'' linewidth 1'
            draw_ops(2) = 'linecolor rgb ''black'' linewidth 3'
            call draw_GP_animated(anim_name,file_name,n_r-1,.false.,&
                &delay=50,draw_ops=draw_ops)
            
            call lvl_ud(-1)
        end subroutine plot_grid_real_GP
        integer function plot_grid_real_HDF5(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)&
            & result(ierr)                                                      ! HDF5 version
            use HDF5_vars, only: open_HDF5_file, add_HDF5_item, &
                &print_HDF5_top, print_HDF5_geom, print_HDF5_3D_data_item, &
                &print_HDF5_grid, close_HDF5_file, &
                &XML_str_type, HDF5_file_type
            
            character(*), parameter :: rout_name = 'plot_grid_real_HDF5'
            
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: plot_title(2)                          ! plot title for flux surface and for field lines
            integer :: id                                                       ! counter
            type(HDF5_file_type) :: file_info                                   ! file info
            type(XML_str_type) :: grids(2)                                      ! grid in spatial collection (flux surface, field line)
            type(XML_str_type), allocatable :: space_col_grids(:)               ! grid with space collection at different times
            type(XML_str_type) :: time_col_grid                                 ! grid with time collection
            type(XML_str_type) :: top                                           ! topology
            type(XML_str_type) :: XYZ(3)                                        ! data items for geometry
            type(XML_str_type) :: geom                                          ! geometry
            integer :: loc_dim(3,2)                                             ! local dimensions (flux surfaces, field lines)
            integer :: n_r                                                      ! nr. of normal points
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Drawing animation with HDF5')
            call lvl_ud(1)
            
            ! set up loc_dim and n_r
            loc_dim(:,1) = [size(X_1,1),size(X_1,2),1]
            loc_dim(:,2) = [size(X_2,1),size(X_2,2),1]
            n_r = size(X_1,3)                                                   ! should be same for all other X_i, Y_i and Z_i
            
            ! set up plot titles
            plot_title(1) = 'Magnetic Flux Surface for alpha job '//&
                &trim(i2str(alpha_job_nr))
            plot_title(2) = 'Magnetic Field Line for alpha job '//&
                &trim(i2str(alpha_job_nr))
            
            ! open HDF5 file
            ierr = open_HDF5_file(file_info,'field_lines_in_flux_surfaces_'//&
                &trim(i2str(alpha_job_nr)),description=anim_name,&
                &ind_plot=.true.)
            CHCKERR('')
            
            ! create grid for time collection
            allocate(space_col_grids(n_r-1))
            
            ! loop over all normal points
            do id = 1,n_r-1
                ! A. start with flux surface
                
                ! print topology
                call print_HDF5_top(top,2,loc_dim(:,1))
                
                ! print data item for X
                ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
                    &'X_surf_'//trim(i2str(id)),X_1(:,:,id:id),loc_dim(:,1))
                CHCKERR('')
                
                ! print data item for Y
                ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Y_surf_'//trim(i2str(id)),Y_1(:,:,id:id),loc_dim(:,1))
                CHCKERR('')
                
                ! print data item for Z
                ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
                    &'Z_surf_'//trim(i2str(id)),Z_1(:,:,id:id),loc_dim(:,1))
                CHCKERR('')
                
                ! print geometry with X, Y and Z data item
                call print_HDF5_geom(geom,2,XYZ,.true.)
                
                ! create a grid with the topology and the geometry
                ierr = print_HDF5_grid(grids(1),plot_title(1),1,&
                    &grid_top=top,grid_geom=geom,reset=.true.)
                CHCKERR('')
                
                ! B. end with magnetic field
                
                ! print topology
                call print_HDF5_top(top,2,loc_dim(:,2))
                
                ! print data item for X
                ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
                    &'X_field_'//trim(i2str(id)),X_2(:,:,id+1:id+1),loc_dim(:,2))
                CHCKERR('')
                
                ! print data item for Y
                ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Y_field_'//trim(i2str(id)),Y_2(:,:,id+1:id+1),loc_dim(:,2))
                CHCKERR('')
                
                ! print data item for Z
                ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
                    &'Z_field_'//trim(i2str(id)),Z_2(:,:,id+1:id+1),loc_dim(:,2))
                CHCKERR('')
                
                ! print geometry with X, Y and Z data item
                call print_HDF5_geom(geom,2,XYZ,.true.)
                
                ! create a grid with the topology and the geometry
                ierr = print_HDF5_grid(grids(2),plot_title(2),1,&
                    &grid_top=top,grid_geom=geom,reset=.true.)
                CHCKERR('')
                
                ! C. merge flux surface and magnetic field in to spatial collection
                ierr = print_HDF5_grid(space_col_grids(id),'spatial collection',&
                    &3,grid_time=id*1._dp,grid_grids=grids,reset=.true.)
                CHCKERR('')
            end do
            
            ! create grid collection from individual grids and reset them
            ierr = print_HDF5_grid(time_col_grid,'time collection',2,&
                &grid_grids=space_col_grids,reset=.true.)
            CHCKERR('')
            
            ! add collection grid to HDF5 file and reset it
            call add_HDF5_item(file_info,time_col_grid,reset=.true.)
            
            ! close HDF5 file
            ierr = close_HDF5_file(file_info)
            CHCKERR('')
            
            call lvl_ud(-1)
        end function plot_grid_real_HDF5
    end function plot_grid_real
end module eq_ops
