!------------------------------------------------------------------------------!
!   Calculates the equilibrium quantities, making use of the metric_ops,       !
!   eq_vars, etc                                                               !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use num_vars, only: pi, dp, max_str_ln
    use message_ops, only: print_ar_2, lvl_ud, writo
    use output_ops, only: print_GP_3D, draw_GP_animated, draw_GP
    use str_ops, only: i2str
    
    implicit none
    private
    public calc_eq

contains
    ! calculate the equilibrium quantities on a grid determined by straight field
    ! lines.
    integer function calc_eq(alpha) result(ierr)
        use eq_vars, only: calc_mesh, calc_flux_q, dealloc_eq, calc_RZL, &
            &normalize_eq_vars, prepare_RZL, check_and_limit_mesh, init_eq, &
            &q_saf_V, q_saf_FD, rot_t_V, rot_t_FD, flux_p_V, theta_V, &
            &flux_p_FD,flux_t_V, flux_t_FD, pres, pres_FD, n_par, zeta_V
        use metric_ops, only: calc_g_C, calc_g_C, calc_T_VC, calc_g_V, &
            &init_metric, calc_T_VF, calc_inv_met, calc_g_F, calc_jac_C, &
            &calc_jac_V, calc_jac_F, calc_f_deriv, dealloc_metric, &
            &normalize_metric_vars, &
            &T_VF, T_FV, g_F, h_F, det_T_VF, det_T_FV, jac_F, g_FD, h_FD, &
            &jac_FD
        use utilities, only: derivs
        use VMEC_vars, only: n_r_eq
        use num_vars, only: max_deriv, ltest, use_pol_flux, plot_grid
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        integer :: id
        
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
            call init_eq
            
            ! initialize metric quantities
            call writo('Initialize metric quantities...')
            call init_metric
            
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
            ! NOT have to coincide with use_pol_flux used by PB3D
            ierr = calc_mesh(alpha)
            CHCKERR('')
            
            ! plot grid if requested
            if (plot_grid) then
                ierr = plot_grid_real(theta_V,zeta_V,1,n_r_eq)
                CHCKERR('')
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
            
            ! calculate  the   cylindrical  variables   R,  Z  and   lambda  and
            ! derivatives
            call writo('Calculate R,Z,L...')
            ierr = prepare_RZL()
            CHCKERR('')
            do id = 0,2
                ierr = calc_RZL(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate flux quantities
            ierr = calc_flux_q()
            CHCKERR('')
            
            ! calculate the metrics in the cylindrical coordinate system
            call writo('Calculate g_C...')                                      ! h_C is not necessary
            do id = 0,1
                ierr = calc_g_C(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the jacobian in the cylindrical coordinate system
            call writo('Calculate jac_C...')
            do id = 0,1
                ierr = calc_jac_C(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
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
            
            ! calculate the transformation matrix V(mec) -> F(lux)
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
                ierr = calc_f_deriv(g_F,T_FV,g_FD,max_deriv-[1,1,1],&
                    &derivs(id))
                CHCKERR('')
                
                ! h_FD
                ierr = calc_f_deriv(h_F,T_FV,h_FD,max_deriv-[1,1,1],&
                    &derivs(id))
                CHCKERR('')
                
                ! jac_FD
                ierr = calc_f_deriv(jac_F,T_FV,jac_FD,max_deriv-[1,1,1],&
                    &derivs(id))
                CHCKERR('')
                
                ! pres_FD
                ierr = calc_f_deriv(pres,T_FV(1,:,2,1,:,0,0),pres_FD(:,id),&
                    &max_deriv(1)-1,id)
                CHCKERR('')
                    
                ! flux_p_FD
                ierr = calc_f_deriv(flux_p_V,T_FV(1,:,2,1,:,0,0),&
                    &flux_p_FD(:,id),max_deriv(1)-1,id)
                CHCKERR('')
                    
                ! flux_t_FD
                ierr = calc_f_deriv(flux_t_V,T_FV(1,:,2,1,:,0,0),&
                    &flux_t_FD(:,id),max_deriv(1)-1,id)
                CHCKERR('')
                flux_t_FD = - flux_t_FD                                         ! conversion LH -> RH coord. system
                
                if (use_pol_flux) then
                    ! q_saf_FD
                    ierr = calc_f_deriv(q_saf_V,T_FV(1,:,2,1,:,0,0),&
                        &q_saf_FD(:,id),max_deriv(1)-1,id)
                    CHCKERR('')
                    q_saf_FD = - q_saf_FD                                       ! conversion LH -> RH coord. system
                else
                    ! rot_t_FD
                    ierr = calc_f_deriv(rot_t_V,T_FV(1,:,2,1,:,0,0),&
                        &rot_t_FD(:,id),max_deriv(1)-1,id)
                    CHCKERR('')
                    rot_t_FD = - rot_t_FD                                       ! conversion LH -> RH coord. system
                end if
            end do
            
            ! normalize the quantities
            call writo('Normalize the equilibrium and metric quantities...')
            call normalize_eq_vars
            call normalize_metric_vars
            
            ! deallocate unused equilibrium quantities
            if (.not.ltest) then
                call writo('Deallocate unused equilibrium and metric quantities...')
                call dealloc_eq
                call dealloc_metric
            end if
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium quantities calculated on equilibrium grid')
        
        call lvl_ud(-1)
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call writo('Done setting up equilibrium quantities')
    end function calc_eq
    
    ! plots the grid in real space
    ! [MPI] Collective call
    integer function plot_grid_real(theta,zeta,min_r,max_r) result(ierr)
        use fourier_ops, only: calc_trigon_factors, fourier2real
        use VMEC_vars, only: R_c, R_s, Z_c, Z_s
        use eq_vars, only: n_par
        
        character(*), parameter :: rout_name = 'plot_grid_real'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:)                                      ! theta in VMEC coordinates
        real(dp), intent(in) :: zeta(:,:)                                       ! zeta in VMEC coordinates
        integer, intent(in) :: min_r, max_r                                     ! minimum and maximum normal index in VMEC output
        
        ! local variables
        real(dp), allocatable :: R(:,:,:), Z(:,:,:)                             ! R and Z in Cylindrical coordinates
        real(dp), allocatable :: X(:,:,:), Y(:,:,:)                             ! X and Y in Axisymmetric coordinates
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: loc_n_par, loc_n_r                                           ! local nr. of parallel and normal points
        real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      ! trigonometric factor cosine for the inverse fourier transf.
        character(len=max_str_ln) :: file_name(2), plot_title(2)                ! file name and title
        integer :: id                                                           ! counter
        real(dp), allocatable :: theta_surf(:,:), zeta_surf(:,:)                ! theta and zeta for flux surface plot
        integer :: n_theta_plot = 20                                            ! nr. of poloidal points in plot
        integer :: n_zeta_plot = 80                                             ! toroidal points in plot
        character(len=max_str_ln) :: draw_ops(2)                ! individual plot command
        
        ! initialize ierr
        ierr = 0
        
        call writo('Plotting magnetic field and flux surfaces')
        
        call lvl_ud(1)
        
        ! some tests
        if (size(theta,1).ne.size(zeta,1) .or. &
            &size(theta,1).ne.n_par .or. &
            &size(theta,2).ne.size(zeta,2) .or. &
            &size(theta,2).ne.max_r-min_r+1) then
            err_msg = 'theta and zeta have to have the correct sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! 1. plot flux surfaces
        
        ! initialize loc_n_r
        loc_n_r = size(theta,2)
        
        ! initialize theta_surf and zeta_surf
        allocate(theta_surf(n_zeta_plot,loc_n_r),&
            &zeta_surf(n_zeta_plot,loc_n_r))
        
        ! initialize R and Z
        allocate(R(n_theta_plot,n_zeta_plot,loc_n_r),&
            &Z(n_theta_plot,n_zeta_plot,loc_n_r))
        allocate(X(n_theta_plot,n_zeta_plot,loc_n_r),&
            &Y(n_theta_plot,n_zeta_plot,loc_n_r))
        
        ! inverse fourier transform with trigonometric factors
        ! set up zeta_surf (same for each poloidal point)
        do id = 1,n_zeta_plot
            zeta_surf(id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
        end do
        do id = 1,n_theta_plot
            ! set up theta_surf (indicates poloidal point)
            theta_surf = (id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)
            
            ierr = calc_trigon_factors(theta_surf,zeta_surf,trigon_factors)
            CHCKERR('')
            ierr = fourier2real(R_c(:,:,min_r:max_r,0),&
                &R_s(:,:,min_r:max_r,0),trigon_factors,&
                &R(id,:,:),[0,0])
            CHCKERR('')
            ierr = fourier2real(Z_c(:,:,min_r:max_r,0),&
                &Z_s(:,:,min_r:max_r,0),trigon_factors,&
                &Z(id,:,:),[0,0])
            CHCKERR('')
            
            ! transform cylindrical to cartesian
            X(id,:,:) = R(id,:,:)*cos(zeta_surf)
            Y(id,:,:) = R(id,:,:)*sin(zeta_surf)
            
            ! deallocate
            deallocate(trigon_factors)
        end do
        
        ! write output files
        plot_title(1) = 'Magnetic Flux Surfaces'
        file_name(1) = 'Flux_surfaces.dat'
        call print_GP_3D(trim(plot_title(1)),trim(file_name(1)),Z,x=X,y=Y,&
            &draw=.false.)
        
        ! deallocate
        deallocate(R,Z,X,Y,theta_surf,zeta_surf)
        
        ! 2. plot field lines
        
        ! initialize loc_n_par and loc_n_r
        loc_n_par = size(theta,1)
        loc_n_r = size(theta,2)
        
        ! initialize R and Z
        allocate(R(1,loc_n_par,loc_n_r),Z(1,loc_n_par,loc_n_r))
        allocate(X(1,loc_n_par,loc_n_r),Y(1,loc_n_par,loc_n_r))
        
        ! inverse fourier transform with trigonometric factors
        ierr = calc_trigon_factors(theta,zeta,trigon_factors)
        CHCKERR('')
        ierr = fourier2real(R_c(:,:,min_r:max_r,0),&
            &R_s(:,:,min_r:max_r,0),trigon_factors,R(1,:,:),[0,0])
        CHCKERR('')
        ierr = fourier2real(Z_c(:,:,min_r:max_r,0),&
            &Z_s(:,:,min_r:max_r,0),trigon_factors,Z(1,:,:),[0,0])
        CHCKERR('')
        
        ! transform cylindrical to cartesian
        X(1,:,:) = R(1,:,:)*cos(zeta)
        Y(1,:,:) = R(1,:,:)*sin(zeta)
        
        ! write output files
        plot_title(2) = 'Magnetic field'
        file_name(2) = 'B_field.dat'
        call print_GP_3D(trim(plot_title(2)),trim(file_name(2)),Z,x=X,y=Y,&
            &draw=.false.)
        
        ! deallocate
        deallocate(R,Z,X,Y)
        deallocate(trigon_factors)
        
        ! 3. plot both output files
        draw_ops(1) = 'linecolor rgb ''#d3d3d3'' linewidth 1'
        draw_ops(2) = 'linecolor rgb ''black'' linewidth 3'
        call draw_GP_animated('Magnetic field in flux surfaces',file_name,&
            &max_r-min_r+1,.false.,delay=50,draw_ops=draw_ops)
        
        call lvl_ud(-1)
    end function plot_grid_real
end module eq_ops
