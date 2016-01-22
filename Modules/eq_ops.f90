!------------------------------------------------------------------------------!
!   Operations on the equilibrium variables                                    !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln
    use grid_vars, only: grid_type, dealloc_grid
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    
    implicit none
    private
    public read_eq, calc_eq, calc_derived_q, calc_normalization_const, &
        &normalize_input, print_output_eq, flux_q_plot
#if ldebug
    public debug_calc_derived_q
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_derived_q = .false.                                   ! plot debug information for calc_derived_q
#endif
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! reads the equilibrium input file
    integer function read_eq() result(ierr)
        use num_vars, only: eq_style, rank, use_pol_flux_E
        use VMEC, only: read_VMEC
        use HELENA, only: read_HEL
        use grid_vars, only: n_r_eq
        
        character(*), parameter :: rout_name = 'read_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only do this for the master
        ! (The other ranks don't yet know what eq_style is!)
        if (rank.eq.0) then
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = read_VMEC(n_r_eq,use_pol_flux_E)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ierr = read_HEL(n_r_eq,use_pol_flux_E)
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
    end function read_eq
    
    ! calculate  the equilibrium  quantities on  a grid  determined by  straight
    ! field lines. This  grid has the dimensions (n_par,loc_n_r)  where n_par is
    ! the number of points taken along  the magnetic field lines and loc_n_r.le.
    ! n_r_eq  is  the normal  extent  in  the  equilibrium  grid of  this  rank,
    ! determined by sharing the workload evenly over the processes.
    integer function calc_eq(grid_eq,eq) result(ierr)
        use eq_vars, only: create_eq
        use num_vars, only: eq_style, max_deriv
        use grid_ops, only: calc_norm_range, setup_grid_eq, calc_ang_grid_eq
        use utilities, only: derivs
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        
        ! local variables
        integer :: id
        integer :: eq_limits(2)                                                 ! min. and max. index of eq. grid of this process
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Start setting up equilibrium quantities')
        
        call lvl_ud(1)
        
        ! split  the  communicator MPI_Comm_world into subcommunicators
        call writo('Calculate normal equilibrium ranges')
        call lvl_ud(1)
        ierr = calc_norm_range(eq_limits=eq_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! setup equilibrium grid
        call writo('Determine the equilibrium grid')
        call lvl_ud(1)
        ierr = setup_grid_eq(grid_eq,eq_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! create equilibrium
        call writo('Initialize equilibrium quantities')
        call lvl_ud(1)
        ierr = create_eq(grid_eq,eq)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! calculate flux quantities and complete equilibrium grid
        call writo('Calculate flux quantities')
        call lvl_ud(1)
        ierr = calc_flux_q(grid_eq,eq)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! calculate angular grid points for equilibrium grid
        call writo('Calculate angular equilibrium grid')
        call lvl_ud(1)
        ierr = calc_ang_grid_eq(grid_eq,eq)
        CHCKERR('')
        call lvl_ud(-1)
        
        call lvl_ud(-1)
        call writo('Done setting up equilibrium quantities')
        
        ! Calcalations depending on equilibrium style being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate  the  cylindrical  variables  R, Z  and  lambda  and
                ! derivatives
                call writo('Calculate R,Z,L...')
                ierr = prepare_RZL(grid_eq)
                CHCKERR('')
                do id = 0,max_deriv+1
                    ierr = calc_RZL(grid_eq,eq,derivs(id))
                    CHCKERR('')
                end do
            case (2)                                                            ! HELENA
                ! do nothing
        end select
    end function calc_eq
    
    ! Calculates derived  equilibrium quantities  in the  Equilibrium coordinate
    ! system [ADD REF]:
    !   - magnetic shear S
    !   - normal curvature kappa_n
    !   - geodesic curvature kappa_g
    !   - parallel current sigma
    subroutine calc_derived_q(grid_eq,eq,met)
        use utilities, only: c
        use eq_vars, only: vac_perm
        use num_vars, only: norm_disc_prec_eq
#if ldebug
        use num_vars, only: use_pol_flux_F
        use utilities, only: calc_deriv
        use grid_vars, only: dealloc_grid
        use grid_ops, only: trim_grid, calc_XYZ_grid
#endif
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:) => null()                                 ! jac
        real(dp), pointer :: D1J(:,:,:) => null()                               ! D_alpha jac
        real(dp), pointer :: D2J(:,:,:) => null()                               ! D_psi jac
        real(dp), pointer :: D3J(:,:,:) => null()                               ! D_theta jac
        ! lower metric factors
        real(dp), pointer :: g13(:,:,:) => null()                               ! g_alpha,theta
        real(dp), pointer :: D2g13(:,:,:) => null()                             ! D_psi g_alpha,theta
        real(dp), pointer :: D3g13(:,:,:) => null()                             ! D_theta g_alpha,theta
        real(dp), pointer :: g23(:,:,:) => null()                               ! g_psi,theta
        real(dp), pointer :: D1g23(:,:,:) => null()                             ! D_alpha g_psi,theta
        real(dp), pointer :: D3g23(:,:,:) => null()                             ! D_theta g_psi,theta
        real(dp), pointer :: g33(:,:,:) => null()                               ! g_theta,theta
        real(dp), pointer :: D1g33(:,:,:) => null()                             ! D_alpha g_theta,theta
        real(dp), pointer :: D2g33(:,:,:) => null()                             ! D_psi g_theta,theta
        real(dp), pointer :: D3g33(:,:,:) => null()                             ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), pointer :: h12(:,:,:) => null()                               ! h^alpha,psi
        real(dp), pointer :: D3h12(:,:,:) => null()                             ! D_theta h^alpha,psi
        real(dp), pointer :: h22(:,:,:) => null()                               ! h^psi,psi
        real(dp), pointer :: D3h22(:,:,:) => null()                             ! D_theta h^psi,psi
        real(dp), pointer :: h23(:,:,:) => null()                               ! h^psi,theta
#if ldebug
        type(grid_type) :: grid_eq_trim                                         ! trimmed equilibrium grid
        real(dp), allocatable :: D3sigma(:,:,:)                                 ! D_theta sigma
        real(dp), allocatable :: D3sigma_ALT(:,:,:)                             ! alternative D_theta sigma
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle theta_F or zeta_F
        real(dp), allocatable :: X_plot(:,:,:)                                  ! X of plot
        real(dp), allocatable :: Y_plot(:,:,:)                                  ! Y of plot
        real(dp), allocatable :: Z_plot(:,:,:)                                  ! Z of plot
        integer :: istat                                                        ! status
        integer :: jd                                                           ! counter
#endif
        
        ! user output
        call writo('Start setting up derived quantities in flux coordinates')
        
        call lvl_ud(1)
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        D1J => met%jac_FD(:,:,:,1,0,0)
        D2J => met%jac_FD(:,:,:,0,1,0)
        D3J => met%jac_FD(:,:,:,0,0,1)
        ! lower metric factors
        g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,0)
        D2g13 => met%g_FD(:,:,:,c([1,3],.true.),0,1,0)
        D3g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,1)
        g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,0)
        D1g23 => met%g_FD(:,:,:,c([2,3],.true.),1,0,0)
        D3g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,1)
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        D1g33 => met%g_FD(:,:,:,c([3,3],.true.),1,0,0)
        D2g33 => met%g_FD(:,:,:,c([3,3],.true.),0,1,0)
        D3g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,1)
        ! upper metric factors
        h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        D3h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,1)
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        D3h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,1)
        h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,0)
        
        ! Calculate the shear S
        eq%S = -(D3h12/h22 - D3h22*h12/h22**2)/J
        
        ! Calculate the normal curvature kappa_n
        do kd = 1,grid_eq%loc_n_r
            eq%kappa_n(:,:,kd) = &
                &vac_perm*J(:,:,kd)**2*eq%pres_FD(kd,1)/g33(:,:,kd) + &
                &1._dp/(2*h22(:,:,kd)) * ( &
                &h12(:,:,kd) * ( D1g33(:,:,kd)/g33(:,:,kd) - &
                &2*D1J(:,:,kd)/J(:,:,kd) ) + &
                &h22(:,:,kd) * ( D2g33(:,:,kd)/g33(:,:,kd) - &
                &2*D2J(:,:,kd)/J(:,:,kd) ) + &
                &h23(:,:,kd) * ( D3g33(:,:,kd)/g33(:,:,kd) - &
                &2*D3J(:,:,kd)/J(:,:,kd) ) )
        end do
        
        ! Calculate the geodesic curvature kappa_g
        eq%kappa_g(:,:,:) = (0.5*D1g33/g33 - D1J/J) - &
            &g13/g33*(0.5*D3g33/g33 - D3J/J)
        
        ! Calculate the parallel current sigma
        write(*,*) 'WHY IS THERE A SPIKE OF SIGMA AT LAST NORMAL POINT ?!?!'
        eq%sigma = 1._dp/vac_perm*&
            &(D1g23/J - g23*D1J/J**2 - D2g13/J + g13*D2J/J**2)
        do kd = 1,grid_eq%loc_n_r
            eq%sigma(:,:,kd) = eq%sigma(:,:,kd) - &
                &eq%pres_FD(kd,1)*J(:,:,kd)*g13(:,:,kd)/g33(:,:,kd)
        end do 
        
#if ldebug
        ! test whether -2 p' J kappa_g = D3sigma and plot kappa components
        if (debug_calc_derived_q) then
            call writo('Testing whether -2 p'' J kappa_g = D3sigma')
            call lvl_ud(1)
            
            ! trim equilibrium grid
            istat = trim_grid(grid_eq,grid_eq_trim,norm_id)
            CHCKSTT
            
            ! allocate variables
            allocate(D3sigma(grid_eq_trim%n(1),grid_eq_trim%n(2),&
                &grid_eq_trim%loc_n_r))
            allocate(D3sigma_ALT(grid_eq_trim%n(1),grid_eq_trim%n(2),&
                &grid_eq_trim%loc_n_r))
            
            ! point parallel angle
            if (use_pol_flux_F) then
                ang_par_F => grid_eq%theta_F
            else
                ang_par_F => grid_eq%zeta_F
            end if
            
            ! get derived sigma
            do kd = norm_id(1),norm_id(2)
                do jd = 1,grid_eq_trim%n(2)
                    istat = calc_deriv(eq%sigma(:,jd,kd),&
                        &D3sigma(:,jd,kd-norm_id(1)+1),ang_par_F(:,jd,kd),1,&
                        &norm_disc_prec_eq)
                    CHCKSTT
                end do
            end do
            
            ! calculate alternatively derived sigma
            do kd = norm_id(1),norm_id(2)
                D3sigma_ALT(:,:,kd-norm_id(1)+1) = &
                    &-2*eq%pres_FD(kd,1)*eq%kappa_g(:,:,kd)*J(:,:,kd)
            end do
            
            ! plot output
            call plot_diff_HDF5(D3sigma,D3sigma_ALT,'TEST_D3sigma',&
                &grid_eq_trim%n,[0,0,grid_eq_trim%i_min-1],&
                &description='To test whether -2 p'' J kappa_g = D3sigma',&
                &output_message=.true.)
            
            ! get X, Y and Z of plot
            allocate(X_plot(grid_eq_trim%n(1),grid_eq_trim%n(2),&
                &grid_eq_trim%loc_n_r))
            allocate(Y_plot(grid_eq_trim%n(1),grid_eq_trim%n(2),&
                &grid_eq_trim%loc_n_r))
            allocate(Z_plot(grid_eq_trim%n(1),grid_eq_trim%n(2),&
                &grid_eq_trim%loc_n_r))
            istat = calc_XYZ_grid(grid_eq_trim,X_plot,Y_plot,Z_plot)
            CHCKSTT
            
            ! plot shear
            call plot_HDF5('shear','TEST_shear',&
                &eq%S(:,:,norm_id(1):norm_id(2)),tot_dim=grid_eq_trim%n,&
                &loc_offset=[0,0,grid_eq_trim%i_min-1],x=X_plot,y=Y_plot,&
                &z=Z_plot)
            
            ! plot sigma
            call plot_HDF5('sigma','TEST_sigma',&
                &eq%sigma(:,:,norm_id(1):norm_id(2)),tot_dim=grid_eq_trim%n,&
                &loc_offset=[0,0,grid_eq_trim%i_min-1],x=X_plot,y=Y_plot,&
                &z=Z_plot)
            
            ! plot kappa_n
            call plot_HDF5('kappa_n','TEST_kappa_n',&
                &eq%kappa_n(:,:,norm_id(1):norm_id(2)),tot_dim=grid_eq_trim%n,&
                &loc_offset=[0,0,grid_eq_trim%i_min-1],x=X_plot,y=Y_plot,&
                &z=Z_plot)
            
            ! plot kappa_g
            call plot_HDF5('kappa_g','TEST_kappa_g',&
                &eq%kappa_g(:,:,norm_id(1):norm_id(2)),tot_dim=grid_eq_trim%n,&
                &loc_offset=[0,0,grid_eq_trim%i_min-1],x=X_plot,y=Y_plot,&
                &z=Z_plot)
            
            ! clean up
            nullify(ang_par_F)
            call dealloc_grid(grid_eq_trim)
            
            call lvl_ud(-1)
            call writo('Testing done')
        end if
#endif
        
        ! clean up
        nullify(J,D1J,D2J,D3J)
        nullify(g13,D2g13,D3g13)
        nullify(g23,D1g23,D3g23)
        nullify(g33,D1g33,D2g33,D3g33)
        nullify(h12,D3h12)
        nullify(h22,D3h22)
        nullify(h23)
        
        call lvl_ud(-1)
        
        call writo('Done setting up derived quantities in flux coordinates')
    end subroutine calc_derived_q
    
    ! Calculates  flux  quantities and  normal  derivatives  in the  Equilibrium
    ! coordinate  system. Also  sets the  normal coordinate  in the  Equilibrium
    ! grid.
    integer function calc_flux_q(grid_eq,eq) result(ierr)
        use num_vars, only: eq_style, max_deriv, use_pol_flux_E, &
            &use_pol_flux_F, rho_style, use_normalization, norm_disc_prec_eq
        use utilities, only: calc_deriv, calc_int
        use eq_vars, only: max_flux_p_E, max_flux_t_E, max_flux_p_F, &
            &max_flux_t_F, rho_0
        
        character(*), parameter :: rout_name = 'calc_flux_q'
        
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_flux_q_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_flux_q_HEL()
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! Calculate particle density rho
        ! choose which density style is being used:
        !   1:  constant, equal to rho_0
        select case (rho_style)
            case (1)                                                            ! arbitrarily constant (normalized value)
                eq%rho = rho_0
            case default
                err_msg = 'No density style associated with '//&
                    &trim(i2str(rho_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        ! normalize rho if necessary
        if (use_normalization) eq%rho = eq%rho/rho_0
    contains
        ! VMEC version
        ! The VMEC normal coord. is  the toroidal (or poloidal) flux, normalized
        ! wrt.  to  the  maximum  flux,  equidistantly,  so  the  step  size  is
        ! 1/(n(3)-1)
        integer function calc_flux_q_VMEC() result(ierr)
            use VMEC, only: rot_t_V, flux_t_V, Dflux_t_V, pres_V
            
            character(*), parameter :: rout_name = 'calc_flux_q_VMEC'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_p_full(:)                            ! version of dflux_p/dp on full normal equilibrium grid (1..n(3))
            real(dp), allocatable :: flux_p_full(:)                             ! version of integrated flux_p on full normal grid (1..n(3))
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate poloidal flux
            allocate(Dflux_p_full(grid_eq%n(3)),flux_p_full(grid_eq%n(3)))
            Dflux_p_full = rot_t_V*Dflux_t_V
            ierr = calc_int(Dflux_p_full,1.0_dp/(grid_eq%n(3)-1.0_dp),&
                &flux_p_full)                                                   ! equidistant grid (0..1) with n(3) points
            CHCKERR('')
            
            ! max flux and normal coord. of eq grid in Equilibrium coordinates
            max_flux_p_E = flux_p_full(grid_eq%n(3))
            max_flux_t_E = flux_t_V(grid_eq%n(3))
            if (use_pol_flux_E) then
                grid_eq%r_E = flux_p_full/flux_p_full(grid_eq%n(3))
                grid_eq%loc_r_E = flux_p_full(grid_eq%i_min:grid_eq%i_max)/&
                    &flux_p_full(grid_eq%n(3))
            else
                grid_eq%r_E = flux_t_V/flux_t_V(grid_eq%n(3))
                grid_eq%loc_r_E = flux_t_V(grid_eq%i_min:grid_eq%i_max)/&
                    &flux_t_V(grid_eq%n(3))
            end if
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            max_flux_p_F = flux_p_full(grid_eq%n(3))
            max_flux_t_F = - flux_t_V(grid_eq%n(3))                             ! conversion VMEC LH -> RH coord. system
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_full/(2*pi)                                ! psi_F = flux_p/2pi
                grid_eq%loc_r_F = flux_p_full(grid_eq%i_min:grid_eq%i_max)/&
                    &(2*pi)                                                     ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = - flux_t_V/(2*pi)                                 ! conversion VMEC LH -> RH coord. system, psi_F = flux_t/2pi
                grid_eq%loc_r_F = - flux_t_V(grid_eq%i_min:grid_eq%i_max)/&
                    &(2*pi)                                                     ! conversion VMEC LH -> RH coord. system, psi_F = flux_t/2pi
            end if
            
            ! poloidal  flux: calculate  using rot_t_V  and flux_t_V,  Dflux_t_V
            ! (easier to use full normal  equilibrium grid flux_p because of the
            ! integral)
            eq%flux_p_E(:,0) = flux_p_full(grid_eq%i_min:grid_eq%i_max)
            eq%flux_p_E(:,1) = Dflux_p_full(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+2
                ierr = calc_deriv(eq%flux_p_E(:,1),eq%flux_p_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd-1,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! toroidal flux: copy from VMEC and derive
            eq%flux_t_E(:,0) = flux_t_V(grid_eq%i_min:grid_eq%i_max)
            eq%flux_t_E(:,1) = Dflux_t_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+2
                ierr = calc_deriv(eq%flux_t_E(:,1),eq%flux_t_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd-1,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! pressure: copy from VMEC and derive
            eq%pres_E(:,0) = pres_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(eq%pres_E(:,0),eq%pres_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! safety factor
            eq%q_saf_E(:,0) = 1.0_dp/rot_t_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%q_saf_E(:,0),eq%q_saf_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! rot. transform
            eq%rot_t_E(:,0) = rot_t_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%rot_t_E(:,0),eq%rot_t_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! deallocate helper variables
            deallocate(Dflux_p_full,flux_p_full)
        end function calc_flux_q_VMEC
        
        ! HELENA version
        ! The HELENA normal coord. is the poloidal flux divided by 2pi
        integer function calc_flux_q_HEL() result(ierr)
            use HELENA, only: qs, flux_p_H, pres_H
            
            character(*), parameter :: rout_name = 'calc_flux_q_HEL'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_t_full(:)                            ! version of dflux_t/dp on full normal equilibrium grid
            real(dp), allocatable :: flux_t_full(:)                             ! version of integrated flux_t on full normal equilibrium grid
            real(dp), allocatable :: Dflux_p_H(:)                               ! normal derivative of flux_p_H
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate toroidal flux
            ! calculate normal derivative of flux_p_H
            allocate(Dflux_p_H(grid_eq%n(3)))
            ierr = calc_deriv(flux_p_H,Dflux_p_H,flux_p_H,1,norm_disc_prec_eq)
            allocate(Dflux_t_full(grid_eq%n(3)),flux_t_full(grid_eq%n(3)))
            Dflux_t_full = qs*Dflux_p_H
            ierr = calc_int(Dflux_t_full,flux_p_H,flux_t_full)
            CHCKERR('')
            
            ! max flux and  normal coord. of eq grid  in Equilibrium coordinates
            ! (uses poloidal flux by default)
            max_flux_p_E = flux_p_H(grid_eq%n(3))
            max_flux_t_E = flux_t_full(grid_eq%n(3))
            grid_eq%r_E = flux_p_H/(2*pi)
            grid_eq%loc_r_E = flux_p_H(grid_eq%i_min:grid_eq%i_max)/(2*pi)
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            max_flux_p_F = flux_p_H(grid_eq%n(3))
            max_flux_t_F = flux_t_full(grid_eq%n(3))
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_H/(2*pi)                                   ! psi_F = flux_p/2pi
                grid_eq%loc_r_F = flux_p_H(grid_eq%i_min:grid_eq%i_max)/(2*pi)  ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = flux_t_full/(2*pi)                                ! psi_F = flux_t/2pi
                grid_eq%loc_r_F = flux_t_full(grid_eq%i_min:grid_eq%i_max)/&
                    &(2*pi)                                                     ! psi_F = flux_t/2pi
            end if
            
            ! poloidal flux: copy from HELENA and derive
            eq%flux_p_E(:,0) = flux_p_H(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%flux_p_E(:,0),eq%flux_p_E(:,kd),&
                    &grid_eq%loc_r_E,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! toroidal flux: calculate using  qs and flux_p_H, Dflux_p_H (easier
            ! to  use  full  normal  equilibrium  grid  flux_t  because  of  the
            ! integral)
            eq%flux_t_E(:,0) = flux_t_full(grid_eq%i_min:grid_eq%i_max)
            eq%flux_t_E(:,1) = Dflux_t_full(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_t_E(:,1),eq%flux_t_E(:,kd),&
                    &grid_eq%loc_r_E,kd-1,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! pressure: copy from HELENA and derive
            eq%pres_E(:,0) = pres_H(grid_eq%i_min:grid_eq%i_max)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(eq%pres_E(:,0),eq%pres_E(:,kd),&
                    &grid_eq%loc_r_E,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! safety factor
            eq%q_saf_E(:,0) = qs(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%q_saf_E(:,0),eq%q_saf_E(:,kd),&
                    &grid_eq%loc_r_E,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! rot. transform
            eq%rot_t_E(:,0) = 1.0_dp/qs(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%rot_t_E(:,0),eq%rot_t_E(:,kd),&
                    &grid_eq%loc_r_E,kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! deallocate helper variables
            deallocate(Dflux_t_full,flux_t_full)
        end function calc_flux_q_HEL
    end function calc_flux_q
    
    ! prepare the cosine  and sine factors that are used  in the inverse Fourier
    ! transformation of R, Z and L and calculate the derivatives.
    integer function prepare_RZL(grid) result(ierr)
        use num_vars, only: max_deriv, norm_disc_prec_eq
        use VMEC, only: calc_trigon_factors, &
            &mpol, ntor, R_V_c, Z_V_s, L_V_s, R_V_s, Z_V_c, L_V_c
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'prepare_RZL'
        
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid for which to prepare the trigoniometric factors
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! normal derivatives of these factors
        ! The VMEC normal coord. is  the toroidal (or poloidal) flux, normalized
        ! wrt.  to  the  maximum  flux,  equidistantly,  so  the  step  size  is
        ! 1/(grid_eq%n(3)-1).
        do kd = 1,max_deriv+1
            do jd = -ntor,ntor
                do id = 0,mpol-1
                    ierr = calc_deriv(R_V_c(id,jd,:,0),R_V_c(id,jd,:,kd),&
                        &grid%n(3)-1._dp,kd,norm_disc_prec_eq)
                    CHCKERR('')
                    ierr = calc_deriv(Z_V_s(id,jd,:,0),Z_V_s(id,jd,:,kd),&
                        &grid%n(3)-1._dp,kd,norm_disc_prec_eq)
                    CHCKERR('')
                    ierr = calc_deriv(L_V_s(id,jd,:,0),&
                        &L_V_s(id,jd,:,kd),grid%n(3)-1._dp,kd,norm_disc_prec_eq)
                    CHCKERR('')
                    !if (lasym) then                                            ! following only needed in assymetric situations
                        ierr = calc_deriv(R_V_s(id,jd,:,0),R_V_s(id,jd,:,kd),&
                            &grid%n(3)-1._dp,kd,norm_disc_prec_eq)
                        CHCKERR('')
                        ierr = calc_deriv(Z_V_c(id,jd,:,0),Z_V_c(id,jd,:,kd),&
                            &grid%n(3)-1._dp,kd,norm_disc_prec_eq)
                        CHCKERR('')
                        ierr = calc_deriv(L_V_c(id,jd,:,0),L_V_c(id,jd,:,kd),&
                            &grid%n(3)-1._dp,kd,norm_disc_prec_eq)
                        CHCKERR('')
                    !end if
                end do
            end do
        end do
        
        ! calculate trigoniometric factors using theta_E and zeta_E
        ierr = calc_trigon_factors(grid%theta_E,grid%zeta_E,grid%trigon_factors)
        CHCKERR('')
    end function prepare_RZL
    
    ! calculate  R, Z  and Lambda  and derivatives  in VMEC  coordinates at  the
    ! grid  points given  by  the  variables theta_E  and  zeta_E (contained  in
    ! trigon_factors) and at  every normal point. The  derivatives are indicated
    ! by the variable "deriv" which has 3 indices
    ! Note: There is no HELENA equivalent  because for HELENA simulations, R and
    ! Z are not necessary for calculation of the metric coefficients, and L does
    ! not exist.
    integer function calc_RZL_ind(grid,eq,deriv) result(ierr)
        use VMEC, only: fourier2real, &
            &R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s
        use utilities, only: check_deriv
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_RZL_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to prepare the trigoniometric factors
        type(eq_type), intent(inout) :: eq                                      ! equilibrium
        integer, intent(in) :: deriv(3)                                         ! derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv+1,'calc_RZL')
        CHCKERR('')
        
        ! calculate the variables R,Z and their angular derivative
        ierr = fourier2real(R_V_c(:,:,grid%i_min:grid%i_max,deriv(1)),&
            &R_V_s(:,:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(Z_V_c(:,:,grid%i_min:grid%i_max,deriv(1)),&
            &Z_V_s(:,:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(L_V_c(:,:,grid%i_min:grid%i_max,deriv(1)),&
            &L_V_s(:,:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%L_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
    end function calc_RZL_ind
    integer function calc_RZL_arr(grid,eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_RZL_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to prepare the trigoniometric factors
        type(eq_type), intent(inout) :: eq                                      ! equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_RZL_ind(grid,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_RZL_arr

    ! plots the flux quantities in the solution grid
    !   safety factor q_saf
    !   rotational transform rot
    !   pressure pres
    !   poloidal flux flux_p
    !   toroidal flux flux_t
    integer function flux_q_plot(grid_eq,eq) result(ierr)
        use num_vars, only: rank, no_plots
        use grid_ops, only: trim_grid
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'flux_q_plot'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! normal grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: id                                                           ! counter
        integer :: n_vars = 5                                                   ! nr. of variables to plot
        character(len=max_str_ln), allocatable :: plot_titles(:)                ! plot titles
        type(grid_type) :: grid_trim                                            ! trimmed grid
        real(dp), allocatable :: X_plot_2D(:,:)                                 ! x values of 2D plot
        real(dp), allocatable :: Y_plot_2D(:,:)                                 ! y values of 2D plot
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! user output
        call writo('Plotting flux quantities')
        
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! initialize plot titles and file names
        allocate(plot_titles(n_vars))
        plot_titles(1) = 'safety factor []'
        plot_titles(2) = 'rotational transform []'
        plot_titles(3) = 'pressure [pa]'
        plot_titles(4) = 'poloidal flux [Tm^2]'
        plot_titles(5) = 'toroidal flux [Tm^2]'
        
        ! plot using HDF5
        ierr = flux_q_plot_HDF5()
        CHCKERR('')
        
        ! plot using GNUPlot
        ierr = flux_q_plot_GP()
        CHCKERR('')
        
        ! clean up
        call dealloc_grid(grid_trim)
        
        call lvl_ud(-1)
    contains
        ! plots the pressure and fluxes in GNUplot
        integer function flux_q_plot_GP() result(ierr)
            use eq_vars, only: max_flux_p_F, max_flux_t_F, pres_0, psi_0
            use num_vars, only: use_pol_flux_F, use_normalization
            
            character(*), parameter :: rout_name = 'flux_q_plot_GP'
            
            ! local variables
            character(len=max_str_ln), allocatable :: file_name(:)              ! file_name
            real(dp), allocatable :: ser_var_loc(:)                             ! local serial var
            
            ! initialize ierr
            ierr = 0
            
            ! allocate variabels if master
            if (rank.eq.0) then
                allocate(X_plot_2D(grid_trim%n(3),n_vars))
                allocate(Y_plot_2D(grid_trim%n(3),n_vars))
            end if
            
            ! fill the 2D version of the plot
            ierr = get_ser_var(eq%q_saf_FD(norm_id(1):norm_id(2),0),ser_var_loc)
            CHCKERR('')
            if (rank.eq.0) Y_plot_2D(:,1) = ser_var_loc
            ierr = get_ser_var(eq%rot_t_FD(norm_id(1):norm_id(2),0),ser_var_loc)
            CHCKERR('')
            if (rank.eq.0) Y_plot_2D(:,2) = ser_var_loc
            ierr = get_ser_var(eq%pres_FD(norm_id(1):norm_id(2),0),ser_var_loc)
            CHCKERR('')
            if (rank.eq.0) Y_plot_2D(:,3) = ser_var_loc
            ierr = get_ser_var(eq%flux_p_FD(norm_id(1):norm_id(2),0),&
                &ser_var_loc)
            CHCKERR('')
            if (rank.eq.0) Y_plot_2D(:,4) = ser_var_loc
            ierr = get_ser_var(eq%flux_t_FD(norm_id(1):norm_id(2),0),&
                &ser_var_loc)
            CHCKERR('')
            if (rank.eq.0) Y_plot_2D(:,5) = ser_var_loc
            
            ! rescale if normalized
            if (use_normalization .and. rank.eq.0) then
                Y_plot_2D(:,3) = Y_plot_2D(:,3)*pres_0                          ! pressure
                Y_plot_2D(:,4) = Y_plot_2D(:,4)*psi_0                           ! flux_p
                Y_plot_2D(:,5) = Y_plot_2D(:,5)*psi_0                           ! flux_t
            end if
            
            ! continue the plot if global master
            if (rank.eq.0) then
                ! deallocate local serial variables
                deallocate(ser_var_loc)
                
                ! set up file names
                allocate(file_name(2))
                file_name(1) = 'pres'
                file_name(2) = 'flux'
                
                ! 2D normal variable (normalized F coordinate)
                if (use_pol_flux_F) then
                    X_plot_2D(:,1) = grid_trim%r_F*2*pi/max_flux_p_F
                else
                    X_plot_2D(:,1) = grid_trim%r_F*2*pi/max_flux_t_F
                end if
                do id = 2,n_vars
                    X_plot_2D(:,id) = X_plot_2D(:,1)
                end do
                
                ! plot the  individual 2D output  of this process  (except q_saf
                ! and rot_t, as they are already in plot_resonance)
                call print_GP_2D(plot_titles(3),file_name(1),&
                    &Y_plot_2D(:,3),X_plot_2D(:,3),draw=.false.)
                ! fluxes
                call print_GP_2D(trim(plot_titles(4))//', '//&
                    &trim(plot_titles(5)),file_name(2),&
                    &Y_plot_2D(:,4:5),X_plot_2D(:,4:5),draw=.false.)
                
                ! draw plot
                call draw_GP(plot_titles(3),file_name(1),file_name(1),1,1,&
                    &.false.)                                                   ! pressure
                call draw_GP(trim(plot_titles(4))//', '//trim(plot_titles(5)),&
                    &file_name(2),file_name(2),2,1,.false.)                     ! fluxes
                
                ! user output
                call writo('Safety factor and rotational transform are plotted &
                    &using "plot_resonance"')
                
                ! clean up
                deallocate(X_plot_2D,Y_plot_2D)
            end if
        end function flux_q_plot_GP
        
        ! plots the flux quantities in HDF5
        integer function flux_q_plot_HDF5() result(ierr)
            use output_ops, only: plot_HDF5
            use grid_ops, only: calc_XYZ_grid, extend_grid_E, trim_grid
            use grid_vars, only: dealloc_grid
            
            character(*), parameter :: rout_name = 'flux_q_plot_HDF5'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: X_plot(:,:,:,:)                            ! x values of total plot
            real(dp), allocatable :: Y_plot(:,:,:,:)                            ! y values of total plot
            real(dp), allocatable :: Z_plot(:,:,:,:)                            ! z values of total plot
            real(dp), allocatable :: f_plot(:,:,:,:)                            ! values of variable of total plot
            integer :: plot_dim(4)                                              ! total plot dimensions
            integer :: plot_offset(4)                                           ! plot offset
            type(grid_type) :: grid_plot                                        ! grid for plotting
            character(len=max_str_ln) :: file_name                              ! file name
            
            ! initialize ierr
            ierr = 0
            
            ! set up file name
            file_name = 'flux_quantities'
            
            ! fill the 2D version of the plot
            allocate(Y_plot_2D(grid_trim%loc_n_r,n_vars))
            
            Y_plot_2D(:,1) = eq%q_saf_FD(norm_id(1):norm_id(2),0)
            Y_plot_2D(:,2) = eq%rot_t_FD(norm_id(1):norm_id(2),0)
            Y_plot_2D(:,3) = eq%pres_FD(norm_id(1):norm_id(2),0)
            Y_plot_2D(:,4) = eq%flux_p_FD(norm_id(1):norm_id(2),0)
            Y_plot_2D(:,5) = eq%flux_t_FD(norm_id(1):norm_id(2),0)
            
            ! extend trimmed equilibrium grid
            ierr = extend_grid_E(grid_trim,grid_plot)
            CHCKERR('')
            
            ! set up plot_dim and plot_offset
            plot_dim = [grid_plot%n(1),grid_plot%n(2),grid_plot%n(3),n_vars]
            plot_offset = [0,0,grid_plot%i_min-1,n_vars]
            
            ! set up total plot variables
            allocate(X_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r,&
                &n_vars))
            allocate(Y_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r,&
                &n_vars))
            allocate(Z_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r,&
                &n_vars))
            allocate(f_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%loc_n_r,&
                &n_vars))
            
            ! calculate 3D X,Y and Z
            ierr = calc_XYZ_grid(grid_plot,X_plot(:,:,:,1),Y_plot(:,:,:,1),&
                &Z_plot(:,:,:,1))
            CHCKERR('')
            do id = 2,n_vars
                X_plot(:,:,:,id) = X_plot(:,:,:,1)
                Y_plot(:,:,:,id) = Y_plot(:,:,:,1)
                Z_plot(:,:,:,id) = Z_plot(:,:,:,1)
            end do
            
            do kd = 1,grid_plot%loc_n_r
                f_plot(:,:,kd,1) = Y_plot_2D(kd,1)                              ! safey factor
                f_plot(:,:,kd,2) = Y_plot_2D(kd,2)                              ! rotational transform
                f_plot(:,:,kd,3) = Y_plot_2D(kd,3)                              ! pressure
                f_plot(:,:,kd,4) = Y_plot_2D(kd,4)                              ! poloidal flux
                f_plot(:,:,kd,5) = Y_plot_2D(kd,5)                              ! toroidal flux
            end do
            
            ! print the output using HDF5
            call plot_HDF5(plot_titles,file_name,f_plot,plot_dim,plot_offset,&
                &X_plot,Y_plot,Z_plot,col=1,description='Flux quantities')
            
            ! deallocate and destroy grid
            deallocate(Y_plot_2D)
            deallocate(X_plot,Y_plot,Z_plot,f_plot)
            call dealloc_grid(grid_plot)
        end function flux_q_plot_HDF5
    end function flux_q_plot
    
    ! sets up normalization constants:
    !   R_0:    major radius (= average R on axis)
    !   rho_0:  mass density on axis (set up through input variable)
    !   pres_0: pressure on axis
    !   B_0:    average magnetic field (= sqrt(mu_0 pres_0))
    !   psi_0:  reference flux (= R_0^2 B_0)
    ! [MPI] only global master
    integer function calc_normalization_const() result(ierr)
        use num_vars, only: rank, eq_style, mu_0_original, use_normalization
        use eq_vars, only: T_0, B_0, pres_0, psi_0, R_0, rho_0
        
        character(*), parameter :: rout_name = 'calc_normalization_const'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        if (use_normalization) then
            ! user output
            call writo('Calculating the normalization constants')
            call lvl_ud(1)
            
            ! calculation
            if (rank.eq.0) then
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        call calc_normalization_const_VMEC
                    case (2)                                                    ! HELENA
                        call calc_normalization_const_HEL
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            end if
            
            ! Alfven velocity
            T_0 = sqrt(mu_0_original*rho_0)*R_0/B_0 
            
            call writo('Major radius    R_0    = '//trim(r2strt(R_0))//' m')
            call writo('Pressure        pres_0 = '//trim(r2strt(pres_0))//' Pa')
            call writo('Mass density    rho_0  = '//trim(r2strt(rho_0))&
                &//' kg/m^3')
            call writo('Magnetic field  B_0    = '//trim(r2strt(B_0))//' T')
            call writo('Magnetic flux   psi_0  = '//trim(r2strt(psi_0))//&
                &' Tm^2')
            call writo('Alfven time     T_0    = '//trim(r2strt(T_0))//' s')
            call writo('Vacuum perm.    mu_0   = '//trim(r2strt(mu_0_original))&
                &//' Tm/A')
            
            ! user output
            call lvl_ud(-1)
            call writo('Normalization constants calculated')
        else
            ! user output
            call writo('Normalization not used')
        end if
    contains 
        ! VMEC version
        subroutine calc_normalization_const_VMEC
            use VMEC, only: R_V_c, pres_V
            
            ! set the major  radius as the average value of  R_E on the magnetic
            ! axis
            R_0 = R_V_c(0,0,1,0)
            
            ! rho_0 is set up through an input variable with the same name
            
            ! set pres_0 as pressure on axis
            pres_0 = pres_V(1)
            
            ! set  the  reference value  for B_0  from B_0  = sqrt(mu_0_original
            ! pres_0)
            B_0 = sqrt(pres_0 * mu_0_original)
            
            ! set reference flux
            psi_0 = R_0**2 * B_0
        end subroutine calc_normalization_const_VMEC
        
        ! HELENA version
        subroutine calc_normalization_const_HEL
            use HELENA, only: R_0_H, B_0_H
            
            ! set the major radius as the HELENA normalization parameter
            R_0 = R_0_H
            
            ! rho_0 is set up through an input variable with the same name
            
            !  set the  reference  value  for B_0  as  the HELENA  normalization
            ! parameter
            B_0 = B_0_H
            
            ! set pres_0 as B_0^2/mu_0_original
            pres_0 = B_0**2/mu_0_original
            
            ! set reference flux
            psi_0 = R_0**2 * B_0
        end subroutine calc_normalization_const_HEL
    end function calc_normalization_const
    
    ! normalize input quantities
    integer function normalize_input() result(ierr)
        use num_vars, only: use_normalization, eq_style, mu_0_original, rank
        use VMEC, only: normalize_VMEC
        use eq_vars, only: vac_perm
        
        character(*), parameter :: rout_name = 'normalize_input'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only normalize if needed
        if (use_normalization .and. rank.eq.0) then
            ! user output
            call writo('Start normalizing the input variables')
            call lvl_ud(1)
            
            ! normalize common variables
            vac_perm = vac_perm/mu_0_original
            
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    call normalize_VMEC
                case (2)                                                        ! HELENA
                    ! HELENA input already normalized
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! user output
            call lvl_ud(-1)
            call writo('Normalization done')
        end if
    end function normalize_input
    
    ! Print equilibrium quantities to an output file:
    !   - grid:     r_F, theta_F, zeta_F
    !   - eq:       pres_FD, q_saf_FD, rot_t_FD, flux_p_FD, flux_t_FD, rho, S,
    !               kappa_n, kappa_g, sigma
    !   - metric:   g_FD, h_FD, jac_FD
    !   - VMEC:     R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s
    !   - HELENA:   R_H, Z_H
    ! Note: The equilibrium quantities are outputted in Flux coordinates.
    integer function print_output_eq(grid_eq,eq,met) result(ierr)
        use num_vars, only: eq_style, PB3D_name, rank
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        use grid_ops, only: trim_grid
        use HELENA, only: chi_H, flux_p_H, R_H, Z_H, nchi
        use VMEC, only:  R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, ntor
        use eq_vars, only: max_flux_p_E, max_flux_t_E, max_flux_p_F, &
            &max_flux_t_F
        
        character(*), parameter :: rout_name = 'print_output_eq'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: norm_id_f(2)                                                 ! untrimmed full normal indices for trimmed grids
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(var_1D_type), allocatable, target :: eq_1D(:)                      ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: eq_1D_loc => null()                       ! local element in eq_1D
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Writing equilibrium variables to output file')
        call lvl_ud(1)
        
        ! user output
        call writo('Preparing variables for writing')
        call lvl_ud(1)
        
        ! trim grids
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set i_min and i_max for variables tabulated on full grid, trimmed
        norm_id_f = norm_id + grid_eq%i_min - 1
        
        ! Set up the 1D equivalents  of the equilibrium variables, with size
        ! depending on equilibrium style
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                allocate(eq_1D(25))
            case (2)                                                            ! HELENA
                allocate(eq_1D(18))
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! Set up common variables eq_1D
        id = 1
        
        ! max_flux
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'max_flux'
        allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
        allocate(eq_1D_loc%loc_i_min(1),eq_1D_loc%loc_i_max(1))
        if (rank.eq.0) then
            eq_1D_loc%loc_i_min = [1]
            eq_1D_loc%loc_i_max = [4]
            allocate(eq_1D_loc%p(4))
            eq_1D_loc%p = [max_flux_p_E,max_flux_t_E,max_flux_p_F,max_flux_t_F]
        else
            eq_1D_loc%loc_i_min = [1]
            eq_1D_loc%loc_i_max = [0]
            allocate(eq_1D_loc%p(0))
        end if
        eq_1D_loc%tot_i_min = [1]
        eq_1D_loc%tot_i_max = [4]
        
        ! pres_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'pres_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%pres_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%pres_FD,2)-1]
        allocate(eq_1D_loc%p(size(eq%pres_FD(norm_id(1):norm_id(2),:))))
        eq_1D_loc%p = reshape(eq%pres_FD(norm_id(1):norm_id(2),:),&
            &[size(eq%pres_FD(norm_id(1):norm_id(2),:))])
        
        ! q_saf_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'q_saf_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%q_saf_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%q_saf_FD,2)-1]
        allocate(eq_1D_loc%p(size(eq%q_saf_FD(norm_id(1):norm_id(2),:))))
        eq_1D_loc%p = reshape(eq%q_saf_FD(norm_id(1):norm_id(2),:),&
            &[size(eq%q_saf_FD(norm_id(1):norm_id(2),:))])
        
        ! rot_t_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'rot_t_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%rot_t_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%rot_t_FD,2)-1]
        allocate(eq_1D_loc%p(size(eq%rot_t_FD(norm_id(1):norm_id(2),:))))
        eq_1D_loc%p = reshape(eq%rot_t_FD(norm_id(1):norm_id(2),:),&
            &[size(eq%rot_t_FD(norm_id(1):norm_id(2),:))])
        
        ! flux_p_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'flux_p_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%flux_p_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%flux_p_FD,2)-1]
        allocate(eq_1D_loc%p(size(eq%flux_p_FD(norm_id(1):norm_id(2),:))))
        eq_1D_loc%p = reshape(eq%flux_p_FD(norm_id(1):norm_id(2),:),&
            &[size(eq%flux_p_FD(norm_id(1):norm_id(2),:))])
        
        ! flux_t_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'flux_t_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%flux_t_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%flux_t_FD,2)-1]
        allocate(eq_1D_loc%p(size(eq%flux_t_FD(norm_id(1):norm_id(2),:))))
        eq_1D_loc%p = reshape(eq%flux_t_FD(norm_id(1):norm_id(2),:),&
            &[size(eq%flux_t_FD(norm_id(1):norm_id(2),:))])
        
        ! rho
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'rho'
        allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
        allocate(eq_1D_loc%loc_i_min(1),eq_1D_loc%loc_i_max(1))
        eq_1D_loc%tot_i_min = 1
        eq_1D_loc%tot_i_max = grid_trim%n(3)
        eq_1D_loc%loc_i_min = grid_trim%i_min
        eq_1D_loc%loc_i_max = grid_trim%i_max
        allocate(eq_1D_loc%p(size(eq%rho(norm_id(1):norm_id(2)))))
        eq_1D_loc%p = eq%rho(norm_id(1):norm_id(2))
        
        ! S
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'S'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = &
            &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
        allocate(eq_1D_loc%p(size(eq%S(:,:,norm_id(1):norm_id(2)))))
        eq_1D_loc%p = reshape(eq%S(:,:,norm_id(1):norm_id(2)),&
            &[size(eq%S(:,:,norm_id(1):norm_id(2)))])
        
        ! kappa_n
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'kappa_n'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = &
            &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
        allocate(eq_1D_loc%p(size(eq%kappa_n(:,:,norm_id(1):norm_id(2)))))
        eq_1D_loc%p = reshape(eq%kappa_n(:,:,norm_id(1):norm_id(2)),&
            &[size(eq%kappa_n(:,:,norm_id(1):norm_id(2)))])
        
        ! kappa_g
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'kappa_g'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = &
            &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
        allocate(eq_1D_loc%p(size(eq%kappa_g(:,:,norm_id(1):norm_id(2)))))
        eq_1D_loc%p = reshape(eq%kappa_g(:,:,norm_id(1):norm_id(2)),&
            &[size(eq%kappa_g(:,:,norm_id(1):norm_id(2)))])
        
        ! sigma
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'sigma'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = [grid_trim%n(1:2),grid_trim%i_max]
        allocate(eq_1D_loc%p(size(eq%sigma(:,:,norm_id(1):norm_id(2)))))
        eq_1D_loc%p = reshape(eq%sigma(:,:,norm_id(1):norm_id(2)),&
            &[size(eq%sigma(:,:,norm_id(1):norm_id(2)))])
        
        ! g_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'g_FD'
        allocate(eq_1D_loc%tot_i_min(7),eq_1D_loc%tot_i_max(7))
        allocate(eq_1D_loc%loc_i_min(7),eq_1D_loc%loc_i_max(7))
        eq_1D_loc%tot_i_min = [1,1,1,1,0,0,0]
        eq_1D_loc%tot_i_max = [grid_trim%n,6,size(met%g_FD,5)-1,&
            &size(met%g_FD,6)-1,size(met%g_FD,7)-1]
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min,1,0,0,0]
        eq_1D_loc%loc_i_max = [grid_trim%n(1),grid_trim%n(2),&
            &grid_trim%i_max,6,size(met%g_FD,5)-1,size(met%g_FD,6)-1,&
            &size(met%g_FD,7)-1]
        allocate(eq_1D_loc%p(size(met%g_FD(:,:,norm_id(1):norm_id(2),:,:,:,:))))
        eq_1D_loc%p = reshape(met%g_FD(:,:,norm_id(1):norm_id(2),:,:,:,:),&
            &[size(met%g_FD(:,:,norm_id(1):norm_id(2),:,:,:,:))])
        
        ! h_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'h_FD'
        allocate(eq_1D_loc%tot_i_min(7),eq_1D_loc%tot_i_max(7))
        allocate(eq_1D_loc%loc_i_min(7),eq_1D_loc%loc_i_max(7))
        eq_1D_loc%tot_i_min = [1,1,1,1,0,0,0]
        eq_1D_loc%tot_i_max = [grid_trim%n,6,size(met%h_FD,5)-1,&
            &size(met%h_FD,6)-1,size(met%h_FD,7)-1]
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min,1,0,0,0]
        eq_1D_loc%loc_i_max = [grid_trim%n(1),grid_trim%n(2),&
            &grid_trim%i_max,6,size(met%h_FD,5)-1,size(met%h_FD,6)-1,&
            &size(met%h_FD,7)-1]
        allocate(eq_1D_loc%p(size(met%h_FD(:,:,norm_id(1):norm_id(2),:,:,:,:))))
        eq_1D_loc%p = reshape(met%h_FD(:,:,norm_id(1):norm_id(2),:,:,:,:),&
            &[size(met%h_FD(:,:,norm_id(1):norm_id(2),:,:,:,:))])
        
        ! jac_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'jac_FD'
        allocate(eq_1D_loc%tot_i_min(6),eq_1D_loc%tot_i_max(6))
        allocate(eq_1D_loc%loc_i_min(6),eq_1D_loc%loc_i_max(6))
        eq_1D_loc%tot_i_min = [1,1,1,0,0,0]
        eq_1D_loc%tot_i_max = [grid_trim%n,size(met%jac_FD,4)-1,&
            &size(met%jac_FD,5)-1,size(met%jac_FD,6)-1]
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min,0,0,0]
        eq_1D_loc%loc_i_max = [grid_trim%n(1),grid_trim%n(2),&
            &grid_trim%i_max,size(met%jac_FD,4)-1,size(met%jac_FD,5)-1,&
            &size(met%jac_FD,6)-1]
        allocate(eq_1D_loc%p(size(met%jac_FD(:,:,norm_id(1):norm_id(2),:,:,:))))
        eq_1D_loc%p = reshape(met%jac_FD(:,:,norm_id(1):norm_id(2),:,:,:),&
            &[size(met%jac_FD(:,:,norm_id(1):norm_id(2),:,:,:))])
        
        ! Set up particular variables, depending on equilibrium style
        !   1: VMEC needs the flux quantities in E coords. and VMEC variables in
        !   order to  calculate the  equilibrium quantities for  different grids
        !   and in order to plot.
        !   2: HELENA  needs HELENA variables  to interpolate the  output tables
        !   and in order to plot.
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! No  need  for  grid variables  as they  are  identical to  the
                ! field-aligned grid variables written above
                
                ! pres_E
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'pres_E'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,0]
                eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%pres_E,2)-1]
                eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%pres_E,2)-1]
                allocate(eq_1D_loc%p(size(eq%pres_E(norm_id(1):norm_id(2),:))))
                eq_1D_loc%p = reshape(eq%pres_E(norm_id(1):norm_id(2),:),&
                    &[size(eq%pres_E(norm_id(1):norm_id(2),:))])
                
                ! q_saf_E
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'q_saf_E'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,0]
                eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%q_saf_E,2)-1]
                eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%q_saf_E,2)-1]
                allocate(eq_1D_loc%p(size(eq%q_saf_E(norm_id(1):norm_id(2),:))))
                eq_1D_loc%p = reshape(eq%q_saf_E(norm_id(1):norm_id(2),:),&
                    &[size(eq%q_saf_E(norm_id(1):norm_id(2),:))])
                
                ! rot_t_E
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'rot_t_E'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,0]
                eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%rot_t_E,2)-1]
                eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%rot_t_E,2)-1]
                allocate(eq_1D_loc%p(size(eq%rot_t_E(norm_id(1):norm_id(2),:))))
                eq_1D_loc%p = reshape(eq%rot_t_E(norm_id(1):norm_id(2),:),&
                    &[size(eq%rot_t_E(norm_id(1):norm_id(2),:))])
                
                ! flux_p_E
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'flux_p_E'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,0]
                eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%flux_p_E,2)-1]
                eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%flux_p_E,2)-1]
                allocate(eq_1D_loc%p(size(eq%flux_p_E(norm_id(1):&
                    &norm_id(2),:))))
                eq_1D_loc%p = reshape(eq%flux_p_E(norm_id(1):norm_id(2),:),&
                    &[size(eq%flux_p_E(norm_id(1):norm_id(2),:))])
                
                ! flux_t_E
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'flux_t_E'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,0]
                eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%flux_t_E,2)-1]
                eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%flux_t_E,2)-1]
                allocate(eq_1D_loc%p(size(eq%flux_t_E(norm_id(1):&
                    &norm_id(2),:))))
                eq_1D_loc%p = reshape(eq%flux_t_E(norm_id(1):norm_id(2),:),&
                    &[size(eq%flux_t_E(norm_id(1):norm_id(2),:))])
                
                ! R_V_c
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'R_V_c'
                allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                allocate(eq_1D_loc%loc_i_min(4),eq_1D_loc%loc_i_max(4))
                eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),&
                    &size(R_V_c,4)-1]
                eq_1D_loc%loc_i_min = [0,-ntor,grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = &
                    &[mpol-1,ntor,grid_trim%i_max,size(R_V_c,4)-1]
                allocate(eq_1D_loc%p(size(R_V_c(:,:,norm_id_f(1):&
                    &norm_id_f(2),:))))
                eq_1D_loc%p = reshape(R_V_c(:,:,norm_id_f(1):norm_id_f(2),:),&
                    &[size(R_V_c(:,:,norm_id_f(1):norm_id_f(2),:))])
                
                ! R_V_s
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'R_V_s'
                allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                allocate(eq_1D_loc%loc_i_min(4),eq_1D_loc%loc_i_max(4))
                eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),&
                    &size(R_V_s,4)-1]
                eq_1D_loc%loc_i_min = [0,-ntor,grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = &
                    &[mpol-1,ntor,grid_trim%i_max,size(R_V_s,4)-1]
                allocate(eq_1D_loc%p(size(R_V_s(:,:,norm_id_f(1):norm_id_f(2),:))))
                eq_1D_loc%p = reshape(R_V_s(:,:,norm_id_f(1):norm_id_f(2),:),&
                    &[size(R_V_s(:,:,norm_id_f(1):norm_id_f(2),:))])
                
                ! Z_V_c
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'Z_V_c'
                allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                allocate(eq_1D_loc%loc_i_min(4),eq_1D_loc%loc_i_max(4))
                eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),&
                    &size(Z_V_c,4)-1]
                eq_1D_loc%loc_i_min = [0,-ntor,grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = &
                    &[mpol-1,ntor,grid_trim%i_max,size(Z_V_c,4)-1]
                allocate(eq_1D_loc%p(size(Z_V_c(:,:,norm_id_f(1):norm_id_f(2),:))))
                eq_1D_loc%p = reshape(Z_V_c(:,:,norm_id_f(1):norm_id_f(2),:),&
                    &[size(Z_V_c(:,:,norm_id_f(1):norm_id_f(2),:))])
                
                ! Z_V_s
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'Z_V_s'
                allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                allocate(eq_1D_loc%loc_i_min(4),eq_1D_loc%loc_i_max(4))
                eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),&
                    &size(Z_V_s,4)-1]
                eq_1D_loc%loc_i_min = [0,-ntor,grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = &
                    &[mpol-1,ntor,grid_trim%i_max,size(Z_V_s,4)-1]
                allocate(eq_1D_loc%p(size(Z_V_s(:,:,norm_id_f(1):&
                    &norm_id_f(2),:))))
                eq_1D_loc%p = reshape(Z_V_s(:,:,norm_id_f(1):norm_id_f(2),:),&
                    &[size(Z_V_s(:,:,norm_id_f(1):norm_id_f(2),:))])
                
                ! L_V_c
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'L_V_c'
                allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                allocate(eq_1D_loc%loc_i_min(4),eq_1D_loc%loc_i_max(4))
                eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),&
                    &size(L_V_c,4)-1]
                eq_1D_loc%loc_i_min = [0,-ntor,grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = &
                    &[mpol-1,ntor,grid_trim%i_max,size(L_V_c,4)-1]
                allocate(eq_1D_loc%p(size(L_V_c(:,:,norm_id_f(1):&
                    &norm_id_f(2),:))))
                eq_1D_loc%p = reshape(L_V_c(:,:,norm_id_f(1):norm_id_f(2),:),&
                    &[size(L_V_c(:,:,norm_id_f(1):norm_id_f(2),:))])
                
                ! L_V_s
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'L_V_s'
                allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                allocate(eq_1D_loc%loc_i_min(4),eq_1D_loc%loc_i_max(4))
                eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),&
                    &size(L_V_s,4)-1]
                eq_1D_loc%loc_i_min = [0,-ntor,grid_trim%i_min,0]
                eq_1D_loc%loc_i_max = &
                    &[mpol-1,ntor,grid_trim%i_max,size(L_V_s,4)-1]
                allocate(eq_1D_loc%p(size(L_V_s(:,:,norm_id_f(1):&
                    &norm_id_f(2),:))))
                eq_1D_loc%p = reshape(L_V_s(:,:,norm_id_f(1):norm_id_f(2),:),&
                    &[size(L_V_s(:,:,norm_id_f(1):norm_id_f(2),:))])
            case (2)                                                            ! HELENA
                ! R_H
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'R_H'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,1]
                eq_1D_loc%tot_i_max = [nchi,grid_trim%n(3)]
                eq_1D_loc%loc_i_min = [1,grid_trim%i_min]
                eq_1D_loc%loc_i_max = [nchi,grid_trim%i_max]
                allocate(eq_1D_loc%p(size(R_H(:,grid_trim%i_min:&
                    &grid_trim%i_max))))
                eq_1D_loc%p = reshape(R_H(:,norm_id_f(1):norm_id_f(2)),&
                    &[size(R_H(:,norm_id_f(1):norm_id_f(2)))])
                
                ! Z_H
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'Z_H'
                allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
                eq_1D_loc%tot_i_min = [1,1]
                eq_1D_loc%tot_i_max = [nchi,grid_trim%n(3)]
                eq_1D_loc%loc_i_min = [1,grid_trim%i_min]
                eq_1D_loc%loc_i_max = [nchi,grid_trim%i_max]
                allocate(eq_1D_loc%p(size(Z_H(:,grid_trim%i_min:&
                    &grid_trim%i_max))))
                eq_1D_loc%p = reshape(Z_H(:,norm_id_f(1):norm_id_f(2)),&
                    &[size(Z_H(:,norm_id_f(1):norm_id_f(2)))])
                
                ! chi_H
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'chi_H'
                allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
                allocate(eq_1D_loc%loc_i_min(1),eq_1D_loc%loc_i_max(1))
                eq_1D_loc%tot_i_min = [1]
                eq_1D_loc%tot_i_max = [nchi]
                eq_1D_loc%loc_i_min = [1]
                if (rank.eq.0) then
                    eq_1D_loc%loc_i_max = [nchi]
                    allocate(eq_1D_loc%p(nchi))
                    eq_1D_loc%p = chi_H
                else
                    eq_1D_loc%loc_i_max = [0]
                    allocate(eq_1D_loc%p(0))
                end if
                
                ! flux_p_H
                eq_1D_loc => eq_1D(id); id = id+1
                eq_1D_loc%var_name = 'flux_p_H'
                allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
                allocate(eq_1D_loc%loc_i_min(1),eq_1D_loc%loc_i_max(1))
                eq_1D_loc%tot_i_min = [1]
                eq_1D_loc%tot_i_max = [grid_trim%n(3)]
                eq_1D_loc%loc_i_min = [grid_trim%i_min]
                eq_1D_loc%loc_i_max = [grid_trim%i_max]
                allocate(eq_1D_loc%p(size(flux_p_H(grid_trim%i_min:&
                    &grid_trim%i_max))))
                eq_1D_loc%p = flux_p_H(norm_id_f(1):norm_id_f(2))
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        call lvl_ud(-1)
        
        ! write
        call writo('Writing using HDF5')
        call lvl_ud(1)
        ierr = print_HDF5_arrs(eq_1D,PB3D_name,'eq')
        CHCKERR('')
        call lvl_ud(-1)
        
        ! deallocate
        deallocate(eq_1D)
        
        ! user output
        call lvl_ud(-1)
        call writo('Equilibrium variables written to output')
    end function print_output_eq
end module eq_ops
