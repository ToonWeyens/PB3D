!------------------------------------------------------------------------------!
!   Operations on the equilibrium variables                                    !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    
    implicit none
    private
    public read_eq, calc_normalization_const, normalize_input, calc_eq, &
        &calc_flux_q, print_output_eq
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! reads the equilibrium input file
    integer function read_eq() result(ierr)
        use num_vars, only: eq_style, glb_rank, use_pol_flux_E
        use VMEC, only: read_VMEC
        use HELENA, only: read_HEL
        use grid_vars, only: n_r_eq
        use input_ops, only: get_log, pause_prog
#if ldebug
        use num_vars, only: ltest
        use HELENA, only: test_metrics_H
#endif
        
        character(*), parameter :: rout_name = 'read_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only do this for the group master
        ! (The other ranks don't yet know what eq_style is!)
        if (glb_rank.eq.0) then
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
#if ldebug
                    if (ltest) then
                        call writo('Test consistency of metric factors?')
                        if(get_log(.false.,ind=.true.)) then
                            ierr = test_metrics_H(n_r_eq)
                            CHCKERR('')
                            call pause_prog(ind=.true.)
                        end if
                    end if
#endif
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
    end function read_eq
    
    ! calculate  the equilibrium  quantities on  a grid  determined by  straight
    ! field lines. This  grid has the dimensions  (n_par,grp_n_r_eq) where n_par
    ! is  the  number  of  points  taken along  the  magnetic  field  lines  and
    ! grp_n_r_eq .le.  n_r_eq is the  normal extent  in the equilibrium  grid of
    ! this rank. It is determined so  that the perturbation quantities that will
    ! be needed in this rank are fully covered, so no communication is necessary
    integer function calc_eq(grid_eq,eq,met) result(ierr)
        use met_ops, only: calc_g_C, calc_T_VC, calc_g_V, calc_T_VF, &
            &calc_inv_met, calc_g_F, calc_jac_C, calc_jac_V, calc_jac_F, &
            &transf_deriv, calc_jac_H, calc_T_HF, calc_h_H
        use met_vars, only: create_met
        use utilities, only: derivs
        use input_ops, only: get_log, pause_prog
        use num_vars, only: max_deriv, eq_style
        use MPI_utilities, only: wait_MPI
#if ldebug
        use num_vars, only: ltest
        use met_ops, only: test_T_EF, test_p, test_jac_F, test_g_V, &
            &test_D12h_H, test_B_F, test_jac_V, test_Dg_E
#endif
        
        use utilities, only: calc_det, calc_inv, calc_mult, c
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(met_type), intent(inout) :: met                                    ! metric variables
        
        ! local variables
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: pmone                                                        ! plus or minus one
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call writo('Start setting up equilibrium quantities in '//&
                    &trim(i2str(grid_eq%n(1)))//' parallel and '//&
                    &trim(i2str(grid_eq%n(3)))//' normal points')
            case (2)                                                            ! HELENA
                call writo('Start setting up equilibrium quantities in '//&
                    &trim(i2str(grid_eq%n(1)))//' poloidal and '//&
                    &trim(i2str(grid_eq%n(3)))//' normal points')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        call lvl_ud(1)
            
            call writo('Calculating quantities on equilibrium grid')
            call lvl_ud(1)
            
                ! initialize metric quantities
                call writo('Initialize metric quantities...')
                ierr = create_met(grid_eq,met)
                CHCKERR('')
                
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! calculate the  cylindrical variables  R, Z  and lambda
                        ! and derivatives
                        call writo('Calculate R,Z,L...')
                        ierr = prepare_RZL(grid_eq)
                        CHCKERR('')
                        do id = 0,max_deriv+1
                            ierr = calc_RZL(grid_eq,eq,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate  the metrics  in the  cylindrical coordinate
                        ! system
                        call writo('Calculate g_C...')                          ! h_C is not necessary
                        do id = 0,max_deriv
                            ierr = calc_g_C(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the  jacobian in the  cylindrical coordinate
                        ! system
                        call writo('Calculate jac_C...')
                        do id = 0,max_deriv
                            ierr = calc_jac_C(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the  transformation matrix  C(ylindrical) ->
                        ! V(MEC)
                        call writo('Calculate T_VC...')
                        do id = 0,max_deriv
                            ierr = calc_T_VC(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate  the metric  factors in the  VMEC coordinate
                        ! system
                        call writo('Calculate g_V...')
                        do id = 0,max_deriv
                            ierr = calc_g_V(met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the jacobian in the VMEC coordinate system
                        call writo('Calculate jac_V...')
                        do id = 0,max_deriv
                            ierr = calc_jac_V(met,derivs(id))
                            CHCKERR('')
                        end do
                        
#if ldebug
                        if (ltest) then
                            call writo('Test calculation of g_V?')
                            if(get_log(.false.)) then
                                ierr = test_g_V(grid_eq,eq,met)
                                CHCKERR('')
                                call pause_prog
                            end if
                            call writo('Test calculation of jac_V?')
                            if(get_log(.false.)) then
                                ierr = test_jac_V(grid_eq,met)
                                CHCKERR('')
                                call pause_prog
                            end if
                        end if
#endif
                        
                        ! calculate the transformation matrix V(MEC) -> F(lux)
                        call writo('Calculate T_VF...')
                        do id = 0,max_deriv
                            ierr = calc_T_VF(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! set up plus minus one
                        pmone = -1                                              ! conversion VMEC LH -> RH coord. system
                    case (2)                                                    ! HELENA
                        ! calculate the jacobian in the HELENA coordinate system
                        call writo('Calculate jac_H...')
                        do id = 0,max_deriv
                            ierr = calc_jac_H(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the metric factors  in the HELENA coordinate
                        ! system
                        call writo('Calculate h_H...')
                        do id = 0,max_deriv
                            ierr = calc_h_H(grid_eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
#if ldebug
                        if (ltest) then
                            call writo('Test calculation of D1 D2 h_H?')
                            if(get_log(.false.)) then
                                ierr = test_D12h_H(grid_eq,met)
                                CHCKERR('')
                                call pause_prog
                            end if
                        end if
#endif
                        
                        ! calculate the inverse g_H of the metric factors h_H
                        call writo('Calculate g_H...')
                        do id = 0,max_deriv
                            ierr = calc_inv_met(met%g_E,met%h_E,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the transformation matrix H(ELENA) -> F(lux)
                        call writo('Calculate T_HF...')
                        do id = 0,max_deriv
                            ierr = calc_T_HF(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! set up plus minus one
                        pmone = 1
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
#if ldebug
                if (ltest) then
                    call writo('Test calculation of the derivatives of g_E?')
                    if(get_log(.false.)) then
                        ierr = test_Dg_E(grid_eq,met)
                        CHCKERR('')
                        call pause_prog
                    end if
                    call writo('Test calculation of T_EF?')
                    if(get_log(.false.)) then
                        ierr = test_T_EF(grid_eq,eq,met)
                        CHCKERR('')
                        call pause_prog
                    end if
                end if
#endif
                
                ! calculate  the  inverse of  the transformation  matrix T_EF
                call writo('Calculate T_FE...')
                do id = 0,max_deriv
                    ierr = calc_inv_met(met%T_FE,met%T_EF,derivs(id))
                    CHCKERR('')
                    ierr = calc_inv_met(met%det_T_FE,met%det_T_EF,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the metric factors in the Flux coordinate system
                call writo('Calculate g_F...')
                do id = 0,max_deriv
                    ierr = calc_g_F(met,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the inverse h_F of the metric factors g_F
                call writo('Calculate h_F...')
                do id = 0,max_deriv
                    ierr = calc_inv_met(met%h_F,met%g_F,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the jacobian in the Flux coordinate system
                call writo('Calculate jac_F...')
                do id = 0,max_deriv
                    ierr = calc_jac_F(met,derivs(id))
                    CHCKERR('')
                end do
                
#if ldebug
                if (ltest) then
                    call writo('Test Jacobian in Flux coordinates?')
                    if(get_log(.false.)) then
                        ierr = test_jac_F(grid_eq,eq,met)
                        CHCKERR('')
                        call pause_prog
                    end if
                    call writo('Test calculation of B_F?')
                    if(get_log(.false.)) then
                        ierr = test_B_F(grid_eq,eq,met)
                        CHCKERR('')
                        call pause_prog
                    end if
                end if
#endif
                
                ! Transform  the derivatives in VMEC  coordinates to derivatives
                ! in the Flux coordinates.
                call writo('Calculate derivatives in Flux coordinates...')
                do id = 0,max_deriv
                    ! g_FD
                    ierr = transf_deriv(met%g_F,met%T_FE,met%g_FD,max_deriv,&
                        &derivs(id))
                    CHCKERR('')
                    
                    ! h_FD
                    ierr = transf_deriv(met%h_F,met%T_FE,met%h_FD,max_deriv,&
                        &derivs(id))
                    CHCKERR('')
                    
                    ! jac_FD
                    ierr = transf_deriv(met%jac_F,met%T_FE,met%jac_FD,&
                        &max_deriv,derivs(id))
                    CHCKERR('')
                    
                    ! pres_FD
                    ierr = transf_deriv(eq%pres_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%pres_FD(:,id),max_deriv,id)
                    CHCKERR('')
                        
                    ! flux_p_FD
                    ierr = transf_deriv(eq%flux_p_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%flux_p_FD(:,id),max_deriv,id)
                    CHCKERR('')
                        
                    ! flux_t_FD
                    ierr = transf_deriv(eq%flux_t_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%flux_t_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%flux_t_FD(:,id) = pmone * eq%flux_t_FD(:,id)             ! multiply by plus minus one
                    
                    ! q_saf_FD
                    ierr = transf_deriv(eq%q_saf_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%q_saf_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%q_saf_FD(:,id) = pmone * eq%q_saf_FD(:,id)               ! multiply by plus minus one
                    
                    ! rot_t_FD
                    ierr = transf_deriv(eq%rot_t_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%rot_t_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%rot_t_FD(:,id) = pmone * eq%rot_t_FD(:,id)               ! multiply by plus minus one
                end do
                
#if ldebug
                if (ltest) then
                    call writo('Test consistency with given pressure?')
                    if(get_log(.false.)) then
                        ierr = test_p(grid_eq,eq,met)
                        CHCKERR('')
                        call pause_prog
                    end if
                end if
#endif
                ! calculate derived equilibrium quantities
                call writo('Calculate derived equilibrium quantities...')
                ierr = calc_derived_q(grid_eq,eq,met)
                CHCKERR('')
                
            call lvl_ud(-1)
            call writo('Equilibrium quantities calculated on equilibrium grid')
            
        call lvl_ud(-1)
        call writo('Done setting up equilibrium quantities')
    end function calc_eq
    
    ! prepare the cosine  and sine factors that are used  in the inverse Fourier
    ! transformation of R, Z and L and derivatives
    integer function prepare_RZL(grid) result(ierr)
        use fourier_ops, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'prepare_RZL'
        
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid for which to prepare the trigoniometric factors
        
        ! initialize ierr
        ierr = 0
        
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
        use fourier_ops, only: fourier2real
        use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s
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
    
    ! Calculates  flux  quantities and  normal  derivatives  in the  Equilibrium
    ! coordinate  system. Also  sets the  normal coordinate  in the  Equilibrium
    ! grid.
    integer function calc_flux_q(eq,grid_eq) result(ierr)
        use num_vars, only: eq_style, max_deriv, grp_nr, use_pol_flux_E, &
            &use_pol_flux_F, plot_flux_q
        use utilities, only: calc_deriv, calc_int
        use eq_vars, only: max_flux_p_E, max_flux_t_E, max_flux_p_F, &
            &max_flux_t_F
        
        character(*), parameter :: rout_name = 'calc_flux_q'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        
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
        
        ! plot flux quantities if requested
        call lvl_ud(1)
        if (plot_flux_q .and. grp_nr.eq.0) then                                 ! only first group because it is the same for all the groups
            ierr = flux_q_plot(eq,grid_eq)
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        call lvl_ud(-1)
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
                grid_eq%grp_r_E = flux_p_full(grid_eq%i_min:grid_eq%i_max)/&
                    &flux_p_full(grid_eq%n(3))
            else
                grid_eq%r_E = flux_t_V/flux_t_V(grid_eq%n(3))
                grid_eq%grp_r_E = flux_t_V(grid_eq%i_min:grid_eq%i_max)/&
                    &flux_t_V(grid_eq%n(3))
            end if
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            max_flux_p_F = flux_p_full(grid_eq%n(3))
            max_flux_t_F = - flux_t_V(grid_eq%n(3))                             ! conversion VMEC LH -> RH coord. system
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_full/(2*pi)                                ! psi_F = flux_p/2pi
                grid_eq%grp_r_F = flux_p_full(grid_eq%i_min:grid_eq%i_max)/&
                    &(2*pi)                                                     ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = - flux_t_V/(2*pi)                                 ! conversion VMEC LH -> RH coord. system, psi_F = flux_t/2pi
                grid_eq%grp_r_F = - flux_t_V(grid_eq%i_min:grid_eq%i_max)/&
                    &(2*pi)                                                     ! conversion VMEC LH -> RH coord. system, psi_F = flux_t/2pi
            end if
            
            ! poloidal  flux: calculate  using rot_t_V  and flux_t_V,  Dflux_t_V
            ! (easier to use full normal  equilibrium grid flux_p because of the
            ! integral)
            eq%flux_p_E(:,0) = flux_p_full(grid_eq%i_min:grid_eq%i_max)
            eq%flux_p_E(:,1) = Dflux_p_full(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_p_E(:,1),eq%flux_p_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd-1,1)
                CHCKERR('')
            end do
            
            ! toroidal flux: copy from VMEC and derive
            eq%flux_t_E(:,0) = flux_t_V(grid_eq%i_min:grid_eq%i_max)
            eq%flux_t_E(:,1) = Dflux_t_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_t_E(:,1),eq%flux_t_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd-1,1)
                CHCKERR('')
            end do
            
            ! pressure: copy from VMEC and derive
            eq%pres_E(:,0) = pres_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(eq%pres_E(:,0),eq%pres_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! safety factor
            eq%q_saf_E(:,0) = 1.0_dp/rot_t_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%q_saf_E(:,0),eq%q_saf_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! rot. transform
            eq%rot_t_E(:,0) = rot_t_V(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%rot_t_E(:,0),eq%rot_t_E(:,kd),&
                    &grid_eq%n(3)-1._dp,kd,1)
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
            real(dp), allocatable :: Dflux_t_full(:)                            ! version of dflux_t/dp on full normal equilibrium grid (1..n(3))
            real(dp), allocatable :: flux_t_full(:)                             ! version of integrated flux_t on full normal equilibrium grid (1..n(3))
            real(dp), allocatable :: Dflux_p_H(:)                               ! normal derivative of flux_p_H
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate toroidal flux
            ! calculate normal derivative of flux_p_H
            allocate(Dflux_p_H(grid_eq%n(3)))
            ierr = calc_deriv(flux_p_H,Dflux_p_H,flux_p_H,1,1)
            allocate(Dflux_t_full(grid_eq%n(3)),flux_t_full(grid_eq%n(3)))
            Dflux_t_full = qs*Dflux_p_H
            ierr = calc_int(Dflux_t_full,flux_p_H,flux_t_full)
            CHCKERR('')
            
            ! max flux and  normal coord. of eq grid  in Equilibrium coordinates
            ! (uses poloidal flux by default)
            max_flux_p_E = flux_p_H(grid_eq%n(3))
            max_flux_t_E = flux_t_full(grid_eq%n(3))
            grid_eq%r_E = flux_p_H/(2*pi)
            grid_eq%grp_r_E = flux_p_H(grid_eq%i_min:grid_eq%i_max)/(2*pi)
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            max_flux_p_F = flux_p_H(grid_eq%n(3))
            max_flux_t_F = flux_t_full(grid_eq%n(3))
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_H/(2*pi)                                   ! psi_F = flux_p/2pi
                grid_eq%grp_r_F = flux_p_H(grid_eq%i_min:grid_eq%i_max)/(2*pi)  ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = flux_t_full/(2*pi)                                ! psi_F = flux_t/2pi
                grid_eq%grp_r_F = flux_t_full(grid_eq%i_min:grid_eq%i_max)/&
                    &(2*pi)                                                     ! psi_F = flux_t/2pi
            end if
            
            ! poloidal flux: copy from HELENA and derive
            eq%flux_p_E(:,0) = flux_p_H(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%flux_p_E(:,0),eq%flux_p_E(:,kd),&
                    &grid_eq%grp_r_E,kd,1)
                CHCKERR('')
            end do
            
            ! toroidal flux: calculate using  qs and flux_p_H, Dflux_p_H (easier
            ! to  use  full  normal  equilibrium  grid  flux_t  because  of  the
            ! integral)
            eq%flux_t_E(:,0) = flux_t_full(grid_eq%i_min:grid_eq%i_max)
            eq%flux_t_E(:,1) = Dflux_t_full(grid_eq%i_min:grid_eq%i_max)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(eq%flux_t_E(:,1),eq%flux_t_E(:,kd),&
                    &grid_eq%grp_r_E,kd-1,1)
                CHCKERR('')
            end do
            
            ! pressure: copy from HELENA and derive
            eq%pres_E(:,0) = pres_H(grid_eq%i_min:grid_eq%i_max)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(eq%pres_E(:,0),eq%pres_E(:,kd),&
                    &grid_eq%grp_r_E,kd,1)
                CHCKERR('')
            end do
            
            ! safety factor
            eq%q_saf_E(:,0) = qs(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%q_saf_E(:,0),eq%q_saf_E(:,kd),&
                    &grid_eq%grp_r_E,kd,1)
                CHCKERR('')
            end do
            
            ! rot. transform
            eq%rot_t_E(:,0) = 1.0_dp/qs(grid_eq%i_min:grid_eq%i_max)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(eq%rot_t_E(:,0),eq%rot_t_E(:,kd),&
                    &grid_eq%grp_r_E,kd,1)
                CHCKERR('')
            end do
            
            ! deallocate helper variables
            deallocate(Dflux_t_full,flux_t_full)
        end function calc_flux_q_HEL
    end function calc_flux_q
    
    ! Calculates derived  equilibrium quantities  in the  Equilibrium coordinate
    ! system [ADD REF]:
    !   - particle density rho
    !   - magnetic shear S
    !   - normal curvature kappa_n
    !   - geodesic curvature kappa_g
    !   - parallel current sigma
    integer function calc_derived_q(grid_eq,eq,met) result(ierr)
        use utilities, only: c
        use num_vars, only: rho_style, use_normalization
        use eq_vars, only: rho_0, vac_perm
        
        character(*), parameter :: rout_name = 'calc_derived_q'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        
        ! local variables
        integer :: kd                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        real(dp), pointer :: D1J(:,:,:)                                         ! D_alpha jac
        real(dp), pointer :: D2J(:,:,:)                                         ! D_psi jac
        real(dp), pointer :: D3J(:,:,:)                                         ! D_theta jac
        ! lower metric factors
        real(dp), pointer :: g13(:,:,:)                                         ! g_alpha,theta
        real(dp), pointer :: D2g13(:,:,:)                                       ! D_psi g_alpha,theta
        real(dp), pointer :: D3g13(:,:,:)                                       ! D_theta g_alpha,theta
        real(dp), pointer :: g23(:,:,:)                                         ! g_psi,theta
        real(dp), pointer :: D1g23(:,:,:)                                       ! D_alpha g_psi,theta
        real(dp), pointer :: D3g23(:,:,:)                                       ! D_theta g_psi,theta
        real(dp), pointer :: g33(:,:,:)                                         ! g_theta,theta
        real(dp), pointer :: D1g33(:,:,:)                                       ! D_alpha g_theta,theta
        real(dp), pointer :: D2g33(:,:,:)                                       ! D_psi g_theta,theta
        real(dp), pointer :: D3g33(:,:,:)                                       ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), pointer :: h12(:,:,:)                                         ! h^alpha,psi
        real(dp), pointer :: D3h12(:,:,:)                                       ! D_theta h^alpha,psi
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        real(dp), pointer :: D3h22(:,:,:)                                       ! D_theta h^psi,psi
        real(dp), pointer :: h23(:,:,:)                                         ! h^psi,theta
        
        ! initialize ierr
        ierr = 0
        
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
        
        ! Calculate particle density rho
        ! choose which density style is being used:
        !   1:  constant, equal to rho_0
        select case (rho_style)
            case (1)                                                            ! HELENA
                eq%rho = rho_0                                                  ! arbitrarily constant (normalized value)
            case default
                err_msg = 'No density style associated with '//&
                    &trim(i2str(rho_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        ! normalize rho if necessary
        if (use_normalization) eq%rho = eq%rho/rho_0
        
        ! Calculate the shear S
        eq%S = -(D3h12/h22 - D3h22*h12/h22**2)/J
        
        ! Calculate the normal curvature kappa_n
        do kd = 1,grid_eq%grp_n_r
            eq%kappa_n(:,:,kd) = J(:,:,kd)**2*eq%pres_FD(kd,1)/g33(:,:,kd) + &
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
        eq%sigma = 1._dp/vac_perm*&
            &(D1g23/J - g23*D1J/J**2 - D2g13/J + g13*D2J/J**2)
        do kd = 1,grid_eq%grp_n_r
            eq%sigma(:,:,kd) = eq%sigma(:,:,kd) - &
                &eq%pres_FD(kd,1)*J(:,:,kd)*g13(:,:,kd)/g33(:,:,kd)
        end do 
    end function calc_derived_q
    
    ! plots the flux quantities in the perturbation grid
    !   safety factor q_saf
    !   rotational transform rot
    !   pressure pres
    !   poloidal flux flux_p
    !   toroidal flux flux_t
    ! [MPI] Only first group
    integer function flux_q_plot(eq,grid_eq) result(ierr)
        use num_vars, only: eq_style, output_style, glb_rank, no_plots
        use grid_ops, only: trim_grid
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'flux_q_plot'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium for this alpha
        type(grid_type), intent(in) :: grid_eq                                  ! normal grid
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_vars = 5                                                   ! nr. of variables to plot
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln), allocatable :: plot_titles(:)                ! plot titles
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: grp_n_r                                                      ! grp_n_r of trimmed grid
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
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set grp_n_r
        grp_n_r = grid_trim%grp_n_r
        
        ! initialize plot titles and file names
        allocate(plot_titles(n_vars))
        plot_titles(1) = 'safety factor []'
        plot_titles(2) = 'rotational transform []'
        plot_titles(3) = 'pressure [pa]'
        plot_titles(4) = 'poloidal flux [Tm^2]'
        plot_titles(5) = 'toroidal flux [Tm^2]'
        
        ! plot the output
        ! choose which output style is being used:
        !   1:  GNUPlot
        !   2:  HDF5
        select case (output_style)
            case (1)                                                        ! GNUPlot
                ierr = flux_q_plot_GP()
                CHCKERR('')
            case (2)                                                        ! HDF5
                ierr = flux_q_plot_HDF5()
                CHCKERR('')
            case default
                err_msg = 'No output style associated with '//&
                    &trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
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
            
            ! allocate variabels if group master
            if (glb_rank.eq.0) then
                allocate(X_plot_2D(grid_trim%n(3),n_vars))
                allocate(Y_plot_2D(grid_trim%n(3),n_vars))
            end if
            
            ! fill the 2D version of the plot
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = get_ser_var(-eq%q_saf_E(1:grp_n_r,0),ser_var_loc)    ! conversion VMEC LH -> RH coord. system
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,1) = ser_var_loc
                    ierr = get_ser_var(-eq%rot_t_E(1:grp_n_r,0),ser_var_loc)    ! conversion VMEC LH -> RH coord. system
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,2) = ser_var_loc
                    ierr = get_ser_var(eq%pres_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,3) = ser_var_loc
                    ierr = get_ser_var(eq%flux_p_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,4) = ser_var_loc
                    ierr = get_ser_var(-eq%flux_t_E(1:grp_n_r,0),ser_var_loc)   ! conversion VMEC LH -> RH coord. system
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,5) = ser_var_loc
                case (2)                                                        ! HELENA
                    ierr = get_ser_var(eq%q_saf_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,1) = ser_var_loc
                    ierr = get_ser_var(eq%rot_t_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,2) = ser_var_loc
                    ierr = get_ser_var(eq%pres_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,3) = ser_var_loc
                    ierr = get_ser_var(eq%flux_p_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,4) = ser_var_loc
                    ierr = get_ser_var(eq%flux_t_E(1:grp_n_r,0),ser_var_loc)
                    CHCKERR('')
                    if (glb_rank.eq.0) Y_plot_2D(:,5) = ser_var_loc
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! rescale if normalized
            if (use_normalization .and. glb_rank.eq.0) then
                Y_plot_2D(:,3) = Y_plot_2D(:,3)*pres_0                          ! pressure
                Y_plot_2D(:,4) = Y_plot_2D(:,4)*psi_0                           ! flux_p
                Y_plot_2D(:,5) = Y_plot_2D(:,5)*psi_0                           ! flux_t
            end if
            
            ! continue the plot if global master
            if (glb_rank.eq.0) then
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
                ! and rot_t, as they are already in plot_jq)
                call writo('The safety factor and rotational transform are not &
                    &plotted here. Instead, use the input variable "plot_jq".')
                call print_GP_2D(plot_titles(3),trim(file_name(1))//'.dat',&
                    &Y_plot_2D(:,3),X_plot_2D(:,3),draw=.false.)
                ! fluxes
                call print_GP_2D(trim(plot_titles(4))//', '//&
                    &trim(plot_titles(5)),trim(file_name(2))//'.dat',&
                    &Y_plot_2D(:,4:5),X_plot_2D(:,4:5),draw=.false.)
                
                ! draw plot
                call draw_GP(plot_titles(3),trim(file_name(1))//'.dat',1,&
                    &.true.,.false.)                                            ! pressure
                call draw_GP(trim(plot_titles(4))//', '//trim(plot_titles(5)),&
                    &trim(file_name(2))//'.dat',2,.true.,.false.)               ! fluxes
                
                ! clean up
                deallocate(X_plot_2D,Y_plot_2D)
            end if
        end function flux_q_plot_GP
        
        ! plots the flux quantities in HDF5
        integer function flux_q_plot_HDF5() result(ierr)
            use output_ops, only: plot_HDF5
            use grid_ops, only: calc_XYZ_grid, trim_grid, extend_grid
            use grid_vars, only: create_grid, dealloc_grid
            
            character(*), parameter :: rout_name = 'flux_q_plot_HDF5'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: X_plot_3D(:,:,:)                           ! x values of 3D plot
            real(dp), allocatable :: Y_plot_3D(:,:,:)                           ! y values of 3D plot
            real(dp), allocatable :: Z_plot_3D(:,:,:)                           ! z values of 3D plot
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
            allocate(Y_plot_2D(grid_trim%grp_n_r,n_vars))
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    Y_plot_2D(:,1) = -eq%q_saf_E(1:grp_n_r,0)                   ! conversion VMEC LH -> RH coord. system
                    Y_plot_2D(:,2) = -eq%rot_t_E(1:grp_n_r,0)                   ! conversion VMEC LH -> RH coord. system
                    Y_plot_2D(:,3) = eq%pres_E(1:grp_n_r,0)
                    Y_plot_2D(:,4) = eq%flux_p_E(1:grp_n_r,0)
                    Y_plot_2D(:,5) = -eq%flux_t_E(1:grp_n_r,0)                  ! conversion VMEC LH -> RH coord. system
                case (2)                                                        ! HELENA
                    Y_plot_2D(:,1) = eq%q_saf_E(1:grp_n_r,0)
                    Y_plot_2D(:,2) = eq%rot_t_E(1:grp_n_r,0)
                    Y_plot_2D(:,3) = eq%pres_E(1:grp_n_r,0)
                    Y_plot_2D(:,4) = eq%flux_p_E(1:grp_n_r,0)
                    Y_plot_2D(:,5) = eq%flux_t_E(1:grp_n_r,0)
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! extend trimmed equilibrium grid
            ierr = extend_grid(grid_trim,grid_plot)
            CHCKERR('')
            
            ! calculate 3D X,Y and Z
            ierr = calc_XYZ_grid(grid_plot,X_plot_3D,Y_plot_3D,Z_plot_3D)
            CHCKERR('')
            
            ! set up plot_dim and plot_offset
            plot_dim = [grid_plot%n(1),grid_plot%n(2),grid_plot%n(3),n_vars]
            plot_offset = [0,0,grid_plot%i_min-1,n_vars]
            
            ! set up total plot variables
            allocate(X_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%grp_n_r,&
                &n_vars))
            allocate(Y_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%grp_n_r,&
                &n_vars))
            allocate(Z_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%grp_n_r,&
                &n_vars))
            allocate(f_plot(grid_plot%n(1),grid_plot%n(2),grid_plot%grp_n_r,&
                &n_vars))
            do id = 1,n_vars
                X_plot(:,:,:,id) = X_plot_3D
                Y_plot(:,:,:,id) = Y_plot_3D
                Z_plot(:,:,:,id) = Z_plot_3D
            end do
            do kd = 1,grid_plot%grp_n_r
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
            deallocate(X_plot_3D,Y_plot_3D,Z_plot_3D)
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
        use num_vars, only: glb_rank, eq_style, mu_0_original, use_normalization
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
            if (glb_rank.eq.0) then
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
        use num_vars, only: use_normalization, eq_style, mu_0_original, glb_rank
        use VMEC, only: normalize_VMEC
        use eq_vars, only: vac_perm
        
        character(*), parameter :: rout_name = 'normalize_input'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only normalize if needed
        if (use_normalization .and. glb_rank.eq.0) then
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
        use num_vars, only: output_style
        
        character(*), parameter :: rout_name = 'print_output_eq'
        
        ! input / output
        type(grid_type) :: grid_eq                                              ! equilibrium grid variables
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Writing equilibrium variables to output file')
        call lvl_ud(1)
        
        ! print according to output_style
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                call writo('WARNING: No possibility to save equilibrium &
                    &results for output style '//trim(i2str(output_style))//&
                    &' implemented')
            case(2)                                                             ! HDF5 output
                ierr = print_output_eq_HDF5(grid_eq,eq,met)
                CHCKERR('')
            case default
                err_msg = 'No style associated with '//&
                    &trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! user output
        call lvl_ud(-1)
        call writo('Equilibrium variables written to output')
    contains
        ! HDF5 version
        integer function print_output_eq_HDF5(grid_eq,eq,met) result(ierr)
            use num_vars, only: eq_style, max_deriv
            use HDF5_ops, only: print_HDF5_arrs, &
                &var_1D
            use HELENA, only: R_H, Z_H, nchi
            use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, ntor
            use grid_ops, only: trim_grid
            
            character(*), parameter :: rout_name = 'print_output_eq_HDF5'
            
            ! input / output
            type(grid_type) :: grid_eq                                          ! equilibrium grid variables
            type(eq_type) :: eq                                                 ! equilibrium variables
            type(met_type) :: met                                               ! metric variables
            
            ! local variables
            type(var_1D), pointer :: eq_1D(:)                                   ! 1D equivalent of eq. variables
            type(var_1D), pointer :: eq_1D_loc                                  ! local element in eq_1D
            character(len=max_str_ln) :: err_msg                                ! error message
            type(grid_type) :: grid_trim                                        ! trimmed grid
            integer :: i_min, i_max                                             ! min. and max. index of variables
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Preparing variables for writing')
            call lvl_ud(1)
            
            ! trim grid
            ierr = trim_grid(grid_eq,grid_trim)
            CHCKERR('')
            
            ! set i_min and i_max
            i_min = 1
            i_max = grid_trim%grp_n_r
            
            ! Set up the 1D equivalents  of the equilibrium variables, with size
            ! depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    allocate(eq_1D(22))
                case (2)                                                        ! HELENA
                    allocate(eq_1D(18))
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! Set up common variables
            
            ! 1. r_F
            eq_1D_loc => eq_1D(1)
            eq_1D_loc%var_name = 'r_F'
            allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
            allocate(eq_1D_loc%grp_i_min(1),eq_1D_loc%grp_i_max(1))
            eq_1D_loc%tot_i_min = [1]
            eq_1D_loc%tot_i_max = [grid_trim%n(3)]
            eq_1D_loc%grp_i_min = [grid_trim%i_min]
            eq_1D_loc%grp_i_max = [grid_trim%i_max]
            allocate(eq_1D_loc%p(size(grid_trim%r_F)))
            eq_1D_loc%p = grid_trim%r_F(i_min:i_max)
            
            ! 2. theta_F
            eq_1D_loc => eq_1D(2)
            eq_1D_loc%var_name = 'theta_F'
            allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
            allocate(eq_1D_loc%grp_i_min(3),eq_1D_loc%grp_i_max(3))
            eq_1D_loc%tot_i_min = [1,1,1]
            eq_1D_loc%tot_i_max = grid_trim%n
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min]
            eq_1D_loc%grp_i_max = &
                &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
            allocate(eq_1D_loc%p(size(grid_trim%theta_F(:,:,i_min:i_max))))
            eq_1D_loc%p = reshape(grid_trim%theta_F(:,:,i_min:i_max),&
                &[size(grid_trim%theta_F(:,:,i_min:i_max))])
            
            ! 3. zeta_F
            eq_1D_loc => eq_1D(3)
            eq_1D_loc%var_name = 'zeta_F'
            allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
            allocate(eq_1D_loc%grp_i_min(3),eq_1D_loc%grp_i_max(3))
            eq_1D_loc%tot_i_min = [1,1,1]
            eq_1D_loc%tot_i_max = grid_trim%n
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min]
            eq_1D_loc%grp_i_max = &
                &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
            allocate(eq_1D_loc%p(size(grid_trim%zeta_F(:,:,i_min:i_max))))
            eq_1D_loc%p = reshape(grid_trim%zeta_F(:,:,i_min:i_max),&
                &[size(grid_trim%zeta_F(:,:,i_min:i_max))])
            
            ! 4. pres_FD
            eq_1D_loc => eq_1D(4)
            eq_1D_loc%var_name = 'pres_FD'
            allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
            allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
            eq_1D_loc%tot_i_min = [1,0]
            eq_1D_loc%tot_i_max = [grid_trim%n(3),max_deriv]
            eq_1D_loc%grp_i_min = [grid_trim%i_min,0]
            eq_1D_loc%grp_i_max = [grid_trim%i_max,max_deriv]
            allocate(eq_1D_loc%p(size(eq%pres_FD(i_min:i_max,:))))
            eq_1D_loc%p = reshape(eq%pres_FD(i_min:i_max,:),&
                &[size(eq%pres_FD(i_min:i_max,:))])
            
            ! 5. q_saf_FD
            eq_1D_loc => eq_1D(5)
            eq_1D_loc%var_name = 'q_saf_FD'
            allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
            allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
            eq_1D_loc%tot_i_min = [1,0]
            eq_1D_loc%tot_i_max = [grid_trim%n(3),max_deriv]
            eq_1D_loc%grp_i_min = [grid_trim%i_min,0]
            eq_1D_loc%grp_i_max = [grid_trim%i_max,max_deriv]
            allocate(eq_1D_loc%p(size(eq%q_saf_FD(i_min:i_max,:))))
            eq_1D_loc%p = reshape(eq%q_saf_FD(i_min:i_max,:),&
                &[size(eq%q_saf_FD(i_min:i_max,:))])
            
            ! 6. rot_t_FD
            eq_1D_loc => eq_1D(6)
            eq_1D_loc%var_name = 'rot_t_FD'
            allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
            allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
            eq_1D_loc%tot_i_min = [1,0]
            eq_1D_loc%tot_i_max = [grid_trim%n(3),max_deriv]
            eq_1D_loc%grp_i_min = [grid_trim%i_min,0]
            eq_1D_loc%grp_i_max = [grid_trim%i_max,max_deriv]
            allocate(eq_1D_loc%p(size(eq%rot_t_FD(i_min:i_max,:))))
            eq_1D_loc%p = reshape(eq%rot_t_FD(i_min:i_max,:),&
                &[size(eq%rot_t_FD(i_min:i_max,:))])
            
            ! 7. flux_p_FD
            eq_1D_loc => eq_1D(7)
            eq_1D_loc%var_name = 'flux_p_FD'
            allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
            allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
            eq_1D_loc%tot_i_min = [1,0]
            eq_1D_loc%tot_i_max = [grid_trim%n(3),max_deriv]
            eq_1D_loc%grp_i_min = [grid_trim%i_min,0]
            eq_1D_loc%grp_i_max = [grid_trim%i_max,max_deriv]
            allocate(eq_1D_loc%p(size(eq%flux_p_FD(i_min:i_max,:))))
            eq_1D_loc%p = reshape(eq%flux_p_FD(i_min:i_max,:),&
                &[size(eq%flux_p_FD(i_min:i_max,:))])
            
            ! 8. flux_t_FD
            eq_1D_loc => eq_1D(8)
            eq_1D_loc%var_name = 'flux_t_FD'
            allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
            allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
            eq_1D_loc%tot_i_min = [1,0]
            eq_1D_loc%tot_i_max = [grid_trim%n(3),max_deriv]
            eq_1D_loc%grp_i_min = [grid_trim%i_min,0]
            eq_1D_loc%grp_i_max = [grid_trim%i_max,max_deriv]
            allocate(eq_1D_loc%p(size(eq%flux_t_FD(i_min:i_max,:))))
            eq_1D_loc%p = reshape(eq%flux_t_FD(i_min:i_max,:),&
                &[size(eq%flux_t_FD(i_min:i_max,:))])
            
            ! 9. rho
            eq_1D_loc => eq_1D(9)
            eq_1D_loc%var_name = 'rho'
            allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
            allocate(eq_1D_loc%grp_i_min(1),eq_1D_loc%grp_i_max(1))
            eq_1D_loc%tot_i_min = 1
            eq_1D_loc%tot_i_max = grid_trim%n(3)
            eq_1D_loc%grp_i_min = grid_trim%i_min
            eq_1D_loc%grp_i_max = grid_trim%i_max
            allocate(eq_1D_loc%p(size(eq%rho(i_min:i_max))))
            eq_1D_loc%p = eq%rho(i_min:i_max)
            
            ! 10. S
            eq_1D_loc => eq_1D(10)
            eq_1D_loc%var_name = 'S'
            allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
            allocate(eq_1D_loc%grp_i_min(3),eq_1D_loc%grp_i_max(3))
            eq_1D_loc%tot_i_min = [1,1,1]
            eq_1D_loc%tot_i_max = grid_trim%n
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min]
            eq_1D_loc%grp_i_max = &
                &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
            allocate(eq_1D_loc%p(size(eq%S(:,:,i_min:i_max))))
            eq_1D_loc%p = reshape(eq%S(:,:,i_min:i_max),&
                &[size(eq%S(:,:,i_min:i_max))])
            
            ! 11. kappa_n
            eq_1D_loc => eq_1D(11)
            eq_1D_loc%var_name = 'kappa_n'
            allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
            allocate(eq_1D_loc%grp_i_min(3),eq_1D_loc%grp_i_max(3))
            eq_1D_loc%tot_i_min = [1,1,1]
            eq_1D_loc%tot_i_max = grid_trim%n
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min]
            eq_1D_loc%grp_i_max = &
                &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
            allocate(eq_1D_loc%p(size(eq%kappa_n(:,:,i_min:i_max))))
            eq_1D_loc%p = reshape(eq%kappa_n(:,:,i_min:i_max),&
                &[size(eq%kappa_n(:,:,i_min:i_max))])
            
            ! 12. kappa_g
            eq_1D_loc => eq_1D(12)
            eq_1D_loc%var_name = 'kappa_g'
            allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
            allocate(eq_1D_loc%grp_i_min(3),eq_1D_loc%grp_i_max(3))
            eq_1D_loc%tot_i_min = [1,1,1]
            eq_1D_loc%tot_i_max = grid_trim%n
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min]
            eq_1D_loc%grp_i_max = &
                &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
            allocate(eq_1D_loc%p(size(eq%kappa_g(:,:,i_min:i_max))))
            eq_1D_loc%p = reshape(eq%kappa_g(:,:,i_min:i_max),&
                &[size(eq%kappa_g(:,:,i_min:i_max))])
            
            ! 13. sigma
            eq_1D_loc => eq_1D(13)
            eq_1D_loc%var_name = 'sigma'
            allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
            allocate(eq_1D_loc%grp_i_min(3),eq_1D_loc%grp_i_max(3))
            eq_1D_loc%tot_i_min = [1,1,1]
            eq_1D_loc%tot_i_max = grid_trim%n
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min]
            eq_1D_loc%grp_i_max = &
                &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
            allocate(eq_1D_loc%p(size(eq%sigma(:,:,i_min:i_max))))
            eq_1D_loc%p = reshape(eq%sigma(:,:,i_min:i_max),&
                &[size(eq%sigma(:,:,i_min:i_max))])
            
            ! 14. g_FD
            eq_1D_loc => eq_1D(14)
            eq_1D_loc%var_name = 'g_FD'
            allocate(eq_1D_loc%tot_i_min(7),eq_1D_loc%tot_i_max(7))
            allocate(eq_1D_loc%grp_i_min(7),eq_1D_loc%grp_i_max(7))
            eq_1D_loc%tot_i_min = [1,1,1,1,0,0,0]
            eq_1D_loc%tot_i_max = [grid_trim%n,6,max_deriv,max_deriv,max_deriv]
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min,1,0,0,0]
            eq_1D_loc%grp_i_max = [grid_trim%n(1),grid_trim%n(2),&
                &grid_trim%i_max,6,max_deriv,max_deriv,max_deriv]
            allocate(eq_1D_loc%p(size(met%g_FD(:,:,i_min:i_max,:,:,:,:))))
            eq_1D_loc%p = reshape(met%g_FD(:,:,i_min:i_max,:,:,:,:),&
                &[size(met%g_FD(:,:,i_min:i_max,:,:,:,:))])
            
            ! 15. h_FD
            eq_1D_loc => eq_1D(15)
            eq_1D_loc%var_name = 'h_FD'
            allocate(eq_1D_loc%tot_i_min(7),eq_1D_loc%tot_i_max(7))
            allocate(eq_1D_loc%grp_i_min(7),eq_1D_loc%grp_i_max(7))
            eq_1D_loc%tot_i_min = [1,1,1,1,0,0,0]
            eq_1D_loc%tot_i_max = [grid_trim%n,6,max_deriv,max_deriv,max_deriv]
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min,1,0,0,0]
            eq_1D_loc%grp_i_max = [grid_trim%n(1),grid_trim%n(2),&
                &grid_trim%i_max,6,max_deriv,max_deriv,max_deriv]
            allocate(eq_1D_loc%p(size(met%h_FD(:,:,i_min:i_max,:,:,:,:))))
            eq_1D_loc%p = reshape(met%h_FD(:,:,i_min:i_max,:,:,:,:),&
                &[size(met%h_FD(:,:,i_min:i_max,:,:,:,:))])
            
            ! 16. jac_FD
            eq_1D_loc => eq_1D(16)
            eq_1D_loc%var_name = 'jac_FD'
            allocate(eq_1D_loc%tot_i_min(6),eq_1D_loc%tot_i_max(6))
            allocate(eq_1D_loc%grp_i_min(6),eq_1D_loc%grp_i_max(6))
            eq_1D_loc%tot_i_min = [1,1,1,0,0,0]
            eq_1D_loc%tot_i_max = [grid_trim%n,max_deriv,max_deriv,max_deriv]
            eq_1D_loc%grp_i_min = [1,1,grid_trim%i_min,0,0,0]
            eq_1D_loc%grp_i_max = [grid_trim%n(1),grid_trim%n(2),&
                &grid_trim%i_max,max_deriv,max_deriv,max_deriv]
            allocate(eq_1D_loc%p(size(met%jac_FD(:,:,i_min:i_max,:,:,:))))
            eq_1D_loc%p = reshape(met%jac_FD(:,:,i_min:i_max,:,:,:),&
                &[size(met%jac_FD(:,:,i_min:i_max,:,:,:))])
            
            ! Set up particular variables, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! 17. R_V_c
                    eq_1D_loc => eq_1D(17)
                    eq_1D_loc%var_name = 'R_V_c'
                    allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                    allocate(eq_1D_loc%grp_i_min(4),eq_1D_loc%grp_i_max(4))
                    eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                    eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),max_deriv]
                    eq_1D_loc%grp_i_min = [0,-ntor,grid_trim%i_min,0]
                    eq_1D_loc%grp_i_max = &
                        &[mpol-1,ntor,grid_trim%i_max,max_deriv]
                    allocate(eq_1D_loc%p(size(R_V_c(:,:,i_min:i_max,:))))
                    eq_1D_loc%p = reshape(R_V_c(:,:,i_min:i_max,:),&
                        &[size(R_V_c(:,:,i_min:i_max,:))])
                    
                    ! 18. R_V_s
                    eq_1D_loc => eq_1D(18)
                    eq_1D_loc%var_name = 'R_V_s'
                    allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                    allocate(eq_1D_loc%grp_i_min(4),eq_1D_loc%grp_i_max(4))
                    eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                    eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),max_deriv]
                    eq_1D_loc%grp_i_min = [0,-ntor,grid_trim%i_min,0]
                    eq_1D_loc%grp_i_max = &
                        &[mpol-1,ntor,grid_trim%i_max,max_deriv]
                    allocate(eq_1D_loc%p(size(R_V_s(:,:,i_min:i_max,:))))
                    eq_1D_loc%p = reshape(R_V_s(:,:,i_min:i_max,:),&
                        &[size(R_V_s(:,:,i_min:i_max,:))])
                    
                    ! 19. Z_V_c
                    eq_1D_loc => eq_1D(19)
                    eq_1D_loc%var_name = 'Z_V_c'
                    allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                    allocate(eq_1D_loc%grp_i_min(4),eq_1D_loc%grp_i_max(4))
                    eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                    eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),max_deriv]
                    eq_1D_loc%grp_i_min = [0,-ntor,grid_trim%i_min,0]
                    eq_1D_loc%grp_i_max = &
                        &[mpol-1,ntor,grid_trim%i_max,max_deriv]
                    allocate(eq_1D_loc%p(size(Z_V_c(:,:,i_min:i_max,:))))
                    eq_1D_loc%p = reshape(Z_V_c(:,:,i_min:i_max,:),&
                        &[size(Z_V_c(:,:,i_min:i_max,:))])
                    
                    ! 20. Z_V_s
                    eq_1D_loc => eq_1D(20)
                    eq_1D_loc%var_name = 'Z_V_s'
                    allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                    allocate(eq_1D_loc%grp_i_min(4),eq_1D_loc%grp_i_max(4))
                    eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                    eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),max_deriv]
                    eq_1D_loc%grp_i_min = [0,-ntor,grid_trim%i_min,0]
                    eq_1D_loc%grp_i_max = &
                        &[mpol-1,ntor,grid_trim%i_max,max_deriv]
                    allocate(eq_1D_loc%p(size(Z_V_s(:,:,i_min:i_max,:))))
                    eq_1D_loc%p = reshape(Z_V_s(:,:,i_min:i_max,:),&
                        &[size(Z_V_s(:,:,i_min:i_max,:))])
                    
                    ! 21. L_V_c
                    eq_1D_loc => eq_1D(21)
                    eq_1D_loc%var_name = 'L_V_c'
                    allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                    allocate(eq_1D_loc%grp_i_min(4),eq_1D_loc%grp_i_max(4))
                    eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                    eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),max_deriv]
                    eq_1D_loc%grp_i_min = [0,-ntor,grid_trim%i_min,0]
                    eq_1D_loc%grp_i_max = &
                        &[mpol-1,ntor,grid_trim%i_max,max_deriv]
                    allocate(eq_1D_loc%p(size(L_V_c(:,:,i_min:i_max,:))))
                    eq_1D_loc%p = reshape(L_V_c(:,:,i_min:i_max,:),&
                        &[size(L_V_c(:,:,i_min:i_max,:))])
                    
                    ! 22. L_V_s
                    eq_1D_loc => eq_1D(22)
                    eq_1D_loc%var_name = 'L_V_s'
                    allocate(eq_1D_loc%tot_i_min(4),eq_1D_loc%tot_i_max(4))
                    allocate(eq_1D_loc%grp_i_min(4),eq_1D_loc%grp_i_max(4))
                    eq_1D_loc%tot_i_min = [0,-ntor,1,0]
                    eq_1D_loc%tot_i_max = [mpol-1,ntor,grid_trim%n(3),max_deriv]
                    eq_1D_loc%grp_i_min = [0,-ntor,grid_trim%i_min,0]
                    eq_1D_loc%grp_i_max = &
                        &[mpol-1,ntor,grid_trim%i_max,max_deriv]
                    allocate(eq_1D_loc%p(size(L_V_s(:,:,i_min:i_max,:))))
                    eq_1D_loc%p = reshape(L_V_s(:,:,i_min:i_max,:),&
                        &[size(L_V_s(:,:,i_min:i_max,:))])
                case (2)                                                        ! HELENA
                    ! 17. R_H
                    eq_1D_loc => eq_1D(17)
                    eq_1D_loc%var_name = 'R_H'
                    allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                    allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
                    eq_1D_loc%tot_i_min = [1,1]
                    eq_1D_loc%tot_i_max = [nchi,grid_trim%n(3)]
                    eq_1D_loc%grp_i_min = [1,grid_trim%i_min]
                    eq_1D_loc%grp_i_max = [nchi,grid_trim%i_max]
                    allocate(eq_1D_loc%p(size(R_H(:,i_min:i_max))))
                    eq_1D_loc%p = reshape(R_H(:,i_min:i_max),&
                        &[size(R_H(:,i_min:i_max))])
                    
                    ! 18. Z_H
                    eq_1D_loc => eq_1D(18)
                    eq_1D_loc%var_name = 'Z_H'
                    allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
                    allocate(eq_1D_loc%grp_i_min(2),eq_1D_loc%grp_i_max(2))
                    eq_1D_loc%tot_i_min = [1,1]
                    eq_1D_loc%tot_i_max = [nchi,grid_trim%n(3)]
                    eq_1D_loc%grp_i_min = [1,grid_trim%i_min]
                    eq_1D_loc%grp_i_max = [nchi,grid_trim%i_max]
                    allocate(eq_1D_loc%p(size(Z_H(:,i_min:i_max))))
                    eq_1D_loc%p = reshape(Z_H(:,i_min:i_max),&
                        &[size(Z_H(:,i_min:i_max))])
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Writing using HDF5')
            call lvl_ud(1)
            
            ! write
            ierr = print_HDF5_arrs(eq_1D)
            CHCKERR('')
            
            ! user output
            call lvl_ud(-1)
            call writo('Writing complete')
        end function print_output_eq_HDF5
    end function print_output_eq
end module eq_ops
