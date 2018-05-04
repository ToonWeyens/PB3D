!------------------------------------------------------------------------------!
!> Operations on the equilibrium variables.
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln, max_deriv
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
#if ldebug
    use num_utilities, only: check_deriv
#endif
    
    implicit none
    private
    public calc_eq, calc_derived_q, calc_normalization_const, normalize_input, &
        &print_output_eq, flux_q_plot, redistribute_output_eq, divide_eq_jobs, &
        &calc_eq_jobs_lims, calc_T_HF, B_plot, J_plot, kappa_plot, delta_r_plot
#if ldebug
    public debug_calc_derived_q, debug_create_VMEC_input, debug_J_plot
#endif
    
    ! global variables
    integer :: fund_n_par                                                       !< fundamental interval width
    logical :: BR_normalization_provided(2)                                     ! used to export HELENA to VMEC
#if ldebug
    !> \ldebug
    logical :: debug_calc_derived_q = .true.                                   !< plot debug information for calc_derived_q()
    !> \ldebug
    logical :: debug_J_plot = .false.                                           !< plot debug information for j_plot()
    !> \ldebug
    logical :: debug_create_VMEC_input = .false.                                !< plot debug information for create_vmec_input()
#endif
    
    ! interfaces
    
    !> \public  Calculate the  equilibrium quantities  on a  grid determined  by
    !! straight field lines.
    !!
    !! This grid has the dimensions \c (n_par,loc_n_r).
    !!
    !! Optionally, for  \c eq_2, the  used variables  can be deallocated  on the
    !! fly, to limit memory usage.
    !!
    !! \return ierr
    interface calc_eq
        !> \public
        module procedure calc_eq_1
        !> \public
        module procedure calc_eq_2
    end interface
    
    !> \public Print equilibrium quantities to an output file:
    !!
    !!  - flux:
    !!      - \c pres_FD,
    !!      - \c q_saf_FD,
    !!      - \c rot_t_FD,
    !!      - \c flux_p_FD,
    !!      - \c flux_t_FD,
    !!      - \c rho,
    !!      - \c S,
    !!      - \c kappa_n,
    !!      - \c kappa_g,
    !!      - \c sigma
    !!  - metric:
    !!      - \c g_FD,
    !!      - \c h_FD,
    !!      - \c jac_FD
    !!
    !! If \c  rich_lvl is provided, \c  "_R_[rich_lvl]" is appended to  the data
    !! name if it is \c >0 (only for \c eq_2).
    !! 
    !! Optionally,  for \c  eq_2, it  can be  specified that  this is  a divided
    !! parallel grid, corresponding  to the variable \c  eq_jobs_lims with index
    !! \c eq_job_nr. In  this case, the total  grid size is adjusted  to the one
    !! specified by \c eq_jobs_lims and the grid is written as a subset.
    !! \note 
    !!  -# The equilibrium quantities are outputted in Flux coordinates.
    !!  -#  The metric  equilibrium quantities  can be  deallocated on  the fly,
    !!  which is useful  if this routine is followed by  a deallocation any way,
    !!  so that memory usage does not almost double.
    !!  -# \c print_output_eq_2  is only used by  HELENA now, as for  VMEC it is
    !!  too slow since there are often multiple VMEC equilibrium jobs, while for
    !!  HELENA this is explicitely forbidden.
    !!
    !! \return ierr
    interface print_output_eq
        !> \public
        module procedure print_output_eq_1
        !> \public
        module procedure print_output_eq_2
    end interface
    
    !> \public  Redistribute  the  equilibrium  variables,  but  only  the  Flux
    !! variables are saved.
    !!
    !! \see redistribute_output_grid()
    !!
    !! \return ierr
    interface redistribute_output_eq
        !> \public
        module procedure redistribute_output_eq_1
        !> \public
        module procedure redistribute_output_eq_2
    end interface

    !> \public Calculate  \f$R\f$, \f$Z\f$  \& \f$\lambda\f$ and  derivatives in
    !! VMEC coordinates.
    !!
    !! This is done at the grid points  given by the variables \c theta_E and \c
    !! zeta_E (contained in \c trigon_factors) and at every normal point.
    !!
    !! The  derivatives are  indicated  by the  variable \c  deriv  which has  3
    !! indices
    !!
    !! \note  There is  no  HELENA equivalent  because  for HELENA  simulations,
    !! \f$R\f$  and \f$Z\f$  are not  necessary  for calculation  of the  metric
    !! coefficients, and \f$\lambda\f$ does not exist.
    !!
    !! \see calc_trigon_factors()
    !!
    !! \return ierr
    interface calc_RZL
        !> \public
        module procedure calc_RZL_ind
        !> \public
        module procedure calc_RZL_arr
    end interface
    
    !> \public  Calculate  the  lower   metric  elements  in  the  C(ylindrical)
    !! coordinate system.
    !!
    !! This is done directly using the formula's in \cite Weyens3D
    !!
    !! \return ierr
    interface calc_g_C
        !> \public
        module procedure calc_g_C_ind
        !> \public
        module procedure calc_g_C_arr
    end interface
    
    !> \public Calculate the lower metric coefficients in the equilibrium V(MEC)
    !! coordinate system.
    !! 
    !! This  is  done  using  the   metric  coefficients  in  the  C(ylindrical)
    !! coordinate system and the transformation matrices
    !!
    !! \note It is assumed that the lower order derivatives have been calculated
    !! already. If not, the results will be incorrect.
    !!
    !! \return ierr
    interface calc_g_V
        !> \public
        module procedure calc_g_V_ind
        !> \public
        module procedure calc_g_V_arr
    end interface
    
    !> \public  Calculate  the  lower  metric coefficients  in  the  equilibrium
    !! H(ELENA) coordinate system.
    !!
    !! This is done using the HELENA output
    !!
    !! \return ierr
    interface calc_g_H
        !> \public
        module procedure calc_g_H_ind
        !> \public
        module procedure calc_g_H_arr
    end interface
    
    !> \public  Calculate  the  metric  coefficients in  the  F(lux)  coordinate
    !! system.
    !!
    !! This is done using the  metric coefficients in the equilibrium coordinate
    !! system and the transformation matrices.
    !!
    !! \return ierr
    interface calc_g_F
        !> \public
        module procedure calc_g_F_ind
        !> \public
        module procedure calc_g_F_arr
    end interface
    
    !> \public  Calculate  \f$\mathcal{J}_\text{V}\f$,  the jacobian  in  V(MEC)
    !! coordinates.
    !!
    !! This is done using VMEC output.
    !!
    !! \return ierr
    interface calc_jac_V
        !> \public
        module procedure calc_jac_V_ind
        !> \public
        module procedure calc_jac_V_arr
    end interface
    
    !> \public  Calculate  \f$\mathcal{J}_\text{H}\f$,  the jacobian  in  HELENA
    !! coordinates.
    !!
    !! This is done directly from
    !!  \f[\mathcal{J}_\text{H} = q \frac{R^2}{F} \f]
    !!
    !! \note It is assumed that the lower order derivatives have been calculated
    !! already. If not, the results will be incorrect.
    !!
    !! \return ierr
    interface calc_jac_H
        !> \public
        module procedure calc_jac_H_ind
        !> \public
        module procedure calc_jac_H_arr
    end interface
    
    !> \public  Calculate  \f$\mathcal{J}_\text{F}\f$,   the  jacobian  in  Flux
    !! coordinates.
    !! 
    !! This is done directly from 
    !!  \f[ \mathcal{J}_\text{F} =
    !!      \text{det}\left(\overline{\text{T}}_\text{F}^\text{E}\right)
    !!      \mathcal{J}_\text{E} \f]
    !!
    !! \note It is assumed that the lower order derivatives have been calculated
    !! already. If not, the results will be incorrect.
    !!
    !! \return ierr
    interface calc_jac_F
        !> \public
        module procedure calc_jac_F_ind
        !> \public
        module procedure calc_jac_F_arr
    end interface
    
    !> \public   Calculate    \f$\overline{\text{T}}_\text{C}^\text{V}\f$,   the
    !! transformation matrix between C(ylindrical) and V(mec) coordinate system.
    !!
    !! This is done directly using the formula's in \cite Weyens3D
    !!
    !! \return ierr
    interface calc_T_VC
        !> \public
        module procedure calc_T_VC_ind
        !> \public
        module procedure calc_T_VC_arr
    end interface
    
    !> \public   Calculate    \f$\overline{\text{T}}_\text{V}^\text{F}\f$,   the
    !! transformation matrix between V(MEC) and F(lux) coordinate systems.
    !!
    !! This is done directly using the formula's in \cite Weyens3D
    !!
    !! \return ierr
    interface calc_T_VF
        !> \public
        module procedure calc_T_VF_ind
        !> \public
        module procedure calc_T_VF_arr
    end interface
    
    !> \public   Calculate    \f$\overline{\text{T}}_\text{H}^\text{F}\f$,   the
    !! transformation matrix between H(ELENA) and F(lux) coordinate systems.
    !!
    !! This is done directly using the formula's in \cite Weyens3D
    !!
    !! \return ierr
    interface calc_T_HF
        !> \public
        module procedure calc_T_HF_ind
        !> \public
        module procedure calc_T_HF_arr
    end interface
    
contains
    !> \private flux version
    integer function calc_eq_1(grid_eq,eq) result(ierr)
        use num_vars, only: eq_style, rho_style, use_normalization, &
            &use_pol_flux_E, use_pol_flux_F
        use eq_vars, only: rho_0
        use num_utilities, only: derivs
        use eq_utilities, only: calc_F_derivs
        
        character(*), parameter :: rout_name = 'calc_eq_1'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        type(eq_1_type), intent(inout) :: eq                                    !< flux equilibrium variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Start setting up flux equilibrium quantities')
        
        call lvl_ud(1)
        
        ! create equilibrium
        call eq%init(grid_eq)
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call calc_flux_q_VMEC()
            case (2)                                                            ! HELENA
                call calc_flux_q_HEL()
        end select
        
        ! take local variables
        grid_eq%loc_r_E = grid_eq%r_E(grid_eq%i_min:grid_eq%i_max)
        grid_eq%loc_r_F = grid_eq%r_F(grid_eq%i_min:grid_eq%i_max)
        
        ! Transform flux equilibrium E into F derivatives
        ierr = calc_F_derivs(grid_eq,eq)
        CHCKERR('')
        
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
        
        call lvl_ud(-1)
        
        call writo('Done setting up flux equilibrium quantities')
    contains
        ! VMEC version
        ! The VMEC normal coord. is  the toroidal (or poloidal) flux, normalized
        ! wrt.  to  the  maximum  flux,  equidistantly,  so  the  step  size  is
        ! 1/(n(3)-1)
        !> \private
        subroutine calc_flux_q_VMEC()
            use VMEC_vars, only: rot_t_V, q_saf_V, flux_t_V, flux_p_V, pres_V
            use eq_vars, only: max_flux_E
            
            ! copy flux variables
            eq%flux_p_E = flux_p_V(grid_eq%i_min:grid_eq%i_max,:)
            eq%flux_t_E = flux_t_V(grid_eq%i_min:grid_eq%i_max,:)
            eq%pres_E = pres_V(grid_eq%i_min:grid_eq%i_max,:)
            eq%q_saf_E = q_saf_V(grid_eq%i_min:grid_eq%i_max,:)
            eq%rot_t_E = rot_t_V(grid_eq%i_min:grid_eq%i_max,:)
            
            ! max flux and  normal coord. of eq grid  in Equilibrium coordinates
            ! (uses poloidal flux by default)
            if (use_pol_flux_E) then
                grid_eq%r_E = flux_p_V(:,0)/max_flux_E
            else
                grid_eq%r_E = flux_t_V(:,0)/max_flux_E
            end if
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_V(:,0)/(2*pi)                              ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = - flux_t_V(:,0)/(2*pi)                            ! psi_F = flux_t/2pi, conversion VMEC LH -> PB3D RH
            end if
        end subroutine calc_flux_q_VMEC
        
        ! HELENA version
        ! The HELENA normal coord. is the poloidal flux divided by 2pi
        !> \private
        subroutine calc_flux_q_HEL()
            use HELENA_vars, only: q_saf_H, rot_t_H, flux_p_H, flux_t_H, pres_H
            use num_utilities, only: calc_int
            
            ! copy flux variables
            eq%flux_p_E = flux_p_H(grid_eq%i_min:grid_eq%i_max,:)
            eq%flux_t_E = flux_t_H(grid_eq%i_min:grid_eq%i_max,:)
            eq%pres_E = pres_H(grid_eq%i_min:grid_eq%i_max,:)
            eq%q_saf_E = q_saf_H(grid_eq%i_min:grid_eq%i_max,:)
            eq%rot_t_E = rot_t_H(grid_eq%i_min:grid_eq%i_max,:)
            
            ! max flux and  normal coord. of eq grid  in Equilibrium coordinates
            ! (uses poloidal flux by default)
            grid_eq%r_E = flux_p_H(:,0)/(2*pi)
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_H(:,0)/(2*pi)                              ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = flux_t_H(:,0)/(2*pi)                              ! psi_F = flux_t/2pi
            end if
        end subroutine calc_flux_q_HEL
    end function calc_eq_1
    !> \private metric version
    integer function calc_eq_2(grid_eq,eq_1,eq_2,dealloc_vars) result(ierr)
        use num_vars, only: eq_style, export_HEL, tol_zero
        use num_utilities, only: derivs, c
        use eq_utilities, only: calc_inv_met, calc_F_derivs
        use VMEC_utilities, only: calc_trigon_factors
#if ldebug
        use num_vars, only: ltest
        use input_utilities, only: get_log, pause_prog
        use HELENA_ops, only: test_metrics_H, test_harm_cont_H
        use HELENA_vars, only: chi_H
#endif
        
        character(*), parameter :: rout_name = 'calc_eq_2'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< metric equilibrium variables
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium variables
        logical, intent(in), optional :: dealloc_vars                           !< deallocate variables on the fly after writing
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: pmone                                                        ! plus or minus one
        logical :: dealloc_vars_loc                                             ! local dealloc_vars
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Start setting up metric equilibrium quantities')
        
        call lvl_ud(1)
        
        ! set up local dealloc_vars
        dealloc_vars_loc = .false.
        if (present(dealloc_vars)) dealloc_vars_loc = dealloc_vars
        
        ! create metric equilibrium variables
        call eq_2%init(grid_eq)
        
        ! do some preparations depending on equilibrium style used
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate  the  cylindrical  variables  R, Z  and  lambda  and
                ! derivatives if not yet done
                call writo('Calculate R,Z,L...')
                if (.not.allocated(grid_eq%trigon_factors)) then
                    ierr = calc_trigon_factors(grid_eq%theta_E,grid_eq%zeta_E,&
                        &grid_eq%trigon_factors)
                    CHCKERR('')
                end if
                do id = 0,max_deriv+1
                    ierr = calc_RZL(grid_eq,eq_2,derivs(id))
                    CHCKERR('')
                end do
            case (2)                                                            ! HELENA
                ! do nothing
        end select
        
        ! Calcalations depending on equilibrium style being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate the metrics in the cylindrical coordinate system
                call writo('Calculate g_C')                                     ! h_C is not necessary
                do id = 0,max_deriv
                    ierr = calc_g_C(eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the transformation matrix C(ylindrical) -> V(MEC)
                call writo('Calculate T_VC')
                do id = 0,max_deriv
                    ierr = calc_T_VC(eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the metric factors in the VMEC coordinate system
                call writo('Calculate g_V')
                do id = 0,max_deriv
                    ierr = calc_g_V(eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the jacobian in the VMEC coordinate system
                call writo('Calculate jac_V')
                do id = 0,max_deriv
                    ierr = calc_jac_V(grid_eq,eq_2,derivs(id))
                    CHCKERR('')
                end do
                
#if ldebug
                if (ltest) then
                    call writo('Test calculation of g_V?')
                    if(get_log(.false.)) then
                        ierr = test_g_V(grid_eq,eq_2)
                        CHCKERR('')
                        call pause_prog
                    end if
                    call writo('Test calculation of jac_V?')
                    if(get_log(.false.)) then
                        ierr = test_jac_V(grid_eq,eq_2)
                        CHCKERR('')
                        call pause_prog
                    end if
                end if
#endif
                
                ! calculate the transformation matrix V(MEC) -> F(lux)
                call writo('Calculate T_VF')
                do id = 0,max_deriv
                    ierr = calc_T_VF(grid_eq,eq_1,eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! set up plus minus one
                pmone = -1                                                      ! conversion VMEC LH -> RH coord. system
                
                ! possibly deallocate
                if (dealloc_vars_loc) then
                    deallocate(eq_2%L_E)
                    deallocate(eq_2%g_C)
                end if
            case (2)                                                            ! HELENA
#if ldebug
                ! Check whether poloidal grid is indeed chi_H
                ! If not,  the boundary  conditions used in  the splines  of the
                ! following routines are not correct.
                ! A  solution to  this problem  would be  to set  up the  HELENA
                ! variables  in  a  full  periodic  grid,  even  when  they  are
                ! top-bottom symmetric.
                do kd = 1,grid_eq%loc_n_r
                    do jd = 1,grid_eq%n(2)
                        do id = 1,grid_eq%n(1)
                            if (abs(grid_eq%theta_E(id,jd,kd)-chi_H(id)).gt.&
                                &tol_zero) then
                                ierr = 1
                                err_msg = 'theta_E is not identical to chi_H'
                                CHCKERR(err_msg)
                            end if
                            if (abs(grid_eq%zeta_E(id,jd,kd)-0._dp).gt.&
                                &tol_zero) then
                                ierr = 1
                                err_msg = 'zeta_E is not identical to 0'
                                CHCKERR(err_msg)
                            end if
                        end do
                    end do
                end do
                
                if (ltest) then
                    call writo('Test consistency of metric factors?')
                    if(get_log(.false.,ind=.true.)) then
                        ierr = test_metrics_H()
                        CHCKERR('')
                        call pause_prog(ind=.true.)
                    end if
                    
                    call writo('Test harmonic content of HELENA output?')
                    if(get_log(.false.,ind=.true.)) then
                        ierr = test_harm_cont_H()
                        CHCKERR('')
                        call pause_prog(ind=.true.)
                    end if
                end if
#endif
                
                ! calculate the jacobian in the HELENA coordinate system        
                call writo('Calculate jac_H')
                do id = 0,max_deriv
                    ierr = calc_jac_H(grid_eq,eq_1,eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the metric factors in the HELENA coordinate system
                call writo('Calculate g_H')
                do id = 0,max_deriv
                    ierr = calc_g_H(grid_eq,eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the inverse h_H of the metric factors g_H
                call writo('Calculate h_H')
                do id = 0,max_deriv
                    ierr = calc_inv_met(eq_2%h_E,eq_2%g_E,derivs(id))
                    CHCKERR('')
                end do
                
#if ldebug
                if (ltest) then
                    call writo('Test calculation of D1 D2 h_H?')
                    if(get_log(.false.)) then
                        ierr = test_D12h_H(grid_eq,eq_2)
                        CHCKERR('')
                        call pause_prog
                    end if
                end if
#endif
                
                ! export for VMEC port
                if (export_HEL) then
                    call writo('Exporting HELENA equilibrium for VMEC porting')
                    call lvl_ud(1)
                    ierr = create_VMEC_input(grid_eq,eq_1)
                    CHCKERR('')
                    call lvl_ud(-1)
                    call writo('Done exporting')
                end if
                
                ! calculate the transformation matrix H(ELENA) -> F(lux)        
                call writo('Calculate T_HF')
                do id = 0,max_deriv
                    ierr = calc_T_HF(grid_eq,eq_1,eq_2,derivs(id))
                    CHCKERR('')
                end do
                
                ! set up plus minus one
                pmone = 1
                
                ! possibly deallocate
                if (dealloc_vars_loc) then
                    deallocate(eq_2%h_E)
                end if
        end select
        
#if ldebug
        if (ltest) then
            call writo('Test calculation of T_EF?')
            if(get_log(.false.)) then
                ierr = test_T_EF(grid_eq,eq_1,eq_2)
                CHCKERR('')
                call pause_prog
            end if
        end if
#endif
        
        ! calculate  the  inverse of  the transformation  matrix T_EF
        call writo('Calculate T_FE')
        do id = 0,max_deriv
            ierr = calc_inv_met(eq_2%T_FE,eq_2%T_EF,derivs(id))
            CHCKERR('')
            ierr = calc_inv_met(eq_2%det_T_FE,eq_2%det_T_EF,derivs(id))
            CHCKERR('')
        end do
        
        ! calculate the metric factors in the Flux coordinate system
        call writo('Calculate g_F')
        do id = 0,max_deriv
            ierr = calc_g_F(eq_2,derivs(id))
            CHCKERR('')
        end do
        
        ! calculate the inverse h_F of the metric factors g_F
        call writo('Calculate h_F')
        do id = 0,max_deriv
            ierr = calc_inv_met(eq_2%h_F,eq_2%g_F,derivs(id))
            CHCKERR('')
        end do
        
        ! calculate the jacobian in the Flux coordinate system
        call writo('Calculate jac_F')
        do id = 0,max_deriv
            ierr = calc_jac_F(eq_2,derivs(id))
            CHCKERR('')
        end do
        
        ! limit Jacobians to small value to avoid infinities
        eq_2%jac_F(:,:,:,0,0,0) = sign(&
            &max(tol_zero,abs(eq_2%jac_F(:,:,:,0,0,0))),eq_2%jac_F(:,:,:,0,0,0))
        eq_2%jac_E(:,:,:,0,0,0) = sign(&
            &max(tol_zero,abs(eq_2%jac_E(:,:,:,0,0,0))),eq_2%jac_E(:,:,:,0,0,0))
        
        ! possibly deallocate
        if (dealloc_vars_loc) then
            deallocate(eq_2%g_E)
        end if
        
        ! Transform metric equilibrium E into F derivatives
        ierr = calc_F_derivs(eq_2)
        CHCKERR('')
        
        ! possibly deallocate
        if (dealloc_vars_loc) then
            deallocate(eq_2%g_F,eq_2%h_F,eq_2%jac_F)
            deallocate(eq_2%T_FE,eq_2%det_T_EF,eq_2%det_T_FE)
        end if
        
        ! Calculate derived metric quantities
        ierr = calc_derived_q(grid_eq,eq_1,eq_2)
        CHCKERR('')
        
#if ldebug
        if (ltest) then
            call writo('Test Jacobian in Flux coordinates?')
            if(get_log(.false.)) then
                ierr = test_jac_F(grid_eq,eq_1,eq_2)
                CHCKERR('')
                call pause_prog
            end if
            call writo('Test calculation of B_F?')
            if(get_log(.false.)) then
                ierr = test_B_F(grid_eq,eq_1,eq_2)
                CHCKERR('')
                call pause_prog
            end if
            call writo('Test consistency with given pressure?')
            if(get_log(.false.)) then
                ierr = test_p(grid_eq,eq_1,eq_2)
                CHCKERR('')
                call pause_prog
            end if
        end if
#endif
        
        call lvl_ud(-1)
        
        call writo('Done setting up metric equilibrium quantities')
    end function calc_eq_2
    
    !> Creates a VMEC input file.
    !!
    !! Optionally, a perturbation  can be added: Either the  displacement of the
    !! plasma position  can be  described (\c  pert_style 1),  or ripple  in the
    !! toroidal magnetic  field (\c  pert_style 2), with  a fixed  toroidal mode
    !! number.
    !!
    !! Both perturbation styles can have various prescription types:
    !!  -# file with Fourier modes in the geometrical angular coordinate
    !!  -# same but manually
    !!  -#  file with  perturbation  from a  2-D map  in  the geometric  angular
    !!  coordinate.
    !!
    !! For  \c pert_style  2,  a file  has  to be  provided  that describes  the
    !! translation between  position perturbation and magnetic  perturbation for
    !! curves of constant  geometrical angle. This file can be  generated for an
    !! already existing  ripple case using POST  with <tt>--compare_tor_pos</tt>
    !! with <tt>n_zeta_plot = 3</tt> and \c min_theta_plot and \c max_theta_plot
    !! indicating half a ripple period.
    !!
    !! The output from  this VMEC run can  then be used to  iteratively create a
    !! new  file  to  translate  toroidal  magnetic  field  ripple  to  position
    !! perturbation.
    !!
    !! \note Meaning of the indices of \c B_F, \c B_F_loc:
    !!  - <tt>(pol modes, cos/sin)</tt> for \c B_F_loc
    !!  - <tt>(tor_modes, pol modes, cos/sin (m theta), R/Z)</tt> for \c B_F
    integer function create_VMEC_input(grid_eq,eq_1) result(ierr)
        use eq_vars, only: pres_0, R_0, psi_0, B_0, vac_perm
        !use eq_vars, only: max_flux_E
        use grid_vars, only: n_r_eq
        use grid_utilities, only: nufft, calc_XYZ_grid
        use HELENA_vars, only: nchi, R_H, Z_H, ias, RBphi_H, flux_p_H, &
            &flux_t_H, chi_H, RMtoG_H, BMtoG_H
        use X_vars, only: min_r_sol, max_r_sol
        use input_utilities, only: pause_prog, get_log, get_int, get_real
        use num_utilities, only: GCD, bubble_sort, order_per_fun, &
            &shift_F, calc_int, spline
        use num_vars, only: eq_name, HEL_pert_i, HEL_export_i, &
            &norm_disc_prec_eq, prop_B_tor_i, tol_zero, mu_0_original, &
            &use_normalization
        use files_utilities, only: skip_comment, count_lines
        use EZspline_obj
        use EZspline
        
        character(*), parameter :: rout_name = 'create_VMEC_input'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid varibles
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium quantities
        
        ! local variables
        type(EZspline2_r8) :: f_spl                                             ! spline object for interpolation
        integer :: id, jd, kd                                                   ! counters
        integer :: n_B                                                          ! nr. of points in Fourier series
        integer :: nfp                                                          ! scale factor for toroidal mode numbers
        integer :: tot_nr_pert                                                  ! total number of perturbations combinations (N,M)
        integer :: nr_n                                                         ! number of different N in perturbation
        integer :: id_n_0                                                       ! index where zero N is situated in B_F
        integer :: n_loc, m_loc                                                 ! local n and m of perturbation
        integer :: pert_map_n_loc                                               ! n_loc for perturbation map
        integer :: n_id                                                         ! index in bundled n_pert
        integer :: n_pert_map(2)                                                ! number of points in perturbation map
        integer :: plot_dim(2)                                                  ! plot dimensions
        integer :: rec_min_m                                                    ! recommended minimum of poloidal modes
        integer :: pert_style                                                   ! style of perturbation (1: r, 2: B_tor)
        integer :: pert_type                                                    ! type of perturbation prescription
        integer :: max_n_B_output                                               ! max. nr. of modes written in output file (constant in VMEC)
        integer :: n_prop_B_tor                                                 ! n of angular poins in prop_B_tor
        integer :: nr_m_max                                                     ! maximum nr. of poloidal mode numbers
        integer :: ncurr                                                        ! ncurr VMEC flag (0: prescribe iota, 1: prescribe tor. current)
        integer :: bcs(2,2)                                                     ! boundary conditions
        integer, allocatable :: n_pert(:)                                       ! tor. mode numbers and pol. mode numbers for each of them
        integer, allocatable :: n_pert_copy(:)                                  ! copy of n_pert, for sorting
        integer, allocatable :: piv(:)                                          ! pivots for sorting
        character(len=1) :: pm(3)                                               ! "+ " or "-"
        character(len=8) :: flux_name(2)                                        ! "poloidal" or "toroidal"
        character(len=6) :: path_prefix = '../../'                              ! prefix of path
        character(len=max_str_ln) :: HEL_pert_file_name                         ! name of perturbation file
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: file_name                                  ! name of file
        character(len=max_str_ln) :: plot_name(2)                               ! name of plot file
        character(len=max_str_ln) :: plot_title(2)                              ! name of plot
        character(len=max_str_ln) :: prop_B_tor_file_name                       ! name of B_tor proportionality file
        real(dp) :: eq_vert_shift                                               ! vertical shift in equilibrium
        real(dp) :: pert_map_vert_shift                                         ! vertical shift in perturbation map
        real(dp) :: max_pert_on_axis                                            ! maximum perturbation on axis (theta = 0)
        real(dp) :: m_tol = 1.E-7_dp                                            ! tolerance for Fourier mode strength
        real(dp) :: delta_loc(2)                                                ! local delta
        real(dp) :: plot_lims(2,2)                                              ! limits of plot dims [pi]
        real(dp) :: norm_B_H                                                    ! normalization for R and Z Fourier modes
        real(dp) :: mult_fac                                                    ! global multiplication factor
        real(dp) :: prop_B_tor_smooth                                           ! smoothing of prop_B_tor_interp via Fourier transform
        real(dp) :: RZ_B_0(2)                                                   ! origin of R and Z of boundary
        real(dp) :: loc_F(0:1,2)                                                ! local trigonometric variable for holding cos and sin
        real(dp) :: s_VJ(99)                                                    ! normal coordinate s for writing
        real(dp) :: pres_VJ(99,0:1)                                             ! pressure for writing
        real(dp) :: rot_T_V(99)                                                 ! rotational transform for writing
        real(dp) :: F_VJ(99)                                                    ! F for writing
        real(dp) :: FFp_VJ(99)                                                  ! FF' for writing
        real(dp) :: flux_J(99)                                                  ! normalized flux for writing
        real(dp) :: q_saf_VJ(99)                                                ! q_saf for writing
        real(dp) :: I_tor_V(99)                                                 ! enclosed toroidal current for VMEC
        real(dp) :: I_tor_J(99)                                                 ! enclosed toroidal current for JOREK
        real(dp) :: norm_J(3)                                                   ! normalization factors for Jorek (R_geo, Z_geo, F0)
        real(dp), allocatable :: psi_t(:)                                       ! normalized toroidal flux
        real(dp), allocatable :: FFp(:)                                         ! FF' in total grid
        real(dp), allocatable :: I_tor(:)                                       ! toroidal current within flux surfaces
        real(dp), allocatable :: norm_transf(:,:)                               ! transformation of normal coordinates between MISHKA and VMEC or JOREK
        real(dp), allocatable :: I_tor_int(:)                                   ! integrated I_tor
        real(dp), allocatable :: RRint(:)                                       ! R^2 integrated poloidally
        real(dp), allocatable :: RRint_loc(:)                                   ! local RRint
        real(dp), allocatable :: R_plot(:,:,:)                                  ! R for plotting of ripple map
        real(dp), allocatable :: Z_plot(:,:,:)                                  ! Z for plotting of ripple map
        real(dp), allocatable :: pert_map_R(:)                                  ! R of perturbation map
        real(dp), allocatable :: pert_map_Z(:)                                  ! Z of perturbation map
        real(dp), allocatable :: pert_map(:,:)                                  ! perturbation map
        real(dp), allocatable :: pert_map_interp(:)                             ! interpolated perturbation from map
        real(dp), allocatable :: pert_map_interp_F(:,:)                         ! Fourier coefficients of pert_map_interp
        real(dp), allocatable :: R_H_loc(:,:)                                   ! local R_H
        real(dp), allocatable :: Z_H_loc(:,:)                                   ! local Z_H
        real(dp), allocatable :: delta(:,:,:)                                   ! amplitudes of perturbations (N,M,c/s)
        real(dp), allocatable :: delta_copy(:,:,:)                              ! copy of delta, for sorting
        real(dp), allocatable :: BH_0(:,:)                                      ! R and Z at unperturbed bounday
        real(dp), allocatable :: B_F(:,:,:,:)                                   ! cos and sin Fourier components
        real(dp), allocatable :: B_F_loc(:,:)                                   ! Local B_F
        real(dp), allocatable :: theta_B(:)                                     ! geometric pol. angle
        real(dp), allocatable :: theta_geo(:,:,:)                               ! geometric pol. angle, for plotting
        real(dp), allocatable :: prop_B_tor(:,:)                                ! proportionality between delta B_tor and delta_norm
        real(dp), allocatable :: prop_B_tor_ord(:,:)                            ! ordered prop_B_tor
        real(dp), allocatable :: prop_B_tor_interp(:)                           ! proportionality between delta B_tor and delta_norm
        real(dp), allocatable :: prop_B_tor_interp_F(:,:)                       ! Fourier modes of prop_B_tor_interp
        real(dp), allocatable :: XYZ_plot(:,:,:,:)                              ! plotting X, Y and Z
        logical :: zero_N_pert                                                  ! there is a perturbation with N = 0
        logical :: pert_eq                                                      ! whether equilibrium is perturbed
        logical :: stel_sym                                                     ! whether there is stellarator symmetry
        logical :: change_max_n_B_output                                        ! whether to change max_n_B_output
        logical :: found                                                        ! whether something was found
#if ldebug
        real(dp), allocatable :: BH_0_ALT(:,:)                                  ! reconstructed R and Z
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! test if full range
        if (min_r_sol.gt.0 .or. max_r_sol.lt.1) then
            ierr = 1
            err_msg = 'This routine should be run with the full normal &
                &range!'
            CHCKERR(err_msg)
        end if
        
        ! test if normalization constants chosen
        if (.not.all(BR_normalization_provided,1)) &
            &call writo('No normalization factors were provided. &
            &Are you sure this is correct?',warning=.true.)
        
        ! set up local R_H and Z_H
        allocate(R_H_loc(nchi,n_r_eq))
        allocate(Z_H_loc(nchi,n_r_eq))
        R_H_loc = R_H
        Z_H_loc = Z_H
        if (use_normalization) then
            R_H_loc = R_H*R_0
            Z_H_loc = Z_H*R_0
        end if
        
        ! ask for vertical shift
        ! Note: HELENA automatically shifts it zo Zvac = 0
        call writo('Was the equilibrium shifted vertically?')
        eq_vert_shift = 0._dp
        if (get_log(.false.)) then
            call writo('How much higher [m] should the &
                &equilibrium be in the real world?')
            eq_vert_shift = get_real()
            call writo('The equilibrium will be shifted by '//&
                &trim(r2strt(eq_vert_shift))//'m to compensate for this')
            Z_H_loc = Z_H_loc + eq_vert_shift
        end if
        
        ! set up F, FF' and <R^2>
        allocate(FFp(n_r_eq))
        allocate(RRint(n_r_eq))
        allocate(RRint_loc(nchi))                                               ! includes overlap, as necessary for integration
        ierr = spline(flux_p_H(:,0)/(2*pi),RBphi_H(:,0)**2*0.5_dp,&
            &flux_p_H(:,0)/(2*pi),FFp,ord=norm_disc_prec_eq,deriv=1)
        CHCKERR('')
        do kd = 1,n_r_eq
            ierr = calc_int(R_H(:,kd)**2,chi_H,RRint_loc)
            CHCKERR('')
            RRint(kd) = RRint_loc(nchi)
        end do
        if (ias.eq.0) RRint = 2*RRint
        
        ! set up toroidal current:
        ! This is the  toroidal current between neighbouring  flux surfaces (see
        ! first commentary in VMEC2013/Sources/General/add_fluxes.f90:
        !   <J^zeta J>  = - [2pi FF'/mu_0 <J/R^2> + p'<J>]
        !               = - [2pi F'q/mu_0 + p' <R^2>q/F] ,
        ! which all  have units  1/m. Here, <.>  means integration  in a
        ! full poloidal period.
        ! As VMEC  wants the toroidal current  within two infintesimally
        ! separated  flux  surface as  a  function  of the  VMEC  normal
        ! coordinate, above has to be multiplied by
        !   dpsi_PB3D/dpsi_V = psi_tor(edge)/(2pi q)
        ! Finally, there is  an extra factor -1 because  the toroidal coordinate
        ! in VMEC is oposite.
        ! For JOREK,  the same  is done  with the poloidal  flux instead  of the
        ! toroidal and without the change of sign.
        allocate(I_tor(n_r_eq))
        allocate(norm_transf(n_r_eq,2))
        norm_transf(:,1) = -flux_t_H(n_r_eq,0)/(2*pi*eq_1%q_saf_E(:,0))
        norm_transf(:,2) = flux_p_H(n_r_eq,0)/(2*pi)
        I_tor = - (eq_1%pres_E(:,1)*RRint*eq_1%q_saf_E(:,0)/RBphi_H(:,0) + &
            &2*pi*FFp*eq_1%q_saf_E(:,0)/(RBphi_H(:,0)*vac_perm))
        
        ! calculate total toroidal current
        allocate(I_tor_int(size(I_tor)))
        ierr = calc_int(I_tor*norm_transf(:,1),&
            &flux_t_H(:,0)/flux_t_H(n_r_eq,0),I_tor_int)
        if (use_normalization) I_tor_int = I_tor_int * R_0*B_0/mu_0_original
        CHCKERR('')
        
        ! user output
        call writo('Initialize boundary')
        call lvl_ud(1)
        
        ! full 0..2pi boundary
        if (ias.eq.0) then                                                      ! symmetric (so nchi is odd)
            n_B = nchi*2-2                                                      ! -2 because of overlapping points at 0, pi
        else                                                                    ! asymmetric (so nchi is aumented by one)
            n_B  = nchi-1                                                       ! -1 because of overlapping points at 0
        end if
        nr_m_max = (n_B-1)/2
        allocate(BH_0(n_B,2))                                                   ! HELENA coords (no overlap at 0)
        allocate(theta_B(n_B))                                                  ! geometrical poloidal angle at boundary (overlap at 0)
        if (ias.eq.0) then                                                      ! symmetric
            BH_0(1:nchi,1) = R_H_loc(:,n_r_eq)
            BH_0(1:nchi,2) = Z_H_loc(:,n_r_eq)
            do id = 1,nchi-2
                BH_0(nchi+id,1) = R_H_loc(nchi-id,n_r_eq)
                BH_0(nchi+id,2) = -Z_H_loc(nchi-id,n_r_eq)
            end do
        else
            BH_0(:,1) = R_H_loc(1:n_B,n_r_eq)
            BH_0(:,2) = Z_H_loc(1:n_B,n_r_eq)
        end if
        RZ_B_0(1) = sum(R_H_loc(:,1))/size(R_H_loc,1)
        RZ_B_0(2) = sum(Z_H_loc(:,1))/size(Z_H_loc,1)
        theta_B = atan2(BH_0(:,2)-RZ_B_0(2),BH_0(:,1)-RZ_B_0(1))
        where (theta_B.lt.0) theta_B = theta_B + 2*pi
        call writo('Magnetic axis used for geometrical coordinates:')
        call lvl_ud(1)
        call writo('('//trim(r2str(RZ_B_0(1)))//'m, '//&
            &trim(r2str(RZ_B_0(2)))//'m)')
        call lvl_ud(-1)
        
        ! plot with HDF5
        allocate(theta_geo(grid_eq%n(1),1,grid_eq%loc_n_r))
        do kd = 1,grid_eq%loc_n_r
            theta_geo(:,1,kd) = atan2(Z_H_loc(:,kd)-RZ_B_0(2),&
                &R_H_loc(:,kd)-RZ_B_0(1))
        end do
        where (theta_geo.lt.0._dp) theta_geo = theta_geo + 2*pi
        allocate(XYZ_plot(grid_eq%n(1),1,grid_eq%loc_n_r,3))
        ierr = calc_XYZ_grid(grid_eq,grid_eq,XYZ_plot(:,:,:,1),&
            &XYZ_plot(:,:,:,2),XYZ_plot(:,:,:,3))
        CHCKERR('')
        call plot_HDF5('theta_geo','theta_geo',theta_geo,&
            &tot_dim=[grid_eq%n(1),1,grid_eq%n(3),3],&
            &loc_offset=[0,0,grid_eq%i_min-1,0],&
            &x=XYZ_plot(:,:,:,1),y=XYZ_plot(:,:,:,2),z=XYZ_plot(:,:,:,3),&
            &descr='geometric poloidal angle used to create the VMEC &
            &input file')
        deallocate(XYZ_plot)
        
        !!! TEMPORARILY
        !!write(*,*) 'TEMPORARILY CIRCULAR TOKAMAK !!!!!!!!!!!'
        !!BH_0(:,1) = 1.5_dp + 0.5*cos(theta_B)
        !!BH_0(:,2) = 0.5*sin(theta_B)
        
        ! plot R and Z
        plot_name(1) = 'RZ'
        plot_title = ['R_H','Z_H']
        call print_ex_2D(plot_title(1:2),plot_name(1),BH_0,&
            &x=reshape(theta_B,[n_B,1])/pi,&
            &draw=.false.)
        call draw_ex(plot_title(1:2),plot_name(1),2,1,.false.)
        
        ! set up auxiliary variable
        flux_name = ['poloidal','toroidal']
        
        call lvl_ud(-1)
        
        ! plot properties
        plot_dim = [100,20]
        plot_lims(:,1) = [0.0_dp,2.0_dp]
        plot_lims(:,2) = [0.5_dp,2.0_dp]
        call writo('Change plot properties from defaults?')
        call lvl_ud(1)
        do id = 1,2
            call writo(trim(i2str(plot_dim(id)))//' geometrical '//&
                &flux_name(id)//' points on range '//&
                &trim(r2strt(plot_lims(1,id)))//' pi .. '//&
                &trim(r2strt(plot_lims(2,id)))//' pi')
        end do
        call lvl_ud(-1)
        if (get_log(.false.)) then
            call lvl_ud(1)
            do id = 1,2
                call writo('number of '//flux_name(id)//' points?')
                plot_dim(id) = get_int(lim_lo=2)
                call writo('lower limit of '//flux_name(id)//&
                    &' range [pi]?')
                plot_lims(1,id) = get_real()
                call writo('upper limit of '//flux_name(id)//&
                    &' range [pi]?')
                plot_lims(2,id) = get_real()
            end do
            call lvl_ud(-1)
        end if
        
        ! get ncurr
        call writo('VMEC current style?')
        call lvl_ud(1)
        call writo('0: prescribe iota')
        call writo('1: prescribe toroidal current')
        ncurr = get_int(lim_lo=0,lim_hi=1)
        call lvl_ud(-1)
        
        ! get input about equilibrium perturbation
        ! Note: negative  M values will  be converted  to negative N  values and
        ! Possibly negative delta's.
        call writo('Do you want to perturb the equilibrium?')
        pert_eq = get_log(.false.)
        
        zero_N_pert = .false.
        if (pert_eq) then
            ! get perturbation style
            call writo('Which perturbation do you want to describe?')
            call lvl_ud(1)
            call writo('1: plasma boundary position')
            call writo('2: B_tor magnetic ripple with fixed N')
            call lvl_ud(-1)
            pert_style = get_int(lim_lo=1,lim_hi=2)
            
            ! for style 2, get proportionality file
            if (pert_style.eq.2) then
                ! find file
                found = .false.
                do while (.not.found)
                    call writo('Proportionality file name?')
                    call lvl_ud(1)
                    read(*,*,IOSTAT=ierr) prop_B_tor_file_name
                    err_msg = 'failed to read prop_B_tor_file_name'
                    CHCKERR(err_msg)
                    
                    ! open
                    do kd = 0,2
                        call writo('Trying "'//path_prefix(1:kd*3)//&
                            &trim(prop_B_tor_file_name)//'"')
                        open(prop_B_tor_i,FILE=path_prefix(1:kd*3)//&
                            &trim(prop_B_tor_file_name),IOSTAT=ierr,&
                            &STATUS='old')
                        if (ierr.eq.0) then
                            call lvl_ud(1)
                            call writo("Success")
                            call lvl_ud(-1)
                            exit
                        end if
                    end do
                    if (ierr.eq.0) then
                        found = .true.
                    else
                        call writo('Could not open file "'//&
                            &trim(prop_B_tor_file_name)//'"')
                    end if
                    
                    call lvl_ud(-1)
                end do
                
                ! read file
                call writo('Analyzing perturbation proportionality file "'&
                    &//trim(prop_B_tor_file_name)//'"')
                call lvl_ud(1)
                n_prop_B_tor = count_lines(prop_B_tor_i)
                allocate(prop_B_tor(n_prop_B_tor,2))
                ierr = skip_comment(prop_B_tor_i,&
                    &file_name=prop_B_tor_file_name)
                CHCKERR('')
                do id = 1,n_prop_B_tor
                    read(prop_B_tor_i,*,IOSTAT=ierr) prop_B_tor(id,:)
                end do
                
                ! user info
                call writo(trim(i2str(n_prop_B_tor))//&
                    &' poloidal points '//&
                    &trim(r2strt(minval(prop_B_tor(:,1))))//'..'//&
                    &trim(r2strt(maxval(prop_B_tor(:,1)))))
                call writo('The proportionality factor should be tabulated &
                    &in a geometrical angle that has the SAME origin as &
                    &the one used here',alert=.true.)
                
                ! interpolate the proportionality  factor on the geometrical
                ! poloidal angle used  here (as opposed to the  one in which
                ! it was tabulated)
                ierr = order_per_fun(prop_B_tor,prop_B_tor_ord,1,&              ! overlap of 1 is enough for splines
                    &tol=0.5_dp*pi/n_prop_B_tor)                                ! set appropriate tolerance: a quarter of a equidistant grid
                CHCKERR('')
                deallocate(prop_B_tor)
                allocate(prop_B_tor_interp(n_B))
                ierr = spline(prop_B_tor_ord(:,1),prop_B_tor_ord(:,2),&
                    &theta_B,prop_B_tor_interp,ord=norm_disc_prec_eq,&
                    &bcs=[-1,-1])
                CHCKERR('')
                deallocate(prop_B_tor_ord)
                if (grid_eq%n(2).ne.1) call writo('There should be only 1 &
                    &geodesical position, but there are '//&
                    &trim(i2str(grid_eq%n(2))),warning=.true.)
                
                ! smooth prop_B_tor_interp
                prop_B_tor_smooth = 1._dp
                call writo('Do you want to smooth the perturbation?')
                if (get_log(.false.)) prop_B_tor_smooth = &
                    &get_real(lim_lo=0._dp,lim_hi=1._dp)
                
                ! calculate NUFFT
                ierr = nufft(theta_B,prop_B_tor_interp,&
                    &prop_B_tor_interp_F)
                CHCKERR('')
                deallocate(prop_B_tor_interp)
                
                call lvl_ud(-1)
            end if
            
            ! get perturbation type
            call writo('How do you want to prescribe the perturbation?')
            call lvl_ud(1)
            call writo('1: through a file with Fourier modes in &
                &geometrical coordinates')
            call lvl_ud(1)
            call writo('It should provide N M delta_cos delta_sin in four &
                &columns')
            call writo('Rows starting with "#" are ignored')
            call lvl_ud(-1)
            call writo('2: same as 1 but interactively')
            call writo('3: through a 2-D map in R, Z for constant &
                &geometrical coordinates')
            call lvl_ud(1)
            call writo('It should provide R, Z and delta')
            call writo('Rows starting with "#" are ignored')
            call lvl_ud(-1)
            pert_type = get_int(lim_lo=1,lim_hi=3)
            call lvl_ud(-1)
            
            select case (pert_type)
                case (1,3)                                                      ! files
                    found = .false.
                    do while (.not.found)
                        ! user input
                        call writo('File name?')
                        call lvl_ud(1)
                        read(*,*,IOSTAT=ierr) HEL_pert_file_name
                        err_msg = 'failed to read HEL_pert_file_name'
                        CHCKERR(err_msg)
                        
                        ! open
                        do kd = 0,2
                            call writo('Trying "'//path_prefix(1:kd*3)//&
                                &trim(HEL_pert_file_name)//'"')
                            open(HEL_pert_i,FILE=path_prefix(1:kd*3)//&
                                &trim(HEL_pert_file_name),IOSTAT=ierr,&
                                &STATUS='old')
                            if (ierr.eq.0) then
                                call lvl_ud(1)
                                call writo("Success")
                                call lvl_ud(-1)
                                exit
                            end if
                        end do
                        if (ierr.eq.0) then
                            found = .true.
                        else
                            call writo('Could not open file "'//&
                                &trim(HEL_pert_file_name)//'"')
                        end if
                        
                        call lvl_ud(-1)
                    end do
                    
                    call writo('Parsing "'//trim(HEL_pert_file_name)//'"')
                    call lvl_ud(1)
                    
                    if (pert_type.eq.1) then                                    ! Fourier modes in geometrical coordinates
                        ! get number of lines
                        tot_nr_pert = count_lines(HEL_pert_i)
                    else if (pert_type.eq.3) then                               ! 2-D map of perturbation
                        ! ask for vertical shift
                        call writo('Was the equilibrium shifted &
                            &vertically?')
                        if (abs(eq_vert_shift).gt.0._dp) then
                            call lvl_ud(1)
                            call writo('Note: If you already shifted just &
                                &now, don''t do it again!',alert=.true.)
                            call lvl_ud(-1)
                        end if
                        pert_map_vert_shift = 0._dp
                        if (get_log(.false.)) then
                            call writo('How much higher [m] should the &
                                &equilibrium be in the real world?')
                            pert_map_vert_shift = get_real()
                            call writo('The perturbation map will be &
                                &shifted by '//&
                                &trim(r2strt(-pert_map_vert_shift))//&
                                &'m to compensate for this')
                        end if
                        
                        ! get map
                        ierr = skip_comment(HEL_pert_i,&
                            &file_name=HEL_pert_file_name)
                        CHCKERR('')
                        read(HEL_pert_i,*) n_pert_map(1)                        ! nr. of points in R
                        allocate(pert_map_R(n_pert_map(1)))
                        do kd = 1,n_pert_map(1)
                            read(HEL_pert_i,*) pert_map_R(kd)
                        end do
                        ierr = skip_comment(HEL_pert_i,&
                            &file_name=HEL_pert_file_name)
                        CHCKERR('')
                        read(HEL_pert_i,*) n_pert_map(2)                        ! nr. of points in R
                        allocate(pert_map_Z(n_pert_map(2)))
                        do kd = 1,n_pert_map(2)
                            read(HEL_pert_i,*) pert_map_Z(kd)
                        end do
                        pert_map_Z = pert_map_Z - pert_map_vert_shift
                        ierr = skip_comment(HEL_pert_i,&
                            &file_name=HEL_pert_file_name)
                        CHCKERR('')
                        allocate(pert_map(n_pert_map(1),n_pert_map(2)))
                        do kd = 1,n_pert_map(1)
                            read(HEL_pert_i,*) pert_map(kd,:)
                        end do
                        
                        ! plot map
                        allocate(R_plot(n_pert_map(1),n_pert_map(2),1))
                        allocate(Z_plot(n_pert_map(1),n_pert_map(2),1))
                        do kd = 1,n_pert_map(2)
                            R_plot(:,kd,1) = pert_map_R
                        end do
                        do kd = 1,n_pert_map(1)
                            Z_plot(kd,:,1) = pert_map_Z
                        end do
                        call plot_HDF5('pert_map','pert_map',&
                            &reshape(pert_map,[n_pert_map,1]),&
                            &x=R_plot,y=Z_plot,descr=&
                            &'perturbation map for '//&
                            &trim(HEL_pert_file_name))
                        deallocate(R_plot,Z_plot)
                        
                        ! test
                        if (minval(R_H_loc).lt.minval(pert_map_R) .or. &
                            &maxval(R_H_loc).gt.maxval(pert_map_R) .or. &
                            &minval(Z_H_loc).lt.minval(pert_map_Z) .or. &
                            &maxval(Z_H_loc).gt.maxval(pert_map_Z)) then
                            ierr = 1
                            call writo('Are you sure you have specified &
                                &a normalization factor R_0?',alert=.true.)
                            err_msg = 'R and Z are not contained in &
                                &perturbation map'
                            CHCKERR(err_msg)
                        end if
                        
                        ! set up 2-D spline
                        allocate(pert_map_interp(n_b))
                        bcs(:,1) = [0,0]                                        ! not a knot
                        bcs(:,2) = [0,0]                                        ! not a knot
                        call EZspline_init(f_spl,n_pert_map(1),n_pert_map(2),&
                            &bcs(:,1),bcs(:,2),ierr)
                        call EZspline_error(ierr)
                        CHCKERR('')
                        
                        ! set grid
                        f_spl%x1 = pert_map_R
                        f_spl%x2 = pert_map_Z
                        
                        ! set up 
                        call EZspline_setup(f_spl,pert_map,ierr,&
                            &exact_dim=.true.)                                  ! match exact dimensions, none of them old Fortran bullsh*t!
                        call EZspline_error(ierr)
                        CHCKERR('')
                        
                        ! interpolate
                        do id = 1,n_b
                            call EZspline_interp(f_spl,BH_0(id,1),BH_0(id,2),&
                                &pert_map_interp(id),ierr)
                            call EZspline_error(ierr)
                            CHCKERR('')
                        end do
                        call print_ex_2D('ripple','pert_map_interp',&
                            &pert_map_interp,x=theta_B,draw=.false.)
                        call draw_ex(['ripple'],'pert_map_interp',&
                            &1,1,.false.)
                        
                        ! calculate NUFFT
                        ierr = nufft(theta_B,pert_map_interp,&
                            &pert_map_interp_F)
                        CHCKERR('')
                        tot_nr_pert = size(pert_map_interp_F,1)*2               ! need both positive and negative n
                        
                        ! get toroidal mode number
                        call writo('toiroidal period for 2-D map?')
                        pert_map_n_loc = get_int(lim_lo=0)
                    end if
                    
                    call lvl_ud(-1)
                case (2)                                                        ! interactively
                    ! user input
                    call writo('Prescribe interactively')
                    call writo('How many different combinations of &
                        &toroidal and poloidal mode numbers (N,M)?')
                    tot_nr_pert = get_int(lim_lo=1)
            end select
            
            ! set up n, m and delta
            allocate(n_pert(tot_nr_pert+1))
            allocate(delta(tot_nr_pert+1,0:nr_m_max,2))
            delta = 0._dp
            nr_n = 0
            do jd = 1,tot_nr_pert
                select case (pert_type)
                    case (1)                                                    ! file with Fourier modes in geometrical coordinates
                        ierr = skip_comment(HEL_pert_i,&
                            &file_name=HEL_pert_file_name)
                        CHCKERR('')
                        read(HEL_pert_i,*,IOSTAT=ierr) n_loc, m_loc, &
                            &delta_loc
                        err_msg = 'Could not read file '//&
                            &trim(HEL_pert_file_name)
                        CHCKERR(err_msg)
                    case (2)                                                    ! same as 1 but interactively
                        call writo('For mode '//trim(i2str(jd))//'/'//&
                            &trim(i2str(tot_nr_pert))//':')
                        call lvl_ud(1)
                        
                        ! input
                        call writo('Toroidal mode number N?')
                        n_loc = get_int()
                        call writo('Poloidal mode number M?')
                        m_loc = get_int()
                        call writo('Perturbation strength ~ cos('//&
                            &trim(i2str(m_loc))//' θ - '//&
                            &trim(i2str(n_loc))//' ζ)?')
                        delta_loc(1) = get_real()
                        call writo('Perturbation strength ~ sin('//&
                            &trim(i2str(m_loc))//' θ - '//&
                            &trim(i2str(n_loc))//' ζ)?')
                        delta_loc(2) = get_real()
                        
                        call lvl_ud(-1)
                    case (3)                                                    ! file with Fourier modes in general coordinates
                        m_loc = (jd-1)/2                                        ! so the output is 0,0,1,1,2,2,3,3,4,4,...
                        delta_loc = pert_map_interp_F(m_loc+1,:)*0.5_dp         ! both + and - helicity, so need to divide amplitude by 2
                        n_loc = pert_map_n_loc
                        if (mod(jd,2).eq.0) n_loc = -n_loc
                end select
                
                ! has N = 0 been prescribed?
                if (n_loc.eq.0) zero_N_pert = .true.
                
                ! possibly convert to positive M
                if (m_loc.lt.0) then
                    n_loc = -n_loc
                    m_loc = -m_loc
                    delta_loc(2) = -delta_loc(2)
                end if
                
                ! bundle into n_pert
                n_id = nr_n+1                                                   ! default new N in next index
                do id = 1,nr_n
                    if (n_loc.eq.n_pert(id)) n_id = id
                end do
                if (n_id.gt.nr_n) then                                          ! new N
                    nr_n = n_id                                                 ! increment number of N
                    n_pert(n_id) = n_loc                                        ! with value n_loc
                end if
                if (m_loc.gt.nr_m_max) then                                     ! enlarge delta
                    allocate(delta_copy(tot_nr_pert+1,0:nr_m_max,2))
                    delta_copy = delta
                    deallocate(delta)
                    allocate(delta(tot_nr_pert+1,0:m_loc,2))
                    delta(:,0:nr_m_max,:) = delta_copy
                    deallocate(delta_copy)
                    nr_m_max = m_loc
                end if
                delta(n_id,m_loc,:) = delta(n_id,m_loc,:) + delta_loc           ! and perturbation amplitude for cos and sin with value delta_loc
            end do
            
            if (pert_type.eq.1 .or. pert_type.eq.3) close(HEL_pert_i)
            
            ! ask for global rescaling
            ! Note: This  is only valid if R(theta=0) is  the place with the    
            ! maximum perturbation
            max_pert_on_axis = sum(delta(:,:,1))
            select case (pert_style)
                case (1)                                                        ! plasma boundary position
                    call writo('Maximum absolute perturbation on axis: '//&
                        &trim(r2strt(100*max_pert_on_axis))//'cm')
                    call writo('Rescale to some value [cm]?')
                case (2)                                                        ! B_tor magnetic ripple with fixed N
                    call writo('Maximum relative perturbation on axis: '//&
                        &trim(r2strt(100*max_pert_on_axis))//'%')
                    call writo('Rescale to some value [%]?')
            end select
            mult_fac = max_pert_on_axis
            if (get_log(.false.)) then
                mult_fac = get_real()
                mult_fac = mult_fac * 0.01_dp
            end if
            delta = delta * mult_fac / max_pert_on_axis
            
            ! add zero to the peturbations if it's not yet there
            ! (for the bubble sort)
            if (.not.zero_N_pert) then
                nr_n = nr_n+1
                n_pert(nr_n) = 0
            end if
        else
            nr_n = 1
        end if
        
        ! user output
        call writo('Set up axisymmetric data')
        call lvl_ud(1)
        
        ! calculate output for unperturbed R and Z
        plot_name(1) = 'R_F'
        plot_name(2) = 'Z_F'
        allocate(B_F(nr_n,0:nr_m_max,2,2))                                      ! (tor modes, pol modes, cos/sin (m theta), R/Z)
        B_F = 0._dp
        do kd = 1,2
            ! user output
            call writo('analyzing '//trim(plot_name(kd)))
            call lvl_ud(1)
            
            ! NUFFT
            ierr = nufft(theta_B,BH_0(:,kd),B_F_loc,plot_name(kd))
            CHCKERR('')
            B_F(1,:,:,kd) = B_F_loc
        
#if ldebug
            if (debug_create_VMEC_input) then
                allocate(BH_0_ALT(size(BH_0,1),size(B_F_loc,1)))
                
                call writo('Comparing R or Z with reconstruction &
                    &through Fourier coefficients')
                BH_0_ALT(:,1) = B_F_loc(1,1)
                do id = 1,size(B_F_loc,1)-1
                    BH_0_ALT(:,id+1) = BH_0_ALT(:,id) + &
                        &B_F_loc(id+1,1)*cos(id*theta_B) + &
                        &B_F_loc(id+1,2)*sin(id*theta_B)
                end do
                call print_ex_2D(['orig BH','alt BH '],'',&
                    &reshape([BH_0(:,kd),BH_0_ALT(:,size(B_F_loc,1))],&
                    &[size(BH_0,1),2]),x=&
                    &reshape([theta_B],[size(BH_0,1),1]))
                
                call writo('Plotting Fourier approximation')
                call print_ex_2D(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                    &'_F_series',BH_0_ALT,&
                    &x=reshape([theta_B],[size(BH_0,1),1]))
                call draw_ex(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                    &'_F_series',size(B_F_loc,1),1,.false.)
                
                call writo('Making animation')
                call print_ex_2D(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                    &'_F_series_anim',BH_0_ALT,&
                    &x=reshape([theta_B],[size(BH_0,1),1]))
                call draw_ex(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                    &'_F_series_anim',size(B_F_loc,1),1,.false.,&
                    &is_animated=.true.)
                
                deallocate(BH_0_ALT)
            end if
#endif
            deallocate(B_F_loc)
            
            call lvl_ud(-1)
        end do
        
        ! plot axisymetric boundary in 3D
        call plot_boundary(B_F(1:1,:,:,:),[0],'last_fs',plot_dim,plot_lims)
        
        call lvl_ud(-1)
        
        ! sort and output
        if (pert_eq) then
            ! user output
            call writo('Set up perturbed data')
            call lvl_ud(1)
            
            ! sort
            allocate(n_pert_copy(nr_n))
            allocate(delta_copy(nr_n,0:nr_m_max,2))
            n_pert_copy = n_pert(1:nr_n)
            delta_copy = delta(1:nr_n,:,:)
            deallocate(n_pert); allocate(n_pert(nr_n)); n_pert = n_pert_copy
            deallocate(delta); allocate(delta(nr_n,0:nr_m_max,2))
            allocate(piv(nr_n))
            call bubble_sort(n_pert,piv)                                        ! sort n
            do jd = 1,nr_n
                delta(jd,:,:) = delta_copy(piv(jd),:,:)
            end do
            deallocate(piv)
            
            ! update the index of N = 0 in B_F
            id_n_0 = minloc(abs(n_pert),1)
            if (id_n_0.gt.1) then
                B_F(id_n_0,:,:,:) = B_F(1,:,:,:)
                B_F(1,:,:,:) = 0._dp
            end if
            
            ! output
            call writo('Summary of perturbation form:')
            call lvl_ud(1)
            do jd = 1,nr_n
                do id = 0,nr_m_max
                    if (maxval(abs(delta(jd,id,:))).lt.tol_zero) cycle          ! this mode has no amplitude
                    if (n_pert(jd).ge.0) then
                        pm(1) = '-'
                    else
                        pm(1) = '+'
                    end if
                    if (delta(jd,id,1).ge.0) then
                        pm(2) = '+'
                    else
                        pm(2) = '-'
                    end if
                    if (jd.eq.1 .and. id.eq.1) pm(2) = ''
                    if (delta(jd,id,2).ge.0) then
                        pm(3) = '+'
                    else
                        pm(3) = '-'
                    end if
                    call writo(&
                        &pm(2)//' '//trim(r2strt(abs(delta(jd,id,1))))//&
                        &' cos('//trim(i2str(id))//' θ '//pm(1)//' '//&
                        &trim(i2str(abs(n_pert(jd))))//' ζ) '//&
                        &pm(3)//' '//trim(r2strt(abs(delta(jd,id,2))))//&
                        &' sin('//trim(i2str(id))//' θ '//pm(1)//' '//&
                        &trim(i2str(abs(n_pert(jd))))//' ζ)')
                end do
            end do
            call lvl_ud(-1)
        else
            id_n_0 = 1
        end if
        
        ! loop  over  all  toroidal  modes   N  and  include  the  corresponding
        ! perturbation, consisting of terms of the form:
        ! delta_c cos(M theta - N zeta) + delta_s sin(M theta - N zeta) = 
        !   delta_c [cos(M theta)cos(N zeta) + sin(M theta) sin(N zeta) ] + 
        !   delta_s [sin(M theta)cos(N zeta) - cos(M theta) sin(N zeta) ].
        ! First, however, the delta_c and delta_s  have to be translated from an
        ! absolute modification of the radius r = sqrt(R^2 + Z^2) to an absolute
        ! modification of R and Z:
        !    ΔR = delta(m) cos(m theta)
        !    ΔZ = delta(m) sin(m theta),
        ! which simply shifts up and down the index m by one.
        ! Note that  for perturbation type 2  (2D map), there is  an extra shift
        ! before.
        pert_N: do jd = 1,nr_n
            if (pert_eq) then
                ! initialize
                allocate(B_F_loc(0:nr_m_max,2))
                
                ! shift due to proportionality map
                if (pert_style.eq.2) then
                    call shift_F(delta(jd,:,:),prop_B_tor_interp_F,&
                        &B_F_loc(:,:))
                    delta(jd,:,:) = B_F_loc
                end if
                
                ! for R, shift due to cos(theta)
                loc_F(0,:) = [0._dp,0._dp]
                loc_F(1,:) = [1._dp,0._dp]
                call shift_F(delta(jd,:,:),loc_F,B_F_loc(:,:))
                B_F(jd,:,:,1) = B_F(jd,:,:,1) + B_F_loc
                
                ! for Z, shift due to sin(theta)
                loc_F(0,:) = [0._dp,0._dp]
                loc_F(1,:) = [0._dp,1._dp]
                call shift_F(delta(jd,:,:),loc_F,B_F_loc(:,:))
                B_F(jd,:,:,2) = B_F(jd,:,:,2) + B_F_loc
                
                deallocate(B_F_loc)
            end if
        end do pert_N
        
        ! plot boundary in 3D
        if (pert_eq) then
            call plot_boundary(B_F,n_pert(1:nr_n),'last_fs_pert',&
                &plot_dim,plot_lims)
            call lvl_ud(-1)
        end if
        
        ! user output
        file_name = "input."//trim(eq_name)
        if (pert_eq) then
            select case (pert_type)
                case (1)                                                        ! Fourier modes in geometrical coordinates
                    file_name = trim(file_name)//'_'//&
                        &trim(HEL_pert_file_name)
                case (2)                                                        ! specified manually
                    do jd = 1,nr_n
                        file_name = trim(file_name)//'_N'//&
                            &trim(i2str(n_pert(jd)))
                        do kd = 0,nr_m_max
                            file_name = trim(file_name)//'M'//&
                                &trim(i2str(kd))
                        end do
                    end do
                case (3)                                                        ! 2-D map
                    file_name = trim(file_name)//'_2DMAP'
            end select
        end if
        call writo('Generate VMEC input file "'//trim(file_name)//'"')
        call writo('This can be used for VMEC porting')
        call lvl_ud(1)
        
        ! set up nfp
        if (pert_eq) then
            ! save greatest common denominator of N into nfp, excluding N=0
            do jd = 2,nr_n
                if (jd.eq.2) then
                    nfp = n_pert(jd)
                else
                    nfp = GCD(nfp,n_pert(jd))
                end if
            end do
            if (nfp.ne.0) then
                n_pert = n_pert/nfp
            else
                nfp = 1
            end if
        else
            nfp = 1
        end if
        
        ! detect whether there is stellarator symmetry
        stel_sym = .false.
        if (maxval(abs(B_F(:,:,2,1)))/maxval(abs(B_F(:,:,1,1))) &
            &.lt. m_tol .and. &                                                 ! R_s << R_c
            maxval(abs(B_F(:,:,1,2)))/maxval(abs(B_F(:,:,2,2))) &
            &.lt. m_tol) &                                                      ! R_c << Z_s
            &stel_sym = .true.
        if (stel_sym) then
            call writo("The equilibrium configuration has stellarator &
                &symmetry")
        else
            call writo("The equilibrium configuration does not have &
                &stellarator symmetry")
        end if
        
        ! find out how many poloidal modes would be necessary
        rec_min_m = 1
        norm_B_H = maxval(abs(B_F))
        do id = 0,nr_m_max
            if (maxval(abs(B_F(:,id,:,1)/norm_B_H)).gt.m_tol .or. &             ! for R
                &maxval(abs(B_F(:,id,:,2)/norm_B_H)).gt.m_tol) &                ! for Z
                &rec_min_m = id
        end do
        call writo("Detected recommended number of poloidal modes: "//&
            &trim(i2str(rec_min_m)))
        
        ! possibly get maximum number of modes to plot
        max_n_B_output = rec_min_m
        call writo('Do you want to change the maximum number of modes &
            &from current '//trim(i2str(max_n_B_output))//'?')
        change_max_n_B_output = get_log(.false.)
        if (change_max_n_B_output) then
            call writo('Maximum number of modes to output?')
            max_n_B_output = get_int(lim_lo=1)
        end if
        
        ! interpolate for VMEC (V) and Jorek (J) output
        s_VJ = [((kd-1._dp)/(size(s_VJ)-1),kd=1,size(s_VJ))]
        allocate(psi_t(grid_eq%n(3)))
        psi_t = eq_1%flux_t_E(:,0)/eq_1%flux_t_E(grid_eq%n(3),0)
        do id = 0,1
            ierr = spline(psi_t,eq_1%pres_E(:,id),s_VJ,pres_VJ(:,id),&
                &ord=norm_disc_prec_eq)
            CHCKERR('')
        end do
        if (use_normalization) then
            pres_VJ = pres_VJ*pres_0
            pres_VJ(:,1) = pres_VJ(:,1) / psi_0
        end if
        ierr = spline(psi_t,RBphi_H(:,0),s_VJ,F_VJ,ord=norm_disc_prec_eq)
        CHCKERR('')
        if (use_normalization) F_VJ = F_VJ * R_0*B_0
        ierr = spline(psi_t,FFp,s_VJ,FFp_VJ,ord=norm_disc_prec_eq)
        CHCKERR('')
        if (use_normalization) then
            FFp_VJ = FFp_VJ * (R_0*B_0)**2/psi_0
        end if
        ierr = spline(psi_t,flux_p_H(:,0),s_VJ,flux_J,ord=norm_disc_prec_eq)
        CHCKERR('')
        if (use_normalization) flux_J = flux_J*psi_0
        ierr = spline(psi_t,eq_1%q_saf_E(:,0),s_VJ,q_saf_VJ,&
            &ord=norm_disc_prec_eq)
        CHCKERR('')
        ierr = spline(psi_t,-eq_1%rot_t_E(:,0),s_VJ,rot_t_V,&
            &ord=norm_disc_prec_eq)
        CHCKERR('')
        ierr = spline(psi_t,I_tor*norm_transf(:,1),s_VJ,I_tor_V,&
            &ord=norm_disc_prec_eq)
        CHCKERR('')
        ierr = spline(psi_t,I_tor*norm_transf(:,2),s_VJ,I_tor_J,&
            &ord=norm_disc_prec_eq)
        CHCKERR('')
        if (use_normalization) I_tor_V = I_tor_V * R_0*B_0/mu_0_original
        if (use_normalization) I_tor_J = I_tor_J * R_0*B_0/mu_0_original
        
        ! output to VMEC input file
        open(HEL_export_i,STATUS='replace',FILE=trim(file_name),&
            &IOSTAT=ierr)
        CHCKERR('Failed to open file')
        
        write(HEL_export_i,"(A)") "!----- General Parameters -----"
        write(HEL_export_i,"(A)") "&INDATA"
        write(HEL_export_i,"(A)") "MGRID_FILE = 'NONE',"
        write(HEL_export_i,"(A)") "PRECON_TYPE = 'NONE'"
        write(HEL_export_i,"(A)") "PREC2D_THRESHOLD = 3.E-8"
        write(HEL_export_i,"(A)") "DELT = 1.00E+00,"
        write(HEL_export_i,"(A)") "NS_ARRAY = 19, 39, 79, 159, 319"
        write(HEL_export_i,"(A)") "LRFP = F"
        write(HEL_export_i,"(A,L1)") "LASYM = ", .not.stel_sym
        write(HEL_export_i,"(A)") "LFREEB = F"
        write(HEL_export_i,"(A)") "NTOR = "//trim(i2str(nr_n-1))                ! -NTOR .. NTOR
        write(HEL_export_i,"(A)") "MPOL = "//trim(i2str(rec_min_m))             ! 0 .. MPOL-1
        write(HEL_export_i,"(A)") "TCON0 = 1"
        write(HEL_export_i,"(A)") "FTOL_ARRAY = 1.E-6, 1.E-6, 1.E-10, 1.E-14, &
            &2.000E-18, "
        write(HEL_export_i,"(A)") "NITER = 100000,"
        write(HEL_export_i,"(A)") "NSTEP = 200,"
        write(HEL_export_i,"(A)") "NFP = "//trim(i2str(nfp))
        if (use_normalization) then
            write(HEL_export_i,"(A)") "PHIEDGE = "//&
                &trim(r2str(-eq_1%flux_t_E(grid_eq%n(3),0)*psi_0))
        else
            write(HEL_export_i,"(A)") "PHIEDGE = "//&
                &trim(r2str(-eq_1%flux_t_E(grid_eq%n(3),0)))
        end if
        write(HEL_export_i,"(A)") "CURTOR = "//&
            &trim(r2str(I_tor_int(size(I_tor_int))))
        
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)") ""
        
        write(HEL_export_i,"(A)") "!----- Pressure Parameters -----"
        write(HEL_export_i,"(A)") "PMASS_TYPE = 'cubic_spline'"
        write(HEL_export_i,"(A)") "GAMMA =  0.00000000000000E+00"
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)",advance="no") "AM_AUX_S ="
        do kd = 1,size(s_VJ)
            write(HEL_export_i,"(A1,ES23.16)") " ", s_VJ(kd)
        end do
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)",advance="no") "AM_AUX_F ="
        do kd = 1,size(pres_VJ,1)
            write(HEL_export_i,"(A1,ES23.16)") " ", pres_VJ(kd,0)
        end do
        
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)") ""
        
        write(HEL_export_i,"(A)") &
            &"!----- Current/Iota Parameters -----"
        write(HEL_export_i,"(A)") "NCURR = "//trim(i2str(ncurr))
        write(HEL_export_i,"(A)") ""
        
        ! iota (is used if ncur = 0)
        write(HEL_export_i,"(A)") "PIOTA_TYPE = 'Cubic_spline'"
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)",advance="no") "AI_AUX_S ="
        do kd = 1,size(s_VJ)
            write(HEL_export_i,"(A1,ES23.16)") " ", s_VJ(kd)
        end do
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)",advance="no") "AI_AUX_F ="
        do kd = 1,size(rot_t_V)
            write(HEL_export_i,"(A1,ES23.16)") " ", rot_t_V(kd)
        end do
        write(HEL_export_i,"(A)") ""
        
        ! I_tor (is used if ncur = 1)
        write(HEL_export_i,"(A)") "PCURR_TYPE = 'cubic_spline_Ip'"
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)",advance="no") "AC_AUX_S ="
        do kd = 1,size(s_VJ)
            write(HEL_export_i,"(A1,ES23.16)") " ", s_VJ(kd)
        end do
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)",advance="no") "AC_AUX_F ="
        do kd = 1,size(rot_t_V)
            write(HEL_export_i,"(A1,ES23.16)") " ", I_tor_V(kd)
        end do
        
        write(HEL_export_i,"(A)") ""
        write(HEL_export_i,"(A)") ""
        
        write(HEL_export_i,"(A)") &
            &"!----- Boundary Shape Parameters -----"
        ierr = print_mode_numbers(HEL_export_i,B_F(id_n_0,:,:,:),0,&
            &max_n_B_output)
        CHCKERR('')
        if (pert_eq) then
            do jd = 1,nr_n
                if (jd.eq.id_n_0) cycle                                         ! skip n = 0 as it is already written
                ierr = print_mode_numbers(HEL_export_i,&
                    &B_F(jd,:,:,:),n_pert(jd),max_n_B_output)
                CHCKERR('')
            end do
        end if
        write(HEL_export_i,"(A,ES23.16)") "RAXIS = ", RZ_B_0(1)
        write(HEL_export_i,"(A,ES23.16)") "ZAXIS = ", RZ_B_0(2)
        write(HEL_export_i,"(A)") "&END"
        
        close(HEL_export_i)
        
        ! output to output file for free-boundary coil fitting in JOREK
        file_name = "flux_and_boundary_quantities_"//trim(eq_name)//'.jorek'
        open(HEL_export_i,STATUS='replace',FILE=trim(file_name),&
            &IOSTAT=ierr)
        CHCKERR('Failed to open file')
        
        ! The normal derivatives in FFp are transformed for JOREK output.
        !   - The difference between the HELENA coordinate system (used in the E
        !   quantities here), and  the JOREK coordinate system is  the fact that
        !   HELENA uses  psi_pol/2pi as  the norma  coordinate while  JOREK uses
        !   psi_pol/psi_pol(s). This step therefore introduces a factor
        !       d/dpsi_J = psi_pol(s)/(2*pi)  d/dpsi_H
        ! However, JOREK  then asks for  a quantity  that has this  exact factor
        ! multiplied back into  it, so that the end results  does not require to
        ! be transformed back. The following lines are therefore commented out.
        !!FFp_VJ = FFp_VJ * max_flux_E/(2*pi)
        !!if (use_normalization) FFp_VJ = FFp_VJ * psi_0                          ! to scale max_flux_E
        if (use_normalization) then
            norm_J = [R_0*RMtoG_H,RZ_B_0(2),R_0*RMtoG_H*B_0*BMtoG_H]
        else
            norm_J = [1._dp,1._dp,1._dp]
        end if
        write(HEL_export_i,"('! ',A)") 'normalization factors:'
        write(HEL_export_i,"('R_geo = ',ES23.16,' ! [m]')") norm_J(1)
        write(HEL_export_i,"('Z_geo = ',ES23.16,' ! [m]')") norm_J(2)
        write(HEL_export_i,"('F0 = ',ES23.16,' ! [Tm]')") norm_J(3)
        write(HEL_export_i,*) ''
        write(HEL_export_i,"('! ',7(A23,' '))") 'pol. flux [Tm^2]', &
            &'normalized pol. flux []', 'mu_0 pressure [T^2]', 'F [Tm]', &
            &'-FFp [T^2 m^2/(Wb/rad)]', 'safety factor []', 'I_tor [A]'
        do kd = 1,size(flux_J)
            write(HEL_export_i,"('  ',7(ES23.16,' '))") flux_J(kd), &
                &flux_J(kd)/flux_J(size(flux_J)), mu_0_original*pres_VJ(kd,0), &
                &F_VJ(kd), -FFp_VJ(kd), q_saf_VJ(kd), I_tor_J(kd)
        end do
        write(HEL_export_i,*) ''
        write(HEL_export_i,*) '! R [m] and Z [m] of boundary at psi [ ]'
        do kd = 1,n_B
            write(HEL_export_i,"('R_boundary(',I4,') = ',ES23.16,&
                &', Z_boundary(',I4,') = ',ES23.16,&
                &', psi_boundary(',I4,') = ',ES23.16)") &
                &kd, BH_0(kd,1), kd, BH_0(kd,2), kd, 1._dp
        end do
        
        close(HEL_export_i)
        
        ! output different files
        file_name = "jorek_ffprime"
        open(HEL_export_i,STATUS='replace',FILE=trim(file_name),&
            &IOSTAT=ierr)
        CHCKERR('Failed to open file')
        !write(HEL_export_i,"('! ',2(A23,' '))") 'psi_norm [ ]', '-FFp [T]'
        do kd = 1,size(flux_J)
            write(HEL_export_i,"('  ',2(ES23.16,' '))") &
                &flux_J(kd)/flux_J(size(flux_J)), -FFp_VJ(kd)
        end do
        close(HEL_export_i)
        file_name = "jorek_temperature"
        open(HEL_export_i,STATUS='replace',FILE=trim(file_name),&
            &IOSTAT=ierr)
        CHCKERR('Failed to open file')
        !write(HEL_export_i,"('! ',2(A23,' '))") 'psi_norm [ ]', &
            !&'mu_0 pressure [T^2]'
        do kd = 1,size(flux_J)
            write(HEL_export_i,"('  ',2(ES23.16,' '))") &
                &flux_J(kd)/flux_J(size(flux_J)), mu_0_original*pres_VJ(kd,0)
        end do
        close(HEL_export_i)
        
        ! user output
        call lvl_ud(-1)
        
        call writo('Done')
        
        ! wait
        call pause_prog
    contains
        ! plots the boundary of a toroidal configuration
        !> \private
        subroutine plot_boundary(B,n,plot_name,plot_dim,plot_lims)
            ! input / output
            real(dp), intent(in) :: B(1:,0:,1:,1:)                              ! cosine and sine of fourier series (tor modes, pol modes, cos/sin (m theta_geo), R/Z)
            integer, intent(in) :: n(:)                                         ! toroidal mode numbers
            character(len=*), intent(in) :: plot_name                           ! name of plot
            integer, intent(in) :: plot_dim(2)                                  ! plot dimensions in geometrical angle
            real(dp), intent(in) :: plot_lims(2,2)                              ! limits of plot dimensions [pi]
            
            ! local variables
            real(dp), allocatable :: ang_plot(:,:,:)                            ! angles of plot (theta,zeta)
            real(dp), allocatable :: XYZ_plot(:,:,:)                            ! coordinates of boundary (X,Y,Z)
            integer :: kd, id                                                   ! counters
            integer :: m_F                                                      ! number of poloidal modes
            
            ! set variables
            allocate(ang_plot(plot_dim(1),plot_dim(2),2))                       ! theta and zeta (geometrical)
            allocate(XYZ_plot(plot_dim(1),plot_dim(2),3))                       ! X, Y and Z
            XYZ_plot = 0._dp
            m_F = size(B,2)-1
            
            ! create grid
            do id = 1,plot_dim(1)
                ang_plot(id,:,1) = pi * (plot_lims(1,1) + (id-1._dp)/&
                    &(plot_dim(1)-1)*(plot_lims(2,1)-plot_lims(1,1)))
            end do
            do kd = 1,plot_dim(2)
                ang_plot(:,kd,2) = pi * (plot_lims(1,2) + (kd-1._dp)/&
                    &(plot_dim(2)-1)*(plot_lims(2,2)-plot_lims(1,2)))
            end do
            
            ! inverse Fourier transform
            do id = 1,size(n)                                                   ! toroidal modes
                do kd = 0,m_F                                                   ! poloidal modes
                    ! R
                    XYZ_plot(:,:,1) = XYZ_plot(:,:,1) + &
                        &B(id,kd,1,1)*cos(kd*ang_plot(:,:,1)-&
                        &n(id)*ang_plot(:,:,2)) + &
                        &B(id,kd,2,1)*sin(kd*ang_plot(:,:,1)-&
                        &n(id)*ang_plot(:,:,2))
                    ! Z
                    XYZ_plot(:,:,3) = XYZ_plot(:,:,3) + &
                        &B(id,kd,1,2)*cos(kd*ang_plot(:,:,1)-&
                        &n(id)*ang_plot(:,:,2)) + &
                        &B(id,kd,2,2)*sin(kd*ang_plot(:,:,1)-&
                        &n(id)*ang_plot(:,:,2))
                end do
            end do
            
            ! print cross-section (R,Z)
            call print_ex_2D(['cross_section_'//trim(plot_name)],&
                &'cross_section_'//trim(plot_name),XYZ_plot(:,:,3),&
                &x=XYZ_plot(:,:,1),draw=.false.)
            call draw_ex(['cross_section_'//trim(plot_name)],&
                &'cross_section_'//trim(plot_name),plot_dim(2),1,.false.)
            call print_ex_2D(['cross_section_'//trim(plot_name)//'_inv'],&
                &'cross_section_'//trim(plot_name)//'_inv',&
                &transpose(XYZ_plot(:,:,3)),x=transpose(XYZ_plot(:,:,1)),&
                &draw=.false.)
            call draw_ex(['cross_section_'//trim(plot_name)//'_inv'],&
                &'cross_section_'//trim(plot_name)//'_inv',plot_dim(1),1,&
                &.false.)
            
            ! convert to X, Y, Z
            XYZ_plot(:,:,2) = XYZ_plot(:,:,1)*sin(ang_plot(:,:,2))
            XYZ_plot(:,:,1) = XYZ_plot(:,:,1)*cos(ang_plot(:,:,2))
            
            ! print
            call print_ex_3D('plasma boundary',trim(plot_name),XYZ_plot(:,:,3),&
                &x=XYZ_plot(:,:,1),y=XYZ_plot(:,:,2),draw=.false.)
            call draw_ex(['plasma boundary'],trim(plot_name),1,2,.false.)
        end subroutine plot_boundary
        
        ! print the mode numbers
        !> \private
        integer function print_mode_numbers(file_i,B,n_pert,max_n_B_output) &
            &result(ierr)
            character(*), parameter :: rout_name = 'print_mode_numbers'
            
            ! input / output
            integer, intent(in) :: file_i                                       ! file to output flux quantities for VMEC export
            real(dp), intent(in) :: B(:,:,:)                                    ! cosine and sine of fourier series for R and Z
            integer, intent(in) :: n_pert                                       ! toroidal mode number
            integer, intent(in) :: max_n_B_output                               ! maximum number of modes to output
            
            ! local variables
            integer :: kd                                                       ! counter
            integer :: m                                                        ! counter
            character :: var_name                                               ! name of variables (R or Z)
            character(len=2*max_str_ln) :: temp_output_str                      ! temporary output string
            
            ! initialize ierr
            ierr = 0
            
            do kd = 1,2                                                         ! R and Z
                select case (kd)
                    case (1)
                        var_name = 'R'
                    case (2)
                        var_name = 'Z'
                    case default
                        var_name = '?'
                end select
                
                ! output to file
                do m = 0,min(size(B,1)-1,max_n_B_output)-1
                    write(temp_output_str,"(' ',A,'BC(',I4,', ',I4,') = ',&
                        &ES23.16,'   ',A,'BS(',I4,', ',I4,') = ',ES23.16)") &
                        &var_name,n_pert,m,B(m+1,1,kd), &
                        &var_name,n_pert,m,B(m+1,2,kd)
                    write(UNIT=file_i,FMT='(A)',IOSTAT=ierr) &
                        &trim(temp_output_str)
                    CHCKERR('Failed to write')
                end do
                write(UNIT=file_i,FMT="(A)",IOSTAT=ierr) ""
                CHCKERR('Failed to write')
            end do
        end function print_mode_numbers
    end function create_VMEC_input
    
    !> \private flux version
    integer function print_output_eq_1(grid_eq,eq,data_name) result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_eq_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid variables
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium variables
        character(len=*), intent(in) :: data_name                               !< name under which to store
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        type(var_1D_type), allocatable, target :: eq_1D(:)                      ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: eq_1D_loc => null()                       ! local element in eq_1D
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: id                                                           ! counter
        integer :: loc_size                                                     ! local size
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write flux equilibrium variables to output file')
        call lvl_ud(1)
        
        ! trim grids
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set up the 1D equivalents of the equilibrium variables
        allocate(eq_1D(max_dim_var_1D))
        
        ! Set up common variables eq_1D
        id = 1
        
        ! pres_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'pres_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%pres_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%pres_FD,2)-1]
        loc_size = size(eq%pres_FD(norm_id(1):norm_id(2),:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%pres_FD(norm_id(1):norm_id(2),:),[loc_size])
        
        ! q_saf_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'q_saf_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%q_saf_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%q_saf_FD,2)-1]
        loc_size = size(eq%q_saf_FD(norm_id(1):norm_id(2),:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%q_saf_FD(norm_id(1):norm_id(2),:),[loc_size])
        
        ! rot_t_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'rot_t_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%rot_t_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%rot_t_FD,2)-1]
        loc_size = size(eq%rot_t_FD(norm_id(1):norm_id(2),:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%rot_t_FD(norm_id(1):norm_id(2),:),[loc_size])
        
        ! flux_p_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'flux_p_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%flux_p_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%flux_p_FD,2)-1]
        loc_size = size(eq%flux_p_FD(norm_id(1):norm_id(2),:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%flux_p_FD(norm_id(1):norm_id(2),:),[loc_size])
        
        ! flux_t_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'flux_t_FD'
        allocate(eq_1D_loc%tot_i_min(2),eq_1D_loc%tot_i_max(2))
        allocate(eq_1D_loc%loc_i_min(2),eq_1D_loc%loc_i_max(2))
        eq_1D_loc%tot_i_min = [1,0]
        eq_1D_loc%tot_i_max = [grid_trim%n(3),size(eq%flux_t_FD,2)-1]
        eq_1D_loc%loc_i_min = [grid_trim%i_min,0]
        eq_1D_loc%loc_i_max = [grid_trim%i_max,size(eq%flux_t_FD,2)-1]
        loc_size = size(eq%flux_t_FD(norm_id(1):norm_id(2),:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%flux_t_FD(norm_id(1):norm_id(2),:),[loc_size])
        
        ! rho
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'rho'
        allocate(eq_1D_loc%tot_i_min(1),eq_1D_loc%tot_i_max(1))
        allocate(eq_1D_loc%loc_i_min(1),eq_1D_loc%loc_i_max(1))
        eq_1D_loc%tot_i_min = 1
        eq_1D_loc%tot_i_max = grid_trim%n(3)
        eq_1D_loc%loc_i_min = grid_trim%i_min
        eq_1D_loc%loc_i_max = grid_trim%i_max
        loc_size = size(eq%rho(norm_id(1):norm_id(2)))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = eq%rho(norm_id(1):norm_id(2))
        
        
        ! write
        ierr = print_HDF5_arrs(eq_1D(1:id-1),PB3D_name,trim(data_name),&
            &ind_print=.not.grid_trim%divided)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        call dealloc_var_1D(eq_1D)
        nullify(eq_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_eq_1
    !> \private metric version
    integer function print_output_eq_2(grid_eq,eq,data_name,rich_lvl,par_div,&
        &dealloc_vars) result(ierr)
        use num_vars, only: PB3D_name, eq_jobs_lims, eq_job_nr
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_eq_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid variables
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium variables
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to print
        logical, intent(in), optional :: par_div                                !< is a parallely divided grid
        logical, intent(in), optional :: dealloc_vars                           !< deallocate variables on the fly after writing
        
        ! local variables
        integer :: n_tot(3)                                                     ! total n
        integer :: par_id(2)                                                    ! local parallel interval
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        type(var_1D_type), allocatable, target :: eq_1D(:)                      ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: eq_1D_loc => null()                       ! local element in eq_1D
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: id                                                           ! counter
        logical :: par_div_loc                                                  ! local par_div
        integer :: loc_size                                                     ! local size
        logical :: dealloc_vars_loc                                             ! local dealloc_vars
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write metric equilibrium variables to output file')
        call lvl_ud(1)
        
        ! trim grids
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set local par_div
        par_div_loc = .false.
        if (present(par_div)) par_div_loc = par_div
        
        ! set total n and parallel interval
        n_tot = grid_trim%n
        par_id = [1,n_tot(1)]
        if (grid_trim%n(1).gt.0 .and. par_div_loc) then                         ! total grid includes all equilibrium jobs
            n_tot(1) = maxval(eq_jobs_lims)-minval(eq_jobs_lims)+1
            par_id = eq_jobs_lims(:,eq_job_nr)
        end if
        
        ! set up local dealloc_vars
        dealloc_vars_loc = .false.
        if (present(dealloc_vars)) dealloc_vars_loc = dealloc_vars
        
        ! set up the 1D equivalents of the equilibrium variables
        allocate(eq_1D(max_dim_var_1D))
        
        ! Set up common variables eq_1D
        id = 1
        
        ! g_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'g_FD'
        allocate(eq_1D_loc%tot_i_min(7),eq_1D_loc%tot_i_max(7))
        allocate(eq_1D_loc%loc_i_min(7),eq_1D_loc%loc_i_max(7))
        eq_1D_loc%tot_i_min = [1,1,1,1,0,0,0]
        eq_1D_loc%tot_i_max = [n_tot,6,size(eq%g_FD,5)-1,size(eq%g_FD,6)-1,&
            &size(eq%g_FD,7)-1]
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,1,0,0,0]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,6,&
            &size(eq%g_FD,5)-1,size(eq%g_FD,6)-1,size(eq%g_FD,7)-1]
        loc_size = size(eq%g_FD(:,:,norm_id(1):norm_id(2),:,:,:,:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%g_FD(:,:,norm_id(1):norm_id(2),:,:,:,:),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%g_FD)
        
        ! h_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'h_FD'
        allocate(eq_1D_loc%tot_i_min(7),eq_1D_loc%tot_i_max(7))
        allocate(eq_1D_loc%loc_i_min(7),eq_1D_loc%loc_i_max(7))
        eq_1D_loc%tot_i_min = [1,1,1,1,0,0,0]
        eq_1D_loc%tot_i_max = [n_tot,6,size(eq%h_FD,5)-1,size(eq%h_FD,6)-1,&
            &size(eq%h_FD,7)-1]
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,1,0,0,0]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,6,&
            &size(eq%h_FD,5)-1,size(eq%h_FD,6)-1,size(eq%h_FD,7)-1]
        loc_size = size(eq%h_FD(:,:,norm_id(1):norm_id(2),:,:,:,:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%h_FD(:,:,norm_id(1):norm_id(2),:,:,:,:),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%h_FD)
        
        ! jac_FD
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'jac_FD'
        allocate(eq_1D_loc%tot_i_min(6),eq_1D_loc%tot_i_max(6))
        allocate(eq_1D_loc%loc_i_min(6),eq_1D_loc%loc_i_max(6))
        eq_1D_loc%tot_i_min = [1,1,1,0,0,0]
        eq_1D_loc%tot_i_max = [n_tot,size(eq%jac_FD,4)-1,size(eq%jac_FD,5)-1,&
            &size(eq%jac_FD,6)-1]
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,0,0,0]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &size(eq%jac_FD,4)-1,size(eq%jac_FD,5)-1,size(eq%jac_FD,6)-1]
        loc_size = size(eq%jac_FD(:,:,norm_id(1):norm_id(2),:,:,:))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%jac_FD(:,:,norm_id(1):norm_id(2),:,:,:),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%jac_FD)
        
        ! S
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'S'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = n_tot
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
        loc_size = size(eq%S(:,:,norm_id(1):norm_id(2)))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%S(:,:,norm_id(1):norm_id(2)),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%S)
        
        ! kappa_n
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'kappa_n'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = n_tot
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
        loc_size = size(eq%kappa_n(:,:,norm_id(1):norm_id(2)))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%kappa_n(:,:,norm_id(1):norm_id(2)),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%kappa_n)
        
        ! kappa_g
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'kappa_g'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = n_tot
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
        loc_size = size(eq%kappa_g(:,:,norm_id(1):norm_id(2)))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%kappa_g(:,:,norm_id(1):norm_id(2)),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%kappa_g)
        
        ! sigma
        eq_1D_loc => eq_1D(id); id = id+1
        eq_1D_loc%var_name = 'sigma'
        allocate(eq_1D_loc%tot_i_min(3),eq_1D_loc%tot_i_max(3))
        allocate(eq_1D_loc%loc_i_min(3),eq_1D_loc%loc_i_max(3))
        eq_1D_loc%tot_i_min = [1,1,1]
        eq_1D_loc%tot_i_max = n_tot
        eq_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max]
        loc_size = size(eq%sigma(:,:,norm_id(1):norm_id(2)))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%sigma(:,:,norm_id(1):norm_id(2)),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%sigma)
        
        ! write
        ierr = print_HDF5_arrs(eq_1D(1:id-1),PB3D_name,trim(data_name),&
            &rich_lvl=rich_lvl,ind_print=.not.grid_trim%divided)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        call dealloc_var_1D(eq_1D)
        nullify(eq_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_eq_2
    
    !> \private flux version
    integer function redistribute_output_eq_1(grid,grid_out,eq,eq_out) &
        &result(ierr)
        use MPI_utilities, only: redistribute_var
        
        character(*), parameter :: rout_name = 'redistribute_output_eq_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid variables
        type(grid_type), intent(in) :: grid_out                                 !< redistributed equilibrium grid variables
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium variables
        type(eq_1_type), intent(inout) :: eq_out                                !< flux equilibrium variables in redistributed grid
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Redistribute flux equilibrium variables')
        call lvl_ud(1)
        
        ! test
        if (grid%n(1).ne.grid_out%n(1) .or. grid%n(2).ne.grid_out%n(2)) then
            ierr = 1
            CHCKERR('grid and grid_out are not compatible')
        end if
        
        ! create redistributed flux equilibrium variables
        call eq_out%init(grid_out,setup_E=.false.,setup_F=.true.)
        
        ! for all derivatives
        do id = 0,max_deriv+1
            ! pres_FD
            ierr = redistribute_var(eq%pres_FD(:,id),eq_out%pres_FD(:,id),&
                &[grid%i_min,grid%i_max],[grid_out%i_min,grid_out%i_max])
            CHCKERR('')
            
            ! q_saf_FD
            ierr = redistribute_var(eq%q_saf_FD(:,id),eq_out%q_saf_FD(:,id),&
                &[grid%i_min,grid%i_max],[grid_out%i_min,grid_out%i_max])
            CHCKERR('')
            
            ! rot_t_FD
            ierr = redistribute_var(eq%rot_t_FD(:,id),eq_out%rot_t_FD(:,id),&
                &[grid%i_min,grid%i_max],[grid_out%i_min,grid_out%i_max])
            CHCKERR('')
            
            ! flux_p_FD
            ierr = redistribute_var(eq%flux_p_FD(:,id),eq_out%flux_p_FD(:,id),&
                &[grid%i_min,grid%i_max],[grid_out%i_min,grid_out%i_max])
            CHCKERR('')
            
            ! flux_t_FD
            ierr = redistribute_var(eq%flux_t_FD(:,id),eq_out%flux_t_FD(:,id),&
                &[grid%i_min,grid%i_max],[grid_out%i_min,grid_out%i_max])
            CHCKERR('')
        end do
        
        ! rho
        ierr = redistribute_var(eq%rho,eq_out%rho,&
            &[grid%i_min,grid%i_max],[grid_out%i_min,grid_out%i_max])
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
    end function redistribute_output_eq_1
    !> \private metric version
    integer function redistribute_output_eq_2(grid,grid_out,eq,eq_out) &
        &result(ierr)
        use MPI_utilities, only: redistribute_var
        
        character(*), parameter :: rout_name = 'redistribute_output_eq_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid variables
        type(grid_type), intent(in) :: grid_out                                 !< redistributed equilibrium grid variables
        type(eq_2_type), intent(in) :: eq                                       !< metric equilibrium variables
        type(eq_2_type), intent(inout) :: eq_out                                !< metric equilibrium variables in redistributed grid
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: lims(2), lims_dis(2)                                         ! limits and distributed limits, taking into account the angular extent
        integer :: siz(3), siz_dis(3)                                           ! size for geometric part of variable
        real(dp), allocatable :: temp_var(:)                                    ! temporary variable
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Redistribute metric equilibrium variables')
        call lvl_ud(1)
        
        ! test
        if (grid%n(1).ne.grid_out%n(1) .or. grid%n(2).ne.grid_out%n(2)) then
            ierr = 1
            CHCKERR('grid and grid_out are not compatible')
        end if
        
        ! create redistributed metric equilibrium variables
        call eq_out%init(grid_out,setup_E=.false.,setup_F=.true.)
        
        ! set up limits taking into account angular extent and temporary var
        lims(1) = product(grid%n(1:2))*(grid%i_min-1)+1
        lims(2) = product(grid%n(1:2))*grid%i_max
        lims_dis(1) = product(grid%n(1:2))*(grid_out%i_min-1)+1
        lims_dis(2) = product(grid%n(1:2))*grid_out%i_max
        siz = [grid%n(1:2),grid%loc_n_r]
        siz_dis = [grid_out%n(1:2),grid_out%loc_n_r]
        allocate(temp_var(product(siz_dis)))
        
        ! for all derivatives and metric factors
        do jd = 0,max_deriv
            do kd = 0,max_deriv
                do ld = 0,max_deriv
                    do id = 1,size(eq%g_FD,4)
                        ! g_FD
                        ierr = redistribute_var(reshape(&
                            &eq%g_FD(:,:,:,id,jd,kd,ld),[product(siz)]),&
                            &temp_var,lims,lims_dis)
                        CHCKERR('')
                        eq_out%g_FD(:,:,:,id,jd,kd,ld) = &
                            &reshape(temp_var,siz_dis)
                        
                        ! h_FD
                        ierr = redistribute_var(reshape(&
                            &eq%h_FD(:,:,:,id,jd,kd,ld),[product(siz)]),&
                            &temp_var,lims,lims_dis)
                        CHCKERR('')
                        eq_out%h_FD(:,:,:,id,jd,kd,ld) = &
                            &reshape(temp_var,siz_dis)
                    end do
                    ! jac_FD
                    ierr = redistribute_var(reshape(&
                        &eq%jac_FD(:,:,:,jd,kd,ld),[product(siz)]),&
                        &temp_var,lims,lims_dis)
                    CHCKERR('')
                    eq_out%jac_FD(:,:,:,jd,kd,ld) = &
                        &reshape(temp_var,siz_dis)
                end do
            end do
        end do
        
        ! S
        ierr = redistribute_var(reshape(eq%S,[product(siz)]),temp_var,&
            &lims,lims_dis)
        CHCKERR('')
        eq_out%S = reshape(temp_var,siz_dis)
        
        ! kappa_n
        ierr = redistribute_var(reshape(eq%kappa_n,[product(siz)]),temp_var,&
            &lims,lims_dis)
        CHCKERR('')
        eq_out%kappa_n = reshape(temp_var,siz_dis)
        
        ! kappa_g
        ierr = redistribute_var(reshape(eq%kappa_g,[product(siz)]),temp_var,&
            &lims,lims_dis)
        CHCKERR('')
        eq_out%kappa_g = reshape(temp_var,siz_dis)
        
        ! sigma
        ierr = redistribute_var(reshape(eq%sigma,[product(siz)]),temp_var,&
            &lims,lims_dis)
        CHCKERR('')
        eq_out%sigma = reshape(temp_var,siz_dis)
        
        ! user output
        call lvl_ud(-1)
    end function redistribute_output_eq_2
    
    !> \private individual version
    integer function calc_RZL_ind(grid_eq,eq,deriv) result(ierr)
        use VMEC_utilities, only: fourier2real
        use VMEC_vars, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, is_asym_V
        
        character(*), parameter :: rout_name = 'calc_RZL_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(3)                                         !< derivatives
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv+1,'calc_RZL')
        CHCKERR('')
#endif
        
        ! calculate the variables R,Z and their derivatives
        ierr = fourier2real(R_V_c(:,grid_eq%i_min:grid_eq%i_max,deriv(1)),&
            &R_V_s(:,grid_eq%i_min:grid_eq%i_max,deriv(1)),&
            &grid_eq%trigon_factors,eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)),&
            &sym=[.true.,is_asym_V],deriv=[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(Z_V_c(:,grid_eq%i_min:grid_eq%i_max,deriv(1)),&
            &Z_V_s(:,grid_eq%i_min:grid_eq%i_max,deriv(1)),&
            &grid_eq%trigon_factors,eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)),&
            &sym=[is_asym_V,.true.],deriv=[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(L_V_c(:,grid_eq%i_min:grid_eq%i_max,deriv(1)),&
            &L_V_s(:,grid_eq%i_min:grid_eq%i_max,deriv(1)),&
            &grid_eq%trigon_factors,eq%L_E(:,:,:,deriv(1),deriv(2),deriv(3)),&
            &sym=[is_asym_V,.true.],deriv=[deriv(2),deriv(3)])
        CHCKERR('')
    end function calc_RZL_ind
    !> \private array version
    integer function calc_RZL_arr(grid_eq,eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_RZL_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_RZL_ind(grid_eq,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_RZL_arr
    
    !> \private individual version
    integer function calc_g_C_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_g_C_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_C')
        CHCKERR('')
#endif
        
        ! initialize
        eq%g_C(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! calculate
        if (sum(deriv).eq.0) then
            eq%g_C(:,:,:,c([1,1],.true.),deriv(1),deriv(2),deriv(3)) = 1.0_dp
            eq%g_C(:,:,:,c([3,3],.true.),deriv(1),deriv(2),deriv(3)) = 1.0_dp
        end if
        ierr = add_arr_mult(eq%R_E,eq%R_E,&
            &eq%g_C(:,:,:,c([2,2],.true.),deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_g_C_ind
    !> \private array version
    integer function calc_g_C_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_C_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_C_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_C_arr
    
    !> \private individual version
    integer function calc_g_V_ind(eq,deriv) result(ierr)
        use eq_utilities, only: calc_g
        
        character(*), parameter :: rout_name = 'calc_g_V_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_V')
        CHCKERR('')
#endif
        
        ierr = calc_g(eq%g_C,eq%T_VC,eq%g_E,deriv,max_deriv)
        CHCKERR('')
    end function calc_g_V_ind
    !> \private  array version
    integer function calc_g_V_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_V_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_V_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_V_arr
    
    !> \private individual version
    integer function calc_g_H_ind(grid_eq,eq,deriv) result(ierr)
        use num_vars, only: norm_disc_prec_eq
        use HELENA_vars, only: ias, nchi, R_H, Z_H, chi_H, h_H_33
        use num_utilities, only: c, spline
        use EZspline_obj
        use EZspline
        
        character(*), parameter :: rout_name = 'calc_g_H_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium 
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! local variables
        type(EZspline2_r8) :: f_spl(2)                                          ! spline object for interpolation, even and odd
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd, ld                                               ! counters
        integer :: bcs(2,3)                                                     ! boundary conditions (theta(even), theta(odd), r)
        integer :: bc_ld(6)                                                     ! boundary condition type for each metric element
        real(dp) :: bcs_val(2,3)                                                ! values for boundary conditions
        real(dp), allocatable :: Rchi(:,:), Rpsi(:,:)                           ! chi and psi derivatives of R
        real(dp), allocatable :: Zchi(:,:), Zpsi(:,:)                           ! chi and psi derivatives of Z
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_H')
        CHCKERR('')
#endif
        
        ! set up boundary conditions
        if (ias.eq.0) then                                                      ! top-bottom symmetric
            bcs(:,1) = [1,1]                                                    ! theta(even): zero first derivative
            bcs(:,2) = [2,2]                                                    ! theta(odd): zero first derivative
        else
            bcs(:,1) = [-1,-1]                                                  ! theta(even): periodic
            bcs(:,2) = [-1,-1]                                                  ! theta(odd): periodic
        end if
        bcs(:,3) = [0,0]                                                        ! r: not-a-knot
        bcs_val = 0._dp
        
        ! boundary condition type for each metric element
        bc_ld = [1,2,0,1,0,1]                                                   ! 1: even, 2: odd, 0: zero
        
        ! initialize
        eq%g_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        if (sum(deriv).eq.0) then                                               ! no derivatives
            ! initialize variables
            allocate(Rchi(nchi,grid_eq%loc_n_r),Rpsi(nchi,grid_eq%loc_n_r))
            allocate(Zchi(nchi,grid_eq%loc_n_r),Zpsi(nchi,grid_eq%loc_n_r))
            
            ! normal derivatives
            do id = 1,grid_eq%n(1)
                ierr = spline(grid_eq%loc_r_E,&
                    &R_H(id,grid_eq%i_min:grid_eq%i_max),grid_eq%loc_r_E,&
                    &Rpsi(id,:),ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,3),&
                    &bcs_val=bcs_val(:,3))
                CHCKERR('')
                ierr = spline(grid_eq%loc_r_E,&
                    &Z_H(id,grid_eq%i_min:grid_eq%i_max),grid_eq%loc_r_E,&
                    &Zpsi(id,:),ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,3),&
                    &bcs_val=bcs_val(:,3))
                CHCKERR('')
            end do
            
            ! poloidal derivatives
            do kd = 1,grid_eq%loc_n_r
                ierr = spline(chi_H,R_H(:,grid_eq%i_min-1+kd),chi_H,&
                    &Rchi(:,kd),ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,1),&
                    &bcs_val=bcs_val(:,1))                                      ! even
                CHCKERR('')
                ierr = spline(chi_H,Z_H(:,grid_eq%i_min-1+kd),chi_H,&
                    &Zchi(:,kd),ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,2),&
                    &bcs_val=bcs_val(:,2))                                      ! odd
                CHCKERR('')
            end do
            
            ! set up g_H
            do jd = 1,grid_eq%n(2)
                eq%g_E(:,jd,:,c([1,1],.true.),0,0,0) = Rpsi*Rpsi+Zpsi*Zpsi
                eq%g_E(:,jd,:,c([1,2],.true.),0,0,0) = Rpsi*Rchi+Zpsi*Zchi
                eq%g_E(:,jd,:,c([2,2],.true.),0,0,0) = Rpsi*Rchi+Zpsi*Zchi
                eq%g_E(:,jd,:,c([2,2],.true.),0,0,0) = Rchi*Rchi+Zchi*Zchi
                eq%g_E(:,jd,:,c([3,3],.true.),0,0,0) = &
                    &1._dp/h_H_33(:,grid_eq%i_min:grid_eq%i_max)
            end do
            
            ! clean up
            deallocate(Zchi,Rchi,Zpsi,Rpsi)
        else if (deriv(3).ne.0) then                                            ! axisymmetry: deriv. in tor. coord. is zero
            !eq%g_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (deriv(1).le.2 .and. deriv(2).le.2) then
            ! initialize 2-D cubic spline for even (1) and odd (2) quantities
            do ld = 1,2
                call EZspline_init(f_spl(ld),grid_eq%n(1),grid_eq%loc_n_r,&
                    &bcs(:,ld),bcs(:,3),ierr) 
                call EZspline_error(ierr)
                CHCKERR('')
                
                ! set grid
                f_spl(ld)%x1 = chi_H
                f_spl(ld)%x2 = grid_eq%loc_r_E
                
                ! set boundary conditions
                f_spl(ld)%bcval1min = bcs_val(1,ld)
                f_spl(ld)%bcval1max = bcs_val(2,ld)
                f_spl(ld)%bcval2min = bcs_val(1,3)
                f_spl(ld)%bcval2max = bcs_val(2,3)
            end do
            
            do ld = 1,6
                do jd = 1,grid_eq%n(2)
                    if (bc_ld(ld).eq.0) cycle                                   ! quantity is zero
                    
                    ! set up 
                    call EZspline_setup(f_spl(bc_ld(ld)),&
                        &eq%g_E(:,jd,:,ld,0,0,0),ierr,exact_dim=.true.)         ! match exact dimensions, none of them old Fortran bullsh*t!
                    call EZspline_error(ierr)
                    CHCKERR('')
                    
                    ! derivate
                    ! Note:  deriv  is the  other  way  around because  poloidal
                    ! indices come before normal
                    call EZspline_derivative(f_spl(bc_ld(ld)),deriv(2),&
                        &deriv(1),grid_eq%n(1),grid_eq%loc_n_r,chi_H,&
                        &grid_eq%loc_r_E,eq%g_E(:,jd,:,ld,deriv(1),deriv(2),0),&
                        &ierr) 
                    call EZspline_error(ierr)
                    CHCKERR('')
                end do
            end do
            
            ! free
            do ld = 1,2
                call EZspline_free(f_spl(ld),ierr)
                call EZspline_error(ierr)
            end do
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_g_H_ind
    !> \private array version
    integer function calc_g_H_arr(grid_eq,eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_H_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibirum grid
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_H_ind(grid_eq,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_H_arr

    !> \private individual version
    integer function calc_g_F_ind(eq,deriv) result(ierr)
        use eq_utilities, only: calc_g
        
        character(*), parameter :: rout_name = 'calc_g_F_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_F')
        CHCKERR('')
#endif
        
        ! calculate g_F
        ierr = calc_g(eq%g_E,eq%T_FE,eq%g_F,deriv,max_deriv)
        CHCKERR('')
    end function calc_g_F_ind
    !> \private array version
    integer function calc_g_F_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_F_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_F_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_F_arr
    
    !> \private individual version
    integer function calc_jac_V_ind(grid,eq,deriv) result(ierr)
        use VMEC_utilities, only: fourier2real
        use VMEC_vars, only: jac_V_c, jac_V_s, is_asym_V
        
        character(*), parameter :: rout_name = 'calc_jac_V_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< grid for which to calculate Jacobian
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_V')
        CHCKERR('')
#endif
        
        ! initialize
        eq%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate the Jacobian and its derivatives
        ierr = fourier2real(jac_V_c(:,grid%i_min:grid%i_max,deriv(1)),&
            &jac_V_s(:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)),&
            &sym=[.true.,is_asym_V],deriv=[deriv(2),deriv(3)])
        CHCKERR('')
    end function calc_jac_V_ind
    !> \private array version
    integer function calc_jac_V_arr(grid,eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_V_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< grid for which to calculate Jacobian
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_V_ind(grid,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_V_arr
    
    !> \private individual version
    integer function calc_jac_H_ind(grid_eq,eq_1,eq_2,deriv) result(ierr)
        use HELENA_vars, only: ias, h_H_33, RBphi_H, chi_H
        use EZspline_obj
        use EZspline
        
        character(*), parameter :: rout_name = 'calc_jac_H_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(EZspline2_r8) :: f_spl                                             ! spline object for interpolation
        integer :: jd, kd                                                       ! counters
        integer :: bcs(2,2)                                                     ! boundary conditions (theta(even), r)
        real(dp) :: bcs_val(2,2)                                                ! values for boundary conditions
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_jac_H')
        CHCKERR('')
#endif
        
        ! set up boundary conditions
        if (ias.eq.0) then                                                      ! top-bottom symmetric
            bcs(:,1) = [1,1]                                                    ! theta(even): zero first derivative
        else
            bcs(:,1) = [-1,-1]                                                  ! periodic
        end if
        bcs(:,2) = [0,0]                                                        ! r: not-a-knot
        bcs_val = 0._dp
        
        ! calculate determinant
        if (sum(deriv).eq.0) then                                               ! no derivatives
            do kd = 1,grid_eq%loc_n_r
                do jd = 1,grid_eq%n(2)
                    eq_2%jac_E(:,jd,kd,0,0,0) = eq_1%q_saf_E(kd,0)/&
                        &(h_H_33(:,grid_eq%i_min-1+kd)*&
                        &RBphi_H(grid_eq%i_min-1+kd,0))
                end do
            end do
        else if (deriv(3).ne.0) then                                            ! axisymmetry: deriv. in tor. coord. is zero
            eq_2%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (deriv(1).le.2 .and. deriv(2).le.2) then
            ! initialize 2-D cubic spline
            call EZspline_init(f_spl,grid_eq%n(1),grid_eq%loc_n_r,&
                &bcs(:,1),bcs(:,2),ierr) 
            call EZspline_error(ierr)
            CHCKERR('')
            
            
            ! set grid
            f_spl%x1 = chi_H
            f_spl%x2 = grid_eq%loc_r_E
            
            ! set boundary conditions
            f_spl%bcval1min = bcs_val(1,1)
            f_spl%bcval1max = bcs_val(2,1)
            f_spl%bcval2min = bcs_val(1,2)
            f_spl%bcval2max = bcs_val(2,2)
            
            do jd = 1,grid_eq%n(2)
                ! set up 
                call EZspline_setup(f_spl,eq_2%jac_E(:,jd,:,0,0,0),ierr,&
                    &exact_dim=.true.)                                          ! match exact dimensions, none of them old Fortran bullsh*t!
                call EZspline_error(ierr)
                CHCKERR('')
                
                ! derivate
                ! Note: deriv is  the other way around  because poloidal indices
                ! come before normal
                call EZspline_derivative(f_spl,deriv(2),deriv(1),grid_eq%n(1),&
                    &grid_eq%loc_n_r,chi_H,grid_eq%loc_r_E,&
                    &eq_2%jac_E(:,jd,:,deriv(1),deriv(2),0),ierr) 
                call EZspline_error(ierr)
                CHCKERR('')
            end do
            
            ! free
            call EZspline_free(f_spl,ierr)
            call EZspline_error(ierr)
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_jac_H_ind
    !> \private array version
    integer function calc_jac_H_arr(grid_eq,eq_1,eq_2,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_H_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrim grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_H_ind(grid_eq,eq_1,eq_2,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_H_arr
    
    !> \private individual version
    integer function calc_jac_F_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_F_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_F')
        CHCKERR('')
#endif
        
        ! initialize
        eq%jac_F(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate determinant
        ierr = add_arr_mult(eq%jac_E,eq%det_T_FE,&
            &eq%jac_F(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_jac_F_ind
    !> \private array version
    integer function calc_jac_F_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_F_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_F_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_F_arr

    !> \private individual version
    integer function calc_T_VC_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_VC_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VC')
        CHCKERR('')
#endif
        
        ! initialize
        eq%T_VC(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! calculate transformation matrix T_V^C
        eq%T_VC(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1)+1,deriv(2),deriv(3))
        !eq%T_VC(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) = 0
        eq%T_VC(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%Z_E(:,:,:,deriv(1)+1,deriv(2),deriv(3))
        eq%T_VC(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1),deriv(2)+1,deriv(3))
        !eq%T_VC(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = 0
        eq%T_VC(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%Z_E(:,:,:,deriv(1),deriv(2)+1,deriv(3))
        eq%T_VC(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)+1)
        if (sum(deriv).eq.0) then
            eq%T_VC(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 1.0_dp
        else
            !eq%T_VC(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 0
        end if
        eq%T_VC(:,:,:,c([3,3],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)+1)
    end function calc_T_VC_ind
    !> \private array version
    integer function calc_T_VC_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VC_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VC_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VC_arr

    !> \private individual version
    integer function calc_T_VF_ind(grid_eq,eq_1,eq_2,deriv) result(ierr)
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_VF_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:,:)                           ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:,:)                            ! - zeta_F and derivatives
        integer :: dims(3)                                                      ! dimensions
        integer :: c1                                                           ! 2D coordinate in met_type storage convention
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VF')
        CHCKERR('')
#endif
        
        ! set up dims
        dims = [grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r]
        
        ! initialize
        eq_2%T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! set up theta_s
        allocate(theta_s(dims(1),dims(2),dims(3),&
            &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
        theta_s = 0.0_dp
        ! start from theta_E
        theta_s(:,:,:,0,0,0) = grid_eq%theta_E
        theta_s(:,:,:,0,1,0) = 1.0_dp
        ! add the deformation described by lambda
        theta_s = theta_s + &
            &eq_2%L_E(:,:,:,0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1)
            
        if (use_pol_flux_F) then
            ! calculate transformation matrix T_V^F
            ! (1,1)
            ierr = add_arr_mult(theta_s,eq_1%q_saf_E(:,1:),&
                &eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ierr = add_arr_mult(eq_2%L_E(:,:,:,1:,0:,0:),eq_1%q_saf_E(:,0:),&
                &eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%T_EF(:,:,kd,c([1,2],.false.),deriv(1),0,0) = &
                        &eq_1%flux_p_E(kd,deriv(1)+1)/(2*pi)
                end do
            !else
                !eq_2%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (1,3)
            eq_2%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                &eq_2%L_E(:,:,:,deriv(1)+1,deriv(2),deriv(3))
            ! (2,1)
            ierr = add_arr_mult(theta_s(:,:,:,0:,1:,0:),eq_1%q_saf_E,&
                &eq_2%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (2,2)
            !eq_2%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,3)
            eq_2%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                &theta_s(:,:,:,deriv(1),deriv(2)+1,deriv(3))
            ! (3,1)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([3,1],.false.),0,0,0) = -1.0_dp
            end if
            ierr = add_arr_mult(eq_2%L_E(:,:,:,0:,0:,1:),eq_1%q_saf_E,&
                &eq_2%T_EF(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (3,2)
            !eq_2%T_EF(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (3,3)
            eq_2%T_EF(:,:,:,c([3,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                &eq_2%L_E(:,:,:,deriv(1),deriv(2),deriv(3)+1)
            
            ! determinant
            eq_2%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ierr = add_arr_mult(-theta_s(:,:,:,0:,1:,0:),&
                &eq_1%flux_p_E(:,1:)/(2*pi),&
                &eq_2%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
        else
            ! set up zeta_s
            allocate(zeta_s(dims(1),dims(2),dims(3),0:deriv(1)+1,&
                &0:deriv(2)+1,0:deriv(3)+1))
            zeta_s = 0.0_dp
            ! start from zeta_E
            zeta_s(:,:,:,0,0,0) = grid_eq%zeta_E
            zeta_s(:,:,:,0,0,1) = 1.0_dp
            
            ! calculate transformation matrix T_V^F
            ! (1,1)
            eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                &-theta_s(:,:,:,deriv(1)+1,deriv(2),deriv(3))
            ierr = add_arr_mult(zeta_s,eq_1%rot_t_E(:,1:),&
                &eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%T_EF(:,:,kd,c([1,2],.false.),deriv(1),0,0) = &
                        &-eq_1%flux_t_E(kd,deriv(1)+1)/(2*pi)
                end do
            !else
                !eq_2%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) 7
                    !&= 0.0_dp
            end if
            ! (1,3)
            !eq_2%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,1)
            eq_2%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                &-theta_s(:,:,:,deriv(1),deriv(2)+1,deriv(3))
            ! (2,2)
            !eq_2%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,3)
            !eq_2%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (3,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%T_EF(:,:,kd,c([3,1],.false.),deriv(1),0,0) = &
                        &eq_1%rot_t_E(kd,deriv(1))
                end do
            end if
            c1 = c([3,1],.false.)                                               ! to avoid compiler crash (as of 26-02-2015)
            eq_2%T_EF(:,:,:,c1,deriv(1),deriv(2),deriv(3)) = &
                &eq_2%T_EF(:,:,:,c1,deriv(1),deriv(2),deriv(3)) &
                &-theta_s(:,:,:,deriv(1),deriv(2),deriv(3)+1)
            ! (3,2)
            !eq_2%T_EF(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (3,3)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([3,3],.false.),0,0,0) = -1.0_dp
            !else
                !eq_2%T_EF(:,:,:,c([3,3],.false.)deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            
            ! determinant
            eq_2%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ierr = add_arr_mult(theta_s(:,:,:,0:,1:,0:),&
                &eq_1%flux_t_E(:,1:)/(2*pi),&
                &eq_2%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
        end if
    end function calc_T_VF_ind
    !> \private array version
    integer function calc_T_VF_arr(grid_eq,eq_1,eq_2,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VF_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VF_ind(grid_eq,eq_1,eq_2,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VF_arr
    
    !> \private individual version
    integer function calc_T_HF_ind(grid_eq,eq_1,eq_2,deriv) result(ierr)
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_HF_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium
        integer, intent(in) :: deriv(:)                                         !< derivatives
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:,:)                           ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:,:)                            ! - zeta_F and derivatives
        integer :: dims(3)                                                      ! dimensions
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_HF')
        CHCKERR('')
#endif
        
        ! set up dims
        dims = [grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r]
        
        ! initialize
        eq_2%T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! set up theta_s
        allocate(theta_s(dims(1),dims(2),dims(3),&
            &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
        theta_s = 0.0_dp
        theta_s(:,:,:,0,0,0) = grid_eq%theta_E
        theta_s(:,:,:,0,1,0) = 1.0_dp
            
        if (use_pol_flux_F) then
            ! calculate transformation matrix T_H^F
            ! (1,1)
            ierr = add_arr_mult(theta_s,-eq_1%q_saf_E(:,1:),&
                &eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([1,2],.false.),0,0,0) = 1._dp
            !else
                !eq_2%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (1,3)
            !eq_2%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%T_EF(:,:,kd,c([2,1],.false.),deriv(1),0,0) = &
                        &-eq_1%q_saf_E(kd,deriv(1))
                end do
            !else
                !eq_2%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (2,2)
            !eq_2%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,3)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([2,3],.false.),0,0,0) = 1.0_dp
            !else
                !eq_2%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (3,1)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([3,1],.false.),0,0,0) = 1.0_dp
            !else
                !eq_2%T_EF(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (3,2)
            !eq_2%T_EF(:,:,:,c([3,2],.false),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (3,3)
            !eq_2%T_EF(:,:,:,c([3,3],.false),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            
            ! determinant
            eq_2%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            if (sum(deriv).eq.0) then
                eq_2%det_T_EF(:,:,:,deriv(1),0,0) = 1._dp
            !else
                !eq_2%det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
        else
            ! set up zeta_s
            allocate(zeta_s(dims(1),dims(2),dims(3),&
                &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
            zeta_s = 0.0_dp
            ! start from zeta_E
            zeta_s(:,:,:,0,0,0) = grid_eq%zeta_E
            zeta_s(:,:,:,0,0,1) = 1.0_dp
            
            ! calculate transformation matrix T_H^F
            ! (1,1)
            ierr = add_arr_mult(zeta_s,eq_1%rot_t_E(:,1:),&
                &eq_2%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%T_EF(:,:,kd,c([1,2],.false.),deriv(1),0,0) = &
                        &eq_1%q_saf_E(kd,deriv(1))
                end do
            !else
                !eq_2%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (1,3)
            !eq_2%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,1)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([2,1],.false.),0,0,0) = -1.0_dp
            !else
                !eq_2%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (2,2)
            !eq_2%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (2,3)
            !eq_2%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (3,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%T_EF(:,:,kd,c([3,1],.false.),deriv(1),0,0) = &
                        &eq_1%rot_t_E(kd,deriv(1))
                end do
            !else
                !eq_2%T_EF(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            ! (3,2)
            !eq_2%T_EF(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                !&0.0_dp
            ! (3,3)
            if (sum(deriv).eq.0) then
                eq_2%T_EF(:,:,:,c([3,3],.false.),deriv(1),0,0) = 1.0_dp
            !else
                !eq_2%T_EF(:,:,:,c([3,3],.false.),deriv(1),deriv(2),deriv(3)) &
                    !&= 0.0_dp
            end if
            
            ! determinant
            eq_2%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    eq_2%det_T_EF(:,:,kd,deriv(1),0,0) = &
                        &eq_1%q_saf_E(kd,deriv(1))
                end do
            !else
                !eq_2%det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
        end if
    end function calc_T_HF_ind
    !> \private array version
    integer function calc_T_HF_arr(grid_eq,eq_1,eq_2,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_HF_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  !< metric equilibrium
        integer, intent(in) :: deriv(:,:)                                       !< derivatives
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_T_HF_ind(grid_eq,eq_1,eq_2,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_HF_arr

    !> Plots the flux quantities in the solution grid.
    !!
    !!  - safety factor \c q_saf
    !!  - rotational transform \c rot_t
    !!  - pressure \c pres
    !!  - poloidal flux \c flux_p
    !!  - toroidal flux \c flux_t
    !!
    !! \return ierr
    integer function flux_q_plot(grid_eq,eq) result(ierr)
        use num_vars, only: rank, no_plots, use_normalization
        use grid_utilities, only: trim_grid
        use MPI_utilities, only: get_ser_var
        use eq_vars, only: pres_0, psi_0
        
        character(*), parameter :: rout_name = 'flux_q_plot'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< normal grid
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium variables
        
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
        
        ! plot using external program
        ierr = flux_q_plot_ex()
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        
        call lvl_ud(-1)
    contains
        ! plots the flux quantities in HDF5
        !> \private
        integer function flux_q_plot_HDF5() result(ierr)
            use num_vars, only: eq_style
            use output_ops, only: plot_HDF5
            use grid_utilities, only: extend_grid_F, trim_grid, calc_XYZ_grid
            use VMEC_utilities, only: calc_trigon_factors
            
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
            ierr = extend_grid_F(grid_trim,grid_plot,grid_eq=grid_eq)
            CHCKERR('')
            
            ! if VMEC, calculate trigonometric factors of plot grid
            if (eq_style.eq.1) then
                ierr = calc_trigon_factors(grid_plot%theta_E,grid_plot%zeta_E,&
                    &grid_plot%trigon_factors)
                CHCKERR('')
            end if
            
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
            ierr = calc_XYZ_grid(grid_eq,grid_plot,X_plot(:,:,:,1),&
                &Y_plot(:,:,:,1),Z_plot(:,:,:,1))
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
            
            ! rescale if normalized
            if (use_normalization) then
                f_plot(:,:,:,3) = f_plot(:,:,:,3)*pres_0                        ! pressure
                f_plot(:,:,:,4) = f_plot(:,:,:,4)*psi_0                         ! flux_p
                f_plot(:,:,:,5) = f_plot(:,:,:,5)*psi_0                         ! flux_t
            end if
            
            ! print the output using HDF5
            call plot_HDF5(plot_titles,file_name,f_plot,plot_dim,plot_offset,&
                &X_plot,Y_plot,Z_plot,col=1,descr='Flux quantities')
            
            ! deallocate and destroy grid
            deallocate(Y_plot_2D)
            deallocate(X_plot,Y_plot,Z_plot,f_plot)
            call grid_plot%dealloc()
        end function flux_q_plot_HDF5
        
        ! plots the pressure and fluxes in external program
        !> \private
        integer function flux_q_plot_ex() result(ierr)
            use eq_vars, only: max_flux_F
            
            character(*), parameter :: rout_name = 'flux_q_plot_ex'
            
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
            !ierr = get_ser_var(eq%q_saf_FD(norm_id(1):norm_id(2),0),ser_var_loc)
            !CHCKERR('')
            !if (rank.eq.0) Y_plot_2D(:,1) = ser_var_loc
            !ierr = get_ser_var(eq%rot_t_FD(norm_id(1):norm_id(2),0),ser_var_loc)
            !CHCKERR('')
            !if (rank.eq.0) Y_plot_2D(:,2) = ser_var_loc
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
            
            ! continue the plot if master
            if (rank.eq.0) then
                ! deallocate local serial variables
                deallocate(ser_var_loc)
                
                ! set up file names
                allocate(file_name(2))
                file_name(1) = 'pres'
                file_name(2) = 'flux'
                
                ! 2D normal variable (normalized F coordinate)
                X_plot_2D(:,1) = grid_trim%r_F*2*pi/max_flux_F
                do id = 2,n_vars
                    X_plot_2D(:,id) = X_plot_2D(:,1)
                end do
                
                ! plot the  individual 2D output  of this process  (except q_saf
                ! and rot_t, as they are already in plot_resonance)
                call print_ex_2D(plot_titles(3),file_name(1),&
                    &Y_plot_2D(:,3),X_plot_2D(:,3),draw=.false.)
                ! fluxes
                call print_ex_2D(plot_titles(4:5),file_name(2),&
                    &Y_plot_2D(:,4:5),X_plot_2D(:,4:5),draw=.false.)
                
                ! draw plot
                call draw_ex(plot_titles(3:3),file_name(1),1,1,.false.)         ! pressure
                call draw_ex(plot_titles(4:5),file_name(2),2,1,.false.)         ! fluxes
                
                ! user output
                call writo('Safety factor and rotational transform are plotted &
                    &using "plot_resonance"')
                
                ! clean up
                deallocate(X_plot_2D,Y_plot_2D)
            end if
        end function flux_q_plot_ex
    end function flux_q_plot
    
    !> Calculates derived equilibrium quantities system.
    !!
    !!  - magnetic shear \c S
    !!  - normal curvature \c kappa_n
    !!  - geodesic curvature \c kappa_g
    !!  - parallel current \c sigma
    !!
    !! This is done using the formula's from \cite Weyens3D
    integer function calc_derived_q(grid_eq,eq_1,eq_2) result(ierr)
        use num_utilities, only: c, spline
        use eq_vars, only: vac_perm
        use num_vars, only: norm_disc_prec_eq, eq_style
        use HELENA_vars, only: RBphi_H, R_H, Z_H, chi_H, q_saf_H, ias, nchi
#if ldebug
        use grid_utilities, only: trim_grid, calc_XYZ_grid
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: calc_int
#endif
        
        character(*), parameter :: rout_name = 'calc_derived_q'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(inout), target :: eq_2                          !< metric equilibrium variables
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: kd_H                                                         ! kd in Helena tables
        integer :: bcs(2,2)                                                     ! boundary conditions (theta(even), theta(odd))
        real(dp) :: bcs_val(2,2)                                                ! values for boundary conditions
        real(dp), allocatable :: Rchi(:,:)                                      ! chi and chi^2 derivatives of R
        real(dp), allocatable :: Zchi(:,:)                                      ! chi and chi^2 derivatives of Z
        real(dp), allocatable :: denom(:)                                       ! dummy denominator
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
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        real(dp), allocatable :: X_plot(:,:,:)                                  ! x values of total plot
        real(dp), allocatable :: Y_plot(:,:,:)                                  ! y values of total plot
        real(dp), allocatable :: Z_plot(:,:,:)                                  ! z values of total plot
        real(dp), pointer :: D13J(:,:,:) => null()                              ! D_alpha,theta jac
        real(dp), pointer :: D23J(:,:,:) => null()                              ! D_psi,theta jac
        real(dp), pointer :: D23g13(:,:,:) => null()                            ! D_psi,theta g_alpha,theta
        real(dp), pointer :: D13g23(:,:,:) => null()                            ! D_alpha,theta g_psi,theta
        real(dp), allocatable :: D3sigma(:,:,:)                                 ! D_theta sigma
        real(dp), allocatable :: D3sigma_ALT(:,:,:)                             ! alternative D_theta sigma
        real(dp), allocatable :: sigma_ALT(:,:,:)                               ! alternative sigma
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle theta_F or zeta_F
        integer :: jd                                                           ! counter
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Set up derived quantities in flux coordinates')
        
        call lvl_ud(1)
        
        ! set up submatrices
        ! jacobian
        J => eq_2%jac_FD(:,:,:,0,0,0)
        D1J => eq_2%jac_FD(:,:,:,1,0,0)
        D2J => eq_2%jac_FD(:,:,:,0,1,0)
        D3J => eq_2%jac_FD(:,:,:,0,0,1)
        ! lower metric factors
        g13 => eq_2%g_FD(:,:,:,c([1,3],.true.),0,0,0)
        D2g13 => eq_2%g_FD(:,:,:,c([1,3],.true.),0,1,0)
        D3g13 => eq_2%g_FD(:,:,:,c([1,3],.true.),0,0,1)
        g23 => eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,0)
        D1g23 => eq_2%g_FD(:,:,:,c([2,3],.true.),1,0,0)
        D3g23 => eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,1)
        g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        D1g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),1,0,0)
        D2g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,1,0)
        D3g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,1)
        ! upper metric factors
        h12 => eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        D3h12 => eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,1)
        h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        D3h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,1)
        h23 => eq_2%h_FD(:,:,:,c([2,3],.true.),0,0,0)
#if ldebug
        D13J => eq_2%jac_FD(:,:,:,1,0,1)
        D23J => eq_2%jac_FD(:,:,:,0,1,1)
        D23g13 => eq_2%g_FD(:,:,:,c([1,3],.true.),0,1,1)
        D13g23 => eq_2%g_FD(:,:,:,c([2,3],.true.),1,0,1)
#endif
        
        ! set up boundary conditions
        if (ias.eq.0) then                                                      ! top-bottom symmetric
            bcs(:,1) = [1,1]                                                    ! theta(even): zero first derivative
            bcs(:,2) = [2,2]                                                    ! theta(odd): zero first derivative
        else
            bcs(:,1) = [-1,-1]                                                  ! theta(even): periodic
            bcs(:,2) = [-1,-1]                                                  ! theta(odd): periodic
        end if
        bcs_val = 0._dp
        
        ! Calculate the shear S
        eq_2%S = -(D3h12/h22 - D3h22*h12/h22**2)/J
        
        ! Calculate the normal curvature kappa_n and geodesic curvature kappa_g
        select case (eq_style)
            case (1)                                                            ! VMEC
                do kd = 1,grid_eq%loc_n_r
                    eq_2%kappa_n(:,:,kd) = vac_perm*&
                        &J(:,:,kd)**2*eq_1%pres_FD(kd,1)/g33(:,:,kd) + &
                        &1._dp/(2*h22(:,:,kd)) * ( &
                        &h12(:,:,kd) * ( D1g33(:,:,kd)/g33(:,:,kd) - &
                        &2*D1J(:,:,kd)/J(:,:,kd) ) + &
                        &h22(:,:,kd) * ( D2g33(:,:,kd)/g33(:,:,kd) - &
                        &2*D2J(:,:,kd)/J(:,:,kd) ) + &
                        &h23(:,:,kd) * ( D3g33(:,:,kd)/g33(:,:,kd) - &
                        &2*D3J(:,:,kd)/J(:,:,kd) ) )
                end do
                eq_2%kappa_g(:,:,:) = (0.5*D1g33/g33 - D1J/J) - &
                    &g13/g33*(0.5*D3g33/g33 - D3J/J)
            case (2)                                                            ! HELENA
                allocate(Rchi(nchi,1:2))
                allocate(Zchi(nchi,1:2))
                allocate(denom(nchi))
                do kd = 1,grid_eq%loc_n_r
                    kd_H = grid_eq%i_min-1+kd
                    do id = 1,2
                        ierr = spline(chi_H,R_H(:,kd_H),chi_H,Rchi(:,id),&
                            &ord=3,deriv=id,bcs=bcs(:,1),bcs_val=bcs_val(:,1))  ! even
                        CHCKERR('')
                        ierr = spline(chi_H,Z_H(:,kd_H),chi_H,Zchi(:,id),&
                            &ord=3,deriv=id,bcs=bcs(:,2),bcs_val=bcs_val(:,2))  ! odd
                        CHCKERR('')
                    end do
                    denom = Rchi(:,1)**2 + Zchi(:,1)**2 + &
                        &(q_saf_H(kd_H,0)*R_H(:,kd_H))**2
                    eq_2%kappa_n(:,1,kd) = &
                        &R_H(:,kd_H)*q_saf_H(kd_H,0)/RBphi_H(kd_H,0) * &
                        &(Zchi(:,1)*Rchi(:,2) - Rchi(:,1)*Zchi(:,2) - &
                        &Zchi(:,1)*q_saf_H(kd_H,0)**2*R_H(:,kd_H)) / &
                        &(denom * (Rchi(:,1)**2 + Zchi(:,1)**2))
                    eq_2%kappa_g(:,1,kd) = &
                        &R_H(:,kd_H)*q_saf_H(kd_H,0) * &
                        &(2*Rchi(:,1) * ( Rchi(:,1)**2 + Zchi(:,1)**2 ) + &
                        &Rchi(:,1) * (q_saf_H(kd_H,0)*R_H(:,kd_H))**2 - &
                        &R_H(:,kd_H) * &
                        &(Rchi(:,1)*Rchi(:,2) + Zchi(:,1)*Zchi(:,2))) / &
                        &(denom**2)
                    if (.not.use_pol_flux_F) then
                        eq_2%kappa_n(:,:,kd) = eq_2%kappa_n(:,:,kd) * &
                            &q_saf_H(kd_H,0)
                        eq_2%kappa_n(:,:,kd) = eq_2%kappa_n(:,:,kd) / &
                            &q_saf_H(kd_H,0)
                    end if
                end do
                do jd = 2,grid_eq%n(2)
                    eq_2%kappa_n(:,jd,:) = eq_2%kappa_n(:,1,:)
                    eq_2%kappa_g(:,jd,:) = eq_2%kappa_g(:,1,:)
                end do
                deallocate(Rchi, Zchi, denom)
        end select
        
        ! Calculate the parallel current sigma
        select case (eq_style)
            case (1)                                                            ! VMEC
                eq_2%sigma = 1._dp/vac_perm*&
                    &(D1g23 - D2g13 - (g23*D1J - g13*D2J)/J)/J
            case (2)                                                            ! HELENA
                do kd = 1,grid_eq%loc_n_r
                    eq_2%sigma(:,:,kd) = -RBphi_H(kd+grid_eq%i_min-1,1)
                end do
        end select
        do kd = 1,grid_eq%loc_n_r
            eq_2%sigma(:,:,kd) = eq_2%sigma(:,:,kd) - &
                &eq_1%pres_FD(kd,1)*J(:,:,kd)*g13(:,:,kd)/g33(:,:,kd)
        end do 
        
#if ldebug
        ! test whether -2 p' J kappa_g = D3sigma and plot kappa components
        if (debug_calc_derived_q) then
            call writo('Testing whether -2 p'' J kappa_g = D3sigma')
            call lvl_ud(1)
            
            ! trim equilibrium grid
            ierr = trim_grid(grid_eq,grid_trim,norm_id)
            CHCKERR('')
            
            ! allocate variables
            allocate(D3sigma(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            allocate(D3sigma_ALT(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            allocate(sigma_ALT(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            allocate(X_plot(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
            allocate(Y_plot(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
            allocate(Z_plot(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
            
            ! calculate grid
            ierr = calc_XYZ_grid(grid_eq,grid_trim,X_plot,Y_plot,Z_plot)
            CHCKERR('')
            
            ! point parallel angle
            if (use_pol_flux_F) then
                ang_par_F => grid_eq%theta_F
            else
                ang_par_F => grid_eq%zeta_F
            end if
            
            ! get derived sigma
            do kd = 1,grid_eq%loc_n_r
                do jd = 1,grid_eq%n(2)
                    ierr = spline(ang_par_F(:,jd,kd),eq_2%sigma(:,jd,kd),&
                        &ang_par_F(:,jd,kd),D3sigma(:,jd,kd),&
                        &ord=norm_disc_prec_eq,deriv=1)
                    CHCKERR('')
                end do
            end do
            
            ! calculate alternatively derived sigma
            do kd = 1,grid_eq%loc_n_r
                D3sigma_ALT(:,:,kd) = &
                    &-2*eq_1%pres_FD(kd,1)*eq_2%kappa_g(:,:,kd)*J(:,:,kd)
            end do
            
            ! calculate alternative sigma by integration
            ! Note:  there  is  an  undetermined constant  of  integration  that
            ! depends on the normal coordinate and the geodesic coordinate.
            do kd = 1,grid_eq%loc_n_r
                do jd = 1,grid_eq%n(2)
                    ierr = calc_int(D3sigma_ALT(:,jd,kd),ang_par_F(:,jd,kd),&
                        &sigma_ALT(:,jd,kd))
                    CHCKERR('')
                end do
                sigma_ALT(:,:,kd) = eq_2%sigma(1,1,kd) + sigma_ALT(:,:,kd)
            end do
            
            ! plot output
            call plot_diff_HDF5(D3sigma,D3sigma_ALT,'TEST_diff_D3sigma',&
                &grid_eq%n,[0,0,grid_eq%i_min-1],&
                &descr='To test whether -2 p'' J kappa_g = D3sigma',&
                &output_message=.true.)
            call plot_diff_HDF5(eq_2%sigma,sigma_ALT,'TEST_diff_sigma',&
                &grid_eq%n,[0,0,grid_eq%i_min-1],descr='To test whether &
                &int(-2 p'' J kappa_g) = sigma',output_message=.true.)
            
            ! plot sigma
            call plot_HDF5('sigma','TEST_sigma',&
                &eq_2%sigma(:,:,norm_id(1):norm_id(2)),&
                &tot_dim=grid_trim%n,loc_offset=[0,0,grid_trim%i_min-1])!,&
                !&X=X_plot,Y=Y_plot,Z=Z_plot)
            
            ! plot shear
            call plot_HDF5('shear','TEST_shear',&
                &eq_2%S(:,:,norm_id(1):norm_id(2)),&
                &tot_dim=grid_trim%n,loc_offset=[0,0,grid_trim%i_min-1])!,&
                !&X=X_plot,Y=Y_plot,Z=Z_plot)
            
            ! plot kappa_n
            call plot_HDF5('kappa_n','TEST_kappa_n',&
                &eq_2%kappa_n(:,:,norm_id(1):norm_id(2)),&
                &tot_dim=grid_trim%n,loc_offset=[0,0,grid_trim%i_min-1])!,&
                !&X=X_plot,Y=Y_plot,Z=Z_plot)
            
            ! plot kappa_g
            call plot_HDF5('kappa_g','TEST_kappa_g',&
                &eq_2%kappa_g(:,:,norm_id(1):norm_id(2)),&
                &tot_dim=grid_trim%n,loc_offset=[0,0,grid_trim%i_min-1])!,&
                !&X=X_plot,Y=Y_plot,Z=Z_plot)
            
            ! clean up
            nullify(ang_par_F)
            
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
    end function calc_derived_q
    
    !> Sets up normalization constants.
    !!
    !!  - VMEC version (\c eq_style=1):\n
    !!      Normalization depends on \c norm_style:
    !!      -# MISHKA
    !!          - R_0:    major radius (= average magnetic axis)
    !!          - B_0:    B on magnetic axis (theta = zeta = 0)
    !!          - pres_0: reference pressure (= B_0^2/mu_0)
    !!          - psi_0:  reference flux (= R_0^2 B_0)
    !!          - rho_0:  reference mass density
    !!      -# COBRA
    !!          - R_0:    major radius (= average geometric axis)
    !!          - pres_0: pressure on magnetic axis
    !!          - B_0:    reference magnetic field (= sqrt(2pres_0mu_0 / beta))
    !!          - psi_0:  reference flux (= R_0^2 B_0 / aspr^2)
    !!          - rho_0:  reference mass density
    !!      where aspr (aspect ratio) and beta are given by VMEC.
    !!  - HELENA version (\c eq_style=2):\n
    !!      Normalization  is  used  by  default  and  does  not  depend  on  \c
    !!      norm_style
    !!
    !! \see read_hel()
    !!
    !! \note \c rho_0  is not given through by the  equilibrium codes and should
    !! be user-supplied
    !!
    !! \return ierr
    integer function calc_normalization_const() result(ierr)
        use num_vars, only: rank, eq_style, mu_0_original, use_normalization, &
            &rich_restart_lvl
        use eq_vars, only: T_0, B_0, pres_0, psi_0, R_0, rho_0
        
        character(*), parameter :: rout_name = 'calc_normalization_const'
        
        ! local variables
        integer :: nr_overridden_const                                          ! nr. of user-overridden constants, to print warning if > 0
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize BR_normalization_provided
        BR_normalization_provided = [.false.,.false.]
        
        if (use_normalization) then
            ! user output
            call writo('Calculating the normalization constants')
            call lvl_ud(1)
            
            ! initialize nr_overridden_const
            nr_overridden_const = 0
            
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
                end select
            end if
            
            ! print warning if user-overridden
            if (nr_overridden_const.gt.0) &
                &call writo(trim(i2str(nr_overridden_const))//&
                &' constants were overridden by user. Consistency is NOT &
                &checked!',warning=.true.)
            
            ! print constants
            call print_normalization_const(R_0,rho_0,B_0,pres_0,psi_0,&
                &mu_0_original,T_0)
            
            ! check whether it is physically consistent
            if (T_0.lt.0._dp) then
                ierr = 1
                err_msg = 'Alfven time is negative. Are you sure the &
                    &equilibrium is consistent?'
                CHCKERR(err_msg)
            end if
            
            ! user output
            call lvl_ud(-1)
            call writo('Normalization constants calculated')
        else if (rich_restart_lvl.eq.1) then                                    ! only for first Richardson leevel
            ! user output
            call writo('Normalization not used')
        end if
    contains 
        ! VMEC version
        !> \private
        subroutine calc_normalization_const_VMEC
            use num_vars, only: norm_style
            use VMEC_vars, only: R_V_c, B_0_V, rmax_surf, rmin_surf, pres_V, &
                &beta_V
            !use VMEC_vars, only: aspr_V
            
            select case (norm_style)
                ! Note that PB3D does not run with the exact COBRA normalization
                ! (for style 2), as it is not a pure nondimensionalization, that
                ! would  modify the  equations by  intrucing extra  factors. The
                ! diffference between the COBRA  normalization used here and the
                ! one used in COBRA is given by:
                !   - B_0 = sqrt(pres_0 mu_0) here
                !       versus
                !     B_0 = sqrt(2 pres_0 mu_0 / beta) in COBRA,
                !   - psi_0 = B_0 R_0^2 here
                !       versus
                !     psi_0 = B_0 (R_0/aspr)^2 in COBRA,
                ! This results in a difference in the factor T_0:
                !   - T_0 = sqrt(rho_0/pres_0) R_0 here
                !       versus
                !     T_0 = sqrt(rho_0/pres_0) R_0 sqrt(beta/2) in COBRA.
                ! Therefore, the Eigenvalues here should be rescaled as
                !   - EV_COBRA = EV_PB3D beta/2. This is done automatically.
                case (1)                                                        ! MISHKA
                    ! user output
                    call writo('Using MISHKA normalization')
                    
                    ! set the  major radius as the  average value of R_V  on the
                    ! magnetic axis
                    if (R_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        R_0 = R_V_c(1,1,0)
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! set B_0 from magnetic field on axis
                    if (B_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        B_0 = B_0_V
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! set reference pres_0 from B_0
                    if (pres_0.ge.huge(1._dp)) then                             ! user did not provide a value
                        pres_0 = B_0**2/mu_0_original
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! set reference flux from R_0 and B_0
                    if (psi_0.ge.huge(1._dp)) then                              ! user did not provide a value
                        psi_0 = R_0**2 * B_0
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                case (2)                                                        ! COBRA
                    ! user output
                    call writo('Using COBRA normalization')
                    
                    ! set the major radius as the average geometric axis
                    if (R_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        R_0 = 0.5_dp*(rmin_surf+rmax_surf)
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! set pres_0 from pressure on axis
                    if (pres_0.ge.huge(1._dp)) then                             ! user did not provide a value
                        pres_0 = pres_V(1,0)
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! set the reference value for B_0 from pres_0 and beta
                    if (B_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        !B_0 = sqrt(2._dp*pres_0*mu_0_original/beta_V)          ! exact COBRA
                        B_0 = sqrt(pres_0*mu_0_original)                        ! pure modified COBRA
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! set reference flux from R_0, B_0 and aspr
                    if (psi_0.ge.huge(1._dp)) then                              ! user did not provide a value
                        !psi_0 = B_0 * (R_0/aspr_V)**2                          ! exact COBRA
                        psi_0 = B_0 * R_0**2                                    ! pure modified COBRA
                    else
                        nr_overridden_const = nr_overridden_const + 1
                    end if
                    
                    ! user output concerning pure modified COBRA
                    call writo('Exact COBRA normalization is substituted by &
                        &"pure", modified version',warning=.true.)
                    call writo('This leaves the equations unmodified',&
                        &alert=.true.)
                    call writo('To translate to exact COBRA, manually multiply &
                        &PB3D Eigenvalues by',alert=.true.)
                    call lvl_ud(1)
                    call writo(trim(r2str(beta_V/2)),alert=.true.)
                    call lvl_ud(-1)
            end select
            
            ! rho_0 is set up through an input variable with the same name
            
            ! set Alfven time
            if (T_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                T_0 = sqrt(mu_0_original*rho_0)*R_0/B_0
            else
                nr_overridden_const = nr_overridden_const + 1
            end if
        end subroutine calc_normalization_const_VMEC
        
        ! HELENA version
        !> \private
        subroutine calc_normalization_const_HEL
            ! user output
            call writo('Using MISHKA normalization')
            
            if (R_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                R_0 = 1._dp
            else
                BR_normalization_provided(1) = .true.
            end if
            if (B_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                B_0 = 1._dp
            else
                BR_normalization_provided(2) = .true.
            end if
            if (pres_0.ge.huge(1._dp)) then                                     ! user did not provide a value
                pres_0 = B_0**2/mu_0_original
            else
                nr_overridden_const = nr_overridden_const + 1
            end if
            if (psi_0.ge.huge(1._dp)) then                                      ! user did not provide a value
                psi_0 = R_0**2 * B_0
            else
                nr_overridden_const = nr_overridden_const + 1
            end if
            if (T_0.ge.huge(1._dp)) T_0 = sqrt(mu_0_original*rho_0)*R_0/B_0     ! only if user did not provide a value
        end subroutine calc_normalization_const_HEL
        
        ! prints the Normalization factors
        !> \private
        subroutine print_normalization_const(R_0,rho_0,B_0,pres_0,psi_0,&
            &mu_0,T_0)
            ! input / output
            real(dp), intent(in), optional :: R_0
            real(dp), intent(in), optional :: rho_0
            real(dp), intent(in), optional :: B_0
            real(dp), intent(in), optional :: pres_0
            real(dp), intent(in), optional :: psi_0
            real(dp), intent(in), optional :: mu_0
            real(dp), intent(in), optional :: T_0
            
            ! user output
            if (present(R_0)) call writo('R_0    = '//trim(r2str(R_0))//' m')
            if (present(rho_0)) call writo('rho_0  = '//trim(r2str(rho_0))//&
                &' kg/m^3')
            if (present(B_0)) call writo('B_0    = '//trim(r2str(B_0))//' T')
            if (present(pres_0)) call writo('pres_0 = '//trim(r2str(pres_0))//&
                &' Pa')
            if (present(psi_0)) call writo('psi_0  = '//trim(r2str(psi_0))//&
                &' Tm^2')
            if (present(mu_0)) call writo('mu_0   = '//&
                &trim(r2str(mu_0))//' Tm/A')
            if (present(T_0)) call writo('T_0    = '//trim(r2str(T_0))//' s')
        end subroutine print_normalization_const
    end function calc_normalization_const
    
    !> Normalize input quantities.
    !!
    !! \see calc_normalization_const()
    subroutine normalize_input()
        use num_vars, only: use_normalization, eq_style, mu_0_original
        use VMEC_ops, only: normalize_VMEC
        use eq_vars, only: vac_perm
        
        ! only normalize if needed
        if (use_normalization) then
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
                    ! other HELENA input already normalized
                    call writo('HELENA input is already normalized with MISHKA &
                        &normalization')
            end select
            
            ! user output
            call lvl_ud(-1)
            call writo('Normalization done')
        end if
    end subroutine normalize_input
    
    !> Plots the  magnetic fields.
    !!
    !! If multiple equilibrium parallel jobs, every  job does its piece, and the
    !! results are joined automatically by plot_HDF5.
    !!
    !! The outputs are  given in contra- and covariant  components and magnitude
    !! in multiple coordinate systems, as indicated in calc_vec_comp().
    !! 
    !! The starting point is the fact that the magnetic field is given by
    !!  \f[\vec{B} = \frac{\vec{e}_{\theta}}{\mathcal{J}}, \f]
    !! in F coordinates. The F covariant components are therefore given by
    !!  \f[B_i = \frac{g_{i3}}{\mathcal{J}} =
    !!      \frac{\vec{e}_i \cdot \vec{e}_3}{\mathcal{J}}, \f]
    !! and the only non-vanishing contravariant component is
    !!  \f[B^3 = \frac{1}{\mathcal{J}}. \f]
    !!
    !! These are then all be transformed to the other coordinate systems.
    !!
    !! \note
    !!  -# Vector plots for different Richardson  levels can be combined to show
    !!  the total grid by just plotting them all individually.
    !!  -# The metric factors and transformation matrices have to be allocated.
    !!
    !! \return ierr
    integer function B_plot(grid_eq,eq_1,eq_2,rich_lvl,plot_fluxes,XYZ) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, use_normalization
        use num_utilities, only: c
        use eq_vars, only: B_0
        
        character(*), parameter :: rout_name = 'B_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        integer, intent(in), optional :: rich_lvl                               !< Richardson level
        logical, intent(in), optional :: plot_fluxes                            !< plot the fluxes
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          !< X, Y and Z of grid
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp), allocatable :: B_com(:,:,:,:,:)                               ! covariant and contravariant components of B (dim1,dim2,dim3,3,2)
        real(dp), allocatable :: B_mag(:,:,:)                                   ! magnitude of B (dim1,dim2,dim3)
        real(dp), allocatable, save :: B_flux_tor(:,:), B_flux_pol(:,:)         ! fluxes
        character(len=10) :: base_name                                          ! base name
        logical :: plot_fluxes_loc                                              ! local plot_fluxes
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (eq_style.eq.1 .and. .not.allocated(grid_eq%trigon_factors)) then
            ierr = 1
            CHCKERR('trigonometric factors not allocated')
        end if
        
        ! set up local plot_fluxes
        plot_fluxes_loc = .false.
        if (present(plot_fluxes)) plot_fluxes_loc = plot_fluxes
        
        ! set up components and magnitude of B
        allocate(B_com(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3,2))
        allocate(B_mag(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        B_com = 0._dp
        B_mag = 0._dp
        do id = 1,3
            B_com(:,:,:,id,1) = eq_2%g_FD(:,:,:,c([3,id],.true.),0,0,0)/&
                &eq_2%jac_FD(:,:,:,0,0,0)
        end do
        B_com(:,:,:,3,2) = 1._dp/eq_2%jac_FD(:,:,:,0,0,0)
        
        ! transform back to unnormalized quantity
        if (use_normalization) B_com = B_com * B_0
        
        ! set plot variables
        base_name = 'B'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! transform coordinates, including the flux
        if (plot_fluxes_loc) then
            ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,B_com,&
                &norm_disc_prec_eq,v_mag=B_mag,base_name=base_name,&
                &v_flux_tor=B_flux_tor,v_flux_pol=B_flux_pol,XYZ=XYZ)
            CHCKERR('')
        else
            ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,B_com,&
                &norm_disc_prec_eq,v_mag=B_mag,base_name=base_name,XYZ=XYZ)
            CHCKERR('')
        end if
    end function B_plot
    
    !> Plots the current.
    !! 
    !! If multiple equilibrium parallel jobs, every  job does its piece, and the
    !! results are joined automatically by plot_HDF5.
    !!
    !! The outputs are  given in contra- and covariant  components and magnitude
    !! in multiple coordinate systems, as indicated in calc_vec_comp().
    !! 
    !! The starting point is the pressure balance
    !!  \f[ \nabla p = \vec{J} \times \vec{B}, \f]
    !! which, using \f$\vec{B} = \frac{\vec{e}_{\theta}}{\mathcal{J}}\f$,
    !! reduces to
    !!  \f[J^\alpha = -p'. \f]
    !! Furthermore, the current has to lie in the magnetic flux surfaces:
    !!  \f[J^\psi = 0. \f]
    !! Finally, the  parallel current \f$\sigma\f$  gives an expression  for the
    !! last contravariant component:
    !!  \f[J^\theta =
    !!      \frac{\sigma}{\mathcal{J}} + p' \frac{B_\alpha}{B_\theta}. \f]
    !! From these, the contravariant components can be calculated as
    !!  \f[J_i = J^\alpha g_{\alpha,i} + J^\theta g_{\theta,i}. \f]
    !! 
    !! These are then all be transformed to the other coordinate systems.
    !! \note
    !!  -# Vector plots for different Richardson  levels can be combined to show
    !!  the total grid by just plotting them all individually.
    !!  -# The metric factors and transformation matrices have to be allocated.
    !!
    !! \return ierr
    integer function J_plot(grid_eq,eq_1,eq_2,rich_lvl,plot_fluxes,XYZ) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, use_normalization
        use num_utilities, only: c
        use eq_vars, only: B_0, R_0, pres_0
#if ldebug
        use eq_vars, only: max_flux_F
        use VMEC_vars, only: J_V_sup_int
        use num_utilities, only: calc_int
#endif
        
        character(*), parameter :: rout_name = 'J_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        integer, intent(in), optional :: rich_lvl                               !< Richardson level
        logical, intent(in), optional :: plot_fluxes                            !< plot the fluxes
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          !< X, Y and Z of grid
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: J_com(:,:,:,:,:)                               ! covariant and contravariant components of J (dim1,dim2,dim3,3,2)
        real(dp), allocatable :: J_mag(:,:,:)                                   ! magnitude of J (dim1,dim2,dim3)
        character(len=10) :: base_name                                          ! base name
        real(dp), allocatable, save :: J_flux_tor(:,:), J_flux_pol(:,:)         ! fluxes
        logical :: plot_fluxes_loc                                              ! local plot_fluxes
        character(len=max_str_ln) :: plot_name                                  ! name of plot
        character(len=max_str_ln) :: plot_titles(2)                             ! titles of plot
#if ldebug
        real(dp), allocatable :: J_V_sup_int2(:,:)                              ! integrated J_V_sup_int
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (eq_style.eq.1 .and. .not.allocated(grid_eq%trigon_factors)) then
            ierr = 1
            CHCKERR('trigonometric factors not allocated')
        end if
        
        ! set up local plot_fluxes
        plot_fluxes_loc = .false.
        if (present(plot_fluxes)) plot_fluxes_loc = plot_fluxes
        
        ! set up components and magnitude of J
        allocate(J_com(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3,2))
        allocate(J_mag(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        J_com = 0._dp
        J_mag = 0._dp
#if ldebug
        if (debug_J_plot) then
            J_com(:,:,:,1,2) = eq_2%g_FD(:,:,:,c([3,3],.true.),0,1,0) - &
                &eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0)*&
                &eq_2%jac_FD(:,:,:,0,1,0)/&
                &eq_2%jac_FD(:,:,:,0,0,0) - &
                &eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,1) + &
                &eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,0)*&
                &eq_2%jac_FD(:,:,:,0,0,1)/&
                &eq_2%jac_FD(:,:,:,0,0,0) 
            J_com(:,:,:,1,2) = J_com(:,:,:,1,2)/eq_2%jac_FD(:,:,:,0,0,0)**2
            J_com(:,:,:,3,2) = eq_2%g_FD(:,:,:,c([2,3],.true.),1,0,0) - &
                &eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,0)*&
                &eq_2%jac_FD(:,:,:,1,0,0)/&
                &eq_2%jac_FD(:,:,:,0,0,0) - &
                &eq_2%g_FD(:,:,:,c([1,3],.true.),0,1,0) + &
                &eq_2%g_FD(:,:,:,c([1,3],.true.),0,0,0)*&
                &eq_2%jac_FD(:,:,:,0,1,0)/&
                &eq_2%jac_FD(:,:,:,0,0,0) 
            J_com(:,:,:,3,2) = J_com(:,:,:,3,2)/eq_2%jac_FD(:,:,:,0,0,0)**2
        else
#endif
            do kd = 1,grid_eq%loc_n_r
                J_com(:,:,kd,1,2) = -eq_1%pres_FD(kd,1)
                J_com(:,:,kd,3,2) = eq_2%sigma(:,:,kd)/&
                    &eq_2%jac_FD(:,:,kd,0,0,0) + eq_1%pres_FD(kd,1)*&
                    &eq_2%g_FD(:,:,kd,c([1,3],.true.),0,0,0) / &
                    &eq_2%g_FD(:,:,kd,c([3,3],.true.),0,0,0)
            end do
#if ldebug
        end if
#endif
        do id = 1,3
            J_com(:,:,:,id,1) = &
                &J_com(:,:,:,1,2)*eq_2%g_FD(:,:,:,c([1,id],.true.),0,0,0) + &
                &J_com(:,:,:,3,2)*eq_2%g_FD(:,:,:,c([3,id],.true.),0,0,0)
        end do
        
        ! transform back to unnormalized quantity
        if (use_normalization) J_com = J_com * pres_0/(R_0*B_0)
        
        ! set plot variables
        base_name = 'J'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! transform coordinates, including the flux
        if (plot_fluxes_loc) then
            ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,J_com,&
                &norm_disc_prec_eq,v_mag=J_mag,base_name=base_name,&
                &v_flux_tor=J_flux_tor,v_flux_pol=J_flux_pol,XYZ=XYZ)
            CHCKERR('')
        else
            ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,J_com,&
                &norm_disc_prec_eq,v_mag=J_mag,base_name=base_name,XYZ=XYZ)
            CHCKERR('')
        end if
        
#if ldebug
        if (eq_style.eq.1) then
            ! transform back to unnormalized quantity
            if (use_normalization) J_V_sup_int = J_V_sup_int * pres_0/(R_0*B_0)
            
            ! plot the fluxes from VMEC
            plot_name = 'J_V_sup_int'
            plot_titles = ['J^theta_V','J^phi_V  ']
            call print_ex_2D(plot_titles,plot_name,J_V_sup_int,&
                &x=reshape(grid_eq%r_F*2*pi/max_flux_F,&
                &[size(grid_eq%r_F),1]),draw=.false.)
            call draw_ex(plot_titles,plot_name,2,1,.false.)
            
            ! integrate
            allocate(J_V_sup_int2(size(J_V_sup_int,1),2))
            do id = 1,2
                ierr = calc_int(J_V_sup_int(:,id),grid_eq%r_E,&
                    &J_V_sup_int2(:,id))
                CHCKERR('')
            end do
            plot_name = 'J_V_sup_int2'
            plot_titles = ['J^theta_V','J^phi_V  ']
            call print_ex_2D(plot_titles,plot_name,J_V_sup_int2,&
                &x=reshape(grid_eq%r_F*2*pi/max_flux_F,&
                &[size(grid_eq%r_F),1]),draw=.false.)
            call draw_ex(plot_titles,plot_name,2,1,.false.)
        end if
#endif
    end function J_plot
    
    !> Plots the curvature.
    !! 
    !! If multiple equilibrium parallel jobs, every  job does its piece, and the
    !! results are joined automatically by plot_HDF5.
    !!
    !! The outputs are  given in contra- and covariant  components and magnitude
    !! in multiple coordinate systems, as indicated in calc_vec_comp().
    !! 
    !! The starting point is the curvature, given by
    !!  \f[\vec{\kappa} =
    !!      \kappa_n \frac{\nabla \psi}{ \left|\nabla \psi\right|^2 } +
    !!      \kappa_g \frac{\nabla \psi \times \vec{B}}{B^2} , \f]
    !! which can be  used to find the covariant and  contravariant components in
    !! Flux coordinates.
    !! 
    !! These are then  transformed to Cartesian coordinates and plotted.
    !!
    !! \note
    !!  -# Vector plots for different Richardson  levels can be combined to show
    !!  the total grid by just plotting them all individually.
    !!  -# The metric factors and transformation matrices have to be allocated.
    !!
    !! \return ierr
    integer function kappa_plot(grid_eq,eq_1,eq_2,rich_lvl,XYZ) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp
        use eq_vars, only: eq_1_type, eq_2_type
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, use_normalization
        use num_utilities, only: c
        use eq_vars, only: R_0
        use VMEC_utilities, only: calc_trigon_factors
#if ldebug
        use num_vars, only: ltest, eq_jobs_lims, eq_job_nr
        use input_utilities, only: get_log
        use grid_utilities, only: trim_grid, calc_XYZ_grid
#endif
        
        character(*), parameter :: rout_name = 'kappa_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< metric equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        integer, intent(in), optional :: rich_lvl                               !< Richardson level
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          !< X, Y and Z of grid
        
        ! local variables
        real(dp), allocatable :: k_com(:,:,:,:,:)                               ! covariant and contravariant components of kappa (dim1,dim2,dim3,3,2)
        real(dp), allocatable :: k_mag(:,:,:)                                   ! magnitude of kappa (dim1,dim2,dim3)
        character(len=15) :: base_name                                          ! base name
#if ldebug
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        integer :: id                                                           ! counter
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_offset(4)                                               ! local offset of plot
        real(dp), allocatable :: XYZ_loc(:,:,:,:,:)                             ! X, Y and Z of surface in cylindrical coordinates, trimmed grid
        real(dp), allocatable :: k_com_inv(:,:,:,:)                             ! inverted cartesian components of curvature
        logical, save :: asked_for_testing = .false.                            ! whether we have been asked to test
        logical, save :: testing = .false.                                      ! whether we are testing
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (eq_style.eq.1 .and. .not.allocated(grid_eq%trigon_factors)) then
            ierr = 1
            CHCKERR('trigonometric factors not allocated')
        end if
        
        ! set up components and magnitude of kappa
        allocate(k_com(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3,2))
        allocate(k_mag(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        k_com = 0._dp
        k_mag = 0._dp
        
        ! covariant components:
        !   kappa . e_alpha = kappa_g 
        !   kappa . e_psi   = kappa_n - kappa_g h12/h22
        !   kappa . e_theta = 0
        ! and similar for Q
        k_com(:,:,:,1,1) = eq_2%kappa_g
        k_com(:,:,:,2,1) = eq_2%kappa_n - &
            &eq_2%kappa_g * &
            &eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,0)/&
            &eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        k_com(:,:,:,3,1) = 0._dp
        
        ! contravariant components:
        !   kappa . nabla alpha = kappa_n h12 + kappa_g g33/(J^2 h22)
        !   kappa . nabla psi   = kappa_n h22
        !   kappa . nabla theta = kappa_n h23 - kappa_g g13/(J^2 h22)
        ! and similar for Q
        k_com(:,:,:,1,2) = eq_2%kappa_n * &
            &eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,0) + &
            &eq_2%kappa_g * &
            &eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0) / &
            &(eq_2%jac_FD(:,:,:,0,0,0)**2*&
            &eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0))
        k_com(:,:,:,2,2) = eq_2%kappa_n * &
            &eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        k_com(:,:,:,3,2) = eq_2%kappa_n * &
            &eq_2%h_FD(:,:,:,c([3,2],.true.),0,0,0) - &
            &eq_2%kappa_g * &
            &eq_2%g_FD(:,:,:,c([3,1],.true.),0,0,0) / &
            &(eq_2%jac_FD(:,:,:,0,0,0)**2*&
            &eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0))
        
        ! transform back to unnormalized quantity
        if (use_normalization) k_com = k_com / R_0
        
        ! set plot variables
        base_name = 'kappa'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! transform coordinates
        ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,k_com,&
            &norm_disc_prec_eq,v_mag=k_mag,base_name=base_name,XYZ=XYZ)
        CHCKERR('')
        
#if ldebug
        if (ltest) then
            if (.not.asked_for_testing) then
                call writo('Do you want to plot the inversion in kappa?')
                call lvl_ud(1)
                testing = get_log(.false.)
                asked_for_testing = .true.
                call lvl_ud(-1)
            end if
            if (testing) then
                call writo('Every point in the grid is transformed by &
                    &displacing it in the direction of the curvature,')
                call writo('by a distance equal to the inverse of the &
                    &curvature. I.e. the point is displaced to the center &
                    &of curvature.')
                call lvl_ud(1)
                
                ! trim equilibrium grid
                ierr = trim_grid(grid_eq,grid_trim,norm_id)
                CHCKERR('')
                
                ! set up plot dimensions and local dimensions
                plot_dim = [grid_trim%n,3]
                plot_offset = [0,0,grid_trim%i_min-1,0]
                
                ! possibly modify if multiple equilibrium parallel jobs
                if (size(eq_jobs_lims,2).gt.1) then
                    plot_dim(1) = eq_jobs_lims(2,size(eq_jobs_lims,2)) - &
                        &eq_jobs_lims(1,1) + 1
                    plot_offset(1) = eq_jobs_lims(1,eq_job_nr) - 1
                end if
                
                ! set up inverted cartesian components
                allocate(k_com_inv(grid_trim%n(1),grid_trim%n(2),&
                    &grid_trim%loc_n_r,3))
                do id = 1,3
                    k_com_inv(:,:,:,id) = &
                        &k_com(:,:,norm_id(1):norm_id(2),id,1)/&
                        &(k_mag(:,:,norm_id(1):norm_id(2))**2)
                end do
                
                ! if VMEC, calculate trigonometric factors of trimmed grid
                if (eq_style.eq.1) then
                    ierr = calc_trigon_factors(grid_trim%theta_E,&
                        &grid_trim%zeta_E,grid_trim%trigon_factors)
                    CHCKERR('')
                end if
                
                ! calculate X, Y and Z
                allocate(XYZ_loc(grid_trim%n(1),grid_trim%n(2),&
                    &grid_trim%loc_n_r,3,3))
                ierr = calc_XYZ_grid(grid_eq,grid_trim,XYZ_loc(:,:,:,1,1),&
                    &XYZ_loc(:,:,:,1,2),XYZ_loc(:,:,:,1,3))
                CHCKERR('')
                
                ! produce a plot of center of curvature for every point
                call plot_HDF5(['cen_of_curv'],'TEST_cen_of_curv_vec',&
                    &k_com_inv,tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=XYZ_loc(:,:,:,:,1),Y=XYZ_loc(:,:,:,:,2),&
                    &Z=XYZ_loc(:,:,:,:,3),col=4,cont_plot=eq_job_nr.gt.1,&
                    &descr='center of curvature')
                
                ! displace the points to the center of curvature
                do id = 1,3
                    XYZ_loc(:,:,:,1,id) = XYZ_loc(:,:,:,1,id) + &
                        &k_com_inv(:,:,:,id)
                end do
                XYZ_loc(:,:,:,2,:) = XYZ_loc(:,:,:,1,:)
                XYZ_loc(:,:,:,3,:) = XYZ_loc(:,:,:,3,:)
                
                ! produce an inverted plot
                call plot_HDF5(['cen_of_curv_inv'],'TEST_cen_of_curv_inv_vec',&
                    &-k_com_inv,tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=XYZ_loc(:,:,:,:,1),Y=XYZ_loc(:,:,:,:,2),&
                    &Z=XYZ_loc(:,:,:,:,3),col=4,cont_plot=eq_job_nr.gt.1,&
                    &descr='center of curvature')
                
                ! clean up
                call grid_trim%dealloc()
                
                call lvl_ud(-1)
            end if
        end if
#endif
    end function kappa_plot
    
    !> Plots  \b HALF  of the  change in  the position  vectors for  2 different
    !! toroidal positions, which can correspond to a ripple.
    !!
    !! Also calculates \b HALF of the relative magnetic perturbation, which also
    !! corresponds to a ripple.
    !!
    !! Finally, if the  output grid contains a  fundamental interval \f$2\pi\f$,
    !! the proportionality between both is written to a file.
    !!
    !! \note  The  metric   factors  and  transformation  matrices  have  to  be
    !! allocated.
    !!
    !! \return ierr
    integer function delta_r_plot(grid_eq,eq_1,eq_2,XYZ,rich_lvl) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp, calc_XYZ_grid, trim_grid, &
            &calc_tor_diff, nufft
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, eq_job_nr, RZ_0, &
            &use_normalization, rank, ex_plot_style, n_procs, prop_B_tor_i, &
            &tol_zero
        use eq_vars, only: B_0
        use num_utilities, only: c, calc_int, order_per_fun
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'delta_r_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        real(dp), intent(in) :: XYZ(:,:,:,:)                                    !< X, Y and Z of grid
        integer, intent(in), optional :: rich_lvl                               !< Richardson level
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: r_lim_2D(2)                                                  ! limits of r to plot in 2-D
        integer :: max_m                                                        ! maximum mode number
        real(dp), allocatable :: XYZ_loc(:,:,:,:)                               ! X, Y and Z of surface in trimmed grid
        real(dp), allocatable :: B_com(:,:,:,:,:)                               ! covariant and contravariant components of B (dim1,dim2,dim3,3,2)
        real(dp), allocatable :: delta_B_tor(:,:,:)                             ! delta_B/B
        real(dp), allocatable :: prop_B_tor_tot(:,:)                            ! delta_r / delta_B/B in total grid
        real(dp), allocatable :: prop_B_tor(:,:,:)                              ! delta_r / delta_B/B
        real(dp), allocatable :: var_tot_loc(:)                                 ! auxilliary variable
        real(dp), allocatable :: x_2D(:,:)                                      ! x for 2-D plotting
        real(dp), allocatable :: delta_r(:,:,:)                                 ! normal displacement
        real(dp), allocatable :: theta_geo(:,:,:)                               ! geometrical poloidal angle for proportionality factor output
        real(dp), allocatable :: r_geo(:,:,:)                                   ! geometrical radius for proportionality factor output
        real(dp), allocatable :: prop_B_tor_plot(:,:)                           ! prop_B_tor and angle at last normal position for plotting
        real(dp), allocatable :: delta_r_F(:,:)                                 ! Fourier components of delta_r at last normal position
        logical :: new_file_found                                               ! name for new file found
        character(len=25) :: base_name                                          ! base name
        character(len=25) :: plot_name                                          ! plot name
        character(len=max_str_ln) :: prop_B_tor_file_name                       ! name of B_tor proportionality file
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (eq_style.eq.1 .and. .not.allocated(grid_eq%trigon_factors)) then
            ierr = 1
            CHCKERR('trigonometric factors not allocated')
        end if
        
        call writo('Calculate normal displacement')
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! calculate local X, Y and Z
        ! (Can be different from provided one for different plot grid styles)
        allocate(XYZ_loc(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,4))
        ierr = calc_XYZ_grid(grid_eq,grid_eq,XYZ_loc(:,:,:,1),XYZ_loc(:,:,:,2),&
            &XYZ_loc(:,:,:,3),R=XYZ_loc(:,:,:,4))
        CHCKERR('')
        
        ! calculate geometrical poloidal angle and radius
        allocate(theta_geo(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(r_geo(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        theta_geo = atan2(XYZ_loc(:,:,:,3)-RZ_0(2),&
            &XYZ_loc(:,:,:,4)-RZ_0(1))
        where (theta_geo.lt.0._dp) theta_geo = theta_geo + 2*pi
        r_geo = sqrt((XYZ_loc(:,:,:,3)-RZ_0(2))**2 + &
            &(XYZ_loc(:,:,:,4)-RZ_0(1))**2)
        
        ! plot
        call plot_HDF5('theta_geo','theta_geo',&
            &theta_geo(:,:,norm_id(1):norm_id(2)),&
            &tot_dim=[grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &x=XYZ(:,:,norm_id(1):norm_id(2),1),&
            &y=XYZ(:,:,norm_id(1):norm_id(2),2),&
            &z=XYZ(:,:,norm_id(1):norm_id(2),3),cont_plot=eq_job_nr.gt.1,&
            &descr='geometrical poloidal angle')
        call plot_HDF5('r_geo','r_geo',&
            &r_geo(:,:,norm_id(1):norm_id(2)),&
            &tot_dim=[grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &x=XYZ(:,:,norm_id(1):norm_id(2),1),&
            &y=XYZ(:,:,norm_id(1):norm_id(2),2),&
            &z=XYZ(:,:,norm_id(1):norm_id(2),3),cont_plot=eq_job_nr.gt.1,&
            &descr='geometrical radius')
        
        ! calculate perturbation delta_r on radius
        ierr = calc_tor_diff(r_geo,theta_geo,norm_disc_prec_eq,&
            &absolute=.true.,r=grid_eq%loc_r_F)
        CHCKERR('')
        allocate(delta_r(grid_eq%n(1),1,grid_eq%loc_n_r))
        delta_r(:,1,:) = r_geo(:,2,:)*0.5_dp                                    ! take factor half as this is absolute
        deallocate(r_geo)
        
        ! set plot variables for delta_r
        base_name = 'delta_r'
        plot_name = 'delta_r'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! plot
        call plot_HDF5(trim(plot_name),trim(base_name),&
            &delta_r(:,:,norm_id(1):norm_id(2)),&
            &tot_dim=[grid_trim%n(1),1,grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &x=XYZ(:,2:2,norm_id(1):norm_id(2),1),&
            &y=XYZ(:,2:2,norm_id(1):norm_id(2),2),&
            &z=XYZ(:,2:2,norm_id(1):norm_id(2),3),&
            &cont_plot=eq_job_nr.gt.1,&
            &descr='plasma position displacement')
        
        ! calculate NUFFT for last point
        if (rank.eq.n_procs-1) then
            ! plot delta_r at last normal position
            call print_ex_2D(trim(plot_name),trim(base_name),&
                &delta_r(1:grid_eq%n(1)-1,1,grid_eq%loc_n_r),&
                &x=theta_geo(1:grid_eq%n(1)-1,2,grid_eq%loc_n_r),&
                &draw=.false.)
            call draw_ex([trim(plot_name)],trim(base_name),1,1,.false.)
            
            ! calculate and plot Fourier modes
            ierr = nufft(theta_geo(1:grid_eq%n(1)-1,2,grid_eq%loc_n_r),&
                &delta_r(1:grid_eq%n(1)-1,1,grid_eq%loc_n_r),delta_r_F)
            CHCKERR('')
            base_name = trim(base_name)//'_F'
            call print_ex_2D(['delta_r cos','delta_r sin'],trim(base_name),&
                &delta_r_F,draw=.false.)
            call draw_ex(['delta_r cos','delta_r sin'],trim(base_name),2,1,&
                &.false.)
            
            ! mode outputs, (p)recycle "prop_B_tor_file_*"
            new_file_found = .false.
            prop_B_tor_file_name = base_name
            kd = 1
            do while (.not.new_file_found)
                open(prop_B_tor_i,FILE=trim(prop_B_tor_file_name)//'_'//&
                    &trim(i2str(kd))//'.dat',IOSTAT=ierr,STATUS='old')
                if (ierr.eq.0) then
                    kd = kd + 1
                else
                    prop_B_tor_file_name = trim(prop_B_tor_file_name)//&
                        &'_'//trim(i2str(kd))//'.dat'
                    new_file_found = .true.
                    open(prop_B_tor_i,FILE=trim(prop_B_tor_file_name),&
                        &IOSTAT=ierr,STATUS='new')
                    CHCKERR('Failed to open file')
                end if
            end do
            
            ! write to output file
            max_m = 20
            write(prop_B_tor_i,'("# ",2(A5," "),2(A23," "))') &
                &'N', 'M', 'delta_c', 'delta_s'
            do id = 1,min(size(delta_r_F,1),max_m+1)
                write(prop_B_tor_i,'("  ",2(I5," "),2(ES23.16," "))') &
                    &18, id-1, delta_r_F(id,:)
            end do
            
            ! close
            close(prop_B_tor_i)
        end if
        
        call lvl_ud(-1)
        
        call writo('Calculate magnetic field ripple')
        call lvl_ud(1)
        
        ! Calculate cylindrical components of B
        allocate(B_com(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3,2))
        B_com = 0._dp
        do id = 1,3
            B_com(:,:,:,id,1) = eq_2%g_FD(:,:,:,c([3,id],.true.),0,0,0)/&
                &eq_2%jac_FD(:,:,:,0,0,0)
        end do
        B_com(:,:,:,3,2) = 1._dp/eq_2%jac_FD(:,:,:,0,0,0)
        if (use_normalization) B_com = B_com * B_0
        ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,B_com,norm_disc_prec_eq,&
            &max_transf=4)
        CHCKERR('')
        
        ! calculate perturbation
        ierr = calc_tor_diff(B_com,theta_geo,norm_disc_prec_eq,&
            &r=grid_eq%loc_r_F)
        CHCKERR('')
        allocate(delta_B_tor(grid_eq%n(1),1,grid_eq%loc_n_r))
        delta_B_tor(:,1,:) = 0.5_dp*sum(B_com(:,2,:,2,:),3)
        deallocate(B_com)
        
        ! set plot variables for delta_B_tor
        base_name = 'delta_B_tor'
        plot_name = 'delta_B_tor'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! plot with HDF5
        call plot_HDF5(trim(plot_name),trim(base_name),&
            &delta_B_tor(:,:,norm_id(1):norm_id(2)),&
            &tot_dim=[grid_trim%n(1),1,grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &x=XYZ(:,2:2,norm_id(1):norm_id(2),1),&
            &y=XYZ(:,2:2,norm_id(1):norm_id(2),2),&
            &z=XYZ(:,2:2,norm_id(1):norm_id(2),3),cont_plot=eq_job_nr.gt.1,&
            &descr='plasma position displacement')
        
        call lvl_ud(-1)
        
        call writo('Calculate proportionality factor')
        call lvl_ud(1)
        
        allocate(prop_B_tor(grid_eq%n(1),1,grid_eq%loc_n_r))
        prop_B_tor = 0.0_dp
        where (abs(delta_B_tor).gt.tol_zero) prop_B_tor = delta_r / delta_B_tor
        
        ! set plot variables for prop_B_tor
        base_name = 'prop_B_tor'
        plot_name = trim(base_name)
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! plot with HDF5
        call plot_HDF5(plot_name,trim(base_name),&
            &prop_B_tor(:,:,norm_id(1):norm_id(2)),&
            &tot_dim=[grid_trim%n(1),1,grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &x=XYZ(:,2:2,norm_id(1):norm_id(2),1),&
            &y=XYZ(:,2:2,norm_id(1):norm_id(2),2),&
            &z=XYZ(:,2:2,norm_id(1):norm_id(2),3),&
            &cont_plot=eq_job_nr.gt.1,&
            &descr='delta_r divided by delta_B_tor')
        
        ! plot in 2-D
        r_lim_2D = [1,grid_trim%n(3)]
        if (ex_plot_style.eq.2) r_lim_2D(1) = &
            &max(r_lim_2D(1),r_lim_2D(2)-127+1)                                 ! Bokeh can only handle 255/2 input arguments
        if (rank.eq.0) then
            allocate(x_2D(grid_trim%n(1),grid_trim%n(3)))
        end if
        allocate(prop_B_tor_tot(grid_trim%n(1),grid_trim%n(3)))
        do id = 1,grid_trim%n(1)
            ! prop_B_tor (all procs)
            ierr = get_ser_var(prop_B_tor(id,1,norm_id(1):norm_id(2)),&
                &var_tot_loc,scatter=.true.)
            CHCKERR('')
            prop_B_tor_tot(id,:) = max(min(var_tot_loc,2._dp),-2._dp)           ! limit to avoid division by small delta_B_tor values
            deallocate(var_tot_loc)
            
            ! theta_geo (only master)
            ierr = get_ser_var(theta_geo(id,2,norm_id(1):norm_id(2)),&
                &var_tot_loc)
            CHCKERR('')
            if (rank.eq.0) x_2D(id,:) = var_tot_loc/pi
            deallocate(var_tot_loc)
        end do
        if (rank.eq.0) then
            call print_ex_2D([trim(plot_name)],trim(base_name),&
                &prop_B_tor_tot(:,r_lim_2D(1):r_lim_2D(2)),&
                &x=x_2D(:,r_lim_2D(1):r_lim_2D(2)),draw=.false.)
            call draw_ex([trim(plot_name)],trim(base_name),&
                &r_lim_2D(2)-r_lim_2D(1)+1,1,.false.)
        end if
        
        call lvl_ud(-1)
        
        call writo('Output to file as function of poloidal flux angle')
        call lvl_ud(1)
        
        ! only last rank outputs to file (it has the last normal position)
        if (rank.eq.n_procs-1) then
            ! find new file name
            new_file_found = .false.
            prop_B_tor_file_name = 'prop_B_tor'
            kd = 1
            do while (.not.new_file_found)
                open(prop_B_tor_i,FILE=trim(prop_B_tor_file_name)//'_'//&
                    &trim(i2str(kd))//'.dat',IOSTAT=ierr,STATUS='old')
                if (ierr.eq.0) then
                    kd = kd + 1
                else
                    prop_B_tor_file_name = trim(prop_B_tor_file_name)//&
                        &'_'//trim(i2str(kd))//'.dat'
                    new_file_found = .true.
                    open(prop_B_tor_i,FILE=trim(prop_B_tor_file_name),&
                        &IOSTAT=ierr,STATUS='new')
                    CHCKERR('Failed to open file')
                end if
            end do
            
            ! user output
            call writo('Save toroidal field proportionality factor in &
                &file "'//trim(prop_B_tor_file_name)//'"',persistent=.true.)
            
            ! order  output, taking away last  point as it is  equivalent to
            ! the first
            ierr = order_per_fun(reshape([&
                &theta_geo(1:grid_trim%n(1)-1,1,grid_trim%loc_n_r),&
                &prop_B_tor_tot(1:grid_trim%n(1)-1,grid_trim%n(3))],&
                &[grid_trim%n(1)-1,2]),prop_B_tor_plot,0)
            CHCKERR('')
            
            ! write to output file
            write(prop_B_tor_i,'("# ",2(A23," "))') &
                &'pol. angle [ ]', 'prop. factor [ ]'
            do id = 1,grid_trim%n(1)-1
                write(prop_B_tor_i,'("  ",2(ES23.16," "))') &
                    &prop_B_tor_plot(id,:)
            end do
            
            ! close
            close(prop_B_tor_i)
        end if
        
        call lvl_ud(-1)
        
        ! clean up
        call grid_trim%dealloc()
    end function delta_r_plot

    !> Divides the equilibrium jobs.
    !!
    !! For PB3D,  the entire  parallel range  has to be  calculated, but  due to
    !! memory limits this  has to be split  up in pieces. Every piece  has to be
    !! able to  contain the equilibrium variables  (see note below), as  well as
    !! the  vectorial  perturbation variables.  These  are  later combined  into
    !! tensorial variables and integrated.
    !! 
    !! The  equilibrium variables  have to  be  operated on  to calculate  them,
    !! which  translates to  a scale  factor \c  mem_scale_fac. However,  in the
    !! perturbation phase,  when they are  just used,  this scale factor  is not
    !! needed.
    !! 
    !! In its most  extreme form, the division in equilibrium  jobs would be the
    !! individual  calculation  on a  fundamental  integration  integral of  the
    !! parallel points:
    !!  - for \c magn_int_style=1 (trapezoidal), this is 1 point,
    !!  - for \c magn_int_style=2 (Simpson 3/8), this is 3 points.
    !!
    !! For  HELENA,  the parallel  derivatives  are  calculated discretely,  the
    !! equilibrium and  vectorial perturbation variables are  tabulated first in
    !! this HELENA  grid. This  happens in  the first  Richardson level.  In all
    !! Richardson levels,  afterwards, these  variables are interpolated  in the
    !! angular directions. In this case, therefore,  there can be no division of
    !! this HELENA output interval for the first Richardson level.
    !!
    !! This procedure  does the  job of  dividing the  grids setting  the global
    !! variables \c eq_jobs_lims.
    !!
    !! The integration of the tensorial perturbation variables is adjusted:
    !!  - If  the first job  of the parallel jobs  and not the  first Richardson
    !!  level: add half  of the integrated tensorial  perturbation quantities of
    !!  the previous level.
    !!  -  If  not the  first  job  of the  parallel  jobs,  add the  integrated
    !!  tensorial perturbation quantities to those of the previous parallel job,
    !!  same Richardson level.
    !!
    !! In fact,  the equilibrium jobs  have much  in common with  the Richardson
    !! levels,  as is  attested by  the existence  of the  routines do_eq()  and
    !! eq_info(), which are equivalent to do_rich() and rich_info().
    !!
    !! In  POST, finally,  the situation  is slightly  different for  HELENA, as
    !! all  the requested  variables  have to  fit,  including the  interpolated
    !! variables, as they are stored whereas  in PB3D they are not. The parallel
    !! range to be  taken is then the  one of the output grid,  including a base
    !! range for the variables tabulated on  the HELENA grid. Also, for extended
    !! output grids,  the size  of the  grid in  the secondary  angle has  to be
    !! included in \c n_par_X (i.e. toroidal when poloidal flux is used and vice
    !! versa). Furthermore, multiple equilibrium jobs are allowed.
    !!
    !! To this  end, optionally, a base  number can be provided  for \c n_par_X,
    !! that is always added to the number of points in the divided \c n_par_X.
    !!
    !! \note For  PB3D, only the  variables \c g_FD, \c  h_FD and \c  jac_FD are
    !! counted, as the equilibrium variables and the transformation matrices are
    !! deleted after use. Also, \c S, \c sigma, \c kappa_n and \c kappa_g can be
    !! neglected  as they  do not  contain  derivatives and  are therefore  much
    !! smaller. in both routines  calc_memory_eq() and calc_memory_x(), however,
    !! a 50% safety factor is used to account for this somewhat.
    !!
    !! \return ierr
    integer function divide_eq_jobs(n_par_X,arr_size,n_div,n_div_max,&
        &n_par_X_base,range_name) result(ierr)
        
        use num_vars, only: max_tot_mem, max_X_mem, magn_int_style, &
            &mem_scale_fac, prog_style
        use rich_vars, only: rich_lvl
        use eq_utilities, only: calc_memory_eq
        use X_utilities, only: calc_memory_X
        use X_vars, only: n_mod_X
        
        character(*), parameter :: rout_name = 'divide_eq_jobs'
        
        ! input / output
        integer, intent(in) :: n_par_X                                          !< number of parallel points to be divided
        integer, intent(in) :: arr_size(2)                                      !< array size (using loc_n_r) for eq_2 and X_1 variables
        integer, intent(inout) :: n_div                                         !< final number of divisions
        integer, intent(in), optional :: n_div_max                              !< maximum n_div
        integer, intent(in), optional :: n_par_X_base                           !< base n_par_X, undivisible
        character(len=*), intent(in), optional :: range_name                    !< name of range
        
        ! local variables
        real(dp) :: mem_size(2)                                                 ! approximation of memory required for eq_2 and X_1 variables
        integer :: n_par_X_base_loc                                             ! local n_par_X_base
        integer :: max_mem_req                                                  ! maximum memory required
        integer :: n_div_max_loc                                                ! maximum n_div
        integer :: n_par_range                                                  ! nr. of points in range
        character(len=max_str_ln) :: range_message                              ! message about how many ranges
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: range_name_loc                             ! local range name
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Dividing the equilibrium jobs')
        call lvl_ud(1)
        
        ! set up width of fundamental interval
        select case (magn_int_style)
            case (1)                                                            ! Trapezoidal rule
                fund_n_par = 1   
            case (2)                                                            ! Simpson's 3/8 rule
                fund_n_par = 3   
        end select
        
        ! set up local n_par_X_base
        n_par_X_base_loc = 0
        if (present(n_par_X_base)) n_par_X_base_loc = n_par_X_base
        
        ! set up local range name
        range_name_loc = 'parallel points'
        if (present(range_name)) range_name_loc = trim(range_name)
        
        ! setup auxilliary variables
        n_div = 0
        mem_size = huge(1._dp)
        if (rich_lvl.eq.1) then
            n_div_max_loc = (n_par_X-1)/fund_n_par
        else
            n_div_max_loc = n_par_X/fund_n_par
        end if
        if (present(n_div_max)) n_div_max_loc = n_div_max
        n_div_max_loc = max(1,n_div_max_loc)                                    ! cannot have less than 1 piece
        
        ! calculate largest possible range of parallel points fitting in memory
        do while (max_tot_mem.lt.sum(mem_size))                                 ! combined mem_size has to be smaller than total memory
            n_div = n_div + 1
            n_par_range = ceiling(n_par_X*1._dp/n_div + n_par_X_base_loc)
            ierr = calc_memory_eq(arr_size(1),n_par_range,mem_size(1))
            CHCKERR('')
            mem_size(1) = mem_size(1)*mem_scale_fac                             ! operations have to be done on eq, it is not just stored.
            ierr = calc_memory_X(1,arr_size(2)*n_par_range,n_mod_X,mem_size(2)) ! all modes have to be stored
            CHCKERR('')
            if (n_div.gt.n_div_max_loc) then                                    ! still not enough memory
                ierr = 1
                err_msg = 'The memory limit is too low, need more than '//&
                    &trim(i2str(max_mem_req))//'MB'
                CHCKERR(err_msg)
            end if
            max_mem_req = ceiling(sum(mem_size))
        end do
        if (n_div.gt.1) then
            range_message = 'The '//trim(i2str(n_par_X))//' '//&
                &trim(range_name_loc)//' are split into '//&
                &trim(i2str(n_div))//' and '//trim(i2str(n_div))//&
                &' collective jobs are done serially'
        else
            range_message = 'The '//trim(i2str(n_par_X))//' '//&
                &trim(range_name_loc)//' can be done without splitting them'
        end if
        call writo(range_message)
        call writo('The total memory for all processes together is estimated &
            &to be about '//trim(i2str(ceiling(sum(mem_size))))//'MB')
        call lvl_ud(1)
        call writo('(maximum: '//trim(i2str(ceiling(max_tot_mem)))//'&
            &MB, user specified)')
        call lvl_ud(-1)
        
        if (prog_style.eq.1) then
            ! calculate max memory available for perturbation calculations
            mem_size(1) = mem_size(1)/mem_scale_fac                             ! rescale by mem_scale_fac as eq is now just stored
            max_X_mem = max_tot_mem - sum(mem_size)
            call writo('In the perturbation phase, the equilibrium variables &
                &are not being operated on:')
            call lvl_ud(1)
            call writo('This translates to a scale factor 1/'//&
                &trim(r2strt(mem_scale_fac)))
            call writo('Therefore, the memory left for the perturbation phase &
                &is '//trim(i2str(ceiling(max_X_mem)))//'MB')
            call lvl_ud(-1)
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Equilibrium jobs divided')
    end function divide_eq_jobs
    
    !> Calculate \c eq_jobs_lims.
    !!
    !! Take  into  account  that  every  job   has  to  start  and  end  at  the
    !! start  and end  of a  fundamental integration  interval, as  discussed in
    !! divide_eq_jobs():
    !!  - for \c magn_int_style=1 (trapezoidal), this is 1 point,
    !!  - for \c magn_int_style=2 (Simpson 3/8), this is 3 points.
    !!
    !! for POST, there are no Richardson levels,  and there has to be overlap of
    !! one always, in order to have correct composite integrals of the different
    !! regions.
    !!
    !! \note The  \c n_par_X passed into  this procedure refers to  the quantity
    !! that is already possibly halved if the Richardson level is higher than 1.
    !! This  information is  then  reflected in  the  eq_jobs_lims, which  refer
    !! to  the local  limits,  i.e.  only the  parallel  points currently  under
    !! consideration.
    !!
    !! \return ierr
    integer function calc_eq_jobs_lims(n_par_X,n_div) result(ierr)
        use num_vars, only: prog_style, eq_jobs_lims, eq_job_nr
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'calc_eq_jobs_lims'
        
        ! input / output
        integer, intent(in) :: n_par_X                                          !< number of parallel points in this Richardson level
        integer, intent(in) :: n_div                                            !< nr. of divisions of parallel ranges
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: ol_width                                                     ! overlap width
        integer, allocatable :: n_par(:)                                        ! number of points per range
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! Also initialize eq_job_nr to 0 as it is incremented by "do_eq".
        eq_job_nr = 0
        
        allocate(n_par(n_div))
        n_par = n_par_X/n_div                                                   ! number of parallel points on this processor
        n_par(1:mod(n_par_X,n_div)) = n_par(1:mod(n_par_X,n_div)) + 1           ! add a point to if there is a remainder
        CHCKERR('')
        
        ! (re)allocate equilibrium jobs limits
        if (allocated(eq_jobs_lims)) deallocate(eq_jobs_lims)
        allocate(eq_jobs_lims(2,n_div))
        
        ! set overlap width
        select case (prog_style)
            case (1)                                                            ! PB3D
                if (rich_lvl.eq.1) then
                    ol_width = 1
                else
                    ol_width = 0
                end if
            case (2)                                                            ! POST
                ol_width = 1
        end select
        
        ! loop over divisions
        do id = 1,n_div
            ! setup  first  guess,  without   taking  into  account  fundamental
            ! integration intervals
            if (id.eq.1) then
                eq_jobs_lims(1,id) = 1
            else
                eq_jobs_lims(1,id) = eq_jobs_lims(2,id-1) + 1 - ol_width
            end if
            eq_jobs_lims(2,id) = sum(n_par(1:id))
            
            ! take into account fundamental interval
            if (n_div.gt.1) then
                eq_jobs_lims(2,id) = eq_jobs_lims(1,id) - 1 + ol_width + &
                    &fund_n_par * max(1,nint(&
                    &(eq_jobs_lims(2,id)-eq_jobs_lims(1,id)+1._dp-ol_width)/&
                    &fund_n_par))
            end if
        end do
        
        ! test whether end coincides with sum of n_par
        if (eq_jobs_lims(2,n_div).ne.sum(n_par)) then
            ierr = 1
            err_msg = 'Limits don''t match, try with more memory or lower &
                &magn_int_style'
            CHCKERR(err_msg)
        end if
    end function calc_eq_jobs_lims
    
#if ldebug
    !> See if \c T_EF it complies with the theory of \cite Weyens3D.
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_T_EF(grid_eq,eq_1,eq_2) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style
        use grid_utilities, only: trim_grid
        use num_utilities, only: c
        use output_ops, only: plot_diff_HDF5
        
        character(*), parameter :: rout_name = 'test_T_EF'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        
        ! local variables
        integer :: id, kd                                                       ! counter
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        real(dp), allocatable :: res(:,:,:,:)                                   ! calculated result
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether T_EF complies with the theory')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,9))
        res = 0.0_dp
        
        ! calc  res,  depending  on flux  used  in F  coordinates  and on  which
        ! equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                if (use_pol_flux_F) then                                        ! using poloidal flux
                    ! calculate T_EF(1,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,1) = &
                            &grid_eq%theta_F(:,:,kd)*eq_1%q_saf_E(kd,1) + &
                            &eq_2%L_E(:,:,kd,1,0,0)*eq_1%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(2,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,2) = &
                            &(1._dp + eq_2%L_E(:,:,kd,0,1,0))*eq_1%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(3,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,3) = -1._dp + &
                            &eq_2%L_E(:,:,kd,0,0,1)*eq_1%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(1,2)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,4) = eq_1%flux_p_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(1,3)
                    res(:,:,:,7) = eq_2%L_E(:,:,norm_id(1):norm_id(2),1,0,0)
                    ! calculate T_EF(2,3)
                    res(:,:,:,8) = 1._dp + &
                        &eq_2%L_E(:,:,norm_id(1):norm_id(2),0,1,0)
                    ! calculate T_EF(3,3)
                    res(:,:,:,9) = eq_2%L_E(:,:,norm_id(1):norm_id(2),0,0,1)
                else                                                            ! using toroidal flux
                    ! calculate T_EF(1,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,1) = - eq_2%L_E(:,:,kd,1,0,0) &
                            &+ grid_eq%zeta_E(:,:,kd)*eq_1%rot_t_E(kd,1)
                    end do
                    ! calculate T_EF(2,1)
                    res(:,:,:,2) = &
                        &- (1._dp + eq_2%L_E(:,:,norm_id(1):norm_id(2),0,1,0))
                    ! calculate T_EF(3,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,3) = &
                            &eq_1%rot_t_E(kd,0) - eq_2%L_E(:,:,kd,0,0,1)
                    end do
                    ! calculate T_EF(1,2)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,4) = &
                            &- eq_1%flux_t_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(3,3)
                    res(:,:,:,9) = -1._dp
                end if
            case (2)                                                            ! HELENA
                if (use_pol_flux_F) then                                        ! using poloidal flux
                    ! calculate T_EF(1,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,1) = &
                            &- eq_1%q_saf_E(kd,1)*grid_eq%theta_E(:,:,kd)
                    end do
                    ! calculate T_EF(2,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,2) = - eq_1%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(3,1)
                    res(:,:,:,3) = 1._dp
                    ! calculate T_EF(1,2)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,4) = eq_1%flux_p_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(2,3)
                    res(:,:,:,8) = 1._dp
                else                                                            ! using toroidal flux
                    ! calculate T_EF(1,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,1) = &
                            &eq_1%rot_t_E(kd,1)*grid_eq%zeta_E(:,:,kd)
                    end do
                    ! calculate T_EF(2,1)
                    res(:,:,:,2) = -1._dp
                    ! calculate T_EF(3,1)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,3) = eq_1%rot_t_E(kd,0)
                    end do
                    ! calculate T_EF(1,2)
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,4) = eq_1%flux_t_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(3,3)
                    res(:,:,:,9) = 1._dp
                end if
        end select
        
        ! set up plot variables for calculated values
        do id = 1,3
            do kd = 1,3
                ! user output
                call writo('Testing T_EF('//trim(i2str(kd))//','//&
                    &trim(i2str(id))//')')
                call lvl_ud(1)
                
                ! set some variables
                file_name = 'TEST_T_EF_'//trim(i2str(kd))//'_'//trim(i2str(id))
                description = 'Testing calculated with given value for T_EF('//&
                    &trim(i2str(kd))//','//trim(i2str(id))//')'
                
                ! plot difference
                call plot_diff_HDF5(res(:,:,:,c([kd,id],.false.)),&
                    &eq_2%T_EF(:,:,norm_id(1):norm_id(2),c([kd,id],.false.),&
                    &0,0,0),file_name,tot_dim,loc_offset,description,&
                    &output_message=.true.)
                
                call lvl_ud(-1)
            end do
        end do
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_T_EF
    
    !> Tests whether \f$ \frac{\partial^2}{\partial u_i \partial u_j} h_\text{H}
    !! \f$ is calculated correctly.
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_D12h_H(grid_eq,eq) result(ierr)
        use grid_utilities, only: trim_grid
        use num_utilities, only: c, spline
        use num_vars, only: norm_disc_prec_eq
        use HELENA_vars, only: ias
        
        character(*), parameter :: rout_name = 'test_D12h_H'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       !< metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: id, jd, kd, ld                                               ! counters
        integer :: bcs(2,2)                                                     ! boundary conditions (theta(even), theta(odd))
        integer :: bc_ld(6)                                                     ! boundary condition type for each metric element
        real(dp) :: bcs_val(2,3)                                                ! values for boundary conditions
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether D1 D2 h_H is calculated correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,6))
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! set up boundary conditions
        if (ias.eq.0) then                                                      ! top-bottom symmetric
            bcs(:,1) = [1,1]                                                    ! theta(even): zero first derivative
            bcs(:,2) = [2,2]                                                    ! theta(odd): zero first derivative
        else
            bcs(:,1) = [-1,-1]                                                  ! theta(even): periodic
            bcs(:,2) = [-1,-1]                                                  ! theta(odd): periodic
        end if
        bcs_val = 0._dp
        
        ! boundary condition type for each metric element
        bc_ld = [1,2,0,1,0,1]                                                   ! 1: even, 2: odd, 0: zero
        
        ! calculate D1 D2 h_H alternatively
        do ld = 1,6
            do kd = norm_id(1),norm_id(2)
                do jd = 1,grid_trim%n(2)
                    if (bc_ld(ld).eq.0) cycle                                   ! quantity is zero
                    
                    ierr = spline(grid_eq%theta_E(:,jd,kd),&
                        &eq%h_E(:,jd,kd,ld,0,1,0),grid_eq%theta_E(:,jd,kd),&
                        &res(:,jd,kd-norm_id(1)+1,ld),ord=norm_disc_prec_eq,&
                        &deriv=1,bcs=bcs(:,bc_ld(ld)),&
                        &bcs_val=bcs_val(:,bc_ld(ld)))
                    CHCKERR('')
                end do
            end do
        end do
        
        ! set up plot variables for calculated values
        do id = 1,3
            do kd = 1,3
                ! user output
                call writo('Testing h_H('//trim(i2str(kd))//','//&
                    &trim(i2str(id))//')')
                call lvl_ud(1)
                
                ! set some variables
                file_name = 'TEST_D12h_H_'//trim(i2str(kd))//'_'//&
                    &trim(i2str(id))
                description = 'Testing calculated with given value for D12h_H('&
                    &//trim(i2str(kd))//','//trim(i2str(id))//')'
                
                ! plot difference
                call plot_diff_HDF5(res(:,:,:,c([kd,id],.true.)),&
                    &eq%h_E(:,:,norm_id(1):norm_id(2),&
                    &c([kd,id],.true.),1,1,0),file_name,tot_dim,loc_offset,&
                    &description,output_message=.true.)
                
                call lvl_ud(-1)
            end do
        end do
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_D12h_H
    
    !> Performs  tests  on \f$  \mathcal{J}_\text{F}\f$.
    !!
    !!  - comparing it with the determinant of \f$g_\text{F}\f$
    !!  - comparing it with the direct formula
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_jac_F(grid_eq,eq_1,eq_2) result(ierr)
        use num_vars, only: eq_style, use_pol_flux_F
        use grid_utilities, only: trim_grid
        use num_utilities, only: calc_det
        use HELENA_vars, only: h_H_33, RBphi_H
        
        character(*), parameter :: rout_name = 'test_jac_F'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in), target :: eq_1                             !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        real(dp), allocatable :: res(:,:,:)                                     ! result variable
        integer :: jd, kd                                                       ! counters
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        real(dp), pointer :: Dflux(:) => null()                                 ! points to D flux_p or D flux_t in E coordinates
        integer :: pmone                                                        ! plus or minus one
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test the calculation of the Jacobian in Flux &
            &coordinates')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! 1. Compare with determinant of g_F
        
        ! calculate Jacobian from determinant of g_F
        ierr = calc_det(res,eq_2%g_F(:,:,norm_id(1):norm_id(2),:,0,0,0),3)
        CHCKERR('')
        
        ! set some variables
        file_name = 'TEST_jac_F_1'
        description = 'Testing whether the Jacobian in Flux coordinates is &
            &consistent with determinant of metric matrix'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &eq_2%jac_F(:,:,norm_id(1):norm_id(2),0,0,0)**2,&
            &file_name,tot_dim,loc_offset,description,output_message=.true.)
        
        ! 2. Compare with explicit formula
        
        ! set up Dflux and pmone
        if (use_pol_flux_F) then                                                ! using poloidal flux
            Dflux => eq_1%flux_p_E(:,1)
            pmone = 1                                                           ! flux_p_V = flux_p_F
        else
            Dflux => eq_1%flux_t_E(:,1)
            pmone = -1                                                          ! flux_t_V = - flux_t_F
        end if
        
        ! calculate jac_F
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                do kd = norm_id(1),norm_id(2)
                    res(:,:,kd-norm_id(1)+1) = pmone * 2*pi * &
                        &eq_2%R_E(:,:,kd,0,0,0)*&
                        &(eq_2%R_E(:,:,kd,1,0,0)*eq_2%Z_E(:,:,kd,0,1,0) - &
                        &eq_2%Z_E(:,:,kd,1,0,0)*eq_2%R_E(:,:,kd,0,1,0)) / &
                        &(Dflux(kd)*(1+eq_2%L_E(:,:,kd,0,1,0)))
                end do
            case (2)                                                            ! HELENA
                do kd = norm_id(1),norm_id(2)
                    do jd = 1,grid_trim%n(2)
                        res(:,jd,kd-norm_id(1)+1) = eq_1%q_saf_E(kd,0)/&
                            &(h_H_33(:,kd+grid_eq%i_min-1)*&
                            &RBphi_H(kd+grid_eq%i_min-1,0))                     ! h_H_33 = 1/R^2 and RBphi_H = F are tabulated in eq. grid
                    end do
                end do
        end select
        
        ! set some variables
        file_name = 'TEST_jac_F_2'
        description = 'Testing whether the Jacobian in Flux coordinates is &
            &consistent with explicit formula'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &eq_2%jac_F(:,:,norm_id(1):norm_id(2),0,0,0),file_name,tot_dim,&
            &loc_offset,description,output_message=.true.)
        
        ! clean up
        nullify(Dflux)
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_jac_F
    
    !> Tests whether \f$g_\text{V}\f$ is calculated correctly.
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_g_V(grid_eq,eq) result(ierr)
        use grid_utilities, only: trim_grid
        use num_utilities, only: c
        
        character(*), parameter :: rout_name = 'test_g_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       !< metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether g_V is calculated correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,6))
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! calculate g_V(1,1)
        res(:,:,:,1) = eq%R_E(:,:,norm_id(1):norm_id(2),1,0,0)**2 + &
            &eq%Z_E(:,:,norm_id(1):norm_id(2),1,0,0)**2
        ! calculate g_V(2,1)
        res(:,:,:,2) = eq%R_E(:,:,norm_id(1):norm_id(2),1,0,0)*&
            &eq%R_E(:,:,norm_id(1):norm_id(2),0,1,0) + &
            &eq%Z_E(:,:,norm_id(1):norm_id(2),1,0,0)*&
            &eq%Z_E(:,:,norm_id(1):norm_id(2),0,1,0)
        ! calculate g_V(3,1)
        res(:,:,:,3) = eq%R_E(:,:,norm_id(1):norm_id(2),1,0,0)*&
            &eq%R_E(:,:,norm_id(1):norm_id(2),0,0,1) + &
            &eq%Z_E(:,:,norm_id(1):norm_id(2),1,0,0)*&
            &eq%Z_E(:,:,norm_id(1):norm_id(2),0,0,1)
        ! calculate g_V(2,2)
        res(:,:,:,4) = eq%R_E(:,:,norm_id(1):norm_id(2),0,1,0)**2 + &
            &eq%Z_E(:,:,norm_id(1):norm_id(2),0,1,0)**2
        ! calculate g_V(3,2)
        res(:,:,:,5) = eq%R_E(:,:,norm_id(1):norm_id(2),0,1,0)*&
            &eq%R_E(:,:,norm_id(1):norm_id(2),0,0,1) + &
            &eq%Z_E(:,:,norm_id(1):norm_id(2),0,1,0)*&
            &eq%Z_E(:,:,norm_id(1):norm_id(2),0,0,1)
        ! calculate g_V(3,3)
        res(:,:,:,6) = eq%R_E(:,:,norm_id(1):norm_id(2),0,0,1)**2 + &
            &eq%Z_E(:,:,norm_id(1):norm_id(2),0,0,1)**2 + &
            &eq%R_E(:,:,norm_id(1):norm_id(2),0,0,0)**2
        
        ! set up plot variables for calculated values
        do id = 1,3
            do kd = 1,3
                ! user output
                call writo('Testing g_V('//trim(i2str(kd))//','//&
                    &trim(i2str(id))//')')
                call lvl_ud(1)
                
                ! set some variables
                file_name = 'TEST_g_V_'//trim(i2str(kd))//'_'//trim(i2str(id))
                description = 'Testing calculated with given value for g_V('//&
                    &trim(i2str(kd))//','//trim(i2str(id))//')'
                
                ! plot difference
                call plot_diff_HDF5(res(:,:,:,c([kd,id],.true.)),&
                    &eq%g_E(:,:,norm_id(1):norm_id(2),c([kd,id],.true.),&
                    &0,0,0),file_name,tot_dim,loc_offset,description,&
                    &output_message=.true.)
                
                call lvl_ud(-1)
            end do
        end do
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_g_V
    
    !> Tests whether \f$\mathcal{J}_\text{V}\f$ is calculated correctly.
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_jac_V(grid_eq,eq) result(ierr)
        use grid_utilities, only: trim_grid
        use num_utilities, only: calc_det
        
        character(*), parameter :: rout_name = 'test_jac_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       !< metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        real(dp), allocatable :: res(:,:,:)                                     ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test the calculation of the Jacobian in the VMEC &
            &coords.')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! calculate Jacobian from determinant of g_V
        ierr = calc_det(res,eq%g_E(:,:,norm_id(1):norm_id(2),:,0,0,0),3)
        CHCKERR('')
        
        ! set some variables
        file_name = 'TEST_jac_V'
        description = 'Testing whether the Jacobian in VMEC coordinates is &
            &consistent with determinant of metric matrix'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &eq%jac_E(:,:,norm_id(1):norm_id(2),0,0,0)**2,&
            &file_name,tot_dim,loc_offset,description,output_message=.true.)
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_jac_V
    
    !> Tests whether \f$\vec{B}_\text{F}\f$ is calculated correctly.
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_B_F(grid_eq,eq_1,eq_2) result(ierr)
        use num_vars, only: eq_style
        use grid_utilities, only: trim_grid
        use num_utilities, only: c
        use VMEC_utilities, only: fourier2real
        use VMEC_vars, only: B_V_sub_s, B_V_sub_c, B_V_c, B_V_s, is_asym_V
        
        character(*), parameter :: rout_name = 'test_B_F'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: norm_id_f(2)                                                 ! norm_id transposed to full grid
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        real(dp), allocatable :: res2(:,:,:,:)                                  ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether the magnetic field is calculated &
            &correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set norm_id_f for quantities tabulated on full grid
        norm_id_f = grid_eq%i_min+norm_id-1
        
        ! set up res and res2
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,4))
        allocate(res2(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,4))
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! get covariant components and magnitude of B_E
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                do id = 1,3
                    ierr = fourier2real(&
                        &B_V_sub_c(:,norm_id_f(1):norm_id_f(2),id),&
                        &B_V_sub_s(:,norm_id_f(1):norm_id_f(2),id),&
                        &grid_eq%trigon_factors(:,:,:,norm_id(1):&
                        &norm_id(2),:),res(:,:,:,id))
                    CHCKERR('')
                end do
                ierr = fourier2real(B_V_c(:,norm_id_f(1):norm_id_f(2)),&
                    &B_V_s(:,norm_id_f(1):norm_id_f(2)),&
                    &grid_eq%trigon_factors(:,:,:,norm_id(1):norm_id(2),:),&
                    &res(:,:,:,4),[.true.,is_asym_V])
                CHCKERR('')
            case (2)                                                            ! HELENA
                res(:,:,:,1) = &
                    &eq_2%g_E(:,:,norm_id(1):norm_id(2),c([1,2],.true.),0,0,0)
                res(:,:,:,2) = &
                    &eq_2%g_E(:,:,norm_id(1):norm_id(2),c([2,2],.true.),0,0,0)
                do kd = norm_id(1),norm_id(2)
                    res(:,:,kd-norm_id(1)+1,3) = eq_1%q_saf_E(kd,0)*&
                        &eq_2%g_E(:,:,kd,c([3,3],.true.),0,0,0)
                    res(:,:,kd-norm_id(1)+1,4) = &
                        &sqrt(eq_2%g_E(:,:,kd,c([2,2],.true.),0,0,0) + &
                        &eq_1%q_saf_E(kd,0)**2*&
                        &eq_2%g_E(:,:,kd,c([3,3],.true.),0,0,0))
                end do
                do id = 1,4
                    do kd = norm_id(1),norm_id(2)
                        res(:,:,kd-norm_id(1)+1,id) = &
                            &res(:,:,kd-norm_id(1)+1,id)/&
                            &eq_2%jac_E(:,:,kd,0,0,0)
                    end do
                end do
        end select
        
        ! transform them to Flux coord. system
        res2 = 0._dp
        do id = 1,3
            do kd = 1,3
                res2(:,:,:,id) = res2(:,:,:,id) + &
                    &res(:,:,:,kd) * eq_2%T_FE(:,:,norm_id(1):&
                    &norm_id(2),c([id,kd],.false.),0,0,0)
            end do
        end do
        res2(:,:,:,4) = res(:,:,:,4)
        
        ! set up plot variables for calculated values
        do id = 1,4
            ! user output
            if (id.lt.4) then
                call writo('Testing B_F_'//trim(i2str(id)))
            else
                call writo('Testing B_F')
            end if
            call lvl_ud(1)
            
            ! set some variables
            if (id.lt.4) then
                file_name = 'TEST_B_F_'//trim(i2str(id))
                description = 'Testing calculated with given value for B_F_'//&
                    &trim(i2str(id))
            else
                file_name = 'TEST_B_F'
                description = 'Testing calculated with given value for B_F'
            end if
            
            ! plot difference
            if (id.lt.4) then
                call plot_diff_HDF5(res2(:,:,:,id),&
                    &eq_2%g_F(:,:,norm_id(1):norm_id(2),c([3,id],.true.),0,0,0)&
                    &/eq_2%jac_F(:,:,norm_id(1):norm_id(2),0,0,0),file_name,&
                    &tot_dim,loc_offset,description,output_message=.true.)
            else
                call plot_diff_HDF5(res2(:,:,:,id),sqrt(&
                    &eq_2%g_F(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,0,0))&
                    &/eq_2%jac_F(:,:,norm_id(1):norm_id(2),0,0,0),file_name,&
                    &tot_dim,loc_offset,description,output_message=.true.)
            end if
            
            call lvl_ud(-1)
        end do
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_B_F
    
    !> Performs tests on pressure balance.
    !!
    !!  \f[\mu_0 \frac{\partial p}{\partial u^2} = \frac{1}{\mathcal{J}}
    !!      \left(\frac{\partial B_2}{\partial u^3} -
    !!      \frac{\partial B_3}{\partial u^2}\right)\f]
    !!  \f[\mu_0 \mathcal{J} \frac{\partial p}{\partial u^3} = 0 \rightarrow
    !!      \left(\frac{\partial B_1}{\partial u^3} =
    !!      \frac{\partial B_3}{\partial u^1}\right), \f]
    !! working in the (modified) Flux coordinates
    !!      \f$\left(\alpha,\psi,\theta\right)_\text{F}\f$
    !!
    !! \note Debug version only
    !!
    !! \return ierr
    integer function test_p(grid_eq,eq_1,eq_2) result(ierr)
        use num_utilities, only: c, spline
        use grid_utilities, only: trim_grid
        use eq_vars, only: vac_perm
        use num_vars, only: eq_style, norm_disc_prec_eq
        use HELENA_vars, only: RBphi_H, h_H_11, h_H_12
        
        character(*), parameter :: rout_name = 'test_p'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: norm_id_f(2)                                                 ! norm_id transposed to full grid
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        integer :: id, jd, kd                                                   ! counters
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test the consistency of the equilibrium variables &
            &with the given pressure')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,2))
        
        ! user output
        call writo('Checking if mu_0 D2 p = 1/J (D3 B_2 - D2_B3)')
        call lvl_ud(1)
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! calculate mu_0 D2p
        res(:,:,:,1) = &
            &(eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,1) - &
            &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,1,0))/&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0) - &
            &(eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,0)*&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,1) - &
            &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,0,0)*&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,1,0))/ &
            &(eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0)**2)
        res(:,:,:,1) = res(:,:,:,1)/eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0)
        !!call plot_HDF5('var','TEST_Dg_FD_23',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,1))
        !!call plot_HDF5('var','TEST_g_FD_23',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,0))
        !!call plot_HDF5('var','TEST_Dg_FD_33',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,1,0))
        !!call plot_HDF5('var','TEST_g_FD_33',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,0,0))
        !!call plot_HDF5('var','TEST_Dg_FD',eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,1))
        !!call plot_HDF5('var','TEST_g_FD',eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0))
        
        ! save mu_0 D2p in res
        do kd = norm_id(1),norm_id(2)
            res(:,:,kd-norm_id(1)+1,2) = vac_perm*eq_1%pres_FD(kd,1)
        end do
        
        ! set some variables
        file_name = 'TEST_D2p'
        description = 'Testing whether mu_0 D2 p = 1/J (D3 B_2 - D2_B3)'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
            &loc_offset,description,output_message=.true.)
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Checking if mu_0 J D1p = 0 => D3 B_1 = D2 B_3')
        call lvl_ud(1)
        
        ! calculate D3 B1
        res(:,:,:,1) = &
            &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([1,3],.true.),0,0,1)/&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0) - &
            &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([1,3],.true.),0,0,0)*&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,1) / &
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0)**2
        
        ! calculate D1 B3
        res(:,:,:,2) = &
            &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),1,0,0)/&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0) - &
            &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,0,0)*&
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),1,0,0) / &
            &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0)**2
        
        ! set some variables
        file_name = 'TEST_D1p'
        description = 'Testing whether D3 B_1 = D1 B_3'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
            &loc_offset,description,output_message=.true.)
        
        call lvl_ud(-1)
        
        ! extra testing if Helena
        if (eq_style.eq.2) then
            ! set norm_id_f for quantities tabulated on full grid
            norm_id_f = grid_eq%i_min+norm_id-1
            
            ! user output
            call writo('Checking if D2 B_3 |F = D1 (qF+qh_11/F) |H')
            call lvl_ud(1)
            
            ! calculate if D2 B_3 = D1 (qF) + D1 (q/F h_11)
            res(:,:,:,1) = &
                &(eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,1,0)/&
                &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0) - &
                &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,0,0)*&
                &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,1,0)/ &
                &(eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0)**2))
            do jd = 1,grid_trim%n(2)
                do id = 1,grid_trim%n(1)
                    ierr = spline(grid_trim%loc_r_E,&
                        &eq_1%q_saf_E(norm_id(1):norm_id(2),0)*&
                        &RBphi_H(norm_id_f(1):norm_id_f(2),0)+&
                        &eq_1%q_saf_E(norm_id(1):norm_id(2),0)*&
                        &h_H_11(id,norm_id_f(1):norm_id_f(2))/&
                        &RBphi_H(norm_id_f(1):norm_id_f(2),0),&
                        &grid_trim%loc_r_E,res(id,jd,:,2),&
                        &ord=norm_disc_prec_eq,deriv=1)
                    CHCKERR('')
                end do
            end do
            
            ! set some variables
            file_name = 'TEST_D2B_3'
            description = 'Testing whether D2 B_3 |F = D1 (qF+qh_11/F) |H'
            
            ! plot difference
            call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
                &loc_offset,description,output_message=.true.)
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Checking if D3 B_2 |F = -q/F D2 h_11 + F q'' |H')
            call lvl_ud(1)
            
            ! calculate if D3 B_2 = -q/F D2 h_11 + F D1 q
            res(:,:,:,1) = &
                &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,1)/&
                &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0) - &
                &eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,0)*&
                &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,1)/ &
                &(eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0)**2)
            do kd = norm_id(1),norm_id(2)
                do jd = 1,grid_trim%n(2)
                    ierr = spline(grid_eq%theta_E(:,jd,kd),&
                        &h_H_12(:,kd+grid_eq%i_min-1),grid_eq%theta_E(:,jd,kd),&
                        &res(:,jd,kd-norm_id(1)+1,2),&
                        &ord=norm_disc_prec_eq,deriv=1)
                    CHCKERR('')
                end do
                res(:,:,kd-norm_id(1)+1,2) = &
                    &-eq_1%q_saf_E(kd,0)/RBphi_H(kd+grid_eq%i_min-1,0) * &
                    &res(:,:,kd-norm_id(1)+1,2) + &
                    &RBphi_H(kd+grid_eq%i_min-1,0) * &
                    &eq_1%q_saf_E(kd,1)
            end do
            
            ! set some variables
            file_name = 'TEST_D3B_2'
            description = 'Testing whether D3 B_2 |F = -D2 (qh_12/F) + f q'' |H'
            
            ! plot difference
            call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
                &loc_offset,description,output_message=.true.)
            
            call lvl_ud(-1)
        end if
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_p
#endif
end module eq_ops
