!------------------------------------------------------------------------------!
!   Operations on the equilibrium variables                                    !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln, max_deriv
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type, eq_2_type
    use num_utilities, only: check_deriv
    
    implicit none
    private
    public calc_eq, calc_derived_q, calc_normalization_const, normalize_input, &
        &print_output_eq, flux_q_plot, redistribute_output_eq, divide_eq_jobs, &
        &calc_eq_jobs_lims, calc_T_HF, B_plot, J_plot, kappa_plot, delta_r_plot
#if ldebug
    public debug_calc_derived_q, debug_create_VMEC_input, debug_J_plot
#endif
    
    ! global variables
    integer :: fund_n_par                                                       ! fundamental interval width
#if ldebug
    logical :: debug_calc_derived_q = .false.                                   ! plot debug information for calc_derived_q
    logical :: debug_J_plot = .false.                                           ! plot debug information for J_plot
    logical :: debug_create_VMEC_input = .false.                                ! plot debug information for create_VMEC_input
#endif
    
    ! interfaces
    interface calc_eq
        module procedure calc_eq_1, calc_eq_2
    end interface
    interface print_output_eq
        module procedure print_output_eq_1, print_output_eq_2
    end interface
    interface redistribute_output_eq
        module procedure redistribute_output_eq_1, redistribute_output_eq_2
    end interface
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface
    interface calc_g_C
        module procedure calc_g_C_ind, calc_g_C_arr
    end interface
    interface calc_g_V
        module procedure calc_g_V_ind, calc_g_V_arr
    end interface
    interface calc_h_H
        module procedure calc_h_H_ind, calc_h_H_arr
    end interface
    interface calc_g_F
        module procedure calc_g_F_ind, calc_g_F_arr
    end interface
    interface calc_jac_C
        module procedure calc_jac_C_ind, calc_jac_C_arr
    end interface
    interface calc_jac_V
        module procedure calc_jac_V_ind, calc_jac_V_arr
    end interface
    interface calc_jac_H
        module procedure calc_jac_H_ind, calc_jac_H_arr
    end interface
    interface calc_jac_F
        module procedure calc_jac_F_ind, calc_jac_F_arr
    end interface
    interface calc_T_VC
        module procedure calc_T_VC_ind, calc_T_VC_arr
    end interface
    interface calc_T_VF
        module procedure calc_T_VF_ind, calc_T_VF_arr
    end interface
    interface calc_T_HF
        module procedure calc_T_HF_ind, calc_T_HF_arr
    end interface

contains
    ! calculate  the equilibrium  quantities on  a grid  determined by  straight
    ! field lines. This grid has the dimensions (n_par,loc_n_r).
    ! Optionally, for eq_2, the used variables can be deallocated on the fly, to
    ! limit memory usage.
    integer function calc_eq_1(grid_eq,eq) result(ierr)                         ! flux version
        use num_vars, only: eq_style, rho_style, use_normalization, &
            &use_pol_flux_E, use_pol_flux_F
        use grid_utilities, only: apply_disc
        use eq_vars, only: rho_0
        use num_utilities, only: derivs
        use eq_utilities, only: calc_F_derivs
        
        character(*), parameter :: rout_name = 'calc_eq_1'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(inout) :: eq                                    ! flux equilibrium variables
        
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
    integer function calc_eq_2(grid_eq,eq_1,eq_2,dealloc_vars) result(ierr)     ! metric version
        use num_vars, only: eq_style, export_HEL, use_normalization
        use num_utilities, only: derivs, c
        use eq_utilities, only: calc_inv_met, calc_F_derivs
        use VMEC_utilities, only: calc_trigon_factors
#if ldebug
        use num_vars, only: ltest
        use input_utilities, only: get_log, pause_prog
        use HELENA_ops, only: test_metrics_H
#endif
        
        character(*), parameter :: rout_name = 'calc_eq_2'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! metric equilibrium variables
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium variables
        logical, intent(in), optional :: dealloc_vars                           ! deallocate variables on the fly after writing
        
        ! local variables
        integer :: id
        integer :: pmone                                                        ! plus or minus one
        logical :: dealloc_vars_loc                                             ! local dealloc_vars
        
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
                ! derivatives
                call writo('Calculate R,Z,L...')
                ierr = calc_trigon_factors(grid_eq%theta_E,grid_eq%zeta_E,&
                    &grid_eq%trigon_factors)
                CHCKERR('')
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
                
                ! calculate the jacobian in the cylindrical coordinate system
                call writo('Calculate jac_C')
                do id = 0,max_deriv
                    ierr = calc_jac_C(eq_2,derivs(id))
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
                    ierr = calc_jac_V(eq_2,derivs(id))
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
                    deallocate(eq_2%R_E,eq_2%Z_E,eq_2%L_E)
                    deallocate(eq_2%g_C,eq_2%jac_C)
                    deallocate(eq_2%T_VC,eq_2%det_T_VC)
                end if
            case (2)                                                            ! HELENA
#if ldebug
                if (ltest) then
                    call writo('Test consistency of metric factors?')
                    if(get_log(.false.,ind=.true.)) then
                        ierr = test_metrics_H()
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
                call writo('Calculate h_H')
                do id = 0,max_deriv
                    ierr = calc_h_H(grid_eq,eq_2,derivs(id))
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
                
                ! calculate the inverse g_H of the metric factors h_H
                call writo('Calculate g_H')
                do id = 0,max_deriv
                    ierr = calc_inv_met(eq_2%g_E,eq_2%h_E,derivs(id))
                    CHCKERR('')
                end do
                
                ! export for VMEC port
                if (export_HEL) then
                    call writo('Exporting HELENA equilibrium for VMEC porting')
                    call lvl_ud(1)
                    ierr = create_VMEC_input(eq_1,eq_2)
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
        
        !!! limit Jacobian to small value to avoid infinities
        !!if (maxval(eq_2%jac_F(:,:,:,0,0,0)).gt.0._dp) then
            !!eq_2%jac_F(:,:,:,0,0,0) = max(eq_2%jac_F(:,:,:,0,0,0),tol_zero)
        !!else
            !!eq_2%jac_F(:,:,:,0,0,0) = min(eq_2%jac_F(:,:,:,0,0,0),-tol_zero)
        !!end if
        
        ! possibly deallocate
        if (dealloc_vars_loc) then
            deallocate(eq_2%g_E,eq_2%jac_E)
        end if
        
        ! Transform metric equilibrium E into F derivatives
        ierr = calc_F_derivs(eq_2)
        CHCKERR('')
        
        ! possibly deallocate
        if (dealloc_vars_loc) then
            deallocate(eq_2%g_F,eq_2%h_F,eq_2%jac_F)
            deallocate(eq_2%T_EF,eq_2%T_FE,eq_2%det_T_EF,eq_2%det_T_FE)
        end if
        
        ! Calculate derived metric quantities
        call calc_derived_q(grid_eq,eq_1,eq_2)
        
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
    contains
        ! Plots flux quantities in file for VMEC port.
        ! Optionally, a  perturbation can be  added: Either the  displacement of
        ! the plasma position can be described  (pert_style 1), or ripple in the
        ! toroidal magnetic  field (pert_style  2), with  a fixed  toroidal mode
        ! number.
        ! Both perturbation styles can have various prescription types:
        !   1. file with Fourier modes in the HELENA angular coordinates
        !   2. same but manually
        !   3. file with perturbation on general points R, Z, phi
        ! For perturbation  style 2, an  additional fourth prescription  type is
        ! through a map of the magnetic perturbation in space.
        ! Also,  for  this  perturbation  style,  a  file  has  to  be  provided
        ! that describes  the translation  between position ripple  and magnetic
        ! perturbation.  This file  can  be generated  for  an already  existing
        ! ripple case using POST with --compare_tor_pos with n_zeta_plot = 2 and
        ! min_theta_plot and max_theta_plot indicating half a ripple period.
        ! After  translating  these to  an  equivalent  position deviation,  the
        ! procedure calculates the  NUFFT in the geometrical angles  in order to
        ! output it to a VMEC input file.
        ! The output from this VMEC run can then be used to iteratively create a
        ! new  file to  translate  toroidal magnetic  field  ripple to  position
        ! perturbation.
        ! A note about the indices of B_F, B_F_pert and B_F_loc:
        !   B_F_loc:  (pol modes, cos/sin)
        !   B_F:      (pol modes, tor modes, cos/sin (m theta), R/Z)
        !   B_F_pert: (pol modes, cos/sin (m theta), R/Z, cos/sin (N zeta))
        integer function create_VMEC_input(eq_1,eq_2) result(ierr)
            use eq_vars, only: pres_0, R_0, psi_0
            use grid_vars, only: n_r_eq
            use grid_utilities, only: setup_interp_data, apply_disc, nufft
            use HELENA_vars, only: nchi, R_H, Z_H, ias, chi_H
            use X_vars, only: min_r_sol, max_r_sol
            use input_utilities, only: pause_prog, get_log, get_int, get_real
            use num_utilities, only: GCD, bubble_sort
            use num_vars, only: eq_name, HEL_pert_i, HEL_export_i, &
                &norm_disc_prec_eq, prop_B_tor_i
            use files_utilities, only: skip_comment
#if ldebug
            use grid_utilities, only: setup_deriv_data
#endif
            
            character(*), parameter :: rout_name = 'create_VMEC_input'
            
            ! input / output
            type(eq_1_type), intent(in) :: eq_1                                 ! flux equilibrium quantities
            type(eq_2_type), intent(in) :: eq_2                                 ! metric equilibrium quantities
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            integer :: jd_min                                                   ! start value of jd, possibly excluding N = 0 from equilibrium
            integer :: n_B                                                      ! nr. of points in Fourier series
            integer :: nfp                                                      ! scale factor for toroidal mode numbers
            integer :: tot_nr_pert                                              ! total number of perturbations combinations (N,M)
            integer :: nr_n                                                     ! number of different N
            integer :: n_loc, m_loc                                             ! local n_pert and m_pert
            integer :: n_id                                                     ! index in bundled n_pert
            integer :: istat                                                    ! status
            integer :: plot_dim(2)                                              ! plot dimensions
            integer :: rec_min_m                                                ! recommended minimum of poloidal modes
            integer :: pert_style                                               ! style of perturbation (1: r, 2: B_tor)
            integer :: pert_type                                                ! type of perturbation prescription
            integer :: r_prop                                                   ! r at which to use proportionality factor
            integer :: max_n_B_output                                           ! max. nr. of modes written in output file (constant in VMEC)
            integer :: n_prop_B_tor(2)                                          ! n of prop_B_tor
            integer :: flux_type_prop_B_tor                                     ! 1 (poloidal) or 2 (toroidal) flux used in prop_B_tor
            integer, allocatable :: n_pert(:), m_pert(:,:)                      ! tor. mode numbers and pol. mode numbers for each of them
            integer, allocatable :: m_pert_copy(:,:)                            ! copy of m_pert, for sorting
            integer, allocatable :: piv(:)                                      ! pivots for sorting
            integer, allocatable :: nr_m(:)                                     ! number of different M for every N
            integer, allocatable :: nr_m_copy(:)                                ! copy of nr_m, for sorting
            character(len=1) :: loc_data_char                                   ! first data character
            character(len=2) :: plus                                            ! "+ " or ""
            character(len=8) :: flux_name(2)                                    ! "poloidal" or "toroidal"
            character(len=6) :: path_prefix = '../../'                          ! prefix of path
            character(len=max_str_ln) :: dummy_string                           ! dummy string
            character(len=max_str_ln) :: HEL_pert_file_name                     ! name of perturbation file
            character(len=max_str_ln) :: err_msg                                ! error message
            character(len=max_str_ln) :: file_name                              ! name of file
            character(len=max_str_ln) :: plot_name(2)                           ! name of plot file
            character(len=max_str_ln) :: plot_title(2)                          ! name of plot
            character(len=max_str_ln) :: prop_B_tor_file_name                   ! name of B_tor proportionality file
            real(dp) :: m_tol = 1.E-7_dp                                        ! tolerance for Fourier mode strength
            real(dp) :: delta_loc(2)                                            ! local delta
            real(dp) :: plot_lims(2,2)                                          ! limits of plot dims [pi]
            real(dp) :: norm_B_H                                                ! normalization for R and Z Fourier modes
            real(dp) :: mult_fac                                                ! global multiplication factor
            real(dp) :: min_theta_prop_B_tor                                    ! starting angle of prop_B_tor
            real(dp) :: RZ_B_0(2)                                               ! origin of R and Z of boundary
            real(dp) :: dummy_input(4)                                          ! dummy variable
            real(dp) :: s_V(99)                                                 ! normal coordinate s for writing
            real(dp) :: pres_V(99)                                              ! pressure for writing
            real(dp) :: rot_T_V(99)                                             ! rotational transform for writing
            real(dp), allocatable :: R_H_loc(:,:)                               ! local R_H
            real(dp), allocatable :: Z_H_loc(:,:)                               ! local Z_H
            real(dp), allocatable :: delta(:,:,:)                               ! amplitudes of perturbations (N,M,c/s)
            real(dp), allocatable :: delta_copy(:,:,:)                          ! copy of m_pert, for sorting
            real(dp), allocatable :: BH_0(:,:)                                  ! R and Z at unperturbed bounday
            real(dp), allocatable :: BH_pert(:,:)                               ! BH perturbed by cos(M theta) or sin(M theta)
            real(dp), allocatable :: BH_deriv(:,:)                              ! theta derivative of R and Z
            real(dp), allocatable :: B_F(:,:,:,:)                               ! cos and sin Fourier components
            real(dp), allocatable :: B_F_pert(:,:,:,:)                          ! cosine and sine of fourier series of perturbation
            real(dp), allocatable :: B_F_loc(:,:)                               ! Local B_F
            real(dp), allocatable :: u_norm(:,:)                                ! normalized unit vector
            real(dp), allocatable :: theta(:,:)                                 ! pol. angle: (geometric, HELENA (equidistant))
            real(dp), allocatable :: prop_B_tor(:,:,:)                          ! proportionality between delta B_tor and delta_norm
            real(dp), allocatable :: prop_B_tor_F_loc(:,:)                      ! local fourier coefficients of prop_B_tor
            real(dp), allocatable :: s_prop_B_tor(:)                            ! normal positions at which prop_B_tor is tabulated
            logical :: zero_N_pert                                              ! there is a perturbation with N = 0
            logical :: pert_eq                                                  ! whether equilibrium is perturbed
            logical :: stel_sym                                                 ! whether there is stellarator symmetry
            logical :: change_max_n_B_output                                    ! whether to change max_n_B_output
            logical :: found                                                    ! whether something was found
            logical :: choose_norm_pos                                          ! whether to choose a different normal position to use
            type(disc_type) :: norm_interp_data                                 ! data for normal interpolation
#if ldebug
            real(dp), allocatable :: BH_0_ALT(:,:)                              ! reconstructed R and Z
            type(disc_type) :: deriv_data                                       ! data for derivatives in theta
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
            
            ! set up local R_H and Z_H
            allocate(R_H_loc(nchi,n_r_eq))
            allocate(Z_H_loc(nchi,n_r_eq))
            R_H_loc = R_H
            Z_H_loc = Z_H
            if (use_normalization) then
                R_H_loc = R_H*R_0
                Z_H_loc = Z_H*R_0
            end if
            
            ! set up auxiliary variable
            flux_name = ['poloidal','toroidal']
            
            ! get input about equilibrium perturbation
            ! Note: negative N values will be converted to negative M values and
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
                        open(prop_B_tor_i,FILE=trim(prop_B_tor_file_name),&
                            &IOSTAT=ierr,STATUS='old')
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
                    ierr = skip_comment(prop_B_tor_i,name=prop_B_tor_file_name)
                    CHCKERR('')
                    read(prop_B_tor_i,*,IOSTAT=ierr) dummy_input
                    err_msg = 'failed to read dummy_input'
                    CHCKERR(err_msg)
                    n_prop_B_tor = nint(dummy_input(1:2))
                    min_theta_prop_B_tor = dummy_input(3)
                    flux_type_prop_B_tor = 1
                    if (dummy_input(4).lt.1) flux_type_prop_B_tor = 2
                    call writo(trim(i2str(n_prop_B_tor(1)))//' poloidal points')
                    call writo(trim(i2str(n_prop_B_tor(2)))//' normal points &
                        &in '//flux_name(flux_type_prop_B_tor)//' flux')
                    call lvl_ud(1)
                    allocate(prop_B_tor(n_prop_B_tor(1),n_prop_B_tor(2),2))
                    allocate(s_prop_B_tor(n_prop_B_tor(2)))
                    ierr = skip_comment(prop_B_tor_i,name=prop_B_tor_file_name)
                    CHCKERR('')
                    do kd = 1,n_prop_B_tor(2)
                        do id = 1,2
                            read(prop_B_tor_i,*,IOSTAT=ierr) dummy_input(id), &
                                &prop_B_tor(:,kd,id)
                            err_msg = 'failed to read dummy_input'
                            CHCKERR(err_msg)
                        end do
                        if (abs(dummy_input(1)-dummy_input(2)).lt.&
                            &10*tiny(1._dp)) then
                            s_prop_B_tor(kd) = dummy_input(1)
                            call writo(trim(i2str(kd))//'/'//&
                                &trim(i2str(n_prop_B_tor(2)))//': s = '//&
                                &trim(r2str(s_prop_B_tor(kd))))
                        else
                            ierr = 1
                            err_msg = 'Check s-values in "'//&
                                &trim(prop_B_tor_file_name)//'"'
                            CHCKERR(err_msg)
                        end if
                    end do
                    call lvl_ud(-1)
                    
                    call writo('Use last normal position?')
                    choose_norm_pos = get_log(.true.)
                    if (.not.choose_norm_pos) then
                        call writo('Normal position to use, in normalized '//&
                            &flux_name(flux_type_prop_B_tor)//' flux')
                        r_prop = get_int(lim_lo=1,lim_hi=n_prop_B_tor(2),&
                            &ind=.true.)
                    else
                        r_prop = n_prop_B_tor(2)
                    end if
                    
                    call lvl_ud(-1)
                    
                    call writo('toroidal modulation mode number?')
                    n_loc = get_int()
                end if
                
                ! get perturbation type
                call writo('How do you want to prescribe the perturbation?')
                call lvl_ud(1)
                call writo('1: through a file with Fourier modes in HELENA &
                    &coordinates')
                call lvl_ud(1)
                call writo('It should provide N M delta_cos delta_sin in four &
                    &columns')
                call writo('Rows starting with "#" are ignored')
                call lvl_ud(-1)
                call writo('2: same as 1 but interactively')
                call writo('3: through a file with Fourier modes in general &
                    &coordinates')
                call lvl_ud(1)
                call writo('It should provide R Z phi delta_r in four columns')
                call writo('Rows starting with "#" are ignored')
                call writo('Note: You can use the output of DESCUR in the &
                    &STELLOPT suite')
                call lvl_ud(-1)
                call writo('4: same as 3 but interactively')
                pert_type = get_int(lim_lo=1,lim_hi=4)
                
                select case (pert_type)
                    case (1,3)                                                  ! file with Fourier modes
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
                        
                        if (pert_type.eq.1) then
                            ! get number of lines
                            istat = 0
                            tot_nr_pert = 0
                            do while (istat.eq.0)
                                ! read first character of data
                                read(HEL_pert_i,*,IOSTAT=istat) loc_data_char
                                if (istat.eq.0 .and. loc_data_char.ne.'#') then ! exclude comment lines
                                    tot_nr_pert = tot_nr_pert + 1
                                end if
                            end do
                            rewind(HEL_pert_i)
                            loc_data_char = '#'
                        else
                            ! set error message
                            err_msg = 'Cannot read the file "'//&
                                &trim(HEL_pert_file_name)//'"'
                            
                            ! forward until 'MB'
                            found = .false.
                            do while (.not.found) 
                                read(HEL_pert_i,*,IOSTAT=ierr) dummy_string
                                if (trim(dummy_string).eq.'MB') then
                                    found = .true.
                                    tot_nr_pert = 0
                                end if
                                CHCKERR(err_msg)
                            end do
                            
                            ! count number of lines
                            found = .false.
                            do while (.not.found)
                                read(HEL_pert_i,'(A)',IOSTAT=ierr) dummy_string
                                if (trim(dummy_string).eq.'') then
                                    found = .true.
                                    write(*,*) ''
                                else
                                    tot_nr_pert = tot_nr_pert + 1
                                end if
                                err_msg = 'Cannot read the file "'//&
                                    &trim(HEL_pert_file_name)//'"'
                                CHCKERR(err_msg)
                            end do
                            do id = 1,tot_nr_pert+1
                                backspace(HEL_pert_i)
                            end do
                        end if
                    case (2,4)                                                  ! interactively
                        ! user input
                        call writo('Prescribe interactively')
                        call writo('How many different combinations of &
                            &toroidal and poloidal mode numbers (N,M)?')
                        tot_nr_pert = get_int(lim_lo=1)
                end select
                
                ! ask for global multiplication factor
                call writo('Globally multiply with some factor?')
                call lvl_ud(1)
                mult_fac = 0.01_dp
                call writo('Current multiplication factor: '//&
                    &trim(r2strt(mult_fac))//'m')
                call lvl_ud(1)
                call writo('minor radius: '//&
                    &trim(r2strt((maxval(R_H_loc)-minval(R_H_loc))/2))//'m')
                call writo('major radius: '//&
                    &trim(r2strt((maxval(R_H_loc)+minval(R_H_loc))/2))//'m')
                call lvl_ud(-1)
                call lvl_ud(-1)
                if (get_log(.false.)) mult_fac = get_real()
                
                ! set up n, m and delta  (use maximum possible size, possibly +1
                ! if N=0 is not a perturbation, and limit afterwards)
                allocate(n_pert(tot_nr_pert+1))
                allocate(m_pert(tot_nr_pert+1,tot_nr_pert+1))
                allocate(delta(tot_nr_pert+1,tot_nr_pert+1,2))
                allocate(nr_m(tot_nr_pert+1))
                nr_n = 0
                nr_m = 0
                do jd = 1,tot_nr_pert
                    select case (pert_type)
                        case (1)                                                ! file with Fourier modes in HELENA coordinates
                            err_msg = 'Could not read file '//&
                                &trim(HEL_pert_file_name)
                            do while (loc_data_char.eq.'#')
                                ! read first character of data
                                read(HEL_pert_i,*,IOSTAT=ierr) &
                                    &loc_data_char
                                CHCKERR(err_msg)
                                if (loc_data_char.ne.'#') then                  ! exclude comment lines
                                    backspace(UNIT=HEL_pert_i)                  ! go back one line
                                    read(HEL_pert_i,*,IOSTAT=ierr) &
                                        &n_loc, m_loc, delta_loc
                                    CHCKERR(err_msg)
                                    delta_loc = delta_loc*mult_fac
                                end if
                            end do
                            loc_data_char = '#'
                            
                            close(HEL_pert_i)
                        case (2)                                                ! same as 1 but interactively
                            call writo('For mode '//trim(i2str(jd))//'/'//&
                                &trim(i2str(tot_nr_pert))//':')
                            call lvl_ud(1)
                            
                            ! input
                            call writo('Toroidal mode number N?')
                            n_loc = get_int()
                            call writo('Poloidal mode number M?')
                            m_loc = get_int()
                            call writo('Perturbation strength ~ cos('//&
                                &trim(i2str(m_loc))//' theta - '//&
                                &trim(i2str(n_loc))//' zeta)?')
                            delta_loc(1) = get_real()
                            call writo('Perturbation strength ~ sin('//&
                                &trim(i2str(m_loc))//' theta - '//&
                                &trim(i2str(n_loc))//' zeta)?')
                            delta_loc(2) = get_real()
                            
                            delta_loc = delta_loc*mult_fac
                            
                            call lvl_ud(-1)
                        case (3)                                                ! file with Fourier modes in general coordinates
                            ! user input
                            ierr = 1
                            CHCKERR('NOT YET IMPLEMENTED')
                        case (4)                                                ! same as 3 but interactively
                            m_loc = jd/2
                            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                            write(*,*) '!!!! THIS IS NOT CORRECT YET !!!!!!!!!!'
                            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                            if (m_loc.ne.(jd-1)/2) m_loc = - m_loc
                            delta_loc = mult_fac/prop_B_tor_F_loc(jd/2+1,:)
                            if (m_loc.ne.0) delta_loc = delta_loc*0.5_dp        ! modes with nonzero m are counted double
                    end select
                    
                    ! has N = 0 been prescribed?
                    if (n_loc.eq.0) zero_N_pert = .true.
                    
                    ! possibly convert to positive N
                    if (n_loc.lt.0) then
                        n_loc = -n_loc
                        m_loc = -m_loc
                        delta_loc(2) = -delta_loc(2)
                    end if
                    
                    ! bundle into n_pert and m_pert
                    n_id = nr_n+1                                               ! default new N in next index
                    do id = 1,nr_n
                        if (n_loc.eq.n_pert(id)) n_id = id
                    end do
                    if (n_id.gt.nr_n) then                                      ! new N
                        nr_n = n_id                                             ! increment number of N
                        n_pert(n_id) = n_loc                                    ! with value n_loc
                    end if
                    nr_m(n_id) = nr_m(n_id) + 1                                 ! one more M in either old or new N
                    m_pert(n_id,nr_m(n_id)) = m_loc                             ! with value m_loc
                    delta(n_id,nr_m(n_id),:) = delta_loc                        ! and perturbation amplitude for cos and sin with value delta_loc
                end do
            else
                nr_n = 0
            end if
            
            ! plot properties
            plot_dim = [100,100]
            plot_lims(:,1) = [0.0_dp,2.0_dp]
            plot_lims(:,2) = [0.5_dp,2.0_dp]
            call writo('Change plot properties from defaults?')
            call lvl_ud(1)
            do id = 1,2
                call writo(trim(i2str(plot_dim(id)))//' '//flux_name(id)//&
                    &' points on range '//trim(r2strt(plot_lims(1,id)))//&
                    &' pi .. '//trim(r2strt(plot_lims(2,id)))//' pi')
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
            
            ! possibly include N = 0 for equilibrium if not already done so
            if (.not.zero_N_pert) then                                          ! N = 0 was not already included
                jd_min = 2
                nr_n = nr_n + 1
                if (pert_eq) then
                    n_pert(nr_n) = 0
                    nr_m(nr_n) = 1
                    m_pert(nr_n,1) = 0
                    delta(nr_n,1,:) = 0
                end if
            else
                jd_min = 1
            end if
            
            ! sort and output
            if (pert_eq) then
                ! sort
                allocate(m_pert_copy(size(m_pert,1),size(m_pert,2)))
                allocate(delta_copy(size(delta,1),size(delta,2),size(delta,3)))
                allocate(nr_m_copy(size(nr_m)))
                m_pert_copy = m_pert
                delta_copy = delta
                nr_m_copy = nr_m
                allocate(piv(nr_n))
                call bubble_sort(n_pert(1:nr_n),piv)
                do jd = 1,nr_n
                    nr_m(jd) = nr_m_copy(piv(jd))
                    m_pert(jd,:) = m_pert_copy(piv(jd),:)
                    delta(jd,:,:) = delta_copy(piv(jd),:,:)
                    !!write(*,*) 'jd', jd, nr_n
                    !!write(*,*) 'nr_m', nr_m(jd)
                    !!write(*,*) ' >> n', n_pert(jd)
                    !!write(*,*) ' >> m', m_pert(jd,1:nr_m(jd))
                end do
                deallocate(m_pert_copy)
                deallocate(delta_copy)
                deallocate(nr_m_copy)
                
                ! output
                call writo('Summary of perturbation form:')
                call lvl_ud(1)
                plus = ''
                do jd = jd_min,nr_n
                    do id = 1,nr_m(jd)
                        call writo(plus//&
                            &trim(r2strt(delta(jd,id,1)))//' cos('//&
                            &trim(i2str(m_pert(jd,id)))//' theta - '//&
                            &trim(i2str(n_pert(jd)))//' zeta) + '//&
                            &trim(r2strt(delta(jd,id,2)))//' sin('//&
                            &trim(i2str(m_pert(jd,id)))//' theta - '//&
                            &trim(i2str(n_pert(jd)))//' zeta)')
                        plus = '+ '
                    end do
                end do
                call lvl_ud(-1)
            end if
            
            ! user output
            call writo('Set up axisymmetric data')
            call lvl_ud(1)
            
            ! Fourier transform of boundary R and Z
            if (ias.eq.0) then                                                  ! symmetric (so nchi is odd)
                n_B = nchi*2-2                                                  ! -2 because of overlapping points at 0, pi
            else                                                                ! asymmetric (so nchi is aumented by one)
                n_B  = nchi-1                                                   ! -1 because of overlapping points at 0
            end if
            allocate(BH_0(n_B,2))                                               ! HELENA coords (no overlap at 0)
            allocate(theta(n_B,2))                                              ! geometrical poloidal angle (overlap at 0)
            if (ias.eq.0) then                                                  ! symmetric
                BH_0(1:nchi,1) = R_H_loc(:,n_r_eq)
                BH_0(1:nchi,2) = Z_H_loc(:,n_r_eq)
                theta(1:nchi,2) = chi_H
                do kd = 1,nchi-2
                    BH_0(nchi+kd,1) = R_H_loc(nchi-kd,n_r_eq)
                    BH_0(nchi+kd,2) = -Z_H_loc(nchi-kd,n_r_eq)
                    theta(nchi+kd,2) = 2*pi-chi_H(nchi-kd)
                end do
            else
                BH_0(:,1) = R_H_loc(1:n_B,n_r_eq)
                BH_0(:,2) = Z_H_loc(1:n_B,n_r_eq)
                theta(:,2) = chi_H(1:n_B)
            end if
            RZ_B_0(1) = sum(R_H_loc(:,1))/size(R_H_loc,1)
            RZ_B_0(2) = sum(Z_H_loc(:,1))/size(Z_H_loc,1)
            theta(:,1) = atan2(BH_0(:,2)-RZ_B_0(2),BH_0(:,1)-RZ_B_0(1))
            
            !!! TEMPORARILY
            !!write(*,*) 'TEMPORARILY CIRCULAR TOKAMAK !!!!!!!!!!!'
            !!theta(:,1) = theta(:,2)
            !!BH_0(:,1) = 1.5_dp + 0.5*cos(theta(:,1))
            !!BH_0(:,2) = 0.5*sin(theta(:,1))
            
            ! plot angles
            plot_name(1) = 'chi'
            plot_title = ['theta_G','chi_H  ']
            call print_ex_2D(plot_title(1:2),plot_name(1),theta/pi,&
                &x=reshape([theta(:,2)],[n_B,1])/pi,&
                &draw=.false.)
            call draw_ex(plot_title(1:2),plot_name(1),2,1,.false.)
            
            ! plot R and Z
            plot_name(1) = 'RZ'
            plot_title = ['R_H','Z_H']
            call print_ex_2D(plot_title(1:2),plot_name(1),BH_0,&
                &x=reshape([theta(:,1)],[n_B,1])/pi,&
                &draw=.false.)
            call draw_ex(plot_title(1:2),plot_name(1),2,1,.false.)
            
            ! calculate output for unperturbed R and Z
            plot_name(1) = 'R_F'
            plot_name(2) = 'Z_F'
            do kd = 1,2
                ! user output
                call writo('analyzing '//trim(plot_name(kd)))
                call lvl_ud(1)
                
                ! NUFFT
                ierr = nufft(theta(:,1),BH_0(:,kd),B_F_loc,plot_name(kd))       ! geometrical pol. Fourier coefficients
                CHCKERR('')
                if (kd.eq.1) then
                    allocate(B_F(size(B_F_loc,1),-(nr_n-1):(nr_n-1),2,2))       ! (pol modes, tor modes, cos/sin (m theta), R/Z)
                    B_F = 0._dp
                end if
                B_F(:,0,:,kd) = B_F_loc
            
#if ldebug
                if (debug_create_VMEC_input) then
                    allocate(BH_0_ALT(size(BH_0,1),size(B_F_loc,1)))
                    
                    call writo('Comparing R or Z with reconstruction &
                        &through Fourier coefficients')
                    BH_0_ALT(:,1) = B_F_loc(1,1)
                    do id = 1,size(B_F_loc,1)-1
                        BH_0_ALT(:,id+1) = BH_0_ALT(:,id) + &
                            &B_F_loc(id+1,1)*cos(id*theta(:,2)) + &
                            &B_F_loc(id+1,2)*sin(id*theta(:,2))
                    end do
                    call print_ex_2D(['orig BH','alt BH '],'',&
                        &reshape([BH_0(:,kd),BH_0_ALT(:,size(B_F_loc,1))],&
                        &[size(BH_0,1),2]),x=&
                        &reshape([theta(:,2)],[size(BH_0,1),1]))
                    
                    call writo('Plotting Fourier approximation')
                    call print_ex_2D(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                        &'_F_series',BH_0_ALT,&
                        &x=reshape([theta(:,2)],[size(BH_0,1),1]))
                    call draw_ex(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                        &'_F_series',size(B_F_loc,1),1,.false.)
                    
                    call writo('Making animation')
                    call print_ex_2D(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                        &'_F_series_anim',BH_0_ALT,&
                        &x=reshape([theta(:,2)],[size(BH_0,1),1]))
                    call draw_ex(['alt BH'],'TEST_'//trim(plot_name(kd))//&
                        &'_F_series_anim',size(B_F_loc,1),1,.false.,&
                        &is_animated=.true.)
                    
                    deallocate(BH_0_ALT)
                end if
#endif
                
                call lvl_ud(-1)
            end do
            
            ! plot axisymetric boundary in 3D
            call plot_boundary(B_F(:,0:0,:,:),[0],'last_fs',plot_dim,plot_lims)
            
            call lvl_ud(-1)
            
            ! user output
            if (pert_eq) then
                ! Calculate the  normal unit vector on  the HELENA (equidistant)
                ! poloidal grid.
                call writo('Modify axisymmetric equilibrium')
                call lvl_ud(1)
                
                ! calculate unit normal vector  by calculating first the Fourier
                ! components in  the Helena poloidal coordinate,  and then using
                ! them to obtain the poloidal derivative  of R and Z. The normal
                ! unit vector is then given by
                !   (Z_theta e_R - R_theta e_Z)/sqrt(Z_theta^2 + R_theta^2)
                allocate(BH_deriv(n_B,2))                                       ! theta derivative of BH
                allocate(u_norm(n_B,2))                                         ! normalized unit vector: (Z_theta e_R - R_theta e_Z)/norm
                BH_deriv = 0._dp
                do kd = 1,2
                    ! NUFFT
                    ierr = nufft(theta(:,2),BH_0(:,kd),B_F_loc)                 ! HELENA pol. Fourier coefficients
                    CHCKERR('')
                    
                    ! calculate poloidal derivative using the Fourier series
                    do id = 0,size(B_F_loc,1)-1
                        BH_deriv(:,kd) = BH_deriv(:,kd) + id*(&
                            &B_F_loc(id+1,2)*cos(id*theta(:,2)) - &
                            &B_F_loc(id+1,1)*sin(id*theta(:,2)))
                    end do
                    
#if ldebug
                    if (debug_create_VMEC_input) then
                        allocate(BH_0_ALT(size(BH_0,1),1))
                        
                        call writo('Comparing derivative of R or Z with &
                            &reconstruction through Fourier coefficients')
                        ierr = setup_deriv_data(theta(:,2),deriv_data,1,1)
                        CHCKERR('')
                        ierr = apply_disc(BH_0(:,kd),deriv_data,BH_0_ALT(:,1))
                        CHCKERR('')
                        call deriv_data%dealloc()
                        call print_ex_2D(['orig BH','alt BH '],'',&
                            &reshape([BH_deriv(:,kd),BH_0_ALT(:,1)],&
                            &[size(BH_deriv,1),2]),x=&
                            &reshape([theta(:,2)],[size(BH_deriv,1),1]))
                        
                        deallocate(BH_0_ALT)
                    end if
#endif
                end do
                u_norm(:,2) = sqrt(BH_deriv(:,1)**2 + BH_deriv(:,2)**2)
                u_norm(:,1) = BH_deriv(:,2)/u_norm(:,2)
                u_norm(:,2) = -BH_deriv(:,1)/u_norm(:,2)
                deallocate(BH_deriv)
                if (debug_create_VMEC_input) then
                    call print_ex_2D(['u_norm'],'',u_norm,x=&
                        &reshape([theta(:,2)],[n_B,1]))
                    call print_ex_2D(['test_normal'],'',&
                        &transpose(reshape([BH_0(:,2)-u_norm(:,2)*0.1,&
                        &BH_0(:,2),BH_0(:,2)+u_norm(:,2)*0.1],&
                        &[size(BH_0,1),3])),x=&
                        &transpose(reshape([BH_0(:,1)-u_norm(:,1)*0.1,&
                        &BH_0(:,1),BH_0(:,1)+u_norm(:,1)*0.1],&
                        &[size(BH_0,1),3])))
                end if
            end if
            
            ! loop  over all  toroidal  modes N  and  include the  corresponding
            ! perturbation, consisting of terms of the form:
            ! delta_c cos(M theta - N zeta) + delta_s sin(M theta - N zeta) = 
            !   delta_c [cos(M theta)cos(N zeta) + sin(M theta) sin(N zeta) ] + 
            !   delta_s [sin(M theta)cos(N zeta) - cos(M theta) sin(N zeta) ],
            ! where theta is in Helena  straight-field line coordinates and zeta
            ! corresponds to the geometrical toroidal angle.
            ! The terms  dependent on theta will  shift up and down  the Fourier
            ! components of the axisymmetric R and Z, but the terms dependent on
            ! zeta will not.
            ! The final output should be on the geometrical poloidal grid, so to
            ! calculate the Fourier components of the  perturbed R and Z on this
            ! geometrical poloidal  grid interpolation from the  HELENA poloidal
            ! grid is necessary. This means  that all multiplications with terms
            ! ~ cos(M theta) and terms ~  sin(M theta) where theta is the HELENA
            ! poloidal angle are performed first, before calculating the Fourier
            ! modes in the geometrical poloidal grid through interpolation.
            ! Note that  the different poloidal  modes for equal  toroidal modes
            ! are combined.
            ! The resulting  cosine and sine  modes in the  geometrical poloidal
            ! angle are  then be combined  with the cosine  and sine modes  of N
            ! zeta. The result  will be a contribution to the  terms cos(m theta
            ! +/- N zeta) and sin(m theta +/- N zeta), described below.
            pert_N: do jd = jd_min,nr_n
                ! user output
                call writo(trim(i2str(jd-jd_min+1))//'/'//&
                    &trim(i2str(nr_n-jd_min+1))//': N = '//&
                    &trim(i2str(n_pert(jd))))
                call lvl_ud(1)
                
                ! calculate the  modal content in geometrical  poloidal angle of
                ! the perturbation.
                allocate(B_F_pert(size(B_F,1),2,2,2))                           ! (pol modes, cos/sin (m theta), R/Z, cos/sin (N zeta))
                
                ! user output
                call writo('Calculating contributions to cos('//&
                    &trim(i2str(n_pert(jd)))//' zeta):')
                call lvl_ud(1)
                
                allocate(BH_pert(n_B,2))                                        ! perturbed BH
                BH_pert = 0._dp
                do kd = 1,nr_m(jd)
                    ! user output
                    call writo('term '//trim(r2strt(delta(jd,kd,1)))//&
                        &' cos('//trim(i2str(m_pert(jd,kd)))//' theta) + '//&
                        &trim(r2strt(delta(jd,kd,2)))//&
                        &' sin('//trim(i2str(m_pert(jd,kd)))//&
                        &' theta) contributes')
                    do id = 1,n_B
                        BH_pert(id,:) = BH_pert(id,:) + u_norm(id,:)*(&
                            &delta(jd,kd,1)*cos(m_pert(jd,kd)*theta(id,2)) + &
                            &delta(jd,kd,2)*sin(m_pert(jd,kd)*theta(id,2)))
                    end do
                end do
                
                ! calculate the Fourier series of perturbed R and Z
                plot_name(1) = 'R_F_pert_cos_'//trim(i2str(jd-jd_min+1))
                plot_name(2) = 'Z_F_pert_cos_'//trim(i2str(jd-jd_min+1))
                do kd = 1,2
                    !!!call print_ex_2D('BH_pert_'//trim(i2str(kd)),&
                        !!!&'BH_pert_'//trim(i2str(kd)),BH_pert(:,kd),&
                        !!!&x=theta(:,1))
                    ! NUFFT
                    ierr = nufft(theta(:,1),BH_pert(:,kd),B_F_loc,&
                        &plot_name(kd))
                    CHCKERR('')
                    B_F_pert(:,:,kd,1) = B_F_loc
                end do
                deallocate(BH_pert)
                
                call lvl_ud(-1)
                
                ! user output
                call writo('Calculating contributions to sin('//&
                    &trim(i2str(n_pert(jd)))//' zeta):')
                call lvl_ud(1)
                
                allocate(BH_pert(n_B,2))                                        ! perturbed BH
                BH_pert = 0._dp
                do kd = 1,nr_m(jd)
                    ! user output
                    call writo('term '//trim(r2strt(delta(jd,kd,1)))//&
                        &' sin('//trim(i2str(m_pert(jd,kd)))//' theta) - '//&
                        &trim(r2strt(delta(jd,kd,2)))//&
                        &' cos('//trim(i2str(m_pert(jd,kd)))//&
                        &' theta) contributes')
                    do id = 1,n_B
                        BH_pert(id,:) = BH_pert(id,:) + u_norm(id,:)*(&
                            &delta(jd,kd,1)*sin(m_pert(jd,kd)*theta(id,2)) - &
                            &delta(jd,kd,2)*cos(m_pert(jd,kd)*theta(id,2)))
                    end do
                end do
                
                ! calculate the Fourier series of perturbed R and Z
                plot_name(1) = 'R_F_pert_sin_'//trim(i2str(jd-jd_min+1))
                plot_name(2) = 'Z_F_pert_sin_'//trim(i2str(jd-jd_min+1))
                do kd = 1,2
                    ! NUFFT
                    ierr = nufft(theta(:,1),BH_pert(:,kd),B_F_loc,&
                        &plot_name(kd))
                    CHCKERR('')
                    B_F_pert(:,:,kd,2) = B_F_loc
                end do
                deallocate(BH_pert)
                
                call lvl_ud(-1)
                
                ! We now  have the cos and  sin factors in front  of cos(N zeta)
                ! and sin(N zeta). The four combinations combine to terms
                !   cos(k theta +/- N zeta) and sin(k theta +/- N zeta)
                ! note that the poloidal modes k are always positive and nonzero
                ! whereas the toroidal mode numbers can be negative or positive.
                ! Therefore, the terms that contribute to a positive N are:
                !   - BC(:,N): 1/2 (alpha^c + beta^s)
                !   - BS(:,N): 1/2 (alpha^s - beta^c)
                !   - BC(:,-N): 1/2 (alpha^c - beta^s)
                !   - BS(:,-N): 1/2 (alpha^s + beta^c)
                ! where  alpha^c and  alpha^s result  from the  calculation with
                ! cos(M theta) and beta^c and beta^s from sin(M theta).
                
                ! user output
                call writo('Convert to contributions to')
                call lvl_ud(1)
                call writo(&
                    &'cos(k theta +/- '//trim(i2str(n_pert(jd)))//&
                    &' zeta) and sin(k theta +/- '//&
                    &trim(i2str(n_pert(jd)))//' zeta)')
                
                RZ: do kd = 1,2
                    ! terms with N
                    B_F_loc(:,1) = &
                        &0.5_dp*(B_F_pert(:,1,kd,1)+B_F_pert(:,2,kd,2))         ! ~ cos, N
                    B_F_loc(:,2) = &
                        &0.5_dp*(B_F_pert(:,2,kd,1)-B_F_pert(:,1,kd,2))         ! ~ sin, N
                    
                    ! update total B_F
                    B_F(:,jd-1,:,kd) = B_F(:,jd-1,:,kd) + B_F_loc
                    
                    ! terms with -N
                    B_F_loc(:,1) = &
                        &0.5_dp*(B_F_pert(:,1,kd,1)-B_F_pert(:,2,kd,2))         ! ~ cos, -N
                    B_F_loc(:,2) = &
                        &0.5_dp*(B_F_pert(:,2,kd,1)+B_F_pert(:,1,kd,2))         ! ~ sin, -N
                    
                    ! update total B_F
                    B_F(:,-(jd-1),:,kd) = B_F(:,-(jd-1),:,kd) + B_F_loc
                end do RZ
                
                deallocate(B_F_pert)
                
                call lvl_ud(-1)
                call lvl_ud(-1)
            end do pert_N
            
            ! plot boundary in 3D
            if (pert_eq) call plot_boundary(&
                &B_F,n_pert(1:nr_n),'last_fs_pert',plot_dim,plot_lims)
            
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
                &.lt. m_tol .and. &                                             ! R_s << R_c
                maxval(abs(B_F(:,:,1,2)))/maxval(abs(B_F(:,:,2,2))) &
                &.lt. m_tol) &                                                  ! R_c << Z_s
                &stel_sym = .true.
            if (stel_sym) then
                call writo("The equilibrium configuration contains &
                    &stellarator symmetry")
            else
                call writo("The equilibrium configuration does not contain &
                    &stellarator symmetry")
            end if
            
            ! find out how many poloidal modes would be necessary
            rec_min_m = 1
            norm_B_H = maxval(abs(B_F))
            do id = 1,size(B_F,1)
                if (maxval(abs(B_F(id,:,:,1)/norm_B_H)).gt.m_tol .or. &         ! for R
                    &maxval(abs(B_F(id,:,:,2)/norm_B_H)).gt.m_tol) &            ! for Z
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
            
            ! user output
            if (pert_eq) then
                call lvl_ud(-1)
            end if
            
            ! user output
            file_name = "input."//trim(eq_name)
            if (pert_eq) then
                select case (pert_type)
                    case (1,2)                                                  ! Fourier modes in HELENA coordinates
                        do jd = jd_min,nr_n
                            file_name = trim(file_name)//'_N'//&
                                &trim(i2str(n_pert(jd)))
                            do id = 1,nr_m(jd)
                                file_name = trim(file_name)//'M'//&
                                    &trim(i2str(m_pert(jd,id)))
                            end do
                        end do
                    case (3,4)                                                  ! Fourier modes in general coordinates
                        ! user input
                        ierr = 1
                        CHCKERR('NOT YET IMPLEMENTED')
                        file_name = trim(file_name)//'_MAP'
                end select
            end if
            call writo('Generate VMEC input file "'//trim(file_name)//'"')
            call writo('This can be used for VMEC porting')
            call lvl_ud(1)
            
            ! interpolate for VMEC output
            s_V = [((kd-1._dp)/(size(s_V)-1),kd=1,size(s_V))]
            ierr = setup_interp_data(&
                &eq_1%flux_t_E(:,0)/eq_1%flux_t_E(grid_eq%n(3),0),s_V,&
                &norm_interp_data,norm_disc_prec_eq)
            CHCKERR('')
            ierr = apply_disc(eq_1%pres_E(:,0),norm_interp_data,pres_V)
            CHCKERR('')
            if (use_normalization) pres_V = pres_V*pres_0
            ierr = apply_disc(-eq_1%rot_t_E(:,0),norm_interp_data,rot_T_V)
            CHCKERR('')
            call norm_interp_data%dealloc()
            
            ! output to VMEC input file
            open(HEL_export_i,STATUS='replace',FILE=trim(file_name),&
                &IOSTAT=ierr)
            CHCKERR('Failed to open file')
            
            write(HEL_export_i,"(A)") "!----- General Parameters -----"
            write(HEL_export_i,"(A)") "&INDATA"
            write(HEL_export_i,"(A)") "MGRID_FILE = 'NONE',"
            write(HEL_export_i,"(A)") "PRECON_TYPE = 'GMRES'"
            write(HEL_export_i,"(A)") "PREC2D_THRESHOLD = 3.E-8"
            write(HEL_export_i,"(A)") "DELT = 1.00E+00,"
            write(HEL_export_i,"(A)") "NS_ARRAY = 19, 39, 79, 159, 319, &
                &639"
            write(HEL_export_i,"(A)") "LRFP = F"
            write(HEL_export_i,"(A,L1)") "LASYM = ", .not.stel_sym
            write(HEL_export_i,"(A)") "LFREEB = F"
            write(HEL_export_i,"(A)") "NTOR = "//trim(i2str(nr_n-1))            ! -NTOR .. NTOR
            write(HEL_export_i,"(A)") "MPOL = "//trim(i2str(rec_min_m))         ! 0 .. MPOL-1
            write(HEL_export_i,"(A)") "TCON0 = 1"
            write(HEL_export_i,"(A)") "FTOL_ARRAY = 1.E-6, 1.E-6, 1.E-6, &
                &1.E-10, 1.E-14, 2.000E-18, "
            write(HEL_export_i,"(A)") "NITER = 20000, NSTEP = 200,"
            write(HEL_export_i,"(A)") "NFP = "//trim(i2str(nfp))
            if (use_normalization) then
                write(HEL_export_i,"(A)") "PHIEDGE = "//&
                    &trim(r2str(-eq_1%flux_t_E(grid_eq%n(3),0)*psi_0))
            else
                write(HEL_export_i,"(A)") "PHIEDGE = "//&
                    &trim(r2str(-eq_1%flux_t_E(grid_eq%n(3),0)))
            end if
            
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)") ""
            
            write(HEL_export_i,"(A)") "!----- Pressure Parameters -----"
            write(HEL_export_i,"(A)") "PMASS_TYPE = 'cubic_spline'"
            write(HEL_export_i,"(A)") "GAMMA =  0.00000000000000E+00"
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)",advance="no") "AM_AUX_S ="
            do kd = 1,size(s_V)
                write(HEL_export_i,"(A1,ES23.16)") " ", s_V(kd)
            end do
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)",advance="no") "AM_AUX_F ="
            do kd = 1,size(pres_V)
                write(HEL_export_i,"(A1,ES23.16)") " ", pres_V(kd)
            end do
            
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)") ""
            
            write(HEL_export_i,"(A)") &
                &"!----- Current/Iota Parameters -----"
            write(HEL_export_i,"(A)") "NCURR = 0"
            write(HEL_export_i,"(A)") "PIOTA_TYPE = 'Cubic_spline'"
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)",advance="no") "AI_AUX_S ="
            do kd = 1,size(s_V)
                write(HEL_export_i,"(A1,ES23.16)") " ", s_V(kd)
            end do
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)",advance="no") "AI_AUX_F ="
            do kd = 1,size(rot_t_V)
                write(HEL_export_i,"(A1,ES23.16)") " ", rot_t_V(kd)
            end do
            
            write(HEL_export_i,"(A)") ""
            write(HEL_export_i,"(A)") ""
            
            write(HEL_export_i,"(A)") &
                &"!----- Boundary Shape Parameters -----"
            ierr = print_mode_numbers(HEL_export_i,B_F(:,0,:,:),0,&
                &max_n_B_output)
            CHCKERR('')
            if (pert_eq) then
                do jd = 2,nr_n
                    ierr = print_mode_numbers(HEL_export_i,&
                        &B_F(:,-(jd-1),:,:),-n_pert(jd),max_n_B_output)
                    CHCKERR('')
                    ierr = print_mode_numbers(HEL_export_i,&
                        &B_F(:,jd-1,:,:),n_pert(jd),max_n_B_output)
                    CHCKERR('')
                end do
            end if
            write(HEL_export_i,"(A,ES23.16)") "RAXIS = ", RZ_B_0(1)
            write(HEL_export_i,"(A,ES23.16)") "ZAXIS = ", RZ_B_0(2)
            write(HEL_export_i,"(A)") "&END"
            
            close(HEL_export_i)
            
            ! user output
            call lvl_ud(-1)
            
            ! wait
            call pause_prog
        end function create_VMEC_input
        
        ! plots the boundary of a toroidal configuration
        subroutine plot_boundary(B,n,plot_name,plot_dim,plot_lims)
            ! input / output
            real(dp), intent(in) :: B(:,:,:,:)                                  ! cosine and sine of fourier series
            integer, intent(in) :: n(:)                                         ! toroidal mode numbers
            character(len=*), intent(in) :: plot_name                           ! name of plot
            integer, intent(in) :: plot_dim(2)                                  ! plot dimensions
            real(dp), intent(in) :: plot_lims(2,2)                              ! limits of plot dimensions [pi]
            
            ! local variables
            real(dp), allocatable :: ang_plot(:,:,:)                            ! angles of plot (theta,zeta)
            real(dp), allocatable :: XYZ_plot(:,:,:)                            ! coordinates of boundary (X,Y,Z)
            integer :: kd, id                                                   ! counters
            integer :: id_loc                                                   ! local id
            integer :: n_F                                                      ! number of toroidal modes (given by n)
            integer :: n_F_loc                                                  ! local n_F
            integer :: m_F                                                      ! number of poloidal modes (including 0)
            
            ! set variables
            allocate(ang_plot(plot_dim(1),plot_dim(2),2))                       ! theta and zeta
            allocate(XYZ_plot(plot_dim(1),plot_dim(2),3))                       ! X, Y and Z
            XYZ_plot = 0._dp
            n_F = size(n)
            m_F = size(B,1)
            
            ! create grid
            do kd = 1,plot_dim(1)
                ang_plot(kd,:,1) = plot_lims(1,1) + (kd-1._dp)/(plot_dim(1)-1)*&
                    &(plot_lims(2,1)-plot_lims(1,1))
            end do
            do kd = 1,plot_dim(2)
                ang_plot(:,kd,2) = plot_lims(1,2) + (kd-1._dp)/(plot_dim(2)-1)*&
                    &(plot_lims(2,2)-plot_lims(1,2))
            end do
            ang_plot = ang_plot*pi
            
            ! inverse Fourier transform
            do id = 1,n_F                                                       ! toroidal modes
                do n_F_loc = -n(id),n(id),max(2*n(id),1)                        ! positive and negative N modes, no duplication at zero
                    id_loc = n_F + (id-1)*n_F_loc/max(n(id),1)
                    do kd = 0,m_F-1                                                 ! poloidal modes
                        ! R
                        XYZ_plot(:,:,1) = XYZ_plot(:,:,1) + &
                            &B(kd+1,id_loc,1,1)*cos(kd*ang_plot(:,:,1)-&
                            &n_F_loc*ang_plot(:,:,2)) + &
                            &B(kd+1,id_loc,2,1)*sin(kd*ang_plot(:,:,1)-&
                            &n_F_loc*ang_plot(:,:,2))
                        ! Z
                        XYZ_plot(:,:,3) = XYZ_plot(:,:,3) + &
                            &B(kd+1,id_loc,1,2)*cos(kd*ang_plot(:,:,1)-&
                            &n_F_loc*ang_plot(:,:,2)) + &
                            &B(kd+1,id_loc,2,2)*sin(kd*ang_plot(:,:,1)-&
                            &n_F_loc*ang_plot(:,:,2))
                    end do
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
    end function calc_eq_2

    ! plots the flux quantities in the solution grid
    !   safety factor q_saf
    !   rotational transform rot
    !   pressure pres
    !   poloidal flux flux_p
    !   toroidal flux flux_t
    integer function flux_q_plot(grid_eq,eq) result(ierr)
        use num_vars, only: rank, no_plots, use_normalization
        use grid_utilities, only: trim_grid
        use MPI_utilities, only: get_ser_var
        use eq_vars, only: pres_0, psi_0
        
        character(*), parameter :: rout_name = 'flux_q_plot'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! normal grid
        type(eq_1_type), intent(in) :: eq                                       ! flux equilibrium variables
        
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
        integer function flux_q_plot_HDF5() result(ierr)
            use num_vars, only: eq_style
            use output_ops, only: plot_HDF5
            use grid_utilities, only: calc_XYZ_grid, extend_grid_F, trim_grid
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
                &X_plot,Y_plot,Z_plot,col=1,description='Flux quantities')
            
            ! deallocate and destroy grid
            deallocate(Y_plot_2D)
            deallocate(X_plot,Y_plot,Z_plot,f_plot)
            call grid_plot%dealloc()
        end function flux_q_plot_HDF5
        
        ! plots the pressure and fluxes in external program
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
    
    ! calculate  R, Z  and Lambda  and derivatives  in VMEC  coordinates at  the
    ! grid  points given  by  the  variables theta_E  and  zeta_E (contained  in
    ! trigon_factors) and at  every normal point. The  derivatives are indicated
    ! by the variable "deriv" which has 3 indices
    ! Note: There is no HELENA equivalent  because for HELENA simulations, R and
    ! Z are not necessary for calculation of the metric coefficients, and L does
    ! not exist.
    integer function calc_RZL_ind(grid,eq,deriv) result(ierr)
        use VMEC_utilities, only: fourier2real
        use VMEC_vars, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, is_asym_V
        use num_utilities, only: check_deriv
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_RZL_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to prepare the trigonometric factors
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(3)                                         ! derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv+1,'calc_RZL')
        CHCKERR('')
        
        ! calculate the variables R,Z and their angular derivative
        ierr = fourier2real(R_V_c(:,grid%i_min:grid%i_max,deriv(1)),&
            &R_V_s(:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)),sym=[.true.,is_asym_V],&
            &deriv=[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(Z_V_c(:,grid%i_min:grid%i_max,deriv(1)),&
            &Z_V_s(:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)),sym=[is_asym_V,.true.],&
            &deriv=[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(L_V_c(:,grid%i_min:grid%i_max,deriv(1)),&
            &L_V_s(:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%L_E(:,:,:,deriv(1),deriv(2),deriv(3)),sym=[is_asym_V,.true.],&
            &deriv=[deriv(2),deriv(3)])
        CHCKERR('')
    end function calc_RZL_ind
    integer function calc_RZL_arr(grid,eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_RZL_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to prepare the trigonometric factors
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_RZL_ind(grid,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_RZL_arr
    
    ! calculate the lower metric elements in the C(ylindrical) coordinate system
    integer function calc_g_C_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_g_C_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_C')
        CHCKERR('')
        
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
    integer function calc_g_C_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_C_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_C_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_C_arr
    
    ! calculate  the metric  coefficients in  the equilibrium  V(MEC) coordinate
    ! system  using  the metric  coefficients  in  the C(ylindrical)  coordinate
    ! system and the transformation matrices
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_g_V_ind(eq,deriv) result(ierr)
        use num_vars, only: max_deriv
        use eq_utilities, only: calc_g
        
        character(*), parameter :: rout_name = 'calc_g_V_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_V')
        CHCKERR('')
        
        ierr = calc_g(eq%g_C,eq%T_VC,eq%g_E,deriv,max_deriv)
        CHCKERR('')
    end function calc_g_V_ind
    integer function calc_g_V_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_V_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_V_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_V_arr
    
    ! calculate the  metric coefficients in the  equilibrium H(ELENA) coordinate
    ! system using the HELENA output
    ! Note: the derivatives are done  using finite differences, which means that
    ! the first  and last points  of the normal grid  have a higher  error. This
    ! could be remedied by choosing the  interval larger. Currently, this is not
    ! done automatically.
    integer function calc_h_H_ind(grid,eq,deriv) result(ierr)
        use num_vars, only: max_deriv, norm_disc_prec_eq
        use HELENA_vars, only: h_H_11, h_H_12, h_H_33
        use num_utilities, only: c
        use grid_utilities, only: setup_deriv_data, apply_disc
        
        character(*), parameter :: rout_name = 'calc_h_H_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: jd, kd                                                       ! counters
        type(disc_type) :: norm_deriv_data                                      ! data for normal derivatives
        type(disc_type) :: ang_deriv_data                                       ! data for angular derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_h_H')
        CHCKERR('')
        
        ! initialize
        eq%h_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        if (deriv(3).ne.0) then                                                 ! axisymmetry: deriv. in tor. coord. is zero
            !eq%h_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (sum(deriv).eq.0) then                                          ! no derivatives
            do jd = 1,grid%n(2)
                eq%h_E(:,jd,:,c([1,1],.true.),0,0,0) = &
                    &h_H_11(:,grid%i_min:grid%i_max)
                eq%h_E(:,jd,:,c([1,2],.true.),0,0,0) = &
                    &h_H_12(:,grid%i_min:grid%i_max)
                eq%h_E(:,jd,:,c([3,3],.true.),0,0,0) = &
                    &h_H_33(:,grid%i_min:grid%i_max)
                eq%h_E(:,jd,:,c([2,2],.true.),0,0,0) = &
                    &1._dp/h_H_11(:,grid%i_min:grid%i_max) * &
                    &(1._dp/(eq%jac_E(:,jd,:,0,0,0)**2*&
                    &h_H_33(:,grid%i_min:grid%i_max)) + &
                    &h_H_12(:,grid%i_min:grid%i_max)**2)
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.0) then                        ! derivative in norm. coord.
            ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,1,&
                &norm_disc_prec_eq+1)                                           ! use extra precision to later calculate mixed derivatives
            CHCKERR('')
            ierr = apply_disc(eq%h_E(:,:,:,:,0,0,0),norm_deriv_data,&
                &eq%h_E(:,:,:,:,1,0,0),3)
            CHCKERR('')
            call norm_deriv_data%dealloc()
        else if (deriv(1).eq.2 .and. deriv(2).eq.0) then                        ! 2nd derivative in norm. coord.
            ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,2,&
                &norm_disc_prec_eq)
            CHCKERR('')
            ierr = apply_disc(eq%h_E(:,:,:,:,0,0,0),norm_deriv_data,&
                &eq%h_E(:,:,:,:,2,0,0),3)
            CHCKERR('')
            call norm_deriv_data%dealloc()
        else if (deriv(1).eq.0 .and. deriv(2).eq.1) then                        ! derivative in pol. coord.
            do kd = 1,grid%loc_n_r
                do jd = 1,grid%n(2)
                    ierr = setup_deriv_data(grid%theta_E(:,jd,kd),&
                        &ang_deriv_data,1,norm_disc_prec_eq+1)                  ! use extra precision to later calculate mixed derivatives
                    CHCKERR('')
                    ierr = apply_disc(eq%h_E(:,jd,kd,:,0,0,0),&
                        &ang_deriv_data,eq%h_E(:,jd,kd,:,0,1,0),1)
                    CHCKERR('')
                end do
            end do
            call ang_deriv_data%dealloc()
        else if (deriv(1).eq.0 .and. deriv(2).eq.2) then                        ! 2nd derivative in pol. coord.
            do kd = 1,grid%loc_n_r
                do jd = 1,grid%n(2)
                    ierr = setup_deriv_data(grid%theta_E(:,jd,kd),&
                        &ang_deriv_data,2,norm_disc_prec_eq)
                    CHCKERR('')
                    ierr = apply_disc(eq%h_E(:,jd,kd,:,0,0,0),&
                        &ang_deriv_data,eq%h_E(:,jd,kd,:,0,2,0),1)
                    CHCKERR('')
                end do
            end do
            call ang_deriv_data%dealloc()
        else if (deriv(1).eq.1 .and. deriv(2).eq.1) then                        ! mixed derivative in norm. and pol. coord.
            ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,1,&
                &norm_disc_prec_eq+1)                                           ! use extra precision to later calculate mixed derivatives
            CHCKERR('')
            ierr = apply_disc(eq%h_E(:,:,:,:,0,1,0),norm_deriv_data,&
                &eq%h_E(:,:,:,:,1,1,0),3) 
            CHCKERR('')
            call norm_deriv_data%dealloc()
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_h_H_ind
    integer function calc_h_H_arr(grid,eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_h_H_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_h_H_ind(grid,eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_h_H_arr

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the  metric coefficients  in  the equilibrium  coordinate  system and  the
    ! trans- formation matrices
    integer function calc_g_F_ind(eq,deriv) result(ierr)
        use num_vars, only: max_deriv
        use eq_utilities, only: calc_g
        
        character(*), parameter :: rout_name = 'calc_g_F_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_F')
        CHCKERR('')
        
        ! calculate g_F
        ierr = calc_g(eq%g_E,eq%T_FE,eq%g_F,deriv,max_deriv)
        CHCKERR('')
    end function calc_g_F_ind
    integer function calc_g_F_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_F_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_g_F_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_F_arr
    
    ! calculate the jacobian in cylindrical coordinates
    integer function calc_jac_C_ind(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_C_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_C')
        CHCKERR('')
        
        eq%jac_C(:,:,:,deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3))
    end function calc_jac_C_ind
    integer function calc_jac_C_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_C_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_C_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_C_arr
    
    ! calculate the jacobian in equilibrium V(MEC) coordinates from 
    !   jac_V = det(T_VC) jac_C
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_V_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_V_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_V')
        CHCKERR('')
        
        ! initialize
        eq%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate determinant
        ierr = add_arr_mult(eq%jac_C,eq%det_T_VC,&
            &eq%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_jac_V_ind
    integer function calc_jac_V_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_V_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_V_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_V_arr
    
    ! calculate the jacobian in HELENA coordinates directly from J = q R^2/F
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_H_ind(grid,eq_1,eq_2,deriv) result(ierr)
        use num_vars, only: norm_disc_prec_eq
        use HELENA_vars, only:  h_H_33, RBphi_H
        use grid_utilities, only: setup_deriv_data, apply_disc
        
        character(*), parameter :: rout_name = 'calc_jac_H_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(disc_type) :: norm_deriv_data                                      ! data for normal derivatives
        type(disc_type) :: ang_deriv_data                                       ! data for angular derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_H')
        CHCKERR('')
        
        ! calculate determinant
        if (deriv(3).ne.0) then                                                 ! axisymmetry: deriv. in tor. coord. is zero
            eq_2%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (sum(deriv).eq.0) then                                          ! no derivatives
            do kd = 1,grid%loc_n_r
                do jd = 1,grid%n(2)
                    eq_2%jac_E(:,jd,kd,0,0,0) = eq_1%q_saf_E(kd,0)/&
                        &(h_H_33(:,grid%i_min-1+kd)*RBphi_H(grid%i_min-1+kd))
                end do
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.0) then                        ! derivative in norm. coord.
            ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,1,&
                &norm_disc_prec_eq+1)                                           ! use extra precision to later calculate mixed derivatives
            CHCKERR('')
            ierr = apply_disc(eq_2%jac_E(:,:,:,0,0,0),norm_deriv_data,&
                &eq_2%jac_E(:,:,:,1,0,0),3)
            CHCKERR('')
            call norm_deriv_data%dealloc()
        else if (deriv(1).eq.2 .and. deriv(2).eq.0) then                        ! 2nd derivative in norm. coord.
            ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,2,&
                &norm_disc_prec_eq)
            CHCKERR('')
            ierr = apply_disc(eq_2%jac_E(:,:,:,0,0,0),norm_deriv_data,&
                &eq_2%jac_E(:,:,:,2,0,0),3)
            CHCKERR('')
            call norm_deriv_data%dealloc()
        else if (deriv(1).eq.0 .and. deriv(2).eq.1) then                        ! derivative in pol. coord.
            do kd = 1,grid%loc_n_r
                do jd = 1,grid%n(2)
                    ierr = setup_deriv_data(grid%theta_E(:,jd,kd),&
                        &ang_deriv_data,1,norm_disc_prec_eq+1)                  ! use extra precision to later calculate mixed derivatives
                    CHCKERR('')
                    ierr = apply_disc(eq_2%jac_E(:,jd,kd,0,0,0),&
                        &ang_deriv_data,eq_2%jac_E(:,jd,kd,0,1,0))
                    CHCKERR('')
                end do
            end do
            call ang_deriv_data%dealloc()
        else if (deriv(1).eq.0 .and. deriv(2).eq.2) then                        ! 2nd derivative in pol. coord.
            do kd = 1,grid%loc_n_r
                do jd = 1,grid%n(2)
                    ierr = setup_deriv_data(grid%theta_E(:,jd,kd),&
                        &ang_deriv_data,2,norm_disc_prec_eq)
                    CHCKERR('')
                    ierr = apply_disc(eq_2%jac_E(:,jd,kd,0,0,0),&
                        &ang_deriv_data,eq_2%jac_E(:,jd,kd,0,2,0))
                    CHCKERR('')
                end do
            end do
            call ang_deriv_data%dealloc()
        else if (deriv(1).eq.1 .and. deriv(2).eq.1) then                        ! mixed derivative in norm. and pol. coord.
            ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,1,&
                &norm_disc_prec_eq+1)                                           ! use extra precision to later calculate mixed derivatives
            CHCKERR('')
            ierr = apply_disc(eq_2%jac_E(:,:,:,0,1,0),norm_deriv_data,&
                &eq_2%jac_E(:,:,:,1,1,0),3) 
            CHCKERR('')
            call norm_deriv_data%dealloc()
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_jac_H_ind
    integer function calc_jac_H_arr(grid,eq_1,eq_2,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_H_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_H_ind(grid,eq_1,eq_2,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_H_arr
    
    ! calculate the jacobian in Flux coordinates from 
    !   jac_F = det(T_FE) jac_E
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_F_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_F_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_F')
        CHCKERR('')
        
        ! initialize
        eq%jac_F(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate determinant
        ierr = add_arr_mult(eq%jac_E,eq%det_T_FE,&
            &eq%jac_F(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_jac_F_ind
    integer function calc_jac_F_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_F_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_F_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_F_arr

    ! calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system
    integer function calc_T_VC_ind(eq,deriv) result(ierr)
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_VC_ind'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VC')
        CHCKERR('')
        
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
        
        ! determinant
        eq%det_T_VC(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = add_arr_mult(eq%R_E(:,:,:,0:,1:,0:),eq%Z_E(:,:,:,1:,0:,0:),&
            &eq%det_T_VC(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ierr = add_arr_mult(-eq%R_E(:,:,:,1:,0:,0:),eq%Z_E(:,:,:,0:,1:,0:),&
            &eq%det_T_VC(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_T_VC_ind
    integer function calc_T_VC_arr(eq,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VC_arr'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VC_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VC_arr

    ! calculate the transformation matrix  between equilibrium V(mec) and F(lux)
    ! oordinate system
    integer function calc_T_VF_ind(grid,eq_1,eq_2,deriv) result(ierr)
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_VF_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium
        integer, intent(in) :: deriv(:)
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:,:)                           ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:,:)                            ! - zeta_F and derivatives
        integer :: dims(3)                                                      ! dimensions
        integer :: c1                                                           ! 2D coordinate in met_type storage convention
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VF')
        CHCKERR('')
        
        ! set up dims
        dims = [grid%n(1),grid%n(2),grid%loc_n_r]
        
        ! initialize
        eq_2%T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! set up theta_s
        allocate(theta_s(dims(1),dims(2),dims(3),&
            &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
        theta_s = 0.0_dp
        ! start from theta_E
        theta_s(:,:,:,0,0,0) = grid%theta_E
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
            zeta_s(:,:,:,0,0,0) = grid%zeta_E
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
    integer function calc_T_VF_arr(grid,eq_1,eq_2,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VF_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VF_ind(grid,eq_1,eq_2,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VF_arr
    
    ! calculate the transformation matrix  between H(ELENA) and F(lux) oordinate
    ! system
    integer function calc_T_HF_ind(grid,eq_1,eq_2,deriv) result(ierr)
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_HF_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium
        integer, intent(in) :: deriv(:)
            
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:,:)                           ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:,:)                            ! - zeta_F and derivatives
        integer :: dims(3)                                                      ! dimensions
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_HF')
        CHCKERR('')
        
        ! set up dims
        dims = [grid%n(1),grid%n(2),grid%loc_n_r]
        
        ! initialize
        eq_2%T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! set up theta_s
        allocate(theta_s(dims(1),dims(2),dims(3),&
            &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
        theta_s = 0.0_dp
        theta_s(:,:,:,0,0,0) = grid%theta_E
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
            zeta_s(:,:,:,0,0,0) = grid%zeta_E
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
    integer function calc_T_HF_arr(grid,eq_1,eq_2,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_HF_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(inout) :: eq_2                                  ! metric equilibrium
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        do id = 1, size(deriv,2)
            ierr = calc_T_HF_ind(grid,eq_1,eq_2,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_HF_arr
    
    ! Calculates derived equilibrium quantities
    ! system [ADD REF]:
    !   - magnetic shear S
    !   - normal curvature kappa_n
    !   - geodesic curvature kappa_g
    !   - parallel current sigma
    subroutine calc_derived_q(grid_eq,eq_1,eq_2)
        use num_utilities, only: c
        use eq_vars, only: vac_perm
        use splines, only: spline3
#if ldebug
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: calc_int
#endif
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(inout), target :: eq_2                          ! metric equilibrium variables
        
        ! local variables
        integer :: kd                                                           ! counter
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
        real(dp), pointer :: D13J(:,:,:) => null()                              ! D_alpha,theta jac
        real(dp), pointer :: D23J(:,:,:) => null()                              ! D_psi,theta jac
        real(dp), pointer :: D23g13(:,:,:) => null()                            ! D_psi,theta g_alpha,theta
        real(dp), pointer :: D13g23(:,:,:) => null()                            ! D_alpha,theta g_psi,theta
        real(dp), allocatable :: D3sigma(:,:,:)                                 ! D_theta sigma
        real(dp), allocatable :: D3sigma_ALT(:,:,:)                             ! alternative D_theta sigma
        real(dp), allocatable :: sigma_ALT(:,:,:)                               ! alternative sigma
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle theta_F or zeta_F
        integer :: istat                                                        ! status
        integer :: jd                                                           ! counter
#endif
        
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
        
        ! Calculate the shear S
        eq_2%S = -(D3h12/h22 - D3h22*h12/h22**2)/J
        
        ! Calculate the normal curvature kappa_n
        do kd = 1,grid_eq%loc_n_r
            eq_2%kappa_n(:,:,kd) = &
                &vac_perm*J(:,:,kd)**2*eq_1%pres_FD(kd,1)/g33(:,:,kd) + &
                &1._dp/(2*h22(:,:,kd)) * ( &
                &h12(:,:,kd) * ( D1g33(:,:,kd)/g33(:,:,kd) - &
                &2*D1J(:,:,kd)/J(:,:,kd) ) + &
                &h22(:,:,kd) * ( D2g33(:,:,kd)/g33(:,:,kd) - &
                &2*D2J(:,:,kd)/J(:,:,kd) ) + &
                &h23(:,:,kd) * ( D3g33(:,:,kd)/g33(:,:,kd) - &
                &2*D3J(:,:,kd)/J(:,:,kd) ) )
        end do
        
        ! Calculate the geodesic curvature kappa_g
        eq_2%kappa_g(:,:,:) = (0.5*D1g33/g33 - D1J/J) - &
            &g13/g33*(0.5*D3g33/g33 - D3J/J)
        
        ! Calculate the parallel current sigma
        eq_2%sigma = 1._dp/vac_perm*&
            &(D1g23 - D2g13 - (g23*D1J - g13*D2J)/J)/J
        do kd = 1,grid_eq%loc_n_r
            eq_2%sigma(:,:,kd) = eq_2%sigma(:,:,kd) - &
                &eq_1%pres_FD(kd,1)*J(:,:,kd)*g13(:,:,kd)/g33(:,:,kd)
        end do 
        
#if ldebug
        ! test whether -2 p' J kappa_g = D3sigma and plot kappa components
        if (debug_calc_derived_q) then
            call writo('Testing whether -2 p'' J kappa_g = D3sigma')
            call lvl_ud(1)
            
            ! allocate variables
            allocate(D3sigma(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            allocate(D3sigma_ALT(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            allocate(sigma_ALT(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
            
            ! point parallel angle
            if (use_pol_flux_F) then
                ang_par_F => grid_eq%theta_F
            else
                ang_par_F => grid_eq%zeta_F
            end if
            
            ! get derived sigma
            do kd = 1,grid_eq%loc_n_r
                do jd = 1,grid_eq%n(2)
                    istat = spline3(ang_par_F(:,jd,kd),eq_2%sigma(:,jd,kd),&
                        &ang_par_F(:,jd,kd),dynew=D3sigma(:,jd,kd))
                    CHCKSTT
                end do
            end do
            !!D3sigma_ALT = -1._dp/vac_perm*&
                !!&(D1g23 - D2g13 - (g23*D1J - g13*D2J)/J)/J**2*D3J + &
                !!&1._dp/vac_perm*&
                !!&(D13g23 - D23g13 - (D3g23*D1J + g23*D13J - &
                !!&D3g13*D2J - g13*D23J)/J + (g23*D1J - g13*D2J)/J**2*D3J)/J
            !!do kd = 1,grid_eq%loc_n_r
                !!D3sigma_ALT(:,:,kd) = D3sigma_ALT(:,:,kd) - &
                    !!&eq_1%pres_FD(kd,1)*(D3J(:,:,kd)*g13(:,:,kd)/g33(:,:,kd) + &
                    !!&J(:,:,kd)*D3g13(:,:,kd)/g33(:,:,kd) - &
                    !!&J(:,:,kd)*g13(:,:,kd)/g33(:,:,kd)**2*D3g33(:,:,kd))
            !!end do 
            
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
                    istat = calc_int(D3sigma_ALT(:,jd,kd),ang_par_F(:,jd,kd),&
                        &sigma_ALT(:,jd,kd))
                    CHCKSTT
                end do
                sigma_ALT(:,:,kd) = eq_2%sigma(1,1,kd) + sigma_ALT(:,:,kd)
            end do
            
            ! plot output
            call plot_diff_HDF5(D3sigma,D3sigma_ALT,'TEST_D3sigma',&
                &grid_eq%n,[0,0,grid_eq%i_min-1],&
                &description='To test whether -2 p'' J kappa_g = D3sigma',&
                &output_message=.true.)
            call plot_diff_HDF5(eq_2%sigma,sigma_ALT,'TEST_sigma',grid_eq%n,&
                &[0,0,grid_eq%i_min-1],description='To test whether &
                &int(-2 p'' J kappa_g) = sigma',output_message=.true.)
            
            ! plot shear
            call plot_HDF5('shear','TEST_shear',eq_2%S,&
                &tot_dim=grid_eq%n,loc_offset=[0,0,grid_eq%i_min-1])
            
            ! plot kappa_n
            call plot_HDF5('kappa_n','TEST_kappa_n',eq_2%kappa_n,&
                &tot_dim=grid_eq%n,loc_offset=[0,0,grid_eq%i_min-1])
            
            ! plot kappa_g
            call plot_HDF5('kappa_g','TEST_kappa_g',eq_2%kappa_g,&
                &tot_dim=grid_eq%n,loc_offset=[0,0,grid_eq%i_min-1])
            
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
    end subroutine calc_derived_q
    
    ! Sets up normalization constants.
    subroutine calc_normalization_const()
        use num_vars, only: rank, eq_style, mu_0_original, use_normalization, &
            &rich_restart_lvl
        use eq_vars, only: T_0, B_0, pres_0, psi_0, R_0, rho_0
        
        ! local variables
        integer :: nr_overriden_const                                           ! nr. of user-overriden constants, to print warning if > 0
        
        if (use_normalization) then
            ! user output
            call writo('Calculating the normalization constants')
            call lvl_ud(1)
            
            ! initialize nr_overriden_const
            nr_overriden_const = 0
            
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
            
            ! print warning if user-overriden
            if (nr_overriden_const.gt.0) &
                &call writo(trim(i2str(nr_overriden_const))//&
                &' constants were overriden by user. Consistency is NOT &
                &checked!',warning=.true.)
            
            ! print constants
            call print_normalization_const(R_0,rho_0,B_0,pres_0,psi_0,&
                &mu_0_original,T_0)
            
            ! user output
            call lvl_ud(-1)
            call writo('Normalization constants calculated')
        else if (rich_restart_lvl.eq.1) then                                    ! only for first Richardson leevel
            ! user output
            call writo('Normalization not used')
        end if
    contains 
        ! VMEC version
        ! Normalization depends on style:
        !   1: MISHKA
        !       R_0:    major radius (= average magnetic axis)
        !       B_0:    B on magnetic axis (theta = zeta = 0)
        !       pres_0: reference pressure (= B_0^2/mu_0)
        !       psi_0:  reference flux (= R_0^2 B_0)
        !       rho_0:  reference mass density
        !   2: COBRA
        !       R_0:    major radius (= average geometric axis)
        !       pres_0: pressure on magnetic axis
        !       B_0:    reference magnetic field (= sqrt(2pres_0mu_0 / beta))
        !       psi_0:  reference flux (= R_0^2 B_0 / aspr^2)
        !       rho_0:  reference mass density
        !   where aspr (aspect ratio) and beta are given by VMEC.
        ! Note that  rho_0 is  not given  through by  the equilibrium  codes and
        ! should be user-supplied
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
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                    
                    ! set B_0 from magnetic field on axis
                    if (B_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        B_0 = B_0_V
                    else
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                    
                    ! set reference pres_0 from B_0
                    if (pres_0.ge.huge(1._dp)) then                             ! user did not provide a value
                        pres_0 = B_0**2/mu_0_original
                    else
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                    
                    ! set reference flux from R_0 and B_0
                    if (psi_0.ge.huge(1._dp)) then                              ! user did not provide a value
                        psi_0 = R_0**2 * B_0
                    else
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                case (2)                                                        ! COBRA
                    ! user output
                    call writo('Using COBRA normalization')
                    
                    ! set the major radius as the average geometric axis
                    if (R_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        R_0 = 0.5_dp*(rmin_surf+rmax_surf)
                    else
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                    
                    ! set pres_0 from pressure on axis
                    if (pres_0.ge.huge(1._dp)) then                             ! user did not provide a value
                        pres_0 = pres_V(1,0)
                    else
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                    
                    ! set the reference value for B_0 from pres_0 and beta
                    if (B_0.ge.huge(1._dp)) then                                ! user did not provide a value
                        !B_0 = sqrt(2._dp*pres_0*mu_0_original/beta_V)          ! exact COBRA
                        B_0 = sqrt(pres_0*mu_0_original)                        ! pure modified COBRA
                    else
                        nr_overriden_const = nr_overriden_const + 1
                    end if
                    
                    ! set reference flux from R_0, B_0 and aspr
                    if (psi_0.ge.huge(1._dp)) then                              ! user did not provide a value
                        !psi_0 = B_0 * (R_0/aspr_V)**2                          ! exact COBRA
                        psi_0 = B_0 * R_0**2                                    ! pure modified COBRA
                    else
                        nr_overriden_const = nr_overriden_const + 1
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
                nr_overriden_const = nr_overriden_const + 1
            end if
        end subroutine calc_normalization_const_VMEC
        
        ! HELENA version
        ! The MISHKA normalization is taken by default, see "read_VMEC" for more
        ! information.
        subroutine calc_normalization_const_HEL
            ! user output
            call writo('Using MISHKA normalization')
            
            if (R_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                R_0 = 1._dp
            end if
            if (B_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                B_0 = 1._dp
            end if
            if (pres_0.ge.huge(1._dp)) then                                     ! user did not provide a value
                pres_0 = B_0**2/mu_0_original
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            if (psi_0.ge.huge(1._dp)) then                                      ! user did not provide a value
                psi_0 = R_0**2 * B_0
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            if (T_0.ge.huge(1._dp)) T_0 = sqrt(mu_0_original*rho_0)*R_0/B_0     ! only if user did not provide a value
        end subroutine calc_normalization_const_HEL
        
        ! prints the Normalization factors
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
    end subroutine calc_normalization_const
    
    ! Normalize input quantities.
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
    
    ! Plots the  magnetic fields. If  multiple equilibrium parallel  jobs, every
    ! job does its piece, and the results are joined automatically by plot_HDF5.
    ! The outputs are given in contra- and covariant components and magnitude in
    ! multiple coordinate systems, as indicated in "calc_vec_comp".
    ! The starting point is the fact that the magnetic field is given by
    !   B = e_theta/J
    ! in F coordinates. The F covariant components are therefore given by
    !   B_i = g_i3/J
    ! and the only non-vanishing contravariant component is
    !   B^3 = 1/J.
    ! These are then all be transformed to the other coordinate systems.
    ! Note that vector plots for different  Richardson levels can be combined to
    ! show the total grid by just plotting them all individually.
    ! Note: The metric factors and transformation matrices have to be allocated.
    integer function B_plot(grid_eq,eq_1,eq_2,rich_lvl,plot_fluxes,XYZ) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, use_normalization
        use num_utilities, only: c
        use eq_vars, only: B_0
        
        character(*), parameter :: rout_name = 'B_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        integer, intent(in), optional :: rich_lvl                               ! Richardson level
        logical, intent(in), optional :: plot_fluxes                            ! plot the fluxes
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          ! X, Y and Z of grid
        
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
    
    ! Plots the current.  If multiple equilibrium parallel jobs,  every job does
    ! its  piece, and  the results  are joined  automatically by  plot_HDF5. The
    ! outputs are  given in  contra- and covariant  components and  magnitude in
    ! multiple coordinate systems, as indicated in "calc_vec_comp".
    ! The starting point is the pressure balance
    !   nabla p = J x B,
    ! which, using B = e_theta/J, reduces to
    !   J^alpha = -p'.
    ! Furthermore, the current has to lie in the magnetic flux surfaces:
    !   J^psi = 0.
    ! Finally,  the parallel  current sigma  gives  an expression  for the  last
    ! contravariant component:
    !   J^theta = sigma/J + p' B_alpha/B_theta.
    ! From these, the contravariant components can be calculated as
    !   J_i = J^alpha g_alpha,i + J^theta g_theta,i.
    ! These are then all be transformed to the other coordinate systems.
    ! Note that vector plots for different  Richardson levels can be combined to
    ! show the total grid by just plotting them all individually.
    ! Note: The metric factors and transformation matrices have to be allocated.
    integer function J_plot(grid_eq,eq_1,eq_2,rich_lvl,plot_fluxes,XYZ) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, use_normalization
        use num_utilities, only: c
        use eq_vars, only: B_0, R_0, pres_0
        
        character(*), parameter :: rout_name = 'J_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        integer, intent(in), optional :: rich_lvl                               ! Richardson level
        logical, intent(in), optional :: plot_fluxes                            ! plot the fluxes
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          ! X, Y and Z of grid
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: J_com(:,:,:,:,:)                               ! covariant and contravariant components of J (dim1,dim2,dim3,3,2)
        real(dp), allocatable :: J_mag(:,:,:)                                   ! magnitude of J (dim1,dim2,dim3)
        character(len=10) :: base_name                                          ! base name
        real(dp), allocatable, save :: J_flux_tor(:,:), J_flux_pol(:,:)         ! fluxes
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
    end function J_plot
    
    ! Plots the curvature. If multiple equilibrium parallel jobs, every job does
    ! its  piece, and  the results  are joined  automatically by  plot_HDF5. The
    ! outputs are  given in  contra- and covariant  components and  magnitude in
    ! multiple coordinate systems, as indicated in "calc_vec_comp".
    ! The starting point is the curvature, given by
    !   kappa = kappa_n nabla psi / |nabla psi|^2 + kappa_g nabla psi x B / B^2
    ! which can  be used to find  the covariant and contravariant  components in
    ! Flux coordinates. These are then  transformed to Cartesian coordinates and
    ! plotted.
    ! Note that vector plots for different  Richardson levels can be combined to
    ! show the total grid by just plotting them all individually.
    ! Note: The metric factors and transformation matrices have to be allocated.
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
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! metric equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        integer, intent(in), optional :: rich_lvl                               ! Richardson level
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          ! X, Y and Z of grid
        
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
                    &description='center of curvature')
                
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
                    &description='center of curvature')
                
                ! clean up
                call grid_trim%dealloc()
                
                call lvl_ud(-1)
            end if
        end if
#endif
    end function kappa_plot
    
    ! Plots the  magnetic fields. If  multiple equilibrium parallel  jobs, every
    ! job does its piece, and the results are joined automatically by plot_HDF5.
    ! The outputs are given in contra- and covariant components and magnitude in
    ! multiple coordinate systems, as indicated in "calc_vec_comp".
    ! The starting point is the fact that the magnetic field is given by
    !   B = e_theta/J
    ! in F coordinates. The F covariant components are therefore given by
    !   B_i = g_i3/J
    ! and the only non-vanishing contravariant component is
    !   B^3 = 1/J.
    ! These are then all be transformed to the other coordinate systems.
    ! Note that vector plots for different  Richardson levels can be combined to
    ! show the total grid by just plotting them all individually.
    ! Note: The metric factors and transformation matrices have to be allocated.
    integer function delta_r_plot(grid_eq,eq_1,eq_2,XYZ,rich_lvl) &
        &result(ierr)
        
        use grid_utilities, only: calc_vec_comp, calc_XYZ_grid, trim_grid
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_style, norm_disc_prec_eq, eq_job_nr, &
            &use_normalization, rank, use_pol_flux_F, ex_plot_style, n_procs, &
            &prop_B_tor_i, min_theta_plot
        use eq_vars, only: B_0, R_0, psi_0, max_flux_F
        use num_utilities, only: c, calc_int
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'delta_r_plot'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        real(dp), intent(in) :: XYZ(:,:,:,:)                                    ! X, Y and Z of grid
        integer, intent(in), optional :: rich_lvl                               ! Richardson level
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: r_lim_2D(2)                                                  ! limits of r to plot in 2-D
        real(dp) :: flux_code                                                   ! positive for poloidal flux, negative for toroidal
        real(dp), allocatable :: XYZ_loc(:,:,:,:)                               ! X, Y and Z of surface
        real(dp), allocatable :: B_com(:,:,:,:,:)                               ! covariant and contravariant components of B (dim1,dim2,dim3,3,2)
        real(dp), allocatable :: delta_B(:,:,:,:)                               ! delta_B/B
        real(dp), allocatable :: prop_B_tor_tot(:,:,:)                          ! delta_r / delta_B/B in total grid
        real(dp), allocatable :: var_tot_loc(:)                                 ! auxilliary variable
        real(dp), allocatable :: x_2D(:,:)                                      ! x for 2-D plotting
        real(dp), allocatable :: delta_r(:,:,:)                                 ! normal displacement
        real(dp), allocatable :: rad(:,:,:)                                     ! minor radius
        real(dp), allocatable :: prev_rad(:)                                    ! rad of previous ranks
        real(dp), allocatable :: arclen(:,:,:)                                  ! length along the arc for the proportionality factor
        real(dp), allocatable :: arclen_tot(:,:,:)                              ! arclen in total grid
        real(dp), allocatable :: u_norm_com(:,:,:,:,:)                          ! covariant and contravariant components of u_norm (dim1,dim2,dim3,3,2)
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle
        logical :: new_file_found                                               ! name for new file found
        character(len=8) :: flux_name                                           ! either poloidal or toroidal
        character(len=25) :: base_name                                          ! base name
        character(len=25) :: plot_names(2)                                      ! plot names
        character(len=max_str_ln) :: prop_B_tor_file_name                       ! name of B_tor proportionality file
        character(len=max_str_ln) :: format_head                                ! header format
        character(len=max_str_ln) :: format_val                                 ! format
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
        allocate(XYZ_loc(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3))
        ierr = calc_XYZ_grid(grid_eq,grid_eq,XYZ_loc(:,:,:,1),XYZ_loc(:,:,:,2),&
            &XYZ_loc(:,:,:,3))
        CHCKERR('')
        
        ! set up components of u_norm and calculate minor radius
        allocate(u_norm_com(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3,2))
        allocate(rad(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
        allocate(prev_rad(n_procs))
        u_norm_com = 0._dp
        rad = 0._dp
        u_norm_com(:,:,:,2,1) = &
            &(eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0))**(-0.5_dp)                ! normal unit vector nabla psi / |nabla psi|
        do id = 1,3
            u_norm_com(:,:,:,id,2) = eq_2%h_FD(:,:,:,c([2,id],.true.),0,0,0)*&
                &u_norm_com(:,:,:,2,1)
        end do
        do jd = 1,grid_trim%n(2)
            do id = 1,grid_trim%n(1)
                ierr = calc_int((eq_2%h_FD(id,jd,norm_id(1):norm_id(2),&
                    &c([2,2],.true.),0,0,0))**(-0.5),grid_trim%loc_r_F,&
                    &rad(id,jd,:))
                CHCKERR('')
                if (use_normalization) rad(id,jd,:) = &
                    &rad(id,jd,:) * psi_0/(R_0*B_0)
                ierr = get_ser_var([rad(id,jd,grid_trim%loc_n_r)],prev_rad,&
                    &scatter=.true.)
                CHCKERR('')
                rad(id,jd,:) = rad(id,jd,:) + sum(prev_rad(1:rank))
            end do
        end do
        call plot_HDF5('rad','rad',rad,&
            &tot_dim=[grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &X=XYZ(:,:,norm_id(1):norm_id(2),1),&
            &Y=XYZ(:,:,norm_id(1):norm_id(2),2),&
            &Z=XYZ(:,:,norm_id(1):norm_id(2),3),&
            &cont_plot=eq_job_nr.gt.1,&
            &description='minor radius')
        
        ! transform coordinates, including the flux
        ierr = calc_vec_comp(grid_eq,grid_eq,eq_1,eq_2,u_norm_com,&
            &norm_disc_prec_eq,XYZ=XYZ,compare_tor_pos=.false.)
        CHCKERR('')
        
        ! dot  position difference  with Cartesian components  of u_norm  to get
        ! delta_r
        allocate(delta_r(grid_eq%n(1),1,grid_eq%loc_n_r))
        delta_r = 0._dp
        do id = 1,3
            delta_r(:,1,:) = delta_r(:,1,:) + 0.5_dp * &
                &(u_norm_com(:,2,:,id,1)+u_norm_com(:,1,:,id,1)) * &
                &(XYZ_loc(:,2,:,id)-XYZ_loc(:,1,:,id))
        end do
        
        ! set plot variables for delta_r
        base_name = 'delta_r'
        plot_names(1) = 'delta_r'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        
        ! plot
        call plot_HDF5(trim(plot_names(1)),trim(base_name),&
            &delta_r(:,:,norm_id(1):norm_id(2)),&
            &tot_dim=[grid_trim%n(1),1,grid_trim%n(3)],&
            &loc_offset=[0,0,grid_trim%i_min-1],&
            &X=0.5_dp*(XYZ(:,1:1,norm_id(1):norm_id(2),1)+&
            &XYZ(:,2:2,norm_id(1):norm_id(2),1)),&
            &Y=0.5_dp*(XYZ(:,1:1,norm_id(1):norm_id(2),2)+&
            &XYZ(:,2:2,norm_id(1):norm_id(2),2)),&
            &Z=0.5_dp*(XYZ(:,1:1,norm_id(1):norm_id(2),3)+&
            &XYZ(:,2:2,norm_id(1):norm_id(2),3)),&
            &cont_plot=eq_job_nr.gt.1,&
            &description='plasma position displacement')
        
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
        
        ! set plot variables for prop_B_tor
        base_name = 'prop_B_tor_from_eq'
        plot_names(1) = trim(base_name)//'_sub'
        plot_names(2) = trim(base_name)//'_sup'
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) base_name = trim(base_name)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        allocate(delta_B(grid_eq%n(1),1,grid_eq%loc_n_r,2))
        delta_B = 2*(B_com(:,2:2,:,2,:)-B_com(:,1:1,:,2,:)) / &
            &(B_com(:,2:2,:,2,:)+B_com(:,1:1,:,2,:))
        
        call lvl_ud(-1)
        
        call writo('Output ripple proportionality factor')
        call lvl_ud(1)
        
        ! plot with HDF5
        call plot_HDF5(plot_names,trim(base_name),&
            &reshape([delta_r(:,:,norm_id(1):norm_id(2)),&
            &delta_r(:,:,norm_id(1):norm_id(2))],&
            &[grid_trim%n(1),1,grid_trim%loc_n_r,2])/&
            &delta_B(:,:,norm_id(1):norm_id(2),:),&
            &tot_dim=[grid_trim%n(1),1,grid_trim%n(3),2],&
            &loc_offset=[0,0,grid_trim%i_min-1,0],&
            &X=0.5_dp*(XYZ(:,1:1,norm_id(1):norm_id(2),1:1)+&
            &XYZ(:,2:2,norm_id(1):norm_id(2),1:1)),&
            &Y=0.5_dp*(XYZ(:,1:1,norm_id(1):norm_id(2),2:2)+&
            &XYZ(:,2:2,norm_id(1):norm_id(2),2:2)),&
            &Z=0.5_dp*(XYZ(:,1:1,norm_id(1):norm_id(2),3:3)+&
            &XYZ(:,2:2,norm_id(1):norm_id(2),3:3)),&
            &cont_plot=eq_job_nr.gt.1,description='delta_r divided by delta_B')
        
        ! plot in 2-D
        r_lim_2D = [1,grid_trim%n(3)]
        if (ex_plot_style.eq.2) r_lim_2D(1) = &
            &max(r_lim_2D(1),r_lim_2D(2)-127+1)                                 ! Bokeh can only handle 255/2 input arguments
        if (rank.eq.0) then
            allocate(prop_B_tor_tot(grid_trim%n(1),grid_trim%n(3),2))
            allocate(x_2D(grid_trim%n(1),grid_trim%n(3)))
        end if
        if (use_pol_flux_F) then
            ang_par_F => grid_trim%theta_F
        else
            ang_par_F => grid_trim%zeta_F
        end if
        do id = 1,grid_trim%n(1)
            do jd = 1,2
                ierr = get_ser_var(delta_r(id,1,norm_id(1):norm_id(2))/&
                    &delta_B(id,1,norm_id(1):norm_id(2),jd),var_tot_loc)
                CHCKERR('')
                if (rank.eq.0) prop_B_tor_tot(id,:,jd) = var_tot_loc
            end do
            ierr = get_ser_var(ang_par_F(id,1,:),var_tot_loc)
            CHCKERR('')
            if (rank.eq.0) x_2D(id,:) = var_tot_loc/pi
        end do
        if (rank.eq.0) then
            do jd = 1,2
                call print_ex_2D([trim(plot_names(jd))],trim(base_name),&
                    &prop_B_tor_tot(:,r_lim_2D(1):r_lim_2D(2),jd),&
                    &x=x_2D(:,r_lim_2D(1):r_lim_2D(2)),draw=.false.)
                call draw_ex([trim(plot_names(jd))],trim(base_name),&
                    &r_lim_2D(2)-r_lim_2D(1)+1,1,.false.)
            end do
        end if
        
        call lvl_ud(-1)
        
        call writo('Output to file as function of fractional arc length')
        call lvl_ud(1)
        
        ! calculate arclength: int dtheta |e_theta|
        allocate(arclen(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        do kd = 1,grid_eq%loc_n_r
            do jd = 1,grid_eq%n(2)
                if (use_pol_flux_F) then
                    ! e_arc = e_theta - q e_alpha
                    ierr = calc_int((&
                        &eq_2%g_FD(:,jd,kd,c([3,3],.true.),0,0,0) - &
                        &eq_2%g_FD(:,jd,kd,c([1,3],.true.),0,0,0) * &
                        &2*eq_1%q_saf_FD(kd,0) + &
                        &eq_2%g_FD(:,jd,kd,c([1,1],.true.),0,0,0) * &
                        &eq_1%q_saf_FD(kd,0)**2)**(0.5),&
                        &grid_eq%theta_F(:,jd,kd),arclen(:,jd,kd))
                    CHCKERR('')
                else
                    ! e_arc = -e_alpha
                    ierr = calc_int(&
                        &(-eq_2%g_FD(:,jd,kd,c([1,1],.true.),0,0,0))**(0.5),&
                        &grid_eq%theta_F(:,jd,kd),arclen(:,jd,kd))
                    CHCKERR('')
                end if
                if (use_normalization) arclen(:,jd,kd) = arclen(:,jd,kd) * R_0  ! not really necessary as we take fractional arc length
                arclen(:,jd,kd) = arclen(:,jd,kd) / arclen(grid_trim%n(1),jd,kd)
            end do
        end do
        if (rank.eq.0) then
            allocate(arclen_tot(grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)))
        end if
        do id = 1,grid_trim%n(1)
            do jd = 1,grid_trim%n(2)
                ierr = get_ser_var(arclen(id,jd,norm_id(1):norm_id(2)),&
                    &var_tot_loc)
                CHCKERR('')
                if (rank.eq.0) arclen_tot(id,jd,:) = var_tot_loc
            end do
        end do
        
        ! only master outputs to file
        if (rank.eq.0) then
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
                    prop_B_tor_file_name = trim(prop_B_tor_file_name)//'_'//&
                        &trim(i2str(kd))//'.dat'
                    new_file_found = .true.
                    open(prop_B_tor_i,FILE=trim(prop_B_tor_file_name),&
                        &IOSTAT=ierr,STATUS='new')
                    CHCKERR('Failed to open file')
                end if
            end do
            format_head = '("# ",'//trim(i2str(size(prop_B_tor_tot,1)+1))//&
                &'(A23," "))'
            format_val = '("  ",'//trim(i2str(size(prop_B_tor_tot,1)+1))//&
                &'(ES23.16," "))'
            
            ! user output
            call writo('Save toroidal field proportionality factor in file "'//&
                &trim(prop_B_tor_file_name)//'"')
            
            ! write to output file
            if (use_normalization) prop_B_tor_tot = prop_B_tor_tot / R_0
            if (use_pol_flux_F) then
                flux_name = 'poloidal'
                flux_code = 1._dp
            else
                flux_name = 'toroidal'
                flux_code = -1._dp
            end if
            write(prop_B_tor_i,'(A)') '# Proportionality factor to be &
                &multiplied by deltaB/B in order to get delta_r/R_0 necessary'
            write(prop_B_tor_i,trim(format_head)) 'n_pol', 'n_r', &
                &'ang_pol_start'
            write(prop_B_tor_i,trim(format_val)) 1._dp*grid_trim%n(1), &
                &1._dp*grid_trim%n(3), min_theta_plot*pi, flux_code
            write(prop_B_tor_i,trim(format_head)) trim(flux_name)//' s', &
                &'poloidal angle'
            write(prop_B_tor_i,trim(format_head)) trim(flux_name)//' s', &
                &'proportionality factor'
            do kd = 1,grid_trim%n(3)
                write(prop_B_tor_i,trim(format_val)) &
                    &grid_trim%r_F(kd)*2*pi/max_flux_F, &
                    &0.5*sum(arclen_tot(:,:,kd),2)
                write(prop_B_tor_i,trim(format_val)) &
                    &grid_trim%r_F(kd)*2*pi/max_flux_F, &
                    &0.5*sum(prop_B_tor_tot(:,kd,:),2)
            end do
            
            ! close
            close(prop_B_tor_i)
        end if
        
        call lvl_ud(-1)
        
        ! clean up
        call grid_trim%dealloc()
        nullify(ang_par_F)
    end function delta_r_plot
    
    ! Print equilibrium quantities to an output file:
    !   - flux:     pres_FD, q_saf_FD, rot_t_FD, flux_p_FD, flux_t_FD, rho, S,
    !               kappa_n, kappa_g, sigma
    !   - metric:   g_FD, h_FD, jac_FD
    ! If "rich_lvl" is  provided, "_R_rich_lvl" is appended to the  data name if
    ! it is > 0 (only for eq_2).
    ! Optionally, for eq_2, it can be  specified that this is a divided parallel
    ! grid, corresponding to the variable "eq_jobs_lims" with index "eq_job_nr".
    ! In this  case, the  total grid size  is adjusted to  the one  specified by
    ! "eq_jobs_lims" and the grid is written as a subset.
    ! Note: The equilibrium quantities are outputted in Flux coordinates.
    ! Note: The  metric equilibrium  quantities can be  deallocated on  the fly,
    ! which is useful if this routine is  followed by a deallocation any way, so
    ! that memory usage does not almost double.
    ! Note: print_output_eq_2 is only used by HELENA  now, as for VMEC it is too
    ! slow  since there  are often  multiple  VMEC equilibrium  jobs, while  for
    ! HELENA this is explicitely forbidden.
    integer function print_output_eq_1(grid_eq,eq,data_name) result(ierr)       ! flux version
        use num_vars, only: PB3D_name
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_eq_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(eq_1_type), intent(in) :: eq                                       ! flux equilibrium variables
        character(len=*), intent(in) :: data_name                               ! name under which to store
        
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
    integer function print_output_eq_2(grid_eq,eq,data_name,rich_lvl,par_div,&
        &dealloc_vars) result(ierr)                                             ! metric version
        use num_vars, only: PB3D_name, eq_jobs_lims, eq_job_nr
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_eq_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium variables
        character(len=*), intent(in) :: data_name                               ! name under which to store
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to print
        logical, intent(in), optional :: par_div                                ! is a parallely divided grid
        logical, intent(in), optional :: dealloc_vars                           ! deallocate variables on the fly after writing
        
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
    
    ! Redistribute  the equilibrium variables,  but only the Flux  variables are
    ! saved. See "redistribute_output_grid" for more information.
    integer function redistribute_output_eq_1(grid,grid_out,eq,eq_out) &
        &result(ierr)                                                           ! flux version
        use MPI_utilities, only: redistribute_var
        
        character(*), parameter :: rout_name = 'redistribute_output_eq_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid variables
        type(grid_type), intent(in) :: grid_out                                 ! redistributed equilibrium grid variables
        type(eq_1_type), intent(inout) :: eq                                    ! flux equilibrium variables
        type(eq_1_type), intent(inout) :: eq_out                                ! flux equilibrium variables in redistributed grid
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Redistribute flux equilibrium variables')
        call lvl_ud(1)
        
        ! create redistributed flux equilibrium variables
        call eq_out%init(grid_out,setup_E=.false.,setup_F=.true.)
        
        ! for all derivatives
        do id = 0,max_deriv
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
    integer function redistribute_output_eq_2(grid,grid_out,eq,eq_out) &
        &result(ierr)                                                           ! metric version
        use MPI_utilities, only: redistribute_var
        
        character(*), parameter :: rout_name = 'redistribute_output_eq_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid variables
        type(grid_type), intent(in) :: grid_out                                 ! redistributed equilibrium grid variables
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium variables
        type(eq_2_type), intent(inout) :: eq_out                                ! metric equilibrium variables in redistributed grid
        
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

    ! Divides the equilibrium jobs.
    ! For  PB3D, the  entire parallel  range has  to be  calculated, but  due to
    ! memory limits  this has to be  split up in  pieces. Every piece has  to be
    ! able to contain the equilibrium variables (see note below), as well as the
    ! vectorial perturbation variables. These  are later combined into tensorial
    ! variables and integrated.
    ! The equilibrium variables have to be  operated on to calculate them, which
    ! translates to a scale factor "mem_scale_fac". However, in the perturbation
    ! phase, when they are just used, this scale factor is not needed.
    ! In its  most extreme form, the  division in equilibrium jobs  would be the
    ! individual  calculation  on  a  fundamental integration  integral  of  the
    ! parallel points:
    !   - for magn_int_style = 1 (trapezoidal), this is 1 point,
    !   - for magn_int_style = 2 (Simpson 3/8), this is 3 points.
    ! For  HELENA,  the  parallel  derivatives are  calculated  discretely,  the
    ! equilibrium and  vectorial perturbation  variables are tabulated  first in
    ! this  HELENA grid.  This happens  in the  first Richardson  level. In  all
    ! Richardson  levels, afterwards,  these variables  are interpolated  in the
    ! angular directions. In  this case, therefore, there can be  no division of
    ! this HELENA output interval for the first Richardson level.
    ! This  procedure does  the job  of dividing  the grids  setting the  global
    ! variables 'eq_jobs_lims'.
    ! The integration of the tensorial perturbation variables is adjusted:
    !   - If  the first job  of the parallel jobs  and not the  first Richardson
    !   level: add half  of the integrated tensorial  perturbation quantities of
    !   the previous level.
    !   -  If  not the  first  job  of the  parallel  jobs,  add the  integrated
    !   tensorial perturbation quantities to those of the previous parallel job,
    !   same Richardson level.
    ! In fact,  the equilibrium  jobs have  much in  common with  the Richardson
    ! levels,  as is  attested  by the  existence of  the  routines "do_eq"  and
    ! "eq_info", which are equivalent to "do_rich" and "rich_info".
    ! In POST, finally,  the situation is slightly different for  HELENA, as all
    ! the requested variables have to fit, including the interpolated variables,
    ! as they are stored whereas in PB3D  they are not. The parallel range to be
    ! taken is then the  one of the output grid, including a  base range for the
    ! variables tabulated on  the HELENA grid. Also, for  extended output grids,
    ! the size of the grid in the  secondary angle has to be included in n_par_X
    ! (i.e. toroidal  when poloidal flux  is used and vice  versa). Furthermore,
    ! multiple equilibrium jobs are allowed.
    ! To this end,  optionally, a base number can be  provided for n_par_X, that
    ! is always added to the number of points in the divided n_par_X.
    ! Note: For PB3D,  only the variables g_FD, h_FD and  jac_FD are counted, as
    ! the  equilibrium variables  and  the transformation  matrices are  deleted
    ! after use. Also, S, sigma, kappa_n and kappa_g can be neglected as they do
    ! not  contain  derivatives   and  are  therefore  much   smaller.  in  both
    ! "calc_memory" routines,  however, a 50%  safety factor is used  to account
    ! for this somewhat.
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
        integer, intent(in) :: n_par_X                                          ! number of parallel points to be divided
        integer, intent(in) :: arr_size(2)                                      ! array size (using loc_n_r) for eq_2 and X_1 variables
        integer, intent(inout) :: n_div                                         ! final number of divisions
        integer, intent(in), optional :: n_div_max                              ! maximum n_div
        integer, intent(in), optional :: n_par_X_base                           ! base n_par_X, undivisible
        character(len=*), intent(in), optional :: range_name                    ! name of range
        
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
    
    ! Calculate eq_jobs_lims.
    ! Take into account that every job has to start and end at the start and end
    ! of a fundamental integration interval, as discussed in above:
    !   - for magn_int_style = 1 (trapezoidal), this is 1 point,
    !   - for magn_int_style = 2 (Simpson 3/8), this is 3 points.
    ! for POST, there are  no Richardson levels, and there has  to be overlap of
    ! one always, in order to have  correct composite integrals of the different
    ! regions.
    integer function calc_eq_jobs_lims(n_par_X,n_div) result(ierr)
        use num_vars, only: prog_style, eq_jobs_lims, eq_job_nr
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'calc_eq_jobs_lims'
        
        ! input / output
        integer, intent(in) :: n_par_X                                          ! number of parallel points in this Richardson level
        integer, intent(in) :: n_div                                            ! nr. of divisions of parallel ranges
        
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
    ! See if T_EF it complies with the theory of [ADD REF]
    integer function test_T_EF(grid_eq,eq_1,eq_2) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style
        use grid_utilities, only: trim_grid
        use num_utilities, only: c
        use output_ops, only: plot_diff_HDF5
        
        character(*), parameter :: rout_name = 'test_T_EF'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium
        
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
    
    ! Tests whether D1 D2 h_H is calculated correctly
    integer function test_D12h_H(grid_eq,eq) result(ierr)
        use grid_utilities, only: trim_grid, setup_deriv_data, apply_disc
        use num_utilities, only: c
        use num_vars, only: norm_disc_prec_eq
        
        character(*), parameter :: rout_name = 'test_D12h_H'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       ! metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        type(disc_type) :: ang_deriv_data                                       ! data for angular derivative
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
        
        ! calculate D1 D2 h_H alternatively
        do kd = norm_id(1),norm_id(2)
            do jd = 1,grid_trim%n(2)
                ierr = setup_deriv_data(grid_eq%theta_E(:,jd,kd),&
                    &ang_deriv_data,1,norm_disc_prec_eq)
                CHCKERR('')
                ierr = apply_disc(eq%h_E(:,jd,kd,:,1,0,0),&
                    &ang_deriv_data,res(:,jd,kd-norm_id(1)+1,:),1)
            end do
        end do
        call ang_deriv_data%dealloc()
        
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
    
    ! performs tests on Jacobian in Flux coordinates
    !   - comparing it with the determinant of g_F
    !   - comparing it with the direct formula
    integer function test_jac_F(grid_eq,eq_1,eq_2) result(ierr)
        use num_vars, only: eq_style, use_pol_flux_F
        use grid_utilities, only: trim_grid
        use num_utilities, only: calc_det
        use HELENA_vars, only: h_H_33, RBphi_H
        
        character(*), parameter :: rout_name = 'test_jac_F'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(in), target :: eq_1                             ! flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium
        
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
                            &RBphi_H(kd+grid_eq%i_min-1))                       ! h_H_33 = 1/R^2 and RBphi_H = F are tabulated in eq. grid
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
    
    ! Tests whether g_V is calculated correctly
    integer function test_g_V(grid_eq,eq) result(ierr)
        use grid_utilities, only: trim_grid
        use num_utilities, only: c
        
        character(*), parameter :: rout_name = 'test_g_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       ! metric equilibrium
        
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
    
    ! Tests whether jac_V is calculated correctly
    integer function test_jac_V(grid_eq,eq) result(ierr)
        use grid_utilities, only: trim_grid, nufft
        use VMEC_utilities, only: fourier2real
        use input_utilities, only: get_int, get_log
        use VMEC_vars, only: jac_V_c, jac_V_s, is_asym_V
        use num_vars, only: rank, n_procs, max_deriv
        
        character(*), parameter :: rout_name = 'test_jac_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       ! metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: norm_id_f(2)                                                 ! norm_id transposed to full grid
        integer :: deriv(2)                                                     ! angular derivatives
        logical :: done                                                         ! done
        logical :: do_NUFFT                                                     ! do the NUFFT
        real(dp), allocatable :: res(:,:,:)                                     ! result variable
        real(dp), allocatable :: F(:,:)                                         ! fourier components
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), loc_offset(3)                                    ! total dimensions and local offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether the Jacobian in the VMEC &
            &coords. jac_V is calculated correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim,norm_id)
        CHCKERR('')
        
        ! set norm_id_f for quantities tabulated on full grid
        norm_id_f = grid_eq%i_min+norm_id-1
        
        ! set total and local dimensions and local offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        loc_offset = [0,0,grid_trim%i_min-1]
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r))
        
        ! set some variables
        description = 'Testing calculated with given value for jac_V'
        
        done = .false.
        do while (.not.done)
            call writo('derivative in theta_E?')
            deriv(1) = get_int(lim_lo=0,lim_hi=max_deriv)
            
            call writo('derivative in zeta_E?')
            deriv(2) = get_int(lim_lo=0,lim_hi=max_deriv)
            
            ! get jac_V from VMEC
            ierr = fourier2real(jac_V_c(:,norm_id_f(1):norm_id_f(2)),&
                &jac_V_s(:,norm_id_f(1):norm_id_f(2)),&
                &grid_eq%trigon_factors(:,:,:,norm_id(1):norm_id(2),:),&
                &res,sym=[.true.,is_asym_V],deriv=deriv)
            CHCKERR('')
            
            ! plot difference
            file_name = 'TEST_jac_V'
            call plot_diff_HDF5(res(:,:,:),&
                &eq%jac_E(:,:,norm_id(1):norm_id(2),0,deriv(1),deriv(2)),&
                &file_name,tot_dim,loc_offset,description,output_message=.true.)
            
            ! calculate NUFFT
            call writo('Calculate NUFFT?')
            do_NUFFT = get_log(.false.)
            if (do_NUFFT .and. rank.eq.n_procs-1) then
                call print_ex_2D('VMEC last surface','VMEC_lf',&
                    &res(:,1,grid_trim%loc_n_r),&
                    &x=grid_trim%theta_E(:,1,grid_trim%loc_n_r))
                call print_ex_2D('PB3D last surface','PB3D_lf',&
                    &eq%jac_E(:,1,grid_eq%loc_n_r,0,0,0),&
                    &x=grid_eq%theta_E(:,1,grid_eq%loc_n_r))
                
                ierr = nufft(grid_trim%theta_E(:,1,grid_trim%loc_n_r),&
                    &res(:,1,grid_trim%loc_n_r),F,plot_name='VMEC')
                CHCKERR('')
                deallocate(F)
                ierr = nufft(grid_eq%theta_E(:,1,grid_eq%loc_n_r),&
                    &eq%jac_E(:,1,grid_eq%loc_n_r,0,0,0),F,plot_name='PB3D')
                CHCKERR('')
                deallocate(F)
            end if
            
            call writo('Repeat?')
            done = .not.get_log(.true.)
        end do
        
        ! clean up
        call grid_trim%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_jac_V
    
    ! Tests whether B_F is calculated correctly
    integer function test_B_F(grid_eq,eq_1,eq_2) result(ierr)
        use num_vars, only: eq_style
        use grid_utilities, only: trim_grid
        use num_utilities, only: c
        use VMEC_utilities, only: fourier2real
        use VMEC_vars, only: B_V_sub_s, B_V_sub_c, B_V_c, B_V_s, is_asym_V
        
        character(*), parameter :: rout_name = 'test_B_F'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium
        
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
    
    ! performs tests on pressure balance
    !   - mu_0 D2p = 1/J (D3 B_2 - D2 B_3)
    !   - mu_0 J D3p = 0 => D3 B_1 = D1 B_3
    ! (working in the (modified) Flux coordinates (alpha,psi,theta)_F)
    integer function test_p(grid_eq,eq_1,eq_2) result(ierr)
        use num_utilities, only: c
        use grid_utilities, only: trim_grid
        use eq_vars, only: vac_perm
        use num_vars, only: eq_style
        use HELENA_vars, only: RBphi_H, h_H_11, h_H_12
        use splines, only: spline3
        
        character(*), parameter :: rout_name = 'test_p'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        
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
                    ierr = spline3(grid_trim%loc_r_E,&
                        &eq_1%q_saf_E(norm_id(1):norm_id(2),0)*&
                        &RBphi_H(norm_id_f(1):norm_id_f(2))+&
                        &eq_1%q_saf_E(norm_id(1):norm_id(2),0)*&
                        &h_H_11(id,norm_id_f(1):norm_id_f(2))/&
                        &RBphi_H(norm_id_f(1):norm_id_f(2)),&
                        &grid_trim%loc_r_E,dynew=res(id,jd,:,2))
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
                    ierr = spline3(grid_eq%theta_E(:,jd,kd),&
                        &h_H_12(:,kd+grid_eq%i_min-1),grid_eq%theta_E(:,jd,kd),&
                        &dynew=res(:,jd,kd-norm_id(1)+1,2))
                    CHCKERR('')
                end do
                res(:,:,kd-norm_id(1)+1,2) = &
                    &-eq_1%q_saf_E(kd,0)/RBphi_H(kd+grid_eq%i_min-1) * &
                    &res(:,:,kd-norm_id(1)+1,2) + &
                    &RBphi_H(kd+grid_eq%i_min-1) * &
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
