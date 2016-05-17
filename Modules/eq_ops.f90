!------------------------------------------------------------------------------!
!   Operations on the equilibrium variables                                    !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln, max_deriv
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type, eq_2_type
    use num_utilities, only: check_deriv
    
    implicit none
    private
    public calc_eq, calc_derived_q, calc_normalization_const, normalize_input, &
        &print_output_eq, flux_q_plot
#if ldebug
    public debug_calc_derived_q
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_derived_q = .false.                                   ! plot debug information for calc_derived_q
#endif
    
    ! interfaces
    interface calc_eq
        module procedure calc_eq_1, calc_eq_2
    end interface
    interface print_output_eq
        module procedure print_output_eq_1, print_output_eq_2
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
            &norm_disc_prec_eq, use_pol_flux_E, use_pol_flux_F
        use grid_utilities, only: setup_deriv_data, apply_disc
        use eq_vars, only: rho_0
        use num_utilities, only: derivs
        use eq_utilities, only: calc_F_derivs
        
        character(*), parameter :: rout_name = 'calc_eq_1'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid
        type(eq_1_type), intent(inout) :: eq                                    ! flux equilibrium variables
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(disc_type), allocatable :: norm_deriv_data(:)                      ! data for normal derivatives
        
        ! variables in full equilibrium grid
        real(dp), allocatable :: flux_p_E_full(:,:)                             ! poloidal flux in full equilibrium coordinates
        real(dp), allocatable :: flux_t_E_full(:,:)                             ! toroidal flux in full equilibrium coordinates
        real(dp), allocatable :: pres_E_full(:,:)                               ! pressure full equilibrium coordinates
        real(dp), allocatable :: q_saf_E_full(:,:)                              ! safety factor in full equilibrium coordinates
        real(dp), allocatable :: rot_t_E_full(:,:)                              ! rotational transform in full equilibrium coordinates
        
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
                ierr = calc_flux_q_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_flux_q_HEL()
                CHCKERR('')
        end select
        
        ! clean up
        do id = 1,size(norm_deriv_data)
            call norm_deriv_data(id)%dealloc()
        end do
        deallocate(norm_deriv_data)
        
        ! take local variables
        eq%flux_p_E = flux_p_E_full(grid_eq%i_min:grid_eq%i_max,:)
        eq%flux_t_E = flux_t_E_full(grid_eq%i_min:grid_eq%i_max,:)
        eq%pres_E = pres_E_full(grid_eq%i_min:grid_eq%i_max,:)
        eq%q_saf_E = q_saf_E_full(grid_eq%i_min:grid_eq%i_max,:)
        eq%rot_t_E = rot_t_E_full(grid_eq%i_min:grid_eq%i_max,:)
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
        integer function calc_flux_q_VMEC() result(ierr)
            use VMEC, only: rot_t_V, flux_t_V, Dflux_t_V, flux_p_V, Dflux_p_V, &
                &pres_V
            use eq_vars, only: max_flux_E
            use grid_vars, only: n_r_in, n_r_eq
            
            character(*), parameter :: rout_name = 'calc_flux_q_VMEC'
            
            ! local variables
            integer :: kd                                                       ! counter
            
            ! initialize ierr
            ierr = 0
            
            ! allocate full variables
            allocate(flux_p_E_full(grid_eq%n(3),0:max_deriv+2))
            allocate(flux_t_E_full(grid_eq%n(3),0:max_deriv+2))
            allocate(pres_E_full(grid_eq%n(3),0:max_deriv+1))
            allocate(q_saf_E_full(grid_eq%n(3),0:max_deriv+1))
            allocate(rot_t_E_full(grid_eq%n(3),0:max_deriv+1))
            
            ! calculate data  for normal derivatives with  equidistant grid with
            ! grid size 1/(n_r_in-1), but n_r_eq points
            allocate(norm_deriv_data(max_deriv+2))
            do kd = 1,max_deriv+2
                ierr = setup_deriv_data(1._dp/(n_r_in-1),n_r_eq,&
                    &norm_deriv_data(kd),kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! calculate derivatives of toroidal flux
            flux_t_E_full(:,0) = flux_t_V
            flux_t_E_full(:,1) = Dflux_t_V
            do kd = 2,max_deriv+2
                ierr = apply_disc(flux_t_E_full(:,1),&
                    &norm_deriv_data(kd-1),flux_t_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! calculate derivatives of poloidal flux
            flux_p_E_full(:,0) = flux_p_V
            flux_p_E_full(:,1) = Dflux_p_V
            do kd = 2,max_deriv+2
                ierr = apply_disc(flux_p_E_full(:,1),&
                    &norm_deriv_data(kd-1),flux_p_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! pressure: copy from VMEC and derive
            pres_E_full(:,0) = pres_V
            do kd = 1, max_deriv+1
                ierr = apply_disc(pres_E_full(:,0),norm_deriv_data(kd),&
                    &pres_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! safety factor
            q_saf_E_full(:,0) = 1._dp/rot_t_V
            do kd = 1, max_deriv+1
                ierr = apply_disc(q_saf_E_full(:,0),norm_deriv_data(kd),&
                    &q_saf_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! rot. transform
            rot_t_E_full(:,0) = rot_t_V
            do kd = 1, max_deriv+1
                ierr = apply_disc(rot_t_E_full(:,0),norm_deriv_data(kd),&
                    &rot_t_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! max flux and  normal coord. of eq grid  in Equilibrium coordinates
            ! (uses poloidal flux by default)
            if (use_pol_flux_E) then
                grid_eq%r_E = flux_p_E_full(:,0)/max_flux_E
            else
                grid_eq%r_E = flux_t_E_full(:,0)/max_flux_E
            end if
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_E_full(:,0)/(2*pi)                         ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = - flux_t_E_full(:,0)/(2*pi)                       ! psi_F = flux_t/2pi, conversion VMEC LH -> PB3D RH
            end if
        end function calc_flux_q_VMEC
        
        ! HELENA version
        ! The HELENA normal coord. is the poloidal flux divided by 2pi
        integer function calc_flux_q_HEL() result(ierr)
            use HELENA_vars, only: qs_H, flux_p_H, pres_H
            use num_utilities, only: calc_int
            
            character(*), parameter :: rout_name = 'calc_flux_q_HEL'
            
            ! local variables
            integer :: kd                                                       ! counter
            
            ! initialize ierr
            ierr = 0
            
            ! allocate full variables
            allocate(flux_p_E_full(grid_eq%n(3),0:max_deriv+1))
            allocate(flux_t_E_full(grid_eq%n(3),0:max_deriv+1))
            allocate(pres_E_full(grid_eq%n(3),0:max_deriv+1))
            allocate(q_saf_E_full(grid_eq%n(3),0:max_deriv+1))
            allocate(rot_t_E_full(grid_eq%n(3),0:max_deriv+1))
            
            ! calculate data for normal derivatives with flux_p_H/2pi
            allocate(norm_deriv_data(max_deriv+1))
            do kd = 1,max_deriv+1
                ierr = setup_deriv_data(flux_p_H/(2*pi),norm_deriv_data(kd),&
                    &kd,norm_disc_prec_eq)
                CHCKERR('')
            end do
            
            ! calculate derivatives of poloidal flux
            flux_p_E_full(:,0) = flux_p_H
            do kd = 1,max_deriv+1
                ierr = apply_disc(flux_p_E_full(:,0),norm_deriv_data(kd),&
                    &flux_p_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! calculate toroidal flux and derivatives
            flux_t_E_full(:,1) = qs_H*flux_p_E_full(:,1)
            ierr = calc_int(flux_t_E_full(:,1),flux_p_H/(2*pi),&
                &flux_t_E_full(:,0))
            CHCKERR('')
            do kd = 2,max_deriv+1
                ierr = apply_disc(flux_t_E_full(:,1),norm_deriv_data(kd-1),&
                    &flux_t_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! pressure: copy from HELENA and derive
            pres_E_full(:,0) = pres_H
            do kd = 1, max_deriv+1
                ierr = apply_disc(pres_E_full(:,0),norm_deriv_data(kd),&
                    &pres_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! safety factor
            q_saf_E_full(:,0) = qs_H
            do kd = 1, max_deriv+1
                ierr = apply_disc(q_saf_E_full(:,0),norm_deriv_data(kd),&
                    &q_saf_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! rot. transform
            rot_t_E_full(:,0) = 1._dp/qs_H
            do kd = 1, max_deriv+1
                ierr = apply_disc(rot_t_E_full(:,0),norm_deriv_data(kd),&
                    &rot_t_E_full(:,kd))
                CHCKERR('')
            end do
            
            ! max flux and  normal coord. of eq grid  in Equilibrium coordinates
            ! (uses poloidal flux by default)
            grid_eq%r_E = flux_p_E_full(:,0)/(2*pi)
            
            ! max flux and normal coord. of eq grid in Flux coordinates
            if (use_pol_flux_F) then
                grid_eq%r_F = flux_p_E_full(:,0)/(2*pi)                         ! psi_F = flux_p/2pi
            else
                grid_eq%r_F = flux_t_E_full(:,0)/(2*pi)                         ! psi_F = flux_t/2pi
            end if
            
            !!!! export for VMEC port
            !!!call write_flux_q_in_file_for_VMEC()
        end function calc_flux_q_HEL
        
        ! plots flux quantities in file for VMEC port
        subroutine write_flux_q_in_file_for_VMEC()
            use eq_vars, only: pres_0
            use files_utilities, only: nextunit
            
            ! local variables
            integer :: kd                                                       ! counter
            integer :: file_i                                                   ! file to output flux quantities for VMEC export
            character(len=max_str_ln) :: file_name                              ! name of file
            
            file_name = 'flux_quantities.dat'
            
            call writo('Writing output to file '//trim(file_name))
            call writo('This can be used for VMEC porting')
            
            open(unit=nextunit(file_i),file=trim(file_name))
            write(file_i,*) '# for export to VMEC'
            write(file_i,*) '# ------------------'
            write(file_i,*) '# total enclosed toroidal flux: ', &
                &flux_t_E_full(grid_eq%n(3),0)
            write(file_i,*) '#'
            write(file_i,*) '# pressure and rotational transform for 2D plot'
            do kd = 1,grid_eq%n(3)
                write(file_i,'(ES23.16,A,ES23.16,A,ES23.16)') &
                    &flux_t_E_full(kd,0)/flux_t_E_full(grid_eq%n(3),0), ' ', &
                    &pres_E_full(kd,0)*pres_0, ' ', rot_t_E_full(kd,0)
            end do
            write(file_i,*) ''
            write(file_i,*) ''
            
            write(file_i,*) '# toroidal flux'
            do kd = 1,grid_eq%n(3)
                write(file_i,'(ES23.16)') &
                    &flux_t_E_full(kd,0)/flux_t_E_full(grid_eq%n(3),0)
            end do
            write(file_i,*) ''
            write(file_i,*) ''
            write(file_i,*) '# pressure'
            do kd = 1,grid_eq%n(3)
                write(file_i,'(ES23.16)') pres_E_full(kd,0)*pres_0
            end do
            write(file_i,*) ''
            write(file_i,*) ''
            write(file_i,*) '# rotational transform'
            do kd = 1,grid_eq%n(3)
                write(file_i,'(ES23.16)') rot_t_E_full(kd,0)
            end do
            write(file_i,*) ''
            write(file_i,*) ''
            
            write(file_i,*) '# toroidal flux once every 4 values'
            do kd = 0,(grid_eq%n(3)-1)/4
                write(file_i,'(ES23.16)') &
                    &flux_t_E_full(1+kd*4,0)/flux_t_E_full(grid_eq%n(3),0)
            end do
            write(file_i,*) ''
            write(file_i,*) ''
            write(file_i,*) '# pressure once every 4 values'
            do kd = 0,(grid_eq%n(3)-1)/4
                write(file_i,'(ES23.16)') pres_E_full(1+kd*4,0)*pres_0
            end do
            write(file_i,*) ''
            write(file_i,*) ''
            write(file_i,*) '# rotational transform once every 4 values'
            do kd = 0,(grid_eq%n(3)-1)/4
                write(file_i,'(ES23.16)') rot_t_E_full(1+kd*4,0)
            end do
            
            ! close file
            close(file_i)
        end subroutine write_flux_q_in_file_for_VMEC
    end function calc_eq_1
    integer function calc_eq_2(grid_eq,eq_1,eq_2,dealloc_vars) result(ierr)     ! metric version
        use num_vars, only: eq_style
        use num_utilities, only: derivs, c
        use eq_utilities, only: calc_inv_met, calc_F_derivs
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
                ierr = prepare_RZL(grid_eq)
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
                        ierr = test_metrics_H(grid_eq%n(3))
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
    end function calc_eq_2

    ! plots the flux quantities in the solution grid
    !   safety factor q_saf
    !   rotational transform rot
    !   pressure pres
    !   poloidal flux flux_p
    !   toroidal flux flux_t
    integer function flux_q_plot(grid_eq,eq) result(ierr)
        use num_vars, only: rank, no_plots
        use grid_utilities, only: trim_grid
        use MPI_utilities, only: get_ser_var
        
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
        
        ! plot using GNUPlot
        ierr = flux_q_plot_GP()
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        
        call lvl_ud(-1)
    contains
        ! plots the flux quantities in HDF5
        integer function flux_q_plot_HDF5() result(ierr)
            use num_vars, only: eq_style
            use output_ops, only: plot_HDF5
            use grid_utilities, only: calc_XYZ_grid, extend_grid_E, trim_grid
            use VMEC, only: calc_trigon_factors
            
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
            
            ! print the output using HDF5
            call plot_HDF5(plot_titles,file_name,f_plot,plot_dim,plot_offset,&
                &X_plot,Y_plot,Z_plot,col=1,description='Flux quantities')
            
            ! deallocate and destroy grid
            deallocate(Y_plot_2D)
            deallocate(X_plot,Y_plot,Z_plot,f_plot)
            call grid_plot%dealloc()
        end function flux_q_plot_HDF5
        
        ! plots the pressure and fluxes in GNUplot
        integer function flux_q_plot_GP() result(ierr)
            use eq_vars, only: max_flux_F, pres_0, psi_0
            use num_vars, only: use_normalization
            
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
    end function flux_q_plot
    
    ! prepare the cosine  and sine factors that are used  in the inverse Fourier
    ! transformation of R, Z and L and calculate the derivatives.
    integer function prepare_RZL(grid) result(ierr)
        use num_vars, only: max_deriv, norm_disc_prec_eq
        use VMEC, only: calc_trigon_factors, &
            &R_V_c, Z_V_s, L_V_s, R_V_s, Z_V_c, L_V_c
        use grid_utilities, only: setup_deriv_data, apply_disc
        use grid_vars, only: n_r_in
        
        character(*), parameter :: rout_name = 'prepare_RZL'
        
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid for which to prepare the trigonometric factors
        
        ! local variables
        integer :: kd                                                           ! counter
        type(disc_type) :: norm_deriv_data                                      ! data for normal derivative
        
        ! initialize ierr
        ierr = 0
        
        ! normal derivatives of these factors
        ! The VMEC normal coord. is  the toroidal (or poloidal) flux, normalized
        ! wrt.  to  the  maximum  flux,  equidistantly,  so  the  step  size  is
        ! 1/(n_r_in-1).
        do kd = 1,max_deriv+1
            ierr = setup_deriv_data(1._dp/(n_r_in-1),grid%n(3),&
                &norm_deriv_data,kd,norm_disc_prec_eq)
            CHCKERR('')
            ierr = apply_disc(R_V_c(:,:,0),norm_deriv_data,R_V_c(:,:,kd),2)
            CHCKERR('')
            ierr = apply_disc(R_V_s(:,:,0),norm_deriv_data,R_V_s(:,:,kd),2)
            CHCKERR('')
            ierr = apply_disc(Z_V_c(:,:,0),norm_deriv_data,Z_V_c(:,:,kd),2)
            CHCKERR('')
            ierr = apply_disc(Z_V_s(:,:,0),norm_deriv_data,Z_V_s(:,:,kd),2)
            CHCKERR('')
            ierr = apply_disc(L_V_c(:,:,0),norm_deriv_data,L_V_c(:,:,kd),2)
            CHCKERR('')
            ierr = apply_disc(L_V_s(:,:,0),norm_deriv_data,L_V_s(:,:,kd),2)
            CHCKERR('')
        end do
        call norm_deriv_data%dealloc()
        
        ! calculate trigonometric factors using theta_E and zeta_E
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
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(Z_V_c(:,grid%i_min:grid%i_max,deriv(1)),&
            &Z_V_s(:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(L_V_c(:,grid%i_min:grid%i_max,deriv(1)),&
            &L_V_s(:,grid%i_min:grid%i_max,deriv(1)),grid%trigon_factors,&
            &eq%L_E(:,:,:,deriv(1),deriv(2),deriv(3)),[deriv(2),deriv(3)])
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
        
        do id = 1, size(deriv,2)
            ierr = calc_g_V_ind(eq,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_V_arr
    
    ! calculate the  metric coefficients in the  equilibrium H(ELENA) coordinate
    ! system using the HELENA output
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
        use num_vars, only: norm_disc_prec_eq, eq_style
#if ldebug
        use num_vars, only: use_pol_flux_F
        use grid_utilities, only: trim_grid, calc_XYZ_grid, setup_deriv_data, &
            &apply_disc
        use VMEC, only: calc_trigon_factors
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
        type(grid_type) :: grid_eq_trim                                         ! trimmed equilibrium grid
        real(dp), allocatable :: D3sigma(:,:,:)                                 ! D_theta sigma
        real(dp), allocatable :: D3sigma_ALT(:,:,:)                             ! alternative D_theta sigma
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle theta_F or zeta_F
        real(dp), allocatable :: X_plot(:,:,:)                                  ! X of plot
        real(dp), allocatable :: Y_plot(:,:,:)                                  ! Y of plot
        real(dp), allocatable :: Z_plot(:,:,:)                                  ! Z of plot
        type(disc_type) :: ang_deriv_data                                      ! data for angular derivatives
        integer :: istat                                                        ! status
        integer :: jd                                                           ! counter
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
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
            &(D1g23/J - g23*D1J/J**2 - D2g13/J + g13*D2J/J**2)
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
            istat = trim_grid(grid_eq,grid_eq_trim,norm_id)
            CHCKSTT
            
            ! if VMEC, calculate trigonometric factors of plot grid
            if (eq_style.eq.1) then
                istat = calc_trigon_factors(grid_eq_trim%theta_E,&
                    &grid_eq_trim%zeta_E,grid_eq_trim%trigon_factors)
                CHCKSTT
            end if
            
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
                    istat = setup_deriv_data(ang_par_F(:,jd,kd),ang_deriv_data,&
                        &1,norm_disc_prec_eq)
                    CHCKSTT
                    istat = apply_disc(eq_2%sigma(:,jd,kd),ang_deriv_data,&
                        &D3sigma(:,jd,kd-norm_id(1)+1))
                    CHCKSTT
                end do
            end do
            call ang_deriv_data%dealloc()
            
            ! calculate alternatively derived sigma
            do kd = norm_id(1),norm_id(2)
                D3sigma_ALT(:,:,kd-norm_id(1)+1) = &
                    &-2*eq_1%pres_FD(kd,1)*eq_2%kappa_g(:,:,kd)*J(:,:,kd)
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
            istat = calc_XYZ_grid(grid_eq_trim,grid_eq_trim,X_plot,Y_plot,&
                &Z_plot)
            CHCKSTT
            
            ! plot shear
            call plot_HDF5('shear','TEST_shear',&
                &eq_2%S(:,:,norm_id(1):norm_id(2)),tot_dim=grid_eq_trim%n,&
                &loc_offset=[0,0,grid_eq_trim%i_min-1],x=X_plot,y=Y_plot,&
                &z=Z_plot)
            
            ! plot sigma
            call plot_HDF5('sigma','TEST_sigma',&
                &eq_2%sigma(:,:,norm_id(1):norm_id(2)),tot_dim=grid_eq_trim%n,&
                &loc_offset=[0,0,grid_eq_trim%i_min-1],x=X_plot,y=Y_plot,&
                &z=Z_plot)
            
            ! plot kappa_n
            call plot_HDF5('kappa_n','TEST_kappa_n',&
                &eq_2%kappa_n(:,:,norm_id(1):norm_id(2)),&
                &tot_dim=grid_eq_trim%n,loc_offset=[0,0,grid_eq_trim%i_min-1],&
                &x=X_plot,y=Y_plot,z=Z_plot)
            
            ! plot kappa_g
            call plot_HDF5('kappa_g','TEST_kappa_g',&
                &eq_2%kappa_g(:,:,norm_id(1):norm_id(2)),&
                &tot_dim=grid_eq_trim%n,loc_offset=[0,0,grid_eq_trim%i_min-1],&
                &x=X_plot,y=Y_plot,z=Z_plot)
            
            ! clean up
            nullify(ang_par_F)
            call grid_eq_trim%dealloc()
            
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
            
            ! user output
            call writo('R_0    = '//trim(r2str(R_0))//' m')
            call writo('rho_0  = '//trim(r2str(rho_0))//' kg/m^3')
            call writo('B_0    = '//trim(r2str(B_0))//' T')
            call writo('pres_0 = '//trim(r2str(pres_0))//' Pa')
            call writo('psi_0  = '//trim(r2str(psi_0))//' Tm^2')
            call writo('mu_0   = '//trim(r2str(mu_0_original))//' Tm/A')
            call writo('T_0    = '//trim(r2str(T_0))//' s')
            
            ! user output
            call lvl_ud(-1)
            call writo('Normalization constants calculated')
        else if (rich_restart_lvl.eq.1) then                                    ! only for first Richardson leevel
            ! user output
            call writo('Normalization not used')
        end if
    contains 
        ! VMEC version
        !   R_0:    major radius (= average R on axis)
        !   B_0:    B on magnetic axis (theta = zeta = 0)
        !   pres_0: reference pressure (= B_0^2/mu_0)
        !   psi_0:  reference flux (= R_0^2 B_0)
        !   rho_0:  reference mass density
        ! Note: The  orthodox way of doing  this is by setting  B_0 the toroidal
        ! field on the magnetic axis, and calculating pres_0 from this, which is
        ! not done here.
        ! Note that  rho_0 is  not given  through by  the equilibrium  codes and
        ! should be user-supplied
        subroutine calc_normalization_const_VMEC
            use VMEC, only: R_V_c, B_0_V
            
            ! set the major  radius as the average value of  R_V on the magnetic
            ! axis
            if (R_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                R_0 = R_V_c(1,1,0)
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            
            ! set the reference value for B_0 = B_0_V
            if (B_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                B_0 = B_0_V
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            
            ! set pres_0 from B_0
            if (pres_0.ge.huge(1._dp)) then                                     ! user did not provide a value
                pres_0 = B_0**2/mu_0_original
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            
            ! set reference flux
            if (psi_0.ge.huge(1._dp)) then                                      ! user did not provide a value
                psi_0 = R_0**2 * B_0
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            
            ! rho_0 is set up through an input variable with the same name
            
            ! set Alfven time
            if (T_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                T_0 = sqrt(mu_0_original*rho_0)*R_0/B_0
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
        end subroutine calc_normalization_const_VMEC
        
        ! HELENA version
        ! The  MISHKA normalization  with R_m  = 1  m,  B_m =  1 T  is taken  by
        ! default, see "read_VMEC" for more information.
        subroutine calc_normalization_const_HEL
            if (R_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                R_0 = 1._dp
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            if (B_0.ge.huge(1._dp)) then                                        ! user did not provide a value
                B_0 = 1._dp
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            if (pres_0.ge.huge(1._dp)) then                                     ! user did not provide a value
                pres_0 = B_0**2/mu_0_original
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            if (psi_0.ge.huge(1._dp)) then                                      ! user did not provide a value
                psi_0 = 1._dp
            else
                nr_overriden_const = nr_overriden_const + 1
            end if
            if (T_0.ge.huge(1._dp)) T_0 = sqrt(mu_0_original*rho_0)*R_0/B_0     ! only if user did not provide a value
        end subroutine calc_normalization_const_HEL
    end subroutine calc_normalization_const
    
    ! Normalize input quantities.
    subroutine normalize_input()
        use num_vars, only: use_normalization, eq_style, mu_0_original
        use VMEC, only: normalize_VMEC
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
    
    ! Print equilibrium quantities to an output file:
    !   - flux:     pres_FD, q_saf_FD, rot_t_FD, flux_p_FD, flux_t_FD, rho, S,
    !               kappa_n, kappa_g, sigma
    !   - metric:   g_FD, h_FD, jac_FD
    ! If "rich_lvl" is  provided, "_R_rich_lvl" is appended to the  data name if
    ! it is > 0 (only for eq_2), and similarly for "eq_job" through "_E_eq_job".
    ! This counts only for eq_2.
    ! Note: The equilibrium quantities are outputted in Flux coordinates.
    ! Note: The  metric equilibrium  quantities can be  deallocated on  the fly,
    ! which is useful if this routine is  followed by a deallocation any way, so
    ! that memory usage does not almost double.
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
        ierr = print_HDF5_arrs(eq_1D(1:id-1),PB3D_name,trim(data_name))
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        call dealloc_var_1D(eq_1D)
        nullify(eq_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_eq_1
    integer function print_output_eq_2(grid_eq,eq,data_name,rich_lvl,eq_job,&
        &dealloc_vars) result(ierr)                                             ! metric version
        use num_vars, only: PB3D_name_eq
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
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to print
        logical, intent(in), optional :: dealloc_vars                           ! deallocate variables on the fly after writing
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        type(var_1D_type), allocatable, target :: eq_1D(:)                      ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: eq_1D_loc => null()                       ! local element in eq_1D
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: id                                                           ! counter
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
        eq_1D_loc%tot_i_max = [grid_trim%n,6,size(eq%g_FD,5)-1,&
            &size(eq%g_FD,6)-1,size(eq%g_FD,7)-1]
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min,1,0,0,0]
        eq_1D_loc%loc_i_max = [grid_trim%n(1),grid_trim%n(2),&
            &grid_trim%i_max,6,size(eq%g_FD,5)-1,size(eq%g_FD,6)-1,&
            &size(eq%g_FD,7)-1]
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
        eq_1D_loc%tot_i_max = [grid_trim%n,6,size(eq%h_FD,5)-1,&
            &size(eq%h_FD,6)-1,size(eq%h_FD,7)-1]
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min,1,0,0,0]
        eq_1D_loc%loc_i_max = [grid_trim%n(1),grid_trim%n(2),&
            &grid_trim%i_max,6,size(eq%h_FD,5)-1,size(eq%h_FD,6)-1,&
            &size(eq%h_FD,7)-1]
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
        eq_1D_loc%tot_i_max = [grid_trim%n,size(eq%jac_FD,4)-1,&
            &size(eq%jac_FD,5)-1,size(eq%jac_FD,6)-1]
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min,0,0,0]
        eq_1D_loc%loc_i_max = [grid_trim%n(1),grid_trim%n(2),&
            &grid_trim%i_max,size(eq%jac_FD,4)-1,size(eq%jac_FD,5)-1,&
            &size(eq%jac_FD,6)-1]
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
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = &
            &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
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
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = &
            &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
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
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = &
            &[grid_trim%n(1),grid_trim%n(2),grid_trim%i_max]
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
        eq_1D_loc%tot_i_max = grid_trim%n
        eq_1D_loc%loc_i_min = [1,1,grid_trim%i_min]
        eq_1D_loc%loc_i_max = [grid_trim%n(1:2),grid_trim%i_max]
        loc_size = size(eq%sigma(:,:,norm_id(1):norm_id(2)))
        allocate(eq_1D_loc%p(loc_size))
        eq_1D_loc%p = reshape(eq%sigma(:,:,norm_id(1):norm_id(2)),&
            &[loc_size])
        if (dealloc_vars_loc) deallocate(eq%sigma)
        
        ! write
        ierr = print_HDF5_arrs(eq_1D(1:id-1),PB3D_name_eq,trim(data_name),&
            &rich_lvl=rich_lvl,eq_job=eq_job)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        call dealloc_var_1D(eq_1D)
        nullify(eq_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_eq_2

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
        use grid_utilities, only: trim_grid
        use VMEC, only: fourier2real, &
            &jac_V_c, jac_V_s
        
        character(*), parameter :: rout_name = 'test_jac_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_2_type), intent(in) :: eq                                       ! metric equilibrium
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: norm_id_f(2)                                                 ! norm_id transposed to full grid
        real(dp), allocatable :: res(:,:,:)                                     ! result variable
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
        
        ! get jac_V from VMEC
        ierr = fourier2real(jac_V_c(:,norm_id_f(1):norm_id_f(2)),&
            &jac_V_s(:,norm_id_f(1):norm_id_f(2)),&
            &grid_eq%trigon_factors(:,:,:,norm_id(1):norm_id(2),:),&
            &res,[0,0])
        CHCKERR('')
        
        ! user output
        call writo('Testing jac_V')
        call lvl_ud(1)
        
        ! set some variables
        file_name = 'TEST_jac_V'
        description = 'Testing calculated with given value for jac_V'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &eq%jac_E(:,:,norm_id(1):norm_id(2),0,0,0),file_name,&
            &tot_dim,loc_offset,description,output_message=.true.)
        
        call lvl_ud(-1)
        
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
        use VMEC, only: fourier2real, &
            &B_V_sub_s, B_V_sub_c, B_V_c, B_V_s
        
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
                        &norm_id(2),:),res(:,:,:,id),[0,0])
                    CHCKERR('')
                end do
                ierr = fourier2real(B_V_c(:,norm_id_f(1):norm_id_f(2)),&
                    &B_V_s(:,norm_id_f(1):norm_id_f(2)),&
                    &grid_eq%trigon_factors(:,:,:,norm_id(1):norm_id(2),:),&
                    &res(:,:,:,4),[0,0])
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
        use grid_utilities, only: trim_grid, setup_deriv_data, apply_disc
        use eq_vars, only: vac_perm
        use num_vars, only: eq_style, norm_disc_prec_eq
        use HELENA_vars, only: RBphi_H, h_H_11, h_H_12
        
        character(*), parameter :: rout_name = 'test_p'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: norm_id_f(2)                                                 ! norm_id transposed to full grid
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        type(disc_type) :: norm_deriv_data                                      ! data for normal derivative
        type(disc_type) :: ang_deriv_data                                       ! data for angular derivative
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
        call plot_HDF5('var','TEST_Dg_FD_23',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,1))
        call plot_HDF5('var','TEST_g_FD_23',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([2,3],.true.),0,0,0))
        call plot_HDF5('var','TEST_Dg_FD_33',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,1,0))
        call plot_HDF5('var','TEST_g_FD_33',eq_2%g_FD(:,:,norm_id(1):norm_id(2),c([3,3],.true.),0,0,0))
        call plot_HDF5('var','TEST_Dg_FD',eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,1))
        call plot_HDF5('var','TEST_g_FD',eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0))
        
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
            ierr = setup_deriv_data(grid_trim%loc_r_E,norm_deriv_data,1,&
                &norm_disc_prec_eq)
            CHCKERR('')
            do jd = 1,grid_trim%n(2)
                do id = 1,grid_trim%n(1)
                    ierr = apply_disc(eq_1%q_saf_E(norm_id(1):norm_id(2),0)*&
                        &RBphi_H(norm_id_f(1):norm_id_f(2))+&
                        &eq_1%q_saf_E(norm_id(1):norm_id(2),0)*&
                        &h_H_11(id,norm_id_f(1):norm_id_f(2))/&
                        &RBphi_H(norm_id_f(1):norm_id_f(2)),&
                        &norm_deriv_data,res(id,jd,:,2))
                    CHCKERR('')
                end do
            end do
            call norm_deriv_data%dealloc()
            
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
                    ierr = setup_deriv_data(grid_eq%theta_E(:,jd,kd),&
                        &ang_deriv_data,1,norm_disc_prec_eq)
                    CHCKERR('')
                    ierr = apply_disc(h_H_12(:,kd+grid_eq%i_min-1),&
                        &ang_deriv_data,res(:,jd,kd-norm_id(1)+1,2))
                end do
                res(:,:,kd-norm_id(1)+1,2) = &
                    &-eq_1%q_saf_E(kd,0)/RBphi_H(kd+grid_eq%i_min-1) * &
                    &res(:,:,kd-norm_id(1)+1,2) + &
                    &RBphi_H(kd+grid_eq%i_min-1) * &
                    &eq_1%q_saf_E(kd,1)
            end do
            call ang_deriv_data%dealloc()
            
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
