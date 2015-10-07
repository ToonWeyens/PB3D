!------------------------------------------------------------------------------!
!   Variables from PB3D output:                                                !
!       - Equilibrium variables                                                !
!       - Perturbation variables                                               !
!------------------------------------------------------------------------------!
module PB3D_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type
    use PB3D_vars, only: PB3D_type
    use HDF5_ops, only: var_1D

    implicit none
    private
    public read_PB3D, reconstruct_PB3D, retrieve_var_1D_id
    
    ! interfaces
    interface conv_1D2ND
        module procedure conv_1D2ND_1D, conv_1D2ND_2D, conv_1D2ND_3D, &
            &conv_1D2ND_4D, &
            !&conv_1D2ND_5D, &
            &conv_1D2ND_6D, conv_1D2ND_7D
    end interface
    
contains
    ! reads PB3D output from user-provided input file
    integer function read_PB3D(read_eq,read_X_1,read_X_2,read_sol) result(ierr)
        use num_vars, only: eq_style, PB3D_name
        use HDF5_ops, only: read_HDF5_arrs
        use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B, vars_1D_X_1, &
            &vars_1D_X_2, vars_1D_sol
        
        character(*), parameter :: rout_name = 'read_PB3D'
        
        ! input / output
        logical, intent(in) :: read_eq, read_X_1, read_X_2, read_sol            ! which files to read
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set common error message
        err_msg = 'Maybe you wanted to use Richardson extrapolation? Not &
            &yet implemented.'
        
        if (read_eq) then
            call writo('Reading equilibrium variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_eq,PB3D_name,'eq')
            CHCKERR('')
            ierr = set_eq_style()
            CHCKERR('')
            if (eq_style.eq.2) then                                             ! only read eq_B for HELENA
                ierr = read_HDF5_arrs(vars_1D_eq_B,PB3D_name,'eq_B')
                CHCKERR('')
            else
                allocate(vars_1D_eq_B(0))                                       ! not needed for other equilibrium styles
            end if
            call lvl_ud(-1)
            call writo('Equilbrium variables read')
        end if
        
        if (read_X_1) then
            call writo('Reading vectorial perturbation variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_X_1,PB3D_name,'X_1')
            CHCKERR(err_msg)
            call lvl_ud(-1)
            call writo('Vectorial perturbation variables read')
        end if
        
        if (read_X_2) then
            call writo('Reading tensorial perturbation variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_X_2,PB3D_name,'X_2')
            CHCKERR(err_msg)
            call lvl_ud(-1)
            call writo('Tensorial perturbation variables read')
        end if
        
        if (read_sol) then
            call writo('Reading solution variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_sol,PB3D_name,'sol')
            CHCKERR(err_msg)
            call lvl_ud(-1)
            call writo('Solution variables read')
        end if
    contains
        ! sets equilibrium style
        integer function set_eq_style() result(ierr)
            use num_vars, only: eq_style, rho_style
            
            character(*), parameter :: rout_name = 'set_eq_style'
            
            ! local variables
            integer :: misc_eq_id                                               ! index of misc_eq
            real(dp), allocatable :: dum_1D(:)                                  ! dummy variable
            
            ! initialize ierr
            ierr = 0
            
            ! retrieve variable
            ierr = retrieve_var_1D_id(vars_1D_eq,'misc_eq',misc_eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D_eq(misc_eq_id),dum_1D)
            eq_style = nint(dum_1D(2))
            rho_style = nint(dum_1D(3))
            select case (eq_style)
                case (1)
                    call writo('This is a VMEC equilibrium')
                case (2)
                    call writo('This is a HELENA equilibrium')
                case default
                    ierr = 1
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    CHCKERR(err_msg)
            end select
        end function set_eq_style
    end function read_PB3D
    
    ! Reconstructs the PB3D variables: eq,  X and/or sol variables, as indicated
    ! through rec_eq, rec_X_1, rec_X_2 and rec_sol.
    ! Additionally,  if  the variables  'eq_limits'  and/or  'X_limits' are  not
    ! provided, the full normal range is taken for every variable.
    ! Note that despite the name, the eq limits are used for eq- and X-variables
    ! and that the X limits for the sol-variables!
    ! Note that  both the  X and the  sol variables need  the eq  variables, but
    ! these necessary variables do not have to be calculated here: They can also
    ! be provided.
    integer function reconstruct_PB3D(rec_eq,rec_X_1,rec_X_2,rec_sol,PB3D,&
        &grid_eq_B,eq_limits,X_limits) result(ierr)
        use PB3D_vars, only: min_PB3D_version
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D'
        
        ! input  / output
        logical, intent(in) :: rec_eq, rec_X_1, rec_X_2, rec_sol                ! which quantities to reconstruct
        type(PB3D_type), intent(inout) :: PB3D                                  ! PB3D for which to do postprocessing
        type(grid_type), intent(inout), optional :: grid_eq_B                   ! optional field-aligned grid for HELENA
        integer, intent(in), optional :: eq_limits(2)                           ! i_limit of eq and X variables
        integer, intent(in), optional :: X_limits(2)                            ! i_limit of sol variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! local variables also used in child functions
        integer :: r_F_X_id, r_E_X_id                                           ! index of perturbation r_F and r_E
        integer :: r_F_eq_B_id, r_E_eq_B_id                                     ! index of field-aligned equilibrium r_F and r_E
        integer :: theta_F_B_id, zeta_F_B_id                                    ! index of field-aligned theta_F and zeta_F
        integer :: theta_E_B_id, zeta_E_B_id                                    ! index of field-aligned theta_E and zeta_E
        integer :: r_F_eq_id, r_E_eq_id                                         ! index of equilibrium r_F and r_E
        integer :: theta_F_id, zeta_F_id                                        ! index of theta_F and zeta_F
        integer :: theta_E_id, zeta_E_id                                        ! index of theta_E and zeta_E
        integer :: pres_FD_id, q_saf_FD_id, rot_t_FD_id                         ! index of pres_FD, q_saf_FD, rot_t_FD
        integer :: flux_p_FD_id, flux_t_FD_id                                   ! index of flux_p_FD, flux_t_FD
        integer :: pres_E_id, q_saf_E_id, rot_t_E_id                            ! index of pres_E, q_saf_E, rot_t_E
        integer :: flux_p_E_id, flux_t_E_id                                     ! index of flux_p_FD, flux_t_FD
        integer :: rho_id, S_id, kappa_n_id, kappa_g_id, sigma_id               ! index of rho, S, kappa_n, kappa_g, sigma
        integer :: R_V_c_id, R_V_s_id                                           ! index of R_V_c and R_V_s
        integer :: Z_V_c_id, Z_V_s_id                                           ! index of Z_V_c and Z_V_s
        integer :: L_V_c_id, L_V_s_id                                           ! index of L_V_c and L_V_s
        integer :: R_H_id, Z_H_id, chi_H_id, flux_p_H_id                        ! index of R_H, Z_H, chi_H and flux_p_H
        integer :: g_FD_id, h_FD_id, jac_FD_id                                  ! index of g_FD, h_FD, jac_FD
        integer :: misc_eq_id, misc_eq_V_id, misc_eq_H_id                       ! index of misc_eq, misc_eq_V, misc_eq_H
        integer :: misc_X_id                                                    ! index of misc_X
        integer :: RE_U_0_id, IM_U_0_id, RE_U_1_id, IM_U_1_id                   ! index of RE_U_0, IM_U_0, RE_U_1, IM_U_1
        integer :: RE_DU_0_id, IM_DU_0_id, RE_DU_1_id, IM_DU_1_id               ! index of RE_DU_0, IM_DU_0, RE_DU_1, IM_DU_1
        integer :: RE_X_val_id, IM_X_val_id                                     ! index of RE_X_val, IM_X_val
        integer :: RE_X_vec_id, IM_X_vec_id                                     ! index of RE_X_vec, IM_X_vec
        real(dp), allocatable :: dum_1D(:), dum_2D(:,:), dum_3D(:,:,:)          ! dummy variables
        real(dp), allocatable :: dum_4D(:,:,:,:)                                ! dummy variables
        !real(dp), allocatable :: dum_5D(:,:,:,:,:)                              ! dummy variables
        real(dp), allocatable :: dum_6D(:,:,:,:,:,:), dum_7D(:,:,:,:,:,:,:)     ! dummy variables
        real(dp), parameter :: tol_version = 1.E-8_dp                           ! tolerance for version control
        integer :: eq_limits_loc(2), X_limits_loc(2)                            ! local versions of eq_limits, X_limits
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Prepare variable indices')
        call lvl_ud(1)
        
        if (rec_eq) then
            ierr = prepare_vars_eq()
            CHCKERR('')
        end if
        
        if (rec_X_1) then
            ierr = prepare_vars_X_1()
            CHCKERR('')
        end if
        
        if (rec_X_2) then
            ierr = prepare_vars_X_2()
            CHCKERR('')
        end if
        
        if (rec_sol) then
            ierr = prepare_vars_sol()
            CHCKERR('')
        end if
        
        call lvl_ud(-1)
        call writo('Indices prepared')
        
        call writo('Reconstructing variables')
        call lvl_ud(1)
        
        call writo('Setting grids')
        call lvl_ud(1)
        
        if (rec_eq) then
            ierr = reconstruct_grid_eq()
            CHCKERR('')
        end if
        
        if (rec_sol) then
            ierr = reconstruct_grid_X()
            CHCKERR('')
        end if
        
        call lvl_ud(-1)
        
        if (rec_eq) then
            ierr = reconstruct_vars_eq()
            CHCKERR('')
        end if
        
        if (rec_X_1) then
            ierr = reconstruct_vars_X_1()
            CHCKERR('')
        end if
        
        if (rec_X_2) then
            ierr = reconstruct_vars_X_2()
            CHCKERR('')
        end if
        
        if (rec_sol) then
            call reconstruct_vars_sol()
        end if
        
        call lvl_ud(-1)
        call writo('Variables reconstructed')
        
        call writo('Running tests')
        call lvl_ud(1)
        
        call writo('PB3D version '//trim(r2strt(PB3D%version)))
        if (PB3D%version.lt.min_PB3D_version*(1-tol_version)) then
            ierr = 1
            err_msg = 'Need at least PB3D version '//&
                &trim(r2strt(min_PB3D_version))
            CHCKERR(err_msg)
        end if
        
        call lvl_ud(-1)
        call writo('Tests done')
    contains
        ! prepare eq vars
        integer function prepare_vars_eq() result(ierr)
            use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'prepare_vars_eq'
            
            ! initialize ierr
            ierr = 0
            
            ! set up 1D indices common for all equilibrium styles
            ierr = retrieve_var_1D_id(vars_1D_eq,'r_F',r_F_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'r_E',r_E_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'theta_F',theta_F_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'zeta_F',zeta_F_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'theta_E',theta_E_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'zeta_E',zeta_E_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'pres_FD',pres_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'q_saf_FD',q_saf_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'rot_t_FD',rot_t_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'flux_p_FD',flux_p_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'flux_t_FD',flux_t_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'rho',rho_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'S',S_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'kappa_n',kappa_n_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'kappa_g',kappa_g_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'sigma',sigma_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'g_FD',g_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'h_FD',h_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'jac_FD',jac_FD_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_eq,'misc_eq',misc_eq_id)
            CHCKERR('')
            
            ! Set up particular 1D eq indices, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = retrieve_var_1D_id(vars_1D_eq,'R_V_c',R_V_c_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'R_V_s',R_V_s_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'Z_V_c',Z_V_c_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'Z_V_s',Z_V_s_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'L_V_c',L_V_c_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'L_V_s',L_V_s_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'misc_eq_V',&
                        &misc_eq_V_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'pres_E',pres_E_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'q_saf_E',q_saf_E_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'rot_t_E',rot_t_E_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'flux_p_E',flux_p_E_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'flux_t_E',flux_t_E_id)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ierr = retrieve_var_1D_id(vars_1D_eq,'R_H',R_H_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'Z_H',Z_H_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'chi_H',chi_H_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'flux_p_H',flux_p_H_id)
                    CHCKERR('')
                    ierr = retrieve_var_1D_id(vars_1D_eq,'misc_eq_H',&
                        &misc_eq_H_id)
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! get 1D eq_B indices if present
            if (present(grid_eq_B)) then 
                ierr = retrieve_var_1D_id(vars_1D_eq_B,'r_F',r_F_eq_B_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq_B,'r_E',r_E_eq_B_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq_B,'theta_F',theta_F_B_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq_B,'zeta_F',zeta_F_B_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq_B,'theta_E',theta_E_B_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq_B,'zeta_E',zeta_E_B_id)
                CHCKERR('')
            end if
        end function prepare_vars_eq
        
        ! prepare vectorial perturbation vars
        integer function prepare_vars_X_1() result(ierr)
            use PB3D_vars, only: vars_1D_X_1
            
            character(*), parameter :: rout_name = 'prepare_vars_X_1'
            
            ! initialize ierr
            ierr = 0
            
            ! get 1D X_1 indices
            ierr = retrieve_var_1D_id(vars_1D_X_1,'misc_X',misc_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'RE_U_0',RE_U_0_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'IM_U_0',IM_U_0_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'RE_U_1',RE_U_1_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'IM_U_1',IM_U_1_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'RE_DU_0',RE_DU_0_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'IM_DU_0',IM_DU_0_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'RE_DU_1',RE_DU_1_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_X_1,'IM_DU_1',IM_DU_1_id)
            CHCKERR('')
        end function prepare_vars_X_1
        
        ! prepare tensorial perturbation vars
        integer function prepare_vars_X_2() result(ierr)
            use PB3D_vars, only: vars_1D_X_2
            
            character(*), parameter :: rout_name = 'prepare_vars_X_2'
            
            ! initialize ierr
            ierr = 0
            
            ! get 1D X_2 indices
            write(*,*) '¡¡¡ NOT YET IMPLEMENTED READING OF X_2 !!!!'
        end function prepare_vars_X_2
        
        ! prepare solution vars
        integer function prepare_vars_sol() result(ierr)
            use PB3D_vars, only: vars_1D_sol
            
            character(*), parameter :: rout_name = 'prepare_vars_sol'
            
            ! initialize ierr
            ierr = 0
            
            ! get 1D sol indices
            ierr = retrieve_var_1D_id(vars_1D_sol,'r_F',r_F_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'r_E',r_E_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'RE_X_val',RE_X_val_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'IM_X_val',IM_X_val_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'RE_X_vec',RE_X_vec_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'IM_X_vec',IM_X_vec_id)
            CHCKERR('')
        end function prepare_vars_sol
        
        ! reconstruct equilibrium grid
        integer function reconstruct_grid_eq() result(ierr)
            use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B
            use grid_vars, only: create_grid
            
            character(*), parameter :: rout_name = 'reconstruct_grid_eq'
            
            ! local variables
            integer :: n_eq(3), n_eq_B(3)                                       ! n of eq and eq_B grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Creating equilibrium grid')
            call lvl_ud(1)
            
            ! set n_eq
            n_eq = vars_1D_eq(theta_F_id)%tot_i_max-&
                &vars_1D_eq(theta_F_id)%tot_i_min+1
            
            ! set up local eq_limits
            eq_limits_loc = [1,n_eq(3)]
            if (present(eq_limits)) eq_limits_loc = eq_limits
            
            ! set grid
            ierr = create_grid(PB3D%grid_eq,n_eq,eq_limits_loc)
            CHCKERR('')
            PB3D%grid_eq%r_F = vars_1D_eq(r_F_eq_id)%p
            PB3D%grid_eq%loc_r_F = &
                &vars_1D_eq(r_F_eq_id)%p(eq_limits_loc(1):eq_limits_loc(2))
            PB3D%grid_eq%r_E = vars_1D_eq(r_E_eq_id)%p
            PB3D%grid_eq%loc_r_E = &
                &vars_1D_eq(r_E_eq_id)%p(eq_limits_loc(1):eq_limits_loc(2))
            call conv_1D2ND(vars_1D_eq(theta_F_id),dum_3D)
            PB3D%grid_eq%theta_F = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(zeta_F_id),dum_3D)
            PB3D%grid_eq%zeta_F = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(theta_E_id),dum_3D)
            PB3D%grid_eq%theta_E = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(zeta_E_id),dum_3D)
            PB3D%grid_eq%zeta_E = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call writo('angular size: ('//trim(i2str(n_eq(1)))//','//&
                &trim(i2str(n_eq(2)))//')')
            call writo('normal size: '//trim(i2str(n_eq(3))))
            call lvl_ud(-1)
            
            if (present(grid_eq_B)) then
                call writo('Creating field-aligned equilibrium grid')
                call lvl_ud(1)
                n_eq_B = vars_1D_eq_B(theta_F_B_id)%tot_i_max-&
                    &vars_1D_eq_B(theta_F_B_id)%tot_i_min+1
                ierr = create_grid(grid_eq_B,n_eq_B,eq_limits_loc)
                CHCKERR('')
                grid_eq_B%r_F = vars_1D_eq_B(r_F_eq_B_id)%p
                grid_eq_B%loc_r_F = vars_1D_eq_B(r_F_eq_B_id)%&
                    &p(eq_limits_loc(1):eq_limits_loc(2))
                grid_eq_B%r_E = vars_1D_eq_B(r_E_eq_B_id)%p
                grid_eq_B%loc_r_E = vars_1D_eq_B(r_E_eq_B_id)%&
                    &p(eq_limits_loc(1):eq_limits_loc(2))
                call conv_1D2ND(vars_1D_eq_B(theta_F_B_id),dum_3D)
                grid_eq_B%theta_F = &
                    &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
                call conv_1D2ND(vars_1D_eq_B(zeta_F_B_id),dum_3D)
                grid_eq_B%zeta_F = &
                    &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
                call conv_1D2ND(vars_1D_eq_B(theta_E_B_id),dum_3D)
                grid_eq_B%theta_E = &
                    &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
                call conv_1D2ND(vars_1D_eq_B(zeta_E_B_id),dum_3D)
                grid_eq_B%zeta_E = &
                    &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
                call writo('angular size: ('//trim(i2str(n_eq_B(1)))//','//&
                    &trim(i2str(n_eq_B(2)))//')')
                call writo('normal size: '//trim(i2str(n_eq_B(3))))
                call lvl_ud(-1)
            end if
        end function reconstruct_grid_eq
        
        ! reconstruct perturbation grid
        integer function reconstruct_grid_X() result(ierr)
            use PB3D_vars, only: vars_1D_sol
            use grid_vars, only: create_grid
            
            character(*), parameter :: rout_name = 'reconstruct_grid_X'
            
            ! local variables
            integer :: n_X                                                      ! n of X grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Creating perturbation grid')
            call lvl_ud(1)
            
            ! set n_X
            n_X = vars_1D_sol(r_F_X_id)%tot_i_max(1)-&
                &vars_1D_sol(r_F_X_id)%tot_i_min(1)+1
            
            ! set up local X_limits
            X_limits_loc = [1,n_X]
            if (present(X_limits)) X_limits_loc = X_limits
            
            ierr = create_grid(PB3D%grid_X,n_X,X_limits_loc)
            CHCKERR('')
            PB3D%grid_X%r_F = vars_1D_sol(r_F_X_id)%p
            PB3D%grid_X%loc_r_F = &
                &vars_1D_sol(r_F_X_id)%p(X_limits_loc(1):X_limits_loc(2))
            PB3D%grid_X%r_E = vars_1D_sol(r_E_X_id)%p
            PB3D%grid_X%loc_r_E = &
                &vars_1D_sol(r_E_X_id)%p(X_limits_loc(1):X_limits_loc(2))
            call writo('normal size: '//trim(i2str(n_X)))
            call lvl_ud(-1)
        end function reconstruct_grid_X
        
        ! reconstruct equilibrium vars
        integer function reconstruct_vars_eq() result(ierr)
            use num_vars, only: eq_style, max_deriv, use_pol_flux_E, &
                &use_pol_flux_F, use_normalization, norm_disc_prec_eq
            use met_vars, only: create_met
            use PB3D_vars, only: vars_1D_eq
            use eq_vars, only: create_eq, R_0, pres_0, B_0, psi_0, rho_0, &
                &T_0, vac_perm, max_flux_p_E, max_flux_t_E, max_flux_p_F, &
                &max_flux_t_F
            use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, &
                &ntor, nfp, lfreeB, lasym
            use HELENA, only: R_H, Z_H, chi_H, flux_p_H, ias, nchi
            
            character(*), parameter :: rout_name = 'reconstruct_vars_eq'
            
            ! local variables
            integer :: n_r_eq                                                   ! n_r of eq grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Setting equilibrium')
            ierr = create_eq(PB3D%grid_eq,PB3D%eq)
            CHCKERR('')
            
            call conv_1D2ND(vars_1D_eq(pres_FD_id),dum_2D)
            PB3D%eq%pres_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(q_saf_FD_id),dum_2D)
            PB3D%eq%q_saf_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(rot_t_FD_id),dum_2D)
            PB3D%eq%rot_t_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(flux_p_FD_id),dum_2D)
            PB3D%eq%flux_p_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(flux_t_FD_id),dum_2D)
            PB3D%eq%flux_t_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(rho_id),dum_1D)
            PB3D%eq%rho = dum_1D(eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_eq(S_id),dum_3D)
            PB3D%eq%S = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(kappa_n_id),dum_3D)
            PB3D%eq%kappa_n = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(kappa_g_id),dum_3D)
            PB3D%eq%kappa_g = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(sigma_id),dum_3D)
            PB3D%eq%sigma = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(g_FD_id),dum_7D)
            
            ! Set up particular eq variables, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    call conv_1D2ND(vars_1D_eq(misc_eq_V_id),dum_1D)
                    lasym = .false.
                    lfreeB = .false.
                    if (dum_1D(1).gt.0) lasym = .true.
                    if (dum_1D(2).gt.0) lfreeB = .true.
                    mpol = nint(dum_1D(3))
                    ntor = nint(dum_1D(4))
                    nfp = nint(dum_1D(5))
                    n_r_eq = vars_1D_eq(R_V_c_id)%tot_i_max(3)-&
                        &vars_1D_eq(R_V_c_id)%tot_i_max(3)+1
                    deallocate(dum_1D)
                    
                    allocate(R_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
                    call conv_1D2ND(vars_1D_eq(R_V_c_id),dum_4D)
                    R_V_c = dum_4D
                    deallocate(dum_4D)
                    allocate(R_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
                    call conv_1D2ND(vars_1D_eq(R_V_s_id),dum_4D)
                    R_V_s = dum_4D
                    deallocate(dum_4D)
                    allocate(Z_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
                    call conv_1D2ND(vars_1D_eq(Z_V_c_id),dum_4D)
                    Z_V_c = dum_4D
                    deallocate(dum_4D)
                    allocate(Z_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
                    call conv_1D2ND(vars_1D_eq(Z_V_s_id),dum_4D)
                    Z_V_s = dum_4D
                    deallocate(dum_4D)
                    allocate(L_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
                    call conv_1D2ND(vars_1D_eq(L_V_c_id),dum_4D)
                    L_V_c = dum_4D
                    deallocate(dum_4D)
                    allocate(L_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
                    call conv_1D2ND(vars_1D_eq(L_V_s_id),dum_4D)
                    L_V_s = dum_4D
                    deallocate(dum_4D)
                    call conv_1D2ND(vars_1D_eq(pres_E_id),dum_2D)
                    PB3D%eq%pres_E = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(q_saf_E_id),dum_2D)
                    PB3D%eq%q_saf_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(rot_t_E_id),dum_2D)
                    PB3D%eq%rot_t_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(flux_p_E_id),dum_2D)
                    PB3D%eq%flux_p_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(flux_t_E_id),dum_2D)
                    PB3D%eq%flux_t_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                case (2)                                                        ! HELENA
                    call conv_1D2ND(vars_1D_eq(misc_eq_H_id),dum_1D)
                    ias = nint(dum_1D(1))
                    nchi = nint(dum_1D(2))
                    n_r_eq = vars_1D_eq(R_H_id)%tot_i_max(2)-&
                        &vars_1D_eq(R_H_id)%tot_i_max(1)+1
                    deallocate(dum_1D)
                    
                    allocate(R_H(1:nchi,1:n_r_eq))
                    call conv_1D2ND(vars_1D_eq(R_H_id),dum_2D)
                    R_H = dum_2D
                    deallocate(dum_2D)
                    allocate(Z_H(1:nchi,1:n_r_eq))
                    call conv_1D2ND(vars_1D_eq(Z_H_id),dum_2D)
                    Z_H = dum_2D
                    deallocate(dum_2D)
                    allocate(flux_p_H(1:n_r_eq))
                    call conv_1D2ND(vars_1D_eq(flux_p_H_id),dum_1D)
                    flux_p_H = dum_1D
                    deallocate(dum_1D)
                    allocate(chi_H(1:nchi))
                    call conv_1D2ND(vars_1D_eq(chi_H_id),dum_1D)
                    chi_H = dum_1D
                    deallocate(dum_1D)
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call writo('Setting metrics')
            ierr = create_met(PB3D%grid_eq,PB3D%met)
            CHCKERR('')
            
            PB3D%met%g_FD = &
                        &dum_7D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:,:)
            deallocate(dum_7D)
            call conv_1D2ND(vars_1D_eq(h_FD_id),dum_7D)
            PB3D%met%h_FD = &
                        &dum_7D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:,:)
            deallocate(dum_7D)
            call conv_1D2ND(vars_1D_eq(jac_FD_id),dum_6D)
            PB3D%met%jac_FD = &
                        &dum_6D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:)
            deallocate(dum_6D)
            call conv_1D2ND(vars_1D_eq(misc_eq_id),dum_1D)
            PB3D%version = dum_1D(1)
            PB3D%alpha = dum_1D(4)
            R_0 = dum_1D(5)
            pres_0 = dum_1D(6)
            B_0 = dum_1D(7)
            psi_0 = dum_1D(8)
            rho_0 = dum_1D(9)
            T_0 = dum_1D(10)
            vac_perm = dum_1D(11)
            max_flux_p_E = dum_1D(12)
            max_flux_t_E = dum_1D(13)
            max_flux_p_F = dum_1D(14)
            max_flux_t_F = dum_1D(15)
            use_pol_flux_E = .false.
            use_pol_flux_F = .false.
            use_normalization = .false.
            if (dum_1D(16).gt.0) use_pol_flux_E = .true.
            if (dum_1D(17).gt.0) use_pol_flux_F = .true.
            if (dum_1D(18).gt.0) use_normalization = .true.
            norm_disc_prec_eq = nint(dum_1D(19))
            deallocate(dum_1D)
        end function reconstruct_vars_eq
        
        ! reconstruct vectorial perturbation vars
        integer function reconstruct_vars_X_1() result(ierr)
            use num_vars, only: norm_disc_prec_X, use_pol_flux_F
            use PB3D_vars, only: vars_1D_X_1
            use X_vars, only: create_X, &
                &min_r_X, max_r_X, min_n_X, max_n_X, min_m_X, max_m_X
            
            character(*), parameter :: rout_name = 'reconstruct_vars_X_1'
            
            ! initialize ierr
            ierr = 0
            
            call writo('Setting vectorial perturbation')
            call lvl_ud(1)
            
            call conv_1D2ND(vars_1D_X_1(misc_X_id),dum_1D)
            min_r_X = dum_1D(1)
            max_r_X = dum_1D(2)
            min_n_X = nint(dum_1D(3))
            max_n_X = nint(dum_1D(4))
            min_m_X = nint(dum_1D(5))
            max_m_X = nint(dum_1D(6))
            norm_disc_prec_X = nint(dum_1D(7))
            deallocate(dum_1D)
            
            ! user output
            call writo('The PB3D output')
            call lvl_ud(1)
            if (use_pol_flux_F) then
                call writo('has modes n = '//trim(i2str(min_n_X))//&
                    &' and m = '//trim(i2str(min_m_X))//'..'//&
                    &trim(i2str(max_m_X)))
                call writo('works with the poloidal flux as normal coordinate')
            else
                call writo('has modes m = '//trim(i2str(min_m_X))//&
                    &' and n='//trim(i2str(min_n_X))//'..'//&
                    &trim(i2str(max_n_X)))
                call writo('works with the toroidal flux as normal coordinate')
            end if
            call writo('and its normal boundaries are '//&
            &trim(r2strt(min_r_X))//' and '//trim(r2strt(max_r_X)))
            call lvl_ud(-1)
            
            ierr = create_X(PB3D%grid_eq,PB3D%X_1)
            CHCKERR('')
            
            ! set up local eq_limits
            eq_limits_loc = [1,PB3D%grid_eq%n(3)]
            if (present(eq_limits)) eq_limits_loc = eq_limits
            
            call conv_1D2ND(vars_1D_X_1(RE_U_0_id),dum_4D)
            PB3D%X_1%U_0 = dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(IM_U_0_id),dum_4D)
            PB3D%X_1%U_0 = PB3D%X_1%U_0 + &
                &iu*dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(RE_U_1_id),dum_4D)
            PB3D%X_1%U_1 = dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(IM_U_1_id),dum_4D)
            PB3D%X_1%U_1 = PB3D%X_1%U_1 + &
                &iu*dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(RE_DU_0_id),dum_4D)
            PB3D%X_1%DU_0 = dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(IM_DU_0_id),dum_4D)
            PB3D%X_1%DU_0 = PB3D%X_1%DU_0 + &
                &iu*dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(RE_DU_1_id),dum_4D)
            PB3D%X_1%DU_1 = dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            call conv_1D2ND(vars_1D_X_1(IM_DU_1_id),dum_4D)
            PB3D%X_1%DU_1 = PB3D%X_1%DU_1 + &
                &iu*dum_4D(:,:,eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_4D)
            
            call lvl_ud(-1)
        end function reconstruct_vars_X_1
        
        ! reconstruct tensorial perturbation vars
        integer function reconstruct_vars_X_2() result(ierr)
            use PB3D_vars, only: vars_1D_X_2
            
            character(*), parameter :: rout_name = 'reconstruct_vars_X_2'
            
            ! initialize ierr
            ierr = 0
            
            call writo('Setting tensorial perturbation')
            call lvl_ud(1)
            write(*,*) '¡¡¡ NOT YET IMPLEMENTED RECONSTRUCTING X_2 !!!!'
            
            call lvl_ud(-1)
        end function reconstruct_vars_X_2
        
        ! reconstruct solution vars
        subroutine reconstruct_vars_sol()
            use PB3D_vars, only: vars_1D_sol
            
            call writo('Setting solution')
            call lvl_ud(1)
            
            ! set up local sol_limits
            X_limits_loc = [1,PB3D%grid_X%n(3)]
            if (present(X_limits)) X_limits_loc = X_limits
            
            allocate(PB3D%sol%val(vars_1D_sol(RE_X_val_id)%tot_i_min(1):&
                &vars_1D_sol(RE_X_val_id)%tot_i_max(1)))
            call conv_1D2ND(vars_1D_sol(RE_X_val_id),dum_1D)
            PB3D%sol%val = dum_1D
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_sol(IM_X_val_id),dum_1D)
            PB3D%sol%val = PB3D%sol%val + iu*dum_1D
            deallocate(dum_1D)
            allocate(PB3D%sol%vec(vars_1D_sol(RE_X_vec_id)%tot_i_min(1):&
                &vars_1D_sol(RE_X_vec_id)%tot_i_max(1),&
                &1:X_limits_loc(2)-X_limits_loc(1)+1,&
                &vars_1D_sol(RE_X_vec_id)%tot_i_min(3):&
                &vars_1D_sol(RE_X_vec_id)%tot_i_max(3)))
            call conv_1D2ND(vars_1D_sol(RE_X_vec_id),dum_3D)
            PB3D%sol%vec = dum_3D(:,X_limits_loc(1):X_limits_loc(2),:)
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_sol(IM_X_vec_id),dum_3D)
            PB3D%sol%vec = PB3D%sol%vec + &
                &iu*dum_3D(:,X_limits_loc(1):X_limits_loc(2),:)
            deallocate(dum_3D)
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_sol
    end function reconstruct_PB3D
    
    ! Retrieves variable index from array 1D equivalents
    integer function retrieve_var_1D_id(vars,var_name,var_id) &
        &result(ierr)
        character(*), parameter :: rout_name = 'retrieve_var_1D'
        
        ! input / output
        type(var_1D), intent(in) :: vars(:)                                     ! array of 1D variables
        character(len=*), intent(in) :: var_name                                ! name of variable to retrieve
        integer, intent(inout) :: var_id                                        ! index of variable
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! look up the variable
        do id = 1,size(vars)
            if (trim(vars(id)%var_name).eq.trim(var_name)) then                 ! found
                var_id = id
                return
            end if
        end do
        
        ! if still here, nothing found
        ierr = 1
        err_msg = 'Variable '//trim(var_name)//' not found'
        CHCKERR(err_msg)
    end function retrieve_var_1D_id
    
    ! Converts 1D to nD variables. The output variable has to be allocatable and
    ! unallocated
    subroutine conv_1D2ND_1D(var_in,var_out)                                    ! 1D version
        ! input / output
        type(var_1D), intent(in) :: var_in                                      ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:)                      ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1])
    end subroutine conv_1D2ND_1D
    subroutine conv_1D2ND_2D(var_in,var_out)                                    ! 2D version
        ! input / output
        type(var_1D), intent(in) :: var_in                                      ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:)                    ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1])
    end subroutine conv_1D2ND_2D
    subroutine conv_1D2ND_3D(var_in,var_out)                                    ! 3D version
        ! input / output
        type(var_1D), intent(in) :: var_in                                      ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:)                  ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1])
    end subroutine conv_1D2ND_3D
    subroutine conv_1D2ND_4D(var_in,var_out)                                    ! 4D version
        ! input / output
        type(var_1D), intent(in) :: var_in                                      ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:)                ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            &var_in%tot_i_max(4)-var_in%tot_i_min(4)+1])
    end subroutine conv_1D2ND_4D
    !subroutine conv_1D2ND_5D(var_in,var_out)                                    ! 5D version
        !! input / output
        !type(var_1D), intent(in) :: var_in                                      ! 1D variable
        !real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:)              ! output variable
        
        !! allocate and copy variable
        !allocate(var_out(&
            !&var_in%tot_i_min(1):var_in%tot_i_max(1),&
            !&var_in%tot_i_min(2):var_in%tot_i_max(2),&
            !&var_in%tot_i_min(3):var_in%tot_i_max(3),&
            !&var_in%tot_i_min(4):var_in%tot_i_max(4),&
            !&var_in%tot_i_min(5):var_in%tot_i_max(5)))
        !var_out = reshape(var_in%p,[&
            !&var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            !&var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            !&var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            !&var_in%tot_i_max(4)-var_in%tot_i_min(4)+1,&
            !&var_in%tot_i_max(5)-var_in%tot_i_min(5)+1])
    !end subroutine conv_1D2ND_5D
    subroutine conv_1D2ND_6D(var_in,var_out)                                    ! 6D version
        ! input / output
        type(var_1D), intent(in) :: var_in                                      ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:)            ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4),&
            &var_in%tot_i_min(5):var_in%tot_i_max(5),&
            &var_in%tot_i_min(6):var_in%tot_i_max(6)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            &var_in%tot_i_max(4)-var_in%tot_i_min(4)+1,&
            &var_in%tot_i_max(5)-var_in%tot_i_min(5)+1,&
            &var_in%tot_i_max(6)-var_in%tot_i_min(6)+1])
    end subroutine conv_1D2ND_6D
    subroutine conv_1D2ND_7D(var_in,var_out)                                    ! 7D version
        ! input / output
        type(var_1D), intent(in) :: var_in                                      ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:,:)          ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4),&
            &var_in%tot_i_min(5):var_in%tot_i_max(5),&
            &var_in%tot_i_min(6):var_in%tot_i_max(6),&
            &var_in%tot_i_min(7):var_in%tot_i_max(7)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            &var_in%tot_i_max(4)-var_in%tot_i_min(4)+1,&
            &var_in%tot_i_max(5)-var_in%tot_i_min(5)+1,&
            &var_in%tot_i_max(6)-var_in%tot_i_min(6)+1,&
            &var_in%tot_i_max(7)-var_in%tot_i_min(7)+1])
    end subroutine conv_1D2ND_7D
end module PB3D_ops
