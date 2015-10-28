!------------------------------------------------------------------------------!
!   Operations on PB3D output.                                                 !
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
    use HDF5_vars, only: var_1D_type

    implicit none
    private
    public init_PB3D_ops, read_PB3D, reconstruct_PB3D, retrieve_var_1D_id
    
    ! interfaces
    interface conv_1D2ND
        module procedure conv_1D2ND_1D, conv_1D2ND_2D, conv_1D2ND_3D, &
            &conv_1D2ND_4D, &
            !&conv_1D2ND_5D, &
            &conv_1D2ND_6D, conv_1D2ND_7D
    end interface
    interface get_full_var_names
        module procedure get_full_var_names_1, get_full_var_names_2
    end interface
    
contains
    ! Initialize the PB3D routines:
    !   - set the equilibrium style
    integer function init_PB3D_ops() result(ierr)
        use num_vars, only: eq_style, PB3D_name, rank
        use PB3D_vars, only: vars_1D_misc
        use HDF5_ops, only: read_HDF5_arrs
        
        character(*), parameter :: rout_name = 'init_PB3D_ops'
        
        ! local variables
        integer :: eq_misc_id                                                   ! indices of miscellaneous variables
        real(dp), allocatable :: dum_1D(:)                                      ! dummy variable
        
        ! initialize ierr
        ierr = 0
        
        if (rank.eq.0) then
            ! user output
            call writo('Initializing PB3D operations')
            call lvl_ud(1)
            
            ! read miscellaneous variables
            ierr = read_HDF5_arrs(vars_1D_misc,PB3D_name,'misc')
            CHCKERR('')
            
            ! get equilibrium style
            ierr = retrieve_var_1D_id(vars_1D_misc,'eq',eq_misc_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D_misc(eq_misc_id),dum_1D)
            eq_style = nint(dum_1D(2))
            
            ! clean up
            deallocate(vars_1D_misc)
            
            ! user output
            call lvl_ud(-1)
            call writo('PB3D operations Initialized')
        end if
    end function init_PB3D_ops
    
    ! Reads PB3D output from user-provided input file.
    ! Optionally, if perturbation  quantities are read, the limits  of m_X (pol.
    ! flux) or  n_X (tor.  flux) can be  passed, instead of  reading all  of the
    ! perturbation quantities.
    integer function read_PB3D(read_misc,read_eq,read_X_1,read_X_2,read_sol,&
        &lim_sec_X_1,lim_sec_X_2) result(ierr)
        use num_vars, only: eq_style, PB3D_name
        use HDF5_ops, only: read_HDF5_arrs
        use X_vars, only: get_suffix, &
            &X_1_var_names, X_2_var_names
        use PB3D_vars, only: vars_1D_misc, vars_1D_eq, vars_1D_eq_B, &
            &vars_1D_X_1, vars_1D_X_2, vars_1D_sol
        
        character(*), parameter :: rout_name = 'read_PB3D'
        
        ! input / output
        logical, intent(in) :: read_misc, read_eq, read_X_1, read_X_2, read_sol ! which files to read
        integer, intent(in), optional :: lim_sec_X_1(2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in), optional :: lim_sec_X_2(2,2)                       ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln), allocatable :: req_var_names(:)              ! requested variable names
        
        ! initialize ierr
        ierr = 0
        
        ! set common error message
        err_msg = 'Maybe you wanted to use Richardson extrapolation? Not &
            &yet implemented.'
        
        if (read_misc) then
            call writo('Reading miscellaneous variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_misc,PB3D_name,'misc')
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Miscellaneous variables read')
        end if
        
        if (read_eq) then
            call writo('Reading equilibrium variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_eq,PB3D_name,'eq')
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
            if (present(lim_sec_X_1)) then
                ! get requested full variable names
                req_var_names = get_full_var_names(X_1_var_names,lim_sec_X_1)
                
                ! read from HDF5 file
                ierr = read_HDF5_arrs(vars_1D_X_1,PB3D_name,'X_1',&
                    &acc_var_names=req_var_names)
                CHCKERR(err_msg)
            else
                ierr = read_HDF5_arrs(vars_1D_X_1,PB3D_name,'X_1')
                CHCKERR(err_msg)
            end if
            call lvl_ud(-1)
            call writo('Vectorial perturbation variables read')
        end if
        
        if (read_X_2) then
            call writo('Reading tensorial perturbation variables')
            call lvl_ud(1)
            if (present(lim_sec_X_2)) then
                ! get requested full variable names
                req_var_names = get_full_var_names(X_2_var_names,&
                    &[.true.,.true.,.false.,.false.,.true.,.true.,&
                    &.true.,.true.,.false.,.false.,.true.,.true.],lim_sec_X_2)
                
                ! read from HDF5 file
                ierr = read_HDF5_arrs(vars_1D_X_2,PB3D_name,'X_2',&
                    &acc_var_names=req_var_names)
                CHCKERR(err_msg)
            else
                ierr = read_HDF5_arrs(vars_1D_X_2,PB3D_name,'X_2')
                CHCKERR(err_msg)
            end if
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
    end function read_PB3D
    
    ! Reconstructs the PB3D variables: eq,  X and/or sol variables, as indicated
    ! through rec_eq, rec_X_1, rec_X_2 and rec_sol.
    ! Additionally,  if the  variables 'eq_limits'  and/or 'sol_limits'  are not
    ! provided, the full normal range is taken for every variable.
    ! Note that both X variables need the equilibrium grid variables.
    ! Optionally, if perturbation  quantities are read, the limits  of m_X (pol.
    ! flux) or  n_X (tor.  flux) can be  passed, instead of  reading all  of the
    ! perturbation quantities.
    integer function reconstruct_PB3D(rec_misc,rec_eq,rec_X_1,rec_X_2,rec_sol,&
        &grid_eq,grid_eq_B,grid_sol,eq,met,X_1,X_2,sol,eq_limits,sol_limits,&
        &lim_sec_X_1,lim_sec_X_2) result(ierr)
        use HDF5_vars, only: dealloc_var_1D
        use PB3D_vars, only: vars_1D_misc, vars_1D_eq, vars_1D_eq_B, &
            &vars_1D_X_1, vars_1D_X_2, vars_1D_sol, min_PB3D_version, &
            &PB3D_version
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D'
        
        ! input  / output
        logical, intent(in) :: rec_misc, rec_eq, rec_X_1, rec_X_2, rec_sol      ! which quantities to reconstruct
        type(grid_type), intent(inout), optional :: grid_eq                     ! equilibrium grid 
        type(grid_type), intent(inout), optional :: grid_eq_B                   ! optional field-aligned grid for HELENA
        type(grid_type), intent(inout), optional :: grid_sol                    ! solution grid 
        type(eq_type), intent(inout), optional :: eq                            ! equilibrium variables
        type(met_type), intent(inout), optional :: met                          ! metric variables
        type(X_1_type), intent(inout), optional :: X_1                          ! vectorial perturbation variables
        type(X_2_type), intent(inout), optional :: X_2                          ! tensorial perturbation variables
        type(sol_type), intent(inout), optional :: sol                          ! solution variables
        integer, intent(in), optional :: eq_limits(2)                           ! i_limit of eq and X variables
        integer, intent(in), optional :: sol_limits(2)                          ! i_limit of sol variables
        integer, intent(in), optional :: lim_sec_X_1(2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in), optional :: lim_sec_X_2(2,2)                       ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! local variables also used in child functions
        integer :: eq_misc_id, X_misc_id, eq_V_misc_id, eq_H_misc_id            ! indices of miscellaneous variables
        integer :: r_F_X_id, r_E_X_id                                           ! index of perturbation r_F and r_E
        integer :: r_F_eq_B_id, r_E_eq_B_id                                     ! index of field-aligned equilibrium r_F and r_E
        integer :: theta_F_B_id, zeta_F_B_id                                    ! index of field-aligned theta_F and zeta_F
        integer :: theta_E_B_id, zeta_E_B_id                                    ! index of field-aligned theta_E and zeta_E
        integer :: r_F_eq_id, r_E_eq_id                                         ! index of equilibrium r_F and r_E
        integer :: theta_F_id, zeta_F_id                                        ! index of theta_F and zeta_F
        integer :: theta_E_id, zeta_E_id                                        ! index of theta_E and zeta_E
        integer :: max_flux_id                                                  ! index of max_flux
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
        integer, allocatable :: RE_U_0_id(:), IM_U_0_id(:)                      ! index of RE_U_0, IM_U_0, RE_U_1, IM_U_1
        integer, allocatable :: RE_U_1_id(:), IM_U_1_id(:)                      ! index of RE_U_1, IM_U_1, RE_U_1, IM_U_1
        integer, allocatable :: RE_DU_0_id(:), IM_DU_0_id(:)                    ! index of RE_DU_0, IM_DU_0, RE_DU_1, IM_DU_1
        integer, allocatable :: RE_DU_1_id(:), IM_DU_1_id(:)                    ! index of RE_DU_1, IM_DU_1, RE_DU_1, IM_DU_1
        integer, allocatable :: RE_PV_int_0_id(:), IM_PV_int_0_id(:)            ! index of RE_PV_int_0, IM_PV_int_0
        integer, allocatable :: RE_PV_int_1_id(:), IM_PV_int_1_id(:)            ! index of RE_PV_int_1, IM_PV_int_1
        integer, allocatable :: RE_PV_int_2_id(:), IM_PV_int_2_id(:)            ! index of RE_PV_int_2, IM_PV_int_2
        integer, allocatable :: RE_KV_int_0_id(:), IM_KV_int_0_id(:)            ! index of RE_KV_int_0, IM_KV_int_0
        integer, allocatable :: RE_KV_int_1_id(:), IM_KV_int_1_id(:)            ! index of RE_KV_int_1, IM_KV_int_1
        integer, allocatable :: RE_KV_int_2_id(:), IM_KV_int_2_id(:)            ! index of RE_KV_int_2, IM_KV_int_2
        integer :: RE_X_val_id, IM_X_val_id                                     ! index of RE_X_val, IM_X_val
        integer :: RE_X_vec_id, IM_X_vec_id                                     ! index of RE_X_vec, IM_X_vec
        real(dp), allocatable :: dum_1D(:), dum_2D(:,:), dum_3D(:,:,:)          ! dummy variables
        real(dp), allocatable :: dum_4D(:,:,:,:)                                ! dummy variables
        !real(dp), allocatable :: dum_5D(:,:,:,:,:)                              ! dummy variables
        real(dp), allocatable :: dum_6D(:,:,:,:,:,:), dum_7D(:,:,:,:,:,:,:)     ! dummy variables
        real(dp), parameter :: tol_version = 1.E-8_dp                           ! tolerance for version control
        integer :: eq_limits_loc(2), sol_limits_loc(2)                          ! local versions of eq_limits, sol_limits
        
        ! initialize ierr
        ierr = 0
        
        ! test whether appropriate variables provided
        if (rec_eq .and. .not.&
            &(present(grid_eq).and.present(eq).and.present(met))) then
            ierr = 1
            err_msg = 'To reconstruct equilibrium need grid_eq, eq and met'
            CHCKERR(err_msg)
        end if
        if (rec_X_1 .and. .not.(present(grid_eq).and.present(X_1))) then
            ierr = 1
            err_msg = 'To reconstruct equilibrium need grid_eq and X_1'
            CHCKERR(err_msg)
        end if
        if (rec_X_2 .and. .not.(present(grid_eq).and.present(X_2))) then
            ierr = 1
            err_msg = 'To reconstruct equilibrium need grid_eq and X_2'
            CHCKERR(err_msg)
        end if
        if (rec_sol .and. .not.(present(grid_sol).and.present(sol))) then
            ierr = 1
            err_msg = 'To reconstruct equilibrium need grid_sol and sol'
            CHCKERR(err_msg)
        end if
        
        if (rec_misc) then
            ! user output
            call writo('Prepare miscellaneous variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_misc()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing miscellaneous variables')
            call lvl_ud(1)
            
            call writo('Setting variables')
            call lvl_ud(1)
            ierr = reconstruct_vars_misc()
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_misc)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Miscellaneous variables reconstructed')
            
            call writo('Running tests')
            call lvl_ud(1)
            
            call writo('PB3D version '//trim(r2strt(PB3D_version)))
            if (PB3D_version.lt.min_PB3D_version*(1-tol_version)) then
                ierr = 1
                err_msg = 'Need at least PB3D version '//&
                    &trim(r2strt(min_PB3D_version))
                CHCKERR(err_msg)
            end if
            
            call lvl_ud(-1)
            call writo('Tests done')
        end if
        
        if (rec_eq) then
            ! user output
            call writo('Prepare equilibrium variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_eq()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing equilibrium variables')
            call lvl_ud(1)
            
            call writo('Setting grids')
            call lvl_ud(1)
            ierr = reconstruct_grid_eq(grid_eq,grid_eq_B)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Setting variables')
            call lvl_ud(1)
            ierr = reconstruct_vars_eq(grid_eq,eq,met)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_eq)
            if (present(grid_eq_B)) then 
                call dealloc_var_1D(vars_1D_eq_B)
            end if
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Equilibrium variables reconstructed')
        end if
        
        if (rec_X_1) then
            ! user output
            call writo('Prepare vectorial perturbation variable indices')
            call lvl_ud(1)
            
            ierr = prepare_vars_X_1(lim_sec_X_1)
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing vectorial perturbation variables')
            call lvl_ud(1)
            
            call writo('Setting variables')
            call lvl_ud(1)
            call reconstruct_vars_X_1(grid_eq,X_1,lim_sec_X_1)
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_X_1)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Vectorial perturbation variables reconstructed')
        end if
        
        if (rec_X_2) then
            ! user output
            call writo('Prepare tensorial perturbation variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_X_2(lim_sec_X_2)
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing tensorial perturbation variables')
            call lvl_ud(1)
            
            call writo('Setting variables')
            call lvl_ud(1)
            call reconstruct_vars_X_2(grid_eq,X_2,lim_sec_X_2)
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_X_2)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('tensorial perturbation variables reconstructed')
        end if
        
        if (rec_sol) then
            ! user output
            call writo('Prepare solution variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_sol()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing solution variables')
            call lvl_ud(1)
            
            call writo('Setting grids')
            call lvl_ud(1)
            ierr = reconstruct_grid_sol(grid_sol)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Setting variables')
            call lvl_ud(1)
            call reconstruct_vars_sol(grid_sol,sol)
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_sol)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Solution variables reconstructed')
        end if
    contains
        ! prepare miscellaneous vars
        integer function prepare_vars_misc() result(ierr)
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'prepare_vars_misc'
            
            ! initialize ierr
            ierr = 0
            
            ! set up 1D indices common for all equilibrium styles
            ierr = retrieve_var_1D_id(vars_1D_misc,'eq',eq_misc_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_misc,'X',X_misc_id)
            CHCKERR('')
            
            ! Set up particular 1D eq indices, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = retrieve_var_1D_id(vars_1D_misc,'eq_V',eq_V_misc_id)
                case (2)                                                        ! HELENA
                    ierr = retrieve_var_1D_id(vars_1D_misc,'eq_H',eq_H_misc_id)
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end function prepare_vars_misc
        
        ! prepare eq vars
        integer function prepare_vars_eq() result(ierr)
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
            ierr = retrieve_var_1D_id(vars_1D_eq,'max_flux',max_flux_id)
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
        integer function prepare_vars_X_1(lim_sec_X) result(ierr)
            use X_vars, only: X_1_var_names
            
            character(*), parameter :: rout_name = 'prepare_vars_X_1'
            
            ! input / output
            integer, intent(in), optional :: lim_sec_X(2)                       ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            character(len=max_str_ln), allocatable :: req_var_names(:)          ! requested variable names
            integer :: id                                                       ! counter
            
            ! initialize ierr
            ierr = 0
            
            ! 1. RE_U_0
            req_var_names = get_full_var_names([X_1_var_names(1)],lim_sec_X)
            allocate(RE_U_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_U_0_id(id))
                CHCKERR('')
            end do
            
            ! 2. IM_U_0
            req_var_names = get_full_var_names([X_1_var_names(2)],lim_sec_X)
            allocate(IM_U_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_U_0_id(id))
                CHCKERR('')
            end do
            
            ! 3. RE_U_1
            req_var_names = get_full_var_names([X_1_var_names(3)],lim_sec_X)
            allocate(RE_U_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_U_1_id(id))
                CHCKERR('')
            end do
            
            ! 4. IM_U_1
            req_var_names = get_full_var_names([X_1_var_names(4)],lim_sec_X)
            allocate(IM_U_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_U_1_id(id))
                CHCKERR('')
            end do
            
            ! 5. RE_DU_0
            req_var_names = get_full_var_names([X_1_var_names(5)],lim_sec_X)
            allocate(RE_DU_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_DU_0_id(id))
                CHCKERR('')
            end do
            
            ! 6. IM_DU_0
            req_var_names = get_full_var_names([X_1_var_names(6)],lim_sec_X)
            allocate(IM_DU_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_DU_0_id(id))
                CHCKERR('')
            end do
            
            ! 7. RE_DU_1
            req_var_names = get_full_var_names([X_1_var_names(7)],lim_sec_X)
            allocate(RE_DU_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_DU_1_id(id))
                CHCKERR('')
            end do
            
            ! 8. IM_DU_1
            req_var_names = get_full_var_names([X_1_var_names(8)],lim_sec_X)
            allocate(IM_DU_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_DU_1_id(id))
                CHCKERR('')
            end do
        end function prepare_vars_X_1
        
        ! prepare tensorial perturbation vars
        integer function prepare_vars_X_2(lim_sec_X) result(ierr)
            use X_vars, only: X_2_var_names
            
            character(*), parameter :: rout_name = 'prepare_vars_X_2'
            
            ! input / output
            integer, intent(in), optional :: lim_sec_X(2,2)                     ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            character(len=max_str_ln), allocatable :: req_var_names(:)          ! requested variable names
            integer :: id                                                       ! counter
            
            ! initialize ierr
            ierr = 0
            
            ! 1. RE_PV_int_0
            req_var_names = get_full_var_names([X_2_var_names(1)],&
                &[.true.],lim_sec_X)
            allocate(RE_PV_int_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_PV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 2. IM_PV_int_0
            req_var_names = get_full_var_names([X_2_var_names(2)],&
                &[.true.],lim_sec_X)
            allocate(IM_PV_int_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_PV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 3. RE_PV_int_1
            req_var_names = get_full_var_names([X_2_var_names(3)],&
                &[.false.],lim_sec_X)
            allocate(RE_PV_int_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_PV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 4. IM_PV_int_1
            req_var_names = get_full_var_names([X_2_var_names(4)],&
                &[.false.],lim_sec_X)
            allocate(IM_PV_int_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_PV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 5. RE_PV_int_2
            req_var_names = get_full_var_names([X_2_var_names(5)],&
                &[.true.],lim_sec_X)
            allocate(RE_PV_int_2_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_PV_int_2_id(id))
                CHCKERR('')
            end do
            
            ! 6. IM_PV_int_2
            req_var_names = get_full_var_names([X_2_var_names(6)],&
                &[.true.],lim_sec_X)
            allocate(IM_PV_int_2_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_PV_int_2_id(id))
                CHCKERR('')
            end do
            
            ! 7. RE_KV_int_0
            req_var_names = get_full_var_names([X_2_var_names(7)],&
                &[.true.],lim_sec_X)
            allocate(RE_KV_int_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_KV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 8. IM_KV_int_0
            req_var_names = get_full_var_names([X_2_var_names(8)],&
                &[.true.],lim_sec_X)
            allocate(IM_KV_int_0_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_KV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 9. RE_KV_int_1
            req_var_names = get_full_var_names([X_2_var_names(9)],&
                &[.false.],lim_sec_X)
            allocate(RE_KV_int_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_KV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 10. IM_KV_int_1
            req_var_names = get_full_var_names([X_2_var_names(10)],&
                &[.false.],lim_sec_X)
            allocate(IM_KV_int_1_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_KV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 11. RE_KV_int_2
            req_var_names = get_full_var_names([X_2_var_names(11)],&
                &[.true.],lim_sec_X)
            allocate(RE_KV_int_2_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_KV_int_2_id(id))
                CHCKERR('')
            end do
            
            ! 12. IM_KV_int_2
            req_var_names = get_full_var_names([X_2_var_names(12)],&
                &[.true.],lim_sec_X)
            allocate(IM_KV_int_2_id(size(req_var_names)))
            do id =1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_KV_int_2_id(id))
                CHCKERR('')
            end do
        end function prepare_vars_X_2
        
        ! prepare solution vars
        integer function prepare_vars_sol() result(ierr)
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
        integer function reconstruct_grid_eq(grid_eq,grid_eq_B) result(ierr)
            use grid_vars, only: create_grid
            
            character(*), parameter :: rout_name = 'reconstruct_grid_eq'
            
            ! input / output
            type(grid_type), intent(inout) :: grid_eq                           ! grid to reconstruct
            type(grid_type), intent(inout), optional :: grid_eq_B               ! field-aligned grid to reconstruct
            
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
            ierr = create_grid(grid_eq,n_eq,eq_limits_loc)
            CHCKERR('')
            grid_eq%r_F = vars_1D_eq(r_F_eq_id)%p
            grid_eq%loc_r_F = &
                &vars_1D_eq(r_F_eq_id)%p(eq_limits_loc(1):eq_limits_loc(2))
            grid_eq%r_E = vars_1D_eq(r_E_eq_id)%p
            grid_eq%loc_r_E = &
                &vars_1D_eq(r_E_eq_id)%p(eq_limits_loc(1):eq_limits_loc(2))
            call conv_1D2ND(vars_1D_eq(theta_F_id),dum_3D)
            grid_eq%theta_F = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(zeta_F_id),dum_3D)
            grid_eq%zeta_F = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(theta_E_id),dum_3D)
            grid_eq%theta_E = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(zeta_E_id),dum_3D)
            grid_eq%zeta_E = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
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
        
        ! reconstruct solution grid
        integer function reconstruct_grid_sol(grid) result(ierr)
            use grid_vars, only: create_grid
            
            character(*), parameter :: rout_name = 'reconstruct_grid_sol'
            
            ! input / output
            type(grid_type), intent(inout) :: grid                              ! grid to reconstruct
            
            ! local variables
            integer :: n_sol                                                    ! n of sol grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Creating solution grid')
            call lvl_ud(1)
            
            ! set n_sol
            n_sol = vars_1D_sol(r_F_X_id)%tot_i_max(1)-&
                &vars_1D_sol(r_F_X_id)%tot_i_min(1)+1
            
            ! set up local sol_limits
            sol_limits_loc = [1,n_sol]
            if (present(sol_limits)) sol_limits_loc = sol_limits
            
            ierr = create_grid(grid,n_sol,sol_limits_loc)
            CHCKERR('')
            grid%r_F = vars_1D_sol(r_F_X_id)%p
            grid%loc_r_F = &
                &vars_1D_sol(r_F_X_id)%p(sol_limits_loc(1):sol_limits_loc(2))
            grid%r_E = vars_1D_sol(r_E_X_id)%p
            grid%loc_r_E = &
                &vars_1D_sol(r_E_X_id)%p(sol_limits_loc(1):sol_limits_loc(2))
            call writo('normal size: '//trim(i2str(n_sol)))
            call lvl_ud(-1)
        end function reconstruct_grid_sol
        
        ! reconstruct miscellaneous vars
        integer function reconstruct_vars_misc() result(ierr)
            use num_vars, only: norm_disc_prec_eq, norm_disc_prec_X, &
                &use_pol_flux_F, eq_style, use_normalization, use_pol_flux_E, &
                &use_pol_flux_F, rho_style
            use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, &
                &T_0, vac_perm
            use grid_vars, only: alpha
            use X_vars, only: min_r_sol, max_r_sol, min_n_X, max_n_X, min_m_X, &
                &max_m_X
            use VMEC, only: lasym, lfreeB, mpol, ntor, nfp
            use HELENA, only: ias, nchi
            
            character(*), parameter :: rout_name = 'reconstruct_vars_misc'
            
            ! initialize ierr
            ierr = 0
            
            call writo('Setting miscellaneous variables')
            call lvl_ud(1)
            
            ! eq
            call conv_1D2ND(vars_1D_misc(eq_misc_id),dum_1D)
            PB3D_version = dum_1D(1)
            ! eq_style has already been set
            rho_style = nint(dum_1D(3))
            alpha = dum_1D(4)
            R_0 = dum_1D(5)
            pres_0 = dum_1D(6)
            B_0 = dum_1D(7)
            psi_0 = dum_1D(8)
            rho_0 = dum_1D(9)
            T_0 = dum_1D(10)
            vac_perm = dum_1D(11)
            use_pol_flux_E = .false.
            use_pol_flux_F = .false.
            use_normalization = .false.
            if (dum_1D(12).gt.0) use_pol_flux_E = .true.
            if (dum_1D(13).gt.0) use_pol_flux_F = .true.
            if (dum_1D(14).gt.0) use_normalization = .true.
            norm_disc_prec_eq = nint(dum_1D(15))
            deallocate(dum_1D)
            
            ! eq_V or eq_H, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    call conv_1D2ND(vars_1D_misc(eq_V_misc_id),dum_1D)
                    lasym = .false.
                    lfreeB = .false.
                    if (dum_1D(1).gt.0) lasym = .true.
                    if (dum_1D(2).gt.0) lfreeB = .true.
                    mpol = nint(dum_1D(3))
                    ntor = nint(dum_1D(4))
                    nfp = nint(dum_1D(5))
                    deallocate(dum_1D)
                case (2)                                                        ! HELENA
                    call conv_1D2ND(vars_1D_misc(eq_H_misc_id),dum_1D)
                    ias = nint(dum_1D(1))
                    nchi = nint(dum_1D(2))
                    deallocate(dum_1D)
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! X
            call conv_1D2ND(vars_1D_misc(X_misc_id),dum_1D)
            min_r_sol = dum_1D(1)
            max_r_sol = dum_1D(2)
            min_n_X = nint(dum_1D(3))
            max_n_X = nint(dum_1D(4))
            min_m_X = nint(dum_1D(5))
            max_m_X = nint(dum_1D(6))
            norm_disc_prec_X = nint(dum_1D(7))
            deallocate(dum_1D)
            
            call lvl_ud(-1)
        end function reconstruct_vars_misc
        
        ! reconstruct equilibrium vars
        integer function reconstruct_vars_eq(grid_eq,eq,met) result(ierr)
            use num_vars, only: eq_style, max_deriv
            use met_vars, only: create_met
            use eq_vars, only: create_eq, max_flux_p_E, max_flux_t_E, &
                &max_flux_p_F, max_flux_t_F
            use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, &
                &ntor, lfreeB
            use HELENA, only: R_H, Z_H, chi_H, flux_p_H, nchi
            
            character(*), parameter :: rout_name = 'reconstruct_vars_eq'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(eq_type), intent(inout) :: eq                                  ! equilibrium variables
            type(met_type), intent(inout) :: met                                ! metric variables
            
            ! local variables
            integer :: n_r_eq                                                   ! n_r of eq grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Setting equilibrium')
            call lvl_ud(1)
            
            ierr = create_eq(grid_eq,eq)
            CHCKERR('')
            
            ! set up n_r_eq, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    n_r_eq = vars_1D_eq(R_V_c_id)%tot_i_max(3)-&
                        &vars_1D_eq(R_V_c_id)%tot_i_max(3)+1
                case (2)                                                        ! HELENA
                    n_r_eq = vars_1D_eq(R_H_id)%tot_i_max(2)-&
                        &vars_1D_eq(R_H_id)%tot_i_max(1)+1
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call conv_1D2ND(vars_1D_eq(max_flux_id),dum_1D)
            max_flux_p_E = dum_1D(1)
            max_flux_t_E = dum_1D(2)
            max_flux_p_F = dum_1D(3)
            max_flux_t_F = dum_1D(4)
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_eq(pres_FD_id),dum_2D)
            eq%pres_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(q_saf_FD_id),dum_2D)
            eq%q_saf_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(rot_t_FD_id),dum_2D)
            eq%rot_t_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(flux_p_FD_id),dum_2D)
            eq%flux_p_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(flux_t_FD_id),dum_2D)
            eq%flux_t_FD = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
            deallocate(dum_2D)
            call conv_1D2ND(vars_1D_eq(rho_id),dum_1D)
            eq%rho = dum_1D(eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_eq(S_id),dum_3D)
            eq%S = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(kappa_n_id),dum_3D)
            eq%kappa_n = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(kappa_g_id),dum_3D)
            eq%kappa_g = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(sigma_id),dum_3D)
            eq%sigma = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq(g_FD_id),dum_7D)
            
            ! Set up particular eq variables, depending on equilibrium style
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
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
                    eq%pres_E = dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(q_saf_E_id),dum_2D)
                    eq%q_saf_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(rot_t_E_id),dum_2D)
                    eq%rot_t_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(flux_p_E_id),dum_2D)
                    eq%flux_p_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                    call conv_1D2ND(vars_1D_eq(flux_t_E_id),dum_2D)
                    eq%flux_t_E = &
                        &dum_2D(eq_limits_loc(1):eq_limits_loc(2),:)
                    deallocate(dum_2D)
                case (2)                                                        ! HELENA
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
            
            call lvl_ud(-1)
            
            call writo('Setting metrics')
            call lvl_ud(1)
            
            ierr = create_met(grid_eq,met)
            CHCKERR('')
            
            met%g_FD = &
                        &dum_7D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:,:)
            deallocate(dum_7D)
            call conv_1D2ND(vars_1D_eq(h_FD_id),dum_7D)
            met%h_FD = &
                        &dum_7D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:,:)
            deallocate(dum_7D)
            call conv_1D2ND(vars_1D_eq(jac_FD_id),dum_6D)
            met%jac_FD = &
                        &dum_6D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:)
            deallocate(dum_6D)
            
            call lvl_ud(-1)
        end function reconstruct_vars_eq
        
        ! reconstruct vectorial perturbation vars
        subroutine reconstruct_vars_X_1(grid_eq,X,lim_sec_X)
            use X_vars, only: create_X
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(X_1_type), intent(inout) :: X                                  ! vectorial perturbation variables
            integer, intent(in), optional :: lim_sec_X(2)                       ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            integer :: id                                                       ! counter
            
            call writo('Setting vectorial perturbation')
            call lvl_ud(1)
            
            call create_X(grid_eq,X,lim_sec_X)
            
            ! set up local eq_limits
            eq_limits_loc = [1,grid_eq%n(3)]
            if (present(eq_limits)) eq_limits_loc = eq_limits
            
            ! RE_U_0
            do id = 1,size(RE_U_0_id)
                call conv_1D2ND(vars_1D_X_1(RE_U_0_id(id)),dum_3D)
                X%U_0(:,:,:,id) = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_U_0
            do id = 1,size(IM_U_0_id)
                call conv_1D2ND(vars_1D_X_1(IM_U_0_id(id)),dum_3D)
                X%U_0(:,:,:,id) = X%U_0(:,:,:,id) + &
                    &iu*dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! RE_U_1
            do id = 1,size(RE_U_1_id)
                call conv_1D2ND(vars_1D_X_1(RE_U_1_id(id)),dum_3D)
                X%U_1(:,:,:,id) = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_U_1
            do id = 1,size(IM_U_1_id)
                call conv_1D2ND(vars_1D_X_1(IM_U_1_id(id)),dum_3D)
                X%U_1(:,:,:,id) = X%U_1(:,:,:,id) + &
                    &iu*dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! RE_DU_0
            do id = 1,size(RE_DU_0_id)
                call conv_1D2ND(vars_1D_X_1(RE_DU_0_id(id)),dum_3D)
                X%DU_0(:,:,:,id) = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_DU_0
            do id = 1,size(IM_DU_0_id)
                call conv_1D2ND(vars_1D_X_1(IM_DU_0_id(id)),dum_3D)
                X%DU_0(:,:,:,id) = X%DU_0(:,:,:,id) + &
                    &iu*dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! RE_DU_1
            do id = 1,size(RE_DU_1_id)
                call conv_1D2ND(vars_1D_X_1(RE_DU_1_id(id)),dum_3D)
                X%DU_1(:,:,:,id) = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_DU_1
            do id = 1,size(IM_DU_1_id)
                call conv_1D2ND(vars_1D_X_1(IM_DU_1_id(id)),dum_3D)
                X%DU_1(:,:,:,id) = X%DU_1(:,:,:,id) + &
                    &iu*dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_X_1
        
        ! reconstruct tensorial perturbation vars
        subroutine reconstruct_vars_X_2(grid_eq,X,lim_sec_X)
            use X_vars, only: create_X
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(X_2_type), intent(inout) :: X                                  ! vectorial perturbation variables
            integer, intent(in), optional :: lim_sec_X(2,2)                     ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            integer :: id                                                       ! counter
            
            call writo('Setting tensorial perturbation')
            call lvl_ud(1)
            
            call create_X(grid_eq,X,lim_sec_X)
            
            ! set up local eq_limits
            eq_limits_loc = [1,grid_eq%n(3)]
            if (present(eq_limits)) eq_limits_loc = eq_limits
            
            ! RE_PV_int_0
            do id = 1,size(RE_PV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(RE_PV_int_0_id(id)),dum_2D)
                X%PV_int_0(id,:,:) = &
                    &dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_PV_int_0
            do id = 1,size(IM_PV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(IM_PV_int_0_id(id)),dum_2D)
                X%PV_int_0(id,:,:) = X%PV_int_0(id,:,:) + &
                    &iu*dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_PV_int_1
            do id = 1,size(RE_PV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(RE_PV_int_1_id(id)),dum_2D)
                X%PV_int_1(id,:,:) = &
                    &dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_PV_int_1
            do id = 1,size(IM_PV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(IM_PV_int_1_id(id)),dum_2D)
                X%PV_int_1(id,:,:) = X%PV_int_1(id,:,:) + &
                    &iu*dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_PV_int_2
            do id = 1,size(RE_PV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(RE_PV_int_2_id(id)),dum_2D)
                X%PV_int_2(id,:,:) = &
                    &dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_PV_int_2
            do id = 1,size(IM_PV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(IM_PV_int_2_id(id)),dum_2D)
                X%PV_int_2(id,:,:) = X%PV_int_2(id,:,:) + &
                    &iu*dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_KV_int_0
            do id = 1,size(RE_KV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(RE_KV_int_0_id(id)),dum_2D)
                X%KV_int_0(id,:,:) = &
                    &dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_KV_int_0
            do id = 1,size(IM_KV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(IM_KV_int_0_id(id)),dum_2D)
                X%KV_int_0(id,:,:) = X%KV_int_0(id,:,:) + &
                    &iu*dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_KV_int_1
            do id = 1,size(RE_KV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(RE_KV_int_1_id(id)),dum_2D)
                X%KV_int_1(id,:,:) = &
                    &dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_KV_int_1
            do id = 1,size(IM_KV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(IM_KV_int_1_id(id)),dum_2D)
                X%KV_int_1(id,:,:) = X%KV_int_1(id,:,:) + &
                    &iu*dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_KV_int_2
            do id = 1,size(RE_KV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(RE_KV_int_2_id(id)),dum_2D)
                X%KV_int_2(id,:,:) = &
                    &dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_KV_int_2
            do id = 1,size(IM_KV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(IM_KV_int_2_id(id)),dum_2D)
                X%KV_int_2(id,:,:) = X%KV_int_2(id,:,:) + &
                    &iu*dum_2D(:,eq_limits_loc(1):eq_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_X_2
        
        ! reconstruct solution vars
        subroutine reconstruct_vars_sol(grid,sol,lim_sec_X)
            use sol_vars, only: create_sol
            
            ! input / output
            type(grid_type), intent(in) :: grid                                 ! solution grid
            type(sol_type), intent(inout) :: sol                                ! solution variables
            integer, intent(in), optional :: lim_sec_X(2)                       ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            integer :: n_EV                                                     ! nr. of Eigenvalues
            
            call writo('Setting solution')
            call lvl_ud(1)
            
            ! set n_EV
            n_EV = vars_1D_sol(RE_X_vec_id)%tot_i_max(3)-&
                &vars_1D_sol(RE_X_vec_id)%tot_i_min(3)+1
            
            ! create solution
            call create_sol(grid_sol,sol,n_EV,lim_sec_X)
            
            ! set up local sol_limits
            sol_limits_loc = [1,grid%n(3)]
            if (present(sol_limits)) sol_limits_loc = sol_limits
            
            call conv_1D2ND(vars_1D_sol(RE_X_val_id),dum_1D)
            sol%val = dum_1D
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_sol(IM_X_val_id),dum_1D)
            sol%val = sol%val + iu*dum_1D
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_sol(RE_X_vec_id),dum_3D)
            sol%vec = dum_3D(:,sol_limits_loc(1):sol_limits_loc(2),:)
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_sol(IM_X_vec_id),dum_3D)
            sol%vec = sol%vec + &
                &iu*dum_3D(:,sol_limits_loc(1):sol_limits_loc(2),:)
            deallocate(dum_3D)
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_sol
    end function reconstruct_PB3D
    
    ! Set  all possible  full  variable names  for given  input  names and  mode
    ! numbers. If no  limits for m_X (pol.  flux) or n_X (tor.  flux) are given,
    ! the  whole range  is taken  from  the global  variables min_n_X,  max_n_X,
    ! min_m_X and max_m_X.
    function get_full_var_names_1(var_names,lim_sec_X) result(full_var_names)   ! vectorial version
        use X_vars, only: get_suffix, &
            &min_n_X, max_n_X, min_m_X, max_m_X
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            ! internal variable names
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits for m_X (pol flux) or n_X (tor flux)
        character(len=max_str_ln), allocatable :: full_var_names(:)             ! full HDF5 variable names
        
        ! local variables
        integer :: n_vars                                                       ! nr. of input variables
        integer :: n_mod                                                        ! nr. of secondary modes
        integer :: id, jd                                                       ! counters
        integer :: id_loc                                                       ! counter
        integer :: lim_sec_X_loc(2)                                             ! local lim_sec_X
        
        ! set local lim_sec_X
        if (use_pol_flux_F) then
            lim_sec_X_loc = [min_m_X,max_m_X]
        else
            lim_sec_X_loc = [min_n_X,max_n_X]
        end if
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set n_vars and n_mod
        n_vars = size(var_names)
        n_mod = lim_sec_X_loc(2) - lim_sec_X_loc(1) + 1
        
        ! allocate full HDF5 variable names
        if (allocated(full_var_names)) deallocate(full_var_names)
        allocate(full_var_names(n_vars*n_mod))
        
        ! loop over all modes
        id_loc = 1
        ! loop over all variables
        do id = 1,n_vars
            do jd = 1,n_mod
                full_var_names(id_loc) = trim(var_names(id))//'_'//&
                    &trim(get_suffix(lim_sec_X_loc,jd))
                id_loc = id_loc+1
            end do
        end do
    end function get_full_var_names_1
    function get_full_var_names_2(var_names,sym,lim_sec_X) &
        &result(full_var_names)                                                 ! tensorial version
        use X_vars, only: set_nn_mod, get_suffix, is_necessary_X, &
            &min_n_X, max_n_X, min_m_X, max_m_X
        use num_vars, only: use_pol_flux_F
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            ! internal variable names
        logical, intent(in) :: sym(:)                                           ! if variable is symmetric
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits for m_X (pol flux) or n_X (tor flux) for both dimensions
        character(len=max_str_ln), allocatable :: full_var_names(:)             ! full HDF5 variable names
        
        ! local variables
        integer :: n_vars                                                       ! nr. of input variables
        integer :: n_mod(2)                                                     ! nr. of secondary modes
        integer :: id, jd, kd                                                   ! counters
        integer :: id_loc                                                       ! counter
        integer :: lim_sec_X_loc(2,2)                                           ! local lim_sec_X
        integer :: nn_mod_tot                                                   ! total nn_mod
        
        ! tests
        if (size(var_names).ne.size(sym)) then
            call writo('WARNING: size of sym has to be equal to size of &
                &var_names')
            return
        end if
        
        ! set local lim_sec_X
        if (use_pol_flux_F) then
            lim_sec_X_loc(:,1) = [min_m_X,max_m_X]
            lim_sec_X_loc(:,2) = [min_m_X,max_m_X]
        else
            lim_sec_X_loc(:,1) = [min_n_X,max_n_X]
            lim_sec_X_loc(:,2) = [min_n_X,max_n_X]
        end if
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set n_vars and n_mod
        n_vars = size(var_names)
        n_mod = lim_sec_X_loc(2,:) - lim_sec_X_loc(1,:) + 1
        
        ! allocate full HDF5 variable names
        if (allocated(full_var_names)) deallocate(full_var_names)
        nn_mod_tot = 0
        do id = 1,n_vars
            if (sym(id)) then
                nn_mod_tot = nn_mod_tot + set_nn_mod(lim_sec_X_loc)
            else
                nn_mod_tot = nn_mod_tot + product(n_mod)
            end if
        end do
        allocate(full_var_names(nn_mod_tot))
        
        ! loop over all modes
        id_loc = 1
        ! loop over all variables
        do id = 1,n_vars
            do kd = 1,n_mod(2)
                do jd = 1,n_mod(1)
                    if (is_necessary_X(lim_sec_X_loc,sym(id),[jd,kd])) then
                        full_var_names(id_loc) = trim(var_names(id))//'_'//&
                            &trim(get_suffix(lim_sec_X_loc,jd,kd))
                        id_loc = id_loc+1
                    end if
                end do
            end do
        end do
    end function get_full_var_names_2
    
    ! Retrieves variable index from array 1D equivalents
    integer function retrieve_var_1D_id(vars,var_name,var_id) &
        &result(ierr)
        character(*), parameter :: rout_name = 'retrieve_var_1D'
        
        ! input / output
        type(var_1D_type), intent(in) :: vars(:)                                ! array of 1D variables
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
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:)                      ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1])
    end subroutine conv_1D2ND_1D
    subroutine conv_1D2ND_2D(var_in,var_out)                                    ! 2D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
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
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
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
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
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
        !type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
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
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
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
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
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
