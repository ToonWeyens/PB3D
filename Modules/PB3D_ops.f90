!------------------------------------------------------------------------------!
!   Operations on PB3D output.                                                 !
!------------------------------------------------------------------------------!
module PB3D_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, max_name_ln, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type
    use HDF5_vars, only: var_1D_type

    implicit none
    private
    public init_PB3D_ops, read_PB3D, reconstruct_PB3D
    
contains
    ! Initialize the PB3D routines:
    !   - set the equilibrium style
    !   - for POST, sets the variables "max_it_rich" and "rich_lvl"
    ! Note: This routine must be called after "open_input" and "read_input"
    integer function init_PB3D_ops() result(ierr)
        use num_vars, only: eq_style, PB3D_name, rank
        use PB3D_vars, only: vars_1D_misc
        use PB3D_utilities, only: retrieve_var_1D_id, conv_1D2ND
        use HDF5_ops, only: read_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D
        
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
            call dealloc_var_1D(vars_1D_misc)
            
            ! user output
            call lvl_ud(-1)
            call writo('PB3D operations Initialized')
        end if
    end function init_PB3D_ops
    
    ! Reads PB3D output from user-provided input file.
    ! Optionally, if perturbation  quantities are read, the limits  of m_X (pol.
    ! flux) or  n_X (tor.  flux) can be  passed, instead of  reading all  of the
    ! perturbation quantities.
    integer function read_PB3D(read_misc,read_grid_eq,read_grid_X,&
        &read_grid_sol,read_eq,read_X_1,read_X_2,read_sol,lim_sec_X_1,&
        &lim_sec_X_2) result(ierr)
        use num_vars, only: eq_style, PB3D_name
        use HDF5_ops, only: read_HDF5_arrs
        use X_vars, only: X_1_var_names, X_2_var_names
        use PB3D_utilities, only: get_full_var_names
        use PB3D_vars, only: vars_1D_misc, vars_1D_grid_eq, vars_1D_grid_eq_B, &
            &vars_1D_grid_X, vars_1D_grid_X_B, vars_1D_grid_sol, vars_1D_eq, &
            &vars_1D_X_1, vars_1D_X_2, vars_1D_sol
        use rich, only: rich_info_short
        
        character(*), parameter :: rout_name = 'read_PB3D'
        
        ! input / output
        logical, intent(in) :: read_misc, read_grid_eq, &
            &read_grid_X, read_grid_sol, read_eq, read_X_1, read_X_2, read_sol  ! which files to read
        integer, intent(in), optional :: lim_sec_X_1(2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in), optional :: lim_sec_X_2(2,2)                       ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_name_ln), allocatable :: req_var_names(:)             ! requested variable names
        
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
        
        if (read_grid_eq) then
            call writo('Reading equilibrium grid variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_grid_eq,PB3D_name,'grid_eq')
            CHCKERR('')
            select case (eq_style)
                case (1)                                                        ! VMEC
                    allocate(vars_1D_grid_eq_B(0))
                case (2)                                                        ! HELENA
                    ierr = read_HDF5_arrs(vars_1D_grid_eq_B,PB3D_name,&
                        &'grid_eq_B')
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            call lvl_ud(-1)
            call writo('Equilbrium grid variables read')
        end if
        
        if (read_grid_X) then
            call writo('Reading perturbation grid variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_grid_X,PB3D_name,&
                &'grid_X'//trim(rich_info_short()))
            CHCKERR('')
            select case (eq_style)
                case (1)                                                        ! VMEC
                    allocate(vars_1D_grid_X_B(0))
                case (2)                                                        ! HELENA
                    ierr = read_HDF5_arrs(vars_1D_grid_X_B,PB3D_name,&
                        &'grid_X_B'//trim(rich_info_short()))
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            call lvl_ud(-1)
            call writo('Perturbation grid variables read')
        end if
        
        if (read_grid_sol) then
            call writo('Reading solution grid variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_grid_sol,PB3D_name,&
                &'grid_sol'//trim(rich_info_short()))
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Solution grid variables read')
        end if
        
        if (read_eq) then
            call writo('Reading equilibrium variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_eq,PB3D_name,'eq')
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Equilbrium variables read')
        end if
        
        if (read_X_1) then
            call writo('Reading vectorial perturbation variables')
            call lvl_ud(1)
            if (present(lim_sec_X_1)) then
                ! get requested full variable names
                call get_full_var_names(X_1_var_names,req_var_names,lim_sec_X_1)
                
                ! read from HDF5 file
                ierr = read_HDF5_arrs(vars_1D_X_1,PB3D_name,&
                &'X_1'//trim(rich_info_short()),acc_var_names=req_var_names)
                CHCKERR(err_msg)
            else
                ierr = read_HDF5_arrs(vars_1D_X_1,PB3D_name,&
                &'X_1'//trim(rich_info_short()))
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
                call get_full_var_names(X_2_var_names,&
                    &[.true.,.true.,.false.,.false.,.true.,.true.,&
                    &.true.,.true.,.false.,.false.,.true.,.true.],&
                    &req_var_names,lim_sec_X_2)
                
                ! read from HDF5 file
                ierr = read_HDF5_arrs(vars_1D_X_2,PB3D_name,&
                &'X_2'//trim(rich_info_short()),acc_var_names=req_var_names)
                CHCKERR(err_msg)
            else
                ierr = read_HDF5_arrs(vars_1D_X_2,PB3D_name,&
                &'X_2'//trim(rich_info_short()))
                CHCKERR(err_msg)
            end if
            call lvl_ud(-1)
            call writo('Tensorial perturbation variables read')
        end if
        
        if (read_sol) then
            call writo('Reading solution variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_sol,PB3D_name,&
                &'sol'//trim(rich_info_short()))
            CHCKERR(err_msg)
            call lvl_ud(-1)
            call writo('Solution variables read')
        end if
    end function read_PB3D
    
    ! Reconstructs the PB3D variables: eq,  X and/or sol variables, as indicated
    ! through rec_eq, rec_X_1, rec_X_2 and rec_sol.
    ! Additionally, if the variables 'eq_limits', 'X_limits' and/or 'sol_limits'
    ! are not provided, the full normal range is taken for every variable.
    ! Optionally, if perturbation  quantities are read, the limits  of m_X (pol.
    ! flux) or  n_X (tor.  flux) can be  passed, instead of  reading all  of the
    ! perturbation quantities.
    ! Note: To reconstruct perturbation variables, mode numbers n_X and m_X have
    ! to be  correctly set up, using  the equilibrium and perturbation  grids as
    ! well as the equilibrium variables. Therefore, if this routine detects that
    ! the  equilibrium  and perturbation  grid  are  reconstructed, as  well  as
    ! equilibrium and perturbation or  solution variables, it automatically also
    ! sets  up  n_X  and  m_X,  a  behavior  that  can  be  switched  off  using
    ! "use_setup_nm_X".
    integer function reconstruct_PB3D(rec_misc,rec_grid_eq,rec_grid_X,&
        &rec_grid_sol,rec_eq,rec_X_1,rec_X_2,rec_sol,grid_eq,grid_eq_B,grid_X,&
        &grid_X_B,grid_sol,eq,met,X_1,X_2,sol,eq_limits,X_limits,sol_limits,&
        &lim_sec_X_1,lim_sec_X_2,use_setup_nm_X) result(ierr)
        
        use HDF5_vars, only: dealloc_var_1D
        use PB3D_utilities, only: get_full_var_names, retrieve_var_1D_id, &
            &conv_1D2ND
        use PB3D_vars, only: vars_1D_misc, vars_1D_grid_eq, vars_1D_grid_eq_B, &
            &vars_1D_grid_X, vars_1D_grid_X_B, vars_1D_grid_sol, vars_1D_eq, &
            &vars_1D_X_1, vars_1D_X_2, vars_1D_sol, min_PB3D_version, &
            &PB3D_version
        use MPI_utilities, only: wait_MPI
        use X_ops, only: setup_nm_X
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D'
        
        ! input  / output
        logical, intent(in) :: rec_misc, rec_grid_eq, rec_grid_X, &
            &rec_grid_sol, rec_eq, rec_X_1, rec_X_2, rec_sol                    ! which quantities to reconstruct
        type(grid_type), intent(inout), optional :: grid_eq                     ! equilibrium grid 
        type(grid_type), intent(inout), optional, pointer :: grid_eq_B          ! optional field-aligned equilibrium grid for HELENA
        type(grid_type), intent(inout), optional :: grid_X                      ! perturbation grid 
        type(grid_type), intent(inout), optional, pointer :: grid_X_B           ! optional field-aligned grid perturbation grid for HELENA
        type(grid_type), intent(inout), optional :: grid_sol                    ! solution grid 
        type(eq_type), intent(inout), optional :: eq                            ! equilibrium variables
        type(met_type), intent(inout), optional :: met                          ! metric variables
        type(X_1_type), intent(inout), optional :: X_1                          ! vectorial perturbation variables
        type(X_2_type), intent(inout), optional :: X_2                          ! tensorial perturbation variables
        type(sol_type), intent(inout), optional :: sol                          ! solution variables
        integer, intent(in), optional :: eq_limits(2)                           ! i_limit of eq variables
        integer, intent(in), optional :: X_limits(2)                            ! i_limit of X variables
        integer, intent(in), optional :: sol_limits(2)                          ! i_limit of sol variables
        integer, intent(in), optional :: lim_sec_X_1(2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in), optional :: lim_sec_X_2(2,2)                       ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        logical, intent(in), optional :: use_setup_nm_X                         ! setup n_X and m_X if possible [.true.]
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: use_setup_nm_X_loc                                           ! local version of use_setup_nm_X
        
        ! local variables also used in child functions
        ! miscellaneous
        integer :: eq_misc_id, X_misc_id, eq_V_misc_id, eq_H_misc_id            ! indices of miscellaneous variables
        ! equilibrium grid
        integer :: r_F_eq_id, r_E_eq_id                                         ! index of equilibrium r_F and r_E
        integer :: theta_F_eq_id, zeta_F_eq_id                                  ! index of equilibrium theta_F and zeta_F
        integer :: theta_E_eq_id, zeta_E_eq_id                                  ! index of equilibrium theta_E and zeta_E
        integer :: r_F_eq_B_id, r_E_eq_B_id                                     ! index of field-aligned equilibrium r_F and r_E
        integer :: theta_F_eq_B_id, zeta_F_eq_B_id                              ! index of field-aligned equilibrium theta_F and zeta_F
        integer :: theta_E_eq_B_id, zeta_E_eq_B_id                              ! index of field-aligned equilibrium theta_E and zeta_E
        ! perturbation grid
        integer :: r_F_X_id, r_E_X_id                                           ! index of perturbation r_F and r_E
        integer :: theta_F_X_id, zeta_F_X_id                                    ! index of perturbation theta_F and zeta_F
        integer :: theta_E_X_id, zeta_E_X_id                                    ! index of perturbation theta_E and zeta_E
        integer :: r_F_X_B_id, r_E_X_B_id                                       ! index of field-aligned perturbation r_F and r_E
        integer :: theta_F_X_B_id, zeta_F_X_B_id                                ! index of field-aligned perturbation theta_F and zeta_F
        integer :: theta_E_X_B_id, zeta_E_X_B_id                                ! index of field-aligned perturbation theta_E and zeta_E
        ! solution grid
        integer :: r_F_sol_id, r_E_sol_id                                       ! index of solution r_F and r_E
        ! equilibrium
        integer :: max_flux_id                                                  ! index of max_flux
        integer :: pres_FD_id, q_saf_FD_id, rot_t_FD_id                         ! index of pres_FD, q_saf_FD, rot_t_FD
        integer :: flux_p_FD_id, flux_t_FD_id                                   ! index of flux_p_FD, flux_t_FD
        integer :: pres_E_id, q_saf_E_id, rot_t_E_id                            ! index of pres_E, q_saf_E, rot_t_E
        integer :: flux_p_E_id, flux_t_E_id                                     ! index of flux_p_FD, flux_t_FD
        integer :: rho_id, S_id, kappa_n_id, kappa_g_id, sigma_id               ! index of rho, S, kappa_n, kappa_g, sigma
        integer :: g_FD_id, h_FD_id, jac_FD_id                                  ! index of g_FD, h_FD, jac_FD
        integer :: R_V_c_id, R_V_s_id                                           ! index of R_V_c and R_V_s
        integer :: Z_V_c_id, Z_V_s_id                                           ! index of Z_V_c and Z_V_s
        integer :: L_V_c_id, L_V_s_id                                           ! index of L_V_c and L_V_s
        integer :: R_H_id, Z_H_id, chi_H_id, flux_p_H_id                        ! index of R_H, Z_H, chi_H and flux_p_H
        ! vectorial perturbation
        integer, allocatable :: RE_U_0_id(:), IM_U_0_id(:)                      ! index of RE_U_0, IM_U_0, RE_U_1, IM_U_1
        integer, allocatable :: RE_U_1_id(:), IM_U_1_id(:)                      ! index of RE_U_1, IM_U_1, RE_U_1, IM_U_1
        integer, allocatable :: RE_DU_0_id(:), IM_DU_0_id(:)                    ! index of RE_DU_0, IM_DU_0, RE_DU_1, IM_DU_1
        integer, allocatable :: RE_DU_1_id(:), IM_DU_1_id(:)                    ! index of RE_DU_1, IM_DU_1, RE_DU_1, IM_DU_1
        ! tensorial perturbation
        integer, allocatable :: RE_PV_int_0_id(:), IM_PV_int_0_id(:)            ! index of RE_PV_int_0, IM_PV_int_0
        integer, allocatable :: RE_PV_int_1_id(:), IM_PV_int_1_id(:)            ! index of RE_PV_int_1, IM_PV_int_1
        integer, allocatable :: RE_PV_int_2_id(:), IM_PV_int_2_id(:)            ! index of RE_PV_int_2, IM_PV_int_2
        integer, allocatable :: RE_KV_int_0_id(:), IM_KV_int_0_id(:)            ! index of RE_KV_int_0, IM_KV_int_0
        integer, allocatable :: RE_KV_int_1_id(:), IM_KV_int_1_id(:)            ! index of RE_KV_int_1, IM_KV_int_1
        integer, allocatable :: RE_KV_int_2_id(:), IM_KV_int_2_id(:)            ! index of RE_KV_int_2, IM_KV_int_2
        ! solution
        integer :: RE_sol_val_id, IM_sol_val_id                                 ! index of RE_sol_val, IM_sol_val
        integer :: RE_sol_vec_id, IM_sol_vec_id                                 ! index of RE_sol_vec, IM_sol_vec
        ! helper variables
        real(dp), allocatable :: dum_1D(:), dum_2D(:,:), dum_3D(:,:,:)          ! dummy variables
        real(dp), allocatable :: dum_4D(:,:,:,:)                                ! dummy variables
        !real(dp), allocatable :: dum_5D(:,:,:,:,:)                              ! dummy variables
        real(dp), allocatable :: dum_6D(:,:,:,:,:,:), dum_7D(:,:,:,:,:,:,:)     ! dummy variables
        real(dp), parameter :: tol_version = 1.E-8_dp                           ! tolerance for version control
        integer :: eq_limits_loc(2), X_limits_loc(2), sol_limits_loc(2)         ! local versions of eq_limits, sol_limits
        
        ! initialize ierr
        ierr = 0
        
        ! test whether appropriate variables provided
        if (rec_eq .and. .not.&
            &(present(grid_eq).and.present(eq).and.present(met))) then
            ierr = 1
            err_msg = 'To reconstruct equilibrium need grid_eq, eq and met'
            CHCKERR(err_msg)
        end if
        if (rec_X_1 .and. .not.(present(grid_X).and.present(X_1))) then
            ierr = 1
            err_msg = 'To reconstruct vectorial perturbation need grid_X &
                &and X_1'
            CHCKERR(err_msg)
        end if
        if (rec_X_2 .and. .not.(present(grid_X).and.present(X_2))) then
            ierr = 1
            err_msg = 'To reconstruct tensorial perturbation need grid_X &
                &and X_2'
            CHCKERR(err_msg)
        end if
        if (rec_sol .and. .not.(present(grid_sol).and.present(sol))) then
            ierr = 1
            err_msg = 'To reconstruct equilibrium need grid_sol and sol'
            CHCKERR(err_msg)
        end if
        
        ! set up local use_setup_nm_X
        if (rec_grid_eq .and. rec_grid_X .and. rec_eq .and. &
            &(rec_X_1.or.rec_X_2.or.rec_sol)) then                              ! use setup_nm_X if not overriden
            use_setup_nm_X_loc = .true.
            if (present(use_setup_nm_X)) use_setup_nm_X_loc = use_setup_nm_X    ! override
        else                                                                    ! not possible to use setup_nm_X
            use_setup_nm_X_loc = .false.
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
        
        if (rec_grid_eq) then
            ! user output
            call writo('Prepare equilibrium grid variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_grid_eq()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing equilibrium grid variables')
            call lvl_ud(1)
            
            call writo('Setting Variables')
            call lvl_ud(1)
            ierr = reconstruct_vars_grid_eq(grid_eq,grid_eq_B)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_grid_eq)
            call dealloc_var_1D(vars_1D_grid_eq_B)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Equilibrium grid variables reconstructed')
        end if
        
        if (rec_grid_X) then
            ! user output
            call writo('Prepare perturbation grid variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_grid_X()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing perturbation grid variables')
            call lvl_ud(1)
            
            call writo('Setting Variables')
            call lvl_ud(1)
            ierr = reconstruct_vars_grid_X(grid_X,grid_X_B)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_grid_X)
            call dealloc_var_1D(vars_1D_grid_X_B)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Perturbation grid variables reconstructed')
        end if
        
        if (rec_grid_sol) then
            ! user output
            call writo('Prepare solution grid variable indices')
            call lvl_ud(1)
            ierr = prepare_vars_grid_sol()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing solution grid variables')
            call lvl_ud(1)
            
            call writo('Setting Variables')
            call lvl_ud(1)
            ierr = reconstruct_vars_grid_sol(grid_sol)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_grid_sol)
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            call writo('Solution grid variables reconstructed')
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
            
            call writo('Setting variables')
            call lvl_ud(1)
            ierr = reconstruct_vars_eq(grid_eq,eq,met)
            CHCKERR('')
            call lvl_ud(-1)
            
            if (use_setup_nm_X_loc) then
                call writo('Setup mode information')
                call lvl_ud(1)
                ierr = setup_nm_X(grid_eq,grid_X,eq)
                CHCKERR('')
                call lvl_ud(-1)
            end if
            
            call writo('Deallocating temporary variables')
            call lvl_ud(1)
            call dealloc_var_1D(vars_1D_eq)
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
            call reconstruct_vars_X_1(grid_X,X_1,lim_sec_X_1)
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
            ierr = wait_MPI()
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Indices prepared')
            
            call writo('Reconstructing tensorial perturbation variables')
            call lvl_ud(1)
            
            call writo('Setting variables')
            call lvl_ud(1)
            call reconstruct_vars_X_2(grid_X,X_2,lim_sec_X_2)
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
        
        ! prepare grid_eq vars
        integer function prepare_vars_grid_eq() result(ierr)
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'prepare_vars_grid_eq'
            
            ! initialize ierr
            ierr = 0
            
            ! set up 1D indices for grid
            ierr = retrieve_var_1D_id(vars_1D_grid_eq,'r_F',r_F_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_eq,'r_E',r_E_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_eq,'theta_F',theta_F_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_eq,'zeta_F',zeta_F_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_eq,'theta_E',theta_E_eq_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_eq,'zeta_E',zeta_E_eq_id)
            CHCKERR('')
            
            ! set up 1D indices for field-aligned grid if present
            if (present(grid_eq_B)) then 
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! do nothing
                    case (2)                                                    ! HELENA
                        ierr = retrieve_var_1D_id(vars_1D_grid_eq_B,'r_F',&
                            &r_F_eq_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_eq_B,'r_E',&
                            &r_E_eq_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_eq_B,'theta_F',&
                            &theta_F_eq_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_eq_B,'zeta_F',&
                            &zeta_F_eq_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_eq_B,'theta_E',&
                            &theta_E_eq_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_eq_B,'zeta_E',&
                            &zeta_E_eq_B_id)
                        CHCKERR('')
                    case default
                        ierr = 1
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        CHCKERR(err_msg)
                end select
            end if
        end function prepare_vars_grid_eq
        
        ! prepare grid_X vars
        integer function prepare_vars_grid_X() result(ierr)
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'prepare_vars_grid_X'
            
            ! initialize ierr
            ierr = 0
            
            ! set up 1D indices for grid
            ierr = retrieve_var_1D_id(vars_1D_grid_X,'r_F',r_F_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_X,'r_E',r_E_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_X,'theta_F',theta_F_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_X,'zeta_F',zeta_F_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_X,'theta_E',theta_E_X_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_X,'zeta_E',zeta_E_X_id)
            CHCKERR('')
            
            ! set up 1D indices for field-aligned grid if present
            if (present(grid_X_B)) then 
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! do nothing
                    case (2)                                                    ! HELENA
                        ierr = retrieve_var_1D_id(vars_1D_grid_X_B,'r_F',&
                            &r_F_X_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_X_B,'r_E',&
                            &r_E_X_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_X_B,'theta_F',&
                            &theta_F_X_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_X_B,'zeta_F',&
                            &zeta_F_X_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_X_B,'theta_E',&
                            &theta_E_X_B_id)
                        CHCKERR('')
                        ierr = retrieve_var_1D_id(vars_1D_grid_X_B,'zeta_E',&
                            &zeta_E_X_B_id)
                        CHCKERR('')
                    case default
                        ierr = 1
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        CHCKERR(err_msg)
                end select
            end if
        end function prepare_vars_grid_X
        
        ! prepare grid_sol vars
        integer function prepare_vars_grid_sol() result(ierr)
            character(*), parameter :: rout_name = 'prepare_vars_grid_sol'
            
            ! initialize ierr
            ierr = 0
            
            ! set up 1D indices for grid
            ierr = retrieve_var_1D_id(vars_1D_grid_sol,'r_F',r_F_sol_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_grid_sol,'r_E',r_E_sol_id)
            CHCKERR('')
        end function prepare_vars_grid_sol
        
        ! prepare eq vars
        integer function prepare_vars_eq() result(ierr)
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'prepare_vars_eq'
            
            ! initialize ierr
            ierr = 0
            
            ! set up 1D indices common for all equilibrium styles
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
        end function prepare_vars_eq
        
        ! prepare vectorial perturbation vars
        integer function prepare_vars_X_1(lim_sec_X) result(ierr)
            use X_vars, only: X_1_var_names, n_mod_X
            
            character(*), parameter :: rout_name = 'prepare_vars_X_1'
            
            ! input / output
            integer, intent(in), optional :: lim_sec_X(2)                       ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            character(len=max_name_ln), allocatable :: req_var_names(:)         ! requested variable names
            integer :: id                                                       ! counter
            integer :: lim_sec_X_loc(2)                                         ! local version of lim_sec_X
            
            ! initialize ierr
            ierr = 0
            
            ! set up local lim_sec_X
            lim_sec_X_loc = [1,n_mod_X]
            if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
            
            ! 1. RE_U_0
            call get_full_var_names([X_1_var_names(1)],req_var_names,&
                &lim_sec_X_loc)
            allocate(RE_U_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_U_0_id(id))
                CHCKERR('')
            end do
            
            ! 2. IM_U_0
            call get_full_var_names([X_1_var_names(2)],req_var_names,&
                &lim_sec_X_loc)
            allocate(IM_U_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_U_0_id(id))
                CHCKERR('')
            end do
            
            ! 3. RE_U_1
            call get_full_var_names([X_1_var_names(3)],req_var_names,&
                &lim_sec_X_loc)
            allocate(RE_U_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_U_1_id(id))
                CHCKERR('')
            end do
            
            ! 4. IM_U_1
            call get_full_var_names([X_1_var_names(4)],req_var_names,&
                &lim_sec_X_loc)
            allocate(IM_U_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_U_1_id(id))
                CHCKERR('')
            end do
            
            ! 5. RE_DU_0
            call get_full_var_names([X_1_var_names(5)],req_var_names,&
                &lim_sec_X_loc)
            allocate(RE_DU_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_DU_0_id(id))
                CHCKERR('')
            end do
            
            ! 6. IM_DU_0
            call get_full_var_names([X_1_var_names(6)],req_var_names,&
                &lim_sec_X_loc)
            allocate(IM_DU_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_DU_0_id(id))
                CHCKERR('')
            end do
            
            ! 7. RE_DU_1
            call get_full_var_names([X_1_var_names(7)],req_var_names,&
                &lim_sec_X_loc)
            allocate(RE_DU_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &RE_DU_1_id(id))
                CHCKERR('')
            end do
            
            ! 8. IM_DU_1
            call get_full_var_names([X_1_var_names(8)],req_var_names,&
                &lim_sec_X_loc)
            allocate(IM_DU_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_1,req_var_names(id),&
                    &IM_DU_1_id(id))
                CHCKERR('')
            end do
        end function prepare_vars_X_1
        
        ! prepare tensorial perturbation vars
        integer function prepare_vars_X_2(lim_sec_X) result(ierr)
            use X_vars, only: X_2_var_names, n_mod_X
            
            character(*), parameter :: rout_name = 'prepare_vars_X_2'
            
            ! input / output
            integer, intent(in), optional :: lim_sec_X(2,2)                     ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            character(len=max_name_ln), allocatable :: req_var_names(:)         ! requested variable names
            integer :: id                                                       ! counter
            integer :: lim_sec_X_loc(2,2)                                       ! local version of lim_sec_X
            
            ! initialize ierr
            ierr = 0
            
            ! set up local lim_sec_X
            lim_sec_X_loc(:,1) = [1,n_mod_X]
            lim_sec_X_loc(:,2) = [1,n_mod_X]
            if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
            
            ! 1. RE_PV_int_0
            call get_full_var_names([X_2_var_names(1)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(RE_PV_int_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_PV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 2. IM_PV_int_0
            call get_full_var_names([X_2_var_names(2)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(IM_PV_int_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_PV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 3. RE_PV_int_1
            call get_full_var_names([X_2_var_names(3)],[.false.],&
                &req_var_names,lim_sec_X_loc)
            allocate(RE_PV_int_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_PV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 4. IM_PV_int_1
            call get_full_var_names([X_2_var_names(4)],[.false.],&
                &req_var_names,lim_sec_X_loc)
            allocate(IM_PV_int_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_PV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 5. RE_PV_int_2
            call get_full_var_names([X_2_var_names(5)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(RE_PV_int_2_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_PV_int_2_id(id))
                CHCKERR('')
            end do
            
            ! 6. IM_PV_int_2
            call get_full_var_names([X_2_var_names(6)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(IM_PV_int_2_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_PV_int_2_id(id))
                CHCKERR('')
            end do
            
            ! 7. RE_KV_int_0
            call get_full_var_names([X_2_var_names(7)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(RE_KV_int_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_KV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 8. IM_KV_int_0
            call get_full_var_names([X_2_var_names(8)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(IM_KV_int_0_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_KV_int_0_id(id))
                CHCKERR('')
            end do
            
            ! 9. RE_KV_int_1
            call get_full_var_names([X_2_var_names(9)],[.false.],&
                &req_var_names,lim_sec_X_loc)
            allocate(RE_KV_int_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_KV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 10. IM_KV_int_1
            call get_full_var_names([X_2_var_names(10)],[.false.],&
                &req_var_names,lim_sec_X_loc)
            allocate(IM_KV_int_1_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &IM_KV_int_1_id(id))
                CHCKERR('')
            end do
            
            ! 11. RE_KV_int_2
            call get_full_var_names([X_2_var_names(11)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(RE_KV_int_2_id(size(req_var_names)))
            do id = 1,size(req_var_names)
                ierr = retrieve_var_1D_id(vars_1D_X_2,req_var_names(id),&
                    &RE_KV_int_2_id(id))
                CHCKERR('')
            end do
            
            ! 12. IM_KV_int_2
            call get_full_var_names([X_2_var_names(12)],[.true.],&
                &req_var_names,lim_sec_X_loc)
            allocate(IM_KV_int_2_id(size(req_var_names)))
            do id = 1,size(req_var_names)
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
            ierr = retrieve_var_1D_id(vars_1D_sol,'RE_sol_val',RE_sol_val_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'IM_sol_val',IM_sol_val_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'RE_sol_vec',RE_sol_vec_id)
            CHCKERR('')
            ierr = retrieve_var_1D_id(vars_1D_sol,'IM_sol_vec',IM_sol_vec_id)
            CHCKERR('')
        end function prepare_vars_sol
        
        ! reconstruct miscellaneous vars
        integer function reconstruct_vars_misc() result(ierr)
            use num_vars, only: norm_disc_prec_eq, norm_disc_prec_X, &
                &use_pol_flux_F, eq_style, use_normalization, use_pol_flux_E, &
                &use_pol_flux_F, rho_style, norm_style, U_style, X_style
            use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, &
                &T_0, vac_perm
            use grid_vars, only: alpha
            use X_vars, only: min_r_sol, max_r_sol, prim_X, min_sec_X, &
                &max_sec_X, n_mod_X
            use VMEC, only: lasym, lfreeB, mpol, ntor, nfp
            use HELENA_vars, only: ias, nchi
            
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
            prim_X = nint(dum_1D(3))
            n_mod_X = nint(dum_1D(4))
            min_sec_X = nint(dum_1D(5))
            max_sec_X = nint(dum_1D(6))
            norm_disc_prec_X = nint(dum_1D(7))
            norm_style = nint(dum_1D(8))
            U_style = nint(dum_1D(9))
            X_style = nint(dum_1D(10))
            deallocate(dum_1D)
            
            call lvl_ud(-1)
        end function reconstruct_vars_misc
        
        ! reconstruct equilibrium grid
        integer function reconstruct_vars_grid_eq(grid_eq,grid_eq_B) &
            &result(ierr)
            use grid_vars, only: create_grid
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'reconstruct_vars_grid_eq'
            
            ! input / output
            type(grid_type), intent(inout), target :: grid_eq                   ! grid to reconstruct
            type(grid_type), intent(inout), optional, pointer :: grid_eq_B      ! field-aligned grid to reconstruct
            
            ! local variables
            integer :: n_eq(3), n_eq_B(3)                                       ! n of eq and eq_B grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Creating equilibrium grid')
            call lvl_ud(1)
            
            ! set n_eq
            n_eq = vars_1D_grid_eq(theta_F_eq_id)%tot_i_max-&
                &vars_1D_grid_eq(theta_F_eq_id)%tot_i_min+1
            
            ! set up local eq_limits
            eq_limits_loc = [1,n_eq(3)]
            if (present(eq_limits)) eq_limits_loc = eq_limits
            
            ! set grid
            ierr = create_grid(grid_eq,n_eq,eq_limits_loc)
            CHCKERR('')
            grid_eq%r_F = vars_1D_grid_eq(r_F_eq_id)%p
            grid_eq%loc_r_F = &
                &vars_1D_grid_eq(r_F_eq_id)%p(eq_limits_loc(1):eq_limits_loc(2))
            grid_eq%r_E = vars_1D_grid_eq(r_E_eq_id)%p
            grid_eq%loc_r_E = &
                &vars_1D_grid_eq(r_E_eq_id)%p(eq_limits_loc(1):eq_limits_loc(2))
            call conv_1D2ND(vars_1D_grid_eq(theta_F_eq_id),dum_3D)
            grid_eq%theta_F = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_grid_eq(zeta_F_eq_id),dum_3D)
            grid_eq%zeta_F = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_grid_eq(theta_E_eq_id),dum_3D)
            grid_eq%theta_E = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_grid_eq(zeta_E_eq_id),dum_3D)
            grid_eq%zeta_E = dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
            deallocate(dum_3D)
            call writo('angular size: ('//trim(i2str(n_eq(1)))//','//&
                &trim(i2str(n_eq(2)))//')')
            call writo('normal size: '//trim(i2str(n_eq(3))))
            call lvl_ud(-1)
            
            if (present(grid_eq_B)) then
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        grid_eq_B => grid_eq
                    case (2)                                                    ! HELENA
                        call writo('Creating field-aligned equilibrium grid')
                        call lvl_ud(1)
                        n_eq_B = vars_1D_grid_eq_B(theta_F_eq_B_id)%tot_i_max-&
                            &vars_1D_grid_eq_B(theta_F_eq_B_id)%tot_i_min+1
                        allocate(grid_eq_B)
                        ierr = create_grid(grid_eq_B,n_eq_B,eq_limits_loc)
                        CHCKERR('')
                        grid_eq_B%r_F = vars_1D_grid_eq_B(r_F_eq_B_id)%p
                        grid_eq_B%loc_r_F = vars_1D_grid_eq_B(r_F_eq_B_id)%&
                            &p(eq_limits_loc(1):eq_limits_loc(2))
                        grid_eq_B%r_E = vars_1D_grid_eq_B(r_E_eq_B_id)%p
                        grid_eq_B%loc_r_E = vars_1D_grid_eq_B(r_E_eq_B_id)%&
                            &p(eq_limits_loc(1):eq_limits_loc(2))
                        call conv_1D2ND(vars_1D_grid_eq_B(theta_F_eq_B_id),&
                            &dum_3D)
                        grid_eq_B%theta_F = &
                            &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                        deallocate(dum_3D)
                        call conv_1D2ND(vars_1D_grid_eq_B(zeta_F_eq_B_id),&
                            &dum_3D)
                        grid_eq_B%zeta_F = &
                            &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                        deallocate(dum_3D)
                        call conv_1D2ND(vars_1D_grid_eq_B(theta_E_eq_B_id),&
                            &dum_3D)
                        grid_eq_B%theta_E = &
                            &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                        deallocate(dum_3D)
                        call conv_1D2ND(vars_1D_grid_eq_B(zeta_E_eq_B_id),&
                            &dum_3D)
                        grid_eq_B%zeta_E = &
                            &dum_3D(:,:,eq_limits_loc(1):eq_limits_loc(2))
                        deallocate(dum_3D)
                        call writo('angular size: ('//trim(i2str(n_eq_B(1)))//&
                            &','//trim(i2str(n_eq_B(2)))//')')
                        call writo('normal size: '//trim(i2str(n_eq_B(3))))
                        call lvl_ud(-1)
                    case default
                        ierr = 1
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        CHCKERR(err_msg)
                end select
            end if
        end function reconstruct_vars_grid_eq
        
        ! reconstruct perturbation grid
        integer function reconstruct_vars_grid_X(grid_X,grid_X_B) result(ierr)
            use grid_vars, only: create_grid
            use num_vars, only: eq_style
            
            character(*), parameter :: rout_name = 'reconstruct_vars_grid_X'
            
            ! input / output
            type(grid_type), intent(inout), target :: grid_X                    ! grid to reconstruct
            type(grid_type), intent(inout), optional, pointer :: grid_X_B       ! field-aligned grid to reconstruct
            
            ! local variables
            integer :: n_X(3), n_X_B(3)                                         ! n of X and X_B grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Creating perturbation grid')
            call lvl_ud(1)
            
            ! set n_X
            n_X = vars_1D_grid_X(theta_F_X_id)%tot_i_max-&
                &vars_1D_grid_X(theta_F_X_id)%tot_i_min+1
            
            ! set up local X_limits
            X_limits_loc = [1,n_X(3)]
            if (present(X_limits)) X_limits_loc = X_limits
            
            ! set grid
            ierr = create_grid(grid_X,n_X,X_limits_loc)
            CHCKERR('')
            grid_X%r_F = vars_1D_grid_X(r_F_X_id)%p
            grid_X%loc_r_F = &
                &vars_1D_grid_X(r_F_X_id)%p(X_limits_loc(1):X_limits_loc(2))
            grid_X%r_E = vars_1D_grid_X(r_E_X_id)%p
            grid_X%loc_r_E = &
                &vars_1D_grid_X(r_E_X_id)%p(X_limits_loc(1):X_limits_loc(2))
            call conv_1D2ND(vars_1D_grid_X(theta_F_X_id),dum_3D)
            grid_X%theta_F = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_grid_X(zeta_F_X_id),dum_3D)
            grid_X%zeta_F = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_grid_X(theta_E_X_id),dum_3D)
            grid_X%theta_E = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_grid_X(zeta_E_X_id),dum_3D)
            grid_X%zeta_E = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
            deallocate(dum_3D)
            call writo('angular size: ('//trim(i2str(n_X(1)))//','//&
                &trim(i2str(n_X(2)))//')')
            call writo('normal size: '//trim(i2str(n_X(3))))
            call lvl_ud(-1)
            
            if (present(grid_X_B)) then
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        grid_X_B => grid_X
                    case (2)                                                    ! HELENA
                        call writo('Creating field-aligned perturbation grid')
                        call lvl_ud(1)
                        n_X_B = vars_1D_grid_X_B(theta_F_X_B_id)%tot_i_max-&
                            &vars_1D_grid_X_B(theta_F_X_B_id)%tot_i_min+1
                        allocate(grid_X_B)
                        ierr = create_grid(grid_X_B,n_X_B,X_limits_loc)
                        CHCKERR('')
                        grid_X_B%r_F = vars_1D_grid_X_B(r_F_X_B_id)%p
                        grid_X_B%loc_r_F = vars_1D_grid_X_B(r_F_X_B_id)%&
                            &p(X_limits_loc(1):X_limits_loc(2))
                        grid_X_B%r_E = vars_1D_grid_X_B(r_E_X_B_id)%p
                        grid_X_B%loc_r_E = vars_1D_grid_X_B(r_E_X_B_id)%&
                            &p(X_limits_loc(1):X_limits_loc(2))
                        call conv_1D2ND(vars_1D_grid_X_B(theta_F_X_B_id),dum_3D)
                        grid_X_B%theta_F = &
                            &dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                        deallocate(dum_3D)
                        call conv_1D2ND(vars_1D_grid_X_B(zeta_F_X_B_id),dum_3D)
                        grid_X_B%zeta_F = &
                            &dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                        deallocate(dum_3D)
                        call conv_1D2ND(vars_1D_grid_X_B(theta_E_X_B_id),dum_3D)
                        grid_X_B%theta_E = &
                            &dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                        deallocate(dum_3D)
                        call conv_1D2ND(vars_1D_grid_X_B(zeta_E_X_B_id),dum_3D)
                        grid_X_B%zeta_E = &
                            &dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                        deallocate(dum_3D)
                        call writo('angular size: ('//trim(i2str(n_X_B(1)))//&
                            &','//trim(i2str(n_X_B(2)))//')')
                        call writo('normal size: '//trim(i2str(n_X_B(3))))
                        call lvl_ud(-1)
                    case default
                        ierr = 1
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        CHCKERR(err_msg)
                end select
            end if
        end function reconstruct_vars_grid_X
        
        ! reconstruct solution grid
        integer function reconstruct_vars_grid_sol(grid) result(ierr)
            use grid_vars, only: create_grid
            
            character(*), parameter :: rout_name = 'reconstruct_vars_grid_sol'
            
            ! input / output
            type(grid_type), intent(inout) :: grid                              ! grid to reconstruct
            
            ! local variables
            integer :: n_sol                                                    ! n of sol grid
            
            ! initialize ierr
            ierr = 0
            
            call writo('Creating solution grid')
            call lvl_ud(1)
            
            ! set n_sol
            n_sol = vars_1D_grid_sol(r_F_sol_id)%tot_i_max(1)-&
                &vars_1D_grid_sol(r_F_sol_id)%tot_i_min(1)+1
            
            ! set up local sol_limits
            sol_limits_loc = [1,n_sol]
            if (present(sol_limits)) sol_limits_loc = sol_limits
            
            ierr = create_grid(grid,n_sol,sol_limits_loc)
            CHCKERR('')
            grid%r_F = vars_1D_grid_sol(r_F_sol_id)%p
            grid%loc_r_F = vars_1D_grid_sol(r_F_sol_id)%&
                &p(sol_limits_loc(1):sol_limits_loc(2))
            grid%r_E = vars_1D_grid_sol(r_E_sol_id)%p
            grid%loc_r_E = vars_1D_grid_sol(r_E_sol_id)%&
                &p(sol_limits_loc(1):sol_limits_loc(2))
            call writo('normal size: '//trim(i2str(n_sol)))
            call lvl_ud(-1)
        end function reconstruct_vars_grid_sol
        
        ! reconstruct equilibrium vars
        integer function reconstruct_vars_eq(grid_eq,eq,met) result(ierr)
            use num_vars, only: eq_style, max_deriv
            use met_vars, only: create_met
            use eq_vars, only: create_eq, max_flux_p_E, max_flux_t_E, &
                &max_flux_p_F, max_flux_t_F
            use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, &
                &ntor, lfreeB
            use HELENA_vars, only: R_H, Z_H, chi_H, flux_p_H, nchi
            
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
            
            met%g_FD = dum_7D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:,:)
            deallocate(dum_7D)
            call conv_1D2ND(vars_1D_eq(h_FD_id),dum_7D)
            met%h_FD = dum_7D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:,:)
            deallocate(dum_7D)
            call conv_1D2ND(vars_1D_eq(jac_FD_id),dum_6D)
            met%jac_FD = dum_6D(:,:,eq_limits_loc(1):eq_limits_loc(2),:,:,:)
            deallocate(dum_6D)
            
            call lvl_ud(-1)
        end function reconstruct_vars_eq
        
        ! reconstruct vectorial perturbation vars
        subroutine reconstruct_vars_X_1(grid_X,X,lim_sec_X)
            use X_vars, only: create_X
            
            ! input / output
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_1_type), intent(inout) :: X                                  ! vectorial perturbation variables
            integer, intent(in), optional :: lim_sec_X(2)                       ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            integer :: id                                                       ! counter
            
            call writo('Setting vectorial perturbation')
            call lvl_ud(1)
            
            call create_X(grid_X,X,lim_sec_X)
            
            ! set up local X_limits
            X_limits_loc = [1,grid_X%n(3)]
            if (present(X_limits)) X_limits_loc = X_limits
            
            ! RE_U_0
            do id = 1,size(RE_U_0_id)
                call conv_1D2ND(vars_1D_X_1(RE_U_0_id(id)),dum_3D)
                X%U_0(:,:,:,id) = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_U_0
            do id = 1,size(IM_U_0_id)
                call conv_1D2ND(vars_1D_X_1(IM_U_0_id(id)),dum_3D)
                X%U_0(:,:,:,id) = X%U_0(:,:,:,id) + &
                    &iu*dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! RE_U_1
            do id = 1,size(RE_U_1_id)
                call conv_1D2ND(vars_1D_X_1(RE_U_1_id(id)),dum_3D)
                X%U_1(:,:,:,id) = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_U_1
            do id = 1,size(IM_U_1_id)
                call conv_1D2ND(vars_1D_X_1(IM_U_1_id(id)),dum_3D)
                X%U_1(:,:,:,id) = X%U_1(:,:,:,id) + &
                    &iu*dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! RE_DU_0
            do id = 1,size(RE_DU_0_id)
                call conv_1D2ND(vars_1D_X_1(RE_DU_0_id(id)),dum_3D)
                X%DU_0(:,:,:,id) = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_DU_0
            do id = 1,size(IM_DU_0_id)
                call conv_1D2ND(vars_1D_X_1(IM_DU_0_id(id)),dum_3D)
                X%DU_0(:,:,:,id) = X%DU_0(:,:,:,id) + &
                    &iu*dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! RE_DU_1
            do id = 1,size(RE_DU_1_id)
                call conv_1D2ND(vars_1D_X_1(RE_DU_1_id(id)),dum_3D)
                X%DU_1(:,:,:,id) = dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            ! IM_DU_1
            do id = 1,size(IM_DU_1_id)
                call conv_1D2ND(vars_1D_X_1(IM_DU_1_id(id)),dum_3D)
                X%DU_1(:,:,:,id) = X%DU_1(:,:,:,id) + &
                    &iu*dum_3D(:,:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_3D)
            enddo
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_X_1
        
        ! reconstruct tensorial perturbation vars
        subroutine reconstruct_vars_X_2(grid_X,X,lim_sec_X)
            use X_vars, only: create_X
            
            ! input / output
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_2_type), intent(inout) :: X                                  ! vectorial perturbation variables
            integer, intent(in), optional :: lim_sec_X(2,2)                     ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            integer :: id                                                       ! counter
            
            call writo('Setting tensorial perturbation')
            call lvl_ud(1)
            
            call create_X(grid_X,X,lim_sec_X)
            
            ! set up local X_limits
            X_limits_loc = [1,grid_X%n(3)]
            if (present(X_limits)) X_limits_loc = X_limits
            
            ! RE_PV_int_0
            do id = 1,size(RE_PV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(RE_PV_int_0_id(id)),dum_2D)
                X%PV_int_0(:,:,id) = &
                    &dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_PV_int_0
            do id = 1,size(IM_PV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(IM_PV_int_0_id(id)),dum_2D)
                X%PV_int_0(:,:,id) = X%PV_int_0(:,:,id) + &
                    &iu*dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_PV_int_1
            do id = 1,size(RE_PV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(RE_PV_int_1_id(id)),dum_2D)
                X%PV_int_1(:,:,id) = &
                    &dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_PV_int_1
            do id = 1,size(IM_PV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(IM_PV_int_1_id(id)),dum_2D)
                X%PV_int_1(:,:,id) = X%PV_int_1(:,:,id) + &
                    &iu*dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_PV_int_2
            do id = 1,size(RE_PV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(RE_PV_int_2_id(id)),dum_2D)
                X%PV_int_2(:,:,id) = &
                    &dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_PV_int_2
            do id = 1,size(IM_PV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(IM_PV_int_2_id(id)),dum_2D)
                X%PV_int_2(:,:,id) = X%PV_int_2(:,:,id) + &
                    &iu*dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_KV_int_0
            do id = 1,size(RE_KV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(RE_KV_int_0_id(id)),dum_2D)
                X%KV_int_0(:,:,id) = &
                    &dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_KV_int_0
            do id = 1,size(IM_KV_int_0_id)
                call conv_1D2ND(vars_1D_X_2(IM_KV_int_0_id(id)),dum_2D)
                X%KV_int_0(:,:,id) = X%KV_int_0(:,:,id) + &
                    &iu*dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_KV_int_1
            do id = 1,size(RE_KV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(RE_KV_int_1_id(id)),dum_2D)
                X%KV_int_1(:,:,id) = &
                    &dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_KV_int_1
            do id = 1,size(IM_KV_int_1_id)
                call conv_1D2ND(vars_1D_X_2(IM_KV_int_1_id(id)),dum_2D)
                X%KV_int_1(:,:,id) = X%KV_int_1(:,:,id) + &
                    &iu*dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! RE_KV_int_2
            do id = 1,size(RE_KV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(RE_KV_int_2_id(id)),dum_2D)
                X%KV_int_2(:,:,id) = &
                    &dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            ! IM_KV_int_2
            do id = 1,size(IM_KV_int_2_id)
                call conv_1D2ND(vars_1D_X_2(IM_KV_int_2_id(id)),dum_2D)
                X%KV_int_2(:,:,id) = X%KV_int_2(:,:,id) + &
                    &iu*dum_2D(:,X_limits_loc(1):X_limits_loc(2))
                deallocate(dum_2D)
            enddo
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_X_2
        
        ! reconstruct solution vars
        subroutine reconstruct_vars_sol(grid_sol,sol,lim_sec_X)
            use sol_vars, only: create_sol
            
            ! input / output
            type(grid_type), intent(in) :: grid_sol                             ! solution grid
            type(sol_type), intent(inout) :: sol                                ! solution variables
            integer, intent(in), optional :: lim_sec_X(2)                       ! limits of m_X (pol. flux) or n_X (tor. flux)
            
            ! local variables
            integer :: n_EV                                                     ! nr. of Eigenvalues
            
            call writo('Setting solution')
            call lvl_ud(1)
            
            ! set n_EV
            n_EV = vars_1D_sol(RE_sol_vec_id)%tot_i_max(3)-&
                &vars_1D_sol(RE_sol_vec_id)%tot_i_min(3)+1
            
            ! create solution
            call create_sol(grid_sol,sol,n_EV,lim_sec_X)
            
            ! set up local sol_limits
            sol_limits_loc = [1,grid_sol%n(3)]
            if (present(sol_limits)) sol_limits_loc = sol_limits
            
            call conv_1D2ND(vars_1D_sol(RE_sol_val_id),dum_1D)
            sol%val = dum_1D
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_sol(IM_sol_val_id),dum_1D)
            sol%val = sol%val + iu*dum_1D
            deallocate(dum_1D)
            call conv_1D2ND(vars_1D_sol(RE_sol_vec_id),dum_3D)
            sol%vec = dum_3D(:,sol_limits_loc(1):sol_limits_loc(2),:)
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_sol(IM_sol_vec_id),dum_3D)
            sol%vec = sol%vec + &
                &iu*dum_3D(:,sol_limits_loc(1):sol_limits_loc(2),:)
            deallocate(dum_3D)
            
            call lvl_ud(-1)
        end subroutine reconstruct_vars_sol
    end function reconstruct_PB3D
end module PB3D_ops
