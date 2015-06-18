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
    use X_vars, only: X_type
    use PB3D_vars, only: PB3D_type
    use HDF5_ops, only: var_1D

    implicit none
    private
    public read_PB3D, reconstruct_PB3D
    
    ! interfaces
    interface conv_1D2ND
        module procedure conv_1D2ND_1D, conv_1D2ND_2D, conv_1D2ND_3D, &
            &conv_1D2ND_4D, &
            !&conv_1D2ND_5D, &
            &conv_1D2ND_6D, conv_1D2ND_7D
    end interface
    
contains
    ! reads PB3D output from user-provided input file
    ! [MPI] only global master
    integer function read_PB3D() result(ierr)
        use num_vars, only: glb_rank, eq_style
        use HDF5_ops, only: read_HDF5_arrs
        use grid_vars, only: create_grid
        use eq_vars, only: create_eq
        use X_vars, only: create_X
        use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B, vars_1D_X
        
        character(*), parameter :: rout_name = 'read_PB3D'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only global master
        if (glb_rank.eq.0) then
            call writo('Reading equilibrium variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_eq,'eq')
            CHCKERR('')
            ierr = set_eq_style()
            CHCKERR('')
            if (eq_style.eq.2) then                                             ! only read eq_B for HELENA
                ierr = read_HDF5_arrs(vars_1D_eq_B,'eq_B')
                CHCKERR('')
            else
                allocate(vars_1D_eq_B(0))                                       ! not needed for other equilibrium styles
            end if
            call lvl_ud(-1)
            call writo('Equilibrium variables read')
            
            call writo('Reading perturbation variables')
            call lvl_ud(1)
            ierr = read_HDF5_arrs(vars_1D_X,'X')
            err_msg = 'Maybe you wanted to use Richardson extrapolation? Not &
                &yet implemented.'
            CHCKERR(err_msg)
            call lvl_ud(-1)
            call writo('Perturbation variables read')
        end if
    contains
        integer function set_eq_style() result(ierr)
            use num_vars, only: eq_style, rho_style
            
            ! local variables
            integer :: misc_eq_id                                               ! index of misc_eq
            real(dp), allocatable :: dum_1D(:)                                  ! dummy variable
            
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
    
    ! Reconstructs the PB3D variables
    ! Note that the equilibrium grid  is different for the different equilibrium
    ! styles:
    !   - VMEC: Field-aligned grid on which  the PB3D coefficients KV_i and PV_i
    !   are calculated.
    !   - HELENA:  Grid  returned  by HELENA  which,  though   technically  also
    !   field-aligned, has to be interpolated to yield the field-aligned grid on
    !   which the PB3D coefficients KV_i and PV_i are calcualted.
    ! This implies that for HELENA the field-aligned grid has to be interpolated
    ! still,  but at  the  same time  every  other grid  can  equally easily  be
    ! interpolated. For VMEC, for every grid one has to start from the basis and
    ! run the routines from calc_eq.
    integer function reconstruct_PB3D(PB3D,grid_eq_B) result(ierr)
        use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B, vars_1D_X, &
            &min_PB3D_version
        use MPI_ops, only: split_MPI_POST
        use grid_vars, only: create_grid
        use met_vars, only: create_met
        use X_vars, only: create_X, &
            &min_r_X, max_r_X, min_n_X, max_n_X, min_m_X, max_m_X
        use eq_vars, only: create_eq, R_0, pres_0, B_0, psi_0, rho_0, T_0, &
            &vac_perm, max_flux_p_E, max_flux_t_E, max_flux_p_F, max_flux_t_F
        use num_vars, only: eq_style, max_deriv, use_pol_flux_E, use_pol_flux_F
        use VMEC, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mpol, ntor, &
            &nfp, lfreeB, lasym
        use HELENA, only: R_H, Z_H, chi_H, flux_p_H, ias, nchi
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D'
        
        ! input  / output
        type(PB3D_type), intent(inout) :: PB3D                                  ! PB3D for which to do postprocessing
        type(grid_type), intent(inout), optional :: grid_eq_B                   ! optional field-aligned grid for HELENA
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
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
        integer :: i_lim_eq(2), i_lim_X(2)                                      ! i_lim of equilibrium and perturbation
        integer :: n_X, n_eq(3), n_eq_B(3)                                      ! n of X, eq and eq_B grid
        integer :: n_r_eq                                                       ! n_r of eq grid
        real(dp), allocatable :: dum_1D(:), dum_2D(:,:), dum_3D(:,:,:)          ! dummy variables
        real(dp), allocatable :: dum_4D(:,:,:,:)                                ! dummy variables
        !real(dp), allocatable :: dum_5D(:,:,:,:,:)                              ! dummy variables
        real(dp), allocatable :: dum_6D(:,:,:,:,:,:), dum_7D(:,:,:,:,:,:,:)     ! dummy variables
        real(dp), parameter :: tol_version = 1.E-8_dp                           ! tolerance for version control
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Prepare variable indices')
        call lvl_ud(1)
        
        ! get 1D X indices
        ierr = retrieve_var_1D_id(vars_1D_X,'r_F',r_F_X_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'r_E',r_E_X_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'misc_X',misc_X_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'RE_U_0',RE_U_0_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'IM_U_0',IM_U_0_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'RE_U_1',RE_U_1_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'IM_U_1',IM_U_1_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'RE_DU_0',RE_DU_0_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'IM_DU_0',IM_DU_0_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'RE_DU_1',RE_DU_1_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'IM_DU_1',IM_DU_1_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'RE_X_val',RE_X_val_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'IM_X_val',IM_X_val_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'RE_X_vec',RE_X_vec_id)
        CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_X,'IM_X_vec',IM_X_vec_id)
        CHCKERR('')
        
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
        
        ! get 1D eq indices
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
        
        ! Set up particular 1D indices, depending on equilibrium style
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
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
                ierr = retrieve_var_1D_id(vars_1D_eq,'misc_eq_V',misc_eq_V_id)
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
            case (2)                                                            ! HELENA
                ierr = retrieve_var_1D_id(vars_1D_eq,'R_H',R_H_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq,'Z_H',Z_H_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq,'chi_H',chi_H_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq,'flux_p_H',flux_p_H_id)
                CHCKERR('')
                ierr = retrieve_var_1D_id(vars_1D_eq,'misc_eq_H',misc_eq_H_id)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        call lvl_ud(-1)
        call writo('Indices prepared')
        
        ! user output
        call writo('Dividing normal range over processes')
        call lvl_ud(1)
        
        ! split normal grids
        ierr = split_MPI_POST(vars_1D_eq(r_F_eq_id)%p,vars_1D_X(r_F_X_id)%p,&
            &i_lim_eq,i_lim_X)
        CHCKERR('')
        
        call lvl_ud(-1)
        call writo('Normal range divided')
        
        call writo('Reconstructing variables')
        call lvl_ud(1)
        
        call writo('Setting grids')
        call lvl_ud(1)
        
        call writo('Creating perturbation grid')
        n_X = vars_1D_X(r_F_X_id)%tot_i_max(1)-&
            &vars_1D_X(r_F_X_id)%tot_i_min(1)+1
        ierr = create_grid(PB3D%grid_X,n_X,i_lim_X)
        CHCKERR('')
        PB3D%grid_X%r_F = vars_1D_X(r_F_X_id)%p
        PB3D%grid_X%grp_r_F = &
            &vars_1D_X(r_F_X_id)%p(i_lim_X(1):i_lim_X(2))
        PB3D%grid_X%r_E = vars_1D_X(r_E_X_id)%p
        PB3D%grid_X%grp_r_E = &
            &vars_1D_X(r_E_X_id)%p(i_lim_X(1):i_lim_X(2))
        
        if (present(grid_eq_B)) then
            call writo('Creating field-aligned equilibrium grid')
            n_eq_B = vars_1D_eq_B(theta_F_B_id)%tot_i_max-&
                &vars_1D_eq_B(theta_F_B_id)%tot_i_min+1
            ierr = create_grid(grid_eq_B,n_eq_B,i_lim_eq)
            CHCKERR('')
            grid_eq_B%r_F = vars_1D_eq_B(r_F_eq_B_id)%p
            grid_eq_B%grp_r_F = &
                &vars_1D_eq_B(r_F_eq_B_id)%p(i_lim_eq(1):i_lim_eq(2))
            grid_eq_B%r_E = vars_1D_eq_B(r_E_eq_B_id)%p
            grid_eq_B%grp_r_E = &
                &vars_1D_eq_B(r_E_eq_B_id)%p(i_lim_eq(1):i_lim_eq(2))
            call conv_1D2ND(vars_1D_eq_B(theta_F_B_id),dum_3D)
            grid_eq_B%theta_F = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq_B(zeta_F_B_id),dum_3D)
            grid_eq_B%zeta_F = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq_B(theta_E_B_id),dum_3D)
            grid_eq_B%theta_E = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
            deallocate(dum_3D)
            call conv_1D2ND(vars_1D_eq_B(zeta_E_B_id),dum_3D)
            grid_eq_B%zeta_E = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
            deallocate(dum_3D)
        end if
        
        call writo('Creating equilibrium grid for output tables')
        n_eq = vars_1D_eq(theta_F_id)%tot_i_max-&
            &vars_1D_eq(theta_F_id)%tot_i_min+1
        ierr = create_grid(PB3D%grid_eq,n_eq,i_lim_eq)
        CHCKERR('')
        PB3D%grid_eq%r_F = vars_1D_eq(r_F_eq_id)%p
        PB3D%grid_eq%grp_r_F = &
            &vars_1D_eq(r_F_eq_id)%p(i_lim_eq(1):i_lim_eq(2))
        PB3D%grid_eq%r_E = vars_1D_eq(r_E_eq_id)%p
        PB3D%grid_eq%grp_r_E = &
            &vars_1D_eq(r_E_eq_id)%p(i_lim_eq(1):i_lim_eq(2))
        call conv_1D2ND(vars_1D_eq(theta_F_id),dum_3D)
        PB3D%grid_eq%theta_F = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(zeta_F_id),dum_3D)
        PB3D%grid_eq%zeta_F = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(theta_E_id),dum_3D)
        PB3D%grid_eq%theta_E = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(zeta_E_id),dum_3D)
        PB3D%grid_eq%zeta_E = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        
        call lvl_ud(-1)
        
        call writo('Setting equilibrium')
        ierr = create_eq(PB3D%grid_eq,PB3D%eq)
        CHCKERR('')
        
        call conv_1D2ND(vars_1D_eq(pres_FD_id),dum_2D)
        PB3D%eq%pres_FD = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_2D)
        call conv_1D2ND(vars_1D_eq(q_saf_FD_id),dum_2D)
        PB3D%eq%q_saf_FD = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_2D)
        call conv_1D2ND(vars_1D_eq(rot_t_FD_id),dum_2D)
        PB3D%eq%rot_t_FD = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_2D)
        call conv_1D2ND(vars_1D_eq(flux_p_FD_id),dum_2D)
        PB3D%eq%flux_p_FD = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_2D)
        call conv_1D2ND(vars_1D_eq(flux_t_FD_id),dum_2D)
        PB3D%eq%flux_t_FD = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_2D)
        call conv_1D2ND(vars_1D_eq(rho_id),dum_1D)
        PB3D%eq%rho = dum_1D(i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_1D)
        call conv_1D2ND(vars_1D_eq(S_id),dum_3D)
        PB3D%eq%S = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(kappa_n_id),dum_3D)
        PB3D%eq%kappa_n = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(kappa_g_id),dum_3D)
        PB3D%eq%kappa_g = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(sigma_id),dum_3D)
        PB3D%eq%sigma = dum_3D(:,:,i_lim_eq(1):i_lim_eq(2))
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_eq(g_FD_id),dum_7D)
        
        ! Set up particular variables, depending on equilibrium style
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
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
                PB3D%eq%pres_E = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
                deallocate(dum_2D)
                call conv_1D2ND(vars_1D_eq(q_saf_E_id),dum_2D)
                PB3D%eq%q_saf_E = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
                deallocate(dum_2D)
                call conv_1D2ND(vars_1D_eq(rot_t_E_id),dum_2D)
                PB3D%eq%rot_t_E = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
                deallocate(dum_2D)
                call conv_1D2ND(vars_1D_eq(flux_p_E_id),dum_2D)
                PB3D%eq%flux_p_E = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
                deallocate(dum_2D)
                call conv_1D2ND(vars_1D_eq(flux_t_E_id),dum_2D)
                PB3D%eq%flux_t_E = dum_2D(i_lim_eq(1):i_lim_eq(2),:)
                deallocate(dum_2D)
            case (2)                                                            ! HELENA
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
        
        PB3D%met%g_FD = dum_7D(:,:,i_lim_eq(1):i_lim_eq(2),:,:,:,:)
        deallocate(dum_7D)
        call conv_1D2ND(vars_1D_eq(h_FD_id),dum_7D)
        PB3D%met%h_FD = dum_7D(:,:,i_lim_eq(1):i_lim_eq(2),:,:,:,:)
        deallocate(dum_7D)
        call conv_1D2ND(vars_1D_eq(jac_FD_id),dum_6D)
        PB3D%met%jac_FD = dum_6D(:,:,i_lim_eq(1):i_lim_eq(2),:,:,:)
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
        if (dum_1D(16).gt.0) use_pol_flux_E = .true.
        if (dum_1D(17).gt.0) use_pol_flux_F = .true.
        deallocate(dum_1D)
        
        call writo('Setting perturbation')
        call lvl_ud(1)
        
        call conv_1D2ND(vars_1D_X(misc_X_id),dum_1D)
        min_r_X = dum_1D(1)
        max_r_X = dum_1D(2)
        min_n_X = nint(dum_1D(3))
        max_n_X = nint(dum_1D(4))
        min_m_X = nint(dum_1D(5))
        max_m_X = nint(dum_1D(6))
        deallocate(dum_1D)
        
        ! user output
        call writo('The PB3D output')
        call lvl_ud(1)
        if (use_pol_flux_F) then
            call writo('has modes n = '//trim(i2str(min_n_X))//&
                &' and m = '//trim(i2str(min_m_X))//'..'//trim(i2str(max_m_X)))
            call writo('works with the poloidal flux as normal coordinate')
        else
            call writo('has modes m = '//trim(i2str(min_m_X))//&
                &' and n='//trim(i2str(min_n_X))//'..'//trim(i2str(max_n_X)))
            call writo('works with the toroidal flux as normal coordinate')
        end if
        call writo('and its normal boundaries are '//trim(r2strt(min_r_X))//&
            &' and '//trim(r2strt(max_r_X)))
        call lvl_ud(-1)
        
        call create_X(PB3D%grid_eq,PB3D%X)
        
        call conv_1D2ND(vars_1D_X(RE_U_0_id),dum_4D)
        PB3D%X%U_0 = dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(IM_U_0_id),dum_4D)
        PB3D%X%U_0 = PB3D%X%U_0 + iu*dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(RE_U_1_id),dum_4D)
        PB3D%X%U_1 = dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(IM_U_1_id),dum_4D)
        PB3D%X%U_1 = PB3D%X%U_1 + iu*dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(RE_DU_0_id),dum_4D)
        PB3D%X%DU_0 = dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(IM_DU_0_id),dum_4D)
        PB3D%X%DU_0 = PB3D%X%DU_0 + iu*dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(RE_DU_1_id),dum_4D)
        PB3D%X%DU_1 = dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        call conv_1D2ND(vars_1D_X(IM_DU_1_id),dum_4D)
        PB3D%X%DU_1 = PB3D%X%DU_1 + iu*dum_4D(:,:,i_lim_eq(1):i_lim_eq(2),:)
        deallocate(dum_4D)
        
        allocate(PB3D%X%val(vars_1D_X(RE_X_val_id)%tot_i_min(1):&
            &vars_1D_X(RE_X_val_id)%tot_i_max(1)))
        call conv_1D2ND(vars_1D_X(RE_X_val_id),dum_1D)
        PB3D%X%val = dum_1D
        deallocate(dum_1D)
        call conv_1D2ND(vars_1D_X(IM_X_val_id),dum_1D)
        PB3D%X%val = PB3D%X%val + iu*dum_1D
        deallocate(dum_1D)
        allocate(PB3D%X%vec(vars_1D_X(RE_X_vec_id)%tot_i_min(1):&
            &vars_1D_X(RE_X_vec_id)%tot_i_max(1),&
            &1:i_lim_X(2)-i_lim_X(1)+1,&
            &vars_1D_X(RE_X_vec_id)%tot_i_min(3):&
            &vars_1D_X(RE_X_vec_id)%tot_i_max(3)))
        call conv_1D2ND(vars_1D_X(RE_X_vec_id),dum_3D)
        PB3D%X%vec = dum_3D(:,i_lim_X(1):i_lim_X(2),:)
        deallocate(dum_3D)
        call conv_1D2ND(vars_1D_X(IM_X_vec_id),dum_3D)
        PB3D%X%vec = PB3D%X%vec + iu*dum_3D(:,i_lim_X(1):i_lim_X(2),:)
        deallocate(dum_3D)
        
        call lvl_ud(-1)
        
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
