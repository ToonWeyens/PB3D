!------------------------------------------------------------------------------!
!   Operations on PB3D output.                                                 !
!   Note: If you have parallel jobs, if the reconstruction routines are called !
!   for  multiple equilibrium  jobs,  this can  only be  done  after the  last !
!   equilibrium job is finished. There are no checks for this.                 !
!------------------------------------------------------------------------------!
module PB3D_ops
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, max_name_ln, iu, min_PB3D_version
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type
    use HDF5_vars, only: dealloc_var_1D, var_1D_type
    use PB3D_utilities, only: setup_rich_id, setup_par_id

    implicit none
    private
    public reconstruct_PB3D_in, reconstruct_PB3D_grid, reconstruct_PB3D_eq_1, &
        &reconstruct_PB3D_eq_2, reconstruct_PB3D_X_1, reconstruct_PB3D_X_2, &
        &reconstruct_PB3D_sol, get_PB3D_grid_size
    
    ! global variables
    real(dp), allocatable :: dum_1D(:)                                          ! dummy variables
    real(dp), allocatable :: dum_2D(:,:)                                        ! dummy variables
    real(dp), allocatable :: dum_3D(:,:,:)                                      ! dummy variables
    real(dp), allocatable :: dum_4D(:,:,:,:)                                    ! dummy variables
    real(dp), allocatable :: dum_6D(:,:,:,:,:,:)                                ! dummy variables
    real(dp), allocatable :: dum_7D(:,:,:,:,:,:,:)                              ! dummy variables
    
contains
    ! Reconstructs the input variables from PB3D output.
    integer function reconstruct_PB3D_in(data_name) result(ierr)
        use num_vars, only: eq_style, rho_style, use_pol_flux_E, max_deriv, &
            &use_pol_flux_F, use_normalization, norm_disc_prec_eq, PB3D_name, &
            &norm_disc_prec_X, norm_style, U_style, X_style, prog_style, &
            &matrix_SLEPC_style, BC_style, EV_style, norm_disc_prec_sol, &
            &EV_BC, magn_int_style, K_style, debug_version
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm, &
            &max_flux_E, max_flux_F
        use grid_vars, onLy: n_r_in, n_r_eq, n_r_sol
        use X_vars, only: min_r_sol, max_r_sol, min_sec_X, max_sec_X, prim_X, &
            &n_mod_X
        use sol_vars, only: alpha
        use HELENA_vars, only: chi_H, flux_p_H, flux_t_H, R_H, Z_H, nchi, ias, &
            &q_saf_H, rot_t_H, pres_H, RBphi_H
        use VMEC_vars, only: is_freeb_V, mnmax_V, mpol_V, ntor_V, is_asym_V, &
            &gam_V, R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mn_V, rot_t_V, &
            &q_saf_V, pres_V, flux_t_V, flux_p_V, nfp_V
        use HELENA_vars, only: h_H_11, h_H_12, h_H_33
#if ldebug
        use VMEC_vars, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, &
            &jac_V_s
#endif
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_in'
            
        ! input / output
        character(len=*), intent(in) :: data_name                               ! name of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(var_1D_type) :: var_1D                                             ! 1D variable
        real(dp), parameter :: tol_version = 1.E-4_dp                           ! tolerance for version control
        real(dp) :: PB3D_version                                                ! version of PB3D variable read
        logical :: debug_version_PB3D                                           ! debug version of in
        character(len=max_str_ln) :: use_debug_str(2)                           ! using debug or not
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Reconstructing input variables from PB3D output')
        call lvl_ud(1)
        
        ! misc_in
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'misc_in')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        PB3D_version = dum_1D(1)
        eq_style = nint(dum_1D(2))
        rho_style = nint(dum_1D(3))
        R_0 = dum_1D(4)
        pres_0 = dum_1D(5)
        B_0 = dum_1D(6)
        psi_0 = dum_1D(7)
        rho_0 = dum_1D(8)
        T_0 = dum_1D(9)
        vac_perm = dum_1D(10)
        use_pol_flux_E = .false.
        use_pol_flux_F = .false.
        use_normalization = .false.
        if (dum_1D(11).gt.0) use_pol_flux_E = .true.
        if (dum_1D(12).gt.0) use_pol_flux_F = .true.
        if (dum_1D(13).gt.0) use_normalization = .true.
        norm_disc_prec_eq = nint(dum_1D(14))
        n_r_in = nint(dum_1D(15))
        n_r_eq = nint(dum_1D(16))
        n_r_sol = nint(dum_1D(17))
        max_flux_E = dum_1D(18)
        max_flux_F = dum_1D(19)
        debug_version_PB3D = .false.
        if (dum_1D(20).gt.0) debug_version_PB3D = .true.
        call dealloc_var_1D(var_1D)
        
        ! tests
        select case (prog_style)
            case (1)                                                            ! PB3D
                ! do nothing
            case (2)                                                            ! POST
                call writo('Run tests')
                call lvl_ud(1)
                
                call writo('PB3D version '//trim(r2strt(PB3D_version)))
                if (PB3D_version.lt.min_PB3D_version*(1-tol_version)) then
                    ierr = 1
                    err_msg = 'Need at least PB3D version '//&
                        &trim(r2strt(min_PB3D_version))
                    CHCKERR(err_msg)
                end if
                
                if (debug_version_PB3D) call writo('debug version')
                
                if (debug_version_PB3D.neqv.debug_version) then
                    ierr = 1
                    if (debug_version_PB3D) then
                        use_debug_str(1) = 'uses debug version'
                    else
                        use_debug_str(1) = 'uses release version'
                    end if
                    if (debug_version) then
                        use_debug_str(2) = 'uses debug version'
                    else
                        use_debug_str(2) = 'uses release version'
                    end if
                    call writo('The PB3D output '//trim(use_debug_str(1))//&
                        &' but POST '//trim(use_debug_str(2)))
                    err_msg = 'Need to use debug version for both, or not for &
                        &both'
                    CHCKERR(err_msg)
                end if
                
                call lvl_ud(-1)
        end select
        
        ! variables depending on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! misc_in_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'misc_in_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                is_asym_V = .false.
                is_freeb_V = .false.
                if (dum_1D(1).gt.0) is_asym_V = .true.
                if (dum_1D(2).gt.0) is_freeb_V = .true.
                mnmax_V = nint(dum_1D(3))
                mpol_V = nint(dum_1D(4))
                ntor_V = nint(dum_1D(5))
                nfp_V = nint(dum_1D(6))
                gam_V = dum_1D(7)
                call dealloc_var_1D(var_1D)
                
                ! flux_t_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_t_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(flux_t_V(n_r_eq,0:max_deriv+2))
                flux_t_V = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! flux_p_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_p_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(flux_p_V(n_r_eq,0:max_deriv+2))
                flux_p_V = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! pres_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'pres_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(pres_V(n_r_eq,0:max_deriv+1))
                pres_V = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! rot_t_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rot_t_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(rot_t_V(n_r_eq,0:max_deriv+1))
                rot_t_V = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! q_saf_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'q_saf_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(q_saf_V(n_r_eq,0:max_deriv+1))
                q_saf_V = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! mn_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'mn_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(mn_V(mnmax_V,2))
                mn_V = nint(dum_2D)
                call dealloc_var_1D(var_1D)
                
                ! RZL_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RZL_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_4D)
                allocate(R_V_c(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(R_V_s(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(Z_V_c(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(Z_V_s(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(L_V_c(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(L_V_s(mnmax_V,n_r_eq,0:max_deriv+1))
                R_V_c = dum_4D(:,:,:,1)
                R_V_s = dum_4D(:,:,:,2)
                Z_V_c = dum_4D(:,:,:,3)
                Z_V_s = dum_4D(:,:,:,4)
                L_V_c = dum_4D(:,:,:,5)
                L_V_s = dum_4D(:,:,:,6)
                call dealloc_var_1D(var_1D)
                
#if ldebug
                ! B_V_sub
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'B_V_sub')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_4D)
                allocate(B_V_sub_c(mnmax_V,n_r_eq,3))
                allocate(B_V_sub_s(mnmax_V,n_r_eq,3))
                B_V_sub_c = dum_4D(:,:,:,1)
                B_V_sub_s = dum_4D(:,:,:,2)
                call dealloc_var_1D(var_1D)
                
                ! B_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'B_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(B_V_c(mnmax_V,n_r_eq))
                allocate(B_V_s(mnmax_V,n_r_eq))
                B_V_c = dum_3D(:,:,1)
                B_V_s = dum_3D(:,:,2)
                call dealloc_var_1D(var_1D)
                
                ! jac_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'jac_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(jac_V_c(mnmax_V,n_r_eq))
                allocate(jac_V_s(mnmax_V,n_r_eq))
                jac_V_c = dum_3D(:,:,1)
                jac_V_s = dum_3D(:,:,2)
                call dealloc_var_1D(var_1D)
#endif
            case (2)                                                            ! HELENA
                ! misc_in_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'misc_in_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                ias = nint(dum_1D(1))
                nchi= nint(dum_1D(2))
                call dealloc_var_1D(var_1D)
                
                ! RZ_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RZ_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(R_H(nchi,n_r_eq))
                allocate(Z_H(nchi,n_r_eq))
                R_H = dum_3D(:,:,1)
                Z_H = dum_3D(:,:,2)
                call dealloc_var_1D(var_1D)
                
                ! chi_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'chi_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(chi_H(nchi))
                chi_H = dum_1D
                call dealloc_var_1D(var_1D)
                
                ! flux_p_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_p_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(flux_p_H(n_r_eq,0:max_deriv+1))
                flux_p_H = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! flux_t_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_t_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(flux_t_H(n_r_eq,0:max_deriv+1))
                flux_t_H = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! q_saf_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'q_saf_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(q_saf_H(n_r_eq,0:max_deriv+1))
                q_saf_H = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! rot_t_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rot_t_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(rot_t_H(n_r_eq,0:max_deriv+1))
                rot_t_H = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! pres_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'pres_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(pres_H(n_r_eq,0:max_deriv+1))
                pres_H = dum_2D
                call dealloc_var_1D(var_1D)
                
                ! RBphi_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RBphi_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(RBphi_H(n_r_eq))
                RBphi_H = dum_1D
                call dealloc_var_1D(var_1D)
                
                ! h_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'h_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(h_H_11(nchi,n_r_eq))
                allocate(h_H_12(nchi,n_r_eq))
                allocate(h_H_33(nchi,n_r_eq))
                h_H_11 = dum_3D(:,:,1)
                h_H_12 = dum_3D(:,:,2)
                h_H_33 = dum_3D(:,:,3)
                call dealloc_var_1D(var_1D)
        end select
        
        ! misc_X
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'misc_X')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        prim_X = nint(dum_1D(1))
        n_mod_X = nint(dum_1D(2))
        min_sec_X = nint(dum_1D(3))
        max_sec_X = nint(dum_1D(4))
        norm_disc_prec_X = nint(dum_1D(5))
        norm_style = nint(dum_1D(6))
        U_style = nint(dum_1D(7))
        X_style = nint(dum_1D(8))
        matrix_SLEPC_style = nint(dum_1D(9))
        select case(prog_style)
            case (1)                                                            ! PB3D
                magn_int_style = nint(dum_1D(10))
            case (2)                                                            ! POST
                magn_int_style = 1                                              ! integration is done in volume, currently only with trapezoidal rule
        end select
        K_style = nint(dum_1D(11))
        call dealloc_var_1D(var_1D)
        
        ! misc_sol
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'misc_sol')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        min_r_sol = dum_1D(1)
        max_r_sol = dum_1D(2)
        alpha = dum_1D(3)
        norm_disc_prec_sol = nint(dum_1D(4))
        BC_style(1) = nint(dum_1D(5))
        BC_style(2) = nint(dum_1D(6))
        EV_style = nint(dum_1D(7))
        EV_BC = dum_1D(8)
        call dealloc_var_1D(var_1D)
        
        ! user output
        call lvl_ud(-1)
        call writo('Input variables from PB3D output reconstructed')
        
        ! clean up
        deallocate(dum_1D)
        deallocate(dum_3D)
    end function reconstruct_PB3D_in
    
    ! Reconstructs grid variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name  if it  is  >  0.
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    ! Furthermore,  using "lim_pos",  you can  obtain a  subset of  the data  by
    ! directly passing its limits to the underlying HDF5 routine. This refers to
    ! the position dimensions only. If provided,  the normal limits of a divided
    ! grid refer to the subset, as in "copy_grid".
    ! Note: "grid_" is added in front the data_name.
    ! Note:  By providing  lim_pos equal  to 0  in the  angular dimensions,  the
    ! angular part can be discarded when reconstructing the grid.
    integer function reconstruct_PB3D_grid(grid,data_name,rich_lvl,tot_rich,&
        &lim_pos,grid_limits) result(ierr)
        
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_grid'
        
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid 
        character(len=*), intent(in) :: data_name                               ! name of grid
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: lim_pos(3,2)                           ! position limits of subset of data
        integer, intent(in), optional :: grid_limits(2)                         ! i_limit of grid
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        integer :: lim_mem(3,2)                                                 ! memory limits for variables
        integer :: n(3)                                                         ! total n
        integer :: id                                                           ! counter
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: par_id_mem(2)                                                ! parallel indices (start, end) in memory (stride 1)
        integer :: par_lim(2)                                                   ! parallel limits
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        
        ! setup rich_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        
        ! get total grid size
        ierr = get_PB3D_grid_size(n,data_name,rich_lvl,tot_rich)
        CHCKERR('')
        
        ! possibly change n to user-specified
        if (present(lim_pos)) then
            if (lim_pos(1,1).ge.0 .and. lim_pos(1,2).ge.0) n(1) = &
                &lim_pos(1,2)-lim_pos(1,1)+1
            if (lim_pos(2,1).ge.0 .and. lim_pos(2,2).ge.0) n(2) = &
                &lim_pos(2,2)-lim_pos(2,1)+1
            if (lim_pos(3,1).ge.0 .and. lim_pos(3,2).ge.0) n(3) = &
                &lim_pos(3,2)-lim_pos(3,1)+1
            where (lim_pos(:,1).eq.lim_pos(:,2) .and. lim_pos(:,1).eq.[0,0,0]) &
                &n = 0
        end if
        
        ! set up parallel limits
        par_lim = [1,n(1)]
        if (present(lim_pos)) then
            where (lim_pos(1,:).ge.0) par_lim = lim_pos(1,:)
        end if
        
        ! create grid
        ierr = grid%init(n,i_lim=grid_limits)
        CHCKERR('')
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id 
            par_id = setup_par_id(grid,rich_lvl_loc,id,tot_rich=tot_rich,&
                &par_lim=par_lim,par_id_mem=par_id_mem)
            
            ! set up local limits for HDF5 reconstruction of full vars
            lim_mem(1,:) = par_id_mem
            lim_mem(2,:) = [1,grid%n(2)]
            lim_mem(3,:) = [1,grid%n(3)]
            if (present(lim_pos)) then
                lim_mem(2,:) = lim_pos(2,:)
                lim_mem(3,:) = lim_pos(3,:)
            end if
            
            ! r_F
            ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//&
                &trim(data_name),'r_F',rich_lvl=id,lim_loc=lim_mem(3:3,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_1D)
            grid%r_F = dum_1D
            call dealloc_var_1D(var_1D)
            
            ! r_E
            ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//&
                &trim(data_name),'r_E',rich_lvl=id,lim_loc=lim_mem(3:3,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_1D)
            grid%r_E = dum_1D
            call dealloc_var_1D(var_1D)
            
            ! loc_r_F
            grid%loc_r_F = grid%r_F(grid%i_min:grid%i_max)
            
            ! loc_r_E
            grid%loc_r_E = grid%r_E(grid%i_min:grid%i_max)
            
            ! overwrite local limits for HDF5 reconstruction of divided vars
            lim_mem(1,:) = par_id_mem
            lim_mem(2,:) = [1,grid%n(2)]
            lim_mem(3,:) = [grid%i_min,grid%i_max]
            if (present(lim_pos)) then
                lim_mem(2,:) = lim_pos(2,:)
                lim_mem(3,:) = lim_mem(3,:) + lim_pos(3,1) - 1                  ! take into account the grid limits (relative to position subset)
            end if
            
            ! only for 3D grids
            if (product(grid%n(1:2)).ne.0) then
                ! theta_F
                ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'theta_F',rich_lvl=id,lim_loc=lim_mem)
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                grid%theta_F(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
                call dealloc_var_1D(var_1D)
                
                ! theta_E
                ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'theta_E',rich_lvl=id,lim_loc=lim_mem)
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                grid%theta_E(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
                call dealloc_var_1D(var_1D)
                
                ! zeta_F
                ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'zeta_F',rich_lvl=id,lim_loc=lim_mem)
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                grid%zeta_F(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
                call dealloc_var_1D(var_1D)
                
                ! zeta_E
                ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'zeta_E',rich_lvl=id,lim_loc=lim_mem)
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                grid%zeta_E(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
                call dealloc_var_1D(var_1D)
            end if
        end do
        
        ! clean up
        deallocate(dum_1D)
        if (product(grid%n(1:2)).ne.0) deallocate(dum_3D)
    end function reconstruct_PB3D_grid
    
    ! Reconstructs the equilibrium variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Furthermore,  using "lim_pos",  you can  obtain a  subset of  the data  by
    ! directly passing its limits to the underlying HDF5 routine. This refers to
    ! the position dimensions only. If provided,  the normal limits of a divided
    ! grid refer to the subset, as in "copy_grid".
    integer function reconstruct_PB3D_eq_1(grid_eq,eq,data_name,lim_pos) &
        &result(ierr)                                                           ! flux version
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_eq_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid 
        type(eq_1_type), intent(inout), optional :: eq                          ! flux equilibrium
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: lim_pos(1,2)                           ! position limits of subset of data
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        integer :: lim_mem(2,2)                                                 ! memory limits for variables
        
        ! initialize ierr
        ierr = 0
        
        ! prepare
        
        ! create equilibrium
        call eq%init(grid_eq,setup_E=.false.,setup_F=.true.)
        
        ! set up local limits for HDF5 reconstruction
        lim_mem(1,:) = [grid_eq%i_min,grid_eq%i_max]
        lim_mem(2,:) = [0,-1]                                                   ! all derivatives, starting from 0
        if (present(lim_pos)) then
            lim_mem(1,:) = lim_mem(1,:) + lim_pos(1,1) - 1                      ! take into account the grid limits (relative to position subset)
        end if
        
        ! restore variables
        
        ! pres_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'pres_FD',&
            &lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%pres_FD = dum_2D
        call dealloc_var_1D(var_1D)
        
        ! q_saf_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'q_saf_FD',&
            &lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%q_saf_FD = dum_2D
        call dealloc_var_1D(var_1D)
        
        ! rot_t_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rot_t_FD',&
            &lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%rot_t_FD = dum_2D
        call dealloc_var_1D(var_1D)
        
        ! flux_p_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'flux_p_FD',&
            &lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%flux_p_FD = dum_2D
        call dealloc_var_1D(var_1D)
        
        ! flux_t_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'flux_t_FD',&
            &lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%flux_t_FD = dum_2D
        call dealloc_var_1D(var_1D)
        
        ! rho
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rho',&
            &lim_loc=lim_mem(1:1,:))
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        eq%rho = dum_1D
        call dealloc_var_1D(var_1D)
        
        ! clean up
        deallocate(dum_1D)
        deallocate(dum_2D)
    end function reconstruct_PB3D_eq_1
    
    ! Reconstructs the equilibrium variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name if it is > 0.
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    ! Furthermore,  using "lim_pos",  you can  obtain a  subset of  the data  by
    ! directly passing its limits to the underlying HDF5 routine. This refers to
    ! the position dimensions only. If provided,  the normal limits of a divided
    ! grid refer to the subset, as in "copy_grid".
    integer function reconstruct_PB3D_eq_2(grid_eq,eq,data_name,rich_lvl,&
        &tot_rich,lim_pos) result(ierr)                                         ! metric version
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_eq_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid 
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: lim_pos(3,2)                           ! position limits of subset of data
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        integer :: lim_mem(7,2)                                                 ! memory limits for variables
        integer :: id                                                           ! counter
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: par_id_mem(2)                                                ! parallel indices (start, end) in memory (stride 1)
        integer :: par_lim(2)                                                   ! parallel limits
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        
        ! setup rich_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        
        ! set up parallel limits
        par_lim = [1,grid_eq%n(1)]
        if (present(lim_pos)) then
            where (lim_pos(1,:).ge.0) par_lim = lim_pos(1,:)
        end if
        
        ! create equilibrium
        call eq%init(grid_eq,setup_E=.false.,setup_F=.true.)
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id
            par_id = setup_par_id(grid_eq,rich_lvl_loc,id,tot_rich=tot_rich,&
                &par_lim=par_lim,par_id_mem=par_id_mem)
            
            ! set up local limits for HDF5 reconstruction
            lim_mem(1,:) = par_id_mem
            lim_mem(2,:) = [1,grid_eq%n(2)]
            lim_mem(3,:) = [grid_eq%i_min,grid_eq%i_max]
            lim_mem(4,:) = [-1,-1]
            lim_mem(5,:) = [-1,-1]
            lim_mem(6,:) = [-1,-1]
            lim_mem(7,:) = [-1,-1]
            if (present(lim_pos)) then
                lim_mem(2,:) = lim_pos(2,:)
                lim_mem(3,:) = lim_mem(3,:) + lim_pos(3,1) - 1                  ! take into account the grid limits (relative to position subset)
            end if
            
            ! g_FD
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'g_FD',&
                &rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_7D)
            eq%g_FD(par_id(1):par_id(2):par_id(3),:,:,:,:,:,:) = dum_7D
            call dealloc_var_1D(var_1D)
            
            ! h_FD
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'h_FD',&
                &rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_7D)
            eq%h_FD(par_id(1):par_id(2):par_id(3),:,:,:,:,:,:) = dum_7D
            call dealloc_var_1D(var_1D)
            
            ! jac_FD
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'jac_FD',&
                &rich_lvl=id,lim_loc=lim_mem(1:6,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_6D)
            eq%jac_FD(par_id(1):par_id(2):par_id(3),:,:,:,:,:) = dum_6D
            call dealloc_var_1D(var_1D)
            
            ! S
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'S',&
                &rich_lvl=id,lim_loc=lim_mem(1:3,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_3D)
            eq%S(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
            call dealloc_var_1D(var_1D)
            
            ! kappa_n
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'kappa_n',rich_lvl=id,lim_loc=lim_mem(1:3,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_3D)
            eq%kappa_n(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
            call dealloc_var_1D(var_1D)
            
            ! kappa_g
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'kappa_g',rich_lvl=id,lim_loc=lim_mem(1:3,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_3D)
            eq%kappa_g(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
            call dealloc_var_1D(var_1D)
            
            ! sigma
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'sigma',&
                &rich_lvl=id,lim_loc=lim_mem(1:3,:))
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_3D)
            eq%sigma(par_id(1):par_id(2):par_id(3),:,:) = dum_3D
            call dealloc_var_1D(var_1D)
        end do
        
        ! clean up
        deallocate(dum_3D)
        deallocate(dum_6D)
        deallocate(dum_7D)
    end function reconstruct_PB3D_eq_2
    
    ! Reconstructs the vectorial perturbation variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name if it is > 0.
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    ! Furthermore,  using "lim_pos",  you can  obtain a  subset of  the data  by
    ! directly passing its limits to the underlying HDF5 routine. This refers to
    ! the position dimensions only. If provided,  the normal limits of a divided
    ! grid refer to the subset, as in "copy_grid".
    integer function reconstruct_PB3D_X_1(grid_X,X,data_name,rich_lvl,&
        &tot_rich,lim_sec_X,lim_pos) result(ierr)
        use num_vars, only: PB3D_name
        use X_vars, only: n_mod_X
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid 
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in), optional :: lim_pos(3,2)                           ! position limits of subset of data
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        integer :: lim_sec_X_loc(2)                                             ! local version of lim_sec_X
        integer :: lim_mem(4,2)                                                 ! memory limits for variables
        integer :: id                                                           ! counter
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: par_id_mem(2)                                                ! parallel indices (start, end) in memory (stride 1)
        integer :: par_lim(2)                                                   ! parallel limits
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        
        ! setup rich_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        
        ! set up local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set up parallel limits
        par_lim = [1,grid_X%n(1)]
        if (present(lim_pos)) then
            where (lim_pos(1,:).ge.0) par_lim = lim_pos(1,:)
        end if
        
        ! create X
        call X%init(grid_X,lim_sec_X)
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id
            par_id = setup_par_id(grid_X,rich_lvl_loc,id,tot_rich=tot_rich,&
                &par_lim=par_lim,par_id_mem=par_id_mem)
            
            ! set up local limits for HDF5 reconstruction
            lim_mem(1,:) = par_id_mem
            lim_mem(2,:) = [1,grid_X%n(2)]
            lim_mem(3,:) = [grid_X%i_min,grid_X%i_max]
            lim_mem(4,:) = lim_sec_X_loc
            if (present(lim_pos)) then
                lim_mem(2,:) = lim_pos(2,:)
                lim_mem(3,:) = lim_mem(3,:) + lim_pos(3,1) - 1                  ! take into account the grid limits (relative to position subset)
            end if
            
            ! RE_U_0
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'RE_U_0',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%U_0(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            call dealloc_var_1D(var_1D)
            
            ! IM_U_0
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'IM_U_0',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%U_0(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%U_0(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            call dealloc_var_1D(var_1D)
            
            ! RE_U_1
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'RE_U_1',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%U_1(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            call dealloc_var_1D(var_1D)
            
            ! IM_U_1
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'IM_U_1',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%U_1(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%U_1(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            call dealloc_var_1D(var_1D)
            
            ! RE_DU_0
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'RE_DU_0',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%DU_0(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            call dealloc_var_1D(var_1D)
            
            ! IM_DU_0
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'IM_DU_0',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%DU_0(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%DU_0(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            call dealloc_var_1D(var_1D)
            
            ! RE_DU_1
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'RE_DU_1',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%DU_1(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            call dealloc_var_1D(var_1D)
            
            ! IM_DU_1
            ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                &'IM_DU_1',rich_lvl=id,lim_loc=lim_mem)
            CHCKERR('')
            call conv_1D2ND(var_1D,dum_4D)
            X%DU_1(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%DU_1(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            call dealloc_var_1D(var_1D)
        end do
        
        ! clean up
        deallocate(dum_4D)
    end function reconstruct_PB3D_X_1
    
    ! Reconstructs the tensorial perturbation variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name if it is > 0.
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    ! Furthermore,  using "lim_pos",  you can  obtain a  subset of  the data  by
    ! directly passing its limits to the underlying HDF5 routine. This refers to
    ! the position dimensions only. If provided,  the normal limits of a divided
    ! grid refer to the subset, as in "copy_grid".
    ! Note: the tensorial perturbation type can  also be used for field- aligned
    ! variables, in  which case the first  index is assumed to  have dimension 1
    ! only. This can be triggered using "is_field_averaged".
    integer function reconstruct_PB3D_X_2(grid_X,X,data_name,rich_lvl,&
        &tot_rich,lim_sec_X,lim_pos,is_field_averaged) result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        use X_utilities, only: get_sec_X_range
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid 
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation vars
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        integer, intent(in), optional :: lim_pos(3,2)                           ! position limits of subset of data
        logical, intent(in), optional :: is_field_averaged                      ! if field-averaged, only one dimension for first index
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D var
        logical :: is_field_averaged_loc                                        ! local is_field_averaged
        integer :: id                                                           ! counter
        integer :: m, k                                                         ! counters
        integer :: sXr_loc(2,2)                                                 ! local secondary X limits for symmetric and asymmetric vars
        integer :: sXr_tot(2,2)                                                 ! total secondary X limits for symmetric and asymmetric vars
        integer :: nn_mod_loc(2)                                                ! local nr. of modes for symmetric and asymmetric vars
        integer :: lim_mem(4,2,2)                                               ! memory limits for variables for symmetric and asymmetric vars
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: par_id_mem(2)                                                ! parallel indices (start, end) in memory (stride 1)
        integer :: par_lim(2)                                                   ! parallel limits
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        logical :: read_this(2)                                                 ! whether symmetric and asymmetric variables need to be read
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        
        ! set up local is_field_averaged
        is_field_averaged_loc = .false.
        if (present(is_field_averaged)) is_field_averaged_loc = &
            &is_field_averaged
        
        ! setup rich_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        
        ! set up parallel limits
        par_lim = [1,grid_X%n(1)]
        if (present(lim_pos)) then
            where (lim_pos(1,:).ge.0) par_lim = lim_pos(1,:)
        end if
        
        ! create X
        call X%init(grid_X,lim_sec_X,is_field_averaged)
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id
            if (is_field_averaged_loc) then
                par_id = [1,1,1]                                                ! only first element
                par_id_mem = [1,1]
            else
                par_id = setup_par_id(grid_X,rich_lvl_loc,id,tot_rich=tot_rich,&
                    &par_lim=par_lim,par_id_mem=par_id_mem)
            end if
        
            ! loop over second dimension (horizontal)
            do m = 1,X%n_mod(2)
                ! get contiguous range of modes of this m
                call get_sec_X_range(sXr_loc(:,1),sXr_tot(:,1),m,.true.,&
                    &lim_sec_X)
                call get_sec_X_range(sXr_loc(:,2),sXr_tot(:,2),m,.false.,&
                    &lim_sec_X)
                nn_mod_loc = sXr_loc(2,:)-sXr_loc(1,:)+1
                read_this = .false.
                do k = 1,2
                    if (sXr_loc(1,k).le.sXr_loc(2,k)) read_this(k) = .true.     ! a bound is found
                end do
                
                ! set up local limits for HDF5 reconstruction of this m
                ! Note:  It  are  the  indices  in  total  matrix  sXr_tot  that
                ! correspond to the local limits.  "tot" just refers the to fact
                ! that they are valid for a  submatrix of the total matrix; They
                ! have been set up using the local grid_X limits as well.
                lim_mem(1,:,1) = par_id_mem
                lim_mem(2,:,1) = [1,grid_X%n(2)]
                lim_mem(3,:,1) = [grid_X%i_min,grid_X%i_max]
                lim_mem(4,:,1) = [sXr_tot(1,1),sXr_tot(2,1)]
                lim_mem(1,:,2) = par_id_mem
                lim_mem(2,:,2) = [1,grid_X%n(2)]
                lim_mem(3,:,2) = [grid_X%i_min,grid_X%i_max]
                lim_mem(4,:,2) = [sXr_tot(1,2),sXr_tot(2,2)]
                if (present(lim_pos)) then
                    lim_mem(2,:,1) = lim_pos(2,:)
                    lim_mem(3,:,1) = lim_mem(3,:,1) + lim_pos(3,1) - 1          ! take into account the grid limits (relative to position subset)
                    lim_mem(2,:,2) = lim_pos(2,:)
                    lim_mem(3,:,2) = lim_mem(3,:,2) + lim_pos(3,1) - 1          ! take into account the grid limits (relative to position subset)
                end if
                
                if (read_this(1)) then
                    ! RE_PV_0
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'RE_PV_0',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%PV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! IM_PV_0
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'IM_PV_0',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%PV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%PV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! RE_PV_2
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'RE_PV_2',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%PV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! IM_PV_2
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'IM_PV_2',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%PV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%PV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! RE_KV_0
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'RE_KV_0',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%KV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! IM_KV_0
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'IM_KV_0',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%KV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%KV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! RE_KV_2
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'RE_KV_2',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%KV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! IM_KV_2
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'IM_KV_2',rich_lvl=id,lim_loc=lim_mem(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%KV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%KV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    call dealloc_var_1D(var_1D)
                end if
                    
                if (read_this(2)) then
                    ! RE_PV_1
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'RE_PV_1',rich_lvl=id,lim_loc=lim_mem(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%PV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! IM_PV_1
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'IM_PV_1',rich_lvl=id,lim_loc=lim_mem(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%PV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = &
                        &X%PV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) + iu*dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! RE_KV_1
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'RE_KV_1',rich_lvl=id,lim_loc=lim_mem(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%KV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = dum_4D
                    call dealloc_var_1D(var_1D)
                    
                    ! IM_KV_1
                    ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                        &'IM_KV_1',rich_lvl=id,lim_loc=lim_mem(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(var_1D,dum_4D)
                    X%KV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = &
                        &X%KV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) + iu*dum_4D
                    call dealloc_var_1D(var_1D)
                end if
            end do
        end do
        
        ! clean up
        deallocate(dum_4D)
    end function reconstruct_PB3D_X_2
    
    ! Reconstructs the solution variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Furthermore,  using "lim_pos",  you can  obtain a  subset of  the data  by
    ! directly passing its limits to the underlying HDF5 routine. This refers to
    ! the position dimensions only. If provided,  the normal limits of a divided
    ! grid refer to the subset, as in "copy_grid".
    integer function reconstruct_PB3D_sol(grid_sol,sol,data_name,rich_lvl,&
        &lim_sec_sol,lim_pos) result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use X_vars, only: n_mod_X
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid 
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: lim_sec_sol(2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer, intent(in), optional :: lim_pos(1,2)                           ! position limits of subset of data
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        integer :: lim_sec_sol_loc(2)                                           ! local version of lim_sec_sol
        integer :: lim_mem(3,2)                                                 ! memory limits for variables
        integer :: n_EV                                                         ! nr. of Eigenvalues
        
        ! initialize ierr
        ierr = 0
        
        ! prepare
        
        ! set n_EV
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RE_sol_vec',&
            &rich_lvl=rich_lvl)
        CHCKERR('')
        n_EV = var_1D%tot_i_max(3)-var_1D%tot_i_min(3)+1
        call dealloc_var_1D(var_1D)
        
        ! set up local lim_sec_sol
        lim_sec_sol_loc = [1,n_mod_X]
        if (present(lim_sec_sol)) lim_sec_sol_loc = lim_sec_sol
        
        ! create solution
        call sol%init(grid_sol,n_EV,lim_sec_sol)
        
        ! set up local limits for HDF5 reconstruction
        lim_mem(1,:) = lim_sec_sol_loc
        lim_mem(2,:) = [grid_sol%i_min,grid_sol%i_max]
        lim_mem(3,:) = [-1,-1]
        if (present(lim_pos)) then
            lim_mem(2,:) = lim_mem(2,:) + lim_pos(1,1) - 1                  ! take into account the grid limits (relative to position subset)
        end if
        
        ! restore variables
        
        ! RE_sol_val
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RE_sol_val',&
            &rich_lvl=rich_lvl)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        sol%val = dum_1D
        deallocate(dum_1D)
        call dealloc_var_1D(var_1D)
        
        ! IM_sol_val
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'IM_sol_val',&
            &rich_lvl=rich_lvl)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        sol%val = sol%val + iu*dum_1D
        deallocate(dum_1D)
        call dealloc_var_1D(var_1D)
        
        ! RE_sol_vec
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RE_sol_vec',&
            &rich_lvl=rich_lvl,lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_3D)
        sol%vec = dum_3D
        deallocate(dum_3D)
        call dealloc_var_1D(var_1D)
        
        ! IM_sol_vec
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'IM_sol_vec',&
            &rich_lvl=rich_lvl,lim_loc=lim_mem)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_3D)
        sol%vec = sol%vec + iu*dum_3D
        deallocate(dum_3D)
        call dealloc_var_1D(var_1D)
    end function reconstruct_PB3D_sol
    
    ! get grid size
    ! Note: "grid_" is added in front the grid_name.
    integer function get_PB3D_grid_size(n,grid_name,rich_lvl,tot_rich) &
        &result(ierr)
        use num_vars, only: PB3D_name
        use PB3D_utilities, only: conv_1D2ND
        use HDF5_ops, only: read_HDF5_arr
        
        character(*), parameter :: rout_name = 'get_PB3D_grid_size'
        
        ! input / output
        integer, intent(inout) :: n(3)                                          ! n of grid
        character(len=*), intent(in) :: grid_name                               ! name of grid
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        real(dp), allocatable :: dum_1D(:)                                      ! dummy variables
        
        ! initialize ierr
        ierr = 0
        
        ! set n from HDF5 for local Richardson level
        n = 0
        ierr = read_HDF5_arr(var_1D,PB3D_name,'grid_'//trim(grid_name),&
            &'n',rich_lvl=rich_lvl)
        CHCKERR('')
        ! loop over all equilibrium jobs
        call conv_1D2ND(var_1D,dum_1D)
        n = nint(dum_1D)
        call dealloc_var_1D(var_1D)
        
        ! possibly only half of the points were saved in local Richardson level
        if (present(tot_rich) .and. present(rich_lvl)) then
            if (tot_rich .and. rich_lvl.gt.1) n(1) = n(1)*2+1
        end if
    end function get_PB3D_grid_size
end module PB3D_ops
