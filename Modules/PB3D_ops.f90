!------------------------------------------------------------------------------!
!   Operations on PB3D output.                                                 !
!------------------------------------------------------------------------------!
module PB3D_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, max_name_ln, iu, min_PB3D_version
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type
    use HDF5_vars, only: dealloc_var_1D, var_1D_type
    use PB3D_utilities, only: setup_rich_id, setup_eq_id, setup_par_id

    implicit none
    private
    public reconstruct_PB3D_in, reconstruct_PB3D_grid, reconstruct_PB3D_eq_1, &
        &reconstruct_PB3D_eq_2, reconstruct_PB3D_X_1, reconstruct_PB3D_X_2, &
        &reconstruct_PB3D_sol
    
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
            &EV_BC, magn_int_style, K_style
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm, &
            &max_flux_E, max_flux_F
        use grid_vars, onLy: n_r_in, n_r_eq, n_r_sol
        use X_vars, only: min_r_sol, max_r_sol, min_sec_X, max_sec_X, prim_X, &
            &n_mod_X
        use sol_vars, only: alpha
        use HELENA_vars, only: chi_H, flux_p_H, R_H, Z_H, nchi, ias, qs_H, &
            &pres_H, RBphi_H
        use VMEC, only: is_freeb_V, mnmax_V, mpol_V, ntor_V, is_asym_V, gam_V, &
            &R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mnmax_V, mn_V, rot_t_V, &
            &pres_V, flux_t_V, Dflux_t_V, flux_p_V, Dflux_p_V, nfp_V
#if ldebug
        use HELENA_vars, only: h_H_11, h_H_12, h_H_33
        use VMEC, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_in'
            
        ! input / output
        character(len=*), intent(in) :: data_name                               ! name of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(var_1D_type) :: var_1D                                             ! 1D variable
        real(dp), parameter :: tol_version = 1.E-4_dp                           ! tolerance for version control
        real(dp) :: PB3D_version                                                ! version of PB3D variable read
        
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
        deallocate(dum_1D)
        call dealloc_var_1D(var_1D)
        
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
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! flux_t_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_t_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(flux_t_V(n_r_eq))
                flux_t_V = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! Dflux_t_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'Dflux_t_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(Dflux_t_V(n_r_eq))
                Dflux_t_V = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! flux_p_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_p_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(flux_p_V(n_r_eq))
                flux_p_V = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! Dflux_p_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'Dflux_p_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(Dflux_p_V(n_r_eq))
                Dflux_p_V = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! pres_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'pres_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(pres_V(n_r_eq))
                pres_V = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! rot_t_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rot_t_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(rot_t_V(n_r_eq))
                rot_t_V = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! mn_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'mn_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_2D)
                allocate(mn_V(mnmax_V,2))
                mn_V = nint(dum_2D)
                deallocate(dum_2D)
                call dealloc_var_1D(var_1D)
                
                ! RZL_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RZL_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(R_V_c(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(R_V_s(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(Z_V_c(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(Z_V_s(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(L_V_c(mnmax_V,n_r_eq,0:max_deriv+1))
                allocate(L_V_s(mnmax_V,n_r_eq,0:max_deriv+1))
                R_V_c(:,:,0) = dum_3D(:,:,1)
                R_V_s(:,:,0) = dum_3D(:,:,2)
                Z_V_c(:,:,0) = dum_3D(:,:,3)
                Z_V_s(:,:,0) = dum_3D(:,:,4)
                L_V_c(:,:,0) = dum_3D(:,:,5)
                L_V_s(:,:,0) = dum_3D(:,:,6)
                deallocate(dum_3D)
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
                deallocate(dum_4D)
                call dealloc_var_1D(var_1D)
                
                ! B_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'B_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(B_V_c(mnmax_V,n_r_eq))
                allocate(B_V_s(mnmax_V,n_r_eq))
                B_V_c = dum_3D(:,:,1)
                B_V_s = dum_3D(:,:,2)
                deallocate(dum_3D)
                call dealloc_var_1D(var_1D)
                
                ! jac_V
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'jac_V')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(jac_V_c(mnmax_V,n_r_eq))
                allocate(jac_V_s(mnmax_V,n_r_eq))
                jac_V_c = dum_3D(:,:,1)
                jac_V_s = dum_3D(:,:,2)
                deallocate(dum_3D)
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
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! RZ_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RZ_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_3D)
                allocate(R_H(nchi,n_r_eq))
                allocate(Z_H(nchi,n_r_eq))
                R_H = dum_3D(:,:,1)
                Z_H = dum_3D(:,:,2)
                deallocate(dum_3D)
                call dealloc_var_1D(var_1D)
                
                ! chi_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'chi_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(chi_H(nchi))
                chi_H = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! flux_p_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),&
                    &'flux_p_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(flux_p_H(n_r_eq))
                flux_p_H = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! qs_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'qs_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(qs_H(n_r_eq))
                qs_H = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! pres_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'pres_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(pres_H(n_r_eq))
                pres_H = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
                ! RBphi_H
                ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'RBphi_H')
                CHCKERR('')
                call conv_1D2ND(var_1D,dum_1D)
                allocate(RBphi_H(n_r_eq))
                RBphi_H = dum_1D
                deallocate(dum_1D)
                call dealloc_var_1D(var_1D)
                
#if ldebug
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
                deallocate(dum_3D)
                call dealloc_var_1D(var_1D)
#endif
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
        magn_int_style = nint(dum_1D(10))
        K_style = nint(dum_1D(11))
        deallocate(dum_1D)
        call dealloc_var_1D(var_1D)
        
        ! misc_sol
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'misc_sol')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        min_r_sol = dum_1D(1)
        max_r_sol = dum_1D(2)
        alpha = dum_1D(3)
        norm_disc_prec_sol = nint(dum_1D(4))
        BC_style = nint(dum_1D(5))
        EV_style = nint(dum_1D(6))
        EV_BC = dum_1D(7)
        deallocate(dum_1D)
        call dealloc_var_1D(var_1D)
        
        ! user output
        call lvl_ud(-1)
        call writo('Input variables from PB3D output reconstructed')
        
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
                
                call lvl_ud(-1)
        end select
    end function reconstruct_PB3D_in
    
    ! Reconstructs grid variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name  if it  is  >  0, and  similarly  for  "eq_job" through  "_E_eq_job".
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    ! Note: "grid_" is added in front the data_name.
    integer function reconstruct_PB3D_grid(grid,data_name,rich_lvl,eq_job,&
        &tot_rich,grid_limits) result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_grid'
        
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid 
        character(len=*), intent(in) :: data_name                               ! name of grid
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to print
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: grid_limits(2)                         ! i_limit of grid
        
        ! local variables
        type(var_1D_type), allocatable :: vars_1D(:)                            ! 1D variables
        integer :: grid_limits_loc(2)                                           ! local versions of grid_limits
        integer :: n(3)                                                         ! total n
        integer :: loc_n_r                                                      ! local n
        integer :: id                                                           ! counter
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        integer :: eq_id(2)                                                     ! equilibrium job indices (start, end)
        logical :: overlap                                                      ! overlap in equilibrium jobs if first Richardson level
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        
        ! setup rich_id and eq_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        eq_id = setup_eq_id('grid_'//trim(data_name),eq_job=eq_job,&
            &rich_lvl=rich_lvl)
        if (minval(eq_id).lt.0) then
            ierr = 1
            err_msg = 'Variable "grid_'//trim(data_name)//'" not found, nor &
                &using _R_ and _E_ suffixes'
            CHCKERR(err_msg)
        end if
        
        ! set n
        n = 0
        ierr = read_HDF5_arr(vars_1D,PB3D_name,'grid_'//trim(data_name),&
            &'n',rich_lvl=rich_lvl_loc,eq_job=eq_id)
        CHCKERR('')
        do id = 1,eq_id(2)-eq_id(1)+1
            call conv_1D2ND(vars_1D(id),dum_1D)
            if (id.eq.1) then
                n = nint(dum_1D)
            else if (n(2).eq.nint(dum_1D(2)) .and. n(3).eq.nint(dum_1D(3))) then
                n(1) = n(1) + nint(dum_1D(1))
                if (rich_lvl_loc.le.1) n(1) = n(1) - 1
            else
                ierr = 1
                err_msg = 'Variables for different equilibrium jobs not &
                    &compatible'
                CHCKERR(err_msg)
            end if
            deallocate(dum_1D)
        end do
        call dealloc_var_1D(vars_1D)
        
        ! possibly only half of the points were saved in last Richardson level
        if (present(tot_rich) .and. rich_lvl_loc.gt.1) then
            if (tot_rich) n(1) = n(1)*2+1
        end if
        
        ! set up local grid_limits
        grid_limits_loc = [1,n(3)]
        if (present(grid_limits)) grid_limits_loc = grid_limits
        
        ! create grid
        ierr = grid%init(n,grid_limits_loc)
        CHCKERR('')
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id, loc_n_r and overlap
            par_id = setup_par_id(grid,rich_lvl_loc,id,tot_rich)
            loc_n_r = (par_id(2)-par_id(1))/par_id(3)+1
            overlap = id.le.1                                                   ! overlap for first Richardson level
            
            ! r_F
            allocate(vars_1D(1))
            ierr = read_HDF5_arr(vars_1D(1),PB3D_name,'grid_'//&
                &trim(data_name),'r_F',rich_lvl=id,eq_job=eq_id(1))
            CHCKERR('')
            call conv_1D2ND(vars_1D(1),dum_1D)                                  ! use first equilibrium job, as it should be invariant
            grid%r_F = dum_1D
            deallocate(dum_1D)
            call dealloc_var_1D(vars_1D)
            
            ! r_E
            allocate(vars_1D(1))
            ierr = read_HDF5_arr(vars_1D(1),PB3D_name,'grid_'//&
                &trim(data_name),'r_E',rich_lvl=id,eq_job=eq_id(1))
            CHCKERR('')
            call conv_1D2ND(vars_1D(1),dum_1D)                                  ! use first equilibrium job, as it should be invariant
            grid%r_E = dum_1D
            deallocate(dum_1D)
            call dealloc_var_1D(vars_1D)
            
            ! loc_r_F
            grid%loc_r_F = grid%r_F(grid_limits_loc(1):grid_limits_loc(2))
            
            ! loc_r_E
            grid%loc_r_E = grid%r_E(grid_limits_loc(1):grid_limits_loc(2))
            
            ! only for 3D grids
            if (product(grid%n(1:2)).ne.0) then
                ! theta_F
                ierr = read_HDF5_arr(vars_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'theta_F',rich_lvl=id,eq_job=eq_id)
                CHCKERR('')
                call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
                grid%theta_F(par_id(1):par_id(2):par_id(3),:,:) = &
                    &dum_3D(:,:,grid_limits_loc(1):grid_limits_loc(2))
                deallocate(dum_3D)
                call dealloc_var_1D(vars_1D)
                
                ! theta_E
                ierr = read_HDF5_arr(vars_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'theta_E',rich_lvl=id,eq_job=eq_id)
                CHCKERR('')
                call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
                grid%theta_E(par_id(1):par_id(2):par_id(3),:,:) = &
                    &dum_3D(:,:,grid_limits_loc(1):grid_limits_loc(2))
                deallocate(dum_3D)
                call dealloc_var_1D(vars_1D)
                
                ! zeta_F
                ierr = read_HDF5_arr(vars_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'zeta_F',rich_lvl=id,eq_job=eq_id)
                CHCKERR('')
                call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
                grid%zeta_F(par_id(1):par_id(2):par_id(3),:,:) = &
                    &dum_3D(:,:,grid_limits_loc(1):grid_limits_loc(2))
                deallocate(dum_3D)
                call dealloc_var_1D(vars_1D)
                
                ! zeta_E
                ierr = read_HDF5_arr(vars_1D,PB3D_name,'grid_'//&
                    &trim(data_name),'zeta_E',rich_lvl=id,eq_job=eq_id)
                CHCKERR('')
                call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
                grid%zeta_E(par_id(1):par_id(2):par_id(3),:,:) = &
                    &dum_3D(:,:,grid_limits_loc(1):grid_limits_loc(2))
                deallocate(dum_3D)
                call dealloc_var_1D(vars_1D)
            end if
        end do
    end function reconstruct_PB3D_grid
    
    ! Reconstructs the equilibrium variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    integer function reconstruct_PB3D_eq_1(grid_eq,eq,data_name) &
        &result(ierr)                                                           ! flux version
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_eq_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid 
        type(eq_1_type), intent(inout), optional :: eq                          ! flux equilibrium
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
        
        ! initialize ierr
        ierr = 0
        
        ! prepare
        
        ! create equilibrium
        call eq%init(grid_eq,setup_E=.false.,setup_F=.true.)
        
        ! restore variables
        
        ! pres_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'pres_FD')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%pres_FD = dum_2D(grid_eq%i_min:grid_eq%i_max,:)
        deallocate(dum_2D)
        call dealloc_var_1D(var_1D)
        
        ! q_saf_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'q_saf_FD')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%q_saf_FD = dum_2D(grid_eq%i_min:grid_eq%i_max,:)
        deallocate(dum_2D)
        call dealloc_var_1D(var_1D)
        
        ! rot_t_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rot_t_FD')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%rot_t_FD = dum_2D(grid_eq%i_min:grid_eq%i_max,:)
        deallocate(dum_2D)
        call dealloc_var_1D(var_1D)
        
        ! flux_p_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'flux_p_FD')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%flux_p_FD = dum_2D(grid_eq%i_min:grid_eq%i_max,:)
        deallocate(dum_2D)
        call dealloc_var_1D(var_1D)
        
        ! flux_t_FD
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'flux_t_FD')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_2D)
        eq%flux_t_FD = dum_2D(grid_eq%i_min:grid_eq%i_max,:)
        deallocate(dum_2D)
        call dealloc_var_1D(var_1D)
        
        ! rho
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'rho')
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_1D)
        eq%rho = dum_1D(grid_eq%i_min:grid_eq%i_max)
        deallocate(dum_1D)
        call dealloc_var_1D(var_1D)
    end function reconstruct_PB3D_eq_1
    
    ! Reconstructs the equilibrium variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name  if it  is  >  0, and  similarly  for  "eq_job" through  "_E_eq_job".
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    integer function reconstruct_PB3D_eq_2(grid_eq,eq,data_name,rich_lvl,&
        &eq_job,tot_rich) result(ierr)                                          ! metric version
        use num_vars, only: PB3D_name_eq
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_eq_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid 
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to print
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        
        ! local variables
        type(var_1D_type), allocatable :: vars_1D(:)                            ! 1D variables
        integer :: loc_n_r                                                      ! local n
        integer :: id                                                           ! counter
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        integer :: eq_id(2)                                                     ! equilibrium job indices (start, end)
        logical :: overlap                                                      ! overlap in equilibrium jobs if first Richardson level
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl and overlap
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        overlap = rich_lvl_loc.le.1
        
        ! setup rich_id and eq_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        eq_id = setup_eq_id(trim(data_name),eq_job=eq_job,rich_lvl=rich_lvl)
        if (minval(eq_id).lt.0) then
            ierr = 1
            err_msg = 'Variable "'//trim(data_name)//'" not found, nor &
                &using _R_ and _E_ suffixes'
            CHCKERR(err_msg)
        end if
        
        ! create equilibrium
        call eq%init(grid_eq,setup_E=.false.,setup_F=.true.)
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id
            par_id = setup_par_id(grid_eq,rich_lvl_loc,id,tot_rich)
            loc_n_r = (par_id(2)-par_id(1))/par_id(3)+1
            
            ! read HDF5 variables
            !ierr = read_HDF5_arrs(vars_1D,PB3D_name_eq,trim(data_name),&
                !&rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            
            ! g_FD
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),'g_FD',&
                &rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_7D)
            eq%g_FD(par_id(1):par_id(2):par_id(3),:,:,:,:,:,:) = &
                &dum_7D(:,:,grid_eq%i_min:grid_eq%i_max,:,:,:,:)
            deallocate(dum_7D)
            call dealloc_var_1D(vars_1D)
            
            ! h_FD
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),'h_FD',&
                &rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_7D)
            eq%h_FD(par_id(1):par_id(2):par_id(3),:,:,:,:,:,:) = &
                &dum_7D(:,:,grid_eq%i_min:grid_eq%i_max,:,:,:,:)
            deallocate(dum_7D)
            call dealloc_var_1D(vars_1D)
            
            ! jac_FD
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),'jac_FD',&
                &rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_6D)
            eq%jac_FD(par_id(1):par_id(2):par_id(3),:,:,:,:,:) = &
                &dum_6D(:,:,grid_eq%i_min:grid_eq%i_max,:,:,:)
            deallocate(dum_6D)
            call dealloc_var_1D(vars_1D)
            
            ! S
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),'S',&
                &rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
            eq%S(par_id(1):par_id(2):par_id(3),:,:) = &
                &dum_3D(:,:,grid_eq%i_min:grid_eq%i_max)
            deallocate(dum_3D)
            call dealloc_var_1D(vars_1D)
            
            ! kappa_n
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'kappa_n',rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
            eq%kappa_n(par_id(1):par_id(2):par_id(3),:,:) = &
                &dum_3D(:,:,grid_eq%i_min:grid_eq%i_max)
            deallocate(dum_3D)
            call dealloc_var_1D(vars_1D)
            
            ! kappa_g
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'kappa_g',rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
            eq%kappa_g(par_id(1):par_id(2):par_id(3),:,:) = &
                &dum_3D(:,:,grid_eq%i_min:grid_eq%i_max)
            deallocate(dum_3D)
            call dealloc_var_1D(vars_1D)
            
            ! sigma
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),'sigma',&
                &rich_lvl=id,eq_job=eq_id)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_3D)
            eq%sigma(par_id(1):par_id(2):par_id(3),:,:) = &
                &dum_3D(:,:,grid_eq%i_min:grid_eq%i_max)
            deallocate(dum_3D)
            call dealloc_var_1D(vars_1D)
        end do
    end function reconstruct_PB3D_eq_2
    
    ! Reconstructs the vectorial perturbation variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name  if it  is  >  0, and  similarly  for  "eq_job" through  "_E_eq_job".
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    integer function reconstruct_PB3D_X_1(grid_X,X,data_name,rich_lvl,eq_job,&
        &tot_rich,lim_sec_X) result(ierr)
        use num_vars, only: PB3D_name_eq
        use X_vars, only: n_mod_X
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid 
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to print
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        type(var_1D_type), allocatable :: vars_1D(:)                            ! 1D variables
        integer :: lim_sec_X_loc(2)                                             ! local version of lim_sec_X
        integer :: lim_loc(4,2)                                                 ! local limits for variables
        integer :: loc_n_r                                                      ! local n
        integer :: id                                                           ! counter
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        integer :: eq_id(2)                                                     ! equilibrium job indices (start, end)
        logical :: overlap                                                      ! overlap in equilibrium jobs if first Richardson level
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl and overlap
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        overlap = rich_lvl_loc.le.1
        
        ! setup rich_id and eq_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        eq_id = setup_eq_id(trim(data_name),eq_job=eq_job,rich_lvl=rich_lvl)
        if (minval(eq_id).lt.0) then
            ierr = 1
            err_msg = 'Variable "'//trim(data_name)//'" not found, nor &
                &using _R_ and _E_ suffixes'
            CHCKERR(err_msg)
        end if
        
        ! set up local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set up local limits for HDF5 reconstruction
        lim_loc(:,1) = [1,1,grid_X%i_min,lim_sec_X_loc(1)]
        lim_loc(:,2) = [grid_X%n(1:2),grid_X%i_max,lim_sec_X_loc(2)]
        
        ! create X
        call X%init(grid_X,lim_sec_X)
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id
            par_id = setup_par_id(grid_X,rich_lvl_loc,id,tot_rich)
            loc_n_r = (par_id(2)-par_id(1))/par_id(3)+1
            
            ! RE_U_0
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'RE_U_0',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%U_0(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! IM_U_0
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'IM_U_0',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%U_0(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%U_0(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! RE_U_1
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'RE_U_1',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%U_1(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! IM_U_1
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'IM_U_1',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%U_1(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%U_1(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! RE_DU_0
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'RE_DU_0',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%DU_0(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! IM_DU_0
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'IM_DU_0',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%DU_0(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%DU_0(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! RE_DU_1
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'RE_DU_1',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%DU_1(par_id(1):par_id(2):par_id(3),:,:,:) = dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
            
            ! IM_DU_1
            ierr = read_HDF5_arr(vars_1D,PB3D_name_eq,trim(data_name),&
                &'IM_DU_1',rich_lvl=id,eq_job=eq_id,lim_loc=lim_loc)
            CHCKERR('')
            call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
            X%DU_1(par_id(1):par_id(2):par_id(3),:,:,:) = &
                X%DU_1(par_id(1):par_id(2):par_id(3),:,:,:) + iu*dum_4D
            deallocate(dum_4D)
            call dealloc_var_1D(vars_1D)
        end do
    end function reconstruct_PB3D_X_1
    
    ! Reconstructs the tensorial perturbation variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the data
    ! name  if it  is  >  0, and  similarly  for  "eq_job" through  "_E_eq_job".
    ! With "tot_rich"  the information  from previous  Richardson levels  can be
    ! combined.
    ! Note: the tensorial perturbation type can  also be used for field- aligned
    ! variables, in  which case the first  index is assumed to  have dimension 1
    ! only. This can be triggered using "is_field_averaged".
    integer function reconstruct_PB3D_X_2(grid_X,X,data_name,rich_lvl,eq_job,&
        &tot_rich,lim_sec_X,is_field_averaged) result(ierr)
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
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to print
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        logical, intent(in), optional :: is_field_averaged                      ! if field-averaged, only one dimension for first index
        
        ! local variables
        type(var_1D_type), allocatable :: vars_1D(:)                            ! 1D vars
        logical :: is_field_averaged_loc                                        ! local is_field_averaged
        integer :: id                                                           ! counter
        integer :: m, k                                                         ! counters
        integer :: sXr_loc(2,2)                                                 ! local secondary X limits for symmetric and asymmetric vars
        integer :: sXr_tot(2,2)                                                 ! total secondary X limits for symmetric and asymmetric vars
        integer :: lim_loc(4,2,2)                                               ! local limits for vars for symmetric and asymmetric vars
        integer :: nn_mod_loc(2)                                                ! local nr. of modes for symmetric and asymmetric vars
        integer :: loc_n_r                                                      ! local n
        integer :: rich_lvl_loc                                                 ! local rich_lvl
        integer :: par_id(3)                                                    ! parallel indices (start, end, stride)
        integer :: par_id_loc(2)                                                ! local parallel indices (start, end)
        integer :: rich_id(2)                                                   ! richardson level indices (start, end)
        integer :: eq_id(2)                                                     ! equilibrium job indices (start, end)
        logical :: read_this(2)                                                 ! whether symmetric and asymmetric variables need to be read
        logical :: overlap                                                      ! overlap in equilibrium jobs if first Richardson level
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up local rich_lvl and overlap
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        overlap = rich_lvl_loc.le.1
        
        ! set up local is_field_averaged
        is_field_averaged_loc = .false.
        if (present(is_field_averaged)) is_field_averaged_loc = &
            &is_field_averaged
        
        ! setup rich_id and eq_id
        rich_id = setup_rich_id(rich_lvl_loc,tot_rich)
        eq_id = setup_eq_id(trim(data_name),eq_job=eq_job,rich_lvl=rich_lvl)
        if (minval(eq_id).lt.0) then
            ierr = 1
            err_msg = 'Variable "'//trim(data_name)//'" not found, nor &
                &using _R_ and _E_ suffixes'
            CHCKERR(err_msg)
        end if
        
        ! create X
        call X%init(grid_X,lim_sec_X,is_field_averaged)
        
        ! restore looping over richardson levels
        do id = rich_id(2),rich_id(1),-1
            ! setup par_id
            if (is_field_averaged_loc) then
                par_id = [1,1,1]                                                ! only first element
                par_id_loc = [1,1]                                              ! only first element
            else
                par_id = setup_par_id(grid_X,rich_lvl_loc,id,tot_rich)
                par_id_loc = [1,(par_id(2)-par_id(1))/par_id(3)+1]              ! the number of elements par_id(1):par_id(2), stride par_id(3)
            end if
            loc_n_r = (par_id(2)-par_id(1))/par_id(3)+1
        
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
                ! correspond to the local limits. "tot" just refers to fact that
                ! they are valid for a submatrix  of the total matrix; They have
                ! been set up using the local grid_X limits as well.
                lim_loc(:,1,1) = [par_id_loc(1),1,grid_X%i_min,sXr_tot(1,1)]    ! lower limits for symmetric vars.
                lim_loc(:,2,1) = [par_id_loc(2),grid_X%n(2),grid_X%i_max,&
                    &sXr_tot(2,1)]                                              ! upper limits for symmetric vars.
                lim_loc(:,1,2) = [par_id_loc(1),1,grid_X%i_min,sXr_tot(1,2)]    ! lower limits for asymmetric vars.
                lim_loc(:,2,2) = [par_id_loc(2),grid_X%n(2),grid_X%i_max,&
                    &sXr_tot(2,2)]                                              ! upper limits for asymmetric vars.
                
                if (read_this(1)) then
                    ! RE_PV_0
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'RE_PV_0',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%PV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! IM_PV_0
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'IM_PV_0',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%PV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%PV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! RE_PV_2
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'RE_PV_2',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%PV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! IM_PV_2
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'IM_PV_2',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%PV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%PV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! RE_KV_0
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'RE_KV_0',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%KV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! IM_KV_0
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'IM_KV_0',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%KV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%KV_0(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! RE_KV_2
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'RE_KV_2',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%KV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! IM_KV_2
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'IM_KV_2',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,1))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%KV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) = &
                        &X%KV_2(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,1):sXr_loc(2,1)) + iu*dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                end if
                    
                if (read_this(2)) then
                    ! RE_PV_1
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'RE_PV_1',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%PV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! IM_PV_1
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'IM_PV_1',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%PV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = &
                        &X%PV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) + iu*dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! RE_KV_1
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'RE_KV_1',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%KV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                    
                    ! IM_KV_1
                    ierr = read_HDF5_arr(vars_1D,PB3D_name,trim(data_name),&
                        &'IM_KV_1',rich_lvl=id,eq_job=eq_id,&
                        &lim_loc=lim_loc(:,:,2))
                    CHCKERR('')
                    call conv_1D2ND(vars_1D,loc_n_r,overlap,dum_4D)
                    X%KV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) = &
                        &X%KV_1(par_id(1):par_id(2):par_id(3),:,:,&
                        &sXr_loc(1,2):sXr_loc(2,2)) + iu*dum_4D
                    deallocate(dum_4D)
                    call dealloc_var_1D(vars_1D)
                end if
            end do
        end do
    end function reconstruct_PB3D_X_2
    
    ! Reconstructs the solution variables from PB3D output.
    ! Optionally, the grid limits can be provided.
    integer function reconstruct_PB3D_sol(grid_sol,sol,data_name,rich_lvl,&
        &lim_sec_sol) result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: read_HDF5_arr
        use PB3D_utilities, only: conv_1D2ND
        
        character(*), parameter :: rout_name = 'reconstruct_PB3D_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid 
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        character(len=*), intent(in) :: data_name                               ! name to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: lim_sec_sol(2)                         ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        type(var_1D_type) :: var_1D                                             ! 1D variable
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
        
        ! create solution
        call sol%init(grid_sol,n_EV,lim_sec_sol)
        
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
            &rich_lvl=rich_lvl)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_3D)
        sol%vec = dum_3D(:,grid_sol%i_min:grid_sol%i_max,:)
        deallocate(dum_3D)
        call dealloc_var_1D(var_1D)
        
        ! IM_sol_vec
        ierr = read_HDF5_arr(var_1D,PB3D_name,trim(data_name),'IM_sol_vec',&
            &rich_lvl=rich_lvl)
        CHCKERR('')
        call conv_1D2ND(var_1D,dum_3D)
        sol%vec = sol%vec + iu*dum_3D(:,grid_sol%i_min:grid_sol%i_max,:)
        deallocate(dum_3D)
        call dealloc_var_1D(var_1D)
    end function reconstruct_PB3D_sol
end module PB3D_ops
