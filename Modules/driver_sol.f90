!------------------------------------------------------------------------------!
!> Driver of the solution part of PB3D.
!------------------------------------------------------------------------------!
module driver_sol
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type, disc_type
    use X_vars, only: X_2_type, &
        &mds_X, mds_sol
    use vac_vars, only: vac_type
    use sol_vars, only: sol_type
    
    implicit none
    private
    public run_driver_sol
    
contains
    !> Main driver of PB3D solution part.
    !!
    !!  - sets up:
    !!      * \c grid_sol (only first Richardson level)
    !!      * \c sol
    !!  - writes to HDF5:
    !!      * \c grid_sol (only first Richardson level)
    !!      * \c sol
    !!  - deallocates:
    !!      * sol before setting up (but after guess)
    !!
    !! \return ierr
    integer function run_driver_sol(grid_eq,grid_X,grid_sol,X,vac,sol) &
        &result(ierr)
        
        use num_vars, only: EV_style, eq_style, rich_restart_lvl, rank, &
            &n_procs, X_grid_style, V_interp_style, norm_disc_prec_X, &
            &use_pol_flux_F
        use grid_vars, only: n_r_sol, n_alpha
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_sol
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_grid_sol, &
            &print_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_lvl
        use rich_ops, only: calc_rich_ex
        use vac_ops, only: calc_vac_res, print_output_vac
        use MPI_utilities, only: get_ser_var
        use grid_utilities, only: trim_grid, setup_interp_data, apply_disc
        use X_ops, only: setup_modes
        use X_vars, only: min_m_X, min_n_X
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in), target :: grid_X                           !< perturbation grid
        type(grid_type), intent(inout) :: grid_sol                              !< solution grid
        type(X_2_type), intent(in) :: X                                         !< integrated tensorial perturbation variables
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        type(sol_type), intent(inout) :: sol                                    !< solution variables
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed solution grid
        type(grid_type) :: grid_X_trim                                          ! trimmed perturbation grid
        type(grid_type) :: grid_X_ser                                           ! serial perturbation grid
        type(grid_type) :: grid_X_sol                                           ! perturbation grid with solution normal part
        type(X_2_type) :: X_ser                                                 ! serial X
        type(X_2_type) :: X_sol                                                 ! interpolated X
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation X->sol
        integer :: pmone                                                        ! plus (increasing r_F) or minus (decreasing r_F) one
        integer :: min_nm_X                                                     ! minimal n (tor. flux) or m (pol. flux)
        integer :: kd, ld                                                       ! counters
        integer :: n_mod_tot                                                    ! total number of modes
        integer :: norm_range_X(2)                                              ! normal range in X grid of total mode
        integer :: norm_range_sol(2)                                            ! normal range in sol grid of total mode
        integer :: n_X(3)                                                       ! n of grid_X_trim
        integer :: n_X_loc                                                      ! local size in grid_X_trim
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed solution grid
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        real(dp), allocatable :: r_F_X(:)                                       ! normal points in X grid
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
        complex(dp), allocatable :: V_X(:,:,:)                                  ! PV or KV for total mode
        complex(dp), allocatable :: ser_var_loc(:)                              ! local serial variable
        logical :: do_vac_ops                                                   ! whether specific calculations for vacuum are necessary
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up solution grid if first level
        if (rich_lvl.eq.rich_restart_lvl) then
            ! Divide solution grid under group processes, calculating the limits
            ! and the normal coordinate.
            allocate(r_F_sol(n_r_sol))
            ierr = calc_norm_range(sol_limits=sol_limits,r_F_sol=r_F_sol)
            CHCKERR('')
            
            if (rich_lvl.eq.1) then
                call writo('Set up solution grid')
                call lvl_ud(1)
                
                call writo('Calculate the grid')
                call lvl_ud(1)
                ierr = setup_grid_sol(grid_X,grid_sol,r_F_sol,sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                call writo('Write to output file')
                call lvl_ud(1)
                ierr = print_output_grid(grid_sol,'solution','sol')
                CHCKERR('')
                call lvl_ud(-1)
                
                ! set up modes
                ierr = setup_modes(mds_sol,grid_eq,grid_sol,plot_nm=.false.)
                CHCKERR('')
                
                call lvl_ud(-1)
            else
                ! restore solution grid and trim it
                ierr = reconstruct_PB3D_grid(grid_sol,'sol',&
                    &grid_limits=sol_limits)
                CHCKERR('')
                ierr = trim_grid(grid_sol,grid_sol_trim)
                CHCKERR('')
                
                ! reconstruct solution on trimmed grid
                ierr = reconstruct_PB3D_sol(mds_sol,grid_sol_trim,sol,'sol',&
                    &rich_lvl=rich_lvl-1)
                CHCKERR('')
                
                ! clean up
                call grid_sol_trim%dealloc()
            end if
            
            deallocate(r_F_sol)
        end if
        
        ! set up  whether Richardson level  has to be  appended to the  name and
        ! whether to do vacuum operations
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
                do_vac_ops = .true.
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
                if (rich_lvl.eq.rich_restart_lvl) then
                    do_vac_ops = .true.
                else
                    do_vac_ops = .false.
                end if
        end select
        write(*,*) '¡¡¡¡¡ NO VACUUM !!!!!'
        do_vac_ops = .false.
        
        if (do_vac_ops) then
            ! calculate vacuum
            ierr = calc_vac_res(mds_sol,vac)
            CHCKERR('')
            
            call writo('Write to output file')
            call lvl_ud(1)
            CHCKERR('')
            if (rank.eq.n_procs-1) then
                ierr = print_output_vac(vac,'vac',rich_lvl=rich_lvl_name)
                CHCKERR('')
            end if
            call lvl_ud(-1)
        end if
        
        ! initialize interpolated X
        ierr = grid_X_sol%init([1,n_alpha,grid_sol%n(3)],&
            &[grid_sol%i_min,grid_sol%i_max],grid_sol%divided)
        call X_sol%init(mds_sol,grid_X_sol,is_field_averaged=.true.)
        call grid_X_sol%dealloc()
        
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! user output
                call writo('Interpolate the perturbation variables to &
                    &solution grid')
                call lvl_ud(1)
                
                ! interpolate depending on style
                select case (V_interp_style)
                    case (1)                                                ! finite differences
                        ! user output
                        call writo('Using finite differences')
                        
                        ! set up normal interpolation factors
                        ierr = setup_interp_data(r_F_X,r_F_sol,&
                            &norm_interp_data,norm_disc_prec_X)
                        
                        ! interpolate
                        ierr = apply_disc(X%PV_0,norm_interp_data,&
                            &X_sol%PV_0,3)
                        CHCKERR('')
                        ierr = apply_disc(X%PV_1,norm_interp_data,&
                            &X_sol%PV_1,3)
                        CHCKERR('')
                        ierr = apply_disc(X%PV_2,norm_interp_data,&
                            &X_sol%PV_2,3)
                        CHCKERR('')
                        ierr = apply_disc(X%KV_0,norm_interp_data,&
                            &X_sol%KV_0,3)
                        CHCKERR('')
                        ierr = apply_disc(X%KV_1,norm_interp_data,&
                            &X_sol%KV_1,3)
                        CHCKERR('')
                        ierr = apply_disc(X%KV_2,norm_interp_data,&
                            &X_sol%KV_2,3)
                        CHCKERR('')
                        
                        ! clean up
                        call norm_interp_data%dealloc()
                    case (2)                                                ! splines
                        ! user output
                        call writo('Using splines')
                        
                        ! set up trimmed grid
                        ierr = trim_grid(grid_X,grid_X_trim,norm_id)
                        CHCKERR('')
                        
                        ! set up undivided grid
                        n_X = grid_X_trim%n
                        n_X(1) = 1                                          ! field-aligned
                        n_X_loc = grid_X_trim%n(2)*grid_X_trim%loc_n_r
                        ierr = grid_X_ser%init(n_X)
                        CHCKERR('')
                        grid_X_ser%r_F = grid_X_trim%r_F
                        call X_ser%init(mds_X,grid_X_ser,&
                            &is_field_averaged=.true.)
                        
                        allocate(ser_var_loc(product(n_X)))
                        if (grid_X_trim%divided) then
                            do ld = 1,size(X%PV_0,4)
                                ! PV_0
                                ierr = get_ser_var(reshape(&
                                    &X%PV_0(:,:,norm_id(1):norm_id(2),ld),&
                                    &[n_X_loc]),ser_var_loc,scatter=.true.)
                                CHCKERR('')
                                X_ser%PV_0(:,:,:,ld) = &
                                    &reshape(ser_var_loc,n_X)
                                
                                ! PV_2
                                ierr = get_ser_var(reshape(&
                                    &X%PV_2(:,:,norm_id(1):norm_id(2),ld),&
                                    &[n_X_loc]),ser_var_loc,scatter=.true.)
                                CHCKERR('')
                                X_ser%PV_2(:,:,:,ld) = &
                                    &reshape(ser_var_loc,n_X)
                                
                                ! KV_0
                                ierr = get_ser_var(reshape(&
                                    &X%KV_0(:,:,norm_id(1):norm_id(2),ld),&
                                    &[n_X_loc]),ser_var_loc,scatter=.true.)
                                CHCKERR('')
                                X_ser%KV_0(:,:,:,ld) = &
                                    &reshape(ser_var_loc,n_X)
                                
                                ! KV_2
                                ierr = get_ser_var(reshape(&
                                    &X%KV_2(:,:,norm_id(1):norm_id(2),ld),&
                                    &[n_X_loc]),ser_var_loc,scatter=.true.)
                                CHCKERR('')
                                X_ser%KV_2(:,:,:,ld) = &
                                    &reshape(ser_var_loc,n_X)
                            end do
                            
                            do ld = 1,size(X%PV_1,4)
                                ! PV_1
                                ierr = get_ser_var(reshape(&
                                    &X%PV_1(:,:,norm_id(1):norm_id(2),ld),&
                                    &[n_X_loc]),ser_var_loc,scatter=.true.)
                                CHCKERR('')
                                X_ser%PV_1(:,:,:,ld) = &
                                    &reshape(ser_var_loc,n_X)
                                
                                ! KV_1
                                ierr = get_ser_var(reshape(&
                                    &X%KV_1(:,:,norm_id(1):norm_id(2),ld),&
                                    &[n_X_loc]),ser_var_loc,scatter=.true.)
                                CHCKERR('')
                                X_ser%KV_1(:,:,:,ld) = &
                                    &reshape(ser_var_loc,n_X)
                            end do
                        else
                           X_ser%PV_0 = X%PV_0
                           X_ser%PV_1 = X%PV_1
                           X_ser%PV_2 = X%PV_2
                           X_ser%KV_0 = X%KV_0
                           X_ser%KV_1 = X%KV_1
                           X_ser%KV_2 = X%KV_2
                        end if
                        
                        ! interpolate
                        ierr = interp_V_spline(X_ser%PV_0,grid_X_ser%r_F,&
                            &X_sol%PV_0,grid_sol%loc_r_F)
                        CHCKERR('')
                        ierr = interp_V_spline(X_ser%PV_1,grid_X_ser%r_F,&
                            &X_sol%PV_1,grid_sol%loc_r_F)
                        CHCKERR('')
                        ierr = interp_V_spline(X_ser%PV_2,grid_X_ser%r_F,&
                            &X_sol%PV_2,grid_sol%loc_r_F)
                        CHCKERR('')
                        ierr = interp_V_spline(X_ser%KV_0,grid_X_ser%r_F,&
                            &X_sol%KV_0,grid_sol%loc_r_F)
                        CHCKERR('')
                        ierr = interp_V_spline(X_ser%KV_1,grid_X_ser%r_F,&
                            &X_sol%KV_1,grid_sol%loc_r_F)
                        CHCKERR('')
                        ierr = interp_V_spline(X_ser%KV_2,grid_X_ser%r_F,&
                            &X_sol%KV_2,grid_sol%loc_r_F)
                        CHCKERR('')
                        
                        ! clean up
                        call grid_X_trim%dealloc()
                        call grid_X_ser%dealloc()
                        call X_ser%dealloc()
                end select
                    
                ! clean up
                deallocate(r_F_X,r_F_sol)
                read(*,*)
            
                call lvl_ud(-1)
            case (2)                                                            ! solution
                ! user output
                call writo('Copy the perturbation variables to solution grid')
                call lvl_ud(1)
                
                ! copy
                call X%copy(mds_sol,grid_sol,X_sol)
                
                call lvl_ud(-1)
        end select
        
        ! solve the system
        call writo('Solving the system')
        call lvl_ud(1)
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(mds_sol,grid_sol,X_sol,vac,sol)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('System solved')
        
        ! write solution variables to output
        ierr = print_output_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl)
        CHCKERR('')
        
        ! calculate Richardson extrapolation factors if necessary
        call calc_rich_ex(sol%val)
        
        
        ! clean up
        call X_sol%dealloc()
    end function run_driver_sol
    
    !> \public  Interpolate  tensorial  perturbation  quantities  in  the  third
    !! dimension.
    !!
    !! The input grid should not be divided, whereas the output grid can be.
    !!
    !! \note The size  of \c r_i has  to be equal to  the size of \c  V_i in the
    !! third dimension,  and similarly for \c  r_o and \c V_o.  Also, the fourth
    !! dimension of \c V_i and \c K_i has to match. This is not checked.
    integer function interp_V_spline(V_i,r_i,V_o,r_o) result(ierr)
        use num_vars, only: norm_disc_prec_sol, iu
        use bspline_sub_module, only: db1ink, db1val, get_status_message
        
        character(*), parameter :: rout_name = 'interp_V_spline'
        
        ! input / output
        complex(dp), intent(in) :: V_i(:,:,:,:)                                 !< tensorial perturbation variable on input grid
        real(dp), intent(in) :: r_i(:)                                          !< points at which \c V_i is tabulated
        complex(dp), intent(inout) :: V_o(:,:,:,:)                              !< interpolated tensorial perturbation variable
        real(dp), intent(in) :: r_o(:)                                          !< points at which \c V_i is interpolated
        
        ! local variables
        integer :: id, jd, kd, ld, md                                           ! counters
        integer :: spline_init(2)                                               ! spline initialization parameter
        integer :: n(3)                                                         ! size of variables
        real(dp) :: V_loc(2)                                                    ! local PV or KV value
        real(dp), allocatable :: spline_knots(:,:)                              ! knots of spline
        real(dp), allocatable :: spline_coeff(:,:)                              ! coefficients of spline
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize
        n = shape(V_i(:,:,:,1))
        allocate(spline_coeff(n(3),2))
        allocate(spline_knots(n(3)+norm_disc_prec_sol,2))
        do ld = 1,size(V_i,4)
            do jd = 1,n(2)
                do id = 1,n(1)
                    call db1ink(r_i,n(3),rp(V_i(id,jd,:,ld)),&
                        &norm_disc_prec_sol,0,spline_knots(:,1),&
                        &spline_coeff(:,1),ierr)
                    err_msg = get_status_message(ierr)
                    CHCKERR(err_msg)
                    call db1ink(r_i,n(3),ip(V_i(id,jd,:,ld)),&
                        &norm_disc_prec_sol,0,spline_knots(:,2),&
                        &spline_coeff(:,2),ierr)
                    err_msg = get_status_message(ierr)
                    CHCKERR(err_msg)
                    spline_init = 1
                    do kd = 1,size(r_o)
                        do md = 1,2
                            call db1val(r_o(kd),0,spline_knots(:,md),n(3),&
                                &norm_disc_prec_sol,spline_coeff(:,md),&
                                &V_loc(md),ierr,spline_init(md))
                            err_msg = get_status_message(ierr)
                            CHCKERR(err_msg)
                        end do
                        V_o(id,jd,kd,ld) = V_loc(1) + iu*V_loc(2)
                    end do
                end do
            end do
        end do
    end function interp_V_spline
end module driver_sol
