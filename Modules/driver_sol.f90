!------------------------------------------------------------------------------!
!   Driver of the solution part of PB3D.                                       !
!------------------------------------------------------------------------------!
module driver_sol
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use X_vars, only: X_2_type
    use sol_vars, only: sol_type
    
    implicit none
    private
    public run_driver_sol
#if ldebug
    public debug_sol_grid
#endif
    
    ! global variables
#if ldebug
    logical :: debug_sol_grid = .false.                                         ! plot debug information for treatment of sol grid
    logical :: debug_X_norm = .false.                                           ! plot debug information X_norm
#endif
    
contains
    ! Main driver of PB3D solution part.
    ! sets up:
    !   - grid_sol (only first Richardson level)
    !   - sol
    ! writes to HDF5:
    !   - grid_sol (only first Richardson level)
    !   - sol
    ! deallocates:
    !   - sol before setting up (but after guess)
    integer function run_driver_sol(grid_X,grid_X_B,grid_sol,X,sol) result(ierr)
        use num_vars, only: EV_style, eq_style, rich_restart_lvl
        use grid_vars, only: n_r_sol
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_sol
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_grid_sol, print_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_lvl
        use rich_ops, only: calc_rich_ex
        !!use num_utilities, only: calc_aux_utilities
#if ldebug
        use num_vars, only: iu, use_pol_flux_F
        use num_utilities, only: c, con
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! input / output
        type(grid_type), intent(in), target :: grid_X                           ! perturbation grid
        type(grid_type), intent(inout), pointer :: grid_X_B                     ! field-aligned perturbation grid
        type(grid_type), intent(inout) :: grid_sol                              ! solution grid
        type(X_2_type), intent(in) :: X                                         ! integrated tensorial perturbation variables
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
#if ldebug
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle theta_F or zeta_F
        complex(dp), allocatable :: X_norm(:,:,:)                               ! |X|^2 or other results to be plotted
        integer :: m, k                                                         ! counters
        integer :: kd                                                           ! counter
        character(len=max_str_ln), allocatable :: var_names(:)                  ! names of variable to be plot
        character(len=max_str_ln) :: file_name                                  ! file name
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set up whether Richardson level has to be appended to the name
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
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
                ierr = setup_grid_sol(grid_X,grid_sol,sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                call writo('Write to output file')
                call lvl_ud(1)
                ierr = print_output_grid(grid_sol,'solution','sol')
                CHCKERR('')
                call lvl_ud(-1)
                
                call lvl_ud(-1)
            else
                ! restore solution grid and previous solution
                ierr = reconstruct_PB3D_grid(grid_sol,'sol',&
                    &grid_limits=sol_limits)
                CHCKERR('')
                ierr = reconstruct_PB3D_sol(grid_sol,sol,'sol',&
                    &rich_lvl=rich_lvl-1)
                CHCKERR('')
            end if
            
            deallocate(r_F_sol)
        end if
        
        ! solve the system
        call writo('Solving the system')
        call lvl_ud(1)
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(grid_X,grid_sol,X,sol)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('System solved')
        
#if ldebug
        ! calculate |X|^2 directly from solution vector
        if (debug_X_norm) then
            call writo('Calculating |X|^2 directly from solution')
            call lvl_ud(1)
            
            ! allocate variables
            allocate(X_norm(grid_X_B%n(1),grid_X_B%n(2),grid_X_B%loc_n_r))
            allocate(var_names(6))
            
            ! set pointers
            if (use_pol_flux_F) then
                ang_par_F => grid_X_B%theta_F
            else
                ang_par_F => grid_X_B%zeta_F
            end if
            
            ! calculate |X|^2
            var_names = 'X_norm'
            X_norm = 0._dp
            do kd = 1,grid_sol%loc_n_r
                do m = 1,sol%n_mod
                    do k = 1,sol%n_mod
                        X_norm(:,:,kd) = X_norm(:,:,kd) + &
                            &conjg(sol%vec(k,kd,1))*sol%vec(m,kd,1)* &
                            &exp(iu*(X%m_1(kd,k)-X%m_2(kd,m))*ang_par_F(:,:,kd))
                    end do
                end do
            end do
            
            ! plot |X|^2
            file_name = 'TEST_X_norm_PB3D'
            call plot_HDF5(var_names(1),file_name,rp(X_norm),&
                &tot_dim=grid_X_B%n,loc_offset=[0,0,grid_X%i_min-1])
            
            ! clean up
            nullify(ang_par_F)
            
            call writo('The output should be compared with the POST-output')
            
            call lvl_ud(-1)
        end if
#endif
        
        ! write solution variables to output
        ierr = print_output_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl)
        CHCKERR('')
        
        ! calculate Richardson extrapolation factors if necessary
        call calc_rich_ex(sol%val)
    end function run_driver_sol
end module driver_sol
