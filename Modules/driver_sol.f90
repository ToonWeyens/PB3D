!------------------------------------------------------------------------------!
!> Driver of the solution part of PB3D.
!------------------------------------------------------------------------------!
module driver_sol
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type
    use X_vars, only: X_2_type, modes_type, &
        &mds_X, mds_sol
    use vac_vars, only: vac_type
    use sol_vars, only: sol_type
    
    implicit none
    private
    public run_driver_sol
    
    ! global variables
#if ldebug
    logical :: debug_run_driver_sol = .false.                                   !< debug information for run_driver_sol \ldebug
#endif
    
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
            &n_procs, X_grid_style, jump_to_sol, BC_style
        use grid_vars, only: n_r_sol, n_alpha
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_sol
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_grid_sol, &
            &print_output_grid, redistribute_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_lvl
        use rich_ops, only: calc_rich_ex
        use vac_ops, only: calc_vac_res, print_output_vac
        use grid_utilities, only: trim_grid
        use X_ops, only: setup_modes, redistribute_output_X
        use sol_utilities, only: interp_V
#if ldebug
        use X_ops, only: print_debug_X_2
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in), target :: grid_X                           !< perturbation grid
        type(grid_type), intent(inout) :: grid_sol                              !< solution grid
        type(X_2_type), intent(inout) :: X                                      !< integrated tensorial perturbation variables, freed once copied
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        type(sol_type), intent(inout) :: sol                                    !< solution variables
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed solution grid
        type(grid_type) :: grid_X_rdst                                          ! redistributed perturbation grid
        type(grid_type) :: grid_X_sol                                           ! perturbation grid with solution normal part
        type(X_2_type) :: X_rdst                                                ! redistributed X
        type(X_2_type) :: X_sol                                                 ! interpolated X
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
        logical :: do_vac_ops                                                   ! whether specific calculations for vacuum are necessary
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up solution grid if first level
        if (rich_lvl.eq.rich_restart_lvl) then
            ! Divide solution grid under group processes, calculating the limits
            ! and the normal coordinate.
            allocate(r_F_sol(n_r_sol))
            ierr = calc_norm_range('PB3D_sol',sol_limits=sol_limits,&
                &r_F_sol=r_F_sol)
            CHCKERR('')
            
            if (rich_lvl.eq.1) then
                call writo('Set up solution grid')
                call lvl_ud(1)
                
                call writo('Calculate the grid')
                call lvl_ud(1)
                ierr = setup_grid_sol(grid_eq,grid_X,grid_sol,r_F_sol,&
                    &sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                call writo('Write to output file')
                call lvl_ud(1)
                ierr = print_output_grid(grid_sol,'solution','sol',&
                    &remove_previous_arrs=(jump_to_sol.and.&
                    &(X_grid_style.eq.1 .or. X_grid_style.eq.3)))
                CHCKERR('')
                call lvl_ud(-1)
                
                ! set up modes
                ierr = setup_modes(mds_sol,grid_eq,grid_sol,plot_name='sol')
                CHCKERR('')
                
                call lvl_ud(-1)
            else
                ! restore solution grid and trim it
                ierr = reconstruct_PB3D_grid(grid_sol,'sol',&
                    &grid_limits=sol_limits,divided=n_procs.gt.1)
                CHCKERR('')
                ierr = trim_grid(grid_sol,grid_sol_trim)
                CHCKERR('')
                
                ! set up modes
                ierr = setup_modes(mds_sol,grid_eq,grid_sol,plot_name='sol')
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
        if (BC_style(2).eq.1) do_vac_ops = .false.                              ! fixed-boundary BC does not use the vacuum response
        
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
        CHCKERR('')
        grid_X_sol%r_F = grid_sol%r_F
        grid_X_sol%r_E = grid_sol%r_E
        
        
        select case (X_grid_style)
            case (1,3)                                                          ! equilibrium or enriched
                ! user output
                call writo('Redistribute and interpolate the perturbation &
                    &variables to solution grid')
                call lvl_ud(1)
                
                ! redistribute grid and X variables
                ierr = redistribute_output_grid(grid_X,grid_X_rdst)
                CHCKERR('')
                ierr = redistribute_output_X(mds_X,grid_X,grid_X_rdst,X,X_rdst)
                CHCKERR('')
                
                ! X is not needed any more: the next Richardson level reads it
                ! from the output file
                call X%dealloc()
                
                ! initialize
                call X_sol%init(mds_sol,grid_X_sol,is_field_averaged=.true.)
                
                ! interpolate
                ierr = interp_V(mds_X,grid_X_rdst,X_rdst,mds_sol,grid_X_sol,&
                    &X_sol)
                CHCKERR('')
                
                ! clean up
                call X_rdst%dealloc()
                call grid_X_rdst%dealloc()
                
                call lvl_ud(-1)
            case (2)                                                            ! solution
                ! user output
                call writo('Copy the perturbation variables to solution grid')
                call lvl_ud(1)
                
                ! copy
                call X%copy(mds_sol,grid_X,X_sol)
                
                ! X is not needed any more: the next Richardson level reads it
                ! from the output file
                call X%dealloc()
                
                call lvl_ud(-1)
        end select
        
#if ldebug
        ! write integrated  field-aligned tensorial  perturbation quantities
        ! to output and plot
        if (debug_run_driver_sol) then
            ierr = print_debug_X_2(mds_sol,grid_X_sol,X_sol)
            CHCKERR('')
        end if
#endif
        
        ! clean up
        call grid_X_sol%dealloc()
        
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
        ierr = print_output_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl,&
            &remove_previous_arrs=&
            &(jump_to_sol.and.(X_grid_style.eq.1 .or. X_grid_style.eq.3)))
        CHCKERR('')
        
        ! calculate Richardson extrapolation factors if necessary
        call calc_rich_ex(sol%val)
        
        ! clean up
        call X_sol%dealloc()
    end function run_driver_sol
end module driver_sol
