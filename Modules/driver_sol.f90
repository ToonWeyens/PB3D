!------------------------------------------------------------------------------!
!   Driver of the solution part of PB3D.                                       !
!------------------------------------------------------------------------------!
module driver_sol
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
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
    integer function run_driver_sol() result(ierr)
        use num_vars, only: EV_style, eq_style
        use grid_vars, only: dealloc_grid, &
            &n_r_sol
        use eq_vars, only: dealloc_eq
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        use num_utilities, only: test_max_memory
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_2, reconstruct_PB3D_sol
        use MPI_utilities, only: wait_MPI
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_grid_sol, print_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_lvl, use_guess
        use rich_ops, only: calc_rich_ex
        use input_utilities, only: dealloc_in
        !!use num_utilities, only: calc_aux_utilities
#if ldebug
        use num_vars, only: iu, use_pol_flux_F
        use num_utilities, only: c, con
        use MPI_utilities, only: get_ser_var
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type), target :: grid_X                                       ! perturbation grid
        type(grid_type) :: grid_sol                                             ! solution grid
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium
        type(eq_2_type) :: eq_2                                                 ! metric equilibrium
        type(X_2_type) :: X                                                     ! field-averged tensorial perturbation variables (only 1 par dim)
        type(sol_type) :: sol                                                   ! solution variables
        type(sol_type) :: sol_prev                                              ! previous solution variables
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
#if ldebug
        type(grid_type), pointer :: grid_X_B                                    ! field-aligned perturbation grid
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle theta_F or zeta_F
        complex(dp), allocatable :: X_norm(:,:,:)                               ! |X|^2 or other results to be plotted
        integer :: m, k                                                         ! counters
        integer :: kd                                                           ! counter
        character(len=max_str_ln), allocatable :: var_names(:)                  ! names of variable to be plot
        character(len=max_str_ln) :: file_name                                  ! file name
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! some preliminary things
        ierr = wait_MPI()
        CHCKERR('')
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
        ! set up whether Richardson level has to be appended to the name
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! Divide solution grid under group processes, calculating the limits and
        ! the normal coordinate.
        allocate(r_F_sol(n_r_sol))
        ierr = calc_norm_range(sol_limits=sol_limits,r_F_sol=r_F_sol)
        CHCKERR('')
        
        ! reconstruct PB3D variables, using sol limits for X grid
        ! user output
        call writo('Reconstructing PB3D output on output grid')
        call lvl_ud(1)
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
            &grid_limits=sol_limits,tot_rich=.true.)
        CHCKERR('')
        if (rich_lvl.gt.1) then                                                 ! also need solution grid
            ierr = reconstruct_PB3D_grid(grid_sol,'sol',grid_limits=sol_limits)
            CHCKERR('')
        end if
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_2(grid_eq,eq_2,'eq_2',&
            &rich_lvl=rich_lvl_name,tot_rich=.true.)
        CHCKERR('')
        ierr = reconstruct_PB3D_X_2(grid_X,X,'X_2_int',rich_lvl=rich_lvl,&
            &X_limits=sol_limits,is_field_averaged=.true.)
        CHCKERR('')
        ! need previous solution as well if guess is used
        if (use_guess .and. rich_lvl.gt.1) then
            ierr = reconstruct_PB3D_sol(grid_X,sol_prev,'sol',&
                &rich_lvl=rich_lvl-1,sol_limits=sol_limits)
            CHCKERR('')
        end if
#if ldebug
        ! need field-aligned perturbation grid as well for debugging
        if (debug_X_norm) then
            select case (eq_style)
                case (1)                                                        ! VMEC
                    grid_X_B => grid_X
                case (2)                                                        ! HELENA
                    allocate(grid_X_B)
                    ierr = reconstruct_PB3D_grid(grid_X_B,'X_B',&
                        &rich_lvl=rich_lvl,grid_limits=sol_limits)
                    CHCKERR('')
            end select
        end if
#endif
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! set up solution grid if first level
        if (rich_lvl.eq.1) then
            call writo('Set up solution grid')
            call lvl_ud(1)
            
            call writo('Calculate the grid')
            call lvl_ud(1)
            ierr = setup_grid_sol(grid_eq,grid_sol,r_F_sol,sol_limits)
            CHCKERR('')
            call lvl_ud(-1)
            
            call writo('Write to output file')
            call lvl_ud(1)
            ierr = print_output_grid(grid_sol,'solution','sol')
            CHCKERR('')
            call lvl_ud(-1)
            
            call lvl_ud(-1)
        end if
        deallocate(r_F_sol)
        
        ! solve the system
        call writo('Solving the system')
        call lvl_ud(1)
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(grid_sol,X,sol,sol_prev=sol_prev)
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
            call plot_HDF5(var_names(1),file_name,realpart(X_norm),&
                &tot_dim=grid_X_B%n,loc_offset=[0,0,grid_X%i_min-1])
            
            ! clean up
            nullify(ang_par_F)
            if (eq_style.eq.2) then
                call dealloc_grid(grid_X_B)
                deallocate(grid_X_B)
            end if
            nullify(grid_X_B)
            
            call writo('The output should be compared with the POST-output')
            
            call lvl_ud(-1)
        end if
#endif
        
        ! write solution variables to output
        ierr = print_output_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl)
        CHCKERR('')
        
        ! calculate Richardson extrapolation factors if necessary
        call calc_rich_ex(sol%val)
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call dealloc_in()
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_grid(grid_sol)
        call dealloc_eq(eq_1)
        call dealloc_eq(eq_2)
        call dealloc_X(X)
        call dealloc_sol(sol)
        if (use_guess .and. rich_lvl.gt.1) call dealloc_sol(sol_prev)
        call lvl_ud(-1)
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
    end function run_driver_sol
end module driver_sol
