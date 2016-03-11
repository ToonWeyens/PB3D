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
    use eq_vars, only: eq_type
    use met_vars, only: met_type
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
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    integer function run_driver_sol() result(ierr)
        use num_vars, only: EV_style, eq_style
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        use utilities, only: test_max_memory
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq, reconstruct_PB3D_X_2
        use MPI_utilities, only: wait_MPI
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_and_calc_grid_sol, &
            &print_output_grid
        use sol_ops, only: print_output_sol
        use rich_vars, only: rich_info_short, &
            &n_r_sol
        use rich_ops, only: calc_rich_ex
        use files_ops, only: dealloc_in
        !!use utilities, only: calc_aux_utilities
#if ldebug
        use num_vars, only: iu, use_pol_flux_F
        use utilities, only: c, con
        use MPI_utilities, only: get_ser_var
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type), target :: grid_X                                       ! perturbation grid
        type(grid_type) :: grid_sol                                             ! solution grid
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_2_type) :: X                                                     ! tensorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
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
        ierr = reconstruct_PB3D_in()                                            ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! Divide solution grid under group processes, calculating the limits
        ! and the normal coordinate.
        allocate(r_F_sol(n_r_sol))
        ierr = calc_norm_range(sol_limits=sol_limits,r_F_sol=r_F_sol)
        CHCKERR('')
        
        ! reconstruct PB3D variables, using sol limits for X grid
        ! user output
        call writo('Reconstructing PB3D output on output grid')
        call lvl_ud(1)
        ierr = reconstruct_PB3D_grid(grid_eq,'grid_eq')
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'grid_X'//trim(rich_info_short()),&
            &grid_limits=sol_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_eq(grid_eq,eq,met)
        CHCKERR('')
        ierr = reconstruct_PB3D_X_2(grid_X,X,X_limits=sol_limits)
        CHCKERR('')
#if ldebug
        ! need field-aligned perturbation grid as well
        select case (eq_style)
            case (1)                                                            ! VMEC
                grid_X_B => grid_X
            case (2)                                                            ! HELENA
                allocate(grid_X_B)
                ierr = reconstruct_PB3D_grid(grid_X_B,'grid_X_B'//&
                    &trim(rich_info_short()),grid_limits=sol_limits)
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
#endif
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! create solution grid with division limits and setup normal
        ! coordinate
        call writo('Setting up solution grid')
        call lvl_ud(1)
        ierr = setup_and_calc_grid_sol(grid_eq,grid_sol,eq,r_F_sol,&
            &sol_limits)
        CHCKERR('')
        deallocate(r_F_sol)
        call lvl_ud(-1)
        
        ! solve the system
        call writo('Solving the system')
        call lvl_ud(1)
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(grid_sol,X,sol)
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
            
            call writo('The output should be compared with the POST-&
                &output')
            
            call lvl_ud(-1)
        end if
#endif
        
        ! write perturbation grid variables to output
        ierr = print_output_grid(grid_sol,'solution',&
            &'sol'//trim(rich_info_short()))
        CHCKERR('')
        
        ! write solution variables to output
        ierr = print_output_sol(grid_sol,sol)
        CHCKERR('')
        
        ! calculate Richardson extrapolation factors if necessary
        call calc_rich_ex(sol%val)
        
        ! clean up
        call writo('Cleaning up')
        call lvl_ud(1)
        ierr = dealloc_in()
        CHCKERR('')
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_grid(grid_sol)
        call dealloc_eq(eq)
        call dealloc_met(met)
        call dealloc_X(X)
        call dealloc_sol(sol)
        call lvl_ud(-1)
        call writo('Clean')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
    end function run_driver_sol
end module driver_sol
