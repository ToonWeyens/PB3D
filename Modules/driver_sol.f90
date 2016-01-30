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
        use rich, only: n_r_sol
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use VMEC, only: dealloc_VMEC
        use HELENA, only: dealloc_HEL
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        use utilities, only: test_max_memory
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use MPI_utilities, only: wait_MPI
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_and_calc_grid_sol, &
            &print_output_grid
        use sol_ops, only: print_output_sol
        use rich, only: rich_info_short, calc_rich_ex
        !!use utilities, only: calc_aux_utilities
#if ldebug
        use num_vars, only: iu, use_pol_flux_F
        use grid_utilities, only: get_norm_interp_data
        use utilities, only: c, con
        use MPI_utilities, only: get_ser_var
#endif
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        type(grid_type) :: grid_sol                                             ! solution grid
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_2_type) :: X                                                     ! tensorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
#if ldebug
        type(grid_type) :: grid_X_B                                             ! field-aligned perturbation grid
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
        
        ! read PB3D output file
        ierr = read_PB3D(.false.,.true.,.true.,.false.,.true.,.false.,.true.,&
            &.false.)                                                           ! read the equilibrium and tensorial perturbation variables
        CHCKERR('')
        
        ! reconstruct PB3D variables, using sol limits for X grid
#if ldebug
        ierr = reconstruct_PB3D(.false.,.true.,.true.,.false.,.true.,.false.,&
            &.true.,.false.,grid_eq=grid_eq,grid_X=grid_X,grid_X_B=grid_X_B,&
            &eq=eq,met=met,X_2=X,X_limits=sol_limits,use_setup_nm_X=.false.)
        CHCKERR('')
#else
        ierr = reconstruct_PB3D(.false.,.true.,.true.,.false.,.true.,.false.,&
            &.true.,.false.,grid_eq=grid_eq,grid_X=grid_X,eq=eq,met=met,X_2=X,&
            &X_limits=sol_limits,use_setup_nm_X=.false.)
        CHCKERR('')
#endif
        
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
        ierr = calc_rich_ex(sol%val)
        CHCKERR('')
        
        ! clean up
        call writo('Cleaning up')
        call lvl_ud(1)
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_grid(grid_sol)
        call dealloc_eq(eq)
        ! deallocate variables depending on equilibrium style
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call dealloc_VMEC
            case (2)                                                            ! HELENA
                call dealloc_HEL
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
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
