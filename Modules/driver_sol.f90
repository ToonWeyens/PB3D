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
#endif
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_driver_sol() result(ierr)
        use num_vars, only: EV_style, rich_lvl_nr, max_it_r, n_sol_requested, &
            &no_guess, rank
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        use utilities, only: test_max_memory
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use MPI_utilities, only: wait_MPI
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_ops, only: calc_norm_range, setup_and_calc_grid_sol
        use sol_ops, only: print_output_sol
        !!use utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_sol'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(grid_type) :: grid_eq, grid_sol                                    ! (field-aligned) equilibrium and solution grid
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_2_type) :: X                                                     ! tensorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
        integer :: n_r_sol                                                      ! nr. of normal points in sol grid for Richardson lvl
        logical :: use_guess                                                    ! whether a guess is formed from previous level of Richardson
        logical :: done_richard                                                 ! is it converged?
        integer :: sol_limits(2)                                                ! min. and max. index of sol grid for this process
        real(dp), allocatable :: x_axis(:,:)                                    ! x axis for plot of Eigenvalues with n_r_sol
        real(dp), allocatable :: r_F_sol(:)                                     ! normal points in solution grid
        complex(dp), allocatable :: X_val_rich(:,:,:)                           ! Richardson array of eigenvalues X val
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! some preliminary things
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('The solution of the generalized Eigenvalue problem is &
            &calculated')
        call lvl_ud(1)
        
        ! read PB3D output file
        ierr = read_PB3D(.false.,.true.,.false.,.true.,.false.)                 ! read the equilibrium and tensorial perturbation variables
        CHCKERR('')
        
        ! reconstruct PB3D variables
        ierr = reconstruct_PB3D(.false.,.true.,.false.,.true.,.false.,&
            &grid_eq=grid_eq,eq=eq,met=met,X_2=X)
        CHCKERR('')
        
        ! deallocate unneeded quantities
        call dealloc_eq(eq)
        call dealloc_met(met)
        
        ! Initalize some variables for Richardson loop
        rich_lvl_nr = 1
        done_richard = .false.
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            allocate(X_val_rich(1:max_it_r,1:max_it_r,1:n_sol_requested))
            X_val_rich = 0.0_dp
            allocate(x_axis(1:max_it_r,1:n_sol_requested))
        end if
        
        ! user output
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Starting perturbation calculation in field-aligned &
                &grid with Richardson extrapolation')
        else
            call writo('Starting perturbation calculation in field-aligned &
                &grid')
        end if
        call lvl_ud(1)                                                          ! before richardson loop
        
        Richard: do while (.not.done_richard .and. rich_lvl_nr.le.max_it_r)
            ! user output
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call writo('Starting Level ' // trim(i2str(rich_lvl_nr)) // &
                    &' of Richardson extrapolation')
                call lvl_ud(1)                                                  ! beginning of one richardson loop
            end if
            
            ! Calculate number of  radial points for the  solution in Richardson
            ! loops and save in n_r_sol.
            call writo('calculating the normal points')
            call lvl_ud(1)
            call calc_n_r_sol(rich_lvl_nr,n_r_sol)
            call lvl_ud(-1)
            
            ! Divide solution grid under group processes, calculating the limits
            ! and the normal coordinate.
            allocate(r_F_sol(n_r_sol))
            ierr = calc_norm_range(sol_limits=sol_limits,r_F_sol=r_F_sol)
            CHCKERR('')
            
            ! create solution grid with division limits and setup normal
            ! coordinate
            call writo('Setting up solution grid')
            call lvl_ud(1)
            ierr = setup_and_calc_grid_sol(grid_eq,grid_sol,eq,r_F_sol,&
                &sol_limits)
            CHCKERR('')
            deallocate(r_F_sol)
            call lvl_ud(-1)
            
            ! set use_guess to .false. if no_guess
            if (no_guess) use_guess = .false.
            
            ! solve the system
            call writo('Solving the system')
            call lvl_ud(1)
            select case (EV_style)
                case(1)                                                         ! SLEPC solver for EV problem
                    ! solve the system
                    ierr = solve_EV_system_SLEPC(grid_eq,grid_sol,X,sol,&
                        &use_guess)
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
            ierr = print_output_sol(grid_sol,sol)
            CHCKERR('')
            
            ! Richardson extrapolation
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                ! update the x axis of the Eigenvalue plot
                x_axis(rich_lvl_nr,:) = 1.0_dp*n_r_sol
                
                ! update  the  variable  X_val_rich  with the  results  of  this
                ! Richardson level
                call writo('updating Richardson extrapolation variables')
                call lvl_ud(1)
                ierr = calc_rich_ex(rich_lvl_nr,sol%val,X_val_rich,&
                    &done_richard,use_guess)
                CHCKERR('')
                call lvl_ud(-1)
            else                                                                ! if not, Richardson is done
                done_richard = .true.
            end if
            
            ! user output
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call lvl_ud(-1)                                                 ! beginning of one richardson loop
                call writo('Completed Level ' // trim(i2str(rich_lvl_nr)) // &
                    &' of Richardson extrapolation')
            end if
            
            ! draw output
            if (max_it_r.gt.1 .and. rank.eq.0) call draw_X_val_rich(X_val_rich)
            
            ! clean up
            call dealloc_grid(grid_eq)
            call dealloc_grid(grid_sol)
            
            ! synchronize MPI
            ierr = wait_MPI()
            CHCKERR('')
        end do Richard
        
        ! user output
        call lvl_ud(-1)                                                         ! done with richardson
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Finished Richardson loop')
        else
            call writo('Finished perturbation calculation')
        end if
        
        ! clean up
        call dealloc_X(X)
        call dealloc_sol(sol)
    contains
        ! Calculates the number  of normal points n_r_sol for  the solution grid
        ! for the various Richardson iterations.
        ! The aim is to  halve the step size, which is given  by dx(n) = 1/(n-1)
        ! or, inverting: n(dx) = 1 + 1/dx.
        ! This yields n(dx/2)/n(dx) = (2+dx)/(1+dx) = (2n(dx)-1)/n(dx)
        ! The recursion formula is therefore: n(dx/2) = 2n(dx) - 1
        subroutine calc_n_r_sol(ir,n_r_sol)
            use X_vars, only: min_n_r_sol
            
            ! input / output
            integer, intent(in) :: ir
            integer, intent(inout) :: n_r_sol
            
            write(*,*) '!!!! TEMPORARILY DISABLED BECAUSE RICH. EXT. NOT &
                &WORKING PROPERLY !!!'
            
            n_r_sol = min_n_r_sol*ir
            !!!if (ir.eq.1) then
                !!!n_r_sol = min_n_r_sol
            !!!else
                !!!n_r_sol = 2 * n_r_sol - 1
            !!!end if
            call writo(trim(i2str(n_r_sol))//' normal points for this level')
        end subroutine calc_n_r_sol
        
        ! Calculates  the  coefficients of  the  Eigenvalues  in the  Richardson
        ! extrapolation.
        ! This is done using the recursive formula
        !   X_val_rich(ir,ir2,:) = X_val_rich(ir,ir2-1,:) +  1/(2^(2ir2) - 1) * 
        !       (X_val_rich(ir,ir2-1,:) - X_val_rich(ir-1,ir2-1,:)),
        ! as described in [ADD REF]
        ! [MPI] All ranks
        integer function calc_rich_ex(ir,X_val,X_val_rich,done_richard,&
                &use_guess_for_next_level) result(ierr)
            use num_vars, only: tol_r
            
            character(*), parameter :: rout_name = 'calc_rich_ex'
            
            ! input / output
            integer, intent(inout) :: ir                                        ! level of Richardson extrapolation (starting at 1)
            complex(dp), intent(in) :: X_val(:)                                 ! EV for this Richardson level
            complex(dp), intent(inout) :: X_val_rich(:,:,:)                     ! extrapolated coefficients
            logical, intent(inout) :: done_richard                              ! if Richardson loop has converged sufficiently
            logical, intent(inout) :: use_guess_for_next_level                  ! if a guessed is used for next Richardson level
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            integer :: ir2                                                      ! counter
            complex(dp), allocatable :: corr(:)                                 ! correction
            real(dp) :: max_corr                                                ! maximum of maximum of correction
            real(dp), save :: prev_max_corr                                     ! max_corr at previous Richardson level
            integer :: loc_max_corr(1)                                          ! location of maximum of correction
            
            ! initialize ierr
            ierr = 0
            write(*,*) '!!! FOR RICHARDSON EXTRAPOLATION YOU HAVE TO FIND &
                &CONSTANT EIGENVALUES !!!'
            
            ! tests
            if (size(X_val_rich,1).ne.size(X_val_rich,2) .or. &
                &ir.gt.size(X_val_rich,1)) then
                ierr = 1
                err_msg = 'X_val_rich has to have correct dimensions'
                CHCKERR(err_msg)
            end if
            
            ! allocate correction
            allocate(corr(size(X_val))); corr = 1.0E15
            
            ! do calculations if ir > 1
            X_val_rich(ir,1,1:size(X_val)) = X_val
            do ir2 = 2,ir
                corr = 1./(2**(2*ir2)-1.) * &
                    &(X_val_rich(ir,ir2-1,1:size(X_val)) - &
                    &X_val_rich(ir-1,ir2-1,1:size(X_val)))
                X_val_rich(ir,ir2,1:size(X_val)) = &
                    &X_val_rich(ir,ir2-1,1:size(X_val)) + corr
            end do
            
            if (ir.gt.1) then                                                   ! only do this if in Richardson level higher than 1
                ! get maximum and location of maximum for relative correction
                max_corr = maxval(abs(corr/X_val_rich(ir,ir,1:size(X_val))))
                loc_max_corr = maxloc(abs(corr/X_val_rich(ir,ir,1:size(X_val))))
                call writo('maximum relative error: '//trim(r2strt(max_corr))//&
                    &' for Eigenvalue '//trim(i2str(loc_max_corr(1))))
                
                ! check whether tolerance has been  reached for last value ir2 =
                ! ir
                if (maxval(abs(corr/X_val_rich(ir,ir,1:size(X_val))))&
                    &.lt.tol_r) then
                    done_richard = .true.
                    call writo('tolerance '//trim(r2strt(tol_r))//&
                        &' reached after '//trim(i2str(ir))//' iterations')
                else
                    call writo('tolerance '//trim(r2strt(tol_r))//' not yet &
                        &reached')
                end if
                
                ! determine use_guess_for_next_level
                use_guess_for_next_level = max_corr.lt.prev_max_corr            ! set guess for next level to true if max_corr is decreasing
                prev_max_corr = max_corr                                        ! save max_corr for next level
            else
                use_guess_for_next_level = .true.                               ! for first Richardson level, set guess for next level to true
            end if
            
            !!! TEMPORARILY !!!
            call writo('!!!!! DO NOT USE GUESS FOR MODIFIED RICHARDSON !!!')
            use_guess_for_next_level = .false.
            
            ! check for convergence
            if (.not.done_richard) then
                if (ir.lt.max_it_r) then                                        ! not yet at maximum Richardson iteration
                    ir = ir + 1
                else                                                            ! maximum nr. of Richardson iterations reached
                    call writo('maximum number of Richardson iterations &
                        &reached')
                    done_richard = .true.
                end if
            end if
        end function calc_rich_ex
        
        ! Draws  the  Eigenvalues   for   the  different  levels  of  Richardson
        ! extrapolation as a function of the number of normal points.
        subroutine draw_X_val_rich(X_val_rich)
            ! input / output
            complex(dp), allocatable :: X_val_rich(:,:,:)                       ! Richardson array of eigenvalues X val
            
            ! local variables
            character(len=max_str_ln) :: plot_title                             ! title for plots
            character(len=max_str_ln) :: plot_name                              ! name of plot
            character(len=max_str_ln) :: draw_ops                               ! optional drawing options
            
            ! user output
            call writo('Plotting Eigenvalues as function of nr. of normal &
                &points in Richardson Extrapolation')
            call lvl_ud(1)
            
            ! set up drawing options
            draw_ops = 'pt 7 ps 0.2'
            
            ! output on screen
            plot_title = 'Eigenvalues as function of nr. of normal points'
            plot_name = 'Eigenvalues_richardson'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(X_val_rich(1:rich_lvl_nr,1,:)),&
                &x=x_axis(1:rich_lvl_nr,:),draw=.false.)
            
            ! same output in file as well
            call draw_GP(plot_title,plot_name,plot_name,&
                &n_sol_requested,1,.false.,draw_ops=draw_ops)
            
            ! user output
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end subroutine draw_X_val_rich
    end function run_driver_sol
end module driver_sol
