!------------------------------------------------------------------------------!
!   Operations concerning Richardson extrapolation                             !
!------------------------------------------------------------------------------!
module rich_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_it_rich, tol_rich
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type
    use sol_vars, only: sol_type
    use rich_vars

    implicit none
    private
    public init_rich, term_rich, start_rich_lvl, stop_rich_lvl, do_rich, &
        &calc_rich_ex, find_max_lvl_rich
    
contains
    integer function init_rich() result(ierr)
        use num_vars, only: n_sol_requested, rich_restart_lvl, eq_style, &
            &rank
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_sol
        use X_ops, only: setup_nm_X
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        use MPI_utilities, only: wait_MPI
        use input_utilities, only: dealloc_in
        
        character(*), parameter :: rout_name = 'init_rich'
        
        ! local variables
        type(grid_type) :: grid_eq_B                                            ! equilibrium grid
        type(grid_type) :: grid_X_B                                             ! field-aligned perturbation grid
        type(eq_1_type) :: eq                                                   ! flux equilibrium variables
        type(sol_type) :: sol                                                   ! solution variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=4) :: grid_eq_name                                        ! name of equilibrium grid
        character(len=4) :: grid_X_name                                         ! name of perturbation grid
        
        ! set variables
        rich_lvl = 1
        rich_conv = .false.
        use_guess = .false.
        
        ! initialize ierr
        ierr = 0
        
        ! user output (only master has variables set up)
        if (rank.eq.0 .and. max_it_rich.gt.1) then                              ! only when more than one level
            if (rich_restart_lvl.eq.1) then
                ! user output
                call writo('Starting Richardson extrapolation loop')
            else
                ! user output
                call writo('Resuming Richardson extrapolation loop at level '&
                    &//trim(i2str(rich_restart_lvl)))
            end if
        end if
        call lvl_ud(1)
        
        ! some preliminary things
        ierr = wait_MPI()
        CHCKERR('')
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! allocate variables
        allocate(sol_val_rich(max_it_rich,max_it_rich,n_sol_requested))
        sol_val_rich = 0.0_dp
        allocate(x_axis_rich(max_it_rich,n_sol_requested))
        allocate(max_rel_err(max_it_rich-1))
        allocate(loc_max_rel_err(max_it_rich-1,1))
        
        ! set variables
        if (max_it_rich.gt.1) then                                              ! only when more than one level
            if (rich_restart_lvl.gt.1) then
                ! user output
                call writo('reconstructing Richardson extrapolation variables')
                call lvl_ud(1)
                
                call writo('Prepare general variables')
                call lvl_ud(1)
                
                ! use guess
                use_guess = .true.
                
                ! set up grid name
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        grid_X_name = 'X'                                       ! grids are already field-aligned
                        grid_eq_name = 'eq'                                     ! grids are already field-aligned
                    case (2)                                                    ! HELENA
                        grid_X_name = 'X_B'                                     ! field-aligned grid is separate
                        grid_eq_name = 'eq_B'                                   ! field-aligned grid is separate
                end select
                
                ! reconstruct field-aligned equilibrium grid (for level 1)
                ierr = reconstruct_PB3D_grid(grid_eq_B,trim(grid_eq_name),&
                    &rich_lvl=1)
                CHCKERR('')
                
                ! reconstruct flux equilibrium variables
                ierr = reconstruct_PB3D_eq_1(grid_eq_B,eq,'eq_1')
                CHCKERR('')
                
                call lvl_ud(-1)
                
                ! loop  over all  previous Richardson  levels using  rich_lvl as
                ! counter
                do rich_lvl = 1,rich_restart_lvl-1
                    call writo('Prepare variables for Richardson level '//&
                        &trim(i2str(rich_lvl)))
                    call lvl_ud(1)
                    
                    ! reconstruct field-aligned perturbation grid
                    ierr = reconstruct_PB3D_grid(grid_X_B,trim(grid_X_name),&
                        &rich_lvl=rich_lvl)
                    CHCKERR('')
                    
                    if (rich_lvl.eq.1) then
                        ! set up n and m variables
                        ierr = setup_nm_X(grid_eq_B,grid_X_B,eq)
                        CHCKERR('')
                    end if
                    
                    ! reconstruct solution
                    ierr = reconstruct_PB3D_sol(grid_X_B,sol,'sol',&
                        &rich_lvl=rich_lvl)
                    CHCKERR('')
                    
                    call lvl_ud(-1)
                    
                    call writo('Recreate Richardson variables')
                    call lvl_ud(1)
                    
                    ! set x_axis_rich
                    if (rich_lvl.eq.1) then
                        ! set up n_par_X for first rich_lvl
                        n_par_X = grid_X_B%n(1)
                        
                        ! check it
                        if (n_par_X.eq.min_n_par_X) then
                            ! set x_axis_rich
                            x_axis_rich(1,:) = min_n_par_X
                        else
                            ierr = 1
                            err_msg = 'Saved variables started with &
                                &min_n_par_X = '//trim(i2str(n_par_X))//&
                                &' whereas here '//trim(i2str(min_n_par_X))//&
                                &' is specified'
                            CHCKERR(err_msg)
                        end if
                    else
                        ! update n_par_X
                        n_par_X = n_par_X*2-1
                        
                        ! set x_axis_rich
                        x_axis_rich(rich_lvl,:) = x_axis_rich(rich_lvl-1,:)*2-1
                    end if
                    
                    ! calculate Richardson extrapolation
                    call calc_rich_ex(sol%val)
                    
                    ! update error variables
                    call check_conv()
                    
                    ! clean up
                    call dealloc_grid(grid_X_B)
                    call dealloc_sol(sol)
                    
                    call lvl_ud(-1)
                end do
                
                ! clean up
                call dealloc_grid(grid_eq_B)
                call dealloc_eq(eq)
                
                ! set rich_lvl
                rich_lvl = rich_restart_lvl
                
                call lvl_ud(-1)
                call writo('Richardson extrapolation variables reconstructed')
            end if
            
            call writo('Maximum number of iterations: '//&
                &trim(i2str(max_it_rich)))
            call writo('Tolerance requested: '//trim(r2strt(tol_rich)))
            
            call writo('')
        end if
        
        ! clean up
        call dealloc_in()
        
        call lvl_ud(-1)
    end function init_rich
    
    subroutine term_rich()
        use num_vars, only: rank,norm_disc_prec_sol
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! user output
        if (max_it_rich.gt.1) then
            ! user output
            call writo('Finishing Richardson extrapolation loop')
            call lvl_ud(1)
            
            if (rich_conv) then
                call writo('Convergence reached in '//trim(i2str(rich_lvl-1))//&
                    &' steps:')
            else
                call writo('After '//trim(i2str(rich_lvl-1))//&
                    &' steps, no convergence reached:')
            end if
            
            call lvl_ud(1)
            
            call writo('Maximum relative error '//&
                &trim(r2str(max_rel_err(rich_lvl-2)))//&
                &' for Eigenvalue '//trim(i2str(loc_max_rel_err(rich_lvl-2,1))))
            call writo('Tolerance requested: '//trim(r2strt(tol_rich)))
            call writo('Resulting best guesses for Eigenvalues:')
            call lvl_ud(1)
            do id = 1,size(sol_val_rich,3)
                call writo('For Eigenvalue '//trim(i2str(id))//': '//&
                    &trim(c2str(sol_val_rich(rich_lvl-1,rich_lvl-1,id)))//',')
            end do
            call lvl_ud(-1)
            call writo('with, theoretically, an error of the order O(Δ^'//&
                &trim(i2str(2*norm_disc_prec_sol*rich_lvl))//'),')
            call writo('compared to O(Δ^'//trim(i2str(2*norm_disc_prec_sol))//&
                &') without Richardson extrapolation')
            call writo('Notice that this is ONLY valid if the problem in the &
                &different Richardson stages is similar:')
            call lvl_ud(1)
            call writo('- It has to be part of the point spectrum')
            call writo('- The Eigenvectors should look similar')
            call lvl_ud(-1)
            
            call lvl_ud(-1)
            
            ! draw X values
            if (rank.eq.0) call draw_sol_val_rich()
            
            call writo('')
            
            call lvl_ud(-1)
        end if
    contains
        ! Draws  the   Eigenvalues  for  the  different   levels  of  Richardson
        ! extrapolation as a function of the number of parallel points.
        subroutine draw_sol_val_rich()
            use num_vars, only: n_sol_requested
            
            ! local variables
            character(len=max_str_ln) :: plot_title                             ! title for plots
            character(len=max_str_ln) :: plot_name                              ! name of plot
            
            ! user output
            call writo('Plotting Eigenvalues as function of nr. of parallel &
                &points')
            call lvl_ud(1)
            
            ! print, using rich_lvl-1 as it has been aumented
            plot_title = 'Eigenvalues as function of nr. of parallel points'
            plot_name = 'Eigenvalues_richardson'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(sol_val_rich(1:rich_lvl-1,1,:)),&
                &x=x_axis_rich(1:rich_lvl-1,:),&
                &draw=.false.)
            
            ! output in file
            call draw_GP(plot_title,plot_name,plot_name,&
                &n_sol_requested,1,.false.)
            
            ! user output
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end subroutine draw_sol_val_rich
    end subroutine term_rich
    
    ! if this Richardson level should be done
    logical function do_rich()
        if (rich_lvl.gt.max_it_rich .or. rich_conv) then
            do_rich = .false.
        else
            do_rich = .true.
        end if
    end function do_rich
    
    ! start a Richardson level
    subroutine start_rich_lvl
        ! Calculate number of parallel points for the solution in Richardson loops
        if (rich_lvl.eq.1) then
            n_par_X = min_n_par_X
        else
            n_par_X = 2 * n_par_X - 1
        end if
        
        ! set use_guess to .false. if user sets no_guess
        if (no_guess) use_guess = .false.
    end subroutine start_rich_lvl
    
    ! stop a Richardson level
    subroutine stop_rich_lvl
        ! update the x axis of the Eigenvalue plot
        x_axis_rich(rich_lvl,:) = 1._dp*n_par_X
        
        ! Richardson extrapolation
        if (max_it_rich.gt.1) then                                              ! only do this if more than 1 Richardson level
            ! user output
            if (rich_lvl.gt.1) call writo('Richardson level '//&
                &trim(i2str(rich_lvl))//' summary')
            call lvl_ud(1)
            
            ! decide whether converged or not
            call check_conv()
            
            ! setup possible guess for next Richardson level
            call set_guess()
            
            if (rich_lvl.gt.1) call writo('')
            call lvl_ud(-1)
        else                                                                    ! if not, Richardson is done
            rich_conv = .true.
        end if
        
        ! increase level
        rich_lvl = rich_lvl + 1
    contains
        ! Decides whether  a guess should be used in  a possible next Richardson
        ! level.
        subroutine set_guess()
            ! decide on guess
            if (rich_lvl.ge.3) then                                             ! start from level 3
                use_guess = max_rel_err(rich_lvl-1).lt.&
                    &max_rel_err(rich_lvl-2)                                    ! set guess for next level to true if max_rel_err is decreasing
            else
                use_guess = .true.                                              ! for first Richardson level, set guess for next level to true
            end if
        end subroutine set_guess
    end subroutine stop_rich_lvl
    
    ! Calculates  the   coefficients  of  the  Eigenvalues   in  the  Richardson
    ! extrapolation. This is done using the recursive formula:
    !   sol_val_rich(lvl,ir,:) = sol_val_rich(lvl,ir-1,:) +  1/(2^(2 p ir) - 1)*
    !       (sol_val_rich(lvl,ir-1,:) - sol_val_rich(lvl-1,ir-1,:)),
    ! where p is norm_disc_prec_sol, the order of the numerical discretization.
    ! as described in [ADD REF]
    subroutine calc_rich_ex(sol_val)
        use num_vars, only: norm_disc_prec_sol
        
        ! input / output
        complex(dp), intent(in) :: sol_val(:)                                   ! EV for this Richardson level
        
        ! local variables
        integer :: ir                                                           ! counter
        
        ! do calculations if rich_lvl > 1
        sol_val_rich(rich_lvl,1,1:size(sol_val)) = sol_val                      ! size(sol_val) can be less than n_sol_requested
        do ir = 2,rich_lvl
            sol_val_rich(rich_lvl,ir,1:size(sol_val)) = &
                &sol_val_rich(rich_lvl,ir-1,1:size(sol_val)) + &
                &1._dp/(2**(2*(ir-1)*norm_disc_prec_sol)-1._dp) * &
                &(sol_val_rich(rich_lvl,ir-1,1:size(sol_val)) - &
                &sol_val_rich(rich_lvl-1,ir-1,1:size(sol_val)))
        end do
    end subroutine calc_rich_ex
    
    ! Decides  whether   the  difference   between  the  approximation   of  the
    ! Eigenvalues in  this Richardson level  and the previous  Richardson level,
    ! with the same order of the error, falls below a threshold
    ! [ADD SOURCE]
    subroutine check_conv()
        ! local variables
        real(dp), allocatable :: corr(:)                                        ! correction for order (rich_lvl-1) of error
        integer :: ir                                                           ! counter
        
        if (rich_lvl.gt.1) then                                                 ! only do this if in Richardson level higher than 1
            ! set relative correction
            allocate(corr(size(sol_val_rich,3)))
            corr = abs(2*(sol_val_rich(rich_lvl,rich_lvl-1,:)-&
                &sol_val_rich(rich_lvl-1,rich_lvl-1,:))/&
                &(sol_val_rich(rich_lvl,rich_lvl-1,:)+&
                &sol_val_rich(rich_lvl-1,rich_lvl-1,:)))
            
            ! get maximum and location of maximum for relative correction
            max_rel_err(rich_lvl-1) = maxval(corr)
            loc_max_rel_err(rich_lvl-1,:) = maxloc(corr)
            
            ! user output
            do ir = 1,rich_lvl-1
                call writo('Richardson level '//trim(i2str(ir))//' -> '//&
                    &trim(i2str(ir+1))//': maximum relative error '//&
                    &trim(r2strt(max_rel_err(ir)))//' for Eigenvalue '//&
                    &trim(i2str(loc_max_rel_err(ir,1))))
            end do
            
            ! check whether tolerance has been reached
            if (max_rel_err(rich_lvl-1).lt.tol_rich) then
                rich_conv = .true.
                call writo('tolerance '//trim(r2strt(tol_rich))//&
                    &' reached after '//trim(i2str(rich_lvl))//&
                    &' iterations')
            else
                call writo('tolerance '//trim(r2strt(tol_rich))//' not yet &
                    &reached')
            end if
        end if
    end subroutine check_conv
    
    ! Probe to find out which Richardson levels are available.
    integer function find_max_lvl_rich(max_lvl_rich_file) result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: probe_HDF5_group
        
        character(*), parameter :: rout_name = 'find_max_lvl_rich'
        
        ! input / output
        integer, intent(inout) :: max_lvl_rich_file                             ! max. Richardson level found in file
        
        ! local variables
        integer :: ir                                                           ! counter
        character(len=max_str_ln) :: group_name                                 ! name of group to probe for
        logical :: group_exists                                                 ! whether probed group exists
        
        ! initialize ierr
        ierr = 0
        
        ! try openining solution without Richardson extrapolation
        group_name = 'sol'
        ierr = probe_HDF5_group(PB3D_name,group_name,group_exists)
        CHCKERR('')
        if (group_exists) then                                                  ! No Richardson extrapolation
            max_lvl_rich_file = 1
        else
            ! try opening solutions for different Richardson level
            ir = 1                                                              ! initialize counter
            group_exists = .true.                                               ! group_exists becomes stopping criterion
            do while (group_exists)
                group_name = 'sol_R_'//trim(i2str(ir))
                ierr = probe_HDF5_group(PB3D_name,group_name,group_exists)
                CHCKERR('')
                ir = ir + 1                                                     ! increment counter
            end do
            max_lvl_rich_file = ir-2                                            ! -2 because there will be one additional iteration
        end if
    end function find_max_lvl_rich
end module rich_ops
