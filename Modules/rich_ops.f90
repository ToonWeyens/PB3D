!------------------------------------------------------------------------------!
!   Operations concerning Richardson extrapolation                             !
!------------------------------------------------------------------------------!
module rich_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_it_rich, tol_rich
    use rich_vars

    implicit none
    private
    public init_rich, term_rich, start_rich_lvl, stop_rich_lvl, do_rich, &
        &calc_rich_ex, find_max_lvl_rich
    
contains
    integer function init_rich() result(ierr)
        use num_vars, only: n_sol_requested, rich_restart_lvl
        use PB3D_ops, only: reconstruct_PB3D_rich
        
        character(*), parameter :: rout_name = 'init_rich'
        
        ! local variables
        integer :: rd                                                           ! counter
        
        ! set variables
        rich_lvl = 1
        rich_conv = .false.
        no_guess = .true.
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        allocate(sol_val_rich(max_it_rich,max_it_rich,n_sol_requested))
        sol_val_rich = 0.0_dp
        allocate(x_axis_rich(max_it_rich,n_sol_requested))
        allocate(max_rel_err(max_it_rich-1))
        allocate(loc_max_rel_err(max_it_rich-1,1))
        
        ! set variables
        if (max_it_rich.gt.1) then                                              ! only when more than one level
            if (rich_restart_lvl.eq.1) then
                ! user output
                call writo('Starting Richardson extrapolation loop')
                call lvl_ud(1)
            else
                ! user output
                call writo('Resuming Richardson extrapolation loop at level '&
                    &//trim(i2str(rich_restart_lvl)))
                call lvl_ud(1)
                
                ! fast-forward counting variables
                rich_lvl = rich_restart_lvl
                if (rich_restart_lvl.ge.2) then
                    n_par_X = min_n_par_X
                    do rd = 2,rich_restart_lvl-1
                        n_par_X = 2 * n_par_X - 1
                    end do
                end if
                
                ! reconstruct Richardson variables
                ierr = reconstruct_PB3D_rich('rich')
                CHCKERR('')
            end if
            
            call writo('Maximum number of iterations: '//&
                &trim(i2str(max_it_rich)))
            call writo('Tolerance requested: '//trim(r2str(tol_rich)))
            
            call writo('')
            call lvl_ud(-1)
        end if
    end function init_rich
    
    integer function term_rich() result(ierr)
        use num_vars, only: rank,norm_disc_prec_sol
        
        character(*), parameter :: rout_name = 'term_rich'
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        if (max_it_rich.gt.1) &
            &call writo('Finishing Richardson extrapolation loop')
        
        call lvl_ud(1)
        
        ! print Richardson variables
        ierr = print_output_rich('rich')
        CHCKERR('')
        
        if (max_it_rich.gt.1) then                                              ! only when more than one level
            ! user output
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
            call writo('Tolerance requested: '//trim(r2str(tol_rich)))
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
        end if
        
        call lvl_ud(-1)
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
    end function term_rich
    
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
        x_axis_rich(rich_lvl,:) = 1.0_dp*n_par_X
        
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
        ! Decides  whether  the  difference  between the  approximation  of  the
        ! Eigenvalues  in  this Richardson  level  and  the previous  Richardson
        ! level, with the same order of the error, falls below a threshold
        ! [ADD SOURCE]
        subroutine check_conv()
            ! local variables
            real(dp), allocatable :: corr(:)                                    ! correction for order (rich_lvl-1) of error
            integer :: ir                                                       ! counter
            
            if (rich_lvl.gt.1) then                                             ! only do this if in Richardson level higher than 1
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
    
    ! Print Richardson variables to an output file:
    !   sol_val_rich, x_axis_rich, max_rel_err, loc_max_rel_err
    integer function print_output_rich(data_name) result(ierr)
        use num_vars, only: PB3D_name, n_sol_requested, rank, no_output
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type, &
            &max_dim_var_1D
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'print_output_rich'
        
        ! input / output
        character(len=*), intent(in) :: data_name                               ! name under which to store
        
        ! local variables
        type(var_1D_type), allocatable, target :: rich_1D(:)                    ! 1D equivalent of X variables
        type(var_1D_type), pointer :: rich_1D_loc => null()                     ! local element in rich_1D
        integer :: id                                                           ! counters
        logical :: no_output_loc                                                ! local copy of no_output
        
        ! initialize ierr
        ierr = 0
        
        ! no output if only one Richardson level
        if (max_it_rich.eq.1) no_output_loc = no_output
        if (max_it_rich.eq.1) no_output = .true.
        
        ! only master does this
        if (rank.eq.0) then
            ! user output
            call writo('Writing richardson variables to output file')
            call lvl_ud(1)
            
            ! Set up the 1D equivalents of the perturbation variables
            allocate(rich_1D(max_dim_var_1D))
            
            ! set up variables rich_1D
            id = 1
            
            ! global_vars
            rich_1D_loc => rich_1D(id); id = id+1
            rich_1D_loc%var_name = 'global_vars'
            allocate(rich_1D_loc%tot_i_min(1),rich_1D_loc%tot_i_max(1))
            allocate(rich_1D_loc%loc_i_min(1),rich_1D_loc%loc_i_max(1))
            rich_1D_loc%tot_i_min = [1]
            rich_1D_loc%tot_i_max = [2]
            rich_1D_loc%loc_i_min = rich_1D_loc%tot_i_min
            rich_1D_loc%loc_i_max = rich_1D_loc%tot_i_max
            allocate(rich_1D_loc%p(2))
            rich_1D_loc%p = [max_it_rich,n_sol_requested]
            
            ! RE_sol_val_rich
            rich_1D_loc => rich_1D(id); id = id+1
            rich_1D_loc%var_name = 'RE_sol_val_rich'
            allocate(rich_1D_loc%tot_i_min(3),rich_1D_loc%tot_i_max(3))
            allocate(rich_1D_loc%loc_i_min(3),rich_1D_loc%loc_i_max(3))
            rich_1D_loc%tot_i_min = [1,1,1]
            rich_1D_loc%tot_i_max = [max_it_rich,max_it_rich,n_sol_requested]
            rich_1D_loc%loc_i_min = rich_1D_loc%tot_i_min
            rich_1D_loc%loc_i_max = rich_1D_loc%tot_i_max
            allocate(rich_1D_loc%p(size(sol_val_rich)))
            rich_1D_loc%p = reshape(realpart(sol_val_rich),[size(sol_val_rich)])
            
            ! IM_sol_val_rich
            rich_1D_loc => rich_1D(id); id = id+1
            rich_1D_loc%var_name = 'IM_sol_val_rich'
            allocate(rich_1D_loc%tot_i_min(3),rich_1D_loc%tot_i_max(3))
            allocate(rich_1D_loc%loc_i_min(3),rich_1D_loc%loc_i_max(3))
            rich_1D_loc%tot_i_min = [1,1,1]
            rich_1D_loc%tot_i_max = [max_it_rich,max_it_rich,n_sol_requested]
            rich_1D_loc%loc_i_min = rich_1D_loc%tot_i_min
            rich_1D_loc%loc_i_max = rich_1D_loc%tot_i_max
            allocate(rich_1D_loc%p(size(sol_val_rich)))
            rich_1D_loc%p = reshape(imagpart(sol_val_rich),[size(sol_val_rich)])
            
            ! x_axis_rich
            rich_1D_loc => rich_1D(id); id = id+1
            rich_1D_loc%var_name = 'x_axis_rich'
            allocate(rich_1D_loc%tot_i_min(2),rich_1D_loc%tot_i_max(2))
            allocate(rich_1D_loc%loc_i_min(2),rich_1D_loc%loc_i_max(2))
            rich_1D_loc%tot_i_min = [1,1]
            rich_1D_loc%tot_i_max = [max_it_rich,n_sol_requested]
            rich_1D_loc%loc_i_min = rich_1D_loc%tot_i_min
            rich_1D_loc%loc_i_max = rich_1D_loc%tot_i_max
            allocate(rich_1D_loc%p(size(x_axis_rich)))
            rich_1D_loc%p = reshape(x_axis_rich,[size(x_axis_rich)])
            
            if (max_it_rich.gt.1) then
                ! max_rel_err
                rich_1D_loc => rich_1D(id); id = id+1
                rich_1D_loc%var_name = 'max_rel_err'
                allocate(rich_1D_loc%tot_i_min(1),rich_1D_loc%tot_i_max(1))
                allocate(rich_1D_loc%loc_i_min(1),rich_1D_loc%loc_i_max(1))
                rich_1D_loc%tot_i_min = [1]
                rich_1D_loc%tot_i_max = [max_it_rich-1]
                rich_1D_loc%loc_i_min = rich_1D_loc%tot_i_min
                rich_1D_loc%loc_i_max = rich_1D_loc%tot_i_max
                allocate(rich_1D_loc%p(size(max_rel_err)))
                rich_1D_loc%p = max_rel_err
                
                ! loc_max_rel_err
                rich_1D_loc => rich_1D(id); id = id+1
                rich_1D_loc%var_name = 'loc_max_rel_err'
                allocate(rich_1D_loc%tot_i_min(2),rich_1D_loc%tot_i_max(2))
                allocate(rich_1D_loc%loc_i_min(2),rich_1D_loc%loc_i_max(2))
                rich_1D_loc%tot_i_min = [1,1]
                rich_1D_loc%tot_i_max = [max_it_rich-1,1]
                rich_1D_loc%loc_i_min = rich_1D_loc%tot_i_min
                rich_1D_loc%loc_i_max = rich_1D_loc%tot_i_max
                allocate(rich_1D_loc%p(size(loc_max_rel_err)))
                rich_1D_loc%p = reshape(loc_max_rel_err*1._dp,&
                    &[size(loc_max_rel_err)])
            end if
            
            ! write
            ierr = print_HDF5_arrs(rich_1D(1:id-1),PB3D_name,trim(data_name))
            CHCKERR('')
            
            ! clean up
            nullify(rich_1D_loc)
            
            ! user output
            call lvl_ud(-1)
            call writo('Richardson variables written to output')
        end if
        
        ! reset no_output
        if (max_it_rich.eq.1) no_output = no_output_loc
        
        ! wait for all processes
        ierr = wait_MPI()
        CHCKERR('')
    end function print_output_rich
end module rich_ops
