!------------------------------------------------------------------------------!
!> Operations concerning Richardson extrapolation.
!------------------------------------------------------------------------------!
module rich_ops
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
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
        &calc_rich_ex, find_max_rich_lvl
    
contains
    !> Initialize Richardson extrapolation system.
    !!
    !! Needs to be called only once.
    !!
    !! \return ierr
    integer function init_rich() result(ierr)
        use num_vars, only: n_sol_requested, rich_restart_lvl, eq_style, rank
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_sol
        use X_ops, only: init_modes, setup_modes
        use X_vars, only: mds_X, mds_sol
        use grid_vars, only: n_alpha, min_alpha, max_alpha, alpha
        use grid_utilities, only: calc_eqd_grid
        
        character(*), parameter :: rout_name = 'init_rich'
        
        ! local variables
        type(grid_type) :: grid_eq_B                                            !< equilibrium grid
        type(grid_type) :: grid_X_B                                             !< field-aligned perturbation grid
        type(grid_type) :: grid_sol                                             !< solution grid
        type(eq_1_type) :: eq                                                   !< flux equilibrium variables
        type(sol_type) :: sol                                                   !< solution variables
        character(len=4) :: grid_eq_name                                        !< name of equilibrium grid
        character(len=4) :: grid_X_name                                         !< name of perturbation grid
        character(len=max_str_ln) :: err_msg                                    !< error message
        
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
                call writo('Richardson extrapolation loop initialization')
            else
                ! user output
                call writo('Richardson extrapolation loop continuation at &
                    &level '//trim(i2str(rich_restart_lvl)))
            end if
        end if
        call lvl_ud(1)
        
        ! some preliminary things
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! allocate variables
        allocate(sol_val_rich(max_it_rich,max_it_rich,n_sol_requested))
        sol_val_rich = 0.0_dp
        allocate(x_axis_rich(max_it_rich,n_sol_requested))
        allocate(max_rel_err(max_it_rich-1))
        allocate(loc_max_rel_err(max_it_rich-1,1))
        
        ! set up alpha
        allocate(alpha(n_alpha))
        ierr = calc_eqd_grid(alpha,min_alpha*pi,max_alpha*pi,excl_last=.true.)
        CHCKERR('')
        
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
                
                ! set up grid names
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
                
                ! reconstruct solution grid
                ierr = reconstruct_PB3D_grid(grid_sol,'sol')
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
                        ! initialize global n and m variables
                        ierr = init_modes(grid_eq_B,eq)
                        CHCKERR('')
                        
                        ! set up n and m variables for field-aligned X grid
                        ierr = setup_modes(mds_X,grid_eq_B,grid_X_B)
                        CHCKERR('')
                        
                        ! set up n and m variables for sol grid
                        ierr = setup_modes(mds_sol,grid_eq_B,grid_sol)
                        CHCKERR('')
                    end if
                    
                    ! reconstruct solution
                    ierr = reconstruct_PB3D_sol(mds_sol,grid_sol,sol,'sol',&
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
                    call grid_X_B%dealloc()
                    call sol%dealloc()
                    
                    call lvl_ud(-1)
                end do
                
                ! clean up
                call grid_eq_B%dealloc()
                call grid_sol%dealloc()
                call eq%dealloc()
                
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
        
        call lvl_ud(-1)
    end function init_rich
    
    !> Terminate the Richardson extrapolation system.
    !!
    !! Needs to be called only once.
    subroutine term_rich()
        use num_vars, only: rank, norm_disc_prec_sol
        use input_utilities, only: dealloc_in
        
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
                call writo('For Eigenvalue '//trim(i2str(id))//':')
                call lvl_ud(1)
                call writo(trim(c2str(sol_val_rich(rich_lvl-1,rich_lvl-1,id))),&
                    &alert=.true.)
                call lvl_ud(-1)
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
        
        ! clean up
        call dealloc_in()
    contains
        ! Draws  the   Eigenvalues  for  the  different   levels  of  Richardson
        ! extrapolation as a function of the number of parallel points.
        !> \private
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
            call print_ex_2D([plot_title],plot_name,&
                &rp(sol_val_rich(1:rich_lvl-1,1,:)),&
                &x=x_axis_rich(1:rich_lvl-1,:),draw=.false.)
            
            ! output in file
            call draw_ex([plot_title],plot_name,n_sol_requested,1,.false.)
            
            ! user output
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end subroutine draw_sol_val_rich
    end subroutine term_rich
    
    !> Tests whether this Richardson level should be done.
    logical function do_rich()
        if (rich_lvl.gt.max_it_rich .or. rich_conv) then
            do_rich = .false.
        else
            do_rich = .true.
        end if
    end function do_rich
    
    !> Start a Richardson level.
    !!
    !! Calculates \c n_par_X for this level.  Then uses this to divide the jobs.
    !! equilibrium Finally, the limits are set up for these jobs.
    !!
    !! \return ierr
    integer function start_rich_lvl() result(ierr)
        use grid_utilities, only: calc_n_par_X_rich
        use grid_vars, only: n_r_eq, n_r_sol, n_alpha
        use eq_ops, only: divide_eq_jobs, calc_eq_jobs_lims
        use HELENA_vars, only: nchi
        use num_vars, only: eq_style, jump_to_sol, rich_restart_lvl
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'start_rich_lvl'
        
        ! local variables
        integer :: n_div                                                        !< number of divisions in equilibrium jobs
        logical :: only_half_grid                                               !< calculate only half grid with even points
        integer :: n_par_X_loc                                                  !< local n_par_X
        integer :: var_size_without_par(2)                                      !< size without parallel dimension for eq_2 and X_1 variables
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate  total  number  of  parallel  points  for  the  solution  in
        ! Richardson loops
        if (rich_lvl.eq.1) then
            n_par_X = min_n_par_X
            only_half_grid = .false.
        else
            n_par_X = 2 * n_par_X - 1
            only_half_grid = .true.
        end if
        
        ! use calc_n_par_X_rich, taking into account possible half grid
        ierr = calc_n_par_X_rich(n_par_X_loc,only_half_grid)
        CHCKERR('')
        
        if (rich_lvl.eq.rich_restart_lvl .and. jump_to_sol) then
            call writo('Jumping straight to solution, no need to divide into &
                &equilibrium jobs')
            n_div = 1
        else
            ! Get local  n_par_X to divide  the equilibrium jobs. Note  that the
            ! average  size of  eq_2  variables for  all  processes together  is
            ! n_r_eq,  times the  number of  field lines,  the size  due to  the
            ! dimensions  corresponding to  the derivatives,  and the  number of
            ! variables (see  subroutine 'calc_memory_eq'  and 'calc_memory_X'),
            ! while it is n_r_sol for the X_1 variables.
            var_size_without_par(1) = n_r_eq
            var_size_without_par(2) = n_r_sol
            var_size_without_par = var_size_without_par * n_alpha
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! divide equilibrium jobs
                    ierr = divide_eq_jobs(n_par_X_loc,var_size_without_par,&
                        &n_div)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ! divide equilibrium jobs
                    ! Note: calculations  for first Richardson level  have to be
                    ! done with a  single job, as this is  the calculation where
                    ! the variables are calculated  that are interpolated in all
                    ! levels.
                    ierr = divide_eq_jobs(nchi,var_size_without_par,n_div,&
                        &n_div_max=1,range_name='poloidal points in the &
                        &axisymmetric HELENA cross-section')                    ! everything is tabulated on nchi poloidal points
                    CHCKERR('')
            end select
        end if
        
        ! calculate equilibrium job limits
        ierr = calc_eq_jobs_lims(n_par_X_loc,n_div)
        CHCKERR('')
        
        ! set use_guess to .false. if user sets no_guess
        if (no_guess) use_guess = .false.
    end function start_rich_lvl
    
    !> Stop a Richardson level.
    !!
    !! It  decides whether  convergence has  been reached,  and possibly  sets a
    !! guess for the next level.
    subroutine stop_rich_lvl()
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
        !> \private
        subroutine set_guess()
            ! local variables
            real(dp), parameter :: tol = 5._dp                                  ! tolerance for increase of error
            
            ! decide on guess
            if (rich_lvl.ge.3) then                                             ! start from level 3
                use_guess = max_rel_err(rich_lvl-1).lt.&
                    &tol*max_rel_err(rich_lvl-2)                                ! set guess for next level to true if max_rel_err is decreasing
            else
                use_guess = .true.                                              ! for first Richardson level, set guess for next level to true
            end if
        end subroutine set_guess
    end subroutine stop_rich_lvl
    
    !> Calculates  the  coefficients  of   the  Eigenvalues  in  the  Richardson
    !! extrapolation.
    !!
    !! This is done using the recursive formula:
    !!  \f[ s(l,i,:) = s(l,i-1,:) +  \frac{1}{2^{2 \ p \ i} - 1}
    !!       (s(l,i-1,:) - s(l-1,i-1,:)),
    !!  \f]
    !! where  \f$s\f$  =  \c  sol_val_rich,  \f$i\f$  =  \c  ir,  \f$l\f$  =  \c
    !! lvl  and \f$p\f$  = \c  norm_disc_prec_sol,  the order  of the  numerical
    !! discretization. as  described in \cite  dahlquist2003numerical, algorithm
    !! 3.2, p. 306.
    subroutine calc_rich_ex(sol_val)
        use num_vars, only: norm_disc_prec_sol
        
        ! input / output
        complex(dp), intent(in) :: sol_val(:)                                   !< EV for this Richardson level
        
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
    
    !> Decides  whether   convergence  has  been  reached   for  the  Richardson
    !! extrapolation.
    !!
    !! This is done by checking whether the difference between the approximation
    !! of the Eigenvalues  in this Richardson level and  the previous Richardson
    !! level, with the same order of the error, falls below a threshold.
    !!
    !! Also increases the next tol_SLEPC if it is too low.
    !!
    !! \see calc_rich_ex()
    subroutine check_conv()
        use num_vars, only: tol_SLEPC
        
        ! local variables
        real(dp), allocatable :: corr(:)                                        ! correction for order (rich_lvl-1) of error
        integer :: ir                                                           ! counter
        real(dp) :: tol_SLEPC_adapt_factor = 0.01_dp                            ! factor to relate tol_SLEPC to maximum relative error
        
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
            
            ! check whether to decrease next tol_SLEPC
            if (rich_lvl.lt.max_it_rich) tol_SLEPC(rich_lvl+1) = min(&
                &tol_SLEPC(rich_lvl+1),&
                &tol_SLEPC_adapt_factor*max_rel_err(rich_lvl-1))
        end if
    end subroutine check_conv
    
    !> Probe to find out which Richardson levels are available.
    !!
    !! \return ierr
    integer function find_max_rich_lvl(lvl_rich_req,max_lvl_rich_file) &
        &result(ierr)
        use num_vars, only: PB3D_name
        use HDF5_utilities, only: probe_HDF5_group
        
        character(*), parameter :: rout_name = 'find_max_rich_lvl'
        
        ! input / output
        integer, intent(inout) :: lvl_rich_req                                  !< requested Richardson level
        integer, intent(inout) :: max_lvl_rich_file                             !< max. Richardson level found in file
        
        ! local variables
        integer :: ir                                                           ! counter
        character(len=max_str_ln) :: group_name                                 ! name of group to probe for
        logical :: group_exists                                                 ! whether probed group exists
        
        ! initialize ierr
        ierr = 0
        
        ! try openining solution for different Richardson extrapolation levels
        group_exists = .true.                                                   ! group_exists becomes stopping criterion
        ir = 1                                                                  ! initialize counter
        if (lvl_rich_req.gt.0 .and. lvl_rich_req.lt.huge(1)) ir = lvl_rich_req  ! move to requested
        do while (group_exists)
            group_name = 'sol_R_'//trim(i2str(ir))
            ierr = probe_HDF5_group(PB3D_name,group_name,group_exists)
            CHCKERR('')
            ir = ir + 1                                                         ! increment counter
        end do
        max_lvl_rich_file = ir-2                                                ! -2 because there will be one additional iteration
    end function find_max_rich_lvl
end module rich_ops
