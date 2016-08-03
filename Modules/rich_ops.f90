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
        &calc_rich_ex, find_max_rich_lvl
    
contains
    integer function init_rich() result(ierr)
        use num_vars, only: n_sol_requested, rich_restart_lvl, eq_style, &
            &rank
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_sol
        use X_ops, only: setup_nm_X
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
                call writo('Richardson extrapolation loop initialization')
            else
                ! user output
                call writo('Richardson extrapolation loop continuation at &
                    &level '//trim(i2str(rich_restart_lvl)))
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
                    call grid_X_B%dealloc()
                    call sol%dealloc()
                    
                    call lvl_ud(-1)
                end do
                
                ! clean up
                call grid_eq_B%dealloc()
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
    integer function start_rich_lvl() result(ierr)
        use grid_utilities, only: calc_n_par_X_rich
        use grid_vars, only: n_r_eq
        use eq_utilities, only: divide_eq_jobs
        use HELENA_vars, only: nchi
        use num_vars, only: eq_jobs_lims, eq_style, rank, n_procs, &
            &minim_output, PB3D_name_eq
        use HDF5_ops, only: create_output_HDF5
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'start_rich_lvl'
        
        ! local variables
        logical :: only_half_grid                                               ! calculate only half grid with even points
        integer :: n_par_X_loc                                                  ! local n_par_X
        real(dp) :: tot_mem_size                                                ! total memory size per process
        real(dp), allocatable :: tot_mem_size_full(:)                           ! tot_mem_size for all processes
        character(len=max_str_ln) :: err_msg                                    ! error message
        
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
        
        ! Get local n_par_X to divide the equilibrium jobs. Note that the 
        ! average size per process is n_r_eeq / n_procs
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! divide equilibrium jobs
                ierr = divide_eq_jobs(n_par_X_loc,n_r_eq/n_procs)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! divide equilibrium jobs
                ierr = divide_eq_jobs(n_par_X_loc,n_r_eq/n_procs,&
                    &n_par_X_base=nchi,tot_mem_size=tot_mem_size)               ! always use nchi points as base
                CHCKERR('')
                
                ! check whether  calculation for  first Richardson level  can be
                ! done,  as this  is  the calculation  where  the variables  are
                ! calculated that are interpolated in all levels.
                ! Note: This could be changed, and HELENA could be calculated by
                ! dividing up  the parallel grid.  However, this is very  low on
                ! the list of priorities as  it would required ghost regions and
                ! all those things. For now: use more memory per process or more
                ! processes.
                if (rich_lvl.eq.1 .and. size(eq_jobs_lims,2).gt.1) then
                    ierr = get_ser_var([tot_mem_size],tot_mem_size_full)
                    CHCKERR('')
                    if (rank.eq.0) then
                        ierr = 1
                        err_msg = 'Need at least '//&
                            &trim(i2str(ceiling(sum(tot_mem_size_full))))//&
                            &'MB for HELENA calculations'
                        CHCKERR(err_msg)
                    end if
                end if
        end select
        
        ! set use_guess to .false. if user sets no_guess
        if (no_guess) use_guess = .false.
        
        ! possibly create separate output for eq
        if (minim_output) then
            call writo('Output file size minimized: non-essential variables &
                &get temporary file')
            call lvl_ud(1)
            ierr = create_output_HDF5(PB3D_name_eq)
            CHCKERR('')
            call lvl_ud(-1)
        end if
    end function start_rich_lvl
    
    ! stop a Richardson level
    subroutine stop_rich_lvl
        use num_vars, only: minim_output, PB3D_name_eq
        use files_utilities, only: delete_file
        
        ! local variables
        integer :: istat                                                        ! status
        
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
        
        ! possibly remove separate output for eq
        if (minim_output) istat = delete_file(PB3D_name_eq)
    contains
        ! Decides whether  a guess should be used in  a possible next Richardson
        ! level.
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
    ! Also increases the next tol_SLEPC if it is too low.
    ! [ADD SOURCE]
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
    
    ! Probe to find out which Richardson levels are available.
    ! Also sets the global variable minim_output if provided.
    integer function find_max_rich_lvl(max_lvl_rich_file,minim_output) &
        &result(ierr)
        use num_vars, only: PB3D_name, eq_style
        use HDF5_utilities, only: probe_HDF5_group
        
        character(*), parameter :: rout_name = 'find_max_rich_lvl'
        
        ! input / output
        integer, intent(inout) :: max_lvl_rich_file                             ! max. Richardson level found in file
        logical, intent(inout), optional :: minim_output                        ! minimal output
        
        ! local variables
        integer :: ir                                                           ! counter
        character(len=max_str_ln) :: group_name                                 ! name of group to probe for
        logical :: group_exists(2)                                              ! whether probed group exists
        
        ! initialize ierr
        ierr = 0
        
        ! try openining solution for different Richardson extrapolation levels
        group_exists(1) = .true.                                                ! group_exists becomes stopping criterion
        ir = 1                                                                  ! initialize counter
        do while (group_exists(1))
            group_name = 'sol_R_'//trim(i2str(ir))
            ierr = probe_HDF5_group(PB3D_name,group_name,group_exists(1))
            CHCKERR('')
            ir = ir + 1                                                         ! increment counter
        end do
        max_lvl_rich_file = ir-2                                                ! -2 because there will be one additional iteration
        
        ! set up minim_output
        if (present(minim_output)) then
            select case(eq_style)
                case (1)                                                        ! VMEC
                    ! try  whether eq_2  and X_1 variables  exist with  first eq
                    ! jobs suffix
                    group_name = 'eq_2_R_'//trim(i2str(max_lvl_rich_file))//&
                        &'_E_1'
                    ierr = probe_HDF5_group(PB3D_name,group_name,&
                        &group_exists(1))
                    CHCKERR('')
                    group_name = 'X_1_R_'//trim(i2str(max_lvl_rich_file))//&
                        &'_E_1'
                    ierr = probe_HDF5_group(PB3D_name,group_name,&
                        &group_exists(2))
                    CHCKERR('')
                    if (group_exists(1) .and. group_exists(2)) then
                        minim_output = .false.
                    else
                        minim_output = .true.
                    end if
                case (2)                                                        ! HELENA
                    minim_output = .false.                                      ! never minimal output
            end select
        end if
    end function find_max_rich_lvl
end module rich_ops
