!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type
    
    implicit none
    private
    public run_rich_driver
#if ldebug
    public debug_X_grid
#endif
    
    ! global variables
#if ldebug
    logical :: debug_X_grid = .false.                                           ! plot debug information for treatment of X grid
#endif
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_rich_driver() result(ierr)
        use num_vars, only: min_alpha, max_alpha, n_alpha, glb_rank, grp_nr, &
            &max_alpha, alpha_job_nr, use_pol_flux_F, eq_style, group_output, &
            &output_name, n_groups, max_it_r, plot_flux_q
        use MPI_ops, only: split_MPI, merge_MPI, get_next_job
        use MPI_utilities, only: wait_MPI
        use eq_vars, only: create_eq, dealloc_eq
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X, min_r_X, &
            &max_r_X, min_n_r_X
        use eq_ops, only: calc_flux_q, flux_q_plot
        use VMEC, only: dealloc_VMEC
        use HELENA, only: dealloc_HEL
        use grid_ops, only: calc_eqd_grid, setup_grid_eq
        use grid_vars, only: dealloc_grid, &
            &n_r_eq, n_par_X, min_par_X, max_par_X
        use HDF5_ops, only: create_output_HDF5
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: eq_limits(2)                                                 ! min. and max. index of eq. grid of this process
        real(dp), allocatable :: alpha(:)                                       ! all field line labes all
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(eq_type) :: eq                                                     ! equilibrium for this alpha
        
        ! initialize ierr
        ierr = 0
        
        ! set up flux_name
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        
        ! user output
        call writo('The calculations will be done')
        call lvl_ud(1)
        
        if (max_it_r.eq.1) then
            call writo('for '//trim(i2str(min_n_r_X))//' values on &
                &normal range '//trim(r2strt(min_r_X))//'..'//&
                &trim(r2strt(max_r_X)))
        else
            call writo('for minimally '//trim(i2str(min_n_r_X))//' values on &
                &normal range '//trim(r2strt(min_r_X))//'..'//&
                &trim(r2strt(max_r_X)))
        end if
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
        if (n_alpha.eq.1) then
            call writo('for alpha = '//trim(r2strt(min_alpha*pi)))
        else
            call writo('for '//trim(i2str(n_alpha))//&
                &' values of alpha '//trim(r2strt(min_alpha*pi))//'..'//&
                &trim(r2strt(max_alpha*pi)))
        end if
        call writo('using the '//trim(flux_name)//' flux as the normal &
            &variable')
        if (use_pol_flux_F) then
            call writo('with toroidal mode number n = '//trim(i2str(min_n_X)))
            call writo('and poloidal mode number m = '//trim(i2str(min_m_X))//&
                &'..'//trim(i2str(max_m_X)))
        else
            call writo('with poloidal mode number m = '//trim(i2str(min_m_X)))
            call writo('and toroidal mode number n = '//trim(i2str(min_n_X))//&
                &'..'//trim(i2str(max_n_X)))
        end if
        
        call lvl_ud(-1)
        
        ! determine the magnetic field lines for which to run the calculations 
        ! (equidistant grid)
        allocate(alpha(n_alpha))
        ierr = calc_eqd_grid(alpha,min_alpha*pi,max_alpha*pi,excl_last=.true.)
        CHCKERR('')
        
        ! split  the  communicator MPI_Comm_world into subcommunicators
        call writo('Set up groups for dynamical load balancing')
        call lvl_ud(1)
        ierr = split_MPI(n_r_eq,eq_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! setup equilibrium grid
        ierr = setup_grid_eq(grid_eq,eq_limits)
        CHCKERR('')
        
        ! create equilibrium
        call writo('Initialize equilibrium quantities')
        ierr = create_eq(grid_eq,eq)
        CHCKERR('')
        
        ! calculate flux quantities and complete equilibrium grid
        call writo('Calculate flux quantities')
        ierr = calc_flux_q(eq,grid_eq)
        CHCKERR('')
        
        ! do the calculations for every  field line, where initially every group
        ! is assigned a field line. On completion of a job, the completing group
        ! inquires  all the other  groups using passive  1-sided MPI to  get the
        ! next group that has to be done
        ! loop over all fied lines
        field_lines: do
            ! non-global group masters cannot output (any more)
            group_output = .false.
            
            ! get next job for current group
            ierr = get_next_job(alpha_job_nr)
            CHCKERR('')
            
            ! Do the calculations for a field line alpha
            if (alpha_job_nr.gt.0) then
                ! non-global group masters can output
                group_output = .true.
                
                ! display message
                call writo('Job '//trim(i2str(alpha_job_nr))//': Calculations &
                    &for field line alpha = '//&
                    &trim(r2strt(alpha(alpha_job_nr)))//&
                    &', allocated to group '//trim(i2str(grp_nr+1)))
                
                call lvl_ud(1)                                                  ! starting calculation for current fied line
                
                ! open HDF5 file for output
                ierr = create_output_HDF5()
                CHCKERR('')
                
                ! calculate for this alpha job
                ierr = run_rich_driver_for_alpha(grid_eq,eq,alpha(alpha_job_nr))
                CHCKERR('')
                
                ! display message
                call lvl_ud(-1)                                                 ! done with calculation for current field line
                call writo('Job '//trim(i2str(alpha_job_nr))//&
                    &': Calculations for field line alpha = '//&
                    &trim(r2strt(alpha(alpha_job_nr)))//&
                    &', completed by group '//trim(i2str(grp_nr+1)))
            else
                exit field_lines
            end if
        end do field_lines
        
        ! plot flux quantities if requested
        if (plot_flux_q .and. grp_nr.eq.0) then                                 ! only first group because it is the same for all the groups
            ierr = flux_q_plot(eq,grid_eq)
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        
        ! deallocate variables
        call dealloc_eq(eq)
        call dealloc_grid(grid_eq)
        
        ! choose which equilibrium style is being used:
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
        
        ! wait on all the groups
        ierr = wait_MPI(2)
        CHCKERR('')
        
        ! merge  the subcommunicator into communicator MPI_Comm_world
        if (glb_rank.eq.0) then
            if (n_alpha.gt.1) then
                call writo('Stability analysis concluded at all '//&
                    &trim(i2str(n_alpha))//' fieldlines')
            else
                call writo('Stability analysis concluded')
            end if
            if (n_groups.eq.2) then
                call writo('The outputs of the other group is saved in "'//&
                    &trim(output_name)//'_2"')
            else if (n_groups.gt.2) then
                call writo('The outputs of the other groups 2..'//&
                    &trim(i2str(n_alpha))//' i are saved in "'//&
                    &trim(output_name)//'_2..'//trim(i2str(n_alpha))//'"')
            end if
            call writo('Merging groups for dynamical load balancing back &
                &together')
        end if
        call lvl_ud(1)
        ierr = merge_MPI()
        CHCKERR('')
        call lvl_ud(-1)
    end function run_rich_driver
    
    ! Runs the  calculations for one  of the alpha's.
    integer function run_rich_driver_for_alpha(grid_eq,eq,alpha) result(ierr)
        use num_vars, only: n_sol_requested, max_it_r, no_guess, &
            &rich_lvl_nr, grp_rank, alpha_job_nr, eq_style
        use X_vars, only: dealloc_X
        use eq_ops, only: calc_eq, print_output_eq
        use X_ops, only: solve_EV_system, calc_magn_ints, prepare_X, &
            &print_output_X, print_output_sol
        use met_vars, only: dealloc_met
        use MPI_ops, only: divide_X_grid
        use grid_vars, only: create_grid, dealloc_grid
        use grid_ops, only: coord_F2E, setup_and_calc_grid_B, calc_ang_grid_eq
        use vac, only: calc_vac
#if ldebug
        use grid_ops, only: trim_grid, untrim_grid
#endif
        
        character(*), parameter :: rout_name = 'run_rich_driver_for_alpha'
        
        ! input / output
        type(grid_type), intent(inout), target :: grid_eq                       ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium for this alpha
        real(dp), intent(in) :: alpha                                           ! alpha at which to run the calculations
        
        ! local variables
        type(grid_type) :: grid_X                                               ! perturbation grid
        type(grid_type), pointer :: grid_eq_B => null()                         ! field-aligned equilibrium grid
        type(met_type) :: met                                                   ! metric variables
        type(X_type) :: X, X_B                                                  ! general and field-aligned perturbation variables
        integer :: n_r_X                                                        ! total number of normal points in X grid for this step in Richardson ex.
        integer :: X_limits(2)                                                  ! min. and max. index of X grid for this process
        real(dp), allocatable :: r_X(:)                                         ! normal points in Flux coords.
        logical :: done_richard                                                 ! is it converged?
        complex(dp), allocatable :: X_val_rich(:,:,:)                           ! Richardson array of eigenvalues X val
        real(dp), allocatable :: x_axis(:,:)                                    ! x axis for plot of Eigenvalues with n_r_X
        logical :: use_guess                                                    ! whether a guess is formed from previous level of Richardson
        integer :: n_sol_found                                                  ! how many solutions found and saved
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! name of plot
        character(len=max_str_ln) :: draw_ops                                   ! optional drawing options
#if ldebug
        type(grid_type) :: grid_eq_trim                                         ! trimmed equilibrium grid
        type(grid_type) :: grid_eq_ghost                                        ! ghosted equilibrium grid
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! calculate angular grid points for equilibrium grid
        ierr = calc_ang_grid_eq(grid_eq,eq,alpha)
        CHCKERR('')
        
        ! calculate the non-flux equilibrium quantities for current alpha
        ierr = calc_eq(grid_eq,eq,met)
        CHCKERR('')
        
        ! prepare matrix elements
        ierr = prepare_X(grid_eq,eq,met,X)
        CHCKERR('')
        
        ! calculate vacuum response
        ierr = calc_vac(X)
        CHCKERR('')
        
        ! set up field-aligned equilibrium grid
        ierr = setup_and_calc_grid_B(grid_eq,grid_eq_B,eq,alpha)
        CHCKERR('')
        
        ! write equilibrium variables to output
        ierr = print_output_eq(grid_eq,grid_eq_B,eq,met,alpha)
        CHCKERR('')
        
        ! deallocate metric variables
        call dealloc_met(met)
        
        ! write X variables to output
        ierr = print_output_X(grid_eq_B,X)
        CHCKERR('')
        
        ! adapt variables to a field-aligned grid
        ierr = adapt_X_to_B(grid_eq,grid_eq_B,X,X_B)
        CHCKERR('')
        
        ! deallocate non-aligned perturbation variables
        call dealloc_X(X)
        
        ! calculate magnetic integrals
        ierr = calc_magn_ints(grid_eq_B,X_B)
        CHCKERR('')
        
        ! Initalize some variables for Richardson loop
        rich_lvl_nr = 1
        done_richard = .false.
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            allocate(X_val_rich(1:max_it_r,1:max_it_r,1:n_sol_requested))
            X_val_rich = 0.0_dp
            allocate(x_axis(1:max_it_r,1:n_sol_requested))
        end if
        
        ! Start Richardson loop
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Starting perturbation calculation in field-aligned &
                &grid with Richardson extrapolation')
        else
            call writo('Starting perturbation calculation in field-aligned &
                &grid')
        end if
        call lvl_ud(1)                                                          ! before richardson loop
        
        ! initialize use_guess for first Richardson loop
        use_guess = .true.
        
        Richard: do while (.not.done_richard .and. rich_lvl_nr.le.max_it_r)
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call writo('Level ' // trim(i2str(rich_lvl_nr)) // &
                    &' of Richardson extrapolation')
                call lvl_ud(1)                                                  ! beginning of one richardson loop
            end if
            
            ! setting up perturbation grid
#if ldebug
            if (.not.debug_X_grid) then
#endif
                ! calculate  number of  radial points  for  the  perturbation in
                ! Richardson loops and save in n_r_X
                call writo('calculating the normal points')
                call lvl_ud(1)
                call calc_n_r_X(rich_lvl_nr,n_r_X)
                call lvl_ud(-1)
                
                ! divide  perturbation grid  under group  processes, calculating
                ! the limits and the normal coordinate
                ierr = divide_X_grid(n_r_X,X_limits,r_X)
                CHCKERR('')
                
                ! create perturbation grid with division limits and setup normal
                ! coordinate
                call writo('creating perturbation grid')
                call lvl_ud(1)
                ierr = create_grid(grid_X,n_r_X,X_limits)
                CHCKERR('')
                grid_X%r_F = r_X
                grid_X%grp_r_F = r_X(X_limits(1):X_limits(2))
                deallocate(r_X)
                ierr = coord_F2E(grid_eq,eq,grid_X%r_F,grid_X%r_E,&
                    &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
                CHCKERR('')
                ierr = coord_F2E(grid_eq,eq,grid_X%grp_r_F,grid_X%grp_r_E,&
                    &r_F_array=grid_eq%r_F,r_E_array=grid_eq%r_E)
                CHCKERR('')
                call lvl_ud(-1)
#if ldebug
            else
                ! user output
                call writo('for debugging, equating the perturbation grid to &
                    &the equilibrium grid with '//trim(i2str(grid_eq%n(3)))//&
                    &' normal points')
                call lvl_ud(1)
                
                ! trim equilibrium grid
                ierr = trim_grid(grid_eq,grid_eq_trim)
                CHCKERR('')
                
                ! untrim grid
                ierr = untrim_grid(grid_eq_trim,grid_eq_ghost,1)
                CHCKERR('')
                
                ! number of radial points given by equilibrium grid
                n_r_X = grid_eq_ghost%n(3)
                X_limits = [grid_eq_ghost%i_min,grid_eq_ghost%i_max]
                
                ! create grid
                ierr = create_grid(grid_X,n_r_X,X_limits)
                CHCKERR('')
                
                ! fill grid
                if (grid_eq_ghost%divided) then                                 ! if input grid divided, grp_r gets priority
                    grid_X%grp_r_E = grid_eq_ghost%grp_r_E
                    grid_X%grp_r_F = grid_eq_ghost%grp_r_F
                end if
                grid_X%r_E = grid_eq_trim%r_E
                grid_X%r_F = grid_eq_trim%r_F
                
                ! clean up
                call dealloc_grid(grid_eq_trim)
                call dealloc_grid(grid_eq_ghost)
                
                call lvl_ud(-1)
            end if
#endif
            
            ! set use_guess to .false. if no_guess
            if (no_guess) use_guess = .false.
            
            ! setup the matrices of the generalized EV system AX = lambda BX and
            ! solve it
            call writo('treating the EV system')
            call lvl_ud(1)
            ierr = solve_EV_system(grid_eq_B,grid_X,X_B,use_guess,n_sol_found)
            CHCKERR('')
            call lvl_ud(-1)
            
            ! write solution variables to output
            ierr = print_output_sol(grid_X,X_B)
            CHCKERR('')
            
            ! Richardson extrapolation
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                ! update the x axis of the Eigenvalue plot
                x_axis(rich_lvl_nr,:) = 1.0_dp*n_r_X
                
                ! update  the  variable  X_val_rich  with the  results  of  this
                ! Richardson level
                call writo('updating Richardson extrapolation variables')
                call lvl_ud(1)
                ierr = calc_rich_ex(rich_lvl_nr,X_B%val,X_val_rich,&
                    &done_richard,use_guess)
                CHCKERR('')
                call lvl_ud(-1)
            else                                                                ! if not, Richardson is done
                done_richard = .true.
            end if
            
            ! destroy grid if not yet done with Richardson Extrapolation
            if (.not.done_richard) call dealloc_grid(grid_X)
            
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call lvl_ud(-1)                                                 ! end of one richardson loop
            end if
        end do Richard
        
        call lvl_ud(-1)                                                         ! done with richardson
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Finished Richardson loop')
        else
            call writo('Finished perturbation calculation')
        end if
        
        ! visualize  the  Richardson  Extrapolation by  plotting Eigenvalues  as
        ! function of nr. of normal points
        if (max_it_r.gt.1 .and. grp_rank.eq.0) then
            call writo('Plotting Eigenvalues as function of nr. of normal &
                &points in Richardson Extrapolation')
            call lvl_ud(1)
            
            ! set up drawing options
            draw_ops = 'pt 7 ps 0.2'
            
            ! output on screen
            plot_title = 'job '//trim(i2str(alpha_job_nr))//' - Eigenvalues &
                &as function of nr. of normal points'
            plot_name = 'Eigenvalues_A'//trim(i2str(alpha_job_nr))//&
                &'_richardson'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(X_val_rich(1:rich_lvl_nr,1,:)),&
                &x=x_axis(1:rich_lvl_nr,:),draw=.false.)
            
            ! same output in file as well
            call draw_GP(plot_title,plot_name,plot_name,&
                &n_sol_requested,1,.false.,draw_ops=draw_ops)
            
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end if
        
        ! deallocate Richardson loop variables
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            deallocate(X_val_rich)
            deallocate(x_axis)
        end if
        
        ! destroy grid
        call dealloc_grid(grid_X)
        
        ! deallocate field-aligned quantities in Flux coords.
        call dealloc_X(X_B)
        if (eq_style.eq.2) call dealloc_grid(grid_eq_B)                         ! only for HELENA
        nullify(grid_eq_B)
    contains
        ! calculates the number of normal  points for the perturbation n_r_X for
        ! the various Richardson iterations
        ! The aim is to  halve the step size, which is given  by dx(n) = 1/(n-1)
        ! or, inverting: n(dx) = 1 + 1/dx.
        ! This yields n(dx/2)/n(dx) = (2+dx)/(1+dx) = (2n(dx)-1)/n(dx)
        ! The recursion formula is therefore: n(dx/2) = 2n(dx) - 1
        subroutine calc_n_r_X(ir,n_r_X)
            use X_vars, only: min_n_r_X
            
            ! input / output
            integer, intent(in) :: ir
            integer, intent(inout) :: n_r_X
            
            write(*,*) '!!!! TEMPORARILY DISABLED BECAUSE RICH. EXT. NOT &
                &WORKING PROPERLY !!!'
            
            n_r_X = min_n_r_X*ir
            !!!if (ir.eq.1) then
                !!!n_r_X = min_n_r_X
            !!!else
                !!!n_r_X = 2 * n_r_X - 1
            !!!end if
            call writo(trim(i2str(n_r_X))//' normal points for this level...')
        end subroutine calc_n_r_X
        
        ! calculates  the  coefficients of  the  Eigenvalues  in the  Richardson
        ! extrapolation
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
    end function run_rich_driver_for_alpha
    
    ! Adapt  X  variables  angularly  to  a  field-aligned  grid,  depending  on
    ! equilibrium type:
    !   X: J_exp_ang_par_F, U_i, DU_i, PV_i, KV_i
    ! To save memory, the original quantities adapted are deallocated.
    ! Note that this, by definition, does not affect and thus doesn not apply to
    ! flux functions.
    integer function adapt_X_to_B(grid_eq,grid_eq_B,X,X_B) &
        &result(ierr)
        use num_vars, only: eq_style
        use HELENA, only: interp_HEL_on_grid
        use X_vars, only: create_X
        
        character(*), parameter :: rout_name = 'adapt_X_to_B'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq, grid_eq_B                       ! general and field-aligned equilibrium grid
        type(X_type), intent(inout) :: X                                        ! general perturbation variables
        type(X_type), intent(inout) :: X_B                                      ! field-aligned perturbation variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! create new perturbation variables
        call create_X(grid_eq_B,X_B)
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! no conversion necessary: already in field-aligned grid
                X_B%J_exp_ang_par_F = X%J_exp_ang_par_F
                deallocate(X%J_exp_ang_par_F)
                X_B%U_0 = X%U_0
                deallocate(X%U_0); nullify(X%U_0)
                X_B%U_1 = X%U_1
                deallocate(X%U_1); nullify(X%U_1)
                X_B%DU_0 = X%DU_0
                deallocate(X%DU_0); nullify(X%DU_0)
                X_B%DU_1 = X%DU_1
                deallocate(X%DU_1); nullify(X%DU_1)
                X_B%PV_0 = X%PV_0
                deallocate(X%PV_0); nullify(X%PV_0)
                X_B%PV_1 = X%PV_1
                deallocate(X%PV_1); nullify(X%PV_1)
                X_B%PV_2 = X%PV_2
                deallocate(X%PV_2); nullify(X%PV_2)
                X_B%KV_0 = X%KV_0
                deallocate(X%KV_0); nullify(X%KV_0)
                X_B%KV_1 = X%KV_1
                deallocate(X%KV_1); nullify(X%KV_1)
                X_B%KV_2 = X%KV_2
                deallocate(X%KV_2); nullify(X%KV_2)
                X_B%vac_res = X%vac_res
                deallocate(X%vac_res)
            case (2)                                                            ! HELENA
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,X,X_B,&
                    &grid_name='field-aligned grid')
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function adapt_X_to_B
end module driver_rich
