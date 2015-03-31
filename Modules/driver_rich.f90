!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
#include <PB3D_macros.h>
    use num_vars, only: max_it_r, dp, pi, max_str_ln
    use str_ops, only: i2str, r2str, r2strt
    use messages, only: writo, print_ar_2, print_ar_1, lvl_ud
    use output_ops, only: print_GP_2D, print_GP_3D, draw_GP, print_HDF5
    
    implicit none
    private
    public run_rich_driver
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_rich_driver() result(ierr)
        use num_vars, only: min_alpha, max_alpha, n_alpha, glb_rank, grp_nr, &
            &max_alpha, alpha_job_nr, use_pol_flux_F, eq_style
        use MPI_ops, only: split_MPI, merge_MPI, get_next_job
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X
        use VMEC, only: dealloc_VMEC_final
        use HELENA, only: dealloc_HEL_final
        use grid_ops, only: calc_eqd_grid
        use grid_vars, only: create_grid, grid_type, &
            &n_r_eq, n_par_X
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: eq_limits(2)                                                 ! min. and max. index of eq. grid of this process
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        real(dp), allocatable :: alpha(:)                                       ! all field line labes all
        
        ! initialize ierr
        ierr = 0
        
        ! set up flux_name
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        
        ! output concerning n_alpha
        call writo('The calculations will be done')
        call lvl_ud(1)
        call writo('using the '//trim(flux_name)//' flux as the normal &
            &variable')
        call writo('for '//trim(i2str(n_alpha))//' values of alpha')
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
        
        ! split  the  communicator MPI_COMM_WORLD into subcommunicators
        call writo('Setting up groups for dynamical load balancing')
        call lvl_ud(1)
        ierr = split_MPI(n_r_eq,eq_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! create equilibrium grid
        call writo('Creating equilibrium grid')
        call lvl_ud(1)
        ierr = create_grid(grid_eq,[n_par_X,1,n_r_eq],eq_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! do the calculations for every  field line, where initially every group
        ! is assigned a field line. On completion of a job, the completing group
        ! inquires  all the other  groups using passive  1-sided MPI to  get the
        ! next group that has to be done
        ! loop over all fied lines
        field_lines: do
            ! get next job for current group
            ierr = get_next_job(alpha_job_nr)
            CHCKERR('')
            
            ! Do the calculations for a field line alpha
            if (alpha_job_nr.gt.0) then
                ! display message
                call writo('Job '//trim(i2str(alpha_job_nr))//': Calculations &
                    &for field line alpha = '//&
                    &trim(r2strt(alpha(alpha_job_nr)))//&
                    &', allocated to group '//trim(i2str(grp_nr)))
                
                call lvl_ud(1)                                                  ! starting calculation for current fied line
                
                ! calculate
                ierr = run_for_alpha(grid_eq,grid_X,alpha(alpha_job_nr))
                CHCKERR('')
                
                ! display message
                call lvl_ud(-1)                                                 ! done with calculation for current field line
                call writo('Job '//trim(i2str(alpha_job_nr))//&
                    &': Calculations for field line alpha = '//&
                    &trim(r2strt(alpha(alpha_job_nr)))//&
                    &', completed by group '//trim(i2str(grp_nr)))
            else
                call writo('Finished all jobs')
                exit field_lines
            end if
        end do field_lines
        
        ! deallocate equilibrium variables
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call dealloc_VMEC_final
            case (2)                                                            ! HELENA
                call dealloc_HEL_final
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! merge  the subcommunicator into communicator MPI_COMM_WORLD
        if (glb_rank.eq.0) then
            call writo('Stability analysis concluded at all '//&
                &trim(i2str(n_alpha))//' fieldlines')
            call writo('Merging groups for dynamical load balancing back &
                &together')
        end if
        call lvl_ud(1)
        ierr = merge_MPI()
        CHCKERR('')
        call lvl_ud(-1)
    end function run_rich_driver
    
    ! Runs the  calculations for one  of the alpha's.
    integer function run_for_alpha(grid_eq,grid_X,alpha) result(ierr)
        use num_vars, only: n_sol_requested, max_it_r, no_guess, &
            &alpha_job_nr, grp_rank
        use eq_vars, only: dealloc_eq_final, eq_type
        use X_vars, only: dealloc_X_final, X_type
        use X_ops, only: solve_EV_system
        use metric_vars, only: dealloc_metric_final, metric_type
        use MPI_ops, only: divide_X_grid
        use grid_vars, only: create_grid, destroy_grid, grid_type
        
        character(*), parameter :: rout_name = 'run_for_alpha'
        
        ! input / output
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        real(dp), intent(in) :: alpha                                           ! alpha at which to run the calculations
        
        ! local variables
        type(eq_type) :: eq                                                     ! equilibrium for this alpha
        type(metric_type) :: met                                                ! metric variables
        type(X_type) :: X                                                       ! perturbation variables
        integer :: ir                                                           ! counter for richardson extrapolation
        integer :: n_r_X                                                        ! total number of normal points in X grid for this step in Richardson ex.
        integer :: X_limits(2)                                                  ! min. and max. index of X grid for this process
        real(dp), allocatable :: grp_r_X(:)                                     ! normal points in Flux coords., globally normalized to (min_r_X..max_r_X)
        logical :: done_richard                                                 ! is it converged?
        complex(dp), allocatable :: X_val_rich(:,:,:)                           ! Richardson array of eigenvalues X val
        real(dp), allocatable :: x_axis(:,:)                                    ! x axis for plot of Eigenvalues with n_r_X
        logical :: use_guess                                                    ! whether a guess is formed from previous level of Richardson
        integer :: n_sol_found                                                  ! how many solutions found and saved
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        
        ! initialize ierr
        ierr = 0
        
        ! calculate necessary quantities on equilibrium grid
        ierr = calculate_vars_in_eq_grid(grid_eq,eq,met,X,alpha)
        CHCKERR('')
        
        ! Initalize some variables for Richardson loop
        ir = 1
        done_richard = .false.
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            allocate(X_val_rich(1:max_it_r,1:max_it_r,1:n_sol_requested))
            X_val_rich = 0.0_dp
            allocate(x_axis(1:max_it_r,1:n_sol_requested))
        end if
        
        ! Start Richardson loop
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Starting perturbation calculation with Richardson &
                &extrapolation')
        else
            call writo('Starting perturbation calculation')
        end if
        call lvl_ud(1)                                                          ! before richardson loop
        
        ! initialize use_guess for first Richardson loop
        use_guess = .true.
        
        Richard: do while (.not.done_richard .and. ir.le.max_it_r)
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call writo('Level ' // trim(i2str(ir)) // &
                    &' of Richardson extrapolation')
                call lvl_ud(1)                                                  ! beginning of one richardson loop
            end if
            
            !  calculate  number  of  radial  points  for  the  perturbation  in
            ! Richardson loops and save in n_r_X
            call writo('calculating the normal points')
            call lvl_ud(1)
            call calc_n_r_X(ir,n_r_X)
            call lvl_ud(-1)
            
            ! divide  perturbation grid  under group processes,  calculating the
            ! limits and the normal coordinate
            ierr = divide_X_grid(n_r_X,X_limits,grp_r_X)
            CHCKERR('')
            
            ! create  perturbation grid  with division  limits and  setup normal
            ! coordinate
            call writo('creating perturbation grid')
            call lvl_ud(1)
            ierr = create_grid(grid_X,n_r_X,X_limits)
            CHCKERR('')
            grid_X%grp_r_F = grp_r_X
            deallocate(grp_r_X)
            call lvl_ud(-1)
            
            ! set use_guess to .false. if no_guess
            if (no_guess) use_guess = .false.
            
            ! setup the matrices of the generalized EV system AX = lambda BX and
            ! solve it
            call writo('treating the EV system')
            call lvl_ud(1)
            ierr = solve_EV_system(grid_eq,grid_X,X,use_guess,n_sol_found)
            CHCKERR('')
            call lvl_ud(-1)
            
            ! Richardson extrapolation
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                ! update the x axis of the Eigenvalue plot
                x_axis(ir,:) = 1.0_dp*n_r_X
                
                ! update  the  variable  X_val_rich  with the  results  of  this
                ! Richardson level
                call writo('updating Richardson extrapolation variables')
                call lvl_ud(1)
                ierr = calc_rich_ex(ir,X%val,X_val_rich,done_richard,use_guess)
                CHCKERR('')
                call lvl_ud(-1)
            else                                                                ! if not, Richardson is done
                done_richard = .true.
            end if
            
            ! destroy grid if not yet done with Richardson Extrapolation
            if (.not.done_richard) call destroy_grid(grid_X)
            
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
            
            ! output on screen
            plot_title = 'job '//trim(i2str(alpha_job_nr))//' - Eigenvalues &
                &as function of nr. of normal points [log]'
            call print_GP_2D(plot_title,'Eigenvalues_'//&
                &trim(i2str(alpha_job_nr))//'_richardson.dat',&
                &log10(abs(realpart(X_val_rich(1:ir,1,:)))),&
                &x=x_axis(1:ir,:))
            
            ! same output in file as well
            call draw_GP(plot_title,'Eigenvalues_'//trim(i2str(alpha_job_nr))//&
                &'_richardson.dat',n_sol_requested,.true.,.false.)
            
            call lvl_ud(-1)
            call writo('Done plotting Eigenvalues')
        end if
            
        ! process output
        call writo('Processing output')
        call lvl_ud(1)
        
        ierr = process_output(grid_eq,eq,grid_X,X,n_sol_found,alpha)
        CHCKERR('')
        
        call lvl_ud(-1)
        call writo('Done processing output')
        
        ! deallocate Richardson loop variables
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            deallocate(X_val_rich)
            deallocate(x_axis)
        end if
        
        ! destroy grid
        call destroy_grid(grid_X)
        
        ! deallocate remaining equilibrium quantities
        call writo('Deallocate remaining quantities')
        call dealloc_eq_final(eq)
        call dealloc_metric_final(met)
        call dealloc_X_final(X)
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
            
            if (ir.eq.1) then
                n_r_X = min_n_r_X
            else
                n_r_X = 2 * n_r_X - 1
            end if
            call writo(trim(i2str(n_r_X))//' normal points for this level')
        end subroutine calc_n_r_X
        
        ! calculates  the  coefficients of  the  Eigenvalues  in the  Richardson
        ! extrapolati                                                         on
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
            
            ! tests
            if (size(X_val_rich,1).ne.size(X_val_rich,2) .or. &
                &size(X_val_rich,3).ne.size(X_val) .or. &
                &ir.gt.size(X_val_rich,1)) then
                ierr = 1
                err_msg = 'X_val_rich has to have correct dimensions'
                CHCKERR(err_msg)
            end if
            
            ! allocate correction
            allocate(corr(size(X_val))); corr = 1.0E15
            
            ! do calculations if ir > 1
            X_val_rich(ir,1,:) = X_val
            do ir2 = 2,ir
                corr = 1./(2**(2*ir2)-1.) * &
                    &(X_val_rich(ir,ir2-1,:) - X_val_rich(ir-1,ir2-1,:))
                X_val_rich(ir,ir2,:) = X_val_rich(ir,ir2-1,:) + corr
            end do
            
            if (ir.gt.1) then                                                   ! only do this if in Richardson level higher than 1
                ! get maximum and location of maximum for relative correction
                max_corr = maxval(abs(corr/X_val_rich(ir,ir,:)))
                loc_max_corr = maxloc(abs(corr/X_val_rich(ir,ir,:)))
                call writo('maximum relative error: '//trim(r2strt(max_corr))//&
                    &' for Eigenvalue '//trim(i2str(loc_max_corr(1))))
                
                ! check whether tolerance has been  reached for last value ir2 =
                ! ir
                if (maxval(abs(corr/X_val_rich(ir,ir,:))) .lt. tol_r) then
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
    end function run_for_alpha
    
    ! sets  up the equilibrium, metric  and some perturbation quantities  on the
    ! equilibrium grid
    integer function calculate_vars_in_eq_grid(grid_eq,eq,met,X,alpha) &
        &result(ierr)
        use grid_vars, only: grid_type
        use eq_vars, only: create_eq, eq_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type, create_X
        use eq_ops, only: calc_eq
        use X_ops, only: prepare_X
        
        character(*), parameter :: rout_name = 'calculate_vars_in_eq_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(metric_type), intent(inout) :: met                                 ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        real(dp), intent(in) :: alpha                                           ! field line label alpha
        
        ! initialize ierr
        ierr = 0
        
        ! 1. create equilibrium
        ierr = create_eq(eq,grid_eq)
        CHCKERR('')
        
        ! 2. create perturbation
        call create_X(grid_eq,X)
        
        ! 3. Calculate the equilibrium quantities for current alpha
        ierr = calc_eq(grid_eq,eq,met,X,alpha)
        CHCKERR('')
        
        ! 4. normalize the quantities
        call normalize(eq,met,X)
        
        ! 5. prepare the matrix elements
        ierr = prepare_X(eq,grid_eq,met,X)
        CHCKERR('')
    contains
        ! normalize quantities
        subroutine normalize(eq,met,X)
            use num_vars, only: use_normalization
            use eq_vars, only: eq_type
            use metric_vars, only: metric_type
            use X_vars, only: X_type
            use eq_ops, only: normalize_eq_vars
            use metric_ops, only: normalize_metric_vars
            use X_ops, only: normalize_X_vars
            
            ! input / output
            type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
            type(metric_type), intent(inout) :: met                                 ! metric variables
            type(X_type), intent(inout) :: X                                        ! perturbation variables
            
            if (use_normalization) then
                call writo('Start normalizing the equilibrium, metric and &
                    &perturbation quantities')
                call lvl_ud(1)
                call normalize_eq_vars(eq)
                call normalize_metric_vars(met)
                call normalize_X_vars(X)
                call lvl_ud(-1)
                call writo('Normalization done')
            end if
        end subroutine normalize
    end function calculate_vars_in_eq_grid
    
    ! Processes the output of the simulations for vizualization, analysis, etc.
    integer function process_output(grid_eq,eq,grid_X,X,n_sol_found,alpha) &
        &result(ierr)
        use num_vars, only: n_theta_plot, n_zeta_plot, no_plots, no_messages
        use grid_vars, only: create_grid, grid_type
        use eq_vars, only: create_eq, eq_type
        use metric_vars, only: metric_type
        use X_vars, only: create_X, X_type
        use X_ops, only: prepare_X
        
        character(*), parameter :: rout_name = 'process_output'
        
        ! input / output
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibirum variables
        type(grid_type) :: grid_X                                               ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: n_sol_found                                      ! how many solutions found and saved
        real(dp), intent(in) :: alpha                                           ! field line label alpha
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        type(grid_type) :: grid_eq_ext                                          ! extended equilibrium grid, covering whole geometry
        type(eq_type) :: eq_ext                                                 ! extended equilibrium variables
        type(metric_type) :: met_ext                                            ! extended metric variables
        type(X_type) :: X_ext                                                   ! extended perturbation variables
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_messages_loc                                              ! local copy of no_messages
        
        ! initialize ierr
        ierr = 0
        
        ! getting range of EV to plot
        call find_plot_ranges(X%val,n_sol_found,min_id,max_id,&
            &last_unstable_id)
        
        ! user output
        call writo('plotting the Eigenvalues')
        
        call lvl_ud(1)
        call plot_X_vals(X,n_sol_found,last_unstable_id)
        call lvl_ud(-1)
        
        ! user output
        call writo('creating plot grid')
        call lvl_ud(1)
        
        ! creating  equilibrium  grid  for  the  output  that covers  the  whole
        ! geometry angularly
        ierr = create_grid(grid_eq_ext,[n_theta_plot,n_zeta_plot,grid_eq%n(3)],&
            &[grid_eq%i_min,grid_eq%i_max])
        CHCKERR('')
        grid_eq_ext%grp_r_F = grid_eq%grp_r_F
        do id = 1,n_theta_plot
            grid_eq_ext%theta_F(id,:,:) = &
                &pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)                      ! starting from pi gives nicer plots
        end do
        do id = 1,n_zeta_plot
            grid_eq_ext%zeta_F(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
        end do
        
        ! user output
        call writo('calculating variables on plot grid')
        
        ! back up no_plots and no_messages and set to .true.
        no_plots_loc = no_plots
        no_messages_loc = no_messages
        no_plots = .true.
        no_messages = .true.
        ! create  equilibrium,  metric   and   some  perturbation  variables  on
        ! equilibrium grid
        ierr = calculate_vars_in_eq_grid(grid_eq_ext,eq_ext,met_ext,X_ext,alpha)
        CHCKERR('')
        ! reset no_plots and no_messages
        no_plots = no_plots_loc
        no_messages = no_messages_loc
        ! copy the Eigenvectors and -values into the extended X
        X_ext%val = X%val
        X_ext%vec = X%vec
        
        !! test whether flux quantities are the same
        !call print_GP_2D('rot_t_FD','',eq%rot_t_FD(:,0))
        !call print_GP_2D('OLD - NEW rot_t_FD','',eq%rot_t_FD(:,0)-eq_ext%rot_t_FD(:,0))
        !call print_GP_2D('q_saf_FD','',eq%q_saf_FD(:,0))
        !call print_GP_2D('OLD - NEW q_saf_FD','',eq%q_saf_FD(:,0)-eq_ext%q_saf_FD(:,0))
        !call print_GP_2D('pres_FD','',eq%pres_FD(:,0))
        !call print_GP_2D('OLD - NEW pres_FD','',eq%pres_FD(:,0)-eq_ext%pres_FD(:,0))
        !call print_GP_2D('flux_p_FD','',eq%flux_p_FD(:,0))
        !call print_GP_2D('OLD - NEW flux_p_FD','',eq%flux_p_FD(:,0)-eq_ext%flux_p_FD(:,0))
        !call print_GP_2D('flux_t_FD','',eq%flux_t_FD(:,0))
        !call print_GP_2D('OLD - NEW flux_t_FD','',eq%flux_t_FD(:,0)-eq_ext%flux_t_FD(:,0))
        
        ! user output
        call writo('plotting the Eigenvectors')
        
        call lvl_ud(1)
        ierr = plot_X_vecs(grid_eq_ext,eq_ext,grid_X,X_ext,n_sol_found,&
            &last_unstable_id,min_id,max_id)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call writo('Decomposing the energy into its terms')
        
        call lvl_ud(1)
        ierr = decompose_energy(grid_eq,eq,grid_X,X,n_sol_found)
        CHCKERR('')
        call lvl_ud(-1)
    contains
        ! finds the plot ranges min_id and max_id
        subroutine find_plot_ranges(X_val,n_sol_found,min_id,max_id,&
            &last_unstable_id)
            use num_vars, only: n_sol_plotted
            
            ! input / output
            complex(dp), intent(in) :: X_val(:)                                 ! Eigenvalues
            integer, intent(in) :: n_sol_found                                  ! how many solutions found and saved
            integer, intent(inout) :: min_id(3), max_id(3)                      ! min. and max. index of range 1, 2 and 3
            integer, intent(inout) :: last_unstable_id                          ! index of last unstable EV
            
            ! local variables
            integer :: id                                                       ! counter
            
            ! find last unstable index (if ends with 0, no unstable EV)
            last_unstable_id = 0
            do id = 1,n_sol_found
                if (realpart(X_val(id)).lt.0._dp) last_unstable_id = id
            end do
            ! set up min. and max. of range 1
            if (last_unstable_id.gt.0) then                                     ! there is an unstable range
                min_id(1) = 1                                                   ! start from most unstable EV
                max_id(1) = min(n_sol_plotted(1),last_unstable_id)              ! end with n_sol_plotted(1) first unstable values if available
            else                                                                ! no unstable range
                min_id(1) = 1                                                   ! no unstable values to plot
                max_id(1) = 0                                                   ! no unstable values to plot
            end if
            ! set up min. and max. of range 2
            if (last_unstable_id.gt.0) then                                     ! there is an unstable range
                min_id(2) = last_unstable_id - n_sol_plotted(2) + 1             ! start from n_sol_plotted(2) last unstable values
                max_id(2) = &
                    &min(last_unstable_id + n_sol_plotted(3),n_sol_found)       ! end with n_sol_plotted(3) first stable values if available
            else                                                                ! no unstable range
                min_id(2) = 1                                                   ! start from first EV (stable)
                max_id(2) = min(n_sol_plotted(3),n_sol_found)                   ! end with n_sol_plotted(3) first stable values if available
            end if
            ! set up min. and max. of range 3
            min_id(3) = n_sol_found - n_sol_plotted(4) + 1                      ! start from n_sol_plotted(4) last stable values
            max_id(3) = n_sol_found                                             ! end with most stable EV
            ! merge ranges 2 and 3 if overlap
            if (min_id(3).le.max_id(2)) then
                max_id(2) = max_id(3)                                           ! range 3 merged with range 2
                min_id(3) = 1                                                   ! range 3 not used any more
                max_id(3) = 0                                                   ! range 3 not used any more
            end if
            ! merge ranges 1 and 2 if overlap
            if (min_id(2).le.max_id(1)) then
                max_id(1) = max_id(2)                                           ! range 2 merged with range 1
                min_id(2) = 1                                                   ! range 2 not used any more
                max_id(2) = 0                                                   ! range 2 not used any more
            end if
        end subroutine find_plot_ranges
        
        ! plots Eigenvalues
        subroutine plot_X_vals(X,n_sol_found,last_unstable_id)
            use num_vars, only: grp_rank, alpha_job_nr
            use X_vars, only: X_type
            
            ! input / output
            type(X_type), intent(in) :: X                                       ! perturbation variables
            integer, intent(in) :: n_sol_found                                  ! how many solutions found and saved
            integer, intent(in) :: last_unstable_id                             ! index of last unstable EV
            
            ! local variables
            character(len=max_str_ln) :: plot_title                             ! title for plots
            integer :: id                                                       ! counter
            
            if (grp_rank.eq.0) then
                ! Last Eigenvalues
                ! output on screen
                plot_title = 'job '//trim(i2str(alpha_job_nr))//' - &
                    &final Eigenvalues omega^2 [log]'
                call print_GP_2D(plot_title,'Eigenvalues_'//&
                    &trim(i2str(alpha_job_nr))//'.dat',&
                    &log10(abs(realpart(X%val(1:n_sol_found)))))
                ! same output in file as well
                call draw_GP(plot_title,'Eigenvalues_'//&
                    &trim(i2str(alpha_job_nr))//'.dat',1,.true.,.false.)
                
                ! Last Eigenvalues: unstable range
                if (last_unstable_id.gt.0) then
                    ! output on screen
                    plot_title = 'job '//trim(i2str(alpha_job_nr))//' - &
                        &final unstable Eigenvalues omega^2'
                    call print_GP_2D(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_unstable.dat',&
                        &realpart(X%val(1:last_unstable_id)),&
                        &x=[(id*1._dp,id=1,last_unstable_id)])
                    ! same output in file as well
                    call draw_GP(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_unstable.dat',1,&
                        &.true.,.false.)
                end if
                
                ! Last Eigenvalues: stable range
                if (last_unstable_id.lt.n_sol_found) then
                    ! output on screen
                    plot_title = 'job '//trim(i2str(alpha_job_nr))//' - &
                        &final stable Eigenvalues omega^2'
                    call print_GP_2D(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_stable.dat',&
                        &realpart(X%val(last_unstable_id+1:n_sol_found)),&
                        &x=[(id*1._dp,id=last_unstable_id+1,n_sol_found)])
                    ! same output in file as well
                    call draw_GP(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_stable.dat',1,&
                        &.true.,.false.)
                end if
            end if
        end subroutine plot_X_vals
        
        ! Plots  Eigenvectors  using  the  angular  part  of  the  the  provided
        ! equilibrium  grid and  the normal  part of  the provided  perturbation
        ! grid.
        ! Note: These don't necessarily have to  coincide with the grids used in
        ! the calculations.
        integer function plot_X_vecs(grid_eq,eq,grid_X,X,n_sol_found,&
            &last_unstable_id,min_id,max_id) result(ierr)
            use num_vars, only: use_pol_flux_F, alpha_job_nr
            use grid_vars, only: grid_type
            use eq_vars, only: eq_type
            use X_vars, only: X_type
            use X_ops, only: plot_X_vec
            use sol_ops, only: calc_real_XUQ
            
            character(*), parameter :: rout_name = 'plot_X_vecs'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(eq_type), intent(in) :: eq                                     ! equilibrium variables
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_type), intent(in) :: X                                       ! perturbation variables
            integer, intent(in) :: n_sol_found                                  ! how many solutions found and saved
            integer, intent(in) :: last_unstable_id                             ! index of last unstable EV
            integer, intent(in) :: min_id(3), max_id(3)                         ! min. and max. index of range 1, 2 and 3
            real(dp), allocatable :: X_F(:,:,:,:)                               ! the normal component of the plasma perturbation
            
            ! local variables
            integer :: id, jd                                                   ! counters
            
            ! initialize ierr
            ierr = 0
            
            ! for all the solutions that are to be saved
            ! three ranges
            do jd = 1,3
                if (min_id(jd).le.max_id(jd)) call &
                    &writo('Plotting results for modes '//&
                    &trim(i2str(min_id(jd)))//'..'//trim(i2str(max_id(jd)))&
                    &//' of range '//trim(i2str(jd)))
                call lvl_ud(1)
                
                ! indices in each range
                do id = min_id(jd),max_id(jd)
                    ! user output
                    call writo('Mode '//trim(i2str(id))//'/'//&
                        &trim(i2str(n_sol_found))//', with eigenvalue '&
                        &//trim(r2strt(realpart(X%val(id))))//' + '//&
                        &trim(r2strt(imagpart(X%val(id))))//' i')
                    
                    call lvl_ud(1)
                    
                    ! plot information about harmonics
                    ierr = plot_harmonics(grid_X,X,id)
                    CHCKERR('')
                    
                    ! calculate the normal component of the perturbation
                    !ierr = calc_real_XUQ_arr(grid_eq,eq,grid_X,X,X_id,style,&
                        !&time,X_F,met,deriv)
                    CHCKERR('')
                    
                    call lvl_ud(-1)
                end do
                
                call lvl_ud(-1)
            end do
        end function plot_X_vecs
        
        ! plots the harmonics and their maximum in 2D
        integer function plot_harmonics(grid_X,X,X_id) result(ierr)
            use MPI_ops, only: wait_MPI, get_ghost_arr, get_ser_var
            use output_ops, only: merge_GP
            use num_vars, only: grp_n_procs, grp_rank, alpha_job_nr
            use grid_vars, only: grid_type
            use X_vars, only: X_type
            
            character(*), parameter :: rout_name = 'plot_harmonics'
            
            ! input / output
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_type), intent(in) :: X                                       ! perturbation variables/
            integer, intent(in) :: X_id                                         ! nr. of Eigenvalue (for output name)
            
            ! local variables
            integer :: id, kd                                                   ! counters
            character(len=max_str_ln) :: file_name                              ! name of file of plots of this proc.
            character(len=max_str_ln), allocatable :: file_names(:)             ! names of file of plots of different procs.
            character(len=max_str_ln) :: plot_title                             ! title for plots
            real(dp), allocatable :: x_plot(:,:)                                ! x values of plot
            complex(dp), allocatable :: X_vec_extended(:,:)                     ! MPI Eigenvector extended with assymetric ghost region
            real(dp), allocatable :: X_vec_max(:)                               ! maximum position index of X_vec of rank
            real(dp), allocatable :: ser_X_vec_max(:)                           ! maximum position index of X_vec of whole group
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Started plot of the harmonics')
            call lvl_ud(1)
            
            ! set up extended X_vec with ghost values
            allocate(X_vec_extended(X%n_mod,size(grid_X%grp_r_F)))
            X_vec_extended(:,1:size(X%vec,2)) = X%vec(:,:,X_id)
            ierr = get_ghost_arr(X_vec_extended,1)
            CHCKERR('')
            
            ! set up x_plot
            allocate(x_plot(size(grid_X%grp_r_F),X%n_mod))
            do kd = 1,X%n_mod
                x_plot(:,kd) = grid_X%grp_r_F
            end do
            
            ! absolute amplitude
            ! set up file name of this rank and plot title
            file_name = 'Eigenvector_'//trim(i2str(alpha_job_nr))//&
                &'_abs.dat'
            plot_title = 'job '//trim(i2str(alpha_job_nr))//' - EV '//&
                &trim(i2str(X_id))//' - absolute value'
            
            ! print amplitude of harmonics of eigenvector for each rank
            call print_GP_2D(trim(plot_title),trim(file_name)//'_'//&
                &trim(i2str(grp_rank)),abs(transpose(X_vec_extended(:,:))),&
                &x=x_plot,draw=.false.)
            
            ! wait for all processes
            ierr = wait_MPI()
            CHCKERR('')
            
            ! plot by group master
            if (grp_rank.eq.0) then
                ! set up file names in array
                allocate(file_names(grp_n_procs))
                do kd = 1,grp_n_procs
                    file_names(kd) = trim(file_name)//'_'//trim(i2str(kd-1))
                end do
                
                ! merge files
                call merge_GP(file_names,file_name,delete=.true.)
                
                ! draw plot
                call draw_GP(trim(plot_title),file_name,X%n_mod,.true.,.false.)
                
                ! deallocate
                deallocate(file_names)
            end if
            
            ! perturbation at midplane theta = zeta = 0
            ! set up file name of this rank and plot title
            file_name = 'Eigenvector_'//trim(i2str(alpha_job_nr))//&
                &'_midplane.dat'
            plot_title = 'job '//trim(i2str(alpha_job_nr))//' - EV '//&
                &trim(i2str(X_id))//' - midplane'
            
            ! print amplitude of harmonics of eigenvector for each rank
            call print_GP_2D(trim(plot_title),trim(file_name)//'_'//&
                &trim(i2str(grp_rank)),&
                &realpart(transpose(X_vec_extended(:,:))),x=x_plot,draw=.false.)
            
            ! wait for all processes
            ierr = wait_MPI()
            CHCKERR('')
            
            ! plot by group master
            if (grp_rank.eq.0) then
                ! set up file names in array
                allocate(file_names(grp_n_procs))
                do kd = 1,grp_n_procs
                    file_names(kd) = trim(file_name)//'_'//trim(i2str(kd-1))
                end do
                
                ! merge files
                call merge_GP(file_names,file_name,delete=.true.)
                
                ! draw plot
                call draw_GP(trim(plot_title),file_name,X%n_mod,.true.,.false.)
                
                ! deallocate
                deallocate(file_names)
            end if
            
            ! maximum of each mode
            allocate(X_vec_max(X%n_mod))
            X_vec_max = 0.0_dp
            do kd = 1,X%n_mod
                X_vec_max(kd) = grid_X%grp_r_F(maxloc(abs(X%vec(kd,:,X_id)),1))
            end do
            
            ! gather all parllel X_vec_max arrays in one serial array
            ierr = get_ser_var(X_vec_max,ser_X_vec_max)
            CHCKERR('')
            
            ! find the maximum of the different ranks and put it in X_vec_max of
            ! group master
            if (grp_rank.eq.0) then
                do kd = 1,X%n_mod
                    X_vec_max(kd) = maxval([(ser_X_vec_max(kd+id*X%n_mod),&
                        &id=0,grp_n_procs-1)])
                end do
                
                ! set up file name of this rank and plot title
                file_name = 'Eigenvector_'//trim(i2str(alpha_job_nr))//'_max.dat'
                plot_title = 'job '//trim(i2str(alpha_job_nr))//' - EV '//&
                    &trim(i2str(X_id))//' - maximum of modes'
                
                ! plot the maximum
                call print_GP_2D(trim(plot_title),trim(file_name),&
                    &[(kd*1._dp,kd=1,X%n_mod)],x=X_vec_max,draw=.false.)
                
                ! draw plot
                call draw_GP(trim(plot_title),file_name,1,.true.,.false.)
            end if
            
            ! deallocate
            deallocate(x_plot,X_vec_extended)
            
            ! user output
            call lvl_ud(-1)
            call writo('Finished plot of the harmonics')
        end function plot_harmonics
        
        ! Decomposes the  energy to see the fraction of  the individual terms in
        ! the total energy on a grid.
        ! The terms are:
        !   - normal line bending term = Qn^2/mu_0 1/|nabla psi|^2
        !   - geodesic line bending term = Qg^2/mu_0 |nabla psi|^2/B^2
        !   - normal ballooning term = -2 X X* p' kappa_n
        !   - geodesic ballooning term = -2 X U* p' kappa_g
        !   - normal kink term = -sigma X* Qg
        !   - geodesic kink term = sigma U* Qn
        integer function decompose_energy(grid_eq,eq,grid_X,X,n_sol_found) &
            &result(ierr)
            use sol_ops, only: calc_real_XUQ, calc_real_U, calc_real_Qn, &
                &calc_real_Qg
            use grid_vars, only: grid_type
            use eq_vars, only: eq_type
            use X_vars, only: X_type
            use grid_ops, only: calc_eqd_grid
            
            character(*), parameter :: rout_name = 'decompose_energy'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(eq_type), intent(in) :: eq                                     ! equilibirum variables
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_type), intent(in) :: X                                       ! perturbation variables
            integer, intent(in) :: n_sol_found                                  ! how many solutions found and saved
            
            ! local variables
            integer :: jd                                                       ! counter
            integer :: n_t = 10                                                 ! nr. of time steps
            real(dp), allocatable :: time(:)                                    ! time at which to calculate decomposition
            real(dp), allocatable :: X_F(:,:,:,:)                               ! normal comp. of perturbation
            real(dp), allocatable :: D3X_F(:,:,:,:)                             ! parallel deriv. of X
            real(dp), allocatable :: U_F(:,:,:,:)                               ! geodesic comp. of perturbation
            real(dp), allocatable :: D3U_F(:,:,:,:)                             ! parallel deriv. of U
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Merging perturbation and equilibrium grid into custom &
                &grid')
            call lvl_ud(1)
            
            ! intialize variables
            !allocate(time(n_t))
            !allocate(X_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            !allocate(D3X_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            !allocate(U_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            !allocate(D3U_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            
            ! user output
            call writo('Calculating X, U, Qn and Qg')
            call lvl_ud(1)
            
            ! set up time
            ierr = calc_eqd_grid(time,0.0_dp,1.0_dp,excl_last=.true.)
            CHCKERR('')
            
            ! loop over all Eigenvalues
            do jd = 1,n_sol_found
                ! calculate X, U, Qn and Qg in real space
                !ierr = calc_real_XUQ(eq,grid_eq,X,jd,time,X_F)
                CHCKERR('')
            end do
            
            call lvl_ud(-1)
        end function decompose_energy
    end function process_output
end module driver_rich
