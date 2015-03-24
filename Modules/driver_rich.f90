!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
#include <PB3D_macros.h>
    use num_vars, only: max_it_r, dp, pi, max_str_ln
    use str_ops, only: i2str, r2str, r2strt
    use messages, only: writo, print_ar_2, print_ar_1, lvl_ud
    use output_ops, only: print_GP_2D, print_GP_3D, print_HDF5
    use grid_vars, only: grid_type
    
    implicit none
    private
    public run_rich_driver, calc_eq

    ! global variables
    real(dp), allocatable :: alpha(:)
    type(grid_type) :: grid_eq                                                  ! equilibrium grid
    type(grid_type) :: grid_X                                                   ! perturbation grid
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_rich_driver() result(ierr)
        use num_vars, only: min_alpha, max_alpha, n_alpha, glb_rank, grp_nr, &
            &max_alpha, alpha_job_nr, use_pol_flux_X, eq_style
        use MPI_ops, only: split_MPI, merge_MPI, get_next_job
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X
        use VMEC, only: dealloc_VMEC_final
        use HELENA, only: dealloc_HEL_final
        use grid_ops, only: calc_eqd_grid
        use grid_vars, only: create_grid, &
            &n_r_eq, n_par_X
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: eq_limits(2)                                                 ! min. and max. index of eq. grid of this process
        
        ! initialize ierr
        ierr = 0
        
        ! set up flux_name
        if (use_pol_flux_X) then
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
        if (use_pol_flux_X) then
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
        ierr = calc_eqd_grid(alpha,n_alpha,min_alpha*pi,max_alpha*pi,&
            &excl_last=.true.)
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
                ierr = run_for_alpha(alpha_job_nr,alpha(alpha_job_nr))
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
    integer function run_for_alpha(job_nr,alpha) result(ierr)
        use num_vars, only: n_sol_requested, max_it_r, grp_rank, no_guess, &
            &alpha_job_nr, use_normalization, use_pol_flux_X
        use output_ops, only: draw_GP
        use eq_vars, only: create_eq, dealloc_eq_final, eq_type
        use eq_ops, only: normalize_eq_vars
        use metric_ops, only: normalize_metric_vars
        use X_vars, only: create_X, dealloc_X_final, X_type
        use X_ops, only: prepare_X, solve_EV_system, plot_X_vec, &
            &normalize_X_vars
        use metric_vars, only: dealloc_metric_final, metric_type
        use MPI_ops, only: divide_X_grid
        use grid_vars, only: create_grid, destroy_grid
        
        character(*), parameter :: rout_name = 'run_for_alpha'
        
        ! input / output
        integer, intent(in) :: job_nr                                           ! job nr.
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
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle in flux coordinates
        
        ! initialize ierr
        ierr = 0
        
        ! create equilibrium
        call writo('Create equilibrium...')
        ierr = create_eq(eq,grid_eq)
        CHCKERR('')
        
        ! create perturbation
        call create_X(grid_eq,X)
        
        ! Calculate the equilibrium quantities for current alpha
        ierr = calc_eq(eq,met,X,alpha)
        CHCKERR('')
        
        ! normalize the quantities
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
        
        ! prepare the  matrix elements by calculating  the integrated magnitudes
        ! KV_int and  PV_int for each of  the n_r equilibrium normal  points and
        ! for the modes (k,m)
        call writo('Preparing perturbation variables')
        call lvl_ud(1)
        ierr = prepare_X(eq,grid_eq,met,X)
        CHCKERR('')
        call lvl_ud(-1)
        call writo('Perturbation variables prepared')
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_X) then
            ang_par_F => grid_eq%theta_F
        else
            ang_par_F => grid_eq%zeta_F
        end if
        
        write(*,*) 'INVESTIGATE HOW IMPROVING THE DERIVATIVES CAN &
            &HELP YOU GET RID OF UNSTABLE SIDE OF SPECTRUM !!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'DOING THIS IN THIS ROUTINE ALREADY HELPED A LOT!!!!!'
        
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
            
            ! divide perturbation grid under group processes
            ierr = divide_X_grid(n_r_X,X_limits,grp_r_X)
            CHCKERR('')
            
            ! create perturbation grid with division limits
            call writo('Creating perturbation grid')
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
            ierr = solve_EV_system(eq,grid_X,X,use_guess,n_sol_found)
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
            
            ! getting range of EV to plot
            call find_plot_ranges(X%val,n_sol_found,min_id,max_id,&
                &last_unstable_id)
            
            ! output the evolution of the  Eigenvalues with the number of normal
            ! points in the perturbation grid  and the largest Eigenfunction for
            ! the last Richardson loop
            if (done_richard .or. ir.eq.max_it_r+1) then
                call writo('plotting the Eigenvalues')
                
                call lvl_ud(1)
                call plot_X_vals(x_axis)
                call lvl_ud(-1)
                
                call writo('plotting the Eigenvectors')
                
                call lvl_ud(1)
                call plot_X_vecs()
                call lvl_ud(-1)
            end if
            
            ! destroy grid
            call destroy_grid(grid_X)
            
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
        
        ! deallocate Richardson loop variables
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            deallocate(X_val_rich)
            deallocate(x_axis)
        end if
        
        ! deallocate remaining equilibrium quantities
        call writo('Deallocate remaining quantities')
        call dealloc_eq_final(eq)
        call dealloc_metric_final(met)
        call dealloc_X_final(X)
    contains
        ! plots the harmonics and their maximum in 2D
        subroutine plot_harmonics(X_vec,X_id)
            use MPI_ops, only: wait_MPI, get_ghost_arr, get_ser_var
            use output_ops, only: merge_GP
            use num_vars, only: grp_n_procs
            
            ! input / output
            complex(dp), intent(in) :: X_vec(:,:)                               ! MPI Eigenvector
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
            
            ! user output
            call writo('Started plot of the harmonics')
            call lvl_ud(1)
            
            ! set up extended X_vec with ghost values
            allocate(X_vec_extended(X%n_mod,size(grid_X%grp_r_F)))
            X_vec_extended(:,1:size(X_vec,2)) = X_vec
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
                X_vec_max(kd) = grid_X%grp_r_F(maxloc(abs(X_vec(kd,:)),1))
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
        end subroutine plot_harmonics
        
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
        subroutine plot_X_vals(x_axis)
            ! input / output
            real(dp), intent(in) :: x_axis(:,:)                                    ! x axis for plot of Eigenvalues with n_r_X
            
            ! local variables
            character(len=max_str_ln) :: plot_title                             ! title for plots
            integer :: id                                                       ! counter
            
            if (grp_rank.eq.0) then
                ! Eigenvalues as function of nr. of normal points
                if (max_it_r.gt.1) then
                    call writo('plotting Eigenvalues as function of nr. of &
                        &normal points')
                    
                    ! output on screen
                    plot_title = 'job '//trim(i2str(alpha_job_nr))//' - &
                        &Eigenvalues as function of nr. of normal points &
                        &[log]'
                    call print_GP_2D(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_richardson.dat',&
                        &log10(abs(realpart(X_val_rich(1:ir,1,:)))),&
                        &x=x_axis(1:ir,:))
                    ! same output in file as well
                    call draw_GP(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_richardson.dat',&
                        &n_sol_requested,.true.,.false.)
                end if
                call writo('plotting final Eigenvalues')
                
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
        
        ! plots Eigenvectors
        subroutine plot_X_vecs()
            ! local variables
            integer :: id, jd                                                   ! counters
            
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
                    call plot_harmonics(X%vec(:,:,id),id)
                    
                    ! plot the vector
                    ierr = plot_X_vec(grid_eq,eq,grid_X,X,id,job_nr,&
                        &[ang_par_F(1,1,1),ang_par_F(grid_eq%n(1),1,1)])
                    CHCKERR('')
                    
                    call lvl_ud(-1)
                end do
                
                call lvl_ud(-1)
            end do
        end subroutine plot_X_vecs
    end function run_for_alpha
    
    ! calculate  the equilibrium  quantities on  a grid  determined by  straight
    ! field lines. This  grid has the dimensions  (n_par,grp_n_r_eq) where n_par
    ! is  the  number  of  points  taken along  the  magnetic  field  lines  and
    ! grp_n_r_eq .le.  n_r_eq is the  normal extent  in the equilibrium  grid of
    ! this rank. It is determined so  that the perturbation quantities that will
    ! be needed in this rank are fully covered, so no communication is necessary
    integer function calc_eq(eq,met,X,alpha) result(ierr)
        use eq_vars, only: dealloc_eq, eq_type
        use eq_ops, only: calc_flux_q, prepare_RZL, calc_RZL, &
            &adapt_HEL_to_eq
        use metric_ops, only: calc_g_C, calc_g_C, calc_T_VC, calc_g_V, &
            &calc_T_VF, calc_inv_met, calc_g_F, calc_jac_C, calc_jac_V, &
            &calc_jac_F, calc_f_deriv, calc_jac_H, calc_T_HF, calc_h_H
        use metric_vars, only: create_metric, dealloc_metric, metric_type
        use utilities, only: derivs
        use num_vars, only: max_deriv, ltest, plot_grid, eq_style
        use grid_ops, only: calc_ang_grid, plot_grid_real
        use HELENA, only: dealloc_HEL
        use X_vars, only: X_type
        
        use utilities, only: calc_det, calc_inv, calc_mult, c
        !use metric_ops, only: plot_info
        !use utilities, only: plot_info_2 => plot_info
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium variables
        type(metric_type), intent(inout) :: met                                 ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        real(dp), intent(in) :: alpha                                           ! field line coordinate of current equilibrium
        
        ! local variables
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: pmone                                                        ! plus or minus one
        
        !! TEMPORARILY, TO PLOT MORE INFO IN ROUTINES
        !plot_info = .true.
        !plot_info_2 = .true.
        
        !real(dp) :: A(1,1,1,16)
        !real(dp) :: B(1,1,1,10)
        !real(dp) :: det(1,1,1)
        !real(dp) :: invA(1,1,1,16)
        !real(dp) :: invB(1,1,1,10)
        !real(dp) :: AinvA(1,1,1,16)
        !real(dp) :: BinvB(1,1,1,10)
        !real(dp) :: AinvB(1,1,1,16)
        !real(dp) :: BinvA(1,1,1,16)
        
        !A = reshape([1,2,3,1,2,1,1,4,3,1,3,2,1,4,2,2],shape(A))
        !B = reshape([1,2,3,1,1,1,4,3,2,2],shape(B))
        !write(*,*) 'A = ', A
        !write(*,*) 'B = ', B
        !ierr = calc_det(det,A,4)
        !CHCKERR('')
        !write(*,*) 'detA = ', det
        !ierr = calc_det(det,B,4)
        !CHCKERR('')
        !write(*,*) 'detB = ', det
        !ierr = calc_inv(invA,A,4)
        !CHCKERR('')
        !write(*,*) 'invA = ', invA
        !ierr = calc_inv(invB,B,4)
        !CHCKERR('')
        !write(*,*) 'invB = ', invB
        !ierr = calc_mult(A,invA,AinvA,4)
        !CHCKERR('')
        !write(*,*) 'AinvA = ', AinvA
        !ierr = calc_mult(B,invB,BinvB,4)
        !CHCKERR('')
        !write(*,*) 'BinvB = ', BinvB
        !ierr = calc_mult(A,invB,AinvB,4)
        !CHCKERR('')
        !write(*,*) 'AinvB = ', AinvB
        !ierr = calc_mult(B,invA,BinvA,4)
        !CHCKERR('')
        !write(*,*) 'BinvA = ', BinvA
        
        ! initialize ierr
        ierr = 0
        
        call writo('Start setting up equilibrium quantities in '//&
            &trim(i2str(grid_eq%n(1)))//' discrete parallel points')
        call lvl_ud(1)
        
            call writo('Start initializing variables')
            call lvl_ud(1)
                
                ! initialize metric quantities
                call writo('Initialize metric quantities...')
                ierr = create_metric(grid_eq,met)
                CHCKERR('')
                
                ! calculate flux quantities
                call writo('Calculate flux quantities...')
                ierr = calc_flux_q(eq,grid_eq,X)
                CHCKERR('')
                
            call lvl_ud(-1)
            call writo('Variables initialized')
            
            call writo('Start determining the equilibrium grid')
            call lvl_ud(1)
            
                ! calculate  angular grid  points (theta,zeta)  that follow  the
                ! magnetic field  line determined  by alpha in  E(quilibrum) and
                ! F(lux) coords.
                ! Note: The  normal grid  is determined  by the  equilibrium, it
                ! can  either  use  the  toroidal  flux  or  the  poloidal  flux
                ! (use_pol_flux_eq).  This  does  not   have  to  coincide  with
                ! use_pol_flux_X used by PB3D.
                ierr = calc_ang_grid(grid_eq,eq,alpha)
                CHCKERR('')
                
                ! adapt the HELENA variables to the equilibrium grid
                if (eq_style.eq.2) then
                    ierr = adapt_HEL_to_eq(grid_eq)
                    CHCKERR('')
                end if
                
                ! plot grid if requested
                if (plot_grid) then
                    ierr = plot_grid_real(grid_eq)
                    CHCKERR('')
                else
                    call writo('Grid plot not requested')
                end if
            
            call lvl_ud(-1)
            call writo('Equilibrium grid determined')
            
            call writo('Calculating equilibrium quantities on equilibrium grid')
            call lvl_ud(1)
            
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! calculate the  cylindrical variables  R, Z  and lambda
                        ! and derivatives
                        call writo('Calculate R,Z,L...')
                        ierr = prepare_RZL(grid_eq)
                        CHCKERR('')
                        do id = 0,max_deriv+1
                            ierr = calc_RZL(grid_eq,eq,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate  the metrics  in the  cylindrical coordinate
                        ! system
                        call writo('Calculate g_C...')                          ! h_C is not necessary
                        do id = 0,max_deriv
                            ierr = calc_g_C(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the  jacobian in the  cylindrical coordinate
                        ! system
                        call writo('Calculate jac_C...')
                        do id = 0,max_deriv
                            ierr = calc_jac_C(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the  transformation matrix  C(ylindrical) ->
                        ! V(MEC)
                        call writo('Calculate T_VC...')
                        do id = 0,max_deriv
                            ierr = calc_T_VC(eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate  the metric  factors in the  VMEC coordinate
                        ! system
                        call writo('Calculate g_V...')
                        do id = 0,max_deriv
                            ierr = calc_g_V(met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the jacobian in the VMEC coordinate system
                        call writo('Calculate jac_V...')
                        do id = 0,max_deriv
                            ierr = calc_jac_V(met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the transformation matrix V(MEC) -> F(lux)
                        call writo('Calculate T_VF...')
                        do id = 0,max_deriv
                            ierr = calc_T_VF(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! set up plus minus one
                        pmone = -1                                              ! conversion VMEC LH -> RH coord. system
                    case (2)                                                    ! HELENA
                        ! calculate the jacobian in the HELENA coordinate system
                        call writo('Calculate jac_H...')
                        do id = 0,max_deriv
                            ierr = calc_jac_H(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the metric factors  in the HELENA coordinate
                        ! system
                        call writo('Calculate h_H...')
                        do id = 0,max_deriv
                            ierr = calc_h_H(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the inverse g_H of the metric factors h_H
                        call writo('Calculate g_H...')
                        do id = 0,max_deriv
                            ierr = calc_inv_met(met%g_E,met%h_E,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! calculate the transformation matrix H(ELENA) -> F(lux)
                        call writo('Calculate T_HF...')
                        do id = 0,max_deriv
                            ierr = calc_T_HF(grid_eq,eq,met,derivs(id))
                            CHCKERR('')
                        end do
                        
                        ! set up plus minus one
                        pmone = 1
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! calculate  the  inverse of  the transformation  matrix T_EF
                call writo('Calculate T_FE...')
                do id = 0,max_deriv
                    ierr = calc_inv_met(met%T_FE,met%T_EF,derivs(id))
                    CHCKERR('')
                    ierr = calc_inv_met(met%det_T_FE,met%det_T_EF,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the metric factors in the Flux coordinate system
                call writo('Calculate g_F...')
                do id = 0,max_deriv
                    ierr = calc_g_F(met,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the inverse h_F of the metric factors g_F
                call writo('Calculate h_F...')
                do id = 0,max_deriv
                    ierr = calc_inv_met(met%h_F,met%g_F,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate the jacobian in the Flux coordinate system
                call writo('Calculate jac_F...')
                do id = 0,max_deriv
                    ierr = calc_jac_F(met,derivs(id))
                    CHCKERR('')
                end do
                
                ! calculate   the  derivatives  in  Flux  coordinates  from  the
                ! derivatives in VMEC coordinates
                call writo('Calculate derivatives in Flux coordinates...')
                do id = 0,max_deriv
                    ! g_FD
                    ierr = calc_f_deriv(met%g_F,met%T_FE,met%g_FD,max_deriv,&
                        &derivs(id))
                    CHCKERR('')
                    
                    ! h_FD
                    ierr = calc_f_deriv(met%h_F,met%T_FE,met%h_FD,max_deriv,&
                        &derivs(id))
                    CHCKERR('')
                    
                    ! jac_FD
                    ierr = calc_f_deriv(met%jac_F,met%T_FE,met%jac_FD,&
                        &max_deriv,derivs(id))
                    CHCKERR('')
                    
                    ! pres_FD
                    ierr = calc_f_deriv(eq%pres_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%pres_FD(:,id),max_deriv,id)
                    CHCKERR('')
                        
                    ! flux_p_FD
                    ierr = calc_f_deriv(eq%flux_p_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%flux_p_FD(:,id),max_deriv,id)
                    CHCKERR('')
                        
                    ! flux_t_FD
                    ierr = calc_f_deriv(eq%flux_t_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%flux_t_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%flux_t_FD(:,id) = pmone * eq%flux_t_FD(:,id)             ! multiply by plus minus one
                    
                    ! q_saf_FD
                    ierr = calc_f_deriv(eq%q_saf_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%q_saf_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%q_saf_FD(:,id) = pmone * eq%q_saf_FD(:,id)               ! multiply by plus minus one
                    
                    ! rot_t_FD
                    ierr = calc_f_deriv(eq%rot_t_E,&
                        &met%T_FE(1,1,:,c([2,1],.false.),:,0,0),&
                        &eq%rot_t_FD(:,id),max_deriv,id)
                    CHCKERR('')
                    eq%rot_t_FD(:,id) = pmone * eq%rot_t_FD(:,id)               ! multiply by plus minus one
                end do
                
                !write(*,*) 'det_T_FE:', 0,0,0
                !call print_HDF5('X','X',met%det_T_FE(:,:,:,0,0,0))
                !read(*,*)
                !do id = 1,9
                    !write(*,*) 'h_FD:', id
                    !call print_HDF5('X','X',met%h_FD(:,:,:,id,0,0,0))
                    !read(*,*)
                !end do
                
                ! deallocate unused equilibrium quantities
                if (.not.ltest) then
                    call writo('Deallocate unused equilibrium and metric &
                        &quantities...')
                    ierr = dealloc_metric(met)
                    CHCKERR('')
                    ! general equilibrium
                    ierr = dealloc_eq(eq)
                    CHCKERR('')
                    ! specific equilibrium
                    ! choose which equilibrium style is being used:
                    !   1:  VMEC
                    !   2:  HELENA
                    select case (eq_style)
                        case (1)                                                ! VMEC
                            ! nothing
                        case (2)                                                ! HELENA
                            call dealloc_HEL
                        case default
                            err_msg = 'No equilibrium style associated with '//&
                                &trim(i2str(eq_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                end if
            
            call lvl_ud(-1)
            call writo('Equilibrium quantities calculated on equilibrium grid')
            
        call lvl_ud(-1)
        call writo('Done setting up equilibrium quantities')
    end function calc_eq
    
    ! calculates the number of normal  points for the perturbation n_r_X for the
    ! various Richardson iterations
    ! The aim is to  halve the step size, which is given by  dx(n) = 1/(n-1) or,
    ! inverting: n(dx) = 1 + 1/dx.
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
    
    ! calculates  the  coefficients  of   the  Eigenvalues  in   the  Richardson
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
        integer, intent(inout) :: ir                                            ! level of Richardson extrapolation (starting at 1)
        complex(dp), intent(in) :: X_val(:)                                     ! EV for this Richardson level
        complex(dp), intent(inout) :: X_val_rich(:,:,:)                         ! extrapolated coefficients
        logical, intent(inout) :: done_richard                                  ! if Richardson loop has converged sufficiently
        logical, intent(inout) :: use_guess_for_next_level                      ! if a guessed is used for next Richardson level
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: ir2                                                          ! counter
        complex(dp), allocatable :: corr(:)                                     ! correction
        real(dp) :: max_corr                                                    ! maximum of maximum of correction
        real(dp), save :: prev_max_corr                                         ! max_corr at previous Richardson level
        integer :: loc_max_corr(1)                                              ! location of maximum of correction
        
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
        
        if (ir.gt.1) then                                                       ! only do this if in Richardson level higher than 1
            ! get maximum and location of maximum for relative correction
            max_corr = maxval(abs(corr/X_val_rich(ir,ir,:)))
            loc_max_corr = maxloc(abs(corr/X_val_rich(ir,ir,:)))
            call writo('maximum relative error: '//trim(r2strt(max_corr))//&
                &' for Eigenvalue '//trim(i2str(loc_max_corr(1))))
            
            ! check whether tolerance has been reached for last value ir2 = ir
            if (maxval(abs(corr/X_val_rich(ir,ir,:))) .lt. tol_r) then
                done_richard = .true.
                call writo('tolerance '//trim(r2strt(tol_r))//&
                    &' reached after '//trim(i2str(ir))//' iterations')
            else
                call writo('tolerance '//trim(r2strt(tol_r))//' not yet &
                    &reached')
            end if
            
            ! determine use_guess_for_next_level
            use_guess_for_next_level = max_corr.lt.prev_max_corr                ! set guess for next level to true if max_corr is decreasing
            prev_max_corr = max_corr                                            ! save max_corr for next level
        else
            use_guess_for_next_level = .true.                                   ! for first Richardson level, set guess for next level to true
        end if
        
        ! check for convergence
        if (.not.done_richard) then
            if (ir.lt.max_it_r) then                                            ! not yet at maximum Richardson iteration
                ir = ir + 1
            else                                                                ! maximum nr. of Richardson iterations reached
                call writo('maximum number of Richardson iterations reached')
                done_richard = .true.
            end if
        end if
    end function calc_rich_ex
end module driver_rich
