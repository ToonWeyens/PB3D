!------------------------------------------------------------------------------!
!   Post-processing for the results of the main driver.                        !
!------------------------------------------------------------------------------!
module post_process
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: max_str_ln, dp
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type
    
    implicit none
    private
    public run_post_processing

contains
    ! the main driver routine
    ! [MPI] All ranks
    integer function run_post_processing() result(ierr)
        character(*), parameter :: rout_name = 'run_post_processing'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        call writo('Processing output')
        call lvl_ud(1)
        
        !ierr = process_output(grid_eq_B,eq,grid_X,X_B,n_sol_found,alpha)
        !CHCKERR('')
        
        call lvl_ud(-1)
        call writo('Done processing output')
    end function run_post_processing
    
    ! Processes the output of the simulations for vizualization, analysis, etc.
    integer function process_output(grid_eq,eq,grid_X,X,n_sol_found,alpha) &
        &result(ierr)
        use num_vars, only: no_plots, no_messages
        use grid_vars, only: create_grid, dealloc_grid
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X, create_X
        use grid_ops, only: coord_E2F, calc_XYZ_grid, extend_grid, calc_XYZ_grid
        use X_ops, only: prepare_X
        use sol_ops, only: plot_X_vecs
        
        character(*), parameter :: rout_name = 'process_output'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibirum variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: n_sol_found                                      ! how many solutions found and saved
        real(dp), intent(in) :: alpha                                           ! field line label alpha
        
        ! local variables
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        type(grid_type) :: grid_eq_ext                                          ! extended equilibrium grid, covering whole geometry
        type(grid_type) :: grid_X_ext                                           ! extended perturbation grid, covering whole geometry
        type(eq_type) :: eq_ext                                                 ! extended equilibrium variables
        type(met_type) :: met_ext                                               ! extended metric variables
        type(X_type) :: X_ext                                                   ! extended perturbation variables
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_messages_loc                                              ! local copy of no_messages
        real(dp), allocatable :: X_ind(:,:,:), Y_ind(:,:,:), Z_ind(:,:,:)       ! individual versions of X, Y and Z
        
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
        call writo('creating extended plot grids')
        call lvl_ud(1)
        
        ierr = extend_grid(grid_eq,grid_eq_ext,grid_eq=grid_eq,eq=eq)           ! extend eq grid and convert to F
        CHCKERR('')
        ierr = extend_grid(grid_X,grid_X_ext,grid_eq=grid_eq,eq=eq)             ! extend X grid and convert to F
        CHCKERR('')
        ierr = calc_XYZ_grid(grid_X_ext,X_ind,Y_ind,Z_ind)                      ! calculate X, Y and Z on extended X grid
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! user output
        call writo('plotting the Eigenvectors')
        call lvl_ud(1)
        
        ! Since the equations have been solved  for a single field line normally
        ! the  equilibrium, metric  and perturbation  quantities that  have been
        ! used have  to be recalculated for  the entire angular range,  which is
        ! done  below.  However,  since  here  only  the  equilibrium  variables
        ! flux_q_FD or q_saf_FD or used, the plotting of the Eigenvectors can be
        ! done already.
        ierr = plot_X_vecs(grid_eq_ext,eq_ext,grid_X,X,&
            &reshape([X_ind,Y_ind,Z_ind],&
            &[grid_X_ext%n(1),grid_X_ext%n(2),grid_X_ext%grp_n_r,3]),&
            &n_sol_found,min_id,max_id)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! user output
        call writo('calculating variables on plot grid')
        call lvl_ud(1)
        
        ! back up no_plots and no_messages and set to .true.
        no_plots_loc = no_plots
        no_messages_loc = no_messages
        write(*,*) 'TEMPORARILIY DO NOT BLOCK MESSAGES'
        !!no_plots = .true.
        !!no_messages = .true.
        ! create extended perturbation
        call create_X(grid_eq_ext,X_ext)
        ! create  equilibrium,  metric   and   some  perturbation  variables  on
        ! equilibrium grid
        !!!ierr = calculate_vars_in_eq_grid(grid_eq_ext,eq_ext,met_ext,X_ext,alpha)
        CHCKERR('')
        ! reset no_plots and no_messages
        no_plots = no_plots_loc
        no_messages = no_messages_loc
        ! copy the Eigenvectors and -values into the extended X
        X_ext%val = X%val
        X_ext%vec = X%vec
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Decomposing the energy into its terms')
        call lvl_ud(1)
        
        !!!ierr = decompose_energy(grid_eq,eq,grid_X,X,n_sol_found)
        !!!CHCKERR('')
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Cleaning up')
        call lvl_ud(1)
        
        deallocate(X_ind,Y_ind,Z_ind)
        call dealloc_grid(grid_eq_ext)
        call dealloc_grid(grid_X_ext)
        ierr = dealloc_eq(eq_ext)
        CHCKERR('')
        ierr = dealloc_met(met_ext)
        CHCKERR('')
        call dealloc_X(X_ext)
        
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
                    &log10(abs(realpart(X%val(1:n_sol_found)))),draw=.false.)
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
                        &x=[(id*1._dp,id=1,last_unstable_id)],draw=.false.)
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
                        &x=[(id*1._dp,id=last_unstable_id+1,n_sol_found)],&
                        &draw=.false.)
                    ! same output in file as well
                    call draw_GP(plot_title,'Eigenvalues_'//&
                        &trim(i2str(alpha_job_nr))//'_stable.dat',1,&
                        &.true.,.false.)
                end if
            end if
        end subroutine plot_X_vals
        
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
            use sol_ops, only: calc_real_XUQ
            use grid_ops, only: calc_eqd_grid
            
            character(*), parameter :: rout_name = 'decompose_energy'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(eq_type), intent(in) :: eq                                     ! equilibirum variables
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_type), intent(in) :: X                                       ! perturbation variables
            integer, intent(in) :: n_sol_found                                  ! how many solutions found and saved
            
            !! local variables
            !integer :: jd                                                       ! counter
            !integer :: n_t = 10                                                 ! nr. of time steps
            !real(dp), allocatable :: time(:)                                    ! time at which to calculate decomposition
            !real(dp), allocatable :: X_F(:,:,:,:)                               ! normal comp. of perturbation
            !real(dp), allocatable :: D3X_F(:,:,:,:)                             ! parallel deriv. of X
            !real(dp), allocatable :: U_F(:,:,:,:)                               ! geodesic comp. of perturbation
            !real(dp), allocatable :: D3U_F(:,:,:,:)                             ! parallel deriv. of U
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Merging perturbation and equilibrium grid into custom &
                &grid')
            call lvl_ud(1)
            
            !! intialize variables
            !!allocate(time(n_t))
            !!allocate(X_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            !!allocate(D3X_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            !!allocate(U_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            !!allocate(D3U_F(grid%n(1),grid%n(2),grid%grp_n_r,n_t))
            
            !! user output
            !call writo('Calculating X, U, Qn and Qg')
            !call lvl_ud(1)
            
            !! set up time
            !ierr = calc_eqd_grid(time,0.0_dp,1.0_dp,excl_last=.true.)
            !CHCKERR('')
            
            !! loop over all Eigenvalues
            !do jd = 1,n_sol_found
                !! calculate X, U, Qn and Qg in real space
                !!ierr = calc_real_XUQ(eq,grid_eq,X,jd,time,X_F)
                !CHCKERR('')
            !end do
            
            call lvl_ud(-1)
        end function decompose_energy
    end function process_output
end module post_process
