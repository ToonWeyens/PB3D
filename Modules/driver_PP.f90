!------------------------------------------------------------------------------!
!   Main driver of PostProcessing of program Peeling Ballooning in 3D          !
!------------------------------------------------------------------------------!
module driver_PP
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: max_str_ln, dp, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type
        use PB3D_vars, only: PB3D_type
    
    implicit none
    private
    public run_driver_PP

contains
    ! The main driver routine for postprocessing
    ! Note  that  the PB3D  output  is given  on different  grids for  different
    ! styles of the equilibrium model:
    !   - VMEC: field-aligned grid on which EV problem has been solved.
    !   - HELENA: output  grid  on  which equilibrium  and metric  variables are
    !   tabulated.
    ! These variables are needed here both in a field-aligned and a plot grid. A
    ! consequence of  the above  is that  for VMEC  the field-aligned  output is
    ! already given, but the  output on the plot grid has  to be calculated from
    ! scratch, while  for HELENA both outputs  have to be interpolated  from the
    ! output tables.
    ! [MPI] All ranks
    integer function run_driver_PP() result(ierr)
        use num_vars, only: no_messages, no_plots, eq_style
        use PB3D_ops, only: reconstruct_PB3D
        use grid_vars, only: create_grid
        use grid_ops, only: calc_XYZ_grid, extend_grid_E
        use eq_vars, only: create_eq
        use met_vars, only: create_met
        use X_vars, only: create_X
        use X_ops, only: prepare_X
        use sol_ops, only: plot_X_vecs, decompose_energy
        use eq_ops, only: calc_eq
        use HELENA, only: interp_HEL_on_grid
        
        character(*), parameter :: rout_name = 'run_driver_PP'
        
        ! local variables
        type(PB3D_type), target :: PB3D                                         ! output PB3D for which to do postprocessing
        type(PB3D_type), pointer :: PB3D_B => null()                            ! PB3D variables on a field-aligned grid
        type(PB3D_type) :: PB3D_plot                                            ! PB3D variables on a plot grid
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_messages_loc                                              ! local copy of no_messages
        real(dp), allocatable :: X_plot(:,:,:), Y_plot(:,:,:), Z_plot(:,:,:)    ! X, Y and Z on plot grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! reconstructing grids depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Reconstruct PB3D output on output grid')
                call lvl_ud(1)
                ! the field-aligned grid is identical to the output grid
                PB3D_B => PB3D
                ! normal call to reconstruct_PB3D
                ierr = reconstruct_PB3D(PB3D)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('The output grid is field-aligned')
            case (2)                                                            ! HELENA
                ! user output
                call writo('Reconstruct PB3D output on output grid')
                call lvl_ud(1)
                ! the field-aligned grid is different form the output grid
                allocate(PB3D_B)
                ! additionally need field-aligned equilibrium grid
                ierr = reconstruct_PB3D(PB3D,PB3D_B%grid_eq)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('Interpolate PB3D output on field-aligned grid')
                call lvl_ud(1)
                
                call writo('Preparing quantities')
                call lvl_ud(1)
                ! setup field-aligned quantities
                ierr = create_eq(PB3D_B%grid_eq,PB3D_B%eq)
                CHCKERR('')
                ierr = create_met(PB3D_B%grid_eq,PB3D_B%met)
                CHCKERR('')
                call create_X(PB3D_B%grid_eq,PB3D_B%X)
                call lvl_ud(-1)
                call writo('Quantities prepared')
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(PB3D%grid_eq,PB3D_B%grid_eq,PB3D%met,&
                    &PB3D_B%met,PB3D%X,PB3D_B%X,eq=PB3D%eq,eq_B=PB3D_B%eq,&
                    &grid_name='field-aligned grid')
                CHCKERR('')
                
                call lvl_ud(-1)
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        
        ! user output
        call writo('Extend PB3D output to plot grid')
        call lvl_ud(1)
        
        call writo('Setting up the grid')
        call lvl_ud(1)
        ierr = extend_grid_E(PB3D%grid_eq,PB3D_plot%grid_eq,&
            &grid_eq=PB3D%grid_eq,eq=PB3D%eq)                                   ! extend eq grid and convert to F
        CHCKERR('')
        ierr = extend_grid_E(PB3D%grid_X,PB3D_plot%grid_X,&
            &grid_eq=PB3D%grid_eq,eq=PB3D%eq)                                   ! extend X grid and convert to F
        CHCKERR('')
        call lvl_ud(-1)
        call writo('Grid set up')
        
        ! depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                call writo('Preparing quantities')
                call lvl_ud(1)
                ! setup plot quantities
                ierr = create_eq(PB3D_plot%grid_eq,PB3D_plot%eq)
                CHCKERR('')
                PB3D_plot%eq%pres_E = PB3D%eq%pres_E
                PB3D_plot%eq%q_saf_E = PB3D%eq%q_saf_E
                PB3D_plot%eq%rot_t_E = PB3D%eq%rot_t_E
                PB3D_plot%eq%flux_p_E = PB3D%eq%flux_p_E
                PB3D_plot%eq%flux_t_E = PB3D%eq%flux_t_E
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                call writo('Calculating quantities on plot grid')
                call lvl_ud(1)
                ! back up no_plots and no_messages and set to .true.
                no_plots_loc = no_plots; no_plots = .true.
                no_messages_loc = no_messages; no_messages = .true.
                ! calculate the non-flux equilibrium quantities for plot grid
                ierr = calc_eq(PB3D_plot%grid_eq,PB3D_plot%eq,PB3D_plot%met)
                CHCKERR('')
                ! prepare matrix elements for plot grid
                ierr = prepare_X(PB3D_plot%grid_eq,PB3D_plot%eq,PB3D_plot%met,&
                    &PB3D_plot%X)
                CHCKERR('')
                ! reset no_plots and no_messages
                no_plots = no_plots_loc
                no_messages = no_messages_loc
                call lvl_ud(-1)
                call writo('Quantities calculated')
            case (2)                                                            ! HELENA
                call writo('Preparing quantities')
                call lvl_ud(1)
                ! setup plot quantities
                ierr = create_eq(PB3D_plot%grid_eq,PB3D_plot%eq)
                CHCKERR('')
                ierr = create_met(PB3D_plot%grid_eq,PB3D_plot%met)
                CHCKERR('')
                call create_X(PB3D_plot%grid_eq,PB3D_plot%X)
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(PB3D%grid_eq,PB3D_plot%grid_eq,&
                    &PB3D%met,PB3D_plot%met,PB3D%X,PB3D_plot%X,eq=PB3D%eq,&
                    &eq_B=PB3D_plot%eq,grid_name='plot grid')
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        
        ! copy the Eigenvectors and -values into the extended X
        allocate(PB3D_plot%X%val(size(PB3D%X%val)))
        allocate(PB3D_plot%X%vec(size(PB3D%X%vec,1),size(PB3D%X%vec,2),&
            &size(PB3D%X%vec,3)))
        PB3D_plot%X%val = PB3D%X%val
        PB3D_plot%X%vec = PB3D%X%vec
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Find stability ranges')
        call lvl_ud(1)
        
        call find_stab_ranges(PB3D%X,min_id,max_id,last_unstable_id)
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Plot the Eigenvalues')
        call lvl_ud(1)
        
        call plot_X_vals(PB3D%X,last_unstable_id)
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Plot the Eigenvectors on plot grid')
        call lvl_ud(1)
        
        ! TEMPORARY: Since  the equations  have been solved  for a  single field
        ! line normally the equilibrium, metric and perturbation quantities that
        ! have been used  have to be recalculated for the  entire angular range,
        ! which  is  done  below.  However,  since  here  only  the  equilibrium
        ! variables  flux_q_FD  or  q_saf_FD  or   used,  the  plotting  of  the
        ! Eigenvectors can be done already.
        ierr = calc_XYZ_grid(PB3D_plot%grid_X,X_plot,Y_plot,Z_plot)             ! calculate X, Y and Z on plot grid
        CHCKERR('')
        ierr = plot_X_vecs(PB3D_plot%grid_eq,PB3D_plot%eq,PB3D_plot%grid_X,&
            &PB3D_plot%X,reshape([X_plot,Y_plot,Z_plot],[PB3D_plot%grid_X%n(1),&
            &PB3D_plot%grid_X%n(2),PB3D_plot%grid_X%grp_n_r,3]),min_id,max_id)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Decompose the energy into its terms')
        call lvl_ud(1)
        
        ierr = decompose_energy(PB3D%grid_eq,PB3D%eq)
        CHCKERR('')
        
        call lvl_ud(-1)
    end function run_driver_PP
    
    ! finds the plot ranges min_id and max_id
    subroutine find_stab_ranges(X,min_id,max_id,last_unstable_id)
        use num_vars, only: n_sol_plotted
        
        ! input / output
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(inout) :: min_id(3), max_id(3)                          ! min. and max. index of range 1, 2 and 3
        integer, intent(inout) :: last_unstable_id                              ! index of last unstable EV
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_sol_found                                                  ! how many solutions found and saved
        
        ! set local variables
        n_sol_found = size(X%val)
        
        ! find last unstable index (if ends with 0, no unstable EV)
        last_unstable_id = 0
        do id = 1,n_sol_found
            if (realpart(X%val(id)).lt.0._dp) last_unstable_id = id
        end do
        ! set up min. and max. of range 1
        if (last_unstable_id.gt.0) then                                         ! there is an unstable range
            min_id(1) = 1                                                       ! start from most unstable EV
            max_id(1) = min(n_sol_plotted(1),last_unstable_id)                  ! end with n_sol_plotted(1) first unstable values if available
        else                                                                    ! no unstable range
            min_id(1) = 1                                                       ! no unstable values to plot
            max_id(1) = 0                                                       ! no unstable values to plot
        end if
        ! set up min. and max. of range 2
        if (last_unstable_id.gt.0) then                                         ! there is an unstable range
            min_id(2) = last_unstable_id - n_sol_plotted(2) + 1                 ! start from n_sol_plotted(2) last unstable values
            max_id(2) = &
                &min(last_unstable_id + n_sol_plotted(3),n_sol_found)           ! end with n_sol_plotted(3) first stable values if available
        else                                                                    ! no unstable range
            min_id(2) = 1                                                       ! start from first EV (stable)
            max_id(2) = min(n_sol_plotted(3),n_sol_found)                       ! end with n_sol_plotted(3) first stable values if available
        end if
        ! set up min. and max. of range 3
        min_id(3) = n_sol_found - n_sol_plotted(4) + 1                          ! start from n_sol_plotted(4) last stable values
        max_id(3) = n_sol_found                                                 ! end with most stable EV
        ! merge ranges 2 and 3 if overlap
        if (min_id(3).le.max_id(2)) then
            max_id(2) = max_id(3)                                               ! range 3 merged with range 2
            min_id(3) = 1                                                       ! range 3 not used any more
            max_id(3) = 0                                                       ! range 3 not used any more
        end if
        ! merge ranges 1 and 2 if overlap
        if (min_id(2).le.max_id(1)) then
            max_id(1) = max_id(2)                                               ! range 2 merged with range 1
            min_id(2) = 1                                                       ! range 2 not used any more
            max_id(2) = 0                                                       ! range 2 not used any more
        end if
    end subroutine find_stab_ranges
    
    ! plots Eigenvalues
    subroutine plot_X_vals(X,last_unstable_id)
        use num_vars, only: grp_rank
        
        ! input / output
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: last_unstable_id                                 ! index of last unstable EV
        
        ! local variables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        integer :: n_sol_found                                                  ! how many solutions found and saved
        integer :: id                                                           ! counter
        
        ! set local variables
        n_sol_found = size(X%val)
        
        ! only let group rank plot
        if (grp_rank.eq.0) then
            ! Last Eigenvalues
            ! output on screen
            plot_title = 'final Eigenvalues omega^2 [log]'
            call print_GP_2D(plot_title,'Eigenvalues.dat',&
                &log10(abs(realpart(X%val(1:n_sol_found)))),draw=.false.)
            ! same output in file as well
            call draw_GP(plot_title,'Eigenvalues.dat',1,.true.,.false.)
            
            ! Last Eigenvalues: unstable range
            if (last_unstable_id.gt.0) then
                ! output on screen
                plot_title = 'final unstable Eigenvalues omega^2'
                call print_GP_2D(plot_title,'Eigenvalues_unstable.dat',&
                    &realpart(X%val(1:last_unstable_id)),&
                    &x=[(id*1._dp,id=1,last_unstable_id)],draw=.false.)
                ! same output in file as well
                call draw_GP(plot_title,'Eigenvalues_unstable.dat',1,&
                    &.true.,.false.)
            end if
            
            ! Last Eigenvalues: stable range
            if (last_unstable_id.lt.n_sol_found) then
                ! output on screen
                plot_title = 'final stable Eigenvalues omega^2'
                call print_GP_2D(plot_title,'Eigenvalues_stable.dat',&
                    &realpart(X%val(last_unstable_id+1:n_sol_found)),&
                    &x=[(id*1._dp,id=last_unstable_id+1,n_sol_found)],&
                    &draw=.false.)
                ! same output in file as well
                call draw_GP(plot_title,'Eigenvalues_stable.dat',1,&
                    &.true.,.false.)
            end if
        end if
    end subroutine plot_X_vals
end module driver_PP
