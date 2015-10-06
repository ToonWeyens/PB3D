!------------------------------------------------------------------------------!
!   Main driver of PostProcessing of program Peeling Ballooning in 3D          !
!------------------------------------------------------------------------------!
module driver_POST
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: max_str_ln, dp, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type
    use PB3D_vars, only: dealloc_PB3D, &
        &PB3D_type
    
    implicit none
    private
    public run_driver_POST

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
    integer function run_driver_POST() result(ierr)
        use num_vars, only: no_messages, no_plots, eq_style, plot_resonance, &
            &plot_flux_q, plot_grid, prog_name, output_name, rank
        use PB3D_ops, only: reconstruct_PB3D
        use grid_vars, only: create_grid
        use eq_vars, only: create_eq
        use met_vars, only: create_met
        use X_vars, only: create_X
        use grid_ops, only: calc_XYZ_grid, extend_grid_E, plot_grid_real
        use eq_ops, only: calc_eq, flux_q_plot
        use X_ops, only: prepare_X, resonance_plot, calc_res_surf
        use sol_ops, only: plot_X_vec, decompose_energy
        use HELENA, only: interp_HEL_on_grid
        use files_utilities, only: nextunit
        use utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_POST'
        
        ! local variables
        type(PB3D_type), target :: PB3D                                         ! output PB3D for which to do postprocessing
        type(PB3D_type), pointer :: PB3D_B => null()                            ! PB3D variables on a field-aligned grid
        type(PB3D_type) :: PB3D_plot                                            ! PB3D variables on a plot grid
        integer :: id, jd                                                       ! counter
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        integer :: output_EN_i                                                  ! file number
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_messages_loc                                              ! local copy of no_messages
        real(dp), allocatable :: X_plot(:,:,:), Y_plot(:,:,:), Z_plot(:,:,:)    ! X, Y and Z on plot grid
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: format_head                                ! header
        
        ! initialize ierr
        ierr = 0
        
        ! calculate auxiliary quantities for utilities
        call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! reconstructing grids depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Reconstruct PB3D output on output grid')
                call lvl_ud(1)
                ! the field-aligned grid is identical to the output grid
                PB3D_B => PB3D
                ! normal call to reconstruct_PB3D
                ierr = reconstruct_PB3D(.true.,.true.,.true.,PB3D)
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
                ierr = reconstruct_PB3D(.true.,.true.,.true.,PB3D,&
                    &PB3D_B%grid_eq)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('Interpolate all PB3D output on field-aligned grid')
                call lvl_ud(1)
                
                call writo('Preparing quantities')
                call lvl_ud(1)
                ! setup field-aligned quantities
                ierr = create_eq(PB3D_B%grid_eq,PB3D_B%eq)
                CHCKERR('')
                ierr = create_met(PB3D_B%grid_eq,PB3D_B%met)
                CHCKERR('')
                ierr = create_grid(PB3D_B%grid_X,PB3D%grid_X%n,&
                    &[PB3D%grid_X%i_min,PB3D%grid_X%i_max])
                CHCKERR('')
                PB3D_B%grid_X%r_F = PB3D%grid_X%r_F
                PB3D_B%grid_X%r_e = PB3D%grid_X%r_E
                PB3D_B%grid_X%loc_r_F = PB3D%grid_X%loc_r_F
                PB3D_B%grid_X%loc_r_e = PB3D%grid_X%loc_r_E
                call create_X(PB3D_B%grid_eq,PB3D_B%X)
                allocate(PB3D_B%X%val(size(PB3D%X%val)))
                allocate(PB3D_B%X%vec(size(PB3D%X%vec,1),size(PB3D%X%vec,2),&
                    &size(PB3D%X%vec,3)))
                PB3D_B%X%val = PB3D%X%val
                PB3D_B%X%vec = PB3D%X%vec
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(PB3D%grid_eq,PB3D_B%grid_eq,PB3D%X,PB3D_B%X,&
                    &met=PB3D%met,met_B=PB3D_B%met,eq=PB3D%eq,eq_B=PB3D_B%eq,&
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
        call writo('Various PB3D outputs')
        call lvl_ud(1)
        
        if (plot_resonance) then
            ierr = resonance_plot(PB3D%eq,PB3D%grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        if (plot_flux_q) then
            ierr = flux_q_plot(PB3D%grid_eq,PB3D%eq)
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        if (plot_grid) then
            ierr = plot_grid_real(PB3D_B%grid_eq)
            CHCKERR('')
        else
            call writo('Magnetic grid plot not requested')
        end if
        
        call lvl_ud(-1)
        
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
        
        ! calculating variables on plot grid depends on equilibrium style
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
                PB3D_plot%eq%rho = PB3D%eq%rho
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                call writo('Calculating quantities on plot grid')
                call lvl_ud(1)
                ! back up no_plots and no_messages and set to .true.
                no_plots_loc = no_plots; no_plots = .true.
                no_messages_loc = no_messages; no_messages = .true.
                ! calculate the non-flux equilibrium quantities for plot grid
                !!ierr = calc_eq(PB3D_plot%grid_eq,PB3D_plot%eq,PB3D_plot%met)
                !!CHCKERR('')
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
                    &PB3D%X,PB3D_plot%X,met=PB3D%met,met_B=PB3D_plot%met,eq=PB3D%eq,&
                    &eq_B=PB3D_plot%eq,grid_name='plot grid')
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        
        ! copy the Eigenvectors and -values into the plot X
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
        call writo('Prepare plots')
        call lvl_ud(1)
        
        call writo('Calculate plot grid')
        call lvl_ud(1)
        allocate(X_plot(PB3D_plot%grid_X%n(1),PB3D_plot%grid_X%n(2),&
            &PB3D_plot%grid_X%loc_n_r))
        allocate(Y_plot(PB3D_plot%grid_X%n(1),PB3D_plot%grid_X%n(2),&
            &PB3D_plot%grid_X%loc_n_r))
        allocate(Z_plot(PB3D_plot%grid_X%n(1),PB3D_plot%grid_X%n(2),&
            &PB3D_plot%grid_X%loc_n_r))
        ierr = calc_XYZ_grid(PB3D_plot%grid_X,X_plot,Y_plot,Z_plot)
        CHCKERR('')
        call lvl_ud(-1)
        
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        ierr = calc_res_surf(PB3D%grid_eq,PB3D%eq,res_surf,info=.false.,&
            &tol_NR=1.E-8_dp,max_it_NR=5000)
        call lvl_ud(-1)
        
        call writo('Open decomposition log file')
        call lvl_ud(1)
        if (rank.eq.0) then
            ! set format strings
            format_head = '("#  ",A23," ",A23," ",A23," ",A23," ",A23," ",A23)'
            ! open output file for the log
            full_output_name = prog_name//'_'//output_name//'_EN.txt'
            open(unit=nextunit(output_EN_i),file=full_output_name,&
                &iostat=ierr)
            CHCKERR('Cannot open EN output file')
            call writo('Log file opened in '//trim(full_output_name))
            write(output_EN_i,'(A)') '# Energy decomposition using the &
                &solution Eigenvectors:'
            write(output_EN_i,format_head) &
                &'RE Eigenvalue          ', 'RE E_pot/E_kin         ', &
                &'RE E_pot               ', 'RE E_kin               '
            write(output_EN_i,format_head) &
                &'IM Eigenvalue          ', 'IM E_pot/E_kin         ', &
                &'IM E_pot               ', 'IM E_kin               '
            write(output_EN_i,format_head) &
                &'RE E_kin_n             ', 'RE E_kin_g             '
            write(output_EN_i,format_head) &
                &'IM E_kin_n             ', 'IM E_kin_g             '
            write(output_EN_i,format_head) &
                &'RE E_pot line_bending_n', 'RE E_pot line_bending_g', &
                &'RE E_pot ballooning_n  ', 'RE E_pot ballooning_g  ', &
                &'RE E_pot kink_n        ', 'RE E_pot kink_g        '
            write(output_EN_i,format_head) &
                &'IM E_pot line_bending_n', 'IM E_pot line_bending_g', &
                &'IM E_pot ballooning_n  ', 'IM E_pot ballooning_g  ', &
                &'IM E_pot kink_n        ', 'IM E_pot kink_g        '
        end if
        call lvl_ud(-1)
        
        call lvl_ud(-1)
        
        ! loop over three ranges
        do jd = 1,3
            if (min_id(jd).le.max_id(jd)) call &
                &writo('RANGE '//trim(i2str(jd))//': modes '//&
                &trim(i2str(min_id(jd)))//'..'//trim(i2str(max_id(jd))))
            call lvl_ud(1)
            
            ! indices in each range
            do id = min_id(jd),max_id(jd)
                ! user output
                call writo('Mode '//trim(i2str(id))//'/'//&
                    &trim(i2str(size(PB3D%X%val)))//', with eigenvalue '&
                    &//trim(c2strt(PB3D%X%val(id))))
                call lvl_ud(1)
                
                call writo('Plot the Eigenvector')
                call lvl_ud(1)
                ierr = plot_X_vec(PB3D_plot%grid_eq,PB3D_plot%eq,&
                    &PB3D_plot%grid_X,PB3D_plot%X,&
                    &reshape([X_plot,Y_plot,Z_plot],[PB3D_plot%grid_X%n(1),&
                    &PB3D_plot%grid_X%n(2),PB3D_plot%grid_X%loc_n_r,3]),id,&
                    &res_surf)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('Decompose the energy into its terms')
                call lvl_ud(1)
                ierr = decompose_energy(PB3D_B,id,output_EN_i,PB3D_plot,&
                    &reshape([X_plot,Y_plot,Z_plot],[PB3D_plot%grid_X%n(1),&
                    &PB3D_plot%grid_X%n(2),PB3D_plot%grid_X%loc_n_r,3]))
                CHCKERR('')
                call lvl_ud(-1)
                
                call lvl_ud(-1)
            end do
            
            call lvl_ud(-1)
        end do
        
        ! clean up 
        call dealloc_PB3D(PB3D)
        call dealloc_PB3D(PB3D_plot)
        if (eq_style.eq.2) call dealloc_PB3D(PB3D_B)                            ! only for HELENA
        nullify(PB3D_B)
        
        ! close output
        if (rank.eq.0) close(output_EN_i)
    end function run_driver_POST
    
    ! finds the plot ranges min_id and max_id
    ! There are three ranges, calculated using n_sol_plotted, which indicates:
    !   1. how many of the first EV's in the unstable range
    !   2. how many of the last EV's in the unstable range
    !   3. how many of the first EV's in the stable range
    !   4. how many of the last EV's in the stable range
    ! have  to be  plotted. This  yields maximally  three different  ranges: One
    ! starting at  the first unstable  EV, one centered  around the zero  of the
    ! EV's and one  ending at the last  EV. These ranges can be  disjoint but do
    ! not have  to be. Also,  it is  possible that a  range does not  exist, for
    ! example if there are no unstable EV's.
    ! Note: A negative value for the elements in n_sol_plotted means "all values
    ! in range":
    !   1. or 2. full unstable range
    !   3. or 4. full stable range
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
            if (n_sol_plotted(1).gt.0) then
                min_id(1) = 1                                                   ! start from most unstable EV
                max_id(1) = min(n_sol_plotted(1),last_unstable_id)              ! end with n_sol_plotted(1) first unstable values if available
            else
                min_id(1) = 1                                                   ! start from most unstable EV
                max_id(1) = last_unstable_id                                    ! end with last unstable EV
            end if
        else                                                                    ! no unstable range
            min_id(1) = 1                                                       ! no unstable values to plot
            max_id(1) = 0                                                       ! no unstable values to plot
        end if
        ! set up min. and max. of range 2
        if (last_unstable_id.gt.0) then                                         ! there is an unstable range
            if (n_sol_plotted(2).gt.0) then
                min_id(2) = last_unstable_id - n_sol_plotted(2) + 1             ! start from n_sol_plotted(2) last unstable values
            else
                min_id(2) = 1                                                   ! start from 1
            end if
            if (n_sol_plotted(3).gt.0) then
                max_id(2) = min(last_unstable_id + n_sol_plotted(3),n_sol_found)! end with n_sol_plotted(3) first stable values if available
            else
                max_id(2) = n_sol_found                                         ! end with last EV
            end if
        else                                                                    ! no unstable range
            min_id(2) = 1                                                       ! start from first EV (stable)
            if (n_sol_plotted(3).gt.0) then
                max_id(2) = min(n_sol_plotted(3),n_sol_found)                   ! end with n_sol_plotted(3) first stable values if available
            else
                max_id(2) = n_sol_found                                         ! end with last EV
            end if
        end if
        ! set up min. and max. of range 3
        if (n_sol_plotted(4).gt.0) then
            min_id(3) = n_sol_found - n_sol_plotted(4) + 1                      ! start from n_sol_plotted(4) last stable values
        else
            min_id(3) = last_unstable_id + 1                                    ! start from n_sol_plotted(4) last stable values
        end if
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
        use num_vars, only: rank
        
        ! input / output
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: last_unstable_id                                 ! index of last unstable EV
        
        ! local variables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        integer :: n_sol_found                                                  ! how many solutions found and saved
        integer :: id                                                           ! counter
        
        ! set local variables
        n_sol_found = size(X%val)
        
        ! only let master plot
        if (rank.eq.0) then
            ! Last Eigenvalues
            ! output on screen
            plot_title = 'final Eigenvalues omega^2 [log]'
            plot_name = 'Eigenvalues'
            call print_GP_2D(plot_title,plot_name,&
                &log10(abs(realpart(X%val(1:n_sol_found)))),draw=.false.)
            ! same output in file as well
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            
            ! Last Eigenvalues: unstable range
            if (last_unstable_id.gt.0) then
                ! output on screen
                plot_title = 'final unstable Eigenvalues omega^2'
                plot_name = 'Eigenvalues_unstable'
                call print_GP_2D(plot_title,plot_name,&
                    &realpart(X%val(1:last_unstable_id)),&
                    &x=[(id*1._dp,id=1,last_unstable_id)],draw=.false.)
                ! same output in file as well
                call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            end if
            
            ! Last Eigenvalues: stable range
            if (last_unstable_id.lt.n_sol_found) then
                ! output on screen
                plot_title = 'final stable Eigenvalues omega^2'
                plot_name = 'Eigenvalues_stable'
                call print_GP_2D(plot_title,plot_name,&
                    &realpart(X%val(last_unstable_id+1:n_sol_found)),&
                    &x=[(id*1._dp,id=last_unstable_id+1,n_sol_found)],&
                    &draw=.false.)
                ! same output in file as well
                call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            end if
        end if
    end subroutine plot_X_vals
end module driver_POST
