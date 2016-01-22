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
    use X_vars, only: X_1_type
    use sol_vars, only: sol_type
    
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
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D, retrieve_var_1D_id
        use PB3D_vars, only: vars_1D_eq, vars_1D_sol
        use grid_vars, only: create_grid, dealloc_grid
        use grid_ops, only: calc_norm_range
        use eq_vars, only: create_eq, dealloc_eq
        use eq_ops, only: calc_eq, flux_q_plot, calc_derived_q
        use met_vars, only: create_met, dealloc_met
        use met_ops, only: calc_met, calc_F_derivs
        use X_vars, only: create_X, dealloc_X
        use sol_vars, only: dealloc_sol
        use grid_ops, only: calc_XYZ_grid, extend_grid_E, plot_grid_real
        use X_ops, only: calc_X, resonance_plot, calc_res_surf
        use sol_ops, only: plot_X_vals, plot_X_vec, decompose_energy
        use HELENA, only: interp_HEL_on_grid
        use files_utilities, only: nextunit
        use utilities, only: calc_aux_utilities
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'run_driver_POST'
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B                                   ! field-aligned equilibrium grid
        type(grid_type) :: grid_eq_plot                                         ! plot equilibrium grid
        type(grid_type), target :: grid_sol                                     ! solution grid
        type(grid_type) :: grid_sol_plot                                        ! plot solution grid
        type(eq_type), target :: eq                                             ! equilbrium variables
        type(eq_type), pointer :: eq_B                                          ! field-aligned equilibrium variables
        type(eq_type) :: eq_plot                                                ! plot equilibrium variables
        type(met_type), target :: met                                           ! metric variables
        type(met_type), pointer :: met_B                                        ! field-aligned metric variables
        type(met_type) :: met_plot                                              ! plot metric variables
        type(X_1_type), target :: X_1                                           ! vectorial perturbation variables
        type(X_1_type), pointer :: X_1_B                                        ! field-aligned vectorial perturbation variables
        type(X_1_type) :: X_1_plot                                              ! plot vectorial perturbation variables
        type(sol_type), target :: sol                                           ! solution variables
        integer :: id, jd                                                       ! counter
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        integer :: output_EN_i                                                  ! file number
        integer :: eq_limits(2)                                                 ! i_limit of eq and X variables
        integer :: sol_limits(2)                                                ! i_limit of sol variables
        integer :: r_F_eq_id, r_F_X_id, r_F_sol_id                              ! index of equilibrium, perturbation and solution r_F
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_messages_loc                                              ! local copy of no_messages
        real(dp), allocatable :: XYZ_plot(:,:,:,:)                              ! X, Y and Z on plot grid
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        complex(dp), allocatable :: X_val_comp(:,:,:)                           ! fraction between total E_pot and E_kin, compared with EV
        complex(dp), allocatable :: X_val_comp_loc(:,:,:)                       ! local X_val_comp
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: format_head                                ! header
        !real(dp), allocatable :: plot_var(:,:,:) 
        
        ! initialize ierr
        ierr = 0
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('Reconstruct PB3D output')
        call lvl_ud(1)
        
        ! read from PB3D  output file equilibrium grid  and variables, vectorial
        ! perturbation grid and variables and solution grid and variables
        ierr = read_PB3D(.false.,.true.,.true.,.true.,.true.,.true.,.false.,&
            &.true.)
        CHCKERR('')
        
        ! set eq and X limits
        ierr = retrieve_var_1D_id(vars_1D_eq,'r_F',r_F_eq_id)
        CHCKERR('')
        !!ierr = retrieve_var_1D_id(vars_1D_X,'r_F',r_F_X_id)
        !!CHCKERR('')
        ierr = retrieve_var_1D_id(vars_1D_sol,'r_F',r_F_sol_id)
        CHCKERR('')
        ierr = calc_norm_range(eq_limits=eq_limits,sol_limits=sol_limits,&
            &r_F_eq=vars_1D_eq(r_F_eq_id)%p,r_F_sol=vars_1D_sol(r_F_sol_id)%p)
        CHCKERR('')
        
        ! reconstructing grids depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Reconstruct PB3D output on output grid')
                call lvl_ud(1)
                ! the field-aligned quantities are already found
                grid_eq_B => grid_eq
                eq_B => eq
                met_B => met
                X_1_B => X_1
                ! normal call to reconstruct PB3D
                ierr = reconstruct_PB3D(.false.,.true.,.true.,.true.,.true.,&
                    &.true.,.false.,.true.,grid_eq=grid_eq,grid_sol=grid_sol,&
                    &eq=eq,met=met,X_1=X_1,sol=sol,eq_limits=eq_limits,&
                    &sol_limits=sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('The output grid is field-aligned')
            case (2)                                                            ! HELENA
                ! user output
                call writo('Reconstruct PB3D output on output grid')
                call lvl_ud(1)
                ! the field-aligned quantities are different
                allocate(grid_eq_B)
                allocate(eq_B)
                allocate(met_B)
                allocate(X_1_B)
                ! additionally need field-aligned equilibrium grid
                ierr = reconstruct_PB3D(.false.,.true.,.true.,.true.,.true.,&
                    &.true.,.false.,.true.,grid_eq=grid_eq,grid_eq_B=grid_eq_B,&
                    &grid_sol=grid_sol,eq=eq,met=met,X_1=X_1,sol=sol,&
                    &eq_limits=eq_limits,sol_limits=sol_limits)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('Interpolate all PB3D output on field-aligned grid')
                call lvl_ud(1)
                
                call writo('Preparing quantities')
                call lvl_ud(1)
                ! setup field-aligned quantities
                ierr = create_eq(grid_eq_B,eq_B)
                CHCKERR('')
                ierr = create_met(grid_eq_B,met_B)
                CHCKERR('')
                call create_X(grid_eq_B,X_1_B)
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,eq=eq,&
                    &eq_out=eq_B,met=met,met_out=met_B,&
                    &X_1=X_1,X_1_out=X_1_B,eq_met=eq,&
                    &grid_name='field-aligned grid')
                CHCKERR('')
                
                !! get X, Y and Z of plot
                !allocate(XYZ_plot(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,3))
                !ierr = calc_XYZ_grid(grid_eq,XYZ_plot(:,:,:,1),&
                    !&XYZ_plot(:,:,:,2),XYZ_plot(:,:,:,3))
                !CHCKERR('')
                !allocate(plot_var(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
                !do id = 1,grid_eq%loc_n_r
                    !plot_var(:,:,id) = -eq%pres_FD(id,1)*eq%kappa_n(:,:,id)
                !end do
                !call plot_HDF5('kappa_n_Dp','TEST_kappa_n_Dp',plot_var,&
                    !&x=XYZ_plot(:,:,:,1),y=XYZ_plot(:,:,:,2),&
                    !&z=XYZ_plot(:,:,:,3))
                !do id = 1,grid_eq%loc_n_r
                    !plot_var(:,:,id) = eq%pres_FD(id,1)
                !end do
                !call plot_HDF5('Dp','TEST_Dp',plot_var,x=XYZ_plot(:,:,:,1),&
                    !&y=XYZ_plot(:,:,:,2),z=XYZ_plot(:,:,:,3))
                !call plot_HDF5('kappa_n','kappa_n',eq%kappa_n,&
                    !&x=XYZ_plot(:,:,:,1),y=XYZ_plot(:,:,:,2),&
                    !&z=XYZ_plot(:,:,:,3))
                !call plot_HDF5('kappa_g','TEST_kappa_g',eq%kappa_g,&
                    !&x=XYZ_plot(:,:,:,1),y=XYZ_plot(:,:,:,2),&
                    !&z=XYZ_plot(:,:,:,3))
                !deallocate(XYZ_plot)
                
                ! user output
                call lvl_ud(-1)
                call writo('PB3D output interpolated on field-aligned grid')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        
        ! user output
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! user output
        call writo('Various PB3D outputs')
        call lvl_ud(1)
        
        !!!!!!!!! OKAY !!!!!!!!!!!!!!!!!
        if (plot_resonance) then
            ierr = resonance_plot(eq,grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        if (plot_flux_q) then
            ierr = flux_q_plot(grid_eq,eq)
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        if (plot_grid) then
            ierr = plot_grid_real(grid_eq_B)
            CHCKERR('')
        else
            call writo('Magnetic grid plot not requested')
        end if
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('PB3D outputs done')
        
        ! user output
        call writo('Plot the Eigenvalues')
        call lvl_ud(1)
        
        call plot_X_vals(sol,last_unstable_id)
        
        ! user output
        call lvl_ud(-1)
        call writo('Eigenvalues plotted')
        !!!!!!!!! END OKAY !!!!!!!!!!!!!
        
        ! user output
        call writo('Extend PB3D output to plot grid')
        call lvl_ud(1)
        
        call writo('Setting up the grid')
        call lvl_ud(1)
        ierr = extend_grid_E(grid_eq,grid_eq_plot,grid_eq=grid_eq,eq=eq)        ! extend eq grid and convert to F
        CHCKERR('')
        ierr = extend_grid_E(grid_sol,grid_sol_plot,grid_eq=grid_eq,eq=eq)      ! extend sol grid and convert to F
        CHCKERR('')
        call lvl_ud(-1)
        call writo('Grid set up')
        
        write(*,*) '!!!!!!!!!!!!! NEED GRID_X TO BE ABLE TO CALCULATE CALC_X !!!!'
        ! calculating variables on plot grid depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                call writo('Preparing quantities')
                call lvl_ud(1)
                ! setup plot quantities
                ierr = create_eq(grid_eq_plot,eq_plot)
                CHCKERR('')
                eq_plot%pres_E = eq%pres_E
                eq_plot%q_saf_E = eq%q_saf_E
                eq_plot%rot_t_E = eq%rot_t_E
                eq_plot%flux_p_E = eq%flux_p_E
                eq_plot%flux_t_E = eq%flux_t_E
                eq_plot%rho = eq%rho
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                call writo('Calculating quantities on plot grid')
                call lvl_ud(1)
                ! back up no_plots and no_messages and set to .true.
                no_plots_loc = no_plots; no_plots = .true.
                no_messages_loc = no_messages; no_messages = .true.
                ! Calculate the equilibrium quantities
                ierr = calc_eq(grid_eq_plot,eq_plot)
                CHCKERR('')
                ! Calculate the metric quantities
                ierr = calc_met(grid_eq_plot,eq_plot,met_plot)
                CHCKERR('')
                ! Transform E into F derivatives
#if ldebug
                ierr = calc_F_derivs(grid_eq_plot,eq_plot,met_plot)
#else
                ierr = calc_F_derivs(eq_plot,met_plot)
#endif
                CHCKERR('')
                ! Calculate derived metric quantities
                call calc_derived_q(grid_eq_plot,eq_plot,met_plot)
                ! calculate X variables, vector phase
                !!ierr = calc_X(grid_eq_plot,eq_plot,met_plot,X_1_plot)
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
                ierr = create_eq(grid_eq_plot,eq_plot)
                CHCKERR('')
                ierr = create_met(grid_eq_plot,met_plot)
                CHCKERR('')
                call create_X(grid_eq_plot,X_1_plot)
                call lvl_ud(-1)
                call writo('Quantities prepared')
                
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_plot,eq=eq,&
                    &eq_out=eq_plot,met=met,met_out=met_plot,&
                    &X_1=X_1,X_1_out=X_1_plot,eq_met=eq,&
                    &grid_name='plot grid')
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        
        ! user output
        call lvl_ud(-1)
        call writo('PB3D output extended to plot grid')
        
        ! user output
        call writo('Find stability ranges')
        call lvl_ud(1)
        
        call find_stab_ranges(sol,min_id,max_id,last_unstable_id)
        
        ! user output
        call lvl_ud(-1)
        call writo('Stability ranges found')
        
        ! user output
        call writo('Prepare plots')
        call lvl_ud(1)
        
        call writo('Calculate plot grid')
        call lvl_ud(1)
        allocate(XYZ_plot(grid_sol_plot%n(1),grid_sol_plot%n(2),&
            &grid_sol_plot%loc_n_r,3))
        ierr = calc_XYZ_grid(grid_sol_plot,XYZ_plot(:,:,:,1),&
            &XYZ_plot(:,:,:,2),XYZ_plot(:,:,:,3))
        CHCKERR('')
        call lvl_ud(-1)
        
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        ierr = calc_res_surf(grid_eq,eq,res_surf,info=.false.,&
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
        
        ! user output
        call lvl_ud(-1)
        call writo('Plots prepared')
        
        ! user output
        call writo('Plotting variables for different ranges')
        call lvl_ud(1)
        
        ! loop over three ranges
        allocate(X_val_comp(2,2,0))
        do jd = 1,3
            if (min_id(jd).le.max_id(jd)) call &
                &writo('RANGE '//trim(i2str(jd))//': modes '//&
                &trim(i2str(min_id(jd)))//'..'//trim(i2str(max_id(jd))))
            call lvl_ud(1)
            
            ! indices in each range
            do id = min_id(jd),max_id(jd)
                ! user output
                call writo('Mode '//trim(i2str(id))//'/'//&
                    &trim(i2str(size(sol%val)))//', with eigenvalue '&
                    &//trim(c2strt(sol%val(id))))
                call lvl_ud(1)
                
                call writo('Plot the Eigenvector')
                call lvl_ud(1)
#if ldebug
                ierr = plot_X_vec(grid_eq_plot,grid_sol_plot,eq_plot,met_plot,&
                    &X_1_plot,sol,XYZ_plot,id,res_surf)
#else
                ierr = plot_X_vec(grid_eq_plot,grid_sol_plot,eq_plot,X_1_plot,&
                    &sol,XYZ_plot,id,res_surf)
#endif
                CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call writo('Decompose the energy into its terms')
                call lvl_ud(1)
                allocate(X_val_comp_loc(2,2,size(X_val_comp,3)+1))
                X_val_comp_loc(:,:,1:size(X_val_comp,3)) = X_val_comp
                ierr = decompose_energy(grid_eq_B,grid_sol,eq_B,met_B,&
                    &X_1_B,sol,id,output_EN_i,.true.,&
                    &X_val_comp=X_val_comp_loc(:,:,size(X_val_comp_loc,3)))
                CHCKERR('')
                deallocate(X_val_comp)
                allocate(X_val_comp(2,2,size(X_val_comp_loc,3)))
                X_val_comp = X_val_comp_loc
                deallocate(X_val_comp_loc)
                write(*,*) '!!!! NOT PLOTTING !!!!!!!!!!'
                !!ierr = decompose_energy(grid_eq_plot,grid_sol,eq_plot,met_plot,&
                    !!&X_1_plot,sol,id,output_EN_i,.false.,XYZ=XYZ_plot)
                !!CHCKERR('')
                call lvl_ud(-1)
                
                ! user output
                call lvl_ud(-1)
                call writo('Mode '//trim(i2str(id))//'/'//&
                    &trim(i2str(size(sol%val)))//' finished')
            end do
            
            if (min_id(jd).le.max_id(jd)) call &
                &writo('RANGE '//trim(i2str(jd))//' plotted')
            call lvl_ud(-1)
        end do
        
        ! print the difference between the Eigenvalue and the energy fraction
        call plot_X_val_comp(X_val_comp)
        
        ! user output
        call lvl_ud(-1)
        call writo('Variables for different ranges plotted')
        
        ! clean up 
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_sol)
        call dealloc_eq(eq)
        call dealloc_met(met)
        call dealloc_X(X_1)
        call dealloc_sol(sol)
        if (eq_style.eq.2) then
            call dealloc_grid(grid_eq_B)
            call dealloc_eq(eq_B)
            call dealloc_met(met_B)
            call dealloc_X(X_1_B)
            deallocate(grid_eq_B)
            deallocate(eq_B)
            deallocate(met_B)
            deallocate(X_1_B)
        end if
        nullify(grid_eq_B)
        nullify(eq_B)
        nullify(met_B)
        nullify(X_1_B)
        call dealloc_grid(grid_eq_plot)
        call dealloc_grid(grid_sol_plot)
        call dealloc_eq(eq_plot)
        call dealloc_met(met_plot)
        call dealloc_X(X_1_plot)
        
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
    subroutine find_stab_ranges(sol,min_id,max_id,last_unstable_id)
        use num_vars, only: n_sol_plotted
        
        ! input / output
        type(sol_type), intent(in) :: sol                                       ! solution variables
        integer, intent(inout) :: min_id(3), max_id(3)                          ! min. and max. index of range 1, 2 and 3
        integer, intent(inout) :: last_unstable_id                              ! index of last unstable EV
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_sol_found                                                  ! how many solutions found and saved
        
        ! set local variables
        n_sol_found = size(sol%val)
        
        ! find last unstable index (if ends with 0, no unstable EV)
        last_unstable_id = 0
        do id = 1,n_sol_found
            if (realpart(sol%val(id)).lt.0._dp) last_unstable_id = id
        end do
        ! set up min. and max. of range 1
        if (last_unstable_id.gt.0) then                                         ! there is an unstable range
            if (n_sol_plotted(1).ge.0) then
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
            if (n_sol_plotted(2).ge.0) then
                min_id(2) = last_unstable_id - n_sol_plotted(2) + 1             ! start from n_sol_plotted(2) last unstable values
            else
                min_id(2) = 1                                                   ! start from 1
            end if
            if (n_sol_plotted(3).ge.0) then
                max_id(2) = min(last_unstable_id + n_sol_plotted(3),n_sol_found)! end with n_sol_plotted(3) first stable values if available
            else
                max_id(2) = n_sol_found                                         ! end with last EV
            end if
        else                                                                    ! no unstable range
            min_id(2) = 1                                                       ! start from first EV (stable)
            if (n_sol_plotted(3).ge.0) then
                max_id(2) = min(n_sol_plotted(3),n_sol_found)                   ! end with n_sol_plotted(3) first stable values if available
            else
                max_id(2) = n_sol_found                                         ! end with last EV
            end if
        end if
        ! set up min. and max. of range 3
        if (n_sol_plotted(4).ge.0) then
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
    
    ! plots difference between Eigenvalues and energy fraction
    subroutine plot_X_val_comp(X_val_comp)
        use num_vars, only: rank
        
        ! input /output
        complex(dp) :: X_val_comp(:,:,:)                                        ! fraction between total E_pot and E_kin, compared with EV
        
        ! local variables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        
        if (rank.eq.0) then
            ! real part
            plot_title = 'RE X_val and E_frac'
            plot_name = 'X_val_comp_RE'
            call print_GP_2D(plot_title,plot_name,realpart(X_val_comp(:,1,:)),&
                &x=realpart(X_val_comp(:,2,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,size(X_val_comp,3),1,&
                &.false.)
            plot_title = 'RE X_val and E_frac rel diff'
            plot_name = 'X_val_comp_RE_rel_diff'
            call print_GP_2D(plot_title,plot_name,realpart(&
                &(X_val_comp(1,2,:)-X_val_comp(2,2,:))/X_val_comp(1,2,:)),&
                &x=realpart(X_val_comp(1,1,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            plot_title = 'RE X_val and E_frac log rel diff'
            plot_name = 'X_val_comp_RE_rel_diff_log'
            call print_GP_2D(plot_title,plot_name,log10(abs(realpart(&
                &(X_val_comp(1,2,:)-X_val_comp(2,2,:))/X_val_comp(1,2,:)))),&
                &x=realpart(X_val_comp(1,1,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            
            ! imaginary part
            plot_title = 'IM X_val and E_frac'
            plot_name = 'X_val_comp_IM'
            call print_GP_2D(plot_title,plot_name,realpart(X_val_comp(:,1,:)),&
                &x=imagpart(X_val_comp(:,2,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,size(X_val_comp,3),1,&
                &.false.)
            plot_title = 'IM X_val and E_frac rel diff'
            plot_name = 'X_val_comp_IM_rel_diff'
            call print_GP_2D(plot_title,plot_name,imagpart(&
                &(X_val_comp(1,2,:)-X_val_comp(2,2,:))/X_val_comp(1,2,:)),&
                &x=realpart(X_val_comp(1,1,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            plot_title = 'IM X_val and E_frac log rel diff'
            plot_name = 'X_val_comp_IM_rel_diff_log'
            call print_GP_2D(plot_title,plot_name,log10(abs(imagpart(&
                &(X_val_comp(1,2,:)-X_val_comp(2,2,:))/X_val_comp(1,2,:)))),&
                &x=realpart(X_val_comp(1,1,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
        end if
    end subroutine plot_X_val_comp
end module driver_POST
