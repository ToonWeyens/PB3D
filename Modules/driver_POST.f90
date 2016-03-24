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
    use eq_vars, only: eq_1_type, eq_2_type
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
        use num_vars, only: no_output, no_plots, eq_style, plot_resonance, &
            &plot_flux_q, plot_grid, prog_name, output_name, rank
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_1, reconstruct_PB3D_sol
        use PB3D_utilities, only: retrieve_var_1D_id
        use grid_vars, only: create_grid, dealloc_grid
        use grid_ops, only: calc_norm_range, plot_grid_real
        use eq_vars, only: create_eq, dealloc_eq
        use eq_ops, only: calc_eq, flux_q_plot, calc_derived_q
        use eq_utilities, only: calc_F_derivs
        use X_vars, only: create_X, dealloc_X
        use sol_vars, only: dealloc_sol
        use grid_utilities, only: calc_XYZ_grid, extend_grid_E
        use X_ops, only: calc_X, resonance_plot, calc_res_surf, setup_nm_X
        use sol_ops, only: plot_sol_vals, plot_sol_vec, decompose_energy
        use HELENA_ops, only: interp_HEL_on_grid
        use files_utilities, only: nextunit
        use input_utilities, only: dealloc_in
        use utilities, only: calc_aux_utilities
        use MPI_utilities, only: wait_MPI
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'run_driver_POST'
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B                                   ! field-aligned equilibrium grid
        type(grid_type) :: grid_eq_plot                                         ! plot equilibrium grid
        type(grid_type), target :: grid_X                                       ! perturbation grid
        type(grid_type), pointer :: grid_X_B                                    ! field-aligned perturbation grid
        type(grid_type) :: grid_X_plot                                          ! plot perturbation grid
        type(grid_type), target :: grid_sol                                     ! solution grid
        type(grid_type) :: grid_sol_plot                                        ! plot solution grid
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium
        type(eq_2_type), target :: eq_2                                         ! metric equilibrium
        type(eq_2_type), pointer :: eq_2_B                                      ! field-aligned metric equilibrium
        type(eq_1_type) :: eq_1_plot                                            ! plot flux equilibrium
        type(eq_2_type) :: eq_2_plot                                            ! plot metric equilibrium
        type(X_1_type), target :: X                                             ! vectorial perturbation variables
        type(X_1_type), pointer :: X_B                                          ! field-aligned vectorial perturbation variables
        type(X_1_type) :: X_plot                                                ! plot vectorial perturbation variables
        type(sol_type), target :: sol                                           ! solution variables
        integer :: id, jd                                                       ! counter
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        integer :: output_EN_i                                                  ! file number
        integer :: eq_limits(2)                                                 ! i_limit of eq variables
        integer :: X_limits(2)                                                  ! i_limit of X variables
        integer :: sol_limits(2)                                                ! i_limit of sol variables
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_output_loc                                                ! local copy of no_output
        real(dp), allocatable :: XYZ_plot(:,:,:,:)                              ! X, Y and Z on plot grid
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        complex(dp), allocatable :: sol_val_comp(:,:,:)                         ! fraction between total E_pot and E_kin, compared with EV
        complex(dp), allocatable :: sol_val_comp_loc(:,:,:)                     ! local sol_val_comp
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: format_head                                ! header
        !real(dp), allocatable :: plot_var(:,:,:) 
        
        ! initialize ierr
        ierr = 0
        
        ! some preliminary things
        ierr = wait_MPI()
        CHCKERR('')
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! set up whether Richardson level has to be appended to the name
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('Setting up grid limits')
        call lvl_ud(1)
        
        ! reconstruct full equilibrium, perturbation and solution grids
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_sol,'sol',rich_lvl=rich_lvl)
        CHCKERR('')
        
        ! set eq and X limits, using r_F of the grids
        ierr = calc_norm_range(eq_limits=eq_limits,X_limits=X_limits,&
            &sol_limits=sol_limits,r_F_eq=grid_eq%r_F,r_F_X=grid_X%r_F,&
            &r_F_sol=grid_sol%r_F)
        CHCKERR('')
        
        ! deallocate the grids
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_grid(grid_sol)
        
        ! user output
        call lvl_ud(-1)
        call writo('Limits set up')
        
        ! user output
        call writo('Reconstruct PB3D output')
        call lvl_ud(1)
        
        ! reconstruct local equilibrium, perturbation and solution grids
        call writo('Reconstructing output grids')
        call lvl_ud(1)
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name,&
            &grid_limits=eq_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
            &grid_limits=X_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_sol,'sol',rich_lvl=rich_lvl,&
            &grid_limits=sol_limits)
        CHCKERR('')
        call lvl_ud(-1)
        call writo('Output grids reconstructed')
        
        ! user output
        call writo('Reconstructing PB3D output on output grid')
        call lvl_ud(1)
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1',eq_limits=eq_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_2(grid_eq,eq_2,'eq_2',&
            &rich_lvl=rich_lvl_name,eq_limits=eq_limits)
        CHCKERR('')
        ierr = setup_nm_X(grid_eq,grid_X,eq_1,plot_nm=.true.)                   ! is necessary for X variables
        CHCKERR('')
        ierr = reconstruct_PB3D_X_1(grid_X,X,'X_1',rich_lvl=rich_lvl_name,&
            &X_limits=X_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl,&
            &sol_limits=sol_limits)
        CHCKERR('')
        ! reconstructing grids depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! the field-aligned quantities are already found
                grid_eq_B => grid_eq
                grid_X_B => grid_X
                eq_2_B => eq_2
                X_B => X
            case (2)                                                            ! HELENA
                ! also need field-aligned quantities
                allocate(grid_eq_B)
                allocate(grid_X_B)
                allocate(eq_2_B)
                allocate(X_B)
                ierr = reconstruct_PB3D_grid(grid_eq_B,'eq_B',&
                    &rich_lvl=rich_lvl,grid_limits=eq_limits)
                CHCKERR('')
                ierr = reconstruct_PB3D_grid(grid_X_B,'X_B',rich_lvl=rich_lvl,&
                    &grid_limits=X_limits)
                CHCKERR('')
                call create_eq(grid_eq_B,eq_2_B)
                call create_X(grid_X_B,X_B)
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,eq_2=eq_2,&
                    &eq_2_out=eq_2_B,eq_1=eq_1,&
                    &grid_name='field-aligned equilibrium grid')
                CHCKERR('')
                ierr = interp_HEL_on_grid(grid_X,grid_X_B,X_1=X,&
                    &X_1_out=X_B,grid_name='field-aligned perturbation grid')
                CHCKERR('')
        end select
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! user output
        call writo('Various PB3D output plots')
        call lvl_ud(1)
        
        if (plot_resonance) then
            ierr = resonance_plot(eq_1,grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        if (plot_flux_q) then
            ierr = flux_q_plot(grid_eq,eq_1)
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
        
        ! user output
        call lvl_ud(-1)
        call writo('PB3D outputs plots done')
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call writo('Extend PB3D output to plot grid')
        call lvl_ud(1)
        
        call writo('Setting up the grids')
        call lvl_ud(1)
        ierr = extend_grid_E(grid_eq,grid_eq_plot,grid_eq=grid_eq)              ! extend eq grid and convert to F
        CHCKERR('')
        ierr = extend_grid_E(grid_X,grid_X_plot,grid_eq=grid_eq)                ! extend X grid and convert to F
        CHCKERR('')
        ierr = extend_grid_E(grid_sol,grid_sol_plot,grid_eq=grid_eq)            ! extend sol grid and convert to F
        CHCKERR('')
        call lvl_ud(-1)
        call writo('Grids set up')
        
        ! calculating variables on plot grid depends on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                call writo('Calculating quantities on plot grid')
                call lvl_ud(1)
                
                ! back up no_plots and no_output and set to .true.
                no_plots_loc = no_plots; no_plots = .true.
                no_output_loc = no_output; no_output = .true.
                
                ! Calculate the flux equilibrium quantities
                ! Note: Only F quantities are saved, but E are needed here
                ierr = calc_eq(grid_eq_plot,eq_1_plot)
                CHCKERR('')
                
                ! Calculate the metric equilibrium quantitities
                ierr = calc_eq(grid_eq_plot,eq_1_plot,eq_2_plot)
                CHCKERR('')
                
                ! Transform E into F derivatives
                ierr = calc_F_derivs(eq_2_plot)
                CHCKERR('')
                
                ! Calculate derived metric quantities
                call calc_derived_q(grid_eq_plot,eq_1_plot,eq_2_plot)
                
                ! calculate X variables, vector phase
                ierr = calc_X(grid_eq_plot,grid_X_plot,eq_1_plot,eq_2_plot,&
                    &X_plot)
                CHCKERR('')
                
                ! no need any more for eq_1_plot E vars so eq_1 will do.
                call dealloc_eq(eq_1_plot)
                
                ! reset no_plots and no_output
                no_plots = no_plots_loc
                no_output = no_output_loc
                
                call lvl_ud(-1)
                call writo('Quantities calculated')
            case (2)                                                            ! HELENA
                call writo('Interpolating quantities on plot grid')
                call lvl_ud(1)
                
                ! setup plot quantities
                call create_eq(grid_eq_plot,eq_2_plot)
                call create_X(grid_X_plot,X_plot)
                
                ! call HELENA grid interpolation
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_plot,eq_2=eq_2,&
                    &eq_2_out=eq_2_plot,eq_1=eq_1,&
                    &grid_name='equilibrium plot grid')
                CHCKERR('')
                ierr = interp_HEL_on_grid(grid_X,grid_X_plot,X_1=X,&
                    &X_1_out=X_plot,grid_name='perturbation plot grid')
                CHCKERR('')
                
                call lvl_ud(-1)
                call writo('Quantities interpolated')
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
        call writo('Plot the Eigenvalues')
        call lvl_ud(1)
        
        ierr = plot_sol_vals(sol,last_unstable_id)
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Eigenvalues plotted')
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call writo('Prepare Eigenvector plots')
        call lvl_ud(1)
        
        call writo('Calculate helper variables')
        call lvl_ud(1)
        allocate(XYZ_plot(grid_sol_plot%n(1),grid_sol_plot%n(2),&
            &grid_sol_plot%loc_n_r,3))
        ierr = calc_XYZ_grid(grid_eq,grid_sol_plot,XYZ_plot(:,:,:,1),&
            &XYZ_plot(:,:,:,2),XYZ_plot(:,:,:,3))
        CHCKERR('')
        call lvl_ud(-1)
        
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        ierr = calc_res_surf(grid_eq,eq_1,res_surf,info=.false.)
        call lvl_ud(-1)
        
        call writo('Open decomposition log file')
        call lvl_ud(1)
        if (rank.eq.0) then
            ! set format strings
            format_head = '("#  ",A23," ",A23," ",A23," ",A23," ",A23," ",A23)'
            ! open output file for the log
            full_output_name = prog_name//'_'//trim(output_name)//'_EN.txt'
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
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call writo('Plotting variables for different ranges')
        call lvl_ud(1)
        
        ! loop over three ranges
        allocate(sol_val_comp(2,2,0))
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
                ierr = plot_sol_vec(grid_eq_plot,grid_X_plot,grid_sol_plot,&
                    &eq_1,eq_2_plot,X_plot,sol,XYZ_plot,id,res_surf)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! synchronize processes
                ierr = wait_MPI()
                CHCKERR('')
                
                ! user output
                call writo('Decompose the energy into its terms')
                call lvl_ud(1)
                allocate(sol_val_comp_loc(2,2,size(sol_val_comp,3)+1))
                sol_val_comp_loc(:,:,1:size(sol_val_comp,3)) = sol_val_comp
                ierr = decompose_energy(grid_eq_B,grid_X_B,grid_sol,eq_1,&
                    &eq_2_B,X_B,sol,id,output_EN_i,.true.,sol_val_comp=&
                    &sol_val_comp_loc(:,:,size(sol_val_comp_loc,3)))
                CHCKERR('')
                deallocate(sol_val_comp)
                allocate(sol_val_comp(2,2,size(sol_val_comp_loc,3)))
                sol_val_comp = sol_val_comp_loc
                deallocate(sol_val_comp_loc)
                ierr = decompose_energy(grid_eq_plot,grid_X_plot,grid_sol,&
                    &eq_1,eq_2_plot,X_plot,sol,id,output_EN_i,.false.,&
                    &XYZ=XYZ_plot)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! synchronize processes
                ierr = wait_MPI()
                CHCKERR('')
                
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
        call plot_sol_val_comp(sol_val_comp)
        
        ! user output
        call lvl_ud(-1)
        call writo('Variables for different ranges plotted')
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call dealloc_in()
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_grid(grid_sol)
        call dealloc_eq(eq_1)
        call dealloc_eq(eq_2)
        call dealloc_X(X)
        call dealloc_sol(sol)
        if (eq_style.eq.2) then
            call dealloc_grid(grid_eq_B)
            call dealloc_grid(grid_X_B)
            call dealloc_eq(eq_2_B)
            call dealloc_X(X_B)
            deallocate(grid_eq_B)
            deallocate(grid_X_B)
            deallocate(eq_2_B)
            deallocate(X_B)
        end if
        nullify(grid_eq_B)
        nullify(grid_X_B)
        nullify(eq_2_B)
        nullify(X_B)
        call dealloc_grid(grid_eq_plot)
        call dealloc_grid(grid_sol_plot)
        call dealloc_eq(eq_2_plot)
        call dealloc_X(X_plot)
        call lvl_ud(-1)
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
        
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
    subroutine plot_sol_val_comp(sol_val_comp)
        use num_vars, only: rank
        
        ! input /output
        complex(dp) :: sol_val_comp(:,:,:)                                      ! fraction between total E_pot and E_kin, compared with EV
        
        ! local variables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        
        if (rank.eq.0) then
            ! real part
            plot_title = 'RE sol_val and E_frac'
            plot_name = 'sol_val_comp_RE'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(sol_val_comp(:,1,:)),&
                &x=realpart(sol_val_comp(:,2,:)),draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,size(sol_val_comp,3),1,&
                &.false.)
            plot_title = 'RE sol_val and E_frac rel diff'
            plot_name = 'sol_val_comp_RE_rel_diff'
            call print_GP_2D(plot_title,plot_name,realpart(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)),x=realpart(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            plot_title = 'RE sol_val and E_frac log rel diff'
            plot_name = 'sol_val_comp_RE_rel_diff_log'
            call print_GP_2D(plot_title,plot_name,log10(abs(realpart(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)))),x=realpart(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            
            ! imaginary part
            plot_title = 'IM sol_val and E_frac'
            plot_name = 'sol_val_comp_IM'
            call print_GP_2D(plot_title,plot_name,&
                &realpart(sol_val_comp(:,1,:)),x=imagpart(sol_val_comp(:,2,:)),&
                &draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,size(sol_val_comp,3),1,&
                &.false.)
            plot_title = 'IM sol_val and E_frac rel diff'
            plot_name = 'sol_val_comp_IM_rel_diff'
            call print_GP_2D(plot_title,plot_name,imagpart(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)),x=realpart(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
            plot_title = 'IM sol_val and E_frac log rel diff'
            plot_name = 'sol_val_comp_IM_rel_diff_log'
            call print_GP_2D(plot_title,plot_name,log10(abs(imagpart(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)))),x=realpart(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_GP(plot_title,plot_name,plot_name,1,1,.false.)
        end if
    end subroutine plot_sol_val_comp
end module driver_POST
