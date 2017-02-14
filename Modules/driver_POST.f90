!------------------------------------------------------------------------------!
!   Main driver of PostProcessing of program Peeling Ballooning in 3D          !
!------------------------------------------------------------------------------!
module driver_POST
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
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
    ! Furthermore,  if Richardson  extrapolation  is used,  the  VMEC grids  and
    ! variables are contained in multiple HDF5 variables.
    ! These variables are needed here both in a field-aligned and a plot grid. A
    ! consequence of  the above  is that  for VMEC  the field-aligned  output is
    ! already given, but the  output on the plot grid has  to be calculated from
    ! scratch, while  for HELENA both outputs  have to be interpolated  from the
    ! output tables.
    integer function run_driver_POST() result(ierr)
        use num_vars, only: no_output, no_plots, eq_style, plot_resonance, &
            &plot_flux_q, plot_magn_grid, rank, norm_disc_prec_X, POST_style, &
            &slab_plots_style, use_pol_flux_F, pi, swap_angles, minim_output, &
            &decomp_i
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_1, reconstruct_PB3D_sol
        use grid_vars, only: disc_type
        use grid_ops, only: calc_norm_range, magn_grid_plot
        use eq_vars, only: max_flux_F
        use eq_ops, only: calc_eq, flux_q_plot, calc_derived_q
        use eq_utilities, only: calc_F_derivs
        use sol_vars, only: alpha
        use grid_utilities, only: calc_XYZ_grid, extend_grid_E, &
            &setup_interp_data, apply_disc
        use X_ops, only: calc_X, resonance_plot, calc_res_surf, setup_nm_X
        use sol_ops, only: plot_sol_vals, plot_harmonics, plot_sol_vec, &
            &decompose_energy
        use HELENA_ops, only: interp_HEL_on_grid
        use VMEC, onLy: calc_trigon_factors
        use input_utilities, only: dealloc_in
        use num_utilities, only: calc_aux_utilities
        use MPI_utilities, only: wait_MPI
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'run_driver_POST'
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_out                                 ! output equilibrium grid
        type(grid_type), target :: grid_X                                       ! perturbation grid
        type(grid_type), pointer :: grid_X_out                                  ! output perturbation grid
        type(grid_type) :: grid_sol                                             ! solution grid
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium
        type(eq_1_type) :: eq_1_out                                             ! output flux equilibrium
        type(eq_2_type), target :: eq_2                                         ! metric equilibrium
        type(eq_2_type), pointer :: eq_2_out                                    ! output metric equilibrium
        type(X_1_type), target :: X                                             ! vectorial perturbation variables
        type(X_1_type), pointer :: X_out                                        ! output vectorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        integer :: id, jd, kd                                                   ! counter
        integer :: min_id(3), max_id(3)                                         ! min. and max. index of range 1, 2 and 3
        integer :: last_unstable_id                                             ! index of last unstable EV
        integer :: eq_limits(2)                                                 ! i_limit of eq variables
        integer :: X_limits(2)                                                  ! i_limit of X variables
        integer :: sol_limits(2)                                                ! i_limit of sol variables
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_output_loc                                                ! local copy of no_output
        logical :: B_aligned                                                    ! whether the grid is field-aligned
        logical :: full_output                                                  ! whether full output is possible
        real(dp), allocatable :: XYZ_out(:,:,:,:)                               ! X, Y and Z on output grid
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        complex(dp), allocatable :: sol_val_comp(:,:,:)                         ! fraction between total E_pot and E_kin, compared with EV
        complex(dp), allocatable :: sol_val_comp_loc(:,:,:)                     ! local sol_val_comp
        character(len=max_str_ln) :: err_msg                                    ! error message
        
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
        
        ! set up whether full output is possible
        full_output = .true.
        if (minim_output .and. POST_style.eq.2) full_output = .false.           ! field-aligned quantities not saved
        if (.not.full_output) call writo('This is a minimal PB3D output: &
            &field-aligned output is not available',alert=.true.)
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('Setting up grid limits')
        call lvl_ud(1)
        
        ! reconstruct full equilibrium, perturbation and solution grids and flux
        ! equilibrium quantities
        ! note: only the normal part is needed, so no need for tot_rich.
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_sol,'sol')
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
        CHCKERR('')
        
        ! set up nm in full grids
        !  Note: This  is needed  as many  routines count  on this,  for example
        ! 'set_nm_X' in X_vars, also used to initialize a solution variable.
        ierr = setup_nm_X(grid_eq,grid_X,eq_1,plot_nm=.true.)                   ! is necessary for X variables
        CHCKERR('')
        
        ! set eq and X limits, using r_F of the grids
        ierr = calc_norm_range(eq_limits=eq_limits,X_limits=X_limits,&
            &sol_limits=sol_limits,r_F_eq=grid_eq%r_F,r_F_X=grid_X%r_F,&
            &r_F_sol=grid_sol%r_F)
        CHCKERR('')
        call writo('grid limits:')
        call lvl_ud(1)
        call writo('proc '//trim(i2str(rank))//' - equilibrium:  '//&
            &trim(i2str(eq_limits(1)))//' .. '//trim(i2str(eq_limits(2))),&
            &persistent=.true.)
        call writo('proc '//trim(i2str(rank))//' - perturbation: '//&
            &trim(i2str(X_limits(1)))//' .. '//trim(i2str(X_limits(2))),&
            &persistent=.true.)
        call writo('proc '//trim(i2str(rank))//' - solution:     '//&
            &trim(i2str(sol_limits(1)))//' .. '//trim(i2str(sol_limits(2))),&
            &persistent=.true.)
        call lvl_ud(-1)
        
        ! deallocate the grids and equilibrium flux quantities
        call grid_eq%dealloc()
        call grid_X%dealloc()
        call grid_sol%dealloc()
        call eq_1%dealloc()
        
        ! user output
        call lvl_ud(-1)
        call writo('Limits set up')
        
        ! user output
        call writo('Reconstructing PB3D output')
        call lvl_ud(1)
        
        ! reconstruct local equilibrium, perturbation and solution grids
        call writo('Reconstructing original PB3D output')
        call lvl_ud(1)
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.,grid_limits=eq_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.,grid_limits=X_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_sol,'sol',grid_limits=sol_limits)
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
        CHCKERR('')
        if (.not.minim_output) then
            ierr = reconstruct_PB3D_eq_2(grid_eq,eq_2,'eq_2',&
                &rich_lvl=rich_lvl_name,tot_rich=.true.)
            CHCKERR('')
            ierr = reconstruct_PB3D_X_1(grid_X,X,'X_1',rich_lvl=rich_lvl_name,&
                &tot_rich=.true.)
            CHCKERR('')
        end if
        ierr = reconstruct_PB3D_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl)
        CHCKERR('')
        call lvl_ud(-1)
        call writo('Original PB3D output reconstructed')
        
        call writo('Set up output grid')
        call lvl_ud(1)
        
        ! set up  output grids grid_eq_out and grid_X_out, depending on POST
        ! style
        select case (POST_style)
            case (1)                                                            ! extended grid
                ! user output
                call writo('POST style 1: Output grid is extended grid')
                
                ! set B-aligned
                B_aligned = .false.
                
                ! allocate grids
                allocate(grid_eq_out)
                allocate(grid_X_out)
                
                ! extend eq grid
                ierr = extend_grid_E(grid_eq,grid_eq_out,grid_eq=grid_eq)       ! extend eq grid and convert to F
                CHCKERR('')
                
                ! interpolate extended eq grid to extended X grid
                call writo('Interpolate to extended perturbation grid')
                call lvl_ud(1)
                ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
                    &norm_interp_data,norm_disc_prec_X)
                CHCKERR('')
                ierr = grid_X_out%init([grid_eq_out%n(1:2),grid_X%n(3)],&
                    &i_lim=X_limits)
                CHCKERR('')
                grid_X_out%r_E = grid_X%r_E
                grid_X_out%r_F = grid_X%r_F
                grid_X_out%loc_r_E = grid_X%loc_r_E
                grid_X_out%loc_r_F = grid_X%loc_r_F
                ierr = apply_disc(grid_eq_out%theta_E,norm_interp_data,&
                    &grid_X_out%theta_E,3)
                CHCKERR('')
                ierr = apply_disc(grid_eq_out%theta_F,norm_interp_data,&
                    &grid_X_out%theta_F,3)
                CHCKERR('')
                ierr = apply_disc(grid_eq_out%zeta_E,norm_interp_data,&
                    &grid_X_out%zeta_E,3)
                CHCKERR('')
                ierr = apply_disc(grid_eq_out%zeta_F,norm_interp_data,&
                    &grid_X_out%zeta_F,3)
                CHCKERR('')
                call norm_interp_data%dealloc()                                 ! clean up
                call lvl_ud(-1)
            case (2)                                                            ! field-aligned grid
                ! user output
                call writo('POST style 2: Output grid is field-aligned &
                    &grid')
                
                ! set B-aligned
                B_aligned = .true.
                
                ! depends on equilibrium style
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! user output
                        call writo('The grids are already field-aligned')
                        
                        ! the grids are already field-aligned
                        grid_eq_out => grid_eq
                        grid_X_out => grid_X
                    case (2)                                                    ! HELENA
                        ! user output
                        call writo('For HELENA, the grids are read from &
                            &the output file')
                        
                        ! allocate grids
                        allocate(grid_eq_out)
                        allocate(grid_X_out)
                        
                        ! reconstruct from output
                        ierr = reconstruct_PB3D_grid(grid_eq_out,'eq_B',&
                            &rich_lvl=rich_lvl,tot_rich=.true.,&
                            &grid_limits=eq_limits)
                        CHCKERR('')
                        ierr = reconstruct_PB3D_grid(grid_X_out,'X_B',&
                            &rich_lvl=rich_lvl,tot_rich=.true.,&
                            &grid_limits=X_limits)
                        CHCKERR('')
                end select
        end select
        
        call lvl_ud(-1)
        call writo('Output grid set up')
        
        if (full_output) then
            call writo('Set up output variables')
            call lvl_ud(1)
            
            ! set up output eq_2 and X_1, depending on equilibrium style
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! reconstruct variables depending on POST style
                    select case (POST_style)
                        case (1)                                                ! extended grid
                            call writo('Recalculating variables on extended &
                                &grid')
                            call lvl_ud(1)
                            
                            ! allocate variables
                            allocate(eq_2_out)
                            allocate(X_out)
                            
                            ! back up no_plots and no_output and set to .true.
                            no_plots_loc = no_plots; no_plots = .true.
                            no_output_loc = no_output; no_output = .true.
                            
                            ! Calculate the flux equilibrium quantities
                            ! Note:  Only  F  quantities  are saved,  but E  are
                            ! needed temporarily now
                            ierr = calc_eq(grid_eq_out,eq_1_out)
                            CHCKERR('')
                            
                            ! Calculate the metric equilibrium quantitities
                            ierr = calc_eq(grid_eq_out,eq_1_out,eq_2_out)
                            CHCKERR('')
                            
                            ! Transform E into F derivatives
                            ierr = calc_F_derivs(eq_2_out)
                            CHCKERR('')
                            
                            ! Calculate derived metric quantities
                            call calc_derived_q(grid_eq_out,eq_1_out,eq_2_out)
                            
                            ! calculate X variables, vector phase
                            ierr = calc_X(grid_eq_out,grid_X_out,eq_1_out,&
                                &eq_2_out,X_out)
                            CHCKERR('')
                            
                            ! no need any more for eq_1_out E vars so eq_1 is ok
                            call eq_1_out%dealloc()
                            
                            ! reset no_plots and no_output
                            no_plots = no_plots_loc
                            no_output = no_output_loc
                            
                            call lvl_ud(-1)
                            call writo('Quantities recalculated')
                        case (2)                                                ! field-aligned grid
                            ! user output
                            call writo('The grids are already field-aligned')
                            
                            ! the grids are already field-aligned
                            eq_2_out => eq_2
                            X_out => X
                    end select
                case (2)                                                        ! HELENA
                    call writo('Interpolate the variables on the output grid')
                    call lvl_ud(1)
                    
                    ! allocate variables
                    allocate(eq_2_out)
                    allocate(X_out)
                    
                    ! interpolate
                    call eq_2_out%init(grid_eq_out)
                    call X_out%init(grid_X_out)
                    ierr = interp_HEL_on_grid(grid_eq,grid_eq_out,eq_2=eq_2,&
                        &eq_2_out=eq_2_out,eq_1=eq_1,&
                        &grid_name='output equilibrium grid')
                    CHCKERR('')
                    ierr = interp_HEL_on_grid(grid_X,grid_X_out,X_1=X,&
                        &X_1_out=X_out,grid_name='output perturbation grid')
                    CHCKERR('')
                    call lvl_ud(-1)
            end select
            
            call lvl_ud(-1)
            call writo('Output variables set up')
        end if
        
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
        
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
        if (plot_magn_grid) then
            if (POST_style.eq.2) then                                           ! output grid is field-aligned
                ierr = magn_grid_plot(grid_eq_out)
                CHCKERR('')
            else
                call writo('Need POST style 2 (field-aligned) for grid plot',&
                    &warning=.true.)
            end if
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
        
        if (full_output) then
            call writo('Calculate helper variables')
            call lvl_ud(1)
            
            ! if VMEC, calculate trigonometric factors of output grid
            if (eq_style.eq.1) then
                ierr = calc_trigon_factors(grid_X_out%theta_E,&
                    &grid_X_out%zeta_E,grid_X_out%trigon_factors)
                CHCKERR('')
            end if
            
            ! set up XYZ
            allocate(XYZ_out(grid_X_out%n(1),grid_X_out%n(2),&
                &grid_X_out%loc_n_r,3))
            select case (slab_plots_style)
                case (0)                                                        ! 3-D geometry
                    ierr = calc_XYZ_grid(grid_eq,grid_X_out,XYZ_out(:,:,:,1),&
                        &XYZ_out(:,:,:,2),XYZ_out(:,:,:,3))
                    CHCKERR('')
                    !!! To plot the cross-section
                    !!call print_ex_2D(['cross_section'],'cross_section',&
                        !!&XYZ_out(:,1,:,3),x=XYZ_out(:,1,:,1),draw=.false.)
                    !!call draw_ex('cross_section','cross_section',&
                        !!&'cross_section',size(XYZ_out,3),1,.false.)
                case (1,2)                                                      ! slab geometry (with or without wrapping to fundamental interval)
                    ! user output
                    call writo('Plots are done in slab geometry')
                    call lvl_ud(1)
                    
                    if (B_aligned .and. slab_plots_style.eq.1) then             ! for slab_plots style 2 both theta and zeta need to be chosen, not alpha
                        if (swap_angles) then
                            call writo('For plotting, the angular coordinates &
                                &are swapped: theta <-> zeta')
                            if (use_pol_flux_F) then
                                XYZ_out(:,:,:,1) = grid_X_out%zeta_F/pi
                            else
                                XYZ_out(:,:,:,1) = grid_X_out%theta_F/pi
                            end if
                        else
                            if (use_pol_flux_F) then
                                XYZ_out(:,:,:,1) = grid_X_out%theta_F/pi
                            else
                                XYZ_out(:,:,:,1) = grid_X_out%zeta_F/pi
                            end if
                        end if
                        XYZ_out(:,:,:,2) = alpha
                    else
                        if (use_pol_flux_F) then
                            XYZ_out(:,:,:,1) = grid_X_out%theta_F/pi
                            XYZ_out(:,:,:,2) = grid_X_out%zeta_F/pi
                        else
                            XYZ_out(:,:,:,1) = grid_X_out%zeta_F/pi
                            XYZ_out(:,:,:,2) = grid_X_out%theta_F/pi
                        end if
                    end if
                    
                    do kd = 1,grid_X_out%loc_n_r
                        XYZ_out(:,:,kd,3) = grid_X_out%loc_r_F(kd)/&
                            &max_flux_F*2*pi
                    end do
                    
                    ! limit to fundamental interval -1..1
                    if (slab_plots_style.eq.2) XYZ_out(:,:,:,1:2) = &
                        &modulo(XYZ_out(:,:,:,1:2)+1._dp,2._dp)-1._dp
                    
                    call lvl_ud(-1)
                case default
                    err_msg = 'No style "'//trim(i2str(slab_plots_style))//&
                        &'" for slab_plots_style'
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! deallocate memory-thirsty trigonometric factors
            if (eq_style.eq.1) deallocate(grid_X_out%trigon_factors)
            
            call lvl_ud(-1)
            
            ! open decomposition log file
            ierr = open_decomp_log()
            CHCKERR('')
        end if
        
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        ierr = calc_res_surf(grid_eq,eq_1,res_surf,info=.false.)
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
                
                call writo('Plot the Harmonics')
                call lvl_ud(1)
                ierr = plot_harmonics(grid_sol,sol,id,res_surf)
                CHCKERR('')
                call lvl_ud(-1)
                
                ! synchronize processes
                ierr = wait_MPI()
                CHCKERR('')
                
                if (full_output) then
                    ! user output
                    call writo('Plot the Eigenvector')
                    call lvl_ud(1)
                    ierr = plot_sol_vec(grid_eq_out,grid_X_out,grid_sol,eq_1,&
                        &eq_2_out,X_out,sol,XYZ_out,id)
                    CHCKERR('')
                    call lvl_ud(-1)
                    
                    ! user output
                    call writo('Decompose the energy into its terms')
                    call lvl_ud(1)
                    allocate(sol_val_comp_loc(2,2,size(sol_val_comp,3)+1))
                    sol_val_comp_loc(:,:,1:size(sol_val_comp,3)) = sol_val_comp
                    ierr = decompose_energy(grid_eq_out,grid_X_out,grid_sol,&
                        &eq_1,eq_2_out,X_out,sol,id,B_aligned,&
                        &log_i=decomp_i,sol_val_comp=&
                        &sol_val_comp_loc(:,:,size(sol_val_comp_loc,3)),&
                        &XYZ=XYZ_out)
                    CHCKERR('')
                    deallocate(sol_val_comp)
                    allocate(sol_val_comp(2,2,size(sol_val_comp_loc,3)))
                    sol_val_comp = sol_val_comp_loc
                    deallocate(sol_val_comp_loc)
                    call lvl_ud(-1)
                    
                    ! synchronize processes
                    ierr = wait_MPI()
                    CHCKERR('')
                end if
                
                ! user output
                call lvl_ud(-1)
                call writo('Mode '//trim(i2str(id))//'/'//&
                    &trim(i2str(size(sol%val)))//' finished')
            end do
            
            call lvl_ud(-1)
            if (min_id(jd).le.max_id(jd)) call &
                &writo('RANGE '//trim(i2str(jd))//' plotted')
        end do
        
        if (full_output) then
            ! print difference between Eigenvalue and energy fraction
            call plot_sol_val_comp(sol_val_comp)
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Variables for different ranges plotted')
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call dealloc_in()
        call grid_eq%dealloc()
        call grid_X%dealloc()
        call grid_sol%dealloc()
        call eq_1%dealloc()
        if (.not.minim_output) then
            call eq_2%dealloc()
            call X%dealloc()
        end if
        call sol%dealloc()
        select case (POST_style)
            case (1)                                                            ! extended grid
                call grid_eq_out%dealloc()
                call grid_X_out%dealloc()
                deallocate(grid_eq_out)
                deallocate(grid_X_out)
                if (full_output) then
                    call eq_2_out%dealloc()
                    call X_out%dealloc()
                    deallocate(eq_2_out)
                    deallocate(X_out)
                end if
            case (2)                                                            ! field-aligned grid
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! do nothing
                    case (2)                                                    ! HELENA
                        call grid_eq_out%dealloc()
                        call grid_X_out%dealloc()
                        deallocate(grid_eq_out)
                        deallocate(grid_X_out)
                        if (full_output) then
                            call eq_2_out%dealloc()
                            call X_out%dealloc()
                            deallocate(eq_2_out)
                            deallocate(X_out)
                        end if
                end select
        end select
        nullify(grid_eq_out)
        nullify(grid_X_out)
        if (.not.minim_output) then
            nullify(eq_2_out)
            nullify(X_out)
        end if
        call lvl_ud(-1)
        
        ! synchronize processes
        ierr = wait_MPI()
        CHCKERR('')
        
        ! close output
        if (rank.eq.0) close(decomp_i)
    end function run_driver_POST
    
    ! opens the decomposition log file
    integer function open_decomp_log() result(ierr)
        use num_vars, only: prog_name, output_name, rank, POST_style, &
            &decomp_i
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'open_decomp_log'
        
        ! local variables
        character(len=max_str_ln) :: format_head                                ! header format
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: grid_name                                  ! name of grid on which calculations are done
        character(len=2*max_str_ln) :: temp_output_str                          ! temporary output string
        
        ! initialize ierr
        ierr = 0
        
        ! set up output name
        select case (POST_style)
            case (1)                                                            ! extended grid
                grid_name = 'extended'
            case (2)                                                            ! field-aligned grid
                grid_name = 'field-aligned'
        end select
        
        call writo('Open decomposition log file')
        call lvl_ud(1)
        if (rank.eq.0) then
            ! set format strings
            format_head = '("#  ",A23," ",A23," ",A23," ",A23," ",A23," ",&
                &A23)'
            ! open output file for the log
            full_output_name = prog_name//'_'//trim(output_name)//'_EN_R'//&
                &trim(i2str(rich_lvl))//'.txt'
            open(UNIT=decomp_i,STATUS='replace',FILE=full_output_name,&
                &IOSTAT=ierr)
            CHCKERR('Cannot open EN output file')
            call writo('Log file opened in '//trim(full_output_name))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &'# Energy decomposition using the solution Eigenvectors'
            CHCKERR('Failed to write')
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &'# (Output on the '//trim(grid_name)//' grid)'
            CHCKERR('Failed to write')
            write(temp_output_str,format_head) &
                &'RE Eigenvalue          ', 'RE E_pot/E_kin         ', &
                &'RE E_pot               ', 'RE E_kin               '
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_head) &
                &'IM Eigenvalue          ', 'IM E_pot/E_kin         ', &
                &'IM E_pot               ', 'IM E_kin               '
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_head) &
                &'RE E_kin_n             ', 'RE E_kin_g             '
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_head) &
                &'IM E_kin_n             ', 'IM E_kin_g             '
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_head) &
                &'RE E_pot line_bending_n', 'RE E_pot line_bending_g', &
                &'RE E_pot ballooning_n  ', 'RE E_pot ballooning_g  ', &
                &'RE E_pot kink_n        ', 'RE E_pot kink_g        '
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_head) &
                &'IM E_pot line_bending_n', 'IM E_pot line_bending_g', &
                &'IM E_pot ballooning_n  ', 'IM E_pot ballooning_g  ', &
                &'IM E_pot kink_n        ', 'IM E_pot kink_g        '
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
        end if
        call lvl_ud(-1)
    end function open_decomp_log
    
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
            if (rp(sol%val(id)).lt.0._dp) last_unstable_id = id
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
        character(len=max_str_ln) :: plot_title(2)                              ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        
        if (rank.eq.0) then
            ! real part
            plot_title = ['RE sol_val','E_frac    ']
            plot_name = 'sol_val_comp_RE'
            call print_ex_2D(plot_title,plot_name,&
                &rp(sol_val_comp(:,1,:)),&
                &x=rp(sol_val_comp(:,2,:)),draw=.false.)
            call draw_ex(plot_title,plot_name,size(sol_val_comp,3),1,.false.)
            plot_title(1) = 'rel diff between RE sol_val and E_frac'
            plot_name = 'sol_val_comp_RE_rel_diff'
            call print_ex_2D(plot_title(1),plot_name,rp(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)),x=rp(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_ex(plot_title(1:1),plot_name,1,1,.false.)
            plot_title(1) = 'rel diff between RE sol_val and E_frac [log]'
            plot_name = 'sol_val_comp_RE_rel_diff_log'
            call print_ex_2D(plot_title(1),plot_name,log10(abs(rp(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)))),x=rp(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_ex(plot_title(1:1),plot_name,1,1,.false.)
            
            ! imaginary part
            plot_title = ['IM sol_val','E_frac    ']
            plot_name = 'sol_val_comp_IM'
            call print_ex_2D(plot_title,plot_name,&
                &rp(sol_val_comp(:,1,:)),x=ip(sol_val_comp(:,2,:)),draw=.false.)
            call draw_ex(plot_title,plot_name,size(sol_val_comp,3),1,.false.)
            plot_title(1) = 'rel diff between IM sol_val and E_frac'
            plot_name = 'sol_val_comp_IM_rel_diff'
            call print_ex_2D(plot_title(1),plot_name,ip(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)),x=rp(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_ex(plot_title(1:1),plot_name,1,1,.false.)
            plot_title(1) = 'rel diff between IM sol_val and E_frac [log]'
            plot_name = 'sol_val_comp_IM_rel_diff_log'
            call print_ex_2D(plot_title(1),plot_name,log10(abs(ip(&
                &(sol_val_comp(1,2,:)-sol_val_comp(2,2,:))/&
                &sol_val_comp(1,2,:)))),x=rp(sol_val_comp(1,1,:)),&
                &draw=.false.)
            call draw_ex(plot_title(1:1),plot_name,1,1,.false.)
        end if
    end subroutine plot_sol_val_comp
end module driver_POST
