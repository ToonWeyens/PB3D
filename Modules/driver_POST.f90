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
    use rich_vars, only: rich_lvl
    
    implicit none
    private
    public run_driver_POST, init_POST, stop_POST
    
    ! local variables
    integer :: lims_norm(2,3)                                                   ! normal i_limit of eq, X and sol variables
    integer :: rich_lvl_name                                                    ! either the Richardson level or zero, to append to names
    integer :: min_id(3), max_id(3)                                             ! min. and max. index of range 1, 2 and 3
    integer :: last_unstable_id                                                 ! index of last unstable EV
    integer :: n_par_tot                                                        ! total n(1) for all equilibrium jobs together
    character(len=4) :: grid_name(3)                                            ! names of grids
    type(sol_type) :: sol                                                       ! solution variables
    type(grid_type) :: grid_eq_XYZ                                              ! eq grid for XYZ calculations
    type(grid_type) :: grid_eq_HEL                                              ! HELENA eq grid
    type(grid_type) :: grid_X_HEL                                               ! HELENA X grid
    type(eq_2_type) :: eq_2_HEL                                                 ! metric equilibrium for HELENA tables
    type(X_1_type) :: X_HEL                                                     ! vectorial perturbation variables for HELENA tables
    complex(dp), allocatable :: E_pot_int(:,:)                                  ! E_pot integrated for requested solutions
    complex(dp), allocatable :: E_kin_int(:,:)                                  ! E_kin integrated for requested solutions

contains
    ! Initializes the POST driver:
    !   - set up preliminary variables
    !       * global variables
    !       * read grids (full)
    !       * eq_1 (full) and n and m (full)
    !   - set up output grids:
    !       * POST_style = 1: extended grid
    !       * POST_style = 2: field-aligned grid
    !   - clean up
    !   - set up final variables
    !       * normal limits
    !       * read grids (divided), eq_1 (divided) and sol (divided) variables
    !   - 1-D output plots
    !       * resonance plot
    !       * flux quantities plot
    !       * magnetic grid plot
    !   - prepare Eigenvalue plots
    !       * calculates resonant surfaces
    !       * plots Eigenvalues
    !       * finds stability ranges for all Eigenvalues to be plot
    !       * plots harmonics for every Eigenvalue
    !       * calculates the parallel ranges of the equilibrium jobs
    !   - clean up
    ! In the  actual driver, more detailed  plots are possibly made  for all the
    ! requested Eigenvalues if full_output.
    integer function init_POST() result(ierr)
        use grid_vars, only: disc_type
        use num_vars, only: POST_style, eq_style, rank, plot_magn_grid, &
            &plot_resonance, plot_flux_q, eq_jobs_lims, plot_grid_style, &
            &n_theta_plot, n_zeta_plot, POST_output_full, POST_output_sol
        use eq_ops, only: flux_q_plot, divide_eq_jobs, calc_eq_jobs_lims
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_1, reconstruct_PB3D_sol, get_PB3D_grid_size
        use X_ops, only: setup_nm_X, resonance_plot, calc_res_surf
        use grid_ops, only: calc_norm_range, magn_grid_plot
        use grid_utilities, only: extend_grid_E, setup_interp_data, &
            &apply_disc, copy_grid
        use HELENA_vars, only: nchi
        use sol_ops, only: plot_sol_vals, plot_harmonics
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'init_POST'
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B                                   ! field-aligned equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        type(grid_type) :: grid_sol                                             ! solution grid
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium
        integer :: n_div                                                        ! number of divisions of parallel points
        integer :: n_out(3,2)                                                   ! total grid sizes for eq and X output grids
        integer :: id, jd                                                       ! counters
        integer :: var_size_without_par(2)                                      ! size without parallel dimension for eq_2 and X_1 variables
        integer :: lim_loc(3,2)                                                 ! grid ranges for local equilibrium job
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Setting up preliminary variables')
        call lvl_ud(1)
        
        ! some preliminary things
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! set up whether Richardson level has to be appended to the name
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        ! set up grid names
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! take the standard grids
                grid_name(1) = 'eq'
                grid_name(2) = 'X'
                grid_name(3) = 'sol'
            case (2)                                                            ! HELENA
                ! take the field-aligned grids
                grid_name(1) = 'eq_B'
                grid_name(2) = 'X_B'
                grid_name(3) = 'sol'
        end select
        
        ! limits to exclude angular part
        lim_loc(1,:) = [0,0]
        lim_loc(2,:) = [0,0]
        lim_loc(3,:) = [1,-1]
        
        ! reconstruct normal part of full equilibrium, perturbation and solution
        ! grids
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.,lim_pos=lim_loc)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.,lim_pos=lim_loc)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_sol,'sol')
        CHCKERR('')
        
        ! set up nm in full grids
        !  Note: This  is needed  as many  routines count  on this,  for example
        ! 'set_nm_X' in X_vars, also used to initialize a solution variable.
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
        CHCKERR('')
        ierr = setup_nm_X(grid_eq,grid_X,eq_1,plot_nm=.true.)                   ! is necessary for X variables
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Preliminary variables set up')
        
        ! user output
        call writo('Dividing grids over MPI processes')
        call lvl_ud(1)
        
        ! set eq and X limits, using r_F of the grids
        ierr = calc_norm_range(eq_limits=lims_norm(:,1),&
            &X_limits=lims_norm(:,2),sol_limits=lims_norm(:,3),&
            &r_F_eq=grid_eq%r_F,r_F_X=grid_X%r_F,r_F_sol=grid_sol%r_F)
        CHCKERR('')
        call writo('normal grid limits:')
        call lvl_ud(1)
        call writo('proc '//trim(i2str(rank))//' - equilibrium:  '//&
            &trim(i2str(lims_norm(1,1)))//' .. '//trim(i2str(lims_norm(2,1))),&
            &persistent=.true.)
        call writo('proc '//trim(i2str(rank))//' - perturbation: '//&
            &trim(i2str(lims_norm(1,2)))//' .. '//trim(i2str(lims_norm(2,2))),&
            &persistent=.true.)
        call writo('proc '//trim(i2str(rank))//' - solution:     '//&
            &trim(i2str(lims_norm(1,3)))//' .. '//trim(i2str(lims_norm(2,3))),&
            &persistent=.true.)
        call lvl_ud(-1)
        
        ! synchronize outputs
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Grids divided over MPI processes')
        
        ! clean up
        call grid_eq%dealloc()
        call grid_X%dealloc()
        call grid_sol%dealloc()
        call eq_1%dealloc()
        
        ! user output
        call writo('Setting up final variables')
        call lvl_ud(1)
        
        ! reconstruct variables
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.,grid_limits=lims_norm(:,1))
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.,grid_limits=lims_norm(:,2))
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_sol,'sol',grid_limits=lims_norm(:,3))
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
        CHCKERR('')
        if (POST_output_sol) then
            ierr = reconstruct_PB3D_sol(grid_sol,sol,'sol',rich_lvl=rich_lvl)
            CHCKERR('')
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Final variables set up')
        
        ! user output
        call writo('1-D PB3D output plots')
        call lvl_ud(1)
        
        if (plot_resonance) then
            ierr = resonance_plot(grid_eq,eq_1)
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
            select case (eq_style)
                case (1)                                                        ! VMEC
                    grid_eq_B => grid_eq
                case (2)                                                        ! HELENA
                    allocate(grid_eq_B)
                    ierr = reconstruct_PB3D_grid(grid_eq_B,'eq_B',&
                        &rich_lvl=rich_lvl,tot_rich=.true.,&
                        &grid_limits=lims_norm(:,1))
                    CHCKERR('')
            end select
            ierr = magn_grid_plot(grid_eq_B)
            CHCKERR('')
            if (eq_style.eq.2) then
                call grid_eq_B%dealloc()
                deallocate(grid_eq_B)
            end if
            nullify(grid_eq_B)
        else
            call writo('Magnetic grid plot not requested')
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('1-D outputs plots done')
        
        ! user output
        call writo('Prepare 3-D plots')
        call lvl_ud(1)
        
        ! calculate resonant surfaces
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        ierr = calc_res_surf(grid_eq,eq_1,res_surf,info=.false.)
        call lvl_ud(-1)
        
        ! plot eigenvalues
        if (POST_output_sol) then
            call writo('Plot the Eigenvalues')
            call lvl_ud(1)
            call plot_sol_vals(sol,last_unstable_id)
            call lvl_ud(-1)
            
            ! find stability regions
            call writo('Find stability ranges')
            call lvl_ud(1)
            call find_stab_ranges(sol,min_id,max_id,last_unstable_id)
            call lvl_ud(-1)
            
            ! plots harmonics for every Eigenvalue, looping over three ranges
            call writo('Plot the Harmonics')
            call lvl_ud(1)
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
                    ierr = plot_harmonics(grid_sol,sol,id,res_surf)
                    CHCKERR('')
                    
                    ! user output
                    call lvl_ud(-1)
                    call writo('Mode '//trim(i2str(id))//'/'//&
                        &trim(i2str(size(sol%val)))//' finished')
                end do
                
                call lvl_ud(-1)
                if (min_id(jd).le.max_id(jd)) call &
                    &writo('RANGE '//trim(i2str(jd))//' plotted')
            end do
            call lvl_ud(-1)
        end if
        
        ! Divide   the  equilibrium   jobs  if   full  output   (see  subroutine
        ! "calc_memory"  inside  "divide_eq_jobs"),  only  necessary  when  full
        ! output. The results are saved in the global variables eq_jobs_lims for
        ! every eq_job.
        ! If full output  not requested, allocate eq_jobs_lims to  zero size, so
        ! that the driver is never run.
        if (POST_output_full) then
            ! set total grid sizes
            do id = 1,2
                ierr = get_PB3D_grid_size(n_out(:,id),grid_name(id),&
                    &rich_lvl=rich_lvl,tot_rich=.true.)
                CHCKERR('')
                if (POST_style.eq.1) n_out(1:2,id) = [n_theta_plot,n_zeta_plot]
            end do
            
            ! test if angular grid sizes are compatible
            if (n_out(1,1).ne.n_out(1,2) .or. n_out(2,1).ne.n_out(2,2)) then
                ierr = 1
                call writo('grid sizes for output grids not compatible. This &
                    &is assumed in the calculation of the equilibrium jobs.')
                call writo('If this has changed, adapt and generalize &
                    &"divide_eq_jobs" appropriately.')
                err_msg = 'Take into account how the ranges transfer from &
                    &grid to grid'
                CHCKERR(err_msg)
            end if
            
            ! set size of all variables, without parallel dimension
            var_size_without_par(1) = n_out(3,1)
            var_size_without_par(2) = n_out(3,2)
            
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! divide equilibrium jobs
                    ierr = divide_eq_jobs(product(n_out(1:2,1)),&
                        &var_size_without_par,n_div)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ! divide equilibrium jobs
                    ierr = divide_eq_jobs(product(n_out(1:2,1)),&
                        &var_size_without_par,n_div,n_par_X_base=nchi)
                    CHCKERR('')
            end select
            
            ! calculate equilibrium job limits
            ierr = calc_eq_jobs_lims(n_out(1,1),n_div)
            CHCKERR('')
            
            ! set up total number of parallel points
            n_par_tot = maxval(eq_jobs_lims(2,:))-minval(eq_jobs_lims(1,:))+1
        else
            allocate(eq_jobs_lims(2,0))
        end if
        
        ! Calculate variables  on HELENA output  grid if full output  and HELENA
        ! for later interpolation
        if (POST_output_full .and. eq_style.eq.2) then
            ! user output
            call writo('Reconstructing HELENA output for later interpolation')
            call lvl_ud(1)
            
            ierr = reconstruct_PB3D_grid(grid_eq_HEL,'eq',&
                &grid_limits=lims_norm(:,1))
            CHCKERR('')
            ierr = reconstruct_PB3D_grid(grid_X_HEL,'X',&
                &grid_limits=lims_norm(:,2))
            CHCKERR('')
            ierr = reconstruct_PB3D_eq_2(grid_eq_HEL,eq_2_HEL,'eq_2')
            CHCKERR('')
            ierr = reconstruct_PB3D_X_1(grid_X_HEL,X_HEL,'X_1')
            CHCKERR('')
            
            ! user output
            call lvl_ud(-1)
            call writo('HELENA output reconstructed for later interpolation')
        end if
        
        ! reconstruct normal part of full eq grid for XYZ reconstruction
        if (POST_output_full .and. plot_grid_style.eq.0 .or. &
            &plot_grid_style.eq.3) then
            ! user output
            call writo('Reconstructing full eq grid for XYZ reconstruction')
            call lvl_ud(1)
            
            ierr = reconstruct_PB3D_grid(grid_eq_XYZ,'eq',&
                &rich_lvl=rich_lvl_name,tot_rich=.true.,lim_pos=lim_loc)
            CHCKERR('')
            
            ! user output
            call lvl_ud(-1)
            call writo('Full eq grid reconstructed for XYZ reconstruction')
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('3-D plots prepared')
        
        ! clean up
        call grid_eq%dealloc()
        call grid_X%dealloc()
        call grid_sol%dealloc()
        call eq_1%dealloc()
    end function init_POST
    
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
    ! The general workflow is as follows:
    !   - take a subset of the output grids for the current equilibrium job.
    !   - for POST_style = 1 (extended grid):
    !       * VMEC: recalculate variables
    !       * HEL:  interpolate variables
    !     for POST_style = 2 (B-aligned grid):
    !       * VMEC: read subset of variables
    !       * HEL:  inerpolate variables
    !   - create helper variables
    !   - create plots and outputs
    integer function run_driver_POST() &
        &result(ierr)
        
        use num_vars, only: no_output, no_plots, eq_style, plot_grid_style, &
            &use_pol_flux_F, pi, swap_angles, eq_jobs_lims, eq_job_nr, &
            &plot_B, plot_J, plot_sol_xi, plot_sol_Q, plot_kappa, &
            &POST_output_sol, POST_style
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_1
        use grid_vars, only: disc_type
        use eq_vars, only: max_flux_F
        use eq_ops, only: calc_eq, calc_derived_q, calc_T_HF, B_plot, J_plot, &
            &kappa_plot
        use eq_utilities, only: calc_F_derivs, calc_inv_met
        use sol_vars, only: alpha
        use grid_utilities, only: calc_XYZ_grid, setup_interp_data, &
            &apply_disc, copy_grid
        use X_ops, only: calc_X, setup_nm_X
        use sol_ops, only: plot_sol_vec, decompose_energy
        use HELENA_ops, only: interp_HEL_on_grid
        use VMEC, onLy: calc_trigon_factors
        use num_utilities, only: calc_aux_utilities
#if ldebug
        use num_vars, only: ltest
#endif
        
        character(*), parameter :: rout_name = 'run_driver_POST'
        
        ! local variables
        type(grid_type) :: grids(3)                                             ! local output eq, X and sol grid, with limited parallel range and divided normal range
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium
        type(eq_2_type) :: eq_2                                                 ! metric equilibrium
        type(X_1_type) :: X                                                     ! vectorial perturbation variables
        integer :: id, jd, kd                                                   ! counter
        integer :: n_EV_out                                                     ! how many EV's are treated
        integer :: i_EV_out                                                     ! counter
        logical :: no_plots_loc                                                 ! local copy of no_plots
        logical :: no_output_loc                                                ! local copy of no_output
        logical :: last_eq_job                                                  ! last parallel job
        real(dp), allocatable :: XYZ(:,:,:,:)                                   ! X, Y and Z on output grid
        real(dp), allocatable :: R(:,:,:)                                       ! R on output grid
        complex(dp), allocatable :: sol_val_comp(:,:,:)                         ! solution eigenvalue for requested solutions
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        logical :: ltest_loc                                                    ! local copy of ltest
#endif
        
        ! initialize ierr
        ierr = 0
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('Parallel range for this job:')
        call lvl_ud(1)
        call writo(trim(i2str(eq_jobs_lims(1,eq_job_nr)))//'..'//&
            &trim(i2str(eq_jobs_lims(2,eq_job_nr)))//' of '//'1..'//&
            &trim(i2str(n_par_tot)))
        call lvl_ud(-1)
        
        ! some variables
        last_eq_job = eq_job_nr.eq.size(eq_jobs_lims,2)
        n_EV_out = sum(max_id-min_id+1)
        
        ! set up restricted output grids for this equilibrium job
        ierr = setup_out_grids(grids)
        CHCKERR('')
        
        ! set up output eq_1, eq_2 and X_1, depending on equilibrium style
        ierr = calc_eq(grids(1),eq_1)
        CHCKERR('')
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Calculate variables on output grid')
                call lvl_ud(1)
                
                ! back up no_plots and no_output and set to .true.
                no_plots_loc = no_plots; no_plots = .true.
                no_output_loc = no_output; no_output = .true.
#if ldebug
                ltest_loc = ltest; ltest = .false.
#endif
                
                ! Calculate the metric equilibrium quantitities
                ierr = calc_eq(grids(1),eq_1,eq_2)
                CHCKERR('')
                
                ! Transform E into F derivatives
                ierr = calc_F_derivs(eq_2)
                CHCKERR('')
                
                ! Calculate derived metric quantities
                call calc_derived_q(grids(1),eq_1,eq_2)
                
                ! calculate X variables, vector phase
                ierr = calc_X(grids(1),grids(2),eq_1,eq_2,X)
                CHCKERR('')
                
                ! reset no_plots and no_output
                no_plots = no_plots_loc
                no_output = no_output_loc
#if ldebug
                ltest = ltest_loc
#endif
                
                ! user output
                call lvl_ud(-1)
                call writo('Variables calculated on output grid')
            case (2)                                                            ! HELENA
                ! user output
                call writo('Interpolate variables on output grid')
                call lvl_ud(1)
                
                call eq_2%init(grids(1),setup_E=.true.,setup_F=.true.)
                call X%init(grids(2))
                ierr = interp_HEL_on_grid(grid_eq_HEL,grids(1),&
                    &eq_2=eq_2_HEL,eq_2_out=eq_2,eq_1=eq_1,&
                    &grid_name='output equilibrium grid')
                CHCKERR('')
                ierr = interp_HEL_on_grid(grid_X_HEL,grids(2),X_1=X_HEL,&
                    &X_1_out=X,grid_name='output perturbation grid')
                CHCKERR('')
                
                ! also set up transformation matrices for plots
                ierr = calc_T_HF(grids(1),eq_1,eq_2,[0,0,0])
                CHCKERR('')
                ierr = calc_inv_met(eq_2%T_FE,eq_2%T_EF,[0,0,0])
                CHCKERR('')
                
                ! user output
                call lvl_ud(-1)
                call writo('Variables interpolated on output grid')
        end select
        
        ! user output
        call writo('Initialize 3-D plots')
        call lvl_ud(1)
        
        ! initialize  integrated energies and  open decompose log file  if first
        ! equilibrium job
        if (POST_output_sol .and. eq_job_nr.eq.1) then
            allocate(E_pot_int(6,n_EV_out))
            allocate(E_kin_int(6,n_EV_out))
            E_pot_int = 0._dp
            E_kin_int = 0._dp
            
            ierr = open_decomp_log()
            CHCKERR('')
        end if
        
        ! if VMEC, calculate trigonometric factors of output grid
        if (eq_style.eq.1) then
            ierr = calc_trigon_factors(grids(2)%theta_E,&
                &grids(2)%zeta_E,grids(2)%trigon_factors)
            CHCKERR('')
        end if
        
        ! set up XYZ
        ! Note:  plot_grid_style is  only applied  for the  quantities that  are
        ! tabulated  on the  perturbation grid,  i.e. the  perturbation and  the
        ! magnetic perturbation, as well as the energy terms that are calculated
        ! from them.
        allocate(XYZ(grids(2)%n(1),grids(2)%n(2),grids(2)%loc_n_r,3))
        select case (plot_grid_style)
            case (0)                                                            ! 3-D geometry
                ! user output
                call writo('Plots are done in 3-D geometry')
                call lvl_ud(1)
                
                ierr = calc_XYZ_grid(grid_eq_XYZ,grids(2),&
                    &XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3))
                CHCKERR('')
                !!! To plot the cross-section
                !!call print_ex_2D(['cross_section'],'cross_section_'//&
                    !!&trim(i2str(eq_job_nr))//'_'//trim(i2str(rank)),&
                    !!&XYZ(:,1,:,3),x=XYZ(:,1,:,1),draw=.false.)
                !!call draw_ex(['cross_section'],'cross_section_'//&
                    !!&trim(i2str(eq_job_nr))//'_'//trim(i2str(rank)),&
                    !!&size(XYZ,3),1,.false.)
                
                call lvl_ud(-1)
            case (3)                                                            ! 3-D geometry with straightened toroidal coordinate
                ! user output
                call writo('Plots are done in an unwrapped torus')
                
                allocate(R(grids(2)%n(1),grids(2)%n(2),grids(2)%loc_n_r))
                ierr = calc_XYZ_grid(grid_eq_XYZ,grids(2),&
                    &XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),R=R)
                CHCKERR('')
                XYZ(:,:,:,1) = R                                                ! set X to R coordinate
                XYZ(:,:,:,2) = grids(2)%zeta_F                                  ! set Y to toroidal coordinate
                !!! To plot the cross-section
                !!call print_ex_2D(['cross_section'],'cross_section_'//&
                    !!&trim(i2str(eq_job_nr))//'_'//trim(i2str(rank)),&
                    !!&XYZ(:,1,:,3),x=XYZ(:,1,:,1),draw=.false.)
                !!call draw_ex(['cross_section'],'cross_section_'//&
                    !!&trim(i2str(eq_job_nr))//'_'//trim(i2str(rank)),&
                    !!&size(XYZ,3),1,.false.)
                
                call lvl_ud(-1)
            case (1,2)                                                          ! slab geometry (with or without wrapping to fundamental interval)
                ! user output
                call writo('Plots are done in slab geometry')
                call lvl_ud(1)
                
                ! set angular coordinates
                select case (POST_style)
                    case (1)                                                    ! extended grid: both theta and zeta need to be chosen, not alpha
                        if (use_pol_flux_F) then
                            XYZ(:,:,:,1) = grids(2)%theta_F/pi
                            XYZ(:,:,:,2) = grids(2)%zeta_F/pi
                        else
                            XYZ(:,:,:,1) = grids(2)%zeta_F/pi
                            XYZ(:,:,:,2) = grids(2)%theta_F/pi
                        end if
                    case (2)                                                    ! field-aligned grid: alpha in combination with parallel angle
                        if (swap_angles) then
                            call writo('For plotting, the angular coordinates &
                                &are swapped: theta <-> zeta')
                            if (use_pol_flux_F) then
                                XYZ(:,:,:,1) = grids(2)%zeta_F/pi
                            else
                                XYZ(:,:,:,1) = grids(2)%theta_F/pi
                            end if
                        else
                            if (use_pol_flux_F) then
                                XYZ(:,:,:,1) = grids(2)%theta_F/pi
                            else
                                XYZ(:,:,:,1) = grids(2)%zeta_F/pi
                            end if
                        end if
                        XYZ(:,:,:,2) = alpha
                end select
                
                ! set normal coordinate
                do kd = 1,grids(2)%loc_n_r
                    XYZ(:,:,kd,3) = grids(2)%loc_r_F(kd)/max_flux_F*2*pi
                end do
                
                ! limit to fundamental interval -1..1 if plot grid style 2
                if (plot_grid_style.eq.2) XYZ(:,:,:,1:2) = &
                    &modulo(XYZ(:,:,:,1:2)+1._dp,2._dp)-1._dp
                
                call lvl_ud(-1)
            case default
                err_msg = 'No style "'//trim(i2str(plot_grid_style))//&
                    &'" for plot_grid_style'
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! user output
        call lvl_ud(-1)
        call writo('3-D plots initialized')
        
        ! user output
        call writo('3-D PB3D output plots')
        call lvl_ud(1)
        
        ! plot magnetic field if requested
        if (plot_B) then
            call writo('Plot the magnetic field')
            call lvl_ud(1)
            ierr = B_plot(grids(1),eq_1,eq_2,plot_fluxes=POST_style.eq.1)
            CHCKERR('')
            call lvl_ud(-1)
        end if
        
        ! plot current if requested
        if (plot_J) then
            call writo('Plot the current')
            call lvl_ud(1)
            ierr = J_plot(grids(1),eq_1,eq_2,plot_fluxes=POST_style.eq.1)
            CHCKERR('')
            call lvl_ud(-1)
        end if
        
        ! plot magnetic field if requested
        if (plot_kappa) then
            call writo('Plot the curvature')
            call lvl_ud(1)
            ierr = kappa_plot(grids(1),eq_1,eq_2)
            CHCKERR('')
            call lvl_ud(-1)
        end if
        
        if (POST_output_sol) then
            ! user output
            call writo('Solution plots for different ranges')
            call lvl_ud(1)
            
            ! loop over three ranges
            i_EV_out = 1
            if (last_eq_job) allocate(sol_val_comp(2,2,n_EV_out))
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
                    
                    ! plot solution vector
                    ierr = plot_sol_vec(grids(1),grids(2),grids(3),&
                        &eq_1,eq_2,X,sol,XYZ,id,[plot_sol_xi,plot_sol_Q])
                    CHCKERR('')
                    
                    ! user output
                    call writo('Decompose the energy into its terms')
                    call lvl_ud(1)
                    ierr = decompose_energy(grids(1),grids(2),grids(3),&
                        &eq_1,eq_2,X,sol,id,B_aligned=POST_style.eq.2,XYZ=XYZ,&
                        &E_pot_int=E_pot_int(:,i_EV_out),&
                        &E_kin_int=E_kin_int(:,i_EV_out))
                    CHCKERR('')
                    call lvl_ud(-1)
                    
                    ! write to  log file and save  difference between Eigenvalue
                    ! and energy fraction if last parallel job
                    if (last_eq_job) then
                        ierr = write_decomp_log(id,E_pot_int(:,i_EV_out),&
                            &E_kin_int(:,i_EV_out))
                        CHCKERR('')
                        
                        sol_val_comp(:,1,i_EV_out) = id*1._dp
                        sol_val_comp(:,2,i_EV_out) = [sol%val(id),&
                            &sum(E_pot_int(:,i_EV_out))/&
                            &sum(E_kin_int(:,i_EV_out))]
                    end if
                    
                    ! increment counter
                    i_EV_out = i_EV_out + 1
                    
                    ! user output
                    call lvl_ud(-1)
                    call writo('Mode '//trim(i2str(id))//'/'//&
                        &trim(i2str(size(sol%val)))//' finished')
                end do
                
                call lvl_ud(-1)
                if (min_id(jd).le.max_id(jd)) call &
                    &writo('RANGE '//trim(i2str(jd))//' plotted')
            end do
            
            ! plot  difference between  Eigenvalue and  energy fraction  if last
            ! parallel job
            if (last_eq_job) then
                call plot_sol_val_comp(sol_val_comp)
            end if
            
            ! user output
            call lvl_ud(-1)
            call writo('Solution for different ranges plotted')
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('3-D output plots done')
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call grids(1)%dealloc()
        call grids(2)%dealloc()
        call grids(3)%dealloc()
        call eq_1%dealloc()
        call eq_2%dealloc()
        call X%dealloc()
        call lvl_ud(-1)
    end function run_driver_POST
    
    ! Cleans up main driver for postprocessing.
    subroutine stop_POST()
        use num_vars, only: plot_grid_style, eq_style, POST_output_full, &
            &POST_output_sol
        use input_utilities, only: dealloc_in
        
        call dealloc_in()
        if (POST_output_sol) call sol%dealloc()
        if (POST_output_full .and. eq_style.eq.2) then
            call grid_eq_HEL%dealloc()
            call grid_X_HEL%dealloc()
            call eq_2_HEL%dealloc()
            call X_HEL%dealloc()
        end if
        if (POST_output_full .and. plot_grid_style.eq.0) then
            call grid_eq_XYZ%dealloc()
        end if
    end subroutine stop_POST
    
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
    
    ! write to decomposition log file
    integer function write_decomp_log(X_id,E_pot_int,E_kin_int) result(ierr)
        use num_vars, only: decomp_i, rank
        
        character(*), parameter :: rout_name = 'write_decomp_log'
        
        ! input / output
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        complex(dp), intent(in) :: E_pot_int(6)                                 ! E_pot integrated for requested solutions
        complex(dp), intent(in) :: E_kin_int(2)                                 ! E_kin integrated for requested solutions
        
        ! local variables
        character(len=max_str_ln) :: format_val                                 ! format
        character(len=2*max_str_ln) :: temp_output_str                          ! temporary output string
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write to log file')
        call lvl_ud(1)
        
        ! only master
        if (rank.eq.0) then
            ! set up format string:
            !   row 1: EV, E_pot/E_kin
            !   row 2: E_pot, E_kin
            !   row 3: E_kin(1), E_kin(2)
            !   row 4: E_pot(1), E_pot(2)
            !   row 5: E_pot(3), E_pot(4)
            !   row 6: E_pot(5), E_pot(6)
            format_val = '("  ",ES23.16," ",ES23.16," ",&
                &ES23.16," ",ES23.16," ",ES23.16," ",ES23.16," ",ES23.16)'
            
            ! write header
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &'# Eigenvalue '//trim(i2str(X_id))
            CHCKERR('Failed to write')
            
            ! write values
            write(temp_output_str,format_val) &
                &rp(sol%val(X_id)),&
                &rp(sum(E_pot_int)/sum(E_kin_int)),&
                &rp(sum(E_pot_int)),&
                &rp(sum(E_kin_int))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_val) &
                &ip(sol%val(X_id)),&
                &ip(sum(E_pot_int)/sum(E_kin_int)), &
                &ip(sum(E_pot_int)),&
                &ip(sum(E_kin_int))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_val) &
                &rp(E_kin_int(1)),&
                &rp(E_kin_int(2))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_val) &
                &ip(E_kin_int(1)),&
                &ip(E_kin_int(2))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_val) &
                &rp(E_pot_int(1)),&
                &rp(E_pot_int(2)),&
                &rp(E_pot_int(3)),&
                &rp(E_pot_int(4)),&
                &rp(E_pot_int(5)),&
                &rp(E_pot_int(6))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &trim(temp_output_str)
            CHCKERR('Failed to write')
            write(temp_output_str,format_val) &
                &ip(E_pot_int(1)),&
                &ip(E_pot_int(2)),&
                &ip(E_pot_int(3)),&
                &ip(E_pot_int(4)),&
                &ip(E_pot_int(5)),&
                &ip(E_pot_int(6))
            write(UNIT=decomp_i,FMT='(A)',IOSTAT=ierr) &
                &trim(temp_output_str)
            CHCKERR('Failed to write')
            
            close(decomp_i)
        end if
        
        call lvl_ud(-1)
    end function write_decomp_log
    
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
    
    ! Sets up the output grid for a particular parallel job
    integer function setup_out_grids(grids_out) result(ierr)
        use num_vars, only: min_theta_plot, max_theta_plot, POST_style, &
            &eq_job_nr, eq_jobs_lims
        use grid_utilities, only: extend_grid_E, copy_grid
        use PB3D_ops, only: reconstruct_PB3D_grid
        
        character(*), parameter :: rout_name = 'setup_out_grids'
        
        ! input / output
        type(grid_type), intent(inout) :: grids_out(3)                          ! eq, X and sol output grids (full parallel and normal range)
        
        ! local variables
        integer :: lim_loc(3,2,3)                                               ! grid ranges for local equilibrium job (last index: eq, X, sol)
        real(dp) :: lim_theta(2)                                                ! theta limits
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        integer :: n_theta                                                      ! nr. of points in theta for extended grid
        ! initialize ierr
        ierr = 0
        
        call writo('Set up output grid')
        call lvl_ud(1)
        
        select case (POST_style)
            case (1)                                                            ! extended grid
                ! limits to exclude angular part
                ! eq
                lim_loc(1,:,1) = [0,0]
                lim_loc(2,:,1) = [0,0]
                lim_loc(3,:,1) = [1,-1]
                ! X
                lim_loc(1,:,2) = [0,0]
                lim_loc(2,:,2) = [0,0]
                lim_loc(3,:,2) = [1,-1]
                ! sol
                lim_loc(1,:,3) = [0,0]
                lim_loc(2,:,3) = [0,0]
                lim_loc(3,:,3) = [1,-1]
                
                ! reconstruct equilibrium, perturbation and solution grids
                ierr = reconstruct_PB3D_grid(grid_eq,'eq',&
                    &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                    &lim_pos=lim_loc(:,:,1),grid_limits=lims_norm(:,1))
                CHCKERR('')
                ierr = reconstruct_PB3D_grid(grid_X,'X',&
                    &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                    &lim_pos=lim_loc(:,:,2),grid_limits=lims_norm(:,2))
                CHCKERR('')
                ierr = reconstruct_PB3D_grid(grids_out(3),'sol',&
                    &lim_pos=lim_loc(:,:,3),grid_limits=lims_norm(:,3))
                CHCKERR('')
                
                ! theta limits
                if (n_par_tot.gt.1) then
                    lim_theta = min_theta_plot + &
                        &(eq_jobs_lims(:,eq_job_nr)-1._dp)*&
                        &(max_theta_plot-min_theta_plot)/(n_par_tot-1._dp)
                else
                    lim_theta = min_theta_plot
                end if
                
                ! extend eq and X grid
                n_theta = eq_jobs_lims(2,eq_job_nr)-eq_jobs_lims(1,eq_job_nr)+1
                ierr = extend_grid_E(grid_eq,grids_out(1),n_theta_plot=n_theta,&
                    &lim_theta_plot=lim_theta,grid_eq=grid_eq)                  ! extend eq grid and convert to F
                CHCKERR('')
                ierr = extend_grid_E(grid_X,grids_out(2),n_theta_plot=n_theta,&
                    &lim_theta_plot=lim_theta,grid_eq=grid_eq)                  ! extend eq grid and convert to F
                CHCKERR('')
                
                ! clean up
                call grid_eq%dealloc()
                call grid_X%dealloc()
            case (2)                                                            ! field-aligned grid
                ! limits to take subset of angular part
                ! eq
                lim_loc(1,:,1) = eq_jobs_lims(:,eq_job_nr)
                lim_loc(2,:,1) = [1,-1]
                lim_loc(3,:,1) = [1,-1]
                ! X
                lim_loc(1,:,2) = eq_jobs_lims(:,eq_job_nr)
                lim_loc(2,:,2) = [1,-1]
                lim_loc(3,:,2) = [1,-1]
                ! sol
                lim_loc(1,:,3) = [0,0]
                lim_loc(2,:,3) = [0,0]
                lim_loc(3,:,3) = [1,-1]
                
                ! reconstruct equilibrium, perturbation and solution grids
                ierr = reconstruct_PB3D_grid(grids_out(1),trim(grid_name(1)),&
                    &rich_lvl=rich_lvl,tot_rich=.true.,&
                    &lim_pos=lim_loc(:,:,1),grid_limits=lims_norm(:,1))
                CHCKERR('')
                ierr = reconstruct_PB3D_grid(grids_out(2),trim(grid_name(2)),&
                    &rich_lvl=rich_lvl,tot_rich=.true.,&
                    &lim_pos=lim_loc(:,:,2),grid_limits=lims_norm(:,2))
                CHCKERR('')
                ierr = reconstruct_PB3D_grid(grids_out(3),trim(grid_name(3)),&
                    &lim_pos=lim_loc(:,:,3),grid_limits=lims_norm(:,3))
                CHCKERR('')
        end select
        
        call lvl_ud(-1)
    end function setup_out_grids
end module driver_POST
