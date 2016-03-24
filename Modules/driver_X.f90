!------------------------------------------------------------------------------!
!   Driver of the perturbation part of PB3D.                                   !
!------------------------------------------------------------------------------!
module driver_X
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type
    
    implicit none
    private
    public run_driver_X
    
    ! global variables
    character(len=8) :: flux_name(2)                                            ! name of flux variable
    character(len=1) :: mode_name(2)                                            ! name of modes
    integer :: rich_lvl_name                                                    ! either the Richardson level or zero, to append to names
    
contains
    ! Main driver of PB3D perturbation part.
    integer function run_driver_X() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, rank, plot_resonance, &
            &X_style, rich_restart_lvl
        use rich_vars, only: n_par_X, rich_lvl
        use MPI_utilities, only: wait_MPI
        use X_vars, only: dealloc_X, &
            &min_sec_X, max_sec_X, prim_X, min_r_sol, max_r_sol, n_mod_X
        use grid_vars, only: dealloc_grid, &
            &n_r_sol, min_par_X, max_par_X
        use eq_vars, only: dealloc_eq
        use grid_ops, only: calc_norm_range, setup_grid_X, &
            &print_output_grid
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_2
        use utilities, only: test_max_memory
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use X_ops, only: calc_X, check_X_modes, resonance_plot, &
            &setup_nm_X, print_output_X
        use input_utilities, only: dealloc_in
        !!use utilities, only: calc_aux_utilities
#if ldebug
        use X_vars, only: create_X
        use eq_vars, only: create_eq
#endif
        
        character(*), parameter :: rout_name = 'run_driver_X'
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B                                   ! field-aligned equilibrium grid
        type(grid_type), target :: grid_X                                       ! perturbation grid
        type(grid_type), pointer :: grid_X_B                                    ! perturbation grid
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium variables
        integer :: X_limits(2)                                                  ! min. and max. index of X grid for this process
        real(dp), allocatable :: r_F_X(:)                                       ! normal points in perturbation grid
        
        ! initialize ierr
        ierr = 0
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
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
        
        ! user output
        call writo('Perturbations are analyzed')
        call lvl_ud(1)
        
        call writo('for '//trim(i2str(n_r_sol))//' values on normal range '//&
            &trim(r2strt(min_r_sol))//'..'//trim(r2strt(max_r_sol)))
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
        
        ! set flux and mode names
        if (use_pol_flux_F) then
            flux_name = ['poloidal','toroidal']
            mode_name = ['m','n']
        else
            flux_name = ['toroidal','poloidal']
            mode_name = ['n','m']
        end if
        
        ! user output
        call writo('with '//flux_name(2)//' mode number '//mode_name(2)//&
            &' = '//trim(i2str(prim_X)))
        select case(X_style)
            case (1)                                                            ! prescribed
                call writo('and '//flux_name(1)//' mode number '//&
                    &mode_name(1)//' = '//trim(i2str(min_sec_X))//'..'//&
                    &trim(i2str(max_sec_X)))
            case (2)                                                            ! fast
                call writo('and the '//trim(i2str(n_mod_X))//&
                    &' '//flux_name(1)//' modes '//mode_name(1)//&
                    &' that resonate most')
        end select
        
        call lvl_ud(-1)
        
        ! reconstruct equilibrium grid and equilibrium output variables
        ! user output
        call writo('Setting up perturbation grid')
        call lvl_ud(1)
        
        call writo('Reconstruct equilibrium grid and variables')
        call lvl_ud(1)
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name)
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
        CHCKERR('')
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! do nothing: grid is already field-aligned
            case (2)                                                            ! HELENA
                ! also reconstruct field-aligned equilibrium grid
                allocate(grid_eq_B)
                ierr = reconstruct_PB3D_grid(grid_eq_B,'eq_B',rich_lvl=rich_lvl)
                CHCKERR('')
        end select
        call lvl_ud(-1)
        
        ! Divide perturbation grid under group processes, calculating the limits
        ! and the normal coordinate.
        call writo('Calculate normal coordinate')
        call lvl_ud(1)
        allocate(r_F_X(n_r_sol))
        ierr = calc_norm_range(X_limits=X_limits,r_F_X=r_F_X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! create and setup perturbation grid with division limits
        call writo('Calculate grid')
        call lvl_ud(1)
        ierr = setup_grid_X(grid_eq,grid_X,r_F_X,X_limits)
        CHCKERR('')
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! do nothing: grid is already field-aligned
            case (2)                                                            ! HELENA
                ! also set up field-aligned perturbation grid
                allocate(grid_X_B)
                ierr = setup_grid_X(grid_eq_B,grid_X_B,r_F_X,X_limits)
                CHCKERR('')
        end select
        call lvl_ud(-1)
        
        ! write to output if master
        if (rank.eq.0) then
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! print output
                    ierr = print_output_grid(grid_X,'perturbation',&
                        &'X',rich_lvl=rich_lvl)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ! print output
                    if (rich_lvl.eq.1) then                                     ! only write once
                        ierr = print_output_grid(grid_X,'perturbation','X')
                        CHCKERR('')
                    end if
                    ! also print field-aligned output
                    ierr = print_output_grid(grid_X_B,'field-aligned &
                        &perturbation','X_B',rich_lvl=rich_lvl)
                    CHCKERR('')
            end select
        end if
        
        call lvl_ud(-1)
        call writo('Perturbation grid set up')
        
        ! setup nm and check the X modes if lowest (or restart) Richardson level
        if (rich_lvl.eq.rich_restart_lvl) then
            ! setup limits of modes
            ierr = setup_nm_X(grid_eq,grid_X,eq_1,plot_nm=.true.)
            CHCKERR('')
            
            ! tests
            ierr = check_X_modes(grid_eq,eq_1)
            CHCKERR('')
        end if
        
        ! plot resonances if requested
        if (plot_resonance .and. rank.eq.0 .and. rich_lvl.eq.1) then
            ierr = resonance_plot(eq_1,grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! Run the  first part of the  driver. For VMEC  this has to be  done for
        ! every Richardson level, but for HELENA only for the first one.
        if (eq_style.eq.1 .or. (eq_style.eq.2 .and. rich_lvl.eq.1)) then
            ierr = run_driver_X_1(grid_eq,grid_X,eq_1)
            CHCKERR('')
        end if
        
        ! Run  the second  part of  the driver. This  has to  be done  for every
        ! Richardson level.
        ierr = run_driver_X_2(eq_1)
        CHCKERR('')
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call dealloc_in()
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_eq(eq_1)
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! do nothing: everything already field-aligned
            case (2)                                                            ! HELENA
                ! clean up
                call dealloc_grid(grid_eq_B)
                call dealloc_grid(grid_X_B)
                deallocate(grid_eq_B)
                deallocate(grid_X_B)
                nullify(grid_eq_B)
                nullify(grid_X_B)
        end select
        call lvl_ud(-1)
    end function run_driver_X
    
    ! Part 1 of driver_X: up to the calculation of magnetic integrals.
    integer function run_driver_X_1(grid_eq,grid_X,eq_1) result(ierr)
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use MPI_utilities, only: wait_MPI
        use num_vars, only: X_job_nr, X_jobs_lims, rank, n_procs
        use PB3D_ops, only: reconstruct_PB3D_eq_2, reconstruct_PB3D_X_1
        use X_ops, only: calc_X, print_output_X
        use eq_vars, only: dealloc_eq
        use X_vars, only: dealloc_X
        use vac, only: calc_vac
        
        character(*), parameter :: rout_name = 'run_driver_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        
        ! local variables
        logical :: dim_reused(2)                                                ! wether dimension is reused
        integer :: id                                                           ! counter
        integer :: arr_size                                                     ! size of array for jobs division
        type(eq_2_type) :: eq_2                                                 ! metric equilibrium variables
        type(X_1_type) :: X_1(2)                                                ! vectorial X variables
        type(X_2_type) :: X_2                                                   ! tensorial X variables
        
        ! initialize ierr
        ierr = 0
        
        !!!!!!!!!!!!!!!!!!!!!
        ! 1. Vectorial jobs !
        !!!!!!!!!!!!!!!!!!!!!
        
        call writo('Setting up vectorial perturbation jobs')
        call lvl_ud(1)
        
        ! divide perturbation jobs
        arr_size = grid_X%loc_n_r*product(grid_X%n(1:2))
        ierr = divide_X_jobs(arr_size,1)
        CHCKERR('')
        
        ! get metric equilibrium quantities
        ierr = reconstruct_PB3D_eq_2(grid_eq,eq_2,'eq_2',rich_lvl=rich_lvl_name)
        CHCKERR('')
        
        call lvl_ud(-1)
        call writo('Vectorial perturbation jobs set up')
        
        ! main loop over vectorial jobs
        X_job_nr = 0
        call writo('Starting vectorial perturbation jobs')
        call lvl_ud(1)
        X_jobs_1: do
            ! get next job
            ierr = get_next_job(X_job_nr)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call print_info_X_1()
            
            ! calc variables
            ierr = calc_X(grid_eq,grid_X,eq_1,eq_2,X_1(1),&
                &lim_sec_X=X_jobs_lims(:,X_job_nr))
            CHCKERR('')
            
            ! write vectorial perturbation variables to output
            ierr = print_output_X(grid_X,X_1(1),'X_1',&
                &rich_lvl=rich_lvl_name,lim_sec_X=X_jobs_lims(:,X_job_nr))
            CHCKERR('')
            
            !!! plot information for comparison between VMEC and HELENA
            !!call plot_info_for_VMEC_HEL_comparison(grid_X,X_1(1))
            
            ! clean up
            call dealloc_X(X_1(1))
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed &
                &by process '//trim(i2str(rank)))
        end do X_jobs_1
        call lvl_ud(-1)
        call writo('Vectorial perturbation jobs finished')
        
        ! user output
        if (n_procs.gt.1) call writo('Waiting for the other processes')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        !!!!!!!!!!!!!!!!!!!!!
        ! 2. Tensorial jobs !
        !!!!!!!!!!!!!!!!!!!!!
        
        call writo('Setting up Tensorial perturbation jobs')
        call lvl_ud(1)
        
        ! divide perturbation jobs, tensor phase
        arr_size = grid_X%loc_n_r*product(grid_X%n(1:2))
        ierr = divide_X_jobs(arr_size,2)
        CHCKERR('')
        
        call lvl_ud(-1)
        call writo('Tensorial perturbation jobs set up')
        
        ! main loop over tensorial jobs
        X_job_nr = 0
        call writo('Starting tensorial perturbation jobs')
        call lvl_ud(1)
        X_jobs_2: do
            ! get next job
            ierr = get_next_job(X_job_nr,dim_reused)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call print_info_X_2()
            
            ! retrieve vectorial perturbation if dimension not reused
            do id = 1,2
                if (dim_reused(id)) then
                    call writo('Vectorial perturbation variables for &
                        &dimension '//trim(i2str(id))//' reused')
                else
                    ! free the variable
                    if (allocated(X_1(id)%n)) call dealloc_X(X_1(id))
                    
                    ! user output
                    call writo('Requesting vectorial perturbation &
                        &variables for dimension '//trim(i2str(id)))
                    call lvl_ud(1)
                    
                    ! reconstruct PB3D X_1 quantities for this dimension
                    ierr = reconstruct_PB3D_X_1(grid_X,X_1(id),'X_1',&
                        &rich_lvl=rich_lvl_name,lim_sec_X=&
                        &X_jobs_lims((id-1)*2+1:(id-1)*2+2,X_job_nr))
                    CHCKERR('')
                    
                    ! user output
                    call lvl_ud(-1)
                    call writo('Vectorial perturbation variables for &
                        &dimension '//trim(i2str(id))//' loaded')
                end if
            end do
            
            ! calculate X variables, tensor phase
            ierr = calc_X(grid_eq,grid_X,eq_1,eq_2,X_1(1),X_1(2),X_2,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! calculate vacuum response
            ierr = calc_vac(X_2)
            CHCKERR('')
            
            ! write tensorial perturbation variables to output file
            ierr = print_output_X(grid_X,X_2,'X_2',rich_lvl=rich_lvl_name,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X_2)
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//&
                &' completed by process '//trim(i2str(rank)))
        end do X_jobs_2
        call lvl_ud(-1)
        call writo('Tensorial perturbation jobs finished')
        
        ! clean up
        call writo('Clean up tensorial perturbation variables')
        call lvl_ud(1)
        call dealloc_eq(eq_2)
        do id = 1,2
            call dealloc_X(X_1(id))
        end do
        call lvl_ud(-1)
        
        ! user output
        if (n_procs.gt.1) call writo('Waiting for the other processes')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
#if ldebug
    contains
        subroutine plot_info_for_VMEC_HEL_comparison(grid_X,X)
            use input_utilities, only: pause_prog, get_int, get_log
            use num_vars, only: use_pol_flux_F
            
            ! input / output
            type(grid_type), intent(in) :: grid_X
            type(X_1_type), intent(in) :: X
            
            ! local variables
            real(dp), allocatable :: x_plot(:,:,:), y_plot(:,:,:), z_plot(:,:,:)
            integer :: id, ld
            logical :: not_ready = .true.
            
            write(*,*) '!!!! PLOTTING INFORMATION FOR COMPARISON BETWEEN &
                &VMEC AND HELENA !!!!'
            do while (not_ready)
                call writo('Plots for which mode number?')
                ld =  get_int(lim_lo=1,lim_hi=X%n_mod)
                
                ! x, y, z
                allocate(x_plot(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(y_plot(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(z_plot(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                if (use_pol_flux_F) then
                    x_plot = grid_X%theta_F
                else
                    x_plot = grid_X%zeta_F
                end if
                y_plot = 0._dp
                do id = 1,grid_X%loc_n_r
                    z_plot(:,:,id) = grid_X%loc_r_F(id)/maxval(grid_X%r_F)
                end do
                call plot_HDF5('x_plot','x_plot',x_plot)
                call plot_HDF5('y_plot','y_plot',y_plot)
                call plot_HDF5('z_plot','z_plot',z_plot)
                
                ! U_0 and U_1
                call plot_HDF5('RE_U_0','RE_U_0',realpart(X%U_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_U_0','IM_U_0',imagpart(X%U_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('RE_U_1','RE_U_1',realpart(X%U_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_U_1','IM_U_1',imagpart(X%U_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                
                ! DU_0 and DU_1
                call plot_HDF5('RE_DU_0','RE_DU_0',realpart(X%DU_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_DU_0','IM_DU_0',imagpart(X%DU_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('RE_DU_1','RE_DU_1',realpart(X%DU_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_DU_1','IM_DU_1',imagpart(X%DU_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                
                ! clean up
                deallocate(x_plot,y_plot,z_plot)
                
                call writo('Want to redo the plotting?')
                not_ready = get_log(.true.)
            end do
            
            write(*,*) '!!!! DONE, PAUSED !!!!'
            call pause_prog()
        end subroutine plot_info_for_VMEC_HEL_comparison
#endif
    end function run_driver_X_1
    
    ! Part 2 of driver_X: calculation of magnetic integrals.
    integer function run_driver_X_2(eq_1) result(ierr)
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use MPI_utilities, only: wait_MPI
        use num_vars, only: X_job_nr, X_jobs_lims, rank, n_procs, eq_style
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_1, reconstruct_PB3D_X_2
        use rich_vars, only: rich_lvl
        use X_ops, only: calc_X, print_output_X, calc_magn_ints
        use grid_vars, only: dealloc_grid
        use eq_vars, only: create_eq, dealloc_eq
        use X_vars, only: dealloc_X
        use HELENA_ops, only: interp_HEL_on_grid
        use vac, only: calc_vac
        
        character(*), parameter :: rout_name = 'run_driver_X_2'
        
        ! input / output
        type(eq_1_type), intent(in) :: eq_1                                     ! vectorial perturbation variables
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B                                   ! field-aligned equilibrium grid
        type(grid_type), target :: grid_X                                       ! perturbation grid
        type(grid_type), pointer :: grid_X_B                                    ! perturbation grid
        type(eq_2_type), target :: eq_2                                         ! tensorial perturbation variables
        type(eq_2_type), pointer :: eq_2_B                                      ! field-aligned tensorial perturbation variables
        type(X_2_type) :: X_2                                                   ! tensorial X variables
        integer :: arr_size                                                     ! size of array for jobs division
        
        ! initialize ierr
        ierr = 0
        
        call writo('Setting up magnetic averaging jobs')
        call lvl_ud(1)
        
        call writo('Reconstruct equilibrium and perturbation grids and &
            &variables')
        call lvl_ud(1)
        
        ! Reconstruct equilibrium  and perturbation grid and  metric equilibrium
        ! variables using all Richardson levels.
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.)
        CHCKERR('')
        ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
            &tot_rich=.true.)
        CHCKERR('')
        ierr = reconstruct_PB3D_eq_2(grid_eq,eq_2,'eq_2',&
            &rich_lvl=rich_lvl_name,tot_rich=.true.)
        CHCKERR('')
        
        ! Setup  field-aligned  equilibrium  and perturbation  grid  and  metric
        ! equilibirum variables depending on equilibrium type.
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! point
                grid_eq_B => grid_eq
                grid_X_B => grid_X
                eq_2_B => eq_2
            case (2)                                                            ! HELENA
                ! allocate
                allocate(grid_eq_B)
                allocate(grid_X_B)
                allocate(eq_2_B)
                
                ! reconstruct field-aligned equilibrium and perturbation grids
                ierr = reconstruct_PB3D_grid(grid_eq_B,'eq_B',&
                    &rich_lvl=rich_lvl)
                CHCKERR('')
                ierr = reconstruct_PB3D_grid(grid_X_B,'X_B',&
                    &rich_lvl=rich_lvl)
                CHCKERR('')
                
                ! interpolate field-aligned metric equilibrium variables
                call create_eq(grid_eq_B,eq_2_B)
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,eq_2=eq_2,&
                    &eq_2_out=eq_2_B,eq_1=eq_1,grid_name='field-aligned &
                    &equilibrium grid')
                CHCKERR('')
        end select
        
        call lvl_ud(-1)
        
        ! divide perturbation jobs, tensor phase
        select case (eq_style)
            case (1)                                                            ! VMEC
                arr_size = grid_X%loc_n_r*product(grid_X_B%n(1:2))
            case (2)                                                            ! HELENA
                arr_size = grid_X%loc_n_r*&
                    &max(product(grid_X%n(1:2)),product(grid_X_B%n(1:2)))
        end select
        ierr = divide_X_jobs(arr_size,2)
        CHCKERR('')
        
        call lvl_ud(-1)
        call writo('Magnetic averaging jobs set up')
        
        ! main loop over tensorial jobs
        X_job_nr = 0
        call writo('Starting magnetic averaging jobs')
        call lvl_ud(1)
        X_jobs_3: do
            ! get next job
            ierr = get_next_job(X_job_nr)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call print_info_X_2()
            
            ! reconstruct the tensorial perturbation variables
            ierr = reconstruct_PB3D_X_2(grid_X,X_2,'X_2',&
                &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! Adapt to field-aligned grid if necessary
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! do nothing
                case (2)                                                        ! HELENA
                    ! interpolate tensorial perturbations to field-aligned grid
                    ierr = interp_HEL_on_grid(grid_X,grid_X_B,X_2=X_2,&
                        &grid_name='field-aligned perturbation grid')
                    CHCKERR('')
            end select
            
            ! integrate magnetic  integrals of tensorial  perturbation variables
            ! over field-aligned grid and write output in same variable X_2.
            ierr = calc_magn_ints(grid_eq_B,grid_X_B,eq_2_B,X_2,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! write field-averaged tensorial perturbation variables to output
            ierr = print_output_X(grid_X_B,X_2,'X_2_int',rich_lvl=rich_lvl,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]),&
                &is_field_averaged=.true.)
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X_2)
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//&
                &' completed by process '//trim(i2str(rank)))
        end do X_jobs_3
        call lvl_ud(-1)
        call writo('Magnetic averaging jobs done')
        
        ! clean up
        call writo('Clean up magnetic averaging variables')
        call lvl_ud(1)
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        call dealloc_eq(eq_2)
        if (eq_style.eq.2) then
            call dealloc_grid(grid_eq_B)
            call dealloc_grid(grid_X_B)
            call dealloc_eq(eq_2_B)
            deallocate(grid_eq_B)
            deallocate(grid_X_B)
            deallocate(eq_2_B)
        end if
        nullify(grid_eq_B)
        nullify(grid_X_B)
        nullify(eq_2_B)
        call lvl_ud(-1)
        
        ! user output
        if (n_procs.gt.1) call writo('Waiting for the other processes')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
    end function run_driver_X_2
    
    ! prints information for vectorial perturbation job
    subroutine print_info_X_1()
        use num_vars, only: rank, X_job_nr, X_jobs_lims
        
        call writo('Job '//trim(i2str(X_job_nr))//' is started by process '&
            &//trim(i2str(rank)))
        call lvl_ud(1)
        
        call writo('The series of '//flux_name(1)//' mode numbers ('//&
            &trim(i2str(X_jobs_lims(1,X_job_nr)))//'..'//&
            &trim(i2str(X_jobs_lims(2,X_job_nr)))//') is calculated')
    end subroutine print_info_X_1
    
    ! prints information for tensorial perturbation job
    subroutine print_info_X_2()
        use num_vars, only: rank, X_job_nr, X_jobs_lims
        
        ! user output
        call writo('Job '//trim(i2str(X_job_nr))//' is started by process '&
            &//trim(i2str(rank)))
        call lvl_ud(1)
        
        call writo('The block of '//flux_name(1)//' mode numbers ('//&
            &trim(i2str(X_jobs_lims(1,X_job_nr)))//'..'//&
            &trim(i2str(X_jobs_lims(2,X_job_nr)))//')x('//&
            &trim(i2str(X_jobs_lims(3,X_job_nr)))//'..'//&
            &trim(i2str(X_jobs_lims(4,X_job_nr)))//') is calculated')
    end subroutine print_info_X_2
end module driver_X
