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
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type, X_2_type
    
    implicit none
    private
    public run_driver_X
#if ldebug
    public debug_run_driver_X
#endif
    
    ! global variables
#if ldebug
    logical :: debug_run_driver_X = .false.                                     ! debug run_driver_X
#endif
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_driver_X() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, rank, plot_resonance, &
            &X_job_nr, X_jobs_lims, n_procs, X_style
        use rich, only: n_r_sol
        use MPI_utilities, only: wait_MPI
        use X_vars, only: dealloc_X, &
            &min_sec_X, max_sec_X, prim_X, min_r_sol, max_r_sol, n_mod_X
        use grid_vars, only: dealloc_grid, &
            &n_par_X, min_par_X, max_par_X
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met, create_met
        use grid_ops, only: calc_norm_range, setup_and_calc_grid_X, &
            &print_output_grid
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use utilities, only: test_max_memory
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use X_ops, only: calc_X, check_X_modes, resonance_plot, &
            &setup_nm_X, print_output_X, calc_magn_ints
        use vac, only: calc_vac
        use HELENA_vars, only: dealloc_HEL
        use HELENA_ops, only: interp_HEL_on_grid
        use VMEC, only: dealloc_VMEC
        use rich, only: rich_info_short
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
        type(eq_type), target :: eq                                             ! equilibrium variables
        type(eq_type), pointer :: eq_B                                          ! field-aligned equilibrium variables
        type(met_type), target :: met                                           ! metric variables
        type(met_type), pointer :: met_B                                        ! field-aligned metric variables
        type(X_1_type) :: X_1(2)                                                ! vectorial X variables
        type(X_2_type) :: X_2                                                   ! tensorial X variables
        integer :: X_limits(2)                                                  ! min. and max. index of X grid for this process
        real(dp), allocatable :: r_F_X(:)                                       ! normal points in perturbation grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: dim_reused(2)                                                ! wether dimension is reused
        integer :: id                                                           ! counter
        character(len=8) :: flux_name                                           ! name of flux variable
#if ldebug
        type(X_1_type) :: X_1_B                                                 ! field-aligned vectorial perturbation variables
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! some preliminary things
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('Perturbations are analyzed')
        call lvl_ud(1)
        
        call writo('for '//trim(i2str(n_r_sol))//' values on normal range '//&
            &trim(r2strt(min_r_sol))//'..'//trim(r2strt(max_r_sol)))
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
        if (use_pol_flux_F) then
            call writo('with toroidal mode number n = '//trim(i2str(prim_X)))
            select case(X_style)
                case (1)                                                        ! prescribed
                    call writo('and poloidal mode number m = '//&
                        &trim(i2str(min_sec_X))//'..'//trim(i2str(max_sec_X)))
                case (2)                                                        ! fast
                    call writo('and the '//trim(i2str(n_mod_X))//&
                        &' poloidal modes m that resonate most')
                case default
                    err_msg = 'No X style associated with '//&
                        &trim(i2str(X_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        else
            call writo('with poloidal mode number m = '//trim(i2str(prim_X)))
            select case(X_style)
                case (1)                                                        ! prescribed
                    call writo('and toroidal mode number n = '//&
                        &trim(i2str(min_sec_X))//'..'//trim(i2str(max_sec_X)))
                case (2)                                                        ! fast
                    call writo('and the '//trim(i2str(n_mod_X))//&
                        &' toroidal modes n that resonate most')
                case default
                    err_msg = 'No X style associated with '//&
                        &trim(i2str(X_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        
        call lvl_ud(-1)
        
        ! read PB3D output file
        ierr = read_PB3D(.false.,.true.,.false.,.false.,.true.,.false.,.false.,&
            &.false.)                                                           ! read the equilibrium grid and variables
        CHCKERR('')
        
        ! reconstructing grids depends on equilibrium style
        ! user output
        call writo('Reconstructing PB3D output on output grid')
        call lvl_ud(1)
        ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.true.,.false.,&
            &.false.,.false.,grid_eq=grid_eq,grid_eq_B=grid_eq_B,eq=eq,met=met)
        CHCKERR('')
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! the  field-aligned metric variables and  equilibrium variables
                ! are identical to output
                met_B => met
#if ldebug
                if (debug_run_driver_X) eq_B => eq
#endif
            case (2)                                                            ! HELENA
                ! also need field-aligined metrics
                allocate(met_B)
                ierr = create_met(grid_eq_B,met_B)
                CHCKERR('')
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,met=met,&
                    &met_out=met_B,eq_met=eq,grid_name='field-aligned &
                    &equilibrium grid')
                CHCKERR('')
#if ldebug
                if (debug_run_driver_X) then
                    ! also field-aligned equilibrium variables
                    allocate(eq_B)
                    ierr = create_eq(grid_eq_B,eq_B)
                    CHCKERR('')
                    if (eq_style.eq.2) then
                        ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,eq=eq,&
                            &eq_out=eq_B,grid_name='field-aligned &
                            &perturbation grid')
                        CHCKERR('')
                    end if
                end if
#endif
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! Divide perturbation grid under group processes, calculating the limits
        ! and the normal coordinate.
        allocate(r_F_X(n_r_sol))
        ierr = calc_norm_range(X_limits=X_limits,r_F_X=r_F_X)
        CHCKERR('')
        
        ! create and setup perturbation grid with division limits
        call writo('Setting up perturbation grid')
        call lvl_ud(1)
        ierr = setup_and_calc_grid_X(grid_eq,grid_X,eq,r_F_X,X_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! setup limits of modes
        ierr = setup_nm_X(grid_eq,grid_X,eq)
        CHCKERR('')
        
        ! tests
        ierr = check_X_modes(grid_eq,eq)
        CHCKERR('')
        
        ! plot resonances if requested
        if (plot_resonance .and. rank.eq.0) then
            ierr = resonance_plot(eq,grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! setup field-aligned perturbation grid
        select case (eq_style)
            case (1)                                                            ! VMEC
                grid_X_B => grid_X
            case (2)                                                            ! HELENA
                call writo('Setting up field-aligned perturbation grid')
                call lvl_ud(1)
                allocate(grid_X_B)
                ierr = setup_and_calc_grid_X(grid_eq_B,grid_X_B,eq,r_F_X,&
                    &X_limits)
                CHCKERR('')
                call lvl_ud(-1)
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        deallocate(r_F_X)
        
        ! write perturbation grid variables to output
        ierr = print_output_grid(grid_X,'perturbation',&
            &'X'//trim(rich_info_short()))
        CHCKERR('')
        if (eq_style.eq.2) then                                                 ! for HELENA also print field-aligned grid
            ierr = print_output_grid(grid_X_B,'field-aligned perturbation',&
                &'X_B'//trim(rich_info_short()))
            CHCKERR('')
        end if
        
        ! divide perturbation jobs
        ierr = divide_X_jobs(grid_X,1)
        CHCKERR('')
        
        ! main loop over vectorial jobs
        X_job_nr = 0
        X_jobs_1: do
            ! get next job
            ierr = get_next_job(X_job_nr)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call writo('Job '//trim(i2str(X_job_nr))//' is started by process '&
                &//trim(i2str(rank)))
            call lvl_ud(1)
            
            call writo('The series of '//flux_name//' mode numbers ('//&
                &trim(i2str(X_jobs_lims(1,X_job_nr)))//'..'//&
                &trim(i2str(X_jobs_lims(2,X_job_nr)))//') is calculated')
            
            ! calculate X variables, vector phase
            ierr = calc_X(grid_eq,grid_X,eq,met,X_1(1),&
                &lim_sec_X=X_jobs_lims(:,X_job_nr))
            CHCKERR('')
            
#if ldebug
            if (debug_run_driver_X) then
                call writo('To debug, the vectorial perturbation quantities &
                    &are interpolated angularly here already')
                call writo('Normally, this is done later, in the tensorial &
                    &step')
                call lvl_ud(1)
                
                ! create perturbation
                call create_X(grid_X_B,X_1_B,lim_sec_X=X_jobs_lims(:,X_job_nr))
                
                ! interpolate on field-aligned grid
                ierr = interp_HEL_on_grid(grid_X,grid_X_B,X_1=X_1(1),&
                    &X_1_out=X_1_B,grid_name='field-aligned grid')
                CHCKERR('')
                
                call lvl_ud(-1)
                
                ! write vectorial perturbation variables to output
                ierr = print_output_X(grid_X_B,X_1_B,lim_sec_X=&
                    &X_jobs_lims(:,X_job_nr))
                CHCKERR('')
                
                !!! plot information for comparison between VMEC and HELENA
                !!call plot_info_for_VMEC_HEL_comparision(grid_X_B,X_1_B)
                
                ! clean up
                call dealloc_X(X_1_B)
            else
#endif
                ! write vectorial perturbation variables to output
                ierr = print_output_X(grid_X,X_1(1),lim_sec_X=&
                    &X_jobs_lims(:,X_job_nr))
                CHCKERR('')
                
                !!! plot information for comparison between VMEC and HELENA
                !!call plot_info_for_VMEC_HEL_comparision(grid_X,X_1(1))
#if ldebug
            end if
#endif
            
            ! clean up
            call dealloc_X(X_1(1))
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed by process '&
                &//trim(i2str(rank)))
        end do X_jobs_1
        
        ! user output
        if (n_procs.gt.1) call writo('Waiting for the other processes')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        ! divide perturbation jobs, tensor phase
        ierr = divide_X_jobs(grid_X,2)
        CHCKERR('')
        
        ! main loop over tensorial jobs
        X_job_nr = 0
        X_jobs_2: do
            ! get next job
            ierr = get_next_job(X_job_nr,dim_reused)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call writo('Job '//trim(i2str(X_job_nr))//' is started by process '&
                &//trim(i2str(rank)))
            call lvl_ud(1)
            
            call writo('The block of '//flux_name//' mode numbers ('//&
                &trim(i2str(X_jobs_lims(1,X_job_nr)))//'..'//&
                &trim(i2str(X_jobs_lims(2,X_job_nr)))//')x('//&
                &trim(i2str(X_jobs_lims(3,X_job_nr)))//'..'//&
                &trim(i2str(X_jobs_lims(4,X_job_nr)))//') is calculated')
            
            ! retrieve vectorial perturbation if dimension not reused
            do id = 1,2
                if (dim_reused(id)) then
                    call writo('Vectorial perturbation variables for &
                        &dimension '//trim(i2str(id))//' reused')
                else
                    ! free the variable
                    if (allocated(X_1(id)%n)) call dealloc_X(X_1(id))
                    
                    ! user output
                    call writo('Requesting vectorial perturbation variables &
                        &for dimension '//trim(i2str(id)))
                    call lvl_ud(1)
                    
                    ! read PB3D output file for this dimension
                    ierr = read_PB3D(.false.,.false.,.false.,.false.,.false.,&
                        &.true.,.false.,.false.,lim_sec_X_1=&
                        &X_jobs_lims((id-1)*2+1:(id-1)*2+2,X_job_nr))
                    CHCKERR('')
                    
                    ! reconstruct PB3D X_1 quantities for this dimension
#if ldebug
                    if (debug_run_driver_X) then
                        ierr = reconstruct_PB3D(.false.,.false.,.false.,&
                            &.false.,.false.,.true.,.false.,.false.,&
                            &grid_X=grid_X_B,X_1=X_1(id),lim_sec_X_1=&
                            &X_jobs_lims((id-1)*2+1:(id-1)*2+2,X_job_nr))
                        CHCKERR('')
                    else
#endif
                        ierr = reconstruct_PB3D(.false.,.false.,.false.,&
                            &.false.,.false.,.true.,.false.,.false.,&
                            &grid_X=grid_X,X_1=X_1(id),lim_sec_X_1=&
                            &X_jobs_lims((id-1)*2+1:(id-1)*2+2,X_job_nr))
                        CHCKERR('')
#if ldebug
                    end if
#endif
                    
                    ! user output
                    call lvl_ud(-1)
                    call writo('Vectorial perturbation variables for &
                        &dimension '//trim(i2str(id))//' loaded')
                end if
            end do
        
#if ldebug
            if (debug_run_driver_X) then
                ! calculate X variables, tensor phase
                ierr = calc_X(grid_eq_B,grid_X_B,eq_B,met_B,X_1(1),X_1(2),X_2,&
                    &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
                CHCKERR('')
            else
#endif
                ! calculate X variables, tensor phase
                ierr = calc_X(grid_eq,grid_X,eq,met,X_1(1),X_1(2),X_2,&
                    &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
                CHCKERR('')
                
                ! adapt tensorial perturbation variables to field-aligned coords
                ! if HELENA
                if (eq_style.eq.2) then
                    ierr = interp_HEL_on_grid(grid_X,grid_X_B,X_2=X_2,&
                        &grid_name='field-aligned grid')
                    CHCKERR('')
                end if
#if ldebug
            end if
#endif
            
            ! calculate vacuum response
            ierr = calc_vac(X_2)
            CHCKERR('')
            
            ! integrate magnetic  integrals of tensorial  perturbation variables
            ! over field-aligned grid
            ierr = calc_magn_ints(grid_eq_B,grid_X_B,met_B,X_2,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! write tensorial perturbation variables to output file
            ierr = print_output_X(grid_X,X_2,lim_sec_X=&
                    &reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X_2)
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed by process '&
                &//trim(i2str(rank)))
        end do X_jobs_2
        
        ! clean up
        do id = 1,2
            call dealloc_X(X_1(id))
        end do
        
        ! cleaning up
        call writo('Cleaning up')
        call lvl_ud(1)
        call dealloc_grid(grid_eq)
        call dealloc_grid(grid_X)
        if (eq_style.eq.2) then
            call dealloc_grid(grid_eq_B)
            deallocate(grid_eq_B)
            call dealloc_grid(grid_X_B)
            deallocate(grid_X_B)
            call dealloc_met(met_B)
            deallocate(met_B)
#if ldebug
            if (debug_run_driver_X) then
                call dealloc_eq(eq_B)
                deallocate(eq_B)
            end if
#endif
        end if
        nullify(grid_eq_B)
        nullify(grid_X_B)
        nullify(met_B)
#if ldebug
            if (debug_run_driver_X) nullify(eq_B)
#endif
        call dealloc_eq(eq)
        call dealloc_met(met)
        ! deallocate variables depending on equilibrium style
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call dealloc_VMEC
            case (2)                                                            ! HELENA
                call dealloc_HEL
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('Clean')
        
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
        subroutine plot_info_for_VMEC_HEL_comparision(grid_X,X)
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
        end subroutine plot_info_for_VMEC_HEL_comparision
#endif
    end function run_driver_X
end module driver_X
