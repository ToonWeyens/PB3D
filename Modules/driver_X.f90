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
        use num_vars, only: use_pol_flux_F, eq_style, max_it_r, rank, &
            &plot_resonance, X_job_nr, X_jobs_lims, n_procs
        use MPI_utilities, only: wait_MPI
        use X_vars, only: dealloc_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X, min_r_sol, max_r_sol, &
            &min_n_r_sol
        use grid_vars, only: dealloc_grid, &
            &n_par_X, min_par_X, max_par_X
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met, create_met
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use utilities, only: test_max_memory
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use X_ops, only: calc_X, check_X_modes, resonance_plot, &
            &print_output_X, calc_magn_ints
        use vac, only: calc_vac
        use HELENA, only: interp_HEL_on_grid, dealloc_HEL
        use VMEC, only: dealloc_VMEC
        !!use utilities, only: calc_aux_utilities
#if ldebug
        use X_vars, only: create_X
        use eq_vars, only: create_eq
#endif
        
        character(*), parameter :: rout_name = 'run_driver_X'
        
        ! local variables
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B                                   ! field-aligned equilibrium grid
        type(eq_type), target :: eq                                             ! equilibrium variables
        type(eq_type), pointer :: eq_B                                          ! field-aligned equilibrium variables
        type(met_type), target :: met                                           ! metric variables
        type(met_type), pointer :: met_B                                        ! field-aligned metric variables
        type(X_1_type) :: X_1(2)                                                ! vectorial X variables
        type(X_2_type) :: X_2                                                   ! tensorial X variables
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
        
        if (max_it_r.eq.1) then
            call writo('for '//trim(i2str(min_n_r_sol))//' values on &
                &normal range '//trim(r2strt(min_r_sol))//'..'//&
                &trim(r2strt(max_r_sol)))
        else
            call writo('for minimally '//trim(i2str(min_n_r_sol))//' values on &
                &normal range '//trim(r2strt(min_r_sol))//'..'//&
                &trim(r2strt(max_r_sol)))
        end if
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
        if (use_pol_flux_F) then
            call writo('with toroidal mode number n = '//trim(i2str(min_n_X)))
            call writo('and poloidal mode number m = '//trim(i2str(min_m_X))//&
                &'..'//trim(i2str(max_m_X)))
        else
            call writo('with poloidal mode number m = '//trim(i2str(min_m_X)))
            call writo('and toroidal mode number n = '//trim(i2str(min_n_X))//&
                &'..'//trim(i2str(max_n_X)))
        end if
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        
        call lvl_ud(-1)
        
        ! read PB3D output file
        ierr = read_PB3D(.false.,.true.,.false.,.false.,.false.)                ! read the equilibrium variables
        CHCKERR('')
        
        ! reconstructing grids depends on equilibrium style
        ! user output
        call writo('Reconstructing PB3D output on output grid')
        call lvl_ud(1)
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! normal call to reconstruct_PB3D
                ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.false.,&
                    &grid_eq=grid_eq,eq=eq,met=met)
                CHCKERR('')
                ! the field-aligned  grid,  metric   variables  and  equilibrium
                ! the variables are identical to output
                grid_eq_B => grid_eq
                met_B => met
#if ldebug
                if (debug_run_driver_X) eq_B => eq
#endif
            case (2)                                                            ! HELENA
                ! additionally need field-aligned equilibrium grid
                allocate(grid_eq_B)
                ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.false.,&
                    &grid_eq=grid_eq,grid_eq_B=grid_eq_B,eq=eq,met=met)
                CHCKERR('')
                ! also need field-aligined metrics
                allocate(met_B)
                ierr = create_met(grid_eq_B,met_B)
                CHCKERR('')
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,met=met,&
                    &met_out=met_B,eq_met=eq,grid_name='field-aligned grid')
                CHCKERR('')
#if ldebug
                if (debug_run_driver_X) then
                    ! also field-aligned equilibrium variables
                    allocate(eq_B)
                    ierr = create_eq(grid_eq_B,eq_B)
                    CHCKERR('')
                    if (eq_style.eq.2) then
                        ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,eq=eq,&
                            &eq_out=eq_B,grid_name='field-aligned grid')
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
        
        ! tests
        ierr = check_X_modes(eq)
        CHCKERR('')
        
        ! plot resonances if requested
        if (plot_resonance .and. rank.eq.0) then
            ierr = resonance_plot(eq,grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! divide perturbation jobs
        ierr = divide_X_jobs(grid_eq,1)
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
            ierr = calc_X(grid_eq,eq,met,X_1(1),&
                &lim_sec_X=X_jobs_lims(:,X_job_nr))
            CHCKERR('')
            
#if ldebug
            if (debug_run_driver_X) then
                call writo('To debug, the vectorial perturbation quantities &
                    &are interpolated here already')
                call lvl_ud(1)
                
                ! create perturbation
                call create_X(grid_eq_B,X_1_B,lim_sec_X=X_jobs_lims(:,X_job_nr))
                
                ! interpolate on field-aligned grid
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,X_1=X_1(1),&
                    &X_1_out=X_1_B,grid_name='field-aligned grid')
                CHCKERR('')
                
                call lvl_ud(-1)
                
                ! write vectorial perturbation variables to output
                ierr = print_output_X(grid_eq_B,X_1_B)
                CHCKERR('')
                
                ! clean up
                call dealloc_X(X_1_B)
            else
#endif
                ! write vectorial perturbation variables to output
                ierr = print_output_X(grid_eq,X_1(1))
                CHCKERR('')
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
        ierr = divide_X_jobs(grid_eq,2)
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
            
            ! calculate vectorial perturbation if dimension not reused
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
                    ierr = read_PB3D(.false.,.false.,.true.,.false.,.false.,&
                        &lim_sec_X_1=X_jobs_lims((id-1)*2+1:(id-1)*2+2,&
                        &X_job_nr))
                    CHCKERR('')
                    
                    ! reconstruct PB3D X_1 quantities for this dimension
#if ldebug
                    if (debug_run_driver_X) then
                        ierr = reconstruct_PB3D(.false.,.false.,.true.,.false.,&
                            &.false.,grid_eq=grid_eq_B,X_1=X_1(id),lim_sec_X_1=&
                            &X_jobs_lims((id-1)*2+1:(id-1)*2+2,X_job_nr))
                        CHCKERR('')
                    else
#endif
                        ierr = reconstruct_PB3D(.false.,.false.,.true.,.false.,&
                            &.false.,grid_eq=grid_eq,X_1=X_1(id),lim_sec_X_1=&
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
                ierr = calc_X(grid_eq_B,eq_B,met_B,X_1(1),X_1(2),X_2,&
                    &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
                CHCKERR('')
            else
#endif
                ! calculate X variables, tensor phase
                ierr = calc_X(grid_eq,eq,met,X_1(1),X_1(2),X_2,&
                    &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
                CHCKERR('')
                ! adapt tensorial perturbation variables to field-aligned coords
                ! if HELENA
                if (eq_style.eq.2) then
                    ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,X_2=X_2,&
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
            call calc_magn_ints(grid_eq_B,met_B,X_2,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            
            call plot_HDF5('PV_int_0_RE','PV_int_0_RE',realpart(X_2%PV_int_0))
            call plot_HDF5('PV_int_1_RE','PV_int_1_RE',realpart(X_2%PV_int_1))
            call plot_HDF5('PV_int_2_RE','PV_int_2_RE',realpart(X_2%PV_int_2))
            call plot_HDF5('KV_int_0_RE','KV_int_0_RE',realpart(X_2%KV_int_0))
            call plot_HDF5('KV_int_1_RE','KV_int_1_RE',realpart(X_2%KV_int_1))
            call plot_HDF5('KV_int_2_RE','KV_int_2_RE',realpart(X_2%KV_int_2))
            call plot_HDF5('PV_int_0_IM','PV_int_0_IM',imagpart(X_2%PV_int_0))
            call plot_HDF5('PV_int_1_IM','PV_int_1_IM',imagpart(X_2%PV_int_1))
            call plot_HDF5('PV_int_2_IM','PV_int_2_IM',imagpart(X_2%PV_int_2))
            call plot_HDF5('KV_int_0_IM','KV_int_0_IM',imagpart(X_2%KV_int_0))
            call plot_HDF5('KV_int_1_IM','KV_int_1_IM',imagpart(X_2%KV_int_1))
            call plot_HDF5('KV_int_2_IM','KV_int_2_IM',imagpart(X_2%KV_int_2))
            
            ! write tensorial perturbation variables to output file
            ierr = print_output_X(grid_eq,X_2)
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
        
        ! user output
        if (n_procs.gt.1) call writo('Waiting for the other processes')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        ! cleaning up
        call dealloc_grid(grid_eq)
        if (eq_style.eq.2) then
            call dealloc_grid(grid_eq_B)
            deallocate(grid_eq_B)
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
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
    end function run_driver_X
end module driver_X
