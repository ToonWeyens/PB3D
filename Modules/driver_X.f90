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
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_driver_X() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, max_it_r, rank, &
            &plot_resonance, X_job_nr, X_jobs_lims
        use MPI_utilities, only: wait_MPI
        use X_vars, only: dealloc_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X, min_r_X, max_r_X, min_n_r_X
        use grid_vars, only: dealloc_grid, &
            &n_par_X, min_par_X, max_par_X
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use PB3D_vars, only: dealloc_PB3D, &
            &PB3D_type
        use utilities, only: test_max_memory
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use X_ops, only: calc_X, check_X_modes, resonance_plot, print_output_X
        use vac, only: calc_vac
        !!use utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_X'
        
        ! local variables
        type(PB3D_type), target :: PB3D                                         ! output PB3D for which to do postprocessing
        type(PB3D_type), pointer :: PB3D_B => null()                            ! PB3D variables on a field-aligned grid
        type(X_1_type) :: X_1(2)                                                ! vectorial X variables
        type(X_2_type) :: X                                                     ! tensorial X variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: dim_reused(2)                                                ! wether dimension is reused
        integer :: id                                                           ! counter
        character(len=8) :: flux_name                                           ! name of flux variable
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
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
            call writo('for '//trim(i2str(min_n_r_X))//' values on &
                &normal range '//trim(r2strt(min_r_X))//'..'//&
                &trim(r2strt(max_r_X)))
        else
            call writo('for minimally '//trim(i2str(min_n_r_X))//' values on &
                &normal range '//trim(r2strt(min_r_X))//'..'//&
                &trim(r2strt(max_r_X)))
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
                ! the field-aligned grid is identical to the output grid
                PB3D_B => PB3D
                ! normal call to reconstruct_PB3D
                ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.false.,&
                    &PB3D)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! the field-aligned grid is different form the output grid
                allocate(PB3D_B)
                ! additionally need field-aligned equilibrium grid
                ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.false.,&
                    &PB3D,PB3D_B%grid_eq)
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! tests
        ierr = check_X_modes(PB3D%eq)
        CHCKERR('')
        
        ! plot resonances if requested
        if (plot_resonance .and. rank.eq.0) then
            ierr = resonance_plot(PB3D%eq,PB3D%grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! divide perturbation jobs
        ierr = divide_X_jobs(PB3D%grid_eq,1)
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
            ierr = calc_X(PB3D%grid_eq,PB3D%eq,PB3D%met,X_1(1),&
                &lim_sec_X=X_jobs_lims(:,X_job_nr))
            CHCKERR('')
            
            ! write vectorial perturbation variables to output
            ierr = print_output_X(PB3D%grid_eq,X_1(1))
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X_1(1))
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed by process '&
                &//trim(i2str(rank)))
        end do X_jobs_1
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        ! divide perturbation jobs, tensor phase
        ierr = divide_X_jobs(PB3D%grid_eq,2)
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
            
            write(*,*) 'TEMPORARILY SLOWING PROC NE 0 DOWDN !!!!!!!!!!!!!!!!!!!'
            call sleep(rank*10)
            
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
                    ierr = reconstruct_PB3D(.false.,.false.,.true.,.false.,&
                        &.false.,PB3D,X_1=X_1(id),lim_sec_X_1=&
                        &X_jobs_lims((id-1)*2+1:(id-1)*2+2,X_job_nr))
                    CHCKERR('')
                    
                    ! user output
                    call lvl_ud(-1)
                    call writo('Vectorial perturbation variables for &
                        &dimension '//trim(i2str(id))//' loaded')
                end if
            end do
            
            ! calculate X variables, tensor phase
            ierr = calc_X(PB3D%grid_eq,PB3D%eq,PB3D%met,X_1(1),X_1(2),X,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! calculate vacuum response
            ierr = calc_vac(X)
            CHCKERR('')
            
            ! integrate  magnetic integrals of tensorial  perturbation variables
            ! over field-aligned grid
            !ierr = calc_magn_ints(grid_eq_B,X_B)
            !CHCKERR('')
            !!!!!!!!!!!!!!!!!
            write(*,*) '!!!! CONTINUE HERE !!!!!!!!!'
            write(*,*) 'NEED TO BE ABLE TO INTEGRATE THE QUANTITIES IN FIELD-ALIGNED GRID WITHOUT COPYING TOO MUCH !!!!!'
            write(*,*) 'BEST DIRECTLY !!!'
            
            ! write tensorial perturbation variables to output file
            ierr = print_output_X(PB3D%grid_eq,X)
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X)
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed by process '&
                &//trim(i2str(rank)))
        end do X_jobs_2
        
        ! clean up
        do id = 1,2
            call dealloc_X(X_1(id))
        end do
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        ! deallocate PB3D variable
        call dealloc_PB3D(PB3D)
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
    end function run_driver_X
end module driver_X
