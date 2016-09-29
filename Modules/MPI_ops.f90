!------------------------------------------------------------------------------!
!   Operations related to MPI                                                  !
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
#include <IO_resilience.h>
    use str_utilities
    use messages
    use MPI
    use num_vars, only: dp, max_str_ln, pi
    use grid_vars, only: grid_type
    use MPI_vars, only: lock_type
    
    implicit none
    private
    public start_MPI, stop_MPI, abort_MPI, broadcast_input_opts, &
        &sudden_stop, get_next_job, print_jobs_info
    
contains
    ! start MPI and gather information
    ! [MPI] Collective call
    integer function start_MPI() result(ierr)
        use num_vars, only: rank, n_procs, time_start
        use MPI_vars, only: HDF5_lock
        use files_utilities, only: nextunit
        use input_utilities, only: pause_prog
#if ldebug
        use grid_vars, only: n_alloc_grids, n_alloc_discs
        use eq_vars, only: n_alloc_eq_1s, n_alloc_eq_2s
        use X_vars, only: n_alloc_X_1s, n_alloc_X_2s
        use sol_vars, only: n_alloc_sols
#endif
#if lrIO
        use num_vars, only: nr_extra_tries_IO
#endif
        
        character(*), parameter :: rout_name = 'start_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! start MPI
        call MPI_init(ierr)                                                     ! initialize MPI
        CHCKERR('MPI init failed')
        call MPI_Comm_rank(MPI_Comm_world,rank,ierr)                            ! rank
        CHCKERR('MPI rank failed')
        call MPI_Comm_size(MPI_Comm_world,n_procs,ierr)                         ! nr. processes
        CHCKERR('MPI size failed')
        
        ! create HDF5 lock
        ierr = HDF5_lock%init(10)
        CHCKERR('')
        
        ! initialize time
        call system_clock(time_start)
        
#if ldebug
        ! set up allocated variable counters and number of extra IO tries
        n_alloc_discs = 0
        n_alloc_grids = 0
        n_alloc_eq_1s = 0
        n_alloc_eq_2s = 0
        n_alloc_X_1s = 0
        n_alloc_X_2s = 0
        n_alloc_sols = 0
#endif
#if lrIO
        nr_extra_tries_IO = 0
#endif
    end function start_MPI
    
    ! stop MPI
    ! [MPI] Collective call
    integer function stop_MPI() result(ierr)
        use MPI_vars, only: dealloc_lock, &
            &HDF5_lock
#if ldebug
        use num_vars, only: rank
        use grid_vars, only: n_alloc_grids, n_alloc_discs
        use eq_vars, only: n_alloc_eq_1s, n_alloc_eq_2s
        use X_vars, only: n_alloc_X_1s, n_alloc_X_2s
        use sol_vars, only: n_alloc_sols
#endif
#if lrIO
        use num_vars, only: rank, nr_extra_tries_IO
#endif
        
        character(*), parameter :: rout_name = 'stop_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Stopping MPI')
        
#if ldebug
        ! information about allocated variables
        if (n_alloc_discs.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_discs = '//trim(i2str(n_alloc_discs)))
        if (n_alloc_grids.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_grids = '//trim(i2str(n_alloc_grids)))
        if (n_alloc_eq_1s.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_eq_1s = '//trim(i2str(n_alloc_eq_1s)))
        if (n_alloc_eq_2s.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_eq_2s = '//trim(i2str(n_alloc_eq_2s)))
        if (n_alloc_X_1s.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_X_1s = '//trim(i2str(n_alloc_X_1s)))
        if (n_alloc_X_2s.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_X_2s = '//trim(i2str(n_alloc_X_2s)))
        if (n_alloc_sols.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', n_alloc_sols = '//trim(i2str(n_alloc_sols)))
#endif
#if lrIO
        if (nr_extra_tries_IO.ne.0) call writo('For rank '//trim(i2str(rank))//&
            &', nr_extra_tries_IO = '//trim(i2str(nr_extra_tries_IO)))
#endif
        
        ! deallocate HDF5 lock
        ierr = dealloc_lock(HDF5_lock)
        CHCKERR('')
        
        ! finalize
        call MPI_finalize(ierr)
        CHCKERR('MPI stop failed')
    end function stop_MPI
    
    ! abort MPI
    ! [MPI] Collective call
    integer function abort_MPI() result(ierr)
        character(*), parameter :: rout_name = 'abort_MPI'
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Abort(MPI_Comm_world,0,ierr)
        CHCKERR('MPI abort failed')
    end function abort_MPI
    
    ! Finds a  suitable next job. If  none are left, set X_job_nr  to a negative
    ! value.
    ! Optionally an array of logicals can be  passed to see whether the modes of
    ! the next job are the same as the previous job.
    integer function get_next_job(X_job_nr,same_modes) result(ierr)
        use num_vars, only: X_jobs_taken, X_jobs_lims, X_jobs_file_name, rank
        use files_utilities, only: nextunit
        use MPI_vars, only: X_jobs_lock
        use MPI_utilities, only: lock_req_acc, lock_return_acc
        
        character(*), parameter :: rout_name = 'get_next_job'
        
        ! input / output
        integer, intent(inout) :: X_job_nr                                      ! perturbation job nr.
        logical, intent(inout), optional :: same_modes(:)                       ! whether the modes are the same
        
        ! local variables
        integer :: current_job                                                  ! current job when calling this routine
        integer :: n_jobs                                                       ! nr. of jobs
        integer :: id, kd                                                       ! counters
        integer :: read_stat                                                    ! read status
        integer :: X_file_i                                                     ! X file number
        integer :: proc                                                         ! process that did a job
        integer :: X_job_done                                                   ! X_job done already
        character :: dummy_char                                                 ! first char of text
        
        ! initialize ierr
        ierr = 0
        
        ! set n_jobs
        n_jobs = size(X_jobs_taken)
        
        ! initialize same modes
        if (present(same_modes)) same_modes = .false.
        
        ! save current job and reset next job to some negative value
        current_job = X_job_nr
        X_job_nr = -1
        
        ! get access to X jobs lock
        ierr = lock_req_acc(X_jobs_lock)
        CHCKERR('')
        rIO(open(nextunit(X_file_i),STATUS='old',ACTION ='readwrite',&
            &FILE=X_jobs_file_name,IOSTAT=ierr),ierr)
        CHCKERR('Failed to open X jobs file')
        
        ! read X_jobs file and update X_jobs_taken
        rIO(read(X_file_i,*,iostat=ierr) dummy_char,ierr)
        if (ierr.eq.0 .and. dummy_char.ne.'#') ierr = 1                         ! even if good read, first char must be #
        CHCKERR('X_job file corrupt')
        read_stat = 0
        do while(read_stat.eq.0) 
            rIO2(read(X_file_i,*,iostat=read_stat) (proc, X_job_done),read_stat)
            if (read_stat.eq.0) X_jobs_taken(X_job_done) = .true.
        end do
        
        ! determine a suitable  next job by finding a job  whose columns or rows
        ! match the previous one
        do id = 1,n_jobs
            if (.not.X_jobs_taken(id) .and. current_job.gt.0) then              ! only if there was a previous job (current job > 0)
                if (matching_job(X_jobs_lims,current_job,id,same_modes)) then
                    X_job_nr = id
                    exit
                end if
            end if
        end do
        
        ! if no easy job found, look for any job
        if (X_job_nr.lt.1) then
            do id = 1,n_jobs
                if (.not.X_jobs_taken(id)) then
                    X_job_nr= id
                    exit
                end if
            end do
        end if
        
        ! if new job found, put back X_jobs_taken
        if (X_job_nr.ge.1) then
            backspace(UNIT=X_file_i)
            rIO(write(UNIT=X_file_i,FMT='(A)',IOSTAT=ierr) &
                &trim(i2str(rank))//' '//trim(i2str(X_job_nr)),ierr)
            CHCKERR('Failed to write')
        end if
        
        ! close X_jobs file and return access to lock
        close(X_file_i,iostat=ierr)
        CHCKERR('Failed to close X jobs file')
        ierr = lock_return_acc(X_jobs_lock)
        CHCKERR('')
    contains
        ! checks whether there is a  job with matching mode numbers. Optionally,
        ! a logical is  returned indicating in which dimension  the matching job
        ! was found.
        logical function matching_job(X_jobs_lims,prev_id,next_id,match_dim) &
            &result(res)
            ! input / output
            integer, intent(in) :: X_jobs_lims(:,:)                             ! limits of X jobs
            integer, intent(in) :: prev_id, next_id                             ! id's of previous and possible next jobs
            logical, intent(inout), optional :: match_dim(:)                    ! dimensions in which match was found
            
            ! local variables
            integer :: id                                                       ! counter
            
            ! initialize result
            res = .false.
            
            ! initialize same modes
            if (present(match_dim)) match_dim = .false.
            
            ! loop over all orders
            do id = 1,size(X_jobs_lims,1)/2
                if (X_jobs_lims((id-1)*2+1,next_id).eq.&
                    &X_jobs_lims((id-1)*2+1,prev_id) .and. &                    ! minimum matches
                    &X_jobs_lims((id-1)*2+2,next_id).eq.&
                    &X_jobs_lims((id-1)*2+2,prev_id)) then                      ! maximum matches
                    res = .true.
                    if (present(match_dim)) match_dim(id) = .true.
                end if
            end do
        end function matching_job
    end function get_next_job
    
    ! outputs job information from other processors.
    integer function print_jobs_info() result(ierr)
        use num_vars, only: X_jobs_file_name, rank
        use files_utilities, only: nextunit
        
        character(*), parameter :: rout_name = 'print_jobs_info'
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: read_stat                                                    ! read status
        integer :: X_file_i                                                     ! X file number
        integer :: proc                                                         ! process that did a job
        integer :: X_job_done                                                   ! X_job done already
        character :: dummy_char                                                 ! first char of text
        
        ! initialize ierr
        ierr = 0
        
        ! only if master
        if (rank.eq.0) then
            ! open X_jobs file
            rIO(open(nextunit(X_file_i),STATUS='old',ACTION='readwrite',&
                &FILE=X_jobs_file_name,IOSTAT=ierr),ierr)
            CHCKERR('Failed to open X jobs file')
            
            ! read  X_jobs  file  and  output message  for  jobs  done by  other
            ! processes
            rIO(read(X_file_i,*,iostat=ierr) dummy_char,ierr)
            if (ierr.eq.0 .and. dummy_char.ne.'#') ierr = 1                     ! even if good read, first char must be #
            CHCKERR('X_job file corrupt')
            read_stat = 0
            do while(read_stat.eq.0) 
                rIO2(read(X_file_i,*,iostat=read_stat) (proc, X_job_done),&
                    &read_stat)
                if (read_stat.eq.0 .and. proc.ne.rank) then
                    call writo('Job '//trim(i2str(X_job_done))//&
                        &' was done by process '//trim(i2str(proc)))
                    call lvl_ud(1)
#if ldebug
                    call writo('The output can be found in ¡¡¡OUTPUT SHOULD &
                        &BE CAPTURED!!!')
#endif
                    call lvl_ud(-1)
                end if
            end do
            
            ! close X_jobs file
            close(X_file_i,status='DELETE',iostat=ierr)
            CHCKERR('Failed to delete X jobs file')
        end if
    end function print_jobs_info
    
    ! Broadcasts options (e.g. user-prescribed) that  are not passed through the
    ! HDF5 output file (i.e. ltest, no_plots, ...).
    integer function broadcast_input_opts() result(ierr)
        use num_vars, only: max_str_ln, ltest, max_it_zero, rank, &
            &max_it_rich, relax_fac_HH, tol_zero, n_procs, n_sol_requested, &
            &tol_rich, max_nr_tries_HH, &
            &retain_all_sol, plot_flux_q, plot_magn_grid, no_plots, &
            &slab_plots, n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &min_theta_plot, max_theta_plot, min_zeta_plot, max_zeta_plot, &
            &swap_angles, plot_resonance, tol_SLEPC, prog_style, POST_style, &
            &minim_output, jump_to_sol, export_HEL, ex_plot_style, &
            &max_it_inv, tol_norm, max_it_slepc, &
            &max_tot_mem_per_proc, max_X_mem_per_proc, plot_size, &
            &do_execute_command_line, print_mem_usage, &
            &rich_restart_lvl, &
            &PB3D_name, PB3D_name_eq
        use grid_vars, only: min_par_X, max_par_X
        use rich_vars, only: no_guess, rich_lvl, min_n_par_X
        
        character(*), parameter :: rout_name = 'broadcast_input_opts'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set err_msg
        err_msg = 'MPI broadcast failed'
        
        if (n_procs.gt.1) then                                                  ! need to broadcast from rank 0 to other processes
            ! print message
            call writo('Broadcasting input option variables')
            call lvl_ud(1)
            
            ! variables that are sent for every program style:
            call MPI_Bcast(ltest,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(no_plots,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(minim_output,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(do_execute_command_line,1,MPI_LOGICAL,0,&
                &MPI_Comm_world,ierr)
            call MPI_Bcast(print_mem_usage,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_flux_q,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_magn_grid,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_resonance,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(relax_fac_HH,1,MPI_DOUBLE_PRECISION,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(tol_zero,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_theta_plot,1,MPI_DOUBLE_PRECISION,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_theta_plot,1,MPI_DOUBLE_PRECISION,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_zeta_plot,1,MPI_DOUBLE_PRECISION,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_zeta_plot,1,MPI_DOUBLE_PRECISION,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_size,2,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_rich,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_zero,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_nr_tries_HH,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(PB3D_name,len(PB3D_name),MPI_CHARACTER,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(PB3D_name_eq,len(PB3D_name_eq),MPI_CHARACTER,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(ex_plot_style,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            
            ! select according to program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    call MPI_Bcast(no_guess,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(jump_to_sol,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(export_HEL,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_n_par_X,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_inv,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_sol_requested,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(retain_all_sol,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_slepc,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(rich_restart_lvl,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_norm,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_par_X,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_par_X,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_tot_mem_per_proc,1,MPI_DOUBLE_PRECISION,&
                        &0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_X_mem_per_proc,1,MPI_DOUBLE_PRECISION,&
                        &0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_rich,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    if (rank.ne.0) allocate(tol_SLEPC(max_it_rich))
                    call MPI_Bcast(tol_SLEPC,max_it_rich,MPI_DOUBLE_PRECISION,&
                        &0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                case(2)                                                         ! POST
                    call MPI_Bcast(slab_plots,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(swap_angles,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_sol_plotted,4,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(rich_lvl,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(POST_style,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call lvl_ud(-1)
            call writo('Input option variables broadcasted')
        end if
    end function broadcast_input_opts
    
    ! Suddenly stops the computations, aborting MPI, etc.
    ! as a special case, if ierr = 66, no error message is printed
    subroutine sudden_stop(ierr)
        use num_vars, only: rank
        
        ! input / output
        integer, intent(in) :: ierr                                             ! error to output
        
        ! local variables
        integer :: ierr_abort                                                   ! error to output
        
        if (ierr.ne.66) then
            call writo('>> calling routine: PB3D (main) of rank '//&
                &trim(i2str(rank)),persistent=.true.)
            call writo('ERROR CODE '//trim(i2str(ierr))//&
                &'. Aborting MPI rank '//trim(i2str(rank)),&
                &persistent=.true.)
            call lvl_ud(1)
            ierr_abort = abort_MPI()
        else
            ierr_abort = stop_MPI()
        end if
        if (ierr_abort.ne.0) then
            call writo('MPI cannot abort...',persistent=.true.)
            call writo('Shutting down',persistent=.true.)
        end if
        call lvl_ud(-1)
        stop
    end subroutine
end module MPI_ops
