!------------------------------------------------------------------------------!
!   Operations related to MPI                                                  !
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use MPI
    use num_vars, only: dp, max_str_ln, pi
    use grid_vars, only: grid_type
    
    implicit none
    private
    public start_MPI, stop_MPI, abort_MPI, broadcast_input_vars, &
        &divide_X_jobs, sudden_stop, get_next_job, print_jobs_info
    
contains
    ! start MPI and gather information
    ! [MPI] Collective call
    integer function start_MPI() result(ierr)
        use num_vars, only: rank, n_procs
        
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
    end function start_MPI
    
    ! stop MPI
    ! [MPI] Collective call
    integer function stop_MPI() result(ierr)
        character(*), parameter :: rout_name = 'stop_MPI'
        
        ! initialize ierr
        ierr = 0
        
        call writo('Stopping MPI')
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
    
    ! divides the perturbation jobs
    ! All,  the (k,m)  pairs have  to be  calculated, so  the smallest  possible
    ! division of work is the calculation of one pair. For a certain mode number
    ! k, therefore, all  the pairs (k,m) can be  calculated sequentially, saving
    ! the values for this k in the process and cycling through m.
    ! This  can be  directly  scaled  up to  a  block of  k  and  m values,  and
    ! ultimately, all the values simultaneously.
    ! Therefore, the whole  load is divided into jobs depending  on the sizes of
    ! the blocks of k  and m values in memory. These jobs  start with the vector
    ! phase  by calculating  U and  DU, followed  in the  tensor phase  by their
    ! combinations in  PV and  KV. Then,  these are  integrated in  the parallel
    ! coordinate,  with  a  possible  interpolation  in  between.  Finally,  the
    ! integrated values are saved and the next jobs stats.
    ! This  function does  the  job of  dividing the  grids  setting the  global
    ! variables 'X_jobs_lims'  and 'X_jobs_taken' for  data of a  certain order,
    ! given by div_ord.  E.g. for the vector  phase, the order is 1  and for the
    ! tensorial phase it is 2.
    integer function divide_X_jobs(grid,div_ord) result(ierr)
        use num_vars, only: max_mem_per_proc, n_procs, X_jobs_lims, rank, &
            &X_jobs_file_name, X_jobs_taken, X_jobs_lock_file_name, &
            &use_pol_flux_F
        use X_vars, only: min_n_X, min_m_X, max_n_X, max_m_X
        use files_utilities, only: nextunit
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'divide_X_jobs'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid on which X vars are tabulated
        integer, intent(in) :: div_ord                                          ! division order
        
        ! local variables
        integer :: arr_size                                                     ! array size
        integer :: n_mod                                                        ! number of Fourier modes
        integer :: n_div                                                        ! factor by which to divide the total size
        real(dp) :: mem_size                                                    ! approximation of memory required for X variables
        integer :: n_mod_block                                                  ! nr. of modes in block
        integer, allocatable :: n_mod_loc(:)                                    ! number of modes per block
        character(len=max_str_ln) :: block_message                              ! message about how many blocks
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: file_i                                                       ! file number
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Dividing the perturbation jobs of order '//&
            &trim(i2str(div_ord)))
        call lvl_ud(1)
        
        ! set nr. of modes
        n_mod = (max_n_X-min_n_X+1)*(max_m_X-min_m_X+1)
        
        ! set arr_size
        arr_size = product(grid%n(1:2))*grid%loc_n_r
        
        ! calculate largest possible block of (k,m) values
        n_div = 0
        mem_size = huge(1._dp)
        do while (mem_size.gt.max_mem_per_proc)
            n_div = n_div + 1
            n_mod_block = ceiling(n_mod*1._dp/n_div)
            ierr = calc_memory(div_ord,arr_size,n_mod_block,mem_size)
            if (n_div.gt.n_mod) then
                ierr = 1
                err_msg = 'The memory limit is too low'
                CHCKERR(err_msg)
            end if
        end do
        if (n_div.gt.1) then
            block_message = 'The '//trim(i2str(n_mod))//&
                &' Fourier modes are split into '//trim(i2str(n_div))//&
                &' and '//trim(i2str(n_div**div_ord))//&
                &' jobs are done separately'
            if (n_procs.lt.n_div**div_ord) then
                block_message = trim(block_message)//', '//&
                    &trim(i2str(n_procs))//' at a time'
            else
                block_message = trim(block_message)//', but simultaneously'
            end if
        else
            block_message = 'The '//trim(i2str(n_mod))//' Fourier modes &
                &can be done without splitting them'
        end if
        call writo(block_message)
        call writo('The memory per process is estimated to be about '//&
            &trim(r2strt(mem_size))//'MB, whereas the maximum was '//&
            &trim(r2strt(max_mem_per_proc))//'MB')
        
        ! set up jobs data as illustrated below for 3 divisions, order 1:
        !   [1,2,3]
        ! or order 2:
        !   [1,4,7]
        !   [2,5,8]
        !   [3,6,9]
        ! etc.
        ! Also initialize the jobs taken to false.
        allocate(n_mod_loc(n_div))
        n_mod_loc = n_mod/n_div                                                 ! number of radial points on this processor
        n_mod_loc(1:mod(n_mod,n_div)) = n_mod_loc(1:mod(n_mod,n_div)) + 1       ! add a mode to if there is a remainder
        X_jobs_lims = calc_X_jobs_lims(n_mod_loc,div_ord)
        if (use_pol_flux_F) then
            X_jobs_lims = X_jobs_lims + min_m_X - 1                             ! scale with minimum poloidal mode number
        else
            X_jobs_lims = X_jobs_lims + min_n_X - 1                             ! scale with minimum toroidal mode number
        end if
        if (allocated(X_jobs_taken)) deallocate(X_jobs_taken)
        allocate(X_jobs_taken(n_div**div_ord))
        X_jobs_taken = .false.
        
        ! initialize file with global variable
        if (rank.eq.0) then
            open(STATUS='REPLACE',unit=nextunit(file_i),file=X_jobs_file_name,&
                &iostat=ierr)
            CHCKERR('Cannot open perturbation jobs file')
            write(file_i,*) '# Process, X job'
            close(file_i,iostat=ierr)
            CHCKERR('Failed to close file')
        end if
        
        ! delete possible lock file
        if (rank.eq.0) then
            open(unit=nextunit(file_i),file=X_jobs_lock_file_name,iostat=ierr)
            CHCKERR('Failed to open lock file')
            close(file_i,status='DELETE',iostat=ierr)
            CHCKERR('Failed to delete lock file')
        end if
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation jobs divided')
    contains
        ! Calculate memory in MB necessary for X variables of a certain order
        !   - order 1: 4x n_par_X x n_geo x loc_n_r x n_mod
        !   - order 2: 2x n_par_X x n_geo x loc_n_r x n_mod^2
        !              4x n_par_X x n_geo x loc_n_r x n_mod^2
        !   - higher order: not used
        ! where n_par_X  x n_geo x  loc_n_r should  be passed as  'arr_size' and
        ! n_mod as well
        integer function calc_memory(ord,arr_size,n_mod,mem_size) result(ierr)
            use ISO_C_BINDING
            character(*), parameter :: rout_name = 'calc_memory'
            
            ! input / output
            integer, intent(in) :: ord                                          ! order of data
            integer, intent(in) :: arr_size                                     ! size of part of X array
            integer, intent(in) :: n_mod                                        ! number of modes
            real(dp), intent(inout) :: mem_size                                 ! total size
            
            ! local variables
            integer(C_SIZE_T) :: dp_size                                        ! size of dp
            real(dp), parameter :: mem_scale_fac = 1.5                          ! scale factor of memory (because only estimation)
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            call lvl_ud(1)
            
            ! get size of complex variable
            dp_size = 2*sizeof(1._dp)                                           ! complex variable
            
            ! calculate memory depending on order
            select case(ord)
                case (1)                                                        ! vectorial data: U, DU
                    ! set memory size
                    mem_size = 4*arr_size
                    mem_size = mem_size*n_mod**ord
                    mem_size = mem_size*dp_size
                case (2)                                                        ! tensorial data: PV, KV
                    ! set memory size
                    mem_size = 6*arr_size
                    mem_size = mem_size*n_mod**ord
                    mem_size = mem_size*dp_size
                case default
                    ierr = 1
                    err_msg = 'Orders > 2 are not implemented'
                    CHCKERR(err_msg)
            end select
            
            ! convert B to MB
            mem_size = mem_size*1.E-6_dp
            
            ! scale memory to account for rough estimation
            mem_size = mem_size*mem_scale_fac
            
            ! test overflow
            if (mem_size.lt.0) then
                ierr = 1
                CHCKERR('Overflow occured')
            end if
            
            call lvl_ud(-1)
        end function calc_memory
        
        ! Calculate X_jobs_lims.
        recursive function calc_X_jobs_lims(n_mod,ord) result(res)
            ! input / output
            integer, intent(in) :: n_mod(:)                                     ! X jobs data
            integer, intent(in) :: ord                                          ! order of data
            integer, allocatable :: res(:,:)                                    ! result
            
            ! local variables
            integer, allocatable :: res_loc(:,:)                                ! local result
            integer :: n_div                                                    ! nr. of divisions of modes
            integer :: id                                                       ! counter
            
            ! set up nr. of divisions
            n_div = size(n_mod)
            
            ! (re)allocate result
            if (allocated(res)) deallocate(res)
            allocate(res(2*ord,n_div**ord))
            
            ! loop over divisions
            do id = 1,n_div
                if (ord.gt.1) then                                              ! call lower order
                    res_loc = calc_X_jobs_lims(n_mod,ord-1)
                    res(1:2*ord-2,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = res_loc
                end if                                                          ! set for own order
                res(2*ord-1,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                    &sum(n_mod(1:id-1))+1
                res(2*ord,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                    &sum(n_mod(1:id))
            end do
        end function calc_X_jobs_lims
    end function divide_X_jobs
    
    ! Finds a  suitable next job. If  none are left, set X_job_nr  to a negative
    ! value.
    ! Optionally an array of logicals can be  passed to see whether the modes of
    ! the next job are the same as the previous job.
    integer function get_next_job(X_job_nr,same_modes) result(ierr)
        use num_vars, only: X_jobs_taken, X_jobs_lims, X_jobs_file_name, rank, &
            &X_jobs_lock_file_name
        use files_utilities, only: nextunit, wait_file
        
        character(*), parameter :: rout_name = 'get_next_job'
        
        ! input / output
        integer, intent(inout) :: X_job_nr                                      ! perturbation job nr.
        logical, intent(inout), optional :: same_modes(:)                       ! whether the modes are the same
        
        ! local variables
        integer :: current_job                                                  ! current job when calling this routine
        integer :: n_jobs                                                       ! nr. of jobs
        integer :: id                                                           ! counter
        integer :: read_stat                                                    ! read status
        integer :: X_file_i                                                     ! X file number
        integer :: lock_file_i                                                  ! lock file number
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
        
        ! open X_jobs file if no lock-file exists
        call wait_file(lock_file_i,X_jobs_lock_file_name)
        open(STATUS='OLD',ACTION = 'READWRITE',unit=nextunit(X_file_i),&
            &file=X_jobs_file_name,iostat=ierr)
        CHCKERR('Failed to open X jobs file')
        
        ! read X_jobs file and update X_jobs_taken
        read(X_file_i,*,iostat=ierr) dummy_char
        if (ierr.eq.0 .and. dummy_char.ne.'#') ierr = 1                         ! even if good read, first char must be #
        CHCKERR('X_job file corrupt')
        read_stat = 0
        do while(read_stat.eq.0) 
            read(X_file_i,*,iostat=read_stat) proc, X_job_done
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
            write(X_file_i,*) rank, X_job_nr
        end if
        
        ! close X_jobs file and lock file
        close(X_file_i,iostat=ierr)
        CHCKERR('Failed to close X jobs file')
        close(lock_file_i,status='DELETE',iostat=ierr)
        CHCKERR('Failed to delete lock file')
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
            open(STATUS='OLD',ACTION = 'READWRITE',unit=nextunit(X_file_i),&
                &file=X_jobs_file_name,iostat=ierr)
            CHCKERR('Failed to open X jobs file')
            
            ! read  X_jobs  file  and  output message  for  jobs  done by  other
            ! processes
            read(X_file_i,*,iostat=ierr) dummy_char
            if (ierr.eq.0 .and. dummy_char.ne.'#') ierr = 1                     ! even if good read, first char must be #
            CHCKERR('X_job file corrupt')
            read_stat = 0
            do while(read_stat.eq.0) 
                read(X_file_i,*,iostat=read_stat) proc, X_job_done
                if (read_stat.eq.0 .and. proc.ne.rank) then
                    call writo('Job '//trim(i2str(X_job_done))//&
                        &' was done by process '//trim(i2str(proc)))
                    call lvl_ud(1)
                    call writo('The output can be found in ¡¡¡OUTPUT SHOULD &
                        &BE CAPTURED!!!')
                    call lvl_ud(-1)
                end if
            end do
            
            ! close X_jobs file
            close(X_file_i,status='DELETE',iostat=ierr)
            CHCKERR('Failed to delete X jobs file')
        end if
    end function print_jobs_info
    
    ! Broadcasts all  the relevant variable that have been  determined so far in
    ! the master process using the inputs to the other processes
    ! [MPI] Collective call
    integer function broadcast_input_vars() result(ierr)
        use num_vars, only: max_str_ln, ltest, EV_style, max_it_NR, &
            &max_it_rich, tol_NR, rank, n_procs, n_sol_requested, &
            &nyq_fac, tol_rich, use_pol_flux_F, use_pol_flux_E, retain_all_sol, &
            &plot_flux_q, plot_grid, no_plots, eq_style, use_normalization, &
            &n_sol_plotted, n_theta_plot, n_zeta_plot, plot_resonance, EV_BC, &
            &tol_SLEPC, rho_style, prog_style, max_it_inv, norm_disc_prec_X, &
            &norm_disc_prec_eq, norm_disc_prec_sol, BC_style, tol_norm, &
            &max_it_slepc, max_mem_per_proc, PB3D_name, norm_style, U_style, &
            &plot_size, test_max_mem
        use VMEC, only: mpol, ntor, lasym, lfreeb, nfp, rot_t_V, gam, R_V_c, &
            &R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, flux_t_V, Dflux_t_V, pres_V
        use HELENA, only: pres_H, qs, flux_p_H, nchi, chi_H, ias, h_H_11, &
            &h_H_12, h_H_33, RBphi, R_H, Z_H
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X, min_r_sol, &
            &max_r_sol
        use grid_vars, only: alpha, n_r_eq, n_par_X, min_par_X, max_par_X
        use HDF5_vars, only: var_1D_type
        use PB3D_vars, only: PB3D_version
        use eq_vars, only: max_flux_p_E, max_flux_t_E, max_flux_p_F, &
            &max_flux_t_F
        use rich, only: min_n_r_sol, no_guess, rich_lvl
#if ldebug
        use VMEC, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif
        
        character(*), parameter :: rout_name = 'broadcast_input_vars'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set err_msg
        err_msg = 'MPI broadcast failed'
        
        if (n_procs.gt.1) then                                                  ! need to broadcast from rank 0 to other processes
            ! print message
            call writo('Broadcasting variables')
            call lvl_ud(1)
            
            ! broadcast eq_style to determine what else to broadcast
            call MPI_Bcast(eq_style,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            
            ! variables that are sent for every program style:
            call MPI_Bcast(ltest,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(use_normalization,1,MPI_LOGICAL,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(no_plots,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(test_max_mem,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(tol_NR,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_r_sol,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_r_sol,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(vac_perm,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(R_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(pres_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(B_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(psi_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(rho_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(T_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(rho_style,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(norm_style,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(U_style,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_flux_q,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_grid,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_resonance,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(use_pol_flux_F,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(use_pol_flux_E,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_size,2,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(norm_disc_prec_eq,1,MPI_INTEGER,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(norm_disc_prec_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(norm_disc_prec_sol,1,MPI_INTEGER,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_m_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_m_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_n_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_n_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_rich,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_NR,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            
            ! select according to program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    call MPI_Bcast(no_guess,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_par_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_r_eq,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(nyq_fac,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_n_r_sol,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_inv,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(BC_style,2,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(EV_style,1,MPI_INTEGER,0,MPI_Comm_world,&
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
                    call MPI_Bcast(gam,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_mem_per_proc,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_rich,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(EV_BC,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_SLEPC,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                case(2)                                                         ! POST
                    call MPI_Bcast(n_sol_plotted,4,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(rich_lvl,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(PB3D_version,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_flux_p_E,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_flux_t_E,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_flux_p_F,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_flux_t_F,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(PB3D_name,len(PB3D_name),MPI_CHARACTER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! For  specific variables, choose  which equilibrium style  is being
            ! used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! variables that are sent for every program style:
                    call MPI_Bcast(lasym,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(lfreeb,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(mpol,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(ntor,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(nfp,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    
                    ! select according to program style
                    select case (prog_style)
                        case(1)                                                 ! PB3D
                            call bcast_size_1_R(flux_t_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(flux_t_V,size(flux_t_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(Dflux_t_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Dflux_t_V,size(Dflux_t_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(rot_t_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(rot_t_V,size(rot_t_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(pres_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(pres_V,size(pres_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(R_V_c,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(R_V_c,size(R_V_c),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(R_V_s,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(R_V_s,size(R_V_s),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(Z_V_c,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Z_V_c,size(Z_V_c),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(Z_V_s,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Z_V_s,size(Z_V_s),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(L_V_c,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(L_V_c,size(L_V_c),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(L_V_s,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(L_V_s,size(L_V_s),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
#if ldebug
                            if (ltest) then                                     ! ltest has already been broadcast
                                call bcast_size_4_R(B_V_sub_c,1)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_sub_c,size(B_V_sub_c),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_4_R(B_V_sub_s,1)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_sub_s,size(B_V_sub_s),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(B_V_c)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_c,size(B_V_c),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(B_V_s)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_s,size(B_V_s),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(jac_V_c)
                                CHCKERR(err_msg)
                                call MPI_Bcast(jac_V_c,size(jac_V_c),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(jac_V_s)
                                CHCKERR(err_msg)
                                call MPI_Bcast(jac_V_s,size(jac_V_s),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                            end if
#endif
                        case(2)                                                 ! POST
                            ! do nothing extra
                        case default
                            err_msg = 'No program style associated with '//&
                                &trim(i2str(prog_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                case (2)                                                        ! HELENA
                    ! variables that are sent for every program style:
                    call MPI_Bcast(ias,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(nchi,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    
                    ! select according to program style
                    select case (prog_style)
                        case(1)                                                 ! PB3D
                            call bcast_size_2_R(R_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(R_H,size(R_H),MPI_DOUBLE_PRECISION,&
                                &0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(Z_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Z_H,size(Z_H),MPI_DOUBLE_PRECISION,&
                                &0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(chi_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(chi_H,size(chi_H),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(flux_p_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(flux_p_H,size(flux_p_H),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(qs)
                            CHCKERR(err_msg)
                            call MPI_Bcast(qs,size(qs),MPI_DOUBLE_PRECISION,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(pres_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(pres_H,size(pres_H),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(RBphi)
                            CHCKERR(err_msg)
                            call MPI_Bcast(RBphi,size(RBphi),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(h_H_11)
                            CHCKERR(err_msg)
                            call MPI_Bcast(h_H_11,size(h_H_11),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(h_H_12)
                            CHCKERR(err_msg)
                            call MPI_Bcast(h_H_12,size(h_H_12),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(h_H_33)
                            CHCKERR(err_msg)
                            call MPI_Bcast(h_H_33,size(h_H_33),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                        case(2)                                                 ! POST
                            ! do nothing extra
                        case default
                            err_msg = 'No program style associated with '//&
                                &trim(i2str(prog_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call lvl_ud(-1)
            call writo('Variables broadcasted')
        end if
    contains
        ! broadcasts the size of an array, so this array can be allocated in the
        ! slave processes
        ! The index of this array is (1:)
        subroutine bcast_size_1_I(arr)                                          ! version with 1 integer argument
            ! input / output
            integer, intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_I
        ! The index of this array is (1:)
        subroutine bcast_size_1_R(arr)                                          ! version with 1 real argument
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_R
        ! The index of this array is (1:)
        subroutine bcast_size_1_var_1D(arr)                                     ! version with 1D var argument (see HDF5_ops)
            ! input / output
            type(var_1D_type), intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_var_1D
        ! The array index is (1:,1:)
        subroutine bcast_size_2_R(arr)                                          ! version with 2 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:)
            
            ! local variables
            integer :: arr_size(2)                                              ! sent ahead so arrays can be allocated
            
            if (rank.eq.0) arr_size = [size(arr,1),size(arr,2)]
            
            call MPI_Bcast(arr_size,2,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (rank.ne.0) allocate(arr(1:arr_size(1),1:arr_size(2)))
        end subroutine bcast_size_2_R
        ! The array index is (0:mpol-1,-ntor:ntor,1:,start_id)
        subroutine bcast_size_4_R(arr,start_id)                                 ! version with 4 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:,:,:)
            integer, intent(in) :: start_id
            
            ! local variables
            integer :: arr_size(4)                                              ! sent ahead so arrays can be allocated
            
            if (rank.eq.0) arr_size = [size(arr,1),size(arr,2),&
                &size(arr,3),size(arr,4)]
            
            call MPI_Bcast(arr_size,4,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (rank.ne.0) allocate(arr(0:arr_size(1)-1,&
                &-(arr_size(2)-1)/2:(arr_size(2)-1)/2,1:arr_size(3),&
                &start_id:arr_size(4)+start_id-1))
        end subroutine bcast_size_4_R
        ! The array index is (0:mpol-1,-ntor:ntor,1:)
        subroutine bcast_size_3_R(arr)                                          ! version with 3 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:,:)
            
            ! local variables
            integer :: arr_size(3)                                              ! sent ahead so arrays can be allocated
            
            if (rank.eq.0) arr_size = [size(arr,1),size(arr,2),size(arr,3)]
            
            call MPI_Bcast(arr_size,3,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (rank.ne.0) allocate(arr(0:arr_size(1)-1,&
                &-(arr_size(2)-1)/2:(arr_size(2)-1)/2,1:arr_size(3)))
        end subroutine bcast_size_3_R
    end function broadcast_input_vars
    
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
