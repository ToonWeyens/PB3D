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
        &divide_X_jobs, sudden_stop, get_next_job
    
contains
    ! start MPI and gather information
    ! [MPI] Collective call
    integer function start_MPI() result(ierr)
        use num_vars, only: rank, n_procs, plt_rank
        
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
        
        ! so set plot rank to rank
        plt_rank = rank
    end function start_MPI
    
    ! stop MPI
    ! [MPI] Collective call
    integer function stop_MPI() result(ierr)
#if lold_MPI
        use num_vars, only: next_job_win
#else
        use num_vars, only: jobs_taken_win
#endif
        
        character(*), parameter :: rout_name = 'stop_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! free window
#if lold_MPI
        call MPI_Win_free(next_job_win,ierr)
#else
        call MPI_Win_free(jobs_taken_win,ierr)
#endif
        CHCKERR('Failed to free the window to jobs_taken')
        
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
    
    !! Finds a  suitable next job. If  none are left, set next_job  to a negative
    !! value
    ! [MPI] Collective call
    integer function get_next_job(X_job_nr) result(ierr)
#if lold_MPI
        use num_vars, only: next_job_win
#else
        use num_vars, only: jobs_taken_win, jobs_taken, jobs_data
#endif
        
        character(*), parameter :: rout_name = 'get_next_job'
        
        ! input / output
        integer, intent(inout) :: X_job_nr                                      ! perturbation job nr.
        
        ! local variables
        integer(kind=MPI_ADDRESS_KIND) :: null_disp = 0                         ! zero target displacement
#if lold_MPI
#else
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: current_job                                                  ! current job when calling this routine
        integer :: n_jobs                                                       ! nr. of jobs
        integer :: id                                                           ! counter
        integer :: get_req                                                      ! get request
        integer :: stat(MPI_STATUS_SIZE)                                        ! MPI status
#endif
        
        ! initialize ierr
        ierr = 0
        
#if lold_MPI
        call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,next_job_win,ierr)
        CHCKERR('Failed to get lock window on master')
        
        call MPI_get(X_job_nr,1,MPI_INTEGER,0,null_disp,1,MPI_INTEGER,&
            &next_job_win,ierr)
        CHCKERR('Failed to get next_job')
        
        call MPI_accumulate(1,1,MPI_INTEGER,0,null_disp,1,MPI_INTEGER,&
            &MPI_SUM,next_job_win,ierr)
        CHCKERR('Failed to increment')
        
        call MPI_Win_unlock(0,next_job_win,ierr)
        CHCKERR('Failed to unlock window')
#else
        ! set n_jobs
        n_jobs = size(jobs_taken)
        
        ! save current job and reset next job to some negative value
        current_job = X_job_nr
        X_job_nr = -1
        
        ! get jobs_taken window
        call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,jobs_taken_win,ierr)
        CHCKERR('Failed to lock window on master')
        call MPI_Rget(jobs_taken,n_jobs,MPI_INTEGER,0,null_disp,n_jobs,&
            &MPI_INTEGER,jobs_taken_win,get_req,ierr)
        CHCKERR('Failed to increment jobs_taken on master')
        call MPI_Wait(get_req,stat,ierr)
        CHCKERR('Failed to wait')
        
        ! determine a suitable  next job by finding a job  whose columns or rows
        ! match the previous one
        do id = 1,n_jobs
            if (jobs_taken(id).eq.0) then
                if (jobs_data(1,id).eq.jobs_data(1,current_job) .and. &
                    &jobs_data(2,id).eq.jobs_data(2,current_job) .or. &         ! columns match
                    &jobs_data(3,id).eq.jobs_data(4,current_job) .and. &
                    &jobs_data(4,id).eq.jobs_data(3,current_job)) then          ! rows match
                    X_job_nr = id
                end if
            end if
        end do
        
        ! if no easy job found, look for any job
        if (X_job_nr.lt.1) then
            do id = 1,n_jobs
                if (jobs_taken(id).eq.0) X_job_nr= id
            end do
        end if
        
        ! if new job found, put back jobs_taken
        if (X_job_nr.ge.1) then
            jobs_taken(X_job_nr) = 1
            call MPI_PUT(jobs_taken,n_jobs,MPI_INTEGER,0,null_disp,&
                &n_jobs,MPI_INTEGER,jobs_taken_win,ierr)
            err_msg = 'Failed to put back the updated jobs_taken'
            CHCKERR(err_msg)
        end if
        
        ! unlock window
        call MPI_Win_unlock(0,jobs_taken_win,ierr)
        CHCKERR('Failed to unlock window on master')
#endif
    end function get_next_job
    
    ! divides the perturbation jobs
    ! In the  end, all the  (k,m) pairs have to  be calculated, so  the smallest
    ! possible  division  of  work  is  the  calculation  of  one  pair.  For  a
    ! certain mode  number k, therefore, all  the pairs (k,m) can  be calculated
    ! sequentially, saving  the values  for this  k in  the process  and cycling
    ! through m.
    ! This  can be  directly  scaled  up to  a  block of  k  and  m values,  and
    ! ultimately, all the values simultaneously.
    ! Therefore, the whole  load is divided into jobs depending  on the sizes of
    ! the blocks of k and m values  in memory. These jobs start by calculating U
    ! and  DU, followed  by their  combinations in  PV and  KV. Then,  these are
    ! integrated in  the parallel coordinate,  with a possible  interpolation in
    ! between. Finally, the integrated values are saved and the next jobs stats.
    integer function divide_X_jobs(grid,n_mod) result(ierr)
        use num_vars, only: max_mem_per_proc, n_procs, jobs_data, rank
#if lold_MPI
        use num_vars, only: next_job_win, next_job
#else
        use num_vars, only: jobs_taken_win, jobs_taken
#endif
        use utilities, only: calc_memory
        
        character(*), parameter :: rout_name = 'divide_X_jobs'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid on which X vars are tabulated
        integer, intent(in) :: n_mod                                            ! number of Fourier modes
        
        ! local variables
        integer :: arr_size                                                     ! array size
        integer :: n_div                                                        ! factor by which to divide the total size
        real(dp) :: mem_size                                                    ! approximation of memory required for X variables
        integer :: n_mod_block                                                  ! nr. of modes in block
        integer, allocatable :: n_mod_loc(:)                                    ! number of modes per block
        character(len=max_str_ln) :: block_message                              ! message about how many blocks
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd                                                       ! counters
        integer(kind=MPI_ADDRESS_KIND) :: intlb, intsize                        ! lower bound and  size of integer
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Dividing the perturbation jobs')
        call lvl_ud(1)
        
        ! set arr_size
        arr_size = product(grid%n(1:2))*grid%loc_n_r
        
        ! calculate largest possible block of (k,m) values
        n_div = 0
        mem_size = huge(1._dp)
        do while (mem_size.gt.max_mem_per_proc)
            n_div = n_div + 1
            n_mod_block = ceiling(n_mod*1._dp/n_div)
            mem_size = calc_memory(arr_size,n_mod_block,block_mem=.true.)
            if (n_div.gt.n_mod) then
                ierr = 1
                err_msg = 'The memory limit is too low - use less grid points &
                    &or aument ''max_mem_per_proc'''
                CHCKERR(err_msg)
            end if
        end do
        if (n_div.gt.1) then
            block_message = 'The '//trim(i2str(n_mod))//&
                &' Fourier modes are split into '//trim(i2str(n_div))//&
                &' and '//trim(i2str(n_div**2))//' jobs are done separately'
            if (n_procs.lt.n_div) then
                block_message = trim(block_message)//', '//&
                    &trim(i2str(n_procs))//' at the same time'
            else
                block_message = trim(block_message)//', but simultaneously'
            end if
        else
            block_message = 'The whole range of Fourier modes can be done &
                &without dividing it'
        end if
        call writo(block_message)
        call writo('The memory per process is estimated to be about '//&
            &trim(r2strt(mem_size))//'MB, whereas the maximum was '//&
            &trim(r2strt(max_mem_per_proc))//'MB')
        
        ! set up jobs as illustrated below for 3 divisions, 9 jobs
        !   [1,4,7]
        !   [2,5,8]
        !   [3,6,9]
        allocate(jobs_data(4,n_div**2))
        allocate(n_mod_loc(n_div**2))
        n_mod_loc = n_mod/n_div                                                 ! number of radial points on this processor
        n_mod_loc(1:mod(n_mod,n_div)) = n_mod_loc(1:mod(n_mod,n_div)) + 1       ! add a mode to if there is a remainder
        do jd = 1,n_div
            do id = 1,n_div
                jobs_data(:,(jd-1)*n_div+id) = &
                    &[sum(n_mod_loc(1:id-1))+1, &                               ! k_min
                    &sum(n_mod_loc(1:id)), &                                    ! k_max
                    &sum(n_mod_loc(1:jd-1))+1, &                                ! m_min
                    &sum(n_mod_loc(1:jd))]                                      ! k_max
            end do
        end do
        
        ! initialize global variable on master and set window
#if lold_MPI
        call MPI_Type_get_extent(MPI_INTEGER,intlb,intsize,ierr)
        err_msg = 'Couldn''t determine the extent of an integer'
        CHCKERR(err_msg)
        if (rank.eq.0) then                                                     ! master
            call MPI_Win_create(next_job,1*intsize,intsize,MPI_INFO_NULL,&
                &MPI_Comm_world,next_job_win,ierr)
        else
            next_job = 0                                                        ! variable only has meaning on master
            call MPI_Win_create(next_job,0*intsize,intsize,MPI_INFO_NULL,&
                &MPI_Comm_world,next_job_win,ierr)
        end if
        CHCKERR('Couldn''t create window to next_job')
        
        ! set fence
        call MPI_Win_fence(0,next_job_win,ierr) 
        CHCKERR('Couldn''t set fence') 
#else
        allocate(jobs_taken(n_div**2))
        call MPI_Type_get_extent(MPI_INTEGER,intlb,intsize,ierr)
        err_msg = 'Couldn''t determine the extent of an integer'
        CHCKERR(err_msg)
        if (rank.eq.0) then                                                     ! master
            jobs_taken = 0                                                      ! no jobs taken yet
            call MPI_Win_create(jobs_taken,size(jobs_taken)*intsize,&
                &intsize,MPI_INFO_NULL,MPI_Comm_world,jobs_taken_win,&
                &ierr)
        else
            call MPI_Win_create(jobs_taken,0,1,MPI_INFO_NULL,&
                &MPI_Comm_world,jobs_taken_win,ierr)
        end if
        CHCKERR('Couldn''t create window to jobs_taken')
        
        ! set fence
        call MPI_Win_fence(0,jobs_taken_win,ierr) 
        CHCKERR('Couldn''t set fence') 
#endif
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation jobs divided')
    end function divide_X_jobs
    
    ! Broadcasts all  the relevant variable that have been  determined so far in
    ! the master process using the inputs to the other processes
    ! [MPI] Collective call
    integer function broadcast_input_vars() result(ierr)
        use num_vars, only: max_str_ln, ltest, EV_style, max_it_NR, max_it_r, &
            &tol_NR, rank, n_procs, no_guess, n_sol_requested, nyq_fac, tol_r, &
            &use_pol_flux_F, use_pol_flux_E, retain_all_sol, plot_flux_q, &
            &plot_grid, no_plots, eq_style, use_normalization, n_sol_plotted, &
            &n_theta_plot, n_zeta_plot, plot_resonance, EV_BC, tol_slepc, &
            &rho_style, prog_style, max_it_inv, norm_disc_prec_X, &
            &norm_disc_prec_eq, norm_disc_prec_sol, BC_style, tol_norm_r, &
            &max_it_slepc, max_mem_per_proc
        use VMEC, only: mpol, ntor, lasym, lfreeb, nfp, rot_t_V, gam, R_V_c, &
            &R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, flux_t_V, Dflux_t_V, pres_V
        use HELENA, only: pres_H, qs, flux_p_H, nchi, chi_H, ias, h_H_11, &
            &h_H_12, h_H_33, RBphi, R_H, Z_H
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X, min_n_r_X, &
            &min_r_X, max_r_X
        use grid_vars, only: alpha, n_r_eq, n_par_X, min_par_X, max_par_X
        use HDF5_ops, only: var_1D
        use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B, vars_1D_X, vars_1D_sol
#if ldebug
        use VMEC, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif
        
        character(*), parameter :: rout_name = 'broadcast_input_vars'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id                                                           ! counter
        
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
            call MPI_Bcast(use_normalization,1,MPI_LOGICAL,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_NR,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(tol_NR,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(vac_perm,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(rho_style,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_flux_q,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_grid,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_resonance,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            
            ! select according to program style
            select case (prog_style)
                case(1)                                                         ! PB3D pre-perturbation
                    call MPI_Bcast(ltest,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(use_pol_flux_F,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(use_pol_flux_E,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(no_guess,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(no_plots,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_par_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_r_eq,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_norm_r,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(nyq_fac,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(norm_disc_prec_eq,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,&
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
                    call MPI_Bcast(R_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(pres_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(B_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(psi_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(rho_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(T_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                case(2)                                                         ! PB3D perturbation
                    call MPI_Bcast(max_mem_per_proc,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_n_r_X,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_r,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_inv,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(BC_style,2,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_r,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(EV_style,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_sol_requested,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(EV_BC,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_slepc,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(retain_all_sol,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(norm_disc_prec_X,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(norm_disc_prec_sol,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_slepc,1,MPI_INTEGER,0,MPI_Comm_world,&
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
                    call bcast_size_1_var_1D(vars_1D_eq)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_eq)
                        call bcast_size_1_R(vars_1D_eq(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%p,size(vars_1D_eq(id)%p),&
                            &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%tot_i_min,&
                            &size(vars_1D_eq(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%tot_i_max,&
                            &size(vars_1D_eq(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                    call bcast_size_1_var_1D(vars_1D_eq_B)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_eq_B)
                        call bcast_size_1_R(vars_1D_eq_B(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%p,&
                            &size(vars_1D_eq_B(id)%p),MPI_DOUBLE_PRECISION,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq_B(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%tot_i_min,&
                            &size(vars_1D_eq_B(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq_B(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%tot_i_max,&
                            &size(vars_1D_eq_B(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                case(3)                                                         ! PB3D_POST
                    call MPI_Bcast(n_sol_plotted,4,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(norm_disc_prec_sol,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call bcast_size_1_var_1D(vars_1D_eq)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_eq)
                        call bcast_size_1_R(vars_1D_eq(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%p,size(vars_1D_eq(id)%p),&
                            &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%tot_i_min,&
                            &size(vars_1D_eq(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%tot_i_max,&
                            &size(vars_1D_eq(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                    call bcast_size_1_var_1D(vars_1D_eq_B)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_eq_B)
                        call bcast_size_1_R(vars_1D_eq_B(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%p,&
                            &size(vars_1D_eq_B(id)%p),MPI_DOUBLE_PRECISION,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq_B(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%tot_i_min,&
                            &size(vars_1D_eq_B(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq_B(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%tot_i_max,&
                            &size(vars_1D_eq_B(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                    call bcast_size_1_var_1D(vars_1D_X)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_X)
                        call bcast_size_1_R(vars_1D_X(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%p,size(vars_1D_X(id)%p),&
                            &MPI_DOUBLE_PRECISION,&
                            &0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_X(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%tot_i_min,&
                            &size(vars_1D_X(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_X(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%tot_i_max,&
                            &size(vars_1D_X(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                    call bcast_size_1_var_1D(vars_1D_sol)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_sol)
                        call bcast_size_1_R(vars_1D_sol(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_sol(id)%p,&
                            &size(vars_1D_sol(id)%p),MPI_DOUBLE_PRECISION,&
                            &0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_sol(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_sol(id)%tot_i_min,&
                            &size(vars_1D_sol(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_sol(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_sol(id)%tot_i_max,&
                            &size(vars_1D_sol(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_sol(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
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
                    ! select according to program style
                    select case (prog_style)
                        case(1)                                                 ! PB3D pre-perturbation
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
                            call MPI_Bcast(lasym,1,MPI_LOGICAL,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(lfreeb,1,MPI_LOGICAL,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(mpol,1,MPI_INTEGER,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(ntor,1,MPI_INTEGER,0,MPI_Comm_world,&
                                &ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(nfp,1,MPI_INTEGER,0,MPI_Comm_world,&
                                    &ierr)
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
                        case(2)                                                 ! PB3D perturbation
                        case(3)                                                 ! PB3D_POST
                            ! do nothing extra
                        case default
                            err_msg = 'No program style associated with '//&
                                &trim(i2str(prog_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                case (2)                                                        ! HELENA
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
                            call MPI_Bcast(nchi,1,MPI_INTEGER,0,MPI_Comm_world,&
                                &ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(ias,1,MPI_INTEGER,0,MPI_Comm_world,&
                                &ierr)
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
                        case(2)                                                 ! PB3D_POST
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
            type(var_1D), intent(inout), allocatable :: arr(:)
            
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
