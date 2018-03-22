!------------------------------------------------------------------------------!
!> Operations related to MPI.
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use MPI
    use num_vars, only: dp, max_str_ln, pi
    use MPI_vars, only: lock_type
    
    implicit none
    private
    public start_MPI, stop_MPI, abort_MPI, broadcast_input_opts, &
        &sudden_stop
    
contains
    !> Start MPI and gather information.
    !!
    !! \return ierr
    integer function start_MPI() result(ierr)
        use num_vars, only: rank, n_procs, time_start
        use MPI_vars, only: HDF5_lock
        use input_utilities, only: pause_prog
#if ldebug
        use grid_vars, only: n_alloc_grids, n_alloc_discs
        use eq_vars, only: n_alloc_eq_1s, n_alloc_eq_2s
        use X_vars, only: n_alloc_X_1s, n_alloc_X_2s
        use sol_vars, only: n_alloc_sols
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
    end function start_MPI
    
    !> Stop MPI.
    !!
    !! Also deallocates:
    !!  - \c grid_eq
    !!  - \c grid_eq_B
    !!  - \c grid_X
    !!  - \c grid_X_B
    !!  - \c grid_sol
    !!  - \c eq_1
    !!  - \c eq_2
    !!  - \c X_1
    !!  - \c X_2
    !!  - \c sol
    !!
    !! \return ierr
    integer function stop_MPI(grid_eq,grid_eq_B,grid_X,grid_X_B,grid_sol,eq_1,&
        &eq_2,X_1,X_2,vac,sol) result(ierr)
        
        use MPI_vars, only: dealloc_lock, &
            &HDF5_lock
        use grid_vars, only: grid_type
        use eq_vars, only: eq_1_type, eq_2_type
        use X_vars, only: X_1_type, X_2_type
        use vac_vars, only: vac_type
        use sol_vars, only: sol_type
        use num_vars, only: eq_style
#if ldebug
        use MPI_utilities, only: get_ser_var, &
            &n_waits
        use num_vars, only: rank, n_procs
        use grid_vars, only: n_alloc_grids, n_alloc_discs
        use eq_vars, only: n_alloc_eq_1s, n_alloc_eq_2s
        use X_vars, only: n_alloc_X_1s, n_alloc_X_2s
        use sol_vars, only: n_alloc_sols
#endif
        
        character(*), parameter :: rout_name = 'stop_MPI'
        
        ! input / output
        type(grid_type), intent(inout), optional :: grid_eq                     !< equilibrium grid
        type(grid_type), intent(inout), pointer, optional :: grid_eq_B          !< field-aligned equilibrium grid
        type(grid_type), intent(inout), optional :: grid_X                      !< perturbation grid
        type(grid_type), intent(inout), pointer, optional :: grid_X_B           !< field-aligned perturbation grid
        type(grid_type), intent(inout), optional :: grid_sol                    !< solution grid
        type(eq_1_type), intent(inout), optional :: eq_1                        !< Flux equilibrium variables
        type(eq_2_type), intent(inout), optional :: eq_2                        !< metric equilibrium variables
        type(X_1_type), intent(inout), optional :: X_1                          !< vectorial perturbation variables
        type(X_2_type), intent(inout), optional :: X_2                          !< integrated tensorial perturbation variables
        type(vac_type), intent(inout), optional :: vac                          !< vacuum variables
        type(sol_type), intent(inout), optional :: sol                          !< solution variables
        
        ! local variables
#if ldebug
        integer :: id                                                           ! counter
        integer, allocatable :: ser_n_waits(:)                                  ! serial n_waits
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Stopping MPI')
        
        ! clean up variables
        if (present(grid_eq)) then
            if (associated(grid_eq%r_F)) call grid_eq%dealloc()
        end if
        if (present(grid_eq_B)) then
            if (associated(grid_eq_B)) then
                if (eq_style.eq.2) then
                    call grid_eq_B%dealloc()
                    deallocate(grid_eq_B)
                end if
                nullify(grid_eq_B)
            end if
        end if
        if (present(grid_X)) then
            if (associated(grid_X%r_F)) call grid_X%dealloc()
        end if
        if (present(grid_X_B)) then
            if (associated(grid_X_B)) then
                if (eq_style.eq.2) then
                    call grid_X_B%dealloc()
                    deallocate(grid_X_B)
                end if
                nullify(grid_X_B)
            end if
        end if
        if (present(grid_sol)) then
            if (associated(grid_sol%r_F)) call grid_sol%dealloc()
        end if
        if (present(eq_1)) then
            if (allocated(eq_1%pres_FD)) call eq_1%dealloc()
        end if
        if (present(eq_2))then
            if (allocated(eq_2%jac_FD)) call eq_2%dealloc()
        end if
        if (present(X_1)) then
            if (allocated(X_1%U_0)) call X_1%dealloc()
        end if
        if (present(X_2)) then
            if (allocated(X_2%PV_0)) call X_2%dealloc()
        end if
        if (present(vac)) then
            if (allocated(vac%res)) call vac%dealloc()
        end if
        if (present(sol)) then
            if (allocated(sol%val)) call sol%dealloc()
        end if
        
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
        
        ! information about synchronization of wait_MPI
        ierr = get_ser_var([n_waits],ser_n_waits) 
        CHCKERR('')
        if (rank.eq.0 .and. (minval(ser_n_waits).ne.n_waits .or. &
            &maxval(ser_n_waits).ne.n_waits)) then
            call writo('The number of wait_MPI commands was not consistent')
            call lvl_ud(1)
            do id = 1,n_procs
                call writo('For rank '//trim(i2str(id-1))//': '//&
                    &trim(i2str(ser_n_waits(id))))
            end do
            call lvl_ud(-1)
        end if
#endif
        
        ! deallocate HDF5 lock
        ierr = dealloc_lock(HDF5_lock)
        CHCKERR('')
        
        ! finalize
        call blacs_exit(1)
        call MPI_finalize(ierr)
        CHCKERR('MPI stop failed')
    end function stop_MPI
    
    !> Abort MPI suddenly.
    !!
    !! \return ierr
    integer function abort_MPI() result(ierr)
        character(*), parameter :: rout_name = 'abort_MPI'
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Abort(MPI_Comm_world,0,ierr)
        CHCKERR('MPI abort failed')
    end function abort_MPI
    
    !> Broadcasts options (e.g. user-prescribed) that are not passed through the
    !! HDF5 output file (i.e. \c ltest, \c no_plots, ...).
    !!
    !! \see read_input_opts()
    !!
    !! \note Some  variables (e.g. \c  eq_style, ...)  are not passed  over MPI.
    !! Every  process should  call  its own  reconstruct_pb3d_in()  in order  to
    !! obtain them.
    !!
    !! \return ierr
    integer function broadcast_input_opts() result(ierr)
        use num_vars, only: max_str_ln, ltest, max_it_zero, rank, &
            &max_it_rich, relax_fac_HH, tol_zero, n_procs, n_sol_requested, &
            &tol_rich, max_nr_backtracks_HH, sol_n_procs, &
            &retain_all_sol, plot_flux_q, plot_magn_grid, plot_B, plot_J, &
            &plot_kappa, plot_sol_xi, plot_sol_Q, plot_E_rec, no_plots, &
            &plot_grid_style, n_sol_plotted, n_theta_plot, n_zeta_plot, &
            &min_theta_plot, max_theta_plot, min_zeta_plot, max_zeta_plot, &
            &min_r_plot, max_r_plot, swap_angles, plot_resonance, tol_SLEPC, &
            &prog_style, POST_style, jump_to_sol, export_HEL, compare_tor_pos, &
            &plot_VMEC_modes, EV_guess, ex_plot_style, solver_SLEPC_style, &
            &pert_mult_factor_POST, POST_output_full, POST_output_sol, &
            &tol_norm, max_it_slepc, plot_vac_pot, plot_size, &
            &max_tot_mem, max_X_mem, &
            &do_execute_command_line, print_mem_usage, &
            &rich_restart_lvl, &
            &PB3D_name, &
            &RZ_0, &
            &min_Rvac_plot, max_Rvac_plot, min_Zvac_plot, max_Zvac_plot, &
            &n_vac_plot
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
            call MPI_Bcast(do_execute_command_line,1,MPI_LOGICAL,0,&
                &MPI_Comm_world,ierr)
            call MPI_Bcast(print_mem_usage,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_flux_q,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_kappa,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_magn_grid,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_B,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_J,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
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
            call MPI_Bcast(min_r_plot,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_r_plot,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_tot_mem,1,MPI_DOUBLE_PRECISION,0,&
                &MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(sol_n_procs,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_rich,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_it_zero,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_nr_backtracks_HH,1,MPI_INTEGER,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(ex_plot_style,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_size,2,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(PB3D_name,len(PB3D_name),MPI_CHARACTER,0,&
                &MPI_Comm_world,ierr)
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
                    call MPI_Bcast(plot_VMEC_modes,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_n_par_X,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(solver_SLEPC_style,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
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
                    call MPI_Bcast(EV_guess,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_norm,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_par_X,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_par_X,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_X_mem,1,MPI_DOUBLE_PRECISION,&
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
                    call MPI_Bcast(swap_angles,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(compare_tor_pos,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(plot_sol_xi,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(plot_sol_Q,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(plot_vac_pot,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(plot_E_rec,1,MPI_LOGICAL,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(POST_output_full,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(POST_output_sol,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_vac_plot,2,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(plot_grid_style,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_sol_plotted,4,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(rich_lvl,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(POST_style,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(pert_mult_factor_POST,1,&
                        &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(RZ_0,2,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_Rvac_plot,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_Rvac_plot,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_Zvac_plot,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_Zvac_plot,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
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
    
    !> Suddenly stops the computations, aborting MPI, etc.
    !!
    !! As a special case, if \c ierr = 66, no error message is printed.
    !!
    !! \return ierr
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
