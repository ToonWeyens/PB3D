!------------------------------------------------------------------------------!
!> Driver of the equilibrium part of PB3D.
!------------------------------------------------------------------------------!
module driver_eq
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use vac_vars, only: vac_type
    
    implicit none
    private
    public run_driver_eq
    
    ! global variables
#if ldebug
    logical :: plot_info = .false.                                              !< plot information for comparison with HELENA \ldebug
#endif
    
contains
    !> Main driver of PB3D equilibrium part.
    !!
    !!  - sets up ([out] means for output):
    !!      * \c grid_eq [out] (for HELENA, only first Richardson level)
    !!      * \c grid_eq_B [out] (for VMEC, equal to grid_eq_out)
    !!      * \c eq_1 [out] (only first Richardson level)
    !!      * \c eq_2 [out] (for HELENA, only first Richardson level)
    !!  where output means
    !!      * on the equilibrium grid if \c X_grid style is 1,3 (no change).
    !!      * on a redistributed grid if \c X_grid_style is 2
    !!
    !!  - writes to HDF5:
    !!      * \c grid_eq (for HELENA, only first Richardson level)
    !!      * \c grid_eq_B (for VMEC, equal to grid_eq)
    !!      * \c eq_1 (only first Richardson level)
    !!      * \c eq_2 (only for HELENA)
    !!
    !!  - deallocates:
    !!      * \c grid_eq [out] before setting up
    !!      * \c grid_B_eq [out] before setting up
    !!      * \c eq_2 [out] before setting up
    !!
    !! \return ierr
    integer function run_driver_eq(grid_eq_out,grid_eq_B_out,eq_1_out,&
        &eq_2_out,vac) result(ierr)
        
        use num_vars, only: use_pol_flux_F, eq_style, plot_flux_q, &
            &plot_magn_grid, plot_B, plot_J, plot_kappa, eq_job_nr, &
            &eq_jobs_lims, jump_to_sol, rich_restart_lvl, ltest, alpha_style, &
            &X_grid_style
        use eq_ops, only: calc_eq, print_output_eq, flux_q_plot, &
            &redistribute_output_eq, B_plot, J_plot, kappa_plot
        use grid_ops, only: setup_grid_eq_B, print_output_grid, &
            &calc_norm_range, setup_grid_eq, calc_ang_grid_eq_B, &
            &magn_grid_plot, redistribute_output_grid
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_eq_1, &
            &reconstruct_PB3D_eq_2
        use grid_vars, only: min_par_X, max_par_X, n_alpha, alpha
        use num_utilities, only: derivs
        use MPI_utilities, only: wait_MPI
        use rich_vars, only: rich_lvl
        use HDF5_ops, only: create_output_HDF5
        use vac_ops, only: store_vac
        !!use num_utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_eq'
        
        ! input / output
        type(grid_type), intent(inout), target :: grid_eq_out                   !< output equilibrium grid
        type(grid_type), intent(inout), pointer :: grid_eq_B_out                !< output field-aligned equilibrium grid
        type(eq_1_type), intent(inout) :: eq_1_out                              !< flux equilibrium variables in output grid
        type(eq_2_type), intent(inout) :: eq_2_out                              !< metric equilibrium variables in output grid
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium variables
        type(eq_2_type) :: eq_2                                                 ! metric equilibrium variables
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_eq_B                                            ! field-aligned equilibrium grid
        integer :: eq_limits(2)                                                 ! min. and max. index of eq. grid of this process
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        logical :: do_eq_1_ops                                                  ! whether specific calculations for eq_1 are necessary
        logical :: do_eq_2_ops                                                  ! whether specific calculations for eq_2 are necessary
        logical :: only_half_grid                                               ! calculate only half grid
        logical :: dealloc_vars                                                 ! whether to deallocate variables to save memory
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        character(len=max_str_ln) :: grid_eq_B_name                             ! name of grid_eq_B
        
        ! initialize ierr
        ierr = 0
        
        ! decide whether to append Richardson level to name
        select case (eq_style)
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        ! Divide equilibrium grid under  group processes, calculating the limits
        ! and the normal coordinate.
        ierr = calc_norm_range('PB3D_eq',eq_limits=eq_limits)
        CHCKERR('')
        
        ! jump to solution if requested
        if (rich_lvl.eq.rich_restart_lvl .and. jump_to_sol) then
            call writo('Prepare to jump to solution')
            call lvl_ud(1)
            
            ! grid_eq_out
            select case (X_grid_style)
                case (1,3)                                                      ! equilibrium or optimized
                    ! perturbation  quantities  are  calculated  on  equilibrium
                    ! normal grid or optimized variant
                    ierr = reconstruct_PB3D_grid(grid_eq_out,'eq',&
                        &rich_lvl=rich_lvl_name,grid_limits=eq_limits)
                    CHCKERR('')
                case (2)                                                        ! solution
                    ! start from equilibrium grid in equilibrium limits and then
                    ! redistribute it
                    ierr = reconstruct_PB3D_grid(grid_eq,'eq',&
                        &rich_lvl=rich_lvl_name,grid_limits=eq_limits)
                    CHCKERR('')
                    ierr = redistribute_output_grid(grid_eq,grid_eq_out)
                    CHCKERR('')
                    call grid_eq%dealloc()
            end select
            
            ! eq_1_out
            ierr = reconstruct_PB3D_eq_1(grid_eq_out,eq_1_out,'eq_1')
            CHCKERR('')
            
            ! eq_2_out
            if (eq_style.eq.2) then
                ierr = reconstruct_PB3D_eq_2(grid_eq_out,eq_2_out,'eq_2')
                CHCKERR('')
            end if
            
            call lvl_ud(-1)
            call writo('Skipping rest to jump to solution')
            return
        end if
        
        ! set up whether to deallocate variables
        dealloc_vars = .true.
#if ldebug
        if (plot_info) dealloc_vars = .false.
        if (ltest) dealloc_vars = .false.
#endif
        if (plot_B .or. plot_J .or. plot_kappa) dealloc_vars = .false.          ! need transformation matrices
        
        ! user output
        call writo('The equilibrium variables are processed')
        call lvl_ud(1)
        
        select case (alpha_style)
            case (1)                                                            ! single field line, multiple turns
                call writo('With a single field-line')
                call lvl_ud(1)
                call writo('with label alpha = '//&
                    &trim(r2strt(alpha(1))))
                call writo('with parallel range '//trim(r2strt(min_par_X))//&
                    &'..'//trim(r2strt(max_par_X)))
                call lvl_ud(-1)
            case (2)                                                            ! multiple field lines, single turns
                call writo('With '//trim(i2str(n_alpha))//' field-lines')
        end select
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        call writo('using the '//trim(flux_name)//' flux as the normal &
            &variable')
        
        call lvl_ud(-1)
        
        ! set up whether half or full parallel grid has to be calculated
        if (rich_lvl.eq.1) then
            only_half_grid = .false.
        else
            only_half_grid = .true.
        end if
        
        ! decide whether to do certain equilibrium calculations
        if (rich_lvl.eq.rich_restart_lvl .and. eq_job_nr.eq.1) then
            do_eq_1_ops = .true.
        else
            do_eq_1_ops = .false.
        end if
        select case (eq_style)
            case (1)                                                            ! VMEC
                do_eq_2_ops = .true.
            case (2)                                                            ! HELENA
                if (rich_lvl.eq.rich_restart_lvl .and. eq_job_nr.eq.1) then
                    do_eq_2_ops = .true.
                else
                    do_eq_2_ops = .false.
                end if
        end select
        
        ! setup equilibrium grid
        call writo('Determine the equilibrium grid')
        call lvl_ud(1)
        ierr = setup_grid_eq(grid_eq,eq_limits)
        call lvl_ud(-1)
        CHCKERR('')
        
        ! Calculate the flux equilibrium quantities.
        ! Note: There  are some E  quantities that are needed  in to set  up the
        ! equilibrium grid and that are not present in eq_1_out.
        ierr = calc_eq(grid_eq,eq_1)
        CHCKERR('')
        
        ! Write to HDF5 only for first calculation
        if (rich_lvl.eq.1 .and. eq_job_nr.eq.1) then
            ! write flux equilibrium variables to output
            ierr = print_output_eq(grid_eq,eq_1,'eq_1')
            CHCKERR('')
            
            ! plot flux quantities if requested
            if (plot_flux_q) then
                ierr = flux_q_plot(grid_eq,eq_1)
                CHCKERR('')
            else
                call writo('Flux quantities plot not requested')
            end if
        end if
        
        ! grid_eq and grid_eq_B
        ! VMEC: The equilibrium grid is  field-aligned but its angular variables
        ! need to be set up still.
        ! HELENA:  The equilibrium  grid  is not  field-aligned  but taken  from
        ! the  equilibrium output.  The  field-aligned  grid has  to  be set  up
        ! separately.
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate angular grid points for equilibrium grid
                call writo('Calculate angular equilibrium grid')
                ierr = calc_ang_grid_eq_B(grid_eq,eq_1,only_half_grid=&
                    &only_half_grid)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! set up field-aligned equilibrium grid
                call writo('Determine the field-aligned equilibrium grid')
                ierr = setup_grid_eq_B(grid_eq,grid_eq_B,eq_1,&
                    &only_half_grid=only_half_grid)
                CHCKERR('')
                
                ! write field-aligned equilibrium grid variables to output
                ierr = print_output_grid(grid_eq_B,'field-aligned equilibrium',&
                    &'eq_B',rich_lvl=rich_lvl,par_div=.true.)
                CHCKERR('')
        end select
        
        ! only do when needed
        if (do_eq_2_ops) then
            ! write equilibrium grid variables to output
            ierr = print_output_grid(grid_eq,'equilibrium','eq',&
                &rich_lvl=rich_lvl_name,par_div=size(eq_jobs_lims,2).gt.1)
            CHCKERR('')
            
            ! Calculate the metric equilibrium quantities
            ierr = calc_eq(grid_eq,eq_1,eq_2,dealloc_vars=dealloc_vars)
            CHCKERR('')
            
#if ldebug
            if (plot_info) then
                ! plot info for comparison between VMEC and HELENA
                call plot_info_for_VMEC_HEL_comparision()
            end if
#endif
            
            if (rich_lvl.eq.1 .and. eq_style.eq.2) then                         ! only if first Richardson level
                ! write metric equilibrium variables to output
                ierr = print_output_eq(grid_eq,eq_2,'eq_2',par_div=.false.)
                CHCKERR('')
            end if
            
            ! store vacuum variables
            ierr = store_vac(grid_eq,eq_1,eq_2,vac)
            CHCKERR('')
        end if
        
        ! set output variables
        select case (X_grid_style)
            case (1,3)                                                          ! equilibrium or enriched
                call writo('Copy the equilibrium grids and variables to &
                    &output grid')
            case (2)                                                            ! solution
                call writo('Redistribute the equilibrium grids and variables &
                    &to output grid')
        end select
        call lvl_ud(1)
        
        ! grid_eq_out
        if (do_eq_2_ops) then
            if (associated(grid_eq_out%r_F)) call grid_eq_out%dealloc()         ! deallocate if necessary
            select case (X_grid_style)
                case (1,3)                                                      ! equilibrium or enriched
                    ierr = grid_eq%copy(grid_eq_out)
                    CHCKERR('')
                case (2)                                                        ! solution
                    ierr = redistribute_output_grid(grid_eq,grid_eq_out)
                    CHCKERR('')
            end select
        end if
        
        ! eq_1_out
        if (do_eq_1_ops) then
            if (allocated(eq_1_out%pres_FD)) call eq_1_out%dealloc()            ! deallocate if necessary
            select case (X_grid_style)
                case (1,3)                                                      ! equilibrium or optimized
                    call eq_1%copy(grid_eq,eq_1_out)
                case (2)                                                        ! solution
                    ierr = redistribute_output_eq(grid_eq,grid_eq_out,eq_1,&
                        &eq_1_out)
                    CHCKERR('')
            end select
        end if
        
        ! eq_2_out
        if (do_eq_2_ops) then
            if (allocated(eq_2_out%jac_FD)) call eq_2_out%dealloc()             ! deallocate if necessary
            select case (X_grid_style)
                case (1,3)                                                      ! equilibrium or optimized
                    call eq_2%copy(grid_eq,eq_2_out)
                case (2)                                                        ! solution
                    ierr = redistribute_output_eq(grid_eq,grid_eq_out,eq_2,&
                        &eq_2_out)
                    CHCKERR('')
            end select
        end if
        
        ! grid_eq_B_out
        select case (eq_style)
            case (1)                                                            ! VMEC
                grid_eq_B_out => grid_eq_out
            case (2)                                                            ! HELENA
                if (associated(grid_eq_B_out)) then
                    call grid_eq_B_out%dealloc()
                    deallocate(grid_eq_B_out)
                    nullify(grid_eq_B_out)
                end if
                allocate(grid_eq_B_out)
                select case (X_grid_style)
                    case (1,3)                                                  ! equilibrium or optimized
                        ierr = grid_eq_B%copy(grid_eq_B_out)
                        CHCKERR('')
                    case (2)                                                    ! solution
                        ierr = redistribute_output_grid(grid_eq_B,grid_eq_B_out)
                        CHCKERR('')
                end select
        end select
        
        ! clean up
        call grid_eq%dealloc()
        call eq_1%dealloc()                                                     ! eq_1 is always calculated because there are some quantities that are not present in eq_1_out
        if (do_eq_2_ops) then
            call eq_2%dealloc()
        end if
        if (eq_style.eq.2) then
            call grid_eq_B%dealloc()
        end if
        
        call lvl_ud(-1)
        
        ! plot magnetic field if requested
        ! (done in parts, for every parallel job)
        if (plot_B) then
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = B_plot(grid_eq_out,eq_1_out,eq_2_out,&
                        &rich_lvl=rich_lvl)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    if (rich_lvl.eq.1) then
                        ierr = B_plot(grid_eq_out,eq_1_out,eq_2_out)
                        CHCKERR('')
                    end if
            end select
        else
            if (eq_job_nr.eq.1) call writo('Magnetic field plot not requested')
        end if
        
        ! plot current if requested
        ! (done in parts, for every parallel job)
        if (plot_J) then
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = J_plot(grid_eq_out,eq_1_out,eq_2_out,&
                        &rich_lvl=rich_lvl)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    if (rich_lvl.eq.1) then
                        ierr = J_plot(grid_eq_out,eq_1_out,eq_2_out)
                        CHCKERR('')
                    end if
            end select
        else
            if (eq_job_nr.eq.1) call writo('Current plot not requested')
        end if
        
        ! plot curvature if requested
        ! (done in parts, for every parallel job)
        if (plot_kappa) then
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = kappa_plot(grid_eq_out,eq_1_out,eq_2_out,&
                        &rich_lvl=rich_lvl)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    if (rich_lvl.eq.1) then
                        ierr = kappa_plot(grid_eq_out,eq_1_out,eq_2_out)
                        CHCKERR('')
                    end if
            end select
        else
            if (eq_job_nr.eq.1) call writo('Curvature plot not requested')
        end if
        
        ! plot full field-aligned grid if requested
        ! Note: As  this needs a total  grid, it cannot be  easiliy created like
        ! plot_B, where the output is updated for every equilibrium job.
        if (plot_magn_grid) then
            if (eq_job_nr.eq.size(eq_jobs_lims,2)) then
                ! set up the name of grid_eq_B
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        grid_eq_B_name = 'eq'                                   ! already field-aligned
                    case (2)                                                    ! HELENA
                        grid_eq_B_name = 'eq_B'                                 ! not already field-aligned
                end select
                
                ! reconstruct the full field-aligned grid
                ierr = reconstruct_PB3D_grid(grid_eq_B,trim(grid_eq_B_name),&
                    &rich_lvl=rich_lvl,tot_rich=.true.,grid_limits=eq_limits)
                CHCKERR('')
                
                ! wait for all processes (sometimes necessary to finish lock)
                ierr = wait_MPI()
                CHCKERR('')
                
                ! plot it
                ierr = magn_grid_plot(grid_eq_B)
                CHCKERR('')
                
                ! deallocate
                call grid_eq_B%dealloc()
            end if
        else
            if (eq_job_nr.eq.1) call writo('Magnetic grid plot not requested')
        end if
#if ldebug
    contains
        ! Plots information for comparison between HELENA and VMEC.
        !> \private
        subroutine plot_info_for_VMEC_HEL_comparision()
            use HELENA_vars, only: R_H, Z_H
            use input_utilities, only: pause_prog, get_int, get_log
            
            ! local variables
            real(dp), allocatable :: x_plot(:,:,:), y_plot(:,:,:), z_plot(:,:,:)
            integer :: id
            integer :: d(3)
            logical :: not_ready            
            
            call writo('Plotting information for comparison between VMEC and &
                &HELENA',alert=.true.)
            not_ready = .true.
            do while (not_ready)
                call writo('derivative in dim 1?')
                d(1) =  get_int(lim_lo=0)
                call writo('derivative in dim 2?')
                d(2) =  get_int(lim_lo=0)
                call writo('derivative in dim 3?')
                d(3) =  get_int(lim_lo=0)
                
                ! x, y, z
                allocate(x_plot(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
                allocate(y_plot(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
                allocate(z_plot(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
                x_plot = grid_eq%theta_F
                y_plot = 0._dp
                do id = 1,grid_eq%loc_n_r
                    z_plot(:,:,id) = grid_eq%loc_r_F(id)/maxval(grid_eq%r_F)
                end do
                call plot_HDF5('x_plot','x_plot',x_plot)
                call plot_HDF5('y_plot','y_plot',y_plot)
                call plot_HDF5('z_plot','z_plot',z_plot)
                
                ! flux quantities
                call print_ex_2D('pres_FD','pres_FD',eq_1%pres_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_ex(['pres_FD'],'pres_FD',1,1,.false.)
                call print_ex_2D('q_saf_FD','q_saf_FD',eq_1%q_saf_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_ex(['q_saf_FD'],'q_saf_FD',1,1,.false.)
                call print_ex_2D('flux_p_FD','flux_p_FD',&
                    &eq_1%flux_p_FD(:,d(2)),x=grid_eq%loc_r_F,draw=.false.)
                call draw_ex(['flux_p_FD'],'flux_p_FD',1,1,.false.)
                call print_ex_2D('flux_t_FD','flux_t_FD',&
                    &eq_1%flux_t_FD(:,d(2)),x=grid_eq%loc_r_F,draw=.false.)
                call draw_ex(['flux_t_FD'],'flux_t_FD',1,1,.false.)
                
                ! R and Z
                select case (eq_style)
                    case(1)
                        call plot_HDF5('R_V','R_V',&
                            &eq_2%R_E(:,:,:,d(1),d(2),d(3)),x=x_plot,y=y_plot,&
                            &z=z_plot)
                        call plot_HDF5('Z_V','Z_V',&
                            &eq_2%Z_E(:,:,:,d(1),d(2),d(3)),x=x_plot,y=y_plot,&
                            &z=z_plot)
                        call plot_HDF5('L_V','L_V',&
                            &eq_2%L_E(:,:,:,d(1),d(2),d(3)),x=x_plot,y=y_plot,&
                            &z=z_plot)
                    case(2)
                        if (sum(d).eq.0) then
                            call plot_HDF5('R_H','R_H',&
                                &reshape(R_H(:,grid_eq%i_min:grid_eq%i_max),&
                                &[grid_eq%n(1:2),grid_eq%loc_n_r]),&
                                &x=x_plot,y=y_plot,z=z_plot)
                            call plot_HDF5('Z_H','Z_H',&
                                &reshape(Z_H(:,grid_eq%i_min:grid_eq%i_max),&
                                &[grid_eq%n(1:2),grid_eq%loc_n_r]),&
                                &x=x_plot,y=y_plot,z=z_plot)
                        else
                            call writo('R and Z can only be plot for d = 0')
                        end if
                end select
                
                ! metric fators
                call plot_HDF5('jac_FD','jac_FD',&
                    &eq_2%jac_FD(:,:,:,d(1),d(2),d(3)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                do id = 1,6
                    call plot_HDF5('g_FD_'//trim(i2str(id)),&
                        &'g_FD_'//trim(i2str(id)),&
                        &eq_2%g_FD(:,:,:,id,d(1),d(2),d(3)),&
                        &x=x_plot,y=y_plot,z=z_plot)
                    call plot_HDF5('h_FD_'//trim(i2str(id)),&
                        &'h_FD_'//trim(i2str(id)),&
                        &eq_2%h_FD(:,:,:,id,d(1),d(2),d(3)),&
                        &x=x_plot,y=y_plot,z=z_plot)
                end do
                
                ! clean up
                deallocate(x_plot,y_plot,z_plot)
                
                call writo('Want to redo the plotting?')
                not_ready = get_log(.true.)
            end do
            
            ! reset not_ready
            not_ready = .true.
            
            call writo('Done, paused',alert=.true.)
            call pause_prog()
        end subroutine plot_info_for_VMEC_HEL_comparision
#endif
    end function run_driver_eq
end module driver_eq
