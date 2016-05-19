!------------------------------------------------------------------------------!
!   Driver of the equilibrium part of PB3D.                                    !
!------------------------------------------------------------------------------!
module driver_eq
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    
    implicit none
    private
    public run_driver_eq
    
    ! global variables
#if ldebug
    logical :: plot_info= .false.                                               ! plot information for comparison with HELENA
#endif
    
contains
    ! Main driver of PB3D equilibrium part.
    integer function run_driver_eq() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, plot_flux_q, &
            &plot_magn_grid, eq_job_nr
        use MPI_utilities, only: wait_MPI
        use eq_ops, only: calc_eq, print_output_eq, flux_q_plot
        use sol_vars, only: alpha
        use grid_ops, only: setup_grid_eq_B, print_output_grid, &
            &calc_norm_range, setup_grid_eq, calc_ang_grid_eq_B, magn_grid_plot
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid
        use num_utilities, only: derivs
        use input_utilities, only: dealloc_in
        use rich_vars, only: rich_lvl
        !!use num_utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_eq'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        type(grid_type), target :: grid_eq                                      ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B => null()                         ! field-aligned equilibrium grid
        type(eq_1_type) :: eq_1                                                 ! equilibrium for
        type(eq_2_type) :: eq_2                                                 ! equilibrium for
        integer :: eq_limits(2)                                                 ! min. and max. index of eq. grid of this process
        logical :: only_half_grid                                               ! calculate only half grid
        logical :: dealloc_vars = .true.                                        ! whether to deallocate variables to save memory
        character(len=max_str_ln) :: grid_eq_B_name                             ! name of grid_eq_B
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        if (plot_info) dealloc_vars = .false.
#endif
        
        ! some preliminary things
        ierr = wait_MPI()
        CHCKERR('')
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('The equilibrium variables are processed')
        call lvl_ud(1)
        
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        call writo('for alpha = '//trim(r2strt(alpha)))
        call writo('using the '//trim(flux_name)//' flux as the normal &
            &variable')
        
        call lvl_ud(-1)
        
        ! Divide equilibrium grid under  group processes, calculating the limits
        ! and the normal coordinate.
        call writo('Calculate normal equilibrium ranges')
        call lvl_ud(1)
        ierr = calc_norm_range(eq_limits=eq_limits)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! set up whether half or full parallel grid has to be calculated
        if (rich_lvl.eq.1) then
            only_half_grid = .false.
        else
            only_half_grid = .true.
        end if
        
        ! setup equilibrium grid
        call writo('Determine the equilibrium grid')
        call lvl_ud(1)
        ierr = setup_grid_eq(grid_eq,eq_limits,only_half_grid=only_half_grid)
        call lvl_ud(-1)
        CHCKERR('')
        
        ! Calculate the flux equilibrium quantities
        ! Note: Only F quantities are saved, but E are needed in every level
        ierr = calc_eq(grid_eq,eq_1)
        CHCKERR('')
        
        ! write  flux equilibrium variables to  output if first level  and first
        ! equilibirum job
        if (rich_lvl.eq.1 .and. eq_job_nr.eq.1) then
            ! print output
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
        
        ! Do actions depending on equilibrium style
        ! VMEC: The  equilibrium grid is  field-aligned and therefore has  to be
        ! recalculated for every  Richardson step. The output names  make use of
        ! the Richardson level.
        ! HELENA: The equilibrium  grid is not field-aligned but  taken from the
        ! equilibrium output.  The variables concerning this  is only calculated
        ! and written for  the first Richardson level.  The variables concerning
        ! the field-aligned grid is written at every Richardson level.
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate angular grid points for equilibrium grid
                call writo('Calculate angular equilibrium grid')
                ierr = calc_ang_grid_eq_B(grid_eq,eq_1,only_half_grid=&
                    &only_half_grid)
                CHCKERR('')
                
                ! write equilibrium grid variables to output
                ierr = print_output_grid(grid_eq,'equilibrium','eq',&
                    &rich_lvl=rich_lvl,eq_job=eq_job_nr)
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
                
                ! write metric equilibrium variables to output
                ierr = print_output_eq(grid_eq,eq_2,'eq_2',rich_lvl=rich_lvl,&
                    &eq_job=eq_job_nr,dealloc_vars=dealloc_vars)
                CHCKERR('')
                
                ! clean up
                call eq_2%dealloc()
            case (2)                                                            ! HELENA
                if (rich_lvl.eq.1) then
                    ! write equilibrium grid variables to output
                    ierr = print_output_grid(grid_eq,'equilibrium','eq')
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
                    
                    ! write metric equilibrium variables to output
                    ierr = print_output_eq(grid_eq,eq_2,'eq_2',&
                        &dealloc_vars=dealloc_vars)
                    CHCKERR('')
                    
                    ! clean up
                    call eq_2%dealloc()
                end if
                
                ! set up field-aligned equilibrium grid
                call writo('Determine the field-aligned equilibrium grid')
                allocate(grid_eq_B)
                ierr = setup_grid_eq_B(grid_eq,grid_eq_B,eq_1,&
                    &only_half_grid=only_half_grid)
                CHCKERR('')
                
                ! write field-aligned equilibrium grid variables to output
                ierr = print_output_grid(grid_eq_B,'field-aligned equilibrium',&
                    &'eq_B',rich_lvl=rich_lvl,eq_job=eq_job_nr)
                CHCKERR('')
        end select
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call grid_eq%dealloc()
        call eq_1%dealloc()
        if (eq_style.eq.2) then
            call grid_eq_B%dealloc()
            deallocate(grid_eq_B)
        end if
        nullify(grid_eq_B)
        call lvl_ud(-1)
        
        ! plot full field-aligned grid if requested
        if (plot_magn_grid) then
            ! allocate
            allocate(grid_eq_B)
            
            ! set up the name of grid_eq_B and rich_lvl
            select case (eq_style)
                case (1)                                                        ! VMEC
                    grid_eq_B_name = 'eq'                                       ! already field-aligned
                    rich_lvl_name = rich_lvl                                    ! append richardson level
                case (2)                                                        ! HELENA
                    grid_eq_B_name = 'eq_B'                                     ! not already field-aligned
                    rich_lvl_name = 0                                           ! do not append name
            end select
            
            ! reconstruct the full field-aligned grid
            ierr = reconstruct_PB3D_grid(grid_eq_B,trim(grid_eq_B_name),&
                &rich_lvl=rich_lvl_name,tot_rich=.true.)
            CHCKERR('')
            
            ! plot it
            ierr = magn_grid_plot(grid_eq_B)
            CHCKERR('')
            
            ! deallocate
            call grid_eq_B%dealloc()
            deallocate(grid_eq_B)
            nullify(grid_eq_B)
        else
            call writo('Magnetic grid plot not requested')
        end if
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        call dealloc_in()
        call lvl_ud(-1)
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
#if ldebug
    contains
        subroutine plot_info_for_VMEC_HEL_comparision()
            use HELENA_vars, only: R_H, Z_H
            use input_utilities, only: pause_prog, get_int, get_log
            
            ! local variables
            real(dp), allocatable :: x_plot(:,:,:), y_plot(:,:,:), z_plot(:,:,:)
            integer :: id
            integer :: d(3)
            logical :: not_ready = .true.
            
            call writo('Plotting information for comparison between VMEC and &
                &HELENA',alert=.true.)
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
                call print_GP_2D('pres_FD','pres_FD',eq_1%pres_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('pres_FD','pres_FD','pres_FD',1,1,.false.)
                call print_GP_2D('q_saf_FD','q_saf_FD',eq_1%q_saf_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('q_saf_FD','q_saf_FD','q_saf_FD',1,1,.false.)
                call print_GP_2D('flux_p_FD','flux_p_FD',&
                    &eq_1%flux_p_FD(:,d(2)),x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('flux_p_FD','flux_p_FD','flux_p_FD',1,1,.false.)
                call print_GP_2D('flux_t_FD','flux_t_FD',&
                    &eq_1%flux_t_FD(:,d(2)),x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('flux_t_FD','flux_t_FD','flux_t_FD',1,1,.false.)
                
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
