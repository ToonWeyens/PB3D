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
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    
    implicit none
    private
    public run_driver_eq
    
contains
    ! Main driver of PB3D_eq.
    ! Note: This is skipped if Richardson restart.
    integer function run_driver_eq() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, plot_flux_q, plot_grid, &
            &rich_restart_lvl
        use MPI_utilities, only: wait_MPI
        use eq_vars, only: dealloc_eq
        use grid_vars, only: dealloc_grid
        use eq_ops, only: calc_eq, calc_derived_q, print_output_eq, flux_q_plot
        use met_ops, only: calc_met, calc_F_derivs
        use met_vars, only: dealloc_met
        use sol_vars, only: alpha
        use grid_ops, only: setup_and_calc_grid_eq_B, print_output_grid, &
            &plot_grid_real
        use grid_ops, only: setup_grid_eq, calc_ang_grid_eq
        use PB3D_ops, only: reconstruct_PB3D_in
        use utilities, only: derivs
        use files_ops, only: dealloc_in
        !!use utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_eq'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B => null()                         ! field-aligned equilibrium grid
        type(eq_type) :: eq                                                     ! equilibrium for
        type(met_type) :: met                                                   ! metric variables
        
        ! initialize ierr
        ierr = 0
        
        ! check on Richardson restart
        if (rich_restart_lvl.eq.0) then
            ! some preliminary things
            ierr = wait_MPI()
            CHCKERR('')
            ierr = reconstruct_PB3D_in()                                        ! reconstruct miscellaneous PB3D output variables
            CHCKERR('')
            
            !!! calculate auxiliary quantities for utilities
            !!call calc_aux_utilities                                             ! calculate auxiliary quantities for utilities
            
            ! user output
            call writo('The equilibrium variables are processed')
            call lvl_ud(1)
            
            if (use_pol_flux_F) then
                flux_name = 'poloidal'
            else
                flux_name = 'toroidal'
            end if
            call writo('for alpha = '//trim(r2strt(alpha*pi)))
            call writo('using the '//trim(flux_name)//' flux as the normal &
                &variable')
            
            call lvl_ud(-1)
            
            ! Calculate the equilibrium quantities
            ierr = calc_eq(grid_eq,eq)
            CHCKERR('')
            
            ! Calculate the metric quantities
            ierr = calc_met(grid_eq,eq,met)
            CHCKERR('')
            
            ! Transform E into F derivatives
#if ldebug
            ierr = calc_F_derivs(grid_eq,eq,met)
#else
            ierr = calc_F_derivs(eq,met)
#endif
            CHCKERR('')
            
            ! plot flux quantities if requested
            if (plot_flux_q) then
                ierr = flux_q_plot(grid_eq,eq)
                CHCKERR('')
            else
                call writo('Flux quantities plot not requested')
            end if
            
            ! Calculate derived metric quantities
            call calc_derived_q(grid_eq,eq,met)
            
            ! set up field-aligned equilibrium grid
            ierr = setup_and_calc_grid_eq_B(grid_eq,grid_eq_B,eq)
            CHCKERR('')
            
            ! plot grid if requested
            if (plot_grid) then
                ierr = plot_grid_real(grid_eq_B)
                CHCKERR('')
            else
                call writo('Magnetic grid plot not requested')
            end if
            
            ! write equilibrium grid variables to output
            ierr = print_output_grid(grid_eq,'equilibrium','eq')
            CHCKERR('')
            if (eq_style.eq.2) then                                             ! for HELENA also print field-aligned grid
                ierr = print_output_grid(grid_eq_B,'field-aligned equilibrium',&
                    &'eq_B')
                CHCKERR('')
            end if
            
            ! write equilibrium variables to output
            ierr = print_output_eq(grid_eq,eq,met)
            CHCKERR('')
            
            !!! plot information for comparison between VMEC and HELENA
            !!call plot_info_for_VMEC_HEL_comparision()
            
            ! cleaning up
            call writo('Cleaning up')
            call lvl_ud(1)
            ierr = dealloc_in()
            CHCKERR('')
            call dealloc_grid(grid_eq)
            call dealloc_eq(eq)
            call dealloc_met(met)
            if (eq_style.eq.2) then
                call dealloc_grid(grid_eq_B)
                deallocate(grid_eq_B)
            end if
            nullify(grid_eq_B)
            call lvl_ud(-1)
            call writo('Clean')
        else
            call writo('Skipping pre-peturbation driver to go straight to the &
                &perturbation phase')
            call lvl_ud(1)
            if (rich_restart_lvl.gt.1) call writo('Restarting at Richardson &
                &level '//trim(i2str(rich_restart_lvl)))
            call lvl_ud(-1)
        end if
        
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
            
            write(*,*) '!!!! PLOTTING INFORMATION FOR COMPARISON BETWEEN &
                &VMEC AND HELENA !!!!'
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
                call print_GP_2D('pres_FD','pres_FD',eq%pres_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('pres_FD','pres_FD','pres_FD',1,1,.false.)
                call print_GP_2D('q_saf_FD','q_saf_FD',eq%q_saf_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('q_saf_FD','q_saf_FD','q_saf_FD',1,1,.false.)
                call print_GP_2D('flux_p_FD','flux_p_FD',eq%flux_p_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('flux_p_FD','flux_p_FD','flux_p_FD',1,1,.false.)
                call print_GP_2D('flux_t_FD','flux_t_FD',eq%flux_t_FD(:,d(2)),&
                    &x=grid_eq%loc_r_F,draw=.false.)
                call draw_GP('flux_t_FD','flux_t_FD','flux_t_FD',1,1,.false.)
                
                ! R and Z
                if (sum(d).eq.0) then
                    select case (eq_style)
                        case(1)
                            call plot_HDF5('R_V','R_V',eq%R_E(:,:,:,0,0,0),&
                                &x=x_plot,y=y_plot,z=z_plot)
                            call plot_HDF5('Z_V','Z_V',eq%Z_E(:,:,:,0,0,0),&
                                &x=x_plot,y=y_plot,z=z_plot)
                        case(2)
                            call plot_HDF5('R_H','R_H',&
                                &reshape(R_H(:,grid_eq%i_min:grid_eq%i_max),&
                                &[grid_eq%n(1:2),grid_eq%loc_n_r]),&
                                &x=x_plot,y=y_plot,z=z_plot)
                            call plot_HDF5('Z_H','Z_H',&
                                &reshape(Z_H(:,grid_eq%i_min:grid_eq%i_max),&
                                &[grid_eq%n(1:2),grid_eq%loc_n_r]),&
                                &x=x_plot,y=y_plot,z=z_plot)
                    end select
                else
                    call writo('R and Z can only be plot for d = 0')
                end if
                
                ! metric fators
                call plot_HDF5('jac_FD','jac_FD',&
                    &met%jac_FD(:,:,:,d(1),d(2),d(3)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                do id = 1,6
                    call plot_HDF5('g_FD_'//trim(i2str(id)),&
                        &'g_FD_'//trim(i2str(id)),&
                        &met%g_FD(:,:,:,id,d(1),d(2),d(3)),&
                        &x=x_plot,y=y_plot,z=z_plot)
                    call plot_HDF5('h_FD_'//trim(i2str(id)),&
                        &'h_FD_'//trim(i2str(id)),&
                        &met%h_FD(:,:,:,id,d(1),d(2),d(3)),&
                        &x=x_plot,y=y_plot,z=z_plot)
                end do
                
                ! clean up
                deallocate(x_plot,y_plot,z_plot)
                
                call writo('Want to redo the plotting?')
                not_ready = get_log(.true.)
            end do
            
            write(*,*) '!!!! DONE, PAUSED !!!!'
            call pause_prog()
        end subroutine plot_info_for_VMEC_HEL_comparision
#endif
    end function run_driver_eq
end module driver_eq
