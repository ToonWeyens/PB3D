!------------------------------------------------------------------------------!
!> Driver of the perturbation part of PB3D.
!------------------------------------------------------------------------------!
module driver_X
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type, modes_type, &
        &mds_X
    
    implicit none
    private
    public run_driver_X
    
    ! global variables
    integer :: X_limits(2)                                                      !< min. and max. index of X grid for this process
    real(dp), allocatable :: r_F_X(:)                                           !< normal points in perturbation grid
    character(len=8) :: flux_name(2)                                            !< name of flux variable
    character(len=1) :: mode_name(2)                                            !< name of modes
    integer :: rich_lvl_name                                                    !< either the Richardson level or zero
#if ldebug
    logical :: debug_run_driver_X_1 = .false.                                   !< debug information for run_driver_X_1 \ldebug
    logical :: debug_run_driver_X_2 = .true.                                   !< debug information for run_driver_X_2 \ldebug
    logical :: plot_info= .false.                                               !< plot information for comparison with HELENA \ldebug
#endif
    
contains
    !> Main driver of PB3D perturbation part.
    !!
    !!  - sets up:
    !!      * [0] \c grid_X (for HELENA, only first Richardson level)
    !!      * [0] \c grid_X_B (for VMEC, equal to grid_X_out)
    !!      * [1] \c X_1 
    !!      * [2] \c X_2 
    !!  - writes to HDF5:
    !!      * [0] \c grid_X (for HELENA, only first Richardson level)
    !!      * [0] \c grid_X_B (for VMEC, equal to grid_X)
    !!      * [1] \c X_1 (only for HELENA, only first Richardson level)
    !!      * [2] \c X_2 
    !!  - deallocates:
    !!      * \c X_1 before setting up
    !!      * \c grid_X before setting up
    !!      * \c grid_X_B before setting up
    !!
    !!  ([x] indicates driver x)
    !!
    !! \note  \c eq_2  needs  to be  intent(inout) because  interp_HEL_on_grid()
    !! requires  this for  generality.  The  variable is  not  modified in  this
    !! driver, though.
    !!
    !! \return ierr
    integer function run_driver_X(grid_eq,grid_eq_B,grid_X,grid_X_B,eq_1,eq_2,&
        &X_1,X_2) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, plot_resonance, &
            &X_style, rich_restart_lvl, jump_to_sol, eq_job_nr, X_grid_style
        use rich_vars, only: n_par_X, rich_lvl
        use PB3D_ops, only: reconstruct_PB3D_grid, reconstruct_PB3D_X_1, &
            &reconstruct_PB3D_X_2
        use X_vars, only: min_sec_X, max_sec_X, prim_X, min_r_sol, max_r_sol, &
            &n_mod_X
        use grid_ops, only: calc_norm_range
        use grid_vars, only: min_par_X, max_par_X, n_r_sol, n_r_eq
        use X_ops, only: calc_X, check_X_modes, resonance_plot, &
            &init_modes, setup_modes
        use files_utilities, only: delete_file
        !!use num_utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_X'
        
        ! input / output
        type(grid_type), intent(in), target :: grid_eq                          !< equilibrium grid (should be in but needs inout for interp_HEL_on_grid)
        type(grid_type), intent(in), pointer :: grid_eq_B                       !< field-aligned equilibrium grid (should be in but needs inout for interp_HEL_on_grid)
        type(grid_type), intent(inout), target :: grid_X                        !< perturbation grid
        type(grid_type), intent(inout), pointer :: grid_X_B                     !< field-aligned perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(inout), target :: eq_2                          !< metric equilibrium variables (should be in but needs inout for interp_HEL_on_grid)
        type(X_1_type), intent(inout) :: X_1                                    !< vectorial perturbation variables
        type(X_2_type), intent(inout) :: X_2                                    !< tensorial perturbation variables
        
        ! local variables
        type(eq_2_type), pointer :: eq_2_B                                      ! field-aligned metric equilibrium variables
        
        ! initialize ierr
        ierr = 0
        
        ! set up whether Richardson level has to be appended to the name
        select case (eq_style) 
            case (1)                                                            ! VMEC
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append
        end select
        
        ! Divide perturbation grid under group processes, calculating the limits
        ! and the normal coordinate.
        if (rich_lvl.eq.rich_restart_lvl .and. eq_job_nr.eq.1) then
            select case (X_grid_style)
                case (1)                                                        ! equilibrium
                    allocate(r_F_X(n_r_eq))
                    ierr = calc_norm_range(eq_limits=X_limits)
                    CHCKERR('')
                    r_F_X = grid_eq%r_F
                case (2)                                                        ! solution
                    allocate(r_F_X(n_r_sol))
                    ierr = calc_norm_range(sol_limits=X_limits,r_F_sol=r_F_X)
                    CHCKERR('')
            end select
        end if
        
        ! jump to solution if requested
        if (rich_lvl.eq.rich_restart_lvl .and. jump_to_sol) then
            call writo('Prepare to jump to solution')
            call lvl_ud(1)
            
            ! grid_X
            ierr = reconstruct_PB3D_grid(grid_X,'X',rich_lvl=rich_lvl_name,&
                &grid_limits=X_limits)
            CHCKERR('')
            
            ! grid_X_B
            select case (eq_style)
                case (1)                                                        ! VMEC
                    grid_X_B => grid_X
                case (2)                                                        ! HELENA
                    allocate(grid_X_B)
                    ierr = reconstruct_PB3D_grid(grid_X_B,'X_B',&
                        &rich_lvl=rich_lvl)
                    CHCKERR('')
            end select
            
            ! initialize modes and set up
            ierr = init_modes(grid_eq,eq_1)
            CHCKERR('')
            ierr = setup_modes(mds_X,grid_eq,grid_X,plot_nm=.false.)
            CHCKERR('')
            
            ! X_1
            if (eq_style.eq.2) then                                             ! HELENA
                ierr = reconstruct_PB3D_X_1(mds_X,grid_X,X_1,'X_1')
                CHCKERR('')
            end if
            
            ! X_2
            ierr = reconstruct_PB3D_X_2(mds_X,grid_X,X_2,'X_2_int',&
                &rich_lvl=rich_lvl,is_field_averaged=.true.)
            CHCKERR('')
            
            call lvl_ud(-1)
            call writo('Skipping rest to jump to solution')
            return
        end if
        
        ! user output
        call writo('Perturbations are analyzed')
        call lvl_ud(1)
        
        call writo('for '//trim(i2str(size(r_F_X)))//&
            &' values on normal range '//trim(r2strt(min_r_sol))//'..'//&
            &trim(r2strt(max_r_sol)))
        if (X_grid_style.eq.1) then
            call lvl_ud(1)
            call writo('interpolation to the solution grid with '//&
                &trim(i2str(n_r_sol))//' will happen in the solution driver')
            call lvl_ud(-1)
        end if
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
        
        ! set flux and mode names
        if (use_pol_flux_F) then
            flux_name = ['poloidal','toroidal']
            mode_name = ['m','n']
        else
            flux_name = ['toroidal','poloidal']
            mode_name = ['n','m']
        end if
        
        ! user output
        call writo('with '//flux_name(2)//' mode number '//mode_name(2)//&
            &' = '//trim(i2str(prim_X)))
        select case(X_style)
            case (1)                                                            ! prescribed
                call writo('and '//flux_name(1)//' mode number '//&
                    &mode_name(1)//' = '//trim(i2str(min_sec_X))//'..'//&
                    &trim(i2str(max_sec_X)))
            case (2)                                                            ! fast
                call writo('and the '//trim(i2str(n_mod_X))//&
                    &' '//flux_name(1)//' modes '//mode_name(1)//&
                    &' that resonate most')
        end select
        
        call lvl_ud(-1)
        
        ! Run grid part of driver
        ierr = run_driver_X_0(grid_eq,grid_eq_B,grid_X,grid_X_B,eq_1,&
            &eq_2,eq_2_B)
        CHCKERR('')
        
        ! setup nm and check the X modes if first Richardson level
        if (rich_lvl.eq.1 .and. eq_job_nr.eq.1) then
            ! initialize modes and set up
            ierr = init_modes(grid_eq,eq_1)
            CHCKERR('')
            ierr = setup_modes(mds_X,grid_eq,grid_X,plot_nm=.true.)
            CHCKERR('')
            
            ! tests
            ierr = check_X_modes(grid_eq,eq_1)
            CHCKERR('')
        end if
        
        ! jump to solution if  requested (assuming that integrated X_2 variables
        ! are present in the HDF5 file)
        if (rich_lvl.eq.rich_restart_lvl .and.  jump_to_sol) then
            call writo('Skipping rest to jump to solution')
        else
            ! plot resonances if requested
            if (plot_resonance .and. rich_lvl.eq.1 .and. eq_job_nr.eq.1) then
                ierr = resonance_plot(mds_X,grid_eq,eq_1)
                CHCKERR('')
            else
                call writo('Resonance plot not requested')
            end if
            
            ! Run vectorial part of driver
            ierr = run_driver_X_1(grid_eq,grid_X,eq_1,eq_2,X_1)
            CHCKERR('')
            
            ! Run Tensorial part of driver
            ierr = run_driver_X_2(grid_eq_B,grid_X,grid_X_B,eq_1,eq_2_B,X_1,X_2)
            CHCKERR('')
        end if
        
        ! clean up
        call writo('Clean up')
        call lvl_ud(1)
        if (eq_style.eq.2) then
            call eq_2_B%dealloc()
            deallocate(eq_2_B)
        end if
        nullify(eq_2_B)
        call lvl_ud(-1)
    end function run_driver_X
    
    !> part  0 of  driver_x:  perturbation  grid as  well  as reconstruction  of
    !! variables.
    integer function run_driver_X_0(grid_eq,grid_eq_B,grid_X,grid_X_B,eq_1,&
        &eq_2,eq_2_B) result(ierr)
        use num_vars, only: eq_style, eq_job_nr, rich_restart_lvl, eq_jobs_lims
        use grid_ops, only: setup_grid_X, print_output_grid
        use rich_vars, only: rich_lvl
        use HELENA_ops, only: interp_HEL_on_grid
        
        character(*), parameter :: rout_name = 'run_driver_X_0'
        
        ! input / output
        type(grid_type), intent(in), target :: grid_eq                          !< equilibrium grid
        type(grid_type), intent(in), pointer :: grid_eq_B                       !< field-aligned equilibrium grid
        type(grid_type), intent(inout), target :: grid_X                        !< perturbation grid
        type(grid_type), intent(inout), pointer :: grid_X_B                     !< field-aligned perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(inout), target :: eq_2                          !< metric equilibrium variables
        type(eq_2_type), intent(inout), pointer :: eq_2_B                       !< field-aligned metric equilibrium variables
        
        ! local variables
        logical :: do_X_2_ops                                                   ! whether specific calculations for X_2 are necessary
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Setting up perturbation grid')
        call lvl_ud(1)
        
        ! decide whether to do certain  perturbation calculations and whether to
        ! append Richardson level to name
        select case (eq_style)
            case (1)                                                            ! VMEC
                do_X_2_ops = .true.
                rich_lvl_name = rich_lvl                                        ! append richardson level
            case (2)                                                            ! HELENA
                if (rich_lvl.eq.rich_restart_lvl .and. eq_job_nr.eq.1) then
                    do_X_2_ops = .true.
                else
                    do_X_2_ops = .false.
                end if
                rich_lvl_name = 0                                               ! do not append name
        end select
        
        call writo('Reconstruct equilibrium variables')
        call lvl_ud(1)
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! point: grid is already field-aligned
                eq_2_B => eq_2
            case (2)                                                            ! HELENA
                ! allocate field-aligned equilibrium variables
                allocate(eq_2_B)
                
                ! interpolate field-aligned metric equilibrium variables
                call eq_2_B%init(grid_eq_B,setup_E=.false.,setup_F=.true.)
                ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,eq_2=eq_2,&
                    &eq_2_out=eq_2_B,eq_1=eq_1,grid_name='field-aligned &
                    &equilibrium grid')
                CHCKERR('')
        end select
        call lvl_ud(-1)
        
        ! create and setup perturbation grid with division limits
        if (do_X_2_ops) then
            call writo('Determine the perturbation grid')
            call lvl_ud(1)
            if (associated(grid_X%r_F)) call grid_X%dealloc()                   ! deallocate if necessary
            ierr = setup_grid_X(grid_eq,grid_X,r_F_X,X_limits)
            CHCKERR('')
            call lvl_ud(-1)
            
            ! print output
            ierr = print_output_grid(grid_X,'perturbation','X',&
                &rich_lvl=rich_lvl_name,par_div=size(eq_jobs_lims,2).gt.1)
            CHCKERR('')
        end if
        
        ! grid_X_B
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! do nothing: grid is already field-aligned
                grid_X_B => grid_X
            case (2)                                                            ! HELENA
                ! set up field-aligned perturbation grid
                if (associated(grid_X_B)) then
                    call grid_X_B%dealloc()
                    deallocate(grid_X_B)
                    nullify(grid_X_B)
                end if
                allocate(grid_X_B)
                ierr = setup_grid_X(grid_eq_B,grid_X_B,r_F_X,X_limits)
                CHCKERR('')
                
                ! print field-aligned output
                ierr = print_output_grid(grid_X_B,'field-aligned &
                    &perturbation','X_B',rich_lvl=rich_lvl,par_div=.true.)
                CHCKERR('')
        end select
        
        call lvl_ud(-1)
        call writo('Perturbation grid set up')
    end function run_driver_X_0
    
    !> Part 1 of driver_x: Vectorial jobs.
    !!
    !! \note  Everything is  done here  in the  original grids,  not necessarily
    !! field-aligned.
    !!
    !! \return ierr
    integer function run_driver_X_1(grid_eq,grid_X,eq_1,eq_2,X_1) result(ierr)
        use num_vars, only: rank, n_procs, eq_style, U_style, eq_job_nr, &
            &rich_restart_lvl
        use X_ops, only: calc_X, print_output_X
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'run_driver_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        type(X_1_type), intent(inout) :: X_1                                    !< vectorial X variables
        
        ! local variables
        logical :: do_X_1_ops                                                   ! whether specific calculations for X_1 are necessary
#if ldebug
        character(len=max_str_ln), allocatable :: var_names(:)                  ! names of variables
        character(len=max_str_ln) :: file_name                                  ! name of file
        integer :: ld                                                           ! counter
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculating U up to order '//trim(i2str(U_style)))
        
        ! decide whether to do certain  perturbation calculations and whether to
        ! append Richardson level to name
        select case (eq_style)
            case (1)                                                            ! VMEC
                do_X_1_ops = .true.
            case (2)                                                            ! HELENA
                if (rich_lvl.eq.rich_restart_lvl .and. eq_job_nr.eq.1) then
                    do_X_1_ops = .true.
                else
                    do_X_1_ops = .false.
                end if
        end select
        
        if (do_X_1_ops) then
            ! possibly deallocate
            if (allocated(X_1%U_0)) call X_1%dealloc()
            
            ! calc variables
            ierr = calc_X(mds_X,grid_eq,grid_X,eq_1,eq_2,X_1)
            CHCKERR('')
        else
            ! For HELENA, everything has been done
            return
        end if
        
        ! write vectorial perturbation variables to output if HELENA
        ! Note: There should be only 1 parallel job
        if (eq_style.eq.2) then                                                 ! only if first Richardson level
            ierr = print_output_X(grid_X,X_1,'X_1',par_div=.false.)
            CHCKERR('')
        end if
        
#if ldebug
        ! write vectorial perturbation quantities to output
        if (debug_run_driver_X_1) then
            ! angles
            call plot_HDF5('theta_F_B','theta_F',grid_X%theta_F)
            call plot_HDF5('zeta_F_B','zeta_F',grid_X%zeta_F)
            ! U_0
            allocate(var_names(size(X_1%U_0,4)))
            if (n_procs.eq.1) then
                file_name = 'U_0'
            else
                file_name = 'U_0'//trim(i2str(rank))
            end if
            do ld = 1,size(var_names)
                var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
            end do
            call plot_HDF5(var_names,'RE_'//trim(file_name),&
                &rp(X_1%U_0),col_id=4,col=1)
            call plot_HDF5(var_names,'IM_'//trim(file_name),&
                &ip(X_1%U_0),col_id=4,col=1)
            deallocate(var_names)
            ! U_1
            allocate(var_names(size(X_1%U_1,4)))
            if (n_procs.eq.1) then
                file_name = 'U_1'
            else
                file_name = 'U_1'//trim(i2str(rank))
            end if
            do ld = 1,size(var_names)
                var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
            end do
            call plot_HDF5(var_names,'RE_'//trim(file_name),&
                &rp(X_1%U_1),col_id=4,col=1)
            call plot_HDF5(var_names,'IM_'//trim(file_name),&
                &ip(X_1%U_1),col_id=4,col=1)
            deallocate(var_names)
            ! DU_0
            allocate(var_names(size(X_1%DU_0,4)))
            if (n_procs.eq.1) then
                file_name = 'DU_0'
            else
                file_name = 'DU_0'//trim(i2str(rank))
            end if
            do ld = 1,size(var_names)
                var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
            end do
            call plot_HDF5(var_names,'RE_'//trim(file_name),&
                &rp(X_1%DU_0),col_id=4,col=1)
            call plot_HDF5(var_names,'IM_'//trim(file_name),&
                &ip(X_1%DU_0),col_id=4,col=1)
            deallocate(var_names)
            ! DU_1
            allocate(var_names(size(X_1%DU_1,4)))
            if (n_procs.eq.1) then
                file_name = 'DU_1'
            else
                file_name = 'DU_1'//trim(i2str(rank))
            end if
            do ld = 1,size(var_names)
                var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
            end do
            call plot_HDF5(var_names,'RE_'//trim(file_name),&
                &rp(X_1%DU_1),col_id=4,col=1)
            call plot_HDF5(var_names,'IM_'//trim(file_name),&
                &ip(X_1%DU_1),col_id=4,col=1)
            deallocate(var_names)
        end if
        
        if (plot_info) then
            ! plot information for comparison between VMEC and HELENA
            call plot_info_for_VMEC_HEL_comparison(grid_X,X_1)
        end if
#endif
#if ldebug
    contains
        !> \private
        subroutine plot_info_for_VMEC_HEL_comparison(grid_X,X)
            use input_utilities, only: pause_prog, get_int, get_log
            use num_vars, only: use_pol_flux_F
            
            ! input / output
            type(grid_type), intent(in) :: grid_X
            type(X_1_type), intent(in) :: X
            
            ! local variables
            real(dp), allocatable :: x_plot(:,:,:), y_plot(:,:,:), z_plot(:,:,:)
            integer :: id, ld
            logical :: not_ready             
            
            call writo('Plotting information for comparison between VMEC and &
                &HELENA',alert=.true.)
            not_ready = .true.
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
                call plot_HDF5('RE_U_0','RE_U_0',rp(X%U_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_U_0','IM_U_0',ip(X%U_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('RE_U_1','RE_U_1',rp(X%U_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_U_1','IM_U_1',ip(X%U_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                
                ! DU_0 and DU_1
                call plot_HDF5('RE_DU_0','RE_DU_0',rp(X%DU_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_DU_0','IM_DU_0',ip(X%DU_0(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('RE_DU_1','RE_DU_1',rp(X%DU_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                call plot_HDF5('IM_DU_1','IM_DU_1',ip(X%DU_1(:,:,:,ld)),&
                    &x=x_plot,y=y_plot,z=z_plot)
                
                ! clean up
                deallocate(x_plot,y_plot,z_plot)
                
                call writo('Want to redo the plotting?')
                not_ready = get_log(.true.)
            end do
            
            call writo('Done, paused',alert=.true.)
            call pause_prog()
        end subroutine plot_info_for_VMEC_HEL_comparison
#endif
    end function run_driver_X_1
    
    !> Part 2 of driver_X: Tensorial jobs.
    !!
    !! \note  Everything  is  done  in the  field-aligned  grids,  where  HELENA
    !! vectorial perturbation variables are first interpolated.
    !!
    !! \return ierr
    integer function run_driver_X_2(grid_eq_B,grid_X,grid_X_B,eq_1,eq_2_B,X_1,&
        &X_2_int) result(ierr)
        use PB3D_ops, only: reconstruct_PB3D_X_2
        use num_vars, only: X_job_nr, X_jobs_lims, rank, n_procs, eq_style, &
            &eq_job_nr, eq_jobs_lims, magn_int_style, no_output, &
            &rich_restart_lvl
        use rich_vars, only: rich_lvl
        use X_ops, only: calc_X, print_output_X, calc_magn_ints, divide_X_jobs
        use HELENA_ops, only: interp_HEL_on_grid
        use X_utilities, only: do_X
#if ldebug
        use num_vars, only: use_pol_flux_F
        use num_utilities, only: c, con
        use X_vars, only: min_m_X, min_n_X, n_mod_X
        use grid_vars, only: alpha, n_alpha
        use grid_utilities, only: trim_grid
        use eq_vars, only: max_flux_F
#endif
        
        character(*), parameter :: rout_name = 'run_driver_X_2'
        
        ! input / output
        type(grid_type), intent(in), pointer :: grid_eq_B                       !< field-aligned equilibrium grid
        type(grid_type), intent(in), target :: grid_X                           !< perturbation grid
        type(grid_type), intent(in), pointer :: grid_X_B                        !< field-aligned perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in), pointer :: eq_2_B                          !< field-aligned metric equilibrium variables
        type(X_1_type), intent(in) :: X_1                                       !< vectorial X variables
        type(X_2_type), intent(inout) :: X_2_int                                !< tensorial X variables
        
        ! local variables
        logical :: no_output_loc                                                ! local copy of no_output
        type(X_1_type) :: X_1_loc(2)                                            ! local X_1
        type(X_2_type) :: X_2                                                   ! tensorial X variables
        integer :: id                                                           ! counter
        integer :: var_size_without_par                                         ! size of array for jobs division
        integer :: prev_style                                                   ! style to treat prev_X
        integer :: lims_loc(2,2)                                                ! local X_jobs_lims
#if ldebug
        character(len=max_str_ln), allocatable :: var_names(:)                  ! names of variables
        character(len=max_str_ln) :: file_name                                  ! name of file
        integer :: k, m                                                         ! local row and column index in flux surface
        integer :: ld, kd, rd, cd                                               ! counters
        integer :: min_nm_X                                                     ! minimal n (tor. flux) or m (pol. flux)
        integer :: n_mod_tot                                                    ! local number of modes
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grid
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_offset(4)                                               ! local offset of plot
        integer :: c_tot(2)                                                     ! total index for symmetric and asymmetric quantities
        complex(dp), allocatable :: PV_int(:,:,:,:,:)                           ! integrated PV_i for all mode combinations
        complex(dp), allocatable :: KV_int(:,:,:,:,:)                           ! integrated KV_i for all mode combinations
        real(dp), allocatable :: X_plot(:,:,:,:)                                ! X of plot
        real(dp), allocatable :: Y_plot(:,:,:,:)                                ! Y of plot
        real(dp), allocatable :: Z_plot(:,:,:,:)                                ! Y of plot
        type(grid_type) :: grid_trim                                            ! trimmed grid
#endif
        
        ! initialize ierr
        ierr = 0
        
        call writo('Setting up Tensorial perturbation jobs')
        call lvl_ud(1)
        
        ! divide perturbation jobs, tensor phase
        var_size_without_par = ceiling(product(grid_X_B%n(1:3))*1._dp)
        ierr = divide_X_jobs(var_size_without_par)
        CHCKERR('')
        
        ! if first equilibrium job of  first Richardson level, set up integrated
        ! tensorial perturbation variables
        if (rich_lvl.eq.rich_restart_lvl .and. eq_job_nr.eq.1) then
            if (rich_lvl.eq.1) then
                call X_2_int%init(mds_X,grid_X_B,is_field_averaged=.true.)
                X_2_int%PV_0 = 0._dp
                X_2_int%PV_1 = 0._dp
                X_2_int%PV_2 = 0._dp
                X_2_int%KV_0 = 0._dp
                X_2_int%KV_1 = 0._dp
                X_2_int%KV_2 = 0._dp
            else
                ierr = reconstruct_PB3D_X_2(mds_X,grid_X,X_2_int,'X_2_int',&
                    &rich_lvl=rich_lvl-1,is_field_averaged=.true.)
                CHCKERR('')
            end if
        end if
        
        ! user output
        call writo('The interpolation method used for the magnetic integrals:')
        call lvl_ud(1)
        select case (magn_int_style)
            case (1)
                call writo('magnetic interpolation style: 1 (trapezoidal rule)')
            case (2)
                call writo('magnetic interpolation style: 2 (Simpson 3/8 rule)')
        end select
        call lvl_ud(-1)
        
        call lvl_ud(-1)
        call writo('Tensorial perturbation jobs set up')
        
        ! main loop over tensorial jobs
        X_job_nr = 0
        call writo('Starting tensorial perturbation jobs')
        call lvl_ud(1)
        X_jobs: do while(do_X())
            ! user output
            call print_info_X_2()
            call lvl_ud(1)
            no_output_loc = no_output
            no_output = .true.
            
            ! set local X_jobs limits
            lims_loc = reshape(X_jobs_lims(:,X_job_nr),[2,2])
            
            ! set local X_1
            do id = 1,2
                ! create X_1_loc
                call X_1_loc(id)%init(mds_X,grid_X,lim_sec_X=lims_loc(:,id))
                
                ! fill X_1_loc
                X_1_loc(id)%U_0 = X_1%U_0(:,:,:,lims_loc(1,id):lims_loc(2,id))
                X_1_loc(id)%U_1 = X_1%U_1(:,:,:,lims_loc(1,id):lims_loc(2,id))
                X_1_loc(id)%DU_0 = X_1%DU_0(:,:,:,lims_loc(1,id):lims_loc(2,id))
                X_1_loc(id)%DU_1 = X_1%DU_1(:,:,:,lims_loc(1,id):lims_loc(2,id))
                
                ! possibly interpolate
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        ! do nothing
                    case (2)                                                    ! HELENA
                        ierr = interp_HEL_on_grid(grid_X,grid_X_B,&
                            &X_1=X_1_loc(id),grid_name=&
                            &'field-aligned perturbation grid')
                        CHCKERR('')
                end select
            end do
            
            ! calculate X variables, tensor phase
            ierr = calc_X(mds_X,grid_eq_B,grid_X_B,eq_1,eq_2_B,&
                &X_1_loc(1),X_1_loc(2),X_2,lim_sec_X=lims_loc)
            CHCKERR('')
            
#if ldebug
            ! write field-aligned tensorial perturbation quantities to output
            if (debug_run_driver_X_2) then
                ! angles
                call plot_HDF5('theta_F_B','theta_F_B',grid_X_B%theta_F)
                call plot_HDF5('zeta_F_B','zeta_F_B',grid_X_B%zeta_F)
                ! PV_0
                allocate(var_names(size(X_2%PV_0,4)))
                if (n_procs.eq.1) then
                    file_name = 'PV_0'
                else
                    file_name = 'PV_0_'//trim(i2str(rank))
                end if
                do ld = 1,size(var_names)
                    var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
                end do
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(X_2%PV_0),col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(X_2%PV_0),col_id=4,col=1)
                deallocate(var_names)
                ! PV_1
                allocate(var_names(size(X_2%PV_1,4)))
                if (n_procs.eq.1) then
                    file_name = 'PV_1'
                else
                    file_name = 'PV_1_'//trim(i2str(rank))
                end if
                do ld = 1,size(var_names)
                    var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
                end do
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(X_2%PV_1),col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(X_2%PV_1),col_id=4,col=1)
                deallocate(var_names)
                ! PV_2
                allocate(var_names(size(X_2%PV_2,4)))
                if (n_procs.eq.1) then
                    file_name = 'PV_2'
                else
                    file_name = 'PV_2_'//trim(i2str(rank))
                end if
                do ld = 1,size(var_names)
                    var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
                end do
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(X_2%PV_2),col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(X_2%PV_2),col_id=4,col=1)
                deallocate(var_names)
                ! KV_0
                allocate(var_names(size(X_2%KV_0,4)))
                if (n_procs.eq.1) then
                    file_name = 'KV_0'
                else
                    file_name = 'KV_0_'//trim(i2str(rank))
                end if
                do ld = 1,size(var_names)
                    var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
                end do
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(X_2%KV_0),col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(X_2%KV_0),col_id=4,col=1)
                deallocate(var_names)
                ! KV_1
                allocate(var_names(size(X_2%KV_1,4)))
                if (n_procs.eq.1) then
                    file_name = 'KV_1'
                else
                    file_name = 'KV_1_'//trim(i2str(rank))
                end if
                do ld = 1,size(var_names)
                    var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
                end do
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(X_2%KV_1),col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(X_2%KV_1),col_id=4,col=1)
                deallocate(var_names)
                ! KV_2
                allocate(var_names(size(X_2%KV_2,4)))
                if (n_procs.eq.1) then
                    file_name = 'KV_2'
                else
                    file_name = 'KV_2_'//trim(i2str(rank))
                end if
                do ld = 1,size(var_names)
                    var_names(ld) = trim(file_name)//'_'//trim(i2str(ld))
                end do
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(X_2%KV_2),col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(X_2%KV_2),col_id=4,col=1)
                deallocate(var_names)
            end if
#endif
            
            ! integrate magnetic  integrals of tensorial  perturbation variables
            ! over  field-aligned  grid.  If   not  first  Richardson  level  or
            ! equilibrium job, previous integrated X_2 is used:
            !   - rich_lvl = 1, eq_job_nr = 1: just calculate integral
            !   - rich_lvl = 1, eq_job_nr > 1: sum to previous integral
            !   - rich_lvl > 1, eq_job_nr = 1: sum to previous integral / 2
            !   - rich_lvl > 1, eq_job_nr > 1: sum to previous integral
            if (rich_lvl.eq.1) then
                if (eq_job_nr.eq.1) then
                    prev_style = 0                                              ! don't add
                else
                    prev_style = 1                                              ! add
                end if
            else
                if (eq_job_nr.eq.1) then
                    prev_style = 2                                              ! divide by 2 and add, change indices
                else
                    prev_style = 3                                              ! add, change indices
                end if
            end if
            ierr = calc_magn_ints(grid_eq_B,grid_X_B,eq_2_B,X_2,X_2_int,&
                &prev_style,lim_sec_X=lims_loc)
            CHCKERR('')
            
            ! clean up
            do id = 1,2
                call X_1_loc(id)%dealloc()
            end do
            call X_2%dealloc()
            
            ! user output
            no_output = no_output_loc
            call lvl_ud(-1)
        end do X_jobs
        
#if ldebug
        ! write  integrated field-aligned  tensorial perturbation  quantities to
        ! output and plot
        if (debug_run_driver_X_2) then
            ! trim grid
            ierr = trim_grid(grid_X,grid_trim,norm_id=norm_id)
            CHCKERR('')
            
            ! set local n_mod and allocate integrated quantities
            n_mod_tot = size(mds_X%sec_ind,2)
            allocate(PV_int(n_mod_tot,n_mod_tot,grid_trim%loc_n_r,&
                &grid_trim%n(2),0:2))
            allocate(KV_int(n_mod_tot,n_mod_tot,grid_trim%loc_n_r,&
                &grid_trim%n(2),0:2))
            allocate(X_plot(n_mod_tot,n_mod_tot,grid_trim%loc_n_r,1))
            allocate(Y_plot(n_mod_tot,n_mod_tot,grid_trim%loc_n_r,1))
            allocate(Z_plot(n_mod_tot,n_mod_tot,grid_trim%loc_n_r,1))
            PV_int = 0._dp
            KV_int = 0._dp
            
            ! minimal index in sec_ind modes
            if (use_pol_flux_F) then
                min_nm_X = minval(min_m_X)
            else
                min_nm_X = minval(min_n_X)
            end if
            
            ! setup X, Y and Z of plot
            ! loop over all normal grid points
            do kd = 1,grid_trim%loc_n_r
                Z_plot(:,:,kd,1) = 100*grid_trim%loc_r_F(kd)/max_flux_F*2*pi
                ! loop over all possible mode combinations
                do cd = 1,n_mod_tot
                    do rd = 1,n_mod_tot
                        Y_plot(rd,cd,kd,1) = min_nm_X+cd-1
                        X_plot(rd,cd,kd,1) = min_nm_X+rd-1
                    end do
                end do
            end do
            
            ! loop over all normal grid points
            do kd = 1,grid_trim%loc_n_r
                ! loop over all mode combinations of this X job
                do m = 1,n_mod_X
                    do k = 1,n_mod_X
                        ! index in tables
                        c_tot(1) = c([k,m],.true.,n_mod_X)
                        c_tot(2) = c([k,m],.false.,n_mod_X)
                        
                        ! total mode combinations indices
                        if (use_pol_flux_F) then
                            rd = mds_X%m(grid_trim%i_min-1+kd,k)
                            cd = mds_X%m(grid_trim%i_min-1+kd,m)
                        else
                            rd = mds_X%n(grid_trim%i_min-1+kd,k)
                            cd = mds_X%n(grid_trim%i_min-1+kd,m)
                        end if
                        cd = cd-min_nm_X+1
                        rd = rd-min_nm_X+1
                        
                        ! symmetric quantities
                        if (k.ge.m) then
                            PV_int(rd,cd,kd,:,0) = con(X_2_int%&
                                &PV_0(1,:,norm_id(1)-1+kd,c_tot(1)),&
                                &[rd,cd],.true.,[n_alpha])
                            PV_int(rd,cd,kd,:,2) = con(X_2_int%&
                                &PV_2(1,:,norm_id(1)-1+kd,c_tot(1)),&
                                &[rd,cd],.true.,[n_alpha])
                            KV_int(rd,cd,kd,:,0) = con(X_2_int%&
                                &KV_0(1,:,norm_id(1)-1+kd,c_tot(1)),&
                                &[rd,cd],.true.,[n_alpha])
                            KV_int(rd,cd,kd,:,2) = con(X_2_int%&
                                &KV_2(1,:,norm_id(1)-1+kd,c_tot(1)),&
                                &[rd,cd],.true.,[n_alpha])
                            PV_int(cd,rd,kd,:,0) = PV_int(rd,cd,kd,:,0)
                            PV_int(cd,rd,kd,:,2) = PV_int(rd,cd,kd,:,2)
                            KV_int(cd,rd,kd,:,0) = KV_int(rd,cd,kd,:,0)
                            KV_int(cd,rd,kd,:,2) = KV_int(rd,cd,kd,:,2)
                        end if
                        
                        ! asymmetric quantities
                        PV_int(rd,cd,kd,:,1) = &
                            &X_2_int%PV_1(1,:,norm_id(1)-1+kd,c_tot(2))
                        KV_int(rd,cd,kd,:,1) = &
                            &X_2_int%KV_1(1,:,norm_id(1)-1+kd,c_tot(2))
                    end do
                end do
            end do
            
            ! plot
            allocate(var_names(n_alpha))
            do ld = 1,size(var_names)
                var_names(ld) = 'alpha = '//trim(r2strt(alpha(ld)))
            end do
            plot_dim = [n_mod_tot,n_mod_tot,grid_trim%n(3),n_alpha]
            plot_offset = [0,0,grid_trim%i_min-1,0]
            do id = 0,2
                file_name = 'PV_'//trim(i2str(id))//'_int_R'//&
                    &trim(i2str(rich_lvl))
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(PV_int(:,:,:,:,id)),x=X_plot,y=Y_plot,z=Z_plot,&
                    &tot_dim=plot_dim, loc_offset=plot_offset,&
                    &col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(PV_int(:,:,:,:,id)),x=X_plot,y=Y_plot,z=Z_plot,&
                    &tot_dim=plot_dim, loc_offset=plot_offset,&
                    &col_id=4,col=1)
                
                file_name = 'KV_'//trim(i2str(id))//'_int_R'//&
                    &trim(i2str(rich_lvl))
                call plot_HDF5(var_names,'RE_'//trim(file_name),&
                    &rp(KV_int(:,:,:,:,id)),x=X_plot,y=Y_plot,z=Z_plot,&
                    &tot_dim=plot_dim, loc_offset=plot_offset,&
                    &col_id=4,col=1)
                call plot_HDF5(var_names,'IM_'//trim(file_name),&
                    &ip(KV_int(:,:,:,:,id)),x=X_plot,y=Y_plot,z=Z_plot,&
                    &tot_dim=plot_dim, loc_offset=plot_offset,&
                    &col_id=4,col=1)
            end do
            deallocate(var_names)
            
            ! clean up
            call grid_trim%dealloc()
        end if
#endif
        
        call lvl_ud(-1)
        call writo('Tensorial perturbation jobs finished')
        
        ! write field-averaged tensorial perturbation variables to output
        if (eq_job_nr.eq.size(eq_jobs_lims,2)) then
            ierr = print_output_X(grid_X_B,X_2_int,'X_2_int',rich_lvl=rich_lvl,&
                &is_field_averaged=.true.)
            CHCKERR('')
        end if
    end function run_driver_X_2
    
    !> Prints information for tensorial perturbation job.
    subroutine print_info_X_2()
        use num_vars, only: X_job_nr, X_jobs_lims
        
        ! user output
        call writo('Job '//trim(i2str(X_job_nr))//'/'//&
            &trim(i2str(size(X_jobs_lims,2)))//': mode numbers ('//&
            &trim(i2str(X_jobs_lims(1,X_job_nr)))//'..'//&
            &trim(i2str(X_jobs_lims(2,X_job_nr)))//')x('//&
            &trim(i2str(X_jobs_lims(3,X_job_nr)))//'..'//&
            &trim(i2str(X_jobs_lims(4,X_job_nr)))//')')
    end subroutine print_info_X_2
end module driver_X
