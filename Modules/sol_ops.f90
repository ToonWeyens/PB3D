!------------------------------------------------------------------------------!
!> Operations on the solution vectors such as decompososing the energy, etc...
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi, rank
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, modes_type
    use sol_vars, only: sol_type
    use vac_vars, only: vac_type

    implicit none
    private
    public plot_harmonics, plot_sol_vec, decompose_energy, print_output_sol, &
        &plot_sol_vals
#if ldebug
    public debug_plot_sol_vec, debug_calc_E, debug_DU, debug_X_norm
#endif
    
    ! global variables
#if ldebug
    logical :: debug_plot_sol_vec = .true.                                     !< plot debug information for plot_sol_vec() \ldebug
    logical :: debug_calc_E = .false.                                           !< plot debug information for calc_E() \ldebug
    logical :: debug_X_norm = .false.                                           !< plot debug information \c X_norm \ldebug
    logical :: debug_DU = .false.                                               !< plot debug information for calculation of \c DU \ldebug
#endif

contains
    !> Plots Eigenvalues.
    subroutine plot_sol_vals(sol,last_unstable_id)
        use num_vars, only: rank
        
        ! input / output
        type(sol_type), intent(in) :: sol                                       !< solution variables
        integer, intent(in) :: last_unstable_id                                 !< index of last unstable EV
        
        ! local variables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        integer :: n_sol_found                                                  ! how many solutions found and saved
        integer :: id                                                           ! counter
        
        ! set local variables
        n_sol_found = size(sol%val)
        
        ! only let master plot
        if (rank.eq.0) then
            ! Last Eigenvalues
            plot_title = 'final Eigenvalues omega^2 [log]'
            plot_name = 'Eigenvalues'
            call print_ex_2D(plot_title,plot_name,&
                &log10(abs(rp(sol%val(1:n_sol_found)))),draw=.false.)
            call draw_ex([plot_title],plot_name,1,1,.false.,ex_plot_style=1)
            
            ! Last Eigenvalues: unstable range
            if (last_unstable_id.gt.0) then
                plot_title = 'final unstable Eigenvalues omega^2'
                plot_name = 'Eigenvalues_unstable'
                call print_ex_2D(plot_title,plot_name,&
                    &rp(sol%val(1:last_unstable_id)),&
                    &x=[(id*1._dp,id=1,last_unstable_id)],draw=.false.)
                call draw_ex([plot_title],plot_name,1,1,.false.,ex_plot_style=1)
            end if
            
            ! Last Eigenvalues: stable range
            if (last_unstable_id.lt.n_sol_found) then
                plot_title = 'final stable Eigenvalues omega^2'
                plot_name = 'Eigenvalues_stable'
                call print_ex_2D(plot_title,plot_name,&
                    &rp(sol%val(last_unstable_id+1:n_sol_found)),&
                    &x=[(id*1._dp,id=last_unstable_id+1,n_sol_found)],&
                    &draw=.false.)
                call draw_ex([plot_title],plot_name,1,1,.false.,ex_plot_style=1)
            end if
        end if
    end subroutine plot_sol_vals
    
    !> Plots Eigenvectors.
    !!
    !! This is done using the angular  part of the the provided equilibrium grid
    !! and the normal part of the provided solution grid.
    !!
    !! The perturbtion grid  is assumed to have the same  angular coordinates as
    !! the equilibrium grid, and the normal coordinates correspond to either the
    !! equilibrium grid (X_grid_style 1) or  the solution grid (X_grid_style 2).
    !! The output grid,  furthermore, has the angular part  corresponding to the
    !! equilibrium grid, and the normal part to the solution grid.
    !!
    !! \note The normalization factors are taken  into account and the output is
    !! transformed back to unnormalized values:
    !!  \f[\begin{aligned}
    !!      \vec{\xi} &\sim \frac{X_0}{R_0 B_0} \ , \\
    !!      \vec{Q}   &\sim \frac{X_0}{R_0^2} \ , \\
    !!  \end{aligned}\f]
    !!
    !! which translates to
    !!  \f[\begin{aligned}
    !!      X  &\sim X_0 \ , \\
    !!      U  &\sim \frac{X_0}{R_0^2 B_0} \ , \\
    !!      Qn &\sim \frac{X_0 B_0}{R_0} \ , \\
    !!      Qg &\sim \frac{X_0}{R_0^3} \ ,
    !!  \end{aligned}\f]
    !! where \f$X_0\f$ is not determined but is common to all factors.
    !!
    !! \return ierr
    integer function plot_sol_vec(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,X,&
        &sol,XYZ,X_id,plot_var) result(ierr)
        use num_vars, only: no_plots, tol_zero, pert_mult_factor_POST, &
            &eq_job_nr, eq_jobs_lims, eq_job_nr, norm_disc_prec_X, &
            &norm_disc_prec_sol, use_normalization, X_grid_style
        use grid_utilities, only: trim_grid, calc_vec_comp, setup_interp_data, &
            &apply_disc
        use sol_utilities, only: calc_XUQ
        use eq_vars, only: R_0, B_0
        use num_utilities, only: c
        use MPI_utilities, only: get_ser_var
#if ldebug
        use num_vars, only: use_pol_flux_F
        use grid_utilities, only: setup_deriv_data
#endif
        
        character(*), parameter :: rout_name = 'plot_sol_vec'
        
        ! input / output
        type(modes_type), intent(in) :: mds_sol                                 !< general modes variables in solution grid
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(in) :: X                                         !< perturbation variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        real(dp), intent(in) :: XYZ(:,:,:,:)                                    !< X, Y and Z of extended eq_grid
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue (for output name)
        logical, intent(in) :: plot_var(2)                                      !< whether variables are plotted (1: plasma perturbation, 2: magnetic perturbation)
        
        ! local variables
        type(grid_type) :: grid_out                                             ! output grid (see description)
        type(grid_type) :: grid_out_trim                                        ! trimmed output grid
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: id, jd, kd, td                                               ! counters
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_offset(4)                                               ! local offset of plot
        integer :: col                                                          ! collection type for HDF5 plot
        logical :: cont_plot                                                    ! continued plot
        logical :: plot_var_loc(2)                                              ! local plot_var
        logical :: XYZ_plot_setup                                               ! XYZ_plot_setup has been set up
        real(dp) :: X_0                                                         ! X at theta = zeta = 0
        real(dp), allocatable :: X_0_ser(:)                                     ! X_0 for all processes
        real(dp), allocatable :: time(:)                                        ! fraction of Alfvén time
        real(dp), allocatable :: XYZ_plot(:,:,:,:,:)                            ! copies of XYZ
        real(dp), allocatable :: XYZ_vec(:,:,:,:,:)                             ! copies of XYZ for vector plot
        real(dp), allocatable :: h12(:,:,:)                                     ! interpolated h_FD(1,2) on solution grid
        real(dp), allocatable :: h22(:,:,:)                                     ! interpolated h_FD(2,2) on solution grid
        real(dp), allocatable :: h23(:,:,:)                                     ! interpolated h_FD(3,2) on solution grid
        real(dp), allocatable :: g13(:,:,:)                                     ! interpolated g_FD(1,3) on solution grid
        real(dp), allocatable :: g33(:,:,:)                                     ! interpolated g_FD(3,3) on solution grid
        real(dp), allocatable :: jac_FD_int(:,:,:)                              ! interpolated jac_FD on solution grid
        real(dp), allocatable :: f_plot_phase(:,:,:,:)                          ! phase of f_plot
        real(dp), allocatable :: ccomp(:,:,:,:,:)                               ! Cart. components of perturbation
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        complex(dp), allocatable :: f_plot(:,:,:,:,:)                           ! the function to plot
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: var_name(2)                                ! name of variable that is plot
        character(len=max_str_ln) :: file_name(2)                               ! name of file
        character(len=max_str_ln) :: description(2)                             ! description
        character(len=2) :: sol_name                                            ! name of solution vector ('xi' or 'Q')
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
#if ldebug
        type(disc_type) :: norm_deriv_data                                      ! normal derivative data
        integer :: nm                                                           ! n (pol. flux) or m (tor. flux)
        integer :: ld                                                           ! counter
        real(dp), allocatable :: Dq_saf(:)                                      ! norm. deriv. of q_saf interpolated on solution grid
        real(dp), allocatable :: Drot_t(:)                                      ! norm. deriv. of rot_t interpolated on solution grid
        complex(dp), allocatable :: U_inf(:,:,:,:)                              ! ideal ballooning U
        complex(dp), allocatable :: U_inf_prop(:,:,:)                           ! proportional part of U_inf
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! return if no variables are plot
        if (count(plot_var).eq.0) return
        plot_var_loc = plot_var
        if (abs(pert_mult_factor_POST).ge.tol_zero) then
            ! set up X_0
            X_0 = 0._dp
            do kd = 1,grid_sol%loc_n_r
                X_0 = max(X_0,abs(sum(sol%vec(:,kd,X_id))))
            end do
            ierr = get_ser_var([X_0],X_0_ser,scatter=.true.)
            CHCKERR('')
            X_0 = maxval(X_0_ser)
            
            ! user output
            call writo('Perturbing the position vector (X,Y,Z) with &
                &multiplicative factor '//trim(r2strt(pert_mult_factor_POST)))
            if (.not.plot_var(1)) then
                call writo('For nonzero "pert_mult_factor_POST", need to also &
                    &plot xi',warning=.true.)
                plot_var_loc(1) = .true.
            end if
        end if
        
        ! user output
        call writo('Plot the solution vector')
        call lvl_ud(1)
        
        ! tests
        if (size(XYZ,4).ne.3) then
            ierr = 1
            err_msg = 'X, Y and Z needed'
            CHCKERR(err_msg)
        end if
        if (grid_eq%n(1).ne.size(XYZ,1) .or. &
            &grid_eq%n(2).ne.size(XYZ,2) .or. &
            &grid_sol%loc_n_r.ne.size(XYZ,3)) then
            ierr = 1
            err_msg = 'XYZ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! set up n_t
        ! if the  Eigenvalue is negative,  the Eigenfunction explodes,  so limit
        ! n_t(2) to 1. If it is positive, the Eigenfunction oscilates, so choose
        ! n_t(2) = 8 for 2 whole periods
        if (rp(sol%val(X_id)).lt.0) then                                        ! exploding, unstable
            n_t(1) = 1                                                          ! 1 point per quarter period
            n_t(2) = 1                                                          ! 1 quarter period
        else                                                                    ! oscillating, stable
            n_t(1) = 2                                                          ! 2 points per quarter period
            n_t(2) = 4                                                          ! 4 quarter periods
        end if
        
        ! set collection type
        if (product(n_t).gt.1) then
            col = 2                                                             ! temporal collection
        else
            col = 1                                                             ! no collection
        end if
        
        ! setup normal interpolation data X->sol
        ierr = setup_interp_data(grid_X%loc_r_F,grid_sol%loc_r_F,&
            &norm_interp_data,norm_disc_prec_X)
        CHCKERR('')
        
        ! set up output grid
        ierr = grid_out%init([grid_X%n(1:2),grid_sol%n(3)],&
            &i_lim=[grid_sol%i_min,grid_sol%i_max],divided=grid_sol%divided)
        CHCKERR('')
        grid_out%r_F = grid_sol%r_F
        grid_out%r_E = grid_sol%r_E
        grid_out%loc_r_F = grid_sol%loc_r_F
        grid_out%loc_r_E = grid_sol%loc_r_E
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! interpolate
                ierr = apply_disc(grid_X%theta_F,norm_interp_data,&
                    grid_out%theta_F,3)
                CHCKERR('')
                ierr = apply_disc(grid_X%theta_E,norm_interp_data,&
                    grid_out%theta_E,3)
                CHCKERR('')
                ierr = apply_disc(grid_X%zeta_F,norm_interp_data,&
                    grid_out%zeta_F,3)
                CHCKERR('')
                ierr = apply_disc(grid_X%zeta_E,norm_interp_data,&
                    grid_out%zeta_E,3)
                CHCKERR('')
            case (2)                                                            ! solution
                ! copy
                grid_out%theta_F = grid_X%theta_F
                grid_out%theta_E = grid_X%theta_E
                grid_out%zeta_F = grid_X%zeta_F
                grid_out%zeta_E = grid_X%zeta_E
        end select
        
        ! trim output grid
        ierr = trim_grid(grid_out,grid_out_trim,norm_id)
        CHCKERR('')
        
        ! set up plot dimensions and local dimensions
        plot_dim = [grid_out_trim%n,product(n_t)]
        plot_offset = [0,0,grid_out_trim%i_min-1,0]
        
        ! possibly modify if multiple equilibrium parallel jobs
        if (size(eq_jobs_lims,2).gt.1) then
            plot_dim(1) = eq_jobs_lims(2,size(eq_jobs_lims,2)) - &
                &eq_jobs_lims(1,1) + 1
            plot_offset(1) = eq_jobs_lims(1,eq_job_nr) - 1
        end if
        
        ! set up continued plot
        cont_plot = eq_job_nr.gt.1
        
        ! For each time step, calculate the time (as fraction of Alfvén time)
        allocate(time(product(n_t)))
        do td = 1,product(n_t)
            if (n_t(1).eq.1) then
                time(td) = (td-1) * 0.25
            else
                time(td) = (mod(td-1,n_t(1))*1._dp/n_t(1) + &
                    &(td-1)/n_t(1)) * 0.25
            end if
        end do
        
        ! set up interpolated metric factors (not trimmed)
        allocate(h12(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
        allocate(h22(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
        allocate(h23(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
        allocate(g13(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
        allocate(g33(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
        allocate(jac_FD_int(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
#if ldebug
        allocate(Dq_saf(grid_out%loc_n_r))
        allocate(Drot_t(grid_out%loc_n_r))
#endif
        
        ! interpolate eq->sol
        ierr = apply_disc(eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,0),&
            &norm_interp_data,h12,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0),&
            &norm_interp_data,h22,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%h_FD(:,:,:,c([2,3],.true.),0,0,0),&
            &norm_interp_data,h23,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%g_FD(:,:,:,c([1,3],.true.),0,0,0),&
            &norm_interp_data,g13,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0),&
            &norm_interp_data,g33,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%jac_FD(:,:,:,0,0,0),&
            &norm_interp_data,jac_FD_int,3)
        CHCKERR('')
#if ldebug
        ierr = apply_disc(eq_1%q_saf_FD(:,1),norm_interp_data,Dq_saf)
        CHCKERR('')
        ierr = apply_disc(eq_1%rot_t_FD(:,1),norm_interp_data,Drot_t)
        CHCKERR('')
#endif
        
        ! calculate omega =  sqrt(sol_val) and make sure to  select the decaying
        ! solution
        omega = sqrt(sol%val(X_id))
        if (ip(omega).gt.0) omega = - omega                                     ! exploding solution, not the decaying one
        
        ! set up f_plot (not trimmed) and XYZ_plot (trimmed)
        ! set up:
        !   * f_plot (not trimmed): X, U, Q
        !   * ccomp (not trimmed): Cartesian components of perturbation
        !   * XYZ_plot (trimmed): copies of XYZ for plot for different times
        !   * XYZ_vec (trimmed): copies of XYZ for vector plot
        allocate(f_plot(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r,&
            &product(n_t),2))
        allocate(ccomp(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r,3,2))
        allocate(XYZ_plot(grid_out%n(1),grid_out%n(2),grid_out_trim%loc_n_r,&
            &size(time),3))
        allocate(XYZ_vec(grid_out%n(1),grid_out%n(2),grid_out_trim%loc_n_r,3,&
            &3))
        XYZ_plot_setup = .false.
        
        ! operations for plasma perturbation and magnetic perturbation
        do jd = 1,2                                                             ! first plasma perturbation, then magnetic perturbation
            if (.not.plot_var_loc(jd)) cycle
            
            ! calculate vector components and  set solution name, variable name,
            ! file name and description
            select case (jd)
                case (1)                                                        ! plasma perturbation
                    ierr = calc_XUQ(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,&
                        &X,sol,X_id,1,time,f_plot(:,:,:,:,1))
                    CHCKERR('')
                    ierr = calc_XUQ(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,&
                        &X,sol,X_id,2,time,f_plot(:,:,:,:,2))
                    CHCKERR('')
                    sol_name = 'xi'
                    file_name(1) = trim(i2str(X_id))//'_sol_X'
                    file_name(2) = trim(i2str(X_id))//'_sol_U'
                case (2)                                                        ! magnetic perturbation
                    ierr = calc_XUQ(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,&
                        &X,sol,X_id,3,time,f_plot(:,:,:,:,1))
                    CHCKERR('')
                    ierr = calc_XUQ(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,&
                        &X,sol,X_id,4,time,f_plot(:,:,:,:,2))
                    CHCKERR('')
                    sol_name = 'Q'
                    file_name(1) = trim(i2str(X_id))//'_sol_Qn'
                    file_name(2) = trim(i2str(X_id))//'_sol_Qg'
            end select
            var_name(1) = 'Normal comp. of '//trim(sol_name)
            var_name(2) = 'Geodesic comp. of '//trim(sol_name)
            description(1) = 'Normal component of vector '//trim(sol_name)//&
                &' for Eigenvalue '//trim(i2str(X_id))//' with omega = '//&
                &trim(r2str(rp(omega)))
            description(2) = 'Geodesic component of vector '//trim(sol_name)//&
                &' for Eigenvalue '//trim(i2str(X_id))//' with omega = '//&
                &trim(r2str(rp(omega)))
            
            ! calculate vector components at every time point
            do td = 1,product(n_t)
                ! covariant components:
                !   xi . e_alpha = U J^2 h22/g33
                !   xi . e_psi   = X 1/h22 - U J^2 h12/g33
                !   xi . e_theta = 0
                ! and similar for Q
                ccomp(:,:,:,1,1) = rp(f_plot(:,:,:,td,2)) * &
                    &jac_FD_int**2*h22/g33
                ccomp(:,:,:,2,1) = rp(f_plot(:,:,:,td,1)) * &
                    &1._dp/h22 - &
                    &rp(f_plot(:,:,:,td,2)) * &
                    &jac_FD_int**2*h12/g33
                ccomp(:,:,:,3,1) = 0._dp
                
                ! contravariant components:
                !   xi . nabla alpha = X h12/h22 + U
                !   xi . nabla psi   = X
                !   xi . nabla theta = X h23/h22 - U g13/g33
                ! and similar for Q
                ccomp(:,:,:,1,2) = rp(f_plot(:,:,:,td,1)) * &
                    &h12/h22 + rp(f_plot(:,:,:,td,2))
                ccomp(:,:,:,2,2) = rp(f_plot(:,:,:,td,1))
                ccomp(:,:,:,3,2) = rp(f_plot(:,:,:,td,1)) * &
                    &h23/h22 - rp(f_plot(:,:,:,td,2)) * g13/g33 
                
                ! multiplication factor if xi
                if (abs(pert_mult_factor_POST).ge.tol_zero .and. jd.eq.1) &
                    &ccomp = ccomp*pert_mult_factor_POST/X_0
                
                ! transform normalized quantities to unnormalized
                if (use_normalization) then
                    select case (jd)
                        case (1)                                                ! plasma perturbation
                            ccomp = ccomp / (R_0*B_0)                           ! xi ~ X_0/(R_0*B_0)
                        case (2)                                                ! magnetic perturbation
                            ccomp = ccomp / (R_0**2)                            ! Q ~ X_0/(R_0^2)
                    end select
                end if
                
                ! transform to Cartesian components
                ierr = calc_vec_comp(grid_out,grid_eq,eq_1,eq_2,ccomp,&
                    &norm_disc_prec_sol)
                CHCKERR('')
                
                ! vector plot
                do id = 1,3
                    XYZ_vec(:,:,:,id,:) = XYZ(:,:,norm_id(1):norm_id(2),:)
                end do
                call plot_HDF5([trim(sol_name)],trim(sol_name)//'_T_'//&
                    &trim(i2str(td)),ccomp(:,:,norm_id(1):norm_id(2),:,1),&
                    &tot_dim=[plot_dim(1:3),3],&
                    &loc_offset=[plot_offset(1:3),0],&
                    &X=XYZ_vec(:,:,:,:,1),Y=XYZ_vec(:,:,:,:,2),&
                    &Z=XYZ_vec(:,:,:,:,3),col=4,cont_plot=cont_plot,&
                    &descr='perturbation for time t = '//&
                    &trim(r2strt(time(td))))
                
                ! perturb position vector if xi;  it does not matter whether
                ! we take covariant or contravariant components
                if (.not.XYZ_plot_setup) then
                    XYZ_plot(:,:,:,td,:) = XYZ(:,:,norm_id(1):norm_id(2),:)
                    if (abs(pert_mult_factor_POST).ge.tol_zero .and. jd.eq.1) &
                        &XYZ_plot(:,:,:,td,:) = XYZ_plot(:,:,:,td,:) + &
                        &ccomp(:,:,norm_id(1):norm_id(2),:,1)
                end if
            end do
            
            ! XYZ_plot has been set up if we get here
            XYZ_plot_setup = .true.
        
#if ldebug
            if (debug_plot_sol_vec .and. jd.eq.1) then
                call writo('Checking how well U approximates the pure &
                    &ballooning result')
                call lvl_ud(1)
                
                ! set up ideal U
                allocate(U_inf(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r,&
                    &product(n_t)))
                
                ! derive the X vector
                ierr = setup_deriv_data(grid_out%loc_r_F,norm_deriv_data,1,&
                    &norm_disc_prec_sol)
                CHCKERR('')
                ierr = apply_disc(f_plot(:,:,:,:,1),norm_deriv_data,U_inf,3)
                CHCKERR('')
                call norm_deriv_data%dealloc()
                
                ! set up dummy variable Theta^alpha + q' theta and nm
                allocate(U_inf_prop(grid_out%n(1),grid_out%n(2),&
                    grid_out%loc_n_r))
                if (use_pol_flux_F) then
                    do kd = 1,grid_out%loc_n_r
                        U_inf_prop(:,:,kd) = h12(:,:,kd)/h22(:,:,kd) + &
                            &Dq_saf(kd)*grid_out%theta_F(:,:,kd)
                    end do
                    nm = sol%n(1,1)
                else
                    do kd = 1,grid_out%loc_n_r
                        U_inf_prop(:,:,kd) = h12(:,:,kd)/h22(:,:,kd) - &
                            &Drot_t(kd)*grid_out%zeta_F(:,:,kd)
                    end do
                    nm = sol%m(1,1)
                end if
                
                call print_ex_2D('X_inf','X_inf',rp(f_plot(50,1,:,1,1)),persistent=.true.)
                call print_ex_2D('U_inf','U_inf',rp(U_inf(50,1,:,1)))
                call plot_HDF5(['U_inf'],'U_inf',&
                    &reshape([rp(iu/nm*U_inf(:,:,:,1)),ip(iu/nm*U_inf(:,:,:,1))],&
                    &[grid_out%n(1:2),grid_out%loc_n_r,2]))
                ! multiply by i/n or i/m and add the proportional part
                do ld = 1,product(n_t)
                    U_inf(:,:,:,ld) = iu/nm * U_inf(:,:,:,ld) - &
                        &f_plot(:,:,:,ld,1) * U_inf_prop
                end do
                deallocate(U_inf_prop)
                
                ! plotting real part
                call plot_HDF5('RE X','TEST_RE_X_'//trim(i2str(X_id)),&
                    &rp(f_plot(:,:,norm_id(1):norm_id(2),1,1)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3))
                call plot_diff_HDF5(rp(f_plot(:,:,norm_id(1):norm_id(2),1,2)),&
                    &rp(U_inf(:,:,norm_id(1):norm_id(2),1)),&
                    &'TEST_RE_U_inf_'//trim(i2str(X_id)),tot_dim=plot_dim(1:3),&
                    &loc_offset=plot_offset(1:3),descr='To test whether U &
                    &approximates the ideal ballooning result',&
                    &output_message=.true.)
                
                ! plotting imaginary part
                call plot_HDF5('IM X','TEST_IM_X_'//trim(i2str(X_id)),&
                    &ip(f_plot(:,:,norm_id(1):norm_id(2),1,1)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3))
                call plot_diff_HDF5(ip(f_plot(:,:,norm_id(1):norm_id(2),1,2)),&
                    &ip(U_inf(:,:,norm_id(1):norm_id(2),1)),&
                    &'TEST_IM_U_inf_'//trim(i2str(X_id)),tot_dim=plot_dim(1:3),&
                    &loc_offset=plot_offset(1:3),descr='To test whether U &
                    &approximates the ideal ballooning result',&
                    &output_message=.true.)
                
                call writo('Checking done')
                call lvl_ud(-1)
            end if
#endif
            
            ! set up temporary variable for phase
            allocate(f_plot_phase(grid_out%n(1),grid_out%n(2),&
                &grid_out_trim%loc_n_r,product(n_t)))
            
            ! transform normalized quantities to unnormalized
            if (use_normalization) then
                select case (jd)
                    case (1)                                                    ! plasma perturbation
                        !f_plot(:,:,:,:,1) = f_plot(:,:,:,:,1)                   ! X
                        f_plot(:,:,:,:,2) = f_plot(:,:,:,:,2) / (R_0**2*B_0)    ! U
                    case (2)                                                    ! magnetic perturbation
                        f_plot(:,:,:,:,1) = f_plot(:,:,:,:,1) * B_0/R_0         ! Qn
                        f_plot(:,:,:,:,2) = f_plot(:,:,:,:,2) / (R_0**3)        ! Qg
                end select
            end if
            
            call plot_HDF5([var_name(1)],'TEST',&
                &rp(f_plot(:,:,norm_id(1):norm_id(2),:,1)),&
                &tot_dim=plot_dim,loc_offset=plot_offset,&
                &col=col,cont_plot=cont_plot,&
                &descr=description(1))
            do kd = 1,2
                call plot_HDF5([var_name(kd)],trim(file_name(kd))//'_RE',&
                    &rp(f_plot(:,:,norm_id(1):norm_id(2),:,kd)),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=XYZ_plot(:,:,:,:,1),Y=XYZ_plot(:,:,:,:,2),&
                    &Z=XYZ_plot(:,:,:,:,3),col=col,cont_plot=cont_plot,&
                    &descr=description(kd))
                call plot_HDF5([var_name(kd)],trim(file_name(kd))//'_IM',&
                    &ip(f_plot(:,:,norm_id(1):norm_id(2),:,kd)),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=XYZ_plot(:,:,:,:,1),Y=XYZ_plot(:,:,:,:,2),&
                    &Z=XYZ_plot(:,:,:,:,3),col=col,cont_plot=cont_plot,&
                    &descr=description(kd))
                f_plot_phase = atan2(&
                    &ip(f_plot(:,:,norm_id(1):norm_id(2),:,kd)),&
                    &rp(f_plot(:,:,norm_id(1):norm_id(2),:,kd)))
                where (f_plot_phase.lt.0) f_plot_phase = f_plot_phase + 2*pi
                call plot_HDF5([var_name(kd)],trim(file_name(kd))//'_PH',&
                    &f_plot_phase,tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=XYZ_plot(:,:,:,:,1),Y=XYZ_plot(:,:,:,:,2),&
                    &Z=XYZ_plot(:,:,:,:,3),col=col,cont_plot=cont_plot,&
                    &descr=description(kd))
            end do
            
            ! clean up
            deallocate(f_plot_phase)
            write(*,*) 'WAITING'
            read(*,*)
        end do
        
        ! clean up
        call norm_interp_data%dealloc()
        call grid_out%dealloc()
        call grid_out_trim%dealloc()
        
        call lvl_ud(-1)
    end function plot_sol_vec
    
    !> Plots the harmonics and their maximum in 2-D.
    !!
    !! \return ierr
    integer function plot_harmonics(mds,grid_sol,sol,X_id,res_surf) result(ierr)
        use MPI_utilities, only: get_ser_var
        use num_vars, only: ex_max_size, rank, no_plots, use_pol_flux_F
        use eq_vars, only: max_flux_F
        use X_vars, only: min_n_X, min_m_X
        use sol_utilities, only: calc_tot_sol_vec
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'plot_harmonics'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        type(sol_type), intent(in) :: sol                                       !< solution variables
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue (for output name)
        real(dp), intent(in) :: res_surf(:,:)                                   !< resonant surfaces
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed sol grid
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: ld                                                           ! counter
        integer :: ld_loc                                                       ! local ld
        integer :: n_mod_tot                                                    ! total number of modes that can resonate
        integer :: n_mod_loc                                                    ! local number of modes that are used in the calculations
        integer :: min_nm_X                                                     ! minimum n or m in total modes
        character(len=max_str_ln) :: file_name                                  ! name of file of plots of this proc.
        character(len=max_str_ln), allocatable :: plot_title(:)                 ! title for plots
        real(dp) :: norm_factor                                                 ! conversion factor max_flux/2pi from flux to normal coordinate
        real(dp), allocatable :: x_plot(:,:)                                    ! x values of plot
        real(dp), allocatable :: y_plot(:,:)                                    ! y values of plot
        complex(dp), allocatable :: sol_vec_ser(:,:)                            ! serial MPI Eigenvector for local number of modes
        complex(dp), allocatable :: sol_vec_ser_tot(:,:)                        ! serial MPI Eigenvector for total number of modes
        complex(dp), allocatable :: sol_vec_ser_loc(:)                          ! local sol_vec_ser
        real(dp), allocatable :: sol_vec_phase(:,:)                             ! phase of sol_vec
        character :: nm_X                                                       ! n or m
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! trim sol grid
        ierr = trim_grid(grid_sol,grid_sol_trim,norm_id)
        CHCKERR('')
        
        ! set up local and total nr.  of modes, which can be different for X
        ! style 2 (see discussion in sol_utilities)
        n_mod_loc = sol%n_mod
        n_mod_tot = size(mds%sec_ind,2)
        
        ! if master, set up plot titles
        if (rank.eq.0) then
        end if
        
        ! set up serial sol_vec on master
        if (rank.eq.0) &
            &allocate(sol_vec_ser(1:n_mod_loc,1:grid_sol_trim%n(3)))
        do ld = 1,n_mod_loc
            ierr = get_ser_var(sol%vec(ld,norm_id(1):norm_id(2),X_id),&
                &sol_vec_ser_loc,scatter=.true.)
            CHCKERR('')
            if (rank.eq.0) sol_vec_ser(ld,:) = sol_vec_ser_loc
            deallocate(sol_vec_ser_loc)
        end do
        
        ! convert to full mode on master and set up normal factor
        if (rank.eq.0) then
            allocate(sol_vec_ser_tot(1:n_mod_tot,1:grid_sol_trim%n(3)))
            ierr = calc_tot_sol_vec(mds,grid_sol%i_min,sol_vec_ser,&
                &sol_vec_ser_tot)                                               ! Note: need to use untrimmed grid for minimal i
            CHCKERR('')
            
            norm_factor = max_flux_F/(2*pi)
        end if
        
        ! master sets up x_plot, y_plot and plot titles
        if (rank.eq.0) then
            ! prepare
            allocate(x_plot(grid_sol_trim%n(3),n_mod_tot))
            allocate(y_plot(grid_sol_trim%n(3),n_mod_tot))
            allocate(plot_title(n_mod_tot))
            if (use_pol_flux_F) then
                nm_X = 'm'
                min_nm_X = minval(min_m_X)
            else
                nm_X = 'n'
                min_nm_X = minval(min_n_X)
            end if
            
            ! identical copies of r_F of solution grid and mode numbers
            do ld = 1,n_mod_tot
                x_plot(:,ld) = grid_sol_trim%r_F
                y_plot(:,ld) = min_nm_X+ld-1
                plot_title(ld) = 'EV - midplane ('//nm_X//' = '//&
                    &trim(i2str(min_nm_X+ld-1))//')'
            end do
            
            ! scale x_plot by normal factor
            x_plot = x_plot/norm_factor
        end if
        
        ! master plots real part at midplane
        if (rank.eq.0) then
            ! set up file name of this rank and plot title
            file_name = trim(i2str(X_id))//'_EV_midplane_RE'
            
            ! print real amplitude of harmonics of eigenvector at midplane
            call print_ex_2D(plot_title,file_name,&
                &rp(transpose(sol_vec_ser_tot)),x=x_plot,draw=.false.)
            call draw_ex(plot_title,file_name,n_mod_tot,1,&
                &.false.,draw_ops=['with lines'],ex_plot_style=1)
            call draw_ex(plot_title,file_name,n_mod_tot,1,&
                &.false.,ex_plot_style=2)
            
            ! plot in file using decoupled 3D if not too big
            if (n_mod_tot*grid_sol_trim%n(3).le.ex_max_size) then
                call draw_ex(plot_title,trim(file_name)//'_3D',&
                    &n_mod_tot,3,.false.,draw_ops=['with lines'],&
                    &data_name=file_name,ex_plot_style=1)
            end if
        end if
        
        ! master plots imaginary part at midplane
        if (rank.eq.0) then
            ! set up file name of this rank and plot title
            file_name = trim(i2str(X_id))//'_EV_midplane_IM'
            
            ! print imag amplitude of harmonics of eigenvector at midplane
            call print_ex_2D(plot_title,file_name,&
                &ip(transpose(sol_vec_ser_tot)),x=x_plot,draw=.false.)
            
            ! plot in file
            call draw_ex(plot_title,file_name,n_mod_tot,1,.false.,&
                &draw_ops=['with lines'],ex_plot_style=1)
            
            ! plot in file using decoupled 3D if not too big
            if (n_mod_tot*grid_sol_trim%n(3).le.ex_max_size) then
                call draw_ex(plot_title,trim(file_name)//'_3D',&
                    &n_mod_tot,3,.false.,draw_ops=['with lines'],&
                    &data_name=file_name,ex_plot_style=1)
            end if
        end if
        
        ! master plots phase at midplane
        if (rank.eq.0) then
            ! set up file name of this rank and plot title
            file_name = trim(i2str(X_id))//'_EV_midplane_PH'
            
            ! set up phase
            allocate(sol_vec_phase(grid_sol_trim%n(3),n_mod_tot))
            sol_vec_phase = atan2(ip(transpose(sol_vec_ser_tot)),&
                &rp(transpose(sol_vec_ser_tot)))
            where (sol_vec_phase.lt.0) sol_vec_phase = sol_vec_phase + 2*pi
            
            ! print imag amplitude of harmonics of eigenvector at midplane
            call print_ex_2D(plot_title,file_name,sol_vec_phase,&
                &x=x_plot,draw=.false.)
            
            ! plot in file
            call draw_ex(plot_title,file_name,n_mod_tot,1,&
                &.false.,draw_ops=['with lines'],ex_plot_style=1)
            
            ! plot in file using decoupled 3D if not too big
            if (n_mod_tot*grid_sol_trim%n(3).le.ex_max_size) then
                call draw_ex(plot_title,trim(file_name)//'_3D',&
                    &n_mod_tot,3,.false.,draw_ops=['with lines'],&
                    &data_name=file_name,ex_plot_style=1)
            end if
            
            ! deallocate variables
            deallocate(plot_title)
        end if
        
        ! masters plots above in HDf5
        if (rank.eq.0) then
            allocate(plot_title(3))
            plot_title(1) = 'real part'
            plot_title(2) = 'imaginary part'
            plot_title(3) = 'phase'
            file_name = trim(i2str(X_id))//'_EV_midplane'
            call plot_HDF5(plot_title,trim(file_name),&
                &reshape([rp(transpose(sol_vec_ser_tot)),&
                &ip(transpose(sol_vec_ser_tot)),&
                &sol_vec_phase],[grid_sol_trim%n(3),n_mod_tot,1,3]),&
                &x=reshape(x_plot,[grid_sol_trim%n(3),n_mod_tot,1,1]),&
                &y=reshape(y_plot,[grid_sol_trim%n(3),n_mod_tot,1,1]))
            
            ! deallocate variables
            deallocate(x_plot)
            deallocate(y_plot)
            deallocate(sol_vec_phase)
            deallocate(plot_title)
        end if
        
        ! only plot maximum if resonant surfaces found
        if (rank.eq.0 .and. size(res_surf,1).gt.0) then
            ! set up file name of this rank and plot title
            file_name = trim(i2str(X_id))//'_EV_max'
            allocate(plot_title(2))
            plot_title(1) = 'mode maximum'
            plot_title(2) = 'resonating surface'
            
            ! set up plot variables
            allocate(x_plot(n_mod_tot,2))
            allocate(y_plot(n_mod_tot,2))
            
            ! set up maximum of each mode at midplane
            ld_loc = 0
            do ld = 1,n_mod_tot
                if (maxval(abs(rp(sol_vec_ser_tot(ld,:)))).gt.0) then           ! mode resonates somewhere
                    ld_loc = ld_loc + 1
                    x_plot(ld_loc,1) = grid_sol_trim%r_F(&
                        &maxloc(abs(rp(sol_vec_ser_tot(ld,:))),1))
                    y_plot(ld_loc,1) = min_nm_X+ld-1
                end if
            end do
            if (ld_loc.eq.0) then
                ierr = 1
                err_msg = 'None of the modes resonates'
                CHCKERR(err_msg)
            end if
            x_plot(ld_loc+1:n_mod_tot,1) = x_plot(ld_loc,1)                     ! identical copies of last point for numerical reasons
            y_plot(ld_loc+1:n_mod_tot,1) = y_plot(ld_loc,1)                     ! identical copies of last point for numerical reasons
            
            ! set up resonating surface for each mode
            x_plot(1:size(res_surf,1),2) = res_surf(:,2)                        ! normal position of resonating mode
            x_plot(size(res_surf,1)+1:n_mod_tot,2) = &
                &res_surf(size(res_surf,1),2)                                   ! identical copies of last point for numerical reasons
            y_plot(1:size(res_surf,1),2) = res_surf(:,1)                        ! total mode index of resonating mode
            y_plot(size(res_surf,1)+1:n_mod_tot,2) = &
                &res_surf(size(res_surf,1),1)                                   ! identical copies of last point for numerical reasons
            
            ! scale x_plot by normal factor
            x_plot = x_plot/norm_factor
            
            ! plot the maximum at midplane
            call print_ex_2D(plot_title,file_name,y_plot,x=x_plot,&
                &draw=.false.)
            
            ! draw plot in file
            call draw_ex(plot_title,file_name,2,1,.false.,extra_ops=&
                &'set xrange ['//&
                &trim(r2str(grid_sol_trim%r_F(1)/norm_factor))//':'//&
                &trim(r2str(grid_sol_trim%r_F(grid_sol_trim%n(3))/norm_factor))&
                &//']',ex_plot_style=1)
        end if
        
        ! clean up
        call grid_sol_trim%dealloc()
    end function plot_harmonics
    
    !> Decomposes  the plasma  potential and  kinetic energy  in its  individual
    !! terms for an individual Eigenvalue.
    !!
    !! Use  is made of  variables representing  the potential and  kinetic energy
    !! \cite Weyens3D.
    !!  - \f$E_\text{pot}\f$:
    !!      - normal line bending term:
    !!          \f$\frac{1}{\mu_0} \frac{Q_n^2}{h^{22}}\f$,
    !!      - geodesic line bending term:
    !!          \f$\frac{1}{\mu_0} \mathcal{J}^2 \frac{h^{22}}{g_{33}} Q_g^2\f$,
    !!      - normal ballooning term:
    !!          \f$-2 p' X^2 \kappa_n\f$,
    !!      - geodesic ballooning term:
    !!          \f$-2 p' X U^* \kappa_g\f$,
    !!      - normal kink term:
    !!          \f$-\sigma X^*Q_g\f$,
    !!      - geodesic kink term:
    !!          \f$\sigma U^*Q_n\f$,
    !!  - \f$E_\text{kin}\f$:
    !!      - normal kinetic term:
    !!          \f$\rho \frac{X^2}{h^{22}}\f$,
    !!      - geodesic kinetic term:
    !!          \f$\rho \mathcal{J}^2 \frac{h^{22}}{g_{33}} U^2\f$.
    !!
    !! The energy terms  are calculated normally on the  sol grid, interpolating
    !! the quantities defined on the eq grid, and angularly in the eq grid.
    !!
    !! Optionally,  the results  can be  plotted  by providing  X, Y  and Z.  By
    !! default, they are instead written to an output file.
    !!
    !! Also, the fraction between potential  and kinetic energy can be returned,
    !! compared with the eigenvalue.
    !!
    !! \return ierr
    integer function decompose_energy(mds_sol,grid_eq,grid_X,grid_sol,eq_1,&
        &eq_2,X,sol,vac,X_id,B_aligned,XYZ,E_pot_int,E_kin_int) result(ierr)
        
        use grid_utilities, only: trim_grid
        use num_vars, only: no_plots, eq_job_nr, eq_jobs_lims, eq_job_nr
        
        character(*), parameter :: rout_name = 'decompose_energy'
        
        ! input / output
        type(modes_type), intent(in) :: mds_sol                                 !< general modes variables in solution grid
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(in) :: X                                         !< perturbation variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue
        logical, intent(in) :: B_aligned                                        !< whether grid is field-aligned
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          !< X, Y and Z for plotting
        complex(dp), intent(inout), optional :: E_pot_int(7)                    !< integrated potential energy
        complex(dp), intent(inout), optional :: E_kin_int(2)                    !< integrated kinetic energy
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed sol grid
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: kd                                                           ! counter
        integer :: loc_dim(3)                                                   ! local dimension
        integer :: plot_dim(3)                                                  ! total dimensions for plot
        integer :: plot_offset(3)                                               ! local offsets for plot
        logical :: cont_plot                                                    ! continued plot
        complex(dp), allocatable, target :: E_pot(:,:,:,:)                      ! potential energy
        complex(dp), allocatable, target :: E_kin(:,:,:,:)                      ! kinetic energy
        complex(dp) :: E_pot_int_loc(7)                                         ! integrated potential energy for this parallel job
        complex(dp) :: E_kin_int_loc(2)                                         ! integrated kinetic energy for this parallel job
        complex(dp), pointer :: E_pot_trim(:,:,:,:) => null()                   ! trimmed part of E_pot
        complex(dp), pointer :: E_kin_trim(:,:,:,:) => null()                   ! trimmed part of E_kin
        real(dp), allocatable, target :: X_tot(:,:,:,:), Y_tot(:,:,:,:), &
            &Z_tot(:,:,:,:)                                                     ! multiple copies of X, Y and Z for collection plot
        real(dp), pointer :: X_tot_trim(:,:,:,:) => null()                      ! trimmed part of X
        real(dp), pointer :: Y_tot_trim(:,:,:,:) => null()                      ! trimmed part of Y
        real(dp), pointer :: Z_tot_trim(:,:,:,:) => null()                      ! trimmed part of Z
        character(len=max_str_ln), allocatable :: var_names_pot(:)              ! name of potential energy variables
        character(len=max_str_ln), allocatable :: var_names_kin(:)              ! name of kinetic energy variables
        character(len=max_str_ln), allocatable :: var_names(:)                  ! name of other variables that are plot
        character(len=max_str_ln) :: file_name                                  ! name of file
        character(len=max_str_ln) :: description                                ! description
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate energy terms')
        call lvl_ud(1)
        
        ! calculate for this parallel job
        ierr = calc_E(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,X,sol,vac,&
            &B_aligned,X_id,E_pot,E_kin,E_pot_int_loc,E_kin_int_loc)
        CHCKERR('')
        
        ! add to totals if requested
        if (present(E_pot_int)) E_pot_int = E_pot_int + E_pot_int_loc
        if (present(E_kin_int)) E_kin_int = E_kin_int + E_kin_int_loc
        
        call lvl_ud(-1)
        
        ! plot if wanted
        if (present(XYZ)) then
            ! bypass plots if no_plots
            if (no_plots) return
            
            ! user output
            call writo('Preparing plots')
            call lvl_ud(1)
            
            ! set up continued plot
            cont_plot = eq_job_nr.gt.1
            
            ! set up potential energy variable names
            allocate(var_names_pot(6))
            var_names_pot(1) = '1. normal line bending term ~ Q_n^2'
            var_names_pot(2) = '2. geodesic line bending term ~ Q_g^2'
            var_names_pot(3) = '3. normal ballooning term ~ -p'' X^2 kappa_n'
            var_names_pot(4) = '4. geodesic ballooning term ~ -p'' X U* kappa_g'
            var_names_pot(5) = '5. normal kink term ~ -sigma X*Q_g'
            var_names_pot(6) = '6. geodesic kink term ~ sigma U*Q_n'
            
            ! set up kinetic energy variable names
            allocate(var_names_kin(2))
            var_names_kin(1) = '1. normal kinetic term ~ rho X^2'
            var_names_kin(2) = '2. geodesic kinetic term ~ rho U^2'
            
            ! trim sol grid
            ierr = trim_grid(grid_sol,grid_sol_trim,norm_id)
            CHCKERR('')
            
            ! set plot variables
            loc_dim = [grid_eq%n(1),grid_eq%n(2),grid_sol_trim%loc_n_r]
            plot_dim = [grid_eq%n(1),grid_eq%n(2),grid_sol_trim%n(3)]
            plot_offset = [0,0,grid_sol_trim%i_min-1]
            allocate(X_tot(loc_dim(1),loc_dim(2),loc_dim(3),6))
            allocate(Y_tot(loc_dim(1),loc_dim(2),loc_dim(3),6))
            allocate(Z_tot(loc_dim(1),loc_dim(2),loc_dim(3),6))
            do kd = 1,6
                X_tot(:,:,:,kd) = XYZ(:,:,norm_id(1):norm_id(2),1)
                Y_tot(:,:,:,kd) = XYZ(:,:,norm_id(1):norm_id(2),2)
                Z_tot(:,:,:,kd) = XYZ(:,:,norm_id(1):norm_id(2),3)
            end do
            allocate(var_names(2))
            
            ! possibly modify if multiple equilibrium parallel jobs
            if (size(eq_jobs_lims,2).gt.1) then
                plot_dim(1) = eq_jobs_lims(2,size(eq_jobs_lims,2)) - &
                    &eq_jobs_lims(1,1) + 1
                plot_offset(1) = eq_jobs_lims(1,eq_job_nr) - 1
            end if
            
            ! point to the trimmed versions of energy
            E_pot_trim => E_pot(:,:,norm_id(1):norm_id(2),:)
            E_kin_trim => E_kin(:,:,norm_id(1):norm_id(2),:)
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Plot energy terms')
            call lvl_ud(1)
            
            ! E_kin
            file_name = trim(i2str(X_id))//'_E_kin'
            description = 'Kinetic energy constituents'
            X_tot_trim => X_tot(:,:,:,1:2)
            Y_tot_trim => Y_tot(:,:,:,1:2)
            Z_tot_trim => Z_tot(:,:,:,1:2)
            call plot_HDF5(var_names_kin,trim(file_name)//'_RE',rp(E_kin_trim),&
                &tot_dim=[plot_dim,2],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,cont_plot=cont_plot,&
                &descr=description)
            call plot_HDF5(var_names_kin,trim(file_name)//'_IM',ip(E_kin_trim),&
                &tot_dim=[plot_dim,2],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,cont_plot=cont_plot,&
                &descr=description)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            
            ! E_pot
            file_name = trim(i2str(X_id))//'_E_pot'
            description = 'Potential energy constituents'
            X_tot_trim => X_tot(:,:,:,1:6)
            Y_tot_trim => Y_tot(:,:,:,1:6)
            Z_tot_trim => Z_tot(:,:,:,1:6)
            call plot_HDF5(var_names_pot,trim(file_name)//'_RE',rp(E_pot_trim),&
                &tot_dim=[plot_dim,6],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,cont_plot=cont_plot,&
                &descr=description)
            call plot_HDF5(var_names_pot,trim(file_name)//'_IM',ip(E_pot_trim),&
                &tot_dim=[plot_dim,6],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,cont_plot=cont_plot,&
                &descr=description)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            
            ! E_stab
            var_names(1) = '1. stabilizing term'
            var_names(2) = '2. destabilizing term'
            file_name = trim(i2str(X_id))//'_E_stab'
            description = '(de)stabilizing energy'
            X_tot_trim => X_tot(:,:,:,1:2)
            Y_tot_trim => Y_tot(:,:,:,1:2)
            Z_tot_trim => Z_tot(:,:,:,1:2)
            call plot_HDF5(var_names,trim(file_name)//'_RE',rp(reshape(&
                &[sum(E_pot_trim(:,:,:,1:2),4),sum(E_pot_trim(:,:,:,3:6),4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[plot_dim,2],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,&
                &cont_plot=cont_plot,descr=description)
            call plot_HDF5(var_names,trim(file_name)//'_IM',ip(reshape(&
                &[sum(E_pot_trim(:,:,:,1:2),4),sum(E_pot_trim(:,:,:,3:6),4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[plot_dim,2],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,&
                &cont_plot=cont_plot,descr=description)
            
            ! E
            var_names(1) = '1. potential energy'
            var_names(2) = '2. kinetic energy'
            file_name = trim(i2str(X_id))//'_E'
            description = 'total potential and kinetic energy'
            call plot_HDF5(var_names,trim(file_name)//'_RE',rp(reshape(&
                &[sum(E_pot_trim,4),sum(E_kin_trim,4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[plot_dim,2],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,cont_plot=cont_plot,&
                &descr=description)
            call plot_HDF5(var_names,trim(file_name)//'_IM',ip(reshape(&
                &[sum(E_pot_trim,4),sum(E_kin_trim,4)],&
                &[loc_dim(1),loc_dim(2),loc_dim(3),2])),&
                &tot_dim=[plot_dim,2],loc_offset=[plot_offset,0],&
                &X=X_tot_trim,Y=Y_tot_trim,Z=Z_tot_trim,cont_plot=cont_plot,&
                &descr=description)
            
            call lvl_ud(-1)
            
            ! clean up
            nullify(E_kin_trim,E_pot_trim)
            nullify(X_tot_trim,Y_tot_trim,Z_tot_trim)
            call grid_sol_trim%dealloc()
        end if
    end function decompose_energy
    
    !> Calculate the energy terms in the energy decomposition.
    !!
    !! \note see explanation of routine in decompose_energy().
    !!
    !! \return ierr
    integer function calc_E(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,X,sol,&
        &vac,B_aligned,X_id,E_pot,E_kin,E_pot_int,E_kin_int) result(ierr)
        
        use num_vars, only: use_pol_flux_F, n_procs, K_style, &
            &norm_disc_prec_sol, rank, eq_job_nr, eq_jobs_lims
        use eq_vars, only: vac_perm
        use num_utilities, only: c
        use grid_utilities, only: calc_int_vol, trim_grid, untrim_grid, &
            &setup_interp_data, apply_disc
        use grid_vars, only: alpha, n_alpha
        use MPI_utilities, only: get_ser_var
        use sol_utilities, only: calc_XUQ
#if ldebug
        use grid_utilities, only: setup_deriv_data, apply_disc
#endif
        
        character(*), parameter :: rout_name = 'calc_E'
        
        ! input / output
        type(modes_type), intent(in) :: mds_sol                                 !< general modes variables in solution grid
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 !< and solution grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(in) :: X                                         !< perturbation variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        type(vac_type), intent(in) :: vac                                       !< vacuum variables
        logical, intent(in) :: B_aligned                                        !< whether grid is field-aligned
        integer, intent(in) :: X_id                                             !< nr. of Eigenvalue
        complex(dp), intent(inout), allocatable :: E_pot(:,:,:,:)               !< potential energy
        complex(dp), intent(inout), allocatable :: E_kin(:,:,:,:)               !< kinetic energy
        complex(dp), intent(inout) :: E_pot_int(7)                              !< integrated potential energy
        complex(dp), intent(inout) :: E_kin_int(2)                              !< integrated kinetic energy
        
        ! local variables
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: jd, kd                                                       ! counter
        integer :: loc_dim(3)                                                   ! local dimension
        type(grid_type) :: grid_sol_trim                                        ! trimmed sol grid
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        real(dp), allocatable :: h22(:,:,:), g33(:,:,:), J(:,:,:)               ! interpolated h_FD(2,2), g_FD(3,3) and J_FD
        real(dp), allocatable :: kappa_n(:,:,:), kappa_g(:,:,:)                 ! interpolated kappa_n and kappa_g
        real(dp), allocatable :: sigma(:,:,:)                                   ! interpolated sigma
        real(dp), allocatable :: D2p(:,:,:), rho(:,:,:)                         ! interpolated Dpres_FD, rho
        real(dp), allocatable :: ang_1(:,:,:), ang_2(:,:,:), norm(:)            ! coordinates of grid in which to calculate total energy
        complex(dp), allocatable :: XUQ(:,:,:,:)                                ! X, U, Q_n and Q_g
        complex(dp), allocatable :: E_int_tot(:)                                ! integrated potential or kinetic energy of all processes
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        real(dp), allocatable :: S(:,:,:)                                       ! interpolated S
        complex(dp), allocatable :: DU(:,:,:)                                   ! D_par U
        complex(dp), allocatable :: DU_ALT(:,:,:)                               ! DU calculated from U
        type(disc_type) :: ang_1_deriv_data                                     ! deriv data for ang_1
#endif
        
        ! set loc_dim
        loc_dim = [grid_X%n(1:2),grid_sol%loc_n_r]                              ! includes ghost regions of width norm_disc_prec_sol
        
        ! allocate local variables
        allocate(h22(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(g33(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(J(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(kappa_n(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(kappa_g(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(sigma(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(D2p(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(rho(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(ang_1(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(ang_2(loc_dim(1),loc_dim(2),loc_dim(3)))
        allocate(norm(loc_dim(3)))
        allocate(XUQ(loc_dim(1),loc_dim(2),loc_dim(3),4))
        allocate(E_pot(loc_dim(1),loc_dim(2),loc_dim(3),6))
        allocate(E_kin(loc_dim(1),loc_dim(2),loc_dim(3),2))
#if ldebug
        if (debug_calc_E .or. debug_DU) then
            allocate(DU(loc_dim(1),loc_dim(2),loc_dim(3)))
            allocate(S(loc_dim(1),loc_dim(2),loc_dim(3)))
            
            ! calculate D_par U
            ierr = calc_XUQ(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,X,sol,&
                &X_id,2,0._dp,DU,deriv=.true.)
            CHCKERR('')
        end if
#endif
        
#if lwith_intel
        ! set J to zero to avoid INTEL compiler bug
        J = 0._dp
#endif
        
        ! setup normal interpolation data for equilibrium grid
        ierr = setup_interp_data(grid_eq%loc_r_F,grid_sol%loc_r_F,&
            &norm_interp_data,norm_disc_prec_sol)
        CHCKERR('')
        
        ! interpolate
        ierr = apply_disc(eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0),&
            &norm_interp_data,h22,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0),&
            &norm_interp_data,g33,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%jac_FD(:,:,:,0,0,0),&
            &norm_interp_data,J,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%kappa_n,norm_interp_data,kappa_n,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%kappa_g,norm_interp_data,kappa_g,3)
        CHCKERR('')
        ierr = apply_disc(eq_2%sigma,norm_interp_data,sigma,3)
        CHCKERR('')
        ierr = apply_disc(eq_1%pres_FD(:,1),norm_interp_data,D2p(1,1,:))
        CHCKERR('')
        ierr = apply_disc(eq_1%rho,norm_interp_data,rho(1,1,:))
        CHCKERR('')
        do kd = 1,loc_dim(3)
            D2p(:,:,kd) = D2p(1,1,kd)
            rho(:,:,kd) = rho(1,1,kd)
        end do
#if ldebug
        if (debug_calc_E) then
            ierr = apply_disc(eq_2%S,norm_interp_data,S,3)
            CHCKERR('')
        end if
#endif
        ! clean up
        call norm_interp_data%dealloc()
        
        ! set angles and norm
        if (B_aligned) then
            if (use_pol_flux_F) then
                ang_1 = grid_X%theta_F
            else
                ang_1 = grid_X%zeta_F
            end if
            do jd = 1,n_alpha
                ang_2(:,jd,:) = alpha(jd)
            end do
        else
            ang_1 = grid_X%theta_F
            ang_2 = grid_X%zeta_F
        end if
        norm = grid_X%loc_r_F
        
        ! calculate X, U, Q_n and Q_g
        do kd = 1,4
            ierr = calc_XUQ(mds_sol,grid_eq,grid_X,grid_sol,eq_1,eq_2,X,sol,&
                &X_id,kd,0._dp,XUQ(:,:,:,kd))
            CHCKERR('')
        end do
        
        ! calc kinetic energy
        E_kin(:,:,:,1) = rho/h22*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
        select case (K_style)
            case (1)                                                            ! normalization of full perpendicular component
                E_kin(:,:,:,2) = &
                    &rho*h22*J**2/g33*XUQ(:,:,:,2)*conjg(XUQ(:,:,:,2))
            case (2)                                                            ! normalization of only normal component
                E_kin(:,:,:,2) = 0._dp
            case default
                err_msg = 'No normalization style associated with '//&
                    &trim(i2str(K_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! calc potential energy
        E_pot(:,:,:,1) = 1._dp/vac_perm*1._dp/h22*XUQ(:,:,:,3)*&
            &conjg(XUQ(:,:,:,3))
        E_pot(:,:,:,2) = 1._dp/vac_perm*h22*J**2/g33*XUQ(:,:,:,4)*&
            &conjg(XUQ(:,:,:,4))
        E_pot(:,:,:,3) = -2*D2p*kappa_n*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
        E_pot(:,:,:,4) = -2*D2p*kappa_g*XUQ(:,:,:,1)*conjg(XUQ(:,:,:,2))
        E_pot(:,:,:,5) = -sigma*XUQ(:,:,:,4)*conjg(XUQ(:,:,:,1))
        E_pot(:,:,:,6) = sigma*XUQ(:,:,:,3)*conjg(XUQ(:,:,:,2))
        
#if ldebug
        if (debug_calc_E) then
            call writo('Testing whether the total potential energy is &
                &equal to the alternative form given in [ADD REF]')
            call lvl_ud(1)
            
            ! alternative formulation for E_pot, always real
            E_pot(:,:,:,3) = (-2*D2p*kappa_n+sigma*S)*&
                &XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))
            E_pot(:,:,:,4) = -sigma/J*2*rp(XUQ(:,:,:,1)*conjg(DU))
            E_pot(:,:,:,5:6) = 0._dp
            
            call lvl_ud(-1)
        end if
        
        if (debug_X_norm) then
            call writo('Plotting the squared amplitude of the normal &
                &component')
            call lvl_ud(1)
            
            call plot_HDF5('X_norm', 'TEST_X_norm_POST_'//trim(i2str(X_id)),&
                &rp(XUQ(:,:,:,1)*conjg(XUQ(:,:,:,1))),&
                &tot_dim=grid_X%n,loc_offset=[0,0,grid_X%i_min-1])
            
            call lvl_ud(-1)
        end if
        
        if (debug_DU .and. B_aligned) then
            call writo('Testing whether DU is indeed the parallel &
                &derivative of U')
            call lvl_ud(1)
            
            allocate(DU_ALT(loc_dim(1),loc_dim(2),loc_dim(3)))
            
            ! derive
            do kd = 1,loc_dim(3)
                do jd = 1,loc_dim(2)
                    ierr = setup_deriv_data(ang_1(:,jd,kd),ang_1_deriv_data,1,&
                        &norm_disc_prec_sol)
                    CHCKERR('')
                    ierr = apply_disc(XUQ(:,jd,kd,2),ang_1_deriv_data,&
                        &DU_ALT(:,jd,kd))
                    CHCKERR('')
                end do
            end do
            call ang_1_deriv_data%dealloc()
            
            ! plot real part
            call plot_HDF5('RE U','TEST_RE_U_'//&
                &trim(i2str(X_id)),rp(XUQ(:,:,:,2)),grid_X%n,&
                &[0,0,grid_sol%i_min-1])
            call plot_diff_HDF5(rp(DU),rp(DU_ALT),&
                &'TEST_RE_DU_'//trim(i2str(X_id)),grid_X%n,&
                &[0,0,grid_sol%i_min-1],descr='To test whether DU is &
                &parallel derivative of U',output_message=.true.)
            
            ! plot imaginary part
            call plot_HDF5('IM U','TEST_IM_U_'//&
                &trim(i2str(X_id)),ip(XUQ(:,:,:,2)),grid_X%n,&
                &[0,0,grid_sol%i_min-1])
            call plot_diff_HDF5(ip(DU),ip(DU_ALT),&
                &'TEST_IM_DU_'//trim(i2str(X_id)),grid_X%n,&
                &[0,0,grid_sol%i_min-1],descr='To test whether DU is &
                &parallel derivative of U',output_message=.true.)
            
            deallocate(DU_ALT)
            
            call lvl_ud(-1)
        end if
#endif
        
        ! trim sol grid
        ierr = trim_grid(grid_sol,grid_sol_trim,norm_id)
        CHCKERR('')
        
        ! add ghost region of width one to the right of the interval
        if (rank.lt.n_procs-1) norm_id(2) = norm_id(2)+1                        ! adjust norm_id as well
        
        ! integrate energy using ghosted variables
        ierr = calc_int_vol(ang_1(:,:,norm_id(1):norm_id(2)),&
            &ang_2(:,:,norm_id(1):norm_id(2)),norm(norm_id(1):norm_id(2)),&
            &J(:,:,norm_id(1):norm_id(2)),&
            &E_kin(:,:,norm_id(1):norm_id(2),:),E_kin_int)
        CHCKERR('')
        ierr = calc_int_vol(ang_1(:,:,norm_id(1):norm_id(2)),&
            &ang_2(:,:,norm_id(1):norm_id(2)),norm(norm_id(1):norm_id(2)),&
            &J(:,:,norm_id(1):norm_id(2)),&
            &E_pot(:,:,norm_id(1):norm_id(2),:),E_pot_int(1:6))
        CHCKERR('')
        
        ! add vacuum energy if last process and first equilibrium job
        E_pot_int(7) = 0._dp
        write(*,*) '¡¡¡¡¡ NO VACUUM !!!!!'
        if (1.eq.2 .and. rank.eq.n_procs-1 .and. eq_job_nr.eq.size(eq_jobs_lims,2)) then
            do jd = 1,sol%n_mod
                do kd = 1,sol%n_mod
                    E_pot_int(7) = E_pot_int(7) + &
                        &conjg(sol%vec(kd,grid_sol%loc_n_r,X_id)) * &
                        &vac%res(kd,jd) * sol%vec(jd,grid_sol%loc_n_r,X_id)
                end do
            end do
        end if
        
        ! bundle all processes
        do kd = 1,2
            ierr = get_ser_var([E_kin_int(kd)],E_int_tot,scatter=.true.)
            CHCKERR('')
            E_kin_int(kd) = sum(E_int_tot)
            deallocate(E_int_tot)
        end do
        do kd = 1,7
            ierr = get_ser_var([E_pot_int(kd)],E_int_tot,scatter=.true.)
            CHCKERR('')
            E_pot_int(kd) = sum(E_int_tot)
            deallocate(E_int_tot)
        end do
        
        ! deallocate variables
        call grid_sol_trim%dealloc()
    end function calc_E
    
    !> Print solution quantities to an output file:
    !!
    !!  - sol:
    !!      - \c val
    !!      - \c vec
    !!
    !! If \c  rich_lvl is provided, \c  "_[R_rich_lvl]" is appended to  the data
    !! name if it is \c >0.
    !!
    !! \return ierr
    integer function print_output_sol(grid,sol,data_name,rich_lvl) result(ierr)
        use num_vars, only: PB3D_name, sol_n_procs, rank
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'print_output_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< solution grid variables
        type(sol_type), intent(in) :: sol                                       !< solution variables
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to print
        
        ! local variables
        type(var_1D_type), allocatable, target :: sol_1D(:)                     ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: sol_1D_loc => null()                      ! local element in sol_1D
        type(grid_type) :: grid_trim                                            ! trimmed solution grid
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! only write to output file if at least one solution
        if (size(sol%val).gt.0) then
            ! user output
            call writo('Write solution variables to output file')
            call lvl_ud(1)
            
            ! trim grids
            ierr = trim_grid(grid,grid_trim)
            CHCKERR('')
            
            ! Set up the 1D equivalents of the solution variables
            allocate(sol_1D(max_dim_var_1D))
            
            id = 1
            
            ! RE_sol_val
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'RE_sol_val'
            allocate(sol_1D_loc%tot_i_min(1),sol_1D_loc%tot_i_max(1))
            allocate(sol_1D_loc%loc_i_min(1),sol_1D_loc%loc_i_max(1))
            sol_1D_loc%tot_i_min = [1]
            sol_1D_loc%tot_i_max = [size(sol%val)]
            sol_1D_loc%loc_i_min = [1]
            sol_1D_loc%loc_i_max = [size(sol%val)]
            allocate(sol_1D_loc%p(size(sol%val)))
            sol_1D_loc%p = rp(sol%val)
            
            ! IM_sol_val
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'IM_sol_val'
            allocate(sol_1D_loc%tot_i_min(1),sol_1D_loc%tot_i_max(1))
            allocate(sol_1D_loc%loc_i_min(1),sol_1D_loc%loc_i_max(1))
            sol_1D_loc%tot_i_min = [1]
            sol_1D_loc%tot_i_max = [size(sol%val)]
            sol_1D_loc%loc_i_min = [1]
            sol_1D_loc%loc_i_max = [size(sol%val)]
            allocate(sol_1D_loc%p(size(sol%val)))
            sol_1D_loc%p = ip(sol%val)
            
            ! RE_sol_vec
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'RE_sol_vec'
            allocate(sol_1D_loc%tot_i_min(3),sol_1D_loc%tot_i_max(3))
            allocate(sol_1D_loc%loc_i_min(3),sol_1D_loc%loc_i_max(3))
            sol_1D_loc%loc_i_min = [1,grid_trim%i_min,1]
            sol_1D_loc%loc_i_max = [sol%n_mod,grid_trim%i_max,size(sol%vec,3)]
            sol_1D_loc%tot_i_min = [1,1,1]
            sol_1D_loc%tot_i_max = [sol%n_mod,grid_trim%n(3),size(sol%vec,3)]
            allocate(sol_1D_loc%p(size(sol%vec)))
            sol_1D_loc%p = reshape(rp(sol%vec),[size(sol%vec)])
            
            ! IM_sol_vec
            sol_1D_loc => sol_1D(id); id = id+1
            sol_1D_loc%var_name = 'IM_sol_vec'
            allocate(sol_1D_loc%tot_i_min(3),sol_1D_loc%tot_i_max(3))
            allocate(sol_1D_loc%loc_i_min(3),sol_1D_loc%loc_i_max(3))
            sol_1D_loc%loc_i_min = [1,grid_trim%i_min,1]
            sol_1D_loc%loc_i_max = [sol%n_mod,grid_trim%i_max,size(sol%vec,3)]
            sol_1D_loc%tot_i_min = [1,1,1]
            sol_1D_loc%tot_i_max = [sol%n_mod,grid_trim%n(3),size(sol%vec,3)]
            allocate(sol_1D_loc%p(size(sol%vec)))
            sol_1D_loc%p = reshape(ip(sol%vec),[size(sol%vec)])
            
            ! write
            ! Note: if processes that are not used  in the solve do not have the
            ! solution  vector and  should therefore  not be  used to  write the
            ! output.
            if (sol_n_procs.gt.1 .or. rank.eq.0) then
                ierr = print_HDF5_arrs(sol_1D(1:id-1),PB3D_name,trim(data_name),&
                    &rich_lvl=rich_lvl)
                CHCKERR('')
            end if
            
            ! clean up
            call grid_trim%dealloc()
            call dealloc_var_1D(sol_1D)
            nullify(sol_1D_loc)
            
            ! user output
            call lvl_ud(-1)
        else
            ! user output
            call writo('No solutions to write')
        end if
    end function print_output_sol
end module sol_ops
