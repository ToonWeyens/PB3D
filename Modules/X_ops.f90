!------------------------------------------------------------------------------!
!   Operations considering perturbation quantities                             !
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    use messages, only: lvl_ud, writo, print_ar_2
    use output_ops, only: print_GP_2D, print_GP_3D, draw_GP, print_HDF5
    use str_ops, only: i2str, r2strt, r2str

    implicit none
    private
    public prepare_X, solve_EV_system, plot_X_vec, calc_PV, calc_KV, calc_U, &
        &calc_extra, normalize_X_vars
    
    !! for testing:
    !character(len=max_str_ln), allocatable :: var_names(:)                      ! names of variables
    !character(len=max_str_ln), allocatable :: file_name                         ! name of file
    !integer :: tot_dim(4), grp_dim(4), grp_offset(4)                            ! total dim., group dim. and group offset

contains
    ! prepare the matrix elements by calculating  KV_i and PV_i, which then will
    ! have to be integrated, with a complex exponential weighting function
    integer function prepare_X(eq,grid,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, plot_jq, grp_nr
        use eq_vars, only: eq_type
        use grid_vars, onlY: grid_type
        use metric_vars, only: metric_type
        use X_vars, only: dealloc_X, X_type
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'prepare_X'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(grid_type) :: grid                                                 ! grid variables
        type(metric_type) :: met                                                ! metric variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables
        integer :: m, k                                                         ! counters
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle in flux coordinates
        character(len=5) :: ang_par_F_name                                      ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Preparing perturbation variables')
        call lvl_ud(1)
        
        ! tests
        ierr = check_modes(eq,X)
        CHCKERR('')
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_F) then
            ang_par_F => grid%theta_F
        else
            ang_par_F => grid%zeta_F
        end if
        
        ! plot resonances if requested
        if (plot_jq) then
            call writo('Resonance plot requested')
            call lvl_ud(1)
            if (grp_nr.eq.0) ierr = resonance_plot(eq,grid,X)
            CHCKERR('')
            call lvl_ud(-1)
            call writo('Resonance plot done')
        end if
        
        ! set exp_ang_par_F
        if (use_pol_flux_F) then
            do m = 1,X%n_mod
                do k = m,X%n_mod
                    X%exp_ang_par_F(:,:,:,c([k,m],.true.,X%n_mod)) = &
                        &exp(iu*(k-m)*grid%theta_F)
                end do
            end do
        else
            do m = 1,X%n_mod
                do k = m,X%n_mod
                    X%exp_ang_par_F(:,:,:,c([k,m],.true.,X%n_mod)) = &
                        &exp(iu*(m-k)*grid%zeta_F)
                end do
            end do
        end if
        
        ! set up ang_par_F_name
        if (use_pol_flux_F) then
            ang_par_F_name = 'theta'
        else
            ang_par_F_name = 'zeta'
        end if
        
        call writo('Calculating table of field line averages &
            &<V e^i[(k-m)'//trim(ang_par_F_name)//']>')
        
        call lvl_ud(1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        ierr = calc_U(eq,grid,met,X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! calculate extra equilibrium quantities
        call writo('Calculating extra equilibrium quantities...')
        call lvl_ud(1)
        ierr = calc_extra(eq,grid,met,X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! Calculate PV_i for all (k,m) pairs and n_r (equilibrium) values of the
        ! normal coordinate
        call writo('Calculating PV...')
        call lvl_ud(1)
        call calc_PV(eq,grid,met,X)
        call lvl_ud(-1)
        
        ! Calculate KV_i for all (k,m) pairs and n_r (equilibrium) values of the
        ! normal coordinate
        call writo('Calculating KV...')
        call lvl_ud(1)
        call calc_KV(eq,grid,met,X)
        call lvl_ud(-1)
        
        ! Calculate PV_int = <PV e^(k-m)ang_par_F>
        call writo('Taking field average of PV')
        call lvl_ud(1)
        ierr = calc_V_int(grid,met,X%exp_ang_par_F,X%n_mod,X%PV_0(:,:,:,:),&
            &X%PV_int_0(:,:,:))
        CHCKERR('')
        ierr = calc_V_int(grid,met,X%exp_ang_par_F,X%n_mod,X%PV_1(:,:,:,:),&
            &X%PV_int_1(:,:,:))
        CHCKERR('')
        ierr = calc_V_int(grid,met,X%exp_ang_par_F,X%n_mod,X%PV_2(:,:,:,:),&
            &X%PV_int_2(:,:,:))
        CHCKERR('')
        call lvl_ud(-1)
        
        ! Calculate KV_int = <KV e^(k-m)ang_par_F>
        call writo('Taking field average of KV')
        call lvl_ud(1)
        ierr = calc_V_int(grid,met,X%exp_ang_par_F,X%n_mod,X%KV_0(:,:,:,:),&
            &X%KV_int_0(:,:,:))
        CHCKERR('')
        ierr = calc_V_int(grid,met,X%exp_ang_par_F,X%n_mod,X%KV_1(:,:,:,:),&
            &X%KV_int_1(:,:,:))
        CHCKERR('')
        ierr = calc_V_int(grid,met,X%exp_ang_par_F,X%n_mod,X%KV_2(:,:,:,:),&
            &X%KV_int_2(:,:,:))
        CHCKERR('')
        call lvl_ud(-1)
        
        ! deallocate equilibrium variables
        call writo('deallocating unused variables')
        call dealloc_X(X)
        
        call lvl_ud(-1)
        
        call writo('Done calculating')
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation variables prepared')
    end function prepare_X
    
    ! plot  q-profile  or iota-profile  in  flux coordinates  with nq-m  = 0  or
    ! n-iotam = 0 indicate if requested
    ! [MPI] Parts by all processes, parts only by global master
    integer function resonance_plot(eq,grid,X) result(ierr)
        use num_vars, only: use_pol_flux_E, use_pol_flux_F, output_style, &
            &grp_n_procs, grp_rank, tol_NR, no_plots
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        use utilities, only: calc_zero_NR, interp_fun
        use MPI_ops, only: get_ser_var
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(grid_type) :: grid                                                 ! grid variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables (also used in child functions)
        real(dp) :: mnfrac_for_function                                         ! fraction m/n or n/m to determine resonant flux surface
        real(dp), allocatable :: jq_for_function(:,:)                           ! rot_t or q_saf
        
        ! local variables (not to be used in child functions)
        integer :: jd, kd                                                       ! counters
        real(dp), allocatable :: x_vars(:,:)                                    ! for plotting
        real(dp), allocatable :: y_vars(:,:)                                    ! for plotting
        real(dp) :: jq_solution                                                 ! solution for q = m/n or iota = n/m
        real(dp) :: jq_solution_transf                                          ! transformed solution to flux coordinates
        integer :: istat                                                        ! status
        character(len=max_str_ln) :: plot_title, file_name                      ! name of plot, of file
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: jq(:,:)                                        ! saf. fac. or rot. transf. in Flxu coords.
        real(dp), allocatable :: jq_loc(:)                                      ! local version of jq
        real(dp), allocatable :: flux_X(:)                                      ! pol. or tor. perturbation flux in Flux coords.
        real(dp), allocatable :: flux_eq(:)                                     ! pol. or tor. equilibrium flux in Flux coords.
        integer, allocatable :: tot_i_min(:)                                    ! i_min of Equilibrium grid of all processes
        integer :: i_lim(2)                                                     ! limits of indices in group arrays
        integer :: n_r                                                          ! total number of normal points
        real(dp) :: tol_NR_old                                                  ! old value of tol_NR
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! get min_i of equilibrium grid
        ierr = get_ser_var([grid%i_min],tot_i_min,scatter=.true.)
        CHCKERR('')
        
        ! set index limits
        i_lim(1) = 1                                                            ! start with first index of this group
        if (grp_rank.lt.grp_n_procs-1) then                                     ! not last process
            i_lim(2) = tot_i_min(grp_rank+2)-tot_i_min(grp_rank+1)              ! end with one before the start of next group
        else                                                                    ! last process
            i_lim(2) = grid%grp_n_r                                             ! end of this last group
        end if
        
        ! get serial version of flux_X and safety factor or rot. transform
        if (use_pol_flux_F) then
            ierr = get_ser_var(eq%flux_p_FD(i_lim(1):i_lim(2),0),flux_X)
            CHCKERR('')
            if (grp_rank.eq.0) allocate(jq(size(flux_X),0:2))
            do jd = 0,2
                ierr = get_ser_var(eq%q_saf_FD(i_lim(1):i_lim(2),jd),jq_loc)
                CHCKERR('')
                if(grp_rank.eq.0) jq(:,jd) = jq_loc*X%max_flux_F**jd
            end do
        else
            ierr = get_ser_var(eq%flux_t_FD(i_lim(1):i_lim(2),0),flux_X)
            CHCKERR('')
            if (grp_rank.eq.0) allocate(jq(size(flux_X),0:2))
            do jd = 0,2
                ierr = get_ser_var(eq%rot_t_FD(i_lim(1):i_lim(2),jd),jq_loc)
                CHCKERR('')
                if(grp_rank.eq.0) jq(:,jd) = jq_loc*X%max_flux_F**jd
            end do
        end if
        if (use_pol_flux_E) then
            ierr = get_ser_var(eq%flux_p_FD(i_lim(1):i_lim(2),0),flux_eq)
            CHCKERR('')
        else
            ierr = get_ser_var(eq%flux_t_FD(i_lim(1):i_lim(2),0),flux_eq)
            CHCKERR('')
        end if
        
        if (use_pol_flux_F) then
            call writo('Plotting safety factor q and resonant surfaces &
                &q = m/n')
        else
            call writo('Plotting rotational transform iota and resonant &
                &surfaces iota = n/m')
        end if
        call lvl_ud(1)
            
        ! the rest is done only by global master
        if (grp_rank.eq.0) then
            call writo('calculating resonant surfaces')
            
            ! initialize variables
            n_r = size(flux_X)
            allocate(x_vars(n_r,X%n_mod+1)); x_vars = 0
            allocate(y_vars(n_r,X%n_mod+1)); y_vars = 0
            allocate(jq_for_function(n_r,0:2))
            
            ! set x_vars and y_vars for first column
            x_vars(:,1) = flux_X/abs(flux_X(n_r))
            y_vars(:,1) = jq(:,0)
            jq_for_function = jq
            
            ! save old tol_NR
            tol_NR_old = tol_NR
            
            ! change tol_NR
            tol_NR = 1.E-4_dp
            
            ! loop over all modes (and shift the index in x and y_vars by 1)
            kd = 2
            do jd = 1, X%n_mod
                ! find place where q = m/n or  iota = n/m in Flux coordinates by
                ! solving q-m/n = 0 or iota-n/m=0, using the functin jq_fun
                call lvl_ud(1)
                
                ! set up mnfrac_for_function
                if (use_pol_flux_F) then
                    mnfrac_for_function = X%m(jd)*1.0_dp/X%n(jd)
                else
                    mnfrac_for_function = X%n(jd)*1.0_dp/X%m(jd)
                end if
                
                ! calculate zero using Newton-Rhapson
                istat = calc_zero_NR(jq_solution,jq_fun,jq_dfun,1.0_dp)
                
                ! intercept error
                if (istat.ne.0) then
                    call writo('Error intercepted: Couldn''t find resonating &
                        &surface for (n,m) = ('//trim(i2str(X%n(jd)))//','//&
                        &trim(i2str(X%m(jd)))//')')
                else if (jq_solution.lt.0.0_dp) then
                    call writo('Mode (n,m) = ('//trim(i2str(X%n(jd)))//','//&
                        &trim(i2str(X%m(jd)))//') does not resonate in plasma')
                else
                    if (jq_solution.gt.1.0_dp) then
                        call writo('Mode (n,m) = ('//trim(i2str(X%n(jd)))//','&
                            &//trim(i2str(X%m(jd)))//') does not resonate &
                            &in plasma')
                        y_vars(n_r,kd) = jq(n_r,0)
                    else
                        ! convert solution to flux coordinates if GNUPlot output
                        if (output_style.eq.1) then
                            istat = interp_fun(jq_solution_transf,&
                                &flux_X/abs(flux_X(n_r)),jq_solution)
                            x_vars(:,kd) = jq_solution_transf
                            call writo('Mode (n,m) = ('//trim(i2str(X%n(jd)))//&
                                &','//trim(i2str(X%m(jd)))//') resonates in &
                                &plasma at normalized flux surface '//&
                                &trim(r2str(jq_solution_transf)))
                        else
                            x_vars(:,kd) = jq_solution
                        end if
                        ! the y axis is always in perturbation grid
                        if (use_pol_flux_F) then
                            y_vars(n_r,kd) = X%m(jd)*1.0_dp/X%n(jd)
                        else
                            y_vars(n_r,kd) = X%n(jd)*1.0_dp/X%m(jd)
                        end if
                    end if
                    kd = kd + 1
                end if
                
                call lvl_ud(-1)
            end do
            
            ! recover old tol_NR
            tol_NR = tol_NR_old
            
            ! check status
            if (istat.ne.0) then
                call writo('WARNING: Failed to produce plot')
                return
            end if
            
            call writo('Plotting results')
            if (use_pol_flux_F) then
                plot_title = 'safety factor q'
                file_name = 'q_saf'
            else
                plot_title = 'rotational transform iota'
                file_name = 'rot_t'
            end if
            select case(output_style)
                case(1)                                                             ! GNUPlot output
                    ! plot on screen
                    call print_GP_2D(plot_title,trim(file_name)//'.dat.',&
                        &y_vars(:,1:kd-1),x=x_vars(:,1:kd-1))
                    
                    ! plot in file as well
                    call draw_GP(plot_title,trim(file_name)//'.dat.',kd-1,&
                        &.true.,.false.)
                case(2)                                                             ! HDF5 output
                    ierr = resonance_plot_HDF5()
                    CHCKERR('')
                case default
                    err_msg = 'No style associated with '//&
                        &trim(i2str(output_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! deallocate local variables
            deallocate(flux_X,flux_eq,jq)
        
        end if
        
        call lvl_ud(-1)
        call writo('Done plotting')
    contains
        ! Returns q-m/n or  iota-n/m in Flux coordinates, used to  solve for q =
        ! m/n or iota = n/m.
        real(dp) function jq_fun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(n_r)                                              ! so interp_fun can be used
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.0) then                                                   ! point requested lower than 0
                ! extrapolate variable from 1
                res = jq_for_function(1,0) - mnfrac_for_function + &
                    &jq_for_function(1,1)*(pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from n_r
                res = jq_for_function(n_r,0) - mnfrac_for_function + &
                    &jq_for_function(n_r,1)*(pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                varin = jq_for_function(:,0) - mnfrac_for_function
                istat = interp_fun(res,varin,pt)
            end if
        end function jq_fun
        
        ! returns d(q-m/n)/dr in Flux coordinates, used to solve for q = m/n
        ! WARNING: This  routine requires that jq_for_function's  derivatives be
        ! calculated up to order 2. This is NOT checked!
        real(dp) function jq_dfun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(n_r)                                           ! so interp_fun can be used
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.0) then                                                   ! point requested lower than 0
                ! extrapolate variable from 1
                res = jq_for_function(1,1) + jq_for_function(1,2)*(pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from n_r
                res = jq_for_function(n_r,1) + jq_for_function(n_r,2)*&
                    &(pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                varin = jq_for_function(:,1)
                istat = interp_fun(res,varin,pt)
            end if
        end function jq_dfun
        
        ! plots the resonance plot in 3D in HDF5 format
        integer function resonance_plot_HDF5() result(ierr)
            use num_vars, only: n_theta_plot, n_zeta_plot
            use output_ops, only: print_HDF5
            use grid_vars, only: create_grid, destroy_grid
            use grid_ops, only: calc_XYZ_grid, calc_eqd_grid, coord_F2E
            
            character(*), parameter :: rout_name = 'resonance_plot_HDF5'
            
            ! local variables
            integer :: id                                                       ! counters
            real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)        ! pol. and tor. angle of plot
            real(dp), allocatable :: r_plot(:)                                  ! normal coordinates of plot
            real(dp), allocatable :: X_plot(:,:,:,:), Y_plot(:,:,:,:), &
                &Z_plot(:,:,:,:)                                                ! X, Y and Z of plot of all surfaces
            real(dp), allocatable :: X_plot_ind(:,:,:), Y_plot_ind(:,:,:), &
                &Z_plot_ind(:,:,:)                                              ! X, Y and Z of plots of individual surfaces
            integer :: plot_dim(4)                                              ! plot dimensions (total = group because only group masters)
            real(dp), allocatable :: vars(:,:,:,:)                              ! variable to plot
            type(grid_type) :: grid_plot                                        ! grid for plotting
            
            ! initialize ierr
            ierr = 0
            
            ! set up pol. and tor. angle for plot
            allocate(theta_plot(n_theta_plot,n_zeta_plot,1))
            ierr = calc_eqd_grid(theta_plot(:,1,1),0._dp,2*pi)
            CHCKERR('')
            do id = 2,n_zeta_plot
                theta_plot(:,id,1) = theta_plot(:,1,1)
            end do
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,1))
            ierr = calc_eqd_grid(zeta_plot(1,:,1),0._dp,2*pi)
            CHCKERR('')
            do id = 2,n_theta_plot
                zeta_plot(id,:,1) = zeta_plot(1,:,1)
            end do
            
            ! set up vars
            allocate(vars(n_theta_plot,n_zeta_plot,1,X%n_mod))
            do id = 1,X%n_mod
                vars(:,:,:,id) = y_vars(n_r,id+1)
            end do
            
            ! set dimensions
            plot_dim = [n_theta_plot,n_zeta_plot,1,X%n_mod]
            
            ! calculate normal coords. in Equilibrium coords.
            allocate(r_plot(X%n_mod))
            ierr = coord_F2E(eq,X,x_vars(n_r,2:X%n_mod+1),r_plot,&
                &r_F_array=flux_X/flux_X(n_r),r_E_array=flux_eq/flux_eq(n_r))
            CHCKERR('')
            
            ! create plot grid
            ierr = create_grid(grid_plot,plot_dim(1:3))
            CHCKERR('')
            grid_plot%theta_E = theta_plot
            grid_plot%zeta_E = zeta_plot
            
            ! set up plot X, Y and Z
            allocate(X_plot(n_theta_plot,n_zeta_plot,1,X%n_mod))
            allocate(Y_plot(n_theta_plot,n_zeta_plot,1,X%n_mod))
            allocate(Z_plot(n_theta_plot,n_zeta_plot,1,X%n_mod))
            
            ! loop over all resonant surfaces to calculate X, Y and Z values
            do id = 1,X%n_mod
                ! set grp_r_E of plot grid
                grid_plot%grp_r_E = r_plot(id)
                
                ierr = calc_XYZ_grid(grid_plot,X_plot_ind,Y_plot_ind,Z_plot_ind)
                CHCKERR('')
                
                ! save the individual variable in the total variables
                X_plot(:,:,:,id) = X_plot_ind
                Y_plot(:,:,:,id) = Y_plot_ind
                Z_plot(:,:,:,id) = Z_plot_ind
                
                ! deallocate the individual variables
                deallocate(X_plot_ind,Y_plot_ind,Z_plot_ind)
            end do
            
            ! print using HDF5
            call print_HDF5([plot_title],file_name,vars,plot_dim,plot_dim,&
                &[0,0,0,0],X=X_plot,Y=Y_plot,Z=Z_plot,col=3,&
                &description='resonant surfaces')
            
            ! deallocate local variables
            deallocate(vars)
            deallocate(theta_plot,zeta_plot,r_plot)
            
            ! delete plot grid
            call destroy_grid(grid_plot)
        end function resonance_plot_HDF5
    end function resonance_plot
    
    ! checks whether  nq-m <<  n or  n-iotam << m  is satisfied  in some  of the
    ! plasma
    ! [MPI] Parts by all processes, parts only by global master
    integer function check_modes(eq,X) result(ierr)
        use MPI_ops, only: get_ser_var
        use num_vars, only: glb_rank, use_pol_flux_F, eq_style
        use eq_vars, only: eq_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'check_modes'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp) :: tol = 0.2                                                   ! tolerance for being out of range of q or iota values
        real(dp) :: min_jq, max_jq                                              ! min. and max. values of q or iota
        integer :: pmone                                                        ! plus or minus one
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: ser_jq(:)                                      ! serial q_saf_E(:,0) or rot_t_E(:,0)
        
        ! initialize ierr
        ierr = 0
        
        ! get serial q_saf_E(:,0) or rot_t_E(:,0)
        ! (There will be some overlap but it does not matter)
        if (use_pol_flux_F) then
            ierr = get_ser_var(eq%q_saf_E(:,0),ser_jq)
            CHCKERR('')
        else
            ierr = get_ser_var(eq%rot_t_E(:,0),ser_jq)
            CHCKERR('')
        end if
        
        if (glb_rank.eq.0) then
            call writo('Checking mode numbers')
            ! set up plus minus one
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    pmone = -1                                                  ! conversion VMEC LH -> RH coord. system
                case (2)                                                        ! HELENA
                    pmone = 1
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! set min_jq and max_jq in flux coordinate system
            min_jq = minval(pmone*ser_jq)
            max_jq = maxval(pmone*ser_jq)
            
            ! multiply by plus minus one and tolerance
            min_jq = min_jq - tol*abs(min_jq)
            max_jq = max_jq + tol*abs(max_jq)
            
            ! for every mode (n,m) check whether  m/n is inside the range of
            ! q values or n/m inside the range of iota values
            do id = 1, X%n_mod
                if (use_pol_flux_F) then
                    if (X%m(id)*1.0/X%n(id) .lt.min_jq .or. &
                        &X%m(id)*1.0/X%n(id) .gt.max_jq) then
                        call writo('for (n,m) = ('//trim(i2str(X%n(id)))//&
                            &','//trim(i2str(X%m(id)))//'), the ratio &
                            &|n q - m|/n is never << 1')
                        ierr = 1
                        err_msg = 'Choose m and n so that |n q - m|/n << 1'
                        CHCKERR(err_msg)
                    end if
                else
                    if (X%n(id)*1.0/X%m(id) .lt.min_jq .or. &
                        &X%n(id)*1.0/X%m(id) .gt.max_jq) then
                        call writo('for (n,m) = ('//trim(i2str(X%n(id)))//&
                            &','//trim(i2str(X%m(id)))//'), the ratio &
                            &|n - iota m|/n is never << 1')
                        ierr = 1
                        err_msg = 'Choose m and n so that &
                            &|n - iota m|/n << 1'
                        CHCKERR(err_msg)
                    end if
                end if
            end do
            call writo('Mode numbers checked')
        end if
    end function check_modes
    
    ! set-up  and solve  the  EV system  by discretizing  the  equations in  the
    ! perturbation  grid,  making  use  of   PV  and  KV,  interpolated  in  the
    ! equilibrium grid.
    integer function solve_EV_system(grid_eq,grid_X,X,use_guess,&
        &n_sol_found) result(ierr)
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use SLEPC_ops, only: solve_EV_system_SLEPC
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        logical, intent(in) :: use_guess                                        ! whether to use a guess or not
        integer, intent(inout) :: n_sol_found                                   ! how many solutions saved
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(grid_eq,grid_X,X,use_guess,&
                    &n_sol_found)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
    
    ! Calculate  rho  from  user input.  This   is  done  here  in stead  of  in
    ! calc_flux_q  because  normalized quantities  are  necessary  to apply  the
    ! adiabatic law:
    !   rho p^gamma = rho_0 p_0^gamma.
    integer function calc_rho(eq,grid) result(ierr)
        use num_vars, only: eq_style, use_normalization
        use VMEC, only: gam
        use eq_vars, only: eq_type, rho_0
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'calc_rho'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(grid_type) :: grid                                                 ! grid
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp) :: expon                                                       ! exponent = 1/gam
        real(dp), parameter :: tol = 1.0E-10_dp                                 ! tolerance for negative pressure
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! allocate rho
        allocate(eq%rho(grid%grp_n_r))
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! set exp
                expon = 1.0_dp/gam
                
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    if (eq%pres_FD(kd,0).gt.0) then
                        eq%rho(kd) = eq%pres_FD(kd,0)**expon
                    else
                        eq%rho(kd) = eq%rho(kd-1)
                        if (eq%pres_FD(kd,0).lt.-tol) &
                        call writo('WARNING: pressure was negative ('//&
                            &trim(r2strt(eq%pres_FD(kd,0)))//') at point '&
                            &//trim(i2str(grid%i_min-1+kd))//'/'//&
                            &trim(i2str(grid%i_max))//&
                            &' so density is set to '//&
                            &trim(r2str(eq%rho(kd))),persistent=.true.)
                    end if
                end do
            case (2)                                                            ! HELENA
                eq%rho = 1.0_dp                                                 ! arbitrarily constant (normalized value)
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! normalize rho
        if (use_normalization) eq%rho = eq%rho/rho_0
    end function calc_rho
    
    ! calculate  ~PV_(k,m)^i  (pol.  flux)  or ~PV_(l,n)^i  (tor.  flux) at  all
    ! eq grp_n_r values
    ! (see [ADD REF] for details)
    subroutine calc_PV(eq,grid,met,X)
        use num_vars, only: use_pol_flux_F, use_normalization, mu_0
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type
        use utilities, only: c
        
        ! use input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:,:)                                 ! common factor |nabla psi|^2/(J^2*B^2)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp) :: mu_0_loc                                                    ! local version of mu_0
        integer :: c1                                                           ! value of c, to avoid compiler hang
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:)                                         ! h^alpha,psi
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^alpha,psi
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        ! lower metric factors
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        ! upper metric factors
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up local mu_0
        if (use_normalization) then
            mu_0_loc = 1._dp
        else
            mu_0_loc = mu_0
        end if
        
        ! set up common factor for PV_i
        allocate(com_fac(grid%n(1),grid%n(2),grid%grp_n_r))
        com_fac = h22/(g33*mu_0_loc)
        
        ! set up fac_n and fac_m
        allocate(fac_n(grid%grp_n_r),fac_m(grid%grp_n_r))
        if (use_pol_flux_F) then
            fac_n = eq%q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = eq%rot_t_FD(:,0)
        end if
        
        ! calculate PV_i (Hermitian)
        do m = 1,X%n_mod
            do k = m,X%n_mod
                ! calculate PV_0
                X%PV_0(:,:,:,c([k,m],.true.,X%n_mod)) = &
                    &com_fac*(X%DU_0(:,:,:,m) - X%extra1 - X%extra2 ) * &
                    &(conjg(X%DU_0(:,:,:,k)) - X%extra1 - X%extra2) - &
                    &mu_0/mu_0_loc * &
                    &(X%mu0sigma/J * (X%extra1 + X%extra2) - X%extra3)
                
                ! add (nq-k)*(nq-m)/(J^2 |nabla psi|^2) to PV_0
                do kd = 1,grid%grp_n_r
                    c1 = c([k,m],.true.,X%n_mod)
                    X%PV_0(:,:,kd,c1) = X%PV_0(:,:,kd,c1) + 1._dp/mu_0_loc * &
                        &(X%n(m)*fac_n(kd)-X%m(m)*fac_m(kd))*&
                        &(X%n(k)*fac_n(kd)-X%m(k)*fac_m(kd)) / &
                        &( J(:,:,kd)**2*h22(:,:,kd) )
                end do
                
                ! calculate PV_2
                X%PV_2(:,:,:,c([k,m],.true.,X%n_mod)) = &
                    &com_fac*X%DU_1(:,:,:,m)*conjg(X%DU_1(:,:,:,k))
                
                !write(*,*) '**k,m = ', k, m
                
                !write(*,*) 'RE PV_0',k,m
                !call print_HDF5('X','X',realpart(X%PV_0(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'IM PV_0',k,m
                !call print_HDF5('X','X',imagpart(X%PV_0(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'RE PV_2',k,m
                !call print_HDF5('X','X',realpart(X%PV_2(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'IM PV_2',k,m
                !call print_HDF5('X','X',imagpart(X%PV_2(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
            end do
        end do
        
        ! PV_1 is not Hermitian
        do m = 1,X%n_mod
            do k = 1,X%n_mod
                ! calculate PV_1
                X%PV_1(:,:,:,c([k,m],.false.,X%n_mod)) = &
                    &com_fac * X%DU_1(:,:,:,m) * &
                    &(conjg(X%DU_0(:,:,:,k)) - X%extra1 - X%extra2)
                
                !write(*,*) 'RE PV_1',k,m
                !call print_HDF5('X','X',realpart(X%PV_1(:,:,:,c([k,m],.false.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'IM PV_1',k,m
                !call print_HDF5('X','X',imagpart(X%PV_1(:,:,:,c([k,m],.false.,X%n_mod))))
                !read(*,*)
            end do
        end do
        
        ! deallocate variables
        nullify(J,g33,h22)
        deallocate(com_fac,fac_m,fac_n)
    end subroutine calc_PV
    
    ! calculate  ~KV_(k,m)^i  (pol.  flux)  or ~KV_(l,n)^i  (tor.  flux) at  all
    ! eq grp_n_r values
    ! (see [ADD REF] for details)
    subroutine calc_KV(eq,grid,met,X)
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type
        use utilities, only: c
        
        ! use input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:,:)                                 ! common factor |nabla psi|^2/(J^2*B^2)
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:)                                         ! h^alpha,psi
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^alpha,psi
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        ! lower metric factors
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        ! upper metric factors
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up common factor
        allocate(com_fac(grid%n(1),grid%n(2),grid%grp_n_r))
        com_fac = J**2*h22/g33
        
        ! for  Hermitian KV_0  and  KV_2, only  half  of the  terms  have to  be
        ! calculated
        do m = 1,X%n_mod
            do k = m,X%n_mod
                ! calculate KV_0
                X%KV_0(:,:,:,c([k,m],.true.,X%n_mod)) = com_fac * &
                    &X%U_0(:,:,:,m) * conjg(X%U_0(:,:,:,k)) + 1._dp/h22
                
                ! calculate KV_2
                X%KV_2(:,:,:,c([k,m],.true.,X%n_mod)) = com_fac * &
                    &X%U_1(:,:,:,m) * conjg(X%U_1(:,:,:,k))
            end do
        end do
        
        ! KV_1 is not Hermitian
        do m = 1,X%n_mod
            do k = 1,X%n_mod
                ! calculate KV_1
                X%KV_1(:,:,:,c([k,m],.false.,X%n_mod)) = com_fac * &
                    &X%U_1(:,:,:,m) * conjg(X%U_0(:,:,:,k))
            end do
        end do
        
        ! multiply by rho
        do kd = 1,grid%grp_n_r
            X%KV_0(:,:,kd,:) = X%KV_0(:,:,kd,:)*eq%rho(kd)
            X%KV_1(:,:,kd,:) = X%KV_1(:,:,kd,:)*eq%rho(kd)
            X%KV_2(:,:,kd,:) = X%KV_2(:,:,kd,:)*eq%rho(kd)
        end do
        
        !do m = 1,X%n_mod
            !do k = m,X%n_mod
                !write(*,*) '**k,m = ', k, m
                
                !write(*,*) 'RE KV_0',k,m
                !call print_HDF5('X','X',realpart(X%KV_0(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'IM KV_0',k,m
                !call print_HDF5('X','X',imagpart(X%KV_0(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'RE KV_2',k,m
                !call print_HDF5('X','X',realpart(X%KV_2(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'IM KV_2',k,m
                !call print_HDF5('X','X',imagpart(X%KV_2(:,:,:,c([k,m],.true.,X%n_mod))))
                !read(*,*)
            !end do
        !end do
        !do m = 1,X%n_mod
            !do k = 1,X%n_mod
                !write(*,*) 'RE KV_1',k,m
                !call print_HDF5('X','X',realpart(X%KV_1(:,:,:,c([k,m],.false.,X%n_mod))))
                !read(*,*)
                !write(*,*) 'IM KV_1',k,m
                !call print_HDF5('X','X',imagpart(X%KV_1(:,:,:,c([k,m],.false.,X%n_mod))))
                !read(*,*)
            !end do
        !end do
    end subroutine calc_KV
    
    ! calculate U_m^0, U_m^1 or U_n^0, U_n^1  at eq grp_n_r values of the normal
    ! coordinate, n_par values  of the parallel coordinate and  size_X values of
    ! the poloidal mode number,  or of the toroidal mode number,  as well as the
    ! poloidal derivatives
    ! Note: The normal derivatives have the  factor i/n or i/m included already,
    ! as opposed to [ADD REF]
    integer function calc_U(eq,grid,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, mu_0, use_normalization, eq_style
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'calc_U'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        real(dp) :: mu_0_loc                                                    ! local version of mu_0
        real(dp) :: n_frac                                                      ! nq-m (pol. flux) or n-iotam (tor. flux)
        real(dp), allocatable :: djq(:)                                         ! either q' (pol. flux) or -iota' (tor. flux)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp), allocatable :: mn(:)                                          ! either n*A_0 (pol. flux) or m (tor.flux)
        complex(dp), allocatable :: U_corr(:,:,:,:)                             ! correction to U for a certain (n,m)
        complex(dp), allocatable :: D3U_corr(:,:,:,:)                           ! D_theta U_corr
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle in flux coordinates
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! helper variables for the correction to U
        real(dp), pointer :: g_frac(:,:,:)                                      ! g_alpha,theta / g_theta,theta
        real(dp), pointer :: T_theta(:,:,:), D1T_theta(:,:,:), &
            &D3T_theta(:,:,:)                                                   ! Theta^theta and derivatives
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        real(dp), pointer :: D3J(:,:,:)                                         ! D_theta jac
        ! lower metric factors
        real(dp), pointer :: g13(:,:,:)                                         ! g_alpha,theta
        real(dp), pointer :: D3g13(:,:,:)                                       ! D_theta g_alpha,theta
        real(dp), pointer :: g23(:,:,:)                                         ! g_psi,theta
        real(dp), pointer :: D3g23(:,:,:)                                       ! D_theta g_psi,theta
        real(dp), pointer :: g33(:,:,:)                                         ! g_theta,theta
        real(dp), pointer :: D3g33(:,:,:)                                       ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), pointer :: h12(:,:,:)                                         ! h^alpha,psi
        real(dp), pointer :: D3h12(:,:,:)                                       ! D_theta h^alpha,psi
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        real(dp), pointer :: D1h22(:,:,:)                                       ! D_alpha h^psi,psi
        real(dp), pointer :: D3h22(:,:,:)                                       ! D_theta h^psi,psi
        real(dp), pointer :: h23(:,:,:)                                         ! h^psi,theta
        real(dp), pointer :: D1h23(:,:,:)                                       ! D_alpha h^psi,theta
        real(dp), pointer :: D3h23(:,:,:)                                       ! D_theta h^psi,theta
        
        ! initialize ierr
        ierr = 0
        
        ! set up local mu_0
        if (use_normalization) then
            mu_0_loc = 1._dp
        else
            mu_0_loc = mu_0
        end if
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_F) then
            ang_par_F => grid%theta_F
        else
            ang_par_F => grid%zeta_F
        end if
        
        ! set up djq, fac_n, fac_m and mn
        allocate(djq(grid%grp_n_r),mn(X%n_mod))
        allocate(fac_n(grid%grp_n_r),fac_m(grid%grp_n_r))
        if (use_pol_flux_F) then
            djq = eq%q_saf_FD(:,1)
            fac_n = eq%q_saf_FD(:,0)
            fac_m = 1.0_dp
            mn = X%n
        else
            djq = -eq%rot_t_FD(:,1)
            fac_n = 1.0_dp
            fac_m = eq%rot_t_FD(:,0)
            mn = X%m
        end if
        
        ! allocate helper variables
        allocate(U_corr(grid%n(1),grid%n(2),grid%grp_n_r,X%n_mod))
        allocate(D3U_corr(grid%n(1),grid%n(2),grid%grp_n_r,X%n_mod))
        allocate(g_frac(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(T_theta(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(D1T_theta(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(D3T_theta(grid%n(1),grid%n(2),grid%grp_n_r))
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        D3J => met%jac_FD(:,:,:,0,0,1)
        ! lower metric factors
        g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,0)
        D3g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,1)
        g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,0)
        D3g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,1)
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        D3g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,1)
        ! upper metric factors
        h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        D3h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,1)
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        D1h22 => met%h_FD(:,:,:,c([2,2],.true.),1,0,0)
        D3h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,1)
        h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,0)
        D1h23 => met%h_FD(:,:,:,c([2,3],.true.),1,0,0)
        D3h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,1)
        
        ! calculate X%U_0, X%U_1, X%DU_0 and X%DU_1
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call calc_U_VMEC
            case (2)                                                            ! HELENA
                ierr = calc_U_HEL()
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! deallocate
        deallocate(djq,mn)
        deallocate(fac_n,fac_m)
        deallocate(U_corr,D3U_corr)
        nullify(J,D3J)
        nullify(g13,D3g13,g23,D3g23,g33,D3g33)
        nullify(h12,D3h12,h22,D1h22,D3h22,h23,D1h23,D3h23)
        nullify(g_frac)
        nullify(T_theta,D1T_theta,D3T_theta)
        
        !do id = 1,X%n_mod
            !write(*,*) 'RE DU_',id
            !call print_HDF5('X','X',realpart(X%DU_0(:,:,:,id)))
            !read(*,*)
            !write(*,*) 'IM DU_0',id
            !call print_HDF5('X','X',imagpart(X%DU_0(:,:,:,id)))
            !read(*,*)
            !write(*,*) 'RE DU_1',id
            !call print_HDF5('X','X',realpart(X%DU_1(:,:,:,id)))
            !read(*,*)
            !write(*,*) 'IM DU_1',id
            !call print_HDF5('X','X',imagpart(X%DU_1(:,:,:,id)))
            !read(*,*)
        !end do 
    contains
        ! VMEC version
        subroutine calc_U_VMEC
            ! local variables
            ! extra helper variables for the correction to U
            real(dp), allocatable :: D3g_frac(:,:,:)                            ! D_theta g_frac
            real(dp), allocatable :: D13T_theta(:,:,:), D33T_theta(:,:,:)       ! Theta^theta derivatives
            ! extra upper metric factors
            real(dp), allocatable :: D13h22(:,:,:)                              ! D^2_alpha,theta h^psi,psi
            real(dp), allocatable :: D33h22(:,:,:)                              ! D^2_theta,theta h^psi,psi
            real(dp), allocatable :: D13h23(:,:,:)                              ! D^2_alpha,theta h^psi,theta
            real(dp), allocatable :: D33h23(:,:,:)                              ! D^2_theta,theta h^psi,theta
            
            ! allocate extra helper variables
            allocate(D3g_frac(grid%n(1),grid%n(2),grid%grp_n_r))
            allocate(D13T_theta(grid%n(1),grid%n(2),grid%grp_n_r))
            allocate(D33T_theta(grid%n(1),grid%n(2),grid%grp_n_r))
            
            ! set up extra submatrices
            ! upper metric factors
            allocate(D13h22(grid%n(1),grid%n(2),grid%grp_n_r))
            D13h22 = met%h_FD(:,:,:,c([2,2],.true.),1,0,1)
            allocate(D33h22(grid%n(1),grid%n(2),grid%grp_n_r))
            D33h22 = met%h_FD(:,:,:,c([2,2],.true.),0,0,2)
            allocate(D13h23(grid%n(1),grid%n(2),grid%grp_n_r))
            D13h23 = met%h_FD(:,:,:,c([2,3],.true.),1,0,1)
            allocate(D33h23(grid%n(1),grid%n(2),grid%grp_n_r))
            D33h23 = met%h_FD(:,:,:,c([2,3],.true.),0,0,2)
            
            ! loop over the M elements of U_X and DU
            do jd = 1,X%n_mod
                ! set up helper variables
                g_frac = g13/g33
                D3g_frac = D3g13/g33 - g13*D3g33/(g33**2)
                T_theta = h23/h22
                D1T_theta = D1h23/h22 - h23*D1h22/(h22**2)
                D3T_theta = D3h23/h22 - h23*D3h22/(h22**2)
                D13T_theta = D13h23/h22 - (D3h23*D1h22+D1h23*D3h22)/(h22**2) &
                    &- h23*D13h22/(h22**2) + 2*h23*D3h22*D1h22/(h22**3)
                D33T_theta = D33h23/h22 - 2*D3h23*D3h22/(h22**2) &
                    &- h23*D33h22/(h22**2) + 2*h23*D3h22**2/(h22**3)
                
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    ! set up correction to U
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*&
                        &(g_frac(:,:,kd)*(D3T_theta(:,:,kd)+&
                        &iu*n_frac*T_theta(:,:,kd))+D1T_theta(:,:,kd))
                    ! set up D_theta U correction
                    D3U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*&
                        &(D3g_frac(:,:,kd)*&
                        &(D3T_theta(:,:,kd)+iu*n_frac*T_theta(:,:,kd))+&
                        &g_frac(:,:,kd)*&
                        &(D33T_theta(:,:,kd)+iu*n_frac*D3T_theta(:,:,kd))+&
                        &D13T_theta(:,:,kd))
                    ! calculate X%U_0 and X%DU_0
                    X%U_0(:,:,kd,jd) = &
                        &-(h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd))&
                        &+ iu/(mn(jd)*g33(:,:,kd)) * (g13(:,:,kd)*djq(kd) + &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*mu_0_loc + iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) - &
                        &g23(:,:,kd) )) + U_corr(:,:,kd,jd)
                    X%DU_0(:,:,kd,jd) = -(D3h12(:,:,kd)/h22(:,:,kd) - &
                        &D3h22(:,:,kd)*h12(:,:,kd)/h22(:,:,kd)**2 + djq(kd)) - &
                        &iu*D3g33(:,:,kd)/(mn(jd)*g33(:,:,kd)**2) * &
                        &(g13(:,:,kd)*djq(kd)+ &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*mu_0_loc + &
                        &iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) &
                        &- g23(:,:,kd) )) + &
                        &iu/(mn(jd)*g33(:,:,kd)) * (D3g13(:,:,kd)*djq(kd) + &
                        &2*D3J(:,:,kd)*J(:,:,kd)*eq%pres_FD(kd,1)*mu_0_loc + &
                        &iu*n_frac &
                        &* ( D3g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) + &
                        &g13(:,:,kd)*djq(kd) - D3g23(:,:,kd) )) + &
                        &iu*n_frac*X%U_0(:,:,kd,jd) &
                        &+ D3U_corr(:,:,kd,jd)
                    ! calculate X%U_1 and X%DU_1
                    X%U_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,:,kd)/g33(:,:,kd))
                    X%DU_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(n_frac/mn(jd) * (D3g13(:,:,kd)/g33(:,:,kd) - &
                        &D3g33(:,:,kd)*g13(:,:,kd)/g33(:,:,kd)**2) + &
                        &iu*n_frac*X%U_1(:,:,kd,jd) )
                end do
            end do
            
            ! deallocate
            deallocate(D3g_frac)
            deallocate(D13T_theta,D33T_theta)
        end subroutine calc_U_VMEC
        
        ! HELENA version
        integer function calc_U_HEL() result(ierr)
            use utilities, only: calc_deriv
            
            character(*), parameter :: rout_name = 'calc_U_HEL'
            
            ! local variablas
            real(dp), allocatable :: D3_var(:,:,:)                              ! derivative of variable
            
            ! allocate extra helper variables
            allocate(D3_var(grid%n(1),grid%n(2),2))
            
            ! initialize ierr
            ierr = 0
            
            ! loop over the M elements of U_X and DU
            do jd = 1,X%n_mod
                ! set up helper variables
                g_frac = g13/g33
                T_theta = h23/h22
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(T_theta(:,id,kd),D3T_theta(:,id,kd),&
                            &grid%theta_F(:,id,kd),1,2)                         ! higher precision because other derivative will be taken later
                    end do
                end do
                
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    ! set up correction to U
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*(g_frac(:,:,kd)*&
                        &(D3T_theta(:,:,kd)+iu*n_frac*T_theta(:,:,kd)))
                    ! calculate X%U_0 and X%DU_0
                    X%U_0(:,:,kd,jd) = &
                        &-(h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd))&
                        &+ iu/(mn(jd)*g33(:,:,kd)) * (g13(:,:,kd)*djq(kd) + &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*mu_0_loc + iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) - &
                        &g23(:,:,kd) )) + U_corr(:,:,kd,jd)
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(realpart(X%U_0(:,id,kd,jd)),&
                            &D3_var(:,id,1),grid%theta_F(:,id,kd),1,1)
                        ierr = calc_deriv(imagpart(X%U_0(:,id,kd,jd)),&
                            &D3_var(:,id,2),grid%theta_F(:,id,kd),1,1)
                        end do
                        CHCKERR('')
                    X%DU_0(:,:,kd,jd) = D3_var(:,:,1) + iu*D3_var(:,:,2) + &
                        &iu*n_frac*X%U_0(:,:,kd,jd)
                    ! calculate X%U_1 and X%DU_1
                    X%U_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,:,kd)/g33(:,:,kd))
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(realpart(X%U_1(:,id,kd,jd)),&
                            &D3_var(:,id,1),grid%theta_F(:,id,kd),1,1)
                        ierr = calc_deriv(imagpart(X%U_1(:,id,kd,jd)),&
                            &D3_var(:,id,2),grid%theta_F(:,id,kd),1,1)
                    end do
                    CHCKERR('')
                    X%DU_1(:,:,kd,jd) = D3_var(:,:,1) + iu*D3_var(:,:,2) + &
                        &iu*n_frac*X%U_1(:,:,kd,jd)
                end do
            end do
            
            ! deallocate
            deallocate(D3_var)
        end function calc_U_HEL
    end function calc_U
    
    ! Calculate mu0sigma, extra1, extra2 and extra3:
    !   extra1 = S*J
    !   extra2 = mu0sigma*J*B^2/h^psi,psi
    !   extra3 = 2*p'*kn
    ! with
    !   shear S = - d Theta^alpha/d_theta 1 / J
    !   parallel current mu0sigma = B . nabla x B / B^2
    !   normal curvature kn = nabla psi / h^psi,psi * nabla (mu0 p + B^2/2)
    integer function calc_extra(eq,grid,met,X) result(ierr)
        use num_vars, only: mu_0
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'calc_extra'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: kd                                                           ! counter
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        real(dp), pointer :: D1J(:,:,:)                                         ! D_alpha jac
        real(dp), pointer :: D2J(:,:,:)                                         ! D_psi jac
        real(dp), pointer :: D3J(:,:,:)                                         ! D_theta jac
        ! lower metric factors
        real(dp), pointer :: g13(:,:,:)                                         ! g_alpha,theta
        real(dp), pointer :: D2g13(:,:,:)                                       ! D_psi g_alpha,theta
        real(dp), pointer :: D3g13(:,:,:)                                       ! D_theta g_alpha,theta
        real(dp), pointer :: g23(:,:,:)                                         ! g_psi,theta
        real(dp), pointer :: D1g23(:,:,:)                                       ! D_alpha g_psi,theta
        real(dp), pointer :: D3g23(:,:,:)                                       ! D_theta g_psi,theta
        real(dp), pointer :: g33(:,:,:)                                         ! g_theta,theta
        real(dp), pointer :: D1g33(:,:,:)                                       ! D_alpha g_theta,theta
        real(dp), pointer :: D2g33(:,:,:)                                       ! D_psi g_theta,theta
        real(dp), pointer :: D3g33(:,:,:)                                       ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), pointer :: h12(:,:,:)                                         ! h^alpha,psi
        real(dp), pointer :: D3h12(:,:,:)                                       ! D_theta h^alpha,psi
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        real(dp), pointer :: D3h22(:,:,:)                                       ! D_theta h^psi,psi
        real(dp), pointer :: h23(:,:,:)                                         ! h^psi,theta
        
        ! initialize ierr
        ierr = 0
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        D1J => met%jac_FD(:,:,:,1,0,0)
        D2J => met%jac_FD(:,:,:,0,1,0)
        D3J => met%jac_FD(:,:,:,0,0,1)
        ! lower metric factors
        g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,0)
        D2g13 => met%g_FD(:,:,:,c([1,3],.true.),0,1,0)
        D3g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,1)
        g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,0)
        D1g23 => met%g_FD(:,:,:,c([2,3],.true.),1,0,0)
        D3g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,1)
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        D1g33 => met%g_FD(:,:,:,c([3,3],.true.),1,0,0)
        D2g33 => met%g_FD(:,:,:,c([3,3],.true.),0,1,0)
        D3g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,1)
        ! upper metric factors
        h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        D3h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,1)
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        D3h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,1)
        h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,0)
        
        ! calculate mu0sigma
        X%mu0sigma = 1._dp/(J*g33) * (g13*(D2g33-D3g23) + g23*(D3g13-D1g33) &
            &+ g33*(D1g23-D2g13))
        
        ! calculate extra1
        X%extra1 = -D3h12/h22 + D3h22*h12/h22**2
        
        ! calculate extra2
        X%extra2 = g33/h22 * X%mu0sigma / J
        
        ! calculate extra3
        do kd = 1,grid%grp_n_r
            X%extra3(:,:,kd) = eq%pres_FD(kd,1) * &
                &( 2*J(:,:,kd)**2*mu_0*eq%pres_FD(kd,1)/g33(:,:,kd) + &
                &1._dp/h22(:,:,kd) * ( &
                &h12(:,:,kd) * ( D1g33(:,:,kd)/g33(:,:,kd) - &
                &2*D1J(:,:,kd)/J(:,:,kd) ) + &
                &h22(:,:,kd) * ( D2g33(:,:,kd)/g33(:,:,kd) - &
                &2*D2J(:,:,kd)/J(:,:,kd) ) + &
                &h23(:,:,kd) * ( D3g33(:,:,kd)/g33(:,:,kd) - &
                &2*D3J(:,:,kd)/J(:,:,kd) ) ) )
        end do
        
        !write(*,*) 'mu0sigma'
        !call print_GP_2D('mu0sigma(:,1,:)','',transpose(X%mu0sigma(:,1,:)))
        !write(*,*) 'extra1'
        !call print_GP_2D('extra1(:,1,:)','',transpose(X%extra1(:,1,:)))
        !write(*,*) 'extra2'
        !call print_GP_2D('extra2(:,1,:)','',transpose(X%extra2(:,1,:)))
        !write(*,*) 'extra3'
        !call print_GP_2D('extra3(:,1,:)','',transpose(X%extra3(:,1,:)))
        
        ! calculate rho from input (best done with normalized pressure already)
        ierr = calc_rho(eq,grid)
        CHCKERR('')
        
        ! deallocate local variables
        nullify(J, D1J, D2J, D3J)
        nullify(g13, D2g13, D3g13)
        nullify(g23, D1g23, D3g23)
        nullify(g33, D1g33, D2g33, D3g33)
        nullify(h12, D3h12)
        nullify(h22, D3h22)
        nullify(h23)
        
        !write(*,*) 'extra1'
        !call print_HDF5('X','X',(X%extra1))
        !read(*,*)
        !write(*,*) 'extra2'
        !call print_HDF5('X','X',(X%extra2))
        !read(*,*)
        !write(*,*) 'extra3'
        !call print_HDF5('X','X',(X%extra3))
        !read(*,*)
    end function calc_extra
    
    ! calculates  magnetic  integral  <V e^[i(k-m)ang_par_F]>,  defined  as  the
    ! matrix  
    !   <V e^[i(k-m)ang_par_F]> = [ oint J V(k,m) e^i(k-m)ang_par_F dang_par_F ]
    ! or
    !   <V e^[i(n-l)ang_par_F]> = [ oint J V(l,n) e^i(n-l)ang_par_F dang_par_F ]
    integer function calc_V_int(grid,met,exp_ang,n_mod,V,V_int) result(ierr)
        use num_vars, only: use_pol_flux_F
        use grid_vars, only: grid_type
        use metric_vars, only: metric_type
        use utilities, only: calc_mult, c, con, is_sym
        
        character(*), parameter :: rout_name = 'calc_V_int'
        
        ! input / output
        type(grid_type) :: grid                                                 ! grid
        type(metric_type) :: met                                                ! metric variables
        complex(dp), intent(in) :: exp_ang(:,:,:,:)                             ! exponential of Flux parallel angle
        integer, intent(in) :: n_mod                                            ! number of 
        complex(dp), intent(in) :: V(:,:,:,:)                                   ! input V(n_par,n_geo,n_r,size_X^2)
        complex(dp), intent(inout) :: V_int(:,:,:)                              ! output <V e^i(k-m)ang_par_F> integrated in parallel Flux coord.
        
        ! local variables
        integer :: k, m, id, jd, kd                                             ! counters
        integer :: nn_mod                                                       ! number of indices for V and V_int
        integer :: k_min                                                        ! minimum k
        logical :: sym                                                          ! whether V and V_int are symmetric
        complex(dp), allocatable :: V_J_e(:,:,:,:)                              ! V*J*exp_ang
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle
        integer :: dims(3)                                                      ! real dimensions
        
        ! initialize ierr
        ierr = 0
        
        ! set nn_mod
        nn_mod = size(V,4)
        
        ! tests
        if (size(V_int,1).ne.nn_mod) then
            ierr = 1
            err_msg = 'V and V_int need to have the same storage convention'
            CHCKERR(err_msg)
        end if
        
        ! set up dims
        dims = [grid%n(1),grid%n(2),grid%grp_n_r]
        
        ! set up V_J_e
        allocate(V_J_e(dims(1),dims(2),dims(3),nn_mod))
        
        ! set up ang_par_F
        if (use_pol_flux_F) then
            ang_par_F => grid%theta_F
        else
            ang_par_F => grid%zeta_F
        end if
        
        ! determine whether matrices are symmetric or not
        ierr = is_sym(n_mod,nn_mod,sym)
        CHCKERR('')
        
        ! set up k_min
        k_min = 1
        
        ! multiply V by Jacobian and exponential
        do m = 1,n_mod
            if (sym) k_min = m
            do k = k_min,n_mod
                V_J_e(:,:,:,c([k,m],sym,n_mod)) = met%jac_FD(:,:,:,0,0,0) * &
                    &con(exp_ang(:,:,:,c([k,m],.true.,n_mod)),&
                    &[k,m],.true.,dims) * &
                    &con(V(:,:,:,c([k,m],sym,n_mod)),[k,m],sym,dims)
                !write(*,*) 'RE V_J', k,m
                !call print_HDF5('X','X',realpart(V_J_e(:,:,:,c([k,m],sym,n_mod))))
                !read(*,*)
                !write(*,*) 'IM V_J', k,m
                !call print_HDF5('X','X',imagpart(V_J_e(:,:,:,c([k,m],sym,n_mod))))
                !read(*,*)
            end do
        end do
        
        ! integrate term  over ang_par_F for  all equilibrium grid  points using
        ! the recursive formula int_1^n f(x) dx
        !   = int_1^(n-1) f(x) dx + (f(n)+f(n-1))*(x(n)-x(n-1))/2
        V_int = 0.0_dp
        ! loop over all geodesic points on this process
        do kd = 1,grid%grp_n_r
            ! loop over all normal points on this process
            do jd = 1,grid%n(2)
                ! parallel integration loop
                do id = 2,grid%n(1)
                    V_int(:,jd,kd) = V_int(:,jd,kd) + &
                        &(V_J_e(id,jd,kd,:)+V_J_e(id-1,jd,kd,:))/2 * &
                        &(ang_par_F(id,jd,kd)-ang_par_F(id-1,jd,kd))
                end do
            end do
        end do
        
        !do m = 1,n_mod
            !do k = 1,n_mod
                !write(*,*) 'INTEGRATED', k,m
                !call print_HDF5('X','X',reshape([realpart(transpose(V_int(c([k,m],sym,n_mod),:,:))),&
                    !&imagpart(transpose(con(V_int(c([k,m],sym,n_mod),:,:),[k,m],sym,[grid%n(2),grid%grp_n_r])))],&
                    !&[grid%grp_n_r,grid%n(2),2]))
                !write(*,*) 'last point:', con(V_int(c([k,m],sym,n_mod),:,grid%grp_n_r-1:grid%grp_n_r),&
                    !&[k,m],sym,[grid%n(2),2])
                !read(*,*)
            !end do
        !end do
        
        ! deallocate local variables
        deallocate(V_J_e)
        nullify(ang_par_F)
    end function calc_V_int
    
    ! Plots an Eigenfunction
    ! [MPI] Collective call
    integer function plot_X_vec(grid_eq,eq,grid_X,X,X_id,job_id,par_range_F) &
        &result(ierr)
        use num_vars, only: output_style
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        real(dp), intent(in), optional :: par_range_F(2)                        ! parallel range [2pi]
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! local variables
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        
        ! initialize ierr
        ierr = 0
        
        ! if the  Eigenvalue is negative,  the Eigenfunction explodes,  so limit
        ! n_t(2) to 1. If it is positive, the Eigenfunction oscilates, so choose
        ! n_t(2) = 8 for 2 whole periods
        if (realpart(X%val(X_id)).lt.0) then                                    ! exploding, unstable
            n_t(1) = 1                                                          ! 1 point per quarter period
            n_t(2) = 1                                                          ! 1 quarter period
        else                                                                    ! oscillating, stable
            n_t(1) = 2                                                          ! 2 points per quarter period
            n_t(2) = 4                                                          ! 4 quarter periods
        end if
        
        ! delegate to child routines
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                ierr = plot_X_vec_GP(eq,grid_X,X,X_id,job_id,n_t,par_range_F)
                CHCKERR('')
            case(2)                                                             ! HDF5 output
                ierr = plot_X_vec_HDF5(grid_eq,eq,grid_X,X,X_id,job_id,n_t,&
                    &.false.)
                CHCKERR('')
            case default
                err_msg = 'No style associated with '//trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function plot_X_vec
    
    ! GNUPlot version of plot_X_vec
    ! plots an  eigenfunction in either 1D  (normal coord. in perturb.  grid) or
    ! 2D  (normal  coord. &  par  coord.),  depending  on whether  the  variable
    ! par_range_F is provided The individual modes are  plotted at t = 0 and t =
    ! (pi/2)/sqrt(X_val) and  an animated  plot is  produced as  well in  a .gif
    ! output file
    integer function plot_X_vec_GP(eq,grid_X,X,X_id,job_id,n_t,par_range_F) &
        &result(ierr)
        use num_vars, only: grp_rank, max_n_plots, grp_n_procs, &
            &use_pol_flux_F, use_pol_flux_E
        use grid_vars, only: grid_type
        use eq_vars, only: eq_type
        use X_vars, only: min_r_X, max_r_X, X_type
        use output_ops, only: draw_GP, draw_GP_animated, merge_GP
        use MPI_ops, only: get_ghost_arr, wait_MPI
        use eq_vars, only: eq_type
        use utilities, only: con2dis, interp_fun
        
        character(*), parameter :: rout_name = 'plot_X_vec_GP'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        integer, intent(in) :: n_t(2)                                           ! nr. of time steps in quarter period, nr. of quarter periods
        real(dp), intent(in), optional :: par_range_F(2)                        ! parallel range [2pi]
        
        ! local variables
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        complex(dp), allocatable :: X_vec_extended(:,:)                         ! MPI Eigenvector extended with assymetric ghost region
        real(dp), allocatable :: X_plot(:,:,:,:)                                ! x values of plot, parallel version
        real(dp), allocatable :: Y_plot(:,:,:,:)                                ! y values of plot, parallel version
        real(dp), allocatable :: Z_plot(:,:,:,:)                                ! z values of plot, parallel version
        real(dp), allocatable :: z_magn_plot(:,:,:)                             ! z values of plot of magnitude, parallel version
        integer :: id, jd, kd, ld                                               ! counter
        integer :: n_par_F = 50                                                 ! how many parallel points
        real(dp) :: kd_loc_eq                                                   ! kd_loc in equilibrium grid
        real(dp) :: kd_loc_i                                                    ! unrounded index corresponding to kd_loc in tables
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: file_name                                  ! file name for plots
        real(dp), allocatable :: ang_par_F_loc(:)                               ! local version of ang_par_F
        complex(dp), allocatable :: X_vec_interp(:)                             ! X_vec at an interpolated normal position
        complex(dp), allocatable :: first_X_vec_interp(:)                       ! first X_vec_interp
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp) :: fac_n_interp, fac_m_interp                                  ! fac_n and fac_m at interpolated normal position
        real(dp) :: first_fac_n_interp, first_fac_m_interp                      ! first fac_n_interp and fac_m_interp
        character(len=max_str_ln), allocatable :: file_names(:)                 ! names of file of plots of different ranks
        real(dp), pointer :: flux_X(:), flux_eq(:)                              ! either pol. or tor. flux
        real(dp), allocatable :: r_plot(:)                                      ! local normal values at which to interpolate for the plot
        integer :: n_norm                                                       ! max nr of normal points to plot
        integer :: grp_n_norm                                                   ! nr. of normal points in plot for this rank
        complex(dp) :: fac_i_interp(1)                                          ! complex copy of fac_n_interp or fac_m_interp
        
        ! initialize ierr
        ierr = 0
        
        ! 1. initialize some quantities
        
        ! calculate  omega = sqrt(X_val)  and make  sure to select  the decaying
        ! solution
        omega = sqrt(X%val(X_id))
        if (imagpart(omega).gt.0) omega = - omega                               ! exploding solution, not the decaying one
        
        ! set up fac_n and fac_m
        allocate(fac_n(size(eq%q_saf_FD,1)),fac_m(size(eq%rot_t_FD,1)))
        if (use_pol_flux_F) then
            fac_n = eq%q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = eq%rot_t_FD(:,0)
        end if
        
        ! set up the angular variable ang_par_F_loc, depending on whether it has
        ! been provided as input or not
        if (present(par_range_F)) then
            allocate(ang_par_F_loc(n_par_F))
            ang_par_F_loc = &
                &[(par_range_F(1) + (jd-1.0_dp)/(n_par_F-1.0_dp) * &
                &(par_range_F(2)-par_range_F(1)),jd=1,n_par_F)]
        else
            n_par_F = 1
            allocate(ang_par_F_loc(n_par_F))
            ang_par_F_loc = 0.0_dp
        end if
        
        ! 2. calculate the x, y and z parts of the plot in parallel
        
        ! calculate starting and ending normal index of array
        call calc_r_plot(r_plot,n_norm,grp_n_norm,.true.)
        
        ! allocate variables
        if (grp_rank.lt.grp_n_procs-1) then
            allocate(X_vec_extended(size(X%vec(:,:,X_id),1),&
                &size(X%vec(:,:,X_id),2)+1))
        else
            allocate(X_vec_extended(size(X%vec(:,:,X_id),1),&
                &size(X%vec(:,:,X_id),2)))
        end if
        allocate(X_plot(grp_n_norm,n_par_F,X%n_mod,product(n_t)))
        allocate(Y_plot(grp_n_norm,n_par_F,X%n_mod,product(n_t)))
        allocate(Z_plot(grp_n_norm,n_par_F,X%n_mod,product(n_t)))
        allocate(z_magn_plot(grp_n_norm,n_par_F,product(n_t)))
        allocate(X_vec_interp(X%n_mod))
        allocate(first_X_vec_interp(X%n_mod))
        
        ! set up extended X_vec
        X_vec_extended(:,1:size(X%vec(:,:,X_id),2)) = X%vec(:,:,X_id)
        ierr = get_ghost_arr(X_vec_extended,1)
        CHCKERR('')
        
        ! set up y-axis in parallel
        if (present(par_range_F)) then
            do id = 1,n_par_F
                Y_plot(:,id,:,:) = ang_par_F_loc(id)
            end do
        else
            Y_plot = 0.0_dp
        end if
        
        ! set up flux_X and flux_eq
        if (use_pol_flux_F) then
            flux_X => eq%flux_p_FD(:,0)
        else
            flux_X => eq%flux_t_FD(:,0)
        end if
        if (use_pol_flux_E) then
            flux_eq => eq%flux_p_FD(:,0)
        else
            flux_eq => eq%flux_t_FD(:,0)
        end if
        
        ! set up x-axis and z-axis values in parallel
        z_magn_plot = 0.0_dp
        ! loop over all time steps
        do id = 1,product(n_t)
            ! loop over all normal points
            do kd = 1,grp_n_norm
                ! set up X_plot
                X_plot(kd,:,:,:) = r_plot(kd)
                
                ! set up the interpolated values of X_vec, fac_n and fac_m
                if (kd.lt.grp_n_norm .and. grp_rank+1.lt.grp_n_procs &
                    &.or. grp_rank+1.eq.grp_n_procs) then                       ! last point is treated differently
                    ! set up interp. X_vec (tabulated in perturbation grid)
                    call con2dis(r_plot(kd),kd_loc_i,grid_X%grp_r_F)            ! perturbation values tabulated at grp_r_F for this group
                    X_vec_interp = X_vec_extended(:,floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(X_vec_extended(:,ceiling(kd_loc_i))-&
                        &X_vec_extended(:,floor(kd_loc_i)))
                    
                    ! set up  interp. fac_n and fac_m  (tabulated in equilibrium
                    ! grid)
                    ierr = interp_fun(kd_loc_eq,flux_eq/eq%max_flux_F,&
                        &r_plot(kd),flux_X/X%max_flux_F)
                    CHCKERR('')
                    call con2dis(kd_loc_eq,kd_loc_i,flux_eq/eq%max_flux_F)      ! equilibrium values tabulated at flux_eq for this group, normalized with max_flux_eq_F
                    fac_n_interp = fac_n(floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(fac_n(ceiling(kd_loc_i))-fac_n(floor(kd_loc_i)))
                    fac_m_interp = fac_m(floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(fac_m(ceiling(kd_loc_i))-fac_m(floor(kd_loc_i)))
                end if
                
                ! save first interpolated X_vec, fac_n and fac_m
                if (kd.eq.1) then
                    first_X_vec_interp = X_vec_interp
                    first_fac_n_interp = fac_n_interp
                    first_fac_m_interp = fac_m_interp
                end if
                
                ! exchange X_vec_interp, fac_n and fac_m at last points
                if (kd.eq.grp_n_norm) then
                    ! X_vec_interp
                    call exch_ghost_values(first_X_vec_interp)
                    if (grp_rank+1.lt.grp_n_procs) &
                        &X_vec_interp = first_X_vec_interp
                    
                    ! fac_n_interp
                    fac_i_interp = first_fac_n_interp
                    call exch_ghost_values(fac_i_interp)
                    if (grp_rank+1.lt.grp_n_procs) &
                        &fac_n_interp = realpart(fac_i_interp(1))
                    
                    ! fac_m_interp
                    fac_i_interp = first_fac_m_interp
                    call exch_ghost_values(fac_i_interp)
                    if (grp_rank+1.lt.grp_n_procs) &
                        &fac_m_interp = realpart(fac_i_interp(1))
                end if
                
                ! loop over all modes and all parallel points to set up Z_plot
                do ld = 1,X%n_mod
                    do jd = 1,n_par_F
                        Z_plot(kd,jd,ld,id) = &
                            &realpart(X_vec_interp(ld) * exp(iu * &
                            &(omega/abs(omega)*0.5*pi*(id-1.0_dp)/&
                            &max(n_t(1)-1,1) + X%n(ld)*fac_n_interp-&
                            &          X%m(ld)*fac_m_interp)&
                            &    * ang_par_F_loc(jd)))
                    end do
                    z_magn_plot(kd,:,id) = z_magn_plot(kd,:,id) + &
                        &Z_plot(kd,:,ld,id)
                end do
            end do
        end do
        
        ! 3. for  every mode, write plot  output data in data files  and plot by
        ! group master
        if (X%n_mod.gt.1) then
            ! plot the individual modes in time
            do ld = 1,min(X%n_mod,max_n_plots)
                ! set plot_title and file name and write data
                plot_title = 'job '//trim(i2str(job_id))//' - EV '//&
                    &trim(i2str(X_id))//' - mode '//trim(i2str(ld))//' of '//&
                    &trim(i2str(X%n_mod))
                file_name = 'mode_'//trim(i2str(X_id))//'_'//&
                    &trim(i2str(ld))//'.dat'
                if (present(par_range_F)) then
                    call print_GP_3D(trim(plot_title),trim(file_name)//'_'//&
                        &trim(i2str(grp_rank)),Z_plot(:,:,ld,:),&
                        &y=Y_plot(:,:,ld,:),x=X_plot(:,:,ld,:),draw=.false.)
                else
                    call print_GP_2D(trim(plot_title),trim(file_name)//'_'//&
                        &trim(i2str(grp_rank)),Z_plot(:,1,ld,:),&
                        &x=X_plot(:,1,ld,:),draw=.false.)
                end if
                    
                ! wait for all processes
                ierr = wait_MPI()
                CHCKERR('')
                
                ! plot by group master
                if (grp_rank.eq.0) then
                    ! set up file names in array
                    allocate(file_names(grp_n_procs))
                    do id = 1,grp_n_procs
                        file_names(id) = trim(file_name)//'_'//trim(i2str(id-1))
                    end do
                    
                    ! merge files
                    call merge_GP(file_names,file_name,delete=.true.)
                    
                    ! draw animation
                    call draw_GP_animated(trim(plot_title),file_name,&
                        &product(n_t),.false.)
                    
                    ! deallocate
                    deallocate(file_names)
                end if
            end do
            if (X%n_mod.eq.max_n_plots+1) then
                call writo('Skipping plot for mode '//&
                    &trim(i2str(max_n_plots+1)))
            else if (X%n_mod.gt.max_n_plots+1) then
                call writo('Skipping plots for modes '//&
                    &trim(i2str(max_n_plots+1))//'..'//trim(i2str(X%n_mod)))
            end if
        end if
        
        ! 4. do the same for the global mode
        ! set plot_title and file name and write data
        plot_title = 'job '//trim(i2str(job_id))//' - EV '//&
            &trim(i2str(X_id))//' - global mode'
        file_name = 'global_mode_'//trim(i2str(X_id))//'.dat'
        if (present(par_range_F)) then
            call print_GP_3D(trim(plot_title),trim(file_name)//&
                &'_'//trim(i2str(grp_rank)),&
                &z_magn_plot(:,:,:),y=Y_plot(:,:,1,:),&
                &x=X_plot(:,:,1,:),draw=.false.)
        else
            call print_GP_2D(trim(plot_title),trim(file_name)//&
                &'_'//trim(i2str(grp_rank)),&
                &z_magn_plot(:,1,:),x=X_plot(:,1,1,:),&
                &draw=.false.)
        end if
        
        ! plot by group master
        if (grp_rank.eq.0) then
            ! set up file names in file
            allocate(file_names(grp_n_procs))
            do id = 1,grp_n_procs
                file_names(id) = trim(file_name)//'_'//trim(i2str(id-1))
            end do
            
            ! merge files
            call merge_GP(file_names,file_name,delete=.true.)
            
            ! draw animation
            call draw_GP_animated(trim(plot_title),file_name,product(n_t),&
                &.false.)
            
            ! deallocate
            deallocate(file_names)
        end if
    contains
        ! calculates the array  of normal points in the  perturbation grid where
        ! to interpolate  to produce the plot.  This assumes that the  vector to
        ! be  plot  has  an  assymetric  ghost  region  of  size  one,  so  that
        ! the  interpolated  plot points  can  be  calculated in  every  process
        ! independently
        ! Additionally,  there can be  an extra  assymetric ghost region  in the
        ! plot grid so that there is overlap (indicated by extra_ghost)
        subroutine calc_r_plot(r_plot,n_norm,grp_n_norm,extra_ghost)
            ! input / output
            real(dp), intent(inout), allocatable :: r_plot(:)                   ! local normal values at which to interpolate for the plot
            integer, intent(inout) :: grp_n_norm                                ! nr. of normal points for this rank
            integer, intent(inout) :: n_norm                                    ! total nr. of normal points
            logical, intent(in), optional :: extra_ghost                        ! .true. if extra ghost region requested (e.g. for GNUPlot)
            
            ! local variables
            integer :: start_i_norm, end_i_norm                                 ! starting and ending index of normal points for this rank
            integer :: max_n_norm = 50                                          ! nr. of normal points
            real(dp) :: inc_norm                                                ! normal increment for plots
            integer :: kd                                                       ! counter
            
            ! set  up how  many  normal  points for  this  rank  n_norm and  the
            ! increment between normal points inc_norm
            ! (taking all the normal points could be way too many)
            n_norm = min(grid_X%n(3),max_n_norm)
            inc_norm = (max_r_X-min_r_X)/(n_norm-1)
            
            ! calculate the  starting and ending index,  first determining which
            ! is the range where normal interpolation can take place
            ! from:
            !   gpr_r_X(1) .le. (min_r_X+(k-1)*inc_norm) .le. max grp_r_F
            ! with k an integer
            start_i_norm = ceiling((grid_X%grp_r_F(1)-min_r_X)/inc_norm + 1)
            end_i_norm = floor((grid_X%grp_r_F(size(grid_X%grp_r_F))-min_r_X)/&
                &inc_norm + 1)                                                  ! grp_r_F of X grid contains ghost region
            
            ! increment end_i_norm with one if not  last rank, so that there can
            ! be overlap in the plots
            if (present(extra_ghost)) then
                if (extra_ghost) then
                    if (grp_rank.lt.grp_n_procs-1) end_i_norm = end_i_norm + 1
                end if
            end if
            
            ! set grp_n_norm
            grp_n_norm = end_i_norm-start_i_norm+1
            
            ! allocate r_plot with the correct size
            allocate(r_plot(grp_n_norm)); 
            
            ! set r_plot
            r_plot = [(min_r_X + (start_i_norm+kd-2)*inc_norm,kd=1,grp_n_norm)]
            if (r_plot(grp_n_norm).gt.1.0_dp) r_plot(grp_n_norm) = 1.0_dp
        end subroutine calc_r_plot
        
        ! exhange the ghost values of a vector
        subroutine exch_ghost_values(vec)
            ! input / output
            complex(dp), intent(inout) :: vec(:)                                ! input vector, later overwritten with exchanged value
            
            ! local variables
            complex(dp), allocatable :: exch_vec(:,:)                           ! to exchange the input vector
            integer :: n_modes                                                  ! number of modes to exchange
            
            ! set up n_modes
            n_modes = size(vec)
            
            ! allocate the exchange vectors
            if (grp_rank.eq.0) then
                allocate(exch_vec(n_modes,1))
            else if (grp_rank+1.eq.grp_n_procs) then
                allocate(exch_vec(n_modes,1))
                exch_vec(:,1) = vec                                             ! fill with vec
            else
                allocate(exch_vec(n_modes,2))
                exch_vec(:,1) = vec                                             ! fill with vec
            end if
            
            ! exchange the last interpolated X_vec
            ierr = get_ghost_arr(exch_vec,1)
            CHCKERR('')
            
            ! copy value to new interpolated X_vec
            vec = exch_vec(:,size(exch_vec,2))
        end subroutine exch_ghost_values
    end function plot_X_vec_GP
    
    ! HDF5 version of plot_X_vec
    ! Plots  an  eigenfunction  in  3D  real space  by  determining  a  grid  in
    ! E(quilibrium)  coordinates   that  covers  the  entire   device  and  then
    ! calculating the  real perturbation on  that grid by inverting  the Fourier
    ! transform that defined the modes for which are solved.
    ! Since the equations have  been solved for a single field  line, a trick is
    ! applied  setting up  a modified  equilibrium and  a modified  perturbation
    ! grid,  which   have  the  same   normal  variables  as   their  unmodified
    ! counterparts, but the  angular variables are modified to  cover the entire
    ! space.  As only  the normal  perturbation  is calculated,  using only  the
    ! equilibrium variables flux_q_FD or q_saf_FD, this is not wrong.
    ! Alternatively, by using the optional argument "follow_B", the modes can be
    ! plot along  the magnetic  field lines  which have been  used to  solve the
    ! system of equations. (TO BE IMPLEMENTED)
    integer function plot_X_vec_HDF5(grid_eq,eq,grid_X,X,X_id,job_id,n_t,&
        &follow_B) result(ierr)
        use eq_vars, only: eq_type
        use grid_vars, only: create_grid, destroy_grid, grid_type
        use X_vars, only: X_type
        use num_vars, only: n_theta_plot, n_zeta_plot, grp_n_procs, grp_rank
        use utilities, only: interp_fun
        use grid_ops, only: calc_XYZ_grid, coord_F2E, coord_E2F
        use sol_ops, only: calc_real_XUQ
        use output_ops, only: print_HDF5
        
        character(*), parameter :: rout_name = 'plot_X_vec_HDF5'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        integer, intent(in) :: n_t(2)                                           ! nr. of time steps in quarter period, nr. of quarter periods
        logical, intent(in), optional :: follow_B                               ! .true. if to be plot only along B
        
        ! local variables
        integer :: id                                                           ! counter
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        type(grid_type) :: grid_X_mod                                           ! modified X grid (see remark above)
        type(grid_type) :: grid_eq_mod                                          ! modified eq grid (see remark above)
        real(dp), allocatable :: X_plot_ind(:,:,:)                              ! individual version of X_plot
        real(dp), allocatable :: Y_plot_ind(:,:,:)                              ! individual version of Y_plot
        real(dp), allocatable :: Z_plot_ind(:,:,:)                              ! individual version of Z_plot
        real(dp), allocatable :: X_plot(:,:,:,:)                                ! x values of plot
        real(dp), allocatable :: Y_plot(:,:,:,:)                                ! y values of plot
        real(dp), allocatable :: Z_plot(:,:,:,:)                                ! z values of plot
        real(dp), allocatable :: f_plot(:,:,:,:)                                ! the function to plot
        integer :: tot_dim(4), grp_dim(4), grp_offset(4)                        ! total dimensions, group dimensions and group offset
        logical :: follow_B_loc                                                 ! local copy of follow_B
        character(len=max_str_ln) :: var_name                                   ! name of variable that is plot
        character(len=max_str_ln) :: file_name                                  ! name of file
        character(len=max_str_ln) :: description                                ! description
        real(dp) :: time_frac                                                   ! fraction of Alfvén time
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Started the plot of the Eigenfunction')
        call lvl_ud(1)
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! calculate  omega = sqrt(X_val)  and make  sure to select  the decaying
        ! solution
        omega = sqrt(X%val(X_id))
        if (imagpart(omega).gt.0) omega = - omega                               ! exploding solution, not the decaying one
        
        ! set up local follow_B
        follow_B_loc = .false.
        if (present(follow_B)) then
            if (follow_B) follow_B_loc = follow_B
        end if
        
        ! Set  up total dimensions,  group dimensions  and group offset  for the
        ! final plot. The normal variables are taken from the perturbation grid.
        tot_dim = [n_theta_plot,n_zeta_plot,grid_X%n(3),product(n_t)]
        grp_dim = [n_theta_plot,n_zeta_plot,grp_n_r_X_loc,product(n_t)]
        grp_offset = [0,0,grid_X%i_min-1,product(n_t)]
        
        ! set up var_name, file_name and description
        var_name = 'Solution vector X_vec'
        file_name = 'X_vec_'//trim(i2str(job_id))//'_'//trim(i2str(X_id))
        description = 'Job '//trim(i2str(job_id))//' - Solution vector X_vec &
            &for Eigenvalue '//trim(i2str(X_id))//' with omega = '&
            &//trim(r2str(realpart(omega)))
        
        ! allocate variables holding X, Y and Z of plot
        allocate(X_plot(n_theta_plot,n_zeta_plot,grp_n_r_X_loc,product(n_t)))
        allocate(Y_plot(n_theta_plot,n_zeta_plot,grp_n_r_X_loc,product(n_t)))
        allocate(Z_plot(n_theta_plot,n_zeta_plot,grp_n_r_X_loc,product(n_t)))
        
        ! create modified perturbation grid (only group quantities needed)
        ierr = create_grid(grid_X_mod,&
            &[n_theta_plot,n_zeta_plot,grid_X%n(3)],&
            &[grid_X%i_min,grid_X%i_min+grp_n_r_X_loc-1])
        CHCKERR('')
        
        ! create modified equilibrium grid (only group quantities needed)
        ierr = create_grid(grid_eq_mod,&
            &[n_theta_plot,n_zeta_plot,grid_eq%n(3)],&
            &[grid_eq%i_min,grid_eq%i_max])
        CHCKERR('')
        
        ! Set modified perturbation and modified equilibrium grid coordinates.
        write(*,*) 'THE GRID SHOULD BE CALCULATED ONLY ONCE AS IT DOESNT CHANGE'
        if (follow_B) then
            ierr = 1
            CHCKERR('NOT YET IMPLEMENTED')
            !!!! SHOULD BE SIMILAR TO THE  GRID PLOT, WITH ADDITIONALLY THE !!!!
            !!!! SOLUTION VECTOR SUPERIMPOSED ALONG THE GRID LINES          !!!!
        else
            ! modified X normal F variable taken from X grp_r_F
            grid_X_mod%grp_r_F = grid_X%grp_r_F(1:grp_n_r_X_loc)
            
            ! modified eq normal F variable taken from eq grp_r_F
            grid_eq_mod%grp_r_F = grid_eq%grp_r_F
            
            ! modified X  theta equidistant and saved in E  variable (used later
            ! to calculate X, Y and Z)
            if (n_theta_plot.eq.1) then
                grid_X_mod%theta_E = 0.0_dp
            else
                do id = 1,n_theta_plot
                    grid_X_mod%theta_E(id,:,:) = &
                        &pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)              ! starting from pi gives nicer plots
                end do
            end if
            
            ! modified X  zeta equidistant and  saved in E  variable (used later
            ! to calculate X, Y and Z)
            if (n_zeta_plot.eq.1) then
                grid_X_mod%zeta_E = 0.0_dp
            else
                do id = 1,n_zeta_plot
                    grid_X_mod%zeta_E(:,id,:) = &
                        &(id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
                end do
            end if
            
            ! convert the X  normal values to eq normal values  and save them in
            ! the E variables (used later to calculate X, Y and Z)
            !!ierr = coord_X2eq(eq,X,grid_X_mod%grp_r_F,grid_X_mod%grp_r_E)
            CHCKERR('')
            
            ! convert all equilibrium coordinates to X coordinates and save them
            ! in the Flux variables (used later to calculate normal component of
            ! the plasma perturbation)
            write(*,*) '!?!?!?! CAN WE USE E eq to F X ?! ?!??!? !? !?'
            write(*,*) '!?!?!?!??????????????????????!!!!!!!!!!!!!!!!!'
            ierr = coord_E2F(eq,grid_eq,X,&
                &grid_X_mod%grp_r_E,grid_X_mod%theta_E,grid_X_mod%zeta_E,&
                &grid_X_mod%grp_r_F,grid_X_mod%theta_F,grid_X_mod%zeta_F)
            CHCKERR('')
            
            !write(*,*) 'theta_E'
            !call print_HDF5('X','X',grid_X_mod%theta_E)
            !read(*,*)
            !write(*,*) 'theta_F'
            !call print_HDF5('X','X',grid_X_mod%theta_F)
            !read(*,*)
            !write(*,*) 'zeta_E'
            !call print_HDF5('X','X',grid_X_mod%zeta_E)
            !read(*,*)
            !write(*,*) 'zeta_F'
            !call print_HDF5('X','X',grid_X_mod%zeta_F)
            !read(*,*)
            !write(*,*) 'grp_r_E = ', grid_X_mod%grp_r_E
            !read(*,*)
            !write(*,*) 'grp_r_F = ', grid_X_mod%grp_r_F
            !read(*,*)
        end if
        
        ! calculate  individual X,Y  and  Z using  the  Equilibrium theta_E  and
        ! zeta_plot, tabulated in the equilibrium normal grid
        ierr = calc_XYZ_grid(grid_X_mod,X_plot_ind,Y_plot_ind,Z_plot_ind)
        CHCKERR('')
        
        !write(*,*) 'X_plot_ind'
        !call print_HDF5('X','X',X_plot_ind)
        !read(*,*)
        !write(*,*) 'Y_plot_ind'
        !call print_HDF5('X','X',Y_plot_ind)
        !read(*,*)
        !write(*,*) 'Z_plot_ind'
        !call print_HDF5('X','X',Z_plot_ind)
        !read(*,*)
        
        ! calculate the function to plot: global mode
        allocate(f_plot(n_theta_plot,n_zeta_plot,grp_n_r_X_loc,product(n_t)))
        
        ! user output
        if (product(n_t).gt.1) call writo('Calculating plot for '//&
            &trim(i2str(product(n_t)))//' time steps')
        
        call lvl_ud(1)
            
        ! For  each time step, produce the 3D plot
        do id = 1,product(n_t)
            ! set time fraction
            if (n_t(1).eq.1) then
                time_frac = (id-1) * 0.25
            else
                time_frac = (mod(id-1,n_t(1))*1._dp/n_t(1) + (id-1)/n_t(1)) * &
                    &0.25
            end if
            
            ! calculate f_plot at this time point
            ! Note: Since for  X (style 1) only the  equilibrium flux quantities
            ! rot_t_FD  or  q_saf_FD  are  needed  in  calc_real_XUQ,  a  little
            ! trick  is  applied  here:  Instead  of  the  equilibrium  grid,  a
            ! modified equilibrium  grid is  passed, which  has the  same normal
            ! characteristics as  the equilibrium grid, but  the desired angular
            ! ones.
            ierr = calc_real_XUQ(grid_eq,eq,grid_X,X,grid_X_mod,X_id,1,&
                &time_frac,f_plot(:,:,:,id))
            CHCKERR('')
            
            !write(*,*) 'f_plot'
            !call print_HDF5('X','X',f_plot(:,:,:,id))
            !read(*,*)
            
            ! user output
            if (product(n_t).gt.1) then
                call writo('Calculated plot for time step '//trim(i2str(id))//&
                    &'/'//trim(i2str(product(n_t))))
            end if
            
            ! set up X, Y and Z of plot
            X_plot(:,:,:,id) = X_plot_ind
            Y_plot(:,:,:,id) = Y_plot_ind
            Z_plot(:,:,:,id) = Z_plot_ind
        end do
        
        call lvl_ud(-1)
        
        ! print using HDF5
        call print_HDF5([var_name],file_name,f_plot,tot_dim,grp_dim=grp_dim,&
            &grp_offset=grp_offset,X=X_plot,Y=Y_plot,Z=Z_plot,col=2,&
            &description=description)
        
        ! deallocate
        deallocate(X_plot_ind,Y_plot_ind,Z_plot_ind)
        deallocate(X_plot,Y_plot,Z_plot,f_plot)
        
        ! destroy grid
        call destroy_grid(grid_X_mod)
        
        ! user output
        call lvl_ud(-1)
        call writo('Finished the plot of the Eigenfunction')
    end function plot_X_vec_HDF5
    
    ! normalizes Perturbation quantities using the normalization constants
    subroutine normalize_X_vars(X)
        use X_vars, only: X_type
        use eq_vars, only: psi_0
        
        ! input / output
        type(X_type), intent(inout) :: X                                        ! perturbation
        
        ! scale the quantities
        X%max_flux_E = X%max_flux_E/psi_0
        X%max_flux_F = X%max_flux_F/psi_0
    end subroutine normalize_X_vars
end module
