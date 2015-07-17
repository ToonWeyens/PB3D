!------------------------------------------------------------------------------!
!   Operations considering perturbation quantities                             !
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    use grid_vars, onlY: grid_type, dealloc_grid
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type

    implicit none
    private
    public prepare_X, solve_EV_system, calc_PV, calc_KV, calc_U, &
        &calc_magn_ints, print_output_X, resonance_plot, calc_res_surf

contains
    ! prepare the matrix elements by calculating  KV_i and PV_i, which then will
    ! have to be integrated, with a complex exponential weighting function
    integer function prepare_X(grid_eq,eq,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, plot_resonance, grp_nr
        use X_vars, only: X_type, create_X
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'prepare_X'
        
        ! input / output
        type(grid_type) :: grid_eq                                              ! equilibrium grid variables
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables
        integer :: m, k                                                         ! counters
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle in flux coordinates
        character(len=10) :: integrand_name                                     ! name of integrand in field-line averages
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Preparing perturbation variables')
        call lvl_ud(1)
        
        ! create perturbation
        call create_X(grid_eq,X)
        
        ! tests
        ierr = check_modes(eq,X)
        CHCKERR('')
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_F) then
            ang_par_F => grid_eq%theta_F
        else
            ang_par_F => grid_eq%zeta_F
        end if
        
        ! plot resonances if requested
        if (plot_resonance .and. grp_nr.eq.0) then
            ierr = resonance_plot(eq,grid_eq,X)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! set exp_ang_par_F
        if (use_pol_flux_F) then
            do m = 1,X%n_mod
                do k = m,X%n_mod
                    X%exp_ang_par_F(:,:,:,c([k,m],.true.,X%n_mod)) = &
                        &exp(iu*(k-m)*ang_par_F)
                end do
            end do
        else
            do m = 1,X%n_mod
                do k = m,X%n_mod
                    X%exp_ang_par_F(:,:,:,c([k,m],.true.,X%n_mod)) = &
                        &exp(iu*(m-k)*ang_par_F)
                end do
            end do
        end if
        
        ! set up integrand_name
        if (use_pol_flux_F) then
            integrand_name = '(k-m)theta'
        else
            integrand_name = '(m-k)zeta'
        end if
        
        call writo('Setting up tables of field line averages &
            &<V e^i['//trim(integrand_name)//']>')
        
        call lvl_ud(1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        ierr = calc_U(eq,grid_eq,met,X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! Calculate PV_i for all (k,m) pairs and n_r (equilibrium) values of the
        ! normal coordinate
        call writo('Calculating PV...')
        call lvl_ud(1)
        call calc_PV(eq,grid_eq,met,X)
        call lvl_ud(-1)
        
        ! Calculate KV_i for all (k,m) pairs and n_r (equilibrium) values of the
        ! normal coordinate
        call writo('Calculating KV...')
        call lvl_ud(1)
        call calc_KV(eq,grid_eq,met,X)
        call lvl_ud(-1)
        
        call lvl_ud(-1)
        
        call writo('Done setting up tables')
        
        ! clean up
        nullify(ang_par_F)
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation variables prepared')
    end function prepare_X
    
    ! plot  q-profile  or iota-profile  in  flux coordinates  with nq-m  = 0  or
    ! n-iotam = 0 indicate if requested
    integer function resonance_plot(eq,grid,X) result(ierr)
        use num_vars, only: use_pol_flux_F, grp_rank, no_plots, n_theta_plot, &
            &n_zeta_plot
        use utilities, only: calc_zero_NR, interp_fun
        use grid_vars, only: create_grid, dealloc_grid
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        use grid_ops, only: trim_grid, calc_XYZ_grid, calc_eqd_grid, coord_F2E
        use MPI_utilities, only: get_ser_var
        use output_ops, only: plot_HDF5
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        
        ! local variables (not to be used in child functions)
        integer :: kd                                                           ! counter
        integer :: n_mod_loc                                                    ! local n_mod
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        real(dp), allocatable :: x_vars(:,:)                                    ! for plotting
        real(dp), allocatable :: y_vars(:,:)                                    ! for plotting
        character(len=max_str_ln) :: plot_title, file_name                      ! name of plot, of file
        real(dp), allocatable :: jq(:)                                          ! saf. fac. or rot. transf. in Flux coords.
        integer :: n_r                                                          ! total number of normal points
        integer :: plot_dim(4)                                                  ! plot dimensions (total = group because only group masters)
        type(grid_type) :: grid_trim                                            ! trimmed version of grid
        type(grid_type) :: grid_plot                                            ! grid for plotting
        real(dp), allocatable :: r_plot_E(:)                                    ! normal E coordinates of plot
        real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)            ! pol. and tor. angle of plot
        real(dp), allocatable :: X_plot(:,:,:,:), Y_plot(:,:,:,:), &
            &Z_plot(:,:,:,:)                                                    ! X, Y and Z of plot of all surfaces
        real(dp), allocatable :: vars(:,:,:,:)                                  ! variable to plot
        character(len=max_str_ln), allocatable :: plot_titles(:)                ! name of plots
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! get trimmed grid
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
        ! initialize variables
        n_r = grid_trim%n(3)
        
        ! print user output
        if (use_pol_flux_F) then
            call writo('Plotting safety factor q and resonant surfaces &
                &q = m/n...')
            plot_title = 'safety factor q'
            file_name = 'q_saf'
        else
            call writo('Plotting rotational transform iota and resonant &
                &surfaces iota = n/m...')
            plot_title = 'rotational transform iota'
            file_name = 'rot_t'
        end if
        
        call lvl_ud(1)
        
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        
        ! find resonating surfaces
        ierr = calc_res_surf(grid,eq,X,res_surf,info=.true.,jq=jq,&
            &tol_NR=1.E-8_dp,max_it_NR=5000)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! only group master
        if (grp_rank.eq.0) then
            ! set local n_mod
            n_mod_loc = size(res_surf,1)
            
            ! initialize x_vars and y_vars
            allocate(x_vars(n_r,n_mod_loc+1)); x_vars = 0
            allocate(y_vars(n_r,n_mod_loc+1)); y_vars = 0
            
            ! set x_vars and y_vars for first column
            x_vars(:,1) = grid_trim%r_F
            y_vars(:,1) = jq(:)
            
            ! set x_vars and y_vars for other columns
            do kd = 1,n_mod_loc
                x_vars(:,kd+1) = res_surf(kd,2)
                y_vars(n_r,kd+1) = res_surf(kd,3)
            end do
            
            ! user message
            call writo('Plot results using HDF5')
            
            call lvl_ud(1)
            
            ! set up pol. and tor. angle for plot
            allocate(theta_plot(n_theta_plot,n_zeta_plot,1))
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,1))
            ierr = calc_eqd_grid(theta_plot,1*pi,3*pi,1)                        ! starting from pi gives nicer results
            CHCKERR('')
            ierr = calc_eqd_grid(zeta_plot,0*pi,2*pi,2)
            CHCKERR('')
            
            ! set up vars
            allocate(vars(n_theta_plot,n_zeta_plot,1,n_mod_loc))
            do kd = 1,n_mod_loc
                vars(:,:,:,kd) = y_vars(n_r,kd+1)
            end do
            
            ! set up plot titles
            allocate(plot_titles(n_mod_loc))
            if (use_pol_flux_F) then                                            ! n is fixed and m = m/n n
                do kd = 1,n_mod_loc
                    plot_titles(kd) = trim(plot_title)//' for m,n = '//&
                        &trim(i2str(nint(res_surf(kd,3)*X%n(1))))//','//&
                        &trim(i2str(X%n(1)))
                end do
            else                                                                ! m is fixed and n = n/m m
                do kd = 1,n_mod_loc
                    plot_titles(kd) = trim(plot_title)//' for m,n = '//&
                        &trim(i2str(X%m(1)))//','//&
                        &trim(i2str(nint(res_surf(kd,3)*X%m(1))))
                end do
            end if
            
            ! set dimensions
            plot_dim = [n_theta_plot,n_zeta_plot,1,n_mod_loc]
            
            ! calculate normal vars in Equilibrium coords.
            allocate(r_plot_E(n_mod_loc))
            ierr = coord_F2E(grid,eq,x_vars(n_r,2:n_mod_loc+1),r_plot_E,&
                &r_F_array=grid%r_F,r_E_array=grid%r_E)
            CHCKERR('')
            
            ! create plot grid
            ierr = create_grid(grid_plot,plot_dim(1:3))
            CHCKERR('')
            grid_plot%theta_E = theta_plot
            grid_plot%zeta_E = zeta_plot
            
            ! set up plot X, Y and Z
            allocate(X_plot(n_theta_plot,n_zeta_plot,1,n_mod_loc))
            allocate(Y_plot(n_theta_plot,n_zeta_plot,1,n_mod_loc))
            allocate(Z_plot(n_theta_plot,n_zeta_plot,1,n_mod_loc))
            
            ! loop over all resonant surfaces to calculate X, Y and Z values
            do kd = 1,n_mod_loc
                ! set grp_r_E of plot grid
                grid_plot%grp_r_E = r_plot_E(kd)
                
                ierr = calc_XYZ_grid(grid_plot,X_plot(:,:,:,kd),&
                    &Y_plot(:,:,:,kd),Z_plot(:,:,:,kd))
                CHCKERR('')
            end do
            
            ! print using HDF5
            call plot_HDF5(plot_titles,file_name,vars,X=X_plot,Y=Y_plot,&
                &Z=Z_plot,col=1,description='resonant surfaces')
            
            ! deallocate local variables
            deallocate(vars)
            deallocate(theta_plot,zeta_plot,r_plot_E)
            deallocate(X_plot,Y_plot,Z_plot)
            
            call lvl_ud(-1)
            
            call writo('Plot results using GNUPlot')
            
            call lvl_ud(1)
            
            ! rescale x_vars by max_flux_F/2*pi
            if (use_pol_flux_F) then
                x_vars = x_vars*2*pi/max_flux_p_F
            else
                x_vars = x_vars*2*pi/max_flux_t_F
            end if
            
            ! print to file
            call print_GP_2D(plot_title,file_name,y_vars,x=x_vars,draw=.false.)
            
            ! plot using GNUPlot
            call draw_GP(plot_title,file_name,file_name,n_mod_loc+1,1,&
                &.false.)
            
            call lvl_ud(-1)
            
            call writo('Plot 2D results using HDF5')
            
            call lvl_ud(1)
            
            ! initialize plot dimensions
            plot_dim = [1,n_mod_loc,1,0]
            
            ! set up X, Y and Z
            allocate(X_plot(plot_dim(1),plot_dim(2),plot_dim(3),1))
            allocate(Y_plot(plot_dim(1),plot_dim(2),plot_dim(3),1))
            allocate(Z_plot(plot_dim(1),plot_dim(2),plot_dim(3),1))
            
            do kd = 1,n_mod_loc
                X_plot(1,kd,1,1) = x_vars(1,kd+1)
                Y_plot(1,kd,1,1) = kd-1
                Z_plot(1,kd,1,1) = 1._dp
            end do
            
            ! plot 2D (to be used with plots of harmonics in sol_ops)
            call plot_HDF5('resonating surfaces','res_surf',&
                &Z_plot(:,:,:,1),x=X_plot(:,:,:,1),y=Y_plot(:,:,:,1))
            
            call lvl_ud(-1)
            
            ! clean up
            call dealloc_grid(grid_plot)
        end if
            
        ! clean up
        call dealloc_grid(grid_trim)
        
        call lvl_ud(-1)
    end function resonance_plot
    
    ! Calculates resonating flux surfaces for the perturbation modes. The output
    ! consists of mode  number, resonating normal position and  the fraction n/m
    ! or m/n for  those modes for which  a solution is found that  is within the
    ! plasma range.
    ! Optionally, user  messages can be  displayed with info  and Newton-Rhapson
    ! tolerance and  max. nr. of  iterations. Also,  the total safety  factor or
    ! rotational transform can be returned to the group master.
    ! Note: Has to be called by all processes
    integer function calc_res_surf(grid,eq,X,res_surf,info,jq,tol_NR,&
        &max_it_NR) result(ierr)
        use num_vars, only: use_pol_flux_F, tol_NR_loc => tol_NR, &
            &max_it_NR_loc => max_it_NR
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        use utilities, only: calc_zero_NR, interp_fun
        use grid_vars, only: dealloc_grid
        use grid_ops, only: trim_grid
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(X_type), intent(in) :: X                                           ! perturbation variables
        real(dp), intent(inout), allocatable :: res_surf(:,:)                   ! resonant surface
        logical, intent(in), optional :: info                                   ! if info is displayed
        real(dp), intent(inout), optional, allocatable :: jq(:)                 ! either safety factor or rotational transform
        real(dp), intent(in), optional :: tol_NR                                ! optional tolerance for Newton-Rhapson iteration
        integer, intent(in), optional :: max_it_NR                              ! optional max. nr. of iterations for Newton-Rhapson iteration
        
        ! local variables (also used in child functions)
        real(dp) :: mnfrac_fun                                                  ! fraction m/n or n/m to determine resonant flux surface
        
        ! local variables (not to be used in child functions)
        integer :: kd                                                           ! counter
        integer :: kd_loc                                                       ! local kd
        logical :: info_loc                                                     ! local version of info
        integer :: istat                                                        ! status
        real(dp), allocatable :: res_surf_loc(:,:)                              ! local copy of res_surf
        real(dp), allocatable :: jq_tot(:,:)                                    ! saf. fac. or rot. transf. and derivs. in Flux coords.
        real(dp), allocatable :: jq_loc(:)                                      ! local version of jq
        real(dp) :: norm_factor                                                 ! normalization factor to make normal coordinate 0..1
        type(grid_type) :: grid_trim                                            ! trimmed version of grid
        real(dp) :: tol_NR_old                                                  ! old value of tol_NR
        integer :: max_it_NR_old                                                ! old value of max_it_NR
        
        ! initialize ierr
        ierr = 0
        
        ! set local info
        info_loc = .false.
        if (present(info)) info_loc = info
        
        ! get trimmed grid
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
        ! get serial version of safety factor or rot. transform
        allocate(jq_tot(grid_trim%n(3),0:2))
        if (use_pol_flux_F) then
            do kd = 0,2
                ierr = get_ser_var(eq%q_saf_FD(1:grid_trim%grp_n_r,kd),jq_loc,&
                    &scatter=.true.)
                CHCKERR('')
                jq_tot(:,kd) = jq_loc
            end do
        else
            do kd = 0,2
                ierr = get_ser_var(eq%rot_t_FD(1:grid_trim%grp_n_r,kd),jq_loc,&
                    &scatter=.true.)
                CHCKERR('')
                jq_tot(:,kd) = jq_loc
            end do
        end if
        if (present(jq)) then
            allocate(jq(grid_trim%n(3)))
            jq = jq_tot(:,0)
        end if
        
        ! initialize local res_surf
        allocate(res_surf_loc(1:X%n_mod,3))
        
        ! save old tol_NR and max_it_NR and modify them if needed
        tol_NR_old = tol_NR_loc
        max_it_NR_old = max_it_NR_loc
        if (present(tol_NR)) tol_NR_loc = tol_NR
        if (present(max_it_NR)) max_it_NR_loc = max_it_NR
        
        ! calculate normalization factor max_flux / 2pi
        if (use_pol_flux_F) then
            norm_factor = max_flux_p_F/(2*pi)
        else
            norm_factor = max_flux_t_F/(2*pi)
        end if
        
        ! loop over all modes (and shift the index in x and y_vars by 1)
        kd_loc = 1
        do kd = 1, X%n_mod
            ! find place  where q  = m/n or  iota = n/m  in Flux  coordinates by
            ! solving q-m/n = 0 or iota-n/m=0, using the functin jq_fun
            
            ! set up mnfrac for function
            if (use_pol_flux_F) then
                mnfrac_fun = 1.0_dp*X%m(kd)/X%n(kd)
            else                             
                mnfrac_fun = 1.0_dp*X%n(kd)/X%m(kd)
            end if
            
            ! calculate zero using Newton-Rhapson
            res_surf_loc(kd_loc,1) = kd
            res_surf_loc(kd_loc,3) = mnfrac_fun
            istat = calc_zero_NR(res_surf_loc(kd_loc,2),jq_fun,jq_dfun,&
                &(maxval(grid_trim%r_F)+minval(grid_trim%r_F))/2)               ! guess halfway between minimum and maximum normal range
            
            ! intercept error
            if (istat.ne.0) then
                call lvl_ud(1)
                if (info_loc) call writo('Error intercepted: Couldn''t &
                    &find resonating surface for (n,m) = ('//&
                    &trim(i2str(X%n(kd)))//','//trim(i2str(X%m(kd)))//')')
                call lvl_ud(-1)
            else if (res_surf_loc(kd_loc,2).lt.minval(grid_trim%r_F) .or. &
                &res_surf_loc(kd_loc,2).gt.maxval(grid_trim%r_F)) then
                if (info_loc) call writo('Mode (n,m) = ('//&
                    &trim(i2str(X%n(kd)))//','//trim(i2str(X%m(kd)))//&
                    &') does not resonate in plasma')
            else
                if (info_loc) call writo('Mode (n,m) = ('//&
                    &trim(i2str(X%n(kd)))//','//trim(i2str(X%m(kd)))//&
                    &') resonates in plasma at normalized flux surface '//&
                    &trim(r2str(res_surf_loc(kd_loc,2)/norm_factor)))
                kd_loc = kd_loc + 1                                             ! advance kd_loc
            end if
        end do
        
        ! set res_surf from local copy
        allocate(res_surf(kd_loc-1,3))
        res_surf = res_surf_loc(1:kd_loc-1,:)
        
        ! recover old tol_NR and max_it_NR
        tol_NR_loc = tol_NR_old
        max_it_NR_loc = max_it_NR_old
        
        ! deallocate local variables
        deallocate(jq_tot)
        call dealloc_grid(grid_trim)
    contains
        ! Returns q-m/n or  iota-n/m in Flux coordinates, used to  solve for q =
        ! m/n or iota = n/m.
        ! WARNING: This routine requires that jq_tot's derivatives be calculated
        ! up to order 1. This is NOT checked!
        real(dp) function jq_fun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            integer :: i_min, i_max                                             ! index of minimum and maximum value of x
            real(dp) :: x_min, x_max                                            ! minimum and maximum value of x
            
            ! initialize res
            res = 0
            
            ! set up min. and max index and value
            x_min = minval(grid_trim%r_F)
            x_max = maxval(grid_trim%r_F)
            i_min = minloc(grid_trim%r_F,1)
            i_max = maxloc(grid_trim%r_F,1)
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.x_min) then                                               ! point requested lower than minimum x
                ! extrapolate variable from minimum value
                res = jq_tot(i_min,0) - mnfrac_fun + jq_tot(i_min,1)*(pt-x_min)
            else if (pt.gt.x_max) then                                          ! point requested higher than maximum x
                ! extrapolate variable from maximum value
                res = jq_tot(i_max,0) - mnfrac_fun + jq_tot(i_max,1)*(pt-x_max)
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                istat = interp_fun(res,jq_tot(:,0)-mnfrac_fun,pt,grid_trim%r_F)
            end if
        end function jq_fun
        
        ! returns d(q-m/n)/dr  in Flux coordinates,  used to  solve for q  = m/n
        ! WARNING: This routine requires that jq_tot's derivatives be calculated
        ! up to order 2. This is NOT checked!
        real(dp) function jq_dfun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            integer :: i_min, i_max                                             ! index of minimum and maximum value of x
            real(dp) :: x_min, x_max                                            ! minimum and maximum value of x
            
            ! initialize res
            res = 0
            
            ! set up min. and max index and value
            x_min = minval(grid_trim%r_F)
            x_max = maxval(grid_trim%r_F)
            i_min = minloc(grid_trim%r_F,1)
            i_max = maxloc(grid_trim%r_F,1)
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.x_min) then                                               ! point requested lower than minimum x
                ! extrapolate variable from minimum value
                res = jq_tot(i_min,1) + jq_tot(i_min,2)*(pt-x_min)
            else if (pt.gt.x_max) then                                          ! point requested higher than maximum x
                ! extrapolate variable from maximum value
                res = jq_tot(i_max,0) + jq_tot(i_max,2)*(pt-x_max)
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                istat = interp_fun(res,jq_tot(:,1),pt,x=grid_trim%r_F)
            end if
        end function jq_dfun
    end function calc_res_surf
    
    ! Checks whether |nq-m|/|n|  < tol and |nq-m|/|m| < tol  (or |q-iotam|/|m| <
    ! tol and |n-iotam|/|n| < tol) is satisfied in some part of the plasma, with
    ! tol << 1.
    ! The condition is determined by the sign of q (or iota) and given by:
    !   max(min_q-tol,min_q/(1+tol)) < m/n < min(max_q+tol,max_q/(1-tol)), q>0
    !   max(min_q-tol,min_q/(1-tol)) < m/n < min(max_q+tol,max_q/(1+tol)), q<0
    ! (or replacing q by iota and m/n by n/m).
    ! [MPI] Parts by all processes, parts only by global master
    integer function check_modes(eq,X) result(ierr)
        use MPI_utilities, only: get_ser_var
        use num_vars, only: glb_rank, use_pol_flux_F, eq_style, tol_norm_r
        
        character(*), parameter :: rout_name = 'check_modes'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp) :: min_jq, max_jq                                              ! min. and max. values of q or iota
        integer :: pmone                                                        ! plus or minus one
        integer :: pmone2                                                       ! plus or minus one
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: jq_name                                    ! either safety factor or rotational transform
        character(len=1) :: mode_name                                           ! either n or m
        real(dp), allocatable :: ser_jq(:)                                      ! serial q_saf_E(:,0) or rot_t_E(:,0)
        real(dp) :: lim_lo, lim_hi                                              ! lower and upper limit on n/m (or m/n)
        real(dp) :: min_sec_X, max_sec_X                                        ! min. and max. of m_X (if pol. flux) or n_X (if tor. flux)
        
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
            call lvl_ud(1)
            
            ! user output
            call writo('The tolerance used is '//trim(r2strt(tol_norm_r))//&
                &'...')
            
            ! set  up  plus  minus  one  to  convert  from Equilibrium  to  Flux
            ! coordinates
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
            
            ! set min_jq and max_jq in Flux coordinate system
            min_jq = minval(pmone*ser_jq)
            max_jq = maxval(pmone*ser_jq)
            
            ! set up jq name
            if (use_pol_flux_F) then
                jq_name = 'safety factor'
                mode_name = 'n'
            else
                jq_name = 'rotational transform'
                mode_name = 'm'
            end if
            
            ! set  up plus  minus one  2, according  to the  sign of  the safety
            ! factor
            if (min_jq.lt.0 .and. max_jq.lt.0) then
                pmone2 = -1
            else if (min_jq.gt.0 .and. max_jq.gt.0) then
                pmone2 = 1
            else
                err_msg = trim(jq_name)//' cannot change sign'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! calculate upper and lower limits
            lim_lo = max(min_jq-tol_norm_r,min_jq/(1+pmone2*tol_norm_r))
            lim_hi = min(max_jq+tol_norm_r,max_jq/(1-pmone2*tol_norm_r))
            
            ! initialize min_sec_X and max_sec_X
            min_sec_X = huge(1._dp)
            max_sec_X = -huge(1._dp)
            
            ! for every mode (n,m) check whether  m/n is inside the range of
            ! q values or n/m inside the range of iota values
            do id = 1, X%n_mod
                ! check if limits are met
                if (use_pol_flux_F) then
                    if (X%m(id)*1.0/X%n(id).lt.lim_lo .or. &
                        &X%m(id)*1.0/X%n(id).gt.lim_hi) then
                        call writo('for (n,m) = ('//trim(i2str(X%n(id)))//&
                            &','//trim(i2str(X%m(id)))//'), there is no range &
                            &the plasma where the ratio |n q - m| << |n|,|m| &
                            &is  met')
                        ierr = 1
                        err_msg = 'Choose m and n so that |n q - m| << |n|,|m|'
                        CHCKERR(err_msg)
                    else
                        if (X%n(id)*lim_lo.lt.min_sec_X) &
                            &min_sec_X = X%n(id)*lim_lo
                        if (X%n(id)*lim_hi.gt.max_sec_X) &
                            &max_sec_X = X%n(id)*lim_hi
                    end if
                else
                    if (X%n(id)*1.0/X%m(id).lt.lim_lo .or. &
                        &X%n(id)*1.0/X%m(id).gt.lim_hi) then
                        call writo('for (n,m) = ('//trim(i2str(X%n(id)))//&
                            &','//trim(i2str(X%m(id)))//'), there is no range &
                            &the plasma where the ratio &
                            &|n - iota m| << |m|,|n| is  met')
                        ierr = 1
                        err_msg = 'Choose m and n so that &
                            &|n - iota m| << |n|,|m|'
                        CHCKERR(err_msg)
                    else
                        if (X%m(id)*lim_lo.lt.min_sec_X) &
                            &min_sec_X = X%m(id)*lim_lo
                        if (X%m(id)*lim_hi.gt.max_sec_X) &
                            &max_sec_X = X%m(id)*lim_hi
                    end if
                end if
            end do
            
            ! output message
            call writo('The modes are all within the allowed range of '//&
                &trim(i2str(ceiling(min_sec_X)))//' < '//mode_name//' < '//&
                &trim(i2str(floor(max_sec_X)))//'...')
            
            call lvl_ud(-1)
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
    
    ! calculate  ~PV_(k,m)^i  (pol.  flux)  or ~PV_(l,n)^i  (tor.  flux) at  all
    ! eq grp_n_r values
    ! (see [ADD REF] for details)
    subroutine calc_PV(eq,grid,met,X)
        use num_vars, only: use_pol_flux_F
        use eq_vars, only: vac_perm
        use utilities, only: c
        
        ! use input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:,:)                                 ! common factor |nabla psi|^2/(J^2*B^2)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        integer :: c1                                                           ! value of c, to avoid compiler hang
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:) => null()                                 ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:) => null()                               ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:) => null()                               ! h^psi,psi
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        ! lower metric factors
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        ! upper metric factors
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up common factor for PV_i
        allocate(com_fac(grid%n(1),grid%n(2),grid%grp_n_r))
        com_fac = h22/(g33*vac_perm)
        
        ! set up fac_n and fac_m
        allocate(fac_n(grid%grp_n_r),fac_m(grid%grp_n_r))
        if (use_pol_flux_F) then
            fac_n = eq%q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = eq%rot_t_FD(:,0)
        end if
        
        ! calculate PV_0 and PV_2 (Hermitian)
        do m = 1,X%n_mod
            do k = m,X%n_mod
                ! calculate PV_0
                X%PV_0(:,:,:,c([k,m],.true.,X%n_mod)) = com_fac*&
                    &(X%DU_0(:,:,:,m) - eq%S*J - eq%sigma/(com_fac*J)) * &
                    &(conjg(&
                    &X%DU_0(:,:,:,k)) - eq%S*J - eq%sigma/(com_fac*J)) - &
                    &eq%sigma/J * (eq%S*J + eq%sigma/(com_fac*J))
                
                ! add (nq-k)*(nq-m)/(mu_0J^2 |nabla psi|^2) - 2p'kappa_n to PV_0
                do kd = 1,grid%grp_n_r
                    c1 = c([k,m],.true.,X%n_mod)
                    X%PV_0(:,:,kd,c1) = X%PV_0(:,:,kd,c1) + &
                        &(X%n(m)*fac_n(kd)-X%m(m)*fac_m(kd))*&
                        &(X%n(k)*fac_n(kd)-X%m(k)*fac_m(kd)) / &
                        &( vac_perm*J(:,:,kd)**2*h22(:,:,kd) ) - &
                        &2*eq%pres_FD(kd,1)*eq%kappa_n(:,:,kd)
                end do
                
                ! calculate PV_2
                X%PV_2(:,:,:,c([k,m],.true.,X%n_mod)) = &
                    &com_fac*X%DU_1(:,:,:,m)*conjg(X%DU_1(:,:,:,k))
            end do
        end do
        
        ! PV_1 is not Hermitian
        do m = 1,X%n_mod
            do k = 1,X%n_mod
                ! calculate PV_1
                X%PV_1(:,:,:,c([k,m],.false.,X%n_mod)) = &
                    &com_fac * X%DU_1(:,:,:,m) * &
                    &(conjg(X%DU_0(:,:,:,k)) - eq%S*J - eq%sigma/(com_fac*J))
            end do
        end do
        
        ! deallocate variables
        nullify(J,g33,h22)
    end subroutine calc_PV
    
    ! calculate  ~KV_(k,m)^i  (pol.  flux)  or ~KV_(l,n)^i  (tor.  flux) at  all
    ! eq grp_n_r values
    ! (see [ADD REF] for details)
    subroutine calc_KV(eq,grid,met,X)
        use utilities, only: c
        
        ! use input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:,:)                                 ! common factor |nabla psi|^2/(J^2*B^2)
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:) => null()                                 ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:) => null()                               ! h^alpha,psi
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:) => null()                               ! h^alpha,psi
        
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
        
        ! deallocate variables
        nullify(J,g33,h22)
    end subroutine calc_KV
    
    ! calculate U_m^0, U_m^1 or U_n^0, U_n^1  at eq grp_n_r values of the normal
    ! coordinate, n_par values  of the parallel coordinate and  size_X values of
    ! the poloidal mode number,  or of the toroidal mode number,  as well as the
    ! poloidal derivatives
    integer function calc_U(eq,grid,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, ltest
        use utilities, only: c
        use input_ops, only: get_log, pause_prog
        use eq_vars, only: vac_perm
        
        character(*), parameter :: rout_name = 'calc_U'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        real(dp) :: n_frac                                                      ! nq-m (pol. flux) or n-iotam (tor. flux)
        real(dp), allocatable :: djq(:)                                         ! either q' (pol. flux) or -iota' (tor. flux)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp), allocatable :: mn(:)                                          ! either n*A_0 (pol. flux) or m (tor.flux)
        complex(dp), allocatable :: U_corr(:,:,:,:)                             ! correction to U for a certain (n,m)
        complex(dp), allocatable :: D3U_corr(:,:,:,:)                           ! D_theta U_corr
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle in flux coordinates
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! helper variables for the correction to U
        real(dp), allocatable :: g_frac(:,:,:)                                  ! g_alpha,theta / g_theta,theta
        real(dp), allocatable :: Theta_3(:,:,:), D1Theta_3(:,:,:), &
            &D3Theta_3(:,:,:)                                                   ! Theta^theta and derivatives
        ! jacobian
        real(dp), pointer :: J(:,:,:) => null()                                 ! jac
        real(dp), pointer :: D3J(:,:,:) => null()                               ! D_theta jac
        ! lower metric factors
        real(dp), pointer :: g13(:,:,:) => null()                               ! g_alpha,theta
        real(dp), pointer :: D3g13(:,:,:) => null()                             ! D_theta g_alpha,theta
        real(dp), pointer :: g23(:,:,:) => null()                               ! g_psi,theta
        real(dp), pointer :: D3g23(:,:,:) => null()                             ! D_theta g_psi,theta
        real(dp), pointer :: g33(:,:,:) => null()                               ! g_theta,theta
        real(dp), pointer :: D3g33(:,:,:) => null()                             ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), pointer :: h12(:,:,:) => null()                               ! h^alpha,psi
        real(dp), pointer :: D3h12(:,:,:) => null()                             ! D_theta h^alpha,psi
        real(dp), pointer :: h22(:,:,:) => null()                               ! h^psi,psi
        real(dp), pointer :: D1h22(:,:,:) => null()                             ! D_alpha h^psi,psi
        real(dp), pointer :: D3h22(:,:,:) => null()                             ! D_theta h^psi,psi
        real(dp), pointer :: h23(:,:,:) => null()                               ! h^psi,theta
        real(dp), pointer :: D1h23(:,:,:) => null()                             ! D_alpha h^psi,theta
        real(dp), pointer :: D3h23(:,:,:)=> null()                              ! D_theta h^psi,theta
        
        ! initialize ierr
        ierr = 0
        
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
        allocate(Theta_3(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(D1Theta_3(grid%n(1),grid%n(2),grid%grp_n_r))
        allocate(D3Theta_3(grid%n(1),grid%n(2),grid%grp_n_r))
        
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
        
#if ldebug
                if (ltest) then
                    call writo('Test calculation of DU')
                    if(get_log(.false.,ind=.true.)) then
                        ierr = test_DU()
                        CHCKERR('')
                        call pause_prog(ind=.true.)
                    end if
                end if
#endif
            case (2)                                                            ! HELENA
                ierr = calc_U_HEL()
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! clean up
        nullify(ang_par_F)
        nullify(J,D3J)
        nullify(g13,D3g13,g23,D3g23,g33,D3g33)
        nullify(h12,D3h12,h22,D1h22,D3h22,h23,D1h23,D3h23)
    contains
        ! VMEC version
        subroutine calc_U_VMEC
            ! local variables
            ! extra helper variables for the correction to U
            real(dp), allocatable :: D3g_frac(:,:,:)                            ! D_theta g_frac
            real(dp), allocatable :: D13Theta_3(:,:,:), D33Theta_3(:,:,:)       ! Theta^theta derivatives
            ! extra upper metric factors
            real(dp), pointer :: D13h22(:,:,:) => null()                        ! D^2_alpha,theta h^psi,psi
            real(dp), pointer :: D33h22(:,:,:) => null()                        ! D^2_theta,theta h^psi,psi
            real(dp), pointer :: D13h23(:,:,:) => null()                        ! D^2_alpha,theta h^psi,theta
            real(dp), pointer :: D33h23(:,:,:) => null()                        ! D^2_theta,theta h^psi,theta
            
            ! allocate extra helper variables
            allocate(D3g_frac(grid%n(1),grid%n(2),grid%grp_n_r))
            allocate(D13Theta_3(grid%n(1),grid%n(2),grid%grp_n_r))
            allocate(D33Theta_3(grid%n(1),grid%n(2),grid%grp_n_r))
            
            ! set up extra submatrices
            ! upper metric factors
            D13h22 => met%h_FD(:,:,:,c([2,2],.true.),1,0,1)
            D33h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,2)
            D13h23 => met%h_FD(:,:,:,c([2,3],.true.),1,0,1)
            D33h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,2)
            
            ! loop over the M elements of U_X and DU
            do jd = 1,X%n_mod
                ! set up helper variables
                g_frac = g13/g33
                D3g_frac = D3g13/g33 - g13*D3g33/(g33**2)
                Theta_3 = h23/h22
                D1Theta_3 = D1h23/h22 - h23*D1h22/(h22**2)
                D3Theta_3 = D3h23/h22 - h23*D3h22/(h22**2)
                D13Theta_3 = D13h23/h22 - (D3h23*D1h22+D1h23*D3h22)/(h22**2) &
                    &- h23*D13h22/(h22**2) + 2*h23*D3h22*D1h22/(h22**3)
                D33Theta_3 = D33h23/h22 - 2*D3h23*D3h22/(h22**2) &
                    &- h23*D33h22/(h22**2) + 2*h23*D3h22**2/(h22**3)
                
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    ! set up correction to U
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*&
                        &(g_frac(:,:,kd)*(D3Theta_3(:,:,kd)+&
                        &iu*n_frac*Theta_3(:,:,kd))+D1Theta_3(:,:,kd))
                    ! set up D_theta U correction
                    D3U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*&
                        &(D3g_frac(:,:,kd)*&
                        &(D3Theta_3(:,:,kd)+iu*n_frac*Theta_3(:,:,kd))+&
                        &g_frac(:,:,kd)*&
                        &(D33Theta_3(:,:,kd)+iu*n_frac*D3Theta_3(:,:,kd))+&
                        &D13Theta_3(:,:,kd))
                    ! calculate X%U_0 and X%DU_0
                    X%U_0(:,:,kd,jd) = &
                        &-(h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd))&
                        &+ iu/(mn(jd)*g33(:,:,kd)) * (g13(:,:,kd)*djq(kd) + &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*vac_perm + iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) - &
                        &g23(:,:,kd) )) + U_corr(:,:,kd,jd)
                    X%DU_0(:,:,kd,jd) = -(D3h12(:,:,kd)/h22(:,:,kd) - &
                        &D3h22(:,:,kd)*h12(:,:,kd)/h22(:,:,kd)**2 + djq(kd)) - &
                        &iu*D3g33(:,:,kd)/(mn(jd)*g33(:,:,kd)**2) * &
                        &(g13(:,:,kd)*djq(kd)+ &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*vac_perm + &
                        &iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) &
                        &- g23(:,:,kd) )) + &
                        &iu/(mn(jd)*g33(:,:,kd)) * (D3g13(:,:,kd)*djq(kd) + &
                        &2*D3J(:,:,kd)*J(:,:,kd)*eq%pres_FD(kd,1)*vac_perm + &
                        &iu*n_frac &
                        &* ( D3g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) + &
                        &g13(:,:,kd)*djq(kd) - D3g23(:,:,kd) )) + &
                        &iu*n_frac*X%U_0(:,:,kd,jd) &
                        &+ D3U_corr(:,:,kd,jd)
                    ! calculate X%U_1 and X%DU_1
                    X%U_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,:,kd)/g33(:,:,kd))
                    X%DU_1(:,:,kd,jd) = iu/mn(jd) * n_frac/mn(jd) * &
                        &( D3g13(:,:,kd)/g33(:,:,kd) - &
                        &  g13(:,:,kd)*D3g33(:,:,kd)/g33(:,:,kd)**2 ) + &
                        &iu*n_frac*X%U_1(:,:,kd,jd) 
                end do
            end do
            
            ! clean up
            nullify(D13h22,D33h22,D13h23,D33h23)
        end subroutine calc_U_VMEC
        
        ! HELENA version
        ! Note:  The parallel  derivatives are  done here  numerically. However,
        ! these  must be  translated to  numerical derivatives  in the  poloidal
        ! coordinate. This depends on the flux on which the normal coordinate is
        ! based:
        ! - poloidal flux: parallel deriv. equal to poloidal deriv,
        ! - toroidal flux: parallel deriv. equal to iota * poloidal deriv.
        integer function calc_U_HEL() result(ierr)
            use utilities, only: calc_deriv
            
            character(*), parameter :: rout_name = 'calc_U_HEL'
            
            ! local variablas
            complex(dp), allocatable :: D3var(:,:)                              ! derivative of variable
            
            ! allocate extra helper variables
            allocate(D3var(grid%n(1),grid%n(2)))
            
            ! initialize ierr
            ierr = 0
            
            ! loop over the M elements of U_X and DU
            do jd = 1,X%n_mod
                ! set up helper variables
                g_frac = g13/g33
                Theta_3 = h23/h22
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(Theta_3(:,id,kd),D3Theta_3(:,id,kd),&
                            &grid%theta_F(:,id,kd),1,2)                         ! higher precision because other derivative will be taken later
                        CHCKERR('')
                    end do
                    if (.not.use_pol_flux_F) D3Theta_3(:,:,kd) = &
                        &D3Theta_3(:,:,kd)*eq%rot_t_FD(kd,0)                    ! parallel deriv. equal to iota * poloidal deriv.
                end do
                
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    ! set up correction to U
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*(g_frac(:,:,kd)*&
                        &(D3Theta_3(:,:,kd)+iu*n_frac*Theta_3(:,:,kd)))
                    ! calculate X%U_0 and X%DU_0
                    X%U_0(:,:,kd,jd) = &
                        &-(h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd))&
                        &+ iu/(mn(jd)*g33(:,:,kd)) * (g13(:,:,kd)*djq(kd) + &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*vac_perm + iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) - &
                        &g23(:,:,kd) )) + U_corr(:,:,kd,jd)
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_0(:,id,kd,jd),D3var(:,id),&
                            &grid%theta_F(:,id,kd),1,2)
                        CHCKERR('')
                    end do
                    if (.not.use_pol_flux_F) D3var = D3var*eq%rot_t_FD(kd,0)    ! parallel deriv. equal to iota * poloidal deriv.
                    X%DU_0(:,:,kd,jd) = D3var + iu*n_frac*X%U_0(:,:,kd,jd)
                    ! calculate X%U_1 and X%DU_1
                    X%U_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,:,kd)/g33(:,:,kd))
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_1(:,id,kd,jd),D3var(:,id),&
                            &grid%theta_F(:,id,kd),1,2)
                        CHCKERR('')
                    end do
                    if (.not.use_pol_flux_F) D3var = D3var*eq%rot_t_FD(kd,0)    ! parallel deriv. equal to iota * poloidal deriv.
                    X%DU_1(:,:,kd,jd) = D3var + iu*n_frac*X%U_1(:,:,kd,jd)
                end do
            end do
            
            ! deallocate
            deallocate(D3var)
        end function calc_U_HEL
        
#if ldebug
        ! test calculation of DU by deriving U numerically
        integer function test_DU() result(ierr)
            use utilities, only: calc_deriv
            use grid_ops, only: trim_grid
            use output_ops, only: plot_diff_HDF5
            
            character(*), parameter :: rout_name = 'test_DU'
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            complex(dp), allocatable :: DU_0(:,:,:)                             ! alternative calculation for DU_0
            complex(dp), allocatable :: DU_1(:,:,:)                             ! alternative calculation for DU_1
            type(grid_type) :: grid_trim                                        ! trimmed equilibrium grid
            integer :: tot_dim(3), grp_offset(3)                                ! total dimensions and group offset
            character(len=max_str_ln) :: file_name                              ! name of plot file
            character(len=max_str_ln) :: description                            ! description of plot
            
            ! initialize ierr
            ierr = 0
            
            ! warning
            call writo('WARNING: This test is only representable if there &
                &are enough parallel points in the grid!')
            
            ! output
            call writo('Going to test whether DU is consistent with U')
            call lvl_ud(1)
            
            ! set up DU_0 and DU_1
            allocate(DU_0(grid%n(1),grid%n(2),grid%grp_n_r))
            allocate(DU_1(grid%n(1),grid%n(2),grid%grp_n_r))
            
            ! trim extended grid into plot grid
            ierr = trim_grid(grid,grid_trim)
            CHCKERR('')
            
            ! set total and group dimensions and group offset
            tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
            grp_offset = [0,0,grid_trim%i_min-1]
            
            ! loop over all modes
            do jd = 1,X%n_mod
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    ! derive numerically
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_0(:,id,kd,jd),DU_0(:,id,kd),&
                            &ang_par_F(:,id,kd),1,2)
                        CHCKERR('')
                        ierr = calc_deriv(X%U_1(:,id,kd,jd),DU_1(:,id,kd),&
                            &ang_par_F(:,id,kd),1,2)
                        CHCKERR('')
                    end do
                    
                    ! add the second part
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    DU_0(:,:,kd) = DU_0(:,:,kd) + iu*n_frac*X%U_0(:,:,kd,jd)
                    DU_1(:,:,kd) = DU_1(:,:,kd) + iu*n_frac*X%U_1(:,:,kd,jd)
                end do
                
                ! set some variables
                file_name = 'TEST_RE_DU_0_'//trim(i2str(jd))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(jd)//', real part')
                
                ! plot difference for RE DU_0
                call plot_diff_HDF5(realpart(DU_0(:,:,1:grid_trim%grp_n_r)),&
                    &realpart(X%DU_0(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_offset,description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_0_'//trim(i2str(jd))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(jd)//', imaginary part')
                
                ! plot difference for IM DU_0
                call plot_diff_HDF5(imagpart(DU_0(:,:,1:grid_trim%grp_n_r)),&
                    &imagpart(X%DU_0(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_offset,description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_RE_DU_1_'//trim(i2str(jd))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(jd)//', real part')
                
                ! plot difference for RE DU_1
                call plot_diff_HDF5(realpart(DU_1(:,:,1:grid_trim%grp_n_r)),&
                    &realpart(X%DU_1(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_offset,description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_1_'//trim(i2str(jd))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(jd)//', imaginary part')
                
                ! plot difference for IM DU_1
                call plot_diff_HDF5(imagpart(DU_1(:,:,1:grid_trim%grp_n_r)),&
                    &imagpart(X%DU_1(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_offset,description,output_message=.true.)
            end do
            
            ! clean up
            call dealloc_grid(grid_trim)
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end function test_DU
#endif
    end function calc_U
    
    ! Calculate the  magnetic integrals  from PV_i and  KV_i. All  the variables
    ! should thus be field-line oriented.
    integer function calc_magn_ints(grid_eq,met,X) result(ierr)
        use grid_ops, only: calc_int_magn
        
        character(*), parameter :: rout_name = 'calc_magn_ints'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculating field-line averages')
        call lvl_ud(1)
        
        ! Calculate PV_int = <PV e^(k-m)ang_par_F>
        call writo('Taking field average of PV')
        call lvl_ud(1)
        ierr = calc_int_magn(grid_eq,met,X%exp_ang_par_F,X%n_mod,&
            &X%PV_0,X%PV_int_0)
        CHCKERR('')
        ierr = calc_int_magn(grid_eq,met,X%exp_ang_par_F,X%n_mod,&
            &X%PV_1,X%PV_int_1)
        CHCKERR('')
        ierr = calc_int_magn(grid_eq,met,X%exp_ang_par_F,X%n_mod,&
            &X%PV_2,X%PV_int_2)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! Calculate KV_int = <KV e^(k-m)ang_par_F>
        call writo('Taking field average of KV')
        call lvl_ud(1)
        ierr = calc_int_magn(grid_eq,met,X%exp_ang_par_F,X%n_mod,&
            &X%KV_0,X%KV_int_0)
        CHCKERR('')
        ierr = calc_int_magn(grid_eq,met,X%exp_ang_par_F,X%n_mod,&
            &X%KV_1,X%KV_int_1)
        CHCKERR('')
        ierr = calc_int_magn(grid_eq,met,X%exp_ang_par_F,X%n_mod,&
            &X%KV_2,X%KV_int_2)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
        call writo('Field-line averages calculated')
    end function calc_magn_ints
    
    ! Print perturbation quantities to an output file:
    !   - grid:     r_F, theta_F, zeta_F
    !   - X:        pres_FD, q_saf_FD, rot_t_FD, flux_p_FD, flux_t_FD, rho, S,
    !               kappa_n, kappa_g, sigma
    ! Note: The equilibrium quantities are outputted in Flux coordinates.
    integer function print_output_X(grid_eq,grid_X,X,X_B) result(ierr)
        use num_vars, only: rich_lvl_nr, max_it_r, grp_rank
        use HDF5_ops, only: print_HDF5_arrs, &
            &var_1D
        use grid_ops, only: trim_grid
        use X_vars, only: min_r_X, max_r_X, min_m_X, max_m_X, min_n_X, max_n_X
        
        character(*), parameter :: rout_name = 'print_output_X'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid variables
        type(X_type), intent(in) :: X                                           ! perturbation variables (for U_0, etc.)
        type(X_type), intent(in) :: X_B                                         ! field-aligned perturbation variables (for val, vec, etc)
        
        ! local variables
        type(var_1D), allocatable, target :: X_1D(:)                            ! 1D equivalent of eq. variables
        type(var_1D), pointer :: X_1D_loc => null()                             ! local element in X_1D
        type(grid_type) :: grid_X_trim                                          ! trimmed X grid
        type(grid_type) :: grid_eq_trim                                         ! trimmed eq grid
        integer :: i_min_X, i_max_X                                             ! min. and max. index of variables
        integer :: i_min_eq, i_max_eq                                           ! min. and max. index of variables
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Writing perturbation variables to output file')
        call lvl_ud(1)
        
        ! user output
        call writo('Preparing variables for writing')
        call lvl_ud(1)
        
        ! trim grids
        ierr = trim_grid(grid_X,grid_X_trim)
        CHCKERR('')
        ierr = trim_grid(grid_eq,grid_eq_trim)
        CHCKERR('')
        
        ! set i_min and i_max
        i_min_X = 1
        i_max_X = grid_X_trim%grp_n_r
        i_min_eq = 1
        i_max_eq = grid_eq_trim%grp_n_r
        
        ! Set up the 1D equivalents  of the perturbation variables
        allocate(X_1D(15))
        id = 1
        
        ! r_F
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'r_F'
        allocate(X_1D_loc%tot_i_min(1),X_1D_loc%tot_i_max(1))
        allocate(X_1D_loc%grp_i_min(1),X_1D_loc%grp_i_max(1))
        X_1D_loc%tot_i_min = [1]
        X_1D_loc%tot_i_max = [grid_X_trim%n(3)]
        X_1D_loc%grp_i_min = [grid_X_trim%i_min]
        X_1D_loc%grp_i_max = [grid_X_trim%i_max]
        allocate(X_1D_loc%p(size(grid_X_trim%grp_r_F(i_min_X:i_max_X))))
        X_1D_loc%p = grid_X_trim%grp_r_F(i_min_X:i_max_X)
        
        ! r_E
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'r_E'
        allocate(X_1D_loc%tot_i_min(1),X_1D_loc%tot_i_max(1))
        allocate(X_1D_loc%grp_i_min(1),X_1D_loc%grp_i_max(1))
        X_1D_loc%tot_i_min = [1]
        X_1D_loc%tot_i_max = [grid_X_trim%n(3)]
        X_1D_loc%grp_i_min = [grid_X_trim%i_min]
        X_1D_loc%grp_i_max = [grid_X_trim%i_max]
        allocate(X_1D_loc%p(size(grid_X_trim%grp_r_E(i_min_X:i_max_X))))
        X_1D_loc%p = grid_X_trim%grp_r_E(i_min_X:i_max_X)
        
        ! RE_U_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_U_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%U_0(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(realpart(X%U_0(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%U_0(:,:,i_min_eq:i_max_eq,:))])
        
        ! IM_U_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_U_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%U_0(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(imagpart(X%U_0(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%U_0(:,:,i_min_eq:i_max_eq,:))])
        
        ! RE_U_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_U_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%U_1(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(realpart(X%U_1(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%U_1(:,:,i_min_eq:i_max_eq,:))])
        
        ! IM_U_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_U_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%U_1(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(imagpart(X%U_1(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%U_1(:,:,i_min_eq:i_max_eq,:))])
        
        ! RE_DU_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_DU_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%DU_0(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(realpart(X%DU_0(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%DU_0(:,:,i_min_eq:i_max_eq,:))])
        
        ! IM_DU_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_DU_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%DU_0(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(imagpart(X%DU_0(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%DU_0(:,:,i_min_eq:i_max_eq,:))])
        
        ! RE_DU_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_DU_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%DU_1(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(realpart(X%DU_1(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%DU_1(:,:,i_min_eq:i_max_eq,:))])
        
        ! IM_DU_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_DU_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%grp_i_min(4),X_1D_loc%grp_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [grid_eq_trim%n,X%n_mod]
        X_1D_loc%grp_i_min = [1,1,grid_eq_trim%i_min,1]
        X_1D_loc%grp_i_max = &
            &[grid_eq_trim%n(1:2),grid_eq_trim%i_max,X%n_mod]
        allocate(X_1D_loc%p(size(X%DU_1(:,:,i_min_eq:i_max_eq,:))))
        X_1D_loc%p = reshape(imagpart(X%DU_1(:,:,i_min_eq:i_max_eq,:)),&
            &[size(X%DU_1(:,:,i_min_eq:i_max_eq,:))])
        
        ! RE_X_val
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_X_val'
        allocate(X_1D_loc%tot_i_min(1),X_1D_loc%tot_i_max(1))
        allocate(X_1D_loc%grp_i_min(1),X_1D_loc%grp_i_max(1))
        X_1D_loc%tot_i_min = [1]
        X_1D_loc%tot_i_max = [size(X_B%val)]
        X_1D_loc%grp_i_min = [1]
        if (grp_rank.eq.0) then
            X_1D_loc%grp_i_max = [size(X_B%val)]
            allocate(X_1D_loc%p(size(X_B%val)))
            X_1D_loc%p = realpart(X_B%val)
        else
            X_1D_loc%grp_i_max = [0]
            allocate(X_1D_loc%p(0))
        end if
        
        ! IM_X_val
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_X_val'
        allocate(X_1D_loc%tot_i_min(1),X_1D_loc%tot_i_max(1))
        allocate(X_1D_loc%grp_i_min(1),X_1D_loc%grp_i_max(1))
        X_1D_loc%tot_i_min = [1]
        X_1D_loc%tot_i_max = [size(X_B%val)]
        X_1D_loc%grp_i_min = [1]
        if (grp_rank.eq.0) then
            X_1D_loc%grp_i_max = [size(X_B%val)]
            allocate(X_1D_loc%p(size(X_B%val)))
            X_1D_loc%p = imagpart(X_B%val)
        else
            X_1D_loc%grp_i_max = [0]
            allocate(X_1D_loc%p(0))
        end if
        
        ! RE_X_vec
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_X_vec'
        allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
        allocate(X_1D_loc%grp_i_min(3),X_1D_loc%grp_i_max(3))
        X_1D_loc%grp_i_min = [1,grid_X_trim%i_min,1]
        X_1D_loc%grp_i_max = [X_B%n_mod,grid_X_trim%i_max,size(X_B%vec,3)]
        X_1D_loc%tot_i_min = [1,1,1]
        X_1D_loc%tot_i_max = [X_B%n_mod,grid_X_trim%n(3),size(X_B%vec,3)]
        allocate(X_1D_loc%p(size(X_B%vec(:,i_min_X:i_max_X,:))))
        X_1D_loc%p = reshape(realpart(X_B%vec(:,i_min_X:i_max_X,:)),&
            &[size(X_B%vec(:,i_min_X:i_max_X,:))])
        
        ! IM_X_vec
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_X_vec'
        allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
        allocate(X_1D_loc%grp_i_min(3),X_1D_loc%grp_i_max(3))
        X_1D_loc%grp_i_min = [1,grid_X_trim%i_min,1]
        X_1D_loc%grp_i_max = [X_B%n_mod,grid_X_trim%i_max,size(X_B%vec,3)]
        X_1D_loc%tot_i_min = [1,1,1]
        X_1D_loc%tot_i_max = [X_B%n_mod,grid_X_trim%n(3),size(X_B%vec,3)]
        allocate(X_1D_loc%p(size(X_B%vec(:,i_min_X:i_max_X,:))))
        X_1D_loc%p = reshape(imagpart(X_B%vec(:,i_min_X:i_max_X,:)),&
            &[size(X_B%vec(:,i_min_X:i_max_X,:))])
        
        ! misc_X
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'misc_X'
        allocate(X_1D_loc%tot_i_min(1),X_1D_loc%tot_i_max(1))
        allocate(X_1D_loc%grp_i_min(1),X_1D_loc%grp_i_max(1))
        if (grp_rank.eq.0) then
            X_1D_loc%grp_i_min = [1]
            X_1D_loc%grp_i_max = [6]
            allocate(X_1D_loc%p(6))
            X_1D_loc%p = [min_r_X,max_r_X,min_n_X*1._dp,max_n_X*1._dp,&
                &min_m_X*1._dp,max_m_X*1._dp]
        else
            X_1D_loc%grp_i_min = [1]
            X_1D_loc%grp_i_max = [0]
            allocate(X_1D_loc%p(0))
        end if
        X_1D_loc%tot_i_min = [1]
        X_1D_loc%tot_i_max = [6]
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Writing using HDF5')
        call lvl_ud(1)
        
        ! write
        if (max_it_r.gt.1) then
            ierr = print_HDF5_arrs(X_1D,'X_R'//trim(i2str(rich_lvl_nr)))
        else
            ierr = print_HDF5_arrs(X_1D,'X')
        end if
        CHCKERR('')
        
        ! clean up
        call dealloc_grid(grid_X_trim)
        call dealloc_grid(grid_eq_trim)
        nullify(X_1D_loc)
        
        ! user output
        call lvl_ud(-1)
        call writo('Writing complete')
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation variables written to output')
    end function print_output_X
end module
