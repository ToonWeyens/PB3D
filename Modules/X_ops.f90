!------------------------------------------------------------------------------!
!   Operations considering perturbation quantities                             !
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    use grid_vars, onlY: grid_type
    use eq_vars, only: eq_type
    use metric_vars, only: metric_type
    use X_vars, only: X_type

    implicit none
    private
    public prepare_X, solve_EV_system, calc_PV, calc_KV, calc_U, calc_extra, &
        &calc_magn_ints

contains
    ! prepare the matrix elements by calculating  KV_i and PV_i, which then will
    ! have to be integrated, with a complex exponential weighting function
    integer function prepare_X(grid,eq,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, plot_jq, grp_nr
        use X_vars, only: dealloc_X, X_type
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'prepare_X'
        
        ! input / output
        type(grid_type) :: grid                                                 ! equilibrium grid variables
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(metric_type) :: met                                                ! metric variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables
        integer :: m, k                                                         ! counters
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle in flux coordinates
        character(len=10) :: integrand_name                                     ! name of integrand in field-line averages
        
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
        use num_vars, only: use_pol_flux_F, output_style, grp_rank, tol_NR, &
            &max_it_NR, no_plots
        use utilities, only: calc_zero_NR, interp_fun
        use grid_vars, only: destroy_grid
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        use grid_ops, only: trim_grid
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(grid_type) :: grid                                                 ! grid variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables (also used in child functions)
        real(dp) :: mnfrac_fun                                                  ! fraction m/n or n/m to determine resonant flux surface
        
        ! local variables (not to be used in child functions)
        integer :: jd, kd                                                       ! counters
        real(dp), allocatable :: x_vars(:,:)                                    ! for plotting
        real(dp), allocatable :: y_vars(:,:)                                    ! for plotting
        real(dp) :: jq_sol                                                      ! solution q or iota in F coords.
        integer :: istat                                                        ! status
        character(len=max_str_ln) :: plot_title, file_name                      ! name of plot, of file
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: jq(:,:)                                        ! saf. fac. or rot. transf. in Flxu coords.
        real(dp), allocatable :: jq_loc(:)                                      ! local version of jq
        integer :: n_r                                                          ! total number of normal points
        real(dp) :: tol_NR_old                                                  ! old value of tol_NR
        integer :: max_it_NR_old                                                ! old value of max_it_NR
        type(grid_type) :: grid_trim                                            ! trimmed version of grid
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! get trimmed grid
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
        ! initialize variables
        n_r = grid_trim%n(3)
        
        ! Get serial version  of safety factor or rot. transform  and print user
        ! message.  As  the normal  coordinate  is  rescaled  to 0..1,  the  ith
        ! derivative has to be rescaled by the factor (max_flux_F/2pi)^i.
        if (use_pol_flux_F) then
            call writo('Plotting safety factor q and resonant surfaces &
                &q = m/n')
            if (grp_rank.eq.0) allocate(jq(n_r,0:2))
            do jd = 0,2
                ierr = get_ser_var(eq%q_saf_FD(1:grid_trim%grp_n_r,jd),jq_loc)
                CHCKERR('')
                if(grp_rank.eq.0) jq(:,jd) = jq_loc*(max_flux_p_F/(2*pi))**jd
            end do
        else
            call writo('Plotting rotational transform iota and resonant &
                &surfaces iota = n/m')
            if (grp_rank.eq.0) allocate(jq(n_r,0:2))
            do jd = 0,2
                ierr = get_ser_var(eq%rot_t_FD(1:grid_trim%grp_n_r,jd),jq_loc)
                CHCKERR('')
                if(grp_rank.eq.0) jq(:,jd) = jq_loc*(max_flux_t_F/(2*pi))**jd
            end do
        end if
        
        ! the rest is done only by global master
        if (grp_rank.eq.0) then
            call writo('calculating resonant surfaces')
            
            ! initialize variables
            allocate(x_vars(n_r,X%n_mod+1)); x_vars = 0
            allocate(y_vars(n_r,X%n_mod+1)); y_vars = 0
            
            ! set x_vars and y_vars for first column
            x_vars(:,1) = grid_trim%r_F/grid_trim%r_F(n_r)
            y_vars(:,1) = jq(:,0)
            
            ! save old tol_NR and max_it_NR
            tol_NR_old = tol_NR
            max_it_NR_old = max_it_NR
            
            ! change tol_NR and max_it_NR
            tol_NR = 1.E-8_dp
            max_it_NR = 5000
            
            ! loop over all modes (and shift the index in x and y_vars by 1)
            kd = 2
            do jd = 1, X%n_mod
                ! find place where q = m/n or  iota = n/m in Flux coordinates by
                ! solving q-m/n = 0 or iota-n/m=0, using the functin jq_fun
                call lvl_ud(1)
                
                ! set up mnfrac for function
                if (use_pol_flux_F) then
                    mnfrac_fun = 1.0_dp*X%m(jd)/X%n(jd)
                else                             
                    mnfrac_fun = 1.0_dp*X%n(jd)/X%m(jd)
                end if
                
                ! calculate zero using Newton-Rhapson
                istat = calc_zero_NR(jq_sol,jq_fun,jq_dfun,1.0_dp)
                
                ! intercept error
                if (istat.ne.0) then
                    call writo('Error intercepted: Couldn''t find resonating &
                        &surface for (n,m) = ('//trim(i2str(X%n(jd)))//','//&
                        &trim(i2str(X%m(jd)))//')')
                else if (jq_sol.lt.0.0_dp) then
                    call writo('Mode (n,m) = ('//trim(i2str(X%n(jd)))//','//&
                        &trim(i2str(X%m(jd)))//') does not resonate in plasma')
                else
                    if (jq_sol.gt.1.0_dp) then
                        call writo('Mode (n,m) = ('//trim(i2str(X%n(jd)))//','&
                            &//trim(i2str(X%m(jd)))//') does not resonate &
                            &in plasma')
                        y_vars(n_r,kd) = jq(n_r,0)
                    else
                        ! set x and y vars
                        x_vars(:,kd) = jq_sol
                        y_vars(n_r,kd) = mnfrac_fun
                        ! GNUplot output
                        if (output_style.eq.1) call writo('Mode (n,m) = ('//&
                            &trim(i2str(X%n(jd)))//','//trim(i2str(X%m(jd)))//&
                            &') resonates in plasma at normalized flux &
                            &surface '//trim(r2str(jq_sol)))
                    end if
                    kd = kd + 1
                end if
                
                call lvl_ud(-1)
            end do
            
            ! recover old tol_NR and max_it_NR
            tol_NR = tol_NR_old
            max_it_NR = max_it_NR_old
            
            ! user message
            call writo('Plotting results')
            if (use_pol_flux_F) then
                plot_title = 'safety factor q'
                file_name = 'q_saf'
            else
                plot_title = 'rotational transform iota'
                file_name = 'rot_t'
            end if
            
            ! plot according to output_style
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
            deallocate(jq)
            call destroy_grid(grid_trim)
        
        end if
        
        call lvl_ud(-1)
        call writo('Done plotting')
    contains
        ! Returns q-m/n or  iota-n/m in Flux coordinates, used to  solve for q =
        ! m/n or iota = n/m.
        ! WARNING: This routine requires that  jq's derivatives be calculated up
        ! to order 1. This is NOT checked!
        real(dp) function jq_fun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.x_vars(1,1)) then                                         ! point requested lower than minimum x
                ! extrapolate variable from 1
                res = jq(1,0) - mnfrac_fun + jq(1,1)*(pt-x_vars(1,1))
            else if (pt.gt.x_vars(n_r,1)) then                                  ! point requested higher than maximum x
                ! extrapolate variable from n_r
                res = jq(n_r,0) - mnfrac_fun + jq(n_r,1)*(pt-x_vars(n_r,1))
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                istat = interp_fun(res,jq(:,0)-mnfrac_fun,pt,x_vars(:,1))
            end if
        end function jq_fun
        
        ! returns d(q-m/n)/dr in Flux coordinates, used to solve for q = m/n
        ! WARNING: This routine requires that  jq's derivatives be calculated up
        ! to order 2. This is NOT checked!
        real(dp) function jq_dfun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.x_vars(1,1)) then                                         ! point requested lower than minimum x
                ! extrapolate variable from 1
                res = jq(1,1) + jq(1,2)*(pt-x_vars(1,1))
            else if (pt.gt.x_vars(n_r,1)) then                                  ! point requested higher than maximum x
                ! extrapolate variable from n_r
                res = jq(n_r,1) + jq(n_r,2)*(pt-x_vars(n_r,1))
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                istat = interp_fun(res,jq(:,1),pt,x=x_vars(:,1))
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
            real(dp), allocatable :: r_plot_E(:)                                ! normal E coordinates of plot
            real(dp), allocatable :: r_plot_F(:)                                ! normal F coordinates of plot
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
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,1))
            ierr = calc_eqd_grid(theta_plot,1*pi,3*pi,1)                        ! starting from pi gives nicer results
            CHCKERR('')
            ierr = calc_eqd_grid(zeta_plot,0*pi,2*pi,2)
            CHCKERR('')
            
            ! set up vars
            allocate(vars(n_theta_plot,n_zeta_plot,1,X%n_mod))
            do id = 1,X%n_mod
                vars(:,:,:,id) = y_vars(n_r,id+1)
            end do
            
            ! set dimensions
            plot_dim = [n_theta_plot,n_zeta_plot,1,X%n_mod]
            
            ! setup normal vars in F coords.
            allocate(r_plot_F(X%n_mod))
            if (use_pol_flux_F) then
                r_plot_F = x_vars(n_r,2:X%n_mod+1)*max_flux_p_F/(2*pi)
            else
                r_plot_F = x_vars(n_r,2:X%n_mod+1)*max_flux_t_F/(2*pi)
            end if
            
            ! calculate normal vars in Equilibrium coords.
            allocate(r_plot_E(X%n_mod))
            ierr = coord_F2E(grid,eq,r_plot_F,r_plot_E)
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
                grid_plot%grp_r_E = r_plot_E(id)
                
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
            deallocate(theta_plot,zeta_plot,r_plot_E,r_plot_F)
            
            ! delete plot grid
            call destroy_grid(grid_plot)
        end function resonance_plot_HDF5
    end function resonance_plot
    
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
        use num_vars, only: glb_rank, use_pol_flux_F, eq_style
        
        character(*), parameter :: rout_name = 'check_modes'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(X_type) :: X                                                       ! perturbation variables
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp) :: tol = 0.1_dp                                                ! tolerance for being out of range of q or iota values
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
            call writo('The tolerance used is '//trim(r2strt(tol)))
            
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
            
            ! initialize min_sec_X and max_sec_X
            min_sec_X = huge(1._dp)
            max_sec_X = -huge(1._dp)
            
            ! for every mode (n,m) check whether  m/n is inside the range of
            ! q values or n/m inside the range of iota values
            do id = 1, X%n_mod
                ! calculate upper and lower limits
                lim_lo = max(min_jq-tol,min_jq/(1+pmone2*tol))
                lim_hi = min(max_jq+tol,max_jq/(1-pmone2*tol))
                
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
                &trim(r2strt(min_sec_X))//' < '//mode_name//' < '//&
                &trim(r2strt(max_sec_X)))
            
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
    
    ! Calculate rho according to the style specified by user input.
    integer function calc_rho(eq,grid) result(ierr)
        use num_vars, only: rho_style, use_normalization
        use eq_vars, only: rho_0
        
        character(*), parameter :: rout_name = 'calc_rho'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(grid_type) :: grid                                                 ! grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! allocate rho
        allocate(eq%rho(grid%grp_n_r))
        
        ! choose which density style is being used:
        !   1:  constant, equal to rho_0
        select case (rho_style)
            case (1)                                                            ! HELENA
                eq%rho = rho_0                                                  ! arbitrarily constant (normalized value)
            case default
                err_msg = 'No density style associated with '//&
                    &trim(i2str(rho_style))
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
        real(dp), pointer :: g33(:,:,:)                                         ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        
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
            end do
        end do
        
        ! PV_1 is not Hermitian
        do m = 1,X%n_mod
            do k = 1,X%n_mod
                ! calculate PV_1
                X%PV_1(:,:,:,c([k,m],.false.,X%n_mod)) = &
                    &com_fac * X%DU_1(:,:,:,m) * &
                    &(conjg(X%DU_0(:,:,:,k)) - X%extra1 - X%extra2)
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
    end subroutine calc_KV
    
    ! calculate U_m^0, U_m^1 or U_n^0, U_n^1  at eq grp_n_r values of the normal
    ! coordinate, n_par values  of the parallel coordinate and  size_X values of
    ! the poloidal mode number,  or of the toroidal mode number,  as well as the
    ! poloidal derivatives
    ! Note: The normal derivatives have the  factor i/n or i/m included already,
    ! as opposed to [ADD REF]
    integer function calc_U(eq,grid,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, mu_0, use_normalization, eq_style, &
            &ltest
        use utilities, only: c
        use input_ops, only: get_log, pause_prog
        
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
        
        ! deallocate
        deallocate(djq,mn)
        deallocate(fac_n,fac_m)
        deallocate(U_corr,D3U_corr)
        nullify(J,D3J)
        nullify(g13,D3g13,g23,D3g23,g33,D3g33)
        nullify(h12,D3h12,h22,D1h22,D3h22,h23,D1h23,D3h23)
        nullify(g_frac)
        nullify(T_theta,D1T_theta,D3T_theta)
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
                    X%DU_1(:,:,kd,jd) = iu/mn(jd) * n_frac/mn(jd) * &
                        &( D3g13(:,:,kd)/g33(:,:,kd) - &
                        &  g13(:,:,kd)*D3g33(:,:,kd)/g33(:,:,kd)**2 ) + &
                        &iu*n_frac*X%U_1(:,:,kd,jd) 
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
            complex(dp), allocatable :: D3_var(:,:)                             ! derivative of variable
            
            ! allocate extra helper variables
            allocate(D3_var(grid%n(1),grid%n(2)))
            
            ! initialize ierr
            ierr = 0
            
            ! loop over the M elements of U_X and DU
            do jd = 1,X%n_mod
            write(*,*) '!!! THIS IS WRONG: YOU HAVE TO TAKE T_ZETA, NOT T_THETA !!!!'
            write(*,*) 'I CHANGED theta_F to ang_par_F, THIS CAUSED INFINITIES!!!'
                ! set up helper variables
                g_frac = g13/g33
                T_theta = h23/h22
                ! loop over all normal points
                do kd = 1,grid%grp_n_r
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(T_theta(:,id,kd),D3T_theta(:,id,kd),&
                            &ang_par_F(:,id,kd),1,2)                            ! higher precision because other derivative will be taken later
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
                        ierr = calc_deriv(X%U_0(:,id,kd,jd),D3_var(:,id),&
                            &ang_par_F(:,id,kd),1,2)
                    end do
                    CHCKERR('')
                    X%DU_0(:,:,kd,jd) = D3_var + iu*n_frac*X%U_0(:,:,kd,jd)
                    ! calculate X%U_1 and X%DU_1
                    X%U_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,:,kd)/g33(:,:,kd))
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_1(:,id,kd,jd),D3_var(:,id),&
                            &ang_par_F(:,id,kd),1,2)
                    end do
                    CHCKERR('')
                    X%DU_1(:,:,kd,jd) = D3_var + iu*n_frac*X%U_1(:,:,kd,jd)
                end do
            end do
            
            ! deallocate
            deallocate(D3_var)
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
            integer :: tot_dim(3), grp_dim(3), grp_offset(3)                    ! total and group dimensions and group offset
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
            grp_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r]
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
                    &tot_dim,grp_dim,grp_offset,description,&
                    &output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_0_'//trim(i2str(jd))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(jd)//', imaginary part')
                
                ! plot difference for IM DU_0
                call plot_diff_HDF5(imagpart(DU_0(:,:,1:grid_trim%grp_n_r)),&
                    &imagpart(X%DU_0(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_dim,grp_offset,description,&
                    &output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_RE_DU_1_'//trim(i2str(jd))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(jd)//', real part')
                
                ! plot difference for RE DU_1
                call plot_diff_HDF5(realpart(DU_1(:,:,1:grid_trim%grp_n_r)),&
                    &realpart(X%DU_1(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_dim,grp_offset,description,&
                    &output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_1_'//trim(i2str(jd))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(jd)//', imaginary part')
                
                ! plot difference for IM DU_1
                call plot_diff_HDF5(imagpart(DU_1(:,:,1:grid_trim%grp_n_r)),&
                    &imagpart(X%DU_1(:,:,1:grid_trim%grp_n_r,jd)),file_name,&
                    &tot_dim,grp_dim,grp_offset,description,&
                    &output_message=.true.)
            end do
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end function test_DU
#endif
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
    end function calc_extra
    
    ! Calculate the  magnetic integrals  from PV_i and  KV_i. All  the variables
    ! should thus be field-line oriented.
    integer function calc_magn_ints(grid_eq,met,X) result(ierr)
        character(*), parameter :: rout_name = 'calc_magn_ints'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
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
    contains
        ! calculates magnetic integral of V, defined as the matrix
        !   <V e^[i(k-m)theta_F]> = [ oint J V(k,m) e^i(k-m)theta_F dtheta_F ],
        ! or
        !   <V e^[i(m-k)zeta_F]> = [ oint J V(k,m) e^i(m-k)zeta_F dzeta_F ],
        ! depending on whether pol. or tor. flux is used as normal coord.
        integer function calc_int_magn(grid,met,exp_ang,n_mod,V,V_int) result(ierr)
            use num_vars, only: use_pol_flux_F
            use utilities, only: calc_mult, c, con, is_sym
            
            character(*), parameter :: rout_name = 'calc_int_magn'
            
            ! input / output
            type(grid_type) :: grid                                             ! grid
            type(metric_type) :: met                                            ! metric variables
            complex(dp), intent(in) :: exp_ang(:,:,:,:)                         ! exponential of Flux parallel angle
            integer, intent(in) :: n_mod                                        ! number of 
            complex(dp), intent(in) :: V(:,:,:,:)                               ! input V(n_par,n_geo,n_r,size_X^2)
            complex(dp), intent(inout) :: V_int(:,:,:)                          ! output <V e^i(k-m)ang_par_F> integrated in parallel Flux coord.
            
            ! local variables
            integer :: k, m, id, jd, kd                                         ! counters
            integer :: nn_mod                                                   ! number of indices for V and V_int
            integer :: k_min                                                    ! minimum k
            logical :: sym                                                      ! whether V and V_int are symmetric
            complex(dp), allocatable :: V_J_e(:,:,:,:)                          ! V*J*exp_ang
            character(len=max_str_ln) :: err_msg                                ! error message
            real(dp), pointer :: ang_par_F(:,:,:)                               ! parallel angle
            integer :: dims(3)                                                  ! real dimensions
            
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
                    V_J_e(:,:,:,c([k,m],sym,n_mod)) = met%jac_FD(:,:,:,0,0,0) &
                        &* con(exp_ang(:,:,:,c([k,m],.true.,n_mod)),&
                        &[k,m],.true.,dims) * &
                        &con(V(:,:,:,c([k,m],sym,n_mod)),[k,m],sym,dims)
                end do
            end do
            
            ! integrate  term over  ang_par_F  for all  equilibrium grid  points
            ! using the recursive formula int_1^n f(x) dx
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
            
            ! deallocate local variables
            deallocate(V_J_e)
            nullify(ang_par_F)
        end function calc_int_magn
    end function calc_magn_ints
end module
