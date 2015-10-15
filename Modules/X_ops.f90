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
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type

    implicit none
    private
    public calc_X, solve_EV_system, calc_magn_ints, print_output_X, &
        &resonance_plot, calc_res_surf, check_X_modes
    
    ! interfaces
    interface calc_X
        module procedure calc_X_1, calc_X_2
    end interface
    interface print_output_X
        module procedure print_output_X_1, print_output_X_2
    end interface

contains
    ! Calculates either vectorial or tensorial perturbation variables.
    ! Optionally, the  secondary mode  numbers can be  specified (m  if poloidal
    ! flux is used and n if toroidal  flux). By default, they are taken from the
    ! global X_vars variables.
    integer function calc_X_1(grid_eq,eq,met,X,lim_sec_X) result(ierr)          ! vectorial version
        use X_vars, only: create_X
        
        character(*), parameter :: rout_name = 'calc_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculating vectorial perturbation variables')
        call lvl_ud(1)
        
        ! create perturbation with modes of current X job
        call create_X(grid_eq,X,lim_sec_X)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        ierr = calc_U(eq,grid_eq,met,X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
        call writo('Vectorial perturbation variables calculated')
    end function calc_X_1
    integer function calc_X_2(grid_eq,eq,met,X_a,X_b,X,lim_sec_X) result(ierr)  ! tensorial version
        use X_vars, only: create_X
        
        character(*), parameter :: rout_name = 'calc_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(inout) :: X_a, X_b                               ! vectorial perturbation variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculating tensorial perturbation variables')
        call lvl_ud(1)
        
        ! create perturbation with modes of current X job
        call create_X(grid_eq,X,lim_sec_X)
        
        ! Calculate  PV_i  for  all (k,m)  pairs  and n_r  (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating PV...')
        call lvl_ud(1)
        ierr = calc_PV(eq,grid_eq,met,X_a,X_b,X,lim_sec_X)
        CHCKERR('')
        call lvl_ud(-1)
        
        !  Calculate KV_i  for  all (k,m)  pairs  and n_r  (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV...')
        call lvl_ud(1)
        ierr = calc_KV(eq,grid_eq,met,X_a,X_b,X,lim_sec_X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
        call writo('Tensorial perturbation variables calculated')
    end function calc_X_2
    
    ! plot  q-profile  or iota-profile  in  flux coordinates  with nq-m  = 0  or
    ! n-iotam = 0 indicate if requested
    integer function resonance_plot(eq,grid) result(ierr)
        use num_vars, only: use_pol_flux_F, no_plots, n_theta_plot, &
            &n_zeta_plot, rank
        use utilities, only: calc_zero_NR, interp_fun
        use grid_vars, only: create_grid, dealloc_grid
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        use grid_ops, only: calc_XYZ_grid, calc_eqd_grid, coord_F2E, trim_grid
        use MPI_utilities, only: get_ser_var
        use X_vars, only: min_n_X, min_m_X
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        
        ! local variables (not to be used in child functions)
        integer :: kd                                                           ! counter
        integer :: n_mod_loc                                                    ! local n_mod
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        real(dp), allocatable :: x_vars(:,:)                                    ! for plotting
        real(dp), allocatable :: y_vars(:,:)                                    ! for plotting
        character(len=max_str_ln) :: plot_title, file_name                      ! name of plot, of file
        real(dp), allocatable :: jq(:)                                          ! saf. fac. or rot. transf. in Flux coords.
        integer :: n_r                                                          ! total number of normal points
        integer :: plot_dim(4)                                                  ! plot dimensions (total = local because only masters)
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
        ierr = calc_res_surf(grid,eq,res_surf,info=.true.,jq=jq,&
            &tol_NR=1.E-8_dp,max_it_NR=5000)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! only master
        if (rank.eq.0) then
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
                        &trim(i2str(nint(res_surf(kd,3)*min_n_X)))//','//&
                        &trim(i2str(min_n_X))
                end do
            else                                                                ! m is fixed and n = n/m m
                do kd = 1,n_mod_loc
                    plot_titles(kd) = trim(plot_title)//' for m,n = '//&
                        &trim(i2str(min_m_X))//','//&
                        &trim(i2str(nint(res_surf(kd,3)*min_m_X)))
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
                ! set loc_r_E of plot grid
                grid_plot%loc_r_E = r_plot_E(kd)
                
                ! calculate X, Y and Z
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
    ! rotational transform can be returned to the master.
    ! Note: Has to be called by all processes
    integer function calc_res_surf(grid,eq,res_surf,info,jq,tol_NR,max_it_NR) &
        &result(ierr)
        use X_vars, only: min_n_X, max_n_X, min_m_X, max_m_X
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
        real(dp), intent(inout), allocatable :: res_surf(:,:)                   ! resonant surface
        logical, intent(in), optional :: info                                   ! if info is displayed
        real(dp), intent(inout), optional, allocatable :: jq(:)                 ! either safety factor or rotational transform
        real(dp), intent(in), optional :: tol_NR                                ! optional tolerance for Newton-Rhapson iteration
        integer, intent(in), optional :: max_it_NR                              ! optional max. nr. of iterations for Newton-Rhapson iteration
        
        ! local variables (also used in child functions)
        real(dp) :: mnfrac_fun                                                  ! fraction m/n or n/m to determine resonant flux surface
        
        ! local variables (not to be used in child functions)
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: kd                                                           ! counter
        integer :: kd_loc                                                       ! local kd
        integer :: n_mod                                                        ! nr. of modes
        logical :: info_loc                                                     ! local version of info
        integer :: istat                                                        ! status
        real(dp), allocatable :: res_surf_loc(:,:)                              ! local copy of res_surf
        real(dp), allocatable :: jq_tot(:,:)                                    ! saf. fac. or rot. transf. and derivs. in Flux coords.
        real(dp), allocatable :: jq_loc(:)                                      ! local version of jq
        real(dp) :: norm_factor                                                 ! normalization factor to make normal coordinate 0..1
        type(grid_type) :: grid_trim                                            ! trimmed version of grid
        real(dp) :: tol_NR_old                                                  ! old value of tol_NR
        integer :: max_it_NR_old                                                ! old value of max_it_NR
        integer :: n_X_loc, m_X_loc                                             ! local n_X and m_X
        
        ! initialize ierr
        ierr = 0
        
        ! set local info
        info_loc = .false.
        if (present(info)) info_loc = info
        
        ! get trimmed grid
        ierr = trim_grid(grid,grid_trim,norm_id)
        CHCKERR('')
        
        ! set n_mod
        n_mod = (max_n_X-min_n_X+1)*(max_m_X-min_m_X+1)
        
        ! get serial version of safety factor or rot. transform
        allocate(jq_tot(grid_trim%n(3),0:2))
        if (use_pol_flux_F) then
            if (grid%divided) then
                do kd = 0,2
                    ierr = get_ser_var(eq%q_saf_FD(norm_id(1):norm_id(2),kd),&
                        &jq_loc,scatter=.true.)
                    CHCKERR('')
                    jq_tot(:,kd) = jq_loc
                end do
            else
                jq_tot = eq%q_saf_FD
            end if
        else
            if (grid%divided) then
                do kd = 0,2
                    ierr = get_ser_var(eq%rot_t_FD(norm_id(1):norm_id(2),kd),&
                        &jq_loc,scatter=.true.)
                    CHCKERR('')
                    jq_tot(:,kd) = jq_loc
                end do
            else
                jq_tot = eq%q_saf_FD
            end if
        end if
        if (present(jq)) then
            allocate(jq(grid_trim%n(3)))
            jq = jq_tot(:,0)
        end if
        
        ! initialize local res_surf
        allocate(res_surf_loc(1:n_mod,3))
        
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
        do kd = 1,n_mod
            ! find place  where q  = m/n or  iota = n/m  in Flux  coordinates by
            ! solving q-m/n = 0 or iota-n/m=0, using the functin jq_fun
            
            ! set up n_X_loc and m_X_loc mnfrac for function
            if (use_pol_flux_F) then
                n_X_loc = min_n_X
                m_X_loc = min_m_X+kd-1
                mnfrac_fun = 1.0_dp*m_X_loc/n_X_loc
            else                             
                m_X_loc = min_m_X
                n_X_loc = min_n_X+kd-1
                mnfrac_fun = 1.0_dp*n_X_loc/m_X_loc
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
                    &trim(i2str(n_X_loc))//','//trim(i2str(m_X_loc))//')')
                call lvl_ud(-1)
            else if (res_surf_loc(kd_loc,2).lt.minval(grid_trim%r_F) .or. &
                &res_surf_loc(kd_loc,2).gt.maxval(grid_trim%r_F)) then
                if (info_loc) call writo('Mode (n,m) = ('//&
                    &trim(i2str(n_X_loc))//','//trim(i2str(m_X_loc))//&
                    &') does not resonate in plasma')
            else
                if (info_loc) call writo('Mode (n,m) = ('//&
                    &trim(i2str(n_X_loc))//','//trim(i2str(m_X_loc))//&
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
    integer function check_X_modes(eq) result(ierr)
        use MPI_utilities, only: get_ser_var
        use X_vars, only: min_n_X, max_n_X, min_m_X, max_m_X
        use num_vars, only: use_pol_flux_F, tol_norm_r
        
        character(*), parameter :: rout_name = 'check_X_modes'
        
        ! input / output
        type(eq_type) :: eq                                                     ! equilibrium variables
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_mod                                                        ! nr. of modes
        real(dp) :: min_jq, max_jq                                              ! min. and max. values of q or iota
        integer :: pmone                                                        ! plus or minus one
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: jq_name                                    ! either safety factor or rotational transform
        character(len=1) :: mode_name                                           ! either n or m
        real(dp) :: lim_lo, lim_hi                                              ! lower and upper limit on n/m (or m/n)
        real(dp) :: min_sec_X, max_sec_X                                        ! min. and max. of m_X (if pol. flux) or n_X (if tor. flux)
        
        ! initialize ierr
        ierr = 0
        
        call writo('Checking mode numbers')
        call lvl_ud(1)
        
        ! setup n_mod
        n_mod = (max_n_X-min_n_X+1)*(max_m_X-min_m_X+1)
        
        ! user output
        call writo('The tolerance used is '//trim(r2strt(tol_norm_r))//&
            &'...')
        
        ! set min_jq and max_jq in Flux coordinate system
        if (use_pol_flux_F) then 
            min_jq = minval(eq%q_saf_FD(:,0))
            max_jq = maxval(eq%q_saf_FD(:,0))
        else
            min_jq = minval(eq%rot_t_FD(:,0))
            max_jq = maxval(eq%rot_t_FD(:,0))
        end if
        
        ! set up jq name
        if (use_pol_flux_F) then
            jq_name = 'safety factor'
            mode_name = 'm'
        else
            jq_name = 'rotational transform'
            mode_name = 'n'
        end if
        
        ! set up plus minus one, according to the sign of the safety factor
        if (min_jq.lt.0 .and. max_jq.lt.0) then
            pmone = -1
        else if (min_jq.gt.0 .and. max_jq.gt.0) then
            pmone = 1
        else
            err_msg = trim(jq_name)//' cannot change sign'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! calculate upper and lower limits
        lim_lo = max(min_jq-tol_norm_r,min_jq/(1+pmone*tol_norm_r))
        lim_hi = min(max_jq+tol_norm_r,max_jq/(1-pmone*tol_norm_r))
        
        ! initialize min_sec_X and max_sec_X
        min_sec_X = huge(1._dp)
        max_sec_X = -huge(1._dp)
        
        ! for every mode (n,m) check whether m/n is inside the range of q values
        ! or n/m inside the range of iota values
        do id = 1, n_mod
            ! check if limits are met
            if (use_pol_flux_F) then
                if ((min_m_X+id-1._dp)/min_n_X.lt.lim_lo .or. &
                    &(min_m_X+id-1._dp)/min_n_X.gt.lim_hi) then
                    call writo('for (n,m) = ('//trim(i2str(min_n_X))//&
                        &','//trim(i2str(min_m_X+id-1))//'), there is no range &
                        &in the plasma where the ratio |n q - m| << |n|,|m| &
                        &is met')
                    ierr = 1
                    err_msg = 'Choose n and m so that |n q - m| << |n|,|m|'
                    CHCKERR(err_msg)
                end if
            else
                if ((min_n_X+id-1._dp)/min_m_X.lt.lim_lo .or. &
                    &(min_n_X+id-1._dp)/min_m_X.gt.lim_hi) then
                    call writo('for (n,m) = ('//trim(i2str(min_n_X+id-1))//&
                        &','//trim(i2str(min_m_X))//'), there is no range &
                        &in the plasma where the ratio |n - iota m| << |m|,|n| &
                        &is met')
                    ierr = 1
                    err_msg = 'Choose n and m so that |n - iota m| << |n|,|m|'
                    CHCKERR(err_msg)
                end if
            end if
        end do
        
        ! set minimum and max secondary mode
        if (use_pol_flux_F) then
            min_sec_X = min_n_X*lim_lo
            max_sec_X = min_n_X*lim_hi
        else
            min_sec_X = min_m_X*lim_lo
            max_sec_X = min_m_X*lim_hi
        end if
        
        ! output message
        call writo('The modes are all within the allowed range of '//&
            &trim(i2str(ceiling(min_sec_X)))//' < '//mode_name//' < '//&
            &trim(i2str(floor(max_sec_X)))//'...')
        
        call lvl_ud(-1)
        call writo('Mode numbers checked')
    end function check_X_modes
    
    ! calculate U_m^0, U_m^1 or U_n^0, U_n^1  at eq loc_n_r values of the normal
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
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        
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
        allocate(djq(grid%loc_n_r),mn(X%n_mod))
        allocate(fac_n(grid%loc_n_r),fac_m(grid%loc_n_r))
        if (use_pol_flux_F) then
            djq = eq%q_saf_FD(:,1)
            fac_n = eq%q_saf_FD(:,0)
            fac_m = 1.0_dp
            mn = X%n(:)
        else
            djq = -eq%rot_t_FD(:,1)
            fac_n = 1.0_dp
            fac_m = eq%rot_t_FD(:,0)
            mn = X%m(:)
        end if
        
        ! allocate helper variables
        allocate(U_corr(grid%n(1),grid%n(2),grid%loc_n_r,X%n_mod))
        allocate(D3U_corr(grid%n(1),grid%n(2),grid%loc_n_r,X%n_mod))
        allocate(g_frac(grid%n(1),grid%n(2),grid%loc_n_r))
        allocate(Theta_3(grid%n(1),grid%n(2),grid%loc_n_r))
        allocate(D1Theta_3(grid%n(1),grid%n(2),grid%loc_n_r))
        allocate(D3Theta_3(grid%n(1),grid%n(2),grid%loc_n_r))
        
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
            allocate(D3g_frac(grid%n(1),grid%n(2),grid%loc_n_r))
            allocate(D13Theta_3(grid%n(1),grid%n(2),grid%loc_n_r))
            allocate(D33Theta_3(grid%n(1),grid%n(2),grid%loc_n_r))
            
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
                do kd = 1,grid%loc_n_r
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
            use num_vars, only: norm_disc_prec_X
            
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
                do kd = 1,grid%loc_n_r
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(Theta_3(:,id,kd),D3Theta_3(:,id,kd),&
                            &grid%theta_F(:,id,kd),1,norm_disc_prec_X+1)        ! higher precision because other derivative will be taken later
                        CHCKERR('')
                    end do
                    if (.not.use_pol_flux_F) D3Theta_3(:,:,kd) = &
                        &D3Theta_3(:,:,kd)*eq%rot_t_FD(kd,0)                    ! parallel deriv. equal to iota * poloidal deriv.
                end do
                
                ! loop over all normal points
                do kd = 1,grid%loc_n_r
                    ! set up correction to U
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*&
                        &(g_frac(:,:,kd)*(D3Theta_3(:,:,kd)+iu*n_frac*&
                        &Theta_3(:,:,kd)))
                    ! calculate X%U_0 and X%DU_0
                    X%U_0(:,:,kd,jd) = &
                        &-(h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd))&
                        &+ iu/(mn(jd)*g33(:,:,kd)) * (g13(:,:,kd)*djq(kd) + &
                        &J(:,:,kd)**2*eq%pres_FD(kd,1)*vac_perm + iu*n_frac * &
                        &( g13(:,:,kd)*djq(kd)*ang_par_F(:,:,kd) - &
                        &g23(:,:,kd) )) + U_corr(:,:,kd,jd)
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_0(:,id,kd,jd),D3var(:,id),&
                            &grid%theta_F(:,id,kd),1,norm_disc_prec_X+1)        ! higher precision because previous derivative
                        CHCKERR('')
                    end do
                    if (.not.use_pol_flux_F) D3var = D3var*eq%rot_t_FD(kd,0)    ! parallel deriv. equal to iota * poloidal deriv.
                    X%DU_0(:,:,kd,jd) = D3var + iu*n_frac*X%U_0(:,:,kd,jd)
                    ! calculate X%U_1 and X%DU_1
                    X%U_1(:,:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,:,kd)/g33(:,:,kd))
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_1(:,id,kd,jd),D3var(:,id),&
                            &grid%theta_F(:,id,kd),1,norm_disc_prec_X+1)
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
            use num_vars, only: norm_disc_prec_X
            
            character(*), parameter :: rout_name = 'test_DU'
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            complex(dp), allocatable :: DU_0(:,:,:)                             ! alternative calculation for DU_0
            complex(dp), allocatable :: DU_1(:,:,:)                             ! alternative calculation for DU_1
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
            allocate(DU_0(grid%n(1),grid%n(2),grid%n(3)))
            allocate(DU_1(grid%n(1),grid%n(2),grid%n(3)))
            
            ! loop over all modes
            do jd = 1,X%n_mod
                ! loop over all normal points
                do kd = 1,grid%n(3)
                    ! derive numerically
                    do id = 1,grid%n(2)
                        ierr = calc_deriv(X%U_0(:,id,kd,jd),DU_0(:,id,kd),&
                            &ang_par_F(:,id,kd),1,norm_disc_prec_X+1)
                        CHCKERR('')
                        ierr = calc_deriv(X%U_1(:,id,kd,jd),DU_1(:,id,kd),&
                            &ang_par_F(:,id,kd),1,norm_disc_prec_X+1)
                        CHCKERR('')
                    end do
                    
                    ! add the second part
                    n_frac = X%n(jd)*fac_n(kd)-X%m(jd)*fac_m(kd)
                    DU_0(:,:,kd) = DU_0(:,:,kd) + iu*n_frac*X%U_0(:,:,kd,jd)
                    DU_1(:,:,kd) = DU_1(:,:,kd) + iu*n_frac*X%U_1(:,:,kd,jd)
                end do
                
                ! set some variables
                file_name = 'TEST_RE_DU_0_'//trim(i2str(X%n(jd)))//&
                        &'_'//trim(i2str(X%m(jd)))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(jd)//', real part')
                
                ! plot difference for RE DU_0
                call plot_diff_HDF5(realpart(DU_0),realpart(X%DU_0(:,:,:,jd)),&
                    &file_name,description=description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_0_'//trim(i2str(X%n(jd)))//&
                        &'_'//trim(i2str(X%m(jd)))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(jd)//', imaginary part')
                
                ! plot difference for IM DU_0
                call plot_diff_HDF5(imagpart(DU_0),imagpart(X%DU_0(:,:,:,jd)),&
                    &file_name,description=description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_RE_DU_1_'//trim(i2str(X%n(jd)))//&
                        &'_'//trim(i2str(X%m(jd)))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(jd)//', real part')
                
                ! plot difference for RE DU_1
                call plot_diff_HDF5(realpart(DU_1),realpart(X%DU_1(:,:,:,jd)),&
                    &file_name,description=description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_1_'//trim(i2str(X%n(jd)))//&
                        &'_'//trim(i2str(X%m(jd)))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(jd)//', imaginary part')
                
                ! plot difference for IM DU_1
                call plot_diff_HDF5(imagpart(DU_1),imagpart(X%DU_1(:,:,:,jd)),&
                    &file_name,description=description,output_message=.true.)
            end do
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end function test_DU
#endif
    end function calc_U
    
    ! calculate  ~PV_(k,m)^i  (pol.  flux)  or ~PV_(l,n)^i  (tor.  flux) at  all
    ! eq loc_n_r values
    ! (see [ADD REF] for details)
    integer function calc_PV(eq,grid,met,X_a,X_b,X,lim_sec_X) result(ierr)
        use num_vars, only: use_pol_flux_F
        use eq_vars, only: vac_perm
        use utilities, only: c
        use X_vars, only: is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        
        character(*), parameter :: rout_name = 'calc_PV'
        
        ! use input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(in) :: X_a, X_b                                  ! vectorial perturbation variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:,:)                                 ! common factor |nabla psi|^2/(J^2*B^2)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:) => null()                                 ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:) => null()                               ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:) => null()                               ! h^psi,psi
        
        ! initialize ierr
        ierr = 0
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        ! lower metric factors
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        ! upper metric factors
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up common factor for PV_i
        allocate(com_fac(grid%n(1),grid%n(2),grid%loc_n_r))
        com_fac = h22/(g33*vac_perm)
        
        ! set up fac_n and fac_m
        allocate(fac_n(grid%loc_n_r),fac_m(grid%loc_n_r))
        if (use_pol_flux_F) then
            fac_n = eq%q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = eq%rot_t_FD(:,0)
        end if
        
        ! set nr. of modes
        n_mod_tot = (max_m_X-min_m_X+1)*(max_n_X-min_n_X+1)
        
        ! loop over all modes
        do m = 1,X_b%n_mod
            do k = 1,X_a%n_mod
                ! check whether modes are correct
                if (X%n_1(k).ne.X_a%n(k) .or. X%n_2(m).ne.X_b%n(m) .or. &
                    &X%m_1(k).ne.X_a%m(k) .or. X%m_2(m).ne.X_b%m(m)) then
                    ierr = 1
                    CHCKERR('Modes do not match')
                end if
                
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(X,.true.,[k,m])
                calc_this(2) = is_necessary_X(X,.false.,[k,m])
                
                ! set up c_loc
                c_loc(1) = c([k,m],.true.,n_mod_tot,lim_sec_X)
                c_loc(2) = c([k,m],.false.,n_mod_tot,lim_sec_X)
                
                ! calculate PV_0
                if (calc_this(1)) then
                    X%PV_0(:,:,:,c_loc(1)) = &
                        &com_fac*(X_b%DU_0(:,:,:,m) - eq%S*J - &
                        &eq%sigma/(com_fac*J)) * (conjg(&
                        &X_a%DU_0(:,:,:,k)) - eq%S*J - eq%sigma/(com_fac*J)) - &
                        &eq%sigma/J * (eq%S*J + eq%sigma/(com_fac*J))
                    
                    ! add (nq-k)*(nq-m)/(mu_0J^2 |nabla  psi|^2) - 2p'kappa_n to
                    ! PV_0
                    do kd = 1,grid%loc_n_r
                        X%PV_0(:,:,kd,c_loc(1)) = X%PV_0(:,:,kd,c_loc(1)) + &
                            &(X_b%n(m)*fac_n(kd)-X_b%m(m)*fac_m(kd))*&
                            &(X_a%n(k)*fac_n(kd)-X_a%m(k)*fac_m(kd)) / &
                            &( vac_perm*J(:,:,kd)**2*h22(:,:,kd) ) - &
                            &2*eq%pres_FD(kd,1)*eq%kappa_n(:,:,kd)
                    end do
                end if
                
                ! calculate PV_1
                if (calc_this(2)) then
                    X%PV_1(:,:,:,c_loc(2)) = com_fac * X_b%DU_1(:,:,:,m) * &
                        &(conjg(X_a%DU_0(:,:,:,k)) - eq%S*J - &
                        &eq%sigma/(com_fac*J))
                end if
                
                ! calculate PV_2
                if (calc_this(1)) then
                    X%PV_2(:,:,:,c_loc(1)) = &
                        &com_fac*X_b%DU_1(:,:,:,m)*conjg(X_a%DU_1(:,:,:,k))
                end if
            end do
        end do
        
        ! deallocate variables
        nullify(J,g33,h22)
    end function calc_PV
    
    ! calculate  ~KV_(k,m)^i  (pol.  flux)  or ~KV_(l,n)^i  (tor.  flux) at  all
    ! eq loc_n_r values
    ! (see [ADD REF] for details)
    integer function calc_KV(eq,grid,met,X_a,X_b,X,lim_sec_X) result(ierr)
        use utilities, only: c
        use X_vars, only: is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        
        character(*), parameter :: rout_name = 'calc_KV'
        
        ! use input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid                                     ! grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(in) :: X_a, X_b                                  ! vectorial perturbation variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:,:)                                 ! common factor |nabla psi|^2/(J^2*B^2)
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        
        ! submatrices
        ! jacobian
        real(dp), pointer :: J(:,:,:) => null()                                 ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:) => null()                               ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:) => null()                               ! h^psi,psi
        
        ! initialize ierr
        ierr = 0
        
        ! set up submatrices
        ! jacobian
        J => met%jac_FD(:,:,:,0,0,0)
        ! lower metric factors
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        ! upper metric factors
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up common factor
        allocate(com_fac(grid%n(1),grid%n(2),grid%loc_n_r))
        com_fac = J**2*h22/g33
        
        ! set nr. of modes
        n_mod_tot = (max_m_X-min_m_X+1)*(max_n_X-min_n_X+1)
        
        ! loop over the modes
        do m = 1,X_b%n_mod
            do k = 1,X_a%n_mod
                ! check whether modes are correct
                if (X%n_1(k).ne.X_a%n(k) .or. X%n_2(m).ne.X_b%n(m) .or. &
                    &X%m_1(k).ne.X_a%m(k) .or. X%m_2(m).ne.X_b%m(m)) then
                    ierr = 1
                    CHCKERR('Modes do not match')
                end if
                
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(X,.true.,[k,m])
                calc_this(2) = is_necessary_X(X,.false.,[k,m])
                
                ! set up c_loc
                c_loc(1) = c([k,m],.true.,n_mod_tot,lim_sec_X)
                c_loc(2) = c([k,m],.false.,n_mod_tot,lim_sec_X)
                
                ! calculate KV_0
                if (calc_this(1)) then
                    X%KV_0(:,:,:,c_loc(1)) = com_fac * &
                        &X_b%U_0(:,:,:,m) * conjg(X_a%U_0(:,:,:,k)) + 1._dp/h22
                end if
                
                ! calculate KV_1
                if (calc_this(2)) then
                    X%KV_1(:,:,:,c_loc(2)) = com_fac * &
                        &X_b%U_1(:,:,:,m) * conjg(X_a%U_0(:,:,:,k))
                end if
                
                ! calculate KV_2
                if (calc_this(1)) then
                    X%KV_2(:,:,:,c_loc(1)) = com_fac * &
                        &X_b%U_1(:,:,:,m) * conjg(X_a%U_1(:,:,:,k))
                end if
            end do
        end do
        
        ! multiply by rho
        do kd = 1,grid%loc_n_r
            X%KV_0(:,:,kd,:) = X%KV_0(:,:,kd,:)*eq%rho(kd)
            X%KV_1(:,:,kd,:) = X%KV_1(:,:,kd,:)*eq%rho(kd)
            X%KV_2(:,:,kd,:) = X%KV_2(:,:,kd,:)*eq%rho(kd)
        end do
        
        ! deallocate variables
        nullify(J,g33,h22)
    end function calc_KV
    
    ! Calculate the  magnetic integrals  from PV_i and  KV_i. All  the variables
    ! should thus be field-line oriented.
    subroutine calc_magn_ints(grid_eq,met,X,lim_sec_X)
        use num_vars, only: use_pol_flux_F
        use X_vars, only: is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        use utilities, only: c
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        integer :: k, m                                                         ! counters
        integer :: id                                                           ! counter
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        complex(dp), allocatable :: J_exp_ang(:,:,:)                            ! J * exponential of Flux parallel angle
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle in flux coordinates
        
        ! user output
        call writo('Calculating field-line averages')
        call lvl_ud(1)
        
        ! set nr. of modes
        n_mod_tot = (max_m_X-min_m_X+1)*(max_n_X-min_n_X+1)
        
        ! allocate J_exp_ang
        allocate(J_exp_ang(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_F) then
            ang_par_F => grid_eq%theta_F
        else
            ang_par_F => grid_eq%zeta_F
        end if
        
        ! initialize integrated quantities
        X%PV_int_0 = 0
        X%PV_int_1 = 0
        X%PV_int_2 = 0
        X%KV_int_0 = 0
        X%KV_int_1 = 0
        X%KV_int_2 = 0
        
        ! loop over all modes
        do m = 1,X%n_mod(2)
            do k = 1,X%n_mod(1)
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(X,.true.,[k,m])
                calc_this(2) = is_necessary_X(X,.false.,[k,m])
                
                ! set up c_loc
                c_loc(1) = c([k,m],.true.,n_mod_tot,lim_sec_X)
                c_loc(2) = c([k,m],.false.,n_mod_tot,lim_sec_X)
                
                ! calculate J_exp_ang
                if (use_pol_flux_F) then
                    J_exp_ang(:,:,:) = met%jac_FD(:,:,:,0,0,0)*&
                        &exp(iu*(X%m_1(k)-X%m_2(m))*ang_par_F)
                else
                    J_exp_ang(:,:,:) = met%jac_FD(:,:,:,0,0,0)*&
                        &exp(-iu*(X%n_1(k)-X%n_2(m))*ang_par_F)
                end if
                
                ! parallel integration loop
                do id = 2,grid_eq%n(1)
                    if (calc_this(1)) then
                        X%PV_int_0(c_loc(1),:,:) = X%PV_int_0(c_loc(1),:,:) + &
                            &(J_exp_ang(id,:,:)*X%PV_0(id,:,:,c_loc(1))&
                            &+J_exp_ang(id-1,:,:)*X%PV_0(id-1,:,:,c_loc(1)))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                        X%PV_int_2(c_loc(1),:,:) = X%PV_int_2(c_loc(1),:,:) + &
                            &(J_exp_ang(id,:,:)*X%PV_2(id,:,:,c_loc(1))&
                            &+J_exp_ang(id-1,:,:)*X%PV_2(id-1,:,:,c_loc(1)))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                        X%KV_int_0(c_loc(1),:,:) = X%KV_int_0(c_loc(1),:,:) + &
                            &(J_exp_ang(id,:,:)*X%KV_0(id,:,:,c_loc(1))&
                            &+J_exp_ang(id-1,:,:)*X%KV_0(id-1,:,:,c_loc(1)))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                        X%KV_int_2(c_loc(1),:,:) = X%KV_int_2(c_loc(1),:,:) + &
                            &(J_exp_ang(id,:,:)*X%KV_2(id,:,:,c_loc(1))&
                            &+J_exp_ang(id-1,:,:)*X%KV_2(id-1,:,:,c_loc(1)))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                    end if
                    if (calc_this(2)) then
                        X%PV_int_1(c_loc(2),:,:) = X%PV_int_1(c_loc(2),:,:) + &
                            &(J_exp_ang(id,:,:)*X%PV_1(id,:,:,c_loc(2))&
                            &+J_exp_ang(id-1,:,:)*X%PV_1(id-1,:,:,c_loc(2)))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                        X%KV_int_1(c_loc(2),:,:) = X%KV_int_1(c_loc(2),:,:) + &
                            &(J_exp_ang(id,:,:)*X%KV_1(id,:,:,c_loc(2))&
                            &+J_exp_ang(id-1,:,:)*X%KV_1(id-1,:,:,c_loc(2)))/2 &
                            &*(ang_par_F(id,:,:)-ang_par_F(id-1,:,:))
                    end if
                end do
            end do
        end do
        
        ! clean up
        nullify(ang_par_F)
        
        ! user output
        call lvl_ud(-1)
        call writo('Field-line averages calculated')
    end subroutine calc_magn_ints
    
    ! set-up  and solve  the  EV system  by discretizing  the  equations in  the
    ! perturbation  grid,  making  use  of   PV  and  KV,  interpolated  in  the
    ! equilibrium grid.
    integer function solve_EV_system(grid_eq,grid_X,X,sol,use_guess,&
        &n_sol_found) result(ierr)
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use SLEPC_ops, only: solve_EV_system_SLEPC
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables
        type(sol_type), intent(inout) :: sol                                    ! tensorial perturbation variables
        logical, intent(in) :: use_guess                                        ! whether to use a guess or not
        integer, intent(inout) :: n_sol_found                                   ! how many solutions saved
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(grid_eq,grid_X,X,sol,use_guess,&
                    &n_sol_found)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
    
    ! Print either  vectorial or tensorial perturbation quantities  of a certain
    ! order to an output file:
    !   - vectorial:    U, DU
    !   - tensorial:    PV_int, KV_int
    !     (the non-integrated variables are heavy and not requested)
    ! Note: Flux coordinates used as normal coordinates
    integer function print_output_X_1(grid_eq,X) result(ierr)                   ! vectorial version
        use num_vars, only: PB3D_name, max_it_r, rich_lvl_nr
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        use X_vars, only: get_suffix
        
        character(*), parameter :: rout_name = 'print_output_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(X_1_type), intent(in) :: X                                         ! vectorial perturbation variables 
        
        ! local variables
        type(var_1D_type), allocatable, target :: X_1D(:)                       ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: X_1D_loc => null()                        ! local element in X_1D
        integer :: id, jd                                                       ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Writing vectorial perturbation variables to output file')
        call lvl_ud(1)
        
        ! user output
        call writo('Preparing variables for writing')
        call lvl_ud(1)
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(X_1D(8*X%n_mod))
        
        ! set up variables X_1D
        id = 1
        
        ! RE_U_0
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_U_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_0(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%U_0(:,:,:,jd)),&
                &[size(X%U_0(:,:,:,jd))])
        end do
        
        ! IM_U_0
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_U_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_0(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%U_0(:,:,:,jd)),&
                &[size(X%U_0(:,:,:,jd))])
        end do
        
        ! RE_U_1
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_U_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_1(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%U_1(:,:,:,jd)),&
                &[size(X%U_1(:,:,:,jd))])
        end do
        
        ! IM_U_1
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_U_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_1(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%U_1(:,:,:,jd)),&
                &[size(X%U_1(:,:,:,jd))])
        end do
        
        ! RE_DU_0
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_DU_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_0(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%DU_0(:,:,:,jd)),&
                &[size(X%DU_0(:,:,:,jd))])
        end do
        
        ! IM_DU_0
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_DU_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_0(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%DU_0(:,:,:,jd)),&
                &[size(X%DU_0(:,:,:,jd))])
        end do
        
        ! RE_DU_1
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_DU_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_1(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%DU_1(:,:,:,jd)),&
                &[size(X%DU_1(:,:,:,jd))])
        end do
        
        ! IM_DU_1
        do jd = 1,X%n_mod
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_DU_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid_eq%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_1(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%DU_1(:,:,:,jd)),&
                &[size(X%DU_1(:,:,:,jd))])
        end do
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Writing using HDF5')
        call lvl_ud(1)
        
        ! write
        if (max_it_r.gt.1) then
            ierr = print_HDF5_arrs(X_1D,PB3D_name,&
                &'X_1_R_'//trim(i2str(rich_lvl_nr)))
        else
            ierr = print_HDF5_arrs(X_1D,PB3D_name,'X_1')
        end if
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        
        ! clean up
        nullify(X_1D_loc)
        
        ! user output
        call lvl_ud(-1)
        call writo('Vectorial perturbation variables written to output')
    end function print_output_X_1
    integer function print_output_X_2(grid_eq,X) result(ierr)                   ! tensorial version
        use num_vars, only: rich_lvl_nr, max_it_r, PB3D_name, use_pol_flux_F
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        use X_vars, only: get_suffix, is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'print_output_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables 
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        type(var_1D_type), allocatable, target :: X_1D(:)                       ! 1D equivalent of eq. variables
        type(var_1D_type), pointer :: X_1D_loc => null()                        ! local element in X_1D
        integer :: id, jd, kd                                                   ! counters
        integer :: lim_submat(2,2)                                              ! limits of submatrix
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: print_this(2)                                                ! whether symmetric and asymmetric variables need to be printed
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Writing tensorial perturbation variables to output file')
        call lvl_ud(1)
        
        ! user output
        call writo('Preparing variables for writing')
        call lvl_ud(1)
        
        ! set nr. of modes
        n_mod_tot = (max_m_X-min_m_X+1)*(max_n_X-min_n_X+1)
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(X_1D(&
            &4*(size(X%PV_int_0,1)+size(X%PV_int_1,1)+size(X%PV_int_2,1))))
        
        ! set limits of submatrix in whole matrix
        if (use_pol_flux_F) then
            lim_submat(:,1) = [minval(X%m_1),maxval(X%m_1)]-min_m_X+1
            lim_submat(:,2) = [minval(X%m_2),maxval(X%m_2)]-min_m_X+1
        else
            lim_submat(:,1) = [minval(X%n_1),maxval(X%n_1)]-min_n_X+1
            lim_submat(:,2) = [minval(X%n_2),maxval(X%n_2)]-min_n_X+1
        end if
        id = 1
        
        do kd = 1,X%n_mod(2)
            do jd = 1,X%n_mod(1)
                ! set local c's
                c_loc(1) = c([jd,kd],.true.,n_mod_tot,lim_submat)
                c_loc(2) = c([jd,kd],.false.,n_mod_tot,lim_submat)
                
                ! determine whether variables need to be printed
                print_this(1) = is_necessary_X(X,.true.,[jd,kd])
                print_this(2) = is_necessary_X(X,.false.,[jd,kd])
                
                ! RE_PV_int_0
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'RE_PV_int_0_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%PV_int_0(1,:,:))))
                    X_1D_loc%p = reshape(realpart(&
                        &X%PV_int_0(c_loc(1),:,:)),[size(X%PV_int_0(1,:,:))])
                end if
                
                ! IM_PV_int_0
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'IM_PV_int_0_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%PV_int_0(1,:,:))))
                    X_1D_loc%p = reshape(imagpart(&
                        &X%PV_int_0(c_loc(1),:,:)),[size(X%PV_int_0(1,:,:))])
                end if
                
                ! RE_PV_int_1
                if (print_this(2)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'RE_PV_int_1_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%PV_int_1(1,:,:))))
                    X_1D_loc%p = reshape(realpart(&
                        &X%PV_int_1(c_loc(2),:,:)),[size(X%PV_int_1(1,:,:))])
                end if
                
                ! IM_PV_int_1
                if (print_this(2)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'IM_PV_int_1_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%PV_int_1(1,:,:))))
                    X_1D_loc%p = reshape(imagpart(&
                        &X%PV_int_1(c_loc(2),:,:)),[size(X%PV_int_1(1,:,:))])
                end if
                
                ! RE_PV_int_2
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'RE_PV_int_2_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%PV_int_2(1,:,:))))
                    X_1D_loc%p = reshape(realpart(&
                        &X%PV_int_2(c_loc(1),:,:)),[size(X%PV_int_2(1,:,:))])
                end if
                
                ! IM_PV_int_2
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'IM_PV_int_2_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%PV_int_2(1,:,:))))
                    X_1D_loc%p = reshape(imagpart(&
                        &X%PV_int_2(c_loc(1),:,:)),[size(X%PV_int_2(1,:,:))])
                end if
                
                ! RE_KV_int_0
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'RE_KV_int_0_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%KV_int_0(1,:,:))))
                    X_1D_loc%p = reshape(realpart(&
                        &X%KV_int_0(c_loc(1),:,:)),[size(X%KV_int_0(1,:,:))])
                end if
                
                ! IM_KV_int_0
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'IM_KV_int_0_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%KV_int_0(1,:,:))))
                    X_1D_loc%p = reshape(imagpart(&
                        &X%KV_int_0(c_loc(1),:,:)),[size(X%KV_int_0(1,:,:))])
                end if
                
                ! RE_KV_int_1
                if (print_this(2)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'RE_KV_int_1_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%KV_int_1(1,:,:))))
                    X_1D_loc%p = reshape(realpart(&
                        &X%KV_int_1(c_loc(2),:,:)),[size(X%KV_int_1(1,:,:))])
                end if
                
                ! IM_KV_int_1
                if (print_this(2)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'IM_KV_int_1_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%KV_int_1(1,:,:))))
                    X_1D_loc%p = reshape(imagpart(&
                        &X%KV_int_1(c_loc(2),:,:)),[size(X%KV_int_1(1,:,:))])
                end if
                
                ! RE_KV_int_2
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'RE_KV_int_2_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%KV_int_2(1,:,:))))
                    X_1D_loc%p = reshape(realpart(&
                        &X%KV_int_2(c_loc(1),:,:)),[size(X%KV_int_2(1,:,:))])
                end if
                
                ! IM_KV_int_2
                if (print_this(1)) then
                    X_1D_loc => X_1D(id); id = id+1
                    X_1D_loc%var_name = 'IM_KV_int_2_'//&
                        &trim(get_suffix(X,jd,kd))
                    allocate(X_1D_loc%tot_i_min(2),X_1D_loc%tot_i_max(2))
                    allocate(X_1D_loc%loc_i_min(2),X_1D_loc%loc_i_max(2))
                    X_1D_loc%tot_i_min = [1,1]
                    X_1D_loc%tot_i_max = grid_eq%n(2:3)
                    X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
                    X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
                    allocate(X_1D_loc%p(size(X%KV_int_2(1,:,:))))
                    X_1D_loc%p = reshape(imagpart(&
                        &X%KV_int_2(c_loc(1),:,:)),[size(X%KV_int_2(1,:,:))])
                end if
            end do
        end do
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Writing using HDF5')
        call lvl_ud(1)
        
        ! write
        if (max_it_r.gt.1) then
            ierr = print_HDF5_arrs(X_1D,PB3D_name,&
                &'X_2_R_'//trim(i2str(rich_lvl_nr)))
        else
            ierr = print_HDF5_arrs(X_1D,PB3D_name,'X_2')
        end if
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        
        ! clean up
        nullify(X_1D_loc)
        
        ! user output
        call lvl_ud(-1)
        call writo('Tensorial perturbation variables written to output')
    end function print_output_X_2
end module
