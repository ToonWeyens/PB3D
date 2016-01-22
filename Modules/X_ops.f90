!------------------------------------------------------------------------------!
!   Operations considering perturbation quantities                             !
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, max_name_ln, pi
    use grid_vars, onlY: grid_type, dealloc_grid
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type, X_2_type

    implicit none
    private
    public calc_X, calc_magn_ints, print_output_X, resonance_plot, &
        &calc_res_surf, check_X_modes
    
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
    integer function calc_X_1(grid_eq,grid_X,eq,met,X,lim_sec_X) result(ierr)   ! vectorial version
        use X_vars, only: create_X
        
        character(*), parameter :: rout_name = 'calc_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid variables
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
        call create_X(grid_X,X,lim_sec_X)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        ierr = calc_U(grid_eq,grid_X,eq,met,X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
        call writo('Vectorial perturbation variables calculated')
    end function calc_X_1
    integer function calc_X_2(grid_eq,grid_X,eq,met,X_a,X_b,X,lim_sec_X) &
        &result(ierr)                                                           ! tensorial version
        use X_vars, only: create_X
        
        character(*), parameter :: rout_name = 'calc_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid variables
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
        call create_X(grid_X,X,lim_sec_X)
        
        ! Calculate  PV_i  for  all (k,m)  pairs  and n_r  (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating PV...')
        call lvl_ud(1)
        ierr = calc_PV(grid_eq,grid_X,eq,met,X_a,X_b,X,lim_sec_X)
        CHCKERR('')
        call lvl_ud(-1)
        
        !  Calculate KV_i  for  all (k,m)  pairs  and n_r  (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV...')
        call lvl_ud(1)
        ierr = calc_KV(grid_eq,grid_X,eq,met,X_a,X_b,X,lim_sec_X)
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
        
        ! only master and only if resonant surfaces
        if (rank.eq.0 .and. size(res_surf,1).gt.0) then
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
                Y_plot(1,kd,1,1) = res_surf(kd,1)-1
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
        use num_vars, only: use_pol_flux_F, tol_norm
        
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
        call writo('The tolerance used is '//trim(r2strt(tol_norm))//&
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
        lim_lo = max(min_jq-tol_norm,min_jq/(1+pmone*tol_norm))
        lim_hi = min(max_jq+tol_norm,max_jq/(1-pmone*tol_norm))
        
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
    ! Use is made  of variables Ti that  are set up in the  equilibrium grid and
    ! are afterwards converted to Ti_X in  the perturbation grid, which needs to
    ! have the same angular coordinates as the equilibrium grid:
    !   - T1 = B_alpha/B_theta
    !   - T2 = Theta^alpha + q' theta
    !   - T3 = B_alpha/B_theta q' + J^2/B_theta mu_0 p'
    !   - T4 = B_alpha/B_theta q' theta - B_psi/B_theta
    !   - T5 = B_alpha/B_theta D3Theta^theta - D1Theta^theta
    !   - T6 = B_alpha/B_theta Theta^theta
    ! which  is valid  for poloidal  Flux  coordinates and  where q'  has to  be
    ! replaced by -iota' and theta by zeta for toroidal Flux coordinates.
    ! The interpolated Ti_X are then used to calculate U:
    !   - U_0 = -T2                                         (order 1)
    !       + (i/n) (T3 + i(nq-m) T4)                       (order 2)
    !       + (i/n)^2 i(nq-m)(-T5 - (nq-m) T6)              (order 3)
    !   - U_1 = i/n                                         (order 1)
    !       + (i/n)^2 i(nq-m) (-T1)                         (order 2)
    ! which is valid for poloidal Flux coordinates and where n is to be replaced
    ! by m and (nq-m) by (n-iota m) for toroidal Flux coordinates.
    ! For VMEC, these factors are also derived in the parallel coordinate.
    integer function calc_U(grid_eq,grid_X,eq,met,X) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, U_style, norm_disc_prec_X
        use utilities, only: c, calc_deriv
        use input_ops, only: get_log, pause_prog
        use eq_vars, only: vac_perm
        use grid_ops, only: get_norm_interp_data
#if ldebug
        use num_vars, only: ltest
#endif
        
        character(*), parameter :: rout_name = 'calc_U'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(inout) :: X                                      ! vectorial perturbation variables
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: loc_r_eq(:)                                    ! unrounded index of eq variables in X grid
        integer :: i_lo, i_hi                                                   ! upper and lower index
        integer :: T_size                                                       ! 2 for VMEC and 1 for HELENA
        
        ! Jacobian
        real(dp), pointer :: J(:,:,:)                                           ! Jacobian
        real(dp), pointer :: D3J(:,:,:)                                         ! D_theta Jacobian
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
        real(dp), pointer :: D13h22(:,:,:)                                      ! D_alpha,theta h^psi,psi
        real(dp), pointer :: D33h22(:,:,:)                                      ! D_theta,theta h^psi,psi
        real(dp), pointer :: h23(:,:,:)                                         ! h^psi,theta
        real(dp), pointer :: D1h23(:,:,:)                                       ! D_alpha h^psi,theta
        real(dp), pointer :: D3h23(:,:,:)                                       ! D_theta h^psi,theta
        real(dp), pointer :: D13h23(:,:,:)                                      ! D_alpha,theta h^psi,theta
        real(dp), pointer :: D33h23(:,:,:)                                      ! D_theta,theta h^psi,theta
        ! helper variables
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! parallel angle in flux coordinates
        real(dp), pointer :: ang_geo_F(:,:,:)                                   ! geodesical angle in flux coordinates
        real(dp), allocatable :: djq(:)                                         ! either q' (pol. flux) or -iota' (tor. flux)
        real(dp), allocatable :: Theta_3(:,:,:), D1Theta_3(:,:,:), &
            &D3Theta_3(:,:,:)                                                   ! Theta^theta and derivatives
        real(dp), allocatable :: D13Theta_3(:,:,:), D33Theta_3(:,:,:)           ! Theta^theta and derivatives
        real(dp), allocatable :: n_frac(:,:)                                    ! nq-m (pol. flux) or n-iotam (tor. flux) for all modes
        real(dp), allocatable :: mn(:)                                          ! either n*A_0 (pol. flux) or m (tor.flux)
        complex(dp), allocatable :: U_corr(:,:,:,:)                             ! correction to U_0 and U_1 for a certain (n,m)
        complex(dp), allocatable :: D3U_corr(:,:,:,:)                           ! D_theta U_corr
        complex(dp), allocatable :: D3var(:)                                    ! derivative of variable
        ! U factors
        real(dp), allocatable :: T1(:,:,:,:)                                    ! B_alpha/B_theta and par. deriv.
        real(dp), allocatable :: T2(:,:,:,:)                                    ! Theta^alpha + q' theta and par. deriv.
        real(dp), allocatable :: T3(:,:,:,:)                                    ! B_alpha/B_theta q' + J^2/B_theta mu_0 p' and par. deriv.
        real(dp), allocatable :: T4(:,:,:,:)                                    ! B_alpha/B_theta q' theta - B_psi/B_theta and par. deriv.
        real(dp), allocatable :: T5(:,:,:,:)                                    ! B_alpha/B_theta D3Theta^theta - D1Theta^theta and par. deriv.
        real(dp), allocatable :: T6(:,:,:,:)                                    ! B_alpha/B_theta Theta^theta and par. deriv.
        real(dp), allocatable :: T1_X(:,:,:,:)                                  ! T1 in X grid and par. deriv.
        real(dp), allocatable :: T2_X(:,:,:,:)                                  ! T2 in X grid and par. deriv.
        real(dp), allocatable :: T3_X(:,:,:,:)                                  ! T3 in X grid and par. deriv.
        real(dp), allocatable :: T4_X(:,:,:,:)                                  ! T4 in X grid and par. deriv.
        real(dp), allocatable :: T5_X(:,:,:,:)                                  ! T5 in X grid and par. deriv.
        real(dp), allocatable :: T6_X(:,:,:,:)                                  ! T6 in X grid and par. deriv.
        
        ! initialize ierr
        ierr = 0
        
        ! message
        call writo('Calculating U up to order '//trim(i2str(U_style)))
        
        ! allocate variables
        ! helper variables
        allocate(mn(X%n_mod))
        allocate(djq(grid_eq%loc_n_r))
        allocate(Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D1Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D3Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D13Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D33Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(n_frac(grid_X%loc_n_r,X%n_mod))
        allocate(U_corr(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,2))
        allocate(D3U_corr(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,2))
        ! U factors
        select case (eq_style)
            case (1)                                                            ! VMEC
                write(*,*) '!!!!!!!!! CALC_U  NOT TESTED FOR VMEC !!!!!!!!!!!!!'
                T_size = 2
            case (2)                                                            ! HELENA
                T_size = 1
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        allocate(T1(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T2(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T4(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T5(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T6(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T1_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        allocate(T2_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        allocate(T3_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        allocate(T4_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        allocate(T5_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        allocate(T6_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        
        ! set pointers
        if (use_pol_flux_F) then
            ang_par_F => grid_eq%theta_F
            ang_geo_F => grid_eq%zeta_F
        else
            ang_par_F => grid_eq%zeta_F
            ang_geo_F => grid_eq%theta_F
        end if
        J => met%jac_FD(:,:,:,0,0,0)
        D3J => met%jac_FD(:,:,:,0,0,1)
        g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,0)
        D3g13 => met%g_FD(:,:,:,c([1,3],.true.),0,0,1)
        g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,0)
        D3g23 => met%g_FD(:,:,:,c([2,3],.true.),0,0,1)
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        D3g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,1)
        h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        D3h12 => met%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        D1h22 => met%h_FD(:,:,:,c([2,2],.true.),1,0,0)
        D3h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,1)
        D13h22 => met%h_FD(:,:,:,c([2,2],.true.),1,0,1)
        D33h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,2)
        h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,0)
        D1h23 => met%h_FD(:,:,:,c([2,3],.true.),1,0,0)
        D3h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,1)
        D13h23 => met%h_FD(:,:,:,c([2,3],.true.),1,0,1)
        D33h23 => met%h_FD(:,:,:,c([2,3],.true.),0,0,2)
        
        ! set up helper variables in eq grid
        if (use_pol_flux_F) then
            mn = X%n(:)
            djq = eq%q_saf_FD(:,1)
        else
            djq = -eq%rot_t_FD(:,1)
            mn = X%m(:)
        end if
        Theta_3 = h23/h22
        select case (eq_style)
            case (1)                                                            ! VMEC
                D1Theta_3 = D1h23/h22 - h23*D1h22/(h22**2)
                D3Theta_3 = D3h23/h22 - h23*D3h22/(h22**2)
                D13Theta_3 = D13h23/h22 - (D3h23*D1h22+D1h23*D3h22)/(h22**2) &
                    &- h23*D13h22/(h22**2) + 2*h23*D3h22*D1h22/(h22**3)
                D33Theta_3 = D33h23/h22 - 2*D3h23*D3h22/(h22**2) &
                    &- h23*D33h22/(h22**2) + 2*h23*D3h22**2/(h22**3)
            case (2)                                                            ! HELENA
                if (grid_eq%n(2).gt.norm_disc_prec_X+3) then                    ! need enough terms
                    do kd = 1,grid_eq%loc_n_r
                        do id = 1,grid_eq%n(1)
                            ierr = calc_deriv(Theta_3(id,:,kd),&
                                &D1Theta_3(id,:,kd),ang_geo_F(id,:,kd),1,&
                                &norm_disc_prec_X+1)                            ! higher precision because other derivative will be taken later
                            CHCKERR('')
                        end do
                    end do
                else
                    D1Theta_3 = 0._dp
                end if
                if (grid_eq%n(1).gt.norm_disc_prec_X+3) then                    ! need enough terms
                    do kd = 1,grid_eq%loc_n_r
                        do jd = 1,grid_eq%n(2)
                            ierr = calc_deriv(Theta_3(:,jd,kd),&
                                &D3Theta_3(:,jd,kd),ang_par_F(:,jd,kd),1,&
                                &norm_disc_prec_X+1)                            ! higher precision because other derivative will be taken later
                            CHCKERR('')
                        end do
                    end do
                else
                    D3Theta_3 = 0._dp
                end if
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! set up U factors in eq grid
        T1(:,:,:,1) = g13/g33
        do kd = 1,grid_eq%loc_n_r
            T2(:,:,kd,1) = h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd)
            T3(:,:,kd,1) = T1(:,:,kd,1)*djq(kd) + &
                &J(:,:,kd)**2*vac_perm*eq%pres_FD(kd,1)/g33(:,:,kd)
            T4(:,:,kd,1) = T1(:,:,kd,1)*djq(kd)*ang_par_F(:,:,kd) - &
                &g23(:,:,kd)/g33(:,:,kd)
        end do
        T5(:,:,:,1) = T1(:,:,:,1)*D3Theta_3 - D1Theta_3
        T6(:,:,:,1) = T1(:,:,:,1)*Theta_3
        
        ! set up parallel derivatives of U in eq grid if VMEC is used
        if (eq_style.eq.1) then
            T1(:,:,:,2) = D3g13/g33 - g13*D3g33/g33**2
            do kd = 1,grid_eq%loc_n_r
                T2(:,:,kd,2) = D3h12(:,:,kd)/h22(:,:,kd) - &
                    &h12(:,:,kd)*D3h22(:,:,kd)/h22(:,:,kd)**2 + djq(kd)
                T3(:,:,kd,2) = T1(:,:,kd,2)*djq(kd) + &
                    &J(:,:,kd)*vac_perm*eq%pres_FD(kd,1)/g33(:,:,kd) * &
                    &(2*D3J(:,:,kd)-D3g33(:,:,kd)*J(:,:,kd)/g33(:,:,kd))
                T4(:,:,kd,1) = (T1(:,:,kd,2)*ang_par_F(:,:,kd)+T1(:,:,kd,1))*&
                    &djq(kd) - D3g23(:,:,kd)/g33(:,:,kd) + &
                    &g23(:,:,kd)*D3g33(:,:,kd)/g33(:,:,kd)**2
            end do
            T5(:,:,:,2) = T1(:,:,:,2)*D3Theta_3 + T1(:,:,:,1)*D33Theta_3 - &
                &D13Theta_3
            T6(:,:,:,2) = T1(:,:,:,2)*Theta_3 + T1(:,:,:,1)*D3Theta_3
        end if
        
        ! get normal interpolation factors
        ierr = get_norm_interp_data(grid_eq,grid_X,loc_r_eq)
        CHCKERR('')
        
        ! interpolate helper variables and U factors
        do kd = 1,grid_X%loc_n_r
            ! set lower and upper index
            i_lo = floor(loc_r_eq(kd))
            i_hi = ceiling(loc_r_eq(kd))
            
            if (use_pol_flux_F) then
                n_frac(kd,:) = X%n*(eq%q_saf_FD(i_lo,0) + (loc_r_eq(kd)-i_lo)*&
                    &(eq%q_saf_FD(i_hi,0)-eq%q_saf_FD(i_lo,0))) - X%m
            else
                n_frac(kd,:) = -X%m*(eq%rot_t_FD(i_lo,0) + (loc_r_eq(kd)-i_lo)*&
                    &(eq%rot_t_FD(i_hi,0)-eq%rot_t_FD(i_lo,0))) + X%n
            end if
            T1_X(:,:,kd,:) = T1(:,:,i_lo,:) + (loc_r_eq(kd)-i_lo)*&
                &(T1(:,:,i_hi,:)-T1(:,:,i_lo,:))
            T2_X(:,:,kd,:) = T2(:,:,i_lo,:) + (loc_r_eq(kd)-i_lo)*&
                &(T2(:,:,i_hi,:)-T2(:,:,i_lo,:))
            T3_X(:,:,kd,:) = T3(:,:,i_lo,:) + (loc_r_eq(kd)-i_lo)*&
                &(T3(:,:,i_hi,:)-T3(:,:,i_lo,:))
            T4_X(:,:,kd,:) = T4(:,:,i_lo,:) + (loc_r_eq(kd)-i_lo)*&
                &(T4(:,:,i_hi,:)-T4(:,:,i_lo,:))
            T5_X(:,:,kd,:) = T5(:,:,i_lo,:) + (loc_r_eq(kd)-i_lo)*&
                &(T5(:,:,i_hi,:)-T5(:,:,i_lo,:))
            T6_X(:,:,kd,:) = T6(:,:,i_lo,:) + (loc_r_eq(kd)-i_lo)*&
                &(T6(:,:,i_hi,:)-T6(:,:,i_lo,:))
        end do
        
        ! calculate U_0 and U_1 for all modes
        do ld = 1,X%n_mod
            ! calculate order 1 of U_0 and U_1
            if (U_style.ge.1) then
                X%U_0(:,:,:,ld) = -T2_X(:,:,:,1)
                X%U_1(:,:,:,ld) = iu/mn(ld)
            end if
            ! include order 2
            if (U_style.ge.2) then
                do kd = 1,grid_X%loc_n_r
                    U_corr(:,:,kd,1) = T3_X(:,:,kd,1) + &
                        &iu*n_frac(kd,ld)*T4_X(:,:,kd,1)
                    U_corr(:,:,kd,2) = n_frac(kd,ld)/mn(ld)*T1_X(:,:,kd,1)
                end do
                X%U_0(:,:,:,ld) = X%U_0(:,:,:,ld) + iu/mn(ld)*&
                    &U_corr(:,:,:,1)
                X%U_1(:,:,:,ld) = X%U_1(:,:,:,ld) + iu/mn(ld)*&
                    &U_corr(:,:,:,2)
            end if
            ! include order 3
            if (U_style.ge.3) then
                do kd = 1,grid_X%loc_n_r
                    U_corr(:,:,kd,1) = iu*n_frac(kd,ld)*&
                        &(-T5_X(:,:,kd,1)-iu*n_frac(kd,ld)*T6_X(:,:,kd,1))
                end do
                X%U_0(:,:,:,ld) = X%U_0(:,:,:,ld) + (iu/mn(ld))**2*&
                    &U_corr(:,:,:,1)
            end if
            if (U_style.ge.4) then
                call writo('WARNING: The geodesic perturbation U is &
                    &implemented up to order 3')
            end if
        end do
        
        ! calculate DU_0 and DU_1
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! calculate order 1 of DU_0 and DU_1
                if (U_style.ge.1) then
                    X%DU_0(:,:,:,ld) = -T2_X(:,:,:,2)
                    X%DU_1(:,:,:,ld) = 0._dp
                end if
                ! include order 2
                if (U_style.ge.2) then
                    do kd = 1,grid_X%loc_n_r
                        U_corr(:,:,kd,1) = T3_X(:,:,kd,2) + &
                            &iu*n_frac(kd,ld)*T4_X(:,:,kd,2)
                        U_corr(:,:,kd,2) = n_frac(kd,ld)/mn(ld)*T1_X(:,:,kd,2)
                    end do
                    X%DU_0(:,:,:,ld) = X%DU_0(:,:,:,ld) + iu/mn(ld)*&
                        &U_corr(:,:,:,1)
                    X%DU_1(:,:,:,ld) = X%DU_1(:,:,:,ld) + iu/mn(ld)*&
                        &U_corr(:,:,:,2)
                end if
                ! include order 3
                if (U_style.ge.3) then
                    do kd = 1,grid_X%loc_n_r
                        U_corr(:,:,kd,2) = iu*n_frac(kd,ld)*&
                            &(-T5_X(:,:,kd,2)-iu*n_frac(kd,ld)*T6_X(:,:,kd,2))
                    end do
                    X%DU_0(:,:,:,ld) = X%DU_0(:,:,:,ld) + (iu/mn(ld))**2*&
                        &U_corr(:,:,:,1)
                end if
                if (U_style.ge.4) then
                    call writo('WARNING: The geodesic perturbation U is &
                        &implemented up to order 3')
                end if
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
                ! reset pointers
                if (use_pol_flux_F) then
                    ang_par_F => grid_X%theta_F
                else
                    ang_par_F => grid_X%zeta_F
                end if
                ! set up helper variable for derivative
                allocate(D3var(grid_X%n(1)))
                ! numerically derive U_0 and U_1
                do ld = 1,X%n_mod
                    do kd = 1,grid_X%loc_n_r
                        ! calculate DX%U_0
                        do jd = 1,grid_X%n(2)
                            ierr = calc_deriv(X%U_0(:,jd,kd,ld),D3var,&
                                &ang_par_F(:,jd,kd),1,norm_disc_prec_X+1)       ! higher precision because previous derivative
                            CHCKERR('')
                            X%DU_0(:,jd,kd,ld) = D3var + &
                                &iu*n_frac(kd,ld)*X%U_0(:,jd,kd,ld)
                        end do
                        ! calculate DX%U_1
                        do jd = 1,grid_X%n(2)
                            ierr = calc_deriv(X%U_1(:,jd,kd,ld),D3var,&
                                &ang_par_F(:,jd,kd),1,norm_disc_prec_X+1)
                            CHCKERR('')
                            X%DU_1(:,jd,kd,ld) = D3var + &
                                &iu*n_frac(kd,ld)*X%U_1(:,jd,kd,ld)
                        end do
                    end do
                end do
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! clean up
        nullify(ang_par_F,ang_geo_F)
        nullify(J,D3J)
        nullify(g13,D3g13)
        nullify(g23,D3g23)
        nullify(g33,D3g33)
        nullify(h12,D3h12)
        nullify(h22,D1h22,D3h22,D13h22,D33h22)
        nullify(h23,D1h23,D3h23,D13h23,D33h23)
#if ldebug
    contains
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
            allocate(DU_0(grid_X%n(1),grid_X%n(2),grid_X%n(3)))
            allocate(DU_1(grid_X%n(1),grid_X%n(2),grid_X%n(3)))
            
            ! loop over all modes
            do jd = 1,X%n_mod
                ! loop over all normal points
                do kd = 1,grid_X%n(3)
                    ! derive numerically
                    do id = 1,grid_X%n(2)
                        ierr = calc_deriv(X%U_0(:,id,kd,jd),DU_0(:,id,kd),&
                            &ang_par_F(:,id,kd),1,norm_disc_prec_X+1)
                        CHCKERR('')
                        ierr = calc_deriv(X%U_1(:,id,kd,jd),DU_1(:,id,kd),&
                            &ang_par_F(:,id,kd),1,norm_disc_prec_X+1)
                        CHCKERR('')
                    end do
                    
                    ! add the second part
                    DU_0(:,:,kd) = DU_0(:,:,kd) + &
                        &iu*n_frac(kd,jd)*X%U_0(:,:,kd,jd)
                    DU_1(:,:,kd) = DU_1(:,:,kd) + &
                        &iu*n_frac(kd,jd)*X%U_1(:,:,kd,jd)
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
    ! eq loc_n_r values.
    ! Like  in calc_U,  use is  made of  variables  Ti that  are set  up in  the
    ! equilibrium grid and are afterwards  converted to Ti_X in the perturbation
    ! grid, which needs to have the  same angular coordinates as the equilibrium
    ! grid:
    !   - T1 = h22/mu_0 g33
    !   - T2 = JS + mu_0 sigma g33/J h22
    !   - T3 = sigma/J T2
    !   - T4 = 1/mu_0 J^2 h22
    !   - T5 = 2 p' kappa_n
    ! The interpolated Ti_X are then used to calculate PV:
    !   - PV_0 = T1(DU_k^0* - T2)(DU_m^0 - T2) + (nq-m)(nq-k)T4 - T5
    !   - PV_1 = T1(DU_k^0* - T2) DU_m^1
    !   - PV_2 = T1 DU_k^1* DU_m^1
    ! which  is valid for  poloidal Flux coordinates and  where (nq-m) is  to be
    ! replaced by (n-iota m) for toroidal Flux coordinates.
    ! (see [ADD REF] for details)
    integer function calc_PV(grid_eq,grid_X,eq,met,X_a,X_b,X,lim_sec_X) &
        &result(ierr)
        use num_vars, only: use_pol_flux_F
        use eq_vars, only: vac_perm
        use utilities, only: c
        use X_vars, only: is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        use grid_ops, only: get_norm_interp_data
        
        character(*), parameter :: rout_name = 'calc_PV'
        
        ! use input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(in) :: X_a, X_b                                  ! vectorial perturbation variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: loc_r_eq(:)                                    ! unrounded index of eq variables in X grid
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        integer :: i_lo, i_hi                                                   ! upper and lower index
        
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:)                                         ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        ! helper variables
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        ! PV factors
        real(dp), allocatable :: T1(:,:,:)                                      ! h22/mu_0 g33
        real(dp), allocatable :: T2(:,:,:)                                      ! JS + mu_0 sigma g33/J h22
        real(dp), allocatable :: T3(:,:,:)                                      ! sigma/J T2
        real(dp), allocatable :: T4(:,:,:)                                      ! 1/mu_0 J^2 h22
        real(dp), allocatable :: T5(:,:,:)                                      ! 2 p' kappa_n
        real(dp), allocatable :: T1_X(:,:,:)                                    ! T1 in X grid and par. deriv.
        real(dp), allocatable :: T2_X(:,:,:)                                    ! T2 in X grid and par. deriv.
        real(dp), allocatable :: T3_X(:,:,:)                                    ! T3 in X grid and par. deriv.
        real(dp), allocatable :: T4_X(:,:,:)                                    ! T4 in X grid and par. deriv.
        real(dp), allocatable :: T5_X(:,:,:)                                    ! T5 in X grid and par. deriv.
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        ! helper variables
        allocate(fac_n(grid_X%loc_n_r),fac_m(grid_X%loc_n_r))
        ! PV factors
        allocate(T1(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T2(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T4(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T5(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T1_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(T2_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(T3_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(T4_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(T5_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        
        ! set pointers
        J => met%jac_FD(:,:,:,0,0,0)
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up PV factors in eq grid
        T1 = h22/(vac_perm*g33)
        T2 = J*eq%S + vac_perm*eq%sigma*g33/(J*h22)
        T3 = eq%sigma/J*T2
        T4 = 1._dp/(vac_perm*J**2*h22)
        do kd = 1,grid_eq%loc_n_r
            T5(:,:,kd) = 2*eq%pres_FD(kd,1)*eq%kappa_n(:,:,kd)
        end do
        
        ! get normal interpolation factors
        ierr = get_norm_interp_data(grid_eq,grid_X,loc_r_eq)
        CHCKERR('')
        
        ! interpolate helper variables and PV factors
        do kd = 1,grid_X%loc_n_r
            ! set lower and upper index
            i_lo = floor(loc_r_eq(kd))
            i_hi = ceiling(loc_r_eq(kd))
            
            if (use_pol_flux_F) then
                fac_n(kd) = eq%q_saf_FD(i_lo,0) + (loc_r_eq(kd)-i_lo)*&
                    &(eq%q_saf_FD(i_hi,0)-eq%q_saf_FD(i_lo,0))
                fac_m(kd) = 1.0_dp
            else
                fac_n(kd) = 1.0_dp
                fac_m(kd) = eq%rot_t_FD(i_lo,0) + (loc_r_eq(kd)-i_lo)*&
                    &(eq%rot_t_FD(i_hi,0)-eq%rot_t_FD(i_lo,0))
            end if
            T1_X(:,:,kd) = T1(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T1(:,:,i_hi)-T1(:,:,i_lo))
            T2_X(:,:,kd) = T2(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T2(:,:,i_hi)-T2(:,:,i_lo))
            T3_X(:,:,kd) = T3(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T3(:,:,i_hi)-T3(:,:,i_lo))
            T4_X(:,:,kd) = T4(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T4(:,:,i_hi)-T4(:,:,i_lo))
            T5_X(:,:,kd) = T5(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T5(:,:,i_hi)-T5(:,:,i_lo))
        end do
        
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
                    do kd = 1,grid_X%loc_n_r
                        X%PV_0(:,:,kd,c_loc(1)) = T1_X(:,:,kd)*&
                            &(X_b%DU_0(:,:,kd,m) - T2_X(:,:,kd)) * &
                            &(conjg(X_a%DU_0(:,:,kd,k)) - T2_X(:,:,kd)) - &
                            &T3_X(:,:,kd) - T5_X(:,:,kd) + T4_X(:,:,kd)*&
                            &(X_b%n(m)*fac_n(kd)-X_b%m(m)*fac_m(kd))*&
                            &(X_a%n(k)*fac_n(kd)-X_a%m(k)*fac_m(kd))
                    end do
                end if
                
                ! calculate PV_1
                if (calc_this(2)) then
                    X%PV_1(:,:,:,c_loc(2)) = T1_X * X_b%DU_1(:,:,:,m) * &
                        &(conjg(X_a%DU_0(:,:,:,k))-T2_X)
                end if
                
                ! calculate PV_2
                if (calc_this(1)) then
                    X%PV_2(:,:,:,c_loc(1)) = T1_X * X_b%DU_1(:,:,:,m) * &
                        &conjg(X_a%DU_1(:,:,:,k))
                end if
            end do
        end do
        
        ! clean up
        nullify(J)
        nullify(g33)
        nullify(h22)
    end function calc_PV
    
    ! calculate  ~KV_(k,m)^i  (pol.  flux)  or ~KV_(l,n)^i  (tor.  flux) at  all
    ! eq loc_n_r values.
    ! Like  in calc_U,  use is  made of  variables  Ti that  are set  up in  the
    ! equilibrium grid and are afterwards  converted to Ti_X in the perturbation
    ! grid, which needs to have the  same angular coordinates as the equilibrium
    ! grid:
    !   - T1 = rho J^2 h22/g33
    !   - T2 = rho/h22
    ! The interpolated Ti_X are then used to calculate KV:
    !   - KV_0 = T1 U_k^0* U_m^0 + T2
    !   - KV_1 = T1 U_k^0* U_m^1
    !   - KV_2 = T1 U_k^1* U_m^1
    ! which  is valid for  poloidal Flux coordinates and  where (nq-m) is  to be
    ! replaced by (n-iota m) for toroidal Flux coordinates.
    ! (see [ADD REF] for details)
    integer function calc_KV(grid_eq,grid_X,eq,met,X_a,X_b,X,lim_sec_X) &
        &result(ierr)
        use num_vars, only: norm_style
        use utilities, only: c
        use X_vars, only: is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        use grid_ops, only: get_norm_interp_data
        
        character(*), parameter :: rout_name = 'calc_KV'
        
        ! use input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_1_type), intent(in) :: X_a, X_b                                  ! vectorial perturbation variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_mod_tot                                                    ! total nr. of modes
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: loc_r_eq(:)                                    ! unrounded index of eq variables in X grid
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        integer :: i_lo, i_hi                                                   ! upper and lower index
        
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:)                                         ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        ! KV factors
        real(dp), allocatable :: T1(:,:,:)                                      ! h22/mu_0 J g33
        real(dp), allocatable :: T2(:,:,:)                                      ! JS + mu_0 sigma g33/h22
        real(dp), allocatable :: T1_X(:,:,:)                                    ! T1 in X grid and par. deriv.
        real(dp), allocatable :: T2_X(:,:,:)                                    ! T2 in X grid and par. deriv.
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        ! KV factors
        allocate(T1(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T2(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T1_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(T2_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        
        ! set pointers
        J => met%jac_FD(:,:,:,0,0,0)
        g33 => met%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        h22 => met%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up KV factors in eq grid
        do kd = 1,grid_eq%loc_n_r
            T1(:,:,kd) = eq%rho(kd)*J(:,:,kd)**2*h22(:,:,kd)/g33(:,:,kd)
            T2(:,:,kd) = eq%rho(kd)/h22(:,:,kd)
        end do
        
        ! get normal interpolation factors
        ierr = get_norm_interp_data(grid_eq,grid_X,loc_r_eq)
        CHCKERR('')
        
        ! interpolate KV factors
        do kd = 1,grid_X%loc_n_r
            ! set lower and upper index
            i_lo = floor(loc_r_eq(kd))
            i_hi = ceiling(loc_r_eq(kd))
            
            T1_X(:,:,kd) = T1(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T1(:,:,i_hi)-T1(:,:,i_lo))
            T2_X(:,:,kd) = T2(:,:,i_lo) + (loc_r_eq(kd)-i_lo)*&
                &(T2(:,:,i_hi)-T2(:,:,i_lo))
        end do
        
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
                
                select case (norm_style)
                    case (1)                                                    ! normalization of full perpendicular component
                        ! calculate KV_0
                        if (calc_this(1)) then
                            X%KV_0(:,:,:,c_loc(1)) = T2_X + T1_X * &
                                &X_b%u_0(:,:,:,m) * conjg(X_a%u_0(:,:,:,k))
                        end if
                        
                        ! calculate KV_1
                        if (calc_thiS(2)) then
                            X%KV_1(:,:,:,c_loc(2)) = T1_X * &
                                &X_b%u_1(:,:,:,m) * conjg(X_a%u_0(:,:,:,k))
                        end if
                        
                        ! calculate kv_2
                        if (calc_this(1)) then
                            X%kv_2(:,:,:,c_loc(1)) = T1_X * &
                                &X_b%U_1(:,:,:,m) * conjg(X_a%U_1(:,:,:,k))
                        end if
                    case (2)                                                    ! normalization of only normal component
                        ! calculate KV_0
                        if (calc_this(1)) then
                            X%KV_0(:,:,:,c_loc(1)) = T2_X
                        end if
                        
                        ! calculate KV_1
                        if (calc_this(2)) then
                            X%KV_1(:,:,:,c_loc(2)) = 0._dp
                        end if
                        
                        ! calculate KV_2
                        if (calc_this(1)) then
                            X%KV_2(:,:,:,c_loc(1)) = 0._dp
                        end if
                    case default
                        err_msg = 'No normalization style associated with '//&
                            &trim(i2str(norm_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            end do
        end do
        
        ! clean up
        nullify(J)
        nullify(g33)
        nullify(h22)
    end function calc_KV
    
    ! Calculate the  magnetic integrals  from PV_i and  KV_i. All  the variables
    ! should thus be field-line oriented.
    integer function calc_magn_ints(grid_eq,grid_X,met,X,lim_sec_X) result(ierr)
        use num_vars, only: use_pol_flux_F
        use X_vars, only: is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        use utilities, only: c
        use grid_ops, only: get_norm_interp_data
        
        character(*), parameter :: rout_name = 'calc_magn_ints'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(met_type), intent(in) :: met                                       ! metric variables
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        integer :: k, m                                                         ! counters
        integer :: id, kd                                                       ! counters
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        complex(dp), allocatable :: J_exp_ang(:,:,:)                            ! J * exponential of Flux parallel angle
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle in flux coordinates
        real(dp), allocatable :: loc_r_eq(:)                                    ! unrounded index of eq variables in X grid
        integer :: i_lo, i_hi                                                   ! upper and lower index
        
        ! jacobian
        real(dp), allocatable :: J(:,:,:)                                       ! jac
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculating field-line averages')
        call lvl_ud(1)
        
        ! allocate variables
        allocate(J(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        allocate(J_exp_ang(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        
        ! get normal interpolation factors
        ierr = get_norm_interp_data(grid_eq,grid_X,loc_r_eq)
        CHCKERR('')
        
        ! interpolate equilibrium variables
        do kd = 1,grid_X%loc_n_r
            ! set lower and upper index
            i_lo = floor(loc_r_eq(kd))
            i_hi = ceiling(loc_r_eq(kd))
            
            
            ! jacobian
            J(:,:,kd) = met%jac_FD(:,:,i_lo,0,0,0) + (loc_r_eq(kd)-i_lo)*&
                &(met%jac_FD(:,:,i_hi,0,0,0)-met%jac_FD(:,:,i_lo,0,0,0))
        end do
        
        ! set nr. of modes
        n_mod_tot = (max_m_X-min_m_X+1)*(max_n_X-min_n_X+1)
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_F) then
            ang_par_F => grid_X%theta_F
        else
            ang_par_F => grid_X%zeta_F
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
                    J_exp_ang = J*exp(iu*(X%m_1(k)-X%m_2(m))*ang_par_F)
                else
                    J_exp_ang = J*exp(-iu*(X%n_1(k)-X%n_2(m))*ang_par_F)
                end if
                
                ! parallel integration loop
                do id = 2,grid_X%n(1)
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
    end function calc_magn_ints
    
    ! Print either  vectorial or tensorial perturbation quantities  of a certain
    ! order to an output file:
    !   - vectorial:    U, DU
    !   - tensorial:    PV_int, KV_int
    !     (the non-integrated variables are heavy and not requested)
    ! Note: Flux coordinates used as normal coordinates
    integer function print_output_X_1(grid,X) result(ierr)                      ! vectorial version
        use num_vars, only: PB3D_name
        use rich, only: rich_info_short
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        use X_vars, only: get_suffix
        
        character(*), parameter :: rout_name = 'print_output_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! perturbation grid variables
        type(X_1_type), intent(in) :: X                                         ! vectorial perturbation variables 
        
        ! local variables
        type(var_1D_type), allocatable, target :: X_1D(:)                       ! 1D equivalent of X variables
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
        
        do jd = 1,X%n_mod
            ! RE_U_0
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_U_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_0(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%U_0(:,:,:,jd)),&
                &[size(X%U_0(:,:,:,jd))])
            
            ! IM_U_0
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_U_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_0(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%U_0(:,:,:,jd)),&
                &[size(X%U_0(:,:,:,jd))])
            
            ! RE_U_1
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_U_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_1(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%U_1(:,:,:,jd)),&
                &[size(X%U_1(:,:,:,jd))])
            
            ! IM_U_1
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_U_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%U_1(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%U_1(:,:,:,jd)),&
                &[size(X%U_1(:,:,:,jd))])
            
            ! RE_DU_0
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_DU_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_0(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%DU_0(:,:,:,jd)),&
                &[size(X%DU_0(:,:,:,jd))])
            
            ! IM_DU_0
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_DU_0_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_0(:,:,:,jd))))
            X_1D_loc%p = reshape(imagpart(X%DU_0(:,:,:,jd)),&
                &[size(X%DU_0(:,:,:,jd))])
            
            ! RE_DU_1
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'RE_DU_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
            X_1D_loc%loc_i_min = X_1D_loc%tot_i_min
            X_1D_loc%loc_i_max = X_1D_loc%tot_i_max
            allocate(X_1D_loc%p(size(X%DU_1(:,:,:,jd))))
            X_1D_loc%p = reshape(realpart(X%DU_1(:,:,:,jd)),&
                &[size(X%DU_1(:,:,:,jd))])
            
            ! IM_DU_1
            X_1D_loc => X_1D(id); id = id+1
            X_1D_loc%var_name = 'IM_DU_1_'//trim(get_suffix(X,jd))
            allocate(X_1D_loc%tot_i_min(3),X_1D_loc%tot_i_max(3))
            allocate(X_1D_loc%loc_i_min(3),X_1D_loc%loc_i_max(3))
            X_1D_loc%tot_i_min = [1,1,1]
            X_1D_loc%tot_i_max = grid%n
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
        ierr = print_HDF5_arrs(X_1D,PB3D_name,'X_1'//trim(rich_info_short()))
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        
        ! clean up
        nullify(X_1D_loc)
        
        ! user output
        call lvl_ud(-1)
        call writo('Vectorial perturbation variables written to output')
    end function print_output_X_1
    integer function print_output_X_2(grid,X) result(ierr)                      ! tensorial version
        use num_vars, only: PB3D_name, use_pol_flux_F
        use rich, only: rich_info_short
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: var_1D_type
        use X_vars, only: get_suffix, is_necessary_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'print_output_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! perturbation grid variables
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables 
        
        ! local variables
        integer :: n_mod_tot                                                    ! total nr. of modes
        type(var_1D_type), allocatable, target :: X_1D(:)                       ! 1D equivalent of X variables
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
                    X_1D_loc%tot_i_max = grid%n(2:3)
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
        ierr = print_HDF5_arrs(X_1D,PB3D_name,'X_2'//trim(rich_info_short()))
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
