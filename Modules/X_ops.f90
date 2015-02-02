!------------------------------------------------------------------------------!
!   Operations considering perturbation quantities                             !
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    use messages, only: lvl_ud, writo, print_ar_2
    use output_ops, only: print_GP_2D, print_GP_3D, draw_GP
    use str_ops, only: i2str, r2strt, r2str

    implicit none
    private
    public init_X, init_m, prepare_X, solve_EV_system, plot_X_vec, calc_PV, &
        &calc_KV, calc_U, calc_extra, calc_V_int
    
    !! for testing:
    !character(len=max_str_ln), allocatable :: var_names(:)                      ! names of variables
    !character(len=max_str_ln), allocatable :: file_name                         ! name of file
    !integer :: tot_dim(4), grp_dim(4), grp_offset(4)                            ! total dim., group dim. and group offset

contains
    ! initialize the variable m and check and/or plot it
    integer function init_m() result(ierr)
        use num_vars, only: plot_jq, use_pol_flux_X, glb_rank
        use X_vars, only: n_X, m_X, min_m_X, max_m_X, min_n_X, max_n_X, size_X
        
        character(*), parameter :: rout_name = 'init_m'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set size_X
        if (use_pol_flux_X) then
            size_X = max_m_X - min_m_X + 1
        else
            size_X = max_n_X - min_n_X + 1
        end if
        
        ! setup n_X and m_X
        allocate(n_X(size_X),m_X(size_X))
        if (use_pol_flux_X) then
            n_X = min_n_X
            m_X = [(id, id = min_m_X, max_m_X)]
        else
            n_X = [(id, id = min_n_X, max_n_X)]
            m_X = min_m_X
        end if
        
        ! plot resonances if requested
        if (plot_jq .and. glb_rank.eq.0) then
            ierr = resonance_plot()
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! check whether m and n are sensible
        ierr = check_m_and_n()
        CHCKERR('')
    contains
        ! checks whether nq-m  << n or n-iotam << m is satisfied  in some of the
        ! plasma
        integer function check_m_and_n() result(ierr)
            use eq_vars, only: q_saf_E_full, rot_t_E_full                       ! full variables are NOT normalized
            use num_vars, only: glb_rank, use_pol_flux_X, eq_style
            
            character(*), parameter :: rout_name = 'check_m_and_n'
            
            ! local variables
            integer :: id                                                       ! counter
            real(dp) :: tol = 0.2                                               ! tolerance for being out of range of q or iota values
            real(dp) :: min_jq, max_jq                                          ! min. and max. values of q or iota
            integer :: pmone                                                    ! plus or minus one
            
            ! initialize ierr
            ierr = 0
            
            if (glb_rank.eq.0) then
                call writo('Checking mode numbers')
                ! set up plus minus one
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        pmone = -1                                              ! conversion VMEC LH -> RH coord. system
                    case (2)                                                    ! HELENA
                        pmone = 1
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! set min_jq and max_jq in flux coordinate system
                if (use_pol_flux_X) then
                    min_jq = minval(pmone*q_saf_E_full(:,0))
                    max_jq = maxval(pmone*q_saf_E_full(:,0))
                else
                    min_jq = minval(pmone*rot_t_E_full(:,0))
                    max_jq = maxval(pmone*rot_t_E_full(:,0))
                end if
                
                ! multiply by plus minus one and tolerance
                min_jq = min_jq - tol*abs(min_jq)
                max_jq = max_jq + tol*abs(max_jq)
                
                ! for every mode (n,m) check whether  m/n is inside the range of
                ! q values or n/m inside the range of iota values
                do id = 1, size_X
                    if (use_pol_flux_X) then
                        if (m_X(id)*1.0/n_X(id) .lt.min_jq .or. &
                            &m_X(id)*1.0/n_X(id) .gt.max_jq) then
                            call writo('for (n,m) = ('//trim(i2str(n_X(id)))//&
                                &','//trim(i2str(m_X(id)))//'), the ratio &
                                &|n q - m|/n is never << 1')
                            ierr = 1
                            err_msg = 'Choose m and n so that |n q - m|/n << 1'
                            CHCKERR(err_msg)
                        end if
                    else
                        if (n_X(id)*1.0/m_X(id) .lt.min_jq .or. &
                            &n_X(id)*1.0/m_X(id) .gt.max_jq) then
                            call writo('for (n,m) = ('//trim(i2str(n_X(id)))//&
                                &','//trim(i2str(m_X(id)))//'), the ratio &
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
        end function check_m_and_n
    end function init_m
    
    ! plot  q-profile  or iota-profile  in  flux coordinates  with nq-m  = 0  or
    ! n-iotam = 0 indicate if requested
    integer function resonance_plot() result(ierr)
        use num_vars, only: use_pol_flux_X, eq_style, output_style
        use utilities, only: calc_zero_NR, interp_fun
        use eq_vars, only: q_saf_E_full, rot_t_E_full, flux_p_E_full, &
            &flux_t_E_full, n_r_eq                                              ! full variables are NOT normalized
        use X_vars, only: size_X, m_X, n_X
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
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
        real(dp), allocatable :: q_saf(:,:), rot_t(:,:)                         ! saf. fac., rot. transf. in Equilibrium coords.
        real(dp), allocatable :: flux_p(:), flux_t(:)                           ! pol. flux, tor. flux in Equilibrium coords.
        integer :: pmone                                                        ! plus or minus one
        
        ! initialize ierr
        ierr = 0
        
        ! set up plus minus one
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                pmone = -1                                                      ! conversion of VMEC LH to Flux RH
            case (2)                                                            ! HELENA
                pmone = 1                                                       ! no conversion of HELENA HH to Flux RH
            case default
                call writo('No equilibrium style associated with '//&
                    &trim(i2str(eq_style)))
                call writo('Aborting')
                return
        end select
        
        if (use_pol_flux_X) then
            allocate(flux_p(size(flux_p_E_full,1)))
            flux_p = flux_p_E_full(:,0)
            allocate(q_saf(size(q_saf_E_full,1),0:2))
            q_saf = q_saf_E_full(:,0:2)
        else
            allocate(flux_t(size(flux_t_E_full,1)))
            flux_t = flux_t_E_full(:,0)
            allocate(rot_t(size(rot_t_E_full,1),0:2))
            rot_t = rot_t_E_full(:,0:2)
        end if
        
        if (use_pol_flux_X) then
            call writo('Plotting safety factor q and resonant surfaces &
                &q = m/n')
        else
            call writo('Plotting rotational transform iota and resonant &
                &surfaces iota = n/m')
        end if
        call lvl_ud(1)
        
        call writo('calculating resonant surfaces')
        call lvl_ud(1)
        
        ! initialize variables
        allocate(x_vars(n_r_eq,size_X+1)); x_vars = 0
        allocate(y_vars(n_r_eq,size_X+1)); y_vars = 0
        allocate(jq_for_function(n_r_eq,0:2))
        
        if (use_pol_flux_X) then
            x_vars(:,1) = flux_p/abs(flux_p(n_r_eq))
            y_vars(:,1) = pmone*q_saf(:,0)
            jq_for_function = q_saf
        else
            x_vars(:,1) = flux_t/abs(flux_t(n_r_eq))
            y_vars(:,1) = pmone*rot_t(:,0)
            jq_for_function = rot_t
        end if
        
        kd = 2
        do jd = 1, size_X
            ! find place  where q  = m/n or  iota = n/m  in VMEC  coordinates by
            ! solving q-m/n = 0 or iota-n/m=0, using the functin jq_fun
            call lvl_ud(1)
            
            ! set up mnfrac_for_function
            if (use_pol_flux_X) then
                mnfrac_for_function = pmone*m_X(jd)*1.0_dp/n_X(jd)
            else
                mnfrac_for_function = pmone*n_X(jd)*1.0_dp/m_X(jd)
            end if
            
            ! calculate zero using Newton-Rhapson
            istat = calc_zero_NR(jq_solution,jq_fun,jq_dfun,1.0_dp)
            call lvl_ud(-1)
            
            ! intercept error
            if (istat.ne.0) then
                call writo('Error intercepted: Couldn''t find resonating &
                    &surface for (n,m) = ('//trim(i2str(n_X(jd)))//','//&
                    &trim(i2str(m_X(jd)))//')')
            else if (jq_solution.lt.0.0_dp) then
                call writo('Mode (n,m) = ('//trim(i2str(n_X(jd)))//','//&
                    &trim(i2str(m_X(jd)))//') does not resonate in plasma')
            else
                if (jq_solution.gt.1.0_dp) then
                    call writo('Mode (n,m) = ('//trim(i2str(n_X(jd)))//','&
                        &//trim(i2str(m_X(jd)))//') does not resonate &
                        &in plasma')
                    if (use_pol_flux_X) then
                        y_vars(n_r_eq,kd) = q_saf(n_r_eq,0)
                    else
                        y_vars(n_r_eq,kd) = rot_t(n_r_eq,0)
                    end if
                else
                    ! convert solution to flux coordinates if GNUPlot output
                    if (output_style.eq.1) then
                        if (use_pol_flux_X) then
                            istat = interp_fun(jq_solution_transf,&
                                &flux_p/abs(flux_p(n_r_eq)),jq_solution)
                        else
                            istat = interp_fun(jq_solution_transf,&
                                &flux_t/abs(flux_t(n_r_eq)),jq_solution)
                        end if
                        x_vars(:,kd) = jq_solution_transf
                        call writo('Mode (n,m) = ('//trim(i2str(n_X(jd)))//&
                            &','//trim(i2str(m_X(jd)))//') resonates in &
                            &plasma at normalized flux surface '//&
                            &trim(r2str(jq_solution_transf)))
                    else                                                        ! HDF5 output needs equilibrium tabulated values
                        x_vars(:,kd) = jq_solution
                    end if
                    ! the y axis is always in perturbation grid
                    if (use_pol_flux_X) then
                        y_vars(n_r_eq,kd) = m_X(jd)*1.0_dp/n_X(jd)
                    else
                        y_vars(n_r_eq,kd) = n_X(jd)*1.0_dp/m_X(jd)
                    end if
                end if
                kd = kd + 1
            end if
        end do
        
        ! check status
        if (istat.ne.0) then
            call writo('WARNING: Failed to produce plot')
            return
        end if
        
        call lvl_ud(-1)
        
        call writo('Plotting results')
        if (use_pol_flux_X) then
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
        
        ! nullify pointers
        if (use_pol_flux_X) then
            deallocate(flux_p,q_saf)
        else
            deallocate(flux_t,rot_t)
        end if
        
        call lvl_ud(-1)
        call writo('Done plotting')
    contains
        ! Returns q-m/n  or iota-n/m in  Equilibrium coordinates, used  to solve
        ! for q = m/n or iota = n/m.
        real(dp) function jq_fun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(n_r_eq)                                           ! so interp_fun can be used
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.0) then                                                   ! point requested lower than 0
                ! extrapolate variable from 1
                res = jq_for_function(1,0) - mnfrac_for_function + &
                    &jq_for_function(1,1)*(pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from n_r_eq
                res = jq_for_function(n_r_eq,0) - mnfrac_for_function + &
                    &jq_for_function(n_r_eq,1)*(pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using interp_fun
                varin = jq_for_function(:,0) - mnfrac_for_function
                istat = interp_fun(res,varin,pt)
            end if
        end function jq_fun
        
        ! returns d(q-m/n)/dr Equilibrium coordinates, used to solve for q = m/n
        ! WARNING: This  routine requires that jq_for_function's  derivatives be
        ! calculated up to order 2. This is NOT checked!
        real(dp) function jq_dfun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(n_r_eq)                                           ! so interp_fun can be used
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.0) then                                                   ! point requested lower than 0
                ! extrapolate variable from 1
                res = jq_for_function(1,1) + jq_for_function(1,2)*(pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from n_r_eq
                res = jq_for_function(n_r_eq,1) + jq_for_function(n_r_eq,2)*&
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
            use grid_ops, only: calc_XYZ_grid, calc_eqd_grid
            
            character(*), parameter :: rout_name = 'resonance_plot_HDF5'
            
            ! local variables
            integer :: id                                                       ! counters
            real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)        ! pol. and tor. angle of plot
            real(dp), allocatable :: X_plot(:,:,:,:), Y_plot(:,:,:,:), &
                &Z_plot(:,:,:,:)                                                ! X, Y and Z of plot of all surfaces
            real(dp), allocatable :: X_plot_ind(:,:,:), Y_plot_ind(:,:,:), &
                &Z_plot_ind(:,:,:)                                              ! X, Y and Z of plots of individual surfaces
            integer :: plot_dim(4)                                              ! plot dimensions (total = group because only group masters)
            real(dp), allocatable :: vars(:,:,:,:)                              ! variable to plot
            
            ! initialize ierr
            ierr = 0
                
            ! set up pol. and tor. angle for plot
            allocate(theta_plot(n_theta_plot,n_zeta_plot,1))
            ierr = calc_eqd_grid(theta_plot(:,1,1),n_theta_plot,0._dp,2*pi)
            CHCKERR('')
            do id = 2,n_zeta_plot
                theta_plot(:,id,1) = theta_plot(:,1,1)
            end do
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,1))
            ierr = calc_eqd_grid(zeta_plot(1,:,1),n_zeta_plot,0._dp,2*pi)
            CHCKERR('')
            do id = 2,n_theta_plot
                zeta_plot(id,:,1) = zeta_plot(1,:,1)
            end do
            
            ! set up vars
            allocate(vars(n_theta_plot,n_zeta_plot,1,size_X))
            do id = 1,size_X
                vars(:,:,:,id) = y_vars(n_r_eq,id+1)
            end do
            
            ! set dimensions
            plot_dim = [n_theta_plot,n_zeta_plot,1,size_X]
            
            ! set up plot X, Y and Z
            allocate(X_plot(n_theta_plot,n_zeta_plot,1,size_X))
            allocate(Y_plot(n_theta_plot,n_zeta_plot,1,size_X))
            allocate(Z_plot(n_theta_plot,n_zeta_plot,1,size_X))
            
            ! loop over all resonant surfaces to calculate X, Y and Z values
            do id = 1,size_X
                ! calculate the X, Y and Z values for this flux surface
                ierr = calc_XYZ_grid([x_vars(n_r_eq,id+1)],theta_plot,&
                    &zeta_plot,X_plot_ind,Y_plot_ind,Z_plot_ind)
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
        end function resonance_plot_HDF5
    end function resonance_plot
    
    ! initialize the perturbation (X) variables
    ! Note: The variables U_X_i, DU_X_i have the following indices:
    !   1:n_par     variable along the magnetic field line (in theta or zeta)
    !   1:grp_n_r_X variable in the normal direction of the perturbation
    !   1:size_X    poloidal or toroidal mode number
    ! and variables PVi,  KVi have an extra  index with the same  range as the
    ! third, but in this case the third index refers to the complex conjugate X*
    ! while the fourth refers to X
    ! Finally, the variables PV_int and KV_int have the following indices:
    !   1:size_X    poloidal or toroidal mode number of X*
    !   1:size_X    poloidal or toroidal mode number of X
    !   1:grp_n_r_X variable in the normal direction of the perturbation
    !   1:3         order of differential in [ADD REF]
    ! ¡¡¡ THIS SHOULD BE  CHANGED TO TAKE INTO ACCOUNT THAT  U_X AND DU_X ARE
    ! HERMITIAN SO ONLY  M*(M+1)/2 ELEMENTS ARE NEEDED INSTEAD  OF M^2. HOWEVER,
    ! TAKE INTO ACCOUNT AS WELL THE MPI STORAGE MATRIX FORMATS !!!
    subroutine init_X
        use eq_vars, only: grp_n_r_eq
        use num_vars, only: use_pol_flux_X
        use X_vars, only: U_X_0, U_X_1, DU_X_0, DU_X_1, PV0, PV1, PV2, KV0, &
            &KV1, KV2, exp_ang_par_F, PV_int, KV_int, size_X, n_par, ang_par_F
        
        ! local variables
        integer :: m, k                                                         ! counter
        
        ! U_X_0
        allocate(U_X_0(n_par,grp_n_r_eq,size_X))
        
        ! U_X_1
        allocate(U_X_1(n_par,grp_n_r_eq,size_X))
        
        ! DU_X_0
        allocate(DU_X_0(n_par,grp_n_r_eq,size_X))
        
        ! DU_X_1
        allocate(DU_X_1(n_par,grp_n_r_eq,size_X))
        
        ! PV0
        allocate(PV0(n_par,grp_n_r_eq,size_X,size_X))
        
        ! PV1
        allocate(PV1(n_par,grp_n_r_eq,size_X,size_X))
        
        ! PV2
        allocate(PV2(n_par,grp_n_r_eq,size_X,size_X))
        
        ! KV0
        allocate(KV0(n_par,grp_n_r_eq,size_X,size_X))
        
        ! KV1
        allocate(KV1(n_par,grp_n_r_eq,size_X,size_X))
        
        ! KV2
        allocate(KV2(n_par,grp_n_r_eq,size_X,size_X))
        
        ! exp_ang_par_F
        allocate(exp_ang_par_F(n_par,grp_n_r_eq,size_X,size_X))
        if (use_pol_flux_X) then
            do m = 1,size_X
                do k = 1,size_X
                    exp_ang_par_F(:,:,k,m) = exp(iu*(k-m)*ang_par_F)
                end do
            end do
        else
            do m = 1,size_X
                do k = 1,size_X
                    exp_ang_par_F(:,:,k,m) = exp(iu*(m-k)*ang_par_F)
                end do
            end do
        end if
        
        ! PV_int
        allocate(PV_int(size_X,size_X,grp_n_r_eq,3))
        
        ! KV_int
        allocate(KV_int(size_X,size_X,grp_n_r_eq,3))
    end subroutine init_X
    
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    integer function prepare_X() result(ierr)
        use X_vars, only: dealloc_X, &
            &PV0, PV1, PV2, KV0, KV1, KV2, PV_int, KV_int
        use num_vars, only: use_pol_flux_X
        
        character(*), parameter :: rout_name = 'prepare_X'
        
        ! local variables
        character(len=5) :: ang_par_F_name                                      ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up ang_par_F_name
        if (use_pol_flux_X) then
            ang_par_F_name = 'theta'
        else
            ang_par_F_name = 'zeta'
        end if
        
        call writo('Calculating table of field line averages &
            &<V e^i[(k-m)'//trim(ang_par_F_name)//']>')
        
        call lvl_ud(1)
        
        ! initialize the variables
        call writo('Initalizing variables...')
        call lvl_ud(1)
        call init_X
        call lvl_ud(-1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        ierr = calc_U()
        CHCKERR('')
        call lvl_ud(-1)
        
        ! calculate extra equilibrium quantities
        call writo('Calculating extra equilibrium quantities...')
        call lvl_ud(1)
        ierr = calc_extra()
        CHCKERR('')
        call lvl_ud(-1)
        
        ! Calculate PV0, PV1  and PV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating PV0, PV1 and PV2...')
        call lvl_ud(1)
        call calc_PV
        call lvl_ud(-1)
        
        ! Calculate KV0, KV1  and KV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV0, KV1 and KV2...')
        call lvl_ud(1)
        call calc_KV
        call lvl_ud(-1)
        
        ! Calculate PV_int = <PVi e^(k-m)ang_par_F>
        call writo('Taking field average of PV')
        call lvl_ud(1)
        call calc_V_int(PV0,PV_int(:,:,:,1))
        call calc_V_int(PV1,PV_int(:,:,:,2))
        call calc_V_int(PV2,PV_int(:,:,:,3))
        call lvl_ud(-1)
        !write(*,*) 'RE min PV = ', minval(realpart(PV_int(:,:,:,1))), &
            !&minval(realpart(PV_int(:,:,:,2))), minval(realpart(PV_int(:,:,:,3)))
        !write(*,*) 'IM min PV = ', minval(imagpart(PV_int(:,:,:,1))), &
            !&minval(imagpart(PV_int(:,:,:,2))), minval(imagpart(PV_int(:,:,:,3)))
        !write(*,*) 'RE max PV = ', maxval(realpart(PV_int(:,:,:,1))), &
            !&maxval(realpart(PV_int(:,:,:,2))), maxval(realpart(PV_int(:,:,:,3)))
        !write(*,*) 'IM max PV = ', maxval(imagpart(PV_int(:,:,:,1))), &
            !&maxval(imagpart(PV_int(:,:,:,2))), maxval(imagpart(PV_int(:,:,:,3)))
        
        ! Calculate KV_int = <KVi e^(k-m)ang_par_F>
        call writo('Taking field average of KV')
        call lvl_ud(1)
        call calc_V_int(KV0,KV_int(:,:,:,1))
        call calc_V_int(KV1,KV_int(:,:,:,2))
        call calc_V_int(KV2,KV_int(:,:,:,3))
        call lvl_ud(-1)
        !write(*,*) 'RE min KV = ', minval(realpart(KV_int(:,:,:,1))), &
            !&minval(realpart(KV_int(:,:,:,2))), minval(realpart(KV_int(:,:,:,3)))
        !write(*,*) 'IM min KV = ', minval(imagpart(KV_int(:,:,:,1))), &
            !&minval(imagpart(KV_int(:,:,:,2))), minval(imagpart(KV_int(:,:,:,3)))
        !write(*,*) 'RE max KV = ', maxval(realpart(KV_int(:,:,:,1))), &
            !&maxval(realpart(KV_int(:,:,:,2))), maxval(realpart(KV_int(:,:,:,3)))
        !write(*,*) 'IM max KV = ', maxval(imagpart(KV_int(:,:,:,1))), &
            !&maxval(imagpart(KV_int(:,:,:,2))), maxval(imagpart(KV_int(:,:,:,3)))
        
        ! deallocate equilibrium variables
        call writo('deallocating unused variables')
        call dealloc_X
        
        call lvl_ud(-1)
        
        call writo('Done calculating')
    end function prepare_X
    
    ! set-up and  solve the  EV system  by discretizing  the equations  in n_r_X
    ! normal points,  making use of  PV0, PV1 and  PV2, interpolated in  the n_r
    ! (equilibrium) values
    integer function solve_EV_system(use_guess,n_sol_found) result(ierr)
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use SLEPC_ops, only: solve_EV_system_SLEPC
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        logical, intent(in) :: use_guess                                        ! whether to use a guess or not
        integer, intent(inout) :: n_sol_found                                   ! how many solutions saved
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! SLEPC solver for EV problem
                ! solve the system
                ierr = solve_EV_system_SLEPC(use_guess,n_sol_found)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
    
    ! calculate rho from user input
    integer function calc_rho() result(ierr)
        use num_vars, only: eq_style, use_normalization
        use eq_vars, only: pres_FD, grp_n_r_eq, grp_min_r_eq, n_r_eq, rho_0, rho
        use VMEC, only: gam
        
        character(*), parameter :: rout_name = 'calc_rho'
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp) :: expon                                                       ! exponent = 1/gam
        real(dp), parameter :: tol = 1.0E-10_dp                                 ! tolerance for negative pressure
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! allocate rho
        allocate(rho(grp_n_r_eq))
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! set exp
                expon = 1.0_dp/gam
                
                ! loop over all normal points
                do kd = 1,grp_n_r_eq
                    if (pres_FD(kd,0).gt.0) then
                        rho(kd) = pres_FD(kd,0)**expon
                    else
                        rho(kd) = rho(kd-1)
                        if (pres_FD(kd,0).lt.-tol) &
                        call writo('WARNING: pressure was negative ('//&
                            &trim(r2strt(pres_FD(kd,0)))//') at point '&
                            &//trim(i2str(grp_min_r_eq-1+kd))//'/'//&
                            &trim(i2str(n_r_eq))//' so density is set to '//&
                            &trim(r2str(rho(kd))),persistent=.true.)
                    end if
                end do
            case (2)                                                            ! HELENA
                rho = 1.0_dp                                                    ! arbitrarily constant (normalized value)
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! normalize rho
        if (use_normalization) rho = rho/rho_0
    end function calc_rho
    
    ! calculate  ~PV_(k,m)^i  (pol.  flux)  or ~PV_(l,n)^i  (tor.  flux) at  all
    ! grp_n_r_eq values
    ! (see [ADD REF] for details)
    subroutine calc_PV
        use metric_vars, only: g_FD, h_FD, jac_FD
        use eq_vars, only: q_saf_FD, rot_t_FD, grp_n_r_eq
        use num_vars, only: use_pol_flux_X, use_normalization, mu_0
        use X_vars, only: PV0, PV1, PV2, DU_X_0, DU_X_1, mu0sigma, extra1, &
            &extra2, extra3, size_X, m_X, n_X, n_par
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:)                                   ! common factor |nabla psi|^2/(J^2*B^2)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp) :: mu_0_loc                                                    ! local version of mu_0
        
        ! submatrices
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        ! lower metric factors
        real(dp), allocatable :: g33(:,:)                                       ! h^alpha,psi
        ! upper metric factors
        real(dp), allocatable :: h22(:,:)                                       ! h^alpha,psi
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,grp_n_r_eq)); J = jac_FD(:,:,0,0,0)
        ! lower metric factors
        allocate(g33(n_par,grp_n_r_eq)); g33 = g_FD(:,:,3,3,0,0,0)
        ! upper metric factors
        allocate(h22(n_par,grp_n_r_eq)); h22 = h_FD(:,:,2,2,0,0,0)
        
        ! set up local mu_0
        if (use_normalization) then
            mu_0_loc = 1._dp
        else
            mu_0_loc = mu_0
        end if
        
        ! set up common factor for PVi
        allocate(com_fac(n_par,grp_n_r_eq))
        com_fac = h22/(g33*mu_0_loc)
        
        ! set up fac_n and fac_m
        allocate(fac_n(grp_n_r_eq),fac_m(grp_n_r_eq))
        if (use_pol_flux_X) then
            fac_n = q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = rot_t_FD(:,0)
        end if
        
        ! calculate PVi
        do m = 1,size_X
            do k = 1,m
                ! calculate PV0
                PV0(:,:,k,m) = com_fac * (DU_X_0(:,:,m) - extra1 - extra2 ) * &
                    &(conjg(DU_X_0(:,:,k)) - extra1 - extra2) - &
                    &mu_0/mu_0_loc * (mu0sigma/J * (extra1 + extra2) - extra3)
                
                ! add (nq-k)*(nq-m)/(J^2 |nabla psi|^2) to PV0
                do kd = 1,grp_n_r_eq
                    PV0(:,kd,k,m) = PV0(:,kd,k,m) + 1._dp/mu_0_loc * &
                        &(n_X(m)*fac_n(kd)-m_X(m)*fac_m(kd))*&
                        &(n_X(k)*fac_n(kd)-m_X(k)*fac_m(kd)) / &
                        &( J(:,kd)**2*h22(:,kd) )
                end do
                
                ! calculate PV2
                PV2(:,:,k,m) = com_fac * DU_X_1(:,:,m) * conjg(DU_X_1(:,:,k))
            end do
        end do
        
        ! fill the Hermitian conjugates
        do m = 1,size_X
            do k = m+1,size_X
                PV0(:,:,k,m) = conjg(PV0(:,:,m,k))
                PV2(:,:,k,m) = conjg(PV2(:,:,m,k))
            end do
        end do
        
        ! PV1 is not Hermitian
        do m = 1,size_X
            do k = 1,size_X
                ! calculate PV1
                PV1(:,:,k,m) = com_fac * DU_X_1(:,:,m) * &
                    &(conjg(DU_X_0(:,:,k)) - extra1 - extra2)
                
                !write(*,*) 'PV1'
                !call print_GP_2D('RE PV2(:,:,'//trim(i2str(k))//','//&
                    !&trim(i2str(m))//')','',realpart(transpose(PV1(:,:,k,m))))
                !call print_GP_2D('IM PV2(:,:,'//trim(i2str(k))//','//&
                    !&trim(i2str(m))//')','',imagpart(transpose(PV1(:,:,k,m))))
            end do
        end do
        
        !! output for test
        !write(*,*) 'non-integrated values'
        !write(*,*) 'RE max PV = ', maxval(realpart(PV0)), &
            !&maxval(realpart(PV1)), maxval(realpart(PV2))
        !write(*,*) 'IM max PV = ', maxval(imagpart(PV0)), &
            !&maxval(imagpart(PV1)), maxval(imagpart(PV2))
        !write(*,*) 'RE min PV = ', minval(realpart(PV0)), &
            !&minval(realpart(PV1)), minval(realpart(PV2))
        !write(*,*) 'IM min PV = ', minval(imagpart(PV0)), &
            !&minval(imagpart(PV1)), minval(imagpart(PV2))
        !allocate(var_names(n_par))
        !do kd = 1,n_par
            !var_names(kd) = 'i_par = '//trim(i2str(kd))
        !end do
        !tot_dim = [n_par,n_r_eq,size_X,size_X]
        !grp_dim = [n_par,grp_n_r_eq,size_X,size_X]
        !grp_offset = [0,grp_min_r_eq-1,0,0]
        !file_name = 'RE PV0'
        !call print_HDF5(var_names,file_name,realpart(PV0),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'IM PV0'
        !call print_HDF5(var_names,file_name,imagpart(PV0),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'RE PV1'
        !call print_HDF5(var_names,file_name,realpart(PV1),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'IM PV1'
        !call print_HDF5(var_names,file_name,imagpart(PV1),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'RE PV2'
        !call print_HDF5(var_names,file_name,realpart(PV2),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'IM PV2'
        !call print_HDF5(var_names,file_name,imagpart(PV2),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !deallocate(var_names)
    end subroutine calc_PV
    
    ! calculate  ~KV_(k,m)^i  (pol.  flux)  or ~KV_(l,n)^i  (tor.  flux) at  all
    ! grp_n_r_eq values
    ! (see [ADD REF] for details)
    subroutine calc_KV
        use metric_vars, only: g_FD, h_FD, jac_FD
        use eq_vars, only: grp_n_r_eq, rho
        use output_ops, only: print_HDF5
        use X_vars, only: KV0, KV1, KV2, U_X_0, U_X_1, size_X, n_par
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:)                                   ! common factor |nabla psi|^2/(J^2*B^2)
        
        ! submatrices
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        ! lower metric factors
        real(dp), allocatable :: g33(:,:)                                       ! h^alpha,psi
        ! upper metric factors
        real(dp), allocatable :: h22(:,:)                                       ! h^alpha,psi
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,grp_n_r_eq)); J = jac_FD(:,:,0,0,0)
        ! lower metric factors
        allocate(g33(n_par,grp_n_r_eq)); g33 = g_FD(:,:,3,3,0,0,0)
        ! upper metric factors
        allocate(h22(n_par,grp_n_r_eq)); h22 = h_FD(:,:,2,2,0,0,0)
        
        ! set up common factor
        allocate(com_fac(n_par,grp_n_r_eq))
        com_fac = J**2*h22/g33
        
        ! for  Hermitian KV0  and  KV2,  only half  of  the  terms  have  to  be
        ! calculated
        do m = 1,size_X
            do k = 1,m
                ! calculate KV0
                KV0(:,:,k,m) = com_fac * U_X_0(:,:,m) * conjg(U_X_0(:,:,k)) + &
                    &1._dp/h22
                
                ! calculate KV2
                KV2(:,:,k,m) = com_fac * U_X_1(:,:,m) * conjg(U_X_1(:,:,k))
            end do
        end do
        
        ! fill the Hermitian conjugates
        do m = 1,size_X
            do k = m+1,size_X
                KV0(:,:,k,m) = conjg(KV0(:,:,m,k))
                KV2(:,:,k,m) = conjg(KV2(:,:,m,k))
            end do
        end do
        
        ! KV1 is not Hermitian
        do m = 1,size_X
            do k = 1,size_X
                ! calculate KV1
                KV1(:,:,k,m) = com_fac * U_X_1(:,:,m) * conjg(U_X_0(:,:,k))
            end do
        end do
        
        ! multiply by rho
        do kd = 1,grp_n_r_eq
            KV0(:,kd,:,:) = KV0(:,kd,:,:)*rho(kd)
            KV1(:,kd,:,:) = KV1(:,kd,:,:)*rho(kd)
            KV2(:,kd,:,:) = KV2(:,kd,:,:)*rho(kd)
        end do
        
        !! output for test
        !write(*,*) 'RE max PV = ', maxval(realpart(PV0)), &
            !&maxval(realpart(PV1)), maxval(realpart(PV2))
        !write(*,*) 'IM max PV = ', maxval(imagpart(PV0)), &
            !&maxval(imagpart(PV1)), maxval(imagpart(PV2))
        !write(*,*) 'RE min PV = ', minval(realpart(PV0)), &
            !&minval(realpart(PV1)), minval(realpart(PV2))
        !write(*,*) 'IM min PV = ', minval(imagpart(PV0)), &
            !&minval(imagpart(PV1)), minval(imagpart(PV2))
        !allocate(var_names(n_par))
        !do kd = 1,n_par
            !var_names(kd) = 'i_par = '//trim(i2str(kd))
        !end do
        !tot_dim = [n_par,n_r_eq,size_X,size_X]
        !grp_dim = [n_par,grp_n_r_eq,size_X,size_X]
        !grp_offset = [0,grp_min_r_eq-1,0,0]
        !file_name = 'RE KV0'
        !call print_HDF5(var_names,file_name,realpart(KV0),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'IM KV0'
        !call print_HDF5(var_names,file_name,imagpart(KV0),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'RE KV1'
        !call print_HDF5(var_names,file_name,realpart(KV1),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'IM KV1'
        !call print_HDF5(var_names,file_name,imagpart(KV1),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'RE KV2'
        !call print_HDF5(var_names,file_name,realpart(KV2),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !file_name = 'IM KV2'
        !call print_HDF5(var_names,file_name,imagpart(KV2),tot_dim,grp_dim,&
            !&grp_offset,col_id=1,col=1)
        !deallocate(var_names)
    end subroutine calc_KV
    
    ! calculate U_m^0, U_m^1 or U_n^0, U_n^1  at grp_n_r_eq values of the normal
    ! coordinate, n_par values  of the parallel coordinate and  size_X values of
    ! the poloidal mode number,  or of the toroidal mode number,  as well as the
    ! poloidal derivatives
    ! Note: The  poloidal  derivatives  have  the  factor  i/n  or i/m  included
    ! already, as opposed to [ADD REF]
    integer function calc_U() result(ierr)
        use num_vars, only: use_pol_flux_X, mu_0, use_normalization, eq_style
        use metric_vars, only: g_FD, h_FD, jac_FD
        use eq_vars, only: rot_t_FD, q_saf_FD, p => pres_FD, &
            &grp_n_r_eq
        use X_vars, only: n_X, m_X, U_X_0, U_X_1, DU_X_0, DU_X_1, size_X, &
            &n_par, ang_par_F
        
        character(*), parameter :: rout_name = 'calc_U'
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        real(dp) :: mu_0_loc                                                    ! local version of mu_0
        real(dp) :: n_frac                                                      ! nq-m (pol. flux) or n-iotam (tor. flux)
        real(dp), allocatable :: djq(:)                                         ! either q' (pol. flux) or -iota' (tor. flux)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp), allocatable :: mn(:)                                          ! either n*A_0 (pol. flux) or m (tor.flux)
        complex(dp), allocatable :: U_corr(:,:,:)                               ! correction to U for a certain (n,m)
        complex(dp), allocatable :: D3U_corr(:,:,:)                             ! D_theta U_corr
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! helper variables for the correction to U
        real(dp), allocatable :: g_frac(:,:)                                    ! g_alpha,theta / g_theta,theta
        real(dp), allocatable :: T_theta(:,:), D1T_theta(:,:), D3T_theta(:,:)   ! Theta^theta and derivatives
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        real(dp), allocatable :: D3J(:,:)                                       ! D_theta jac
        ! lower metric factors
        real(dp), allocatable :: g13(:,:)                                       ! g_alpha,theta
        real(dp), allocatable :: D3g13(:,:)                                     ! D_theta g_alpha,theta
        real(dp), allocatable :: g23(:,:)                                       ! g_psi,theta
        real(dp), allocatable :: D3g23(:,:)                                     ! D_theta g_psi,theta
        real(dp), allocatable :: g33(:,:)                                       ! g_theta,theta
        real(dp), allocatable :: D3g33(:,:)                                     ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), allocatable :: h12(:,:)                                       ! h^alpha,psi
        real(dp), allocatable :: D3h12(:,:)                                     ! D_theta h^alpha,psi
        real(dp), allocatable :: h22(:,:)                                       ! h^psi,psi
        real(dp), allocatable :: D1h22(:,:)                                     ! D_alpha h^psi,psi
        real(dp), allocatable :: D3h22(:,:)                                     ! D_theta h^psi,psi
        real(dp), allocatable :: h23(:,:)                                       ! h^psi,theta
        real(dp), allocatable :: D1h23(:,:)                                     ! D_alpha h^psi,theta
        real(dp), allocatable :: D3h23(:,:)                                     ! D_theta h^psi,theta
        
        ! initialize ierr
        ierr = 0
        
        ! set up local mu_0
        if (use_normalization) then
            mu_0_loc = 1._dp
        else
            mu_0_loc = mu_0
        end if
        
        ! set up djq, fac_n, fac_m and mn
        allocate(djq(grp_n_r_eq),fac_n(grp_n_r_eq),fac_m(grp_n_r_eq),mn(size_X))
        if (use_pol_flux_X) then
            djq = q_saf_FD(:,1)
            fac_n = q_saf_FD(:,0)
            fac_m = 1.0_dp
            mn = n_X
        else
            djq = -rot_t_FD(:,1)
            fac_n = 1.0_dp
            fac_m = rot_t_FD(:,0)
            mn = m_X
        end if
        
        ! allocate helper variables
        allocate(U_corr(n_par,grp_n_r_eq,size(m_X)))
        allocate(D3U_corr(n_par,grp_n_r_eq,size(m_X)))
        allocate(g_frac(n_par,grp_n_r_eq))
        allocate(T_theta(n_par,grp_n_r_eq))
        allocate(D1T_theta(n_par,grp_n_r_eq),D3T_theta(n_par,grp_n_r_eq))
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,grp_n_r_eq)); J = jac_FD(:,:,0,0,0)
        allocate(D3J(n_par,grp_n_r_eq)); D3J = jac_FD(:,:,0,0,1)
        ! lower metric factors
        allocate(g13(n_par,grp_n_r_eq)); g13 = g_FD(:,:,1,3,0,0,0)
        allocate(D3g13(n_par,grp_n_r_eq)); D3g13 = g_FD(:,:,1,3,0,0,1)
        allocate(g23(n_par,grp_n_r_eq)); g23 = g_FD(:,:,2,3,0,0,0)
        allocate(D3g23(n_par,grp_n_r_eq)); D3g23 = g_FD(:,:,2,3,0,0,1)
        allocate(g33(n_par,grp_n_r_eq)); g33 = g_FD(:,:,3,3,0,0,0)
        allocate(D3g33(n_par,grp_n_r_eq)); D3g33 = g_FD(:,:,3,3,0,0,1)
        ! upper metric factors
        allocate(h12(n_par,grp_n_r_eq)); h12 = h_FD(:,:,1,2,0,0,0)
        allocate(D3h12(n_par,grp_n_r_eq)); D3h12 = h_FD(:,:,1,2,0,0,1)
        allocate(h22(n_par,grp_n_r_eq)); h22 = h_FD(:,:,2,2,0,0,0)
        allocate(D1h22(n_par,grp_n_r_eq)); D1h22 = h_FD(:,:,2,2,1,0,0)
        allocate(D3h22(n_par,grp_n_r_eq)); D3h22 = h_FD(:,:,2,2,0,0,1)
        allocate(h23(n_par,grp_n_r_eq)); h23 = h_FD(:,:,2,3,0,0,0)
        allocate(D1h23(n_par,grp_n_r_eq)); D1h23 = h_FD(:,:,2,3,1,0,0)
        allocate(D3h23(n_par,grp_n_r_eq)); D3h23 = h_FD(:,:,2,3,0,0,1)
        
        ! calculate U_X_0, U_X_1, DU_X_0 and DU_X_1
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
        
        !do jd = 1,size(m_X)
            !write(*,*) 'U0'
            !write(*,*) 'RE U0'
            !call print_GP_2D('RE U0(:,:,'//trim(i2str(jd))//')','',&
                !&transpose(realpart(U_X_0(:,:,jd))))
            !write(*,*) 'RE U_corr'
            !call print_GP_2D('RE U_corr(:,:,'//trim(i2str(jd))//')','',&
                !&transpose(realpart(U_corr(:,:,jd))))
            !write(*,*) 'FRACTION RE U_corr'
            !call print_GP_2D('RE U_corr(:,:,'//trim(i2str(jd))//')','',&
                !&transpose(realpart(U_corr(:,:,jd)/U_X_0(:,:,jd))))
            !write(*,*) 'IM U0'
            !call print_GP_2D('IM U0(:,:,'//trim(i2str(jd))//')','',&
                !&transpose(imagpart(U_X_0(:,:,jd))))
            !write(*,*) 'IM U_corr'
            !call print_GP_2D('IM RE U_corr(:,:,'//trim(i2str(jd))//')','',&
                !&transpose(imagpart(U_corr(:,:,jd))))
            !write(*,*) 'FRACTION IM U_corr'
            !call print_GP_2D('IM U_corr(:,:,'//trim(i2str(jd))//')','',&
                !&transpose(imagpart(U_corr(:,:,jd)/U_X_0(:,:,jd))))
        !end do
        
        ! deallocate
        deallocate(U_corr,D3U_corr)
        deallocate(J,D3J)
        deallocate(g13,D3g13,g23,D3g23,g33,D3g33)
        deallocate(h12,D3h12,h22,D1h22,D3h22,h23,D1h23,D3h23)
        deallocate(g_frac)
        deallocate(T_theta,D1T_theta,D3T_theta)
    contains
        ! VMEC version
        subroutine calc_U_VMEC
            ! local variables
            ! extra helper variables for the correction to U
            real(dp), allocatable :: D3g_frac(:,:)                              ! D_theta g_frac
            real(dp), allocatable :: D13T_theta(:,:), D33T_theta(:,:)           ! Theta^theta derivatives
            ! extra upper metric factors
            real(dp), allocatable :: D13h22(:,:)                                ! D^2_alpha,theta h^psi,psi
            real(dp), allocatable :: D33h22(:,:)                                ! D^2_theta,theta h^psi,psi
            real(dp), allocatable :: D13h23(:,:)                                ! D^2_alpha,theta h^psi,theta
            real(dp), allocatable :: D33h23(:,:)                                ! D^2_theta,theta h^psi,theta
            
            ! allocate extra helper variables
            allocate(D3g_frac(n_par,grp_n_r_eq))
            allocate(D13T_theta(n_par,grp_n_r_eq),D33T_theta(n_par,grp_n_r_eq))
            
            ! set up extra submatrices
            ! upper metric factors
            allocate(D13h22(n_par,grp_n_r_eq)); D13h22 = h_FD(:,:,2,2,1,0,1)
            allocate(D33h22(n_par,grp_n_r_eq)); D33h22 = h_FD(:,:,2,2,0,0,2)
            allocate(D13h23(n_par,grp_n_r_eq)); D13h23 = h_FD(:,:,2,3,1,0,1)
            allocate(D33h23(n_par,grp_n_r_eq)); D33h23 = h_FD(:,:,2,3,0,0,2)
            
            ! loop over the M elements of U_X and DU
            do jd = 1,size(m_X)
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
                do kd = 1,grp_n_r_eq
                    ! set up correction to U
                    n_frac = n_X(jd)*fac_n(kd)-m_X(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*(g_frac(:,kd)*&
                        &(D3T_theta(:,kd)+iu*n_frac*T_theta(:,kd))+&
                        &D1T_theta(:,kd))
                    ! set up D_theta U correction
                    D3U_corr(:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*&
                        &(D3g_frac(:,kd)*&
                        &(D3T_theta(:,kd)+iu*n_frac*T_theta(:,kd))+&
                        &g_frac(:,kd)*&
                        &(D33T_theta(:,kd)+iu*n_frac*D3T_theta(:,kd))+&
                        &D13T_theta(:,kd))
                    ! calculate U_X_0 and DU_X_0
                    U_X_0(:,kd,jd) = &
                        &-(h12(:,kd)/h22(:,kd) + djq(kd)*ang_par_F(:,kd))&
                        &+ iu/(mn(jd)*g33(:,kd)) * (g13(:,kd)*djq(kd) + &
                        &J(:,kd)**2*p(kd,1)*mu_0_loc + iu*n_frac * &
                        &( g13(:,kd)*djq(kd)*ang_par_F(:,kd) - g23(:,kd) )) &
                        &+ U_corr(:,kd,jd)
                    DU_X_0(:,kd,jd) = -(D3h12(:,kd)/h22(:,kd) - &
                        &D3h22(:,kd)*h12(:,kd)/h22(:,kd)**2 + djq(kd)) - &
                        &iu*D3g33(:,kd)/(mn(jd)*g33(:,kd)**2) * &
                        &(g13(:,kd)*djq(kd)+ J(:,kd)**2*p(kd,1)*mu_0_loc + &
                        &iu*n_frac * &
                        &( g13(:,kd)*djq(kd)*ang_par_F(:,kd) - g23(:,kd) )) + &
                        &iu/(mn(jd)*g33(:,kd)) * (D3g13(:,kd)*djq(kd) + &
                        &2*D3J(:,kd)*J(:,kd)*p(kd,1)*mu_0_loc + iu*n_frac * &
                        &( D3g13(:,kd)*djq(kd)*ang_par_F(:,kd) + &
                        &g13(:,kd)*djq(kd) - D3g23(:,kd) )) + &
                        &iu*n_frac*U_X_0(:,kd,jd) &
                        &+ D3U_corr(:,kd,jd)
                    ! calculate U_X_1 and DU_X_1
                    U_X_1(:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,kd)/g33(:,kd))
                    DU_X_1(:,kd,jd) = iu/mn(jd) * &
                        &(n_frac/mn(jd) * (D3g13(:,kd)/g33(:,kd) - &
                        &D3g33(:,kd)*g13(:,kd)/g33(:,kd)**2) + &
                        &iu*n_frac*U_X_1(:,kd,jd) )
                end do
            end do
            
            ! deallocate
            deallocate(D3g_frac)
            deallocate(D13T_theta,D33T_theta)
        end subroutine calc_U_VMEC
        
        ! HELENA version
        integer function calc_U_HEL() result(ierr)
            use eq_vars, only: theta_E
            use utilities, only: calc_deriv
            
            character(*), parameter :: rout_name = 'calc_U_HEL'
            
            ! local variablas
            real(dp), allocatable :: D3_var(:,:)                                ! derivative of variable
            
            ! allocate extra helper variables
            allocate(D3_var(n_par,2))
            
            ! initialize ierr
            ierr = 0
            
            write(*,*) 'INVESTIGATE HOW IMPROVING THE DERIVATIVES CAN &
                &HELP YOU GET RID OF UNSTABLE SIDE OF SPECTRUM !!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) 'DOING THIS IN THIS ROUTINE ALREADY HELPED A LOT!!!!!'
                
            ! loop over the M elements of U_X and DU
            do jd = 1,size(m_X)
                ! set up helper variables
                g_frac = g13/g33
                T_theta = h23/h22
                ! loop over all normal points
                do kd = 1,grp_n_r_eq
                    ierr = calc_deriv(T_theta(:,kd),D3T_theta(:,kd),&
                        &theta_E(:,kd),1,2)                                     ! higher precision because other derivative will be taken later
                end do
                
                ! loop over all normal points
                do kd = 1,grp_n_r_eq
                    ! set up correction to U
                    n_frac = n_X(jd)*fac_n(kd)-m_X(jd)*fac_m(kd)
                    ! set up U correction
                    U_corr(:,kd,jd) = iu/mn(jd)*n_frac/mn(jd)*(g_frac(:,kd)*&
                        &(D3T_theta(:,kd)+iu*n_frac*T_theta(:,kd)))
                    ! calculate U_X_0 and DU_X_0
                    U_X_0(:,kd,jd) = &
                        &-(h12(:,kd)/h22(:,kd) + djq(kd)*ang_par_F(:,kd))&
                        &+ iu/(mn(jd)*g33(:,kd)) * (g13(:,kd)*djq(kd) + &
                        &J(:,kd)**2*p(kd,1)*mu_0_loc + iu*n_frac * &
                        &( g13(:,kd)*djq(kd)*ang_par_F(:,kd) - g23(:,kd) )) &
                        &+ U_corr(:,kd,jd)
                    ierr = calc_deriv(realpart(U_X_0(:,kd,jd)),D3_var(:,1),&
                        &theta_E(:,kd),1,1)
                    ierr = calc_deriv(imagpart(U_X_0(:,kd,jd)),D3_var(:,2),&
                        &theta_E(:,kd),1,1)
                    CHCKERR('')
                    DU_X_0(:,kd,jd) = D3_var(:,1) + iu*D3_var(:,2) + &
                        &iu*n_frac*U_X_0(:,kd,jd)
                    ! calculate U_X_1 and DU_X_1
                    U_X_1(:,kd,jd) = iu/mn(jd) * &
                        &(1 + n_frac/mn(jd) * g13(:,kd)/g33(:,kd))
                    ierr = calc_deriv(realpart(U_X_1(:,kd,jd)),D3_var(:,1),&
                        &theta_E(:,kd),1,1)
                    ierr = calc_deriv(imagpart(U_X_1(:,kd,jd)),D3_var(:,2),&
                        &theta_E(:,kd),1,1)
                    CHCKERR('')
                    DU_X_1(:,kd,jd) = D3_var(:,1) + iu*D3_var(:,2) + &
                        &iu*n_frac*U_X_1(:,kd,jd)
                end do
            end do
            
            ! deallocate
            deallocate(D3_var)
        end function calc_U_HEL
    end function calc_U
    
    ! Calculate extra1, extra2 and extra3 in non-normalized coordinates:
    !   extra1 = S*J
    !   extra2 = mu0sigma*J*B^2/h^psi,psi
    !   extra3 = 2*p'*kn
    ! with
    !   shear S = - d Theta^alpha/d_theta 1 / J
    !   parallel current mu0sigma = B . nabla x B / B^2
    !   normal curvature kn = nabla psi / h^psi,psi * nabla (mu0 p + B^2/2)
    integer function calc_extra() result(ierr)
        use metric_vars, only: g_FD, h_FD, jac_FD
        use eq_vars, only: p => pres_FD, grp_n_r_eq
        use num_vars, only: mu_0
        use X_vars, only: mu0sigma, extra1, extra2, extra3, n_par
        
        character(*), parameter :: rout_name = 'calc_extra'
        
        ! local variables
        integer :: kd                                                           ! counter
        
        ! submatrices
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        real(dp), allocatable :: D1J(:,:)                                       ! D_alpha jac
        real(dp), allocatable :: D2J(:,:)                                       ! D_psi jac
        real(dp), allocatable :: D3J(:,:)                                       ! D_theta jac
        ! lower metric factors
        real(dp), allocatable :: g13(:,:)                                       ! g_alpha,theta
        real(dp), allocatable :: D2g13(:,:)                                     ! D_psi g_alpha,theta
        real(dp), allocatable :: D3g13(:,:)                                     ! D_theta g_alpha,theta
        real(dp), allocatable :: g23(:,:)                                       ! g_psi,theta
        real(dp), allocatable :: D1g23(:,:)                                     ! D_alpha g_psi,theta
        real(dp), allocatable :: D3g23(:,:)                                     ! D_theta g_psi,theta
        real(dp), allocatable :: g33(:,:)                                       ! g_theta,theta
        real(dp), allocatable :: D1g33(:,:)                                     ! D_alpha g_theta,theta
        real(dp), allocatable :: D2g33(:,:)                                     ! D_psi g_theta,theta
        real(dp), allocatable :: D3g33(:,:)                                     ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), allocatable :: h12(:,:)                                       ! h^alpha,psi
        real(dp), allocatable :: D3h12(:,:)                                     ! D_theta h^alpha,psi
        real(dp), allocatable :: h22(:,:)                                       ! h^psi,psi
        real(dp), allocatable :: D3h22(:,:)                                     ! D_theta h^psi,psi
        real(dp), allocatable :: h23(:,:)                                       ! h^psi,theta
        
        ! initialize ierr
        ierr = 0
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,grp_n_r_eq)); J = jac_FD(:,:,0,0,0)
        allocate(D1J(n_par,grp_n_r_eq)); D1J = jac_FD(:,:,1,0,0)
        allocate(D2J(n_par,grp_n_r_eq)); D2J = jac_FD(:,:,0,1,0)
        allocate(D3J(n_par,grp_n_r_eq)); D3J = jac_FD(:,:,0,0,1)
        ! lower metric factors
        allocate(g13(n_par,grp_n_r_eq)); g13 = g_FD(:,:,1,3,0,0,0)
        allocate(D2g13(n_par,grp_n_r_eq)); D2g13 = g_FD(:,:,1,3,0,1,0)
        allocate(D3g13(n_par,grp_n_r_eq)); D3g13 = g_FD(:,:,1,3,0,0,1)
        allocate(g23(n_par,grp_n_r_eq)); g23 = g_FD(:,:,2,3,0,0,0)
        allocate(D1g23(n_par,grp_n_r_eq)); D1g23 = g_FD(:,:,2,3,1,0,0)
        allocate(D3g23(n_par,grp_n_r_eq)); D3g23 = g_FD(:,:,2,3,0,0,1)
        allocate(g33(n_par,grp_n_r_eq)); g33 = g_FD(:,:,3,3,0,0,0)
        allocate(D1g33(n_par,grp_n_r_eq)); D1g33 = g_FD(:,:,3,3,1,0,0)
        allocate(D2g33(n_par,grp_n_r_eq)); D2g33 = g_FD(:,:,3,3,0,1,0)
        allocate(D3g33(n_par,grp_n_r_eq)); D3g33 = g_FD(:,:,3,3,0,0,1)
        ! upper metric factors
        allocate(h12(n_par,grp_n_r_eq)); h12 = h_FD(:,:,1,2,0,0,0)
        allocate(D3h12(n_par,grp_n_r_eq)); D3h12 = h_FD(:,:,1,2,0,0,1)
        allocate(h22(n_par,grp_n_r_eq)); h22 = h_FD(:,:,2,2,0,0,0)
        allocate(D3h22(n_par,grp_n_r_eq)); D3h22 = h_FD(:,:,2,2,0,0,1)
        allocate(h23(n_par,grp_n_r_eq)); h23 = h_FD(:,:,2,3,0,0,0)
        
        ! allocate  mu0sigma  and extra  terms.  Later they  are deallocated  in
        ! calc_PV
        allocate(mu0sigma(n_par,grp_n_r_eq))
        allocate(extra1(n_par,grp_n_r_eq))
        allocate(extra2(n_par,grp_n_r_eq))
        allocate(extra3(n_par,grp_n_r_eq))
        
        ! calculate mu0sigma
        mu0sigma = 1._dp/(J*g33) * (g13*(D2g33-D3g23) + g23*(D3g13-D1g33) &
            &+ g33*(D1g23-D2g13))
        
        ! calculate extra1
        extra1 = -D3h12/h22 + D3h22*h12/h22**2
        
        ! calculate extra2
        extra2 = g33/h22 * mu0sigma / J
        
        ! calculate extra3
        do kd = 1,grp_n_r_eq
            extra3(:,kd) = p(kd,1) * ( 2*J(:,kd)**2*mu_0*p(kd,1)/g33(:,kd) + &
                &1._dp/h22(:,kd) * ( &
                &h12(:,kd) * ( D1g33(:,kd)/g33(:,kd) - 2*D1J(:,kd)/J(:,kd) ) + &
                &h22(:,kd) * ( D2g33(:,kd)/g33(:,kd) - 2*D2J(:,kd)/J(:,kd) ) + &
                &h23(:,kd) * ( D3g33(:,kd)/g33(:,kd) - 2*D3J(:,kd)/J(:,kd) ) ) )
        end do
        
        ! calculate rho from input
        ierr = calc_rho()
        CHCKERR('')
    end function calc_extra
    
    ! calculates  magnetic  integral  <V e^[i(k-m)ang_par_F]>,  defined  as  the
    ! matrix  
    !   <V e^[i(k-m)ang_par_F]> = [ oint J V(k,m) e^i(k-m)ang_par_F dang_par_F ]
    ! or
    !   <V e^[i(n-l)ang_par_F]> = [ oint J V(l,n) e^i(n-l)ang_par_F dang_par_F ]
    ! Note: V is changed on exit
    subroutine calc_V_int(V,V_int)
        use eq_vars, only: grp_n_r_eq
        use metric_vars, only: jac_FD
        use X_vars, only: size_X, exp_ang_par_F, n_par, ang_par_F
        
        ! input / output
        complex(dp), intent(inout) :: V(n_par,grp_n_r_eq,size_X,size_X)         ! input V(n_par,n_r,size_X,size_X)
        complex(dp), intent(inout) :: V_int(size_X,size_X,grp_n_r_eq)           ! output <V e^i(k-m)ang_par_F> at normal equilibrium grid points
        
        ! local variables
        integer :: k, m, jd, kd                                                 ! counters
        
        do m = 1,size_X
            do k = 1,size_X
                ! multiply V by Jacobian and exponential
                V(:,:,k,m) = exp_ang_par_F(:,:,k,m) * jac_FD(:,:,0,0,0) * &
                    &V(:,:,k,m)
            end do
        end do
        
        ! integrate term  over ang_par_F for  all equilibrium grid  points using
        ! the recursive formula int_1^n f(x) dx
        !   = int_1^(n-1) f(x) dx + (f(n)+f(n-1))*(x(n)-x(n-1))/2
        V_int = 0.0_dp
        ! loop over all normal points on this process
        do kd = 1,grp_n_r_eq
            ! parallel integration loop
            do jd = 2, n_par
                V_int(:,:,kd) = V_int(:,:,kd) + (V(jd,kd,:,:)+V(jd-1,kd,:,:))/2&
                    & * (ang_par_F(jd,kd)-ang_par_F(jd-1,kd))
            end do
        end do
    end subroutine calc_V_int
    
    ! Plots an Eigenfunction
    ! [MPI] Collective call
    integer function plot_X_vec(X_vec,X_val,X_id,job_id,par_range_F) &
        &result(ierr)
        use num_vars, only: output_style
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        real(dp), intent(in), optional :: par_range_F(2)                        ! parallel range [2pi]
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! local variables
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        
        ! initialize ierr
        ierr = 0
        
        ! if omega  is predominantly  imaginary, the Eigenfunction  explodes, so
        ! limit  n_t(2) to  1. If  it is  predominantly real,  the Eigenfunction
        ! oscilates, so choose n_t(2) = 8 for 2 whole periods
        omega = sqrt(X_val)
        if (abs(realpart(omega)/imagpart(omega)).lt.1) then                     ! exploding, unstable
            if (imagpart(omega).gt.0) omega = - omega                           ! exploding solution, not the decaying one
            n_t(1) = 1                                                          ! 1 point per quarter period
            n_t(2) = 1                                                          ! 1 quarter period
        else                                                                    ! oscillating, stable
            n_t(1) = 2                                                          ! 2 points per quarter period
            n_t(2) = 4                                                          ! 4 quarter periods
        end if
        
        ! delegate to child routines
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                ierr = plot_X_vec_GP(X_vec,omega,X_id,job_id,n_t,par_range_F)
                CHCKERR('')
            case(2)                                                             ! HDF5 output
                ierr = plot_X_vec_HDF5(X_vec,omega,X_id,job_id,n_t,.false.)
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
    ! (pi/2)/omega and  an animated plot  is produced as  well in a  .gif output
    ! file
    integer function plot_X_vec_GP(X_vec,omega,X_id,job_id,n_t,par_range_F) &
        &result(ierr)
        use num_vars, only: grp_rank, max_n_plots, grp_n_procs, &
            &use_pol_flux_X, use_pol_flux_eq
        use output_ops, only: draw_GP, draw_GP_animated, merge_GP
        use X_vars, only: n_r_X, size_X, n_X, m_X, grp_r_X, min_r_X, max_r_X
        use MPI_ops, only: get_ghost_X_vec, wait_MPI
        use eq_vars, only: q_saf_FD, rot_t_FD, grp_n_r_eq, max_flux_X_F, &
            &max_flux_eq_F, flux_p_FD, flux_t_FD
        use utilities, only: con2dis, interp_fun
        
        character(*), parameter :: rout_name = 'plot_X_vec_GP'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: omega                                        ! sqrt of Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        integer, intent(in) :: n_t(2)                                           ! nr. of time steps in quarter period, nr. of quarter periods
        real(dp), intent(in), optional :: par_range_F(2)                        ! parallel range [2pi]
        
        ! local variables
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
        
        ! set up fac_n and fac_m
        allocate(fac_n(grp_n_r_eq),fac_m(grp_n_r_eq))
        if (use_pol_flux_X) then
            fac_n = q_saf_FD(:,0)
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = rot_t_FD(:,0)
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
            allocate(X_vec_extended(size(X_vec,1),size(X_vec,2)+1))
        else
            allocate(X_vec_extended(size(X_vec,1),size(X_vec,2)))
        end if
        allocate(X_plot(grp_n_norm,n_par_F,size_X,product(n_t)))
        allocate(Y_plot(grp_n_norm,n_par_F,size_X,product(n_t)))
        allocate(Z_plot(grp_n_norm,n_par_F,size_X,product(n_t)))
        allocate(z_magn_plot(grp_n_norm,n_par_F,product(n_t)))
        allocate(X_vec_interp(size_X))
        allocate(first_X_vec_interp(size_X))
        
        ! set up extended X_vec
        X_vec_extended(:,1:size(X_vec,2)) = X_vec
        ierr = get_ghost_X_vec(X_vec_extended,1)
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
        if (use_pol_flux_X) then
            flux_X => flux_p_FD(:,0)
        else
            flux_X => flux_t_FD(:,0)
        end if
        if (use_pol_flux_eq) then
            flux_eq => flux_p_FD(:,0)
        else
            flux_eq => flux_t_FD(:,0)
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
                    call con2dis(r_plot(kd),kd_loc_i,grp_r_X)                   ! perturbation values tabulated at grp_r_X for this group
                    X_vec_interp = X_vec_extended(:,floor(kd_loc_i)) + &
                        &(kd_loc_i-floor(kd_loc_i)) * &
                        &(X_vec_extended(:,ceiling(kd_loc_i))-&
                        &X_vec_extended(:,floor(kd_loc_i)))
                    
                    ! set up  interp. fac_n and fac_m  (tabulated in equilibrium
                    ! grid)
                    ierr = interp_fun(kd_loc_eq,flux_eq/max_flux_eq_F,&
                        &r_plot(kd),flux_X/max_flux_X_F)
                    CHCKERR('')
                    call con2dis(kd_loc_eq,kd_loc_i,flux_eq/max_flux_eq_F)      ! equilibrium values tabulated at flux_eq for this group, normalized with max_flux_eq_F
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
                do ld = 1,size_X
                    do jd = 1,n_par_F
                        Z_plot(kd,jd,ld,id) = &
                            &realpart(X_vec_interp(ld) * exp(iu * &
                            &(omega/abs(omega)*0.5*pi*(id-1.0_dp)/&
                            &max(n_t(1)-1,1) + n_X(ld)*fac_n_interp-&
                            &          m_X(ld)*fac_m_interp)&
                            &    * ang_par_F_loc(jd)))
                    end do
                    z_magn_plot(kd,:,id) = z_magn_plot(kd,:,id) + &
                        &Z_plot(kd,:,ld,id)
                end do
            end do
        end do
        
        ! 3. for  every mode, write plot  output data in data files  and plot by
        ! group master
        if (size_X.gt.1) then
            ! plot the individual modes in time
            do ld = 1,min(size_X,max_n_plots)
                ! set plot_title and file name and write data
                plot_title = 'job '//trim(i2str(job_id))//' - EV '//&
                    &trim(i2str(X_id))//' - mode '//trim(i2str(ld))//' of '//&
                    &trim(i2str(size_X))
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
            if (size_X.eq.max_n_plots+1) then
                call writo('Skipping plot for mode '//&
                    &trim(i2str(max_n_plots+1)))
            else if (size_X.gt.max_n_plots+1) then
                call writo('Skipping plots for modes '//&
                    &trim(i2str(max_n_plots+1))//'..'//trim(i2str(size_X)))
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
            n_norm = min(n_r_X,max_n_norm)
            inc_norm = (max_r_X-min_r_X)/(n_norm-1)
            
            ! calculate the  starting and ending index,  first determining which
            ! is the range where normal interpolation can take place
            ! from:
            !   gpr_r_X(1) .le. (min_r_X+(k-1)*inc_norm) .le. grp_r_X(grp_n_r_X)
            ! with k an integer
            start_i_norm = ceiling((grp_r_X(1)-min_r_X)/inc_norm + 1)
            end_i_norm = floor((grp_r_X(size(grp_r_X))-min_r_X)/inc_norm + 1)   ! size(grp_r_X) can be grp_n_r_X or grp_n_r_X + 1 (ghost region in X grid)
            
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
            ierr = get_ghost_X_vec(exch_vec,1)
            CHCKERR('')
            
            ! copy value to new interpolated X_vec
            vec = exch_vec(:,size(exch_vec,2))
        end subroutine exch_ghost_values
    end function plot_X_vec_GP
    
    ! HDF5 version of plot_X_vec
    ! Plots an  eigenfunction in  3D real  space by determining  a grid  in VMEC
    ! coordinates that covers  the entire device (But they are  tabulated in the
    ! perturbation normal  grid) and then  calculating the real  perturbation on
    ! that grid  by inverting the Fourier  transform that defined the  modes for
    ! which are  solved.
    ! Alternatively, by using the optional argument "follow_B", the modes can be
    ! plot along  the magnetic  field lines  which have been  used to  solve the
    ! system of equations. (TO BE IMPLEMENTED)
    integer function plot_X_vec_HDF5(X_vec,omega,X_id,job_id,n_t,follow_B) &
        &result(ierr)
        use X_vars, only: grp_min_r_X, grp_n_r_X, n_r_X, grp_r_X, grp_n_r_X
        use num_vars, only: n_theta_plot, n_zeta_plot
        use utilities, only: interp_fun
        use grid_ops, only: calc_XYZ_grid, coord_F2E, coord_E2F
        use sol_ops, only: calc_real_X
        use output_ops, only: print_HDF5
        
        character(*), parameter :: rout_name = 'plot_X_vec_HDF5'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: omega                                        ! sqrt of Eigenvalue
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue (for output name)
        integer, intent(in) :: job_id                                           ! nr. of alpha job
        integer, intent(in) :: n_t(2)                                           ! nr. of time steps in quarter period, nr. of quarter periods
        logical, intent(in), optional :: follow_B                               ! .true. if to be plot only along B
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp), allocatable :: theta_E(:,:,:), zeta_E(:,:,:)                  ! theta_E and zeta_E for flux surface plot
        real(dp), allocatable :: r_E(:)                                         ! r for plots in pert. coords.
        real(dp), allocatable :: theta_F(:,:,:), zeta_F(:,:,:)                  ! theta_F and zeta_F for flux surface plot
        real(dp), allocatable :: r_F(:)                                         ! r for plots in eq. coords.
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
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Started the plot of the Eigenfunction')
        call lvl_ud(1)
        
        ! set up local follow_B
        follow_B_loc = .false.
        if (present(follow_B)) then
            if (follow_B) follow_B_loc = follow_B
        end if
        
        ! set up total dimensions, group dimensions and group offset
        tot_dim = [n_theta_plot,n_zeta_plot,n_r_X,product(n_t)]
        grp_dim = [n_theta_plot,n_zeta_plot,grp_n_r_X,product(n_t)]
        grp_offset = [0,0,grp_min_r_X-1,product(n_t)]
        
        ! set up var_name, file_name and description
        var_name = 'Solution vector X_vec'
        file_name = 'X_vec_'//trim(i2str(job_id))//'_'//trim(i2str(X_id))
        description = 'Job '//trim(i2str(job_id))//' - Solution vector X_vec &
            &for Eigenvalue '//trim(i2str(X_id))//' with omega = '&
            &//trim(r2str(realpart(omega)))
        
        ! allocate variables holding X, Y and Z of plot
        allocate(X_plot(n_theta_plot,n_zeta_plot,grp_n_r_X,product(n_t)))
        allocate(Y_plot(n_theta_plot,n_zeta_plot,grp_n_r_X,product(n_t)))
        allocate(Z_plot(n_theta_plot,n_zeta_plot,grp_n_r_X,product(n_t)))
        
        ! initialize theta_E and zeta_E
        if (follow_B) then
            ierr = 1
            CHCKERR('NOT YET IMPLEMENTED')
            !!!! SHOULD BE SIMILAR TO THE  GRID PLOT, WITH ADDITIONALLY THE !!!!
            !!!! SOLUTION VECTOR SUPERIMPOSED ALONG THE GRID LINES          !!!!
        else
            ! normal variable taken from grp_r_X
            allocate(r_F(grp_n_r_X))
            r_F = grp_r_X(1:grp_n_r_X)
            ! theta equidistant
            allocate(theta_E(n_theta_plot,n_zeta_plot,grp_n_r_X))
            if (n_theta_plot.eq.1) then
                theta_E = 0.0_dp
            else
                do id = 1,n_theta_plot
                    theta_E(id,:,:) = &
                        &pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)              ! starting from pi gives nicer plots
                end do
            end if
            ! zeta equidistant
            allocate(zeta_E(n_theta_plot,n_zeta_plot,grp_n_r_X))
            if (n_zeta_plot.eq.1) then
                zeta_E = 0.0_dp
            else
                do id = 1,n_zeta_plot
                    zeta_E(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
                end do
            end if
            ! convert  the perturbation  normal  values  grp_r_X to  equilibrium
            ! normal values
            allocate(r_E(grp_n_r_X))
            ierr = coord_F2E(grp_r_X(1:grp_n_r_X),r_E)
            CHCKERR('')
            ! convert all equilibrium coordinates to flux coordinates
            allocate(theta_F(n_theta_plot,n_zeta_plot,grp_n_r_X))
            allocate(zeta_F(n_theta_plot,n_zeta_plot,grp_n_r_X))
            ierr = coord_E2F(r_E,theta_E,zeta_E,r_F,theta_F,zeta_F)
            CHCKERR('')
        end if
        
        ! calculate  individual X,Y  and  Z using  the  Equilibrium theta_E  and
        ! zeta_plot, tabulated in the equilibrium normal grid
        ierr = calc_XYZ_grid(r_E,theta_E,zeta_E,X_plot_ind,Y_plot_ind,&
            &Z_plot_ind)
        CHCKERR('')
        
        ! calculate the function to plot: global mode
        allocate(f_plot(n_theta_plot,n_zeta_plot,grp_n_r_X,product(n_t)))
        
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
            ierr = calc_real_X(X_vec,omega**2,r_F,theta_F,zeta_F,time_frac,&
                &f_plot(:,:,:,id))
            CHCKERR('')
            
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
        deallocate(theta_E,zeta_E,r_E)
        deallocate(theta_F,zeta_F,r_F)
        
        ! user output
        call lvl_ud(-1)
        call writo('Finished the plot of the Eigenfunction')
    end function plot_X_vec_HDF5
end module
