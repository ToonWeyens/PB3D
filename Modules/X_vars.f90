!------------------------------------------------------------------------------!
!   Holds variables pertaining to the perturbation quantities
!------------------------------------------------------------------------------!
module X_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, max_str_ln, iu
    use message_ops, only: lvl_ud, writo, print_ar_1
    use output_ops, only: print_GP_2D, draw_GP
    use str_ops, only: r2strt, i2str
    
    implicit none
    
    private
    public init_X, calc_PV, calc_KV, calc_U, calc_extra, &
        &dealloc_X, dealloc_X_final, init_m, calc_V_int, &
        &rho, n_x, m_X, n_r_X, U_X_0, U_X_1, DU_X_0, DU_X_1, extra1, grp_r_X, &
        &extra2, extra3, PV0, PV1, PV2, KV0, KV1, KV2, X_vec, X_val, min_m_X, &
        &max_m_X, grp_n_r_X, size_X, PV_int, KV_int, grp_min_r_X, grp_max_r_X, &
        &min_n_X, max_n_X
    
    ! global variables
    real(dp), allocatable :: rho(:)                                             ! density
    integer :: size_X                                                           ! size of n_X and m_X (nr. of modes)
    integer :: min_n_X                                                          ! lowest poloidal mode number m_X
    integer :: max_n_X                                                          ! highest poloidal mode number m_X
    integer, allocatable :: n_X(:)                                              ! vector of poloidal mode numbers
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    integer, allocatable :: m_X(:)                                              ! vector of poloidal mode numbers
    integer :: n_r_X                                                            ! number of normal points for perturbations
    integer :: grp_n_r_X                                                        ! number of normal points for perturbations on this process
    integer :: grp_min_r_X, grp_max_r_X                                         ! min. and max. r range of this process in alpha group
    real(dp), allocatable :: grp_r_X(:)                                         ! normal points in Flux coords., globally normalized to (min_r_X..max_r_X)
    complex(dp), allocatable :: U_X_0(:,:,:), U_X_1(:,:,:)                      ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
    complex(dp), allocatable :: DU_X_0(:,:,:), DU_X_1(:,:,:)                    ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
    real(dp), allocatable :: mu0sigma(:,:)                                      ! parallel current
    real(dp), allocatable :: extra1(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra2(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra3(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    complex(dp), allocatable :: PV0(:,:,:,:)                                    ! ~PV^0 coefficient
    complex(dp), allocatable :: PV1(:,:,:,:)                                    ! ~PV^1 coefficient
    complex(dp), allocatable :: PV2(:,:,:,:)                                    ! ~PV^2 coefficient
    complex(dp), allocatable :: KV0(:,:,:,:)                                    ! ~KV^0 coefficient
    complex(dp), allocatable :: KV1(:,:,:,:)                                    ! ~KV^1 coefficient
    complex(dp), allocatable :: KV2(:,:,:,:)                                    ! ~KV^2 coefficient
    complex(dp), allocatable :: PV_int(:,:,:,:)                                 ! <~PV e^i(k-m)ang_par_F> coefficient
    complex(dp), allocatable :: KV_int(:,:,:,:)                                 ! <~KV e^i(k-m)ang_par_F> coefficient
    complex(dp), allocatable :: X_vec(:,:,:)                                    ! Eigenvector solution, with ghost region
    complex(dp), allocatable :: X_val(:)                                        ! Eigenvalue solution
    complex(dp), allocatable :: exp_ang_par_F(:,:,:,:)                          ! exp(i (k-m) ang_par_F)
    
contains
    ! initialize the variable m and check and/or plot it
    integer function init_m() result(ierr)
        use num_vars, only: plot_jq, use_pol_flux
        
        character(*), parameter :: rout_name = 'init_m'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set size_X
        if (use_pol_flux) then
            size_X = max_m_X - min_m_X + 1
        else
            size_X = max_n_X - min_n_X + 1
        end if
        
        ! setup n_X and m_X
        allocate(n_X(size_X),m_X(size_X))
        if (use_pol_flux) then
            n_X = min_n_X
            m_X = [(id, id = min_m_X, max_m_X)]
        else
            n_X = [(id, id = min_n_X, max_n_X)]
            m_X = min_m_X
        end if
        
        ! plot resonances if requested
        if (plot_jq) then
            call resonance_plot
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
            use eq_vars, only: q_saf_V_full, rot_t_V_full                       ! q_saf_V_full and rot_t_V_full are NOT normalized
            use num_vars, only: glb_rank, use_pol_flux
            
            character(*), parameter :: rout_name = 'check_m_and_n'
            
            ! local variables
            integer :: id                                                       ! counter
            real(dp) :: tol = 0.1                                               ! tolerance for being out of range of q or iota values
            real(dp) :: min_jq, max_jq                                          ! min. and max. values of q or iota
            
            ! initialize ierr
            ierr = 0
            
            if (glb_rank.eq.0) then
                call writo('Checking mode numbers')
                ! set min_jq and max_jq in flux coordinate system
                if (use_pol_flux) then
                    min_jq = minval(-q_saf_V_full(:,0))
                    max_jq = maxval(-q_saf_V_full(:,0))
                else
                    min_jq = minval(-rot_t_V_full(:,0))
                    max_jq = maxval(-rot_t_V_full(:,0))
                end if
                
                min_jq = min_jq - tol*abs(min_jq)
                max_jq = max_jq + tol*abs(max_jq)
                
                ! for every mode (n,m) check whether  m/n is inside the range of
                ! q values or n/m inside the range of iota values
                do id = 1, size_X
                    if (use_pol_flux) then
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
    
    ! initialize the variables of this module
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
        use eq_vars, only: n_par, grp_n_r_eq, ang_par_F
        use num_vars, only: use_pol_flux
        
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
        if (use_pol_flux) then
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
    
    ! plot  q-profile  or iota-profile  in  flux coordinates  with nq-m  = 0  or
    ! n-iotam = 0 indicate if requested
    subroutine resonance_plot
        use num_vars, only: glb_rank, tol_NR, use_pol_flux
        use utilities, only: calc_zero_NR, interp_fun_1D
        use eq_vars, only: q_saf_V_full, rot_t_V_full, &                        ! q_saf_V_full and rot_t_V_full are NOT normalized
            &flux_p_V_full, flux_t_V_full
        use VMEC_vars, only: n_r_eq
        
        ! local variables (also used in child functions)
        real(dp) :: mnfrac_for_function                                         ! fraction m/n or n/m to determine resonant flux surface
        real(dp), allocatable :: jq_for_function(:,:)                           ! rot_t or q_saf
        
        ! local variables (not to be used in child functions)
        integer :: jd, kd                                                       ! counters
        real(dp), allocatable :: x_vars(:,:)                                    ! for plotting
        real(dp), allocatable :: y_vars(:,:)                                    ! for plotting
        real(dp) :: jq_solution                                                 ! solution for q = m/n or iota = n/m
        real(dp) :: jq_solution_transf                                          ! transformed solution to flux coordinates
        real(dp) :: old_tol_NR                                                  ! to backup tol_NR
        integer :: istat                                                        ! status
        character(len=max_str_ln) :: plot_title, file_name                      ! name of plot, of file
        
        if (glb_rank.eq.0) then
            ! backup tol_NR
            old_tol_NR = tol_NR
            
            ! set  lower  tolerance  for  Newton-Rhapson  because  working  with
            ! discrete data that is not perfectly continous
            tol_NR = 5E-3_dp
            
            if (use_pol_flux) then
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
            
            if (use_pol_flux) then
                x_vars(:,1) = flux_p_V_full(:,0)
                y_vars(:,1) = -q_saf_V_full(:,0)                                ! convert to flux coordinates
                jq_for_function = -q_saf_V_full(:,0:2)                          ! (need -1 to convert form Flux to VMEC)
            else
                x_vars(:,1) = flux_t_V_full(:,0)
                y_vars(:,1) = -rot_t_V_full(:,0)                                ! convert to flux coordinates
                jq_for_function = -rot_t_V_full(:,0:2)                          ! (need -1 to convert form Flux to VMEC)
            end if
            
            kd = 2
            do jd = 1, size_X
                ! find place where q = m/n or  iota = n/m in VMEC coordinates by
                ! solving q-m/n = 0 or iota-n/m=0, using the functin jq_fun
                call lvl_ud(1)
                
                ! set up mnfrac_for_function
                if (use_pol_flux) then
                    mnfrac_for_function = m_X(jd)*1.0_dp/n_X(jd)
                else
                    mnfrac_for_function = n_X(jd)*1.0_dp/m_X(jd)
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
                        if (use_pol_flux) then
                            y_vars(n_r_eq,kd) = -q_saf_V_full(n_r_eq,0)         ! convert to flux coordinates
                        else
                            y_vars(n_r_eq,kd) = -rot_t_V_full(n_r_eq,0)         ! convert to flux coordinates
                        end if
                    else
                        call writo('Mode (n,m) = ('//trim(i2str(n_X(jd)))//','&
                            &//trim(i2str(m_X(jd)))//') resonates in plasma &
                            &at normalized flux surface '//&
                            &trim(r2strt(jq_solution)))
                        if (use_pol_flux) then
                            y_vars(n_r_eq,kd) = m_X(jd)*1.0_dp/n_X(jd)
                        else
                            y_vars(n_r_eq,kd) = n_X(jd)*1.0_dp/m_X(jd)
                        end if
                    end if
                    
                    ! convert solution to flux coordinates
                    if (use_pol_flux) then
                        istat = interp_fun_1D(jq_solution_transf,&
                            &flux_p_V_full(:,0),jq_solution)
                    else
                        istat = interp_fun_1D(jq_solution_transf,&
                            &flux_t_V_full(:,0),jq_solution)
                    end if
                    
                    x_vars(:,kd) = jq_solution_transf
                    kd = kd + 1
                end if
            end do
            
            ! check status
            if (istat.ne.0) then
                call writo('WARNING: Failed to produce plot')
                return
            end if
            
            ! rescale x axis
            if (use_pol_flux) then
                x_vars = x_vars/flux_p_V_full(n_r_eq,0)
            else
                x_vars = x_vars/flux_t_V_full(n_r_eq,0)
            end if
            
            call lvl_ud(-1)
            
            call writo('Plotting results')
            if (use_pol_flux) then
                plot_title = 'safety factor q'
                file_name = 'q_saf.dat'
            else
                plot_title = 'rotational transform iota'
                file_name = 'rot_t.dat'
            end if
            ! plot on screen
            call print_GP_2D(plot_title,file_name,y_vars(:,1:kd-1),&
                &x=x_vars(:,1:kd-1))
            ! plot in file as well
            call draw_GP(plot_title,file_name,kd-1,.true.,.false.)
            
            call lvl_ud(-1)
            call writo('Done plotting')
            
            tol_NR = old_tol_NR
        end if
    contains
        ! returns q-m/n or  iota-n/m in VMEC coordinates, used to  solve for q =
        ! m/n or iota = n/m
        real(dp) function jq_fun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(n_r_eq)                                           ! so interp_fun_1D can be used
            
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
                ! interpolate using interp_fun_1D
                varin = jq_for_function(:,0) - mnfrac_for_function
                istat = interp_fun_1D(res,varin,pt)
            end if
        end function jq_fun
        
        ! returns d(q-m/n)/dr VMEC coordinates, used to solve for q = m/n
        ! (the toroidal mode number in VMEC coordinates is equal to -n_X)
        ! WARNING: This  routine requires that jq_for_function's  derivatives be
        ! calculated up to order 2. This is NOT checked!
        real(dp) function jq_dfun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(n_r_eq)                                           ! so interp_fun_1D can be used
            
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
                ! interpolate using interp_fun_1D
                varin = jq_for_function(:,1)
                istat = interp_fun_1D(res,varin,pt)
            end if
        end function jq_dfun
    end subroutine resonance_plot
    
    ! calculate rho from user input
    subroutine calc_rho
        use eq_vars, only: pres_V, grp_n_r_eq, grp_min_r_eq
        use VMEC_vars, only: gam
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp) :: expon                                                       ! exponent = 1/gam
        
        ! allocate rho
        allocate(rho(grp_n_r_eq))
        
        ! set exp
        expon = 1.0_dp/gam
        
        ! loop over all normal points
        do kd = 1,grp_n_r_eq
            if (pres_V(kd,0).gt.0) then
                rho(kd) = pres_V(kd,0)**expon
            else
                call writo('WARNING: pressure was negative at point '//&
                    &trim(i2str(grp_min_r_eq+kd))//'/'//&
                    &trim(i2str(grp_min_r_eq+grp_n_r_eq)),persistent=.true.)
                rho(kd) = rho(kd-1)
            end if
        end do
    end subroutine calc_rho
    
    ! calculate  ~PV_(k,m)^i  (pol.  flux)  or ~PV_(l,n)^i  (tor.  flux) at  all
    ! grp_n_r_eq values
    ! (see [ADD REF] for details)
    subroutine calc_PV
        use metric_ops, only: g_FD, h_FD, jac_FD
        use eq_vars, only: n_par, q_saf_FD, rot_t_FD, grp_n_r_eq
        use num_vars, only: use_pol_flux
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:)                                   ! common factor |nabla psi|^2/(J^2*B^2)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        
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
        
        ! set up common factor for PVi
        allocate(com_fac(n_par,grp_n_r_eq))
        com_fac = h22/g33
        
        ! set up fac_n and fac_m
        allocate(fac_n(grp_n_r_eq),fac_m(grp_n_r_eq))
        if (use_pol_flux) then
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
                    &mu0sigma/J * (extra1 + extra2) - extra3
                
                ! add (nq-k)*(nq-m)/(J^2 |nabla psi|^2) to PV0
                do kd = 1,grp_n_r_eq
                    PV0(:,kd,k,m) = PV0(:,kd,k,m) + &
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
        
        ! deallocate mu0sigma and the extra's
        deallocate(mu0sigma,extra1,extra2,extra3)
    end subroutine calc_PV
    
    ! calculate  ~KV_(k,m)^i  (pol.  flux)  or ~KV_(l,n)^i  (tor.  flux) at  all
    ! grp_n_r_eq values
    ! (see [ADD REF] for details)
    subroutine calc_KV
        use metric_ops, only: g_FD, h_FD, jac_FD
        use eq_vars, only: n_par, grp_n_r_eq
        
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
                    &1/h22
                
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
            KV2(:,kd,:,:) = KV2(:,kd,:,:)*rho(kd)
            KV1(:,kd,:,:) = KV1(:,kd,:,:)*rho(kd)
        end do
    end subroutine calc_KV
    
    ! calculate U_m^0, U_m^1 or U_n^0, U_n^1  at grp_n_r_eq values of the normal
    ! coordinate, n_par values  of the parallel coordinate and  size_X values of
    ! the poloidal mode number,  or of the toroidal mode number,  as well as the
    ! poloidal derivatives
    ! Note: The  poloidal  derivatives  have  the  factor  i/n  or i/m  included
    ! already, as opposed to [ADD REF]
    subroutine calc_U
        use num_vars, only: use_pol_flux
        use metric_ops, only: g_FD, h_FD, jac_FD
        use eq_vars, only: rot_t_FD, q_saf_FD, ang_par_F, n_par, p => pres_FD, &
            &grp_n_r_eq
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        real(dp) :: n_frac                                                      ! nq-m (pol. flux) or n-iotam (tor. flux)
        real(dp), allocatable :: djq(:)                                         ! either q' (pol. flux) or -iota' (tor. flux)
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        real(dp), allocatable :: mn(:)                                          ! either n*A_0 (pol. flux) or m (tor.flux)
        
        ! submatrices
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
        real(dp), allocatable :: D3h22(:,:)                                     ! D_theta h^psi,psi
        
        ! set up djq, fac_n, fac_m and mn
        allocate(djq(grp_n_r_eq),fac_n(grp_n_r_eq),fac_m(grp_n_r_eq),mn(size_X))
        if (use_pol_flux) then
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
        allocate(D3h22(n_par,grp_n_r_eq)); D3h22 = h_FD(:,:,2,2,0,0,1)
        
        ! loop over the M elements of U_X and DU
        do jd = 1,size(m_X)
            ! loop over all normal points
            do kd = 1,grp_n_r_eq
                n_frac = n_X(jd)*fac_n(kd)-m_X(jd)*fac_m(kd)
                U_X_0(:,kd,jd) = &
                    &-(h12(:,kd)/h22(:,kd) + djq(kd)*ang_par_F(:,kd))&
                    &+ iu/(mn(jd)*g33(:,kd)) * (g13(:,kd)*djq(kd) + &
                    &J(:,kd)**2*p(kd,1) + iu*n_frac * &
                    &( g13(:,kd)*djq(kd)*ang_par_F(:,kd) - g23(:,kd) ))
                DU_X_0(:,kd,jd) = -(D3h12(:,kd)/h22(:,kd) - &
                    &D3h22(:,kd)*h12(:,kd)/h22(:,kd)**2 + djq(kd)) - &
                    &iu*D3g33(:,kd)/(mn(jd)*g33(:,kd)**2) * (g13(:,kd)*djq(kd)&
                    &+ J(:,kd)**2*p(kd,1) + iu*n_frac * &
                    &( g13(:,kd)*djq(kd)*ang_par_F(:,kd) - g23(:,kd) )) + &
                    &iu/(mn(jd)*g33(:,kd)) * (D3g13(:,kd)*djq(kd) + &
                    &2*D3J(:,kd)*J(:,kd)*p(kd,1) + iu*n_frac * &
                    &( D3g13(:,kd)*djq(kd)*ang_par_F(:,kd) + g13(:,kd)*djq(kd) &
                    &- D3g23(:,kd) )) + iu*n_frac*U_X_0(:,kd,jd)
                U_X_1(:,kd,jd) = iu/mn(jd) * &
                    &(1 + n_frac/mn(jd) * g13(:,kd)/g33(:,kd))
                DU_X_1(:,kd,jd) = iu/mn(jd) * &
                    &(n_frac/mn(jd) * (D3g13(:,kd)/g33(:,kd) - &
                    &D3g33(:,kd)*g13(:,kd)/g33(:,kd)**2) + &
                    &iu*n_frac*U_X_1(:,kd,jd) )
            end do
            
            !write(*,*) 'U1'
            !call print_GP_2D('RE U1(:,:,'//trim(i2str(jd))//')','',&
                !&realpart(transpose(U_X_1(:,:,jd))))
            !call print_GP_2D('IM U1(:,:,'//trim(i2str(jd))//')','',&
                !&imagpart(transpose(U_X_1(:,:,jd))))
            
            !write(*,*) 'DU1'
            !call print_GP_2D('RE DU1(:,:,'//trim(i2str(jd))//')','',&
                !&realpart(transpose(DU_X_1(:,:,jd))))
            !call print_GP_2D('IM DU1(:,:,'//trim(i2str(jd))//')','',&
                !&imagpart(transpose(DU_X_1(:,:,jd))))
        end do
    end subroutine calc_U
    
    ! calculate extra1, extra2 and extra3
    !   extra1 = S*J
    !   extra2 = mu0sigma*J*B^2/h^psi,psi
    !   extra3 = 2*p'*kn
    ! with
    !   shear S = - d Theta^alpha/d_theta 1 / J
    !   parallel current mu0sigma = B . nabla x B / B^2
    !   normal curvature kn = nabla psi / h^psi,psi * nabla (mu0 p + B^2/2)
    subroutine calc_extra
        use metric_ops, only: g_FD, h_FD, jac_FD
        use eq_vars, only: n_par, p => pres_FD, grp_n_r_eq
        
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
        
        ! allocate mu0sigma and extra terms. Later they are deallocated in calc_PV
        allocate(mu0sigma(n_par,grp_n_r_eq))
        allocate(extra1(n_par,grp_n_r_eq))
        allocate(extra2(n_par,grp_n_r_eq))
        allocate(extra3(n_par,grp_n_r_eq))
        
        ! calculate mu0sigma
        mu0sigma = 1/(J*g33) * (g13*(D2g33-D3g23) + g23*(D3g13-D1g33) &
            &+ g33*(D1g23-D2g13))
        
        ! calculate extra1
        extra1 = -D3h12/h22 + D3h22*h12/h22**2
        
        ! calculate extra2
        extra2 = g33/h22 * mu0sigma / J
        
        ! calculate extra3
        do kd = 1,grp_n_r_eq
            extra3(:,kd) = p(kd,1) * ( 2*J(:,kd)**2*p(kd,1)/g33(:,kd) + &
                &1/h22(:,kd) * ( &
                &h12(:,kd) * ( D1g33(:,kd)/g33(:,kd) - 2*D1J(:,kd)/J(:,kd) ) + &
                &h22(:,kd) * ( D2g33(:,kd)/g33(:,kd) - 2*D2J(:,kd)/J(:,kd) ) + &
                &h23(:,kd) * ( D3g33(:,kd)/g33(:,kd) - 2*D3J(:,kd)/J(:,kd) ) ) )
        end do
        
        ! calculate rho from input
        call calc_rho
    end subroutine calc_extra
    
    ! calculates  magnetic  integral  <V e^[i(k-m)ang_par_F]>,  defined  as  the
    ! matrix  
    !   <V e^[i(k-m)ang_par_F]> = [ oint J V(k,m) e^i(k-m)ang_par_F dang_par_F ]
    ! or
    !   <V e^[i(n-l)ang_par_F]> = [ oint J V(l,n) e^i(n-l)ang_par_F dang_par_F ]
    ! Note: V is changed on exit
    subroutine calc_V_int(V,V_int)
        use eq_vars, only: n_par, grp_n_r_eq, ang_par_F
        use metric_ops, only: jac_FD
        use utilities, only: calc_int
        
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
    
    ! deallocates  perturbation quantities that  are not used anymore  after the
    ! equilibrium phase
    subroutine dealloc_X
        deallocate(rho)
        deallocate(U_X_0,U_X_1,DU_X_0,DU_X_1)
        deallocate(PV0,PV1,PV2,KV0,KV1,KV2)
        deallocate(exp_ang_par_F)
    end subroutine dealloc_X
    
    ! deallocates  perturbation quantities that  are not used anymore  after the
    ! calculations for a certain alpha
    subroutine dealloc_X_final
        deallocate(n_X,m_X)
        deallocate(X_vec,X_val)
        deallocate(PV_int,KV_int)
        deallocate(grp_r_X)
    end subroutine dealloc_X_final
end module
