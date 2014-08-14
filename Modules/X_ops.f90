!------------------------------------------------------------------------------!
!   Calculate the matrix  elements due to the Plasma of  the system of equations
!   that has to be solved as described in [ADD REF]
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, mu_0, max_str_ln
    use output_ops, only: lvl_ud, writo, print_ar_2, print_GP_2D
    use str_ops, only: i2str, r2strt

    implicit none
    private
    public prepare_matrix_X, solve_EV_system

contains
    ! initialize the variables of this module
    ! ¡¡¡ THIS SHOULD BE  CHANGED TO TAKE INTO ACCOUNT THAT  U_X AND DU_X ARE
    ! HERMITIAN SO ONLY  M*(M+1)/2 ELEMENTS ARE NEEDED INSTEAD  OF M^2. HOWEVER,
    ! TAKE INTO ACCOUNT AS WELL THE MPI STORAGE MATRIX FORMATS !!!
    integer function init_X_ops() result(ierr)
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        use X_vars, only: U_X_0, U_X_1, DU_X_0, DU_X_1, sigma, extra1, extra2, &
            &extra3, PV0, PV1, PV2, KV0, KV1, KV2, min_m_X, max_m_X, m_X
        use num_vars, only: plot_q
        
        character(*), parameter :: rout_name = 'init_X_ops'
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_m_X                                                        ! number of poloidal modes m
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        
        ! initialize ierr
        ierr = 0
        
        ! initialize
        n_m_X = max_m_X - min_m_X + 1
        if (n_m_X.lt.1) then
            ierr = 1
            err_msg = 'max_m_X has to be larger or equal to min_m_X'
            CHCKERR(err_msg)
        end if
        
        ! setup m_X
        if (allocated(m_X)) deallocate(m_X)
        allocate(m_X(n_m_X))
        m_X = [(id, id = min_m_X, max_m_X)]
        if (plot_q) then
            call resonance_plot
        else
            call writo('resonance plot not requested')
        end if
        
        ! U_X_0
        if (allocated(U_X_0)) deallocate(U_X_0)
        allocate(U_X_0(n_par,n_r,n_m_X))
        
        ! U_X_1
        if (allocated(U_X_1)) deallocate(U_X_1)
        allocate(U_X_1(n_par,n_r,n_m_X))
        
        ! DU_X_0
        if (allocated(DU_X_0)) deallocate(DU_X_0)
        allocate(DU_X_0(n_par,n_r,n_m_X))
        
        ! DU_X_1
        if (allocated(DU_X_1)) deallocate(DU_X_1)
        allocate(DU_X_1(n_par,n_r,n_m_X))
        
        ! sigma
        if (allocated(sigma)) deallocate(sigma)
        allocate(sigma(n_par,n_r))
        
        ! extra terms (see calc_extra)
        if (allocated(extra1)) deallocate(extra1)
        allocate(extra1(n_par,n_r))
        if (allocated(extra2)) deallocate(extra2)
        allocate(extra2(n_par,n_r))
        if (allocated(extra3)) deallocate(extra3)
        allocate(extra3(n_par,n_r))
        
        ! PV0
        if (allocated(PV0)) deallocate(PV0)
        allocate(PV0(n_par,n_r,n_m_X,n_m_X))
        
        ! PV1
        if (allocated(PV1)) deallocate(PV1)
        allocate(PV1(n_par,n_r,n_m_X,n_m_X))
        
        ! PV2
        if (allocated(PV2)) deallocate(PV2)
        allocate(PV2(n_par,n_r,n_m_X,n_m_X))
        
        ! KV0
        if (allocated(KV0)) deallocate(KV0)
        allocate(KV0(n_par,n_r,n_m_X,n_m_X))
        
        ! KV1
        if (allocated(KV1)) deallocate(KV1)
        allocate(KV1(n_par,n_r,n_m_X,n_m_X))
        
        ! KV2
        if (allocated(KV2)) deallocate(KV2)
        allocate(KV2(n_par,n_r,n_m_X,n_m_X))
    end function init_X_ops
    
    ! plot q-profile with nq-m = 0 indicated if requested
    subroutine resonance_plot
        use num_vars, only: glob_rank, tol_NR
        use utilities, only: calc_zero_NR, calc_interp
        use eq_vars, only: q_saf
        use VMEC_vars, only: n_r
        use X_vars, only: min_m_X, max_m_X, m_X, n_X
        
        ! local variables (also used in child functions)
        integer :: m_X_for_function                                             ! used in q_fun to determine flux surface at which q = m/n
        
        ! local variables (not to be used in child functions)
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: x_vars(:,:)                                    ! for plotting
        real(dp), allocatable :: y_vars(:,:)                                    ! for plotting
        real(dp) :: q_solution                                                  ! solution for q = m/n
        real(dp) :: old_tol_NR                                                  ! to backup tol_NR
        integer :: istat                                                        ! status
        integer :: n_m_X                                                        ! number of poloidal modes m
        
        if (glob_rank.eq.0) then
            ! initialize
            n_m_X = max_m_X - min_m_X + 1
            
            ! backup tol_NR
            old_tol_NR = tol_NR
            
            ! set  lower  tolerance  for  Newton-Rhapson  because  working  with
            ! discrete data that is not perfectly continous
            tol_NR = 5E-3_dp
            
            call writo('plotting safety factor q and resonant surfaces q = m/n')
            call lvl_ud(1)
            
            call writo('calculating resonant surfaces q = m/n')
            call lvl_ud(1)
            
            allocate(x_vars(n_r,n_m_X+1))
            allocate(y_vars(n_r,n_m_X+1))
            x_vars(:,1) = [((id-1.0_dp)/(n_r-1.0_dp), id = 1,n_r)]
            !y_vars = minval(q_saf(:,0))
            y_vars = 0
            y_vars(:,1) = q_saf(:,0)
            
            kd = 2
            do jd = 1, n_m_X
                ! find place where q = m/n by solving q-m/n = 0, using
                ! the functin q_fun with m_X_for_function = m_X(jd)
                m_X_for_function = m_X(jd)
                call lvl_ud(1)
                istat = calc_zero_NR(q_solution,q_fun,q_dfun,1.0_dp)
                call lvl_ud(-1)
                ! intercept error
                if (istat.ne.0) then
                    call writo('Error intercepted: Couldn''t find resonating &
                        &surface for poloidal mode '//trim(i2str(m_X(jd))))
                else
                    if (q_solution.gt.1.0_dp) then
                        call writo('Poloidal mode '//trim(i2str(m_X(jd)))//&
                            &' does not resonate in plasma')
                        y_vars(n_r,kd) = maxval(q_saf(:,0))
                    else
                        call writo('Poloidal mode '//trim(i2str(m_X(jd)))//&
                            &' resonates in plasma at normalized &
                            &radius '//trim(r2strt(q_solution)))
                        y_vars(n_r,kd) = m_X(jd)*1.0_dp/n_X
                    end if
                    x_vars(:,kd) = q_solution
                    kd = kd + 1
                end if
            end do
            
            call lvl_ud(-1)
            
            call writo('Plotting results')
            call print_GP_2D('safety factor q','',y_vars(:,1:kd-1),&
                &x=x_vars(:,1:kd-1))
            
            call lvl_ud(-1)
            
            tol_NR = old_tol_NR
        end if
    contains
        ! returns q-m/n, used to solve for q = m/n
        real(dp) function q_fun(pt) result(res)
            use utilities, only: calc_interp
            
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(1,n_r,1,1)                                        ! so calc_interp can be used
            real(dp) :: varout(1,1,1)                                           ! so calc_interp can be used
            real(dp) :: pt_copy                                                 ! copy of pt to use with calc_interp
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.0) then                                                   ! point requested lower than 0
                ! extrapolate variable from 0
                res = q_saf(1,0) - m_X_for_function*1.0_dp/n_X + &
                    &q_saf(0,1) * (pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from 1
                res = q_saf(n_r,0) - m_X_for_function*1.0_dp/n_X + &
                    &q_saf(n_r,1) * (pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using calc_interp
                pt_copy = pt
                
                varin(1,:,1,1) = q_saf(:,0) - m_X_for_function*1.0_dp/n_X
                istat = calc_interp(varin,varout,pt_copy)
                
                res = varout(1,1,1)
            end if
        end function q_fun
        
        ! returns d(q-m/n)/dr, used to solve for q = m/n
        real(dp) function q_dfun(pt) result(res)
            use utilities, only: calc_interp
            
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            real(dp) :: varin(1,n_r,1,1)                                        ! so calc_interp can be used
            real(dp) :: varout(1,1,1)                                           ! so calc_interp can be used
            real(dp) :: pt_copy                                                 ! copy of pt to use with calc_interp
            
            ! initialize res
            res = 0
            
            ! check whether to interpolate or extrapolate
            if (pt.lt.0) then                                                   ! point requested lower than 0
                ! extrapolate variable from 0
                res = q_saf(1,1) + q_saf(0,2) * (pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from 1
                res = q_saf(n_r,1) + q_saf(n_r,2) * (pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using calc_interp
                pt_copy = pt
                
                varin(1,:,1,1) = q_saf(:,1)
                istat = calc_interp(varin,varout,pt_copy)
                
                res = varout(1,1,1)
            end if
        end function q_dfun
    end subroutine resonance_plot
        
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    integer function prepare_matrix_X() result(ierr)
        character(*), parameter :: rout_name = 'prepare_matrix_X'
        
        ! initialize ierr
        ierr = 0
        
        ! initialize the variables
        call writo('Initalizing variables...')
        call lvl_ud(1)
        ierr = init_X_ops()
        CHCKERR('')
        call lvl_ud(-1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        call calc_U
        call lvl_ud(-1)
        
        ! calculate extra equilibrium quantities
        call writo('Calculating extra equilibrium quantities...')
        call lvl_ud(1)
        call calc_extra
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
    end function prepare_matrix_X
    
    ! calculate ~PV_(k,m)^i at all n_r values (see [ADD REF] for details)
    subroutine calc_PV
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par, q => q_saf_FD
        use VMEC_vars, only: n_r
        use X_vars, only: DU_X_0, DU_X_1, sigma, extra1, extra2, extra3, PV0, &
            &PV1, PV2, n_X, m_X
        
        ! local variables
        integer :: n_m_X                                                        ! number of poloidal modes m
        integer :: m, k, kd                                                     ! counters
        real(dp), allocatable :: com_fac(:,:)                                   ! common factor |nabla psi|^2/(mu_0*J^2*B^2)
        
        ! submatrices
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        ! lower metric factors
        real(dp), allocatable :: g33(:,:)                                       ! h^alpha,psi
        ! upper metric factors
        real(dp), allocatable :: h22(:,:)                                       ! h^alpha,psi
        
        ! initialize
        n_m_X = size(m_X)
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,n_r)); J = jac_F_FD(:,:,0,0,0)
        ! lower metric factors
        allocate(g33(n_par,n_r)); g33 = g_F_FD(:,:,3,3,0,0,0)
        ! upper metric factors
        allocate(h22(n_par,n_r)); h22 = h_F_FD(:,:,2,2,0,0,0)
        
        ! set up common factor for PVi
        allocate(com_fac(n_par,n_r))
        com_fac = h22/(mu_0*J*g33)
        
        ! calculate PVi
        do m = 1,n_m_X
            do k = 1,m
                ! calculate PV0
                PV0(:,:,k,m) = com_fac * (DU_X_0(:,:,m) - extra1 - extra2 ) * &
                    &(conjg(DU_X_0(:,:,k)) - extra1 - extra2) - &
                    &sigma/J * (extra1 + extra2) - extra3
                
                ! add (nq-k)*(nq-m)/(mu_0 J^2 |nabla psi|^2) to PV0
                do kd = 1,n_r
                    PV0(:,kd,k,m) = PV0(:,kd,k,m) + &
                        &(n_X*q(kd,0)-m_X(m))*(n_X*q(kd,0)-m_X(k)) / &
                        &( mu_0*J(:,kd)**2*h22(:,kd) )
                end do
                
                ! calculate PV2
                PV2(:,:,k,m) = com_fac * DU_X_1(:,:,m) * conjg(DU_X_1(:,:,k))
            end do
        end do
        
        ! fill the Hermitian conjugates
        do m = 1,n_m_X
            do k = m+1,n_m_X
                PV0(:,:,k,m) = conjg(PV0(:,:,m,k))
                PV2(:,:,k,m) = conjg(PV2(:,:,m,k))
            end do
        end do
        
        ! PV1 is not Hermitian
        do m = 1,n_m_X
            do k = 1,n_m_X
                ! calculate PV1
                PV1(:,:,k,m) = com_fac * DU_X_1(:,:,m) * &
                    &(conjg(DU_X_0(:,:,k)) - extra1 - extra2)
            end do
        end do
    end subroutine
    
    ! calculate ~KV_(k,m)^i at all n_r values (see [ADD REF] for details)
    subroutine calc_KV
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        use X_vars, only: rho, U_X_0, U_X_1, KV0, KV1, KV2, m_X
        
        
        ! local variables
        integer :: m, k                                                         ! counters
        real(dp), allocatable :: com_fac(:,:)                                   ! common factor |nabla psi|^2/(mu_0*J^2*B^2)
        integer :: n_m_X                                                        ! number of poloidal modes m
        
        ! submatrices
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        ! lower metric factors
        real(dp), allocatable :: g33(:,:)                                       ! h^alpha,psi
        ! upper metric factors
        real(dp), allocatable :: h22(:,:)                                       ! h^alpha,psi
        
        ! initialize
        n_m_X = size(m_X)
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,n_r)); J = jac_F_FD(:,:,0,0,0)
        ! lower metric factors
        allocate(g33(n_par,n_r)); g33 = g_F_FD(:,:,3,3,0,0,0)
        ! upper metric factors
        allocate(h22(n_par,n_r)); h22 = h_F_FD(:,:,2,2,0,0,0)
        
        ! set up common factor
        allocate(com_fac(n_par,n_r))
        com_fac = J*h22/g33
        
        ! for  Hermitian KV0  and  KV2,  only half  of  the  terms  have  to  be
        ! calculated
        do m = 1,n_m_X
            do k = 1,m
                ! calculate KV0
                KV0(:,:,k,m) = com_fac * U_X_0(:,:,m) * conjg(U_X_0(:,:,k)) + &
                    &1/h22
                
                ! calculate KV2
                KV2(:,:,k,m) = com_fac * U_X_1(:,:,m) * conjg(U_X_1(:,:,k))
                
                ! multiply by rho
                KV0(:,:,k,m) = KV0(:,:,k,m)*rho
                KV2(:,:,k,m) = KV2(:,:,k,m)*rho
            end do
        end do
        ! fill the Hermitian conjugates
        do m = 1,n_m_X
            do k = m+1,n_m_X
                KV0(:,:,k,m) = conjg(KV0(:,:,m,k))
                KV2(:,:,k,m) = conjg(KV2(:,:,m,k))
            end do
        end do
        
        ! KV1 is not Hermitian
        do m = 1,n_m_X
            do k = 1,n_m_X
                ! calculate KV1
                KV1(:,:,k,m) = com_fac * U_X_1(:,:,m) * conjg(U_X_0(:,:,k))
                
                ! multiply by rho
                KV1(:,:,k,m) = KV1(:,:,k,m)*rho
            end do
        end do
    end subroutine
    
    ! set-up and  solve the  EV system  by discretizing  the equations  in n_r_X
    ! normal points,  making use of  PV0, PV1 and  PV2, interpolated in  the n_r
    ! (equilibrium) values
    integer function solve_EV_system() result(ierr)
        use X_vars, only: n_r_X, m_X, n_X, X_vec
        use VMEC_vars, only: n_r
        use num_vars, only: EV_style, group_rank
        use str_ops, only: i2str
        use slepc_ops, only: solve_EV_system_slepc
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! slepc solver for EV problem
                ierr = solve_EV_system_slepc(m_X,n_r_X,n_X,n_r-1._dp)
                CHCKERR('')
                ! TEMPORARY!!
                write(*,*) 'plotting results'
                if (group_rank.eq.0) then
                    call print_GP_2D('norm of solution','',&
                        &transpose(sqrt(realpart(X_vec(:,:,1))**2 + &
                        &imagpart(X_vec(:,:,1))**2)))
                end if
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
    
    ! calculate  U_m^0, U_m^1  at n_r  values  of the  normal coordinate,  n_par
    ! values  of the  parallel  coordinate and  M values  of  the poloidal  mode
    ! number, with m = size(m_X)
    subroutine calc_U
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: q => q_saf_FD, theta, n_par, p => pres_FD
        use VMEC_vars, only: n_r
        use X_vars, only: U_X_0, U_X_1, DU_X_0, DU_X_1, n_X, m_X
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        real(dp) :: n_frac                                                      ! = (nq-m)/n, to make notation easier
        
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
        
        ! set up submatrices
        ! jacobian
        allocate(J(n_par,n_r)); J = jac_F_FD(:,:,0,0,0)
        allocate(D3J(n_par,n_r)); D3J = jac_F_FD(:,:,0,0,1)
        ! lower metric factors
        allocate(g13(n_par,n_r)); g13 = g_F_FD(:,:,1,3,0,0,0)
        allocate(D3g13(n_par,n_r)); D3g13 = g_F_FD(:,:,1,3,0,0,1)
        allocate(g23(n_par,n_r)); g23 = g_F_FD(:,:,2,3,0,0,0)
        allocate(D3g23(n_par,n_r)); D3g23 = g_F_FD(:,:,2,3,0,0,1)
        allocate(g33(n_par,n_r)); g33 = g_F_FD(:,:,3,3,0,0,0)
        allocate(D3g33(n_par,n_r)); D3g33 = g_F_FD(:,:,3,3,0,0,1)
        ! upper metric factors
        allocate(h12(n_par,n_r)); h12 = h_F_FD(:,:,1,2,0,0,0)
        allocate(D3h12(n_par,n_r)); D3h12 = h_F_FD(:,:,1,2,0,0,1)
        allocate(h22(n_par,n_r)); h22 = h_F_FD(:,:,2,2,0,0,0)
        allocate(D3h22(n_par,n_r)); D3h22 = h_F_FD(:,:,2,2,0,0,1)
        
        ! loop over the M elements of U_X and DU
        do jd = 1,size(m_X)
            ! loop over all normal points
            do kd = 1,n_r
                n_frac = (n_X*q(kd,0)-m_X(jd))/n_X
                U_X_1(:,kd,jd) = 1 + n_frac * g13(:,kd)/g33(:,kd)
                DU_X_1(:,kd,jd) = n_frac * D3g13(:,kd)/g33(:,kd) - &
                    &D3g33(:,kd)*g13(:,kd)/g33(:,kd)**2 + &
                    &iu*n_frac*n_X*U_X_1(:,kd,jd)
                U_X_0(:,kd,jd) = -(h12(:,kd)/h22(:,kd) + q(kd,1)*theta(:,kd)) +&
                    &iu/(n_X*g33(:,kd)) * (g13(:,kd)*q(kd,1) + &
                    &J(:,kd)**2*mu_0*p(kd,1) + iu*n_frac*n_X * &
                    &( g13(:,kd)*q(kd,1)*theta(:,kd) - g23(:,kd) ))
                DU_X_0(:,kd,jd) = -(D3h12(:,kd)/h22(:,kd) - &
                    &D3h22(:,kd)*h12(:,kd)/h22(:,kd)**2 + q(kd,1)) - &
                    &iu*D3g33(:,kd)/(n_X*g33(:,kd)**2) * (g13(:,kd)*q(kd,1) + &
                    &J(:,kd)**2*mu_0*p(kd,1) + iu*n_frac*n_X * &
                    &( g13(:,kd)*q(kd,1)*theta(:,kd) - g23(:,kd) )) + &
                    &iu/(n_X*g33(:,kd)) * (D3g13(:,kd)*q(kd,1) + &
                    &2*D3J(:,kd)*J(:,kd)*mu_0*p(kd,1) + iu*n_frac*n_X * &
                    &( D3g13(:,kd)*q(kd,1)*theta(:,kd) + g13(:,kd)*q(kd,1) &
                    &- D3g23(:,kd) )) * iu*n_frac*n_X*U_X_0(:,kd,jd)
            end do
        end do
    end subroutine
    
    ! calculate extra1, extra2 and extra3
    !   extra1 = S*J
    !   extra2 = mu_0*sigma*J*B^2/h^psi,psi
    !   extra3 = 2*p'*kn
    ! with
    !   shear S = - d Theta^alpha/d_theta
    !   parallel current sigma = B . nabla x B / B^2
    !   normal curvature kn = nabla psi / h^psi,psi * nabla (p + B^2/2)
    subroutine calc_extra
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par, p => pres_FD
        use VMEC_vars, only: n_r
        use X_vars, only: calc_rho, sigma, extra1, extra2, extra3
        
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
        allocate(J(n_par,n_r)); J = jac_F_FD(:,:,0,0,0)
        allocate(D1J(n_par,n_r)); D1J = jac_F_FD(:,:,1,0,0)
        allocate(D2J(n_par,n_r)); D2J = jac_F_FD(:,:,0,1,0)
        allocate(D3J(n_par,n_r)); D3J = jac_F_FD(:,:,0,0,1)
        ! lower metric factors
        allocate(g13(n_par,n_r)); g13 = g_F_FD(:,:,1,3,0,0,0)
        allocate(D2g13(n_par,n_r)); D2g13 = g_F_FD(:,:,1,3,0,1,0)
        allocate(D3g13(n_par,n_r)); D3g13 = g_F_FD(:,:,1,3,0,0,1)
        allocate(g23(n_par,n_r)); g23 = g_F_FD(:,:,2,3,0,0,0)
        allocate(D1g23(n_par,n_r)); D1g23 = g_F_FD(:,:,2,3,1,0,0)
        allocate(D3g23(n_par,n_r)); D3g23 = g_F_FD(:,:,2,3,0,0,1)
        allocate(g33(n_par,n_r)); g33 = g_F_FD(:,:,3,3,0,0,0)
        allocate(D1g33(n_par,n_r)); D1g33 = g_F_FD(:,:,3,3,1,0,0)
        allocate(D2g33(n_par,n_r)); D2g33 = g_F_FD(:,:,3,3,0,1,0)
        allocate(D3g33(n_par,n_r)); D3g33 = g_F_FD(:,:,3,3,0,0,1)
        ! upper metric factors
        allocate(h12(n_par,n_r)); h12 = h_F_FD(:,:,1,2,0,0,0)
        allocate(D3h12(n_par,n_r)); D3h12 = h_F_FD(:,:,1,2,0,0,1)
        allocate(h22(n_par,n_r)); h22 = h_F_FD(:,:,2,2,0,0,0)
        allocate(D3h22(n_par,n_r)); D3h22 = h_F_FD(:,:,2,2,0,0,1)
        allocate(h23(n_par,n_r)); h23 = h_F_FD(:,:,2,3,0,0,0)
        
        ! calculate sigma
        sigma = 1/(mu_0*J*g33) * (g13*(D2g33-D3g23) + g23*(D3g13-D1g33) &
            &+ g33*(D1g23-D2g13))
        
        ! calculate extra1
        extra1 = -D3h12/h22 + D3h22*h12/h22**2
        
        ! calculate extra2
        extra2 = g33/h22 * mu_0*sigma
        
        ! calculate extra3
        do kd = 1,n_r
            extra3(:,kd) = p(kd,1) * ( J(:,kd)**2*mu_0*p(kd,1)/g33(:,kd) + &
                &1/h22(:,kd) * ( &
                &h12(:,kd) * ( D1g33(:,kd)/g33(:,kd) - 2*D1J(:,kd)/J(:,kd) ) + &
                &h22(:,kd) * ( D2g33(:,kd)/g33(:,kd) - 2*D2J(:,kd)/J(:,kd) ) + &
                &h23(:,kd) * ( D3g33(:,kd)/g33(:,kd) - 2*D3J(:,kd)/J(:,kd) ) ) )
        end do
        
        ! calculate rho from input
        call calc_rho(n_par,n_r)
    end subroutine
end module
