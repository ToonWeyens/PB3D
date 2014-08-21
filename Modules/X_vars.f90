!------------------------------------------------------------------------------!
!   Holds variables pertaining to the perturbation quantities
!------------------------------------------------------------------------------!
module X_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, mu_0, max_str_ln, iu
    use output_ops, only: lvl_ud, writo, print_GP_2D
    use str_ops, only: r2strt, i2str
    
    implicit none
    
    private
    public calc_rho, init_X_ops, calc_PV, calc_KV, calc_U, calc_extra, &
        &dealloc_X_vars, &
        &rho, n_x, min_r, max_r, m_X, n_r_X, U_X_0, U_X_1, DU_X_0, DU_X_1, &
        &sigma, extra1, extra2, extra3, PV0, PV1, PV2, KV0, KV1, KV2, X_vec, &
        &X_val, min_m_X, max_m_X

    ! global variables
    real(dp), allocatable :: rho(:,:)                                           ! density
    integer :: n_X                                                              ! toroidal mode number
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    integer, allocatable :: m_X(:)                                              ! vector of poloidal mode numbers
    integer :: n_r_X                                                            ! number of normal points for perturbations
    real(dp) :: min_r, max_r                                                    ! mimimum and maximum radius
    complex(dp), allocatable :: U_X_0(:,:,:), U_X_1(:,:,:)                      ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
    complex(dp), allocatable :: DU_X_0(:,:,:), DU_X_1(:,:,:)                    ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
    real(dp), allocatable :: sigma(:,:)                                         ! parallel current
    real(dp), allocatable :: extra1(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra2(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra3(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    complex(dp), allocatable :: PV0(:,:,:,:)                                    ! ~PV^0 coefficient
    complex(dp), allocatable :: PV1(:,:,:,:)                                    ! ~PV^1 coefficient
    complex(dp), allocatable :: PV2(:,:,:,:)                                    ! ~PV^2 coefficient
    complex(dp), allocatable :: KV0(:,:,:,:)                                    ! ~KV^0 coefficient
    complex(dp), allocatable :: KV1(:,:,:,:)                                    ! ~KV^1 coefficient
    complex(dp), allocatable :: KV2(:,:,:,:)                                    ! ~KV^2 coefficient
    complex(dp), allocatable :: X_vec(:,:,:)                                    ! Eigenvector solution
    complex(dp), allocatable :: X_val(:)                                        ! EIgenvalue solution
    
contains
    ! initialize the variables of this module
    ! ¡¡¡ THIS SHOULD BE  CHANGED TO TAKE INTO ACCOUNT THAT  U_X AND DU_X ARE
    ! HERMITIAN SO ONLY  M*(M+1)/2 ELEMENTS ARE NEEDED INSTEAD  OF M^2. HOWEVER,
    ! TAKE INTO ACCOUNT AS WELL THE MPI STORAGE MATRIX FORMATS !!!
    integer function init_X_ops() result(ierr)
        use eq_vars, only: n_par, n_r
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
        allocate(m_X(n_m_X))
        m_X = [(id, id = min_m_X, max_m_X)]
        if (plot_q) then
            call resonance_plot
        else
            call writo('resonance plot not requested')
        end if
        ierr = check_m_and_n()
        CHCKERR('')
        
        ! U_X_0
        allocate(U_X_0(n_par,n_r,n_m_X))
        
        ! U_X_1
        allocate(U_X_1(n_par,n_r,n_m_X))
        
        ! DU_X_0
        allocate(DU_X_0(n_par,n_r,n_m_X))
        
        ! DU_X_1
        allocate(DU_X_1(n_par,n_r,n_m_X))
        
        ! PV0
        allocate(PV0(n_par,n_r,n_m_X,n_m_X))
        
        ! PV1
        allocate(PV1(n_par,n_r,n_m_X,n_m_X))
        
        ! PV2
        allocate(PV2(n_par,n_r,n_m_X,n_m_X))
        
        ! KV0
        allocate(KV0(n_par,n_r,n_m_X,n_m_X))
        
        ! KV1
        allocate(KV1(n_par,n_r,n_m_X,n_m_X))
        
        ! KV2
        allocate(KV2(n_par,n_r,n_m_X,n_m_X))
    contains
        ! checks whether nq-m is << n is satisfied in some of the plasma
        integer function check_m_and_n() result(ierr)
            use eq_vars, only: q_saf_FD
            
            character(*), parameter :: rout_name = 'check_m_and_n'
            
            ! local variables
            integer :: id                                                       ! counter
            real(dp) :: tol = 0.1                                               ! tolerance to be out of range of q values
            real(dp) :: min_q, max_q                                            ! min. and max. values of q
            
            ! initialize ierr
            ierr = 0
            
            ! set min_q and max_q
            min_q = minval(q_saf_FD(:,0))
            max_q = maxval(q_saf_FD(:,0))
            
            ! for every  poloidal mode number m whether m/n  is inside the range
            ! of q values
            do id = 1, n_m_X
                if (m_X(id)*1.0/n_X .lt.(1-tol)*min_q .or. &
                    &m_X(id)*1.0/n_X .gt.(1+tol)*max_q) then
                    call writo('for m_X = '//trim(i2str(m_X(id)))//', the ratio &
                        &|nq-m|/n is never << 1')
                    ierr = 1
                    err_msg = 'Choose m and n so that (n*q-m)/n << 1'
                    CHCKERR(err_msg)
                end if
            end do
        end function check_m_and_n
    end function init_X_ops
    
    ! plot q-profile with nq-m = 0 indicated if requested
    subroutine resonance_plot
        use num_vars, only: glob_rank, tol_NR
        use utilities, only: calc_zero_NR, calc_interp
        use eq_vars, only: q_saf, n_r
        
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
                ! extrapolate variable from 1
                res = q_saf(1,0) - m_X_for_function*1.0_dp/n_X + &
                    &q_saf(1,1) * (pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from n_r
                res = q_saf(n_r,0) - m_X_for_function*1.0_dp/n_X + &
                    &q_saf(n_r,1) * (pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using calc_interp
                pt_copy = pt
                
                varin(1,:,1,1) = q_saf(:,0) - m_X_for_function*1.0_dp/n_X
                istat = calc_interp(varin,[1,n_r],varout,pt_copy)
                
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
                ! extrapolate variable from 1
                res = q_saf(1,1) + q_saf(1,2) * (pt-0.0_dp)
            else if (pt.gt.1) then                                              ! point requested higher than 1
                ! extrapolate variable from n_r
                res = q_saf(n_r,1) + q_saf(n_r,2) * (pt-1.0_dp)
            else                                                                ! point requested between 0 and 1
                ! interpolate using calc_interp
                pt_copy = pt
                
                varin(1,:,1,1) = q_saf(:,1)
                istat = calc_interp(varin,[1,n_r],varout,pt_copy)
                
                res = varout(1,1,1)
            end if
        end function q_dfun
    end subroutine resonance_plot
    
    ! calculate rho from user input
    ! TO BE IMPLEMENTED. TEMPORARILY SET TO 1 EVERYWHERE !!!
    subroutine calc_rho(n_par,n_r)
        ! input variables
        integer, intent(in) :: n_par, n_r                                       ! dimensions of the size to be allocated for the rho matrix
        
        call writo('calc_rho NOT YET IMPLEMENTED!!!')
        
        ! allocate rho
        allocate(rho(n_par,n_r))
        
        ! TEMPORAL REPLACEMENT !!!
        rho = 1.0E-7_dp
        call writo('TEMPORALLY SETTING rho to '//trim(r2strt(rho(1,1)))//'!!!')
    end subroutine
    
    ! calculate ~PV_(k,m)^i at all n_r values (see [ADD REF] for details)
    subroutine calc_PV
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par, q => q_saf_FD, n_r
        
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
        com_fac = h22/(mu_0*g33)
        
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
                
                !write(*,*) 'PV1'
                !call print_GP_2D('RE PV2(:,:,'//trim(i2str(k))//','//&
                    !&trim(i2str(m))//')','',realpart(transpose(PV1(:,:,k,m))))
                !call print_GP_2D('IM PV2(:,:,'//trim(i2str(k))//','//&
                    !&trim(i2str(m))//')','',imagpart(transpose(PV1(:,:,k,m))))
            end do
        end do
        
        ! deallocate sigma and the extra's
        deallocate(sigma,extra1,extra2,extra3)
    end subroutine calc_PV
    
    ! calculate ~KV_(k,m)^i at all n_r values (see [ADD REF] for details)
    subroutine calc_KV
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par, n_r
        
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
    
    ! calculate  U_m^0, U_m^1  at n_r  values  of the  normal coordinate,  n_par
    ! values  of the  parallel  coordinate and  M values  of  the poloidal  mode
    ! number, with m = size(m_X)
    subroutine calc_U
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: q => q_saf_FD, theta, n_par, p => pres_FD, n_r
        
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
                DU_X_1(:,kd,jd) = n_frac * (D3g13(:,kd)/g33(:,kd) - &
                    &D3g33(:,kd)*g13(:,kd)/g33(:,kd)**2) + &
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
                    &- D3g23(:,kd) )) + iu*n_frac*n_X*U_X_0(:,kd,jd)
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
    end subroutine
    
    ! calculate extra1, extra2 and extra3
    !   extra1 = S*J
    !   extra2 = mu_0*sigma*J*B^2/h^psi,psi
    !   extra3 = 2*p'*kn
    ! with
    !   shear S = - d Theta^alpha/d_theta 1 / J
    !   parallel current mu0 sigma = B . nabla x B / B^2
    !   normal curvature kn = nabla psi / h^psi,psi * nabla (mu0 p + B^2/2)
    subroutine calc_extra
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par, p => pres_FD, n_r
        
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
        
        ! allocatee sigma and extra terms. Later they are deallocated in calc_PV
        allocate(sigma(n_par,n_r))
        allocate(extra1(n_par,n_r))
        allocate(extra2(n_par,n_r))
        allocate(extra3(n_par,n_r))
        
        ! calculate sigma
        sigma = 1/(mu_0*J*g33) * (g13*(D2g33-D3g23) + g23*(D3g13-D1g33) &
            &+ g33*(D1g23-D2g13))
        
        ! calculate extra1
        extra1 = -D3h12/h22 + D3h22*h12/h22**2
        
        ! calculate extra2
        extra2 = g33/h22 * mu_0*(sigma*J)
        
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
    
    ! deallocates  perturbation quantities that  are not used anymore  after the
    ! equilibrium phase
    subroutine dealloc_X_vars
        use num_vars, only: group_rank
        
        deallocate(rho)
        deallocate(m_X)
        deallocate(U_X_0,U_X_1,DU_X_0,DU_X_1)
        deallocate(PV0,PV1,PV2,KV0,KV1,KV2)
        if (group_rank.eq.0) then
            deallocate(X_vec,X_val)
        end if
    end subroutine dealloc_X_vars
end module
