!------------------------------------------------------------------------------!
!   Calculate the matrix  elements due to the Plasma of  the system of equations
!   that has to be solved as described in [ADD REF]
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, mu_0, max_str_ln
    use output_ops, only: lvl_ud, writo, print_ar_2

    implicit none
    private
    public prepare_matrix_X, solve_EV_system

contains
    ! initialize the variables of this module
    ! ¡¡¡ THIS SHOULD BE  CHANGED TO TAKE INTO ACCOUNT THAT  U_X AND DU_X ARE
    ! HERMITIAN SO ONLY  M*(M+1)/2 ELEMENTS ARE NEEDED INSTEAD  OF M^2. HOWEVER,
    ! TAKE INTO ACCOUNT AS WELL THE MPI STORAGE MATRIX FORMATS !!!
    subroutine init_X_ops(n_m_X)
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        use X_vars, only: U_X_0, U_X_1, DU_X_0, DU_X_1, sigma, extra1, extra2, &
            &extra3, PV0, PV1, PV2, KV0, KV1, KV2
        
        ! input / output
        integer :: n_m_X
        
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
    end subroutine init_X_ops
    
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    subroutine prepare_matrix_X
        use X_vars, only: n_X, m_X
        
        ! initialize the variables
        call writo('Initalizing variables...')
        call lvl_ud(1)
        call init_X_ops(size(m_X))
        call lvl_ud(-1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        call calc_U(n_X,m_X)
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
        call calc_PV(n_X,m_X)
        call lvl_ud(-1)
        
        ! Calculate KV0, KV1  and KV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV0, KV1 and KV2...')
        call lvl_ud(1)
        call calc_KV(size(m_X))
        call lvl_ud(-1)
    end subroutine
    
    ! calculate ~PV_(k,m)^i at all n_r values (see [ADD REF] for details)
    subroutine calc_PV(n_X,m_X)
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par, q => q_saf_FD
        use VMEC_vars, only: n_r
        use X_vars, only: DU_X_0, DU_X_1, sigma, extra1, extra2, extra3, PV0, PV1, PV2
        
        ! input / output
        integer, intent(in) :: n_X                                              ! tor. mode nr n
        integer, intent(in) :: m_X(:)                                           ! vector containing pol. mode nrs m
        
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
    subroutine calc_KV(n_m_X)
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        use X_vars, only: rho, U_X_0, U_X_1, KV0, KV1, KV2
        
        ! input / output
        integer, intent(in) :: n_m_X                                            ! number of poloidal modes m
        
        ! local variables
        integer :: m, k                                                         ! counters
        real(dp), allocatable :: com_fac(:,:)                                   ! common factor |nabla psi|^2/(mu_0*J^2*B^2)
        
        ! submatrices
        ! jacobian
        real(dp), allocatable :: J(:,:)                                         ! jac
        ! lower metric factors
        real(dp), allocatable :: g33(:,:)                                       ! h^alpha,psi
        ! upper metric factors
        real(dp), allocatable :: h22(:,:)                                       ! h^alpha,psi
        
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
    subroutine solve_EV_system(ierr)
        use X_vars, only: n_r_X, m_X, n_X
        use VMEC_vars, only: n_r
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use slepc_ops, only: solve_EV_system_slepc
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        integer, intent(inout) :: ierr                                          ! error
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! slepc solver for EV problem
                call solve_EV_system_slepc(m_X,n_r_X,n_X,n_r-1._dp,ierr)
                CHCKERR('')
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end subroutine
    
    ! calculate  U_m^0, U_m^1  at n_r  values  of the  normal coordinate,  n_par
    ! values  of the  parallel  coordinate and  M values  of  the poloidal  mode
    ! number, with m = size(m_X)
    subroutine calc_U(n_X,m_X)
        use metric_ops, only: g_F_FD, h_F_FD, jac_F_FD
        use eq_vars, only: q => q_saf_FD, theta, n_par, p => pres_FD
        use VMEC_vars, only: n_r
        use X_vars, only: U_X_0, U_X_1, DU_X_0, DU_X_1
        
        ! input / output
        integer, intent(in) :: n_X                                              ! tor. mode nr n
        integer, intent(in) :: m_X(:)                                           ! vector containing pol. mode nrs m
        
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
