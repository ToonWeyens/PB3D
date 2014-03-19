!-------------------------------------------------------
!   Variables, subroutines and  functions that have to do with the metric elements
!-------------------------------------------------------
module metric_ops 
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud
    use str_ops, only: r2str, i2str
    use VMEC_vars, only: n_r
    use eq_vars, only: n_par, R, Z
    
    implicit none
    private
    public metric_C2V, metric_V, metric_C, metric_V2F, metric_F, &
        &C2V_up, C2V_dn, jac_V, h_V, g_V, V2F_up, V2F_dn, jac_F, h_F, g_F

    ! upper (h) and lower (g) metric factors
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: g_C(:,:,:,:), h_C(:,:,:,:)                         ! in the C(ylindrical) coordinate system
    real(dp), allocatable :: g_V(:,:,:,:), h_V(:,:,:,:)                         ! in the V(MEC) coordinate system
    real(dp), allocatable :: g_F(:,:,:,:), h_F(:,:,:,:)                         ! in the F(lux) coordinate system
    ! upper and lower transformation matrices
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: C2V_up(:,:,:,:), C2V_dn(:,:,:,:)                   ! C(ylindrical) to V(MEC)
    real(dp), allocatable :: V2F_up(:,:,:,:), V2F_dn(:,:,:,:)                   ! V(MEC) to F(lux)
    real(dp), allocatable :: jac_V(:,:)                                         ! jacobian of V(MEC) coordinate system
    real(dp), allocatable :: jac_F(:,:)                                         ! jacobian of F(lux) coordinate system

contains

    ! calculate the trivial metric elements in the C(ylindrical) coordinate system
    subroutine metric_C
        ! local variables
        integer :: kd
        
        ! checking for previous allocation
        if (allocated(g_C)) deallocate(g_C)
        if (allocated(h_C)) deallocate(h_C)
        
        allocate(g_C(3,3,n_par,n_r)); g_C = 0.0_dp
        allocate(h_C(3,3,n_par,n_r)); h_C = 0.0_dp
        
        do kd = 1,n_r
            g_C(1,1,:,:) = 1.0_dp
            g_C(2,2,:,:) = R(:,:,1)**2
            g_C(3,3,:,:) = 1.0_dp
            h_C(1,1,:,:) = 1.0_dp
            h_C(2,2,:,:) = 1.0_dp/R(:,:,1)**2
            h_C(3,3,:,:) = 1.0_dp
        end do
    end subroutine

    ! Calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system from the VMEC file.
    subroutine metric_C2V
        ! local variables
        real(dp) :: cf(n_par)                                                  ! common factor used in calculating the elements of the matrix
        integer :: kd
        
        ! deallocate if allocated
        if (allocated(C2V_dn)) deallocate(C2V_dn)
        if (allocated(C2V_up)) deallocate(C2V_up)
        if (allocated(jac_V)) deallocate(jac_V)
        ! reallocate
        allocate(C2V_dn(3,3,n_par,n_r)); C2V_dn = 0.0_dp
        allocate(C2V_up(3,3,n_par,n_r)); C2V_up = 0.0_dp
        allocate(jac_V(n_par,n_r)); jac_V = 0.0_dp
        
        ! calculate transformation matrix for contravariant coordinates
        do kd = 1,n_r
            jac_V(:,kd) = R(:,kd,1)*(R(:,kd,2)*Z(:,kd,3)-R(:,kd,3)*Z(:,kd,2))
            cf = R(:,kd,1)/jac_V(:,kd)
            C2V_up(1,1,:,kd) = cf*Z(:,kd,3)
            C2V_up(1,2,:,kd) = cf*(R(:,kd,4)*Z(:,kd,3)-R(:,kd,3)*Z(:,kd,4))
            C2V_up(1,3,:,kd) = -cf*R(:,kd,3)
            C2V_up(2,1,:,kd) = -cf*Z(:,kd,2)
            C2V_up(2,2,:,kd) = cf*(R(:,kd,2)*Z(:,kd,4)-R(:,kd,4)*Z(:,kd,2))
            C2V_up(2,3,:,kd) = cf*R(:,kd,2)
            C2V_up(3,1,:,kd) = 0.0_dp
            C2V_up(3,2,:,kd) = cf*(R(:,kd,3)*Z(:,kd,2)-R(:,kd,2)*Z(:,kd,3))
            C2V_up(3,3,:,kd) = 0.0_dp
        end do
        
        ! calculate transformation matrix for covariant coordinates
        do kd = 1,n_r
            C2V_dn(1,1,:,kd) = R(:,kd,2)
            C2V_dn(1,2,:,kd) = 0.0_dp
            C2V_dn(1,3,:,kd) = Z(:,kd,2)
            C2V_dn(2,1,:,kd) = R(:,kd,3)
            C2V_dn(2,2,:,kd) = 0.0_dp
            C2V_dn(2,3,:,kd) = Z(:,kd,3)
            C2V_dn(3,1,:,kd) = R(:,kd,4)
            C2V_dn(3,2,:,kd) = -1
            C2V_dn(3,3,:,kd) = Z(:,kd,4)
        end do
    end subroutine

    ! Calculate  the transformation  matrix between  the V(mec)  and the  F(lux)
    ! coordinate system
    subroutine metric_V2F
        use eq_vars, only: q_saf, theta, n_par, lam, flux_p
        use num_vars, only: pi
        
        ! local variables
        integer :: id, kd
        real(dp) :: theta_s                                                     ! = theta + lambda
        
        ! deallocate if allocated
        if (allocated(V2F_dn)) deallocate(V2F_dn)
        if (allocated(V2F_up)) deallocate(V2F_up)
        if (allocated(jac_F)) deallocate(jac_F)
        ! reallocate
        allocate(V2F_dn(3,3,n_par,n_r)); V2F_dn = 0.0_dp
        allocate(V2F_up(3,3,n_par,n_r)); V2F_up = 0.0_dp
        allocate(jac_F(n_par,n_r)); jac_F = 0.0_dp
        
        do kd = 1,n_r
            do id = 1, n_par
                theta_s = theta(id,kd)+lam(id,kd,1)
                V2F_up(1,1,id,kd) = -q_saf(kd,2)*theta_s &
                    &- q_saf(kd,1)*lam(id,kd,2)
                V2F_up(1,2,id,kd) = -q_saf(kd,1)*(1+lam(id,kd,3))
                V2F_up(1,3,id,kd) = 1 - q_saf(kd,1)*lam(id,kd,4)
                V2F_up(2,1,id,kd) = flux_p(kd,2)/(2*pi)
                V2F_up(2,2,id,kd) = 0.0_dp
                V2F_up(2,3,id,kd) = 0.0_dp
                V2F_up(3,1,id,kd) = lam(id,kd,2)
                V2F_up(3,2,id,kd) = 1 + lam(id,kd,3)
                V2F_up(3,3,id,kd) = lam(id,kd,4)
                jac_F(id,kd) = (flux_p(kd,2)/(2*pi)*(1+lam(id,kd,3)))**(-1)*&
                    &jac_V(id,kd)
            end do
        end do
        
        do kd = 1,n_r
            do id = 1, n_par
                theta_s = theta(id,kd)+lam(id,kd,1)
                V2F_dn(1,1,id,kd) = 0.0_dp
                V2F_dn(1,2,id,kd) = -lam(id,kd,4)/(1+lam(id,kd,3))
                V2F_dn(1,3,id,kd) = 1.0_dp
                V2F_dn(2,1,id,kd) = 2*pi/flux_p(kd,2)
                V2F_dn(2,2,id,kd) = -2*pi/flux_p(kd,2) * (lam(id,kd,2)+&
                    &lam(id,kd,4)*q_saf(kd,2)*theta_s)/(1+lam(id,kd,3))
                V2F_dn(2,3,id,kd) = 2*pi/flux_p(kd,2)*q_saf(kd,2)*theta_s
                V2F_dn(3,1,id,kd) = 0.0_dp
                V2F_dn(3,2,id,kd) = (1-q_saf(kd,1)*lam(id,kd,4))/&
                    &(1+lam(id,kd,3))
                V2F_dn(3,3,id,kd) = q_saf(kd,1)
            end do
        end do
    end subroutine

    ! calculate the  metric coefficients in  the V(MEC) coordinate  system using
    ! the metric  coefficients in  the C(ylindrical)  coordinate system  and the
    ! trans- formation matrices
    subroutine metric_V
        ! local variables
        integer :: id, kd
        
        ! deallocate if allocated
        if (allocated(g_V)) deallocate(g_V)
        if (allocated(h_V)) deallocate(h_V)
        ! reallocate
        allocate(g_V(3,3,n_par,n_r))
        allocate(h_V(3,3,n_par,n_r))
            
        do kd = 1,n_r
            do id = 1,n_par
                g_V(:,:,id,kd) = &
                    &metric_calc(g_C(:,:,id,kd), C2V_dn(:,:,id,kd))
                h_V(:,:,id,kd) = &
                    &metric_calc(h_C(:,:,id,kd), C2V_up(:,:,id,kd))
            end do
        end do
    end subroutine

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the metric  coefficients in  the V(MEC) coordinate  system and  the trans-
    ! formation matrices
    subroutine metric_F
        ! local variables
        integer :: id, kd
        
        ! deallocate if allocated
        if (allocated(g_F)) deallocate(g_F)
        if (allocated(h_F)) deallocate(h_F)
        ! reallocate
        allocate(g_F(3,3,n_par,n_r))
        allocate(h_F(3,3,n_par,n_r))
            
        do kd = 1,n_r
            do id = 1,n_par
                g_F(:,:,id,kd) = &
                    &metric_calc(g_V(:,:,id,kd), V2F_dn(:,:,id,kd))
                h_F(:,:,id,kd) = &
                    &metric_calc(h_V(:,:,id,kd), V2F_up(:,:,id,kd))
            end do
        end do
    end subroutine
    
    ! calculate the metric coefficients from the transformation matrices  in the 
    ! coordinate system B using 
    ! (lower) g_B = TC2V_dn*g_A*TC2V_dn^T
    ! (upper) h_B = TC2V_up*h_A*T2CV_up^T
    ! where g_A and h_A are the metric coefficients of the lower (covariant) and
    ! upper (contravariant) unit vectors in the coordinate system A
    function metric_calc(gh_A, TC2V)
        use var_ops, only: mat_mult
        
        ! input / output
        real(dp), intent(in) :: gh_A(3,3)
        real(dp), intent(in) :: TC2V(3,3)
        real(dp) :: metric_calc(3,3)
        
        ! local variables
        real(dp) :: temp_mat(3,3)
        
        temp_mat = mat_mult(gh_A,transpose(TC2V))
        metric_calc = mat_mult(TC2V,temp_mat)
    end function metric_calc

end module metric_ops
