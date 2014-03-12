module metric_ops
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud
    use var_ops, only: r2str, i2str
    use VMEC_vars, only: n_r
    use grid_vars, only: n_theta, n_zeta, R, Z
    
    implicit none
    private
    public metric_C2V, metric_V, metric_C, &
        &C2V_up, C2V_dn, jac_V, h_V, g_V

    ! upper (h) and lower (g) metric factors
    ! (index 1,2,3: theta, zeta, r, index 4,5: 3x3 matrix)
    real(dp), allocatable :: g_C(:,:,:,:,:), h_C(:,:,:,:,:)                     ! in the C(ylindrical) coordinate system
    real(dp), allocatable :: g_V(:,:,:,:,:), h_V(:,:,:,:,:)                     ! in the V(MEC) coordinate system
    ! upper and lower transformation matrices
    ! (index 1,2,3: theta, zeta, r, index 4,5: 3x3 matrix)
    real(dp), allocatable :: C2V_up(:,:,:,:,:), C2V_dn(:,:,:,:,:)
    real(dp), allocatable :: jac_V(:,:,:)                                       ! jacobian of V(MEC) coordinate system

contains

    ! calculate the trivial metric elements in the C(ylindrical) coordinate system
    subroutine metric_C
        use VMEC_vars, only: n_r
        use grid_vars, only: R, n_theta, n_zeta
        
        ! local variables
        integer :: kd
        
        ! checking for previous allocation
        if (allocated(g_C)) deallocate(g_C)
        if (allocated(h_C)) deallocate(h_C)
        
        allocate(g_C(n_theta,n_zeta,n_r,3,3)); g_C = 0.0_dp
        allocate(h_C(n_theta,n_zeta,n_r,3,3)); h_C = 0.0_dp
        
        do kd = 1,n_r
            g_C(:,:,:,1,1) = 1.0_dp
            g_C(:,:,:,2,2) = R(:,:,:,1)**2
            g_C(:,:,:,3,3) = 1.0_dp
            h_C(:,:,:,1,1) = 1.0_dp
            h_C(:,:,:,2,2) = 1.0_dp/R(:,:,:,1)**2
            h_C(:,:,:,3,3) = 1.0_dp
        end do
    end subroutine

    ! Calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system from the VMEC file.
    subroutine metric_C2V
        ! local variables
        real(dp) :: cf(n_theta,n_zeta)                                          ! common factor used in calculating the elements of the matrix
        integer :: kd
        
        ! deallocate if allocated
        if (allocated(C2V_dn)) deallocate(C2V_dn)
        if (allocated(C2V_up)) deallocate(C2V_up)
        if (allocated(jac_V)) deallocate(jac_V)
        ! reallocate
        allocate(C2V_dn(n_theta,n_zeta,n_r,3,3))
        allocate(C2V_up(n_theta,n_zeta,n_r,3,3))
        allocate(jac_V(n_theta,n_zeta,n_r))
        
        ! calculate transformation matrix for contravariant coordinates
        do kd = 1,n_r
            jac_V(:,:,kd) = R(:,:,kd,1)*(R(:,:,kd,2)*Z(:,:,kd,3)&
                &-R(:,:,kd,3)*Z(:,:,kd,2))
            cf = R(:,:,kd,1)/jac_V(:,:,kd)
            C2V_up(:,:,kd,1,1) = cf*Z(:,:,kd,3)
            C2V_up(:,:,kd,1,2) = cf*(R(:,:,kd,4)*Z(:,:,kd,3)&
                &-R(:,:,kd,3)*Z(:,:,kd,4))
            C2V_up(:,:,kd,1,3) = -cf*R(:,:,kd,3)
            C2V_up(:,:,kd,2,1) = -cf*Z(:,:,kd,2)
            C2V_up(:,:,kd,2,2) = cf*(R(:,:,kd,2)*Z(:,:,kd,4)&
                &-R(:,:,kd,4)*Z(:,:,kd,2))
            C2V_up(:,:,kd,2,3) = cf*R(:,:,kd,2)
            C2V_up(:,:,kd,3,1) = 0.0_dp
            C2V_up(:,:,kd,3,2) = cf*(R(:,:,kd,3)*Z(:,:,kd,2)&
                &-R(:,:,kd,2)*Z(:,:,kd,3))
            C2V_up(:,:,kd,3,3) = 0.0_dp
        end do
        
        ! calculate transformation matrix for covariant coordinates
        do kd = 1,n_r
            C2V_dn(:,:,kd,1,1) = R(:,:,kd,2)
            C2V_dn(:,:,kd,1,2) = 0.0_dp
            C2V_dn(:,:,kd,1,3) = Z(:,:,kd,2)
            C2V_dn(:,:,kd,2,1) = R(:,:,kd,3)
            C2V_dn(:,:,kd,2,2) = 0.0_dp
            C2V_dn(:,:,kd,2,3) = Z(:,:,kd,3)
            C2V_dn(:,:,kd,3,1) = R(:,:,kd,4)
            C2V_dn(:,:,kd,3,2) = -1
            C2V_dn(:,:,kd,3,3) = Z(:,:,kd,4)
        end do
    end subroutine

    ! calculate the  metric coefficients in  the V(MEC) coordinate  system using
    ! the metric  coefficients in  the C(ylindrical)  coordinate system  and the
    ! trans- formation matrices
    subroutine metric_V
        ! local variables
        integer :: id, jd, kd
        
        ! deallocate if allocated
        if (allocated(g_V)) deallocate(g_V)
        if (allocated(h_V)) deallocate(h_V)
        ! reallocate
        allocate(g_V(n_theta,n_zeta,n_r,3,3))
        allocate(h_V(n_theta,n_zeta,n_r,3,3))
            
        do kd = 1,n_r
            do jd = 1,n_zeta
                do id = 1,n_theta
                    g_V(id,jd,kd,:,:) = metric_calc(g_C(id,jd,kd,:,:),&
                        &C2V_dn(id,jd,kd,:,:))
                    h_V(id,jd,kd,:,:) = metric_calc(h_C(id,jd,kd,:,:),&
                        &C2V_up(id,jd,kd,:,:))
                end do
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
