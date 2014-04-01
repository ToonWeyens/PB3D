!------------------------------------------------------------------------------!
!   Variables, subroutines and  functions that have to do with the metric      !
!   elements                                                                   !
!------------------------------------------------------------------------------!
module metric_ops 
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud, write_out
    use str_ops, only: r2str, i2str
    
    implicit none
    private
    public metric_C2V, metric_V, metric_C, metric_V2F, metric_F, &
        &C2V_up, C2V_dn, C2V_up_H, C2V_dn_H, jac_V, h_V, g_V, h_V_H, g_V_H, &
        &V2F_up, V2F_dn, V2F_up_H, V2F_dn_H, jac_F, jac_F_H, &
        &h_F, g_F, h_F_H, g_F_H

    ! upper (h) and lower (g) metric factors
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: g_C(:,:,:,:), h_C(:,:,:,:)                         ! in the C(ylindrical) coordinate system (FM)
    real(dp), allocatable :: g_C_H(:,:,:,:), h_C_H(:,:,:,:)                     ! in the C(ylindrical) coordinate system (HM)
    real(dp), allocatable :: g_V(:,:,:,:), h_V(:,:,:,:)                         ! in the V(MEC) coordinate system (FM)
    real(dp), allocatable :: g_V_H(:,:,:,:), h_V_H(:,:,:,:)                     ! in the V(MEC) coordinate system (HM)
    real(dp), allocatable :: g_F(:,:,:,:), h_F(:,:,:,:)                         ! in the F(lux) coordinate system (FM)
    real(dp), allocatable :: g_F_H(:,:,:,:), h_F_H(:,:,:,:)                     ! in the F(lux) coordinate system (HM)
    ! upper and lower transformation matrices
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: C2V_up(:,:,:,:), C2V_dn(:,:,:,:)                   ! C(ylindrical) to V(MEC) (FM)
    real(dp), allocatable :: C2V_up_H(:,:,:,:), C2V_dn_H(:,:,:,:)               ! C(ylindrical) to V(MEC) (HM)
    real(dp), allocatable :: V2F_up(:,:,:,:), V2F_dn(:,:,:,:)                   ! V(MEC) to F(lux) (FM)
    real(dp), allocatable :: V2F_up_H(:,:,:,:), V2F_dn_H(:,:,:,:)               ! V(MEC) to F(lux) (HM)
    real(dp), allocatable :: jac_V(:,:)                                         ! jacobian of V(MEC) coordinate system (FM)
    real(dp), allocatable :: jac_V_H(:,:)                                       ! jacobian, of V(MEC) coordinate system (HM)
    real(dp), allocatable :: jac_F(:,:)                                         ! jacobian of F(lux) coordinate system (FM)
    real(dp), allocatable :: jac_F_H(:,:)                                       ! jacobian of F(lux) coordinate system (HM)

contains

    ! calculate the trivial metric elements in the C(ylindrical) coordinate system
    subroutine metric_C
        use eq_vars, only: f2h, &
            &n_par, R
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: jd
        real(dp) :: R_H(n_r)                                                    ! HM version of R
        
        ! checking for previous allocation
        if (allocated(g_C)) deallocate(g_C)
        if (allocated(h_C)) deallocate(h_C)
        if (allocated(g_C_H)) deallocate(g_C_H)
        if (allocated(h_C_H)) deallocate(h_C_H)
        
        allocate(g_C(3,3,n_par,n_r)); g_C = 0.0_dp
        allocate(h_C(3,3,n_par,n_r)); h_C = 0.0_dp
        allocate(g_C_H(3,3,n_par,n_r)); g_C_H = 0.0_dp
        allocate(h_C_H(3,3,n_par,n_r)); h_C_H = 0.0_dp
        
        do jd = 1,n_par
            ! FM
            g_C(1,1,jd,:) = 1.0_dp
            g_C(2,2,jd,:) = R(jd,:,1)**2
            g_C(3,3,jd,:) = 1.0_dp
            h_C(1,1,jd,:) = 1.0_dp
            h_C(2,2,jd,:) = 1.0_dp/R(jd,:,1)**2
            h_C(3,3,jd,:) = 1.0_dp
            
            ! HM
            R_H = f2h(R(jd,:,1))
            g_C_H(1,1,jd,:) = 1.0_dp
            g_C_H(2,2,jd,:) = R_H**2
            g_C_H(3,3,jd,:) = 1.0_dp
            h_C_H(1,1,jd,:) = 1.0_dp
            h_C_H(2,2,jd,:) = 1.0_dp/R_H**2
            h_C_H(3,3,jd,:) = 1.0_dp
        end do
    end subroutine

    ! calculate the  metric coefficients in  the V(MEC) coordinate  system using
    ! the metric  coefficients in  the C(ylindrical)  coordinate system  and the
    ! trans- formation matrices
    subroutine metric_V
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: id, kd
        
        ! deallocate if allocated
        if (allocated(g_V)) deallocate(g_V)
        if (allocated(h_V)) deallocate(h_V)
        if (allocated(g_V_H)) deallocate(g_V_H)
        if (allocated(h_V_H)) deallocate(h_V_H)
        ! reallocate
        allocate(g_V(3,3,n_par,n_r))
        allocate(h_V(3,3,n_par,n_r))
        allocate(g_V_H(3,3,n_par,n_r))
        allocate(h_V_H(3,3,n_par,n_r))
            
        do kd = 1,n_r
            do id = 1,n_par
                g_V(:,:,id,kd) = &
                    &metric_calc(g_C(:,:,id,kd), C2V_dn(:,:,id,kd))
                h_V(:,:,id,kd) = &
                    &metric_calc(h_C(:,:,id,kd), C2V_up(:,:,id,kd))
                g_V_H(:,:,id,kd) = &
                    &metric_calc(g_C_H(:,:,id,kd), C2V_dn_H(:,:,id,kd))
                h_V_H(:,:,id,kd) = &
                    &metric_calc(h_C_H(:,:,id,kd), C2V_up_H(:,:,id,kd))
            end do
        end do
    end subroutine

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the metric  coefficients in  the V(MEC) coordinate  system and  the trans-
    ! formation matrices
    subroutine metric_F
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: id, kd
        
        ! deallocate if allocated
        if (allocated(g_F)) deallocate(g_F)
        if (allocated(h_F)) deallocate(h_F)
        if (allocated(g_F_H)) deallocate(g_F_H)
        if (allocated(h_F_H)) deallocate(h_F_H)
        ! reallocate
        allocate(g_F(3,3,n_par,n_r))
        allocate(h_F(3,3,n_par,n_r))
        allocate(g_F_H(3,3,n_par,n_r))
        allocate(h_F_H(3,3,n_par,n_r))
            
        do kd = 1,n_r
            do id = 1,n_par
                g_F(:,:,id,kd) = &
                    &metric_calc(g_V(:,:,id,kd), V2F_dn(:,:,id,kd))
                h_F(:,:,id,kd) = &
                    &metric_calc(h_V(:,:,id,kd), V2F_up(:,:,id,kd))
                g_F_H(:,:,id,kd) = &
                    &metric_calc(g_V_H(:,:,id,kd), V2F_dn_H(:,:,id,kd))
                h_F_H(:,:,id,kd) = &
                    &metric_calc(h_V_H(:,:,id,kd), V2F_up_H(:,:,id,kd))
            end do
        end do
    end subroutine
    
    ! Calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system from the VMEC file.
    subroutine metric_C2V
        use VMEC_vars, only: n_r, jac_V_H_c, jac_V_H_s, mpol, ntor
        use eq_vars, only: h2f, f2h, calc_norm_deriv, &
            &theta_H, zeta_H, n_par, R, Z
        use fourier_ops, only: f2r, mesh_cs
        
        ! local variables
        real(dp) :: cf(n_r)                                                     ! common factor used in calculating the elements of the matrix
        integer :: id, jd, kd
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        real(dp) :: jac_V_H_alt(n_par, n_r)                                     ! Jacobian from VMEC (HM)
        real(dp) :: jac_V_H_dff(n_par, n_r)
        real(dp) :: R_H(n_r,4), Z_H(n_r,4)                                      ! HM versions of R, Z
        
        ! deallocate if allocated
        if (allocated(C2V_dn)) deallocate(C2V_dn)
        if (allocated(C2V_up)) deallocate(C2V_up)
        if (allocated(C2V_dn_H)) deallocate(C2V_dn_H)
        if (allocated(C2V_up_H)) deallocate(C2V_up_H)
        if (allocated(jac_V)) deallocate(jac_V)
        if (allocated(jac_V_H)) deallocate(jac_V_H)
        ! reallocate
        allocate(C2V_dn(3,3,n_par,n_r)); C2V_dn = 0.0_dp
        allocate(C2V_up(3,3,n_par,n_r)); C2V_up = 0.0_dp
        allocate(C2V_dn_H(3,3,n_par,n_r)); C2V_dn_H = 0.0_dp
        allocate(C2V_up_H(3,3,n_par,n_r)); C2V_up_H = 0.0_dp
        allocate(jac_V(n_par,n_r)); jac_V = 0.0_dp
        allocate(jac_V_H(n_par,n_r)); jac_V_H = 0.0_dp
        
        ! transform the VMEC jacobian to real space, using FM variants
        par: do jd = 1, n_par                                                   ! parallel: along the magnetic field line
            ! import jac_V on HM from VMEC
            perp: do kd = 1, n_r                                                ! perpendicular: normal to the flux surfaces
                ! calculate the (co)sines at the current mesh points
                cs = mesh_cs(mpol,ntor,theta_H(jd,kd),zeta_H(jd,kd))
                jac_V_H_alt(jd,kd) = &
                    &-f2r(jac_V_H_c(:,:,kd),jac_V_H_s(:,:,kd),cs,mpol,ntor,1)
            end do perp
        end do par
        
        ! calculate transformation matrix for C -> V (HM), making use of of R_H,
        ! Z_H, and derivatives
        R_H = 0.0_dp
        Z_H = 0.0_dp
        do jd = 1,n_par
            ! calculate R_H and Z_H
            R_H(:,1) = f2h(R(jd,:,1))
            R_H(:,2) = calc_norm_deriv(R(jd,:,1),n_r-1._dp,.true.,.false.)
            R_H(:,3) = f2h(R(jd,:,3))
            R_H(:,4) = f2h(R(jd,:,4))
            Z_H(:,1) = f2h(Z(jd,:,1))
            Z_H(:,2) = calc_norm_deriv(Z(jd,:,1),n_r-1._dp,.true.,.false.)
            Z_H(:,3) = f2h(Z(jd,:,3))
            Z_H(:,4) = f2h(Z(jd,:,4))
            
            !call write_out(2,n_r/10,transpose(reshape([(1.0_dp*kd/(n_r-1),&
                !&kd=1,n_r/10),R(jd,2:n_r/10+1,1)],[n_r/10,2])),&
                !&'R_H at id = '//trim(i2str(id)))
            !call write_out(2,n_r/10,transpose(reshape([(1.0_dp*kd/(n_r-1),&
                !&kd=1,n_r/10),R_H(2:n_r/10+1,2)],[n_r/10,2])),&
                !&'Rr at id = '//trim(i2str(jd)))
            
            jac_V_H(jd,:) = R_H(:,1)*(R_H(:,2)*Z_H(:,3)-R_H(:,3)*Z_H(:,2))
            
            ! upper metric factors
            cf = R_H(:,1)/jac_V_H(jd,:); cf(1) = 0.0_dp
            C2V_up_H(1,1,jd,:) = cf*Z_H(:,3)
            C2V_up_H(1,2,jd,:) = cf*(R_H(:,4)*Z_H(:,3)-R_H(:,3)*Z_H(:,4))
            C2V_up_H(1,3,jd,:) = -cf*R_H(:,3)
            C2V_up_H(2,1,jd,:) = -cf*Z_H(:,2)
            C2V_up_H(2,2,jd,:) = cf*(R_H(:,2)*Z_H(:,4)-R_H(:,4)*Z_H(:,2))
            C2V_up_H(2,3,jd,:) = cf*R_H(:,2)
            C2V_up_H(3,1,jd,:) = 0.0_dp
            C2V_up_H(3,2,jd,:) = -1.0_dp
            C2V_up_H(3,3,jd,:) = 0.0_dp
            
            ! lower metric factors
            C2V_dn_H(1,1,jd,:) = R_H(:,2)
            C2V_dn_H(1,2,jd,:) = 0.0_dp
            C2V_dn_H(1,3,jd,:) = Z_H(:,2)
            C2V_dn_H(2,1,jd,:) = R_H(:,3)
            C2V_dn_H(2,2,jd,:) = 0.0_dp
            C2V_dn_H(2,3,jd,:) = Z_H(:,3)
            C2V_dn_H(3,1,jd,:) = R_H(:,4)
            C2V_dn_H(3,2,jd,:) = -1.0_dp
            C2V_dn_H(3,3,jd,:) = Z_H(:,4)
            
            do kd = 1, n_r
                jac_V_H_dff(jd,kd) = 2*(jac_V_H_alt(jd,kd)-jac_V_H(jd,kd))/&
                    &(jac_V_H_alt(jd,kd)+jac_V_H(jd,kd))
            end do
            
            ! convert to full mesh
            ! jacobian
            jac_V(jd,:) = h2f(jac_V_H(jd,:))
            
            ! transformation matrices
            do id = 1,3
                do kd = 1,3
                    C2V_up(kd,id,jd,:) = h2f(C2V_up_H(kd,id,jd,:))
                    C2V_dn(kd,id,jd,:) = h2f(C2V_dn_H(kd,id,jd,:))
                end do
            end do
        end do
        !write(*,*) 'Jac_V_H = '
        !call print_ar_2(jac_V_H)
        !write(*,*) 'Jac_V_H_alt = '
        !call print_ar_2(jac_V_H_alt)
        !write(*,*) 'Jac_V_H_dff = '
        !call print_ar_2(jac_V_H_dff)
        !write(*,*) 'max diff = ', maxval(abs(jac_V_H_dff))
        !do id = 1,n_par
            !call write_out(2,n_r-1,transpose(reshape([((jd-0.5_dp)/(n_r-1),&
                !&jd=1,n_r-1),jac_V_H(id,2:n_r)],[n_r-1,2])),&
                !&'jac at id = '//trim(i2str(id)))
        !end do
        !read(*,*)
    end subroutine

    ! Calculate  the transformation  matrix between  the V(mec)  and the  F(lux)
    ! coordinate system
    subroutine metric_V2F
        use eq_vars, only: h2f, &
            &q_saf_H, theta_H, n_par, lam_H, flux_p_H
        use num_vars, only: pi
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: id, jd, kd
        real(dp) :: theta_s                                                     ! = theta + lambda
        
        ! deallocate if allocated
        if (allocated(V2F_dn)) deallocate(V2F_dn)
        if (allocated(V2F_up)) deallocate(V2F_up)
        if (allocated(V2F_dn_H)) deallocate(V2F_dn_H)
        if (allocated(V2F_up_H)) deallocate(V2F_up_H)
        if (allocated(jac_F)) deallocate(jac_F)
        if (allocated(jac_F_H)) deallocate(jac_F_H)
        ! reallocate
        allocate(V2F_dn(3,3,n_par,n_r)); V2F_dn = 0.0_dp
        allocate(V2F_up(3,3,n_par,n_r)); V2F_up = 0.0_dp
        allocate(V2F_dn_H(3,3,n_par,n_r)); V2F_dn_H = 0.0_dp
        allocate(V2F_up_H(3,3,n_par,n_r)); V2F_up_H = 0.0_dp
        allocate(jac_F(n_par,n_r)); jac_F = 0.0_dp
        allocate(jac_F_H(n_par,n_r)); jac_F_H = 0.0_dp
        
        ! calculate transformation matrix for V -> F (HM)
        do kd = 1,n_r
            do id = 1, n_par
                theta_s = theta_H(id,kd)+lam_H(id,kd,1)
                V2F_up_H(1,1,id,kd) = -q_saf_H(kd,2)*theta_s &
                    &- q_saf_H(kd,1)*lam_H(id,kd,2)
                V2F_up_H(1,2,id,kd) = -q_saf_H(kd,1)*(1+lam_H(id,kd,3))
                V2F_up_H(1,3,id,kd) = 1 - q_saf_H(kd,1)*lam_H(id,kd,4)
                V2F_up_H(2,1,id,kd) = flux_p_H(kd,2)/(2*pi)
                V2F_up_H(2,2,id,kd) = 0.0_dp
                V2F_up_H(2,3,id,kd) = 0.0_dp
                V2F_up_H(3,1,id,kd) = lam_H(id,kd,2)
                V2F_up_H(3,2,id,kd) = 1 + lam_H(id,kd,3)
                V2F_up_H(3,3,id,kd) = lam_H(id,kd,4)
                jac_F_H(id,kd) = (flux_p_H(kd,2)/(2*pi)*&
                    &(1+lam_H(id,kd,3)))**(-1)*jac_V_H(id,kd)
                
                V2F_dn_H(1,1,id,kd) = 0.0_dp
                V2F_dn_H(1,2,id,kd) = -lam_H(id,kd,4)/(1+lam_H(id,kd,3))
                V2F_dn_H(1,3,id,kd) = 1.0_dp
                V2F_dn_H(2,1,id,kd) = 2*pi/flux_p_H(kd,2)
                V2F_dn_H(2,2,id,kd) = -2*pi/flux_p_H(kd,2) * (lam_H(id,kd,2)+&
                    &lam_H(id,kd,4)*q_saf_H(kd,2)*theta_s)/(1+lam_H(id,kd,3))
                V2F_dn_H(2,3,id,kd) = 2*pi/flux_p_H(kd,2)*q_saf_H(kd,2)*theta_s
                V2F_dn_H(3,1,id,kd) = 0.0_dp
                V2F_dn_H(3,2,id,kd) = (1-q_saf_H(kd,1)*lam_H(id,kd,4))/&
                    &(1+lam_H(id,kd,3))
                V2F_dn_H(3,3,id,kd) = q_saf_H(kd,1)
            end do
        end do
        
        ! convert to full mesh
        do jd = 1,n_par
            ! jacobian
            jac_F(jd,:) = h2f(jac_F_H(jd,:))
            
            ! transformation matrices
            do id = 1,3
                do kd = 1,3
                    V2F_up(kd,id,jd,:) = h2f(V2F_up_H(kd,id,jd,:))
                    V2F_dn(kd,id,jd,:) = h2f(V2F_dn_H(kd,id,jd,:))
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
