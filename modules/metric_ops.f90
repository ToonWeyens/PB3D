!------------------------------------------------------------------------------!
!   Variables, subroutines and  functions that have to do with the metric      !
!   elements                                                                   !
!------------------------------------------------------------------------------!
module metric_ops 
    use num_vars, only: dp, max_deriv
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud, write_out
    use str_ops, only: r2str, i2str
    use VMEC_vars, only: n_r
    use utilities, only: check_deriv
    
    implicit none
    private
    public calc_T_VC, calc_g_C, calc_jac_C, calc_g_V, T_VC, calc_jac_V, &
        &init_metric, calc_T_VF, calc_inv_met, calc_g_F, calc_jac_F, &
        &calc_f_deriv, &
        &g_V, jac_F, h_F, g_F, g_C, T_VF, T_FV, jac_V, det_T_VF, det_T_FV

    ! upper (h) and lower (g) metric factors
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: g_C(:,:,:,:,:,:,:)                                 ! in the C(ylindrical) coordinate system (FM)
    real(dp), allocatable :: g_V(:,:,:,:,:,:,:)                                 ! in the V(MEC) coordinate system (FM)
    real(dp), allocatable :: g_F(:,:,:,:,:,:,:), h_F(:,:,:,:,:,:,:)             ! in the F(lux) coordinate system (FM)
    ! upper and lower transformation matrices
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: T_VC(:,:,:,:,:,:,:)                                ! C(ylindrical) to V(MEC) (FM) (lower)
    real(dp), allocatable :: T_VF(:,:,:,:,:,:,:)                                ! V(MEC) to F(lux) (FM) (upper)
    real(dp), allocatable :: T_FV(:,:,:,:,:,:,:)                                ! V(MEC) to F(lux) (FM) (lower)
    real(dp), allocatable :: det_T_VC(:,:,:,:,:)                                ! determinant of T_VC
    real(dp), allocatable :: det_T_VF(:,:,:,:,:)                                ! determinant of T_VF
    real(dp), allocatable :: det_T_FV(:,:,:,:,:)                                ! determinant of T_FV
    real(dp), allocatable :: jac_C(:,:,:,:,:)                                   ! jacobian of C(ylindrical) coordinate system (FM)
    real(dp), allocatable :: jac_V(:,:,:,:,:)                                   ! jacobian of V(MEC) coordinate system (FM)
    real(dp), allocatable :: jac_F(:,:,:,:,:)                                   ! jacobian of F(lux) coordinate system (FM)
        
    interface calc_g_C
        module procedure calc_g_C_ind, calc_g_C_arr
    end interface
    interface calc_g_V
        module procedure calc_g_V_ind, calc_g_V_arr
    end interface
    interface calc_g_F
        module procedure calc_g_F_ind, calc_g_F_arr
    end interface
    interface calc_jac_C
        module procedure calc_jac_C_ind, calc_jac_C_arr
    end interface
    interface calc_jac_V
        module procedure calc_jac_V_ind, calc_jac_V_arr
    end interface
    interface calc_jac_F
        module procedure calc_jac_F_ind, calc_jac_F_arr
    end interface
    interface calc_T_VC
        module procedure calc_T_VC_ind, calc_T_VC_arr
    end interface
    interface calc_T_VF
        module procedure calc_T_VF_ind, calc_T_VF_arr
    end interface
    interface calc_inv_met
        module procedure calc_inv_met_ind, calc_inv_met_arr, &
            &calc_inv_met_ind_0D, calc_inv_met_arr_0D
    end interface
    interface calc_f_deriv
        module procedure calc_f_deriv_ind, calc_f_deriv_arr, &
            &calc_f_deriv_arr_2D
    end interface

contains
    ! initialize metric variables
    subroutine init_metric
        use VMEC_vars, only: n_r
        use eq_vars, only: n_par
        
        ! g_C
        if (allocated(g_C)) deallocate(g_C)
        allocate(g_C(n_par,n_r,3,3,0:max_deriv(1),0:max_deriv(2),&
            &0:max_deriv(3))); g_C = 0.0_dp
        
        ! g_V
        if (allocated(g_V)) deallocate(g_V)
        allocate(g_V(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! g_F
        if (allocated(g_F)) deallocate(g_F)
        allocate(g_F(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! h_F
        if (allocated(h_F)) deallocate(h_F)
        allocate(h_F(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! T_VC
        if (allocated(T_VC)) deallocate(T_VC)
        allocate(T_VC(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1)); T_VC = 0.0_dp
        
        ! T_VF
        if (allocated(T_VF)) deallocate(T_VF)
        allocate(T_VF(n_par,n_r,3,3,0:max_deriv(1)-2,0:max_deriv(2)-2,&
            &0:max_deriv(3)-2)); T_VC = 0.0_dp
        
        ! T_FV
        if (allocated(T_FV)) deallocate(T_FV)
        allocate(T_FV(n_par,n_r,3,3,0:max_deriv(1)-2,0:max_deriv(2)-2,&
            &0:max_deriv(3)-2)); T_VC = 0.0_dp
        
        ! det_T_VC
        if (allocated(det_T_VC)) deallocate(det_T_VC)
        allocate(det_T_VC(n_par,n_r,0:max_deriv(1),0:max_deriv(2),&
            &0:max_deriv(3)))
        
        ! det_T_VF
        if (allocated(det_T_VF)) deallocate(det_T_VF)
        allocate(det_T_VF(n_par,n_r,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! det_T_FV
        if (allocated(det_T_FV)) deallocate(det_T_FV)
        allocate(det_T_FV(n_par,n_r,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! jac_C
        if (allocated(jac_C)) deallocate(jac_C)
        allocate(jac_C(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! jac_V
        if (allocated(jac_V)) deallocate(jac_V)
        allocate(jac_V(n_par,n_r,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! jac_F
        if (allocated(jac_F)) deallocate(jac_F)
        allocate(jac_F(n_par,n_r,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
    end subroutine
    
    ! calculate the lower metric elements in the C(ylindrical) coordinate system
    subroutine calc_g_C_ind(deriv)
        use eq_vars, only: VMEC_R
        use utilities, only: arr_mult
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv,'calc_g_C')
        
        g_C(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        if (sum(deriv).eq.0) then
            g_C(:,:,1,1,deriv(1),deriv(2),deriv(3)) = 1.0_dp
            g_C(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 1.0_dp
        else
            g_C(:,:,1,1,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            g_C(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        end if
        call arr_mult(VMEC_R,VMEC_R,g_C(:,:,2,2,deriv(1),deriv(2),deriv(3)),&
            &deriv)
    end subroutine
    subroutine calc_g_C_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_g_C_ind(deriv(:,id))
        end do
    end subroutine

    ! calculate the  metric coefficients in  the V(MEC) coordinate  system using
    ! the metric  coefficients in  the C(ylindrical)  coordinate system  and the
    ! transformation matrices
    subroutine calc_g_V_ind(deriv)
        use num_vars, only: max_deriv
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv-[1,1,1],'calc_g_V')
        
        call calc_g(T_VC,g_C,g_V,deriv,max_deriv-[1,1,1])
    end subroutine
    subroutine calc_g_V_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_g_V_ind(deriv(:,id))
        end do
    end subroutine

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the metric  coefficients in  the V(MEC) coordinate  system and  the trans-
    ! formation matrices
    subroutine calc_g_F_ind(deriv)
        use num_vars, only: max_deriv
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv-[1,1,1],'calc_g_F')
        
        call calc_g(T_FV,g_V,g_F,deriv,max_deriv-[1,1,1])
    end subroutine
    subroutine calc_g_F_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_g_F_ind(deriv(:,id))
        end do
    end subroutine
    
    ! Calculate  the metric  coefficients in  a  coordinate system  B using  the
    ! metric coefficients in a coordinate system A and the transformation matrix
    ! T_BA, general treatment using the formula described in [ADD REF TO DOC]
    ! D_1^m1 D_2^m2 D_3^m3 g_B = sum_1 sum_2 sum_3 x
    !   C1 C2 C3 D^(i1,j1,k1) T_BA D^(i2,j2,k2) G_A D^(i3,j3,k3) (T_BA)^T
    ! with
    !   kl = ml-il-jl
    ! and with the sum_i  double summations, with the first one  going from 0 to
    ! and mi the second one from 0 to m-j
    ! the coefficients Ci are calculated as mi!/(i!j!(m-i-j)!) x 
    !   1   for odd values of m-i
    !   1/2 for even values of m-i
    subroutine calc_g(T_BA,g_A,g_B,deriv,max_deriv)
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        
        ! input / output
        integer, intent(in) :: deriv(3)
        integer, intent(in) :: max_deriv(3)
        real(dp), intent(in) :: T_BA(:,:,:,:,0:,0:,0:)
        real(dp), intent(in) :: g_A(:,:,:,:,0:,0:,0:)
        real(dp), intent(inout) :: g_B(:,:,:,:,0:,0:,0:)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: i1, j1, i2, j2, i3, j3, k1, k2, k3                           ! counters
        integer :: m1, m2, m3                                                   ! alias for deriv
        real(dp), allocatable :: C1(:), C2(:), C3(:)                            ! coeff. for (il)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv,'calc_g')
        
        ! set ml
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! initialize the requested derivative to zero
        g_B(:,:,:,:,m1,m2,m3) = 0.0_dp
        
        ! calculate the terms in the summation
        d1: do j1 = 0,m1                                                        ! derivatives in first coordinate
            call calc_C(m1,j1,C1)                                               ! calculate coeff. C1
            do i1 = m1-j1,0,-1
                d2: do j2 = 0,m2                                                ! derivatives in second coordinate
                    call calc_C(m2,j2,C2)                                       ! calculate coeff. C2
                    do i2 = m2-j2,0,-1
                        d3: do j3 = 0,m3                                        ! derivatives in third coordinate
                            call calc_C(m3,j3,C3)                               ! calculate coeff. C3
                            do i3 = m3-j3,0,-1
                                do kd = 1,n_r                                   ! all normal points
                                    do id = 1,n_par                             ! all parallel points
                                        k1 = m1 - j1 - i1
                                        k2 = m2 - j2 - i2
                                        k3 = m3 - j3 - i3
            ! add this term to the summation
            g_B(id,kd,:,:,m1,m2,m3) = g_B(id,kd,:,:,m1,m2,m3) + &
                &C1(i1)*C2(i2)*C3(i3) * matmul(T_BA(id,kd,:,:,i1,i2,i3),&       ! (i i i)(j j j)(k k k)
                &matmul(g_A(id,kd,:,:,j1,j2,j3),&
                &transpose(T_BA(id,kd,:,:,k1,k2,k3))))
                                    end do
                                end do
                            end do
                        end do d3
                    end do
                end do d2
            end do
        end do d1
    contains
        ! Calculate the coeff. C  at the all i values for  the current value for
        ! j, optionally making use of the coeff. C at previous value for j:
        ! - If it is the very first value (0,m), then it is just equal to 1
        ! - If  it is  the first value  in the current  i-summation, then  it is
        ! calculated from the first coeff. of the previous i-summation:
        !   C(j,m) = C(j-1,m) * (m-j+1)/j
        ! - If it is not the first  value in the current i-summation, then it is
        ! calculated from the previous value of the current i-summation
        !   C(j,i) = C(j,i+1) * (i+1)/(m-i-j)
        ! - If  it is  the last  value in  the current  i-summation, then  it is
        ! divided by 2 if (m-j) is even (see [ADD REF])
        subroutine calc_C(m,j,C)
            ! input / output
            integer, intent(in) :: m, j
            real(dp), intent(inout), allocatable :: C(:)
            
            ! local variables
            integer :: i                                                        ! counter
            real(dp) :: Cm_prev
            
            ! if j>0, save the value Cm_prev 
            if (j.gt.0) Cm_prev = C(m-j+1)
            
            ! allocate C
            if (allocated(C)) deallocate (C)
            allocate(C(0:m-j))
            
            ! first value i = m-j
            if (j.eq.0) then                                                    ! first coeff., for j = 0
                C(m) = 1.0_dp
            else                                                                ! first coeff., for j > 0 -> need value Cm_prev from previous array
                C(m-j) = (m-j+1.0_dp)/j * Cm_prev
            end if
            
            ! all other values i = m-1 .. 0
            do i = m-j-1,0,-1
                C(i) = (i+1.0_dp)/(m-i-j) * C(i+1)
            end do
        end subroutine
    end subroutine
    
    ! calculate the jacobian in cylindrical coordinates
    subroutine calc_jac_C_ind(deriv)
        use eq_vars, only: VMEC_R
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv,'calc_J_C')
        
        jac_C(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        jac_C(:,:,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_R(:,:,deriv(1),deriv(2),deriv(3))
    end subroutine
    subroutine calc_jac_C_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_jac_C_ind(deriv(:,id))
        end do
    end subroutine
    
    ! calculate the jacobian in VMEC coordinates from 
    !   jac_V = T_VC jac_C
    subroutine calc_jac_V_ind(deriv)
        use utilities, only: arr_mult
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv-[1,1,1],'calc_J_V')
        
        ! calculate determinant
        jac_V(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        call arr_mult(jac_C,det_T_VC,jac_V(:,:,deriv(1),deriv(2),deriv(3)),&
            &deriv)
    end subroutine
    subroutine calc_jac_V_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_jac_V_ind(deriv(:,id))
        end do
    end subroutine
    
    ! calculate the jacobian in Flux coordinates from 
    !   jac_F = T_FV jac_V
    subroutine calc_jac_F_ind(deriv)
        use utilities, only: arr_mult
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv-[1,1,1],'calc_J_V')
        
        ! calculate determinant
        jac_F(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        call arr_mult(jac_V,det_T_FV,jac_F(:,:,deriv(1),deriv(2),deriv(3)),&
            &deriv)
    end subroutine
    subroutine calc_jac_F_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_jac_F_ind(deriv(:,id))
        end do
    end subroutine

    ! calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system
    subroutine calc_T_VC_ind(deriv)
        use utilities, only: arr_mult
        use eq_vars, only: VMEC_R, VMEC_Z
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv-[1,1,1],'calc_T_VC')
        
        ! initialize
        T_VC(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate transformation matrix T_V^C (FM)
        T_VC(:,:,1,1,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_R(:,:,deriv(1)+1,deriv(2),deriv(3))
        T_VC(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0
        T_VC(:,:,1,3,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_Z(:,:,deriv(1)+1,deriv(2),deriv(3))
        T_VC(:,:,2,1,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_R(:,:,deriv(1),deriv(2)+1,deriv(3))
        T_VC(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0
        T_VC(:,:,2,3,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_Z(:,:,deriv(1),deriv(2)+1,deriv(3))
        T_VC(:,:,3,1,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_R(:,:,deriv(1),deriv(2),deriv(3)+1)
        if (sum(deriv).eq.0) then
            T_VC(:,:,3,2,deriv(1),deriv(2),deriv(3)) = -1.0_dp
        else
            T_VC(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0
        end if
        T_VC(:,:,3,3,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_Z(:,:,deriv(1),deriv(2),deriv(3)+1)
        
        ! determinant
        det_T_VC(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        call arr_mult(VMEC_R(:,:,1:,0:,0:),VMEC_Z(:,:,0:,1:,0:),&
            &det_T_VC(:,:,deriv(1),deriv(2),deriv(3)),deriv)
        call arr_mult(-VMEC_R(:,:,0:,1:,0:),VMEC_Z(:,:,1:,0:,0:),&
            &det_T_VC(:,:,deriv(1),deriv(2),deriv(3)),deriv)
    end subroutine
    subroutine calc_T_VC_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_T_VC_ind(deriv(:,id))
        end do
    end subroutine

    ! calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system
    subroutine calc_T_VF_ind(deriv)
        use num_vars, only: pi
        use eq_vars, only: VMEC_L, q_saf, n_par, theta, flux_p
        use utilities, only: arr_mult
        use VMEC_vars, only: n_r
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, jd, kd                                                   ! counter
        real(dp) :: theta_s(n_par,n_r,0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1)   ! straight field line coordinate theta
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv-[2,2,2],'calc_T_VF')
        
        ! initialize
        theta_s = 0.0_dp
        T_VF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate theta_s
        ! start from theta
        theta_s(:,:,0,0,0) = theta                                              ! underived
        theta_s(:,:,0,1,0) = 1.0_dp                                             ! first derivative in theta
        ! add the deformation described by lambda
        do kd = 0,deriv(3)+1                                                    ! zeta derivatives
            do jd = 0,deriv(2)+1                                                ! theta derivatives
                do id = 0,deriv(1)+1                                            ! r derivatives
                    theta_s(:,:,id,jd,kd) = &
                        &theta_s(:,:,id,jd,kd) + VMEC_L(:,:,id,jd,kd)
                end do
            end do
        end do
        
        ! calculate transformation matrix T_V^F (FM)
        ! (1,1)
        call arr_mult(theta_s,-q_saf(:,1:),&
            &T_VF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
        call arr_mult(VMEC_L(:,:,1:,0:,0:),-q_saf,&
            &T_VF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
        ! (1,2)
        if (deriv(2).eq.0 .and. deriv(3).eq.0) then
            do id = 1,n_par
                T_VF(id,:,1,2,deriv(1),0,0) = flux_p(:,deriv(1)+1)/(2*pi)
            end do
        else
            T_VF(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        end if
        ! (1,3)
        T_VF(:,:,1,3,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_L(:,:,deriv(1)+1,deriv(2),deriv(3))
        ! (2,1)
        call arr_mult(theta_s(:,:,0:,1:,0:),-q_saf,&
            &T_VF(:,:,2,1,deriv(1),deriv(2),deriv(3)),deriv)
        ! (2,2)
        T_VF(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ! (2,3)
        T_VF(:,:,2,3,deriv(1),deriv(2),deriv(3)) = &
            &theta_s(:,:,deriv(1),deriv(2)+1,deriv(3))
        ! (3,1)
        if (sum(deriv).eq.0) then
            T_VF(:,:,3,1,0,0,0) = 1.0_dp
        end if
        call arr_mult(VMEC_L(:,:,0:,0:,1:),-q_saf,&
            &T_VF(:,:,3,1,deriv(1),deriv(2),deriv(3)),deriv)
        ! (3,2)
        T_VF(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ! (3,3)
        T_VF(:,:,3,3,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_L(:,:,deriv(1),deriv(2),deriv(3)+1)
        
        ! determinant
        det_T_VF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        call arr_mult(theta_s(:,:,0:,1:,0:)/(2*pi),flux_p(:,1:),&
            &det_T_VF(:,:,deriv(1),deriv(2),deriv(3)),deriv)
    end subroutine
    subroutine calc_T_VF_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_T_VF_ind(deriv(:,id))
        end do
    end subroutine
    
    ! calculate D_1^m1 D_2^m2 D_3^m3 X from D_1^i1 D_2^i2 D_3^3 X and 
    ! D_1^j1 D_2^j2 D_3^j3 Y where XY = 1, i,j = 0..m, according to [ADD REF]
    subroutine calc_inv_met_ind(X,Y,deriv)                                      ! matrix version
        use utilities, only: inv
        use VMEC_vars, only: n_r
        use eq_vars, only: n_par
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: r, t, z                                                      ! counters for derivatives
        integer :: m1, m2, m3                                                   ! alias for deriv(i)
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! initialize the requested derivative
        X(:,:,:,:,m1,m2,m3) = 0.0_dp
        
        ! calculate terms
        if (m1.eq.0 .and. m2.eq.0 .and. m3.eq.0) then                           ! direct inverse
            X(:,:,:,:,0,0,0) = inv(Y(:,:,:,:,0,0,0))
        else                                                                    ! calculate using the other inverses
            do z = 0,m3                                                         ! derivs in third coord
                if (z.eq.0) then                                                ! first term in sum
                    bin_fac(3) = 1.0_dp
                else
                    bin_fac(3) = bin_fac(3)*(m3-(z-1))/z
                end if
                do t = 0,m2                                                     ! derivs in second coord
                    if (t.eq.0) then                                            ! first term in sum
                        bin_fac(2) = 1.0_dp
                    else
                        bin_fac(2) = bin_fac(2)*(m2-(t-1))/t
                    end if
                    do r = 0,m1                                                 ! derivs in first coord
                        if (r.eq.0) then                                        ! first term in sum
                            bin_fac(1) = 1.0_dp
                        else
                            bin_fac(1) = bin_fac(1)*(m1-(r-1))/r
                        end if
                        if (z+t+r .lt. m1+m2+m3) then                           ! only add if not all r,t,z are equal to m1,m2,m3
                            do kd = 1,n_r
                                do id = 1,n_par
                                    X(id,kd,:,:,m1,m2,m3) = &
                                        &X(id,kd,:,:,m1,m2,m3)-&
                                        &bin_fac(1)*bin_fac(2)*bin_fac(3)* &
                                        &matmul(X(id,kd,:,:,r,t,z),&
                                        &Y(id,kd,:,:,m1-r,m2-t,m3-z))
                                end do
                            end do
                        end if
                    end do
                end do
            end do
            
            ! right-multiply by X(0,0,0)
            do kd = 1,n_r
                do id = 1,n_par
                    X(id,kd,:,:,m1,m2,m3) = matmul(X(id,kd,:,:,m1,m2,m3),&
                        &X(id,kd,:,:,0,0,0))
                end do
            end do
        end if
    end subroutine
    subroutine calc_inv_met_ind_0D(X,Y,deriv)                                   ! scalar version
        use utilities, only: inv
        use VMEC_vars, only: n_r
        use eq_vars, only: n_par
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: r, t, z                                                      ! counters for derivatives
        integer :: m1, m2, m3                                                   ! alias for deriv(i)
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! initialize the requested derivative
        X(:,:,m1,m2,m3) = 0.0_dp
        
        ! calculate terms
        if (m1.eq.0 .and. m2.eq.0 .and. m3.eq.0) then                           ! direct inverse
            X(:,:,0,0,0) = 1/Y(:,:,0,0,0)
        else                                                                    ! calculate using the other inverses
            do z = 0,m3                                                         ! derivs in third coord
                if (z.eq.0) then                                                ! first term in sum
                    bin_fac(3) = 1.0_dp
                else
                    bin_fac(3) = bin_fac(3)*(m3-(z-1))/z
                end if
                do t = 0,m2                                                     ! derivs in second coord
                    if (t.eq.0) then                                            ! first term in sum
                        bin_fac(2) = 1.0_dp
                    else
                        bin_fac(2) = bin_fac(2)*(m2-(t-1))/t
                    end if
                    do r = 0,m1                                                 ! derivs in first coord
                        if (r.eq.0) then                                        ! first term in sum
                            bin_fac(1) = 1.0_dp
                        else
                            bin_fac(1) = bin_fac(1)*(m1-(r-1))/r
                        end if
                        if (z+t+r .lt. m1+m2+m3) then                           ! only add if not all r,t,z are equal to m1,m2,m3
                            do kd = 1,n_r
                                do id = 1,n_par
                                    X(id,kd,m1,m2,m3) = X(id,kd,m1,m2,m3)-&
                                        &bin_fac(1)*bin_fac(2)*bin_fac(3)* &
                                        &X(id,kd,r,t,z)*Y(id,kd,m1-r,m2-t,m3-z)
                                end do
                            end do
                        end if
                    end do
                end do
            end do
            
            ! right-multiply by X(0,0,0)
            do kd = 1,n_r
                do id = 1,n_par
                    X(id,kd,m1,m2,m3) = X(id,kd,m1,m2,m3)*X(id,kd,0,0,0)
                end do
            end do
        end if
    end subroutine
    subroutine calc_inv_met_arr(X,Y,deriv)
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_inv_met_ind(X,Y,deriv(:,id))
        end do
    end subroutine
    subroutine calc_inv_met_arr_0D(X,Y,deriv)
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_inv_met_ind_0D(X,Y,deriv(:,id))
        end do
    end subroutine
    
    ! calculate the derivatives in the flux coordinates from derivatives in VMEC
    ! coordinates. These are replaced
    subroutine calc_f_deriv_ind(X,deriv)
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:)
        
        ! CONTINUE HERE 
        ! !!!!!!!!!!!!!!!
    end subroutine
    subroutine calc_f_deriv_arr(X,deriv)                                        ! scalar version
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_f_deriv_ind(X,deriv(:,id))
        end do
        write(*,*) 'THIS IS NOT YET IMPLEMENTED. CONTINUE!!!'
    end subroutine
    subroutine calc_f_deriv_arr_2D(X,deriv)                                     ! matrix version
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id, jd, kd
        
        do kd = 1,size(X,4)
            do jd = 1,size(X,3)
                do id = 1, size(deriv,2)
                    call calc_f_deriv_ind(X(:,:,jd,kd,:,:,:),deriv(:,id))
                end do
            end do
        end do
        write(*,*) 'THIS IS NOT YET IMPLEMENTED. CONTINUE!!!'
    end subroutine
end module metric_ops
    
    !! UNUSED ALTERNATIVE FOR calc_g_ALT: USING IMIN TECHNIQUE
    !! Calculate the  metric coefficients in  the V(MEC) coordinate  system using
    !! the metric  coefficients in  the C(ylindrical)  coordinate system  and the
    !! trans- formation  matrices, general treatment using  the formula described
    !! in [ADD REF TO DOC]
    !! D_1^m1 D_2^m2 D_3^m3 g_B = sum_1 sum_2 sum_3 x
    !!   C1 C2 C3 D^(i1,j1,k1) T_BA D^(i2,j2,k2) G_A D^(i3,j3,k3) (T_BA)^T
    !! with
    !!   kl = ml-il-jl
    !! and with the sum_i  double summations, with the first one  going from 0 to
    !! and mi the second one from floor((m-j+1)/2) to m-j
    !! the coefficients Ci are calculated as mi!/(i!j!(m-i-j)!) x 
    !!   1   for odd values of m-i
    !!   1/2 for even values of m-i
    !subroutine calc_g_ALT(T_BA,g_A,g_B,deriv,max_deriv)
        !use eq_vars, only: n_par
        !use VMEC_vars, only: n_r
        
        !! input / output
        !integer, intent(in) :: deriv(3)
        !integer, intent(in) :: max_deriv(3)
        !real(dp), intent(in) :: T_BA(:,:,:,:,0:,0:,0:)
        !real(dp), intent(in) :: g_A(:,:,:,:,0:,0:,0:)
        !real(dp), intent(inout) :: g_B(:,:,:,:,0:,0:,0:)
        
        !! local variables
        !integer :: id, kd                                                       ! counters
        !integer :: i1, j1, i2, j2, i3, j3, k1, k2, k3                           ! counters
        !integer :: imin1, imin2, imin3                                          ! min. for i variable in sum
        !integer :: m1, m2, m3                                                   ! alias for deriv
        !real(dp), allocatable :: C1(:), C2(:), C3(:)                            ! coeff. for (il)
        
        !! check the derivatives requested
        !call check_deriv(deriv,max_deriv,'calc_g')
        
        !! set ml
        !m1 = deriv(1)
        !m2 = deriv(2)
        !m3 = deriv(3)
        
        !! initialize the requested derivative to zero
        !g_B(:,:,:,:,m1,m2,m3) = 0.0_dp
        
        !! calculate the terms in the summation
        !d1: do j1 = 0,m1                                                        ! derivatives in first coordinate
            !call calc_C(m1,j1,imin1,C1)                                         ! calculate coeff. C1
            !do i1 = m1-j1,imin1,-1
                !d2: do j2 = 0,m2                                                ! derivatives in second coordinate
                    !call calc_C(m2,j2,imin2,C2)                                 ! calculate coeff. C2
                    !do i2 = m2-j2,imin2,-1
                        !d3: do j3 = 0,m3                                        ! derivatives in third coordinate
                            !call calc_C(m3,j3,imin3,C3)                         ! calculate coeff. C3
                            !do i3 = m3-j3,imin3,-1
                                !do kd = 1,n_r                                   ! all normal points
                                    !do id = 1,n_par                             ! all parallel points
                                        !k1 = m1 - j1 - i1
                                        !k2 = m2 - j2 - i2
                                        !k3 = m3 - j3 - i3
            !! add this term to the summation
            !g_B(id,kd,:,:,m1,m2,m3) = g_B(id,kd,:,:,m1,m2,m3) + &
                !&C1(i1)*C2(i2)*C3(i3) * matmul(T_BA(id,kd,:,:,i1,i2,i3),&       ! (i i i)(j j j)(k k k)
                !&matmul(g_A(id,kd,:,:,j1,j2,j3),&
                !&transpose(T_BA(id,kd,:,:,k1,k2,k3)))) + &
                !&C1(i1)*C2(i2)*C3(i3) * matmul(T_BA(id,kd,:,:,i1,i2,k3),&       ! (i i k)(j j j)(k k i)
                !&matmul(g_A(id,kd,:,:,j1,j2,j3),&
                !&transpose(T_BA(id,kd,:,:,k1,k2,i3)))) + &
                !&C1(i1)*C2(i2)*C3(i3) * matmul(T_BA(id,kd,:,:,k1,i2,k3),&       ! (k i k)(j j j)(i k i)
                !&matmul(g_A(id,kd,:,:,j1,j2,j3),&
                !&transpose(T_BA(id,kd,:,:,i1,k2,i3)))) + &
                !&C1(i1)*C2(i2)*C3(i3) * matmul(T_BA(id,kd,:,:,k1,i2,i3),&       ! (k i i)(j j j)(i i k)
                !&matmul(g_A(id,kd,:,:,j1,j2,j3),&
                !&transpose(T_BA(id,kd,:,:,i1,k2,k3))))
                                    !end do
                                !end do
                            !end do
                        !end do d3
                    !end do
                !end do d2
            !end do
        !end do d1
        
        !! add the transpose of what is already there
        !do kd = 1,n_r
            !do id = 1,n_par
                !g_B(id,kd,:,:,m1,m2,m3) = g_B(id,kd,:,:,m1,m2,m3) + &
                    !&transpose(g_B(id,kd,:,:,m1,m2,m3))
            !end do
        !end do
    !contains
        !! Calculate imin and  the coeff. C at  the all i values  for the current
        !! value for j,  optionally making use of the coeff.  C at previous value
        !! for j:
        !! - If it is the very first value (0,m), then it is just equal to 1
        !! - If  it is  the first value  in the current  i-summation, then  it is
        !! calculated from the first coeff. of the previous i-summation:
        !!   C(j,m) = C(j-1,m) * (m-j+1)/j
        !! - If it is not the first  value in the current i-summation, then it is
        !! calculated from the previous value of the current i-summation
        !!   C(j,i) = C(j,i+1) * (i+1)/(m-i-j)
        !! - If  it is  the last  value in  the current  i-summation, then  it is
        !! divided by 2 if (m-j) is even (see [ADD REF])
        !subroutine calc_C(m,j,imin,C)
            !! input / output
            !integer, intent(in) :: m, j
            !integer, intent(inout) :: imin
            !real(dp), intent(inout), allocatable :: C(:)
            
            !! local variables
            !integer :: i                                                        ! counter
            !real(dp) :: Cm_prev
            
            !! calculate minimum index for i
            !imin = floor((m-j+1.0_dp)/2)
        
            !! if j>0, save the value Cm_prev 
            !if (j.gt.0) Cm_prev = C(m-j+1)
            
            !! allocate C
            !if (allocated(C)) deallocate (C)
            !allocate(C(imin:m-j))
            
            !! first value i = m-j
            !if (j.eq.0) then                                                    ! first coeff., for j = 0
                !C(m) = 1.0_dp
            !else                                                                ! first coeff., for j > 0 -> need value Cm_prev from previous array
                !C(m-j) = (m-j+1.0_dp)/j * Cm_prev
            !end if
            
            !! all other values i = m-1 .. imin
            !do i = m-j-1,imin,-1
                !C(i) = (i+1.0_dp)/(m-i-j) * C(i+1)
            !end do
            
            !! divide coeff. by 2 if m-j even
            !if (mod(m-j,2).eq.0) then
                !C(imin) = C(imin)/2.0_dp
            !end if
        !end subroutine
    !end subroutine

