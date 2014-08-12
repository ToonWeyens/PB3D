!------------------------------------------------------------------------------!
!   Variables, subroutines and  functions that have to do with the metric      !
!   elements                                                                   !
!------------------------------------------------------------------------------!
module metric_ops 
#include <PB3D_macros.h>
    use num_vars, only: dp, max_deriv, max_str_ln
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud
    use str_ops, only: r2str, i2str
    use VMEC_vars, only: n_r
    use utilities, only: check_deriv
    
    implicit none
    private
    public calc_T_VC, calc_g_C, calc_jac_C, calc_g_V, T_VC, calc_jac_V, &
        &init_metric, calc_T_VF, calc_inv_met, calc_g_F, calc_jac_F, &
        &calc_f_deriv, &
        &g_V, jac_F, h_F, g_F, g_C, T_VF, T_FV, jac_V, det_T_VF, det_T_FV, &
        &g_F_FD, h_F_FD, jac_F_FD

    ! upper (h) and lower (g) metric factors
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: g_C(:,:,:,:,:,:,:)                                 ! in the C(ylindrical) coordinate system (FM)
    real(dp), allocatable :: g_V(:,:,:,:,:,:,:)                                 ! in the V(MEC) coordinate system (FM)
    real(dp), allocatable :: g_F(:,:,:,:,:,:,:), h_F(:,:,:,:,:,:,:)             ! in the F(lux) coordinate system (FM) with derivatves in the V(MEC) system
    real(dp), allocatable :: g_F_FD(:,:,:,:,:,:,:), h_F_FD(:,:,:,:,:,:,:)       ! in the F(lux) coordinate system (FM) with derivatives in the F(lux) system
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
    real(dp), allocatable :: jac_F(:,:,:,:,:)                                   ! jacobian of F(lux) coordinate system (FM) with derivatives in the V(MEC) system
    real(dp), allocatable :: jac_F_FD(:,:,:,:,:)                                ! jacobian of F(lux) coordinate system (FM) with derivatives in the F(lux) system
        
    ! interfaces
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
        module procedure calc_f_deriv_3_ind, calc_f_deriv_3_arr, &
            &calc_f_deriv_3_arr_2D, calc_f_deriv_1_ind
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
        
        ! g_F_FD
        if (allocated(g_F_FD)) deallocate(g_F_FD)
        allocate(g_F_FD(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! h_F_FD
        if (allocated(h_F_FD)) deallocate(h_F_FD)
        allocate(h_F_FD(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
        
        ! T_VC
        if (allocated(T_VC)) deallocate(T_VC)
        allocate(T_VC(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1)); T_VC = 0.0_dp
        
        ! T_VF
        if (allocated(T_VF)) deallocate(T_VF)
        allocate(T_VF(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1)); T_VC = 0.0_dp
        
        ! T_FV
        if (allocated(T_FV)) deallocate(T_FV)
        allocate(T_FV(n_par,n_r,3,3,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1)); T_VC = 0.0_dp
        
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
        
        ! jac_F_FD
        if (allocated(jac_F_FD)) deallocate(jac_F_FD)
        allocate(jac_F_FD(n_par,n_r,0:max_deriv(1)-1,0:max_deriv(2)-1,&
            &0:max_deriv(3)-1))
    end subroutine
    
    ! calculate the lower metric elements in the C(ylindrical) coordinate system
    integer function calc_g_C_ind(deriv) result(ierr)
        use eq_vars, only: VMEC_R
        use utilities, only: arr_mult
        
        character(*), parameter :: rout_name = 'calc_g_C_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_C')
        CHCKERR('')
        
        g_C(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        if (sum(deriv).eq.0) then
            g_C(:,:,1,1,deriv(1),deriv(2),deriv(3)) = 1.0_dp
            g_C(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 1.0_dp
        else
            g_C(:,:,1,1,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            g_C(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        end if
        ierr = arr_mult(VMEC_R,VMEC_R,g_C(:,:,2,2,deriv(1),deriv(2),deriv(3)),&
            &deriv)
        CHCKERR('')
    end function calc_g_C_ind
    integer function calc_g_C_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_C_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_g_C_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_C_arr

    ! calculate the  metric coefficients in  the V(MEC) coordinate  system using
    ! the metric  coefficients in  the C(ylindrical)  coordinate system  and the
    ! transformation matrices
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_g_V_ind(deriv) result(ierr)
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_g_V_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv-[1,1,1],'calc_g_V')
        CHCKERR('')
        
        ierr = calc_g(g_C,T_VC,g_V,deriv,max_deriv-[1,1,1])
        CHCKERR('')
    end function calc_g_V_ind
    integer function calc_g_V_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_V_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_g_V_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_V_arr

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the metric  coefficients in  the V(MEC) coordinate  system and  the trans-
    ! formation matrices
    integer function calc_g_F_ind(deriv) result(ierr)
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_g_F_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv-[1,1,1],'calc_g_F')
        CHCKERR('')
        
        ierr = calc_g(g_V,T_FV,g_F,deriv,max_deriv-[1,1,1])
        CHCKERR('')
    end function calc_g_F_ind
    integer function calc_g_F_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_F_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_g_F_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_F_arr
    
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
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_g(g_A,T_BA,g_B,deriv,max_deriv) result(ierr)
        use eq_vars, only: n_par
        use VMEC_vars, only: n_r
        
        character(*), parameter :: rout_name = 'calc_g'
        
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
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g')
        CHCKERR('')
        
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
    end function calc_g
    
    ! calculate the jacobian in cylindrical coordinates
    integer function calc_jac_C_ind(deriv) result(ierr)
        use eq_vars, only: VMEC_R
        
        character(*), parameter :: rout_name = 'calc_jac_C_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_C')
        CHCKERR('')
        
        jac_C(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        jac_C(:,:,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_R(:,:,deriv(1),deriv(2),deriv(3))
    end function calc_jac_C_ind
    integer function calc_jac_C_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_C_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_C_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_C_arr
    
    ! calculate the jacobian in VMEC coordinates from 
    !   jac_V = T_VC jac_C
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_V_ind(deriv) result(ierr)
        use utilities, only: arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_V_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv-[1,1,1],'calc_J_V')
        CHCKERR('')
        
        ! calculate determinant
        jac_V(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = arr_mult(jac_C,det_T_VC,jac_V(:,:,deriv(1),deriv(2),deriv(3)),&
            &deriv)
        CHCKERR('')
    end function calc_jac_V_ind
    integer function calc_jac_V_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_V_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_V_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_V_arr
    
    ! calculate the jacobian in Flux coordinates from 
    !   jac_F = T_FV jac_V
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_F_ind(deriv) result(ierr)
        use utilities, only: arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_F_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv-[1,1,1],'calc_J_V')
        CHCKERR('')
        
        ! calculate determinant
        jac_F(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = arr_mult(jac_V,det_T_FV,jac_F(:,:,deriv(1),deriv(2),deriv(3)),&
            &deriv)
        CHCKERR('')
    end function calc_jac_F_ind
    integer function calc_jac_F_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_F_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_F_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_F_arr

    ! calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system
    integer function calc_T_VC_ind(deriv) result(ierr)
        use utilities, only: arr_mult
        use eq_vars, only: VMEC_R, VMEC_Z
        
        character(*), parameter :: rout_name = 'calc_T_VC_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv-[1,1,1],'calc_T_VC')
        CHCKERR('')
        
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
        ierr = arr_mult(VMEC_R(:,:,1:,0:,0:),VMEC_Z(:,:,0:,1:,0:),&
            &det_T_VC(:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ierr = arr_mult(-VMEC_R(:,:,0:,1:,0:),VMEC_Z(:,:,1:,0:,0:),&
            &det_T_VC(:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_T_VC_ind
    integer function calc_T_VC_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VC_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VC_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VC_arr

    ! calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system
    integer function calc_T_VF_ind(deriv) result(ierr)
        use num_vars, only: pi
        use eq_vars, only: VMEC_L, q_saf, n_par, theta, flux_p
        use utilities, only: arr_mult
        use VMEC_vars, only: n_r
        
        character(*), parameter :: rout_name = 'calc_T_VF_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp) :: theta_s(n_par,n_r,0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1)   ! straight field line coordinate theta
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv-[1,1,1],'calc_T_VF')
        CHCKERR('')
        
        ! initialize
        theta_s = 0.0_dp
        T_VF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate theta_s
        ! start from theta
        theta_s(:,:,0,0,0) = theta                                              ! underived
        theta_s(:,:,0,1,0) = 1.0_dp                                             ! first derivative in theta
        ! add the deformation described by lambda
        theta_s = theta_s + VMEC_L(:,:,0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1)
        
        ! calculate transformation matrix T_V^F (FM)
        ! (1,1)
        ierr = arr_mult(theta_s,-q_saf(:,1:),&
            &T_VF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ierr = arr_mult(VMEC_L(:,:,1:,0:,0:),-q_saf,&
            &T_VF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
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
        ierr = arr_mult(theta_s(:,:,0:,1:,0:),-q_saf,&
            &T_VF(:,:,2,1,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ! (2,2)
        T_VF(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ! (2,3)
        T_VF(:,:,2,3,deriv(1),deriv(2),deriv(3)) = &
            &theta_s(:,:,deriv(1),deriv(2)+1,deriv(3))
        ! (3,1)
        if (sum(deriv).eq.0) then
            T_VF(:,:,3,1,0,0,0) = 1.0_dp
        end if
        ierr = arr_mult(VMEC_L(:,:,0:,0:,1:),-q_saf,&
            &T_VF(:,:,3,1,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ! (3,2)
        T_VF(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ! (3,3)
        T_VF(:,:,3,3,deriv(1),deriv(2),deriv(3)) = &
            &VMEC_L(:,:,deriv(1),deriv(2),deriv(3)+1)
        
        ! determinant
        det_T_VF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = arr_mult(theta_s(:,:,0:,1:,0:)/(2*pi),flux_p(:,1:),&
            &det_T_VF(:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_T_VF_ind
    integer function calc_T_VF_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VF_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VF_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VF_arr
    
    ! calculate D_1^m1 D_2^m2 D_3^m3 X from D_1^i1 D_2^i2 D_3^3 X and 
    ! D_1^j1 D_2^j2 D_3^j3 Y where XY = 1, i,j = 0..m, according to [ADD REF]
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_inv_met_ind(X,Y,deriv) result(ierr)                   ! matrix version
        use utilities, only: calc_inv
        use VMEC_vars, only: n_r
        use eq_vars, only: n_par
        
        character(*), parameter :: rout_name = 'calc_inv_met_ind'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: r, t, z                                                      ! counters for derivatives
        integer :: m1, m2, m3                                                   ! alias for deriv(i)
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        
        ! initialize ierr
        ierr = 0
        
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! initialize the requested derivative
        X(:,:,:,:,m1,m2,m3) = 0.0_dp
        
        ! calculate terms
        if (m1.eq.0 .and. m2.eq.0 .and. m3.eq.0) then                           ! direct inverse
            ierr = calc_inv(X(:,:,:,:,0,0,0),Y(:,:,:,:,0,0,0))
            CHCKERR('')
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
    end function calc_inv_met_ind
    integer function calc_inv_met_ind_0D(X,Y,deriv) result(ierr)                ! scalar version
        use VMEC_vars, only: n_r
        use eq_vars, only: n_par
        character(*), parameter :: rout_name = 'calc_inv_met_ind_0D'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: r, t, z                                                      ! counters for derivatives
        integer :: m1, m2, m3                                                   ! alias for deriv(i)
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        
        ! initialize ierr
        ierr = 0
        
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
                    ! ADD SOME ERROR CHECKING HERE !!!!
                    CHCKERR('')
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
    end function calc_inv_met_ind_0D
    integer function calc_inv_met_arr(X,Y,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_met_arr'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_inv_met_ind(X,Y,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_inv_met_arr
    integer function calc_inv_met_arr_0D(X,Y,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_met_arr_0D'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_inv_met_ind_0D(X,Y,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_inv_met_arr_0D
    
    ! calculate the derivatives in the flux coordinates from derivatives in VMEC
    ! coordinates.
    ! The routine works by exchangeing the  derivatives in the coordinates B for
    ! derivatives in coordinates A using the formula
    !   D^m_B X = T_BA D_A (D^m-1_B X)
    ! This is done for the derivatives in each of the coordinates B until degree
    ! 0  is  reached. Furthermore,  each  of  these  degrees of  derivatives  in
    ! coordinates  B  has  to  be   derived  optionally  also  in  the  original
    ! coordinates A, which yields the formula:
    !   D^p_A (D^m_B X) = sum_q binom(p,q) D^q_A (T_BA) D^p-q_A D_A (D^m-1_B X)
    ! This way, ultimately  the desired derivatives in the coordinates  B can be
    ! obtained recursively from the lower orders in the coordinates B and higher
    ! orders in the coordinates A
    ! see [ADD REF] for more detailed information
    integer recursive function calc_f_deriv_3_ind(X_A,T_BA,X_B,max_deriv,&
        &deriv_B,deriv_A_input) result(ierr)                                    ! normal variable version
        use utilities, only: arr_mult
        
        character(*), parameter :: rout_name = 'calc_f_deriv_3_ind'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,0:,0:,0:)                             ! variable and derivs. in coord. system A
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        real(dp), intent(inout) :: X_B(1:,1:)                                   ! requested derivs. of variable in coord. system B
        integer, intent(in) :: max_deriv(:)                                     ! maximum degrees of derivs.
        integer, intent(in) :: deriv_B(:)                                       ! derivs. in coord. system B
        integer, intent(in), optional :: deriv_A_input(:)                       ! derivs. in coord. system A (optional)
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: deriv_id_B                                                   ! holds the deriv. currently calculated
        integer, allocatable :: deriv_A(:)                                      ! holds either deriv_A or 0
        real(dp), allocatable :: X_B_x(:,:,:,:,:)                               ! X_B and derivs. D^p-q_A with extra 1 exchanged deriv. D_A
        integer, allocatable :: deriv_A_x(:)                                    ! holds A derivs. for exchanged X_B_x
        integer, allocatable :: deriv_B_x(:)                                    ! holds B derivs. for exchanged X_B_x
        
        ! initialize ierr
        ierr = 0
        
        ! setup deriv_A
        allocate(deriv_A(size(deriv_B)))
        if (present(deriv_A_input)) then
            deriv_A = deriv_A_input
        else
            deriv_A = 0
        end if
        
        ! check the derivatives requested
        ! (every B deriv. needs all the A derivs. -> sum(deriv_B) needed)
        ierr = check_deriv(deriv_A +sum(deriv_B),max_deriv,'calc_f_deriv')
        CHCKERR('')
        
        ! detect first deriv. in the B coord. system that needs to be exchanged,
        ! with derivs in the A coord. system going from B coord. 1 to B coord. 2
        ! to B coord. 3
        ! if equal to  zero on termination, this means that  no more exchange is
        ! necessary
        deriv_id_B = 0
        do id = 1,3
            if (deriv_B(id).gt.0) then
                deriv_id_B = id
                exit
            end if
        end do
        
        ! calculate the  derivative in coord.  deriv_id_B of coord. system  B of
        ! one order lower than requested here
        ! check if we have reached the zeroth derivative in coord. B
        if (deriv_id_B.eq.0) then                                               ! return the function with its requierd A derivs.
            X_B = X_A(:,:,deriv_A(1),deriv_A(2),deriv_A(3))
        else                                                                    ! apply the formula to obtain X_B using exchanged derivs.
            ! allocate  the  exchanged  X_B_x  and  the  derivs.  in  A  and  B,
            ! corresponding to D^m-1_B X
            allocate(X_B_x(size(X_B,1),size(X_B,2),0:deriv_A(1),0:deriv_A(2),&
                &0:deriv_A(3)))
            allocate(deriv_A_x(size(deriv_A)))
            allocate(deriv_B_x(size(deriv_B)))
            
            ! initialize X_B
            X_B = 0.0_dp
            
            ! loop over the 3 terms in the sum due to the transf. mat.
            do id = 1,3
                ! calculate X_B_x for orders 0 up to deriv_A
                do jd = 0,deriv_A(1)
                    do kd = 0,deriv_A(2)
                        do ld = 0,deriv_A(3)
                            ! calculate the  derivs. in  A for X_B_x:  for every
                            ! component in the sum due  to the transf. mat., the
                            ! deriv. is one order  higher than [jd,kd,ld] at the
                            ! corresponding coord. id
                            deriv_A_x = [jd,kd,ld]
                            deriv_A_x(id) = deriv_A_x(id) + 1
                            
                            ! calculate the derivs. in B  for X_B_x: this is one
                            ! order lower in the coordinate deriv_id_B
                            deriv_B_x = deriv_B
                            deriv_B_x(deriv_id_B) = deriv_B_x(deriv_id_B) - 1
                            
                            ! recursively call the subroutine again to calculate
                            ! X_B_x for the current derivatives
                            ierr = calc_f_deriv_3_ind(X_A,T_BA,&
                                &X_B_x(:,:,jd,kd,ld),max_deriv,deriv_B_x,&
                                &deriv_A_x)
                            CHCKERR('')
                        end do
                    end do
                end do
                
                ! use X_B_x at this term in summation due to the transf. mat. to
                ! update X_B using the formula
                ierr = arr_mult(T_BA(:,:,deriv_id_B,id,0:,0:,0:),&
                    &X_B_x(:,:,0:,0:,0:),X_B,deriv_A)
                CHCKERR('')
            end do
        end if
    end function calc_f_deriv_3_ind
    integer recursive function calc_f_deriv_1_ind(X_A,T_BA,X_B,max_deriv,&
        &deriv_B,deriv_A_input) result(ierr)                                    ! flux variable version
        use utilities, only: arr_mult
        
        character(*), parameter :: rout_name = 'calc_f_deriv_1_ind'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,0:)                                      ! variable and derivs. in coord. system A
        real(dp), intent(inout) :: X_B(1:)                                      ! requested derivs. of variable in coord. system B
        real(dp), intent(in) :: T_BA(1:,0:)                                     ! transf. mat. and derivs. between coord. systems A and B
        integer, intent(in), optional :: deriv_A_input                          ! derivs. in coord. system A (optional)
        integer, intent(in) :: deriv_B                                          ! derivs. in coord. system B
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
        
        ! local variables
        integer :: jd                                                           ! counters
        integer :: deriv_A                                                      ! holds either deriv_A or 0
        real(dp), allocatable :: X_B_x(:,:)                                     ! X_B and derivs. D^p-q_A with extra 1 exchanged deriv. D_A
        integer :: deriv_A_x                                                    ! holds A derivs. for exchanged X_B_x
        integer :: deriv_B_x                                                    ! holds B derivs. for exchanged X_B_x
        
        ! initialize ierr
        ierr = 0
        
        ! setup deriv_A
        if (present(deriv_A_input)) then
            deriv_A = deriv_A_input
        else
            deriv_A = 0
        end if
        
        ! check the derivatives requested
        ! (every B deriv. needs all the A derivs. -> sum(deriv_B) needed)
        ierr = check_deriv([deriv_A+deriv_B,0,0],[max_deriv,0,0],'calc_f_deriv')
        CHCKERR('')
        
        ! calculate the  derivative in coord.  deriv_id_B of coord. system  B of
        ! one order lower than requested here
        ! check if we have reached the zeroth derivative in coord. B
        if (deriv_B.eq.0) then                                                  ! return the function with its requierd A derivs.
            X_B = X_A(:,deriv_A)
        else                                                                    ! apply the formula to obtain X_B using exchanged derivs.
            ! allocate  the  exchanged  X_B_x  and  the  derivs.  in  A  and  B,
            ! corresponding to D^m-1_B X
            allocate(X_B_x(size(X_B,1),0:deriv_A))
            
            ! initialize X_B
            X_B = 0.0_dp
            
            ! calculate X_B_x for orders 0 up to deriv_A
            do jd = 0,deriv_A
                ! calculate the derivs.  in A for X_B_x: for  every component in
                ! the  sum due  to the  transf. mat.,  the deriv.  is one  order
                ! higher than jd
                deriv_A_x = jd+1
                
                ! calculate the derivs. in B for X_B_x: this is one order lower
                deriv_B_x = deriv_B-1
                
                ! recursively call  the subroutine again to  calculate X_B_x for
                ! the current derivatives
                ierr = calc_f_deriv_1_ind(X_A,T_BA,X_B_x(:,jd),max_deriv,&
                    &deriv_B_x,deriv_A_input=deriv_A_x)
                CHCKERR('')
            end do
            
            ! use X_B_x to calculate X_B using the formula
            ierr = arr_mult(T_BA(:,0:),X_B_x(:,0:),X_B,[deriv_A,0,0])
            CHCKERR('')
        end if
    end function calc_f_deriv_1_ind
    integer function calc_f_deriv_3_arr(X_A,T_BA,X_B,max_deriv,derivs) &
        &result(ierr)                                                           ! matrix version
        character(*), parameter :: rout_name = 'calc_f_deriv_3_arr'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,0:,0:,0:)                             ! variable and derivs. in coord. system A
        real(dp), intent(inout) :: X_B(1:,1:,0:,0:,0:)                          ! requested derivs. of variable in coord. system B
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        integer, intent(in) :: derivs(:,:)                                      ! series of derivs. (in coordinate system B)
        integer, intent(in) :: max_deriv(:)                                     ! maximum degrees of derivs.
        
        ! local variables
        integer :: id                                                           ! counter
        
        do id = 1, size(derivs,2)
            ierr = calc_f_deriv_3_ind(X_A,T_BA,&
                &X_B(:,:,derivs(1,id),derivs(2,id),derivs(3,id)),&
                &max_deriv,derivs(:,id))
            CHCKERR('')
        end do
    end function calc_f_deriv_3_arr
    integer function calc_f_deriv_3_arr_2D(X_A,T_BA,X_B,max_deriv,derivs) &
        &result(ierr)                                                           ! matrix version
        character(*), parameter :: rout_name = 'calc_f_deriv_3_arr_2D'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,1:,1:,0:,0:,0:)                       ! variable and derivs. in coord. system A
        real(dp), intent(inout) :: X_B(1:,1:,1:,1:,0:,0:,0:)                    ! requested derivs. of variable in coord. system B
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        integer, intent(in) :: derivs(:,:)                                      ! series of derivs. (in coordinate system B)
        integer, intent(in) :: max_deriv(:)                                     ! maximum degrees of derivs.
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! test whether the dimensions of X_A and X_B agree
        if (size(X_A,3).ne.size(X_B,3) .or. size(X_A,4).ne.size(X_B,4)) then
            err_msg = 'X_A and X_B need to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        do kd = 1,size(X_A,4)
            do jd = 1,size(X_A,3)
                do id = 1, size(derivs,2)
                    ierr = calc_f_deriv_3_ind(X_A(:,:,jd,kd,:,:,:),T_BA,&
                        &X_B(:,:,jd,kd,derivs(1,id),derivs(2,id),derivs(3,id)),&
                        &max_deriv,derivs(:,id))
                    CHCKERR('')
                end do
            end do
        end do
    end function calc_f_deriv_3_arr_2D
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
        !ierr = check_deriv(deriv,max_deriv,'calc_g')
        
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

