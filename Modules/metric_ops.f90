!------------------------------------------------------------------------------!
!   Operations that have to do with the metric elements                        !
!------------------------------------------------------------------------------!
module metric_ops 
#include <PB3D_macros.h>
    use num_vars, only: dp, max_deriv, max_str_ln
    use messages, only: writo, print_ar_2, print_ar_1, lvl_ud
    use str_ops, only: r2str, i2str
    use eq_vars, only: grp_n_r_eq
    use utilities, only: check_deriv
    use output_ops, only: print_GP_3D, print_GP_2D
    
    implicit none
    private
    public calc_T_VC, calc_g_C, calc_jac_C, calc_g_V, calc_jac_V, &
        &calc_jac_H, calc_T_HF, init_metric, calc_T_VF, calc_h_H, &
        &calc_inv_met, calc_g_F, calc_jac_F, normalize_metric_vars, &
        &calc_f_deriv
    
    ! interfaces
    interface calc_g_C
        module procedure calc_g_C_ind, calc_g_C_arr
    end interface
    interface calc_g_V
        module procedure calc_g_V_ind, calc_g_V_arr
    end interface
    interface calc_h_H
        module procedure calc_h_H_ind, calc_h_H_arr
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
    interface calc_jac_H
        module procedure calc_jac_H_ind, calc_jac_H_arr
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
    interface calc_T_HF
        module procedure calc_T_HF_ind, calc_T_HF_arr
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
    integer function init_metric() result(ierr)
        use num_vars, only: eq_style
        use eq_vars, only: grp_n_r_eq
        use X_vars, only: n_par
        use metric_vars, only: g_E, h_E, g_F, h_F, g_FD, h_FD, jac_E, jac_F, &
            &jac_FD, T_EF, T_FE, det_T_EF, det_T_FE, g_C, T_VC, det_T_VC, jac_C
        
        character(*), parameter :: rout_name = 'init_metric'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize variables that are used for all equilibrium styles
        ! g_E
        allocate(g_E(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! h_E
        allocate(h_E(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! g_F
        allocate(g_F(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! h_F
        allocate(h_F(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! g_FD
        allocate(g_FD(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! h_FD
        allocate(h_FD(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! jac_E
        allocate(jac_E(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! jac_F
        allocate(jac_F(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! jac_FD
        allocate(jac_FD(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,0:max_deriv))
        
        ! T_EF
        allocate(T_EF(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,&
            &0:max_deriv)); T_EF = 0.0_dp
        
        ! T_FE
        allocate(T_FE(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,&
            &0:max_deriv)); T_FE = 0.0_dp
        
        ! det_T_EF
        allocate(det_T_EF(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,&
            &0:max_deriv))
        
        ! det_T_FE
        allocate(det_T_FE(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,&
            &0:max_deriv))
        
        ! initialize variables that are  specifici to which equilibrium style is
        ! being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! g_C
                allocate(g_C(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,&
                    &0:max_deriv)); g_C = 0.0_dp
                
                ! T_VC
                allocate(T_VC(n_par,grp_n_r_eq,3,3,0:max_deriv,0:max_deriv,&
                    &0:max_deriv)); T_VC = 0.0_dp
                
                ! det_T_VC
                allocate(det_T_VC(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,&
                    &0:max_deriv))
                
                ! jac_C
                allocate(jac_C(n_par,grp_n_r_eq,0:max_deriv,0:max_deriv,&
                    &0:max_deriv))
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function init_metric
    
    ! calculate the lower metric elements in the C(ylindrical) coordinate system
    integer function calc_g_C_ind(deriv) result(ierr)
        use eq_vars, only: R_E
        use utilities, only: add_arr_mult
        use metric_vars, only: g_C
        
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
        end if
        ierr = add_arr_mult(R_E,R_E,g_C(:,:,2,2,deriv(1),deriv(2),deriv(3)),&
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

    ! calculate  the metric  coefficients in  the equilibrium  V(MEC) coordinate
    ! system  using  the metric  coefficients  in  the C(ylindrical)  coordinate
    ! system and the transformation matrices
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_g_V_ind(deriv) result(ierr)
        use num_vars, only: max_deriv
        use metric_vars, only: g_C, T_VC, g_E
        
        character(*), parameter :: rout_name = 'calc_g_V_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_V')
        CHCKERR('')
        
        ierr = calc_g(g_C,T_VC,g_E,deriv,max_deriv)
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
    
    ! calculate the  metric coefficients in the  equilibrium H(ELENA) coordinate
    ! system using the HELENA output
    integer function calc_h_H_ind(deriv) result(ierr)
        use num_vars, only: max_deriv
        use HELENA, only: h_H_11, h_H_12, h_H_33, flux_H
        use eq_vars, only: grp_min_r_eq, grp_max_r_eq, theta_E
        use X_vars, only: n_par
        use utilities, only: calc_deriv, calc_det
        use metric_vars, only: h_E, jac_E
        
        character(*), parameter :: rout_name = 'calc_h_H_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, kd                                                       ! counters
        !real(dp), allocatable :: jac_E_alt(:,:)                                 ! jac_E calculated alternatively
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_h_H')
        CHCKERR('')
        
        ! initialize h_E
        h_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        if (deriv(3).ne.0) then                                                 ! axisymmetry: deriv. in tor. coord. is zero
            !h_H(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (sum(deriv).eq.0) then                                          ! no derivatives
            h_E(:,:,1,1,0,0,0) = h_H_11
            h_E(:,:,1,2,0,0,0) = h_H_12
            h_E(:,:,2,1,0,0,0) = h_E(:,:,1,2,0,0,0)
            h_E(:,:,3,3,0,0,0) = h_H_33
            h_E(:,:,2,2,0,0,0) = 1._dp/h_H_11 * &
                &(1._dp/(jac_E(:,:,0,0,0)**2*h_H_33) + h_H_12**2)
            !! check if 1/jac^2 is recovered as determinant
            !allocate(jac_E_alt(size(jac_E,1),size(jac_E,2)))
            !ierr = calc_det(jac_E_alt,h_E(:,:,:,:,0,0,0))
            !CHCKERR('')
            !call print_GP_3D('jac_E','',reshape([1._dp/(jac_E(:,:,0,0,0)**2),&
                !&jac_E_alt],[shape(jac_E_alt),2]))
            !call print_GP_3D('jac_E (diff)','',1._dp/(jac_E(:,:,0,0,0)**2)-&
                !&jac_E_alt)
        else if (deriv(1).eq.1 .and. deriv(2).eq.0) then                        ! derivative in norm. coord.
            do id = 1,n_par
                ierr = calc_deriv(h_E(id,:,1,1,0,0,0),h_E(id,:,1,1,1,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,2)                     ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,1,2,0,0,0),h_E(id,:,1,2,1,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,2)                     ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,2,2,0,0,0),h_E(id,:,2,2,1,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,2)                     ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,3,3,0,0,0),h_E(id,:,3,3,1,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,2)                     ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
            end do
            h_E(:,:,2,1,1,0,0) = h_E(:,:,1,2,1,0,0)
        else if (deriv(1).eq.2 .and. deriv(2).eq.0) then                        ! 2nd derivative in norm. coord.
            do id = 1,n_par
                ierr = calc_deriv(h_E(id,:,1,1,0,0,0),h_E(id,:,1,1,2,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),2,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,1,2,0,0,0),h_E(id,:,1,2,2,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),2,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,2,2,0,0,0),h_E(id,:,2,2,2,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),2,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,3,3,0,0,0),h_E(id,:,3,3,2,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),2,1)
                CHCKERR('')
            end do
            h_E(:,:,2,1,2,0,0) = h_E(:,:,1,2,2,0,0)
        else if (deriv(1).eq.0 .and. deriv(2).eq.1) then                        ! derivative in pol. coord.
            do kd = 1,grp_n_r_eq
                ierr = calc_deriv(h_E(:,kd,1,1,0,0,0),h_E(:,kd,1,1,0,1,0),&
                    &theta_E(:,kd),1,2)                                         ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
                ierr = calc_deriv(h_E(:,kd,1,2,0,0,0),h_E(:,kd,1,2,0,1,0),&
                    &theta_E(:,kd),1,2)                                         ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
                ierr = calc_deriv(h_E(:,kd,2,2,0,0,0),h_E(:,kd,2,2,0,1,0),&
                    &theta_E(:,kd),1,2)                                         ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
                ierr = calc_deriv(h_E(:,kd,3,3,0,0,0),h_E(:,kd,3,3,0,1,0),&
                    &theta_E(:,kd),1,2)                                         ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
            end do
            h_E(:,:,2,1,0,1,0) = h_E(:,:,1,2,0,1,0)
        else if (deriv(1).eq.0 .and. deriv(2).eq.2) then                        ! 2nd derivative in pol. coord.
            do kd = 1,grp_n_r_eq
                ierr = calc_deriv(h_E(:,kd,1,1,0,0,0),h_E(:,kd,1,1,0,2,0),&
                    &theta_E(:,kd),2,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(:,kd,1,2,0,0,0),h_E(:,kd,1,2,0,2,0),&
                    &theta_E(:,kd),2,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(:,kd,2,2,0,0,0),h_E(:,kd,2,2,0,2,0),&
                    &theta_E(:,kd),2,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(:,kd,3,3,0,0,0),h_E(:,kd,3,3,0,2,0),&
                    &theta_E(:,kd),2,1)
                CHCKERR('')
            end do
            h_E(:,:,2,1,0,2,0) = h_E(:,:,1,2,0,2,0)
        else if (deriv(1).eq.1 .and. deriv(2).eq.1) then                        ! mixed derivative in norm. and pol. coord.
            do id = 1,n_par
                ierr = calc_deriv(h_E(id,:,1,1,0,1,0),h_E(id,:,1,1,1,1,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,1,2,0,1,0),h_E(id,:,1,2,1,1,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,2,2,0,1,0),h_E(id,:,2,2,1,1,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,1)
                CHCKERR('')
                ierr = calc_deriv(h_E(id,:,3,3,0,1,0),h_E(id,:,3,3,1,1,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,1)
                CHCKERR('')
            end do
            h_E(:,:,2,1,1,1,0) = h_E(:,:,1,2,1,1,0)
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_h_H_ind
    integer function calc_h_H_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_h_H_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_h_H_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_h_H_arr

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the  metric coefficients  in  the equilibrium  coordinate  system and  the
    ! trans- formation matrices
    integer function calc_g_F_ind(deriv) result(ierr)
        use num_vars, only: max_deriv
        use metric_vars, only: g_E, g_F, T_FE
        
        character(*), parameter :: rout_name = 'calc_g_F_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_F')
        CHCKERR('')
        
        ! calculate g_F
        ierr = calc_g(g_E,T_FE,g_F,deriv,max_deriv)
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
    ! mi and the second one from m-j to 0
    ! the coefficients Ci are calculated as mi!/(i!j!(m-i-j)!)
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect. This is not checked!
    integer function calc_g(g_A,T_BA,g_B,deriv,max_deriv) result(ierr)
        use eq_vars, only: grp_n_r_eq
        use X_vars, only: n_par
        
        character(*), parameter :: rout_name = 'calc_g'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        integer, intent(in) :: max_deriv
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
                                do kd = 1,grp_n_r_eq                            ! all normal points
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
        use eq_vars, only: R_E
        use metric_vars, only: jac_C
        
        character(*), parameter :: rout_name = 'calc_jac_C_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_C')
        CHCKERR('')
        
        jac_C(:,:,deriv(1),deriv(2),deriv(3)) = &
            &R_E(:,:,deriv(1),deriv(2),deriv(3))
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
    
    ! calculate the jacobian in equilibrium V(MEC) coordinates from 
    !   jac_V = det(T_VC) jac_C
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_V_ind(deriv) result(ierr)
        use utilities, only: add_arr_mult
        use metric_vars, only: jac_E, det_T_VC, jac_C
        
        character(*), parameter :: rout_name = 'calc_jac_V_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_V')
        CHCKERR('')
        
        ! calculate determinant
        jac_E(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = add_arr_mult(jac_C,det_T_VC,&
            &jac_E(:,:,deriv(1),deriv(2),deriv(3)),deriv)
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
    
    ! calculate the jacobian in HELENA coordinates directly from J = q R^2/F
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_H_ind(deriv) result(ierr)
        use HELENA, only:  qs, h_H_33, RBphi, flux_H
        use utilities, only: calc_deriv
        use X_vars, only: n_par
        use eq_vars, only: grp_n_r_eq, grp_min_r_eq, theta_E, &
            &grp_min_r_eq, grp_max_r_eq
        use metric_vars, only: jac_E
        
        character(*), parameter :: rout_name = 'calc_jac_H_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_H')
        CHCKERR('')
        
        ! calculate determinant
        if (deriv(3).ne.0) then                                                 ! axisymmetry: deriv. in tor. coord. is zero
            jac_E(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (sum(deriv).eq.0) then                                          ! no derivatives
            do kd = 1,grp_n_r_eq
                jac_E(:,kd,0,0,0) = 1.0_dp/(h_H_33(:,kd)*&
                    &qs(grp_min_r_eq-1+kd)*RBphi(grp_min_r_eq-1+kd))
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.0) then                        ! derivative in norm. coord.
            do id = 1,n_par
                ierr = calc_deriv(jac_E(id,:,0,0,0),jac_E(id,:,1,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),1,2)                     ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
            end do
        else if (deriv(1).eq.2 .and. deriv(2).eq.0) then                        ! 2nd derivative in norm. coord.
            do id = 1,n_par
                ierr = calc_deriv(jac_E(id,:,0,0,0),jac_E(id,:,2,0,0),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),2,1)
                CHCKERR('')
            end do
        else if (deriv(1).eq.0 .and. deriv(2).eq.1) then                        ! derivative in pol. coord.
            do kd = 1,grp_n_r_eq
                ierr = calc_deriv(jac_E(:,kd,0,0,0),jac_E(:,kd,0,1,0),&
                    &theta_E(:,kd),1,2)                                         ! use extra precision to later calculate mixed derivatives
                CHCKERR('')
            end do
        else if (deriv(1).eq.0 .and. deriv(2).eq.2) then                        ! 2nd derivative in pol. coord.
            do kd = 1,grp_n_r_eq
                ierr = calc_deriv(jac_E(:,kd,0,0,0),jac_E(:,kd,0,2,0),&
                    &theta_E(:,kd),2,1)
                CHCKERR('')
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.1) then                        ! mixed derivative in norm. and pol. coord.
            do kd = 1,grp_n_r_eq
                ierr = calc_deriv(jac_E(:,kd,1,0,0),jac_E(:,kd,1,1,0),&
                    &theta_E(:,kd),1,1)     
                CHCKERR('')
            end do
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_jac_H_ind
    integer function calc_jac_H_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_H_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_H_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_H_arr
    
    ! calculate the jacobian in Flux coordinates from 
    !   jac_F = det(T_FE) jac_E
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_F_ind(deriv) result(ierr)
        use utilities, only: add_arr_mult
        use metric_vars, only: jac_F, det_T_FE, jac_E
        
        character(*), parameter :: rout_name = 'calc_jac_F_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_V')
        CHCKERR('')
        
        ! initialize jac_F
        jac_F(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate jac_F
        ierr = add_arr_mult(jac_E,det_T_FE,&
            &jac_F(:,:,deriv(1),deriv(2),deriv(3)),deriv)
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
        use utilities, only: add_arr_mult
        use eq_vars, only: R_E, Z_E
        use metric_vars, only: T_VC, det_T_VC
        
        character(*), parameter :: rout_name = 'calc_T_VC_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VC')
        CHCKERR('')
        
        ! initialize
        T_VC(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate transformation matrix T_V^C
        T_VC(:,:,1,1,deriv(1),deriv(2),deriv(3)) = &
            &R_E(:,:,deriv(1)+1,deriv(2),deriv(3))
        !T_VC(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0
        T_VC(:,:,1,3,deriv(1),deriv(2),deriv(3)) = &
            &Z_E(:,:,deriv(1)+1,deriv(2),deriv(3))
        T_VC(:,:,2,1,deriv(1),deriv(2),deriv(3)) = &
            &R_E(:,:,deriv(1),deriv(2)+1,deriv(3))
        !T_VC(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0
        T_VC(:,:,2,3,deriv(1),deriv(2),deriv(3)) = &
            &Z_E(:,:,deriv(1),deriv(2)+1,deriv(3))
        T_VC(:,:,3,1,deriv(1),deriv(2),deriv(3)) = &
            &R_E(:,:,deriv(1),deriv(2),deriv(3)+1)
        if (sum(deriv).eq.0) then
            T_VC(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 1.0_dp
        !else
            !T_VC(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0
        end if
        T_VC(:,:,3,3,deriv(1),deriv(2),deriv(3)) = &
            &Z_E(:,:,deriv(1),deriv(2),deriv(3)+1)
        
        ! determinant
        det_T_VC(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = add_arr_mult(R_E(:,:,0:,1:,0:),Z_E(:,:,1:,0:,0:),&
            &det_T_VC(:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ierr = add_arr_mult(-R_E(:,:,1:,0:,0:),Z_E(:,:,0:,1:,0:),&
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

    ! calculate the transformation matrix  between equilibrium V(mec) and F(lux)
    ! oordinate system
    integer function calc_T_VF_ind(deriv) result(ierr)
        use num_vars, only: pi, use_pol_flux_X
        use X_vars, only: n_par
        use eq_vars, only: L_E, q_saf_E, rot_t_E, theta_E, zeta_E, flux_p_E, &
            &flux_t_E, grp_n_r_eq
        use utilities, only: add_arr_mult
        use metric_vars, only: T_EF, det_T_EF
        
        character(*), parameter :: rout_name = 'calc_T_VF_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:)                             ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:)                              ! - zeta_F and derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VF')
        CHCKERR('')
        
        ! initialize T_EF
        T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! set up theta_s
        allocate(theta_s(n_par,grp_n_r_eq,0:deriv(1)+1,0:deriv(2)+1,&
            &0:deriv(3)+1))
        theta_s = 0.0_dp
        ! start from theta_E
        theta_s(:,:,0,0,0) = theta_E
        theta_s(:,:,0,1,0) = 1.0_dp
        ! add the deformation described by lambda
        theta_s = theta_s + L_E(:,:,0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1)
            
        if (use_pol_flux_X) then
            ! calculate transformation matrix T_V^F
            ! (1,1)
            ierr = add_arr_mult(theta_s,q_saf_E(:,1:),&
                &T_EF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ierr = add_arr_mult(L_E(:,:,1:,0:,0:),q_saf_E,&
                &T_EF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    T_EF(id,:,1,2,deriv(1),0,0) = flux_p_E(:,deriv(1)+1)/(2*pi)
                end do
            !else
                !T_EF(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (1,3)
            T_EF(:,:,1,3,deriv(1),deriv(2),deriv(3)) = &
                &L_E(:,:,deriv(1)+1,deriv(2),deriv(3))
            ! (2,1)
            ierr = add_arr_mult(theta_s(:,:,0:,1:,0:),q_saf_E,&
                &T_EF(:,:,2,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ! (2,2)
            !T_EF(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            T_EF(:,:,2,3,deriv(1),deriv(2),deriv(3)) = &
                &theta_s(:,:,deriv(1),deriv(2)+1,deriv(3))
            ! (3,1)
            if (sum(deriv).eq.0) then
                T_EF(:,:,3,1,0,0,0) = -1.0_dp
            end if
            ierr = add_arr_mult(L_E(:,:,0:,0:,1:),q_saf_E,&
                &T_EF(:,:,3,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ! (3,2)
            !T_EF(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            T_EF(:,:,3,3,deriv(1),deriv(2),deriv(3)) = &
                &L_E(:,:,deriv(1),deriv(2),deriv(3)+1)
            
            ! determinant
            det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ierr = add_arr_mult(-theta_s(:,:,0:,1:,0:),flux_p_E(:,1:)/(2*pi),&
                &det_T_EF(:,:,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
        else
            ! set up zeta_s
            allocate(zeta_s(n_par,grp_n_r_eq,0:deriv(1)+1,0:deriv(2)+1,&
                &0:deriv(3)+1))
            zeta_s = 0.0_dp
            ! start from theta_E
            zeta_s(:,:,0,0,0) = zeta_E
            zeta_s(:,:,0,0,1) = 1.0_dp
            
            ! calculate transformation matrix T_V^F
            ! (1,1)
            T_EF(:,:,1,1,deriv(1),deriv(2),deriv(3)) = &
                &-theta_s(:,:,deriv(1)+1,deriv(2),deriv(3))
            ierr = add_arr_mult(zeta_s,rot_t_E(:,1:),&
                &T_EF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    T_EF(id,:,1,2,deriv(1),0,0) = -flux_t_E(:,deriv(1)+1)/(2*pi)
                end do
            !else
                !T_EF(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (1,3)
            !T_EF(:,:,1,3,deriv(1),deriv(2),deriv(3)) = 0
            ! (2,1)
            T_EF(:,:,2,1,deriv(1),deriv(2),deriv(3)) = &
                &-theta_s(:,:,deriv(1),deriv(2)+1,deriv(3))
            ! (2,2)
            !T_EF(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            !T_EF(:,:,2,3,deriv(1),deriv(2),deriv(3)) = 0
            ! (3,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    T_EF(id,:,3,1,deriv(1),0,0) = rot_t_E(:,deriv(1))
                end do
            end if
            T_EF(:,:,3,1,deriv(1),deriv(2),deriv(3)) = &
                &T_EF(:,:,3,1,deriv(1),deriv(2),deriv(3)) &
                &-theta_s(:,:,deriv(1),deriv(2),deriv(3)+1)
            ! (3,2)
            !T_EF(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            if (sum(deriv).eq.0) then
                T_EF(:,:,3,3,0,0,0) = -1.0_dp
            !else
                !T_EF(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            
            ! determinant
            det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ierr = add_arr_mult(theta_s(:,:,0:,1:,0:),flux_t_E(:,1:)/(2*pi),&
                &det_T_EF(:,:,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
        end if
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
    
    ! calculate the transformation matrix  between H(ELENA) and F(lux) oordinate
    ! system
    integer function calc_T_HF_ind(deriv) result(ierr)
        use num_vars, only: pi, use_pol_flux_X
        use X_vars, only: n_par
        use eq_vars, only: q_saf_E, rot_t_E, theta_E, zeta_E, flux_t_E, &
            &grp_n_r_eq
        use utilities, only: add_arr_mult
        use metric_vars, only: T_EF, det_T_EF
        
        character(*), parameter :: rout_name = 'calc_T_HF_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
            
        ! local variables
        integer :: id                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:)                             ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:)                              ! - zeta_F and derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_HF')
        CHCKERR('')
        
        ! initialize T_EF
        T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! set up theta_s
        allocate(theta_s(n_par,grp_n_r_eq,0:deriv(1)+1,0:deriv(2)+1,&
            &0:deriv(3)+1))
        theta_s = 0.0_dp
        theta_s(:,:,0,0,0) = theta_E
        theta_s(:,:,0,1,0) = 1.0_dp
            
        if (use_pol_flux_X) then
            ! calculate transformation matrix T_H^F
            ! (1,1)
            ierr = add_arr_mult(theta_s,-q_saf_E(:,1:),&
                &T_EF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ! (1,2)
            if (sum(deriv).eq.0) then
                T_EF(:,:,1,2,0,0,0) = 1.0_dp
            !else
                !T_EF(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (1,3)
            !T_EF(:,:,1,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    T_EF(id,:,2,1,deriv(1),0,0) = -q_saf_E(:,deriv(1))
                end do
            !else
                !T_EF(:,:,2,1,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (2,2)
            !T_EF(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            if (sum(deriv).eq.0) then
                T_EF(:,:,2,3,0,0,0) = 1.0_dp
            !else
                !T_EF(:,:,2,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (3,1)
            if (sum(deriv).eq.0) then
                T_EF(:,:,3,1,0,0,0) = 1.0_dp
            !else
                !T_EF(:,:,3,1,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (3,2)
            !T_EF(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            !T_EF(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            
            ! determinant
            det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            if (sum(deriv).eq.0) then
                det_T_EF(:,:,0,0,0) = 1.0_dp
            !else
                !det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
        else
            ! set up zeta_s
            allocate(zeta_s(n_par,grp_n_r_eq,0:deriv(1)+1,0:deriv(2)+1,&
                &0:deriv(3)+1))
            zeta_s = 0.0_dp
            ! start from theta_E
            zeta_s(:,:,0,0,0) = zeta_E
            zeta_s(:,:,0,0,1) = 1.0_dp
            
            ! calculate transformation matrix T_H^F
            ! (1,1)
            ierr = add_arr_mult(zeta_s,rot_t_E(:,1:),&
                &T_EF(:,:,1,1,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    T_EF(id,:,1,2,deriv(1),0,0) = flux_t_E(:,deriv(1)+1)/(2*pi)
                end do
            !else
                !T_EF(:,:,1,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (1,3)
            !T_EF(:,:,1,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,1)
            if (sum(deriv).eq.0) then
                T_EF(:,:,2,1,deriv(1),0,0) = -1.0_dp
            !else
                !T_EF(:,:,2,1,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (2,2)
            !T_EF(:,:,2,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            !T_EF(:,:,2,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    T_EF(id,:,3,1,deriv(1),0,0) = rot_t_E(:,deriv(1))
                end do
            !else
                !T_EF(:,:,3,1,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            ! (3,2)
            !T_EF(:,:,3,2,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            if (sum(deriv).eq.0) then
                T_EF(:,:,3,3,deriv(1),0,0) = 1.0_dp
            !else
                !T_EF(:,:,3,3,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
            
            ! determinant
            det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do id = 1,n_par
                    det_T_EF(id,:,deriv(1),0,0) = flux_t_E(:,deriv(1)+1)/(2*pi)
                end do
            !else
                !det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
        end if
    end function calc_T_HF_ind
    integer function calc_T_HF_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_HF_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_T_HF_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_HF_arr
    
    ! calculate D_1^m1 D_2^m2 D_3^m3 X from D_1^i1 D_2^i2 D_3^3 X and 
    ! D_1^j1 D_2^j2 D_3^j3 Y where XY = 1, i,j = 0..m, according to [ADD REF]
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_inv_met_ind(X,Y,deriv) result(ierr)                   ! matrix version
        use utilities, only: calc_inv
        use X_vars, only: n_par
        use eq_vars, only: grp_n_r_eq
        
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
                            do kd = 1,grp_n_r_eq
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
            do kd = 1,grp_n_r_eq
                do id = 1,n_par
                    X(id,kd,:,:,m1,m2,m3) = matmul(X(id,kd,:,:,m1,m2,m3),&
                        &X(id,kd,:,:,0,0,0))
                end do
            end do
        end if
    end function calc_inv_met_ind
    integer function calc_inv_met_ind_0D(X,Y,deriv) result(ierr)                ! scalar version
        use eq_vars, only: grp_n_r_eq
        use X_vars, only: n_par
        
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
            X(:,:,0,0,0) = 1._dp/Y(:,:,0,0,0)
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
                            do kd = 1,grp_n_r_eq
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
            do kd = 1,grp_n_r_eq
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
        use utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'calc_f_deriv_3_ind'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,0:,0:,0:)                             ! variable and derivs. in coord. system A
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        real(dp), intent(inout) :: X_B(1:,1:)                                   ! requested derivs. of variable in coord. system B
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
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
        ierr = check_deriv(deriv_A + sum(deriv_B),max_deriv,'calc_f_deriv')
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
                ierr = add_arr_mult(T_BA(:,:,deriv_id_B,id,0:,0:,0:),&
                    &X_B_x(:,:,0:,0:,0:),X_B,deriv_A)
                CHCKERR('')
            end do
        end if
    end function calc_f_deriv_3_ind
    integer recursive function calc_f_deriv_1_ind(X_A,T_BA,X_B,max_deriv,&
        &deriv_B,deriv_A_input) result(ierr)                                    ! flux variable version
        use utilities, only: add_arr_mult
        
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
        ierr = check_deriv([deriv_A+deriv_B,0,0],max_deriv,'calc_f_deriv')
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
            ierr = add_arr_mult(T_BA(:,0:),X_B_x(:,0:),X_B,[deriv_A,0,0])
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
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
        
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
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
        
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
    
    ! normalizes metric coefficients  and jacobian in flux  coordinates, as well
    ! as the  poloidal and toroidal flux
    subroutine normalize_metric_vars
        use eq_vars, only: R_0, B_0, psi_0
        use metric_vars, only: g_FD, h_FD, jac_FD
        
        ! local variables
        real(dp) :: g_0(3,3)                                                    ! normalization factor for the covariant metric factors
        integer :: jd, d1, d2                                                   ! counter
        
        ! set up g_0
        g_0(1,1) = R_0 * R_0
        g_0(1,2) = 1._dp/B_0
        g_0(1,3) = R_0 * R_0
        g_0(2,1) = g_0(1,2)
        g_0(2,2) = 1._dp/(R_0 * B_0) * 1/(R_0 * B_0)
        g_0(2,3) = 1._dp/B_0
        g_0(3,1) = g_0(1,3)
        g_0(3,2) = g_0(2,3)
        g_0(3,3) = R_0 * R_0
        
        ! normalize the metric coefficients
        do d2 = 1,3
            do d1 = 1,3
                g_FD(:,:,d1,d2,:,:,:) = g_FD(:,:,d1,d2,:,:,:) / g_0(d1,d2)
                h_FD(:,:,d1,d2,:,:,:) = h_FD(:,:,d1,d2,:,:,:) * g_0(d1,d2)
            end do
        end do
        do jd = 1,size(g_FD,6)-1                                                ! derivatives in alpha are scaled by A
            g_FD(:,:,:,:,:,jd,:) = g_FD(:,:,:,:,:,jd,:) * psi_0**jd 
            h_FD(:,:,:,:,:,jd,:) = h_FD(:,:,:,:,:,jd,:) * psi_0**jd
        end do
        
        ! normalize the Jacobian
        jac_FD = jac_FD * B_0/R_0
        do jd = 1,size(jac_FD,4)-1                                              ! derivatives in alpha are scaled by A
            jac_FD(:,:,:,jd,:) = jac_FD(:,:,:,jd,:) * psi_0**jd
        end do
    end subroutine normalize_metric_vars
end module metric_ops
