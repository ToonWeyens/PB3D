!------------------------------------------------------------------------------!
!   Operations that have to do with the metric elements                        !
!------------------------------------------------------------------------------!
module metric_ops 
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, max_deriv, max_str_ln, pi, grp_rank
    use utilities, only: check_deriv
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use metric_vars, only: metric_type
    
    implicit none
    private
    public calc_T_VC, calc_g_C, calc_jac_C, calc_g_V, calc_jac_V, &
        &calc_jac_H, calc_T_HF, calc_T_VF, calc_h_H, &
        &calc_inv_met, calc_g_F, calc_jac_F, transf_deriv
#if ldebug
    public test_T_EF, test_p, test_jac_F, test_g_V, test_B_F, test_D12h_H, &
        &test_jac_V, &
        &debug_calc_inv_met_ind
#endif
    
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
    interface transf_deriv
        module procedure transf_deriv_3_ind, transf_deriv_3_arr, &
            &transf_deriv_3_arr_2D, transf_deriv_1_ind
    end interface
    
    ! global variables
#if ldebug
    logical :: debug_calc_inv_met_ind = .false.                                 ! plot debug information for calc_inv_met_ind
#endif

contains
    ! calculate the lower metric elements in the C(ylindrical) coordinate system
    integer function calc_g_C_ind(eq,met,deriv) result(ierr)
        use utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_g_C_ind'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_C')
        CHCKERR('')
        
        ! initialize
        met%g_C(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! calculate
        if (sum(deriv).eq.0) then
            met%g_C(:,:,:,c([1,1],.true.),deriv(1),deriv(2),deriv(3)) = 1.0_dp
            met%g_C(:,:,:,c([3,3],.true.),deriv(1),deriv(2),deriv(3)) = 1.0_dp
        end if
        ierr = add_arr_mult(eq%R_E,eq%R_E,&
            &met%g_C(:,:,:,c([2,2],.true.),deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_g_C_ind
    integer function calc_g_C_arr(eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_C_arr'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_g_C_ind(eq,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_C_arr

    ! calculate  the metric  coefficients in  the equilibrium  V(MEC) coordinate
    ! system  using  the metric  coefficients  in  the C(ylindrical)  coordinate
    ! system and the transformation matrices
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_g_V_ind(met,deriv) result(ierr)
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_g_V_ind'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_V')
        CHCKERR('')
        
        ierr = calc_g(met%g_C,met%T_VC,met%g_E,deriv,max_deriv)
        CHCKERR('')
    end function calc_g_V_ind
    integer function calc_g_V_arr(met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_V_arr'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_g_V_ind(met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_g_V_arr
    
    ! calculate the  metric coefficients in the  equilibrium H(ELENA) coordinate
    ! system using the HELENA output
    integer function calc_h_H_ind(grid,met,deriv) result(ierr)
        use num_vars, only: max_deriv
        use HELENA, only: h_H_11, h_H_12, h_H_33
        use utilities, only: calc_deriv, c
        
        character(*), parameter :: rout_name = 'calc_h_H_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd, ld                                               ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_h_H')
        CHCKERR('')
        
        ! initialize
        met%h_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        if (deriv(3).ne.0) then                                                 ! axisymmetry: deriv. in tor. coord. is zero
            !met%h_E(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (sum(deriv).eq.0) then                                          ! no derivatives
            do jd = 1,grid%n(2)
                met%h_E(:,jd,:,c([1,1],.true.),0,0,0) = &
                    &h_H_11(:,grid%i_min:grid%i_max)
                met%h_E(:,jd,:,c([1,2],.true.),0,0,0) = &
                    &h_H_12(:,grid%i_min:grid%i_max)
                met%h_E(:,jd,:,c([3,3],.true.),0,0,0) = &
                    &h_H_33(:,grid%i_min:grid%i_max)
                met%h_E(:,jd,:,c([2,2],.true.),0,0,0) = &
                    &1._dp/h_H_11(:,grid%i_min:grid%i_max) * &
                    &(1._dp/(met%jac_E(:,jd,:,0,0,0)**2*&
                    &h_H_33(:,grid%i_min:grid%i_max)) + &
                    &h_H_12(:,grid%i_min:grid%i_max)**2)
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.0) then                        ! derivative in norm. coord.
            do ld = 1,6
                do jd = 1,grid%n(2)
                    do id = 1,grid%n(1)
                        ierr = calc_deriv(met%h_E(id,jd,:,ld,0,0,0),&
                            &met%h_E(id,jd,:,ld,1,0,0),grid%grp_r_E,1,2)        ! use extra precision to later calculate mixed derivatives
                        CHCKERR('')
                    end do
                end do
            end do
        else if (deriv(1).eq.2 .and. deriv(2).eq.0) then                        ! 2nd derivative in norm. coord.
            do ld = 1,6
                do jd = 1,grid%n(2)
                    do id = 1,grid%n(1)
                        ierr = calc_deriv(met%h_E(id,jd,:,ld,0,0,0),&
                            &met%h_E(id,jd,:,ld,2,0,0),grid%grp_r_E,2,1)
                        CHCKERR('')
                    end do
                end do
            end do
        else if (deriv(1).eq.0 .and. deriv(2).eq.1) then                        ! derivative in pol. coord.
            do ld = 1,6
                do kd = 1,grid%grp_n_r
                    do jd = 1,grid%n(2)
                        ierr = calc_deriv(met%h_E(:,jd,kd,ld,0,0,0),&
                            &met%h_E(:,jd,kd,ld,0,1,0),&
                            &grid%theta_E(:,jd,kd),1,2)                         ! use extra precision to later calculate mixed derivatives
                        CHCKERR('')
                    end do
                end do
            end do
        else if (deriv(1).eq.0 .and. deriv(2).eq.2) then                        ! 2nd derivative in pol. coord.
            do ld = 1,6
                do kd = 1,grid%grp_n_r
                    do jd = 1,grid%n(2)
                        ierr = calc_deriv(met%h_E(:,jd,kd,ld,0,0,0),&
                            &met%h_E(:,jd,kd,ld,0,2,0),&
                            &grid%theta_E(:,jd,kd),2,1)
                        CHCKERR('')
                    end do
                end do
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.1) then                        ! mixed derivative in norm. and pol. coord.
            do ld = 1,6
                do jd = 1,grid%n(2)
                    do id = 1,grid%n(1)
                        ierr = calc_deriv(met%h_E(id,jd,:,ld,0,1,0),&
                            &met%h_E(id,jd,:,ld,1,1,0),grid%grp_r_E,1,1)
                        CHCKERR('')
                    end do
                end do
            end do
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_h_H_ind
    integer function calc_h_H_arr(grid,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_h_H_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_h_H_ind(grid,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_h_H_arr

    ! calculate the  metric coefficients in  the F(lux) coordinate  system using
    ! the  metric coefficients  in  the equilibrium  coordinate  system and  the
    ! trans- formation matrices
    integer function calc_g_F_ind(met,deriv) result(ierr)
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_g_F_ind'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g_F')
        CHCKERR('')
        
        ! calculate g_F
        ierr = calc_g(met%g_E,met%T_FE,met%g_F,deriv,max_deriv)
        CHCKERR('')
    end function calc_g_F_ind
    integer function calc_g_F_arr(met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_g_F_arr'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_g_F_ind(met,deriv(:,id))
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
        use utilities, only: calc_mult, conv_sym
        
        character(*), parameter :: rout_name = 'calc_g'
        
        ! input / output
        real(dp), intent(in) :: g_A(:,:,:,:,0:,0:,0:)
        real(dp), intent(in) :: T_BA(:,:,:,:,0:,0:,0:)
        real(dp), intent(inout) :: g_B(:,:,:,:,0:,0:,0:)
        integer, intent(in) :: deriv(3)
        integer, intent(in) :: max_deriv
        
        ! local variables
        integer :: i1, j1, i2, j2, i3, j3, k1, k2, k3                           ! counters
        integer :: m1, m2, m3                                                   ! alias for deriv
        real(dp), allocatable :: C1(:), C2(:), C3(:)                            ! coeff. for (il)
        real(dp), allocatable :: dum1(:,:,:,:)                                  ! dummy variable representing metric matrix for some derivatives
        real(dp), allocatable :: dum2(:,:,:,:)                                  ! dummy variable representing metric matrix for some derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_g')
        CHCKERR('')
        
        ! set ml
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! set up dummies
        allocate(dum1(size(g_A,1),size(g_A,2),size(g_A,3),9))
        allocate(dum2(size(g_A,1),size(g_A,2),size(g_A,3),9))
        
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
                                ! set up ki
                                k1 = m1 - j1 - i1
                                k2 = m2 - j2 - i2
                                k3 = m3 - j3 - i3
                                ! calculate term to add to the summation
                                ierr = calc_mult(g_A(:,:,:,:,j1,j2,j3),&
                                    &T_BA(:,:,:,:,k1,k2,k3),dum1,3,&
                                    &transp=[.false.,.true.])
                                CHCKERR('')
                                ierr = calc_mult(T_BA(:,:,:,:,i1,i2,i3),&
                                    &dum1,dum2,3)
                                CHCKERR('')
                                ! convert result to lower-diagonal storage
                                ierr = conv_sym(dum2,dum1(:,:,:,1:6),3)
                                CHCKERR('')
                                g_B(:,:,:,:,m1,m2,m3) = g_B(:,:,:,:,m1,m2,m3) &
                                    &+C1(i1)*C2(i2)*C3(i3)*dum1(:,:,:,1:6)
                            end do
                        end do d3
                    end do
                end do d2
            end do
        end do d1
        
        ! deallocate local variables
        deallocate(dum1,dum2)
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
    integer function calc_jac_C_ind(eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_C_ind'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_C')
        CHCKERR('')
        
        met%jac_C(:,:,:,deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3))
    end function calc_jac_C_ind
    integer function calc_jac_C_arr(eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_C_arr'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_C_ind(eq,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_C_arr
    
    ! calculate the jacobian in equilibrium V(MEC) coordinates from 
    !   jac_V = det(T_VC) jac_C
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_V_ind(met,deriv) result(ierr)
        use utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_V_ind'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_V')
        CHCKERR('')
        
        ! initialize
        met%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate determinant
        ierr = add_arr_mult(met%jac_C,met%det_T_VC,&
            &met%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_jac_V_ind
    integer function calc_jac_V_arr(met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_V_arr'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_V_ind(met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_V_arr
    
    ! calculate the jacobian in HELENA coordinates directly from J = q R^2/F
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_H_ind(grid,eq,met,deriv) result(ierr)
        use HELENA, only:  h_H_33, RBphi
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'calc_jac_H_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_H')
        CHCKERR('')
        
        ! calculate determinant
        if (deriv(3).ne.0) then                                                 ! axisymmetry: deriv. in tor. coord. is zero
            met%jac_E(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        else if (sum(deriv).eq.0) then                                          ! no derivatives
            do kd = 1,grid%grp_n_r
                do jd = 1,grid%n(2)
                    met%jac_E(:,jd,kd,0,0,0) = eq%q_saf_E(kd,0)/&
                        &(h_H_33(:,grid%i_min-1+kd)*RBphi(grid%i_min-1+kd))
                end do
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.0) then                        ! derivative in norm. coord.
            do jd = 1,grid%n(2)
                do id = 1,grid%n(1)
                    ierr = calc_deriv(met%jac_E(id,jd,:,0,0,0),&
                        &met%jac_E(id,jd,:,1,0,0),grid%grp_r_E,1,2)             ! use extra precision to later calculate mixed derivatives
                    CHCKERR('')
                end do
            end do
        else if (deriv(1).eq.2 .and. deriv(2).eq.0) then                        ! 2nd derivative in norm. coord.
            do jd = 1,grid%n(2)
                do id = 1,grid%n(1)
                    ierr = calc_deriv(met%jac_E(id,jd,:,0,0,0),&
                        &met%jac_E(id,jd,:,2,0,0),grid%grp_r_E,2,1)
                    CHCKERR('')
                end do
            end do
        else if (deriv(1).eq.0 .and. deriv(2).eq.1) then                        ! derivative in pol. coord.
            do kd = 1,grid%grp_n_r
                do jd = 1,grid%n(2)
                    ierr = calc_deriv(met%jac_E(:,jd,kd,0,0,0),&
                        &met%jac_E(:,jd,kd,0,1,0),grid%theta_E(:,jd,kd),1,2)    ! use extra precision to later calculate mixed derivatives
                    CHCKERR('')
                end do
            end do
        else if (deriv(1).eq.0 .and. deriv(2).eq.2) then                        ! 2nd derivative in pol. coord.
            do kd = 1,grid%grp_n_r
                do jd = 1,grid%n(2)
                    ierr = calc_deriv(met%jac_E(:,jd,kd,0,0,0),&
                        &met%jac_E(:,jd,kd,0,2,0),grid%theta_E(:,jd,kd),2,1)
                    CHCKERR('')
                end do
            end do
        else if (deriv(1).eq.1 .and. deriv(2).eq.1) then                        ! mixed derivative in norm. and pol. coord.
            do kd = 1,grid%grp_n_r
                do jd = 1,grid%n(2)
                    ierr = calc_deriv(met%jac_E(:,jd,kd,1,0,0),&
                        &met%jac_E(:,jd,kd,1,1,0),grid%theta_E(:,jd,kd),1,1)
                    CHCKERR('')
                end do
            end do
        else
            ierr = 1
            err_msg = 'Derivative of order ('//trim(i2str(deriv(1)))//','//&
                &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                &) not supported'
            CHCKERR(err_msg)
        end if
    end function calc_jac_H_ind
    integer function calc_jac_H_arr(grid,eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_H_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_H_ind(grid,eq,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_H_arr
    
    ! calculate the jacobian in Flux coordinates from 
    !   jac_F = det(T_FE) jac_E
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_jac_F_ind(met,deriv) result(ierr)
        use utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'calc_jac_F_ind'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_J_F')
        CHCKERR('')
        
        ! initialize
        met%jac_F(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        
        ! calculate determinant
        ierr = add_arr_mult(met%jac_E,met%det_T_FE,&
            &met%jac_F(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_jac_F_ind
    integer function calc_jac_F_arr(met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_jac_F_arr'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_jac_F_ind(met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_jac_F_arr

    ! calculate  the transformation  matrix  between  C(ylindrical) and  V(mec)
    ! coordinate system
    integer function calc_T_VC_ind(eq,met,deriv) result(ierr)
        use utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_VC_ind'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VC')
        CHCKERR('')
        
        ! initialize
        met%T_VC(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! calculate transformation matrix T_V^C
        met%T_VC(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1)+1,deriv(2),deriv(3))
        !met%T_VC(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) = 0
        met%T_VC(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%Z_E(:,:,:,deriv(1)+1,deriv(2),deriv(3))
        met%T_VC(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1),deriv(2)+1,deriv(3))
        !met%T_VC(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = 0
        met%T_VC(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%Z_E(:,:,:,deriv(1),deriv(2)+1,deriv(3))
        met%T_VC(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%R_E(:,:,:,deriv(1),deriv(2),deriv(3)+1)
        if (sum(deriv).eq.0) then
            met%T_VC(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 1.0_dp
        else
            !met%T_VC(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 0
        end if
        met%T_VC(:,:,:,c([3,3],.false.),deriv(1),deriv(2),deriv(3)) = &
            &eq%Z_E(:,:,:,deriv(1),deriv(2),deriv(3)+1)
        
        ! determinant
        met%det_T_VC(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
        ierr = add_arr_mult(eq%R_E(:,:,:,0:,1:,0:),eq%Z_E(:,:,:,1:,0:,0:),&
            &met%det_T_VC(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
        ierr = add_arr_mult(-eq%R_E(:,:,:,1:,0:,0:),eq%Z_E(:,:,:,0:,1:,0:),&
            &met%det_T_VC(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
        CHCKERR('')
    end function calc_T_VC_ind
    integer function calc_T_VC_arr(eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VC_arr'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VC_ind(eq,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VC_arr

    ! calculate the transformation matrix  between equilibrium V(mec) and F(lux)
    ! oordinate system
    integer function calc_T_VF_ind(grid,eq,met,deriv) result(ierr)
        use num_vars, only: use_pol_flux_F
        use utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_VF_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:,:)                           ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:,:)                            ! - zeta_F and derivatives
        integer :: dims(3)                                                      ! dimensions
        integer :: c1                                                           ! 2D coordinate in metric_type storage convention
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_VF')
        CHCKERR('')
        
        ! set up dims
        dims = [grid%n(1),grid%n(2),grid%grp_n_r]
        
        ! initialize
        met%T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! set up theta_s
        allocate(theta_s(dims(1),dims(2),dims(3),&
            &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
        theta_s = 0.0_dp
        ! start from theta_E
        theta_s(:,:,:,0,0,0) = grid%theta_E
        theta_s(:,:,:,0,1,0) = 1.0_dp
        ! add the deformation described by lambda
        theta_s = theta_s + eq%L_E(:,:,:,0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1)
            
        if (use_pol_flux_F) then
            ! calculate transformation matrix T_V^F
            ! (1,1)
            ierr = add_arr_mult(theta_s,eq%q_saf_E(:,1:),&
                &met%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ierr = add_arr_mult(eq%L_E(:,:,:,1:,0:,0:),eq%q_saf_E(:,0:),&
                &met%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%T_EF(:,:,kd,c([1,2],.false.),deriv(1),0,0) = &
                        &eq%flux_p_E(kd,deriv(1)+1)/(2*pi)
                end do
            !else
                !met%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (1,3)
            met%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                &eq%L_E(:,:,:,deriv(1)+1,deriv(2),deriv(3))
            ! (2,1)
            ierr = add_arr_mult(theta_s(:,:,:,0:,1:,0:),eq%q_saf_E,&
                &met%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (2,2)
            !met%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            met%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                &theta_s(:,:,:,deriv(1),deriv(2)+1,deriv(3))
            ! (3,1)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([3,1],.false.),0,0,0) = -1.0_dp
            end if
            ierr = add_arr_mult(eq%L_E(:,:,:,0:,0:,1:),eq%q_saf_E,&
                &met%T_EF(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (3,2)
            !met%T_EF(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            met%T_EF(:,:,:,c([3,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                &eq%L_E(:,:,:,deriv(1),deriv(2),deriv(3)+1)
            
            ! determinant
            met%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ierr = add_arr_mult(-theta_s(:,:,:,0:,1:,0:),&
                &eq%flux_p_E(:,1:)/(2*pi),&
                &met%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
        else
            ! set up zeta_s
            allocate(zeta_s(dims(1),dims(2),dims(3),0:deriv(1)+1,&
                &0:deriv(2)+1,0:deriv(3)+1))
            zeta_s = 0.0_dp
            ! start from zeta_E
            zeta_s(:,:,:,0,0,0) = grid%zeta_E
            zeta_s(:,:,:,0,0,1) = 1.0_dp
            
            ! calculate transformation matrix T_V^F
            ! (1,1)
            met%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                &-theta_s(:,:,:,deriv(1)+1,deriv(2),deriv(3))
            ierr = add_arr_mult(zeta_s,eq%rot_t_E(:,1:),&
                &met%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%T_EF(:,:,kd,c([1,2],.false.),deriv(1),0,0) = &
                        &-eq%flux_t_E(kd,deriv(1)+1)/(2*pi)
                end do
            !else
                !met%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (1,3)
            !met%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,1)
            met%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                &-theta_s(:,:,:,deriv(1),deriv(2)+1,deriv(3))
            ! (2,2)
            !met%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            !met%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%T_EF(:,:,kd,c([3,1],.false.),deriv(1),0,0) = &
                        &eq%rot_t_E(kd,deriv(1))
                end do
            end if
            c1 = c([3,1],.false.)                                               ! to avoid compiler crash (as of 26-02-2015)
            met%T_EF(:,:,:,c1,deriv(1),deriv(2),deriv(3)) = &
                &met%T_EF(:,:,:,c1,deriv(1),deriv(2),deriv(3)) &
                &-theta_s(:,:,:,deriv(1),deriv(2),deriv(3)+1)
            ! (3,2)
            !met%T_EF(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([3,3],.false.),0,0,0) = -1.0_dp
            !else
                !met%T_EF(:,:,:,c([3,3],.false.)deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            
            ! determinant
            met%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ierr = add_arr_mult(theta_s(:,:,:,0:,1:,0:),&
                &eq%flux_t_E(:,1:)/(2*pi),&
                &met%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)),deriv)
            CHCKERR('')
        end if
    end function calc_T_VF_ind
    integer function calc_T_VF_arr(grid,eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_VF_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_T_VF_ind(grid,eq,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_VF_arr
    
    ! calculate the transformation matrix  between H(ELENA) and F(lux) oordinate
    ! system
    integer function calc_T_HF_ind(grid,eq,met,deriv) result(ierr)
        use num_vars, only: use_pol_flux_F
        use utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'calc_T_HF_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:)
            
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: theta_s(:,:,:,:,:,:)                           ! theta_F and derivatives
        real(dp), allocatable :: zeta_s(:,:,:,:,:,:)                            ! - zeta_F and derivatives
        integer :: dims(3)                                                      ! dimensions
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_T_HF')
        CHCKERR('')
        
        ! set up dims
        dims = [grid%n(1),grid%n(2),grid%grp_n_r]
        
        ! initialize
        met%T_EF(:,:,:,:,deriv(1),deriv(2),deriv(3)) = 0._dp
        
        ! set up theta_s
        allocate(theta_s(dims(1),dims(2),dims(3),&
            &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
        theta_s = 0.0_dp
        theta_s(:,:,:,0,0,0) = grid%theta_E
        theta_s(:,:,:,0,1,0) = 1.0_dp
            
        if (use_pol_flux_F) then
            ! calculate transformation matrix T_H^F
            ! (1,1)
            ierr = add_arr_mult(theta_s,-eq%q_saf_E(:,1:),&
                &met%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([1,2],.false.),0,0,0) = 1._dp
            !else
                !met%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (1,3)
            !met%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%T_EF(:,:,kd,c([2,1],.false.),deriv(1),0,0) = &
                        &-eq%q_saf_E(kd,deriv(1))
                end do
            !else
                !met%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (2,2)
            !met%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([2,3],.false.),0,0,0) = 1.0_dp
            !else
                !met%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (3,1)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([3,1],.false.),0,0,0) = 1.0_dp
            !else
                !met%T_EF(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (3,2)
            !met%T_EF(:,:,:,c([3,2],.false),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            !met%T_EF(:,:,:,c([3,3],.false),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            
            ! determinant
            met%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                met%det_T_EF(:,:,:,deriv(1),0,0) = 1._dp
            !else
                !met%det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
        else
            ! set up zeta_s
            allocate(zeta_s(dims(1),dims(2),dims(3),&
                &0:deriv(1)+1,0:deriv(2)+1,0:deriv(3)+1))
            zeta_s = 0.0_dp
            ! start from zeta_E
            zeta_s(:,:,:,0,0,0) = grid%zeta_E
            zeta_s(:,:,:,0,0,1) = 1.0_dp
            
            ! calculate transformation matrix T_H^F
            ! (1,1)
            ierr = add_arr_mult(zeta_s,eq%rot_t_E(:,1:),&
                &met%T_EF(:,:,:,c([1,1],.false.),deriv(1),deriv(2),deriv(3)),&
                &deriv)
            CHCKERR('')
            ! (1,2)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%T_EF(:,:,kd,c([1,2],.false.),deriv(1),0,0) = &
                        &eq%q_saf_E(kd,deriv(1))
                end do
            !else
                !met%T_EF(:,:,:,c([1,2],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (1,3)
            !met%T_EF(:,:,:,c([1,3],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,1)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([2,1],.false.),0,0,0) = -1.0_dp
            !else
                !met%T_EF(:,:,:,c([2,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (2,2)
            !met%T_EF(:,:,:,c([2,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (2,3)
            !met%T_EF(:,:,:,c([2,3],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,1)
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%T_EF(:,:,kd,c([3,1],.false.),deriv(1),0,0) = &
                        &eq%rot_t_E(kd,deriv(1))
                end do
            !else
                !met%T_EF(:,:,:,c([3,1],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            ! (3,2)
            !met%T_EF(:,:,:,c([3,2],.false.),deriv(1),deriv(2),deriv(3)) = 0.0_dp
            ! (3,3)
            if (sum(deriv).eq.0) then
                met%T_EF(:,:,:,c([3,3],.false.),deriv(1),0,0) = 1.0_dp
            !else
                !met%T_EF(:,:,:,c([3,3],.false.),deriv(1),deriv(2),deriv(3)) = &
                    !&0.0_dp
            end if
            
            ! determinant
            met%det_T_EF(:,:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            if (deriv(2).eq.0 .and. deriv(3).eq.0) then
                do kd = 1,dims(3)
                    met%det_T_EF(:,:,kd,deriv(1),0,0) = eq%q_saf_E(kd,deriv(1))
                end do
            !else
                !met%det_T_EF(:,:,deriv(1),deriv(2),deriv(3)) = 0.0_dp
            end if
        end if
    end function calc_T_HF_ind
    integer function calc_T_HF_arr(grid,eq,met,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_T_HF_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_T_HF_ind(grid,eq,met,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_T_HF_arr
    
    ! calculate D_1^m1 D_2^m2 D_3^m3 X from D_1^i1 D_2^i2 D_3^3 X and 
    ! D_1^j1 D_2^j2 D_3^j3 Y where XY = 1, i,j = 0..m, according to [ADD REF]
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_inv_met_ind(X,Y,deriv) result(ierr)                   ! matrix version
        use utilities, only: calc_inv, calc_mult, c, conv_sym
#if ldebug
        use num_vars, only: grp_n_procs
#endif
        
        character(*), parameter :: rout_name = 'calc_inv_met_ind'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)                      ! X
        real(dp), intent(in) :: Y(1:,1:,1:,1:,0:,0:,0:)                         ! Y
        integer, intent(in) :: deriv(:)
        
        ! local variables
        integer :: r, t, z                                                      ! counters for derivatives
        integer :: m1, m2, m3                                                   ! alias for deriv(i)
        integer :: dims(3)                                                      ! dimensions of coordinates
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        real(dp), allocatable :: dum1(:,:,:,:)                                  ! dummy variable representing metric matrix for some derivatives
        real(dp), allocatable :: dum2(:,:,:,:)                                  ! dummy variable representing metric matrix for some derivatives
        integer :: nn                                                           ! size of matrices
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        real(dp), allocatable :: XY(:,:,:,:)                                    ! to check if XY is indeed unity
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set nn
        nn = size(X,4)
        
        ! tests
        if (size(Y,4).ne.nn) then
            ierr = 1
            err_msg = 'Both matrices X and Y have to be either in full or &
                &lower diagonal storage'
            CHCKERR(err_msg)
        end if
        
        ! set mi
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! set dims
        dims = [size(X,1),size(X,2),size(X,3)]
        
        ! initialize the requested derivative to zero
        X(:,:,:,:,m1,m2,m3) = 0._dp
        
        ! set up dummy variables
        allocate(dum1(dims(1),dims(2),dims(3),9))
        allocate(dum2(dims(1),dims(2),dims(3),9))
        
        ! initialize dum2
        dum2 = 0._dp
        
        ! calculate terms
        if (m1.eq.0 .and. m2.eq.0 .and. m3.eq.0) then                           ! direct inverse
            ierr = calc_inv(X(:,:,:,:,0,0,0),Y(:,:,:,:,0,0,0),3)
            CHCKERR('')
            
#if ldebug
            ! check if X*Y is unity indeed if debugging
            if (debug_calc_inv_met_ind) then
                if (grp_n_procs.gt.1) call writo('WARNING: Debugging of &
                    &calc_inv_met_ind should be done with one processor only')
                call writo('checking if X*Y = 1')
                call lvl_ud(1)
                allocate(XY(dims(1),dims(2),dims(3),9))
                ierr = calc_mult(X(:,:,:,:,0,0,0),Y(:,:,:,:,0,0,0),XY,3)
                CHCKERR('')
                do r = 1,3
                    do t = 1,3
                        call writo(&
                            &trim(r2strt(minval(XY(:,1,:,c([r,t],.false.)))))//&
                            &' < element ('//trim(i2str(t))//','//&
                            &trim(i2str(r))//') < '//&
                            &trim(r2strt(maxval(XY(:,1,:,c([r,t],.false.))))))
                    end do
                end do
                deallocate(XY)
                call lvl_ud(-1)
            end if
#endif
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
                            ! multiply X(r,t,z) and Y(m1-r,m2-t,m3-z)
                            ierr = calc_mult(X(:,:,:,:,r,t,z),&
                                &Y(:,:,:,:,m1-r,m2-t,m3-z),dum1,3)
                            CHCKERR('')
                            ! add result, multiplied by bin_fac, to dum2
                            dum2 = dum2 - product(bin_fac)*dum1
                        end if
                    end do
                end do
            end do
            
            ! right-multiply by X(0,0,0) and save in dum1
            ierr = calc_mult(dum2,X(:,:,:,:,0,0,0),dum1,3)
            CHCKERR('')
            ! save result in X(m1,m2,m3), converting if necessary
            ierr = conv_sym(dum1,X(:,:,:,:,m1,m2,m3),3)
            CHCKERR('')
            
            ! deallocate local variables
            deallocate(dum1,dum2)
        end if
    end function calc_inv_met_ind
    integer function calc_inv_met_ind_0D(X,Y,deriv) result(ierr)                ! scalar version
        character(*), parameter :: rout_name = 'calc_inv_met_ind_0D'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:)
        
        ! local variables
        integer :: r, t, z                                                      ! counters for derivatives
        integer :: m1, m2, m3                                                   ! alias for deriv(i)
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        
        ! initialize ierr
        ierr = 0
        
        m1 = deriv(1)
        m2 = deriv(2)
        m3 = deriv(3)
        
        ! initialize the requested derivative
        X(:,:,:,m1,m2,m3) = 0.0_dp
        
        ! calculate terms
        if (m1.eq.0 .and. m2.eq.0 .and. m3.eq.0) then                           ! direct inverse
            X(:,:,:,0,0,0) = 1._dp/Y(:,:,:,0,0,0)
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
                            X(:,:,:,m1,m2,m3) = X(:,:,:,m1,m2,m3)-&
                                &bin_fac(1)*bin_fac(2)*bin_fac(3)* &
                                &X(:,:,:,r,t,z)*Y(:,:,:,m1-r,m2-t,m3-z)
                        end if
                    end do
                end do
            end do
            
            ! right-multiply by X(0,0,0)
            X(:,:,:,m1,m2,m3) = X(:,:,:,m1,m2,m3)*X(:,:,:,0,0,0)
        end if
    end function calc_inv_met_ind_0D
    integer function calc_inv_met_arr(X,Y,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_met_arr'
        
        ! input / output
        real(dp), intent(inout) :: X(1:,1:,1:,1:,0:,0:,0:)                      ! X
        real(dp), intent(in) :: Y(1:,1:,1:,1:,0:,0:,0:)                         ! Y
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
        real(dp), intent(inout) :: X(1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: Y(1:,1:,1:,0:,0:,0:)
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_inv_met_ind_0D(X,Y,deriv(:,id))
            CHCKERR('')
        end do
    end function calc_inv_met_arr_0D
    
    ! Calculates  derivatives in  a coordinate  system B  from derivatives  in a
    ! coordinates system A, making use of the transformation matrix T_BA.
    ! The routine works  by exchanging the derivatives in the  coordinates B for
    ! derivatives in coordinates A using the formula
    !   D^m_B X = T_BA D_A (D^m-1_B X)
    ! This is done for the derivatives in each of the coordinates B until degree
    ! 0  is  reached. Furthermore,  each  of  these  degrees of  derivatives  in
    ! coordinates  B has to be derived optionally, also in the original A, which
    ! coordinates yields the formula:
    !   D^p_A (D^m_B X) = sum_q binom(p,q) D^q_A (T_BA) D^p-q_A D_A (D^m-1_B X)
    ! This way, ultimately  the desired derivatives in the coordinates  B can be
    ! obtained recursively from the lower orders in the coordinates B and higher
    ! orders in the coordinates A
    ! see [ADD REF] for more detailed information
    integer recursive function transf_deriv_3_ind(X_A,T_BA,X_B,max_deriv,&
        &deriv_B,deriv_A_input) result(ierr)                                    ! normal variable version
        use utilities, only: add_arr_mult, c
        
        character(*), parameter :: rout_name = 'transf_deriv_3_ind'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,1:,0:,0:,0:)                          ! variable and derivs. in coord. system A
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        real(dp), intent(inout) :: X_B(1:,1:,1:)                                ! requested derivs. of variable in coord. system B
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
        integer, intent(in) :: deriv_B(:)                                       ! derivs. in coord. system B
        integer, intent(in), optional :: deriv_A_input(:)                       ! derivs. in coord. system A (optional)
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: deriv_id_B                                                   ! holds the deriv. currently calculated
        integer, allocatable :: deriv_A(:)                                      ! holds either deriv_A or 0
        real(dp), allocatable :: X_B_x(:,:,:,:,:,:)                             ! X_B and derivs. D^p-q_A with extra 1 exchanged deriv. D_A
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
        ierr = check_deriv(deriv_A + sum(deriv_B),max_deriv,'transf_deriv')
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
            X_B = X_A(:,:,:,deriv_A(1),deriv_A(2),deriv_A(3))
        else                                                                    ! apply the formula to obtain X_B using exchanged derivs.
            ! allocate  the  exchanged  X_B_x  and  the  derivs.  in  A  and  B,
            ! corresponding to D^m-1_B X
            allocate(X_B_x(size(X_B,1),size(X_B,2),size(X_B,3),&
                &0:deriv_A(1),0:deriv_A(2),0:deriv_A(3)))
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
                            ierr = transf_deriv_3_ind(X_A,T_BA,&
                                &X_B_x(:,:,:,jd,kd,ld),max_deriv,deriv_B_x,&
                                &deriv_A_x)
                            CHCKERR('')
                        end do
                    end do
                end do
                
                ! use X_B_x at this term in summation due to the transf. mat. to
                ! update X_B using the formula
                ierr = add_arr_mult(&
                    &T_BA(:,:,:,c([deriv_id_B,id],.false.),0:,0:,0:),&
                    &X_B_x(:,:,:,0:,0:,0:),X_B,deriv_A)
                CHCKERR('')
            end do
        end if
    end function transf_deriv_3_ind
    integer function transf_deriv_3_arr(X_A,T_BA,X_B,max_deriv,derivs) &
        &result(ierr)                                                           ! matrix version
        character(*), parameter :: rout_name = 'transf_deriv_3_arr'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,1:,0:,0:,0:)                             ! variable and derivs. in coord. system A
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        real(dp), intent(inout) :: X_B(1:,1:,1:,0:,0:,0:)                          ! requested derivs. of variable in coord. system B
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
        integer, intent(in) :: derivs(:,:)                                      ! series of derivs. (in coordinate system B)
        
        ! local variables
        integer :: id                                                           ! counter
        
        do id = 1, size(derivs,2)
            ierr = transf_deriv_3_ind(X_A,T_BA,&
                &X_B(:,:,:,derivs(1,id),derivs(2,id),derivs(3,id)),&
                &max_deriv,derivs(:,id))
            CHCKERR('')
        end do
    end function transf_deriv_3_arr
    integer function transf_deriv_3_arr_2D(X_A,T_BA,X_B,max_deriv,derivs) &
        &result(ierr)                                                           ! matrix version
        use utilities, only: is_sym, c
        
        character(*), parameter :: rout_name = 'transf_deriv_3_arr_2D'
        
        ! input / output
        real(dp), intent(in) :: X_A(1:,1:,1:,1:,0:,0:,0:)                       ! variable and derivs. in coord. system A
        real(dp), intent(in) :: T_BA(1:,1:,1:,1:,0:,0:,0:)                      ! transf. mat. and derivs. between coord. systems A and B
        real(dp), intent(inout) :: X_B(1:,1:,1:,1:,0:,0:,0:)                    ! requested derivs. of variable in coord. system B
        integer, intent(in) :: max_deriv                                        ! maximum degrees of derivs.
        integer, intent(in) :: derivs(:,:)                                      ! series of derivs. (in coordinate system B)
        
        ! local variables
        integer :: id, kd                                                       ! counters
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: sym                                                          ! sym
        
        ! test whether the dimensions of X_A and X_B agree
        if (size(X_A,4).ne.size(X_B,4)) then
            err_msg = 'X_A and X_B need to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ierr = is_sym(3,size(X_A,4),sym)
        CHCKERR('')
        
        do kd = 1,size(X_A,4)
            do id = 1, size(derivs,2)
                ierr = transf_deriv_3_ind(X_A(:,:,:,kd,:,:,:),T_BA,&
                    &X_B(:,:,:,kd,derivs(1,id),derivs(2,id),derivs(3,id)),&
                    &max_deriv,derivs(:,id))
                CHCKERR('')
            end do
        end do
    end function transf_deriv_3_arr_2D
    integer recursive function transf_deriv_1_ind(X_A,T_BA,X_B,max_deriv,&
        &deriv_B,deriv_A_input) result(ierr)                                    ! flux variable version
        use utilities, only: add_arr_mult
        
        character(*), parameter :: rout_name = 'transf_deriv_1_ind'
        
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
        ierr = check_deriv([deriv_A+deriv_B,0,0],max_deriv,'transf_deriv')
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
                ierr = transf_deriv_1_ind(X_A,T_BA,X_B_x(:,jd),max_deriv,&
                    &deriv_B_x,deriv_A_input=deriv_A_x)
                CHCKERR('')
            end do
            
            ! use X_B_x to calculate X_B using the formula
            ierr = add_arr_mult(T_BA(:,0:),X_B_x(:,0:),X_B,[deriv_A,0,0])
            CHCKERR('')
        end if
    end function transf_deriv_1_ind
    
#if ldebug
    ! See if T_EF it complies with the theory of [ADD REF]
    integer function test_T_EF(grid_eq,eq,met) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style
        use grid_vars, only: destroy_grid
        use grid_ops, only: trim_grid
        use utilities, only: c
        use output_ops, only: plot_diff_HDF5
        
        character(*), parameter :: rout_name = 'test_T_EF'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(metric_type), intent(in) :: met                                    ! metric variables
        
        ! local variables
        integer :: id, kd                                                       ! counter
        real(dp), allocatable :: res(:,:,:,:)                                   ! calculated result
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether T_EF complies with the theory')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r,9))
        res = 0.0_dp
        
        ! calc  res,  depending  on flux  used  in F  coordinates  and on  which
        ! equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                if (use_pol_flux_F) then                                        ! using poloidal flux
                    ! calculate T_EF(1,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,1) = &
                            &grid_trim%theta_F(:,:,kd)*eq%q_saf_E(kd,1) + &
                            &eq%L_E(:,:,kd,1,0,0)*eq%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(2,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,2) = &
                            &(1._dp + eq%L_E(:,:,kd,0,1,0))*eq%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(3,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,3) = -1._dp + &
                            &eq%L_E(:,:,kd,0,0,1)*eq%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(1,2)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,4) = eq%flux_p_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(1,3)
                    res(:,:,:,7) = eq%L_E(:,:,1:grid_trim%grp_n_r,1,0,0)
                    ! calculate T_EF(2,3)
                    res(:,:,:,8) = 1._dp + &
                        &eq%L_E(:,:,1:grid_trim%grp_n_r,0,1,0)
                    ! calculate T_EF(3,3)
                    res(:,:,:,9) = eq%L_E(:,:,1:grid_trim%grp_n_r,0,0,1)
                else                                                            ! using toroidal flux
                    ! calculate T_EF(1,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,1) = - eq%L_E(:,:,kd,1,0,0) + &
                            &grid_trim%zeta_E(:,:,kd)*eq%rot_t_E(kd,1)
                    end do
                    ! calculate T_EF(2,1)
                    res(:,:,:,2) = &
                        &- (1._dp + eq%L_E(:,:,1:grid_trim%grp_n_r,0,1,0))
                    ! calculate T_EF(3,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,3) = eq%rot_t_E(kd,0) - eq%L_E(:,:,kd,0,0,1)
                    end do
                    ! calculate T_EF(1,2)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,4) = - eq%flux_t_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(3,3)
                    res(:,:,:,9) = -1._dp
                end if
            case (2)                                                            ! HELENA
                if (use_pol_flux_F) then                                        ! using poloidal flux
                    ! calculate T_EF(1,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,1) = &
                            &- eq%q_saf_E(kd,1)*grid_trim%theta_E(:,:,kd)
                    end do
                    ! calculate T_EF(2,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,2) = - eq%q_saf_E(kd,0)
                    end do
                    ! calculate T_EF(3,1)
                    res(:,:,:,3) = 1._dp
                    ! calculate T_EF(1,2)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,4) = eq%flux_p_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(2,3)
                    res(:,:,:,8) = 1._dp
                else                                                            ! using toroidal flux
                    ! calculate T_EF(1,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,1) = &
                            &eq%rot_t_E(kd,1)*grid_trim%zeta_E(:,:,kd)
                    end do
                    ! calculate T_EF(2,1)
                    res(:,:,:,2) = -1._dp
                    ! calculate T_EF(3,1)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,3) = eq%rot_t_E(kd,0)
                    end do
                    ! calculate T_EF(1,2)
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,4) = eq%flux_t_E(kd,1)/(2*pi)
                    end do
                    ! calculate T_EF(3,3)
                    res(:,:,:,9) = 1._dp
                end if
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! set up plot variables for calculated values
        do id = 1,3
            do kd = 1,3
                ! user output
                call writo('Testing T_EF('//trim(i2str(kd))//','//&
                    &trim(i2str(id))//')')
                call lvl_ud(1)
                
                ! set some variables
                file_name = 'TEST_T_EF_'//trim(i2str(kd))//'_'//trim(i2str(id))
                description = 'Testing calculated with given value for T_EF('//&
                    &trim(i2str(kd))//','//trim(i2str(id))//')'
                
                ! plot difference
                call plot_diff_HDF5(res(:,:,:,c([kd,id],.false.)),&
                    &met%T_EF(:,:,1:grid_trim%grp_n_r,c([kd,id],.false.),&
                    &0,0,0),file_name,tot_dim,grp_offset,description,&
                    &output_message=.true.)
                
                call lvl_ud(-1)
            end do
        end do
        
        ! clean up
        call destroy_grid(grid_trim)
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_T_EF
    
    ! performs tests on pressure balance
    !   - mu_0 D2p = 1/J (D3 B_2 - D2 B_3)
    !   - mu_0 J D3p = 0 => D3 B_1 = D1 B_3
    ! (working in the (modified) Flux coordinates (alpha,psi,theta)_F)
    integer function test_p(grid_eq,eq,met) result(ierr)
        use utilities, only: c, calc_deriv
        use grid_ops, only: trim_grid
        use eq_vars, only: vac_perm
        use num_vars, only: eq_style
        use HELENA, only: RBphi, h_H_11, h_H_12
        
        character(*), parameter :: rout_name = 'test_p'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        
        ! local variables
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        integer :: id, jd, kd                                                   ! counters
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        integer :: n_H(2)                                                       ! indices in Helena coords.
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test the consistency of the metric variables with &
            &the given pressure')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r,2))
        
        ! user output
        call writo('Checking if mu_0 D2 p = 1/J (D3 B_2 - D2_B3)')
        call lvl_ud(1)
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! calculate mu_0 D2p
        res(:,:,:,1) = &
            &(met%g_FD(:,:,1:grid_trim%grp_n_r,c([2,3],.true.),0,0,1) - &
            &met%g_FD(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),0,1,0))/&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) - &
            &(met%g_FD(:,:,1:grid_trim%grp_n_r,c([2,3],.true.),0,0,0)*&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,1) - &
            &met%g_FD(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),0,0,0)*&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,1,0))/ &
            &(met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0)**2)
        res(:,:,:,1) = res(:,:,:,1)/met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0)
        
        ! save mu_0 D2p in res
        do kd = 1,grid_trim%grp_n_r
            res(:,:,kd,2) = vac_perm*eq%pres_FD(kd,1)
        end do
        
        ! set some variables
        file_name = 'TEST_D2p'
        description = 'Testing whether mu_0 D2 p = 1/J (D3 B_2 - D2_B3)'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
            &grp_offset,description,output_message=.true.)
        
        call lvl_ud(-1)
        
        ! user output
        call writo('Checking if mu_0 J D1p = 0 => D3 B_1 = D2 B_3')
        call lvl_ud(1)
        
        ! calculate D3 B1
        res(:,:,:,1) = met%g_FD(:,:,1:grid_trim%grp_n_r,c([1,3],.true.),0,0,1)/&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) - &
            &met%g_FD(:,:,1:grid_trim%grp_n_r,c([1,3],.true.),0,0,0)*&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,1) / &
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0)**2
        
        ! calculate D1 B3
        res(:,:,:,2) = met%g_FD(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),1,0,0)/&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) - &
            &met%g_FD(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),0,0,0)*&
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,1,0,0) / &
            &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0)**2
        
        ! set some variables
        file_name = 'TEST_D1p'
        description = 'Testing whether D3 B_1 = D1 B_3'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
            &grp_offset,description,output_message=.true.)
        
        call lvl_ud(-1)
        
        ! extra testing if Helena
        if (eq_style.eq.2) then
            ! set n_H
            n_H = [grid_eq%i_min,grid_eq%i_min+grid_trim%grp_n_r-1]
            
            ! user output
            call writo('Checking if D2 B_3 |F = D1 (qF+qh_11/F) |H')
            call lvl_ud(1)
            
            ! calculate if D2 B_3 = D1 (qF) + D1 (q/F h_11)
            res(:,:,:,1) = 1._dp/met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) * &
                &(met%g_FD(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),0,1,0)/&
                &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) - &
                &met%g_FD(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),0,0,0)*&
                &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,1,0)/ &
                &(met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0)**2))
            do jd = 1,grid_trim%n(2)
                do id = 1,grid_trim%n(1)
                    ierr = calc_deriv(eq%q_saf_E(:,0)*RBphi(n_H(1):n_H(2))+&
                        &eq%q_saf_E(:,0)*h_H_11(id,n_H(1):n_H(2))/&
                        &RBphi(n_H(1):n_H(2)),res(id,jd,:,2),&
                        &grid_trim%grp_r_E,1,1)
                    CHCKERR('')
                end do
            end do
            
            ! set some variables
            file_name = 'TEST_D2B_3'
            description = 'Testing whether D2 B_3 |F = D1 (qF+qh_11/F) |H'
            
            ! plot difference
            call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
                &grp_offset,description,output_message=.true.)
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Checking if D3 B_2 |F = D2 -(qh_11/F) + q''F |H')
            call lvl_ud(1)
            
            ! calculate if D2 B_3 = D1 (q/F h_11) + q'F
            res(:,:,:,1) = 1._dp/met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) * &
                &(met%g_FD(:,:,1:grid_trim%grp_n_r,c([2,3],.true.),0,0,1)/&
                &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0) - &
                &met%g_FD(:,:,1:grid_trim%grp_n_r,c([2,3],.true.),0,0,0)*&
                &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,1)/ &
                &(met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0)**2))
            do jd = 1,grid_trim%n(2)
                do id = 1,grid_trim%n(1)
                    ierr = calc_deriv(eq%q_saf_E(:,0)*h_H_12(id,n_H(1):n_H(2))/&
                        &RBphi(n_H(1):n_H(2))+&
                        &eq%q_saf_E(:,1)*RBphi(n_H(1):n_H(2)),res(id,jd,:,2),&
                        &grid_trim%grp_r_E,1,1)
                    CHCKERR('')
                end do
            end do
            
            ! set some variables
            file_name = 'TEST_D3B_2'
            description = 'Testing whether D3 B_2 |F = -D2 (qh_12/F) +q''F |H'
            
            ! plot difference
            call plot_diff_HDF5(res(:,:,:,1),res(:,:,:,2),file_name,tot_dim,&
                &grp_offset,description,output_message=.true.)
            
            call lvl_ud(-1)
            
            call print_HDF5('B_2_TEST','B_2_TEST',&
                &met%g_FD(:,:,1:grid_trim%grp_n_r,c([2,3],.true.),0,0,0)/&
                &met%jac_FD(:,:,1:grid_trim%grp_n_r,0,0,0),tot_dim=tot_dim,&
                &grp_offset=grp_offset)
            
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_p
    
    ! performs tests on Jacobian in Flux coordinates
    !   - comparing it with the determinant of g_F
    !   - comparing it with the direct formula
    integer function test_jac_F(grid_eq,eq,met) result(ierr)
        use num_vars, only: eq_style, use_pol_flux_F
        use grid_ops, only: trim_grid
        use utilities, only: calc_det
        use HELENA, only: h_H_33, RBphi
        
        character(*), parameter :: rout_name = 'test_jac_F'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        
        ! local variables
        real(dp), allocatable :: res(:,:,:)                                     ! result variable
        integer :: jd, kd                                                       ! counters
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), pointer :: Dflux(:)                                           ! points to D flux_p or D flux_t in E coordinates
        integer :: pmone                                                        ! plus or minus one
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test the calculation of the Jacobian in Flux &
            &coordinates')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r))
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! 1. Compare with determinant of g_F
        
        ! calculate Jacobian from determinant of g_F
        ierr = calc_det(res,met%g_F(:,:,1:grid_trim%grp_n_r,:,0,0,0),3)
        CHCKERR('')
        
        ! set some variables
        file_name = 'TEST_jac_F_1'
        description = 'Testing whether the Jacobian in Flux coordinates is &
            &consistent with determinant of metric matrix'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &met%jac_F(:,:,1:grid_trim%grp_n_r,0,0,0)**2,&
            &file_name,tot_dim,grp_offset,description,output_message=.true.)
        
        ! 2. Compare with explicit formula
        
        ! set up Dflux and pmone
        if (use_pol_flux_F) then                                                ! using poloidal flux
            Dflux => eq%flux_p_E(:,1)
            pmone = 1                                                           ! flux_p_V = flux_p_F
        else
            Dflux => eq%flux_t_E(:,1)
            pmone = -1                                                          ! flux_t_V = - flux_t_F
        end if
        
        ! calculate jac_F
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                do kd = 1,grid_trim%grp_n_r
                    res(:,:,kd) = pmone * 2*pi * eq%R_E(:,:,kd,0,0,0)*&
                        &(eq%R_E(:,:,kd,1,0,0)*eq%Z_E(:,:,kd,0,1,0) - &
                        &eq%Z_E(:,:,kd,1,0,0)*eq%R_E(:,:,kd,0,1,0)) / &
                        &(Dflux(kd)*(1+eq%L_E(:,:,kd,0,1,0)))
                end do
            case (2)                                                            ! HELENA
                do kd = 1,grid_trim%grp_n_r
                    do jd = 1,grid_trim%n(2)
                        res(:,jd,kd) = eq%q_saf_E(kd,0)/&
                            &(h_H_33(:,grid_eq%i_min-1+kd)*&
                            &RBphi(grid_eq%i_min-1+kd))                         ! h_H_33 = 1/R^2 and RBphi = F are tabulated in eq. grid
                    end do
                end do
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! set some variables
        file_name = 'TEST_jac_F_2'
        description = 'Testing whether the Jacobian in Flux coordinates is &
            &consistent with explicit formula'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &met%jac_F(:,:,1:grid_trim%grp_n_r,0,0,0),file_name,tot_dim,&
            &grp_offset,description,output_message=.true.)
        
        ! clean up
        nullify(Dflux)
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_jac_F
    
    ! Tests whether g_V is calculated correctly
    integer function test_g_V(grid_eq,eq,met) result(ierr)
        use grid_ops, only: trim_grid
        use utilities, only: c
        
        character(*), parameter :: rout_name = 'test_g_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether g_V is calculated correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r,6))
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! calculate g_V(1,1)
        res(:,:,:,1) = eq%R_E(:,:,1:grid_trim%grp_n_r,1,0,0)**2 + &
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,1,0,0)**2
        ! calculate g_V(2,1)
        res(:,:,:,2) = eq%R_E(:,:,1:grid_trim%grp_n_r,1,0,0)*&
            &eq%R_E(:,:,1:grid_trim%grp_n_r,0,1,0) + &
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,1,0,0)*&
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,0,1,0)
        ! calculate g_V(3,1)
        res(:,:,:,3) = eq%R_E(:,:,1:grid_trim%grp_n_r,1,0,0)*&
            &eq%R_E(:,:,1:grid_trim%grp_n_r,0,0,1) + &
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,1,0,0)*&
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,0,0,1)
        ! calculate g_V(2,2)
        res(:,:,:,4) = eq%R_E(:,:,1:grid_trim%grp_n_r,0,1,0)**2 + &
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,0,1,0)**2
        ! calculate g_V(3,2)
        res(:,:,:,5) = eq%R_E(:,:,1:grid_trim%grp_n_r,0,1,0)*&
            &eq%R_E(:,:,1:grid_trim%grp_n_r,0,0,1) + &
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,0,1,0)*&
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,0,0,1)
        ! calculate g_V(3,3)
        res(:,:,:,6) = eq%R_E(:,:,1:grid_trim%grp_n_r,0,0,1)**2 + &
            &eq%Z_E(:,:,1:grid_trim%grp_n_r,0,0,1)**2 + &
            &eq%R_E(:,:,1:grid_trim%grp_n_r,0,0,0)**2
        
        ! set up plot variables for calculated values
        do id = 1,3
            do kd = 1,3
                ! user output
                call writo('Testing g_V('//trim(i2str(kd))//','//&
                    &trim(i2str(id))//')')
                call lvl_ud(1)
                
                ! set some variables
                file_name = 'TEST_g_V_'//trim(i2str(kd))//'_'//trim(i2str(id))
                description = 'Testing calculated with given value for g_V('//&
                    &trim(i2str(kd))//','//trim(i2str(id))//')'
                
                ! plot difference
                call plot_diff_HDF5(res(:,:,:,c([kd,id],.true.)),&
                    &met%g_E(:,:,1:grid_trim%grp_n_r,c([kd,id],.true.),&
                    &0,0,0),file_name,tot_dim,grp_offset,description,&
                    &output_message=.true.)
                
                call lvl_ud(-1)
            end do
        end do
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_g_V
    
    ! Tests whether jac_V is calculated correctly
    integer function test_jac_V(grid_eq,met) result(ierr)
        use grid_ops, only: trim_grid
        use fourier_ops, only: fourier2real
        use VMEC, only: jac_V_c, jac_V_s
        
        character(*), parameter :: rout_name = 'test_jac_V'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        
        ! local variables
        real(dp), allocatable :: res(:,:,:)                                     ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether the Jacobian in the VMEC &
            &coords. jac_V is calculated correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r))
        
        ! get jac_V from VMEC
        ierr = fourier2real(jac_V_c(:,:,grid_trim%i_min:grid_trim%i_max),&
            &jac_V_s(:,:,grid_trim%i_min:grid_trim%i_max),&
            &grid_eq%trigon_factors(:,:,:,:,1:grid_trim%grp_n_r,:),&
            &res,[0,0])
        CHCKERR('')
        write(*,*) 'done calculating fourier'
        
        ! user output
        call writo('Testing jac_V')
        call lvl_ud(1)
        
        ! set some variables
        file_name = 'TEST_jac_V'
        description = 'Testing calculated with given value for jac_V'
        
        ! plot difference
        call plot_diff_HDF5(res(:,:,:),&
            &met%jac_E(:,:,1:grid_trim%grp_n_r,0,0,0),file_name,&
            &tot_dim,grp_offset,description,output_message=.true.)
        
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_jac_V
    
    ! Tests whether B_F is calculated correctly
    integer function test_B_F(grid_eq,eq,met) result(ierr)
        use num_vars, only: eq_style
        use grid_ops, only: trim_grid
        use utilities, only: c
        use fourier_ops, only: fourier2real
        use VMEC, only: B_V_sub_s, B_V_sub_c, B_V_c, B_V_s
        
        character(*), parameter :: rout_name = 'test_B_F'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        real(dp), allocatable :: res2(:,:,:,:)                                  ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether the magnetic field is calculated &
            &correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set up res and res2
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r,4))
        allocate(res2(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r,4))
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! get covariant components and magnitude of B_E
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                do id = 1,3
                    ierr = fourier2real(&
                        &B_V_sub_c(:,:,grid_trim%i_min:grid_trim%i_max,id),&
                        &B_V_sub_s(:,:,grid_trim%i_min:grid_trim%i_max,id),&
                        &grid_eq%trigon_factors(:,:,:,:,1:grid_trim%grp_n_r,:),&
                        &res(:,:,:,id),[0,0])
                    CHCKERR('')
                end do
                ierr = fourier2real(B_V_c(:,:,grid_trim%i_min:grid_trim%i_max),&
                    &B_V_s(:,:,grid_trim%i_min:grid_trim%i_max),&
                    &grid_eq%trigon_factors(:,:,:,:,1:grid_trim%grp_n_r,:),&
                    &res(:,:,:,4),[0,0])
                CHCKERR('')
            case (2)                                                            ! HELENA
                res(:,:,:,1) = &
                    &met%g_E(:,:,1:grid_trim%grp_n_r,c([1,2],.true.),0,0,0)
                res(:,:,:,2) = &
                    &met%g_E(:,:,1:grid_trim%grp_n_r,c([2,2],.true.),0,0,0)
                do kd = 1,grid_trim%grp_n_r
                    res(:,:,kd,3) = eq%q_saf_E(kd,0)*&
                        &met%g_E(:,:,kd,c([3,3],.true.),0,0,0)
                    res(:,:,kd,4) = &
                        &sqrt(met%g_E(:,:,kd,c([2,2],.true.),0,0,0) + &
                        &eq%q_saf_E(kd,0)**2*&
                        &met%g_E(:,:,kd,c([3,3],.true.),0,0,0))
                end do
                do id = 1,4
                    do kd = 1,grid_trim%grp_n_r
                        res(:,:,kd,id) = res(:,:,kd,id)/met%jac_E(:,:,kd,0,0,0)
                    end do
                end do
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! transform them to Flux coord. system
        res2 = 0._dp
        do id = 1,3
            do kd = 1,3
                res2(:,:,:,id) = res(:,:,:,kd) * &
                    &met%T_FE(:,:,1:grid_trim%grp_n_r,c([id,kd],.false.),0,0,0)
            end do
        end do
        res2(:,:,:,4) = res(:,:,:,4)
        
        ! set up plot variables for calculated values
        do id = 1,4
            ! user output
            if (id.lt.4) then
                call writo('Testing B_F_'//trim(i2str(id)))
            else
                call writo('Testing B_F')
            end if
            call lvl_ud(1)
            
            ! set some variables
            if (id.lt.4) then
                file_name = 'TEST_B_F_'//trim(i2str(id))
                description = 'Testing calculated with given value for B_F_'//&
                    &trim(i2str(id))
            else
                file_name = 'TEST_B_F'
                description = 'Testing calculated with given value for B_F'
            end if
            
            ! plot difference
            if (id.lt.4) then
                call plot_diff_HDF5(res2(:,:,:,id),&
                    &met%g_F(:,:,1:grid_trim%grp_n_r,c([3,id],.true.),0,0,0)/&
                    &met%jac_F(:,:,1:grid_trim%grp_n_r,0,0,0),file_name,&
                    &tot_dim,grp_offset,description,output_message=.true.)
            else
                call plot_diff_HDF5(res2(:,:,:,id),sqrt(&
                    &met%g_F(:,:,1:grid_trim%grp_n_r,c([3,3],.true.),0,0,0))/&
                    &met%jac_F(:,:,1:grid_trim%grp_n_r,0,0,0),file_name,&
                    &tot_dim,grp_offset,description,output_message=.true.)
            end if
            
            call lvl_ud(-1)
        end do
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_B_F
    
    ! Tests whether D1 D2 h_H is calculated correctly
    integer function test_D12h_H(grid_eq,met) result(ierr)
        use grid_ops, only: trim_grid
        use utilities, only: c, calc_deriv
        
        character(*), parameter :: rout_name = 'test_D12h_H'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        real(dp), allocatable :: res(:,:,:,:)                                   ! result variable
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: tot_dim(3), grp_offset(3)                                    ! total dimensions and group offset
        type(grid_type) :: grid_trim                                            ! trimmed equilibrium grid
        
        ! initialize ierr
        ierr = 0
        
        ! output
        call writo('Going to test whether D1 D2 h_H is calculated correctly')
        call lvl_ud(1)
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_eq,grid_trim)
        CHCKERR('')
        
        ! set up res
        allocate(res(grid_trim%n(1),grid_trim%n(2),grid_trim%grp_n_r,6))
        
        ! set total and group dimensions and group offset
        tot_dim = [grid_trim%n(1),grid_trim%n(2),grid_trim%n(3)]
        grp_offset = [0,0,grid_trim%i_min-1]
        
        ! calculate D1 D2 h_H alternatively
        do ld = 1,6
            do kd = 1,grid_trim%grp_n_r
                do jd = 1,grid_trim%n(2)
                    ierr = calc_deriv(met%h_E(:,jd,kd,ld,1,0,0),&
                        &res(:,jd,kd,ld),grid_trim%theta_E(:,jd,kd),1,1)
                    CHCKERR('')
                end do
            end do
        end do
        
        ! set up plot variables for calculated values
        do id = 1,3
            do kd = 1,3
                ! user output
                call writo('Testing h_H('//trim(i2str(kd))//','//&
                    &trim(i2str(id))//')')
                call lvl_ud(1)
                
                ! set some variables
                file_name = 'TEST_D12h_H_'//trim(i2str(kd))//'_'//&
                    &trim(i2str(id))
                description = 'Testing calculated with given value for D12h_H('&
                    &//trim(i2str(kd))//','//trim(i2str(id))//')'
                
                ! plot difference
                call plot_diff_HDF5(res(:,:,:,c([kd,id],.true.)),&
                    &met%h_E(:,:,grid_trim%i_min:grid_trim%i_max,&
                    &c([kd,id],.true.),1,1,0),file_name,tot_dim,grp_offset,&
                    &description,output_message=.true.)
                
                call lvl_ud(-1)
            end do
        end do
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_D12h_H
#endif
end module metric_ops
