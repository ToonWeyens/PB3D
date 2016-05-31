!------------------------------------------------------------------------------!
!   Numerical utilities related to equilibrium variables                       !
!------------------------------------------------------------------------------!
module eq_utilities
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln, max_deriv
    use grid_vars, only: grid_type, disc_type
    use eq_vars, only: eq_1_type, eq_2_type
    use num_utilities, only: check_deriv
    
    implicit none
    private
    public calc_inv_met, calc_g, transf_deriv, calc_F_derivs, divide_eq_jobs, &
        &do_eq, eq_info, print_info_eq
#if ldebug
    public debug_calc_inv_met_ind
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_inv_met_ind = .false.                                 ! plot debug information for calc_inv_met_ind
#endif
    
    ! interfaces
    interface calc_F_derivs
        module procedure calc_F_derivs_1, calc_F_derivs_2
    end interface
    interface calc_inv_met
        module procedure calc_inv_met_ind, calc_inv_met_arr, &
            &calc_inv_met_ind_0D, calc_inv_met_arr_0D
    end interface
    interface transf_deriv
        module procedure transf_deriv_3_ind, transf_deriv_3_arr, &
            &transf_deriv_3_arr_2D, transf_deriv_1_ind
    end interface
    
contains
    ! calculate D_1^m1 D_2^m2 D_3^m3 X from D_1^i1 D_2^i2 D_3^3 X and 
    ! D_1^j1 D_2^j2 D_3^j3 Y where XY = 1, i,j = 0..m, according to [ADD REF]
    ! NOTE: It is assumed that the  lower order derivatives have been calculated
    !       already. If not, the results will be incorrect!
    integer function calc_inv_met_ind(X,Y,deriv) result(ierr)                   ! matrix version
        use num_utilities, only: calc_inv, calc_mult, c, conv_mat
#if ldebug
        use num_vars, only: n_procs
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
                if (n_procs.gt.1) call writo('Debugging of calc_inv_met_ind &
                    &should be done with one processor only',warning=.true.)
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
            ierr = conv_mat(dum1,X(:,:,:,:,m1,m2,m3),3)
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
        use num_utilities, only: calc_mult, conv_mat
        
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
                                ierr = conv_mat(dum2,dum1(:,:,:,1:6),3)
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
        use num_utilities, only: add_arr_mult, c
        
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
        if (deriv_id_B.eq.0) then                                               ! return the function with its required A derivs.
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
        use num_utilities, only: is_sym, c
        
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
        use num_utilities, only: add_arr_mult
        
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
        if (deriv_B.eq.0) then                                                  ! return the function with its required A derivs.
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
    
    ! Transforms  derivatives of  the  equilibrium quantities in  E
    ! coordinates to derivatives in the F coordinates.
    integer function calc_F_derivs_1(grid_eq,eq) result(ierr)                   ! flux version
        use num_vars, only: eq_style, max_deriv, use_pol_flux_F
        use num_utilities, only: derivs, c, fac
        use eq_vars, only: max_flux_E
        
        character(*), parameter :: rout_name = 'calc_F_derivs_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_1_type), intent(inout) :: eq                                    ! flux equilibrium variables
        
        ! local variables
        integer :: id
        real(dp), allocatable :: T_FE_loc(:,:)                                  ! T_FE(2,1)
        
        ! initialize ierr
        ierr = 0
        
        ! Transform depending on equilibrium style being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! user output
                call writo('Transform flux equilibrium quantities to flux &
                    &coordinates')
                
                call lvl_ud(1)
                
                ! set up local T_FE
                allocate(T_FE_loc(grid_eq%loc_n_r,0:max_deriv+1))
                if (use_pol_flux_F) then
                    do id = 0,max_deriv+1
                        T_FE_loc(:,id) = 2*pi/max_flux_E*eq%q_saf_E(:,id)
                    end do
                else
                    T_FE_loc(:,0) = -2*pi/max_flux_E
                    T_FE_loc(:,1:max_deriv+1) = 0._dp
                end if
                
                ! Transform the  derivatives in E coordinates  to derivatives in
                ! the F coordinates.
                do id = 0,max_deriv+1
                    ! pres_FD
                    ierr = transf_deriv(eq%pres_E,T_FE_loc,eq%pres_FD(:,id),&
                        &max_deriv+1,id)
                    CHCKERR('')
                    
                    ! flux_p_FD
                    ierr = transf_deriv(eq%flux_p_E,T_FE_loc,&
                        &eq%flux_p_FD(:,id),max_deriv+1,id)
                    CHCKERR('')
                    
                    ! flux_t_FD
                    ierr = transf_deriv(eq%flux_t_E,T_FE_loc,&
                        &eq%flux_t_FD(:,id),max_deriv+1,id)
                    CHCKERR('')
                    eq%flux_t_FD(:,id) = -eq%flux_t_FD(:,id)                    ! multiply by plus minus one
                    
                    ! q_saf_FD
                    ierr = transf_deriv(eq%q_saf_E,T_FE_loc,eq%q_saf_FD(:,id),&
                        &max_deriv+1,id)
                    CHCKERR('')
                    eq%q_saf_FD(:,id) = -eq%q_saf_FD(:,id)                      ! multiply by plus minus one
                    
                    ! rot_t_FD
                    ierr = transf_deriv(eq%rot_t_E,T_FE_loc,eq%rot_t_FD(:,id),&
                        &max_deriv+1,id)
                    CHCKERR('')
                    eq%rot_t_FD(:,id) = -eq%rot_t_FD(:,id)                      ! multiply by plus minus one
                end do
                
                call lvl_ud(-1)
            case (2)                                                            ! HELENA
                ! E derivatives are equal to F derivatives
                eq%flux_p_FD = eq%flux_p_E
                eq%flux_t_FD = eq%flux_t_E
                eq%pres_FD = eq%pres_E
                eq%q_saf_FD = eq%q_saf_E
                eq%rot_t_FD = eq%rot_t_E
        end select
    end function calc_F_derivs_1
    integer function calc_F_derivs_2(eq) result(ierr)                           ! metric version
        use num_vars, only: eq_style
        use num_utilities, only: derivs, c
        
        character(*), parameter :: rout_name = 'calc_F_derivs_2'
        
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! metric equilibrium variables
        
        ! local variables
        integer :: id
        integer :: pmone                                                        ! plus or minus one
        
        ! initialize ierr
        ierr = 0
        
        ! Set up pmone, depending on equilibrium style being used
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                pmone = -1                                                      ! conversion VMEC LH -> RH coord. system
            case (2)                                                            ! HELENA
                pmone = 1
        end select
        
        ! user output
        call writo('Transform metric equilibrium quantities to flux &
            &coordinates')
        
        call lvl_ud(1)
        
        ! Transform the  derivatives in E coordinates  to derivatives in        
        ! the F coordinates.
        do id = 0,max_deriv
            ! g_FD
            ierr = transf_deriv(eq%g_F,eq%T_FE,eq%g_FD,max_deriv,derivs(id))
            CHCKERR('')
            
            ! h_FD
            ierr = transf_deriv(eq%h_F,eq%T_FE,eq%h_FD,max_deriv,derivs(id))
            CHCKERR('')
            
            ! jac_FD
            ierr = transf_deriv(eq%jac_F,eq%T_FE,eq%jac_FD,max_deriv,&
                &derivs(id))
            CHCKERR('')
        end do
        
        call lvl_ud(-1)
    end function calc_F_derivs_2
    
    ! Divides the equilibrium jobs.
    ! The entire parallel  range has to be calculated, but  due to memory limits
    ! this has to be split up in pieces. In its most extreme form, this would be
    ! the individual  calculation on a  fundamental integration integral  of the
    ! parallel points:
    !   - for magn_int_style = 1 (trapezoidal), this is 1 point,
    !   - for magn_int_style = 2 (Simpson 3/8), this is 3 points.
    ! For HELENA, as the parallel derivatives are calculated discretely, and the
    ! calculations for  the first Richardson  level happen on the  HELENA output
    ! grid, followed  by an interpolation  that depends on the  Richardson level
    ! and an  integral of the interpolated  X_2 values, the splitting  up of the
    ! parallel range must  be understood in terms of splitting  up the different
    ! parts of the  interpolated and integrated X_2 variables, though  it has no
    ! direct influence on the equilibrium variables.
    ! Therefore, the whole  load is divided into jobs depending  on the sizes of
    ! the sizes of the pieces in  memory. The division of equilibrium jobs forms
    ! an  "outer  loop"  compared  to  the  "inner  loop"  of  the  division  of
    ! perturbation jobs.
    ! This  function does  the  job of  dividing the  grids  setting the  global
    ! variables 'eq_jobs_lims'.  In contrast  with "divide_X_jobs", there  is no
    ! need for a variable such as "eq_jobs_taken". See note below.
    ! Optionally, a base n_par_X can be set, which is undivisible. 
    ! Also, optionally the total memory size per process is also returned.
    ! Note: a  very important difference  between the equilibrium  parallel jobs
    ! and the perturbation jobs is that  the parallel jobs are done serially, by
    ! all processes. The perturbations jobs are individual asynchronous jobs for
    ! each of the processes. Furthermore, perturbation jobs only comprise either
    ! vectorial or  tensorial perturbation jobs  while parallel jobs  comprise a
    ! whole equilibrium  and perturbation calculation, where  the integration of
    ! the tensorial perturbation variables is possibly adjusted:
    !   - If  the first job  of the parallel jobs  and not the  first Richardson
    !   level: add half  of the integrated tensorial  perturbation quantities of
    !   the previous level.
    !   -  If  not the  first  job  of the  parallel  jobs,  add the  integrated
    !   tensorial perturbation quantities to those of the previous parallel job,
    !   same Richardson level.
    ! In fact,  the equilibrium  jobs have  much in  common with  the Richardson
    ! levels,  as is  attested  by the  existence of  the  routines "do_eq"  and
    ! "eq_info", which are equivalent to "do_rich" and "rich_info".
    integer function divide_eq_jobs(n_par_X_rich,arr_size,n_par_X_base,&
        &tot_mem_size) result(ierr)
        
        use num_vars, only: max_tot_mem_per_proc, max_X_mem_per_proc, &
            &eq_jobs_lims, mem_scale_fac, eq_job_nr, magn_int_style
        use files_utilities, only: nextunit
        use rich_vars, only: rich_lvl
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'divide_eq_jobs'
        
        ! input / output
        integer, intent(in) :: n_par_X_rich                                     ! number of parallel points in this Richardson level
        integer, intent(in) :: arr_size                                         ! array size (using loc_n_r)
        integer, intent(in), optional :: n_par_X_base                           ! base n_par_X, undivisible
        real(dp), intent(inout), optional :: tot_mem_size                       ! total memory size
        
        ! local variables
        real(dp) :: mem_size                                                    ! approximation of memory required for X variables
        integer :: n_div                                                        ! factor by which to divide the total size
        integer :: n_div_max                                                    ! maximum n_div
        integer :: n_par_range                                                  ! nr. of points in range
        integer :: fund_n_par                                                   ! fundamental interval width
        integer, allocatable :: n_par_loc(:)                                    ! number of points per range
        character(len=max_str_ln) :: range_message                              ! message about how many ranges
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Dividing the equilibrium jobs')
        call lvl_ud(1)
        
        ! set up width of fundamental interval
        select case (magn_int_style)
            case (1)                                                            ! Trapezoidal rule
                fund_n_par = 1   
            case (2)                                                            ! Simpson's 3/8 rule
                fund_n_par = 3   
        end select
        
        ! calculate largest possible range of parallel points
        n_div = 0
        mem_size = huge(1._dp)
        if (rich_lvl.eq.1) then
            n_div_max = (n_par_X_rich-1)/fund_n_par                             ! through adapt_min_n_par_X: n_par_X_rich = 1 + fund_n_par(1+k)
        else
            n_div_max = n_par_X_rich/fund_n_par                                 ! next levels: n_par_X_rich - 1 = fund_n_par(1+k)
        end if
        do while (mem_size.gt.max_tot_mem_per_proc)
            n_div = n_div + 1
            n_par_range = ceiling(n_par_X_rich*1._dp/n_div)
            if (present(n_par_X_base)) n_par_range = n_par_range + n_par_X_base
            ierr = calc_memory(arr_size,n_par_range,mem_size)
            CHCKERR('')
            if (n_div.eq.n_div_max) then                                        ! reached the maximum n_div
                if (mem_size.gt.max_tot_mem_per_proc) then                      ! still not enough memory
                    ierr = 1
                    err_msg = 'The memory limit is too low, need more than '//&
                        &trim(i2str(ceiling(mem_size)))//'MB'
                    CHCKERR(err_msg)
                end if
            end if
        end do
        if (n_div.gt.1) then
            range_message = 'The '//trim(i2str(n_par_X_rich))//' parallel &
                &points are split into '//trim(i2str(n_div))//' and '//&
                &trim(i2str(n_div))//' collective jobs are done serially'
        else
            range_message = 'The '//trim(i2str(n_par_X_rich))//' parallel &
                &points can be done without splitting them'
        end if
        call writo(range_message)
        call writo('The memory per process is estimated to be about '//&
            &trim(i2str(ceiling(mem_size)))//'MB (maximum: '//&
            &trim(i2str(ceiling(max_tot_mem_per_proc)))//'MB)')
        
        ! calculate max memory available for perturbation calculations
        max_X_mem_per_proc = max_tot_mem_per_proc - mem_size/mem_scale_fac      ! don't need scale factor as no operations are done on eq vars
        call writo('In the perturbation phase, the equilibrium variables are &
            &not being operated on:')
        call lvl_ud(1)
        call writo('This translates to a scale factor 1/'//&
            &trim(r2strt(mem_scale_fac)))
        call writo('Therefore, the memory left for the perturbation phase is '&
            &//trim(i2str(ceiling(max_X_mem_per_proc)))//'MB')
        call lvl_ud(-1)
        
        ! set total memory size if requested
        if (present(tot_mem_size)) then
            n_par_range = n_par_X_rich
            if (present(n_par_X_base)) n_par_range = n_par_range + n_par_X_base
            ierr = calc_memory(arr_size,n_par_range,tot_mem_size)
            CHCKERR('')
            tot_mem_size = 1._dp*ceiling(tot_mem_size)                          ! round up
        end if
        
        ! set up jobs data as illustrated below for 3 divisions:
        !   [1,2,3]
        ! Also initialize eq_job_nr to 0 as it is incremented by "do_eq".
        allocate(n_par_loc(n_div))
        n_par_loc = n_par_X_rich/n_div                                          ! number of radial points on this processor
        n_par_loc(1:mod(n_par_X_rich,n_div)) = &
            &n_par_loc(1:mod(n_par_X_rich,n_div)) + 1                           ! add a point to if there is a remainder
        ierr = calc_eq_jobs_lims(n_par_loc,eq_jobs_lims)
        CHCKERR('')
        eq_job_nr = 0
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Equilibrium jobs divided')
    contains
        ! Calculate memory in MB necessary for eq variables:
        !   (2*6+1) x n_geo x loc_n_r x n_par x (max_deriv + 1)^3
        ! where  n_geo x  loc_n_r should  be passed as  'arr_size' and  n_par as
        ! well.
        ! Note:  only g_FD,  h_FD and  jac_FD  are counted,  as the  equilibrium
        ! variables and the transformation matrices are deleted after use. Also,
        ! S, sigma, kappa_n and kappa_g can  be neglected as they do not contain
        ! derivatives.
        integer function calc_memory(arr_size,n_par,mem_size) result(ierr)
            use ISO_C_BINDING
            use num_vars, only: max_deriv
            
            character(*), parameter :: rout_name = 'calc_memory'
            
            ! input / output
            integer, intent(in) :: arr_size                                     ! size of part of X array
            integer, intent(in) :: n_par                                        ! number of parallel points
            real(dp), intent(inout) :: mem_size                                 ! total size
            
            ! local variables
            integer(C_SIZE_T) :: dp_size                                        ! size of dp
            
            ! initialize ierr
            ierr = 0
            
            call lvl_ud(1)
            
            ! get size of real variable
            dp_size = sizeof(1._dp)
            
            ! set memory size
            mem_size = 13*arr_size
            mem_size = mem_size*n_par
            mem_size = mem_size*(max_deriv+1)**3
            mem_size = mem_size*dp_size
            
            ! convert B to MB
            mem_size = mem_size*1.E-6_dp
            
            ! scale memory to account for rough estimation
            mem_size = mem_size*mem_scale_fac
            
            ! test overflow
            if (mem_size.lt.0) then
                ierr = 1
                CHCKERR('Overflow occured')
            end if
            
            call lvl_ud(-1)
        end function calc_memory
        
        ! Calculate eq_jobs_lims.
        ! Take into account that every job has to start and end at the start and
        ! end of a fundamental integration interval, as discussed in above:
        !   - for magn_int_style = 1 (trapezoidal), this is 1 point,
        !   - for magn_int_style = 2 (Simpson 3/8), this is 3 points.
        integer function calc_eq_jobs_lims(n_par,res) result(ierr)
            character(*), parameter :: rout_name = 'calc_eq_jobs_lims'
            
            ! input / output
            integer, intent(in) :: n_par(:)                                     ! eq jobs data
            integer, allocatable :: res(:,:)                                    ! result
            
            ! local variables
            integer :: n_div                                                    ! nr. of divisions of parallel ranges
            integer :: id                                                       ! counter
            integer :: ol_width                                                 ! overlap width
            
            ! initialize ierr
            ierr = 0
            
            ! set up nr. of divisions
            n_div = size(n_par)
            
            ! (re)allocate result
            if (allocated(res)) deallocate(res)
            allocate(res(2,n_div))
            
            ! set overlap width
            if (rich_lvl.eq.1) then
                ol_width = 1
            else
                ol_width = 0
            end if
            
            ! loop over divisions
            do id = 1,n_div
                ! setup first guess, without taking into account fundamental
                ! integration intervals
                if (id.eq.1) then
                    res(1,id) = 1
                else
                    res(1,id) = res(2,id-1) + 1 - ol_width
                end if
                res(2,id) = sum(n_par(1:id))
                
                ! take into account fundamental interval
                res(2,id) = res(1,id) - 1 + ol_width + fund_n_par * max(1,&
                    &nint((res(2,id)-res(1,id)+1._dp-ol_width)/fund_n_par))
            end do
            
            ! test whether end coincides with sum of n_par
            if (res(2,n_div).ne.sum(n_par)) then
                ierr = 1
                err_msg = 'Limits don''t match, try with more memory or lower &
                    &magn_int_style'
                CHCKERR(err_msg)
            end if
        end function calc_eq_jobs_lims
    end function divide_eq_jobs
    
    ! if this equilibrium job should be done, also increment eq_job_nr
    logical function do_eq()
        use num_vars, only: eq_jobs_lims, eq_job_nr, rich_restart_lvl, &
            &jump_to_sol
        use rich_vars, only: rich_lvl
        
        ! increment equilibrium job nr.
        eq_job_nr = eq_job_nr + 1
        
        if (eq_job_nr.le.size(eq_jobs_lims,2)) then
            do_eq = .true.
        else
            do_eq = .false.
        end if
        
        ! possibly skip equilibrium jobs for first Richardson level
        if (rich_lvl.le.rich_restart_lvl .and. jump_to_sol) do_eq = .false.
    end function do_eq
    
    ! Possible  extension with  equilibrium job  as well  as  parallel job  or
    ! nothing if only one level and one parallel job.
    elemental character(len=max_str_ln) function eq_info()
        use num_vars, only: eq_jobs_lims, eq_job_nr
        
        if (size(eq_jobs_lims,2).gt.1) then
            eq_info = ' for Equilibrium job '//trim(i2str(eq_job_nr))//&
                &' of '//trim(i2str(size(eq_jobs_lims,2)))
        else
            eq_info = ''
        end if
    end function eq_info
    
    ! prints information for equilibrium parallel job
    subroutine print_info_eq(n_par_X_rich)
        use num_vars, only: eq_job_nr, eq_jobs_lims
        
        ! input / output
        integer, intent(in) :: n_par_X_rich                                     ! number of parallel points in this Richardson level
        
        ! user output
        if (size(eq_jobs_lims,2).gt.1) then
            call writo('but parallel limits for this job only ['//&
                &trim(i2str(eq_jobs_lims(1,eq_job_nr)))//'..'//&
                &trim(i2str(eq_jobs_lims(2,eq_job_nr)))//&
                &'] of [1..'//trim(i2str(n_par_X_rich))//']')
        end if
    end subroutine print_info_eq
end module eq_utilities
