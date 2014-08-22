!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module utilities
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln
    use output_ops, only: writo, print_ar_1, print_ar_2
    use str_ops, only: i2str, r2strt, r2str
    
    implicit none
    private
    public calc_zero_NR, calc_ext_var, calc_det, calc_int, arr_mult, &
        &VMEC_norm_deriv, VMEC_conv_FHM, check_deriv, calc_inv, &
        &init_utilities, derivs, calc_interp, con2dis, dis2con
    
    ! the possible derivatives of order i
    integer, allocatable :: derivs_0(:,:)                                       ! all possible derivatives of order 0
    integer, allocatable :: derivs_1(:,:)                                       ! all possible derivatives of order 1
    integer, allocatable :: derivs_2(:,:)                                       ! all possible derivatives of order 2
    integer, allocatable :: derivs_3(:,:)                                       ! all possible derivatives of order 3
    integer, allocatable :: derivs_4(:,:)                                       ! all possible derivatives of order 3

    ! interfaces
    interface arr_mult
        module procedure arr_mult_3_3, arr_mult_3_1, arr_mult_1_1
    end interface
    interface calc_det
        module procedure calc_det_0D, calc_det_2D
    end interface
    interface calc_inv
        module procedure calc_inv_0D, calc_inv_2D
    end interface
    interface calc_int
        module procedure calc_int_real, calc_int_complex
    end interface
    interface calc_interp
        module procedure calc_interp_real, calc_interp_complex
    end interface
contains
    ! initialize utilities:
    ! calculate all possible combinations of derivatives of a certain order
    subroutine init_utilities
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: ci, cj, ck, cl                                               ! counters
        
        ! allocate
        allocate(derivs_0(3,1))
        allocate(derivs_1(3,3))
        allocate(derivs_2(3,6))
        allocate(derivs_3(3,10))
        allocate(derivs_4(3,15))
        
        ci = 1
        cj = 1
        ck = 1
        cl = 1
        
        derivs_0 = 0
        derivs_1 = 0
        derivs_2 = 0
        derivs_3 = 0
        derivs_4 = 0
        
        do id = 1,3
            derivs_1(id,ci) = derivs_1(id,ci) + 1
            ci = ci+1
            do jd = 1,id
                derivs_2(id,cj) = derivs_2(id,cj) + 1
                derivs_2(jd,cj) = derivs_2(jd,cj) + 1
                cj = cj+1
                do kd = 1,jd
                    derivs_3(id,ck) = derivs_3(id,ck) + 1
                    derivs_3(jd,ck) = derivs_3(jd,ck) + 1
                    derivs_3(kd,ck) = derivs_3(kd,ck) + 1
                    ck = ck+1
                    do ld = 1,kd
                        derivs_4(id,cl) = derivs_4(id,cl) + 1
                        derivs_4(jd,cl) = derivs_4(jd,cl) + 1
                        derivs_4(kd,cl) = derivs_4(kd,cl) + 1
                        derivs_4(ld,cl) = derivs_4(ld,cl) + 1
                        cl = cl+1
                    end do
                end do
            end do
        end do
    end subroutine
    
    function derivs(order)
        ! input / output
        integer, intent(in) :: order
        integer, allocatable :: derivs(:,:)
        
        select case (order)
            case (0)
                allocate(derivs(3,size(derivs_0,2)))
                derivs = derivs_0
            case (1)
                allocate(derivs(3,size(derivs_1,2)))
                derivs = derivs_1
            case (2)
                allocate(derivs(3,size(derivs_2,2)))
                derivs = derivs_2
            case (3)
                allocate(derivs(3,size(derivs_3,2)))
                derivs = derivs_3
            case (4)
                allocate(derivs(3,size(derivs_4,2)))
                derivs = derivs_4
            case default
        end select
    end function
    
    ! numerically derives a  function whose values are given on  a regular mesh,
    ! specified by  the inverse  step size to  an order specified  by ord  and a
    ! precision specified by prec (which is the  power of the step size to which
    ! the result is still correct. E.g.:  for forward differences, prec = 0, and
    ! for central differences prec=1)
    integer function VMEC_norm_deriv(var,dvar,inv_step,ord,prec) result(ierr)
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'VMEC_norm_deriv'
        
        ! input / output
        real(dp), intent(in) :: var(:), inv_step
        real(dp), intent(inout) :: dvar(:)
        integer, intent(in) :: ord, prec
        
        ! local variables
        integer :: max_n
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        max_n = size(var)
        
        ! tests
        if (size(dvar).ne.max_n) then
            err_msg = 'Derived vector has to be of the same length as input &
                &vector'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (ord.lt.1 .or. ord.gt.max_deriv(3)) then
            err_msg = 'VMEC_norm_deriv can only derive from order 1 &
                &up to order '//trim(i2str(max_deriv(3)))
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! choose correct precision
        select case (prec)
            case (1)
                call prec1
            case default
                err_msg = 'Precision of order '//trim(i2str(prec))//&
                    &' not implemented'
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        subroutine prec1
            ! test whether enough points are given
            if (max_n-2.lt.ord) then
                err_msg = 'For a derivative of order '//&
                    &trim(i2str(ord))//', with precision '//trim(i2str(prec))&
                    &//', at least '//trim(i2str(ord+2))//&
                    &' input values have to be passed'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! apply derivation rules precise up to order 1
            select case (ord)
                case(1)                                                         ! first derivative
                    ! first point
                    dvar(1) = (-3*var(1)+4*var(2)-var(3))*inv_step*0.5
                    ! middle points
                    dvar(2:max_n-1) = (-var(1:max_n-2)+var(3:max_n))*inv_step*0.5
                    ! last point
                    dvar(max_n) = (var(max_n-2) - 4*var(max_n-1) &
                        &+ 3*var(max_n))* inv_step*0.5
                case(2)                                                         ! second derivative
                    ! first point
                    dvar(1) = (2*var(1)-5*var(2)+4*var(3)-var(4))*inv_step**2
                    ! middle points
                    dvar(2:max_n-1) = (var(1:max_n-2)-2*var(2:max_n-1)&
                        &+var(3:max_n))*inv_step**2
                    ! last point
                    dvar(max_n) = (-var(max_n-3)+4*var(max_n-2)-5*var(max_n-1)&
                        &+2*var(max_n))*inv_step**2
                case(3)                                                         ! third derivative
                    ! first point
                    dvar(1) = (-5*var(1)+18*var(2)-24*var(3)+14*var(4)&
                        &-3*var(5))* inv_step**3*0.5
                    ! second point
                    dvar(2) = (-3*var(1)+10*var(2)-12*var(3)+6*var(4)-var(5))* &
                        &inv_step**3*0.5
                    ! middle points
                    dvar(3:max_n-2) = (-var(1:max_n-4)+2*var(2:max_n-3)-&
                        &2*var(4:max_n-1)+var(5:max_n))*inv_step**3*0.5
                    ! next-to-last point
                    dvar(max_n-1) = (var(max_n-4)-6*var(max_n-3)&
                        &+12*var(max_n-2)-10*var(max_n-1)+3*var(max_n))&
                        &*inv_step**3*0.5
                    ! last point
                    dvar(max_n) = (3*var(max_n-4)-14*var(max_n-3)&
                        &+24*var(max_n-2)-18*var(max_n-1)+5*var(max_n))&
                        &*inv_step**3/2
                case(4)                                                         ! fourth derivative
                    ! first point
                    dvar(1) = (3*var(1)-14*var(2)+26*var(3)-24*var(4)+11*var(5)&
                        &-2*var(6))*inv_step**4
                    ! second point
                    dvar(2) = (2*var(1)-9*var(2)+16*var(3)-14*var(4)+6*var(5)&
                        &-var(6))*inv_step**4
                    ! middle points
                    dvar(3:max_n-2) = (var(1:max_n-4)-4*var(2:max_n-3)&
                        &+6*var(3:max_n-2)-4*var(4:max_n-1)+var(5:max_n))&
                        &*inv_step**4
                    ! next-to-last point
                    dvar(max_n-1) = (-var(max_n-5)+6*var(max_n-4)&
                        &-14*var(max_n-3)+16*var(max_n-2)-9*var(max_n-1)&
                        &+2*var(max_n))*inv_step**4
                    ! last point
                    dvar(max_n) = (-2*var(max_n-5)+11*var(max_n-4)-&
                        &24*var(max_n-3)+26*var(max_n-2)-14*var(max_n-1)&
                        &+3*var(max_n))*inv_step**4
                case(5)                                                         ! fifth derivative
                    ! first point
                    dvar(1) = (-7*var(1)+40*var(2)-95*var(3)+120*var(4)&
                        &-85*var(5)+32*var(6)-5*var(7))*inv_step**5*0.5
                    ! second point
                    dvar(2) = (-5*var(1)+28*var(2)-65*var(3)+80*var(4)&
                        &-55*var(5)+20*var(6)-3*var(7))*inv_step**5*0.5
                    ! third point
                    dvar(3) = (-3*var(1)+16*var(2)-35*var(3)+40*var(4)&
                        &-25*var(5)+8*var(6)-var(7))*inv_step**5*0.5
                    ! middle points
                    dvar(4:max_n-3) = (-var(1:max_n-6)+4*var(2:max_n-5)&
                        &-5*var(3:max_n-4)+5*var(5:max_n-2)-4*var(6:max_n-1)&
                        &+var(7:max_n))*inv_step**5*0.5
                    ! next-to-next-to-last point
                    dvar(max_n-2) = (var(max_n-6)-8*var(max_n-5)&
                        &+25*var(max_n-4)-40*var(max_n-3)+35*var(max_n-2)&
                        &-16*var(max_n-1)+3*var(max_n))*inv_step**5*0.5
                    ! next-to-last point
                    dvar(max_n-1) = (3*var(max_n-6)-20*var(max_n-5)&
                        &+55*var(max_n-4)-80*var(max_n-3)+65*var(max_n-2)&
                        &-28*var(max_n-1)+5*var(max_n))*inv_step**5*0.5
                    ! last point
                    dvar(max_n) = (5*var(max_n-6)-32*var(max_n-5)&
                        &+85*var(max_n-4)-120*var(max_n-3)+95*var(max_n-2)&
                        &-40*var(max_n-1)+7*var(max_n))*inv_step**5*0.5
                case default
                    ! This you should never reach!
                    err_msg = 'Derivation of order '//trim(i2str(ord))//&
                        &' not possible in VMEC_norm_deriv'
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end subroutine
    end function VMEC_norm_deriv

    ! numerically interpolates a function that is given on either FM to HM or to
    ! FM. If FM2HM is .true., the starting variable is FM, if .false., it is HM
    integer function VMEC_conv_FHM(var,cvar,FM2HM) result(ierr)
        character(*), parameter :: rout_name = 'VMEC_conv_FHM'
        
        ! input / output
        real(dp), intent(in) :: var(:)
        real(dp), intent(inout) :: cvar(:)
        logical, intent(in) :: FM2HM
        
        ! local variables
        integer :: max_n
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        max_n = size(var)
        
        ! tests
        if (size(cvar).ne.max_n) then
            err_msg = 'The converted vector has to be of the same length as &
                &the input vector'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        if (FM2HM) then                                                         ! FM to HM
            cvar(1) = 0.0_dp
            cvar(2:max_n) = (var(1:max_n-1)+var(2:max_n))/2.0_dp
        else                                                                    ! HM to FM
            cvar(1) = (3.0_dp*var(2)-var(3))/2.0_dp
            cvar(2:max_n-1) = (var(2:max_n-1)+var(3:max_n))/2.0_dp
            cvar(max_n) = (-var(max_n-1)+3.0_dp*var(max_n))/2.0_dp
        end if
    end function
    
    ! checks  whether the  derivatives requested  for a  certain subroutine  are
    ! valid
    integer function check_deriv(deriv,max_deriv,sr_name) result(ierr)
        character(*), parameter :: rout_name = 'check_deriv'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        integer, intent(in) :: max_deriv(3)
        character(len=*), intent(in) :: sr_name
        
        ! local variables
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test the derivatives
        do id = 1, 3
            if (deriv(id).gt.max_deriv(id)) then
                err_msg = 'Asking '//trim(sr_name)//' for a derivative&
                    & in the '//trim(i2str(id))//'th VMEC coordinate of order '&
                    &//trim(i2str(deriv(id)))//', while the maximum order is '&
                    &//trim(i2str(max_deriv(id)))
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do
    end function

    ! integrates a function using the trapezoidal rule using constant step size:
    !   int_1^n f(x) dx = (x(2)-x(1))/2 f(a) + (x(n)-x(n-1))/2 f(n) 
    !                     + sum_k=2^(n-1){f(k)*(x(k+1)-x(k-1))/2},
    ! with n the number of points. So, n  points have to be specified as well as
    ! n values  for the function  to be interpolated. They  have to be  given in
    ! ascending order but the step size does not have to be constant
    ! this yields the following difference formula:
    !   int_1^n f(x) dx = int_1^(n-1) f(x) dx + (f(n)+f(n+1))*(x(n+1)-x(n))/2,
    ! which is used here
    integer function calc_int_real(var,x,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_real'
        
        ! input / output
        real(dp) :: var_int(:)
        real(dp) :: var(:)
        real(dp) :: x(:)
        
        ! local variables
        integer :: n_max
        integer :: kd
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        n_max = size(var)
        
        ! tests
        if (size(x).ne.n_max) then
            err_msg = 'The arrays of points and values are not of the same &
                &length'
            ierr = 1
            CHCKERR(err_msg)
        else if (n_max.lt.2) then
            err_msg = 'Asking calc_int to integrate with '//trim(i2str(n_max))&
                &//' points. Need at least 2'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! calculate integral for all points
        var_int = 0.0_dp
        
        do kd = 2, n_max
            var_int(kd) = var_int(kd-1) + &
                &(var(kd)+var(kd-1))/2 * (x(kd)-x(kd-1))
        end do
    end function calc_int_real
    integer function calc_int_complex(var,x,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_complex'
        
        ! input / output
        complex(dp) :: var_int(:)
        complex(dp) :: var(:)
        real(dp) :: x(:)
        
        ! local variables
        integer :: n_max
        integer :: kd
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        n_max = size(var)
        
        ! tests
        if (size(x).ne.n_max) then
            err_msg = 'The arrays of points and values are not of the same &
                &length'
            CHCKERR(err_msg)
        else if (n_max.lt.2) then
            err_msg = 'Asking calc_int to integrate with '//trim(i2str(n_max))&
                &//' points. Need at least 2'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! calculate integral for all points
        var_int = 0.0_dp
        
        do kd = 2, n_max
            var_int(kd) = var_int(kd-1) + &
                &(var(kd)+var(kd-1))/2 * (x(kd)-x(kd-1))
        end do
    end function calc_int_complex
    
    ! extrapolates  a   function,  using  linear  or   quadratic  interpolation,
    ! depending on  the number of  points and values  given. The data  should be
    ! given sorted, in ascending order, without duplicate points in var_points.
    ! It uses the following solution: A b = c
    ! | x_1^0 x_1^1 x_1^2 ... | | a_0 |   | var(1) |
    ! | x_2^0 x_2^1 x_2^2 ... | | a_1 | = | var(2) |
    ! | x_3^0 x_3^1 x_3^2 ... | | a_2 |   | var(3) |
    ! | ...   ...   ...   ... | | ... |   |  ...   |
    ! for the interpolating polynomial ext_var = a_0 + a_1 x + a_2 x^2 + ...
    ! There is an optional flag to look  for the k^th derivative of the function
    ! instead of the function itself, where k should be lower than the degree of
    ! the polynomial
    integer function calc_ext_var(ext_var,var,var_points,ext_point,deriv_in) &
        &result(ierr)
        character(*), parameter :: rout_name = 'ext_var'
        
        ! input / output
        real(dp), intent(inout) :: ext_var                                      ! output
        real(dp), intent(in) :: var(:), var_points(:)                           ! input function and points at which tabulated
        real(dp), intent(in) :: ext_point                                       ! point to which extrapolate
        integer, intent(in), optional :: deriv_in                               ! specifies an optional derivative
        
        ! local variables
        integer :: id, jd
        integer :: istat
        integer :: pol_deg
        integer :: deriv
        real(dp), allocatable :: A(:,:), LU(:,:)
        real(dp), allocatable :: a_i(:)
        real(dp) :: fact                                                        ! multiplicative factor, see below
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! determine the degree of the polynomial ext_var
        pol_deg = size(var)
        
        ! tests
        if (size(var).ne.size(var_points)) then
            err_msg = 'The size of the arrays containing the variable and &
                &the points do not match'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (present(deriv_in)) then 
            deriv = deriv_in
            if (deriv.ge.pol_deg) then                                          ! order of derivative too hight
                err_msg = 'Asking ext_var for '//trim(i2str(deriv_in))&
                    &//'th derivative, while using a polynomial of degree '&
                    &//trim(i2str(pol_deg))
                ierr = 1
                CHCKERR(err_msg)
            else if (deriv.lt.0) then
                err_msg = 'Asking ext_var for a derivative of order '&
                    &//trim(i2str(deriv_in))
                ierr = 1
                CHCKERR(err_msg)
            end if
        else
            deriv = 0
        end if
        
        ! fill matrix
        allocate(A(pol_deg,pol_deg))
        allocate(LU(pol_deg,pol_deg))
        A(:,1) = 1.0_dp
        col: do jd = 2, pol_deg
            row: do id = 1, pol_deg
                A(id,jd) = var_points(id)**(jd-1)
            end do row
        end do col
        
        ! solve the system for a_i
        allocate(a_i(pol_deg))
        a_i = var
        LU = A
        call dgesv( pol_deg, 1, LU, pol_deg, [(id,id=1,pol_deg)], &
            &a_i, pol_deg, istat)
        
        
        if (istat.ne.0) then
            err_msg = 'The solution could not be found. Are the extrapolating &
                &points independent?'
            ierr = 1
            if (ierr.ne.0) then
                call writo('The matrix equation to be solved was Ax=b, A = ')
                call print_ar_2(A)
                call writo('and b = ')
                call print_ar_1(var)
            end if
            CHCKERR(err_msg)
        end if
        
        ext_var = 0.0_dp
        ! for deriv of order i: 
        ! d^i/dx^i f = sum_j = i ^ pol_deg ( j*(j-1)*...*(j-i) a_j x^(j-i) )
        polynom: do id = 1+deriv, pol_deg
            ! construct the factor 
            fact = 1.0_dp
            d_dx: do jd = 1,deriv
                fact = fact * float((id-1)-(jd-1))
            end do d_dx
            ! update ext_var with current polynomial term a_i(id)
            ext_var = ext_var + fact*a_i(id)*ext_point**(id-1-deriv)
        end do polynom
    end function calc_ext_var
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    ! VERSION with arr_1 and arr_2 in three coords.
    integer function arr_mult_3_3(arr_1,arr_2,arr_3,deriv) result(ierr)
        character(*), parameter :: rout_name = 'arr_mult_3_3'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: arr_2(1:,1:,0:,0:,0:)
        real(dp), intent(out) :: arr_3(1:,1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        integer :: r, m, n                                                      ! current degree of norm., pol., tor. derivatives
        integer :: kd                                                           ! normal counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(arr_1,1).ne.size(arr_2,1) .or. size(arr_1,2).ne.size(arr_2,2)) then
            err_msg = 'Arrays 1 and 2 need to have the same size'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(arr_1,1).ne.size(arr_3,1) .or. size(arr_1,2).ne.size(arr_3,2)) then
            err_msg = 'Arrays 1 and 2 need to have the same size as the &
                &resulting array 3'
            ierr = 1
            CHCKERR(err_msg)
        end if
        do kd = 1,3
            if (size(arr_1,2+kd).lt.deriv(kd)+1 .or. &
                &size(arr_2,2+kd).lt.deriv(kd)+1) then
                err_msg = 'Arrays 1 and 2 do not provide the necessary &
                    &orders of derivatives to calculate a derivative of order &
                    &['//trim(i2str(deriv(1)))//','//trim(i2str(deriv(2)))//&
                    &','//trim(i2str(deriv(3)))//', at least in coordinate '&
                    &//trim(i2str(kd))
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do
        
        ! distribute normal and angular derivatives using binomial theorem
        ! normal derivatives
        do r = 0,deriv(1)
            if (r.eq.0) then                                                    ! first term in sum
                bin_fac(1) = 1.0_dp
            else
                bin_fac(1) = bin_fac(1)*(deriv(1)-(r-1))/r
            end if
            ! poloidal derivatives
            do m = 0,deriv(2)
                if (m.eq.0) then                                                ! first term in sum
                    bin_fac(2) = 1.0_dp
                else
                    bin_fac(2) = bin_fac(2)*(deriv(2)-(m-1))/m
                end if
                ! toroidal derivatives
                do n = 0,deriv(3)
                    if (n.eq.0) then                                            ! first term in sum
                        bin_fac(3) = 1.0_dp
                    else
                        bin_fac(3) = bin_fac(3)*(deriv(3)-(n-1))/n
                    end if
                    
                    ! current term in the tripple summation
                    do kd = 1, size(arr_1,2)
                        arr_3(:,kd) = arr_3(:,kd) + &
                            &bin_fac(1)*bin_fac(2)*bin_fac(3) &
                            &* arr_1(:,kd,r,m,n)&
                            &* arr_2(:,kd,deriv(1)-r,deriv(2)-m,deriv(3)-n)
                    end do
                end do
            end do
        end do
    end function arr_mult_3_3
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    ! VERSION with arr_1 in three coords and arr_2 only in the flux coord.
    integer function arr_mult_3_1(arr_1,arr_2,arr_3,deriv) result(ierr)
        character(*), parameter :: rout_name = 'arr_mult_3_1'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: arr_2(1:,0:)
        real(dp), intent(out) :: arr_3(1:,1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac                                                     ! binomial factor for norm. sum
        integer :: r                                                            ! current degree of norm. derivative
        integer :: kd                                                           ! normal counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(arr_1,2).ne.size(arr_2,1)) then
            err_msg = 'Arrays 1 and 2 need to have the same size in the flux &
                &coord.'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(arr_1,1).ne.size(arr_3,1) .or. &
            &size(arr_1,2).ne.size(arr_3,2)) then
            err_msg = 'Array 1 needs to have the same size as the resulting &
                &array 3'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! distribute   normal  derivatives   using  binomial   theorem,  angular
        ! derivatives only operate on the first array
        ! normal derivatives
        do r = 0,deriv(1)
            if (r.eq.0) then                                                    ! first term in sum
                bin_fac = 1.0_dp
            else
                bin_fac = bin_fac*(deriv(1)-(r-1))/r
            end if
            
            ! current term in the tripple summation
            do kd = 1, size(arr_1,2)
                arr_3(:,kd) = arr_3(:,kd) + bin_fac * &
                    &arr_1(:,kd,r,deriv(2),deriv(3))* arr_2(kd,deriv(1)-r)
            end do
        end do
    end function arr_mult_3_1
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    ! VERSION with arr_1 and arr_2 only in the flux coord.
    integer function arr_mult_1_1(arr_1,arr_2,arr_3,deriv) result(ierr)
        character(*), parameter :: rout_name = 'arr_mult_1_1'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,0:)
        real(dp), intent(in) :: arr_2(1:,0:)
        real(dp), intent(out) :: arr_3(1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac                                                     ! binomial factor for norm. sum
        integer :: r                                                            ! current degree of norm. derivative
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(arr_1,1).ne.size(arr_2,1)) then
            err_msg = 'Arrays 1 and 2 need to have the same size in the &
                &flux coord.'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(arr_1,1).ne.size(arr_3,1)) then
            err_msg = 'Array 1 needs to have the same size as the resulting &
                &array 3'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        if (deriv(2).gt.0 .or. deriv(3).gt.0) then                              ! no addition to arr_3
            return
        else
            ! distribute  normal  derivatives  using binomial  theorem,  angular
            ! derivatives only operate on the first array
            ! normal derivatives
            do r = 0,deriv(1)
                if (r.eq.0) then                                                    ! first term in sum
                    bin_fac = 1.0_dp
                else
                    bin_fac = bin_fac*(deriv(1)-(r-1))/r
                end if
                
                ! current term in the tripple summation
                arr_3 = arr_3 + bin_fac * arr_1(:,r)* arr_2(:,deriv(1)-r)
            end do
        end if
    end function arr_mult_1_1

    ! calculate determinant of a matrix which is defined on a 2D grid
    ! the size should be small, as the direct formula is used.
    integer recursive function calc_det_2D(detA,A) result (ierr)
        character(*), parameter :: rout_name = 'calc_det_2D'
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! input
        real(dp), intent(inout) :: detA(:,:)                                    ! output
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n                                                            ! holds size of A
        real(dp) :: sgn                                                         ! holds either plus or minus one
        integer, allocatable :: slct(:)                                         ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: work(:,:)                                      ! work array
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,3).ne.size(A,4)) then
            err_msg = 'The matrix A has to be square'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(A,3).gt.4) then
            err_msg = 'This should only be used for small matrices'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(detA,1).ne.size(A,1) .or. size(detA,2).ne.size(A,2)) then
            err_msg = 'The output and input matrix have to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! intialize
        n = size(A,3)
        sgn = 1.0_dp
        allocate (idx(n))
        idx = [(id,id=1,n)]
        allocate (slct(n))
        slct = 1
        allocate (work(size(A,1),size(A,2)))
        
        detA = 0.0_dp
        
        if (n.eq.1) then                                                        ! shouldn't be used
            detA = A(:,:,1,1)
        else if (n.eq.2) then
            detA = A(:,:,1,1)*A(:,:,2,2)-A(:,:,1,2)*A(:,:,2,1)
        else
            do id = 1, n
                slct(id) = 0
                ierr = calc_det_2D(work,A(:,:,2:n,pack(idx,slct.gt.0)))
                CHCKERR('')
                detA = detA + A(:,:,1,id)*sgn*work
                sgn = -sgn
                slct(id) = 1
            end do
        end if
    end function calc_det_2D

    ! calculate determinant of a constant matrix
    ! (adapted from http://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/)
    ! ¡¡¡WARNING: the matrix A is changed on output!!!
    integer function calc_det_0D(det_0D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_det_0D'
        
        ! input / output
        real(dp), intent(inout) :: A(:,:)                                       ! input
        real(dp), intent(inout) :: det_0D                                       ! output
        
        ! local variables
        integer :: n 
        real(dp) :: sgn
        integer :: id, info
        integer, allocatable :: ipiv(:)
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,1).ne.size(A,2)) then
            err_msg = 'The matrix A has to be square'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initialize
        n = size(A,1)
        allocate(ipiv(n))
        ipiv = 0
        
        call dgetrf(n, n, A, n, ipiv, info)
        
        det_0D = 1.0_dp
        do id = 1, n
            det_0D = det_0D*A(id, id)
        end do
        
        sgn = 1.0_dp
        do id = 1, n
            if(ipiv(id).ne.id) then
                sgn = -sgn
            end if
        end do
        det_0D = sgn*det_0D
    end function calc_det_0D
    
    ! calculate inverse of  square matrix A which has  elements depending on
    ! 2D mesh (first two coords)
    ! this method uses direct inversion using  Cramer's rule, since the matrix A
    ! is supposed to be  very small (i.e. 3x3 or 1x1) and  since the inverse has
    ! to be calculated at  each of the points in the 2D mesh,  this can be quite
    ! fast
    integer function calc_inv_2D(inv_2D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_2D'
        
        ! input / output
        real(dp), intent(inout) :: inv_2D(:,:,:,:)                              ! output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! input
        
        ! local variables
        real(dp), allocatable :: detA(:,:)
        integer :: id, kd                                                       ! counters
        integer :: n
        integer, allocatable :: slct(:,:)                                       ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        character(len=max_str_ln) :: err_msg                                    ! error message
       
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,3).ne.size(A,4)) then
            err_msg = 'The input matrix has to be square'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(inv_2D,1).ne.size(A,1) .or. size(inv_2D,2).ne.size(A,2) .or. &
            &size(inv_2D,3).ne.size(A,3) .or. size(inv_2D,4).ne.size(A,4)) then
            err_msg = 'The output and input matrix have to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initializing
        n = size(A,3)
        allocate (idx(n))
        idx = [(id,id=1,n)]
        allocate (slct(n,2))
        slct = 1
        allocate(detA(size(A,1),size(A,2)))
        detA = 0.0_dp
        
        ! calculate cofactor matrix, replacing original elements
        ! use is made of detA
        do kd = 1,n
            do id = 1,n
                slct(kd,1) = 0                                                  ! take out row kd
                slct(id,2) = 0                                                  ! take out column id
                ierr = calc_det(detA,&
                    &A(:,:,pack(idx,slct(:,1).gt.0),pack(idx,slct(:,2).gt.0)))
                CHCKERR('')
                inv_2D(:,:,id,kd) = (-1.0_dp)**(id+kd)*detA
                slct(kd,1) = 1
                slct(id,2) = 1
            end do
        end do
        
        ! calculate determinant in detA
        ierr = calc_det(detA,A)
        CHCKERR('')
        
        ! divide by determinant
        do kd = 1,n
            do id = 1,n
                inv_2D(:,:,id,kd) = inv_2D(:,:,id,kd) / detA
            end do
        end do
    end function calc_inv_2D
    ! calculate inverse of a constant square matrix
    integer function calc_inv_0D(inv_0D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_0D'
        
        ! input / output
        real(dp), intent(inout) :: inv_0D(:,:)                                  ! output
        real(dp), intent(in) :: A(:,:)                                          ! input
        
        ! local variables
        integer :: n 
        integer, allocatable :: ipiv(:)                                         ! pivot variable, used by Lapack
        real(dp), allocatable :: w(:)                                           ! work variable
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,1).ne.size(A,1)) then
            err_msg = 'The input matrix A has to be square'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(inv_0D,1).ne.size(A,1) .or. size(inv_0D,2).ne.size(A,2)) then
            err_msg = 'The output and input matrix have to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initialize
        n = size(A,1)
        allocate(ipiv(n))
        ipiv = 0
        allocate(w(n))
        
        inv_0D = A
        
        call dgetrf(n, n, inv_0D, n, ipiv, ierr)                                ! LU factorization
        err_msg = 'Lapack couldn''t find the LU factorization'
        CHCKERR(err_msg)
        
        call dgetri(n, inv_0D, n, ipiv, w, n, ierr)                             ! inverse of LU
        CHCKERR('Lapack couldn''t compute the inverse')
    end function calc_inv_0D

    ! finds the zero of a function using Newton-Rhapson iteration
    integer function calc_zero_NR(zero_NR,fun,dfun,guess) result(ierr)
        use num_vars, only: max_it_NR, tol_NR
        
        character(*), parameter :: rout_name = 'calc_zero_NR'
        
        ! input / output
        real(dp), intent(inout) :: zero_NR                                      ! output
        real(dp), intent(in) :: guess                                           ! first guess
        interface
            function fun(x)                                                     ! the function
                use num_vars, only: dp
                real(dp) :: fun
                real(dp), intent(in) :: x
            end function fun
            function dfun(x)                                                    ! derivative of the function
                use num_vars, only: dp
                real(dp) :: dfun
                real(dp), intent(in) :: x
            end function dfun
        end interface
        
        ! local variables
        integer :: jd
        real(dp) :: corr
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        zero_NR = guess
        
        NR: do jd = 1,max_it_NR
            ! correction to theta_NR
            corr = -fun(zero_NR)/dfun(zero_NR)
            !write(*,*) 'jd, zero_NR = ', jd, zero_NR
            !write(*,*) 'fun, dfun = ', fun(zero_NR), dfun(zero_NR)
            zero_NR = zero_NR + corr
            
            ! check for convergence
            ! CHANGE THIS TO RELATIVE
            if (abs(corr).lt.tol_NR) then
                exit
            else if (jd .eq. max_it_NR) then
                err_msg = 'Not converged after '//trim(i2str(jd))//&
                    &' iterations, with residual = '//trim(r2str(corr))
                zero_NR = 0.0_dp
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do NR
    end function calc_zero_NR

    ! convert  between points  on  a continuous  grid  (.e.g. perturbation  grid
    ! (min_r_X..max_r_X))  and  on  a  discrete  grid  (.e.g.  equilibrium  grid
    ! (1..n_r)),  taking  into  account  the   limits  on  the  continuous  grid
    ! (lim_cont) and the limits on the discrete grid (lim_dis)
    ! The formula used is:
    !    pt_c - min_c     pt_d - min_d
    !   ------------- =  ------------- ,
    !   max_c - min_c    max_d - min_d
    ! note: for con2dis,  the limits of the discrete grid  are given as discrete
    ! values, while the output  pt_d is given back as a real  value, which is to
    ! be rounded by the calling subroutine
    subroutine con2dis(pt_c,lim_c,pt_d,lim_d)
        ! input / output
        real(dp), intent(in) :: pt_c                                            ! point on continous grid
        real(dp), intent(in) :: lim_c(2)                                        ! [min_c,max_c]
        real(dp), intent(inout) :: pt_d                                         ! point on discrete grid
        integer, intent(in) :: lim_d(2)                                         ! [min_d,max_d]
        
        ! calculate
        pt_d = lim_d(1) + (lim_d(2)-lim_d(1)) * (pt_c-lim_c(1)) / &
            &(lim_c(2)-lim_c(1))
    end subroutine con2dis
    subroutine dis2con(pt_d,lim_d,pt_c,lim_c)
        ! input / output
        integer, intent(in) :: pt_d                                             ! point on discrete grid
        integer, intent(in) :: lim_d(2)                                         ! [min_d,max_d]
        real(dp), intent(inout) :: pt_c                                         ! point on continous grid
        real(dp), intent(in) :: lim_c(2)                                        ! [min_c,max_c]
        
        ! calculate
        pt_c = lim_c(1) + (lim_c(2)-lim_c(1)) * (pt_d-lim_d(1)) / &
            &(lim_d(2)-lim_d(1))
    end subroutine dis2con
    
    ! Simple  linear interpolation  of  matrix varin,  tabulated at  equilibrium
    ! grid,  at point  ptout where  ptout =  0..1 with  result stored  in matrix
    ! varout
    ! pt_arr is  the index in  the array corresponding  to ptout, which  is then
    ! rounded above  and below to  yield the  indices ind_lo and  ind_hi between
    ! which to interpolate
    ! The interpolation between ind_lo and ind_hi is done linearly:
    !   varout = varin(ind_lo) + (pt_arr-ind_lo) * (varin(ind_hi)-varin(ind_lo))
    ! because ind_hi - ind_lo = 1
    ! note: this routine is also correct if ind_hi = ind_lo (for the last point)
    integer function calc_interp_real(varin,limin,varout,ptout,r_offset) &
        &result(ierr)
        character(*), parameter :: rout_name = 'calc_interp_real'
        
        ! input / output
        real(dp), intent(in) :: varin(:,:,:,:)                                  ! input variable, not interpolated
        integer, intent(in) :: limin(2)                                         ! start and end points at which variable is tabulated
        real(dp), intent(inout) :: varout(:,:,:)                                ! output variable, interpolated
        real(dp), intent(in) :: ptout                                           ! point at which to interpolate (0...1)
        integer, intent(in), optional :: r_offset                               ! offset in the tables in the normal variable, wrt. 1
        
        ! local variables
        real(dp) :: ptout_loc                                                   ! local copy of ptout
        integer :: ind_lo, ind_hi                                               ! lower and upper index
        real(dp) :: pt_arr                                                      ! point in array referring to ptout
        real(dp) :: margin = 1E-4_dp                                            ! margin for the comparisons
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set local copy of ptout
        ptout_loc = ptout
        
        ! tests
        if (size(varin,1).ne.size(varout,1) .or. size(varin,3).ne.&
            &size(varout,2) .or. size(varin,4).ne.size(varout,3)) then
            err_msg = 'The sizes of varin and varout have to match!'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (ptout_loc.lt.0) then
            if (ptout_loc.lt.0-margin) then
                err_msg = 'ptout has to lie within 0..1, yet it is equal &
                    &to '//trim(r2str(ptout))
                ierr = 1
                CHCKERR(err_msg)
            else
                ptout_loc = 0.0_dp
            end if
        end if
        if (ptout_loc.gt.1) then
            if (ptout_loc.gt.1+margin) then
                err_msg = 'ptout has to lie within 0..1, yet it is equal &
                    &to '//trim(r2str(ptout))//'...'
                ierr = 1
                CHCKERR(err_msg)
            else
                ptout_loc = 1.0_dp
            end if
        end if
        
        ! find interpolation points
        pt_arr = limin(1) +(limin(2)-limin(1))*ptout_loc                        ! corresponding array index, not rounded
        ind_lo = floor(pt_arr)                                                  ! lower bound
        ind_hi = ceiling(pt_arr)                                                ! upper bound
        
        ! include offset
        if (present(r_offset)) then
            ind_lo = ind_lo - r_offset
            ind_hi = ind_hi - r_offset
        end if
        
        ! interpolate
        varout = varin(:,ind_lo,:,:) + (pt_arr-ind_lo) * &
            &(varin(:,ind_hi,:,:)-varin(:,ind_lo,:,:))
    end function calc_interp_real
    integer function calc_interp_complex(varin,limin,varout,ptout,r_offset) &
        &result(ierr)
        character(*), parameter :: rout_name = 'calc_interp_complex'
        
        ! input / output
        complex(dp), intent(in) :: varin(:,:,:,:)                               ! input variable, not interpolated
        integer, intent(in) :: limin(2)                                         ! start and end points at which variable is tabulated
        complex(dp), intent(inout) :: varout(:,:,:)                             ! output variable, interpolated
        real(dp), intent(in) :: ptout                                           ! point at which to interpolate (0...1)
        integer, intent(in), optional :: r_offset                               ! offset in the tables in the normal variable, wrt. 1
        
        ! local variables
        real(dp) :: ptout_loc                                                   ! local copy of ptout
        integer :: ind_lo, ind_hi                                               ! lower and upper index
        real(dp) :: pt_arr                                                      ! point in array referring to ptout
        real(dp) :: margin = 1E-4_dp                                            ! margin for the comparisons
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set local copy of ptout
        ptout_loc = ptout
        
        ! tests
        if (size(varin,1).ne.size(varout,1) .or. size(varin,3).ne.&
            &size(varout,2) .or. size(varin,4).ne.size(varout,3)) then
            err_msg = 'The sizes of varin and varout have to match!'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (ptout_loc.lt.0) then
            if (ptout_loc.lt.0-margin) then
                err_msg = 'ptout has to lie within 0..1, yet it is equal &
                    &to '//trim(r2str(ptout))
                ierr = 1
                CHCKERR(err_msg)
            else
                ptout_loc = 0.0_dp
            end if
        end if
        if (ptout_loc.gt.1) then
            if (ptout_loc.gt.1+margin) then
                err_msg = 'ptout has to lie within 0..1, yet it is equal &
                    &to '//trim(r2str(ptout))//'...'
                ierr = 1
                CHCKERR(err_msg)
            else
                ptout_loc = 1.0_dp
            end if
        end if
        
        ! find interpolation points
        pt_arr = limin(1) +(limin(2)-limin(1))*ptout_loc                        ! corresponding array index, not rounded
        ind_lo = floor(pt_arr)                                                  ! lower bound
        ind_hi = ceiling(pt_arr)                                                ! upper bound
        
        ! interpolate
        if (present(r_offset)) then                                             ! include offset
            varout = varin(:,ind_lo-r_offset,:,:) + (pt_arr-ind_lo) * &
                &(varin(:,ind_hi-r_offset,:,:)-varin(:,ind_lo-r_offset,:,:))
        else
            varout = varin(:,ind_lo,:,:) + (pt_arr-ind_lo) * &
                &(varin(:,ind_hi,:,:)-varin(:,ind_lo,:,:))
        end if
    end function calc_interp_complex
end module utilities
