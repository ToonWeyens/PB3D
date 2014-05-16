!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module utilities
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_1, print_ar_2
    use str_ops, only: i2str, r2strt, r2str
    
    implicit none
    private
    public zero_NR, ext_var, det, calc_int, arr_mult, find_VMEC_norm_coord, &
        &VMEC_norm_deriv, VMEC_conv_FHM, check_deriv, inv

    interface arr_mult
        module procedure arr_mult_3_3, arr_mult_3_1, arr_mult_1_1
    end interface
    interface det
        module procedure det_0D, det_2D
    end interface
    interface inv
        module procedure inv_0D, inv_2D
    end interface
contains
    ! numerically derives a  function whose values are given on  a regular mesh,
    ! specified by  the inverse  step size to  an order specified  by ord  and a
    ! precision specified by prec (which is the  power of the step size to which
    ! the result is still correct. E.g.:  for forward differences, prec = 0, and
    ! for central differences prec=1)
    subroutine VMEC_norm_deriv(var,dvar,inv_step,ord,prec)
        use num_vars, only: max_deriv
        
        ! input / output
        real(dp), intent(in) :: var(:), inv_step
        real(dp), intent(inout) :: dvar(:)
        integer, intent(in) :: ord, prec
        
        ! local variables
        integer :: max_n
        
        max_n = size(var)
        
        ! tests
        if (size(dvar).ne.max_n) then
            call writo('ERROR: In VMEC_norm_deriv, the derived vector has to &
                &be of the same length as the input vector')
            stop
        end if
        if (ord.lt.1 .or. ord.gt.max_deriv(3)) then
            call writo('ERROR: VMEC_norm_deriv can only derive from order 1 &
                &up to order '//trim(i2str(max_deriv(3))))
            stop
        end if
        
        ! choose correct precision
        select case (prec)
            case (1)
                call prec1
            case default
                call writo('ERROR: precision of order '//trim(i2str(prec))//&
                    &' not implemented')
                stop
        end select
    contains
        subroutine prec1
            ! test whether enough points are given
            if (max_n-2.lt.ord) then
                call writo('ERROR: for a derivative of order '//&
                    &trim(i2str(ord))//', with precision '//trim(i2str(prec))&
                    &//', at least '//trim(i2str(ord+2))//&
                    &' input values have to be passed')
                stop
            end if
            
            ! apply derivation rules precise up to order 1
            select case (ord)
                case(1)                                                         ! first derivative
                    ! first point
                    dvar(1) = (-3*var(1)+4*var(2)-var(3))*inv_step/2
                    ! middle points
                    dvar(2:max_n-1) = (-var(1:max_n-2)+var(3:max_n))*inv_step/2
                    ! last point
                    dvar(max_n) = (var(max_n-2) - 4*var(max_n-1) &
                        &+ 3*var(max_n))* inv_step/2
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
                        &-3*var(5))* inv_step**3/2
                    ! second point
                    dvar(2) = (-3*var(1)+10*var(2)-12*var(3)+6*var(4)-var(5))* &
                        &inv_step**3/2
                    ! middle points
                    dvar(3:max_n-2) = (-var(1:max_n-4)+2*var(2:max_n-3)-&
                        &2*var(4:max_n-1)+var(5:max_n))*inv_step**3/2
                    ! next-to-last point
                    dvar(max_n-1) = (var(max_n-4)-6*var(max_n-3)&
                        &+12*var(max_n-2)-10*var(max_n-1)+3*var(max_n))&
                        &*inv_step**3/2
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
                        &-85*var(5)+32*var(6)-5*var(7))*inv_step**5/2
                    ! second point
                    dvar(2) = (-5*var(1)+28*var(2)-65*var(3)+80*var(4)&
                        &-55*var(5)+20*var(6)-3*var(7))*inv_step**5/2
                    ! third point
                    dvar(3) = (-3*var(1)+16*var(2)-35*var(3)+40*var(4)&
                        &-25*var(5)+8*var(6)-var(7))*inv_step**5/2
                    ! middle points
                    dvar(4:max_n-3) = (-var(1:max_n-6)+4*var(2:max_n-5)&
                        &-5*var(3:max_n-4)+5*var(5:max_n-2)-4*var(6:max_n-1)&
                        &+var(7:max_n))*inv_step**5/2
                    ! next-to-next-to-last point
                    dvar(max_n-2) = (var(max_n-6)-8*var(max_n-5)&
                        &+25*var(max_n-4)-40*var(max_n-3)+35*var(max_n-2)&
                        &-16*var(max_n-1)+3*var(max_n))*inv_step**5/2
                    ! next-to-last point
                    dvar(max_n-1) = (3*var(max_n-6)-20*var(max_n-5)&
                        &+55*var(max_n-4)-80*var(max_n-3)+65*var(max_n-2)&
                        &-28*var(max_n-1)+5*var(max_n))*inv_step**5/2
                    ! last point
                    dvar(max_n) = (5*var(max_n-6)-32*var(max_n-5)&
                        &+85*var(max_n-4)-120*var(max_n-3)+95*var(max_n-2)&
                        &-40*var(max_n-1)+7*var(max_n))*inv_step**5/2
                case default
                    ! This you should never reach!
                    call writo('ERROR: Derivation of order '//&
                        &trim(i2str(ord))//' not possible in VMEC_norm_deriv')
            end select
        end subroutine
    end subroutine

    ! numerically interpolates a function that is given on either FM to HM or to
    ! FM. If FM2HM is .true., the starting variable is FM, if .false., it is HM
    subroutine VMEC_conv_FHM(var,cvar,FM2HM)
        ! input / output
        real(dp), intent(in) :: var(:)
        real(dp), intent(inout) :: cvar(:)
        logical, intent(in) :: FM2HM
        
        ! local variables
        integer :: max_n
        
        max_n = size(var)
        
        ! tests
        if (size(cvar).ne.max_n) then
            call writo('ERROR: In VMEC_conv_FHM, the converted vector has to &
                &be of the same length as the input vector')
            stop
        end if
        
        if (FM2HM) then                                                         ! FM to HM
            cvar(1) = 0.0_dp
            cvar(2:max_n) = (var(1:max_n-1)+var(2:max_n))/2.0_dp
        else                                                                    ! HM to FM
            cvar(1) = (3.0_dp*var(2)-var(3))/2.0_dp
            cvar(2:max_n-1) = (var(2:max_n-1)+var(3:max_n))/2.0_dp
            cvar(max_n) = (-var(max_n-1)+3.0_dp*var(max_n))/2.0_dp
        end if
    end subroutine
    
    ! checks  whether the  derivatives requested  for a  certain subroutine  are
    ! valid
    subroutine check_deriv(deriv,max_deriv,sr_name)
        ! input / output
        integer, intent(in) :: deriv(3)
        integer, intent(in) :: max_deriv(3)
        character(len=*), intent(in) :: sr_name
        
        ! local variables
        integer :: id
        
        ! test the derivatives
        do id = 1, 3
            if (deriv(id).gt.max_deriv(id)) then
                call writo('ERROR: Asking '//trim(sr_name)//' for a derivative&
                    & in the '//trim(i2str(id))//'th VMEC coordinate of order '&
                    &//trim(i2str(deriv(id)))//', while the maximum order is '&
                    &//trim(i2str(max_deriv(id))))
                stop
            end if
        end do
    end subroutine

    ! integrates a function using the trapezoidal rule using constant step size:
    !   int_1^n f(x) dx = (x(2)-x(1))/2 f(a) + (x(n)-x(n-1))/2 f(n) 
    !                     + sum_k=2^(n-1){f(k)*(x(k+1)-x(k-1))/2},
    ! with n the number of points. So, n  points have to be specified as well as
    ! n values  for the function  to be interpolated. They  have to be  given in
    ! ascending order but the step size does not have to be constant
    ! this yields the following difference formula:
    !   int_1^n f(x) dx = int_1^(n-1) f(x) dx + (f(n)+f(n+1))*(x(n+1)-x(n))/2,
    ! which is used here
    function calc_int(var,x)
        ! input / output
        real(dp), allocatable :: calc_int(:)
        real(dp) :: var(:)
        real(dp) :: x(:)
        
        ! local variables
        integer :: n_max
        integer :: kd
        
        n_max = size(var)
        
        ! tests
        if (size(x).ne.n_max) then
            call writo('ERROR: in calc_int, the arrays of points and values &
                &are not of the same length')
            stop
        else if (n_max.lt.2) then
            call writo('ERROR: asking calc_int to integrate with '&
                &//trim(i2str(n_max))//' points. Need at least 2')
            stop
        end if
        
        ! allocate output
        allocate(calc_int(n_max)); calc_int = 0.0_dp
        
        ! calculate integral for all points
        do kd = 2, n_max
            calc_int(kd) = calc_int(kd-1) + &
                &(var(kd)+var(kd-1))/2 * (x(kd)-x(kd-1))
        end do
    end function calc_int
    
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
    ! instead of  the function  itself. Of  course, k should  be lower  than the
    ! degree of the polynomial
    function ext_var(var,var_points,ext_point,deriv_in)
        ! input / output
        real(dp) :: ext_var                                                     ! output
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
        
        
        ! determine the degree of the polynomial ext_var
        pol_deg = size(var)
        
        ! tests
        if (size(var).ne.size(var_points)) then
            call writo('ERROR: in ext_var, the size of the arrays containing &
                &the variable and the points do not match')
            stop
        end if
        if (present(deriv_in)) then 
            deriv = deriv_in
            if (deriv.ge.pol_deg) then                                          ! order of derivative too hight
                call writo('ERROR: asking ext_var for '//trim(i2str(deriv_in))&
                    &//'th derivative, while using a polynomial of degree '&
                    &//trim(i2str(pol_deg)))
                stop
            else if (deriv.lt.0) then
                call writo('ERROR: asking ext_var for a derivative of order '&
                    &//trim(i2str(deriv_in)))
                stop
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
            call writo('ERROR: in ext_var, the solution could not be found.')
            call writo(' -> Are you sure that you are extrapolating independent&
                & points?')
            call writo(' -> The matrix equation to be solved was Ax=b, A = ')
            call print_ar_2(A)
            call writo('    b = ')
            call print_ar_1(var)
            stop
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
    end function ext_var
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    ! VERSION with arr_1 and arr_2 in three coords.
    subroutine arr_mult_3_3(arr_1,arr_2,arr_3,deriv)
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: arr_2(1:,1:,0:,0:,0:)
        real(dp), intent(out) :: arr_3(1:,1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        integer :: r, m, n                                                      ! current degree of norm., pol., tor. derivatives
        integer :: kd                                                           ! normal counter
        
        ! tests
        if (size(arr_1).ne.size(arr_2)) then
            call writo('ERROR: In arr_mult, arrays 1 and 2 need to have the &
                &same size')
            stop
        end if
        if (size(arr_1,1).ne.size(arr_3,1) .or. size(arr_1,2).ne.size(arr_3,2)) then
            call writo('ERROR: In arr_mult, arrays 1 and 2 need to have the &
                &same size as the resulting array 3')
            stop
        end if
        
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
    end subroutine
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    ! VERSION with arr_1 in three coords and arr_2 only in the flux coord.
    subroutine arr_mult_3_1(arr_1,arr_2,arr_3,deriv)
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,0:,0:,0:)
        real(dp), intent(in) :: arr_2(1:,0:)
        real(dp), intent(out) :: arr_3(1:,1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac                                                     ! binomial factor for norm. sum
        integer :: r                                                            ! current degree of norm. derivative
        integer :: kd                                                           ! normal counter
        
        ! tests
        if (size(arr_1,2).ne.size(arr_2,1)) then
            call writo('ERROR: In arr_mult, arrays 1 and 2 need to have the &
                &same size in the flux coord.')
            stop
        end if
        if (size(arr_1,1).ne.size(arr_3,1) .or. size(arr_1,2).ne.size(arr_3,2)) then
            call writo('ERROR: In arr_mult, array 1 needs to have the same &
                &size as the resulting array 3')
            stop
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
    end subroutine
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    ! VERSION with arr_1 and arr_2 only in the flux coord.
    subroutine arr_mult_1_1(arr_1,arr_2,arr_3,deriv)
        ! input / output
        real(dp), intent(in) :: arr_1(1:,0:)
        real(dp), intent(in) :: arr_2(1:,0:)
        real(dp), intent(out) :: arr_3(1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac                                                     ! binomial factor for norm. sum
        integer :: r                                                            ! current degree of norm. derivative
        
        ! tests
        if (size(arr_1,1).ne.size(arr_2,1)) then
            call writo('ERROR: In arr_mult, arrays 1 and 2 need to have the &
                &same size in the flux coord.')
            stop
        end if
        if (size(arr_1,1).ne.size(arr_3,1)) then
            call writo('ERROR: In arr_mult, array 1 needs to have the same &
                &size as the resulting array 3')
            stop
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
    end subroutine

    ! calculate determinant of a matrix which is defined on a 2D grid
    ! the size should be small, as the direct formula is used.
    recursive function det_2D(A) result (detA)
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)
        real(dp), allocatable :: detA(:,:)
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n                                                            ! holds size of A
        real(dp) :: sgn                                                         ! holds either plus or minus one
        integer, allocatable :: slct(:)                                         ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        
        ! tests
        if (size(A,3).ne.size(A,4)) then
            call writo('ERROR: In det(A), the matrix A has to be square')
            stop
        end if
        if (size(A,3).gt.4) then
            call writo('ERROR: det(A) should only be used for small matrices')
            stop
        end if
        
        ! intialize
        n = size(A,3)
        sgn = 1.0_dp
        allocate (idx(n))
        idx = [(id,id=1,n)]
        allocate (slct(n))
        slct = 1
        allocate (detA(size(A,1),size(A,2)))
        detA = 0.0_dp
        
        if (n.eq.1) then                                                        ! shouldn't be used
            detA = A(:,:,1,1)
        else if (n.eq.2) then
            detA = A(:,:,1,1)*A(:,:,2,2)-A(:,:,1,2)*A(:,:,2,1)
        else
            do id = 1, n
                slct(id) = 0
                detA = detA + &
                    &A(:,:,1,id)*sgn*det_2D(A(:,:,2:n,pack(idx,slct.gt.0)))
                sgn = -sgn
                slct(id) = 1
            end do
        end if
    end function

    ! calculate determinant of a constant matrix
    ! (adapted from http://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/)
    ! ¡¡¡WARNING: the matrix A is changed on output!!!
    function det_0D(A)
        ! input / output
        real(dp), intent(inout) :: A(:,:)
        real(dp) :: det_0D
        
        ! local variables
        integer :: n 
        real(dp) :: sgn
        integer :: id, info
        integer, allocatable :: ipiv(:)
        
        ! tests
        if (size(A,1).ne.size(A,2)) then
            call writo('ERROR: In det(A), the matrix A has to be square')
            stop
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
    end function det_0D
    
    ! calculate inverse of  square matrix A which has  elements depending on
    ! 2D mesh (first two coords)
    ! this method uses direct inversion using  Cramer's rule, since the matrix A
    ! is supposed to be  very small (i.e. 3x3 or 1x1) and  since the inverse has
    ! to be calculated at  each of the points in the 2D mesh,  this can be quite
    ! fast
    function inv_2D(A)
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)
        real(dp), allocatable :: inv_2D(:,:,:,:)
        
        ! local variables
        real(dp), allocatable :: detA(:,:)
        integer :: id, kd                                                       ! counters
        integer :: n
        integer, allocatable :: slct(:,:)                                       ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        
        ! tests
        if (size(A,3).ne.size(A,4)) then
            call writo('ERROR: In inv_2D, the matrix has to be square')
            stop
        end if
        
        ! initializing
        n = size(A,3)
        allocate(detA(size(A,1),size(A,2)))
        detA = 0.0_dp
        allocate (idx(n))
        idx = [(id,id=1,n)]
        allocate (slct(n,2))
        slct = 1
        allocate(inv_2D(size(A,1),size(A,2),n,n))
        
        ! calculate determinant
        detA = det(A)
        
        ! calculate cofactor matrix, replacing original elements
        do kd = 1,n
            do id = 1,n
                slct(kd,1) = 0
                slct(id,2) = 0
                inv_2D(:,:,id,kd) = (-1.0_dp)**(id+kd)*det(A(:,:,&
                    &pack(idx,slct(:,1).gt.0),pack(idx,slct(:,2).gt.0)))
                slct(kd,1) = 1
                slct(id,2) = 1
            end do
        end do
        
        ! divide by determinant
        do kd = 1,n
            do id = 1,n
                inv_2D(:,:,id,kd) = inv_2D(:,:,id,kd) / detA
            end do
        end do
    end function inv_2D
    
    ! calculate inverse of a constant square matrix
    function inv_0D(A)
        ! input / output
        real(dp), intent(in) :: A(:,:)
        real(dp), allocatable :: inv_0D(:,:)
        
        ! local variables
        integer :: n 
        integer :: info
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: w(:)                                           ! work variable
        
        ! tests
        if (size(A,1).ne.size(A,1)) then
            call writo('ERROR: In inv_0D(A), the matrix A has to be square')
            stop
        end if
        
        ! initialize
        n = size(A,1)
        allocate(ipiv(n))
        ipiv = 0
        allocate(w(n))
        allocate(inv_0D(n,n))
        
        inv_0D = A
        
        call dgetrf(n, n, inv_0D, n, ipiv, info)                                ! LU factorization
        
        call dgetri(n, inv_0D, n, ipiv, w, n, info)                             ! inverse of LU
    end function inv_0D

    

 
    ! finds the zero of a function using Newton-Rhapson iteration
    function zero_NR(fun,dfun,guess)
        use num_vars, only: max_it_NR, tol_NR
        
        ! input / output
        real(dp) :: zero_NR
        real(dp) :: guess
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
        
        zero_NR = guess
        
        NR: do jd = 1,max_it_NR
            ! correction to theta_NR
            corr = -fun(zero_NR)/dfun(zero_NR)
            zero_NR = zero_NR + corr
            
            ! check for convergence
            ! CHANGE THIS TO RELATIVE
            if (abs(corr).lt.tol_NR) then
                exit
            else if (jd .eq. max_it_NR) then
                call writo('ERROR: Newton-Rhapson method not converged after '&
                    &//trim(i2str(jd))//' iterations. Try increasing max_it_NR&
                    & in input file?')
                call writo('(the residual was equal to '//&
                    &trim(r2str(corr)))
                zero_NR = 0.0_dp
                stop
            end if
        end do NR
    end function zero_NR

    ! calculates the radial mesh point for VMEC variables
    subroutine find_VMEC_norm_coord(fun_name,pt_r,FM,kd,n_max)
        ! input / output
        real(dp) :: pt_r
        logical :: FM
        integer :: kd
        character(len=*) :: fun_name
        integer, intent(in) :: n_max
        
        ! local variables
        real(dp), parameter :: prec = 100*epsilon(1.0_dp)                       ! precision which we want
        
        ! determine whether to use FM or HM and which radial value, test
        if (mod(pt_r,1.0_dp).lt.prec) then                                     ! FM
            FM = .true.
            kd = nint(pt_r)
            if (kd.lt.1 .or. kd.gt.n_max) then
                call writo('ERROR: FM r-range for '//trim(fun_name)//' is &
                    &discrete values between 1 and '//trim(i2str(n_max))//&
                    &', yet asking for '//trim(i2str(kd)))
                stop
            end if
        else if (mod(pt_r,0.5_dp).lt.prec) then                                ! HM
            FM = .false.
            kd = nint(pt_r+0.5_dp)
            if (kd.lt.2 .or. kd.gt.n_max) then
                call writo('ERROR: HM r-range for '//trim(fun_name)//' is &
                    &discrete values between 2 and '//trim(i2str(n_max))//&
                    &', yet asking for '//trim(i2str(kd)))
                stop
            end if
        else 
            call writo('ERROR: for '//trim(fun_name)//', use integer values &
                &for FM quantities or half-integer values for HM. The &
                &provided value deviated too much from this.')
            stop
        end if 
    end subroutine
end module utilities
