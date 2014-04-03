!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module utilities
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_1, print_ar_2
    use str_ops, only: i2str, r2strt, r2str
    
    implicit none
    private
    public zero_NR, ext_var, calc_norm_deriv, f2h, h2f, mat_mult, det, &
        &matvec_mult, calc_int

contains
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
        use var_ops, only: matvec_mult
        
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

    ! transform half-mesh to full-mesh quantities
    ! Adapted from SIESTA routines
    function h2f(var)
        ! input / output
        real(dp), allocatable :: h2f(:)
        real(dp), intent(in) :: var(:)
        
        ! local variables
        integer :: n_max
        integer :: kd
        
        n_max = size(var)
        allocate(h2f(n_max))
        
        ! interpolate the inner quantities
        h2f(2:n_max-1) = 0.5_dp*(var(2:n_max-1) + var(3:n_max))
        ! extrapolate the outer quantities using 3 points
        h2f(n_max) = ext_var([(var(kd),kd=n_max-2,n_max)],&
           &[(1.0_dp*(kd-0.5)/(n_max-1),kd=n_max-3,n_max-1)],1.0_dp,0)
        h2f(1) = ext_var([(var(kd),kd=2,4)],&
           &[(1.0_dp*(kd-0.5)/(n_max-1),kd=1,3)],0.0_dp,0)
    end function

    ! transform  half-mesh to  full-mesh  quantities based  on  h2f. The  result
    ! is  reversible   (h2f*f2h  =   h)  only  if   the  second   derivative  of
    ! the  function  is  zero  at  all  points,  because  the  condition  h_i  =
    ! 1/4*(h_(i-1)+2*h_i+h_(i+1)) can  be derived. This is  logical, because the
    ! linear interpolation is only exact for first order polynomials.
    function f2h(var)
        ! input / output
        real(dp), allocatable :: f2h(:)
        real(dp), intent(in) :: var(:)
        
        ! local variables
        integer :: n_max
        
        n_max = size(var)
        allocate(f2h(n_max))
        
        ! interpolate the inner quantities
        f2h(2:n_max) = 0.5_dp*(var(1:n_max-1) + var(2:n_max))                   ! Only last values must come from B.C.
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
    
    ! calculates normal  derivatives of a  FM or  HM quantity. The  inverse step
    ! size  has to  be specified,  normalized  so that  the whole  range of  the
    ! variable is between 0 and 1 (both inclusive)
    ! Both input and output can be FM or HM:
    !   +-----------------------------------+
    !   | I\O |      HM      |      FM      |
    !   | ----------------------------------|
    !   | HM  | centr. diff. |  (i+1)-(i)   |
    !   | ----------------------------------|
    !   | FM  |   (i)-(i-1)  | centr. diff. |
    !   +-----------------------------------+
    function calc_norm_deriv(var,inv_step,FM_i,FM_o)
        ! input / output
        real(dp), allocatable :: calc_norm_deriv(:)                             ! output variable
        real(dp), intent(in) :: var(:)                                          ! input variable
        logical, intent(in) :: FM_i, FM_o                                       ! whether or not Full Mesh for input and output (true: FM, false: HM)
        real(dp), intent(in) :: inv_step                                        ! inverse of step size in var
        
        ! local variables
        real(dp), allocatable :: varout(:)                                      ! temporary variable to hold output
        integer :: kd                                                           ! counters
        integer :: n_max                                                        ! maximum index of array
        
        n_max = size(var)
        allocate(calc_norm_deriv(n_max))
        allocate(varout(n_max))
        
        if (FM_i) then                                                          ! FM input
            if (FM_o) then                                                      ! FM output
                ! (2,2) FM_i, FM_o: CENTRAL DIFFERENCES
                ! first normal point: extrapolate using 4 points
                varout(1) = ext_var([(var(kd),kd=1,4)],&
                    &[(kd/inv_step,kd=0,3)],0.0_dp,1)
                ! internal points
                do kd = 3, n_max
                    varout(kd-1) = (var(kd)-var(kd-2)) * inv_step / 2.0_dp      ! centered difference
                end do
                ! last normal point: extrapolate using 4 points
                varout(n_max) = ext_var([(var(kd),kd=n_max-3,n_max)],&
                    &[(kd/inv_step,kd=n_max-4,n_max-1)],1.0_dp,1)
            else                                                                ! HM output
                ! (2,1) FM_i, HM_o: (i)-(i-1)
                varout(1) = 0.0_dp
                do kd = 2, n_max
                    varout(kd) = (var(kd)-var(kd-1)) * inv_step                 ! centered difference
                end do
                !varout(2) = ext_var([(var(jd),jd=1,10)],&
                    !&[(jd/inv_step,jd=0,9)],0.5_dp/inv_step,1)
                !varout(3) = ext_var([(var(jd),jd=1,10)],&
                    !&[(jd/inv_step,jd=0,9)],1.5_dp/inv_step,1)
                !varout(4) = ext_var([(var(jd),jd=1,10)],&
                    !&[(jd/inv_step,jd=0,9)],2.5_dp/inv_step,1)
                !varout(5) = ext_var([(var(jd),jd=1,10)],&
                    !&[(jd/inv_step,jd=0,9)],3.5_dp/inv_step,1)
            end if
        else                                                                    ! HM input
            if (FM_o) then                                                      ! FM output
                ! (1,2) HM_i, FM_o: (i+1)-(i)
                ! first normal point: extrapolate using 4 points
                varout(1) = ext_var([(var(kd),kd=2,5)],&
                    &[((kd-0.5)/inv_step,kd=1,4)],0.0_dp,1)
                ! internal points
                do kd = 2, n_max-1
                    varout(kd) = (var(kd+1)-var(kd)) * inv_step                 ! centered difference
                end do
                ! last normal point: extrapolate using 4 points
                varout(n_max) = ext_var([(var(kd),kd=n_max-4,n_max)],&
                    &[((kd-0.5)/inv_step,kd=n_max-5,n_max-1)],1.0_dp,1)
            else                                                                ! HM output
                ! (1,1) HM_i, HM_o: CENTRAL DIFFERENCES
                ! first normal point: extrapolate using 4 points
                varout(1) = 0.0_dp
                varout(2) = ext_var([(var(kd),kd=2,5)],&
                    &[((kd-0.5)/inv_step,kd=1,4)],0.5_dp/inv_step,1)
                ! internal points
                do kd = 4, n_max
                    varout(kd-1) = (var(kd)-var(kd-2)) * inv_step / 2.0_dp      ! centered difference
                end do
                ! last normal point: extrapolate using 4 points
                varout(n_max) = ext_var([(var(kd),kd=n_max-4,n_max)],&
                    &[((kd-0.5)/inv_step,kd=n_max-5,n_max-1)],&
                    &(n_max-1.5_dp)/inv_step,1)
            end if
        end if
        
        calc_norm_deriv = varout
    end function calc_norm_deriv
    
    ! multipy two matrices
    function mat_mult(A,B)
        real(dp) :: A(:,:), B(:,:)
        real(dp), allocatable :: mat_mult(:,:)
        integer :: id, jd, kd
        
        if (size(A,2).ne.size(B,1)) then
            call writo('ERROR: matrices A and B not compatible')
            stop
        end if
        
        allocate(mat_mult(size(A,1),size(B,2)))
        mat_mult = 0.0_dp
        do jd = 1,size(B,2)
            do id = 1,size(A,1)
                do kd = 1, size(A,2)
                    mat_mult(id,jd) = mat_mult(id,jd) + A(id,kd)*B(kd,jd)
                end do
            end do
        end do
    end function mat_mult
    
    ! multiply a matrix with a vector
    function matvec_mult(A,b)
        real(dp) :: A(:,:), b(:)
        real(dp), allocatable :: matvec_mult(:)
        integer :: id, kd
        
        if (size(A,2).ne.size(b)) then
            call writo('ERROR: Matrix A and vector b not compatible')
            stop
        end if
        
        allocate(matvec_mult(size(A,1)))
        matvec_mult = 0.0_dp
        do id = 1, size(A,1)
            do kd = 1, size(A,2)
                matvec_mult(id) = matvec_mult(id) + A(id,kd)*b(kd)
            end do
        end do
    end function matvec_mult

    ! subtract two matrices
    function mat_sub(A,B)
        real(dp) :: A(:,:), B(:,:)
        real(dp), allocatable :: mat_sub(:,:)
        integer :: id, jd
        
        if (size(A,1).ne.size(B,1) .or. size(A,2).ne.size(B,2)) then
            call writo('ERROR: matrices A and B not compatible')
            stop
        end if
        
        allocate(mat_sub(size(A,1),size(A,2)))
        mat_sub = 0.0_dp
        do jd = 1,size(B,2)
            do id = 1,size(A,1)
                mat_sub(id,jd) = A(id,jd)-B(id,jd)
            end do
        end do
    end function mat_sub

    ! calculate determinant of a matrix
    ! (adapted from http://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/)
    real(dp) function det(N, mat)
        integer, intent(in) :: N 
        real(dp), intent(inout) :: mat(:,:)
        integer :: i, info
        integer, allocatable :: ipiv(:)
        real(dp) :: sgn
        
        allocate(ipiv(N))
        
        ipiv = 0
        
        call dgetrf(N, N, mat, N, ipiv, info)
        
        det = 1.0_dp
        do i = 1, N
            det = det*mat(i, i)
        end do
        
        sgn = 1.0_dp
        do i = 1, N
            if(ipiv(i) /= i) then
                sgn = -sgn
            end if
        end do
        det = sgn*det   
    end function det
end module utilities
