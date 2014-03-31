!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module utilities
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_1, print_ar_2
    use str_ops, only: i2str, r2strt, r2str
    
    implicit none
    private
    public zero_NR, ext_var

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
end module utilities
