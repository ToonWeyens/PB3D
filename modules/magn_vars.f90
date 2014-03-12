module magn_vars
    use num_vars, only: dp
    implicit none
    private
    public magn_theta

contains
    ! Calculates the  poloidal angle theta as  a function of the  toroidal angle
    ! zeta following a particular magnetic field line alpha.
    ! This is done using a Newton-Rhapson scheme that calculates the zero's of the
    ! function f(theta) = zeta - q (theta + lambda(theta)) - alpha
    function magn_theta(alpha,zeta,lambda)
        use num_vars, only: max_it_NR
        
        ! input / output
        real(dp) :: magn_theta
        real(dp), intent(in) :: alpha, zeta, lambda
        
        ! local variables
        integer :: kd
        
        do kd = 1,max_it_NR
        end do
    end function magn_theta
end module magn_vars
