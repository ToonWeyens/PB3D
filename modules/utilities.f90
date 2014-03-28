!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module utilities
    use num_vars, only: dp
    use output_ops, only: writo
    use str_ops, only: i2str, r2strt, r2str
    
    implicit none
    private
    public zero_NR

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

end module utilities
