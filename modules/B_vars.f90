!-------------------------------------------------------
!   Variables, subroutines and functions that have to do with the magnetic field
!-------------------------------------------------------
module B_vars
    use num_vars, only: dp
    use output_ops, only: writo, print_ar_2
    use var_ops, only: r2str, r2strt, i2str

    implicit none
    private
    public theta_B

contains
    ! Calculates the  poloidal angle theta as  a function of the  toroidal angle
    ! zeta following a particular magnetic field line alpha.
    ! This is done using a Newton-Rhapson scheme that calculates the zero's of the
    ! function f(theta) = zeta - q (theta + lambda(theta)) - alpha
    function theta_B(alpha_in,zeta_in,theta_in)
        use num_vars, only: max_it_NR, tol_NR
        use VMEC_vars, only: mpol, ntor, n_r, l_c, l_s, iotaf
        use fourier_ops, only: mesh_cs, f2r
        
        ! input / output
        real(dp) :: theta_B(n_r)                                                ! theta(zeta)
        real(dp), intent(in) :: alpha_in, zeta_in                               ! alpha, zeta
        real(dp), optional :: theta_in(n_r)                                     ! optional input (guess) for theta
        
        ! local variables
        integer :: jd,kd
        real(dp) :: lam(4)
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)
        real(dp) :: f, f_theta
        real(dp) :: theta_NR
        
        ! for all radial points
        rad: do kd = 1, n_r
            ! initialization
            f = 0.0_dp
            f_theta = 0.0_dp
            
            ! if first guess for theta is given
            if (present(theta_in)) then
                theta_NR = theta_in(kd)
            else if (kd.eq.1) then
                theta_NR = zeta_in                                ! take the provided zeta
            else                                                                ! take solution for previous flux surface
                theta_NR = theta_B(kd-1)
            end if
            
            ! Newton-Rhapson loop
            NR: do jd = 1,max_it_NR
                ! calculate the (co)sines
                cs = mesh_cs(mpol,ntor,theta_NR,zeta_in)
                !write(*,*) 'kd (rad), jd (NR) = ', kd, jd
                !write(*,*) 'mpol, ntor, theta_NR, zeta_in = ', mpol, ntor, &
                    !&theta_NR, zeta_in
                
                ! calculate lambda and angular derivatives
                lam = f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor)
                !write(*,*) 'lambda = ', lam
                
                ! calculate the factors f and f_theta
                f = zeta_in - (theta_NR+lam(1))/iotaf(kd) - alpha_in
                f_theta = -(1.0_dp + lam(3))/iotaf(kd)
                !write(*,*) 'f, f_theta = ', f, f_theta
                
                ! correction to theta_NR
                theta_NR = theta_NR - f/f_theta
                
                ! check for convergence
                if (abs(f/f_theta).lt.tol_NR) then
                    theta_B(kd) = theta_NR
                    exit
                else if (jd .eq. max_it_NR) then
                    call writo('ERROR: Newton-Rhapson method to find &
                        &theta(zeta) not converged after '//trim(i2str(jd))//&
                        &' iterations. Try increasing max_it_NR in input file?')
                    call writo('(the residual was equal to '//&
                        &trim(r2str(f/f_theta)))
                    call writo(' with f = '//trim(r2strt(f))//' and f_theta = '&
                        &//trim(r2strt(f_theta))//')')
                    stop
                end if
            end do NR
        end do rad
    end function theta_B
end module B_vars
