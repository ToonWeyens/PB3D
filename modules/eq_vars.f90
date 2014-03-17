!-------------------------------------------------------
!   Variables,  subroutines  and functions  that  have  to do  with  equilibrium
!   quantities and the mesh used in the calculations
!-------------------------------------------------------
module eq_vars
    use num_vars, only: dp, pi
    use output_ops, only: writo, lvl_ud, print_ar_2
    use str_ops, only: r2strt, i2str, r2str

    implicit none
    private
    public eqd_mesh, tor_mesh, pol_mesh, calc_RZl, theta_B, &
        &theta, zeta, n_zeta, R, Z, lam, min_zeta, max_zeta 

    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative)
    real(dp), allocatable :: R(:,:,:,:), Z(:,:,:,:), lam(:,:,:,:)
    real(dp), allocatable :: theta(:,:), zeta(:,:)                                ! grid points
    real(dp) :: min_zeta, max_zeta
    integer :: n_zeta

contains
    ! calculate the toroidal mesh, filling the global variable zeta. This is done
    ! for every flux surface. 
    ! ¡¡¡SHOULD BE DONE ADAPTIVELY!!!
    ! ¡¡¡NEED A CHECK FOR THE FINENESS!!!
    subroutine tor_mesh
        use VMEC_vars, only: n_r
        
        ! local variables
        integer :: kd
        
        if (allocated(zeta)) deallocate(zeta)
        allocate(zeta(n_zeta, n_r)); zeta = 0.0_dp
        do kd = 1,n_r
            zeta(:,kd) = eqd_mesh(n_zeta, min_zeta, max_zeta)
        end do
    end subroutine tor_mesh

    ! calculate the poloidal mesh, in which the magnetic field line is straight,
    ! for a particular alpha
    subroutine pol_mesh(alpha)
        use VMEC_vars, only: n_r
        
        ! input / output
        real(dp) :: alpha
        real(dp) :: zeta_in(n_r)                                                ! to avoid annoying warnings on array creation
        
        ! local variables
        integer :: jd
           
        if (allocated(theta)) deallocate(theta)
        allocate(theta(n_zeta, n_r)); theta = 0.0_dp
        do jd = 1,n_zeta
            zeta_in = zeta(jd,:)
            theta(jd,:) = theta_B(alpha,zeta_in)
        end do
    end subroutine pol_mesh
    
    ! calculate mesh of equidistant points
    function eqd_mesh(n_ang, min_ang, max_ang)
        ! input and output
        real(dp), allocatable :: eqd_mesh(:)
        real(dp), intent(in) :: min_ang, max_ang
        integer, intent(in) :: n_ang
        
        ! local variables
        integer :: id
        real(dp) :: delta_ang
        
        ! test some values
        if (n_ang.lt.1) then
            call writo('ERROR: the angular array has to have a length of &
                &at least 1')
            stop
        end if
        if (min_ang.gt.max_ang) then
            call writo('ERROR: the minimum angle has to be smaller than &
                &the maximum angle')
            stop
        end if
        
        ! initialize output vector
        allocate(eqd_mesh(n_ang)); eqd_mesh = 0.0_dp
        ! There are (n_ang-1) pieces in the total interval but the last one 
        ! is not needed as the functions are all periodic
        
        delta_ang = (max_ang-min_ang)/(n_ang)
        
        eqd_mesh(1) = min_ang
        do id = 2,n_ang
            eqd_mesh(id) = eqd_mesh(id-1) + delta_ang
        end do
    end function eqd_mesh
    
    ! transform half-mesh to full-mesh quantities
    function h2f(var)
        use VMEC_vars, only: n_r
        real(dp) :: h2f(n_zeta,n_zeta,n_r,4)
        real(dp), intent(in) :: var(n_zeta,n_zeta,n_r,4)
        
        ! interpolate the inner quantities
        h2f(:,:,2:n_r-1,:) = 0.5_dp*(var(:,:,2:n_r-1,:) + var(:,:,3:n_r,:))          ! Only last values must come from B.C.
        h2f(:,:,n_r,:) = 2*var(:,:,n_r,:) - h2f(:,:,n_r-1,:)                       ! Extrapolate first, last points
        h2f(:,:,1,:) = 2*var(:,:,2,:)-h2f(:,:,2,:)
    end function

    ! calculate the coordinates R  and Z and l in real  space from their Fourier
    ! decomposition using the grid points currently stored in the variables 
    ! theta, zeta
    subroutine calc_RZl
        use fourier_ops, only: mesh_cs, f2r
        use VMEC_vars, only: R_c, R_s, Z_c, Z_s, l_c, l_s, n_r, rmax_surf, &
            &rmin_surf, zmax_surf, ntor, mpol
        
        ! local variables
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        real(dp) :: delta_r                                                     ! normal step size
        integer :: id, jd, kd
        
        ! deallocate if allocated
        if (allocated(R)) deallocate(R)
        if (allocated(Z)) deallocate(Z)
        if (allocated(lam)) deallocate(lam)
        ! reallocate
        allocate(R(n_zeta,n_zeta,n_r,4)); R = 0.0_dp
        allocate(Z(n_zeta,n_zeta,n_r,4)); Z = 0.0_dp
        allocate(lam(n_zeta,n_zeta,n_r,4)); lam = 0.0_dp
        
        ! do calculations for all angular points
        tor: do jd = 1, n_zeta
            pol: do id = 1,n_zeta
                ! calculate the variables R, Z, and their angular derivatives for all normal points and current angular point
                do kd = 1, n_r
                    ! calculate the (co)sines at the current mesh points
                    cs = mesh_cs(mpol,ntor,theta(id,kd),zeta(jd,kd))
                    
                    R(id,jd,kd,:) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor)
                    Z(id,jd,kd,:) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor)
                    lam(id,jd,kd,:) = f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor)
                end do
                
                ! numerically calculate normal derivatives at the currrent angular points
                ! first normal point
                ! FULL MESH quantities R and Z
                delta_r = 1.0/(n_r-1)                                           ! step size between first points
                R(id,jd,1,2) = (R(id,jd,2,1)-R(id,jd,1,1))/delta_r              ! forward difference
                Z(id,jd,1,2) = (Z(id,jd,2,1)-Z(id,jd,1,1))/delta_r              ! forward difference
                do kd = 3, n_r
                    delta_r = 2.0/(n_r-1)                                       ! intermediate step size: centered differences
                    R(id,jd,kd-1,2) = (R(id,jd,kd,1)-R(id,jd,kd-2,1))/delta_r   ! centered difference
                    Z(id,jd,kd-1,2) = (Z(id,jd,kd,1)-Z(id,jd,kd-2,1))/delta_r   ! centered difference
                end do
                ! last normal point
                delta_r = 1.0/(n_r-1)                                           ! step size between last points
                R(id,jd,n_r,2) = (R(id,jd,n_r,1)-R(id,jd,n_r-1,1))/delta_r
                Z(id,jd,n_r,2) = (Z(id,jd,n_r,1)-Z(id,jd,n_r-1,1))/delta_r
                ! HALF MESH quantity lambda
                delta_r = 1.0/(n_r-1)                                           ! step size between first points
                lam(id,jd,2,2) = (lam(id,jd,3,1)-lam(id,jd,2,1))/delta_r              ! forward difference
                do kd = 4, n_r
                    delta_r = 2.0/(n_r-1)                                       ! intermediate step size: centered differences
                    lam(id,jd,kd-1,2) = (lam(id,jd,kd,1)-lam(id,jd,kd-2,1))&
                        &/delta_r                                               ! centered difference
                end do
                ! last normal point
                delta_r = 1.0/(n_r-1)                                           ! step size between last points
                lam(id,jd,n_r,2) = (lam(id,jd,n_r,1)-lam(id,jd,n_r-1,1))/delta_r
            end do pol
        end do tor
        
        ! output a message if the found R  and Z values are not within the VMEC-
        ! provided bounds
        call writo('Checking the bounds of R')
        call lvl_ud(1)
        call within_bounds(R(:,:,:,1),rmin_surf,rmax_surf)
        call lvl_ud(-1)
        call writo('Checking the bounds of Z')
        call lvl_ud(1)
        call within_bounds(Z(:,:,:,1),-zmax_surf,zmax_surf)
        call lvl_ud(-1)
        ! ̉¿¿¿¿¿ IS THERE SOME CRITERION FOR THE MAXIMUM AND / OR MINIMUM OF LAMBDA ???
    contains
        ! display whether the calculated table for R or Z fits within the limits
        ! outputted by VMEC
        subroutine within_bounds(var,min_VMEC,max_VMEC)
            ! input and output
            real(dp), intent(in) :: var(:,:,:)
            real(dp), intent(in) :: min_VMEC, max_VMEC
            
            ! local variables
            real(dp) :: margin
            real(dp) :: min_frac , max_frac
            
            margin = 1.0E-2_dp                                                  ! 0.1% for comparison 
            
            ! minimum and maximum of variable
            min_frac = 2*abs((minval(var)-min_VMEC)/(minval(var)+min_VMEC))
            max_frac = 2*abs((maxval(var)-max_VMEC)/(maxval(var)+max_VMEC))
            
            if (min_frac.lt.margin .and. max_frac.lt.margin) then
                return
            else if (min_frac.gt.margin) then                                   ! too low minimum
                call writo('WARNING: difference between minimum of variable in &
                    &real angular space and VMEC provided minimum is '//&
                    &trim(r2strt(100*min_frac))//'%...')
                call writo(' -> Maybe use an even number of poloidal points, &
                    &or more angular mesh points in general?')
            else if (max_frac.ge.margin) then                                   ! too high maximum
                call writo('WARNING: difference between maximum of variable in &
                    &real angular space and VMEC provided maximum is '//&
                    &trim(r2strt(100*max_frac))//'%...')
                call writo(' -> Maybe use more angular mesh points')
            end if
        end subroutine
    end subroutine calc_RZl
    
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
        real(dp), intent(in) :: alpha_in, zeta_in(n_r)                          ! alpha, zeta
        real(dp), optional :: theta_in(n_r)                                     ! optional input (guess) for theta
        
        ! local variables
        integer :: jd,kd
        real(dp) :: lam(4)                                                      ! lambda, interpolated from current and next mesh points
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)
        real(dp) :: f, f_theta
        real(dp) :: theta_NR                                                    ! temporary solution for a given r, iteration
        
        ! for all normal points
        norm: do kd = 1, n_r
            ! initialization
            f = 0.0_dp
            f_theta = 0.0_dp
            
            ! if first guess for theta is given
            if (present(theta_in)) then
                theta_NR = theta_in(kd)
            else if (kd.eq.1) then
                theta_NR = pi                                                   ! take pi, because it is in the middle of 0...2pi
            else                                                                ! take solution for previous flux surface
                theta_NR = theta_B(kd-1)
            end if
            
            ! Newton-Rhapson loop
            NR: do jd = 1,max_it_NR
                ! transform lambda from Fourier space to real space
                ! calculate the (co)sines
                cs = mesh_cs(mpol,ntor,theta_NR,zeta_in(kd))
                
                ! calculate lambda and angular derivatives, converted to FM
                if (kd.eq.1) then                                               ! first point not defined on HM -> extrapolate
                    lam = 3./2. * f2r(l_c(:,:,2),l_s(:,:,2),cs,mpol,ntor) - &
                        &1./2. * f2r(l_c(:,:,3),l_s(:,:,3),cs,mpol,ntor)
                else if (kd.eq.n_r) then                                        ! last point -> extrapolate as well
                    lam = 3./2. * f2r(l_c(:,:,n_r),l_s(:,:,n_r),cs,mpol,ntor)-&
                        &1./2. * f2r(l_c(:,:,n_r-1),l_s(:,:,n_r-1),cs,mpol,ntor)
                else                                                            ! intermediate points -> interpolate
                    lam = 1./2. * f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor) + &
                        &1./2. * f2r(l_c(:,:,kd+1),l_s(:,:,kd+1),cs,mpol,ntor)
                end if
                
                ! calculate the factors f and f_theta
                f = zeta_in(kd) - (theta_NR+lam(1))/iotaf(kd) - alpha_in
                f_theta = -(1.0_dp + lam(3))/iotaf(kd)
                
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
        end do norm
    end function theta_B
end module eq_vars
