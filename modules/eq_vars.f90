!-------------------------------------------------------
!   Variables,  subroutines  and functions  that  have  to do  with  equilibrium
!   quantities and the mesh used in the calculations
!   These are called in the subroutine calc_eq in the module eq_ops
!-------------------------------------------------------
module eq_vars
    use num_vars, only: dp, pi
    use output_ops, only: writo, lvl_ud, print_ar_2, print_ar_1, write_out
    use str_ops, only: r2strt, i2str, r2str

    implicit none
    private
    public eqd_mesh, calc_mesh, calc_RZl, ang_B, calc_flux_q, h2f, &
        &check_mesh, &
        &theta, zeta, n_par, R, Z, lam, min_par, max_par, q_saf, &
        &flux_p, flux_t

    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative)
    real(dp), allocatable :: R(:,:,:), Z(:,:,:), lam(:,:,:)
    real(dp), allocatable :: theta(:,:), zeta(:,:)                              ! grid points
    real(dp), allocatable :: q_saf(:,:), flux_p(:,:), flux_t(:,:)               ! safety factor, pol. flux, tor. flux, and normal derivative
    real(dp) :: min_par, max_par
    integer :: n_par

contains
    ! calculate the angular mesh, filling the global variables (theta, zeta). 
    ! This is done for every flux surface. 
    ! The variable  theta_var_along_B determines  whether theta  is used  as the
    ! base variable. If .false., zeta is used.
    ! ¡¡¡SHOULD BE DONE ADAPTIVELY!!!
    ! ¡¡¡NEED A CHECK FOR THE FINENESS!!!
    subroutine calc_mesh(alpha)
        use VMEC_vars, only: n_r
        use num_vars, only: theta_var_along_B
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        integer :: kd
        real(dp) :: var_in(n_r)
        
        if (allocated(zeta)) deallocate(zeta)
        allocate(zeta(n_par, n_r)); zeta = 0.0_dp
        if (allocated(theta)) deallocate(theta)
        allocate(theta(n_par, n_r)); theta = 0.0_dp
        
        if (theta_var_along_B) then                                             ! first calculate theta
            do kd = 1,n_r
                theta(:,kd) = eqd_mesh(n_par, min_par, max_par)
            end do
            do kd = 1,n_par
                var_in = theta(kd,:)
                zeta(kd,:) = ang_B(.false.,alpha,var_in)
            end do
        else                                                                    ! first calculate zeta
            do kd = 1,n_r
                zeta(:,kd) = eqd_mesh(n_par, min_par, max_par)
            end do
            do kd = 1,n_par
                var_in = zeta(kd,:)
                theta(kd,:) = ang_B(.true.,alpha,var_in)
            end do
        end if
    end subroutine calc_mesh
    
    ! check  whether   the  straight  field   line  mesh  has   been  calculated
    ! satisfactorily,  by  calculating  again  the zeta's  from  the  calculated
    ! theta's
    subroutine check_mesh(alpha)
        use VMEC_vars, only: n_r
        use num_vars, only: tol_NR, theta_var_along_B, max_str_ln
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        real(dp) :: var_calc(n_par,n_r)
        real(dp) :: var_diff(n_par,n_r)
        real(dp) :: var_in(n_r)
        integer :: id
        character(len=max_str_ln) :: par_ang, dep_ang                           ! parallel angle and dependant angle, for output message
        
        if (theta_var_along_B) then                                             ! calculate theta again from zeta
            do id = 1,n_par
                var_in = zeta(id,:)
                var_calc(id,:) = ang_B(.true.,alpha,var_in)
            end do
            par_ang = 'poloidal'; dep_ang = 'toroidal'
            var_diff = theta - var_calc
        else                                                                    ! calculate zeta again from theta
            do id = 1,n_par
                var_in = theta(id,:)
                var_calc(id,:) = ang_B(.false.,alpha,var_in)
            end do
            par_ang = 'toroidal'; dep_ang = 'poloidal'
            var_diff = zeta - var_calc
        end if
        
        if (maxval(abs(var_diff)).gt.tol_NR*100) then
            call writo('ERROR: Calculating again the '//trim(par_ang)//&
                &' grid points from the '//trim(dep_ang)//' grid points, &
                &along the magnetic fields, does not result in the same &
                &values as initally given.')
            call writo('The maximum error is equal to '//&
                &trim(r2strt(100*maxval(abs(var_diff))))//'%')
            stop
        end if
    end subroutine check_mesh

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
    ! Adapted from SIESTA routines
    function h2f(var,nx,ny)
        use VMEC_vars, only: n_r
        integer, intent(in) :: nx, ny                                           ! size of first two dimensins (third is r)
        real(dp) :: h2f(nx,ny,n_r)
        real(dp), intent(in) :: var(nx,ny,n_r)
        
        ! interpolate the inner quantities
        h2f(:,:,2:n_r-1) = 0.5_dp*(var(:,:,2:n_r-1) + var(:,:,3:n_r))           ! Only last values must come from B.C.
        h2f(:,:,n_r) = 2*var(:,:,n_r) - h2f(:,:,n_r-1)                          ! Extrapolate first, last points
        h2f(:,:,1) = 2*var(:,:,2)-h2f(:,:,2)
    end function

    ! calculate the coordinates R  and Z and l in real  space from their Fourier
    ! decomposition  using the  grid points  currently stored  in the  variables
    ! theta, zeta.
    ! the output is all FULL MESH (FM)
    subroutine calc_RZl
        use fourier_ops, only: mesh_cs, f2r
        use VMEC_vars, only: R_c, R_s, Z_c, Z_s, l_c, l_s, n_r, rmax_surf, &
            &rmin_surf, zmax_surf, ntor, mpol
        
        ! local variables
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        integer :: id, kd
        real(dp) :: tempvar(n_r)
        real(dp) :: l_c_F(0:mpol-1,-ntor:ntor,1:n_r)                            ! FM version of HM l_c
        real(dp) :: l_s_F(0:mpol-1,-ntor:ntor,1:n_r)                            ! FM version of HM l_s
        real(dp) :: lam_H(n_r)
        
        ! deallocate if allocated
        if (allocated(R)) deallocate(R)
        if (allocated(Z)) deallocate(Z)
        if (allocated(lam)) deallocate(lam)
        ! reallocate
        allocate(R(n_par,n_r,4)); R = 0.0_dp
        allocate(Z(n_par,n_r,4)); Z = 0.0_dp
        allocate(lam(n_par,n_r,4)); lam = 0.0_dp
        
        ! convert l_c and l_s to FM
        l_c_F = 0.0_dp; l_c_F = h2f(l_c, mpol, 2*ntor+1)
        l_s_F = 0.0_dp; l_s_F = h2f(l_s, mpol, 2*ntor+1)
        
        ! do calculations for all angular points
        par: do id = 1, n_par                                                  ! parallel: along the magnetic field line
            ! calculate the variables R, Z, and their angular derivatives for all normal points and current angular point
            perp: do kd = 1, n_r                                                ! perpendicular: normal to the flux surfaces
                ! calculate the (co)sines at the current mesh points
                cs = mesh_cs(mpol,ntor,theta(id,kd),zeta(id,kd))
                
                R(id,kd,1) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,1)
                R(id,kd,3) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,3)
                R(id,kd,4) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,4)
                Z(id,kd,1) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,1)
                Z(id,kd,3) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,3)
                Z(id,kd,4) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,4)
                lam(id,kd,1) = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,1)
                lam(id,kd,3) = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,3)
                lam(id,kd,4) = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,4)
                lam_H(kd) = f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor,1)             ! HM needed for the calculation of the normal derivative
            end do perp
            
            ! numerically calculate normal derivatives at the currrent angular points
            tempvar = R(id,:,1)
            R(id,:,2) = calc_norm_deriv(tempvar,.true.)
            tempvar = Z(id,:,1)
            Z(id,:,2) = calc_norm_deriv(tempvar,.true.)
            lam(id,:,2) = calc_norm_deriv(lam_H,.false.)
        end do par
        
        ! output a message if the found R  and Z values are not within the VMEC-
        ! provided bounds
        call writo('Checking the bounds of R')
        call lvl_ud(1)
        call within_bounds(R(:,:,1),rmin_surf,rmax_surf)
        call lvl_ud(-1)
        call writo('Checking the bounds of Z')
        call lvl_ud(1)
        call within_bounds(Z(:,:,1),-zmax_surf,zmax_surf)
        call lvl_ud(-1)
        
    contains
        ! display whether the calculated table for R or Z fits within the limits
        ! outputted by VMEC
        ! Since the calculations are done following a magnetic field line, it is
        ! actually quite  probable that the  minimum or maximum is  not reached.
        ! Therefore, only hard  warnings are given, when the  minimum or maximum
        ! calculated value is grossly under or above the VMEC minimum or maximum
        subroutine within_bounds(var,min_VMEC,max_VMEC)
            ! input and output
            real(dp), intent(in) :: var(:,:)
            real(dp), intent(in) :: min_VMEC, max_VMEC
            
            ! local variables
            real(dp) :: margin
            real(dp) :: min_frac , max_frac
            
            margin = 5.0E-2_dp                                                  ! 1% for comparison 
            
            ! minimum and maximum of variable
            min_frac = 2*(minval(var)-min_VMEC)/abs(minval(var)+min_VMEC)       ! positive if within bounds
            max_frac = 2*(maxval(var)-max_VMEC)/abs(maxval(var)+max_VMEC)       ! positive if out of bounds
            
            if (min_frac.lt.-margin .and. max_frac.gt.margin) then
                return
            else if (min_frac.lt.-margin) then                                  ! too low minimum
                call writo('WARNING: minimum of variable in real angular space &
                    & is lower than VMEC provided minimum by '//&
                    &trim(r2strt(100*min_frac))//'%...')
                call writo(' -> Maybe use an even number of poloidal points, &
                    &or more angular mesh points in general?')
                call writo(' -> Maybe run VMEC with more accuracy?')
            else if (max_frac.gt.margin) then                                   ! too high maximum
                call writo('WARNING: maximum of variable in real angular space &
                    & is greater than VMEC provided minimum by '//&
                    &trim(r2strt(100*max_frac))//'%...')
                call writo(' -> Maybe use more angular mesh points')
                call writo(' -> Maybe run VMEC with more accuracy?')
            end if
        end subroutine
    end subroutine calc_RZl

    ! calculates normal derivatives in FM and HM
    function calc_norm_deriv(var,FM)
        use VMEC_vars, only: n_r
        
        ! input / output
        real(dp) :: calc_norm_deriv(n_r)
        real(dp), intent(in) :: var(n_r)
        logical :: FM                                                           ! whether or not Full Mesh (true: FM, false: HM)
        
        ! local variables
        real(dp) :: delta_r                                                     ! step size
        real(dp) :: varout(n_r)                                                 ! temporary variable to hold output
        integer :: kd
        
        if (FM) then                                                            ! full mesh calculation
            ! first normal point
            delta_r = 1.0/(n_r-1)                                               ! step size between first points
            varout(1) = (var(2)-var(1))/delta_r                                 ! forward difference
            ! internal points
            do kd = 3, n_r
                delta_r = 2.0/(n_r-1)                                           ! intermediate step size: centered differences
                varout(kd-1) = (var(kd)-var(kd-2))/delta_r                      ! centered difference
            end do
            ! last normal point
            delta_r = 1.0/(n_r-1)                                               ! step size between last points
            varout(n_r) = (var(n_r)-var(n_r-1))/delta_r
        else                                                                    ! half mesh calculation
            ! first normal point: choose equal to derivative at second point (linear approx.)
            delta_r = 1.0/(n_r-1)                                               ! step size between first points
            varout(1) = (var(3)-var(2))/delta_r                                 ! centered difference
            ! internal points
            do kd = 2, n_r-1
                delta_r = 1.0/(n_r-1)                                           ! intermediate step size: centered differences
                varout(kd) = (var(kd+1)-var(kd))/delta_r                        ! centered difference
            end do
            ! last normal point: choose equal to derivative at next to last point (linear approx.)
            delta_r = 1.0/(n_r-1)                                               ! step size between last points
            varout(n_r) = (var(n_r)-var(n_r-1))/delta_r
        end if
        
        calc_norm_deriv = varout
    end function calc_norm_deriv
    
    ! Calculates flux quantities in FM
    subroutine calc_flux_q
        use VMEC_vars, only: &
            &iotaf, n_r, phi, phipf
        
        ! local variables
        integer :: kd
        
        ! reallocate
        if (allocated(q_saf)) deallocate(q_saf)
        allocate(q_saf(n_r,2))
        if (allocated(flux_p)) deallocate(flux_p)
        allocate(flux_p(n_r,2))
        if (allocated(flux_t)) deallocate(flux_t)
        allocate(flux_t(n_r,2))
        
        ! safety factor q_saf: invert iotaf
        do kd = 1,n_r
            q_saf(kd,1) = 1/iotaf(kd)
        end do
        q_saf(:,2) = calc_norm_deriv(q_saf(:,1),.true.)
        
        ! toroidal flux: copy from VMEC
        flux_t(:,1) = phi
        flux_t(:,2) = phipf
        
        ! poloidal flux: calculate using iotaf and phi, phipf
        do kd = 1,n_r
            flux_p(kd,1) = iotaf(kd)*phi(kd)
        end do
        flux_p(:,2) = calc_norm_deriv(flux_p(:,1),.true.)
    end subroutine
    
    ! Calculates the  poloidal/toroidal angle theta(zeta)  as a function  of the
    ! toroidal/poloidal angle  zeta/theta following a particular  magnetic field
    ! line alpha.
    ! This is done using a Newton-Rhapson scheme that calculates the zero's of the
    ! function f(theta) = zeta - q (theta + lambda(theta)) - alpha
    ! the logical find_theta determines whether theta is sought or zeta:
    !   find_theta = .true. : look for theta as a function of zeta
    !   find_theta = .false. : look for zeta as a function of theta
    function ang_B(find_theta,alpha_in,input_ang,guess)
        use num_vars, only: max_it_NR, tol_NR
        use VMEC_vars, only: mpol, ntor, n_r, l_c, l_s, iotaf
        use fourier_ops, only: mesh_cs, f2r
        
        ! input / output
        real(dp) :: ang_B(n_r)                                                  ! theta(zeta)/zeta(theta)
        logical :: find_theta                                                   ! whether theta or zeta is sought
        real(dp), intent(in) :: alpha_in                                        ! alpha
        real(dp), intent(in) :: input_ang(n_r)                                  ! the input angle zeta/theta
        real(dp), intent(in), optional :: guess(n_r)                            ! optional input (guess) for theta/zeta
        
        ! local variables
        integer :: jd,kd
        real(dp) :: lam, dlam                                                   ! lambda and poloidal/toroidal derivative
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! factors of cosine and sine for inverse fourier transform
        real(dp) :: f, df                                                       ! function f and poloidal/toroidal derivative
        real(dp) :: ang_NR                                                      ! temporary solution for a given r, iteration
        real(dp) :: l_c_F(0:mpol-1,-ntor:ntor,1:n_r)                            ! FM version of HM l_c
        real(dp) :: l_s_F(0:mpol-1,-ntor:ntor,1:n_r)                            ! FM version of HM l_s
        real(dp) :: alpha_calc                                                  ! calculated alpha, to check with the given alpha
        real(dp) :: corr                                                        ! correction calculated
        
        ! convert l_c and l_s to FM
        l_c_F = 0.0_dp; l_c_F = h2f(l_c, mpol, 2*ntor+1)
        l_s_F = 0.0_dp; l_s_F = h2f(l_s, mpol, 2*ntor+1)
        
        ! for all normal points
        norm: do kd = 1, n_r
            ! initialization
            f = 0.0_dp
            df = 0.0_dp
            
            ! if first guess for theta is given
            if (present(guess)) then
                ang_NR = guess(kd)
            else if (kd.eq.1) then
                ang_NR = pi                                                     ! take pi, because it is in the middle of 0...2pi
            else                                                                ! take solution for previous flux surface
                ang_NR = ang_B(kd-1)
            end if
            
            ! Newton-Rhapson loop
            NR: do jd = 1,max_it_NR
                ! transform lambda from Fourier space to real space
                ! calculate the (co)sines
                if (find_theta) then                                            ! looking for theta
                    cs = mesh_cs(mpol,ntor,ang_NR,input_ang(kd))
                else                                                            ! looking for zeta
                    cs = mesh_cs(mpol,ntor,input_ang(kd),ang_NR)
                end if
                
                ! calculate lambda and angular derivatives, using FM coeff.
                lam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,1)
                if (find_theta) then                                            ! looking for theta
                    dlam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,3)
                else                                                            ! looking for zeta
                    dlam = f2r(l_c_F(:,:,kd),l_s_F(:,:,kd),cs,mpol,ntor,4)
                end if
                
                ! calculate the factors f and f_t
                if (find_theta) then                                            ! looking for theta
                    f = input_ang(kd) - (ang_NR+lam)/iotaf(kd) - alpha_in
                    df = -(1.0_dp + dlam)/iotaf(kd)
                else                                                            ! looking for zeta
                    f = ang_NR - (input_ang(kd)+lam)/iotaf(kd) - alpha_in
                    df = 1.0_dp - dlam/iotaf(kd)
                end if
                
                ! correction to theta_NR
                corr = -f/df
                ang_NR = ang_NR + corr
                
                ! check for convergence
                if (abs(corr).lt.tol_NR) then
                    ang_B(kd) = ang_NR
                    exit
                else if (jd .eq. max_it_NR) then
                    call writo('ERROR: Newton-Rhapson method to find &
                        &theta(zeta) not converged after '//trim(i2str(jd))//&
                        &' iterations. Try increasing max_it_NR in input file?')
                    call writo('(the residual was equal to '//&
                        &trim(r2str(corr)))
                    call writo(' with f = '//trim(r2strt(f))//' and df = '&
                        &//trim(r2strt(df))//')')
                    stop
                end if
            end do NR
            
            ! do a check whether the result is indeed alpha
            if (find_theta) then                                                ! looking for theta
                alpha_calc = input_ang(kd) - (ang_B(kd) + lam)/iotaf(kd)
            else                                                                ! looking for zeta
                alpha_calc = ang_B(kd) - (input_ang(kd) + lam)/iotaf(kd)
            end if
            if (alpha_calc-alpha_in.gt.tol_NR*100) then
                call writo('ERROR: In theta_B, calculating alpha as a check,&
                    & using the theta that is the solution of alpha '&
                    &//trim(r2strt(alpha_in))//', yields alpha_calc that &
                    &deviates '//trim(r2strt(100*abs((alpha_calc-alpha_in))))&
                    &//'% from the original alpha_in')
                stop
            end if
        end do norm
    end function ang_B
end module eq_vars
