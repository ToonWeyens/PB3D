!------------------------------------------------------------------------------!
!   Variables,  subroutines  and functions  that  have  to do  with            !
!   equilibrium quantities and the mesh used in the calculations               !
!   These are called in the subroutine calc_eq in the module eq_ops            !
!------------------------------------------------------------------------------!
module eq_vars
    use num_vars, only: dp, pi
    use output_ops, only: writo, lvl_ud, print_ar_2, print_ar_1, write_out
    use str_ops, only: r2strt, i2str, r2str

    implicit none
    private
    public eqd_mesh, calc_mesh, calc_RZl, ang_B, calc_flux_q, h2f, f2h, &
        &check_mesh, calc_norm_deriv, &
        &theta, zeta, theta_H, zeta_H, n_par, R, Z, lam_H, min_par, max_par, &
        &q_saf, q_saf_H, flux_p, flux_p_H, flux_t, flux_t_H

    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative)
    real(dp), allocatable :: R(:,:,:), Z(:,:,:)                                 ! R, Z (FM)
    real(dp), allocatable :: lam_H(:,:,:)                                       ! lambda in (HM)
    real(dp), allocatable :: theta(:,:), zeta(:,:)                              ! grid points (n_par, n_r, 2) (FM)
    real(dp), allocatable :: theta_H(:,:), zeta_H(:,:)                          ! grid points (n_par, n_r, 2) (HM)
    real(dp), allocatable :: q_saf(:,:), q_saf_H(:,:)                           ! safety factor (FM) and (HM)
    real(dp), allocatable :: flux_p(:,:), flux_t(:,:), pres(:,:)                ! pol. flux, tor. flux and pressure, and normal derivative (FM)
    real(dp), allocatable :: flux_p_H(:,:), flux_t_H(:,:), pres_H(:,:)          ! pol. flux, tor.flux and pressure (HM)
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
        integer :: jd, kd
        real(dp) :: var_in(n_r)
        
        if (allocated(zeta)) deallocate(zeta)
        allocate(zeta(n_par,n_r)); zeta = 0.0_dp
        if (allocated(zeta_H)) deallocate(zeta_H)
        allocate(zeta_H(n_par,n_r)); zeta_H = 0.0_dp
        if (allocated(theta)) deallocate(theta)
        allocate(theta(n_par,n_r)); theta = 0.0_dp
        if (allocated(theta_H)) deallocate(theta_H)
        allocate(theta_H(n_par,n_r)); theta_H = 0.0_dp
        
        if (theta_var_along_B) then                                             ! first calculate theta
            ! ¡¡¡¡¡¡! TEMPORARILY REPLACED !!!!!
            !do jd = 1,n_r
                !theta(:,jd) = eqd_mesh(n_par, min_par, max_par)
                !theta_H(:,jd) = eqd_mesh(n_par, min_par, max_par)
            !end do
            !do kd = 1,n_par
                !var_in = theta(kd,:)
                !zeta(kd,:) = ang_B(.false.,.true.,alpha,var_in)
                !var_in = theta_H(kd,:)
                !zeta-H(kd,:) = ang_B(.false.,.false.,alpha,var_in)
            !end do
            ! TEMPORARY REPLACEMENT:
            theta = 1.*pi/2
            theta_H = 1.*pi/2
            do jd = 1,n_r
                zeta(:,jd) = eqd_mesh(n_par, min_par, max_par)
                zeta_H(:,jd) = eqd_mesh(n_par, min_par, max_par)
            end do
        else                                                                    ! first calculate zeta
            do jd = 1,n_r
                zeta(:,jd) = eqd_mesh(n_par, min_par, max_par)
                zeta_H(:,jd) = eqd_mesh(n_par, min_par, max_par)
            end do
            do kd = 1,n_par
                var_in = zeta(kd,:)
                theta(kd,:) = ang_B(.true.,.true.,alpha,var_in)
                var_in = zeta(kd,:)
                theta_H(kd,:) = ang_B(.true.,.false.,alpha,var_in)
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
        character(len=max_str_ln) :: par_ang, dep_ang                           ! parallel angle and dependent angle, for output message
        
        if (theta_var_along_B) then                                             ! calculate theta again from zeta
            do id = 1,n_par
                var_in = zeta(id,:)
                var_calc(id,:) = ang_B(.true.,.true.,alpha,var_in)
            end do
            par_ang = 'poloidal'; dep_ang = 'toroidal'
            var_diff = theta(:,:) - var_calc
        else                                                                    ! calculate zeta again from theta
            do id = 1,n_par
                var_in = theta(id,:)
                var_calc(id,:) = ang_B(.false.,.true.,alpha,var_in)
            end do
            par_ang = 'toroidal'; dep_ang = 'poloidal'
            var_diff = zeta(:,:) - var_calc
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
    function h2f(var)
        ! input / output
        real(dp), allocatable :: h2f(:)
        real(dp), intent(in) :: var(:)
        
        ! local variables
        integer :: n_max
        
        n_max = size(var)
        allocate(h2f(n_max))
        
        ! interpolate the inner quantities
        h2f(2:n_max-1) = 0.5_dp*(var(2:n_max-1) + var(3:n_max))                 ! Only last values must come from B.C.
        h2f(n_max) = 2*var(n_max) - h2f(n_max-1)                                ! Extrapolate first, last points
        h2f(1) = 2*var(2)-h2f(2)
    end function

    ! transform  half-mesh to  full-mesh  quantities based  on  h2f. The  result
    ! is  reversible   (h2f*f2h  =   h)  only  if   the  second   derivative  of
    ! the  function  is  zero  at  all  points,  because  the  condition  h_i  =
    ! 1/4*(h_(i-1)+2*h_i+h_(i+1)) can  be derived. This is  logical, because the
    ! linear interpolation is only exact for first order polynomials.
    ! ¡¡¡THEREFORE,  THE  WHOLE  INTERPOLATION SHOULD  BE
    ! SUBSTITUTED BY SOMETHING ELSE, SUCH AS SPLINES!!!
    function f2h(var)
        real(dp), allocatable :: f2h(:)
        real(dp), intent(in) :: var(:)
        
        ! local variables
        integer :: n_max
        
        n_max = size(var)
        allocate(f2h(n_max))
        
        ! interpolate the inner quantities
        f2h(2:n_max) = 0.5_dp*(var(1:n_max-1) + var(2:n_max))                   ! Only last values must come from B.C.
        f2h(1) = 0.0_dp
    end function

    ! calculate the coordinates R  and Z and l in real  space from their Fourier
    ! decomposition  using the  grid points  currently stored  in the  variables
    ! theta, zeta.
    subroutine calc_RZl
        use fourier_ops, only: mesh_cs, f2r
        use VMEC_vars, only: R_c, R_s, Z_c, Z_s, lam_c, lam_s, n_r, &
            &rmax_surf, rmin_surf, zmax_surf, ntor, mpol
        
        ! local variables
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        integer :: id, kd
        
        ! deallocate if allocated
        if (allocated(R)) deallocate(R)
        if (allocated(Z)) deallocate(Z)
        if (allocated(lam_H)) deallocate(lam_H)
        ! reallocate
        allocate(R(n_par,n_r,4)); R = 0.0_dp
        allocate(Z(n_par,n_r,4)); Z = 0.0_dp
        allocate(lam_H(n_par,n_r,4)); lam_H = 0.0_dp
        
        ! do calculations for all angular points
        par: do id = 1, n_par                                                   ! parallel: along the magnetic field line
            ! calculate the  variables R, Z,  and their angular  derivatives for
            ! all normal points and current angular point
            perp: do kd = 1, n_r                                                ! perpendicular: normal to the flux surfaces
                ! FM quantities
                cs = mesh_cs(mpol,ntor,theta(id,kd),zeta(id,kd))
                R(id,kd,1) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,1)
                R(id,kd,3) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,3)
                R(id,kd,4) = f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,4)
                Z(id,kd,1) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,1)
                Z(id,kd,3) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,3)
                Z(id,kd,4) = f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,4)
                
                ! HM quantities
                cs = mesh_cs(mpol,ntor,theta_H(id,kd),zeta_H(id,kd))
                lam_H(id,kd,1) = f2r(lam_c(:,:,kd),lam_s(:,:,kd),cs,mpol,ntor,1)
                lam_H(id,kd,3) = f2r(lam_c(:,:,kd),lam_s(:,:,kd),cs,mpol,ntor,3)
                lam_H(id,kd,4) = f2r(lam_c(:,:,kd),lam_s(:,:,kd),cs,mpol,ntor,4)
            end do perp
            
            ! numerically calculate  normal derivatives at the  currrent angular
            ! points
            R(id,:,2) = calc_norm_deriv(R(id,:,1),dfloat(n_r-1),.true.,.true.)
            Z(id,:,2) = calc_norm_deriv(Z(id,:,1),dfloat(n_r-1),.true.,.true.)
            lam_H(id,:,2) = calc_norm_deriv(lam_H(id,:,1),dfloat(n_r-1),.false.&
                &,.false.)
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

    ! calculates normal derivatives of a FM or HM quantity
    ! Both input and output can be FM or HM:
    !   +-----------------------------------+
    !   | I\O |      HM      |      FM      |
    !   | ----------------------------------|
    !   | HM  | centr. diff. |  (i+1)-(i)   |
    !   | ----------------------------------|
    !   | FM  |   (i)-(i-1)  | centr. diff. |
    !   +-----------------------------------+
    ! ¡¡¡¡ TRY TO IMPROVE THIS BY USING HIGHER ORDER EXTRAPOLATION FOR FIRST, LAST POINTS!!!!
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
                ! first normal point
                varout(1) = (var(2)-var(1)) * inv_step                          ! forward difference
                ! internal points
                do kd = 3, n_max
                    varout(kd-1) = (var(kd)-var(kd-2)) * inv_step / 2.0_dp      ! centered difference
                end do
                ! last normal point
                varout(n_max) = (var(n_max)-var(n_max-1)) * inv_step
            else                                                                ! HM output
                ! (2,1) FM_i, HM_o: (i)-(i-1)
                do kd = 2, n_max
                    varout(kd) = (var(kd)-var(kd-1)) * inv_step                 ! centered difference
                end do
            end if
        else                                                                    ! HM input
            if (FM_o) then                                                      ! FM output
                ! (1,2) HM_i, FM_o: (i+1)-(i)
                ! first normal point: choose equal to derivative at second point (linear approx.)
                varout(1) = (var(3)-var(2)) * inv_step                          ! centered difference
                ! internal points
                do kd = 2, n_max-1
                    varout(kd) = (var(kd+1)-var(kd)) * inv_step                 ! centered difference
                end do
                ! last normal point: choose equal to derivative at next to last point (linear approx.)
                varout(n_max) = (var(n_max)-var(n_max-1)) * inv_step
            else                                                                ! HM output
                ! (1,1) HM_i, HM_o: CENTRAL DIFFERENCES
                ! first normal point
                varout(2) = (var(3)-var(2)) * inv_step                          ! forward difference
                ! internal points
                do kd = 4, n_max
                    varout(kd-1) = (var(kd)-var(kd-2)) * inv_step / 2.0_dp      ! centered difference
                end do
                ! last normal point
                varout(n_max) = (var(n_max)-var(n_max-1)) * inv_step
            end if
        end if
        
        calc_norm_deriv = varout
    end function calc_norm_deriv
    
    ! Calculates flux quantities
    subroutine calc_flux_q
        use VMEC_vars, only: &
            &iotaf, iotah, n_r, phi, phi_r, phi_r_H, presf, presh
        
        ! local variables
        integer :: kd
        
        ! reallocate
        if (allocated(q_saf)) deallocate(q_saf)
        allocate(q_saf(n_r,2))
        if (allocated(q_saf_H)) deallocate(q_saf_H)
        allocate(q_saf_H(n_r,2))
        if (allocated(flux_p)) deallocate(flux_p)
        allocate(flux_p(n_r,2))
        if (allocated(flux_p_H)) deallocate(flux_p_H)
        allocate(flux_p_H(n_r,2))
        if (allocated(flux_t)) deallocate(flux_t)
        allocate(flux_t(n_r,2))
        if (allocated(flux_t_H)) deallocate(flux_t_H)
        allocate(flux_t_H(n_r,2))
        if (allocated(pres)) deallocate(pres)
        allocate(pres(n_r,2))
        if (allocated(pres_H)) deallocate(pres_H)
        allocate(pres_H(n_r,2))
        
        ! safety factor q_saf: invert iota and derivate
        do kd = 1,n_r
            q_saf(kd,1) = 1/iotaf(kd)
        end do
        do kd = 2,n_r
            q_saf_H(kd,1) = 1/iotah(kd)
        end do
        q_saf(:,2) = calc_norm_deriv(q_saf_H(:,1),dfloat(n_r-1),.false.,.true.)
        q_saf_H(:,2) = calc_norm_deriv(q_saf(:,1),dfloat(n_r-1),.true.,.false.)
        
        ! toroidal flux: copy from VMEC
        flux_t(:,1) = phi
        flux_t(:,2) = phi_r
        flux_t_H(:,1) = phi
        flux_t_H(:,2) = phi_r_H
        
        ! poloidal flux: calculate using iotaf and phi, phi_r
        ! THIS IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do kd = 1,n_r
            flux_p(kd,1) = iotaf(kd)*phi(kd)
            flux_p_H(kd,1) = iotah(kd)*phi(kd)
        end do
        ! THIS IS CORRECT:
        flux_p(:,2) = iotaf*phi_r
        flux_p_H(:,2) = iotah*phi_r_H
        
        ! pressure: copy from VMEC and derivate
        pres(:,1) = presf
        pres_H(:,1) = presh
        pres(:,2) = calc_norm_deriv(pres_H(:,1),dfloat(n_r-1),.true.,.false.)
        pres_H(:,2) = calc_norm_deriv(pres(:,1),dfloat(n_r-1),.false.,.true.)
    end subroutine

    ! Calculates the  poloidal/toroidal angle theta(zeta)  as a function  of the
    ! toroidal/poloidal angle  zeta/theta following a particular  magnetic field
    ! line alpha.
    ! This is done using a Newton-Rhapson scheme that calculates the zero's of the
    ! function f(theta) = zeta - q (theta + lambda(theta)) - alpha
    ! the logical find_theta determines whether theta is sought or zeta:
    !   find_theta = .true. : look for theta as a function of zeta
    !   find_theta = .false. : look for zeta as a function of theta
    function ang_B(find_theta,fm ,alpha_in,input_ang,guess)
        use num_vars, only: tol_NR
        use VMEC_vars, only: mpol, ntor, n_r, lam_c, lam_s, iotaf, iotah
        use fourier_ops, only: mesh_cs, f2r
        use utilities, only: zero_NR
        
        ! input / output
        real(dp) :: ang_B(n_r)                                                  ! theta(zeta)/zeta(theta)
        logical, intent(in) :: find_theta                                       ! whether theta or zeta is sought
        logical, intent(in) :: FM                                               ! whether full or half mesh quantity is sought
        real(dp), intent(in) :: alpha_in                                        ! alpha
        real(dp), intent(in) :: input_ang(n_r)                                  ! the input angle zeta/theta
        real(dp), intent(in), optional :: guess(n_r)                            ! optional input (guess) for theta/zeta
        
        ! local variables (also used in child functions)
        integer :: jd, kd
        real(dp) :: lam_c_loc(0:mpol-1,-ntor:ntor,1:n_r)                        ! local version of HM l_c (either FM or HM)
        real(dp) :: lam_s_loc(0:mpol-1,-ntor:ntor,1:n_r)                        ! local version of HM l_s (either FM or HM)
        real(dp) :: q(n_r)                                                      ! local version of 1/iota (either FM or HM)
        real(dp) :: lam                                                         ! lambda
        real(dp) :: dlam                                                        ! angular derivative of lambda
        real(dp) :: tempcoeff(n_r)                                              ! temporary holds a coefficient
        
        ! local variables (not to be used in child functions)
        real(dp) :: alpha_calc                                                  ! calculated alpha, to check with the given alpha
        real(dp) :: ang_NR                                                      ! temporary solution for a given r, iteration
        integer :: norm_start                                                   ! either 1 (FM) or 2 (HM)
        
        ! convert lam_c and lam_s to FM if FM is .true., select correct q
        if (FM) then
            do jd = -ntor, ntor
                do kd = 0, mpol-1
                    tempcoeff = lam_c(kd,jd,:)
                    lam_c_loc(kd,jd,:) = h2f(tempcoeff)
                    tempcoeff = lam_s(kd,jd,:)
                    lam_s_loc(kd,jd,:) = h2f(tempcoeff)
                end do
            end do
            q = iotaf
        else 
            lam_c_loc = lam_c
            lam_s_loc = lam_s
            q = iotah
        end if
        
        ! for all normal points
        if (fm) then
            norm_start = 1
        else
            norm_start = 2
        end if
        ang_B = 0.0_dp
        norm: do kd = 1, n_r
            ! if first guess for theta is given
            if (present(guess)) then
                ang_NR = guess(kd)
            else if (kd.eq.1) then
                ang_NR = pi                                                     ! take pi, because it is in the middle of 0...2pi
            else                                                                ! take solution for previous flux surface
                ang_NR = ang_B(kd-1)
            end if
            
            ! Newton-Rhapson loop for current normal point
            ang_B(kd) = zero_NR(fun_ang_B,dfun_ang_B,ang_NR)
            
            ! do a check  whether the result is indeed alpha,  making use of the
            ! last lam and dlam that have been calculated in the child functions
            if (find_theta) then                                                ! looking for theta
                alpha_calc = input_ang(kd) - (ang_B(kd) + lam)*q(kd)
            else                                                                ! looking for zeta
                alpha_calc = ang_B(kd) - (input_ang(kd) + lam)*q(kd)
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
        
    contains
        ! function that returns  f = alpha -  alpha_0. It uses kd  from the main
        ! loop in the  parent function as the normal position  where to evaluate
        ! the quantitites, and input_ang(kd) as the angle for which the magnetic
        ! match is sought, as well as q, alpha_in, lam_s_loc and lam_c_loc
        function fun_ang_B(ang)
            ! input / output
            real(dp) :: fun_ang_B
            real(dp), intent(in) :: ang
            
            ! local variables
            real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                               ! factors of cosine and sine for inverse fourier transform
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            if (find_theta) then                                                ! looking for theta
                cs = mesh_cs(mpol,ntor,ang,input_ang(kd))
            else                                                                ! looking for zeta
                cs = mesh_cs(mpol,ntor,input_ang(kd),ang)
            end if
            
            ! calculate lambda
            lam = f2r(lam_c_loc(:,:,kd),lam_s_loc(:,:,kd),cs,mpol,ntor,1)
            
            ! calculate the output function
            if (find_theta) then                                                ! looking for theta
                fun_ang_B = input_ang(kd) - (ang+lam)*q(kd) - alpha_in
            else                                                                ! looking for zeta
                fun_ang_B = ang - (input_ang(kd)+lam)*q(kd) - alpha_in
            end if
        end function fun_ang_B
        
        ! function that returns  df/d(ang) = d(alpha -  alpha_0)/d(ang). It uses
        ! kd from  the main loop in  the parent function as  the normal position
        ! where to evaluate the quantitites,  and input_ang(kd) as the angle for
        ! which the magnetic match is sought,  as well as q, alpha_in, lam_s_loc
        ! and lam_c_loc
        function dfun_ang_B(ang)
            ! input / output
            real(dp) :: dfun_ang_B
            real(dp), intent(in) :: ang
            
            ! local variables
            real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                               ! factors of cosine and sine for inverse fourier transform
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            if (find_theta) then                                                ! looking for theta
                cs = mesh_cs(mpol,ntor,ang,input_ang(kd))
            else                                                                ! looking for zeta
                cs = mesh_cs(mpol,ntor,input_ang(kd),ang)
            end if
            
            ! calculate angular derivatives of lambda
            if (find_theta) then                                                ! looking for theta
                dlam = f2r(lam_c_loc(:,:,kd),lam_s_loc(:,:,kd),cs,mpol,ntor,3)
            else                                                                ! looking for zeta
                dlam = f2r(lam_c_loc(:,:,kd),lam_s_loc(:,:,kd),cs,mpol,ntor,4)
            end if
            
            ! calculate the output function
            if (find_theta) then                                                ! looking for theta
                dfun_ang_B = -(1.0_dp + dlam)*q(kd)
            else                                                                ! looking for zeta
                dfun_ang_B = 1.0_dp - dlam*q(kd)
            end if
        end function dfun_ang_B
    end function ang_B
end module eq_vars
