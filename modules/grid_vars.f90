module grid_vars
    use num_vars, only: dp
    use output_ops, only: writo, lvl_ud, print_ar_2

    implicit none
    private
    public calc_ang_mesh, calc_RZl, &
        &theta, zeta, rad, n_theta, n_zeta, R, Z, lam

    real(dp), allocatable :: theta(:), zeta(:), rad(:)                          ! grid points
    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative)
    real(dp), allocatable :: R(:,:,:,:), Z(:,:,:,:), lam(:,:,:,:)
    integer :: n_theta, n_zeta     

contains
    ! calculate the angular points in the mesh
    ! currently trivial implementation, but can be extended for mesh adaptation
    ! or, for example, following a field line
    function calc_ang_mesh(n_ang,min_ang, max_ang)
        ! input and output
        real(dp), allocatable :: calc_ang_mesh(:)
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
        
        allocate(calc_ang_mesh(n_ang))
        ! There are (n_ang-1) pieces in the total interval but the last one 
        ! is not needed as the functions are all periodic
        delta_ang = (max_ang-min_ang)/(n_ang)
        calc_ang_mesh = 0.0_dp
        
        calc_ang_mesh(1) = min_ang
        do id = 2,n_ang
            calc_ang_mesh(id) = calc_ang_mesh(id-1) + delta_ang
        end do
    end function calc_ang_mesh

    ! calculate the coordinates R  and Z and l in real  space from their Fourier
    ! decomposition using the grid points currently stored in the variables 
    ! theta, zeta, rad
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
        allocate(R(n_theta,n_zeta,n_r,4)); R = 0.0_dp
        allocate(Z(n_theta,n_zeta,n_r,4)); Z = 0.0_dp
        allocate(lam(n_theta,n_zeta,n_r,4)); lam = 0.0_dp
        
        ! do calculations for all angular points
        tor: do jd = 1, n_zeta
            pol: do id = 1,n_theta
                ! calculate the (co)sines at the current mesh points
                cs = mesh_cs(mpol,ntor,theta(id),zeta(jd))
                
                ! calculate the variables R, Z, and their angular derivatives for all normal points and current angular point
                do kd = 1, n_r
                    R(id,jd,kd,:) = &
                        &f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor)
                    Z(id,jd,kd,:) = &
                        &f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor)
                    lam(id,jd,kd,:) = &
                        &f2r(l_c(:,:,kd),l_s(:,:,kd),cs,mpol,ntor)
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
        write(*,*) 'maximum value of lambda', maxval(lam(:,:,:,1))
        write(*,*) 'minimum value of lambda', minval(lam(:,:,:,1))
        
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
            use var_ops, only: r2strt
            
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
end module grid_vars
