module plasma_vars
    use num_vars, only: dp
    use var_ops, only: r2str
    use output_ops, only: lvl_ud, writo
    use file_ops, only: &
        &VMEC_name
    use read_wout_mod , only: read_wout_file, read_wout_deallocate, &
        &iasym, version_, lfreeb, &
        &n_r => ns, mpol, ntor, xn, xm, mnmax, nfp, &
        &phi, iotaf, &                                                          ! toroidal flux (FM), iota (FM)
        &presf, gmns, gmnc, &                                                   ! pressure (FM), jacobian (HM)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf                                        ! max and min values of R, Z
    implicit none
    private
    public read_VMEC, calc_metric, &
        &mnmax, rmnc, mpol, ntor, nfp, n_r, &
        &hrr, hzz, htt, htz, hrt, hrz, jac, grr, gzz, gtt, gtz, grt, grz

    real(dp), allocatable :: hrr(:,:,:), hzz(:,:,:), htt(:,:,:), &              ! upper metrical coefficients (r: normal, t: theta, z: zeta) & jacobian
        &htz(:,:,:), hrt(:,:,:), hrz(:,:,:), jac(:,:,:)
    real(dp), allocatable :: grr(:,:,:), gzz(:,:,:), gtt(:,:,:), &              ! lower metrical coefficients (r: normal, t: theta, z: zeta)
        &gtz(:,:,:), grt(:,:,:), grz(:,:,:)

contains
    ! Reads the VMEC equilibrium data
    subroutine read_VMEC()
        integer :: istat                                                        ! holds status of read_wout_file

        call writo('Reading data from VMEC output "' &
            &// trim(VMEC_name) // '"')
        call lvl_ud(1)
        call read_wout_file(VMEC_name, istat)
        if (istat.ne.0) then                                                    ! not able to read from the VMEC file
            call writo("ERROR: can't read the VMEC file")
            stop
        end if
        
        call writo('VMEC version is ' // trim(r2str(version_)))
        if (lfreeb) then
            call writo('Free boundary VMEC')
        else
            call writo('Fixed boundary VMEC')
        end if

        call writo('Running some tests')
        call lvl_ud(1)
        if (mnmax.ne.((ntor+1)+(2*ntor+1)*(mpol-1))) then                       ! are ntor, mpol and mnmax what is expected?
            call writo('ERROR: inconsistency in ntor, mpol and mnmax')
            stop
        end if
        call lvl_ud(-1)

        call writo('Reading the plasma variables')
        call lvl_ud(-1)
        call writo('Data from VMEC output successfully read')

        call writo('Reading the grid parameters')
        call lvl_ud(1)
        call calc_metric
        call lvl_ud(-1)
    end subroutine

    ! Calculate the metric matrix from the VMEC file
    subroutine calc_metric()
        use fourier_ops, only: repack, mesh_cs, f2r
        use output_ops, only: writo, print_ar_2, print_ar_1
        use num_vars, only: min_theta, max_theta, min_zeta, max_zeta, n_theta, &
            &n_zeta
            
        real(dp), allocatable :: R_c(:,:,:), R_s(:,:,:), Z_c(:,:,:), &          ! Coeff. of R, Z, lambda in (co)sine series
            &Z_s(:,:,:), l_c(:,:,:), l_s(:,:,:)
        real(dp), allocatable :: R(:,:,:,:)                                     ! R and derivatives in real (as opposed to Fourier) space (see below)
        real(dp), allocatable :: Z(:,:,:,:)                                     ! Z and derivatives in real (as opposed to Fourier) space (see below)
        real(dp), allocatable :: l(:,:,:,:)                                     ! lambda and derivatives in real (as opposed to Fourier) space (see below)
        real(dp) :: delta_r, common_factor                                      ! normal step size, common factor for lower metric coefficients
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        real(dp) :: theta(n_theta), zeta(n_zeta)                                ! array of theta's and zeta's where metrics are to be calculated
        integer :: id, jd, kd
        real(dp), allocatable :: jac_low(:,:,:), jac_diff                       ! jacobian calculated using det. of lower metrical factors and max diff
        real(dp), allocatable :: max_R_an(:)                                    ! max of R, calculated analytically
        real(dp) :: testarray(2,3)
        
        testarray = 1_dp
        testarray(1,2) = 2_dp
        testarray(2,3) = 3_dp
        call writo('Calculating the matrix of metric elements from VMEC &
            &coordinates to flux coordinates')
        ! ¿¿¿ MAYBE TAKE INTO ACCOUNT THE ADDITIONAL POSSIBLE ASSYMETRY FOR FASTER CODE???

        ! Allocate and repack the Fourier coefficients of the variables R, Z and
        ! lambda to translate them for use in this code
        allocate(Z_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(Z_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(l_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(l_s(0:mpol-1,-ntor:ntor,1:n_r))
        Z_c = repack(zmnc,mnmax,n_r,mpol,ntor,xm,xn/nfp)
        Z_s = repack(zmns,mnmax,n_r,mpol,ntor,xm,xn/nfp)
        R_c = repack(rmnc,mnmax,n_r,mpol,ntor,xm,xn/nfp)
        R_s = repack(rmns,mnmax,n_r,mpol,ntor,xm,xn/nfp)
        l_c = repack(lmnc,mnmax,n_r,mpol,ntor,xm,xn/nfp)
        l_s = repack(lmns,mnmax,n_r,mpol,ntor,xm,xn/nfp)
        
        ! calc array of angular mesh points (for now, VMEC's mesh is used)
        theta = 0.0_dp
        zeta = 0.0_dp
        theta = ang_mesh(n_theta,min_theta,max_theta)
        zeta = ang_mesh(n_zeta,min_zeta,max_zeta)
        
        ! Allocate the variables R, Z and lambda in real space
        ! index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative
        allocate(R(n_theta,n_zeta,n_r,4))
        R = 0.0_dp
        allocate(Z(n_theta,n_zeta,n_r,4))
        Z = 0.0_dp
        allocate(l(n_theta,n_zeta,n_r,4))
        l = 0.0_dp
        
        ! for every mesh point, calculate the metrics
        do id = 1,n_theta
            do jd = 1,n_zeta
                ! calculate the (co)sines at the current mesh points
                cs = mesh_cs(mpol,ntor,nfp,theta(id),zeta(jd))
                ! calculate the variables R, Z, and their angular derivatives for all normal points
                do kd = 1, n_r
                    R(id,jd,kd,:) = &
                        &f2r(R_c(:,:,kd),R_s(:,:,kd),cs,mpol,ntor,nfp)
                    Z(id,jd,kd,:) = &
                        &f2r(Z_c(:,:,kd),Z_s(:,:,kd),cs,mpol,ntor,nfp)
                end do
                ! numerically calculate normal derivatives at the currrent angular points
                ! first normal point
                delta_r = phi(2)-phi(1)                                         ! step size between first points
                R(id,jd,1,2) = (R(id,jd,2,1)-R(id,jd,1,1))/delta_r              ! forward difference
                Z(id,jd,1,2) = (Z(id,jd,2,1)-Z(id,jd,1,1))/delta_r              ! forward difference
                do kd = 3, n_r
                    delta_r = phi(kd)-phi(kd-2)
                    R(id,jd,kd-1,2) = (R(id,jd,kd,1)-R(id,jd,kd-2,1))/delta_r
                    Z(id,jd,kd-1,2) = (Z(id,jd,kd,1)-Z(id,jd,kd-2,1))/delta_r
                end do
                ! last normal point
                delta_r = phi(n_r)-phi(n_r-1)
                R(id,jd,n_r,2) = (R(id,jd,n_r,1)-R(id,jd,n_r-1,1))/delta_r
                Z(id,jd,n_r,2) = (Z(id,jd,n_r,1)-Z(id,jd,n_r-1,1))/delta_r
            end do
        end do
        ! analytically calculate maximum R by summing the modes (theta = 0)
        allocate(max_R_an(n_r))
        max_R_an = 0.0_dp
        do kd = 1,n_r
            max_R_an(kd) = sum(R_c(:,:,kd))
        end do
        
        ! output a message if the found R  and Z values are not within the VMEC-
        ! provided bounds
        call writo('Testing whether R is in agreement with VMEC')
        call lvl_ud(1)
        call within_bounds(R(:,:,n_r,1),rmin_surf,rmax_surf)
        call lvl_ud(-1)
        call writo('Testing whether Z is in agreement with VMEC')
        call lvl_ud(1)
        call within_bounds(Z(:,:,n_r,1),-zmax_surf,zmax_surf)
        call lvl_ud(-1)
        
        ! allocate the upper metric coefficients
        allocate(hrr(n_theta,n_zeta,n_r))
        allocate(hzz(n_theta,n_zeta,n_r))
        allocate(htt(n_theta,n_zeta,n_r))
        allocate(htz(n_theta,n_zeta,n_r))
        allocate(hrt(n_theta,n_zeta,n_r))
        allocate(hrz(n_theta,n_zeta,n_r))
        allocate(jac(n_theta,n_zeta,n_r))
        hrr = 0.0_dp
        hzz = 0.0_dp
        htt = 0.0_dp
        htz = 0.0_dp
        hrt = 0.0_dp
        hrz = 0.0_dp
        jac = 0.0_dp
        
        ! calculate upper metric coefficients and Jacobian using formula's in [1]
        do id = 1,n_theta
            do jd = 1,n_zeta
                do kd = 1, n_r
                    jac(id,jd,kd) = R(id,jd,kd,1)*(R(id,jd,kd,3)*Z(id,jd,kd,2)-&
                        &R(id,jd,kd,2)*Z(id,jd,kd,3))
                    hzz(id,jd,kd) = 1.0_dp/(R(id,jd,kd,1)**2)
                    htz(id,jd,kd) = (R(id,jd,kd,2)*Z(id,jd,kd,4)-R(id,jd,kd,4)&
                        &*Z(id,jd,kd,2))/(R(id,jd,kd,1)*jac(id,jd,kd))
                    hrr(id,jd,kd) = (R(id,jd,kd,1)/jac(id,jd,kd))**2 * &
                        &(Z(id,jd,kd,3)**2 + R(id,jd,kd,3)**2)
                    htt(id,jd,kd) = (R(id,jd,kd,1)/jac(id,jd,kd))**2 * &
                        &(Z(id,jd,kd,2)**2 + R(id,jd,kd,2)**2) + &
                        &(R(id,jd,kd,1)*htz(id,jd,kd))**2
                    hrt(id,jd,kd) = -(R(id,jd,kd,1)/jac(id,jd,kd))**2* &
                        &(Z(id,jd,kd,3)*Z(id,jd,kd,2)+&
                        &R(id,jd,kd,3)*R(id,jd,kd,2))
                    hrz(id,jd,kd) = 0.0_dp
                end do
            end do
        end do
        
        ! allocate the lower metric coefficients
        allocate(grr(n_theta,n_zeta,n_r))
        allocate(gzz(n_theta,n_zeta,n_r))
        allocate(gtt(n_theta,n_zeta,n_r))
        allocate(gtz(n_theta,n_zeta,n_r))
        allocate(grt(n_theta,n_zeta,n_r))
        allocate(grz(n_theta,n_zeta,n_r))
        allocate(jac_low(n_theta,n_zeta,n_r))
        grr = 0.0_dp
        gzz = 0.0_dp
        gtt = 0.0_dp
        gtz = 0.0_dp
        grt = 0.0_dp
        grz = 0.0_dp
        jac_low = 0.0_dp
        ! calculate lower metric coefficients and inverse Jacobian using formula's in [1]
        do id = 1,n_theta
            do jd = 1,n_zeta
                do kd = 1, n_r
                    grr(id,jd,kd) = R(id,jd,kd,2)**2 + Z(id,jd,kd,2)**2
                    gtt(id,jd,kd) = R(id,jd,kd,3)**2 + Z(id,jd,kd,3)**2
                    grt(id,jd,kd) = R(id,jd,kd,2)*R(id,jd,kd,3) &
                        &+ Z(id,jd,kd,2)*Z(id,jd,kd,3)
                    common_factor = - R(id,jd,kd,1)**2 * htz(id,jd,kd)
                    gzz(id,jd,kd) = R(id,jd,kd,1)**2 + gtt(id,jd,kd) * &
                        &common_factor**2
                    gtz(id,jd,kd) = gtt(id,jd,kd)*common_factor
                    grz(id,jd,kd) = grt(id,jd,kd)*common_factor
                    jac_low(id,jd,kd) = R(id,jd,kd,1)**2 * (grr(id,jd,kd)* &
                        &gtt(id,jd,kd) - grt(id,jd,kd)**2)
                end do
            end do
        end do
        
        ! test whether jac_low == jac
        call writo('Test whether Jacobians calculated from upper and lower &
            &metrics are in agreement')
        call lvl_ud(1)
        jac_diff = maxval(abs(jac(:,:,5)**2-jac_low(:,:,5)))
        if (jac_diff.gt.maxval(jac(:,:,:))/100) then                            ! difference is too big -> error
            call writo('ERROR: maximum difference between Jacobian &
                &calculated from lower and upper metrics is ' // &
                &trim(r2str(jac_diff/maxval(jac(:,:,:))*100)) // '%')
        else if (jac_diff.gt.1E5*epsilon(1.0_dp)) then                          ! difference is substantial -> warning
            call writo('WARNING: maximum difference between Jacobian &
                &calculated from lower and upper metrics is not completely&
                & negligible')
        end if
        call lvl_ud(-1)

    contains
        ! calculate the angular points in the mesh
        ! currently trivial implementation, but can be extended for mesh
        ! adaptation
        function ang_mesh(n_ang,min_ang, max_ang)
            real(dp), allocatable :: ang_mesh(:)
            real(dp), intent(in) :: min_ang, max_ang
            integer, intent(in) :: n_ang
            
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
            
            allocate(ang_mesh(n_ang))
            ! There are (n_ang-1) pieces in the total interval but the last one 
            ! is not needed as the functions are all periodic
            delta_ang = (max_ang-min_ang)/(n_ang)
            ang_mesh = 0.0_dp
            
            ang_mesh(1) = min_ang
            do id = 2,n_ang
                ang_mesh(id) = ang_mesh(id-1) + delta_ang
            end do
        end function ang_mesh
            
        ! display whether the calculated table for R or Z fits within the limits
        ! outputted by VMEC
        subroutine within_bounds(var,min_VMEC,max_VMEC)
            use var_ops, only: r2strt
            
            real(dp) :: var(:,:)
            real(dp) :: min_VMEC, max_VMEC
            
            real(dp) :: min_frac, max_frac, margin
            
            margin = 1E-2                                                       ! 0.1% for comparison 
            min_frac = abs(2*(minval(var)-min_VMEC)/(minval(var)+min_VMEC))
            max_frac = abs(2*(maxval(var)-max_VMEC)/(maxval(var)+max_VMEC))
            
            if (min_frac.lt.margin .and. max_frac.lt.margin) then
                return
            else if (min_frac.ge.margin) then                                   ! too low minimum
                call writo('WARNING: minimum of variable in &
                &real angular space deviates '// trim(r2strt(100*min_frac)) &
                &//'% from VMEC provided minimum... -> Maybe use an even number&
                & of poloidal points, or more angular mesh points in general?')
            else if (max_frac.ge.margin) then                                   ! too high maximum
                call writo('WARNING: maximum of variable in &
                &real angular space deviates '// trim(r2strt(100*max_frac)) &
                &//'% from VMEC provided maximum... &
                &-> Maybe use more angular mesh points?')
            end if
        end subroutine
    end subroutine
end module plasma_vars
