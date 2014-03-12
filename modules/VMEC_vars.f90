module VMEC_vars
    use num_vars, only: dp
    use var_ops, only: r2str, i2str
    use output_ops, only: lvl_ud, writo
    use file_ops, only: &
        &VMEC_name
    use read_wout_mod , only: read_wout_file, read_wout_deallocate, &
        &iasym, version_, lfreeb, &
        &n_r => ns, mpol, ntor, xn, xm, mnmax, nfp, &                           ! mpol, ntor = # modes
        &phi, iotaf, &                                                          ! toroidal flux (FM), iota (FM)
        &presf, gmns, gmnc, &                                                   ! pressure (FM), jacobian (HM)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf                                        ! max and min values of R, Z
    implicit none
    private
    public read_VMEC, &
        &mnmax, rmnc, mpol, ntor, n_r, R_c, R_s, Z_c, Z_s, l_c, l_s, &
        &rmax_surf, rmin_surf, zmax_surf, iotaf

    real(dp), allocatable :: R_c(:,:,:), R_s(:,:,:), Z_c(:,:,:), &              ! Coeff. of R, Z, lambda in (co)sine series
        &Z_s(:,:,:), l_c(:,:,:), l_s(:,:,:)

contains
    ! Reads the VMEC equilibrium data
    subroutine read_VMEC()
        use fourier_ops, only: repack
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

        call writo('Updating some variables')
        ! make ntor already include the factor nfp for easier handling in this code
        ntor = ntor*nfp

        call lvl_ud(-1)
        call writo('Data from VMEC output successfully read')

        call writo('Reading the grid parameters')
        call lvl_ud(1)

        ! Allocate and repack the Fourier coefficients of the variables R, Z and
        ! lambda to translate them for use in this code
        allocate(Z_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(Z_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(l_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(l_s(0:mpol-1,-ntor:ntor,1:n_r))
        R_c = repack(rmnc,mnmax,n_r,mpol,ntor,xm,xn)
        R_s = repack(rmns,mnmax,n_r,mpol,ntor,xm,xn)
        Z_c = repack(zmnc,mnmax,n_r,mpol,ntor,xm,xn)                            ! xn includes the nfp factor from VMEC, which will not be looked at in this code
        Z_s = repack(zmns,mnmax,n_r,mpol,ntor,xm,xn)
        l_c = repack(lmnc,mnmax,n_r,mpol,ntor,xm,xn)
        l_s = repack(lmns,mnmax,n_r,mpol,ntor,xm,xn)

        call lvl_ud(-1)
    end subroutine

end module VMEC_vars
