module VMEC_vars
    use num_vars, only: &
        &dp, max_str_ln
    use str_ops, only: r2str, i2str
    use output_ops, only: lvl_ud, writo
    use read_wout_mod , only: read_wout_file, read_wout_deallocate, &           ! from LIBSTELL
        &iasym, version_, lfreeb, &                                             ! whether it is symmetric, version number, free boundary or not
        &n_r => ns, mpol, ntor, xn, xm, mnmax, nfp, &                           ! mpol, ntor = # modes
        &phi, phi_r => phipf, &                                                 ! toroidal flux (FM), norm. deriv. of toroidal flux (FM)
        &iotaf, iotah => iotas, &                                               ! iota = 1/q (FM) and (HM)
        &presf, presh => pres, gmns, gmnc, &                                    ! pressure (FM) and (HM), jacobian (HM)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &bmns, bmnc, &                                                          ! magnitude of B (HM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf                                        ! max and min values of R, Z
    implicit none
    private
    public read_VMEC, &
        &mnmax, rmnc, mpol, ntor, n_r, R_c, R_s, Z_c, Z_s, lam_c, lam_s, &
        &presf, presh, rmax_surf, rmin_surf, zmax_surf, iotaf, iotah, &
        &VMEC_name, phi, phi_r, B_V_sub_s_M, B_V_sub_c_M, B_V_c_H, B_V_s_H, &
        &jac_V_c, jac_V_s

    real(dp), allocatable :: R_c(:,:,:), R_s(:,:,:), Z_c(:,:,:), Z_s(:,:,:), &  ! Coeff. of R and Z in (co)sine series (FM)
        &lam_c(:,:,:), lam_s(:,:,:), &                                          ! Coeff. of lambda in (co)sine series (HM)
        &B_V_sub_c_M(:,:,:,:), B_V_sub_s_M(:,:,:,:), &                          ! Coeff. of B_i in (co)sine series (last index: r,theta,phi) (FM, HM, HM)
        &B_V_c_H(:,:,:), B_V_s_H(:,:,:), &                                      ! Coeff. of magnitude of B (HM)
        &jac_V_c(:,:,:), jac_V_s(:,:,:)                                         ! Jacobian in VMEC coordinates (HM)
    
    character(len=max_str_ln) :: VMEC_name                                      ! will hold name of the VMEC input file
    

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
        
        ! Allocate and repack the Fourier coefficients to translate them for use
        ! in this code
        allocate(Z_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(Z_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(lam_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(lam_s(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(B_V_sub_c_M(0:mpol-1,-ntor:ntor,1:n_r,3))
        allocate(B_V_sub_s_M(0:mpol-1,-ntor:ntor,1:n_r,3))
        allocate(B_V_c_H(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(B_V_s_H(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(jac_V_c(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(jac_V_s(0:mpol-1,-ntor:ntor,1:n_r))
        
        R_c = repack(rmnc,mnmax,n_r,mpol,ntor,xm,xn)
        R_s = repack(rmns,mnmax,n_r,mpol,ntor,xm,xn)
        Z_c = repack(zmnc,mnmax,n_r,mpol,ntor,xm,xn)                            ! xn includes the nfp factor from VMEC, which will not be looked at in this code
        Z_s = repack(zmns,mnmax,n_r,mpol,ntor,xm,xn)
        lam_c = repack(lmnc,mnmax,n_r,mpol,ntor,xm,xn)
        lam_s = repack(lmns,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_sub_c_M(:,:,:,1) = repack(bsubsmnc,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_sub_s_M(:,:,:,1) = repack(bsubsmns,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_sub_c_M(:,:,:,2) = repack(bsubumnc,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_sub_s_M(:,:,:,2) = repack(bsubumns,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_sub_c_M(:,:,:,3) = repack(bsubvmnc,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_sub_s_M(:,:,:,3) = repack(bsubvmns,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_c_H(:,:,:) = repack(bmnc,mnmax,n_r,mpol,ntor,xm,xn)
        B_V_s_H(:,:,:) = repack(bmns,mnmax,n_r,mpol,ntor,xm,xn)
        jac_V_c(:,:,:) = repack(gmnc,mnmax,n_r,mpol,ntor,xm,xn)
        jac_V_s(:,:,:) = repack(gmns,mnmax,n_r,mpol,ntor,xm,xn)
        
        call lvl_ud(-1)
    end subroutine

end module VMEC_vars
