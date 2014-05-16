!------------------------------------------------------------------------------!
!   functions and routines that concern variables given by VMEC                !
!------------------------------------------------------------------------------!
module VMEC_vars
    use num_vars, only: &
        &dp, max_str_ln, pi
    use str_ops, only: r2str, i2str
    use output_ops, only: lvl_ud, writo, print_ar_1, print_ar_2, write_out
    use read_wout_mod , only: read_wout_file, read_wout_deallocate, &           ! from LIBSTELL
        &iasym, version_, lfreeb, &                                             ! whether it is symmetric, version number, free boundary or not
        &n_r => ns, mpol, ntor, xn, xm, mnmax, nfp, &                           ! mpol, ntor = # modes
        &phi, phi_r => phipf, &                                                 ! toroidal flux (FM), norm. deriv. of toroidal flux (FM)
        &phi_r_H => phip, &                                                     ! norm. deriv. of toroidal flux (HM) / 2 pi
        &iotaf, iotah => iotas, &                                               ! iota = 1/q (FM) and (HM)
        &presf, gmns, gmnc, &                                                   ! pressure (FM) jacobian (HM)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &bmns, bmnc, &                                                          ! magnitude of B (HM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf                                        ! max and min values of R, Z
    implicit none
    private
    public read_VMEC, &
        &mnmax, rmnc, mpol, ntor, n_r, R_c, R_s, Z_c, Z_s, L_c, L_s, &
        &presf, rmax_surf, rmin_surf, zmax_surf, iotaf, iotah, &
        &VMEC_name, phi, phi_r, &
        &phi_r_H, B_V_sub_s_M, B_V_sub_c_M, B_V_c_H, B_V_s_H, &
        &jac_V_H_c, jac_V_H_s

    real(dp), allocatable :: &
        &R_c(:,:,:,:), R_s(:,:,:,:), &                                          ! Coeff. of R in (co)sine series (FM) and norm. deriv.
        &Z_c(:,:,:,:), Z_s(:,:,:,:), &                                          ! Coeff. of Z in (co)sine series (FM) and norm. deriv.
        &L_c(:,:,:,:), L_s(:,:,:,:)                                             ! Coeff. of lambda in (co)sine series (HM) and norm. deriv.

    real(dp), allocatable :: &
        &B_V_sub_c_M(:,:,:,:), B_V_sub_s_M(:,:,:,:), &                          ! Coeff. of B_i in (co)sine series (last index: r,theta,phi) (FM, HM, HM)
        &B_V_c_H(:,:,:), B_V_s_H(:,:,:), &                                      ! Coeff. of magnitude of B (HM)
        &jac_V_H_c(:,:,:), jac_V_H_s(:,:,:)                                     ! Jacobian in VMEC coordinates (HM)
    
    character(len=max_str_ln) :: VMEC_name                                      ! will hold name of the VMEC input file
    

contains
    ! Reads the VMEC equilibrium data
    subroutine read_VMEC()
        use fourier_ops, only: repack
        use utilities, only: VMEC_norm_deriv, VMEC_conv_FHM
        use num_vars, only: ltest, max_deriv
        
        ! local variables
        integer :: istat                                                        ! holds status of read_wout_file
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: L_c_H(:,:,:,:)                                 ! temporary HM variable
        real(dp), allocatable :: L_s_H(:,:,:,:)                                 ! temporary HM variable
        
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
        
        call writo('Running tests')
        call lvl_ud(1)
        if (mnmax.ne.((ntor+1)+(2*ntor+1)*(mpol-1))) then                       ! are ntor, mpol and mnmax what is expected?
            call writo('ERROR: inconsistency in ntor, mpol and mnmax')
            stop
        end if
        call lvl_ud(-1)
        
        call writo('Updating variables')
        ! make ntor already include the factor nfp for easier handling in this code
        ntor = ntor*nfp
        
        call lvl_ud(-1)
        call writo('Data from VMEC output successfully read')
        
        call writo('Reading the grid parameters')
        call lvl_ud(1)
        
        ! rescale phi_r_H
        phi_r_H = -phi_r_H*2*pi
        
        ! Allocate and repack the Fourier coefficients to translate them for use
        ! in this code
        allocate(R_c(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(R_s(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(Z_c(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(Z_s(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(L_c(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(L_s(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(L_c_H(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        allocate(L_s_H(0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3)))
        
        ! factors R_c,s; Z_c,s and L_C,s and HM varieties
        R_c(:,:,:,0) = repack(rmnc,mnmax,n_r,mpol,ntor,xm,xn)
        R_s(:,:,:,0) = repack(rmns,mnmax,n_r,mpol,ntor,xm,xn)
        Z_c(:,:,:,0) = repack(zmnc,mnmax,n_r,mpol,ntor,xm,xn)
        Z_s(:,:,:,0) = repack(zmns,mnmax,n_r,mpol,ntor,xm,xn)
        L_c_H(:,:,:,0) = repack(lmnc,mnmax,n_r,mpol,ntor,xm,xn)
        L_s_H(:,:,:,0) = repack(lmns,mnmax,n_r,mpol,ntor,xm,xn)
        
        ! normal derivatives of these factors
        do kd = 1,max_deriv(1)
            do jd = -ntor,ntor
                do id = 0,mpol-1
                    call VMEC_norm_deriv(R_c(id,jd,:,0),R_c(id,jd,:,kd),&
                        &n_r-1._dp,kd,1)
                    call VMEC_norm_deriv(R_s(id,jd,:,0),R_s(id,jd,:,kd),&
                        &n_r-1._dp,kd,1)
                    call VMEC_norm_deriv(Z_c(id,jd,:,0),Z_c(id,jd,:,kd),&
                        &n_r-1._dp,kd,1)
                    call VMEC_norm_deriv(Z_s(id,jd,:,0),Z_s(id,jd,:,kd),&
                        &n_r-1._dp,kd,1)
                    call VMEC_norm_deriv(L_c_H(id,jd,:,0),L_c_H(id,jd,:,kd),&
                        &n_r-1._dp,kd,1)
                    call VMEC_norm_deriv(L_s_H(id,jd,:,0),L_s_H(id,jd,:,kd),&
                        &n_r-1._dp,kd,1)
                end do
            end do
        end do
        
        ! conversion HM -> FM (L)
        do kd = 0,max_deriv(3)
            do jd = -ntor,ntor
                do id = 0,mpol-1
                    call VMEC_conv_FHM(L_c_H(id,jd,:,kd),L_c(id,jd,:,kd),&
                        &.false.)
                    call VMEC_conv_FHM(L_s_H(id,jd,:,kd),L_s(id,jd,:,kd),&
                        &.false.)
                end do
            end do
        end do
        
        ! for tests
        if (ltest) then
            allocate(B_V_sub_c_M(0:mpol-1,-ntor:ntor,1:n_r,3))
            allocate(B_V_sub_s_M(0:mpol-1,-ntor:ntor,1:n_r,3))
            allocate(B_V_c_H(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(B_V_s_H(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(jac_V_H_c(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(jac_V_H_s(0:mpol-1,-ntor:ntor,1:n_r))
            
            B_V_sub_c_M(:,:,:,1) = repack(bsubsmnc,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_sub_s_M(:,:,:,1) = repack(bsubsmns,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_sub_c_M(:,:,:,2) = repack(bsubumnc,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_sub_s_M(:,:,:,2) = repack(bsubumns,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_sub_c_M(:,:,:,3) = repack(bsubvmnc,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_sub_s_M(:,:,:,3) = repack(bsubvmns,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_c_H(:,:,:) = repack(bmnc,mnmax,n_r,mpol,ntor,xm,xn)
            B_V_s_H(:,:,:) = repack(bmns,mnmax,n_r,mpol,ntor,xm,xn)
            jac_V_H_c(:,:,:) = -repack(gmnc,mnmax,n_r,mpol,ntor,xm,xn)
            jac_V_H_s(:,:,:) = -repack(gmns,mnmax,n_r,mpol,ntor,xm,xn)
        end if
        
        call lvl_ud(-1)
        call writo('Grid parameters successfully read')
    end subroutine

end module VMEC_vars
