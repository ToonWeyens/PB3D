!------------------------------------------------------------------------------!
!   functions and routines that concern variables given by VMEC                !
!------------------------------------------------------------------------------!
module VMEC_vars
#include <PB3D_macros.h>
    use num_vars, only: &
        &dp, max_str_ln, pi
    use str_ops, only: r2str, i2str
    use message_ops, only: lvl_ud, writo, print_ar_1, print_ar_2
    use read_wout_mod, only: read_wout_file, read_wout_deallocate, &            ! from LIBSTELL
        &lasym, VMEC_version => version_, lfreeb, &                             ! stellerator symmetry, version number, free boundary or not
        &n_r_eq => ns, mpol, ntor, xn, xm, mnmax, nfp, &                        ! mpol, ntor = # modes
        &phi, phi_r => phipf, &                                                 ! toroidal flux (FM), norm. deriv. of toroidal flux (FM)
        &iotaf, &                                                               ! iota = 1/q (FM)
        &presf, gmns, gmnc, &                                                   ! pressure (FM) jacobian (HM)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &bmns, bmnc, &                                                          ! magnitude of B (HM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf, &                                     ! max and min values of R, Z
        &VMEC_use_pol_flux => lrfp, &                                           ! whether or not the poloidal flux is used as radial variable
        &gam => gamma                                                           ! gamma in adiabatic law
    implicit none
    private
    public read_VMEC, dealloc_VMEC, repack, &
        &mnmax, rmnc, mpol, ntor, nfp, n_r_eq, R_c, R_s, Z_c, Z_s, L_c, L_s, &
        &presf, rmax_surf, rmin_surf, zmax_surf, iotaf, lasym, &
        &VMEC_use_pol_flux, lfreeb, phi, phi_r, VMEC_version, gam, &
        &B_V_sub_s_M, B_V_sub_c_M, B_V_c_H, B_V_s_H, &
        &jac_V_c_H, jac_V_s_H

    real(dp), allocatable :: &
        &R_c(:,:,:,:), R_s(:,:,:,:), &                                          ! Coeff. of R in (co)sine series (FM) and norm. deriv.
        &Z_c(:,:,:,:), Z_s(:,:,:,:), &                                          ! Coeff. of Z in (co)sine series (FM) and norm. deriv.
        &L_c(:,:,:,:), L_s(:,:,:,:)                                             ! Coeff. of lambda in (co)sine series (HM) and norm. deriv.

    real(dp), allocatable :: &
        &B_V_sub_c_M(:,:,:,:), B_V_sub_s_M(:,:,:,:), &                          ! Coeff. of B_i in (co)sine series (last index: r,theta,phi) (FM, HM, HM)
        &B_V_c_H(:,:,:), B_V_s_H(:,:,:), &                                      ! Coeff. of magnitude of B (HM)
        &jac_V_c_H(:,:,:), jac_V_s_H(:,:,:)                                     ! Jacobian in VMEC coordinates (HM)

contains
    ! Reads the VMEC equilibrium data
    ! [MPI] only global master
    integer function read_VMEC() result(ierr)
        use utilities, only: VMEC_norm_deriv, VMEC_conv_FHM
        use num_vars, only: max_deriv, eq_i, glb_rank, eq_name
#if ldebug
        use num_vars, only: ltest
#endif
        character(*), parameter :: rout_name = 'read_VMEC'
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: L_c_H(:,:,:,:)                                 ! temporary HM variable
        real(dp), allocatable :: L_s_H(:,:,:,:)                                 ! temporary HM variable
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then                                                 ! only global master
            call writo('Reading data from VMEC output "' &
                &// trim(eq_name) // '"')
            call lvl_ud(1)
            
            ! read VMEC output using LIBSTELL
            !call read_wout_file(eq_name, ierr)                                ! DOESN'T WORK WITH VMEC .TXT OUTPUT -> use file number
            call read_wout_file(eq_i, ierr)                                     ! read the VMEC file number
            CHCKERR('Can''t read the VMEC file')
            close(eq_i)                                                         ! close the VMEC file
            
            call writo('VMEC version is ' // trim(r2str(VMEC_version)))
            if (lfreeb) then
                call writo('Free boundary VMEC')
                err_msg = 'Free boundary VMEC is not yet supported by PB3D...'
                ierr = 1
                CHCKERR(err_msg)
            else
                call writo('Fixed boundary VMEC')
            end if
            if (lasym) then
                call writo('No stellerator symmetry')
            else
                call writo('Stellerator symmetry applicable')
                call writo('¡¡¡ ITS USAGE COULD BE IMPLEMENTED !!!')
            end if
            call writo('VMEC has '//trim(i2str(mpol))//' poloidal and '&
                &//trim(i2str(ntor))//' toroidal modes, defined on '&
                &//trim(i2str(n_r_eq))//' flux surfaces')
            
            call writo('Running tests')
            call lvl_ud(1)
            if (mnmax.ne.((ntor+1)+(2*ntor+1)*(mpol-1))) then                   ! are ntor, mpol and mnmax what is expected?
                err_msg = 'Inconsistency in ntor, mpol and mnmax'
                ierr = 1
                CHCKERR(err_msg)
            end if
            call lvl_ud(-1)
            
            call writo('Updating variables')
            xn = xn/nfp                                                         ! so we have xn excluding nfp
            
            call lvl_ud(-1)
            call writo('Data from VMEC output successfully read')
            
            call writo('Reading the grid parameters')
            call lvl_ud(1)
            
            ! Allocate and repack the Fourier coefficients to translate them for
            ! use in this code
            allocate(R_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
            allocate(Z_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
            allocate(L_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
            allocate(L_s_H(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
            !if (lasym) then                                                     ! following only needed in assymetric situations
                allocate(R_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
                allocate(Z_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
                allocate(L_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
                allocate(L_c_H(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv(3)))
            !end if
            
            ! factors R_c,s; Z_c,s and L_C,s and HM varieties
            R_c(:,:,:,0) = repack(rmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            Z_s(:,:,:,0) = repack(zmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            L_s_H(:,:,:,0) = repack(lmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            !if (lasym) then                                                     ! following only needed in assymetric situations
                R_s(:,:,:,0) = repack(rmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
                Z_c(:,:,:,0) = repack(zmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
                L_c_H(:,:,:,0) = repack(lmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            !end if
            
            ! normal derivatives of these factors
            do kd = 1,max_deriv(1)
                do jd = -ntor,ntor
                    do id = 0,mpol-1
                        ierr = VMEC_norm_deriv(R_c(id,jd,:,0),R_c(id,jd,:,kd),&
                            &n_r_eq-1._dp,kd,1)
                        CHCKERR('')
                        ierr = VMEC_norm_deriv(Z_s(id,jd,:,0),Z_s(id,jd,:,kd),&
                            &n_r_eq-1._dp,kd,1)
                        CHCKERR('')
                        ierr = VMEC_norm_deriv(L_s_H(id,jd,:,0),&
                            &L_s_H(id,jd,:,kd),n_r_eq-1._dp,kd,1)
                        CHCKERR('')
                        !if (lasym) then                                         ! following only needed in assymetric situations
                            ierr = VMEC_norm_deriv(R_s(id,jd,:,0),&
                                &R_s(id,jd,:,kd),n_r_eq-1._dp,kd,1)
                            CHCKERR('')
                            ierr = VMEC_norm_deriv(Z_c(id,jd,:,0),&
                                &Z_c(id,jd,:,kd),n_r_eq-1._dp,kd,1)
                            CHCKERR('')
                            ierr = VMEC_norm_deriv(L_c_H(id,jd,:,0),&
                                &L_c_H(id,jd,:,kd),n_r_eq-1._dp,kd,1)
                            CHCKERR('')
                        !end if
                    end do
                end do
            end do
            
            ! conversion HM -> FM (L)
            do kd = 0,max_deriv(1)
                do jd = -ntor,ntor
                    do id = 0,mpol-1
                        ierr = VMEC_conv_FHM(L_s_H(id,jd,:,kd),L_s(id,jd,:,kd),&
                            &.false.)
                        CHCKERR('')
                        !if (lasym) then                                         ! following only needed in assymetric situations
                            ierr = VMEC_conv_FHM(L_c_H(id,jd,:,kd),&
                                &L_c(id,jd,:,kd),.false.)
                            CHCKERR('')
                        !end if
                    end do
                end do
            end do
            
#if ldebug
            ! for tests
            if (ltest) then
                allocate(B_V_sub_c_M(0:mpol-1,-ntor:ntor,1:n_r_eq,3))
                allocate(B_V_sub_s_M(0:mpol-1,-ntor:ntor,1:n_r_eq,3))
                allocate(B_V_c_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
                allocate(B_V_s_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
                allocate(jac_V_c_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
                allocate(jac_V_s_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
                
                B_V_sub_c_M(:,:,:,1) = repack(bsubsmnc,mnmax,n_r_eq,mpol,ntor,&
                    &xm,xn)
                B_V_sub_s_M(:,:,:,1) = repack(bsubsmns,mnmax,n_r_eq,mpol,ntor,&
                    &xm,xn)
                B_V_sub_c_M(:,:,:,2) = repack(bsubumnc,mnmax,n_r_eq,mpol,ntor,&
                    &xm,xn)
                B_V_sub_s_M(:,:,:,2) = repack(bsubumns,mnmax,n_r_eq,mpol,ntor,&
                    &xm,xn)
                B_V_sub_c_M(:,:,:,3) = repack(bsubvmnc,mnmax,n_r_eq,mpol,ntor,&
                    &xm,xn)
                B_V_sub_s_M(:,:,:,3) = repack(bsubvmns,mnmax,n_r_eq,mpol,ntor,&
                    &xm,xn)
                B_V_c_H(:,:,:) = repack(bmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
                B_V_s_H(:,:,:) = repack(bmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
                jac_V_c_H(:,:,:) = -repack(gmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
                jac_V_s_H(:,:,:) = -repack(gmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            end if
#endif
            
            call lvl_ud(-1)
            call writo('Grid parameters successfully read')
        end if
    end function read_VMEC
    
    ! deallocates VMEC quantities that are not used anymore
    subroutine dealloc_VMEC
        deallocate(phi,phi_r)
        deallocate(iotaf)
        deallocate(presf)
        if (allocated(gmns)) deallocate(gmns)
        if (allocated(gmnc)) deallocate(gmnc)
        if (allocated(bsubumns)) deallocate(bsubumns)
        if (allocated(bsubumnc)) deallocate(bsubumnc)
        if (allocated(bsubvmns)) deallocate(bsubvmns)
        if (allocated(bsubvmnc)) deallocate(bsubvmnc)
        if (allocated(bsubsmns)) deallocate(bsubsmns)
        if (allocated(bsubsmnc)) deallocate(bsubsmnc)
        if (allocated(bmns)) deallocate(bmns)
        if (allocated(bmnc)) deallocate(bmnc)
        if (allocated(lmns)) deallocate(lmns)
        if (allocated(lmnc)) deallocate(lmnc)
        if (allocated(rmns)) deallocate(rmns)
        if (allocated(rmnc)) deallocate(rmnc)
        if (allocated(zmns)) deallocate(zmns)
        if (allocated(zmnc)) deallocate(zmnc)
        if (allocated(R_c)) deallocate(R_c)
        if (allocated(R_s)) deallocate(R_s)
        if (allocated(Z_c)) deallocate(Z_c)
        if (allocated(Z_s)) deallocate(Z_s)
        if (allocated(L_c)) deallocate(L_c)
        if (allocated(L_s)) deallocate(L_s)
        if (allocated(B_V_sub_c_M)) deallocate(B_V_sub_c_M)
        if (allocated(B_V_sub_s_M)) deallocate(B_V_sub_s_M)
        if (allocated(B_V_c_H)) deallocate(B_V_c_H)
        if (allocated(B_V_s_H)) deallocate(B_V_s_H)
        if (allocated(jac_V_c_H)) deallocate(jac_V_c_H)
        if (allocated(jac_V_s_H)) deallocate(jac_V_s_H)
    end subroutine dealloc_VMEC

    ! Repack  variables  representing the  Fourier  composition  such as  R,  Z,
    ! lambda, ...  In VMEC these  are stored as  (1:mnmax, 1:ns) with  mnmax the
    ! total  number of  all modes.  Here  they are  to be  stored as  (0:mpol-1,
    ! -ntor:ntor, 1:ns), which  is valid due to the symmetry  of the modes (only
    ! one of either theta_V  or zeta_V has to be able to  change sign because of
    ! the (anti)-symmetry of the (co)sine.
    ! It is  possible that the  input variable is  not allocated. In  this case,
    ! output zeros
    ! [MPI] only global master
    !       (this is a precaution: only the global master should use it)
    function repack(var_VMEC,mnmax,n_r,mpol,ntor,xm,xn)
        use num_vars, only: glb_rank
        
        ! input / output
        integer, intent(in) :: mnmax, n_r, mpol, ntor
        real(dp), intent(in) :: xm(mnmax), xn(mnmax)
        real(dp), allocatable :: var_VMEC(:,:)
        real(dp) :: repack(0:mpol-1,-ntor:ntor,1:n_r)
            
        ! local variables
        integer :: mode, m, n
        
        if (allocated(var_VMEC) .and. glb_rank.eq.0) then                       ! only global rank
            ! check if the  values in xm and xn don't  exceed the maximum number
            ! of poloidal and toroidal modes (xm  and xn are of length mnmax and
            ! contain the pol/tor mode number)
            if (maxval(xm).gt.mpol .or. maxval(abs(xn)).gt.ntor) then
                call writo('WARNING: In repack, less modes are used than in the&
                    & VMEC format')
            end if
                
            repack = 0.0_dp
            ! copy the VMEC modes using the PB3D format
            do mode = 1,mnmax
                m = nint(xm(mode))
                n = nint(xn(mode))
                ! if the modes don't fit, cycle
                if (m.gt.mpol .or. abs(n).gt.ntor) then
                    call writo('WARNING: In repack, m > mpol or n > ntor!')
                    cycle
                end if
                
                repack(m,n,:) = var_VMEC(mode,:)
            end do
        else
            repack = 0.0_dp
        end if
    end function repack
end module VMEC_vars
