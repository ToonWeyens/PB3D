!------------------------------------------------------------------------------!
!> Operations that concern the output of VMEC.
!------------------------------------------------------------------------------!
module VMEC_ops
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: &
        &dp, max_str_ln, pi
    use read_wout_mod, only: read_wout_file, read_wout_deallocate, &            ! from LIBSTELL
        &n_r_VMEC => ns, &                                                      ! n_r
        &xn, xm, &                                                              ! xn, xm
        &lrfp, &                                                                ! whether or not the poloidal flux is used as radial variable
        &phi, phipf, &                                                          ! toroidal flux (FM), norm. deriv. of toroidal flux (FM)
        &iotaf, &                                                               ! rot. transf. (tor. flux) or saf. fac. (pol. flux) (FM)
        &presf, &                                                               ! pressure (FM)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &jcuru, jcurv, &                                                        ! poloidal (FM) and toroidal (FM) current
        &bmns, bmnc, &                                                          ! magnitude of B (HM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &gmnc, gmns                                                             ! Jacobian in VMEC coordinates
    use VMEC_vars
    
    implicit none
    private
    public read_VMEC, normalize_VMEC

contains
    !> Reads the VMEC equilibrium data.
    !!
    !! \return ierr
    integer function read_VMEC(n_r_in,use_pol_flux_V) result(ierr)
        use num_utilities, only: calc_int, spline
        use num_vars, only: eq_name, max_deriv, norm_disc_prec_eq
        
        character(*), parameter :: rout_name = 'read_VMEC'
        
        ! input / output
        integer, intent(inout) :: n_r_in                                        !< nr. of normal points in input grid
        logical, intent(inout) :: use_pol_flux_V                                !< .true. if VMEC equilibrium is based on poloidal flux
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: r_V(:)                                         ! normal coordinate
        real(dp), allocatable :: L_c_H(:,:,:)                                   ! temporary HM variable
        real(dp), allocatable :: L_s_H(:,:,:)                                   ! temporary HM variable
        real(dp), allocatable :: jac_c_H(:,:,:)                                 ! temporary HM variable
        real(dp), allocatable :: jac_s_H(:,:,:)                                 ! temporary HM variable
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=8) :: flux_name                                           ! either poloidal or toroidal
        real(dp), allocatable :: B_V_sub_c_M(:,:,:), B_V_sub_s_M(:,:,:)         ! Coeff. of B_i in (co)sine series (r,theta,phi) (FM, HM, HM)
#if ldebug
        real(dp), allocatable :: B_V_c_H(:,:), B_V_s_H(:,:)                     ! Coeff. of magnitude of B (HM)
#endif
        
        ! initialize ierr
        ierr = 0
        
        call writo('Reading data from VMEC output "' &
            &// trim(eq_name) // '"')
        call lvl_ud(1)
        
        ! read VMEC output using LIBSTELL
        if (scan(eq_name, ".") .ne. scan(eq_name, ".", .true.)) &
            &call writo('There seems to be a dot "."  in "'//trim(eq_name)//&
            &'". This can create problems.', alert=.true.)
        call read_wout_file(eq_name,ierr)                                       ! read the VMEC file
        CHCKERR('Failed to read the VMEC file')
        
        ! set some variables
        n_r_in = n_r_VMEC
        allocate(flux_t_V(n_r_in,0:max_deriv+1))
        allocate(flux_p_V(n_r_in,0:max_deriv+1))
        allocate(rot_t_V(n_r_in,0:max_deriv+1))
        allocate(q_saf_V(n_r_in,0:max_deriv+1))
        allocate(pres_V(n_r_in,0:max_deriv+1))
        use_pol_flux_V = lrfp
        flux_t_V(:,0) = phi
        flux_t_V(:,1) = phipf
        flux_p_V(:,1) = iotaf*phipf
        ierr = calc_int(flux_p_V(:,1),1.0_dp/(n_r_in-1.0_dp),flux_p_V(:,0))
        CHCKERR('')
        pres_V(:,0) = presf
        rot_t_V(:,0) = iotaf
        q_saf_V(:,0) = 1._dp/iotaf
        if (use_pol_flux_V) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        B_0_V = sum(1.5*bmnc(:,2)-0.5_dp*bmnc(:,3))                             ! convert to full mesh
        
        call writo('VMEC version is ' // trim(r2str(VMEC_version)))
        if (is_freeb_V) then
            call writo('Free boundary VMEC')
        else
            call writo('Fixed boundary VMEC')
        end if
        if (is_asym_V) then
            call writo('No stellarator symmetry')
        else
            call writo('Stellarator symmetry')
        end if
        call writo('VMEC has '//trim(i2str(mpol_V))//' poloidal and '&
            &//trim(i2str(ntor_V))//' toroidal modes,')
        if (ntor_V.gt.0) call writo('(basic toroidal field period NFP = '//&
            &trim(i2str(nfp_V))//')')
        call writo('defined on '//trim(i2str(n_r_in))//' '//flux_name//&
            &' flux surfaces')
        
        call writo('Running tests')
        call lvl_ud(1)
        if (mnmax_V.ne.((ntor_V+1)+(2*ntor_V+1)*(mpol_V-1))) then               ! for mpol_V = 0, only ntor_V+1 values needed
            err_msg = 'Inconsistency in ntor_V, mpol_V and mnmax_V'
            ierr = 1
            CHCKERR(err_msg)
        end if
        call lvl_ud(-1)
        
        call writo('Updating variables')
        allocate(mn_V(mnmax_V,2))
        mn_V(:,1) = nint(xm(1:mnmax_V))
        mn_V(:,2) = nint(xn(1:mnmax_V))
        
        call lvl_ud(-1)
        call writo('Data from VMEC output successfully read')
        
        call writo('Saving VMEC output to PB3D')
        call lvl_ud(1)
        
        ! Allocate the Fourier coefficients
        allocate(R_V_c(mnmax_V,n_r_in,0:max_deriv+1))
        allocate(R_V_s(mnmax_V,n_r_in,0:max_deriv+1))
        allocate(Z_V_c(mnmax_V,n_r_in,0:max_deriv+1))
        allocate(Z_V_s(mnmax_V,n_r_in,0:max_deriv+1))
        allocate(L_V_c(mnmax_V,n_r_in,0:max_deriv+1))
        allocate(L_V_s(mnmax_V,n_r_in,0:max_deriv+1))
        allocate(jac_V_c(mnmax_V,n_r_in,0:max_deriv))
        allocate(jac_V_s(mnmax_V,n_r_in,0:max_deriv))
        allocate(L_c_H(mnmax_V,n_r_in,0:0))
        allocate(L_s_H(mnmax_V,n_r_in,0:0))
        allocate(jac_c_H(mnmax_V,n_r_in,0:0))
        allocate(jac_s_H(mnmax_V,n_r_in,0:0))
        
        ! factors R_V_c,s; Z_V_c,s and L_C,s and HM varieties
        R_V_c(:,:,0) = rmnc(1:mnmax_V,:)
        Z_V_s(:,:,0) = zmns(1:mnmax_V,:)
        L_s_H(:,:,0) = lmns(1:mnmax_V,:)
        jac_c_H(:,:,0) = gmnc(1:mnmax_V,:)
        if (is_asym_V) then                                                     ! following only needed in asymmetric situations
            R_V_s(:,:,0) = rmns(1:mnmax_V,:)
            Z_V_c(:,:,0) = zmnc(1:mnmax_V,:)
            L_c_H(:,:,0) = lmnc(1:mnmax_V,:)
            jac_s_H(:,:,0) = gmns(1:mnmax_V,:)
        else
            R_V_s = 0._dp
            Z_V_c = 0._dp
            L_c_H = 0._dp
            jac_s_H = 0._dp
        end if
        
        ! calculate data for normal derivatives  with the toroidal (or poloidal)
        ! flux, normalized wrt. to the  maximum flux, equidistantly, so the step
        ! size is 1/(n_r_in-1).
        allocate(r_V(n_r_in))
        r_V = [((kd-1._dp)/(n_r_in-1),kd=1,n_r_in)]
        
        ! up to order 2, full mesh
        do kd = 1,2
            ierr = spline(r_V,q_saf_V(:,0),r_V,q_saf_V(:,kd),&
                &ord=norm_disc_prec_eq,deriv=kd)
            CHCKERR('')
            ierr = spline(r_V,rot_t_V(:,0),r_V,rot_t_V(:,kd),&
                &ord=norm_disc_prec_eq,deriv=kd)
            CHCKERR('')
            ierr = spline(r_V,pres_V(:,0),r_V,pres_V(:,kd),&
                &ord=norm_disc_prec_eq,deriv=kd)
            CHCKERR('')
            ierr = spline(r_V,flux_t_V(:,0),r_V,flux_t_V(:,kd),&
                &ord=norm_disc_prec_eq,deriv=kd)
            CHCKERR('')
            ierr = spline(r_V,flux_p_V(:,0),r_V,flux_p_V(:,kd),&
                &ord=norm_disc_prec_eq,deriv=kd)
            CHCKERR('')
        end do
        
        ! up to order 2, half mesh: also extrapolate
        do kd = 0,2
            do id = 1,mnmax_V
                ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                    &jac_c_H(id,2:n_r_in,0),r_V,jac_V_c(id,:,kd),&
                    &ord=norm_disc_prec_eq,deriv=kd,extrap=.true.)
                CHCKERR('')
                if (is_asym_V) then
                    ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                        &jac_s_H(id,2:n_r_in,0),r_V,jac_V_s(id,:,kd),&
                        &ord=norm_disc_prec_eq,deriv=kd,extrap=.true.)
                    CHCKERR('')
                end if
            end do
        end do
        
        ! up to order 2, full mesh
        do kd = 1,3
            do id = 1,mnmax_V
                ! full mesh
                ierr = spline(r_V,R_V_c(id,:,0),r_V,R_V_c(id,:,kd),&
                    &ord=norm_disc_prec_eq,deriv=kd)
                CHCKERR('')
                ierr = spline(r_V,Z_V_s(id,:,0),r_V,Z_V_s(id,:,kd),&
                    &ord=norm_disc_prec_eq,deriv=kd)
                CHCKERR('')
                if (is_asym_V) then
                    ierr = spline(r_V,R_V_s(id,:,0),r_V,R_V_s(id,:,kd),&
                        &ord=norm_disc_prec_eq,deriv=kd)
                    CHCKERR('')
                    ierr = spline(r_V,Z_V_c(id,:,0),r_V,Z_V_c(id,:,kd),&
                        &ord=norm_disc_prec_eq,deriv=kd)
                    CHCKERR('')
                end if
            end do
        end do
        
        ! up to order 2, half mesh: also extrapolate
        do kd = 0,3
            do id = 1,mnmax_V
                ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                    &L_s_H(id,2:n_r_in,0),r_V,L_V_s(id,:,kd),&
                    &ord=norm_disc_prec_eq,deriv=kd,extrap=.true.)
                CHCKERR('')
                if (is_asym_V) then
                    ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                        &L_c_H(id,2:n_r_in,0),r_V,L_V_c(id,:,kd),&
                        &ord=norm_disc_prec_eq,deriv=kd,extrap=.true.)
                    CHCKERR('')
                end if
            end do
        end do
        
        ! These are only necessary if the compiler searches for uninitialized
        ! memory
        flux_t_V(:,3:) = 0._dp
        flux_p_V(:,3:) = 0._dp
        q_saf_V(:,3:) = 0._dp
        rot_t_V(:,3:) = 0._dp
        pres_V(:,3:) = 0._dp
        
        ! allocate helper variables
        allocate(B_V_sub_c_M(mnmax_V,n_r_in,3)); B_V_sub_c_M = 0._dp
        allocate(B_V_sub_s_M(mnmax_V,n_r_in,3)); B_V_sub_s_M = 0._dp
#if ldebug
        allocate(B_V_c_H(mnmax_V,n_r_in)); B_V_c_H = 0._dp
        allocate(B_V_s_H(mnmax_V,n_r_in)); B_V_s_H = 0._dp
#endif
        
        ! store in helper variables
        B_V_sub_s_M(:,:,1) = bsubsmns(1:mnmax_V,:)
        B_V_sub_c_M(:,:,2) = bsubumnc(1:mnmax_V,:)
        B_V_sub_c_M(:,:,3) = bsubvmnc(1:mnmax_V,:)
        if (is_asym_V) then                                                     ! following only needed in asymmetric situations
            B_V_sub_c_M(:,:,1) = bsubsmnc(1:mnmax_V,:)
            B_V_sub_s_M(:,:,2) = bsubumns(1:mnmax_V,:)
            B_V_sub_s_M(:,:,3) = bsubvmns(1:mnmax_V,:)
        else
            B_V_sub_c_M(:,:,1) = 0._dp
            B_V_sub_s_M(:,:,2) = 0._dp
            B_V_sub_s_M(:,:,3) = 0._dp
        end if
#if ldebug
        B_V_c_H(:,:) = bmnc(1:mnmax_V,:)
        if (is_asym_V) then                                                     ! following only needed in asymmetric situations
            B_V_s_H(:,:) = bmns(1:mnmax_V,:)
        else
            B_V_s_H(:,:) = 0._dp
        end if
#endif
        
        ! allocate FM variables
        allocate(B_V_sub_c(mnmax_V,n_r_in,3))
        allocate(B_V_sub_s(mnmax_V,n_r_in,3))
#if ldebug
        allocate(B_V_c(mnmax_V,n_r_in))
        allocate(B_V_s(mnmax_V,n_r_in))
        allocate(J_V_sup_int(n_r_in,2))
#endif
        
        ! half mesh: extrapolate
        do id = 1,mnmax_V
            do kd = 2,3
                ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                    &B_V_sub_c_M(id,2:n_r_in,kd),r_V,B_V_sub_c(id,:,kd),&
                    &ord=norm_disc_prec_eq,deriv=0,extrap=.true.)
                CHCKERR('')
                ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                    &B_V_sub_s_M(id,2:n_r_in,kd),r_V,B_V_sub_s(id,:,kd),&
                    &ord=norm_disc_prec_eq,deriv=0,extrap=.true.)
                CHCKERR('')
            end do
        end do
#if ldebug
        do id = 1,mnmax_V
            ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                &B_V_c_H(id,2:n_r_in),r_V,B_V_c(id,:),&
                &ord=norm_disc_prec_eq,deriv=0,extrap=.true.)
            CHCKERR('')
            ierr = spline(-0.5_dp/n_r_in+r_V(2:n_r_in),&
                &B_V_s_H(id,2:n_r_in),r_V,B_V_s(id,:),&
                &ord=norm_disc_prec_eq,deriv=0,extrap=.true.)
            CHCKERR('')
        end do
        J_V_sup_int(:,1) = jcuru
        J_V_sup_int(:,2) = jcurv
#endif
        
        ! deallocate repacked variables
        call read_wout_deallocate()
        
        call lvl_ud(-1)
        call writo('Conversion complete')
    end function read_VMEC
    
    !> Normalizes VMEC input.
    !!
    !! \note  The  normal  VMEC  coordinate  runs from  0  to  1,  whatever  the
    !! normalization.
    subroutine normalize_VMEC
        use eq_vars, only: pres_0, psi_0, R_0, B_0
        
        ! scale the VMEC quantities
        pres_V = pres_V/pres_0
        flux_t_V = flux_t_V/psi_0
        flux_p_V = flux_p_V/psi_0
        R_V_c = R_V_c/R_0
        R_V_s = R_V_s/R_0
        Z_V_c = Z_V_c/R_0
        Z_V_s = Z_V_s/R_0
        L_V_c = L_V_c
        L_V_s = L_V_s
        jac_V_c = jac_V_c/(R_0**3)
        jac_V_s = jac_V_s/(R_0**3)
        B_V_sub_s = B_V_sub_s/(R_0*B_0)
        B_V_sub_c = B_V_sub_c/(R_0*B_0)
#if ldebug
        B_V_c = B_V_c/B_0
        B_V_s = B_V_s/B_0
        J_V_sup_int = J_V_sup_int*psi_0/pres_0
#endif
    end subroutine normalize_VMEC
end module VMEC_ops
