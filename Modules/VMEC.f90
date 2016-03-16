!------------------------------------------------------------------------------!
!   Operations and variables that concern the output of VMEC                   !
!------------------------------------------------------------------------------!
module VMEC
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use grid_vars, only: grid_type
    use num_vars, only: &
        &dp, max_str_ln, pi
    use read_wout_mod, only: read_wout_file, read_wout_deallocate, &            ! from LIBSTELL
        &is_asym_V => lasym, VMEC_version => version_, is_freeb_V => lfreeb, &  ! stellerator symmetry, version number, free boundary or not
        &n_r_VMEC => ns, mpol_V => mpol, ntor_V => ntor, xn, xm, &              ! n_r, mpol, ntor, xn, xm
        &mnmax_V => mnmax, nfp_V => nfp, &                                      ! mnmax, nfp
        &phi, phipf, &                                                          ! toroidal flux (FM), norm. deriv. of toroidal flux (FM)
        &iotaf, &                                                               ! rot. transf. (tor. flux) or saf. fac. (pol. flux) (FM)
        &presf, &                                                               ! pressure (FM)
        &lrfp, &                                                                ! whether or not the poloidal flux is used as radial variable
        &gam_V => gamma, &                                                      ! gamma in adiabatic law (not important here, incompressibility)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &bmns, bmnc, &                                                          ! magnitude of B (HM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf, &                                     ! max and min values of R, Z
        &gmnc, gmns                                                             ! Jacobian in VMEC coordinates
    
    implicit none
    private
    public read_VMEC, dealloc_VMEC, normalize_VMEC, &
        &calc_trigon_factors, fourier2real, &
        &R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mnmax_V, mpol_V, ntor_V, &
        &mn_V, pres_V, rot_t_V, is_asym_V, flux_t_V, Dflux_t_V, flux_p_V, &
        &Dflux_p_V, VMEC_version, gam_V, is_freeb_V, nfp_V, B_0_V
#if ldebug
    public B_V_sub_s, B_V_sub_c, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif

    ! local variables
    integer, allocatable :: mn_V(:,:)                                           ! m and n of modes
    real(dp) :: B_0_V                                                           ! the magnitude of B at the magnetic axis, theta = zeta = 0
    real(dp), allocatable :: flux_t_V(:), Dflux_t_V(:)                          ! toroidal flux and derivative
    real(dp), allocatable :: flux_p_V(:), Dflux_p_V(:)                          ! poloidal flux and derivative
    real(dp), allocatable :: pres_V(:)                                          ! pressure
    real(dp), allocatable :: rot_t_V(:)                                         ! rotational transform
    real(dp), allocatable :: R_V_c(:,:,:), R_V_s(:,:,:)                         ! Coeff. of R in (co)sine series (FM) and norm. deriv.
    real(dp), allocatable :: Z_V_c(:,:,:), Z_V_s(:,:,:)                         ! Coeff. of Z in (co)sine series (FM) and norm. deriv.
    real(dp), allocatable :: L_V_c(:,:,:), L_V_s(:,:,:)                         ! Coeff. of lambda in (co)sine series (HM) and norm. deriv.
    real(dp), allocatable :: B_V_sub_c(:,:,:), B_V_sub_s(:,:,:)                 ! Coeff. of B_i in (co)sine series (r,theta,phi) (FM)
    real(dp), allocatable :: B_V_c(:,:), B_V_s(:,:)                             ! Coeff. of magnitude of B (HM and FM)
    real(dp), allocatable :: jac_V_c(:,:), jac_V_s(:,:)                         ! Jacobian in VMEC coordinates (HM and FM)

contains
    ! Reads the VMEC equilibrium data
    ! [MPI] only master
    integer function read_VMEC(n_r_in,use_pol_flux_V) result(ierr)
        use utilities, only: conv_FHM, calc_int
        use num_vars, only: eq_name
        
        character(*), parameter :: rout_name = 'read_VMEC'
        
        ! input / output
        integer, intent(inout) :: n_r_in                                        ! nr. of normal points in input grid
        logical, intent(inout) :: use_pol_flux_V                                ! .true. if VMEC equilibrium is based on poloidal flux
        
        ! local variables
        integer :: id                                                           ! counters
        real(dp), allocatable :: L_c_H(:,:,:)                                   ! temporary HM variable
        real(dp), allocatable :: L_s_H(:,:,:)                                   ! temporary HM variable
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=8) :: flux_name                                           ! either poloidal or toroidal
#if ldebug
        integer :: kd                                                           ! counter
        real(dp), allocatable :: B_V_sub_c_M(:,:,:), B_V_sub_s_M(:,:,:)         ! Coeff. of B_i in (co)sine series (r,theta,phi) (FM, HM, HM)
        real(dp), allocatable :: B_V_c_H(:,:), B_V_s_H(:,:)                     ! Coeff. of magnitude of B (HM)
        real(dp), allocatable :: jac_V_c_H(:,:), jac_V_s_H(:,:)                 ! Jacobian in VMEC coordinates (HM)
#endif
        
        ! initialize ierr
        ierr = 0
        
        call writo('Reading data from VMEC output "' &
            &// trim(eq_name) // '"')
        call lvl_ud(1)
        
        ! read VMEC output using LIBSTELL
        call read_wout_file(eq_name,ierr)                                       ! read the VMEC file
        CHCKERR('Failed to read the VMEC file')
        
        ! set some variables
        n_r_in = n_r_VMEC
        allocate(flux_t_V(n_r_in))
        allocate(Dflux_t_V(n_r_in))
        allocate(flux_p_V(n_r_in))
        allocate(Dflux_p_V(n_r_in))
        allocate(rot_t_V(n_r_in))
        allocate(pres_V(n_r_in))
        use_pol_flux_V = lrfp
        flux_t_V = phi
        Dflux_t_V = phipf
        Dflux_p_V = iotaf*phipf
        ierr = calc_int(Dflux_p_V,1.0_dp/(n_r_in-1.0_dp),flux_p_V)              ! needs to be calculated here, on full input grid
        CHCKERR('')
        pres_V = presf
        rot_t_V = iotaf
        if (use_pol_flux_V) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        B_0_V = sum(bmnc(:,2))
        
        call writo('VMEC version is ' // trim(r2str(VMEC_version)))
        if (is_freeb_V) then
            call writo('Free boundary VMEC')
            err_msg = 'Free boundary VMEC is not yet supported by PB3D...'
            ierr = 1
            CHCKERR(err_msg)
        else
            call writo('Fixed boundary VMEC')
        end if
        if (is_asym_V) then
            call writo('No stellerator symmetry')
        else
            call writo('Stellerator symmetry')
        end if
        call writo('VMEC has '//trim(i2str(mpol_V))//' poloidal and '&
            &//trim(i2str(ntor_V))//' toroidal modes,')
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
        mn_V(:,1) = nint(xm)
        mn_V(:,2) = nint(xn)
        
        call lvl_ud(-1)
        call writo('Data from VMEC output successfully read')
        
        call writo('Saving VMEC output to PB3D')
        call lvl_ud(1)
        
        ! Allocate the Fourier coefficients and store the underived ones
        ! (only allocate for the underived indices 0)
        allocate(R_V_c(mnmax_V,n_r_in,0:0))
        allocate(R_V_s(mnmax_V,n_r_in,0:0))
        allocate(Z_V_c(mnmax_V,n_r_in,0:0))
        allocate(Z_V_s(mnmax_V,n_r_in,0:0))
        allocate(L_V_c(mnmax_V,n_r_in,0:0))
        allocate(L_V_s(mnmax_V,n_r_in,0:0))
        allocate(L_c_H(mnmax_V,n_r_in,0:0))
        allocate(L_s_H(mnmax_V,n_r_in,0:0))
        
        ! factors R_V_c,s; Z_V_c,s and L_C,s and HM varieties
        R_V_c(:,:,0) = rmnc
        Z_V_s(:,:,0) = zmns
        L_s_H(:,:,0) = lmns
        if (is_asym_V) then                                                     ! following only needed in asymmetric situations
            R_V_s(:,:,0) = rmns
            Z_V_c(:,:,0) = zmnc
            L_c_H(:,:,0) = lmnc
        else
            R_V_s(:,:,0) = 0._dp
            Z_V_c(:,:,0) = 0._dp
            L_c_H(:,:,0) = 0._dp
        end if
        
        ! conversion HM -> FM (L)
        do id = 1,mnmax_V
            ierr = conv_FHM(L_s_H(id,:,0),L_V_s(id,:,0),.false.)
            CHCKERR('')
            ierr = conv_FHM(L_c_H(id,:,0),L_V_c(id,:,0),.false.)
            CHCKERR('')
        end do
        
#if ldebug
        ! allocate helper variables
        allocate(B_V_sub_c_M(mnmax_V,n_r_in,3)); B_V_sub_c_M = 0._dp
        allocate(B_V_sub_s_M(mnmax_V,n_r_in,3)); B_V_sub_s_M = 0._dp
        allocate(B_V_c_H(mnmax_V,n_r_in)); B_V_c_H = 0._dp
        allocate(B_V_s_H(mnmax_V,n_r_in)); B_V_s_H = 0._dp
        allocate(jac_V_c_H(mnmax_V,n_r_in)); jac_V_c_H = 0._dp
        allocate(jac_V_s_H(mnmax_V,n_r_in))
        
        ! store in helper variables
        B_V_sub_s_M(:,:,1) = bsubsmns
        B_V_sub_c_M(:,:,2) = bsubumnc
        B_V_sub_c_M(:,:,3) = bsubvmnc
        B_V_c_H(:,:) = bmnc
        jac_V_c_H(:,:) = gmnc
        if (is_asym_V) then                                                 ! following only needed in asymmetric situations
            B_V_sub_c_M(:,:,1) = bsubsmnc
            B_V_sub_s_M(:,:,2) = bsubumns
            B_V_sub_s_M(:,:,3) = bsubvmns
            B_V_s_H(:,:) = bmns
            jac_V_s_H(:,:) = gmns
        else
            B_V_sub_c_M(:,:,1) = 0._dp
            B_V_sub_s_M(:,:,2) = 0._dp
            B_V_sub_s_M(:,:,3) = 0._dp
            B_V_s_H(:,:) = 0._dp
            jac_V_s_H(:,:) = 0._dp
        end if
        
        ! allocate FM variables
        allocate(B_V_sub_c(mnmax_V,n_r_in,3))
        allocate(B_V_sub_s(mnmax_V,n_r_in,3))
        allocate(B_V_c(mnmax_V,n_r_in))
        allocate(B_V_s(mnmax_V,n_r_in))
        allocate(jac_V_c(mnmax_V,n_r_in))
        allocate(jac_V_s(mnmax_V,n_r_in))
        
        ! conversion HM -> FM (B_V_sub, B_V, jac_V)
        do id = 1,mnmax_V
            do kd = 1,3
                ierr = conv_FHM(B_V_sub_c_M(id,:,kd),B_V_sub_c(id,:,kd),&
                    &.false.)
                CHCKERR('')
                ierr = conv_FHM(B_V_sub_s_M(id,:,kd),B_V_sub_s(id,:,kd),&
                    &.false.)
                CHCKERR('')
            end do
            ierr = conv_FHM(B_V_c_H(id,:),B_V_c(id,:),.false.)
            CHCKERR('')
            ierr = conv_FHM(B_V_s_H(id,:),B_V_s(id,:),.false.)
            CHCKERR('')
            ierr = conv_FHM(jac_V_c_H(id,:),jac_V_c(id,:),.false.)
            CHCKERR('')
            ierr = conv_FHM(jac_V_s_H(id,:),jac_V_s(id,:),.false.)
            CHCKERR('')
        end do
        
        ! deallocate helper variables
        deallocate(B_V_sub_c_M,B_V_sub_s_M)
        deallocate(B_V_c_H,B_V_s_H)
        deallocate(jac_V_c_H,jac_V_s_H)
#endif
        
        ! deallocate repacked variables
        call read_wout_deallocate()
        
        call lvl_ud(-1)
        call writo('Conversion complete')
    end function read_VMEC
    
    ! deallocates VMEC quantities that are not used anymore
    subroutine dealloc_VMEC
        deallocate(rot_t_V)
        deallocate(pres_V)
        deallocate(flux_t_V)
        deallocate(Dflux_t_V)
        deallocate(flux_p_V)
        deallocate(Dflux_p_V)
        deallocate(mn_V)
        deallocate(R_V_c)
        deallocate(R_V_s)
        deallocate(Z_V_c)
        deallocate(Z_V_s)
        deallocate(L_V_c)
        deallocate(L_V_s)
#if ldebug
        deallocate(B_V_sub_c)
        deallocate(B_V_sub_s)
        deallocate(B_V_c)
        deallocate(B_V_s)
        deallocate(jac_V_c)
        deallocate(jac_V_s)
#endif
    end subroutine dealloc_VMEC
    
    ! Calculate the trigoniometric cosine and sine factors on a grid (1:mnmax_V)
    ! at given 3D arrays for the (VMEC) E(quilibrium) angles theta_E and zeta_E.
    integer function calc_trigon_factors(theta,zeta,trigon_factors) &
        &result(ierr)
        
        character(*), parameter :: rout_name = 'calc_trigon_factors'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:,:)                                    ! poloidal angles in equilibrium coords.
        real(dp), intent(in) :: zeta(:,:,:)                                     ! toroidal angles in equilibrium coords.
        real(dp), intent(inout), allocatable :: trigon_factors(:,:,:,:,:)       ! trigonometric factor cosine and sine at these angles
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_ang_1, n_ang_2, n_r                                        ! sizes of 3D real output array
        real(dp), allocatable :: cos_theta(:,:,:,:)                             ! cos(m theta) for all m
        real(dp), allocatable :: sin_theta(:,:,:,:)                             ! sin(m theta) for all m
        real(dp), allocatable :: cos_zeta(:,:,:,:)                              ! cos(n theta) for all n
        real(dp), allocatable :: sin_zeta(:,:,:,:)                              ! sin(n theta) for all n
        integer :: id, m, n                                                     ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! set n_ang_1 and n_r
        n_ang_1 = size(theta,1)
        n_ang_2 = size(theta,2)
        n_r = size(theta,3)
        
        ! tests
        if (size(zeta,1).ne.n_ang_1 .or. size(zeta,2).ne.n_ang_2 .or. &
            &size(zeta,3).ne.n_r) then
            ierr = 1
            err_msg = 'theta and zeta need to have the same size'
            CHCKERR(err_msg)
        end if
        
        ! setup cos_theta, sin_theta, cos_zeta and sin_zeta
        allocate(cos_theta(0:mpol_V-1,n_ang_1,n_ang_2,n_r))
        allocate(sin_theta(0:mpol_V-1,n_ang_1,n_ang_2,n_r))
        allocate(cos_zeta(-ntor_V:ntor_V,n_ang_1,n_ang_2,n_r))
        allocate(sin_zeta(-ntor_V:ntor_V,n_ang_1,n_ang_2,n_r))
        
        do m = 0,mpol_V-1
            cos_theta(m,:,:,:) = cos(m*theta)
            sin_theta(m,:,:,:) = sin(m*theta)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! MAYBE IN CALC_TRIGON_FACTORS, MAYBE YOU HAVE TO USE - ZETA???????  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n = -ntor_V,ntor_V
            cos_zeta(n,:,:,:) = cos(n*zeta)
            sin_zeta(n,:,:,:) = sin(n*zeta)
        end do
        
        ! initialize trigon_factors
        allocate(trigon_factors(mnmax_V,n_ang_1,n_ang_2,n_r,2))
        trigon_factors = 0.0_dp
        
        ! calculate cos(m theta - n zeta) = cos(m theta) cos(n zeta) + 
        !   sin(m theta) sin(n zeta) and
        ! sin(m theta - n zeta) = sin(m theta) cos(n zeta) -
        !   cos(m theta) sin(n zeta)
        do id = 1,mnmax_V
            trigon_factors(id,:,:,:,1) = &
                &cos_theta(mn_V(id,1),:,:,:)*cos_zeta(mn_V(id,2)/nfp_V,:,:,:) + &
                &sin_theta(mn_V(id,1),:,:,:)*sin_zeta(mn_V(id,2)/nfp_V,:,:,:)
            trigon_factors(id,:,:,:,2) = &
                &sin_theta(mn_V(id,1),:,:,:)*cos_zeta(mn_V(id,2)/nfp_V,:,:,:) - &
                &cos_theta(mn_V(id,1),:,:,:)*sin_zeta(mn_V(id,2)/nfp_V,:,:,:)
        end do
        
        ! deallocate variables
        deallocate(cos_theta,sin_theta)
        deallocate(cos_zeta,sin_zeta)
    end function calc_trigon_factors
    
    ! Inverse Fourier transformation, from VMEC. Also calculates the poloidal or
    ! toroidal  derivatives  in  VMEC  coords., as  indicated  by  the  variable
    ! deriv(2).
    ! (Normal derivative  is done on the variables in  Fourier space, and should
    ! be provided here in var_fourier_i if needed).
    integer function fourier2real(var_fourier_c,var_fourier_s,trigon_factors,&
        &var_real,deriv) result(ierr)
        
        character(*), parameter :: rout_name = 'fourier2real'
        
        ! input / output
        real(dp), intent(in) :: var_fourier_c(:,:)                              ! cos factor of variable in Fourier space
        real(dp), intent(in) :: var_fourier_s(:,:)                              ! sin factor of variable in Fourier space
        real(dp), intent(in) :: trigon_factors(:,:,:,:,:)                       ! trigonometric factor cosine and sine at these angles
        real(dp), intent(inout) :: var_real(:,:,:)                              ! variable in real space
        integer, intent(in), optional :: deriv(2)                               ! optional derivatives in angular coordinates
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_ang_1, n_r, n_ang_2                                        ! sizes of 3D real output array
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: fac_cos(:), fac_sin(:)                         ! factor in front of cos and sin, after taking derivatives
        real(dp), allocatable :: fac_trigon_temp(:)                             ! temporary variable that holds fac_cos or fac_sin
        
        ! initialize ierr
        ierr = 0
        
        ! set n_ang_1 and n_r
        n_ang_1 = size(trigon_factors,2)
        n_ang_2 = size(trigon_factors,3)
        n_r = size(trigon_factors,4)
        
        ! tests
        if (size(trigon_factors,5).ne.2) then
            ierr = 1
            err_msg = 'trigon_factors needs to contain sines and cosines'
            CHCKERR(err_msg)
        end if
        if (size(trigon_factors,1).ne.mnmax_V) then
            ierr = 1
            err_msg = 'trigon_factors needs to be defined for the right number &
                &of modes'
            CHCKERR(err_msg)
        end if
        if (size(var_fourier_c,2).ne.n_r .or. &
            &size(var_fourier_s,2).ne.n_r) then
            ierr = 1
            err_msg = 'var_fourier_c and _s need to have the right number of &
                &normal points'
            CHCKERR(err_msg)
        end if
        
        ! initialize fac_cos, fac_sin and fac_trigon_temp
        allocate(fac_cos(n_r),fac_sin(n_r))                                     ! factor in front of cos and sin, after taking derivatives
        allocate(fac_trigon_temp(n_r))                                          ! temporary variable that holds fac_cos or fac_sin
        
        ! initialize
        var_real = 0.0_dp
        
        ! sum over modes
        do id = 1,mnmax_V
            ! initialize factors in front of cos and sin
            fac_cos = var_fourier_c(id,:)
            fac_sin = var_fourier_s(id,:)
            
            ! angular derivatives
            if (present(deriv)) then 
                ! apply possible poloidal derivatives
                do kd = 1,deriv(1)
                    fac_trigon_temp = - mn_V(id,1) * fac_cos
                    fac_cos = mn_V(id,1) * fac_sin
                    fac_sin = fac_trigon_temp
                end do
                ! apply possible toroidal derivatives
                do kd = 1,deriv(2)
                    fac_trigon_temp = mn_V(id,2) * fac_cos
                    fac_cos = - mn_V(id,2) * fac_sin
                    fac_sin = fac_trigon_temp
                end do
            end if
            
            ! sum
            do kd = 1,n_r
                var_real(:,:,kd) = var_real(:,:,kd) + &
                    &fac_cos(kd)*trigon_factors(id,:,:,kd,1) + &
                    &fac_sin(kd)*trigon_factors(id,:,:,kd,2)
            end do
        end do
    end function fourier2real
    
    ! Normalizes VMEC input
    ! Note  that  the normal  VMEC coordinate  runs from  0 to  1, whatever  the
    ! normalization.
    subroutine normalize_VMEC
        use eq_vars, only: pres_0, psi_0, R_0
#if ldebug
        use  eq_vars, only: B_0
#endif
        
        ! scale the VMEC quantities
        pres_V = pres_V/pres_0
        flux_t_V = flux_t_V/psi_0
        Dflux_t_V = Dflux_t_V/psi_0
        flux_p_V = flux_p_V/psi_0
        Dflux_p_V = Dflux_p_V/psi_0
        R_V_c = R_V_c/R_0
        R_V_s = R_V_s/R_0
        Z_V_c = Z_V_c/R_0
        Z_V_s = Z_V_s/R_0
        L_V_c = L_V_c
        L_V_s = L_V_s
#if ldebug
        B_V_sub_s = B_V_sub_s/(R_0*B_0)
        B_V_sub_c = B_V_sub_c/(R_0*B_0)
        B_V_c = B_V_c/B_0
        B_V_s = B_V_s/B_0
        jac_V_c = jac_V_c/(R_0**3)
        jac_V_s = jac_V_s/(R_0**3)
#endif
    end subroutine normalize_VMEC
end module VMEC
