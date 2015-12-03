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
        &lasym, VMEC_version => version_, lfreeb, &                             ! stellerator symmetry, version number, free boundary or not
        &n_r_VMEC => ns, mpol, ntor, xn, xm, mnmax, nfp, &                      ! mpol, ntor = # modes
        &flux_t_V => phi, Dflux_t_V => phipf, &                                 ! toroidal flux (FM), norm. deriv. of toroidal flux (FM)
        &rot_t_V => iotaf, &                                                    ! rotational transform = 1/q (FM)
        &pres_V => presf, &                                                     ! pressure (FM)
        &lrfp, &                                                                ! whether or not the poloidal flux is used as radial variable
        &gam => gamma, &                                                        ! gamma in adiabatic law (NOT important here because of incompressibility)
        &bsubumns, bsubumnc, bsubvmns, bsubvmnc, bsubsmns, bsubsmnc, &          ! B_theta (HM), B_zeta (HM), B_r (FM)
        &bmns, bmnc, &                                                          ! magnitude of B (HM)
        &lmns, lmnc, rmns, rmnc, zmns, zmnc, &                                  ! lambda (HM), R (FM), Z(FM)
        &rmax_surf, rmin_surf, zmax_surf, &                                     ! max and min values of R, Z
        &gmnc, gmns                                                             ! Jacobian in VMEC coordinates
    
    implicit none
    private
    public read_VMEC, dealloc_VMEC, repack, normalize_VMEC, &
        &calc_trigon_factors, fourier2real, &
        &mpol, ntor, nfp, R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, &
        &pres_V, rot_t_V, lasym, flux_t_V, Dflux_t_V, VMEC_version, gam, lfreeB
#if ldebug
    public B_V_sub_s, B_V_sub_c, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif

    real(dp), allocatable :: R_V_c(:,:,:,:), R_V_s(:,:,:,:)                     ! Coeff. of R in (co)sine series (FM) and norm. deriv.
    real(dp), allocatable :: Z_V_c(:,:,:,:), Z_V_s(:,:,:,:)                     ! Coeff. of Z in (co)sine series (FM) and norm. deriv.
    real(dp), allocatable :: L_V_c(:,:,:,:), L_V_s(:,:,:,:)                     ! Coeff. of lambda in (co)sine series (HM) and norm. deriv.
    real(dp), allocatable :: B_V_sub_c(:,:,:,:), B_V_sub_s(:,:,:,:)             ! Coeff. of B_i in (co)sine series (r,theta,phi) (FM)
    real(dp), allocatable :: B_V_c(:,:,:), B_V_s(:,:,:)                         ! Coeff. of magnitude of B (HM and FM)
    real(dp), allocatable :: jac_V_c(:,:,:), jac_V_s(:,:,:)                     ! Jacobian in VMEC coordinates (HM and FM)

contains
    ! Reads the VMEC equilibrium data
    ! [MPI] only global master
    integer function read_VMEC(n_r_eq,use_pol_flux_V) result(ierr)
        use utilities, only: conv_FHM
        use num_vars, only: max_deriv, eq_i, eq_name
#if ldebug
        use num_vars, only: ltest
#endif
        character(*), parameter :: rout_name = 'read_VMEC'
        
        ! input / output
        integer, intent(inout) :: n_r_eq                                        ! nr. of normal points in equilibrium grid
        logical, intent(inout) :: use_pol_flux_V                                ! .true. if VMEC equilibrium is based on poloidal flux
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: L_c_H(:,:,:,:)                                 ! temporary HM variable
        real(dp), allocatable :: L_s_H(:,:,:,:)                                 ! temporary HM variable
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: B_V_sub_c_M(:,:,:,:), B_V_sub_s_M(:,:,:,:)     ! Coeff. of B_i in (co)sine series (r,theta,phi) (FM, HM, HM)
        real(dp), allocatable :: B_V_c_H(:,:,:), B_V_s_H(:,:,:)                 ! Coeff. of magnitude of B (HM)
        real(dp), allocatable :: jac_V_c_H(:,:,:), jac_V_s_H(:,:,:)             ! Jacobian in VMEC coordinates (HM)
        
        ! initialize ierr
        ierr = 0
        
        call writo('Reading data from VMEC output "' &
            &// trim(eq_name) // '"')
        call lvl_ud(1)
        
        ! read VMEC output using LIBSTELL
        call read_wout_file(eq_i,ierr)                                          ! read the VMEC file number
        CHCKERR('Failed to read the VMEC file')
        
        ! close the VMEC output file
        close(eq_i)
        
        ! set some variables
        n_r_eq = n_r_VMEC
        use_pol_flux_V = lrfp
        
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
        if (mnmax.ne.((ntor+1)+(2*ntor+1)*(mpol-1))) then                       ! are ntor, mpol and mnmax what is expected?
            err_msg = 'Inconsistency in ntor, mpol and mnmax'
            ierr = 1
            CHCKERR(err_msg)
        end if
        call lvl_ud(-1)
        
        call writo('Updating variables')
        xn = xn/nfp                                                             ! so we have xn excluding nfp
        
        call lvl_ud(-1)
        call writo('Data from VMEC output successfully read')
        
        call writo('Converting VMEC output to PB3D format')
        call lvl_ud(1)
        
        ! Allocate and repack the Fourier coefficients to translate them for use
        ! in this code
        allocate(R_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
        allocate(Z_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
        allocate(L_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
        allocate(L_s_H(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
        !if (lasym) then                                                         ! following only needed in asymmetric situations
            allocate(R_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
            allocate(Z_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
            allocate(L_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
            allocate(L_c_H(0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv+1))
        !end if
        
        ! factors R_V_c,s; Z_V_c,s and L_C,s and HM varieties
        R_V_c(:,:,:,0) = repack(rmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
        Z_V_s(:,:,:,0) = repack(zmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
        L_s_H(:,:,:,0) = repack(lmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
        !if (lasym) then                                                         ! following only needed in asymmetric situations
            R_V_s(:,:,:,0) = repack(rmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            Z_V_c(:,:,:,0) = repack(zmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            L_c_H(:,:,:,0) = repack(lmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
        !end if
        
        ! conversion HM -> FM (L)
        do jd = -ntor,ntor
            do id = 0,mpol-1
                ierr = conv_FHM(L_s_H(id,jd,:,0),L_V_s(id,jd,:,0),.false.)
                CHCKERR('')
                !if (lasym) then                                                ! following only needed in asymmetric situations
                    ierr = conv_FHM(L_c_H(id,jd,:,0),L_V_c(id,jd,:,0),&
                        &.false.)
                    CHCKERR('')
                !end if
            end do
        end do
        
#if ldebug
        ! for tests
        if (ltest) then
            ! allocate helper variables
            allocate(B_V_sub_c_M(0:mpol-1,-ntor:ntor,1:n_r_eq,3))
            allocate(B_V_sub_s_M(0:mpol-1,-ntor:ntor,1:n_r_eq,3))
            allocate(B_V_c_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
            allocate(B_V_s_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
            allocate(jac_V_c_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
            allocate(jac_V_s_H(0:mpol-1,-ntor:ntor,1:n_r_eq))
            
            ! store in helper variables
            B_V_sub_c_M(:,:,:,1) = repack(bsubsmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_sub_s_M(:,:,:,1) = repack(bsubsmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_sub_c_M(:,:,:,2) = repack(bsubumnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_sub_s_M(:,:,:,2) = repack(bsubumns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_sub_c_M(:,:,:,3) = repack(bsubvmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_sub_s_M(:,:,:,3) = repack(bsubvmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_c_H(:,:,:) = repack(bmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            B_V_s_H(:,:,:) = repack(bmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            jac_V_c_H(:,:,:) = repack(gmnc,mnmax,n_r_eq,mpol,ntor,xm,xn)
            jac_V_s_H(:,:,:) = repack(gmns,mnmax,n_r_eq,mpol,ntor,xm,xn)
            
            ! allocate FM variables
            allocate(B_V_sub_c(0:mpol-1,-ntor:ntor,1:n_r_eq,3))
            allocate(B_V_sub_s(0:mpol-1,-ntor:ntor,1:n_r_eq,3))
            allocate(B_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq))
            allocate(B_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq))
            allocate(jac_V_c(0:mpol-1,-ntor:ntor,1:n_r_eq))
            allocate(jac_V_s(0:mpol-1,-ntor:ntor,1:n_r_eq))
            
            ! conversion HM -> FM (B_V_sub, B_V, jac_V)
            do jd = -ntor,ntor
                do id = 0,mpol-1
                    do kd = 1,3
                        ierr = conv_FHM(B_V_sub_c_M(id,jd,:,kd),&
                            &B_V_sub_c(id,jd,:,kd),.false.)
                        CHCKERR('')
                        ierr = conv_FHM(B_V_sub_s_M(id,jd,:,kd),&
                            &B_V_sub_s(id,jd,:,kd),.false.)
                        CHCKERR('')
                    end do
                    ierr = conv_FHM(B_V_c_H(id,jd,:),B_V_c(id,jd,:),.false.)
                    CHCKERR('')
                    ierr = conv_FHM(B_V_s_H(id,jd,:),B_V_s(id,jd,:),.false.)
                    CHCKERR('')
                    ierr = conv_FHM(jac_V_c_H(id,jd,:),jac_V_c(id,jd,:),.false.)
                    CHCKERR('')
                    ierr = conv_FHM(jac_V_s_H(id,jd,:),jac_V_s(id,jd,:),.false.)
                    CHCKERR('')
                end do
            end do
            
            ! deallocate helper variables
            deallocate(B_V_sub_c_M,B_V_sub_s_M)
            deallocate(B_V_c_H,B_V_s_H)
            deallocate(jac_V_c_H,jac_V_s_H)
        end if
#endif
        
        ! deallocate repacked variables
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
        
        call lvl_ud(-1)
        call writo('Conversion complete')
    end function read_VMEC
    
    ! deallocates VMEC quantities that are not used anymore
    subroutine dealloc_VMEC
        deallocate(flux_t_V,Dflux_t_V)
        deallocate(rot_t_V)
        deallocate(pres_V)
        if (allocated(R_V_c)) deallocate(R_V_c)
        if (allocated(R_V_s)) deallocate(R_V_s)
        if (allocated(Z_V_c)) deallocate(Z_V_c)
        if (allocated(Z_V_s)) deallocate(Z_V_s)
        if (allocated(L_V_c)) deallocate(L_V_c)
        if (allocated(L_V_s)) deallocate(L_V_s)
#if ldebug
        if (allocated(B_V_sub_c)) deallocate(B_V_sub_c)
        if (allocated(B_V_sub_s)) deallocate(B_V_sub_s)
        if (allocated(B_V_c)) deallocate(B_V_c)
        if (allocated(B_V_s)) deallocate(B_V_s)
        if (allocated(jac_V_c)) deallocate(jac_V_c)
        if (allocated(jac_V_s)) deallocate(jac_V_s)
#endif
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
        use num_vars, only: rank
        
        ! input / output
        integer, intent(in) :: mnmax, n_r, mpol, ntor
        real(dp), intent(in) :: xm(mnmax), xn(mnmax)
        real(dp), allocatable :: var_VMEC(:,:)
        real(dp) :: repack(0:mpol-1,-ntor:ntor,1:n_r)
            
        ! local variables
        integer :: mode, m, n
        
        if (allocated(var_VMEC) .and. rank.eq.0) then                           ! only global rank
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
    
    ! Calculate the trigoniometric cosine and  sine factors on a grid (0:mpol-1,
    ! -ntor:ntor) at given 3D arrays for the (VMEC) E(quilibrium) angles theta_E
    ! and zeta_E.
    integer function calc_trigon_factors(theta,zeta,trigon_factors) &
        &result(ierr)
        
        character(*), parameter :: rout_name = 'calc_trigon_factors'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:,:)                                    ! poloidal angles in equilibrium coords.
        real(dp), intent(in) :: zeta(:,:,:)                                     ! toroidal angles in equilibrium coords.
        real(dp), intent(inout), allocatable :: trigon_factors(:,:,:,:,:,:)     ! trigonometric factor cosine and sine at these angles
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_ang_1, n_ang_2, n_r                                        ! sizes of 3D real output array
        real(dp), allocatable :: cos_theta(:,:,:,:)                             ! cos(m theta) for all m
        real(dp), allocatable :: sin_theta(:,:,:,:)                             ! sin(m theta) for all m
        real(dp), allocatable :: cos_zeta(:,:,:,:)                              ! cos(n nfp theta) for all n
        real(dp), allocatable :: sin_zeta(:,:,:,:)                              ! sin(n nfp theta) for all n
        integer :: n, m                                                         ! counters
        
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
        allocate(cos_theta(mpol,n_ang_1,n_ang_2,n_r))
        allocate(sin_theta(mpol,n_ang_1,n_ang_2,n_r))
        allocate(cos_zeta(2*ntor+1,n_ang_1,n_ang_2,n_r))
        allocate(sin_zeta(2*ntor+1,n_ang_1,n_ang_2,n_r))
        
        do m = 0,mpol-1
            cos_theta(m+1,:,:,:) = cos(m*theta)
            sin_theta(m+1,:,:,:) = sin(m*theta)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! MAYBE IN CALC_TRIGON_FACTORS, MAYBE YOU HAVE TO USE - ZETA???????  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do n = -ntor,ntor
            cos_zeta(n+ntor+1,:,:,:) = cos(n*nfp*zeta)
            sin_zeta(n+ntor+1,:,:,:) = sin(n*nfp*zeta)
        end do
        
        ! initialize trigon_factors
        allocate(trigon_factors(mpol,2*ntor+1,n_ang_1,n_ang_2,n_r,2))
        trigon_factors = 0.0_dp
        
        ! calculate cos(m theta - n nfp zeta) = cos(m theta) cos(n nfp zeta) + 
        !   sin(m theta) sin(n nfp zeta) and
        ! sin(m theta - n nfp zeta) = sin(m theta) cos(n nfp zeta) -
        !   cos(m theta) sin(n nfp zeta)
        do n = -ntor,ntor
            do m = 0,mpol-1
                trigon_factors(m+1,n+ntor+1,:,:,:,1) = &
                    &cos_theta(m+1,:,:,:)*cos_zeta(n+ntor+1,:,:,:) + &
                    &sin_theta(m+1,:,:,:)*sin_zeta(n+ntor+1,:,:,:)
                trigon_factors(m+1,n+ntor+1,:,:,:,2) = &
                    &sin_theta(m+1,:,:,:)*cos_zeta(n+ntor+1,:,:,:) - &
                    &cos_theta(m+1,:,:,:)*sin_zeta(n+ntor+1,:,:,:)
            end do
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
        real(dp), intent(in) :: var_fourier_c(:,:,:)                            ! cos factor of variable in Fourier space
        real(dp), intent(in) :: var_fourier_s(:,:,:)                            ! sin factor of variable in Fourier space
        real(dp), intent(in) :: trigon_factors(:,:,:,:,:,:)                     ! trigonometric factor cosine and sine at these angles
        real(dp), intent(inout) :: var_real(:,:,:)                              ! variable in real space
        integer, intent(in), optional :: deriv(2)                               ! optional derivatives in angular coordinates
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_ang_1, n_r, n_ang_2                                        ! sizes of 3D real output array
        integer :: m,n                                                          ! counters for mode numbers
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: fac_cos(:), fac_sin(:)                         ! factor in front of cos and sin, after taking derivatives
        real(dp), allocatable :: fac_trigon_temp(:)                             ! temporary variable that holds fac_cos or fac_sin
        
        ! initialize ierr
        ierr = 0
        
        ! set n_ang_1 and n_r
        n_ang_1 = size(trigon_factors,3)
        n_ang_2 = size(trigon_factors,4)
        n_r = size(trigon_factors,5)
        
        ! tests
        if (size(trigon_factors,6).ne.2) then
            ierr = 1
            err_msg = 'trigon_factors needs to contain sines and cosines'
            CHCKERR(err_msg)
        end if
        if (size(trigon_factors,1).ne.mpol .or. &
            &size(trigon_factors,2).ne.(2*ntor+1)) then
            ierr = 1
            err_msg = 'trigon_factors needs to be defined for the right number &
                &of modes'
            CHCKERR(err_msg)
        end if
        if (size(var_fourier_c,3).ne.n_r .or. &
            &size(var_fourier_s,3).ne.n_r) then
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
        do m = 0,mpol-1
            do n = -ntor,ntor
                ! initialize factors in front of cos and sin
                fac_cos = var_fourier_c(m+1,n+ntor+1,:)
                fac_sin = var_fourier_s(m+1,n+ntor+1,:)
                
                ! angular derivatives
                if (present(deriv)) then 
                    ! apply possible poloidal derivatives
                    do id = 1,deriv(1)
                        fac_trigon_temp = - m * fac_cos
                        fac_cos = m * fac_sin
                        fac_sin = fac_trigon_temp
                    end do
                    ! apply possible toroidal derivatives
                    do id = 1,deriv(2)
                        fac_trigon_temp = n * nfp * fac_cos
                        fac_cos = - n * nfp * fac_sin
                        fac_sin = fac_trigon_temp
                    end do
                end if
                
                ! sum
                do kd = 1,n_r
                    var_real(:,:,kd) = var_real(:,:,kd) + &
                        &fac_cos(kd)*trigon_factors(m+1,n+ntor+1,:,:,kd,1) + &
                        &fac_sin(kd)*trigon_factors(m+1,n+ntor+1,:,:,kd,2)
                end do
            end do
        end do
    end function fourier2real
    
    ! Normalizes VMEC input
    ! Note  that  the normal  VMEC coordinate  runs from  0 to  1, whatever  the
    ! normalization.
    subroutine normalize_VMEC
        use num_vars, only: ltest
        use eq_vars, only: pres_0, psi_0, R_0, B_0
        
        ! scale the VMEC quantities
        pres_V = pres_V/pres_0
        flux_t_V = flux_t_V/psi_0
        Dflux_t_V = Dflux_t_V/psi_0
        R_V_c = R_V_c/R_0
        R_V_s = R_V_s/R_0
        Z_V_c = Z_V_c/R_0
        Z_V_s = Z_V_s/R_0
        L_V_c = L_V_c
        L_V_s = L_V_s
#if ldebug
        if (ltest) then
            B_V_sub_s = B_V_sub_s/(R_0*B_0)
            B_V_sub_c = B_V_sub_c/(R_0*B_0)
            B_V_c = B_V_c/B_0
            B_V_s = B_V_s/B_0
            jac_V_c = jac_V_c/(R_0**3)
            jac_V_s = jac_V_s/(R_0**3)
        end if
#endif
    end subroutine normalize_VMEC
end module VMEC
