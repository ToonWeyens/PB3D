!------------------------------------------------------------------------------!
!> Numerical utilities related to the output of VMEC.
!------------------------------------------------------------------------------!
module VMEC_utilities
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: &
        &dp, max_str_ln, pi
    use VMEC_vars
    
    implicit none
    private
    public calc_trigon_factors, fourier2real
#if ldebug
    public debug_calc_trigon_factors
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_trigon_factors = .false.                              !< plot debug information for calc_trigon_factors() \ldebug
#endif

    ! interfaces
    
    !> \public Inverse Fourier transformation, from VMEC.
    !!
    !! Also calculates the poloidal or  toroidal derivatives in VMEC coords., as
    !! indicated by the variable \c deriv(2).
    !!
    !! (The normal derivative  is done on the variables in  Fourier space, and should
    !! be provided here in \c varf_i if needed).
    !!
    !! There are two variants:
    !!  -# version using trigon_factors, which is  useful when the grid on which
    !!  the trigonometric  factors are defined  is not regular and  ideally when
    !!  they are reused multiple times.
    !!  -# version using \f$\theta\f$ and  \f$\zeta\f$ directly, which is useful
    !!  for small, unique calculations.
    !!
    !! Both  these  versions  make  use  of  a  factor  that  represents  angular
    !! derivatives. For <tt>deriv = [j,k]</tt>, this is:
    !!  - \f$m^j (-n)^k (-1)^{\frac{j+k+1}{2}}\f$   for \c varf_c,
    !!  - \f$m^j (-n)^k (-1)^{\frac{j+k}{2}}\f$     for \c varf_s,
    !!
    !! where  the  divisions  have  to  be  done  using  integers,  i.e.  without
    !! remainder. The  first two factors  are straightforward, and the  third one
    !! originates in  the change of  sign when deriving a  cosine, but not  for a
    !! sine.
    !!
    !! Finally, depending  on whether  \f$(j+k)\f$ is even  or odd,  the correct
    !! \f$\cos\f$ or \f$\sin\f$ is chosen.
    !!
    !! \return ierr
    interface fourier2real
        !> \public
        module procedure fourier2real_1
        !> \public
        module procedure fourier2real_2
    end interface

contains
    !> \private version with trigonometric factors
    !!
    !! The trigonometric  factors are in  basis form (see  calc_trigon_factors()):
    !! the cosine and  sine of \f$m \theta\f$ and of  \f$n N_\text{fp} \zeta\f$,
    !! from  which   \f$\cos(m \theta  - n  \zeta)\f$  and  \f$\sin(m  \theta -  n
    !! \zeta)\f$ are formed  per mode. Any mode  list within the basis  can thus be
    !! evaluated: the main VMEC set \c mn_V by default, or the Nyquist set
    !! \c mn_V_nyq of the Jacobian, \f$|B|\f$ and \f$B_i\f$ through \c mn.
    integer function fourier2real_1(varf_c,varf_s,trigon_factors,varr,sym,&
        &deriv,mn) result(ierr)
        
        character(*), parameter :: rout_name = 'fourier2real_1'
        
        ! input / output
        real(dp), intent(in) :: varf_c(:,:)                                     !< \f$\cos\f$ factor of variable in Fourier space
        real(dp), intent(in) :: varf_s(:,:)                                     !< \f$\sin\f$ factor of variable in Fourier space
        real(dp), intent(in) :: trigon_factors(:,:,:,:,:)                       !< trigonometric basis at these angles (see calc_trigon_factors())
        real(dp), intent(inout) :: varr(:,:,:)                                  !< variable in real space
        logical, intent(in), optional :: sym(2)                                 !< whether to use \c varf_c (1) and / or \c varf_s (2)
        integer, intent(in), optional :: deriv(2)                               !< optional derivatives in angular coordinates
        integer, intent(in), optional :: mn(:,:)                                !< mode numbers (m, n) of the Fourier coefficients (default: \c mn_V)
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: dims(3)                                                      ! dimensions of varr
        integer :: id, kd                                                       ! counters
        integer :: n_mod                                                        ! number of modes
        integer, allocatable :: i_m(:), i_n(:)                                  ! indices of m and n of every mode in the trigonometric basis
        integer :: i_c, i_s                                                     ! which of cos (1) or sin (2) multiplies varf_c and varf_s
        logical :: need(2)                                                      ! whether cos (1) and sin (2) of (m theta - n zeta) are needed
        integer :: deriv_loc(2)                                                 ! local derivative
        logical :: sym_loc(2)                                                   ! local sym
        integer, allocatable :: mn_loc(:,:)                                     ! local mn
        real(dp), allocatable :: fac_c(:), fac_s(:)                             ! factors due to derivatives for varf_c and varf_s, per mode
        real(dp), allocatable :: cs(:,:,:)                                      ! cos and sin of (m theta - n zeta) at one normal point
        
        ! initialize ierr
        ierr = 0
        
        ! set local deriv and sym
        deriv_loc = [0,0]
        sym_loc = [.true.,.true.]
        if (present(deriv)) deriv_loc = deriv
        if (present(sym)) sym_loc = sym
        
        ! set local mode list
        n_mod = size(varf_c,1)
        allocate(mn_loc(n_mod,2))
        if (present(mn)) then
            if (size(mn,1).ne.n_mod .or. size(mn,2).ne.2) then
                ierr = 1
                err_msg = 'Mode list does not match the Fourier coefficients'
                CHCKERR(err_msg)
            end if
            mn_loc = mn
        else
            if (n_mod.ne.mnmax_V) then
                ierr = 1
                err_msg = 'Fourier coefficients are not on the main VMEC mode &
                    &set; pass their mode list through mn'
                CHCKERR(err_msg)
            end if
            mn_loc = mn_V
        end if
        
        ! set dimensions
        dims = shape(varr)
        
        ! tests
        if (.not.sym_loc(1) .and. .not.sym_loc(2)) then
            ierr = 1
            err_msg = 'Need at least the cosine or the sine factor'
            CHCKERR(err_msg)
        end if
        if (size(varf_s,1).ne.n_mod .or. size(varf_c,2).ne.dims(3) .or. &
            &size(varf_s,2).ne.dims(3)) then
            ierr = 1
            err_msg = 'Fourier coefficients do not match the real-space &
                &variable'
            CHCKERR(err_msg)
        end if
        if (size(trigon_factors,1).ne.dims(1) .or. &
            &size(trigon_factors,2).ne.dims(2) .or. &
            &size(trigon_factors,3).ne.dims(3) .or. &
            &size(trigon_factors,4).ne.mpol_tf+2*ntor_tf+1 .or. &
            &size(trigon_factors,5).ne.2) then
            ierr = 1
            err_msg = 'Trigonometric factors are not set up'
            CHCKERR(err_msg)
        end if
        
        ! an even number of angular derivatives keeps cos with varf_c and sin
        ! with varf_s, an odd number swaps them
        if (mod(sum(deriv_loc),2).eq.0) then
            i_c = 1
            i_s = 2
        else
            i_c = 2
            i_s = 1
        end if
        
        ! a coefficient set that is zero everywhere (e.g. the sine coefficients
        ! of a stellarator-symmetric equilibrium) costs nothing
        if (sym_loc(1)) sym_loc(1) = any(abs(varf_c).gt.0._dp)
        if (sym_loc(2)) sym_loc(2) = any(abs(varf_s).gt.0._dp)
        need(1) = (sym_loc(1) .and. i_c.eq.1) .or. (sym_loc(2) .and. i_s.eq.1)
        need(2) = (sym_loc(1) .and. i_c.eq.2) .or. (sym_loc(2) .and. i_s.eq.2)
        
        ! initialize
        varr = 0.0_dp
        if (.not.sym_loc(1) .and. .not.sym_loc(2)) return
        
        ! per mode: indices in the trigonometric basis (n in units of nfp_V)
        ! and the derivative factors
        allocate(i_m(n_mod),i_n(n_mod),fac_c(n_mod),fac_s(n_mod))
        do id = 1,n_mod
            if (mn_loc(id,1).lt.0 .or. mn_loc(id,1).gt.mpol_tf-1 .or. &
                &abs(mn_loc(id,2)).gt.ntor_tf*nfp_V) then
                ierr = 1
                err_msg = 'Mode ('//trim(i2str(mn_loc(id,1)))//','//&
                    &trim(i2str(mn_loc(id,2)))//') lies outside the &
                    &trigonometric basis'
                CHCKERR(err_msg)
            end if
            i_m(id) = mn_loc(id,1)+1
            i_n(id) = mpol_tf+ntor_tf+1+mn_loc(id,2)/nfp_V
            fac_c(id) = mn_loc(id,1)**deriv_loc(1)*&
                &(-mn_loc(id,2))**deriv_loc(2)*(-1)**((sum(deriv_loc)+1)/2)
            fac_s(id) = mn_loc(id,1)**deriv_loc(1)*&
                &(-mn_loc(id,2))**deriv_loc(2)*(-1)**(sum(deriv_loc)/2)
        end do
        
        ! sum over modes, one normal point at a time so that the basis of that
        ! point stays in cache while the modes stream over it
        allocate(cs(dims(1),dims(2),2))
        do kd = 1,dims(3)
            do id = 1,n_mod
                if (abs(varf_c(id,kd)).le.0._dp .and. abs(varf_s(id,kd)).le.0._dp) cycle
                
                ! cos(m theta - n zeta) =
                !   cos(m theta) cos(n zeta) + sin(m theta) sin(n zeta)
                if (need(1)) cs(:,:,1) = &
                    &trigon_factors(:,:,kd,i_m(id),1)*&
                    &trigon_factors(:,:,kd,i_n(id),1) + &
                    &trigon_factors(:,:,kd,i_m(id),2)*&
                    &trigon_factors(:,:,kd,i_n(id),2)
                ! sin(m theta - n zeta) =
                !   sin(m theta) cos(n zeta) - cos(m theta) sin(n zeta)
                if (need(2)) cs(:,:,2) = &
                    &trigon_factors(:,:,kd,i_m(id),2)*&
                    &trigon_factors(:,:,kd,i_n(id),1) - &
                    &trigon_factors(:,:,kd,i_m(id),1)*&
                    &trigon_factors(:,:,kd,i_n(id),2)
                
                ! add terms ~ varf_c and ~ varf_s
                if (sym_loc(1)) varr(:,:,kd) = varr(:,:,kd) + &
                    &varf_c(id,kd)*fac_c(id)*cs(:,:,i_c)
                if (sym_loc(2)) varr(:,:,kd) = varr(:,:,kd) + &
                    &varf_s(id,kd)*fac_s(id)*cs(:,:,i_s)
            end do
        end do
    end function fourier2real_1
    !> \private Angle version: the trigonometric basis is built here for the
    !! angles given (calc_trigon_factors()) and the table version does the
    !! rest, so both versions share one code path.
    integer function fourier2real_2(varf_c,varf_s,theta,zeta,varr,sym,deriv,&
        &mn) result(ierr)
        
        character(*), parameter :: rout_name = 'fourier2real_2'
        
        ! input / output
        real(dp), intent(in) :: varf_c(:,:)                                     !< \f$\cos\f$ factor of variable in Fourier space
        real(dp), intent(in) :: varf_s(:,:)                                     !< \f$\sin\f$ factor of variable in Fourier space
        real(dp), intent(in) :: theta(:,:,:)                                    !< \f$\theta_\text{E}\f$
        real(dp), intent(in) :: zeta(:,:,:)                                     !< \f$\zeta_\text{E}\f$
        real(dp), intent(inout) :: varr(:,:,:)                                  !< variable in real space
        logical, intent(in), optional :: sym(2)                                 !< whether to use \c varf_c (1) and / or \c varf_s (2)
        integer, intent(in), optional :: deriv(2)                               !< optional derivatives in angular coordinates
        integer, intent(in), optional :: mn(:,:)                                !< mode numbers (m, n) of the Fourier coefficients (default: \c mn_V)
        
        ! local variables
        real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      ! trigonometric basis at these angles
        
        ierr = calc_trigon_factors(theta,zeta,trigon_factors)
        CHCKERR('')
        ierr = fourier2real_1(varf_c,varf_s,trigon_factors,varr,sym,deriv,mn)
        CHCKERR('')
    end function fourier2real_2
    
    !> Calculate the trigonometric cosine and sine basis.
    !! 
    !! The basis holds, at given 3-D arrays for the (VMEC) E(quilibrium) angles
    !! \f$\theta_\text{E}\f$ and \f$\zeta_\text{E}\f$,
    !!  - \f$\cos(m \theta)\f$ and \f$\sin(m \theta)\f$ for \f$m = 0 \ldots
    !!  \text{mpol_tf}-1\f$ at index \f$m+1\f$,
    !!  - \f$\cos(n N_\text{fp} \zeta)\f$  and \f$\sin(n N_\text{fp} \zeta)\f$ for
    !!  \f$n = -\text{ntor_tf} \ldots \text{ntor_tf}\f$ at index \f$\text{mpol_tf}
    !!  + \text{ntor_tf} + 1 + n\f$,
    !!
    !! so that the dimensions of the output array are
    !!
    !! <tt>(1:n_ang_1,1:n_ang_2,1:n_r,1:mpol_tf+2*ntor_tf+1,1:2)</tt>
    !!
    !! where the last index selects cosine  (1) or sine (2) and \c n_r  is the
    !! total number of normal points; the  points come first so that one basis
    !! function is contiguous over the grid. fourier2real() forms \f$\cos(m \theta  - n
    !! \zeta)\f$ and \f$\sin(m \theta - n \zeta)\f$ from the basis per mode, so
    !! the same  basis serves the  main VMEC mode set  and the Nyquist  set of the
    !! Jacobian,  \f$|B|\f$ and \f$B_i\f$:  \c mpol_tf and \c ntor_tf cover both
    !! (see VMEC_ops.read_VMEC()). Previously the  factors were stored per mode of
    !! the main set only, which is why  the Nyquist modes beyond it were dropped
    !! (issue #4, step 9).
    !!
    !! \see See grid_vars.grid_type for a discussion on \c ang_1 and \c ang_2.
    !!
    !! \return ierr
    integer function calc_trigon_factors(theta,zeta,trigon_factors) &
        &result(ierr)
        
        character(*), parameter :: rout_name = 'calc_trigon_factors'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:,:)                                    !< poloidal angles in equilibrium coords.
        real(dp), intent(in) :: zeta(:,:,:)                                     !< toroidal angles in equilibrium coords.
        real(dp), intent(inout), allocatable :: trigon_factors(:,:,:,:,:)       !< trigonometric basis at these angles
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_ang_1, n_ang_2, n_r                                        ! sizes of 3D real output array
        integer :: m, n                                                         ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! set n_ang_1 and n_r
        n_ang_1 = size(theta,1)
        n_ang_2 = size(theta,2)
        n_r = size(theta,3)

#if ldebug
        if (debug_calc_trigon_factors) then
            call writo('Calculate trigonometric basis for grid of size ['//&
                &trim(i2str(n_ang_1))//','//trim(i2str(n_ang_2))//','//&
                &trim(i2str(n_r))//'], m = 0..'//trim(i2str(mpol_tf-1))//&
                &', n = '//trim(i2str(-ntor_tf))//'..'//trim(i2str(ntor_tf)))
        end if
#endif
        
        ! tests
        if (size(zeta,1).ne.n_ang_1 .or. size(zeta,2).ne.n_ang_2 .or. &
            &size(zeta,3).ne.n_r) then
            ierr = 1
            err_msg = 'theta and zeta need to have the same size'
            CHCKERR(err_msg)
        end if
        if (mpol_tf.lt.1 .or. ntor_tf.lt.0) then
            ierr = 1
            err_msg = 'The extent of the trigonometric basis is not set'
            CHCKERR(err_msg)
        end if
        
        ! initialize trigon_factors
        allocate(trigon_factors(n_ang_1,n_ang_2,n_r,mpol_tf+2*ntor_tf+1,2),&
            &STAT=ierr)
        CHCKERR('Failed to allocate trigonometric factors')
        
        ! poloidal part: cos(m theta), sin(m theta)
        do m = 0,mpol_tf-1
            trigon_factors(:,:,:,m+1,1) = cos(m*theta)
            trigon_factors(:,:,:,m+1,2) = sin(m*theta)
        end do
        
        ! toroidal part: cos(n nfp zeta), sin(n nfp zeta)
        do n = -ntor_tf,ntor_tf
            trigon_factors(:,:,:,mpol_tf+ntor_tf+1+n,1) = cos(n*nfp_V*zeta)
            trigon_factors(:,:,:,mpol_tf+ntor_tf+1+n,2) = sin(n*nfp_V*zeta)
        end do
    end function calc_trigon_factors
end module VMEC_utilities
