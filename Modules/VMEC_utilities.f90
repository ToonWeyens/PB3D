!------------------------------------------------------------------------------!
!   Numerical utilities related to the output of VMEC                          !
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
    logical :: debug_calc_trigon_factors = .false.                              ! plot debug information for calc_trigon_factors
#endif

    interface fourier2real
        module procedure fourier2real_1, fourier2real_2
    end interface

contains
    ! Calculate the trigonometric cosine and  sine factors on a grid (1:mnmax_V)
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

#if ldebug
        if (debug_calc_trigon_factors) then
            call writo('Calculate trigonometric factors for grid of size ['//&
                &trim(i2str(n_ang_1))//','//trim(i2str(n_ang_2))//','//&
                &trim(i2str(n_r))//']')
        end if
#endif
        
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
#if ldebug
            if (debug_calc_trigon_factors) then
                call writo('Calculating for m = '//trim(i2str(m))//&
                    &' of [0..'//trim(i2str(mpol_V-1))//']')
            end if
#endif
            cos_theta(m,:,:,:) = cos(m*theta)
            sin_theta(m,:,:,:) = sin(m*theta)
        end do
        do n = -ntor_V,ntor_V
#if ldebug
            if (debug_calc_trigon_factors) then
                call writo('Calculating for n = '//trim(i2str(n))//&
                    &' of ['//trim(i2str(-ntor_V))//'..'//&
                    &trim(i2str(ntor_V))//']')
            end if
#endif
            cos_zeta(n,:,:,:) = cos(n*nfp_V*zeta)
            sin_zeta(n,:,:,:) = sin(n*nfp_V*zeta)
        end do
        
        ! initialize trigon_factors
        allocate(trigon_factors(mnmax_V,n_ang_1,n_ang_2,n_r,2),STAT=ierr)
        CHCKERR('Failed to allocate trigonometric factors')
        trigon_factors = 0.0_dp
        
        ! calculate cos(m theta - n zeta) =
        !   cos(m theta) cos(n zeta) + sin(m theta) sin(n zeta)
        ! and sin(m theta - n zeta) =
        !   sin(m theta) cos(n zeta) - cos(m theta) sin(n zeta)
        ! Note: need to scale the indices mn_V(:,2) by nfp_V.
        do id = 1,mnmax_V
#if ldebug
            if (debug_calc_trigon_factors) then
                call writo('Calculating for id = '//trim(i2str(id))//&
                        &' of [1..'//trim(i2str(mnmax_V))//']')
            end if
#endif
            trigon_factors(id,:,:,:,1) = &
                &cos_theta(mn_V(id,1),:,:,:)*&
                &cos_zeta(mn_V(id,2)/nfp_V,:,:,:) + &
                &sin_theta(mn_V(id,1),:,:,:)*&
                &sin_zeta(mn_V(id,2)/nfp_V,:,:,:)
            trigon_factors(id,:,:,:,2) = &
                &sin_theta(mn_V(id,1),:,:,:)*&
                &cos_zeta(mn_V(id,2)/nfp_V,:,:,:) - &
                &cos_theta(mn_V(id,1),:,:,:)*&
                &sin_zeta(mn_V(id,2)/nfp_V,:,:,:)
        end do
        
        ! deallocate variables
        deallocate(cos_theta,sin_theta)
        deallocate(cos_zeta,sin_zeta)
    end function calc_trigon_factors
    
    ! Inverse Fourier transformation, from VMEC. Also calculates the poloidal or
    ! toroidal  derivatives  in  VMEC  coords., as  indicated  by  the  variable
    ! deriv(2).
    ! (Normal derivative  is done on the variables in  Fourier space, and should
    ! be provided here in varf_i if needed).
    ! There are two variants:
    !   1: version using trigon_factors, which is  useful when the grid on which
    !   the trigonometric  factors are defined  is not regular and  ideally when
    !   they are reused multiple times.
    !   2: version  using theta and  zeta directly,  which is useful  for small,
    !   unique calculations.
    ! Both  these  versions  make  use  of  a  factor  that  represents  angular
    ! derivatives. For deriv = [j,k], this is:
    !   m^j (-n)^k (-1)^((j+k+1)/2)             for varf_c,
    !   m^j (-n)^k (-1)^((j+k)/2)               for varf_c,
    ! where  the  divisions  have  to  be  done  using  integers,  i.e.  without
    ! remainder. The  first two factors  are straightforward, and the  third one
    ! originates in  the change of  sign when deriving a  cosine, but not  for a
    ! sine.
    ! Finally,  depending on (j+k)  is even  or odd, the  correct cos or  sin is
    ! chosen.
    integer function fourier2real_1(varf_c,varf_s,trigon_factors,varr,sym,&
        &deriv) result(ierr)
        
        character(*), parameter :: rout_name = 'fourier2real_1'
        
        ! input / output
        real(dp), intent(in) :: varf_c(:,:)                                     ! cos factor of variable in Fourier space
        real(dp), intent(in) :: varf_s(:,:)                                     ! sin factor of variable in Fourier space
        real(dp), intent(in) :: trigon_factors(:,:,:,:,:)                       ! trigonometric factor cosine and sine at these angles
        real(dp), intent(inout) :: varr(:,:,:)                                  ! variable in real space
        logical, intent(in), optional :: sym(2)                                 ! whether to use varf_c (1) and / or varf_s (2)
        integer, intent(in), optional :: deriv(2)                               ! optional derivatives in angular coordinates
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: dims(3)                                                      ! dimensions of varr
        integer :: id, kd                                                       ! counters
        integer :: deriv_loc(2)                                                 ! local derivative
        logical :: sym_loc(2)                                                   ! local sym
        real(dp) :: deriv_fac                                                   ! factur due to derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! set local deriv and sym
        deriv_loc = [0,0]
        sym_loc = [.true.,.true.]
        if (present(deriv)) deriv_loc = deriv
        if (present(sym)) sym_loc = sym
        
        ! set dimensions
        dims = shape(varr)
        
        ! test
        if (.not.sym_loc(1) .and. .not.sym_loc(2)) then
            ierr = 1
            err_msg = 'Need at least the cosine or the sine factor'
            CHCKERR(err_msg)
        end if
        if (size(trigon_factors,1).ne.mnmax_V .and. &
            &size(trigon_factors,2).ne.dims(1) .and. &
            &size(trigon_factors,3).ne.dims(2) .and. &
            &size(trigon_factors,4).ne.dims(3) .and. &
            &size(trigon_factors,5).ne.2) then
            ierr = 1
            err_msg = 'Trigonometric factors are not set up'
            CHCKERR(err_msg)
        end if
        
        ! initialize
        varr = 0.0_dp
        
        ! sum over modes
        do id = 1,mnmax_V
            ! setup derivative factor for varf_c
            deriv_fac = mn_V(id,1)**deriv_loc(1)*(-mn_V(id,2))**deriv_loc(2)*&
                &(-1)**((sum(deriv_loc)+1)/2)
            
            ! add terms ~ varf_c
            if (sym_loc(1)) then
                if (mod(sum(deriv_loc),2).eq.0) then                            ! even number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_c(id,kd) * deriv_fac * &
                            &trigon_factors(id,:,:,kd,1)
                    end do
                else                                                            ! odd number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_c(id,kd) * deriv_fac * &
                            &trigon_factors(id,:,:,kd,2)
                    end do
                end if
            end if
            
            ! setup derivative factor for varf_s
            deriv_fac = mn_V(id,1)**deriv_loc(1)*(-mn_V(id,2))**deriv_loc(2)*&
                &(-1)**(sum(deriv_loc)/2)
            
            ! add terms ~ varf_s
            if (sym_loc(2)) then
                if (mod(sum(deriv_loc),2).eq.0) then                            ! even number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_s(id,kd) * deriv_fac * &
                            &trigon_factors(id,:,:,kd,2)
                    end do
                else                                                            ! odd number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_s(id,kd) * deriv_fac * &
                            &trigon_factors(id,:,:,kd,1)
                    end do
                end if
            end if
        end do
    end function fourier2real_1
    integer function fourier2real_2(varf_c,varf_s,theta,zeta,varr,sym,deriv) &
        &result(ierr)
        
        character(*), parameter :: rout_name = 'fourier2real_2'
        
        ! input / output
        real(dp), intent(in) :: varf_c(:,:)                                     ! cos factor of variable in Fourier space
        real(dp), intent(in) :: varf_s(:,:)                                     ! sin factor of variable in Fourier space
        real(dp), intent(in) :: theta(:,:,:)                                    ! theta
        real(dp), intent(in) :: zeta(:,:,:)                                     ! zeta
        real(dp), intent(inout) :: varr(:,:,:)                                  ! variable in real space
        logical, intent(in), optional :: sym(2)                                 ! whether to use varf_c (1) and / or varf_s (2)
        integer, intent(in), optional :: deriv(2)                               ! optional derivatives in angular coordinates
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: dims(3)                                                      ! dimensions of varr
        integer :: id, kd                                                       ! counters
        integer :: deriv_loc(2)                                                 ! local derivative
        logical :: sym_loc(2)                                                   ! local sym
        real(dp) :: deriv_fac                                                   ! factur due to derivatives
        real(dp), allocatable :: ang(:,:,:)                                     ! angle of trigonometric functions
        
        ! initialize ierr
        ierr = 0
        
        ! set local deriv and sym
        deriv_loc = [0,0]
        sym_loc = [.true.,.true.]
        if (present(deriv)) deriv_loc = deriv
        if (present(sym)) sym_loc = sym
        
        ! test
        if (.not.sym_loc(1) .and. .not.sym_loc(2)) then
            ierr = 1
            err_msg = 'Need at least the cosine or the sine factor'
            CHCKERR(err_msg)
        end if
        
        ! set dimensions
        dims = shape(varr)
        
        ! initialize angle
        allocate(ang(dims(1),dims(2),dims(3)))
        
        ! initialize output
        varr = 0.0_dp
        
        ! sum over modes
        do id = 1,mnmax_V
            ! set angle
            ang = mn_V(id,1)*theta - mn_V(id,2)*zeta
            
            ! setup derivative factor for varf_c
            deriv_fac = mn_V(id,1)**deriv_loc(1)*(-mn_V(id,2))**deriv_loc(2)*&
                &(-1)**((sum(deriv_loc)+1)/2)
            
            ! add terms ~ varf_c
            if (sym_loc(1)) then
                if (mod(sum(deriv_loc),2).eq.0) then                            ! even number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_c(id,kd) * deriv_fac * cos(ang(:,:,kd))
                    end do
                else                                                            ! odd number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_c(id,kd) * deriv_fac * sin(ang(:,:,kd))
                    end do
                end if
            end if
            
            ! setup derivative factor for varf_s
            deriv_fac = mn_V(id,1)**deriv_loc(1)*(-mn_V(id,2))**deriv_loc(2)*&
                &(-1)**(sum(deriv_loc)/2)
            
            ! add terms ~ varf_s
            if (sym_loc(2)) then
                if (mod(sum(deriv_loc),2).eq.0) then                            ! even number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_s(id,kd) * deriv_fac * sin(ang(:,:,kd))
                    end do
                else                                                            ! odd number of derivatives
                    do kd = 1,dims(3)
                        varr(:,:,kd) = varr(:,:,kd) + &
                            &varf_s(id,kd) * deriv_fac * cos(ang(:,:,kd))
                    end do
                end if
            end if
        end do
    end function fourier2real_2
end module VMEC_utilities
