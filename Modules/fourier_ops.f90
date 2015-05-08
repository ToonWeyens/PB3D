!------------------------------------------------------------------------------!
!   Operations that  have to do with  the fourier representation used  and its !
!   relations with the real space                                              !
!------------------------------------------------------------------------------!
module fourier_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln

    implicit none
    private
    public calc_trigon_factors, fourier2real

contains
    ! Calculate the trigoniometric cosine and  sine factors on a grid (0:mpol-1,
    ! -ntor:ntor) at given 3D arrays for the (VMEC) E(quilibrium) angles theta_E
    ! and zeta_E.
    integer function calc_trigon_factors(theta,zeta,trigon_factors) &
        &result(ierr)
        use VMEC, only: mpol, ntor, nfp
        
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
        use VMEC, only: mpol, ntor, nfp
        
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
end module fourier_ops
