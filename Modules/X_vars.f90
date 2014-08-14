!------------------------------------------------------------------------------!
!   Holds variables pertaining to the perturbation quantities
!------------------------------------------------------------------------------!
module X_vars
    use num_vars, only: dp
    use output_ops, only: lvl_ud, writo
    use str_ops, only: r2strt
    
    implicit none
    
    private
    public calc_rho, &
        &rho, n_x, min_r, max_r, m_X, n_r_X, U_X_0, U_X_1, DU_X_0, &
        &DU_X_1, sigma, extra1, extra2, extra3, PV0, PV1, PV2, KV0, KV1, KV2, &
        &X_vec, X_val, min_m_X, max_m_X

    ! global variables
    real(dp), allocatable :: rho(:,:)                                           ! density
    integer :: n_X                                                              ! toroidal mode number
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    integer, allocatable :: m_X(:)                                              ! vector of poloidal mode numbers
    integer :: n_r_X                                                            ! number of normal points for perturbations
    real(dp) :: min_r, max_r                                                    ! mimimum and maximum radius
    complex(dp), allocatable :: U_X_0(:,:,:), U_X_1(:,:,:)                      ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
    complex(dp), allocatable :: DU_X_0(:,:,:), DU_X_1(:,:,:)                    ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
    real(dp), allocatable :: sigma(:,:)                                         ! parallel current
    real(dp), allocatable :: extra1(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra2(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra3(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    complex(dp), allocatable :: PV0(:,:,:,:)                                    ! ~PV^0 coefficient
    complex(dp), allocatable :: PV1(:,:,:,:)                                    ! ~PV^1 coefficient
    complex(dp), allocatable :: PV2(:,:,:,:)                                    ! ~PV^2 coefficient
    complex(dp), allocatable :: KV0(:,:,:,:)                                    ! ~KV^0 coefficient
    complex(dp), allocatable :: KV1(:,:,:,:)                                    ! ~KV^1 coefficient
    complex(dp), allocatable :: KV2(:,:,:,:)                                    ! ~KV^2 coefficient
    complex(dp), allocatable :: X_vec(:,:,:)                                    ! Eigenvector solution
    complex(dp), allocatable :: X_val(:)                                        ! EIgenvalue solution
    
contains
    ! calculate rho from user input
    ! TO BE IMPLEMENTED. TEMPORARILY SET TO 1 EVERYWHERE !!!
    subroutine calc_rho(n_par,n_r)
        ! input variables
        integer, intent(in) :: n_par, n_r                                       ! dimensions of the size to be allocated for the rho matrix
        
        call writo('calc_rho NOT YET IMPLEMENTED!!!')
        
        ! allocate rho
        if (allocated(rho)) deallocate(rho)
        allocate(rho(n_par,n_r))
        
        ! TEMPORAL REPLACEMENT !!!
        rho = 1.0E-7_dp
        call writo('TEMPORALLY SETTING rho to '//trim(r2strt(rho(1,1)))//'!!!')
    end subroutine
end module
