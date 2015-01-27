!------------------------------------------------------------------------------!
!   Holds variables pertaining to the perturbation quantities
!------------------------------------------------------------------------------!
module X_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, max_str_ln, iu
    
    implicit none
    
    private
    public dealloc_X, dealloc_X_final, &
        &n_x, m_X, n_r_X, U_X_0, U_X_1, DU_X_0, DU_X_1, extra1, grp_r_X, &
        &extra2, extra3, PV0, PV1, PV2, KV0, KV1, KV2, X_vec, X_val, min_m_X, &
        &max_m_X, grp_n_r_X, size_X, PV_int, KV_int, grp_min_r_X, grp_max_r_X, &
        &min_n_X, max_n_X, exp_ang_par_F, mu0sigma
    
    ! global variables
    integer :: size_X                                                           ! size of n_X and m_X (nr. of modes)
    integer :: min_n_X                                                          ! lowest poloidal mode number m_X
    integer :: max_n_X                                                          ! highest poloidal mode number m_X
    integer, allocatable :: n_X(:)                                              ! vector of poloidal mode numbers
    integer :: min_m_X                                                          ! lowest poloidal mode number m_X
    integer :: max_m_X                                                          ! highest poloidal mode number m_X
    integer, allocatable :: m_X(:)                                              ! vector of poloidal mode numbers
    integer :: n_r_X                                                            ! number of normal points for perturbations
    integer :: grp_n_r_X                                                        ! number of normal points for perturbations on this process
    integer :: grp_min_r_X, grp_max_r_X                                         ! min. and max. r range of this process in alpha group
    real(dp), allocatable :: grp_r_X(:)                                         ! normal points in Flux coords., globally normalized to (min_r_X..max_r_X)
    complex(dp), allocatable :: U_X_0(:,:,:), U_X_1(:,:,:)                      ! U_m(X_m) = [ U_m^0 + U_m^1 i/n d/dx] (X_m)
    complex(dp), allocatable :: DU_X_0(:,:,:), DU_X_1(:,:,:)                    ! d(U_m(X_m))/dtheta = [ DU_m^0 + DU_m^1 i/n d/dx] (X_m)
    real(dp), allocatable :: mu0sigma(:,:)                                      ! parallel current
    real(dp), allocatable :: extra1(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra2(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    real(dp), allocatable :: extra3(:,:)                                        ! extra terms in PV_0 (see routine calc_extra)
    complex(dp), allocatable :: PV0(:,:,:,:)                                    ! ~PV^0 coefficient
    complex(dp), allocatable :: PV1(:,:,:,:)                                    ! ~PV^1 coefficient
    complex(dp), allocatable :: PV2(:,:,:,:)                                    ! ~PV^2 coefficient
    complex(dp), allocatable :: KV0(:,:,:,:)                                    ! ~KV^0 coefficient
    complex(dp), allocatable :: KV1(:,:,:,:)                                    ! ~KV^1 coefficient
    complex(dp), allocatable :: KV2(:,:,:,:)                                    ! ~KV^2 coefficient
    complex(dp), allocatable :: PV_int(:,:,:,:)                                 ! <~PV e^i(k-m)ang_par_F> coefficient
    complex(dp), allocatable :: KV_int(:,:,:,:)                                 ! <~KV e^i(k-m)ang_par_F> coefficient
    complex(dp), allocatable :: X_vec(:,:,:)                                    ! Eigenvector solution, with ghost region
    complex(dp), allocatable :: X_val(:)                                        ! Eigenvalue solution
    complex(dp), allocatable :: exp_ang_par_F(:,:,:,:)                          ! exp(i (k-m) ang_par_F)
    
contains
    ! deallocates  perturbation quantities that  are not used anymore  after the
    ! equilibrium phase
    subroutine dealloc_X
        deallocate(U_X_0,U_X_1,DU_X_0,DU_X_1)
        deallocate(PV0,PV1,PV2,KV0,KV1,KV2)
        deallocate(exp_ang_par_F)
        deallocate(mu0sigma,extra1,extra2,extra3)
    end subroutine dealloc_X
    
    ! deallocates  perturbation quantities that  are not used anymore  after the
    ! calculations for a certain alpha
    subroutine dealloc_X_final
        deallocate(n_X,m_X)
        deallocate(X_vec,X_val)
        deallocate(PV_int,KV_int)
        deallocate(grp_r_X)
    end subroutine dealloc_X_final
end module
