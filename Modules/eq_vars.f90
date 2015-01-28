!------------------------------------------------------------------------------!
!   Variables that have to do with equilibrium quantities and the grid used in !
!   the calculations.                                                          !
!------------------------------------------------------------------------------!
module eq_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, pi, max_str_ln
    use messages, only: writo, lvl_ud, print_ar_2, print_ar_1
    use str_ops, only: i2str

    implicit none
    private
    public dealloc_eq, dealloc_eq_final, &
        &grp_n_r_eq, n_r_eq, q_saf_E, q_saf_E_full, &
        &flux_p_E, flux_t_E, R_E, Z_E, L_E, pres_E, q_saf_FD, &
        &flux_p_FD, flux_t_FD, pres_FD, grp_min_r_eq, grp_max_r_eq, rho, R_0, &
        &pres_0, B_0, psi_0, rho_0, T_0, rot_t_FD, rot_t_E, flux_p_E_full, &
        &flux_t_E_full, rot_t_E_full, max_flux, max_flux_eq, &
        &max_flux_F, max_flux_eq_F, theta_E, zeta_E, trigon_factors

    ! global variables
    ! Note: The indices in [derivatives] are:
    !   [r,theta,zeta]
    real(dp), allocatable :: theta_E(:,:), zeta_E(:,:)                          ! grid points (n_par, grp_n_r_eq, 2) in equilibrium coords.
    real(dp), allocatable :: trigon_factors(:,:,:,:,:)                          ! trigonometric factor cosine for the inverse fourier transf.
    real(dp), allocatable :: R_E(:,:,:,:,:)                                     ! R in Equilibrium coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: Z_E(:,:,:,:,:)                                     ! Z in Equilibrium coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: L_E(:,:,:,:,:)                                     ! L(ambda) in Equilibrium coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: rho(:)                                             ! density (in all coord. systems)
    real(dp), allocatable :: pres_E(:,:)                                        ! pressure, and norm. Deriv. in equilibrium coords.
    real(dp), allocatable :: q_saf_E(:,:)                                       ! safety factor in equilibrium coordinates
    real(dp), allocatable :: q_saf_E_full(:,:)                                  ! safety factor in full normal grid in equilibrium coordinates
    real(dp), allocatable :: rot_t_E(:,:)                                       ! rot. transform in equilibrium coordinates
    real(dp), allocatable :: rot_t_E_full(:,:)                                  ! rot. transform in full normal grid in equilibrium coordinates
    real(dp), allocatable, target :: flux_p_E(:,:), flux_t_E(:,:)               ! pol. and tor. flux and norm. Deriv. in equilibrium coords.
    real(dp), allocatable :: flux_p_E_full(:,:), flux_t_E_full(:,:)             ! pol. flux, tor. flux, and norm. Deriv. values and Derivs. in equilibrium coords.
    real(dp), allocatable :: pres_FD(:,:)                                       ! pressure, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable :: q_saf_FD(:,:)                                      ! safety factor, Deriv. in Flux coords.
    real(dp), allocatable :: rot_t_FD(:,:)                                      ! rot. transform, Deriv. in Flux coords.
    real(dp), allocatable, target :: flux_p_FD(:,:), flux_t_FD(:,:)             ! pol. and tor. flux, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp) :: max_flux, max_flux_eq                                           ! max. flux (pol. or tor.) in Equilibrium coordinates
    real(dp) :: max_flux_F, max_flux_eq_F                                       ! max. flux (pol. or tor.) in Flux coordinates
    real(dp) :: R_0, pres_0, rho_0                                              ! independent normalization constants for nondimensionalization
    real(dp) :: B_0, psi_0, T_0                                                 ! derived normalization constants for nondimensionalization
    integer :: n_r_eq                                                           ! total nr. of normal points in equilibrim grid
    integer :: grp_n_r_eq                                                       ! nr. of normal points in this process in alpha group
    integer :: grp_min_r_eq, grp_max_r_eq                                       ! min. and max. r range of this process in alpha group

contains
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! equilibrium phase
    integer function dealloc_eq() result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'dealloc_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! deallocate general variables
        deallocate(flux_p_E)
        deallocate(flux_t_E)
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                deallocate(R_E,Z_E,L_E)
                deallocate(trigon_factors)
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function dealloc_eq
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! calculation for a certain alpha
    subroutine dealloc_eq_final
        use num_vars, only: grp_rank, use_pol_flux_X
        
        ! deallocate general variables
        deallocate(pres_FD)
        deallocate(flux_p_FD,flux_t_FD)
        if (use_pol_flux_X) then
            deallocate(q_saf_FD)
        else
            deallocate(rot_t_FD)
        end if
        
        deallocate(pres_E)
        if (use_pol_flux_X) then
            deallocate(q_saf_E)
        else
            deallocate(rot_t_E)
        end if
        
        deallocate(rho)
        
        ! group masters
        if (grp_rank.eq.0) then
            deallocate(flux_p_E_full)
            deallocate(flux_t_E_full)
            deallocate(q_saf_E_full)
            deallocate(rot_t_E_full)
        end if
        deallocate(zeta_E,theta_E)
    end subroutine dealloc_eq_final
end module eq_vars
