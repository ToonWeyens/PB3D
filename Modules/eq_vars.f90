!------------------------------------------------------------------------------!
!   Variables,  subroutines  and functions  that  have  to do  with            !
!   equilibrium quantities and the mesh used in the calculations               !
!   These are called in the subroutine calc_eq in the module eq_ops            !
!------------------------------------------------------------------------------!
module eq_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, pi, max_str_ln
    use message_ops, only: writo, lvl_ud, print_ar_2, print_ar_1
    use output_ops, only: print_GP_3D, draw_GP, draw_GP_animated, print_GP_2D
    use str_ops, only: r2strt, i2str, r2str

    implicit none
    private
    public calc_eqd_mesh, calc_mesh, calc_ang_B, calc_flux_q, prepare_RZL, &
        &check_and_limit_mesh, init_eq, calc_RZL, dealloc_eq, calc_XYZ_grid, &
        &dealloc_eq_final, normalize_eq_vars, calc_norm_const, &
        &theta_V, zeta_V, n_par, grp_n_r_eq, lam_H, min_par, max_par, n_r_eq, &
        &q_saf_V, q_saf_V_full, flux_p_V, flux_t_V, VMEC_R, VMEC_Z, VMEC_L, &
        &pres_V, q_saf_FD, flux_p_FD, flux_t_FD, pres_FD, grp_min_r_eq, &
        &grp_max_r_eq, R_0, pres_0, B_0, psi_0, rho_0, ang_par_F, &
        &rot_t_FD, rot_t_V, flux_p_V_full, flux_t_V_full, rot_t_V_full, &
        &max_flux, max_flux_eq, eq_use_pol_flux, q_saf_H, rot_t_H, theta_H, &
        &zeta_H, flux_p_H, flux_t_H

    ! global variables
    real(dp), allocatable :: trigon_factors(:,:,:,:,:)                          ! trigonometric factor cosine for the inverse fourier transf.
    real(dp), allocatable :: VMEC_R(:,:,:,:,:)                                  ! R in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: VMEC_Z(:,:,:,:,:)                                  ! Z in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: VMEC_L(:,:,:,:,:)                                  ! L(ambda) in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: lam_H(:,:,:)                                       ! lambda in (HM)
    real(dp), allocatable :: theta_V(:,:), zeta_V(:,:)                          ! grid points (n_par, grp_n_r_eq, 2) in VMEC coords
    real(dp), allocatable :: theta_H(:,:), zeta_H(:,:)                          ! grid points (n_par, grp_n_r_eq, 2) in HELENA coords
    real(dp), allocatable :: ang_par_F(:,:)                                     ! parallel angle in flux coordinates (either zeta_F or theta_F)
    real(dp), allocatable :: q_saf_V(:,:)                                       ! safety factor in VMEC coordinates
    real(dp), allocatable :: q_saf_H(:,:)                                       ! safety factor in HELENA coordinates
    real(dp), allocatable :: q_saf_V_full(:,:)                                  ! safety factor in full normal mesh in VMEC coordinates
    real(dp), allocatable :: q_saf_H_full(:,:)                                  ! safety factor in full normal mesh in HELENA coordinates
    real(dp), allocatable :: q_saf_FD(:,:)                                      ! safety factor, Deriv. in Flux coords.
    real(dp), allocatable :: rot_t_V(:,:)                                       ! rot. transform in VMEC coordinates
    real(dp), allocatable :: rot_t_H(:,:)                                       ! rot. transform in HELENA coordinates
    real(dp), allocatable :: rot_t_V_full(:,:)                                  ! rot. transform in full normal mesh in VMEC coordinates
    real(dp), allocatable :: rot_t_H_full(:,:)                                  ! rot. transform in full normal mesh in HELENA coordinates
    real(dp), allocatable :: rot_t_FD(:,:)                                      ! rot. transform, Deriv. in Flux coords.
    real(dp), allocatable :: flux_p_V(:,:), flux_t_V(:,:), pres_V(:,:)          ! pol. flux, tor. flux and pressure, and norm. Deriv. in VMEC coords.
    real(dp), allocatable :: flux_p_H(:,:), flux_t_H(:,:), pres_H(:,:)          ! pol. flux, tor. flux and pressure, and norm. Deriv. in HELENA coords.
    real(dp), allocatable :: pres_FD(:,:)                                       ! pressure, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable, target :: flux_p_FD(:,:), flux_t_FD(:,:)             ! pol. and tor. flux, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable :: flux_p_V_full(:,:), flux_t_V_full(:,:)             ! pol. flux, tor. flux, and norm. Deriv. values and Derivs. in VMEC coords.
    real(dp), allocatable :: flux_p_H_full(:,:), flux_t_H_full(:,:)             ! pol. flux, tor. flux, and norm. Deriv. values and Derivs. in HELENA coords.
    real(dp) :: max_flux, max_flux_eq                                           ! max. flux (pol. or tor.) (min.flux is trivially equal to 0)
    real(dp) :: min_par, max_par                                                ! min. and max. of parallel coordinate [pi]
    real(dp) :: R_0, pres_0, rho_0                                              ! independent normalization constants for nondimensionalization
    real(dp) :: B_0, psi_0                                                      ! derived normalization constants for nondimensionalization
    integer :: n_par, n_r_eq                                                    ! total nr. of parallel and normal points in equilibrim grid
    integer :: grp_n_r_eq                                                       ! nr. of normal points in this process in alpha group
    integer :: grp_min_r_eq, grp_max_r_eq                                       ! min. and max. r range of this process in alpha group
    logical :: eq_use_pol_flux                                                  ! .true. if equilibrium uses pol. flux and .false. if toroidal flux
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! initialize the equilibrium variables
    integer function init_eq() result(ierr)
        use num_vars, only: max_deriv, grp_rank, use_pol_flux, eq_style
        
        character(*), parameter :: rout_name = 'init_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize variables that are used for all equilibrium styles
        ! calculate grp_n_r_eq
        grp_n_r_eq = grp_max_r_eq - grp_min_r_eq + 1
        
        ! pres_FD
        allocate(pres_FD(grp_n_r_eq,0:max_deriv(1)))
        
        ! flux_p_FD
        allocate(flux_p_FD(grp_n_r_eq,0:max_deriv(1)))
        
        ! flux_t_FD
        allocate(flux_t_FD(grp_n_r_eq,0:max_deriv(1)))
        
        if (use_pol_flux) then
            ! q_saf_FD
            allocate(q_saf_FD(grp_n_r_eq,0:max_deriv(1)))
        else
            ! rot_t_FD
            allocate(rot_t_FD(grp_n_r_eq,0:max_deriv(1)))
        end if
        
        ! initialize variables that are  specifici to which equilibrium style is
        ! being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! R
                allocate(VMEC_R(n_par,grp_n_r_eq,0:max_deriv(1),0:max_deriv(2),&
                    &0:max_deriv(3)))
                
                ! Z
                allocate(VMEC_Z(n_par,grp_n_r_eq,0:max_deriv(1),0:max_deriv(2),&
                    &0:max_deriv(3)))
                
                ! lambda
                allocate(VMEC_L(n_par,grp_n_r_eq,0:max_deriv(1),0:max_deriv(2),&
                    &0:max_deriv(3)))
                
                ! pres_V
                allocate(pres_V(grp_n_r_eq,0:max_deriv(1)))
                
                ! flux_p_V
                allocate(flux_p_V(grp_n_r_eq,0:max_deriv(1)))
                
                ! flux_t_V
                allocate(flux_t_V(grp_n_r_eq,0:max_deriv(1)))
                
                if (use_pol_flux) then
                    ! q_saf_V
                    allocate(q_saf_V(grp_n_r_eq,0:max_deriv(1)))
                else
                    ! rot_t_V
                    allocate(rot_t_V(grp_n_r_eq,0:max_deriv(1)))
                end if
                
                ! full variables for group masters
                if (grp_rank.eq.0) then
                    ! flux_p_V_full
                    allocate(flux_p_V_full(n_r_eq,0:max_deriv(1)))
                    
                    ! flux_t_V_full
                    allocate(flux_t_V_full(n_r_eq,0:max_deriv(1)))
                    
                    ! q_saf_V_full
                    allocate(q_saf_V_full(n_r_eq,0:max_deriv(1)))
                    
                    ! rot_t_V_full
                    allocate(rot_t_V_full(n_r_eq,0:max_deriv(1)))
                end if
            case (2)                                                            ! HELENA
                ! pres_H
                allocate(pres_H(grp_n_r_eq,0:max_deriv(1)))
                
                ! flux_p_H
                allocate(flux_p_H(grp_n_r_eq,0:max_deriv(1)))
                
                ! flux_t_H
                allocate(flux_t_H(grp_n_r_eq,0:max_deriv(1)))
                
                if (use_pol_flux) then
                    ! q_saf_H
                    allocate(q_saf_H(grp_n_r_eq,0:max_deriv(1)))
                else
                    ! rot_t_H
                    allocate(rot_t_H(grp_n_r_eq,0:max_deriv(1)))
                end if
                
                ! full variables for group masters
                if (grp_rank.eq.0) then
                    ! flux_p_H_full
                    allocate(flux_p_H_full(n_r_eq,0:max_deriv(1)))
                    
                    ! flux_t_H_full
                    allocate(flux_t_H_full(n_r_eq,0:max_deriv(1)))
                    
                    ! q_saf_H_full
                    allocate(q_saf_H_full(n_r_eq,0:max_deriv(1)))
                    
                    ! rot_t_H_full
                    allocate(rot_t_H_full(n_r_eq,0:max_deriv(1)))
                end if
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function init_eq
    
    ! prepare the cosine  and sine factors that are used  in the inverse Fourier
    ! transformation of R, Z and L and derivatives
    integer function prepare_RZL() result(ierr)
        use fourier_ops, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'prepare_RZL'
        
        ! initialize ierr
        ierr = 0
        
        ierr = calc_trigon_factors(theta_V,zeta_V,trigon_factors)
        CHCKERR('')
    end function prepare_RZL
    
    ! calculate R, Z and Lambda and  derivatives in VMEC coordinates at the mesh
    ! points given by the variables VMEC_theta and VMEC_zeta and at every normal
    ! point. The derivatives  are indicated by the variable "deriv"  which has 3
    ! indices
    integer function calc_RZL_ind(deriv) result(ierr)
        use fourier_ops, only: fourier2real
        use VMEC_vars, only: R_c, R_s, Z_c, Z_s, L_c, L_s
        use utilities, only: check_deriv
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_RZL_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_RZL')
        CHCKERR('')
        
        ! calculate the variables R,Z and their angular derivative
        ierr = fourier2real(R_c(:,:,grp_min_r_eq:grp_max_r_eq,deriv(1)),&
            &R_s(:,:,grp_min_r_eq:grp_max_r_eq,deriv(1)),&
            &trigon_factors,VMEC_R(:,:,deriv(1),deriv(2),deriv(3)),&
            &[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(Z_c(:,:,grp_min_r_eq:grp_max_r_eq,deriv(1)),&
            &Z_s(:,:,grp_min_r_eq:grp_max_r_eq,deriv(1)),&
            &trigon_factors,VMEC_Z(:,:,deriv(1),deriv(2),deriv(3)),&
            &[deriv(2),deriv(3)])
        CHCKERR('')
        ierr = fourier2real(L_c(:,:,grp_min_r_eq:grp_max_r_eq,deriv(1)),&
            &L_s(:,:,grp_min_r_eq:grp_max_r_eq,deriv(1)),&
            &trigon_factors,VMEC_L(:,:,deriv(1),deriv(2),deriv(3)),&
            &[deriv(2),deriv(3)])
        CHCKERR('')
    end function calc_RZL_ind
    integer function calc_RZL_arr(deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_RZL_arr'
        
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            ierr = calc_RZL_ind(deriv(:,id))
            CHCKERR('')
        end do
    end function calc_RZL_arr

    ! calculates flux quantities  and normal derivatives in  the VMEC coordinate
    ! system
    integer function calc_flux_q() result(ierr)
        use num_vars, only: eq_style, max_deriv, grp_rank, use_pol_flux
        use utilities, only: norm_deriv, calc_int
        
        character(*), parameter :: rout_name = 'calc_flux_q_VMEC'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_flux_q_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_flux_q_HEL()
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! VMEC version
        integer function calc_flux_q_VMEC() result(ierr)
            use VMEC_vars, only: iotaf, phi, phi_r, presf
            
            character(*), parameter :: rout_name = 'calc_flux_q_VMEC'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_p_full(:)                            ! version of dflux_p/dp on full normal mesh (1..n_r_eq)
            real(dp), allocatable :: flux_p_int_full(:)                         ! version of integrated flux_p on full normal mesh (1..n_r_eq)
            
            ! initialize ierr
            ierr = 0
            
            ! pressure: copy from VMEC and derive
            pres_V(:,0) = presf(grp_min_r_eq:grp_max_r_eq)
            do kd = 1, max_deriv(1)
                ierr = norm_deriv(pres_V(:,0),pres_V(:,kd),n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! set up helper variables to calculate poloidal flux
            allocate(Dflux_p_full(n_r_eq),flux_p_int_full(n_r_eq))
            Dflux_p_full = iotaf*phi_r
            ierr = calc_int(Dflux_p_full,&
                &[(kd*1.0_dp/(n_r_eq-1.0_dp),kd=0,n_r_eq-1)],flux_p_int_full)
            CHCKERR('')
            
            ! poloidal flux: calculate using iotaf and phi, phi_r
            ! easier to use full normal mesh flux_p because of the integral
            flux_p_V(:,1) = Dflux_p_full(grp_min_r_eq:grp_max_r_eq)
            flux_p_V(:,0) = flux_p_int_full(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv(1)
                ierr = norm_deriv(flux_p_V(:,1),flux_p_V(:,kd),&
                    &n_r_eq-1._dp,kd-1,1)
                CHCKERR('')
            end do
                
            ! toroidal flux: copy from VMEC and derive
            flux_t_V(:,0) = phi(grp_min_r_eq:grp_max_r_eq)
            flux_t_V(:,1) = phi_r(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv(1)
                ierr = norm_deriv(flux_t_V(:,1),flux_t_V(:,kd),n_r_eq-1._dp,&
                    &kd-1,1)
                CHCKERR('')
            end do
            
            if (use_pol_flux) then
                ! safety factor
                q_saf_V(:,0) = 1.0_dp/iotaf(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(q_saf_V(:,0),q_saf_V(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_p_int_full(n_r_eq)
            else
                ! rot. transform
                rot_t_V(:,0) = iotaf(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(rot_t_V(:,0),rot_t_V(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = phi(n_r_eq)
            end if
            
            ! max_flux_eq
            if (eq_use_pol_flux) then
                max_flux_eq = flux_p_int_full(n_r_eq)
            else
                max_flux_eq = phi(n_r_eq)
            end if
                
            ! the global master needs flux_p_V_full, flux_t_V_full, q_saf_V_full
            ! and rot_t_V_full  for the resonance plot  and checking of m  and n
            ! and  the  group masters  need  q_saf_V  on  full normal  mesh  for
            ! plot_X_vec
            if (grp_rank.eq.0) then
                ! flux_t_V_full
                flux_t_V_full(:,0) = phi
                flux_t_V_full(:,1) = phi_r
                do kd = 2,max_deriv(1)
                    ierr = norm_deriv(flux_t_V_full(:,1),flux_t_V_full(:,kd),&
                        &n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! flux_p_V_full
                flux_p_V_full(:,0) = flux_p_int_full
                flux_p_V_full(:,1) = Dflux_p_full
                do kd = 2,max_deriv(1)
                    ierr = norm_deriv(flux_p_V_full(:,1),&
                        &flux_p_V_full(:,kd),n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! q_saf_V_full
                q_saf_V_full(:,0) = 1.0_dp/iotaf
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(q_saf_V_full(:,0),&
                        &q_saf_V_full(:,kd),n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
                
                ! rot_t_V_full
                rot_t_V_full(:,0) = iotaf
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(rot_t_V_full(:,0),&
                        &rot_t_V_full(:,kd),n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
            end if
            
            ! deallocate helper variables
            if (use_pol_flux .or. grp_rank.eq.0) then
                deallocate(Dflux_p_full,flux_p_int_full)
            end if
        end function calc_flux_q_VMEC
        
        ! HELENA version
        integer function calc_flux_q_HEL() result(ierr)
            use HEL_vars, only: qs, flux_H, p0
            
            character(*), parameter :: rout_name = 'calc_flux_q_HEL'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_t_full(:)                            ! version of dflux_t/dp on full normal mesh (1..n_r_eq)
            real(dp), allocatable :: flux_t_int_full(:)                         ! version of integrated flux_t on full normal mesh (1..n_r_eq)
            real(dp), allocatable :: flux_H_r(:)                                ! normal derivative of flux_H
            
            ! initialize ierr
            ierr = 0
            
            ! pressure: copy from HELENA and derive
            pres_H(:,0) = p0(grp_min_r_eq:grp_max_r_eq)
            do kd = 1, max_deriv(1)
                ierr = norm_deriv(pres_H(:,0),pres_H(:,kd),n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! set up helper variables to calculate toroidal flux
            ! calculate normal derivative of flux_H
            allocate(flux_H_r(n_r_eq))
            ierr = norm_deriv(flux_H,flux_H_r,n_r_eq-1._dp,1,1)
            CHCKERR('')
            allocate(Dflux_t_full(n_r_eq),flux_t_int_full(n_r_eq))
            Dflux_t_full = qs*flux_H_r
            ierr = calc_int(Dflux_t_full,&
                &[(kd*1.0_dp/(n_r_eq-1.0_dp),kd=0,n_r_eq-1)],flux_t_int_full)
            CHCKERR('')
            
            ! toroidal flux: calculate using qs and flux_H, flux_H_r
            ! easier to use full normal mesh flux_t because of the integral
            flux_t_H(:,1) = Dflux_t_full(grp_min_r_eq:grp_max_r_eq)
            flux_t_H(:,0) = flux_t_int_full(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv(1)
                ierr = norm_deriv(flux_t_H(:,1),flux_t_H(:,kd),&
                    &n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
                
            ! poloidal flux: copy from HELENA and derive
            flux_p_H(:,0) = flux_H(grp_min_r_eq:grp_max_r_eq)
            do kd = 1,max_deriv(1)
                ierr = norm_deriv(flux_p_H(:,1),flux_p_H(:,kd),n_r_eq-1._dp,&
                    &kd,1)
                CHCKERR('')
            end do
            
            if (use_pol_flux) then
                ! safety factor
                q_saf_H(:,0) = qs(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(q_saf_H(:,0),q_saf_H(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_H(n_r_eq)
            else
                ! rot. transform
                rot_t_H(:,0) = 1.0_dp/qs(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(rot_t_H(:,0),rot_t_H(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_t_int_full(n_r_eq)
            end if
            
            ! max_flux_eq
            max_flux_eq = flux_H(n_r_eq)
            
            ! the global master needs flux_p_H_full, flux_t_H_full, q_saf_H_full
            ! and rot_t_H_full  for the resonance plot  and checking of m  and n
            ! and  the  group masters  need  q_saf_H  on  full normal  mesh  for
            ! plot_X_vec
            if (grp_rank.eq.0) then
                ! flux_t_H_full
                flux_t_H_full(:,0) = flux_t_int_full
                flux_t_H_full(:,1) = Dflux_t_full
                do kd = 2,max_deriv(1)
                    ierr = norm_deriv(flux_t_H_full(:,1),&
                        &flux_t_H_full(:,kd),n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! flux_p_H_full
                flux_p_H_full(:,0) = flux_H
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(flux_p_H_full(:,1),flux_p_H_full(:,kd),&
                        &n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
                
                ! q_saf_H_full
                q_saf_H_full(:,0) = qs
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(q_saf_H_full(:,0),&
                        &q_saf_H_full(:,kd),n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
                
                ! rot_t_H_full
                rot_t_H_full(:,0) = 1.0_dp/qs
                do kd = 1,max_deriv(1)
                    ierr = norm_deriv(rot_t_H_full(:,0),&
                        &rot_t_H_full(:,kd),n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
            end if
            
            ! deallocate helper variables
            if (use_pol_flux .or. grp_rank.eq.0) then
                deallocate(Dflux_t_full,flux_t_int_full)
            end if
        end function calc_flux_q_HEL
    end function calc_flux_q

    ! Calculate the angular mesh.
    ! The variable  use_pol_flux determines  whether theta (.true.)  or zeta
    ! (.false.) is used as the parallel variable.
    integer function calc_mesh(alpha) result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'calc_mesh'
        
        ! input / output
        real(dp), intent(in) :: alpha
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize some variables
        allocate(ang_par_F(n_par,n_r_eq)); ang_par_F = 0.0_dp
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_mesh_VMEC(alpha)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_mesh_HEL(alpha)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! VMEC  Version.
        ! For  the normal  mesh type  (following the  magnetic field),  theta_V,
        ! zeta_V and ang_par_F are defined on  the entire normal mesh. Later, in
        ! check_and_limit_mesh, the  normal extent is  limited if the  tests are
        ! positive. 
        ! Aside from the  normal computational mesh, there are  other mesh types
        ! and  the global  variable calc_mesh_style  can optionally  force other
        ! types of meshes, which is useful for tests concerning VMEC.
        integer function calc_mesh_VMEC(alpha) result(ierr)
            use num_vars, only: calc_mesh_style, use_pol_flux
            use VMEC_vars, only: iotaf
            
            character(*), parameter :: rout_name = 'calc_mesh_VMEC'
            
            ! input / output
            real(dp), intent(in) :: alpha
            
            ! local variables
            integer :: kd, id                                                   ! counters
            real(dp), allocatable :: var_in(:)                                  ! input of calc_ang_B
            real(dp), allocatable :: var_out(:)                                 ! output of calc_ang_B
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! initialize variables
            allocate(zeta_V(n_par,n_r_eq)); zeta_V = 0.0_dp
            allocate(theta_V(n_par,n_r_eq)); theta_V = 0.0_dp
            
            allocate(var_in(n_r_eq))
            allocate(var_out(n_r_eq))
            
            ! calculate the mesh
            select case (calc_mesh_style)
                ! grid along the magnetic field lines
                case (0)
                    ! set up  parallel angle in flux  coordinates on equidistant
                    ! mesh
                    ierr = calc_eqd_mesh(ang_par_F(:,1),n_par,min_par,max_par)
                    CHCKERR('')
                    do kd = 2,n_r_eq
                        ang_par_F(:,kd) = ang_par_F(:,1)
                    end do
                    ! calculate zeta_V from alpha and ang_par_F
                    if (use_pol_flux) then                                      ! theta_F is parallel coordinate
                        do id = 1,n_par
                            zeta_V(id,:) = ang_par_F(id,:)/iotaf
                        end do
                        zeta_V = zeta_V - alpha
                    else                                                        ! zeta_F is parallel coordinate
                        zeta_V = - ang_par_F
                    end if
                    ! calculate theta_V from alpha and zeta_V
                    do id = 1,n_par
                        var_in = zeta_V(id,:)
                        ierr = calc_ang_B(var_out,.true.,alpha,var_in)
                        theta_V(id,:) = var_out
                        CHCKERR('')
                    end do
                ! grid with constant theta_V, equidistant zeta_V
                case (1)
                    theta_V = alpha
                    do kd = 1,n_r_eq
                        ierr = calc_eqd_mesh(zeta_V(:,kd),n_par,min_par,max_par)
                        CHCKERR('')
                end do
                    call writo('Defining grid with theta_V = '//&
                        &trim(r2strt(theta_V(1,1)))//&
                        &' and zeta_V equidistant from '//&
                        &trim(r2strt(min_par))//' pi to '//&
                        &trim(r2strt(max_par))//' pi')
                ! grid with constant zeta_V, equidistant theta_V
                case (2)
                    zeta_V = alpha
                    do kd = 1,n_r_eq
                        ierr = calc_eqd_mesh(theta_V(:,kd),n_par,min_par,&
                            &max_par)
                        CHCKERR('')
                    end do
                    call writo('Defining grid with zeta_V = '//&
                        &trim(r2strt(zeta_V(1,1)))//&
                        &' and theta_V equidistant from '//&
                        &trim(r2strt(min_par))//' pi to '//&
                        &trim(r2strt(max_par))//' pi')
                ! grid style error
                case default
                    err_msg = 'No style is associated with '&
                        &//trim(i2str(calc_mesh_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end function calc_mesh_VMEC
        
        ! HELENA  Version.
        ! Only the normal mesh type (following the magnetic field) is used, with
        ! theta_H,  zeta_H and  ang_par_F  defined on  the  entire normal  mesh.
        ! Later, in  check_and_limit_mesh, the normal  extent is limited  if the
        ! tests are positive.
        integer function calc_mesh_HEL(alpha) result(ierr)
            use HEL_vars, only: qs
            
            character(*), parameter :: rout_name = 'calc_mesh_HEL'
            
            ! input / output
            real(dp), intent(in) :: alpha
            
            ! local variables
            integer :: kd, id                                                   ! counters
            
            ! initialize ierr
            ierr = 0
            
            ! initialize variables
            allocate(zeta_H(n_par,n_r_eq)); zeta_H = 0.0_dp
            allocate(theta_H(n_par,n_r_eq)); theta_H = 0.0_dp
            
            ! set  up parallel  angle in  flux coordinates  on equidistant  mesh
            ! Note: this includes  chi, or half the array chi,  except for maybe
            ! the last point if HELENA is top-bottom symmetric (see read_HEL)
            ierr = calc_eqd_mesh(ang_par_F(:,1),n_par,min_par,max_par)
            CHCKERR('')
            do kd = 2,n_r_eq
                ang_par_F(:,kd) = ang_par_F(:,1)
            end do
            
            ! calculate theta_H and zeta_H from alpha and ang_par_F
            theta_H = ang_par_F                                                 ! theta_F is parallel coordinate
            do id = 1,n_par
                zeta_H(id,:) = ang_par_F(id,:)*qs
            end do
            zeta_H = zeta_H + alpha
        end function calc_mesh_HEL
    end function calc_mesh
    
    ! check  whether   the  straight  field   line  mesh  has   been  calculated
    ! satisfactorily,  by calculating  again the  theta_V's from  the calculated
    ! zeta_V's
    integer function check_and_limit_mesh(alpha) result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'check_and_limit_mesh'
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = check_and_limit_mesh_VMEC(alpha)
                CHCKERR('')
            case (2)                                                            ! HELENA
                call check_and_limit_mesh_HEL
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! VMEC version
        integer function check_and_limit_mesh_VMEC(alpha) result(ierr)
            use num_vars, only: tol_NR, calc_mesh_style, use_pol_flux
            
            character(*), parameter :: rout_name = 'check_and_limit_mesh_VMEC'
            
            ! input / output
            real(dp) :: alpha
            
            ! local variables
            real(dp) :: var_calc(n_par,n_r_eq)
            real(dp) :: var_diff(n_par,n_r_eq)
            real(dp) :: var_in(n_r_eq)
            integer :: id
            real(dp), allocatable :: work(:)                                    ! work array
            character(len=5) :: par_c                                           ! name of parallel coord.
            character(len=8) :: nor_c                                           ! name of normal coord.
            
            ! initialize ierr
            ierr = 0
            
            if (calc_mesh_style.eq.0) then                                      ! only do the check if the magnetic field lines are followed
                ! allocate work variable
                allocate(work(n_r_eq))
                
                do id = 1,n_par
                    var_in = theta_V(id,:)
                    ierr = calc_ang_B(work,.false.,alpha,var_in)
                    var_calc(id,:) =  work
                    CHCKERR('')
                end do
                var_diff = zeta_V - var_calc
                
                ! check the difference
                if (maxval(abs(var_diff)).gt.tol_NR*100) then                   ! difference too large
                    err_msg = 'Calculating again the toroidal grid points from &
                        &the poloidal grid points, along the magnetic fields, &
                        &does not result in the same values as initally &
                        &given... The maximum error is equal to '//&
                        &trim(r2strt(100*maxval(abs(var_diff))))//'%'
                    ierr = 1
                    CHCKERR(err_msg)
                else                                                            ! difference within tolerance
                    ! limit the normal size of theta_V and zeta_V
                    call limit_normal_size(theta_V)
                    call limit_normal_size(zeta_V)
                    call limit_normal_size(ang_par_F)
                end if
                
                ! deallocate work variable
                deallocate(work)
                
                ! user message
                if (use_pol_flux) then
                    par_c = 'theta'
                else
                    par_c = 'zeta'
                end if
                if (eq_use_pol_flux) then
                    nor_c = 'poloidal'
                else
                    nor_c = 'toroidal'
                end if
                call writo('angular grid set up with coordinate '//&
                    &trim(par_c)//' oriented along the magnetic field')
                call writo('normal grid set with '//trim(nor_c)//&
                    &' flux as normal coordinate')
                if (use_pol_flux.neqv.eq_use_pol_flux) write(*,*) &
                    &'!!!!!!!!!!!!!!! ADD TEST TO CHECK WHETHER WE CAN INDEED &
                    &USE DIFFERENT FLUX THAN VMEC!!!!!!!!!!!'
            else
                ! limit the normal size of theta_V and zeta_V
                call limit_normal_size(theta_V)
                call limit_normal_size(zeta_V)
                call limit_normal_size(ang_par_F)
            end if
        end function check_and_limit_mesh_VMEC
        
        ! HELENA version: the coordinates are already magnetic
        subroutine check_and_limit_mesh_HEL()
            ! limit the normal size of theta_H and zeta_H
            call limit_normal_size(theta_H)
            call limit_normal_size(zeta_H)
            call limit_normal_size(ang_par_F)
            
            call writo('angular grid set up with coordinate theta oriented &
                &along the magnetic field')
            call writo('normal grid set with poloidal flux as normal &
                &coordinate')
        end subroutine check_and_limit_mesh_HEL
        
        ! limits the normal size of theta_V and zeta_V
        subroutine limit_normal_size(angle)
            ! input / output
            real(dp), intent(inout), allocatable :: angle(:,:)
            
            ! local variables
            real(dp), allocatable :: work2(:,:)                                 ! work array
            
            ! allocate work array
            allocate(work2(n_par,n_r_eq))
            
            work2 = angle
            deallocate(angle); allocate(angle(n_par,grp_n_r_eq))
            angle = work2(:,grp_min_r_eq:grp_max_r_eq)
        end subroutine limit_normal_size
    end function check_and_limit_mesh

    ! calculate mesh of equidistant points
    ! Note: input is given in units of pi
    integer function calc_eqd_mesh(eqd_mesh,n_ang, min_ang, max_ang) result(ierr)
        character(*), parameter :: rout_name = 'eqd_mesh'
        
        ! input and output
        real(dp), intent(inout) :: eqd_mesh(:)                                  ! output
        real(dp), intent(in) :: min_ang, max_ang                                ! min. and max. of angles [pi]
        integer, intent(in) :: n_ang                                            ! nr. of points
        
        ! local variables
        integer :: id
        real(dp) :: delta_ang
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test some values
        if (n_ang.lt.1) then
            err_msg = 'The angular array has to have a length of at &
                &least 1'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (min_ang.gt.max_ang) then
            err_msg = 'The minimum angle has to be smaller than the &
                &maximum angle'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initialize output vector
        eqd_mesh = 0.0_dp
        
        ! There are (n_ang-1)  pieces in the total interval but  the last one is
        ! not included
        
        delta_ang = (max_ang-min_ang)/(n_ang) * pi
        
        eqd_mesh(1) = min_ang*pi
        do id = 2,n_ang
            eqd_mesh(id) = eqd_mesh(id-1) + delta_ang
        end do
    end function calc_eqd_mesh
    
    ! Calculates the  VMEC angle  theta_V or  zeta_V as a  function of  the VMEC
    ! angle zeta_V or theta_V following a particular magnetic field line alpha.
    ! This is done  using a Newton-Rhapson scheme that calculates  the zero's of
    ! either the function
    !   f = - zeta_V + q_V (theta_V + lambda(theta_V,zeta_V)) - alpha_F
    ! if the poloidal flux is used as normal coordinate, or
    !   f = iota_V zeta_V - (theta_V + lambda(theta_V,zeta_V)) - alpha_T
    ! if the toroidal flux is used as normal coordinate
    ! The  logical find_theta determines  whether theta_V (.true.) is  sought or
    ! zeta_V (.false.)
    integer function calc_ang_B(ang_B,find_theta,alpha_in,input_ang,guess) &
        &result(ierr)
        use num_vars, only: tol_NR, use_pol_flux
        use VMEC_vars, only: L_c, L_s, iotaf
        use fourier_ops, only: calc_trigon_factors, fourier2real
        use utilities, only: calc_zero_NR
        
        character(*), parameter :: rout_name = 'ang_B'
        
        ! input / output
        real(dp) :: ang_B(n_r_eq)                                               ! theta_F(zeta_F)/zeta_F(theta_F)
        logical, intent(in) :: find_theta                                       ! find theta_F or zeta_F
        real(dp), intent(in) :: alpha_in                                        ! alpha
        real(dp), intent(in) :: input_ang(n_r_eq)                               ! the input angle zeta_F/theta_F
        real(dp), intent(in), optional :: guess(n_r_eq)                         ! optional input (guess) for theta_F/zeta_F
        
        ! local variables (also used in child functions)
        integer :: kd                                                           ! counters
        real(dp) :: lam(1,1)                                                    ! lambda (needs to be 2D matrix)
        real(dp) :: dlam(1,1)                                                   ! angular derivative of lambda (needs to be 2D matrix)
        real(dp) :: theta_loc(1,1),zeta_loc(1,1)                                ! theta_V and zeta_V (needs to be 2D matrix)
        
        ! local variables (not to be used in child functions)
        real(dp) :: ang_NR                                                      ! temporary solution for a given r, iteration
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! for all normal points
        ang_B = 0.0_dp
        norm: do kd = 1, n_r_eq
            ! if first guess for theta_F is given
            if (present(guess)) then
                ang_NR = guess(kd)
            else if (kd.eq.1) then                                              ! take pi, because it is in the middle of 0...2pi
                ang_NR = pi
            else                                                                ! take solution for previous flux surface
                ang_NR = ang_B(kd-1)
            end if
            
            ! Newton-Rhapson loop for current normal point
            ierr = calc_zero_NR(ang_B(kd),fun_ang_B,dfun_ang_B,ang_NR)
            CHCKERR('')
            
            ! do a check  whether the result is indeed alpha,  making use of the
            if (abs(fun_ang_B(ang_B(kd))).gt.tol_NR*100) then
                err_msg = 'In theta_B, calculating alpha as a check, using &
                    &the theta_V that is the solution of alpha = '&
                    &//trim(r2strt(alpha_in))//', yields a calcualted alpha &
                    &that deviates '//&
                    &trim(r2strt(100*abs(fun_ang_B(ang_B(kd)))))&
                    &//'% from the original alpha_in'
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do norm
    contains
        ! function that returns  f = alpha -  alpha_0. It uses kd  from the main
        ! loop in the  parent function as the normal position  where to evaluate
        ! the quantitites,  lam and dlam to  contain the variable lambda  or its
        ! derivative, theta_loc and zeta_loc as  the local angular variables and
        ! input_ang(kd) as the angle for which  the magnetic match is sought, as
        ! well as iotaf, alpha_in, L_s and L_c
        function fun_ang_B(ang)
            character(*), parameter :: rout_name = 'fun_ang_B'
            
            ! input / output
            real(dp) :: fun_ang_B
            real(dp), intent(in) :: ang
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:)              ! trigonometric factor cosine for the inverse fourier transf.
            
            ! initialize fun_ang_B
            fun_ang_B = 0.0_dp
            
            ! transform lambda from Fourier space to real space
            ! set up local theta_V and zeta_V
            if (find_theta) then                                                ! looking for theta_F
                theta_loc = ang
                zeta_loc = input_ang(kd)
            else                                                                ! looking for zeta_F
                theta_loc = input_ang(kd)
                zeta_loc = ang
            end if
            ! calculate the (co)sines
            ierr = calc_trigon_factors(theta_loc,zeta_loc,trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_c(:,:,kd:kd,0),L_s(:,:,kd:kd,0),&
                &trigon_factors_loc,lam)
            CHCKERR('')
            
            ! calculate the output function
            if (use_pol_flux) then                                              ! poloidal flux as normal coordinate
                fun_ang_B = - zeta_loc(1,1) + &
                    &(theta_loc(1,1)+lam(1,1))/iotaf(kd) - alpha_in
            else                                                                ! toroidal flux as normal coordinate
                fun_ang_B = iotaf(kd)*zeta_loc(1,1) - &
                    &(theta_loc(1,1)+lam(1,1)) - alpha_in
            end if
        end function fun_ang_B
        
        ! function that returns  df/d(ang) = d(alpha -  alpha_0)/d(ang). It uses
        ! kd from  the main loop in  the parent function as  the normal position
        ! where  to  evaluate the  quantitites,  lam  and  dlam to  contain  the
        ! variable lambda or its derivative, theta_loc and zeta_loc as the local
        ! angular  variables  and  input_ang(kd)  as the  angle  for  which  the
        ! magnetic match is sought, as well as iotaf, alpha_in, L_s and L_c
        function dfun_ang_B(ang)
            character(*), parameter :: rout_name = 'dfun_ang_B'
            
            ! input / output
            real(dp) :: dfun_ang_B
            real(dp), intent(in) :: ang
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:)              ! trigonometric factor cosine for the inverse fourier transf.
            integer :: derivs_loc(2)                                            ! derivatives in theta_V or zeta_V
            
            ! initialize dfun_ang_B
            dfun_ang_B = 0.0_dp
            
            ! transform lambda from Fourier space to real space
            ! set up local theta_V and zeta_V
            if (find_theta) then                                                ! looking for theta_F
                theta_loc = ang
                zeta_loc = input_ang(kd)
            else                                                                ! looking for zeta_F
                theta_loc = input_ang(kd)
                zeta_loc = ang
            end if
            ! calculate the (co)sines
            ierr = calc_trigon_factors(theta_loc,zeta_loc,trigon_factors_loc)
            CHCKERR('')
            
            ! calculate angular derivatives of lambda
            if (find_theta) then                                                ! looking for theta_F
                derivs_loc = [1,0]
            else                                                                ! looking for zeta_F
                derivs_loc = [0,1]
            end if
            
            ierr = fourier2real(L_c(:,:,kd:kd,0),L_s(:,:,kd:kd,0),&
                &trigon_factors_loc,dlam,derivs_loc)
            CHCKERR('')
            
            ! calculate the output function
            if (find_theta) then                                                ! looking for theta_F
                !dfun_ang_B = -q(kd)
            else                                                                ! looking for zeta_F
                !dfun_ang_B = 1.0_dp
            end if
            if (use_pol_flux) then                                              ! poloidal flux as normal coordinate
                if (find_theta) then                                            ! derivatives in theta_V
                    dfun_ang_B = (1+dlam(1,1))/iotaf(kd)
                else                                                            ! derivatives in zeta_V
                    dfun_ang_B = -1 + dlam(1,1)/iotaf(kd)
                end if
            else                                                                ! toroidal flux as normal coordinate
                if (find_theta) then                                            ! derivatives in theta_V
                    dfun_ang_B = - (1+dlam(1,1))
                else                                                            ! derivatives in zeta_V
                    dfun_ang_B = iotaf(kd) - dlam(1,1)
                end if
            end if
        end function dfun_ang_B
    end function calc_ang_B
    
    ! normalizes equilibrium quantities pres_FD, q_saf_FD or rot_t_FD, flux_p_FD
    ! or flux_t_FD, max_flux and pres using the normalization constants
    subroutine normalize_eq_vars
        use num_vars, only: use_pol_flux
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! scale the quantities
        pres_FD = pres_FD/pres_0
        flux_p_FD = flux_p_FD/psi_0
        flux_t_FD = flux_t_FD/psi_0
        max_flux = max_flux/psi_0
        max_flux_eq = max_flux_eq/psi_0
        
        ! scale  the  derivatives  by  psi_p_0
        do id = 1,size(pres_FD,2)-1
            pres_FD(:,id) = pres_FD(:,id) * psi_0**(id)
            flux_p_FD(:,id) = flux_p_FD(:,id) * psi_0**(id)
            flux_t_FD(:,id) = flux_t_FD(:,id) * psi_0**(id)
            if (use_pol_flux) then
                q_saf_FD(:,id) = q_saf_FD(:,id) * psi_0**(id)
            else
                rot_t_FD(:,id) = rot_t_FD(:,id) * psi_0**(id)
            end if
        end do
    end subroutine normalize_eq_vars
    
    ! sets up normalization constants:
    !   R_0:    major radius (= average R on axis)
    !   rho_0:  mass density on axis (set up through input variable)
    !   pres_0: pressure on axis
    !   B_0:    average magnetic field (= sqrt(mu_0 pres_0))
    !   psi_0:  reference flux (= R_0^2 B_0)
    ! [MPI] only global master
    integer function calc_norm_const() result(ierr)
        use num_vars, only: glb_rank, eq_style, mu_0
        
        character(*), parameter :: rout_name = 'calc_norm_const'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    call calc_norm_const_VMEC
                case (2)                                                        ! HELENA
                    call calc_norm_const_HEL
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
    contains 
        ! VMEC version
        subroutine calc_norm_const_VMEC
            use VMEC_vars, only: R_c, presf
            
            ! set  the major  radius  as  the average  value  of  VMEC_R on  the
            ! magnetic axis
            R_0 = R_c(0,0,1,0)
            
            ! rho_0 is set up through an input variable with the same name
            
            ! set pres_0 as pressure on axis
            pres_0 = presf(1)
            
            ! set the reference value for B_0 from B_0 = sqrt(mu_0 pres_0)
            B_0 = sqrt(pres_0 * mu_0)
            
            ! set reference flux
            psi_0 = (R_0**2 * B_0)
        end subroutine calc_norm_const_VMEC
        
        ! HELENA version
        subroutine calc_norm_const_HEL
            use HEL_vars, only: R_H, p0
            
            ! set the major  radius as the average value of  VMEC_R on the
            ! magnetic axis
            R_0 = R_H
            
            ! rho_0 is set up through an input variable with the same name
            
            ! set pres_0 as pressure on axis
            pres_0 = p0(1)
            
            ! set the reference value for B_0 from B_0 = sqrt(mu_0 pres_0)
            B_0 = sqrt(pres_0 * mu_0)
            
            ! set reference flux
            psi_0 = (R_0**2 * B_0)
        end subroutine calc_norm_const_HEL
    end function calc_norm_const
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! equilibrium phase
    subroutine dealloc_eq
        deallocate(VMEC_R,VMEC_Z,VMEC_L)
        deallocate(flux_p_V)
        deallocate(flux_t_V)
        deallocate(trigon_factors)
    end subroutine dealloc_eq
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! calculation for a certain alpha
    subroutine dealloc_eq_final
        use num_vars, only: grp_rank, use_pol_flux
        
        deallocate(pres_V,pres_FD)
        if (use_pol_flux) then
            deallocate(q_saf_V,q_saf_FD)
        else
            deallocate(rot_t_V,rot_t_FD)
        end if
        deallocate(flux_p_FD)
        deallocate(flux_t_FD)
        
        ! group masters
        if (grp_rank.eq.0) then
            deallocate(flux_p_V_full)
            deallocate(flux_t_V_full)
            deallocate(q_saf_V_full)
            deallocate(rot_t_V_full)
        end if
        deallocate(zeta_V,theta_V,ang_par_F)
    end subroutine dealloc_eq_final
    
    ! calculates X,Y  and Z on  a grid in the  VMEC poloidal and  toroidal angle
    ! theta_V and  zeta_V, for every normal  point in either the  equilibrium or
    ! the perturbation grid,  indicated by the variable  eq_grid. The dimensions
    ! are (n_theta,n_zeta,n_r)
    ! This  routine also optionally  calculates lambda on  the grid, as  this is
    ! also needed some times.
    ! Note:  If  eq_grid  =  .false.,  the flux  arrays  have  to  be  correctly
    ! initialized for this routine to work
    integer function calc_XYZ_grid(theta,zeta,r_range,eq_grid,X,Y,Z,L) &
        &result(ierr)
        use fourier_ops, only: calc_trigon_factors, fourier2real
        use VMEC_vars, only: R_c, R_s, Z_c, Z_s, L_c, L_s, mpol, ntor
        use num_vars, only: use_pol_flux
        use utilities, only: interp_fun_1D
        
        character(*), parameter :: rout_name = 'calc_XYZ_grid'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:,:), zeta(:,:,:)                       ! points at which to calculate the grid
        real(dp), intent(in) :: r_range(2)                                      ! min. and max. normal range in eq. range or pert. range
        logical, intent(in) :: eq_grid                                          ! .true. if eq. grid is used, .false. if pert. grid
        real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:)    ! X, Y and Z of grid
        real(dp), intent(inout), allocatable, optional :: L(:,:,:)              ! lambda of grid
        
        ! local variables
        real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      ! trigonometric factor cosine for the inverse fourier transf.
        real(dp), allocatable :: R(:,:,:)                                       ! R in Cylindrical coordinates
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_theta, n_zeta, n_r                                         ! dimensions of the grid
        integer :: id                                                           ! counters
        real(dp), allocatable :: R_c_int(:,:,:), R_s_int(:,:,:)                 ! interpolated version of R_c and R_s
        real(dp), allocatable :: Z_c_int(:,:,:), Z_s_int(:,:,:)                 ! interpolated version of Z_c and Z_s
        real(dp), allocatable :: L_c_int(:,:,:), L_s_int(:,:,:)                 ! interpolated version of L_c and L_s
        real(dp) :: r_loc                                                       ! current r value in range specified by eq_grid
        real(dp) :: r_loc_eq                                                    ! current r value in equilibrium range
        real(dp), pointer :: flux(:), flux_VMEC(:)                              ! either pol. or tor. flux
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (size(theta,1).ne.size(zeta,1) .or. size(theta,2).ne.size(zeta,2) &
            &.or. size(theta,3).ne.size(zeta,3)) then
            ierr = 1
            err_msg =  'Theta and Zeta need to have the same dimensions'
            CHCKERR(err_msg)
        end if
        
        ! set up flux and flux_VMEC
        if (.not.eq_grid) then
            if (use_pol_flux) then
                flux => flux_p_FD(:,0)
            else
                flux => flux_t_FD(:,0)
            end if
            if (eq_use_pol_flux) then
                flux_VMEC => flux_p_FD(:,0)
            else
                flux_VMEC => flux_t_FD(:,0)
            end if
        end if
        
        ! set up n_theta, n_zeta and n_r
        n_theta = size(theta,1)
        n_zeta = size(theta,2)
        n_r = size(theta,3)
        
        ! initialize R and Z
        allocate(R(n_theta,n_zeta,n_r),Z(n_theta,n_zeta,n_r))
        allocate(X(n_theta,n_zeta,n_r),Y(n_theta,n_zeta,n_r))
        if (present(L)) allocate(L(n_theta,n_zeta,n_r))
        
        ! set up interpolated R_c_int, ..
        allocate(R_c_int(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(R_s_int(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(Z_c_int(0:mpol-1,-ntor:ntor,1:n_r))
        allocate(Z_s_int(0:mpol-1,-ntor:ntor,1:n_r))
        if (present(L)) then
            allocate(L_c_int(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(L_s_int(0:mpol-1,-ntor:ntor,1:n_r))
        end if
        
        ! interpolate for every requested normal point
        r_loc = r_range(1)
        do id = 1,n_r                                                           ! loop over all normal points
            ! convert r_loc to r_loc_eq if necessary
            if (eq_grid) then
                r_loc_eq = r_loc
            else
                ierr = interp_fun_1D(r_loc_eq,flux_VMEC/max_flux_eq,&
                    &r_loc,flux/max_flux)
                CHCKERR('')
            end if
            ! interpolate R_c, ... tables to fill R_c_int, ...
            call interp_VMEC_table(R_c(:,:,:,0),R_c_int(:,:,id),r_loc_eq)
            call interp_VMEC_table(R_s(:,:,:,0),R_s_int(:,:,id),r_loc_eq)
            call interp_VMEC_table(Z_c(:,:,:,0),Z_c_int(:,:,id),r_loc_eq)
            call interp_VMEC_table(Z_s(:,:,:,0),Z_s_int(:,:,id),r_loc_eq)
            if (present(L)) then
                call interp_VMEC_table(L_c(:,:,:,0),L_c_int(:,:,id),r_loc_eq)
                call interp_VMEC_table(L_s(:,:,:,0),L_s_int(:,:,id),r_loc_eq)
            end if
            ! increment r_loc
            r_loc = r_loc + 1._dp*(r_range(2)-r_range(1))/(n_r-1)
        end do
        
        ! inverse fourier transform with trigonometric factors
        do id = 1,n_theta
            ! calculate trigonometric factors
            ierr = calc_trigon_factors(theta(id,:,:),zeta(id,:,:),&
                &trigon_factors)
            CHCKERR('')
            ierr = fourier2real(R_c_int,R_s_int,trigon_factors,R(id,:,:),[0,0])
            CHCKERR('')
            ierr = fourier2real(Z_c_int,Z_s_int,trigon_factors,Z(id,:,:),[0,0])
            CHCKERR('')
            if (present(L)) then
                ierr = fourier2real(L_c_int,L_s_int,trigon_factors,L(id,:,:),&
                    &[0,0])
                CHCKERR('')
            end if
            
            ! transform cylindrical to cartesian
            ! (the geometrical zeta is the inverse of VMEC zeta)
            X(id,:,:) = R(id,:,:)*cos(-zeta(id,:,:))
            Y(id,:,:) = R(id,:,:)*sin(-zeta(id,:,:))
            
            ! deallocate
            deallocate(trigon_factors)
        end do
        
        ! deallocate
        deallocate(R_c_int,R_s_int,Z_c_int,Z_s_int)
        if (present(L)) deallocate(L_c_int,L_s_int)
        deallocate(R)
    contains
        ! interpolates an equilibrium table (R_c, Z_s, ...)
        subroutine interp_VMEC_table(var_in,var_out,pt_in)
            use utilities, only: con2dis
            
            ! input / output
            real(dp), intent(in) :: var_in(:,:,:)                               ! R_c, Z_s, ... VMEC table
            real(dp), intent(inout) :: var_out(:,:)                             ! interpolated version of R_c, Z_s, ...
            real(dp), intent(in) :: pt_in                                       ! point at which to interpolate (in equilibrium range)
            
            ! local variables
            integer :: id_lo, id_hi                                             ! low and high indices for the interpolation
            real(dp) :: pt_dis                                                  ! discrete, unrounded equivalent of pt_in in eq. grid
            
            ! convert pt_in to discrete equilibrium grid, unrounded
            call con2dis(pt_in,[0._dp,1._dp],pt_dis,[1,n_r_eq])
            
            ! round up and down
            id_lo = floor(pt_dis)
            id_hi = ceiling(pt_dis)
            
            ! limit both to the total range of the VMEC tables
            ! (they can fall outside the range due to numerical errors)
            if (id_lo.lt.1) id_lo = 1
            if (id_hi.lt.1) id_hi = 1
            if (id_lo.gt.n_r_eq) id_lo = n_r_eq
            if (id_hi.gt.n_r_eq) id_hi = n_r_eq
            
            ! interpolate var_in 
            var_out = var_in(:,:,id_lo) + (pt_dis-id_lo) * &
                &(var_in(:,:,id_hi)-var_in(:,:,id_lo))                          ! because id_hi - id_lo = 1
        end subroutine interp_VMEC_table
    end function calc_XYZ_grid
end module eq_vars
