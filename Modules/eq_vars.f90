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
        &adapt_HEL_to_eq, &
        &theta_V, zeta_V, n_par, grp_n_r_eq, lam_H, min_par, max_par, n_r_eq, &
        &q_saf_V, q_saf_V_full, flux_p_V, flux_t_V, VMEC_R, VMEC_Z, VMEC_L, &
        &pres_V, q_saf_FD, flux_p_FD, flux_t_FD, pres_FD, grp_min_r_eq, &
        &grp_max_r_eq, R_0, pres_0, B_0, psi_0, rho_0, T_0, ang_par_F, &
        &rot_t_FD, rot_t_V, flux_p_V_full, flux_t_V_full, rot_t_V_full, &
        &max_flux, max_flux_eq, eq_use_pol_flux, q_saf_H, rot_t_H, theta_H, &
        &zeta_H, pres_H, flux_p_H, flux_t_H, q_saf_H_full, rot_t_H_full, &
        &flux_p_H_full, flux_t_H_full, max_flux_F, max_flux_eq_F

    ! global variables
    ! Note: The indices in [derivatives] are:
    !   [VMEC_r,VMEC_theta,VMEC_zeta]   for VMEC variables
    !   [r,theta_F,zeta_F]              for F(lux) variables
    !   [r_H,theta_H,zeta_H]            for H(ELENA) variables
    real(dp), allocatable :: trigon_factors(:,:,:,:,:)                          ! trigonometric factor cosine for the inverse fourier transf.
    real(dp), allocatable :: VMEC_R(:,:,:,:,:)                                  ! R in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: VMEC_Z(:,:,:,:,:)                                  ! Z in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: VMEC_L(:,:,:,:,:)                                  ! L(ambda) in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: lam_H(:,:,:)                                       ! lambda in (HM)
    real(dp), allocatable :: theta_V(:,:), zeta_V(:,:)                          ! grid points (n_par, grp_n_r_eq, 2) in VMEC coords
    real(dp), allocatable :: theta_H(:,:), zeta_H(:,:)                          ! grid points (n_par, grp_n_r_eq, 2) in HELENA coords
    real(dp), allocatable :: ang_par_F(:,:)                                     ! parallel angle in flux coordinates (either zeta_F or theta_F)
    real(dp), allocatable, target :: q_saf_V(:,:)                               ! safety factor in VMEC coordinates
    real(dp), allocatable, target :: q_saf_H(:,:)                               ! safety factor in HELENA coordinates
    real(dp), allocatable :: q_saf_V_full(:,:)                                  ! safety factor in full normal mesh in VMEC coordinates
    real(dp), allocatable :: q_saf_H_full(:,:)                                  ! safety factor in full normal mesh in HELENA coordinates
    real(dp), allocatable :: q_saf_FD(:,:)                                      ! safety factor, Deriv. in Flux coords.
    real(dp), allocatable, target :: rot_t_V(:,:)                               ! rot. transform in VMEC coordinates
    real(dp), allocatable, target :: rot_t_H(:,:)                               ! rot. transform in HELENA coordinates
    real(dp), allocatable :: rot_t_V_full(:,:)                                  ! rot. transform in full normal mesh in VMEC coordinates
    real(dp), allocatable :: rot_t_H_full(:,:)                                  ! rot. transform in full normal mesh in HELENA coordinates
    real(dp), allocatable :: rot_t_FD(:,:)                                      ! rot. transform, Deriv. in Flux coords.
    real(dp), allocatable, target :: flux_p_V(:,:), flux_t_V(:,:), pres_V(:,:)  ! pol. flux, tor. flux and pressure, and norm. Deriv. in VMEC coords.
    real(dp), allocatable, target :: flux_p_H(:,:), flux_t_H(:,:), pres_H(:,:)  ! pol. flux, tor. flux and pressure, and norm. Deriv. in HELENA coords.
    real(dp), allocatable :: pres_FD(:,:)                                       ! pressure, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable, target :: flux_p_FD(:,:), flux_t_FD(:,:)             ! pol. and tor. flux, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable :: flux_p_V_full(:,:), flux_t_V_full(:,:)             ! pol. flux, tor. flux, and norm. Deriv. values and Derivs. in VMEC coords.
    real(dp), allocatable :: flux_p_H_full(:,:), flux_t_H_full(:,:)             ! pol. flux, tor. flux, and norm. Deriv. values and Derivs. in HELENA coords.
    real(dp) :: max_flux, max_flux_eq                                           ! max. flux (pol. or tor.) in Equilibrium coordinates
    real(dp) :: max_flux_F, max_flux_eq_F                                       ! max. flux (pol. or tor.) in Flux coordinates
    real(dp) :: min_par, max_par                                                ! min. and max. of parallel coordinate [pi]
    real(dp) :: R_0, pres_0, rho_0                                              ! independent normalization constants for nondimensionalization
    real(dp) :: B_0, psi_0, T_0                                                 ! derived normalization constants for nondimensionalization
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
        allocate(pres_FD(grp_n_r_eq,0:max_deriv))
        
        ! flux_p_FD
        allocate(flux_p_FD(grp_n_r_eq,0:max_deriv))
        
        ! flux_t_FD
        allocate(flux_t_FD(grp_n_r_eq,0:max_deriv))
        
        if (use_pol_flux) then
            ! q_saf_FD
            allocate(q_saf_FD(grp_n_r_eq,0:max_deriv))
        else
            ! rot_t_FD
            allocate(rot_t_FD(grp_n_r_eq,0:max_deriv))
        end if
        
        ! initialize variables that are  specifici to which equilibrium style is
        ! being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! R
                allocate(VMEC_R(n_par,grp_n_r_eq,0:max_deriv+1,0:max_deriv+1,&
                    &0:max_deriv+1))
                
                ! Z
                allocate(VMEC_Z(n_par,grp_n_r_eq,0:max_deriv+1,0:max_deriv+1,&
                    &0:max_deriv+1))
                
                ! lambda
                allocate(VMEC_L(n_par,grp_n_r_eq,0:max_deriv+1,0:max_deriv+1,&
                    &0:max_deriv+1))
                
                ! pres_V
                allocate(pres_V(grp_n_r_eq,0:max_deriv+1))
                
                ! flux_p_V
                allocate(flux_p_V(grp_n_r_eq,0:max_deriv+1))
                
                ! flux_t_V
                allocate(flux_t_V(grp_n_r_eq,0:max_deriv+1))
                
                if (use_pol_flux) then
                    ! q_saf_V
                    allocate(q_saf_V(grp_n_r_eq,0:max_deriv+1))
                else
                    ! rot_t_V
                    allocate(rot_t_V(grp_n_r_eq,0:max_deriv+1))
                end if
                
                ! full variables for group masters
                if (grp_rank.eq.0) then
                    ! flux_p_V_full
                    allocate(flux_p_V_full(n_r_eq,0:max_deriv+1))
                    
                    ! flux_t_V_full
                    allocate(flux_t_V_full(n_r_eq,0:max_deriv+1))
                    
                    ! q_saf_V_full
                    allocate(q_saf_V_full(n_r_eq,0:max_deriv+1))
                    
                    ! rot_t_V_full
                    allocate(rot_t_V_full(n_r_eq,0:max_deriv+1))
                end if
            case (2)                                                            ! HELENA
                ! pres_H
                allocate(pres_H(grp_n_r_eq,0:max_deriv+1))
                
                ! flux_p_H
                allocate(flux_p_H(grp_n_r_eq,0:max_deriv+1))
                
                ! flux_t_H
                allocate(flux_t_H(grp_n_r_eq,0:max_deriv+1))
                
                if (use_pol_flux) then
                    ! q_saf_H
                    allocate(q_saf_H(grp_n_r_eq,0:max_deriv+1))
                else
                    ! rot_t_H
                    allocate(rot_t_H(grp_n_r_eq,0:max_deriv+1))
                end if
                
                ! full variables for group masters
                if (grp_rank.eq.0) then
                    ! flux_p_H_full
                    allocate(flux_p_H_full(n_r_eq,0:max_deriv+1))
                    
                    ! flux_t_H_full
                    allocate(flux_t_H_full(n_r_eq,0:max_deriv+1))
                    
                    ! q_saf_H_full
                    allocate(q_saf_H_full(n_r_eq,0:max_deriv+1))
                    
                    ! rot_t_H_full
                    allocate(rot_t_H_full(n_r_eq,0:max_deriv+1))
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
        ierr = check_deriv(deriv,max_deriv+1,'calc_RZL')
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

    ! Adapt the HELENA  quantities for the equilibrium parallel  grid taking the
    ! correct poloidal HELENA range and possibly multiples.
    integer function adapt_HEL_to_eq() result(ierr)
        use num_vars, only: pi
        use HEL_vars, only: h_H_11, h_H_12, h_H_33, ias, chi_H
        use utilities, only: interp_fun_1D
        
        character(*), parameter :: rout_name = 'adapt_HEL_to_eq'
        
        ! local variables
        real(dp), allocatable :: old_h_H_11(:,:)                                ! upper metric factor 11 (gem11)
        real(dp), allocatable :: old_h_H_12(:,:)                                ! upper metric factor 12 (gem12)
        real(dp), allocatable :: old_h_H_33(:,:)                                ! upper metric factor 33 (1/gem33)
        integer :: id, kd                                                       ! counters
        real(dp) :: par_loc                                                     ! local parallel (= poloidal) point
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'HELENA SHOULD BE ADAPTED FOR EVERY ALPHA, BUT THE HELENA OUTPUT SHOULD BE SAVED !!!!'
        
        ! set old arrays
        allocate(old_h_H_11(size(h_H_11,1),size(h_H_11,2)))
        old_h_H_11 = h_H_11
        allocate(old_h_H_12(size(h_H_12,1),size(h_H_12,2)))
        old_h_H_12 = h_H_12
        allocate(old_h_H_33(size(h_H_33,1),size(h_H_33,2)))
        old_h_H_33 = h_H_33
        
        ! reallocate the new arrays
        deallocate(h_H_11); allocate(h_H_11(n_par,n_r_eq))
        deallocate(h_H_12); allocate(h_H_12(n_par,n_r_eq))
        deallocate(h_H_33); allocate(h_H_33(n_par,n_r_eq))
        
        ! For every poloidal point, check  which half poloidal circle it belongs
        ! to.  If this  is a  bottom part  and HELENA  is symmetric  (ias =  0),
        ! the  quantities have  to  be taken  from  their symmetric  counterpart
        ! (2pi-theta) and  the metric factors  h_H_12 carry an  additional minus
        ! sign.
        ! Note that min_par and max_par have units of [pi]
        do id = 1,n_par
            ! loop over all normal points
            do kd = 1,n_r_eq
                ! set the local poloidal point from theta_H
                par_loc = theta_H(id,kd)
                ! add or subtract 2pi to the parallel angle until it is at least
                ! 0 to get principal range 0..2pi
                if (par_loc.lt.0._dp) then
                    do while (par_loc.lt.0._dp)
                        par_loc = par_loc + 2*pi
                    end do
                else if (par_loc.gt.2*pi) then
                    do while (par_loc.gt.2._dp*pi)
                        par_loc = par_loc - 2*pi
                    end do
                end if
                ! Interpolate  the  HELENA  variables  poloidally,  taking  into
                ! account the possible symmetry
                if (ias.eq.0 .and. par_loc.gt.pi) then
                    ierr = interp_fun_1D(h_H_11(id,kd),old_h_H_11(:,kd),&
                        &2*pi-par_loc,x=chi_H)
                    CHCKERR('')
                    ierr = interp_fun_1D(h_H_12(id,kd),-old_h_H_12(:,kd),&
                        &2*pi-par_loc,x=chi_H)                                  ! change of sign
                    CHCKERR('')
                    ierr = interp_fun_1D(h_H_33(id,kd),old_h_H_33(:,kd),&
                        &2*pi-par_loc,x=chi_H)
                    CHCKERR('')
                else
                    ierr = interp_fun_1D(h_H_11(id,kd),old_h_H_11(:,kd),&
                        &par_loc,x=chi_H)
                    CHCKERR('')
                    ierr = interp_fun_1D(h_H_12(id,kd),old_h_H_12(:,kd),&
                        &par_loc,x=chi_H)
                    CHCKERR('')
                    ierr = interp_fun_1D(h_H_33(id,kd),old_h_H_33(:,kd),&
                        &par_loc,x=chi_H)
                    CHCKERR('')
                end if
            end do
        end do
        
        ! deallocate old arrays
        deallocate(old_h_H_11,old_h_H_12,old_h_H_33)
    end function adapt_HEL_to_eq
    
    ! calculates flux quantities  and normal derivatives in  the VMEC coordinate
    ! system
    integer function calc_flux_q() result(ierr)
        use num_vars, only: eq_style, max_deriv, grp_rank, use_pol_flux, &
            &glb_rank, plot_flux_q
        use utilities, only: calc_deriv, calc_int
        
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
        
        ! plot flux quantities if requested
        if (plot_flux_q .and. glb_rank.eq.0) then
            ierr = flux_q_plot()
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
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
            do kd = 1, max_deriv+1
                ierr = calc_deriv(pres_V(:,0),pres_V(:,kd),n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! set up helper variables to calculate poloidal flux
            allocate(Dflux_p_full(n_r_eq),flux_p_int_full(n_r_eq))
            Dflux_p_full = iotaf*phi_r
            ierr = calc_int(Dflux_p_full,1.0_dp/(n_r_eq-1.0_dp),flux_p_int_full)
            CHCKERR('')
            
            ! poloidal flux: calculate using iotaf and phi, phi_r
            ! easier to use full normal mesh flux_p because of the integral
            flux_p_V(:,1) = Dflux_p_full(grp_min_r_eq:grp_max_r_eq)
            flux_p_V(:,0) = flux_p_int_full(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(flux_p_V(:,1),flux_p_V(:,kd),&
                    &n_r_eq-1._dp,kd-1,1)
                CHCKERR('')
            end do
                
            ! toroidal flux: copy from VMEC and derive
            flux_t_V(:,0) = phi(grp_min_r_eq:grp_max_r_eq)
            flux_t_V(:,1) = phi_r(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(flux_t_V(:,1),flux_t_V(:,kd),n_r_eq-1._dp,&
                    &kd-1,1)
                CHCKERR('')
            end do
            
            if (use_pol_flux) then
                ! safety factor
                q_saf_V(:,0) = 1.0_dp/iotaf(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(q_saf_V(:,0),q_saf_V(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_p_int_full(n_r_eq)
                max_flux_F = max_flux
            else
                ! rot. transform
                rot_t_V(:,0) = iotaf(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(rot_t_V(:,0),rot_t_V(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = phi(n_r_eq)
                max_flux_F = - max_flux                                         ! conversion VMEC LH -> RH coord. system
            end if
            
            ! max_flux_eq
            if (eq_use_pol_flux) then
                max_flux_eq = flux_p_int_full(n_r_eq)
                max_flux_eq_F = max_flux_eq
            else
                max_flux_eq = phi(n_r_eq)
                max_flux_eq_F = - max_flux_eq                                   ! conversion VMEC LH -> RH coord. system
            end if
                
            ! the global master needs flux_p_V_full, flux_t_V_full, q_saf_V_full
            ! and rot_t_V_full  for the resonance plot  and checking of m  and n
            ! and  the  group masters  need  q_saf_V  on  full normal  mesh  for
            ! plot_X_vec
            if (grp_rank.eq.0) then
                ! flux_t_V_full
                flux_t_V_full(:,0) = phi
                flux_t_V_full(:,1) = phi_r
                do kd = 2,max_deriv+1
                    ierr = calc_deriv(flux_t_V_full(:,1),flux_t_V_full(:,kd),&
                        &n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! flux_p_V_full
                flux_p_V_full(:,0) = flux_p_int_full
                flux_p_V_full(:,1) = Dflux_p_full
                do kd = 2,max_deriv+1
                    ierr = calc_deriv(flux_p_V_full(:,1),&
                        &flux_p_V_full(:,kd),n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! q_saf_V_full
                q_saf_V_full(:,0) = 1.0_dp/iotaf
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(q_saf_V_full(:,0),&
                        &q_saf_V_full(:,kd),n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
                
                ! rot_t_V_full
                rot_t_V_full(:,0) = iotaf
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(rot_t_V_full(:,0),&
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
            do kd = 1, max_deriv+1
                ierr = calc_deriv(pres_H(:,0),pres_H(:,kd),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),kd,1)
                CHCKERR('')
            end do
            
            ! set up helper variables to calculate toroidal flux
            ! calculate normal derivative of flux_H
            allocate(flux_H_r(n_r_eq))
            ierr = calc_deriv(flux_H,flux_H_r,flux_H,1,1)
            CHCKERR('')
            allocate(Dflux_t_full(n_r_eq),flux_t_int_full(n_r_eq))
            Dflux_t_full = qs*flux_H_r
            ierr = calc_int(Dflux_t_full,flux_H,flux_t_int_full)
            CHCKERR('')
            
            ! toroidal flux: calculate using qs and flux_H, flux_H_r
            ! easier to use full normal mesh flux_t because of the integral
            flux_t_H(:,1) = Dflux_t_full(grp_min_r_eq:grp_max_r_eq)
            flux_t_H(:,0) = flux_t_int_full(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(flux_t_H(:,1),flux_t_H(:,kd),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),kd-1,1)
                CHCKERR('')
            end do
                
            ! poloidal flux: copy from HELENA and derive
            flux_p_H(:,0) = flux_H(grp_min_r_eq:grp_max_r_eq)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(flux_p_H(:,0),flux_p_H(:,kd),&
                    &flux_H(grp_min_r_eq:grp_max_r_eq),kd,1)
                CHCKERR('')
            end do
            
            if (use_pol_flux) then
                ! safety factor
                q_saf_H(:,0) = qs(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(q_saf_H(:,0),q_saf_H(:,kd),&
                        &flux_H(grp_min_r_eq:grp_max_r_eq),kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_H(n_r_eq)
            else
                ! rot. transform
                rot_t_H(:,0) = 1.0_dp/qs(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(rot_t_H(:,0),rot_t_H(:,kd),&
                        &flux_H(grp_min_r_eq:grp_max_r_eq),kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_t_int_full(n_r_eq)
            end if
            max_flux_F = max_flux
            
            ! max_flux_eq
            max_flux_eq = flux_H(n_r_eq)
            max_flux_eq_F = max_flux_eq
            
            ! the global master needs flux_p_H_full, flux_t_H_full, q_saf_H_full
            ! and rot_t_H_full  for the resonance plot  and checking of m  and n
            ! and  the  group masters  need  q_saf_H  on  full normal  mesh  for
            ! plot_X_vec
            if (grp_rank.eq.0) then
                ! flux_t_H_full
                flux_t_H_full(:,0) = flux_t_int_full
                flux_t_H_full(:,1) = Dflux_t_full
                do kd = 2,max_deriv
                    ierr = calc_deriv(flux_t_H_full(:,1),flux_t_H_full(:,kd),&
                        &flux_H,kd-1,1)
                    CHCKERR('')
                end do
                
                ! flux_p_H_full
                flux_p_H_full(:,0) = flux_H
                do kd = 1,max_deriv
                    ierr = calc_deriv(flux_p_H_full(:,0),flux_p_H_full(:,kd),&
                        &flux_H,kd,1)
                    CHCKERR('')
                end do
                
                ! q_saf_H_full
                q_saf_H_full(:,0) = qs(:)
                do kd = 1,max_deriv
                    ierr = calc_deriv(q_saf_H_full(:,0),q_saf_H_full(:,kd),&
                        &flux_H,kd,1)
                    CHCKERR('')
                end do
                
                ! rot_t_H_full
                rot_t_H_full(:,0) = 1.0_dp/qs(:)
                do kd = 1,max_deriv
                    ierr = calc_deriv(rot_t_H_full(:,0),rot_t_H_full(:,kd),&
                        &flux_H,kd,1)
                    CHCKERR('')
                end do
            end if
            
            ! deallocate helper variables
            if (use_pol_flux .or. grp_rank.eq.0) then
                deallocate(Dflux_t_full,flux_t_int_full)
            end if
        end function calc_flux_q_HEL
    end function calc_flux_q
    
    ! plots the flux quantities in the perturbation grid
    !   safety factor q_saf
    !   rotational transform rot
    !   pressure pres
    !   poloidal flux flux_p
    !   toroidal flux flux_t
    integer function flux_q_plot() result(ierr)
        use num_vars, only: eq_style, use_pol_flux, output_style
        use VMEC_vars, only: presf
        use HEL_vars, only: p0, flux_H
        
        character(*), parameter :: rout_name = 'flux_q_plot'
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_vars = 5                                                   ! nr. of variables to plot
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: x_plot_2D(:,:)                                 ! x values of 2D plot
        real(dp), allocatable :: y_plot_2D(:,:)                                 ! y values of 2D plot
        character(len=max_str_ln), allocatable :: plot_titles(:)                ! plot titles
        character(len=max_str_ln), allocatable :: file_names(:)                 ! file_names
        real(dp), allocatable :: r_plot(:)                                      ! normal r for plot
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Plotting flux quantities')
        
        call lvl_ud(1)
        
        ! initialize x_plot_2D and y_plot_2D
        allocate(x_plot_2D(n_r_eq,n_vars))
        allocate(y_plot_2D(n_r_eq,n_vars))
        
        ! set up plot titles and file names
        allocate(plot_titles(n_vars))
        allocate(file_names(2))
        plot_titles(1) = 'safety factor []'
        plot_titles(2) = 'rotational transform []'
        plot_titles(3) = 'pressure [pa]'
        plot_titles(4) = 'poloidal flux [Tm^2]'
        plot_titles(5) = 'toroidal flux [Tm^2]'
        file_names(1) = 'pres'
        file_names(2) = 'flux'
        
        ! fill the 2D version of the plot
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                y_plot_2D(:,1) = -q_saf_V_full(:,0)                             ! conversion VMEC LH -> RH coord. system
                y_plot_2D(:,2) = -rot_t_V_full(:,0)                             ! conversion VMEC LH -> RH coord. system
                y_plot_2D(:,3) = presf
                y_plot_2D(:,4) = flux_p_V_full(:,0)
                y_plot_2D(:,5) = -flux_t_V_full(:,0)                            ! conversion VMEC LH -> RH coord. system
                ! 2D normal variable (y_plot_2D tabulated in eq. grid)
                if (use_pol_flux) then
                    x_plot_2D(:,1) = flux_p_V_full(:,0)/max_flux
                else
                    x_plot_2D(:,1) = flux_t_V_full(:,0)/max_flux
                end if
                do id = 2,n_vars
                    x_plot_2D(:,id) = x_plot_2D(:,1)
                end do
            case (2)                                                            ! HELENA
                y_plot_2D(:,1) = q_saf_H_full(:,0)
                y_plot_2D(:,2) = rot_t_H_full(:,0)
                y_plot_2D(:,3) = p0
                y_plot_2D(:,4) = flux_p_H_full(:,0)
                y_plot_2D(:,5) = flux_t_H_full(:,0)
                if (use_pol_flux) then
                    x_plot_2D(:,1) = flux_p_H_full(:,0)/max_flux
                else
                    x_plot_2D(:,1) = flux_t_H_full(:,0)/max_flux
                end if
                do id = 2,n_vars
                    x_plot_2D(:,id) = x_plot_2D(:,1)
                end do
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! plot the output
        ! choose which output style is being used:
        !   1:  GNUPlot
        !   2:  HDF5
        select case (output_style)
            case (1)                                                            ! GNUPlot
                ! plot  the  2D output  (except  q_saf and  rot_t,  as they  are
                ! plotalready ted in plot_jq)
                call writo('The safety factor and rotational transform are not &
                    &plotted here. Instead, use the input variable "plot_jq".')
                call print_GP_2D(plot_titles(3),file_names(1)//'.dat',&
                    &y_plot_2D(:,3),x_plot_2D(:,3),draw=.false.)
                call draw_GP(plot_titles(3),file_names(1),1,.true.,.false.)
                call print_GP_2D(trim(plot_titles(4))//', '//&
                    &trim(plot_titles(5)),file_names(2),y_plot_2D(:,4:5),&
                    &x_plot_2D(:,4:5),draw=.false.)
                call draw_GP(trim(plot_titles(4))//', '//trim(plot_titles(5)),&
                    &file_names(2),2,.true.,.false.)
            case (2)                                                            ! HDF5
                ! set up r_plot
                allocate(r_plot(n_r_eq))
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        r_plot = [((id-1._dp)/(n_r_eq-1),id=1,n_r_eq)]          ! normal variable is equidistant normal poloidal flux
                    case (2)                                                    ! HELENA
                        r_plot = flux_H/flux_H(n_r_eq)                          ! output from HELENA
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ierr = flux_q_plot_HDF5(r_plot)
                CHCKERR('')
            case default
                err_msg = 'No output style associated with '//&
                    &trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! deallocate
        deallocate(x_plot_2D,y_plot_2D)
        
        call lvl_ud(-1)
    contains
        ! convert 2D plot to real plot in 3D and output in HDF5
        integer function flux_q_plot_HDF5(r_plot) result(ierr)
            use output_ops, only: print_HDF5
            use num_vars, only: n_theta_plot, n_zeta_plot
            
            character(*), parameter :: rout_name = 'flux_q_plot_HDF5'
            
            ! input / output
            real(dp), intent(in) :: r_plot(:)                                   ! normal r for 3D plot
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)        ! theta and zeta for 3D plot
            real(dp), allocatable :: x_plot_3D(:,:,:)                           ! x values of 3D plot
            real(dp), allocatable :: y_plot_3D(:,:,:)                           ! y values of 3D plot
            real(dp), allocatable :: z_plot_3D(:,:,:)                           ! z values of 3D plot
            real(dp), allocatable :: x_plot(:,:,:,:)                            ! x values of total plot
            real(dp), allocatable :: y_plot(:,:,:,:)                            ! y values of total plot
            real(dp), allocatable :: z_plot(:,:,:,:)                            ! z values of total plot
            real(dp), allocatable :: f_plot(:,:,:,:)                            ! values of variable of total plot
            integer :: n_r_plot                                                 ! how many normal points
            integer :: plot_dim(4)                                              ! total plot dimensions
            
            ! initialize ierr
            ierr = 0
            
            ! intitialize n_r_plot
            n_r_plot = size(r_plot)
            
            ! initialize theta_plot and zeta_plot
            allocate(theta_plot(n_theta_plot,n_zeta_plot,n_r_plot))
            if (n_theta_plot.eq.1) then
                theta_plot = 0.0_dp
            else
                do id = 1,n_theta_plot
                    theta_plot(id,:,:) = &
                        &pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)              ! starting from pi gives nicer plots
                end do
            end if
            ! zeta equidistant
            allocate(zeta_plot(n_theta_plot,n_zeta_plot,n_r_plot))
            if (n_zeta_plot.eq.1) then
                zeta_plot = 0.0_dp
            else
                do id = 1,n_zeta_plot
                    zeta_plot(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
                end do
            end if
            
            ! calculate X,Y and Z
            ierr = calc_XYZ_grid(theta_plot,zeta_plot,r_plot,&
                &x_plot_3D,y_plot_3D,z_plot_3D)
            CHCKERR('')
            
            ! set up plot_dim
            plot_dim = [n_theta_plot,n_zeta_plot,n_r_plot,n_vars]
            
            ! set up total plot variables
            allocate(x_plot(n_theta_plot,n_zeta_plot,n_r_plot,n_vars))
            allocate(y_plot(n_theta_plot,n_zeta_plot,n_r_plot,n_vars))
            allocate(z_plot(n_theta_plot,n_zeta_plot,n_r_plot,n_vars))
            allocate(f_plot(n_theta_plot,n_zeta_plot,n_r_plot,n_vars))
            do id = 1,n_vars
                x_plot(:,:,:,id) = x_plot_3D
                y_plot(:,:,:,id) = y_plot_3D
                z_plot(:,:,:,id) = z_plot_3D
            end do
            do kd = 1,n_r_plot
                f_plot(:,:,kd,1) = y_plot_2D(kd,1)                              ! safey factor
                f_plot(:,:,kd,2) = y_plot_2D(kd,2)                              ! rotational transform
                f_plot(:,:,kd,3) = y_plot_2D(kd,3)                              ! pressure
                f_plot(:,:,kd,4) = y_plot_2D(kd,4)                              ! poloidal flux
                f_plot(:,:,kd,5) = y_plot_2D(kd,5)                              ! toroidal flux
            end do
            
            ! print the output using HDF5
            call print_HDF5(plot_titles,'flux_quantities',f_plot,plot_dim,&
                &plot_dim,[0,0,0,0],x_plot,y_plot,z_plot,col=1,&
                &description='Flux quantities')
            
            ! deallocate
            deallocate(theta_plot,zeta_plot)
            deallocate(x_plot_3D,y_plot_3D,z_plot_3D)
            deallocate(x_plot,y_plot,z_plot,f_plot)
        end function flux_q_plot_HDF5
    end function flux_q_plot

    ! Calculate the angular mesh with  the normal variable being the equilibrium
    ! normal variable.
    ! Later, in check_and_limit_mesh, the normal  extent is limited if the tests
    ! are positive.
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
        ! zeta_V and ang_par_F are defined on  the entire normal mesh.
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
        integer function calc_mesh_HEL(alpha) result(ierr)
            use HEL_vars, only: qs
            use num_vars, only: use_pol_flux
            
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
            ! Note: this  includes chi_H,  or half the  array chi_H,  except for
            ! maybe  the  last point  if  HELENA  is top-bottom  symmetric  (see
            ! read_HEL)
            ierr = calc_eqd_mesh(ang_par_F(:,1),n_par,min_par,max_par)
            CHCKERR('')
            do kd = 2,n_r_eq
                ang_par_F(:,kd) = ang_par_F(:,1)
            end do
            
            ! calculate theta_H and zeta_H from alpha and ang_par_F
            if (use_pol_flux) then                                              ! theta_H is parallel coordinate
                theta_H = ang_par_F
                do id = 1,n_par
                    zeta_H(id,:) = ang_par_F(id,:)*qs(:)
                end do
                zeta_H = zeta_H + alpha
            else                                                                ! zeta_H is parallel coordinate
                zeta_H = ang_par_F
                do id = 1,n_par
                    theta_H(id,:) = ang_par_F(id,:)/qs(:)
                end do
                theta_H = theta_H - alpha
            end if
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

    ! calculate mesh of equidistant points,  where optionally the last point can
    ! be excluded
    ! Note: input is given in units of pi
    integer function calc_eqd_mesh(eqd_mesh,n_ang,min_ang,max_ang,excl_last) &
        &result(ierr)
        character(*), parameter :: rout_name = 'eqd_mesh'
        
        ! input and output
        real(dp), intent(inout) :: eqd_mesh(:)                                  ! output
        real(dp), intent(in) :: min_ang, max_ang                                ! min. and max. of angles [pi]
        integer, intent(in) :: n_ang                                            ! nr. of points
        logical, intent(in), optional :: excl_last                              ! .true. if last point excluded
        
        ! local variables
        integer :: id
        real(dp) :: delta_ang
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: excl_last_loc                                                ! local copy of excl_last
        
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
        
        ! set up local excl_last
        excl_last_loc = .false.
        if (present(excl_last)) excl_last_loc = excl_last
        
        ! initialize output vector
        eqd_mesh = 0.0_dp
        
        ! There are (n_ang-1) pieces in the total interval but if excl_last, the
        ! last one is not included
        if (excl_last_loc) then
            delta_ang = (max_ang-min_ang)/(n_ang) * pi
        else
            delta_ang = (max_ang-min_ang)/(n_ang-1) * pi
        end if
        
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
            
            ! do a check whether the result is indeed alpha
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
        max_flux_F = max_flux_F/psi_0
        max_flux_eq = max_flux_eq/psi_0
        max_flux_eq_F = max_flux_eq_F/psi_0
        
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
        use num_vars, only: glb_rank, eq_style, mu_0, use_normalization
        
        character(*), parameter :: rout_name = 'calc_norm_const'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        if (use_normalization) then
            ! user output
            call writo('Calculating the normalization constants')
            
            call lvl_ud(1)
            
            ! calculation
            if (glb_rank.eq.0) then
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        call calc_norm_const_VMEC
                    case (2)                                                    ! HELENA
                        call calc_norm_const_HEL
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            end if
            
            ! Alfven velocity
            T_0 = sqrt(mu_0*rho_0)*R_0/B_0 
            
            call writo('Major radius    R_0 = '//trim(r2strt(R_0))//' m')
            call writo('Pressure        pres_0 = '//trim(r2strt(pres_0))//' Pa')
            call writo('Mass density    rho_0 = '//trim(r2strt(rho_0))&
                &//' kg/m^3')
            call writo('Magnetic field  B_0 = '//trim(r2strt(B_0))//' T')
            call writo('Magnetic flux   psi_0 = '//trim(r2strt(psi_0))//' Tm^2')
            call writo('Alfven time     T_0 = '//trim(r2strt(T_0))//' s')
            
            call lvl_ud(-1)
            
            ! user output
            call writo('Normalization constants calculated')
        else
            ! user output
            call writo('Normalization not used')
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
            psi_0 = R_0**2 * B_0
        end subroutine calc_norm_const_VMEC
        
        ! HELENA version
        subroutine calc_norm_const_HEL
            use HEL_vars, only: R_0_H, B_0_H
            
            ! set the major radius as the HELENA normalization parameter
            R_0 = R_0_H
            
            ! rho_0 is set up through an input variable with the same name
            
            !  set the  reference  value  for B_0  as  the HELENA  normalization
            ! parameter
            B_0 = B_0_H
            
            ! set pres_0 as B_0^2/mu_0
            pres_0 = B_0**2/mu_0
            
            ! set reference flux
            psi_0 = R_0**2 * B_0
        end subroutine calc_norm_const_HEL
    end function calc_norm_const
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! equilibrium phase
    integer function dealloc_eq() result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'dealloc_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                deallocate(VMEC_R,VMEC_Z,VMEC_L)
                deallocate(flux_p_V)
                deallocate(flux_t_V)
                deallocate(trigon_factors)
            case (2)                                                            ! HELENA
                deallocate(flux_p_H)
                deallocate(flux_t_H)
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function dealloc_eq
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! calculation for a certain alpha
    integer function dealloc_eq_final() result(ierr)
        use num_vars, only: grp_rank, use_pol_flux, eq_style
        
        character(*), parameter :: rout_name = 'dealloc_eq_final'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! deallocate general variables
        deallocate(pres_FD)
        deallocate(flux_p_FD,flux_t_FD)
        if (use_pol_flux) then
            deallocate(q_saf_FD)
        else
            deallocate(rot_t_FD)
        end if
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                deallocate(pres_V)
                if (use_pol_flux) then
                    deallocate(q_saf_V)
                else
                    deallocate(rot_t_V)
                end if
                
                ! group masters
                if (grp_rank.eq.0) then
                    deallocate(flux_p_V_full)
                    deallocate(flux_t_V_full)
                    deallocate(q_saf_V_full)
                    deallocate(rot_t_V_full)
                end if
                deallocate(zeta_V,theta_V,ang_par_F)
            case (2)                                                            ! HELENA
                deallocate(pres_H)
                if (use_pol_flux) then
                    deallocate(q_saf_H)
                else
                    deallocate(rot_t_H)
                end if
                
                ! group masters
                if (grp_rank.eq.0) then
                    deallocate(flux_p_H_full)
                    deallocate(flux_t_H_full)
                    deallocate(q_saf_H_full)
                    deallocate(rot_t_H_full)
                end if
                deallocate(zeta_H,theta_H,ang_par_F)
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function dealloc_eq_final
    
    ! calculates X,Y  and Z on a  grid in the equilibrium  poloidal and toroidal
    ! angle theta and zeta, for every normal point in the equilibrium grid.
    ! The dimensions are (n_theta,n_zeta,n_r)
    ! If VMEC is the equilibrium  model, this routine also optionally calculates
    ! lambda on the grid, as this is  also needed some times. If HELENA is used,
    ! this variable is not allocated.
    ! Note: The variables X, Y, Z and optionally L have to be unallocated.
    integer function calc_XYZ_grid(theta,zeta,norm_r,X,Y,Z,L) &
        &result(ierr)
        use num_vars, only: eq_style
        use utilities, only: interp_fun_1D, round_with_tol
        
        character(*), parameter :: rout_name = 'calc_XYZ_grid'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:,:), zeta(:,:,:), norm_r(:)            ! points at which to calculate the grid
        real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:)    ! X, Y and Z of grid
        real(dp), intent(inout), allocatable, optional :: L(:,:,:)              ! lambda of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (size(theta,1).ne.size(zeta,1) .or. &
            &size(theta,2).ne.size(zeta,2) .or. &
            &size(theta,3).ne.size(zeta,3) .or. &
            &size(theta,3).ne.size(norm_r)) then
            ierr = 1
            err_msg =  'theta, zeta and norm_r need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_XYZ_grid_VMEC(theta,zeta,norm_r,X,Y,Z,L)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_XYZ_grid_HEL(theta,zeta,norm_r,X,Y,Z)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! VMEC version
        integer function calc_XYZ_grid_VMEC(theta,zeta,norm_r,X,Y,Z,L) &
            &result(ierr)
            use VMEC_vars, only: R_c, R_s, Z_c, Z_s, L_c, L_s, mpol, ntor
            use fourier_ops, only: calc_trigon_factors, fourier2real
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_VMEC'
            
            ! input / output
            real(dp), intent(in) :: theta(:,:,:), zeta(:,:,:)                   ! points at which to calculate the grid
            real(dp), intent(in) :: norm_r(:)                                   ! r values in equilibrium range
            real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), &
                &Z(:,:,:)                                                       ! X, Y and Z of grid
            real(dp), intent(inout), allocatable, optional :: L(:,:,:)          ! lambda of grid
            
            ! local variables
            integer :: id, kd, m, n                                             ! counters
            integer :: n_theta, n_zeta, n_r                                     ! dimensions of the grid
            real(dp), allocatable :: R_c_int(:,:,:), R_s_int(:,:,:)             ! interpolated version of R_c and R_s
            real(dp), allocatable :: Z_c_int(:,:,:), Z_s_int(:,:,:)             ! interpolated version of Z_c and Z_s
            real(dp), allocatable :: L_c_int(:,:,:), L_s_int(:,:,:)             ! interpolated version of L_c and L_s
            real(dp), allocatable :: trigon_factors(:,:,:,:,:)                  ! trigonometric factor cosine for the inverse fourier transf.
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            
            ! initialize ierr
            ierr = 0
            
            ! set up n_theta, n_zeta and n_r
            n_theta = size(theta,1)
            n_zeta = size(theta,2)
            n_r = size(theta,3)
            
            ! initialize R, Z, X and Y
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
            
            ! interpolate VMEC tables for every requested normal point
            do kd = 1,n_r                                                       ! loop over all normal points
                ! interpolate the VMEC tables in normal direction
                ! Note: VMEC uses an equidistant grid, here normalized from 0 to
                ! 1
                do n = -ntor,ntor
                    do m = 0,mpol-1
                        ierr = interp_fun_1D(R_c_int(m,n,kd),R_c(m,n,:,0),&
                            &norm_r(kd))
                        CHCKERR('')
                        ierr = interp_fun_1D(R_s_int(m,n,kd),R_s(m,n,:,0),&
                            &norm_r(kd))
                        CHCKERR('')
                        ierr = interp_fun_1D(Z_c_int(m,n,kd),Z_c(m,n,:,0),&
                            &norm_r(kd))
                        CHCKERR('')
                        ierr = interp_fun_1D(Z_s_int(m,n,kd),Z_s(m,n,:,0),&
                            &norm_r(kd))
                        CHCKERR('')
                        if (present(L)) then
                            ierr = interp_fun_1D(L_c_int(m,n,kd),L_c(m,n,:,0),&
                                &norm_r(kd))
                            CHCKERR('')
                            ierr = interp_fun_1D(L_s_int(m,n,kd),L_s(m,n,:,0),&
                                &norm_r(kd))
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! inverse fourier transform with trigonometric factors
            do id = 1,n_theta
                ! calculate trigonometric factors
                ierr = calc_trigon_factors(theta(id,:,:),zeta(id,:,:),&
                    &trigon_factors)
                CHCKERR('')
                ierr = fourier2real(R_c_int,R_s_int,trigon_factors,R(id,:,:),&
                        &[0,0])
                CHCKERR('')
                ierr = fourier2real(Z_c_int,Z_s_int,trigon_factors,Z(id,:,:),&
                        &[0,0])
                CHCKERR('')
                if (present(L)) then
                    ierr = fourier2real(L_c_int,L_s_int,trigon_factors,&
                        &L(id,:,:),[0,0])
                    CHCKERR('')
                end if
                
                ! deallocate
                deallocate(trigon_factors)
            end do
            
            ! transform cylindrical to cartesian
            ! (the geometrical zeta is the inverse of VMEC zeta)
            X = R*cos(-zeta)
            Y = R*sin(-zeta)
            
            ! deallocate
            deallocate(R_c_int,R_s_int,Z_c_int,Z_s_int)
            if (present(L)) deallocate(L_c_int,L_s_int)
            deallocate(R)
        end function calc_XYZ_grid_VMEC
        
        ! HELENA version
        integer function calc_XYZ_grid_HEL(theta,zeta,norm_r,X,Y,Z) result(ierr)
            use HEL_vars, only: R_H, Z_H, chi_H, ias, flux_H
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_HEL'
            
            ! input / output
            real(dp), intent(in) :: theta(:,:,:), zeta(:,:,:)                   ! points at which to calculate the grid
            real(dp), intent(in) :: norm_r(:)                                   ! r values in equilibrium range
            real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), &
                &Z(:,:,:)                                                       ! X, Y and Z of grid
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            integer :: n_theta, n_zeta, n_r                                     ! dimensions of the grid
            real(dp), allocatable :: R_H_int(:), Z_H_int(:)                     ! R and Z at interpolated normal value
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            real(dp) :: theta_loc                                               ! local copy of theta
            
            ! initialize ierr
            ierr = 0
            
            ! set up n_theta, n_zeta and n_r
            n_theta = size(theta,1)
            n_zeta = size(theta,2)
            n_r = size(theta,3)
            
            ! initialize R, Z X and Y
            allocate(R(n_theta,n_zeta,n_r),Z(n_theta,n_zeta,n_r))
            allocate(X(n_theta,n_zeta,n_r),Y(n_theta,n_zeta,n_r))
            
            ! set up interpolated R and Z
            allocate(R_H_int(size(R_H,1)),Z_H_int(size(Z_H,1)))
            
            ! interpolate HELENA output  R_H and Z_H for  every requested normal
            ! point
            do kd = 1,n_r                                                       ! loop over all normal points
                ! interpolate the HELENA tables in normal direction
                ! Note:  HELENA   uses  a regular,  non-equidistant  grid,  here
                ! normalized from 0 to 1 in flux_H/flux_H(n_r_eq)
                do id = 1,size(R_H,1)
                    ierr = interp_fun_1D(R_H_int(id),R_H(id,:),norm_r(kd),&
                        &x=flux_H/flux_H(n_r_eq))
                    CHCKERR('')
                end do
                do id = 1,size(Z_H,1)
                    ierr = interp_fun_1D(Z_H_int(id),Z_H(id,:),norm_r(kd),&
                        &x=flux_H/flux_H(n_r_eq))
                    CHCKERR('')
                end do
                ! loop over toroidal points
                do jd = 1,n_zeta
                    ! interpolate at the requested poloidal points
                    ! Note: R_H  and Z_H are  not adapted to the  parallel grid,
                    ! but tabulated in the original HELENA poloidal grid
                    do id = 1,n_theta
                        theta_loc = theta(id,jd,kd)
                        ! add or subtract 2pi to  the parallel angle until it is
                        ! at least 0 to get principal range 0..2pi
                        if (theta_loc.lt.0._dp) then
                            do while (theta_loc.lt.0._dp)
                                theta_loc = theta_loc + 2*pi
                            end do
                        else if (theta_loc.gt.2*pi) then
                            do while (theta_loc.gt.2._dp*pi)
                                theta_loc = theta_loc - 2*pi
                            end do
                        end if
                        ! Interpolate  the HELENA  variables poloidally,  taking
                        ! into account the possible symmetry
                        if (ias.eq.0 .and. theta_loc.gt.pi) then
                            ierr = interp_fun_1D(R(id,jd,kd),R_H_int,&
                                &2*pi-theta_loc,x=chi_H)
                            CHCKERR('')
                            ierr = interp_fun_1D(Z(id,jd,kd),-Z_H_int,&
                                &2*pi-theta_loc,x=chi_H)                        ! Z at second half of period is inverted
                            CHCKERR('')
                        else
                            ierr = interp_fun_1D(R(id,jd,kd),R_H_int,theta_loc,&
                                &x=chi_H)
                            CHCKERR('')
                            ierr = interp_fun_1D(Z(id,jd,kd),Z_H_int,theta_loc,&
                                &x=chi_H)
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! calculate X and Y, transforming cylindrical to cartesian
            X = R*cos(zeta)
            Y = R*sin(zeta)
            
            ! deallocate
            deallocate(R)
        end function calc_XYZ_grid_HEL
    end function calc_XYZ_grid
end module eq_vars
