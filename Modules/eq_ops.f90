!------------------------------------------------------------------------------!
!   Calculates the equilibrium quantities, making use of the metric_ops,       !
!   eq_vars, etc                                                               !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use num_vars, only: pi, dp, max_str_ln
    use message_ops, only: print_ar_2, lvl_ud, writo
    use output_ops, only: print_GP_3D, draw_GP_animated, draw_GP, &
        &print_GP_2D
    use str_ops, only: i2str, r2strt
    
    implicit none
    private
    public init_eq, read_eq, calc_flux_q, prepare_RZL, calc_RZL, &
        &normalize_eq_vars, calc_norm_const, adapt_HEL_to_eq
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! initialize the equilibrium variables
    integer function init_eq() result(ierr)
        use num_vars, only: max_deriv, grp_rank, use_pol_flux, eq_style
        use eq_vars, only: pres_FD, flux_p_FD, flux_t_FD, q_saf_FD, rot_t_FD, &
            &pres_E, flux_p_E, flux_t_E, q_saf_E, rot_t_E, flux_p_E_full, &
            &flux_t_E_full, q_saf_E_full, rot_t_E_full, VMEC_R, VMEC_Z, &
            &VMEC_L, grp_min_r_eq, grp_max_r_eq, grp_n_r_eq, n_par, n_r_eq
        
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
        
        ! pres_E
        allocate(pres_E(grp_n_r_eq,0:max_deriv+1))
        
        ! flux_p_E
        allocate(flux_p_E(grp_n_r_eq,0:max_deriv+1))
        
        ! flux_t_E
        allocate(flux_t_E(grp_n_r_eq,0:max_deriv+1))
        
        if (use_pol_flux) then
            ! q_saf_E
            allocate(q_saf_E(grp_n_r_eq,0:max_deriv+1))
        else
            ! rot_t_E
            allocate(rot_t_E(grp_n_r_eq,0:max_deriv+1))
        end if
        
        ! full variables for group masters
        if (grp_rank.eq.0) then
            ! flux_p_E_full
            allocate(flux_p_E_full(n_r_eq,0:max_deriv+1))
            
            ! flux_t_E_full
            allocate(flux_t_E_full(n_r_eq,0:max_deriv+1))
            
            ! q_saf_E_full
            allocate(q_saf_E_full(n_r_eq,0:max_deriv+1))
            
            ! rot_t_E_full
            allocate(rot_t_E_full(n_r_eq,0:max_deriv+1))
        end if
        
        ! initialize variables that are specificic to which equilibrium style is
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
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function init_eq
    
    ! reads the equilibrium input file
    integer function read_eq() result(ierr)
        use num_vars, only: eq_style, glb_rank
        use VMEC_ops, only: read_VMEC
        use HEL_ops, only: read_HEL
        use eq_vars, only: n_r_eq, eq_use_pol_flux
        
        character(*), parameter :: rout_name = 'read_eq'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! only do this for the group master
        if (glb_rank.eq.0) then
            ! choose which equilibrium style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ierr = read_VMEC(n_r_eq,eq_use_pol_flux)
                    CHCKERR('')
                case (2)                                                        ! HELENA
                    ierr = read_HEL(n_r_eq,eq_use_pol_flux)
                    CHCKERR('')
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end if
    end function read_eq
    
    ! prepare the cosine  and sine factors that are used  in the inverse Fourier
    ! transformation of R, Z and L and derivatives
    integer function prepare_RZL() result(ierr)
        use fourier_ops, only: calc_trigon_factors
        use eq_vars, only: trigon_factors, theta_E, zeta_E
        
        character(*), parameter :: rout_name = 'prepare_RZL'
        
        ! initialize ierr
        ierr = 0
        
        ierr = calc_trigon_factors(theta_E,zeta_E,trigon_factors)
        CHCKERR('')
    end function prepare_RZL
    
    ! calculate R, Z and Lambda and  derivatives in VMEC coordinates at the grid
    ! points given by the variables VMEC_theta and VMEC_zeta and at every normal
    ! point. The derivatives  are indicated by the variable "deriv"  which has 3
    ! indices
    integer function calc_RZL_ind(deriv) result(ierr)
        use fourier_ops, only: fourier2real
        use VMEC_ops, only: R_c, R_s, Z_c, Z_s, L_c, L_s
        use utilities, only: check_deriv
        use num_vars, only: max_deriv
        use eq_vars, only: VMEC_R, VMEC_Z, VMEC_L, grp_min_r_eq, grp_max_r_eq, &
            &trigon_factors
        
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
    ! correct poloidal HELENA range (and possibly multiples).
    integer function adapt_HEL_to_eq() result(ierr)
        use num_vars, only: pi
        use HEL_ops, only: h_H_11, h_H_12, h_H_33, ias, chi_H
        use utilities, only: interp_fun
        use eq_vars, only: grp_min_r_eq, n_par, grp_n_r_eq, theta_E
        
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
        deallocate(h_H_11); allocate(h_H_11(n_par,grp_n_r_eq))
        deallocate(h_H_12); allocate(h_H_12(n_par,grp_n_r_eq))
        deallocate(h_H_33); allocate(h_H_33(n_par,grp_n_r_eq))
        
        ! For every poloidal point, check  which half poloidal circle it belongs
        ! to.  If this  is a  bottom part  and HELENA  is symmetric  (ias =  0),
        ! the  quantities have  to  be taken  from  their symmetric  counterpart
        ! (2pi-theta) and  the metric factors  h_H_12 carry an  additional minus
        ! sign.
        ! Note that min_par and max_par have units of [pi]
        do id = 1,n_par
            ! loop over all normal points of this rank
            do kd = 1,grp_n_r_eq
                ! set the local poloidal point from theta_H
                par_loc = theta_E(id,kd)
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
                    ierr = interp_fun(h_H_11(id,kd),&
                        &old_h_H_11(:,grp_min_r_eq-1+kd),2*pi-par_loc,x=chi_H)
                    CHCKERR('')
                    ierr = interp_fun(h_H_12(id,kd),&
                        &-old_h_H_12(:,grp_min_r_eq-1+kd),2*pi-par_loc,x=chi_H) ! change of sign
                    CHCKERR('')
                    ierr = interp_fun(h_H_33(id,kd),&
                        &old_h_H_33(:,grp_min_r_eq-1+kd),2*pi-par_loc,x=chi_H)
                    CHCKERR('')
                else
                    ierr = interp_fun(h_H_11(id,kd),&
                        &old_h_H_11(:,grp_min_r_eq-1+kd),par_loc,x=chi_H)
                    CHCKERR('')
                    ierr = interp_fun(h_H_12(id,kd),&
                        &old_h_H_12(:,grp_min_r_eq-1+kd),par_loc,x=chi_H)
                    CHCKERR('')
                    ierr = interp_fun(h_H_33(id,kd),&
                        &old_h_H_33(:,grp_min_r_eq-1+kd),par_loc,x=chi_H)
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
        use eq_vars, only: flux_p_E, flux_t_E, pres_E, q_saf_E, rot_t_E, &
            &flux_t_E_full, flux_p_E_full, q_saf_E_full, rot_t_E_full, &
            &grp_min_r_eq, grp_max_r_eq, max_flux, max_flux_eq, max_flux_F, &
            &max_flux_eq_F, n_r_eq, eq_use_pol_flux
        
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
        
        call lvl_ud(1)
        
        ! plot flux quantities if requested
        if (plot_flux_q .and. glb_rank.eq.0) then
            ierr = flux_q_plot()
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        
        call lvl_ud(-1)
    contains
        ! VMEC version
        integer function calc_flux_q_VMEC() result(ierr)
            use VMEC_ops, only: iotaf, phi, phi_r, presf
            
            character(*), parameter :: rout_name = 'calc_flux_q_VMEC'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_p_full(:)                            ! version of dflux_p/dp on full normal grid (1..n_r_eq)
            real(dp), allocatable :: flux_p_int_full(:)                         ! version of integrated flux_p on full normal grid (1..n_r_eq)
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate poloidal flux
            allocate(Dflux_p_full(n_r_eq),flux_p_int_full(n_r_eq))
            Dflux_p_full = iotaf*phi_r
            ierr = calc_int(Dflux_p_full,1.0_dp/(n_r_eq-1.0_dp),flux_p_int_full)
            CHCKERR('')
            
            ! poloidal flux: calculate using iotaf and phi, phi_r
            ! easier to use full normal grid flux_p because of the integral
            flux_p_E(:,1) = Dflux_p_full(grp_min_r_eq:grp_max_r_eq)
            flux_p_E(:,0) = flux_p_int_full(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(flux_p_E(:,1),flux_p_E(:,kd),&
                    &n_r_eq-1._dp,kd-1,1)
                CHCKERR('')
            end do
                
            ! toroidal flux: copy from VMEC and derive
            flux_t_E(:,0) = phi(grp_min_r_eq:grp_max_r_eq)
            flux_t_E(:,1) = phi_r(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(flux_t_E(:,1),flux_t_E(:,kd),n_r_eq-1._dp,&
                    &kd-1,1)
                CHCKERR('')
            end do
            
            ! pressure: copy from VMEC and derive
            pres_E(:,0) = presf(grp_min_r_eq:grp_max_r_eq)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(pres_E(:,0),pres_E(:,kd),n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
            
            if (use_pol_flux) then
                ! safety factor
                q_saf_E(:,0) = 1.0_dp/iotaf(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(q_saf_E(:,0),q_saf_E(:,kd),n_r_eq-1._dp,&
                        &kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_p_int_full(n_r_eq)
                max_flux_F = max_flux
            else
                ! rot. transform
                rot_t_E(:,0) = iotaf(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(rot_t_E(:,0),rot_t_E(:,kd),n_r_eq-1._dp,&
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
                
            ! the global master needs flux_p_E_full, flux_t_E_full, q_saf_E_full
            ! and rot_t_E_full  for the resonance plot  and checking of m  and n
            ! and  the  group masters  need  q_saf_E  on  full normal  grid  for
            ! plot_X_vec
            if (grp_rank.eq.0) then
                ! flux_t_E_full
                flux_t_E_full(:,0) = phi
                flux_t_E_full(:,1) = phi_r
                do kd = 2,max_deriv+1
                    ierr = calc_deriv(flux_t_E_full(:,1),flux_t_E_full(:,kd),&
                        &n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! flux_p_E_full
                flux_p_E_full(:,0) = flux_p_int_full
                flux_p_E_full(:,1) = Dflux_p_full
                do kd = 2,max_deriv+1
                    ierr = calc_deriv(flux_p_E_full(:,1),&
                        &flux_p_E_full(:,kd),n_r_eq-1._dp,kd-1,1)
                    CHCKERR('')
                end do
                
                ! q_saf_E_full
                q_saf_E_full(:,0) = 1.0_dp/iotaf
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(q_saf_E_full(:,0),&
                        &q_saf_E_full(:,kd),n_r_eq-1._dp,kd,1)
                    CHCKERR('')
                end do
                
                ! rot_t_E_full
                rot_t_E_full(:,0) = iotaf
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(rot_t_E_full(:,0),&
                        &rot_t_E_full(:,kd),n_r_eq-1._dp,kd,1)
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
            use HEL_ops, only: qs, flux_H, p0
            
            character(*), parameter :: rout_name = 'calc_flux_q_HEL'
            
            ! local variables
            integer :: kd                                                       ! counter
            real(dp), allocatable :: Dflux_t_full(:)                            ! version of dflux_t/dp on full normal grid (1..n_r_eq)
            real(dp), allocatable :: flux_t_int_full(:)                         ! version of integrated flux_t on full normal grid (1..n_r_eq)
            real(dp), allocatable :: flux_H_r(:)                                ! normal derivative of flux_H
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables to calculate toroidal flux
            ! calculate normal derivative of flux_H
            allocate(flux_H_r(n_r_eq))
            ierr = calc_deriv(flux_H,flux_H_r,flux_H,1,1)
            CHCKERR('')
            allocate(Dflux_t_full(n_r_eq),flux_t_int_full(n_r_eq))
            Dflux_t_full = qs*flux_H_r
            ierr = calc_int(Dflux_t_full,flux_H,flux_t_int_full)
            CHCKERR('')
                
            ! poloidal flux: copy from HELENA and derive
            flux_p_E(:,0) = flux_H(grp_min_r_eq:grp_max_r_eq)
            do kd = 1,max_deriv+1
                ierr = calc_deriv(flux_p_E(:,0),flux_p_E(:,kd),&
                    &flux_p_E(:,0),kd,1)
                CHCKERR('')
            end do
            
            ! toroidal flux: calculate using qs and flux_H, flux_H_r
            ! easier to use full normal grid flux_t because of the integral
            flux_t_E(:,1) = Dflux_t_full(grp_min_r_eq:grp_max_r_eq)
            flux_t_E(:,0) = flux_t_int_full(grp_min_r_eq:grp_max_r_eq)
            do kd = 2,max_deriv+1
                ierr = calc_deriv(flux_t_E(:,1),flux_t_E(:,kd),&
                    &flux_p_E(:,0),kd-1,1)
                CHCKERR('')
            end do
            
            ! pressure: copy from HELENA and derive
            pres_E(:,0) = p0(grp_min_r_eq:grp_max_r_eq)
            do kd = 1, max_deriv+1
                ierr = calc_deriv(pres_E(:,0),pres_E(:,kd),&
                    &flux_p_E(:,0),kd,1)
                CHCKERR('')
            end do
            
            if (use_pol_flux) then
                ! safety factor
                q_saf_E(:,0) = qs(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(q_saf_E(:,0),q_saf_E(:,kd),&
                        &flux_p_E(:,0),kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_H(n_r_eq)
            else
                ! rot. transform
                rot_t_E(:,0) = 1.0_dp/qs(grp_min_r_eq:grp_max_r_eq)
                do kd = 1,max_deriv+1
                    ierr = calc_deriv(rot_t_E(:,0),rot_t_E(:,kd),&
                        &flux_p_E(:,0),kd,1)
                    CHCKERR('')
                end do
                
                ! set max_flux
                max_flux = flux_t_int_full(n_r_eq)
            end if
            max_flux_F = max_flux
            
            ! max_flux_eq
            max_flux_eq = flux_H(n_r_eq)
            max_flux_eq_F = max_flux_eq
            
            ! the global master needs flux_p_E_full, flux_t_E_full, q_saf_E_full
            ! and rot_t_E_full  for the resonance plot  and checking of m  and n
            ! and  the  group masters  need  q_saf_E  on  full normal  grid  for
            ! plot_X_vec
            if (grp_rank.eq.0) then
                ! flux_p_E_full
                flux_p_E_full(:,0) = flux_H
                do kd = 1,max_deriv
                    ierr = calc_deriv(flux_p_E_full(:,0),flux_p_E_full(:,kd),&
                        &flux_p_E_full(:,0),kd,1)
                    CHCKERR('')
                end do
                
                ! flux_t_E_full
                flux_t_E_full(:,0) = flux_t_int_full
                flux_t_E_full(:,1) = Dflux_t_full
                do kd = 2,max_deriv
                    ierr = calc_deriv(flux_t_E_full(:,1),flux_t_E_full(:,kd),&
                        &flux_p_E_full(:,0),kd-1,1)
                    CHCKERR('')
                end do
                
                ! q_saf_E_full
                q_saf_E_full(:,0) = qs(:)
                do kd = 1,max_deriv
                    ierr = calc_deriv(q_saf_E_full(:,0),q_saf_E_full(:,kd),&
                        &flux_p_E_full(:,0),kd,1)
                    CHCKERR('')
                end do
                
                ! rot_t_E_full
                rot_t_E_full(:,0) = 1.0_dp/qs(:)
                do kd = 1,max_deriv
                    ierr = calc_deriv(rot_t_E_full(:,0),rot_t_E_full(:,kd),&
                        &flux_p_E_full(:,0),kd,1)
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
        use VMEC_ops, only: presf
        use HEL_ops, only: p0
        use eq_vars, only: q_saf_E_full, rot_t_E_full, flux_p_E_full, &
            &flux_t_E_full, max_flux, max_flux_eq, n_r_eq
        
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
                y_plot_2D(:,1) = -q_saf_E_full(:,0)                             ! conversion VMEC LH -> RH coord. system
                y_plot_2D(:,2) = -rot_t_E_full(:,0)                             ! conversion VMEC LH -> RH coord. system
                y_plot_2D(:,3) = presf
                y_plot_2D(:,4) = flux_p_E_full(:,0)
                y_plot_2D(:,5) = -flux_t_E_full(:,0)                            ! conversion VMEC LH -> RH coord. system
                ! 2D normal variable (y_plot_2D tabulated in eq. grid)
                if (use_pol_flux) then
                    x_plot_2D(:,1) = flux_p_E_full(:,0)/max_flux
                else
                    x_plot_2D(:,1) = flux_t_E_full(:,0)/max_flux
                end if
                do id = 2,n_vars
                    x_plot_2D(:,id) = x_plot_2D(:,1)
                end do
            case (2)                                                            ! HELENA
                y_plot_2D(:,1) = q_saf_E_full(:,0)
                y_plot_2D(:,2) = rot_t_E_full(:,0)
                y_plot_2D(:,3) = p0
                y_plot_2D(:,4) = flux_p_E_full(:,0)
                y_plot_2D(:,5) = flux_t_E_full(:,0)
                if (use_pol_flux) then
                    x_plot_2D(:,1) = flux_p_E_full(:,0)/max_flux
                else
                    x_plot_2D(:,1) = flux_t_E_full(:,0)/max_flux
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
                        r_plot = flux_p_E_full(:,0)/max_flux_eq                 ! output from HELENA
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
            use coord_ops, only: calc_XYZ_grid
            
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
            ierr = calc_XYZ_grid(r_plot,theta_plot,zeta_plot,&
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
    
    ! normalizes equilibrium quantities pres_FD, q_saf_FD or rot_t_FD, flux_p_FD
    ! or flux_t_FD, max_flux and pres using the normalization constants
    subroutine normalize_eq_vars
        use num_vars, only: use_pol_flux
        use eq_vars, only: pres_FD, flux_p_FD, flux_t_FD, q_saf_FD, rot_t_FD, &
            &max_flux, max_flux_eq, max_flux_F, max_flux_eq_F, pres_0, psi_0
        
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
        use eq_vars, only: T_0, B_0, pres_0, psi_0, R_0, rho_0
        
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
            use VMEC_ops, only: R_c, presf
            
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
            use HEL_ops, only: R_0_H, B_0_H
            
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
end module eq_ops
