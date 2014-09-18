!------------------------------------------------------------------------------!
!   Variables,  subroutines  and functions  that  have  to do  with            !
!   equilibrium quantities and the mesh used in the calculations               !
!   These are called in the subroutine calc_eq in the module eq_ops            !
!------------------------------------------------------------------------------!
module eq_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, pi, max_str_ln
    use output_ops, only: writo, lvl_ud, print_ar_2, print_ar_1,  print_GP_2D
    use str_ops, only: r2strt, i2str, r2str

    implicit none
    private
    public calc_eqd_mesh, calc_mesh, calc_ang_B, calc_flux_q, prepare_RZL, &
        &check_and_limit_mesh, init_eq, calc_RZL, dealloc_eq, &
        &dealloc_eq_final, normalize_eq_vars, calc_norm_const, &
        &theta_V, zeta_V, n_par, grp_n_r_eq, lam_H, min_par, max_par, &
        &q_saf_V, q_saf_V_full, flux_p_V, flux_t_V, VMEC_R, VMEC_Z, VMEC_L, &
        &pres, q_saf_FD, flux_p_FD, flux_t_FD, pres_FD, grp_min_r_eq, &
        &grp_max_r_eq, R_0, A_0, pres_0, B_0, psi_p_0, rho_0, ang_par_F, &
        &rot_t_FD, rot_t_V, flux_p_V_full, flux_t_V_full, rot_t_V_full, &
        &max_flux, max_flux_VMEC

    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r deriv., 2: theta_V deriv., 3: zeta_V deriv.)
    real(dp), allocatable :: trigon_factors(:,:,:,:,:)                          ! trigonometric factor cosine for the inverse fourier transf.
    real(dp), allocatable :: VMEC_R(:,:,:,:,:)                                  ! R in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: VMEC_Z(:,:,:,:,:)                                  ! Z in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: VMEC_L(:,:,:,:,:)                                  ! L(ambda) in VMEC coordinates (n_par, grp_n_r_eq, [derivatives])
    real(dp), allocatable :: lam_H(:,:,:)                                       ! lambda in (HM)
    real(dp), allocatable :: theta_V(:,:), zeta_V(:,:)                          ! grid points (n_par, grp_n_r_eq, 2) in VMEC coords
    real(dp), allocatable :: ang_par_F(:,:)                                     ! parallel angle in flux coordinates (either zeta_F or theta_F)
    real(dp), allocatable :: q_saf_V(:,:)                                       ! safety factor in VMEC coordinates
    real(dp), allocatable :: q_saf_V_full(:,:)                                  ! safety factor in full normal mesh in flux coordinates
    real(dp), allocatable :: q_saf_FD(:,:)                                      ! safety factor, Deriv. in Flux coords.
    real(dp), allocatable :: rot_t_V(:,:)                                       ! rot. transform in VMEC coordinates
    real(dp), allocatable :: rot_t_V_full(:,:)                                  ! rot. transform in full normal mesh in flux coordinates
    real(dp), allocatable :: rot_t_FD(:,:)                                      ! rot. transform, Deriv. in Flux coords.
    real(dp), allocatable :: flux_p_V(:,:), flux_t_V(:,:), pres(:,:)            ! pol. flux, tor. flux and pressure, and norm. Deriv. in VMEC coords.
    real(dp), allocatable :: pres_FD(:,:)                                       ! pressure, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable, target :: flux_p_FD(:,:), flux_t_FD(:,:)             ! pol. and tor. flux, and norm. Deriv. with values and Derivs. in flux coords.
    real(dp), allocatable :: flux_p_V_full(:,:), flux_t_V_full(:,:)             ! pol. flux, tor. flux, and norm. Deriv. values and Derivs. in VMEC coords.
    real(dp) :: max_flux, max_flux_VMEC                                         ! max. flux (pol. or tor.) (min.flux is trivially equal to 0)
    real(dp) :: min_par, max_par                                                ! min. and max. of parallel coordinate
    real(dp) :: R_0, A_0, pres_0, B_0, psi_p_0, rho_0                           ! normalization constants for nondimensionalization
    integer :: n_par, grp_n_r_eq                                                ! nr. of parallel and normal points in this process in alpha group
    integer :: grp_min_r_eq, grp_max_r_eq                                       ! min. and max. r range of this process in alpha group
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! initialize the equilibrium variables
    subroutine init_eq
        use num_vars, only: max_deriv, grp_rank, use_pol_flux
        use VMEC_vars, only: n_r_eq
        
        ! calculate grp_n_r_eq
        grp_n_r_eq = grp_max_r_eq - grp_min_r_eq + 1
        
        ! R
        allocate(VMEC_R(n_par,grp_n_r_eq,0:max_deriv(1),0:max_deriv(2),&
            &0:max_deriv(3)))
        
        ! Z
        allocate(VMEC_Z(n_par,grp_n_r_eq,0:max_deriv(1),0:max_deriv(2),&
            &0:max_deriv(3)))
        
        ! lambda
        allocate(VMEC_L(n_par,grp_n_r_eq,0:max_deriv(1),0:max_deriv(2),&
            &0:max_deriv(3)))
        
        ! pres
        allocate(pres(grp_n_r_eq,0:max_deriv(1)))
        
        ! pres_FD
        allocate(pres_FD(grp_n_r_eq,0:max_deriv(1)))
        
        ! flux_p_V
        allocate(flux_p_V(grp_n_r_eq,0:max_deriv(1)))
        
        ! flux_p_FD
        allocate(flux_p_FD(grp_n_r_eq,0:max_deriv(1)))
            
        ! flux_t_V
        allocate(flux_t_V(grp_n_r_eq,0:max_deriv(1)))
        
        ! flux_t_FD
        allocate(flux_t_FD(grp_n_r_eq,0:max_deriv(1)))
        
        if (use_pol_flux) then
            ! q_saf_V
            allocate(q_saf_V(grp_n_r_eq,0:max_deriv(1)))
            
            ! q_saf_FD
            allocate(q_saf_FD(grp_n_r_eq,0:max_deriv(1)))
        else
            ! rot_t_V
            allocate(rot_t_V(grp_n_r_eq,0:max_deriv(1)))
            
            ! rot_t_FD
            allocate(rot_t_FD(grp_n_r_eq,0:max_deriv(1)))
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
    end subroutine init_eq
    
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
        use fourier_ops, only: calc_trigon_factors, fourier2real
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
        use VMEC_vars, only: iotaf, phi, phi_r, presf, n_r_eq, VMEC_use_pol_flux
        use utilities, only: VMEC_norm_deriv, calc_int
        use num_vars, only: max_deriv, grp_rank, use_pol_flux
        
        character(*), parameter :: rout_name = 'calc_flux_q'
        
        ! local variables
        integer :: kd                                                           ! counter
        real(dp), allocatable :: Dflux_p_full(:)                                ! version of flux_dp on full normal mesh (1..n_r_eq)
        real(dp), allocatable :: flux_p_int_full(:)                             ! version of integrated flux_dp on full normal mesh (1..n_r_eq)
        
        ! initialize ierr
        ierr = 0
        
        ! pressure: copy from VMEC and derive
        pres(:,0) = presf(grp_min_r_eq:grp_max_r_eq)
        do kd = 1, max_deriv(1)
            ierr = VMEC_norm_deriv(pres(:,0),pres(:,kd),n_r_eq-1._dp,kd,1)
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
            ierr = VMEC_norm_deriv(flux_p_V(:,1),flux_p_V(:,kd),&
                &n_r_eq-1._dp,kd-1,1)
            CHCKERR('')
        end do
            
        ! toroidal flux: copy from VMEC and derive
        flux_t_V(:,0) = phi(grp_min_r_eq:grp_max_r_eq)
        flux_t_V(:,1) = phi_r(grp_min_r_eq:grp_max_r_eq)
        do kd = 2,max_deriv(1)
            ierr = VMEC_norm_deriv(flux_t_V(:,1),flux_t_V(:,kd),n_r_eq-1._dp,&
                &kd-1,1)
            CHCKERR('')
        end do
        
        if (use_pol_flux) then
            ! safety factor
            q_saf_V(:,0) = 1.0_dp/iotaf(grp_min_r_eq:grp_max_r_eq)
            do kd = 1,max_deriv(1)
                ierr = VMEC_norm_deriv(q_saf_V(:,0),q_saf_V(:,kd),n_r_eq-1._dp,&
                    &kd,1)
                CHCKERR('')
            end do
            
            ! set max_flux
            max_flux = flux_p_int_full(n_r_eq)
        else
            ! rot. transform
            rot_t_V(:,0) = iotaf(grp_min_r_eq:grp_max_r_eq)
            do kd = 1,max_deriv(1)
                ierr = VMEC_norm_deriv(rot_t_V(:,0),rot_t_V(:,kd),n_r_eq-1._dp,&
                    &kd,1)
                CHCKERR('')
            end do
            
            ! set max_flux
            max_flux = phi(n_r_eq)
        end if
        
        ! max_flux_VMEC
        if (VMEC_use_pol_flux) then
            max_flux_VMEC = flux_p_int_full(n_r_eq)
        else
            max_flux_VMEC = phi(n_r_eq)
        end if
            
        ! the global master needs flux_p_V_full, flux_t_V_full, q_saf_V_full and
        ! rot_t_V_full for  the resonance plot and  checking of m and  n and the
        ! group masters need q_saf_V on full normal mesh for plot_X_vec
        if (grp_rank.eq.0) then
            ! flux_t_V_full
            flux_t_V_full(:,0) = phi
            flux_t_V_full(:,1) = phi_r
            do kd = 2,max_deriv(1)
                ierr = VMEC_norm_deriv(flux_t_V_full(:,1),flux_t_V_full(:,kd),&
                    &n_r_eq-1._dp,kd-1,1)
                CHCKERR('')
            end do
            
            ! flux_p_V_full
            flux_p_V_full(:,0) = flux_p_int_full
            flux_p_V_full(:,1) = Dflux_p_full
            do kd = 2,max_deriv(1)
                ierr = VMEC_norm_deriv(flux_p_V_full(:,1),&
                    &flux_p_V_full(:,kd),n_r_eq-1._dp,kd-1,1)
                CHCKERR('')
            end do
            
            ! q_saf_V_full
            q_saf_V_full(:,0) = 1.0_dp/iotaf
            do kd = 1,max_deriv(1)
                ierr = VMEC_norm_deriv(q_saf_V_full(:,0),&
                    &q_saf_V_full(:,kd),n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
            
            ! rot_t_V_full
            rot_t_V_full(:,0) = iotaf
            do kd = 1,max_deriv(1)
                ierr = VMEC_norm_deriv(rot_t_V_full(:,0),&
                    &rot_t_V_full(:,kd),n_r_eq-1._dp,kd,1)
                CHCKERR('')
            end do
        end if
        
        ! deallocate helper variables
        if (use_pol_flux .or. grp_rank.eq.0) then
            deallocate(Dflux_p_full,flux_p_int_full)
        end if
    end function calc_flux_q

    ! calculate   the    angular   mesh,   filling   the    global   variables. 
    ! (theta_V,zeta_V) and ang_par_F, for every flux surface      . 
    ! The  variable  use_pol_flux  determines  whether theta  (.true.)  or  zeta
    ! (.false.) is used as the parallel variable
    ! for  the  normal  mesh  type  (following  the  magnetic  field),  theta_V,
    ! zeta_V and  ang_par_F are  defined on  the entire  normal mesh.  Later, in
    ! check_and_limit_mesh,  the  normal extent  is  limited  if the  tests  are
    ! positive
    ! the global  variable calc_mesh_style can  optionally force other  types of
    ! meshes, which is useful for tests
    integer function calc_mesh(alpha) result(ierr)
        use num_vars, only: calc_mesh_style, use_pol_flux
        use VMEC_vars, only: n_r_eq, iotaf
        
        character(*), parameter :: rout_name = 'calc_mesh'
        
        ! input / output
        real(dp), intent(in) :: alpha
        
        ! local variables
        integer :: kd, id                                                       ! counters
        real(dp), allocatable :: var_in(:)                                      ! input of calc_ang_B
        real(dp), allocatable :: var_out(:)                                     ! output of calc_ang_B
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize variables
        allocate(zeta_V(n_par,n_r_eq)); zeta_V = 0.0_dp
        allocate(theta_V(n_par,n_r_eq)); theta_V = 0.0_dp
        allocate(ang_par_F(n_par,n_r_eq)); ang_par_F = 0.0_dp
        
        allocate(var_in(n_r_eq))
        allocate(var_out(n_r_eq))
        
        ! calculate the mesh
        select case (calc_mesh_style)
            ! grid along the magnetic field lines
            case (0)
                ! set up parallel angle in flux coordinates on equidistant mesh
                do kd = 1,n_r_eq
                    ierr = calc_eqd_mesh(ang_par_F(:,kd),n_par,min_par,max_par)
                    CHCKERR('')
                end do
                ! calculate zeta_V from alpha and ang_par_F
                if (use_pol_flux) then                                          ! theta_F is parallel coordinate
                    do id = 1,n_par
                        zeta_V(id,:) = ang_par_F(id,:)/iotaf
                    end do
                    zeta_V = zeta_V - alpha
                else                                                            ! zeta_F is parallel coordinate
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
                    ierr = calc_eqd_mesh(zeta_V(:,kd),n_par, min_par, max_par)
                    CHCKERR('')
            end do
                call writo('Defining grid with theta_V = '//&
                    &trim(r2strt(theta_V(1,1)))//&
                    &' and zeta_V equidistant from '//&
                    &trim(r2strt(min_par))//' to '//trim(r2strt(max_par)))
            ! grid with constant zeta_V, equidistant theta_V
            case (2)
                zeta_V = alpha
                do kd = 1,n_r_eq
                    ierr = calc_eqd_mesh(theta_V(:,kd),n_par,min_par,max_par)
                    CHCKERR('')
                end do
                call writo('Defining grid with zeta_V = '//&
                    &trim(r2strt(zeta_V(1,1)))//&
                    &' and theta_V equidistant from '//&
                    &trim(r2strt(min_par))//' to '//trim(r2strt(max_par)))
            ! grid style error
            case default
                err_msg = 'No style is associated with '&
                    &//trim(i2str(calc_mesh_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function calc_mesh
    
    ! check  whether   the  straight  field   line  mesh  has   been  calculated
    ! satisfactorily,  by calculating  again the  theta_V's from  the calculated
    ! zeta_V's
    integer function check_and_limit_mesh(alpha) result(ierr)
        use num_vars, only: tol_NR, max_str_ln, calc_mesh_style, use_pol_flux
        use VMEC_vars, only: n_r_eq, VMEC_use_pol_flux
        
        character(*), parameter :: rout_name = 'check_and_limit_mesh'
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        real(dp) :: var_calc(n_par,n_r_eq)
        real(dp) :: var_diff(n_par,n_r_eq)
        real(dp) :: var_in(n_r_eq)
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: work(:)                                        ! work array
        character(len=5) :: par_c                                               ! name of parallel coord.
        character(len=8) :: nor_c                                               ! name of normal coord.
        
        ! initialize ierr
        ierr = 0
        
        if (calc_mesh_style.eq.0) then                                          ! only do the check if the magnetic field lines are followed
            ! allocate work variable
            allocate(work(n_r_eq))
            
            do id = 1,n_par
                var_in = theta_V(id,:)
                ierr = calc_ang_B(work,.false.,alpha,var_in)
                var_calc(id,:) =  work
                CHCKERR('')
            end do
            var_diff = zeta_V(:,:) - var_calc
            
            ! check the difference
            if (maxval(abs(var_diff)).gt.tol_NR*100) then                       ! difference too large
                err_msg = 'Calculating again the toroidal grid points from &
                    &the poloidal grid points, along the magnetic fields, &
                    &does not result in the same values as initally given... &
                    &The maximum error is equal to '//&
                    &trim(r2strt(100*maxval(abs(var_diff))))//'%'
                ierr = 1
                CHCKERR(err_msg)
            else                                                                ! difference within tolerance
                ! limit the normal size of theta_V and zeta_V
                call limit_normal_size(theta_V)
                call limit_normal_size(zeta_V)
                call limit_normal_size(ang_par_F)
            end if
            
            ! deallocate work variable
            deallocate(work)
            
            if (use_pol_flux) then
                par_c = 'theta'
            else
                par_c = 'zeta'
            end if
            if (VMEC_use_pol_flux) then
                nor_c = 'poloidal'
            else
                nor_c = 'toroidal'
            end if
            call writo('angular grid set up with coordinate '//trim(par_c)//&
                &' oriented along the magnetic field')
            call writo('normal grid set with '//trim(nor_c)//' flux as normal &
                &coordinate')
            if (use_pol_flux.neqv.VMEC_use_pol_flux) write(*,*) &
                &'!!!!!!!!!!!!!!! ADD TEST TO CHECK WHETHER WE CAN INDEED &
                &USE DIFFERENT FLUX THAN VMEC!!!!!!!!!!!'
        else
            ! limit the normal size of theta_V and zeta_V
            call limit_normal_size(theta_V)
            call limit_normal_size(zeta_V)
            call limit_normal_size(ang_par_F)
        end if
    contains
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
    integer function calc_eqd_mesh(eqd_mesh,n_ang, min_ang, max_ang) result(ierr)
        character(*), parameter :: rout_name = 'eqd_mesh'
        
        ! input and output
        real(dp), intent(inout) :: eqd_mesh(:)                                  ! output
        real(dp), intent(in) :: min_ang, max_ang                                ! min. and max. of angles
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
        
        ! There are (n_ang-1) pieces in the total interval but the last one 
        ! is not needed as the functions are all periodic
        
        delta_ang = (max_ang-min_ang)/(n_ang)
        
        eqd_mesh(1) = min_ang
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
        use VMEC_vars, only: L_c, L_s, iotaf, n_r_eq
        use fourier_ops, only: calc_trigon_factors, fourier2real
        use utilities, only: calc_zero_NR, VMEC_conv_FHM
        
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
        use VMEC_vars, only: VMEC_use_pol_flux
        
        ! local variables
        integer :: id                                                           ! counter
        real(dp) :: psi_0                                                       ! either psi_p_0 (pol. flux) or psi_p_0*A_0 (tor. flux)
        
        ! scale the quantities
        pres_FD = pres_FD/pres_0
        flux_p_FD = flux_p_FD/psi_p_0
        flux_t_FD = flux_t_FD/(psi_p_0*A_0)
        if (use_pol_flux) then
            q_saf_FD = q_saf_FD/A_0
            max_flux = max_flux/psi_p_0
        else
            rot_t_FD = rot_t_FD 
            max_flux = max_flux/(psi_p_0*A_0)
        end if
        if (VMEC_use_pol_flux) then
            max_flux_VMEC = max_flux_VMEC/psi_p_0
        else
            max_flux_VMEC = max_flux_VMEC/(psi_p_0*A_0)
        end if
        
        ! set up psi_0
        if (use_pol_flux) then
            psi_0 = psi_p_0
        else
            psi_0 = psi_p_0 * A_0
        end if
        
        ! scale  the  derivatives  by  psi_p_0  if poloidal  flux  is  used,  or
        ! psi_p_0*A_0 if toroidal flux is used
        if (use_pol_flux) then
            do id = 1,size(pres_FD,2)-1
                pres_FD(:,id) = pres_FD(:,id) * psi_0**(id)
                q_saf_FD(:,id) = q_saf_FD(:,id) * psi_0**(id)
                flux_p_FD(:,id) = flux_p_FD(:,id) * psi_0**(id)
            end do
        else
            do id = 1,size(pres_FD,2)-1
                pres_FD(:,id) = pres_FD(:,id) * psi_0**(id)
                rot_t_FD(:,id) = rot_t_FD(:,id) * psi_0**(id)
                flux_t_FD(:,id) = flux_t_FD(:,id) * psi_0**(id)
            end do
        end if
    end subroutine normalize_eq_vars
    
    ! sets up normalization constants:
    !   R_0:     major radius (= average R on axis)
    !   A_0:     major / minor radius
    !   pres_0:  pressure on axis
    !   B_0:     average poloidal field on axis
    !   psi_p_0:   reference poloidal flux
    ! [MPI] only global master
    subroutine calc_norm_const
        use num_vars, only: mu_0, glb_rank
        use VMEC_vars, only: R_c, zmax_surf, presf
        
        if (glb_rank.eq.0) then
            ! calculate  the  major radius  as the  average value  of VMEC_R  on the
            ! magnetic axis
            R_0 = R_c(0,0,1,0)
            
            ! calculate  the aspect ratio as  the major radius divided  by the minor
            ! radius
            A_0 = R_0 / zmax_surf
            
            ! calculate pres_0 as pressure on axis
            pres_0 = presf(1)
            
            ! calculate the reference value for B_0 from B_0 = sqrt(mu_0 pres_0)
            B_0 = sqrt(pres_0 * mu_0)
            
            ! calculate reference pol. flux
            ! (reference tor. flux is multiplied by A_0)
            psi_p_0 = (R_0**2 * B_0) / A_0
        end if
    end subroutine calc_norm_const
    
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
        
        deallocate(pres,pres_FD)
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
end module eq_vars
