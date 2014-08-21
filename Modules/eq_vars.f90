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
    public calc_eqd_mesh, calc_mesh, calc_ang_B, calc_flux_q, check_mesh, &
        &init_eq, calc_RZL, dealloc_eq_vars, dealloc_eq_vars_final, &
        &theta, zeta, n_par, n_r, R, Z, lam_H, min_par, max_par, &
        &q_saf, flux_p, flux_t, VMEC_R, VMEC_Z, VMEC_L, pres, &
        &q_saf_FD, flux_p_FD, flux_t_FD, pres_FD

    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative)
    real(dp), allocatable :: R(:,:,:), Z(:,:,:)                                 ! R, Z (FM)
    real(dp), allocatable :: VMEC_R(:,:,:,:,:)                                  ! R in VMEC coordinates (n_par, n_r, [derivatives])
    real(dp), allocatable :: VMEC_Z(:,:,:,:,:)                                  ! Z in VMEC coordinates (n_par, n_r, [derivatives])
    real(dp), allocatable :: VMEC_L(:,:,:,:,:)                                  ! L(ambda) in VMEC coordinates (n_par, n_r, [derivatives])
    real(dp), allocatable :: lam_H(:,:,:)                                       ! lambda in (HM)
    real(dp), allocatable :: theta(:,:), zeta(:,:)                              ! grid points (n_par, n_r, 2) (FM)
    real(dp), allocatable :: q_saf(:,:)                                         ! safety factor (FM)
    real(dp), allocatable :: q_saf_FD(:,:)                                      ! safety factor (FM), derivatives in Flux coords.
    real(dp), allocatable :: flux_p(:,:), flux_t(:,:), pres(:,:)                ! pol. flux, tor. flux and pressure, and normal derivative (FM)
    real(dp), allocatable :: flux_p_FD(:,:), flux_t_FD(:,:), pres_FD(:,:)       ! pol. flux, tor. flux and pressure, and normal derivative (FM), derivs. in flux coords.
    real(dp) :: min_par, max_par
    integer :: n_par, n_r
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! initialize the variables n_r, VMEC_R, VMEC_Z and VMEC_L
    subroutine init_eq
        use num_vars, only: max_deriv
        use VMEC_vars, only: ns
        
        ! n_r
        n_r = ns
        
        ! R
        allocate(VMEC_R(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! Z
        allocate(VMEC_Z(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! lambda
        allocate(VMEC_L(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! q_saf
        allocate(q_saf(n_r,0:max_deriv(1)))
        
        ! q_saf_FD
        allocate(q_saf_FD(n_r,0:max_deriv(1)))
        
        ! flux_p
        allocate(flux_p(n_r,0:max_deriv(1)))
        
        ! flux_p_FD
        allocate(flux_p_FD(n_r,0:max_deriv(1)))
        
        ! flux_t
        allocate(flux_t(n_r,0:max_deriv(1)))
        
        ! flux_t_FD
        allocate(flux_t_FD(n_r,0:max_deriv(1)))
        
        ! pres
        allocate(pres(n_r,0:max_deriv(1)))
        
        ! pres_FD
        allocate(pres_FD(n_r,0:max_deriv(1)))
    end subroutine
    
    ! calculate R, Z and Lambda and  derivatives in VMEC coordinates at the mesh
    ! points given by the variables VMEC_theta and VMEC_zeta and at every normal
    ! point. The derivatives  are indicated by the variable "deriv"  which has 3
    ! indices
    integer function calc_RZL_ind(deriv) result(ierr)
        use fourier_ops, only: calc_mesh_cs, f2r
        use VMEC_vars, only: mpol, ntor, nfp, R_c, R_s, Z_c, Z_s, L_c, L_s
        use utilities, only: check_deriv
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_RZL_ind'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        
        ! initialize ierr
        ierr = 0
        
        ! check the derivatives requested
        ierr = check_deriv(deriv,max_deriv,'calc_RZL')
        CHCKERR('')
        
        ! calculate the  variables R, Z,  and their angular derivatives  for all
        ! normal points and current angular point
        perp: do kd = 1, n_r                                                    ! perpendicular: normal to the flux surfaces
            par: do id = 1, n_par                                               ! parallel: along the magnetic field line
                ierr = calc_mesh_cs(cs,mpol,ntor,nfp,theta(id,kd),zeta(id,kd))
                CHCKERR('')
                VMEC_R(id,kd,deriv(1),deriv(2),deriv(3)) = &
                    &f2r(R_c(:,:,kd,deriv(1)),R_s(:,:,kd,deriv(1)),cs&
                    &,mpol,ntor,nfp,[deriv(2),deriv(3)],ierr)
                CHCKERR('')
                VMEC_Z(id,kd,deriv(1),deriv(2),deriv(3)) = &
                    &f2r(Z_c(:,:,kd,deriv(1)),Z_s(:,:,kd,deriv(1)),cs&
                    &,mpol,ntor,nfp,[deriv(2),deriv(3)],ierr)
                CHCKERR('')
                VMEC_L(id,kd,deriv(1),deriv(2),deriv(3)) = &
                    &f2r(L_c(:,:,kd,deriv(1)),L_s(:,:,kd,deriv(1)),cs&
                    &,mpol,ntor,nfp,[deriv(2),deriv(3)],ierr)
                CHCKERR('')
            end do par
        end do perp
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

    ! Calculates flux quantities and normal derivatives
    integer function calc_flux_q() result(ierr)
        use VMEC_vars, only: iotaf, phi, phi_r, presf
        use utilities, only: VMEC_norm_deriv, calc_int
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_flux_q'
        
        ! local variables
        integer :: kd
        
        ! initialize ierr
        ierr = 0
        
        ! safety factor q_saf: invert iota and derive
        do kd = 1,n_r
            q_saf(kd,0) = 1/iotaf(kd)
        end do
        do kd = 1,max_deriv(1)
            ierr = VMEC_norm_deriv(q_saf(:,0),q_saf(:,kd),n_r-1._dp,kd,1)
            CHCKERR('')
        end do
        
        ! toroidal flux: copy from VMEC and derive
        flux_t(:,0) = phi
        flux_t(:,1) = phi_r
        do kd = 2,max_deriv(1)
            ierr = VMEC_norm_deriv(flux_t(:,1),flux_t(:,kd),n_r-1._dp,kd-1,1)
            CHCKERR('')
        end do
        
        ! poloidal flux: calculate using iotaf and phi, phi_r
        flux_p = 0.0_dp
        flux_p(:,1) = iotaf*phi_r
        ierr = calc_int(flux_p(:,1),[(kd/(n_r-1.0_dp),kd=0,n_r-1)],flux_p(:,0))
        CHCKERR('')
        do kd = 2,max_deriv(1)
            ierr = VMEC_norm_deriv(flux_p(:,1),flux_p(:,kd),n_r-1._dp,kd-1,1)
            CHCKERR('')
        end do
        
        ! pressure: copy from VMEC and derive
        pres(:,0) = presf
        do kd = 1, max_deriv(1)
            ierr = VMEC_norm_deriv(pres(:,0),pres(:,kd),n_r-1._dp,kd,1)
            CHCKERR('')
        end do
        !call print_GP_2D('pres','',presf,draw=.true.)
        !do kd = 0,max_deriv(1)
            !call print_GP_2D('pres^'//trim(i2str(kd)),'',pres(:,kd),draw=.true.)
        !end do
    end function calc_flux_q

    ! calculate the angular mesh, filling the global variables (theta, zeta). 
    ! This is done for every flux surface. 
    ! The variable  theta_var_along_B determines  whether theta  is used  as the
    ! base variable. If .false., zeta is used.
    ! the global  variable calc_mesh_style can  optionally force other  types of
    ! meshes, which is useful for tests
    integer function calc_mesh(alpha) result(ierr)
        use num_vars, only: theta_var_along_B, calc_mesh_style
        
        character(*), parameter :: rout_name = 'calc_mesh'
        
        ! input / output
        real(dp), intent(in) :: alpha
        
        ! local variables
        integer :: jd, kd
        real(dp) :: var_in(n_r)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: work(:)                                        ! work array
        
        ! initialize ierr
        ierr = 0
        
        allocate(zeta(n_par,n_r)); zeta = 0.0_dp
        allocate(theta(n_par,n_r)); theta = 0.0_dp
        
        allocate(work(n_r))
        
        ! calculate the mesh
        select case (calc_mesh_style)
            case (0)                                                            ! along the magnetic field lines
                ierr = check_periodicity(min_par,max_par)
                CHCKERR('')
                
                if (theta_var_along_B) then                                     ! first calculate theta
                    do jd = 1,n_r
                        ierr = calc_eqd_mesh(theta(:,jd),n_par,min_par,max_par)
                        CHCKERR('')
                    end do
                    do kd = 1,n_par
                        var_in = theta(kd,:)
                        ierr = calc_ang_B(work,.false.,alpha,var_in)
                        zeta(kd,:) = work
                        CHCKERR('')
                    end do
                else                                                            ! first calculate zeta
                    do jd = 1,n_r
                        ierr = calc_eqd_mesh(zeta(:,jd),n_par,min_par,max_par)
                        CHCKERR('')
                    end do
                    do kd = 1,n_par
                        var_in = zeta(kd,:)
                        ierr = calc_ang_B(work,.true.,alpha,var_in)
                        theta(kd,:) = work
                        CHCKERR('')
                    end do
                end if
            case (1)                                                            ! constant theta, equidistant zeta
                theta = alpha
                do jd = 1,n_r
                    ierr = calc_eqd_mesh(zeta(:,jd),n_par, min_par, max_par)
                    CHCKERR('')
            end do
                call writo('Defining grid with theta = '//&
                    &trim(r2strt(theta(1,1)))//' and zeta equidistant from '//&
                    &trim(r2strt(min_par))//' to '//trim(r2strt(max_par)))
            case (2)                                                            ! constant zeta, equidistant theta
                zeta = alpha
                do jd = 1,n_r
                    ierr = calc_eqd_mesh(theta(:,jd),n_par,min_par,max_par)
                    CHCKERR('')
                end do
                call writo('Defining grid with zeta = '//&
                    &trim(r2strt(zeta(1,1)))//' and theta equidistant from '//&
                    &trim(r2strt(min_par))//' to '//trim(r2strt(max_par)))
            case default
                err_msg = 'No style is associated with '&
                    &//trim(i2str(calc_mesh_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! checks whether max_par - min_par indeed is a multiple of 2pi, which is
        ! a requisite for the correct functioning  of the code, as the integrals
        ! which are to be calculated in the parallel direction have to be closed
        ! loop intervals in the parallel direction
        integer function check_periodicity(min_par,max_par) result(ierr)
            character(*), parameter :: rout_name = 'check_periodicity'
            
            ! input / output
            real(dp), intent(in) :: min_par, max_par                            ! min. and max. of parallel angle
            
            ! local variables
            real(dp) :: tol = 1E-5                                              ! tolerance
            real(dp) :: modulus
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! check whether max_par > min_par
            if (max_par.lt.min_par+2*pi*(1-tol)) then
                ierr = 1
                err_msg = 'max_par has to be at least min_par + 2 pi'
                CHCKERR(err_msg)
            end if
            
            ! calculate modulus
            modulus = mod(max_par-min_par,2*pi)
            modulus = min(modulus,2*pi-modulus)
            
            ! check whether tolerance is reached
            if (modulus/(2*pi).gt.tol) then
                ierr = 1
                err_msg = 'max_par - min_par has to be a multiple of 2 pi'
                CHCKERR(err_msg)
            end if
        end function check_periodicity
    end function calc_mesh
    
    ! check  whether   the  straight  field   line  mesh  has   been  calculated
    ! satisfactorily,  by  calculating  again  the zeta's  from  the  calculated
    ! theta's
    integer function check_mesh(alpha) result(ierr)
        use num_vars, only: tol_NR, theta_var_along_B, max_str_ln, &
            &calc_mesh_style
        
        character(*), parameter :: rout_name = 'check_mesh'
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        real(dp) :: var_calc(n_par,n_r)
        real(dp) :: var_diff(n_par,n_r)
        real(dp) :: var_in(n_r)
        integer :: id
        character(len=max_str_ln) :: par_ang, dep_ang                           ! parallel angle and dependent angle, for output message
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: work(:)                                        ! work array
        
        ! initialize ierr
        ierr = 0
        
        allocate(work(n_r))
        
        if (calc_mesh_style.eq.0) then                                          ! only do the check if the magnetic field lines are followed
            if (theta_var_along_B) then                                             ! calculate theta again from zeta
                do id = 1,n_par
                    var_in = zeta(id,:)
                    ierr = calc_ang_B(work,.true.,alpha,var_in)
                    var_calc(id,:) = work
                    CHCKERR('')
                end do
                par_ang = 'poloidal'; dep_ang = 'toroidal'
                var_diff = theta(:,:) - var_calc
            else                                                                    ! calculate zeta again from theta
                do id = 1,n_par
                    var_in = theta(id,:)
                    ierr = calc_ang_B(work,.false.,alpha,var_in)
                    var_calc(id,:) =  work
                    CHCKERR('')
                end do
                par_ang = 'toroidal'; dep_ang = 'poloidal'
                var_diff = zeta(:,:) - var_calc
            end if
            
            if (maxval(abs(var_diff)).gt.tol_NR*100) then
                err_msg = 'Calculating again the '//trim(par_ang)//&
                    &' grid points from the '//trim(dep_ang)//' grid points, &
                    &along the magnetic fields, does not result in the same &
                    &values as initally given... The maximum error is equal &
                    &to '//trim(r2strt(100*maxval(abs(var_diff))))//'%'
                ierr = 1
                CHCKERR(err_msg)
            end if
        end if
    end function check_mesh

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
    
    ! Calculates the  poloidal/toroidal angle theta(zeta)  as a function  of the
    ! toroidal/poloidal angle  zeta/theta following a particular  magnetic field
    ! line alpha.
    ! This is done using a Newton-Rhapson scheme that calculates the zero's of the
    ! function f(theta) = zeta - q (theta + lambda(theta)) - alpha
    ! the logical find_theta determines whether theta is sought or zeta:
    !   find_theta = .true. : look for theta as a function of zeta
    !   find_theta = .false. : look for zeta as a function of theta
    integer function calc_ang_B(ang_B,find_theta,alpha_in,input_ang,guess) result(ierr)
        use num_vars, only: tol_NR
        use VMEC_vars, only: mpol, ntor, L_c, L_s, iotaf, nfp
        use fourier_ops, only: calc_mesh_cs, f2r
        use utilities, only: calc_zero_NR, VMEC_conv_FHM
        
        character(*), parameter :: rout_name = 'ang_B'
        
        ! input / output
        real(dp) :: ang_B(n_r)                                                  ! theta(zeta)/zeta(theta)
        logical, intent(in) :: find_theta                                       ! whether theta or zeta is sought
        real(dp), intent(in) :: alpha_in                                        ! alpha
        real(dp), intent(in) :: input_ang(n_r)                                  ! the input angle zeta/theta
        real(dp), intent(in), optional :: guess(n_r)                            ! optional input (guess) for theta/zeta
        
        ! local variables (also used in child functions)
        integer :: kd                                                           ! counters
        real(dp) :: q(n_r)                                                      ! local version of 1/iota
        real(dp) :: lam                                                         ! lambda
        real(dp) :: dlam                                                        ! angular derivative of lambda
        
        ! local variables (not to be used in child functions)
        real(dp) :: alpha_calc                                                  ! calculated alpha, to check with the given alpha
        real(dp) :: ang_NR                                                      ! temporary solution for a given r, iteration
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! setup q
        q = 1/iotaf
        
        ! for all normal points
        ang_B = 0.0_dp
        norm: do kd = 1, n_r
            ! if first guess for theta is given
            if (present(guess)) then
                ang_NR = guess(kd)
            else if (kd.eq.1) then
                ang_NR = pi                                                     ! take pi, because it is in the middle of 0...2pi
            else                                                                ! take solution for previous flux surface
                ang_NR = ang_B(kd-1)
            end if
            
            ! Newton-Rhapson loop for current normal point
            ierr = calc_zero_NR(ang_B(kd),fun_ang_B,dfun_ang_B,ang_NR)
            CHCKERR('')
            
            ! do a check  whether the result is indeed alpha,  making use of the
            ! last lam and dlam that have been calculated in the child functions
            if (find_theta) then                                                ! looking for theta
                alpha_calc = input_ang(kd) - (ang_B(kd) + lam)*q(kd)
            else                                                                ! looking for zeta
                alpha_calc = ang_B(kd) - (input_ang(kd) + lam)*q(kd)
            end if
            if (alpha_calc-alpha_in.gt.tol_NR*100) then
                err_msg = 'In theta_B, calculating alpha as a check,&
                    & using the theta that is the solution of alpha '&
                    &//trim(r2strt(alpha_in))//', yields alpha_calc that &
                    &deviates '//trim(r2strt(100*abs((alpha_calc-alpha_in))))&
                    &//'% from the original alpha_in'
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do norm
    contains
        ! function that returns  f = alpha -  alpha_0. It uses kd  from the main
        ! loop in the  parent function as the normal position  where to evaluate
        ! the quantitites, and input_ang(kd) as the angle for which the magnetic
        ! match is sought, as well as q, alpha_in, L_s and L_c
        function fun_ang_B(ang)
            character(*), parameter :: rout_name = 'fun_ang_B'
            
            ! input / output
            real(dp) :: fun_ang_B
            real(dp), intent(in) :: ang
            
            ! local variables
            real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                               ! factors of cosine and sine for inverse fourier transform
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            if (find_theta) then                                                ! looking for theta
                ierr = calc_mesh_cs(cs,mpol,ntor,nfp,ang,input_ang(kd))
                CHCKERR('')
            else                                                                ! looking for zeta
                ierr = calc_mesh_cs(cs,mpol,ntor,nfp,input_ang(kd),ang)
                CHCKERR('')
            end if
            
            ! calculate lambda
            lam = f2r(L_c(:,:,kd,0),L_s(:,:,kd,0),cs,mpol,ntor,nfp,ierr=ierr)
            CHCKERR('')
            
            ! calculate the output function
            if (find_theta) then                                                ! looking for theta
                fun_ang_B = input_ang(kd) - (ang+lam)*q(kd) - alpha_in
            else                                                                ! looking for zeta
                fun_ang_B = ang - (input_ang(kd)+lam)*q(kd) - alpha_in
            end if
        end function fun_ang_B
        
        ! function that returns  df/d(ang) = d(alpha -  alpha_0)/d(ang). It uses
        ! kd from  the main loop in  the parent function as  the normal position
        ! where to evaluate the quantitites,  and input_ang(kd) as the angle for
        ! which the  magnetic match is sought,  as well as q,  alpha_in, L_s and
        ! L_c
        function dfun_ang_B(ang)
            character(*), parameter :: rout_name = 'dfun_ang_B'
            
            ! input / output
            real(dp) :: dfun_ang_B
            real(dp), intent(in) :: ang
            
            ! local variables
            real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                               ! factors of cosine and sine for inverse fourier transform
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            if (find_theta) then                                                ! looking for theta
                ierr = calc_mesh_cs(cs,mpol,ntor,nfp,ang,input_ang(kd))
                CHCKERR('')
            else                                                                ! looking for zeta
                ierr = calc_mesh_cs(cs,mpol,ntor,nfp,input_ang(kd),ang)
                CHCKERR('')
            end if
            
            ! calculate angular derivatives of lambda
            if (find_theta) then                                                ! looking for theta
                dlam = f2r(L_c(:,:,kd,0),L_s(:,:,kd,0),cs,mpol,ntor,nfp,&
                    &[1,0],ierr)
                CHCKERR('')
            else                                                                ! looking for zeta
                dlam = f2r(L_c(:,:,kd,0),L_s(:,:,kd,0),cs,mpol,ntor,nfp,&
                    &[0,1],ierr)
                CHCKERR('')
            end if
            
            ! calculate the output function
            if (find_theta) then                                                ! looking for theta
                dfun_ang_B = -(1.0_dp + dlam)*q(kd)
            else                                                                ! looking for zeta
                dfun_ang_B = 1.0_dp - dlam*q(kd)
            end if
        end function dfun_ang_B
    end function calc_ang_B
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! equilibrium phase
    subroutine dealloc_eq_vars
        deallocate(VMEC_R,VMEC_Z,VMEC_L)
        deallocate(flux_p,flux_p_FD,flux_t,flux_t_FD)
    end subroutine dealloc_eq_vars
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! calculation for a certain alpha
    subroutine dealloc_eq_vars_final
        deallocate(q_saf,q_saf_FD)
        deallocate(pres,pres_FD)
        deallocate(zeta,theta)
    end subroutine dealloc_eq_vars_final
end module eq_vars
