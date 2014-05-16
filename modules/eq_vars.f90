!------------------------------------------------------------------------------!
!   Variables,  subroutines  and functions  that  have  to do  with            !
!   equilibrium quantities and the mesh used in the calculations               !
!   These are called in the subroutine calc_eq in the module eq_ops            !
!------------------------------------------------------------------------------!
module eq_vars
    use num_vars, only: dp, pi
    use output_ops, only: writo, lvl_ud, print_ar_2, print_ar_1, write_out
    use str_ops, only: r2strt, i2str, r2str

    implicit none
    private
    public eqd_mesh, calc_mesh, ang_B, calc_flux_q, check_mesh, &
        &init_eq, calc_RZL, &
        &theta, zeta, theta_H, zeta_H, n_par, R, Z, lam_H, min_par, max_par, &
        &q_saf, flux_p, flux_t, VMEC_R, VMEC_Z, VMEC_L, pres

    ! R and Z and derivatives in real (as opposed to Fourier) space (see below)
    ! (index 1: variable, 2: r derivative, 2: theta derivative, 3: zeta derivative)
    real(dp), allocatable :: R(:,:,:), Z(:,:,:)                                 ! R, Z (FM)
    real(dp), allocatable :: VMEC_R(:,:,:,:,:)                                  ! R in VMEC coordinates (n_par, n_r, [derivatives])
    real(dp), allocatable :: VMEC_Z(:,:,:,:,:)                                  ! Z in VMEC coordinates (n_par, n_r, [derivatives])
    real(dp), allocatable :: VMEC_L(:,:,:,:,:)                                  ! L(ambda) in VMEC coordinates (n_par, n_r, [derivatives])
    real(dp), allocatable :: lam_H(:,:,:)                                       ! lambda in (HM)
    real(dp), allocatable :: theta(:,:), zeta(:,:)                              ! grid points (n_par, n_r, 2) (FM)
    real(dp), allocatable :: theta_H(:,:), zeta_H(:,:)                          ! grid points (n_par, n_r, 2) (HM)
    real(dp), allocatable :: q_saf(:,:)                                         ! safety factor (FM)
    real(dp), allocatable :: flux_p(:,:), flux_t(:,:), pres(:,:)                ! pol. flux, tor. flux and pressure, and normal derivative (FM)
    real(dp) :: min_par, max_par
    integer :: n_par
    
    interface calc_RZL
        module procedure calc_RZL_ind, calc_RZL_arr
    end interface

contains
    ! initialize the variables VMEC_R, VMEC_Z and VMEC_L
    subroutine init_eq
        use num_vars, only: max_deriv
        use VMEC_vars, only: n_r
        
        ! R
        if (allocated(VMEC_R)) deallocate(VMEC_R)
        allocate(VMEC_R(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! Z
        if (allocated(VMEC_Z)) deallocate(VMEC_Z)
        allocate(VMEC_Z(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! lambda
        if (allocated(VMEC_L)) deallocate(VMEC_L)
        allocate(VMEC_L(n_par,n_r,0:max_deriv(1),0:max_deriv(2),0:max_deriv(3)))
        
        ! q_saf
        if (allocated(q_saf)) deallocate(q_saf)
        allocate(q_saf(n_r,0:max_deriv(1)))
        
        ! flux_p
        if (allocated(flux_p)) deallocate(flux_p)
        allocate(flux_p(n_r,0:max_deriv(1)))
        
        ! flux_t
        if (allocated(flux_t)) deallocate(flux_t)
        allocate(flux_t(n_r,0:max_deriv(1)))
        
        ! pres
        if (allocated(pres)) deallocate(pres)
        allocate(pres(n_r,0:max_deriv(1)))
    end subroutine
    
    ! calculate R, Z and Lambda and  derivatives in VMEC coordinates at the mesh
    ! points given by the variables VMEC_theta and VMEC_zeta and at every normal
    ! point. The derivatives  are indicated by the variable "deriv"  which has 3
    ! indices
    subroutine calc_RZL_ind(deriv)
        use fourier_ops, only: mesh_cs, f2r
        use VMEC_vars, only: mpol, ntor, n_r, R_c, R_s, Z_c, Z_s, L_c, L_s
        use utilities, only: check_deriv
        use num_vars, only: max_deriv
        
        ! input / output
        integer, intent(in) :: deriv(3)
        
        ! local variables
        integer :: id, kd
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        
        ! check the derivatives requested
        call check_deriv(deriv,max_deriv,'calc_RZL')
        
        ! calculate the  variables R, Z,  and their angular derivatives  for all
        ! normal points and current angular point
        perp: do kd = 1, n_r                                                    ! perpendicular: normal to the flux surfaces
            par: do id = 1, n_par                                               ! parallel: along the magnetic field line
                cs = mesh_cs(mpol,ntor,theta(id,kd),zeta(id,kd))
                VMEC_R(id,kd,deriv(1),deriv(2),deriv(3)) = &
                    &f2r(R_c(:,:,kd,deriv(1)),R_s(:,:,kd,deriv(1)),cs&
                    &,mpol,ntor,[deriv(2),deriv(3)])
                VMEC_Z(id,kd,deriv(1),deriv(2),deriv(3)) = &
                    &f2r(Z_c(:,:,kd,deriv(1)),Z_s(:,:,kd,deriv(1)),cs&
                    &,mpol,ntor,[deriv(2),deriv(3)])
                VMEC_L(id,kd,deriv(1),deriv(2),deriv(3)) = &
                    &f2r(L_c(:,:,kd,deriv(1)),L_s(:,:,kd,deriv(1)),cs&
                    &,mpol,ntor,[deriv(2),deriv(3)])
            end do par
        end do perp
    end subroutine
    subroutine calc_RZL_arr(deriv)
        ! input / output
        integer, intent(in) :: deriv(:,:)
        
        ! local variables
        integer :: id
        
        do id = 1, size(deriv,2)
            call calc_RZL_ind(deriv(:,id))
        end do
    end subroutine

    ! Calculates flux quantities and normal derivatives
    subroutine calc_flux_q
        use VMEC_vars, only: &
            &iotaf, n_r, phi, phi_r, presf
        use utilities, only: VMEC_norm_deriv, calc_int
        use num_vars, only: max_deriv
        
        ! local variables
        integer :: kd
        
        ! safety factor q_saf: invert iota and derivate
        do kd = 1,n_r
            q_saf(kd,0) = 1/iotaf(kd)
        end do
        do kd = 1,max_deriv(1)
            call VMEC_norm_deriv(q_saf(:,0),q_saf(:,kd),n_r-1._dp,kd,1)
        end do
        
        ! toroidal flux: copy from VMEC and interpolate
        flux_t(:,0) = phi
        flux_t(:,1) = phi_r
        do kd = 2,max_deriv(1)
            call VMEC_norm_deriv(flux_t(:,1),flux_t(:,kd),n_r-1._dp,kd-1,1)
        end do
        
        ! poloidal flux: calculate using iotaf and phi, phi_r
        flux_p(:,1) = iotaf*phi_r
        flux_p(:,0) = calc_int(flux_p(:,1),[(kd/(n_r-1.0_dp),kd=0,n_r-1)])
        do kd = 2,max_deriv(1)
            call VMEC_norm_deriv(flux_p(:,1),flux_p(:,kd),n_r-1._dp,kd-1,1)
        end do
        
        ! pressure: copy from VMEC and derivate
        pres(:,0) = presf
        do kd = 1, max_deriv(1)
            call VMEC_norm_deriv(pres(:,0),pres(:,kd),n_r-1._dp,kd,1)
        end do
        !do kd = 0,max_deriv(1)
            !call write_out(1,n_r,pres(:,kd),'pres^'//trim(i2str(kd)))
        !end do
    end subroutine



    





    ! calculate the angular mesh, filling the global variables (theta, zeta). 
    ! This is done for every flux surface. 
    ! The variable  theta_var_along_B determines  whether theta  is used  as the
    ! base variable. If .false., zeta is used.
    ! ¡¡¡SHOULD BE DONE ADAPTIVELY!!!
    ! ¡¡¡NEED A CHECK FOR THE FINENESS!!!
    subroutine calc_mesh(alpha)
        use VMEC_vars, only: n_r
        use num_vars, only: theta_var_along_B
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        integer :: jd, kd
        real(dp) :: var_in(n_r)
        
        if (allocated(zeta)) deallocate(zeta)
        allocate(zeta(n_par,n_r)); zeta = 0.0_dp
        if (allocated(zeta_H)) deallocate(zeta_H)
        allocate(zeta_H(n_par,n_r)); zeta_H = 0.0_dp
        if (allocated(theta)) deallocate(theta)
        allocate(theta(n_par,n_r)); theta = 0.0_dp
        if (allocated(theta_H)) deallocate(theta_H)
        allocate(theta_H(n_par,n_r)); theta_H = 0.0_dp
        
        if (theta_var_along_B) then                                             ! first calculate theta
            ! ¡¡¡¡¡¡¡¡¡¡ TO CALCULATE RADIAL DERIVATIVES YOU NEED COORDINATES THETA
            !            AND ZETA (VMEC) TO BE CONSTANT IN DIFFERENT FLUX SURFACES
            !            -> YOU CANNOT FOLLOW A MAGNETIC FIELD JUST LIKE THIS!!!
            !            -> THIS WOULD YIELD THE NORMAL DERIVATIVE WITH ALPHA = CONSTANT
            !               WHICH IS IN THE FLUX COORDINATE SYSTEM
            ! !¡¡¡¡¡¡¡ ANOTHER PROBLEM: YOU SHOULD FIND THIS AS A FUNCTION OF THETA_F, 
            !          WITHIN A CERTAIN RANGE, NOT OF THETA_V!
            !          -> FIND EQUIDISTANT MESH IN THETA_F AND LOOK FOR THETA_V FIRST, 
            !             FROM WHICH YOU CAN FIND ZETA_V
            ! ¡¡¡¡¡¡! TEMPORARILY REPLACED !!!!!
            !do jd = 1,n_r
                !theta(:,jd) = eqd_mesh(n_par, min_par, max_par)
                !theta_H(:,jd) = eqd_mesh(n_par, min_par, max_par)
            !end do
            !do kd = 1,n_par
                !var_in = theta(kd,:)
                !zeta(kd,:) = ang_B(.false.,.true.,alpha,var_in)
                !var_in = theta_H(kd,:)
                !zeta-H(kd,:) = ang_B(.false.,.false.,alpha,var_in)
            !end do
            ! TEMPORARY REPLACEMENT:
            theta = 0.4*pi/2
            theta_H = 0.4*pi/2
            do jd = 1,n_r
                zeta(:,jd) = eqd_mesh(n_par, min_par, max_par)
                zeta_H(:,jd) = eqd_mesh(n_par, min_par, max_par)
            end do
        else                                                                    ! first calculate zeta
            do jd = 1,n_r
                zeta(:,jd) = eqd_mesh(n_par, min_par, max_par)
                zeta_H(:,jd) = eqd_mesh(n_par, min_par, max_par)
            end do
            do kd = 1,n_par
                var_in = zeta(kd,:)
                theta(kd,:) = ang_B(.true.,.true.,alpha,var_in)
                var_in = zeta(kd,:)
                theta_H(kd,:) = ang_B(.true.,.false.,alpha,var_in)
            end do
        end if
    end subroutine calc_mesh
    
    ! check  whether   the  straight  field   line  mesh  has   been  calculated
    ! satisfactorily,  by  calculating  again  the zeta's  from  the  calculated
    ! theta's
    subroutine check_mesh(alpha)
        use VMEC_vars, only: n_r
        use num_vars, only: tol_NR, theta_var_along_B, max_str_ln
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        real(dp) :: var_calc(n_par,n_r)
        real(dp) :: var_diff(n_par,n_r)
        real(dp) :: var_in(n_r)
        integer :: id
        character(len=max_str_ln) :: par_ang, dep_ang                           ! parallel angle and dependent angle, for output message
        
        if (theta_var_along_B) then                                             ! calculate theta again from zeta
            do id = 1,n_par
                var_in = zeta(id,:)
                var_calc(id,:) = ang_B(.true.,.true.,alpha,var_in)
            end do
            par_ang = 'poloidal'; dep_ang = 'toroidal'
            var_diff = theta(:,:) - var_calc
        else                                                                    ! calculate zeta again from theta
            do id = 1,n_par
                var_in = theta(id,:)
                var_calc(id,:) = ang_B(.false.,.true.,alpha,var_in)
            end do
            par_ang = 'toroidal'; dep_ang = 'poloidal'
            var_diff = zeta(:,:) - var_calc
        end if
        
        if (maxval(abs(var_diff)).gt.tol_NR*100) then
            call writo('ERROR: Calculating again the '//trim(par_ang)//&
                &' grid points from the '//trim(dep_ang)//' grid points, &
                &along the magnetic fields, does not result in the same &
                &values as initally given.')
            call writo('The maximum error is equal to '//&
                &trim(r2strt(100*maxval(abs(var_diff))))//'%')
            stop
        end if
    end subroutine check_mesh

    ! calculate mesh of equidistant points
    function eqd_mesh(n_ang, min_ang, max_ang)
        ! input and output
        real(dp), allocatable :: eqd_mesh(:)
        real(dp), intent(in) :: min_ang, max_ang
        integer, intent(in) :: n_ang
        
        ! local variables
        integer :: id
        real(dp) :: delta_ang
        
        ! test some values
        if (n_ang.lt.1) then
            call writo('ERROR: the angular array has to have a length of &
                &at least 1')
            stop
        end if
        if (min_ang.gt.max_ang) then
            call writo('ERROR: the minimum angle has to be smaller than &
                &the maximum angle')
            stop
        end if
        
        ! initialize output vector
        allocate(eqd_mesh(n_ang)); eqd_mesh = 0.0_dp
        ! There are (n_ang-1) pieces in the total interval but the last one 
        ! is not needed as the functions are all periodic
        
        delta_ang = (max_ang-min_ang)/(n_ang)
        
        eqd_mesh(1) = min_ang
        do id = 2,n_ang
            eqd_mesh(id) = eqd_mesh(id-1) + delta_ang
        end do
    end function eqd_mesh
    
    ! Calculates the  poloidal/toroidal angle theta(zeta)  as a function  of the
    ! toroidal/poloidal angle  zeta/theta following a particular  magnetic field
    ! line alpha.
    ! This is done using a Newton-Rhapson scheme that calculates the zero's of the
    ! function f(theta) = zeta - q (theta + lambda(theta)) - alpha
    ! the logical find_theta determines whether theta is sought or zeta:
    !   find_theta = .true. : look for theta as a function of zeta
    !   find_theta = .false. : look for zeta as a function of theta
    function ang_B(find_theta,fm ,alpha_in,input_ang,guess)
        use num_vars, only: tol_NR
        use VMEC_vars, only: mpol, ntor, n_r, L_c, L_s, iotaf, iotah
        use fourier_ops, only: mesh_cs, f2r
        use utilities, only: zero_NR
        
        ! input / output
        real(dp) :: ang_B(n_r)                                                  ! theta(zeta)/zeta(theta)
        logical, intent(in) :: find_theta                                       ! whether theta or zeta is sought
        logical, intent(in) :: FM                                               ! whether full or half mesh quantity is sought
        real(dp), intent(in) :: alpha_in                                        ! alpha
        real(dp), intent(in) :: input_ang(n_r)                                  ! the input angle zeta/theta
        real(dp), intent(in), optional :: guess(n_r)                            ! optional input (guess) for theta/zeta
        
        !! local variables (also used in child functions)
        !integer :: jd, kd
        !real(dp) :: L_c_loc(0:mpol-1,-ntor:ntor,1:n_r)                        ! local version of HM l_c (either FM or HM)
        !real(dp) :: L_s_loc(0:mpol-1,-ntor:ntor,1:n_r)                        ! local version of HM l_s (either FM or HM)
        !real(dp) :: q(n_r)                                                      ! local version of 1/iota (either FM or HM)
        !real(dp) :: lam                                                         ! lambda
        !real(dp) :: dlam                                                        ! angular derivative of lambda
        !real(dp) :: tempcoeff(n_r)                                              ! temporary holds a coefficient
        
        !! local variables (not to be used in child functions)
        !real(dp) :: alpha_calc                                                  ! calculated alpha, to check with the given alpha
        !real(dp) :: ang_NR                                                      ! temporary solution for a given r, iteration
        !integer :: norm_start                                                   ! either 1 (FM) or 2 (HM)
        
        !! convert L_c and L_s to FM if FM is .true., select correct q
        !if (FM) then
            !do jd = -ntor, ntor
                !do kd = 0, mpol-1
                    !tempcoeff = L_c(kd,jd,:,0)
                    !L_c_loc(kd,jd,:) = h2f(tempcoeff)
                    !tempcoeff = L_s(kd,jd,:,0)
                    !L_s_loc(kd,jd,:) = h2f(tempcoeff)
                !end do
            !end do
            !q = iotaf
        !else 
            !L_c_loc = L_c(:,:,:,0)
            !L_s_loc = L_s(:,:,:,0)
            !q = iotah
        !end if
        
        !! for all normal points
        !if (fm) then
            !norm_start = 1
        !else
            !norm_start = 2
        !end if
        !ang_B = 0.0_dp
        !norm: do kd = 1, n_r
            !! if first guess for theta is given
            !if (present(guess)) then
                !ang_NR = guess(kd)
            !else if (kd.eq.1) then
                !ang_NR = pi                                                     ! take pi, because it is in the middle of 0...2pi
            !else                                                                ! take solution for previous flux surface
                !ang_NR = ang_B(kd-1)
            !end if
            
            !! Newton-Rhapson loop for current normal point
            !ang_B(kd) = zero_NR(fun_ang_B,dfun_ang_B,ang_NR)
            
            !! do a check  whether the result is indeed alpha,  making use of the
            !! last lam and dlam that have been calculated in the child functions
            !if (find_theta) then                                                ! looking for theta
                !alpha_calc = input_ang(kd) - (ang_B(kd) + lam)*q(kd)
            !else                                                                ! looking for zeta
                !alpha_calc = ang_B(kd) - (input_ang(kd) + lam)*q(kd)
            !end if
            !if (alpha_calc-alpha_in.gt.tol_NR*100) then
                !call writo('ERROR: In theta_B, calculating alpha as a check,&
                    !& using the theta that is the solution of alpha '&
                    !&//trim(r2strt(alpha_in))//', yields alpha_calc that &
                    !&deviates '//trim(r2strt(100*abs((alpha_calc-alpha_in))))&
                    !&//'% from the original alpha_in')
                !stop
            !end if
        !end do norm
        
    !contains
        !! function that returns  f = alpha -  alpha_0. It uses kd  from the main
        !! loop in the  parent function as the normal position  where to evaluate
        !! the quantitites, and input_ang(kd) as the angle for which the magnetic
        !! match is sought, as well as q, alpha_in, L_s_loc and L_c_loc
        !function fun_ang_B(ang)
            !! input / output
            !real(dp) :: fun_ang_B
            !real(dp), intent(in) :: ang
            
            !! local variables
            !real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                               ! factors of cosine and sine for inverse fourier transform
            
            !! transform lambda from Fourier space to real space
            !! calculate the (co)sines
            !if (find_theta) then                                                ! looking for theta
                !cs = mesh_cs(mpol,ntor,ang,input_ang(kd))
            !else                                                                ! looking for zeta
                !cs = mesh_cs(mpol,ntor,input_ang(kd),ang)
            !end if
            
            !! calculate lambda
            !lam = f2r(L_c_loc(:,:,kd),L_s_loc(:,:,kd),cs,mpol,ntor)
            
            !! calculate the output function
            !if (find_theta) then                                                ! looking for theta
                !fun_ang_B = input_ang(kd) - (ang+lam)*q(kd) - alpha_in
            !else                                                                ! looking for zeta
                !fun_ang_B = ang - (input_ang(kd)+lam)*q(kd) - alpha_in
            !end if
        !end function fun_ang_B
        
        !! function that returns  df/d(ang) = d(alpha -  alpha_0)/d(ang). It uses
        !! kd from  the main loop in  the parent function as  the normal position
        !! where to evaluate the quantitites,  and input_ang(kd) as the angle for
        !! which the  magnetic match is sought,  as well as q,  alpha_in, L_s_loc
        !! and L_c_loc
        !function dfun_ang_B(ang)
            !! input / output
            !real(dp) :: dfun_ang_B
            !real(dp), intent(in) :: ang
            
            !! local variables
            !real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                               ! factors of cosine and sine for inverse fourier transform
            
            !! transform lambda from Fourier space to real space
            !! calculate the (co)sines
            !if (find_theta) then                                                ! looking for theta
                !cs = mesh_cs(mpol,ntor,ang,input_ang(kd))
            !else                                                                ! looking for zeta
                !cs = mesh_cs(mpol,ntor,input_ang(kd),ang)
            !end if
            
            !! calculate angular derivatives of lambda
            !if (find_theta) then                                                ! looking for theta
                !dlam = f2r(L_c_loc(:,:,kd),L_s_loc(:,:,kd),cs,mpol,ntor,[1,0])
            !else                                                                ! looking for zeta
                !dlam = f2r(L_c_loc(:,:,kd),L_s_loc(:,:,kd),cs,mpol,ntor,[0,1])
            !end if
            
            !! calculate the output function
            !if (find_theta) then                                                ! looking for theta
                !dfun_ang_B = -(1.0_dp + dlam)*q(kd)
            !else                                                                ! looking for zeta
                !dfun_ang_B = 1.0_dp - dlam*q(kd)
            !end if
        !end function dfun_ang_B
    end function ang_B
end module eq_vars
