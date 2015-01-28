!------------------------------------------------------------------------------!
!   Operations that have to do with the different coordinate systems           !
!------------------------------------------------------------------------------!
module coord_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, pi, max_str_ln
    use str_ops, only: r2strt, i2str, r2str
    use output_ops, only: print_GP_2D, print_GP_3D, draw_GP_animated
    use messages, only: lvl_ud, writo

    implicit none
    private
    public coord_F2E, coord_E2F, calc_XYZ_grid, calc_eqd_grid, calc_ang_grid, &
        &plot_grid_real
    
contains
    ! Calculate  the angular grid  in the Flux (perturbation)  angular variable,
    ! but tabulated in the Equilibrium (!) normal variable.
    ! The variable  use_pol_flux determines  whether theta (.true.)  or zeta
    ! (.false.) is used as the parallel variable.
    integer function calc_ang_grid(alpha) result(ierr)
        use num_vars, only: grid_style, use_pol_flux_X, use_pol_flux_eq, &
            &eq_style
        use X_vars, only: n_par, min_par, max_par, ang_par_F
        use eq_vars, only: theta_E, zeta_E, flux_p_E, flux_t_E, &
            &grp_n_r_eq, q_saf_E, rot_t_E
        
        character(*), parameter :: rout_name = 'calc_ang_grid'
        
        ! input / output
        real(dp), intent(in) :: alpha                                           ! field line coordinate
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: ang_comp_F(:,:)                                ! complementary Flux angle to ang_par_F
        integer :: kd                                                           ! counter
        real(dp), allocatable :: r_E_loc(:)                                     ! flux in Equilibrium coords.
        real(dp), allocatable :: theta_F(:,:,:), zeta_F(:,:,:)                  ! theta_F and zeta_F
        real(dp), allocatable :: theta_E_loc(:,:,:), zeta_E_loc(:,:,:)          ! local copy of theta_E and zeta_E
        real(dp), pointer :: flux(:), flux_eq(:)                                ! either pol. or tor. flux
        integer :: pmone                                                        ! plus or minus one
        
        ! initialize ierr
        ierr = 0
        
        ! calculate the grid
        ! grid style
        !   1:  grid along magnetic field lines
        select case (grid_style)
            ! grid along the magnetic field lines
            case (1)
                ! set up parallel and complementary Flux coordinate
                allocate(ang_par_F(n_par,grp_n_r_eq)); ang_par_F = 0.0_dp
                allocate(ang_comp_F(n_par,grp_n_r_eq)); ang_comp_F = 0.0_dp
                
                ! set up  parallel angle in Flux  coordinates on equidistant
                ! grid
                ierr = calc_eqd_grid(ang_par_F(:,1),n_par,min_par,max_par)
                CHCKERR('')
                do kd = 2,grp_n_r_eq
                    ang_par_F(:,kd) = ang_par_F(:,1)
                end do
                
                ! set up plus minus one
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        pmone = -1                                              ! conversion VMEC LH -> RH coord. system
                    case (2)                                                    ! HELENA
                        pmone = 1
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! calculate the other Flux angle from alpha
                if (use_pol_flux_X) then
                    do kd = 1,grp_n_r_eq
                        ang_comp_F(:,kd) = pmone*q_saf_E(kd,0)*ang_par_F(:,kd)
                    end do
                    ang_comp_F = ang_comp_F + alpha
                else
                    do kd = 1,grp_n_r_eq
                        ang_comp_F(:,kd) = pmone*rot_t_E(kd,0)*ang_par_F(:,kd)
                    end do
                    ang_comp_F = ang_comp_F - alpha
                end if
                
                ! allocate flux coordinates
                allocate(theta_F(1,n_par,grp_n_r_eq),zeta_F(1,n_par,grp_n_r_eq))
                
                ! set up local equilibrium coordinates
                allocate(r_E_loc(grp_n_r_eq))
                allocate(theta_E_loc(1,n_par,grp_n_r_eq))
                allocate(zeta_E_loc(1,n_par,grp_n_r_eq))
                
                ! set up theta_F and zeta_F
                if (use_pol_flux_X) then
                    theta_F(1,:,:) = ang_par_F
                    zeta_F(1,:,:) = ang_comp_F
                else
                    theta_F(1,:,:) = ang_comp_F
                    zeta_F(1,:,:) = ang_par_F
                end if
                
                ! set up flux and flux_eq
                if (use_pol_flux_X) then
                    flux => flux_p_E(:,0)
                else
                    flux => flux_t_E(:,0)
                end if
                if (use_pol_flux_eq) then
                    flux_eq => flux_p_E(:,0)
                else
                    flux_eq => flux_t_E(:,0)
                end if
                
                ! convert Flux coordinates to equilibrium coordinates
                ierr = coord_F2E(flux,theta_F,zeta_F,&
                    &r_E_loc,theta_E_loc,zeta_E_loc,&
                    &r_F_array=flux,r_E_array=flux_eq)
                CHCKERR('')
                
                ! update theta_E and zeta_E with their local copies
                theta_E = theta_E_loc(1,:,:)
                zeta_E = zeta_E_loc(1,:,:)
                
                ! deallocate local variables
                deallocate(ang_comp_F)
                deallocate(theta_E_loc,zeta_E_loc,r_E_loc)
                deallocate(theta_F,zeta_F)
                nullify(flux,flux_eq)
            ! grid style error
            case default
                err_msg = 'No grid style is associated with '&
                    &//trim(i2str(grid_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function calc_ang_grid
    
    ! Converts  Flux  coordinates  (r,theta,zeta)_F to  Equilibrium  coordinates
    ! (r,theta,zeta)_E. Optionally,  two arrays  r_F_array and r_E_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, the normalized (wrt. 1) flux variables are used.
    integer function coord_F2E(r_F,theta_F,zeta_F,r_E,theta_E,zeta_E,&
        &r_F_array,r_E_array) result(ierr)
        use num_vars, only: tol_NR, use_pol_flux_X, use_pol_flux_eq, eq_style
        use fourier_ops, only: fourier2real, calc_trigon_factors
        use VMEC, only: L_c, L_s
        use utilities, only: interp_fun
        use eq_vars, only: flux_p_FD, flux_t_FD, max_flux_F, max_flux_eq_F
        
        character(*), parameter :: rout_name = 'coord_F2E'
        
        ! input / output
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        real(dp), intent(inout) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)        ! equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        real(dp), allocatable :: L_c_loc(:,:,:)                                 ! local version of L_c
        real(dp), allocatable :: L_s_loc(:,:,:)                                 ! local version of L_s
        real(dp) :: theta_V_loc(1,1), zeta_V_loc(1,1)                           ! local version of theta_V, zeta_V (needs to be 2D matrix)
        real(dp) :: lam(1,1)                                                    ! lambda (needs to be 2D matrix)
        real(dp) :: dlam(1,1)                                                   ! angular derivative of lambda (needs to be 2D matrix)
        real(dp), pointer :: flux(:), flux_eq(:)                                ! either pol. or tor. flux
        integer :: id, jd, kd                                                   ! counters
        real(dp) :: r_F_factor, r_E_factor                                      ! mult. factors for r_F and r_E
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_theta = size(theta_E,1)
        n_zeta = size(theta_E,2)
        n_r = size(theta_E,3)
        
        ! tests
        if (size(theta_F,1).ne.size(zeta_F,1) .or. &
            &size(theta_F,2).ne.size(zeta_F,2) .or. &
            &size(theta_F,3).ne.size(zeta_F,3) .or. &
            &size(theta_F,3).ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(zeta_E,1) .or. n_zeta.ne.size(zeta_E,2) .or. &
            &n_r.ne.size(zeta_E,3) .or. n_r.ne.size(r_E)) then
            ierr = 1
            err_msg = 'theta_E, zeta_E and r_E need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_F_array).neqv.present(r_E_array)) then
            ierr = 1
            err_msg = 'both r_F_array and r_E_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! set up flux, flux_eq and multiplicative factors r_F_factor, r_E_factor
        if (present(r_F_array)) then
            flux => r_F_array
            flux_eq => r_E_array
            r_F_factor = 1._dp
            r_E_factor = 1._dp
        else
            if (use_pol_flux_X) then
                flux => flux_p_FD(:,0)
            else
                flux => flux_t_FD(:,0)
            end if
            if (use_pol_flux_eq) then
                flux_eq => flux_p_FD(:,0)
            else
                flux_eq => flux_t_FD(:,0)
            end if
            r_F_factor = max_flux_F
            r_E_factor = max_flux_eq_F
        end if
        
        ! convert normal position
        do kd = 1,n_r
            ierr = interp_fun(r_E(kd),flux_eq,r_F(kd)*r_F_factor,flux)
            CHCKERR('')
            r_E(kd) = r_E(kd)*r_E_factor
        end do
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = coord_F2E_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! trivial HELENA uses flux coordinates
                theta_E = theta_F
                zeta_E = zeta_F
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        integer function coord_F2E_VMEC() result(ierr)
            use utilities, only: calc_zero_NR
            use VMEC, only: mpol, ntor
            use eq_vars, only: grp_min_r_eq, grp_max_r_eq
            
            character(*), parameter :: rout_name = 'coord_F2E_VMEC'
            
            ! local variables
            real(dp), allocatable :: L_c_loc_ind(:,:), L_s_loc_ind(:,:)         ! individual versions of L_c_loc and L_s_loc
            
            ! initialize ierr
            ierr = 0
            
            ! the toroidal angle is trivial
            zeta_E = - zeta_F                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_c and L_s
            allocate(L_c_loc(0:mpol-1,-ntor:ntor,1:1))
            allocate(L_s_loc(0:mpol-1,-ntor:ntor,1:1))
            allocate(L_c_loc_ind(0:mpol-1,-ntor:ntor))
            allocate(L_s_loc_ind(0:mpol-1,-ntor:ntor))
            
            ! loop over all normal points
            do kd = 1,n_r
                ! interpolate L_c and L_s at requested normal point r_F
                ierr = interp_fun(L_c_loc_ind,&
                    &L_c(:,:,grp_min_r_eq:grp_max_r_eq,0),&
                    &r_F(kd)*r_F_factor,flux)
                CHCKERR('')
                ierr = interp_fun(L_s_loc_ind,&
                    &L_s(:,:,grp_min_r_eq:grp_max_r_eq,0),&
                    &r_F(kd)*r_F_factor,flux)
                CHCKERR('')
                
                ! copy individual to array version
                L_c_loc(:,:,1) = L_c_loc_ind
                L_s_loc(:,:,1) = L_s_loc_ind
                
                ! the poloidal angle has to be found as the zero of
                !   f = theta_F - theta_E - lambda
                ! loop over all angular points
                do jd = 1,n_zeta
                    do id = 1,n_theta
                        ! calculate zero of f
                        ierr = calc_zero_NR(theta_E(id,jd,kd),fun_pol,dfun_pol,&
                            &theta_F(id,jd,kd))                                 ! use theta_F as guess for theta_E
                        CHCKERR('')
                        
                        ! do a check whether the result is indeed zero
                        if (abs(fun_pol(theta_E(id,jd,kd))).gt.tol_NR*100) then
                            err_msg = 'In coord_F2E_VMEC, calculating f as a &
                                &check, using the theta_E that is the solution &
                                &of f = 0, yields a calculated f that &
                                &deviates from 0 by '//trim(r2strt(&
                                &100*abs(fun_pol(theta_E(id,jd,kd)))))//'%'
                            ierr = 1
                            CHCKERR(err_msg)
                        end if
                    end do
                end do
            end do
            
            ! deallocate variables
            deallocate(L_c_loc,L_s_loc)
            deallocate(L_c_loc_ind,L_s_loc_ind)
        end function coord_F2E_VMEC
        
        ! function that returns f = theta_F  - theta_V - lambda. It uses theta_F
        ! and  zeta_E (=  zeta_V)  from the  parent function,  lam  and dlam  to
        ! contain the  variable lambda  or its derivative,  as well  as L_s_loc,
        ! L_c_loc, theta_V_loc, zeta_V_loc, id, jd and kd.
        function fun_pol(theta_E_in)
            character(*), parameter :: rout_name = 'fun_pol'
            
            ! input / output
            real(dp) :: fun_pol
            real(dp), intent(in) :: theta_E_in
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:)              ! trigonometric factor cosine for the inverse fourier transf.
            
            ! initialize fun_pol
            fun_pol = 0.0_dp
            
            ! set up local copies of theta_V_loc, zeta_V_loc
            theta_V_loc = theta_E_in
            zeta_V_loc = zeta_E(id,jd,kd)
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            ierr = calc_trigon_factors(theta_V_loc,zeta_V_loc,&
                &trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_c_loc,L_s_loc,trigon_factors_loc,lam)
            CHCKERR('')
            
            ! calculate the output function
            fun_pol = theta_F(id,jd,kd) - theta_E_in - lam(1,1)
            
            ! deallocate trigoniometric factors
            deallocate(trigon_factors_loc)
        end function fun_pol
        
        ! function that  returns df/dtheta_V  = -1  - dlambda/dtheta_V.  It uses
        ! theta_F and zeta_E  (= zeta_V) from the parent function,  lam and dlam
        ! to contain the variable lambda or  its derivative, as well as L_s_loc,
        ! L_c_loc, theta_V_loc, zeta_V_loc, id, jd and kd.
        function dfun_pol(theta_E_in)
            character(*), parameter :: rout_name = 'dfun_pol'
            
            ! input / output
            real(dp) :: dfun_pol
            real(dp), intent(in) :: theta_E_in
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:)              ! trigonometric factor cosine for the inverse fourier transf.
            
            ! initialize dfun_pol
            dfun_pol = 0.0_dp
            
            ! set up local copies of theta_V_loc, zeta_V_loc
            theta_V_loc = theta_E_in
            zeta_V_loc = zeta_E(id,jd,kd)
            
            ! transform lambda from Fourier space to real space
            ! calculate the (co)sines
            ierr = calc_trigon_factors(theta_V_loc,zeta_V_loc,&
                &trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_c_loc,L_s_loc,trigon_factors_loc,dlam,[1,0])
            CHCKERR('')
            
            ! calculate the output function
            dfun_pol = -1._dp - dlam(1,1)
            
            ! deallocate trigoniometric factors
            deallocate(trigon_factors_loc)
        end function dfun_pol
    end function coord_F2E
    
    ! Converts  Equilibrium  coordinates  (r,theta,zeta)_E to  Flux  coordinates
    ! (r,theta,zeta)_F. Optionally,  two arrays  r_E_array and r_F_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, the normalized (wrt. 1) flux variables are used.
    integer function coord_E2F(r_E,theta_E,zeta_E,r_F,theta_F,zeta_F,&
        &r_E_array,r_F_array) result(ierr)
        use num_vars, only: eq_style, use_pol_flux_X, use_pol_flux_eq
        use utilities, only: interp_fun
        use eq_vars, only: flux_p_FD, flux_t_FD, max_flux_F, max_flux_eq_F
        
        character(*), parameter :: rout_name = 'coord_E2F'
        
        ! input / output
        real(dp), intent(in) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)           ! equilibrium coords.
        real(dp), intent(inout) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)        ! Flux (perturbation) coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), pointer :: flux(:), flux_eq(:)                                ! either pol. or tor. flux
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: id, jd, kd                                                   ! counters
        real(dp) :: r_F_factor, r_E_factor                                      ! mult. factors for r_F and r_E
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_theta = size(theta_E,1)
        n_zeta = size(theta_E,2)
        n_r = size(theta_E,3)
        
        ! tests
        if (size(theta_F,1).ne.size(zeta_F,1) .or. &
            &size(theta_F,2).ne.size(zeta_F,2) .or. &
            &size(theta_F,3).ne.size(zeta_F,3) .or. &
            &size(theta_F,3).ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(zeta_E,1) .or. n_zeta.ne.size(zeta_E,2) .or. &
            &n_r.ne.size(zeta_E,3) .or. n_r.ne.size(r_E)) then
            ierr = 1
            err_msg = 'theta_E, zeta_E and r_E need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_F_array).neqv.present(r_E_array)) then
            ierr = 1
            err_msg = 'both r_F_array and r_E_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! set up flux, flux_eq and multiplicative factors r_F_factor, r_E_factor
        if (present(r_F_array)) then
            flux => r_F_array
            flux_eq => r_E_array
            r_F_factor = 1._dp
            r_E_factor = 1._dp
        else
            if (use_pol_flux_X) then
                flux => flux_p_FD(:,0)
            else
                flux => flux_t_FD(:,0)
            end if
            if (use_pol_flux_eq) then
                flux_eq => flux_p_FD(:,0)
            else
                flux_eq => flux_t_FD(:,0)
            end if
            r_F_factor = max_flux_F
            r_E_factor = max_flux_eq_F
        end if
        
        ! convert normal position
        do kd = 1,n_r
            ierr = interp_fun(r_F(kd),flux,r_E(kd)*r_E_factor,flux_eq)
            CHCKERR('')
            r_F(kd) = r_F(kd)*r_F_factor
        end do
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = coord_E2F_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! trivial HELENA uses flux coordinates
                theta_F = theta_E
                zeta_F = zeta_E
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        integer function coord_E2F_VMEC() result(ierr)
            use VMEC, only: mpol, ntor, L_c, L_s
            use fourier_ops, only: calc_trigon_factors, fourier2real
            use eq_vars, only: grp_min_r_eq, grp_max_r_eq
            
            ! local variables
            real(dp), allocatable :: L_c_loc(:,:,:), L_s_loc(:,:,:)             ! local version of L_c and L_s
            real(dp), allocatable :: L_c_loc_ind(:,:), L_s_loc_ind(:,:)         ! individual versions of L_c_loc and L_s_loc
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:)              ! trigonometric factor cosine for the inverse fourier transf.
            real(dp) :: lam(1,1)                                                ! lambda (needs to be 2D matrix)
            real(dp) :: theta_V_loc(1,1), zeta_V_loc(1,1)                       ! local version of theta_V, zeta_V (needs to be 2D matrix)
            
            ! initialize ierr
            ierr = 0
            
            ! the toroidal angle is trivial
            zeta_F = - zeta_E                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_c and L_s
            allocate(L_c_loc(0:mpol-1,-ntor:ntor,1:1))
            allocate(L_s_loc(0:mpol-1,-ntor:ntor,1:1))
            allocate(L_c_loc_ind(0:mpol-1,-ntor:ntor))
            allocate(L_s_loc_ind(0:mpol-1,-ntor:ntor))
            
            ! loop over all normal points
            do kd = 1,n_r
                ! interpolate L_c and L_s at requested normal point r_E
                ierr = interp_fun(L_c_loc_ind,&
                    &L_c(:,:,grp_min_r_eq:grp_max_r_eq,0),&
                    &r_E(kd)*r_E_factor,flux_eq)
                CHCKERR('')
                ierr = interp_fun(L_s_loc_ind,&
                    &L_s(:,:,grp_min_r_eq:grp_max_r_eq,0),&
                    &r_E(kd)*r_E_factor,flux_eq)
                CHCKERR('')
                
                ! copy individual to array version
                L_c_loc(:,:,1) = L_c_loc_ind
                L_s_loc(:,:,1) = L_s_loc_ind
                
                ! loop over all angular points
                do jd = 1,n_zeta
                    do id = 1,n_theta
                        ! set up local copies of theta_V_loc, zeta_V_loc
                        theta_V_loc = theta_E(id,jd,kd)
                        zeta_V_loc = zeta_E(id,jd,kd)
                        
                        ! transform lambda from Fourier space to real space
                        ! calculate the (co)sines
                        ierr = calc_trigon_factors(theta_V_loc,zeta_V_loc,&
                            &trigon_factors_loc)
                        CHCKERR('')
                        
                        ! calculate lambda
                        ierr = fourier2real(L_c_loc,L_s_loc,trigon_factors_loc,&
                            &lam)
                        CHCKERR('')
                        
                        ! the poloidal angle has to be found as
                        !   theta_F = theta_E + lambda
                        theta_F(id,jd,kd) = theta_E(id,jd,kd) + lam(1,1)
                        
                        ! deallocate trigoniometric factors
                        deallocate(trigon_factors_loc)
                    end do
                end do
            end do
            
            ! deallocate variables
            deallocate(L_c_loc,L_s_loc)
            deallocate(L_c_loc_ind,L_s_loc_ind)
        end function coord_E2F_VMEC
    end function coord_E2F
    
    ! calculates X,Y  and Z on a  grid in the equilibrium  poloidal and toroidal
    ! angle theta and zeta, for every normal point in the equilibrium grid.
    ! The dimensions are (n_theta,n_zeta,n_r)
    ! If VMEC is the equilibrium  model, this routine also optionally calculates
    ! lambda on the grid, as this is  also needed some times. If HELENA is used,
    ! this variable is not allocated.
    ! Note: The variables X, Y, Z and optionally L have to be unallocated.
    integer function calc_XYZ_grid(r_E,theta_E,zeta_E,X,Y,Z,L) result(ierr)
        use num_vars, only: eq_style
        use utilities, only: interp_fun, round_with_tol
        
        character(*), parameter :: rout_name = 'calc_XYZ_grid'
        
        ! input / output
        real(dp), intent(in) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)           ! points at which to calculate the grid
        real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:)    ! X, Y and Z of grid
        real(dp), intent(inout), allocatable, optional :: L(:,:,:)              ! lambda of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (size(theta_E,1).ne.size(zeta_E,1) .or. &
            &size(theta_E,2).ne.size(zeta_E,2) .or. &
            &size(theta_E,3).ne.size(zeta_E,3) .or. &
            &size(theta_E,3).ne.size(r_E)) then
            ierr = 1
            err_msg =  'theta_E, zeta_E and r_E need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_XYZ_grid_VMEC(r_E,theta_E,zeta_E,X,Y,Z,L)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_XYZ_grid_HEL(r_E,theta_E,zeta_E,X,Y,Z)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! VMEC version
        integer function calc_XYZ_grid_VMEC(r_E,theta_E,zeta_E,X,Y,Z,L) &
            &result(ierr)
            use VMEC, only: R_c, R_s, Z_c, Z_s, L_c, L_s, mpol, ntor
            use fourier_ops, only: calc_trigon_factors, fourier2real
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_VMEC'
            
            ! input / output
            real(dp), intent(in) :: r_E(:)                                      ! r values in equilibrium range
            real(dp), intent(in) :: theta_E(:,:,:), zeta_E(:,:,:)               ! points at which to calculate the grid
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
            n_theta = size(theta_E,1)
            n_zeta = size(theta_E,2)
            n_r = size(theta_E,3)
            
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
                        ierr = interp_fun(R_c_int(m,n,kd),R_c(m,n,:,0),r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(R_s_int(m,n,kd),R_s(m,n,:,0),r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(Z_c_int(m,n,kd),Z_c(m,n,:,0),r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(Z_s_int(m,n,kd),Z_s(m,n,:,0),r_E(kd))
                        CHCKERR('')
                        if (present(L)) then
                            ierr = interp_fun(L_c_int(m,n,kd),L_c(m,n,:,0),&
                                &r_E(kd))
                            CHCKERR('')
                            ierr = interp_fun(L_s_int(m,n,kd),L_s(m,n,:,0),&
                                &r_E(kd))
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! inverse fourier transform with trigonometric factors
            do id = 1,n_theta
                ! calculate trigonometric factors
                ierr = calc_trigon_factors(theta_E(id,:,:),zeta_E(id,:,:),&
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
            X = R*cos(-zeta_E)
            Y = R*sin(-zeta_E)
            
            ! deallocate
            deallocate(R_c_int,R_s_int,Z_c_int,Z_s_int)
            if (present(L)) deallocate(L_c_int,L_s_int)
            deallocate(R)
        end function calc_XYZ_grid_VMEC
        
        ! HELENA version
        integer function calc_XYZ_grid_HEL(r_E,theta_E,zeta_E,X,Y,Z) &
            &result(ierr)
            use HELENA, only: R_H, Z_H, chi_H, ias, flux_H
            use eq_vars, only: n_r_eq
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_HEL'
            
            ! input / output
            real(dp), intent(in) :: r_E(:)                                      ! r values in equilibrium range
            real(dp), intent(in) :: theta_E(:,:,:), zeta_E(:,:,:)               ! points at which to calculate the grid
            real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), &
                &Z(:,:,:)                                                       ! X, Y and Z of grid
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            integer :: n_theta, n_zeta, n_r                                     ! dimensions of the grid
            real(dp), allocatable :: R_H_int(:), Z_H_int(:)                     ! R and Z at interpolated normal value
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            real(dp) :: theta_loc                                               ! local copy of theta_E
            
            ! initialize ierr
            ierr = 0
            
            ! set up n_theta, n_zeta and n_r
            n_theta = size(theta_E,1)
            n_zeta = size(theta_E,2)
            n_r = size(theta_E,3)
            
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
                    ierr = interp_fun(R_H_int(id),R_H(id,:),r_E(kd),&
                        &x=flux_H/flux_H(n_r_eq))
                    CHCKERR('')
                end do
                do id = 1,size(Z_H,1)
                    ierr = interp_fun(Z_H_int(id),Z_H(id,:),r_E(kd),&
                        &x=flux_H/flux_H(n_r_eq))
                    CHCKERR('')
                end do
                ! loop over toroidal points
                do jd = 1,n_zeta
                    ! interpolate at the requested poloidal points
                    ! Note: R_H  and Z_H are  not adapted to the  parallel grid,
                    ! but tabulated in the original HELENA poloidal grid
                    do id = 1,n_theta
                        theta_loc = theta_E(id,jd,kd)
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
                            ierr = interp_fun(R(id,jd,kd),R_H_int,&
                                &2*pi-theta_loc,x=chi_H)
                            CHCKERR('')
                            ierr = interp_fun(Z(id,jd,kd),-Z_H_int,&
                                &2*pi-theta_loc,x=chi_H)                        ! Z at second half of period is inverted
                            CHCKERR('')
                        else
                            ierr = interp_fun(R(id,jd,kd),R_H_int,theta_loc,&
                                &x=chi_H)
                            CHCKERR('')
                            ierr = interp_fun(Z(id,jd,kd),Z_H_int,theta_loc,&
                                &x=chi_H)
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! calculate X and Y, transforming cylindrical to cartesian
            X = R*cos(zeta_E)
            Y = R*sin(zeta_E)
            
            ! deallocate
            deallocate(R)
        end function calc_XYZ_grid_HEL
    end function calc_XYZ_grid

    ! calculate grid of equidistant points,  where optionally the last point can
    ! be excluded
    ! Note: input is given in units of pi
    integer function calc_eqd_grid(eqd_grid,n_ang,min_ang,max_ang,excl_last) &
        &result(ierr)
        character(*), parameter :: rout_name = 'eqd_grid'
        
        ! input and output
        real(dp), intent(inout) :: eqd_grid(:)                                  ! output
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
        eqd_grid = 0.0_dp
        
        ! There are (n_ang-1) pieces in the total interval but if excl_last, the
        ! last one is not included
        if (excl_last_loc) then
            delta_ang = (max_ang-min_ang)/(n_ang) * pi
        else
            delta_ang = (max_ang-min_ang)/(n_ang-1) * pi
        end if
        
        eqd_grid(1) = min_ang*pi
        do id = 2,n_ang
            eqd_grid(id) = eqd_grid(id-1) + delta_ang
        end do
    end function calc_eqd_grid
    
    ! plots the grid in real space
    ! The  variables theta  and  zeta follow  the magnetic  field  lines in  the
    ! equilibrium grid  and are tabulated along  the magnetic field lines  for a
    ! series  of flux  surfaces. The  start  and end  index of  the normal  flux
    ! surfaces in the equilibrium grid has  to be provided with min_r and max_r.
    ! Ideally, this routine should be run  with the full equilibrium grid normal
    ! extent, which implies for example, min_r = 0, max_r = 1.
    ! Note: This routine has its  own n_theta_plot and n_zeta_plot, since it has
    ! to be 3D also in the axisymmetric case.
    ! [MPI] Collective call
    integer function plot_grid_real(theta,zeta,min_r,max_r) &
        &result(ierr)
        use X_vars, only: n_par
        use num_vars, only: output_style, alpha_job_nr
        
        character(*), parameter :: rout_name = 'plot_grid_real'
        
        ! input / output
        real(dp), intent(in) :: theta(:,:)                                      ! theta in VMEC or HELENA coordinates
        real(dp), intent(in) :: zeta(:,:)                                       ! zeta in VMEC or HELENA coordinates
        real(dp), intent(in) :: min_r, max_r                                    ! minimum and maximum normal range
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=*), parameter :: anim_name = &
            &'Magnetic field in flux surfaces'                                  ! name of animation
        real(dp), allocatable :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)             ! X, Y and Z of surface in Axisymmetric coordinates
        real(dp), allocatable :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)             ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
        real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)            ! theta and zeta for plots
        real(dp), allocatable :: r_plot(:)                                      ! r for plots
        integer :: id                                                           ! counter
        integer :: loc_n_par, loc_n_r                                           ! local nr. of parallel and normal points
        integer :: n_theta_plot = 40                                            ! nr. of poloidal points in plot
        integer :: n_zeta_plot = 160                                            ! nr. of toroidal points in plot
        
        ! initialize ierr
        ierr = 0
        
        call writo('Plotting magnetic field and flux surfaces')
        
        call lvl_ud(1)
        
        ! some tests
        if (size(theta,1).ne.size(zeta,1) .or. &
            &size(theta,1).ne.n_par .or. &
            &size(theta,2).ne.size(zeta,2)) then
            err_msg = 'theta and zeta have to have the correct sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! 1. plot flux surfaces
        call writo('writing flux surfaces')
        
        ! initialize loc_n_par and loc_n_r
        loc_n_par = size(theta,1)
        loc_n_r = size(theta,2)
        
        ! initialize theta_plot, zeta_plot and r_plot
        allocate(theta_plot(n_theta_plot,n_zeta_plot,loc_n_r))
        do id = 1,n_theta_plot
            theta_plot(id,:,:) = pi+(id-1.0_dp)*2*pi/(n_theta_plot-1.0_dp)      ! better to start from pi for the plot
        end do
        allocate(zeta_plot(n_theta_plot,n_zeta_plot,loc_n_r))
        do id = 1,n_zeta_plot
            zeta_plot(:,id,:) = (id-1.0_dp)*2*pi/(n_zeta_plot-1.0_dp)
        end do
        allocate(r_plot(loc_n_r))
        do id = 1,loc_n_r
            r_plot(id) = min_r + (id-1._dp)/(loc_n_r-1)*(max_r-min_r)
        end do
        
        ! calculate X,Y and Z
        ierr = calc_XYZ_grid(r_plot,theta_plot,zeta_plot,X_1,Y_1,Z_1)
        CHCKERR('')
        
        ! deallocate
        deallocate(theta_plot,zeta_plot)
        
        ! 2. plot field lines
        call writo('writing field lines')
        
        ! initialize theta_plot and zeta_plot
        allocate(theta_plot(1,loc_n_par,loc_n_r),zeta_plot(1,loc_n_par,loc_n_r))
        do id = 1,loc_n_par
            theta_plot(1,id,:) = theta(id,:)
            zeta_plot(1,id,:) = zeta(id,:)
        end do
        
        ! calculate X,Y and Z
        ierr = calc_XYZ_grid(r_plot,theta_plot,zeta_plot,X_2,Y_2,Z_2)
        CHCKERR('')
        
        ! deallocate
        deallocate(theta_plot,zeta_plot)
        
        ! delegate to child routines
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                call plot_grid_real_GP(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)
            case(2)                                                             ! HDF5 output
                ierr = plot_grid_real_HDF5(x_1,x_2,y_1,y_2,z_1,z_2,anim_name)
                CHCKERR('')
            case default
                err_msg = 'No style associated with '//trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! deallocate
        deallocate(X_1,Y_1,Z_1)
        deallocate(X_2,Y_2,Z_2)
        
        call lvl_ud(-1)
        
        call writo('Done plotting magnetic field and flux surfaces')
    contains
        subroutine plot_grid_real_GP(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)         ! GNUPlot version
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: file_name(2), plot_title(2)            ! file name and title
            character(len=max_str_ln) :: draw_ops(2)                            ! individual plot command
            integer :: n_r                                                      ! total number of normal points
           
            ! set up n_r
            n_r = size(X_1,3)                                                   ! should be same for all other X_i, Y_i and Z_i
            
            ! user output
            call writo('Drawing animation with GNUPlot')
            call lvl_ud(1)
            
            ! set names
            plot_title(1) = 'Magnetic Flux Surface for alpha job '//&
                &trim(i2str(alpha_job_nr))
            file_name(1) = 'Flux_surfaces_'//trim(i2str(alpha_job_nr))//'.dat'
            plot_title(2) = 'Magnetic Field Line for alpha job '//&
                &trim(i2str(alpha_job_nr))
            file_name(2) = 'B_field_'//trim(i2str(alpha_job_nr))//'.dat'
            
            ! write flux surfaces
            call print_GP_3D(trim(plot_title(1)),trim(file_name(1)),&
                &Z_1(:,:,1:n_r-1),x=X_1(:,:,1:n_r-1),&
                &y=Y_1(:,:,1:n_r-1),draw=.false.)                               ! plot all but last points
            
            ! initialize n_r
            n_r = size(X_2,3)
            
            ! write magnetic field lines
            call print_GP_3D(trim(plot_title(2)),trim(file_name(2)),&
                &Z_2(:,:,2:n_r),x=X_2(:,:,2:n_r),y=Y_2(:,:,2:n_r),&
                &draw=.false.)                                                  ! plot all bit first points
            
            ! plot both output files
            call writo('creating GNUPlot animation')
            draw_ops(1) = 'linecolor rgb ''#d3d3d3'' linewidth 1'
            draw_ops(2) = 'linecolor rgb ''black'' linewidth 3'
            call draw_GP_animated(anim_name,file_name,n_r-1,.false.,&
                &delay=50,draw_ops=draw_ops)
            
            call lvl_ud(-1)
        end subroutine plot_grid_real_GP
        integer function plot_grid_real_HDF5(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)&
            & result(ierr)                                                      ! HDF5 version
            use HDF5_ops, only: open_HDF5_file, add_HDF5_item, &
                &print_HDF5_top, print_HDF5_geom, print_HDF5_3D_data_item, &
                &print_HDF5_grid, close_HDF5_file, &
                &XML_str_type, HDF5_file_type
            
            character(*), parameter :: rout_name = 'plot_grid_real_HDF5'
            
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: plot_title(2)                          ! plot title for flux surface and for field lines
            integer :: id                                                       ! counter
            type(HDF5_file_type) :: file_info                                   ! file info
            type(XML_str_type) :: grids(2)                                      ! grid in spatial collection (flux surface, field line)
            type(XML_str_type), allocatable :: space_col_grids(:)               ! grid with space collection at different times
            type(XML_str_type) :: time_col_grid                                 ! grid with time collection
            type(XML_str_type) :: top                                           ! topology
            type(XML_str_type) :: XYZ(3)                                        ! data items for geometry
            type(XML_str_type) :: geom                                          ! geometry
            integer :: loc_dim(3,2)                                             ! local dimensions (flux surfaces, field lines)
            integer :: n_r                                                      ! nr. of normal points
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Drawing animation with HDF5')
            call lvl_ud(1)
            
            ! set up loc_dim and n_r
            loc_dim(:,1) = [size(X_1,1),size(X_1,2),1]
            loc_dim(:,2) = [size(X_2,1),size(X_2,2),1]
            n_r = size(X_1,3)                                                   ! should be same for all other X_i, Y_i and Z_i
            
            ! set up plot titles
            plot_title(1) = 'Magnetic Flux Surface for alpha job '//&
                &trim(i2str(alpha_job_nr))
            plot_title(2) = 'Magnetic Field Line for alpha job '//&
                &trim(i2str(alpha_job_nr))
            
            ! open HDF5 file
            ierr = open_HDF5_file(file_info,'field_lines_in_flux_surfaces_'//&
                &trim(i2str(alpha_job_nr)),description=anim_name,&
                &ind_plot=.true.)
            CHCKERR('')
            
            ! create grid for time collection
            allocate(space_col_grids(n_r-1))
            
            ! loop over all normal points
            do id = 1,n_r-1
                ! A. start with flux surface
                
                ! print topology
                call print_HDF5_top(top,2,loc_dim(:,1))
                
                ! print data item for X
                ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
                    &'X_surf_'//trim(i2str(id)),X_1(:,:,id:id),loc_dim(:,1))
                CHCKERR('')
                
                ! print data item for Y
                ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Y_surf_'//trim(i2str(id)),Y_1(:,:,id:id),loc_dim(:,1))
                CHCKERR('')
                
                ! print data item for Z
                ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
                    &'Z_surf_'//trim(i2str(id)),Z_1(:,:,id:id),loc_dim(:,1))
                CHCKERR('')
                
                ! print geometry with X, Y and Z data item
                call print_HDF5_geom(geom,2,XYZ,.true.)
                
                ! create a grid with the topology and the geometry
                ierr = print_HDF5_grid(grids(1),plot_title(1),1,&
                    &grid_top=top,grid_geom=geom,reset=.true.)
                CHCKERR('')
                
                ! B. end with magnetic field
                
                ! print topology
                call print_HDF5_top(top,2,loc_dim(:,2))
                
                ! print data item for X
                ierr = print_HDF5_3D_data_item(XYZ(1),file_info,&
                    &'X_field_'//trim(i2str(id)),X_2(:,:,id+1:id+1),loc_dim(:,2))
                CHCKERR('')
                
                ! print data item for Y
                ierr = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Y_field_'//trim(i2str(id)),Y_2(:,:,id+1:id+1),loc_dim(:,2))
                CHCKERR('')
                
                ! print data item for Z
                ierr = print_HDF5_3D_data_item(XYZ(3),file_info,&
                    &'Z_field_'//trim(i2str(id)),Z_2(:,:,id+1:id+1),loc_dim(:,2))
                CHCKERR('')
                
                ! print geometry with X, Y and Z data item
                call print_HDF5_geom(geom,2,XYZ,.true.)
                
                ! create a grid with the topology and the geometry
                ierr = print_HDF5_grid(grids(2),plot_title(2),1,&
                    &grid_top=top,grid_geom=geom,reset=.true.)
                CHCKERR('')
                
                ! C. merge flux surface and magnetic field in to spatial collection
                ierr = print_HDF5_grid(space_col_grids(id),'spatial collection',&
                    &3,grid_time=id*1._dp,grid_grids=grids,reset=.true.)
                CHCKERR('')
            end do
            
            ! create grid collection from individual grids and reset them
            ierr = print_HDF5_grid(time_col_grid,'time collection',2,&
                &grid_grids=space_col_grids,reset=.true.)
            CHCKERR('')
            
            ! add collection grid to HDF5 file and reset it
            call add_HDF5_item(file_info,time_col_grid,reset=.true.)
            
            ! close HDF5 file
            ierr = close_HDF5_file(file_info)
            CHCKERR('')
            
            call lvl_ud(-1)
        end function plot_grid_real_HDF5
    end function plot_grid_real
end module coord_ops
