!------------------------------------------------------------------------------!
!   Operations that have to do with the grids and different coordinate systems !
!------------------------------------------------------------------------------!
module grid_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type

    implicit none
    private
    public coord_F2E, coord_E2F, calc_XYZ_grid, calc_eqd_grid, calc_ang_grid, &
        &plot_grid_real, trim_grid, extend_grid
#if ldebug
    public debug_calc_ang_grid
#endif
    
    interface coord_F2E
        module procedure coord_F2E_r, coord_F2E_rtz
    end interface
    interface coord_E2F
        module procedure coord_E2F_r, coord_E2F_rtz
    end interface
    interface calc_eqd_grid
        module procedure calc_eqd_grid_1D, calc_eqd_grid_3D
    end interface
    
    ! global variables
#if ldebug
    logical :: debug_calc_ang_grid = .false.                                    ! plot debug information for calc_ang_grid
#endif
    
contains
    ! Calculate the angular equilibrium grid  in the Flux (perturbation) angular
    ! variable, but tabulated in the Equilibrium (!) normal variable.
    ! The variable  use_pol_flux determines  whether theta (.true.)  or zeta
    ! (.false.) is used as the parallel variable.
    integer function calc_ang_grid(grid_eq,eq,alpha) result(ierr)
        use num_vars, only: grid_style, use_pol_flux_F, use_pol_flux_E, &
            &eq_style, tol_NR
        use grid_vars, only: min_par_X, max_par_X
        use eq_vars, only: max_flux_p_E, max_flux_t_E
        
        character(*), parameter :: rout_name = 'calc_ang_grid'
        
        ! input / output
        type(grid_type), intent(inout) :: grid_eq                               ! equilibrium grid of which to calculate angular part
        type(eq_type), intent(in) :: eq                                         ! equilibrium containing the angular grid
        real(dp), intent(in) :: alpha                                           ! field line coordinate
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: r_E_loc(:)                                     ! flux in Equilibrium coords.
        real(dp), pointer :: flux_F(:), flux_E(:)                               ! flux that the F and E use as normal coord.
        real(dp) :: r_F_factor, r_E_factor                                      ! mult. factors for r_F and r_E
        integer :: pmone                                                        ! plus or minus one
        integer :: kd                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (grid_eq%n(2).ne.1) then
            ierr = 1
            err_msg = 'The angular grid_eq must be for a single field line, so &
                &it can only have one geodesic value.'
            CHCKERR(err_msg)
        end if
        
        ! calculate the grid
        ! grid style
        !   1:  grid along magnetic field lines
        select case (grid_style)
            ! grid along the magnetic field lines
            case (1)
                ! set up flux_E and plus minus one
                ! Note: this routine is similar  to calc_grp_r, but that routine
                ! cannot  be used  here because  the FD  quantities are  not yet
                ! defined.
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        if (use_pol_flux_E) then                                ! normalized poloidal flux
                            flux_E => eq%flux_p_E(:,0)
                            r_E_factor = max_flux_p_E
                        else                                                    ! normalized toroidal flux
                            flux_E => eq%flux_t_E(:,0)
                            r_E_factor = max_flux_t_E
                        end if
                        pmone = -1                                              ! conversion VMEC LH -> RH coord. system
                    case (2)                                                    ! HELENA
                        flux_E => eq%flux_p_E(:,0)
                        r_E_factor = 1._dp
                        pmone = 1
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! set up flux_F
                if (use_pol_flux_F) then                                        ! poloidal flux / 2pi
                    flux_F => eq%flux_p_E(:,0)
                    r_F_factor = 2*pi
                else                                                            ! toroidal flux / 2pi
                    flux_F => eq%flux_t_E(:,0)
                    r_F_factor = pmone*2*pi                                     ! possible conversion VMEC LH -> RH coord. system
                end if
                
                ! set up parallel angle in  Flux coordinates on equidistant grid
                ! and use this to calculate the other angle as well
                if (use_pol_flux_F) then                                        ! parallel angle theta
                    ierr = calc_eqd_grid(grid_eq%theta_F,min_par_X*pi,&
                        &max_par_X*pi,1)                                        ! first index corresponds to parallel angle
                    CHCKERR('')
                    do kd = 1,grid_eq%grp_n_r
                        grid_eq%zeta_F(:,:,kd) = pmone*eq%q_saf_E(kd,0)*&
                            &grid_eq%theta_F(:,:,kd)
                    end do
                    grid_eq%zeta_F = grid_eq%zeta_F + alpha
                else                                                            ! parallel angle zeta
                    ierr = calc_eqd_grid(grid_eq%zeta_F,min_par_X*pi,&
                        &max_par_X*pi,1)                                        ! first index corresponds to parallel angle
                    CHCKERR('')
                    do kd = 1,grid_eq%grp_n_r
                        grid_eq%theta_F(:,:,kd) = pmone*eq%rot_t_E(kd,0)*&
                            &grid_eq%zeta_F(:,:,kd)
                    end do
                    grid_eq%theta_F = grid_eq%theta_F - alpha
                end if
                
                ! allocate local r_E
                allocate(r_E_loc(size(flux_F)))
                
                ! convert  Flux  coordinates  to  Equilibrium  coordinates  (use
                ! custom flux_E and  flux_F because the Flux  quantities are not
                ! yet calculated)
                ierr = coord_F2E(eq,grid_eq,&
                    &grid_eq%grp_r_F,grid_eq%theta_F,grid_eq%zeta_F,&
                    &r_E_loc,grid_eq%theta_E,grid_eq%zeta_E,&
                    &r_F_array=flux_F/r_F_factor,r_E_array=flux_E/r_E_factor)
                CHCKERR('')
                
                ! test whether r_E_loc indeed corresponds to grp_r_E of eq. grid
                if (maxval(abs(grid_eq%grp_r_E-r_E_loc)).gt.10*tol_NR) then
                    ierr = 1
                    err_msg = 'grp_r_E of equilibrium grid is not recovered'
                    CHCKERR(err_msg)
                end if
                
                ! deallocate local variables
                deallocate(r_E_loc)
                nullify(flux_F,flux_E)
                
#if ldebug
                if (debug_calc_ang_grid) then
                    call writo('Plotting theta_E, theta_F, zeta_E and zeta_F')
                    call print_HDF5('TEST_theta_E','TEST_theta_E',grid_eq%theta_E)
                    call print_HDF5('TEST_theta_F','TEST_theta_F',grid_eq%theta_F)
                    call print_HDF5('TEST_zeta_E','TEST_zeta_E',grid_eq%zeta_E)
                    call print_HDF5('TEST_zeta_F','TEST_zeta_F',grid_eq%zeta_F)
                end if
#endif
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
    ! Standard, for E, the poloidal or  toroidal normalized flux is used and for
    ! F, the poloidal or toroidal flux in F coordinates, divided by 2pi.
    integer function coord_F2E_rtz(eq,grid_eq,r_F,theta_F,zeta_F,r_E,&
        &theta_E,zeta_E,r_F_array,r_E_array) result(ierr)                       ! version with r, theta and zeta
        use num_vars, only: tol_NR, eq_style
        use fourier_ops, only: fourier2real, calc_trigon_factors
        use VMEC, only: L_c, L_s
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_F2E_rtz'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal group limits)
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux coords.
        real(dp), intent(inout) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)        ! Equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r, n_par, n_geo                                            ! dimensions of the grid
        real(dp), allocatable :: L_c_loc(:,:,:)                                 ! local version of L_c
        real(dp), allocatable :: L_s_loc(:,:,:)                                 ! local version of L_s
        real(dp) :: theta_V_loc(1,1,1), zeta_V_loc(1,1,1)                       ! local version of theta_V, zeta_V (needs to be 3D matrix)
        real(dp) :: lam(1,1,1)                                                  ! lambda (needs to be 3D matrix)
        real(dp) :: dlam(1,1,1)                                                 ! angular derivative of lambda (needs to be 3D matrix)
        integer :: id, jd, kd                                                   ! counters
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_par = size(theta_E,1)
        n_geo = size(theta_E,2)
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
        if (n_par.ne.size(zeta_E,1) .or. n_geo.ne.size(zeta_E,2) .or. &
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
        
        ! convert normal position
        ierr = coord_F2E_r(grid_eq,eq,r_F,r_E,r_F_array,r_E_array)
        CHCKERR('')
        
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
            
            character(*), parameter :: rout_name = 'coord_F2E_VMEC'
            
            ! local variables
            real(dp), allocatable :: L_c_loc_ind(:,:), L_s_loc_ind(:,:)         ! individual versions of L_c_loc and L_s_loc
            real(dp), allocatable :: grp_r_F(:)                                 ! grp_r in F coords.
            
            ! initialize ierr
            ierr = 0
            
            ! set up grp_r_F
            if (present(r_F_array)) then
                allocate(grp_r_F(grid_eq%grp_n_r))
                grp_r_F = r_F_array
            else
                ierr = calc_grp_r(grid_eq,eq,grp_r_F,style=2)
                CHCKERR('')
            end if
            
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
                    &L_c(:,:,grid_eq%i_min:grid_eq%i_max,0),r_F(kd),grp_r_F)
                CHCKERR('')
                ierr = interp_fun(L_s_loc_ind,&
                    &L_s(:,:,grid_eq%i_min:grid_eq%i_max,0),r_F(kd),grp_r_F)
                CHCKERR('')
                
                ! copy individual to array version
                L_c_loc(:,:,1) = L_c_loc_ind
                L_s_loc(:,:,1) = L_s_loc_ind
                
                ! the poloidal angle has to be found as the zero of
                !   f = theta_F - theta_E - lambda
                ! loop over all angular points
                do jd = 1,n_geo
                    do id = 1,n_par
                        ! calculate zero of f
                        ierr = calc_zero_NR(theta_E(id,jd,kd),fun_pol,&
                            &dfun_pol,theta_F(id,jd,kd))                        ! use theta_F as guess for theta_E
                        CHCKERR('')
                        
                        ! do a check whether the result is indeed zero
                        if (abs(fun_pol(theta_E(id,jd,kd))).gt.tol_NR*100) then
                            err_msg = 'In coord_F2E_VMEC, calculating f as a &
                                &check, using the theta_E that is the solution&
                                & of f = 0, yields a calculated f that &
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
        ! and  zeta_E (=  zeta_V) from  the parent  function, lam  and dlam  to
        ! contain the  variable lambda  or its derivative,  as well  as L_s_loc,
        ! L_c_loc, theta_V_loc, zeta_V_loc, id, jd and kd.
        function fun_pol(theta_E_in)
            character(*), parameter :: rout_name = 'fun_pol'
            
            ! input / output
            real(dp) :: fun_pol
            real(dp), intent(in) :: theta_E_in
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:,:)            ! trigonometric factor cosine for the inverse fourier transf.
            
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
            fun_pol = theta_F(id,jd,kd) - theta_E_in - lam(1,1,1)
            
            ! deallocate trigoniometric factors
            deallocate(trigon_factors_loc)
        end function fun_pol
        
        ! function that  returns df/dtheta_V  = -1  - dlambda/dtheta_V.  It uses
        ! theta_F and zeta_E (= zeta_V) from  the parent function, lam and dlam
        ! to contain the variable lambda or  its derivative, as well as L_s_loc,
        ! L_c_loc, theta_V_loc, zeta_V_loc, id, jd and kd.
        function dfun_pol(theta_E_in)
            character(*), parameter :: rout_name = 'dfun_pol'
            
            ! input / output
            real(dp) :: dfun_pol
            real(dp), intent(in) :: theta_E_in
            
            ! local variables
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:,:)            ! trigonometric factor cosine for the inverse fourier transf.
            
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
            dfun_pol = -1._dp - dlam(1,1,1)
            
            ! deallocate trigoniometric factors
            deallocate(trigon_factors_loc)
        end function dfun_pol
    end function coord_F2E_rtz
    integer function coord_F2E_r(grid_eq,eq,r_F,r_E,r_F_array,r_E_array) &
        &result(ierr)                                                           ! version with only r
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_F2E_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal group limits)
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        real(dp), intent(in) :: r_F(:)                                          ! perturbation coords.
        real(dp), intent(inout) :: r_E(:)                                       ! Equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r                                                          ! dimension of the grid
        integer :: kd                                                           ! counter
        real(dp), allocatable :: grp_r_E(:), grp_r_F(:)                         ! grp_r in E and F coords.
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_r = size(r_F)
        
        ! tests
        if (n_r.ne.size(r_E)) then
            ierr = 1
            err_msg = 'r_F and r_E need to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_F_array).neqv.present(r_E_array)) then
            ierr = 1
            err_msg = 'both r_F_array and r_E_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! set up grp_r_E and grp_r_F
        if (present(r_F_array)) then
            allocate(grp_r_E(grid_eq%grp_n_r))
            allocate(grp_r_F(grid_eq%grp_n_r))
            grp_r_E = r_E_array
            grp_r_F = r_F_array
        else
            ierr = calc_grp_r(grid_eq,eq,grp_r_E,style=1)
            CHCKERR('')
            ierr = calc_grp_r(grid_eq,eq,grp_r_F,style=2)
            CHCKERR('')
        end if
        
        ! convert normal position
        do kd = 1,n_r
            ierr = interp_fun(r_E(kd),grp_r_E,r_F(kd),grp_r_F)
            CHCKERR('')
            r_E(kd) = r_E(kd)
        end do
    end function coord_F2E_r
    
    ! Converts  Equilibrium  coordinates  (r,theta,zeta)_E to  Flux  coordinates
    ! (r,theta,zeta)_F. Optionally,  two arrays  r_E_array and r_F_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, for E, the poloidal or  toroidal normalized flux is used and for
    ! F, the poloidal or toroidal flux in F coordinates, divided by 2pi.
    integer function coord_E2F_rtz(eq,grid_eq,r_E,theta_E,zeta_E,r_F,&
        &theta_F,zeta_F,r_E_array,r_F_array) result(ierr)                       ! version with r, theta and zeta
        use num_vars, only: eq_style
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_E2F_rtz'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal group limits)
        real(dp), intent(in) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)           ! Equilibrium coords.
        real(dp), intent(inout) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)        ! Flux coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r, n_par, n_geo                                            ! dimensions of the grid
        integer :: kd                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_par = size(theta_E,1)
        n_geo = size(theta_E,2)
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
        if (n_par.ne.size(zeta_E,1) .or. n_geo.ne.size(zeta_E,2) .or. &
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
        
        ! convert normal position
        ierr = coord_E2F_r(grid_eq,eq,r_E,r_F,r_E_array,r_F_array)
        CHCKERR('')
        
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
            
            character(*), parameter :: rout_name = 'coord_E2F_VMEC'
            
            ! local variables
            real(dp), allocatable :: L_c_loc(:,:,:), L_s_loc(:,:,:)             ! local version of L_c and L_s
            real(dp), allocatable :: trigon_factors_loc(:,:,:,:,:,:)            ! trigonometric factor cosine for the inverse fourier transf.
            real(dp), allocatable :: lam(:,:,:)                                 ! lambda
            real(dp), allocatable :: grp_r_E(:)                                 ! grp_r in E coords.
            
            ! initialize ierr
            ierr = 0
            
            ! set up grp_r_E
            if (present(r_E_array)) then
                allocate(grp_r_E(grid_eq%grp_n_r))
                grp_r_E = r_E_array
            else
                ierr = calc_grp_r(grid_eq,eq,grp_r_E,style=1)
                CHCKERR('')
            end if
            
            ! the toroidal angle is trivial
            zeta_F = - zeta_E                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_c and L_s and lambda
            allocate(L_c_loc(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(L_s_loc(0:mpol-1,-ntor:ntor,1:n_r))
            allocate(lam(n_par,n_geo,n_r))
            
            ! interpolate L_c and L_s at requested normal point r_E
            ! loop over all normal points
            do kd = 1,n_r
                ierr = interp_fun(L_c_loc(:,:,kd),&
                    &L_c(:,:,grid_eq%i_min:grid_eq%i_max,0),r_E(kd),grp_r_E)
                CHCKERR('')
                ierr = interp_fun(L_s_loc(:,:,kd),&
                    &L_s(:,:,grid_eq%i_min:grid_eq%i_max,0),r_E(kd),grp_r_E)
                CHCKERR('')
            end do
            
            ! calculate the (co)sines to transform lambda to real space
            ierr = calc_trigon_factors(theta_E,zeta_E,trigon_factors_loc)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_c_loc,L_s_loc,trigon_factors_loc,lam)
            CHCKERR('')
            
            ! the poloidal angle has to be found as
            !   theta_F = theta_E + lambda
            theta_F = theta_E + lam
            
            ! deallocate variables
            deallocate(trigon_factors_loc)
            deallocate(L_c_loc,L_s_loc)
            deallocate(lam)
        end function coord_E2F_VMEC
    end function coord_E2F_rtz
    integer function coord_E2F_r(grid_eq,eq,r_E,r_F,r_E_array,r_F_array) &
        &result(ierr)                                                           ! version with only r
        use utilities, only: interp_fun
        
        character(*), parameter :: rout_name = 'coord_E2F_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal group limits)
        type(eq_type), intent(in) :: eq                                         ! equilibrium in which to convert variables
        real(dp), intent(in) :: r_E(:)                                          ! Equilibrium coords.
        real(dp), intent(inout) :: r_F(:)                                       ! Flux coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r                                                          ! dimension of the grid
        integer :: kd                                                           ! counter
        real(dp), allocatable :: grp_r_E(:), grp_r_F(:)                         ! grp_r in E and F coords.
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        n_r = size(r_E)
        
        ! tests
        if (n_r.ne.size(r_F)) then
            ierr = 1
            err_msg = 'r_E and r_F need to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(r_E_array).neqv.present(r_F_array)) then
            ierr = 1
            err_msg = 'both r_E_array and r_F_array have to be provided'
            CHCKERR(err_msg)
        end if
        
        ! set up grp_r_E and grp_r_F
        if (present(r_F_array)) then
            allocate(grp_r_E(grid_eq%grp_n_r))
            allocate(grp_r_F(grid_eq%grp_n_r))
            grp_r_E = r_E_array
            grp_r_F = r_F_array
        else
            ierr = calc_grp_r(grid_eq,eq,grp_r_E,style=1)
            CHCKERR('')
            ierr = calc_grp_r(grid_eq,eq,grp_r_F,style=2)
            CHCKERR('')
        end if
        
        ! convert normal position
        do kd = 1,n_r
            ierr = interp_fun(r_F(kd),grp_r_F,r_E(kd),grp_r_E)
            CHCKERR('')
        end do
    end function coord_E2F_r
    
    ! Calculates X,Y  and Z on a  grid, which should have  the group equilibrium
    ! grid angles  set up in E(quilibrium)  coordinates. The total r_E,  and the
    ! F(lux) variables are ignored.
    ! If VMEC is the equilibrium  model, this routine also optionally calculates
    ! lambda on the grid, as this is  also needed some times. If HELENA is used,
    ! this variable is not used.
    ! X, Y, Z and  optionally L need to have the correct  dimensions if they are
    ! allocated.  They  can  also  be  passed unallocated,  in  which  case  the
    ! allocation happens automatically.
    integer function calc_XYZ_grid(grid,X,Y,Z,L) result(ierr)
        use num_vars, only: eq_style
        use utilities, only: interp_fun, round_with_tol
        
        character(*), parameter :: rout_name = 'calc_XYZ_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid for which to calculate X, Y, Z and optionally L
        real(dp), intent(inout), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:)    ! X, Y and Z of grid
        real(dp), intent(inout), allocatable, optional :: L(:,:,:)              ! lambda of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (allocated(X)) then
            if (size(X,1).ne.grid%n(1) .or. size(X,2).ne.grid%n(2) .or. &
                &size(X,3).ne.grid%grp_n_r) then
                ierr = 1
                err_msg =  'X needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(X(grid%n(1),grid%n(2),grid%grp_n_r))
        end if
        if (allocated(Y)) then
            if (size(Y,1).ne.grid%n(1) .or. size(Y,2).ne.grid%n(2) .or. &
                &size(Y,3).ne.grid%grp_n_r) then
                ierr = 1
                err_msg =  'Y needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(Y(grid%n(1),grid%n(2),grid%grp_n_r))
        end if
        if (allocated(Z)) then
            if (size(Z,1).ne.grid%n(1) .or. size(Z,2).ne.grid%n(2) .or. &
                &size(Z,3).ne.grid%grp_n_r) then
                ierr = 1
                err_msg =  'Z needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(Z(grid%n(1),grid%n(2),grid%grp_n_r))
        end if
        if (present(L)) then
            if (allocated(L)) then
                if (size(L,1).ne.grid%n(1) .or. size(L,2).ne.grid%n(2) .or. &
                    &size(L,3).ne.grid%grp_n_r) then
                    ierr = 1
                    err_msg =  'L needs to have the correct dimensions'
                    CHCKERR(err_msg)
                end if
            else
                allocate(L(grid%n(1),grid%n(2),grid%grp_n_r))
            end if
        end if
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_XYZ_grid_VMEC(grid,X,Y,Z,L)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_XYZ_grid_HEL(grid,X,Y,Z)
                CHCKERR('')
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        ! VMEC version
        integer function calc_XYZ_grid_VMEC(grid,X,Y,Z,L) result(ierr)
            use VMEC, only: R_c, R_s, Z_c, Z_s, L_c, L_s, mpol, ntor
            use fourier_ops, only: calc_trigon_factors, fourier2real
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_VMEC'
            
            ! input / output
            type(grid_type), intent(in) :: grid                                 ! grid for which to calculate X, Y, Z and optionally L
            real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)             ! X, Y and Z of grid
            real(dp), intent(inout), optional :: L(:,:,:)                       ! lambda of grid
            
            ! local variables
            integer :: kd, m, n                                                 ! counters
            real(dp), allocatable :: R_c_int(:,:,:), R_s_int(:,:,:)             ! interpolated version of R_c and R_s
            real(dp), allocatable :: Z_c_int(:,:,:), Z_s_int(:,:,:)             ! interpolated version of Z_c and Z_s
            real(dp), allocatable :: L_c_int(:,:,:), L_s_int(:,:,:)             ! interpolated version of L_c and L_s
            real(dp), allocatable :: trigon_factors(:,:,:,:,:,:)                ! trigonometric factor cosine for the inverse fourier transf.
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            
            ! initialize ierr
            ierr = 0
            
            ! set up interpolated R_c_int, ..
            allocate(R_c_int(0:mpol-1,-ntor:ntor,1:grid%grp_n_r))
            allocate(R_s_int(0:mpol-1,-ntor:ntor,1:grid%grp_n_r))
            allocate(Z_c_int(0:mpol-1,-ntor:ntor,1:grid%grp_n_r))
            allocate(Z_s_int(0:mpol-1,-ntor:ntor,1:grid%grp_n_r))
            if (present(L)) then
                allocate(L_c_int(0:mpol-1,-ntor:ntor,1:grid%grp_n_r))
                allocate(L_s_int(0:mpol-1,-ntor:ntor,1:grid%grp_n_r))
            end if
            
            ! interpolate VMEC tables for every requested normal point
            do kd = 1,grid%grp_n_r                                              ! loop over all group normal points
                ! interpolate the VMEC tables in normal direction
                ! Note: VMEC uses an equidistant grid, here normalized from 0 to
                ! 1
                do n = -ntor,ntor
                    do m = 0,mpol-1
                        ierr = interp_fun(R_c_int(m,n,kd),R_c(m,n,:,0),&
                            &grid%grp_r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(R_s_int(m,n,kd),R_s(m,n,:,0),&
                            &grid%grp_r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(Z_c_int(m,n,kd),Z_c(m,n,:,0),&
                            &grid%grp_r_E(kd))
                        CHCKERR('')
                        ierr = interp_fun(Z_s_int(m,n,kd),Z_s(m,n,:,0),&
                            &grid%grp_r_E(kd))
                        CHCKERR('')
                        if (present(L)) then
                            ierr = interp_fun(L_c_int(m,n,kd),L_c(m,n,:,0),&
                                &grid%grp_r_E(kd))
                            CHCKERR('')
                            ierr = interp_fun(L_s_int(m,n,kd),L_s(m,n,:,0),&
                                &grid%grp_r_E(kd))
                            CHCKERR('')
                        end if
                    end do
                end do
            end do
            
            ! calculate trigonometric factors
            ierr = calc_trigon_factors(grid%theta_E,grid%zeta_E,trigon_factors)
            CHCKERR('')
            
            ! allocate R
            allocate(R(grid%n(1),grid%n(2),grid%grp_n_r))
            
            ! inverse fourier transform with trigonometric factors
            ierr = fourier2real(R_c_int,R_s_int,trigon_factors,R,[0,0])
            CHCKERR('')
            ierr = fourier2real(Z_c_int,Z_s_int,trigon_factors,Z,[0,0])
            CHCKERR('')
            if (present(L)) then
                ierr = fourier2real(L_c_int,L_s_int,trigon_factors,L,[0,0])
                CHCKERR('')
            end if
            
            ! deallocate
            deallocate(trigon_factors)
            
            ! transform cylindrical to cartesian
            ! (the geometrical zeta is the inverse of VMEC zeta)
            X = R*cos(-grid%zeta_E)
            Y = R*sin(-grid%zeta_E)
            
            ! deallocate
            deallocate(R_c_int,R_s_int,Z_c_int,Z_s_int)
            if (present(L)) deallocate(L_c_int,L_s_int)
            deallocate(R)
        end function calc_XYZ_grid_VMEC
        
        ! HELENA version
        integer function calc_XYZ_grid_HEL(grid,X,Y,Z) result(ierr)
            use HELENA, only: R_H, Z_H, chi_H, ias, flux_H
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_HEL'
            
            ! input / output
            type(grid_type), intent(in) :: grid                                 ! grid for which to calculate X, Y, Z and optionally L
            real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)             ! X, Y and Z of grid
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            real(dp), allocatable :: R_H_int(:), Z_H_int(:)                     ! R and Z at interpolated normal value
            real(dp), allocatable :: R(:,:,:)                                   ! R in Cylindrical coordinates
            real(dp) :: theta_loc                                               ! local copy of theta_E
            
            ! initialize ierr
            ierr = 0
            
            ! allocate R
            allocate(R(grid%n(1),grid%n(2),grid%grp_n_r))
            
            ! set up interpolated R and Z
            allocate(R_H_int(size(R_H,1)),Z_H_int(size(Z_H,1)))
            
            ! interpolate HELENA output  R_H and Z_H for  every requested normal
            ! point
            ! Note:  R_H and  Z_H  are not  adapted to  the  parallel grid,  but
            ! tabulated in the original HELENA poloidal grid
            ! interpolate the HELENA tables in normal direction
            do kd = 1,grid%grp_n_r                                              ! loop over all normal points
                do id = 1,size(R_H,1)
                    ierr = interp_fun(R_H_int(id),R_H(id,:),grid%grp_r_E(kd),&
                        &x=flux_H)
                    CHCKERR('')
                end do
                do id = 1,size(Z_H,1)
                    ierr = interp_fun(Z_H_int(id),Z_H(id,:),grid%grp_r_E(kd),&
                        &x=flux_H)
                    CHCKERR('')
                end do
                ! loop over toroidal points
                do jd = 1,grid%n(2)
                    ! interpolate at the requested poloidal points
                    do id = 1,grid%n(1)
                        theta_loc = grid%theta_E(id,jd,kd)
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
            X = R*cos(grid%zeta_E)
            Y = R*sin(grid%zeta_E)
            
            ! deallocate
            deallocate(R)
        end function calc_XYZ_grid_HEL
    end function calc_XYZ_grid

    ! calculate grid of equidistant points,  where optionally the last point can
    ! be excluded
    integer function calc_eqd_grid_3D(var,min_grid,max_grid,grid_dim,&
        &excl_last) result(ierr)                                                ! 3D version
        character(*), parameter :: rout_name = 'calc_eqd_grid_3D'
        
        ! input and output
        real(dp), intent(inout) :: var(:,:,:)                                   ! output
        real(dp), intent(in) :: min_grid, max_grid                              ! min. and max. of angles [pi]
        integer, intent(in) :: grid_dim                                         ! in which dimension to create the grid
        logical, intent(in), optional :: excl_last                              ! .true. if last point excluded
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: excl_last_loc                                                ! local copy of excl_last
        integer :: grid_size                                                    ! nr. of points
        integer :: grid_size_mod                                                ! grid_size or grid_size + 1
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (grid_dim.lt.1 .or. grid_dim.gt.3) then
            ierr = 1
            err_msg = 'grid_dim has to point to a dimension going from 1 to 3'
            CHCKERR(err_msg)
        end if
        
        ! set up grid_size
        grid_size = size(var,grid_dim)
        
        ! test some values
        if (grid_size.lt.1) then
            err_msg = 'The angular array has to have a length of at &
                &least 1'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (min_grid.gt.max_grid) then
            err_msg = 'The minimum angle has to be smaller than the &
                &maximum angle'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set up local excl_last
        excl_last_loc = .false.
        if (present(excl_last)) excl_last_loc = excl_last
        
        ! add 1 to modified grid size if last point is to be excluded
        if (excl_last_loc) then
            grid_size_mod = grid_size + 1
        else
            grid_size_mod = grid_size
        end if
        
        ! initialize output vector
        var = 0.0_dp
        
        ! calculate grid points
        if (grid_size.eq.1) then
            var = min_grid                                                      ! = max_grid
        else
            if (grid_dim.eq.1) then
                do id = 1,grid_size
                    var(id,:,:) = min_grid + &
                        &(max_grid-min_grid)*(id-1)/(grid_size_mod-1)
                end do
            else if (grid_dim.eq.2) then
                do id = 1,grid_size
                    var(:,id,:) = min_grid + &
                        &(max_grid-min_grid)*(id-1)/(grid_size_mod-1)
                end do
            else
                do id = 1,grid_size
                    var(:,:,id) = min_grid + &
                        &(max_grid-min_grid)*(id-1)/(grid_size_mod-1)
                end do
            end if
        end if
    end function calc_eqd_grid_3D
    integer function calc_eqd_grid_1D(var,min_grid,max_grid,excl_last) &
        &result(ierr)                                                           ! 1D version
        character(*), parameter :: rout_name = 'calc_eqd_grid_1D'
        
        ! input and output
        real(dp), intent(inout) :: var(:)                                       ! output
        real(dp), intent(in) :: min_grid, max_grid                              ! min. and max. of angles [pi]
        logical, intent(in), optional :: excl_last                              ! .true. if last point excluded
        
        ! local variables
        real(dp), allocatable :: var_3D(:,:,:)                                  ! 3D version of var
        
        ! initialize ierr
        ierr = 0
        
        ! set up var_3D
        allocate(var_3D(size(var),1,1))
        
        ! call 3D version
        ierr = calc_eqd_grid_3D(var_3D,min_grid,max_grid,1,excl_last)
        CHCKERR('')
        
        ! update var
        var = var_3D(:,1,1)
    end function calc_eqd_grid_1D
    
    ! plots the grid in real space
    ! The  equilibrium grid  should contain  the fieldline-oriented  angles with
    ! ang_1 the parallel angle and ang_2 the field line label.
    ! Note:  This  routine  does  not  use  n_theta_plot  and  n_zeta_plot  from
    ! num_vars, but instead  temporarily overwrites them with its  own, since it
    ! has to be 3D also in the axisymmetric case.
    ! [MPI] Collective call
    integer function plot_grid_real(grid) result(ierr)
        use num_vars, only: output_style, alpha_job_nr, grp_rank, grp_n_procs, &
            &no_plots, n_theta_plot, n_zeta_plot
        use grid_vars, only: create_grid, destroy_grid
        
        character(*), parameter :: rout_name = 'plot_grid_real'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! fieldline-oriented grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)             ! X, Y and Z of surface in Axisymmetric coordinates
        real(dp), allocatable :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)             ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
        real(dp), pointer :: X_1_tot(:,:,:), Y_1_tot(:,:,:), Z_1_tot(:,:,:)     ! total X, Y and Z
        real(dp), pointer :: X_2_tot(:,:,:), Y_2_tot(:,:,:), Z_2_tot(:,:,:)     ! total X, Y and Z
        type(grid_type) :: grid_ext                                             ! angularly extended grid
        type(grid_type) :: grid_plot                                            ! grid for plotting
        integer :: id, jd                                                       ! counters
        integer :: n_theta_plot_old                                             ! backup of n_theta_plot
        integer :: n_zeta_plot_old                                              ! backup of n_zeta_plot
        character(len=*), parameter :: anim_name = &
            &'Magnetic field in flux surfaces'                                  ! name of animation
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        call writo('Plotting magnetic field and flux surfaces')
        
        call lvl_ud(1)
        
        ! 0. set up variables
        
        ! save n_theta_plot and n_zeta_plot and change them
        n_theta_plot_old = n_theta_plot
        n_zeta_plot_old = n_zeta_plot
        n_theta_plot = 40
        n_zeta_plot = 160
        
        ! extend grid
        ierr = extend_grid(grid,grid_ext)
        CHCKERR('')
        
        ! restore n_theta_plot and n_zeta_plot
        n_theta_plot = n_theta_plot_old
        n_zeta_plot = n_zeta_plot_old
        
        ! trim extended grid into plot grid
        ierr = trim_grid(grid_ext,grid_plot)
        CHCKERR('')
        
        ! 1. plot flux surfaces
        call writo('writing flux surfaces')
        
        ! calculate X_1,Y_1 and Z_1
        ierr = calc_XYZ_grid(grid_plot,X_1,Y_1,Z_1)
        CHCKERR('')
        
        ! destroy grid
        call destroy_grid(grid_plot)
        
        ! 2. plot field lines
        call writo('writing field lines')
        
        ! trim grid into plot grid
        ierr = trim_grid(grid,grid_plot)
        CHCKERR('')
        
        ! calculate X_2,Y_2 and Z_2
        ierr = calc_XYZ_grid(grid_plot,X_2,Y_2,Z_2)
        CHCKERR('')
        
        ! destroy grids
        call destroy_grid(grid_plot)
        call destroy_grid(grid_ext)
        
        ! get pointers to full X, Y and Z
        ! The reason for  this is that the  plot is not as simple  as usual, and
        ! also  efficiency  is not  the  biggest  priority. Therefore,  all  the
        ! plotting of the group is handled by a single process, the master.
        call get_full_XYZ(X_1,Y_1,Z_1,X_1_tot,Y_1_tot,Z_1_tot,'flux surfaces')
        call get_full_XYZ(X_2,Y_2,Z_2,X_2_tot,Y_2_tot,Z_2_tot,'field lines')
        
        ! delegate to child routines
        select case(output_style)
            case(1)                                                             ! GNUPlot output
                call plot_grid_real_GP(X_1_tot,X_2_tot,Y_1_tot,Y_2_tot,Z_1_tot,&
                    &Z_2_tot,anim_name)
            case(2)                                                             ! HDF5 output
                ierr = plot_grid_real_HDF5(X_1_tot,X_2_tot,Y_1_tot,Y_2_tot,&
                    &Z_1_tot,Z_2_tot,anim_name)
                CHCKERR('')
            case default
                err_msg = 'No style associated with '//trim(i2str(output_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! deallocate
        deallocate(X_1,Y_1,Z_1)
        deallocate(X_2,Y_2,Z_2)
        nullify(X_1_tot,Y_1_tot,Z_1_tot)
        nullify(X_2_tot,Y_2_tot,Z_2_tot)
        
        call lvl_ud(-1)
        
        call writo('Done plotting magnetic field and flux surfaces')
    contains
        ! get pointer to full plotting variables X, Y and Z
        subroutine get_full_XYZ(X,Y,Z,X_tot,Y_tot,Z_tot,merge_name)
            use MPI_ops, only: get_ser_var
            
            ! input / output
            real(dp), intent(in), target :: X(:,:,:), Y(:,:,:), Z(:,:,:)        ! X, Y and Z of either flux surfaces or magnetic field lines
            real(dp), intent(inout), pointer :: X_tot(:,:,:), Y_tot(:,:,:), &
                &Z_tot(:,:,:)                                                   ! pointer to full X, Y and Z
            character(len=*) :: merge_name                                      ! name of variable to be merged
            
            ! local variables
            real(dp), allocatable :: loc_XYZ(:)                                 ! local copy of X, Y or Z in an angular point
            real(dp), allocatable :: ser_loc_XYZ(:)                             ! serial copy of loc_XYZ
            integer, allocatable :: tot_dim(:)                                  ! total dimensions for plot 
            
            ! merge plots for flux surfaces if more than one process
            if (grp_n_procs.gt.1) then                                          ! merge group plots
                ! user output
                call writo('Merging group plots for '//merge_name)
                
                ! get total dimension
                ierr = get_ser_var([size(X,3)],tot_dim)
                CHCKERR('')
                
                ! allocate total X, Y and Z (should have same sizes)
                if (grp_rank.eq.0) then
                    allocate(X_tot(size(X,1),size(X,2),sum(tot_dim)))
                    allocate(Y_tot(size(X,1),size(X,2),sum(tot_dim)))
                    allocate(Z_tot(size(X,1),size(X,2),sum(tot_dim)))
                end if
                
                allocate(loc_XYZ(sum(tot_dim)))
                do jd = 1,size(X,2)
                    do id = 1,size(X,1)
                        loc_XYZ = X(id,jd,:)
                        ierr = get_ser_var(loc_XYZ,ser_loc_XYZ)
                        CHCKERR('')
                        if (grp_rank.eq.0) X_tot(id,jd,:) = ser_loc_XYZ
                        loc_XYZ = Y(id,jd,:)
                        ierr = get_ser_var(loc_XYZ,ser_loc_XYZ)
                        CHCKERR('')
                        if (grp_rank.eq.0) Y_tot(id,jd,:) = ser_loc_XYZ
                        loc_XYZ = Z(id,jd,:)
                        ierr = get_ser_var(loc_XYZ,ser_loc_XYZ)
                        CHCKERR('')
                        if (grp_rank.eq.0) Z_tot(id,jd,:) = ser_loc_XYZ
                    end do
                end do
            else                                                                ! just point
                X_tot => X
                Y_tot => Y
                Z_tot => Z
            end if
        end subroutine get_full_XYZ
        
        ! plot with GNUPlot
        subroutine plot_grid_real_GP(X_1,X_2,Y_1,Y_2,Z_1,Z_2,anim_name)
            ! input / output
            real(dp), intent(in) :: X_1(:,:,:), Y_1(:,:,:), Z_1(:,:,:)          ! X, Y and Z of surface in Axisymmetric coordinates
            real(dp), intent(in) :: X_2(:,:,:), Y_2(:,:,:), Z_2(:,:,:)          ! X, Y and Z of magnetic field lines in Axisymmetric coordinates
            character(len=*), intent(in) :: anim_name                           ! name of animation
            
            ! local variables
            character(len=max_str_ln) :: file_name(2), plot_title(2)            ! file name and title
            character(len=max_str_ln) :: draw_ops(2)                            ! individual plot command
            integer :: n_r                                                      ! number of normal points
            
            ! initialize ierr
            ierr = 0
           
            ! only group masters
            if (grp_rank.eq.0) then
                ! user output
                call writo('Drawing animation with GNUPlot')
                call lvl_ud(1)
                
                ! set n_r
                n_r = size(X_1,3)-1
                
                ! set names
                plot_title(1) = 'Magnetic Flux Surface for alpha job '//&
                    &trim(i2str(alpha_job_nr))
                file_name(1) = 'Flux_surfaces_'//trim(i2str(alpha_job_nr))
                plot_title(2) = 'Magnetic Field Line for alpha job '//&
                    &trim(i2str(alpha_job_nr))
                file_name(2) = 'B_field_'//trim(i2str(alpha_job_nr))
                
                ! write flux surfaces of this process
                call print_GP_3D(trim(plot_title(1)),trim(file_name(1))//&
                    &'.dat',Z_1(:,:,1:n_r),x=X_1(:,:,1:n_r),y=Y_1(:,:,1:n_r),&
                    &draw=.false.)
                
                ! write magnetic field lines
                call print_GP_3D(trim(plot_title(2)),trim(file_name(2))//&
                    &'.dat',Z_2(:,:,2:n_r+1),x=X_2(:,:,2:n_r+1),&
                    &y=Y_2(:,:,2:n_r+1),draw=.false.)
                
                ! add '.dat' to file names
                do id = 1,2
                    file_name(id) = trim(file_name(id))//'.dat'
                end do
                
                ! draw both files
                draw_ops(1) = 'linecolor rgb ''#d3d3d3'' linewidth 1'
                draw_ops(2) = 'linecolor rgb ''black'' linewidth 3'
                call draw_GP_animated(anim_name,file_name,n_r,.false.,&
                    &delay=50,draw_ops=draw_ops)
                
                call lvl_ud(-1)
            end if
        end subroutine plot_grid_real_GP
        
        ! Plot with HDF5
        integer function plot_grid_real_HDF5(X_1,X_2,Y_1,Y_2,Z_1,Z_2,&
            &anim_name) result(ierr)
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
            
            ! only group masters
            if (grp_rank.eq.0) then
                ! user output
                call writo('Drawing animation with HDF5')
                call lvl_ud(1)
                
                ! set up loc_dim and n_r
                loc_dim(:,1) = [size(X_1,1),size(X_1,2),1]
                loc_dim(:,2) = [size(X_2,1),size(X_2,2),1]
                n_r = size(X_1,3) - 1                                           ! should be same for all other X_i, Y_i and Z_i
                
                ! set up plot titles
                plot_title(1) = 'Magnetic Flux Surface for alpha job '//&
                    &trim(i2str(alpha_job_nr))
                plot_title(2) = 'Magnetic Field Line for alpha job '//&
                    &trim(i2str(alpha_job_nr))
                
                ! open HDF5 file
                ierr = open_HDF5_file(file_info,&
                    &'field_lines_in_flux_surfaces_'//&
                    &trim(i2str(alpha_job_nr)),description=anim_name,&
                    &ind_plot=.true.)
                CHCKERR('')
                
                ! create grid for time collection
                allocate(space_col_grids(n_r))
                
                ! loop over all normal points
                do id = 1,n_r
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
                    ierr = print_HDF5_3D_data_item(XYZ(1),file_info,'X_field_'&
                        &//trim(i2str(id)),X_2(:,:,id+1:id+1),loc_dim(:,2))
                    CHCKERR('')
                    
                    ! print data item for Y
                    ierr = print_HDF5_3D_data_item(XYZ(2),file_info,'Y_field_'&
                        &//trim(i2str(id)),Y_2(:,:,id+1:id+1),loc_dim(:,2))
                    CHCKERR('')
                    
                    ! print data item for Z
                    ierr = print_HDF5_3D_data_item(XYZ(3),file_info,'Z_field_'&
                        &//trim(i2str(id)),Z_2(:,:,id+1:id+1),loc_dim(:,2))
                    CHCKERR('')
                    
                    ! print geometry with X, Y and Z data item
                    call print_HDF5_geom(geom,2,XYZ,.true.)
                    
                    ! create a grid with the topology and the geometry
                    ierr = print_HDF5_grid(grids(2),plot_title(2),1,&
                        &grid_top=top,grid_geom=geom,reset=.true.)
                    CHCKERR('')
                    
                    ! C.  merge flux  surface and magnetic  field in  to spatial
                    ! collection
                    ierr = print_HDF5_grid(space_col_grids(id),&
                        &'spatial collection',3,grid_time=id*1._dp,&
                        &grid_grids=grids,reset=.true.)
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
            end if
        end function plot_grid_real_HDF5
    end function plot_grid_real
    
    ! Trim a grid, removing any overlap between the different regions.
    integer function trim_grid(grid_in,grid_out) result(ierr)
        use grid_vars, only: create_grid
        use num_vars, only: grp_n_procs, grp_rank
        use MPI_ops, only: get_ser_var
        
        character(*), parameter :: rout_name = 'trim_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! input grid
        type(grid_type), intent(inout) :: grid_out                              ! trimmed grid
        
        ! local variables
        integer, allocatable :: tot_i_min(:)                                    ! i_min of Equilibrium grid of all processes
        integer, allocatable :: tot_i_max(:)                                    ! i_max of Equilibrium grid of all processes
        integer :: i_lim_out(2)                                                 ! i_lim of output grid
        integer :: n_out(3)                                                     ! n of output grid
        
        ! initialize ierr
        ierr = 0
        
        ! get min_i's of the grid_in
        ierr = get_ser_var([grid_in%i_min],tot_i_min,scatter=.true.)
        CHCKERR('')
        
        ! get max_i's of the grid_in
        ierr = get_ser_var([grid_in%i_max],tot_i_max,scatter=.true.)
        CHCKERR('')
        
        ! set i_lim of output grid
        i_lim_out(1) = grid_in%i_min-tot_i_min(1)+1                             ! subtract i_min of first process
        if (grp_rank.lt.grp_n_procs-1) then                                     ! not last process
            i_lim_out(2) = tot_i_min(grp_rank+2)-tot_i_min(1)                   ! take i_min of next process and subtract i_min of first
        else                                                                    ! last process
            i_lim_out(2) = grid_in%i_max-tot_i_min(1)+1                         ! end of this last group
        end if
        
        ! set n of output grid
        n_out(1) = grid_in%n(1)
        n_out(2) = grid_in%n(2)
        n_out(3) = tot_i_max(grp_n_procs)-tot_i_min(1)+1
        
        ! create new grid
        ierr = create_grid(grid_out,n_out,i_lim_out)
        CHCKERR('')
        
        ! copy arrays
        grid_out%theta_E = grid_in%theta_E(:,:,1:grid_out%grp_n_r)
        grid_out%zeta_E = grid_in%zeta_E(:,:,1:grid_out%grp_n_r)
        grid_out%theta_F = grid_in%theta_F(:,:,1:grid_out%grp_n_r)
        grid_out%zeta_F = grid_in%zeta_F(:,:,1:grid_out%grp_n_r)
        grid_out%r_E = grid_in%r_E(tot_i_min(1):tot_i_max(grp_n_procs))         ! copy total r
        grid_out%r_F = grid_in%r_F(tot_i_min(1):tot_i_max(grp_n_procs))
        if (grid_in%divided) then                                               ! but if input grid divided, grp_r gets priority
            grid_out%grp_r_E = grid_in%grp_r_E(1:grid_out%grp_n_r)
            grid_out%grp_r_F = grid_in%grp_r_F(1:grid_out%grp_n_r)
        end if
    end function trim_grid
    
    ! Extend a  grid angularly using  equidistant variables of  n_theta_plot and
    ! n_zeta_plot angular and  own grp_n_r points in  E coordinates. Optionally,
    ! the grid can  also be converted to  F coordinates if equilibrium  grid and
    ! the variables are provided.
    integer function extend_grid(grid_in,grid_ext,grid_eq,eq) result(ierr)
        use num_vars, only: n_theta_plot, n_zeta_plot
        use grid_vars, only: create_grid
        
        character(*), parameter :: rout_name = 'extend_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! grid to be extended
        type(grid_type), intent(inout) :: grid_ext                              ! extended grid
        type(grid_type), intent(in), optional :: grid_eq                        ! equilibrium grid
        type(eq_type), intent(in), optional :: eq                               ! equilibirum variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (present(grid_eq).neqv.present(eq)) then
            ierr = 1
            err_msg = 'When converting to Flux coordinates, also equilibrium &
                &grid and variables needed'
            CHCKERR(err_msg)
        end if
        
        ! creating  equilibrium  grid  for  the output  that  covers  the  whole
        ! geometry angularly in E coordinates
        ierr = create_grid(grid_ext,&
            &[n_theta_plot,n_zeta_plot,grid_in%n(3)],&
            &[grid_in%i_min,grid_in%i_max])
        CHCKERR('')
        grid_ext%grp_r_E = grid_in%grp_r_E
        ierr = calc_eqd_grid(grid_ext%theta_E,1*pi,3*pi,1)                      ! starting from pi gives nicer plots
        CHCKERR('')
        ierr = calc_eqd_grid(grid_ext%zeta_E,0*pi,2*pi,2)
        CHCKERR('')
        
        ! convert all E coordinates to F coordinates if requested
        if (present(grid_eq)) then
            ierr = coord_E2F(eq,grid_eq,&
                &grid_ext%grp_r_E,grid_ext%theta_E,grid_ext%zeta_E,&
                &grid_ext%grp_r_F,grid_ext%theta_F,grid_ext%zeta_F)
            CHCKERR('')
        end if
    end function extend_grid
    
    ! calculates grp_r_E (style 1) or grp_r_F (style 2)
    integer function calc_grp_r(grid_eq,eq,grp_r,style) result (ierr)
        use num_vars, only: use_pol_flux_E, use_pol_flux_F, eq_style
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        
        character(*), parameter :: rout_name = 'calc_grp_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        integer, intent(in) :: style                                            ! whether to calculate in E or F
        real(dp), intent(inout), allocatable :: grp_r(:)                        ! grp_r
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! allocate grp_r if necessary
        if (allocated(grp_r)) then
            if (size(grp_r).ne.grid_eq%grp_n_r) then
                ierr = 1
                err_msg = 'grp_r needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(grp_r(grid_eq%grp_n_r))
        end if
        
        ! choose which style is being used:
        select case (style)
            case (1)                                                            ! E coords.
                ! choose which equilibrium style is being used:
                !   1:  VMEC
                !   2:  HELENA
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        if (use_pol_flux_E) then                                ! normalized poloidal flux
                            grp_r = eq%flux_p_FD(:,0)/max_flux_p_F
                        else                                                    ! normalized toroidal flux
                            grp_r = eq%flux_t_FD(:,0)/max_flux_t_F
                        end if
                    case (2)                                                    ! HELENA
                        grp_r = eq%flux_p_FD(:,0)                               ! poloidal flux
                    case default
                        err_msg = 'No equilibrium style associated with '//&
                            &trim(i2str(eq_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            case (2)                                                            ! F coords.
                if (use_pol_flux_F) then                                        ! poloidal flux / 2pi
                    grp_r = eq%flux_p_FD(:,0)/(2*pi)
                else                                                            ! toroidal flux / 2pi
                    grp_r = eq%flux_t_FD(:,0)/(2*pi)
                end if
            case default
                err_msg = 'No style associated with '//trim(i2str(style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function calc_grp_r
end module grid_ops
