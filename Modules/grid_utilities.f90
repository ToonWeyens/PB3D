!------------------------------------------------------------------------------!
!   Numerical utilities related to the grids and different coordinate systems  !
!------------------------------------------------------------------------------!
module grid_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type, disc_type

    implicit none
    private
    public coord_F2E, coord_E2F, calc_XYZ_grid, calc_eqd_grid, extend_grid_E, &
        &calc_int_vol, copy_grid, trim_grid, untrim_grid, setup_deriv_data, &
        &setup_interp_data, apply_disc, calc_n_par_X_rich, nufft
#if ldebug
    public debug_calc_int_vol
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_int_vol = .false.                                     ! plot debug information for calc_int_vol
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
    interface setup_deriv_data
        module procedure setup_deriv_data_eqd, setup_deriv_data_reg
    end interface
    interface apply_disc
        module procedure &
            &apply_disc_4D_real, apply_disc_4D_complex, &
            &apply_disc_3D_real, apply_disc_3D_complex, &
            &apply_disc_2D_real, apply_disc_2D_complex, &
            &apply_disc_1D_real, apply_disc_1D_complex
    end interface
    
contains
    ! Converts  Flux  coordinates  (r,theta,zeta)_F to  Equilibrium  coordinates
    ! (r,theta,zeta)_E. Optionally,  two arrays  r_F_array and r_E_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, for E, the poloidal or  toroidal normalized flux is used and for
    ! F, the poloidal or toroidal flux in F coordinates, divided by 2pi.
    integer function coord_F2E_rtz(grid_eq,r_F,theta_F,zeta_F,r_E,&
        &theta_E,zeta_E,r_F_array,r_E_array) result(ierr)                       ! version with r, theta and zeta
        use num_vars, only: tol_zero, eq_style
        use VMEC, only: fourier2real, L_V_c, L_V_s, is_asym_V
        
        character(*), parameter :: rout_name = 'coord_F2E_rtz'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux coords.
        real(dp), intent(inout) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)        ! Equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: dims(3)                                                      ! dimensions of the grid
        real(dp), allocatable :: L_V_c_loc(:,:)                                 ! local version of L_V_c
        real(dp), allocatable :: L_V_s_loc(:,:)                                 ! local version of L_V_s
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        dims(1:2) = grid_eq%n(1:2)
        dims(3) = grid_eq%loc_n_r
        
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
        if (dims(1).ne.size(zeta_E,1) .or. dims(2).ne.size(zeta_E,2) .or. &
            &dims(3).ne.size(zeta_E,3) .or. dims(3).ne.size(r_E)) then
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
        ierr = coord_F2E_r(grid_eq,r_F,r_E,r_F_array,r_E_array)
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
        end select
    contains
        integer function coord_F2E_VMEC() result(ierr)
            use num_vars, only: norm_disc_prec_eq
            use num_ops, only: calc_zero_HH
            use VMEC, only: mnmax_V
            
            character(*), parameter :: rout_name = 'coord_F2E_VMEC'
            
            ! local variables
            real(dp), pointer :: loc_r_F(:) => null()                           ! loc_r in F coords.
            type(disc_type) :: norm_interp_data                                 ! data for normal interpolation
            real(dp), allocatable :: theta_E_guess(:,:)                         ! guess for theta_E for a flux surface
            real(dp), allocatable :: theta_E_guess_3D(:,:,:)                    ! guess for theta_E
            
            ! initialize ierr
            ierr = 0
            
            ! set up loc_r_F
            if (present(r_F_array)) then
                loc_r_F => r_F_array
            else
                loc_r_F => grid_eq%loc_r_F
            end if
            
            ! the toroidal angle is trivial
            zeta_E = - zeta_F                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_V_c and L_V_s
            allocate(L_V_c_loc(mnmax_V,dims(3)))
            allocate(L_V_s_loc(mnmax_V,dims(3)))
            
            ! set up guess
            allocate(theta_E_guess(dims(1),dims(2)))
            allocate(theta_E_guess_3D(dims(1),dims(2),dims(3)))
            theta_E_guess_3D = theta_F
            
            ! set up interpolation data
            ierr = setup_interp_data(loc_r_F,r_F,norm_interp_data,&
                &norm_disc_prec_eq)
            CHCKERR('')
            
            ! interpolate L_V_c and L_V_s at requested normal points r_F
            ierr = apply_disc(L_V_c(:,grid_eq%i_min:grid_eq%i_max,0),&
                &norm_interp_data,L_V_c_loc,2)
            CHCKERR('')
            ierr = apply_disc(L_V_s(:,grid_eq%i_min:grid_eq%i_max,0),&
                &norm_interp_data,L_V_s_loc,2)
            CHCKERR('')
            
            ! the poloidal angle has to be found as the zero of
            !   f = theta_F - theta_E - lambda
            err_msg = calc_zero_HH(dims,theta_E,fun_pol,3,&
                &theta_E_guess_3D,relax_fac=1.0_dp)
            if (err_msg.ne.'') then
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! do a check whether the result is indeed zero
            if (maxval(abs(fun_pol(dims,theta_E,0))).gt.tol_zero*100) then
                err_msg = 'In coord_F2E_VMEC, calculating whether f=0 as a &
                    &check, yields a deviation from 0 by '//trim(r2strt(&
                    &100*maxval(abs(fun_pol(dims,theta_E,0)))))//'%'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! clean up
            nullify(loc_r_F)
            call norm_interp_data%dealloc()
        end function coord_F2E_VMEC
        
        ! function  that  returns  f  =  theta_F  -  theta_V  -  lambda  or  its
        ! derivatives in the poloidal direction.  It uses zeta_E (= zeta_V) from
        ! the parent function.
        function fun_pol(dims,theta_E_in,dpol)
            character(*), parameter :: rout_name = 'fun_pol'
            
            ! input / output
            integer, intent(in) :: dims(3)
            real(dp), intent(in) :: theta_E_in(dims(1),dims(2),dims(3))
            integer, intent(in) :: dpol
            real(dp) :: fun_pol(dims(1),dims(2),dims(3))
            
            ! local variables
            real(dp) :: lam(dims(1),dims(2),dims(3))
            
            ! initialize fun_pol
            fun_pol = 0.0_dp
            
            ! calculate lambda
            ierr = fourier2real(L_V_c_loc,L_V_s_loc,&
                &theta_E_in,zeta_E,lam,sym=[is_asym_V,.true.],deriv=[dpol,0])
            CHCKERR('')
            
            ! calculate the output function
            if (dpol.eq.0) then
                fun_pol = theta_F - theta_E_in - lam
            else if (dpol.eq.1) then
                fun_pol = -1 - lam
            else if (dpol.gt.1) then
                fun_pol = - lam
            end if
        end function fun_pol
    end function coord_F2E_rtz
    integer function coord_F2E_r(grid_eq,r_F,r_E,r_F_array,r_E_array) &
        &result(ierr)                                                           ! version with only r
        use num_vars, only: norm_disc_prec_eq
        
        character(*), parameter :: rout_name = 'coord_F2E_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        real(dp), intent(in) :: r_F(:)                                          ! Flux coords.
        real(dp), intent(inout) :: r_E(:)                                       ! Equilibrium coords.
        real(dp), intent(in), optional, target :: r_F_array(:), r_E_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r                                                          ! dimension of the grid
        real(dp), pointer :: loc_r_E(:) => null()                               ! loc_r in E coords.
        real(dp), pointer :: loc_r_F(:) => null()                               ! loc_r in F coords.
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
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
        
        ! set up loc_r_E and loc_r_F
        if (present(r_F_array)) then
            loc_r_E => r_E_array
            loc_r_F => r_F_array
        else
            loc_r_E => grid_eq%loc_r_E
            loc_r_F => grid_eq%loc_r_F
        end if
        
        ! set up interpolation data
        ierr = setup_interp_data(loc_r_F,r_F,norm_interp_data,&
            &norm_disc_prec_eq)
        CHCKERR('')
        
        ! convert normal position
        ierr = apply_disc(loc_r_E,norm_interp_data,r_E)
        CHCKERR('')
        
        ! clean up
        nullify(loc_r_E,loc_r_F)
        call norm_interp_data%dealloc()
    end function coord_F2E_r
    
    ! Converts  Equilibrium  coordinates  (r,theta,zeta)_E to  Flux  coordinates
    ! (r,theta,zeta)_F. Optionally,  two arrays  r_E_array and r_F_array  can be
    ! provided,  that define  the mapping  between the  both coordinate  system.
    ! Standard, for E, the poloidal or  toroidal normalized flux is used and for
    ! F, the poloidal or toroidal flux in F coordinates, divided by 2pi.
    integer function coord_E2F_rtz(grid_eq,r_E,theta_E,zeta_E,r_F,&
        &theta_F,zeta_F,r_E_array,r_F_array) result(ierr)                       ! version with r, theta and zeta
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'coord_E2F_rtz'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        real(dp), intent(in) :: r_E(:), theta_E(:,:,:), zeta_E(:,:,:)           ! Equilibrium coords.
        real(dp), intent(inout) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)        ! Flux coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables (also used in child functions)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: dims(3)                                                      ! dimensions of the grid
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes
        dims = shape(theta_E)
        
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
        if (dims(1).ne.size(zeta_E,1) .or. dims(2).ne.size(zeta_E,2) .or. &
            &dims(3).ne.size(zeta_E,3) .or. dims(3).ne.size(r_E)) then
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
        ierr = coord_E2F_r(grid_eq,r_E,r_F,r_E_array,r_F_array)
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
        end select
    contains
        integer function coord_E2F_VMEC() result(ierr)
            use num_vars, only: norm_disc_prec_eq
            use VMEC, only: fourier2real, mnmax_V, L_V_c, L_V_s, is_asym_V
            
            character(*), parameter :: rout_name = 'coord_E2F_VMEC'
            
            ! local variables
            real(dp), allocatable :: L_V_c_loc(:,:), L_V_s_loc(:,:)             ! local version of L_V_c and L_V_s
            real(dp), allocatable :: lam(:,:,:)                                 ! lambda
            real(dp), pointer :: loc_r_E(:) => null()                           ! loc_r in E coords.
            type(disc_type) :: norm_interp_data                                 ! data for normal interpolation
            
            ! initialize ierr
            ierr = 0
            
            ! set up loc_r_E
            if (present(r_E_array)) then
                loc_r_E => r_E_array
            else
                loc_r_E => grid_eq%loc_r_E
            end if
            
            ! the toroidal angle is trivial
            zeta_F = - zeta_E                                                   ! conversion VMEC LH -> RH coord. system
            
            ! allocate local copies of L_V_c and L_V_s and lambda
            allocate(L_V_c_loc(mnmax_V,dims(3)))
            allocate(L_V_s_loc(mnmax_V,dims(3)))
            allocate(lam(dims(1),dims(2),dims(3)))
            
            ! set up interpolation data
            ierr = setup_interp_data(loc_r_E,r_E,norm_interp_data,&
                &norm_disc_prec_eq)
            CHCKERR('')
            
            ! interpolate L_V_c and L_V_s at requested normal points r_E
            ierr = apply_disc(L_V_c(:,grid_eq%i_min:grid_eq%i_max,0),&
                &norm_interp_data,L_V_c_loc,2)
            CHCKERR('')
            ierr = apply_disc(L_V_s(:,grid_eq%i_min:grid_eq%i_max,0),&
                &norm_interp_data,L_V_s_loc,2)
            CHCKERR('')
            
            ! calculate lambda
            ierr = fourier2real(L_V_c_loc,L_V_s_loc,theta_E,zeta_E,lam,&
                &sym=[is_asym_V,.true.])
            CHCKERR('')
            
            ! the poloidal angle has to be found as
            !   theta_F = theta_E + lambda
            theta_F = theta_E + lam
            
            ! clean up
            nullify(loc_r_E)
            call norm_interp_data%dealloc()
        end function coord_E2F_VMEC
    end function coord_E2F_rtz
    integer function coord_E2F_r(grid_eq,r_E,r_F,r_E_array,r_F_array) &
        &result(ierr)                                                           ! version with only r
        use num_vars, only: norm_disc_prec_eq
        
        character(*), parameter :: rout_name = 'coord_E2F_r'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid (for normal local limits)
        real(dp), intent(in) :: r_E(:)                                          ! Equilibrium coords.
        real(dp), intent(inout) :: r_F(:)                                       ! Flux coords.
        real(dp), intent(in), optional, target :: r_E_array(:), r_F_array(:)    ! optional arrays that define mapping between two coord. systems
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_r                                                          ! dimension of the grid
        real(dp), pointer :: loc_r_E(:) => null()                               ! loc_r in E coords.
        real(dp), pointer :: loc_r_F(:) => null()                               ! loc_r in F coords.
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
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
        
        ! set up loc_r_E and loc_r_F
        if (present(r_F_array)) then
            loc_r_E => r_E_array
            loc_r_F => r_F_array
        else
            loc_r_E => grid_eq%loc_r_E
            loc_r_F => grid_eq%loc_r_F
        end if
        
        ! set up interpolation data
        ierr = setup_interp_data(loc_r_E,r_E,norm_interp_data,&
            &norm_disc_prec_eq)
        CHCKERR('')
        
        ! convert normal position
        ierr = apply_disc(loc_r_F,norm_interp_data,r_F)
        CHCKERR('')
        
        ! clean up
        nullify(loc_r_E,loc_r_F)
        call norm_interp_data%dealloc()
    end function coord_E2F_r
    
    ! Calculates  X,Y  and  Z  on   a  grid  grid_XYZ,  determined  through  its
    ! E(quilibrium) coordinates.
    ! Furthermore, a grid  grid_eq must be provided, which is  the grid in which
    ! the variables concerning R and Z  are tabulated, i.e. the full equilibrium
    ! grid  in E(quilibrium)  coordinates. Of  this grid,  however, only  r_E is
    ! used,  and the  rest ignored.  It can  therefore be  provided without  the
    ! angular part, i.e. by reconstructing it with a subset.
    ! If VMEC is the equilibrium  model, this routine also optionally calculates
    ! lambda on the grid, as this is  also needed some times. If HELENA is used,
    ! this variable is not used.
    ! Note: For  VMEC, the trigonometric factors of grid_XYZ  must be calculated
    ! beforehand.
    integer function calc_XYZ_grid(grid_eq,grid_XYZ,X,Y,Z,L,R) result(ierr)
        use num_vars, only: eq_style, use_normalization
        use num_utilities, only: round_with_tol
        use eq_vars, only: R_0
        
        character(*), parameter :: rout_name = 'calc_XYZ_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_XYZ                                 ! grid for which to calculate X, Y, Z and optionally L
        real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)                 ! X, Y and Z of grid
        real(dp), intent(inout), optional :: L(:,:,:)                           ! lambda of grid
        real(dp), intent(inout), optional :: R(:,:,:)                           ! R of grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (size(X,1).ne.grid_XYZ%n(1) .or. size(X,2).ne.grid_XYZ%n(2) .or. &
            &size(X,3).ne.grid_XYZ%loc_n_r) then
            ierr = 1
            err_msg =  'X needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (size(Y,1).ne.grid_XYZ%n(1) .or. size(Y,2).ne.grid_XYZ%n(2) .or. &
            &size(Y,3).ne.grid_XYZ%loc_n_r) then
            ierr = 1
            err_msg =  'Y needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (size(Z,1).ne.grid_XYZ%n(1) .or. size(Z,2).ne.grid_XYZ%n(2) .or. &
            &size(Z,3).ne.grid_XYZ%loc_n_r) then
            ierr = 1
            err_msg =  'Z needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(L)) then
            if (size(L,1).ne.grid_XYZ%n(1) .or. size(L,2).ne.grid_XYZ%n(2) &
                &.or. size(L,3).ne.grid_XYZ%loc_n_r) then
                ierr = 1
                err_msg =  'L needs to have the correct dimensions'
                CHCKERR(err_msg)
            end if
        end if
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = calc_XYZ_grid_VMEC(grid_eq,grid_XYZ,X,Y,Z,L=L,R=R)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = calc_XYZ_grid_HEL(grid_eq,grid_XYZ,X,Y,Z,R=R)
                CHCKERR('')
        end select
        
        ! if normalized, return to physical value
        if (use_normalization) then
            X = X*R_0
            Y = Y*R_0
            Z = Z*R_0
            if (present(R)) R = R*R_0
        end if
    contains
        ! VMEC version
        integer function calc_XYZ_grid_VMEC(grid_eq,grid_XYZ,X,Y,Z,L,R) &
            &result(ierr)
            use num_vars, only: norm_disc_prec_eq
            use VMEC, only: fourier2real, &
                &R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, mnmax_V, is_asym_V
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_VMEC'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(grid_type), intent(in) :: grid_XYZ                             ! grid for which to calculate X, Y, Z and optionally L
            real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)             ! X, Y and Z of grid
            real(dp), intent(inout), optional :: L(:,:,:)                       ! lambda of grid
            real(dp), intent(inout), optional :: R(:,:,:)                       ! R of grid
            
            ! local variables
            real(dp), allocatable :: R_V_c_int(:,:), R_V_s_int(:,:)             ! interpolated version of R_V_c and R_V_s
            real(dp), allocatable :: Z_V_c_int(:,:), Z_V_s_int(:,:)             ! interpolated version of Z_V_c and Z_V_s
            real(dp), allocatable :: L_V_c_int(:,:), L_V_s_int(:,:)             ! interpolated version of L_V_c and L_V_s
            real(dp), allocatable :: R_loc(:,:,:)                               ! R in Cylindrical coordinates
            type(disc_type) :: norm_interp_data                                 ! data for normal interpolation
            
            ! initialize ierr
            ierr = 0
            
            ! set up interpolated R_V_c_int, ..
            allocate(R_V_c_int(mnmax_V,grid_XYZ%loc_n_r))
            allocate(R_V_s_int(mnmax_V,grid_XYZ%loc_n_r))
            allocate(Z_V_c_int(mnmax_V,grid_XYZ%loc_n_r))
            allocate(Z_V_s_int(mnmax_V,grid_XYZ%loc_n_r))
            if (present(L)) then
                allocate(L_V_c_int(mnmax_V,grid_XYZ%loc_n_r))
                allocate(L_V_s_int(mnmax_V,grid_XYZ%loc_n_r))
            end if
            
            ! setup normal interpolation data
            ierr = setup_interp_data(grid_eq%r_E,grid_XYZ%loc_r_E,&
                &norm_interp_data,norm_disc_prec_eq)
            CHCKERR('')
            
            ! interpolate VMEC tables
            ierr = apply_disc(R_V_c(:,:,0),norm_interp_data,R_V_c_int,2)
            CHCKERR('')
            ierr = apply_disc(R_V_s(:,:,0),norm_interp_data,R_V_s_int,2)
            CHCKERR('')
            ierr = apply_disc(Z_V_c(:,:,0),norm_interp_data,Z_V_c_int,2)
            CHCKERR('')
            ierr = apply_disc(Z_V_s(:,:,0),norm_interp_data,Z_V_s_int,2)
            CHCKERR('')
            if (present(L)) then
                ierr = apply_disc(L_V_c(:,:,0),norm_interp_data,L_V_c_int,2)
                CHCKERR('')
                ierr = apply_disc(L_V_s(:,:,0),norm_interp_data,L_V_s_int,2)
                CHCKERR('')
            end if
            
            ! clean up
            call norm_interp_data%dealloc()
            
            ! allocate local R
            allocate(R_loc(grid_XYZ%n(1),grid_XYZ%n(2),grid_XYZ%loc_n_r))
            
            ! inverse fourier transform with trigonometric factors
            ierr = fourier2real(R_V_c_int,R_V_s_int,grid_XYZ%trigon_factors,&
                &R_loc,sym=[.true.,is_asym_V])
            CHCKERR('')
            ierr = fourier2real(Z_V_c_int,Z_V_s_int,grid_XYZ%trigon_factors,Z,&
                &sym=[is_asym_V,.true.])
            CHCKERR('')
            if (present(L)) then
                ierr = fourier2real(L_V_c_int,L_V_s_int,&
                    &grid_XYZ%trigon_factors,L,sym=[is_asym_V,.true.])
                CHCKERR('')
            end if
            if (present(R)) R = R_loc
            
            ! transform cylindrical to cartesian
            X = R_loc*cos(grid_XYZ%zeta_E)
            Y = R_loc*sin(grid_XYZ%zeta_E)
            
            ! deallocate
            deallocate(R_V_c_int,R_V_s_int,Z_V_c_int,Z_V_s_int)
            if (present(L)) deallocate(L_V_c_int,L_V_s_int)
            deallocate(R_loc)
        end function calc_XYZ_grid_VMEC
        
        ! HELENA version
        integer function calc_XYZ_grid_HEL(grid_eq,grid_XYZ,X,Y,Z,R) &
            &result(ierr)
            use HELENA_vars, only: R_H, Z_H, chi_H, ias, nchi
            use num_vars, only: norm_disc_prec_eq
            
            character(*), parameter :: rout_name = 'calc_XYZ_grid_HEL'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(grid_type), intent(in) :: grid_XYZ                             ! grid for which to calculate X, Y, Z and optionally L
            real(dp), intent(inout) :: X(:,:,:), Y(:,:,:), Z(:,:,:)             ! X, Y and Z of grid
            real(dp), intent(inout), optional :: R(:,:,:)                       ! R of grid
            
            ! local variables
            integer :: id, jd, kd                                               ! counters
            integer :: pmone                                                    ! plus or minus one
            real(dp), allocatable :: R_H_int(:,:), Z_H_int(:,:)                 ! R and Z at interpolated normal value
            real(dp), allocatable :: R_loc(:,:,:)                               ! R in Cylindrical coordinates
            real(dp) :: theta_loc                                               ! local copy of theta_E
            type(disc_type) :: norm_interp_data                                 ! data for normal interpolation
            type(disc_type) :: pol_interp_data                                  ! data for poloidal interpolation
            
            ! initialize ierr
            ierr = 0
            
            ! allocate local R
            allocate(R_loc(grid_XYZ%n(1),grid_XYZ%n(2),grid_XYZ%loc_n_r))
            
            ! set up interpolated R and Z
            allocate(R_H_int(nchi,grid_XYZ%loc_n_r))
            allocate(Z_H_int(nchi,grid_XYZ%loc_n_r))
            
            ! set up interpolation data
            ierr = setup_interp_data(grid_eq%r_E,grid_XYZ%loc_r_E,&
                &norm_interp_data,norm_disc_prec_eq)
            CHCKERR('')
            
            ! interpolate HELENA output  R_H and Z_H for  every requested normal
            ! point
            ierr = apply_disc(R_H,norm_interp_data,R_H_int,2)
            CHCKERR('')
            ierr = apply_disc(Z_H,norm_interp_data,Z_H_int,2)
            CHCKERR('')
            
            ! Note:  R_H and  Z_H  are not  adapted to  the  parallel grid,  but
            ! tabulated in the original HELENA poloidal grid.
            ! loop over normal points
            do kd = 1,grid_XYZ%loc_n_r                                          ! loop over all normal points
                ! loop over toroidal points
                do jd = 1,grid_XYZ%n(2)
                    ! interpolate at the requested poloidal points
                    do id = 1,grid_XYZ%n(1)
                        ! initialize local theta
                        theta_loc = grid_XYZ%theta_E(id,jd,kd)
                        
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
                        
                        ! take into account possible symmetry
                        if (ias.eq.0 .and. theta_loc.gt.pi) then
                            theta_loc = 2*pi-theta_loc
                            pmone = -1
                        else
                            pmone = 1
                        end if
                        
                        ! set up interpolation data
                        ierr = setup_interp_data(chi_H,[theta_loc],&
                            &pol_interp_data,norm_disc_prec_eq)                 ! use same precision as normal discretization
                        CHCKERR('')
                        
                        ! interpolate  the HELENA  variables poloidally
                        ierr = apply_disc(R_H_int(:,kd),pol_interp_data,&
                            &R_loc(id:id,jd,kd))
                        CHCKERR('')
                        ierr = apply_disc(pmone*Z_H_int(:,kd),pol_interp_data,&
                            &Z(id:id,jd,kd))
                        CHCKERR('')
                    end do
                end do
            end do
            if (present(R)) R = R_loc
            
            ! calculate X and Y, transforming cylindrical to cartesian
            ! (the geometrical zeta is the inverse of HELENA zeta)
            X = R_loc*cos(-grid_XYZ%zeta_E)
            Y = R_loc*sin(-grid_XYZ%zeta_E)
            
            ! deallocate
            deallocate(R_loc)
            call norm_interp_data%dealloc()
            call pol_interp_data%dealloc()
        end function calc_XYZ_grid_HEL
    end function calc_XYZ_grid

    ! Calculate grid of equidistant points,  where optionally the last point can
    ! be excluded.
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
    
    ! Extend a  grid angularly using  equidistant variables of  n_theta_plot and
    ! n_zeta_plot angular and own loc_n_r points in E coordinates.
    ! Optionally,  the grid  can  also  be converted  to  F  coordinates if  the
    ! equilibrium grid is provided. Though this is required for many operations,
    ! it  is  not   required  for  the  calculation  of  X,   Y  and  Z  through
    ! calc_XYZ_grid.
    integer function extend_grid_E(grid_in,grid_ext,grid_eq,n_theta_plot,&
        &n_zeta_plot,lim_theta_plot,lim_zeta_plot) result(ierr)
        use num_vars, only: n_theta => n_theta_plot, n_zeta => n_zeta_plot, &
            &min_theta_plot, max_theta_plot, min_zeta_plot, max_zeta_plot
        
        character(*), parameter :: rout_name = 'extend_grid_E'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! grid to be extended
        type(grid_type), intent(inout) :: grid_ext                              ! extended grid
        type(grid_type), intent(in), optional :: grid_eq                        ! equilibrium grid
        integer, intent(in), optional :: n_theta_plot                           ! number of poins in theta direction
        integer, intent(in), optional :: n_zeta_plot                            ! number of poins in zeta direction
        real(dp), intent(in), optional :: lim_theta_plot(2)                     ! limits in theta
        real(dp), intent(in), optional :: lim_zeta_plot(2)                      ! limits in zeta
        
        ! local variables
        integer :: n_theta_plot_loc                                             ! local n_theta_plot
        integer :: n_zeta_plot_loc                                              ! local n_zeta_plot
        real(dp) :: lim_theta_plot_loc(2)                                       ! local limits of theta_plot
        real(dp) :: lim_zeta_plot_loc(2)                                        ! local limits of zeta_plot
        
        ! initialize ierr
        ierr = 0
        
        ! set up local variables
        n_theta_plot_loc = n_theta
        if (present(n_theta_plot)) n_theta_plot_loc = n_theta_plot
        n_zeta_plot_loc = n_zeta
        if (present(n_zeta_plot)) n_zeta_plot_loc = n_zeta_plot
        lim_theta_plot_loc = [min_theta_plot,max_theta_plot]
        if (present(lim_theta_plot)) lim_theta_plot_loc = lim_theta_plot
        lim_zeta_plot_loc = [min_zeta_plot,max_zeta_plot]
        if (present(lim_zeta_plot)) lim_zeta_plot_loc = lim_zeta_plot
        
        ! creating  equilibrium  grid  for  the output  that  covers  the  whole
        ! geometry angularly in E coordinates
        ierr = grid_ext%init([n_theta_plot_loc,n_zeta_plot_loc,grid_in%n(3)],&
            &[grid_in%i_min,grid_in%i_max])
        CHCKERR('')
        grid_ext%r_E = grid_in%r_E
        grid_ext%loc_r_E = grid_in%loc_r_E
        ierr = calc_eqd_grid(grid_ext%theta_E,lim_theta_plot_loc(1)*pi,&
            &lim_theta_plot_loc(2)*pi,1)
        CHCKERR('')
        ierr = calc_eqd_grid(grid_ext%zeta_E,lim_zeta_plot_loc(1)*pi,&
            &lim_zeta_plot_loc(2)*pi,2)
        CHCKERR('')
        
        ! convert all E coordinates to F coordinates if requested
        if (present(grid_eq)) then
            grid_ext%r_F = grid_in%r_F
            ierr = coord_E2F(grid_eq,&
                &grid_ext%loc_r_E,grid_ext%theta_E,grid_ext%zeta_E,&
                &grid_ext%loc_r_F,grid_ext%theta_F,grid_ext%zeta_F)
            CHCKERR('')
        end if
    end function extend_grid_E
    
    ! Copy a grid A to a new grid B, that was not yet initialized. This new grid
    ! can contain  a subsection of the  previous grid in all  dimensions. It can
    ! also be a divided grid, by providing the limits
    ! Note: these normal limits for the divided grid should be given in terms of
    ! the normal dimension of grid B.
    ! Note: if the grids are not purely normal, the procedure can currently only
    ! handle the copying of a grid_A that is not divided.
    integer function copy_grid(grid_A,grid_B,lims_B,i_lim) result(ierr)
        character(*), parameter :: rout_name = 'copy_grid'
        
        ! input / output
        class(grid_type), intent(in) :: grid_A                                  ! grid to be initialized
        class(grid_type), intent(inout) :: grid_B                               ! grid to be initialized
        integer, intent(in), optional :: lims_B(3,2)                            ! ranges for subset of grid
        integer, intent(in), optional :: i_lim(2)                               ! min. and max. local normal index
        
        ! local variables
        integer :: n_B(3)                                                       ! n of grid B
        integer :: lims_B_loc(3,2)                                              ! local lims_B
        integer :: i_lim_loc(2)                                                 ! local i_lim
        
        ! initialize ierr
        ierr = 0
        
        ! set up dimensions and of B
        lims_B_loc(1,:) = [1,grid_A%n(1)]
        lims_B_loc(2,:) = [1,grid_A%n(2)]
        lims_B_loc(3,:) = [1,grid_A%n(3)]
        if (present(lims_B)) lims_B_loc = lims_B
        n_B = lims_B_loc(:,2)-lims_B_loc(:,1)+1
        where (lims_B_loc(:,1).eq.lims_B_loc(:,2) &
            &.and. lims_B_loc(:,1).eq.[0,0,0]) n_B = 0
        i_lim_loc = lims_B_loc(3,:)
        if (present(i_lim)) i_lim_loc = lims_B_loc(3,1) - 1 + i_lim
        ! tests
        if (n_B(1).gt.grid_A%n(1) .or. n_B(2).gt.grid_A%n(2) .or. &
            &n_B(3).gt.grid_A%n(3)) then
            write(*,*) 'n_B' ,n_B
            write(*,*) 'grid_A', grid_A%n
            ierr = 1
            CHCKERR('lims_B is too large')
        end if
        if (i_lim_loc(1).gt.grid_A%n(3) .or. i_lim_loc(2).gt.grid_A%n(3)) then
            ierr = 1
            CHCKERR('normal limits too wide')
        end if
        
        ! allocate the new grid
        ierr = grid_B%init(n_B,i_lim)
        CHCKERR('')
        
        ! copy arrays, possibly subset
        grid_B%r_F = grid_A%r_F(lims_B_loc(3,1):lims_B_loc(3,2))
        grid_B%r_E = grid_A%r_E(lims_B_loc(3,1):lims_B_loc(3,2))
        grid_B%loc_r_F = grid_A%r_F(i_lim_loc(1):i_lim_loc(2))
        grid_B%loc_r_E = grid_A%r_E(i_lim_loc(1):i_lim_loc(2))
        if (n_B(1).ne.0 .and. n_B(2).ne.0) then
            if (grid_A%divided) then
                ierr = 1
                CHCKERR('grid_A cannot be divided')
            end if
            grid_B%theta_F = grid_A%theta_F(lims_B_loc(1,1):lims_B_loc(1,2),&
                &lims_B_loc(2,1):lims_B_loc(2,2),i_lim_loc(1):i_lim_loc(2))
            grid_B%theta_E = grid_A%theta_E(lims_B_loc(1,1):lims_B_loc(1,2),&
                &lims_B_loc(2,1):lims_B_loc(2,2),i_lim_loc(1):i_lim_loc(2))
            grid_B%zeta_F = grid_A%zeta_F(lims_B_loc(1,1):lims_B_loc(1,2),&
                &lims_B_loc(2,1):lims_B_loc(2,2),i_lim_loc(1):i_lim_loc(2))
            grid_B%zeta_E = grid_A%zeta_E(lims_B_loc(1,1):lims_B_loc(1,2),&
                &lims_B_loc(2,1):lims_B_loc(2,2),i_lim_loc(1):i_lim_loc(2))
        end if
    end function copy_grid
    
    ! Calculates volume integral on a 3D grid.
    ! Two angular and  one normal variable has  to be provided on a the grid. If
    ! the  i'th dimension  of the  grid  is equal  to  one, the  function to  be
    ! integrated is assumed not to vary in this dimension.
    ! Furthermore, if i is 1 or  2, the corresponding i'th (angular) variable is
    ! the only  variable that is  assumed to vary  in that dimension.  The other
    ! angular variable as well as the normal variable are assumed to be constant
    ! like the  function itself. However,  if i is 3,  an error is  displayed as
    ! this does not represent a physical situation.
    ! A common  case through which to  understand this is the  axisymmetric case
    ! where the first  angular variable theta varies in the  dimensions 1 and 3,
    ! the second angular  variable zeta varies only in the  dimension 2, and the
    ! normal variable only varies in the dimension 3.
    ! Alternatively,  there is  the case  of a  grid-aligned set  of coordinates
    ! theta and  alpha, where the  first dimension corresponds to  the direction
    ! along the magnetic field line, the  second to the geodesical direction and
    ! the third to the normal direction. If the calculations for different field
    ! lines are  decoupled, the variation in  the second dimension is  not taken
    ! into account and no integration happens along it.
    ! Internally, the angular  variables and the normal variable  are related to
    ! the coordinates (x,y,z) that correspond to the three dimensions. They thus
    ! form a computational orthogonal grid to which the original coordinates are
    ! related through the transformation of Jacobians:
    !   Jxyz = J r_F_z (ang_1_x ang_2_y-ang_1_y ang_2_x) ,
    ! so that the integral becomes
    !   sum_xyz f(x,y,z) J(x,y,z) r_F_z (ang_1_x ang_2_y-ang_1_y ang_2_x) dxdydz
    ! where dx, dy and dz are all trivially equal to 1.
    ! The integrand has to be evaluated at the intermediate positions inside the
    ! cells. This is  done by taking the  average of the 2^3=8 points  for fJ as
    ! well as the transformation of the Jacobian.
    ! Note: if the coordinates are independent, this method is equivalent to the
    ! repeated numerical integration using the trapezoidal method, NOT Simpson's
    ! 3/8 rule!
    ! Note: by  setting debug_calc_int_vol, this  method can be compared  to the
    ! trapezoidal and simple  method for independent coordinates,  again NOT for
    ! Simpson's 3/8 rule!
    ! The  Simpson's  3/8  rule  could  be developed  but  it  is not  of  great
    ! importance.
    integer function calc_int_vol(ang_1,ang_2,norm,J,f,f_int) result(ierr)
#if ldebug
        use num_vars, only: rank, n_procs
#endif
        
        character(*), parameter :: rout_name = 'calc_int_vol'
        
        ! input / output
        real(dp), intent(in) :: ang_1(:,:,:), ang_2(:,:,:), norm(:)             ! coordinate variables
        real(dp), intent(in) :: J(:,:,:)                                        ! Jacobian
        complex(dp), intent(in) :: f(:,:,:,:)                                   ! input f(n_par,n_geo,n_r,size_X^2)
        complex(dp), intent(inout) :: f_int(:)                                  ! output integrated f
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: nn_mod                                                       ! number of indices for V and V_int
        integer :: k_min                                                        ! minimum k
        integer :: dims(3)                                                      ! real dimensions
        real(dp), allocatable :: transf_J(:,:,:,:)                              ! comps. of transf. between J and Jxyz: theta_x, zeta_y, theta_y, zeta_x and r_F_z
        real(dp), allocatable :: transf_J_tot(:,:,:)                            ! transf. between J and Jxyz: (theta_x*zeta_y-theta_y*zeta_x)*r_F_z
        complex(dp), allocatable :: Jf(:,:,:,:)                                 ! J*f
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: i_a(3), i_b(3)                                               ! limits of building blocks of intermediary variables
        logical :: dim_1(2)                                                     ! whether the dimension is equal to one
#if ldebug
        character(len=max_str_ln), allocatable :: var_names(:,:)                ! names of variables
        complex(dp), allocatable :: f_int_ALT(:)                                ! alternative calculation for output
        complex(dp), allocatable :: f_int_ALT_1D(:,:)                           ! intermediary step in f_int_ALT
        complex(dp), allocatable :: f_int_ALT_2D(:,:,:)                         ! intermediary step in f_int_ALT
        complex(dp), allocatable :: f_int_ALT_ALT(:)                            ! another alternative calculation for output
        integer :: loc_norm(2)                                                  ! local normal index
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set nn_mod
        nn_mod = size(f,4)
        
        ! tests
        if (size(f_int).ne.nn_mod) then
            ierr = 1
            err_msg = 'f and f_int need to have the same storage convention'
            CHCKERR(err_msg)
        end if
        if (size(ang_1).ne.size(ang_2) .or. size(ang_1,3).ne.size(norm)) then
            ierr = 1
            err_msg = 'The angular variables have to be compatible'
            CHCKERR(err_msg)
        end if
        if (size(norm).eq.1) then
            ierr = 1
            err_msg = 'The normal grid size has to be greater than one'
            CHCKERR(err_msg)
        end if
        
        ! set up dims
        dims = shape(ang_1)
        
        ! set up dim_1
        dim_1 = .false.
        if (size(ang_1,1).eq.1) dim_1(1) = .true.
        if (size(ang_2,2).eq.1) dim_1(2) = .true.
        
        ! set up Jf and transf_J
        allocate(Jf(max(dims(1)-1,1),max(dims(2)-1,1),dims(3)-1,nn_mod))
        allocate(transf_J(max(dims(1)-1,1),max(dims(2)-1,1),dims(3)-1,5))
        allocate(transf_J_tot(max(dims(1)-1,1),max(dims(2)-1,1),dims(3)-1))
        
        ! set up k_min
        k_min = 1
        
        ! get intermediate Jf and components of transf_J:
        !   1: ang_1_x
        !   2: ang_2_y
        !   3: ang_1_y
        !   4: ang_2_x
        !   5: r_F_z
        ! Note: the are always 8 terms in the summation, hence the factor 0.125
        Jf = 0._dp
        transf_J = 0._dp
        do kd = 1,2
            if (kd.eq.1) then
                i_a(3) = 1
                i_b(3) = dims(3)-1
            else
                i_a(3) = 2
                i_b(3) = dims(3)
            end if
            do jd = 1,2
                if (dim_1(2)) then
                        i_a(2) = 1
                        i_b(2) = dims(2)
                else
                    if (jd.eq.1) then
                        i_a(2) = 1
                        i_b(2) = dims(2)-1
                    else
                        i_a(2) = 2
                        i_b(2) = dims(2)
                    end if
                end if
                do id = 1,2
                    if (dim_1(1)) then
                            i_a(1) = 1
                            i_b(1) = dims(1)
                    else
                        if (id.eq.1) then
                            i_a(1) = 1
                            i_b(1) = dims(1)-1
                        else
                            i_a(1) = 2
                            i_b(1) = dims(1)
                        end if
                    end if
                    do ld = 1,nn_mod
                        ! Jf
                        Jf(:,:,:,ld) = Jf(:,:,:,ld) + 0.125_dp * &
                            &J(i_a(1):i_b(1),i_a(2):i_b(2),i_a(3):i_b(3)) * &
                            &f(i_a(1):i_b(1),i_a(2):i_b(2),i_a(3):i_b(3),ld)
                    end do
                    if (dim_1(1)) then
                        ! transf_J(1): ang_1_x
                        transf_J(:,:,:,1) = 1._dp
                        ! transf_J(4): ang_2_x
                        transf_J(:,:,:,4) = 0._dp
                    else
                        do ld = 1,dims(1)-1
                            ! transf_J(1): ang_1_x
                            transf_J(ld,:,:,1) = transf_J(ld,:,:,1) + &
                                &0.125_dp * ( &
                                &ang_1(ld+1,i_a(2):i_b(2),i_a(3):i_b(3)) - &
                                &ang_1(ld,i_a(2):i_b(2),i_a(3):i_b(3)) )
                            ! transf_J(4): ang_2_x
                            transf_J(ld,:,:,4) = transf_J(ld,:,:,4) + &
                                &0.125_dp * ( &
                                &ang_2(ld+1,i_a(2):i_b(2),i_a(3):i_b(3)) - &
                                &ang_2(ld,i_a(2):i_b(2),i_a(3):i_b(3)) )
                        end do
                    end if
                    if (dim_1(2)) then
                        ! transf_J(2): ang_2_y
                        transf_J(:,:,:,2) = 1._dp
                        ! transf_J(3): ang_1_y
                        transf_J(:,:,:,3) = 0._dp
                    else
                        do ld = 1,dims(2)-1
                            ! transf_J(2): ang_2_y
                            transf_J(:,ld,:,2) = transf_J(:,ld,:,2) + &
                                &0.125_dp * ( &
                                &ang_2(i_a(1):i_b(1),ld+1,i_a(3):i_b(3)) - &
                                &ang_2(i_a(1):i_b(1),ld,i_a(3):i_b(3)) )
                            ! transf_J(3): ang_1_y
                            transf_J(:,ld,:,3) = transf_J(:,ld,:,3) + &
                                &0.125_dp * ( &
                                &ang_1(i_a(1):i_b(1),ld+1,i_a(3):i_b(3)) - &
                                &ang_1(i_a(1):i_b(1),ld,i_a(3):i_b(3)) )
                        end do
                    end if
                    do ld = 1,dims(3)-1
                        ! transf_J(5): r_F_z
                        transf_J(:,:,ld,5) = transf_J(:,:,ld,5) + 0.125_dp * ( &
                            &norm(ld+1) - norm(ld) )
                    end do
                end do
            end do
        end do
        
        ! set up total transf_J
        transf_J_tot = (transf_J(:,:,:,1)*transf_J(:,:,:,2)-&
            &transf_J(:,:,:,3)*transf_J(:,:,:,4))*transf_J(:,:,:,5)
        
        ! integrate term  over ang_par_F for all equilibrium  grid points, using
        ! dx = dy = dz = 1
        do ld = 1,nn_mod
            f_int(ld) = sum(Jf(:,:,:,ld)*transf_J_tot)
        end do
        
#if ldebug
        if (debug_calc_int_vol) then
            call writo('Testing whether volume integral is calculated &
                &correctly')
            call lvl_ud(1)
            
            allocate(var_names(nn_mod,2))
            var_names(:,1) = 'real part of integrand J f'
            var_names(:,2) = 'imaginary part of integrand J f'
            do ld = 1,nn_mod
                var_names(ld,1) = trim(var_names(ld,1))//' '//trim(i2str(ld))
                var_names(ld,2) = trim(var_names(ld,2))//' '//trim(i2str(ld))
            end do
            
            call plot_HDF5(var_names(:,1),'TEST_RE_Jf_'//trim(i2str(rank)),&
                &rp(Jf))
            call plot_HDF5(var_names(:,2),'TEST_IM_Jf_'//trim(i2str(rank)),&
                &ip(Jf))
            call plot_HDF5('transformation of Jacobians',&
                &'TEST_transf_J_'//trim(i2str(rank)),transf_J_tot)
            
            ! do alternative calculations
            ! assuming that coordinate i varies only in dimensions i!!!
            allocate(f_int_ALT(nn_mod))
            allocate(f_int_ALT_1D(size(f,3),nn_mod))
            allocate(f_int_ALT_2D(size(f,2),size(f,3),nn_mod))
            allocate(f_int_ALT_ALT(nn_mod))
            f_int_ALT = 0._dp
            f_int_ALT_1D = 0._dp
            f_int_ALT_2D = 0._dp
            f_int_ALT_ALT = 0._dp
            
            ! integrate in first coordinate
            do kd = 1,size(f,3)
                do jd = 1,size(f,2)
                    if (size(f,1).gt.1) then
                        do id = 2,size(f,1)
                            f_int_ALT_2D(jd,kd,:) = f_int_ALT_2D(jd,kd,:) + &
                                &0.5*(J(id,jd,kd)*f(id,jd,kd,:)+&
                                &J(id-1,jd,kd)*f(id-1,jd,kd,:))*&
                                &(ang_1(id,jd,kd)-ang_1(id-1,jd,kd))
                        end do
                    else
                        f_int_ALT_2D(jd,kd,:) = J(1,jd,kd)*f(1,jd,kd,:)
                    end if
                end do
            end do
            
            ! integrate in second coordinate
            do kd = 1,size(f,3)
                if (size(f,2).gt.1) then
                    do jd = 2,size(f,2)
                        f_int_ALT_1D(kd,:) = f_int_ALT_1D(kd,:) + &
                            &0.5*(f_int_ALT_2D(jd,kd,:)+&
                            &f_int_ALT_2D(jd-1,kd,:))*&
                            &(ang_2(1,jd,kd)-ang_2(1,jd-1,kd))                  ! assuming that ang_2 is not dependent on dimension 1
                    end do
                else
                    f_int_ALT_1D(kd,:) = f_int_ALT_2D(1,kd,:)
                end if
            end do
            
            ! integrate in third coordinate
            if (size(f,3).gt.1) then
                do kd = 2,size(f,3)
                    f_int_ALT(:) = f_int_ALT(:) + &
                        &0.5*(f_int_ALT_1D(kd,:)+f_int_ALT_1D(kd-1,:))*&
                        &(norm(kd)-norm(kd-1))
                end do
            else
                f_int_ALT = f_int_ALT_1D(1,:)
            end if
            
            ! second alternative, first order integration
            loc_norm = [1,size(f,3)]
            if (rank.lt.n_procs-1) loc_norm(2) = loc_norm(2)-1                  ! no ghost region needed for this method
            do ld = 1,size(f,4)
                f_int_ALT_ALT(ld) = sum(J(:,:,loc_norm(1):loc_norm(2))*&
                    &f(:,:,loc_norm(1):loc_norm(2),ld))
            end do
            if (size(f,1).gt.1) f_int_ALT_ALT = f_int_ALT_ALT*&
                &(ang_1(2,1,1)-ang_1(1,1,1))
            if (size(f,2).gt.1) f_int_ALT_ALT = f_int_ALT_ALT*&
                &(ang_2(1,2,1)-ang_2(1,1,1))
            if (size(f,3).gt.1) f_int_ALT_ALT = f_int_ALT_ALT*&
                &(norm(2)-norm(1))
            
            call writo('Resulting integral: ',persistent=.true.)
            call lvl_ud(1)
            do ld = 1,nn_mod
                call writo('rank '//trim(i2str(rank))//' integral      '//&
                    &trim(i2str(ld))//': '//trim(c2str(f_int(ld))),&
                    &persistent=.true.)
                call writo('rank '//trim(i2str(rank))//'   alternative '//&
                    &trim(i2str(ld))//': '//trim(c2str(f_int_ALT(ld))),&
                    &persistent=.true.)
                call writo('rank '//trim(i2str(rank))//'   other alt.  '//&
                    &trim(i2str(ld))//': '//trim(c2str(f_int_ALT_ALT(ld))),&
                    &persistent=.true.)
            end do
            call lvl_ud(-1)
            
            call writo('(Note that alternative calculations are only valid if &
                &independent coordinates!)')
            
            call lvl_ud(-1)
        end if
#endif
        
        ! deallocate local variables
        deallocate(Jf,transf_J,transf_J_tot)
    end function calc_int_vol
    
    ! Set up  the factors for the  derivative calculation in a  matrix "A" using
    ! finite differences of order "ord" with precision "prec", by which is meant
    ! that the order of the error is at least ~ Delta^(prec+1).
    ! Afterwards, the  matrix A  has to  be multiplied with  the variable  to be
    ! derived to obtain the requested derivative.
    ! For a particular order and precision,  the number of different points that
    ! have to be combined is ord+prec, however, as it is better to use symmetric
    ! expressions, this number is possibly aumented by one. This number is saved
    ! in "n_loc".
    ! The index kd is then defined to go from -(n_loc-1)/2 .. (n_loc-1)/2, which
    ! is  translated to  local coordinates  by adding  (n_loc+1)/2 and  to total
    ! coordinates by adding generally the index in total coordinates i where the
    ! local  problem  is  to be  set  up,  but  capping  it by  (n_loc+1)/2  and
    ! n-(n_loc-1)/2, so by using
    !   k_tot = kd + max((n_loc+1)/2,min(id,n-(n_loc-1)/2)),
    ! so that these never go out of bounds 1 .. n.
    ! The local matrix element "mat_loc" is then set up as follows:
    !   [        1               1       ...       1               1        ]
    !   [     D^i-2_i         D^i-1,i     0     D^i+1_i         D^i+2_i     ]
    !   [ (D^i-2_i)^2/2!  (D^i-1,i)^2/2!  0  (D^i+1_i)^2/2!  (D^i+2_i)^2/2! ]
    !   [ (D^i-2_i)^3/3!  (D^i-1,i)^3/3!  0  (D^i+1_i)^3/3!  (D^i+2_i)^3/3! ]
    !   [       ...             ...      ...      ...             ...       ]
    ! for bulk matrices ((n_loc+1)/2 <= i <= n-(n_loc-1)/2), with
    !   D^j_i = x(j)-x(i). 
    ! For other elements, this is shifted, e.g.:
    !   mat_loc(j,k) = (x(k_tot)-X(i))^(j-1)
    ! where  j = 1..n_loc  and k  = -(n_mod-1)/2..(n_mod-1)/2. k_tot  is defined
    ! above.
    ! The solution A_i of  mat_loc A_i = rhs_loc is then  saved in an individual
    ! row i in the total matrix  A, where rhs_loc = [0,0,..,0,1,0,..,0] with the
    ! unit indicating the order of the derivative.
    ! However, instead of saving the entire matrix A, only the elements from A_i
    ! are saved, omitting the zero's.
    ! For equidistant  grids, the  situation becomes  easier as  D^j_i is  not a
    ! function of id and the rows of local matrix are displaced copies while the
    ! first  rows are  the mirror  image  of the  last  rows so  that the  local
    ! matrices become symmetric. Also, the size of the total variables has to be
    ! passed  in  n,  in contrast  to  the  regular  version,  where it  is  set
    ! automatically.
    integer function setup_deriv_data_eqd(step,n,A,ord,prec) result(ierr)       ! equidistant version
        use num_utilities, only: fac
        
        character(*), parameter :: rout_name = 'setup_deriv_data_eqd'
        
        ! input / output
        real(dp), intent(in) :: step                                            ! step size
        integer, intent(in) :: n                                                ! problem size
        type(disc_type), intent(inout) :: A                                     ! derivation data
        integer, intent(in) :: ord                                              ! order of derivative
        integer, intent(in) :: prec                                             ! precision
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: kd_tot                                                       ! kd in total index
        integer :: n_loc                                                        ! local size of problem to solve
        integer, allocatable :: ipiv(:)                                         ! pivot variable, used by lapack
        real(dp), allocatable :: mat_loc(:,:)                                   ! local matrix
        real(dp), allocatable :: rhs_loc(:)                                     ! local right-hand side
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set n_loc
        n_loc = prec+ord+1
        if (mod(n_loc,2).eq.0) n_loc = n_loc + 1                                ! add one point if even
        
        ! tests
        if (prec.lt.1) then
            ierr = 1
            err_msg = 'precision has to be at least 1'
            CHCKERR(err_msg)
        end if
        if (ord.lt.1) then
            ierr = 1
            err_msg = 'order has to be at least 1'
            CHCKERR(err_msg)
        end if
        if (n.lt.2*n_loc) then
            ierr = 1
            err_msg = 'need at least '//trim(i2str(2*n_loc))//' points in grid'
            CHCKERR(err_msg)
        end if
        
        ! set variables
        allocate(mat_loc(n_loc,n_loc))
        allocate(rhs_loc(n_loc))
        allocate(ipiv(n_loc))
        ipiv = 0
        ierr = A%init(n,n_loc)
        CHCKERR('')
        
        ! iterate over all x values
        do id = 1,n
            ! for bulk of matrix, do calculation only once
            if (id.le.(n_loc+1)/2 .or. id.gt.n-(n_loc-1)/2) then                ! first, last, or first of bulk
                ! calculate elements of local matrix
                mat_loc(1,:) = 1._dp
                do kd = -(n_loc-1)/2,(n_loc-1)/2
                    kd_tot = kd+max((n_loc+1)/2,min(id,n-(n_loc-1)/2))
                    mat_loc(2,kd+(n_loc+1)/2) = (kd_tot-id)*step
                    do jd = 3,n_loc
                        mat_loc(jd,kd+(n_loc+1)/2) = 1._dp/fac(jd-1) * &
                            &mat_loc(2,kd+(n_loc+1)/2)**(jd-1)
                    end do
                end do
                
                ! calculate rhs
                rhs_loc = 0._dp
                rhs_loc(ord+1) = 1._dp                                          ! looking for derivative of this order
                
                ! solve with lapack, making use of lu factorization
                call dgetrf(n_loc,n_loc,mat_loc,n_loc,ipiv,ierr)                ! lu factorization
                err_msg = 'lapack couldn''t find the lu factorization'
                CHCKERR(err_msg)
                call dgetrs('n',n_loc,1,mat_loc,n_loc,ipiv,rhs_loc,n_loc,ierr)  ! solve
                err_msg = 'lapack couldn''t compute the inverse'
                CHCKERR(err_msg)
                
                ! save in A
                A%dat(id,:) = rhs_loc
                A%id_start(id) = max(1,min(id-(n_loc-1)/2,n-n_loc+1))
            else                                                                ! bulk of matrix
                ! copy from first bulk element in matrix A
                A%dat(id,:) = A%dat(id-1,:)
                A%id_start(id) = A%id_start(id-1)+1
            end if
        end do
    end function setup_deriv_data_eqd
    integer function setup_deriv_data_reg(x,A,ord,prec) result(ierr)            ! regular version
        use num_utilities, only: fac
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'setup_deriv_data_reg'
        
        ! input / output
        real(dp), intent(in) :: x(:)                                            ! independent variable
        type(disc_type), intent(inout) :: A                                     ! discretization mat
        integer, intent(in) :: ord                                              ! order of derivative
        integer, intent(in) :: prec                                             ! precision
        
        ! local variables
        integer :: n                                                            ! size of x
        integer :: id, jd, kd                                                   ! counters
        integer :: kd_tot                                                       ! kd in total index
        integer :: n_loc                                                        ! local size of problem to solve
        integer, allocatable :: ipiv(:)                                         ! pivot variable, used by lapack
        real(dp), allocatable :: mat_loc(:,:)                                   ! local matrix
        real(dp), allocatable :: rhs_loc(:)                                     ! local right-hand side
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set n_loc and n
        n_loc = prec+ord+1
        if (mod(n_loc,2).eq.0) n_loc = n_loc + 1                                ! add one point if even
        n = size(x)
        
        ! tests
        if (prec.lt.1) then
            ierr = 1
            err_msg = 'precision has to be at least 1'
            CHCKERR(err_msg)
        end if
        if (ord.lt.1) then
            ierr = 1
            err_msg = 'order has to be at least 1'
            CHCKERR(err_msg)
        end if
        if (n.lt.2*n_loc) then
            ierr = 1
            err_msg = 'need at least '//trim(i2str(2*n_loc))//' points in grid'
            CHCKERR(err_msg)
        end if
        
        ! set variables
        allocate(mat_loc(n_loc,n_loc))
        allocate(rhs_loc(n_loc))
        allocate(ipiv(n_loc))
        ipiv = 0
        ierr = A%init(n,n_loc)
        CHCKERR('')
        
        ! iterate over all x values
        do id = 1,n
            ! calculate elements of local matrix
            mat_loc(1,:) = 1._dp
            do kd = -(n_loc-1)/2,(n_loc-1)/2
                kd_tot = kd+max((n_loc+1)/2,min(id,n-(n_loc-1)/2))
                mat_loc(2,kd+(n_loc+1)/2) = x(kd_tot)-x(id)
                do jd = 3,n_loc
                    mat_loc(jd,kd+(n_loc+1)/2) = 1._dp/fac(jd-1) * &
                        &mat_loc(2,kd+(n_loc+1)/2)**(jd-1)
                end do
            end do
            
            ! calculate local rhs
            rhs_loc = 0._dp
            rhs_loc(ord+1) = 1._dp                                              ! looking for derivative of this order
            
            ! solve with lapack, making use of lu factorization
            call dgetrf(n_loc,n_loc,mat_loc,n_loc,ipiv,ierr)                    ! lu factorization
            err_msg = 'lapack couldn''t find the lu factorization'
            CHCKERR(err_msg)
            call dgetrs('n',n_loc,1,mat_loc,n_loc,ipiv,rhs_loc,n_loc,ierr)      ! solve
            err_msg = 'lapack couldn''t compute the inverse'
            CHCKERR(err_msg)
            
            ! save in total matrix A
            A%dat(id,:) = rhs_loc
            A%id_start(id) = max(1,min(id-(n_loc-1)/2,n-n_loc+1))
        end do
    end function setup_deriv_data_reg
    
    ! Set up the factors for the interpolation calculation in a matrix "A" using
    ! a polynomial order "ord".
    ! This  is  done  using  Barycentric Lagrangian  polynomials,  which  are  a
    ! computationally  more  advantageous  way   of  expressing  the  Lagrangian
    ! interpolating polynomial:
    !   L(x) = sum_j=0^N (w_j f_j /(x-x_j)) / sum_j=0^N (w_j/(x-x_j)) , with 
    !   w_j = 1 / l'(x_j) , where
    !   l(x) = product_i=0^N (x-x_i) , so that
    !   l'(x_j) = product_j.ne.i (x_j-x_i)
    ! [https://people.maths.ox.ac.uk/trefethen/barycentric.pdf, 02/05/2016]
    !  The  optimal range  of grid  points  "x" around  the interpolation  point
    ! "x_interp" is first  calculated, after which the weight  functions w_j are
    ! calculated  for  each of  the  interpolation  points,  multiplied by  1  /
    ! (x-x_j). The result is scaled by the sum of all, and this is tabulated.
    ! Note: To avoid numerical problems, the  factors (x_j-x_i) are divided by a
    ! characteristic  length len,  equal to  the average  length of  the working
    ! interval. This  scales both the  enumerator and the denominator  by len^N,
    ! which  therefore cancel.  Also,  if the  interpolation  point x_interp  is
    ! within a tiny tolerance of the grid points x, the machinery is bypassed to
    ! avoid unstabilities.
    ! Optionally,  this routine  can do  trigonometric interpolation  instead of
    ! polynomial, which results in using sin((x-x_j)/2) instead of (x-x_j).
    ! [https://en.wikipedia.org/wiki/Trigonometric_interpolation, 01/09/2016]
    ! Currently, this only works  for an odd number of points,  as there is also
    ! an odd number of constraints.
    ! Note that this can  be used to resample a periodic  function that has been
    ! sampled on an irregular grid, in  order to calculate the fourier transform
    ! with a FFT, which requires a regular grid. In this case, a full 2pi period
    ! is assumed.
    ! The  fact  that  this is  can  be  done  can  be seen  by  realizing  that
    ! trigonometric  interpolation can  also be  written in  Lagrange form,  and
    ! subsequently in barycentric  Lagrange form. There is  an additional factor
    ! 1/2 in l'(x_j),  due to the derivation of sin((x-x_k)/2)  but this cancels
    ! as it is present in both enumerator and denominator.
    ! Optionally  also  a normalization  length  can  be  provided so  that  the
    ! individual terms in  the interpolation do not blow up.  This can be useful
    ! for example for interpolation with many points.
    integer function setup_interp_data(x,x_interp,A,ord,is_trigon,norm_len) &
        &result(ierr)
        use grid_vars, only: disc_type
        use num_utilities, only: con2dis
        
        character(*), parameter :: rout_name = 'setup_interp_data'
        
        ! input / output
        real(dp), intent(in) :: x(:)                                            ! independent variable of tabulated magnitude
        real(dp), intent(in) :: x_interp(:)                                     ! values of x at which to interpolate
        type(disc_type), intent(inout) :: A                                     ! interpolation mat
        integer, intent(in) :: ord                                              ! order of derivative
        logical, intent(in), optional :: is_trigon                              ! trigonometric interpolation
        real(dp), intent(in), optional :: norm_len                              ! custom normalization length
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: n                                                            ! size of x_interp
        integer :: n_loc                                                        ! nr. of points used for interpolation locally
        real(dp) :: len                                                         ! interval length
        real(dp) :: x_interp_disc                                               ! discrete index of current x_interp
        real(dp), parameter :: tol = 1.0E-8                                     ! tolerance
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: weight(:)                                      ! weights w_j
        real(dp) :: weight_loc                                                  ! local weight factor
        logical :: is_trigon_loc                                                ! local is_trigon
        
        ! initialize ierr
        ierr = 0
        
        ! set local is_trigon
        is_trigon_loc = .false.
        if (present(is_trigon)) is_trigon_loc = is_trigon
        
        ! tests
        if (ord.lt.1) then
            ierr = 1
            err_msg = 'order has to be at least 1'
            CHCKERR(err_msg)
        else if (ord+1.gt.size(x)) then
            ierr = 1
            err_msg = 'order can be at most size(x)-1 = '//&
                &trim(i2str(size(x)-1))
            CHCKERR(err_msg)
        end if
        if (is_trigon_loc .and. mod(ord,2).ne.0) then
            ierr = 1
            err_msg = 'Need an even order, so that there is an odd number of &
                &points'
            CHCKERR(err_msg)
        end if
        
        ! set n_loc and n
        n_loc = ord+1
        n = size(x_interp)
        
        ! set variables
        allocate(weight(n_loc))
        ierr = A%init(n,n_loc)
        CHCKERR('')
        
        ! iterate over all x_interp values
        do id = 1,n
            ! find the base of the interpolation value
            ierr = con2dis(x_interp(id),x_interp_disc,x)
            CHCKERR('')
            
            ! set up start id
            if (n_loc.lt.size(x)) then
                A%id_start(id) = ceiling(x_interp_disc-n_loc*0.5_dp)            ! get minimum of bounding indices
                A%id_start(id) = max(1,min(A%id_start(id),size(x)-n_loc+1))     ! limit to lie within x
            else
                A%id_start(id) = 1                                              ! all values of x are used
            end if
            
            ! check for (near) exact match
            if (mod(x_interp_disc,1._dp).lt.tol) then
                ! directly set the correct index to 1
                A%dat(id,:) = 0._dp
                A%dat(id,nint(x_interp_disc)-A%id_start(id)+1) = 1._dp
            else
                ! calculate w_j/(x-x_j)
                len = (x(A%id_start(id)+n_loc-1) - x(A%id_start(id)))/n_loc
                if (is_trigon_loc) len = sin(len/2)                             ! take sines
                if (present(norm_len)) len = norm_len
                
                ! set up the interp. coeffs. due to each of the points used
                weight = 1._dp
                do jd = 1,n_loc
                    do kd = 1,n_loc
                        ! basis of local weight
                        if (kd.ne.jd) then                                      ! skip k = j
                            ! calculate l'(x_j) = product_k.ne.j (x_j-x_k)
                            weight_loc = &
                                &x(A%id_start(id)-1+jd)-x(A%id_start(id)-1+kd)
                        else
                            ! multiply additionally by (x_interp - x_j)
                            weight_loc = x_interp(id)-x(A%id_start(id)-1+jd)
                        end if
                        
                        ! modification of local weight
                        if (is_trigon_loc) weight_loc = sin(weight_loc/2)       ! take sines
                        weight_loc = weight_loc/len                             ! normalize
                        
                        ! including modified local weight
                        weight(jd) = weight(jd) * weight_loc
                    end do
                end do
                
                ! invert
                weight = 1._dp/weight
                
                ! scale by sum and save in A
                A%dat(id,:) = weight/sum(weight)
            end if
        end do
    end function setup_interp_data
    
    ! Applies  the   discretization  data  calculated  in   setup_deriv_data  or
    ! setup_interp_data  to  calculate  the  derivative or  interpolation  in  a
    ! dimension "disc_dim"
    integer function apply_disc_4D_real(var,disc_data,dvar,disc_dim) &
        &result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_4D_real'
        
        ! input / output
        real(dp), intent(in) :: var(:,:,:,:)                                    ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        real(dp), intent(inout) :: dvar(:,:,:,:)                                ! operated variable
        integer :: disc_dim                                                     ! dimension in which to discretization operation
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n                                                            ! size of variable
        integer :: id, jd                                                       ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (disc_dim.lt.1 .or. disc_dim.gt.4) then
            err_msg = 'invalid dimension of derivative '//trim(i2str(disc_dim))
            CHCKERR(err_msg)
        end if
        n = size(dvar,disc_dim)
        if (n.ne.disc_data%n) then
            err_msg = 'variables are not compatible'
            CHCKERR(err_msg)
        end if
        
        ! initialize output
        dvar = 0._dp
        
        ! multiply the matrices
        select case (disc_dim)
            case (1)
                do id = 1,n
                    do jd = 1,disc_data%n_loc
                        dvar(id,:,:,:) = dvar(id,:,:,:) + &
                            &disc_data%dat(id,jd)*&
                            &var(jd+disc_data%id_start(id)-1,:,:,:)
                    end do
                end do
            case (2)
                do id = 1,n
                    do jd = 1,disc_data%n_loc
                        dvar(:,id,:,:) = dvar(:,id,:,:) + &
                            &disc_data%dat(id,jd)*&
                            &var(:,jd+disc_data%id_start(id)-1,:,:)
                    end do
                end do
            case (3)
                do id = 1,n
                    do jd = 1,disc_data%n_loc
                        dvar(:,:,id,:) = dvar(:,:,id,:) + &
                            &disc_data%dat(id,jd)*&
                            &var(:,:,jd+disc_data%id_start(id)-1,:)
                    end do
                end do
            case (4)
                do id = 1,n
                    do jd = 1,disc_data%n_loc
                        dvar(:,:,:,id) = dvar(:,:,:,id) + &
                            &disc_data%dat(id,jd)*&
                            &var(:,:,:,jd+disc_data%id_start(id)-1)
                    end do
                end do
            case default
                ierr = 1
                err_msg = 'dimension of discretization operation '//&
                    &trim(i2str(disc_dim))//'impossible'
                CHCKERR(err_msg)
        end select
    end function apply_disc_4D_real
    integer function apply_disc_4D_complex(var,disc_data,dvar,disc_dim) &
        &result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_4D_complex'
        
        ! input / output
        complex(dp), intent(in) :: var(:,:,:,:)                                 ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        complex(dp), intent(inout) :: dvar(:,:,:,:)                             ! operated variable
        integer :: disc_dim                                                     ! dimension in which to discretization operation
        
        ! local variables
        real(dp), allocatable :: dvar_loc(:,:,:,:)                              ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),size(dvar,2),size(dvar,3),size(dvar,4)))
        
        ! call the real version for the real and imaginary part separately
        ierr = apply_disc_4D_real(rp(var),disc_data,dvar_loc,disc_dim)
        CHCKERR('')
        dvar = dvar_loc
        ierr = apply_disc_4D_real(ip(var),disc_data,dvar_loc,disc_dim)
        CHCKERR('')
        dvar = dvar + iu*dvar_loc
    end function apply_disc_4D_complex
    integer function apply_disc_3D_real(var,disc_data,dvar,disc_dim) &
        &result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_3D_real'
        
        ! input / output
        real(dp), intent(in) :: var(:,:,:)                                      ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        real(dp), intent(inout) :: dvar(:,:,:)                                  ! operated variable
        integer :: disc_dim                                                     ! dimension in which to discretization operation
        
        ! local variables
        real(dp), allocatable :: dvar_loc(:,:,:,:)                              ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),size(dvar,2),size(dvar,3),1))
        
        ! call the 4D version
        ierr = apply_disc_4D_real(reshape(var,[shape(var),1]),disc_data,&
            &dvar_loc,disc_dim)
        CHCKERR('')
        dvar = dvar_loc(:,:,:,1)
    end function apply_disc_3D_real
    integer function apply_disc_3D_complex(var,disc_data,dvar,disc_dim) &
        &result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_3D_complex'
        
        ! input / output
        complex(dp), intent(in) :: var(:,:,:)                                   ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        complex(dp), intent(inout) :: dvar(:,:,:)                               ! operated variable
        integer :: disc_dim                                                     ! dimension in which to discretization operation
        
        ! local variables
        complex(dp), allocatable :: dvar_loc(:,:,:,:)                           ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),size(dvar,2),size(dvar,3),1))
        
        ! call the 4D version
        ierr = apply_disc_4D_complex(reshape(var,[shape(var),1]),disc_data,&
            &dvar_loc,disc_dim)
        CHCKERR('')
        dvar = dvar_loc(:,:,:,1)
    end function apply_disc_3D_complex
    integer function apply_disc_2D_real(var,disc_data,dvar,disc_dim) &
        &result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_2D_real'
        
        ! input / output
        real(dp), intent(in) :: var(:,:)                                        ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        real(dp), intent(inout) :: dvar(:,:)                                    ! operated variable
        integer :: disc_dim                                                     ! dimension in which to discretization operation
        
        ! local variables
        real(dp), allocatable :: dvar_loc(:,:,:,:)                              ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),size(dvar,2),1,1))
        
        ! call the 4D version
        ierr = apply_disc_4D_real(reshape(var,[shape(var),1,1]),disc_data,&
            &dvar_loc,disc_dim)
        CHCKERR('')
        dvar = dvar_loc(:,:,1,1)
    end function apply_disc_2D_real
    integer function apply_disc_2D_complex(var,disc_data,dvar,disc_dim) &
        &result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_2D_complex'
        
        ! input / output
        complex(dp), intent(in) :: var(:,:)                                     ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        complex(dp), intent(inout) :: dvar(:,:)                                 ! operated variable
        integer :: disc_dim                                                     ! dimension in which to discretization operation
        
        ! local variables
        complex(dp), allocatable :: dvar_loc(:,:,:,:)                           ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),size(dvar,2),1,1))
        
        ! call the 4D version
        ierr = apply_disc_4D_complex(reshape(var,[shape(var),1,1]),disc_data,&
            &dvar_loc,disc_dim)
        CHCKERR('')
        dvar = dvar_loc(:,:,1,1)
    end function apply_disc_2D_complex
    integer function apply_disc_1D_real(var,disc_data,dvar) result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_1D_real'
        
        ! input / output
        real(dp), intent(in) :: var(:)                                          ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        real(dp), intent(inout) :: dvar(:)                                      ! operated variable
        
        ! local variables
        real(dp), allocatable :: dvar_loc(:,:,:,:)                              ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),1,1,1))
        
        ! call the 4D version
        ierr = apply_disc_4D_real(reshape(var,[shape(var),1,1,1]),disc_data,&
            &dvar_loc,1)
        CHCKERR('')
        dvar = dvar_loc(:,1,1,1)
    end function apply_disc_1D_real
    integer function apply_disc_1D_complex(var,disc_data,dvar) result(ierr)
        use grid_vars, only: disc_type
        
        character(*), parameter :: rout_name = 'apply_disc_1D_complex'
        
        ! input / output
        complex(dp), intent(in) :: var(:)                                       ! variable to be operated on
        type(disc_type), intent(in) :: disc_data                                ! disc_data calculated in setup_deriv_data or setup_interp_data
        complex(dp), intent(inout) :: dvar(:)                                   ! operated variable
        
        ! local variables
        complex(dp), allocatable :: dvar_loc(:,:,:,:)                           ! local dvar
        
        ! initialize ierr
        ierr = 0
        
        ! set up local dvar
        allocate(dvar_loc(size(dvar,1),1,1,1))
        
        ! call the 4D version
        ierr = apply_disc_4D_complex(reshape(var,[shape(var),1,1,1]),&
            &disc_data,dvar_loc,1)
        CHCKERR('')
        dvar = dvar_loc(:,1,1,1)
    end function apply_disc_1D_complex
    
    ! Trim a grid, removing any overlap between the different regions.
    ! by default, the routine assumes a  symmetric ghost region and cuts as many
    ! grid points from the end of the  previous process as from the beginning of
    ! the next process, but if the number of overlapping grid points is odd, the
    ! previous process looses one more point.
    ! optionally, the trimmed indices in the normal direction can be provided in
    ! "norm_id", i.e. the indices in the  old, untrimmed grid that correspond to
    ! the start and end indices of the trimmed grid. E.g. if
    !   - proc 0:  3 ... 25
    !   - proc 1: 20 ... 50
    ! then the trimmed grid will be:
    !   - proc 0:  3 ... 22
    !   - proc 1: 23 ... 50
    ! which is shifted down by 2 to 
    !   - proc 0:  1 ... 20
    !   - proc 1: 21 ... 48
    ! in the trimmed grid. The indices of the previous step (3 & 22 and 23 & 50)
    ! are saved in norm_id.
    integer function trim_grid(grid_in,grid_out,norm_id) result(ierr)
        use num_vars, only: n_procs, rank
        use mpi_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'trim_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! input grid
        type(grid_type), intent(inout) :: grid_out                              ! trimmed grid
        integer, intent(inout), optional :: norm_id(2)                          ! normal indices corresponding to trimmed part
        
        ! local variables
        integer, allocatable :: tot_i_min(:)                                    ! i_min of grid of all processes
        integer, allocatable :: tot_i_max(:)                                    ! i_max of grid of all processes
        integer :: i_lim_out(2)                                                 ! i_lim of output grid
        integer :: n_out(3)                                                     ! n of output grid
        
        ! initialize ierr
        ierr = 0
        
        ! detect whether grid divided
        if (grid_in%divided) then
            ! get min_i's of the grid_in
            ierr = get_ser_var([grid_in%i_min],tot_i_min,scatter=.true.)
            CHCKERR('')
            
            ! get max_i's of the grid_in
            ierr = get_ser_var([grid_in%i_max],tot_i_max,scatter=.true.)
            CHCKERR('')
            
            ! set  i_lim of trimmed output  grid (not yet shifted  by first proc
            ! min)
            if (rank.gt.0) then
                i_lim_out(1) = max(tot_i_min(1),tot_i_min(rank+1)+&
                    &floor((tot_i_max(rank)-tot_i_min(rank+1)+1._dp)/2))
            else
                i_lim_out(1) = tot_i_min(1)
            end if
            if (rank.lt.n_procs-1) then
                i_lim_out(2) = min(tot_i_max(n_procs),tot_i_max(rank+1)-&
                    &ceiling((tot_i_max(rank+1)-tot_i_min(rank+2)+1._dp)/2))
            else
                i_lim_out(2) = tot_i_max(n_procs)
            end if
            
            ! get  min_i's and max_i's  of the grid_out,  not shifted by  min of
            ! first process
            ierr = get_ser_var([i_lim_out(1)],tot_i_min,scatter=.true.)
            CHCKERR('')
            ierr = get_ser_var([i_lim_out(2)],tot_i_max,scatter=.true.)
            CHCKERR('')
            
            ! set n of output grid
            n_out(1) = grid_in%n(1)
            n_out(2) = grid_in%n(2)
            n_out(3) = sum(tot_i_max-tot_i_min+1)
            
            ! create new grid
            ierr = grid_out%init(n_out,i_lim_out-tot_i_min(1)+1)                ! limits shifted by min of first process
            CHCKERR('')
            
            ! recycle  i_lim_out  for  shifted  array  indices, set  norm_id  if
            ! requested
            i_lim_out = i_lim_out - grid_in%i_min + 1
            if (present(norm_id)) norm_id = i_lim_out
        else
            ! set n of output grid
            n_out = grid_in%n
            
            ! set unshifted i_lim_out and norm_id if requested
            i_lim_out = [grid_in%i_min,grid_in%i_max]
            if (present(norm_id)) norm_id = i_lim_out
            
            ! shift i_lim_out by min.
            i_lim_out = i_lim_out-i_lim_out(1)+1
            
            ! create new grid
            ierr = grid_out%init(n_out,i_lim_out)                               ! grid not divided
            CHCKERR('')
        end if
        
        ! copy local arrays
        if (grid_in%n(1).ne.0 .and. grid_in%n(2).ne.0) then                     ! only if 3d grid
            grid_out%theta_e = grid_in%theta_e(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%zeta_e = grid_in%zeta_e(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%theta_f = grid_in%theta_f(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%zeta_f = grid_in%zeta_f(:,:,i_lim_out(1):i_lim_out(2))
        end if
        if (grid_in%divided) then                                               ! but if input grid divided, loc_r gets priority
            grid_out%loc_r_E = grid_in%loc_r_E(i_lim_out(1):i_lim_out(2))
            grid_out%loc_r_F = grid_in%loc_r_F(i_lim_out(1):i_lim_out(2))
        end if
        
        ! if divided, set total arrays
        if (grid_in%divided) then
            grid_out%r_E = grid_in%r_E(tot_i_min(1):tot_i_max(n_procs))
            grid_out%r_F = grid_in%r_F(tot_i_min(1):tot_i_max(n_procs))
        else
            grid_out%r_E = grid_in%r_E
            grid_out%r_F = grid_in%r_F
        end if
    end function trim_grid
    
    ! untrims a trimmed  grid by introducing an assymetric ghost  regions at the
    ! right. the width of the ghost region has to be provided.
    ! note: the ghosted grid should be deallocated (with dealloc_grid).
    ! note: the input grid has to be trimmed!
    integer function untrim_grid(grid_in,grid_out,size_ghost) result(ierr)
        use num_vars, only: n_procs, rank
        use mpi_utilities, only: get_ghost_arr, get_ser_var
        
        character(*), parameter :: rout_name = 'untrim_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  ! input grid
        type(grid_type), intent(inout) :: grid_out                              ! ghosted grid
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: i_lim_in(2)                                                  ! limits of input grid
        integer :: i_lim_out(2)                                                 ! limits of ghosted grid
        integer, allocatable :: tot_loc_n_r(:)                                  ! loc_n_r of all processes
        integer :: size_ghost_loc                                               ! local size_ghost
        
        ! initialize ierr
        ierr = 0
        
        ! get array sizes of all processes
        ierr = get_ser_var([grid_in%loc_n_r],tot_loc_n_r,scatter=.true.)
        CHCKERR('')
        
        ! set local size_ghost
        size_ghost_loc = min(size_ghost,minval(tot_loc_n_r)-1)
        
        ! set i_lim_in and i_lim_out
        i_lim_in = [grid_in%i_min,grid_in%i_max]
        i_lim_out = i_lim_in
        if (rank.lt.n_procs-1) i_lim_out(2) = i_lim_out(2)+size_ghost_loc
        
        ! create grid
        ierr = grid_out%init(grid_in%n,i_lim_out)
        CHCKERR('')
        
        ! set ghosted variables
        if (grid_in%n(1).ne.0 .and. grid_in%n(2).ne.0) then                     ! only if 3d grid
            grid_out%theta_e(:,:,1:grid_in%loc_n_r) = grid_in%theta_e
            grid_out%zeta_e(:,:,1:grid_in%loc_n_r) = grid_in%zeta_e
            grid_out%theta_f(:,:,1:grid_in%loc_n_r) = grid_in%theta_f
            grid_out%zeta_f(:,:,1:grid_in%loc_n_r) = grid_in%zeta_f
            ierr = get_ghost_arr(grid_out%theta_e,size_ghost_loc)
            CHCKERR('')
            ierr = get_ghost_arr(grid_out%zeta_e,size_ghost_loc)
            CHCKERR('')
            ierr = get_ghost_arr(grid_out%theta_f,size_ghost_loc)
            CHCKERR('')
            ierr = get_ghost_arr(grid_out%zeta_f,size_ghost_loc)
            CHCKERR('')
        end if
        if (grid_in%divided) then                                               ! but if input grid divided, loc_r gets priority
            grid_out%loc_r_e = grid_in%r_e(i_lim_out(1):i_lim_out(2))
            grid_out%loc_r_f = grid_in%r_f(i_lim_out(1):i_lim_out(2))
        end if
        grid_out%r_e = grid_in%r_e
        grid_out%r_f = grid_in%r_f
    end function untrim_grid
    
    ! Calculates the  local number of  parallel grid points for  this Richardson
    ! level, taking into account that it ould be half the actual number.
    integer function calc_n_par_X_rich(n_par_X_rich,only_half_grid) result(ierr)
        use rich_vars, only: n_par_X
        
        character(*), parameter :: rout_name = 'calc_n_par_X_rich'
        
        ! input / output
        integer, intent(inout) :: n_par_X_rich                                  ! n_par_X for this Richardson level
        logical, intent(in), optional :: only_half_grid                         ! calculate only half grid with even points
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! possibly divide n_par_X by 2
        n_par_X_rich = n_par_X
        if (present(only_half_grid)) then
            if (only_half_grid) then
                if (mod(n_par_X,2).eq.1) then
                    n_par_X_rich = (n_par_X-1)/2
                else
                    ierr = 1
                    err_msg = 'Need odd number of points'
                    CHCKERR(err_msg)
                end if
            end if
        end if
    end function calc_n_par_X_rich
    
    ! calculates the  cosine and sine  mode numbers of  a function defined  on a
    ! non-regular grid.
    ! If a plot name is provided, the modes are plotted.
    ! Note that the fundamental interval is  assumed to be 0..2pi but that there
    ! should be no overlap between the first and last point.
    integer function nufft(x,f,f_F,plot_name) result(ierr)
        character(*), parameter :: rout_name = 'nufft'
        
        ! input / output
        real(dp), intent(in) :: x(:)                                            ! coordinate values
        real(dp), intent(in) :: f(:)                                            ! function values
        real(dp), intent(inout), allocatable :: f_F(:,:)                        ! Fourier modes for cos and sin
        character(len=*), intent(in), optional :: plot_name                     ! name of possible plot
        
        ! local variables
        integer :: m_F                                                          ! nr. of modes
        integer :: n_x                                                          ! nr. of points
        integer :: interp_ord = 4                                               ! order of interpolation
        integer :: id                                                           ! counter
        real(dp), allocatable :: work(:)                                        ! work array
        real(dp), allocatable :: f_loc(:)                                       ! local copy of f
        type(disc_type) :: trigon_interp_data                                   ! data for non-equidistant sampling fourier coefficients
        character(len=max_str_ln) :: plot_title(2)                              ! name of plot
        logical :: print_log = .false.                                          ! print log plot as well
        
        ! tests
        if (size(x).ne.size(f)) then
            ierr = 1
            CHCKERR('x and f must have same size')
        end if
        
        ! set local f and interpolate
        n_x = size(x)
        allocate(f_loc(n_x))
        ierr = setup_interp_data([x,2*pi],[((id-1._dp)/n_x*2*pi,id=1,n_x)],&
            &trigon_interp_data,interp_ord)
        ierr = apply_disc([f,f(1)],trigon_interp_data,f_loc)
        CHCKERR('')
        call trigon_interp_data%dealloc()
        
        !!do id = 1,n_x
            !!f_loc(id) = -4._dp+&
                !!&3*cos((id-1._dp)/n_x*2*pi*1)+&
                !!&0.5*cos((id-1._dp)/n_x*2*pi*100)+&
                !!&1.5*cos((id-1._dp)/n_x*2*pi*101)+&
                !!&4*sin((id-1._dp)/n_x*2*pi*1)+&
                !!&2*sin((id-1._dp)/n_x*2*pi*50)+&
                !!&4.5*sin((id-1._dp)/n_x*2*pi*55)
        !!end do
        
        ! set up variables for fft
        m_F = (n_x-1)/2                                                         ! nr. of modes
        allocate(work(3*n_x+15))                                                ! from fftpack manual
        
        ! calculate fft
        call dffti(n_x,work)
        call dfftf(n_x,f_loc,work)
        deallocate(work)
        
        ! rescale
        f_loc(:) = f_loc(:)*2/n_x
        f_loc(1) = f_loc(1)/2
        
        ! separate cos and sine
        if (allocated(f_F)) deallocate(f_F)
        allocate(f_F(m_F+1,2))
        f_F(1,1) = f_loc(1)
        f_F(1,2) = 0._dp
        f_F(2:m_F+1,1) = f_loc(2:2*m_F:2)
        f_F(2:m_F+1,2) = -f_loc(3:2*m_F+1:2)                                    ! routine returns - sine factors
        
        ! output in plot if requested
        if (present(plot_name)) then
            plot_title = ['cos','sin']
            call print_ex_2D(plot_title,plot_name,f_F,draw=.false.)
            call draw_ex(plot_title,plot_name,2,1,.false.)
            if (print_log) then
                plot_title = ['cos [log]','sin [log]']
                call print_ex_2D(plot_title,trim(plot_name)//'_log',&
                    &log10(abs(f_F)),draw=.false.)
                call draw_ex(plot_title,trim(plot_name)//'_log',2,1,.false.)
            end if
        end if
    end function nufft
end module grid_utilities
