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
    public coord_F2E, coord_E2F, calc_XYZ_grid, calc_eqd_grid, extend_grid_F, &
        &calc_int_vol, copy_grid, trim_grid, untrim_grid, setup_deriv_data, &
        &setup_interp_data, apply_disc, calc_n_par_X_rich, calc_vec_comp, &
        &nufft, find_compr_range, calc_arc_angle
#if ldebug
    public debug_calc_int_vol, debug_calc_vec_comp
#endif
    
    ! global variables
#if ldebug
    logical :: debug_calc_int_vol = .false.                                     ! plot debug information for calc_int_vol
    logical :: debug_calc_vec_comp = .false.                                    ! plot debug information for calc_vec_comp
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
        use VMEC_utilities, only: fourier2real
        use VMEC_vars, only: L_V_c, L_V_s, is_asym_V
        
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
        dims = shape(theta_E)
        
        ! tests
        if (dims(3).ne.size(r_F) .or. dims(3).ne.size(r_E)) then
            ierr = 1
            err_msg = 'r_F and r_E need to have the correct dimensions'
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
            use VMEC_vars, only: mnmax_V
            
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
            
            ! set up guess
            allocate(theta_E_guess(dims(1),dims(2)))
            allocate(theta_E_guess_3D(dims(1),dims(2),dims(3)))
            
            ! set up guess: try theta_F - lambda(theta_F)
            ! (lambda(theta_V) would be the exact solution)
            ierr = fourier2real(L_V_c_loc,L_V_s_loc,theta_F,zeta_E,&
                &theta_E_guess_3D,sym=[is_asym_V,.true.],deriv=[0,0])
            CHCKERR('')
            theta_E_guess_3D = theta_F - theta_E_guess_3D
            
            ! the poloidal angle has to be found as the zero of
            !   f = theta_F - theta_E - lambda
            err_msg = calc_zero_HH(dims,theta_E,fun_pol,3,&
                &theta_E_guess_3D,relax_fac=0.99_dp)                            ! relax factor 1.0 gave problems sometimes
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
        dims = shape(theta_F)
        
        ! tests
        if (dims(3).ne.size(r_F) .or. dims(3).ne.size(r_E)) then
            ierr = 1
            err_msg = 'r_F and r_E need to have the correct dimensions'
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
            use VMEC_utilities, only: fourier2real
            use VMEC_vars, only: mnmax_V, L_V_c, L_V_s, is_asym_V
            
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
    ! Note: The  normalization factor R_0 for  length is taken into  account and
    ! the output is transformed back to unnormalized values:
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
            use VMEC_utilities, only: fourier2real
            use VMEC_vars, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, &
                &mnmax_V, is_asym_V
            
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
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            ! test
            if (.not.allocated(grid_XYZ%trigon_factors)) then
                ierr = 1
                err_msg = 'trigon factors of grid_XYZ need to be allocated'
                CHCKERR(err_msg)
            end if
            
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
            ! (the geometrical zeta is equal to the VMEC zeta)
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
    ! n_zeta_plot angular and own loc_n_r points in F coordinates.
    ! Optionally,  the grid  can  also  be converted  to  E  coordinates if  the
    ! equilibrium grid is  provided. This is required for  many operations, such
    ! as the calculation of X, Y and Z through calc_XYZ_grid.
    integer function extend_grid_F(grid_in,grid_ext,grid_eq,n_theta_plot,&
        &n_zeta_plot,lim_theta_plot,lim_zeta_plot) result(ierr)
        use num_vars, only: n_theta => n_theta_plot, n_zeta => n_zeta_plot, &
            &min_theta_plot, max_theta_plot, min_zeta_plot, max_zeta_plot
        
        character(*), parameter :: rout_name = 'extend_grid_F'
        
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
        
        ! user ouput
#if ldebug
        call writo('Theta limits: '//trim(r2strt(lim_theta_plot_loc(1)))//&
            &'pi .. '//trim(r2strt(lim_theta_plot_loc(2)))//'pi')
        call writo('Zeta limits: '//trim(r2strt(lim_zeta_plot_loc(1)))//&
            &'pi .. '//trim(r2strt(lim_zeta_plot_loc(2)))//'pi')
#endif
        
        ! creating  equilibrium  grid  for  the output  that  covers  the  whole
        ! geometry angularly in F coordinates
        ierr = grid_ext%init([n_theta_plot_loc,n_zeta_plot_loc,grid_in%n(3)],&
            &[grid_in%i_min,grid_in%i_max])
        CHCKERR('')
        grid_ext%r_F = grid_in%r_F
        grid_ext%loc_r_F = grid_in%loc_r_F
        ierr = calc_eqd_grid(grid_ext%theta_F,lim_theta_plot_loc(1)*pi,&
            &lim_theta_plot_loc(2)*pi,1)
        CHCKERR('')
        ierr = calc_eqd_grid(grid_ext%zeta_F,lim_zeta_plot_loc(1)*pi,&
            &lim_zeta_plot_loc(2)*pi,2)
        CHCKERR('')
        
        ! convert all F coordinates to E coordinates if requested
        if (present(grid_eq)) then
            grid_ext%r_F = grid_in%r_F
            ierr = coord_F2E(grid_eq,&
                &grid_ext%loc_r_F,grid_ext%theta_F,grid_ext%zeta_F,&
                &grid_ext%loc_r_E,grid_ext%theta_E,grid_ext%zeta_E)
            CHCKERR('')
        end if
    end function extend_grid_F
    
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
    ! next process looses one more point.
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
                    &ceiling((tot_i_max(rank)-tot_i_min(rank+1)+1._dp)/2))
            else
                i_lim_out(1) = tot_i_min(1)
            end if
            if (rank.lt.n_procs-1) then
                i_lim_out(2) = min(tot_i_max(n_procs),tot_i_max(rank+1)-&
                    &floor((tot_i_max(rank+1)-tot_i_min(rank+2)+1._dp)/2))
            else
                i_lim_out(2) = tot_i_max(n_procs)
            end if
            
            ! limit to input limits (necessary if a process has only one point)
            i_lim_out(1) = max(min(i_lim_out(1),grid_in%i_max),grid_in%i_min)
            i_lim_out(2) = max(min(i_lim_out(2),grid_in%i_max),grid_in%i_min)
            
            ! get  min_i's and max_i's  of the grid_out,  not shifted by  min of
            ! first process
            ierr = get_ser_var([i_lim_out(1)],tot_i_min,scatter=.true.)
            CHCKERR('')
            ierr = get_ser_var([i_lim_out(2)],tot_i_max,scatter=.true.)
            CHCKERR('')
            
            ! set n of output grid
            n_out(1) = grid_in%n(1)
            n_out(2) = grid_in%n(2)
            n_out(3) = tot_i_max(n_procs)-tot_i_min(1)+1
            
            ! create new grid
            ierr = grid_out%init(n_out,i_lim_out-tot_i_min(1)+1,divided=.true.) ! limits shifted by min of first process
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
        if (grid_in%n(1).ne.0 .and. grid_in%n(2).ne.0) then                     ! only if 3D grid
            grid_out%theta_E = grid_in%theta_E(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%zeta_E = grid_in%zeta_E(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%theta_F = grid_in%theta_F(:,:,i_lim_out(1):i_lim_out(2))
            grid_out%zeta_F = grid_in%zeta_F(:,:,i_lim_out(1):i_lim_out(2))
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
    
    ! Calculates  contra-  and covariant  components  of  a vector  in  multiple
    ! coordinate systems:
    !   1. F: (alpha,psi,theta),
    !   2. M: (psi,theta,zeta),
    !   3. - H: (psi,theta,phi) if HELENA is used,
    !      - V: (r,theta,zeta) if VMEC is used,
    !   4. C: (R,phi,Z),
    !   5. X: (X,Y,Z),
    ! as well  as the magnitude,  starting from  the input in  Flux coordinates.
    ! Note that the components in X, Y and Z can be used to plot the real vector
    ! and that covariant components should  be equal to contravariant components
    ! in this coordinate system.
    ! Also, the fluxes can be calculated and plot.
    ! By  default,  the cartesian  components  are  returned,  but this  can  be
    ! indicated differently.  Furthermore, the results for  different coordinate
    ! systems can be plotted by providing a base name.
    ! Note  that plots for different  Richardson levels can be  combined to show
    ! the total grid by just plotting them all individually.
    ! Note: The metric factors and transformation matrices have to be allocated.
    ! They  can be  calculated  using  the routines  from  eq_ops,  for deriv  =
    ! [0,0,0].
    ! Note: For VMEC,  the trigonometric factors of grid_XYZ  must be calculated
    ! beforehand. Optionally,  by providing X, Y  and Z, the ones  calculated in
    ! this  routine are  overwritten. This  is  usefull for,  for example,  slab
    ! geometries.
    ! Note: The normalization  factors are taken into account and  the output is
    ! transformed back to unnormalized values.
    integer function calc_vec_comp(grid,grid_eq,eq_1,eq_2,v_com,norm_disc_prec,&
        &v_mag,base_name,max_transf,v_flux_tor,v_flux_pol,XYZ,compare_tor_pos) &
        &result(ierr)
        
        use grid_vars, only: disc_type
        use eq_vars, only: eq_1_type, eq_2_type, max_flux_F
        use eq_utilities, only: calc_inv_met
        use num_vars, only: eq_jobs_lims, eq_job_nr, use_pol_flux_F, eq_style, &
            &use_normalization, rank, tol_zero, &
            &compare_tor_pos_glob => compare_tor_pos
        use num_utilities, only: c, calc_int
        use eq_vars, only: R_0, B_0, psi_0
        use VMEC_utilities, only: calc_trigon_factors
        use mpi_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'calc_vec_comp'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid on which vector is calculated
        type(grid_type), intent(in) :: grid_eq                                  ! grid on which equilibrium variables are calculated
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium quantities
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        real(dp), intent(inout) :: v_com(:,:,:,:,:)                             ! covariant and contravariant components of v (dim1,dim2,dim3,3,2)
        integer, intent(in) :: norm_disc_prec                                   ! precision for normal derivatives
        real(dp), intent(inout), optional :: v_mag(:,:,:)                       ! magnitude of v (dim1,dim2,dim3)
        character(len=*), intent(in), optional :: base_name                     ! base name for output plotting
        integer, intent(in), optional :: max_transf                             ! maximum transformation level (2: Magnetic, 3: Equilibrium, 4: Cylindrical, 5: Cartesian [def])
        real(dp), intent(inout), allocatable, optional :: v_flux_pol(:,:)       ! poloidal flux of v as function of normal coordinate for all toroidal positions
        real(dp), intent(inout), allocatable, optional :: v_flux_tor(:,:)       ! toroidal flux of v as function of normal coordinate for all poloidal positions
        real(dp), intent(in), optional :: XYZ(:,:,:,:)                          ! X, Y and Z of grid
        logical, intent(in), optional :: compare_tor_pos                        ! compare toroidal positions
        
        ! local variables
        type(grid_type) :: grid_trim                                            ! trimmed plot grid
        integer :: max_transf_loc                                               ! local max_transf
        integer :: id, jd, kd                                                   ! counter
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_offset(4)                                               ! local offset of plot
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: c_loc                                                        ! local c
        integer :: tor_id(2)                                                    ! toroidal indices
        logical :: cont_plot                                                    ! continued plot
        logical :: do_plot                                                      ! perform plotting
        logical :: compare_tor_pos_loc                                          ! local compare_tor_pos
        character(len=max_str_ln) :: description(3)                             ! description of plots
        character(len=max_str_ln) :: file_names(3)                              ! plot file names
        character(len=max_str_ln) :: var_names(3,2)                             ! variable names
        character(len=5) :: coord_names(3)                                      ! name of coordinates
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp) :: norm_len                                                    ! normalization factor for lengths, to cancel the one introduced in "calc_XYZ_grid"
        real(dp), allocatable :: q_saf(:,:)                                     ! interpolated q_saf_FD and derivative
        real(dp), allocatable :: jac(:,:,:)                                     ! interpolated jac_FD
        real(dp), allocatable :: XYZR(:,:,:,:)                                  ! X, Y, Z and R of surface in cylindrical coordinates, untrimmed grid
        real(dp), allocatable :: X(:,:,:,:), Y(:,:,:,:), Z(:,:,:,:)             ! copy of X, Y and Z, trimmed grid
        real(dp), allocatable :: v_temp(:,:,:,:,:)                              ! temporary variable for v
        real(dp), allocatable :: v_ser_temp(:)                                  ! temporary serial variable
        real(dp), allocatable :: v_ser_temp_int(:)                              ! temporary integrated serial variable
        real(dp), allocatable :: T_BA(:,:,:,:,:,:,:), T_AB(:,:,:,:,:,:,:)       ! transformation matrices from A to B
        real(dp), allocatable :: D1R(:,:,:), D2R(:,:,:)                         ! dR/dpsi, dR/dtheta in HELENA coords.
        real(dp), allocatable :: D1Z(:,:,:), D2Z(:,:,:)                         ! dZ/dpsi, dZ/dtheta in HELENA coords.
        type(disc_type) :: norm_deriv_data                                      ! data for normal derivation
        type(disc_type) :: pol_deriv_data                                       ! data for poloidal derivation
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Prepare calculation of vector components')
        call lvl_ud(1)
        
        ! set up maximum level up to which to transform and whether to plot
        max_transf_loc = 5
        if (present(max_transf)) then
            if (max_transf.ge.1 .and. max_transf.le.5) &
                &max_transf_loc = max_transf
        end if
        do_plot = .false.
        if (present(base_name)) then
            if (trim(base_name).ne.'') do_plot = .true.
        end if
        
        ! set up continued plot
        cont_plot = eq_job_nr.gt.1
        
        ! set up normalization factor
        norm_len = 1._dp
        if (use_normalization) norm_len = R_0
        
        ! trim plot grid
        ierr = trim_grid(grid,grid_trim,norm_id)
        CHCKERR('')
        
        ! if r starts at 0, take away first points as it is singular
        if (abs(grid_trim%r_F(1)).lt.tol_zero) then
            if (rank.eq.0) then
                norm_id(1) = 1+norm_disc_prec
                grid_trim%loc_n_r = grid_trim%loc_n_r-norm_disc_prec
            else
                grid_trim%i_min = grid_trim%i_min-norm_disc_prec
            end if
            grid_trim%n(3) = grid_trim%n(3)-norm_disc_prec
            grid_trim%i_max = grid_trim%i_max-norm_disc_prec
        end if
        
        ! setup normal interpolation data for equilibrium grid
        ierr = setup_interp_data(grid_eq%loc_r_F,grid%loc_r_F,&
            &norm_interp_data,norm_disc_prec)
        CHCKERR('')
        
        ! set up X, Y, Z and R in grid
        allocate(XYZR(grid%n(1),grid%n(2),grid%loc_n_r,4))
        ierr = calc_XYZ_grid(grid_eq,grid,XYZR(:,:,:,1),XYZR(:,:,:,2),&
            &XYZR(:,:,:,3),R=XYZR(:,:,:,4))
        CHCKERR('')
        
        ! set up local compare_tor_pos
        compare_tor_pos_loc = compare_tor_pos_glob
        if (present(compare_tor_pos)) compare_tor_pos_loc = compare_tor_pos
        
        ! if plotting is required
        if (do_plot) then
            ! set up plot dimensions and local dimensions
            plot_dim = [grid_trim%n,3]
            plot_offset = [0,0,grid_trim%i_min-1,0]
            tor_id = [1,size(v_com,2)]
            if (compare_tor_pos_loc) then
                if (plot_dim(2).ne.2) then
                    ierr = 1
                    err_msg = 'When comparing toroidal positions, need to &
                        &have 2 toroidal points'
                    CHCKERR(err_msg)
                end if
                plot_dim(2) = 1
                tor_id(2) = 1
            end if
            
            ! possibly modify if multiple equilibrium parallel jobs
            if (size(eq_jobs_lims,2).gt.1) then
                plot_dim(1) = eq_jobs_lims(2,size(eq_jobs_lims,2)) - &
                    &eq_jobs_lims(1,1) + 1
                plot_offset(1) = eq_jobs_lims(1,eq_job_nr) - 1
                if (compare_tor_pos_loc) then
                    ierr = 1
                    err_msg = 'When comparing toroidal positions, cannot have &
                        &multiple equilibrium jobs'
                    CHCKERR(err_msg)
                end if
            end if
            
            ! tests
            if (eq_style.eq.1 .and. .not.allocated(grid%trigon_factors)) then
                ierr = 1
                err_msg = 'trigonometric factors not allocated'
                CHCKERR(err_msg)
            end if
            
            ! user output
            if (compare_tor_pos_loc) call writo('Comparing toroidal positions')
            
            ! copy X, Y and Z to trimmed grid copies
            allocate(X(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,1))
            allocate(Y(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,1))
            allocate(Z(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,1))
            if (compare_tor_pos_loc) then
                if (present(XYZ)) then
                    X(:,1,:,1) = 0.5_dp*&
                        &(XYZ(:,1,norm_id(1):norm_id(2),1) + &
                        &XYZ(:,2,norm_id(1):norm_id(2),1))
                    Y(:,1,:,1) = 0.5_dp*&
                        &(XYZ(:,1,norm_id(1):norm_id(2),2) + &
                        &XYZ(:,2,norm_id(1):norm_id(2),2))
                    Z(:,1,:,1) = 0.5_dp*&
                        &(XYZ(:,1,norm_id(1):norm_id(2),3) + &
                        &XYZ(:,2,norm_id(1):norm_id(2),3))
                else
                    X(:,1,:,1) = 0.5_dp*&
                        &(XYZR(:,1,norm_id(1):norm_id(2),1) + &
                        &XYZR(:,2,norm_id(1):norm_id(2),1))
                    Y(:,1,:,1) = 0.5_dp*&
                        &(XYZR(:,1,norm_id(1):norm_id(2),2) + &
                        &XYZR(:,2,norm_id(1):norm_id(2),2))
                    Z(:,1,:,1) = 0.5_dp*&
                        &(XYZR(:,1,norm_id(1):norm_id(2),3) + &
                        &XYZR(:,2,norm_id(1):norm_id(2),3))
                end if
            else
                if (present(XYZ)) then
                    X(:,:,:,1) = XYZ(:,:,norm_id(1):norm_id(2),1)
                    Y(:,:,:,1) = XYZ(:,:,norm_id(1):norm_id(2),2)
                    Z(:,:,:,1) = XYZ(:,:,norm_id(1):norm_id(2),3)
                else
                    X(:,:,:,1) = XYZR(:,:,norm_id(1):norm_id(2),1)
                    Y(:,:,:,1) = XYZR(:,:,norm_id(1):norm_id(2),2)
                    Z(:,:,:,1) = XYZR(:,:,norm_id(1):norm_id(2),3)
                end if
            end if
        end if
        
        ! set up temporal interpolated copy of  v, T_BA and T_AB
        allocate(v_temp(grid%n(1),grid%n(2),grid%loc_n_r,3,2))
        allocate(T_BA(grid%n(1),grid%n(2),grid%loc_n_r,9,0:0,0:0,0:0))
        allocate(T_AB(grid%n(1),grid%n(2),grid%loc_n_r,9,0:0,0:0,0:0))
        allocate(q_saf(grid%loc_n_r,0:1))
        if ((present(v_flux_tor) .or. present(v_flux_pol)) .and. &
            &.not.compare_tor_pos_loc) then
            allocate(jac(grid%n(1),grid%n(2),grid%loc_n_r))
            if (rank.eq.0 .and. .not.cont_plot) then
                if (present(v_flux_tor)) then
                    allocate(v_flux_tor(grid_trim%n(3),plot_dim(2)))
                    v_flux_tor = 0._dp
                end if
                if (present(v_flux_pol)) then
                    allocate(v_flux_pol(grid_trim%n(3),plot_dim(1)))
                    v_flux_pol = 0._dp
                end if
            end if
        end if
        
        call lvl_ud(-1)
        
        ! 1. Flux coordinates (alpha,psi,theta)
        call writo('Flux coordinates')
        call lvl_ud(1)
        if (present(v_mag)) then
            v_mag = 0._dp
            do id = 1,3
                v_mag = v_mag + v_com(:,:,:,id,1)*v_com(:,:,:,id,2)
            end do
            v_mag = sqrt(v_mag)
        end if
        
        ! save temporary copy, normalized
        v_temp = v_com
        
        ! set up plot variables
        if (do_plot) then
            coord_names(1) = 'alpha'
            coord_names(2) = 'psi'
            if (use_pol_flux_F) then
                coord_names(3) = 'theta'
            else
                coord_names(3) = 'zeta'
            end if
            var_names = trim(base_name)
            do id = 1,3
                var_names(id,1) = trim(var_names(id,1))//'_sub_'//&
                    &trim(coord_names(id))
                var_names(id,2) = trim(var_names(id,2))//'_sup_'//&
                    &trim(coord_names(id))
            end do
            description(1) = 'covariant components of the magnetic field in &
                &Flux coordinates'
            description(2) = 'contravariant components of the magnetic field &
                &in Flux coordinates'
            description(3) = 'magnitude of the magnetic field in Flux &
                &coordinates'
            file_names(1) = trim(base_name)//'_F_sub'
            file_names(2) = trim(base_name)//'_F_sup'
            file_names(3) = trim(base_name)//'_F_mag'
            if (compare_tor_pos_loc) then
                do id = 1,3
                    file_names(id) = trim(file_names(id))//'_COMP'
                end do
            end if
            if (use_normalization) then
                v_com(:,:,:,1,1) = v_com(:,:,:,1,1) * R_0                       ! norm factor for e_alpha
                v_com(:,:,:,2,1) = v_com(:,:,:,2,1) / (R_0*B_0)                 ! norm factor for e_psi
                v_com(:,:,:,3,1) = v_com(:,:,:,3,1) * R_0                       ! norm factor for e_theta
                v_com(:,:,:,1,2) = v_com(:,:,:,1,2) / R_0                       ! norm factor for e^alpha
                v_com(:,:,:,2,2) = v_com(:,:,:,2,2) * (R_0*B_0)                 ! norm factor for e^psi
                v_com(:,:,:,3,2) = v_com(:,:,:,3,2) / R_0                       ! norm factor for e^theta
            end if
            if (compare_tor_pos_loc) v_com(:,1,:,:,:) = 2._dp*&
                &(v_com(:,2,:,:,:) - v_com(:,1,:,:,:))/&
                &(v_com(:,2,:,:,:) + v_com(:,1,:,:,:))
            do id = 1,2
                call plot_HDF5(var_names(:,id),trim(file_names(id)),&
                    &v_com(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2),:,id),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=X(:,tor_id(1):tor_id(2),:,:),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,:),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,:),&
                    &cont_plot=cont_plot,description=description(id))
            end do
            if (present(v_mag)) then
                if (compare_tor_pos_loc) v_mag(:,1,:) = 2._dp*&
                    &(v_mag(:,2,:) - v_mag(:,1,:))/&
                    &(v_mag(:,2,:) + v_mag(:,1,:))
                call plot_HDF5(trim(base_name),trim(file_names(3)),&
                    &v_mag(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3),&
                    &X=X(:,tor_id(1):tor_id(2),:,1),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,1),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,1),&
                    &cont_plot=cont_plot,description=description(3))
            end if
        end if
        
        call lvl_ud(-1)
        
        ! 2.   Magnetic  coordinates   (phi,theta,zeta)   by  converting   using
        ! transformation matrices:
        !   (v_i)_M = T_M^F (v_i)_F
        !   (v^i)_M = (v^i)_F T_F^M
        ! where
        !             ( -q' theta 1 0 )
        !   T_MF    = (    -q     0 1 )
        !             (     1     0 0 )
        ! if the poloidal flux is used, or
        !             ( -q'/q zeta q 0 )
        !   T_MF    = (    -1      0 0 )
        !             (    1/q     0 1 )
        ! if it is the toroidal flux. Its inverse can be calculated as well.
        call writo('magnetic coordinates')
        call lvl_ud(1)
        T_BA = 0._dp
        T_AB = 0._dp
        ierr = apply_disc(eq_1%q_saf_FD(:,0:1),norm_interp_data,q_saf,1)
        CHCKERR('')
        c_loc = c([1,1],.false.)
        if (use_pol_flux_F) then
            T_BA(:,:,:,c_loc,0,0,0) = -grid%theta_F
            do kd = 1,grid%loc_n_r
                T_BA(:,:,kd,c_loc,0,0,0) = T_BA(:,:,kd,c_loc,0,0,0)*q_saf(kd,1)
                T_BA(:,:,kd,c([2,1],.false.),0,0,0) = -q_saf(kd,0)
            end do
            T_BA(:,:,:,c([3,1],.false.),0,0,0) = 1._dp
            T_BA(:,:,:,c([1,2],.false.),0,0,0) = 1._dp
            T_BA(:,:,:,c([2,3],.false.),0,0,0) = 1._dp
        else
            T_BA(:,:,:,c_loc,0,0,0) = -grid%zeta_F
            do kd = 1,grid%loc_n_r
                T_BA(:,:,kd,c_loc,0,0,0) = &
                    &T_BA(:,:,kd,c_loc,0,0,0)*q_saf(kd,1)/q_saf(kd,0)
                T_BA(:,:,kd,c([2,1],.false.),0,0,0) = 1._dp/q_saf(kd,0)
                T_BA(:,:,kd,c([1,2],.false.),0,0,0) = q_saf(kd,0)
            end do
            T_BA(:,:,:,c([2,1],.false.),0,0,0) = -1._dp
            T_BA(:,:,:,c([3,3],.false.),0,0,0) = 1._dp
        end if
        ierr = calc_inv_met(T_AB,T_BA,[0,0,0])
        CHCKERR('')
        v_com = 0._dp
        do jd = 1,3
            do id = 1,3
                v_com(:,:,:,jd,1) = v_com(:,:,:,jd,1) + T_BA(:,:,:,&
                    &c([jd,id],.false.),0,0,0) * v_temp(:,:,:,id,1)
                v_com(:,:,:,jd,2) = v_com(:,:,:,jd,2) + T_AB(:,:,:,&
                    &c([id,jd],.false.),0,0,0) * v_temp(:,:,:,id,2)
            end do
        end do
        if (present(v_mag)) then
            v_mag = 0._dp
            do id = 1,3
                v_mag = v_mag + v_com(:,:,:,id,1)*v_com(:,:,:,id,2)
            end do
            v_mag = sqrt(v_mag)
        end if
        
        ! set up plot variables
        if (do_plot) then
            coord_names(1) = 'psi'
            coord_names(2) = 'theta'
            coord_names(3) = 'zeta'
            var_names = trim(base_name)
            do id = 1,3
                var_names(id,1) = trim(var_names(id,1))//'_sub_'//&
                    &trim(coord_names(id))
                var_names(id,2) = trim(var_names(id,2))//'_sup_'//&
                    &trim(coord_names(id))
            end do
            description(1) = 'covariant components of the magnetic field in &
                &Magnetic coordinates'
            description(2) = 'contravariant components of the magnetic field &
                &in Magnetic coordinates'
            description(3) = 'magnitude of the magnetic field in Magnetic &
                &coordinates'
            file_names(1) = trim(base_name)//'_M_sub'
            file_names(2) = trim(base_name)//'_M_sup'
            file_names(3) = trim(base_name)//'_M_mag'
            if (compare_tor_pos_loc) then
                do id = 1,3
                    file_names(id) = trim(file_names(id))//'_COMP'
                end do
            end if
            if (use_normalization) then
                v_com(:,:,:,1,1) = v_com(:,:,:,1,1) / (R_0*B_0)                 ! norm factor for e_psi
                v_com(:,:,:,2,1) = v_com(:,:,:,2,1) * R_0                       ! norm factor for e_theta
                v_com(:,:,:,3,1) = v_com(:,:,:,3,1) * R_0                       ! norm factor for e_zeta
                v_com(:,:,:,1,2) = v_com(:,:,:,1,2) * (R_0*B_0)                 ! norm factor for e^psi
                v_com(:,:,:,2,2) = v_com(:,:,:,2,2) / R_0                       ! norm factor for e^theta
                v_com(:,:,:,3,2) = v_com(:,:,:,3,2) / R_0                       ! norm factor for e^zeta
            end if
            if (compare_tor_pos_loc) v_com(:,1,:,:,:) = 2._dp*&
                &(v_com(:,2,:,:,:) - v_com(:,1,:,:,:))/&
                &(v_com(:,2,:,:,:) + v_com(:,1,:,:,:))
            do id = 1,2
                call plot_HDF5(var_names(:,id),trim(file_names(id)),&
                    &v_com(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2),:,id),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=X(:,tor_id(1):tor_id(2),:,:),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,:),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,:),&
                    &cont_plot=cont_plot,description=description(id))
            end do
            if (present(v_mag)) then
                if (compare_tor_pos_loc) v_mag(:,1,:) = 2._dp*&
                    &(v_mag(:,2,:) - v_mag(:,1,:))/&
                    &(v_mag(:,2,:) + v_mag(:,1,:))
                call plot_HDF5(trim(base_name),trim(file_names(3)),&
                    &v_mag(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3),&
                    &X=X(:,tor_id(1):tor_id(2),:,1),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,1),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,1),&
                    &cont_plot=cont_plot,description=description(3))
            end if
        end if
        
        if ((present(v_flux_tor) .or. present(v_flux_pol)) .and. &
            &.not.compare_tor_pos_loc) then
            ! tests
            if (maxval(grid%theta_F(grid%n(1),:,:)).lt.&
                &minval(grid%theta_F(1,:,:))) then
                call writo('theta of the grid monotomically decreases in Flux &
                    &quantities.',alert=.true.)
                call lvl_ud(1)
                call writo('This inverts the sign of the toroidal flux.')
                call writo('Remember that the grid limits are prescribed &
                    &in Flux quantities.')
                call lvl_ud(-1)
            end if
            if (maxval(grid%zeta_F(:,grid%n(2),:)).lt.&
                &minval(grid%zeta_F(:,1,:))) then
                call writo('zeta of the grid monotomically decreases in Flux &
                    &quantities.',alert=.true.)
                call lvl_ud(1)
                call writo('This inverts the sign of the poloidal flux.')
                call writo('Remember that the grid limits are prescribed &
                    &in Flux quantities.')
                if (eq_style.eq.1) call writo('For VMEC, these are inverse.')
                call lvl_ud(-1)
            end if
            if (grid_trim%r_F(grid_trim%n(3)).lt.grid_trim%r_F(1)) then
                call writo('r of the grid monotomically decreases in Flux &
                    &quantities.',alert=.true.)
                call lvl_ud(1)
                call writo('This inverts the sign of the fluxes.')
                call writo('Remember that the grid limits are prescribed &
                    &in Flux quantities.')
                call lvl_ud(-1)
            end if
            if (abs(minval(grid_trim%r_F)).gt.tol_zero) then
                call writo('r of the grid does not start at zero.',alert=.true.)
                call lvl_ud(1)
                call writo('This leaves out part of the fluxes.')
                call lvl_ud(-1)
            else
                call writo('r of the grid starts at zero.', alert=.true.)
                call lvl_ud(1)
                call writo('To avoid singularities, the first '//&
                    &trim(i2str(norm_disc_prec))//' normal points are left out')
                call writo('This has some effect on the total integrated &
                    &fluxes')
                call lvl_ud(-1)
            end if
            
            ! initialize temporary variable
            allocate(v_ser_temp_int(grid_trim%n(3)))
            
            ! set up plot variables
            var_names(1,2) = 'integrated poloidal flux of '//trim(base_name)
            var_names(2,2) = 'integrated toroidal flux of '//trim(base_name)
            file_names(1) = trim(base_name)//'_flux_pol_int'
            file_names(2) = trim(base_name)//'_flux_tor_int'
            
#if ldebug
            if (debug_calc_vec_comp) then
                ! user output
                call writo('Instead of calculating fluxes, returning &
                    &integrals in the angular variables.',alert=.true.)
                call lvl_ud(1)
                    call writo('Useful to check whether Maxwell holds.')
                    call writo('i.e. whether loop integral of J returns &
                        &B flux')
                    call writo('Note that the J-variables now refer to B and &
                        &vice-versa')
                    call writo('Don''t forget the contribution of the toroidal &
                        &field B_zeta on axis, times 2piR.')
                    call writo('Don''t forget the vacuum permeability!')
                call lvl_ud(-1)
            else
#endif
                ! multiply angular contravariant coordinates by Jacobian
                ierr = apply_disc(eq_2%jac_FD(:,:,:,0,0,0),norm_interp_data,&
                    &jac,3)
                CHCKERR('')
                if (use_normalization) jac = jac*R_0/B_0
                v_com(:,:,:,2,2) = v_com(:,:,:,2,2)*jac
                v_com(:,:,:,3,2) = v_com(:,:,:,3,2)*jac
                deallocate(jac)
#if ldebug
            end if
#endif
            
            ! integrate toroidal flux
            if (present(v_flux_tor) .and. grid_trim%n(1).gt.1 .and. &
                &.not.compare_tor_pos_loc) then                                 ! integrate poloidally and normally
                ! loop over all toroidal points
                do jd = 1,grid_trim%n(2)
#if ldebug
                    if (debug_calc_vec_comp) then
                        ! for  all  normal   points,  integrate  poloidally  the
                        ! covariant poloidal quantity
                        do kd = norm_id(1),norm_id(2)
                            ierr = calc_int(v_com(:,jd,kd,2,1),&
                                &grid%theta_F(:,jd,kd),v_com(:,jd,kd,2,2))      ! save in contravariant variable
                            CHCKERR('')
                        end do
                        
                        ! gather on master
                        ierr = get_ser_var(v_com(grid_trim%n(1),jd,&
                            &norm_id(1):norm_id(2),2,2),v_ser_temp)
                        CHCKERR('')
                        
                        ! update integrals if master
                        if (rank.eq.0) then
                            if (use_pol_flux_F) then
                                v_flux_tor(:,jd) = v_flux_tor(:,jd) + v_ser_temp
                            else
                                v_flux_tor(:,jd+plot_offset(1)) = v_ser_temp
                            end if
                        end if
                    else
#endif
                        ! for all normal points, integrate poloidally
                        do kd = norm_id(1),norm_id(2)
                            ierr = calc_int(v_com(:,jd,kd,3,2),&
                                &grid%theta_F(:,jd,kd),v_com(:,jd,kd,3,1))      ! save in covariant variable
                            CHCKERR('')
                        end do
                        
                        ! gather on master
                        ierr = get_ser_var(v_com(grid_trim%n(1),jd,&
                            &norm_id(1):norm_id(2),3,1),v_ser_temp)
                        CHCKERR('')
                        
                        ! integrate normally and update integrals if master
                        if (rank.eq.0) then
                            ierr = calc_int(v_ser_temp,&
                                &grid_trim%r_F(norm_id(1):),v_ser_temp_int)
                            CHCKERR('')
                            if (use_normalization) v_ser_temp_int = &
                                &v_ser_temp_int*psi_0
                            if (use_pol_flux_F) then
                                v_flux_tor(:,jd) = v_flux_tor(:,jd) + &
                                    &v_ser_temp_int
                            else
                                v_flux_tor(:,jd+plot_offset(1)) = v_ser_temp_int
                            end if
                        end if
#if ldebug
                    end if
#endif
                end do
                
                ! plot output
                if (rank.eq.0 .and. do_plot .and. &
                    &eq_job_nr.eq.size(eq_jobs_lims,2)) then
                    call print_ex_2D([var_names(2,2)],&
                        &trim(file_names(2)),v_flux_tor,&
                        &x=reshape(grid_trim%r_F(norm_id(1):)*2*pi/max_flux_F,&
                        &[grid_trim%n(3),1]),draw=.false.)
                    call draw_ex([var_names(2,2)],trim(file_names(2)),&
                        &plot_dim(2),1,.false.)
                end if
            end if
            
            ! integrate poloidal flux
            if (present(v_flux_pol) .and. grid_trim%n(2).gt.1 .and. &
                &.not.compare_tor_pos_loc) then                                 ! integrate toroidally and normally
                ! loop over all poloidal points
                do id = 1,grid_trim%n(1)
#if ldebug
                    if (debug_calc_vec_comp) then
                        ! for  all  normal   points,  integrate  toroidally  the
                        ! covariant toroidal quantity
                        do kd = norm_id(1),norm_id(2)
                            ierr = calc_int(v_com(id,:,kd,3,1),&
                                &grid%zeta_F(id,:,kd),v_com(id,:,kd,3,2))       ! save in contravariant variable
                            CHCKERR('')
                        end do
                        
                        ! gather on master
                        ierr = get_ser_var(v_com(id,grid_trim%n(2),&
                            &norm_id(1):norm_id(2),3,2),v_ser_temp)
                        CHCKERR('')
                        
                        ! update integrals if master
                        if (rank.eq.0) then
                            if (use_pol_flux_F) then
                                v_flux_pol(:,id+plot_offset(1)) = v_ser_temp
                            else
                                v_flux_pol(:,id) = v_flux_pol(:,id) + v_ser_temp
                            end if
                        end if
                    else
#endif
                        ! for all normal points, integrate toroidally
                        do kd = norm_id(1),norm_id(2)
                            ierr = calc_int(v_com(id,:,kd,2,2),&
                                &grid%zeta_F(id,:,kd),v_com(id,:,kd,2,1))       ! save in covariant variable
                            CHCKERR('')
                        end do
                        
                        ! gather on master
                        ierr = get_ser_var(v_com(id,grid_trim%n(2),&
                            &norm_id(1):norm_id(2),2,1),v_ser_temp)
                        CHCKERR('')
                        
                        ! integrate normally and update integrals if master
                        if (rank.eq.0) then
                            ierr = calc_int(v_ser_temp,&
                                &grid_trim%r_F(norm_id(1):),v_ser_temp_int)
                            CHCKERR('')
                            if (use_normalization) v_ser_temp_int = &
                                &v_ser_temp_int*psi_0
                            if (use_pol_flux_F) then
                                v_flux_pol(:,id+plot_offset(1)) = v_ser_temp_int
                            else
                                v_flux_pol(:,id) = v_flux_pol(:,id) + &
                                    &v_ser_temp_int
                            end if
                        end if
#if ldebug
                    end if
#endif
                end do
                
                ! plot output
                if (rank.eq.0 .and. do_plot .and. &
                    &eq_job_nr.eq.size(eq_jobs_lims,2)) then
                    call print_ex_2D([var_names(1,2)],&
                        &trim(file_names(1)),v_flux_pol,&
                        &x=reshape(grid_trim%r_F(norm_id(1):)*2*pi/max_flux_F,&
                        &[grid_trim%n(3),1]),draw=.false.)
                    call draw_ex([var_names(1,2)],trim(file_names(1)),&
                        &plot_dim(1),1,.false.)
                end if
            end if
            
            ! clean up
            deallocate(v_ser_temp)
            if (rank.eq.0) deallocate(v_ser_temp_int)
        end if
        
        call lvl_ud(-1)
        
        if (max_transf_loc.eq.2) return
        
        ! 3.   HELENA   coordinates   (psi,theta,zeta)   or   VMEC   coordinates
        ! (r,theta,phi), by converting using transformation matrices
        !   (v_i)_H = T_H^F (v_i)_F
        !   (v^i)_H = (v^i)_F T_F^H
        ! Note: This is  done directly from step 1; The results  from step 2 are
        ! skipped, i.e. v_temp is not overwritten after step 2.
        call writo('Equilibrium coordinates')
        call lvl_ud(1)
        ierr = apply_disc(eq_2%T_FE(:,:,:,:,0,0,0),norm_interp_data,&
            &T_AB(:,:,:,:,0,0,0),3)
        CHCKERR('')
        ierr = apply_disc(eq_2%T_EF(:,:,:,:,0,0,0),norm_interp_data,&
            &T_BA(:,:,:,:,0,0,0),3)
        CHCKERR('')
        v_com = 0._dp
        do jd = 1,3
            do id = 1,3
                v_com(:,:,:,jd,1) = v_com(:,:,:,jd,1) + &
                    &T_BA(:,:,:,c([jd,id],.false.),0,0,0) * &
                    &v_temp(:,:,:,id,1)
                v_com(:,:,:,jd,2) = v_com(:,:,:,jd,2) + &
                    &T_AB(:,:,:,c([id,jd],.false.),0,0,0) * &
                    &v_temp(:,:,:,id,2)
            end do
        end do
        if (present(v_mag)) then
            v_mag = 0._dp
            do id = 1,3
                v_mag = v_mag + v_com(:,:,:,id,1)*v_com(:,:,:,id,2)
            end do
            v_mag = sqrt(v_mag)
        end if
        
        ! save temporary copy, normalized
        v_temp = v_com
        
        ! set up plot variables
        if (do_plot) then
            select case (eq_style)
                case (1)                                                        ! VMEC
                    coord_names(1) = 'r'
                    coord_names(2) = 'theta'
                    coord_names(3) = 'phi'
                    description(1) = 'covariant components of the magnetic &
                        &field in VMEC coordinates'
                    description(2) = 'contravariant components of the magnetic &
                        &field in VMEC coordinates'
                    description(3) = 'magnitude of the magnetic field in VMEC &
                        &coordinates'
                    file_names(1) = trim(base_name)//'_V_sub'
                    file_names(2) = trim(base_name)//'_V_sup'
                    file_names(3) = trim(base_name)//'_V_mag'
                case (2)                                                        ! HELENA
                    coord_names(1) = 'psi'
                    coord_names(2) = 'theta'
                    coord_names(3) = 'phi'
                    description(1) = 'covariant components of the magnetic &
                        &field in HELENA coordinates'
                    description(2) = 'contravariant components of the magnetic &
                        &field in HELENA coordinates'
                    description(3) = 'magnitude of the magnetic field in &
                        &HELENA coordinates'
                    file_names(1) = trim(base_name)//'_H_sub'
                    file_names(2) = trim(base_name)//'_H_sup'
                    file_names(3) = trim(base_name)//'_H_mag'
            end select
            if (compare_tor_pos_loc) then
                do id = 1,3
                    file_names(id) = trim(file_names(id))//'_COMP'
                end do
            end if
            var_names = trim(base_name)
            do id = 1,3
                var_names(id,1) = trim(var_names(id,1))//'_sub_'//&
                    &trim(coord_names(id))
                var_names(id,2) = trim(var_names(id,2))//'_sup_'//&
                    &trim(coord_names(id))
            end do
            if (use_normalization) then
                v_com(:,:,:,1,1) = v_com(:,:,:,1,1) * R_0                       ! norm factor for e_r or e_psi
                v_com(:,:,:,2,1) = v_com(:,:,:,2,1) * R_0                       ! norm factor for e_theta
                v_com(:,:,:,3,1) = v_com(:,:,:,3,1) * R_0                       ! norm factor for e_zeta or e_phi
                v_com(:,:,:,1,2) = v_com(:,:,:,1,2) / R_0                       ! norm factor for e^r or e^psi
                v_com(:,:,:,2,2) = v_com(:,:,:,2,2) / R_0                       ! norm factor for e^theta
                v_com(:,:,:,3,2) = v_com(:,:,:,3,2) / R_0                       ! norm factor for e^zeta or e^phi
            end if
            if (compare_tor_pos_loc) v_com(:,1,:,:,:) = 2._dp*&
                &(v_com(:,2,:,:,:) - v_com(:,1,:,:,:))/&
                &(v_com(:,2,:,:,:) + v_com(:,1,:,:,:))
            do id = 1,2
                call plot_HDF5(var_names(:,id),trim(file_names(id)),&
                    &v_com(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2),:,id),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=X(:,tor_id(1):tor_id(2),:,:),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,:),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,:),&
                    &cont_plot=cont_plot,description=description(id))
            end do
            if (present(v_mag)) then
                if (compare_tor_pos_loc) v_mag(:,1,:) = 2._dp*&
                    &(v_mag(:,2,:) - v_mag(:,1,:))/&
                    &(v_mag(:,2,:) + v_mag(:,1,:))
                call plot_HDF5(trim(base_name),trim(file_names(3)),&
                    &v_mag(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3),&
                    &X=X(:,tor_id(1):tor_id(2),:,1),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,1),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,1),&
                    &cont_plot=cont_plot,description=description(3))
            end if
        end if
        
        call lvl_ud(-1)
        
        if (max_transf_loc.eq.3) return
        
        ! 4.   cylindrical    coordinates   (R,phi,Z)   by    converting   using
        ! transformation matrices:
        !   (v_i)_C = T_C^E (v_i)_E
        !   (v^i)_C = (v^i)_E T_E^C
        ! where for VMEC (E=V), the transformation matrix T_VC is tabulated, and
        ! its inverse  is calculate,  and for  HELENA (E=H),  the transformation
        ! matrix T_HC is given by
        !             (    R'    0    Z'   )
        !   T_HC    = ( R_theta  0 Z_theta )
        !             (    0    -1    0    )
        ! and its inverse can be calculated as well.
        call writo('Cylindrical coordinates')
        call lvl_ud(1)
        T_BA = 0._dp
        T_AB = 0._dp
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = apply_disc(eq_2%T_VC(:,:,:,:,0,0,0),norm_interp_data,&
                    &T_AB(:,:,:,:,0,0,0),3)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! setup D1R and D2R, and same thing for Z
                ! (use untrimmed grid, with ghost region)
                allocate(D1R(grid%n(1),grid%n(2),grid%loc_n_r))
                allocate(D2R(grid%n(1),grid%n(2),grid%loc_n_r))
                allocate(D1Z(grid%n(1),grid%n(2),grid%loc_n_r))
                allocate(D2Z(grid%n(1),grid%n(2),grid%loc_n_r))
                ierr = setup_deriv_data(grid%loc_r_E,norm_deriv_data,1,&
                    &norm_disc_prec)
                CHCKERR('')
                ierr = apply_disc(XYZR(:,:,:,4)/norm_len,norm_deriv_data,D1R,3)
                CHCKERR('')
                ierr = apply_disc(XYZR(:,:,:,3)/norm_len,norm_deriv_data,D1Z,3)
                CHCKERR('')
                call norm_deriv_data%dealloc()
                do kd = 1,grid%loc_n_r
                    do jd = 1,grid%n(2)
                        ierr = setup_deriv_data(grid%theta_E(:,jd,kd),&
                            &pol_deriv_data,1,norm_disc_prec)
                        CHCKERR('')
                        ierr = apply_disc(XYZR(:,jd,kd,4)/norm_len,&
                            &pol_deriv_data,D2R(:,jd,kd))
                        CHCKERR('')
                        ierr = apply_disc(XYZR(:,jd,kd,3)/norm_len,&
                            &pol_deriv_data,D2Z(:,jd,kd))
                        CHCKERR('')
                        call pol_deriv_data%dealloc()
                    end do
                end do
                T_AB(:,:,:,c([1,1],.false.),0,0,0) = D1R
                T_AB(:,:,:,c([1,3],.false.),0,0,0) = D1Z
                T_AB(:,:,:,c([2,1],.false.),0,0,0) = D2R
                T_AB(:,:,:,c([2,3],.false.),0,0,0) = D2Z
                T_AB(:,:,:,c([3,2],.false.),0,0,0) = -1._dp
                deallocate(D1R,D2R,D1Z,D2Z)
        end select
        ierr = calc_inv_met(T_BA,T_AB,[0,0,0])
        CHCKERR('')
        v_com = 0._dp
        do jd = 1,3
            do id = 1,3
                v_com(:,:,:,jd,1) = v_com(:,:,:,jd,1) + T_BA(:,:,:,&
                    &c([jd,id],.false.),0,0,0) * v_temp(:,:,:,id,1)
                v_com(:,:,:,jd,2) = v_com(:,:,:,jd,2) + T_AB(:,:,:,&
                    &c([id,jd],.false.),0,0,0) * v_temp(:,:,:,id,2)
            end do
        end do
        if (present(v_mag)) then
            v_mag = 0._dp
            do id = 1,3
                v_mag = v_mag + v_com(:,:,:,id,1)*v_com(:,:,:,id,2)
            end do
            v_mag = sqrt(v_mag)
        end if
        
        ! save temporary copy, normalized
        v_temp = v_com
        
        ! set up plot variables
        if (do_plot) then
            description(1) = 'covariant components of the magnetic field in &
                &cylindrical coordinates'
            description(2) = 'contravariant components of the magnetic field &
                &in cylindrical coordinates'
            description(3) = 'magnitude of the magnetic field in cylindrical &
                &coordinates'
            file_names(1) = trim(base_name)//'_C_sub'
            file_names(2) = trim(base_name)//'_C_sup'
            file_names(3) = trim(base_name)//'_C_mag'
            if (compare_tor_pos_loc) then
                do id = 1,3
                    file_names(id) = trim(file_names(id))//'_COMP'
                end do
            end if
            coord_names(1) = 'R'
            coord_names(2) = 'phi'
            coord_names(3) = 'Z'
            var_names = trim(base_name)
            do id = 1,3
                var_names(id,1) = trim(var_names(id,1))//'_sub_'//&
                    &trim(coord_names(id))
                var_names(id,2) = trim(var_names(id,2))//'_sup_'//&
                    &trim(coord_names(id))
            end do
            if (use_normalization) then
                !v_com(:,:,:,1,1) = v_com(:,:,:,1,1)                             ! norm factor for e_R
                v_com(:,:,:,2,1) = v_com(:,:,:,2,1) * R_0                       ! norm factor for e_phi
                !v_com(:,:,:,3,1) = v_com(:,:,:,3,1)                             ! norm factor for e_Z
                !v_com(:,:,:,1,2) = v_com(:,:,:,1,2)                             ! norm factor for e^R
                v_com(:,:,:,2,2) = v_com(:,:,:,2,2) / R_0                       ! norm factor for e^phi
                !v_com(:,:,:,3,2) = v_com(:,:,:,3,2)                             ! norm factor for e^Z
            end if
            if (compare_tor_pos_loc) v_com(:,1,:,:,:) = 2._dp*&
                &(v_com(:,2,:,:,:) - v_com(:,1,:,:,:))/&
                &(v_com(:,2,:,:,:) + v_com(:,1,:,:,:))
            do id = 1,2
                call plot_HDF5(var_names(:,id),trim(file_names(id)),&
                    &v_com(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2),:,id),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=X(:,tor_id(1):tor_id(2),:,:),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,:),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,:),&
                    &cont_plot=cont_plot,description=description(id))
            end do
            if (present(v_mag)) then
                if (compare_tor_pos_loc) v_mag(:,1,:) = 2._dp*&
                    &(v_mag(:,2,:) - v_mag(:,1,:))/&
                    &(v_mag(:,2,:) + v_mag(:,1,:))
                call plot_HDF5(trim(base_name),trim(file_names(3)),&
                    &v_mag(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3),&
                    &X=X(:,tor_id(1):tor_id(2),:,1),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,1),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,1),&
                    &cont_plot=cont_plot,description=description(3))
            end if
        end if
        
        call lvl_ud(-1)
        
        if (max_transf_loc.eq.4) return
        
        ! 5. Cartesian  coordinates (X,Y,Z)  by converting  using transformation
        ! matrices
        !    (e_R)    (   cos(phi)   sin(phi)  0 ) (e_X)
        !   (e_phi) = ( -R sin(phi) R cos(phi) 0 ) (e_Y)
        !    (e_Z)    (     0          0       1 ) (e_Z)
        ! with  phi  =  -  zeta_F.   For  the  transformation  of  contravariant
        ! components, the factor R becomes 1/R and the transpose is taken.
        call writo('Cartesian coordinates')
        call lvl_ud(1)
        T_BA = 0._dp
        T_AB = 0._dp
        T_AB(:,:,:,c([1,1],.false.),0,0,0) = cos(-grid%zeta_F)
        T_AB(:,:,:,c([1,2],.false.),0,0,0) = sin(-grid%zeta_F)
        T_AB(:,:,:,c([2,1],.false.),0,0,0) = -XYZR(:,:,:,4)/norm_len*&
            &sin(-grid%zeta_F)
        T_AB(:,:,:,c([2,2],.false.),0,0,0) = XYZR(:,:,:,4)/norm_len*&
            &cos(-grid%zeta_F)
        T_AB(:,:,:,c([3,3],.false.),0,0,0) = 1._dp
        ierr = calc_inv_met(T_BA,T_AB,[0,0,0])
        CHCKERR('')
        v_com = 0._dp
        do jd = 1,3
            do id = 1,3
                v_com(:,:,:,jd,1) = v_com(:,:,:,jd,1) + T_BA(:,:,:,&
                    &c([jd,id],.false.),0,0,0) * v_temp(:,:,:,id,1)
                v_com(:,:,:,jd,2) = v_com(:,:,:,jd,2) + T_AB(:,:,:,&
                    &c([id,jd],.false.),0,0,0) * v_temp(:,:,:,id,2)
            end do
        end do
        if (present(v_mag)) then
            v_mag = 0._dp
            do id = 1,3
                v_mag = v_mag + v_com(:,:,:,id,1)*v_com(:,:,:,id,2)
            end do
            v_mag = sqrt(v_mag)
        end if
        
        ! set up plot variables
        if (do_plot) then
            description(1) = 'covariant components of the magnetic field in &
                &cartesian coordinates'
            description(2) = 'contravariant components of the magnetic field &
                &in cartesian coordinates'
            description(3) = 'magnitude of the magnetic field in cartesian &
                &coordinates'
            file_names(1) = trim(base_name)//'_X_sub'
            file_names(2) = trim(base_name)//'_X_sup'
            file_names(3) = trim(base_name)//'_X_mag'
            if (compare_tor_pos_loc) then
                do id = 1,3
                    file_names(id) = trim(file_names(id))//'_COMP'
                end do
            end if
            coord_names(1) = 'X'
            coord_names(2) = 'Y'
            coord_names(3) = 'Z'
            var_names = trim(base_name)
            do id = 1,3
                var_names(id,1) = trim(var_names(id,1))//'_sub_'//&
                    &trim(coord_names(id))
                var_names(id,2) = trim(var_names(id,2))//'_sup_'//&
                    &trim(coord_names(id))
            end do
            if (compare_tor_pos_loc) v_com(:,1,:,:,:) = 2._dp*&
                &(v_com(:,2,:,:,:) - v_com(:,1,:,:,:))/&
                &(v_com(:,2,:,:,:) + v_com(:,1,:,:,:))
            do id = 1,2
                call plot_HDF5(var_names(:,id),trim(file_names(id)),&
                    &v_com(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2),:,id),&
                    &tot_dim=plot_dim,loc_offset=plot_offset,&
                    &X=X(:,tor_id(1):tor_id(2),:,:),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,:),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,:),&
                    &cont_plot=cont_plot,description=description(id))
            end do
            if (present(v_mag)) then
                if (compare_tor_pos_loc) v_mag(:,1,:) = 2._dp*&
                    &(v_mag(:,2,:) - v_mag(:,1,:))/&
                    &(v_mag(:,2,:) + v_mag(:,1,:))
                call plot_HDF5(trim(base_name),trim(file_names(3)),&
                    &v_mag(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2)),&
                    &tot_dim=plot_dim(1:3),loc_offset=plot_offset(1:3),&
                    &X=X(:,tor_id(1):tor_id(2),:,1),&
                    &Y=Y(:,tor_id(1):tor_id(2),:,1),&
                    &Z=Z(:,tor_id(1):tor_id(2),:,1),&
                    &cont_plot=cont_plot,description=description(3))
            end if
            
            ! plot vector as well
            call plot_HDF5([trim(base_name)],trim(base_name)//'_vec',&
                &v_com(:,tor_id(1):tor_id(2),norm_id(1):norm_id(2),:,1),&
                &tot_dim=plot_dim,loc_offset=plot_offset,&
                &X=X(:,tor_id(1):tor_id(2),:,:),Y=Y(:,tor_id(1):tor_id(2),:,:),&
                &Z=Z(:,tor_id(1):tor_id(2),:,:),col=4,cont_plot=cont_plot,&
                &description='magnetic field vector')
        end if
        
        call lvl_ud(-1)
        
        ! clean up
        call norm_interp_data%dealloc()
        call grid_trim%dealloc()
    end function  calc_vec_comp
    
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
        integer :: ml_x                                                         ! location of maximum of x
        integer :: lim(2)                                                       ! limits in x
        integer :: lim_loc(2)                                                   ! limits in loc_x
        real(dp), allocatable :: x_fund(:)                                      ! x on fundamental interval 0..2pi
        real(dp), allocatable :: x_loc(:)                                       ! local x, extended
        real(dp), allocatable :: work(:)                                        ! work array
        real(dp), allocatable :: f_loc(:)                                       ! local f, extended
        real(dp), allocatable :: f_int(:)                                       ! interpolated f, later Fourier modes
        type(disc_type) :: trigon_interp_data                                   ! data for non-equidistant sampling fourier coefficients
        character(len=max_str_ln) :: plot_title(2)                              ! name of plot
        logical :: print_log = .false.                                          ! print log plot as well
        
        ! tests
        if (size(x).ne.size(f)) then
            ierr = 1
            CHCKERR('x and f must have same size')
        end if
        if (maxval(x)-minval(x).gt.2*pi) then
            ierr = 1
            CHCKERR('x is not a fundamental interval')
        end if
        
        ! create local x and  f that go from <0 to <2pi,  with enough points for
        ! full order precision
        allocate(x_fund(size(x)))
        x_fund = mod(x,2*pi)
        where(x_fund.lt.0._dp) x_fund = x_fund+2*pi
        allocate(x_loc(size(x)+2*interp_ord))
        allocate(f_loc(size(x)+2*interp_ord))
        ml_x = maxloc(x_fund,1)
        lim_loc = [1,interp_ord]
        lim = [ml_x-interp_ord+1,ml_x]
        x_loc(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))-2*pi
        f_loc(lim_loc(1):lim_loc(2)) = f(lim(1):lim(2))
        lim_loc = [interp_ord+1,interp_ord+size(x)-ml_x]
        lim = [ml_x+1,size(x)]
        x_loc(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))
        f_loc(lim_loc(1):lim_loc(2)) = f(lim(1):lim(2))
        lim_loc = [interp_ord+size(x)-ml_x+1,size(x)+interp_ord]
        lim = [1,ml_x]
        x_loc(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))
        f_loc(lim_loc(1):lim_loc(2)) = f(lim(1):lim(2))
        lim_loc(1) = lim_loc(2)+1
        do id = ml_x+1,min(size(x),ml_x+interp_ord)
            x_loc(lim_loc(1)) = x_fund(id)+2*pi
            f_loc(lim_loc(1)) = f(id)
            lim_loc(1) = lim_loc(1)+1
        end do
        do id = 1,ml_x+interp_ord-size(x)
            x_loc(lim_loc(1)) = x_fund(id)+2*pi
            f_loc(lim_loc(1)) = f(id)
            lim_loc(1) = lim_loc(1)+1
        end do
        
        ! set local f and interpolate
        n_x = size(x)
        allocate(f_int(n_x))
        ierr = setup_interp_data(x_loc,[((id-1._dp)/n_x*2*pi,id=1,n_x)],&
            &trigon_interp_data,interp_ord,is_trigon=.true.)
        ierr = apply_disc(f_loc,trigon_interp_data,f_int)
        CHCKERR('')
        
        ! clean up
        deallocate(x_fund)
        deallocate(x_loc)
        deallocate(f_loc)
        call trigon_interp_data%dealloc()
        
        !!do id = 1,n_x
            !!f_int(id) = -4._dp+&
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
        call dfftf(n_x,f_int,work)
        deallocate(work)
        
        ! rescale
        f_int(:) = f_int(:)*2/n_x
        f_int(1) = f_int(1)/2
        
        ! separate cos and sine
        if (allocated(f_F)) deallocate(f_F)
        allocate(f_F(m_F+1,2))
        f_F(1,1) = f_int(1)
        f_F(1,2) = 0._dp
        f_F(2:m_F+1,1) = f_int(2:2*m_F:2)
        f_F(2:m_F+1,2) = -f_int(3:2*m_F+1:2)                                    ! routine returns - sine factors
        
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
    
    ! finds smallest range that comprises a minimum and maximum value
    subroutine find_compr_range(r_F,lim_r,lim_id)
        ! input / output
        real(dp), intent(in) :: r_F(:)                                          ! all values of coordinate
        real(dp), intent(in) :: lim_r(2)                                        ! limiting range
        integer, intent(inout) :: lim_id(2)                                     ! limiting indices
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_r                                                          ! number of points in coordinate
        real(dp), parameter :: tol = 1.E-6                                      ! tolerance for grids
        
        ! set n_r
        n_r = size(r_F)
        
        lim_id = [1,n_r]                                                        ! initialize with full range
        if (r_F(1).lt.r_F(n_r)) then                                            ! ascending r_F
            do id = 1,n_r
                if (r_F(id).le.lim_r(1)-tol) lim_id(1) = id                     ! move lower limit up
                if (r_F(n_r-id+1).ge.lim_r(2)+tol) &
                    &lim_id(2) = n_r-id+1                                       ! move upper limit down
            end do
        else                                                                    ! descending r_F
            do id = 1,n_r
                if (r_F(id).ge.lim_r(2)+tol) lim_id(1) = id                     ! move lower limit up
                if (r_F(n_r-id+1).le.lim_r(1)-tol) &
                    &lim_id(2) = n_r-id+1                                       ! move upper limit down
            end do
        end if
    end subroutine find_compr_range
    
    ! Calculate arclength  angle, based  on calculating  the arclength  from the
    ! start of the grid, and then normalizing it from min(theta)..max(theta)
    ! By  default, the  Flux variables are  used, but this  can be  changed with
    ! "use_E".
    ! Note:  For field-aligned  grids the  projection is  taken on  the poloidal
    ! cross-section. For non-axisymmetric equilibria, this makes little sense.
    integer function calc_arc_angle(grid,eq_1,eq_2,arc,use_E) result(ierr)
        use eq_vars, only: eq_1_type, eq_2_type, &
            &R_0
        use num_vars, only: use_normalization, use_pol_flux_F
        use num_utilities, only: c, calc_int
        
        character(*), parameter :: rout_name = 'calc_arc_angle'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     ! flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     ! metric equilibrium variables
        real(dp), intent(inout), allocatable :: arc(:,:,:)                      ! arclength angle
        logical, intent(in), optional :: use_E                                  ! use E variables instead of F
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        logical :: use_E_loc                                                    ! local use_E
        
        ! initialize ierr
        ierr = 0
        
        ! set local use_E
        use_E_loc = .false.
        if (present(use_E)) use_E_loc = use_E
        
        ! calculate arclength: int dtheta |e_theta|
        allocate(arc(grid%n(1),grid%n(2),grid%loc_n_r))
        do kd = 1,grid%loc_n_r
            do jd = 1,grid%n(2)
                if (use_E_loc) then
                    ierr = calc_int(&
                        &eq_2%g_E(:,jd,kd,c([2,2],.true.),0,0,0)**(0.5),&
                        &grid%theta_E(:,jd,kd),arc(:,jd,kd))
                    CHCKERR('')
                    if (use_normalization) arc(:,jd,kd) = arc(:,jd,kd) * R_0    ! not really necessary as we take fractional arc length
                    arc(:,jd,kd) = grid%theta_E(1,jd,kd) + &
                        &arc(:,jd,kd) / arc(grid%n(1),jd,kd) * &
                        &(grid%theta_E(grid%n(1),jd,kd)-grid%theta_E(1,jd,kd))
                else
                    if (use_pol_flux_F) then
                        ! e_arc = e_theta - q e_alpha
                        ierr = calc_int((&
                            &eq_2%g_FD(:,jd,kd,c([3,3],.true.),0,0,0) - &
                            &eq_2%g_FD(:,jd,kd,c([1,3],.true.),0,0,0) * &
                            &2*eq_1%q_saf_FD(kd,0) + &
                            &eq_2%g_FD(:,jd,kd,c([1,1],.true.),0,0,0) * &
                            &eq_1%q_saf_FD(kd,0)**2)**(0.5),&
                            &grid%theta_F(:,jd,kd),arc(:,jd,kd))
                        CHCKERR('')
                    else
                        ! e_arc = -e_alpha
                        ierr = calc_int(&
                            &(-eq_2%g_FD(:,jd,kd,c([1,1],.true.),0,0,0))&
                            &**(0.5),grid%theta_F(:,jd,kd),arc(:,jd,kd))
                        CHCKERR('')
                    end if
                    if (use_normalization) arc(:,jd,kd) = arc(:,jd,kd) * R_0    ! not really necessary as we take fractional arc length
                    arc(:,jd,kd) = grid%theta_F(1,jd,kd) + &
                        &arc(:,jd,kd) / arc(grid%n(1),jd,kd) * &
                        &(grid%theta_F(grid%n(1),jd,kd)-grid%theta_F(1,jd,kd))
                end if
            end do
        end do
    end function calc_arc_angle
end module grid_utilities
