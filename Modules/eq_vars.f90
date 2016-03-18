!------------------------------------------------------------------------------!
!   Variables that have to do with equilibrium quantities and the grid used in !
!   the calculations:                                                          !
!       - The  equilibrium  variables  are  comprised  of  the  variables that !
!       result from the equilibrium  calculation, such as pressure, rotational !
!       transform, etc.                                                        !
!       - The flux variables are tabulated on a 1D grid.                       !
!       - The metric variables are tabulated on a 3D grid                      !
!           + In the normal coordinate,  they are tabulated in the equilibrium !
!           grid.                                                              !
!           + In the  angular coordinates, they are tabulated  in the solution !
!           grid (VMEC), or in the  equilibrium grid followed by an adaptation !
!           to the solution grid (HELENA).                                     !
!       - However, this  is not necessarily always the case,  as it depends on !
!       the grid angles being aligned with  the grid (e.g. Providing theta and !
!       zeta in  the equilibrium  grid so  that the  magnetic field  lines are !
!       followed.                                                              !
!   Note: In general in PB3D, there are two kinds of variables, differing from !
!   one another in the way in which they are tabulated:                        !
!       - variables tabulated on the full output grid of the equilibrium code  !
!       - variables tabulated in an internal grid of this code                 !
!   In many places  in the code a  range in the normal  coordinate is selected !
!   for each of the variables on different processes. This selection has to be !
!   done correctly  and things  can get  a little  bit complicated  if trimmed !
!   grids  are  used (grids  that  have  no  overlap between  processes).  See !
!   grid_ops.                                                                  !
!------------------------------------------------------------------------------!
module eq_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, mu_0_original
    use grid_vars, only: grid_type

    implicit none
    private
    public create_eq, dealloc_eq, &
        &eq_1_type, eq_2_type, &
        &R_0, pres_0, B_0, psi_0, rho_0, T_0, max_flux_E, max_flux_F, vac_perm
    
    ! global variables
    real(dp) :: R_0, pres_0, rho_0                                              ! independent normalization constants for nondimensionalization
    real(dp) :: B_0, psi_0, T_0                                                 ! derived normalization constants for nondimensionalization
    real(dp) :: vac_perm = mu_0_original                                        ! either usual mu_0 (default) or normalized
    real(dp) :: max_flux_E                                                      ! max. flux in Equilibrium coordinates
    real(dp) :: max_flux_F                                                      ! max. flux in Flux coordinates
    
    ! flux equilibrium type
    ! The arrays here are of the form
    !   - (r)                               for vars. without derivs.
    !   - (r,Dr)                            for vars. with derivs.
    type :: eq_1_type
        real(dp), allocatable :: pres_E(:,:)                                    ! pressure, and norm. deriv.
        real(dp), allocatable :: q_saf_E(:,:)                                   ! safety factor
        real(dp), allocatable :: rot_t_E(:,:)                                   ! rot. transform
        real(dp), pointer :: flux_p_E(:,:) => null()                            ! poloidal flux and norm. deriv.
        real(dp), pointer :: flux_t_E(:,:) => null()                            ! toroidal flux and norm. deriv.
        real(dp), allocatable :: pres_FD(:,:)                                   ! pressure, and norm. deriv.
        real(dp), allocatable :: q_saf_FD(:,:)                                  ! safety factor
        real(dp), allocatable :: rot_t_FD(:,:)                                  ! rot. transform
        real(dp), pointer :: flux_p_FD(:,:) => null()                           ! poloidal flux and norm. deriv.
        real(dp), pointer :: flux_t_FD(:,:) => null()                           ! toroidal flux and norm. deriv.
        real(dp), allocatable :: rho(:)                                         ! density
    end type
    
    ! metric equilibrium type
    ! The arrays here are of the form
    !   - (angle_1,angle_2,r)               for scalar vars. without derivs.
    !   - (angle_1,angle_2,r,D123)          for scalar vars. with derivs.
    !   - (angle_1,angle_2,r,6/9,D1,D2,D3)  for tensorial vars. with derivs.
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    ! The last index refers  to the derivatives in coordinate 1,  2 and 3, which
    ! refer to the coordinates described in [ADD REF]:
    !   - For E(quilibrium) coordinates, they are (r,theta,zeta)_E
    !   - For F(lux) coordinates, they are (alpha,psi,ang_par)_F, where
    !       + alpha = zeta - q theta and ang_par = theta for pol. flux,
    !       + alpha = -theta + iota zeta and ang_par = zeta for tor. flux.
    ! The  fourth  index for  tensorial  variables  correspond  to  the 9  or  6
    ! (symmetric) different values:
    !   (1 4 7)      (1    )
    !   (2 5 8)  or  (2 4  )
    !   (3 6 9)      (3 5 6)
    ! Note that Fortran only allows for 7 dimensions in arrays.
    type :: eq_2_type
        ! coordinate variables R, Z and L (VMEC)
        real(dp), allocatable :: R_E(:,:,:,:,:,:)                               ! R in E(quilibrium) coord
        real(dp), allocatable :: Z_E(:,:,:,:,:,:)                               ! Z in E(quilibrium) coords.
        real(dp), allocatable :: L_E(:,:,:,:,:,:)                               ! L(ambda) in E(quilibrium) coords.
        ! upper (h) and lower (g) metric factors
        real(dp), allocatable :: g_C(:,:,:,:,:,:,:)                             ! in the C(ylindrical) coord. system
        real(dp), allocatable :: g_E(:,:,:,:,:,:,:)                             ! in the E(quilibrium) coord. system
        real(dp), allocatable :: h_E(:,:,:,:,:,:,:)                             ! in the E(quilibrium) coord. system
        real(dp), allocatable :: g_F(:,:,:,:,:,:,:)                             ! in the F(lux) coord. system with derivs. in V(MEC) system
        real(dp), allocatable :: h_F(:,:,:,:,:,:,:)                             ! in the F(lux) coord. system with derivs. in V(MEC) system
        real(dp), pointer :: g_FD(:,:,:,:,:,:,:) => null()                      ! in the F(lux) coord. system with derivs in F(lux) system
        real(dp), pointer :: h_FD(:,:,:,:,:,:,:) => null()                      ! in the F(lux) coord. system with derivs in F(lux) system
        ! transformation matrices
        real(dp), allocatable :: T_VC(:,:,:,:,:,:,:)                            ! C(ylindrical) to V(MEC)
        real(dp), allocatable :: T_EF(:,:,:,:,:,:,:)                            ! E(quilibrium) to F(lux)
        real(dp), allocatable :: T_FE(:,:,:,:,:,:,:)                            ! F(lux) to E(quilibrium)
        ! determinants of transformation matrices
        real(dp), allocatable :: det_T_VC(:,:,:,:,:,:)                          ! determinant of T_VC
        real(dp), allocatable :: det_T_EF(:,:,:,:,:,:)                          ! determinant of T_EF
        real(dp), allocatable :: det_T_FE(:,:,:,:,:,:)                          ! determinant of T_FE
        ! Jacobians
        real(dp), allocatable :: jac_C(:,:,:,:,:,:)                             ! jacobian of C(ylindrical) coord. system
        real(dp), allocatable :: jac_E(:,:,:,:,:,:)                             ! jacobian of E(quilibrium) coord. system
        real(dp), allocatable :: jac_F(:,:,:,:,:,:)                             ! jacobian of F(lux) coord. system with derivs. in V(MEC) system
        real(dp), pointer :: jac_FD(:,:,:,:,:,:) => null()                      ! jacobian of F(lux) coord. system with derivs. in F(lux) system
        ! derived variables
        real(dp), allocatable :: S(:,:,:)                                       ! magnetic shear
        real(dp), allocatable :: kappa_n(:,:,:)                                 ! normal curvature
        real(dp), allocatable :: kappa_g(:,:,:)                                 ! geodesic curvature
        real(dp), allocatable :: sigma(:,:,:)                                   ! parallel current
    end type
    
    ! interfaces
    interface create_eq
        module procedure create_eq_1, create_eq_2
    end interface
    interface dealloc_eq
        module procedure dealloc_eq_1, dealloc_eq_2
    end interface

contains
    ! creates new equilibrium
    ! The normal and angular grid can be  in any coord. system, as only the grid
    ! sizes are used, not the coordinate values.
    ! Optionally, it can be chosen individually whether the E or F(D) quantities
    ! are allocated. The rationale behind this is that the E quantities are only
    ! used in  the pre-perturbation  phase. Note  that they  are not  written in
    ! print_output_eq either.
    ! Note:  The quantities  that  do not  have a  derivative  are considered  F
    ! quantities. Alternatively, all quantities that  have only one version, are
    ! considered F quantities, such as rho, kappa_n, ...
    subroutine create_eq_1(grid,eq,setup_E,setup_F)                             ! flux version
        use num_vars, only: max_deriv, eq_style
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(eq_1_type), intent(inout) :: eq                                    ! equilibrium to be created
        logical, intent(in), optional :: setup_E                                ! whether to set up E
        logical, intent(in), optional :: setup_F                                ! whether to set up F
        
        ! local variables
        integer :: loc_n_r, n                                                   ! local and total nr. of normal points
        logical :: setup_E_loc, setup_F_loc                                     ! local versions of setup_E and setup_F
        
        ! setup local setup_E and setup_F
        setup_E_loc = .true.
        if (present(setup_E)) setup_E_loc = setup_E
        setup_F_loc = .true.
        if (present(setup_F)) setup_F_loc = setup_F
        
        ! set local variables
        loc_n_r = grid%loc_n_r
        n = grid%n(3)
        
        if (setup_E_loc) then
            ! pres_E
            allocate(eq%pres_E(loc_n_r,0:max_deriv+1))
            
            ! q_saf_E
            allocate(eq%q_saf_E(loc_n_r,0:max_deriv+1))
            
            ! rot_t_E
            allocate(eq%rot_t_E(loc_n_r,0:max_deriv+1))
            
            ! initialize  variables that  are  specificic  to which  equilibrium
            ! style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! flux_p_E
                    allocate(eq%flux_p_E(loc_n_r,0:max_deriv+2))                ! Need extra order because used in transformation of flux q.
                    
                    ! flux_t_E
                    allocate(eq%flux_t_E(loc_n_r,0:max_deriv+2))                ! Need extra order because used in transformation of flux q.
                case (2)                                                        ! HELENA
                    ! flux_p_E
                    allocate(eq%flux_p_E(loc_n_r,0:max_deriv+1))
                    
                    ! flux_t_E
                    allocate(eq%flux_t_E(loc_n_r,0:max_deriv+1))
            end select
        end if
        
        if (setup_F_loc) then
            ! pres_FD
            allocate(eq%pres_FD(loc_n_r,0:max_deriv+1))
            
            ! flux_p_FD
            allocate(eq%flux_p_FD(loc_n_r,0:max_deriv+1))
            
            ! flux_t_FD
            allocate(eq%flux_t_FD(loc_n_r,0:max_deriv+1))
            
            ! q_saf_FD
            allocate(eq%q_saf_FD(loc_n_r,0:max_deriv+1))
            
            ! rot_t_FD
            allocate(eq%rot_t_FD(loc_n_r,0:max_deriv+1))
            
            ! rho
            allocate(eq%rho(grid%loc_n_r))
        end if
    end subroutine create_eq_1
    subroutine create_eq_2(grid,eq,setup_E,setup_F)                             ! metric version
        use num_vars, only: max_deriv, eq_style
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(eq_2_type), intent(inout) :: eq                                    ! equilibrium to be created
        logical, intent(in), optional :: setup_E                                ! whether to set up E
        logical, intent(in), optional :: setup_F                                ! whether to set up F
        
        ! local variables
        logical :: setup_E_loc, setup_F_loc                                     ! local versions of setup_E and setup_F
        integer :: dims(3)                                                      ! dimensions
        
        ! set up local variables
        dims = [grid%n(1),grid%n(2),grid%loc_n_r]
        
        ! setup local setup_E and setup_F
        setup_E_loc = .true.
        if (present(setup_E)) setup_E_loc = setup_E
        setup_F_loc = .true.
        if (present(setup_F)) setup_F_loc = setup_F
        
        if (setup_E_loc) then
            ! initialize variables that are used for all equilibrium styles
            ! g_E
            allocate(eq%g_E(dims(1),dims(2),dims(3),6,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! h_E
            allocate(eq%h_E(dims(1),dims(2),dims(3),6,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! g_F
            allocate(eq%g_F(dims(1),dims(2),dims(3),6,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! h_F
            allocate(eq%h_F(dims(1),dims(2),dims(3),6,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! jac_E
            allocate(eq%jac_E(dims(1),dims(2),dims(3),&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! jac_F
            allocate(eq%jac_F(dims(1),dims(2),dims(3),&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! T_EF
            allocate(eq%T_EF(dims(1),dims(2),dims(3),9,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! T_FE
            allocate(eq%T_FE(dims(1),dims(2),dims(3),9,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! det_T_EF
            allocate(eq%det_T_EF(dims(1),dims(2),dims(3),&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! det_T_FE
            allocate(eq%det_T_FE(dims(1),dims(2),dims(3),&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! initialize  variables that  are  specificic  to which  equilibrium
            ! style is being used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! g_C
                    allocate(eq%g_C(dims(1),dims(2),dims(3),6,&
                        &0:max_deriv,0:max_deriv,0:max_deriv))
                    
                    ! T_VC
                    allocate(eq%T_VC(dims(1),dims(2),dims(3),9,&
                        &0:max_deriv,0:max_deriv,0:max_deriv))
                    
                    ! det_T_VC
                    allocate(eq%det_T_VC(dims(1),dims(2),dims(3),&
                        &0:max_deriv,0:max_deriv,0:max_deriv))
                    
                    ! jac_C
                    allocate(eq%jac_C(dims(1),dims(2),dims(3),&
                        &0:max_deriv,0:max_deriv,0:max_deriv))
                    
                    ! R
                    allocate(eq%R_E(dims(1),dims(2),dims(3),&
                        &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                    
                    ! Z
                    allocate(eq%Z_E(dims(1),dims(2),dims(3),&
                        &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                    
                    ! lambda
                    allocate(eq%L_E(dims(1),dims(2),dims(3),&
                        &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                case (2)                                                        ! HELENA
                    ! do nothing
            end select
        end if
        
        if (setup_F_loc) then
            ! initialize variables that are used for all equilibrium styles
            ! g_FD
            allocate(eq%g_FD(dims(1),dims(2),dims(3),6,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! h_FD
            allocate(eq%h_FD(dims(1),dims(2),dims(3),6,&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! jac_FD
            allocate(eq%jac_FD(dims(1),dims(2),dims(3),&
                &0:max_deriv,0:max_deriv,0:max_deriv))
            
            ! magnetic shear
            allocate(eq%S(dims(1),dims(2),dims(3)))
            
            ! normal curvature
            allocate(eq%kappa_n(dims(1),dims(2),dims(3)))
            
            ! geodesic curvature
            allocate(eq%kappa_g(dims(1),dims(2),dims(3)))
            
            ! parallel current
            allocate(eq%sigma(dims(1),dims(2),dims(3)))
        end if
    end subroutine create_eq_2
    
    ! deallocates equilibrium quantities
    subroutine dealloc_eq_1(eq)                                                 ! flux version
        ! input / output
        type(eq_1_type), intent(inout) :: eq                                    ! equilibrium to be deallocated
        
        ! deallocate and nullify allocated pointers
        if (associated(eq%flux_p_E)) then
            deallocate(eq%flux_p_E)
            nullify(eq%flux_p_E)
        end if
        if (associated(eq%flux_t_E)) then
            deallocate(eq%flux_t_E)
            nullify(eq%flux_t_E)
        end if
        if (associated(eq%flux_p_FD)) then
            deallocate(eq%flux_p_FD)
            nullify(eq%flux_p_FD)
        end if
        if (associated(eq%flux_t_FD)) then
            deallocate(eq%flux_t_FD)
            nullify(eq%flux_t_FD)
        end if
        
        ! deallocate allocatable variables
        call dealloc_eq_1_final(eq)
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_eq_1_final(eq)
            ! input / output
            type(eq_1_type), intent(out) :: eq                                  ! equilibrium to be deallocated
        end subroutine dealloc_eq_1_final
    end subroutine dealloc_eq_1
    subroutine dealloc_eq_2(eq)                                                 ! metric version
        ! input / output
        type(eq_2_type), intent(inout) :: eq                                    ! equilibrium to be deallocated
        
        ! deallocate and nullify allocated pointers
        if (associated(eq%g_FD)) then
            deallocate(eq%g_FD)
            nullify(eq%g_FD)
        end if
        if (associated(eq%h_FD)) then
            deallocate(eq%h_FD)
            nullify(eq%h_FD)
        end if
        if (associated(eq%jac_FD)) then
            deallocate(eq%jac_FD)
            nullify(eq%jac_FD)
        end if
        
        ! deallocate allocatable variables
        call dealloc_eq_2_final(eq)
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_eq_2_final(eq)
            ! input / output
            type(eq_2_type), intent(out) :: eq                                  ! equilibrium to be deallocated
        end subroutine dealloc_eq_2_final
    end subroutine dealloc_eq_2
end module eq_vars
