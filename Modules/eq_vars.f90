!------------------------------------------------------------------------------!
!> Variables that have to do with equilibrium quantities and the grid used in
!! the calculations:
!!  - The equilibrium variables are comprised  of the variables that result from
!!  the equilibrium calculation, such as pressure, rotational transform, etc.
!!  - The flux variables are tabulated on a 1D grid.
!!  - The metric variables are tabulated on a 3D grid
!!      + In the normal coordinate, they are tabulated in the equilibrium grid.
!!      + In  the angular coordinates, they  are tabulated in the  solution grid
!!      (VMEC), or  in the  equilibrium grid  followed by  an adaptation  to the
!!      solution grid (HELENA).
!!  - However, this  is not necessarily always  the case, as it  depends on the 
!!  grid angles being  aligned with the grid (e.g. Providing  theta and zeta in 
!!  the equilibrium grid so that the magnetic field lines are followed.         
!!
!! \note In general in PB3D, there are two kinds of variables, differing from 
!! one another in the way in which they are tabulated:
!!  - variables tabulated on the full output grid of the equilibrium code
!!  - variables tabulated in an internal grid of this code
!! \note In many places in the code a range in the normal coordinate is selected
!! for each  of the variables on  different processes. This selection  has to be
!! done correctly and  things can get a little bit  complicated if trimmed grids
!! are used (grids that have no overlap between processes).
!! 
!! \see grid_ops()
!------------------------------------------------------------------------------!
module eq_vars
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, pi, max_str_ln, mu_0_original, weight_dp
    use grid_vars, only: grid_type

    implicit none
    private
    public R_0, pres_0, B_0, psi_0, rho_0, T_0, max_flux_E, max_flux_F, vac_perm
#if ldebug
    public n_alloc_eq_1s, n_alloc_eq_2s
#endif
    
    ! global variables
    real(dp) :: R_0                                                             !< independent normalization constant for nondimensionalization
    real(dp) :: pres_0                                                          !< independent normalization constant for nondimensionalization
    real(dp) :: rho_0                                                           !< independent normalization constant for nondimensionalization
    real(dp) :: B_0                                                             !< derived normalization constant for nondimensionalization
    real(dp) :: psi_0                                                           !< derived normalization constant for nondimensionalization
    real(dp) :: T_0                                                             !< derived normalization constant for nondimensionalization
    real(dp) :: vac_perm = mu_0_original                                        !< either usual mu_0 (default) or normalized
    real(dp) :: max_flux_E                                                      !< max. flux in Equilibrium coordinates, set in calc_norm_range_PB3D_in
    real(dp) :: max_flux_F                                                      !< max. flux in Flux coordinates, set in calc_norm_range_PB3D_in
#if ldebug
    !> \ldebug
    integer :: n_alloc_eq_1s                                                    !< nr. of allocated \c eq_1 variables
    !> \ldebug
    integer :: n_alloc_eq_2s                                                    !< nr. of allocated \c eq_2 variables
#endif
    
    !> flux equilibrium type
    !!
    !! The arrays here are of the form
    !!  - <tt>(r)</tt>                          for variables without derivs.
    !!  - <tt>(r,Dr)</tt>                       for variables with derivs.
    type, public :: eq_1_type
        real(dp), allocatable :: pres_E(:,:)                                    !< pressure, and norm. deriv.
        real(dp), allocatable :: q_saf_E(:,:)                                   !< safety factor
        real(dp), allocatable :: rot_t_E(:,:)                                   !< rot. transform
        real(dp), allocatable :: flux_p_E(:,:)                                  !< poloidal flux and norm. deriv.
        real(dp), allocatable :: flux_t_E(:,:)                                  !< toroidal flux and norm. deriv.
        real(dp), allocatable :: pres_FD(:,:)                                   !< pressure, and norm. deriv.
        real(dp), allocatable :: q_saf_FD(:,:)                                  !< safety factor
        real(dp), allocatable :: rot_t_FD(:,:)                                  !< rot. transform
        real(dp), allocatable :: flux_p_FD(:,:)                                 !< poloidal flux and norm. deriv.
        real(dp), allocatable :: flux_t_FD(:,:)                                 !< toroidal flux and norm. deriv.
        real(dp), allocatable :: rho(:)                                         !< density
#if ldebug
        real(dp) :: estim_mem_usage(2)                                          !< expected memory usage \ldebug
#endif
    contains
        !> initialize
        procedure :: init => init_eq_1
        !> copy
        procedure :: copy => copy_eq_1
        !> deallocate
        procedure :: dealloc => dealloc_eq_1
    end type
    
    !> metric equilibrium type
    !! 
    !! The arrays here are of the form
    !!  - <tt>(angle_1,angle_2,r)</tt>              for scalars without derivs.
    !!  - <tt>(angle_1,angle_2,r,D123)</tt>         for scalars with derivs.
    !!  - <tt>(angle_1,angle_2,r,6/9,D1,D2,D3)</tt> for tensors with derivs.
    !! \n  where it  is  refered to  the  discussion  of the  grid  type for  an
    !! explanation of the angles \c angle_1 and \c angle_2.
    !!
    !! The last index refers to the derivatives  in coordinate 1, 2 and 3, which
    !! refer to the coordinates described in \cite Weyens3D.
    !!  - For E(quilibrium) coordinates, they are
    !!      \f$\left(r,\theta,\zeta\right)_\text{E}\f$.
    !!  - For F(lux) coordinates, they are
    !!      \f$\left(\alpha,\psi,\gamma\right)_\text{F}\f$, where
    !!      + \f$\alpha = \zeta - q \theta\f$ and \f$\gamma = \theta\f$
    !!          for pol. flux,
    !!      + \f$\alpha = -\theta + \iota \zeta\f$ and \f$\gamma = \zeta\f$
    !!         Nfor tor. flux.
    !!
    !! The  fourth index  for  tensorial  variables correspond  to  the  9 or  6
    !! (symmetric) different values:
    !!  \f[ \left( \begin{array}{ccc}1&4&7\\2&5&8\\3&6&9\end{array} \right) \
    !!  \text{or} \quad
    !!  \left( \begin{array}{ccc}1& & \\2&4& \\3&5&6\end{array} \right) \f]
    !!
    !! \note Fortran only allows for 7 dimensions in arrays.
    type, public :: eq_2_type
        ! coordinate variables R, Z and L (VMEC)
        real(dp), allocatable :: R_E(:,:,:,:,:,:)                               !< R in E(quilibrium) coord
        real(dp), allocatable :: Z_E(:,:,:,:,:,:)                               !< Z in E(quilibrium) coords.
        real(dp), allocatable :: L_E(:,:,:,:,:,:)                               !< L(ambda) in E(quilibrium) coords.
        ! upper (h) and lower (g) metric factors
        real(dp), allocatable :: g_C(:,:,:,:,:,:,:)                             !< in the C(ylindrical) coord. system
        real(dp), allocatable :: g_E(:,:,:,:,:,:,:)                             !< in the E(quilibrium) coord. system
        real(dp), allocatable :: h_E(:,:,:,:,:,:,:)                             !< in the E(quilibrium) coord. system
        real(dp), allocatable :: g_F(:,:,:,:,:,:,:)                             !< in the F(lux) coord. system with derivs. in V(MEC) system
        real(dp), allocatable :: h_F(:,:,:,:,:,:,:)                             !< in the F(lux) coord. system with derivs. in V(MEC) system
        real(dp), allocatable :: g_FD(:,:,:,:,:,:,:)                            !< in the F(lux) coord. system with derivs in F(lux) system
        real(dp), allocatable :: h_FD(:,:,:,:,:,:,:)                            !< in the F(lux) coord. system with derivs in F(lux) system
        ! transformation matrices
        real(dp), allocatable :: T_VC(:,:,:,:,:,:,:)                            !< C(ylindrical) to V(MEC)
        real(dp), allocatable :: T_EF(:,:,:,:,:,:,:)                            !< E(quilibrium) to F(lux)
        real(dp), allocatable :: T_FE(:,:,:,:,:,:,:)                            !< F(lux) to E(quilibrium)
        ! determinants of transformation matrices
        real(dp), allocatable :: det_T_EF(:,:,:,:,:,:)                          !< determinant of T_EF
        real(dp), allocatable :: det_T_FE(:,:,:,:,:,:)                          !< determinant of T_FE
        ! Jacobians
        real(dp), allocatable :: jac_E(:,:,:,:,:,:)                             !< jacobian of E(quilibrium) coord. system
        real(dp), allocatable :: jac_F(:,:,:,:,:,:)                             !< jacobian of F(lux) coord. system with derivs. in V(MEC) system
        real(dp), allocatable :: jac_FD(:,:,:,:,:,:)                            !< jacobian of F(lux) coord. system with derivs. in F(lux) system
        ! derived variables
        real(dp), allocatable :: S(:,:,:)                                       !< magnetic shear
        real(dp), allocatable :: kappa_n(:,:,:)                                 !< normal curvature
        real(dp), allocatable :: kappa_g(:,:,:)                                 !< geodesic curvature
        real(dp), allocatable :: sigma(:,:,:)                                   !< parallel current
#if ldebug
        real(dp) :: estim_mem_usage(2)                                          !< expected memory usage \ldebug
#endif
    contains
        !> initialize
        procedure :: init => init_eq_2
        !> copy
        procedure :: copy => copy_eq_2
        !> deallocate
        procedure :: dealloc => dealloc_eq_2
    end type

contains
    !> \public Initializes new flux equilibrium.
    !!
    !! The normal and angular grid can be in any coord. system, as only the grid
    !! sizes are used, not the coordinate values.
    !!
    !! Optionally,  it  can  be  chosen  individually  whether  the  E  or  F(D)
    !! quantities are allocated. D means that  the derivatives are in F as well.
    !! The rationale behind this  is that the E quantities are  only used in the
    !! pre-perturbation phase.
    !! \note
    !!  -# The E quantities are not written in eq_ops.print_output_eq().
    !!  -#  The quantities  that  do  not have  a  derivative  are considered  F
    !!  quantities. Alternatively,  all quantities  that have only  one version,
    !!  are considered F quantities, such as \c rho, \c kappa_n, ...
    subroutine init_eq_1(eq,grid,setup_E,setup_F)
        use num_vars, only: max_deriv, eq_style
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        
        ! input / output
        class(eq_1_type), intent(inout) :: eq                                   !< equilibrium to be initialized
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid
        logical, intent(in), optional :: setup_E                                !< whether to set up E
        logical, intent(in), optional :: setup_F                                !< whether to set up F
        
        ! local variables
        integer :: loc_n_r, n                                                   ! local and total nr. of normal points
        logical :: setup_E_loc, setup_F_loc                                     ! local versions of setup_E and setup_F
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) eq%estim_mem_usage = 0._dp
#endif
        
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
            
#if ldebug
            ! set estimated memory usage
            if (print_mem_usage) eq%estim_mem_usage(1) = &
                &eq%estim_mem_usage(1) + loc_n_r*(max_deriv+2)*3
#endif
            
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
                    
#if ldebug
                    ! set estimated memory usage
                    if (print_mem_usage) eq%estim_mem_usage(1) = &
                        &eq%estim_mem_usage(1) + loc_n_r*(max_deriv+3)*2
#endif
                case (2)                                                        ! HELENA
                    ! flux_p_E
                    allocate(eq%flux_p_E(loc_n_r,0:max_deriv+1))
                    
                    ! flux_t_E
                    allocate(eq%flux_t_E(loc_n_r,0:max_deriv+1))
                    
#if ldebug
                    ! set estimated memory usage
                    if (print_mem_usage) eq%estim_mem_usage(1) = &
                        &eq%estim_mem_usage(1) + loc_n_r*(max_deriv+2)*2
#endif
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
            
#if ldebug
            ! set estimated memory usage
            if (print_mem_usage) eq%estim_mem_usage(2) = &
                &eq%estim_mem_usage(2) + loc_n_r*((max_deriv+2)*5+1)
#endif
        end if
        
#if ldebug
        ! increment n_alloc_eq_1s
        n_alloc_eq_1s = n_alloc_eq_1s + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of eq_1: '//&
            &trim(r2strt(eq%estim_mem_usage(1)*0.008))//' kB for E and '//&
            &trim(r2strt(eq%estim_mem_usage(2)*0.008))//' kB for E]',&
            &alert=.true.)
#endif
    end subroutine init_eq_1
    
    !> \public Initializes new metric equilibrium.
    !!
    !! \see For explanation see init_eq_1().
    subroutine init_eq_2(eq,grid,setup_E,setup_F)                               ! metric version
        use num_vars, only: max_deriv, eq_style
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        
        ! input / output
        class(eq_2_type), intent(inout) :: eq                                   !< equilibrium to be initialized
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid
        logical, intent(in), optional :: setup_E                                !< whether to set up E
        logical, intent(in), optional :: setup_F                                !< whether to set up F
        
        ! local variables
        logical :: setup_E_loc, setup_F_loc                                     ! local versions of setup_E and setup_F
        integer :: dims(3)                                                      ! dimensions
        
        ! set up local variables
        dims = [grid%n(1),grid%n(2),grid%loc_n_r]
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) eq%estim_mem_usage = 0._dp
#endif
        
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
            
#if ldebug
            ! set estimated memory usage
            if (print_mem_usage) eq%estim_mem_usage(1) = &
                &eq%estim_mem_usage(1) + product(dims)*&
                &((max_deriv+1)**3*(4*6+2*9+4))
#endif
            
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
                    
                    ! R
                    allocate(eq%R_E(dims(1),dims(2),dims(3),&
                        &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                    
                    ! Z
                    allocate(eq%Z_E(dims(1),dims(2),dims(3),&
                        &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                    
                    ! lambda
                    allocate(eq%L_E(dims(1),dims(2),dims(3),&
                        &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                    
#if ldebug
                    ! set estimated memory usage
                    if (print_mem_usage) eq%estim_mem_usage(1) = &
                        &eq%estim_mem_usage(1) + product(dims)*(&
                        &(max_deriv+1)**3*(1*6+1*9+2) + &
                        &(max_deriv+2)**3*(3))
#endif
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
            
#if ldebug
            ! set estimated memory usage
            if (print_mem_usage) eq%estim_mem_usage(2) = &
                &eq%estim_mem_usage(2) + product(dims)*(&
                &(max_deriv+1)**3*(2*6+1) + &
                &4)
#endif
        end if
        
#if ldebug
        ! increment n_alloc_eq_2s
        n_alloc_eq_2s = n_alloc_eq_2s + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of eq_2: '//&
            &trim(r2strt(eq%estim_mem_usage(1)*0.008))//' kB for E and '//&
            &trim(r2strt(eq%estim_mem_usage(2)*0.008))//' kB for E]',&
            &alert=.true.)
#endif
    end subroutine init_eq_2
    
    !> \public Deep copy of flux equilibrium variables.
    subroutine copy_eq_1(eq_i,grid_i,eq_o)
        use grid_vars, only: grid_type
        
        ! input / output
        class(eq_1_type), intent(in) :: eq_i                                    !< eq_1 to be copied
        type(grid_type), intent(in) :: grid_i                                   !< grid of eq_i
        type(eq_1_type), intent(inout) :: eq_o                                  !< copied eq_1
        
        ! local variables
        logical :: setup_E
        logical :: setup_F
        
        setup_E = allocated(eq_i%pres_E)
        setup_F = allocated(eq_i%pres_FD)
        
        call eq_o%init(grid_i,setup_E=setup_E,setup_F=setup_F)
        
        if (setup_E) then
            if (allocated(eq_i%pres_E)) eq_o%pres_E = eq_i%pres_E
            if (allocated(eq_i%q_saf_E)) eq_o%q_saf_E = eq_i%q_saf_E
            if (allocated(eq_i%rot_t_E)) eq_o%rot_t_E = eq_i%rot_t_E
            if (allocated(eq_i%flux_p_E)) eq_o%flux_p_E = eq_i%flux_p_E
            if (allocated(eq_i%flux_t_E)) eq_o%flux_t_E = eq_i%flux_t_E
        end if
        
        if (setup_F) then
            if (allocated(eq_i%pres_FD)) eq_o%pres_FD = eq_i%pres_FD
            if (allocated(eq_i%flux_p_FD)) eq_o%flux_p_FD = eq_i%flux_p_FD
            if (allocated(eq_i%flux_t_FD)) eq_o%flux_t_FD = eq_i%flux_t_FD
            if (allocated(eq_i%q_saf_FD)) eq_o%q_saf_FD = eq_i%q_saf_FD
            if (allocated(eq_i%rot_t_FD)) eq_o%rot_t_FD = eq_i%rot_t_FD
            if (allocated(eq_i%rho)) eq_o%rho = eq_i%rho
        end if
    end subroutine copy_eq_1
    
    !> \public Deep copy of metric equilibrium variables.
    !!
    !! \return ierr
    subroutine copy_eq_2(eq_i,grid_i,eq_o)
        use num_vars, only: eq_style
        use grid_vars, only: grid_type
        
        ! input / output
        class(eq_2_type), intent(in) :: eq_i                                    !< eq_2 to be copied
        type(grid_type), intent(in) :: grid_i                                   !< grid of eq_i
        type(eq_2_type), intent(inout) :: eq_o                                  !< copied eq_1
        
        ! local variables
        logical :: setup_E
        logical :: setup_F
        
        setup_E = allocated(eq_i%jac_E)
        setup_F = allocated(eq_i%jac_FD)
        
        call eq_o%init(grid_i,setup_E=setup_E,setup_F=setup_F)
        
        if (setup_E) then
            if (allocated(eq_i%g_E)) eq_o%g_E = eq_i%g_E
            if (allocated(eq_i%h_E)) eq_o%h_E = eq_i%h_E
            if (allocated(eq_i%g_F)) eq_o%g_F = eq_i%g_F
            if (allocated(eq_i%h_F)) eq_o%h_F = eq_i%h_F
            if (allocated(eq_i%jac_E)) eq_o%jac_E = eq_i%jac_E
            if (allocated(eq_i%jac_F)) eq_o%jac_F = eq_i%jac_F
            if (allocated(eq_i%T_EF)) eq_o%T_EF = eq_i%T_EF
            if (allocated(eq_i%T_FE)) eq_o%T_FE = eq_i%T_FE
            if (allocated(eq_i%det_T_EF)) eq_o%det_T_EF = eq_i%det_T_EF
            if (allocated(eq_i%det_T_FE)) eq_o%det_T_FE = eq_i%det_T_FE
            
            select case (eq_style)
                case (1)                                                        ! VMEC
                    if (allocated(eq_i%g_C)) eq_o%g_C = eq_i%g_C
                    if (allocated(eq_i%T_VC)) eq_o%T_VC = eq_i%T_VC
                    if (allocated(eq_i%R_E)) eq_o%R_E = eq_i%R_E
                    if (allocated(eq_i%Z_E)) eq_o%Z_E = eq_i%Z_E
                    if (allocated(eq_i%l_E)) eq_o%l_E = eq_i%l_E
                case (2)                                                        ! HELENA
                    ! do nothing
            end select
        end if
        
        if (setup_F) then
            if (allocated(eq_i%g_FD)) eq_o%g_FD = eq_i%g_FD
            if (allocated(eq_i%h_FD)) eq_o%h_FD = eq_i%h_FD
            if (allocated(eq_i%jac_FD)) eq_o%jac_FD = eq_i%jac_FD
            if (allocated(eq_i%S)) eq_o%S = eq_i%S
            if (allocated(eq_i%kappa_n)) eq_o%kappa_n = eq_i%kappa_n
            if (allocated(eq_i%kappa_g)) eq_o%kappa_g = eq_i%kappa_g
            if (allocated(eq_i%sigma)) eq_o%sigma = eq_i%sigma
        end if
    end subroutine copy_eq_2
    
    !> \public Deallocates flux equilibrium quantities.
    subroutine dealloc_eq_1(eq)
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        
        ! input / output
        class(eq_1_type), intent(inout) :: eq                                   !< equilibrium to be deallocated
        
#if ldebug
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = sum(eq%estim_mem_usage)
        end if
#endif
        
        ! deallocate allocatable variables
        call dealloc_eq_1_final(eq)
        
#if ldebug
        ! decrement n_alloc_eq_1s
        n_alloc_eq_1s = n_alloc_eq_1s - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating eq_1 ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        !> \private
        subroutine dealloc_eq_1_final(eq)
            ! input / output
            type(eq_1_type), intent(out) :: eq                                  ! equilibrium to be deallocated
        end subroutine dealloc_eq_1_final
    end subroutine dealloc_eq_1

    !> \public Deallocates metric equilibrium quantities.
    subroutine dealloc_eq_2(eq)
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        ! input / output
        class(eq_2_type), intent(inout) :: eq                                   !< equilibrium to be deallocated
        
#if ldebug
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = sum(eq%estim_mem_usage)
        end if
#endif
        
        ! deallocate allocatable variables
        call dealloc_eq_2_final(eq)
        
#if ldebug
        ! decrement n_alloc_eq_2s
        n_alloc_eq_2s = n_alloc_eq_2s - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating eq_2 ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        !> \private
        subroutine dealloc_eq_2_final(eq)
            ! input / output
            type(eq_2_type), intent(out) :: eq                                  ! equilibrium to be deallocated
        end subroutine dealloc_eq_2_final
    end subroutine dealloc_eq_2
end module eq_vars
