!------------------------------------------------------------------------------!
!   Variables that have to do with equilibrium quantities and the grid used in !
!   the calculations:                                                          !
!       -  The  equilibrium  variables  are comprised  of  the  variable  that !
!       result from the equilibrium  calculation, such as pressure, rotational !
!       transform, etc.                                                        !
!       - These variables  are tabulated on a 3D grid,  following the magnetic !
!       field:                                                                 !
!           + In the normal coordinate,  they are tabulated in the equilibrium !
!           grid.                                                              !
!           +  In  the   angular  coordinates,  they  are   tabulated  in  the !
!           perturbation grid (VMEC),  or in the equilibrium  grid followed by !
!           an adaptation to the perturbation grid (HELENA).                   !
!       - However, this  is not necessarily always the case,  as it depends on !
!       the grid angles being aligned with  the grid (e.g. Providing theta and !
!       zeta in  the equilibrium  grid so  that the  magnetic field  lines are !
!       followed.                                                              !
!   Note: In general in PB3D, there are two kinds of variables, differing from !
!   one another in the way in which they are tabulated:                        !
!       - variables tabulated on the full output grid of the equilibrium code  !
!       - variables tabulated in an internal grid of this code                 !
!   In many places  in the code a  range in the normal  coordinate is selected !
!   for each of  the variables on different processors. This  selection has to !
!   be done correctly  and things can get a little  bit complicated if trimmed !
!   grids are used (grids that have no overlap between processes):             !
!       -  Variables  on  full  output  grid  need to  keep  in  mind  that  a !
!       grid  trimmed  internal  can  start   at  a  position  different  from !
!       the  grid  starting position  of  the  full equilibrium  output  grid. !
!       Therefore,  the   correct  way  to   indicate  the  normal   range  of !
!       variables  tabulated in  the  full equilibrium  output  grid would  be !
!       [grid%i_min:grid%i_min+grid_trim%grp_n_r], where grid is the untrimmed !
!       and grid_trim the trimmed internal grid.                               !
!       -  Internal  grid by  default  are  trimmed  keeping the  lower  range !
!       unchanged,  so   the  usage   of  [1:grid_trim%grp_n_r]   is  allowed. !
!       However,  using shift_grid  shifting the  grid  by an  amount a,  both !
!       a+[1:grid_trim%grp_n_r]  or  [grid%i_min:grid%i_min+grid_trim%grp_n_r] !
!       can be used.                                                           !
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
        &eq_type, &
        &R_0, pres_0, B_0, psi_0, rho_0, T_0, max_flux_p_E, max_flux_p_F, &
        &max_flux_t_E, max_flux_t_F, vac_perm
    
    ! global variables
    real(dp) :: R_0, pres_0, rho_0                                              ! independent normalization constants for nondimensionalization
    real(dp) :: B_0, psi_0, T_0                                                 ! derived normalization constants for nondimensionalization
    real(dp) :: vac_perm = mu_0_original                                        ! either usual mu_0 (default) or normalized
    real(dp) :: max_flux_p_E, max_flux_p_F                                      ! max. pol. flux in Equilibrium and Flux coordinates
    real(dp) :: max_flux_t_E, max_flux_t_F                                      ! max. tor. flux in Equilibrium and Flux coordinates
    
    ! equilibrium type
    ! The arrays here are of the form (except for rho):
    !   - (r,Dr)                        for flux quantities
    !   - (angle_1,angle_2,r,D1,D2,D3)  for normal quantities
    ! where it is refered to the discussion  of the grid type for an explanation
    ! of the angles angle_1 and angle_2.
    ! The last three indices refer to the  derivatives in coordinate 1, 2 and 3,
    ! which refer to the coordinates described in [ADD REF]: 
    !   - For E(quilibrium) coordinates, they are (r,theta,zeta)_E
    !   - For F(lux) coordinates, they are (alpha,psi,ang_par)_F, where
    !       + alpha = zeta - q theta and ang_par = theta for pol. flux,
    !       + alpha = -theta + iota zeta and ang_par = zeta for tor. flux.
    ! This order  of variables is  important when setting up  the transformation
    ! matrices between the E and F coordinate systems.
    type :: eq_type
        real(dp), allocatable :: R_E(:,:,:,:,:,:)                               ! R in E(quilibrium) coord
        real(dp), allocatable :: Z_E(:,:,:,:,:,:)                               ! Z in E(quilibrium) coords.
        real(dp), allocatable :: L_E(:,:,:,:,:,:)                               ! L(ambda) in E(quilibrium) coords.
        real(dp), allocatable :: pres_E(:,:)                                    ! pressure, and norm. Deriv. in E(equilibrium) coords.
        real(dp), allocatable :: q_saf_E(:,:)                                   ! safety factor in E(equilibrium) coordinates
        real(dp), allocatable :: rot_t_E(:,:)                                   ! rot. transform in E(equilibrium) coordinates
        real(dp), pointer :: flux_p_E(:,:) => null()                            ! poloidal flux and norm. Deriv. in E(equilibrium) coords.
        real(dp), pointer :: flux_t_E(:,:) => null()                            ! toroidal flux and norm. Deriv. in E(equilibrium) coords.
        real(dp), allocatable :: pres_FD(:,:)                                   ! pressure, and norm. Deriv. with values and Derivs. in F(lux) coords.
        real(dp), allocatable :: q_saf_FD(:,:)                                  ! safety factor, Deriv. in F(lux) coords.
        real(dp), allocatable :: rot_t_FD(:,:)                                  ! rot. transform, Deriv. in F(lux) coords.
        real(dp), pointer :: flux_p_FD(:,:) => null()                           ! poloidal flux, and norm. Deriv. with values and Derivs. in F(lux) coords.
        real(dp), pointer :: flux_t_FD(:,:) => null()                           ! toroidal flux, and norm. Deriv. with values and Derivs. in F(lux) coords.
        real(dp), allocatable :: rho(:)                                         ! density (in all coord. systems)
        real(dp), allocatable :: S(:,:,:)                                       ! magnetic shear
        real(dp), allocatable :: kappa_n(:,:,:)                                 ! normal curvature
        real(dp), allocatable :: kappa_g(:,:,:)                                 ! geodesic curvature
        real(dp), allocatable :: sigma(:,:,:)                                   ! parallel current
    end type

contains
    ! creates new equilibrium
    ! The normal and angular grid can be  in any coord. system, as only the grid
    ! sizes are used, not the coordinate values.
    ! Note: intent(out) automatically deallocates the variable
    integer function create_eq(grid,eq) result(ierr)
        use num_vars, only: max_deriv, eq_style
        
        character(*), parameter :: rout_name = 'create_eq'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        type(eq_type), intent(out) :: eq                                        ! equilibrium to be created
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: grp_n_r, n                                                   ! group and total nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        
        ! initialize ierr
        ierr = 0
        
        ! set local variables
        grp_n_r = grid%grp_n_r
        n_par = grid%n(1)
        n_geo = grid%n(2)
        n = grid%n(3)
        
        ! pres_FD
        allocate(eq%pres_FD(grp_n_r,0:max_deriv))
        
        ! flux_p_FD
        allocate(eq%flux_p_FD(grp_n_r,0:max_deriv))
        
        ! flux_t_FD
        allocate(eq%flux_t_FD(grp_n_r,0:max_deriv))
        
        ! q_saf_FD
        allocate(eq%q_saf_FD(grp_n_r,0:max_deriv))
        
        ! rot_t_FD
        allocate(eq%rot_t_FD(grp_n_r,0:max_deriv))
        
        ! pres_E
        allocate(eq%pres_E(grp_n_r,0:max_deriv+1))
        
        ! flux_p_E
        allocate(eq%flux_p_E(grp_n_r,0:max_deriv+1))
        
        ! flux_t_E
        allocate(eq%flux_t_E(grp_n_r,0:max_deriv+1))
        
        ! q_saf_E
        allocate(eq%q_saf_E(grp_n_r,0:max_deriv+1))
        
        ! rot_t_E
        allocate(eq%rot_t_E(grp_n_r,0:max_deriv+1))
        
        ! rho
        allocate(eq%rho(grid%grp_n_r))
        
        ! magnetic shear
        allocate(eq%S(n_par,n_geo,grp_n_r))
        
        ! normal curvature
        allocate(eq%kappa_n(n_par,n_geo,grp_n_r))
        
        ! geodesic curvature
        allocate(eq%kappa_g(n_par,n_geo,grp_n_r))
        
        ! parallel current
        allocate(eq%sigma(n_par,n_geo,grp_n_r))
        
        ! initialize variables that are specificic to which equilibrium style is
        ! being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! R
                allocate(eq%R_E(n_par,n_geo,grp_n_r,&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                
                ! Z
                allocate(eq%Z_E(n_par,n_geo,grp_n_r,&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                
                ! lambda
                allocate(eq%L_E(n_par,n_geo,grp_n_r,&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function create_eq
    
    ! deallocates equilibrium quantities
    subroutine dealloc_eq(eq)
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium to be deallocated
        
        ! deallocate allocated pointers
        deallocate(eq%flux_p_E,eq%flux_t_E)
        deallocate(eq%flux_p_FD,eq%flux_t_FD)
        
        ! nullify pointers
        nullify(eq%flux_p_E,eq%flux_t_E)
        nullify(eq%flux_p_FD,eq%flux_t_FD)
        
        ! deallocate allocatable variables
        call dealloc_eq_final(eq)
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_eq_final(eq)
            ! input / output
            type(eq_type), intent(out) :: eq                                    ! equilibrium to be deallocated
        end subroutine dealloc_eq_final
    end subroutine dealloc_eq
end module eq_vars
