!------------------------------------------------------------------------------!
!   Variables that have to do with equilibrium quantities and the grid used in !
!   the calculations.                                                          !
!   Notes about XDMF:                                                          !
!       -  The  equilibrium  variables  are comprised  of  the  variable  that !
!       result from the equilibrium  calculation, such as pressure, rotational !
!       transform, etc.                                                        !
!       -  These variables  are tabulated  on  a 2D  grid, in  the normal  and !
!       parallel coordinate (thus following the magnetic field lines):         !
!           + In the normal coordinate,  they are tabulated in the equilibrium !
!           grid.                                                              !
!           +  In  the   parallel  coordinate,  they  are   tabulated  in  the !
!           perturbation grid (VMEC),  or in the equilibrium  grid followed by !
!           an adaptation to the perturbation grid (HELENA).                   !
!       - However, this  is not necessarily always the case,  as it depends on !
!       the grid angles being aligned with  the grid (e.g. Providing theta and !
!       zeta in  the equilibrium  grid so  that the  magnetic field  lines are !
!       followed.                                                              !
!------------------------------------------------------------------------------!
module eq_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, pi, max_str_ln
    use messages, only: writo, lvl_ud, print_ar_2, print_ar_1
    use str_ops, only: i2str
    use grid_vars, only: grid_type

    implicit none
    private
    public create_eq, dealloc_eq, dealloc_eq_final, &
        &eq_type, &
        &R_0, pres_0, B_0, psi_0, rho_0, T_0, max_flux_p_E, max_flux_p_F, &
        &max_flux_t_E, max_flux_t_F
    
    ! global variables
    real(dp) :: R_0, pres_0, rho_0                                              ! independent normalization constants for nondimensionalization
    real(dp) :: B_0, psi_0, T_0                                                 ! derived normalization constants for nondimensionalization
    real(dp) :: max_flux_p_E, max_flux_p_F                                      ! max. pol. flux in Equilibrium and Flux coordinates
    real(dp) :: max_flux_t_E, max_flux_t_F                                      ! max. tor. flux in Equilibrium and Flux coordinates
    
    ! equilibrium type
    ! The arrays here are of the form (except for rho):
    !   - (r,Dr)                                    for flux quantities
    !   - (angle_1,angle_2,r,Dr,Dangle_1,Dangle_2)  for normal quantities
    ! where angle_1 and angle_2 can be any angle that completely describe a flux
    ! surface. For example, they  can refer to a grid of  theta and zeta values,
    ! but they  can also refer  to (Modified)  flux coordinates with  a parallel
    ! angle and  a field line coordinate  (alpha). By choosing angle_1  equal to
    ! the  parallel coordinate  and angle_2  to the  field line  coordinate, the
    ! second dimension of  the matrices is then  chosen to be of size  1 for the
    ! calculations  on a  single  field line.  At the  same  time, the  parallel
    ! coordinate, in which integrations will have  to be done, ocupies the first
    ! index. This is good for numerical efficiency.
    ! An  equilibrium type  should be  complemented  by grid  type, which  gives
    ! information about the grid points  at which the equilibrium quantities are
    ! tabulated.
    type :: eq_type
        real(dp), allocatable :: R_E(:,:,:,:,:,:)                               ! R in E(quilibrium) coord
        real(dp), allocatable :: Z_E(:,:,:,:,:,:)                               ! Z in E(quilibrium) coords.
        real(dp), allocatable :: L_E(:,:,:,:,:,:)                               ! L(ambda) in E(quilibrium) coords.
        real(dp), allocatable :: rho(:)                                         ! density (in all coord. systems)
        real(dp), allocatable :: pres_E(:,:)                                    ! pressure, and norm. Deriv. in E(equilibrium) coords.
        real(dp), allocatable :: q_saf_E(:,:)                                   ! safety factor in E(equilibrium) coordinates
        real(dp), allocatable :: rot_t_E(:,:)                                   ! rot. transform in E(equilibrium) coordinates
        real(dp), pointer :: flux_p_E(:,:), flux_t_E(:,:)                       ! pol. and tor. flux and norm. Deriv. in E(equilibrium) coords.
        real(dp), allocatable :: pres_FD(:,:)                                   ! pressure, and norm. Deriv. with values and Derivs. in F(lux) coords.
        real(dp), allocatable :: q_saf_FD(:,:)                                  ! safety factor, Deriv. in F(lux) coords.
        real(dp), allocatable :: rot_t_FD(:,:)                                  ! rot. transform, Deriv. in F(lux) coords.
        real(dp), pointer :: flux_p_FD(:,:), flux_t_FD(:,:)                     ! pol. and tor. flux, and norm. Deriv. with values and Derivs. in F(lux) coords.
    end type

contains
    ! creates new equilibrium
    ! The normal and angular grid can be  in any coord. system, as only the grid
    ! sizes are used, not the coordinate values.
    integer function create_eq(eq,grid) result(ierr)
        use num_vars, only: max_deriv, eq_style
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'create_eq'
        
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium to be created
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: grp_n_r, n                                                   ! group and total nr. of normal points
        integer :: n_par, n_geo                                                 ! tot. nr. of angular points in parallel and geodesic direction
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Create equilibrium...')
        call lvl_ud(1)
        
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
        
        call lvl_ud(-1)
    end function create_eq
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! equilibrium phase
    integer function dealloc_eq(eq) result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'dealloc_eq'
        
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium to be deallocated
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! deallocate general variables
        deallocate(eq%flux_p_E)
        deallocate(eq%flux_t_E)
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                deallocate(eq%R_E,eq%Z_E,eq%L_E)
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function dealloc_eq
    
    ! deallocates  equilibrium quantities  that are not  used anymore  after the
    ! calculation for a certain alpha
    subroutine dealloc_eq_final(eq)
        ! input / output
        type(eq_type), intent(inout) :: eq                                      ! equilibrium to be deallocated
        
        ! deallocate general variables
        deallocate(eq%pres_FD)
        deallocate(eq%flux_p_FD,eq%flux_t_FD)
        deallocate(eq%q_saf_FD)
        deallocate(eq%rot_t_FD)
        
        deallocate(eq%pres_E)
        deallocate(eq%q_saf_E)
        deallocate(eq%rot_t_E)
        
        deallocate(eq%rho)
    end subroutine dealloc_eq_final
end module eq_vars
