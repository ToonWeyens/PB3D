!------------------------------------------------------------------------------!
!   Variables from PB3D output:                                                !
!       - Equilibrium variables                                                !
!       - Perturbation variables                                               !
!------------------------------------------------------------------------------!
module PB3D_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type
    use HDF5_ops, only: var_1D

    implicit none
    private
    public dealloc_PB3D, &
        &PB3D_type, vars_1D_eq, vars_1D_eq_B, vars_1D_X, min_PB3D_version
    
    ! PB3D type
    type :: PB3D_type
        real(dp) :: version                                                     ! version
        real(dp) :: alpha                                                       ! field line label alpha
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_type) :: X                                                       ! perturbation variables
    end type PB3D_type
    
    ! global variables
    type(var_1D), allocatable :: vars_1D_eq(:), vars_1D_eq_B(:), vars_1D_X(:)   ! 1D variables
    real(dp), parameter :: min_PB3D_version = 0.80_dp                           ! minimum PB3D version
    
contains
    ! deallocates PB3D quantities
    subroutine dealloc_PB3D(PB3D)
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X
        
        ! input / output
        type(PB3D_type), intent(inout) :: PB3D                                  ! PB3D variables to be deallocated
        
        ! deallocate individual components
        call dealloc_grid(PB3D%grid_eq)
        call dealloc_grid(PB3D%grid_X)
        call dealloc_eq(PB3D%eq)
        call dealloc_met(PB3D%met)
        call dealloc_X(PB3D%X)
    end subroutine dealloc_PB3D
end module PB3D_vars
