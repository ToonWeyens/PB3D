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
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type
    use HDF5_vars, only: var_1D_type

    implicit none
    private
    public dealloc_PB3D, &
        &PB3D_type, vars_1D_misc, vars_1D_eq, vars_1D_eq_B, vars_1D_X_1, &
        &vars_1D_X_2, vars_1D_sol, min_PB3D_version
    
    ! PB3D type
    type :: PB3D_type
        real(dp) :: version                                                     ! version
        real(dp) :: alpha                                                       ! field line label alpha
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type) :: grid_X                                               ! perturbation grid
        type(eq_type) :: eq                                                     ! equilibrium variables
        type(met_type) :: met                                                   ! metric variables
        type(X_1_type) :: X_1                                                   ! vectorial perturbation variables
        type(X_2_type) :: X_2                                                   ! tensorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
    end type PB3D_type
    
    ! global variables
    type(var_1D_type), allocatable :: vars_1D_misc(:)                           ! 1D miscellaneous variables
    type(var_1D_type), allocatable :: vars_1D_eq(:)                             ! 1D equilibrium variables
    type(var_1D_type), allocatable :: vars_1D_eq_B(:)                           ! 1D field-aligned equilbrium variables
    type(var_1D_type), allocatable :: vars_1D_X_1(:)                            ! 1D vectorial perturbation variables
    type(var_1D_type), allocatable :: vars_1D_X_2(:)                            ! 1D tensorial perturbation variables
    type(var_1D_type), allocatable :: vars_1D_sol(:)                            ! 1D solution variables
    real(dp), parameter :: min_PB3D_version = 0.94_dp                           ! minimum PB3D version
    
contains
    ! deallocates PB3D quantities
    subroutine dealloc_PB3D(PB3D)
        use grid_vars, only: dealloc_grid
        use eq_vars, only: dealloc_eq
        use met_vars, only: dealloc_met
        use X_vars, only: dealloc_X
        use sol_vars, only: dealloc_sol
        
        ! input / output
        type(PB3D_type), intent(inout) :: PB3D                                  ! PB3D variables to be deallocated
        
        ! deallocate individual components
        call dealloc_grid(PB3D%grid_eq)
        call dealloc_grid(PB3D%grid_X)
        call dealloc_eq(PB3D%eq)
        call dealloc_met(PB3D%met)
        call dealloc_X(PB3D%X_1)
        call dealloc_X(PB3D%X_2)
        call dealloc_sol(PB3D%sol)
    end subroutine dealloc_PB3D
end module PB3D_vars
