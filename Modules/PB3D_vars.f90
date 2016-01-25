!------------------------------------------------------------------------------!
!   Variables from PB3D output:                                                !
!       - Miscellaneous variables                                              !
!       - Equilibrium variables                                                !
!       - Perturbation variables                                               !
!       - Solution variables                                                   !
!------------------------------------------------------------------------------!
module PB3D_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use HDF5_vars, only: var_1D_type

    implicit none
    private
    public vars_1D_misc, vars_1D_grid_eq, vars_1D_grid_eq_B, vars_1D_grid_X, &
        &vars_1D_grid_X_B, vars_1D_grid_sol, vars_1D_eq, vars_1D_X_1, &
        &vars_1D_X_2, vars_1D_sol, min_PB3D_version, PB3D_version
    
    ! global variables
    type(var_1D_type), allocatable :: vars_1D_misc(:)                           ! 1D miscellaneous variables
    type(var_1D_type), allocatable :: vars_1D_grid_eq(:)                        ! 1D equilibrium grid variables
    type(var_1D_type), allocatable :: vars_1D_grid_eq_B(:)                      ! 1D field-aligned equilibrium grid variables
    type(var_1D_type), allocatable :: vars_1D_grid_X(:)                         ! 1D perturbation grid variables
    type(var_1D_type), allocatable :: vars_1D_grid_X_B(:)                       ! 1D field-aligned perturbation grid variables
    type(var_1D_type), allocatable :: vars_1D_grid_sol(:)                       ! 1D solution grid variables
    type(var_1D_type), allocatable :: vars_1D_eq(:)                             ! 1D equilibrium variables
    type(var_1D_type), allocatable :: vars_1D_X_1(:)                            ! 1D vectorial perturbation variables
    type(var_1D_type), allocatable :: vars_1D_X_2(:)                            ! 1D tensorial perturbation variables
    type(var_1D_type), allocatable :: vars_1D_sol(:)                            ! 1D solution variables
    real(dp), parameter :: min_PB3D_version = 1.06_dp                           ! minimum PB3D version
    real(dp) :: PB3D_version                                                    ! version of PB3D variable read
end module PB3D_vars
