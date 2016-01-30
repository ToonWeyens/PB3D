!------------------------------------------------------------------------------!
!   Variables pertaining to the solution quantities                            !
!------------------------------------------------------------------------------!
module sol_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln, iu
    use grid_vars, only: grid_type
    use output_ops
    
    implicit none
    
    private
    public dealloc_sol, create_sol, &
        &sol_type
    
    ! solution type
    ! The arrays here are of the form:
    !   - val: 1:n_EV
    !   - vec: (1:n_mod,1:loc_n_r,1:n_EV)
    type :: sol_type
        integer :: n_mod                                                        ! size of n and m (nr. of modes)
        integer, allocatable :: n(:,:)                                          ! vector of poloidal mode numbers
        integer, allocatable :: m(:,:)                                          ! vector of poloidal mode numbers
        complex(dp), allocatable :: vec(:,:,:)                                  ! Eigenvector solution
        complex(dp), allocatable :: val(:)                                      ! Eigenvalue solution
    end type
    
contains
    ! Create a solution type and allocate  the variables. The number of modes as
    ! well as n and m set up.
    ! Optionally, the secondary mode number can be specified (m if poloidal flux
    ! is used  and n  if toroidal  flux). By  default, they  are taken  from the
    ! global X_vars variables.
    ! Note: The lowest limits of the grid  need to be 1; e.g. grid_sol%i_min = 1
    ! for first process.
    subroutine create_sol(grid_sol,sol,n_EV,lim_sec_X)
        use X_vars, only: set_nm_X
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        integer, intent(in) :: n_EV                                             ! nr. of Eigenvalues
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! set mode numbers
        call set_nm_X(grid_sol,sol%n,sol%m,lim_sec_X)
        
        ! set n_mod
        sol%n_mod = size(sol%n,2)
        
        ! allocate val
        allocate(sol%val(1:n_EV))
        
        ! allocate vec
        allocate(sol%vec(1:sol%n_mod,1:grid_sol%loc_n_r,1:n_EV))
    end subroutine create_sol
    
    ! deallocates solution variables
    ! Note: intent(out) automatically deallocates the variable
    subroutine dealloc_sol(sol)
        ! input / output
        type(sol_type), intent(out) :: sol                                      ! solution variables to be deallocated
    end subroutine dealloc_sol
end module sol_vars
