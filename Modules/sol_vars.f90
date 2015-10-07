!------------------------------------------------------------------------------!
!   Variables pertaining to the solution quantities                            !
!------------------------------------------------------------------------------!
module sol_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln, iu
    use grid_vars, only: grid_type
    
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
        integer, allocatable :: n(:)                                            ! vector of poloidal mode numbers
        integer, allocatable :: m(:)                                            ! vector of poloidal mode numbers
        complex(dp), allocatable :: vec(:,:,:)                                  ! Eigenvector solution
        complex(dp), allocatable :: val(:)                                      ! Eigenvalue solution
    end type
    
contains
    ! Create a solution type and allocate  the variables. The number of modes as
    ! well as n and m set up.
    ! Optionally, the secondary mode number can be specified (m if poloidal flux
    ! is used  and n  if toroidal  flux). By  default, they  are taken  from the
    ! global X_vars variables.
    integer function create_sol(grid_X,sol,n_EV,lim_sec_X) result(ierr)
        use num_vars, only: use_pol_flux_F
        use X_vars, only: min_n_X, max_n_X, min_m_X, max_m_X
        
        character(*), parameter :: rout_name = 'create_sol'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        integer, intent(in) :: n_EV                                             ! nr. of Eigenvalues
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
        ! local variables
        integer :: id                                                           ! counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: lim_n_X(2)                                                   ! min. and max. of n_X
        integer :: lim_m_X(2)                                                   ! min. and max. of m_X
        
        ! initialize ierr
        ierr = 0
        
        ! set lim_n_X and lim_m_X
        if (use_pol_flux_F) then
            lim_n_X = [min_n_X,max_n_X]
        else
            if (present(lim_sec_X)) then
                lim_m_X = lim_sec_X
            else
                lim_m_X = [min_m_X,max_m_X]
            end if
        end if
        
        ! set n_mod
        sol%n_mod = (lim_m_X(2)-lim_m_X(1)+1)*(lim_n_X(2)-lim_n_X(1)+1)
        
        ! set n and m
        allocate(sol%n(sol%n_mod),sol%m(sol%n_mod))
        if (use_pol_flux_F) then
            sol%n = lim_n_X(1)
            sol%m = [(id, id = lim_m_X(1), lim_m_X(2))]
        else
            sol%n = [(id, id = lim_n_X(1), lim_n_X(2))]
            sol%m = lim_m_X(1)
        end if
        
        ! allocate val
        allocate(sol%val(1:n_EV))
        
        ! allocate vec
        allocate(sol%vec(1:sol%n_mod,1:grid_X%loc_n_r,1:n_EV))
    end function create_sol
    
    ! deallocates solution variables
    ! Note: intent(out) automatically deallocates the variable
    subroutine dealloc_sol(sol)
        ! input / output
        type(sol_type), intent(out) :: sol                                      ! solution variables to be deallocated
    end subroutine dealloc_sol
end module sol_vars
