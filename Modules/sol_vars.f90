!------------------------------------------------------------------------------!
!> Variables pertaining to the solution quantities.
!------------------------------------------------------------------------------!
module sol_vars
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, max_str_ln, iu, weight_dp
    use grid_vars, only: grid_type
    use output_ops
    
    implicit none
    
    private
    public alpha
#if ldebug
    public n_alloc_sols
#endif
    
    ! global variables
    real(dp) :: alpha                                                           !< field line label alpha \ldebug
#if ldebug
    integer :: n_alloc_sols                                                     !< nr. of allocated grids \ldebug
#endif
    
    !> solution type
    !!
    !! The arrays here are of the form:
    !!  - \c val: <tt>(1:n_EV)</tt>
    !!  - \c vec: <tt>(1:n_mod,1:loc_n_r,1:n_EV)</tt>
    type, public :: sol_type
        integer :: n_mod                                                        !< size of n and m (nr. of modes)
        integer, allocatable :: n(:,:)                                          !< vector of toroidal mode numbers
        integer, allocatable :: m(:,:)                                          !< vector of poloidal mode numbers
        complex(dp), allocatable :: vec(:,:,:)                                  !< Eigenvector solution
        complex(dp), allocatable :: val(:)                                      !< Eigenvalue solution
#if ldebug
        real(dp) :: estim_mem_usage                                             !< estimated memory usage \ldebug
#endif
    contains
        !> initialize
        procedure :: init => init_sol
        !> deallocate
        procedure :: dealloc => dealloc_sol
    end type
    
contains
    !> \public Initialize a solution type and allocate the variables.
    !!
    !! The number of modes as well as \c n and \c m are also set up.
    !!
    !! Optionally, the secondary mode number can  be specified (\c m if poloidal
    !! flux is used as normal coordinate and  c n if toroidal flux). By default,
    !! they are taken from the global \c X_vars variables.
    !!
    !! \note If the lowest limits of  the grid is not 1 (e.g. <tt>grid_sol%i_min
    !! = 1</tt> for first process), the input variable \c i_min should be set to
    !! set correctly. For a full grid, it should be set to 1.
    subroutine init_sol(sol,grid_sol,n_EV,lim_sec_X)
        use X_vars, only: set_nm_X
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        
        ! input / output
        class(sol_type), intent(inout) :: sol                                   !< solution variables
        type(grid_type), intent(in) :: grid_sol                                 !< solution grid
        integer, intent(in) :: n_EV                                             !< nr. of Eigenvalues
        integer, intent(in), optional :: lim_sec_X(2)                           !< limits of \c m_X (pol. flux) or \c n_X (tor. flux)
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) sol%estim_mem_usage = 0._dp
#endif
        
        ! set mode numbers
        call set_nm_X(grid_sol,sol%n,sol%m,lim_sec_X)
        
        ! set n_mod
        sol%n_mod = size(sol%n,2)
        
        ! allocate val
        allocate(sol%val(1:n_EV))
        
        ! allocate vec
        allocate(sol%vec(1:sol%n_mod,1:grid_sol%loc_n_r,1:n_EV))
        
#if ldebug
        ! set estimated memory usage
        if (print_mem_usage) sol%estim_mem_usage = sol%estim_mem_usage + &
            &n_EV * (1 + sol%n_mod*grid_sol%loc_n_r)
        
        ! increment n_alloc_sols
        n_alloc_sols = n_alloc_sols + 1
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of sol: '//&
            &trim(r2strt(sol%estim_mem_usage*weight_dp*2))//' kB]',alert=.true.)
#endif
    end subroutine init_sol
    
    !> \public Deallocates solution variables.
    !!
    !! \note intent(out) automatically deallocates the variable
    subroutine dealloc_sol(sol)
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        
        ! input / output
        class(sol_type), intent(inout) :: sol                                   !< solution variables to be deallocated
        
#if ldebug
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = sol%estim_mem_usage
        end if
#endif
        
        ! deallocate allocatable variables
        call dealloc_sol_final(sol)
        
#if ldebug
        ! decrement n_alloc_sols
        n_alloc_sols = n_alloc_sols - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating sol ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp*2))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        !> \private
        subroutine dealloc_sol_final(sol)
            ! input / output
            type(sol_type), intent(out) :: sol                                      ! solution to be deallocated
        end subroutine dealloc_sol_final
    end subroutine dealloc_sol
end module sol_vars
