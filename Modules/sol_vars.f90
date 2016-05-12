!------------------------------------------------------------------------------!
!   Variables pertaining to the solution quantities                            !
!------------------------------------------------------------------------------!
module sol_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln, iu, weight_dp
    use grid_vars, only: grid_type
    use output_ops
    
    implicit none
    
    private
    public sol_type, &
        &alpha
#if ldebug
    public n_alloc_sols
#endif
    
    ! global variables
    real(dp) :: alpha                                                           ! field line label alpha
#if ldebug
    integer :: n_alloc_sols                                                     ! nr. of allocated grids
#endif
    
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
#if ldebug
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
#endif
    contains
        procedure :: init => init_sol
        procedure :: dealloc => dealloc_sol
    end type
    
contains
    ! Initialize a solution type and allocate the variables. The number of modes
    ! as well as n and m set up.
    ! Optionally, the secondary mode number can be specified (m if poloidal flux
    ! is used  and n  if toroidal  flux). By  default, they  are taken  from the
    ! global X_vars variables.
    ! Note: The lowest limits of the grid  need to be 1; e.g. grid_sol%i_min = 1
    ! for first process.
    subroutine init_sol(sol,grid_sol,n_EV,lim_sec_X)
        use X_vars, only: set_nm_X
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        
        ! input / output
        class(sol_type), intent(inout) :: sol                                   ! solution variables
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        integer, intent(in) :: n_EV                                             ! nr. of Eigenvalues
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        
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
    
    ! deallocates solution variables
    ! Note: intent(out) automatically deallocates the variable
    subroutine dealloc_sol(sol)
#if ldebug
        use messages, only: get_mem_usage
        use num_vars, only: rank, print_mem_usage
#endif
        
        ! input / output
        class(sol_type), intent(inout) :: sol                                   ! solution variables to be deallocated
        
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
        subroutine dealloc_sol_final(sol)
            ! input / output
            type(sol_type), intent(out) :: sol                                      ! solution to be deallocated
        end subroutine dealloc_sol_final
    end subroutine dealloc_sol
end module sol_vars
