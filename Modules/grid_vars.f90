!------------------------------------------------------------------------------!
!> Variables pertaining to the different grids used.
!------------------------------------------------------------------------------!
module grid_vars
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, pi, max_str_ln, weight_dp

    implicit none
    private
    public n_r_in, n_r_eq, n_r_X, n_r_sol, min_par_X, max_par_X, n_alpha, &
        &min_alpha, max_alpha, alpha
#if ldebug
    public n_alloc_grids, n_alloc_discs
#endif
    
    ! global variables
    integer :: n_r_in                                                           !< nr. of normal points in input grid
    integer :: n_r_eq                                                           !< nr. of normal points in equilibrium grid
    integer :: n_r_X                                                            !< nr. of normal points in perturbation grid
    integer :: n_r_sol                                                          !< nr. of normal points in solution grid
    integer :: n_alpha                                                          !< nr. of field-lines
    real(dp) :: min_par_X                                                       !< min. of parallel coordinate [\f$\pi\f$] in field-aligned grid
    real(dp) :: max_par_X                                                       !< max. of parallel coordinate [\f$\pi\f$] in field-aligned grid
    real(dp) :: min_alpha                                                       !< min. of field-line label [\f$\pi\f$] in field-aligned grid
    real(dp) :: max_alpha                                                       !< max. of field-line label [\f$\pi\f$] in field-aligned grid
    real(dp), allocatable :: alpha(:)                                           !< field line label alpha
#if ldebug
    integer :: n_alloc_grids                                                    !< nr. of allocated grids \ldebug
    integer :: n_alloc_discs                                                    !< nr. of allocated discretizations \ldebug
#endif

    !> Type for grids
    !!
    !! The     grids     are     saved      in     the     following     format:
    !! <tt>(angle_1,angle_2,r)</tt>, where \c angle_1 and  \c angle_2 can be any
    !! angle that completely describe a flux surface.
    !! 
    !! For example,  they can refer  to a  grid of \f$\theta\f$  and \f$\zeta\f$
    !! values, but  they can also  refer to  (Modified) flux coordinates  with a
    !! parallel angle and a field line coordinate (\f$\alpha\f$).
    !!
    !! For  field-aligned grids,  \c angle_1  is generally  chosen equal  to the
    !! parallel coordinate in  the Flux coordinate system, and  \c angle_2 equal
    !! to the field line label (so the  second dimension of the matrices is then
    !! chosen to be of size 1 for the calculations on a single field line).
    !!
    !! At the  same time,  the parallel coordinate,  in which  integrations will
    !! have to  be done,  ocupies the  first index. This  is good  for numerical
    !! efficiency.
    !!
    !! For specific equilibrium  grids, such as the case  for HELENA equilibria,
    !! \c angle_1  corresponds to  \f$\theta\f$ and \c  angle_2 to  the symmetry
    !! angle \f$\zeta\f$.
    !!
    !! It is important that this order of the coordinates in space is consistent
    !! among all the variables.
    type, public :: grid_type
        integer :: n(3)                                                         !< tot nr. of points
        integer :: loc_n_r                                                      !< local nr. of normal points
        integer :: i_min                                                        !< min. normal index of this process in full arrays
        integer :: i_max                                                        !< max. normal index of this process in full arrays
        logical :: divided                                                      !< whether several processes each hold a part of the grid
        real(dp), pointer :: r_E(:) => null()                                   !< E(quilibrium) coord. values at n points 
        real(dp), pointer :: r_F(:) => null()                                   !< F(lux) coord. values at n points 
        real(dp), pointer :: loc_r_E(:) => null()                               !< E(quilibrium) coord. values of this process: r_E(i_min:i_max)
        real(dp), pointer :: loc_r_F(:) => null()                               !< F(lux) coord. values of this process: r_F(i_min:i_max)
        real(dp), pointer :: theta_E(:,:,:) => null()                           !< E(quilibrium) coord. values of first angle at n points in 3-D
        real(dp), pointer :: theta_F(:,:,:) => null()                           !< F(lux) coord. values of first angle at n points in 3-D
        real(dp), pointer :: zeta_E(:,:,:) => null()                            !< E(quilibrium) coord. values of second angle at n points in 3-D
        real(dp), pointer :: zeta_F(:,:,:) => null()                            !< F(lux) coord. values of second angle at n points in 3-D
        real(dp), allocatable :: trigon_factors(:,:,:,:,:)                      !< trigonometric basis for the inverse Fourier transform: cos and sin of m theta and of n nfp zeta (see VMEC_utilities.calc_trigon_factors())
#if ldebug
        real(dp) :: estim_mem_usage                                             !< estimated memory usage \ldebug
#endif
    contains
        !> initialize
        procedure :: init => init_grid
        !> copy
        procedure :: copy => copy_grid
        !> deallocate
        procedure :: dealloc => dealloc_grid
    end type
    
contains
    !> \public Initializes a new grid.
    !!
    !! Optionally, the local limits can be provided, if this process holds only
    !! a part of the grid.
    !!
    !! The normal coordinates \c r_E and \c r_F of the whole grid are allocated
    !! on every process. \c loc_r_E  and \c loc_r_F are pointers to their part
    !! between the  local limits,  so that  setting the  coordinates of the whole
    !! grid sets the local ones as well.
    !!
    !! \c divided states whether several processes each hold a part of the grid,
    !! so that the processes  have to communicate to treat the whole grid. It has
    !! to be the same on every process and can therefore not be derived from the
    !! local limits. The limits and \c divided are independent:
    !!  - whole grid, not divided: an ordinary grid of one process, or a grid
    !!  that every process holds completely for itself.
    !!  - part of the grid, divided: the usual case on several processes.
    !!  - whole grid, divided: a process whose part, with its ghost regions, is
    !!  the whole grid. It takes part in the communication as the others do.
    !!  - part of the grid, not divided: a grid restricted to a part of which
    !!  no other process holds the rest, as POST on one process with a reduced
    !!  normal range. trim_grid() returns that part.
    !!
    !! \return ierr
    integer function init_grid(grid,n,i_lim,divided) result(ierr)
#if ldebug
        use num_vars, only: print_mem_usage, rank
#endif
        character(*), parameter :: rout_name = 'init_grid'
        
        ! input / output
        class(grid_type), intent(inout) :: grid                                 !< grid to be initialized
        integer, intent(in) :: n(3)                                             !< tot. nr. of points (par,r,alpha)
        integer, intent(in), optional :: i_lim(2)                               !< min. and max. local normal index
        logical, intent(in) :: divided                                          !< whether several processes each hold a part of the grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (present(i_lim)) then
            if (i_lim(2)-i_lim(1)+1.gt.n(3)) then
                ierr = 1
                err_msg = 'The local nr. of normal points cannot be higher &
                    &than the total nr. of normal points'
                CHCKERR(err_msg)
            end if
            if (i_lim(1).lt.1 .or. i_lim(2).gt.n(3)) then
                ierr = 1
                err_msg = 'The local normal limits '//trim(i2str(i_lim(1)))//&
                    &'..'//trim(i2str(i_lim(2)))//' do not lie in the grid &
                    &of '//trim(i2str(n(3)))//' normal points'
                CHCKERR(err_msg)
            end if
        end if
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) grid%estim_mem_usage = 0._dp
#endif
        
        ! set local limits and loc_n_r
        if (present(i_lim)) then                                                ! part of the grid
            if (i_lim(2).lt.i_lim(1)) then
                ierr = 1
                write(*,*) 'i_lim = ', i_lim
                err_msg = 'faulty i_lim'
                CHCKERR(err_msg)
            end if
            grid%i_min = i_lim(1)
            grid%i_max = i_lim(2)
            grid%loc_n_r = i_lim(2)-i_lim(1)+1
        else                                                                    ! whole grid
            grid%i_min = 1
            grid%i_max = n(3)
            grid%loc_n_r = n(3)
        end if
        grid%divided = divided
        
        ! allocate normal variables, the local ones point at the local part
        grid%n = n
        allocate(grid%r_E(n(3)))
        allocate(grid%r_F(n(3)))
        grid%loc_r_E => grid%r_E(grid%i_min:grid%i_max)
        grid%loc_r_F => grid%r_F(grid%i_min:grid%i_max)
        
#if ldebug
        ! set estimated memory usage
        if (print_mem_usage) grid%estim_mem_usage = grid%estim_mem_usage + &
            &grid%n(3)*2
#endif
        
        ! allocate angular variables
        if (n(1).ne.0 .and. n(2).ne.0) then
            allocate(grid%theta_E(n(1),n(2),grid%loc_n_r))
            allocate(grid%zeta_E(n(1),n(2),grid%loc_n_r))
            allocate(grid%theta_F(n(1),n(2),grid%loc_n_r))
            allocate(grid%zeta_F(n(1),n(2),grid%loc_n_r))
            
#if ldebug
            ! set estimated memory usage
            if (print_mem_usage) grid%estim_mem_usage = grid%estim_mem_usage + &
                &grid%loc_n_r*n(1)*n(2)*4
#endif
        end if
        
#if ldebug
        ! increment n_alloc_grids
        n_alloc_grids = n_alloc_grids + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of grid: '//&
            &trim(r2strt(grid%estim_mem_usage*weight_dp))//' kB]',&
            &alert=.true.)
#endif
    end function init_grid
    
    !> \public Deallocates a grid.
    subroutine dealloc_grid(grid)
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        ! input / output
        class(grid_type), intent(inout) :: grid                                 !< grid to be deallocated
        
        ! local variables
#if ldebug
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = grid%estim_mem_usage
        end if
#endif
        
        ! deallocate allocated pointers
        deallocate(grid%r_E,grid%r_F)
        if (grid%n(1).ne.0 .and. grid%n(2).ne.0) then                           ! 3D grid
            deallocate(grid%theta_E,grid%zeta_E)
            deallocate(grid%theta_F,grid%zeta_F)
        end if
        
        ! nullify pointers
        nullify(grid%r_E,grid%r_F)
        nullify(grid%theta_E,grid%zeta_E)
        nullify(grid%theta_F,grid%zeta_F)
        nullify(grid%loc_r_E,grid%loc_r_F)
        
        ! deallocate allocatable variables
        call dealloc_grid_final(grid)
        
#if ldebug
        ! decrement n_alloc_grids
        n_alloc_grids = n_alloc_grids - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating grid ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        !> \private
        subroutine dealloc_grid_final(grid)
            ! input / output
            type(grid_type), intent(out) :: grid                                ! grid to be deallocated
        end subroutine dealloc_grid_final
    end subroutine dealloc_grid
    
    !> \public Deep copy of a grid.
    !!
    !! \note Does not copy possible trigoniometric factors.
    !!
    !! \return ierr
    integer function copy_grid(grid_i,grid_o) result(ierr)
        character(*), parameter :: rout_name = 'copy_grid'
        
        ! input / output
        class(grid_type), intent(in) :: grid_i                                  !< grid to be copied
        type(grid_type), intent(inout) :: grid_o                                !< copied grid
        
        ! initialize ierr
        ierr = 0
        
        ierr = init_grid(grid_o,grid_i%n,i_lim=[grid_i%i_min,grid_i%i_max],&
            &divided=grid_i%divided)
        CHCKERR('')
        grid_o%r_E = grid_i%r_E
        grid_o%r_F = grid_i%r_F
        if (associated(grid_i%theta_E)) grid_o%theta_E = grid_i%theta_E
        if (associated(grid_i%theta_F)) grid_o%theta_F = grid_i%theta_F
        if (associated(grid_i%zeta_E)) grid_o%zeta_E = grid_i%zeta_E
        if (associated(grid_i%zeta_F)) grid_o%zeta_F = grid_i%zeta_F
        if (allocated(grid_i%trigon_factors)) then
            allocate(grid_o%trigon_factors(size(grid_i%trigon_factors,1),&
                &size(grid_i%trigon_factors,2),size(grid_i%trigon_factors,3),&
                &size(grid_i%trigon_factors,4),size(grid_i%trigon_factors,5)))
            grid_o%trigon_factors = grid_i%trigon_factors
        end if
    end function copy_grid
end module grid_vars
