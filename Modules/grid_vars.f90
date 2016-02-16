!------------------------------------------------------------------------------!
!   Variables pertaining to the different grids used                           !
!------------------------------------------------------------------------------!
module grid_vars
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln

    implicit none
    private
    public create_grid, dealloc_grid, grid_type, create_disc, dealloc_disc, &
        &disc_type, &
        &alpha, n_r_eq, n_par_X, min_par_X, max_par_X

    ! type for grids
    ! The grids  are saved in  the following format:  (angle_1,angle_2,r), where
    ! angle_1  and angle_2  can be  any angle  that completely  describe a  flux
    ! surface. For example, they  can refer to a grid of  theta and zeta values,
    ! but they  can also refer  to (Modified)  flux coordinates with  a parallel
    ! angle and a field line coordinate (alpha).
    ! Normally, for field-aligned grids, angle_1 is chosen equal to the parallel
    ! coordinate in the  Flux coordinate system, and angle_2 equal  to the field
    ! line label (so the  second dimension of the matrices is  then chosen to be
    ! of size 1 for the calculations on  a single field line). At the same time,
    ! the  parallel coordinate,  in which  integrations  will have  to be  done,
    ! ocupies the first index. This is good for numerical efficiency.
    ! For specific  equilibrium grids, such  as the case for  HELENA equilibria,
    ! angle_1 corresponds  normally to theta  and angle_2 to the  symmetry angle
    ! zeta. 
    ! It is important that this order  of the coordinates in space is consistent
    ! among all the variables.
    type :: grid_type
        integer :: n(3)                                                         ! tot nr. of points
        integer :: loc_n_r                                                      ! local nr. of normal points
        integer :: i_min, i_max                                                 ! min. and max. normal index of this process in full arrays
        logical :: divided                                                      ! whether the grid is split up among the processes
        real(dp), pointer :: r_E(:) => null()                                   ! E(quilibrium) coord. values at n points (should be intialized to n)
        real(dp), pointer :: r_F(:) => null()                                   ! F(lux) coord. values at n points (should be intialized to n)
        real(dp), pointer :: loc_r_E(:) => null()                               ! E(quilibrium) coord. values at n points (should be intialized to n)
        real(dp), pointer :: loc_r_F(:) => null()                               ! F(lux) coord. values at n points (should be intialized to n)
        real(dp), pointer :: theta_E(:,:,:) => null()                           ! E(quilibrium) coord. values of first angle at n points in 3D
        real(dp), pointer :: theta_F(:,:,:) => null()                           ! F(lux) coord. values of first angle at n points in 3D
        real(dp), pointer :: zeta_E(:,:,:) => null()                            ! E(quilibrium) coord. values of second angle at n points in 3D
        real(dp), pointer :: zeta_F(:,:,:) => null()                            ! F(lux) coord. values of second angle at n points in 3D
        real(dp), allocatable :: trigon_factors(:,:,:,:,:,:)                    ! trigonometric factor cosine for the inverse fourier transf.
    end type
    
    ! type for data of discretization operations:
    !   - derivatives
    !   - interpolation
    ! See routines setup_deriv_data and setup_interp_data in grid_utilities where
    ! discretization data is setup and apply_disc where it is used.
    type :: disc_type
        integer :: n, n_loc                                                     ! total and local size of discretization variables
        real(dp), allocatable :: dat(:,:)                                       ! nonzero elements of matrix corresponding to discretization
        integer, allocatable :: id_start(:)                                     ! start index of data in dat
    end type
    
    ! global variables
    ! (These  variables  should  be   used   only  until  the  grids  have  been
    ! established. They are put here for lack of a better place.)
    real(dp) :: alpha                                                           ! field line label alpha
    integer :: n_r_eq                                                           ! nr. of normal points in equilibrium grid
    integer :: n_par_X                                                          ! nr. of parallel points in field-aligned grid
    real(dp) :: min_par_X, max_par_X                                            ! min. and max. of parallel coordinate [pi] in field-aligned grid
    
    ! interfaces
    interface create_grid
        module procedure create_grid_3D, create_grid_1D
    end interface
    
contains
    ! creates a new grid
    ! Optionally, the local limits can be provided for a divided grid.
    ! Note: intent(out) automatically deallocates the variable
    integer function create_grid_3D(grid,n,i_lim) result(ierr)                  ! 3D version
        character(*), parameter :: rout_name = 'create_grid_3D'
        
        ! input / output
        type(grid_type), intent(out) :: grid                                    ! grid to be created
        integer :: n(3)                                                         ! tot. nr. of points (par,r,alpha)
        integer, optional :: i_lim(2)                                           ! min. and max. local normal index
        
        ! local Variables
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
        end if
        
        ! check if 1D version has to be called
        if (n(1).eq.0 .and. n(2).eq.0) then
            ierr = create_grid_1D(grid,n(3),i_lim)
            CHCKERR('')
        else
            ! set divided and loc_n_r
            grid%divided = .false.
            if (present(i_lim)) then                                            ! might be divided grid
                grid%i_min = i_lim(1)
                grid%i_max = i_lim(2)
                grid%loc_n_r = i_lim(2)-i_lim(1)+1
                if (i_lim(2)-i_lim(1)+1.lt.n(3)) grid%divided = .true.          ! only divided if difference between local and total
            else                                                                ! certainly not divided grid
                grid%loc_n_r = n(3)
            end if
            
            ! copy the values
            grid%n = n
            allocate(grid%theta_E(n(1),n(2),grid%loc_n_r))
            allocate(grid%zeta_E(n(1),n(2),grid%loc_n_r))
            allocate(grid%theta_F(n(1),n(2),grid%loc_n_r))
            allocate(grid%zeta_F(n(1),n(2),grid%loc_n_r))
            allocate(grid%r_E(n(3)))
            allocate(grid%r_F(n(3)))
            if (grid%divided) then
                allocate(grid%loc_r_E(grid%loc_n_r))
                allocate(grid%loc_r_F(grid%loc_n_r))
            else
                grid%loc_r_E => grid%r_E
                grid%loc_r_F => grid%r_F
            end if
        end if
    end function create_grid_3D
    integer function create_grid_1D(grid,n,i_lim) result(ierr)                  ! 1D version
        character(*), parameter :: rout_name = 'create_grid_1D'
        
        ! input / output
        type(grid_type), intent(out) :: grid                                    ! grid to be created
        integer :: n                                                            ! tot. nr. of points (par,r,alpha)
        integer, optional :: i_lim(2)                                           ! min. and max. local normal index
        
        ! local Variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (present(i_lim)) then
            if (i_lim(2)-i_lim(1)+1.gt.n) then
                ierr = 1
                err_msg = 'The local nr. of normal points cannot be higher &
                    &than the total nr. of normal points'
                CHCKERR(err_msg)
            end if
        end if
        
        ! set divided and loc_n_r
        grid%divided = .false.
        if (present(i_lim)) then                                                ! might be divided grid
            grid%i_min = i_lim(1)
            grid%i_max = i_lim(2)
            grid%loc_n_r = i_lim(2)-i_lim(1)+1
            if (i_lim(2)-i_lim(1)+1.lt.n) grid%divided = .true.                 ! only divided if difference between local and total
        else                                                                    ! certainly not divided grid
            grid%loc_n_r = n
        end if
        
        ! copy the values
        grid%n = [0,0,n]
        allocate(grid%r_E(n))
        allocate(grid%r_F(n))
        if (grid%divided) then
            allocate(grid%loc_r_E(grid%loc_n_r))
            allocate(grid%loc_r_F(grid%loc_n_r))
        else
            grid%loc_r_E => grid%r_E
            grid%loc_r_F => grid%r_F
        end if
    end function create_grid_1D
    
    ! deallocates a grid
    subroutine dealloc_grid(grid)
        ! input / output
        type(grid_type), intent(inout) :: grid                                  ! grid to be deallocated
        
        ! deallocate allocated pointers
        deallocate(grid%r_E,grid%r_F)
        if (grid%n(1).ne.0 .and. grid%n(2).ne.0) then                           ! 3D grid
            deallocate(grid%theta_E,grid%zeta_E)
            deallocate(grid%theta_F,grid%zeta_F)
        end if
        if (grid%divided) deallocate(grid%loc_r_E,grid%loc_r_F)
        
        ! nullify pointers
        nullify(grid%r_E,grid%r_F)
        nullify(grid%theta_E,grid%zeta_E)
        nullify(grid%theta_F,grid%zeta_F)
        nullify(grid%loc_r_E,grid%loc_r_F)
        
        ! deallocate allocatable variables
        call dealloc_grid_final(grid)
    contains
        ! Note: intent(out) automatically deallocates the variable
        subroutine dealloc_grid_final(grid)
            ! input / output
            type(grid_type), intent(out) :: grid                                ! grid to be deallocated
        end subroutine dealloc_grid_final
    end subroutine dealloc_grid
    
    ! Create discretization variable, overwriting.
    integer function create_disc(deriv,n,n_loc) result(ierr)
        character(*), parameter :: rout_name = 'create_disc'
        
        ! input / output
        type(disc_type), intent(inout) :: deriv                                 ! discretization variable
        integer, intent(in) :: n, n_loc                                         ! total and local size of discretization
        
        ! local Variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (n.lt.1 .or. n_loc.lt.1) then
            ierr = 1
            err_msg = 'n and n_loc need to be > 1'
            CHCKERR(err_msg)
        end if
        
        ! (re)allocate
        if (allocated(deriv%dat)) deallocate(deriv%dat)
        if (allocated(deriv%id_start)) deallocate(deriv%id_start)
        allocate(deriv%dat(n,n_loc))
        allocate(deriv%id_start(n))
        deriv%dat = 0._dp
        deriv%id_start = 0
        
        ! set n, n_loc
        deriv%n = n
        deriv%n_loc = n_loc
    end function create_disc
    
    ! Deallocate discretization variable type
    subroutine dealloc_disc(disc_data)
        ! input / output
        type(disc_type), intent(out) :: disc_data                               ! discretization variable to be deallocated
    end subroutine dealloc_disc
end module grid_vars
