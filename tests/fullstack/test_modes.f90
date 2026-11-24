!------------------------------------------------------------------------------!
!> Full-stack tests of the mode tables of X_style 2 (X_ops.setup_modes).
!!
!! The window of n_mod_X poloidal modes is prescribed on the equilibrium
!! points (min_m_X), as init_modes would set it from n q, and setup_modes
!! builds the tables n, m and sec on the equilibrium grid and on other grids.
!! The tests assert, independently of the construction in setup_modes:
!!
!!  - at every grid point, the rows of sec that contain the point hold each
!!    mode of the window exactly once, and every one of them sits in its own
!!    slot of m;
!!  - every grid has the same rows of sec as the equilibrium grid, in the same
!!    order, which interp_V and calc_XUQ rely on to pair them;
!!  - no row is empty on the equilibrium grid.
!!
!! The windows move up, down, up and then down (a minimum of the safety
!! factor), down and then up, and once up and back down at one equilibrium
!! point only, which a coarse grid does not see. A grid on which the window
!! moves over its full width between two points is refused.
!------------------------------------------------------------------------------!
module test_modes
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, use_pol_flux_F
    use str_utilities, only: i2str
    use grid_vars, only: grid_type
    use X_vars, only: modes_type, n_mod_X, prim_X, min_n_X, max_n_X, &
        &min_m_X, max_m_X
    use X_ops, only: setup_modes

    implicit none
    private
    public collect_modes

contains
    !> Collect the tests.
    subroutine collect_modes(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("window_up", test_up), &
            &new_unittest("window_down", test_down), &
            &new_unittest("window_up_then_down", test_up_down), &
            &new_unittest("window_down_then_up", test_down_up), &
            &new_unittest("excursion_between_points", test_excursion), &
            &new_unittest("move_over_full_width_refused", test_refused) &
            &]
    end subroutine collect_modes

    !--------------------------------------------------------------------------
    ! helpers
    !--------------------------------------------------------------------------

    !> A grid with only normal points.
    integer function make_grid(grid,r) result(ierr)
        type(grid_type), intent(inout) :: grid                                  !< grid
        real(dp), intent(in) :: r(:)                                            !< normal coordinate

        ierr = grid%init([1,1,size(r)],divided=.false.)
        if (ierr.ne.0) return
        grid%r_F = r
        grid%r_E = r
    end function make_grid

    !> Equally spaced points from a to b.
    function eqd(n,a,b) result(r)
        integer, intent(in) :: n                                                !< number of points
        real(dp), intent(in) :: a, b                                            !< first and last point
        real(dp) :: r(n)
        integer :: kd                                                           ! counter

        r = [(a+(b-a)*(kd-1)/(n-1), kd=1,n)]
    end function eqd

    !> Prescribe the window on the equilibrium points, as init_modes does for
    !! X_style 2 and poloidal flux.
    subroutine set_window(w_eq,n_mod)
        integer, intent(in) :: w_eq(:)                                          !< lowest mode at the equilibrium points
        integer, intent(in) :: n_mod                                            !< width of the window

        use_pol_flux_F = .true.
        n_mod_X = n_mod
        prim_X = 10
        if (allocated(min_n_X)) deallocate(min_n_X,max_n_X,min_m_X,max_m_X)
        allocate(min_n_X(size(w_eq)),max_n_X(size(w_eq)))
        min_n_X = prim_X
        max_n_X = prim_X
        min_m_X = w_eq
        max_m_X = w_eq + n_mod - 1
    end subroutine set_window

    !> Set up the tables on the equilibrium grid and on the other grids, and
    !! check them.
    subroutine check_tables(error,w_eq,n_mod,r_eq,grids)
        type(error_type), allocatable, intent(out) :: error
        integer, intent(in) :: w_eq(:)                                          !< lowest mode at the equilibrium points
        integer, intent(in) :: n_mod                                            !< width of the window
        real(dp), intent(in) :: r_eq(:)                                         !< equilibrium points
        real(dp), intent(in) :: grids(:,:)                                      !< other grids, one per column (normal points from 0 to 1)

        type(grid_type) :: grid_eq, grid                                        ! equilibrium and other grid
        type(modes_type) :: mds_eq, mds                                         ! their modes
        integer :: ierr, id                                                     ! error status, counter

        call set_window(w_eq,n_mod)
        ierr = make_grid(grid_eq,r_eq)
        call check(error, ierr, 0, 'equilibrium grid')
        if (allocated(error)) return
        ierr = setup_modes(mds_eq,grid_eq,grid_eq)
        call check(error, ierr, 0, 'setup_modes on the equilibrium grid')
        if (allocated(error)) return
        call check_grid(error,mds_eq,grid_eq,r_eq,w_eq,'equilibrium grid')
        if (allocated(error)) return
        call check(error, all(mds_eq%sec(:,3).ge.mds_eq%sec(:,2)), &
            &'empty row on the equilibrium grid')
        if (allocated(error)) return

        do id = 1,size(grids,2)
            ierr = make_grid(grid,grids(:,id))
            call check(error, ierr, 0, 'grid '//trim(i2str(id)))
            if (allocated(error)) return
            ierr = setup_modes(mds,grid_eq,grid)
            call check(error, ierr, 0, 'setup_modes on grid '//trim(i2str(id)))
            if (allocated(error)) return
            call check_grid(error,mds,grid,r_eq,w_eq,'grid '//trim(i2str(id)))
            if (allocated(error)) return
            call check(error, size(mds%sec,1).eq.size(mds_eq%sec,1), &
                &'grid '//trim(i2str(id))//': '//trim(i2str(size(mds%sec,1)))//&
                &' rows against '//trim(i2str(size(mds_eq%sec,1))))
            if (allocated(error)) return
            call check(error, all(mds%sec(:,1).eq.mds_eq%sec(:,1)), &
                &'grid '//trim(i2str(id))//': rows differ from the &
                &equilibrium grid')
            if (allocated(error)) return
            call mds%dealloc()
            call grid%dealloc()
        end do

        call mds_eq%dealloc()
        call grid_eq%dealloc()
    end subroutine check_tables

    !> At every point of a grid, the rows that contain the point hold the
    !! window of that point, each mode once and in its own slot.
    subroutine check_grid(error,mds,grid,r_eq,w_eq,name)
        type(error_type), allocatable, intent(out) :: error
        type(modes_type), intent(in) :: mds                                     !< modes of the grid
        type(grid_type), intent(in) :: grid                                     !< grid
        real(dp), intent(in) :: r_eq(:)                                         !< equilibrium points
        integer, intent(in) :: w_eq(:)                                          !< lowest mode at the equilibrium points
        character(*), intent(in) :: name                                        !< name of the grid

        integer :: kd, ld, md                                                   ! counters
        integer :: w                                                            ! lowest mode of the window at a point
        integer :: n_in                                                         ! number of rows that contain a mode at a point
        real(dp) :: v                                                           ! interpolated window

        do kd = 1,grid%n(3)
            ! window at this point: linear interpolation of w_eq, rounded
            do ld = 1,size(r_eq)-2
                if (grid%r_F(kd).le.r_eq(ld+1)) exit
            end do
            v = w_eq(ld) + (w_eq(ld+1)-w_eq(ld))*(grid%r_F(kd)-r_eq(ld))/&
                &(r_eq(ld+1)-r_eq(ld))
            w = nint(v)
            do md = w,w+n_mod_X-1
                n_in = count(mds%sec(:,1).eq.md .and. mds%sec(:,2).le.kd .and. &
                    &mds%sec(:,3).ge.kd)
                call check(error, n_in, 1, trim(name)//', point '//&
                    &trim(i2str(kd))//': mode '//trim(i2str(md))//' in '//&
                    &trim(i2str(n_in))//' rows')
                if (allocated(error)) return
            end do
            call check(error, count(mds%sec(:,2).le.kd .and. &
                &mds%sec(:,3).ge.kd), n_mod_X, trim(name)//', point '//&
                &trim(i2str(kd))//': rows outside the window')
            if (allocated(error)) return
            do ld = 1,size(mds%sec,1)
                if (mds%sec(ld,2).le.kd .and. mds%sec(ld,3).ge.kd) then
                    call check(error, mds%m(kd,mds%sec(ld,4)), mds%sec(ld,1), &
                        &trim(name)//', point '//trim(i2str(kd))//&
                        &': wrong slot for mode '//trim(i2str(mds%sec(ld,1))))
                    if (allocated(error)) return
                end if
            end do
        end do
    end subroutine check_grid

    !--------------------------------------------------------------------------
    ! tests
    !--------------------------------------------------------------------------

    !> The window moves up, by two at once on the coarse grid.
    subroutine test_up(error)
        type(error_type), allocatable, intent(out) :: error

        call check_tables(error,[3,3,4,5,5,6,8],3,eqd(7,0._dp,1._dp),&
            &reshape([eqd(51,0._dp,1._dp),eqd(51,0.4_dp,0.9_dp),&
            &eqd(51,0._dp,0.5_dp)],[51,3]))
    end subroutine test_up

    !> The window moves down; the second grid starts after the first move.
    subroutine test_down(error)
        type(error_type), allocatable, intent(out) :: error

        call check_tables(error,[9,9,8,7,7,6,5],4,eqd(7,0._dp,1._dp),&
            &reshape([eqd(40,0._dp,1._dp),eqd(40,0.3_dp,1._dp)],[40,2]))
    end subroutine test_down

    !> The window moves up and then down, as for a minimum of q, and the
    !! second grid has four points only.
    subroutine test_up_down(error)
        type(error_type), allocatable, intent(out) :: error

        call check_tables(error,[3,4,5,5,4,3,3],3,eqd(7,0._dp,1._dp),&
            &reshape([eqd(4,0._dp,1._dp),eqd(4,0.1_dp,0.95_dp)],[4,2]))
        if (allocated(error)) return
        call check_tables(error,[3,4,5,5,4,3,3],3,eqd(7,0._dp,1._dp),&
            &reshape([eqd(61,0._dp,1._dp)],[61,1]))
    end subroutine test_up_down

    !> The window moves down and then up, as for a maximum of q.
    subroutine test_down_up(error)
        type(error_type), allocatable, intent(out) :: error

        call check_tables(error,[7,6,5,4,4,5,6,7],5,eqd(8,0._dp,1._dp),&
            &reshape([eqd(5,0._dp,1._dp),eqd(33,0._dp,1._dp)],[5,2]))
    end subroutine test_down_up

    !> The window moves up and back down at one equilibrium point. A grid
    !! without a point near it has empty rows for that excursion and still
    !! the same rows as the equilibrium grid.
    subroutine test_excursion(error)
        type(error_type), allocatable, intent(out) :: error

        call check_tables(error,[3,3,3,4,3,3,3],3,eqd(7,0._dp,1._dp),&
            &reshape([0._dp,0.3_dp,0.7_dp,1._dp,eqd(4,0._dp,1._dp)],[4,2]))
    end subroutine test_excursion

    !> A grid on which the window moves by more than its width between two
    !! points is refused.
    subroutine test_refused(error)
        type(error_type), allocatable, intent(out) :: error

        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(modes_type) :: mds                                                 ! modes
        integer :: ierr                                                         ! error status

        call set_window([3,7],3)
        ierr = make_grid(grid_eq,[0._dp,1._dp])
        call check(error, ierr, 0, 'equilibrium grid')
        if (allocated(error)) return
        ierr = setup_modes(mds,grid_eq,grid_eq)
        call check(error, ierr.ne.0, 'a move by 4 with 3 modes was accepted')
        call grid_eq%dealloc()
    end subroutine test_refused
end module test_modes
