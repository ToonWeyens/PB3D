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
!! point only, which a coarse grid does not see. One case has a normal
!! coordinate that decreases, as for a VMEC file in which iota and the
!! toroidal flux have opposite signs. A grid on which the window moves over
!! its full width between two points is refused.
!!
!! Two more tests are about the routines that use the tables:
!!
!!  - interp_V (driver_sol) between two grids that start in different
!!    windows, so that a mode sits in another slot in the output grid than in
!!    the input grid: every element of the interpolated matrices must be the
!!    one of its two modes, also where the two slots are in the other order;
!!  - calc_tot_sol_vec and calc_loc_sol_vec (sol_utilities) on a grid with
!!    rows that are empty before its first point, between two points and
!!    after its last point: an empty row may not touch any point.
!------------------------------------------------------------------------------!
module test_modes
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, pi, use_pol_flux_F, X_style, norm_disc_prec_sol
    use eq_vars, only: max_flux_F
    use num_utilities, only: c
    use str_utilities, only: i2str
    use grid_vars, only: grid_type
    use X_vars, only: modes_type, X_2_type, n_mod_X, prim_X, min_n_X, &
        &max_n_X, min_m_X, max_m_X
    use X_ops, only: setup_modes
    use driver_sol, only: interp_V
    use sol_utilities, only: calc_tot_sol_vec, calc_loc_sol_vec

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
            &new_unittest("decreasing_coordinate", test_decreasing), &
            &new_unittest("excursion_between_points", test_excursion), &
            &new_unittest("move_over_full_width_refused", test_refused), &
            &new_unittest("interp_V_other_slots", test_interp_V), &
            &new_unittest("sol_vec_empty_rows", test_sol_vec) &
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
        real(dp) :: r_dir                                                       ! 1 if the normal coordinate increases, -1 if it decreases

        r_dir = sign(1._dp,r_eq(size(r_eq))-r_eq(1))
        do kd = 1,grid%n(3)
            ! window at this point: linear interpolation of w_eq, rounded
            do ld = 1,size(r_eq)-2
                if (r_dir*grid%r_F(kd).le.r_dir*r_eq(ld+1)) exit
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

    !> A normal coordinate that decreases from the first to the last point,
    !! as the poloidal flux of a VMEC file in which iota and the toroidal flux
    !! have opposite signs, with a window that moves up and then down.
    subroutine test_decreasing(error)
        type(error_type), allocatable, intent(out) :: error

        call check_tables(error,[3,4,5,5,4,3,3],3,-eqd(7,0._dp,1._dp),&
            &reshape([-eqd(4,0._dp,1._dp),-eqd(4,0.1_dp,0.95_dp)],[4,2]))
        if (allocated(error)) return
        call check_tables(error,[3,4,5,5,4,3,3],3,-eqd(7,0._dp,1._dp),&
            &reshape([-eqd(61,0._dp,1._dp)],[61,1]))
    end subroutine test_decreasing

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

    !> interp_V between an input grid and an output grid that starts one
    !! window further, with a window that moves up and then down. The input
    !! matrices hold a known function of the two mode numbers, Hermitian for
    !! the symmetric quantities and without symmetry for the others, and do
    !! not depend on the normal coordinate, so that the interpolation itself
    !! is exact.
    !!
    !! \note In the Debug build interp_V also plots the normal ranges. Without
    !! a folder Data this prints an error of HDF5, which does not affect the
    !! result.
    subroutine test_interp_V(error)
        type(error_type), allocatable, intent(out) :: error

        type(grid_type) :: grid_i, grid_o                                       ! input (equilibrium) and output grid
        type(modes_type) :: mds_i, mds_o                                        ! their modes
        type(X_2_type) :: X_i, X_o                                              ! matrices on the two grids
        integer :: ierr, kd, a, b                                               ! error status, counters
        integer :: n_flip                                                       ! pairs of modes with their slots in the other order
        integer, parameter :: n_mod = 4                                         ! width of the window
        integer, parameter :: w_eq(9) = [3,3,4,4,5,5,4,4,3]                     ! window at the input points
        real(dp), parameter :: tol = 1.E-12_dp                                  ! tolerance

        max_flux_F = 2*pi                                                       ! the Debug build of interp_V plots normal ranges divided by it
        call set_window(w_eq,n_mod)
        ierr = make_grid(grid_i,eqd(9,0._dp,1._dp))
        call check(error, ierr, 0, 'input grid')
        if (allocated(error)) return
        ierr = make_grid(grid_o,eqd(12,0.2_dp,0.9_dp))
        call check(error, ierr, 0, 'output grid')
        if (allocated(error)) return
        ierr = setup_modes(mds_i,grid_i,grid_i)
        call check(error, ierr, 0, 'setup_modes on the input grid')
        if (allocated(error)) return
        ierr = setup_modes(mds_o,grid_i,grid_o)
        call check(error, ierr, 0, 'setup_modes on the output grid')
        if (allocated(error)) return
        call check(error, mds_o%m(1,1).ne.mds_i%m(1,1), &
            &'the two grids start in the same window')
        if (allocated(error)) return

        ! input matrices
        call X_i%init(mds_i,grid_i,is_field_averaged=.true.)
        call X_o%init(mds_o,grid_o,is_field_averaged=.true.)
        do kd = 1,grid_i%n(3)
            do b = 1,n_mod
                do a = 1,n_mod
                    if (a.ge.b) then
                        X_i%PV_0(1,1,kd,c([a,b],.true.,n_mod)) = &
                            &herm(mds_i%m(kd,a),mds_i%m(kd,b))
                        X_i%PV_2(1,1,kd,c([a,b],.true.,n_mod)) = &
                            &2*herm(mds_i%m(kd,a),mds_i%m(kd,b))
                        X_i%KV_0(1,1,kd,c([a,b],.true.,n_mod)) = &
                            &3*herm(mds_i%m(kd,a),mds_i%m(kd,b))
                        X_i%KV_2(1,1,kd,c([a,b],.true.,n_mod)) = &
                            &4*herm(mds_i%m(kd,a),mds_i%m(kd,b))
                    end if
                    X_i%PV_1(1,1,kd,c([a,b],.false.,n_mod)) = &
                        &asym(mds_i%m(kd,a),mds_i%m(kd,b))
                    X_i%KV_1(1,1,kd,c([a,b],.false.,n_mod)) = &
                        &2*asym(mds_i%m(kd,a),mds_i%m(kd,b))
                end do
            end do
        end do

        ierr = interp_V(mds_i,grid_i,X_i,mds_o,grid_o,X_o)
        call check(error, ierr, 0, 'interp_V')
        if (allocated(error)) return

        ! every element in the output grid belongs to its two modes
        do kd = 1,grid_o%n(3)
            do b = 1,n_mod
                do a = 1,n_mod
                    if (a.ge.b) then
                        call check(error, &
                            &abs(X_o%PV_0(1,1,kd,c([a,b],.true.,n_mod))-&
                            &herm(mds_o%m(kd,a),mds_o%m(kd,b))).lt.tol .and. &
                            &abs(X_o%PV_2(1,1,kd,c([a,b],.true.,n_mod))-&
                            &2*herm(mds_o%m(kd,a),mds_o%m(kd,b))).lt.tol .and. &
                            &abs(X_o%KV_0(1,1,kd,c([a,b],.true.,n_mod))-&
                            &3*herm(mds_o%m(kd,a),mds_o%m(kd,b))).lt.tol .and. &
                            &abs(X_o%KV_2(1,1,kd,c([a,b],.true.,n_mod))-&
                            &4*herm(mds_o%m(kd,a),mds_o%m(kd,b))).lt.tol, &
                            &'symmetric quantity at point '//trim(i2str(kd))//&
                            &' for modes '//trim(i2str(mds_o%m(kd,a)))//&
                            &' and '//trim(i2str(mds_o%m(kd,b))))
                        if (allocated(error)) return
                    end if
                    call check(error, &
                        &abs(X_o%PV_1(1,1,kd,c([a,b],.false.,n_mod))-&
                        &asym(mds_o%m(kd,a),mds_o%m(kd,b))).lt.tol .and. &
                        &abs(X_o%KV_1(1,1,kd,c([a,b],.false.,n_mod))-&
                        &2*asym(mds_o%m(kd,a),mds_o%m(kd,b))).lt.tol, &
                        &'asymmetric quantity at point '//trim(i2str(kd))//&
                        &' for modes '//trim(i2str(mds_o%m(kd,a)))//&
                        &' and '//trim(i2str(mds_o%m(kd,b))))
                    if (allocated(error)) return
                end do
            end do
        end do

        ! the case is only a test if some pair of modes has its slots in
        ! the other order in the two grids
        n_flip = 0
        do a = 1,size(mds_i%sec,1)
            do b = 1,size(mds_i%sec,1)
                if ((mds_i%sec(a,4).lt.mds_i%sec(b,4)) .neqv. &
                    &(mds_o%sec(a,4).lt.mds_o%sec(b,4))) n_flip = n_flip + 1
            end do
        end do
        call check(error, n_flip.gt.0, 'no pair of modes with its slots in &
            &the other order')

        call X_i%dealloc()
        call X_o%dealloc()
        call mds_i%dealloc()
        call mds_o%dealloc()
        call grid_i%dealloc()
        call grid_o%dealloc()
    contains
        ! Hermitian in the two mode numbers
        complex(dp) function herm(m1,m2)
            integer, intent(in) :: m1, m2
            herm = cmplx(m1+m2,m1-m2,dp)
        end function herm
        ! without symmetry
        complex(dp) function asym(m1,m2)
            integer, intent(in) :: m1, m2
            asym = cmplx(m1,10*m2,dp)
        end function asym
    end subroutine test_interp_V

    !> The solution vector per slot and per mode on a grid with rows that are
    !! empty: the mode of the first equilibrium window has left before the
    !! first point, an excursion of the window falls between two points, and
    !! a mode enters after the last point. Every slot holds its mode times
    !! 1 + r, so that the value and the derivative with respect to r are
    !! known.
    subroutine test_sol_vec(error)
        type(error_type), allocatable, intent(out) :: error

        type(grid_type) :: grid_eq, grid                                        ! equilibrium and solution grid
        type(modes_type) :: mds                                                 ! modes of the solution grid
        complex(dp), allocatable :: vec_loc(:,:), vec_tot(:,:), vec_back(:,:)   ! solution vector per slot, per mode, and back
        integer :: ierr, kd, ld, id                                             ! error status, counters
        integer :: md                                                           ! mode number
        integer :: min_mode                                                     ! lowest mode number of the table
        integer, parameter :: n_mod = 3                                         ! width of the window
        integer, parameter :: w_eq(10) = [3,4,4,4,5,4,4,4,4,5]                  ! window at the equilibrium points
        real(dp), parameter :: r(5) = [0.15_dp,0.3_dp,0.6_dp,0.75_dp,0.85_dp]   ! solution grid
        real(dp), parameter :: tol = 1.E-12_dp                                  ! tolerance
        complex(dp) :: ref                                                      ! reference value

        X_style = 2
        norm_disc_prec_sol = 1
        call set_window(w_eq,n_mod)
        ierr = make_grid(grid_eq,eqd(10,0._dp,1._dp))
        call check(error, ierr, 0, 'equilibrium grid')
        if (allocated(error)) return
        ierr = make_grid(grid,r)
        call check(error, ierr, 0, 'solution grid')
        if (allocated(error)) return
        ierr = setup_modes(mds,grid_eq,grid)
        call check(error, ierr, 0, 'setup_modes')
        if (allocated(error)) return
        call check(error, count(mds%sec(:,3).lt.1).gt.0 .and. &
            &count(mds%sec(:,2).gt.size(r)).gt.0 .and. &
            &count(mds%sec(:,3).lt.mds%sec(:,2)).gt.2, &
            &'no empty rows before, between and after the points')
        if (allocated(error)) return
        min_mode = minval(mds%sec(:,1))

        ! derivative id = 0 and 1
        allocate(vec_loc(n_mod,size(r)))
        do id = 0,1
            do kd = 1,size(r)
                vec_loc(:,kd) = mds%m(kd,:)*(1._dp+r(kd))
            end do
            ierr = calc_tot_sol_vec(mds,grid,vec_loc,vec_tot,deriv=id)
            call check(error, ierr, 0, 'calc_tot_sol_vec')
            if (allocated(error)) return
            do kd = 1,size(r)
                do ld = 1,size(vec_tot,1)
                    md = min_mode+ld-1
                    ref = 0._dp
                    if (any(mds%m(kd,:).eq.md)) ref = md*(1-id)*(1._dp+r(kd)) + &
                        &md*id
                    call check(error, abs(vec_tot(ld,kd)-ref).lt.tol, &
                        &'derivative '//trim(i2str(id))//': mode '//&
                        &trim(i2str(md))//' at point '//trim(i2str(kd)))
                    if (allocated(error)) return
                end do
            end do
            ierr = calc_loc_sol_vec(mds,grid%i_min,vec_tot,vec_back)
            call check(error, ierr, 0, 'calc_loc_sol_vec')
            if (allocated(error)) return
            do kd = 1,size(r)
                do ld = 1,n_mod
                    ref = mds%m(kd,ld)*(1-id)*(1._dp+r(kd)) + mds%m(kd,ld)*id
                    call check(error, abs(vec_back(ld,kd)-ref).lt.tol, &
                        &'derivative '//trim(i2str(id))//': slot '//&
                        &trim(i2str(ld))//' at point '//trim(i2str(kd)))
                    if (allocated(error)) return
                end do
            end do
            deallocate(vec_tot,vec_back)
        end do

        call mds%dealloc()
        call grid%dealloc()
        call grid_eq%dealloc()
    end subroutine test_sol_vec
end module test_modes
