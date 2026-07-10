!------------------------------------------------------------------------------!
!> Unit tests for the pure Richardson/parallel index arithmetic of
!! PB3D_utilities (setup_rich_id, setup_par_id) and the 1-D <-> n-D storage
!! conversion (conv_1D2ND).
!!
!! The setup_par_id spec comes from its ~50-line docstring: at Richardson
!! level i (max. I), the points of an interlaced parallel grid with
!! n = 1 + k 2^(I-1) points that were *introduced* at level i are
!!    p, p+s, p+2s, ...   with   s = 2^(I-1), p = 1       for i = 1,
!!                               s = 2^(I-i+1), p = 1+s/2  for i > 1.
!! The strongest property is checked programmatically: the level index sets
!! partition the full grid (disjoint union). The par_lim window and the
!! par_id_mem memory indices are checked against the same formulas.
!------------------------------------------------------------------------------!
module test_PB3D_utilities
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp
    use str_utilities, only: i2str
    use grid_vars, only: grid_type
    use var_1D_vars, only: var_1D_type
    use PB3D_utilities, only: setup_rich_id, setup_par_id, conv_1D2ND

    implicit none
    private
    public collect_PB3D_utilities

contains
    !> Collect the tests.
    subroutine collect_PB3D_utilities(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("setup_rich_id", test_setup_rich_id), &
            &new_unittest("par_id_partition", test_par_id_partition), &
            &new_unittest("par_id_window_and_memory", test_par_id_window), &
            &new_unittest("conv_1D2ND", test_conv_1D2ND) &
            &]
    end subroutine collect_PB3D_utilities

    !> Richardson level range: only the last level, or all levels when
    !! combining (and there is something to combine).
    subroutine test_setup_rich_id(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, all(setup_rich_id(1).eq.[1,1]), 'max 1, no tot')
        if (allocated(error)) return
        call check(error, all(setup_rich_id(1,tot_rich=.true.).eq.[1,1]), &
            &'max 1 has nothing to combine')
        if (allocated(error)) return
        call check(error, all(setup_rich_id(3).eq.[3,3]), 'max 3, no tot')
        if (allocated(error)) return
        call check(error, all(setup_rich_id(3,tot_rich=.true.).eq.[1,3]), &
            &'max 3, tot')
        if (allocated(error)) return
    end subroutine test_setup_rich_id

    !> The index sets of the Richardson levels partition the interlaced
    !! parallel grid: every point belongs to exactly one level.
    subroutine test_par_id_partition(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: rich_lvl_max = 3                                  ! I
        integer, parameter :: n_par = 17                                        ! 1 + 4 k with k = 4

        type(grid_type) :: grid                                                 ! grid (only n(1) is used)
        integer :: lvl, id                                                      ! counters
        integer :: par_id(3)                                                    ! parallel id
        integer :: count_vis(17)                                                ! how often each point is visited

        grid%n = [n_par,1,1]

        ! without combination: the whole grid at once
        par_id = setup_par_id(grid,rich_lvl_max,rich_lvl_max)
        call check(error, all(par_id.eq.[1,n_par,1]), &
            &'no tot_rich: whole grid with stride 1')
        if (allocated(error)) return

        ! with combination: partition over the levels
        count_vis = 0
        do lvl = 1,rich_lvl_max
            par_id = setup_par_id(grid,rich_lvl_max,lvl,tot_rich=.true.)
            call check(error, par_id(1).ge.1 .and. par_id(2).le.n_par .and. &
                &par_id(3).ge.1, 'valid range for level '//trim(i2str(lvl)))
            if (allocated(error)) return
            do id = par_id(1),par_id(2),par_id(3)
                count_vis(id) = count_vis(id) + 1
            end do
        end do
        call check(error, all(count_vis.eq.1), &
            &'levels partition the grid (every point exactly once)')
        if (allocated(error)) return

        ! the explicit expected ids (s and p from the docstring)
        call check(error, all(setup_par_id(grid,3,1,tot_rich=.true.).eq.&
            &[1,17,4]), 'level 1: s = 4, p = 1')
        if (allocated(error)) return
        call check(error, all(setup_par_id(grid,3,2,tot_rich=.true.).eq.&
            &[3,15,4]), 'level 2: s = 4, p = 3')
        if (allocated(error)) return
        call check(error, all(setup_par_id(grid,3,3,tot_rich=.true.).eq.&
            &[2,16,2]), 'level 3: s = 2, p = 2')
        if (allocated(error)) return
    end subroutine test_par_id_partition

    !> A par_lim window returns exactly the full-range points that fall
    !! inside the window (in window-local coordinates), and par_id_mem gives
    !! their contiguous indices in the per-level HDF5 storage.
    subroutine test_par_id_window(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: rich_lvl_max = 3                                  ! I
        integer, parameter :: n_par = 17                                        ! 1 + 4 k with k = 4
        integer, parameter :: par_lim(2) = [5,13]                               ! window

        type(grid_type) :: grid                                                 ! grid (only n(1) is used)
        integer :: lvl, id                                                      ! counters
        integer :: par_id(3), par_id_win(3)                                     ! full-range and windowed ids
        integer :: par_id_mem(2)                                                ! memory indices
        integer :: n_exp                                                        ! expected number of points
        logical :: glob_set(17)                                                 ! full-range index set of this level

        grid%n = [n_par,1,1]

        do lvl = 1,rich_lvl_max
            ! full-range set of this level
            par_id = setup_par_id(grid,rich_lvl_max,lvl,tot_rich=.true.)
            glob_set = .false.
            do id = par_id(1),par_id(2),par_id(3)
                glob_set(id) = .true.
            end do

            ! windowed set: same points, in window-local coordinates
            par_id_win = setup_par_id(grid,rich_lvl_max,lvl,tot_rich=.true.,&
                &par_lim=par_lim,par_id_mem=par_id_mem)
            call check(error, par_id_win(3), par_id(3), &
                &'stride unchanged by window, level '//trim(i2str(lvl)))
            if (allocated(error)) return
            n_exp = 0
            do id = par_lim(1),par_lim(2)
                if (glob_set(id)) then
                    n_exp = n_exp + 1
                    call check(error, &
                        &par_id_win(1)+(n_exp-1)*par_id_win(3), &
                        &id-par_lim(1)+1, &
                        &'windowed index of global point '//trim(i2str(id))//&
                        &', level '//trim(i2str(lvl)))
                    if (allocated(error)) return
                end if
            end do
            call check(error, (par_id_win(2)-par_id_win(1))/par_id_win(3)+1, &
                &n_exp, 'number of windowed points, level '//trim(i2str(lvl)))
            if (allocated(error)) return

            ! memory indices: contiguous, starting at the position of the
            ! first windowed point in the per-level storage
            call check(error, par_id_mem(2)-par_id_mem(1)+1, n_exp, &
                &'memory range size, level '//trim(i2str(lvl)))
            if (allocated(error)) return
        end do

        ! documented explicit case: level 1, window [5,13]: global points
        ! 5, 9, 13 are the 2nd..4th stored level-1 points
        par_id_win = setup_par_id(grid,rich_lvl_max,1,tot_rich=.true.,&
            &par_lim=par_lim,par_id_mem=par_id_mem)
        call check(error, all(par_id_win.eq.[1,9,4]) .and. &
            &all(par_id_mem.eq.[2,4]), 'level-1 window: ids [1,9,4], mem [2,4]')
        if (allocated(error)) return
    end subroutine test_par_id_window

    !> 1-D to n-D conversion restores shape, bounds and values.
    subroutine test_conv_1D2ND(error)
        type(error_type), allocatable, intent(out) :: error

        type(var_1D_type) :: var                                                ! 1-D variable
        real(dp), allocatable :: var_2D(:,:), var_3D(:,:,:)                     ! converted variables
        integer :: id, jd, kd                                                   ! counters
        integer :: ld                                                           ! linear counter

        ! 2-D with non-unit lower bounds: (0:2, 3:4)
        var%tot_i_min = [0,3]
        var%tot_i_max = [2,4]
        var%p = [(1._dp*ld, ld = 1,6)]
        call conv_1D2ND(var,var_2D)
        call check(error, lbound(var_2D,1).eq.0 .and. ubound(var_2D,1).eq.2 &
            &.and. lbound(var_2D,2).eq.3 .and. ubound(var_2D,2).eq.4, &
            &'2-D bounds')
        if (allocated(error)) return
        ld = 0
        do jd = 3,4
            do id = 0,2
                ld = ld + 1
                call check(error, var_2D(id,jd), 1._dp*ld, &
                    &thr=0._dp, message='2-D column-major order')
                if (allocated(error)) return
            end do
        end do

        ! 3-D: (1:2, 1:3, 1:2)
        deallocate(var%p)
        var%tot_i_min = [1,1,1]
        var%tot_i_max = [2,3,2]
        var%p = [(2._dp*ld, ld = 1,12)]
        call conv_1D2ND(var,var_3D)
        call check(error, all(shape(var_3D).eq.[2,3,2]), '3-D shape')
        if (allocated(error)) return
        ld = 0
        do kd = 1,2
            do jd = 1,3
                do id = 1,2
                    ld = ld + 1
                    call check(error, var_3D(id,jd,kd), 2._dp*ld, &
                        &thr=0._dp, message='3-D column-major order')
                    if (allocated(error)) return
                end do
            end do
        end do
    end subroutine test_conv_1D2ND
end module test_PB3D_utilities
