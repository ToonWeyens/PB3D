!------------------------------------------------------------------------------!
!> Unit tests for the pure mode-index logic of X_utilities: the local-to-
!! total secondary-index translation, the symmetric-storage necessity
!! criterion, the contiguous HDF5 ranges of tensorial perturbation
!! variables (pinning the worked example in the get_sec_X_range
!! documentation), and the input/output mode-range trimming.
!------------------------------------------------------------------------------!
module test_X_utilities
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use str_utilities, only: i2str
    use X_vars, only: n_mod_X, modes_type
    use X_utilities, only: sec_ind_loc2tot, is_necessary_X, get_sec_X_range, &
        &check_same_modes

    implicit none
    private
    public collect_X_utilities

contains
    !> Collect the tests.
    subroutine collect_X_utilities(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("sec_ind_loc2tot", test_sec_ind_loc2tot), &
            &new_unittest("is_necessary_X", test_is_necessary_X), &
            &new_unittest("get_sec_X_range_doc_example", &
            &   test_get_sec_X_range), &
            &new_unittest("check_same_modes", test_check_same_modes) &
            &]
    end subroutine collect_X_utilities

    !> Local mode indices translate to total ones by the offset of the
    !! lower limit.
    subroutine test_sec_ind_loc2tot(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: res2(2)                                                      ! tensorial result

        n_mod_X = 5

        ! vectorial: default limits [1, n_mod_X]
        call check(error, sec_ind_loc2tot(3), 3, 'default limits')
        if (allocated(error)) return
        call check(error, sec_ind_loc2tot(3,lim_sec_X=[2,4]), 4, &
            &'shifted limits')
        if (allocated(error)) return

        ! tensorial: independent limits per dimension
        res2 = sec_ind_loc2tot(2,3,lim_sec_X=&
            &reshape([2,4,3,5],[2,2]))                                          ! dim 1: [2,4], dim 2: [3,5]
        call check(error, res2(1).eq.3 .and. res2(2).eq.5, &
            &'tensorial: got ['//trim(i2str(res2(1)))//','//&
            &trim(i2str(res2(2)))//']')
        if (allocated(error)) return
    end subroutine test_sec_ind_loc2tot

    !> Asymmetric variables are always necessary; symmetric ones only on or
    !! below the diagonal in *total* indices.
    subroutine test_is_necessary_X(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: kd, md                                                       ! counters
        integer, parameter :: lims(2,2) = reshape([2,3,2,5],[2,2])              ! dim 1: [2,3], dim 2: [2,5]
        logical :: res, expected                                                ! result and reference

        n_mod_X = 5

        ! asymmetric: always
        call check(error, is_necessary_X(.false.,[1,5]), &
            &'asymmetric always necessary')
        if (allocated(error)) return

        ! symmetric, full range: on or below the diagonal
        do kd = 1,5
            do md = 1,5
                res = is_necessary_X(.true.,[kd,md])
                expected = kd.ge.md
                call check(error, res.eqv.expected, &
                    &'symmetric full range at ('//trim(i2str(kd))//','//&
                    &trim(i2str(md))//')')
                if (allocated(error)) return
            end do
        end do

        ! symmetric, subrange [2:3, 2:5]: diagonal in total indices
        do kd = 1,2                                                             ! local dim 1: total 2..3
            do md = 1,4                                                         ! local dim 2: total 2..5
                res = is_necessary_X(.true.,[kd,md],lim_sec_X=lims)
                expected = lims(1,1)+kd.ge.lims(1,2)+md
                call check(error, res.eqv.expected, &
                    &'symmetric subrange at ('//trim(i2str(kd))//','//&
                    &trim(i2str(md))//')')
                if (allocated(error)) return
            end do
        end do
    end subroutine test_is_necessary_X

    !> The worked example in the documentation: subrange [2:3, 2:5] of a
    !! total range [1:5, 1:5].
    subroutine test_get_sec_X_range(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: lims(2,2) = reshape([2,3,2,5],[2,2])              ! dim 1: [2,3], dim 2: [2,5]
        integer, parameter :: loc_asym(2,4) = &
            &reshape([1,2,3,4,5,6,7,8],[2,4])                                   ! documented local ranges, asymmetric
        integer, parameter :: tot_asym(2,4) = &
            &reshape([7,8,12,13,17,18,22,23],[2,4])                             ! documented total ranges, asymmetric
        integer, parameter :: loc_sym(2,2) = reshape([1,2,3,3],[2,2])           ! documented local ranges, symmetric (nonempty part)
        integer, parameter :: tot_sym(2,2) = reshape([6,7,10,10],[2,2])         ! documented total ranges, symmetric (nonempty part)

        integer :: md                                                           ! counter
        integer :: range_loc(2), range_tot(2)                                   ! local and total ranges

        n_mod_X = 5

        ! asymmetric: all four columns
        do md = 1,4
            call get_sec_X_range(range_loc,range_tot,md,.false.,&
                &lim_sec_X=lims)
            call check(error, all(range_loc.eq.loc_asym(:,md)) .and. &
                &all(range_tot.eq.tot_asym(:,md)), &
                &'asymmetric ranges for m = '//trim(i2str(md))//': loc ['//&
                &trim(i2str(range_loc(1)))//':'//trim(i2str(range_loc(2)))//&
                &'], tot ['//trim(i2str(range_tot(1)))//':'//&
                &trim(i2str(range_tot(2)))//']')
            if (allocated(error)) return
        end do

        ! symmetric: two nonempty columns, then empty ranges
        do md = 1,2
            call get_sec_X_range(range_loc,range_tot,md,.true.,&
                &lim_sec_X=lims)
            call check(error, all(range_loc.eq.loc_sym(:,md)) .and. &
                &all(range_tot.eq.tot_sym(:,md)), &
                &'symmetric ranges for m = '//trim(i2str(md))//': loc ['//&
                &trim(i2str(range_loc(1)))//':'//trim(i2str(range_loc(2)))//&
                &'], tot ['//trim(i2str(range_tot(1)))//':'//&
                &trim(i2str(range_tot(2)))//']')
            if (allocated(error)) return
        end do
        do md = 3,4
            call get_sec_X_range(range_loc,range_tot,md,.true.,&
                &lim_sec_X=lims)
            call check(error, range_loc(1).gt.range_loc(2), &
                &'symmetric range for m = '//trim(i2str(md))//' not empty')
            if (allocated(error)) return
        end do
    end subroutine test_get_sec_X_range

    !> Two mode tables pass only with the same mode in every row, also when
    !! a mode has more than one row, as for a safety factor with an extremum.
    subroutine test_check_same_modes(error)
        type(error_type), allocatable, intent(out) :: error

        type(modes_type) :: mds_i, mds_o                                        ! input and output modes
        integer :: ierr                                                         ! error status
        integer, parameter :: modes(6) = [4,5,6,7,4,3]                          ! the window moves up and down

        allocate(mds_i%sec(6,4),mds_o%sec(6,4))
        mds_i%sec = 0
        mds_o%sec = 1                                                           ! ranges and slots may differ
        mds_i%sec(:,1) = modes
        mds_o%sec(:,1) = modes
        ierr = check_same_modes(mds_i,mds_o)
        call check(error, ierr, 0, 'same rows refused')
        if (allocated(error)) return

        ! one mode differs
        mds_o%sec(5,1) = 8
        ierr = check_same_modes(mds_i,mds_o)
        call check(error, ierr.ne.0, 'a different mode was accepted')
        if (allocated(error)) return
        call mds_o%dealloc()

        ! one row fewer
        allocate(mds_o%sec(5,4))
        mds_o%sec(:,1) = modes(1:5)
        ierr = check_same_modes(mds_i,mds_o)
        call check(error, ierr.ne.0, 'a missing row was accepted')
        if (allocated(error)) return

        call mds_i%dealloc()
        call mds_o%dealloc()
    end subroutine test_check_same_modes
end module test_X_utilities
