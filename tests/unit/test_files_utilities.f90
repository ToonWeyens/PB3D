!------------------------------------------------------------------------------!
!> Unit tests for files_utilities.
!!
!! Temporary files are created in the current working directory (the CTest
!! working dir) and removed afterwards.
!------------------------------------------------------------------------------!
module test_files_utilities
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp, max_str_ln
    use files_utilities

    implicit none
    private
    public collect_files_utilities

    character(len=*), parameter :: tmp_name = 'pb3d_test_tmp.txt'               ! scratch file name

contains
    !> Collect all tests of this suite.
    subroutine collect_files_utilities(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("get_full_PB3D_name", test_get_full_PB3D_name), &
            &new_unittest("nextunit", test_nextunit), &
            &new_unittest("count_lines_and_skip_comment", test_count_skip), &
            &new_unittest("delete_file", test_delete_file) &
            &]
    end subroutine collect_files_utilities

    !> get_full_PB3D_name: with and without Richardson level appendix
    subroutine test_get_full_PB3D_name(error)
        use num_vars, only: PB3D_name

        type(error_type), allocatable, intent(out) :: error

        PB3D_name = 'PB3D_out'
        call check(error, trim(get_full_PB3D_name()), 'PB3D_out.h5')
        if (allocated(error)) return
        call check(error, trim(get_full_PB3D_name(rich_lvl=2)), &
            &'PB3D_out_R_2.h5')
        if (allocated(error)) return
        ! non-positive level is ignored
        call check(error, trim(get_full_PB3D_name(rich_lvl=0)), 'PB3D_out.h5')
    end subroutine test_get_full_PB3D_name

    !> nextunit: returns a free unit in [70,1000], skips occupied ones
    subroutine test_nextunit(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: u1, u2                                                       ! units
        integer :: u1_out                                                       ! output argument version

        u1 = nextunit(u1_out)
        call check(error, u1.ge.70 .and. u1.le.1000, &
            &'nextunit outside [70,1000]')
        if (allocated(error)) return
        call check(error, u1, u1_out, 'function result and out argument differ')
        if (allocated(error)) return

        ! occupy u1, so the next call has to return something else
        open(UNIT=u1,FILE=tmp_name,STATUS='replace')
        u2 = nextunit()
        call check(error, u1.ne.u2, 'nextunit returned an occupied unit')
        close(UNIT=u1,STATUS='delete')
    end subroutine test_nextunit

    !> count_lines counts non-comment lines; skip_comment positions the file
    !! at the first non-comment line
    subroutine test_count_skip(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: file_i                                                       ! file unit
        integer :: ierr                                                         ! error status
        integer :: val                                                          ! value read

        ! set up scratch file: 3 data lines, 2 comment lines
        open(UNIT=nextunit(file_i),FILE=tmp_name,STATUS='replace')
        write(file_i,'(A)') '# leading comment'
        write(file_i,'(A)') '1'
        write(file_i,'(A)') '# interior comment'
        write(file_i,'(A)') '2'
        write(file_i,'(A)') '3'
        rewind(file_i)

        call check(error, count_lines(file_i), 3, 'count_lines miscounted')
        if (allocated(error)) then
            close(UNIT=file_i,STATUS='delete')
            return
        end if

        ! skip_comment must leave the file positioned at the first data line
        ierr = skip_comment(file_i,tmp_name)
        call check(error, ierr, 0, 'skip_comment failed')
        if (allocated(error)) then
            close(UNIT=file_i,STATUS='delete')
            return
        end if
        read(file_i,*) val
        call check(error, val, 1, 'skip_comment not at first data line')

        close(UNIT=file_i,STATUS='delete')
    end subroutine test_count_skip

    !> delete_file removes an existing file
    subroutine test_delete_file(error)
        type(error_type), allocatable, intent(out) :: error

        integer :: file_i                                                       ! file unit
        integer :: istat                                                        ! status
        logical :: exists                                                       ! file existence
        character(len=max_str_ln) :: file_name                                  ! mutable file name

        ! create scratch file
        open(UNIT=nextunit(file_i),FILE=tmp_name,STATUS='replace')
        write(file_i,'(A)') 'to be deleted'
        close(file_i)

        file_name = tmp_name
        istat = delete_file(file_name)
        call check(error, istat, 0, 'delete_file failed')
        if (allocated(error)) return

        inquire(FILE=tmp_name,EXIST=exists)
        call check(error, .not.exists, 'file still exists after delete_file')
    end subroutine test_delete_file
end module test_files_utilities
