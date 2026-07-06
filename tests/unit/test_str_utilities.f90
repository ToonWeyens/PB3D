!------------------------------------------------------------------------------!
!> Unit tests for str_utilities.
!------------------------------------------------------------------------------!
module test_str_utilities
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use num_vars, only: dp
    use str_utilities

    implicit none
    private
    public collect_str_utilities

contains
    !> Collect all tests of this suite.
    subroutine collect_str_utilities(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            &new_unittest("i2str", test_i2str), &
            &new_unittest("ii2str", test_ii2str), &
            &new_unittest("r2str", test_r2str), &
            &new_unittest("r2strt", test_r2strt), &
            &new_unittest("c2str", test_c2str), &
            &new_unittest("case_conversion", test_case_conversion), &
            &new_unittest("merge_strings", test_merge_strings) &
            &]
    end subroutine collect_str_utilities

    !> i2str: integer to string
    subroutine test_i2str(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, trim(i2str(42)), "42")
        if (allocated(error)) return
        call check(error, trim(i2str(-7)), "-7")
        if (allocated(error)) return
        call check(error, trim(i2str(0)), "0")
    end subroutine test_i2str

    !> ii2str: kind-8 integer to string
    subroutine test_ii2str(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, trim(ii2str(123456789012_8)), "123456789012")
    end subroutine test_ii2str

    !> r2str: full-precision real to string (ES23.16)
    subroutine test_r2str(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, trim(r2str(1.5_dp)), "1.5000000000000000E+00")
        if (allocated(error)) return
        call check(error, trim(r2str(-0.03125_dp)), "-3.1250000000000000E-02")
    end subroutine test_r2str

    !> r2strt: truncated real to string (ES9.2)
    subroutine test_r2strt(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, trim(r2strt(1.5_dp)), "1.50E+00")
        if (allocated(error)) return
        call check(error, trim(r2strt(-12345.6789_dp)), "-1.23E+04")
    end subroutine test_r2strt

    !> c2str: complex to string, sign handling of imaginary part
    !!
    !! \note The double space after the sign is long-standing behavior: the
    !! imaginary part keeps the leading blank of its ES edit descriptor.
    subroutine test_c2str(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, trim(c2str((1.5_dp,-2.25_dp))), &
            &"1.5000000000000000E+00 -  2.2500000000000000E+00")
        if (allocated(error)) return
        call check(error, trim(c2strt((1.5_dp,2.25_dp))), &
            &"1.50E+00 +  2.25E+00 i")
    end subroutine test_c2str

    !> strh2l / strl2h: case conversion, non-letters untouched
    subroutine test_case_conversion(error)
        type(error_type), allocatable, intent(out) :: error

        call check(error, strh2l("PB3D Rocks!"), "pb3d rocks!")
        if (allocated(error)) return
        call check(error, strl2h("pb3d rocks!"), "PB3D ROCKS!")
        if (allocated(error)) return
        ! round trip
        call check(error, strl2h(strh2l("MiXeD 123")), "MIXED 123")
    end subroutine test_case_conversion

    !> merge_strings: comma-separated concatenation
    subroutine test_merge_strings(error)
        type(error_type), allocatable, intent(out) :: error

        character(len=5) :: strs(3)

        strs = ["one  ", "two  ", "three"]
        call check(error, trim(merge_strings(strs)), "one, two, three")
        if (allocated(error)) return
        call check(error, trim(merge_strings(strs(1:1))), "one")
    end subroutine test_merge_strings
end module test_str_utilities
