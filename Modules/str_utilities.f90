!------------------------------------------------------------------------------!
!> Operations on strings.
!------------------------------------------------------------------------------!
module str_utilities
#include <wrappers.h>
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public i2str, ii2str, r2str, r2strt, c2str, c2strt, strh2l, strl2h, &
        &merge_strings

contains
    !> Convert an integer to string.
    !!
    !! \see from
    !! <http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran>
    elemental character(len=max_str_ln) function i2str(k)
        ! input / output
        integer, intent(in) :: k                                                !< integer to convert
        
        write (i2str, *) k
        i2str = adjustl(i2str)
    end function i2str
    !> Convert an integer to string.
    !!
    !! Version with kind 8 integers.
    !!
    !! \see See i2str().
    elemental character(len=max_str_ln) function ii2str(k)
        ! input / output
        integer(kind=8), intent(in) :: k                                        !< integer to convert
        
        write (ii2str, *) k
        ii2str = adjustl(ii2str)
    end function ii2str

    !> Convert a real (double) to string
    !!
    !! \note See  <http://www.fortran90.org/src/best-practices.html> how  to not
    !! lose precision.
    elemental character(len=max_str_ln) function r2str(k)
        ! input / output
        real(dp), intent(in) :: k                                               !< real to convert
        
        write (r2str, '(ES23.16)') k
        r2str = adjustl(r2str)
    end function r2str
    !> Convert a real (double) to string
    !!
    !! Version with less precise output.
    !!
    !! \see See r2str().
    elemental character(len=max_str_ln) function r2strt(k)
        ! input / output
        real(dp), intent(in) :: k                                               !< real to convert
        
        write (r2strt, '(ES9.2)') k
        r2strt = adjustl(r2strt)
    end function r2strt
    
    !> Convert a complex (double) to string
    !!
    !! \note See  <http://www.fortran90.org/src/best-practices.html> how  to not
    !! lose precision.
    elemental character(len=max_str_ln) function c2str(k)
        ! input / output
        complex(dp), intent(in) :: k                                            !< complex to convert
        
        ! local variables
        character(len=max_str_ln) :: dum_str                                    ! dummy string
        
        write (c2str, '(ES23.16)') rp(k)
        write (dum_str, '(ES23.16)') abs(imag(k))
        if (imag(k).lt.0) then
            c2str = trim(c2str)//' -'
        else
            c2str = trim(c2str)//' +'
        end if
        c2str = trim(c2str)//' '//dum_str
        c2str = adjustl(c2str)
    end function c2str
    !> Convert a complex (double) to string
    !!
    !! Version with less precise output.
    !!
    !! \see See c2str().
    elemental character(len=max_str_ln) function c2strt(k)
        ! input / output
        complex(dp), intent(in) :: k                                            !< complex to convert
        
        ! local variables
        character(len=max_str_ln) :: dum_str                                    ! dummy string
        
        write (c2strt, '(ES9.2)') rp(k)
        write (dum_str, '(ES9.2)') abs(imag(k))
        if (imag(k).lt.0) then
            c2strt = trim(c2strt)//' -'
        else
            c2strt = trim(c2strt)//' +'
        end if
        c2strt = trim(c2strt)//' '//trim(dum_str)//' i'
        c2strt = adjustl(c2strt)
    end function c2strt

    !> Convert a string to lowercase.
    !!
    !! \see from \cite RedwineF90, figure 3.5B, pg 80. 
    function strh2l(input_string) result(output_string)
        ! input / output
        character(*), intent(in)     :: input_string                            !< input string
        character(len(input_string)) :: output_string                           !< lowercase version
        
        ! local variables
        integer :: i, n
        character(*), parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
        character(*), parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        
        ! copy input string
        output_string = input_string
        
        ! convert case character by character
        do i = 1, len(output_string)
          n = index(upper_case, output_string(i:i))
          if ( n /= 0 ) output_string(i:i) = lower_case(n:n)
        end do
    end function strh2l
    !> convert a string to uppercase.
    !!
    !! \see See strh2l()
    function strl2h(input_string) result(output_string)
        ! input / output
        character(*), intent(in)     :: input_string                            !< input string
        character(len(input_string)) :: output_string                           !< uppercase version
        
        ! local variables
        integer :: i, n
        character(*), parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
        character(*), parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        
        ! copy input string
        output_string = input_string
        
        ! convert case character by character
        do i = 1, len(output_string)
          n = index(lower_case, output_string(i:i))
          if ( n /= 0 ) output_string(i:i) = upper_case(n:n)
        end do
    end function strl2h
    
    !> Merge array of strings.
    function merge_strings(input_strings)
        ! input / output
        character(*), intent(in) :: input_strings(:)                            !< array of strings
        character((len(input_strings)+2)*size(input_strings)) :: merge_strings  !< merged string
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! start with first string
        if (size(input_strings).gt.0) then
            merge_strings = trim(input_strings(1))
            
            ! loop over next strings
            do id = 2,size(input_strings)
                merge_strings = trim(merge_strings)//', '//&
                    &trim(input_strings(id))
            end do
        else
            merge_strings = ''
        end if
    end function merge_strings
end module str_utilities
