!------------------------------------------------------------------------------!
!   This module contains operations on strings                                 !
!------------------------------------------------------------------------------!
module str_ops
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public i2str, r2str, r2strt, strh2l

contains
    ! Convert an integer to string 
    ! (from http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran)
    character(len=max_str_ln) function i2str(k)
        integer, intent(in) :: k
        write (i2str, *) k
        i2str = adjustl(i2str)
    end function i2str

    ! Convert a real (double) to string
    character(len=max_str_ln) function r2str(k)
        real(dp), intent(in) :: k
        write (r2str, *) k
        r2str = adjustl(r2str)
    end function r2str

    ! Convert a real (double) to string and truncate it
    character(len=max_str_ln) function r2strt(k)
        real(dp), intent(in) :: k
        write (r2strt, '(ES9.2)') k
        r2strt = adjustl(r2strt)
    end function r2strt

    ! Convert a string to lowercase 
    ! (from Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
    ! 1995 Springer-Verlag, New York. ) 
    function strh2l(input_string) result(output_string)
        character(*), intent(in)     :: input_string
        character(len(input_string)) :: output_string
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
end module str_ops
