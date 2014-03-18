!-------------------------------------------------------
!   This module contains operations concerning giving input
!-------------------------------------------------------
module input_ops
    use str_ops, only: strh2l
    use num_vars, only: max_str_ln
    implicit none
    private
    public yes_no

contains
    ! queries for yes or no, depending on the flag yes:
    !   yes = .true.: yes is default answer
    !   yes = .false.: no is default answer
    logical function yes_no(yes)
        use output_ops, only: &
            &lvl_sep, lvl
        use time, only: start_time, stop_time
        
        ! input / output
        character(len=max_str_ln) :: answer_str
        logical :: yes
        
        ! local variables
        integer :: id 

        do id = 1,lvl                                                           ! to get the response on the right collumn
            write(*,'(A)',advance='no') lvl_sep
        end do
        if (yes) then
            write(*,'(A)',advance='no') 'y(es)/n(o) [yes]: '
            yes_no = .true.
        else
            write(*,'(A)',advance='no') 'y(es)/n(o) [no]: '
            yes_no = .false.
        end if
        call stop_time
        read (*, '(A)') answer_str
        call start_time

        select case (strh2l(trim(answer_str)))
            !case ('y','Y','yes','Yes','YEs','YES','yEs','yES','yeS','YeS')
            case ('y','yes')
                yes_no = .true.
            case ('n','N','no','No','NO','nO') 
                yes_no = .false.
            case default 
        end select
    end function yes_no
end module input_ops
