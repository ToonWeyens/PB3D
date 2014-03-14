! This module contains operations concerning giving output, on the screen as 
! well as in a file
module output_ops
    use netcdf
    use str_ops, only: i2str, r2str
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public init_output_ops, lvl_ud, writo, write_out, print_ar_1, print_ar_2, &
        &lvl, lvl_sep, format_out

    ! global variables
    integer :: lvl                                                              ! lvl determines the indenting. higher lvl = more indenting
    character(len=2) :: lvl_sep = ''                                            ! characters that separate different levels of output
    integer :: format_out


contains
    ! initialize the variables for the module
    subroutine init_output_ops
        lvl = 1
    end subroutine

    ! increases/decreases lvl of output
    subroutine lvl_ud(inc)                                                      ! ud from up / down
        integer :: inc
        if (lvl+inc.lt.1) then
            write(*,*) 'WARNING: cannot go below lowest level.'
            lvl = 1
        else
            lvl = lvl + inc
        end if
    end subroutine
 
    ! write output using method indicated by format_out
    subroutine write_out(nx, ny, fun, fun_name, alt_output_i, comment)
        use num_vars, only: output_i
        
        integer, intent(in) :: nx, ny
        character(len=*) :: fun_name
        real(dp) :: fun(1:nx,1:ny)
        integer, intent(in), optional :: alt_output_i
        character(len=*), optional :: comment

        integer :: fin_output_i

        if (present(alt_output_i)) then 
            fin_output_i = alt_output_i
        else
            fin_output_i = output_i
        end if
        
        ! choose correct format depending on format_out
        select case (format_out)
            case (1)                                                            ! NETCDF
                call writo('WARNING: NETCDF output not yet implemented')
            case (2)                                                            ! matlab
                call write_out_matlab(fin_output_i, nx, ny, fun, fun_name, &
                    &comment)
            case default 
                call writo('WARNING: No output format associated with ' &
                    &// i2str(format_out) )
        end select
    contains
        ! writes a 2D array into matlab format
        ! Called only by write_out!
        subroutine write_out_matlab(output_i, nx, ny, fun, fun_name, comment)
            integer, intent(in) :: output_i, nx, ny
            character(len=*) :: fun_name
            real(dp) :: fun(1:nx,1:ny)
            character(len=*), optional :: comment
            
            integer :: ix, iy
            character(len=max_str_ln) :: form_str
            
            form_str = ''
            form_str = '(A,' // trim(i2str(ny)) // '(17x,i3,16x))'
            
            if (present(comment)) then
                write(output_i,*) '% comment: ' // trim(comment)
            end if
            write(output_i,*) trim(fun_name) // ' = ['
            write(output_i,form_str) '%', (iy,iy=1,ny)
            form_str = ''
            form_str = '(' // trim(i2str(ny)) // &
                &'(es36.22e4),1x,A,1x,i3,1x)'
            do ix=1,nx
                write(output_i,form_str) (fun(ix,iy), iy = 1, ny), '%', ix 
            enddo 
            write(output_i,*) '];'
        end subroutine
    end subroutine

    ! write output to optional file number 'file_i' using the correct 
    ! indentation for the level ('lvl_loc') of the output
    subroutine writo(input_str,lvl_input,file_i)
        character(len=*), intent(in) :: input_str                               ! the name that is searched for
        integer, optional, intent(in) :: lvl_input                              ! optinally set the output lvl
        integer, optional, intent(in) :: file_i                                 ! optionally set the number of output file
        
        character(len=max_str_ln) :: output_str                                 ! the name that is searched for
        
        character(len=max_str_ln) :: header_str                                 ! the name that is searched for
        integer :: lvl_loc                                                      ! holds either lvl_loc or lvl
        integer :: id, i_part, max_len_part, num_parts, st_part, en_part
        
        ! assign either lvl_input or lvl to lvl_loc
        if (present(lvl_input)) then
            lvl_loc = lvl_input
        else
            lvl_loc = lvl
        end if
        
        ! Divide the input string length by the max_str_ln and loop over the different parts
        max_len_part = max_str_ln-(lvl_loc-1)*len(lvl_sep)                      ! max length of a part
        num_parts = (len(trim(input_str))-1)/(max_len_part) + 1                 ! how many parts there are
        do i_part = 1, num_parts
            ! construct input string for the appropriate level
            st_part = (i_part-1)*max_len_part+1                                 ! index of start of this part
            if (i_part.lt.num_parts) then                                       ! index of end of this part
                en_part = i_part*max_len_part 
            else                                                                ! last part is shorter
                en_part = len(trim(input_str))
            end if
            output_str = &
                &input_str(st_part:en_part)
            do id = 1,lvl_loc-1                                                 ! start with lvl_loc 1
                output_str = lvl_sep // trim(output_str)
            end do
            
            ! construct header string of equal length as output strength
            header_str = ''
            do id = 1, len(trim(output_str))
                header_str =  trim(header_str) // '-'
            end do
            
            if (present(file_i)) then
                if (lvl_loc.eq.1) write(file_i,*) header_str                    ! first level gets extra lines
                write(file_i,*) output_str
                if (lvl_loc.eq.1) write(file_i,*) header_str
            else
                if (lvl_loc.eq.1) write(*,*) header_str                         ! first level gets extra lines
                write(*,*) output_str
                if (lvl_loc.eq.1) write(*,*) header_str
            end if
        end do
    end subroutine

    ! print an array of dimension 2 on the screen
    subroutine print_ar_2(arr)
        real(dp) :: arr(:,:)
        
        integer :: id
        
        do id=1,size(arr,1)
            call print_ar_1(arr(id,:))
        end do
    end subroutine

    ! print an array of dimension 1 on the screen
    subroutine print_ar_1(arr)
        real(dp) :: arr(:)
        
        integer :: id
        character(len=2*max_str_ln) :: output_str                               ! holds string for a range of values
        character(len=14) :: var_str                                             ! holds string for one value, for last value
        integer :: vlen                                                         ! space for one variable + leading space
        integer :: n_free                                                       ! how much space free
        
        logical :: str_full                                                     ! whether the output_str is full
        
        output_str = '|'
        str_full = .false.
        id = 1
        vlen = 10                                                               ! 9 for variable and 1 for space
        
        do while(.not.str_full .and. id.le.size(arr))
            var_str = ''
            n_free = len(output_str) - len(trim(output_str))                    ! how many characters free
            
            if (n_free.ge.2*vlen+6) then                                        ! >2 fit
                continue
            else if (2*vlen+2.le.n_free .and. n_free.lt.2*vlen+6 .and. &        ! special case: maximum 2 left -> all fit
                &id+1.ge.size(arr)) then
                    continue
            else if (vlen+6.le.n_free .and. n_free.lt.2*vlen+6) then            ! 1 fit, rest truncated
                if (id.ge.size(arr)) then                                       ! only 1 left -> all fit
                    continue
                else                                                            ! truncate
                    str_full = .true.
                end if
            else if (vlen+2.le.n_free .and. n_free.lt.vlen+6 .and. &            ! special case: maximum 1 left -> all fit
                &id.eq.size(arr)) then
                continue
            else                                                                ! this should never be reached
                call writo('WARNING: not enough room to display variable. &
                    &Something went wrong')
                return
            end if
            
            
            if (str_full) then                                                  ! truncate
                write (var_str, '(ES9.2)') arr(size(arr))                            ! last variable
                var_str = ' ... ' // trim(var_str)
            else 
                write (var_str, '(ES9.2)') arr(id)                                   ! current variable
                var_str = ' ' // trim(var_str)
            end if
            output_str = trim(output_str) // trim(var_str)
            
            id = id+1
        end do
        output_str = trim(output_str) // ' |'
        write(*,*) output_str
    end subroutine
end module output_ops
