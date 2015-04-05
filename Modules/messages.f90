!------------------------------------------------------------------------------!
!   Operations and variables  concerning giving output, on the  screen as well !
!   as in output files                                                         !
!------------------------------------------------------------------------------!
module messages
    use str_ops
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public init_messages, lvl_ud, writo, print_ar_1, print_ar_2, &
        &print_err_msg, init_time, start_time, stop_time, passed_time, &
        &print_hello, print_goodbye, &
        &temp_output_id, max_len_temp_output, temp_output, lvl, lvl_sep, &
        &temp_output_active, temp_output_omitted, time_sep

    ! global variables
    integer :: lvl                                                              ! lvl determines the indenting. higher lvl = more indenting
    character(len=2) :: lvl_sep = ''                                            ! characters that separate different levels of output
    character(len=10) :: time_sep = ''                                          ! defines the length of time part of output
    real(dp) :: deltat                                                          ! length of time interval
    real(dp) :: t1, t2                                                          ! end points of time interval
    logical :: running                                                          ! whether the timer is running
    logical :: temp_output_active                                               ! true if temporary output is to be written in temp_output
    integer :: temp_output_omitted                                              ! how many temporary outputs omitted
    integer :: temp_output_id                                                   ! counts the entries in temp_output
    integer, parameter :: max_len_temp_output = 70                              ! max. length of array of chars temp_output
    character(len=max_str_ln), allocatable :: temp_output(:)                    ! temporary output, before output file is opened
    
contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    subroutine init_messages
        ! output level
        lvl = 1
        
        ! time
        deltat = 0
        t1 = 0
        t2 = 0
        running = .false. 
        
        ! temporary output
        temp_output_active = .true.
        temp_output_omitted = 0
        temp_output_id = 1
        allocate(temp_output(max_len_temp_output))
    end subroutine init_messages
    
    ! prints first message
    subroutine print_hello
        use num_vars, only: glb_rank, prog_name, prog_version
        
        if (glb_rank.eq.0) then
            write(*,*) 'Simulation started on '//get_date()//&
                &', at '//get_clock()
            write(*,*) trim(prog_name)//' version: '//trim(prog_version)
            write(*,*) ''
        end if
    end subroutine print_hello

    ! prints last message
    subroutine print_goodbye
        use num_vars, only: glb_rank
        
        if (glb_rank.eq.0) then
            write(*,*) ''
            write(*,*) 'Simulation finished on '//get_date()//&
                &', at '//get_clock()
        end if
    end subroutine print_goodbye
    
    ! intialize the time passed to 0
    ! [MPI] All ranks
    subroutine init_time
        deltat = 0
        t1 = 0
        t2 = 0
        running = .false. 
    end subroutine

    ! start a timer
    ! [MPI] All ranks
    subroutine start_time
        if (running) then
            call writo('WARNING: Tried to start timer, but was already running')
        else
            call second0(t1)
            running = .true.
        end if
    end subroutine

    ! stop a timer
    subroutine stop_time
        if (running) then
            call second0(t2)
            
            ! increase deltat
            deltat = deltat+t2-t1

            ! set t1 and t2 back to zero
            t1 = 0
            t2 = 0
            running = .false.
        else
            call writo('WARNING: Tried to stop timer, but was already stopped')
        end if
    end subroutine

    ! display the time that has passed between t1 and t2
    ! automatically stops time and resets everything to zero
    subroutine passed_time
        ! local variables
        character(len=max_str_ln) :: begin_str, end_str
        integer :: time_in_units(4)                                             ! time in days, hours, minutes, seconds
        integer :: total_time
        character(len=6) :: time_units(4) = ['day   ','hour  ','minute',&
            &'second']
        integer :: id                                                           ! counter

        ! stop at current time if running
        if (running) call stop_time

        begin_str = '(this took'
        if (deltat.lt.1) then
            end_str = ' less than 1 second)'
        else
            ! get days, hours, minutes, seconds
            total_time = floor(deltat)                                          ! get total time in seconds
            time_in_units(1) = total_time/(60*60*24)                            ! get days
            total_time = mod(total_time,60*60*24)                               ! subtract days from total time
            time_in_units(2) = total_time/(60*60)                               ! get hours
            total_time = mod(total_time,60*60)                                  ! subtract hours from total time
            time_in_units(3) = total_time/(60)                                  ! get minutes
            total_time = mod(total_time,60)                                     ! subtract minutes from total time
            time_in_units(4) = total_time                                       ! get seconds
            
            ! set up end_str
            end_str = ''
            do id = 1,4
                if (time_in_units(id).gt.0) then
                    if (trim(end_str).ne.'') end_str = trim(end_str)//','
                    end_str = trim(end_str)//' '//&
                        &trim(i2str(time_in_units(id)))//' '//&
                        &trim(time_units(id))
                    if (time_in_units(id).gt.1) end_str = trim(end_str)//'s'
                end if
            end do
            end_str = trim(end_str)//')'
        end if
        call writo(trim(begin_str) // trim(end_str))
        
        ! restart deltat
        call init_time
    end subroutine
    
    ! returns the date
    ! (from http://infohost.nmt.edu/tcc/help/lang/fortran/date.html)
    function get_date() result(date)
        ! input / output
        character(len=10) :: date                                               ! date
        
        ! local variables
        integer :: today(3)
        
        call idate(today)                                                       ! today(1)=day, (2)=month, (3)=year
        
        write (date,'(i2.2,"/",i2.2,"/",i4.4)')  today(2), today(1), today(3)
    end function get_date
    
    ! returns the time
    ! (from http://infohost.nmt.edu/tcc/help/lang/fortran/date.html)
    function get_clock() result(time)
        ! input / output
        character(len=8) :: time                                                ! time
        
        ! local variables
        integer :: now(3)
        
        call itime(now)                                                         ! now(1)=hour, (2)=minute, (3)=second
        
        write (time,'(i2.2,":",i2.2,":",i2.2)')  now
    end function get_clock
    
    ! prints an error  message that is either user-provided, or  the name of the
    ! calling routine
    subroutine print_err_msg(err_msg,routine_name)
        use num_vars, only: glb_rank
        character(len=*) :: err_msg, routine_name
        
        if (trim(err_msg).eq.'') then
            lvl = 2
            call writo('>> calling routine: '//trim(routine_name)//' of rank '&
                &//trim(i2str(glb_rank)),persistent=.true.)
        else
            call writo('ERROR in '//trim(routine_name)//': '//trim(err_msg),&
                &persistent=.true.)
        end if
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
    
    ! write output to file identified by output_i, using the correct indentation
    ! for the level ('lvl_loc') of the output. If first alpha group, also output
    ! to screen
    ! [MPI] Only masters of groups of alpha and the global master call these
    !       The  global master  outputs to  the  master output  file, while  the
    !       masters  of the  groups of  alpha  write their  output to  different
    !       files, which  are then read  by the  global master when  the group's
    !       work is done
    !       Before the groups are created, the group rank is set to be identical
    !       to  the global  rank. This  way,  the group  rank always  determines
    !       whether a  process outputs  or not,  also when  there are  no groups
    !       (yet)
    subroutine writo(input_str,persistent)
        use num_vars, only: grp_rank, output_i, no_messages
        
        ! input / output
        character(len=*), intent(in) :: input_str                               ! the name that is searched for
        logical, optional, intent(in) :: persistent                             ! output even if not group master
        
        ! local variables
        character(len=max_str_ln) :: output_str                                 ! the name that is searched for
        character(len=max_str_ln) :: header_str                                 ! the name that is searched for
        character(len=max_str_ln) :: temp_output_err_str(2)                     ! error if temp_output is full
        integer :: id, i_part                                                   ! counters
        integer :: max_len_part, num_parts, st_part, en_part                    ! variables controlling strings
        logical :: ignore                                                       ! normally, everybody but group master is ignored
        
        ! bypass messages if no_messages
        if (no_messages) return
        
        ! setup ignore
        ignore = .true.                                                         ! ignore by default
        if (grp_rank.eq.0) ignore = .false.                                     ! group masters can output
        if (present(persistent)) ignore = .not.persistent                       ! persistent can override this
        
        if (.not.ignore) then                                                   ! only group master (= global master if no groups) or persistent
            ! Divide the input string length by the max_str_ln and loop over the
            ! different parts
            max_len_part = max_str_ln-(lvl-1)*len(lvl_sep) - len(time_sep)      ! max length of a part including time part
            num_parts = (len(trim(input_str))-1)/(max_len_part) + 1             ! how many parts there are
            do i_part = 1, num_parts
                ! construct input string for the appropriate level
                st_part = (i_part-1)*max_len_part+1                             ! index of start of this part
                if (i_part.lt.num_parts) then                                   ! index of end of this part
                    en_part = i_part*max_len_part 
                else                                                            ! last part is shorter
                    en_part = len(trim(input_str))
                end if
                output_str = input_str(st_part:en_part)
                call format_str(lvl,output_str)
                
                ! construct header string of equal length as output strength
                header_str = ''
                do id = 1, len(trim(output_str)) - len(time_sep)                ! including time part
                    header_str =  trim(header_str) // '-'
                end do
                header_str = '          '//trim(header_str)
                
                ! error string for temporary output
                temp_output_err_str(1) = 'WARNING: not enough space in &
                    &temp_output, start omiting output from log file'           ! error message
                temp_output_err_str(2) = 'Increase the variable &
                    &"max_len_temp_output" in output_ops'                       ! error message
                
                ! write output to file output_i or to temporary output
                if (temp_output_active) then                                    ! temporary output to internal variable temp_output
                    if (temp_output_id.le.max_len_temp_output-3) then           ! still room for at least 3 lines of temporary output
                        if (lvl.eq.1) then                                      ! first level
                            temp_output(temp_output_id) = header_str            ! first level gets extra lines
                            temp_output_id = temp_output_id + 1
                            temp_output(temp_output_id) = output_str
                            temp_output_id = temp_output_id + 1
                            temp_output(temp_output_id) = header_str
                            temp_output_id = temp_output_id + 1                 ! increase counter by 3
                        else                                                    ! other levels only need one line (but error needs 2)
                            temp_output(temp_output_id) = output_str
                            temp_output_id = temp_output_id + 1                 ! increase counter by 1
                        end if
                    else                                                        ! maximum 2 lines left for temporary output: for error
                        if (temp_output_omitted.eq.0) then                      ! up to now success: first error
                            ! error message
                            do id = 1,2
                                call format_str(lvl,temp_output_err_str(id))    ! format error message
                                temp_output(temp_output_id) = &
                                    &temp_output_err_str(id)                    ! save error message as final temporary output
                                write(*,*) temp_output_err_str(id)              ! write error message on screen
                                temp_output_id = temp_output_id + 1             ! 2 last entries in log file
                            end do
                        end if
                        if (lvl.eq.1) then                                      ! first level
                            temp_output_omitted = temp_output_omitted + 3       ! 3 output lines omitted
                        else                                                    ! other levels
                            temp_output_omitted = temp_output_omitted + 1       ! 1 output lines omitted
                        end if
                    end if
                else                                                            ! normal output to file output_i
                    if (lvl.eq.1) write(output_i,*) header_str                  ! first level gets extra lines
                    write(output_i,*) output_str
                    if (lvl.eq.1) write(output_i,*) header_str
                end if
                
                ! also write output to screen
                if (lvl.eq.1) write(*,*) header_str                             ! first level gets extra lines
                write(*,*) output_str
                if (lvl.eq.1) write(*,*) header_str
            end do
        end if
    contains
        subroutine format_str(lvl,str)
            ! input / output
            character(len=*), intent(inout) :: str                              ! string that is to be formatted for the lvl
            integer, intent(in) :: lvl                                          ! lvl of indentation
            
            ! local variables
            integer :: id                                                       ! counter
            
            do id = 1,lvl-1                                                     ! start with lvl 1
                str = lvl_sep // trim(str)
            end do
            str = get_clock()//': '//trim(str)
        end subroutine format_str
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
end module messages
