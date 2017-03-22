!------------------------------------------------------------------------------!
!   Numerical utilities related to giving output                               !
!------------------------------------------------------------------------------!
module messages
    use str_utilities
    use num_vars, only: dp, max_str_ln
    use foul
    
    implicit none
    private
    public init_output, lvl_ud, writo, print_ar_1, print_ar_2, &
        &print_err_msg, init_time, start_time, stop_time, passed_time, &
        &print_hello, print_goodbye, &
        &temp_output, lvl, lvl_sep, temp_output_active, time_sep
#if ldebug
    public get_mem_usage
#endif

    ! global variables
    integer :: lvl                                                              ! lvl determines the indenting. higher lvl = more indenting
    character(len=2) :: lvl_sep = ''                                            ! characters that separate different levels of output
    character(len=10) :: time_sep = ''                                          ! defines the length of time part of output
    real(dp) :: deltat                                                          ! length of time interval
    real(dp) :: t1, t2                                                          ! end points of time interval
    logical :: running                                                          ! whether the timer is running
    logical :: temp_output_active                                               ! true if temporary output is to be written in temp_output
    character(len=max_str_ln), allocatable :: temp_output(:)                    ! temporary output, before output file is opened
    
contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    subroutine init_output
        use num_vars, only: rank
#if ldebug
        use num_vars, only: mem_usage_count
#endif
        
        ! output level
        lvl = 1
        
        ! time
        deltat = 0
        t1 = 0
        t2 = 0
        running = .false. 
        
        ! temporary output for master
        if (rank.eq.0) then
            temp_output_active = .true.
            allocate(temp_output(0))
        else
            temp_output_active = .false.
        end if
        
#if ldebug
        mem_usage_count = 0
#endif
    end subroutine init_output
    
    ! prints first message
    subroutine print_hello
        use num_vars, only: rank, prog_name, prog_version, n_procs
        
        ! local variables
        integer :: istat                                                        ! status
        
        if (rank.eq.0) then
            call write_formatted(' Simulation started on '//get_date()//', at '&
                &//get_clock(),'italic')
            call write_formatted(' '//prog_name//' version: '//&
                &trim(r2strt(prog_version)),'italic')
            if (n_procs.eq.1) then
                call write_formatted(' 1 MPI process','italic')
            else
                call write_formatted(' '//trim(i2str(n_procs))//&
                    &' MPI processes','italic')
            end if
            write(*,*,IOSTAT=istat) ''
        end if
    end subroutine print_hello

    ! prints last message
    subroutine print_goodbye
        use num_vars, only: rank
        
        ! local variables
        integer :: istat                                                        ! status
        
        if (rank.eq.0) then
            write(*,*,IOSTAT=istat) ''
            call write_formatted(' Simulation finished on '//get_date()//&
                &', at '//get_clock(),'italic')
        end if
    end subroutine print_goodbye
    
    ! intialize the time passed to 0
    ! [MPI] All ranks
    subroutine init_time
        deltat = 0
        t1 = 0
        t2 = 0
        running = .false. 
    end subroutine init_time

    ! start a timer
    ! [MPI] All ranks
    subroutine start_time
        if (running) then
            call writo('Tried to start timer, but was already running',&
                &warning=.true.)
        else
            call cpu_time(t1)
            running = .true.
        end if
    end subroutine start_time

    ! stop a timer
    subroutine stop_time
        if (running) then
            call cpu_time(t2)
            
            ! increase deltat
            deltat = deltat+t2-t1
            
            ! set t1 and t2 back to zero
            t1 = 0
            t2 = 0
            running = .false.
        else
            call writo('Tried to stop timer, but was already stopped',&
                &warning=.true.)
        end if
    end subroutine stop_time

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
    end subroutine passed_time
    
    ! returns the date
    ! (from http://infohost.nmt.edu/tcc/help/lang/fortran/date.html)
    function get_date() result(now)
#if (lwith_intel && !lwith_gnu)
        use IFPORT
#endif
        ! input / output
        character(len=10) :: now                                                ! date
        
        ! local variables
        integer :: today(3)
        
        call idate(today)                                                       ! today(1)=day, (2)=month, (3)=year
        
        write (now,'(i2.2,"/",i2.2,"/",i4.4)')  today(2), today(1), today(3)
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
        use num_vars, only: rank
        
        ! input / output
        character(len=*), intent(in) :: err_msg, routine_name
        
        if (trim(err_msg).eq.'') then
            lvl = 2
            call writo('>> calling routine: '//trim(routine_name)//' of rank '&
                &//trim(i2str(rank)),persistent=.true.)
        else
            call writo('ERROR in '//trim(routine_name)//': '//trim(err_msg),&
                &persistent=.true.,error=.true.)
        end if
    end subroutine print_err_msg
    
    ! increases/decreases lvl of output
    subroutine lvl_ud(inc)                                                      ! ud from up / down
        integer :: inc
        if (lvl+inc.lt.1) then
            lvl = 1
            call writo('cannot go below lowest level',warning=.true.)
        else
            lvl = lvl + inc
        end if
    end subroutine lvl_ud
    
    ! write output to file identified by output_i, using the correct indentation
    ! for the level ('lvl_loc') of the output. If first alpha group, also output
    ! to screen
    ! Optionally, special formatting for error, warning or alert can be chosen.
    ! [MPI] Only masters of groups of alpha and the global master call these
    !       The  global master  outputs to  the  master output  file, while  the
    !       masters  of the  groups of  alpha  write their  output to  different
    !       files, which  are then read  by the  global master when  the group's
    !       work is done
    subroutine writo(input_str,persistent,error,warning,alert)
        use num_vars, only: rank, output_i, no_output, max_tot_mem, max_X_mem
#if ldebug
        use MPI
        use num_vars, only: print_mem_usage, prog_name, mem_usage_name, &
            &mem_usage_i, mem_usage_count, time_start
#endif
        
        ! input / output
        character(len=*), intent(in) :: input_str                               ! the name that is searched for
        logical, intent(in), optional :: persistent                             ! output even if not group master
        logical, intent(in), optional :: error                                  ! error message
        logical, intent(in), optional :: warning                                ! warning message
        logical, intent(in), optional :: alert                                  ! alert message
        
        ! local variables
        character(len=max_str_ln) :: output_str                                 ! output string
        character(len=max_str_ln) :: time_str                                   ! time string
        character(len=max_str_ln) :: header_str                                 ! header string
        character(len=max_str_ln) :: input_str_loc                              ! local input string
        character(len=max_str_ln), allocatable :: temp_output_loc(:)            ! local temporary output
        integer :: id, i_part                                                   ! counters
        integer :: max_len_part, num_parts, st_part, en_part                    ! variables controlling strings
        logical :: ignore                                                       ! normally, everybody but group master is ignored
        logical :: error_loc                                                    ! local error
        logical :: warning_loc                                                  ! local warning
        logical :: alert_loc                                                    ! local alert
        integer :: istat                                                        ! status
#if ldebug
        integer :: mem_usage                                                    ! memory usage
        integer(kind=8) :: clock                                                ! current clock
#endif
        
        ! bypass output if no_output
        if (no_output) return
        
        ! setup ignore, error, warning and alert
        ignore = .true.                                                         ! ignore by default
        if (rank.eq.0) ignore = .false.                                         ! master process can output
        if (present(persistent)) ignore = .not.persistent                       ! persistent can override this
        error_loc = .false.
        if (present(error)) error_loc = error
        warning_loc = .false.
        if (present(warning)) warning_loc = warning
        alert_loc = .false.
        if (present(alert)) alert_loc = alert
        
        ! set local input string
        input_str_loc = input_str
        
        ! prepend "WARNING: " if warning
        if (warning_loc) input_str_loc = 'WARNING: '//trim(input_str_loc)
        
#if ldebug
        ! memory usage
        if (print_mem_usage) then
            ! increment counter
            mem_usage_count = mem_usage_count + 1
            
            ! get memory usage
            mem_usage = get_mem_usage()
            
            ! append count and memory usage in MegaBytes
            input_str_loc = trim(input_str_loc)//' - ['//&
                &trim(i2str(mem_usage_count))//': '//&
                &trim(i2str(mem_usage))//'kB]'
            
            ! get clock
            call system_clock(clock)
            
            ! write rank, count, time, memory usage to file if not temp_output
            if (.not.temp_output_active) then
                open(UNIT=mem_usage_i,FILE=prog_name//'_'//&
                    &trim(mem_usage_name)//'.dat',STATUS='old',&
                    &POSITION='append',IOSTAT=istat)
                if (istat.eq.0) then
                    write(mem_usage_i,"(1X,2I10,I21,I10,2ES23.16)",&
                        &IOSTAT=istat) &
                        &rank, mem_usage_count, clock-time_start, mem_usage, &
                        &max_tot_mem*1000, max_X_mem*1000
                    close(UNIT=mem_usage_i,IOSTAT=istat)
                end if
            end if
        end if
#endif
        
        if (.not.ignore) then                                                   ! only group master (= global master if no groups) or persistent
            ! set local error
            
            ! Divide the input string length by the max_str_ln and loop over the
            ! different parts
            max_len_part = max_str_ln-(lvl-1)*len(lvl_sep) - len(time_sep)      ! max length of a part including time part
            num_parts = (len(trim(input_str_loc))-1)/(max_len_part) + 1         ! how many parts there are
            do i_part = 1, num_parts
                ! construct input string for the appropriate level
                st_part = (i_part-1)*max_len_part+1                             ! index of start of this part
                if (i_part.lt.num_parts) then                                   ! index of end of this part
                    en_part = i_part*max_len_part 
                else                                                            ! last part is shorter
                    en_part = len(trim(input_str_loc))
                end if
                output_str = input_str_loc(st_part:en_part)
                call format_str(lvl,output_str)
                call get_time_str(time_str)
                
                ! construct header string of equal length as output strength
                header_str = ''
                do id = 1, len(trim(output_str)) + len(trim(time_str)) + 1 - &
                    &len(time_sep)                                              ! not including time part, 1 for the space
                    header_str =  trim(header_str) // '-'
                end do
                header_str =  '          '//trim(header_str)                    ! number of spaces matches time string
                
                ! write output to file output_i or to temporary output
                if (temp_output_active) then                                    ! temporary output to internal variable temp_output
                    ! back up previous temporary output in local variable
                    allocate(temp_output_loc(size(temp_output)))
                    temp_output_loc = temp_output
                    
                    ! reallocate temp_output
                    deallocate(temp_output)
                    if (lvl.eq.1) then                                          ! first level
                        allocate(temp_output(size(temp_output_loc)+3))
                        temp_output(1:size(temp_output_loc)) = temp_output_loc
                        temp_output(size(temp_output_loc)+1) = trim(header_str) ! first level gets extra lines
                        temp_output(size(temp_output_loc)+2) = trim(time_str)//&
                            &' '//trim(output_str)
                        temp_output(size(temp_output_loc)+3) = trim(header_str)
                    else                                                        ! other levels only need one line
                        allocate(temp_output(size(temp_output_loc)+1))
                        temp_output(1:size(temp_output_loc)) = temp_output_loc
                        temp_output(size(temp_output_loc)+1) = trim(time_str)//&
                            &' '//output_str
                    end if
                    
                    ! deallocate local variable
                    deallocate(temp_output_loc)
                else                                                            ! normal output to file output_i
                    if (output_i.ne.0) then
                        if (lvl.eq.1) write(output_i,"(1X,A)",IOSTAT=istat) &
                            &trim(header_str)                                   ! first level gets extra lines
                        write(output_i,"(1X,A)",IOSTAT=istat) &
                            &trim(time_str)//' '//trim(output_str)
                        if (lvl.eq.1) write(output_i,"(1X,A)",IOSTAT=istat) &
                            &trim(header_str)
                    end if
                end if
                
                ! also write output to screen
                if (error_loc) then
                    call write_formatted(' '//trim(time_str)//' ',&
                        &'background_white',trim(output_str),'italic underline')
                else if (warning_loc .or. alert_loc) then
                    call write_formatted(' '//trim(time_str)//' ',&
                        &'background_white',trim(output_str),'')
                else if (lvl.eq.1) then
                    call write_formatted(' '//trim(header_str),&
                        &'bright red')
                    call write_formatted(' '//trim(time_str)//' ','',&
                        &trim(output_str),'bright red')
                    call write_formatted(' '//trim(header_str),&
                        &'bright red')
                else if (lvl.eq.2) then
                    call write_formatted(' '//trim(time_str)//' ','',&
                        &trim(output_str),'bright blue')
                else
                    write(*,"(1X,A)",IOSTAT=istat) &
                        &trim(time_str)//' '//trim(output_str)
                end if
            end do
        end if
    contains
        ! formats the string, merging time information with the string.
        subroutine format_str(lvl,str)
            ! input / output
            integer, intent(in) :: lvl                                          ! lvl of indentation
            character(len=*), intent(inout) :: str                              ! string that is to be formatted for the lvl
            
            ! local variables
            integer :: id                                                       ! counter
            
            do id = 1,lvl-1                                                     ! start with lvl 1
                str = lvl_sep // trim(str)
            end do
        end subroutine format_str
        
        ! prints a string with time information
        subroutine get_time_str(time_str)
            ! input / output
            character(len=*), intent(inout) :: time_str                         ! string with time information
            
            time_str = get_clock()//':'
        end subroutine get_time_str
    end subroutine writo

    ! print an array of dimension 2 on the screen
    subroutine print_ar_2(arr)
        real(dp) :: arr(:,:)
        
        integer :: id
        
        do id=1,size(arr,1)
            call print_ar_1(arr(id,:))
        end do
    end subroutine print_ar_2

    ! print an array of dimension 1 on the screen
    subroutine print_ar_1(arr)
        real(dp) :: arr(:)
        
        ! local variables
        integer :: id
        character(len=2*max_str_ln) :: output_str                               ! holds string for a range of values
        character(len=14) :: var_str                                             ! holds string for one value, for last value
        integer :: vlen                                                         ! space for one variable + leading space
        integer :: n_free                                                       ! how much space free
        logical :: str_full                                                     ! whether the output_str is full
        integer :: istat                                                        ! status
        
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
                call writo('not enough room to display variable. Something &
                    &went wrong',warning=.true.)
                return
            end if
            
            
            if (str_full) then                                                  ! truncate
                write (var_str, '(ES9.2)') arr(size(arr))                       ! last variable
                var_str = ' ... ' // trim(var_str)
            else 
                write (var_str, '(ES9.2)') arr(id)                              ! current variable
                var_str = ' ' // trim(var_str)
            end if
            output_str = trim(output_str) // trim(var_str)
            
            id = id+1
        end do
        output_str = trim(output_str) // ' |'
        write(*,'(1X,A)',IOSTAT=istat) output_str
    end subroutine print_ar_1
    
#if ldebug
    ! Returns the memory usage in kilobytes.
    ! (based on http://stackoverflow.com/questions/22028571/
    !  track-memory-usage-in-fortran-90)
    ! Note: only works under linux
    integer function get_mem_usage() result(mem)
        use num_vars, only: mem_usage_i
#if ( lwith_intel && !lwith_gnu)
        use IFPORT
#endif
        
        ! local variables
        character(len=200):: filename=' '                                       ! name of file where memory stored
        character(len=80) :: line                                               ! line of memory file
        character(len=8)  :: pid_char=' '                                       ! process ID in string
        logical :: exists                                                       ! whether memory file exists
        integer :: pid                                                          ! process ID
        integer :: istat                                                        ! status
        
        ! initiazlie mem to negative number
        mem = -1
        
        ! get process ID and write it in string
        pid = getpid()
        write(pid_char,'(I8)') pid
        
        ! set up file name
        filename='/proc/'//trim(adjustl(pid_char))//'/status'
        
        ! inquire about memory file
        inquire (file=filename,exist=exists)
        if (.not.exists) return
        
        ! open and read memory file
        open(UNIT=mem_usage_i-1,FILE=filename,ACTION='read',IOSTAT=istat)
        if (istat.ne.0) return
        do
            read (mem_usage_i-1,'(A)',IOSTAT=istat) line
            if (istat.ne.0) return
            if (line(1:6).eq.'VmRSS:') then
                read (line(7:),*) mem
                exit
            endif
        end do
        
        ! close memory file
        close(mem_usage_i-1)
    end function get_mem_usage
#endif
end module messages
