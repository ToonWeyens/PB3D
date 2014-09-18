!------------------------------------------------------------------------------!
!   This module contains operations concerning giving output, on the screen as !
!   well as in output files                                                    !
!------------------------------------------------------------------------------!
module output_ops
    use str_ops, only: i2str, r2strt, r2str
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public init_output_ops, lvl_ud, writo, print_GP_2D, print_GP_3D, &
        &print_ar_1, print_ar_2, draw_GP, print_err_msg, init_time, &
        &start_time, stop_time, passed_time, print_hello, print_goodbye, &
        &draw_GP_animated, temp_output_id, max_len_temp_output, temp_output, &
        &lvl, lvl_sep, format_out, no_plots, temp_output_active, &
        &temp_output_omitted

    ! global variables
    integer :: lvl                                                              ! lvl determines the indenting. higher lvl = more indenting
    character(len=2) :: lvl_sep = ''                                            ! characters that separate different levels of output
    integer :: format_out                                                       ! format of output
    character(len=5) :: plot_dir = 'Plots'                                      ! directory where to save plots
    character(len=7) :: script_dir = 'Scripts'                                  ! directory where to save scripts for plots
    character(len=4) :: data_dir = 'Data'                                       ! directory where to save data for plots
    real(dp) :: deltat                                                          ! length of time interval
    real(dp) :: t1, t2                                                          ! end points of time interval
    logical :: running                                                          ! whether the timer is running
    logical :: no_plots = .false.                                               ! true if no plots should be made
    logical :: temp_output_active                                               ! true if temporary output is to be written in temp_output
    integer :: temp_output_omitted                                              ! how many temporary outputs omitted
    integer :: temp_output_id                                                   ! counts the entries in temp_output
    integer, parameter :: max_len_temp_output = 70                              ! max. length of array of chars temp_output
    character(len=max_str_ln), allocatable :: temp_output(:)                    ! temporary output, before output file is opened

    ! interfaces
    interface print_GP_2D
        module procedure print_GP_2D_ind, print_GP_2D_arr
    end interface
    interface print_GP_3D
        module procedure print_GP_3D_ind, print_GP_3D_arr
    end interface
    interface draw_GP
        module procedure draw_GP_ind, draw_GP_arr
    end interface
    interface draw_GP_animated
        module procedure draw_GP_animated_ind, draw_GP_animated_arr
    end interface
    
contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    subroutine init_output_ops
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
    end subroutine
    
    ! prints first message
    subroutine print_hello
        use num_vars, only: glb_rank
        
        if (glb_rank.eq.0) then
            write(*,*) 'Simulation started on '//get_date()//&
                &', at '//get_clock()
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
                &//trim(i2str(glb_rank)))
        else
            call writo('ERROR in '//trim(routine_name)//': '//trim(err_msg))
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

    ! print 2D output on a file
    ! The variables fun_name and file_name hold the  name of the plot and of the
    ! file in  which the plot  data is  to be saved,  respectively. y is  the an
    ! array containing the function which is  stored and x is an optional vector
    ! with the  x-values. The  logical draw=.false. optionally  disables calling
    ! the  GNUPlot drawing  procedure for  output on  screen [default],  without
    ! modifying the plot file
    ! The first index of y (and x) contains the points of a current plot
    ! The second index indicates various plots (one or more)
    subroutine print_GP_2D_ind(fun_name,file_name,y,x,draw)                     ! individual plot
        ! input / output
        character(len=*), intent(in) :: fun_name
        character(len=*), intent(in) :: file_name
        real(dp), intent(in) :: y(1:)
        real(dp), intent(in), optional :: x(1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: npoints
        
        ! set npoints, ny
        npoints = size(y)
        
        ! call multiplot version
        if (present(draw)) then
            if (present(x)) then
                call print_GP_2D_arr(fun_name,file_name,reshape(y,[npoints,1]),&
                    &x=reshape(x,[npoints,1]),draw=draw)
            else
                call print_GP_2D_arr(fun_name,file_name,reshape(y,[npoints,1]),&
                    &draw=draw)
            end if
        else
            if (present(x)) then
                call print_GP_2D_arr(fun_name,file_name,reshape(y,[npoints,1]),&
                    &x=reshape(x,[npoints,1]))
            else
                call print_GP_2D_arr(fun_name,file_name,reshape(y,[npoints,1]))
            end if
        end if
    end subroutine print_GP_2D_ind
    subroutine print_GP_2D_arr(fun_name,file_name_i,y,x,draw)                   ! multiple plots
        use num_vars, only: n_seq_0, grp_rank
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: fun_name
        character(len=*), intent(in) :: file_name_i
        real(dp), intent(in) :: y(1:,1:)
        real(dp), intent(in), optional :: x(1:,1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: file_i
        integer :: iplt, ipnt, nplt, npnt
        real(dp), allocatable :: x_fin(:,:)
        integer :: istat
        character(len=max_str_ln) :: file_name
        
        ! return if no_plots
        if (no_plots) then
            call writo('WARNING: plot ignored because no_plots is on')
            return
        end if
        
        ! set nplt, npnt
        npnt = size(y,1)
        nplt = size(y,2)
        if (present(x)) then
            if (size(x,1).ne.size(y,1) .or. size(x,2).ne.size(y,2)) then
                call writo('WARNING: In print_GP_2D, the size of x and y has &
                    &to be the same... Skipping plot')
                return
            end if
        end if
        
        ! set x
        if (present (x)) then
            x_fin = x
        else
            allocate(x_fin(size(y,1),size(y,2)))
            do ipnt = 1,npnt
                x_fin(ipnt,:) = ipnt
            end do
        end if
        
        ! set default file name if empty
        if (trim(file_name_i).eq.'') then
            file_name = 'temp_data_print_GP_2D_'//trim(i2str(grp_rank))//'.dat'
        else
            file_name = trim(file_name_i)
        end if
        
        ! open output file
        file_i = n_seq_0
        call safe_open(file_i,istat,data_dir//'/'//trim(file_name),'replace',&
            &'formatted',delim_in='none')
        
        ! write to output file
        write(file_i,*) '# '//trim(fun_name)//':'
        do ipnt = 1,npnt
            write(file_i,*) (x_fin(ipnt,iplt), iplt = 1,nplt), &
                &(y(ipnt,iplt), iplt = 1,nplt)
        enddo 
        write(file_i,*) ''
        
        ! close output file
        close(file_i)
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(fun_name,trim(file_name),nplt,.true.,.true.)
        
        if (trim(file_name_i).eq.'') then
            call execute_command_line('rm '//data_dir//'/'//trim(file_name))
        end if
    end subroutine print_GP_2D_arr
    
    ! print 3D output on a file or on the screen, using GNUPlot
    ! The variables fun_name and file_name hold the  name of the plot and of the
    ! file in  which the plot  data is  to be saved,  respectively. z is  the an
    ! array containing  the function which  is stored and  x and y  are optional
    ! vectors  with the  x  and y-values.  The  logical draw=.false.  optionally
    ! disables  calling  the GNUPlot  drawing  procedure  for output  on  screen
    ! [default], without modifying the plot file
    ! The first index of z (and x, y) contains the points of a current plot
    ! The second index indicates various plots (one or more)
    subroutine print_GP_3D_ind(fun_name,file_name,z,y,x,draw)                   ! individual plot
        ! input / output
        character(len=*), intent(in) :: fun_name
        character(len=*), intent(in) :: file_name
        real(dp), intent(in) :: z(1:,1:)
        real(dp), intent(in), optional :: y(1:,1:)
        real(dp), intent(in), optional :: x(1:,1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: npoints(2)
        
        ! set npoints, ny
        npoints = [size(z,1),size(z,2)]
        
        ! call multiplot version
        if (present(draw)) then
            if (present(y)) then
                if (present(x)) then
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),y=reshape(y,[npoints,1]),&
                        &x=reshape(x,[npoints,1]),draw=draw)
                else
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),y=reshape(y,[npoints,1]),&
                        &draw=draw)
                end if
            else
                if (present(x)) then
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),x=reshape(x,[npoints,1]),&
                        &draw=draw)
                else
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),draw=draw)
                end if
            end if
        else
            if (present(y)) then
                if (present(x)) then
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),y=reshape(y,[npoints,1]),&
                        &x=reshape(x,[npoints,1]))
                else
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),y=reshape(y,[npoints,1]))
                end if
            else
                if (present(x)) then
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]),x=reshape(x,[npoints,1]))
                else
                    call print_GP_3D_arr(fun_name,file_name,&
                        &reshape(z,[npoints,1]))
                end if
            end if
        end if
    end subroutine print_GP_3D_ind
    subroutine print_GP_3D_arr(fun_name,file_name_i,z,y,x,draw)                 ! multiple plots
        use num_vars, only: n_seq_0, grp_rank
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: fun_name
        character(len=*), intent(in) :: file_name_i
        real(dp), intent(in) :: z(1:,1:,1:)
        real(dp), intent(in), optional :: y(1:,1:,1:)
        real(dp), intent(in), optional :: x(1:,1:,1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: file_i
        integer :: iplt, ipntx, ipnty, nplt, npntx, npnty
        real(dp), allocatable :: x_fin(:,:,:)
        real(dp), allocatable :: y_fin(:,:,:)
        integer :: istat
        character(len=max_str_ln) :: file_name
        
        ! return if no_plots
        if (no_plots) then
            call writo('WARNING: plot ignored because no_plots is on')
            return
        end if
        
        ! set nplt, npnt
        npntx = size(z,1)
        npnty = size(z,2)
        nplt = size(z,3)
        if (present(x)) then
            if (size(x,1).ne.size(z,1) .or. size(x,2).ne.size(z,2) .or. &
                &size(x,3).ne.size(z,3)) then
                call writo('WARNING: In print_GP, the size of x and z has to &
                    &be the same... Skipping plot')
                return
            end if
        end if
        if (present(y)) then
            if (size(y,1).ne.size(z,1) .or. size(y,2).ne.size(z,2) .or. &
                &size(y,3).ne.size(z,3)) then
                call writo('WARNING: In print_GP, the size of y and z has to &
                    &be the same... Skipping plot')
                return
            end if
        end if
        
        ! set x
        if (present (x)) then
            x_fin = x
        else
            allocate(x_fin(size(z,1),size(z,2),size(z,3)))
            do ipntx = 1,npntx
                x_fin(ipntx,:,:) = ipntx
            end do
        end if
        ! set y
        if (present (y)) then
            y_fin = y
        else
            allocate(y_fin(size(z,1),size(z,2),size(z,3)))
            do ipnty = 1,npnty
                y_fin(:,ipnty,:) = ipnty
            end do
        end if
        
        ! set default file name if empty
        if (trim(file_name_i).eq.'') then
            file_name = 'temp_data_print_GP_3D_'//&
                &trim(i2str(grp_rank))//'.dat'
        else
            file_name = trim(file_name_i)
        end if
        
        ! open output file
        file_i = n_seq_0
        call safe_open(file_i,istat,data_dir//'/'//trim(file_name),'replace',&
            &'formatted',delim_in='none')
        
        ! write to output file
        write(file_i,*) '# '//trim(fun_name)//':'
        do ipntx = 1,npntx
            do ipnty = 1,npnty
                write(file_i,*) (x_fin(ipntx,ipnty,iplt), iplt = 1,nplt), &
                    &(y_fin(ipntx,ipnty,iplt), iplt = 1,nplt), &
                    &(z(ipntx,ipnty,iplt), iplt = 1,nplt)
            end do
            write(file_i,*) ''
        enddo 
        write(file_i,*) ''
        
        ! close output file
        close(file_i)
        
        if (trim(file_name_i).eq.'') then
            call execute_command_line('rm '//data_dir//'/'//trim(file_name))
        end if
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(fun_name,trim(file_name),nplt,.false.,.true.)
    end subroutine print_GP_3D_arr
    
    ! use GNUPlot to draw a plot
    ! The variables  fun_name and file_name  hold the name  of the plot  and the
    ! name of the file in which the plot is to be found. is_2D indicates whether
    ! the plot is 2D  or 3D and plot_on_screen indicates whether  the plot is to
    ! be shown on the screen, or to be saved in a file called [fun_name].pdf
    subroutine draw_GP_ind(fun_name,file_name,nplt,is_2D,plot_on_screen)
        use safe_open_mod, only: safe_open
        use num_vars, only: grp_rank, n_seq_0
        
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        
        ! local variables
        integer :: iplt                                                         ! counter
        integer :: istat                                                        ! status of opening a file
        integer :: cmd_i                                                        ! file number at which gnuplot_cmd opened
        character(len=max_str_ln) :: script_name                                ! name of script, including directory
        
        ! return if no_plots
        if (no_plots) then
            call writo('WARNING: plot ignored because no_plots is on')
            return
        end if
        
        ! create the GNUPlot command
        if (plot_on_screen) then
            script_name = trim(script_dir)//'/'//'temp_script_draw_GP_'//&
                &trim(i2str(grp_rank))//'.cmd'
        else
            script_name = trim(script_dir)//'/'//trim(fun_name)//'.cmd'
        end if
        
        cmd_i = n_seq_0
        call safe_open(cmd_i,istat,trim(script_name),'replace','formatted',&
            &delim_in='none')
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command')
            return
        end if
        
        ! create the GNUPlot command
        if (plot_on_screen) then
            ! basic GNUPlot commands
            !gnuplot_cmd = 'gnuplot -persist -e "set grid; set border 4095 &
                !&front linetype -1 linewidth 1.000; set terminal x11;'
            write(cmd_i,*) 'set grid; set border 4095 front linetype -1 &
                &linewidth 1.000; set terminal wxt;'
            ! no legend if too many plots
            if (nplt.gt.10) write(cmd_i,*) 'set nokey;'
            ! individual plots
            if (is_2D) then
                write(cmd_i,*) 'plot \'
                do iplt = 1,nplt
                    write(cmd_i,*) ' '''//trim(data_dir)//'/'//trim(file_name)&
                        &//''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                        &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                        &')'' with lines, \'
                end do
            else
                write(cmd_i,*) 'splot \'
                do iplt = 1,nplt
                    write(cmd_i,*) ' '''//trim(data_dir)//'/'//trim(file_name)&
                        &//''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//':'//trim(i2str(2*nplt+iplt))&
                        &//' title '''//trim(fun_name)//' ('//trim(i2str(iplt))&
                        &//'/'//trim(i2str(nplt))//')'' with lines, \'
                end do
            end if
            ! finishing the GNUPlot command
            write(cmd_i,*) ''
            write(cmd_i,*) 'pause -1'
        else
            ! basic GNUPlot commands
            write(cmd_i,*) 'set grid; set border 4095 front linetype -1 &
                &linewidth 1.000; set terminal pdf; set output '''//&
                &trim(plot_dir)//'/'//trim(fun_name)//'.pdf'';'
            ! no legend if too many plots
            if (nplt.gt.10) write(cmd_i,*) 'set nokey;'
            ! individual plots
            if (is_2D) then
                write(cmd_i,*) ' plot \'
                do iplt = 1,nplt
                    write(cmd_i,*) ' '''//trim(data_dir)//'/'//trim(file_name)&
                        &//''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                        &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                        &')'' with lines, \'
                end do
            else
                write(cmd_i,*) ' splot \'
                do iplt = 1,nplt
                    write(cmd_i,*) ' '''//trim(data_dir)//'/'//trim(file_name)&
                        &//''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//':'//trim(i2str(2*nplt+iplt))&
                        &//' title '''//trim(fun_name)//' ('//trim(i2str(iplt))&
                        &//'/'//trim(i2str(nplt))//')'' with lines, \'
                end do
            end if
            ! finishing the GNUPlot command
            write(cmd_i,*) ''
        end if
        
        ! closing the GNUPlot command
        close(cmd_i)
        
        ! stop timer
        call stop_time
        
        ! call GNUPlot
        call execute_command_line('gnuplot "'//trim(script_name)//&
            &'" 2> /dev/null',EXITSTAT=istat)
        
        if (plot_on_screen) then
            if (istat.ne.0) then
                call writo('Failed to plot '//trim(fun_name))
            end if
            call execute_command_line('rm '//trim(script_name))
        else
            if (istat.eq.0) then
                call writo('Created plot in output file '''//&
                    &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
            else
                call writo('Failed to create plot in output file '''//&
                    &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
            end if
        end if
        
        ! start timer
        call start_time
    end subroutine draw_GP_ind
    subroutine draw_GP_arr(fun_name,file_names,nplt,is_2D,plot_on_screen)
        use num_vars, only: grp_rank
        
        ! input / output
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        
        ! local variables
        character(len=max_str_ln) :: file_name                                  ! holds file name
        character(len=size(file_names)*max_str_ln) :: shell_cmd                 ! shell cmd
        integer :: istat                                                        ! status of opening a file
        integer :: id                                                           ! counter
        
        ! concatenate the files indicated by file_names in a temporary file
        ! set file name
        file_name = 'temp_file_draw_GP_'//trim(i2str(grp_rank))
        shell_cmd = 'cat'
        do id = 1,size(file_names)
            shell_cmd = trim(shell_cmd)//' "'//trim(data_dir)//'/'//&
                &trim(file_names(id))//'"'
        end do
        shell_cmd = trim(shell_cmd)//' > '//trim(data_dir)//'/'//trim(file_name)
        call execute_command_line(trim(shell_cmd),EXITSTAT=istat)
        
        ! call draw_GP_ind
        call draw_GP_ind(fun_name,file_name,nplt,is_2D,plot_on_screen)
        
        ! delete temporary file
        call execute_command_line('rm -f "'//trim(data_dir)//'/'&
            &//trim(file_name)//'"',EXITSTAT=istat)
        
        if (istat.ne.0) then
            call writo('Failed to create plot in output file '''//&
                &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
        end if
    end subroutine draw_GP_arr
    
    ! use GNUPlot to draw an animated plot in an animated gif
    ! The variables  fun_name and file_name  hold the name  of the plot  and the
    ! name of the file in which the plot is to be found. is_2D indicates whether
    ! the plot is 2D or 3D. The plot is saved in a file called [fun_name].gif
    subroutine draw_GP_animated_ind(fun_name,file_name,nplt,is_2D,ranges)
        use num_vars, only: n_seq_0
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        
        ! local variables
        integer :: iplt                                                         ! counter
        integer :: istat                                                        ! status of opening a file
        integer :: cmd_i                                                        ! file number at which gnuplot_cmd opened
        character(len=max_str_ln) :: script_name                                ! name of script, including directory
        real(dp), allocatable :: ranges_loc(:,:)                                ! local copy of ranges
        
        ! return if no_plots
        if (no_plots) then
            call writo('WARNING: plot ignored because no_plots is on')
            return
        end if
        
        ! create the GNUPlot command
        script_name = ''//trim(script_dir)//'/'//trim(fun_name)//'.cmd'
        cmd_i = n_seq_0
        call safe_open(cmd_i,istat,trim(script_name),'replace','formatted',&
            &delim_in='none')
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command')
            return
        end if
        
        ! basic GNUPlot commands
        write(cmd_i,*) 'set grid; set border 4095 front &
            &linetype -1 linewidth 1.000; set terminal gif animate delay '//&
            &trim(i2str(max(250/nplt,1)))//' size 1280, 720; &
            &set output '''//trim(plot_dir)//'/'//trim(fun_name)//&
            &'.gif'';'
        
        ! individual plots
        if (is_2D) then
            ! find ranges
            allocate(ranges_loc(2,2))
            if (present(ranges)) then
                if (size(ranges,1).eq.2 .and. size(ranges,2).eq.2) then
                    ranges_loc = ranges
                else
                    call writo('WARNING: invalid ranges given to &
                        &draw_GP_animated')
                end if
            else
                call get_ranges(ranges_loc)
            end if
            
            ! set ranges
            write(cmd_i,*) 'set xrange ['//trim(r2str(ranges_loc(1,1)))//':'//&
                &trim(r2str(ranges_loc(1,2)))//'];'
            write(cmd_i,*) 'set yrange ['//trim(r2str(ranges_loc(2,1)))//':'//&
                &trim(r2str(ranges_loc(2,2)))//'];'
            
            ! plot
            do iplt = 1,nplt
                write(cmd_i,*) ' plot'
                write(cmd_i,*) ' '''//trim(data_dir)//'/'//trim(file_name)//&
                    &''' using '//trim(i2str(iplt))//':'//&
                    &trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                    &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                    &')'' with lines'
            end do
        else
            ! find ranges
            allocate(ranges_loc(3,2))
            if (present(ranges)) then
                if (size(ranges,1).eq.3 .and. size(ranges,2).eq.2) then
                    ranges_loc = ranges
                else
                    call writo('WARNING: invalid ranges given to &
                        &draw_GP_animated')
                end if
            else
                call get_ranges(ranges_loc)
            end if
            
            ! set ranges
            write(cmd_i,*) ' set xrange ['//trim(r2str(ranges_loc(1,1)))//':'//&
                &trim(r2str(ranges_loc(1,2)))//'];'
            ! y range
            write(cmd_i,*) ' set yrange ['//trim(r2str(ranges_loc(2,1)))//':'//&
                &trim(r2str(ranges_loc(2,2)))//'];'
            ! z range
            write(cmd_i,*) ' set zrange ['//trim(r2str(ranges_loc(3,1)))//':'//&
                &trim(r2str(ranges_loc(3,2)))//'];'
            ! color scale
            write(cmd_i,*) ' set cbrange ['//trim(r2str(ranges_loc(3,1)))//':'//&
                &trim(r2str(ranges_loc(3,2)))//'];'
            
            ! plot
            write(cmd_i,*) &
                &'set style line 1 linecolor rgb ''#cccccc''; set pm3d at b;'
            write(cmd_i,*) 'set view 45,45;'
            do iplt = 1,nplt
                write(cmd_i,*) 'splot '''//trim(data_dir)//'/'//trim(file_name)//&
                    &''' using '//trim(i2str(iplt))//':'//&
                    &trim(i2str(nplt+iplt))//':'//&
                    &trim(i2str(2*nplt+iplt))//' title '''//trim(fun_name)//&
                    &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                    &')'' with lines linestyle 1'
            end do
        end if
        ! finishing and closing the GNUPlot command
        write(cmd_i,*) 'set output'
        close(cmd_i)
        
        ! no need to stop the time
        
        ! call GNUPlot
        call execute_command_line('gnuplot "'//trim(script_name)//&
            &'" 2> /dev/null',EXITSTAT=istat)
        
        if (istat.eq.0) then
            call writo('Created animated plot in output file '''//&
                &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
        else
            call writo('Failed to create animated plot in output file '''//&
                &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
        end if
    contains
        !subroutine set_range(rng,var,cmd_i)
            !! input / output
            !integer, intent(in) :: rng(2)
            !character(len=*), intent(in) :: var
            !integer, intent(in) :: cmd_i
            
            !! update gnuplot_cmd
            !write(cmd_i,*) 'do for [t='//trim(i2str(rng(1)))//':'//&
                !&trim(i2str(rng(2)))//'] {'
            !write(cmd_i,*) ' stats "'//trim(data_dir)//'/'//trim(file_name)//&
                !&'" using t nooutput;'
            !write(cmd_i,*) ' if (t=='//trim(i2str(rng(1)))//&
                !&') {min_'//trim(var)//'=STATS_min; max_'//trim(var)//&
                !&'=STATS_max;} else { if (STATS_min < min_'//trim(var)//') &
                !&{ min_'//trim(var)//' = STATS_min} &
                !&if (STATS_max > max_'//trim(var)//') &
                !&{ max_'//trim(var)//' = STATS_max} } }'
            !write(cmd_i,*) 'set '//trim(var)//'range [min_'//trim(var)//&
                !&':max_'//trim(var)//'];'
        !end subroutine set_range
        
        subroutine get_ranges(ranges)
            use safe_open_mod, only: safe_open
            use num_vars, only: n_seq_0
            
            ! input / output
            real(dp), intent(inout) :: ranges(:,:)                              ! x and y range, and z range (if 3D) of plot
            
            ! local variables
            integer :: data_i                                                   ! file number of data file
            real(dp), allocatable :: loc_data(:)                                ! one line of data
            character(len=1) :: loc_data_char
            
            
            ! open data file
            data_i = n_seq_0
            call safe_open(data_i,istat,trim(data_dir)//'/'//trim(file_name),&
                &'old','formatted',delim_in='none')
            
            ! set up loc_data
            if (is_2D) then
                allocate(loc_data(2*nplt))
            else
                allocate(loc_data(3*nplt))
            end if
            
            ! initialize ranges
            ranges(:,1) = 1.E14_dp                                              ! minimum value
            ranges(:,2) = -1.E14_dp                                             ! maximum value
            
            ! read the data file
            istat = 0
            do while (istat.eq.0)
                read(data_i,*,IOSTAT=istat) loc_data_char                       ! read first character of data
                if (istat.eq.0) then                                            ! read succesful
                    if (loc_data_char.ne.'#') then                              ! exclude comment lines
                        backspace(UNIT=data_I)                                  ! go back one line
                        read(data_i,*,IOSTAT=istat) loc_data                    ! read data again, but now ful
                        ranges(1,1) = &
                            &min(ranges(1,1),minval(loc_data(1:nplt)))
                        ranges(1,2) = &
                            &max(ranges(1,2),maxval(loc_data(1:nplt)))
                        ranges(2,1) = &
                            &min(ranges(2,1),minval(loc_data(nplt+1:2*nplt)))
                        ranges(2,2) = &
                            &max(ranges(2,2),maxval(loc_data(nplt+1:2*nplt)))
                        if (.not.is_2D) then
                            ranges(3,1) = min(ranges(3,1),&
                                &minval(loc_data(2*nplt+1:3*nplt)))
                            ranges(3,2) = max(ranges(3,2),&
                                &maxval(loc_data(2*nplt+1:3*nplt)))
                        end if
                    end if
                end if
            end do
            !write(*,*) 'ranges(1,:) = ', ranges(1,:)
            !write(*,*) 'ranges(2,:) = ', ranges(2,:)
            !write(*,*) 'ranges(3,:) = ', ranges(3,:)
            
            ! close data file
            close(data_i)
        end subroutine get_ranges
    end subroutine draw_GP_animated_ind
    subroutine draw_GP_animated_arr(fun_name,file_names,nplt,is_2D,ranges)
        use num_vars, only: grp_rank
        
        ! input / output
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        
        ! local variables
        character(len=max_str_ln) :: file_name                                  ! holds file name
        character(len=size(file_names)*max_str_ln) :: shell_cmd                 ! shell cmd
        integer :: istat                                                        ! status of opening a file
        integer :: id                                                           ! counter
        
        ! concatenate the files indicated by file_names in a temporary file
        ! set file name
        file_name = 'temp_file_draw_GP_animated_'//trim(i2str(grp_rank))
        shell_cmd = 'cat'
        do id = 1,size(file_names)
            shell_cmd = trim(shell_cmd)//' "'//trim(data_dir)//'/'//&
                &trim(file_names(id))//'"'
        end do
        shell_cmd = trim(shell_cmd)//' > '//trim(data_dir)//'/'//trim(file_name)
        call execute_command_line(trim(shell_cmd),EXITSTAT=istat)
        
        ! call draw_GP_animated_ind
        call draw_GP_animated_ind(fun_name,file_name,nplt,is_2D,ranges)
        
        ! delete temporary file
        call execute_command_line('rm -f "'//trim(data_dir)//'/'//&
            &trim(file_name)//'"',EXITSTAT=istat)
        
        if (istat.ne.0) then
            call writo('Failed to create plot in output file '''//&
                &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
        end if
    end subroutine draw_GP_animated_arr
    
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
        use num_vars, only: grp_rank, glb_rank, output_i
        
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
        
        ! setup ignore
        ignore = .true.                                                         ! ignore by default
        if (grp_rank.eq.0) ignore = .false.                                     ! group masters can output
        if (present(persistent)) ignore = .not.persistent                       ! persistent can override this
        
        if (.not.ignore) then                                                   ! only group master (= global master if no groups) or persistent
            ! Divide the input string length by the max_str_ln and loop over the
            ! different parts
            max_len_part = max_str_ln-(lvl-1)*len(lvl_sep) - 10                 ! max length of a part (10 for time)
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
                do id = 1, len(trim(output_str)) - 10                           ! 10 for time
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
                                call format_str(lvl,temp_output_err_str(id))   ! format error message
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
                
                ! if first group, also write output to screen
                if (glb_rank.eq.0) then                                         ! if output_i = 0, output has already been written to screen
                    if (lvl.eq.1) write(*,*) header_str                         ! first level gets extra lines
                    write(*,*) output_str
                    if (lvl.eq.1) write(*,*) header_str
                end if
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
end module output_ops

    !! PREVIOUS OUTPUT ROUTINE, REPLACED BY PRINT_GP
    !! write output using method indicated by format_out
    !!     format_out = 1: NETCDF output
    !!     format_out = 2: matlab output
    !!     format_out = 3: DISLIN output
    !!     format_out = 4: GNUplot output
    !! They can all be called for 2D or 3D plots
    !! For  2D plots  either  x-axis values  can  be specified,  or  they can  be
    !! automatically be inserted as 1 to nx
    !! For 3D plots there is currently no option to specify the axis
    !! A special case applies for GNUplot output. By applying comment='multiplot',
    !! a 3D array is interpreted as a multiplot of 2D plots
    !subroutine write_out(nx, ny, fun, fun_name, output_i_in, comment)
        !use num_vars, only: output_i, n_seq_0
        !use safe_open_mod, only: safe_open
        
        !! input / output
        !integer, intent(in) :: nx, ny
        !character(len=*) :: fun_name
        !real(dp) :: fun(1:nx,1:ny)                                              ! explicit shape is used so first index can be left out if equal to 1
        !integer, intent(in), optional :: output_i_in
        !character(len=*), optional :: comment
        
        !! local variables
        !integer :: fin_output_i
        !real(dp), allocatable :: fun_alt(:,:)
        !integer :: id
        !integer :: istat
        !logical :: write_x11
        
        !! initialize
        !write_x11 = .false.
        
        !if (present(output_i_in)) then 
            !fin_output_i = output_i_in
        !else
            !fin_output_i = output_i
        !end if
        
        !! choose correct format depending on format_out
        !select case (format_out)
            !case (1)                                                            ! NETCDF
                !call writo('WARNING: NETCDF output not yet implemented')
            !case (2)                                                            ! matlab
                !call write_out_matlab(fin_output_i, nx, ny, fun, fun_name, &
                    !&comment)
            !case (3)
                !if (nx.eq.1) then                                               ! DISLIN 2D without x-vector
                    !! allocate a new variable fun_alt that holds a standard axis
                    !! and the original data in the y-axis
                    !allocate(fun_alt(2,size(fun,2))); fun_alt = 0.0_dp
                    !fun_alt(1,:) = (/ ( id , id = 1,size(fun,2)) /)
                    !fun_alt(2,:) = fun(1,:)
                    !call write_out_dislin_2D(ny, fun_alt, fun_name, comment)
                !else if (nx.eq.2) then                                          ! DISLIN 2D
                    !call write_out_dislin_2D(ny, fun, fun_name, comment)
                !else if (nx.gt.2) then                                          ! DISLIN 3D
                    !call write_out_dislin_3D(nx, ny, fun, fun_name, comment)
                !else
                    !call writo('WARNING: Dimension 1 has to be at least 1')
                !end if
            !case (4)                                                            ! GNUplot
                !! if output_i is equal to zero, the plots are not done to a file
                !! but directly to the X11 screen :open a new temporary file
                !if (fin_output_i.eq.0) then
                    !fin_output_i = n_seq_0
                    !call safe_open(fin_output_i,istat,'tempoutput.dat',&
                        !&'replace','formatted',delim_in='none')
                    !write_x11 = .true.
                !end if
                !if (nx.eq.1) then
                    !! allocate a new variable fun_alt that holds a standard axis
                    !! and the original data in the y-axis
                    !allocate(fun_alt(2,size(fun,2))); fun_alt = 0.0_dp
                    !fun_alt(1,:) = (/ ( id , id = 1,size(fun,2)) /)
                    !fun_alt(2,:) = fun(1,:)
                    !call write_out_gnuplot_2D(fin_output_i, ny, fun_alt, &
                        !&fun_name, comment)
                !else if (nx.eq.2) then
                    !call write_out_gnuplot_2D(fin_output_i, ny, fun, &
                        !&fun_name, comment)
                !else
                    !if (trim(comment).eq.'multiplot') then
                        !call write_out_gnuplot_2D_multi(fin_output_i, ny, &
                            !&fun, fun_name)
                    !else
                        !call write_out_gnuplot_3D(fin_output_i, nx, ny, fun, &
                            !&fun_name, comment)
                    !end if
                !end if
                !!gnuplot -e "splot 'results.txt' with lines; set contour; set cntrparam levels auto 10; pause -1"
                !close(fin_output_i)
                !if (write_x11) call execute_command_line('rm tempoutput.dat')
            !case (5)                                                            ! matplotlib
                !call writo('WARNING: matplotlib output not yet implemented')
                !!!! MAKE A PLOT LIKE IN CONTOURF3D_DEMO2.PY !!!
            !case default 
                !call writo('WARNING: No output format associated with ' &
                    !&// i2str(format_out) )
        !end select
    !contains
        !! writes a 2D array into matlab format
        !subroutine write_out_matlab(output_i, nx, ny, fun, fun_name, comment)
            !! input / output
            !integer, intent(in) :: output_i, nx, ny
            !character(len=*) :: fun_name
            !real(dp) :: fun(1:nx,1:ny)
            !character(len=*), optional :: comment
            
            !! local variables
            !integer :: ix, iy
            !character(len=max_str_ln) :: form_str
            
            !form_str = ''
            !form_str = '(A,' // trim(i2str(ny)) // '(17x,i3,16x))'
            
            !if (present(comment)) then
                !write(output_i,*) '% comment: ' // trim(comment)
            !end if
            !write(output_i,*) trim(fun_name) // ' = ['
            !write(output_i,form_str) '%', (iy,iy=1,ny)
            !form_str = ''
            !form_str = '(' // trim(i2str(ny)) // &
                !&'(es36.22e4),1x,A,1x,i3,1x)'
            !do ix=1,nx
                !write(output_i,form_str) (fun(ix,iy), iy = 1, ny), '%', ix 
            !enddo 
            !write(output_i,*) '];'
        !end subroutine
        
        !! writes a 2D array onto the screen using DISLIN
        !subroutine write_out_dislin_2D(np, fun, fun_name, comment)
            !!use dislin
            
            !! input / output
            !integer, intent(in) :: np
            !character(len=*) :: fun_name
            !real(dp) :: fun(2,np)
            !character(len=*), optional :: comment
            
            !!! local variables
            !!integer :: ic
            !!real(dp) :: t_axis(4), f_axis(4)                                    ! axis that envelops minimum, maximum of data
            !!real(dp) :: x_val(np), y_val(np)                                    ! to avoid annoying compiler performance warnings
            
            !!call metafl('xwin')                                                 ! xWin terminal
            !!!call window(0,0,1024,512)                                           ! set window size
            !!call setpag('DA4L')                                                 ! A4 landscape
            !!!call sclfac(0.4_dp)                                                 ! scale factor
            !!call disini                                                         ! start DISFIN
            !!call winkey('return')                                               ! return key also exits
            !!call errmod('all','off')                                            ! disable all messages
            !!call errmod('warnings','on')                                        ! turn warnings on
            
            !!call pagfll(255)                                                    ! white background
            !!call setrgb(0.0_dp,0.0_dp,0.0_dp)                                   ! black foreground for rest
            !!call pagera()                                                       ! plots page border
            !!call complx()                                                       ! sets complex font
            !!!call axspos(451,1800)                                               ! position of lower left corner of axis
            !!!call axslen(2200,1200)                                              ! length of axis
            
            !!call name('x','x')                                                  ! name and label of x axis
            !!call name('f','y')                                                  ! name and label of y axis
            !!call ticks(5,'x')
            !!call ticks(5,'y')
            
            !!call titlin ('f(x) = '//fun_name , 1)                               ! main title
            !!call titlin (comment, 3)                                            ! subtitle
            
            !!ic=intrgb(0.95_dp,0.95_dp,0.95_dp)                                  ! light grey in RGB
            !!call axsbgd(ic)                                                     ! background color for axis
            
            !!! min and max of axis
            !!t_axis = det_axis(fun(1,:),10)
            !!f_axis = det_axis(fun(2,:),10)
            !!if (abs(f_axis(1)-f_axis(2)).lt.1E-5_dp) f_axis = &
                !!&[f_axis(1)-1,f_axis(1)+1,f_axis(1)-1,0.5_dp]
            !!call labdig(max(0,-floor(log10(t_axis(4)))),'x')                    ! number of digits after decimal point in lables
            !!call labdig(max(0,-floor(log10(f_axis(4)))),'y')                    ! number of digits after decimal point in lables
            
            !!call graf(t_axis(1),t_axis(2),t_axis(3),t_axis(4),&
                !!&f_axis(1),f_axis(2),f_axis(3),f_axis(4))
            !!call title
            
            !!call color('red')                                                   ! red color
            !!call thkcrv(5)
            !!x_val = fun(1,:); y_val = fun(2,:)
            !!call curve(x_val,y_val,np)                                          ! plot 2D curve
            !!call thkcrv(1)
            
            !!call color('black')                                                 ! black color
            !!!call axgit                                                          ! plot axis
            !!call dash                                                           ! dashed line style
            !!call xaxgit                                                         ! plot only x axis
            !!call solid                                                          ! solid line style
            !!!call yaxgit                                                         ! plot only y axis
            
            !!call setrgb(0.5_dp,0.5_dp,0.5_dp)                                   ! dark grey foreground for grid
            !!call grid(1,1)                                                      ! 1 line per 'ticks' for x-axis and also for y
            !!call setrgb(0.0_dp,0.0_dp,0.0_dp)                                   ! black foreground for rest
            
            !!call disfin                                                         ! terminate DISFIN
        !end subroutine
        
        !! writes a surfrace plot onto the screen using DISLIN
        !subroutine write_out_dislin_3D(nx, ny, fun, fun_name, comment)
            !!use dislin
            
            !! input / output
            !integer, intent(in) :: nx, ny
            !character(len=*) :: fun_name
            !real(dp) :: fun(1:nx,1:ny)
            !character(len=*), optional :: comment
            
            !!! local variables
            !!integer :: ic, id
            !!real(dp) :: f_axis(4)                                               ! axis that envelops minimum, maximum of data
            !!real(dp) :: x_arr(nx), y_arr(ny)
            !!real(dp) :: minx, maxx, miny, maxy
            !!real(dp) :: fun_joined(1:nx*ny)
            !!!real(dp) :: zlev                                                    ! for surface plots
            
            !!call metafl('xwin')                                                 ! xWin terminal
            !!!call window(0,0,1024,512)                                           ! set window size
            !!call setpag('DA4L')                                                 ! A4 landscape
            !!!call sclfac(0.4_dp)                                                 ! scale factor
            !!call disini                                                         ! start DISFIN
            !!call winkey('return')                                               ! return key also exits
            !!call errmod('all','off')                                            ! disable all messages
            !!call errmod('warnings','on')                                        ! turn warnings on
            
            !!call pagfll(255)                                                    ! white background
            !!call setrgb(0.0_dp,0.0_dp,0.0_dp)                                   ! black foreground for rest
            !!call pagera()                                                       ! plots page border
            !!call complx()                                                       ! sets complex font
            !!!call axspos(451,1800)                                               ! position of lower left corner of axis
            !!!call axslen(2200,1200)                                              ! length of axis

            !!call name('x','x')                                                  ! name and label of x axis
            !!call name('y','y')                                                  ! name and label of y axis
            !!call name('f','z')                                                  ! name and label of z axis
            !!call ticks(5,'xyz')

            !!call titlin ('f = '//fun_name , 1)                                  ! main title
            !!call titlin (comment, 3)                                            ! subtitle
            
            !!ic=intrgb(0.95_dp,0.95_dp,0.95_dp)                                  ! light grey in RGB
            !!call axsbgd(ic)                                                     ! background color for axis
            
            !!! min and max of axis
            !!minx = 1; maxx = size(fun,1)
            !!miny = 1; maxy = size(fun,2)
            !!do id = 1,ny
                !!fun_joined((id-1)*nx+1:id*nx) = fun(:,id)
            !!end do
            !!f_axis = det_axis(fun_joined,10)
            !!call labdig(-1,'xy')
            !!call labdig(max(0,-floor(log10(f_axis(4)))),'z')                    ! number of digits after decimal point in lables
            
            !!call graf3D(minx,maxx,minx,(maxx-minx)/6,miny,maxy,miny,&           ! (shaded) surface plot
                !!&(maxy-miny)/6,f_axis(1),f_axis(2),f_axis(3),f_axis(4))
            !!!call graf(minx,maxx,minx,(maxx-minx)/6,miny,maxy,miny,&             ! contour plot
                 !!!&(maxy-miny)/6)
            !!call box3d                                                          ! 3D box
            !!call title
            
            !!x_arr = (/( id, id = 1,nx)/)
            !!y_arr = (/( id, id = 1,ny)/)
            !!call color('red')                                                   ! red color
            !!call thkcrv(5)
            !!call surmat(fun,nx,ny,1,1)                                          ! surface plot
            !!!call shdmod('smooth','sufrace')                                    ! shaded surface plot
            !!!call surshd(x_arr,nx,y_arr,ny,fun)
            !!!do id = 1,9                                                         ! contour plot
                !!!zlev = minf + (id-1)*(maxf-minf)/10
                !!!call setclr(id*25)
                !!!if(id.eq.5) then
                  !!!call labels('none','contur')
                !!!else
                  !!!call labels('float','contur')
                !!!end if
                !!!call contur(x_arr,nx,y_arr,ny,fun,zlev)
            !!!end do
            !!call thkcrv(1)
            !!call color('black')                                                 ! black color
            
            !!call disfin                                                         ! terminate DISFIN
        !end subroutine
        
        !! determine the best axis for a data range of the form:
        !! (from http://www2.mps.mpg.de/dislin/kap4.html)
        !!   det_axis(1) = minimum
        !!   det_axis(2) = maximum
        !!   det_axis(3) = where to start the first label
        !!   det_axis(4) = step betwen labels
        !! (minf,maxf,minf,(maxf-minf)/6)
        !function det_axis(var, n_points)
            !! input / output
            !integer :: n_points
            !real(dp) :: var(:)
            !real(dp) :: det_axis(4)
            
            !! local variables
            !real(dp) :: delta_var                                               ! delta from data
            !integer :: base_exp_x2                                              ! exponent of 2X delta from data
            !real(dp) :: aux                                                     ! auxiliary variable
            !integer :: aux2                                                     ! auxiliary variable
            !real(dp) :: delta_ax                                                ! delta used in axis
            !real(dp) :: l_b, u_b                                                ! lower and upper bound
            !real(dp) :: margin                                                  ! margin for comparison
            
            !! initialize output
            !det_axis = 0.0_dp
            
            !! find the distance between labels on the axis by dividing the range
            !! of  the input  data  (var) by  the number  of  points wanted,  and
            !! rounding to the nearest multiple of 1Ej * 0.5 or 1Ej * 1.0
            !! First find the value of j. E.g.: if delta_var = 0.036 and n_points
            !! = 6,  the value delta_var/n_points  = 0.006  has to be  rounded to
            !! 0.005, which means that j = 2
            !! Then, the nearest minimum and  maximum that are multiples of 0.005
            !! are found, that encompass the data in var
            !l_b = minval(var)
            !u_b = maxval(var)
            !delta_var = u_b - l_b
            !aux = 2.0_dp*delta_var/n_points
            !base_exp_x2 = floor(log10(aux)) + 1                                 ! so that we get a number of the form 0.abcd
            !aux2 = nint(aux*10.0_dp**(-base_exp_x2))
            !if (aux2.eq.0) then
                !delta_ax = 10.0_dp**(base_exp_x2)*0.1_dp                        ! so that we get 0.1 from aux2/2
            !else
                !delta_ax = 10.0_dp**(base_exp_x2)*aux2*0.5_dp                   ! so that we get 0.5 from aux2/2
            !end if
            
            !! determine the lower and upper  bounds by checking whether there is
            !! a remainder when  dividing the lowest, highest  value by delta_ax.
            !margin = 1.0E-5_dp                                                  ! visually indistinguishable
            !if (abs(dmod(l_b,delta_ax)).gt.margin) then                         ! take rounded times delta
                !det_axis(1) = ( floor(l_b/delta_ax) ) * delta_ax
            !else
                !det_axis(1) = l_b
            !end if
            !if (abs(dmod(u_b,delta_ax)).gt.margin) then                         ! take rounded times delta
                !det_axis(2) = ( ceiling(u_b/delta_ax) ) * delta_ax
            !else
                !det_axis(2) = u_b
            !end if
            
            !! starting point of axis
            !det_axis(3) = det_axis(1)
            
            !! determine the step size that is a multiple of delta_ax, closest to
            !! delta_var/n_points = aux/2 -> step size * delta_ax
            !det_axis(4) = nint(0.5_dp*aux/delta_ax)*delta_ax
        !end function det_axis
        
        !! writes a 2D array into GNUplot format
        !subroutine write_out_gnuplot_2D(output_i, ny, fun, fun_name, comment)
            !! input / output
            !integer, intent(in) :: output_i, ny
            !character(len=*) :: fun_name
            !real(dp) :: fun(1:,1:)
            !character(len=*), optional :: comment
            
            !! local variables
            !integer :: iy
            !character(len=2*max_str_ln) :: gnuplot_cmd
            
            !! test
            !if (size(fun,1).ne.2 .or. size(fun,2).ne.ny) then
                !call writo('WARNING: Write_out_gnuplot_2D has to be called &
                    !&with an array of dimension [2,ny]')
                !return
            !end if
            
            !if (present(comment)) then
                !write(output_i,*) '# comment: ' // trim(comment)
            !end if
            !write(output_i,*) '# '//trim(fun_name)//':'
            !do iy=1,ny
                !write(output_i,*) fun(1,iy), fun(2,iy)
            !enddo 
            !write(output_i,*) ''
            
            !if (write_x11) then
                !gnuplot_cmd = 'gnuplot -e "plot ''tempoutput.dat'' using 1:2 &
                    !&title '''//trim(fun_name)//' ('//trim(comment)//' )'' &
                    !&with lines; set grid; replot; pause -1"'
                !call execute_command_line(gnuplot_cmd)
            !end if
        !end subroutine
        
        !! multiplot version of 2D plot
        !subroutine write_out_gnuplot_2D_multi(output_i, ny, fun, fun_name)
            !! input / output
            !integer, intent(in) :: output_i, ny
            !character(len=*) :: fun_name
            !real(dp) :: fun(1:,1:)
            
            !! local variables
            !integer :: iy, ix
            !character(len=50*max_str_ln) :: gnuplot_cmd
            
            !if (present(comment)) then
                !write(output_i,*) '# comment: ' // trim(comment)
            !end if
            !write(output_i,*) '# '//trim(fun_name)//':'
            !do iy=1,ny
                !write(output_i,*) iy, (fun(ix,iy), ix=1,size(fun,1))
            !enddo 
            !write(output_i,*) ''
            
            !if (write_x11) then
                !gnuplot_cmd = 'gnuplot -e "plot'
                !do ix = 2,nx+1
                    !gnuplot_cmd = trim(gnuplot_cmd) // ' ''tempoutput.dat'' &
                        !&using 1:'//trim(i2str(ix))//' title '''//&
                        !&trim(fun_name)//'_'//trim(i2str(ix))//''' with lines'
                    !if (ix.lt.nx+1) gnuplot_cmd = trim(gnuplot_cmd)//','
                !end do
                !gnuplot_cmd = trim(gnuplot_cmd)//'; set grid; replot; pause -1"'
                !call execute_command_line(gnuplot_cmd)
            !end if
        !end subroutine
        
        !! writes a 3D array into GNUplot format
        !subroutine write_out_gnuplot_3D(output_i, nx, ny, fun, fun_name, comment)
            !! input / output
            !integer, intent(in) :: output_i, nx, ny
            !character(len=*) :: fun_name
            !real(dp) :: fun(1:nx,1:ny)
            !character(len=*), optional :: comment
            
            !! local variables
            !integer :: ix, iy
            !character(len=2*max_str_ln) :: gnuplot_cmd
            
            !if (present(comment)) then
                !write(output_i,*) '# comment: ' // trim(comment)
            !end if
            !write(output_i,*) '# "'//trim(fun_name)//'":'
            !do ix=1,nx
                !do iy = 1,ny
                    !write(output_i,*) ix, iy, fun(ix,iy)
                !end do
                !write(output_i,*) ''
            !end do 
            
            !if (write_x11) then
                !gnuplot_cmd = 'gnuplot -e "splot'
                !gnuplot_cmd = trim(gnuplot_cmd)//" 'tempoutput.dat'"
                !gnuplot_cmd = trim(gnuplot_cmd)//' title '''//trim(fun_name)&
                    !&//'_'//trim(i2str(ix))//'''with lines; set contour; &
                    !&set cntrparam levels auto 10; set pm3d; replot; pause -1"'
                !call execute_command_line(gnuplot_cmd)
            !end if
        !end subroutine
    !end subroutine
