!------------------------------------------------------------------------------!
!   Operations concerning  giving output, on the  screen as well as  in output !
!   files                                                                      !
!------------------------------------------------------------------------------!
module output_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln, no_plots, iu, plot_dir, data_dir, &
        &script_dir, plot_size
    use files_utilities, only: nextunit
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    
    
    implicit none
    private
    public print_GP_2D, print_GP_3D, draw_GP, draw_GP_animated, merge_GP, &
        &plot_HDF5, plot_diff_HDF5
    
    ! interfaces
    interface print_GP_2D
        module procedure print_GP_2D_ind, print_GP_2D_arr
    end interface
    interface print_GP_3D
        module procedure print_GP_3D_ind, print_GP_3D_arr
    end interface
    interface plot_HDF5
        module procedure plot_HDF5_ind, plot_HDF5_arr
    end interface
    interface draw_GP
        module procedure draw_GP_ind, draw_GP_arr
    end interface
    interface draw_GP_animated
        module procedure draw_GP_animated_ind, draw_GP_animated_arr
    end interface
    
    ! global variables
    character(len=9) :: line_clrs(9) = ["'#000090'","'#000fff'","'#0090ff'",&
        &"'#0fffee'", "'#90ff70'", "'#ffee00'", "'#ff7000'", "'#ee0000'", &
        &"'#7f0000'"]                                                           ! color specifications for plotting
    !character(len=max_str_ln) :: line_style = 'lt 1 lw 1 pt 7 pi -1 ps 0.5; &
        !&set pointintervalbox 0.75;'                                            ! with little space in line connecting points
    character(len=max_str_ln) :: line_style = 'lt 1 lw 1 pt 7 ps 0.5;'          ! without little space in line connecting points
#if ldebug
    character(len=0) :: err_output_str = ''                                     ! string with error output
#else
    character(len=14) :: err_output_str = ' 2> /dev/null'                       ! string with error output (/dev/null)
#endif
    
contains
    ! print 2D output on a file
    ! The variables var_name and file_name hold the  name of the plot and of the
    ! file in  which the plot  data is  to be saved,  respectively. y is  the an
    ! array containing the function which is  stored and x is an optional vector
    ! with the  x-values. The  logical draw=.false. optionally  disables calling
    ! the  GNUPlot drawing  procedure for  output on  screen [default],  without
    ! modifying the plot file
    ! The first index of y (and x) contains the points of a current plot
    ! The second index indicates various plots (one or more)
    subroutine print_GP_2D_ind(var_name,file_name,y,x,draw)                     ! individual plot
        ! input / output
        character(len=*), intent(in) :: var_name
        character(len=*), intent(in) :: file_name
        real(dp), intent(in) :: y(1:)
        real(dp), intent(in), optional :: x(1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: npoints
        
        ! set npoints, ny
        npoints = size(y)
        
        ! call multiplot version
        if (present(x)) then
            call print_GP_2D_arr(var_name,file_name,reshape(y,[npoints,1]),&
                &x=reshape(x,[npoints,1]),draw=draw)
        else
            call print_GP_2D_arr(var_name,file_name,reshape(y,[npoints,1]),&
                &draw=draw)
        end if
    end subroutine print_GP_2D_ind
    subroutine print_GP_2D_arr(var_name,file_name_i,y,x,draw)                   ! multiple plots
        use num_vars, only: rank
        
        ! input / output
        character(len=*), intent(in) :: var_name
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
        
        ! set nplt, npnt
        npnt = size(y,1)
        nplt = size(y,2)
        if (present(x)) then
            if (size(x,1).ne.size(y,1) .or. size(x,2).ne.size(y,2)) then
                call writo('WARNING: In print_GP_2D, the size of x and y has &
                    &to be the same... Skipping plot',persistent=.true.)
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
            file_name = 'temp_data_print_GP_2D_'//trim(i2str(rank))
        else
            file_name = trim(file_name_i)
        end if
        
        ! open output file
        open(unit=nextunit(file_i),file=data_dir//'/'//trim(file_name)//'.dat',&
            &iostat=istat)
        
        ! write to output file
        write(file_i,*) '# '//trim(var_name)//':'
        do ipnt = 1,npnt
            write(file_i,*) (x_fin(ipnt,iplt), iplt = 1,nplt), &
                &(y(ipnt,iplt), iplt = 1,nplt)
        enddo 
        write(file_i,*) ''
        
        ! close output file
        close(file_i)
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(var_name,file_name,file_name,nplt,1,.true.)
        
        if (trim(file_name_i).eq.'') then
            call use_execute_command_line('rm '//data_dir//'/'//&
                &trim(file_name)//'.dat')
        end if
    end subroutine print_GP_2D_arr
    
    ! print 3D output on a file or on the screen, using GNUPlot
    ! The variables var_name and file_name hold the  name of the plot and of the
    ! file in  which the plot  data is  to be saved,  respectively. z is  the an
    ! array containing  the function which  is stored and  x and y  are optional
    ! vectors  with the  x  and y-values.  The  logical draw=.false.  optionally
    ! disables  calling  the GNUPlot  drawing  procedure  for output  on  screen
    ! [default], without modifying the plot file
    ! The first index of z (and x, y) contains the points of a current plot
    ! The second index indicates various plots (one or more)
    subroutine print_GP_3D_ind(var_name,file_name,z,y,x,draw)                   ! individual plot
        ! input / output
        character(len=*), intent(in) :: var_name
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
        if (present(y)) then
            if (present(x)) then
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[size(z,1),size(z,2),1]),&
                    &y=reshape(y,[size(z,1),size(z,2),1]),&
                    &x=reshape(x,[size(z,1),size(z,2),1]),draw=draw)
            else
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[size(z,1),size(z,2),1]),&
                    &y=reshape(y,[size(z,1),size(z,2),1]),&
                    &draw=draw)
            end if
        else
            if (present(x)) then
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[size(z,1),size(z,2),1]),&
                    &x=reshape(x,[size(z,1),size(z,2),1]),&
                    &draw=draw)
            else
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[size(z,1),size(z,2),1]),draw=draw)
            end if
        end if
    end subroutine print_GP_3D_ind
    subroutine print_GP_3D_arr(var_name,file_name_i,z,y,x,draw)                 ! multiple plots
        use num_vars, only: rank
        
        ! input / output
        character(len=*), intent(in) :: var_name
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
        
        ! set nplt, npnt
        npntx = size(z,1)
        npnty = size(z,2)
        nplt = size(z,3)
        if (present(x)) then
            if (size(x,1).ne.size(z,1) .or. size(x,2).ne.size(z,2) .or. &
                &size(x,3).ne.size(z,3)) then
                call writo('WARNING: In print_GP, the size of x and z has to &
                    &be the same... Skipping plot',persistent=.true.)
                return
            end if
        end if
        if (present(y)) then
            if (size(y,1).ne.size(z,1) .or. size(y,2).ne.size(z,2) .or. &
                &size(y,3).ne.size(z,3)) then
                call writo('WARNING: In print_GP, the size of y and z has to &
                    &be the same... Skipping plot',persistent=.true.)
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
                &trim(i2str(rank))//'.dat'
        else
            file_name = trim(file_name_i)
        end if
        
        ! open output file
        open(unit=nextunit(file_i),file=data_dir//'/'//trim(file_name)//'.dat',&
            &iostat=istat)
        
        ! write to output file
        write(file_i,*) '# '//trim(var_name)//':'
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
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(var_name,file_name,file_name,nplt,2,.true.)
        
        if (trim(file_name_i).eq.'') then
            call use_execute_command_line('rm '//data_dir//'/'//&
                &trim(file_name)//'.dat')
        end if
    end subroutine print_GP_3D_arr
    
    ! use GNUPlot to draw a plot
    ! The variable file_name(s) holds the file(s)  which are to be plot. Each of
    ! them should contain nplt plots, arranged in columns. Furthermore, draw_dim
    ! determines  whether  the  plot  is  to  be 2D,  3D  or  decoupled  3D  and
    ! plot_on_screen whether the plot is be shown  on the screen, or to be saved
    ! in a file  called [var_name].pdf. Finally, for each of  the file names, an
    ! optional command draw_ops  can be provided, that specifies  the line style
    ! for the plots from the file.
    subroutine draw_GP_ind(var_name,file_name,draw_name,nplt,draw_dim,&
        &plot_on_screen,draw_ops,extra_ops)
        
        ! input / output
        character(len=*), intent(in) :: draw_name                               ! name of drawing
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: var_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        integer, intent(in) :: draw_dim                                         ! 1: 2D, 2: 3D, 3: decoupled 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        character(len=*), intent(in), optional :: draw_ops                      ! drawing option
        character(len=*), intent(in), optional :: extra_ops                     ! extra option
        
        ! call the array version
        if (present(draw_ops)) then
            call draw_GP_arr(var_name,[file_name],draw_name,nplt,draw_dim,&
                &plot_on_screen,draw_ops=[draw_ops],extra_ops=extra_ops)
        else
            call draw_GP_arr(var_name,[file_name],draw_name,nplt,draw_dim,&
                &plot_on_screen,extra_ops=extra_ops)
        end if
    end subroutine draw_GP_ind
    subroutine draw_GP_arr(var_name,file_names,draw_name,nplt,draw_dim,&
        &plot_on_screen,draw_ops,extra_ops)
        use num_vars, only: rank
        
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of function
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: draw_name                               ! name of drawing
        integer, intent(in) :: nplt                                             ! number of plots
        integer, intent(in) :: draw_dim                                         ! 1: 2D, 2: 3D, 3: decoupled 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        character(len=*), intent(in), optional :: draw_ops(:)                   ! extra drawing options (one per file)
        character(len=*), intent(in), optional :: extra_ops                     ! extra option
        
        ! local variables
        character(len=3*size(file_names)*max_str_ln) :: plot_cmd                ! individual plot command
        character(len=max_str_ln) :: script_name                                ! name of script, including path
        character(len=max_str_ln) :: loc_draw_op                                ! local draw option
        integer :: istat                                                        ! status of opening a file
        integer :: iplt, ifl                                                    ! counters
        integer :: nfl                                                          ! number of files to plot
        integer :: cmd_i                                                        ! file number for script file
        integer :: plt_count                                                    ! counts the number of plots
        integer :: cmdstat                                                      ! status of C system command
        character(len=5*max_str_ln) :: cmdmsg                                   ! error message of C system command
        
        ! set up nfl
        nfl = size(file_names)
        
        ! tests
        if (present(draw_ops)) then
            if (nfl.ne.size(draw_ops)) then
                call writo('WARNING: in draw_GP, one value for draw_ops has to &
                    &be provided for each file to plot',persistent=.true.)
                return
            end if
        end if
        
        ! create the GNUPlot command
        if (plot_on_screen) then
            script_name = trim(script_dir)//'/'//'temp_script_draw_GP_'//&
                &trim(i2str(rank))//'.gnu'
        else
            script_name = trim(script_dir)//'/'//trim(draw_name)//'.gnu'
        end if
        
        ! open script file
        open(unit=nextunit(cmd_i),file=trim(script_name),iostat=istat)
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command',&
                &persistent=.true.)
            return
        end if
        
        ! initialize the GNUPlot script
        if (plot_on_screen) then                                                ! wxt terminal
            write(cmd_i,*) 'set grid; set border 4095 front linetype -1 &
                &linewidth 1.000; set terminal wxt;'
        else                                                                    ! pdf terminal
            write(cmd_i,*) 'set grid; set border 4095 front linetype -1 &
                &linewidth 1.000; set terminal pdf size '//&
                &trim(i2str(plot_size(1)))//','//trim(i2str(plot_size(2)))//&
                &'; set output '''//trim(plot_dir)//'/'//trim(draw_name)//&
                &'.pdf'';'
        end if
        
        ! no legend if too many plots
        if (nplt.gt.10) write(cmd_i,*) 'set nokey;'
        
        ! write extra options
        if (present(extra_ops)) write(cmd_i,*) trim(extra_ops)
        
        ! set up line styles
        if (.not.present(draw_ops)) then
            do plt_count = 1,min(nplt*nfl,size(line_clrs))
                write(cmd_i,*) 'set style line '//trim(i2str(plt_count))//&
                    &' lc rgb '//line_clrs(plt_count)//' '//trim(line_style)
            end do
        end if
        
        ! set up plt_count
        plt_count = 1
        
        ! individual plots
        loc_draw_op = ''
        select case (draw_dim)
            case (1)                                                            ! 2D
                write(cmd_i,*) 'plot \'
                do iplt = 1,nplt
                    plot_cmd = ''
                    do ifl = 1,nfl
                        if (present(draw_ops)) then
                            loc_draw_op = trim(draw_ops(ifl))
                        else
                            loc_draw_op = 'with linespoints linestyle '//&
                                &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                        end if
                        plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                        &trim(file_names(ifl))//'.dat'' using '//&
                        &trim(i2str(iplt))//':'//trim(i2str(nplt+iplt))//&
                        &' title '''//trim(var_name)//' ('//trim(i2str(iplt))//&
                        &'/'//trim(i2str(nplt))//')'' '//trim(loc_draw_op)//','
                        if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
                        plt_count = plt_count + 1
                    end do
                    write(cmd_i,*) trim(plot_cmd)
                end do
            case (2)                                                            ! 3D
                write(cmd_i,*) 'splot \'
                do iplt = 1,nplt
                    plot_cmd = ''
                    do ifl = 1,nfl
                        if (present(draw_ops)) then
                            loc_draw_op = trim(draw_ops(ifl))
                        else
                            loc_draw_op = 'with linespoints linestyle '//&
                                &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                        end if
                        plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                        &trim(file_names(ifl))//'.dat'' using '//&
                        &trim(i2str(iplt))//':'//trim(i2str(nplt+iplt))//':'//&
                        &trim(i2str(2*nplt+iplt))//' title '''//trim(var_name)&
                        &//' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                        &')'' '//trim(loc_draw_op)//','
                        if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
                        plt_count = plt_count + 1
                    end do
                    write(cmd_i,*) trim(plot_cmd)
                end do
            case (3)
                write(cmd_i,*) 'splot \'
                do iplt = 1,nplt
                    plot_cmd = ''
                    do ifl = 1,nfl
                        if (present(draw_ops)) then
                            loc_draw_op = trim(draw_ops(ifl))
                        else
                            loc_draw_op = 'with linespoints linestyle '//&
                                &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                        end if
                        plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                        &trim(file_names(ifl))//'.dat'' using ('//&
                        &trim(i2str(iplt))//'):'//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//' title '''//trim(var_name)//&
                        &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                        &')'' '//trim(loc_draw_op)//','
                        if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
                        plt_count = plt_count + 1
                    end do
                    write(cmd_i,*) trim(plot_cmd)
                end do
            case default
                call writo('No draw_dim associated with '//&
                    &trim(i2str(draw_dim)),persistent=.true.)
                istat = 1
                CHCKSTT
        end select
        
        ! finishing the GNUPlot command
        write(cmd_i,*) ''
        if (plot_on_screen) write(cmd_i,*) 'pause -1'
        
        ! closing the GNUPlot command
        close(cmd_i)
        
        ! stop timer
        call stop_time
        
        ! bypass plots if no_plots
        if (.not.no_plots) then
            ! call GNUPlot
            call use_execute_command_line('gnuplot "'//trim(script_name)//'"'//&
                &err_output_str,exitstat=istat,cmdstat=cmdstat,cmdmsg=cmdmsg)
            
            if (istat.ne.0) then
                call writo('Failed to plot '//trim(draw_name)//'.pdf',&
                    &persistent=.true.)
            else
                if (cmdstat.ne.0) then
                    call writo('Failed to plot '//trim(draw_name)//'.pdf',&
                        &persistent=.true.)
                    call lvl_ud(1)
                    if (cmdstat.ne.0) call writo('System message: "'//&
                        &trim(cmdmsg)//'"',persistent=.true.)
                    if (.not.plot_on_screen) call writo(&
                        &'Try running "gnuplot "'//trim(script_name)//'"'//&
                        &'" manually',persistent=.true.)
                    call lvl_ud(-1)
                else
                    if (plot_on_screen) then
                        call use_execute_command_line('rm '//trim(script_name),&
                            &exitstat=istat,cmdstat=cmdstat,cmdmsg=cmdmsg)
                        ! ignore errors
                    else
                        call writo('Created plot in output file '''//&
                            &trim(plot_dir)//'/'//trim(draw_name)//'.pdf''',&
                            &persistent=.true.)
                    end if
                end if
            end if
        end if
        
        ! start timer
        call start_time
    end subroutine draw_GP_arr
    
    ! The variable file_name(s) holds the file(s)  which are to be plot. Each of
    ! them should contain nplt plots, arranged in columns. Furthermore, draw_dim
    ! determines whether the plot  is to be 2D, 3D or decoupled  3D. The plot is
    ! saved in  a file  called [var_name].pdf.  For each of  the file  names, an
    ! optional command draw_ops  can be provided, that specifies  the line style
    ! for  the plots  from the  file. Finally,  also optional  commands for  the
    ! ranges of the figure and the delay between the frames can be specified.
    subroutine draw_GP_animated_ind(var_name,file_name,draw_name,nplt,draw_dim,&
        &ranges,delay,draw_ops,extra_ops)
        
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of function
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: draw_name                               ! name of drawing
        integer, intent(in) :: nplt                                             ! number of plots
        integer, intent(in) :: draw_dim                                         ! 1: 2D, 2: 3D, 3: decoupled 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        integer, intent(in), optional :: delay                                  ! time delay between plots
        character(len=*), intent(in), optional :: draw_ops                      ! extra commands
        character(len=*), intent(in), optional :: extra_ops                     ! extra option
        
        ! call the array version
        if (present(draw_ops)) then
            call draw_GP_animated_arr(var_name,[file_name],draw_name,nplt,&
                &draw_dim,ranges=ranges,delay=delay,draw_ops=[draw_ops],&
                &extra_ops=extra_ops)
        else
            call draw_GP_animated_arr(var_name,[file_name],draw_name,nplt,&
                &draw_dim,delay=delay,ranges=ranges,extra_ops=extra_ops)
        end if
    end subroutine draw_GP_animated_ind
    subroutine draw_GP_animated_arr(var_name,file_names,draw_name,nplt,&
        &draw_dim,ranges,delay,draw_ops,extra_ops)
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of function
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: draw_name                               ! name of drawing
        integer, intent(in) :: nplt                                             ! number of plots
        integer, intent(in) :: draw_dim                                         ! 1: 2D, 2: 3D, 3: decoupled 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        integer, intent(in), optional :: delay                                  ! time delay between plots
        character(len=*), intent(in), optional :: draw_ops(:)                   ! extra commands
        character(len=*), intent(in), optional :: extra_ops                     ! extra option
        
        ! local variables
        character(len=3*size(file_names)*max_str_ln) :: plot_cmd                ! individual plot command
        character(len=max_str_ln) :: script_name                                ! name of script, including directory
        character(len=max_str_ln) :: loc_draw_op                                ! local draw option
        integer :: istat                                                        ! status of opening a file
        integer :: iplt, ifl                                                    ! counters
        integer :: nfl                                                          ! number of files to plot
        integer :: cmd_i                                                        ! file number for script file
        integer :: delay_loc                                                    ! local copy of delay
        real(dp), allocatable :: ranges_loc(:,:)                                ! local copy of ranges
        integer :: plt_count                                                    ! counts the number of plots
        integer :: cmdstat                                                      ! status of C system command
        character(len=5*max_str_ln) :: cmdmsg                                   ! error message of C system command
        
        ! set up nfl
        nfl = size(file_names)
        
        ! tests
        if (present(draw_ops)) then
            if (nfl.ne.size(draw_ops)) then
                call writo('WARNING: in draw_GP_animated, one value for &
                    &draw_ops has to be provided for each file to plot',&
                    &persistent=.true.)
                return
            end if
        end if
        
        ! calculate delay [1/100 s]
        if (present(delay)) then
            delay_loc = delay
        else
            delay_loc = max(250/nplt,10)
        end if
        
        ! create the GNUPlot command
        script_name = ''//trim(script_dir)//'/'//trim(draw_name)//'.gnu'
        
        ! open script file
        open(unit=nextunit(cmd_i),file=trim(script_name),iostat=istat)
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command',&
                &persistent=.true.)
            return
        end if
        
        ! initialize the GNUPlot script
        write(cmd_i,*) 'set grid; set border 4095 front &
            &linetype -1 linewidth 1.000; set terminal gif animate delay '//&
            &trim(i2str(delay_loc))//' size 1280, 720; &
            &set output '''//trim(plot_dir)//'/'//trim(draw_name)//&
            &'.gif'';'
        
        ! find and set ranges
        select case (draw_dim)
            case (1,3)                                                          ! 2D or decoupled 3D
                ! find ranges
                allocate(ranges_loc(2,2))
                if (present(ranges)) then
                    if (size(ranges,1).eq.2 .and. size(ranges,2).eq.2) then
                        ranges_loc = ranges
                    else
                        call writo('WARNING: invalid ranges given to &
                            &draw_GP_animated',persistent=.true.)
                    end if
                else
                    ! initialize ranges
                    ranges_loc(:,1) = huge(1._dp)                               ! minimum value
                    ranges_loc(:,2) = -huge(1._dp)                              ! maximum value
                    
                    do ifl = 1,nfl
                        call get_ranges(file_names(ifl),ranges_loc)
                    end do
                end if
                
                ! set ranges
                write(cmd_i,*) 'set xrange ['//trim(r2str(ranges_loc(1,1)))//&
                    &':'//trim(r2str(ranges_loc(1,2)))//'];'
                write(cmd_i,*) 'set yrange ['//trim(r2str(ranges_loc(2,1)))//&
                    &':'//trim(r2str(ranges_loc(2,2)))//'];'
            case (2)
                ! find ranges
                allocate(ranges_loc(3,2))
                if (present(ranges)) then
                    if (size(ranges,1).eq.3 .and. size(ranges,2).eq.2) then
                        ranges_loc = ranges
                    else
                        call writo('WARNING: invalid ranges given to &
                            &draw_GP_animated',persistent=.true.)
                    end if
                else
                    ranges_loc(:,1) = 1.E14_dp                                  ! minimum value
                    ranges_loc(:,2) = -1.E14_dp                                 ! maximum value
                    
                    do ifl = 1,nfl
                        call get_ranges(file_names(ifl),ranges_loc)
                    end do
                end if
                
                ! set ranges
                write(cmd_i,*) ' set xrange ['//trim(r2str(ranges_loc(1,1)))//&
                    &':'//trim(r2str(ranges_loc(1,2)))//'];'
                ! y range
                write(cmd_i,*) ' set yrange ['//trim(r2str(ranges_loc(2,1)))//&
                    &':'//trim(r2str(ranges_loc(2,2)))//'];'
                ! z range
                write(cmd_i,*) ' set zrange ['//trim(r2str(ranges_loc(3,1)))//&
                    &':'//trim(r2str(ranges_loc(3,2)))//'];'
                ! color scale
                write(cmd_i,*) ' set cbrange ['//trim(r2str(ranges_loc(3,1)))//&
                    &':'//trim(r2str(ranges_loc(3,2)))//'];'
                
                ! other definitions for 3D plots
                write(cmd_i,*) 'set view 45,45;'
                write(cmd_i,*) 'set hidden3d offset 0'
            case default
                call writo('No draw_dim associated with '//&
                    &trim(i2str(draw_dim)),persistent=.true.)
                istat = 1
                CHCKSTT
        end select
        
        ! no legend if too many plots
        !if (nfl.gt.10) write(cmd_i,*) 'set nokey;'
        
        ! write extra options
        if (present(extra_ops)) write(cmd_i,*) trim(extra_ops)
        
        ! set up line styles
        if (.not.present(draw_ops)) then
            do plt_count = 1,min(nplt*nfl,size(line_clrs))
                write(cmd_i,*) 'set style line '//trim(i2str(plt_count))//&
                    &' lc rgb '//line_clrs(plt_count)//' '//trim(line_style)
            end do
        end if
        
        ! set up plt_count
        plt_count = 1
        
        ! individual plots
        loc_draw_op = ''
        select case (draw_dim)
            case (1)                                                            ! 2D
                do iplt = 1,nplt
                    plot_cmd = 'plot'
                    do ifl = 1,nfl
                        if (present(draw_ops)) then
                            loc_draw_op = trim(draw_ops(ifl))
                        else
                            loc_draw_op = 'with linespoints linestyle '//&
                                &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                        end if
                        plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                        &trim(file_names(ifl))//'.dat'' using '//&
                        &trim(i2str(iplt))//':'//trim(i2str(nplt+iplt))//&
                        &' title '''//trim(var_name)//' ('//trim(i2str(iplt))//&
                        &'/'//trim(i2str(nplt))//')'' '//trim(loc_draw_op)
                        if (ifl.ne.nfl) plot_cmd = trim(plot_cmd)//', '
                        plt_count = plt_count + 1
                    end do
                    write(cmd_i,*) trim(plot_cmd)
                end do
            case (2)                                                            ! 2D
                do iplt = 1,nplt
                    plot_cmd = 'splot '
                    do ifl = 1,nfl
                        if (present(draw_ops)) then
                            loc_draw_op = trim(draw_ops(ifl))
                        else
                            loc_draw_op = 'with linespoints linestyle '//&
                                &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                        end if
                        plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                        &trim(file_names(ifl))//'.dat'' using '//&
                        &trim(i2str(iplt))//':'//trim(i2str(nplt+iplt))//':'//&
                        &trim(i2str(2*nplt+iplt))//' title '''//trim(var_name)&
                        &//' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                        &')'' '//trim(loc_draw_op)
                        if (ifl.ne.nfl) plot_cmd = trim(plot_cmd)//', '
                        plt_count = plt_count + 1
                    end do
                    write(cmd_i,*) trim(plot_cmd)
                end do
            case (3)
                do iplt = 1,nplt
                    plot_cmd = 'splot'
                    do ifl = 1,nfl
                        if (present(draw_ops)) then
                            loc_draw_op = trim(draw_ops(ifl))
                        else
                            loc_draw_op = 'with linespoints linestyle '//&
                                &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                        end if
                        plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                        &trim(file_names(ifl))//'.dat'' using ('//&
                        &trim(i2str(iplt))//'):'//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//' title '''//trim(var_name)//&
                        &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                        &')'' '//trim(loc_draw_op)
                        if (ifl.ne.nfl) plot_cmd = trim(plot_cmd)//', '
                        plt_count = plt_count + 1
                    end do
                    write(cmd_i,*) trim(plot_cmd)
                end do
            case default
                call writo('No draw_dim associated with '//&
                    &trim(i2str(draw_dim)),persistent=.true.)
                istat = 1
                CHCKSTT
        end select
        
        ! finishing and closing the GNUPlot command
        write(cmd_i,*) 'set output'
        close(cmd_i)
        
        ! no need to stop the time
        
        ! bypass plots if no_plots
        if (.not.no_plots) then
            ! call GNUPlot
            call use_execute_command_line('gnuplot "'//trim(script_name)//'"'//&
                &err_output_str,exitstat=istat,cmdstat=cmdstat,&
                &cmdmsg=cmdmsg)
            
            if (istat.ne.0) then
                call writo('Failed to plot '//trim(draw_name)//'.pdf',&
                    &persistent=.true.)
            else
                if (cmdstat.ne.0) then
                    call writo('Failed to plot '//trim(draw_name)//'.pdf',&
                    &persistent=.true.)
                    call lvl_ud(1)
                    if (cmdstat.ne.0) call writo('System message: "'//&
                        &trim(cmdmsg)//'"',persistent=.true.)
                    call lvl_ud(-1)
                else
                    call writo('Created animated plot in output file '''//&
                        &trim(plot_dir)//'/'//trim(draw_name)//'.pdf''',&
                        &persistent=.true.)
                end if
            end if
        end if
    contains
        subroutine get_ranges(file_name,ranges)
            ! input / output
            character(len=*), intent(in) :: file_name                           ! name of file
            real(dp), intent(inout) :: ranges(:,:)                              ! x and y range, and z range (if 3D) of plot
            
            ! local variables
            integer :: data_i                                                   ! file number of data file
            real(dp), allocatable :: loc_data(:)                                ! one line of data
            character(len=1) :: loc_data_char
            
            
            ! open data file
            open(unit=nextunit(data_i),file=data_dir//'/'//trim(file_name),&
                &status='old',iostat=istat)
            
            ! set up loc_data
            select case (draw_dim)
                case (1,3)                                                          ! 2D or decoupled 3D
                    allocate(loc_data(2*nplt))
                case (2)
                    allocate(loc_data(3*nplt))
                case default
                    call writo('No draw_dim associated with '//&
                        &trim(i2str(draw_dim)),persistent=.true.)
                    istat = 1
                    CHCKSTT
            end select
            
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
                        if (draw_dim.eq.2) then
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
    end subroutine draw_GP_animated_arr
    
    ! merges GNUPlot data files
    subroutine merge_GP(file_names_in,file_name_out,delete)
        ! input / output
        character(len=*), intent(in) :: file_names_in(:)                        ! name of input file
        character(len=*), intent(in) :: file_name_out                           ! name of output file
        logical, intent(in), optional :: delete                                 ! .true. if input files are to be deleted
        
        ! local variables
        integer :: istat                                                        ! status
        character(len=size(file_names_in)*max_str_ln) :: shell_cmd              ! shell command
        integer :: id                                                           ! counter
        logical :: delete_loc                                                   ! local copy of delete
        
        ! initialize istat
        istat = 0
        
        ! set up concatenation command
        shell_cmd = 'cat'
        do id = 1,size(file_names_in)
            shell_cmd = trim(shell_cmd)//' "'//trim(data_dir)//'/'//&
                &trim(file_names_in(id))//'"'
        end do
        shell_cmd = trim(shell_cmd)//' > '//trim(data_dir)//'/'//&
            &trim(file_name_out)
        
        ! concatenate using shell
        call use_execute_command_line(trim(shell_cmd),exitstat=istat)
        if (istat.ne.0) then
            call writo('WARNING: merge_GP failed to merge',persistent=.true.)
            return
        end if
        
        ! set up delete
        if (present(delete)) then
            delete_loc = delete
        else
            delete_loc = .false.
        end if
        
        if (delete_loc) then
            ! set up delete command
            shell_cmd = 'rm -f'
            do id = 1,size(file_names_in)
                shell_cmd = trim(shell_cmd)//' "'//trim(data_dir)//'/'//&
                    &trim(file_names_in(id))//'"'
            end do
            
            ! delete using shell
            call use_execute_command_line(trim(shell_cmd),exitstat=istat)
            if (istat.ne.0) then
                call writo('WARNING: merge_GP failed to delete',&
                &persistent=.true.)
                return
            end if
        end if
    end subroutine

    ! Prints variables  "vars" with names "var_names"  in a HDF5 file  with name
    ! "file_name" and  accompanying XDMF file.  For collections, only  the first
    ! value for var_names is used, so it should have a size of one.
    ! The plot is generally  3D, but if one of the  dimensions provided is equal
    ! to 1, it is checked whether  there is poloidal or toroidal axisymmetry and
    ! if so, the plot becomes 2D.
    ! Optionally, the  (curvilinear) grid can  be provided through "X",  "Y" and
    ! "Z". If  not, the grid  is assumed to  be Cartesian with  discrete indices
    ! where X corresponds to the first dimensions,  Y to the second and Z to the
    ! third.
    ! Additionally, the  total grid  size and  local offset  can be  provided in
    ! "tot_dim" and "loc_offset" to run this routine in parallel automatically.
    ! Optionally, one of the dimensions (col_id, default 4) can be associated to
    ! a collection dimension using "col" different from 1:
    !   col = 1: no collection, just plots of different variables
    !   col = 2: time collection
    !   col = 3: spatial collection
    ! Note: To plot this with VisIt, use:
    !   - for temporal collections: pseudocolor using the variable name (other 
    !       names are ignored).
    !   - for spatial collections:  subset of blocks or pseudocolor using the 
    !       variable name (other names are ignored).
    !   - for without collections:  pseudocolor using the variable names.
    ! Note: Currently all of possibly multiple processes that write data have to
    ! cover the  entire range  of the  variables in  the dimension  indicated by
    ! col_id. (This could  be implemented by changing how n_plot  is defined and
    ! selectively  letting  each  processer  write  in  the  main  loop  at  its
    ! corresponding indices.)
    ! Note:  This  routine is  written  with  toroidal configurations  in  mind.
    ! Therefore, if X or Y is negative,  and there is poloidal symmetry, this is
    ! ignored and X and Y are combined into R = sqrt(X^2+Y^2).
    subroutine plot_HDF5_arr(var_names,file_name,vars,tot_dim,loc_offset,&
        &X,Y,Z,col_id,col,description)                                          ! array version
        use HDF5_ops, only: open_HDF5_file, add_HDF5_item, print_HDF5_top, &
            &print_HDF5_geom, print_HDF5_3D_data_item, print_HDF5_att, &
            &print_HDF5_grid, close_HDF5_file
        use HDF5_vars, only: XML_str_type, HDF5_file_type
        use num_vars, only: n_procs
        use MPI_utilities, only: get_ser_var
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            ! names of variable to be plot
        character(len=*), intent(in) :: file_name                               ! file name
        real(dp), intent(in), target :: vars(:,:,:,:)                           ! variables to plot
        integer, intent(in), optional :: tot_dim(4)                             ! total dimensions of the arrays
        integer, intent(in), optional :: loc_offset(4)                          ! offset of local dimensions
        real(dp), intent(in), target, optional :: X(:,:,:,:)                    ! curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:,:)                    ! curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:,:)                    ! curvlinear grid Z points
        integer, intent(in), optional :: col_id                                 ! index of time dimension
        integer, intent(in), optional :: col                                    ! whether a collection is made
        character(len=*), intent(in), optional :: description                   ! description
        
        ! local variables
        type(HDF5_file_type) :: file_info                                       ! file info
        integer :: istat                                                        ! status
        integer :: col_id_loc                                                   ! local copy of col_id
        integer :: col_loc                                                      ! local copy of col
        integer :: n_plot                                                       ! nr. of plots
        integer :: id, jd                                                       ! counter
        integer :: sym_type                                                     ! type of symmetry (1: no symmetry, 2: toroidal, 3: poloidal)
        integer :: sym_pol, sym_tor                                             ! used to determine symmetry
        integer, allocatable :: tot_sym_pol(:), tot_sym_tor(:)                  ! sym_pol and sym_tor for all processes
        integer :: tot_dim_loc(4)                                               ! local copy of tot_dim
        integer :: loc_offset_loc(4)                                            ! local copy of loc_offset
        integer :: tot_dim_3D(3)                                                ! tot_dim except collection
        integer :: loc_dim_3D(3)                                                ! loc_dim except collection
        integer :: loc_offset_3D(3)                                             ! loc_offset except collection
        type(XML_str_type) :: col_grid                                          ! grid with collection
        type(XML_str_type), allocatable :: grids(:)                             ! the grids in the time collection
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type), allocatable :: XYZ(:)                               ! data items for geometry
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        logical :: col_mask(4)                                                  ! to select out the collection dimension
        logical :: ind_plot                                                     ! individual plot
        real(dp), allocatable :: sym_ang(:,:,:)                                 ! angle to be checked for symmetry
        real(dp) :: tol_sym = 1.E-8_dp                                          ! tolerance for symmetry determination
        real(dp), pointer :: var_3D(:,:,:) => null()                            ! pointer to vars
        real(dp), pointer :: X_3D(:,:,:) => null()                              ! pointer to X
        real(dp), pointer :: Y_3D(:,:,:) => null()                              ! pointer to Y 
        real(dp), pointer :: Z_3D(:,:,:) => null()                              ! pointer to Z
        character(len=max_str_ln), allocatable :: grd_names(:)                  ! grid names
        character(len=max_str_ln), allocatable :: att_names(:)                  ! attribute names
        
        ! set up local col_id and col
        col_id_loc = 4                                                          ! default collection dimension: last index
        if (present(col_id)) col_id_loc = col_id
        col_loc = 1                                                             ! default no spatial collection
        if (present(col)) col_loc = col
        
        ! set up nr. of plots
        n_plot = size(vars,col_id_loc)
        
        ! set up local tot_dim and loc_offset
        tot_dim_loc = shape(vars)
        if (present(tot_dim)) tot_dim_loc = tot_dim
        loc_offset_loc = 0
        if (present(loc_offset)) loc_offset_loc = loc_offset
        
        ! tests
        if (tot_dim_loc(col_id_loc).ne.n_plot) then
            istat = 1
            call writo('WARNING: In plot_HDF5, all the processes need to have &
                &the full range in the dimension given by col_id',&
                &persistent=.true.)
            CHCKSTT
        end if
        if (n_plot.eq.1 .and. col_loc.ne.1) then
            istat = 1
            call writo('WARNING: In plot_HDF5, if single plot, the collection &
                &type needs to be one',persistent=.true.)
            CHCKSTT
        end if
        
        ! set 3D dimensions
        col_mask = .false.
        col_mask(col_id_loc) = .true.
        tot_dim_3D = pack(tot_dim_loc,.not.col_mask)
        loc_dim_3D = pack(shape(vars),.not.col_mask)
        loc_offset_3D = pack(loc_offset_loc,.not.col_mask)
        
        ! set up individual plot
        if (tot_dim_3D(1).eq.loc_dim_3D(1) .and. &
            &tot_dim_3D(2).eq.loc_dim_3D(2) .and. &
            &tot_dim_3D(3).eq.loc_dim_3D(3)) then
            ind_plot = .true.
        else
            ind_plot = .false.
        end if
        
        ! default symmetry type
        sym_type = 1
        
        ! Find  symmetry type  by  checking whether  Y/X  is constant  (toroidal
        ! symmetry) or  Z^2/(X^2+Y^2) is  constant (poloidal symmetry),  for all
        ! plots.
        if (minval(tot_dim_3D).eq.1) then                                       ! possibly symmetry
            ! allocate helper variable
            allocate(sym_ang(loc_dim_3D(1),loc_dim_3D(2),loc_dim_3D(3)))
            ! initialize sym_pol and sym_tor
            sym_pol = 0
            sym_tor = 0
            ! loop over all plots
            do id = 1,n_plot
                ! assign pointers
                call assign_pointers(id)
                ! check poloidal angle
                sym_ang = atan(Y_3D/X_3D)
                if (maxval(sym_ang)-minval(sym_ang).lt.tol_sym) &
                    &sym_pol = sym_pol+1                                        ! poloidal symmetry for this plot
                ! check toroidal angle
                sym_ang = atan(sqrt(Z_3D**2/(X_3D**2+Y_3D**2)))
                if (maxval(sym_ang)-minval(sym_ang).lt.tol_sym) &
                    &sym_tor = sym_tor+1                                        ! toroidal symmetry for this plot
            end do
            ! get total sym_pol and sym_tor
            if (ind_plot) then                                                  ! so that below test succeeds
                allocate(tot_sym_pol(n_procs),tot_sym_tor(n_procs))
                tot_sym_pol = sym_pol
                tot_sym_tor = sym_tor
            else                                                                ! get from all the processes
                istat = get_ser_var([sym_pol],tot_sym_pol,scatter=.true.)
                CHCKSTT
                istat = get_ser_var([sym_tor],tot_sym_tor,scatter=.true.)
                CHCKSTT
            end if
            
            ! check total results
            if (sum(tot_sym_pol).eq.n_procs*n_plot) then                        ! poloidal symmetry for all plots
                sym_type = 2
            else if (sum(tot_sym_tor).eq.n_procs*n_plot) then                   ! toroidal symmetry for all plots
                sym_type = 3
            end if
            ! deallocate helper variables
            deallocate(sym_ang)
        end if
        
        ! set up local var_names
        allocate(grd_names(n_plot))
        allocate(att_names(n_plot))
        if (col_loc.eq.1) then                                                  ! without collection
            if (n_plot.eq.1) then                                               ! just one plot: attribute name is important
                if (size(var_names).eq.n_plot) then                             ! the right number of variable names provided
                    att_names = var_names
                else if (size(var_names).gt.n_plot) then                        ! too many variable names provided
                    att_names = var_names(1:n_plot)
                    call writo('WARNING: Too many variable names provided',&
                &persistent=.true.)
                else                                                            ! not enough variable names provided
                    att_names(1:size(var_names)) = var_names
                    do id = size(var_names)+1,n_plot
                        att_names(id) = 'unnamed variable '//trim(i2str(id))
                    end do
                    call writo('WARNING: Not enough variable names provided',&
                &persistent=.true.)
                end if
                grd_names = 'default_grid_name'
            else                                                                ! multiple plots: grid name is important
                if (size(var_names).eq.n_plot) then                             ! the right number of variable names provided
                    grd_names = var_names
                else if (size(var_names).gt.n_plot) then                        ! too many variable names provided
                    grd_names = var_names(1:n_plot)
                    call writo('WARNING: Too many variable names provided',&
                &persistent=.true.)
                else                                                            ! not enough variable names provided
                    grd_names(1:size(var_names)) = var_names
                    do id = size(var_names)+1,n_plot
                        grd_names(id) = 'unnamed variable '//trim(i2str(id))
                    end do
                    call writo('WARNING: Not enough variable names provided',&
                &persistent=.true.)
                end if
                att_names = 'default_att_name'
            end if
        else                                                                    ! collections: attribute name is important
            att_names = var_names(1)
            if (size(var_names).gt.1) call writo('WARNING: For collections, &
                &only the first variable name is used',persistent=.true.)
            grd_names = 'default_grid_name'
        end if
        
        ! open HDF5 file
        istat = open_HDF5_file(file_info,file_name,description,&
            &ind_plot=ind_plot)
        CHCKSTT
        
        ! create grid for collection
        allocate(grids(n_plot))
        
        ! allocate geometry arrays
        if (sym_type.eq.1) then                                                 ! 3D grid
            allocate(XYZ(3))
        else                                                                    ! 2D grid
            allocate(XYZ(2))
        end if
        
        ! loop over all plots
        do id = 1,n_plot
            ! print topology
            if (sym_type.eq.1) then                                             ! 3D grid
                call print_HDF5_top(top,2,tot_dim_3D,ind_plot=ind_plot)
            else                                                                ! 2D grid
                call print_HDF5_top(top,1,tot_dim_3D,ind_plot=ind_plot)
            end if
            
            ! assign pointers
            call assign_pointers(id)
            
            ! print data  item for X, Y  and Z (no symmetry),  R = sqrt(X^2+Y^2)
            ! and Z (poloidal symmetry) or X and Y (toroidal symmetry)
            select case (sym_type)
                case (1)                                                        ! no symmetry
                    istat = print_HDF5_3D_data_item(XYZ(1),file_info,&
                        &'X_'//trim(i2str(id)),X_3D,tot_dim_3D,loc_dim_3D,&
                        &loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                    istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                        &'Y_'//trim(i2str(id)),Y_3D,tot_dim_3D,loc_dim_3D,&
                        &loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                    istat = print_HDF5_3D_data_item(XYZ(3),file_info,&
                        &'Z_'//trim(i2str(id)),Z_3D,tot_dim_3D,loc_dim_3D,&
                        &loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                case (2)                                                        ! poloidal symmetry
                    istat = print_HDF5_3D_data_item(XYZ(1),file_info,&
                        &'R_'//trim(i2str(id)),sqrt(X_3D**2+Y_3D**2),&
                        &tot_dim_3D,loc_dim_3D,loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                    istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                        &'Z_'//trim(i2str(id)),Z_3D,tot_dim_3D,loc_dim_3D,&
                        &loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                case (3)                                                        ! toroidal symmetry
                    istat = print_HDF5_3D_data_item(XYZ(1),file_info,&
                        &'X_'//trim(i2str(id)),X_3D,tot_dim_3D,loc_dim_3D,&
                        &loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                    istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                        &'Y_'//trim(i2str(id)),Y_3D,tot_dim_3D,loc_dim_3D,&
                        &loc_offset_3D,ind_plot=ind_plot)
                    CHCKSTT
                case default                                                    ! no symmetry
                    istat = 1
                    call writo('WARNING: symmetry type '//&
                        &trim(i2str(sym_type))//' not recognized',&
                        &persistent=.true.)
                    CHCKSTT
            end select
            
            ! print geometry with X, Y and Z data item
            if (sym_type.eq.1) then                                             ! no symmetry so 3D geometry
                call print_HDF5_geom(geom,2,XYZ,reset=.true.,ind_plot=ind_plot)
            else                                                                ! symmetry so 2D geometry
                call print_HDF5_geom(geom,1,XYZ,reset=.true.,ind_plot=ind_plot)
            end if
            
            ! print data item for plot variable
            istat = print_HDF5_3D_data_item(XYZ(1),file_info,'var_'//&
                &trim(i2str(id)),var_3D,tot_dim_3D,loc_dim_3D,loc_offset_3D,&
                &ind_plot=ind_plot)                                             ! reuse XYZ(1)
            CHCKSTT
            
            ! print attribute with this data item
            call print_HDF5_att(att(1),XYZ(1),att_names(id),1,reset=.true.,&
                &ind_plot=ind_plot)
            
            ! create a  grid with the topology, the geometry,  the attribute and
            ! time if time collection
            if (col_loc.eq.2) then                                              ! time collection
                istat = print_HDF5_grid(grids(id),grd_names(id),1,&
                    &grid_time=id*1._dp,grid_top=top,grid_geom=geom,&
                    &grid_atts=att,reset=.true.,ind_plot=ind_plot)
                CHCKSTT
            else                                                                ! no time collection
                istat = print_HDF5_grid(grids(id),grd_names(id),1,&
                    &grid_top=top,grid_geom=geom,grid_atts=att,reset=.true.,&
                    &ind_plot=ind_plot)
                CHCKSTT
            end if
        end do
        
        ! either create collection or just use individual grids
        if (col_loc.eq.1) then
            ! add individual grids to HDF5 file and reset them
            do id = 1,n_plot
                call add_HDF5_item(file_info,grids(id),reset=.true.,&
                &ind_plot=ind_plot)
            end do
        else
            ! create grid collection from individual grids and reset them
            istat = print_HDF5_grid(col_grid,'domain of mesh',col_loc,&
                &grid_grids=grids,reset=.true.,ind_plot=ind_plot)
            CHCKSTT
            
            ! add collection grid to HDF5 file and reset it
            call add_HDF5_item(file_info,col_grid,reset=.true.,&
                &ind_plot=ind_plot)
        end if
        
        ! close HDF5 file
        istat = close_HDF5_file(file_info,ind_plot=ind_plot)
        CHCKSTT
        
        ! clean up
        nullify(var_3D)
        nullify(X_3D,Y_3D,Z_3D)
    contains
        ! assigns the 3D subarray pointer variables
        subroutine assign_pointers(id)
            ! input / output
            integer :: id                                                       ! index at which to assing pointer
            
            ! X
            if (present(X)) then
                if (col_id_loc.eq.1) then
                    X_3D => X(id,:,:,:)
                else if (col_id_loc.eq.2) then
                    X_3D => X(:,id,:,:)
                else if (col_id_loc.eq.3) then
                    X_3D => X(:,:,id,:)
                else if (col_id_loc.eq.4) then
                    X_3D => X(:,:,:,id)
                else
                    istat = 1
                    CHCKSTT
                end if
            else
                allocate(X_3D(loc_dim_3D(1),loc_dim_3D(2),loc_dim_3D(3)))
                do jd = 1,loc_dim_3D(1)
                    X_3D(jd,:,:) = loc_offset_3D(1) + jd - 1
                end do
            end if
            ! Y
            if (present(Y)) then
                if (col_id_loc.eq.1) then
                    Y_3D => Y(id,:,:,:)
                else if (col_id_loc.eq.2) then
                    Y_3D => Y(:,id,:,:)
                else if (col_id_loc.eq.3) then
                    Y_3D => Y(:,:,id,:)
                else if (col_id_loc.eq.4) then
                    Y_3D => Y(:,:,:,id)
                else
                    istat = 1
                    CHCKSTT
                end if
            else
                allocate(Y_3D(loc_dim_3D(1),loc_dim_3D(2),loc_dim_3D(3)))
                do jd = 1,loc_dim_3D(2)
                    Y_3D(:,jd,:) = loc_offset_3D(2) + jd - 1
                end do
            end if
            ! X
            if (present(Z)) then
                if (col_id_loc.eq.1) then
                    Z_3D => Z(id,:,:,:)
                else if (col_id_loc.eq.2) then
                    Z_3D => Z(:,id,:,:)
                else if (col_id_loc.eq.3) then
                    Z_3D => Z(:,:,id,:)
                else if (col_id_loc.eq.4) then
                    Z_3D => Z(:,:,:,id)
                else
                    istat = 1
                    CHCKSTT
                end if
            else
                allocate(Z_3D(loc_dim_3D(1),loc_dim_3D(2),loc_dim_3D(3)))
                do jd = 1,loc_dim_3D(3)
                    Z_3D(:,:,jd) = loc_offset_3D(3) + jd - 1
                end do
            end if
            ! variable
            if (col_id_loc.eq.1) then
                var_3D => vars(id,:,:,:)
            else if (col_id_loc.eq.2) then
                var_3D => vars(:,id,:,:)
            else if (col_id_loc.eq.3) then
                var_3D => vars(:,:,id,:)
            else if (col_id_loc.eq.4) then
                var_3D => vars(:,:,:,id)
            else
                istat = 1
                CHCKSTT
            end if
        end subroutine assign_pointers
    end subroutine plot_HDF5_arr
    subroutine plot_HDF5_ind(var_name,file_name,var,tot_dim,loc_offset,&
        &X,Y,Z,description)                                                     ! individual version
        
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of variable to be plot
        character(len=*), intent(in) :: file_name                               ! file name
        real(dp), intent(in) :: var(:,:,:)                                      ! variable to plot
        integer, intent(in), optional :: tot_dim(3)                             ! total dimensions of the arrays
        integer, intent(in), optional :: loc_offset(3)                          ! offset of local dimensions
        real(dp), intent(in), target, optional :: X(:,:,:)                      ! curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:)                      ! curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:)                      ! curvlinear grid Z points
        character(len=*), intent(in), optional :: description                   ! description
        
        ! local variables
        integer :: tot_dim_loc(4), loc_offset_loc(4)                            ! local versions of total dimensions and local offset
        
        ! set up local tot_dim and loc_offset
        tot_dim_loc = [shape(var),1]
        if (present(tot_dim)) tot_dim_loc = [tot_dim,1]
        loc_offset_loc = 0
        if (present(loc_offset)) loc_offset_loc = [loc_offset,0]
        
        if (present(X)) then
            if (present(Y)) then
                if (present(Z)) then
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &Y=reshape(Y,[size(Y,1),size(Y,2),size(Y,3),1]),&
                        &Z=reshape(Z,[size(Z,1),size(Z,2),size(Z,3),1]),&
                        &col=1,description=description)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &Y=reshape(Y,[size(Y,1),size(Y,2),size(Y,3),1]),&
                        &col=1,description=description)
                end if
            else
                if (present(Z)) then
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &Z=reshape(Z,[size(Z,1),size(Z,2),size(Z,3),1]),&
                        &col=1,description=description)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &col=1,description=description)
                end if
            end if
        else
            if (present(Y)) then
                if (present(Z)) then
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &Y=reshape(Y,[size(Y,1),size(Y,2),size(Y,3),1]),&
                        &Z=reshape(Z,[size(Z,1),size(Z,2),size(Z,3),1]),&
                        &col=1,description=description)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &Y=reshape(Y,[size(Y,1),size(Y,2),size(Y,3),1]),&
                        &col=1,description=description)
                end if
            else
                if (present(Z)) then
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &col=1,Z=reshape(Z,[size(Z,1),size(Z,2),size(Z,3),1]),&
                        &description=description)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &col=1,description=description)
                end if
            end if
        end if
    end subroutine plot_HDF5_ind
    
    ! Takes  two input  vectors and  plots  these as  well as  the relative  and
    ! absolute difference in  a HDF5 file, similar to  plot_HDF5. Optionally, an
    ! output message  can be displayed on  screen with the maximum  relative and
    ! absolute error.
    subroutine plot_diff_HDF5(A,B,file_name,tot_dim,loc_offset,description,&
        &output_message)
        use messages, only: lvl_ud
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:), B(:,:,:)                              ! vectors A and B
        character(len=*), intent(in) :: file_name                               ! name of plot
        integer, intent(in), optional :: tot_dim(3)                             ! total dimensions of the arrays
        integer, intent(in), optional :: loc_offset(3)                          ! offset of local dimensions
        character(len=*), intent(in), optional :: description                   ! description
        logical, intent(in), optional :: output_message                         ! whether to display a message or not
        
        ! local variables
        real(dp), allocatable :: plot_var(:,:,:,:)                              ! variable containing plot
        integer :: tot_dim_loc(4)                                               ! local version of tot_dim
        integer :: loc_offset_loc(4)                                            ! local version of loc_offset
        character(len=max_str_ln) :: var_names(5)                               ! names of variables in plot
        logical :: output_message_loc                                           ! local version of output_message
        real(dp) :: lim_lo                                                      ! lower limit of errors
        real(dp) :: lim_hi                                                      ! upper limit of errors
        real(dp) :: err_av                                                      ! average error
        real(dp), allocatable :: tot_lim(:)                                     ! total lower or upper limit
        real(dp), allocatable :: tot_err(:)                                     ! total sum of errors
        integer :: istat                                                        ! status
        logical :: ind_plot                                                     ! individual plot or not
        
        ! set up local tot_dim and loc_offset
        tot_dim_loc = [shape(A),5]
        if (present(tot_dim)) tot_dim_loc = [tot_dim,5]
        loc_offset_loc = 0
        if (present(loc_offset)) loc_offset_loc = [loc_offset,5]
        
        ! tests
        if (size(A,1).ne.size(B,1) .or. size(A,2).ne.size(B,2) .or. &
            &size(A,3).ne.size(B,3)) then
            call writo('WARNING: in plot_diff_HDF5, A and B need to have the &
                &correct size',persistent=.true.)
            return
        end if
        
        ! set up local ouput message
        output_message_loc = .false.
        if (present(output_message)) output_message_loc = output_message
        
        ind_plot = .false.
        if (tot_dim_loc(1).eq.size(A,1) .and. tot_dim_loc(2).eq.size(A,2) &
            &.and. tot_dim_loc(3).eq.size(A,3)) ind_plot = .true.
        
        ! set up plot_var
        allocate(plot_var(size(A,1),size(A,2),size(A,3),5))
        plot_var(:,:,:,1) = A
        plot_var(:,:,:,2) = B
        plot_var(:,:,:,3) = diff(A,B,shape(A),rel=.true.)                       ! rel. diff.
        plot_var(:,:,:,4) = log10(abs(diff(A,B,shape(A),rel=.true.)))           ! log of abs. value of rel. diff.
        plot_var(:,:,:,5) = diff(A,B,shape(A),rel=.false.)                      ! abs. diff.
        
        ! set up var_names
        var_names(1) = 'v1'
        var_names(2) = 'v2'
        var_names(3) = 'rel v1 - v2'
        var_names(4) = 'log abs rel v1 - v2'
        var_names(5) = 'abs v1 - v2'
        
        ! plot
        call plot_HDF5_arr(var_names,file_name,plot_var,tot_dim=tot_dim_loc,&
            &loc_offset=loc_offset_loc,col=1,description=description)
        
        ! output message if requested
        if (output_message_loc) then
            call writo('Information about errors:',persistent=.true.)
            call lvl_ud(1)
            ! relative error
            call stats(tot_dim_loc(1:3),plot_var(:,:,:,3),lim_lo,lim_hi,err_av)
            call writo(trim(r2strt(lim_lo))//' < rel. err. < '//&
                &trim(r2strt(lim_hi))//', average value: '//&
                &trim(r2strt(err_av)),persistent=.true.)
            ! log of absolute relative error
            call stats(tot_dim_loc(1:3),plot_var(:,:,:,4),lim_lo,lim_hi,err_av)
            call writo(trim(r2strt(lim_lo))//' < log(abs(rel. err.)) < '//&
                &trim(r2strt(lim_hi))//', average value: '//&
                &trim(r2strt(err_av)),persistent=.true.)
            ! absolute error
            call stats(tot_dim_loc(1:3),plot_var(:,:,:,5),lim_lo,lim_hi,err_av)
            call writo(trim(r2strt(lim_lo))//' < abs. err. < '//&
                &trim(r2strt(lim_hi))//', average value: '//&
                &trim(r2strt(err_av)),persistent=.true.)
            call lvl_ud(-1)
        end if
    contains
        ! returns relative or absolute difference between inputs A and B
        function diff(A,B,dims,rel) result(C)
            ! local variables
            real(dp) :: max_diff = 1.E10                                        ! maximum absolute difference
            
            ! input / output
            real(dp), intent(in) :: A(:,:,:)                                    ! input A
            real(dp), intent(in) :: B(:,:,:)                                    ! input B
            integer, intent(in) :: dims(3)                                      ! dimensions of A and B
            logical, intent(in) :: rel                                          ! .true. if relative and .false. if absolute error
            real(dp) :: C(dims(1),dims(2),dims(3))                              ! output C
            
            ! return output
            if (rel) then
                C = min(max_diff,max(-max_diff,(A-B)/(A+B)))
            else
                C = abs(A-B)
            end if
        end function diff
        
        ! returns limits and average value
        subroutine stats(tot_dims,var,lim_lo,lim_hi,err_av)
            use MPI_utilities, only: get_ser_var
            
            ! input / output
            integer, intent(in) :: tot_dims(3)                                  ! total dimensions of grid
            real(dp), intent(in) :: var(:,:,:)                                  ! input variable for which to calculate statistics
            real(dp), intent(inout) :: lim_lo                                   ! lower limit of errors
            real(dp), intent(inout) :: lim_hi                                   ! upper limit of errors
            real(dp), intent(inout) :: err_av                                   ! average error
            
            ! local variables
            real(dp) :: sum_err                                                 ! sum of errors
            
            ! calculate limits on absolute error
            lim_lo = minval(var)
            lim_hi = maxval(var)
            sum_err = sum(var)
            
            ! get most stringent limits from all processes
            if (.not.ind_plot) then
                istat = get_ser_var([lim_lo],tot_lim)
                CHCKSTT
                lim_lo = minval(tot_lim)
                istat = get_ser_var([lim_hi],tot_lim)
                CHCKSTT
                lim_hi = maxval(tot_lim)
                istat = get_ser_var([sum_err],tot_err)
                CHCKSTT
                sum_err = sum(tot_err)
            end if
            
            ! calculate average error
            err_av = sum_err/product(tot_dims)
        end subroutine
    end subroutine plot_diff_HDF5
    
    ! Executes command line, or displays a message if disabled.
    subroutine use_execute_command_line(command,exitstat,cmdstat,cmdmsg)
        use num_vars, only: no_execute_command_line
        
        ! input / output
        character(len=*), intent(in) :: command                                 ! command to execute
        integer, intent(inout), optional :: exitstat                            ! exit status
        integer, intent(inout), optional :: cmdstat                             ! command status
        character(len=*), intent(inout), optional :: cmdmsg                     ! command message
        
        if (no_execute_command_line) then
            call writo('WARNING: Not executing command')
            call lvl_ud(1)
            call writo(command)
            call lvl_ud(-1)
            call writo('This can be run manually')
            if (present(exitstat)) exitstat = 1
            if (present(cmdstat)) cmdstat = 0
            if (present(cmdmsg)) cmdmsg = ''
        else
            call execute_command_line(command,EXITSTAT=exitstat,CMDSTAT=cmdstat,&
                &CMDMSG=cmdmsg)
        end if
    end subroutine
end module output_ops
