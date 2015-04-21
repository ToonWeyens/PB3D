!------------------------------------------------------------------------------!
!   Operations concerning  giving output, on the  screen as well as  in output !
!   files                                                                      !
!------------------------------------------------------------------------------!
module output_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln, no_plots, iu, plot_dir, data_dir, &
        &script_dir
    use files, only: nextunit
    
    implicit none
    private
    public print_GP_2D, print_GP_3D, draw_GP, draw_GP_animated, merge_GP, &
        &print_HDF5, plot_diff_HDF5
    
    ! interfaces
    interface print_GP_2D
        module procedure print_GP_2D_ind, print_GP_2D_arr
    end interface
    interface print_GP_3D
        module procedure print_GP_3D_ind, print_GP_3D_arr
    end interface
    interface print_HDF5
        module procedure print_HDF5_ind, print_HDF5_arr
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
        use num_vars, only: grp_rank
        
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
        open(unit=nextunit(file_i),file=data_dir//'/'//trim(file_name),&
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
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(var_name,trim(file_name),nplt,.true.,.true.)
        
        if (trim(file_name_i).eq.'') then
            call execute_command_line('rm '//data_dir//'/'//trim(file_name))
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
                    &reshape(z,[npoints,1]),y=reshape(y,[npoints,1]),&
                    &x=reshape(x,[npoints,1]),draw=draw)
            else
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[npoints,1]),y=reshape(y,[npoints,1]),&
                    &draw=draw)
            end if
        else
            if (present(x)) then
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[npoints,1]),x=reshape(x,[npoints,1]),&
                    &draw=draw)
            else
                call print_GP_3D_arr(var_name,file_name,&
                    &reshape(z,[npoints,1]),draw=draw)
            end if
        end if
    end subroutine print_GP_3D_ind
    subroutine print_GP_3D_arr(var_name,file_name_i,z,y,x,draw)                 ! multiple plots
        use num_vars, only: grp_rank
        
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
        open(unit=nextunit(file_i),file=data_dir//'/'//trim(file_name),&
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
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(var_name,trim(file_name),nplt,.false.,.true.)
        
        if (trim(file_name_i).eq.'') then
            call execute_command_line('rm '//data_dir//'/'//trim(file_name))
        end if
    end subroutine print_GP_3D_arr
    
    ! use GNUPlot to draw a plot
    ! The variable file_name(s) holds the file(s)  which are to be plot. Each of
    ! them should  contain nplt plots,  arranged in columns.  Furthermore, is_2D
    ! determines whether the  plot is to be 2D or  3D and plot_on_screen whether
    ! the plot  is be  shown on  the screen,  or to  be saved  in a  file called
    ! [var_name].pdf. Finally, for  each of the file names,  an optional command
    ! draw_ops can be provided, that specifies the line style for the plots from
    ! the file.
    subroutine draw_GP_ind(var_name,file_name,nplt,is_2D,plot_on_screen,&
        &draw_ops)
        
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: var_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        character(len=*), intent(in), optional :: draw_ops                      ! extra drawing option
        
        ! call the array version
        if (present(draw_ops)) then
            call draw_GP_arr(var_name,[file_name],nplt,is_2D,&
                &plot_on_screen,draw_ops=[draw_ops])
        else
            call draw_GP_arr(var_name,[file_name],nplt,is_2D,&
                &plot_on_screen)
        end if
    end subroutine draw_GP_ind
    subroutine draw_GP_arr(var_name,file_names,nplt,is_2D,plot_on_screen,&
        &draw_ops)
        use num_vars, only: grp_rank
        
        ! input / output
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: var_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        character(len=*), intent(in), optional :: draw_ops(:)                   ! extra drawing options (one per file)
        
        ! local variables
        character(len=3*size(file_names)*max_str_ln) :: plot_cmd                ! individual plot command
        character(len=max_str_ln) :: script_name                                ! name of script, including path
        character(len=max_str_ln) :: loc_draw_op                                ! local draw option
        integer :: istat                                                        ! status of opening a file
        integer :: iplt, ifl                                                    ! counters
        integer :: nfl                                                          ! number of files to plot
        integer :: cmd_i                                                        ! file number for script file
        integer :: plt_count                                                    ! counts the number of plots
        
        ! return if no_plots
        if (no_plots) then
            call writo('WARNING: plot ignored because no_plots is on')
            return
        end if
        
        ! set up nfl
        nfl = size(file_names)
        
        ! tests
        if (present(draw_ops)) then
            if (nfl.ne.size(draw_ops)) then
                call writo('WARNING: in draw_GP, one value for draw_ops has to &
                    &be provided for each file to plot')
                return
            end if
        end if
        
        ! create the GNUPlot command
        if (plot_on_screen) then
            script_name = trim(script_dir)//'/'//'temp_script_draw_GP_'//&
                &trim(i2str(grp_rank))//'.gnu'
        else
            script_name = trim(script_dir)//'/'//trim(var_name)//'.gnu'
        end if
        
        ! open script file
        open(unit=nextunit(cmd_i),file=trim(script_name),iostat=istat)
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command')
            return
        end if
        
        ! initialize the GNUPlot script
        if (plot_on_screen) then                                                ! wxt terminal
            write(cmd_i,*) 'set grid; set border 4095 front linetype -1 &
                &linewidth 1.000; set terminal wxt;'
        else                                                                    ! pdf terminal
            write(cmd_i,*) 'set grid; set border 4095 front linetype -1 &
                &linewidth 1.000; set terminal pdf; set output '''//&
                &trim(plot_dir)//'/'//trim(var_name)//'.pdf'';'
        end if
        
        ! no legend if too many plots
        if (nplt.gt.10) write(cmd_i,*) 'set nokey;'
        
        ! set up line styles
        if (present(draw_ops)) then
            do ifl = 1,nfl
                write(cmd_i,*) trim(draw_ops(ifl))//';'
            end do
        else
            do plt_count = 1,min(nplt*nfl,size(line_clrs))
                write(cmd_i,*) 'set style line '//trim(i2str(plt_count))//&
                    &' lc rgb '//line_clrs(plt_count)//' '//trim(line_style)
            end do
        end if
        
        ! set up plt_count
        plt_count = 1
        
        ! individual plots
        loc_draw_op = ''
        if (is_2D) then
            write(cmd_i,*) 'plot \'
            do iplt = 1,nplt
                plot_cmd = ''
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = trim(i2str(ifl))
                    else
                        loc_draw_op = 'with linespoints linestyle '//&
                            &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//' title '''//trim(var_name)//&
                    &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                    &')'' '//trim(loc_draw_op)//','
                    if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
                    plt_count = plt_count + 1
                end do
                write(cmd_i,*) trim(plot_cmd)
            end do
        else
            write(cmd_i,*) 'splot \'
            do iplt = 1,nplt
                plot_cmd = ''
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = trim(i2str(ifl))
                    else
                        loc_draw_op = 'with linespoints linestyle '//&
                            &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//':'//trim(i2str(2*nplt+iplt))&
                    &//' title '''//trim(var_name)//' ('//trim(i2str(iplt))&
                    &//'/'//trim(i2str(nplt))//')'' '//trim(loc_draw_op)//','
                    if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
                    plt_count = plt_count + 1
                end do
                write(cmd_i,*) trim(plot_cmd)
            end do
        end if
        
        ! finishing the GNUPlot command
        write(cmd_i,*) ''
        if (plot_on_screen) write(cmd_i,*) 'pause -1'
        
        ! closing the GNUPlot command
        close(cmd_i)
        
        ! stop timer
        call stop_time
        
        ! call GNUPlot
        call execute_command_line('gnuplot "'//trim(script_name)//&
            &'" 2> /dev/null',EXITSTAT=istat)
        
        if (plot_on_screen) then
            if (istat.ne.0) then
                call writo('Failed to plot '//trim(var_name))
            end if
            call execute_command_line('rm '//trim(script_name))
        else
            if (istat.eq.0) then
                call writo('Created plot in output file '''//&
                    &trim(plot_dir)//'/'//trim(var_name)//'.pdf''')
            else
                call writo('Failed to create plot in output file '''//&
                    &trim(plot_dir)//'/'//trim(var_name)//'.pdf''')
            end if
        end if
        
        ! start timer
        call start_time
    end subroutine draw_GP_arr
    
    ! use GNUPlot to draw an animated (gif) plot
    ! The variable file_name(s) holds the file(s)  which are to be plot. Each of
    ! them should contain nplt plots, arranged  in columns. The plot is saved in
    ! a file called [var_name].pdf. Subsequently, for each of the file names, an
    ! optional command draw_ops  can be provided, that specifies  the line style
    ! for  the plots  from the  file. Finally,  also optional  commands for  the
    ! ranges of the figure and the delay between the frames can be specified
    subroutine draw_GP_animated_ind(var_name,file_name,nplt,is_2D,ranges,delay,&
        &draw_ops)
        
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: var_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        integer, intent(in), optional :: delay                                  ! time delay between plots
        character(len=*), intent(in), optional :: draw_ops                      ! extra commands
        
        ! call the array version
        if (present(draw_ops)) then
            call draw_GP_animated_arr(var_name,[file_name],nplt,is_2D,&
                &ranges=ranges,delay=delay,draw_ops=[draw_ops])
        else
            call draw_GP_animated_arr(var_name,[file_name],nplt,is_2D,&
                &delay=delay,ranges=ranges)
        end if
    end subroutine draw_GP_animated_ind
    subroutine draw_GP_animated_arr(var_name,file_names,nplt,is_2D,ranges,&
        &delay,draw_ops)
        ! input / output
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: var_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        integer, intent(in), optional :: delay                                  ! time delay between plots
        character(len=*), intent(in), optional :: draw_ops(:)                   ! extra commands
        
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
        
        ! return if no_plots
        if (no_plots) then
            call writo('WARNING: plot ignored because no_plots is on')
            return
        end if
        
        ! set up nfl
        nfl = size(file_names)
        
        ! tests
        if (present(draw_ops)) then
            if (nfl.ne.size(draw_ops)) then
                call writo('WARNING: in draw_GP_animated, one value for &
                    &draw_ops has to be provided for each file to plot')
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
        script_name = ''//trim(script_dir)//'/'//trim(var_name)//'.gnu'
        
        ! open script file
        open(unit=nextunit(cmd_i),file=trim(script_name),iostat=istat)
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command')
            return
        end if
        
        ! initialize the GNUPlot script
        write(cmd_i,*) 'set grid; set border 4095 front &
            &linetype -1 linewidth 1.000; set terminal gif animate delay '//&
            &trim(i2str(delay_loc))//' size 1280, 720; &
            &set output '''//trim(plot_dir)//'/'//trim(var_name)//&
            &'.gif'';'
        
        ! find and set ranges
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
                ! initialize ranges
                ranges_loc(:,1) = huge(1._dp)                                   ! minimum value
                ranges_loc(:,2) = -huge(1._dp)                                  ! maximum value
                
                do ifl = 1,nfl
                    call get_ranges(file_names(ifl),ranges_loc)
                end do
            end if
            
            ! set ranges
            write(cmd_i,*) 'set xrange ['//trim(r2str(ranges_loc(1,1)))//':'//&
                &trim(r2str(ranges_loc(1,2)))//'];'
            write(cmd_i,*) 'set yrange ['//trim(r2str(ranges_loc(2,1)))//':'//&
                &trim(r2str(ranges_loc(2,2)))//'];'
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
                ranges_loc(:,1) = 1.E14_dp                                      ! minimum value
                ranges_loc(:,2) = -1.E14_dp                                     ! maximum value
                
                do ifl = 1,nfl
                    call get_ranges(file_names(ifl),ranges_loc)
                end do
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
            
            ! other definitions for 3D plots
            write(cmd_i,*) 'set view 45,45;'
            write(cmd_i,*) 'set hidden3d offset 0'
        end if
        
        ! no legend if too many plots
        !if (nfl.gt.10) write(cmd_i,*) 'set nokey;'
        
        ! set up line styles
        if (present(draw_ops)) then
            do ifl = 1,nfl
                write(cmd_i,*) 'set style line '//trim(i2str(ifl))//' '//&
                    &trim(draw_ops(ifl))//';'
            end do
        else
            do plt_count = 1,min(nplt*nfl,size(line_clrs))
                write(cmd_i,*) 'set style line '//trim(i2str(plt_count))//&
                    &' lc rgb '//line_clrs(plt_count)//' '//trim(line_style)
            end do
        end if
        
        ! individual plots
        loc_draw_op = ''
        if (is_2D) then
            do iplt = 1,nplt
                plot_cmd = 'plot'
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = trim(i2str(ifl))
                    else
                        loc_draw_op = 'with linespoints linestyle '//&
                            &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//' title '''//trim(var_name)//&
                    &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                    &')'' '//trim(loc_draw_op)
                    if (ifl.ne.nfl) plot_cmd = trim(plot_cmd)//', '
                end do
                write(cmd_i,*) trim(plot_cmd)
            end do
        else
            do iplt = 1,nplt
                plot_cmd = 'splot '
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = trim(i2str(ifl))
                    else
                        loc_draw_op = 'with linespoints linestyle '//&
                            &trim(i2str(mod(plt_count-1,size(line_clrs))+1))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//':'//trim(i2str(2*nplt+iplt))&
                    &//' title '''//trim(var_name)//' ('//trim(i2str(iplt))&
                    &//'/'//trim(i2str(nplt))//')'' '//trim(loc_draw_op)
                    if (ifl.ne.nfl) plot_cmd = trim(plot_cmd)//', '
                end do
                write(cmd_i,*) trim(plot_cmd)
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
                &trim(plot_dir)//'/'//trim(var_name)//'.gif''')
        else
            call writo('Failed to create animated plot in output file '''//&
                &trim(plot_dir)//'/'//trim(var_name)//'.gif''')
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
            if (is_2D) then
                allocate(loc_data(2*nplt))
            else
                allocate(loc_data(3*nplt))
            end if
            
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
        call execute_command_line(trim(shell_cmd),EXITSTAT=istat)
        if (istat.ne.0) then
            call writo('WARNING: merge_GP failed to merge')
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
            call execute_command_line(trim(shell_cmd),EXITSTAT=istat)
            if (istat.ne.0) then
                call writo('WARNING: merge_GP failed to delete')
                return
            end if
        end if
    end subroutine

    ! Prints variables  "vars" with names "var_names"  in a HDF5 file  with name
    ! "file_name" and  accompanying XDMF file.  For collections, only  the first
    ! value for var_names is used, so it should have a size of one.
    ! The plot is generally  3D, but if one of the  dimensions provided is equal
    ! to 1, the plot becomes 2D. THIS SHOULD BE DONE ONLY IF IT IS INDEED SO!!!
    ! The axes of the 2D plot correpond to X-Z if the second dimension (zeta) is
    ! singular,  or X-Y  if  it is  the  third (theta).  This  corresponds to  a
    ! poloidal or a toroidal cross section, respectively.
    ! Optionally, the  (curvilinear) grid can  be provided through "X",  "Y" and
    ! "Z". If so, "Y" or "Z" may be neglected if it is a 2D plot.
    ! Additionally, the total grid size has  to be provided in "tot_dim", and if
    ! the routine  is called in  parallel, the  group dimensions and  offsets as
    ! well, in "grp_dim" and "grp_offset". 
    ! Optionally, one of the dimensions (col_id, default 4) can be associated to
    ! a collection dimension using "col" different from 1:
    !   col = 1: no collection, just plots of different variables
    !   col = 2: time collection
    !   col = 3: spatial collection
    ! Note: To plot this with VisIt, use:
    !   - for temporal collections: pseudocolor using the variable name (other 
    !       names are ignored)
    !   - for spatial collections:  subset of blocks or pseudocolor using the 
    !       variable name (other names are ignored)
    !   - for without collections:  pseudocolor using the variables name_i
    ! Note: Currently all of possibly multiple processes that write data have to
    ! cover the  entire range  of the  variables in  the dimension  indicated by
    ! col_id.  (This could  be implemented  by  changing n_plot  is defined  and
    ! selectively  letting  each  processer  write  in  the  main  loop  at  its
    ! corresponding indices.)
    subroutine print_HDF5_arr(var_names,file_name,vars,tot_dim,grp_dim,&
        &grp_offset,X,Y,Z,col_id,col,description)                               ! array version
        use HDF5_ops, only: open_HDF5_file, add_HDF5_item, print_HDF5_top, &
            &print_HDF5_geom, print_HDF5_3D_data_item, print_HDF5_att, &
            &print_HDF5_grid, close_HDF5_file, &
            &XML_str_type, HDF5_file_type
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            ! names of variable to be plot
        character(len=*), intent(in) :: file_name                               ! file name
        real(dp), intent(in), target :: vars(:,:,:,:)                           ! variables to plot
        integer, intent(in), optional :: tot_dim(4)                             ! total dimensions of the arrays
        integer, intent(in), optional :: grp_dim(4)                             ! group dimensions of the arrays
        integer, intent(in), optional :: grp_offset(4)                          ! offset of group dimensions
        real(dp), intent(in), target, optional :: X(:,:,:,:)                    ! curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:,:)                    ! curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:,:)                    ! curvlinear grid Z points
        integer, intent(in), optional :: col_id                                 ! index of time dimension
        integer, intent(in), optional :: col                                    ! whether a collection is made
        character(len=*), intent(in), optional :: description                   ! description
        
        ! local variables
        integer :: istat                                                        ! status
        integer :: col_id_loc                                                   ! local copy of col_id
        integer :: col_loc                                                      ! local copy of col
        type(HDF5_file_type) :: file_info                                       ! file info
        integer :: n_plot                                                       ! nr. of plots
        integer :: id, jd                                                       ! counter
        integer :: tot_dim_loc(4)                                               ! local copy of tot_dim
        integer :: tot_dim_3D(3)                                                ! tot_dim except collection
        integer :: grp_dim_3D(3)                                                ! grp_dim except collection
        integer :: grp_offset_3D(3)                                             ! grp_offset except collection
        type(XML_str_type) :: col_grid                                          ! grid with collection
        type(XML_str_type), allocatable :: grids(:)                             ! the grids in the time collection
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type), allocatable :: XYZ(:)                               ! data items for geometry
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        logical :: col_mask(4)                                                  ! to select out the collection dimension
        real(dp), pointer :: var_3D(:,:,:)                                      ! pointer to vars
        real(dp), pointer :: X_3D(:,:,:)                                        ! pointer to X
        real(dp), pointer :: Y_3D(:,:,:)                                        ! pointer to Y 
        real(dp), pointer :: Z_3D(:,:,:)                                        ! pointer to Z
        character(len=max_str_ln), allocatable :: grd_names(:)                  ! grid names
        character(len=max_str_ln), allocatable :: att_names(:)                  ! attribute names
        
        ! set up local col_id and col
        col_id_loc = 4                                                          ! default collection dimension: last index
        if (present(col_id)) col_id_loc = col_id
        col_loc = 1                                                             ! default no spatial collection
        if (present(col)) col_loc = col
        
        ! set up local tot_dim
        tot_dim_loc = shape(vars)
        if (present(tot_dim)) tot_dim_loc = tot_dim
        
        ! tests
        if (present(grp_dim)) then
            if (tot_dim_loc(col_id_loc).ne.grp_dim(col_id_loc)) then
                istat = 1
                call writo('WARNING: In print_HDF5, all the processes need to &
                    &have the full range in the dimension given by col_id')
                CHCKSTT
            end if
        end if
        
        ! set real dimensions
        col_mask = .false.
        col_mask(col_id_loc) = .true.
        tot_dim_3D = pack(tot_dim_loc,.not.col_mask)
        if (present(grp_dim)) then
            grp_dim_3D = pack(grp_dim,.not.col_mask)
        else
            grp_dim_3D = tot_dim_3D
        end if
        if (present(grp_offset)) then
            grp_offset_3D = pack(grp_offset,.not.col_mask)
        else
            grp_offset_3D = 0
        end if
        
        ! set up nr. of plots
        n_plot = size(vars,col_id_loc)
        
        ! if nr of plots is equal to one, individual version should be used
        ! (otherwise an error occurs in VisIt)
        if (n_plot.eq.1) then
            ! test
            if (size(var_names).ne.1) then
                call writo('WARNING: var_names in print_HDF5 has to have &
                    &length 1 for individual plot')
            end if
            ! assign pointers
            call assign_pointers(1)
            ! call individual version with first name
            call print_HDF5_ind(var_names(1),file_name,var_3D,tot_dim_3D,&
                &grp_dim_3D,grp_offset_3D,X_3D,Y_3D,Z_3D,description)
            ! return
            return
        end if
        
        ! set up local var_names
        allocate(grd_names(n_plot))
        allocate(att_names(n_plot))
        if (col_loc.eq.1) then                                                  ! without collection: grid name is important
            if (size(var_names).eq.n_plot) then                                 ! the right number of variable names provided
                grd_names = var_names
            else if (size(var_names).gt.n_plot) then                            ! too many variable names provided
                grd_names = var_names(1:n_plot)
                call writo('WARNING: Too many variable names provided')
            else                                                                ! not enough variable names provided
                grd_names(1:size(var_names)) = var_names
                do id = size(var_names)+1,n_plot
                    grd_names(id) = 'variable '//trim(i2str(id))
                end do
                call writo('WARNING: Not enough variable names provided')
            end if
            att_names = 'var'
        else                                                                    ! collections: attribute name is important
            att_names = var_names(1)
            if (size(var_names).gt.1) then
                call writo('WARNING: For collections, only the first variable &
                    &name is used')
            end if
            grd_names = 'grid'
        end if
        
        ! open HDF5 file
        if (tot_dim_3D(1).eq.grp_dim_3D(1) .and. &
            &tot_dim_3D(2).eq.grp_dim_3D(2) .and. &
            &tot_dim_3D(3).eq.grp_dim_3D(3)) then
            istat = open_HDF5_file(file_info,file_name,description,&
                &ind_plot=.true.)                                               ! individual plot
        else
            istat = open_HDF5_file(file_info,file_name,description)             ! group plot
        end if
        CHCKSTT
        
        ! create grid for collection
        allocate(grids(n_plot))
        
        ! allocate geometry arrays
        ! (first two indices correspond to angular dimensions)
        if (tot_dim_3D(1).eq.1 .or. tot_dim_3D(2).eq.1) then                    ! 2D grid
            allocate(XYZ(2))
        else                                                                    ! 3D grid
            allocate(XYZ(3))
        end if
        
        ! loop over all plots
        do id = 1,n_plot
            ! print topology
            if (tot_dim_3D(1).eq.1 .or. tot_dim_3D(2).eq.1) then                ! 2D grid
                call print_HDF5_top(top,1,tot_dim_3D)
            else                                                                ! 3D grid
                call print_HDF5_top(top,2,tot_dim_3D)
            end if
            
            ! assign pointers
            call assign_pointers(id)
            
            ! print data item for X
            istat = print_HDF5_3D_data_item(XYZ(1),file_info,&
                &'X_'//trim(i2str(id)),X_3D,tot_dim_3D,grp_dim_3D,&
                &grp_offset_3D)
            CHCKSTT
            
            ! print data item for Y and / or Z
            if (tot_dim_3D(2).eq.1) then                                        ! if toroidally symmetric, no Y axis
                istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Z_'//trim(i2str(id)),Z_3D,tot_dim_3D,grp_dim_3D,&
                    &grp_offset_3D)
                CHCKSTT
            else if (tot_dim_3D(1).eq.1) then                                   ! if poloidally symmetric, no Z axis
                istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Y_'//trim(i2str(id)),Y_3D,tot_dim_3D,grp_dim_3D,&
                    &grp_offset_3D)
                CHCKSTT
            else
                istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                    &'Y_'//trim(i2str(id)),Y_3D,tot_dim_3D,grp_dim_3D,&
                    &grp_offset_3D)
                CHCKSTT
                istat = print_HDF5_3D_data_item(XYZ(3),file_info,&
                    &'Z_'//trim(i2str(id)),Z_3D,tot_dim_3D,grp_dim_3D,&
                    &grp_offset_3D)
                CHCKSTT
            end if
            
            ! print geometry with X, Y and Z data item
            if (tot_dim_3D(2).eq.1 .or. tot_dim_3D(1).eq.1) then                ! if symmetry 2D geometry
                call print_HDF5_geom(geom,1,XYZ,.true.)
            else                                                                ! if no symmetry 3D geometry
                call print_HDF5_geom(geom,2,XYZ,.true.)
            end if
            
            ! print data item for plot variable
            istat = print_HDF5_3D_data_item(XYZ(1),file_info,'var_'//&
                &trim(i2str(id)),var_3D,tot_dim_3D,grp_dim_3D,grp_offset_3D)    ! reuse XYZ(1)
            CHCKSTT
            
            ! print attribute with this data item
            call print_HDF5_att(att(1),XYZ(1),att_names(id),1,.true.)
            
            ! create a  grid with the topology, the geometry,  the attribute and
            ! time if time collection
            if (col_loc.eq.2) then                                              ! time collection
                istat = print_HDF5_grid(grids(id),grd_names(id),1,&
                    &grid_time=id*1._dp,grid_top=top,grid_geom=geom,&
                    &grid_atts=att,reset=.true.)
                CHCKSTT
            else                                                                ! no time collection
                istat = print_HDF5_grid(grids(id),grd_names(id),1,&
                    &grid_top=top,grid_geom=geom,grid_atts=att,reset=.true.)
                CHCKSTT
            end if
        end do
        
        ! either create collection or just use individual grids
        if (col_loc.eq.1) then
            ! add individual grids to HDF5 file and reset them
            do id = 1,n_plot
                call add_HDF5_item(file_info,grids(id),reset=.true.)
            end do
        else
            ! create grid collection from individual grids and reset them
            istat = print_HDF5_grid(col_grid,'domain of mesh',col_loc,&
                &grid_grids=grids,reset=.true.)
            CHCKSTT
            
            ! add collection grid to HDF5 file and reset it
            call add_HDF5_item(file_info,col_grid,reset=.true.)
        end if
        
        ! close HDF5 file
        istat = close_HDF5_file(file_info)
        CHCKSTT
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
                allocate(X_3D(grp_dim_3D(1),grp_dim_3D(2),grp_dim_3D(3)))
                do jd = 1,grp_dim_3D(1)
                    X_3D(jd,:,:) = grp_offset_3D(1) + jd
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
                allocate(Y_3D(grp_dim_3D(1),grp_dim_3D(2),grp_dim_3D(3)))
                do jd = 1,grp_dim_3D(2)
                    Y_3D(:,jd,:) = grp_offset_3D(2) + jd
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
                allocate(Z_3D(grp_dim_3D(1),grp_dim_3D(2),grp_dim_3D(3)))
                do jd = 1,grp_dim_3D(3)
                    Z_3D(:,:,jd) = grp_offset_3D(3) + jd
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
    end subroutine print_HDF5_arr
    subroutine print_HDF5_ind(var_name,file_name,var,tot_dim,grp_dim,&
        &grp_offset,X,Y,Z,description)                                          ! individual version
        use HDF5_ops, only: open_HDF5_file, add_HDF5_item, &
            &XML_str_type, HDF5_file_type, print_HDF5_top, print_HDF5_geom, &
            &print_HDF5_3D_data_item, print_HDF5_att, print_HDF5_grid, &
            &close_HDF5_file
        
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of variable to be plot
        character(len=*), intent(in) :: file_name                               ! file name
        real(dp), intent(in), target :: var(:,:,:)                              ! variable to plot
        integer, intent(in), optional :: tot_dim(3)                             ! total dimensions of the arrays
        integer, intent(in), optional :: grp_dim(3)                             ! group dimensions of the arrays
        integer, intent(in), optional :: grp_offset(3)                          ! offset of group dimensions
        real(dp), intent(in), target, optional :: X(:,:,:)                      ! curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:)                      ! curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:)                      ! curvlinear grid Z points
        character(len=*), intent(in), optional :: description                   ! description
        
        ! local variables
        integer :: istat                                                        ! status
        type(HDF5_file_type) :: file_info                                       ! file info
        integer :: id                                                           ! counter
        integer :: tot_dim_loc(3)                                               ! local version of tot_dim
        integer :: grp_dim_loc(3)                                               ! local version of grp_dim
        integer :: grp_offset_loc(3)                                            ! local version of grp_offset
        type(XML_str_type) :: grid                                              ! grid
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type), allocatable :: XYZ(:)                               ! data items for geometry
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        real(dp), pointer :: X_ptr(:,:,:)                                       ! pointer to X
        real(dp), pointer :: Y_ptr(:,:,:)                                       ! pointer to Y
        real(dp), pointer :: Z_ptr(:,:,:)                                       ! pointer to Z
        
        ! set up file info
        file_info%name = file_name
        
        ! set up local tot_dim
        tot_dim_loc = shape(var)
        if (present(tot_dim)) tot_dim_loc = tot_dim
        
        ! set up local grp_dim and grp_offset
        if (present(grp_dim)) then
            grp_dim_loc = grp_dim
        else
            grp_dim_loc = tot_dim_loc
        end if
        if (present(grp_offset)) then
            grp_offset_loc = grp_offset
        else
            grp_offset_loc = [0,0,0]
        end if
        
        ! open HDF5 file
        if (tot_dim_loc(1).eq.grp_dim_loc(1) .and. &
            &tot_dim_loc(2).eq.grp_dim_loc(2) .and. &
            &tot_dim_loc(3).eq.grp_dim_loc(3)) then
            istat = open_HDF5_file(file_info,file_name,description,&
                &ind_plot=.true.)                                               ! individual plot
        else
            istat = open_HDF5_file(file_info,file_name,description)             ! group plot
        end if
        CHCKSTT
        
        ! allocate geometry arrays
        if (tot_dim_loc(1).eq.1 .or. tot_dim_loc(2).eq.1) then                  ! 2D grid
            allocate(XYZ(2))
        else                                                                    ! 3D grid
            allocate(XYZ(3))
        end if
        
        ! print topology
        if (tot_dim_loc(1).eq.1 .or. tot_dim_loc(2).eq.1) then                  ! 2D grid
            call print_HDF5_top(top,1,tot_dim_loc)
        else                                                                    ! 3D grid
            call print_HDF5_top(top,2,tot_dim_loc)
        end if
            
        ! calculate geometry arrays
        if (present(X)) then
            X_ptr => X
        else
            allocate(X_ptr(size(var,1),size(var,2),size(var,3)))
            do id = 1,size(var,1)
                X_ptr(id,:,:) = grp_offset_loc(1) + id
            end do
        end if
        if (present(Y)) then
            Y_ptr => Y
        else
            allocate(Y_ptr(size(var,1),size(var,2),size(var,3)))
            do id = 1,size(var,2)
                Y_ptr(:,id,:) =  grp_offset_loc(2) +id
            end do
        end if
        if (present(Z)) then
            Z_ptr => Z
        else
            allocate(Z_ptr(size(var,1),size(var,2),size(var,3)))
            do id = 1,size(var,3)
                Z_ptr(:,:,id) =  grp_offset_loc(3) +id
            end do
        end if
        
        ! print data item for X
        istat = print_HDF5_3D_data_item(XYZ(1),file_info,'X',X_ptr,&
            &tot_dim_loc,grp_dim=grp_dim_loc,grp_offset=grp_offset_loc)
        CHCKSTT
        
        ! print data item for Y and / or Z
        if (tot_dim_loc(2).eq.1) then                                           ! if toroidally symmetric, no Y axis
            istat = print_HDF5_3D_data_item(XYZ(2),file_info,'Z',Z_ptr,&
                &tot_dim_loc,grp_dim=grp_dim_loc,grp_offset=grp_offset_loc)
            CHCKSTT
        else if (tot_dim_loc(1).eq.1) then                                      ! if poloidally symmetric, no Z axis
            istat = print_HDF5_3D_data_item(XYZ(2),file_info,'Y',Y_ptr,&
                &tot_dim_loc,grp_dim=grp_dim_loc,grp_offset=grp_offset_loc)
            CHCKSTT
        else
            istat = print_HDF5_3D_data_item(XYZ(2),file_info,'Y',Y_ptr,&
                &tot_dim_loc,grp_dim=grp_dim_loc,grp_offset=grp_offset_loc)
            CHCKSTT
            istat = print_HDF5_3D_data_item(XYZ(3),file_info,'Z',Z_ptr,&
                &tot_dim_loc,grp_dim=grp_dim_loc,grp_offset=grp_offset_loc)
            CHCKSTT
        end if
        
        ! print geometry with X, Y and Z data item
        if (tot_dim_loc(2).eq.1 .or. tot_dim_loc(1).eq.1) then                  ! if symmetry 2D geometry
            call print_HDF5_geom(geom,1,XYZ,.true.)
        else                                                                    ! if no symmetry 3D geometry
            call print_HDF5_geom(geom,2,XYZ,.true.)
        end if
        
        ! print data item for plot variable
        istat = print_HDF5_3D_data_item(XYZ(1),file_info,'var',&
            &var,tot_dim_loc,grp_dim=grp_dim_loc,grp_offset=grp_offset_loc)     ! reuse XYZ(1)
        CHCKSTT
        
        ! print attribute with this data item
        call print_HDF5_att(att(1),XYZ(1),var_name,1,.true.)
        
        ! create a grid with the topology, the geometry and the attribute
        istat = print_HDF5_grid(grid,'grid',1,grid_top=top,&
            &grid_geom=geom,grid_atts=att,reset=.true.)
        CHCKSTT
        
        ! add grid to HDF5 file and reset it
        call add_HDF5_item(file_info,grid,reset=.true.)
        
        ! close HDF5 file
        istat = close_HDF5_file(file_info)
        CHCKSTT
    end subroutine print_HDF5_ind
    
    ! Takes  two input  vectors and  plots  these as  well as  the relative  and
    ! absolute difference in a HDF5  file, similar to print_HDF5. Optionally, an
    ! output message  can be displayed on  screen with the maximum  relative and
    ! absolute error.
    subroutine plot_diff_HDF5(A,B,file_name,tot_dim,grp_dim,grp_offset,&
        &description,output_message)
        use messages, only: lvl_ud
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:), B(:,:,:)                              ! vectors A and B
        character(len=*), intent(in) :: file_name                               ! name of plot
        integer, intent(in), optional :: tot_dim(3)                             ! total dimensions of the arrays
        integer, intent(in), optional :: grp_dim(3)                             ! group dimensions of the arrays
        integer, intent(in), optional :: grp_offset(3)                          ! offset of group dimensions
        character(len=*), intent(in), optional :: description                   ! description
        logical, intent(in), optional :: output_message                         ! whether to display a message or not
        
        ! local variables
        real(dp), allocatable :: plot_var(:,:,:,:)                              ! variable containing plot
        integer :: tot_dim_loc(4)                                               ! local version of tot_dim
        integer :: grp_dim_loc(4)                                               ! local version of grp_dim
        integer :: grp_offset_loc(4)                                            ! local version of grp_offset
        character(len=max_str_ln) :: var_names(5)                               ! names of variables in plot
        logical :: output_message_loc                                           ! local version of output_message
        real(dp) :: lim_lo                                                      ! lower limit of errors
        real(dp) :: lim_hi                                                      ! upper limit of errors
        real(dp) :: err_av                                                      ! average error
        real(dp), allocatable :: tot_lim(:)                                     ! total lower or upper limit
        real(dp), allocatable :: tot_err(:)                                     ! total sum of errors
        integer :: istat                                                        ! status
        logical :: ind_plot                                                     ! individual plot or not
        
        ! set up local tot_dim
        tot_dim_loc = [shape(A),5]
        if (present(tot_dim)) tot_dim_loc = [tot_dim,5]
        
        ! set up local grp_dim and grp_offset
        if (present(grp_dim)) then
            grp_dim_loc = [grp_dim,5]
        else
            grp_dim_loc = tot_dim_loc
        end if
        if (present(grp_offset)) then
            grp_offset_loc = [grp_offset,5]
        else
            grp_offset_loc = [0,0,0,0]
        end if
        
        ! tests
        if (grp_dim_loc(1).ne.size(B,1) .or. grp_dim_loc(2).ne.size(B,2) .or. &
            &grp_dim_loc(3).ne.size(B,3) .or. grp_dim_loc(1).ne.size(A,1) .or. &
            &grp_dim_loc(2).ne.size(A,2) .or. grp_dim_loc(3).ne.size(A,3)) then
            call writo('WARNING: in plot_diff_HDF5, A and B need to have the &
                &correct size')
            return
        end if
        
        ! set up local ouput message
        output_message_loc = .false.
        if (present(output_message)) output_message_loc = output_message
        
        ind_plot = .false.
        if (tot_dim_loc(1).eq.grp_dim_loc(1) .and. &
            &tot_dim_loc(2).eq.grp_dim_loc(2) .and. &
            &tot_dim_loc(3).eq.grp_dim_loc(3)) ind_plot = .true.
        
        ! set up plot_var
        allocate(plot_var(grp_dim_loc(1),grp_dim_loc(2),grp_dim_loc(3),5))
        plot_var(:,:,:,1) = A
        plot_var(:,:,:,2) = B
        plot_var(:,:,:,3) = diff(A,B,grp_dim_loc,rel=.true.)                    ! rel. diff.
        plot_var(:,:,:,4) = log10(abs(diff(A,B,grp_dim_loc,rel=.true.)))        ! log of abs. value of rel. diff.
        plot_var(:,:,:,5) = diff(A,B,grp_dim_loc,rel=.false.)                   ! abs. diff.
        
        ! set up var_names
        var_names(1) = 'v1'
        var_names(2) = 'v2'
        var_names(3) = 'rel v1 - v2'
        var_names(4) = 'log abs rel v1 - v2'
        var_names(5) = 'abs v1 - v2'
        
        ! plot
        call print_HDF5_arr(var_names,file_name,plot_var,tot_dim=tot_dim_loc,&
            &grp_dim=grp_dim_loc,grp_offset=grp_offset_loc,col=1,&
            &description=description)
        
        ! output message if requested
        if (output_message_loc) then
            call writo('Information about errors:')
            call lvl_ud(1)
            ! relative error
            call stats(tot_dim_loc(1:3),plot_var(:,:,:,3),lim_lo,lim_hi,err_av)
            call writo(trim(r2strt(lim_lo))//' < rel. err. < '//&
                &trim(r2strt(lim_hi))//', average value: '//&
                &trim(r2strt(err_av)))
            ! log of absolute relative error
            call stats(tot_dim_loc(1:3),plot_var(:,:,:,4),lim_lo,lim_hi,err_av)
            call writo(trim(r2strt(lim_lo))//' < log(abs(rel. err.)) < '//&
                &trim(r2strt(lim_hi))//', average value: '//&
                &trim(r2strt(err_av)))
            ! absolute error
            call stats(tot_dim_loc(1:3),plot_var(:,:,:,5),lim_lo,lim_hi,err_av)
            call writo(trim(r2strt(lim_lo))//' < abs. err. < '//&
                &trim(r2strt(lim_hi))//', average value: '//&
                &trim(r2strt(err_av)))
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
end module output_ops
