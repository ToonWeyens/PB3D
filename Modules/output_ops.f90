!------------------------------------------------------------------------------!
!   This module contains operations concerning giving output, on the screen as !
!   well as in output files                                                    !
!------------------------------------------------------------------------------!
module output_ops
#include <PB3D_macros.h>
    use str_ops, only: i2str, r2strt, r2str
    use num_vars, only: dp, max_str_ln, no_plots, iu, plot_dir, data_dir, &
        &script_dir
    use message_ops, only: writo, stop_time, start_time, lvl_ud
    
    implicit none
    private
    public print_GP_2D, print_GP_3D, draw_GP, draw_GP_animated, merge_GP, &
        &print_HDF5_3D
    
    ! interfaces
    interface print_GP_2D
        module procedure print_GP_2D_ind, print_GP_2D_arr
    end interface
    interface print_GP_3D
        module procedure print_GP_3D_ind, print_GP_3D_arr
    end interface
    interface print_HDF5_3D
        module procedure print_HDF5_3D_ind, print_HDF5_3D_arr
    end interface
    interface draw_GP
        module procedure draw_GP_ind, draw_GP_arr
    end interface
    interface draw_GP_animated
        module procedure draw_GP_animated_ind, draw_GP_animated_arr
    end interface
    
contains
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
        if (present(x)) then
            call print_GP_2D_arr(fun_name,file_name,reshape(y,[npoints,1]),&
                &x=reshape(x,[npoints,1]),draw=draw)
        else
            call print_GP_2D_arr(fun_name,file_name,reshape(y,[npoints,1]),&
                &draw=draw)
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
        
        ! if draw is present and equal to .false., cancel calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(fun_name,trim(file_name),nplt,.false.,.true.)
        
        if (trim(file_name_i).eq.'') then
            call execute_command_line('rm '//data_dir//'/'//trim(file_name))
        end if
    end subroutine print_GP_3D_arr
    
    ! use GNUPlot to draw a plot
    ! The variable file_name(s) holds the file(s)  which are to be plot. Each of
    ! them should  contain nplt plots,  arranged in columns.  Furthermore, is_2D
    ! determines whether the  plot is to be 2D or  3D and plot_on_screen whether
    ! the plot  is be  shown on  the screen,  or to  be saved  in a  file called
    ! [fun_name].pdf. Finally, for  each of the file names,  an optional command
    ! draw_ops can be provided, that specifies the line style for the plots from
    ! the file.
    subroutine draw_GP_ind(fun_name,file_name,nplt,is_2D,plot_on_screen,&
        &draw_ops)
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        character(len=*), intent(in), optional :: draw_ops                      ! extra drawing option
        
        ! call the array version
        if (present(draw_ops)) then
            call draw_GP_arr(fun_name,[file_name],nplt,is_2D,&
                &plot_on_screen,draw_ops=[draw_ops])
        else
            call draw_GP_arr(fun_name,[file_name],nplt,is_2D,&
                &plot_on_screen)
        end if
    end subroutine draw_GP_ind
    subroutine draw_GP_arr(fun_name,file_names,nplt,is_2D,plot_on_screen,&
        &draw_ops)
        use safe_open_mod, only: safe_open
        use num_vars, only: grp_rank, n_seq_0
        
        ! input / output
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
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
            script_name = trim(script_dir)//'/'//trim(fun_name)//'.gnu'
        end if
        
        ! open script file
        cmd_i = n_seq_0
        call safe_open(cmd_i,istat,trim(script_name),'replace','formatted',&
            &delim_in='none')
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
                &trim(plot_dir)//'/'//trim(fun_name)//'.pdf'';'
        end if
        
        ! no legend if too many plots
        if (nplt.gt.10) write(cmd_i,*) 'set nokey;'
        
        ! set up line styles
        if (present(draw_ops)) then
            do ifl = 1,nfl
                write(cmd_i,*) 'set style line '//trim(i2str(ifl))//' '//&
                    &trim(draw_ops(ifl))//';'
            end do
        end if
        
        ! individual plots
        loc_draw_op = ''
        if (is_2D) then
            write(cmd_i,*) 'plot \'
            do iplt = 1,nplt
                plot_cmd = ''
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = 'linestyle '//trim(i2str(ifl))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                    &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                    &')'' with lines '//trim(loc_draw_op)//','
                    if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
                end do
                write(cmd_i,*) trim(plot_cmd)
            end do
        else
            write(cmd_i,*) 'splot \'
            do iplt = 1,nplt
                plot_cmd = ''
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = 'linestyle '//trim(i2str(ifl))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//':'//trim(i2str(2*nplt+iplt))&
                    &//' title '''//trim(fun_name)//' ('//trim(i2str(iplt))&
                    &//'/'//trim(i2str(nplt))//')'' with lines '//&
                    &trim(loc_draw_op)//','
                    if (ifl.eq.nfl) plot_cmd = trim(plot_cmd)//' \'
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
    end subroutine draw_GP_arr
    
    ! use GNUPlot to draw an animated (gif) plot
    ! The variable file_name(s) holds the file(s)  which are to be plot. Each of
    ! them should contain nplt plots, arranged  in columns. The plot is saved in
    ! a file called [fun_name].pdf. Subsequently, for each of the file names, an
    ! optional command draw_ops  can be provided, that specifies  the line style
    ! for  the plots  from the  file. Finally,  also optional  commands for  the
    ! ranges of the figure and the delay between the frames can be specified
    subroutine draw_GP_animated_ind(fun_name,file_name,nplt,is_2D,ranges,delay,&
        &draw_ops)
        
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
        integer, intent(in) :: nplt                                             ! number of plots
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        real(dp), intent(in), optional :: ranges(:,:)                           ! x and y range, and z range (if 3D) of plot
        integer, intent(in), optional :: delay                                  ! time delay between plots
        character(len=*), intent(in), optional :: draw_ops                      ! extra commands
        
        ! call the array version
        if (present(draw_ops)) then
            call draw_GP_animated_arr(fun_name,[file_name],nplt,is_2D,&
                &ranges=ranges,delay=delay,draw_ops=[draw_ops])
        else
            call draw_GP_animated_arr(fun_name,[file_name],nplt,is_2D,&
                &delay=delay,ranges=ranges)
        end if
    end subroutine draw_GP_animated_ind
    subroutine draw_GP_animated_arr(fun_name,file_names,nplt,is_2D,ranges,&
        &delay,draw_ops)
        use num_vars, only: n_seq_0
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: file_names(:)                           ! name of file
        character(len=*), intent(in) :: fun_name                                ! name of function
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
        script_name = ''//trim(script_dir)//'/'//trim(fun_name)//'.gnu'
        
        ! open script file
        cmd_i = n_seq_0
        call safe_open(cmd_i,istat,trim(script_name),'replace','formatted',&
            &delim_in='none')
        if (istat.ne.0) then
            call writo('WARNING: Could not open file for gnuplot command')
            return
        end if
        
        ! initialize the GNUPlot script
        write(cmd_i,*) 'set grid; set border 4095 front &
            &linetype -1 linewidth 1.000; set terminal gif animate delay '//&
            &trim(i2str(delay_loc))//' size 1280, 720; &
            &set output '''//trim(plot_dir)//'/'//trim(fun_name)//&
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
        end if
        
        ! individual plots
        loc_draw_op = ''
        if (is_2D) then
            do iplt = 1,nplt
                plot_cmd = 'plot'
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = 'linestyle '//trim(i2str(ifl))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                    &' ('//trim(i2str(iplt))//'/'//trim(i2str(nplt))//&
                    &')'' with lines '//trim(loc_draw_op)
                    if (ifl.ne.nfl) plot_cmd = trim(plot_cmd)//', '
                end do
                write(cmd_i,*) trim(plot_cmd)
            end do
        else
            do iplt = 1,nplt
                plot_cmd = 'splot '
                do ifl = 1,nfl
                    if (present(draw_ops)) then
                        loc_draw_op = 'linestyle '//trim(i2str(ifl))
                    end if
                    plot_cmd = trim(plot_cmd)//' '''//trim(data_dir)//'/'//&
                    &trim(file_names(ifl))//''' using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//':'//trim(i2str(2*nplt+iplt))&
                    &//' title '''//trim(fun_name)//' ('//trim(i2str(iplt))&
                    &//'/'//trim(i2str(nplt))//')'' with lines '//&
                    &trim(loc_draw_op)
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
                &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
        else
            call writo('Failed to create animated plot in output file '''//&
                &trim(plot_dir)//'/'//trim(fun_name)//'.gif''')
        end if
    contains
        subroutine get_ranges(file_name,ranges)
            use safe_open_mod, only: safe_open
            use num_vars, only: n_seq_0
            
            ! input / output
            character(len=*), intent(in) :: file_name                           ! name of file
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

    ! Prints  variables "vars" with names  "var_names" in a HDF5  file with name
    ! "file_name" and accompanying XDMF file.
    ! Optionally,  the (curvilinear) grid can  be provided through "X",  "Y" and
    ! "Z".
    ! Additionally, the total grid size has  to be provided in "tot_dim", and if
    ! the routine  is called in  parallel, the  group dimensions and  offsets as
    ! well, in "grp_dim" and "grp_offset". 
    ! Optionally,  the  dimension which  corresponds  to  time can  be  provided
    ! (default: 4)  in "time_id"  (ignored if outside  range (0..1)).  Also, the
    ! time series  can be made into  an animation with "anim"  and a description
    ! can be provided.
    subroutine print_HDF5_3D_arr(var_name,file_name,vars,tot_dim,grp_dim,&
        &grp_offset,X,Y,Z,time_id,anim,description)                             ! array version
        use HDF5_vars, only: open_HDF5_file, add_HDF5_item, print_HDF5_top, &
            &print_HDF5_geom, print_HDF5_3D_data_item, print_HDF5_att, &
            &print_HDF5_grid, close_HDF5_file, &
            &XML_str_type, HDF5_file_type
        
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of variable to be plot
        character(len=*), intent(in) :: file_name                               ! file name
        real(dp), intent(in), target :: vars(:,:,:,:)                           ! variables to plot
        integer, intent(in) :: tot_dim(4)                                       ! total dimensions of the arrays
        integer, intent(in), optional :: grp_dim(4)                             ! group dimensions of the arrays
        integer, intent(in), optional :: grp_offset(4)                          ! offset of group dimensions
        real(dp), intent(in), target, optional :: X(:,:,:,:)                    ! curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:,:)                    ! curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:,:)                    ! curvlinear grid Z points
        integer, intent(in), optional :: time_id                                ! index of time dimension
        logical, intent(in), optional :: anim                                   ! .true. if animation
        character(len=*), intent(in), optional :: description                   ! description
        
        ! local variables
        integer :: istat                                                        ! status
        integer :: time_id_loc                                                  ! local copy of time_id
        integer :: anim_loc                                                     ! 2 if animation, 3 if not
        type(HDF5_file_type) :: file_info                                       ! file info
        integer :: n_time                                                       ! nr. of points in time
        integer :: id, jd                                                       ! counter
        integer :: tot_dim_3D(3)                                                ! tot_dim except time
        integer :: grp_dim_3D(3)                                                ! grp_dim except time
        integer :: grp_offset_3D(3)                                             ! grp_offset except time
        type(XML_str_type) :: time_col_grid                                     ! grid with time collection
        type(XML_str_type), allocatable :: grids(:)                             ! the grids in the time collection
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type) :: XYZ(3)                                            ! data items for geometry
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        logical :: time_mask(4) = .false.                                       ! to select out the time dimension
        real(dp), pointer :: var_ptr(:,:,:)                                     ! pointer to vars, X, Y or z
        character(len=max_str_ln) :: full_var_name                              ! full variable name
        
        
        ! set up local time_id and anim
        if (present(time_id)) then
            time_id_loc = time_id
        else
            time_id_loc = 4                                                     ! default time dimension: last index
        end if
        anim_loc = 3                                                            ! default spatial collection (no animation)
        if (present(anim)) then
            if(anim) anim_loc = 2                                               ! time collection
        end if
        
        ! set real dimensions
        time_mask(time_id_loc) = .true.
        tot_dim_3D = pack(tot_dim,.not.time_mask)
        if (present(grp_dim)) then
            grp_dim_3D = pack(grp_dim,.not.time_mask)
        else
            grp_dim_3D = tot_dim_3D
        end if
        if (present(grp_offset)) then
            grp_offset_3D = pack(grp_offset,.not.time_mask)
        else
            grp_offset_3D = 0
        end if
        
        ! set up nr. of plots
        n_time = size(vars,time_id_loc)
        
        ! open HDF5 file
        istat = open_HDF5_file(file_info,file_name,description)
        CHCKSTT
        
        ! create grid for time collection
        allocate(grids(n_time))
        
        ! loop over all points in time
        do id = 1,n_time
            ! print topology
            call print_HDF5_top(top,2,tot_dim_3D)
            
            ! print data item for X
            if (present(X)) then
                if (time_id_loc.eq.1) then
                    var_ptr => X(id,:,:,:)
                else if (time_id_loc.eq.2) then
                    var_ptr => X(:,id,:,:)
                else if (time_id_loc.eq.3) then
                    var_ptr => X(:,:,id,:)
                else if (time_id_loc.eq.4) then
                    var_ptr => X(:,:,:,id)
                else
                    istat = 1
                    CHCKSTT
                end if
            else
                allocate(var_ptr(grp_dim_3D(1),grp_dim_3D(2),grp_dim_3D(3)))
                do jd = 1,grp_dim_3D(1)
                    var_ptr(jd,:,:) = jd
                end do
            end if
            istat = print_HDF5_3D_data_item(XYZ(1),file_info,&
                &'X_'//trim(i2str(id)),var_ptr,tot_dim_3D,grp_dim_3D,&
                &grp_offset_3D)
            CHCKSTT
            
            ! print data item for Y
            if (present(Y)) then
                if (time_id_loc.eq.1) then
                    var_ptr => Y(id,:,:,:)
                else if (time_id_loc.eq.2) then
                    var_ptr => Y(:,id,:,:)
                else if (time_id_loc.eq.3) then
                    var_ptr => Y(:,:,id,:)
                else if (time_id_loc.eq.4) then
                    var_ptr => Y(:,:,:,id)
                else
                    istat = 1
                    CHCKSTT
                end if
            else
                allocate(var_ptr(grp_dim_3D(1),grp_dim_3D(2),grp_dim_3D(3)))
                do jd = 1,grp_dim_3D(2)
                    var_ptr(:,jd,:) = jd
                end do
            end if
            istat = print_HDF5_3D_data_item(XYZ(2),file_info,&
                &'Y_'//trim(i2str(id)),var_ptr,tot_dim_3D,grp_dim_3D,&
                &grp_offset_3D)
            CHCKSTT
            
            ! print data item for Z
            if (present(Z)) then
                if (time_id_loc.eq.1) then
                    var_ptr => Z(id,:,:,:)
                else if (time_id_loc.eq.2) then
                    var_ptr => Z(:,id,:,:)
                else if (time_id_loc.eq.3) then
                    var_ptr => Z(:,:,id,:)
                else if (time_id_loc.eq.4) then
                    var_ptr => Z(:,:,:,id)
                else
                    istat = 1
                    CHCKSTT
                end if
            else
                allocate(var_ptr(grp_dim_3D(1),grp_dim_3D(2),grp_dim_3D(3)))
                do jd = 1,grp_dim_3D(3)
                    var_ptr(:,:,jd) = jd
                end do
            end if
            istat = print_HDF5_3D_data_item(XYZ(3),file_info,&
                &'Z_'//trim(i2str(id)),var_ptr,tot_dim_3D,grp_dim_3D,&
                &grp_offset_3D)
            CHCKSTT
            
            ! print geometry with X, Y and Z data item
            call print_HDF5_geom(geom,2,XYZ,.true.)
            
            ! print data item for plot variable
            if (time_id_loc.eq.1) then
                var_ptr => vars(id,:,:,:)
            else if (time_id_loc.eq.2) then
                var_ptr => vars(:,id,:,:)
            else if (time_id_loc.eq.3) then
                var_ptr => vars(:,:,id,:)
            else if (time_id_loc.eq.4) then
                var_ptr => vars(:,:,:,id)
            else
                istat = 1
                CHCKSTT
            end if
            istat = print_HDF5_3D_data_item(XYZ(1),file_info,'var_'//&
                &trim(i2str(id)),var_ptr,tot_dim_3D,grp_dim_3D,grp_offset_3D)   ! reuse XYZ(1)
            CHCKSTT
            
            ! print attribute with this data item
            call print_HDF5_att(att(1),XYZ(1),var_name,1,.true.)
            
            ! create a grid with the topology, the geometry and the attribute
            if (anim_loc.eq.2) then                                             ! time collection
                full_var_name = var_name
            else                                                                ! spatial collection: need different names for grids
                full_var_name = var_name//'_'//trim(i2str(id))
            end if
            istat = print_HDF5_grid(grids(id),full_var_name,1,&
                &grid_time=id*1._dp,grid_top=top,grid_geom=geom,&
                &grid_atts=att,reset=.true.)
            CHCKSTT
        end do
        
        ! create grid collection from individual grids and reset them
        istat = print_HDF5_grid(time_col_grid,'time collection',anim_loc,&
            &grid_grids=grids,reset=.true.)
        CHCKSTT
        
        ! add collection grid to HDF5 file and reset it
        call add_HDF5_item(file_info,time_col_grid,reset=.true.)
        
        ! close HDF5 file
        istat = close_HDF5_file(file_info)
        CHCKSTT
    end subroutine print_HDF5_3D_arr
    subroutine print_HDF5_3D_ind(var_name,file_name,var,tot_dim,grp_dim,&
        &grp_offset,X,Y,Z,description)                                          ! individual version
        use HDF5_vars, only: open_HDF5_file, add_HDF5_item, &
            &XML_str_type, HDF5_file_type, print_HDF5_top, print_HDF5_geom, &
            &print_HDF5_3D_data_item, print_HDF5_att, print_HDF5_grid, &
            &close_HDF5_file
        
        ! input / output
        character(len=*), intent(in) :: var_name                                ! name of variable to be plot
        character(len=*), intent(in) :: file_name                               ! file name
        real(dp), intent(in), target :: var(:,:,:)                              ! variable to plot
        integer, intent(in) :: tot_dim(3)                                       ! total dimensions of the arrays
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
        type(XML_str_type) :: grid                                              ! grid
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type) :: XYZ(3)                                            ! data items for geometry
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        real(dp), pointer :: var_ptr(:,:,:)                                     ! pointer to vars, X, Y or z
        
        ! set up file info
        file_info%name = file_name
        
        ! open HDF5 file
        istat = open_HDF5_file(file_info,description)
        CHCKSTT
        
        ! print topology
        call print_HDF5_top(top,2,tot_dim)
            
        ! print data item for X
        if (present(X)) then
            var_ptr => X
        else
            allocate(var_ptr(size(var,1),size(var,2),size(var,3)))
            do id = 1,size(var,1)
                var_ptr(id,:,:) = id
            end do
        end if
        istat = print_HDF5_3D_data_item(XYZ(1),file_info,'X',var_ptr,&
            &tot_dim,grp_dim=grp_dim,grp_offset=grp_offset)
        CHCKSTT
        
        ! print data item for Y
        if (present(Y)) then
            var_ptr => Y
        else
            allocate(var_ptr(size(var,1),size(var,2),size(var,3)))
            do id = 1,size(var,1)
                var_ptr(:,id,:) = id
            end do
        end if
        istat = print_HDF5_3D_data_item(XYZ(2),file_info,'Y',var_ptr,&
            &tot_dim,grp_dim=grp_dim,grp_offset=grp_offset)
        CHCKSTT
        
        ! print data item for Z
        if (present(Z)) then
            var_ptr => Z
        else
            allocate(var_ptr(size(var,1),size(var,2),size(var,3)))
            do id = 1,size(var,1)
                var_ptr(:,:,id) = id
            end do
        end if
        istat = print_HDF5_3D_data_item(XYZ(3),file_info,'Z',var_ptr,&
            &tot_dim,grp_dim=grp_dim,grp_offset=grp_offset)
        CHCKSTT
        
        ! print geometry with X, Y and Z data item
        call print_HDF5_geom(geom,2,XYZ,.true.)
        
        ! print data item for plot variable
        istat = print_HDF5_3D_data_item(XYZ(1),file_info,var_name,&
            &var,tot_dim,grp_dim=grp_dim,grp_offset=grp_offset)             ! reuse XYZ(1)
        CHCKSTT
        
        ! print attribute with this data item
        call print_HDF5_att(att(1),XYZ(1),var_name,1,.true.)
        
        ! create a grid with the topology, the geometry and the attribute
        istat = print_HDF5_grid(grid,var_name,1,grid_top=top,&
            &grid_geom=geom,grid_atts=att,reset=.true.)
        CHCKSTT
        
        ! add grid to HDF5 file and reset it
        call add_HDF5_item(file_info,grid,reset=.true.)
        
        ! close HDF5 file
        istat = close_HDF5_file(file_info)
        CHCKSTT
    end subroutine print_HDF5_3D_ind
end module output_ops
