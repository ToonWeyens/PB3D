!------------------------------------------------------------------------------!
!   This module contains operations concerning giving output, on the screen as !
!   well as in output files                                                    !
!------------------------------------------------------------------------------!
module output_ops
    use str_ops, only: i2str, r2str
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public init_output_ops, lvl_ud, writo, print_GP_2D, print_GP_3D, &
        &print_ar_1, print_ar_2, draw_GP, print_err_msg, &
        &lvl, lvl_sep, format_out

    ! global variables
    integer :: lvl                                                              ! lvl determines the indenting. higher lvl = more indenting
    character(len=2) :: lvl_sep = ''                                            ! characters that separate different levels of output
    integer :: format_out
    character(len=5) :: plot_dir = 'Plots'

    ! interfaces
    interface print_GP_2D
        module procedure print_GP_2D_ind, print_GP_2D_arr
    end interface
    interface print_GP_3D
        module procedure print_GP_3D_ind, print_GP_3D_arr
    end interface
    
contains
    ! initialize the variables for the module
    ! [MPI] All ranks
    subroutine init_output_ops
        use num_vars, only: glob_rank
        
        lvl = 1
        
        ! print date and time
        if (glob_rank.eq.0) then
            write(*,*) 'Simulation started on '//get_date()//&
                &', at '//get_time()
            write(*,*) ''
        end if
    end subroutine

    ! prints an error  message that is either user-provided, or  the name of the
    ! calling routine
    subroutine print_err_msg(err_msg,routine_name)
        use num_vars, only: glob_rank
        character(len=*) :: err_msg, routine_name
        
        if (trim(err_msg).eq.'') then
            lvl = 2
            call writo('>> calling routine: '//trim(routine_name)//' of rank '&
                &//trim(i2str(glob_rank)))
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
    end subroutine
    subroutine print_GP_2D_arr(fun_name,file_name_i,y,x,draw)                   ! multiple plots
        use num_vars, only: n_seq_0
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: fun_name
        character(len=*), intent(in) :: file_name_i
        real(dp), intent(in) :: y(1:,1:)
        real(dp), intent(in), optional :: x(1:,1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: nr
        integer :: iplt, ipnt, nplt, npnt
        real(dp), allocatable :: x_fin(:,:)
        integer :: ostat
        character(len=max_str_ln) :: file_name
        
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
            file_name = 'plot_output.dat'
        else
            file_name = file_name_i
        end if
        
        ! open output file
        nr = n_seq_0
        call safe_open(nr,ostat,trim(file_name),'replace','formatted',&
            &delim_in='none')
        
        ! write to output file
        write(nr,*) '# '//trim(fun_name)//':'
        do ipnt = 1,npnt
            write(nr,*) (x_fin(ipnt,iplt), iplt = 1,nplt), &
                &(y(ipnt,iplt), iplt = 1,nplt)
        enddo 
        write(nr,*) ''
        
        ! close output file
        close(nr)
        
        ! if draw is present and equal to .false., cancell calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(fun_name,trim(file_name),nplt,.true.,.true.)
    end subroutine
    
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
    end subroutine
    subroutine print_GP_3D_arr(fun_name,file_name_i,z,y,x,draw)                 ! multiple plots
        use num_vars, only: n_seq_0
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: fun_name
        character(len=*), intent(in) :: file_name_i
        real(dp), intent(in) :: z(1:,1:,1:)
        real(dp), intent(in), optional :: y(1:,1:,1:)
        real(dp), intent(in), optional :: x(1:,1:,1:)
        logical, intent(in), optional :: draw
        
        ! local variables
        integer :: nr
        integer :: iplt, ipntx, ipnty, nplt, npntx, npnty
        real(dp), allocatable :: x_fin(:,:,:)
        real(dp), allocatable :: y_fin(:,:,:)
        integer :: ostat
        character(len=max_str_ln) :: file_name
        
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
            file_name = 'plot_output.dat'
        else
            file_name = file_name_i
        end if
        
        ! open output file
        nr = n_seq_0
        call safe_open(nr,ostat,trim(file_name),'replace','formatted',&
            &delim_in='none')
        
        ! write to output file
        write(nr,*) '# '//trim(fun_name)//':'
        do ipntx = 1,npntx
            do ipnty = 1,npnty
                write(nr,*) (x_fin(ipntx,ipnty,iplt), iplt = 1,nplt), &
                    &(y_fin(ipntx,ipnty,iplt), iplt = 1,nplt), &
                    &(z(ipntx,ipnty,iplt), iplt = 1,nplt)
            end do
            write(nr,*) ''
        enddo 
        write(nr,*) ''
        
        ! close output file
        close(nr)
        
        ! if draw is present and equal to .false., cancell calling draw_GP
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_GP(fun_name,trim(file_name),nplt,.false.,.true.)
    end subroutine
    
    ! use GNUPlot to draw a plot
    ! The variables  fun_name and file_name  hold the name  of the plot  and the
    ! name of the file in which the plot is to be found. is_2D indicates whether
    ! the plot is 2D  or 3D and plot_on_screen indicates whether  the plot is to
    ! be shown on the screen, or to be saved in a file called [fun_name].pdf
    subroutine draw_GP(fun_name,file_name,nplt,is_2D,plot_on_screen)
        use safe_open_mod, only: safe_open
        
        ! input / output
        character(len=*), intent(in) :: file_name
        character(len=*), intent(in) :: fun_name
        integer, intent(in) :: nplt
        logical, intent(in) :: is_2D                                            ! True if 2D, false if 3D
        logical, intent(in) :: plot_on_screen                                   ! True if on screen, false if in file
        
        ! local variables
        character(len=1000*max_str_ln) :: gnuplot_cmd
        integer :: iplt
        
        ! create the GNUPlot command
        if (plot_on_screen) then
            ! basic GNUPlot commands
            !gnuplot_cmd = 'gnuplot -persist -e "set grid; set border 4095 &
                !&front linetype -1 linewidth 1.000; set terminal x11;'
            gnuplot_cmd = 'gnuplot -persist -e "set grid; set border 4095 &
                &front linetype -1 linewidth 1.000; set terminal wxt;'
            ! no legend if too many plots
            if (nplt.gt.10) gnuplot_cmd = trim(gnuplot_cmd) // ' set nokey;'
            ! individual plots
            if (is_2D) then
                gnuplot_cmd = trim(gnuplot_cmd) // ' plot'
                do iplt = 1,nplt
                    gnuplot_cmd = trim(gnuplot_cmd) // ' '''//trim(file_name)//&
                        &''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                        &'_'//trim(i2str(iplt))//''' with lines'
                    if (iplt.lt.nplt) gnuplot_cmd = trim(gnuplot_cmd)//','
                end do
            else
                gnuplot_cmd = trim(gnuplot_cmd) // ' splot'
                do iplt = 1,nplt
                    gnuplot_cmd = trim(gnuplot_cmd) // ' '''//trim(file_name)//&
                        &''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//':'//&
                        &trim(i2str(2*nplt+iplt))//' title '''//&
                        &trim(fun_name)//'_'//trim(i2str(iplt))//''' with lines'
                    if (iplt.lt.nplt) gnuplot_cmd = trim(gnuplot_cmd)//','
                end do
            end if
            ! finishing the GNUPlot command
            gnuplot_cmd = trim(gnuplot_cmd)//'; pause -1"'
        else
            ! basic GNUPlot commands
            gnuplot_cmd = 'gnuplot -e "set grid; set border 4095 front &
                &linetype -1 linewidth 1.000; set terminal pdf; &
                &set output '''//trim(plot_dir)//'/'//trim(fun_name)//&
                &'.pdf'';'
            ! no legend if too many plots
            if (nplt.gt.10) gnuplot_cmd = trim(gnuplot_cmd) // ' set nokey;'
            ! individual plots
            if (is_2D) then
                gnuplot_cmd = trim(gnuplot_cmd) // ' plot'
                do iplt = 1,nplt
                    gnuplot_cmd = trim(gnuplot_cmd) // ' '''//trim(file_name)//&
                        &''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//' title '''//trim(fun_name)//&
                        &'_'//trim(i2str(iplt))//''' with lines'
                    if (iplt.lt.nplt) gnuplot_cmd = trim(gnuplot_cmd)//','
                end do
            else
                gnuplot_cmd = trim(gnuplot_cmd) // ' splot'
                do iplt = 1,nplt
                    gnuplot_cmd = trim(gnuplot_cmd) // ' '''//trim(file_name)//&
                        &''' using '//trim(i2str(iplt))//':'//&
                        &trim(i2str(nplt+iplt))//':'//&
                        &trim(i2str(2*nplt+iplt))//' title '''//&
                        &trim(fun_name)//'_'//trim(i2str(iplt))//''' with lines'
                    if (iplt.lt.nplt) gnuplot_cmd = trim(gnuplot_cmd)//','
                end do
            end if
            ! finishing the GNUPlot command
            gnuplot_cmd = trim(gnuplot_cmd)//';"'
        end if
        !write(*,*) 'gnuplot_cmd = ', trim(gnuplot_cmd)
        
        ! call GNUPlot
        call system(gnuplot_cmd)
    end subroutine
    
    ! write output to optional file number 'file_i' using the correct 
    ! indentation for the level ('lvl_loc') of the output
    ! [MPI] Only masters of groups of alpha and the global master call these
    !       The  global master  outputs to  the  master output  file, while  the
    !       masters  of the  groups of  alpha  write their  output to  different
    !       files, which  are then read  by the  global master when  the group's
    !       work is done
    !       Before the groups are created, the group rank is set to be identical
    !       to  the global  rank. This  way,  the group  rank always  determines
    !       whether a  process outputs  or not,  also when  there are  no groups
    !       (yet)
    subroutine writo(input_str,file_i)
        use num_vars, only: group_rank, glob_rank, output_i
        
        ! input / output
        character(len=*), intent(in) :: input_str                               ! the name that is searched for
        integer, optional, intent(in) :: file_i                                 ! optionally set the number of output file
        
        ! local variables
        character(len=max_str_ln) :: output_str                                 ! the name that is searched for
        character(len=max_str_ln) :: header_str                                 ! the name that is searched for
        integer :: id, i_part, max_len_part, num_parts, st_part, en_part
        integer :: loc_file_i
        
        if (group_rank.eq.0) then                                               ! only group master (= global master if no groups)
            ! set  loc_file_i,  depending on  whether  it is  given, or  whether
            ! non-global group master calling
            if (present(file_i)) then
                loc_file_i = file_i
            else
                if (glob_rank.ne.0) then
                    loc_file_i = output_i
                else
                    loc_file_i = 0
                end if
            end if
            
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
                do id = 1,lvl-1                                                 ! start with lvl 1
                    output_str = lvl_sep // trim(output_str)
                end do
                output_str = get_time()//': '//trim(output_str)
                
                ! construct header string of equal length as output strength
                header_str = ''
                do id = 1, len(trim(output_str)) - 10                           ! 10 for time
                    header_str =  trim(header_str) // '-'
                end do
                header_str = '          '//trim(header_str)
                
                if (lvl.eq.1) write(loc_file_i,*) header_str                    ! first level gets extra lines
                write(loc_file_i,*) output_str
                if (lvl.eq.1) write(loc_file_i,*) header_str
            end do
        end if
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
    function get_time() result(time)
        ! input / output
        character(len=8) :: time                                                ! time
        
        ! local variables
        integer :: now(3)
        
        call itime(now)                                                         ! now(1)=hour, (2)=minute, (3)=second
        
        write (time,'(i2.2,":",i2.2,":",i2.2)')  now
    end function get_time
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
        !integer :: ostat
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
                    !call safe_open(fin_output_i,ostat,'tempoutput.dat',&
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
                !if (write_x11) call system('rm tempoutput.dat')
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
                !call system(gnuplot_cmd)
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
                !call system(gnuplot_cmd)
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
                !call system(gnuplot_cmd)
            !end if
        !end subroutine
    !end subroutine
