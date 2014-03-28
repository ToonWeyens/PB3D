!------------------------------------------------------------------------------!
!   This module contains operations concerning giving output, on the screen as !
!   well as in output files                                                    !
!------------------------------------------------------------------------------!
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
        
        ! input / output
        integer, intent(in) :: nx, ny
        character(len=*) :: fun_name
        real(dp) :: fun(1:nx,1:ny)
        integer, intent(in), optional :: alt_output_i
        character(len=*), optional :: comment
        
        ! local variables
        integer :: fin_output_i
        real(dp), allocatable :: fun_alt(:,:)
        integer :: id
        
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
            case (3)
                if (nx.eq.1) then                                               ! DISLIN 2D without x-vector
                    ! allocate a new variable fun_alt that holds a standard-axis
                    ! ax nd the original data in the y-axis
                    allocate(fun_alt(2,size(fun,2))); fun_alt = 0.0_dp
                    fun_alt(1,:) = (/ ( id , id = 1,size(fun,2)) /)
                    fun_alt(2,:) = fun(1,:)
                    call write_out_dislin_2D(ny, fun_alt, fun_name, comment)
                else if (nx.eq.2) then                                          ! DISLIN 2D
                    call write_out_dislin_2D(ny, fun, fun_name, comment)
                else if (nx.gt.2) then                                          ! DISLIN 3D
                    call write_out_dislin_3D(nx, ny, fun, fun_name, comment)
                else
                    call writo('WARNING: Dimension 1 has to be at least 1')
                end if
            case default 
                call writo('WARNING: No output format associated with ' &
                    &// i2str(format_out) )
        end select
    contains
        ! writes a 2D array into matlab format
        ! Called only by write_out!
        subroutine write_out_matlab(output_i, nx, ny, fun, fun_name, comment)
            ! input / output
            integer, intent(in) :: output_i, nx, ny
            character(len=*) :: fun_name
            real(dp) :: fun(1:nx,1:ny)
            character(len=*), optional :: comment
            
            ! local variables
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
        
        ! writes a 2D array onto the screen using DISLIN
        subroutine write_out_dislin_2D(np, fun, fun_name, comment)
            use dislin
            
            ! input / output
            integer, intent(in) :: np
            character(len=*) :: fun_name
            real(dp) :: fun(2,np)
            character(len=*), optional :: comment
            
            ! local variables
            integer :: ic
            real(dp) :: t_axis(4), f_axis(4)                                    ! axis that envelops minimum, maximum of data
            real(dp) :: x_val(np), y_val(np)                                    ! to avoid annoying compiler performance warnings
            
            call metafl('xwin')                                                 ! xWin terminal
            !call window(0,0,1024,512)                                           ! set window size
            call setpag('DA4L')                                                 ! A4 landscape
            !call sclfac(0.4_dp)                                                 ! scale factor
            call disini                                                         ! start DISFIN
            call winkey('return')                                               ! return key also exits
            call errmod('all','off')                                            ! disable all messages
            call errmod('warnings','on')                                        ! turn warnings on
            
            call pagfll(255)                                                    ! white background
            call setrgb(0.0_dp,0.0_dp,0.0_dp)                                   ! black foreground for rest
            call pagera()                                                       ! plots page border
            call complx()                                                       ! sets complex font
            !call axspos(451,1800)                                               ! position of lower left corner of axis
            !call axslen(2200,1200)                                              ! length of axis

            call name('x','x')                                                  ! name and label of x axis
            call name('f','y')                                                  ! name and label of y axis
            call ticks(5,'x')
            call ticks(5,'y')

            call titlin ('f(x) = '//fun_name , 1)                               ! main title
            call titlin (comment, 3)                                            ! subtitle
            
            ic=intrgb(0.95_dp,0.95_dp,0.95_dp)                                  ! light grey in RGB
            call axsbgd(ic)                                                     ! background color for axis
            
            ! min and max of axis
            t_axis = det_axis(fun(1,:),10)
            f_axis = det_axis(fun(2,:),10)
            call labdig(max(0,-floor(log10(t_axis(4)))),'x')                    ! number of digits after decimal point in lables
            call labdig(max(0,-floor(log10(f_axis(4)))),'y')                    ! number of digits after decimal point in lables
            
            call graf(t_axis(1),t_axis(2),t_axis(3),t_axis(4),&
                &f_axis(1),f_axis(2),f_axis(3),f_axis(4))
            call title
        
            call color('red')                                                   ! red color
            call thkcrv(5)
            x_val = fun(1,:); y_val = fun(2,:)
            call curve(x_val,y_val,np)                                          ! plot 2D curve
            call thkcrv(1)
            
            call color('black')                                                 ! black color
            !call axgit                                                          ! plot axis
            call dash                                                           ! dashed line style
            call xaxgit                                                         ! plot only x axis
            call solid                                                          ! solid line style
            !call yaxgit                                                         ! plot only y axis
            
            call setrgb(0.5_dp,0.5_dp,0.5_dp)                                   ! dark grey foreground for grid
            call grid(1,1)                                                      ! 1 line per 'ticks' for x-axis and also for y
            call setrgb(0.0_dp,0.0_dp,0.0_dp)                                   ! black foreground for rest
            
            call disfin                                                         ! terminate DISFIN
        end subroutine
        
        ! writes a surfrace plot onto the screen using DISLIN
        subroutine write_out_dislin_3D(nx, ny, fun, fun_name, comment)
            use dislin
            
            ! input / output
            integer, intent(in) :: nx, ny
            character(len=*) :: fun_name
            real(dp) :: fun(1:nx,1:ny)
            character(len=*), optional :: comment
            
            ! local variables
            integer :: ic, id
            real(dp) :: f_axis(4)                                               ! axis that envelops minimum, maximum of data
            real(dp) :: x_arr(nx), y_arr(ny)
            real(dp) :: minx, maxx, miny, maxy
            real(dp) :: fun_joined(1:nx*ny)
            !real(dp) :: zlev                                                    ! for surface plots
        
            
            call metafl('xwin')                                                 ! xWin terminal
            !call window(0,0,1024,512)                                           ! set window size
            call setpag('DA4L')                                                 ! A4 landscape
            !call sclfac(0.4_dp)                                                 ! scale factor
            call disini                                                         ! start DISFIN
            call winkey('return')                                               ! return key also exits
            call errmod('all','off')                                            ! disable all messages
            call errmod('warnings','on')                                        ! turn warnings on
            
            call pagfll(255)                                                    ! white background
            call setrgb(0.0_dp,0.0_dp,0.0_dp)                                   ! black foreground for rest
            call pagera()                                                       ! plots page border
            call complx()                                                       ! sets complex font
            !call axspos(451,1800)                                               ! position of lower left corner of axis
            !call axslen(2200,1200)                                              ! length of axis

            call name('x','x')                                                  ! name and label of x axis
            call name('y','y')                                                  ! name and label of y axis
            call name('f','z')                                                  ! name and label of z axis
            call ticks(5,'xyz')

            call titlin ('f = '//fun_name , 1)                                  ! main title
            call titlin (comment, 3)                                            ! subtitle
            
            ic=intrgb(0.95_dp,0.95_dp,0.95_dp)                                  ! light grey in RGB
            call axsbgd(ic)                                                     ! background color for axis
            
            ! min and max of axis
            minx = 1; maxx = size(fun,1)
            miny = 1; maxy = size(fun,2)
            do id = 1,ny
                fun_joined((id-1)*nx+1:id*nx) = fun(:,id)
            end do
            f_axis = det_axis(fun_joined,10)
            call labdig(-1,'xy')
            call labdig(max(0,-floor(log10(f_axis(4)))),'z')                    ! number of digits after decimal point in lables
            
            call graf3D(minx,maxx,minx,(maxx-minx)/6,miny,maxy,miny,&           ! (shaded) surface plot
                &(maxy-miny)/6,f_axis(1),f_axis(2),f_axis(3),f_axis(4))
            !call graf(minx,maxx,minx,(maxx-minx)/6,miny,maxy,miny,&             ! contour plot
                 !&(maxy-miny)/6)
            call box3d                                                          ! 3D box
            call title
            
            x_arr = (/( id, id = 1,nx)/)
            y_arr = (/( id, id = 1,ny)/)
            call color('red')                                                   ! red color
            call thkcrv(5)
            call surmat(fun,nx,ny,1,1)                                          ! surface plot
            !call shdmod('smooth','sufrace')                                    ! shaded surface plot
            !call surshd(x_arr,nx,y_arr,ny,fun)
            !do id = 1,9                                                         ! contour plot
                !zlev = minf + (id-1)*(maxf-minf)/10
                !call setclr(id*25)
                !if(id.eq.5) then
                  !call labels('none','contur')
                !else
                  !call labels('float','contur')
                !end if
                !call contur(x_arr,nx,y_arr,ny,fun,zlev)
            !end do
            call thkcrv(1)
            call color('black')                                                 ! black color
            
            call disfin                                                         ! terminate DISFIN
        end subroutine
        
        ! determine the best axis for a data range of the form:
        ! (from http://www2.mps.mpg.de/dislin/kap4.html)
        !   det_axis(1) = minimum
        !   det_axis(2) = maximum
        !   det_axis(3) = where to start the first label
        !   det_axis(4) = step betwen labels
        ! (minf,maxf,minf,(maxf-minf)/6)
        function det_axis(var, n_points)
            ! input / output
            integer :: n_points
            real(dp) :: var(:)
            real(dp) :: det_axis(4)
            
            ! local variables
            real(dp) :: delta_var                                               ! delta from data
            integer :: base_exp_x2                                              ! exponent of 2X delta from data
            real(dp) :: aux                                                     ! auxiliary variable
            integer :: aux2                                                     ! auxiliary variable
            real(dp) :: delta_ax                                                ! delta used in axis
            real(dp) :: l_b, u_b                                                ! lower and upper bound
            real(dp) :: margin                                                  ! margin for comparison
            
            ! initialize output
            det_axis = 0.0_dp
            
            ! find the distance between labels on the axis by dividing the range
            ! of  the input  data  (var) by  the number  of  points wanted,  and
            ! rounding to the nearest multiple of 1Ej * 0.5 or 1Ej * 1.0
            ! First find the value of j. E.g.: if delta_var = 0.036 and n_points
            ! = 6,  the value delta_var/n_points  = 0.006  has to be  rounded to
            ! 0.005, which means that j = 2
            ! Then, the nearest minimum and  maximum that are multiples of 0.005
            ! are found, that encompass the data in var
            l_b = minval(var)
            u_b = maxval(var)
            delta_var = u_b - l_b
            aux = 2.0_dp*delta_var/n_points
            base_exp_x2 = floor(log10(aux)) + 1                                 ! so that we get a number of the form 0.abcd
            aux2 = nint(aux*10.0_dp**(-base_exp_x2))
            if (aux2.eq.0) then
                delta_ax = 10.0_dp**(base_exp_x2)*0.1_dp                        ! so that we get 0.1 from aux2/2
            else
                delta_ax = 10.0_dp**(base_exp_x2)*aux2*0.5_dp                   ! so that we get 0.5 from aux2/2
            end if
            
            ! determine the lower and upper  bounds by checking whether there is
            ! a remainder when  dividing the lowest, highest  value by delta_ax.
            margin = 1.0E-5_dp                                                  ! visually indistinguishable
            if (abs(dmod(l_b,delta_ax)).gt.margin) then                         ! take rounded times delta
                det_axis(1) = ( floor(l_b/delta_ax) ) * delta_ax
            else
                det_axis(1) = l_b
            end if
            if (abs(dmod(u_b,delta_ax)).gt.margin) then                         ! take rounded times delta
                det_axis(2) = ( ceiling(u_b/delta_ax) ) * delta_ax
            else
                det_axis(2) = u_b
            end if
            
            ! starting point of axis
            det_axis(3) = det_axis(1)
            
            ! determine the step size that is a multiple of delta_ax, closest to
            ! delta_var/n_points = aux/2 -> step size * delta_ax
            det_axis(4) = nint(0.5_dp*aux/delta_ax)*delta_ax
        end function det_axis
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
