!------------------------------------------------------------------------------!
!> Operations  concerning giving  output, on  the screen  as well  as in  output
!! files.
!------------------------------------------------------------------------------!
module output_ops
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, max_str_ln, no_plots, iu, plot_dir, data_dir, &
        &script_dir, plot_size
    use files_utilities, only: nextunit
    use grid_vars, only: grid_type
    
    
    implicit none
    private
    public print_ex_2D, print_ex_3D, draw_ex, plot_HDF5, plot_diff_HDF5
    
    ! global variables
    integer :: temp_id(2) = 1                                                   ! will be appended to temporary data and script files
    character(len=9) :: line_clrs(20) = [&
        &'"#7297E6"','"#67EB84"','"#F97A6D"','"#F9C96D"','"#1D4599"',&
        &'"#11AD34"','"#E62B17"','"#E69F17"','"#2F3F60"','"#2F6C3D"',&
        &'"#8F463F"','"#8F743F"','"#031A49"','"#025214"','"#6D0D03"',&
        &'"#6D4903"','"#A9BDE6"','"#A6EBB5"','"#F9B7B0"','"#F9E0B0"']           ! line colors, from <https://github.com/Gnuplotting/gnuplot-palettes>
    character(len=max_str_ln) :: line_style = 'lt 1 lw 1 pt 7 ps 0.5;'          ! line style without little space in line connecting points
#if ldebug
    character(len=0) :: err_output_str = ''                                     ! string with error output
#else
    character(len=14) :: err_output_str = ' 2> /dev/null'                       ! string with error output (/dev/null)
#endif
    
    ! interfaces
    
    !> \public Print 2-D output on a file.
    !!
    !! The variables \c var_name and \c file_name  hold the name of the plot and
    !! of the file in which the plot data  is to be saved, respectively. \c y is
    !! the an  array containing  the function  which is  stored and  \c x  is an
    !! optional vector  with the  x-values. The logical  \c draw  can optionally
    !! disable  calling the  external  drawing procedure  for  output on  screen
    !! [default], without modifying the plot file.
    !!
    !! The first index of \c y (and \c x) contains the points of a current plot.
    !!
    !! The second index indicates various plots (one or more)
    interface print_ex_2D
        !> \public
        module procedure print_ex_2D_ind
        !> \public
        module procedure print_ex_2D_arr
    end interface
    
    !> \public Print 3-D output on a file.
    !!
    !! The variables \c var_name and \c file_name  hold the name of the plot and
    !! of the file in which the plot data  is to be saved, respectively. \c z is
    !! the an array containing the function which is stored and \c x and\c y are
    !! optional  vectors with  the  x  and y-values.  The  logical  \c draw  can
    !! optionally disable calling  the external drawing procedure  for output on
    !! screen [default], without modifying the plot file.
    !!
    !! The first index of \c z (and \c x, \c y) contains the points of a current
    !! The plot second index indicates various plots (one or more)
    interface print_ex_3D
        !> \public
        module procedure print_ex_3D_ind
        !> \public
        module procedure print_ex_3D_arr
    end interface
    

    !> \public Prints variables \c vars with  names \c var_names in an HDF5 file
    !! with name c file_name and accompanying XDMF file.
    !!
    !! For XDMF collections (\cite xdmf), only  the first value for \c var_names
    !! is used, so it should have a size of one.
    !!
    !! The plot is generally 3-D, but if one of the dimensions provided is equal
    !! to 1, it is checked whether there is poloidal or toroidal axisymmetry and
    !! if so, the plot becomes 2-D. This  can be forced using the optional input
    !! argument \c sym_type.
    !!
    !! Optionally, the (curvilinear) grid can be provided through \c X, \c Y and
    !! \c Z. If not,  the grid is assumed to be  Cartesian with discrete indices
    !! where \c X corresponds to the first dimensions, \c Y to the second and \c
    !! Z to the third.
    !!
    !! Additionally, the total grid size and  local offset can be provided in \c
    !! tot_dim and \c loc_offset to  run this routine in parallel automatically.
    !! Optionally,  one  of  the  dimensions  (\c  col_id,  default  4)  can  be
    !! associated to a collection dimension using \c col different from 1:
    !!  - \c col = 1: no collection, just plots of different variables
    !!  - \c col = 2: time collection
    !!  - \c col = 3: spatial collection
    !!  - \c col = 4: vector field
    !!
    !! Furthermore,   using  the   variable  \c   cont_plot,  a   plot  can   be
    !! (over-)written in multiple writes. By this  is meant that there should be
    !! an initial  plot, with collection type  1, 2, 3  or 4, which can  then be
    !! followed by an arbitrary number of additional writes. As these additional
    !! writes  currently  cannot  modify  the   plot  structure,  nor  the  XDMF
    !! variables, their collection dimension should be complete from the start.
    !!
    !! This has no  implications for single plots but means  that for collection
    !! types 2  to 4  all the  elements in  the collection  have to  be present,
    !! though they do not necessary need  to have been completely written in the
    !! other dimensions.
    !!
    !! Subsequent writes with  \c cont_plot can then, for  instance, write parts
    !! of the data  that had not yet  been written, or overwrite  ones that had.
    !! This can be useful for  post-processing where the memory requirements are
    !! large so that the work has to be split.
    !!
    !! \note
    !!  -# For a vector field, the number of variables has to be 2 for 2-D plots
    !!  or 3  for 3-D  plots. This should  be rewritten in  the future,  so that
    !!  collections can be used for vectors as well.
    !!  -# In order to merge collections in their collection dimension, the XDMF
    !!  files can always easily be joined.
    !!  -# If  necessary, a lock system  should be used when  multiple processes
    !!  are writing the same file, including continued writes.
    !!  -# To plot this with VisIt, use:
    !!      -  for temporal  collections:  pseudocolor using  the variable  name
    !!      (other names are ignored).
    !!      - for spatial collections: subset of blocks or pseudocolor using the
    !!      variable name (other names are ignored).
    !!      - for vector plot: Vector plot.
    !!      - for without collections: pseudocolor using the variable names.
    !!  -# Currently all of possibly multiple  processes that write data have to
    !!  cover the entire range of the variables in the dimension indicated by \c
    !!  col_id. (This could be implemented by changing how n_plot is defined and
    !!  selectively  letting  each processer  write  in  the  main loop  at  its
    !!  corresponding indices.)
    !!  -# To project  the data to 2-D  in VisIt, use the  projection tool under
    !!  Operators > Transform
    interface plot_HDF5
        !> \public
        module procedure plot_HDF5_ind
        !> \public
        module procedure plot_HDF5_arr
    end interface
    
contains
    !> \private individual version
    subroutine print_ex_2D_ind(var_name,file_name_i,y,x,draw,persistent)
        ! input / output
        character(len=*), intent(in) :: var_name                                !< name of variable in legend
        character(len=*), intent(in) :: file_name_i                             !< name of input file
        real(dp), intent(in) :: y(1:)                                           !< ordinate
        real(dp), intent(in), optional :: x(1:)                                 !< absicca
        logical, intent(in), optional :: draw                                   !< whether to draw the plot as well
        logical, intent(in), optional :: persistent                             !< keep on-screen plot open
        
        ! local variables
        integer :: npoints
        
        ! set npoints, ny
        npoints = size(y)
        
        ! call multiplot version
        if (present(x)) then
            call print_ex_2D_arr([var_name],file_name_i,reshape(y,[npoints,1]),&
                &x=reshape(x,[npoints,1]),draw=draw,persistent=persistent)
        else
            call print_ex_2D_arr([var_name],file_name_i,reshape(y,[npoints,1]),&
                &draw=draw,persistent=persistent)
        end if
    end subroutine print_ex_2D_ind
    !> \private array version
    subroutine print_ex_2D_arr(var_names,file_name_i,y,x,draw,persistent)
        use num_vars, only: rank
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            !< names of variables in legend
        character(len=*), intent(in) :: file_name_i                             !< name of input file
        real(dp), intent(in) :: y(1:,1:)                                        !< ordinate
        real(dp), intent(in), optional :: x(1:,1:)                              !< absicca
        logical, intent(in), optional :: draw                                   !< whether to draw the plot as well
        logical, intent(in), optional :: persistent                             !< keep on-screen plot open
        
        ! local variables
        integer :: file_i
        integer :: iplt, ipnt, nplt, npnt
        real(dp), allocatable :: x_fin(:,:)
        integer :: istat
        character(len=max_str_ln) :: file_name
        character(len=max_str_ln) :: write_format
        logical :: persistent_loc                                               ! local persistent
        
        ! set nplt, npnt
        npnt = size(y,1)
        nplt = size(y,2)
        if (present(x)) then
            if (size(x,1).ne.size(y,1) .or. &
                &(size(x,2).ne.1 .and. size(x,2).ne.size(y,2))) then
                write(*,*) 'shape(x)', shape(x)
                write(*,*) 'shape(y)', shape(y)
                call writo('In print_ex_2D, the size of x and y has to be the &
                    &same... Skipping plot',persistent=.true.,warning=.true.)
            end if
        end if
        
        ! set other variables
        persistent_loc = .false.
        if (present(persistent)) persistent_loc = persistent
        
        ! set x
        allocate(x_fin(npnt,nplt))
        if (present(x)) then
            if (size(x,2).eq.1) then                                            ! assume copy
                do iplt = 1,nplt
                    x_fin(:,iplt) = x(:,1)
                end do
            else
                x_fin = x
            end if
        else
            do ipnt = 1,npnt
                x_fin(ipnt,:) = ipnt
            end do
        end if
        
        ! set default file name if empty
        if (trim(file_name_i).eq.'') then
            file_name = 'temp_data_print_ex_2D_'//trim(i2str(rank))//'_'//&
                trim(i2str(temp_id(1)))
            temp_id(1) = temp_id(1)+1
        else
            file_name = trim(file_name_i)
        end if
        
        ! open output file
        open(UNIT=nextunit(file_i),FILE=data_dir//'/'//trim(file_name)//'.dat',&
            &IOSTAT=istat)
        CHCKSTT
        
        ! write to output file
        write(file_i,'(1X,A)',IOSTAT=istat) &
            &'# '//trim(merge_strings(var_names))//':'
        write_format="("//trim(i2str(nplt*2))//"(ES23.16,' '))"
        do ipnt = 1,npnt
            write(file_i,FMT=trim(write_format),IOSTAT=istat) &
                &(x_fin(ipnt,iplt), iplt = 1,nplt), &
                &(y(ipnt,iplt), iplt = 1,nplt)
        enddo 
        write(file_i,'(1X,A)',IOSTAT=istat) ''
        
        ! close output file
        close(file_i)
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! if draw is present and equal to .false., cancel calling draw_ex
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_ex(var_names,file_name,nplt,1,.true.,persistent=persistent)
        
        if (trim(file_name_i).eq.'' .and. .not.persistent_loc) then
            call use_execute_command_line('rm '//data_dir//'/'//&
                &trim(file_name)//'.dat')
        end if
    end subroutine print_ex_2D_arr
    
    !> \private individual version
    subroutine print_ex_3D_ind(var_name,file_name_i,z,y,x,draw)
        ! input / output
        character(len=*), intent(in) :: var_name                                !< name of variable in legend
        character(len=*), intent(in) :: file_name_i                             !< name of input file
        real(dp), intent(in) :: z(1:,1:)                                        !< ordinate
        real(dp), intent(in), optional :: x(1:,1:)                              !< absicca
        real(dp), intent(in), optional :: y(1:,1:)                              !< absicca
        logical, intent(in), optional :: draw                                   !< whether to draw the plot as well
        
        ! local variables
        integer :: npoints(2)
        
        ! set npoints, ny
        npoints = [size(z,1),size(z,2)]
        
        ! call multiplot version
        if (present(y)) then
            if (present(x)) then
                call print_ex_3D_arr([var_name],file_name_i,&
                    &reshape(z,[size(z,1),size(z,2),1]),&
                    &y=reshape(y,[size(z,1),size(z,2),1]),&
                    &x=reshape(x,[size(z,1),size(z,2),1]),draw=draw)
            else
                call print_ex_3D_arr([var_name],file_name_i,&
                    &reshape(z,[size(z,1),size(z,2),1]),&
                    &y=reshape(y,[size(z,1),size(z,2),1]),&
                    &draw=draw)
            end if
        else
            if (present(x)) then
                call print_ex_3D_arr([var_name],file_name_i,&
                    &reshape(z,[size(z,1),size(z,2),1]),&
                    &x=reshape(x,[size(z,1),size(z,2),1]),&
                    &draw=draw)
            else
                call print_ex_3D_arr([var_name],file_name_i,&
                    &reshape(z,[size(z,1),size(z,2),1]),draw=draw)
            end if
        end if
    end subroutine print_ex_3D_ind
    !> \private array version
    subroutine print_ex_3D_arr(var_names,file_name_i,z,x,y,draw)
        use num_vars, only: rank
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            !< names of variables in legend
        character(len=*), intent(in) :: file_name_i                             !< name of input file
        real(dp), intent(in) :: z(1:,1:,1:)                                     !< ordinate
        real(dp), intent(in), optional :: x(1:,1:,1:)                           !< absicca
        real(dp), intent(in), optional :: y(1:,1:,1:)                           !< absicca
        logical, intent(in), optional :: draw                                   !< whether to draw the plot as well
        
        ! local variables
        integer :: file_i
        integer :: iplt, ipntx, ipnty, nplt, npntx, npnty
        real(dp), allocatable :: x_fin(:,:,:)
        real(dp), allocatable :: y_fin(:,:,:)
        integer :: istat
        character(len=max_str_ln) :: file_name
        character(len=max_str_ln) :: write_format
        
        ! set nplt, npnt
        npntx = size(z,1)
        npnty = size(z,2)
        nplt = size(z,3)
        if (present(x)) then
            if (size(x,1).ne.size(z,1) .or. size(x,2).ne.size(z,2) .or. &
                &(size(x,3).ne.1 .and. size(x,3).ne.size(z,3))) then
                write(*,*) 'shape(x)', shape(x)
                write(*,*) 'shape(z)', shape(z)
                call writo('In print_ex, the size of x and z has to be the &
                    &same... Skipping plot',persistent=.true.,warning=.true.)
                return
            end if
        end if
        if (present(y)) then
            if (size(y,1).ne.size(z,1) .or. size(y,2).ne.size(z,2) .or. &
                &(size(y,3).ne.1 .and. size(y,3).ne.size(z,3))) then
                write(*,*) 'shape(y)', shape(y)
                write(*,*) 'shape(z)', shape(z)
                call writo('In print_ex, the size of y and z has to be the &
                    &same... Skipping plot',persistent=.true.,warning=.true.)
                return
            end if
        end if
        
        ! set x
        allocate(x_fin(npntx,npnty,nplt))
        if (present (x)) then
            if (size(x,3).eq.1) then                                            ! assume copy
                do iplt = 1,nplt
                    x_fin(:,:,iplt) = x(:,:,1)
                end do
            else
                x_fin = x
            end if
        else
            do ipntx = 1,npntx
                x_fin(ipntx,:,:) = ipntx
            end do
        end if
        ! set y
        allocate(y_fin(npntx,npnty,nplt))
        if (present (y)) then
            if (size(y,3).eq.1) then                                            ! assume copy
                do iplt = 1,nplt
                    y_fin(:,:,iplt) = y(:,:,1)
                end do
            else
                y_fin = y
            end if
        else
            do ipnty = 1,npnty
                y_fin(:,ipnty,:) = ipnty
            end do
        end if
        
        ! set default file name if empty
        if (trim(file_name_i).eq.'') then
            file_name = 'temp_data_print_ex_3D_'//trim(i2str(rank))//'_'//&
                trim(i2str(temp_id(1)))
            temp_id(1) = temp_id(1)+1
        else
            file_name = trim(file_name_i)
        end if
        
        ! open output file
        open(nextunit(file_i),FILE=data_dir//'/'//trim(file_name)//'.dat',&
            &IOSTAT=istat)
        CHCKSTT
        
        ! write to output file
        write(file_i,'(A)',IOSTAT=istat) &
            &'# '//trim(merge_strings(var_names))//':'
        write_format="("//trim(i2str(nplt*3))//"(ES23.16,' '))"
        do ipntx = 1,npntx
            do ipnty = 1,npnty
                write(file_i,FMT=trim(write_format),IOSTAT=istat) &
                    &(x_fin(ipntx,ipnty,iplt), iplt = 1,nplt), &
                    &(y_fin(ipntx,ipnty,iplt), iplt = 1,nplt), &
                    &(z(ipntx,ipnty,iplt), iplt = 1,nplt)
            end do
            write(file_i,'(A)',IOSTAT=istat) ''
        enddo 
        write(file_i,'(A)',IOSTAT=istat) ''
        
        ! close output file
        close(file_i)
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! if draw is present and equal to .false., cancel calling draw_ex
        if (present(draw)) then
            if (.not.draw) return
        end if
        call draw_ex(var_names,file_name,nplt,2,.true.)
        
        if (trim(file_name_i).eq.'') then
            call use_execute_command_line('rm '//data_dir//'/'//&
                &trim(file_name)//'.dat')
        end if
    end subroutine print_ex_3D_arr
    
    !> \private array version
    subroutine plot_HDF5_arr(var_names,file_name,vars,tot_dim,loc_offset,&
        &X,Y,Z,col_id,col,sym_type,cont_plot,descr)
        use HDF5_ops, only: open_HDF5_file, add_HDF5_item, print_HDF5_top, &
            &print_HDF5_geom, print_HDF5_3D_data_item, print_HDF5_att, &
            &print_HDF5_grid, close_HDF5_file, merge_HDF5_3D_data_items
        use HDF5_vars, only: dealloc_XML_str, &
            &XML_str_type, HDF5_file_type
        use num_vars, only: n_procs
        use MPI_utilities, only: get_ser_var
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            !< names of variable to be plot
        character(len=*), intent(in) :: file_name                               !< file name
        real(dp), intent(in), target :: vars(:,:,:,:)                           !< variables to plot
        integer, intent(in), optional :: tot_dim(4)                             !< total dimensions of the arrays
        integer, intent(in), optional :: loc_offset(4)                          !< offset of local dimensions
        real(dp), intent(in), target, optional :: X(:,:,:,:)                    !< curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:,:)                    !< curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:,:)                    !< curvlinear grid Z points
        integer, intent(in), optional :: col_id                                 !< index of time dimension
        integer, intent(in), optional :: col                                    !< whether a collection is made
        integer, intent(in), optional :: sym_type                               !< type of symmetry (1: no symmetry, 2: toroidal, 3: poloidal)
        logical, intent(in), optional :: cont_plot                              !< continued plot
        character(len=*), intent(in), optional :: descr                         !< description
        
        ! local variables
        type(HDF5_file_type) :: file_info                                       ! file info
        integer :: istat                                                        ! status
        integer :: col_id_loc                                                   ! local copy of col_id
        integer :: col_loc                                                      ! local copy of col
        integer :: n_plot                                                       ! nr. of plots
        integer :: id, jd                                                       ! counter
        integer :: sym_type_loc                                                 ! local sym_type
        integer :: sym_pol, sym_tor                                             ! used to determine symmetry
        integer :: att_id                                                       ! index in att
        integer, allocatable :: tot_sym_pol(:), tot_sym_tor(:)                  ! sym_pol and sym_tor for all processes
        integer :: tot_dim_loc(4)                                               ! local copy of tot_dim
        integer :: loc_offset_loc(4)                                            ! local copy of loc_offset
        integer :: tot_dim_3D(3)                                                ! tot_dim except collection
        integer :: loc_dim_3D(3)                                                ! loc_dim except collection
        integer :: loc_offset_3D(3)                                             ! loc_offset except collection
        type(XML_str_type) :: col_grid                                          ! grid with collection
        type(XML_str_type), allocatable :: grids(:)                             ! the grids in the time collection
        type(XML_str_type) :: top                                               ! topology
        type(XML_str_type), allocatable :: dat(:)                               ! data items for geometry and attribute
        type(XML_str_type) :: dat_vec                                           ! data items for merged attribute
        type(XML_str_type) :: geom                                              ! geometry
        type(XML_str_type) :: att(1)                                            ! attribute
        logical :: col_mask(4)                                                  ! to select out the collection dimension
        logical :: ind_plot                                                     ! individual plot
        logical :: cont_plot_loc                                                ! local cont_plot
        logical :: first_vc, last_vc                                            ! first or last vector components (for the others, some XDMF operations can be skipped)
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
            call writo('In plot_HDF5, all the processes need to have the full &
                &range in the dimension given by col_id',persistent=.true.,&
                &warning=.true.)
            CHCKSTT
        end if
        if (n_plot.eq.1 .and. col_loc.ne.1) then
            istat = 1
            call writo('In plot_HDF5, if single plot, the collection type &
                &needs to be one',persistent=.true.,warning=.true.)
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
        
        ! set up continued plot
        cont_plot_loc = .false.
        if (present(cont_plot)) cont_plot_loc = cont_plot
        
        ! default symmetry type
        sym_type_loc = 1
        
        ! Find  symmetry type  by  checking whether  Y/X  is constant  (toroidal
        ! symmetry) or  Z^2/(X^2+Y^2) is  constant (poloidal symmetry),  for all
        ! plots.
        if (present(sym_type)) then
            sym_type_loc = sym_type
        else if (minval(tot_dim_3D).eq.1) then                                  ! possibly symmetry
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
                sym_ang = atan2(Y_3D,X_3D)
                if (maxval(sym_ang)-minval(sym_ang).lt.tol_sym .and. &
                    &(maxval(X_3D).ge.0._dp .neqv. &
                    &minval(X_3D).lt.0._dp).and.&                               ! X has to be either positive or negative
                    &(maxval(Y_3D).ge.0._dp .neqv. &
                    &minval(Y_3D).lt.0._dp)) &                                  ! Y has to be either positive or negative
                    &sym_pol = sym_pol+1                                        ! poloidal symmetry for this plot
                ! check toroidal angle
                sym_ang = atan2(sqrt(Z_3D**2),sqrt(X_3D**2+Y_3D**2))
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
                sym_type_loc = 2
            else if (sum(tot_sym_tor).eq.n_procs*n_plot) then                   ! toroidal symmetry for all plots
                sym_type_loc = 3
            end if
            ! deallocate helper variables
            deallocate(sym_ang)
        end if
        
        ! tests
        if (col_loc.eq.4) then
            if (sym_type_loc.eq.1) then
                if (n_plot.ne.3) then
                    istat = 1
                    call writo('For vector field plots, need 3 dimensions',&
                        &warning=.true.)
                    CHCKSTT
                end if
            else
                if (n_plot.ne.2) then
                    istat = 1
                    call writo('For symmetric vector field plots, need 2 &
                        &dimensions',warning=.true.)
                    CHCKSTT
                end if
            end if
        end if
        
        ! set up and local var names if not a continued plot
        if (.not.cont_plot_loc) then
            ! set up local var_names
            allocate(grd_names(n_plot))
            allocate(att_names(n_plot))
            if (col_loc.eq.1) then                                              ! without collection
                if (n_plot.eq.1) then                                           ! just one plot: attribute name is important
                    if (size(var_names).eq.n_plot) then                         ! the right number of variable names provided
                        att_names = var_names
                    else if (size(var_names).gt.n_plot) then                    ! too many variable names provided
                        att_names = var_names(1:n_plot)
                        call writo('Too many variable names provided',&
                            &persistent=.true.,warning=.true.)
                    else                                                        ! not enough variable names provided
                        att_names(1:size(var_names)) = var_names
                        do id = size(var_names)+1,n_plot
                            att_names(id) = 'unnamed variable '//trim(i2str(id))
                        end do
                        call writo('Not enough variable names provided',&
                            &persistent=.true.,warning=.true.)
                    end if
                    grd_names = 'default_grid_name'
                else                                                            ! multiple plots: grid name is important
                    if (size(var_names).eq.n_plot) then                         ! the right number of variable names provided
                        grd_names = var_names
                    else if (size(var_names).gt.n_plot) then                    ! too many variable names provided
                        grd_names = var_names(1:n_plot)
                        call writo('Too many variable names provided',&
                            &persistent=.true.,warning=.true.)
                    else                                                        ! not enough variable names provided
                        grd_names(1:size(var_names)) = var_names
                        do id = size(var_names)+1,n_plot
                            grd_names(id) = 'unnamed variable '//trim(i2str(id))
                        end do
                        call writo('Not enough variable names provided',&
                            &persistent=.true.,warning=.true.)
                    end if
                    att_names = 'default_att_name'
                end if
            else                                                                ! collections: attribute name is important
                att_names = var_names(1)
                if (size(var_names).gt.1) call writo('For collections, only &
                    &the first variable name is used',persistent=.true.,&
                    &warning=.true.)
                grd_names = 'default_grid_name'
            end if
        end if
        
        ! open HDF5 file
        istat = open_HDF5_file(file_info,file_name,sym_type=sym_type_loc,&
            &descr=descr,ind_plot=ind_plot,cont_plot=cont_plot_loc)
        CHCKSTT
        
        ! create grid for collection if not continued plot
        if (.not.cont_plot_loc) allocate(grids(n_plot))
        
        ! allocate data item arrays
        if (sym_type_loc.eq.1) then                                             ! 3D grid
            allocate(dat(3))
        else                                                                    ! 2D grid
            allocate(dat(2))
        end if
        
        ! loop over all plots
        do id = 1,n_plot
            ! set up whether last vector component
            ! (for scalar plots, this is always true)
            first_vc = (id.eq.1 .or. col_loc.ne.4)
            last_vc = (id.eq.n_plot .or. col_loc.ne.4)
            
            ! print topology if not continued plot and last vector components
            if (.not.cont_plot_loc .and. last_vc) then
                if (sym_type_loc.eq.1) then                                     ! 3D grid
                    call print_HDF5_top(top,2,tot_dim_3D,ind_plot=ind_plot)
                else                                                            ! 2D grid
                    call print_HDF5_top(top,1,tot_dim_3D,ind_plot=ind_plot)
                end if
            end if
            
            ! assign pointers
            call assign_pointers(id)
            
            ! print data  item for X, Y  and Z (no symmetry),  R = sqrt(X^2+Y^2)
            ! and Z (poloidal symmetry) or X  and Y (toroidal symmetry) if first
            ! vector component (afterwards, dat is overwritten)
            if (first_vc) then
                select case (sym_type_loc)
                    case (1)                                                    ! no symmetry
                        istat = print_HDF5_3D_data_item(dat(1),file_info,&
                            &'X_'//trim(i2str(id)),X_3D,tot_dim_3D,loc_dim_3D,&
                            &loc_offset_3D,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        CHCKSTT
                        istat = print_HDF5_3D_data_item(dat(2),file_info,&
                            &'Y_'//trim(i2str(id)),Y_3D,tot_dim_3D,loc_dim_3D,&
                            &loc_offset_3D,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        CHCKSTT
                        istat = print_HDF5_3D_data_item(dat(3),file_info,&
                            &'Z_'//trim(i2str(id)),Z_3D,tot_dim_3D,loc_dim_3D,&
                            &loc_offset_3D,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        CHCKSTT
                    case (2)                                                    ! poloidal symmetry
                        istat = print_HDF5_3D_data_item(dat(1),file_info,&
                            &'R_'//trim(i2str(id)),sqrt(X_3D**2+Y_3D**2),&
                            &tot_dim_3D,loc_dim_3D,loc_offset_3D,&
                            &ind_plot=ind_plot,cont_plot=cont_plot_loc)
                        CHCKSTT
                        istat = print_HDF5_3D_data_item(dat(2),file_info,&
                            &'Z_'//trim(i2str(id)),Z_3D,tot_dim_3D,loc_dim_3D,&
                            &loc_offset_3D,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        CHCKSTT
                    case (3)                                                    ! toroidal symmetry
                        istat = print_HDF5_3D_data_item(dat(1),file_info,&
                            &'X_'//trim(i2str(id)),X_3D,tot_dim_3D,loc_dim_3D,&
                            &loc_offset_3D,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        CHCKSTT
                        istat = print_HDF5_3D_data_item(dat(2),file_info,&
                            &'Y_'//trim(i2str(id)),Y_3D,tot_dim_3D,loc_dim_3D,&
                            &loc_offset_3D,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        CHCKSTT
                    case default                                                ! no symmetry
                        istat = 1
                        call writo('symmetry type '//&
                            &trim(i2str(sym_type_loc))//' not recognized',&
                            &persistent=.true.,warning=.true.)
                        CHCKSTT
                end select
            end if
            
            ! print geometry with X, Y and Z data item if not continued plot and
            ! first vector component
            if (.not.cont_plot_loc .and. first_vc) then
                if (sym_type_loc.eq.1) then                                     ! no symmetry so 3D geometry
                    call print_HDF5_geom(geom,2,dat,reset=.true.,&
                        &ind_plot=ind_plot)
                else                                                            ! symmetry so 2D geometry
                    call print_HDF5_geom(geom,1,dat,reset=.true.,&
                        &ind_plot=ind_plot)
                end if
            end if
            
            ! print data item for plot variable
            if (col_loc.eq.4) then
                att_id = id
            else
                att_id = 1
            end if
            istat = print_HDF5_3D_data_item(dat(att_id),file_info,'var_'//&
                &trim(i2str(id)),var_3D,tot_dim_3D,loc_dim_3D,&
                &loc_offset_3D,ind_plot=ind_plot,cont_plot=cont_plot_loc)
            CHCKSTT
            
            ! print attribute with this data item if not continued plot
            if (.not.cont_plot_loc) then
                if (col_loc.eq.4) then
                    if (last_vc) then
                        call merge_HDF5_3D_data_items(dat_vec,dat,'var_vec',&
                            &tot_dim_3D,reset=.true.,ind_plot=ind_plot,&
                            &cont_plot=cont_plot_loc)
                        call print_HDF5_att(att(1),dat_vec,&
                            &att_names(1),1,2,reset=.true.,ind_plot=ind_plot)
                    end if
                else
                    call print_HDF5_att(att(1),dat(1),&
                        &att_names(id),1,1,reset=.true.,ind_plot=ind_plot)
                end if
            end if
            
            ! create a grid  with the topology, the geometry,  the attribute and
            ! time if  time collection,  if not continued  plot and  last vector
            ! component
            if (.not.cont_plot_loc .and. last_vc) then
                if (col_loc.eq.2) then                                          ! time collection
                    istat = print_HDF5_grid(grids(id),grd_names(id),1,&
                        &grid_time=id*1._dp,grid_top=top,grid_geom=geom,&
                        &grid_atts=att,reset=.true.,ind_plot=ind_plot)
                    CHCKSTT
                else                                                            ! no time collection
                    istat = print_HDF5_grid(grids(id),grd_names(id),1,&
                        &grid_top=top,grid_geom=geom,grid_atts=att,&
                        &reset=.true.,ind_plot=ind_plot)
                    CHCKSTT
                end if
            end if
        end do
        
        ! either create collection or just use individual grids if not continued
        ! plot
        if (.not.cont_plot_loc) then
            if (col_loc.eq.1 .or. col_loc.eq.4) then
                ! add individual grids to HDF5 file and reset them
                do id = 1,n_plot
                    if (col_loc.eq.4 .and. id.lt.n_plot) cycle                  ! only plot last for vector
                    istat = add_HDF5_item(file_info,grids(id),reset=.true.,&
                    &ind_plot=ind_plot)
                    CHCKSTT
                end do
            else
                ! create grid collection from individual grids and reset them
                istat = print_HDF5_grid(col_grid,'domain of mesh',col_loc,&
                    &grid_grids=grids,reset=.true.,ind_plot=ind_plot)
                CHCKSTT
                
                ! add collection grid to HDF5 file and reset it
                istat = add_HDF5_item(file_info,col_grid,reset=.true.,&
                    &ind_plot=ind_plot)
                CHCKSTT
            end if
        end if
        
        ! close HDF5 file
        istat = close_HDF5_file(file_info,ind_plot=ind_plot,&
            &cont_plot=cont_plot_loc)
        CHCKSTT
        
        ! clean up
        nullify(var_3D)
        nullify(X_3D,Y_3D,Z_3D)
        call dealloc_XML_str(dat)
        if (.not.cont_plot_loc) then
            call dealloc_XML_str(grids)
            call dealloc_XML_str(top)
            call dealloc_XML_str(geom)
            call dealloc_XML_str(att(1))
            call dealloc_XML_str(col_grid)
        end if
    contains
        !> \private assigns the 3D subarray pointer variables
        subroutine assign_pointers(id)
            ! input / output
            integer :: id                                                       ! index at which to assing pointer
            
            ! local variables
            integer :: id_loc                                                   ! local id
            
            ! X
            if (present(X)) then
                id_loc = id
                if (col_id_loc.eq.1) then
                    if (size(X,1).eq.1) id_loc = 1
                    X_3D => X(id_loc,:,:,:)
                else if (col_id_loc.eq.2) then
                    if (size(X,2).eq.1) id_loc = 1
                    X_3D => X(:,id_loc,:,:)
                else if (col_id_loc.eq.3) then
                    if (size(X,3).eq.1) id_loc = 1
                    X_3D => X(:,:,id_loc,:)
                else if (col_id_loc.eq.4) then
                    if (size(X,4).eq.1) id_loc = 1
                    X_3D => X(:,:,:,id_loc)
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
                id_loc = id
                if (col_id_loc.eq.1) then
                    if (size(Y,1).eq.1) id_loc = 1
                    Y_3D => Y(id_loc,:,:,:)
                else if (col_id_loc.eq.2) then
                    if (size(Y,2).eq.1) id_loc = 1
                    Y_3D => Y(:,id_loc,:,:)
                else if (col_id_loc.eq.3) then
                    if (size(Y,3).eq.1) id_loc = 1
                    Y_3D => Y(:,:,id_loc,:)
                else if (col_id_loc.eq.4) then
                    if (size(Y,4).eq.1) id_loc = 1
                    Y_3D => Y(:,:,:,id_loc)
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
            ! Z
            if (present(Z)) then
                id_loc = id
                if (col_id_loc.eq.1) then
                    if (size(Z,1).eq.1) id_loc = 1
                    Z_3D => Z(id_loc,:,:,:)
                else if (col_id_loc.eq.2) then
                    if (size(Z,2).eq.1) id_loc = 1
                    Z_3D => Z(:,id_loc,:,:)
                else if (col_id_loc.eq.3) then
                    if (size(Z,3).eq.1) id_loc = 1
                    Z_3D => Z(:,:,id_loc,:)
                else if (col_id_loc.eq.4) then
                    if (size(Z,4).eq.1) id_loc = 1
                    Z_3D => Z(:,:,:,id_loc)
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
    !> \private individual version
    subroutine plot_HDF5_ind(var_name,file_name,var,tot_dim,loc_offset,&
        &X,Y,Z,cont_plot,descr)
        
        ! input / output
        character(len=*), intent(in) :: var_name                                !< name of variable to be plot
        character(len=*), intent(in) :: file_name                               !< file name
        real(dp), intent(in) :: var(:,:,:)                                      !< variable to plot
        integer, intent(in), optional :: tot_dim(3)                             !< total dimensions of the arrays
        integer, intent(in), optional :: loc_offset(3)                          !< offset of local dimensions
        real(dp), intent(in), target, optional :: X(:,:,:)                      !< curvlinear grid X points
        real(dp), intent(in), target, optional :: Y(:,:,:)                      !< curvlinear grid Y points
        real(dp), intent(in), target, optional :: Z(:,:,:)                      !< curvlinear grid Z points
        logical, intent(in), optional :: cont_plot                              !< continued plot
        character(len=*), intent(in), optional :: descr                         !< description
        
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
                        &col=1,cont_plot=cont_plot,descr=descr)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &Y=reshape(Y,[size(Y,1),size(Y,2),size(Y,3),1]),&
                        &col=1,cont_plot=cont_plot,descr=descr)
                end if
            else
                if (present(Z)) then
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &Z=reshape(Z,[size(Z,1),size(Z,2),size(Z,3),1]),&
                        &col=1,cont_plot=cont_plot,descr=descr)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &X=reshape(X,[size(X,1),size(X,2),size(X,3),1]),&
                        &col=1,cont_plot=cont_plot,descr=descr)
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
                        &col=1,cont_plot=cont_plot,descr=descr)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &Y=reshape(Y,[size(Y,1),size(Y,2),size(Y,3),1]),&
                        &col=1,cont_plot=cont_plot,descr=descr)
                end if
            else
                if (present(Z)) then
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &col=1,Z=reshape(Z,[size(Z,1),size(Z,2),size(Z,3),1]),&
                        &cont_plot=cont_plot,descr=descr)
                else
                    call plot_HDF5_arr([var_name],file_name,&
                        &reshape(var,[size(var,1),size(var,2),size(var,3),1]),&
                        &tot_dim_loc,loc_offset_loc,&
                        &col=1,cont_plot=cont_plot,descr=descr)
                end if
            end if
        end if
    end subroutine plot_HDF5_ind
    
    !> Use external program to draw a plot.
    !!
    !!
    !! The external programs used are determined by \c ex_plot_style:
    !!  - \c ex_plot_style = 1: Use GnuPlot for 2-D and 3-D
    !!  - \c ex_plot_style = 1: Use Bokeh in 2-D and Mayavi in 3-D
    !!
    !! The  output  is saved  in  a  file <tt>[draw_name].[ext]</tt>  where  the
    !! extension depends  on the  plotting style  and the  dimension, or  on the
    !! screen, depending on \c plot_on_screen.
    !!
    !! If  not specified  otherwise through  data_name, it  is assumed  that the
    !! datafile to be plot is situated in <tt>[draw_name].dat</tt>.
    !!
    !! The title(s)  of the plot(s)  are provided  through \c var_names  and the
    !! number of plots through \c nplt. If  there are less names than the number
    !! of plots, the first one is taken and appended by the plot number.
    !!
    !! Furthermore, \c draw_dim determines whether the plot is to be 2-D, 3-D or
    !! 2-D slices in 3-D.
    !!
    !! Also, an optional command \c extra_ops can be provided, to provide overal
    !! options, as well as draw_ops that  specifies the line style for the plots
    !! from the file. If \c less draw_ops  are provided than plots, they will be
    !! cycled.
    !!
    !! Finally,  there is  an  option to  provide animated  plots,  but is  only
    !! available in 2-D.  Also, for external plot style 1,  no on-screen view is
    !! possible. Animations can take optional ranges arguments as well as delay.
    !! \note \c About draw_dim:
    !!  -  \c draw_dim  = 1:  2-D  plot; should  be  called with  the output  of
    !!  output_ops.print_ex_2d(), \c nplt should be correctly set.
    !!  -  \c draw_dim  = 2:  3-D  plot; should  be  called with  the output  of
    !!  print_ex_3D. For GNUPlot, \c nplt is  ignored, but not for Mayavi, as it
    !!  needs this information to be able to  reconstruct 2D arrays for X, Y and
    !!  Z.
    !!  - \c  draw_dim = 3: 2-D  plot in 3-D  slices; should be called  with the
    !!  output of output_ops.print_ex_2d(), \c nplt should be correctly set.
    subroutine draw_ex(var_names,draw_name,nplt,draw_dim,plot_on_screen,&
        &ex_plot_style,data_name,draw_ops,extra_ops,is_animated,ranges,delay,&
        &persistent)
        
        use num_vars, only: ex_plot_style_global => ex_plot_style, rank
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            !< name of variables
        character(len=*), intent(in) :: draw_name                               !< name of drawing
        integer, intent(in) :: nplt                                             !< number of plots
        integer, intent(in) :: draw_dim                                         !< 1: 2-D, 2: 3-D, 3: decoupled 3-D
        logical, intent(in) :: plot_on_screen                                   !< True if on screen, false if in file
        integer, intent(in), optional :: ex_plot_style                          !< alternative external plot style
        character(len=*), intent(in), optional :: data_name                     !< name of data file
        character(len=*), intent(in), optional :: draw_ops(:)                   !< drawing options
        character(len=*), intent(in), optional :: extra_ops                     !< extra options
        logical, intent(in), optional :: is_animated                            !< plot is animated
        real(dp), intent(in), optional :: ranges(:,:)                           !< x and y range of animated plot
        integer, intent(in), optional :: delay                                  !< time delay between animated plot frames
        logical, intent(in), optional :: persistent                             !< keep on-screen plot open
        
        ! local variables
        character(len=4) :: ex_ext                                              ! output extension of external program
        character(len=5*max_str_ln) :: cmdmsg                                   ! error message of C system command
        character(len=max_str_ln) :: script_name                                ! name of script, including path
        character(len=max_str_ln) :: ex_prog_name                               ! name of external program
        character(len=max_str_ln) :: data_name_loc                              ! local data name
        character(len=max_str_ln), allocatable :: var_names_loc(:)              ! local variables names
        logical :: is_animated_loc                                              ! local is_animated
        logical :: run_shell_necessary                                          ! whether shell needs to be run
        logical :: persistent_loc                                               ! local persistent
        real(dp), allocatable :: ranges_loc(:,:)                                ! local copy of ranges
        integer :: delay_loc                                                    ! local copy of delay
        integer :: n_draw_ops                                                   ! number of draw options provided
        integer :: istat                                                        ! status of opening a file
        integer :: iplt                                                         ! counter
        integer :: cmd_i                                                        ! file number for script file
        integer :: cmdstat                                                      ! status of C system command
        integer :: ex_plot_style_loc                                            ! local ex_plot_style
        
        ! skip if no plots
        if (no_plots) return
        
        ! set up local data_name, variable names, is_animated and ex_plot_style
        data_name_loc = draw_name
        if (present(data_name)) data_name_loc = data_name
        allocate(var_names_loc(nplt))
        if (size(var_names).eq.nplt) then
            var_names_loc = var_names
        else
            do iplt = 1,nplt
                var_names_loc(iplt) = trim(var_names(1))//' ('//&
                    &trim(i2str(iplt))//' / '//trim(i2str(nplt))//')'
            end do
        end if
        is_animated_loc = .false.
        if (present(is_animated)) is_animated_loc = is_animated
        ex_plot_style_loc = ex_plot_style_global
        if (present(ex_plot_style)) ex_plot_style_loc = ex_plot_style
        persistent_loc = .false.
        if (present(persistent)) persistent_loc = persistent
        
        ! run shell to produce the plot by default
        run_shell_necessary = .true.
        
        ! test whether animation is in 2D and set some variables
        if (is_animated_loc) then
            if (draw_dim.ne.1) then
                call writo('Animations are only possible in 2D',warning=.true.,&
                    &persistent=.true.)
                return
            else
                ! calculate delay [1/100 s]
                if (present(delay)) then
                    delay_loc = delay
                else
                    delay_loc = max(250/nplt,10)
                end if
            end if
        end if
        
        ! create the external program command
        if (plot_on_screen .and. &
            &(ex_plot_style_loc.ne.1 .and. is_animated_loc)) then               ! not for GNUPlot animation
            script_name = trim(script_dir)//'/'//'temp_script_draw_ex_'//&
                &trim(i2str(rank))//'_'//trim(i2str(temp_id(2)))
            temp_id(2) = temp_id(2)+1
        else
            script_name = trim(script_dir)//'/'//trim(draw_name)
        end if
        select case (ex_plot_style_loc)
            case (1)                                                            ! GNUPlot
                script_name = trim(script_name)//'.gnu'
                ex_ext = 'pdf'
                ex_prog_name = 'gnuplot'
            case (2)                                                            ! Bokeh or Mayavi
                script_name = trim(script_name)//'.py'
                if (draw_dim.eq.1) then                                         ! 2D: Bokeh
                    ex_ext = 'html'
                else                                                            ! 3D and 2D slices in 3D: Mayavi
                    ex_ext = 'png'
                end if
                ex_prog_name = 'python'
        end select
        
        ! open script file
        open(nextunit(cmd_i),FILE=trim(script_name),IOSTAT=istat)
        CHCKSTT
        
        ! set number of draw ops
        if (present(draw_ops)) then
            n_draw_ops = size(draw_ops)
        else
            n_draw_ops = 0
        end if
        
        select case (ex_plot_style_loc)
            case (1)                                                            ! GNUPlot
                if (is_animated_loc) then
                    call draw_ex_animated_GNUPlot()
                else
                    call draw_ex_GNUPlot()
                end if
            case (2)                                                            ! Bokeh
                select case (draw_dim)
                    case (1)                                                    ! 2D
                        !!run_shell_necessary = plot_on_screen                    ! only need to run if plot on screen (?)
                        if (is_animated_loc) then
                            call draw_ex_animated_Bokeh()
                        else
                            call draw_ex_Bokeh()
                        end if
                    case (2:3)                                                  ! 3D, 2D slices in 3D
                        call draw_ex_Mayavi()
                end select
        end select
        
        ! closing the command
        close(cmd_i)
        
        ! run shell if necessary
        istat = 0
        cmdstat = 0
        if (run_shell_necessary) &
            &call use_execute_command_line(trim(ex_prog_name)//' "'//&
                &trim(script_name)//'"'//err_output_str,exitstat=istat,&
                &cmdstat=cmdstat,cmdmsg=cmdmsg)
        
        if (istat.ne.0) then
            call writo('Failed to plot '//trim(draw_name)//'.'//&
                &trim(ex_ext),persistent=.true.,warning=.true.)
        else
            if (cmdstat.ne.0) then
                call writo('Failed to plot '//trim(draw_name)//'.'//&
                    &trim(ex_ext),persistent=.true.,warning=.true.)
                call lvl_ud(1)
                call writo('System message: "'//trim(cmdmsg)//'"',&
                    &persistent=.true.)
                if (.not.plot_on_screen) call writo(&
                    &'Try running "'//trim(ex_prog_name)//' "'//&
                    &trim(script_name)//'"'//'" manually',persistent=.true.)
                call lvl_ud(-1)
            else
                if (plot_on_screen .and. run_shell_necessary .and. &
                    &.not.persistent_loc) then
                    call use_execute_command_line('rm '//&
                        &trim(script_name),exitstat=istat,&
                        &cmdstat=cmdstat,cmdmsg=cmdmsg)
                    ! ignore errors
                else
                    call writo('Plot in output file "'//&
                        &trim(plot_dir)//'/'//trim(draw_name)//&
                        &'.'//trim(ex_ext)//'"',persistent=.true.)
                end if
            end if
        end if
    contains
        !> \private GNUPlot version: wxt terminal or pdf output
        subroutine draw_ex_GNUPlot
            ! initialize the script
            write(cmd_i,"(A)",IOSTAT=istat) 'set grid'
            write(cmd_i,"(A)",IOSTAT=istat) 'set border 4095 front linetype -1 &
                &linewidth 1.0'
            if (plot_on_screen) then
                if (persistent_loc) then
                    write(cmd_i,"(A)",IOSTAT=istat) 'set terminal wxt persist'
                else
                    write(cmd_i,"(A)",IOSTAT=istat) 'set terminal wxt'
                end if
            else
                write(cmd_i,"(A)",IOSTAT=istat) 'set terminal pdf size '//&
                    &trim(i2str(plot_size(1)))//','//trim(i2str(plot_size(2)))
                write(cmd_i,"(A)",IOSTAT=istat) 'set output "'//&
                    &trim(plot_dir)//'/'//trim(draw_name)//'.pdf"'
            end if
            
            ! no legend if too many plots
            if (nplt.gt.10) write(cmd_i,"(A)",IOSTAT=istat) 'set nokey;'
            
            ! write extra options
            if (present(extra_ops)) write(cmd_i,"(A)",IOSTAT=istat) &
                &trim(extra_ops)
            
            ! set up line styles
            if (n_draw_ops.eq.0) then
                do iplt = 1,min(nplt,size(line_clrs))
                    write(cmd_i,"(A)",IOSTAT=istat) 'set style line '//&
                        &trim(i2str(iplt))//' lc rgb '//&
                        &line_clrs(iplt)//' '//trim(line_style)
                end do
            end if
            
            ! individual plots
            select case (draw_dim)
                case (1)                                                        ! 2D
                    write(cmd_i,"(A)",IOSTAT=istat) 'plot \'
                    do iplt = 1,nplt
                        write(cmd_i,"(A)",IOSTAT=istat,ADVANCE="no") &
                            &' "'//trim(data_dir)//'/'//trim(data_name_loc)//&
                            &'.dat" using '//trim(i2str(iplt))//':'//&
                            &trim(i2str(nplt+iplt))//' title "'//&
                            &trim(var_names_loc(iplt))//'" '//&
                            &trim(loc_draw_op())
                        if (iplt.eq.nplt) then
                            write(cmd_i,"(A)",IOSTAT=istat) ''
                        else
                            write(cmd_i,"(A)",IOSTAT=istat) ', \'
                        end if
                    end do
                case (2)                                                        ! 3D
                    write(cmd_i,"(A)",IOSTAT=istat) 'splot \'
                    do iplt = 1,nplt
                        write(cmd_i,"(A)",IOSTAT=istat,ADVANCE="no") &
                            &' "'//trim(data_dir)//'/'//trim(data_name_loc)//&
                            &'.dat" using '//trim(i2str(iplt))//':'//&
                            &trim(i2str(nplt+iplt))//':'//&
                            &trim(i2str(2*nplt+iplt))//' title "'//&
                            &trim(var_names_loc(iplt))//&
                            &'" '//trim(loc_draw_op())
                        if (iplt.eq.nplt) then
                            write(cmd_i,"(A)",IOSTAT=istat) ''
                        else
                            write(cmd_i,"(A)",IOSTAT=istat) ', \'
                        end if
                    end do
                case (3)                                                        ! 2D slices in 3D
                    write(cmd_i,"(A)",IOSTAT=istat) 'splot \'
                    do iplt = 1,nplt
                        write(cmd_i,"(A)",IOSTAT=istat,ADVANCE="no") &
                            &' "'//trim(data_dir)//'/'//trim(data_name_loc)//&
                            &'.dat" using ('//trim(i2str(iplt))//'):'//&
                            &trim(i2str(iplt))//':'//trim(i2str(nplt+iplt))//&
                            &' title "'//trim(var_names_loc(iplt))//'" '//&
                            &trim(loc_draw_op())
                        if (iplt.eq.nplt) then
                            write(cmd_i,"(A)",IOSTAT=istat) ''
                        else
                            write(cmd_i,"(A)",IOSTAT=istat) ', \'
                        end if
                    end do
                case default
                    call writo('No draw_dim associated with '//&
                        &trim(i2str(draw_dim)),persistent=.true.)
                    istat = 1
                    CHCKSTT
            end select
            
            ! finishing the command
            write(cmd_i,"(A)",IOSTAT=istat) ''
            if (plot_on_screen) write(cmd_i,"(A)",IOSTAT=istat) 'pause -1'
        end subroutine draw_ex_GNUPlot
        
        !> \private Bokeh version: 2D html output
        ! Interactive checkbox from https://github.com/bokeh/bokeh/issues/3715.
        ! (Could be replaced by http://bokeh.pydata.org/en/latest/docs/
        !  user_guide/interaction/legends.html#userguide-interaction-legends)
        subroutine draw_ex_Bokeh
            ! initialize the script
            write(cmd_i,"(A)",IOSTAT=istat) 'from numpy import genfromtxt'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.plotting import &
                &figure, output_file, save, show'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.layouts import &
                &widgetbox, row'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.models.widgets import &
                &CheckboxGroup'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.models import CustomJS'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'TOOLS="hover,pan,wheel_zoom,&
                &box_zoom,reset,save"'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'output_file("'//trim(plot_dir)//&
                &'/'//trim(draw_name)//'.html", title="'//trim(draw_name)//&
                &'", mode="cdn") '
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'data = genfromtxt("'//&
                &trim(data_dir)//'/'//trim(data_name_loc)//'.dat")'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'p = figure(toolbar_location=&
                &"above",tools=TOOLS,active_drag="box_zoom",active_scroll=&
                &"wheel_zoom",plot_width=600, plot_height=600)'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            
            ! write extra options
            if (present(extra_ops)) write(cmd_i,"(A)",IOSTAT=istat) &
                &trim(extra_ops)
            
            ! individual plots
            do iplt = 1,nplt
                ! plot the lines
                write(cmd_i,"(A)",IOSTAT=istat) 'l'//trim(i2str(iplt))//&
                    &' = p.line(data[:,'//trim(i2str(iplt-1))//'],data[:,'//&
                    &trim(i2str(nplt+iplt-1))//'],'//trim(loc_draw_op(1))//')'
                ! plot the circles
                write(cmd_i,"(A)",IOSTAT=istat) 'c'//trim(i2str(iplt))//&
                    &' = p.circle(data[:,'//&
                    &trim(i2str(iplt-1))//'],data[:,'//&
                    &trim(i2str(nplt+iplt-1))//'],'//trim(loc_draw_op(2))//')'
            end do
            write(cmd_i,"(A)",IOSTAT=istat) ''
            
            ! create checkbox
            write(cmd_i,"(A)",IOSTAT=istat) 'labels=list()'
            do iplt = 1,nplt
                write(cmd_i,"(A)",IOSTAT=istat) 'labels.append("'//&
                    &trim(var_names_loc(iplt))//'")'
            end do
            write(cmd_i,"(A)",IOSTAT=istat) 'active=list()'
            do iplt = 1,nplt
                write(cmd_i,"(A)",IOSTAT=istat) 'active.append('//&
                    &trim(i2str(iplt-1))//')'
            end do
            write(cmd_i,"(A)",IOSTAT=istat) 'checkbox = CheckboxGroup(&
                &labels=labels, active=active)'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            
            ! create callback
            write(cmd_i,"(A)",IOSTAT=istat) 'args=dict()'
            do iplt = 1,nplt
                write(cmd_i,"(A)",IOSTAT=istat) 'args["l'//&
                    &trim(i2str(iplt))//'"] = l'//trim(i2str(iplt))
                write(cmd_i,"(A)",IOSTAT=istat) 'args["c'//&
                    &trim(i2str(iplt))//'"] = c'//trim(i2str(iplt))
            end do
            
            write(cmd_i,"(A)",IOSTAT=istat) 'code="""'
            write(cmd_i,"(A)",IOSTAT=istat) '//console.log(cb_obj.active);'
            do iplt = 1,nplt
                write(cmd_i,"(A)",IOSTAT=istat) 'l'//trim(i2str(iplt))//&
                    &'.visible = false;'
                write(cmd_i,"(A)",IOSTAT=istat) 'c'//trim(i2str(iplt))//&
                    &'.visible = false;'
            end do
            write(cmd_i,"(A)",IOSTAT=istat) 'for (i in cb_obj.active) {'
            write(cmd_i,"(A)",IOSTAT=istat) &
                &'    //console.log(cb_obj.active[i]);'
            do iplt = 1,nplt
                write(cmd_i,"(A)",IOSTAT=istat) '    if (cb_obj.active[i] == '&
                    &//trim(i2str(iplt-1))//') {'
                write(cmd_i,"(A)",IOSTAT=istat) '    l'//trim(i2str(iplt))//&
                    &'.visible = true;'
                write(cmd_i,"(A)",IOSTAT=istat) '    c'//trim(i2str(iplt))//&
                    &'.visible = true;'
                write(cmd_i,"(A)",IOSTAT=istat) '    }'
            end do
            write(cmd_i,"(A)",IOSTAT=istat) '}'
            write(cmd_i,"(A)",IOSTAT=istat) '"""'
            
            write(cmd_i,"(A)",IOSTAT=istat) 'checkbox.callback = CustomJS(&
                &args=args, code=code)'
            
            ! create layout
            write(cmd_i,"(A)",IOSTAT=istat) &
                &'layout = row(p, widgetbox(checkbox))'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            
            ! finishing the command
            write(cmd_i,"(A)",IOSTAT=istat) ''
            if (plot_on_screen) then
                write(cmd_i,"(A)",IOSTAT=istat) 'show(layout)'
            else
                write(cmd_i,"(A)",IOSTAT=istat) 'save(layout)'
            end if
        end subroutine draw_ex_Bokeh
        
        !> \private Mayavi version: 3D png output
        subroutine draw_ex_Mayavi
            ! initialize the script
            write(cmd_i,"(A)",IOSTAT=istat) 'from numpy import genfromtxt, &
                &array, size, zeros, amin, amax'
            write(cmd_i,"(A)",IOSTAT=istat) 'from mayavi.mlab import outline, &
                &mesh, points3d, savefig, colorbar, axes, savefig, show'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'nplt = '//trim(i2str(nplt))
            
            ! calculate the number of y points
            write(cmd_i,"(A)",IOSTAT=istat) 'npnty = -1'
            write(cmd_i,"(A)",IOSTAT=istat) 'with open("'//trim(data_dir)//&
                &'/'//trim(data_name_loc)//'.dat") as f:'
            write(cmd_i,"(A)",IOSTAT=istat) '    for line in f:'
            write(cmd_i,"(A)",IOSTAT=istat) '        if (npnty < 0):'
            write(cmd_i,"(A)",IOSTAT=istat) '            if &
                &(line.strip()[0] == "#"):'
            write(cmd_i,"(A)",IOSTAT=istat) '                npnty = 0'
            write(cmd_i,"(A)",IOSTAT=istat) '        else:'
            write(cmd_i,"(A)",IOSTAT=istat) '            if (not line.strip()):'
            write(cmd_i,"(A)",IOSTAT=istat) '                break'
            write(cmd_i,"(A)",IOSTAT=istat) '            else:'
            write(cmd_i,"(A)",IOSTAT=istat) '                npnty += 1'
            
            ! get data
            write(cmd_i,"(A)",IOSTAT=istat) 'data = genfromtxt("'//&
                &trim(data_dir)//'/'//trim(data_name_loc)//'.dat")'
            write(cmd_i,"(A)",IOSTAT=istat) 'dims = &
                &array([size(data,0)/npnty,npnty,nplt])'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'X = zeros(dims)'
            write(cmd_i,"(A)",IOSTAT=istat) 'Y = zeros(dims)'
            write(cmd_i,"(A)",IOSTAT=istat) 'Z = zeros(dims)'
            write(cmd_i,"(A)",IOSTAT=istat) 'for i in range(0,dims[0]):'
            write(cmd_i,"(A)",IOSTAT=istat) '    for j in range(0,dims[2]):'
            write(cmd_i,"(A)",IOSTAT=istat) '        X[i,:,j] = &
                &data[npnty*i:npnty*(i+1),j]'
            write(cmd_i,"(A)",IOSTAT=istat) '        Y[i,:,j] = &
                &data[npnty*i:npnty*(i+1),nplt+j]'
            write(cmd_i,"(A)",IOSTAT=istat) '        Z[i,:,j] = &
                &data[npnty*i:npnty*(i+1),nplt*2+j]'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'minZ = amin(Z)'
            write(cmd_i,"(A)",IOSTAT=istat) 'maxZ = amax(Z)'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            
            ! write extra options
            if (present(extra_ops)) write(cmd_i,"(A)",IOSTAT=istat) &
                &trim(extra_ops)
            
            ! plot
            select case (draw_dim)
                case (2)                                                        ! 3D
                    write(cmd_i,"(A)",IOSTAT=istat) 'for j in range(0,dims[2]):'
                    write(cmd_i,"(A)",IOSTAT=istat) '    &
                        &mesh(X[:,:,j],Y[:,:,j],Z[:,:,j],vmin=minZ,vmax=maxZ)'
                case (3)                                                        ! 2D slices in 3D
                    write(cmd_i,"(A)",IOSTAT=istat) 'for j in range(0,dims[2]):'
                    write(cmd_i,"(A)",IOSTAT=istat) '    &
                        &points3d(X[:,:,j],Y[:,:,j],Z[:,:,j],mode="sphere",&
                        &scale_factor=0.05,vmin=minZ,vmax=maxZ)'
            end select
            write(cmd_i,"(A)",IOSTAT=istat) 'outline()'
            write(cmd_i,"(A)",IOSTAT=istat) 'colorbar()'
            write(cmd_i,"(A)",IOSTAT=istat) 'axes()'
            
            ! pause on screen or output to file
            if (plot_on_screen) then
                write(cmd_i,"(A)",IOSTAT=istat) ''
                write(cmd_i,"(A)",IOSTAT=istat) 'show()'
            else
                write(cmd_i,"(A)",IOSTAT=istat) ''
                write(cmd_i,"(A)",IOSTAT=istat) '#show()'
                write(cmd_i,"(A)",IOSTAT=istat) 'savefig("'//trim(plot_dir)//&
                    &'/'//trim(draw_name)//'.png"'//',magnification=4)'
            end if
        end subroutine draw_ex_Mayavi
        
        !> \private GNUPlot animated version: gif output
        subroutine draw_ex_animated_GNUPlot
            ! initialize the script
            write(cmd_i,"(A)",IOSTAT=istat) 'set grid'
            write(cmd_i,"(A)",IOSTAT=istat) 'set border 4095 front linetype -1 &
                &linewidth 1.0'
            write(cmd_i,"(A)",IOSTAT=istat) 'set terminal gif animate &
                &delay '//trim(i2str(delay_loc))//' size '//&
                &trim(i2str(128*plot_size(1)))//','//&
                &trim(i2str(128*plot_size(2)))
            write(cmd_i,"(A)",IOSTAT=istat) 'set output "'//trim(plot_dir)//&
                &'/'//trim(draw_name)//'.gif"'
            
            ! find and set ranges
            allocate(ranges_loc(2,2))
            if (present(ranges)) then
                if (size(ranges,1).eq.2 .and. size(ranges,2).eq.2) then
                    ranges_loc = ranges
                else
                    call writo('invalid ranges given &
                        &to draw_ex_animated',persistent=.true.,&
                        &warning=.true.)
                end if
            else
                ! initialize ranges
                ranges_loc(:,1) = huge(1._dp)                                   ! minimum value
                ranges_loc(:,2) = -huge(1._dp)                                  ! maximum value
                
                call get_ranges(ranges_loc)
            end if
            
            ! set ranges
            write(cmd_i,"(A)",IOSTAT=istat) 'set xrange ['//&
                &trim(r2str(ranges_loc(1,1)))//':'//&
                &trim(r2str(ranges_loc(1,2)))//'];'
            write(cmd_i,"(A)",IOSTAT=istat) 'set yrange ['//&
                &trim(r2str(ranges_loc(2,1)))//':'//&
                &trim(r2str(ranges_loc(2,2)))//'];'
            
            ! write extra options
            if (present(extra_ops)) write(cmd_i,"(A)",IOSTAT=istat) &
                &trim(extra_ops)
            
            ! set up line styles
            if (n_draw_ops.eq.0) then
                do iplt = 1,min(nplt,size(line_clrs))
                    write(cmd_i,"(A)",IOSTAT=istat) 'set style line '//&
                        &trim(i2str(iplt))//' lc rgb '//&
                        &line_clrs(iplt)//' '//trim(line_style)
                end do
            end if
            
            ! individual plots
            do iplt = 1,nplt
                write(cmd_i,"(A)",IOSTAT=istat) 'plot "'//trim(data_dir)//'/'//&
                    &trim(data_name_loc)//'.dat" using '//trim(i2str(iplt))//&
                    &':'//trim(i2str(nplt+iplt))//' title "'//&
                    &trim(var_names_loc(iplt))//'" '//&
                    &trim(loc_draw_op())
            end do
            
            ! finishing the command
            write(cmd_i,"(A)",IOSTAT=istat) ''
        end subroutine draw_ex_animated_GNUPlot
        
        !> \private Bokeh animated version: html output
        subroutine draw_ex_animated_Bokeh
            ! initialize the script
            write(cmd_i,"(A)",IOSTAT=istat) '# Note: it is necessary to first &
                &run "bokeh serve" before calling python'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'from numpy import genfromtxt'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.plotting import &
                &figure, output_file, save, show, curdoc'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.client import &
                &push_session'
            write(cmd_i,"(A)",IOSTAT=istat) 'from bokeh.driving import repeat'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'TOOLS="crosshair,pan,wheel_zoom,&
                &box_zoom,reset,save"'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'output_file("'//trim(plot_dir)//&
                &'/'//trim(draw_name)//'.html", title="'//trim(draw_name)//&
                &'", mode="cdn") '
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'data = genfromtxt("'//&
                &trim(data_dir)//'/'//trim(data_name_loc)//'.dat")'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'p = figure(toolbar_location=&
                &"above",'//'tools=TOOLS,active_drag="pan",active_scroll=&
                &"wheel_zoom",'//'plot_width=600, plot_height=600)'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            
            ! write extra options
            if (present(extra_ops)) write(cmd_i,"(A)",IOSTAT=istat) &
                &trim(extra_ops)
            
            ! individual plot for first plot
            write(cmd_i,"(A)",IOSTAT=istat) 'line = p.line(data[:,0],data[:,'//&
                &trim(i2str(nplt))//'],'//trim(loc_draw_op(1))//',legend="'//&
                &trim(var_names_loc(1))//'")'
            ! plot the circles
            write(cmd_i,"(A)",IOSTAT=istat) 'circle = p.circle(data[:,'//&
                &trim(i2str(0))//'],data[:,'//trim(i2str(nplt))//'],'//&
                &trim(loc_draw_op(2))//')'
            
            ! finishing the command (needs to use show)
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'session = push_session(curdoc())'
            write(cmd_i,"(A)",IOSTAT=istat) 'session.show(p)'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'seq = range(0,'//&
                &trim(i2str(nplt))//')'
            write(cmd_i,"(A)",IOSTAT=istat) '@repeat(seq)'
            write(cmd_i,"(A)",IOSTAT=istat) 'def update(i):'
            write(cmd_i,"(A)",IOSTAT=istat) '    line.data_source.data["x"] = &
                &data[:,i]'
            write(cmd_i,"(A)",IOSTAT=istat) '    line.data_source.data["y"] = &
                &data[:,'//trim(i2str(nplt))//'+i]'
            write(cmd_i,"(A)",IOSTAT=istat) '    circle.data_source.data["x"] = &
                &data[:,i]'
            write(cmd_i,"(A)",IOSTAT=istat) '    circle.data_source.data["y"] = &
                &data[:,'//&
                &trim(i2str(nplt))//'+i]'
            write(cmd_i,"(A)",IOSTAT=istat) ''
            write(cmd_i,"(A)",IOSTAT=istat) 'curdoc().add_periodic_callback(&
                &update, '//trim(i2str(delay_loc*10))//')'
            write(cmd_i,"(A)",IOSTAT=istat) 'session.loop_until_closed()'
        end subroutine draw_ex_animated_Bokeh
        
        !> \private gets ranges for animated plot
        subroutine get_ranges(ranges)
            ! input / output
            real(dp), intent(inout) :: ranges(:,:)                              ! x and y range, and z range (if 3D) of plot
            
            ! local variables
            integer :: data_i                                                   ! file number of data file
            real(dp), allocatable :: loc_data(:)                                ! one line of data
            character(len=1) :: loc_data_char                                   ! local first char
            
            
            ! open data file
            open(nextunit(data_i),FILE=data_dir//'/'//trim(data_name_loc)//&
                &'.dat',STATUS='old',IOSTAT=istat)
            
            if (istat.eq.0) then
                ! set up loc_data
                allocate(loc_data(2*nplt))
                
                ! read the data file
                istat = 0
                do while (istat.eq.0)
                    read(data_i,*,IOSTAT=istat) loc_data_char                   ! read first character of data
                    if (istat.eq.0) then                                        ! read succesful
                        if (loc_data_char.ne.'#') then                          ! exclude comment lines
                            backspace(UNIT=data_i)                              ! go back one line
                            read(data_i,*,IOSTAT=istat) loc_data                ! read data again, but now ful
                            ranges(1,1) = &
                                &min(ranges(1,1),minval(loc_data(1:nplt)))
                            ranges(1,2) = &
                                &max(ranges(1,2),maxval(loc_data(1:nplt)))
                            ranges(2,1) = &
                                &min(ranges(2,1),&
                                &minval(loc_data(nplt+1:2*nplt)))
                            ranges(2,2) = &
                                &max(ranges(2,2),&
                                &maxval(loc_data(nplt+1:2*nplt)))
                        end if
                    end if
                end do
                !write(*,*) 'ranges(1,:) = ', ranges(1,:)
                !write(*,*) 'ranges(2,:) = ', ranges(2,:)
                !write(*,*) 'ranges(3,:) = ', ranges(3,:)
                
                ! close data file
                close(data_i,IOSTAT=istat)
            end if
        end subroutine get_ranges
        
        !> \private sets local draw options either from pre-defined line styles or through
        ! user specified option, for GNUPlot or Bokeh
        ! As for  Bokeh, there are  two default  draw options possible  (for the
        ! lines and for the points), there is an option to select between them.
        function loc_draw_op(Bokeh_style)
            ! input / output
            character(len=max_str_ln) :: loc_draw_op                            ! local drawing option
            integer, intent(in), optional :: Bokeh_style                        ! lines(1) or points (2)
            
            select case (ex_plot_style_loc)
                case (1)                                                        ! GNUPlot
                    if (n_draw_ops.gt.0) then
                        loc_draw_op = trim(draw_ops(mod(iplt-1,n_draw_ops)+1))
                    else
                        loc_draw_op = 'with linespoints linestyle '//&
                            &trim(i2str(mod(iplt-1,size(line_clrs))+1))
                    end if
                case (2)                                                        ! Bokeh / Mayavi
                    if (n_draw_ops.gt.0) then
                        loc_draw_op = trim(draw_ops(mod(iplt-1,n_draw_ops)+1))
                    else
                        select case (Bokeh_style)
                            case (1)                                            ! lines
                                loc_draw_op = 'line_color='//&
                                    &line_clrs(mod(iplt-1,size(line_clrs))+1)//&
                                    &',line_width=2'
                            case (2)                                            ! points
                                loc_draw_op = 'color='//&
                                    &line_clrs(mod(iplt-1,size(line_clrs))+1)//&
                                    &',size=5'
                            case default
                                loc_draw_op = ''
                        end select
                    end if
            end select
        end function loc_draw_op
    end subroutine draw_ex
    
    !> Takes  two input  vectors and  plots these  as well  as the  relative and
    !! absolute difference in a HDF5 file.
    !!
    !! This is similar to a basic version of output_ops.plot_hdf5().
    !! 
    !! Optionally, an output message can be displayed on screen with the maximum
    !! relative and absolute error.
    subroutine plot_diff_HDF5(A,B,file_name,tot_dim,loc_offset,descr,&
        &output_message)
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:)                                        !< vector A
        real(dp), intent(in) :: B(:,:,:)                                        !< vector B
        character(len=*), intent(in) :: file_name                               !< name of plot
        integer, intent(in), optional :: tot_dim(3)                             !< total dimensions of the arrays
        integer, intent(in), optional :: loc_offset(3)                          !< offset of local dimensions
        character(len=*), intent(in), optional :: descr                         !< description
        logical, intent(in), optional :: output_message                         !< whether to display a message or not
        
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
            call writo('in plot_diff_HDF5, A and B need to have the correct &
                &size',persistent=.true.,warning=.true.)
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
            &loc_offset=loc_offset_loc,col=1,descr=descr)
        
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
        !> \private
        function diff(A,B,dims,rel) result(C)
            use num_vars, only: tol_zero
            
            ! input / output
            real(dp), intent(in) :: A(:,:,:)                                    ! input A
            real(dp), intent(in) :: B(:,:,:)                                    ! input B
            integer, intent(in) :: dims(3)                                      ! dimensions of A and B
            logical, intent(in) :: rel                                          ! .true. if relative and .false. if absolute error
            real(dp) :: C(dims(1),dims(2),dims(3))                              ! output C
            
            ! return output
            if (rel) then
                C = 2*(A-B)/max(tol_zero,(abs(A)+abs(B)))
            else
                C = abs(A-B)
            end if
        end function diff
        
        ! returns limits and average value
        !> \private
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
    
    !> Executes command line, or displays a message if disabled.
    !!
    !! It also keeps a log of all shell commands executed.
    subroutine use_execute_command_line(command,exitstat,cmdstat,cmdmsg)
        use num_vars, only: do_execute_command_line, prog_name, &
            &shell_commands_name, shell_commands_i
#if ( lwith_intel && !lwith_gnu)
        use IFPORT
#endif
        
        ! input / output
        character(len=*), intent(in) :: command                                 !< command to execute
        integer, intent(inout), optional :: exitstat                            !< exit status
        integer, intent(inout), optional :: cmdstat                             !< command status
        character(len=*), intent(inout), optional :: cmdmsg                     !< command message
        
        ! local variables
        integer :: istat                                                        ! status
        character(len=max_str_ln) :: full_name                                  ! full name
        
        ! set up the full shell commands file name
        full_name = prog_name//'_'//trim(shell_commands_name)//'.sh'
        
        ! write the command to the file
        open(shell_commands_i,FILE=trim(full_name),STATUS='old',&
            &POSITION='append',IOSTAT=istat)
        if (istat.eq.0) then
            write(shell_commands_i,'(A)',IOSTAT=istat) trim(command)
            close(shell_commands_i,IOSTAT=istat)
        else
            call writo('Failed to write to shell commands log file "'&
                &//trim(full_name)//'"',warning=.true.)
        end if
        
        ! initialize stati
        if (present(exitstat)) exitstat = 0
        if (present(cmdstat)) cmdstat = 0
        if (present(cmdmsg)) cmdmsg = ''
        
        ! execute command line
        if (do_execute_command_line) then
#if ( lwith_intel && !lwith_gnu)
            exitstat = system(command)
#else
            call execute_command_line(command,EXITSTAT=exitstat,&
                &CMDSTAT=cmdstat,CMDMSG=cmdmsg)
#endif
        else
            call writo('Command added to log file "'//trim(full_name)//'":')
            call lvl_ud(1)
            call writo(command)
            call lvl_ud(-1)
            call writo('Run with "--do_execute_command_line" to execute &
                &command line instead')
        end if
    end subroutine
end module output_ops
