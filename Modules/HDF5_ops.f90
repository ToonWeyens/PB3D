!------------------------------------------------------------------------------!
!   Variables pertaining to HDF5 and XDMF                                      !
!   Notes about XDMF:                                                          !
!       - Collections can be spatial or temporal. If a variable is to be       !
!         evolved in time or if its domain is decomposed (e.g. the same        !
!         physical variable is defined in multiple HDF5 variables), then the   !
!         attribute name of this variable has to be the same for all the grids !
!         in the collection.                                                   !
!       - If no collection is used, but just multiple grids, the attribute     !
!         can be different but does not have to be, as the different grid are  !
!         distinguished by their different grid names, not by the attribute    !
!         of the variables they contain.                                       !
!       - The selection of hyperslabs is used, as described here:              !
!         http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/12_Dataspaces.html        !
!       - Chunking is used for partial I/O, as described here:                 !
!         https://www.hdfgroup.org/HDF5/doc/Advanced/Chunking/                 !
!         As  no  chunk  cache  is  reused,   the  w0  setting  should  be  1. !
!         Furthermore, the number of slots is chosen to be equal to the number !
!         of chunks. And Finally,  the chunk size is chosen to  be the size of !
!         the previous dimensions, if it does not exceed 4GB                   !
!------------------------------------------------------------------------------!
module HDF5_ops
#include <PB3D_macros.h>
    use num_vars, only: max_str_ln, dp, plot_dir, data_dir, script_dir
    use messages, only: writo, lvl_ud
    use str_ops, only: i2str, r2str, r2strt
    use HDF5_vars
    use HDF5
    
    implicit none
    private
    public print_HDF5_grid, print_HDF5_geom, print_HDF5_top, &
        &print_HDF5_att, print_HDF5_3D_data_item, open_HDF5_file, &
        &close_HDF5_file, reset_HDF5_item, add_HDF5_item, create_output_HDF5, &
        &print_HDF5_arrs, read_HDF5_arr
#if ldebug
    public debug_HDF5_ops, debug_print_HDF5_arrs
#endif
    
#if ldebug
    ! global variables
    logical :: debug_HDF5_ops = .false.                                         ! set to true to debug general information
    logical :: debug_print_HDF5_arrs = .false.                                  ! set to true to debug print_HDF5_arrs
#endif
    
    ! interfaces
    interface reset_HDF5_item
        module procedure reset_HDF5_item_ind, reset_HDF5_item_arr
    end interface
    interface read_HDF5_arr
        module procedure read_HDF5_arr_ind, read_HDF5_arr_range
    end interface
    
contains
    ! Opens an HDF5 file and accompanying xmf file and returns the handles.
    ! Optionally, a description of the file  can be provided. Also, the plot can
    ! be done for only one process, setting the variable "ind_plot" to .true.
    ! [MPI] Parts by all processes, parts only by group master
    integer function open_HDF5_file(file_info,file_name,description,&
        &ind_plot) result(ierr)
        use num_vars, only: rank
        use MPI
        use files_utilities, only: nextunit
        
        character(*), parameter :: rout_name = 'open_HDF5_file'
        
        ! input / output
        type(HDF5_file_type), intent(inout) :: file_info                        ! info about HDF5 file
        character(len=*), intent(in) :: file_name                               ! name of HDF5 file
        character(len=*), intent(in), optional  :: description                  ! description of file
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        character(len=max_str_ln) :: full_file_name                             ! full file name
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: plist_id                                              ! property list identifier 
        integer :: MPI_Comm                                                     ! MPI Communicator used
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! initialize ierr
        ierr = 0
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! set up MPI Communicator
        if (ind_plot_loc) then
            MPI_Comm = MPI_Comm_self                                            ! individual plot
        else
            MPI_Comm = MPI_Comm_world                                           ! default world communicator
        end if
        
        ! set up full file name
        full_file_name = data_dir//'/'//trim(file_name)
        
        ! initialize FORTRAN predefined datatypes
        call H5Open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! setup file access property list with parallel I/O access if needed
        call H5Pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
        CHCKERR('Failed to create property list')
        call H5Pset_fapl_mpio_f(plist_id,MPI_Comm,MPI_INFO_NULL,ierr)
        CHCKERR('Failed to set file access property')
        
        ! create the file collectively.
        call H5Fcreate_f(trim(full_file_name)//'.h5',H5F_ACC_TRUNC_F,HDF5_i,&
            &ierr,access_prp=plist_id)
        CHCKERR('Failed to create file')
        call H5Pclose_f(plist_id,ierr)
        CHCKERR('Failed to close property list')
        
        ! update file info, converting integers and setting name
        file_info%HDF5_i = HDF5_i
        file_info%name = file_name
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! open accompanying xmf file
            open(unit=nextunit(file_info%XDMF_i),&
                &file=trim(full_file_name)//'.xmf',iostat=ierr)
            CHCKERR('Failed to open xmf file')
            
            ! write header if group master
            write(file_info%XDMF_i,xmf_fmt) '<?xml version="1.0" ?>'
            write(file_info%XDMF_i,xmf_fmt) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" &
                &[]>'
            write(file_info%XDMF_i,xmf_fmt) '<Xdmf Version="2.0">'
            write(file_info%XDMF_i,xmf_fmt) '<Domain>'
            if (present(description)) then
                write(file_info%XDMF_i,xmf_fmt) &
                    &'<Information Name="Description">'
                write(file_info%XDMF_i,xmf_fmt) trim(description)
                write(file_info%XDMF_i,xmf_fmt) '</Information>'
            end if
        end if
    end function open_HDF5_file
    
    ! Closes an HDF5 file and writes the accompanying xmf file
    ! [MPI] Parts by all processes, parts only by group master
    integer function close_HDF5_file(file_info,ind_plot) result(ierr)
        use num_vars, only: rank
        
        character(*), parameter :: rout_name = 'close_HDF5_file'
        
        ! input / output
        type(HDF5_file_type), intent(inout) :: file_info                        ! info about HDF5 file
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: full_file_name                             ! full file name
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! initialize ierr
        ierr = 0
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! set full file name and HDF5_i (converting integers)
        full_file_name = data_dir//'/'//trim(file_info%name)
        HDF5_i = file_info%HDF5_i
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! close header
            write(file_info%XDMF_i,xmf_fmt) '</Domain>'
            write(file_info%XDMF_i,xmf_fmt) '</Xdmf>'
            
            ! close accompanying xmf file
            close(file_info%XDMF_i,IOSTAT=ierr)
            CHCKERR('Failed to close xmf file')
            
            ! user output
            call writo('Created HDF5/XMF plot in output file "'//&
                &trim(full_file_name)//'.xmf''')
        end if
    end function close_HDF5_file
    
    ! Add an XDMF item to a HDF5 file
    ! Note:  This  should  only  be  used  with  grids  (or for  topologies  and
    ! geometries that are used throughout)
    ! [MPI] Only group master
    subroutine add_HDF5_item(file_info,XDMF_item,reset,ind_plot)
        use num_vars, only: rank
        
        ! input / output
        type(HDF5_file_type), intent(inout) :: file_info                        ! info about HDF5 file
        type(XML_str_type), intent(inout) :: XDMF_item                          ! XDMF item to add
        logical, intent(in), optional :: reset                                  ! if .true., XDMF_item is reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: item_len                                                     ! length of item
        logical :: reset_loc                                                    ! local copy of reset
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! set item_len
            item_len = size(XDMF_item%xml_str)
            
            ! set local reset
            if (present(reset)) then
                reset_loc = reset
            else
                reset_loc = .false.
            end if
            
            ! write the strings to the file
            do id = 1,item_len
                write(file_info%XDMF_i,xmf_fmt) trim(XDMF_item%xml_str(id))
            end do
            
            ! reset if requested
            if (reset_loc) call reset_HDF5_item(XDMF_item,ind_plot_loc)
        end if
    end subroutine add_HDF5_item
    
    ! prints an HDF5 data item
    ! If this is a parallel data item, the group dimension and offset have to be
    ! specified as well.
    ! [MPI] Parts by all processes, parts only by group master
    integer function print_HDF5_3D_data_item(dataitem_id,file_info,var_name,&
        &var,dim_tot,loc_dim,loc_offset,init_val,ind_plot) result(ierr)
        use num_vars, only: rank
        
        character(*), parameter :: rout_name = 'print_HDF5_3D_data_item'
        
        ! input / output
        type(XML_str_type), intent(inout) :: dataitem_id                        ! ID of data item
        type(HDF5_file_type), intent(in) :: file_info                           ! info about HDF5 file
        character(len=*), intent(in) :: var_name                                ! name of variable
        real(dp), intent(in) :: var(:,:,:)                                      ! variable to write
        integer, intent(in) :: dim_tot(3)                                       ! total dimensions of variable
        integer, intent(in), optional :: loc_dim(3)                             ! dimensions in this group
        integer, intent(in), optional :: loc_offset(3)                          ! offset in this group
        real(dp), intent(in), optional :: init_val                              ! initial fill value
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: loc_dim_loc(3)                                               ! local copy of loc_dim
        integer :: loc_offset_loc(3)                                            ! local copy of loc_offset
        integer(HSIZE_T) :: dimsf(3)                                            ! total dataset dimensions
        integer(HSIZE_T) :: dimsm(3)                                            ! group dataset dimensions
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HID_T) :: memspace                                              ! Dataspace identifier in memory
        integer(HID_T) :: plist_id                                              ! property list identifier 
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        integer(HSIZE_T) :: mem_offset(3)                                       ! offset in memory
        integer(HSIZE_T) :: mem_block(3)                                        ! block size in memory
        integer(HSIZE_T) :: mem_stride(3)                                       ! stride in memory
        integer(HSIZE_T) :: mem_count(3)                                        ! nr. of repetitions of block in memory
        character(len=max_str_ln) :: dim_str                                    ! string with dimensions
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        integer(HID_T) :: HDF5_kind_64                                          ! HDF5 type corresponding to dp
        
        ! initialize ierr
        ierr = 0
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! set HDF5 type corresponding to dp
        HDF5_kind_64 = H5Kind_to_type(dp,H5_REAL_KIND)
        
        ! set the group dimensions and offset
        ierr = check_for_parallel_3D(dim_tot,loc_dim_loc,loc_offset_loc,&
            &loc_dim,loc_offset)
        CHCKERR('')
        
        ! create the data spaces for the dataset. 
        dimsf = dim_tot
        dimsm = loc_dim_loc
        call H5Screate_simple_f(size(dimsf),dimsf,filespace,ierr)
        CHCKERR('Failed to create file space')
        call H5Screate_simple_f(size(dimsm),dimsm,memspace,ierr)
        CHCKERR('Failed to create memory space')
        
        ! set initial fill value if provided and create data set
        if (present(init_val)) then
            call H5Pcreate_f(H5P_DATASET_CREATE_F,plist_id,ierr)
            err_msg = 'Failed to create property list'
            CHCKERR(err_msg)
            call H5Pset_fill_value_f(plist_id,HDF5_kind_64,init_val,ierr) 
            err_msg = 'Failed to set default fill value property'
            CHCKERR(err_msg)
            call H5Dcreate_f(file_info%HDF5_i,trim(var_name),HDF5_kind_64,&
                &filespace,dset_id,ierr,plist_id)
            CHCKERR('Failed to create file data set')
            call H5Pclose_f(plist_id,ierr)
            CHCKERR('Failed to close property list')
        else
            call H5Dcreate_f(file_info%HDF5_i,trim(var_name),HDF5_kind_64,&
                &filespace,dset_id,ierr)
            CHCKERR('Failed to create file data set')
        end if
        call H5Sclose_f(filespace,ierr)
        CHCKERR('Failed to close file space')
        
        ! select hyperslab in the file.
        mem_count = 1
        mem_stride = 1
        mem_block = loc_dim_loc
        mem_offset = loc_offset_loc
        call H5Dget_space_f(dset_id,filespace,ierr)
        CHCKERR('Failed to get file space')
        call H5Sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,mem_offset,&
            &mem_count,ierr,mem_stride,mem_block)
        CHCKERR('Failed to select hyperslab')
        
        ! create property list for collective dataset write
        call H5Pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr) 
        CHCKERR('Failed to create property list')
        call H5Pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
        CHCKERR('Failed to set parallel property')
        
        ! write the dataset collectively. 
        call H5Dwrite_f(dset_id,HDF5_kind_64,var,dimsf,ierr,&
            &file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
        CHCKERR('Failed to write data set')
        call H5Pclose_f(plist_id,ierr)
        CHCKERR('Failed to close property list')
        
        ! close dataspaces.
        call H5Sclose_f(filespace,ierr)
        CHCKERR('Unable to close file space')
        call H5Sclose_f(memspace,ierr)
        CHCKERR('Unable to close memory space')
        
        ! close the dataset.
        call H5Dclose_f(dset_id,ierr)
        CHCKERR('Failed to close data set')
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! set up dimensions string
            dim_str = ''
            do id = 1,size(dim_tot)
                if (dim_tot(id).gt.1) then
                    dim_str = trim(dim_str)//' '//trim(i2str(dim_tot(id)))
                end if
            end do
            dataitem_id%name = 'DataItem - '//trim(var_name)
            allocate(dataitem_id%xml_str(3))
            dataitem_id%xml_str(1) = '<DataItem Dimensions="'//trim(dim_str)//&
                &'" NumberType="Float" Precision="8" Format="HDF">'
            dataitem_id%xml_str(2) = trim(file_info%name)//'.h5:/'//&
                &trim(var_name)
            dataitem_id%xml_str(3) = '</DataItem>'
            
            if (debug_HDF5_ops) write(*,*) 'created data item "'//&
                &trim(dataitem_id%name)//'"'
        end if
    contains
        ! check whether the  dimensions provided are a  valid parallel indicator
        ! or not
        integer function check_for_parallel_3D(dim_tot,loc_dim_out,&
            &loc_offset_out,loc_dim_in,loc_offset_in) result(ierr)
            character(*), parameter :: rout_name = 'check_for_parallel_3D'
            
            ! input / output
            integer, intent(in) :: dim_tot(3)                                   ! total dimension
            integer, intent(inout) :: loc_dim_out(3), loc_offset_out(3)         ! output group dimension and offset
            integer, intent(in), optional :: loc_dim_in(3), loc_offset_in(3)    ! input group dimension and offset
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            integer :: id                                                       ! counter
            
            ! initialize ierr
            ierr = 0
            
            if (present(loc_dim_in)) then
                if (.not.present(loc_offset_in)) then
                    err_msg = 'Need to specify offset as well as group &
                        &dimensions'
                    ierr = 1
                    CHCKERR(err_msg)
                else
                    loc_dim_out = loc_dim_in
                    loc_offset_out = loc_offset_in
                end if
                do id = 1,size(dim_tot)
                    if (loc_dim_in(id).gt.dim_tot(id)) then                     ! error in parallel file
                        err_msg = 'Total dimension '//trim(i2str(id))//&
                            &' cannot be smaller than group dimension'
                        ierr = 1
                        CHCKERR(err_msg)
                    end if
                end do
            else
                loc_dim_out = dim_tot
                loc_offset_out = 0
            end if
        end function check_for_parallel_3D
    end function print_HDF5_3D_data_item
    
    ! prints an HDF5 attribute
    ! [MPI] Only group master
    subroutine print_HDF5_att(att_id,att_dataitem,att_name,att_center,reset,&
        &ind_plot)
        use num_vars, only: rank
        
        ! input / output
        type(XML_str_type), intent(inout) :: att_id                             ! ID of attribute
        type(XML_str_type), intent(inout) :: att_dataitem                       ! dataitem of attribute
        character(len=*), intent(in) :: att_name                                ! name of attribute
        integer, intent(in) :: att_center                                       ! center of attribute
        logical, intent(in), optional :: reset                                  ! if .true., data items are reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: dataitem_len                                                 ! length of data item
        integer :: id                                                           ! counter
        logical :: reset_loc                                                    ! local copy of reset
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! set dataitem_len
            dataitem_len = size(att_dataitem%xml_str)
            
            ! set local reset
            if (present(reset)) then
                reset_loc = reset
            else
                reset_loc = .false.
            end if
            
            ! set XDMF attribute ID
            att_id%name = 'Attribute - '//trim(att_name)
            allocate(att_id%xml_str(dataitem_len+2))
            att_id%xml_str(1) = '<Attribute Name="'//trim(att_name)//&
                &'" AttributeType="Scalar" Center="'//&
                &trim(XDMF_center_types(att_center))//'">'
            do id = 1,dataitem_len
                att_id%xml_str(id+1) = att_dataitem%xml_str(id)
            end do
            att_id%xml_str(dataitem_len+2) = '</Attribute>'
            
            if (debug_HDF5_ops) write(*,*) 'created attribute "'//&
                &trim(att_id%name)//'"'
            
            ! reset if requested
            if (reset_loc) call reset_HDF5_item(att_dataitem,ind_plot_loc)
        end if
    end subroutine print_HDF5_att
    
    ! prints an HDF5 topology
    ! Note: currently only structured grids supported
    ! [MPI] Only group master
    subroutine print_HDF5_top(top_id,top_type,top_n_elem,ind_plot)
        use num_vars, only: rank
        
        ! input / output
        type(XML_str_type), intent(inout) ::  top_id                            ! ID of topology
        integer, intent(in) :: top_type                                         ! type
        integer, intent(in) :: top_n_elem(:)                                    ! nr. of elements
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_dims                                                       ! nr. of dimensions
        character(len=max_str_ln) :: work_str                                   ! work string
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! set n_dims
            n_dims = size(top_n_elem)
            
            ! fill work string
            work_str = ''
            do id = n_dims,1,-1
                if (top_n_elem(id).gt.1) work_str = &
                    &trim(work_str)//' '//trim(i2str(top_n_elem(id)))
            end do
            
            ! set XDMF topology ID
            top_id%name = 'Topology'
            allocate(top_id%xml_str(1))
            top_id%xml_str(1) = '<Topology TopologyType="'//&
                &trim(XDMF_top_types(top_type))//'" NumberOfElements="'//&
                &trim(work_str)//'"/>'
            
            if (debug_HDF5_ops) write(*,*) 'created topology "'//&
                &trim(top_id%name)//'"'
        end if
    end subroutine print_HDF5_top
    
    ! prints an HDF5 geometry
    ! [MPI] Only group master
    subroutine print_HDF5_geom(geom_id,geom_type,geom_dataitems,reset,ind_plot)
        use num_vars, only: rank
        
        ! input / output
        type(XML_str_type), intent(inout) ::  geom_id                           ! ID of geometry
        integer, intent(in) :: geom_type                                        ! type of geometry
        type(XML_str_type), intent(inout) :: geom_dataitems(:)                  ! data items of geometry
        logical, intent(in), optional :: reset                                  ! if .true., data items are reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: id_sum                                                       ! counter
        integer, allocatable :: dataitem_len(:)                                 ! length of data item
        integer :: n_dataitems                                                  ! nr. of data items
        logical :: reset_loc                                                    ! local copy of reset
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! set n_dataitems
            n_dataitems = size(geom_dataitems)
            
            ! set dataitem_len for every data item
            allocate(dataitem_len(n_dataitems))
            do id = 1,n_dataitems
                dataitem_len(id) = size(geom_dataitems(id)%xml_str)
            end do
            
            ! set local reset
            if (present(reset)) then
                reset_loc = reset
            else
                reset_loc = .false.
            end if
            
            ! set XDMF geometry ID
            geom_id%name = 'Geometry'
            allocate(geom_id%xml_str(sum(dataitem_len)+2))
            geom_id%xml_str(1) = '<Geometry GeometryType="'//&
                &trim(XDMF_geom_types(geom_type))//'">'
            id_sum = 2                                                          ! starting index
            do id = 1,n_dataitems
                do jd = 1,dataitem_len(id)
                    geom_id%xml_str(id_sum) = geom_dataitems(id)%xml_str(jd)
                    id_sum = id_sum + 1
                end do
            end do
            geom_id%xml_str(id_sum) = '</Geometry>'
            
            if (debug_HDF5_ops) write(*,*) 'created geometry "'//&
                &trim(geom_id%name)//'"'
            
            ! reset if requested
            if (reset_loc) call reset_HDF5_item(geom_dataitems,ind_plot_loc)
        end if
    end subroutine print_HDF5_geom
    
    ! prints an HDF5 grid
    ! If  this is  is  a uniform  grid,  the  geometry and  topology  has to  be
    ! specified, or else  it will be assumed  that it is already  present in the
    ! XDMF domain, and  reused. If the grid  is a collection grid,  the grids in
    ! the  collection have  to  be specified.  Optionally, also  a  time can  be
    ! specified (for the grids in a collection grid).
    ! [MPI] Only group master
    integer function print_HDF5_grid(grid_id,grid_name,grid_type,grid_time,&
        &grid_top,grid_geom,grid_atts,grid_grids,reset,ind_plot) result(ierr)
        use num_vars, only: rank
        
        character(*), parameter :: rout_name = 'print_HDF5_grid'
        
        ! input / output
        type(XML_str_type), intent(inout) :: grid_id                            ! ID of grid
        character(len=*), intent(in) :: grid_name                               ! name
        integer, intent(in) :: grid_type                                        ! type
        real(dp), intent(in), optional :: grid_time                             ! time of grid
        type(XML_str_type), optional :: grid_top                                ! topology
        type(XML_str_type), optional :: grid_geom                               ! geometry
        type(XML_str_type), optional :: grid_atts(:)                            ! attributes
        type(XML_str_type), optional :: grid_grids(:)                           ! grids
        logical, intent(in), optional :: reset                                  ! if .true., top, geom and atts or grids are reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: id_sum                                                       ! counter
        integer :: n_grids                                                      ! nr. of grids
        integer :: n_atts                                                       ! nr. of attributes
        integer :: time_len                                                     ! length of time
        integer :: top_len                                                      ! length of topology
        integer :: geom_len                                                     ! length of geometry
        integer, allocatable :: atts_len(:)                                     ! lengths of attributes
        integer, allocatable :: grids_len(:)                                    ! lengths of grids
        logical :: reset_loc                                                    ! local copy of reset
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! initialize ierr
        ierr = 0
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            ! test whether the correct arguments are provided
            if (grid_type.eq.1) then                                            ! uniform grid
                ! no  requirements: topology and/or  geometry can be  defined in
                ! domain for use throughout
            else if (grid_type.eq.2 .or. grid_type.eq.3) then                   ! collection grid
                if (.not.present(grid_grids)) then
                    ierr = 1
                    err_msg = 'For grid collections, the grids in the &
                        &collection have to be specified'
                    CHCKERR(err_msg)
                end if
            else
                ierr = 1
                err_msg = 'Grid type '//trim(i2str(grid_type))//' not supported'
                CHCKERR(err_msg)
            end if
            
            ! set geom_len, atts_len and grids_len
            time_len = 0
            if (present(grid_time)) time_len = 1
            top_len = 0
            if (present(grid_top)) then
                top_len = size(grid_top%xml_str)
            else
                if (grid_type.eq.1) top_len = 1                                 ! use of general domain topology
            end if
            geom_len = 0
            if (present(grid_geom)) then
                geom_len = size(grid_geom%xml_str)
            else
                if (grid_type.eq.1) geom_len = 1                                ! use of general domain geometry
            end if
            n_atts = 0
            if (present(grid_atts)) then
                n_atts = size(grid_atts)
                allocate(atts_len(n_atts))
                atts_len = 0
                do id = 1,n_atts
                    atts_len(id) = size(grid_atts(id)%xml_str)
                end do
            else
                allocate(atts_len(0))
                atts_len = 0
            end if
            n_grids = 0
            if (present(grid_grids)) then
                n_grids = size(grid_grids)
                allocate(grids_len(n_grids))
                grids_len = 0
                do id = 1,n_grids
                    grids_len(id) = size(grid_grids(id)%xml_str)
                end do
            else
                allocate(grids_len(0))
                grids_len = 0
            end if
            
            ! set local reset
            if (present(reset)) then
                reset_loc = reset
            else
                reset_loc = .false.
            end if
            
            ! set XDMF grid ID
            id_sum = 2                                                          ! starting index
            if (grid_type.eq.1) then                                            ! uniform
                grid_id%name = 'Uniform Grid - '//trim(grid_name)
                allocate(&
                    &grid_id%xml_str(time_len+top_len+geom_len+sum(atts_len)+2))
                grid_id%xml_str(1) = '<Grid Name="'//trim(grid_name)//&
                    &'" GridType="Uniform">'
                if (present(grid_time)) then                                    ! time
                    grid_id%xml_str(id_sum) = '<Time Value="'//&
                        &trim(r2str(grid_time))//'" />'
                    id_sum = id_sum + 1
                end if
                if (present(grid_top)) then
                    do jd = 1,top_len                                           ! topology
                        grid_id%xml_str(id_sum) = grid_top%xml_str(jd)
                        id_sum = id_sum + 1
                    end do
                else
                    grid_id%xml_str(id_sum) = '<Topology Reference="/Xdmf/&
                        &Domain/Topology[1]"/>'
                    id_sum = id_sum + 1
                end if
                if (present(grid_geom)) then
                    do jd = 1,geom_len                                          ! geometry
                        grid_id%xml_str(id_sum) = grid_geom%xml_str(jd)
                        id_sum = id_sum + 1
                    end do
                else
                    grid_id%xml_str(id_sum) = '<Geometry Reference="/Xdmf/&
                        &Domain/Geometry[1]"/>'
                    id_sum = id_sum + 1
                end if
                do id = 1,n_atts                                                ! attributes
                    do jd = 1,atts_len(id)
                        grid_id%xml_str(id_sum) = grid_atts(id)%xml_str(jd)
                        id_sum = id_sum + 1
                    end do
                end do
            else                                                                ! collection
                grid_id%name = 'Collection Grid - '//trim(grid_name)
                allocate(grid_id%xml_str(time_len+sum(grids_len)+2))
                grid_id%xml_str(1) = '<Grid Name="'//trim(grid_name)//&
                    &'" GridType="Collection" CollectionType="'//&
                    &trim(XDMF_grid_types(grid_type))//'">'
                if (present(grid_time)) then                                    ! time
                    grid_id%xml_str(id_sum) = '<Time Value="'//&
                        &trim(r2str(grid_time))//'" />'
                    id_sum = id_sum + 1
                end if
                do id = 1,n_grids                                               ! grids
                    do jd = 1,grids_len(id)
                        grid_id%xml_str(id_sum) = grid_grids(id)%xml_str(jd)
                        id_sum = id_sum + 1
                    end do
                end do
            end if
            grid_id%xml_str(id_sum) = '</Grid>'
            
            if (debug_HDF5_ops) write(*,*) 'created grid "'//&
                &trim(grid_id%name)//'"'
            
            ! reset if requested
            if (reset_loc) then
                if (grid_type.eq.1)  then
                    if (present(grid_top)) &
                        &call reset_HDF5_item(grid_top,ind_plot_loc)
                    if (present(grid_geom)) &
                        &call reset_HDF5_item(grid_geom,ind_plot_loc)
                    if (present(grid_atts)) &
                        &call reset_HDF5_item(grid_atts,ind_plot_loc)
                else if (grid_type.eq.2 .or. grid_type.eq.3) then
                    if (present(grid_grids)) &
                        &call reset_HDF5_item(grid_grids,ind_plot_loc)
                end if
            end if
        end if
    end function print_HDF5_grid
    
    ! creates an HDF5 output file
    integer function create_output_HDF5(HDF5_name) result(ierr)
        character(*), parameter :: rout_name = 'create_output_HDF5'
        
        ! input / output
        character(len=*), intent(in) :: HDF5_name                               ! name of HDF5 file
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! initialize FORTRAN predefined datatypes
        call H5Open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! create the file by group master
        call H5Fcreate_f(trim(HDF5_name),H5F_ACC_TRUNC_F,HDF5_i,&
            &ierr)
        CHCKERR('Failed to create file')
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
        
        call lvl_ud(-1)
        
        ! user output
        call writo('HDF5 output file '//trim(HDF5_name)//' created')
    end function create_output_HDF5
    
    ! Prints a series of arrays, in the form of an array of pointers, to an HDF5
    ! file.
    ! Optionally, output can be given about the variable being written.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the head
    ! name if  it is  > 0.  Similarly, if "eq_job"  is provided,  "_E_eq_job" is
    ! appended.
    ! Note: See https://www.hdfgroup.org/HDF5/doc/UG/12_Dataspaces.html, 7.4.2.3
    ! for an explanation of the selection of the dataspaces.
    integer function print_HDF5_arrs(vars,PB3D_name,head_name,rich_lvl,eq_job,&
        &disp_info,ind_print) result(ierr)
        use num_vars, only: n_procs, rank, HDF5_lock_file_name
        use messages, only: lvl_ud
        use files_utilities, only: wait_file
        use MPI
        use HDF5_utilities, only: set_1D_vars
        
        character(*), parameter :: rout_name = 'print_HDF5_arrs'
        
        ! input / output
        type(var_1D_type), intent(in) :: vars(:)                                ! variables to write
        character(len=*), intent(in) :: PB3D_name                               ! name of PB3D file
        character(len=*), intent(in) :: head_name                               ! head name of variables
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to append to file name
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to append to file name
        logical, intent(in), optional :: disp_info                              ! display additional information about variable being read
        logical, intent(in), optional :: ind_print                              ! individual write (possibly partial I/O)
        
        ! local variables
        integer :: id, jd                                                       ! counter
        integer :: MPI_Comm                                                     ! MPI Communicator used
        integer :: istat                                                        ! status
        integer :: n_dims                                                       ! nr. of dimensions
        integer :: lock_file_i                                                  ! lock file number
        integer, allocatable :: lim_tot(:,:)                                    ! total limits
        integer, allocatable :: lim_loc(:,:)                                    ! local limits
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: head_name_loc                              ! local head_name
        logical :: ind_print_loc                                                ! local ind_print
        logical :: disp_info_loc                                                ! local disp_info
        integer(HID_T) :: a_plist_id                                            ! access property list identifier 
        integer(HID_T) :: chunk_a_plist_id                                      ! chunk access property list identifier 
        integer(HID_T) :: chunk_c_plist_id                                      ! chunk create property list identifier 
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: HDF5_kind_64                                          ! HDF5 type corresponding to dp
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HID_T) :: memspace                                              ! Dataspace identifier in memory
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        integer(HID_T) :: group_id                                              ! group identifier
        integer(HID_T) :: head_group_id                                         ! head group identifier
        integer(HSIZE_T) :: dimsf(1)                                            ! total dataset dimensions
        integer(HSIZE_T) :: dimsm(1)                                            ! group dataset dimensions
#if ldebug
        real :: rdcc_w0                                                         ! Preemption policy.
        integer(SIZE_T) :: rdcc_nslots                                          ! Number of chunk slots in the raw data chunk  cache hash table.
        integer(SIZE_T) :: rdcc_nbytes                                          ! Total size of the raw data chunk cache, in bytes. 
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set local disp_info
        disp_info_loc = .false.
        if (present(disp_info)) disp_info_loc = disp_info
        
        ! set up local head_name
        head_name_loc = head_name
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) head_name_loc = trim(head_name_loc)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        if (present(eq_job)) then
            if (eq_job.gt.0) head_name_loc = trim(head_name_loc)//'_E_'//&
                &trim(i2str(eq_job))
        end if
        
        ! detect whether individual or collective print
        call detect_ind_print(ind_print_loc)
        if (present(ind_print)) ind_print_loc = ind_print
        
        ! set up MPI Communicator
        if (ind_print_loc) then
            MPI_Comm = MPI_Comm_self                                            ! individual print
        else
            MPI_Comm = MPI_Comm_world                                           ! default world communicator
        end if
        
        ! initialize FORTRAN predefined datatypes
        call H5Open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! setup file access property list with parallel I/O access if needed
        call H5Pcreate_f(H5P_FILE_ACCESS_F,a_plist_id,ierr)
        CHCKERR('Failed to create property list')
        call H5Pset_fapl_mpio_f(a_plist_id,MPI_Comm,MPI_INFO_NULL,ierr)
        CHCKERR('Failed to set file access property')
        
        ! wait for file access if individual print and multiple processes
        if (n_procs.gt.1 .and. ind_print_loc) then
            call wait_file(lock_file_i,HDF5_lock_file_name)
        end if
        
        ! open the file
        call H5Fopen_f(trim(PB3D_name),H5F_ACC_RDWR_F,HDF5_i,ierr,&
            &access_prp=a_plist_id)
        CHCKERR('Failed to open file')
        call H5Pclose_f(a_plist_id,ierr)
        CHCKERR('Failed to close property list')
        
        ! set HDF5 type corresponding to dp
        HDF5_kind_64 = H5Kind_to_type(dp,H5_REAL_KIND)
        
        ! check if head group exists and if not create it
        call H5Eset_auto_f(0,ierr)
        CHCKERR('Failed to disable error printing')
        call H5Gopen_f(HDF5_i,trim(head_name_loc),head_group_id,istat)
        if (istat.ne.0) then                                                    ! group does not exist
            call H5Gcreate_f(HDF5_i,trim(head_name_loc),head_group_id,ierr)
            CHCKERR('Failed to create group')
        end if
        call H5Eset_auto_f(1,ierr)
        CHCKERR('Failed to enable error printing')
        
        ! user output
        call writo('Write data to PB3D output "'//trim(PB3D_name)//'/'//&
            &trim(head_name_loc)//'"')
        call lvl_ud(1)
        
        ! loop over all elements in vars
        do id = 1,size(vars)
            ! user output
            if (disp_info_loc) then
                call writo('Writing '//trim(vars(id)%var_name))
            end if
            
            ! check if group exists and if not create it
            call H5Eset_auto_f(0,ierr)
            CHCKERR('Failed to disable error printing')
            call H5Gopen_f(head_group_id,trim(vars(id)%var_name),group_id,&
                &istat)
            call H5Eset_auto_f(1,ierr)
            CHCKERR('Failed to enable error printing')
            if (istat.ne.0) then                                                ! group does not exist
                ! create group
                call H5Gcreate_f(head_group_id,trim(vars(id)%var_name),&
                    &group_id,ierr)
                CHCKERR('Failed to create group')
            end if
            
            ! 1. Data itself
            
            ! set up local and total limits
            n_dims = size(vars(id)%tot_i_min)
            lim_tot = &
                &reshape([vars(id)%tot_i_min,vars(id)%tot_i_max],[n_dims,2])
            lim_loc = &
                &reshape([vars(id)%loc_i_min,vars(id)%loc_i_max],[n_dims,2])
            
            ! set dimsf
            dimsf = product(lim_tot(:,2)-lim_tot(:,1)+1)                        ! file space has total dimensions
            
            ! create or open data set with filespace
            if (istat.ne.0) then                                                ! group did not exist
                ! create file data space
                call H5Screate_simple_f(1,dimsf,filespace,ierr)
                CHCKERR('Failed to create file space')
                
                ! set up chunk creation and access property list
                ierr = set_1D_vars(lim_tot,lim_loc,c_plist_id=chunk_c_plist_id,&
                    &a_plist_id=chunk_a_plist_id)
                CHCKERR('')

                ! create file data set in group
                call H5Dcreate_f(group_id,'var',HDF5_kind_64,filespace,dset_id,&
                    &ierr,dcpl_id=chunk_c_plist_id,dapl_id=chunk_a_plist_id)
                CHCKERR('Failed to create file data set')
                
                ! close chunk property list
                call H5Pclose_f(chunk_c_plist_id,ierr)
                CHCKERR('Failed to close property list')
            else                                                                ! group already exists
                ! set up chunk access property list
                ierr = set_1D_vars(lim_tot,lim_loc,a_plist_id=chunk_a_plist_id)
                CHCKERR('')

                ! open file data set
                call H5Dopen_f(group_id,'var',dset_id,ierr,dapl_id=&
                    &chunk_a_plist_id)
                CHCKERR('Failed to open file data set')
                
                ! get file data space
                call H5Dget_space_f(dset_id,filespace,ierr) 
                CHCKERR('Failed to get file space')
            end if
            
            ! close data access property list
            call H5Pclose_f(chunk_a_plist_id,ierr)
            CHCKERR('Failed to close property list')
            
#if ldebug
            if (debug_print_HDF5_arrs) then
                ! get data access property list
                call H5Dget_access_plist_f(dset_id,a_plist_id,ierr)
                CHCKERR('Failed to get property list')
                
                ! get chunk cache
                call H5PGet_chunk_cache_f(a_plist_id,rdcc_nslots,&
                    &rdcc_nbytes,rdcc_w0,ierr)
                write(*,*) 'Number of chunk slots in the raw data chunk cache &
                    &hash table:', rdcc_nslots
                write(*,*) 'Total size of the raw data chunk cache, in &
                    &Mbytes:', rdcc_nbytes*1.E-6_dp
                write(*,*) 'Preemption Policy:', rdcc_w0
                CHCKERR('Failed to get chunk cache')
                
                ! close data access property list
                call H5Pclose_f(a_plist_id,ierr)
                CHCKERR('Failed to close property list')
            end if
#endif
            
            ! set 1D filespace hyperslab selection
            ierr = set_1D_vars(lim_tot,lim_loc,space_id=filespace)
            CHCKERR('')
            
            ! create property list for collective dataset write
            call H5Pcreate_f(H5P_DATASET_XFER_F,a_plist_id,ierr) 
            CHCKERR('Failed to create property list')
            call H5Pset_dxpl_mpio_f(a_plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
            CHCKERR('Failed to set parallel property')
            
            ! set dimsm
            dimsm = product(lim_loc(:,2)-lim_loc(:,1)+1)                        ! memory space has local dimensions
            
            ! create memory data space
            call H5Screate_simple_f(1,dimsm,memspace,ierr)
            CHCKERR('Failed to create memory space')
            
            ! write the dataset collectively. 
            call H5Dwrite_f(dset_id,HDF5_kind_64,vars(id)%p,&
                &dimsf,ierr,file_space_id=filespace,&
                &mem_space_id=memspace,xfer_prp=a_plist_id)
            if (ierr.ne.0) call writo('Did you increase max_tot_mem_per_proc &
                &while restarting Richardson? If so, must start from 1...',&
                &alert=.true.)
            CHCKERR('Failed to write data data set')
            call H5Pclose_f(a_plist_id,ierr)
            CHCKERR('Failed to close property list')
            
            ! close dataspaces
            call H5Sclose_f(filespace,ierr)
            CHCKERR('Unable to close file space')
            call H5Sclose_f(memspace,ierr)
            CHCKERR('Unable to close memory space')
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! 2. Limit information
            
            ! set dimsf
            dimsf = size(lim_tot)                                               ! min. and max. limit value per dimension
            
            ! create or open data set with filespace
            if (istat.ne.0) then                                                ! group did not exist
                ! create file data space
                call H5Screate_simple_f(1,dimsf,filespace,ierr)
                CHCKERR('Failed to create file space')
                
                ! create file data set in group
                call H5Dcreate_f(group_id,'lim',H5T_NATIVE_INTEGER,filespace,&
                    &dset_id,ierr)
                CHCKERR('Failed to create file data set')
                
                ! close dataspace
                call H5Sclose_f(filespace,ierr)
                CHCKERR('Unable to close file space')
            else                                                                ! group already exists
                ! open file data set
                call H5Dopen_f(group_id,'lim',dset_id,ierr)
                CHCKERR('Failed to open file data set')
            end if
            
            ! only leader writes
            if (rank.eq.0 .or. ind_print_loc) then
                ! write the dataset
                call H5Dwrite_f(dset_id,H5T_NATIVE_INTEGER,&
                    &[vars(id)%tot_i_min,vars(id)%tot_i_max],dimsf,ierr)
                CHCKERR('Failed to write limit data set')
            end if
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! close the group
            call H5gclose_f(group_id,ierr)
            CHCKERR('Failed to close group')
            
            ! deallocate limits
            deallocate(lim_tot,lim_loc)
        end do
        
        call lvl_ud(-1)
        
        ! close the head group
        call H5gclose_f(head_group_id,ierr)
        CHCKERR('Failed to close group')
        
        ! close the HDF5 and possibly lock file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        if (n_procs.gt.1 .and. ind_print_loc) then
            close(lock_file_i,status='DELETE',iostat=ierr)
            CHCKERR('Failed to delete lock file')
        end if
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
    contains
        ! detects whether individual print
        subroutine detect_ind_print(ind_print)
            ! input / output
            logical, intent(inout) :: ind_print                                 ! whether an individual or collective print
            
            ! initialie ind_print to true
            ind_print = .true.
            
            ! one process ensures individual print
            if (n_procs.eq.1) return
            
            ! loop over all elements in vars
            do id = 1,size(vars)
                do jd = 1,size(vars(id)%tot_i_min)
                    if (vars(id)%tot_i_min(jd).ne.vars(id)%loc_i_min(jd) .or. &
                        &vars(id)%tot_i_max(jd).ne.vars(id)%loc_i_max(jd)) then
                        ind_print = .false.
                        return
                    end if
                end do
            end do
        end subroutine
    end function print_HDF5_arrs
    
    ! Reads a  PB3D output file in  HDF5 format. This happens  in a non-parallel
    ! way. By  default, all  variables are  read, but an  array of  strings with
    ! acceptable variable names can be passed.
    ! Optionally, output can be given about the variable being read.
    ! Also, if  "rich_lvl" is  provided, "_R_rich_lvl" is  appended to  the head
    ! name if it is > 0, and similarly for "eq_job" in "_E_eq_job".
    integer function read_HDF5_arr_ind(var,PB3D_name,head_name,var_name,&
        &rich_lvl,eq_job,disp_info,lim_loc) result(ierr)                        ! individual version
        use HDF5_utilities, only: list_all_vars_in_group, set_1D_vars
        
        character(*), parameter :: rout_name = 'read_HDF5_arr_ind'
        
        ! input / output
        type(var_1D_type), intent(inout) :: var                                 ! variable to read
        character(len=*), intent(in) :: PB3D_name                               ! name of PB3D file
        character(len=*), intent(in) :: head_name                               ! head name of variables
        character(len=*), intent(in) :: var_name                                ! name of variable
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to append to file name
        logical, intent(in), optional :: disp_info                              ! display additional information about variable being read
        integer, intent(in), optional :: lim_loc(:,:)                           ! local limits
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        integer(HID_T) :: HDF5_kind_64                                          ! HDF5 type corresponding to dp
        integer(HID_T) :: group_id                                              ! group identifier
        integer(HID_T) :: head_group_id                                         ! head group identifier
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HID_T) :: memspace                                              ! Dataspace identifier in memory
        integer(HID_T) :: chunk_a_plist_id                                      ! chunk access property list identifier 
        integer(HSIZE_T) :: id                                                  ! counter
        integer(HSIZE_T) :: data_size                                           ! size of data set
        integer(HSIZE_T) :: n_dims                                              ! nr. of dimensions
        integer(SIZE_T) :: name_len                                             ! length of name of group
        integer :: storage_type                                                 ! type of storage used in HDF5 file
        integer :: nr_lnks_head                                                 ! number of links in head group
        integer :: max_corder                                                   ! current maximum creation order value for group
        integer, allocatable :: lim_tot(:,:)                                    ! total limits
        integer, allocatable :: lim_loc_loc(:,:)                                ! local lim_loc
        character(len=max_str_ln) :: name_len_loc                               ! local copy of name_len
        character(len=max_str_ln) :: group_name                                 ! name of group
        character(len=max_str_ln) :: head_name_loc                              ! local head_name
        logical :: disp_info_loc                                                ! local disp_info
        
        ! initialize ierr
        ierr = 0
        
        ! set local disp_info
        disp_info_loc = .false.
        if (present(disp_info)) disp_info_loc = disp_info
        
        ! set up local head_name
        head_name_loc = head_name
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) head_name_loc = trim(head_name_loc)//'_R_'//&
                &trim(i2str(rich_lvl))
        end if
        if (present(eq_job)) then
            if (eq_job.gt.0) head_name_loc = trim(head_name_loc)//'_E_'//&
                &trim(i2str(eq_job))
        end if
        
        ! user output
        if (debug_HDF5_ops) then
            write(*,*) 'Reading data from PB3D output "'//trim(PB3D_name)//&
                &'/'//trim(head_name_loc)//'/'//trim(var_name)//'"'
        end if
        
        ! initialize FORTRAN predefined datatypes
        call H5open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! preparation
        HDF5_kind_64 = H5kind_to_type(dp,H5_REAL_KIND)
        
        ! user output
        if (disp_info_loc) then
            call writo('Opening file '//trim(PB3D_name))
        end if
        
        ! open the file
        call H5Fopen_f(trim(PB3D_name),H5F_ACC_RDONLY_F,HDF5_i,ierr)
        CHCKERR('Failed to open file')
        
        ! user output
        if (disp_info_loc) then
            call writo('Opening variable "'//trim(head_name_loc)//'"')
        end if
        
        ! open head group
        call H5Gopen_f(HDF5_i,trim(head_name_loc),head_group_id,ierr)
        CHCKERR('Failed to open head group')
        
#if ldebug
        ! display all variables in group
        if (disp_info_loc) then
            ierr = list_all_vars_in_group(head_group_id)
            CHCKERR('')
        end if
#endif
        
        ! get number of objects in group to allocate vars
        call H5Gget_info_f(head_group_id,storage_type,nr_lnks_head,max_corder,&
            &ierr)
        CHCKERR('Failed to get group info')
        
        ! iterate over all elements in group to save all acceptable variables
        do id = 1, nr_lnks_head
            call H5Lget_name_by_idx_f(head_group_id,'.',H5_INDEX_NAME_F,&
                &H5_ITER_NATIVE_F,id-1,group_name,ierr,size=name_len)
            CHCKERR('Failed to get name')
            
#if ldebug
            ! check length
            if (name_len.gt.max_str_ln) then
                write(name_len_loc,*) name_len
                ierr = 1
                err_msg = 'Recompile with max_str_ln > '//trim(name_len_loc)
                CHCKERR(err_msg)
            end if
#endif
            
            if (trim(group_name).ne.trim(var_name)) then                        ! not found
                cycle
            end if
            
            ! user output
            if (disp_info_loc) then
                call writo('Reading '//trim(group_name))
            end if
            
            ! set name
            var%var_name = trim(group_name)
            
            ! open group
            call H5Gopen_f(head_group_id,trim(group_name),group_id,ierr)
            CHCKERR('Failed to open group')
            
            ! 1. limit information
            
            ! open limit information dataset
            call h5dopen_f(group_id,'lim',dset_id,ierr)
            CHCKERR('Failed to open dataset')
            
            ! get dataspace
            call H5Dget_space_f(dset_id,filespace,ierr)
            CHCKERR('Failed to get file space')
            
            ! get size to allocate the 1D variable
            call H5Sget_simple_extent_npoints_f(filespace,data_size,ierr)
            n_dims = data_size/2
            CHCKERR('Failed to get storage size')
            allocate(var%tot_i_min(n_dims))
            allocate(var%tot_i_max(n_dims))
            
            ! close dataspace
            call H5Sclose_f(filespace,ierr)
            CHCKERR('Failed to close file space')
            
            ! allocate helper variables
            allocate(lim_tot(n_dims,2))
            allocate(lim_loc_loc(n_dims,2))
            
            ! read into 1D variable
            call H5Dread_f(dset_id,H5T_NATIVE_INTEGER,lim_tot,&
                &[2*n_dims],ierr)
            CHCKERR('Failed to read dataset')
            
            ! set up local limits
            if (present(lim_loc)) then
                lim_loc_loc = lim_loc
            else
                lim_loc_loc = lim_tot
            end if
            
            ! copy into limits of var
            ! Note: The local  limits in the HDF5 file are  the total limits
            ! in the output.
            var%tot_i_min = lim_loc_loc(:,1)
            var%tot_i_max = lim_loc_loc(:,2)
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! 2. variable
            
            ! set up local data size to allocate the 1D variable
            data_size = product(lim_loc_loc(:,2)-lim_loc_loc(:,1)+1)
            allocate(var%p(data_size))
            
            ! set up chunk access property list
            ierr = set_1D_vars(lim_tot,lim_loc_loc,&
                &a_plist_id=chunk_a_plist_id)
            CHCKERR('')
            
            ! open variable dataset
            call H5Dopen_f(group_id,'var',dset_id,ierr,&
                &dapl_id=chunk_a_plist_id)
            CHCKERR('Failed to open dataset')
            
            ! close data access property list
            call H5Pclose_f(chunk_a_plist_id,ierr)
            CHCKERR('Failed to close property list')
            
            ! get dataspace
            call H5Dget_space_f(dset_id,filespace,ierr)
            CHCKERR('Failed to get file space')
            
            ! set 1D filespace hyperslab selection
            ierr = set_1D_vars(lim_tot,lim_loc_loc,space_id=filespace)
            CHCKERR('')
            
            ! create memory data space
            call H5Screate_simple_f(1,[data_size],memspace,ierr)
            CHCKERR('Failed to create memory space')
            
            ! read into 1D variable
            call H5Dread_f(dset_id,HDF5_kind_64,var%p,[data_size],&
                &ierr,mem_space_id=memspace,file_space_id=filespace)
            CHCKERR('Failed to read dataset')
            
            ! close dataspaces
            call H5Sclose_f(filespace,ierr)
            CHCKERR('Failed to close file space')
            call H5Sclose_f(memspace,ierr)
            CHCKERR('Unable to close memory space')
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! end of group
            
            ! deallocate limits
            deallocate(lim_tot,lim_loc_loc)
            
            ! close the group
            call H5Gclose_f(group_id,ierr)
            CHCKERR('Failed to close group')
            
            ! exit loop
            exit
        end do
        
        ! close the head group
        call H5Gclose_f(head_group_id,ierr)
        CHCKERR('Failed to close head group')
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
    end function read_HDF5_arr_ind
    integer function read_HDF5_arr_range(vars,PB3D_name,head_name,var_name,&
        &rich_lvl,eq_job,disp_info,lim_loc) result(ierr)                        ! range version
        character(*), parameter :: rout_name = 'read_HDF5_arr_range'
        
        ! input / output
        type(var_1D_type), intent(inout), allocatable :: vars(:)                ! variables to read
        character(len=*), intent(in) :: PB3D_name                               ! name of PB3D file
        character(len=*), intent(in) :: head_name                               ! head name of variables
        character(len=*), intent(in) :: var_name                                ! name of variable
        integer, intent(in), optional :: rich_lvl                               ! Richardson level to reconstruct
        integer, intent(in), optional :: eq_job(2)                              ! equilibrium job to append to file name
        logical, intent(in), optional :: disp_info                              ! display additional information about variable being read
        integer, intent(in), optional :: lim_loc(:,:)                           ! local limits
        
        ! local variables
        integer :: id                                                           ! counter
        type(var_1D_type) :: var_loc                                            ! local vars
        integer, allocatable :: eq_job_loc(:)                                   ! local eq_job
        
        ! set up local eq_job
        allocate(eq_job_loc(2))
        if (present(eq_job)) then
            eq_job_loc = eq_job
        else
            eq_job_loc = 0
        end if
        
        ! set up vars
        allocate(vars(eq_job(2)-eq_job(1)+1))
        
        ! call individual version for every local eq_job
        do id = eq_job_loc(1),eq_job_loc(2)
            ierr = read_HDF5_arr_ind(var_loc,PB3D_name,head_name,var_name,&
                &rich_lvl=rich_lvl,eq_job=id,disp_info=disp_info,&
                &lim_loc=lim_loc)
            CHCKERR('')
            vars(id-eq_job_loc(1)+1) = var_loc
            call dealloc_var_1D(var_loc)
        end do
    end function read_HDF5_arr_range
    
    ! resets an HDF5 XDMF item
    ! Note: individual version cannot make use of array version because then the
    ! deallocation does not work properly.
    ! [MPI] Only group master
    subroutine reset_HDF5_item_arr(XDMF_items,ind_plot)                         ! array version
        use num_vars, only: rank
        
        ! input / output
        type(XML_str_type), intent(inout) :: XDMF_items(:)                      ! XDMF items to reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_items                                                      ! nr. of items
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! set n_items
        n_items = size(XDMF_items)
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            do id = 1,n_items
                if (.not.allocated(XDMF_items(id)%xml_str)) then
                    call writo('Could not reset HDF5 XDMF item "'&
                        &//trim(XDMF_items(id)%name)//'"',warning=.true.)
                else
                    if (debug_HDF5_ops) write(*,*) 'reset "'//&
                        &trim(XDMF_items(id)%name)//'"'
                    XDMF_items(id)%name = ''
                    deallocate(XDMF_items(id)%xml_str)
                end if
            end do
        end if
    end subroutine reset_HDF5_item_arr
    subroutine reset_HDF5_item_ind(XDMF_item,ind_plot)                          ! individual version
        use num_vars, only: rank
        
        ! input / output
        type(XML_str_type), intent(inout) :: XDMF_item                          ! XDMF item to reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.rank.eq.0) then
            if (.not.allocated(XDMF_item%xml_str)) then
                call writo('Could not reset HDF5 XDMF item "'&
                    &//trim(XDMF_item%name)//'"',warning=.true.)
            else
                if (debug_HDF5_ops) write(*,*) 'reset "'//trim(XDMF_item%name)&
                    &//'"'
                XDMF_item%name = ''
                deallocate(XDMF_item%xml_str)
            end if
        end if
    end subroutine reset_HDF5_item_ind
end module HDF5_ops
