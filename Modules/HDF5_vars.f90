!------------------------------------------------------------------------------!
!   variables and routines pertaining to HDF5 and XDMF                         !
!------------------------------------------------------------------------------!
module HDF5_vars
#include <PB3D_macros.h>
    use num_vars, only: max_str_ln, dp, plot_dir, data_dir, script_dir
    use message_ops, only: writo
    use str_ops, only: i2str, r2str, r2strt
    use HDF5
    
    implicit none
    private
    public init_HDF5, print_HDF5_grid, print_HDF5_geom, print_HDF5_top, &
        &print_HDF5_att, print_HDF5_3D_data_item, open_HDF5_file, &
        &close_HDF5_file, &
        &XML_str_type, HDF5_file_type, add_HDF5_item, reset_HDF5_item
    
    ! global variables
    integer, parameter :: max_xml_ln = 300                                      ! max. length of xml string
    character(len=6) :: xmf_fmt = '(999A)'                                      ! format to write the xmf file
    logical :: debug = .false.                                                   ! set to true to debug information
    
    ! XML strings used in XDMF
    type :: XML_str_type
        character(len=max_str_ln) :: name                                       ! name of this item
        integer :: max_xml_ln = 300                                             ! max. length of xml string
        character(len=max_xml_ln), allocatable :: xml_str(:)                    ! XML string
    end type XML_str_type
    
    ! HDF5 data tyipe
    type :: HDF5_file_type                                                      ! type containing the information about HDF5 files
        integer :: HDF5_i                                                       ! HDF5 file handle
        integer :: XDMF_i                                                       ! XDMF file handle
        character(len=max_str_ln) :: name                                       ! name of files (without extensions ".h5" and ".xmf")
    end type HDF5_file_type
    
    ! XDMF possibilities
    character(len=max_str_ln) :: XDMF_num_types(2)                              ! possible XDMF number types
    character(len=max_str_ln) :: XDMF_format_types(2)                           ! possible XDMF format types
    character(len=max_str_ln) :: XDMF_geom_types(2)                             ! possible XDMF geometry types
    character(len=max_str_ln) :: XDMF_top_types(2)                              ! possible XDMF topology types
    character(len=max_str_ln) :: XDMF_att_types(1)                              ! possible XDMF attribute types
    character(len=max_str_ln) :: XDMF_center_types(2)                           ! possible XDMF attribute center types
    character(len=max_str_ln) :: XDMF_grid_types(3)                             ! possible XDMF grid types
        
    ! interfaces
    interface reset_HDF5_item
        module procedure reset_HDF5_item_ind, reset_HDF5_item_arr
    end interface
    
contains
    ! Initializes the HDF5 types.
    subroutine init_HDF5
        ! XDMF_num_types
        XDMF_num_types(1) = "Int"
        XDMF_num_types(2) = "Float"
        
        ! XDMF_format_types
        XDMF_format_types(1) = "XML"
        XDMF_format_types(2) = "HD5"
        
        ! XDMF_geom_types
        XDMF_geom_types(1) = "X_Y"
        XDMF_geom_types(2) = "X_Y_Z"
        
        ! XDMF_top_types
        XDMF_top_types(1) = "2DSMesh"
        XDMF_top_types(2) = "3DSMesh"
        
        ! XDMF_att_types
        XDMF_att_types(1) = "Scalar"
        
        ! XDMF_center_types
        XDMF_center_types(1) = "Node"
        XDMF_center_types(2) = "Cell"
        
        ! XDMF_grid_types
        XDMF_grid_types(1) = "None"
        XDMF_grid_types(2) = "Temporal"
        XDMF_grid_types(3) = "Spatial"
    end subroutine init_HDF5
    
    ! Opens an HDF5 file and accompanying xmf file and returns the handles.
    ! Optionally, a description of the file  can be provided. Also, the plot can
    ! be done for only one process, setting the variable "ind_plot" to .true.
    ! [MPI] Parts by all processes, parts only by group master
    integer function open_HDF5_file(file_info,file_name,description,&
        &ind_plot) result(ierr)
        use num_vars, only: MPI_Comm_groups, grp_rank
        use safe_open_mod, only: safe_open
        use MPI
        
        character(*), parameter :: rout_name = 'open_HDF5_file'
        
        ! input / output
        type(HDF5_file_type) :: file_info                                       ! info about HDF5 file
        character(len=*), intent(in) :: file_name                               ! name of HDF5 file
        character(len=*), intent(in), optional  :: description                  ! description of file
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        character(len=max_str_ln) :: full_file_name                             ! full file name
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: plist_id                                              ! property list identifier 
        integer :: MPI_Comm                                                     ! MPI Communicator used
        
        ! initialize ierr
        ierr = 0
        
        ! set up MPI Communicator
        MPI_Comm = MPI_Comm_groups                                              ! default group communicator
        if (present(ind_plot)) then
            if (ind_plot) MPI_Comm = MPI_Comm_self                              ! individual plot
        end if
        
        ! set up full file name
        full_file_name = data_dir//'/'//trim(file_name)//'.h5'
        
        ! initialize FORTRAN predefined datatypes
        call H5open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! setup file access property list with parallel I/O access if needed
        call H5Pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
        CHCKERR('Failed to create property list')
        call H5Pset_fapl_mpio_f(plist_id,MPI_Comm,MPI_INFO_NULL,ierr)
        CHCKERR('Failed to set file access property')
        
        ! create the file collectively.
        call H5Fcreate_f(trim(full_file_name),H5F_ACC_TRUNC_F,HDF5_i,ierr,&
            &access_prp=plist_id)
        CHCKERR('Failed to create file')
        call H5Pclose_f(plist_id,ierr)
        CHCKERR('Failed to close property list')
        
        ! update file info, converting integers and setting name
        file_info%HDF5_i = HDF5_i
        file_info%name = file_name
        
        ! user output
        call writo('HDF5 file "'//trim(full_file_name)//'" created')
            
        ! only if group master
        if (grp_rank.eq.0) then
            ! open accompanying xmf file
            full_file_name = data_dir//'/'//trim(file_name)//'.xmf'
            call safe_open(file_info%XDMF_i,ierr,full_file_name,'replace',&
                &'formatted',delim_in='none')
            CHCKERR('Failed to open xmf file')
            
            ! user output
            call writo('XDMF file "'//trim(full_file_name)//'" created')
            
            ! write header if group master
            write(file_info%XDMF_i,xmf_fmt) '<?xml version="1.0" ?>'
            write(file_info%XDMF_i,xmf_fmt) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" &
                &[]>'
            write(file_info%XDMF_i,xmf_fmt) '<Xdmf Version="2.0">'
            write(file_info%XDMF_i,xmf_fmt) '<Domain>'
            if (present(description)) then
                write(file_info%XDMF_i,xmf_fmt) &
                    &'<Information Name="Description">'
                write(file_info%XDMF_i,xmf_fmt) description
                write(file_info%XDMF_i,xmf_fmt) '</Information>'
            end if
        end if
    end function open_HDF5_file
    
    ! Closes an HDF5 file and writes the accompanying xmf file
    ! [MPI] Parts by all processes, parts only by group master
    integer function close_HDF5_file(file_info) result(ierr)
        use num_vars, only: grp_rank
        
        character(*), parameter :: rout_name = 'close_HDF5_file'
        
        ! input / output
        type(HDF5_file_type) :: file_info                                       ! info about HDF5 file
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: full_file_name                             ! full file name
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        
        ! initialize ierr
        ierr = 0
        
        ! set full file name and HDF5_i (converting integers)
        full_file_name = data_dir//'/'//trim(file_info%name)//'.h5'
        HDF5_i = file_info%HDF5_i
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
        
        ! user output
        call writo('HDF5 file "'//trim(full_file_name)//'" closed')
            
        ! only if group master
        if (grp_rank.eq.0) then
            ! close header
            if (grp_rank.eq.0) then
                write(file_info%XDMF_i,xmf_fmt) '</Domain>'
                write(file_info%XDMF_i,xmf_fmt) '</Xdmf>'
            end if
            
            ! set full file name
            full_file_name = data_dir//'/'//trim(file_info%name)//'.xmf'
            
            ! close accompanying xmf file
            close(file_info%XDMF_i,IOSTAT=ierr)
            CHCKERR('Failed to close xmf file')
            
            ! user output
            call writo('XDMF file "'//trim(full_file_name)//'" closed')
        end if
    end function close_HDF5_file
    
    ! Add an XDMF item to a HDF5 file
    ! Note:  This  should  only  be  used  with  grids  (or for  topologies  and
    ! geometries that are used throughout)
    ! [MPI] Only group master
    subroutine add_HDF5_item(file_info,XDMF_item,reset)
        use num_vars, only: grp_rank
        
        ! input / output
        type(HDF5_file_type) :: file_info                                       ! info about HDF5 file
        type(XML_str_type) :: XDMF_item                                         ! XDMF item to add
        logical, intent(in), optional :: reset                                  ! if .true., XDMF_item is reset
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: item_len                                                     ! length of item
        logical :: reset_loc                                                    ! local copy of reset
        
        ! only if group master
        if (grp_rank.eq.0) then
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
                write(file_info%XDMF_i,xmf_fmt) XDMF_item%xml_str(id)
            end do
            
            ! reset if requested
            if (reset_loc) call reset_HDF5_item(XDMF_item)
        end if
    end subroutine add_HDF5_item
    
    ! resets an HDF5 XDMF item
    ! Note: individual version cannot make use of array version because then the
    ! deallocation does not work properly>
    ! [MPI] Only group master
    subroutine reset_HDF5_item_arr(XDMF_items)                                  ! array vesion
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) :: XDMF_items(:)                                     ! XDMF items to reset
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_items                                                      ! nr. of items
        
        ! set n_items
        n_items = size(XDMF_items)
        
        ! only if group master
        if (grp_rank.eq.0) then
            do id = 1,n_items
                if (.not.allocated(XDMF_items(id)%xml_str)) then
                    call writo('WARNING: Could not reset HDF5 XDMF item "'&
                        &//trim(XDMF_items(id)%name)//'"')
                else
                    if (debug) write(*,*) 'reset "'//&
                        &trim(XDMF_items(id)%name)//'"'
                    XDMF_items(id)%name = ''
                    deallocate(XDMF_items(id)%xml_str)
                end if
            end do
        end if
    end subroutine reset_HDF5_item_arr
    subroutine reset_HDF5_item_ind(XDMF_item)                                   ! individual vesion
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) :: XDMF_item                                         ! XDMF item to reset
        
        ! only if group master
        if (grp_rank.eq.0) then
            if (.not.allocated(XDMF_item%xml_str)) then
                call writo('WARNING: Could not reset HDF5 XDMF item "'&
                    &//trim(XDMF_item%name)//'"')
            else
                if (debug) write(*,*) 'reset "'//trim(XDMF_item%name)//'"'
                XDMF_item%name = ''
                deallocate(XDMF_item%xml_str)
            end if
        end if
    end subroutine reset_HDF5_item_ind
    
    ! prints an HDF5 data item
    ! If this is a parallel data item, the group dimension and offset have to be
    ! specified as well.
    ! [MPI] Parts by all processes, parts only by group master
    integer function print_HDF5_3D_data_item(dataitem_id,file_info,var_name,&
        &var,tot_dim,grp_dim,grp_offset) result(ierr)
        use num_vars, only: grp_rank
        
        character(*), parameter :: rout_name = 'print_HDF5_3D_data_item'
        
        ! input / output
        type(XML_str_type) :: dataitem_id                                       ! ID of data item
        type(HDF5_file_type) :: file_info                                       ! info about HDF5 file
        character(len=*), intent(in) :: var_name                                ! name of variable
        real(dp), intent(in) :: var(:,:,:)                                      ! variable to write
        integer, intent(in) :: tot_dim(3)                                       ! total dimensions of variable
        integer, intent(in), optional :: grp_dim(3)                             ! dimensions in this group
        integer, intent(in), optional :: grp_offset(3)                          ! offset in this group
        
        ! local variables
        integer :: grp_dim_loc(3)                                               ! local copy of grp_dim
        integer :: grp_offset_loc(3)                                            ! local copy of grp_offset
        integer(HSIZE_T) :: dimsf(3)                                            ! total dataset dimensions
        integer(HSIZE_T) :: dimsm(3)                                            ! group dataset dimensions
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HID_T) :: memspace                                              ! Dataspace identifier in memory
        integer(HID_T) :: plist_id                                              ! property list identifier 
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        integer(HSSIZE_T) :: mem_offset(3)                                      ! offset in memory
        integer(HSIZE_T) :: mem_block(3)                                        ! block size in memory
        integer(HSIZE_T) :: mem_stride(3)                                       ! stride in memory
        integer(HSIZE_T) :: mem_count(3)                                        ! nr. of repetitions of block in memory
        
        ! initialize ierr
        ierr = 0
        
        ! set the group dimensions and offset
        ierr = check_for_parallel_3D(tot_dim,grp_dim_loc,grp_offset_loc,&
            &grp_dim,grp_offset)
        CHCKERR('')
        
        ! find out if var_name contains a group name ("/")
        !!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!
        
        ! create the data spaces for the dataset. 
        dimsf = tot_dim
        dimsm = grp_dim_loc
        call H5Screate_simple_f(size(dimsf),dimsf,filespace,ierr)
        CHCKERR('Failed to create file space')
        call H5Screate_simple_f(size(dimsm),dimsm,memspace,ierr)
        CHCKERR('Failed to create memory space')
        
        ! create chunked dataset.
        !call H5Pcreate_f(H5P_DATASET_CREATE_F,plist_id,ierr)
        !CHCKERR('Failed to create property list')
        !dimsm = file_info%grp_dim
        !call H5Pset_chunk_f(plist_id,size(dimsm),dimsm,ierr)
        !CHCKERR('Failed to set chunk property')
        !call H5Dcreate_f(file_info%HDF5_i,trim(var_name),H5T_NATIVE_DOUBLE,&
            !&filespace,dset_id,ierr,plist_id)
        call H5Dcreate_f(file_info%HDF5_i,trim(var_name),H5T_NATIVE_DOUBLE,&
            &filespace,dset_id,ierr)                                            ! DOESN'T SEEM TO WORK WITH CHUNKED PROPERTY!!!
        CHCKERR('Failed to create file data set')
        call H5Sclose_f(filespace,ierr)
        CHCKERR('Failed to close file space')
        !call h5pclose_f(plist_id,ierr)
        !CHCKERR('Failed to close property list')
        
        ! select hyperslab in the file.
        mem_count = 1
        mem_stride = 1
        mem_block = grp_dim_loc
        mem_offset = grp_offset_loc
        call H5Dget_space_f(dset_id,filespace,ierr)
        CHCKERR('Failed to get file space')
        call H5Sselect_hyperslab_f (filespace,H5S_SELECT_SET_F,mem_offset,&
            &mem_count,ierr,mem_stride,mem_block)
        CHCKERR('Failed to select chunk hyperslab')
        
        ! create property list for collective dataset write
        call H5Pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr) 
        CHCKERR('Failed to create property list')
        call H5Pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
        CHCKERR('Failed to set parallel property')
        
        ! write the dataset collectively. 
        dimsf = tot_dim
        call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,dimsf,ierr,&
            &file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
        CHCKERR('Failed to write data set')
        call H5Pclose_f(plist_id,ierr)
        CHCKERR('Failed to create property list')
        
        ! close dataspaces.
        call H5Sclose_f(filespace,ierr)
        CHCKERR('Unable to close file space')
        call H5Sclose_f(memspace,ierr)
        CHCKERR('Unable to close memory space')
        
        ! close the dataset.
        call H5Dclose_f(dset_id,ierr)
        CHCKERR('Failed to close data set')
        
        ! set XDMF dataitem ID if group master
        if (grp_rank.eq.0) then
            dataitem_id%name = 'DataItem - '//trim(var_name)
            allocate(dataitem_id%xml_str(3))
            dataitem_id%xml_str(1) = '<DataItem Dimensions="'//&
                &trim(i2str(tot_dim(3)))//' '//&
                &trim(i2str(tot_dim(2)))//' '//&
                &trim(i2str(tot_dim(1)))//&
                &'" NumberType="Float" Precision="8" Format="HDF">'
            dataitem_id%xml_str(2) = trim(file_info%name)//'.h5:/'//&
                &trim(var_name)
            dataitem_id%xml_str(3) = '</DataItem>'
            
            if (debug) write(*,*) 'created data item "'//&
                &trim(dataitem_id%name)//'"'
        end if
    contains
        ! check whether the  dimensions provided are a  valid parallel indicator
        ! or not
        integer function check_for_parallel_3D(tot_dim,grp_dim_out,&
            &grp_offset_out,grp_dim_in,grp_offset_in) result(ierr)
            character(*), parameter :: rout_name = 'check_for_parallel_3D'
            
            ! input / output
            integer, intent(in) :: tot_dim(3)                                   ! total dimension
            integer, intent(inout) :: grp_dim_out(3), grp_offset_out(3)         ! output group dimension and offset
            integer, intent(in), optional :: grp_dim_in(3), grp_offset_in(3)    ! input group dimension and offset
            
            ! local variables
            character(len=max_str_ln) :: err_msg                                ! error message
            integer :: id                                                       ! counter
            
            ! initialize ierr
            ierr = 0
            
            if (present(grp_dim_in)) then
                if (.not.present(grp_offset_in)) then
                    err_msg = 'Need to specify offset as well as group &
                        &dimensions'
                    ierr = 1
                    CHCKERR(err_msg)
                else
                    grp_dim_out = grp_dim_in
                    grp_offset_out = grp_offset_in
                end if
                do id = 1,size(tot_dim)
                    if (grp_dim_in(id).gt.tot_dim(id)) then                     ! error in parallel file
                        err_msg = 'Total dimension '//trim(i2str(id))//&
                            &' cannot be smaller than group dimension'
                        ierr = 1
                        CHCKERR(err_msg)
                    end if
                end do
            else
                grp_dim_out = tot_dim
                grp_offset_out = 0
            end if
        end function check_for_parallel_3D
    end function print_HDF5_3D_data_item
    
    ! prints an HDF5 attribute
    ! [MPI] Only group master
    subroutine print_HDF5_att(att_id,att_dataitem,att_name,att_center,reset)
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) :: att_id                                            ! ID of attribute
        type(XML_str_type) :: att_dataitem                                      ! dataitem of attribute
        character(len=*), intent(in) :: att_name                                ! name of attribute
        integer, intent(in) :: att_center                                       ! center of attribute
        logical, intent(in), optional :: reset                                  ! if .true., data items are reset
        
        ! local variables
        integer :: dataitem_len                                                 ! length of data item
        integer :: id                                                           ! counter
        logical :: reset_loc                                                    ! local copy of reset
        
        ! only if group master
        if (grp_rank.eq.0) then
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
            
            if (debug) write(*,*) 'created attribute "'//trim(att_id%name)//'"'
            
            ! reset if requested
            if (reset_loc) call reset_HDF5_item(att_dataitem)
        end if
    end subroutine print_HDF5_att
    
    ! prints an HDF5 topology
    ! Note: currently only structured grids supported
    ! [MPI] Only group master
    subroutine print_HDF5_top(top_id,top_type,top_n_elem)
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) ::  top_id                                           ! ID of topology
        integer, intent(in) :: top_type                                         ! type
        integer, intent(in) :: top_n_elem(:)                                    ! nr. of elements
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_dims                                                       ! nr. of dimensions
        character(len=max_str_ln) :: work_str                                   ! work string
        
        ! only if group master
        if (grp_rank.eq.0) then
            ! set n_dims
            n_dims = size(top_n_elem)
            
            ! fill work string
            work_str = ''
            do id = n_dims,1,-1
                work_str = trim(work_str)//' '//trim(i2str(top_n_elem(id)))
            end do
            
            ! set XDMF topology ID
            top_id%name = 'Topology'
            allocate(top_id%xml_str(1))
            top_id%xml_str(1) = '<Topology TopologyType="'//&
                &trim(XDMF_top_types(top_type))//'" NumberOfElements="'//&
                &trim(work_str)//'"/>'
            
            if (debug) write(*,*) 'created topology "'//trim(top_id%name)//'"'
        end if
    end subroutine print_HDF5_top
    
    ! prints an HDF5 geometry
    ! [MPI] Only group master
    subroutine print_HDF5_geom(geom_id,geom_type,geom_dataitems,reset)
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) ::  geom_id                                          ! ID of geometry
        integer, intent(in) :: geom_type                                        ! type of geometry
        type(XML_str_type) :: geom_dataitems(:)                                 ! data items of geometry
        logical, intent(in), optional :: reset                                  ! if .true., data items are reset
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: id_sum                                                       ! counter
        integer, allocatable :: dataitem_len(:)                                 ! length of data item
        integer :: n_dataitems                                                  ! nr. of data items
        logical :: reset_loc                                                    ! local copy of reset
        
        ! only if group master
        if (grp_rank.eq.0) then
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
            
            if (debug) write(*,*) 'created geometry "'//trim(geom_id%name)//'"'
            
            ! reset if requested
            if (reset_loc) call reset_HDF5_item(geom_dataitems)
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
        &grid_top,grid_geom,grid_atts,grid_grids,reset) result(ierr)
        use num_vars, only: grp_rank
        
        character(*), parameter :: rout_name = 'print_HDF5_grid'
        
        ! input / output
        type(XML_str_type) :: grid_id                                           ! ID of grid
        character(len=*), intent(in) :: grid_name                               ! name
        integer, intent(in) :: grid_type                                        ! type
        real(dp), intent(in), optional :: grid_time                             ! time of grid
        type(XML_str_type), optional :: grid_top                                ! topology
        type(XML_str_type), optional :: grid_geom                               ! geometry
        type(XML_str_type), optional :: grid_atts(:)                            ! attributes
        type(XML_str_type), optional :: grid_grids(:)                           ! grids
        logical, intent(in), optional :: reset                                  ! if .true., top, geom and atts or grids are reset
        
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
        
        ! initialize ierr
        ierr = 0
        
        ! only if group master
        if (grp_rank.eq.0) then
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
            
            if (debug) write(*,*) 'created grid "'//trim(grid_id%name)//'"'
            
            ! reset if requested
            if (reset_loc) then
                if (grid_type.eq.1)  then
                    if (present(grid_top)) call reset_HDF5_item(grid_top)
                    if (present(grid_geom)) call reset_HDF5_item(grid_geom)
                    if (present(grid_atts)) call reset_HDF5_item(grid_atts)
                else if (grid_type.eq.2 .or. grid_type.eq.3) then
                    if (present(grid_grids)) call reset_HDF5_item(grid_grids)
                end if
            end if
        end if
    end function print_HDF5_grid
end module HDF5_vars
