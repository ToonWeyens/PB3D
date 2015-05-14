!------------------------------------------------------------------------------!
!   Operations and variables pertaining to HDF5 and XDMF                       !
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
!------------------------------------------------------------------------------!
module HDF5_ops
#include <PB3D_macros.h>
    use num_vars, only: max_str_ln, dp, plot_dir, data_dir, script_dir
    use messages, only: writo, lvl_ud
    use str_ops, only: i2str, r2str, r2strt
    use HDF5
    
    implicit none
    private
    public init_HDF5, print_HDF5_grid, print_HDF5_geom, print_HDF5_top, &
        &print_HDF5_att, print_HDF5_3D_data_item, open_HDF5_file, &
        &close_HDF5_file, reset_HDF5_item, add_HDF5_item, create_output_HDF5, &
        &print_HDF5_arrs, read_HDF5_arrs, &
        &XML_str_type, HDF5_file_type, var_1D
    
    ! global variables
    integer, parameter :: max_xml_ln = 300                                      ! max. length of xml string
    character(len=6) :: xmf_fmt = '(999A)'                                      ! format to write the xmf file
    logical :: debug_HDF5_ops = .false.                                         ! set to true to debug information
    
    ! XML strings used in XDMF
    type :: XML_str_type
        character(len=max_str_ln) :: name                                       ! name of this item
        integer :: max_xml_ln = 300                                             ! max. length of xml string
        character(len=max_xml_ln), allocatable :: xml_str(:)                    ! XML string
    end type XML_str_type
    
    ! HDF5 data tipe
    type :: HDF5_file_type                                                      ! type containing the information about HDF5 files
        integer :: HDF5_i                                                       ! HDF5 file handle
        integer :: XDMF_i                                                       ! XDMF file handle
        character(len=max_str_ln) :: name                                       ! name of files (without extensions ".h5" and ".xmf")
    end type HDF5_file_type
    
    ! 1D equivalent of multidimensional variables
    type var_1D
        real(dp), allocatable :: p(:)                                           ! 1D equivalent of data of variable
        integer, allocatable :: tot_i_min(:), tot_i_max(:)                      ! total min. and max. of indices of variable
        integer, allocatable :: grp_i_min(:), grp_i_max(:)                      ! group min. and max. of indices of variable
        character(len=max_str_ln) :: var_name                                   ! name of variable
    end type var_1D
    
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
        use MPI
        use files_utilities, only: nextunit
        
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
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! initialize ierr
        ierr = 0
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! set up MPI Communicator
        if (ind_plot_loc) then
            MPI_Comm = MPI_Comm_self                                            ! individual plot
        else
            MPI_Comm = MPI_Comm_groups                                          ! default group communicator
        end if
        
        ! set up full file name
        full_file_name = data_dir//'/'//trim(file_name)
        
        ! initialize FORTRAN predefined datatypes
        call H5open_f(ierr) 
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
            ! open accompanying xmf file
            open(unit=nextunit(file_info%XDMF_i),&
                &file=trim(full_file_name)//'.xmf',iostat=ierr)
            CHCKERR('Failed to open xmf file')
            
            ! user output
            call writo('HDF5 file "'//trim(full_file_name)//'.h5" created',&
                &persistent=.true.)
            call writo('XDMF file "'//trim(full_file_name)//'.xmf" created',&
                &persistent=.true.)
            
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
        use num_vars, only: grp_rank
        
        character(*), parameter :: rout_name = 'close_HDF5_file'
        
        ! input / output
        type(HDF5_file_type) :: file_info                                       ! info about HDF5 file
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
            ! close header
            write(file_info%XDMF_i,xmf_fmt) '</Domain>'
            write(file_info%XDMF_i,xmf_fmt) '</Xdmf>'
            
            ! close accompanying xmf file
            close(file_info%XDMF_i,IOSTAT=ierr)
            CHCKERR('Failed to close xmf file')
            
            ! user output
            call writo('HDF5 file "'//trim(full_file_name)//'.h5" closed',&
                &persistent=.true.)
            call writo('XDMF file "'//trim(full_file_name)//'.xmf" closed',&
                &persistent=.true.)
        end if
    end function close_HDF5_file
    
    ! Add an XDMF item to a HDF5 file
    ! Note:  This  should  only  be  used  with  grids  (or for  topologies  and
    ! geometries that are used throughout)
    ! [MPI] Only group master
    subroutine add_HDF5_item(file_info,XDMF_item,reset,ind_plot)
        use num_vars, only: grp_rank
        
        ! input / output
        type(HDF5_file_type) :: file_info                                       ! info about HDF5 file
        type(XML_str_type) :: XDMF_item                                         ! XDMF item to add
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
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
            if (reset_loc) call reset_HDF5_item(XDMF_item,ind_plot_loc)
        end if
    end subroutine add_HDF5_item
    
    ! resets an HDF5 XDMF item
    ! Note: individual version cannot make use of array version because then the
    ! deallocation does not work properly>
    ! [MPI] Only group master
    subroutine reset_HDF5_item_arr(XDMF_items,ind_plot)                         ! array version
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) :: XDMF_items(:)                                     ! XDMF items to reset
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
            do id = 1,n_items
                if (.not.allocated(XDMF_items(id)%xml_str)) then
                    call writo('WARNING: Could not reset HDF5 XDMF item "'&
                        &//trim(XDMF_items(id)%name)//'"')
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
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) :: XDMF_item                                         ! XDMF item to reset
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        logical :: ind_plot_loc = .false.                                       ! local version of ind_plot
        
        ! set up local ind_plot
        if (present(ind_plot)) ind_plot_loc = ind_plot
        
        ! only group master if parallel plot or current rank if individual plot
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
            if (.not.allocated(XDMF_item%xml_str)) then
                call writo('WARNING: Could not reset HDF5 XDMF item "'&
                    &//trim(XDMF_item%name)//'"')
            else
                if (debug_HDF5_ops) write(*,*) 'reset "'//trim(XDMF_item%name)&
                    &//'"'
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
        &var,tot_dim,grp_dim,grp_offset,init_val,ind_plot) result(ierr)
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
        real(dp), intent(in), optional :: init_val                              ! initial fill value
        logical, intent(in), optional :: ind_plot                               ! .true. if not a collective plot
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: grp_dim_loc(3)                                               ! local copy of grp_dim
        integer :: grp_offset_loc(3)                                            ! local copy of grp_offset
        integer(HSIZE_T) :: dimsf(3)                                            ! total dataset dimensions
        integer(HSIZE_T) :: dimsm(3)                                            ! group dataset dimensions
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HID_T) :: memspace                                              ! Dataspace identifier in memory
        integer(HID_T) :: plist_id                                              ! property list identifier 
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        !!!integer(HSSIZE_T) :: mem_offset(3)                                      ! offset in memory
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
        HDF5_kind_64 = H5kind_to_type(dp,H5_REAL_KIND)
        
        ! set the group dimensions and offset
        ierr = check_for_parallel_3D(tot_dim,grp_dim_loc,grp_offset_loc,&
            &grp_dim,grp_offset)
        CHCKERR('')
        
        ! create the data spaces for the dataset. 
        dimsf = tot_dim
        dimsm = grp_dim_loc
        call H5Screate_simple_f(size(dimsf),dimsf,filespace,ierr)
        CHCKERR('Failed to create file space')
        call H5Screate_simple_f(size(dimsm),dimsm,memspace,ierr)
        CHCKERR('Failed to create memory space')
        
        !! create chunked dataset.
        !dimsm = file_info%grp_dim
        !call H5Pset_chunk_f(plist_id,size(dimsm),dimsm,ierr)
        !CHCKERR('Failed to set chunk property')                                 ! DOESN'T SEEM TO WORK WITH CHUNKED PROPERTY!!!
        
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
        else
            call H5Dcreate_f(file_info%HDF5_i,trim(var_name),HDF5_kind_64,&
                &filespace,dset_id,ierr)
            CHCKERR('Failed to create file data set')
        end if
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
            ! set up dimensions string
            dim_str = ''
            do id = 1,size(tot_dim)
                if (tot_dim(id).gt.1) then
                    dim_str = trim(dim_str)//' '//trim(i2str(tot_dim(id)))
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
    subroutine print_HDF5_att(att_id,att_dataitem,att_name,att_center,reset,&
        &ind_plot)
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) :: att_id                                            ! ID of attribute
        type(XML_str_type) :: att_dataitem                                      ! dataitem of attribute
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
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
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) ::  top_id                                           ! ID of topology
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
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
        use num_vars, only: grp_rank
        
        ! input / output
        type(XML_str_type) ::  geom_id                                          ! ID of geometry
        integer, intent(in) :: geom_type                                        ! type of geometry
        type(XML_str_type) :: geom_dataitems(:)                                 ! data items of geometry
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
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
        if (ind_plot_loc .or. .not.ind_plot_loc.and.grp_rank.eq.0) then
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
    integer function create_output_HDF5() result(ierr)
        use num_vars, only: grp_rank, output_name, alpha_job_nr, n_alpha
        
        character(*), parameter :: rout_name = 'create_output_HDF5'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        character(len=max_str_ln) :: full_output_name                           ! full name
        
        ! initialize ierr
        ierr = 0
        
        ! set full output name
        full_output_name = trim(output_name)
        if (n_alpha.gt.1) full_output_name = &
            &trim(full_output_name)//'_A'//trim(i2str(alpha_job_nr))            ! append alpha job number
        full_output_name = trim(full_output_name)//'.h5'
        
        ! only for group master
        if (grp_rank.eq.0) then
            ! initialize FORTRAN predefined datatypes
            call H5open_f(ierr) 
            CHCKERR('Failed to initialize HDF5')
            
            ! create the file by group master
            call H5Fcreate_f(full_output_name,H5F_ACC_TRUNC_F,HDF5_i,&
                &ierr)
            CHCKERR('Failed to create file')
            
            ! close the HDF5 file
            call H5Fclose_f(HDF5_i,ierr)
            CHCKERR('failed to close HDF5 file')
            
            ! close FORTRAN interfaces and HDF5 library.
            call H5Close_f(ierr)
            err_msg = 'Failed to close FORTRAN HDF5 interface'
            CHCKERR(err_msg)
        end if
        
        ! user output
        call writo('HDF5 output file '//trim(full_output_name)//' created')
    end function create_output_HDF5
    
    ! Prints a series of arrays, in the form of an array of pointers, to an HDF5
    ! file.
    ! Note: See https://www.hdfgroup.org/HDF5/doc/UG/12_Dataspaces.html, 7.4.2.3
    ! for an explanation of the selection of the dataspaces.
    integer function print_HDF5_arrs(vars,head_name) result(ierr)
        use num_vars, only: MPI_Comm_groups, grp_n_procs, output_name, &
            &alpha_job_nr, n_alpha, grp_rank
        use messages, only: lvl_ud
        use MPI
        
        character(*), parameter :: rout_name = 'print_HDF5_arrs'
        
        ! input / output
        type(var_1D), intent(in) :: vars(:)                                     ! variables to write
        character(len=*), intent(in) :: head_name                               ! head name of variables
        
        ! local variables
        integer :: id, jd                                                       ! counter
        integer :: MPI_Comm                                                     ! MPI Communicator used
        integer :: prev_dims                                                    ! total dimensions previous to divided dimension
        integer :: next_dims                                                    ! total dimensions next to divided dimension
        integer :: div_dim                                                      ! index of divided dimension
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: full_output_name                           ! full name
        integer(HID_T) :: plist_id                                              ! property list identifier 
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: HDF5_kind_64                                          ! HDF5 type corresponding to dp
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HID_T) :: memspace                                              ! Dataspace identifier in memory
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        integer(HID_T) :: group_id                                              ! group identifier
        integer(HID_T) :: head_group_id                                         ! head group identifier
        integer(HSIZE_T) :: dimsf(1)                                            ! total dataset dimensions
        integer(HSIZE_T) :: dimsm(1)                                            ! group dataset dimensions
        integer(HSIZE_T) :: mem_block(1)                                        ! block size in memory
        integer(HSIZE_T) :: file_block(1)                                       ! block size in file
        integer(HSIZE_T) :: mem_stride(1)                                       ! stride in memory
        integer(HSIZE_T) :: mem_count(1)                                        ! nr. of repetitions of block in memory
        integer(HSIZE_T) :: mem_offset(1)                                       ! offset in memory
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Preparing...')
        
        ! set up full output name
        full_output_name = trim(output_name)
        if (n_alpha.gt.1) full_output_name = &
            &trim(full_output_name)//'_A'//trim(i2str(alpha_job_nr))            ! append alpha job number
        full_output_name = trim(full_output_name)//'.h5'
        
        ! set up MPI Communicator
        if (grp_n_procs.eq.1) then
            MPI_Comm = MPI_Comm_self                                            ! individual plot
        else
            MPI_Comm = MPI_Comm_groups                                          ! default group communicator
        end if
        
        ! initialize FORTRAN predefined datatypes
        call H5open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! setup file access property list with parallel I/O access if needed
        call H5Pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
        CHCKERR('Failed to create property list')
        call H5Pset_fapl_mpio_f(plist_id,MPI_Comm,MPI_INFO_NULL,ierr)
        CHCKERR('Failed to set file access property')
        
        ! open the file collectively
        call H5Fopen_f(trim(full_output_name),H5F_ACC_RDWR_F,HDF5_i,ierr,&
            &access_prp=plist_id)
        CHCKERR('Failed to open file')
        call H5Pclose_f(plist_id,ierr)
        CHCKERR('Failed to close property list')
        ! set HDF5 type corresponding to dp
        HDF5_kind_64 = H5kind_to_type(dp,H5_REAL_KIND)
        
        ! create head group
        call H5gcreate_f(HDF5_i,trim(head_name),head_group_id,ierr)
        CHCKERR('Failed to create group')
        
        ! user output
        call writo('Writing variables '//trim(head_name)//'...')
        call lvl_ud(1)
        
        ! loop over all elements in vars
        do id = 1,size(vars)
            ! user output
            call writo('Writing '//trim(vars(id)%var_name))
            
            ! create group
            call H5gcreate_f(head_group_id,trim(vars(id)%var_name),group_id,&
                &ierr)
            CHCKERR('Failed to create group')
            
            ! 1. Data itself
            
            ! calculate memory variables block, offset, count and stride
            call calc_bocs
            
            ! create file and memory data space for the variable dataset
            dimsf = file_block                                                  ! file space has total dimensions
            call H5Screate_simple_f(1,dimsf*next_dims,filespace,ierr)
            CHCKERR('Failed to create file space')
            dimsm = mem_block                                                   ! memory space has group dimensions
            call H5Screate_simple_f(1,dimsm*next_dims,memspace,ierr)
            CHCKERR('Failed to create memory space')
            
            ! create file data set in group
            call H5Dcreate_f(group_id,'var',HDF5_kind_64,filespace,dset_id,&
                &ierr)
            CHCKERR('Failed to create file data set')
            
            ! select first hyperslab in the file
            call H5Sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,&
                &mem_offset,mem_count,ierr,mem_stride,mem_block)
            CHCKERR('Failed to select hyperslab')
            
            ! iterate over all the next dimensions to sum the other hyperslabs
            do jd = 2,next_dims
                ! select hyperslab in the file.
                call H5Sselect_hyperslab_f(filespace,H5S_SELECT_OR_F,&
                    &mem_offset+(jd-1)*dimsf,mem_count,ierr,mem_stride,&
                    &mem_block)
                CHCKERR('Failed to select hyperslab')
            end do
            
            ! create property list for collective dataset write
            call H5Pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr) 
            CHCKERR('Failed to create property list')
            call H5Pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
            CHCKERR('Failed to set parallel property')
            
            ! write the dataset collectively. 
            call H5Dwrite_f(dset_id,HDF5_kind_64,vars(id)%p,dimsf*next_dims,&
                &ierr,file_space_id=filespace,mem_space_id=memspace,&
                &xfer_prp=plist_id)
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
            
            ! 2. Limit information
            
            ! create file data space for the limit dataset
            dimsf = 2*size(vars(id)%tot_i_min)                                  ! min. and max. limit value per dimension
            call H5Screate_simple_f(1,dimsf,filespace,ierr)
            CHCKERR('Failed to create file space')
            
            ! create file data set in group
            call H5Dcreate_f(group_id,'lim',H5T_NATIVE_INTEGER,filespace,&
                &dset_id,ierr)
            CHCKERR('Failed to create file data set')
            
            ! only group leader writes
            if (grp_rank.eq.0) then
                ! write the dataset
                call H5Dwrite_f(dset_id,H5T_NATIVE_INTEGER,&
                    &[vars(id)%tot_i_min,vars(id)%tot_i_max],dimsf,ierr)
                CHCKERR('Failed to write data set')
            end if
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! close the group
            call H5gclose_f(group_id,ierr)
            CHCKERR('Failed to close group')
        end do
        
        call lvl_ud(-1)
        
        ! close the head group
        call H5gclose_f(head_group_id,ierr)
        CHCKERR('Failed to close group')
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
        
        ! user output
        call writo('Variables written in '//trim(full_output_name)//'...')
    contains
        ! Calculate the 1D memory block, offset, count and stride:
        !   -  block corresponds  to the  group size  in the  divided dimension,
        !   multiplied by the total sizes of the previous dimensions
        !   - offset corresponds to i_min in the divided dimension multiplied by
        !   the total sizes of the previous dimensions
        !   - count is one
        !   - stride is one
        subroutine calc_bocs
            ! local variables
            integer :: kd                                                       ! counter
            
            ! intialize local variables
            prev_dims = 1
            next_dims = 1
            div_dim = 0
            
            ! get divided dimension
            if (grp_n_procs.gt.1) then
                ! find divided dimension
                do kd = 1,size(vars(id)%tot_i_min)
                    if (vars(id)%grp_i_max(kd)-vars(id)%grp_i_min(kd) .lt. &
                        &vars(id)%tot_i_max(kd)-vars(id)%tot_i_min(kd)) then    ! dimension is divided
                        if (div_dim.le.0) then                                  ! first divided dimension
                            div_dim = kd
                        else                                                    ! already divided dimension
                            ierr = 1
                            err_msg = 'Only one divided dimension possible'
                            CHCKERR(err_msg)
                        end if
                    end if
                end do
                
                ! check if divided dimension has been found
                if (div_dim.le.0) div_dim = size(vars(id)%tot_i_min)            ! assume that this is an individual print, so set div_dim to last index
            else                                                                ! equivalent with div_dim at last index
                div_dim = size(vars(id)%tot_i_min)
            end if
            
            ! set prev_dims
            do kd = 1,div_dim-1
                prev_dims = prev_dims*&
                    &(vars(id)%tot_i_max(kd)-vars(id)%tot_i_min(kd)+1)
            end do
            
            ! set next_dims
            do kd = div_dim+1,size(vars(id)%tot_i_min)
                next_dims = next_dims*&
                    &(vars(id)%tot_i_max(kd)-vars(id)%tot_i_min(kd)+1)
            end do
            
            ! set memory variables
            mem_block = prev_dims*&
                &(vars(id)%grp_i_max(div_dim)-vars(id)%grp_i_min(div_dim)+1)
            mem_offset = prev_dims*&
                &(vars(id)%grp_i_min(div_dim)-vars(id)%tot_i_min(div_dim))
            mem_count = 1
            mem_stride = 1
            
            ! set file variables
            file_block = prev_dims*&
                &(vars(id)%tot_i_max(div_dim)-vars(id)%tot_i_min(div_dim)+1)
        end subroutine calc_bocs
    end function print_HDF5_arrs
    
    ! reads a PB3D output file in HDF5 format
    integer function read_HDF5_arrs(vars,head_name) result(ierr)
        use num_vars, only: PB3D_name
        
        character(*), parameter :: rout_name = 'read_HDF5_arrs'
        
        ! input / output
        type(var_1D), intent(inout), allocatable :: vars(:)                     ! variables to write
        character(len=*), intent(in) :: head_name                               ! head name of variables
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: dset_id                                               ! dataset identifier 
        integer(HID_T) :: HDF5_kind_64                                          ! HDF5 type corresponding to dp
        integer(HID_T) :: group_id                                              ! group identifier
        integer(HID_T) :: head_group_id                                         ! head group identifier
        integer(HID_T) :: filespace                                             ! dataspace identifier in file 
        integer(HSIZE_T) :: id                                                  ! counter
        integer(HSIZE_T) :: data_size                                           ! size of data set
        integer(SIZE_T) :: name_len                                             ! length of name of group
        integer :: storage_type                                                 ! type of storage used in HDF5 file
        integer :: nr_lnks_head                                                 ! number of links in head group
        integer :: max_corder                                                   ! current maximum creation order value for group
        integer, allocatable :: lim_loc(:)                                      ! local copy of limits
        character(len=max_str_ln) :: name_len_loc                               ! local copy of name_len
        character(len=max_str_ln) :: group_name                                 ! name of group
        
        ! initialize ierr
        ierr = 0
        
        call writo('Reading data from PB3D output "'//trim(PB3D_name)//'/'//&
            &trim(head_name)//'"')
        call lvl_ud(1)
        
        ! initialize FORTRAN predefined datatypes
        call H5open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! preparation
        HDF5_kind_64 = H5kind_to_type(dp,H5_REAL_KIND)
        
        ! open the file
        call H5Fopen_f(trim(PB3D_name),H5F_ACC_RDONLY_F,HDF5_i,ierr)
        CHCKERR('Failed to open file')
        
        ! open head group
        call H5Gopen_f(HDF5_i,trim(head_name),head_group_id,ierr)
        CHCKERR('Failed to open group')
        
        ! get number of objects in group to allocate vars
        call H5Gget_info_f(head_group_id,storage_type,nr_lnks_head,max_corder,&
            &ierr)
        CHCKERR('Failed to get group info')
        allocate(vars(nr_lnks_head))
        
        ! iterate over all elements in group
        do id = 1, nr_lnks_head
            call H5Lget_name_by_idx_f(head_group_id,'.',H5_INDEX_NAME_F,&
                &H5_ITER_INC_F,id-1,group_name,ierr,size=name_len)
            CHCKERR('Failed to get name')
            
            ! check length
            if (name_len.gt.max_str_ln) then
                write(name_len_loc,*) name_len
                ierr = 1
                err_msg = 'Recompile with max_str_ln > '//trim(name_len_loc)
                CHCKERR(err_msg)
            end if
            
            ! user output and set name
            call writo('Reading '//trim(group_name))
            vars(id)%var_name = trim(group_name)
            
            ! open group
            call H5Gopen_f(head_group_id,trim(group_name),group_id,ierr)
            CHCKERR('Failed to open group')
            
            ! 1. variable
            
            ! open variable dataset
            call h5dopen_f(group_id,'var',dset_id,ierr)
            CHCKERR('Failed to open dataset')
            
            ! get dataspace
            call H5Dget_space_f(dset_id,filespace,ierr)
            CHCKERR('Failed to get file space')
            
            ! get size to allocate the 1D variable
            call H5Sget_simple_extent_npoints_f(filespace,data_size,ierr)
            CHCKERR('Failed to get storage size')
            allocate(vars(id)%p(data_size))
            
            ! read into 1D variable
            call H5Dread_f(dset_id,HDF5_kind_64,vars(id)%p,[data_size],ierr)
            CHCKERR('Failed to read dataset')
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! 2. limit information
            
            ! open limit information dataset
            call h5dopen_f(group_id,'lim',dset_id,ierr)
            CHCKERR('Failed to open dataset')
            
            ! get dataspace
            call H5Dget_space_f(dset_id,filespace,ierr)
            CHCKERR('Failed to get file space')
            
            ! get size to allocate the 1D variable
            call H5Sget_simple_extent_npoints_f(filespace,data_size,ierr)
            CHCKERR('Failed to get storage size')
            allocate(vars(id)%tot_i_min(data_size/2))
            allocate(vars(id)%tot_i_max(data_size/2))
            
            ! allocate local variables
            allocate(lim_loc(data_size))
            
            ! read into 1D variable
            call H5Dread_f(dset_id,H5T_NATIVE_INTEGER,lim_loc,[data_size],ierr)
            CHCKERR('Failed to read dataset')
            
            ! copy into tot_i_min and max
            vars(id)%tot_i_min = lim_loc(1:data_size/2)
            vars(id)%tot_i_max = lim_loc(data_size/2+1:data_size)
            
            ! close the dataset.
            call H5Dclose_f(dset_id,ierr)
            CHCKERR('Failed to close data set')
            
            ! close the group
            call H5gclose_f(group_id,ierr)
            CHCKERR('Failed to close group')
            
            ! deallocate local variables
            deallocate(lim_loc)
        end do
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        err_msg = 'Failed to close FORTRAN HDF5 interface'
        CHCKERR(err_msg)
        
        call lvl_ud(-1)
        call writo('Data read')
    end function read_HDF5_arrs
end module HDF5_ops
