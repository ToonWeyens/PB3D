program test
#define CHCKERR(s) if(ierr.ne.0) then; write(*,*) s; end if
    use HDF5
    use ISO_FORTRAN_ENV
    
    implicit none
    
    ! variables
    integer :: ierr
    integer, parameter :: dp = REAL64
    integer, parameter :: max_str_ln = 100
    character(len=max_str_ln) :: err_msg
    character(len=max_str_ln) :: file_name
    character(len=max_str_ln) :: group_name
    real(dp), allocatable :: buf(:)
    integer :: dims
    integer(HID_T) :: HDF5_kind_64
    integer(HID_T) :: HDF5_i
    integer(HID_T) :: group_id
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: mem_offset(1)
    integer(HSIZE_T) :: mem_stride(1)
    integer(HSIZE_T) :: mem_count(1)
    integer(HSIZE_T) :: mem_block(1)
    
    write(*,*) 'TESTING WRITING AND READING OF HDF5'
    
    ! set buffer
    dims = 10000
    allocate(buf(dims))
    buf = 1._dp
    
    ! set memory information
    mem_offset = 0
    mem_count = 1
    mem_stride = 1
    mem_block = dims
    
    ! initialize FORTRAN predefined datatypes
    call H5open_f(ierr) 
    CHCKERR('Failed to initialize HDF5')
    
    ! set HDF5 type corresponding to dp
    HDF5_kind_64 = H5kind_to_type(dp,H5_REAL_KIND)
    
    ! create the file
    file_name = 'testfile.h5'
    call H5Fcreate_f(file_name,H5F_ACC_TRUNC_F,HDF5_i,ierr)
    CHCKERR('Failed to create file')
    
    ! create group
    group_name = 'group'
    call H5Gcreate_f(HDF5_i,trim(group_name),group_id,ierr)
    CHCKERR('Failed to create group')
    
    ! create file data space
    call H5Screate_simple_f(1,mem_block,filespace,ierr)
    CHCKERR('Failed to create file space')
    
    ! create file data set in group
    call H5Dcreate_f(group_id,'var',HDF5_kind_64,filespace,dset_id,&
        &ierr)
    CHCKERR('Failed to create file data set')
    
    ! select hyperslab in the file
    call H5Sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,mem_offset,&
        &mem_count,ierr,mem_stride,mem_block)
    CHCKERR('Failed to select hyperslab')
    
    ! create memory data space
    call H5Screate_simple_f(1,mem_block,memspace,ierr)
    CHCKERR('Failed to create memory space')
    
    ! write the dataset collectively. 
    call H5Dwrite_f(dset_id,HDF5_kind_64,buf,mem_block,ierr,&
        &file_space_id=filespace,mem_space_id=memspace)
    CHCKERR('Failed to write data data set')
    
    ! close dataspaces.
    call H5Sclose_f(filespace,ierr)
    CHCKERR('Unable to close file space')
    call H5Sclose_f(memspace,ierr)
    CHCKERR('Unable to close memory space')
    
    ! close the dataset.
    call H5Dclose_f(dset_id,ierr)
    CHCKERR('Failed to close data set')
    
    ! close the group
    call H5gclose_f(group_id,ierr)
    CHCKERR('Failed to close group')
    
    ! close the HDF5 file
    call H5Fclose_f(HDF5_i,ierr)
    CHCKERR('failed to close HDF5 file')
    
    ! close FORTRAN interfaces and HDF5 library.
    call H5Close_f(ierr)
    err_msg = 'Failed to close FORTRAN HDF5 interface'
    CHCKERR(err_msg)
    
    write(*,*) 'TESTS DONE'
end program
