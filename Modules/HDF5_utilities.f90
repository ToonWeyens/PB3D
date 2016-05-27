!------------------------------------------------------------------------------!
!   Utilities pertaining to HDF5 and XDMF                                      !
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
module HDF5_utilities
#include <PB3D_macros.h>
    use num_vars, only: max_str_ln, dp, plot_dir, data_dir, script_dir
    use messages, only: writo, lvl_ud
    use str_ops, only: i2str, r2str, r2strt
    use HDF5_vars
    use HDF5
    
    implicit none
    private
    public probe_HDF5_group, set_1D_vars
#if ldebug
    public list_all_vars_in_group, &
        &debug_set_1D_vars
#endif
    
#if ldebug
    ! global variables
    logical :: debug_set_1D_vars = .false.                                      ! set to true to debug set_1D_vars
#endif
    
contains
    ! Sets up the  chunk property list and/or the 1D  hyperslabs that correspond
    ! to a local subset  of lim_tot given by the limits  lim_loc. An example for
    ! the hyperslab is given in 3 dimensions with size n = 5,3,3:
    !   for dim = 1, limits: 2..4
    !   for dim = 2, limits: 3..3
    !   for dim = 3, limits: 2..3.
    ! The 1D equivalent can be represented as:
    ! [.....][.....][.....] | [.....][.....][.....] | [.....][.....][.....]
    ! dimension 1 will only allow for the following elements (given by "-"):
    ! [.---.][.---.][.---.] | [.---.][.---.][.---.] | [.---.][.---.][.---.]
    ! dimension 2 will only allow for the following elements:
    ! [.....][.....][-----] | [.....][.....][-----] | [.....][.....][-----]
    ! dimension 3 will only allow for the following elements:
    ! [.....][.....][.....] | [-----][-----][-----] | [-----][-----][-----]
    ! so that the total is given by:
    ! [.....][.....][.....] | [.....][.....][.....] | [.....][.....][.---.]
    ! In practice, it is easiest to work  with a full selection, where for every
    ! dimension the ranges that do not correspond to the range of that dimension
    ! are taken away.  Clearly, if a dimension has the  full range, nothing will
    ! be taken out:
    !   block_i     = n_1*n_2*..*n_(i-1) * (b_i-a_i+1) ,
    !   stride_i    = n_1*n_2*..*n_(i-1) * (B_i-A_i+1) ,
    !   offset_i    = n_1*n_2*..*n_(i-1) * (a_i-A_i) ,
    !   count_i     = n_(i+1)*n_(i+2)*..*n_N ,
    ! where a_i and b_i  represent the local limits of dimension  i, A_i and B_i
    ! the total ones and the number of dimensions is N, as can be verified.
    ! The chunk property list for creation and  access can be set up so that its
    ! size is equal to the largest dimensions  of full range, or an integer part
    ! of that if it exceeds 4GB, and the number of elements in the hash table is
    ! set up optimally.  Furthermore, since the variables don't need  to be used
    ! more than once, the w0 factor is set to 1.
    integer function set_1D_vars(lim_tot,lim_loc,space_id,c_plist_id,&
        &a_plist_id) result(ierr)
        character(*), parameter :: rout_name = 'set_1D_vars'
        
        ! input / output
        integer, intent(in) :: lim_tot(:,:)                                     ! total limits
        integer, intent(in) :: lim_loc(:,:)                                     ! local limits
        integer(HID_T), intent(in), optional :: space_id                        ! dataspace identifier
        integer(HID_T), intent(inout), optional :: c_plist_id                   ! chunk creation property list identifier 
        integer(HID_T), intent(inout), optional :: a_plist_id                   ! chunk access property list identifier 
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: n_dims                                                       ! nr. of dimensions
        integer :: n_prod(2)                                                    ! n_prod(1) = n_1*n_2*..*n_(i-1), n_prod(2) = n_prod(1)*n_i
        integer :: div_dim                                                      ! first divided dimension
        integer :: chunk_size                                                   ! size of chunk
        integer(HSIZE_T) :: block(1)                                            ! block size in memory
        integer(HSIZE_T) :: stride(1)                                           ! stride in memory
        integer(HSIZE_T) :: count(1)                                            ! nr. of repetitions of block in memory
        integer(HSIZE_T) :: offset(1)                                           ! offset in memory
        integer(HSIZE_T) :: chunk_dims(1)                                       ! chunk dimensions
        integer(SIZE_T) :: chunk_nslots                                         ! number of chunk slots in the raw data chunk  cache hash table.
        integer(SIZE_T) :: chunk_nbytes                                         ! total size of the raw data chunk cache, in bytes. 
        real(dp) :: max_chunk_size = 4.E9_dp / sizeof(1._dp)                    ! maximum 4 GB (from manual)
        real :: chunk_w0                                                        ! preemption policy.
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set n_dims
        n_dims = size(lim_tot,1)

#if ldebug
        if (size(lim_loc,1).ne.n_dims) then
            ierr = 1
            err_msg = 'lim_tot and lim_loc are not compatible'
            CHCKERR(err_msg)
        end if
        if (size(lim_tot,2).ne.2 .or. size(lim_loc,2).ne.2) then
            ierr = 1
            err_msg = 'lim_tot and_lim_loc need to contain 2 columns'
            CHCKERR(err_msg)
        end if
        do id = 1,n_dims
            if (lim_tot(id,1).gt.lim_loc(id,1) .or. &
                &lim_tot(id,2).lt.lim_loc(id,2)) then
                ierr = 1
                err_msg = 'Total limits must comprise local ones'
                CHCKERR(err_msg)
            end if
        end do
        
        if (debug_set_1D_vars) then
            write(*,*) 'in set_1D_vars:'
            write(*,*) 'lim_tot_lo = ', lim_tot(:,1)
            write(*,*) 'lim_tot_hi = ', lim_tot(:,2)
            write(*,*) 'lim_loc_lo = ', lim_loc(:,1)
            write(*,*) 'lim_loc_hi = ', lim_loc(:,2)
        end if
#endif
        
        ! initialize divided dimension
        div_dim = n_dims+1
        
        ! set hyperslab to be everything if requested
        if (present(space_id)) then
            call H5Sselect_all_f(space_id,ierr) 
            CHCKERR('Failed to select hyperslab')
        end if
        
        ! loop over dimensions to possibly restrict hyperslab and set div_dim
        do id = 1,n_dims
            ! only restrict hyperslab if local range differ from total one
            if (lim_loc(id,1).gt.lim_tot(id,1) .or. &
                &lim_loc(id,2).lt.lim_tot(id,2)) then
                ! set divided dimension if first time
                if (div_dim.eq.n_dims+1) div_dim = id
                
                ! calculations if hyperslab restriction requested
                if (present(space_id)) then
                    ! set auxiliary variables
                    n_prod(1) = product(lim_tot(1:id-1,2)-lim_tot(1:id-1,1)+1)
                    n_prod(2) = n_prod(1)*(lim_tot(id,2)-lim_tot(id,1)+1)
                    
#if ldebug
                    if (debug_set_1D_vars) then
                        write(*,*) 'dimension', id, 'of', n_dims
                        write(*,*) 'n_prod = ', n_prod
                    end if
#endif
                    
                    block = (lim_loc(id,2)-lim_loc(id,1)+1) * n_prod(1)
                    stride = n_prod(2)
                    offset = (lim_loc(id,1)-lim_tot(id,1)) * n_prod(1)
                    count = 1
                    do kd = id+1,n_dims
                        count = count*(lim_tot(kd,2)-lim_tot(kd,1)+1)
                    end do
                    
#if ldebug
                    if (debug_set_1D_vars) then
                        write(*,*) 'block, offset, stride, count = ', block, &
                            &offset, stride, count
                    end if
#endif
                    
                    call H5Sselect_hyperslab_f(space_id,H5S_SELECT_AND_F,&
                        &offset,count,ierr,stride=stride,block=block)
                    CHCKERR('Failed to select hyperslab')
                end if
            end if
        end do
        
        ! set chunk property list if requested
        if (present(c_plist_id) .or. present(a_plist_id)) then
            ! set up variables
            chunk_size = &
                &product(lim_tot(1:div_dim-1,2)-lim_tot(1:div_dim-1,1)+1)       ! all dimensions up to divided one
            chunk_size = chunk_size / (1 + floor(chunk_size/max_chunk_size))    ! limit to max chunk size
            
            ! set creation property list if provided
            if (present(c_plist_id)) then
                chunk_dims = chunk_size
                call H5Pcreate_f(H5P_DATASET_CREATE_F,c_plist_id,ierr) 
                CHCKERR('Failed to create property list')
                call H5Pset_chunk_f(c_plist_id,1,chunk_dims,ierr)               ! one dimension for 1D chunk
                CHCKERR('Failed to create chunk')
            end if
            
            ! set access property list if provided
            if (present(a_plist_id)) then
                chunk_nbytes = chunk_size * sizeof(1._dp)                       ! each chunk has a size in memory
                chunk_nslots = &
                    &ceiling(1._dp*product(lim_tot(:,2)-lim_tot(:,1)+1) /&
                    &chunk_size)                                                ! divide total number of elements by elements in chunk
                chunk_w0 = 1.0                                                  ! cache is never reused
                call H5Pcreate_f(H5P_DATASET_ACCESS_F,a_plist_id,ierr)
                CHCKERR('Failed to create property list')
                call H5Pset_chunk_cache_f(a_plist_id,chunk_nslots,&
                    &chunk_nbytes,chunk_w0,ierr)
                CHCKERR('Failed to set chunk cache')
            end if
        end if
    end function set_1D_vars
    
    ! Probe HDF5 file for group existence.
    integer function probe_HDF5_group(HDF5_name,group_name,group_exists) &
        &result(ierr)
        character(*), parameter :: rout_name = 'probe_HDF5_group'
        
        ! input / output
        character(len=*), intent(in) :: HDF5_name                               ! name of HDF5 file
        character(len=*), intent(in) :: group_name                              ! name of group to probe for
        logical, intent(inout) :: group_exists                                  ! whether group exists
        
        ! local variables
        integer :: istat                                                        ! status
        integer(HID_T) :: HDF5_i                                                ! file identifier 
        integer(HID_T) :: group_id                                              ! group identifier
        
        ! initialize ierr
        ierr = 0
        
        ! initialize FORTRAN predefined datatypes
        call H5open_f(ierr) 
        CHCKERR('Failed to initialize HDF5')
        
        ! open the file
        call H5Fopen_f(HDF5_name,H5F_ACC_RDONLY_F,HDF5_i,ierr)
        CHCKERR('Failed to open file. Does it exist?')
        
        ! disable error messages
        call h5eset_auto_f(0,ierr)
        CHCKERR('Failed to disable error printing')
        
        ! try to open group
        call H5Gopen_f(HDF5_i,group_name,group_id,istat)
        group_exists = istat.eq.0
        
        ! reenable error messages
        call h5eset_auto_f(1,ierr)
        CHCKERR('Failed to enable error printing')
        
        ! close group
        if (group_exists) then
            call H5gclose_f(group_id,ierr)
            CHCKERR('Failed to close head group')
        end if
        
        ! close the HDF5 file
        call H5Fclose_f(HDF5_i,ierr)
        CHCKERR('failed to close HDF5 file')
        
        ! close FORTRAN interfaces and HDF5 library.
        call H5Close_f(ierr)
        CHCKERR('Failed to close FORTRAN HDF5 interface')
    end function probe_HDF5_group

#if ldebug
    integer function list_all_vars_in_group(group_id) result(ierr)
        character(*), parameter :: rout_name = 'list_all_vars_in_group'
        
        ! input / output
        integer(HID_T), intent(in) :: group_id                                  ! group identifier
        
        ! local variables
        integer :: storage_type                                                 ! type of storage used in HDF5 file
        integer :: nr_lnks                                                      ! number of links in group
        integer :: max_corder                                                   ! current maximum creation order value for group
        integer(HSIZE_T) :: id                                                  ! counter
        integer :: id_loc                                                       ! local id
        integer(SIZE_T) :: name_len                                             ! length of name of group
        character(len=max_str_ln) :: var_name                                   ! name of variable
        
        ! get number of objects in group
        call H5Gget_info_f(group_id,storage_type,nr_lnks,max_corder,&
            &ierr)
        CHCKERR('Failed to get group info')
        
        ! user output
        call writo('The group has '//trim(i2str(nr_lnks))//' elements')
        
        ! loop over all elements
        do id = 1,nr_lnks
            ! get variable name
            call H5Lget_name_by_idx_f(group_id,'.',H5_INDEX_NAME_F,&
                &H5_ITER_NATIVE_F,id-1,var_name,ierr,size=name_len)
            CHCKERR('Failed to get name')
            
            ! user output
            id_loc = int(id,1)
            call writo('Element '//trim(i2str(id_loc))//'/'//&
                &trim(i2str(nr_lnks))//': '//trim(var_name))
        end do
    end function list_all_vars_in_group
#endif
end module HDF5_utilities
