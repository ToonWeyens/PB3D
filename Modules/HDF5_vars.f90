!------------------------------------------------------------------------------!
!   Variables pertaining to HDF5 and XDMF                                      !
!------------------------------------------------------------------------------!
module HDF5_vars
#include <PB3D_macros.h>
    use num_vars, only: max_str_ln, dp, plot_dir, data_dir, script_dir
    use messages, only: writo, lvl_ud
    use str_ops, only: i2str, r2str, r2strt
    use HDF5
    
    implicit none
    private
    public init_HDF5, dealloc_var_1D, &
        &XML_str_type, HDF5_file_type, var_1D_type, xmf_fmt, XDMF_num_types, &
        &XDMF_format_types, XDMF_geom_types, XDMF_top_types, XDMF_att_types, &
        &XDMF_center_types, XDMF_grid_types, max_dim_var_1D
    
    ! global variables
    integer, parameter :: max_xml_ln = 300                                      ! max. length of xml string
    character(len=6) :: xmf_fmt = '(999A)'                                      ! format to write the xmf file
    
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
    type var_1D_type
        real(dp), allocatable :: p(:)                                           ! 1D equivalent of data of variable
        integer, allocatable :: tot_i_min(:), tot_i_max(:)                      ! total min. and max. of indices of variable
        integer, allocatable :: loc_i_min(:), loc_i_max(:)                      ! group min. and max. of indices of variable
        character(len=max_str_ln) :: var_name                                   ! name of variable
    end type var_1D_type
    
    ! XDMF possibilities
    character(len=max_str_ln) :: XDMF_num_types(2)                              ! possible XDMF number types
    character(len=max_str_ln) :: XDMF_format_types(2)                           ! possible XDMF format types
    character(len=max_str_ln) :: XDMF_geom_types(2)                             ! possible XDMF geometry types
    character(len=max_str_ln) :: XDMF_top_types(2)                              ! possible XDMF topology types
    character(len=max_str_ln) :: XDMF_att_types(1)                              ! possible XDMF attribute types
    character(len=max_str_ln) :: XDMF_center_types(2)                           ! possible XDMF attribute center types
    character(len=max_str_ln) :: XDMF_grid_types(3)                             ! possible XDMF grid types
    
    ! global variables
    integer, parameter :: max_dim_var_1D = 1000                                 ! maximum dimension of var_1D
    
    ! interfaces
    interface dealloc_var_1D
        module procedure dealloc_var_1D_ind, dealloc_var_1D_arr
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
    
    ! deallocates 1D variable
    subroutine dealloc_var_1D_arr(var_1D)                                       ! array version
        ! input / output
        type(var_1D_type), intent(out), allocatable :: var_1D(:)                ! arra of 1D variables to be deallocated
    end subroutine dealloc_var_1D_arr
    subroutine dealloc_var_1D_ind(var_1D)                                       ! individual version
        ! input / output
        type(var_1D_type), intent(out) :: var_1D                                ! 1D variable to be deallocated
    end subroutine dealloc_var_1D_ind
end module HDF5_vars
