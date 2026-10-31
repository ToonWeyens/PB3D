!------------------------------------------------------------------------------!
!> The 1-D flattened variable type used for internal HDF5 storage.
!!
!! This lives in its own dependency-light module (split off from HDF5_vars,
!! which re-exports it for compatibility) so that modules whose only HDF5
!! connection is this type - such as PB3D_utilities - can be compiled and
!! unit-tested without the HDF5 library.
!------------------------------------------------------------------------------!
module var_1D_vars
    use num_vars, only: dp, max_str_ln

    implicit none
    private
    public dealloc_var_1D, max_dim_var_1D

    ! global variables
    integer, parameter :: max_dim_var_1D = 100000                               !< maximum dimension of var_1D

    !> 1D  equivalent  of multidimensional  variables,  used  for internal  HDF5
    !! storage.
    type, public :: var_1D_type
        real(dp), allocatable :: p(:)                                           !< 1D equivalent of data of variable
        integer, allocatable :: tot_i_min(:)                                    !< total min.of indices of variable
        integer, allocatable :: tot_i_max(:)                                    !< total max.of indices of variable
        integer, allocatable :: loc_i_min(:)                                    !< group min.of indices of variable
        integer, allocatable :: loc_i_max(:)                                    !< group max.of indices of variable
        character(len=max_str_ln) :: var_name                                   !< name of variable
    end type var_1D_type

    ! interfaces

    !> \public Deallocates 1-D variables.
    interface dealloc_var_1D
        !> \public
        module procedure dealloc_var_1D_ind
        !> \public
        module procedure dealloc_var_1D_arr
        !> \public
        module procedure dealloc_var_1D_arr_2
    end interface

contains
    !> \private rank 2 array version
    subroutine dealloc_var_1D_arr_2(var_1D)
        ! input / output
        type(var_1D_type), intent(inout), allocatable :: var_1D(:,:)            !< array of 1D variables to be deallocated

        ! local variables
        integer :: id, jd                                                       ! counters

        ! deallocate individual arrays
        do jd = 1,size(var_1D,2)
            do id = 1,size(var_1D,1)
                call dealloc_var_1D_ind(var_1D(id,jd))
            end do
        end do

        ! deallocate the array
        deallocate(var_1D)
    end subroutine dealloc_var_1D_arr_2
    !> \private array version
    subroutine dealloc_var_1D_arr(var_1D)
        ! input / output
        type(var_1D_type), intent(inout), allocatable :: var_1D(:)              !< array of 1D variables to be deallocated

        ! local variables
        integer :: id                                                           ! counter

        ! deallocate individual arrays
        do id = 1,size(var_1D)
            call dealloc_var_1D_ind(var_1D(id))
        end do

        ! deallocate the array
        deallocate(var_1D)
    end subroutine dealloc_var_1D_arr
    !> \private individual version
    subroutine dealloc_var_1D_ind(var_1D)
        ! input / output
        type(var_1D_type), intent(out) :: var_1D                                !< 1D variable to be deallocated
    end subroutine dealloc_var_1D_ind
end module var_1D_vars
