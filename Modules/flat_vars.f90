!------------------------------------------------------------------------------!
!> The flattened variable type used for internal HDF5 storage.
!!
!! It has its own module, without HDF5, so that modules whose only HDF5
!! connection is this type, such as PB3D_utilities, can be compiled and
!! unit-tested without the HDF5 library.
!------------------------------------------------------------------------------!
module flat_vars
    use num_vars, only: dp, max_str_ln

    implicit none
    private
    public dealloc_var_flat, max_dim_var_flat

    ! global variables
    integer, parameter :: max_dim_var_flat = 100000                             !< maximum dimension of var_flat

    !> 1D  equivalent  of multidimensional  variables,  used  for internal  HDF5
    !! storage.
    type, public :: var_flat_type
        real(dp), allocatable :: p(:)                                           !< 1D equivalent of data of variable
        integer, allocatable :: tot_i_min(:)                                    !< total min.of indices of variable
        integer, allocatable :: tot_i_max(:)                                    !< total max.of indices of variable
        integer, allocatable :: loc_i_min(:)                                    !< group min.of indices of variable
        integer, allocatable :: loc_i_max(:)                                    !< group max.of indices of variable
        character(len=max_str_ln) :: var_name                                   !< name of variable
    end type var_flat_type

    ! interfaces

    !> \public Deallocates 1-D variables.
    interface dealloc_var_flat
        !> \public
        module procedure dealloc_var_flat_ind
        !> \public
        module procedure dealloc_var_flat_arr
        !> \public
        module procedure dealloc_var_flat_arr_2
    end interface

contains
    !> \private rank 2 array version
    subroutine dealloc_var_flat_arr_2(var_flat)
        ! input / output
        type(var_flat_type), intent(inout), allocatable :: var_flat(:,:)        !< array of 1D variables to be deallocated

        ! local variables
        integer :: id, jd                                                       ! counters

        ! deallocate individual arrays
        do jd = 1,size(var_flat,2)
            do id = 1,size(var_flat,1)
                call dealloc_var_flat_ind(var_flat(id,jd))
            end do
        end do

        ! deallocate the array
        deallocate(var_flat)
    end subroutine dealloc_var_flat_arr_2
    !> \private array version
    subroutine dealloc_var_flat_arr(var_flat)
        ! input / output
        type(var_flat_type), intent(inout), allocatable :: var_flat(:)          !< array of 1D variables to be deallocated

        ! local variables
        integer :: id                                                           ! counter

        ! deallocate individual arrays
        do id = 1,size(var_flat)
            call dealloc_var_flat_ind(var_flat(id))
        end do

        ! deallocate the array
        deallocate(var_flat)
    end subroutine dealloc_var_flat_arr
    !> \private individual version
    subroutine dealloc_var_flat_ind(var_flat)
        ! input / output
        type(var_flat_type), intent(out) :: var_flat                            !< 1D variable to be deallocated
    end subroutine dealloc_var_flat_ind
end module flat_vars
