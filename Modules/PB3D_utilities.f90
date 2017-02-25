!------------------------------------------------------------------------------!
!   Numerical utilities related to PB3D operations                             !
!------------------------------------------------------------------------------!
module PB3D_utilities
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_name_ln
    use HDF5_vars, only: var_1D_type

    implicit none
    private
    public conv_1D2ND, setup_rich_id, setup_par_id
    
    ! interfaces
    interface conv_1D2ND
        module procedure &
            &conv_1D2ND_1D, conv_1D2ND_2D, conv_1D2ND_3D, conv_1D2ND_4D, &
            &conv_1D2ND_6D, conv_1D2ND_7D
    end interface
    
contains
    ! Converts 1D to nD variables.
    subroutine conv_1D2ND_1D(var_in,var_out)                                    ! 1D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:)                      ! output variable
        
        ! allocate and copy variable
        if (allocated(var_out)) deallocate(var_out)
        allocate(var_out(var_in%tot_i_min(1):var_in%tot_i_max(1)))
        var_out = reshape(var_in%p,shape(var_out))
    end subroutine conv_1D2ND_1D
    subroutine conv_1D2ND_2D(var_in,var_out)                                    ! 2D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:)                    ! output variable
        
        ! allocate and copy variable
        if (allocated(var_out)) deallocate(var_out)
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1])
    end subroutine conv_1D2ND_2D
    subroutine conv_1D2ND_3D(var_in,var_out)                                    ! 3D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:)                  ! output variable
        
        ! allocate and copy variable
        if (allocated(var_out)) deallocate(var_out)
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3)))
        var_out = reshape(var_in%p,shape(var_out))
    end subroutine conv_1D2ND_3D
    subroutine conv_1D2ND_4D(var_in,var_out)                                    ! 4D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:)                ! output variable
        
        ! allocate and copy variable
        if (allocated(var_out)) deallocate(var_out)
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4)))
        var_out = reshape(var_in%p,shape(var_out))
    end subroutine conv_1D2ND_4D
    subroutine conv_1D2ND_6D(var_in,var_out)                                    ! 6D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:)            ! output variable
        
        ! allocate and copy variable
        if (allocated(var_out)) deallocate(var_out)
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4),&
            &var_in%tot_i_min(5):var_in%tot_i_max(5),&
            &var_in%tot_i_min(6):var_in%tot_i_max(6)))
        var_out = reshape(var_in%p,shape(var_out))
    end subroutine conv_1D2ND_6D
    subroutine conv_1D2ND_7D(var_in,var_out)                                    ! 7D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:,:)          ! output variable
        
        ! allocate and copy variable
        if (allocated(var_out)) deallocate(var_out)
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4),&
            &var_in%tot_i_min(5):var_in%tot_i_max(5),&
            &var_in%tot_i_min(6):var_in%tot_i_max(6),&
            &var_in%tot_i_min(7):var_in%tot_i_max(7)))
        var_out = reshape(var_in%p,shape(var_out))
    end subroutine conv_1D2ND_7D
    
    ! setup parallel id:
    !   par_id(1): start index
    !   par_id(2): end index
    !   par_id(3): stride
    ! These  are set  up by  considering the  transformation between  a point  r
    ! (1..n) in the local variable,  with corresponding indices (a..b) indicated
    ! by par_lim in the HDF5 variable in memory, so that n = b-a+1:
    !   p + k s = a + r - 1,
    ! where k is  an integer, s is the stride for Richardson  level i (with max.
    ! level I):
    !   s = 2^(I-1)     for i = 1
    !     = 2^(I-i+1)   for i > 1
    ! and where p is the starting point in the HDF5 variables:
    !   p = 1           for i = 1
    !     = 1 + s/2     for i > 1
    ! Therefore, the minimum possible r lies in (0..s-1):
    !   mod(r-1,s) = r_min-1,
    ! which leads to
    !   r_min = 1 + mod(p-a,s).
    ! The maximum possible r, then, lies in (n-s+1..n), which leads to:
    !   r_max = n - s + 1 + mod(p-b-1,s).
    ! These limits and the strides are saved in par_id.
    ! If the optional indices a and b are not given they are assumed to be 1 and
    ! n, with n = 1+ks, which simplifies the equations to:
    !   r_min = p
    !   r_max = n-p+1.
    ! Optionally, the indices in the HDF5 arrays can also be returned. These are
    ! equal to k-1, where k is the integer refered to above:
    !   par_id_mem = 1 + (a-1-p+par_id)/s.
    function setup_par_id(grid,rich_lvl_max,rich_lvl_loc,tot_rich,par_lim,&
        &par_id_mem) result(par_id)
        use grid_vars, only: grid_type
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        integer, intent(in) :: rich_lvl_max                                     ! maximum Richardson level
        integer, intent(in) :: rich_lvl_loc                                     ! local Richardson level
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer, intent(in), optional :: par_lim(2)                             ! limits (a..b) of variable requested
        integer, intent(inout), optional :: par_id_mem(2)                       ! parallel id in memory
        integer :: par_id(3)                                                    ! parallel id
        
        ! local variables
        integer :: s                                                            ! stride
        integer :: p                                                            ! lower limit in memory
        integer :: n                                                            ! b-a+1
        logical :: tot_rich_loc                                                 ! local tot_rich
        integer :: par_lim_loc(2)                                               ! local par_lim
        
        ! set up local tot_rich and par_lim
        tot_rich_loc = .false.
        if (present(tot_rich) .and. rich_lvl_max.gt.1) tot_rich_loc = tot_rich  ! only for higher Richardson levels
        par_lim_loc = [1,grid%n(1)]
        if (present(par_lim)) par_lim_loc = par_lim
        
        ! set up parallel indices
        n = par_lim_loc(2)-par_lim_loc(1)+1
        if (tot_rich_loc) then
            if (rich_lvl_loc.eq.1) then
                s = 2**(rich_lvl_max-1)
                p = 1
            else
                s = 2**(rich_lvl_max+1-rich_lvl_loc)
                p = 1 + s/2
            end if
        else
            s = 1
            p = 1
        end if
        par_id(1) = 1 + modulo(p-par_lim_loc(1),s)
        par_id(2) = n - s + 1 + modulo(p-par_lim_loc(2)-1,s)
        par_id(3) = s
        
        ! set par_id in memory if requested
        if (present(par_id_mem)) par_id_mem = &
            &1 + (par_lim_loc(1)-1-p+par_id(1:2))/s
    end function setup_par_id
    
    ! setup richardson id:
    !   rich_id(1): start Richardson level
    !   rich_id(2): end Richardson level
    function setup_rich_id(rich_lvl_max,tot_rich) result(rich_id)
        ! input / output
        integer, intent(in) :: rich_lvl_max                                     ! maximum Richardson level
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer :: rich_id(2)                                                   ! Richardson id
        
        ! set up rich_id
        rich_id = [rich_lvl_max,rich_lvl_max]
        if (present(tot_rich) .and. rich_lvl_max.gt.1) then                     ! only for higher Richardson levels
            if (tot_rich) rich_id = [1,rich_lvl_max]
        end if
    end function setup_rich_id
end module PB3D_utilities
