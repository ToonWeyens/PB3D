!------------------------------------------------------------------------------!
!   Numerical utilities related to PB3D operations                             !
!------------------------------------------------------------------------------!
module PB3D_utilities
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_name_ln
    use HDF5_vars, only: var_1D_type

    implicit none
    private
    public conv_1D2ND, setup_rich_id, setup_eq_id, &
        &setup_par_id
    
    ! interfaces
    interface conv_1D2ND
        module procedure &
            &conv_1D2ND_1D_ind, conv_1D2ND_2D_ind, conv_1D2ND_3D_ind, &
            &conv_1D2ND_3D_arr, conv_1D2ND_4D_ind, conv_1D2ND_6D_ind, &
            &conv_1D2ND_7D_ind, conv_1D2ND_4D_arr, conv_1D2ND_6D_arr, &
            &conv_1D2ND_7D_arr
    end interface
    
contains
    ! Converts 1D to nD variables. The output variable has to be allocatable and
    ! unallocated.
    ! The  array versions,  which exist  for  the 3D,  4D, 6D  and 7D  versions,
    ! perform this for multiple indices in  position 1, as this generally agrees
    ! with  the  parallel index.  If  this  is  not  the desired  behavior,  the
    ! individual  versions have  to be  called separately.  Furthermore, in  the
    ! array routines, the total size of the first index has to be provided.
    ! Note that the  output variables are allocatables, which  has the important
    ! consequence that  lower and upper bounds  are passed as well.  This can be
    ! seen in the reconstruction routines "reconstruct_PB3D_..".
    subroutine conv_1D2ND_1D_ind(var_in,var_out)                                ! 1D individual version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:)                      ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1])
    end subroutine conv_1D2ND_1D_ind
    subroutine conv_1D2ND_2D_ind(var_in,var_out)                                ! 2D individual version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:)                    ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1])
    end subroutine conv_1D2ND_2D_ind
    subroutine conv_1D2ND_3D_ind(var_in,var_out)                                ! 3D individual version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:)                  ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1])
    end subroutine conv_1D2ND_3D_ind
    subroutine conv_1D2ND_3D_arr(var_in,dim_1,overlap,var_out)                  ! 3D array version
        use output_ops
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:)                              ! 2D variable
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:)                  ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:)                             ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_in)
            call conv_1D2ND_3D_ind(var_in(id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,&
                &lbound(var_out_loc,2):ubound(var_out_loc,2),&
                &lbound(var_out_loc,3):ubound(var_out_loc,3)))
            var_out(id_loc:id_loc+size(var_out_loc,1)-1,:,:) = var_out_loc
            id_loc = id_loc + size(var_out_loc,1)
            if (overlap) id_loc = id_loc - 1
            deallocate(var_out_loc)
        end do
    end subroutine conv_1D2ND_3D_arr
    subroutine conv_1D2ND_4D_ind(var_in,var_out)                                ! 4D individual version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:)                ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            &var_in%tot_i_max(4)-var_in%tot_i_min(4)+1])
    end subroutine conv_1D2ND_4D_ind
    subroutine conv_1D2ND_4D_arr(var_in,dim_1,overlap,var_out)                  ! 4D array version
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:)                              ! 2D variable
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:)                ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:,:)                           ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_in)
            call conv_1D2ND_4D_ind(var_in(id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,&
                &lbound(var_out_loc,2):ubound(var_out_loc,2),&
                &lbound(var_out_loc,3):ubound(var_out_loc,3),&
                &lbound(var_out_loc,4):ubound(var_out_loc,4)))
            var_out(id_loc:id_loc+size(var_out_loc,1)-1,:,:,:) = var_out_loc
            id_loc = id_loc + size(var_out_loc,1)
            if (overlap) id_loc = id_loc - 1
            deallocate(var_out_loc)
        end do
    end subroutine conv_1D2ND_4D_arr
    subroutine conv_1D2ND_6D_ind(var_in,var_out)                                ! 6D individual version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:)            ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4),&
            &var_in%tot_i_min(5):var_in%tot_i_max(5),&
            &var_in%tot_i_min(6):var_in%tot_i_max(6)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            &var_in%tot_i_max(4)-var_in%tot_i_min(4)+1,&
            &var_in%tot_i_max(5)-var_in%tot_i_min(5)+1,&
            &var_in%tot_i_max(6)-var_in%tot_i_min(6)+1])
    end subroutine conv_1D2ND_6D_ind
    subroutine conv_1D2ND_6D_arr(var_in,dim_1,overlap,var_out)                  ! 6D array version
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:)                              ! 2D variable
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:)            ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:,:,:,:)                       ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_in)
            call conv_1D2ND_6D_ind(var_in(id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,&
                &lbound(var_out_loc,2):ubound(var_out_loc,2),&
                &lbound(var_out_loc,3):ubound(var_out_loc,3),&
                &lbound(var_out_loc,4):ubound(var_out_loc,4),&
                &lbound(var_out_loc,5):ubound(var_out_loc,5),&
                &lbound(var_out_loc,6):ubound(var_out_loc,6)))
            var_out(id_loc:id_loc+size(var_out_loc,1)-1,:,:,:,:,:) = var_out_loc
            id_loc = id_loc + size(var_out_loc,1)
            if (overlap) id_loc = id_loc - 1
            deallocate(var_out_loc)
        end do
    end subroutine conv_1D2ND_6D_arr
    subroutine conv_1D2ND_7D_ind(var_in,var_out)                                ! 7D individual version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:,:)          ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1),&
            &var_in%tot_i_min(2):var_in%tot_i_max(2),&
            &var_in%tot_i_min(3):var_in%tot_i_max(3),&
            &var_in%tot_i_min(4):var_in%tot_i_max(4),&
            &var_in%tot_i_min(5):var_in%tot_i_max(5),&
            &var_in%tot_i_min(6):var_in%tot_i_max(6),&
            &var_in%tot_i_min(7):var_in%tot_i_max(7)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            &var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            &var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            &var_in%tot_i_max(4)-var_in%tot_i_min(4)+1,&
            &var_in%tot_i_max(5)-var_in%tot_i_min(5)+1,&
            &var_in%tot_i_max(6)-var_in%tot_i_min(6)+1,&
            &var_in%tot_i_max(7)-var_in%tot_i_min(7)+1])
    end subroutine conv_1D2ND_7D_ind
    subroutine conv_1D2ND_7D_arr(var_in,dim_1,overlap,var_out)                  ! 7D array version
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:)                              ! 2D variable
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:,:)          ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:,:,:,:,:)                     ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_in)
            call conv_1D2ND_7D_ind(var_in(id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,&
                &lbound(var_out_loc,2):ubound(var_out_loc,2),&
                &lbound(var_out_loc,3):ubound(var_out_loc,3),&
                &lbound(var_out_loc,4):ubound(var_out_loc,4),&
                &lbound(var_out_loc,5):ubound(var_out_loc,5),&
                &lbound(var_out_loc,6):ubound(var_out_loc,6),&
                &lbound(var_out_loc,7):ubound(var_out_loc,7)))
            var_out(id_loc:id_loc+size(var_out_loc,1)-1,:,:,:,:,:,:) = &
                &var_out_loc
            id_loc = id_loc + size(var_out_loc,1)
            if (overlap) id_loc = id_loc - 1
            deallocate(var_out_loc)
        end do
    end subroutine conv_1D2ND_7D_arr
    
    ! setup parallel id:
    !   par_id(1): start index
    !   par_id(2): end index
    !   par_id(3): stride
    function setup_par_id(grid,rich_lvl_max,rich_lvl_loc,tot_rich) &
        &result(par_id)
        use grid_vars, only: grid_type
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     ! grid
        integer, intent(in) :: rich_lvl_max                                     ! maximum Richardson level
        integer, intent(in) :: rich_lvl_loc                                     ! local Richardson level
        logical, intent(in), optional :: tot_rich                               ! whether to combine with previous Richardson levels
        integer :: par_id(3)                                                    ! parallel id
        
        ! local variables
        logical :: tot_rich_loc                                                 ! local tot_rich
        
        ! set up local tot_rich
        tot_rich_loc = .false.
        if (present(tot_rich) .and. rich_lvl_max.gt.1) tot_rich_loc = tot_rich  ! only for higher Richardson levels
        
        ! set up parallel indices
        if (tot_rich_loc) then
            if (rich_lvl_loc.eq.1) then
                par_id(1) = 1                                                   ! start at 1
                par_id(2) = grid%n(1)                                           ! end at n(1)
                par_id(3) = 2**(rich_lvl_max-1)                                 ! stride 2^(rich_lvl_max-1), same as for rich_lvl_loc = 2
            else
                par_id(1) = 1+2**(rich_lvl_max-rich_lvl_loc)                    ! start at 1+2^(rich_lvl_max-rich_lvl_loc)
                par_id(2) = grid%n(1)-2**(rich_lvl_max-rich_lvl_loc)            ! end at n(1)-2^(rich_lvl_max-rich_lvl_loc)
                par_id(3) = 2**(rich_lvl_max+1-rich_lvl_loc)                    ! stride 2^(rich_lvl_max+1-rich_lvl_loc)
            end if
        else
            par_id = [1,grid%n(1),1]                                            ! start at 1, end at n(1), stride 1
        end if
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
    
    ! setup equilibrium job id:
    !   eq_id(1): start equilibrium job (always 1)
    !   eq_id(2): end equilibrium job
    ! This  procedure   first  seeks  whether   a  specific  eq_job  is   to  be
    ! reconstructed. If not,  it checks whether there is a  variable without the
    ! equilibrium job suffix,  in which case the value 0  is returned for eq_id.
    ! If neither of the  previous two cases, it checks whether  there is a valid
    ! range for equilibrium job indices, which will be returned.
    ! If nothing found, negative numbers are returned.
    function setup_eq_id(group_name,eq_job,rich_lvl) result(eq_id)
        use num_vars, only: PB3D_name
        use HDF5_utilities, only: probe_HDF5_group
        
        ! input / output
        character(len=*), intent(in) :: group_name                              ! name of variable to find
        integer, intent(in), optional :: eq_job                                 ! equilibrium job to reconstruct
        integer, intent(in), optional :: rich_lvl                               ! Richardson level of variable
        integer :: eq_id(2)                                                     ! equilibrium id
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: istat                                                        ! status
        integer :: eq_job_loc                                                   ! local eq_job
        logical :: group_exists                                                 ! whether group exists
        character(len=max_str_ln) :: group_name_loc                             ! local group name
        
        ! set local eq_job
        eq_job_loc = 0
        if (present(eq_job)) eq_job_loc = eq_job
        
        ! set local group name
        group_name_loc = group_name
        
        ! possibly append Richardson level
        if (present(rich_lvl)) then
            if (rich_lvl.gt.0) then
                group_name_loc = trim(group_name_loc)//'_R_'//&
                    &trim(i2str(rich_lvl))
            end if
        end if
        
        ! 1. check for specific equilibrium job
        if (eq_job_loc.gt.0) then
            eq_id = eq_job_loc
            return
        end if
        
        ! 2. check whether there is a variable without eq_job suffix
        istat = probe_HDF5_group(PB3D_name,group_name_loc,group_exists)
        CHCKSTT
        if (group_exists) then
            eq_id = 0
            return
        end if
        
        ! 3. check whether there is a range of eq_job indices, starting with 1
        eq_id = [1,0]
        id = 1
        group_exists = .true.
        do while (group_exists)
            istat = probe_HDF5_group(PB3D_name,trim(group_name_loc)//'_E_'//&
                &trim(i2str(id)),group_exists)
            CHCKSTT
            if (group_exists) eq_id(2) = id
            id = id+1
        end do
        if (eq_id(2).gt.0) return
        
        ! nothing found: return negative number
        eq_id = -1
    end function setup_eq_id
end module PB3D_utilities
