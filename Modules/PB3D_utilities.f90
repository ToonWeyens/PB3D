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
    public get_full_var_names, retrieve_var_1D_id, conv_1D2ND, setup_rich_id, &
        &setup_eq_id, setup_par_id
    
    ! interfaces
    interface get_full_var_names
        module procedure get_full_var_names_1, get_full_var_names_2
    end interface
    interface conv_1D2ND
        module procedure &
            &conv_1D2ND_1D_ind, conv_1D2ND_2D_ind, conv_1D2ND_3D_ind, &
            &conv_1D2ND_3D_arr, conv_1D2ND_4D_ind, conv_1D2ND_6D_ind, &
            &conv_1D2ND_7D_ind, conv_1D2ND_4D_arr, conv_1D2ND_6D_arr, &
            &conv_1D2ND_7D_arr
    end interface
    interface retrieve_var_1D_id
        module procedure retrieve_var_1D_id_ind, retrieve_var_1D_id_arr
    end interface
    
contains
    ! Set  all possible  full  variable names  for given  input  names and  mode
    ! numbers.
    ! Note that  due to  strange FORTRAN behavior,  "lim_sec_X" is  not optional
    ! here: If it  is defined so, for  some reason the allocatable  array is not
    ! passed correctly. This is the only place where consistency is broken.
    subroutine get_full_var_names_1(var_names,full_var_names,lim_sec_X)         ! vectorial version
        use X_utilities, only: get_suffix
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            ! internal variable names
        character(len=max_name_ln), intent(inout), allocatable :: &
            &full_var_names(:)                                                  ! full HDF5 variable names
        integer, intent(in) :: lim_sec_X(2)                                     ! limits for m_X (pol flux) or n_X (tor flux)
        
        ! local variables
        integer :: n_vars                                                       ! nr. of input variables
        integer :: n_mod                                                        ! nr. of secondary modes
        integer :: id, jd                                                       ! counters
        integer :: id_loc                                                       ! counter
        integer :: suffix                                                       ! suffix
        
        ! set n_vars and n_mod
        n_vars = size(var_names)
        n_mod = lim_sec_X(2) - lim_sec_X(1) + 1
        
        ! allocate full HDF5 variable names
        if (allocated(full_var_names)) deallocate(full_var_names)
        allocate(full_var_names(n_vars*n_mod))
        
        ! loop over all modes
        id_loc = 1
        ! loop over all variables
        do id = 1,n_vars
            do jd = 1,n_mod
                suffix = get_suffix(jd,lim_sec_X)
                full_var_names(id_loc) = trim(var_names(id))//'_'//&
                    &trim(i2str(suffix))
                id_loc = id_loc+1
            end do
        end do
    end subroutine get_full_var_names_1
    subroutine get_full_var_names_2(var_names,sym,full_var_names,lim_sec_X)     ! tensorial version
        use X_vars, only: set_nn_mod
        use X_utilities, only: get_suffix, is_necessary_X
        
        ! input / output
        character(len=*), intent(in) :: var_names(:)                            ! internal variable names
        logical, intent(in) :: sym(:)                                           ! if variable is symmetric
        character(len=max_name_ln), intent(inout), allocatable :: &
            &full_var_names(:)                                                  ! full HDF5 variable names
        integer, intent(in) :: lim_sec_X(2,2)                                   ! limits for m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: n_vars                                                       ! nr. of input variables
        integer :: n_mod(2)                                                     ! nr. of secondary modes
        integer :: id, jd, kd                                                   ! counters
        integer :: id_loc                                                       ! counter
        integer :: nn_mod_tot                                                   ! total nn_mod
        integer :: suffix(2)                                                    ! suffix
        
        ! tests
        if (size(var_names).ne.size(sym)) then
            call writo('size of sym has to be equal to size of var_names',&
                &warning=.true.)
            return
        end if
        
        ! set n_vars and n_mod
        n_vars = size(var_names)
        n_mod = lim_sec_X(2,:) - lim_sec_X(1,:) + 1
        
        ! allocate full HDF5 variable names
        if (allocated(full_var_names)) deallocate(full_var_names)
        nn_mod_tot = 0
        do id = 1,n_vars
            if (sym(id)) then
                nn_mod_tot = nn_mod_tot + set_nn_mod(lim_sec_X)
            else
                nn_mod_tot = nn_mod_tot + product(n_mod)
            end if
        end do
        allocate(full_var_names(nn_mod_tot))
        
        ! loop over all modes
        id_loc = 1
        ! loop over all variables
        do id = 1,n_vars
            do kd = 1,n_mod(2)
                do jd = 1,n_mod(1)
                    suffix = get_suffix(jd,kd,lim_sec_X)
                    if (is_necessary_X(sym(id),[jd,kd],lim_sec_X)) then
                        full_var_names(id_loc) = trim(var_names(id))//'_'//&
                            &trim(i2str(suffix(1)))//'_'//trim(i2str(suffix(2)))
                        id_loc = id_loc+1
                    end if
                end do
            end do
        end do
    end subroutine get_full_var_names_2
    
    ! Retrieves variable index from array 1D equivalents
    integer function retrieve_var_1D_id_ind(vars,var_name,var_id) result(ierr)  ! individual version
        character(*), parameter :: rout_name = 'retrieve_var_1D_ind'
        
        ! input / output
        type(var_1D_type), intent(in) :: vars(:)                                ! array of 1D variables
        character(len=*), intent(in) :: var_name                                ! name of variable to retrieve
        integer, intent(inout) :: var_id                                        ! index of variable
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! look up the variable
        do id = 1,size(vars)
            if (trim(vars(id)%var_name).eq.trim(var_name)) then                 ! found
                var_id = id
                return
            end if
        end do
        
        ! if still here, nothing found
        ierr = 1
        err_msg = 'Variable '//trim(var_name)//' not found'
        CHCKERR(err_msg)
    end function retrieve_var_1D_id_ind
    integer function retrieve_var_1D_id_arr(vars,var_name,var_id) result(ierr)  ! array version
        character(*), parameter :: rout_name = 'retrieve_var_1D_arr'
        
        ! input / output
        type(var_1D_type), intent(in) :: vars(:,:)                              ! array of 1D variables
        character(len=*), intent(in) :: var_name                                ! name of variable to retrieve
        integer, intent(inout), allocatable :: var_id(:)                        ! index of variable
        
        ! local variables
        integer :: id                                                           ! counter
        
        ! allocate and test
        if (allocated(var_id)) deallocate(var_id)
        allocate(var_id(size(vars,2)))
        
        ! call individual version
        do id = 1,size(vars,2)
            ierr = retrieve_var_1D_id_ind(vars(:,id),var_name,var_id(id))
            CHCKERR('')
        end do
    end function retrieve_var_1D_id_arr
    
    ! Converts 1D to nD variables. The output variable has to be allocatable and
    ! unallocated.
    ! The  array versions,  which exist  for  the 3D,  4D, 6D  and 7D  versions,
    ! perform this for multiple indices in  position 1, as this generally agrees
    ! with  the  parallel index.  If  this  is  not  the desired  behavior,  the
    ! individual  versions have  to be  called separately.  Furthermore, in  the
    ! array routines, the total size of the first index has to be provided.
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
    subroutine conv_1D2ND_3D_arr(var_in,var_1D_id,dim_1,overlap,var_out)        ! 3D array version
        use output_ops
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:,:)                            ! 2D variable
        integer, intent(in) :: var_1D_id(:)                                     ! indices in var_in
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:)                  ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:)                             ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_1D_id)
            call conv_1D2ND_3D_ind(var_in(var_1D_id(id),id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,size(var_out_loc,2),&
                &size(var_out_loc,3)))
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
    subroutine conv_1D2ND_4D_arr(var_in,var_1D_id,dim_1,overlap,var_out)        ! 4D array version
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:,:)                            ! 2D variable
        integer, intent(in) :: var_1D_id(:)                                     ! indices in var_in
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:)                ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:,:)                           ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_1D_id)
            call conv_1D2ND_4D_ind(var_in(var_1D_id(id),id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,size(var_out_loc,2),&
                &size(var_out_loc,3),size(var_out_loc,4)))
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
    subroutine conv_1D2ND_6D_arr(var_in,var_1D_id,dim_1,overlap,var_out)        ! 6D array version
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:,:)                            ! 2D variable
        integer, intent(in) :: var_1D_id(:)                                     ! indices in var_in
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:)            ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:,:,:,:)                       ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_1D_id)
            call conv_1D2ND_6D_ind(var_in(var_1D_id(id),id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,size(var_out_loc,2),&
                &size(var_out_loc,3),size(var_out_loc,4),size(var_out_loc,5),&
                &size(var_out_loc,6)))
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
    subroutine conv_1D2ND_7D_arr(var_in,var_1D_id,dim_1,overlap,var_out)        ! 7D array version
        ! input / output
        type(var_1D_type), intent(in) :: var_in(:,:)                            ! 2D variable
        integer, intent(in) :: var_1D_id(:)                                     ! indices in var_in
        integer :: dim_1                                                        ! size of dimension 1
        logical :: overlap                                                      ! whether overlap of width 1 for equilibrium jobs
        real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:,:,:)          ! output variable
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: id_loc                                                       ! local id
        real(dp), allocatable :: var_out_loc(:,:,:,:,:,:,:)                     ! local var_out
        
        ! loop over all indices
        id_loc = 1
        do id = 1,size(var_1D_id)
            call conv_1D2ND_7D_ind(var_in(var_1D_id(id),id),var_out_loc)
            if (id.eq.1) allocate(var_out(dim_1,size(var_out_loc,2),&
                &size(var_out_loc,3),size(var_out_loc,4),size(var_out_loc,5),&
                &size(var_out_loc,6),size(var_out_loc,7)))
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
    !   eq_id(1): start equilibirum job (always 1)
    !   eq_id(2): end equilibrium job
    ! This  procedure   first  seeks  whether   a  specific  eq_job  is   to  be
    ! reconstructed. If not,  it checks whether there is a  variable without the
    ! equilibrium job suffix,  in which case the value 0  is returned for eq_id.
    ! If neither of the  previous two cases, it checks whether  there is a valid
    ! range for equilibrium job indices, which will be returned.
    ! If nothing found, negative numbers are returned.
    function setup_eq_id(group_name,eq_job,rich_lvl) result(eq_id)
        use num_vars, only: PB3D_name
        use HDF5_ops, only: probe_HDF5_group
        
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
