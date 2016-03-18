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
    public get_full_var_names, retrieve_var_1D_id, conv_1D2ND
    
    ! interfaces
    interface get_full_var_names
        module procedure get_full_var_names_1, get_full_var_names_2
    end interface
    interface conv_1D2ND
        module procedure conv_1D2ND_1D, conv_1D2ND_2D, conv_1D2ND_3D, &
            &conv_1D2ND_4D, &
            !&conv_1D2ND_5D, &
            &conv_1D2ND_6D, conv_1D2ND_7D
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
            call writo('WARNING: size of sym has to be equal to size of &
                &var_names')
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
    integer function retrieve_var_1D_id(vars,var_name,var_id) result(ierr)
        character(*), parameter :: rout_name = 'retrieve_var_1D'
        
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
    end function retrieve_var_1D_id
    
    ! Converts 1D to nD variables. The output variable has to be allocatable and
    ! unallocated
    subroutine conv_1D2ND_1D(var_in,var_out)                                    ! 1D version
        ! input / output
        type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        real(dp), intent(inout), allocatable :: var_out(:)                      ! output variable
        
        ! allocate and copy variable
        allocate(var_out(&
            &var_in%tot_i_min(1):var_in%tot_i_max(1)))
        var_out = reshape(var_in%p,[&
            &var_in%tot_i_max(1)-var_in%tot_i_min(1)+1])
    end subroutine conv_1D2ND_1D
    subroutine conv_1D2ND_2D(var_in,var_out)                                    ! 2D version
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
    end subroutine conv_1D2ND_2D
    subroutine conv_1D2ND_3D(var_in,var_out)                                    ! 3D version
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
    end subroutine conv_1D2ND_3D
    subroutine conv_1D2ND_4D(var_in,var_out)                                    ! 4D version
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
    end subroutine conv_1D2ND_4D
    !subroutine conv_1D2ND_5D(var_in,var_out)                                    ! 5D version
        !! input / output
        !type(var_1D_type), intent(in) :: var_in                                 ! 1D variable
        !real(dp), intent(inout), allocatable :: var_out(:,:,:,:,:)              ! output variable
        
        !! allocate and copy variable
        !allocate(var_out(&
            !&var_in%tot_i_min(1):var_in%tot_i_max(1),&
            !&var_in%tot_i_min(2):var_in%tot_i_max(2),&
            !&var_in%tot_i_min(3):var_in%tot_i_max(3),&
            !&var_in%tot_i_min(4):var_in%tot_i_max(4),&
            !&var_in%tot_i_min(5):var_in%tot_i_max(5)))
        !var_out = reshape(var_in%p,[&
            !&var_in%tot_i_max(1)-var_in%tot_i_min(1)+1,&
            !&var_in%tot_i_max(2)-var_in%tot_i_min(2)+1,&
            !&var_in%tot_i_max(3)-var_in%tot_i_min(3)+1,&
            !&var_in%tot_i_max(4)-var_in%tot_i_min(4)+1,&
            !&var_in%tot_i_max(5)-var_in%tot_i_min(5)+1])
    !end subroutine conv_1D2ND_5D
    subroutine conv_1D2ND_6D(var_in,var_out)                                    ! 6D version
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
    end subroutine conv_1D2ND_6D
    subroutine conv_1D2ND_7D(var_in,var_out)                                    ! 7D version
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
    end subroutine conv_1D2ND_7D
end module PB3D_utilities
