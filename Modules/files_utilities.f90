!------------------------------------------------------------------------------!
!   Numerical utilities related to files                                       !
!------------------------------------------------------------------------------!
module files_utilities
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public search_file, nextunit, wait_file, get_file_info

contains
    ! looks for the full name of a file and tests for its existence
    ! output:   full name of file, empty string if non-existent
    subroutine search_file(i_unit,file_name)
        character(len=*), intent(inout) :: file_name                            ! the name that is searched for
        integer, intent(out) :: i_unit                                          ! will hold the file handle
        
        character(len=max_str_ln) :: mod_file_name                               ! modified file name
        integer :: istat
        
        ! try to open the given name
        mod_file_name = file_name                                               ! copy file_name to mod_file_name
        open(unit=nextunit(i_unit),file=file_name,status='old',iostat=istat)
        if (istat.eq.0) return
        
        ! no matches, empty file_name and i_unit = 0 returned
        file_name = ""
        i_unit = 0
    end subroutine search_file
    
    ! Search for available  new unit where lun_min and lun_max  define the range
    ! of possible luns to check. The unit value is returned by the function, and
    ! also  by the  optional  argument.  This allows  the  function  to be  used
    ! directly in an  open statement, and optionally save the  result in a local
    ! variable. If no units are available, -1 is returned.
    ! (Adapted from: "http://fortranwiki.org/fortran/show/newunit")
    integer function nextunit(unit)
        ! input / output
        integer, intent(out), optional :: unit
        
        ! local variables
        integer, parameter :: lun_min=10, lun_max=1000
        logical :: opened
        integer :: lun
        
        ! iterate over permitted luns until available one found
        nextunit=-1
        do lun=lun_min,lun_max
            inquire(unit=lun,opened=opened)
            if (.not. opened) then
                nextunit=lun
                exit
            end if
        end do
        
        ! return unit number if present
        if (present(unit)) unit=nextunit
    end function nextunit
    
    ! Waits for file acces through a lock file.
    subroutine wait_file(lock_file_i,lock_file_name)
        ! input / output
        integer, intent(inout) :: lock_file_i                                   ! lock file nr.
        character(len=*), intent(in) :: lock_file_name                          ! name of lock file
        
        ! local variables
        logical :: file_exists                                                  ! file exists status
        integer :: open_stat                                                    ! file open status
        
        ! open X_jobs file if no lock-file exists
        file_exists = .true.
        do while (file_exists)
            open(STATUS='NEW',unit=nextunit(lock_file_i),&
                &file=lock_file_name,iostat=open_stat)
            if (open_stat.eq.0) then
                file_exists = .false.
            !else
                !call sleep(1)
            end if
        end do
    end subroutine wait_file
    
    ! Gets file information
    ! The  time informations  can be  converted to  strings using  the intrinsic
    ! function "ctime".
    subroutine get_file_info(file_name,file_size,acc_time,mod_time)
        ! input / output
        character(len=*), intent(in) :: file_name                               ! name of file
        integer, intent(inout), optional :: file_size                           ! file size
        integer, intent(inout), optional :: acc_time                            ! file access time
        integer, intent(inout), optional :: mod_time                            ! file modification time
        
        ! local variables
        integer :: vals(13)                                                     ! values
        integer :: istat                                                        ! status
        
        ! call stat
        call stat(trim(file_name),vals,istat)
        
        ! check status
        if (istat.ne.0) then
            call writo('WARNING: call to stat failed')
            if (present(file_size)) file_size = 0
            if (present(acc_time)) acc_time = 0
            if (present(mod_time)) mod_time = 0
        else
            if (present(file_size)) file_size = vals(8)
            if (present(acc_time)) acc_time = vals(9)
            if (present(mod_time)) mod_time = vals(10)
        end if
    end subroutine get_file_info
end module files_utilities
