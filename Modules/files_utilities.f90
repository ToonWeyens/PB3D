!------------------------------------------------------------------------------!
!   Numerical utilities related to files                                       !
!------------------------------------------------------------------------------!
module files_utilities
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public nextunit, get_file_info, get_full_PB3D_name, delete_file, &
        &skip_comment, count_lines

contains
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
        integer, parameter :: lun_min=70, lun_max=1000
        logical :: file_open
        integer :: lun
        
        ! iterate over permitted luns until available one found
        nextunit=-1
        do lun=lun_min,lun_max
            inquire(UNIT=lun,OPENED=file_open)
            if (.not.file_open) then
                nextunit=lun
                exit
            end if
        end do
        
        ! return unit number if present
        if (present(unit)) unit=nextunit
    end function nextunit
    
    ! Skips comment.
    integer function skip_comment(file_i,file_name) result(ierr)
        character(*), parameter :: rout_name = 'skip_comment'
        
        ! input / output
        integer, intent(in) :: file_i                                           ! file identifier
        character(len=*), intent(in), optional :: file_name                     ! file name
        
        ! local variables
        character(len=1) :: loc_data_char                                       ! first data character
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        loc_data_char = '#'
        do while (loc_data_char.eq.'#')
            read(file_i,*,IOSTAT=ierr) loc_data_char
            if (present(file_name)) then
                err_msg = 'failed to read "'//trim(file_name)//'"'
            else
                err_msg = 'failed to read file'
            end if
            CHCKERR(err_msg)
        end do
        backspace(file_i)
    end function skip_comment
    
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
            call writo('call to stat failed',warning=.true.)
            if (present(file_size)) file_size = 0
            if (present(acc_time)) acc_time = 0
            if (present(mod_time)) mod_time = 0
        else
            if (present(file_size)) file_size = vals(8)
            if (present(acc_time)) acc_time = vals(9)
            if (present(mod_time)) mod_time = vals(10)
        end if
    end subroutine get_file_info
    
    ! Returns the name of the PB3D output file, including the current Richardson
    ! level. Optionally, this can be overridden If not positive, it is ignored.
    character(len=max_str_ln) function get_full_PB3D_name(rich_lvl) &
        &result(full_PB3D_name)
        use num_vars, only: PB3D_name
        
        ! input / output
        integer, intent(in), optional :: rich_lvl                               ! optional richardson level to be appended
        
        ! local variables
        integer :: rich_lvl_loc                                                 ! local richardson level
        
        ! set local rich_lvl
        rich_lvl_loc = 0
        if (present(rich_lvl)) rich_lvl_loc = rich_lvl
        
        ! set ouput
        if (rich_lvl_loc.gt.0) then
            full_PB3D_name = trim(PB3D_name)//'_R_'//&
                &trim(i2str(rich_lvl_loc))//'.h5'
        else
            full_PB3D_name = trim(PB3D_name)//'.h5'
        end if
    end function get_full_PB3D_name

    ! removes a file
    integer function delete_file(file_name) result(istat)
        ! input / output
        character(len=*), intent(inout) :: file_name                            ! the name that is deleted
        
        ! local variables
        integer :: file_i                                                       ! file unit
        logical :: file_open                                                    ! whether file is open
        
        ! initialize istat
        istat = 0 
        
        ! check if opened and open if not
        inquire(FILE=trim(file_name),OPENED=file_open,NUMBER=file_i,&
            &IOSTAT=istat)
        CHCKSTT
        if (.not.file_open) then
            open(UNIT=nextunit(file_i),FILE=trim(file_name),IOSTAT=istat)
        end if
        CHCKSTT
        
        ! remove open file
        close(UNIT=file_i,STATUS='delete',IOSTAT=istat)
        CHCKSTT
    end function delete_file
    
    ! count non-comment lines
    integer function count_lines(file_i) result(nr_lines)
        ! input / output
        integer, intent(in) :: file_i                                           ! file identifier
        
        ! local variables
        integer :: istat                                                        ! status
        character(len=1) :: loc_data_char                                       ! first data character
        
        ! initialize istat
        istat = 0
        
        nr_lines = 0
        rewind(file_i)
        do while (istat.eq.0)
            ! read first character of data
            read(file_i,*,IOSTAT=istat) loc_data_char
            if (istat.eq.0 .and. loc_data_char.ne.'#') then                     ! exclude comment lines
                nr_lines = nr_lines + 1
            end if
        end do
        rewind(file_i)
    end function count_lines
end module files_utilities
