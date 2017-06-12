!------------------------------------------------------------------------------!
!   Numerical utilities related to input                                       !
!------------------------------------------------------------------------------!
module input_utilities
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    
    implicit none
    private
    public pause_prog, get_real, get_int, get_log, dealloc_in
    
contains
    ! Queries for a logical value yes or no, where the default answer is also to
    ! be provided.
    ! [MPI] All ranks, but only master can give input
    function get_log(yes,ind) result(val)
        use num_vars, only: rank
        use MPI_utilities, only: broadcast_var
        
        ! input / output
        logical :: val                                                          ! output value
        logical :: yes
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        character(len=max_str_ln) :: answer_str                                 ! string with answer
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! only master can receive input
        if (rank.eq.0) then
            write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
            if (yes) then
                write(*,'(A)',advance='no') 'y(es)/n(o) [yes]: '
                val = .true.
            else
                write(*,'(A)',advance='no') 'y(es)/n(o) [no]: '
                val = .false.
            end if
            call stop_time
            read (*, '(A)') answer_str
            call start_time
            
            select case (strh2l(trim(answer_str)))
                !case ('y','yes')
                case ('y','Y','yes','Yes','YEs','YES','yEs','yES','yeS','YeS')
                    val = .true.
                case ('n','N','no','No','NO','nO') 
                    val = .false.
            end select
        end if
        
        ! if not individual, broadcast result
        if (.not.ind_loc) then
            istat = broadcast_var(val)
            if (istat.ne.0) call writo('In get_log, something went wrong. &
                &Default used.',warning=.true.)
        end if
    end function get_log
    
    ! Queries  for user  input for a  real value, where  allowable range  can be
    ! provided as well.
    ! [MPI] All ranks, but only global rank can give input
    function get_real(lim_lo,lim_hi,ind) result(val)
        use num_vars, only: rank
        use MPI_utilities, only: broadcast_var
        
        ! input / output
        real(dp) :: val                                                         ! output value
        real(dp), intent(in), optional :: lim_lo                                ! upper and lower limit
        real(dp), intent(in), optional :: lim_hi                                ! upper and lower limit
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        real(dp) :: lims_loc(2)                                                 ! local version of limits
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! initialize val
        val = 0._dp
        
        ! only master can receive input
        if (rank.eq.0) then
            ! set up local limits
            lims_loc = [-huge(1._dp),huge(1._dp)]
            if (present(lim_lo)) lims_loc(1) = lim_lo
            if (present(lim_hi)) lims_loc(2) = lim_hi
            
            ! get input
            do
                write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
                write(*,'(A)',advance='no') 'Input a value'
                if (present(lim_lo).or.present(lim_hi)) then
                    write(*,'(A)',advance='no') ' ['
                    if (present(lim_lo)) write(*,'(A)',advance='no') &
                        &trim(r2strt(lim_lo))
                    write(*,'(A)',advance='no') '..'
                    if (present(lim_hi)) write(*,'(A)',advance='no') &
                        &trim(r2strt(lim_hi))
                    write(*,'(A)',advance='no') ']'
                end if
                write(*,'(A)',advance='no') ': '
                call stop_time
                read (*,*,iostat=istat) val
                call start_time
                
                if (istat.ne.0 .or. val.lt.lims_loc(1) .or. val.gt. &
                    &lims_loc(2)) then
                    write(*,'(A)',advance='no') empty_str                       ! first print empty string so that output is visible
                    write(*,'(A)',advance='no') 'Choose a value between '//&
                        &trim(r2strt(lims_loc(1)))//' and '//&
                        &trim(r2strt(lims_loc(2)))
                    write(*,*) ''
                    cycle
                else
                    exit
                end if
            end do
        end if
        
        ! if not individual, broadcast result
        if (.not.ind_loc) then
            istat = broadcast_var(val)
            if (istat.ne.0) call writo('In get_real, something went &
                &wrong. Default of zero used.',warning=.true.)
        end if
    end function get_real
    
    ! Queries for user input for an  integer value, where allowable range can be
    ! provided as well.
    ! [MPI] All ranks, but only global rank can give input
    function get_int(lim_lo,lim_hi,ind) result(val)
        use num_vars, only: rank
        use MPI_utilities, only: broadcast_var
        
        ! input / output
        integer :: val                                                          ! output value
        integer, intent(in), optional :: lim_lo                                 ! upper and lower limit
        integer, intent(in), optional :: lim_hi                                 ! upper and lower limit
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        integer :: lims_loc(2)                                                ! local version of limits
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! initialize val
        val = 0
        
        ! only master can receive input
        if (rank.eq.0) then
            ! set up local limits
            lims_loc = [-huge(1),huge(1)]
            if (present(lim_lo)) lims_loc(1) = lim_lo
            if (present(lim_hi)) lims_loc(2) = lim_hi
            
            ! get input
            do
                write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
                write(*,'(A)',advance='no') 'Input a value'
                if (present(lim_lo).or.present(lim_hi)) then
                    write(*,'(A)',advance='no') ' ['
                    if (present(lim_lo)) write(*,'(A)',advance='no') &
                        &trim(i2str(lim_lo))
                    write(*,'(A)',advance='no') '..'
                    if (present(lim_hi)) write(*,'(A)',advance='no') &
                        &trim(i2str(lim_hi))
                    write(*,'(A)',advance='no') ']'
                end if
                write(*,'(A)',advance='no') ': '
                call stop_time
                read (*,*,iostat=istat) val
                call start_time
                
                if (istat.ne.0 .or. val.lt.lims_loc(1) .or. val.gt. &
                    &lims_loc(2)) then
                    write(*,'(A)',advance='no') empty_str                       ! first print empty string so that output is visible
                    write(*,'(A)',advance='no') 'Choose a value between '//&
                        &trim(i2str(lims_loc(1)))//' and '//&
                        &trim(i2str(lims_loc(2)))
                    write(*,*) ''
                    cycle
                else
                    exit
                end if
            end do
        end if
        
        ! if not individual, broadcast result
        if (.not.ind_loc) then
            istat = broadcast_var(val)
            if (istat.ne.0) call writo('In get_int, something went &
                &wrong. Default of zero used.',warning=.true.)
        end if
    end function get_int
    
    ! pauses the running of the program
    ! [MPI] All ranks or, optinally, only current rank
    subroutine pause_prog(ind)
        use MPI_utilities, only: wait_MPI
        use num_vars, only: rank, rank
        
        ! input / output
        logical, intent(in), optional :: ind                                    ! individual pause or not
        
        ! local variables
        character(len=11) :: empty_str = ''                                     ! empty string
        integer :: istat                                                        ! status
        logical :: ind_loc                                                      ! local version of ind
        character(len=max_str_ln) :: hidden_msg                                 ! hidden message
        
        ! output message
        if (rank.eq.0) then
            write(*,'(A)',advance='no') empty_str                               ! first print empty string so that output is visible
            write(*,'(A)',advance='no') 'Paused. Press enter...'
        end if
        
        ! set local ind
        ind_loc = .false.
        if (present(ind)) ind_loc = ind
        
        ! only master can receive input
        if (rank.eq.0) then
            call stop_time
            read (*,'(A)') hidden_msg
            call start_time
        end if
        
        ! wait for MPI
        if (.not.ind_loc) then
            istat = wait_MPI()
            if (istat.ne.0) call writo('In pause_prog, something went &
                &wrong. Continuing.',warning=.true.)
        end if
        
        ! hidden message
        if (trim(hidden_msg).eq.'stop') stop 0
    end subroutine pause_prog
    
    ! Cleans up input from equilibrium codes.
    subroutine dealloc_in()
        use num_vars, only: eq_style
        use VMEC_vars, only: dealloc_VMEC
        use HELENA_vars, only: dealloc_HEL
        
        ! deallocate depending on equilibrium style
        select case (eq_style)
            case (1)                                                            ! VMEC
                call dealloc_VMEC
            case (2)                                                            ! HELENA
                call dealloc_HEL
        end select
    end subroutine dealloc_in
end module input_utilities
