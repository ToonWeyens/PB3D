!------------------------------------------------------------------------------!
!   routines and functions that have to do with the computation time           !
!------------------------------------------------------------------------------!
module time
    use num_vars, only: max_str_ln, dp
    use str_ops, only: r2strt
    use output_ops, only: writo

    implicit none
    private
    public init_time, start_time, stop_time, passed_time

    real(dp) :: deltat                                                          ! length of time interval
    real(dp) :: t1, t2                                                          ! end points of time interval
    logical :: running                                                          ! whether the timer is running

contains
    ! intialize the time passed to 0
    subroutine init_time
        deltat = 0
        t1 = 0
        t2 = 0
        running = .false. 
    end subroutine

    ! start a timer
    subroutine start_time
        if (running) then
            call writo('WARNING: Tried to start timer, but was already running')
        else
            call second0(t1)
            running = .true.
        end if
    end subroutine

    ! stop a timer
    subroutine stop_time
        if (running) then
            call second0(t2)
            
            ! increase deltat
            deltat = deltat+t2-t1

            ! set t1 and t2 back to zero
            t1 = 0
            t2 = 0
            running = .false.
        else
            call writo('WARNING: Tried to stop timer, but was already stopped')
        end if
    end subroutine

    ! display the time that has passed between t1 and t2
    ! automatically stops time and resets everything to zero
    subroutine passed_time
        character(len=max_str_ln) :: begin_str, end_str

        ! stop at current time if running
        if (running) call stop_time

        begin_str = '(this took'
        if (deltat.lt.0.1) then
            end_str = ' less than 0.1 seconds)'
        else
            end_str = ' ' // trim(r2strt(deltat)) // ' seconds)'
        end if
        call writo(trim(begin_str) // trim(end_str))

        ! restart deltat
        call init_time
    end subroutine
end module time
