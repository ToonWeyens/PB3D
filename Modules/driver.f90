!------------------------------------------------------------------------------!
!   Main driver of program Peeling Ballooning in 3D, chooses between available !
!   drivers, according to the value of the parameter style:                    !
!       1. driver_rich, using Richardson's extrapolation and normal            !
!           discretization                                                     !
!------------------------------------------------------------------------------!
module driver
#include <PB3D_macros.h>
    use driver_rich, only: run_rich_driver
    use num_vars, only: style, max_str_ln
    use str_ops, only: i2str
    use output_ops, only: writo, lvl_ud
    implicit none
    private
    public run_driver

contains
    ! the main driver routine
    ! [MPI] All ranks
    integer function run_driver() result(ierr)
        character(*), parameter :: rout_name = 'run_driver'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! run the appropriate driver depending on "style"
        select case (style)
            ! Richardson's exptrapolation and normal discretization
            case (1)
                call writo('Numerical method chosen:')
                call lvl_ud(1)
                call writo('Richardson''s extrapolation and &
                    &normal discretization, generalized eigenvalue problem')
                call lvl_ud(-1)
                ierr = run_rich_driver()
                CHCKERR('')
            case default
                err_msg = 'Style "' // trim(i2str(style)) // '" is not&
                    & valid'
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function run_driver
end module driver