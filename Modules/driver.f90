!------------------------------------------------------------------------------!
!   Main driver of program Peeling Ballooning in 3D, chooses between available !
!   drivers, according to the value of the parameter style:                    !
!       1. driver_rich, using Richardson's extrapolation and normal            !
!           discretization                                                     !
!------------------------------------------------------------------------------!
module driver
    use driver_rich, only: run_rich_driver
    use num_vars, only: style
    use str_ops, only: i2str
    use output_ops, only: writo, lvl_ud
    implicit none
    private
    public run_driver

contains
    subroutine run_driver()
        select case (style)
            ! Richardson's exptrapolation and normal discretization
            case (1)
                call writo('Numerical method chosen:')
                call lvl_ud(1)
                call writo('Richardson''s extrapolation and &
                    &normal discretization, generalized eigenvalue problem')
                call run_rich_driver
                call lvl_ud(-1)
            case default
                call writo('ERROR: style "' // trim(i2str(style)) // '" is not&
                    & valid')
                stop
        end select
    end subroutine
end module driver
