!-------------------------------------------------------
!   Driver employing Richardson's extrapolation and normal discretization of 
!   the ODE's.
!-------------------------------------------------------
module driver_rich
    use num_vars, only: max_r
    use var_ops, only: i2str
    use output_ops, only: writo
    implicit none
    private
    public run_rich_driver
contains
    subroutine run_rich_driver()
        integer :: ir                                                           ! Richardson's level
        logical :: converged                                                    ! is it converged?
        
        ! Initalize some variables
        ir = 1
        converged = .false.
        
        Richard: do while (.not.converged .and. ir.le.max_r)
            call writo('Level ' // trim(i2str(ir)) // ' of Richardson''s &
                &extrapolation')
            ir = ir + 1
        end do Richard
    end subroutine
end module driver_rich
