!-------------------------------------------------------
!   Driver employing Richardson's extrapolation and normal discretization of 
!   the ODE's.
!-------------------------------------------------------
module driver_rich
    use num_vars, only: max_r, dp, pi
    use var_ops, only: i2str
    use output_ops, only: writo, print_ar_2, lvl_ud
    implicit none
    private
    public run_rich_driver

contains
    subroutine run_rich_driver()
        use grid_vars, only: calc_ang_mesh, calc_RZl, &
            &theta, zeta, n_theta, n_zeta
        use metric_ops, only: metric_C, metric_C2V, metric_V
        
        integer :: ir                                                           ! Richardson's level
        logical :: converged                                                    ! is it converged?
        real(dp) :: min_theta, max_theta, min_zeta, max_zeta
        
        ! Initalize some variables
        ir = 1
        converged = .false.
        
        Richard: do while (.not.converged .and. ir.le.max_r)
            call writo('Level ' // trim(i2str(ir)) // ' of Richardson''s &
                &extrapolation')
            call lvl_ud(1)
            ir = ir + 1
            
            call writo('Calculate cylindrical metrics')
            ! calculate the grid points in current Richardson level
            n_theta = 10; min_theta = 0; max_theta = 2*pi
            n_zeta = 10; min_zeta = 0; max_zeta = 2*pi
            theta = calc_ang_mesh(n_theta, min_theta, max_theta)
            zeta = calc_ang_mesh(n_zeta, min_zeta, max_zeta)
            
            ! calculate the cylindrical variables R and Z and derivatives
            call calc_RZl
            
            ! calculate the metrics in the cylindrical coordinate system
            call metric_C
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call metric_C2V
            
            ! calculate  the  metric  factors in the VMEC coordinate system 
            call metric_V
            
            call lvl_ud(-1)
        end do Richard
    end subroutine
end module driver_rich
