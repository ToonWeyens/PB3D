!-------------------------------------------------------
!   Calculates the equilibrium quantities, making use of the metric_ops, eq_vars, etc
!-------------------------------------------------------
module eq_ops
    use num_vars, only: pi, dp
    use output_ops, only: print_ar_2, lvl_ud, writo
    
    implicit none
    private
    public calc_eq

contains
    ! calculate the equilibrium quantities on a grid determined by straight field
    ! lines.
    subroutine calc_eq(alpha)
        use eq_vars, only: eqd_mesh, calc_mesh, calc_RZl, calc_flux_q, &
            &check_mesh, flux_brkdwn
        use metric_ops, only: metric_C, metric_C2V, metric_V, metric_V2F
        
        real(dp) :: alpha
        
        call writo('Start setting up equilibrium quantities')
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(1)
        
        call writo('Start determining the equilibrium grid')
        call lvl_ud(1)
            
            ! calculate mesh points (theta, zeta) that follow the magnetic field
            ! line
            call calc_mesh(alpha)
            
            ! check whether the mesh has been calculated correctly
            call check_mesh(alpha)
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium grid determined')
            
            call writo('Calculating equilibrium quantities on equilibrium grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! calculate the cylindrical variables R, Z and lambda and derivatives
            call calc_RZl
            
            ! calculate flux quantities
            call calc_flux_q
            
            ! find out where using the poloidal flux as normal coordinate breaks
            ! down and calculate the transformation points
            call flux_brkdwn
            
            ! calculate the metrics in the cylindrical coordinate system
            call metric_C
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call metric_C2V
            
            ! calculate  the  metric  factors in the VMEC coordinate system 
            call metric_V
            
            ! calculate the transformation matrix V(mec) -> F(lux)
            call metric_V2F
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium quantities calculated on equilibrium grid')
        
        call lvl_ud(-1)
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call writo('Done setting up equilibrium quantities')
    end subroutine
end module eq_ops
