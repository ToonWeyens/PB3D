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
        use eq_vars, only: eqd_mesh, tor_mesh, pol_mesh, calc_RZl
        use metric_ops, only: metric_C, metric_C2V, metric_V
        
        real(dp) :: alpha
        
        call writo('Start setting up equilibrium quantities')
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(1)
        
        call writo('Start determining the grid')
        call lvl_ud(1)
        
            ! determine the toroidal mesh points
            call tor_mesh
            
            ! calculate poloidal mesh points that follow the magnetic field line
            ! cfor urrent toroidal mesh points and field line (alpha)
            call pol_mesh(alpha)
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Grid determined')
            
            call writo('Calculating equilibrium quantities on grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! calculate the cylindrical variables R, Z and lambda and derivatives
            call calc_RZl
            
            ! calculate the metrics in the cylindrical coordinate system
            call metric_C
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call metric_C2V
            
            ! calculate  the  metric  factors in the VMEC coordinate system 
            call metric_V
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium quantities calculated on grid')
        
        call lvl_ud(-1)
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call writo('Done setting up equilibrium quantities')
    end subroutine
end module eq_ops
