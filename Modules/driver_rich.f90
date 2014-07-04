!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
    use num_vars, only: max_it_r, dp, pi
    use str_ops, only: i2str, r2str, r2strt
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud
    implicit none
    private
    public run_rich_driver

    ! global variables
    real(dp), allocatable :: alpha(:)
    
contains
    subroutine run_rich_driver()
        use num_vars, only: min_alpha, max_alpha, n_alpha
        use eq_vars, only: eqd_mesh
        use eq_ops, only: calc_eq
        use X_vars, only: n_r_X
        use matrix_X, only: calc_matrix_X
        
        integer :: ir, ia                                                       ! counters
        logical :: converged                                                    ! is it converged?
        
        ! determine the magnetic field lines for which to run the calculations 
        ! (equidistant mesh)
        if (allocated(alpha)) deallocate(alpha)
        allocate(alpha(n_alpha))
        alpha = eqd_mesh(n_alpha, min_alpha, max_alpha)                    ! just evenly spread them over 0..2*pi
        
        call writo('The calculations will be done for '//trim(i2str(n_alpha))&
            &//' values of alpha')
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(1)
        
        ! Do the calculations for every field line
        field_lines: do ia = 1, n_alpha
            ! Display message
            call writo(trim(i2str(ia))//'/'//trim(i2str(n_alpha))//&
                &': Calculations for field line alpha = '&
                &//trim(r2strt(alpha(ia)))//':')
            ! 2----------------------------------------------------------------
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
                
                ! Calculate the equilibrium quantities for current alpha
                call calc_eq(alpha(ia))
                
                ! Initalize some variables
                ir = 1
                converged = .false.
                
                ! Start Richardson loop
                call writo('Starting Richardson loop')
                ! 3------------------------------------------------------------
                call lvl_ud(1)
                
                    Richard: do while (.not.converged .and. ir.le.max_it_r)
                        call writo('Level ' // trim(i2str(ir)) // &
                            &' of Richardson''s extrapolation')
                        call lvl_ud(1)
                        ir = ir + 1
                        
                        !  calculate   number   of   radial   points   for   the
                        ! perturbation in Richardson loops and save in n_r_X
                        call calc_n_r_X
                        
                        ! calculate the plasma matrix elements mat_P
                        call writo('calculating the matrix elements due to &
                            &the plasma potential energy')
                        call lvl_ud(1)
                        call calc_matrix_X(n_r_X)
                        call lvl_ud(-1)
                        
                        call lvl_ud(-1)
                    end do Richard
                ! 3------------------------------------------------------------
                call lvl_ud(-1)
            
            ! 2----------------------------------------------------------------
            ! 2----------------------------------------------------------------
            call lvl_ud(-1)
        end do field_lines
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(-1)
    end subroutine
    
    ! calculates the number of normal  points for the perturbation n_r_X for the
    ! various Richardson iterations
    subroutine calc_n_r_X
        use X_vars, only: n_r_X
        
        call writo('TEMPORALLY SETTING n_r_X to 10 !!!')
        n_r_X = 10
    end subroutine
end module driver_rich
