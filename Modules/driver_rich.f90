!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
#include <PB3D_macros.h>
    use num_vars, only: max_it_r, dp, pi
    use str_ops, only: i2str, r2str, r2strt
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud
    implicit none
    private
    public run_rich_driver

    ! global variables
    real(dp), allocatable :: alpha(:)
    
contains
    subroutine run_rich_driver(ierr)
        use num_vars, only: min_alpha, max_alpha, n_alpha
        use eq_vars, only: eqd_mesh
        use eq_ops, only: calc_eq
        use X_ops, only: prepare_matrix_X, solve_EV_system
        use MPI_ops, only: split_MPI
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! input / output
        integer, intent(inout) :: ierr                                          ! error
        
        ! local variables
        integer :: ir, ia                                                       ! counters
        logical :: converged                                                    ! is it converged?
        
        ! initialize ierr
        ierr = 0
        
        ! split  the  communicator MPI_COMM_WORLD into subcommunicators
        call writo('Setting up groups for dynamical load balancing')
        call lvl_ud(1)
        call split_MPI(ierr)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! determine the magnetic field lines for which to run the calculations 
        ! (equidistant mesh)
        if (allocated(alpha)) deallocate(alpha)
        allocate(alpha(n_alpha))
        alpha = eqd_mesh(n_alpha,min_alpha,max_alpha,ierr)                      ! just evenly spread them over 0..2*pi
        CHCKERR('')
        
        call writo('The calculations will be done for '//trim(i2str(n_alpha))&
            &//' values of alpha')
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(1)
        
        ! Do the calculations for every field line
        
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!! DIVIDE DYNAMICALLY WITHIN THE GROUPS !!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        field_lines: do ia = 1, n_alpha
            ! Display message
            call writo(trim(i2str(ia))//'/'//trim(i2str(n_alpha))//&
                &': Calculations for field line alpha = '&
                &//trim(r2strt(alpha(ia)))//':')
            ! 2----------------------------------------------------------------
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
                
                ! Calculate the equilibrium quantities for current alpha
                call calc_eq(alpha(ia),ierr)
                CHCKERR('')
                
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
                        call writo('calculating the normal points')
                        call lvl_ud(1)
                        call calc_n_r_X
                        call lvl_ud(-1)
                        
                        ! prepare  the   matrix  elements  by   calculating  the
                        ! magnitudes KV and  PV for each of  the n_r equilibrium
                        ! normal points and for the modes (k,m)
                        call writo('calculating magnitudes KV and PV at these &
                            &normal points')
                        call lvl_ud(1)
                        call prepare_matrix_X
                        call lvl_ud(-1)
                        
                        ! setup the matrices of the generalized EV system
                        ! AX = lambda BX and solve it
                        call writo('treating the EV system')
                        call lvl_ud(1)
                        call solve_EV_system(ierr)
                        CHCKERR('')
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
        
        n_r_X = 50
        call writo('TEMPORALLY SETTING n_r_X to '//trim(i2str(n_r_X))//'!!!')
    end subroutine
end module driver_rich
