!------------------------------------------------------------------------------!
!   Driver of the equilibrium part of PB3D.                                    !
!------------------------------------------------------------------------------!
module driver_eq
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    
    implicit none
    private
    public run_driver_eq
    
contains
    ! Main driver of PB3D_eq.
    integer function run_driver_eq() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, plot_flux_q
        use MPI_utilities, only: wait_MPI
        use eq_vars, only: dealloc_eq
        use VMEC, only: dealloc_VMEC
        use HELENA, only: dealloc_HEL
        use grid_vars, only: dealloc_grid, &
            &alpha
        use eq_ops, only: calc_eq, calc_derived_q, print_output_eq, &
            &flux_q_plot
        use met_ops, only: calc_met, calc_F_derivs
        use met_vars, only: dealloc_met
        use grid_ops, only: setup_and_calc_grid_B
        !!use utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver_eq'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        character(len=max_str_ln) :: err_msg                                    ! error message
        type(grid_type) :: grid_eq                                              ! equilibrium grid
        type(grid_type), pointer :: grid_eq_B => null()                         ! field-aligned equilibrium grid
        type(eq_type) :: eq                                                     ! equilibrium for
        type(met_type) :: met                                                   ! metric variables
        
        ! initialize ierr
        ierr = 0
        
        ! some preliminary things
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('The equilibrium variables are processed')
        call lvl_ud(1)
        
        if (use_pol_flux_F) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        call writo('for alpha = '//trim(r2strt(alpha*pi)))
        call writo('using the '//trim(flux_name)//' flux as the normal &
            &variable')
        
        call lvl_ud(-1)
        
        ! Calculate the equilibrium quantities
        ierr = calc_eq(grid_eq,eq)
        CHCKERR('')
        
        ! Calculate the metric quantities
        ierr = calc_met(grid_eq,eq,met)
        CHCKERR('')
        
        ! Transform E into F derivatives
        ierr = calc_F_derivs(grid_eq,eq,met)
        CHCKERR('')
        
        ! plot flux quantities if requested
        if (plot_flux_q) then
            ierr = flux_q_plot(grid_eq,eq)
            CHCKERR('')
        else
            call writo('Flux quantities plot not requested')
        end if
        
        ! Calculate derived metric quantities
        call calc_derived_q(grid_eq,eq,met)
        
        ! set up field-aligned equilibrium grid
        ierr = setup_and_calc_grid_B(grid_eq,grid_eq_B,eq)
        CHCKERR('')
        
        ! write equilibrium variables to output
        ierr = print_output_eq(grid_eq,grid_eq_B,eq,met)
        CHCKERR('')
        
        ! deallocate variables
        call writo('Deallocating variables')
        call lvl_ud(1)
        call dealloc_grid(grid_eq)
        nullify(grid_eq_B)
        call dealloc_eq(eq)
        call dealloc_met(met)
        ! deallocate variables depending on equilibrium style
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                call dealloc_VMEC
            case (2)                                                            ! HELENA
                call dealloc_HEL
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('Variables deallocated')
    end function run_driver_eq
end module driver_eq
