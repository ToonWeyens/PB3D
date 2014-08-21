!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
#include <PB3D_macros.h>
    use num_vars, only: max_it_r, dp, pi
    use str_ops, only: i2str, r2str, r2strt
    use output_ops, only: writo, print_ar_2, print_ar_1, lvl_ud, print_GP_2D
    implicit none
    private
    public run_rich_driver

    ! global variables
    real(dp), allocatable :: alpha(:)
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_rich_driver() result(ierr)
        use num_vars, only: min_alpha, max_alpha, n_alpha, glob_rank, &
            &group_nr, max_alpha
        use eq_vars, only: calc_eqd_mesh
        use MPI_ops, only: split_MPI, merge_MPI, get_next_job
        use X_vars, only: n_X, min_m_X, max_m_X
        use VMEC_vars, only: dealloc_VMEC_vars
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! local variables
        integer :: next_job                                                     ! holds next job of current group, filled with get_next_job
        
        ! initialize ierr
        ierr = 0
        
        ! output concerning n_alpha
        call writo('The calculations will be done for '//&
            &trim(i2str(n_alpha))//' values of alpha with:')
        call lvl_ud(1)
        call writo('toroidal mode number n = '//trim(i2str(n_X)))
        call writo('poloidal mode number m = '//trim(i2str(min_m_X))//'..'//&
            &trim(i2str(max_m_X)))
        call lvl_ud(-1)
        
        ! determine the magnetic field lines for which to run the calculations 
        ! (equidistant mesh)
        allocate(alpha(n_alpha))
        ierr = calc_eqd_mesh(alpha,n_alpha,min_alpha,max_alpha)                 ! just evenly spread them over 0..2*pi
        CHCKERR('')
        
        ! split  the  communicator MPI_COMM_WORLD into subcommunicators
        call writo('Setting up groups for dynamical load balancing')
        call lvl_ud(1)
        ierr = split_MPI()
        CHCKERR('')
        call lvl_ud(-1)
        
        ! do the calculations for every  field line, where initially every group
        ! is assigned a field line. On completion of a job, the completing group
        ! inquires  all the other  groups using passive  1-sided MPI to  get the
        ! next group that has to be done
        ! loop over all fied lines
        field_lines: do
            ! get next job for current group
            ierr = get_next_job(next_job)
            CHCKERR('')
            
            ! Do the calculations for a field line alpha
            if (next_job.gt.0) then
                ! display message
                call writo('Job '//trim(i2str(next_job))//': Calculations for &
                    &field line alpha = '//trim(r2strt(alpha(next_job)))//&
                    &', allocated to group '//trim(i2str(group_nr)))
                
                call lvl_ud(1)                                                  ! starting calculation for current fied line
                
                ! calculate
                ierr = run_for_alpha(alpha(next_job))
                CHCKERR('')
                
                ! display message
                call lvl_ud(-1)                                                 ! done with calculation for current field line
                call writo('Job '//trim(i2str(next_job))//': Calculations for &
                    &field line alpha = '//trim(r2strt(alpha(next_job)))//&
                    &', completed by group '//trim(i2str(group_nr)))
            else
                call writo('Finished all jobs')
                exit field_lines
            end if
        end do field_lines
        
        ! deallocate VMEC variables
        call dealloc_VMEC_vars
        
        ! merge  the subcommunicator into communicator MPI_COMM_WORLD
        if (glob_rank.eq.0) then
            call writo('Stability analysis concluded at all '//&
                &trim(i2str(n_alpha))//' fieldlines')
            call writo('Merging groups for dynamical load balancing back &
                &together')
        end if
        call lvl_ud(1)
        ierr = merge_MPI()
        CHCKERR('')
        call lvl_ud(-1)
    end function run_rich_driver
    
    ! runs the calculations for one of the alpha's
    integer function run_for_alpha(alpha) result(ierr)
        use num_vars, only: n_sol_requested, max_it_r, group_rank
        use eq_ops, only: calc_eq
        use eq_vars, only: dealloc_eq_vars_final
        use X_ops, only: prepare_matrix_X, solve_EV_system
        use X_vars, only: dealloc_X_vars, X_val
        use metric_ops, only: dealloc_metric_vars_final
        
        character(*), parameter :: rout_name = 'run_for_alpha'
        
        ! input / output
        real(dp), intent(in) :: alpha                                           ! alpha at which to run the calculations
        
        ! local variables
        integer :: ir                                                           ! counter for richardson extrapolation
        integer :: id                                                           ! counter
        logical :: converged                                                    ! is it converged?
        complex(dp), allocatable :: X_val_rich(:,:)                             ! Richardson array of eigenvalues X_val
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate the equilibrium quantities for current alpha
        ierr = calc_eq(alpha)
        CHCKERR('')
        
        ! Initalize some variables for Richardson loop
        ir = 1
        converged = .false.
        allocate(X_val_rich(1:n_sol_requested,1:max_it_r))
        
        ! Start Richardson loop
        call writo('Starting Richardson loop')
        
        call lvl_ud(1)                                                          ! before richardson loop
        
        Richard: do while (.not.converged .and. ir.le.max_it_r)
            call writo('Level ' // trim(i2str(ir)) // &
                &' of Richardson''s extrapolation')
            call lvl_ud(1)                                                      ! beginning of one richardson loop
            
            !  calculate   number   of   radial   points   for   the
            ! perturbation in Richardson loops and save in n_r_X
            call writo('calculating the normal points')
            call lvl_ud(1)
            call calc_n_r_X(ir)
            call lvl_ud(-1)
            
            ! prepare  the   matrix  elements  by   calculating  the
            ! magnitudes KV and  PV for each of  the n_r equilibrium
            ! normal points and for the modes (k,m)
            call writo('calculating magnitudes KV and PV at these &
                &normal points')
            call lvl_ud(1)
            ierr = prepare_matrix_X()
            CHCKERR('')
            call lvl_ud(-1)
            
            ! setup the matrices of the generalized EV system
            ! AX = lambda BX and solve it
            call writo('treating the EV system')
            call lvl_ud(1)
            ierr = solve_EV_system()
            CHCKERR('')
            call lvl_ud(-1)
            
            ! save the output eigenvalue for this Richardson level
            if (group_rank.eq.0) then
                do id = 1,n_sol_requested
                    X_val_rich(id,ir) = X_val(id)
                end do
            end if
            
            !!!!! PERFORM EXTRAPOLATION OF EV TO GET A GOOD APPROX !!!!
            write(*,*) '!!! Extrapolate EV !!!!'
            
            ! deallocate perturbation variables
            call writo('Deallocating perturbation variables')
            call dealloc_X_vars
            
            ! increment counter
            ir = ir + 1
            
            call lvl_ud(-1)                                                     ! end of one richardson loop
        end do Richard
        
        call lvl_ud(-1)                                                         ! done with richardson
        call writo('Finished Richardson loop')
        
        ! output
        call writo('Plot output')
        if (group_rank.eq.0) then
            call print_GP_2D('RE X_val(1)','',realpart(X_val_rich(1,:)))
            call print_GP_2D('IM X_val(1)','',imagpart(X_val_rich(1,:)))
        end if
        
        ! deallocate remaining equilibrium quantities
        call writo('Deallocate remaining quantities')
        call dealloc_eq_vars_final
        call dealloc_metric_vars_final
    end function run_for_alpha
    
    ! calculates the number of normal  points for the perturbation n_r_X for the
    ! various Richardson iterations
    ! The aim is to  halve the step size, which is given by  dx(n) = 1/(n-1) or,
    ! inverting: n(dx) = 1 + 1/dx.
    ! This yields n(dx/2)/n(dx) = (2+dx)/(1+dx) = (2n(dx)-1)/n(dx)
    ! The recursion formula is therefore: n(dx/2) = 2n(dx) - 1
    subroutine calc_n_r_X(ir)
        use X_vars, only: n_r_X
        
        ! input / output
        integer, intent(in) :: ir
        
        if (ir.eq.1) then
            n_r_X = 10
        else
            n_r_X = 2 * n_r_X - 1
        end if
        call writo(trim(i2str(n_r_X))//' normal points for this level')
        write(*,*) 'UPDATE THE MATRICES!!!!! DO NOT JUST DELETE THEM! THEY CAN BE REUSED!!!!!'
    end subroutine
end module driver_rich
