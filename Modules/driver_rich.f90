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
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_rich_driver() result(ierr)
        use num_vars, only: min_alpha, max_alpha, n_alpha, glob_rank, &
            &n_groups, group_nr, group_rank, max_alpha, next_job, next_job_win
        use eq_vars, only: calc_eqd_mesh
        use MPI_ops, only: split_MPI, merge_MPI, get_next_job
        use MPI
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! initialize ierr
        ierr = 0
        
        ! determine the magnetic field lines for which to run the calculations 
        ! (equidistant mesh)
        if (allocated(alpha)) deallocate(alpha)
        allocate(alpha(n_alpha))
        ierr = calc_eqd_mesh(alpha,n_alpha,min_alpha,max_alpha)                 ! just evenly spread them over 0..2*pi
        CHCKERR('')
        
        ! split  the  communicator MPI_COMM_WORLD into subcommunicators
        call writo('Setting up groups for dynamical load balancing')
        call lvl_ud(1)
        ierr = split_MPI()
        CHCKERR('')
        call lvl_ud(-1)
        
        if (glob_rank.eq.0) then                                                ! only global rank
            call writo('The calculations will be done for '//&
                &trim(i2str(n_alpha))//' values of alpha')
        end if
        
        call lvl_ud(1)
            
        ! initialize next_job to 0
        next_job = 0
        
        !call MPI_barrier(MPI_COMM_WORLD,ierr)
        !write(*,*) 'stopping rank ', glob_rank
        !stop
        
        ! do the calculations for every  field line, where initially every group
        ! is assigned a field line. On completion of a job, the completing group
        ! inquires  all the other  groups using passive  1-sided MPI to  get the
        ! next group that has to be done
        ! loop over all fied lines
        field_lines: do
            ! get next job for current group
            ierr = get_next_job()
            CHCKERR('')
            
            ! output
            if (next_job.gt.0) then
                call writo('Job '//trim(i2str(next_job))//' assigned to group '//&
                    &trim(i2str(group_nr)))
            else
                call writo('Finished all jobs')
                exit
            end if
            
            ! Do the calculations for a field line alpha
            if (next_job.gt.0) then
                !ierr = run_for_alpha(alpha(next_job),next_job)
                CHCKERR('')
            end if
        end do field_lines
        
        call lvl_ud(-1)
        
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
    integer function run_for_alpha(alpha,i_alpha) result(ierr)
        use num_vars, only: group_nr, group_rank, n_alpha
        use eq_ops, only: calc_eq
        use X_ops, only: prepare_matrix_X, solve_EV_system
        
        character(*), parameter :: rout_name = 'run_for_alpha'
        
        ! input / output
        real(dp), intent(in) :: alpha                                           ! alpha at which to run the calculations
        integer, intent(in) :: i_alpha                                          ! which of the alpha's this is
        
        ! local variables
        integer :: ir                                                           ! counter for richardson extrapolation
        logical :: converged                                                    ! is it converged?
        
        ! initialize ierr
        ierr = 0
        
        ! initialize ierr
        ierr = 0
        
        ! Display message
        call writo(trim(i2str(i_alpha))//'/'//trim(i2str(n_alpha))//&
            &': Calculations for field line alpha = '//trim(r2strt(alpha))//&
            &', allocated to group '//trim(i2str(group_nr)))
        
        call lvl_ud(1)                                                          ! starting calculation for current fied line
                
        ! Calculate the equilibrium quantities for current alpha
        ierr = calc_eq(alpha)
        CHCKERR('')
        
        ! Initalize some variables
        ir = 1
        converged = .false.
        
        ! Start Richardson loop
        call writo('Starting Richardson loop')
        
        call lvl_ud(1)                                                          ! before richardson loop
        
        Richard: do while (.not.converged .and. ir.le.max_it_r)
            call writo('Level ' // trim(i2str(ir)) // &
                &' of Richardson''s extrapolation')
            call lvl_ud(1)                                                      ! beginning of one richardson loop
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
            ierr = solve_EV_system()
            CHCKERR('')
            call lvl_ud(-1)
            
            call lvl_ud(-1)                                                     ! end of one richardson loop
        end do Richard
            
        call lvl_ud(-1)                                                         ! done with richardson
        call writo('Finished Richardson loop')
            
        call lvl_ud(-1)                                                         ! done with calculation for current field line
        call writo(trim(i2str(i_alpha))//'/'//trim(i2str(n_alpha))//&
            &': Finished calculations for field line alpha = '//&
            &trim(r2strt(alpha)))
    end function run_for_alpha
    
    ! calculates the number of normal  points for the perturbation n_r_X for the
    ! various Richardson iterations
    subroutine calc_n_r_X
        use X_vars, only: n_r_X
        
        n_r_X = 50
        call writo('TEMPORALLY SETTING n_r_X to '//trim(i2str(n_r_X))//'!!!')
    end subroutine
end module driver_rich
