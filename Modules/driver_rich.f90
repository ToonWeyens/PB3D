!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_rich
#include <PB3D_macros.h>
    use num_vars, only: max_it_r, dp, pi, max_str_ln
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
        use num_vars, only: min_alpha, max_alpha, n_alpha, glb_rank, grp_nr, &
            &max_alpha, alpha_job_nr, use_pol_flux
        use eq_vars, only: calc_eqd_mesh
        use MPI_ops, only: split_MPI, merge_MPI, get_next_job
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X
        use VMEC_vars, only: dealloc_VMEC
        
        character(*), parameter :: rout_name = 'run_rich_driver'
        
        ! local variables
        character(len=8) :: flux_name                                           ! toroidal or poloidal
        
        ! initialize ierr
        ierr = 0
        
        ! set up flux_name
        if (use_pol_flux) then
            flux_name = 'poloidal'
        else
            flux_name = 'toroidal'
        end if
        
        ! output concerning n_alpha
        call writo('The calculations will be done')
        call lvl_ud(1)
        call writo('using the '//trim(flux_name)//' flux as the normal &
            &variable')
        call writo('for '//trim(i2str(n_alpha))//' values of alpha')
        if (use_pol_flux) then
            call writo('with toroidal mode number n = '//trim(i2str(min_n_X)))
            call writo('and poloidal mode number m = '//trim(i2str(min_m_X))//&
                &'..'//trim(i2str(max_m_X)))
        else
            call writo('with poloidal mode number m = '//trim(i2str(min_m_X)))
            call writo('and toroidal mode number n = '//trim(i2str(min_n_X))//&
                &'..'//trim(i2str(max_n_X)))
        end if
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
            ierr = get_next_job(alpha_job_nr)
            CHCKERR('')
            
            ! Do the calculations for a field line alpha
            if (alpha_job_nr.gt.0) then
                ! display message
                call writo('Job '//trim(i2str(alpha_job_nr))//': Calculations &
                    &for field line alpha = '//&
                    &trim(r2strt(alpha(alpha_job_nr)))//&
                    &', allocated to group '//trim(i2str(grp_nr)))
                
                call lvl_ud(1)                                                  ! starting calculation for current fied line
                
                ! calculate
                ierr = run_for_alpha(alpha(alpha_job_nr))
                CHCKERR('')
                
                ! display message
                call lvl_ud(-1)                                                 ! done with calculation for current field line
                call writo('Job '//trim(i2str(alpha_job_nr))//&
                    &': Calculations for field line alpha = '//&
                    &trim(r2strt(alpha(alpha_job_nr)))//&
                    &', completed by group '//trim(i2str(grp_nr)))
            else
                call writo('Finished all jobs')
                exit field_lines
            end if
        end do field_lines
        
        ! deallocate VMEC variables
        call dealloc_VMEC
        
        ! merge  the subcommunicator into communicator MPI_COMM_WORLD
        if (glb_rank.eq.0) then
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
        use num_vars, only: n_sol_requested, max_it_r, grp_rank, no_guess, &
            &alpha_job_nr
        use eq_ops, only: calc_eq
        use eq_vars, only: ang_par_F, n_par
        use output_ops, only: draw_GP
        use eq_vars, only: dealloc_eq_final
        use X_ops, only: prepare_X, solve_EV_system, plot_X_vec
        use X_vars, only: X_vec, X_val, init_m, dealloc_X_final, n_r_X
        use metric_ops, only: dealloc_metric_final
        
        character(*), parameter :: rout_name = 'run_for_alpha'
        
        ! input / output
        real(dp), intent(in) :: alpha                                           ! alpha at which to run the calculations
        
        ! local variables
        integer :: ir                                                           ! counter for richardson extrapolation
        integer :: id                                                           ! counter
        logical :: done_richard                                                 ! is it converged?
        complex(dp), allocatable :: X_val_rich(:,:,:)                           ! Richardson array of eigenvalues X_val
        real(dp), allocatable :: x_axis(:,:)                                    ! x axis for plot of Eigenvalues with n_r_X
        logical :: use_guess                                                    ! whether a guess is formed from previous level of Richardson
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate the equilibrium quantities for current alpha
        ierr = calc_eq(alpha)
        CHCKERR('')
        
        ! initialize m
        ierr = init_m()
        CHCKERR('')
        
        ! prepare the  matrix elements by calculating  the integrated magnitudes
        ! KV_int and  PV_int for each of  the n_r equilibrium normal  points and
        ! for the modes (k,m)
        call prepare_X
        
        ! Initalize some variables for Richardson loop
        ir = 1
        done_richard = .false.
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            allocate(X_val_rich(1:max_it_r,1:max_it_r,1:n_sol_requested))
            X_val_rich = 0.0_dp
            allocate(x_axis(1:max_it_r,1:n_sol_requested))
        end if
        
        ! Start Richardson loop
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Starting perturbation calculation with Richardson &
                &extrapolation')
        else
            call writo('Starting perturbation calculation')
        end if
        call lvl_ud(1)                                                          ! before richardson loop
        
        ! initialize use_guess for first Richardson loop
        use_guess = .true.
        
        Richard: do while (.not.done_richard .and. ir.le.max_it_r)
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call writo('Level ' // trim(i2str(ir)) // &
                    &' of Richardson extrapolation')
                call lvl_ud(1)                                                  ! beginning of one richardson loop
            end if
            
            !  calculate  number  of  radial  points  for  the  perturbation  in
            ! Richardson loops and save in n_r_X
            call writo('calculating the normal points')
            call lvl_ud(1)
            call calc_n_r_X(ir)
            call lvl_ud(-1)
            
            ! set use_guess to .false. if no_guess
            if (no_guess) use_guess = .false.
            
            ! setup the matrices of the generalized EV system AX = lambda BX and
            ! solve it
            call writo('treating the EV system')
            call lvl_ud(1)
            ierr = solve_EV_system(use_guess)
            CHCKERR('')
            call lvl_ud(-1)
            
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                ! update the x axis of the Eigenvalue plot
                x_axis(ir,:) = 1.0_dp*n_r_X
                
                ! update  the  variable  X_val_rich  with the  results  of  this
                ! Richardson level
                ierr = calc_rich_ex(ir,X_val,X_val_rich,done_richard,use_guess)
                CHCKERR('')
                call writo('updating Richardson extrapolation variables')
                call lvl_ud(1)
                call lvl_ud(-1)
            else                                                                ! if not, Richardson is done
                done_richard = .true.
            end if
            
            ! output the evolution of the  Eigenvalues with the number of normal
            ! points in the perturbation grid  and the largest Eigenfunction for
            ! the last Richardson loop
            if (done_richard .or. ir.eq.max_it_r+1) then
                if (grp_rank.eq.0) then
                    if (max_it_r.gt.1) then
                        call writo('plotting the Eigenvalues')
                        
                        ! output on screen
                        plot_title = 'JOB '//trim(i2str(alpha_job_nr))//' - &
                            &Eigenvalues as function of nr. of normal points'
                        call print_GP_2D(plot_title,'Eigenvalues_'//&
                            &trim(i2str(alpha_job_nr))//'.dat',realpart(&
                            &X_val_rich(1:ir,1,:)),x=x_axis(1:ir,:))
                        ! same output in file as well
                        call draw_GP(plot_title,'Eigenvalues_'//&
                            &trim(i2str(alpha_job_nr))//'.dat',&
                            &n_sol_requested,.true.,.false.)
                    end if
                end if
                
                call writo('plotting the Eigenvectors')
                
                call lvl_ud(1)
                
                do id = 1,n_sol_requested
                    call writo('plotting results for mode '//trim(i2str(id))//&
                        &'/'//trim(i2str(n_sol_requested))//&
                        &', with eigenvalue '&
                        &//trim(r2strt(realpart(X_val(id))))//' + '//&
                        &trim(r2strt(imagpart(X_val(id))))//' i')
                    
                    ierr = plot_X_vec(X_vec(:,:,id),X_val(id),id,&
                        &[ang_par_F(1,1),ang_par_F(n_par,1)])
                    CHCKERR('')
                end do
                
                call lvl_ud(-1)
            end if
            
            if (max_it_r.gt.1) then                                             ! only do this if more than 1 Richardson level
                call lvl_ud(-1)                                                 ! end of one richardson loop
            end if
        end do Richard
        
        call lvl_ud(-1)                                                         ! done with richardson
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            call writo('Finished Richardson loop')
        else
            call writo('Finished perturbation calculation')
        end if
        
        ! deallocate Richardson loop variables
        if (max_it_r.gt.1) then                                                 ! only do this if more than 1 Richardson level
            deallocate(X_val_rich)
        end if
        
        ! deallocate remaining equilibrium quantities
        call writo('Deallocate remaining quantities')
        call dealloc_eq_final
        call dealloc_metric_final
        call dealloc_X_final
    end function run_for_alpha
    
    ! calculates the number of normal  points for the perturbation n_r_X for the
    ! various Richardson iterations
    ! The aim is to  halve the step size, which is given by  dx(n) = 1/(n-1) or,
    ! inverting: n(dx) = 1 + 1/dx.
    ! This yields n(dx/2)/n(dx) = (2+dx)/(1+dx) = (2n(dx)-1)/n(dx)
    ! The recursion formula is therefore: n(dx/2) = 2n(dx) - 1
    subroutine calc_n_r_X(ir)
        use X_vars, only: n_r_X
        use num_vars, only: min_n_r_X
        
        ! input / output
        integer, intent(in) :: ir
        
        if (ir.eq.1) then
            n_r_X = min_n_r_X
        else
            n_r_X = 2 * n_r_X - 1
        end if
        call writo(trim(i2str(n_r_X))//' normal points for this level')
    end subroutine
    
    ! calculates  the  coefficients  of   the  Eigenvalues  in   the  Richardson
    ! extrapolation
    ! This is done using the recursive formula
    !   X_val_rich(ir,ir2,:) = X_val_rich(ir,ir2-1,:) +  1/(2^(2ir2) - 1) * 
    !       (X_val_rich(ir,ir2-1,:) - X_val_rich(ir-1,ir2-1,:)),
    ! as described in [ADD REF]
    ! [MPI] All ranks
    integer function calc_rich_ex(ir,X_val,X_val_rich,done_richard,&
            &use_guess_for_next_level) result(ierr)
        use num_vars, only: tol_r
        
        character(*), parameter :: rout_name = 'calc_rich_ex'
        
        ! input / output
        integer, intent(inout) :: ir                                            ! level of Richardson extrapolation (starting at 1)
        complex(dp), intent(in) :: X_val(:)                                     ! EV for this Richardson level
        complex(dp), intent(inout) :: X_val_rich(:,:,:)                         ! extrapolated coefficients
        logical, intent(inout) :: done_richard                                  ! if Richardson loop has converged sufficiently
        logical, intent(inout) :: use_guess_for_next_level                      ! if a guessed is used for next Richardson level
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: ir2                                                          ! counter
        complex(dp), allocatable :: corr(:)                                     ! correction
        real(dp) :: max_corr                                                    ! maximum of maximum of correction
        real(dp), save :: prev_max_corr                                         ! max_corr at previous Richardson level
        integer :: loc_max_corr(1)                                              ! location of maximum of correction
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(X_val_rich,1).ne.size(X_val_rich,2) .or. &
            &size(X_val_rich,3).ne.size(X_val) .or. &
            &ir.gt.size(X_val_rich,1)) then
            ierr = 1
            err_msg = 'X_val_rich has to have correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! allocate correction
        allocate(corr(size(X_val))); corr = 1.0E15
        
        ! do calculations if ir > 1
        X_val_rich(ir,1,:) = X_val
        do ir2 = 2,ir
            corr = 1./(2**(2*ir2)-1.) * &
                &(X_val_rich(ir,ir2-1,:) - X_val_rich(ir-1,ir2-1,:))
            X_val_rich(ir,ir2,:) = X_val_rich(ir,ir2-1,:) + corr
        end do
        
        if (ir.gt.1) then                                                       ! only do this if in Richardson level higher than 1
            ! get maximum and location of maximum for relative correction
            max_corr = maxval(abs(corr/X_val_rich(ir,ir,:)))
            loc_max_corr = maxloc(abs(corr/X_val_rich(ir,ir,:)))
            call writo('maximum relative error: '//trim(r2strt(max_corr))//&
                &' for Eigenvalue '//trim(i2str(loc_max_corr(1))))
            
            ! check whether tolerance has been reached for last value ir2 = ir
            if (maxval(abs(corr/X_val_rich(ir,ir,:))) .lt. tol_r) then
                done_richard = .true.
                call writo('tolerance '//trim(r2strt(tol_r))//&
                    &' reached after '//trim(i2str(ir))//' iterations')
            else
                call writo('tolerance '//trim(r2strt(tol_r))//' not yet &
                    &reached')
            end if
            
            ! determine use_guess_for_next_level
            use_guess_for_next_level = max_corr.lt.prev_max_corr                ! set guess for next level to true if max_corr is decreasing
            prev_max_corr = max_corr                                            ! save max_corr for next level
        else
            use_guess_for_next_level = .true.                                   ! for first Richardson level, set guess for next level to true
        end if
        
        ! check for convergence
        if (.not.done_richard) then
            if (ir.lt.max_it_r) then                                            ! not yet at maximum Richardson iteration
                ir = ir + 1
            else                                                                ! maximum nr. of Richardson iterations reached
                call writo('maximum number of Richardson iterations reached')
                done_richard = .true.
            end if
        end if
    end function calc_rich_ex
end module driver_rich
