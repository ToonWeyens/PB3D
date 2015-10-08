!------------------------------------------------------------------------------!
!   Driver employing Richardson's extrapolation and normal discretization of   !
!   the ODE's.                                                                 !
!------------------------------------------------------------------------------!
module driver_PERT
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type, X_2_type
    
    implicit none
    private
    public run_driver_PERT
#if ldebug
    public debug_X_grid
#endif
    
    ! global variables
#if ldebug
    logical :: debug_X_grid = .false.                                           ! plot debug information for treatment of X grid
#endif
    
contains
    ! Implementation   of  the  driver,  using  Richardson's  extrapolation  and
    ! discretization in the normal direction
    ! [MPI] All ranks, parts are done  by global master only, other parts by the
    !       other ranks only
    !       Here, the tasks are split and allocated dynamically by the master to 
    !       the other groups
    integer function run_driver_PERT() result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, max_it_r, rank, &
            &plot_resonance, X_job_nr, X_jobs_lims
        use MPI_utilities, only: wait_MPI
        use X_vars, only: dealloc_X, &
            &min_m_X, max_m_X, min_n_X, max_n_X, min_r_X, max_r_X, min_n_r_X
        use grid_vars, only: dealloc_grid, &
            &n_par_X, min_par_X, max_par_X
        use PB3D_ops, only: read_PB3D, reconstruct_PB3D
        use PB3D_vars, only: dealloc_PB3D, &
            &PB3D_type
        use utilities, only: test_max_memory
        use MPI_ops, only: divide_X_jobs, get_next_job, print_jobs_info
        use X_ops, only: calc_X, check_X_modes, resonance_plot, print_output_X
        use vac, only: calc_vac
        !!use utilities, only: calc_aux_utilities
        
        character(*), parameter :: rout_name = 'run_driver'
        
        ! local variables
        type(PB3D_type), target :: PB3D                                         ! output PB3D for which to do postprocessing
        type(PB3D_type), pointer :: PB3D_B => null()                            ! PB3D variables on a field-aligned grid
        type(X_1_type) :: X_a, X_b                                              ! vectorial X variables
        type(X_2_type) :: X                                                     ! tensorial X variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! some preliminary things
        
        ! test maximum memory
        ierr = test_max_memory()
        CHCKERR('')
        
        !!! calculate auxiliary quantities for utilities
        !!call calc_aux_utilities                                                 ! calculate auxiliary quantities for utilities
        
        ! user output
        call writo('Perturbations are analyzed')
        call lvl_ud(1)
        
        if (max_it_r.eq.1) then
            call writo('for '//trim(i2str(min_n_r_X))//' values on &
                &normal range '//trim(r2strt(min_r_X))//'..'//&
                &trim(r2strt(max_r_X)))
        else
            call writo('for minimally '//trim(i2str(min_n_r_X))//' values on &
                &normal range '//trim(r2strt(min_r_X))//'..'//&
                &trim(r2strt(max_r_X)))
        end if
        call writo('for '//trim(i2str(n_par_X))//' values on parallel &
            &range '//trim(r2strt(min_par_X*pi))//'..'//&
            &trim(r2strt(max_par_X*pi)))
        if (use_pol_flux_F) then
            call writo('with toroidal mode number n = '//trim(i2str(min_n_X)))
            call writo('and poloidal mode number m = '//trim(i2str(min_m_X))//&
                &'..'//trim(i2str(max_m_X)))
        else
            call writo('with poloidal mode number m = '//trim(i2str(min_m_X)))
            call writo('and toroidal mode number n = '//trim(i2str(min_n_X))//&
                &'..'//trim(i2str(max_n_X)))
        end if
        
        call lvl_ud(-1)
        
        ! read PB3D output file
        ierr = read_PB3D(.false.,.true.,.false.,.false.,.false.)                ! read the equilibrium files
        CHCKERR('')
        
        ! reconstructing grids depends on equilibrium style
        ! user output
        call writo('Reconstructing PB3D output on output grid')
        call lvl_ud(1)
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! the field-aligned grid is identical to the output grid
                PB3D_B => PB3D
                ! normal call to reconstruct_PB3D
                ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.false.,&
                    &PB3D)
                CHCKERR('')
            case (2)                                                            ! HELENA
                ! the field-aligned grid is different form the output grid
                allocate(PB3D_B)
                ! additionally need field-aligned equilibrium grid
                ierr = reconstruct_PB3D(.false.,.true.,.false.,.false.,.false.,&
                    &PB3D,PB3D_B%grid_eq)
                CHCKERR('')
            case default
                ierr = 1
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        call writo('PB3D output reconstructed')
        
        ! tests
        ierr = check_X_modes(PB3D%eq)
        CHCKERR('')
        
        ! plot resonances if requested
        if (plot_resonance .and. rank.eq.0) then
            ierr = resonance_plot(PB3D%eq,PB3D%grid_eq)
            CHCKERR('')
        else
            call writo('Resonance plot not requested')
        end if
        
        ! divide perturbation jobs
        ierr = divide_X_jobs(PB3D%grid_eq,1)
        CHCKERR('')
        
        ! main loop over jobs or order 1
        X_job_nr = 0
        X_jobs_1: do
            ! get next job
            ierr = get_next_job(X_job_nr)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call writo('Job '//trim(i2str(X_job_nr))//' is started by process '&
                &//trim(i2str(rank)))
            call lvl_ud(1)
            
            ! calculate X variables, vector phase
            ierr = calc_X(PB3D%grid_eq,PB3D%eq,PB3D%met,X_a,&
                &lim_sec_X=X_jobs_lims(:,X_job_nr))
            CHCKERR('')
            
            ! write vectorial perturbation variables to output
            ierr = print_output_X(PB3D%grid_eq,X_a)
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X_a)
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed by process '&
                &//trim(i2str(rank)))
        end do X_jobs_1
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        ! divide perturbation jobs, tensor phase
        ierr = divide_X_jobs(PB3D%grid_eq,2)
        CHCKERR('')
        
        ! main loop over jobs or order 2
        X_job_nr = 0
        X_jobs_2: do
            ! get next job
            ierr = get_next_job(X_job_nr)
            CHCKERR('')
            if (X_job_nr.lt.0) exit
            
            ! user output
            call writo('Job '//trim(i2str(X_job_nr))//' is started by process '&
                &//trim(i2str(rank)))
            call lvl_ud(1)
            
            write(*,*) 'TEMPORARILY SLOWING PROC NE 0 DOWDN !!!!!!!!!!!!!!!!!!!'
            call sleep(rank*10)
            
            ! user output
            call writo('Requesting vectorial perturbation variables for &
                &dimension 1')
            call lvl_ud(1)
            
            ! read PB3D output file for dimension 1
            ierr = read_PB3D(.false.,.false.,.true.,.false.,.false.,&
                &lim_sec_X_1=X_jobs_lims(1:2,X_job_nr))
            CHCKERR('')
            
            ! normal call to reconstruct_PB3D X_1 quantities for dimension 1
            ierr = reconstruct_PB3D(.false.,.false.,.true.,.false.,.false.,&
                &PB3D,X_1=X_a,lim_sec_X_1=X_jobs_lims(1:2,X_job_nr))
            CHCKERR('')
            
            ! user output
            call lvl_ud(-1)
            call writo('Vectorial perturbation variables for dimension 1 &
                &loaded')
            
            ! user output
            call writo('Requesting vectorial perturbation variables for &
                &dimension 2')
            call lvl_ud(1)
            
            ! read PB3D output file for dimension 2
            ierr = read_PB3D(.false.,.false.,.true.,.false.,.false.,&
                &lim_sec_X_1=X_jobs_lims(3:4,X_job_nr))
            CHCKERR('')
            
            ! normal call to reconstruct_PB3D X_1 quantities for dimension 2
            ierr = reconstruct_PB3D(.false.,.false.,.true.,.false.,.false.,&
                &PB3D,X_1=X_b,lim_sec_X_1=X_jobs_lims(3:4,X_job_nr))
            CHCKERR('')
            
            ! user output
            call lvl_ud(-1)
            call writo('Vectorial perturbation variables for dimension 2 &
                &loaded')
            
            ! calculate X variables, tensor phase
            ierr = calc_X(PB3D%grid_eq,PB3D%eq,PB3D%met,X_a,X_b,X,&
                &lim_sec_X=reshape(X_jobs_lims(:,X_job_nr),[2,2]))
            CHCKERR('')
            
            ! calculate vacuum response
            ierr = calc_vac(X)
            CHCKERR('')
            
            ! integrate tensorial perturbation variables over field-aligned grid
            !!!!!!!!!!!!!!!!!
            
            ! write tensorial perturbation variables to output file
            ierr = print_output_X(PB3D%grid_eq,X)
            CHCKERR('')
            
            ! clean up
            call dealloc_X(X_a)
            call dealloc_X(X_b)
            call dealloc_X(X)
            
            ! user output
            call lvl_ud(-1)
            call writo('Job '//trim(i2str(X_job_nr))//' completed by process '&
                &//trim(i2str(rank)))
            
            write(*,*) 'waiting'
            if (rank.eq.0) read(*,*) 
        end do X_jobs_2
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! print jobs information from other processes
        ierr = print_jobs_info()
        CHCKERR('')
        
        ! deallocate PB3D variable
        call dealloc_PB3D(PB3D)
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
    end function run_driver_PERT
    
    
    !contains
        !! calculates the number of normal  points for the perturbation n_r_X for
        !! the various Richardson iterations
        !! The aim is to  halve the step size, which is given  by dx(n) = 1/(n-1)
        !! or, inverting: n(dx) = 1 + 1/dx.
        !! This yields n(dx/2)/n(dx) = (2+dx)/(1+dx) = (2n(dx)-1)/n(dx)
        !! The recursion formula is therefore: n(dx/2) = 2n(dx) - 1
        !subroutine calc_n_r_X(ir,n_r_X)
            !use X_vars, only: min_n_r_X
            
            !! input / output
            !integer, intent(in) :: ir
            !integer, intent(inout) :: n_r_X
            
            !write(*,*) '!!!! TEMPORARILY DISABLED BECAUSE RICH. EXT. NOT &
                !&WORKING PROPERLY !!!'
            
            !n_r_X = min_n_r_X*ir
            !!!!if (ir.eq.1) then
                !!!!n_r_X = min_n_r_X
            !!!!else
                !!!!n_r_X = 2 * n_r_X - 1
            !!!!end if
            !call writo(trim(i2str(n_r_X))//' normal points for this level...')
        !end subroutine calc_n_r_X
        
        !! calculates  the  coefficients of  the  Eigenvalues  in the  Richardson
        !! extrapolation
        !! This is done using the recursive formula
        !!   X_val_rich(ir,ir2,:) = X_val_rich(ir,ir2-1,:) +  1/(2^(2ir2) - 1) * 
        !!       (X_val_rich(ir,ir2-1,:) - X_val_rich(ir-1,ir2-1,:)),
        !! as described in [ADD REF]
        !! [MPI] All ranks
        !integer function calc_rich_ex(ir,X_val,X_val_rich,done_richard,&
                !&use_guess_for_next_level) result(ierr)
            !use num_vars, only: tol_r
            
            !character(*), parameter :: rout_name = 'calc_rich_ex'
            
            !! input / output
            !integer, intent(inout) :: ir                                        ! level of Richardson extrapolation (starting at 1)
            !complex(dp), intent(in) :: X_val(:)                                 ! EV for this Richardson level
            !complex(dp), intent(inout) :: X_val_rich(:,:,:)                     ! extrapolated coefficients
            !logical, intent(inout) :: done_richard                              ! if Richardson loop has converged sufficiently
            !logical, intent(inout) :: use_guess_for_next_level                  ! if a guessed is used for next Richardson level
            
            !! local variables
            !character(len=max_str_ln) :: err_msg                                ! error message
            !integer :: ir2                                                      ! counter
            !complex(dp), allocatable :: corr(:)                                 ! correction
            !real(dp) :: max_corr                                                ! maximum of maximum of correction
            !real(dp), save :: prev_max_corr                                     ! max_corr at previous Richardson level
            !integer :: loc_max_corr(1)                                          ! location of maximum of correction
            
            !! initialize ierr
            !ierr = 0
            !write(*,*) '!!! FOR RICHARDSON EXTRAPOLATION YOU HAVE TO FIND &
                !&CONSTANT EIGENVALUES !!!'
            
            !! tests
            !if (size(X_val_rich,1).ne.size(X_val_rich,2) .or. &
                !&ir.gt.size(X_val_rich,1)) then
                !ierr = 1
                !err_msg = 'X_val_rich has to have correct dimensions'
                !CHCKERR(err_msg)
            !end if
            
            !! allocate correction
            !allocate(corr(size(X_val))); corr = 1.0E15
            
            !! do calculations if ir > 1
            !X_val_rich(ir,1,1:size(X_val)) = X_val
            !do ir2 = 2,ir
                !corr = 1./(2**(2*ir2)-1.) * &
                    !&(X_val_rich(ir,ir2-1,1:size(X_val)) - &
                    !&X_val_rich(ir-1,ir2-1,1:size(X_val)))
                !X_val_rich(ir,ir2,1:size(X_val)) = &
                    !&X_val_rich(ir,ir2-1,1:size(X_val)) + corr
            !end do
            
            !if (ir.gt.1) then                                                   ! only do this if in Richardson level higher than 1
                !! get maximum and location of maximum for relative correction
                !max_corr = maxval(abs(corr/X_val_rich(ir,ir,1:size(X_val))))
                !loc_max_corr = maxloc(abs(corr/X_val_rich(ir,ir,1:size(X_val))))
                !call writo('maximum relative error: '//trim(r2strt(max_corr))//&
                    !&' for Eigenvalue '//trim(i2str(loc_max_corr(1))))
                
                !! check whether tolerance has been  reached for last value ir2 =
                !! ir
                !if (maxval(abs(corr/X_val_rich(ir,ir,1:size(X_val))))&
                    !&.lt.tol_r) then
                    !done_richard = .true.
                    !call writo('tolerance '//trim(r2strt(tol_r))//&
                        !&' reached after '//trim(i2str(ir))//' iterations')
                !else
                    !call writo('tolerance '//trim(r2strt(tol_r))//' not yet &
                        !&reached')
                !end if
                
                !! determine use_guess_for_next_level
                !use_guess_for_next_level = max_corr.lt.prev_max_corr            ! set guess for next level to true if max_corr is decreasing
                !prev_max_corr = max_corr                                        ! save max_corr for next level
            !else
                !use_guess_for_next_level = .true.                               ! for first Richardson level, set guess for next level to true
            !end if
            
            !!!! TEMPORARILY !!!
            !call writo('!!!!! DO NOT USE GUESS FOR MODIFIED RICHARDSON !!!')
            !use_guess_for_next_level = .false.
            
            !! check for convergence
            !if (.not.done_richard) then
                !if (ir.lt.max_it_r) then                                        ! not yet at maximum Richardson iteration
                    !ir = ir + 1
                !else                                                            ! maximum nr. of Richardson iterations reached
                    !call writo('maximum number of Richardson iterations &
                        !&reached')
                    !done_richard = .true.
                !end if
            !end if
        !end function calc_rich_ex
    !end function run_driver_for_alpha
    
    !! Adapt  X  variables  angularly  to  a  field-aligned  grid,  depending  on
    !! equilibrium type:
    !!   X: J_exp_ang_par_F, U_i, DU_i, PV_i, KV_i
    !! To save memory, the original quantities adapted are deallocated.
    !! Note that this, by definition, does not affect and thus doesn not apply to
    !! flux functions.
    !integer function adapt_X_to_B(grid_eq,grid_eq_B,X,X_B) &
        !&result(ierr)
        !use num_vars, only: eq_style
        !use HELENA, only: interp_HEL_on_grid
        !use X_vars, only: create_X
        
        !character(*), parameter :: rout_name = 'adapt_X_to_B'
        
        !! input / output
        !type(grid_type), intent(in) :: grid_eq, grid_eq_B                       ! general and field-aligned equilibrium grid
        !type(X_type), intent(inout) :: X                                        ! general perturbation variables
        !type(X_type), intent(inout) :: X_B                                      ! field-aligned perturbation variables
        
        !! local variables
        !character(len=max_str_ln) :: err_msg                                    ! error message
        
        !! initialize ierr
        !ierr = 0
        
        !! create new perturbation variables
        !!!!call create_X(grid_eq_B,X_B)
        !write(*,*) 'NOT CREATING X !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        
        !! choose which equilibrium style is being used:
        !!   1:  VMEC
        !!   2:  HELENA
        !select case (eq_style)
            !case (1)                                                            ! VMEC
                !! no conversion necessary: already in field-aligned grid
                !X_B%J_exp_ang_par_F = X%J_exp_ang_par_F
                !X_B%U_0 = X%U_0
                !X_B%U_1 = X%U_1
                !X_B%DU_0 = X%DU_0
                !X_B%DU_1 = X%DU_1
                !X_B%PV_0 = X%PV_0
                !X_B%PV_1 = X%PV_1
                !X_B%PV_2 = X%PV_2
                !X_B%KV_0 = X%KV_0
                !X_B%KV_1 = X%KV_1
                !X_B%KV_2 = X%KV_2
                !X_B%vac_res = X%vac_res
            !case (2)                                                            ! HELENA
                !! call HELENA grid interpolation
                !ierr = interp_HEL_on_grid(grid_eq,grid_eq_B,X,X_B,&
                    !&grid_name='field-aligned grid')
                !CHCKERR('')
            !case default
                !err_msg = 'No equilibrium style associated with '//&
                    !&trim(i2str(eq_style))
                !ierr = 1
                !CHCKERR(err_msg)
        !end select
    !end function adapt_X_to_B
end module driver_PERT
