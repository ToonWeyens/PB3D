!------------------------------------------------------------------------------!
!   Operations that use slepc (and petsc) routines
!------------------------------------------------------------------------------!
module slepc_ops
#include <PB3D_macros.h>
#include <finclude/slepcepsdef.h>
!#include <finclude/petscsys.h>
    use slepceps
    use num_vars, only: iu, dp, max_str_ln
    use output_ops, only: lvl_ud, writo, print_ar_2, print_GP_2D, print_GP_3D
    use str_ops, only: r2strt, r2str, i2str

    implicit none
    private
    public solve_EV_system_slepc
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in [ADD REF] and solves them using the slepc suite
    
    integer function solve_EV_system_slepc(use_guess) result(ierr)
        use X_vars, only: X_vec, X_val, n_r_X, grp_n_r_X, n_m_X,  n_r_X, &
            &grp_n_r_X, PV_int, KV_int, grp_min_r_X
        use num_vars, only: MPI_comm_groups, grp_n_procs, n_sol_requested, &
            &grp_rank, min_r_X, max_r_X
        use file_ops, only: opt_args
        use MPI_ops, only: divide_grid
        use num_vars, only: mu_0
        use eq_vars, only: B_0, R_0, rho_0
        
        character(*), parameter :: rout_name = 'solve_EV_system_slepc'
        
        ! input / output
        PetscBool, intent(in) :: use_guess                                      ! whether to use a guess or not
        
        ! local variables
        ! petsc / MPI variables
        EPS :: solver                                                           ! EV solver
        PetscInt :: n_it                                                        ! nr. of iterations
        PetscInt :: n_conv                                                      ! nr. of converged solutions
        PetscInt :: n_ev, ncv, mpd                                              ! nr. of requested EV, max. dim of subspace and for projected problem
        EPSType :: EPS_type                                                     ! string with name of type of EPS solver
        PetscReal :: error                                                      ! error of EPS solver
        PetscReal :: tol                                                        ! tolerance of EPS solver
        PetscInt ::  max_it                                                     ! maximum number of iterations
        character(len=max_str_ln) :: option_name                                ! options
        PetscBool :: flg                                                        ! flag to catch options
        
        ! perturbation quantities
        Mat :: A                                                                ! matrix A in EV problem A X = lambda B X
        Mat :: B                                                                ! matrix B in EV problem A X = lambda B X
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        
        ! guess   
        Vec, allocatable :: guess_vec(:)                                        ! guess to solution EV parallel vector
        PetscInt :: n_prev_guess                                                ! nr. elements in previous vector to put in guess
        PetscInt, allocatable :: guess_id(:)                                    ! indices of elements in guess
        PetscInt, save :: guess_start_id = -10                                  ! start of index of previous vector, saved for next iteration
        PetscInt, save :: prev_n_EV                                             ! nr. of solutions of previous vector
        PetscInt :: istat                                                       ! status for commands
        !PetscViewer :: guess_viewer                                             ! viewer for guess object
        
        ! solution
        Vec :: sol_vec                                                          ! solution EV parallel vector
        PetscInt :: n_sol                                                       ! how many solutions can be requested (normally n_sol_requested)
        PetscInt :: one = 1                                                     ! one
        
        ! other variables
        PetscInt :: id, jd, kd                                                  ! counters
        PetscInt :: max_n_EV                                                    ! nr. of EV's saved, up to n_sol_requested
        PetscReal :: step_size                                                  ! step size of perturbation grid
        
        !! for tests
        !Mat :: A_t                                                              ! Hermitian transpose of A
        !Mat :: B_t                                                              ! Hermitian transpose of B
        !PetscScalar :: one = 1.0                                                ! one
        
        ! initialize ierr
        ierr = 0
        
        ! initialize slepc
        call writo('initialize slepc...')
        call lvl_ud(1)
        
        ! divide perturbation grid under group processes
        ierr = divide_grid()
        CHCKERR('')
        
        ! use MPI_Comm_groups for PETSC_COMM_WORLD
        PETSC_COMM_WORLD = MPI_Comm_groups
        if (grp_n_procs.gt.n_r_X) then                                          ! too many processors
            call writo('WARNING: using too many processors per group: '&
                &//trim(i2str(grp_n_procs))//', while beyond '//&
                &trim(i2str(n_r_X))//' does not bring improvement')
            call writo(' -> consider setting "n_procs_per_alpha" to lower &
                &values, increasing n_r_X or increasing number of field &
                &lines n_alpha')
        end if
        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)                         ! initialize slepc
        CHCKERR('slepc failed to initialize')
        
        ! output
        call writo('slepc started with '//trim(i2str(grp_n_procs))&
            &//' processors')
        
        call lvl_ud(-1)
        
        ! checking for complex numbers
        call writo('run tests...')
#if defined(PETSC_USE_COMPLEX)
#else
        err_msg = 'Petsc and slepc have to be configured and compiled &
            &with the option "--with-scalar-type=complex"'
        call SlepcFinalize(ierr)
        ierr = 1
        ERR(err_msg)
#endif
        
        ! catch options so they don't give a Petsc warning
        do id = 1,size(opt_args)
            option_name = opt_args(id)
            if (trim(option_name).ne.'') call PetscOptionsHasName(&
                &PETSC_NULL_CHARACTER,trim(option_name),flg,ierr)
        end do
        
        ! setting variables
        call writo('set variables...')
        ! The step size is either the (normalized) poloidal or toroidal flux
        step_size = (max_r_X-min_r_X)/(n_r_X - 1.0)
        
        ! set up the matrix
        call writo('set up matrices...')
        
        call lvl_ud(1)
        
        ! create a  matrix A and B  with the appropriate number  of preallocated
        ! entries
        ! the numbers of non-zeros in the diagonal and off-diagonal parts
        allocate(d_nz(grp_n_r_X*n_m_X)); d_nz = 0
        allocate(o_nz(grp_n_r_X*n_m_X)); o_nz = 0
        if (grp_n_r_X.eq.1) then                                                ! 1 normal point in the process
            d_nz = n_m_X
            o_nz = 2*n_m_X
        else                                                                    ! more than 1 normal point in the process
            d_nz(1:n_m_X) = 2*n_m_X
            d_nz(n_m_X+1:(grp_n_r_X-1)*n_m_X) = 3*n_m_X
            d_nz((grp_n_r_X-1)*n_m_X+1:grp_n_r_X*n_m_X) = 2*n_m_X
            if(grp_rank.ne.0) then
                o_nz(1:n_m_X) = n_m_X
            end if
            if (grp_rank.ne.grp_n_procs-1) then
                o_nz((grp_n_r_X-1)*n_m_X+1:grp_n_r_X*n_m_X) = n_m_X
            end if
        end if
        ! create matrix A
        call MatCreateAIJ(PETSC_COMM_WORLD,grp_n_r_X*n_m_X,grp_n_r_X*n_m_X,&
            &n_r_X*n_m_X,n_r_X*n_m_X,PETSC_NULL_INTEGER,d_nz,&
            &PETSC_NULL_INTEGER,o_nz,A,ierr)
        CHCKERR('MatCreateAIJ failed for matrix A')
        ! deallocate d_nz and o_nz
        deallocate(d_nz,o_nz)
        ! fill the matrix A
        ierr = fill_mat(step_size,PV_int,A)
        CHCKERR('')
        call writo('matrix A set up')
        
        ! duplicate A into B
        ! (the  advantage of  letting communication  and calculation  overlap by
        ! calculating matrix B while assembling A is offset by the extra cost of
        ! not being able to use MatDuplicate. This is also easier)
        call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)                   ! B has same structure as A
        CHCKERR('failed to duplicate A into B')
        ! fill the matrix B
        ierr = fill_mat(step_size,KV_int,B)
        CHCKERR('')
        call writo('matrix B set up')
        
        call lvl_ud(-1)
        
        !! test if A and B hermitian
        !call MatHermitianTranspose(A,MAT_INITIAL_MATRIX,A_t,ierr)
        !CHCKERR('Hermitian transpose of A failed')
        !call MatAXPY(A_t,-one,A,SAME_NONZERO_PATTERN,ierr)
        !CHCKERR('A-A_t failed')
        !call MatHermitianTranspose(B,MAT_INITIAL_MATRIX,B_t,ierr)
        !CHCKERR('Hermitian transpose of B failed')
        !call MatAXPY(B_t,-one,B,SAME_NONZERO_PATTERN,ierr)
        !CHCKERR('B-B_t failed')
        
        !! visualize the matrices
        !call PetscOptionsSetValue('-draw_pause','-1',ierr)
        !write(*,*) 'A ='
        !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !!call MatView(A,PETSC_VIEWER_DRAW_WORLD,ierr)
        !write(*,*) 'B ='
        !call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !!call MatView(B,PETSC_VIEWER_DRAW_WORLD,ierr)
        !write(*,*) 'A_t ='
        !call MatView(A_t,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !write(*,*) 'B_t ='
        !call MatView(B_t,PETSC_VIEWER_STDOUT_WORLD,ierr)
        
        !! destroy matrices
        !call MatDestroy(A_t,ierr)
        !CHCKERR('Failed to destroy matrix A_t')
        !call MatDestroy(B_t,ierr)
        !CHCKERR('Failed to destroy matrix B_t')
        
        ! solve EV problem
        call writo('solve the EV problem...')
        
        call lvl_ud(1)
        
        !call PetscOptionsSetValue('-eps_view','-1',ierr)
        call EPSCreate(PETSC_COMM_WORLD,solver,ierr)
        CHCKERR('EPSCreate failed')
        
        call EPSSetOperators(solver,A,B,ierr)                                   ! generalized EV problem A X = lambda B X
        !call EPSSetOperators(solver,A,PETSC_NULL_OBJECT,ierr)                   ! normal EV problem A X = lambda X
        CHCKERR('EPSetOperators failed')
        
        !call EPSSetProblemType(solver,EPS_GHEP,ierr)
        !CHCKERR('Set problem type failed')
        
        !call EPSSetType(solver,EPSLAPACK,ierr)
        !CHCKERR('Failed to set type to LAPACK')
        
        ! search for Eigenvalue with largest real value
        call EPSSetWhichEigenpairs(solver,EPS_LARGEST_REAL,ierr)
        CHCKERR('Failed to set which eigenpairs')
        
        ! request n_sol_requested Eigenpairs
        if (n_sol_requested.gt.n_r_X*n_m_X) then
            call writo('WARNING: max. nr. of solutions requested capped to &
                &problem dimension ('//trim(i2str(n_r_X*n_m_X))//')')
            call writo('Increase either min_n_r_X or number of pol. modes or &
                &decrease n_sol_requested')
            n_sol = n_r_X*n_m_X
        else
            n_sol = n_sol_requested
        end if
        call EPSSetDimensions(solver,n_sol,PETSC_DECIDE,&
            &PETSC_DECIDE,ierr)
        
        ! set guess for EV if X_vec is allocated
        if (allocated(X_vec)) then
            if (use_guess) then
                ! allocate guess vectors
                allocate(guess_vec(prev_n_EV))
                
                ! create vecctor guess_vec
                do kd = 1,prev_n_EV
                    call MatGetVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,istat)    ! get compatible parallel vectors to matrix A
                    call VecSetOption(guess_vec(kd),&
                        &VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,istat)          ! ignore negative indices for compacter notation
                end do
                
                ! set the indices and nr. of values to set in guess_vec
                n_prev_guess = size(X_vec,2)
                allocate(guess_id(n_prev_guess*n_m_X))
                do id = 0,n_prev_guess-1
                    do jd = 1,n_m_X
                        guess_id(id*n_m_X+jd) = (guess_start_id-1)*n_m_X*2 + &
                            &id*n_m_X*2 + jd
                    end do
                end do
                
                do kd = 1,prev_n_EV
                    ! set even values of guess_vec
                    call VecSetValues(guess_vec(kd),n_prev_guess*n_m_X,&
                        &guess_id-1,X_vec(:,:,kd),INSERT_VALUES,istat)
                    ! set odd values of guess_vec
                    call VecSetValues(guess_vec(kd),n_prev_guess*n_m_X,&
                        &guess_id-1-n_m_X,X_vec(:,:,kd),INSERT_VALUES,istat)
                        
                    ! assemble the guess vector
                    call VecAssemblyBegin(guess_vec(kd),istat)
                    call VecAssemblyEnd(guess_vec(kd),istat)
                end do
                    
                ! set guess
                call EPSSetInitialSpace(solver,prev_n_EV,guess_vec,istat)
                
                ! destroy guess vector
                call VecDestroy(guess_vec,istat)                                ! destroy guess vector
                
                ! check the status
                if (istat.ne.0) then
                    call writo('WARNING: unable to set guess for Eigenvector, &
                        &working with random guess')
                else
                    call writo('initial guess set up from previous Richardson &
                        &level')
                end if
            end if
            
            ! deallocate X_vec and X_val
            deallocate(X_vec,X_val)
        end if
        
        call writo('run slepc solver')
        
        ! X_vec and X_val
        allocate(X_vec(1:n_m_X,1:grp_n_r_X,1:n_sol_requested))
        X_vec = 0.0_dp
        allocate(X_val(1:n_sol_requested))
        X_val = 0.0_dp
        
        ! set run-time options
        call EPSSetFromOptions(solver,ierr)
        CHCKERR('EPSetFromOptions failed')
        
        ! solve EV problem
        call EPSSolve(solver,ierr) 
        CHCKERR('EPS couldn''t find a solution')
        
        call lvl_ud(-1)
        
        ! output
        call writo('summarize solution...')
        
        call lvl_ud(1)
        
        call EPSGetType(solver,EPS_type,ierr)                                   ! EPS solver type
        CHCKERR('EPSGetType failed')
        call EPSGetTolerances(solver,tol,max_it);                               ! tolerance and max. nr. of iterations
        CHCKERR('EPSGetTolerances failed')
        call writo(trim(EPS_type)//' solver with tolerance '//&
            &trim(r2strt(tol))//' and maximum '//trim(i2str(max_it))//&
            &' iterations')
        call EPSGetConverged(solver,n_conv,ierr)                                ! nr. of converged solutions
        CHCKERR('EPSGetConverged failed')
        call EPSGetIterationNumber(solver,n_it,ierr)                            ! nr. of iterations
        CHCKERR('EPSGetIterationNumber failed')
        call EPSGetDimensions(solver,n_ev,ncv,mpd,ierr)                         ! nr. of requested EV, ...
        CHCKERR('EPSGetDimensions failed')
        call writo(trim(i2str(n_conv))//' converged solutions &
            &after '//trim(i2str(n_it))//' iterations, with '//&
            &trim(i2str(n_ev))//' requested solution')
        
        call lvl_ud(-1)
        
        ! store results
        ! set maximum nr of solutions to be saved
        if (n_sol.gt.n_conv) then
            call writo('WARNING: max. nr. of solutions found only '//&
                &trim(i2str(n_conv)))
            max_n_EV = n_conv
        else
            max_n_EV = n_sol
        end if
        
        ! print info
        call writo('storing results for '//trim(i2str(max_n_EV))//' highest &
            &Eigenvalues...')
        
        call lvl_ud(1)
        
        !call MatGetVecs(A,sol_vec,PETSC_NULL_OBJECT,ierr)                       ! get compatible parallel vector to matrix A
        !CHCKERR('MatGetVecs failed')
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,grp_n_r_X*n_m_X,&
            &n_r_X*n_m_X,PETSC_NULL_SCALAR,sol_vec,ierr)
        CHCKERR('Failed to create MPI vector with arrays')
        
        ! store them
        do id = 1,max_n_EV
            ! get EV solution in vector X_vec
            call VecPlaceArray(sol_vec,X_vec(:,:,id),ierr)
            CHCKERR('Failed to place array')
            call EPSGetEigenpair(solver,id-1,X_val(id),PETSC_NULL_OBJECT,&
                &sol_vec,PETSC_NULL_OBJECT,ierr)                                ! get solution EV vector and value (starts at index 0)
            CHCKERR('EPSGetEigenpair failed')
            call EPSComputeRelativeError(solver,id-1,error,ierr)                ! get error (starts at index 0)
            CHCKERR('EPSComputeRelativeError failed')
            call writo('for solution '//trim(i2str(id))//&
                &'/'//trim(i2str(n_conv))//':')
            
            !! visualize guess
            !if (allocated(X_vec) .and. id.le.prev_n_EV) then
                !call PetscViewerDrawOpen(0,PETSC_NULL_CHARACTER,&
                    !&'guess vector '//trim(i2str(id))//' with '//&
                    !&trim(i2str(n_r_X))//' points',0,0,1000,1000,&
                    !&guess_viewer,istat)
                !call VecView(guess_vec(id),guess_viewer,istat)
                !call VecDestroy(guess_vec(id),istat)                            ! destroy guess vector
            !end if
            !! visualize solution
            !call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,&
                !&'solution '//trim(i2str(id))//' vector with '//&
                !&trim(i2str(n_r_X))//' points',0,&
                !&0,1000,1000,guess_viewer,istat)
            !call VecView(sol_vec,guess_viewer,istat)
            
            ! go back to  physical values from normalization  by multiplying the
            ! Eigenvalues by 1/T0^2 = B_0^2 / (mu_0 rho_0 R_0^2)
            X_val(id) = X_val(id) * B_0**2 / (mu_0*rho_0*R_0**2)
            
            call lvl_ud(1)
            
            ! print output
            call writo('eigenvalue: '//trim(r2strt(realpart(X_val(id))))//&
                &' + '//trim(r2strt(imagpart(X_val(id))))//' i')
            call writo('with rel. error: '//trim(r2str(error)))
            !call EPSPrintSolution(solver,PETSC_NULL_OBJECT,ierr)
            
            !call writo('Eigenvector = ')
            !call VecView(sol_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !CHCKERR('Cannot view vector')
            
            call lvl_ud(-1)
            
            ! reset vector
            call VecResetArray(sol_vec,ierr)
            CHCKERR('Failed to reset array')
        end do
        
        ! set guess_start_id and prev_n_EV
        guess_start_id = grp_min_r_X
        prev_n_EV = max_n_EV
        
        call VecDestroy(sol_vec,ierr)                                           ! destroy compatible vector
        CHCKERR('Failed to destroy sol_vec')
        
        call lvl_ud(-1)
        
        ! destroy and finalize
        call writo('finalize slepc...')
        call EPSDestroy(solver,ierr)
        CHCKERR('Failed to destroy EPS solver')
        call MatDestroy(A,ierr)
        CHCKERR('Failed to destroy matrix A')
        call MatDestroy(B,ierr)
        CHCKERR('Failed to destroy matrix B')
        call SlepcFinalize(ierr)
        CHCKERR('Failed to Finalize slepc')
    end function solve_EV_system_slepc
    
    ! fills a matrix according  to [ADD REF]. Is used for both  the matrix A and
    ! B, corresponding  to the plasma  potential energy and  the (perpendicular)
    ! kinetic energy
    ! The procedure is as follows:
    !   1. tables are set up for (k,m) pairs, in the (par,r) equilibrium grid
    !   2. at the first normal point  in the perturbation grid, belonging to the
    !   current process, the tables are interpolated as a start
    !   3. at  every normal  point in  the perturbation  grid, belonging  to the
    !   current process,  the tables are interpolated  at the next point  in the
    !   corresponding  perturbation  grid.  This  information, as  well  as  the
    !   interpolated value  of the previous  point allow for the  calculation of
    !   every quantity
    ! !!!! THE BOUNDARY CONDITIONS ARE STILL MISSING, SO IT IS TREATED AS HERMITIAN !!!
    integer function fill_mat(step_size,V_int_tab,mat) result(ierr)
        use num_vars, only: min_r_X, max_r_X
        use eq_vars, only: grp_n_r_eq, A_0, flux_p_FD, flux_t_FD
        use utilities, only: dis2con
        use X_vars, only: grp_min_r_X, grp_max_r_X, n_X, n_r_X, n_m_X
        use VMEC_vars, only: use_pol_flux
        
        character(*), parameter :: rout_name = 'fill_mat'
        
        ! input / output
        PetscScalar, intent(in) :: V_int_tab(:,:,:,:)                           ! either PV_int or KV_int tables in grp_n_r_eq normal points
        Mat, intent(inout) :: mat                                               ! either A or B
        PetscReal, intent(in) :: step_size                                      ! step size of current Richardson level
        
        ! local variables (not to be used in child routines)
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscScalar, allocatable :: loc_block(:,:)                              ! (n_m_X x n_m_X) block matrix for 1 normal point
        PetscScalar, allocatable :: V_interp(:,:,:)                             ! interpolated V_int
        PetscScalar, allocatable :: V_interp_next(:,:,:)                        ! interpolated V_int at next normal point
        PetscInt, allocatable :: loc_k(:), loc_m(:)                             ! the locations at which to add the blocks to the matrices
        PetscReal :: step_size_FD                                               ! step size in flux coordinates
        PetscInt :: id, jd                                                      ! counters
        
        ! for tests
        PetscInt :: r_X_start, r_X_end                                          ! start block and end block
        
        ! local variables (to be used in child routines)
        PetscReal :: norm_pt                                                    ! current normal point (min_r_X...max_r_X)
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        allocate(V_interp(n_m_X,n_m_X,3))
        allocate(V_interp_next(n_m_X,n_m_X,3))
        allocate(loc_block(n_m_X,n_m_X))
        allocate(loc_k(n_m_X))
        allocate(loc_m(n_m_X))
        
        !  test  whether  the  matrix   range  coincides  with  grp_min_r_X  and
        ! grp_max_r_X
        call MatGetOwnershipRange(mat,r_X_start,r_X_end,ierr)                   ! starting and ending row r_X_start and r_X_end
        CHCKERR('Couldn''t get ownership range of matrix')
        r_X_start = r_X_start/n_m_X                                             ! count per block
        r_X_end = r_X_end/n_m_X                                                 ! count per block
        if (grp_min_r_X.ne.r_X_start+1) then
            ierr = 1
            err_msg = 'start of matrix in this process does not coincide with &
                &tabulated value'
            CHCKERR(err_msg)
        end if
        if (grp_max_r_X.ne.r_X_end) then
            ierr = 1
            err_msg = 'end of matrix in this process does not coincide with &
                &tabulated value'
            CHCKERR(err_msg)
        end if
        
        ! get first interpolated point, corresponding to r_X_start
        call dis2con(grp_min_r_X,[1,n_r_X],norm_pt,[min_r_X,max_r_X])
        
        ! set up step_size_FD
        if (use_pol_flux) then
            step_size_FD = step_size * (flux_p_FD(grp_n_r_eq,0)-flux_p_FD(1,0)) / &
                &(grp_n_r_eq-1)
        else
            step_size_FD = step_size * (flux_t_FD(grp_n_r_eq,0)-flux_t_FD(1,0)) / &
                &(grp_n_r_eq-1)
        end if
        
        ! get the interpolated terms in V_interp
        ierr = get_V_interp(V_int_tab,norm_pt,V_interp_next)
        CHCKERR('')
        
        ! iterate over all rows of this rank
        do id = grp_min_r_X-1, grp_max_r_X-1                                    ! (indices start with 0 here)
            ! fill V_interp with V_interp_next of previous step
            V_interp = V_interp_next
            
            ! ---------------!
            ! BLOCKS (id,id) !
            ! ---------------!
            loc_block = 0
            ! part ~ V0(i)
            loc_block = loc_block + V_interp(:,:,1)
            ! part ~ V2(i)
            loc_block = loc_block + V_interp(:,:,3) * &
                &2.0/(step_size_FD*n_X*A_0)**2
            
            ! add block to the matrix A
            loc_k = [(jd, jd = 0,n_m_X-1)] + id*n_m_X
            loc_m = loc_k
            call MatSetValues(mat,n_m_X,loc_k,n_m_X,loc_m,loc_block,&
                &INSERT_VALUES,ierr)
            CHCKERR('Couldn''t add values to matrix')
            
            ! -----------------------------------------!
            ! BLOCKS (id,id+1) and hermitian conjugate !
            ! -----------------------------------------!
            if (id.ne.n_r_X-1) then                                             ! don't do the right block for last normal point
                ! calculate next norm_pt by adding  the step size 1/(n_r_X+1) to
                ! the previous norm_pt
                norm_pt = norm_pt + step_size
                
                ! get the interpolated terms in V_interp
                ierr = get_V_interp(V_int_tab,norm_pt,V_interp_next)
                CHCKERR('')
                
                loc_block = 0
                ! part ~ V1*(i+1)
                loc_block = loc_block + conjg(transpose(V_interp_next(:,:,2)))&
                    & * PETSC_i/(2*step_size_FD*n_X*A_0)
                ! part ~ V1(i)
                loc_block = loc_block + V_interp(:,:,2) * &
                    &PETSC_i/(2*step_size_FD*n_X*A_0)
                ! part ~ V2(i+1/2)
                loc_block = loc_block - (V_interp(:,:,3) + &
                    &V_interp_next(:,:,3))/2 * 1.0/(step_size_FD*n_X*A_0)**2
                
                ! add block to the matrix
                loc_k = [(jd, jd = 0,n_m_X-1)] + id*n_m_X
                loc_m = loc_k+n_m_X
                call MatSetValues(mat,n_m_X,loc_k,n_m_X,loc_m,loc_block,&
                    &INSERT_VALUES,ierr)
                CHCKERR('Couldn''t add values to matrix')
                
                ! add Hermitian conjugate
                loc_block = conjg(transpose(loc_block))
                call MatSetValues(mat,n_m_X,loc_m,n_m_X,loc_k,loc_block,&
                    &INSERT_VALUES,ierr)
                CHCKERR('Couldn''t add values to matrix')
            end if
        end do
        
        ! assemble the matrix and view it
        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t begin assembly of matrix')
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t end assembly of matrix')
        
        !call MatSetOption(mat,MAT_HERMITIAN,PETSC_TRUE,ierr)                    ! Hermitian to a first approximation
        !CHCKERR('Coulnd''t set option Hermitian')

        ! deallocate variables
        deallocate(V_interp,V_interp_next)
        deallocate(loc_block)
        deallocate(loc_k,loc_m)
    contains
        ! decides  whether V_interp  has to  be interpolated  from V_int_tab  or
        ! whether it  can be taken  from the  input V_interp_in. Also  saves the
        ! calculated terms in V_interp_out for use in next Richardson iteration
        ! note: uses found_id, norm_pt_found from parent routine
        integer function get_V_interp(V_int_tab,norm_pt,V_interp) result(ierr)
            use utilities, only: con2dis, round_with_tol
            use VMEC_vars, only: n_r_eq
            use eq_vars, only: grp_min_r_eq
            
            character(*), parameter :: rout_name = 'get_V_interp'
            
            ! input / output
            PetscScalar, intent(in) :: V_int_tab(n_m_X,n_m_X,grp_n_r_eq,3)      ! table of (k,m) integrated values at (1..grp_n_r_eq)
            PetscReal, intent(inout) :: norm_pt                                 ! point at which to interpolate
            PetscScalar, intent(inout) :: V_interp(n_m_X,n_m_X,3)               ! output pairs of (k,m), interpolated and integrated at norm_pt
            
            ! local variables
            PetscInt :: id                                                      ! counter
            PetscReal :: norm_pt_dis                                            ! equivalent of norm_pt in table of equilibrium normal points
            PetscInt :: norm_pt_lo, norm_pt_hi                                  ! rounded up or down norm_pt_dis
            PetscInt :: offset                                                  ! offset in the table V_int
            
            ! initialize ierr
            ierr = 0
            
            ! loop over all terms 0,1,2 (with index 1,2,3 in variables)
            do id = 1,3
                ! round norm_pt if necessary
                ierr = round_with_tol(norm_pt,0.0_dp,1.0_dp)
                CHCKERR('')
                
                ! get the table indices between which to interpolate. The tables
                ! are set up in the equilibrium grid, which uses the same normal
                ! coordinate as  the discretization (determined by  VMEC through
                ! use_pol_flux)
                call con2dis(norm_pt,[0._dp,1._dp],norm_pt_dis,[1,n_r_eq])      ! convert norm_pt to the index in the array V_int
                
                ! round up and down and set offset
                norm_pt_lo = floor(norm_pt_dis)
                norm_pt_hi = ceiling(norm_pt_dis)
                offset =  grp_min_r_eq - 1
                
                ! interpolate, also correct if norm_pt_hi = norm_pt_lo
                V_interp(:,:,id) = V_int_tab(:,:,norm_pt_lo-offset,id) + &
                    &(V_int_tab(:,:,norm_pt_hi-offset,id)-&
                    &V_int_tab(:,:,norm_pt_lo-offset,id))*&
                    &(norm_pt_dis-norm_pt_lo)                                   ! because norm_pt_hi - norm_pt_lo = 1
            end do
        end function get_V_interp
    end function fill_mat
end module
