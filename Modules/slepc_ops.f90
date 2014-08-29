!------------------------------------------------------------------------------!
!   Operations that use slepc (and petsc) routines
!------------------------------------------------------------------------------!
module slepc_ops
#include <PB3D_macros.h>
#include <finclude/slepcepsdef.h>
!#include <finclude/petscsys.h>
    use slepceps
    use num_vars, only: iu, mu_0, dp, max_str_ln
    use output_ops, only: lvl_ud, writo, print_ar_2, print_GP_2D
    use str_ops, only: r2strt, r2str, i2str

    implicit none
    private
    public solve_EV_system_slepc
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in [ADD REF] and solves them using the slepc suite
    
    integer function solve_EV_system_slepc(A_terms,B_terms,A_info,B_info) &
        &result(ierr)
        use X_vars, only: PV0, PV1, PV2, KV0, KV1, KV2, X_vec, X_val, m_X, &
            &n_r_X, n_X, group_n_r_X
        use num_vars, only: MPI_comm_groups, group_n_procs, n_sol_requested, &
            &group_rank, min_r_X, max_r_X
        use file_ops, only: opt_args
        
        character(*), parameter :: rout_name = 'solve_EV_system_slepc'
        
        ! input / output
        PetscScalar, intent(inout), allocatable, optional :: A_terms(:,:,:,:)   ! termj_int of matrix A from previous Richardson loop
        PetscScalar, intent(inout), allocatable, optional :: B_terms(:,:,:,:)   ! termj_int of matrix B from previous Richardson loop
        PetscReal, intent(inout), optional :: A_info(2)                         ! info about A_terms: min of r and step_size
        PetscReal, intent(inout), optional :: B_info(2)                         ! info about B_terms: min of r and step_size
        
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
        PetscInt :: n_m_X                                                       ! nr. of poloidal modes
        Mat :: A                                                                ! matrix A in EV problem A X = lambda B X
        Mat :: B                                                                ! matrix B in EV problem A X = lambda B X
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        
        ! solution
        Vec :: sol_vec                                                          ! solution EV parallel vector
        PetscInt :: n_sol                                                       ! how many solutions can be requested (normally n_sol_requested)
        PetscInt :: one = 1                                                     ! one
        
        ! other variables
        PetscInt :: id                                                          ! counter
        PetscInt :: max_n_EV                                                    ! nr. of EV's saved, up to n_sol_requested
        PetscReal :: step_size                                                  ! step size of perturbation grid
        
        !! for tests
        !Mat :: A_t                                                              ! Hermitian transpose of A
        !Mat :: B_t                                                              ! Hermitian transpose of B
        !PetscScalar :: one = 1.0                                                ! one
        
        ! initialize slepc
        call writo('initialize slepc...')
        call lvl_ud(1)
        
        ! use MPI_Comm_groups for PETSC_COMM_WORLD
        PETSC_COMM_WORLD = MPI_Comm_groups
        if (group_n_procs.gt.n_r_X) then                                        ! too many processors
            call writo('WARNING: using too many processors per group: '&
                &//trim(i2str(group_n_procs))//', while beyond '//&
                &trim(i2str(n_r_X))//' does not bring improvement')
            call writo(' -> consider setting "n_procs_per_alpha" to lower &
                &values, increasing n_r_X or increasing number of field &
                &lines n_alpha')
        end if
        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)                         ! initialize slepc
        CHCKERR('slepc failed to initialize')
        
        ! output
        call writo('slepc started with '//trim(i2str(group_n_procs))&
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
        n_m_X = size(m_X)                                                       ! number of poloidal modes (= size of one block)
        step_size = (max_r_X-min_r_X)/(n_r_X - 1.0)
        
        ! X_vec and X_val
        allocate(X_vec(1:n_m_X,1:group_n_r_X,1:n_sol_requested))
        X_vec = 0.0_dp
        allocate(X_val(1:n_sol_requested))
        X_val = 0.0_dp
        
        ! set up the matrix
        call writo('set up matrices...')
        
        call lvl_ud(1)
        
        ! create a  matrix A and B  with the appropriate number  of preallocated
        ! entries
        ! the numbers of non-zeros in the diagonal and off-diagonal parts
        allocate(d_nz(group_n_r_X*n_m_X)); d_nz = 0
        d_nz(1:n_m_X) = 2*n_m_X
        d_nz(n_m_X+1:(group_n_r_X-1)*n_m_X) = 3*n_m_X
        d_nz((group_n_r_X-1)*n_m_X+1:group_n_r_X*n_m_X) = 2*n_m_X
        allocate(o_nz(group_n_r_X*n_m_X)); o_nz = 0
        if(group_rank.ne.0) then
            o_nz(1:n_m_X) = n_m_X
        end if
        if (group_rank.ne.group_n_procs-1) then
            o_nz((group_n_r_X-1)*n_m_X+1:group_n_r_X*n_m_X) = n_m_X
        end if
        ! create matrix A
        call MatCreateAIJ(PETSC_COMM_WORLD,group_n_r_X*n_m_X,group_n_r_X*n_m_X,&
            &n_r_X*n_m_X,n_r_X*n_m_X,PETSC_NULL_INTEGER,d_nz,&
            &PETSC_NULL_INTEGER,o_nz,A,ierr)
        CHCKERR('MatCreateAIJ failed for matrix A')
        ! deallocate d_nz and o_nz
        deallocate(d_nz,o_nz)
        ! fill the matrix A
        ierr = fill_mat(n_X,n_r_X,n_m_X,step_size,PV0,PV1,PV2,A,A_terms,A_info)
        CHCKERR('')
        call writo('Matrix A set up')
        
        ! duplicate A into B
        ! (the  advantage of  letting communication  and calculation  overlap by
        ! calculating matrix B while assembling A is offset by the extra cost of
        ! not being able to use MatDuplicate. This is also easier)
        call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)                   ! B has same structure as A
        CHCKERR('failed to duplicate A into B')
        ! fill the matrix B
        ierr = fill_mat(n_X,n_r_X,n_m_X,step_size,KV0,KV1,KV2,B,B_terms,B_info)
        CHCKERR('')
        call writo('Matrix B set up')
        
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
       ! 
        ! print info
        call writo('storing results for '//trim(i2str(max_n_EV))//' highest &
            &Eigenvalues...')
        
        call lvl_ud(1)
        
        !call MatGetVecs(A,sol_vec,PETSC_NULL_OBJECT,ierr)                       ! get compatible parallel vector to matrix A
        !CHCKERR('MatGetVecs failed')
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,group_n_r_X*n_m_X,&
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
    integer function fill_mat(n_X,n_r_X,n_m_X,step_size,V0,V1,V2,mat,mat_terms,&
        &mat_info) result(ierr)
        use metric_ops, only: jac_F_FD
        use num_vars, only: min_r_X, max_r_X
        use eq_vars, only: n_par, theta, n_r_eq
        use utilities, only: dis2con
        
        character(*), parameter :: rout_name = 'fill_mat'
        
        ! input / output
        PetscInt, intent(in) :: n_X                                             ! toroidal mode number
        PetscInt, intent(in) :: n_r_X                                           ! nr. of normal points for perturbation quantities
        PetscInt, intent(in) :: n_m_X                                           ! nr. of poloidal modes
        PetscScalar, intent(in) :: V0(:,:,:,:), V1(:,:,:,:), V2(:,:,:,:)        ! either PVi or KVi
        Mat, intent(inout) :: mat                                               ! either A or B
        PetscReal, intent(in) :: step_size                                      ! step size
        PetscScalar, intent(inout), allocatable, optional :: mat_terms(:,:,:,:) ! termj_int of matrix mat from previous Richardson loop
        PetscReal, intent(inout), optional :: mat_info(2)                       ! info about mat_terms: min of r and step_size
        
        ! local variables (not to be used in child routines)
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscScalar, allocatable :: loc_block(:,:)                              ! block matrix for 1 normal point
        PetscScalar, allocatable :: term0(:,:,:,:)                              ! holds e^i(k-m) J V0 for (k,m) functions of (par,r)
        PetscScalar, allocatable :: term1(:,:,:,:)                              ! holds e^i(k-m) J V1 for (k,m) functions of (par,r)
        PetscScalar, allocatable :: term2(:,:,:,:)                              ! holds e^i(k-m) J V2 for (k,m) functions of (par,r)
        PetscScalar, allocatable :: term0_int(:,:)                              ! interpolated and integrated version of term0
        PetscScalar, allocatable :: term1_int(:,:)                              ! interpolated and integrated version of term1
        PetscScalar, allocatable :: term2_int(:,:)                              ! interpolated and integrated version of term2
        PetscScalar, allocatable :: term0_int_next(:,:)                         ! interpolated and integrated version of term0, at next point
        PetscScalar, allocatable :: term1_int_next(:,:)                         ! interpolated and integrated version of term1, at next point
        PetscScalar, allocatable :: term2_int_next(:,:)                         ! interpolated and integrated version of term2, at next point
        PetscScalar, allocatable :: exp_theta(:,:,:,:)                          ! exp(i (k-m) theta)
        PetscInt, allocatable :: loc_k(:), loc_m(:)                             ! the locations at which to add the blocks to the matrices
        PetscInt :: id, jd, m, k                                                ! counters
        PetscInt :: n_start, n_end                                              ! start row and end row
        PetscInt :: r_X_start, r_X_end                                          ! start block and end block (= n_start/n_m_X and equiv.)
        
        ! local variables (to be used in child routines)
        PetscReal :: norm_pt                                                    ! current normal point (min_r_X...max_r_X)
        PetscReal :: mat_info_in(2)                                             ! input of mat_info
        PetscScalar, allocatable :: mat_terms_in(:,:,:,:)                       ! input of mat_terms
        PetscBool :: norm_pt_found                                              ! if norm_pt is found in previous Richardson level
        PetscInt :: found_id                                                    ! in which place norm_pt is found
        
        !! for output
        !real(dp), allocatable :: x_plot(:,:)                                    ! x_axis of plot
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (present(mat_terms)) then                                            ! mat terms present
            if (.not.present(mat_info)) then                                    ! mat_info not present
                ierr = 1
                err_msg = 'mat_terms and info_terms have to be provided both'
                CHCKERR(err_msg)
            end if
        else                                                                    ! mat terms not present
            if (present(mat_info)) then                                         ! mat_terms present
                ierr = 1
                err_msg = 'mat_terms and info_terms have to be provided both'
                CHCKERR(err_msg)
            end if
        end if
        
        ! save input mat_terms and mat_info
        if (present(mat_terms)) then
            if (allocated(mat_terms)) then                                      ! check for allocation of mat_terms
                allocate(mat_terms_in(size(mat_terms,1),size(mat_terms,2),&
                    &size(mat_terms,3),size(mat_terms,4)))                      ! allocate mat_terms_in with size of mat_terms
                mat_terms_in = mat_terms                                        ! save mat_terms in mat_terms_in
                mat_info_in = mat_info                                          ! save mat_info_in in mat_info
                deallocate(mat_terms)                                           ! deallocate mat_terms to reallocate
            end if
        end if
        
        ! allocate variables
        allocate(term0(n_par,n_r_eq,n_m_X,n_m_X))
        allocate(term1(n_par,n_r_eq,n_m_X,n_m_X))
        allocate(term2(n_par,n_r_eq,n_m_X,n_m_X))
        allocate(term0_int(n_m_X,n_m_X))
        allocate(term1_int(n_m_X,n_m_X))
        allocate(term2_int(n_m_X,n_m_X))
        allocate(term0_int_next(n_m_X,n_m_X))
        allocate(term1_int_next(n_m_X,n_m_X))
        allocate(term2_int_next(n_m_X,n_m_X))
        allocate(loc_block(n_m_X,n_m_X))
        allocate(loc_k(n_m_X))
        allocate(loc_m(n_m_X))
        allocate(exp_theta(n_par,n_r_eq,n_m_X,n_m_X))
        
        ! set exp_theta for multiple use
        do m = 1,n_m_X
            do k = 1,n_m_X
                exp_theta(:,:,k,m) = exp(PETSC_i*(k-m)*theta)
            end do
        end do
        
        ! set term1, term2 and term3
        do m = 1,n_m_X
            do k = 1,n_m_X
                term0(:,:,k,m) = exp_theta(:,:,k,m) * jac_F_FD(:,:,0,0,0) * &
                    &V0(:,:,k,m)
                term1(:,:,k,m) = exp_theta(:,:,k,m) * jac_F_FD(:,:,0,0,0) * &
                    &V1(:,:,k,m)
                term2(:,:,k,m) = exp_theta(:,:,k,m) * jac_F_FD(:,:,0,0,0) * &
                    &V2(:,:,k,m)
            end do
        end do
        
        ! deallocate variables
        deallocate(exp_theta)
        
        !! plot for debugging
        !allocate(x_plot(1:n_r_eq,1:n_par))
        !do id = 1,n_r_eq
            !x_plot(id,:) = (id-1.0)/(n_r_eq-1)
        !end do
        !write(*,*) 'term0 = '
        !call print_GP_2D('RE term0','',realpart(transpose(term0(:,:,1,1))),x=x_plot)
        !call print_GP_2D('IM term0','',imagpart(transpose(term0(:,:,1,1))),x=x_plot)
        !write(*,*) 'term1 = '
        !call print_GP_2D('RE term1','',realpart(transpose(term1(:,:,1,1))),x=x_plot)
        !call print_GP_2D('IM term1','',imagpart(transpose(term1(:,:,1,1))),x=x_plot)
        !write(*,*) 'term2 = '
        !call print_GP_2D('RE term2','',realpart(transpose(term2(:,:,1,1))),x=x_plot)
        !call print_GP_2D('IM term2','',imagpart(transpose(term2(:,:,1,1))),x=x_plot)
        !deallocate(x_plot)
        
        ! get ownership of the rows
        call MatGetOwnershipRange(mat,n_start,n_end,ierr)                       ! starting and ending row n_start and n_end
        CHCKERR('Couldn''t get ownership range of matrix')
        r_X_start = n_start/n_m_X                                               ! count per block
        r_X_end = n_end/n_m_X                                                   ! count per block
        
        ! get first interpolated point, corresponding to r_X_start
        call dis2con(r_X_start+1,[1,n_r_X],norm_pt,[min_r_X,max_r_X])
        
        ! (re)allocate  mat_terms and set  mat_info with first normal  point and
        ! step size
        if (present(mat_terms)) then
            allocate(mat_terms(n_m_X,n_m_X,r_X_end-r_X_start+1,3))
            mat_info = [norm_pt,step_size]
        end if
        
        ! check if this point is provided by the previous Richardson level
        call find_norm_pt
        ! calculate termj_int_next
        ierr = get_term_int(term0,norm_pt,term0_int_next,0,1)
        CHCKERR('')
        ierr = get_term_int(term1,norm_pt,term1_int_next,1,1)
        CHCKERR('')
        ierr = get_term_int(term2,norm_pt,term2_int_next,2,1)
        CHCKERR('')
        
        ! iterate over all rows of this rank
        do id = r_X_start, r_X_end-1
            ! fill termi_int with termi_int_next of previous step
            term0_int = term0_int_next
            term1_int = term1_int_next
            term2_int = term2_int_next
            
            ! ---------------!
            ! BLOCKS (id,id) !
            ! ---------------!
            loc_block = 0
            ! part ~ V0(i)
            loc_block = loc_block + term0_int
            ! part ~ V2(i)
            loc_block = loc_block + term2_int * 2.0/(step_size*n_X)**2
            
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
                
                ! check if this point is provided by the previous Richardson level
                call find_norm_pt
                ! calculate termj_int_next
                ierr = get_term_int(term0,norm_pt,term0_int_next,0,&
                    &id-r_X_start+2)
                CHCKERR('')
                ierr = get_term_int(term1,norm_pt,term1_int_next,1,&
                    &id-r_X_start+2)
                CHCKERR('')
                ierr = get_term_int(term2,norm_pt,term2_int_next,2,&
                    &id-r_X_start+2)
                CHCKERR('')
                
                loc_block = 0
                ! part ~ V1*(i+1)
                loc_block = loc_block + conjg(transpose(term1_int_next)) * &
                    &PETSC_i/(2*step_size*n_X)
                ! part ~ V1(i)
                loc_block = loc_block + term1_int * PETSC_i/(2*step_size*n_X)
                ! part ~ V2(i+1/2)
                loc_block = loc_block - (term2_int+term2_int_next)/2 * &
                    &1.0/(step_size*n_X)**2
                
                ! add block to the matrix A
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
        deallocate(term0,term1,term2)
        deallocate(term0_int,term1_int,term2_int)
        deallocate(term0_int_next,term1_int_next,term2_int_next)
        deallocate(loc_block)
        deallocate(loc_k,loc_m)
    contains
        ! checks if a termj can be found in the previous Richardson level
        ! uses  norm_pt, mat_info_in,  mat_terms_in, norm_pt_found  and found_id
        ! from parent routine
        subroutine find_norm_pt
            ! local variables
            PetscInt :: mat_size                                                ! size of mat_terms_in
            PetscInt :: id                                                      ! counter
            PetscReal :: tol = 1.E-5                                            ! tolerance
            PetscReal :: norm_pt_terms                                          ! current normal point of mat_terms_in
            
            ! initialize norm_pt_found
            norm_pt_found = .false.
            
            ! check if mat_terms_in is allocated
            if (allocated(mat_terms_in)) then
                ! set mat_size and initialize found_id and norm_pt_terms
                mat_size = size(mat_terms_in,3)
                found_id = -1
                norm_pt_terms = mat_info_in(1)
                
                ! find the column
                do id = 1,mat_size
                    if (abs(norm_pt-norm_pt_terms).lt.tol) then                 ! norm_pt found
                        found_id = id
                        norm_pt_found = .true.
                        return
                    else                                                        ! norm_pt not found
                        norm_pt_terms = norm_pt_terms + mat_info_in(2)          ! increase norm_pt_terms
                    end if
                end do
            end if
        end subroutine find_norm_pt
        
        ! decides  whether  termj has  to  be  interpolated and  integrated,  or
        ! whether it can be taken from the input mat_terms_in
        ! uses mat_terms_in, norm_pt_found and found_id from parent routine
        integer function get_term_int(term,norm_pt,term_int,which_term,&
            &mat_terms_index) result(ierr)
            character(*), parameter :: rout_name = 'get_term_int'
            
            ! input / output
            PetscScalar, intent(in) :: term(n_par,n_r_eq,n_m_X,n_m_X)           ! input pairs of (k,m) in equilibrium table (par,r)
            PetscReal, intent(in) :: norm_pt                                    ! point at which to interpolate
            PetscScalar, intent(inout) :: term_int(n_m_X,n_m_X)                 ! output pairs of (k,m), interpolated and integrated
            PetscInt, intent(in) :: which_term                                  ! which input of mat_terms 0, 1 or 2
            PetscInt, intent(in) :: mat_terms_index                             ! 1..size of terms
            
            ! initialize ierr
            ierr = 0
            
            if (norm_pt_found) then
                term_int = mat_terms_in(:,:,found_id,which_term+1)
            else
                ierr = interp_and_int(term,norm_pt,term_int)
                CHCKERR('')
            end if
            
            ! update mat_terms with first point
            if (present(mat_terms)) then
                mat_terms(:,:,mat_terms_index,which_term+1) = term_int
            end if
        end function get_term_int
        
        ! interpolates  at normal point  norm_pt and integrates in  the parallel
        ! direction
        integer function interp_and_int(term,norm_pt,term_int) result(ierr)
            use utilities, only: calc_interp, calc_int
            use eq_vars, only: min_r_eq, n_r_eq
            use VMEC_vars, only: n_r
            
            character(*), parameter :: rout_name = 'interp_and_int'
            
            ! input / output
            PetscScalar, intent(in) :: term(n_par,n_r_eq,n_m_X,n_m_X)           ! input pairs of (k,m) in equilibrium table (par,r)
            PetscReal, intent(in) :: norm_pt                                    ! point at which to interpolate
            PetscScalar, intent(inout) :: term_int(n_m_X,n_m_X)                 ! output pairs of (k,m), interpolated and integrated
            
            ! local variables
            PetscScalar :: term_interp(n_par,n_m_X,n_m_X)                       ! term at inteprolated position
            PetscScalar :: term_int_full(n_par,n_m_X,n_m_X)                     ! term, interpolated and integrated but not yet evaluated
            PetscReal, allocatable :: theta_interp(:,:,:)                       ! theta at interpolated position
            PetscInt :: k, m                                                    ! counters
            
            ! initialize ierr
            ierr = 0
            
            ! interpolate term at norm_pt
            ierr = calc_interp(term,[1,n_r],term_interp,norm_pt,&
                &r_offset=min_r_eq-1)
            CHCKERR('')
            
            ! interpolate theta at norm_pt
            allocate(theta_interp(n_par,1,1))                                   ! dimension 3 to be able to use the routine from utilities
            ierr = calc_interp(reshape(theta,[n_par,n_r_eq,1,1]),[1,n_r],&
                &theta_interp,norm_pt,r_offset=min_r_eq-1)
            CHCKERR('')
            
            ! integrate term over theta for all pairs (k,m)
            do m = 1,n_m_X
                do k = 1,n_m_X
                    ierr = calc_int(term_interp(:,k,m),theta_interp(:,1,1),&
                        &term_int_full(:,k,m))
                    CHCKERR('')
                end do
            end do
            
            ! evaluate integral at last point
            term_int = term_int_full(n_par,:,:)
        end function interp_and_int
    end function fill_mat
end module
