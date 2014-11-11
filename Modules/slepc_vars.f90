!------------------------------------------------------------------------------!
!   Holds variables and routines that use slepc (and petsc) routines
!------------------------------------------------------------------------------!
module slepc_vars
#include <PB3D_macros.h>
#include <finclude/slepcepsdef.h>
!#include <finclude/petscsys.h>
    use slepceps
    use num_vars, only: iu, dp, max_str_ln
    use message_ops, only: lvl_ud, writo, print_ar_2
    use str_ops, only: r2strt, r2str, i2str

    implicit none
    private
    public start_slepc, stop_slepc, setup_matrices, setup_solver, setup_guess, &
        &get_solution, summarize_solution, store_results
    
contains
    !  This  subroutine starts  petsc  and  slepc  with  the correct  number  of
    ! processors
    integer function start_slepc() result(ierr)
        use num_vars, only: MPI_Comm_groups, grp_n_procs
        use X_vars, only: n_r_X
        use file_ops, only: opt_args
        use MPI_ops, only: divide_grid
        
        character(*), parameter :: rout_name = 'start_slepc'
        
        ! local variables
        PetscBool :: flg                                                        ! flag to catch options
        PetscInt :: id                                                          ! counters
        character(len=max_str_ln) :: option_name                                ! options
        
        ! initialize ierr
        ierr = 0
        
        ! user message
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
    end function start_slepc
    
    ! sets up the matrices A and B in the EV problem A X = lambda B X
    integer function setup_matrices(A,B) result(ierr)
        use num_vars, only: grp_n_procs, grp_rank
        use X_vars, only: grp_n_r_X, size_X, PV_int, KV_int, n_r_X
        
        character(*), parameter :: rout_name = 'setup_matrices'
        
        ! input / output
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        
        ! local variables
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        PetscReal, allocatable :: grp_r_eq(:)                                   ! unrounded index in tables V_int
        PetscInt :: offset                                                      ! offset in the table V_int
        
        !! for tests
        !Mat :: A_t                                                              ! Hermitian transpose of A
        !Mat :: B_t                                                              ! Hermitian transpose of B
        !PetscScalar :: one = 1.0                                                ! one
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set  up arrays  grp_r_eq and  offset that  determine how  to fill  the
        ! matrices
        ierr = get_interp_data(grp_r_eq,offset)
        CHCKERR('')
        
        ! create a  matrix A and B  with the appropriate number  of preallocated
        ! entries
        ! the numbers of non-zeros in the diagonal and off-diagonal parts
        allocate(d_nz(grp_n_r_X*size_X)); d_nz = 0
        allocate(o_nz(grp_n_r_X*size_X)); o_nz = 0
        if (grp_n_r_X.eq.1) then                                                ! 1 normal point in the process
            d_nz = size_X
            o_nz = 2*size_X
        else                                                                    ! more than 1 normal point in the process
            d_nz(1:size_X) = 2*size_X
            d_nz(size_X+1:(grp_n_r_X-1)*size_X) = 3*size_X
            d_nz((grp_n_r_X-1)*size_X+1:grp_n_r_X*size_X) = 2*size_X
            if(grp_rank.ne.0) then
                o_nz(1:size_X) = size_X
            end if
            if (grp_rank.ne.grp_n_procs-1) then
                o_nz((grp_n_r_X-1)*size_X+1:grp_n_r_X*size_X) = size_X
            end if
        end if
        
        ! create matrix A
        call MatCreateAIJ(PETSC_COMM_WORLD,grp_n_r_X*size_X,grp_n_r_X*size_X,&
            &n_r_X*size_X,n_r_X*size_X,PETSC_NULL_INTEGER,d_nz,&
            &PETSC_NULL_INTEGER,o_nz,A,ierr)
        CHCKERR('MatCreateAIJ failed for matrix A')
        
        ! deallocate d_nz and o_nz
        deallocate(d_nz,o_nz)
        
        ! fill the matrix A
        ierr = fill_mat(PV_int,A,grp_r_eq,offset)
        CHCKERR('')
        call writo('matrix A set up')
        
        ! duplicate A into B
        ! (the  advantage of  letting communication  and calculation  overlap by
        ! calculating matrix B while assembling A is offset by the extra cost of
        ! not being able to use MatDuplicate. This is also easier)
        call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)                   ! B has same structure as A
        CHCKERR('failed to duplicate A into B')
        ! fill the matrix B
        ierr = fill_mat(KV_int,B,grp_r_eq,offset)
        CHCKERR('')
        call writo('matrix B set up')
        
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
        
        ! deallocate quantities
        deallocate(grp_r_eq)
        
        call lvl_ud(-1)
    contains
        ! fills a  matrix according to  [ADD REF]. Is  used for both  the matrix
        ! A  and  B,  corresponding  to  the plasma  potential  energy  and  the
        ! (perpendicular) kinetic energy
        ! The procedure is as follows:
        !   1. tables  are set  up for  (k,m) (pol. flux)  or (l,n)  pairs (tor.
        !   flux), in the (par,r) equilibrium grid
        !   2. at the first normal point  in the perturbation grid, belonging to
        !   the current process, the tables are interpolated as a start
        !   3. at every normal point in  the perturbation grid, belonging to the
        !   current process,  the tables are  interpolated at the next  point in
        !   the  corresponding  perturbation  grid. This  information,  as  well
        !   as  the interpolated  value  of  the previous  point  allow for  the
        !   calculation of every quantity
        ! Note: the factors i/n or i/m are already included in V_int_tab
        ! !!!! THE BOUNDARY CONDITIONS ARE STILL MISSING !!!!!!
        integer function fill_mat(V_int_tab,mat,grp_r_eq,offset) result(ierr)
            use num_vars, only: min_r_X, max_r_X
            use eq_vars, only: max_flux
            use utilities, only: dis2con
            use X_vars, only: grp_min_r_X, grp_max_r_X, n_r_X, size_X
            
            character(*), parameter :: rout_name = 'fill_mat'
            
            ! input / output
            PetscScalar, intent(in) :: V_int_tab(:,:,:,:)                       ! either PV_int or KV_int tables in grp_n_r_eq normal points
            Mat, intent(inout) :: mat                                           ! either A or B
            PetscReal, intent(in) :: grp_r_eq(:)                                ! unrounded index in tables V_int
            PetscInt, intent(in) :: offset                                      ! offset in the table V_int
            
            ! local variables (not to be used in child routines)
            character(len=max_str_ln) :: err_msg                                ! error message
            PetscScalar, allocatable :: loc_block(:,:)                          ! (size_X x size_X) block matrix for 1 normal point
            PetscScalar, allocatable :: V_interp(:,:,:)                         ! interpolated V_int
            PetscScalar, allocatable :: V_interp_next(:,:,:)                    ! interpolated V_int at next normal point
            PetscScalar, allocatable :: V_interp_next_OLD(:,:,:)                ! interpolated V_int at next normal point
            PetscInt, allocatable :: loc_k(:), loc_m(:)                         ! the locations at which to add the blocks to the matrices
            PetscReal :: step_size                                              ! step size in flux coordinates
            PetscInt :: id, jd                                                  ! counters
            PetscInt :: grp_r_eq_lo, grp_r_eq_hi                                ! rounded up or down grp_r_eq
            
            ! for tests
            PetscInt :: r_X_start, r_X_end                                      ! start block and end block
            
            ! initialize ierr
            ierr = 0
            
            ! allocate variables
            allocate(V_interp(size_X,size_X,3))
            allocate(V_interp_next(size_X,size_X,3))
            allocate(V_interp_next_OLD(size_X,size_X,3))
            allocate(loc_block(size_X,size_X))
            allocate(loc_k(size_X))
            allocate(loc_m(size_X))
            
            !  test  whether the  matrix  range coincides  with grp_min_r_X  and
            ! grp_max_r_X
            call MatGetOwnershipRange(mat,r_X_start,r_X_end,ierr)               ! starting and ending row r_X_start and r_X_end
            err_msg = 'Couldn''t get ownership range of matrix'
            CHCKERR(err_msg)
            r_X_start = r_X_start/size_X                                        ! count per block
            r_X_end = r_X_end/size_X                                            ! count per block
            if (grp_min_r_X.ne.r_X_start+1) then
                ierr = 1
                err_msg = 'start of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            if (grp_max_r_X.ne.r_X_end) then
                ierr = 1
                err_msg = 'end of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            
            ! set up step_size
            step_size = (max_r_X-min_r_X)/(n_r_X-1.0) * max_flux
            
            ! get  the   interpolated  terms   in  V_interp  (also   correct  if
            ! grp_r_eq_hi = grp_r_eq_lo)
            grp_r_eq_lo = floor(grp_r_eq(1))
            grp_r_eq_hi = ceiling(grp_r_eq(1))
            do jd = 1,3
                V_interp_next(:,:,jd) = V_int_tab(:,:,grp_r_eq_lo-offset,jd) + &
                    &(V_int_tab(:,:,grp_r_eq_hi-offset,jd)-&
                    &V_int_tab(:,:,grp_r_eq_lo-offset,jd))*&
                    &(grp_r_eq(1)-grp_r_eq_lo)                                  ! because grp_r_eq_hi - grp_r_eq_lo = 1
            end do
            
            ! iterate over all rows of this rank
            do id = grp_min_r_X-1, grp_max_r_X-1                                ! (indices start with 0 here)
                ! fill V_interp with V_interp_next of previous step
                V_interp = V_interp_next
                
                ! ---------------!
                ! BLOCKS (id,id) !
                ! ---------------!
                loc_block = 0
                ! part ~ V0(i)
                loc_block = loc_block + V_interp(:,:,1)
                ! part ~ V2(i)
                loc_block = loc_block - V_interp(:,:,3) * &
                    &2.0_dp/(step_size)**2
                
                ! add block to the matrix A
                loc_k = [(jd, jd = 0,size_X-1)] + id*size_X
                loc_m = loc_k
                call MatSetValues(mat,size_X,loc_k,size_X,loc_m,loc_block,&
                    &INSERT_VALUES,ierr)
                CHCKERR('Couldn''t add values to matrix')
                
                ! -----------------------------------------!
                ! BLOCKS (id,id+1) and hermitian conjugate !
                ! -----------------------------------------!
                if (id.ne.n_r_X-1) then                                         ! don't do the right block for last normal point
                    ! get the interpolated terms in V_interp
                    grp_r_eq_lo = floor(grp_r_eq(id-grp_min_r_X+3))
                    grp_r_eq_hi = ceiling(grp_r_eq(id-grp_min_r_X+3))
                    do jd = 1,3
                        V_interp_next(:,:,jd) = &
                            &V_int_tab(:,:,grp_r_eq_lo-offset,jd) + &
                            &(V_int_tab(:,:,grp_r_eq_hi-offset,jd)-&
                            &V_int_tab(:,:,grp_r_eq_lo-offset,jd))*&
                            &(grp_r_eq(id-grp_min_r_X+3)-grp_r_eq_lo)           ! because grp_r_eq_hi - grp_r_eq_lo = 1
                    end do
                    
                    loc_block = 0
                    ! part ~ V1*(i+1)
                    loc_block = loc_block + 1.0_dp/(2*step_size) * &
                        &conjg(transpose(V_interp_next(:,:,2)))
                    ! part ~ V1(i)
                    loc_block = loc_block + 1.0_dp/(2*step_size) * &
                        &V_interp(:,:,2)
                    ! part ~ V2(i+1/2)
                    loc_block = loc_block + (V_interp(:,:,3) + &
                        &V_interp_next(:,:,3))/2 * 1.0/(step_size)**2
                    
                    ! add block to the matrix
                    loc_k = [(jd, jd = 0,size_X-1)] + id*size_X
                    loc_m = loc_k+size_X
                    call MatSetValues(mat,size_X,loc_k,size_X,loc_m,loc_block,&
                        &INSERT_VALUES,ierr)
                    err_msg = 'Couldn''t add values to matrix'
                    CHCKERR(err_msg)
                    
                    ! add Hermitian conjugate
                    loc_block = conjg(transpose(loc_block))
                    call MatSetValues(mat,size_X,loc_m,size_X,loc_k,loc_block,&
                        &INSERT_VALUES,ierr)
                    err_msg = 'Couldn''t add values to matrix'
                    CHCKERR(err_msg)
                end if
            end do
            
            ! assemble the matrix and view it
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t begin assembly of matrix')
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t end assembly of matrix')
            
            !call MatSetOption(mat,MAT_HERMITIAN,PETSC_TRUE,ierr)                ! Hermitian to a first approximation
            !CHCKERR('Coulnd''t set option Hermitian')
            
            ! deallocate variables
            deallocate(V_interp,V_interp_next)
            deallocate(loc_block)
            deallocate(loc_k,loc_m)
        end function fill_mat
        
        ! calculates the variables grp_r_eq and  offset, which are later used to
        ! calculate V_interp from the tabulated  values in the equilibrium grid,
        ! the oordinate of which is determined by the variable VMEC_use_pol_flux
        integer function get_interp_data(grp_r_eq,offset) &
                &result(ierr)
            use utilities, only: con2dis, round_with_tol, interp_fun_1D
            use VMEC_vars, only: n_r_eq, VMEC_use_pol_flux
            use eq_vars, only: grp_min_r_eq, max_flux, max_flux_VMEC, &
                &flux_p_FD, flux_t_FD
            use X_vars, only: grp_r_X
            use num_vars, only: use_pol_flux
            
            character(*), parameter :: rout_name = 'get_interp_data'
            
            ! input / output
            PetscReal, allocatable, intent(inout) :: grp_r_eq(:)                ! unrounded index in tables V_int
            PetscInt, intent(inout) :: offset                                   ! offset in the table V_int
            
            ! local variables
            PetscReal :: r_X_loc                                                ! local copy of grp_r_X
            PetscReal :: grp_r_eq_eq_con                                        ! equivalent of grp_r_eq in table of equilibrium normal points
            PetscReal, pointer :: flux(:), flux_VMEC(:)                         ! either pol. or tor. flux
            PetscInt :: kd                                                      ! counter
            
            ! initialize ierr
            ierr = 0
            
            ! set up flux and flux_VMEC
            if (use_pol_flux) then
                flux => flux_p_FD(:,0)
            else
                flux => flux_t_FD(:,0)
            end if
            if (VMEC_use_pol_flux) then
                flux_VMEC => flux_p_FD(:,0)
            else
                flux_VMEC => flux_t_FD(:,0)
            end if
            
            ! allocate grp_r_eq_lo and hi and offset
            allocate(grp_r_eq(size(grp_r_X)))
            
            ! loop over all normal points
            do kd = 1,size(grp_r_X)
                ! set local copy of grp_r_X
                r_X_loc = grp_r_X(kd)
                
                ! round r_X_loc if necessary
                ierr = round_with_tol(r_X_loc,0.0_dp,1.0_dp)
                CHCKERR('')
                
                ! get the table indices between which to interpolate. The tables
                ! are set up in the equilibrium grid, which not necessarily uses
                ! the same  normal coordinate as the  discretization. Therefore,
                ! conversion is necessary
                ! 1. continuous equilibrium grid (0..1)
                ierr = interp_fun_1D(grp_r_eq_eq_con,flux_VMEC/max_flux_VMEC,&
                    &r_X_loc,flux/max_flux)
                CHCKERR('')
                ! 2. discrete equilibrium grid, unrounded
                call con2dis(grp_r_eq_eq_con,[0._dp,1._dp],grp_r_eq(kd),&
                    &[1,n_r_eq])
                
                ! set offset
                offset =  grp_min_r_eq - 1
            end do
        end function get_interp_data
    end function setup_matrices
    
    ! sets up EV solver
    integer function setup_solver(A,B,solver) result(ierr)
        use num_vars, only: n_sol_requested
        use X_vars, only: n_r_X, size_X
        
        character(*), parameter :: rout_name = 'setup_solver'
        
        ! input / output
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! local variables
        PetscInt :: n_sol                                                       ! how many solutions can be requested (normally n_sol_requested)
        
        ! initialize ierr
        ierr = 0
        
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
        if (n_sol_requested.gt.n_r_X*size_X) then
            call writo('WARNING: max. nr. of solutions requested capped to &
                &problem dimension ('//trim(i2str(n_r_X*size_X))//')')
            call writo('Increase either min_n_r_X or number of pol. modes or &
                &decrease n_sol_requested')
            n_sol = n_r_X*size_X
        else
            n_sol = n_sol_requested
        end if
        call EPSSetDimensions(solver,n_sol,PETSC_DECIDE,&
            &PETSC_DECIDE,ierr)
        
        call lvl_ud(-1)
    end function setup_solver
    
    ! sets up guess in solver
    subroutine setup_guess(A,solver,start_id,prev_n_EV)
        use X_vars, only: X_vec, size_X
        
        ! input / output
        Mat, intent(in) :: A                                                    ! matrix A (or B)
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(in) :: start_id                                        ! start of index of previous vector, saved for next iteration
        PetscInt, intent(in) :: prev_n_EV                                       ! nr. of solutions of previous vector
        
        ! local variables
        Vec, allocatable :: guess_vec(:)                                        ! guess to solution EV parallel vector
        PetscInt, allocatable :: guess_id(:)                                    ! indices of elements in guess
        PetscInt :: id, jd, kd                                                  ! counters
        PetscInt :: istat                                                       ! status for commands
        PetscInt :: n_prev_guess                                                ! nr. elements in previous vector to put in guess
        !PetscViewer :: guess_viewer                                             ! viewer for guess object
        
        call lvl_ud(1)
        
        ! set guess for EV if X_vec is allocated and use_guess is true
        if (allocated(X_vec)) then
            ! allocate guess vectors
            allocate(guess_vec(prev_n_EV))
            
            ! create vecctor guess_vec
            do kd = 1,prev_n_EV
                call MatGetVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,istat)        ! get compatible parallel vectors to matrix A
                call VecSetOption(guess_vec(kd),&
                    &VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,istat)              ! ignore negative indices for compacter notation
            end do
            
            ! set the indices and nr. of values to set in guess_vec
            n_prev_guess = size(X_vec,2)
            allocate(guess_id(n_prev_guess*size_X))
            do id = 0,n_prev_guess-1
                do jd = 1,size_X
                    guess_id(id*size_X+jd) = (start_id-1)*size_X*2 + &
                        &id*size_X*2 + jd
                end do
            end do
            
            do kd = 1,prev_n_EV
                ! set even values of guess_vec
                call VecSetValues(guess_vec(kd),n_prev_guess*size_X,&
                    &guess_id-1,X_vec(:,:,kd),INSERT_VALUES,istat)
                ! set odd values of guess_vec
                call VecSetValues(guess_vec(kd),n_prev_guess*size_X,&
                    &guess_id-1-size_X,X_vec(:,:,kd),INSERT_VALUES,istat)
                    
                ! assemble the guess vector
                call VecAssemblyBegin(guess_vec(kd),istat)
                call VecAssemblyEnd(guess_vec(kd),istat)
            end do
                
            ! set guess
            call EPSSetInitialSpace(solver,prev_n_EV,guess_vec,istat)
            
            !! visualize guess
            !if (allocated(X_vec) .and. id.le.prev_n_EV) then
                !call PetscViewerDrawOpen(0,PETSC_NULL_CHARACTER,&
                    !&'guess vector '//trim(i2str(id))//' with '//&
                    !&trim(i2str(n_r_X))//' points',0,0,1000,1000,&
                    !&guess_viewer,istat)
                !call VecView(guess_vec(id),guess_viewer,istat)
                !call VecDestroy(guess_vec(id),istat)                            ! destroy guess vector
            !end if
            
            ! destroy guess vector
            call VecDestroy(guess_vec,istat)                                    ! destroy guess vector
            
            ! check the status
            if (istat.ne.0) then
                call writo('WARNING: unable to set guess for Eigenvector, &
                    &working with random guess')
            else
                call writo('initial guess set up from previous Richardson &
                    &level')
            end if
        end if
        
        call lvl_ud(-1)
    end subroutine setup_guess
    
    ! get the solution vectors
    integer function get_solution(solver) result(ierr)
        use num_vars, only: n_sol_requested
        use X_vars, only: X_vec, X_val, grp_n_r_X, size_X
        
        character(*), parameter :: rout_name = 'get_solution'
        
        ! input / output
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        if (allocated(X_vec)) then
            ! deallocate X_vec and X_val
            deallocate(X_vec,X_val)
        end if
        
        ! X_vec and X_val
        allocate(X_vec(1:size_X,1:grp_n_r_X,1:n_sol_requested))
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
    end function get_solution
    
    ! summarizes solution
    integer function summarize_solution(solver,max_n_EV) result (ierr)
        use num_vars, only: n_sol_requested
        
        character(*), parameter :: rout_name = 'summarize_solution'
        
        ! input / output
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(inout) :: max_n_EV                                     ! nr. of EV's saved, up to n_conv
        
        ! local variables
        PetscInt :: n_it                                                        ! nr. of iterations
        PetscInt :: n_conv                                                      ! nr. of converged solutions
        PetscInt :: n_ev, ncv, mpd                                              ! nr. of requested EV, max. dim of subspace and for projected problem
        EPSType :: EPS_type                                                     ! string with name of type of EPS solver
        PetscReal :: tol                                                        ! tolerance of EPS solver
        PetscInt ::  max_it                                                     ! maximum number of iterations
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! user output
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
        
        ! set maximum nr of solutions to be saved
        if (n_sol_requested.gt.n_conv) then
            call writo('WARNING: max. nr. of solutions found only '//&
                &trim(i2str(n_conv)))
            max_n_EV = n_conv
        else
            max_n_EV = n_sol_requested
        end if
        
        call lvl_ud(-1)
    end function summarize_solution
    
    ! stores the results
    integer function store_results(solver,max_n_EV) result(ierr)
        use X_vars, only: size_X, grp_n_r_X, n_r_X, X_vec, X_val
        use eq_vars, only: B_0, R_0, rho_0
        use num_vars, only: mu_0
        
        character(*), parameter :: rout_name = 'store_results'
        
        ! input / output
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(inout) :: max_n_EV                                     ! nr. of EV's saved, up to n_conv
        
        ! local variables
        PetscInt :: id                                                          ! counter
        Vec :: sol_vec                                                          ! solution EV parallel vector
        PetscReal :: error                                                      ! error of EPS solver
        PetscInt :: one = 1                                                     ! one
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! create solution vector
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,grp_n_r_X*size_X,&
            &n_r_X*size_X,PETSC_NULL_SCALAR,sol_vec,ierr)
        CHCKERR('Failed to create MPI vector with arrays')
        
        ! store them
        do id = 1,max_n_EV
            ! get EV solution in vector X_vec
            call VecPlaceArray(sol_vec,X_vec(:,1:grp_n_r_X,id),ierr)            ! place array sol_vec in X_vec
            CHCKERR('Failed to place array')
            call EPSGetEigenpair(solver,id-1,X_val(id),PETSC_NULL_OBJECT,&
                &sol_vec,PETSC_NULL_OBJECT,ierr)                                ! get solution EV vector and value (starts at index 0)
            CHCKERR('EPSGetEigenpair failed')
            !ierr = get_ghost_X_vec(X_vec(:,:,id),1)
            !CHCKERR('')
            call EPSComputeRelativeError(solver,id-1,error,ierr)                ! get error (starts at index 0)
            CHCKERR('EPSComputeRelativeError failed')
            call writo('for solution '//trim(i2str(id))//&
                &'/'//trim(i2str(max_n_EV))//':')
            
            !! visualize solution
            !call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,&
                !&'solution '//trim(i2str(id))//' vector with '//&
                !&trim(i2str(n_r_X))//' points',0,&
                !&0,1000,1000,guess_viewer,istat)
            !call VecView(sol_vec,guess_viewer,istat)
            
            ! go back to  physical values from normalization  by multiplying the
            ! Eigenvalues by 1/T0^2 = B_0^2 / (mu_0 rho_0 R_0^2)
            !X_val(id) = X_val(id) * B_0**2 / (mu_0*rho_0*R_0**2)
            
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
    end function store_results
    
    ! stop petsc and slepc
    ! [MPI] Collective call
    integer function stop_slepc(A,B,solver,guess_start_id,prev_n_EV,max_n_EV) &
        &result(ierr)
        use X_vars, only: grp_min_r_X
        
        character(*), parameter :: rout_name = 'stop_slepc'
        
        ! input / output
        Mat, intent(in) :: A, B                                                 ! matrices A and B in EV problem A X = lambda B X
        EPS, intent(in) :: solver                                               ! EV solver
        PetscInt, intent(inout):: guess_start_id                                ! start of index of vector, saved for next iteration
        PetscInt, intent(inout) :: prev_n_EV                                    ! nr. of solutions of previous vector, saved for next iteration
        PetscInt, intent(in) :: max_n_EV                                        ! nr. of EV's saved
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set guess_start_id and prev_n_EV
        guess_start_id = grp_min_r_X
        prev_n_EV = max_n_EV
        
        ! destroy objects
        call EPSDestroy(solver,ierr)
        CHCKERR('Failed to destroy EPS solver')
        call MatDestroy(A,ierr)
        CHCKERR('Failed to destroy matrix A')
        call MatDestroy(B,ierr)
        CHCKERR('Failed to destroy matrix B')
        
        ! stop slepc
        call SlepcFinalize(ierr)
        CHCKERR('Failed to Finalize slepc')
        
        call lvl_ud(-1)
    end function stop_slepc
end module
