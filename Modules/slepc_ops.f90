!------------------------------------------------------------------------------!
!   Operations that use slepc (and petsc) routines
!------------------------------------------------------------------------------!
module slepc_ops
#include <PB3D_macros.h>
#include <finclude/slepcepsdef.h>
    use slepceps
    use num_vars, only: iu, mu_0, dp, max_str_ln
    use output_ops, only: lvl_ud, writo, print_ar_2
    use str_ops, only: r2strt, r2str, i2str

    implicit none
    private
    public solve_EV_system_slepc
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in [ADD REF] and solves them using the slepc suite
    
    integer function solve_EV_system_slepc(m_X,n_r_X,n_X,inv_step) result(ierr)
        use X_vars, only: PV0, PV1, PV2, KV0, KV1, KV2, X_vec, X_val, min_m_X, &
            &max_m_X
        use num_vars, only: MPI_comm_groups, group_n_procs, n_sol_requested
        use file_ops, only: opt_args
        
        character(*), parameter :: rout_name = 'solve_EV_system_slepc'
        
        ! input / output
        PetscInt, intent(in) :: n_r_X                                           ! nr. of normal points for perturbation quantities
        PetscInt, intent(in) :: m_X(:)                                          ! number of poloidal modes m
        PetscInt, intent(in) :: n_X                                             ! toroidal mode number
        PetscReal, intent(in) :: inv_step                                       ! inverse step size
        
        ! local variables
        ! petsc / MPI variables
        PetscMPIInt :: Petsc_rank, n_procs                                      ! rank and nr. of processors
        EPS :: solver                                                           ! EV solver
        PetscInt :: n_it                                                        ! nr. of iterations
        PetscInt :: n_conv                                                      ! nr. of converged solutions
        PetscInt :: n_ev, ncv, mpd                                              ! nr. of requested EV, max. dim of subspace and for projected problem
        EPSType :: EPS_type                                                     ! string with name of type of EPS solver
        PetscReal :: error                                                      ! error of EPS solver
        PetscReal :: tol                                                        ! tolerance of EPS solver
        PetscInt ::  max_it                                                     ! maximum number of iterations
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: option_name                                ! options
        PetscBool :: flg                                                        ! flag to catch options
        
        ! perturbation quantities
        PetscInt :: n_m_X                                                       ! nr. of poloidal modes
        PetscInt :: n_r_X_loc                                                   ! nr. of poloidal modes treated in this processor (=rows)
        Mat :: A                                                                ! matrix A in EV problem A X = lambda B X
        Mat :: B                                                                ! matrix B in EV problem A X = lambda B X
        
        ! solution
        Vec :: vec_par                                                          ! solution EV parallel vector
        Vec :: vec_seq                                                          ! solution EV sequential vector
        PetscReal :: val_r, val_i                                               ! real and imag. part of solution EV value
        VecScatter :: scatter_ctx                                               ! context for scattering of solution vector to group master
        PetscScalar, pointer :: vec_pointer(:)                                  ! pointer to the solution vectors
        
        ! other variables
        PetscInt :: id, jd, kd                                                  ! counters
        PetscInt :: max_n_EV                                                    ! nr. of EV's saved, up to n_sol_requested
        PetscInt :: vec_index                                                   ! vector index to export out of Petsc
        
        ! initialize slepc
        call writo('initialize slepc...')
        call lvl_ud(1)
        ! set PETSC_COMM_WORLD
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
        call MPI_Comm_rank(PETSC_COMM_WORLD,Petsc_rank,ierr)                    ! rank
        CHCKERR('MPI rank failed')
        call MPI_Comm_size(PETSC_COMM_WORLD,n_procs,ierr)                       ! number of processors
        CHCKERR('MPI size failed')
        call writo('slepc started with '//trim(i2str(n_procs))&
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
        n_r_X_loc = n_r_X/n_procs                                               ! number of radial points on this processor
        if (mod(n_r_X,n_procs).gt.0) then                                       ! check if there is a remainder
            if (mod(n_r_X,n_procs).gt.Petsc_rank) n_r_X_loc = n_r_X_loc + 1           ! add a point to the first ranks if there is a remainder
        end if
        
        ! set up the matrix
        call writo('set up matrices...')
        
        call lvl_ud(1)
        
        ! create a  matrix A and B  with the appropriate number  of preallocated
        ! entries
        call MatCreateAIJ(PETSC_COMM_WORLD,n_r_X_loc*n_m_X,n_r_X_loc*n_m_X,&    ! Hermitian matrix as a first approximation
            &n_r_X*n_m_X,n_r_X*n_m_X,2*n_m_X,PETSC_NULL_INTEGER,&
            &2*n_m_X,PETSC_NULL_INTEGER,A,ierr)
        CHCKERR('MatCreateAIJ failed for matrix A')
        call MatCreateAIJ(PETSC_COMM_WORLD,n_r_X_loc*n_m_X,n_r_X_loc*n_m_X,&    ! Hermitian matrix as a first approximation
            &n_r_X*n_m_X,n_r_X*n_m_X,2*n_m_X,PETSC_NULL_INTEGER,&
            &2*n_m_X,PETSC_NULL_INTEGER,B,ierr)
        CHCKERR('MatCreateAIJ failed for matrix B')
        
        ! fill the elements of A and B
        ierr = fill_mat(n_X,n_r_X,n_m_X,inv_step,PV0,PV1,PV2,A)
        CHCKERR('')
        call writo('Matrix A set up')
        ierr = fill_mat(n_X,n_r_X,n_m_X,inv_step,KV0,KV1,KV2,B)
        CHCKERR('')
        call writo('Matrix B set up')
        
        call lvl_ud(-1)
        
        ! visualize the matrices
        !call PetscOptionsSetValue('-draw_pause','-1',ierr)
        write(*,*) 'A ='
        call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !call MatView(A,PETSC_VIEWER_DRAW_WORLD,ierr)
        read(*,*)
        write(*,*) 'B ='
        call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !call MatView(B,PETSC_VIEWER_DRAW_WORLD,ierr)
        read(*,*)
        
        ! solve EV problem
        call writo('solve the EV problem...')
        !call PetscOptionsSetValue('-eps_view','-1',ierr)
        call EPSCreate(PETSC_COMM_WORLD,solver,ierr)
        CHCKERR('EPSCreate failed')
        call EPSSetOperators(solver,A,B,ierr)                                   ! generalized EV problem A X = lambda B X
        CHCKERR('EPSetOperators failed')
        !call EPSSetOperators(solver,A,PETSC_NULL_OBJECT,ierr)                   ! normal EV problem A X = lambda X
        !call EPSSetProblemType(solver,EPS_HEP,ierr)
        call EPSSetFromOptions(solver,ierr)
        CHCKERR('EPSetFromOptions failed')
        call EPSSolve(solver,ierr) 
        CHCKERR('EPS couldn''t find a solution')
        
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
        if (n_sol_requested.gt.n_conv) then
            max_n_EV = n_conv
        else
            max_n_EV = n_sol_requested
        end if
        
        ! allocate X_vec and X_val if group master if first time
        if (Petsc_rank.eq.0) then
            if (allocated(X_vec)) deallocate(X_vec)
            allocate(X_vec(min_m_X:max_m_X,1:n_r_X,1:max_n_EV))
            if (allocated(X_val)) deallocate(X_val)
            allocate(X_val(1:max_n_EV))
        end if
        
        ! print info
        call writo('storing results with '//trim(i2str(max_n_EV))//' highest &
            &Eigenvalues...')
        
        call lvl_ud(1)
        
        call MatGetVecs(A,vec_par,PETSC_NULL_OBJECT,ierr)                       ! get compatible parallel vector to matrix A
        CHCKERR('MatGetVecs failed')
        
        ! store them
        do id = 1,max_n_EV
            ! get EV solution
            call EPSGetEigenpair(solver,id-1,val_r,val_i,vec_par,&
                &PETSC_NULL_OBJECT,ierr)                                        ! get solution EV vector and value (starts at index 0)
            CHCKERR('EPSGetEigenpair failed')
            call EPSComputeRelativeError(solver,id-1,error,ierr)                ! get error (starts at index 0)
            CHCKERR('EPSComputeRelativeError failed')
            call writo('for solution '//trim(i2str(id))//&
                &'/'//trim(i2str(n_conv))//':')
            
            call lvl_ud(1)
            
            ! print output
            call writo('eigenvalue: '//trim(r2str(val_r))//' + '//&
                &trim(r2str(val_i))//' i')
            call writo('with rel. error: '//trim(r2str(error)))
            !call writo('Eigenvector = ')
            !call VecView(vec_par,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !CHCKERR('Cannot view vector')
            
            call lvl_ud(-1)
            
            ! get solution vector to group master
            call VecScatterCreateToZero(vec_par,scatter_ctx,vec_seq,ierr)       ! create scatter to group master
            err_msg = 'Failed to create scatter to all context'
            CHCKERR(err_msg)
            call VecScatterBegin(scatter_ctx,vec_par,vec_seq,INSERT_VALUES,&
                &SCATTER_FORWARD,ierr)                                          ! begin scatter to group master
            err_msg = 'Failed to begin scatter to group master'
            CHCKERR(err_msg)
            call VecScatterEnd(scatter_ctx,vec_par,vec_seq,INSERT_VALUES,&
                &SCATTER_FORWARD,ierr)                                          ! end scatter to group master
            err_msg = 'Failed to end scatter to group master'
            CHCKERR(err_msg)
            call VecScatterDestroy(scatter_ctx,ierr)                            ! destroy scatter to group master
            err_msg = 'Failed to destroy scatter to all context'
            CHCKERR(err_msg)
            
            ! export solution vector to X_vec
            call VecGetArrayF90(vec_seq,vec_pointer,ierr)
            err_msg = 'Failed to get pointer to solution vector'
            CHCKERR(err_msg)
            if (Petsc_rank.eq.0) then
                do kd = 1,n_r_X
                    do jd = min_m_X,max_m_X
                        vec_index = (kd-1) * n_m_X + (jd - min_m_X) + 1
                        X_vec(jd,kd,id) = vec_pointer(vec_index)
                    end do
                end do
                X_val(id) = val_r + val_i * PETSC_i
            end if
        end do
        
        call VecDestroy(vec_par,ierr)                                           ! destroy compatible vector
        CHCKERR('Failed to destroy vec_par')
        
        call lvl_ud(-1)
        
        ! destroy and finalize
        call writo('finalize petsc...')
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
    ! !!!! THE BOUNDARY CONDITIONS ARE STILL MISSING, SO IT IS TREATED AS HERMITIAN !!!
    integer function fill_mat(n_X,n_r_X,n_m_X,inv_step,V0,V1,V2,mat) result(ierr)
        use metric_ops, only: jac_F_FD
        use X_vars, only: min_r, max_r
        use VMEC_vars, only: n_r
        use eq_vars, only: n_par, theta
        use utilities, only: calc_interp, calc_int
        
        character(*), parameter :: rout_name = 'fill_mat'
        
        ! input / output
        PetscInt, intent(in) :: n_X                                             ! toroidal mode number
        PetscInt, intent(in) :: n_r_X                                           ! nr. of normal points for perturbation quantities
        PetscInt, intent(in) :: n_m_X                                           ! nr. of poloidal modes
        PetscScalar, intent(in) :: V0(:,:,:,:), V1(:,:,:,:), V2(:,:,:,:)        ! either PVi or KVi
        Mat, intent(inout) :: mat                                               ! either A or B
        PetscReal, intent(in) :: inv_step                                       ! inverse step size
        
        ! local variables
        PetscScalar, allocatable :: loc_block(:,:)                              ! block matrix for 1 normal point
        PetscScalar, allocatable :: intgrnd(:,:,:,:)                            ! constituents of loc_bloc, not interpolated
        PetscScalar, allocatable :: intgrnd_int(:,:,:)                          ! constituents of loc_bloc, interpolated
        PetscScalar, allocatable :: exp_theta(:,:,:,:)                          ! exp(i (k-m) theta)
        PetscReal, allocatable :: theta_int(:,:,:)                              ! theta at interpolated position
        PetscScalar, allocatable :: integral(:)                                 ! the integral of the integrand
        PetscInt, allocatable :: loc_k(:), loc_m(:)                             ! the locations at which to add the blocks to the matrices
        PetscReal :: norm_pt, norm_pt_5, norm_pt_1                              ! current normal point (0...1) and pt+1/2, pt+1
        PetscInt :: id, jd, m, k                                                ! counters
        PetscInt :: n_start, n_end                                              ! start row and end row
        PetscInt :: r_X_start, r_X_end                                          ! start block and end block (= n_start/n_m_X and equiv.)
        
        ! initialize ierr
        ierr = 0
        
        allocate(intgrnd(n_par,n_r,n_m_X,n_m_X))
        allocate(intgrnd_int(n_par,n_m_X,n_m_X))
        allocate(loc_block(n_m_X,n_m_X))
        allocate(theta_int(n_par,1,1))                                          ! dimension 3 to be able to use the routine from utilities
        allocate(integral(n_par))
        allocate(exp_theta(n_par,n_r,n_m_X,n_m_X))
        do m = 1,n_m_X
            do k = 1,n_m_X
                exp_theta(:,:,k,m) = exp(PETSC_i*(k-m)*theta)
            end do
        end do
        allocate(loc_k(n_m_X))
        allocate(loc_m(n_m_X))
        
        ! get ownership of the rows
        call MatGetOwnershipRange(mat,n_start,n_end,ierr)                       ! starting and ending row n_start and n_end
        CHCKERR('Couldn''t get ownership range of matrix')
        r_X_start = n_start/n_m_X
        r_X_end = n_end/n_m_X
        
        ! iterate over all rows of this rank
        do id = r_X_start, r_X_end-1
            ! calculate the  normal points at which  a row is to  be added. This
            ! should go from min_r to max_r,  making use of n_r_X points (so the
            ! global block  row number id goes  from 0 to n_r_X-1).  This yields
            ! the equation:
            !   norm_pt - min_r         id
            !   ---------------  =   --------- ,
            !    max_r - min_r       n_r_X - 1
            ! which gives a formula for norm_pt
            norm_pt   = min_r + (max_r-min_r) * (id-0.0)/(n_r_X-1.0)            ! id
            norm_pt_5 = min_r + (max_r-min_r) * (id+0.5)/(n_r_X-1.0)            ! id + 0.5
            norm_pt_1 = min_r + (max_r-min_r) * (id+1.0)/(n_r_X-1.0)            ! id + 1
            
            ! ---------------!
            ! BLOCKS (id,id) !
            ! ---------------!
            loc_block = 0
            ! part ~ V0(i)
            do m = 1,n_m_X
                do k = 1,n_m_X
                    intgrnd(:,:,k,m) = exp_theta(:,:,k,m) * &
                        &jac_F_FD(:,:,0,0,0) * V0(:,:,k,m)                      ! integrand, not integrated
                    ierr = calc_interp(intgrnd,intgrnd_int,norm_pt)             ! integrand, integrated
                    CHCKERR('')
                    ierr = calc_interp(reshape(theta,[n_par,n_r,1,1]),&
                        &theta_int,norm_pt)                                     ! theta array at interpolated position
                    CHCKERR('')
                    ierr = calc_int(intgrnd_int(:,k,m),theta_int(:,1,1),&
                        &integral)                                              ! integral of entire function
                    CHCKERR('')
                    loc_block(k,m) = loc_block(k,m) + integral(n_par)           ! add total integral (at pos n_par)
                end do
            end do
            ! part ~ V2(i)
            do m = 1,n_m_X
                do k = 1,n_m_X
                    intgrnd(:,:,k,m) = 2 * exp_theta(:,:,k,m) * &
                        &(inv_step/n_X)**2 * jac_F_FD(:,:,0,0,0) * &
                        &V2(:,:,k,m)                                            ! integrand, not integrated
                    ierr = calc_interp(intgrnd,intgrnd_int,norm_pt)             ! integrand, integrated
                    CHCKERR('')
                    ierr = calc_interp(reshape(theta,[n_par,n_r,1,1]),&
                        &theta_int,norm_pt)                                     ! theta array at interpolated position
                    CHCKERR('')
                    ierr = calc_int(intgrnd_int(:,k,m),theta_int(:,1,1),&
                        &integral)                                              ! integral of entire function
                    CHCKERR('')
                    loc_block(k,m) = loc_block(k,m) + integral(n_par)           ! add total integral (at pos n_par)
                end do
            end do
            
            ! add block to the matrix A
            loc_k = [(jd, jd = 0,n_m_X-1)] + id*n_m_X
            loc_m = loc_k
            call MatSetValues(mat,n_m_X,loc_k,n_m_X,loc_m,loc_block,&
                &INSERT_VALUES,ierr)
            CHCKERR('Couldn''t add values to matrix')
            
            ! -----------------!
            ! BLOCKS (id,id+1) !
            ! -----------------!
            if (id.ne.n_r_X-1) then                                             ! don't do the right block for last normal point
                loc_block = 0
                ! part ~ V1*(i+1)
                do m = 1,n_m_X
                    do k = 1,n_m_X
                        intgrnd(:,:,k,m) = exp_theta(:,:,k,m) * &
                            &PETSC_i*inv_step/(2*n_X) * &
                            &jac_F_FD(:,:,0,0,0) * conjg(V1(:,:,m,k))           ! integrand, not integrated
                        ierr = calc_interp(intgrnd,intgrnd_int,norm_pt_1)       ! integrand, integrated
                        CHCKERR('')
                        ierr = calc_interp(reshape(theta,[n_par,n_r,1,1]),&
                            &theta_int,norm_pt_1)                               ! theta array at interpolated position
                        CHCKERR('')
                        ierr = calc_int(intgrnd_int(:,k,m),theta_int(:,1,1),&
                            &integral)                                          ! integral of entire function
                        CHCKERR('')
                        loc_block(k,m) = loc_block(k,m) + integral(n_par)       ! add total integral (at pos n_par)
                    end do
                end do
                ! part ~ V1(i)
                do m = 1,n_m_X
                    do k = 1,n_m_X
                        intgrnd(:,:,k,m) = exp_theta(:,:,k,m) * &
                            &PETSC_i*inv_step/(2*n_X) * &
                            &jac_F_FD(:,:,0,0,0) * V1(:,:,k,m)                  ! integrand, not integrated
                        ierr = calc_interp(intgrnd,intgrnd_int,norm_pt)         ! integrand, integrated
                        CHCKERR('')
                        ierr = calc_interp(reshape(theta,[n_par,n_r,1,1]),&
                            &theta_int,norm_pt)                                 ! theta array at interpolated position
                        CHCKERR('')
                        ierr = calc_int(intgrnd_int(:,k,m),theta_int(:,1,1),&
                            &integral)                                          ! integral of entire function
                        CHCKERR('')
                        loc_block(k,m) = loc_block(k,m) + integral(n_par)       ! add total integral (at pos n_par)
                    end do
                end do
                ! part ~ V2(i+1/2)
                do m = 1,n_m_X
                    do k = 1,n_m_X
                        intgrnd(:,:,k,m) = - exp_theta(:,:,k,m) * &
                            &(inv_step/n_X)**2 * jac_F_FD(:,:,0,0,0) * &
                            &V2(:,:,k,m)                                        ! integrand, not integrated
                        ierr = calc_interp(intgrnd,intgrnd_int,norm_pt_5)       ! integrand, integrated
                        CHCKERR('')
                        ierr = calc_interp(reshape(theta,[n_par,n_r,1,1]),&
                            &theta_int,norm_pt_5)                               ! theta array at interpolated position
                        CHCKERR('')
                        ierr = calc_int(intgrnd_int(:,k,m),theta_int(:,1,1),&
                            &integral)                                          ! integral of entire function
                        CHCKERR('')
                        loc_block(k,m) = loc_block(k,m) + integral(n_par)       ! add total integral (at pos n_par)
                    end do
                end do
                
                ! add block to the matrix A
                loc_k = [(jd, jd = 0,n_m_X-1)] + id*n_m_X
                loc_m = loc_k+n_m_X
                call MatSetValues(mat,n_m_X,loc_k,n_m_X,loc_m,loc_block,&
                    &INSERT_VALUES,ierr)
                CHCKERR('Couldn''t add values to matrix')
            end if
        end do
        
        ! assemble the matrix and view it
        call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t begin assembly of matrix')
        call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t end assembly of matrix')
        
        call MatSetOption(mat,MAT_HERMITIAN,PETSC_TRUE,ierr)                    ! Hermitian to a first approximation
        CHCKERR('Coulnd''t set option Hermitian')
    end function fill_mat
end module
