!------------------------------------------------------------------------------!
!   Operations that use SLEPC (and PETSC) routines
!------------------------------------------------------------------------------!
module SLEPC_ops
#include <PB3D_macros.h>
#include <finclude/slepcepsdef.h>
!#include <finclude/petscsys.h>
    use str_ops
    use output_ops
    use messages
    use slepceps
    use num_vars, only: iu, dp, max_str_ln
    use grid_vars, only: grid_type
    use X_vars, only: X_type

    implicit none
    private
    public solve_EV_system_SLEPC
#if ldebug
    public debug_setup_matrices, debug_fill_mat
#endif
    
    ! global variables
#if ldebug
    logical :: debug_setup_matrices = .false.                                   ! plot debug information for setup_matrices
    logical :: debug_fill_mat = .false.                                         ! plot debug information for fill_mat
#endif
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in  [ADD REF] and  solves them using  the SLEPC suite.  The most
    ! unstable solutions  are obtained, where the  variable "max_n_EV" indicates
    ! how many.
    integer function solve_EV_system_SLEPC(grid_eq,grid_X,X,use_guess,&
        &max_n_EV) result(ierr)
        character(*), parameter :: rout_name = 'solve_EV_system_SLEPC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        PetscBool, intent(in) :: use_guess                                      ! whether to use a guess or not
        PetscInt, intent(inout) :: max_n_EV                                     ! how many solutions saved
        
        ! local variables
        Mat :: A                                                                ! matrix A in EV problem A X = lambda B X
        Mat :: B                                                                ! matrix B in EV problem A X = lambda B X
        EPS :: solver                                                           ! EV solver
        PetscInt, save :: guess_start_id = -10                                  ! start of index of previous vector, saved for next iteration
        PetscInt, save :: prev_n_EV                                             ! nr. of solutions of previous vector
        
        ! initialize ierr
        ierr = 0
        
        ! start SLEPC
        ierr = start_SLEPC(grid_X%n(3))
        CHCKERR('')
        
        ! set up the matrix
        call writo('set up matrices...')
        
        ierr = setup_matrices(grid_eq,grid_X,X,A,B)
        CHCKERR('')
        
        ! set up solver
        call writo('set up EV solver...')
        
        ierr = setup_solver(grid_X,X,A,B,solver)
        CHCKERR('')
        
        ! set up guess
        call writo('set up guess...')
        
        if (use_guess) call setup_guess(X,A,solver,guess_start_id,prev_n_EV)
        
        ! get solution
        call writo('get solution...')
        
        ierr = get_solution(solver)
        CHCKERR('')
        
        ! summarize solution
        call writo('summarize solution...')
        
        ierr = summarize_solution(solver,max_n_EV)
        CHCKERR('')
        
        ! store results
        call writo('storing results for '//trim(i2str(max_n_EV))//' highest &
            &Eigenvalues...')
        
        ierr = store_results(grid_X,X,solver,max_n_EV)
        CHCKERR('')
        
        ! finalize
        call writo('finalize SLEPC...')
        
        ierr = stop_SLEPC(grid_X,A,B,solver,guess_start_id,prev_n_EV,max_n_EV)
        CHCKERR('')
    end function solve_EV_system_SLEPC
        
    !  This  subroutine starts  PETSC  and  SLEPC  with  the correct  number  of
    ! processors
    integer function start_SLEPC(n_r_X) result(ierr)
        use num_vars, only: MPI_Comm_groups, grp_n_procs
        use files, only: opt_args
        
        character(*), parameter :: rout_name = 'start_SLEPC'
        
        ! input / output
        integer, intent(in) :: n_r_X                                            ! nr. of normal points in perturbation grid
        
        ! local variables
        PetscBool :: flg                                                        ! flag to catch options
        PetscInt :: id                                                          ! counters
        character(len=max_str_ln) :: option_name                                ! options
        
        ! initialize ierr
        ierr = 0
        
        ! user message
        call writo('initialize SLEPC...')
        call lvl_ud(1)
        
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
        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)                         ! initialize SLEPC
        CHCKERR('SLEPC failed to initialize')
        
        ! output
        call writo('slepc started with '//trim(i2str(grp_n_procs))&
            &//' processors')
        
        call lvl_ud(-1)
        
        ! checking for complex numbers
        call writo('run tests...')
#if defined(PETSC_USE_COMPLEX)
        ! OK
#else
        err_msg = 'PETSC and SLEPC have to be configured and compiled &
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
    end function start_SLEPC
    
    ! Sets  up the  matrices  A and  B in  the  EV problem  A  X =  lambda B  X.
    ! Optionally the  geodesic index of  the perturbation variables at  which to
    ! perform the calculations can be changed  from its default value of 1 using
    ! the variable i_geo.
    ! Note: For normal usage, i_geo should be 1, or not present.
    integer function setup_matrices(grid_eq,grid_X,X,A,B,i_geo) result(ierr)
        use num_vars, only: grp_n_procs, grp_rank
        use grid_ops, only: get_norm_interp_data
        
        character(*), parameter :: rout_name = 'setup_matrices'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        integer, intent(in), optional :: i_geo                                  ! at which geodesic index to perform the calculations
        
        ! local variables
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        PetscReal, allocatable :: grp_r_eq(:)                                   ! unrounded index in tables V_int
        integer :: i_geo_loc                                                    ! local copy of i_geo
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        
#if ldebug
        Mat :: A_t                                                              ! Hermitian transpose of A
        Mat :: B_t                                                              ! Hermitian transpose of B
        PetscScalar :: one = 1.0                                                ! one
#endif
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set up local i_geo
        i_geo_loc = 1
        if (present(i_geo)) i_geo_loc = i_geo
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! set up arrays grp_r_eq that determine how to fill the matrices
        ierr = get_norm_interp_data(grid_eq,grid_X,grp_r_eq)
        CHCKERR('')
        
        ! create a  matrix A and B  with the appropriate number  of preallocated
        ! entries (excluding ghost regions)
        ! the numbers of non-zeros in the diagonal and off-diagonal parts
        allocate(d_nz(grp_n_r_X_loc*X%n_mod)); d_nz = 0
        allocate(o_nz(grp_n_r_X_loc*X%n_mod)); o_nz = 0
        if (grp_n_r_X_loc.eq.1) then                                            ! 1 normal point in the process
            d_nz = X%n_mod
            o_nz = 2*X%n_mod
        else                                                                    ! more than 1 normal point in the process
            d_nz(1:X%n_mod) = 2*X%n_mod
            d_nz(X%n_mod+1:(grp_n_r_X_loc-1)*X%n_mod) = 3*X%n_mod
            d_nz((grp_n_r_X_loc-1)*X%n_mod+1:grp_n_r_X_loc*X%n_mod) = &
                &2*X%n_mod
            if(grp_rank.ne.0) then
                o_nz(1:X%n_mod) = X%n_mod
            end if
            if (grp_rank.ne.grp_n_procs-1) then
                o_nz((grp_n_r_X_loc-1)*X%n_mod+1:grp_n_r_X_loc*X%n_mod) = &
                    &X%n_mod
            end if
        end if
        
        ! create matrix A
        call MatCreateAIJ(PETSC_COMM_WORLD,grp_n_r_X_loc*X%n_mod,&
            &grp_n_r_X_loc*X%n_mod,grid_X%n(3)*X%n_mod,grid_X%n(3)*X%n_mod,&
            &PETSC_NULL_INTEGER,d_nz,PETSC_NULL_INTEGER,o_nz,A,ierr)
        CHCKERR('MatCreateAIJ failed for matrix A')
        
        ! deallocate d_nz and o_nz
        deallocate(d_nz,o_nz)
        
        ! fill the matrix A
        ierr = fill_mat(X%PV_int_0(:,i_geo_loc,:),X%PV_int_1(:,i_geo_loc,:),&
            &X%PV_int_2(:,i_geo_loc,:),A)
        CHCKERR('')
        call writo('matrix A set up')
        
        ! duplicate A into B
        ! (the  advantage of  letting communication  and calculation  overlap by
        ! calculating matrix B while assembling A is offset by the extra cost of
        ! not being able to use MatDuplicate. This is also easier)
        call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)                   ! B has same structure as A
        CHCKERR('failed to duplicate A into B')
        
        ! fill the matrix B
        ierr = fill_mat(X%KV_int_0(:,i_geo_loc,:),X%KV_int_1(:,i_geo_loc,:),&
            &X%KV_int_2(:,i_geo_loc,:),B)
        CHCKERR('')
        call writo('matrix B set up')
        
#if ldebug
        if (debug_setup_matrices) then
            ! test if A and B hermitian
            call MatHermitianTranspose(A,MAT_INITIAL_MATRIX,A_t,ierr)
            CHCKERR('Hermitian transpose of A failed')
            call MatAXPY(A_t,-one,A,SAME_NONZERO_PATTERN,ierr)
            CHCKERR('A-A_t failed')
            call MatHermitianTranspose(B,MAT_INITIAL_MATRIX,B_t,ierr)
            CHCKERR('Hermitian transpose of B failed')
            call MatAXPY(B_t,-one,B,SAME_NONZERO_PATTERN,ierr)
            CHCKERR('B-B_t failed')
            
            ! visualize the matrices
            call PetscOptionsSetValue('-draw_pause','-1',ierr)
            write(*,*) 'A ='
            call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !call MatView(A,PETSC_VIEWER_DRAW_WORLD,ierr)
            write(*,*) 'B ='
            call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !call MatView(B,PETSC_VIEWER_DRAW_WORLD,ierr)
            write(*,*) 'A_t ='
            call MatView(A_t,PETSC_VIEWER_STDOUT_WORLD,ierr)
            write(*,*) 'B_t ='
            call MatView(B_t,PETSC_VIEWER_STDOUT_WORLD,ierr)
            
            ! destroy matrices
            call MatDestroy(A_t,ierr)
            CHCKERR('Failed to destroy matrix A_t')
            call MatDestroy(B_t,ierr)
            CHCKERR('Failed to destroy matrix B_t')
        end if
#endif
        
        ! set boundary conditions
        ierr = set_mat_BC(A,B)
        CHCKERR('')
        
#if ldebug
        if (debug_setup_matrices) then
            ! visualize the matrices
            call PetscOptionsSetValue('-draw_pause','-1',ierr)
            write(*,*) 'A ='
            call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !call MatView(A,PETSC_VIEWER_DRAW_WORLD,ierr)
            write(*,*) 'B ='
            call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !call MatView(B,PETSC_VIEWER_DRAW_WORLD,ierr)
        end if
#endif
        
        ! deallocate quantities
        deallocate(grp_r_eq)
        
        call lvl_ud(-1)
    contains
        ! fills a  matrix according to  [ADD REF]. Is  used for both  the matrix
        ! A  and  B,  corresponding  to  the plasma  potential  energy  and  the
        ! (perpendicular) kinetic energy
        ! The procedure is as follows:
        !   1. Tables  are set  up for  (k,m) pairs in the equilibrium grid.
        !   2. At the first normal point  in the perturbation grid, belonging to
        !   the current process, the tables are interpolated as a start.
        !   3. At every normal point in  the perturbation grid, belonging to the
        !   current process,  the tables are  interpolated at the next  point in
        !   the  corresponding  perturbation  grid. This  information,  as  well
        !   as  the interpolated  value  of  the previous  point  allow for  the
        !   calculation of every quantity.
        ! Note: the factors i/n or i/m are already included.
        ! !!!! THE REAL BOUNDARY CONDITIONS ARE STILL MISSING !!!!!!
        integer function fill_mat(V_0,V_1,V_2,mat) result(ierr)
            use num_vars, only: use_pol_flux_F
            use X_vars, only: min_r_X, max_r_X
            use utilities, only: c, con
            use eq_vars, only: max_flux_p_F, max_flux_t_F
            
            character(*), parameter :: rout_name = 'fill_mat'
            
            ! input / output
            PetscScalar, intent(in) :: V_0(:,:)                                 ! either PV_int_0 or KV_int_0 in equilibrium normal grid
            PetscScalar, intent(in) :: V_1(:,:)                                 ! either PV_int_1 or KV_int_1 in equilibrium normal grid
            PetscScalar, intent(in) :: V_2(:,:)                                 ! either PV_int_2 or KV_int_2 in equilibrium normal grid
            Mat, intent(inout) :: mat                                           ! either A or B
            
            ! local variables (not to be used in child routines)
            character(len=max_str_ln) :: err_msg                                ! error message
            PetscInt :: nn_mod_1, nn_mod_2                                      ! number of indices for a quantity that is symmetric or not
            PetscScalar, allocatable :: loc_block(:,:)                          ! (n_mod x n_mod) block matrix for 1 normal point
            PetscScalar, allocatable :: V_int_0(:)                              ! interpolated V_0
            PetscScalar, allocatable :: V_int_1(:)                              ! interpolated V_1
            PetscScalar, allocatable :: V_int_2(:)                              ! interpolated V_2
            PetscScalar, allocatable :: V_int_0_next(:)                         ! interpolated V_0 at next normal point
            PetscScalar, allocatable :: V_int_1_next(:)                         ! interpolated V_1 at next normal point
            PetscScalar, allocatable :: V_int_2_next(:)                         ! interpolated V_2 at next normal point
            PetscInt, allocatable :: loc_k(:), loc_m(:)                         ! the locations at which to add the blocks to the matrices
            PetscReal :: step_size                                              ! step size in flux coordinates
            PetscInt :: id, jd                                                  ! counters
            PetscInt :: k, m                                                    ! counters
            integer :: i_max_X_loc                                              ! local copy of i_max of X grid
            PetscReal :: max_flux_F                                             ! maximum flux in Flux coordinates
            
            ! for tests
            PetscInt :: r_X_start, r_X_end                                      ! start block and end block
            
            ! initialize ierr
            ierr = 0
            
            ! set up max_flux_F
            if (use_pol_flux_F) then
                max_flux_F = max_flux_p_F
            else
                max_flux_F = max_flux_t_F
            end if
            
            ! set nn_mod_1 and nn_mod_2
            nn_mod_1 = X%n_mod**2
            nn_mod_2 = X%n_mod*(X%n_mod+1)/2
            
            ! set up local grp_n_r of X grid
            if (grp_rank.lt.grp_n_procs-1) then
                i_max_X_loc = grid_X%i_max - 1                                  ! i_max of X grid includes ghost region
            else
                i_max_X_loc = grid_X%i_max
            end if
            
            ! allocate variables
            allocate(V_int_0(nn_mod_2))
            allocate(V_int_1(nn_mod_1))
            allocate(V_int_2(nn_mod_2))
            allocate(V_int_0_next(nn_mod_2))
            allocate(V_int_1_next(nn_mod_1))
            allocate(V_int_2_next(nn_mod_2))
            allocate(loc_block(X%n_mod,X%n_mod))
            allocate(loc_k(X%n_mod))
            allocate(loc_m(X%n_mod))
            
            ! test whether the matrix range coincides with i_min and i_max
            call MatGetOwnershipRange(mat,r_X_start,r_X_end,ierr)               ! starting and ending row r_X_start and r_X_end
            err_msg = 'Couldn''t get ownership range of matrix'
            CHCKERR(err_msg)
            r_X_start = r_X_start/X%n_mod                                       ! count per block
            r_X_end = r_X_end/X%n_mod                                           ! count per block
            if (grid_X%i_min.ne.r_X_start+1) then
                ierr = 1
                err_msg = 'start of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            if (i_max_X_loc.ne.r_X_end) then
                ierr = 1
                err_msg = 'end of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            
            ! set up step_size
            step_size = (max_r_X-min_r_X)/(grid_X%n(3)-1.0) * max_flux_F        ! equidistant perturbation grid in perturb. coords.
            
            ! get the interpolated terms in V_int_i_next
            call interp_V(V_0,grp_r_eq(1),V_int_0_next)
            call interp_V(V_1,grp_r_eq(1),V_int_1_next)
            call interp_V(V_2,grp_r_eq(1),V_int_2_next)
            
            ! iterate over all rows of this rank
            do id = grid_X%i_min-1, i_max_X_loc-1                               ! (indices start with 0 here)
                ! fill V_int_i with V_int_i_next of previous step
                V_int_0 = V_int_0_next
                V_int_1 = V_int_1_next
                V_int_2 = V_int_2_next
                
                ! ---------------!
                ! BLOCKS (id,id) !
                ! ---------------!
                loc_block = 0
                do m = 1,X%n_mod
                    do k = 1,X%n_mod
                        ! part ~ V0(i)
                        loc_block(k,m) = loc_block(k,m) + &
                            &con(V_int_0(c([k,m],.true.,X%n_mod)),[k,m],.true.)
                        ! part ~ V2(i)
                        loc_block(k,m) = loc_block(k,m) - &
                            &con(V_int_2(c([k,m],.true.,X%n_mod)),&
                            &[k,m],.true.) * 2.0_dp/(step_size)**2
                    end do
                end do
                
                ! add block to the matrix A
                loc_k = [(jd, jd = 0,X%n_mod-1)] + id*X%n_mod
                loc_m = loc_k
#if ldebug
                if (debug_fill_mat) then
                    call writo('at (k,m) = ')
                    write(*,*) loc_k
                    write(*,*) loc_m
                    call writo('following local block is added: ')
                    write(*,*) 'Re ='
                    call print_ar_2(realpart(loc_block))
                    write(*,*) 'Im ='
                    call print_ar_2(imagpart(loc_block))
                endif
#endif
                call MatSetValues(mat,X%n_mod,loc_k,X%n_mod,loc_m,loc_block,&
                    &INSERT_VALUES,ierr)
                CHCKERR('Couldn''t add values to matrix')
                
                ! -----------------------------------------!
                ! BLOCKS (id,id+1) and hermitian conjugate !
                ! -----------------------------------------!
                if (id.ne.grid_X%n(3)-1) then                                   ! don't do the right block for last normal point
                    ! get the interpolated terms in V_int_i_next
                    call interp_V(V_0,grp_r_eq(id-grid_X%i_min+3),V_int_0_next)
                    call interp_V(V_1,grp_r_eq(id-grid_X%i_min+3),V_int_1_next)
                    call interp_V(V_2,grp_r_eq(id-grid_X%i_min+3),V_int_2_next)
                    
                    loc_block = 0
                    do m = 1,X%n_mod
                        do k = 1,X%n_mod
                            ! part ~ V1*(i+1)
                            loc_block(k,m) = loc_block(k,m) + &
                                &1.0_dp/(2*step_size) * &
                                &conjg(V_int_1_next(c([m,k],.false.,X%n_mod)))  ! transpose, assymetric matrices don't need con()
                            ! part ~ V1(i)
                            loc_block(k,m) = loc_block(k,m) + &
                                &1.0_dp/(2*step_size) * &
                                &V_int_1(c([k,m],.false.,X%n_mod))              ! assymetric matrices don't need con()
                            ! part ~ V2(i+1/2)
                            loc_block(k,m) = loc_block(k,m) + &
                                &con(V_int_2(c([k,m],.true.,X%n_mod)) + &
                                &V_int_2_next(c([k,m],.true.,X%n_mod)),&
                                &[k,m],.true.)/2 * 1.0/(step_size)**2
                        end do
                    end do
                    
                    ! add block to the matrix
                    loc_k = [(jd, jd = 0,X%n_mod-1)] + id*X%n_mod
                    loc_m = loc_k+X%n_mod
                    call MatSetValues(mat,X%n_mod,loc_k,X%n_mod,loc_m,&
                        &loc_block,INSERT_VALUES,ierr)
                    err_msg = 'Couldn''t add values to matrix'
                    CHCKERR(err_msg)
                    
                    ! add Hermitian conjugate
                    loc_block = conjg(transpose(loc_block))
                    call MatSetValues(mat,X%n_mod,loc_m,X%n_mod,loc_k,&
                        &loc_block,INSERT_VALUES,ierr)
                    err_msg = 'Couldn''t add values to matrix'
                    CHCKERR(err_msg)
                end if
            end do
            
            ! assemble the matrix
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t begin assembly of matrix')
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t end assembly of matrix')
            
            !call MatSetOption(mat,MAT_HERMITIAN,PETSC_TRUE,ierr)                ! Hermitian to a first approximation
            !CHCKERR('Coulnd''t set option Hermitian')
            
            ! deallocate variables
            deallocate(V_int_0,V_int_1,V_int_2)
            deallocate(V_int_0_next,V_int_1_next,V_int_2_next)
            deallocate(loc_block)
            deallocate(loc_k,loc_m)
        end function fill_mat
        
        ! Interpolates an array V at position x.
        ! (also  correct  if i_hi = i_lo)
        subroutine interp_V(V,x,V_int)
            ! input / output
            PetscScalar, intent(in) :: V(:,:)                                   ! V
            PetscReal, intent(in) :: x                                          ! coordinate at which to interpolate
            PetscScalar, intent(inout) :: V_int(:)                              ! interpolated V
            
            ! local variables
            PetscInt :: i_lo, i_hi                                              ! upper and lower index
            
            ! set up i_lo and i_hi
            i_lo = floor(x)
            i_hi = ceiling(x)
            
            V_int = V(:,i_lo) + (V(:,i_hi)-V(:,i_lo))*(x-i_lo)                  ! because i_hi - i_lo = 1
        end subroutine interp_V
        
        ! Sets the  boundary conditions by overwriting  the Hermitian components
        ! written in fill_mat.  An artificial Eigenvalue EV_BC  is introduced by
        ! setting the diagonal components  of A to EV_BC and of B  to 1, and the
        ! off-diagonal elements to zero.
        ! The boundary condition at the plasma surface IS STILL AN OPEN QUESTION
        integer function set_mat_BC(A,B) result(ierr)
            use num_vars, only: EV_BC
            
            character(*), parameter :: rout_name = 'set_mat_BC'
            
            ! input / output
            Mat, intent(inout) :: A, B                                          ! Matrices A and B from A X = lambda B X
            
            ! local variables
            PetscScalar, allocatable :: loc_block_one(:,:)                      ! (n_mod,n_mod) block matrix for 1 normal point, one
            PetscScalar, allocatable :: loc_block_zero(:,:)                     ! (n_mod,n_mod) block matrix for 1 normal point, zero
            PetscInt, allocatable :: loc_k(:), loc_m(:)                         ! the locations at which to add the blocks to the matrices
            PetscInt :: id                                                      ! counter
            
            ! initialize ierr
            ierr = 0
            
            ! initialize local blocks
            allocate(loc_block_one(X%n_mod,X%n_mod))
            loc_block_one = 0.0_dp
            do id = 1,X%n_mod
                loc_block_one(id,id) = 1.0_dp
            end do
            allocate(loc_block_zero(X%n_mod,X%n_mod))
            loc_block_zero = 0.0_dp
            do id = 1,X%n_mod
                loc_block_zero(id,id) = 0._dp
            end do
            
            if (grp_rank.eq.0) then
                ! -------------!
                ! BLOCKS (0,0) !
                ! -------------!
                ! set indices where to insert the block
                loc_k = [(id, id = 0,X%n_mod-1)]
                loc_m = loc_k
                
                ! set the values of A (EV_BC)
                call MatSetValues(A,X%n_mod,loc_k,X%n_mod,loc_m,&
                    &EV_BC*loc_block_one,INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! set the values of B (one)
                call MatSetValues(B,X%n_mod,loc_k,X%n_mod,loc_m,loc_block_one,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! -------------!
                ! BLOCKS (0,1) !
                ! -------------!
                ! set indices where to insert the block
                loc_k = [(id, id = 0,X%n_mod-1)]
                loc_m = loc_k + X%n_mod
                
                ! set the values of A (zero)
                call MatSetValues(A,X%n_mod,loc_k,X%n_mod,loc_m,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! set the values of B (zero)
                call MatSetValues(B,X%n_mod,loc_k,X%n_mod,loc_m,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! -------------!
                ! BLOCKS (1,0) !
                ! -------------!
                ! set the values of A (zero)
                call MatSetValues(A,X%n_mod,loc_m,X%n_mod,loc_k,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! set the values of B (zero)
                call MatSetValues(B,X%n_mod,loc_m,X%n_mod,loc_k,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
            end if
            
            if (grp_rank.eq.grp_n_procs-1) then
                ! -------------------------!
                ! BLOCKS (n-1,n-1) !
                ! -------------------------!
                ! set indices where to insert the block
                loc_k = [(id, id = 0,X%n_mod-1)] + (grid_X%n(3)-1)*X%n_mod
                loc_m = loc_k
                
                ! set the values of A (EV_BC)
                call MatSetValues(A,X%n_mod,loc_k,X%n_mod,loc_m,&
                    &EV_BC*loc_block_one,INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! set the values of B (one)
                call MatSetValues(B,X%n_mod,loc_k,X%n_mod,loc_m,loc_block_one,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! -------------------------!
                ! BLOCKS (n-1,n-2) !
                ! -------------------------!
                ! set indices where to insert the block
                loc_k = [(id, id = 0,X%n_mod-1)] + (grid_X%n(3)-1)*X%n_mod
                loc_m = loc_k - X%n_mod
                
                ! set the values of A (zero)
                call MatSetValues(A,X%n_mod,loc_k,X%n_mod,loc_m,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! set the values of B (zero)
                call MatSetValues(B,X%n_mod,loc_k,X%n_mod,loc_m,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! -------------------------!
                ! BLOCKS (n-2,n-1) !
                ! -------------------------!
                ! set the values of A (zero)
                call MatSetValues(A,X%n_mod,loc_m,X%n_mod,loc_k,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
                
                ! set the values of B (zero)
                call MatSetValues(B,X%n_mod,loc_m,X%n_mod,loc_k,loc_block_zero,&
                    &INSERT_VALUES,ierr)
                CHCKERR('')
            end if
            
            ! assemble the matrices
            call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t begin assembly of matrix')
            call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t begin assembly of matrix')
            call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t end assembly of matrix')
            call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t end assembly of matrix')
        end function set_mat_BC
    end function setup_matrices
    
    ! sets up EV solver
    integer function setup_solver(grid_X,X,A,B,solver) result(ierr)
        use num_vars, only: n_sol_requested
        
        character(*), parameter :: rout_name = 'setup_solver'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
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
        
        ! search for  Eigenvalue with smallest real value (so  imaginary part of
        ! sqrt(EV) = omega is largest)
        call EPSSetWhichEigenpairs(solver,EPS_SMALLEST_REAL,ierr)
        CHCKERR('Failed to set which eigenpairs')
        
        ! request n_sol_requested Eigenpairs
        if (n_sol_requested.gt.grid_X%n(3)*X%n_mod) then
            call writo('WARNING: max. nr. of solutions requested capped to &
                &problem dimension ('//trim(i2str(grid_X%n(3)*X%n_mod))//')')
            call writo('Increase either min_n_r_X or number of pol. modes or &
                &decrease n_sol_requested')
            n_sol = grid_X%n(3)*X%n_mod
        else
            n_sol = n_sol_requested
        end if
        call EPSSetDimensions(solver,n_sol,PETSC_DECIDE,&
            &PETSC_DECIDE,ierr)
        
        call lvl_ud(-1)
    end function setup_solver
    
    ! sets up guess in solver
    subroutine setup_guess(X,A,solver,start_id,prev_n_EV)
        ! input / output
        type(X_type), intent(in) :: X                                           ! perturbation variables
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
        
        ! set guess for EV if X vec is allocated and use_guess is true
        if (allocated(X%vec)) then
            ! allocate guess vectors
            allocate(guess_vec(prev_n_EV))
            
            ! create vecctor guess_vec
            do kd = 1,prev_n_EV
                call MatGetVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,istat)        ! get compatible parallel vectors to matrix A
                call VecSetOption(guess_vec(kd),&
                    &VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,istat)              ! ignore negative indices for compacter notation
            end do
            
            ! set the indices and nr. of values to set in guess_vec
            n_prev_guess = size(X%vec,2)
            allocate(guess_id(n_prev_guess*X%n_mod))
            do id = 0,n_prev_guess-1
                do jd = 1,X%n_mod
                    guess_id(id*X%n_mod+jd) = (start_id-1)*X%n_mod*2 + &
                        &id*X%n_mod*2 + jd
                end do
            end do
            
            do kd = 1,prev_n_EV
                ! set even values of guess_vec
                call VecSetValues(guess_vec(kd),n_prev_guess*X%n_mod,&
                    &guess_id-1,X%vec(:,:,kd),INSERT_VALUES,istat)
                ! set odd values of guess_vec
                call VecSetValues(guess_vec(kd),n_prev_guess*X%n_mod,&
                    &guess_id-1-X%n_mod,X%vec(:,:,kd),INSERT_VALUES,istat)
                    
                ! assemble the guess vector
                call VecAssemblyBegin(guess_vec(kd),istat)
                call VecAssemblyEnd(guess_vec(kd),istat)
            end do
                
            ! set guess
            call EPSSetInitialSpace(solver,prev_n_EV,guess_vec,istat)
            
            !! visualize guess
            !if (allocated(X%vec) .and. id.le.prev_n_EV) then
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
        character(*), parameter :: rout_name = 'get_solution'
        
        ! input / output
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
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
    integer function store_results(grid_X,X,solver,max_n_EV) result(ierr)
        use eq_vars, only: T_0
        use num_vars, only: use_normalization, grp_n_procs, grp_rank, EV_BC
        
        character(*), parameter :: rout_name = 'store_results'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(inout) :: max_n_EV                                     ! nr. of EV's saved, up to n_conv
        
        ! local variables
        PetscInt :: id, id_tot                                                  ! counters
        Vec :: sol_vec                                                          ! solution EV parallel vector
        PetscReal :: error                                                      ! error of EPS solver
        PetscInt :: one = 1                                                     ! one
        PetscReal, parameter :: infinity = 1.E40_dp                             ! beyond this value, modes are not saved
        PetscInt :: grp_n_r_X_loc                                               ! local copy of grp_n_r of X grid
        PetscScalar :: EV_loc                                                   ! local copy of EV
        PetscReal :: tol_complex = 1.E-5_dp                                     ! tolerance for an EV to be considered complex
        PetscReal :: tol_EV_BC = 1.E-3_dp                                       ! tolerance for an EV to be considered due to the BC
        PetscReal :: EV_BC_loc                                                  ! local copy of EV_BC
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! create solution vector
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,grp_n_r_X_loc*X%n_mod,&
            &grid_X%n(3)*X%n_mod,PETSC_NULL_SCALAR,sol_vec,ierr)
        CHCKERR('Failed to create MPI vector with arrays')
        
        ! set local EV_BC and user output
        if (use_normalization) then
            call writo('inverse normalization factor:')
            call lvl_ud(1)
            call writo('omega_0 = 1/T_0^2 = '//trim(r2strt(1/(T_0**2)))//&
                &' s^-2')
            EV_BC_loc = EV_BC/(T_0**2)
            call lvl_ud(-1)
        else
            EV_BC_loc = EV_BC
        end if
        
        ! allocate vec
        if (allocated(X%vec)) deallocate(X%vec)
        allocate(X%vec(1:X%n_mod,1:grp_n_r_X_loc,1:max_n_EV))
        X%vec = 0.0_dp
        
        ! allocate val
        if (allocated(X%val)) deallocate(X%val)
        allocate(X%val(1:max_n_EV))
        X%val = 0.0_dp
        
        ! start id
        id = 1
        id_tot = 1
        
        ! store them
        do while (id.le.max_n_EV)
            ! get EV solution in vector X vec
            call VecPlaceArray(sol_vec,X%vec(:,:,id),ierr)                      ! place array sol_vec in X vec
            CHCKERR('Failed to place array')
            call EPSGetEigenpair(solver,id_tot-1,X%val(id),PETSC_NULL_OBJECT,&
                &sol_vec,PETSC_NULL_OBJECT,ierr)                                ! get solution EV vector and value (starts at index 0)
            CHCKERR('EPSGetEigenpair failed')
            call EPSComputeRelativeError(solver,id_tot-1,error,ierr)            ! get error (starts at index 0)
            CHCKERR('EPSComputeRelativeError failed')
            
            ! user output
            call writo(trim(i2str(id))//'/'//trim(i2str(max_n_EV))//':')
            call lvl_ud(1)
            
            ! go back to  physical values from normalization  by multiplying the
            ! Eigenvalues by 1/T_0^2
            if (use_normalization) then
                EV_loc = X%val(id)
                X%val(id) = X%val(id)/(T_0**2)
            end if
            
            ! print output
            call writo('eigenvalue: '//trim(r2strt(realpart(X%val(id))))//&
                &' + '//trim(r2strt(imagpart(X%val(id))))//' i')
            if (use_normalization) call writo('(normalized: '//&
                &trim(r2strt(realpart(EV_loc)))//' + '//&
                &trim(r2strt(imagpart(EV_loc)))//' i)')
            call writo('with rel. error: '//trim(r2str(error)))
            
            ! tests
            if (abs(imagpart(X%val(id))/realpart(X%val(id))).gt.tol_complex) &
                &then                                                           ! test for unphysical complex solution
                call writo('WARNING: Unphysical complex Eigenvalue!')
                call remove_EV(id,max_n_EV)
            else if (abs(realpart(X%val(id)-EV_BC_loc)/EV_BC_loc).lt.&
                &tol_EV_BC) then                                                ! test for artificial EV due to BC
                call writo('WARNING: Eigenvalue probably due to BC''s &
                    &(artifically set to '//trim(r2strt(EV_BC))//')')
                call remove_EV(id,max_n_EV)
            else if (abs(X%val(id)).gt.infinity) then                           ! test for infinity
                call writo('WARNING: The next Eigenvalues are larger than '//&
                    &trim(r2strt(infinity))//' and are discarded')
                call remove_EV(id,max_n_EV,remove_next=.true.)
                exit
            else                                                                ! tests passed
                ! increment id
                id = id + 1
            end if
            
            !! visualize solution
            !call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,&
                !&'solution '//trim(i2str(id-1))//' vector with '//&
                !&trim(i2str(grid_X%n(3)))//' points',0,&
                !&0,1000,1000,guess_viewer,istat)
            !call VecView(sol_vec,guess_viewer,istat)
            
            !call writo('Eigenvector = ')
            !call VecView(sol_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !CHCKERR('Cannot view vector')
            
            call lvl_ud(-1)
            
            ! reset vector
            call VecResetArray(sol_vec,ierr)
            CHCKERR('Failed to reset array')
            
            ! increment total id
            id_tot = id_tot+1
        end do
        
        !call writo('Summary:')
        !call EPSPrintSolution(solver,PETSC_NULL_OBJECT,ierr)
        
        call VecDestroy(sol_vec,ierr)                                           ! destroy compatible vector
        CHCKERR('Failed to destroy sol_vec')
        
        call lvl_ud(-1)
    contains
        ! Removes id'th  Eigenvalue and  -vector. Optionally, all  the following
        ! can be removed as well.
        subroutine remove_EV(id,max_id,remove_next)
            use num_vars, only: retain_all_sol
            
            ! input / output
            PetscInt, intent(inout) :: id                                       ! id of faulty values
            PetscInt, intent(inout) :: max_id                                   ! maximum id
            PetscBool, intent(in), optional :: remove_next                      ! whether all next values have to be removed as well
            
            ! local variables
            PetscScalar, allocatable :: X_val_loc(:)                            ! local copy of X_val
            PetscScalar, allocatable :: X_vec_loc(:,:,:)                        ! local copy of X_vec
            PetscBool :: remove_next_loc = .false.                              ! local copy of remove_next
            
            ! only remove if faulty solutions are not optionally retained
            if (retain_all_sol) then
                id = id+1                                                       ! increment the solution
            else
                ! set up local remove_next
                if (present(remove_next)) remove_next_loc = remove_next
                
                ! user output
                call lvl_ud(1)
                if (remove_next_loc) then
                    call writo('All next solutions are removed')
                else
                    call writo('Current solution is removed')
                end if
                call writo('(Override using "retain_all_sol")')
                
                ! save old arrays
                allocate(X_val_loc(max_id))
                allocate(X_vec_loc(X%n_mod,grp_n_r_X_loc,max_id))
                X_val_loc = X%val
                X_vec_loc = X%vec
                
                ! reallocate arrays
                deallocate(X%val,X%vec)
                if (remove_next_loc) then
                    allocate(X%val(id-1))
                    allocate(X%vec(X%n_mod,grp_n_r_X_loc,id-1))
                else
                    allocate(X%val(max_id-1))
                    allocate(X%vec(X%n_mod,grp_n_r_X_loc,max_id-1))
                end if
                
                ! copy other values
                X%val(1:id-1) = X_val_loc(1:id-1)
                X%vec(:,:,1:id-1) = X_vec_loc(:,:,1:id-1)
                if (.not.remove_next_loc) then
                    X%val(id:max_id-1) = X_val_loc(id+1:max_id)
                    X%vec(:,:,id:max_id-1) = X_vec_loc(:,:,id+1:max_id)
                end if
                
                ! adapt max_id
                if (remove_next_loc) then
                    max_id = id-1
                else
                    max_id = max_id-1
                end if
                
                ! clean up
                deallocate(X_val_loc,X_vec_loc)
                
                call lvl_ud(-1)
            end if
        end subroutine remove_EV
    end function store_results
    
    ! stop PETSC and SLEPC
    ! [MPI] Collective call
    integer function stop_SLEPC(grid_X,A,B,solver,guess_start_id,prev_n_EV,&
        &max_n_EV) result(ierr)
        character(*), parameter :: rout_name = 'stop_SLEPC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        Mat, intent(in) :: A, B                                                 ! matrices A and B in EV problem A X = lambda B X
        EPS, intent(in) :: solver                                               ! EV solver
        PetscInt, intent(inout):: guess_start_id                                ! start of index of vector, saved for next iteration
        PetscInt, intent(inout) :: prev_n_EV                                    ! nr. of solutions of previous vector, saved for next iteration
        PetscInt, intent(in) :: max_n_EV                                        ! nr. of EV's saved
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set guess_start_id and prev_n_EV
        guess_start_id = grid_X%i_min
        prev_n_EV = max_n_EV
        
        ! destroy objects
        call EPSDestroy(solver,ierr)
        CHCKERR('Failed to destroy EPS solver')
        call MatDestroy(A,ierr)
        CHCKERR('Failed to destroy matrix A')
        call MatDestroy(B,ierr)
        CHCKERR('Failed to destroy matrix B')
        
        ! stop SLEPC
        call SlepcFinalize(ierr)
        CHCKERR('Failed to Finalize SLEPC')
        
        call lvl_ud(-1)
    end function stop_SLEPC
end module
