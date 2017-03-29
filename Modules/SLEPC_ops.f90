!------------------------------------------------------------------------------!
!   Operations that use SLEPC (and PETSC) routines                             !
!   Note that the routines here require a TRIMMED solution grid!               !
!------------------------------------------------------------------------------!
! References:                                                                  !
! [1]   http://slepc.upv.es/documentation/current/docs/manualpages/EPS/        !
!       EPSComputeRelativeError.html                                           !
!------------------------------------------------------------------------------!
module SLEPC_ops
#include <PB3D_macros.h>
#include <wrappers.h>
! for slepc 3.6.0:
#include <slepc/finclude/slepcepsdef.h>
! for slepc 3.5.3:
!#include <finclude/slepcepsdef.h>
    use str_utilities
    use output_ops
    use messages
    use slepceps
    use num_vars, only: iu, dp, max_str_ln
    use grid_vars, only: grid_type
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type

    implicit none
    private
    public solve_EV_system_SLEPC
#if ldebug
    public debug_setup_mats, debug_set_BC
#endif
    
    ! global variables
    PetscInt :: n_r                                                             ! n_r of solution
    PetscInt :: loc_n_r                                                         ! loc_n_r of solution
    integer, allocatable :: c_tot(:,:,:)                                        ! total c corresponding to symmetric and asymmetric tensorial variables
    logical :: use_hermitian = .false.                                          ! use hermitian matrices (does not work well currently, see http://lists.mcs.anl.gov/pipermail/petsc-users/2016-October/030781.html)
#if ldebug
    logical :: debug_setup_mats = .false.                                       ! plot debug information for setup_mats
    logical :: debug_set_BC = .false.                                           ! plot debug information for set_BC
    logical :: debug_calc_V_0_mod = .false.                                     ! plot debug information for calc_V_0_mod
    logical :: test_diff = .false.                                              ! test introduction of numerical diff
    real(dp) :: diff_coeff                                                      ! diff coefficient
#endif
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in  [ADD REF] and  solves them using  the SLEPC suite.  The most
    ! unstable solutions  are obtained.
    ! Optionally the  geodesic index of  the perturbation variables at  which to
    ! perform the calculations can be changed  from its default value of 1 using
    ! the variable i_geo.
    ! Note: the  perturbation grid needs  to be the  same as the  solution grid,
    ! if not,  this module has to  be adapted to interpolate  them. A deprecated
    ! technique,  using "get_norm_interp_data"  and  "interp_V" can  be seen  in
    ! builds previous to 1.06. These have been superseded by "setup_interp_data"
    ! followed by "apply_disc".
    integer function solve_EV_system_SLEPC(grid_X,grid_sol,X,sol,i_geo) &
        &result(ierr)
        
        use num_vars, only: max_it_inv, ndps => norm_disc_prec_sol, &
            &matrix_SLEPC_style
        use num_utilities, only: calc_coeff_fin_diff
        use rich_vars, only: use_guess
#if ldebug
        use num_vars, only: ltest
        use input_utilities, only: get_real, get_log
#endif
        
        character(*), parameter :: rout_name = 'solve_EV_system_SLEPC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        integer, intent(in), optional :: i_geo                                  ! at which geodesic index to perform the calculations
        
        ! local variables
        Mat, target :: A                                                        ! matrix A in EV problem A X = lambda B X
        Mat, target :: B                                                        ! matrix B in EV problem A X = lambda B X
        EPS :: solver                                                           ! EV solver
        PetscInt :: max_n_EV                                                    ! how many solutions saved
        PetscInt :: inv_lvl_nr                                                  ! level of inverse calculation
        PetscInt :: i_geo_loc                                                   ! local copy of i_geo
        PetscReal, allocatable :: norm_disc_coeff(:)                            ! discretization coefficients
        PetscReal :: step_size                                                  ! step size in flux coordinates
        PetscBool :: done_inverse                                               ! is it converged?
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialization
        call writo('Initialize generalized Eigenvalue problem')
        
        ! set up local i_geo
        i_geo_loc = 1
        if (present(i_geo)) i_geo_loc = i_geo
        
        ! start SLEPC
        ierr = start_SLEPC(grid_sol%n(3))
        CHCKERR('')
        
        ! set up step size
        step_size = (grid_sol%r_F(grid_sol%n(3))-grid_sol%r_F(1))/&
            &(grid_sol%n(3)-1)
        
        ! get discretization coefficients
        call writo('Get discretization coefficients')
        call lvl_ud(1)
        ierr = calc_coeff_fin_diff(1,ndps,norm_disc_coeff)
        CHCKERR('')
        norm_disc_coeff = norm_disc_coeff/step_size                             ! scale by step size
        call lvl_ud(-1)
        
        ! set up the matrix
        call writo('Set up matrices')
        call lvl_ud(1)
        ierr = setup_mats(grid_X,grid_sol,X,A,B,i_geo_loc,norm_disc_coeff)
        CHCKERR('')
        call lvl_ud(-1)
        
#if ldebug
        if (debug_setup_mats) then
            call writo('Testing if A and B are Hermitian by multiplying with &
                &their Hermitian Transpose and subtracting 1')
            call lvl_ud(1)
                
            ierr = test_mat_hermiticity(A,'A')
            CHCKERR('')
            
            ierr = test_mat_hermiticity(B,'B')
            CHCKERR('')
            
            call lvl_ud(-1)
        end if
#endif
        
        ! iterate over matrix inverse
        if (max_it_inv.gt.1) then
            call writo('Iterating over the inverse')
            call lvl_ud(1)
        end if
        
        done_inverse = .false.
        inv_lvl_nr = 1
        Inverse: do while (.not.done_inverse .and. inv_lvl_nr.le.max_it_inv)
            if (max_it_inv.gt.1) then
                call writo('Iteration '//trim(i2str(inv_lvl_nr))//&
                    &' of calculation of inverse')
                call lvl_ud(1)
            end if
            
            ! set boundary conditions
            call writo('Set up boundary conditions')
            call lvl_ud(1)
            
            select case (matrix_SLEPC_style)
                case (1)                                                        ! sparse
                    ierr = set_BC(grid_X,X,A,B,i_geo_loc,grid_sol%n(3),&
                        &norm_disc_coeff)
                    CHCKERR('')
            
#if ldebug
                    if (debug_set_BC) then
                        call writo('Testing if AFTER INSERTING BC''s, A and B &
                            &are Hermitian by multiplying with their &
                            &Hermitian Transpose and subtracting 1')
                        call lvl_ud(1)
                        
                        ierr = test_mat_hermiticity(A,'A_BC')
                        CHCKERR('')
                        
                        ierr = test_mat_hermiticity(B,'B_BC')
                        CHCKERR('')
                        
                        call lvl_ud(-1)
                    end if
#endif
                case (2)                                                        ! shell
                    ierr = 1
                    err_msg = 'NOT YET IMPLEMENTED FOR SHELL MATRICES'
                    CHCKERR(err_msg)
            end select
            
            call lvl_ud(-1)
            
            ! set up solver
            call writo('Set up EV solver')
#if ldebug
            if (ltest) then
                call writo('Test spectrum of A or B instead of solving &
                    &generalized Eigenvalue problem?')
                if (get_log(.false.)) then
                    call writo('Spectrum of A (true) or B (false)?')
                    if (get_log(.true.)) then                                   ! A
                        ierr = setup_solver(X,A,PETSC_NULL_OBJECT,solver)
                        CHCKERR('')
                    else                                                        ! B
                        ierr = setup_solver(X,B,PETSC_NULL_OBJECT,solver)
                        CHCKERR('')
                    end if
                else
                    ierr = setup_solver(X,A,B,solver)
                    CHCKERR('')
                end if
            else
#endif
                ierr = setup_solver(X,A,B,solver)
                CHCKERR('')
#if ldebug
            end if
#endif
            
            ! set up guess
            if (use_guess) then
                call writo('Set up guess')
                call lvl_ud(1)
                ierr = setup_guess(sol,A,solver)
                CHCKERR('')
                call lvl_ud(-1)
            end if
            
            ! get solution
            call writo('Get solution')
            
            ierr = get_solution(solver)
            CHCKERR('')
            
            ! check for convergence
            inv_lvl_nr = inv_lvl_nr+1
            !!! TO BE IMPLEMENTED !!!
            !!! ALSO, GUESS SHOULD BE ADAPTED TO PREVIOUS GUESS PROBABLY ???
            
            if (max_it_inv.gt.1) call lvl_ud(-1)
        end do Inverse
        
        if (max_it_inv.gt.1) then
            call lvl_ud(-1)
            call writo('Done iterating over the inverse')
        end if
        
        ! summarize solution
        call writo('Summarize solution')
        
        ierr = summarize_solution(solver,max_n_EV)
        CHCKERR('')
        
        ! store results
        call writo('Store results for '//trim(i2str(max_n_EV))//' least &
            &stable Eigenvalues')
        
        ierr = store_results(grid_sol,sol,solver,max_n_EV,A,B,step_size)
        CHCKERR('')
        
        ! clean up
        deallocate(norm_disc_coeff)
        
        ! finalize
        call writo('Finalize SLEPC')
        
        ! stop SLEPC
        ierr = stop_SLEPC(A,B,solver)
        CHCKERR('')
        
#if ldebug
    contains
        integer function test_mat_hermiticity(mat,mat_name) result(ierr)
            use num_vars, only: data_dir
            !use num_vars, only: rank
            
            character(*), parameter :: rout_name = 'test_mat_hermiticity'
            
            ! input / output
            Mat, intent(in) :: mat                                              ! matrix to test
            character(len=*) :: mat_name                                        ! name of matrix
            
            ! local variables
            Mat :: mat_t                                                        ! Hermitian transpose of mat
            Mat :: mat_loc                                                      ! copy of mat, but without block storage (Hermitian transpose does not work for block)
            PetscScalar :: one = 1.0                                            ! one
            PetscViewer :: file_viewer                                          ! viewer to write matrix to file
            character(len=max_str_ln) :: file_name                              ! file name
            
            ! initialize ierr
            ierr = 0
            
            ! duplicate mat into mat_loc
            call MatConvert(mat,MATMPIAIJ,MAT_INITIAL_MATRIX,mat_loc,ierr)
            CHCKERR('Failed to convert mat into mat_loc')
            
            ! test if mat is hermitian
            call MatHermitianTranspose(mat_loc,MAT_INITIAL_MATRIX,mat_t,ierr)
            CHCKERR('Hermitian transpose of mat failed')
            call MatAXPY(mat_t,-one,mat_loc,SAME_NONZERO_PATTERN,ierr)
            CHCKERR('mat-mat_t failed')
            
            ! visualize the matrices
            call PetscOptionsSetValue('-draw_pause','-1',ierr)
            call writo(trim(mat_name)//' =')
            call MatView(mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !call MatView(mat,PETSC_VIEWER_DRAW_WORLD,ierr)
            call writo(trim(mat_name)//' - '//trim(mat_name)//'^*T =')
            call MatView(mat_t,PETSC_VIEWER_STDOUT_WORLD,ierr)
            
            ! destroy matrices 
            call MatDestroy(mat_t,ierr)
            CHCKERR('Failed to destroy matrix mat_t')
            call MatDestroy(mat_loc,ierr)
            CHCKERR('Failed to destroy matrix mat_loc')
            
            !! duplicate mat into mat_loc
            !call MatDuplicate(mat,MAT_SHARE_NONZERO_PATTERN,mat_loc,ierr)
            !CHCKERR('failed to duplicate mat into mat_loc')
            
            ! write to file
            !file_name = trim(data_dir)//'/'//trim(mat_name)//'_RE'
            file_name = trim(data_dir)//'/'//trim(mat_name)
            call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(file_name),&
                &file_viewer,ierr)
            CHCKERR('Unable to open file viewer')
            !call PetscViewerSetFormat(file_viewer,PETSC_VIEWER_ASCII_DENSE,ierr)
            call PetscViewerSetFormat(file_viewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
            CHCKERR('Unable to set format')
            !call MatView(mat_loc,file_viewer,ierr)
            call MatView(mat,file_viewer,ierr)
            CHCKERR('Unable to write matrix to file')
            call PetscViewerDestroy(file_viewer,ierr)
            CHCKERR('Unable to destroy file viewer')
            !if (rank.eq.0) then
                !call execute_command_line('tail -n +3 '//trim(file_name)//&
                    !&' > tempfile && mv tempfile '//trim(file_name),&
                    !&EXITSTAT=ierr)
                !CHCKERR('Failed to execute command')
            !end if
            
            !! write imaginary part to file
            !file_name = trim(data_dir)//'/'//trim(mat_name)//'_IM'
            !call MatCopy(mat,mat_loc,MAT_SHARE_NONZERO_PATTERN,ierr)            ! mat_loc has same structure as mat
            !CHCKERR('Failed to copy mat into mat_loc')
            !call MatImaginaryPart(mat_loc,ierr)
            !CHCKERR('Failed to get imaginary part')
            !call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(file_name),&
                !&file_viewer,ierr)
            !CHCKERR('Unable to open file viewer')
            !call PetscViewerSetFormat(file_viewer,PETSC_VIEWER_ASCII_DENSE,ierr)
            !CHCKERR('Unable to set format')
            !call MatView(mat_loc,file_viewer,ierr)
            !CHCKERR('Unable to write matrix to file')
            !call PetscViewerDestroy(file_viewer,ierr)
            !CHCKERR('Unable to destroy file viewer')
            !if (rank.eq.0) then
                !call execute_command_line('tail -n +3 '//trim(file_name)//&
                    !&' > tempfile && mv tempfile '//trim(file_name),&
                    !&EXITSTAT=ierr)
                !CHCKERR('Failed to execute command')
            !end if
            
            !! destroy matrices
            !call MatDestroy(mat_loc,ierr)
            !CHCKERR('Failed to destroy matrix mat_loc')
        end function test_mat_hermiticity
#endif
    end function solve_EV_system_SLEPC
    
    !  This  subroutine starts  PETSC  and  SLEPC  with  the correct  number  of
    ! processors
    integer function start_SLEPC(n_r_sol) result(ierr)
        use num_vars, only: n_procs, sol_n_procs
        use files_ops, only: opt_args
        
        character(*), parameter :: rout_name = 'start_SLEPC'
        
        ! input / output
        integer, intent(in) :: n_r_sol                                          ! nr. of normal points in solution grid
        
        ! local variables
        PetscBool :: flg                                                        ! flag to catch options
        PetscInt :: id                                                          ! counters
        character(len=max_str_ln) :: option_name                                ! options
#if !defined(PETSC_USE_COMPLEX)
        character(len=max_str_ln) :: err_msg                                    ! error message
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! user message
        call writo('initialize SLEPC...')
        call lvl_ud(1)
        
        ! use MPI_Comm_world for PETSC_COMM_WORLD
        PETSC_COMM_WORLD = MPI_Comm_world
        if (n_procs.gt.n_r_sol) then                                            ! too many processors
            call writo('using too many processors: '//trim(i2str(n_procs))//&
                &', while beyond '//trim(i2str(n_r_sol))//&
                &' does not bring improvement',warning=.true.)
        end if
        call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)                         ! initialize SLEPC
        CHCKERR('SLEPC failed to initialize')
        
        ! output
        call writo('slepc started with '//trim(i2str(n_procs))&
            &//' processes')
        if (sol_n_procs.lt.n_procs) call writo('(but only '//&
            &trim(i2str(sol_n_procs))//' processes will be used)')
        
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
    ! The  geodesical index  at  which to  perform the  calculations  has to  be
    ! provided as well.
    ! The matrices are set up in the  solution grid, so some processes might not
    ! have any part in the storage  thereof (if sol_n_procs < n_procs), but each
    ! process sets up the part of the  grid that corresponds to their own normal
    ! range in the perturbation grid.
    ! Note: For normal usage, i_geo should be 1, or not present.
    integer function setup_mats(grid_X,grid_sol,X,A,B,i_geo,norm_disc_coeff) &
        &result(ierr)
        
        use num_vars, only: ndps => norm_disc_prec_sol, matrix_SLEPC_style, &
            &rank, sol_n_procs
        use grid_utilities, only: trim_grid
        use SLEPC_utilities, only: insert_block_mat
        use X_vars, only: n_mod_X
#if ldebug
        use num_vars, only: ltest
        use input_utilities, only: get_real, get_log
#endif
        
        use num_utilities, only: c
        
        character(*), parameter :: rout_name = 'setup_mats'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in), target :: X                                 ! field-averaged perturbation variables (so only first index)
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        integer, intent(in) :: i_geo                                            ! at which geodesic index to perform the calculations
        PetscReal, intent(in) :: norm_disc_coeff(:)                             ! discretization coefficients for normal derivatives
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed solution grid
        type(grid_type) :: grid_X_trim                                          ! trimmed perturbation grid
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed perturbation grid
        PetscInt :: kd, id                                                      ! counter
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! local variables also used in child routines
        PetscInt :: bulk_i_lim(2)                                               ! absolute limits of bulk matrix (excluding the BC's)
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Normal discretization with central finite differences of &
            &order '//trim(i2str(ndps))//', stencil width '//&
            &trim(i2str(4*ndps+1)))
        
#if ldebug
        if (ltest) then
            call writo('Test introduction of numerical friction?')
            if (get_log(.false.)) then
                test_diff = .true.
                call writo('The friction coefficient will be multiplied with &
                    &the central block, and inserted in the off-diagonal parts')
                call writo('How large should it be?')
                if (get_log(.false.)) diff_coeff = get_real(lim_lo=0._dp)
            end if
        end if
#endif
        
        ! check nr. of modes
        if (n_mod_X.ne.X%n_mod(1) .or. n_mod_X.ne.X%n_mod(2)) then
            ierr = 1
            err_msg = 'Need square matrix of size [n_mod_X:n_mod_X]'
            CHCKERR(err_msg)
        end if
        
        ! setup total c if not yet done
        if (.not.allocated(c_tot)) then
            allocate(c_tot(n_mod_X,n_mod_X,2))
            do kd = 1,n_mod_X
                do id = 1,n_mod_X
                    c_tot(id,kd,1) = c([id,kd],.true.,n_mod_X)
                    c_tot(id,kd,2) = c([id,kd],.false.,n_mod_X)
                end do
            end do
        end if
        
        ! trim grids
        ierr = trim_grid(grid_sol,grid_sol_trim)
        CHCKERR('')
        ierr = trim_grid(grid_X,grid_X_trim,norm_id)
        CHCKERR('')
        
        ! set up loc_n_r and n_r
        if (rank.lt.sol_n_procs) then
            loc_n_r = grid_sol_trim%loc_n_r
        else
            loc_n_r = 0
        end if
        n_r = grid_sol_trim%n(3)
        
        ! set up bulk matrix absolute limits
        ierr = set_bulk_lims(grid_X_trim,bulk_i_lim)
        CHCKERR('')
        
        ! setup matrix A and B
        select case (matrix_SLEPC_style)
            case (1)                                                            ! sparse
                ! calculate nonzeros
                call set_nonzeros()
                
                ! create matrix A
                call MatCreateBAIJ(PETSC_COMM_WORLD,n_mod_X,n_mod_X*loc_n_r,&
                    &n_mod_X*loc_n_r,n_mod_X*n_r,n_mod_X*n_r,&
                    &PETSC_NULL_INTEGER,d_nz,PETSC_NULL_INTEGER,o_nz,A,ierr)
                err_msg = 'MatCreateAIJ failed for matrix A'
                CHCKERR(err_msg)
                
                ! deallocate tot_nz, d_nz and o_nz
                deallocate(d_nz,o_nz)
                
                ! fill the matrix A
                ierr = fill_mat(X%PV_0(1,i_geo,norm_id(1):norm_id(2),:),&
                    &X%PV_1(1,i_geo,norm_id(1):norm_id(2),:),&
                    &X%PV_2(1,i_geo,norm_id(1):norm_id(2),:),bulk_i_lim,A)
                CHCKERR('')
                
                call writo('matrix A set up:')
                ierr = disp_mat_info(A)
                CHCKERR('')
                
                ! duplicate A into B
                call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)           ! B has same structure as A
                CHCKERR('failed to duplicate A into B')
                
                ! fill the matrix B
                ierr = fill_mat(X%KV_0(1,i_geo,norm_id(1):norm_id(2),:),&
                    &X%KV_1(1,i_geo,norm_id(1):norm_id(2),:),&
                    &X%KV_2(1,i_geo,norm_id(1):norm_id(2),:),bulk_i_lim,B)
                CHCKERR('')
                
                call writo('matrix B set up:')
                ierr = disp_mat_info(B)
                CHCKERR('')
            case (2)                                                            ! shell
                ierr = 1
                err_msg = 'NOT YET IMPLEMENTED FOR SHELL MATRICES'
                CHCKERR(err_msg)
        end select
        
        ! clean up
        call grid_X_trim%dealloc()
        call grid_sol_trim%dealloc()
    contains
        ! Sets the limits of the indices of the bulk matrix, depending on the BC
        ! style.
        integer function set_bulk_lims(grid_X,i_lim) result(ierr)
            use num_vars, only: BC_style
            
            character(*), parameter :: rout_name = 'set_bulk_lims'
            
            ! input / output
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            integer, intent(inout) :: i_lim(2)                                  ! min and max of bulk limits
            
            ! initialize ierr
            ierr = 0
            
            ! set up i_min
            select case (BC_style(1))
                case (1)
                    i_lim(1) = grid_X%i_min                                     ! will be overwritten
                case (2)
                    i_lim(1) = max(grid_X%i_min,1+ndps)                        ! first norm_disc_prec_sol rows write left BC's
                case (3)
                    err_msg = 'Left BC''s cannot have BC type 3'
                    ierr = 1
                    CHCKERR(err_msg)
                case default
                    err_msg = 'No BC style associated with '//&
                        &trim(i2str(BC_style(1)))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! set up i_max
            select case (BC_style(2))
                case (1)
                    i_lim(2) = grid_X%i_max                                     ! will be overwritten
                case (2)
                    i_lim(2) = min(grid_X%i_max,n_r-ndps)                       ! last norm_disc_prec_sol rows write right BC's
                case (3)
                    i_lim(2) = min(grid_X%i_max,n_r-1)                          ! last row writes right BC's
                case default
                    err_msg = 'No BC style associated with '//&
                        &trim(i2str(BC_style(1)))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end function set_bulk_lims
        
        ! Fills a  matrix according to norm_disc_prec_sol [ADD REF]:
        !   1. first order accuracy for all terms
        !   2. higher order accuracy for internal first order derivatives
        ! it is used  for both the matrix  A and B, corresponding  to the plasma
        ! potential energy and the (perpendicular) kinetic energy
        ! The procedure is as follows:
        !   1. Tables  are set  up for  (k,m) pairs in the equilibrium grid.
        !   2. At the first normal point  in the solution grid, belonging to the
        !   current process, the tables are interpolated as a start.
        !   3. At  every normal  point in  the solution  grid, belonging  to the
        !   current process,  the tables are  interpolated at the next  point in
        !   the  corresponding  perturbation  grid. This  information,  as  well
        !   as  the interpolated  value  of  the previous  point  allow for  the
        !   calculation of every quantity.
        ! Makes use of n_r, norm_disc_coeff, grid_X_trim and  grid_sol_trim
        integer function fill_mat(V_0,V_1,V_2,bulk_i_lim,mat) result(ierr)
            use num_utilities, only: con
            
            character(*), parameter :: rout_name = 'fill_mat'
            
            ! input / output
            PetscScalar, intent(in) :: V_0(:,:)                                 ! either PV_int_0 or KV_int_0 in equilibrium normal grid
            PetscScalar, intent(in) :: V_1(:,:)                                 ! either PV_int_1 or KV_int_1 in equilibrium normal grid
            PetscScalar, intent(in) :: V_2(:,:)                                 ! either PV_int_2 or KV_int_2 in equilibrium normal grid
            PetscInt, intent(in) :: bulk_i_lim(2)                               ! absolute limits of bulk matrix (excluding the BC's)
            Mat, intent(inout) :: mat                                           ! either A or B
            
            ! local variables (not to be used in child routines)
            character(len=max_str_ln) :: err_msg                                ! error message
            PetscScalar, allocatable :: loc_block(:,:)                          ! [n_mod_X:n_mod_X] block matrix for 1 normal point
            PetscInt :: id, jd, kd                                              ! counters
            PetscInt :: kd_loc                                                  ! kd in local variables
            PetscInt :: k, m                                                    ! counters
#if ldebug
            PetscScalar, allocatable :: loc_block_0_backup(:,:)                 ! backed up V_0 local block
#endif
            
            ! for tests
            PetscInt :: r_sol_start, r_sol_end                                  ! start block and end block
            
            ! initialize ierr
            ierr = 0
            
            ! test whether the matrix range coincides with i_min and i_max
            call MatGetOwnershipRange(mat,r_sol_start,r_sol_end,ierr)           ! starting and ending row r_sol_start and r_sol_end
            err_msg = 'Couldn''t get ownership range of matrix'
            CHCKERR(err_msg)
            r_sol_start = r_sol_start/n_mod_X                                   ! count per block
            r_sol_end = r_sol_end/n_mod_X                                       ! count per block
            if (rank.lt.sol_n_procs) then
                if (grid_sol_trim%i_min.ne.r_sol_start+1) then
                    ierr = 1
                    err_msg = 'start of matrix in this process does not &
                        &coincide with tabulated value'
                    CHCKERR(err_msg)
                end if
                if (grid_sol_trim%i_max.ne.r_sol_end) then
                    ierr = 1
                    err_msg = 'end of matrix in this process does not &
                        &coincide with tabulated value'
                    CHCKERR(err_msg)
                end if
            end if
            
            ! allocate local block
            allocate(loc_block(n_mod_X,n_mod_X))
            
            ! iterate over all rows of this rank
            do kd = bulk_i_lim(1)-1,bulk_i_lim(2)-1                             ! (indices start with 0 here)
                ! set up local kd
                kd_loc = kd+2-grid_X_trim%i_min
                
                ! fill the blocks
                
                ! -------------!
                ! BLOCKS ~ V_0 !
                ! -------------!
                ! fill local block
                do m = 1,n_mod_X
                    do k = 1,n_mod_X
                        loc_block(k,m) = &
                            &con(V_0(kd_loc,c_tot(k,m,1)),&
                            &[k,m],.true.)                                      ! symmetric matrices need con()
                    end do
                end do

#if ldebug
                if (test_diff) then
                    ! back up the V_0 block
                    allocate(loc_block_0_backup(n_mod_X,n_mod_X))
                    loc_block_0_backup = loc_block
                end if
#endif

                ! add block to kd + (0,0)
                ierr = insert_block_mat(loc_block,mat,kd,[0,0],n_r)
                CHCKERR('')
                
                ! -------------!
                ! BLOCKS ~ V_1 !
                ! -------------!
                ! fill local block
                do m = 1,n_mod_X
                    do k = 1,n_mod_X
                        loc_block(k,m) = V_1(kd_loc,c_tot(k,m,2))               ! asymetric matrices don't need con()
                    end do
                end do
                
                ! add block to kd + (0,-p..p) + Hermitian conjugate
                do jd = -ndps,ndps
                    ierr = insert_block_mat(loc_block*&
                        &norm_disc_coeff(jd+ndps+1),mat,kd,&
                        &[0,jd],n_r,transp=.true.)
                    CHCKERR('')
                end do
                
                ! -------------!
                ! BLOCKS ~ V_2 !
                ! -------------!
                ! fill local block
                do m = 1,n_mod_X
                    do k = 1,n_mod_X
                        loc_block(k,m) = &
                            &con(V_2(kd_loc,c_tot(k,m,1)),[k,m],.true.)         ! symmetric matrices need con()
                    end do
                end do
                
#if ldebug
                if (test_diff) then
                    loc_block = loc_block + diff_coeff*loc_block_0_backup
                    deallocate(loc_block_0_backup)
                end if
#endif
            
                ! add block to kd + (-p..p,-p..p)
                do jd = -ndps,ndps
                    do id = -ndps,ndps
                        ierr = insert_block_mat(loc_block*&
                            &norm_disc_coeff(id+ndps+1)*&
                            &norm_disc_coeff(jd+ndps+1),mat,&
                            &kd,[id,jd],n_r)
                        CHCKERR('')
                    end do
                end do
            end do
            
            ! assemble the matrix
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t begin assembly of matrix')
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t end assembly of matrix')
            
            ! Hermitian matrix
            if (use_hermitian) then
                call MatSetOption(mat,MAT_HERMITIAN,PETSC_TRUE,ierr)
                CHCKERR('Coulnd''t set option Hermitian')
            end if
            
            ! clean up
            deallocate(loc_block)
        end function fill_mat
        
        ! Display information about matrix.
        integer function disp_mat_info(mat) result(ierr)
            character(*), parameter :: rout_name = 'disp_mat_info'
            
            ! input / output
            Mat, intent(inout) :: mat                                           ! either A or B
            
            ! local variables
            real(dp) :: mat_info(MAT_INFO_SIZE)                                 ! information about matrix
            
            ! initialize ierr
            ierr = 0
            
            call lvl_ud(1)
            
            select case (matrix_SLEPC_style)
                case (1)                                                        ! sparse
                    call MatGetInfo(mat,MAT_GLOBAL_SUM,mat_info,ierr)
                    CHCKERR('')
                    call writo('memory usage: '//&
                        &trim(r2strt(mat_info(MAT_INFO_MEMORY)*1.E-6_dp))//&
                        &' MB')
                    call writo('nonzero''s allocated: '//&
                        &trim(r2strt(mat_info(MAT_INFO_NZ_ALLOCATED))))
                    if (mat_info(MAT_INFO_NZ_UNNEEDED).gt.0._dp) then
                        call writo('of which unused: '//&
                            &trim(r2strt(mat_info(MAT_INFO_NZ_UNNEEDED))))
                    end if
                case (2)                                                        ! shell
                    call writo('shell matrix')
            end select
            
            call lvl_ud(-1)
        end function disp_mat_info
        
        ! Sets nonzero elements d_nz and o_nz.
        ! Makes use of ndps and grid_sol_trim.
        subroutine set_nonzeros()
            ! local variables
            PetscInt, allocatable :: tot_nz(:)                                  ! nr. of total non-zeros
            PetscInt :: st_size                                                 ! half stencil size
            
            ! initialize the numbers of non-zeros in diagonal and off-diagonal
            allocate(tot_nz(loc_n_r)); tot_nz = 0
            allocate(d_nz(loc_n_r)); d_nz = 0
            allocate(o_nz(loc_n_r)); o_nz = 0
            
            if (rank.lt.sol_n_procs) then
                ! set (half) stencil size
                st_size = 2*ndps
                
                ! calculate number of total nonzero entries
                tot_nz = 1+2*st_size
                do kd = 1,st_size
                    ! limit due to left BC
                    if (kd.ge.grid_sol_trim%i_min .and. &
                        &kd.le.grid_sol_trim%i_max) &
                        &tot_nz(kd-grid_sol_trim%i_min+1) = &
                        &tot_nz(kd-grid_sol_trim%i_min+1) - (st_size+1-kd)
                end do
                do kd = n_r,n_r-st_size+1,-1
                    ! limit due to right BC
                    if (kd.ge.grid_sol_trim%i_min .and. &
                        &kd.le.grid_sol_trim%i_max) &
                        &tot_nz(kd-grid_sol_trim%i_min+1) = &
                        &tot_nz(kd-grid_sol_trim%i_min+1) - (kd-n_r+st_size)
                end do
                
                ! calculate number of nonzero diagonal entries
                if (loc_n_r.le.(st_size+1)) then                                ! up to (st_size+1) normal points in this process
                    d_nz = loc_n_r
                else                                                            ! more than (st_size+1) normal points in this process
                    do kd = 1,st_size
                        d_nz(kd) = st_size+kd
                        d_nz(loc_n_r-kd+1) = st_size+kd
                    end do
                    d_nz(st_size+1:loc_n_r-st_size) = 2*st_size+1
                end if
                d_nz = min(d_nz,loc_n_r)                                        ! limit to loc_n_r
                
                ! calculate number of nonzero off-diagonal entries
                do kd = 1,loc_n_r
                    o_nz(kd) = tot_nz(kd)-d_nz(kd)
                end do
            else
                d_nz = 0
                o_nz = 0
            end if
            
            ! clean up
            deallocate(tot_nz)
        end subroutine set_nonzeros
    end function setup_mats
    
    ! Sets the  boundary conditions:
    ! Deep  in the  plasma,  an  artificial Eigenvalue  EV_BC  is introduced  by
    ! setting the  diagonal components  of A  to EV_BC and  of B  to 1,  and the
    ! off-diagonal elements to zero.
    ! At the plasma surface, the surface energy is minimized as in [ADD REF].
    integer function set_BC(grid_X,X,A,B,i_geo,n_sol,norm_disc_coeff) &
        &result(ierr)
        use num_vars, only: ndps => norm_disc_prec_sol, BC_style, EV_BC
        use X_vars, only: n_mod_X
        use MPI_utilities, only: get_ser_var, wait_MPI
        use num_utilities, only: con, c
        use SLEPC_utilities, only: insert_block_mat
        
        character(*), parameter :: rout_name = 'set_BC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        Mat, intent(inout) :: A, B                                              ! Matrices A and B from A X = lambda B X
        integer, intent(in) :: i_geo                                            ! at which geodesic index to perform the calculations
        integer, intent(in) :: n_sol                                            ! number of grid points of solution grid
        PetscReal, intent(in) :: norm_disc_coeff(:)                             ! discretization coefficients for normal derivatives
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscInt :: n_min, n_max                                                ! absolute limits excluding the BC's
        
        ! local variables
        PetscInt :: kd                                                          ! counter
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        call writo('Preparing variables')
        call lvl_ud(1)
        
        ! check nr. of modes
        if (n_mod_X.ne.X%n_mod(1) .or. n_mod_X.ne.X%n_mod(2)) then
            ierr = 1
            err_msg = 'Need square matrix of size [n_mod_X:n_mod_X]'
            CHCKERR(err_msg)
        end if
        
        call lvl_ud(-1)
        
        ! set up n_min, depending on BC style
        select case (BC_style(1))
            case (1)
                n_min = 2*ndps                                                  ! dirichlet BC requires half stencil
            case (2)
                n_min = ndps                                                    ! mixed BC requires only one fourth of stencil
            case (3)
                err_msg = 'Left BC''s cannot have BC type 3'
                ierr = 1
                CHCKERR(err_msg)
            case default
                err_msg = 'No BC style associated with '//&
                    &trim(i2str(BC_style(1)))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        ! set up n_max, depending on BC style
        select case (BC_style(2))
            case (1)
                n_max = 2*ndps                                                  ! dirichlet BC requires half stencil
            case (2)
                n_max = ndps                                                    ! mixed BC requires only one fourth of stencil
            case (3)
                n_max = 1                                                       ! only last element carries BC
            case default
                err_msg = 'No BC style associated with '//&
                    &trim(i2str(BC_style(1)))
                ierr = 1
                CHCKERR(err_msg)
        end select
        
        call writo('Setting BC deep within plasma')
        call lvl_ud(1)
        call writo('Using artificial eigenvalue EV_BC = '//trim(r2strt(EV_BC)))
        
        ! iterate over all positions where to set left BC
        do kd = 1,n_min
            if (grid_X%i_min.le.kd .and. grid_X%i_max.ge.kd) then               ! decide which process does the setting
                select case (BC_style(1))
                    case (1)
                        ierr = set_BC_1(kd-1,A,B,.false.)                       ! indices start at 0
                        CHCKERR('')
                    case (2)
                        ierr = set_BC_2(kd-1,kd-grid_X%i_min+1,X,A,B,&
                            &i_geo,grid_X%n(3),norm_disc_coeff,.false.)       ! indices start at 0
                        CHCKERR('')
                    case (3)
                        err_msg = 'Left BC''s cannot have BC type 3'
                        ierr = 1
                        CHCKERR(err_msg)
                    case default
                        err_msg = 'No BC style associated with '//&
                            &trim(i2str(BC_style(1)))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            end if
        end do
        
        ! assemble the matrices (cannot mix overwrite and not)
        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t begin assembly of matrix')
        call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t begin assembly of matrix')
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t end assembly of matrix')
        call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t end assembly of matrix')
        
        call lvl_ud(-1)
        
        ! wait to get messages consistent
        ierr = wait_MPI()
        CHCKERR('')
        
        call writo('Setting BC at plasma edge')
        call lvl_ud(1)
        
        ! iterate over all positions where to set right BC
        do kd = grid_X%n(3),grid_X%n(3)-n_max+1,-1
            if (grid_X%i_min.le.kd .and. grid_X%i_max.ge.kd) then               ! decide which process does the setting
                select case (BC_style(2))
                    case (1)
                        ierr = set_BC_1(kd-1,A,B,.true.)                        ! indices start at 0
                        CHCKERR('')
                    case (2)
                        ierr = set_BC_2(kd-1,kd-grid_X%i_min+1,X,A,B,i_geo,&
                            &grid_X%n(3),norm_disc_coeff,.true.)                ! indices start at 0
                        CHCKERR('')
                    case (3)
                        if (kd.eq.grid_X%n(3)) then                             ! only for last point, irrespective of discretization order
                            ierr = set_BC_3(kd-1,X,A)                           ! indices start at 0
                            CHCKERR('')
                        end if
                    case default
                        err_msg = 'No BC style associated with '//&
                            &trim(i2str(BC_style(1)))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
            end if
        end do
        
        ! assemble the matrices
        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t begin assembly of matrix')
        call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t begin assembly of matrix')
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t end assembly of matrix')
        call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
        CHCKERR('Coulnd''t end assembly of matrix')
        
        call lvl_ud(-1)
        
        call lvl_ud(-1)
    contains
        ! set BC style 1:
        !   A(ind,ind) = EV_BC, B(ind,ind) = 1,
        !   A(ind+1..ind+p,ind+1..ind+p) = 0, B(ind+1..ind+p,ind+1..ind+p) = 0,
        ! where ind indicates the row where the BC is centered.
        integer function set_BC_1(r_id,A,B,BC_right) result(ierr)
            character(*), parameter :: rout_name = 'set_BC_1'
            
            ! input / output
            integer, intent(in) :: r_id                                         ! position at which to set BC
            Mat, intent(inout) :: A, B                                          ! Matrices A and B from A X = lambda B X
            logical :: BC_right                                                 ! if BC is at right (so there are vacuum terms)
            
            ! local variables
            integer :: kd                                                       ! counter
            PetscScalar, allocatable :: loc_block(:,:)                          ! [n_mod_X:n_mod_X] block matrix for 1 normal point
            PetscInt :: pmone                                                   ! plus or minus one
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(r_id+1))//&
                &': Eigenvector set to zero',persistent=.true.)
            
            ! initialize local blocks
            allocate(loc_block(n_mod_X,n_mod_X))
            loc_block = 0.0_dp
            do kd = 1,n_mod_X
                loc_block(kd,kd) = 1.0_dp
            end do
            
            ! set up plus minus one
            if (BC_right) then
                pmone = 1
            else
                pmone = -1
            end if
            
            ! set block r_id + (0,0)
            ierr = insert_block_mat(EV_BC*loc_block,A,r_id,[0,0],n_sol,&
                &overwrite=.true.)
            CHCKERR('')
            ierr = insert_block_mat(loc_block,B,r_id,[0,0],n_sol,&
                &overwrite=.true.)
            CHCKERR('')
            
            ! iterate over range 2
            do kd = 1,2*ndps
                ! set block r_id +/- (0,kd) and Hermitian conjugate
                ierr = insert_block_mat(0*loc_block,A,r_id,[0,-pmone*kd],&
                    &n_sol,overwrite=.true.,transp=.true.)
                CHCKERR('')
                ierr = insert_block_mat(0*loc_block,B,r_id,[0,-pmone*kd],&
                    &n_sol,overwrite=.true.,transp=.true.)
                CHCKERR('')
            end do
            
            ! clean up
            deallocate(loc_block)
        end function set_BC_1
        
        ! set BC style 2:
        !   minimization of surface energy term (see [ADD REF])
        integer function set_BC_2(r_id,r_id_loc,X,A,B,i_geo,n_sol,&
            &norm_disc_coeff,BC_right)                                          result(ierr)
            character(*), parameter :: rout_name = 'set_BC_2'
            
            ! input / output
            integer, intent(in) :: r_id                                         ! global position at which to set BC (starting at 0)
            integer, intent(in) :: r_id_loc                                     ! local index in perturbation tables
            type(X_2_type), intent(in) :: X                                     ! field-averaged perturbation variables (so only first index)
            Mat, intent(inout) :: A, B                                          ! Matrices A and B from A X = lambda B X
            integer, intent(in) :: i_geo                                        ! at which geodesic index to perform the calculations
            integer, intent(in) :: n_sol                                        ! number of grid points of solution grid
            PetscReal, intent(in) :: norm_disc_coeff(:)                         ! discretization coefficients for normal derivatives
            logical :: BC_right                                                 ! if BC is at right (so there are vacuum terms)
            
            ! local variables
            PetscInt :: jd                                                      ! counter
            PetscScalar, allocatable :: V_int_0_mod(:,:,:)                      ! V_0 + (V_1+delta) V_2 (V_1+delta)^*T
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(r_id+1))//&
                &': Minimization of surface energy',persistent=.true.)
            
            ! -----------------!
            ! BLOCKS ~ V_0,mod !
            ! -----------------!
            ! calculate modified terms V_int_0_mod
            allocate(V_int_0_mod(n_mod_X,n_mod_X,2))
            ierr = calc_V_0_mod(X%PV_0(1,i_geo,r_id_loc,:),&
                &X%KV_0(1,i_geo,r_id_loc,:),&
                &X%PV_1(1,i_geo,r_id_loc,:),&
                &X%KV_1(1,i_geo,r_id_loc,:),&
                &X%KV_2(1,i_geo,r_id_loc,:),V_int_0_mod)
            CHCKERR('')
            
            ! add block to r_id + (0,0)
            ierr = insert_block_mat(V_int_0_mod(:,:,1),A,r_id,[0,0],n_sol)
            CHCKERR('')
            ierr = insert_block_mat(V_int_0_mod(:,:,2),B,r_id,[0,0],n_sol)
            CHCKERR('')
            
            ! deallocate modified V_0
            deallocate(V_int_0_mod)
            
            ! -------------!
            ! BLOCKS ~ vac !
            ! -------------!
            ! add block to r_id + (0,-p..p) + Hermitian conjugate
            if (BC_right) then
                do jd = -ndps,ndps
                    if (r_id.lt.n_sol .and. r_id+jd.lt.n_sol) then
                        ierr = insert_block_mat(-X%vac_res*&
                            &norm_disc_coeff(jd+ndps+1),A,&
                            &r_id,[0,jd],n_sol,transp=.true.)
                        CHCKERR('')
                    end if
                end do
            end if
        end function set_BC_2
        
        ! set BC style 3:
        !   minimization of vacuum energy (see [ADD REF])
        integer function set_BC_3(r_id,X,A) result(ierr)
            character(*), parameter :: rout_name = 'set_BC_3'
            
            ! input / output
            integer, intent(in) :: r_id                                         ! position at which to set BC
            type(X_2_type), intent(in) :: X                                     ! field-averaged perturbation variables (so only first index)
            Mat, intent(inout) :: A                                             ! Matrices A from A X = lambda B X
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(r_id+1))//&
                &': Minimization of vacuum energy',persistent=.true.)
            
            ! -------------!
            ! BLOCKS ~ vac !
            ! -------------!
            ! add block to r_id + (0,0)
            ierr = insert_block_mat(X%vac_res,A,r_id,[0,0],n_sol)
            CHCKERR('')
        end function set_BC_3
        
        ! calculates V_0 + (V_1+delta) V_2^-1 (V_1+delta)^*T
        integer function calc_V_0_mod(PV_0,KV_0,PV_1,KV_1,KV_2,V_0_mod) &
            &result(ierr)
            character(*), parameter :: rout_name = 'calc_V_0_mod'
            
            ! input / output
            complex(dp), intent(in) :: PV_0(:), KV_0(:)                         ! PV_0 and KV_0
            complex(dp), intent(in) :: PV_1(:), KV_1(:)                         ! PV_1 and KV_1
            complex(dp), intent(in) :: KV_2(:)                                  ! KV_2
            complex(dp), intent(inout) :: V_0_mod(:,:,:)                        ! output V_0_mod
            
            ! local variables
            complex(dp), allocatable :: KV_2_inv(:,:)                           ! inverse of KV_2
            complex(dp), allocatable :: KV_1_loc(:,:)                           ! local copy of KV_1
            complex(dp), allocatable :: PV_1_loc(:,:)                           ! local copy of PV_1
            complex(dp), allocatable :: V_triple(:,:)                           ! product of 3 V's
            character :: uplo                                                   ! upper or lower Cholesky factorization used to invert KV_2
            PetscInt :: k, m                                                    ! counters
#if ldebug
            integer :: istat                                                    ! status
#endif
            
            ! initialize ierr
            ierr = 0
            
            ! initialize variables
            allocate(KV_2_inv(n_mod_X,n_mod_X))
            allocate(V_triple(n_mod_X,n_mod_X))
            allocate(KV_1_loc(n_mod_X,n_mod_X))
            allocate(PV_1_loc(n_mod_X,n_mod_X))
            
            ! make local copy of KV_2 into the inverse array
            KV_2_inv = 0._dp
            do m = 1,n_mod_X
                do k = 1,m                                                      ! only save upper diagonal part
                    KV_2_inv(k,m) = con(KV_2(c([k,m],.true.,n_mod_X)),&
                        &[k,m],.true.)                                          ! symmetric matrices need con()
                end do
            end do
            
#if ldebug
            if (debug_calc_V_0_mod) then
                write(*,*,IOSTAT=istat) 'KV_2 = '
                
                write(*,*,IOSTAT=istat) 'real part:'
                call print_ar_2(rp(KV_2_inv))
                write(*,*,IOSTAT=istat) 'imaginary part:'
                call print_ar_2(ip(KV_2_inv))
            end if
#endif
            
            ! invert matrix KV_2 using Lapack
            uplo = 'U'                                                          ! only lower diagonal matters, but is arbitrary
            call zpotrf(uplo,n_mod_X,KV_2_inv,n_mod_X,ierr)
            CHCKERR('Failed to decompose KV_2')
            call zpotri(uplo,n_mod_X,KV_2_inv,n_mod_X,ierr)
            CHCKERR('Failed to invert KV_2')
            do k = 1,n_mod_X
                do m = 1,k-1
                    KV_2_inv(k,m) = conjg(KV_2_inv(m,k))
                end do
            end do
            
            ! save local copy of KV_1 and PV_1
            do m = 1,n_mod_X
                do k = 1,n_mod_X
                    KV_1_loc(k,m) = KV_1(c([k,m],.false.,n_mod_X))              ! asymetric matrices don't need con()
                    PV_1_loc(k,m) = PV_1(c([k,m],.false.,n_mod_X))              ! asymetric matrices don't need con()
                end do
            end do
            
            ! -------------!
            ! BLOCKS in PV !
            ! -------------!
            
            ! multiply PV_1 + delta^vac with inverse KV_2 and save in V_triple
            V_triple = matmul(PV_1_loc+X%vac_res,KV_2_inv)
            
            ! multiply with KV_1^*T and save in V_triple
            V_triple = matmul(V_triple,conjg(transpose(KV_1_loc)))
            
            ! fill local block
            do m = 1,n_mod_X
                do k = 1,n_mod_X
                    V_0_mod(k,m,1) = &
                        &con(PV_0(c([k,m],.true.,n_mod_X)),[k,m],.true.)        ! symmetric matrices need con()
                end do
            end do
            V_0_mod(:,:,1) = V_0_mod(:,:,1) + &
                &V_triple - conjg(transpose(V_triple))
            
            ! -------------!
            ! BLOCKS in KV !
            ! -------------!
            
            ! multiply KV_1 with inverse KV_2 and save in V_triple
            V_triple = matmul(KV_1_loc,KV_2_inv)
            
            ! multiply with KV_1^*T and save in V_triple
            V_triple = matmul(V_triple,conjg(transpose(KV_1_loc)))
            
            ! fill V_0_mod(1)
            do m = 1,n_mod_X
                do k = 1,n_mod_X
                    V_0_mod(k,m,2) = &
                        &con(KV_0(c([k,m],.true.,n_mod_X)),[k,m],.true.)        ! symmetric matrices need con()
                end do
            end do
            V_0_mod(:,:,2) = V_0_mod(:,:,2) - V_triple
#if ldebug
        if (debug_calc_V_0_mod) then
            write(*,*,IOSTAT=istat) 'KV_2^-1 = '
            
            write(*,*,IOSTAT=istat) 'real part:'
            call print_ar_2(rp(KV_2_inv))
            write(*,*,IOSTAT=istat) 'imaginary part:'
            call print_ar_2(ip(KV_2_inv))
            
            write(*,*,IOSTAT=istat) 'PV_1 = '
            
            write(*,*,IOSTAT=istat) 'real part:'
            call print_ar_2(rp(PV_1_loc))
            write(*,*,IOSTAT=istat) 'imaginary part:'
            call print_ar_2(ip(PV_1_loc))
            
            write(*,*,IOSTAT=istat) 'vac = '
            
            write(*,*,IOSTAT=istat) 'real part:'
            call print_ar_2(rp(X%vac_res))
            write(*,*,IOSTAT=istat) 'imaginary part:'
            call print_ar_2(ip(X%vac_res))
            
            write(*,*,IOSTAT=istat) 'KV_1 = '
            
            write(*,*,IOSTAT=istat) 'real part:'
            call print_ar_2(rp(KV_1_loc))
            write(*,*,IOSTAT=istat) 'imaginary part:'
            call print_ar_2(ip(KV_1_loc))
            
            write(*,*,IOSTAT=istat) 'PV_0_mod = '
            
            write(*,*,IOSTAT=istat) 'real part:'
            call print_ar_2(rp(V_0_mod(:,:,1)))
            write(*,*,IOSTAT=istat) 'imaginary part:'
            call print_ar_2(ip(V_0_mod(:,:,1)))
            
            write(*,*,IOSTAT=istat) 'KV_0_mod = '
            
            write(*,*,IOSTAT=istat) 'real part:'
            call print_ar_2(rp(V_0_mod(:,:,2)))
            write(*,*,IOSTAT=istat) 'imaginary part:'
            call print_ar_2(ip(V_0_mod(:,:,2)))
        end if
#endif
            
            ! clean up
            deallocate(KV_2_inv)
            deallocate(KV_1_loc)
            deallocate(PV_1_loc)
            deallocate(V_triple)
        end function calc_V_0_mod
    end function set_BC
    
    ! sets up EV solver
    integer function setup_solver(X,A,B,solver) result(ierr)
        use num_vars, only: n_sol_requested, tol_SLEPC, max_it_slepc
        use rich_vars, only: rich_lvl
        use X_vars, only: n_mod_X
        use grid_vars, only: n_r_sol
        
        character(*), parameter :: rout_name = 'setup_solver'
        
        ! input / output
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! local variables
        PetscInt :: n_sol                                                       ! how many solutions can be requested (normally n_sol_requested)
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! check nr. of modes
        if (n_mod_X.ne.X%n_mod(1) .or. n_mod_X.ne.X%n_mod(2)) then
            ierr = 1
            err_msg = 'Need square matrix of size [n_mod_X:n_mod_X]'
            CHCKERR(err_msg)
        end if
        
        !call PetscOptionsSetValue('-eps_view','-1',ierr)
        call EPSCreate(PETSC_COMM_WORLD,solver,ierr)
        CHCKERR('EPSCreate failed')
        
        call EPSSetOperators(solver,A,B,ierr)                                   ! generalized EV problem A X = lambda B X
        CHCKERR('EPSetOperators failed')
        
        if (use_hermitian) then
            call EPSSetProblemType(solver,EPS_GHEP,ierr)
            CHCKERR('Set problem type failed')
        end if
        
        call EPSSetType(solver,EPSGD,ierr)
        err_msg = 'Failed to set type to generalized davidson'
        CHCKERR(err_msg)
        
        ! search for  Eigenvalue with smallest real value (so  imaginary part of
        ! sqrt(EV) = omega is largest)
        call EPSSetWhichEigenpairs(solver,EPS_SMALLEST_REAL,ierr)
        CHCKERR('Failed to set which eigenpairs')
        
        ! request n_sol_requested Eigenpairs
        if (n_sol_requested.gt.n_r_sol*n_mod_X) then
            call writo('max. nr. of solutions requested capped to problem &
                &dimension ('//trim(i2str(n_r_sol*n_mod_X))//')',&
                &warning=.true.)
            call writo('Increase either n_r_sol or number of pol. modes or &
                &decrease n_sol_requested')
            n_sol = n_r_sol*n_mod_X
        else
            n_sol = n_sol_requested
        end if
        call EPSSetDimensions(solver,n_sol,PETSC_DECIDE,&
            &PETSC_DECIDE,ierr)
        CHCKERR('EPSSetDimensions failed')
        call EPSSetTolerances(solver,tol_SLEPC(rich_lvl),max_it_slepc,ierr)
        CHCKERR('EPSSetTolerances failed')
        
        ! user output
        call writo('tolerance: '//trim(r2str(tol_SLEPC(rich_lvl))))
        call writo('maximum nr. of iterations: '//trim(i2str(max_it_slepc)))
        
        call lvl_ud(-1)
    end function setup_solver
    
    ! sets up guess in solver
    integer function setup_guess(sol,A,solver) result(ierr)
        character(*), parameter :: rout_name = 'setup_guess'
        
        ! input / output
        type(sol_type), intent(in) :: sol                                       ! solution variables
        Mat, intent(in) :: A                                                    ! matrix A (or B)
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! local variables
        Vec, allocatable :: guess_vec(:)                                        ! guess to solution EV parallel vector
        PetscInt :: kd                                                          ! counter
        PetscInt :: n_EV_prev                                                   ! nr. of previous solutions
        PetscScalar, pointer ::  guess_vec_ptr(:)                               ! pointer to local guess_vec
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set guess for EV if sol vec is allocated
        if (allocated(sol%val)) then
            ! set n_EV_prev
            n_EV_prev = size(sol%val)
            
            ! allocate guess vectors
            allocate(guess_vec(n_EV_prev))
            
            ! create vecctor guess_vec and set values
            do kd = 1,n_EV_prev
                ! create the vectors
                !call MatGetVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,ierr)         ! get compatible parallel vectors to matrix A (petsc 3.5.3)
                call MatCreateVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,ierr)      ! get compatible parallel vectors to matrix A (petsc 3.6.1)
                CHCKERR('Failed to create vector')
                
                ! get pointer
                call VecGetArrayF90(guess_vec(kd),guess_vec_ptr,ierr)
                CHCKERR('Failed to get pointer')
                
                ! copy the values
                guess_vec_ptr = reshape([sol%vec(:,:,kd)],&
                    &[size(sol%vec(:,:,kd))])
                
                ! return pointer
                call VecRestoreArrayF90(guess_vec(kd),guess_vec_ptr,ierr)
                CHCKERR('Failed to restore pointer')
                
                !! visualize guess
                !call VecView(guess_vec(kd),PETSC_VIEWER_STDOUT_WORLD,ierr)
                !CHCKERR('Cannot view vector')
            end do
                
            ! set guess
            call EPSSetInitialSpace(solver,n_EV_prev,guess_vec,ierr)
            CHCKERR('Failed to set guess')
            
            ! destroy guess vector
            call VecDestroy(guess_vec,ierr)                                     ! destroy guess vector
            CHCKERR('Failed to destroy vector')
        end if
        
        call lvl_ud(-1)
    end function setup_guess
    
    ! get the solution vectors
    integer function get_solution(solver) result(ierr)
        character(*), parameter :: rout_name = 'get_solution'
        
        ! input / output
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set run-time options
        call EPSSetFromOptions(solver,ierr)
        CHCKERR('EPSetFromOptions failed')
        
        ! solve EV problem
        call EPSSolve(solver,ierr) 
        err_msg = 'EPS couldn''t find a solution. Maybe you should increase &
            &the number of parallel points.'
        CHCKERR(err_msg)
        
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
        call EPSGetTolerances(solver,tol,max_it,ierr)                           ! tolerance and max. nr. of iterations
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
        call writo('number of iterations: '//trim(i2str(n_it)))
        call writo('number of converged solutions: '//trim(i2str(n_conv)))
        call writo('number of requested solutions: '//trim(i2str(n_ev)))
        call writo('maximum dimension of the subspace to be used by solver: '//&
            &trim(i2str(ncv)))
        call writo('maximum dimension allowed for projected problem : '//&
            &trim(i2str(mpd)))
        
        ! set maximum nr of solutions to be saved
        if (n_sol_requested.gt.n_conv) then
            if (n_conv.eq.0) then
                call writo('no solutions found',warning=.true.)
            else
                call writo('max. nr. of solutions found only '//&
                    &trim(i2str(n_conv)),warning=.true.)
            end if
            max_n_EV = n_conv
        else
            max_n_EV = n_sol_requested
        end if
        
        ! check whether solutions found
        if (max_n_EV.le.0) then
            ierr = 1
            CHCKERR('No solutions found')
        end if
        
        call lvl_ud(-1)
    end function summarize_solution
    
    ! stores the results
    integer function store_results(grid_sol,sol,solver,max_n_EV,A,B,&
        &step_size) result(ierr)
        
        use eq_vars, only: T_0
        use X_vars, only: n_mod_X
        use num_vars, only: use_normalization, EV_BC, prog_name, output_name, &
            &n_procs, rank, tol_SLEPC, eq_style, output_EV_i
        use MPI_utilities, only: get_ser_var
        use rich_vars, only: rich_lvl
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'store_results'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(inout) :: max_n_EV                                     ! nr. of EV's saved, up to n_conv
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        PetscReal, intent(in) :: step_size                                      ! step size in flux coordinates
        
        ! local variables
        type(grid_type) :: grid_sol_trim                                        ! trimmed solution grid
        Mat :: err_mat                                                          ! A - omega^2 B
        Vec :: sol_vec                                                          ! solution EV parallel vector
        Vec :: err_vec                                                          ! AX - omega^2 BX
        Vec :: E_vec(2)                                                         ! AX and BX
        PetscInt :: id                                                          ! counters
        PetscInt :: id_tot                                                      ! counters
        PetscInt :: one = 1                                                     ! one
        PetscInt, allocatable :: sol_vec_max_loc(:)                             ! location of sol_vec_max of this process
        PetscReal :: error                                                      ! error of EPS solver
        PetscReal :: error_est                                                  ! error estimate of EPS solver
        PetscReal :: err_norm                                                   ! norm of error
        PetscReal :: tol_complex                                                ! tolerance for an EV to be considered complex
        PetscReal :: tol_EV_BC = 1.E-3_dp                                       ! tolerance for an EV to be considered due to the BC
        PetscReal, parameter :: infinity = 1.E40_dp                             ! beyond this value, modes are not saved
        PetscScalar :: sol_val_loc                                              ! local X val
        PetscScalar :: err_val                                                  ! X*(A-omega^2B)X
        PetscScalar :: E_val(2)                                                 ! X*AX and X*BX
        PetscScalar, allocatable :: sol_vec_max(:)                              ! max of sol_vec of all processes, including phase at max
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: format_val                                 ! format
        character(len=max_str_ln) :: format_head                                ! header
        character(len=max_str_ln) :: EV_err_str                                 ! String with information about EV error
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=2*max_str_ln) :: temp_output_str                          ! temporary output string
        integer :: n_digits                                                     ! nr. of digits for the integer number
        integer :: n_err(3)                                                     ! how many errors there were
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid_sol,grid_sol_trim)
        CHCKERR('')
        
        ! create solution variables
        if (allocated(sol%val)) call sol%dealloc()
        call sol%init(grid_sol_trim,max_n_EV)
        
        ! create solution vector
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,loc_n_r*n_mod_X,&
            &n_r*n_mod_X,PETSC_NULL_SCALAR,sol_vec,ierr)
        CHCKERR('Failed to create MPI vector with arrays')
        
        ! set up EV error string and format string:
        !   1: index of EV
        !   2: Real part
        !   3: Imaginary part
        !   4: relative error ||Ax - EV Bx||_2/||EV x||_2 [1]
        n_err = 0
        EV_err_str = ''
        if (rank.eq.0) then
            n_digits = min(ceiling(log10(1._dp*max_n_EV)),16)
            format_val = '(A2,I'//trim(i2str(n_digits))//'," ",ES23.16," ",&
                &ES23.16," ",ES23.16)'
            format_head = '(A'//trim(i2str(n_digits+3))//'," ",A23," ",A23," ",&
                &A23)'
        end if
        
        ! open output file for the log
        if (rank.eq.0) then
            full_output_name = prog_name//'_'//trim(output_name)//'_EV_R_'//&
                &trim(i2str(rich_lvl))//'.txt'
            open(output_EV_i,FILE=full_output_name,STATUS='replace',&
                &IOSTAT=ierr)
            CHCKERR('Cannot open EV output file')
        end if
        
        ! output to file
        if (rank.eq.0) then
            if (use_normalization) then
                write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) &
                    &'# Eigenvalues normalized to the squared Alfven &
                    &frequency omega_A^2 = '
                CHCKERR('Failed to write')
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) &
                            &'#     ('//trim(r2str(1._dp/T_0))//' Hz)^2 = '//&
                            &trim(r2str(1._dp/T_0**2))//' Hz^2'
                        CHCKERR('Failed to write')
                    case (2)                                                    ! HELENA
                        write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) &
                            &'#     (HELENA normalization)'
                        CHCKERR('Failed to write')
                end select
            else
                write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) '# Eigenvalues'
                CHCKERR('Failed to write')
            end if
            write(temp_output_str,format_head) &
                &'#  I                            ', &
                &'real part               ', 'imaginary part          ', &
                &'relative precision      '
            write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) trim(temp_output_str)
            CHCKERR('Failed to write')
        end if
        
        ! initialize helper variables
        allocate(sol_vec_max(n_procs))
        allocate(sol_vec_max_loc(2))
        
        ! start id
        id = 1
        id_tot = 1
        
        ! store them
        do while (id.le.max_n_EV)
            ! get EV solution in vector X vec
            call VecPlaceArray(sol_vec,sol%vec(:,:,id),ierr)                    ! place array sol_vec in solution vec
            CHCKERR('Failed to place array')
            call EPSGetEigenpair(solver,id_tot-1,sol%val(id),PETSC_NULL_OBJECT,&
                &sol_vec,PETSC_NULL_OBJECT,ierr)                                ! get solution EV vector and value (starts at index 0)
            CHCKERR('EPSGetEigenpair failed')
            call EPSComputeError(solver,id_tot-1,EPS_ERROR_RELATIVE,error,ierr) ! get error (starts at index 0) (petsc 3.6.1)
            !call EPSComputeRelativeError(solver,id_tot-1,error,ierr)            ! get error (starts at index 0) (petsc 3.5.3)
            CHCKERR('EPSComputeError failed')
            
            ! set up local solution val
            sol_val_loc = sol%val(id)
            
            ! tests
            tol_complex = tol_SLEPC(rich_lvl)
            if (abs(ip(sol%val(id))/rp(sol%val(id))).gt.tol_complex) then       ! test for unphysical complex solution
                EV_err_str = '# WARNING: Unphysical complex Eigenvalue!'
                n_err(1) = n_err(1)+1
                call remove_EV(id,max_n_EV)
            else if (abs(rp(sol%val(id)-EV_BC)/EV_BC).lt.tol_EV_BC) then        ! test for artificial EV due to BC
                EV_err_str = '# WARNING: Eigenvalue probably due to BC''s &
                    &(artifically set to '//trim(r2strt(EV_BC))//')'
                n_err(2) = n_err(2)+1
                call remove_EV(id,max_n_EV)
            else if (abs(sol%val(id)).gt.infinity) then                         ! test for infinity
                EV_err_str = '# WARNING: The next Eigenvalues are larger &
                    &than '//trim(r2strt(infinity))//' and are discarded'
                n_err(3) = 1
                call remove_EV(id,max_n_EV,remove_next=.true.)
                exit
            end if
            
            if (rank.eq.0) then
                ! output EV to file
                if (EV_err_str.ne.'') then
                    write(temp_output_str,format_val) '#', id,rp(sol_val_loc),&
                        &ip(sol_val_loc),error                                  ! if error, comment the next line
                else
                    write(temp_output_str,format_val) '', id,rp(sol_val_loc),&
                        &ip(sol_val_loc),error
                end if
                write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) &
                    &trim(temp_output_str)
                CHCKERR('Failed to write')
                
                ! if error, print explanation
                if (EV_err_str.ne.'') then
                    write(UNIT=output_EV_i,FMT='(A)',IOSTAT=ierr) &
                        &trim(EV_err_str)
                    CHCKERR('Failed to write')
                end if
            end if
            
            ! normalize Eigenvectors to make output more easily comparable
            if (EV_err_str.eq.'') then
                ! find local maximum
                sol_vec_max_loc = maxloc(abs(sol%vec(:,:,id)))
                ierr = get_ser_var(&
                    &[sol%vec(sol_vec_max_loc(1),sol_vec_max_loc(2),id)],&
                    &sol_vec_max,scatter=.true.)
                CHCKERR('')
                ! find global maximum
                sol_vec_max_loc(1) = maxloc(abs(sol_vec_max),1)
                sol%vec(:,:,id) = sol%vec(:,:,id) / &
                    &sol_vec_max(sol_vec_max_loc(1))
                
                call EPSGetErrorEstimate(solver,id_tot-1,error_est,ierr)        ! get error estimate
                CHCKERR('EPSGetErrorEstimate failed')
                ! user message
                call writo('Checking whether A x - omega^2 B x = 0 for EV '//&
                    &trim(i2str(id))//': '//trim(c2strt(sol%val(id))))
                call lvl_ud(1)
                
                ! set up error matrix A - omega^2 B
                call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,err_mat,ierr)
                err_msg = 'failed to duplicate mat into err_mat'
                CHCKERR(err_msg)
                call MatCopy(A,err_mat,MAT_SHARE_NONZERO_PATTERN,ierr)          ! err_mat has same structure as A
                CHCKERR('Failed to copy mat into mat_loc')
                call MatAXPY(err_mat,-sol%val(id),B,DIFFERENT_NONZERO_PATTERN,&
                    &ierr)                                                      ! for some reason, SAME_NONZERO_PATTERN does not work
                CHCKERR('Failed to perform AXPY')
                
                ! set up error vector
                call VecDuplicate(sol_vec,err_vec,ierr)
                CHCKERR('Failed to duplicate vector')
                
                ! calculate Ax-lambdaBx
                call MatMult(err_mat,sol_vec,err_vec,ierr)
                CHCKERR('Failed to multiply')
                
                ! multiply with X* and give output
                call VecDot(err_vec,sol_vec,err_val,ierr)
                CHCKERR('Failed to do dot product')
                call writo('X*(A-omega^2B)X = '//trim(c2strt(err_val)))
                
                ! calculate absolute norm
                call VecNorm(err_vec,NORM_2,err_norm,ierr)
                CHCKERR('Failed to calculate norm')
                
                ! get relative norm
                err_norm = err_norm/abs(sol%val(id))
                
                ! visualize solution and error
                call writo('error: '//trim(r2str(err_norm))//', given: '//&
                    &trim(r2str(error))//', estimate: '//trim(r2str(error_est)))
                !!call VecView(err_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
                !!CHCKERR('Cannot view vector')
                
                ! set up Energy vector
                call VecDuplicate(sol_vec,E_vec(1),ierr)
                CHCKERR('Failed to duplicate vector')
                call VecDuplicate(sol_vec,E_vec(2),ierr)
                CHCKERR('Failed to duplicate vector')
                
                ! calculate Ax and Bx
                call MatMult(A,sol_vec,E_vec(1),ierr)
                CHCKERR('Failed to multiply')
                call MatMult(B,sol_vec,E_vec(2),ierr)
                CHCKERR('Failed to multiply')
                
                ! multiply with X* and give output
                call VecDot(E_vec(1),sol_vec,E_val(1),ierr)
                CHCKERR('Failed to do dot product')
                call writo('step_size = '//trim(r2str(step_size)))
                call writo('E_pot = X*AX*step_size = '//&
                    &trim(c2str(E_val(1)*step_size)))
                call VecDot(E_vec(2),sol_vec,E_val(2),ierr)
                CHCKERR('Failed to do dot product')
                call writo('E_kin = X*BX*step_size = '//&
                    &trim(c2str(E_val(2)*step_size)))
                call writo('X*AX/X*BX = '//trim(c2str(E_val(1)/E_val(2))))
                call writo('omega^2   = '//trim(c2str(sol%val(id))))
                
                ! if stable EV, give warning
                if (rp(sol%val(id)).gt.0._dp) then
                    call writo('The eigenvalue is positive, so the mode could &
                        &be stable.',warning=.true.)
                    call lvl_ud(1)
                    call writo('But this could also be due to value of "ncv" &
                        &that is too low.')
                    call writo('Have a look at the POST output and/or try &
                        &increasing "ncv".')
                    call lvl_ud(-1)
                end if
                
                ! increment counter
                id = id + 1
                
                call lvl_ud(-1)
                
                ! clean up
                call VecDestroy(err_vec,ierr)                                   ! destroy error vector
                CHCKERR('Failed to destroy err_vec')
                call MatDestroy(err_mat,ierr)                                   ! destroy error matrix
                CHCKERR('Failed to destroy err_mat')
                call VecDestroy(E_vec,ierr)                                     ! destroy energy vector
                CHCKERR('Failed to destroy E_vec')
            else
                ! reinitialize error string if error
                EV_err_str = ''
            end if
            
            ! reset vector
            call VecResetArray(sol_vec,ierr)
            CHCKERR('Failed to reset array')
            
            ! increment total id
            id_tot = id_tot+1
        end do
        
        ! deallocate helper variables
        deallocate(sol_vec_max_loc)
        deallocate(sol_vec_max)
        
        ! close output file if master
        if (rank.eq.0) close(output_EV_i)
        
        ! user output
        if (rank.eq.0) call writo(trim(i2str(max_n_EV))//&
            &' Eigenvalues were written in the file "'//&
            &trim(full_output_name)//'"')
        if (sum(n_err).ne.0) then
            call writo('there were some Eigenvalues that were removed:')
            call lvl_ud(1)
            
            if (n_err(1).gt.0) call writo(trim(i2str(n_err(1)))//&
                &' unphysical complex Eigenvalues')
            if (n_err(2).gt.0) call writo(trim(i2str(n_err(2)))//&
                &' Eigenvalues due to the boundary conditions')
            if (n_err(3).gt.0) call writo('The final Eigenvalues became &
                &infinite')
            
            call lvl_ud(-1)
            call writo('(Override this behavior using "retain_all_sol")')
        end if
        
        if (max_n_EV.gt.0) then
            call writo('basic statistics:')
            call lvl_ud(1)
            
            call writo('min: '//trim(c2strt(sol%val(1))))
            call writo('max: '//trim(c2strt(sol%val(max_n_EV))))
            
            call lvl_ud(-1)
        end if
        
        !call EPSPrintSolution(solver,PETSC_NULL_OBJECT,ierr)
        
        call VecDestroy(sol_vec,ierr)                                           ! destroy solution vector
        CHCKERR('Failed to destroy sol_vec')
        
        call lvl_ud(-1)
        
        ! clean up
        call grid_sol_trim%dealloc()
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
            PetscScalar, allocatable :: sol_val_loc(:)                          ! local copy of sol_val
            PetscScalar, allocatable :: sol_vec_loc(:,:,:)                      ! local copy of sol_vec
            PetscBool :: remove_next_loc = .false.                              ! local copy of remove_next
            
            ! only remove if faulty solutions are not optionally retained
            if (retain_all_sol) then
                id = id+1                                                       ! increment the solution
            else
                ! set up local remove_next
                if (present(remove_next)) remove_next_loc = remove_next
                
                ! save old arrays
                allocate(sol_val_loc(max_id))
                allocate(sol_vec_loc(n_mod_X,grid_sol_trim%loc_n_r,max_id))
                sol_val_loc = sol%val
                sol_vec_loc = sol%vec
                
                ! reallocate arrays
                deallocate(sol%val,sol%vec)
                if (remove_next_loc) then
                    allocate(sol%val(id-1))
                    allocate(sol%vec(n_mod_X,grid_sol_trim%loc_n_r,id-1))
                else
                    allocate(sol%val(max_id-1))
                    allocate(sol%vec(n_mod_X,grid_sol_trim%loc_n_r,max_id-1))
                end if
                
                ! copy other values
                sol%val(1:id-1) = sol_val_loc(1:id-1)
                sol%vec(:,:,1:id-1) = sol_vec_loc(:,:,1:id-1)
                if (.not.remove_next_loc) then
                    sol%val(id:max_id-1) = sol_val_loc(id+1:max_id)
                    sol%vec(:,:,id:max_id-1) = sol_vec_loc(:,:,id+1:max_id)
                end if
                
                ! adapt max_id
                if (remove_next_loc) then
                    max_id = id-1
                else
                    max_id = max_id-1
                end if
                
                ! clean up
                deallocate(sol_val_loc,sol_vec_loc)
            end if
        end subroutine remove_EV
    end function store_results
    
    ! stop PETSC and SLEPC
    ! [MPI] Collective call
    integer function stop_SLEPC(A,B,solver) result(ierr)
        character(*), parameter :: rout_name = 'stop_SLEPC'
        
        ! input / output
        Mat, intent(in) :: A, B                                                 ! matrices A and B in EV problem A X = lambda B X
        EPS, intent(in) :: solver                                               ! EV solver
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
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
