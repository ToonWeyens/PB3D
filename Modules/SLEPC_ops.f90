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
! for slepc 3.6.0:
#include <slepc/finclude/slepcepsdef.h>
! for slepc 3.5.3:
!#include <finclude/slepcepsdef.h>
    use str_ops
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
    integer function solve_EV_system_SLEPC(grid_sol,X,sol,sol_prev,i_geo) &
        &result(ierr)
        use num_vars, only: max_it_inv, norm_disc_prec_sol, matrix_SLEPC_style
        use num_utilities, only: calc_coeff_fin_diff
        use rich_vars, only: use_guess
#if ldebug
        use num_vars, only: ltest
        use input_utilities, only: get_real, get_log
#endif
        
        character(*), parameter :: rout_name = 'solve_EV_system_SLEPC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        type(sol_type), intent(in), optional :: sol_prev                        ! previous solution variables
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
        ierr = calc_coeff_fin_diff(1,norm_disc_prec_sol,norm_disc_coeff)
        CHCKERR('')
        norm_disc_coeff = norm_disc_coeff/step_size                             ! scale by step size
        call lvl_ud(-1)
        
        ! set up the matrix
        call writo('Set up matrices')
        call lvl_ud(1)
        ierr = setup_mats(grid_sol,X,A,B,i_geo_loc,norm_disc_coeff)
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
            
            select case (matrix_SLEPC_style)
                case (1)                                                        ! sparse
                    ierr = set_BC(grid_sol,X,A,B,i_geo_loc,grid_sol%n(3),&
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
                    write(*,*) '!!!! DISABLED !!!!'
                    !!ierr = 1
                    !!err_msg = 'NOT YET IMPLEMENTED FOR SHELL MATRICES'
                    !!CHCKERR(err_msg)
            end select
            
            ! set up solver
            call writo('Set up EV solver')
#if ldebug
            if (ltest) then
                call writo('Test spectrum of A or B instead of solving &
                    &generalized Eigenvalue problem?')
                if (get_log(.false.)) then
                    call writo('Spectrum of A (true) or B (false)?')
                    if (get_log(.true.)) then                                   ! A
                        ierr = setup_solver(grid_sol,X,A,&
                            &PETSC_NULL_OBJECT,solver)
                        CHCKERR('')
                    else                                                        ! B
                        ierr = setup_solver(grid_sol,X,B,&
                            &PETSC_NULL_OBJECT,solver)
                        CHCKERR('')
                    end if
                else
                    ierr = setup_solver(grid_sol,X,A,B,solver)
                    CHCKERR('')
                end if
            else
#endif
                ierr = setup_solver(grid_sol,X,A,B,solver)
                CHCKERR('')
#if ldebug
            end if
#endif
            
            ! set up guess
            if (present(sol_prev) .and. use_guess) then
                call writo('Set up guess')
                call lvl_ud(1)
                ierr = setup_guess(sol_prev,A,solver)
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
        
#if ldebug
        ierr = store_results(grid_sol,X,sol,solver,max_n_EV,A,B,step_size)
        CHCKERR('')
#else
        ierr = store_results(grid_sol,X,sol,solver,max_n_EV)
        CHCKERR('')
#endif
        
        ! finalize
        call writo('Finalize SLEPC')
        
        ! stop SLEPC
        ierr = stop_SLEPC(grid_sol,A,B,solver)
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
            !Mat :: mat_loc                                                      ! local version of mat
            PetscScalar :: one = 1.0                                            ! one
            PetscViewer :: file_viewer                                          ! viewer to write matrix to file
            character(len=max_str_ln) :: file_name                              ! file name
            
            ! initialize ierr
            ierr = 0
            
            ! test if mat is hermitian
            call MatHermitianTranspose(mat,MAT_INITIAL_MATRIX,mat_t,ierr)
            CHCKERR('Hermitian transpose of mat failed')
            call MatAXPY(mat_t,-one,mat,SAME_NONZERO_PATTERN,ierr)
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
        use num_vars, only: n_procs
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
    ! The  geodesical index  at  which to  perform the  calculations  has to  be
    ! provided as well.
    ! Note: For normal usage, i_geo should be 1, or not present.
    integer function setup_mats(grid_sol,X,A,B,i_geo,norm_disc_coeff) &
        &result(ierr)
        use num_vars, only: norm_disc_prec_sol, matrix_SLEPC_style
        use SLEPC_utilities, only: insert_block_mat
#if ldebug
        use num_vars, only: ltest
        use input_utilities, only: get_real, get_log
#endif
        !!! TEMPORARY !!!
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'setup_mats'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        integer, intent(in) :: i_geo                                            ! at which geodesic index to perform the calculations
        PetscReal, intent(in) :: norm_disc_coeff(:)                             ! discretization coefficients for normal derivatives
        
        ! local variables
        integer :: n_mod                                                        ! nr. of modes
        PetscInt :: kd                                                          ! counter
        PetscInt :: n_r                                                         ! n_r of trimmed sol grid
        PetscInt :: loc_n_r                                                     ! loc_n_r of trimmed sol grid
        PetscInt :: st_size                                                     ! half stencil size
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        !!! TEMPORARY !!!
        Vec :: vec_in
        Vec :: vec_out
        PetscScalar, pointer :: vec_in_loc(:)
        
        ! local variables also used in child routines
        PetscInt :: bulk_i_lim(2)                                               ! absolute limits of bulk matrix (excluding the BC's)
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Normal discretization with central finite differences of &
            &order '//trim(i2str(norm_disc_prec_sol))//', stencil width '//&
            &trim(i2str(4*norm_disc_prec_sol+1)))
        
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
        
        ! set nr. of modes
        n_mod = X%n_mod(1)
        if (n_mod.ne.X%n_mod(2)) then
            ierr = 1
            CHCKERR('Need square matrix')
        end if
        
        ! set up loc_n_r and n_r
        loc_n_r = grid_sol%loc_n_r
        n_r = grid_sol%n(3)
        
        ! set (half) stencil size
        st_size = 2*norm_disc_prec_sol
        
        ! set up bulk matrix absolute limits
        ierr = set_bulk_lims(grid_sol,bulk_i_lim)
        CHCKERR('')
        
        ! setup matrix A and B
        select case (matrix_SLEPC_style)
            case (1)                                                            ! sparse
                ! calculate nonzeros
                call set_nonzeros()
                
                ! create matrix A
                call MatCreateAIJ(PETSC_COMM_WORLD,loc_n_r*n_mod,loc_n_r*n_mod,&
                    &n_r*n_mod,n_r*n_mod,PETSC_NULL_INTEGER,d_nz*n_mod,&
                    &PETSC_NULL_INTEGER,o_nz*n_mod,A,ierr)
                err_msg = 'MatCreateAIJ failed for matrix A'
                CHCKERR(err_msg)
                
                ! deallocate tot_nz, d_nz and o_nz
                deallocate(d_nz,o_nz)
                
                ! fill the matrix A
                ierr = fill_mat(X%PV_0(1,i_geo,:,:),X%PV_1(1,i_geo,:,:),&
                    &X%PV_2(1,i_geo,:,:),n_r,norm_disc_coeff,bulk_i_lim,A)
                CHCKERR('')
                
                call writo('matrix A set up:')
                ierr = disp_mat_info(A)
                CHCKERR('')
                
                ! duplicate A into B
                ! (the  advantage  of   letting  communication  and  calculation
                ! overlap by calculating  matrix B while assembling  A is offset
                ! by the extra cost of not  being able to use MatDuplicate. This
                ! is also easier)
                call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)           ! B has same structure as A
                CHCKERR('failed to duplicate A into B')
                
                ! fill the matrix B
                ierr = fill_mat(X%KV_0(1,i_geo,:,:),X%KV_1(1,i_geo,:,:),&
                    &X%KV_2(1,i_geo,:,:),n_r,norm_disc_coeff,bulk_i_lim,B)
                CHCKERR('')
                
                call writo('matrix B set up:')
                ierr = disp_mat_info(B)
                CHCKERR('')
            case (2)                                                            ! shell
                ! create matrix A
                call MatCreateShell(PETSC_COMM_WORLD,loc_n_r*n_mod,&
                    &loc_n_r*n_mod,n_r*n_mod,n_r*n_mod,PETSC_NULL_OBJECT,A,ierr)
                err_msg = 'MatCreateShell failed for matrix A'
                CHCKERR(err_msg)
                
                ! set operations for A
                err_msg = 'MatShellSetOperation failed for matrix A'
                call MatShellSetOperation(A,MATOP_MULT,shell_A_product,ierr)
                CHCKERR(err_msg)
                call MatShellSetOperation(A,MATOP_GET_DIAGONAL,&
                    &shell_A_diagonal,ierr)
                CHCKERR(err_msg)
                
                !!!!!!!! TEMPORARY !!!!!!!!!!
                ! test multiplication of A
                call MatCreateVecs(A,vec_in,PETSC_NULL_OBJECT,ierr)        ! get compatible parallel vectors to matrix A (petsc 3.6.1)
                CHCKERR('Failed to create vector')
                call MatCreateVecs(A,vec_out,PETSC_NULL_OBJECT,ierr)        ! get compatible parallel vectors to matrix A (petsc 3.6.1)
                CHCKERR('Failed to create vector')
                call VecGetArrayF90(vec_in,vec_in_loc,ierr)
                CHCKERR('Failed to get pointer')
                vec_in_loc = [(kd*1._dp,kd=1,size(vec_in_loc))]
                call MatMult(A,vec_in,vec_out,ierr)
                CHCKERR('Failed to multiply')
                write(*,*) 'vec_in = '
                call VecView(vec_in,PETSC_VIEWER_STDOUT_WORLD,ierr)
                CHCKERR('Cannot view vector')
                write(*,*) 'vec_out = '
                call VecView(vec_out,PETSC_VIEWER_STDOUT_WORLD,ierr)
                CHCKERR('Cannot view vector')
                call VecRestoreArrayReadF90(vec_in,vec_in_loc,ierr)
                CHCKERR('Failed to restore pointer')
                ierr = wait_MPI()
                CHCKERR('')
                !!!!!!!! TEMPORARY !!!!!!!!!!
                
                ! create matrix B
                call MatCreateShell(PETSC_COMM_WORLD,loc_n_r*n_mod,&
                    &loc_n_r*n_mod,n_r*n_mod,n_r*n_mod,PETSC_NULL_OBJECT,B,ierr)
                err_msg = 'MatCreateShell failed for matrix B'
                CHCKERR(err_msg)
                
                ! set operations for B
                err_msg = 'MatShellSetOperation failed for matrix B'
                call MatShellSetOperation(B,MATOP_MULT,shell_B_product,ierr)
                CHCKERR(err_msg)
                call MatShellSetOperation(B,MATOP_GET_DIAGONAL,&
                    &shell_B_diagonal,ierr)
                CHCKERR(err_msg)
        end select
    contains
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
        integer function fill_mat(V_0,V_1,V_2,n_sol,norm_disc_coeff,bulk_i_lim,&
            &mat) result(ierr)
            use num_utilities, only: c, con
            
            character(*), parameter :: rout_name = 'fill_mat'
            
            ! input / output
            PetscScalar, intent(in) :: V_0(:,:)                                 ! either PV_int_0 or KV_int_0 in equilibrium normal grid
            PetscScalar, intent(in) :: V_1(:,:)                                 ! either PV_int_1 or KV_int_1 in equilibrium normal grid
            PetscScalar, intent(in) :: V_2(:,:)                                 ! either PV_int_2 or KV_int_2 in equilibrium normal grid
            integer, intent(in) :: n_sol                                        ! number of grid points of solution grid
            PetscReal, intent(in) :: norm_disc_coeff(:)                         ! discretization coefficients for normal derivatives
            PetscInt, intent(in) :: bulk_i_lim(2)                               ! absolute limits of bulk matrix (excluding the BC's)
            Mat, intent(inout) :: mat                                           ! either A or B
            
            ! local variables (not to be used in child routines)
            character(len=max_str_ln) :: err_msg                                ! error message
            PetscScalar, allocatable :: loc_block(:,:)                          ! (n_mod x n_mod) block matrix for 1 normal point
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
            r_sol_start = r_sol_start/n_mod                                     ! count per block
            r_sol_end = r_sol_end/n_mod                                         ! count per block
            if (grid_sol%i_min.ne.r_sol_start+1) then
                ierr = 1
                err_msg = 'start of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            if (grid_sol%i_max.ne.r_sol_end) then
                ierr = 1
                err_msg = 'end of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            
            ! allocate local block
            allocate(loc_block(n_mod,n_mod))
            
            ! iterate over all rows of this rank
            do kd = bulk_i_lim(1)-1,bulk_i_lim(2)-1                             ! (indices start with 0 here)
                ! set up local kd
                kd_loc = kd+2-grid_sol%i_min
                
                ! fill the blocks
                
                ! -------------!
                ! BLOCKS ~ V_0 !
                ! -------------!
                ! fill local block
                do m = 1,n_mod
                    do k = 1,n_mod
                        loc_block(k,m) = con(V_0(kd_loc,c([k,m],.true.,n_mod)),&
                            &[k,m],.true.)                                      ! symmetric matrices need con()
                    end do
                end do

#if ldebug
                if (test_diff) then
                    ! back up the V_0 block
                    allocate(loc_block_0_backup(n_mod,n_mod))
                    loc_block_0_backup = loc_block
                end if
#endif

                ! add block to kd + (0,0)
                ierr = insert_block_mat(loc_block,mat,kd,[0,0],n_sol)
                CHCKERR('')
                
                ! -------------!
                ! BLOCKS ~ V_1 !
                ! -------------!
                ! fill local block
                do m = 1,n_mod
                    do k = 1,n_mod
                        loc_block(k,m) = V_1(kd_loc,c([k,m],.false.,n_mod))     ! asymetric matrices don't need con()
                    end do
                end do
                
                ! add block to kd + (0,-p..p) + Hermitian conjugate
                do jd = -norm_disc_prec_sol,norm_disc_prec_sol
                    ierr = insert_block_mat(loc_block*&
                        &norm_disc_coeff(jd+norm_disc_prec_sol+1),mat,kd,&
                        &[0,jd],n_sol,transp=.true.)
                    CHCKERR('')
                end do
                
                ! -------------!
                ! BLOCKS ~ V_2 !
                ! -------------!
                ! fill local block
                do m = 1,n_mod
                    do k = 1,n_mod
                        loc_block(k,m) = con(V_2(kd_loc,c([k,m],.true.,n_mod)),&
                            &[k,m],.true.)                                      ! symmetric matrices need con()
                    end do
                end do
                
#if ldebug
                if (test_diff) then
                    loc_block = loc_block + diff_coeff*loc_block_0_backup
                    deallocate(loc_block_0_backup)
                end if
#endif
            
                ! add block to kd + (-p..p,-p..p)
                do jd = -norm_disc_prec_sol,norm_disc_prec_sol
                    do id = -norm_disc_prec_sol,norm_disc_prec_sol
                        ierr = insert_block_mat(loc_block*&
                            &norm_disc_coeff(id+norm_disc_prec_sol+1)*&
                            &norm_disc_coeff(jd+norm_disc_prec_sol+1),mat,&
                            &kd,[id,jd],n_sol)
                        CHCKERR('')
                    end do
                end do
            end do
            
            ! assemble the matrix
            call MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t begin assembly of matrix')
            call MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY,ierr)
            CHCKERR('Coulnd''t end assembly of matrix')
            
            !call MatSetOption(mat,MAT_HERMITIAN,PETSC_TRUE,ierr)                ! Hermitian to a first approximation
            !CHCKERR('Coulnd''t set option Hermitian')
            
            ! deallocate variables
            deallocate(loc_block)
        end function fill_mat
        
        ! Sets the limits of the indices of the bulk matrix, depending on the BC
        ! style.
        integer function set_bulk_lims(grid_sol,i_lim) result(ierr)
            use num_vars, only: BC_style
            
            character(*), parameter :: rout_name = 'set_bulk_lims'
            
            ! input / output
            type(grid_type), intent(in) :: grid_sol                             ! solution grid
            integer, intent(inout) :: i_lim(2)                                  ! min and max of bulk limits
            
            ! initialize ierr
            ierr = 0
            
            ! set up i_min
            select case (BC_style(1))
                case (1)
                    i_lim(1) = grid_sol%i_min                                   ! will be overwritten
                case (2)
                    i_lim(1) = max(grid_sol%i_min,1+norm_disc_prec_sol)         ! first norm_disc_prec_sol rows write left BC's
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
                    i_lim(2) = grid_sol%i_max                                   ! will be overwritten
                case (2)
                    i_lim(2) = min(grid_sol%i_max,n_r-norm_disc_prec_sol)       ! last norm_disc_prec_sol rows write right BC's
                case (3)
                    i_lim(2) = min(grid_sol%i_max,n_r-1)                        ! last row writes right BC's
                case default
                    err_msg = 'No BC style associated with '//&
                        &trim(i2str(BC_style(1)))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end function set_bulk_lims
        
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
            
            call MatGetInfo(mat,MAT_GLOBAL_SUM,mat_info,ierr)
            CHCKERR('')
            call writo('memory usage: '//&
                &trim(r2strt(mat_info(MAT_INFO_MEMORY)*1.E-6_dp))//' MB')
            call writo('nonzero''s allocated: '//&
                &trim(r2strt(mat_info(MAT_INFO_NZ_ALLOCATED))))
            if (mat_info(MAT_INFO_NZ_UNNEEDED).gt.0._dp) then
                call writo('of which unused: '//&
                    &trim(r2strt(mat_info(MAT_INFO_NZ_UNNEEDED))))
            end if
            
            call lvl_ud(-1)
        end function disp_mat_info
        
        ! Sets nonzero elements d_nz and o_nz.
        subroutine set_nonzeros()
            ! local variables
            PetscInt, allocatable :: tot_nz(:)                                  ! nr. of total non-zeros
            
            ! initialize the numbers of non-zeros in diagonal and off-diagonal
            allocate(tot_nz(loc_n_r*n_mod)); tot_nz = 0
            allocate(d_nz(loc_n_r*n_mod)); d_nz = 0
            allocate(o_nz(loc_n_r*n_mod)); o_nz = 0
            
            ! calculate number of total nonzero entries
            tot_nz = 1+2*st_size
            do kd = 1,st_size
                if (kd.ge.grid_sol%i_min .and. kd.le.grid_sol%i_max) &          ! limit due to left BC
                    &tot_nz((kd-grid_sol%i_min)*n_mod+1:&
                    &(kd-grid_sol%i_min+1)*n_mod) = &
                    &tot_nz((kd-grid_sol%i_min)*n_mod+1:&
                    &(kd-grid_sol%i_min+1)*n_mod) - &
                    &(st_size+1-kd)
            end do
            do kd = n_r,n_r-st_size+1,-1
                if (kd.ge.grid_sol%i_min .and. kd.le.grid_sol%i_max) &          ! limit due to right BC
                    &tot_nz((kd-grid_sol%i_min)*n_mod+1:&
                    &(kd-grid_sol%i_min+1)*n_mod) = &
                    &tot_nz((kd-grid_sol%i_min)*n_mod+1:&
                    &(kd-grid_sol%i_min+1)*n_mod) - &
                    &(kd-n_r+st_size)
            end do
            
            ! calculate number of nonzero diagonal entries
            if (loc_n_r.le.(st_size+1)) then                                    ! up to (st_size+1) normal points in this process
                d_nz = loc_n_r
            else                                                                ! more than (st_size+1) normal points in this process
                do kd = 1,st_size
                    d_nz((kd-1)*n_mod+1:kd*n_mod) = &
                        &st_size+kd
                    d_nz((loc_n_r-kd)*n_mod+1:(loc_n_r-kd+1)*n_mod) = &
                        &st_size+kd
                end do
                d_nz(st_size*n_mod+1:(loc_n_r-st_size)*n_mod) = &
                    &2*st_size+1
            end if
            d_nz = min(d_nz,loc_n_r)                                            ! limit to loc_n_r
            
            ! calculate number of nonzero off-diagonal entries
            do kd = 1,loc_n_r*n_mod
                o_nz(kd) = tot_nz(kd)-d_nz(kd)
            end do
        end subroutine set_nonzeros
        
        ! Defines matrix product for shell matrices: Y = A*X and Y = B*X
        ! Uses the variables A, B, bulk_i_lim, grid_sol from parent routine.
        subroutine shell_A_product(mat,vec_in,vec_out,ierr)
            character(*), parameter :: rout_name = 'shell_A_product'
            
            ! input / output
            Mat :: mat                                                          ! matrix (unused)
            Vec :: vec_in, vec_out                                              ! input and output vector
            PetscInt :: ierr                                                    ! error variable
            
            ! initialize ierr
            ierr = 0
            
            ! call with PV
            ierr = shell_mat_product(X%PV_0(1,i_geo,:,:),&
                &X%PV_1(1,i_geo,:,:),X%PV_2(1,i_geo,:,:),vec_in,vec_out)
            CHCKERR('')
        end subroutine shell_A_product
        subroutine shell_B_product(mat,vec_in,vec_out,ierr)
            character(*), parameter :: rout_name = 'shell_B_product'
            
            ! input / output
            Mat :: mat                                                          ! matrix (unused)
            Vec :: vec_in, vec_out                                              ! input and output vector
            PetscInt :: ierr                                                    ! error variable
            
            ! initialize ierr
            ierr = 0
            
            ! call with KV
            ierr = shell_mat_product(X%KV_0(1,i_geo,:,:),&
                &X%KV_1(1,i_geo,:,:),X%KV_2(1,i_geo,:,:),vec_in,vec_out)
            CHCKERR('')
        end subroutine shell_B_product
        integer function shell_mat_product(V_0,V_1,V_2,vec_in,vec_out) &
            &result(ierr)
            character(*), parameter :: rout_name = 'shell_mat_product'
            
            ! input / output
            PetscScalar, intent(in) :: V_0(:,:)                                 ! either PV_int_0 or KV_int_0 in equilibrium normal grid
            PetscScalar, intent(in) :: V_1(:,:)                                 ! either PV_int_1 or KV_int_1 in equilibrium normal grid
            PetscScalar, intent(in) :: V_2(:,:)                                 ! either PV_int_2 or KV_int_2 in equilibrium normal grid
            Vec, intent(in) :: vec_in                                           ! input vector
            Vec, intent(inout) :: vec_out                                       ! output vector
            
            ! local variables
            PetscScalar, pointer :: vec_in_loc(:)                               ! local pointer to vec_in
            PetscScalar, pointer :: vec_out_loc(:)                              ! local pointer to vec_out
            PetscInt :: kd                                                      ! counter
            PetscInt :: kd_loc                                                  ! kd in local variables
            
            write(*,*) '!!! GOT INTO USER DEFINED SHELL MATRIX PRODUCT !!!!!!!!'
            
            ! initialize ierr
            ierr = 0
            
            ! get vectors
            err_msg = 'Failed to get array'
            call VecGetArrayReadF90(vec_in,vec_in_loc,ierr)
            CHCKERR(err_msg)
            call VecGetArrayF90(vec_out,vec_out_loc,ierr)
            CHCKERR(err_msg)
            
            ! loop over all the diagonal elements of the matrix
            do kd = bulk_i_lim(1)-1,bulk_i_lim(2)-1                             ! (indices start with 0 here)
                ! set up local kd
                kd_loc = kd+2-grid_sol%i_min
                
                ! set diagonal
                !diag_loc(kd) = V_0(
            end do
            !!! TEMPORARY !!!!
            vec_out_loc = 3._dp
            
            ! restore vectors
            err_msg = 'Failed to restore array'
            call VecRestoreArrayReadF90(vec_in,vec_in_loc,ierr)
            CHCKERR(err_msg)
            call VecRestoreArrayF90(vec_out,vec_out_loc,ierr)
            CHCKERR(err_msg)
        end function shell_mat_product
        
        ! Gets diagonal of shell matrix A or B.
        subroutine shell_A_diagonal(mat,diag,ierr)
            character(*), parameter :: rout_name = 'shell_A_diagonal'
            
            ! input / output
            Mat :: mat                                                          ! matrix (unused)
            Vec :: diag                                                         ! diagonal
            PetscInt :: ierr                                                    ! error variable
            
            ! call with PV
            ierr = shell_mat_diagonal(X%PV_0(1,i_geo,:,:),diag)
            CHCKERR('')
        end subroutine shell_A_diagonal
        subroutine shell_B_diagonal(mat,diag,ierr)
            character(*), parameter :: rout_name = 'shell_B_diagonal'
            
            ! input / output
            Mat :: mat                                                          ! matrix (unused)
            Vec :: diag                                                         ! diagonal
            PetscInt :: ierr                                                    ! error variable
            
            ! initialize ierr
            ierr = 0
            
            ! call with KV
            ierr = shell_mat_diagonal(X%KV_0(1,i_geo,:,:),diag)
            CHCKERR('')
        end subroutine shell_B_diagonal
        integer function shell_mat_diagonal(V_0,diag) result(ierr)
            character(*), parameter :: rout_name = 'shell_mat_diagonal'
            
            ! input / output
            PetscScalar, intent(in) :: V_0(:,:)                                 ! either PV_int_0 or KV_int_0 in equilibrium normal grid
            Vec, intent(inout) :: diag                                          ! diagonal
            
            ! local variables
            PetscInt :: kd                                                      ! counter
            PetscInt :: kd_loc                                                  ! kd in local variables
            PetscScalar, pointer :: diag_loc(:)                                 ! local pointer to diag
            
            ! initialize ierr
            ierr = 0
            
            write(*,*) '!!! GOT INTO USER DEFINED SHELL MATRIX PRODUCT !!!!!!!!'
            
            ! get vectors
            err_msg = 'Failed to get array'
            call VecGetArrayF90(diag,diag_loc,ierr)
            CHCKERR(err_msg)
            
            ! loop over all the diagonal elements of the matrix
            do kd = bulk_i_lim(1)-1,bulk_i_lim(2)-1                             ! (indices start with 0 here)
                ! set up local kd
                kd_loc = kd+2-grid_sol%i_min
                
                ! set diagonal
                !diag_loc(kd) = V_0(
            end do
            
            ! restore vectors
            err_msg = 'Failed to restore array'
            call VecRestoreArrayF90(diag,diag_loc,ierr)
            CHCKERR(err_msg)
        end function shell_mat_diagonal
    end function setup_mats
    
    ! Sets the  boundary conditions:
    ! Deep  in the  plasma,  an  artificial Eigenvalue  EV_BC  is introduced  by
    ! setting the  diagonal components  of A  to EV_BC and  of B  to 1,  and the
    ! off-diagonal elements to zero.
    ! At the plasma surface, the surface energy is minimized as in [ADD REF].
    integer function set_BC(grid_sol,X,A,B,i_geo,n_sol,norm_disc_coeff) &
        &result(ierr)
        use num_vars, only: norm_disc_prec_sol, BC_style
        use MPI_utilities, only: get_ser_var, wait_MPI
        use num_utilities, only: con, c
        use SLEPC_utilities, only: insert_block_mat
        
        character(*), parameter :: rout_name = 'set_BC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        Mat, intent(inout) :: A, B                                              ! Matrices A and B from A X = lambda B X
        integer, intent(in) :: i_geo                                            ! at which geodesic index to perform the calculations
        integer, intent(in) :: n_sol                                            ! number of grid points of solution grid
        PetscReal, intent(in) :: norm_disc_coeff(:)                             ! discretization coefficients for normal derivatives
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscInt :: n_min, n_max                                                ! absolute limits excluding the BC's
        
        ! local variables
        integer :: n_mod                                                        ! nr. of modes
        PetscInt :: kd                                                          ! counter
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        call writo('Preparing variables')
        call lvl_ud(1)
        
        ! set nr. of modes
        n_mod = X%n_mod(1)
        if (n_mod.ne.X%n_mod(2)) then
            ierr = 1
            CHCKERR('Need square matrix')
        end if
        
        call lvl_ud(-1)
        
        ! set up n_min, depending on BC style
        select case (BC_style(1))
            case (1)
                n_min = 2*norm_disc_prec_sol                                    ! dirichlet BC requires half stencil
            case (2)
                n_min = norm_disc_prec_sol                                      ! mixed BC requires only one fourth of stencil
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
                n_max = 2*norm_disc_prec_sol                                    ! dirichlet BC requires half stencil
            case (2)
                n_max = norm_disc_prec_sol                                      ! mixed BC requires only one fourth of stencil
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
        
        ! iterate over all positions where to set left BC
        do kd = 1,n_min
            if (grid_sol%i_min.le.kd .and. grid_sol%i_max.ge.kd) then           ! decide which process does the setting
                select case (BC_style(1))
                    case (1)
                        ierr = set_BC_1(kd-1,A,B,.false.)                       ! indices start at 0
                        CHCKERR('')
                    case (2)
                        ierr = set_BC_2(kd-1,kd-grid_sol%i_min+1,X,A,B,&
                            &i_geo,grid_sol%n(3),norm_disc_coeff,.false.)       ! indices start at 0
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
        do kd = grid_sol%n(3),grid_sol%n(3)-n_max+1,-1
            if (grid_sol%i_min.le.kd .and. grid_sol%i_max.ge.kd) then           ! decide which process does the setting
                select case (BC_style(2))
                    case (1)
                        ierr = set_BC_1(kd-1,A,B,.true.)                        ! indices start at 0
                        CHCKERR('')
                    case (2)
                        ierr = set_BC_2(kd-1,kd-grid_sol%i_min+1,X,A,B,i_geo,&
                            &grid_sol%n(3),norm_disc_coeff,.true.)              ! indices start at 0
                        CHCKERR('')
                    case (3)
                        if (kd.eq.grid_sol%n(3)) then                           ! only for last point, irrespective of discretization order
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
            use num_vars, only: EV_BC, norm_disc_prec_sol
            
            character(*), parameter :: rout_name = 'set_BC_1'
            
            ! input / output
            integer, intent(in) :: r_id                                         ! position at which to set BC
            Mat, intent(inout) :: A, B                                          ! Matrices A and B from A X = lambda B X
            logical :: BC_right                                                 ! if BC is at right (so there are vacuum terms)
            
            ! local variables
            integer :: kd                                                       ! counter
            PetscScalar, allocatable :: loc_block(:,:)                          ! (n_mod,n_mod) block matrix for 1 normal point
            PetscInt :: pmone                                                   ! plus or minus one
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(r_id+1))//&
                &': Eigenvector set to zero',persistent=.true.)
            
            ! initialize local blocks
            allocate(loc_block(n_mod,n_mod))
            loc_block = 0.0_dp
            do kd = 1,n_mod
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
            do kd = 1,2*norm_disc_prec_sol
                ! set block r_id +/- (0,kd) and Hermitian conjugate
                ierr = insert_block_mat(0*loc_block,A,r_id,[0,-pmone*kd],&
                    &n_sol,overwrite=.true.,transp=.true.)
                CHCKERR('')
                ierr = insert_block_mat(0*loc_block,B,r_id,[0,-pmone*kd],&
                    &n_sol,overwrite=.true.,transp=.true.)
                CHCKERR('')
            end do
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
            allocate(V_int_0_mod(n_mod,n_mod,2))
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
                do jd = -norm_disc_prec_sol,norm_disc_prec_sol
                    if (r_id.lt.n_sol .and. r_id+jd.lt.n_sol) then
                        ierr = insert_block_mat(-X%vac_res*&
                            &norm_disc_coeff(jd+norm_disc_prec_sol+1),A,&
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
            
            ! initialize ierr
            ierr = 0
            
            ! initialize variables
            allocate(KV_2_inv(n_mod,n_mod))
            allocate(V_triple(n_mod,n_mod))
            allocate(KV_1_loc(n_mod,n_mod))
            allocate(PV_1_loc(n_mod,n_mod))
            
            ! make local copy of KV_2 into the inverse array
            KV_2_inv = 0._dp
            do m = 1,n_mod
                do k = 1,m                                                      ! only save upper diagonal part
                    KV_2_inv(k,m) = con(KV_2(c([k,m],.true.,n_mod)),&
                        &[k,m],.true.)                                          ! symmetric matrices need con()
                end do
            end do
            
#if ldebug
            if (debug_calc_V_0_mod) then
                write(*,*) 'KV_2 = '
                
                write(*,*) 'real part:'
                call print_ar_2(realpart(KV_2_inv))
                write(*,*) 'imaginary part:'
                call print_ar_2(imagpart(KV_2_inv))
            end if
#endif
            
            ! invert matrix KV_2 using Lapack
            uplo = 'U'                                                          ! only lower diagonal matters, but is arbitrary
            call zpotrf(uplo,n_mod,KV_2_inv,n_mod,ierr)
            CHCKERR('Failed to decompose KV_2')
            call zpotri(uplo,n_mod,KV_2_inv,n_mod,ierr)
            CHCKERR('Failed to invert KV_2')
            do k = 1,n_mod
                do m = 1,k-1
                    KV_2_inv(k,m) = conjg(KV_2_inv(m,k))
                end do
            end do
            
            ! save local copy of KV_1 and PV_1
            do m = 1,n_mod
                do k = 1,n_mod
                    KV_1_loc(k,m) = KV_1(c([k,m],.false.,n_mod))                ! asymetric matrices don't need con()
                    PV_1_loc(k,m) = PV_1(c([k,m],.false.,n_mod))                ! asymetric matrices don't need con()
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
            do m = 1,n_mod
                do k = 1,n_mod
                    V_0_mod(k,m,1) = &
                        &con(PV_0(c([k,m],.true.,n_mod)),[k,m],.true.)          ! symmetric matrices need con()
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
            do m = 1,n_mod
                do k = 1,n_mod
                    V_0_mod(k,m,2) = &
                        &con(KV_0(c([k,m],.true.,n_mod)),[k,m],.true.)          ! symmetric matrices need con()
                end do
            end do
            V_0_mod(:,:,2) = V_0_mod(:,:,2) - V_triple
#if ldebug
        if (debug_calc_V_0_mod) then
            write(*,*) 'KV_2^-1 = '
            
            write(*,*) 'real part:'
            call print_ar_2(realpart(KV_2_inv))
            write(*,*) 'imaginary part:'
            call print_ar_2(imagpart(KV_2_inv))
            
            write(*,*) 'PV_1 = '
            
            write(*,*) 'real part:'
            call print_ar_2(realpart(PV_1_loc))
            write(*,*) 'imaginary part:'
            call print_ar_2(imagpart(PV_1_loc))
            
            write(*,*) 'vac = '
            
            write(*,*) 'real part:'
            call print_ar_2(realpart(X%vac_res))
            write(*,*) 'imaginary part:'
            call print_ar_2(imagpart(X%vac_res))
            
            write(*,*) 'KV_1 = '
            
            write(*,*) 'real part:'
            call print_ar_2(realpart(KV_1_loc))
            write(*,*) 'imaginary part:'
            call print_ar_2(imagpart(KV_1_loc))
            
            write(*,*) 'PV_0_mod = '
            
            write(*,*) 'real part:'
            call print_ar_2(realpart(V_0_mod(:,:,1)))
            write(*,*) 'imaginary part:'
            call print_ar_2(imagpart(V_0_mod(:,:,1)))
            
            write(*,*) 'KV_0_mod = '
            
            write(*,*) 'real part:'
            call print_ar_2(realpart(V_0_mod(:,:,2)))
            write(*,*) 'imaginary part:'
            call print_ar_2(imagpart(V_0_mod(:,:,2)))
        end if
#endif
        end function calc_V_0_mod
    end function set_BC
    
    ! sets up EV solver
    integer function setup_solver(grid_sol,X,A,B,solver) result(ierr)
        use num_vars, only: n_sol_requested, tol_SLEPC, max_it_slepc
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'setup_solver'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        EPS, intent(inout) :: solver                                            ! EV solver
        
        ! local variables
        integer :: n_mod                                                        ! nr. of modes
        PetscInt :: n_sol                                                       ! how many solutions can be requested (normally n_sol_requested)
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! set nr. of modes
        n_mod = X%n_mod(1)
        if (n_mod.ne.X%n_mod(2)) then
            ierr = 1
            CHCKERR('Need square matrix')
        end if
        
        !call PetscOptionsSetValue('-eps_view','-1',ierr)
        call EPSCreate(PETSC_COMM_WORLD,solver,ierr)
        CHCKERR('EPSCreate failed')
        
        call EPSSetOperators(solver,A,B,ierr)                                   ! generalized EV problem A X = lambda B X
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
        if (n_sol_requested.gt.grid_sol%n(3)*n_mod) then
            call writo('max. nr. of solutions requested capped to problem &
                &dimension ('//trim(i2str(grid_sol%n(3)*n_mod))//')',&
                &warning=.true.)
            call writo('Increase either n_r_sol or number of pol. modes or &
                &decrease n_sol_requested')
            n_sol = grid_sol%n(3)*n_mod
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
        call writo('number of iterations: '//trim(i2str(n_it)))
        call writo('number of converged solutions: '//trim(i2str(n_conv)))
        call writo('number of requested solutions: '//trim(i2str(n_ev)))
        call writo('maximum dimension of the subspace to be used by solver: '//&
            &trim(i2str(ncv)))
        call writo('maximum dimension allowed for projected problem : '//&
            &trim(i2str(mpd)))
        
        ! set maximum nr of solutions to be saved
        if (n_sol_requested.gt.n_conv) then
            call writo('max. nr. of solutions found only '//&
                &trim(i2str(n_conv)),warning=.true.)
            max_n_EV = n_conv
        else
            max_n_EV = n_sol_requested
        end if
        
        call lvl_ud(-1)
    end function summarize_solution
    
    ! stores the results
#if ldebug
    integer function store_results(grid_sol,X,sol,solver,max_n_EV,A,B,&
        &step_size) result(ierr)
#else
    integer function store_results(grid_sol,X,sol,solver,max_n_EV) result(ierr)
#endif
        use eq_vars, only: T_0
        use num_vars, only: use_normalization, EV_BC, prog_name, output_name, &
            &n_procs, rank, tol_SLEPC, eq_style
        use files_utilities, only: nextunit
        use MPI_utilities, only: get_ser_var
        use sol_vars, only: create_sol
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'store_results'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
        type(X_2_type), intent(in) :: X                                         ! field-averaged perturbation variables (so only first index)
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(inout) :: max_n_EV                                     ! nr. of EV's saved, up to n_conv
#if ldebug
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        PetscReal, intent(in) :: step_size                                      ! step size in flux coordinates
#endif
        
        ! local variables
        PetscInt :: id, kd                                                      ! counters
        PetscInt :: id_tot                                                      ! counters
        Vec :: sol_vec                                                          ! solution EV parallel vector
        PetscInt :: one = 1                                                     ! one
        PetscInt, allocatable :: sol_vec_max_loc(:)                             ! location of sol_vec_max of this process
        PetscReal :: error                                                      ! error of EPS solver
        PetscReal :: tol_complex                                                ! tolerance for an EV to be considered complex
        PetscReal :: tol_EV_BC = 1.E-3_dp                                       ! tolerance for an EV to be considered due to the BC
        PetscReal, parameter :: infinity = 1.E40_dp                             ! beyond this value, modes are not saved
        PetscScalar :: sol_val_loc                                              ! local X val
        PetscScalar, allocatable :: sol_vec_max(:)                              ! max of sol_vec of all processes, including phase at max
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: format_val                                 ! format
        character(len=max_str_ln) :: format_head                                ! header
        character(len=max_str_ln) :: EV_err_str                                 ! String with information about EV error
        integer :: output_EV_i                                                  ! file number
        integer :: n_digits                                                     ! nr. of digits for the integer number
        integer :: n_err(3)                                                     ! how many errors there were
        character(len=max_str_ln) :: err_msg                                    ! error message
#if ldebug
        Mat :: err_mat                                                          ! A - omega^2 B
        Vec :: err_vec                                                          ! AX - omega^2 BX
        Vec :: E_vec(2)                                                         ! AX and BX
        PetscReal :: err_norm                                                   ! norm of error
        PetscReal :: error_est                                                  ! error estimate of EPS solver
        PetscScalar :: err_val                                                  ! X*(A-omega^2B)X
        PetscScalar :: E_val(2)                                                 ! X*AX and X*BX
#endif
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! tests
        if (X%n_mod(1).ne.X%n_mod(2)) then
            ierr = 1
            err_msg = 'Complete matrix with all the Fourier modes needed'
            CHCKERR(err_msg)
        end if
        do id = 1,X%n_mod(1)
            do kd = 1,grid_sol%loc_n_r
                if (X%n_1(kd,id).ne.X%n_2(kd,id) .or. &
                    &X%m_1(kd,id).ne.X%m_2(kd,id)) then
                    ierr = 1
                    err_msg = 'Fourier modes do not match in both dimensions'
                    CHCKERR(err_msg)
                end if
            end do
        end do
        
        ! create solution variables
        call create_sol(grid_sol,sol,max_n_EV)
        
        ! create solution vector
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,&
            &grid_sol%loc_n_r*sol%n_mod,grid_sol%n(3)*sol%n_mod,&
            &PETSC_NULL_SCALAR,sol_vec,ierr)
        CHCKERR('Failed to create MPI vector with arrays')
        
        ! set up EV error string and format string:
        !   1: index of EV
        !   2: Real part
        !   3: Imaginary part
        !   4: relative error ||Ax - EV Bx||_2/||EV x||_2 [1]
        n_err = 0
        EV_err_str = ''
        if (rank.eq.0) then
            n_digits = ceiling(log10(1._dp*max_n_EV))
            format_val = '(I'//trim(i2str(n_digits))//'," ",ES23.16," ",&
                &ES23.16," ",ES23.16)'
            format_head = '(A'//trim(i2str(n_digits+3))//'," ",A23," ",A23," ",&
                &A23)'
        end if
        
        ! open output file for the log
        if (rank.eq.0) then
            full_output_name = prog_name//'_'//trim(output_name)//'_EV_R_'//&
                &trim(i2str(rich_lvl))//'.txt'
            open(unit=nextunit(output_EV_i),file=full_output_name,iostat=ierr)
            CHCKERR('Cannot open EV output file')
        end if
        
        ! output to file
        if (rank.eq.0) then
            if (use_normalization) then
                write(output_EV_i,'(A)') '# Eigenvalues normalized to the &
                    &squared  Alfven frequency omega_A^2 = '
                select case (eq_style)
                    case (1)                                                    ! VMEC
                        write(output_EV_i,'(A)') '#     ('//&
                            &trim(r2str(1._dp/T_0))//' Hz)^2 = '//&
                            &trim(r2str(1._dp/T_0**2))//' Hz^2'
                    case (2)                                                    ! HELENA
                        write(output_EV_i,'(A)') '#     (HELENA normalization)'
                end select
            else
                write(output_EV_i,'(A)') '# Eigenvalues'
            end if
            write(output_EV_i,format_head) '#  I                            ', &
                &'real part               ', 'imaginary part          ', &
                &'relative precision      '
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
            call EPSGetErrorEstimate(solver,id_tot-1,error_est)                 ! get error estimate
            CHCKERR('EPSGetErrorEstimate failed')
            
            ! set up local solution val
            sol_val_loc = sol%val(id)
            
            ! tests
            tol_complex = tol_SLEPC(rich_lvl)
            if (abs(imagpart(sol%val(id))/realpart(sol%val(id))).gt.&
                &tol_complex) then                                              ! test for unphysical complex solution
                EV_err_str = '# WARNING: Unphysical complex Eigenvalue!'
                n_err(1) = n_err(1)+1
                call remove_EV(id,max_n_EV)
            else if (abs(realpart(sol%val(id)-EV_BC)/EV_BC).lt.tol_EV_BC) then  ! test for artificial EV due to BC
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
            
            ! if error, comment the next line
            if (EV_err_str.ne.'') then
                if (rank.eq.0) write(output_EV_i,'(A)',advance='no') '# '
            else
                if (rank.eq.0) write(output_EV_i,'(A)',advance='no') '  '
            end if
            
            if (rank.eq.0) then
                ! output EV to file
                write(output_EV_i,format_val) id,realpart(sol_val_loc),&
                    &imagpart(sol_val_loc),error
                
                ! if error, print explanation
                if (EV_err_str.ne.'') write(output_EV_i,'(A)') trim(EV_err_str)
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
            end if
            
            if (EV_err_str.eq.'') then
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
                call MatAXPY(err_mat,-sol%val(id),B,SAME_NONZERO_PATTERN,ierr)
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
                
                call lvl_ud(-1)
                
                ! clean up
                call VecDestroy(err_vec,ierr)                                   ! destroy error vector
                CHCKERR('Failed to destroy err_vec')
            end if
            
            ! reinitialize error string if error and increment counter if not
            if (EV_err_str.ne.'') then
                EV_err_str = ''
            else                                                                ! tests passed
                id = id + 1
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
                allocate(sol_vec_loc(sol%n_mod,grid_sol%loc_n_r,max_id))
                sol_val_loc = sol%val
                sol_vec_loc = sol%vec
                
                ! reallocate arrays
                deallocate(sol%val,sol%vec)
                if (remove_next_loc) then
                    allocate(sol%val(id-1))
                    allocate(sol%vec(sol%n_mod,grid_sol%loc_n_r,id-1))
                else
                    allocate(sol%val(max_id-1))
                    allocate(sol%vec(sol%n_mod,grid_sol%loc_n_r,max_id-1))
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
    integer function stop_SLEPC(grid_sol,A,B,solver) result(ierr)
        character(*), parameter :: rout_name = 'stop_SLEPC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_sol                                 ! solution grid
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
