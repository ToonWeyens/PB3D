!------------------------------------------------------------------------------!
!   Operations that use SLEPC (and PETSC) routines                             !
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
    use grid_vars, only: grid_type, dealloc_grid
    use X_vars, only: X_1_type, X_2_type
    use sol_vars, only: sol_type

    implicit none
    private
    public solve_EV_system_SLEPC
#if ldebug
    public debug_setup_mats, debug_set_BC, debug_insert_block_mat
#endif
    
    ! global variables
#if ldebug
    logical :: debug_setup_mats = .false.                                       ! plot debug information for setup_mats
    logical :: debug_set_BC = .false.                                           ! plot debug information for set_BC
    logical :: debug_insert_block_mat = .false.                                 ! plot debug information for insert_block_mat
    logical :: debug_calc_V_0_mod = .false.                                     ! plot debug information for calc_V_0_mod
    logical :: debug_store_results = .false.                                    ! plot debug information for store_results
    logical :: test_diff = .false.                                              ! test introduction of numerical diff
    real(dp) :: diff_coeff                                                      ! diff coefficient
#endif
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in  [ADD REF] and  solves them using  the SLEPC suite.  The most
    ! unstable solutions  are obtained, where the  variable "max_n_EV" indicates
    ! how many.
    ! Optionally the  geodesic index of  the perturbation variables at  which to
    ! perform the calculations can be changed  from its default value of 1 using
    ! the variable i_geo.
    integer function solve_EV_system_SLEPC(grid_eq,grid_X,X,sol,use_guess,&
        &max_n_EV,i_geo) result(ierr)
        use num_vars, only: max_it_inv, norm_disc_prec_X
        use grid_ops, only: get_norm_interp_data, trim_grid
        use grid_vars, only: dealloc_grid
        use utilities, only: calc_coeff_fin_diff
#if ldebug
        use num_vars, only: ltest
        use input_ops, only: get_real, get_log
#endif
        
        character(*), parameter :: rout_name = 'solve_EV_system_SLEPC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(in) :: X                                         ! perturbation variables
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        PetscBool, intent(in) :: use_guess                                      ! whether to use a guess or not
        PetscInt, intent(inout) :: max_n_EV                                     ! how many solutions saved
        integer, intent(in), optional :: i_geo                                  ! at which geodesic index to perform the calculations
        
        ! local variables
        type(grid_type) :: grid_X_trim                                          ! trimmed X grid
        Mat, target :: A                                                        ! matrix A in EV problem A X = lambda B X
        Mat, target :: B                                                        ! matrix B in EV problem A X = lambda B X
        EPS :: solver                                                           ! EV solver
        PetscInt, save :: guess_start_id = -10                                  ! start of index of previous vector, saved for next iteration
        PetscInt, save :: prev_n_EV                                             ! nr. of solutions of previous vector
        PetscInt :: inv_lvl_nr                                                  ! level of inverse calculation
        PetscInt :: i_geo_loc                                                   ! local copy of i_geo
        PetscReal, allocatable :: loc_r_eq(:)                                   ! unrounded index in tables V_int
        PetscReal, allocatable :: norm_disc_coeff(:)                            ! discretization coefficients
        PetscReal :: step_size                                                  ! step size in flux coordinates
        PetscBool :: done_inverse                                               ! is it converged?
        
        ! initialize ierr
        ierr = 0
        
        ! initialization
        call writo('initialization...')
        
        ! trim X grid
        ierr = trim_grid(grid_X,grid_X_trim)
        CHCKERR('')
        ! set up local i_geo
        i_geo_loc = 1
        if (present(i_geo)) i_geo_loc = i_geo
        ! set up arrays loc_r_eq that determine how to fill the matrices
        ierr = get_norm_interp_data(grid_eq,grid_X_trim,loc_r_eq)
        CHCKERR('')
        ! start SLEPC
        ierr = start_SLEPC(grid_X_trim%n(3))
        CHCKERR('')
        ! set up step size
        step_size = (grid_X%r_F(grid_X%n(3))-grid_X%r_F(1))/(grid_X%n(3)-1)
        
        ! get discretization coefficients
        call writo('get discretization coefficients...')
        call lvl_ud(1)
        ierr = calc_coeff_fin_diff(1,norm_disc_prec_X,norm_disc_coeff)
        CHCKERR('')
        norm_disc_coeff = norm_disc_coeff/step_size                             ! scale by step size
        call lvl_ud(-1)
        
        ! set up the matrix
        call writo('set up matrices...')
        
        ierr = setup_mats(grid_X_trim,X,A,B,loc_r_eq,i_geo_loc,norm_disc_coeff)
        CHCKERR('')
        
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
        
        ! set up guess
        call writo('set up guess...')
        
        if (use_guess) call setup_guess(sol,A,solver,guess_start_id,prev_n_EV)
        
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
            call writo('set up boundary conditions...')
            
            ierr = set_BC(grid_X_trim,X,A,B,loc_r_eq,i_geo_loc,grid_X%n(3),&
                &norm_disc_coeff)
            CHCKERR('')
            
#if ldebug
            if (debug_set_BC) then
                call writo('Testing if AFTER INSERTING BC''s, A and B are &
                    &Hermitian by multiplying with their Hermitian Transpose &
                    &and subtracting 1')
                call lvl_ud(1)
                
                ierr = test_mat_hermiticity(A,'A_BC')
                CHCKERR('')
                
                ierr = test_mat_hermiticity(B,'B_BC')
                CHCKERR('')
                
                call lvl_ud(-1)
            end if
#endif
            
            ! set up solver
            call writo('set up EV solver...')
#if ldebug
            if (ltest) then
                call writo('Test spectrum of A or B instead of solving &
                    &generalized Eigenvalue problem?')
                if (get_log(.false.)) then
                    call writo('Spectrum of A (true) or B (false)?')
                    if (get_log(.true.)) then                                   ! A
                        ierr = setup_solver(grid_X_trim,X,A,PETSC_NULL_OBJECT,&
                            &solver)
                        CHCKERR('')
                    else                                                        ! B
                        ierr = setup_solver(grid_X_trim,X,B,PETSC_NULL_OBJECT,&
                            &solver)
                        CHCKERR('')
                    end if
                else
                    ierr = setup_solver(grid_X_trim,X,A,B,solver)
                    CHCKERR('')
                end if
            else
#endif
                ierr = setup_solver(grid_X_trim,X,A,B,solver)
                CHCKERR('')
#if ldebug
            end if
#endif
            
            ! get solution
            call writo('get solution...')
            
            ierr = get_solution(solver)
            CHCKERR('')
            
            ! check for convergence
            inv_lvl_nr = inv_lvl_nr+1
            !!! TO BE IMPLEMENTED !!!
            
            if (max_it_inv.gt.1) call lvl_ud(-1)
        end do Inverse
        
        if (max_it_inv.gt.1) then
            call lvl_ud(-1)
            call writo('Done iterating over the inverse')
        end if
        
        ! summarize solution
        call writo('summarize solution...')
        
        ierr = summarize_solution(solver,max_n_EV)
        CHCKERR('')
        
        ! store results
        call writo('storing results for '//trim(i2str(max_n_EV))//' least &
            &stable Eigenvalues...')
        
#if ldebug
        ierr = store_results(grid_X_trim,sol,solver,max_n_EV,A,B)
        CHCKERR('')
#else
        ierr = store_results(grid_X_trim,sol,solver,max_n_EV)
        CHCKERR('')
#endif
        
        ! finalize
        call writo('finalize SLEPC...')
        
        ! stop SLEPC
        ierr = stop_SLEPC(grid_X_trim,A,B,solver,guess_start_id,prev_n_EV,&
            &max_n_EV)
        CHCKERR('')
        
        ! clean up
        call dealloc_grid(grid_X_trim)
        deallocate(loc_r_eq)
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
    integer function start_SLEPC(n_r_X) result(ierr)
        use num_vars, only: n_procs
        use files_ops, only: opt_args
        
        character(*), parameter :: rout_name = 'start_SLEPC'
        
        ! input / output
        integer, intent(in) :: n_r_X                                            ! nr. of normal points in perturbation grid
        
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
        if (n_procs.gt.n_r_X) then                                              ! too many processors
            call writo('WARNING: using too many processors: '&
                &//trim(i2str(n_procs))//', while beyond '//&
                &trim(i2str(n_r_X))//' does not bring improvement')
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
    integer function setup_mats(grid_X,X,A,B,loc_r_eq,i_geo,norm_disc_coeff) &
        &result(ierr)
        use num_vars, only: norm_disc_prec_X
#if ldebug
        use num_vars, only: ltest
        use input_ops, only: get_real, get_log
#endif
        
        character(*), parameter :: rout_name = 'setup_mats'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(in) :: X                                         ! perturbation variables
        Mat, intent(inout) :: A, B                                              ! matrix A and B
        PetscReal, intent(in) :: loc_r_eq(:)                                    ! unrounded index in tables V_int
        integer, intent(in) :: i_geo                                            ! at which geodesic index to perform the calculations
        PetscReal, intent(in) :: norm_disc_coeff(:)                             ! discretization coefficients for normal derivatives
        
        ! local variables
        PetscInt :: kd                                                          ! counter
        PetscInt :: n_r                                                         ! n_r of trimmed X grid
        PetscInt :: loc_n_r                                                     ! loc_n_r of trimmed X grid
        PetscInt :: st_size                                                     ! half stencil size
        PetscInt, allocatable :: tot_nz(:)                                      ! nr. of total non-zeros
        PetscInt, allocatable :: d_nz(:)                                        ! nr. of diagonal non-zeros
        PetscInt, allocatable :: o_nz(:)                                        ! nr. of off-diagonal non-zeros
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) '!!!!!!! SETUP MATS HAS TO BE REWRITTEN !!!!!!!!!!!!!!!!!!!!'
        call lvl_ud(1)
        
        ! user output
        call writo('Normal discretization with central finite differences of &
            &order '//trim(i2str(norm_disc_prec_X))//', stencil width '//&
            &trim(i2str(4*norm_disc_prec_X+1)))
        
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
        
        ! set up loc_n_r and n_r
        loc_n_r = grid_X%loc_n_r
        n_r = grid_X%n(3)
        
        ! set (half) stencil size
        st_size = 2*norm_disc_prec_X
        
        ! create a  matrix A and B  with the appropriate number  of preallocated
        ! entries (excluding ghost regions)
        ! initialize the numbers of non-zeros in diagonal and off-diagonal
        allocate(tot_nz(loc_n_r*X%n_mod)); tot_nz = 0
        allocate(d_nz(loc_n_r*X%n_mod)); d_nz = 0
        allocate(o_nz(loc_n_r*X%n_mod)); o_nz = 0
        
        ! calculate number of total nonzero entries
        tot_nz = 1+2*st_size
        do kd = 1,st_size
            if (kd.ge.grid_X%i_min .and. kd.le.grid_X%i_max) &                  ! limit due to left BC
                &tot_nz((kd-grid_X%i_min)*X%n_mod+1:&
                &(kd-grid_X%i_min+1)*X%n_mod) = &
                &tot_nz((kd-grid_X%i_min)*X%n_mod+1:&
                &(kd-grid_X%i_min+1)*X%n_mod) - &
                &(st_size+1-kd)
        end do
        do kd = n_r,n_r-st_size+1,-1
            if (kd.ge.grid_X%i_min .and. kd.le.grid_X%i_max) &                  ! limit due to right BC
                &tot_nz((kd-grid_X%i_min)*X%n_mod+1:&
                &(kd-grid_X%i_min+1)*X%n_mod) = &
                &tot_nz((kd-grid_X%i_min)*X%n_mod+1:&
                &(kd-grid_X%i_min+1)*X%n_mod) - &
                &(kd-n_r+st_size)
        end do
        
        ! calculate number of nonzero diagonal entries
        if (loc_n_r.le.(st_size+1)) then                                        ! up to (st_size+1) normal points in this process
            d_nz = loc_n_r
        else                                                                    ! more than (st_size+1) normal points in this process
            do kd = 1,st_size
                d_nz((kd-1)*X%n_mod+1:kd*X%n_mod) = &
                    &st_size+kd
                d_nz((loc_n_r-kd)*X%n_mod+1:(loc_n_r-kd+1)*X%n_mod) = &
                    &st_size+kd
            end do
            d_nz(st_size*X%n_mod+1:(loc_n_r-st_size)*X%n_mod) = &
                &2*st_size+1
        end if
        d_nz = min(d_nz,loc_n_r)                                                ! limit to loc_n_r
        
        ! calculate number of nonzero off-diagonal entries
        do kd = 1,loc_n_r*X%n_mod
            o_nz(kd) = tot_nz(kd)-d_nz(kd)
        end do
        
        ! create matrix A
        call MatCreateAIJ(PETSC_COMM_WORLD,loc_n_r*X%n_mod,loc_n_r*X%n_mod,&
            &grid_X%n(3)*X%n_mod,grid_X%n(3)*X%n_mod,PETSC_NULL_INTEGER,&
            &d_nz*X%n_mod,PETSC_NULL_INTEGER,o_nz*X%n_mod,A,ierr)
        CHCKERR('MatCreateAIJ failed for matrix A')
        
        ! deallocate tot_nz, d_nz and o_nz
        deallocate(tot_nz,d_nz,o_nz)
        
        ! fill the matrix A
        ierr = fill_mat(X%PV_int_0(:,i_geo,:),X%PV_int_1(:,i_geo,:),&
            &X%PV_int_2(:,i_geo,:),n_r,norm_disc_coeff,A)
        CHCKERR('')
        call writo('matrix A set up:')
        
        call lvl_ud(1)
        
        ierr = disp_mat_inf(A)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! duplicate A into B
        ! (the  advantage of  letting communication  and calculation  overlap by
        ! calculating matrix B while assembling A is offset by the extra cost of
        ! not being able to use MatDuplicate. This is also easier)
        call MatDuplicate(A,MAT_SHARE_NONZERO_PATTERN,B,ierr)                   ! B has same structure as A
        CHCKERR('failed to duplicate A into B')
        
        ! fill the matrix B
        ierr = fill_mat(X%KV_int_0(:,i_geo,:),X%KV_int_1(:,i_geo,:),&
            &X%KV_int_2(:,i_geo,:),n_r,norm_disc_coeff,B)
        CHCKERR('')
        call writo('matrix B set up:')
        
        call lvl_ud(1)
        
        ierr = disp_mat_inf(B)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        call lvl_ud(-1)
    contains
        ! fills a  matrix according to norm_disc_prec_X [ADD REF]:
        !   1. first order accuracy for all terms
        !   2. higher order accuracy for internal first order derivatives
        ! it is used  for both the matrix  A and B, corresponding  to the plasma
        ! potential energy and the (perpendicular) kinetic energy
        ! The procedure is as follows:
        !   1. Tables  are set  up for  (k,m) pairs in the equilibrium grid.
        !   2. At the first normal point  in the perturbation grid, belonging to
        !   the current process, the tables are interpolated as a start.
        !   3. At every normal point in  the perturbation grid, belonging to the
        !   current process,  the tables are  interpolated at the next  point in
        !   the  corresponding  perturbation  grid. This  information,  as  well
        !   as  the interpolated  value  of  the previous  point  allow for  the
        !   calculation of every quantity.
        integer function fill_mat(V_0,V_1,V_2,n_X,norm_disc_coeff,mat) &
            &result(ierr)
            use num_vars, only: norm_disc_prec_X, BC_style
            use utilities, only: c, con
            
            character(*), parameter :: rout_name = 'fill_mat'
            
            ! input / output
            PetscScalar, intent(in) :: V_0(:,:)                                 ! either PV_int_0 or KV_int_0 in equilibrium normal grid
            PetscScalar, intent(in) :: V_1(:,:)                                 ! either PV_int_1 or KV_int_1 in equilibrium normal grid
            PetscScalar, intent(in) :: V_2(:,:)                                 ! either PV_int_2 or KV_int_2 in equilibrium normal grid
            integer, intent(in) :: n_X                                          ! number of grid points of perturbation grid
            PetscReal, intent(in) :: norm_disc_coeff(:)                         ! discretization coefficients for normal derivatives
            Mat, intent(inout) :: mat                                           ! either A or B
            
            ! local variables (not to be used in child routines)
            character(len=max_str_ln) :: err_msg                                ! error message
            PetscScalar, allocatable :: loc_block(:,:)                          ! (n_mod x n_mod) block matrix for 1 normal point
            PetscScalar, allocatable :: V_int_0(:)                              ! interpolated V_0 at different normal positions
            PetscScalar, allocatable :: V_int_1(:)                              ! interpolated V_1 at different normal positions
            PetscScalar, allocatable :: V_int_2(:)                              ! interpolated V_2 at different normal positions
            PetscInt :: id, jd, kd                                              ! counters
            PetscInt :: k, m                                                    ! counters
            PetscInt :: i_min, i_max                                            ! absolute limits excluding the BC's
#if ldebug
            PetscScalar, allocatable :: loc_block_0_backup(:,:)                 ! backed up V_0 local block
#endif
            
            ! for tests
            PetscInt :: r_X_start, r_X_end                                      ! start block and end block
            
            ! initialize ierr
            ierr = 0
            
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
            if (grid_X%i_max.ne.r_X_end) then
                ierr = 1
                err_msg = 'end of matrix in this process does not coincide &
                    &with tabulated value'
                CHCKERR(err_msg)
            end if
            
            ! test whether numerical coefficients matrix has right size
            if (size(norm_disc_coeff).ne.2*norm_disc_prec_X+1) then
                ierr = 1
                err_msg = 'Array of discretization coefficients needs to have &
                    &the correct size 2*norm_disc_prec_X+1'
                CHCKERR(err_msg)
            end if
            
            ! set up i_min, depending on BC style
            select case (BC_style(1))
                case (1)
                    i_min = grid_X%i_min                                        ! will be overwritten
                case (2)
                    i_min = max(grid_X%i_min,1+norm_disc_prec_X)                ! first norm_disc_prec_X rows write left BC's
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
            
            ! set up i_max, depending on BC style
            select case (BC_style(2))
                case (1)
                    i_max = grid_X%i_max                                        ! will be overwritten
                case (2)
                    i_max = min(grid_X%i_max,grid_X%n(3)-norm_disc_prec_X)      ! last norm_disc_prec_X rows write right BC's
                case (3)
                    i_max = min(grid_X%i_max,grid_X%n(3)-1)                     ! last row writes right BC's
                case default
                    err_msg = 'No BC style associated with '//&
                        &trim(i2str(BC_style(1)))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! allocate interpolated V_int_i and local block
            allocate(V_int_0(X%n_mod**2))
            allocate(V_int_1(X%n_mod**2))
            allocate(V_int_2(X%n_mod**2))
            allocate(loc_block(X%n_mod,X%n_mod))
            
            ! iterate over all rows of this rank
            do kd = i_min-1,i_max-1                                             ! (indices start with 0 here)
                ! get interpolated terms in V_int_i
                call interp_V(V_0,loc_r_eq(kd+2-grid_X%i_min),V_int_0)
                call interp_V(V_1,loc_r_eq(kd+2-grid_X%i_min),V_int_1)
                call interp_V(V_2,loc_r_eq(kd+2-grid_X%i_min),V_int_2)
                
                ! fill the blocks
                
                ! -------------!
                ! BLOCKS ~ V_0 !
                ! -------------!
                ! fill local block
                do m = 1,X%n_mod
                    do k = 1,X%n_mod
                        loc_block(k,m) = con(V_int_0(c([k,m],.true.,X%n_mod)),&
                            &[k,m],.true.)                                      ! symmetric matrices need con()
                    end do
                end do

# if ldebug
                if (test_diff) then
                    ! back up the V_0 block
                    allocate(loc_block_0_backup(X%n_mod,X%n_mod))
                    loc_block_0_backup = loc_block
                end if
#endif

                ! add block to kd + (0,0)
                ierr = insert_block_mat(loc_block,mat,[kd,kd],n_X)
                CHCKERR('')
                
                ! -------------!
                ! BLOCKS ~ V_1 !
                ! -------------!
                ! fill local block
                do m = 1,X%n_mod
                    do k = 1,X%n_mod
                        loc_block(k,m) = V_int_1(c([k,m],.false.,X%n_mod))      ! asymetric matrices don't need con()
                    end do
                end do
                
                ! add block to kd + (0,-p..p) + Hermitian conjugate
                do jd = -norm_disc_prec_X,norm_disc_prec_X
                    ierr = insert_block_mat(loc_block*&
                        &norm_disc_coeff(jd+norm_disc_prec_X+1),mat,&
                        &[kd,kd+jd],n_X,transp=.true.)
                    CHCKERR('')
                end do
                
                ! -------------!
                ! BLOCKS ~ V_2 !
                ! -------------!
                ! fill local block
                do m = 1,X%n_mod
                    do k = 1,X%n_mod
                        loc_block(k,m) = con(V_int_2(c([k,m],.true.,X%n_mod)),&
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
                do jd = -norm_disc_prec_X,norm_disc_prec_X
                    do id = -norm_disc_prec_X,norm_disc_prec_X
                        ierr = insert_block_mat(loc_block*&
                            &norm_disc_coeff(id+norm_disc_prec_X+1)*&
                            &norm_disc_coeff(jd+norm_disc_prec_X+1),mat,&
                            &[kd+id,kd+jd],n_X)
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
            deallocate(V_int_0,V_int_1,V_int_2)
            deallocate(loc_block)
        end function fill_mat
        
        ! display information about matrix
        integer function disp_mat_inf(mat) result(ierr)
            character(*), parameter :: rout_name = 'disp_mat_inf'
            
            ! input / output
            Mat, intent(inout) :: mat                                           ! either A or B
            
            ! local variables
            real(dp) :: mat_info(MAT_INFO_SIZE)                                 ! information about matrix
            
            ! initialize ierr
            ierr = 0
            
            call MatGetInfo(mat,MAT_GLOBAL_SUM,mat_info,ierr)
            CHCKERR('')
            call writo('memory usage: '//&
                &trim(r2strt(mat_info(MAT_INFO_MEMORY)/1000))//' MB')
            call writo('nonzero''s allocated: '//&
                &trim(r2strt(mat_info(MAT_INFO_NZ_ALLOCATED))))
            call writo('nonzero''s used: '//&
                &trim(r2strt(mat_info(MAT_INFO_NZ_USED))))
            call writo('nonzero''s unused: '//&
                &trim(r2strt(mat_info(MAT_INFO_NZ_UNNEEDED))))
        end function disp_mat_inf
    end function setup_mats
    
    ! Sets the  boundary conditions:
    ! Deep  in the  plasma,  an  artificial Eigenvalue  EV_BC  is introduced  by
    ! setting the  diagonal components  of A  to EV_BC and  of B  to 1,  and the
    ! off-diagonal elements to zero.
    ! At the plasma surface, the surface energy is minimized as in [ADD REF].
    integer function set_BC(grid_X,X,A,B,loc_r_eq,i_geo,n_X,norm_disc_coeff) &
        &result(ierr)
        use num_vars, only: norm_disc_prec_X, BC_style
        use MPI_utilities, only: get_ser_var
        use utilities, only: con, c
        use grid_ops, only: trim_grid
        
        character(*), parameter :: rout_name = 'set_BC'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables
        Mat, intent(inout) :: A, B                                              ! Matrices A and B from A X = lambda B X
        PetscReal, intent(in) :: loc_r_eq(:)                                    ! unrounded index in tables V_int
        integer, intent(in) :: i_geo                                            ! at which geodesic index to perform the calculations
        integer, intent(in) :: n_X                                              ! number of grid points of perturbation grid
        PetscReal, intent(in) :: norm_disc_coeff(:)                             ! discretization coefficients for normal derivatives
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscInt :: n_min, n_max                                                ! absolute limits excluding the BC's
        
        ! local variables
        type(grid_type) :: grid_X_trim                                          ! trimmed perturbation grid
        PetscInt :: kd                                                          ! counter
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        call writo('Preparing variables')
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid_X,grid_X_trim)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! set up n_min, depending on BC style
        select case (BC_style(1))
            case (1)
                n_min = 2*norm_disc_prec_X                                      ! dirichlet BC requieres half stencil
            case (2)
                n_min = norm_disc_prec_X                                        ! mixed BC requires only one fourth of stencil
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
                n_max = 2*norm_disc_prec_X                                      ! dirichlet BC requieres half stencil
            case (2)
                n_max = norm_disc_prec_X                                        ! mixed BC requires only one fourth of stencil
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
            if (grid_X_trim%i_min.le.kd .and. grid_X_trim%i_max.ge.kd) then     ! decide which process does the setting
                select case (BC_style(1))
                    case (1)
                        ierr = set_BC_1(kd-1,A,B,.false.)                       ! indices start at 0
                        CHCKERR('')
                    case (2)
                        ierr = set_BC_2(kd-1,X,A,B,loc_r_eq(kd-grid_X%i_min+1),&
                            &i_geo,grid_X%n(3),norm_disc_coeff,.false.)         ! indices start at 0
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
        
        call writo('Setting BC at plasma edge')
        call lvl_ud(1)
        
        ! iterate over all positions where to set right BC
        do kd = grid_X%n(3),grid_X%n(3)-n_max+1,-1
            if (grid_X_trim%i_min.le.kd .and. grid_X_trim%i_max.ge.kd) then     ! decide which process does the setting
                select case (BC_style(2))
                    case (1)
                        ierr = set_BC_1(kd-1,A,B,.true.)                        ! indices start at 0
                        CHCKERR('')
                    case (2)
                        ierr = set_BC_2(kd-1,X,A,B,loc_r_eq(kd-grid_X%i_min+1),&
                            &i_geo,grid_X%n(3),norm_disc_coeff,.true.)          ! indices start at 0
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
        
        ! clean up
        call dealloc_grid(grid_X_trim)
        
        call lvl_ud(-1)
        
        call lvl_ud(-1)
    contains
        ! set BC style 1:
        !   A(ind,ind) = EV_BC, B(ind,ind) = 1,
        !   A(ind+1..ind+p,ind+1..ind+p) = 0, B(ind+1..ind+p,ind+1..ind+p) = 0,
        ! where ind indicates the row where the BC is centered.
        integer function set_BC_1(ind,A,B,BC_right) result(ierr)
            use num_vars, only: EV_BC, norm_disc_prec_X
            
            character(*), parameter :: rout_name = 'set_BC_1'
            
            ! input / output
            integer, intent(in) :: ind                                          ! position at which to set BC
            Mat, intent(inout) :: A, B                                          ! Matrices A and B from A X = lambda B X
            logical :: BC_right                                                 ! if BC is at right (so there are vacuum terms)
            
            ! local variables
            integer :: kd                                                       ! counter
            PetscScalar, allocatable :: loc_block(:,:)                          ! (n_mod,n_mod) block matrix for 1 normal point
            PetscInt :: pmone                                                   ! plus or minus one
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(ind+1))//&
                &': Eigenvector set to zero',persistent=.true.)
            
            ! initialize local blocks
            allocate(loc_block(X%n_mod,X%n_mod))
            loc_block = 0.0_dp
            do kd = 1,X%n_mod
                loc_block(kd,kd) = 1.0_dp
            end do
            
            ! set up plus minus one
            if (BC_right) then
                pmone = 1
            else
                pmone = -1
            end if
            
            ! set block (ind,ind) and Hermitian conjugate
            ierr = insert_block_mat(EV_BC*loc_block,A,[ind,ind],n_X,&
                &overwrite=.true.)
            CHCKERR('')
            ierr = insert_block_mat(loc_block,B,[ind,ind],n_X,overwrite=.true.)
            CHCKERR('')
            
            ! iterate over range 2
            do kd = 1,2*norm_disc_prec_X
                ! set block (ind,ind+/-kd) and Hermitian conjugate
                ierr = insert_block_mat(0*loc_block,A,[ind,ind-pmone*kd],n_X,&
                    &overwrite=.true.,transp=.true.)
                CHCKERR('')
                ierr = insert_block_mat(0*loc_block,B,[ind,ind-pmone*kd],n_X,&
                    &overwrite=.true.,transp=.true.)
                CHCKERR('')
            end do
        end function set_BC_1
        
        ! set BC style 2:
        !   minimization of surface energy term (see [ADD REF])
        integer function set_BC_2(ind,X,A,B,loc_r_eq,i_geo,n_X,norm_disc_coeff,&
            &BC_right) result(ierr)
            character(*), parameter :: rout_name = 'set_BC_2'
            
            ! input / output
            integer, intent(in) :: ind                                          ! position at which to set BC
            type(X_2_type), intent(in) :: X                                     ! tensorial perturbation variables
            Mat, intent(inout) :: A, B                                          ! Matrices A and B from A X = lambda B X
            PetscReal, intent(in) :: loc_r_eq                                   ! unrounded index in tables V_int
            integer, intent(in) :: i_geo                                        ! at which geodesic index to perform the calculations
            integer, intent(in) :: n_X                                          ! number of grid points of perturbation grid
            PetscReal, intent(in) :: norm_disc_coeff(:)                         ! discretization coefficients for normal derivatives
            logical :: BC_right                                                 ! if BC is at right (so there are vacuum terms)
            
            ! local variables
            PetscInt :: jd                                                      ! counter
            PetscScalar, allocatable :: PV_int_0(:)                             ! PV_int_0 in equilibrium normal grid
            PetscScalar, allocatable :: PV_int_1(:)                             ! PV_int_1 in equilibrium normal grid
            PetscScalar, allocatable :: PV_int_2(:)                             ! PV_int_2 in equilibrium normal grid
            PetscScalar, allocatable :: KV_int_0(:)                             ! KV_int_0 in equilibrium normal grid
            PetscScalar, allocatable :: KV_int_1(:)                             ! KV_int_1 in equilibrium normal grid
            PetscScalar, allocatable :: KV_int_2(:)                             ! KV_int_2 in equilibrium normal grid
            PetscScalar, allocatable :: V_int_0_mod(:,:,:)                      ! V_0 + (V_1+delta) V_2 (V_1+delta)^*T
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(ind+1))//&
                &': Minimization of surface energy',persistent=.true.)
            
            ! initialize PV_i and KV_i
            allocate(PV_int_0(X%n_mod**2))
            allocate(PV_int_1(X%n_mod**2))
            allocate(PV_int_2(X%n_mod**2))
            allocate(KV_int_0(X%n_mod**2))
            allocate(KV_int_1(X%n_mod**2))
            allocate(KV_int_2(X%n_mod**2))
           
            ! get interpolated terms in V_int_i
            call interp_V(X%PV_int_0(:,i_geo,:),loc_r_eq,PV_int_0)
            call interp_V(X%PV_int_1(:,i_geo,:),loc_r_eq,PV_int_1)
            call interp_V(X%PV_int_2(:,i_geo,:),loc_r_eq,PV_int_2)
            call interp_V(X%KV_int_0(:,i_geo,:),loc_r_eq,KV_int_0)
            call interp_V(X%KV_int_1(:,i_geo,:),loc_r_eq,KV_int_1)
            call interp_V(X%KV_int_2(:,i_geo,:),loc_r_eq,KV_int_2)
            
            ! -----------------!
            ! BLOCKS ~ V_0,mod !
            ! -----------------!
            ! calculate modified terms V_int_0_mod
            allocate(V_int_0_mod(X%n_mod,X%n_mod,2))
            ierr = calc_V_0_mod(PV_int_0,KV_int_0,PV_int_1,KV_int_1,KV_int_2,&
                &V_int_0_mod)
            CHCKERR('')
            
            ! add block to ind + (0,0)
            ierr = insert_block_mat(V_int_0_mod(:,:,1),A,[ind,ind],n_X)
            CHCKERR('')
            ierr = insert_block_mat(V_int_0_mod(:,:,2),B,[ind,ind],n_X)
            CHCKERR('')
            
            ! deallocate modified V_0
            deallocate(V_int_0_mod)
            
            ! -------------!
            ! BLOCKS ~ vac !
            ! -------------!
            ! add block to ind + (0,-p..p) + Hermitian conjugate
            if (BC_right) then
                do jd = -norm_disc_prec_X,norm_disc_prec_X
                    if (ind.lt.n_X .and. ind+jd.lt.n_X) then
                        ierr = insert_block_mat(-X%vac_res*&
                            &norm_disc_coeff(jd+norm_disc_prec_X+1),A,&
                            &[ind,ind+jd],n_X,transp=.true.)
                        CHCKERR('')
                    end if
                end do
            end if
        end function set_BC_2
        
        ! set BC style 3:
        !   minimization of vacuum energy (see [ADD REF])
        integer function set_BC_3(ind,X,A) result(ierr)
            character(*), parameter :: rout_name = 'set_BC_2'
            
            ! input / output
            integer, intent(in) :: ind                                          ! position at which to set BC
            type(X_2_type), intent(in) :: X                                     ! tensorial perturbation variables
            Mat, intent(inout) :: A                                             ! Matrices A from A X = lambda B X
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Boundary style at row '//trim(i2str(ind+1))//&
                &': Minimization of vacuum energy',persistent=.true.)
            
            ! -------------!
            ! BLOCKS ~ vac !
            ! -------------!
            ! add block to ind + (0,0)
            ierr = insert_block_mat(X%vac_res,A,[ind,ind],n_X)
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
            allocate(KV_2_inv(X%n_mod,X%n_mod))
            allocate(V_triple(X%n_mod,X%n_mod))
            allocate(KV_1_loc(X%n_mod,X%n_mod))
            allocate(PV_1_loc(X%n_mod,X%n_mod))
            
            ! make local copy of KV_2 into the inverse array
            KV_2_inv = 0._dp
            do m = 1,X%n_mod
                do k = 1,m                                                      ! only save upper diagonal part
                    KV_2_inv(k,m) = con(KV_2(c([k,m],.true.,X%n_mod)),&
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
            call zpotrf(uplo,X%n_mod,KV_2_inv,X%n_mod,ierr)
            CHCKERR('Failed to decompose KV_2')
            call zpotri(uplo,X%n_mod,KV_2_inv,X%n_mod,ierr)
            CHCKERR('Failed to invert KV_2')
            do k = 1,X%n_mod
                do m = 1,k-1
                    KV_2_inv(k,m) = conjg(KV_2_inv(m,k))
                end do
            end do
            
            ! save local copy of KV_1 and PV_1
            do m = 1,X%n_mod
                do k = 1,X%n_mod
                    KV_1_loc(k,m) = KV_1(c([k,m],.false.,X%n_mod))              ! asymetric matrices don't need con()
                    PV_1_loc(k,m) = PV_1(c([k,m],.false.,X%n_mod))              ! asymetric matrices don't need con()
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
            do m = 1,X%n_mod
                do k = 1,X%n_mod
                    V_0_mod(k,m,1) = &
                        &con(PV_0(c([k,m],.true.,X%n_mod)),[k,m],.true.)        ! symmetric matrices need con()
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
            do m = 1,X%n_mod
                do k = 1,X%n_mod
                    V_0_mod(k,m,2) = &
                        &con(KV_0(c([k,m],.true.,X%n_mod)),[k,m],.true.)        ! symmetric matrices need con()
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
    integer function setup_solver(grid_X,X,A,B,solver) result(ierr)
        use num_vars, only: n_sol_requested, tol_slepc, max_it_slepc
        
        character(*), parameter :: rout_name = 'setup_solver'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_2_type), intent(in) :: X                                         ! tensorial perturbation variables
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
        CHCKERR('EPSSetDimensions failed')
        call EPSSetTolerances(solver,tol_slepc,max_it_slepc,ierr)
        CHCKERR('EPSSetTolerances failed')
        
        ! user output
        call writo('tolerance: '//trim(r2str(tol_slepc)))
        call writo('maximum nr. of iterations: '//trim(i2str(max_it_slepc)))
        
        call lvl_ud(-1)
    end function setup_solver
    
    ! sets up guess in solver
    subroutine setup_guess(sol,A,solver,start_id,prev_n_EV)
        ! input / output
        type(sol_type), intent(in) :: sol                                       ! solution variables
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
        
        ! set guess for EV if sol vec is allocated and use_guess is true
        if (allocated(sol%vec)) then
            ! allocate guess vectors
            allocate(guess_vec(prev_n_EV))
            
            ! create vecctor guess_vec
            do kd = 1,prev_n_EV
                !call MatGetVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,istat)        ! get compatible parallel vectors to matrix A (petsc 3.5.3)
                call MatCreateVecs(A,guess_vec(kd),PETSC_NULL_OBJECT,istat)     ! get compatible parallel vectors to matrix A (petsc 3.6.1)
                call VecSetOption(guess_vec(kd),&
                    &VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE,istat)              ! ignore negative indices for compacter notation
            end do
            
            ! set the indices and nr. of values to set in guess_vec
            n_prev_guess = size(sol%vec,2)
            allocate(guess_id(n_prev_guess*sol%n_mod))
            do id = 0,n_prev_guess-1
                do jd = 1,sol%n_mod
                    guess_id(id*sol%n_mod+jd) = (start_id-1)*sol%n_mod*2 + &
                        &id*sol%n_mod*2 + jd
                end do
            end do
            
            do kd = 1,prev_n_EV
                ! set even values of guess_vec
                call VecSetValues(guess_vec(kd),n_prev_guess*sol%n_mod,&
                    &guess_id-1,sol%vec(:,:,kd),INSERT_VALUES,istat)
                ! set odd values of guess_vec
                call VecSetValues(guess_vec(kd),n_prev_guess*sol%n_mod,&
                    &guess_id-1-sol%n_mod,sol%vec(:,:,kd),INSERT_VALUES,istat)
                    
                ! assemble the guess vector
                call VecAssemblyBegin(guess_vec(kd),istat)
                call VecAssemblyEnd(guess_vec(kd),istat)
            end do
                
            ! set guess
            call EPSSetInitialSpace(solver,prev_n_EV,guess_vec,istat)
            
            !! visualize guess
            !if (allocated(sol%vec) .and. id.le.prev_n_EV) then
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
        call writo('number of iterations: '//trim(i2str(n_it)))
        call writo('number of converged solutions: '//trim(i2str(n_conv)))
        call writo('number of requested solutions: '//trim(i2str(n_ev)))
        call writo('maximum dimension of the subspace to be used by solver: '//&
            &trim(i2str(ncv)))
        call writo('maximum dimension allowed for projected problem : '//&
            &trim(i2str(mpd)))
        
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
#if ldebug
    integer function store_results(grid_X,sol,solver,max_n_EV,A,B) result(ierr)
#else
    integer function store_results(grid_X,sol,solver,max_n_EV) result(ierr)
#endif
        use eq_vars, only: T_0
        use num_vars, only: use_normalization, EV_BC, prog_name, &
            &output_name, rich_lvl_nr, max_it_r, n_procs, rank
        use files_utilities, only: nextunit
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'store_results'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(sol_type), intent(inout) :: sol                                    ! solution variables
        EPS, intent(inout) :: solver                                            ! EV solver
        PetscInt, intent(inout) :: max_n_EV                                     ! nr. of EV's saved, up to n_conv
#if ldebug
        Mat, intent(inout) :: A, B                                              ! matrix A and B
#endif
        
        ! local variables
        PetscInt :: id, id_tot                                                  ! counters
        Vec :: sol_vec                                                          ! solution EV parallel vector
        PetscInt :: one = 1                                                     ! one
        PetscInt, allocatable :: X_vec_max_loc(:)                               ! location of X_vec_max of this process
        PetscReal :: error                                                      ! error of EPS solver
        PetscReal :: tol_complex = 1.E-5_dp                                     ! tolerance for an EV to be considered complex
        PetscReal :: tol_EV_BC = 1.E-3_dp                                       ! tolerance for an EV to be considered due to the BC
        PetscReal, parameter :: infinity = 1.E40_dp                             ! beyond this value, modes are not saved
        PetscScalar :: X_val_loc                                                ! local X val
        PetscScalar, allocatable :: X_vec_max(:)                                ! max of X_vec of all processes, including phase at max
        character(len=max_str_ln) :: full_output_name                           ! full name
        character(len=max_str_ln) :: format_val                                 ! format
        character(len=max_str_ln) :: format_head                                ! header
        character(len=max_str_ln) :: EV_err_str                                 ! String with information about EV error
        integer :: output_EV_i                                                  ! file number
        integer :: n_digits                                                     ! nr. of digits for the integer number
        integer :: n_err(3)                                                     ! how many errors there were
#if ldebug
        character(len=max_str_ln) :: err_msg                                    ! error message
        Mat :: err_mat                                                          ! A - omega^2 B
        Vec :: err_vec                                                          ! Ax - omega^2 BX
        PetscReal :: err_norm                                                   ! norm of error
        PetscReal :: error_est                                                  ! error estimate of EPS solver
#endif
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! create solution vector
        call VecCreateMPIWithArray(PETSC_COMM_WORLD,one,&
            &grid_X%loc_n_r*sol%n_mod,grid_X%n(3)*sol%n_mod,PETSC_NULL_SCALAR,&
            &sol_vec,ierr)
        CHCKERR('Failed to create MPI vector with arrays')
        
        write(*,*) '!!!!!!!!!!!!!!!!!!! ALLOCATE VEC AND VAL USING THE CREATE_SOL ROUTINE !!!!!!!!!!'
        ! allocate vec
        if (allocated(sol%vec)) deallocate(sol%vec)
        allocate(sol%vec(1:sol%n_mod,1:grid_X%loc_n_r,1:max_n_EV))
        sol%vec = 0.0_dp
        
        ! allocate val
        if (allocated(sol%val)) deallocate(sol%val)
        allocate(sol%val(1:max_n_EV))
        sol%val = 0.0_dp
        
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
            full_output_name = prog_name//'_'//output_name//'_EV'
            if (max_it_r.gt.1) full_output_name = &
                &trim(full_output_name)//'_R'//trim(i2str(rich_lvl_nr))         ! append Richardson level
            full_output_name = trim(full_output_name)//'.txt'
            open(unit=nextunit(output_EV_i),file=full_output_name,iostat=ierr)
            CHCKERR('Cannot open EV output file')
        end if
        
        ! output to file
        if (rank.eq.0) then
            if (use_normalization) then
                write(output_EV_i,'(A)') '# Eigenvalues normalized to the &
                    &squared  Alfven frequency omega_A^2 = '
                write(output_EV_i,'(A)') '#     ('//trim(r2str(1._dp/T_0))//&
                    &' Hz)^2 = '//trim(r2str(1._dp/T_0**2))//' Hz^2'
            else
                write(output_EV_i,'(A)') '# Eigenvalues'
            end if
            write(output_EV_i,format_head) '#  I                            ', &
                &'real part               ', 'imaginary part          ', &
                &'relative precision      '
        end if
        
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
#if ldebug
            if (debug_store_results) then
                call EPSGetErrorEstimate(solver,id_tot-1,error_est)             ! get error estimate
                CHCKERR('EPSGetErrorEstimate failed')
            end if
#endif
            
            ! set up local solution val
            X_val_loc = sol%val(id)
            
            ! tests
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
                write(output_EV_i,format_val) id,realpart(X_val_loc),&
                    &imagpart(X_val_loc),error
                
                ! if error, print explanation
                if (EV_err_str.ne.'') write(output_EV_i,'(A)') trim(EV_err_str)
            end if
            
            !! visualize solution
            !call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,&
                !&'solution '//trim(i2str(id))//' vector with '//&
                !&trim(i2str(grid_X%n(3)))//' points',0,&
                !&0,1000,1000,guess_viewer,ierr)
            !call VecView(sol_vec,guess_viewer,ierr)
            !CHCKERR('Cannot view vector')
            
            !call writo('Eigenvector = ')
            !call VecView(sol_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
            !CHCKERR('Cannot view vector')
            
#if ldebug
            if (debug_store_results) then
                ! user message
                call writo('Testing whether A x - omega^2 B x = 0 for EV '//&
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
                
                ! calculate absolute norm
                call VecNorm(err_vec,NORM_2,err_norm,ierr)
                CHCKERR('Failed to calculate norm')
                
                ! get relative norm
                err_norm = err_norm/abs(sol%val(id))
                
                ! visualize solution
                call writo('error: '//trim(r2str(err_norm))//', given :'//&
                    &trim(r2str(error))//', estimate: '//trim(r2str(error_est)))
                !call VecView(err_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
                !CHCKERR('Cannot view vector')
                
                call lvl_ud(-1)
                
                ! clean up
                call VecDestroy(err_vec,ierr)                                   ! destroy error vector
                CHCKERR('Failed to destroy err_vec')
            end if
#endif
            
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
        
        ! close output file if master
        if (rank.eq.0) close(output_EV_i)
        
        ! normalize Eigenvectors to make output more easily comparable
        allocate(X_vec_max(n_procs))
        allocate(X_vec_max_loc(2))
        do id = 1,size(sol%vec,3)
            ! find local maximum
            X_vec_max_loc = maxloc(abs(sol%vec(:,:,id)))
            ierr = get_ser_var([sol%vec(X_vec_max_loc(1),X_vec_max_loc(2),id)],&
                &X_vec_max,scatter=.true.)
            CHCKERR('')
            ! find global maximum
            X_vec_max_loc(1) = maxloc(abs(X_vec_max),1)
            sol%vec(:,:,id) = sol%vec(:,:,id) / X_vec_max(X_vec_max_loc(1))
        end do
        deallocate(X_vec_max_loc)
        deallocate(X_vec_max)
        
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
        
        call writo('basic statistics:')
        call lvl_ud(1)
        
        call writo('min: '//trim(c2strt(sol%val(1))))
        call writo('max: '//trim(c2strt(sol%val(max_n_EV))))
        
        call lvl_ud(-1)
        
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
            PetscScalar, allocatable :: X_val_loc(:)                            ! local copy of X_val
            PetscScalar, allocatable :: X_vec_loc(:,:,:)                        ! local copy of X_vec
            PetscBool :: remove_next_loc = .false.                              ! local copy of remove_next
            
            ! only remove if faulty solutions are not optionally retained
            if (retain_all_sol) then
                id = id+1                                                       ! increment the solution
            else
                ! set up local remove_next
                if (present(remove_next)) remove_next_loc = remove_next
                
                ! save old arrays
                allocate(X_val_loc(max_id))
                allocate(X_vec_loc(sol%n_mod,grid_X%loc_n_r,max_id))
                X_val_loc = sol%val
                X_vec_loc = sol%vec
                
                ! reallocate arrays
                deallocate(sol%val,sol%vec)
                if (remove_next_loc) then
                    allocate(sol%val(id-1))
                    allocate(sol%vec(sol%n_mod,grid_X%loc_n_r,id-1))
                else
                    allocate(sol%val(max_id-1))
                    allocate(sol%vec(sol%n_mod,grid_X%loc_n_r,max_id-1))
                end if
                
                ! copy other values
                sol%val(1:id-1) = X_val_loc(1:id-1)
                sol%vec(:,:,1:id-1) = X_vec_loc(:,:,1:id-1)
                if (.not.remove_next_loc) then
                    sol%val(id:max_id-1) = X_val_loc(id+1:max_id)
                    sol%vec(:,:,id:max_id-1) = X_vec_loc(:,:,id+1:max_id)
                end if
                
                ! adapt max_id
                if (remove_next_loc) then
                    max_id = id-1
                else
                    max_id = max_id-1
                end if
                
                ! clean up
                deallocate(X_val_loc,X_vec_loc)
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
    
    ! Interpolates an array V at position x.
    ! (also  correct  if i_hi = i_lo)
    subroutine interp_V(V,x,V_int)
        ! input / output
        PetscScalar, intent(in) :: V(:,:)                                       ! V
        PetscReal, intent(in) :: x                                              ! coordinate at which to interpolate
        PetscScalar, intent(inout) :: V_int(:)                                  ! interpolated V
        
        ! local variables
        PetscInt :: i_lo, i_hi                                                  ! upper and lower index
        
        ! set up i_lo and i_hi
        i_lo = floor(x)
        i_hi = ceiling(x)
        
        V_int = V(:,i_lo) + (V(:,i_hi)-V(:,i_lo))*(x-i_lo)                      ! because i_hi - i_lo = 1
    end subroutine interp_V
    
    ! Insert  a block pertaining  to 1 normal  point into a  matrix. Optionally,
    ! also set the Hermitian transpose and / or overwrite instead of add value.
    integer function insert_block_mat(block,mat,ind,n_X,transp,overwrite) &
        &result(ierr)
        use MPI_utilities, only: wait_MPI
#if ldebug
        use num_vars, only: rank
#endif
        
        character(*), parameter :: rout_name = 'insert_block_mat'
        
        ! input / output
        PetscScalar :: block(:,:)                                               ! (n_mod x n_mod) block matrix for 1 normal point
        Mat, intent(inout) :: mat                                               ! matrix in which to insert block
        PetscInt, intent(in) :: ind(2)                                          ! 2D index in matrix
        PetscInt, intent(in) :: n_X                                             ! number of grid points of perturbation grid
        PetscBool, intent(in), optional :: transp                               ! also set Hermitian transpose
        PetscBool, intent(in), optional :: overwrite                            ! overwrite
        
        ! local variables
        PetscInt, allocatable :: loc_k(:), loc_m(:)                             ! the locations at which to add the blocks to the matrices
        PetscInt :: kd                                                          ! counter
        PetscBool :: transp_loc                                                 ! local copy of transp
        PetscBool :: overwrite_loc                                              ! local copy of overwrite
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscInt :: operation                                                   ! either ADD_VALUES or INSERT_VALUES
        
        ! initialize ierr
        ierr = 0
        
        ! set up local transp and overwrite
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        overwrite_loc = .false.
        if (present(overwrite)) overwrite_loc = overwrite
        
        ! set up local k and m
        allocate(loc_k(size(block,1)))
        allocate(loc_m(size(block,2)))
        loc_k = [(kd, kd = 0,size(block,1)-1)] + ind(1)*size(block,1)
        loc_m = [(kd, kd = 0,size(block,2)-1)] + ind(2)*size(block,2)
        
        ! set operation
        if (overwrite_loc) then
            operation = INSERT_VALUES
        else
            operation = ADD_VALUES
        end if
        
#if ldebug
        ! user output
        if (debug_insert_block_mat) then
            call writo('at (k,m) = ')
            write(*,*) loc_k
            write(*,*) loc_m
        end if
#endif
        
        ! only set values if within matrix range
        if (minval(ind).ge.0 .and. maxval(ind).lt.n_X) then
#if ldebug
            if (debug_insert_block_mat) then
                call writo('following local block is going to be added: ')
                write(*,*) 'Re ='
                call print_ar_2(realpart(block))
                write(*,*) 'Im ='
                call print_ar_2(imagpart(block))
                if (overwrite_loc) then
                    call writo('(with operation INSERT_VALUES)')
                else
                    call writo('(with operation INSERT_VALUES)')
                end if
            endif
#endif
            
            ! set error message
            err_msg = 'Couldn''t add values to matrix'
            
            ! set values
            call MatSetValues(mat,size(block,1),loc_k,size(block,2),loc_m,&
                &block,operation,ierr)
            CHCKERR(err_msg)
            
            ! set transpose if wanted
            if (transp_loc) then
                call MatSetValues(mat,size(block,2),loc_m,size(block,1),loc_k,&
                    &conjg(transpose(block)),operation,ierr)
                CHCKERR(err_msg)
            end if
        else
#if ldebug
            if (debug_insert_block_mat) &
                &call writo('Due to out of range no local block is added')
#endif
        end if
        
#if ldebug
            if (debug_insert_block_mat) then
                if (rank.eq.0) read(*,*)
                ierr = wait_MPI()
                CHCKERR('')
            endif
#endif
    end function insert_block_mat
end module
