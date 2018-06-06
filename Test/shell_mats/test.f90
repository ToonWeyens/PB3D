program shell_mats
#include <slepc/finclude/slepceps.h>
    use slepceps
    
    implicit none
    
    ! variables
    PetscInt :: ierr
    Mat :: A, B
    PetscInt :: n
    EPS :: solver
    PetscInt :: rank
    PetscInt :: its, nev, maxit
    PetscReal :: tol
    
    PETSC_COMM_WORLD = MPI_Comm_world
    
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    
    n = 10
    call MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,PETSC_NULL_OBJECT,A,ierr)
    call MatShellSetOperation(A,MATOP_MULT,shell_mat_mult_left,ierr)
    call MatShellSetOperation(A,MATOP_GET_DIAGONAL,shell_mat_get_diagonal_left,&
        &ierr)
    
    call MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,PETSC_NULL_OBJECT,B,ierr)
    call MatShellSetOperation(B,MATOP_MULT,shell_mat_mult_right,ierr)
    call MatShellSetOperation(B,MATOP_GET_DIAGONAL,&
        &shell_mat_get_diagonal_right,ierr)
    
    call EPSCreate(PETSC_COMM_WORLD,solver,ierr)
    !call EPSSetOperators(solver,A,PETSC_NULL_OBJECT,ierr)
    call EPSSetOperators(solver,A,B,ierr)
    call EPSSetFromOptions(solver,ierr)
    call EPSSolve(solver,ierr) 
    
    call EPSGetIterationNumber(solver,its,ierr)
    call EPSGetDimensions(solver,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
    call EPSGetTolerances(solver,tol,maxit,ierr)
    if (rank .eq. 0) then
        write(*,'(A,I4)') ' Number of iterations of the method: ', its
        write(*,'(A,I2)') ' Number of requested eigenvalues:', nev
        write(*,'(A,1PE10.4,A,I6)') ' Stopping condition: tol=', tol,', &
            &maxit=', maxit
    endif
    
    call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,&
        &PETSC_VIEWER_ASCII_INFO_DETAIL,ierr)
    call EPSReasonView(solver,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call EPSErrorView(solver,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr)
    call EPSDestroy(solver,ierr)
    call MatDestroy(A,ierr)
    
    call SlepcFinalize(ierr)
contains
    subroutine shell_mat_mult_left(mat,X,Y,ierr)
        ! input / output
        Mat :: mat
        Vec :: X, Y
        PetscInt :: ierr
        
        ! local variables
        PetscInt :: kd
        PetscScalar, pointer :: X_loc(:)                                        ! local pointer to X
        PetscScalar, pointer :: Y_loc(:)                                        ! local pointer to Y
        
        ! get vectors
        call VecGetArrayReadF90(X,X_loc,ierr)
        call VecGetArrayF90(Y,Y_loc,ierr)
        
        ! multiply
        do kd = 1,size(X_loc)
            Y_loc(kd) = kd*X_loc(kd)
        end do
        
        ! restore vectors
        call VecRestoreArrayReadF90(X,X_loc,ierr)
        call VecRestoreArrayF90(Y,Y_loc,ierr)
    end subroutine shell_mat_mult_left
    
    subroutine shell_mat_mult_right(mat,X,Y,ierr)
        ! input / output
        Mat :: mat
        Vec :: X, Y
        PetscInt :: ierr
        
        ! local variables
        PetscInt :: kd
        PetscScalar, pointer :: X_loc(:)                                        ! local pointer to X
        PetscScalar, pointer :: Y_loc(:)                                        ! local pointer to Y
        
        ! get vectors
        call VecGetArrayReadF90(X,X_loc,ierr)
        call VecGetArrayF90(Y,Y_loc,ierr)
        
        ! multiply
        do kd = 1,size(X_loc)
            Y_loc(kd) = kd*X_loc(kd)
        end do
        
        ! restore vectors
        call VecRestoreArrayReadF90(X,X_loc,ierr)
        call VecRestoreArrayF90(Y,Y_loc,ierr)
    end subroutine shell_mat_mult_right
    
    subroutine shell_mat_get_diagonal_left(mat,diag,ierr)
        ! input / output
        Mat :: mat
        Vec :: diag
        PetscInt :: ierr
        
        ! local variables
        PetscScalar, pointer :: diag_loc(:)
        PetscInt :: kd
        
        ! get vectors
        call VecGetArrayF90(diag,diag_loc,ierr)
        
        ! multiply
        do kd = 1,size(diag_loc)
            diag_loc(kd) = kd
        end do
        
        ! restore vectors
        call VecRestoreArrayF90(diag,diag_loc,ierr)
    end subroutine shell_mat_get_diagonal_left
    
    subroutine shell_mat_get_diagonal_right(mat,diag,ierr)
        ! input / output
        Mat :: mat
        Vec :: diag
        PetscInt :: ierr
        
        ! local variables
        PetscScalar, pointer :: diag_loc(:)
        PetscInt :: kd
        
        ! get vectors
        call VecGetArrayF90(diag,diag_loc,ierr)
        
        ! multiply
        do kd = 1,size(diag_loc)
            diag_loc(kd) = kd
        end do
        
        ! restore vectors
        call VecRestoreArrayF90(diag,diag_loc,ierr)
    end subroutine shell_mat_get_diagonal_right
end program
