!------------------------------------------------------------------------------!
!   Operations that use slepc (and petsc) routines
!------------------------------------------------------------------------------!
module slepc_ops
#include <PB3D_macros.h>
#include <finclude/slepcepsdef.h>
!#include <finclude/petscsys.h>
    use slepceps
    use num_vars, only: iu, dp, max_str_ln
    use message_ops, only: writo, print_ar_2
    use str_ops, only: r2strt, r2str, i2str

    implicit none
    private
    public solve_EV_system_slepc
    
contains
    ! This subroutine sets up  the matrices A ad B of  the generalized EV system
    ! described in [ADD REF] and solves them using the slepc suite
    integer function solve_EV_system_slepc(use_guess) result(ierr)
        use slepc_vars, only: start_slepc, stop_slepc, setup_matrices, &
            &setup_solver, setup_guess, get_solution, summarize_solution, &
            &store_results
        
        character(*), parameter :: rout_name = 'solve_EV_system_slepc'
        
        ! input / output
        PetscBool, intent(in) :: use_guess                                      ! whether to use a guess or not
        
        ! local variables
        Mat :: A                                                                ! matrix A in EV problem A X = lambda B X
        Mat :: B                                                                ! matrix B in EV problem A X = lambda B X
        EPS :: solver                                                           ! EV solver
        PetscInt :: max_n_EV                                                    ! nr. of EV's saved
        PetscInt, save :: guess_start_id = -10                                  ! start of index of previous vector, saved for next iteration
        PetscInt, save :: prev_n_EV                                             ! nr. of solutions of previous vector
        
        ! initialize ierr
        ierr = 0
        
        ! start slepc
        ierr = start_slepc()
        CHCKERR('')
        
        ! set up the matrix
        call writo('set up matrices...')
        
        ierr = setup_matrices(A,B)
        CHCKERR('')
        
        ! set up solver
        call writo('set up EV solver...')
        
        ierr = setup_solver(A,B,solver)
        CHCKERR('')
        
        ! set up guess
        call writo('set up guess...')
        
        if (use_guess) call setup_guess(A,solver,guess_start_id,prev_n_EV)
        
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
        
        ierr = store_results(solver,max_n_EV)
        CHCKERR('')
        
        ! finalize
        call writo('finalize slepc...')
        
        ierr = stop_slepc(A,B,solver,guess_start_id,prev_n_EV,max_n_EV)
        CHCKERR('')
    end function solve_EV_system_slepc
end module
