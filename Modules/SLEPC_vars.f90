!------------------------------------------------------------------------------!
!   Variables pertaining to SLEPC (and PETSC) operations                       !
!------------------------------------------------------------------------------!
! References:                                                                  !
! [1]   http://slepc.upv.es/documentation/current/docs/manualpages/EPS/        !
!       EPSComputeRelativeError.html                                           !
!------------------------------------------------------------------------------!
module SLEPC_vars
#include <PB3D_macros.h>
#include <wrappers.h>
! for slepc 3.6.0:
#include <slepc/finclude/slepcepsdef.h>
! for slepc 3.5.3:
!#include <finclude/slepcepsdef.h>
    use slepceps
    use num_vars, only: dp

    implicit none
    private
    public shell_ctx, MatShellGetContext, init_shell_ctx, dealloc_shell_ctx
    
    ! shell matrix type context
    type shell_ctx
        PetscScalar, pointer :: V_0(:,:,:,:)                                    ! integrated tensorial variables of order 0
        PetscScalar, pointer :: V_1(:,:,:,:)                                    ! integrated tensorial variables of order 1
        PetscScalar, pointer :: V_2(:,:,:,:)                                    ! integrated tensorial variables of order 2
        PetscScalar :: sigma                                                    ! shift of diagonal
    end type shell_ctx
    
    ! explicit interface  to Petsc routine so that custom  type shell_ctx can be
    ! used (see http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/
    ! MatShellGetContext.html, as well as http://lists.mcs.anl.gov/pipermail/
    ! petsc-users/2015-October/027236.html)
    Interface MatShellGetContext                                                
        Subroutine MatShellGetContext(M,ctx,ierr)
            Import :: shell_ctx
            Mat :: M
            type(shell_ctx), pointer :: ctx
            PetscErrorCode :: ierr
        End Subroutine
    End Interface MatShellGetContext   
contains
    ! Initializes a shell matrix context.
    subroutine init_shell_ctx(M_ctx,V_0,V_1,V_2)
        type(shell_ctx), intent(inout), pointer :: M_ctx                        ! context for shell matrix M
        PetscScalar, intent(in), allocatable, target :: V_0(:,:,:,:)            ! integrated tensorial variables of order 0
        PetscScalar, intent(in), allocatable, target :: V_1(:,:,:,:)            ! integrated tensorial variables of order 1
        PetscScalar, intent(in), allocatable, target :: V_2(:,:,:,:)            ! integrated tensorial variables of order 2
        
        allocate(M_ctx)
        M_ctx%V_0 => V_0
        M_ctx%V_1 => V_1
        M_ctx%V_2 => V_2
        M_ctx%sigma = 0._dp
    end subroutine init_shell_ctx
    
    ! deallocates a shell matrix context
    subroutine dealloc_shell_ctx(M_ctx)
        type(shell_ctx), intent(inout), pointer :: M_ctx                        ! context for shell matrix M
        
        nullify(M_ctx%V_0)
        nullify(M_ctx%V_1)
        nullify(M_ctx%V_2)
        deallocate(M_ctx)
        nullify(M_ctx)
    end subroutine dealloc_shell_ctx
end module SLEPC_vars
