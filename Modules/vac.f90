!------------------------------------------------------------------------------!
!   Operations and variables pertaining to the vacuum response                 !
!------------------------------------------------------------------------------!
module vac
#include <PB3D_macros.h>
    use str_ops
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type

    implicit none
    private
    public calc_vac

contains
    ! calculates the vacuum response according to [ADD REF].
    integer function calc_vac(X) result(ierr)
        ! input / output
        type(X_type), intent(inout) :: X                                        ! perturbation variables
        
        ! initialize ierr
        ierr = 0
        
        ! user message
        call writo('Calculating vacuum response')
        call lvl_ud(1)
        
        call writo('NOT YET IMPLEMENTED: SETTING TO ZERO')
        X%vac_res = 0._dp
        
        call lvl_ud(-1)
        call writo('Vacuum response calculated')
    end function calc_vac
end module vac
