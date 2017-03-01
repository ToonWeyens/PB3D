!------------------------------------------------------------------------------!
!   Operations and variables pertaining to the vacuum response                 !
!------------------------------------------------------------------------------!
module vac
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use X_vars, only: X_2_type

    implicit none
    private
    public calc_vac

contains
    ! calculates the vacuum response according to [ADD REF].
    integer function calc_vac(X) result(ierr)
        ! input / output
        type(X_2_type), intent(inout) :: X                                      ! tensorial perturbation variables
        
        ! initialize ierr
        ierr = 0
        
        ! user message
        call writo('Calculate vacuum response')
        call lvl_ud(1)
        
        X%vac_res = 0._dp
        
        call lvl_ud(-1)
    end function calc_vac
end module vac
