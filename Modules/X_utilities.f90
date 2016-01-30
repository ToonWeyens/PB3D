!------------------------------------------------------------------------------!
!   Numerical utilities related to perturbation operations                     !
!------------------------------------------------------------------------------!
module X_utilities
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_name_ln, iu
    use grid_vars, only: grid_type
    use X_vars, only: n_mod_X
    
    implicit none
    
    private
    public get_suffix, is_necessary_X
    
    ! interfaces
    interface get_suffix
        module procedure get_suffix_1, get_suffix_2
    end interface
    
contains
    ! Sets the suffix used to refer to a perturbation quantity.
    function get_suffix_1(id,lim_sec_X) result(res)                             ! vectorial version
        ! input / output
        integer, intent(in) :: id                                               ! mode index
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer :: res                                                          ! output
        
        ! local variables
        integer :: lim_sec_X_loc(2)                                             ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set suffix
        res = lim_sec_X_loc(1)-1+id
    end function get_suffix_1
    function get_suffix_2(id,jd,lim_sec_X) result(res)                          ! tensorial version
        ! input / output
        integer, intent(in) :: id, jd                                           ! mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        integer :: res(2)                                                       ! output
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set suffix
        res = [lim_sec_X_loc(1,1)-1+id,lim_sec_X_loc(1,2)-1+jd]
    end function get_suffix_2
    
    ! Determines whether a variable needs to be  considered: Only if it is on or
    ! below the diagonal for symmetric quantities.
    logical function is_necessary_X(sym,sec_X_id,lim_sec_X) result(res)
        ! input / output
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in) :: sec_X_id(2)                                      ! mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local version of lim_sec_X
        
        ! initialize res
        res = .true.
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! modify res depending on symmetry
        if (sym) then
            if (lim_sec_X_loc(1,1)+sec_X_id(1).lt.&
                &lim_sec_X_loc(1,2)+sec_X_id(2)) res = .false.
        end if
    end function is_necessary_X
end module X_utilities
