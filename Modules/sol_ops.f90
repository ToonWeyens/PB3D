!------------------------------------------------------------------------------!
!   Operations  on the  solution  vectors  such as  returning  X,  U or  their !
!   derivates, etc.                                                            !
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    !use messages, only: lvl_ud, writo, print_ar_2
    !use output_ops, only: print_GP_2D, print_GP_3D
    !use str_ops, only: i2str, r2strt, r2str

    implicit none
    private
    public calc_real_X, calc_real_U

contains
    ! calculates  the  normal   component  of  the  perturbation,   or  a  first
    ! derivative.  The   input  is   given  for   a  range   (r,theta,zeta)_F in
    ! Flux coordinates
    integer function calc_real_X(r_F,theta_F,zeta_F,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_real_X'
        
        ! input / output
        real(dp), intent(in) :: r_F(:), theta_F(:), zeta_F(:)                   ! range of normal coord, theta and zeta
        integer, intent(in), optional :: deriv(3)                               ! optional derivatives in angular coordinates
        
        ! local variables
        integer :: n_pts                                                        ! nr. of points
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_pts
        n_pts = size(r_F)
        
        ! tests
        if (size(theta_F).ne.n_pts .or. size(zeta_F).ne.n_pts) then
            ierr = 1
            err_msg = 'The coordinate arrays need to have the same size'
            CHCKERR(err_msg)
        end if
        if (maxval(deriv).gt.1) then
            ierr = 1
            err_msg = 'Maximum derivative is of order 1'
            CHCKERR(err_msg)
        end if
    end function calc_real_X
    
    ! calculates  the  geodesic  component  of  the  perturbation,  or  a  first
    ! derivative.  The   input  is   given  for   a  range   (r,theta,zeta),  or
    ! alternatively, (r,theta,alpha) or (r,zeta,alpha)  (depending on which flux
    ! is used as normal coordinate), selected by use_alpha
    integer function calc_real_U() result(ierr)
        ! initialize ierr
        ierr = 0
    end function calc_real_U
end module sol_ops
