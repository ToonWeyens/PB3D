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
    integer function calc_real_X(r_F,theta_F,zeta_F,time,X_F,deriv) result(ierr)
        character(*), parameter :: rout_name = 'calc_real_X'
        
        ! input / output
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        real(dp), intent(in) :: time(:)                                         ! time range
        complex(dp), intent(inout) :: X_F(:,:,:)                                ! normal component of perturbation
        integer, intent(in), optional :: deriv(3)                               ! optional derivatives in angular coordinates
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(theta_F,1)
        n_zeta = size(theta_F,2)
        n_r = size(theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(zeta_F,1) .or. n_zeta.ne.size(zeta_F,2) .or. &
            &n_r.ne.size(zeta_F,3) .or. n_r.ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(X_F,1) .or. n_zeta.ne.size(X_F,2) .or. &
            &n_r.ne.size(X_F,3)) then
            ierr = 1
            err_msg = 'X_F needs to have the correct dimensions'
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
