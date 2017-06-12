!------------------------------------------------------------------------------!
!   Splines, simplified interface to bspline-fortran from                      !
!   https://github.com/jacobwilliams/bspline-fortran.git                       !
!------------------------------------------------------------------------------!
module splines
#include <PB3D_macros.h>
#include <wrappers.h>
    use num_vars, only: dp, iu, max_str_ln
    use str_utilities
    
    implicit none
    private
    public spline3

    interface spline3
        module procedure spline3_real, spline3_complex
    end interface
    
contains
    ! This procedure  makes use of the  bspline library, but makes  it easier to
    ! use for 1-D applications where speed is not the main priority.
    integer function spline3_real(ord,x,y,xnew,ynew,dynew,d2ynew,extrap) &
        &result(ierr)                                                           ! real version
        use bspline_sub_module, only: db1ink, db1val, get_status_message
        
        character(*), parameter :: rout_name = 'spline3_real'
        
        ! input / output
        integer, intent(in) :: ord                                              ! order
        real(dp), intent(in) :: x(:), xnew(:)                                   ! coordinates
        real(dp), intent(in) :: y(:)                                            ! function value
        real(dp), intent(out), optional :: ynew(:), dynew(:), d2ynew(:)         ! function values and derivatives
        logical, intent(in), optional :: extrap                                 ! whether extrapolation is allowed
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: n                                                            ! size of x, y, ...
        integer :: spline_init                                                  ! spline initialization parameter
        real(dp), allocatable :: spline_knots(:)                                ! knots of spline
        real(dp), allocatable :: spline_coeff(:)                                ! coefficients of spline
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize
        n = size(x)
        allocate(spline_coeff(n))
        allocate(spline_knots(n+ord))
        
        ! calculate coefficients
        call db1ink(x,n,y,ord,0,spline_knots,spline_coeff,ierr)
        err_msg = get_status_message(ierr)
        CHCKERR(err_msg)
        spline_init = 1
        do kd = 1,size(xnew)
            if (present(ynew)) then
                call db1val(xnew(kd),0,spline_knots,n,ord,&
                    &spline_coeff,ynew(kd),ierr,spline_init,extrap=extrap)
                err_msg = get_status_message(ierr)
                CHCKERR(err_msg)
            end if
            if (present(dynew)) then
                call db1val(xnew(kd),1,spline_knots,n,ord,&
                    &spline_coeff,dynew(kd),ierr,spline_init,extrap=extrap)
                err_msg = get_status_message(ierr)
                CHCKERR(err_msg)
            end if
            if (present(d2ynew)) then
                call db1val(xnew(kd),2,spline_knots,n,ord,&
                    &spline_coeff,d2ynew(kd),ierr,spline_init,extrap=extrap)
                err_msg = get_status_message(ierr)
                CHCKERR(err_msg)
            end if
        end do
    end function spline3_real
    integer function spline3_complex(ord,x,y,xnew,ynew,dynew,d2ynew,extrap) &
        &result(ierr)                                                           ! complex version
        use bspline_sub_module, only: db1ink, db1val, get_status_message
        
        character(*), parameter :: rout_name = 'spline3_complex'
        
        ! input / output
        integer, intent(in) :: ord                                              ! order
        real(dp), intent(in) :: x(:), xnew(:)                                   ! coordinates
        complex(dp), intent(in) :: y(:)                                         ! function value
        complex(dp), intent(out), optional :: ynew(:), dynew(:), d2ynew(:)      ! function values and derivatives
        logical, intent(in), optional :: extrap                                 ! whether extrapolation is allowed
        
        ! local variables
        integer :: kd, id                                                       ! counters
        integer :: n                                                            ! size of x, y, ...
        integer :: spline_init                                                  ! spline initialization parameter
        real(dp) :: dummy_var(2)                                                ! dummy variable
        real(dp), allocatable :: spline_knots(:,:)                              ! knots of spline
        real(dp), allocatable :: spline_coeff(:,:)                              ! coefficients of spline
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! initialize
        n = size(x)
        allocate(spline_coeff(n,2))
        allocate(spline_knots(n+ord,2))
        
        ! calculate coefficients for real part and complex part
        call db1ink(x,n,rp(y),ord,0,spline_knots(:,1),spline_coeff(:,1),ierr)
        err_msg = get_status_message(ierr)
        CHCKERR(err_msg)
        call db1ink(x,n,ip(y),ord,0,spline_knots(:,2),spline_coeff(:,2),ierr)
        err_msg = get_status_message(ierr)
        CHCKERR(err_msg)
        spline_init = 1
        do kd = 1,size(xnew)
            if (present(ynew)) then
                do id = 1,2
                    call db1val(xnew(kd),0,spline_knots(:,id),n,ord,&
                        &spline_coeff(:,id),dummy_var(id),ierr,spline_init,&
                        &extrap=extrap)
                    err_msg = get_status_message(ierr)
                    CHCKERR(err_msg)
                end do
                ynew(kd) = dummy_var(1) + iu*dummy_var(2)
            end if
            if (present(dynew)) then
                do id = 1,2
                    call db1val(xnew(kd),1,spline_knots(:,id),n,ord,&
                        &spline_coeff(:,id),dummy_var(id),ierr,spline_init,&
                        &extrap=extrap)
                    err_msg = get_status_message(ierr)
                    CHCKERR(err_msg)
                end do
                dynew(kd) = dummy_var(1) + iu*dummy_var(2)
            end if
            if (present(d2ynew)) then
                do id = 1,2
                    call db1val(xnew(kd),2,spline_knots(:,id),n,ord,&
                        &spline_coeff(:,id),dummy_var(id),ierr,spline_init,&
                        &extrap=extrap)
                    err_msg = get_status_message(ierr)
                    CHCKERR(err_msg)
                end do
                d2ynew(kd) = dummy_var(1) + iu*dummy_var(2)
            end if
        end do
    end function spline3_complex
end module splines
