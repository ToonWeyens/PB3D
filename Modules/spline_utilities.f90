!------------------------------------------------------------------------------!
!> Spline interpolation, wrapping the PSPLINE / EZspline library.
!!
!! Split out of num_utilities so that the (otherwise pure) numerical
!! utilities can be compiled and unit-tested without the PSPLINE external
!! dependency.
!------------------------------------------------------------------------------!
module spline_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use messages
    use num_vars, only: dp, iu, max_str_ln

    implicit none
    private
    public spline

    interface spline
        !> \public
        module procedure spline_real
        !> \public
        module procedure spline_complex
    end interface

contains
    
    !> \private real version
    integer function spline_real(x,y,xnew,ynew,ord,deriv,bcs,bcs_val,extrap) &
        &result(ierr)
        
        use EZspline_obj
        use EZspline
        
        character(*), parameter :: rout_name = 'spline_real'
        
        ! input / output
        real(dp), intent(in), target :: x(:)                                    !< coordinates
        real(dp), intent(in) :: y(:)                                            !< function value
        real(dp), intent(in), target :: xnew(:)                                 !< new coordinates
        real(dp), intent(out) :: ynew(:)                                        !< new function values
        integer, intent(in), optional :: ord                                    !< order [def 3]
        integer, intent(in), optional :: deriv                                  !< derivative [def 0]
        integer, intent(in), optional :: bcs(2)                                 !< boundary conditions [def 0]
        real(dp), intent(in), optional :: bcs_val(2)                            !< boundary conditions [no def]
        logical, intent(in), optional :: extrap                                 !< whether extrapolation is allowed [def .false.]
        
        ! local variables
        type(EZspline1_r8) :: f_spl                                             ! spline object
        integer :: ord_loc                                                      ! local order
        integer :: deriv_loc                                                    ! local deriv
        integer :: bcs_loc(2)                                                   ! local bcs
        integer :: kd                                                           ! counter
        integer :: n                                                            ! size of x, y
        integer :: nnew                                                         ! size of output x, y
        integer :: nnew_interp                                                  ! size of interpolated x, y, rest is extrapolated
        integer :: il(2)                                                        ! limits of interp., what falls outside is extrapolated
        real(dp), pointer :: x_loc(:)                                           ! x or -x
        real(dp), pointer :: xnew_loc(:)                                        ! xnew or -xnew
        real(dp) :: bcs_val_loc(2)                                              ! local bcs_val
        real(dp) :: lim_vals(0:3)                                               ! limit values at boundaries for extrapolation
        real(dp), allocatable :: xnew_EZ(:)                                     ! xnew in EZ spline doubles
        real(dp), allocatable :: ynew_EZ(:)                                     ! ynew in EZ spline doubles
        character(len=max_str_ln) :: err_msg                                    ! error message
        logical :: flip_x                                                       ! whether x axis is flipped
        logical :: extrap_loc                                                   ! local extrap
        
        ! initialize ierr
        ierr = 0
        
        ! set local variables
        ord_loc = 3
        if (present(ord)) ord_loc = ord
        deriv_loc = 0
        if (present(deriv)) deriv_loc = deriv
        bcs_loc = [0,0]
        if (present(bcs)) bcs_loc = bcs
        if (present(bcs_val)) bcs_val_loc = bcs_val
        extrap_loc = .false.
        if (present(extrap)) extrap_loc = extrap
        
        ! set other variables
        n = size(x)
        nnew = size(xnew)
        
        ! set local x and xnew, possibly flipped if negative
        flip_x = (x(2) .lt. x(1))
        if (flip_x) then
            allocate(x_loc(n))
            allocate(xnew_loc(nnew))
            x_loc = -x
            xnew_loc = -xnew
        else
            x_loc => x
            xnew_loc => xnew
        end if
        
#if ldebug
        ! check monotony
        do kd = 2,n
            if (x_loc(kd).le.x_loc(kd-1)) then
                ierr = 1
                err_msg = '|x| is not monotonously increasing'
                CHCKERR(err_msg)
            end if
        end do
        
        ! check array sizes
        if (size(y).ne.n) then
            ierr = 1
            err_msg = 'x and y need to have the same size'
            CHCKERR(err_msg)
        end if
        if (size(ynew).ne.nnew) then
            ierr = 1
            err_msg = 'xnew and ynew need to have the same size'
            CHCKERR(err_msg)
        end if
        
        ! check order
        if (ord_loc.lt.1 .or. ord_loc.gt.3) then
            ierr = 1
            err_msg = 'Only order 1 (linear), 2 (akima hermite) or 3 (cubic) &
                &possible'
            CHCKERR(err_msg)
        end if
        
        ! check derviative
        select case (deriv_loc)
            case (:-1)
                ierr = 1
                err_msg = 'only nonnegative degrees of derivative are possible'
                CHCKERR(err_msg)
            case (0:1)
                ! do nothing
            case (2:3)
                if (ord_loc.eq.1 .or. ord_loc.eq.2) then
                    ierr = 1
                    err_msg = 'Derivative of degree '//trim(i2str(deriv_loc))//&
                        &' not possible for order '//trim(i2str(ord_loc))
                    CHCKERR(err_msg)
                end if
            case (4:)
                ierr = 1
                err_msg = 'maximum degree of derivative is 2 for order 4'
                CHCKERR(err_msg)
        end select
        
        ! check boundary conditions
        select case (ord_loc)
            case (1)
                if (present(bcs) .or. present(bcs_val)) then
                    call writo('for order 1, no boundary conditions can be &
                        &prescribed',alert=.true.)
                end if
            case (2:3)
                ! check possibility
                do kd = 1,2
                    if (bcs_loc(kd).lt.-1 .or. bcs_loc(kd).gt.ord_loc-1) then
                        ierr = 1
                        err_msg = 'For order '//trim(i2str(ord_loc))//&
                            &' bcs has to be -1..'//trim(i2str(ord_loc-1))
                        CHCKERR(err_msg)
                    end if
                end do
                
                ! check whether derivatives are provided
                if ((any(bcs_loc.eq.1) .or. any(bcs_loc.eq.2)) &
                    &.and. .not.present(bcs_val)) then
                    ierr = 1
                    err_msg = 'When prescribing first or second derviatives, &
                        &need to provide bcs_val'
                    CHCKERR(err_msg)
                end if
        end select
#endif
        
        ! set up interpolation limits
        do kd = 1,nnew
            if (xnew_loc(kd).ge.x_loc(1)) exit
        end do
        il(1) = kd
        do kd = nnew,1,-1
            if (xnew_loc(kd).le.x_loc(n)) exit
        end do
        il(2) = kd
        nnew_interp = il(2)-il(1)+1
        
        ! check for extrapolation
        if ((il(1).ne.1 .or. il(2).ne.nnew) .and. &
            &.not.extrap_loc) then
            ierr = 1
            call writo('xnew = ['//trim(r2str(minval(xnew)))//'..'//&
                &trim(r2str(maxval(xnew)))//']')
            call writo('but x = ['//trim(r2str(minval(x)))//'..'//&
                &trim(r2str(maxval(x)))//']')
            err_msg = 'Extrapolation needed, but not allowed'
            CHCKERR(err_msg)
        end if
        
        ! initialize
        if (ord_loc.eq.1) then
            call EZlinear_init(f_spl,n,ierr)
            call EZspline_error(ierr)
            CHCKERR('')
        else
            call EZspline_init(f_spl,n,bcs_loc,ierr)
            call EZspline_error(ierr)
            CHCKERR('')
            if (ord_loc.eq.2) f_spl%isHermite = 1
        end if
        
        ! set grid
        f_spl%x1 = x_loc
        
        ! set boundary condition
        if (present(bcs_val)) then
            f_spl%bcval1min = bcs_val_loc(1)
            f_spl%bcval1max = bcs_val_loc(2)
            if (flip_x) then                                                    ! if x-axis flipped, also first derivative flipped
                if (bcs_loc(1) .eq. 1) f_spl%bcval1min = -f_spl%bcval1min
                if (bcs_loc(2) .eq. 1) f_spl%bcval1max = -f_spl%bcval1max
            end if
        end if
        
        ! set up 
        call EZspline_setup(f_spl,y,ierr,exact_dim=.true.)                      ! match exact dimensions
        call EZspline_error(ierr)
        CHCKERR('')
        ! interpolated part
        if (il(1).le.il(2)) then
            allocate(xnew_EZ(nnew_interp))
            allocate(ynew_EZ(nnew_interp))
            xnew_EZ = xnew_loc(il(1):il(2))
            
            if (deriv_loc.eq.0) then
                ! interpolate
                call EZspline_interp(f_spl,nnew_interp,xnew_EZ,ynew_EZ,ierr)
                call EZspline_error(ierr)
                CHCKERR('')
            else
                call EZspline_derivative(f_spl,deriv_loc,nnew_interp,xnew_EZ,&
                    &ynew_EZ,ierr) 
                call EZspline_error(ierr)
                CHCKERR('')
            end if
            
            ynew(il(1):il(2)) = ynew_EZ
            deallocate(xnew_EZ,ynew_EZ)
        end if
        
        ! extrapolated part
        if (extrap_loc) then
            ! extrapolate left
            if (il(1).gt.1) then
                ! set up limit values at first point
                ierr = setup_lim_vals(x_loc(1),lim_vals)
                CHCKERR('')
                
                ! calculate extrapolation
                call calc_extrap(xnew_loc(1:il(1)-1),x_loc(1),lim_vals,&
                    &ynew(1:il(1)-1))
            end if
            
            ! extrapolate right
            if (il(2).lt.nnew) then
                ! set up limit values at last point
                ierr = setup_lim_vals(x_loc(n),lim_vals)
                CHCKERR('')
                
                ! calculate extrapolation
                call calc_extrap(xnew_loc(il(2)+1:nnew),x_loc(n),lim_vals,&
                    &ynew(il(2)+1:nnew))
            end if
        end if
        
        ! free
        call EZspline_free(f_spl,ierr)
        CHCKERR('')
        call EZspline_error(ierr)
        CHCKERR('')
        if (flip_x) then
            deallocate(x_loc)
            deallocate(xnew_loc)
        end if
        nullify(x_loc)
        nullify(xnew_loc)
    contains
        !> /private set up limit values
        !!
        !! Makes use of f_spl
        integer function setup_lim_vals(xb,lim_vals) result(ierr)
            character(*), parameter :: rout_name = 'setup_lim_vals'
            
            ! input / output
            real(dp), intent(in) :: xb                                          ! x of boundary
            real(dp), intent(out) :: lim_vals(0:3)                              ! limit values
            
            ! local variables
            integer :: kd                                                       ! counter
            integer :: max_deriv                                                ! maximum degree of derivative
            
            ! initialize ierr
            ierr = 0
            
            ! initialize
            lim_vals = 0._dp
            select case(ord_loc)
                case (1:2)
                    max_deriv = 1
                case (3)
                    max_deriv = 3
            end select
            
            ! interpolation
            call EZspline_interp(f_spl,xb,lim_vals(0),ierr)
            call EZspline_error(ierr)
            CHCKERR('')
            
            ! derivatives
            do kd = 1,max_deriv
                call EZspline_derivative(f_spl,kd,xb,lim_vals(kd),ierr)
                call EZspline_error(ierr)
                CHCKERR('')
            end do
        end function setup_lim_vals
        
        !> /private calculate extrapolation
        subroutine calc_extrap(xnew,xb,lim_vals,ynew)
            ! input / output
            real(dp), intent(in) :: xnew(:)                                     ! xnew where to extrapolate
            real(dp), intent(in) :: xb                                          ! x of boundary
            real(dp), intent(in) :: lim_vals(0:2)                               ! limit values
            real(dp), intent(out) :: ynew(:)                                    ! y at xnew
            
            ! local variables
            real(dp), allocatable :: xdel(:)                                    ! delta x for extrapolation
            
            allocate(xdel(size(xnew)))
            xdel = xnew-xb
            
            select case (deriv_loc)
                case (0)
                    ynew = lim_vals(0) + &
                        &lim_vals(1)*xdel + &
                        &lim_vals(2)*xdel**2*0.5_dp
                case (1)
                    ynew = lim_vals(1) + &
                        &lim_vals(2)*xdel
                case (2)
                    ynew = lim_vals(2)
            end select
        end subroutine calc_extrap
    end function spline_real
    !> \private complex version
    integer function spline_complex(x,y,xnew,ynew,ord,deriv,bcs,bcs_val,&
        &extrap) result(ierr)
        
        use EZspline_obj
        use EZspline
        
        character(*), parameter :: rout_name = 'spline_complex'
        
        ! input / output
        real(dp), intent(in) :: x(:)                                            !< coordinates
        complex(dp), intent(in) :: y(:)                                         !< function value
        real(dp), intent(in) :: xnew(:)                                         !< new coordinates
        complex(dp), intent(out) :: ynew(:)                                     !< new function values
        integer, intent(in), optional :: ord                                    !< order [def 3]
        integer, intent(in), optional :: deriv                                  !< derivative [def 0]
        integer, intent(in), optional :: bcs(2)                                 !< boundary conditions [def 0]
        complex(dp), intent(in), optional :: bcs_val(2)                         !< boundary conditions [no def]
        logical, intent(in), optional :: extrap                                 !< whether extrapolation is allowed [def .false.]
        
        ! local variables
        real(dp), allocatable :: ynew_loc(:,:)                                  ! local ynew
        
        ! set up local variables
        allocate(ynew_loc(size(ynew),2))
        
        ! call real version for real part
        if (present(bcs_val)) then
            ierr = spline_real(x,rp(y),xnew,ynew_loc(:,1),ord,deriv,bcs,&
                &bcs_val=rp(bcs_val),extrap=extrap)
            CHCKERR('')
        else
            ierr = spline_real(x,rp(y),xnew,ynew_loc(:,1),ord,deriv,bcs,&
                &extrap=extrap)
            CHCKERR('')
        end if
        
        ! call real version for complex part
        if (present(bcs_val)) then
            ierr = spline_real(x,ip(y),xnew,ynew_loc(:,2),ord,deriv,bcs,&
                &bcs_val=ip(bcs_val),extrap=extrap)
            CHCKERR('')
        else
            ierr = spline_real(x,ip(y),xnew,ynew_loc(:,2),ord,deriv,bcs,&
                &extrap=extrap)
            CHCKERR('')
        end if
        
        ! save
        ynew = ynew_loc(:,1) + iu*ynew_loc(:,2)
    end function spline_complex
end module spline_utilities
