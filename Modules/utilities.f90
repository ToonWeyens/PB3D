!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module utilities
#include <PB3D_macros.h>
    use num_vars, only: dp, qp, iu, max_str_ln
    use messages, only: writo, print_ar_1, print_ar_2
    use str_ops, only: i2str, r2strt, r2str
    use output_ops, only: print_HDF5
    
    implicit none
    private
    public calc_zero_NR, calc_ext_var, calc_det, calc_int, add_arr_mult, c, &
        &calc_deriv, conv_FHM, check_deriv, calc_inv, interp_fun, calc_mult, &
        &init_utilities, derivs, con2dis, dis2con, round_with_tol, conv_sym, &
        &is_sym, calc_spline_3, con, &
        &plot_info
    
    ! the possible derivatives of order i
    integer, allocatable :: derivs_0(:,:)                                       ! all possible derivatives of order 0
    integer, allocatable :: derivs_1(:,:)                                       ! all possible derivatives of order 1
    integer, allocatable :: derivs_2(:,:)                                       ! all possible derivatives of order 2
    integer, allocatable :: derivs_3(:,:)                                       ! all possible derivatives of order 3
    integer, allocatable :: derivs_4(:,:)                                       ! all possible derivatives of order 3

    ! interfaces
    interface add_arr_mult
        module procedure add_arr_mult_3_3, add_arr_mult_3_1, add_arr_mult_1_1
    end interface
    interface calc_det
        module procedure calc_det_0D, calc_det_2D
    end interface
    interface calc_inv
        module procedure calc_inv_0D, calc_inv_2D
    end interface
    interface calc_mult
        module procedure calc_mult_0D_real, calc_mult_2D_real, &
            &calc_mult_2D_complex
    end interface
    interface calc_deriv
        module procedure calc_deriv_equidistant_real, &
            &calc_deriv_equidistant_complex, calc_deriv_regular_real, &
            &calc_deriv_regular_complex
    end interface
    interface calc_int
        module procedure calc_int_equidistant, calc_int_regular
    end interface
    interface round_with_tol
        module procedure round_with_tol_arr, round_with_tol_ind
    end interface
    interface con2dis
        module procedure con2dis_equidistant, con2dis_regular
    end interface
    interface dis2con
        module procedure dis2con_equidistant, dis2con_regular
    end interface
    interface interp_fun
        module procedure interp_fun_0D, interp_fun_1D, interp_fun_2D
    end interface
    interface con
        module procedure con_3D, con_2D, con_1D, con_0D
    end interface
    
    ! global variables for debugging
    logical :: plot_info = .false.                                              ! to plot information if wanted
    
contains
    ! initialize utilities:
    ! calculate all possible combinations of derivatives of a certain order
    subroutine init_utilities
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: ci, cj, ck, cl                                               ! counters
        
        ! allocate
        allocate(derivs_0(3,1))
        allocate(derivs_1(3,3))
        allocate(derivs_2(3,6))
        allocate(derivs_3(3,10))
        allocate(derivs_4(3,15))
        
        ci = 1
        cj = 1
        ck = 1
        cl = 1
        
        derivs_0 = 0
        derivs_1 = 0
        derivs_2 = 0
        derivs_3 = 0
        derivs_4 = 0
        
        do id = 1,3
            derivs_1(id,ci) = derivs_1(id,ci) + 1
            ci = ci+1
            do jd = 1,id
                derivs_2(id,cj) = derivs_2(id,cj) + 1
                derivs_2(jd,cj) = derivs_2(jd,cj) + 1
                cj = cj+1
                do kd = 1,jd
                    derivs_3(id,ck) = derivs_3(id,ck) + 1
                    derivs_3(jd,ck) = derivs_3(jd,ck) + 1
                    derivs_3(kd,ck) = derivs_3(kd,ck) + 1
                    ck = ck+1
                    do ld = 1,kd
                        derivs_4(id,cl) = derivs_4(id,cl) + 1
                        derivs_4(jd,cl) = derivs_4(jd,cl) + 1
                        derivs_4(kd,cl) = derivs_4(kd,cl) + 1
                        derivs_4(ld,cl) = derivs_4(ld,cl) + 1
                        cl = cl+1
                    end do
                end do
            end do
        end do
    end subroutine
    
    function derivs(order)
        ! input / output
        integer, intent(in) :: order
        integer, allocatable :: derivs(:,:)
        
        select case (order)
            case (0)
                allocate(derivs(3,size(derivs_0,2)))
                derivs = derivs_0
            case (1)
                allocate(derivs(3,size(derivs_1,2)))
                derivs = derivs_1
            case (2)
                allocate(derivs(3,size(derivs_2,2)))
                derivs = derivs_2
            case (3)
                allocate(derivs(3,size(derivs_3,2)))
                derivs = derivs_3
            case (4)
                allocate(derivs(3,size(derivs_4,2)))
                derivs = derivs_4
            case default
        end select
    end function
    
    ! numerically derives  a function  whose values are  given on  a equidistant
    ! grid, specified by the inverse step size  to an order specified by ord and
    ! a precision  specified by  prec (which is  the power of  the step  size to
    ! which the result  is still correct. E.g.: for forward  differences, prec =
    ! 0, and for central differences prec=1)
    integer function calc_deriv_equidistant_real(var,dvar,inv_step,ord,prec) &
        &result(ierr)                                                           ! equidistant version
        
        character(*), parameter :: rout_name = 'calc_deriv_equidistant_real'
        
        ! input / output
        real(dp), intent(in) :: var(:), inv_step
        real(dp), intent(inout) :: dvar(:)
        integer, intent(in) :: ord, prec
        
        ! local variables
        integer :: max_n
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: max_order(2) = [5,1]                                         ! maximum orders for precicions
        integer :: max_prec = 2                                                 ! maximum precision
        
        ! initialize ierr
        ierr = 0
        
        max_n = size(var)
        
        ! tests
        if (size(dvar).ne.max_n) then
            err_msg = 'Derived vector has to be of the same length as input &
                &vector'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (prec.lt.1 .or. prec.gt.max_prec) then
            err_msg = 'Precision '//trim(i2str(prec))//' not implemented'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (ord.lt.1 .or. ord.gt.max_order(prec)) then
            err_msg = 'For precision '//trim(i2str(prec))//&
                &', can only derive from order 1 up to order '//&
                &trim(i2str(max_order(prec)))//', not '//trim(i2str(ord))
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! choose correct precision
        select case (prec)
            case (1)
                call prec1
            case (2)
                call prec2
            case default
                err_msg = 'Precision of order '//trim(i2str(prec))//&
                    &' not implemented'
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        subroutine prec1
            ! test whether enough points are given
            if (max_n-2.lt.ord) then
                err_msg = 'For a derivative of order '//&
                    &trim(i2str(ord))//', with precision '//trim(i2str(prec))&
                    &//', at least '//trim(i2str(ord+2))//&
                    &' input values have to be passed'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! apply derivation rules precise up to order 1
            select case (ord)
                case(1)                                                         ! first derivative
                    ! first point
                    dvar(1) = (-3*var(1)+4*var(2)-var(3))*inv_step*0.5
                    ! middle points
                    dvar(2:max_n-1) = (-var(1:max_n-2)+var(3:max_n))*inv_step*0.5
                    ! last point
                    dvar(max_n) = (var(max_n-2) - 4*var(max_n-1) &
                        &+ 3*var(max_n))* inv_step*0.5
                case(2)                                                         ! second derivative
                    ! first point
                    dvar(1) = (2*var(1)-5*var(2)+4*var(3)-var(4))*inv_step**2
                    ! middle points
                    dvar(2:max_n-1) = (var(1:max_n-2)-2*var(2:max_n-1)&
                        &+var(3:max_n))*inv_step**2
                    ! last point
                    dvar(max_n) = (-var(max_n-3)+4*var(max_n-2)-5*var(max_n-1)&
                        &+2*var(max_n))*inv_step**2
                case(3)                                                         ! third derivative
                    ! first point
                    dvar(1) = (-5*var(1)+18*var(2)-24*var(3)+14*var(4)&
                        &-3*var(5))* inv_step**3*0.5
                    ! second point
                    dvar(2) = (-3*var(1)+10*var(2)-12*var(3)+6*var(4)-var(5))* &
                        &inv_step**3*0.5
                    ! middle points
                    dvar(3:max_n-2) = (-var(1:max_n-4)+2*var(2:max_n-3)-&
                        &2*var(4:max_n-1)+var(5:max_n))*inv_step**3*0.5
                    ! next-to-last point
                    dvar(max_n-1) = (var(max_n-4)-6*var(max_n-3)&
                        &+12*var(max_n-2)-10*var(max_n-1)+3*var(max_n))&
                        &*inv_step**3*0.5
                    ! last point
                    dvar(max_n) = (3*var(max_n-4)-14*var(max_n-3)&
                        &+24*var(max_n-2)-18*var(max_n-1)+5*var(max_n))&
                        &*inv_step**3/2
                case(4)                                                         ! fourth derivative
                    ! first point
                    dvar(1) = (3*var(1)-14*var(2)+26*var(3)-24*var(4)+11*var(5)&
                        &-2*var(6))*inv_step**4
                    ! second point
                    dvar(2) = (2*var(1)-9*var(2)+16*var(3)-14*var(4)+6*var(5)&
                        &-var(6))*inv_step**4
                    ! middle points
                    dvar(3:max_n-2) = (var(1:max_n-4)-4*var(2:max_n-3)&
                        &+6*var(3:max_n-2)-4*var(4:max_n-1)+var(5:max_n))&
                        &*inv_step**4
                    ! next-to-last point
                    dvar(max_n-1) = (-var(max_n-5)+6*var(max_n-4)&
                        &-14*var(max_n-3)+16*var(max_n-2)-9*var(max_n-1)&
                        &+2*var(max_n))*inv_step**4
                    ! last point
                    dvar(max_n) = (-2*var(max_n-5)+11*var(max_n-4)-&
                        &24*var(max_n-3)+26*var(max_n-2)-14*var(max_n-1)&
                        &+3*var(max_n))*inv_step**4
                case(5)                                                         ! fifth derivative
                    ! first point
                    dvar(1) = (-7*var(1)+40*var(2)-95*var(3)+120*var(4)&
                        &-85*var(5)+32*var(6)-5*var(7))*inv_step**5*0.5
                    ! second point
                    dvar(2) = (-5*var(1)+28*var(2)-65*var(3)+80*var(4)&
                        &-55*var(5)+20*var(6)-3*var(7))*inv_step**5*0.5
                    ! third point
                    dvar(3) = (-3*var(1)+16*var(2)-35*var(3)+40*var(4)&
                        &-25*var(5)+8*var(6)-var(7))*inv_step**5*0.5
                    ! middle points
                    dvar(4:max_n-3) = (-var(1:max_n-6)+4*var(2:max_n-5)&
                        &-5*var(3:max_n-4)+5*var(5:max_n-2)-4*var(6:max_n-1)&
                        &+var(7:max_n))*inv_step**5*0.5
                    ! next-to-next-to-last point
                    dvar(max_n-2) = (var(max_n-6)-8*var(max_n-5)&
                        &+25*var(max_n-4)-40*var(max_n-3)+35*var(max_n-2)&
                        &-16*var(max_n-1)+3*var(max_n))*inv_step**5*0.5
                    ! next-to-last point
                    dvar(max_n-1) = (3*var(max_n-6)-20*var(max_n-5)&
                        &+55*var(max_n-4)-80*var(max_n-3)+65*var(max_n-2)&
                        &-28*var(max_n-1)+5*var(max_n))*inv_step**5*0.5
                    ! last point
                    dvar(max_n) = (5*var(max_n-6)-32*var(max_n-5)&
                        &+85*var(max_n-4)-120*var(max_n-3)+95*var(max_n-2)&
                        &-40*var(max_n-1)+7*var(max_n))*inv_step**5*0.5
                case default
                    ! This you should never reach!
                    err_msg = 'Derivation of order '//trim(i2str(ord))//&
                        &' not possible in calc_deriv_equidistant'
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end subroutine
        
        subroutine prec2
            ! test whether enough points are given
            if (max_n-3.lt.ord) then
                err_msg = 'For a derivative of order '//&
                    &trim(i2str(ord))//', with precision '//trim(i2str(prec))&
                    &//', at least '//trim(i2str(ord+3))//&
                    &' input values have to be passed'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! apply derivation rules precise up to order 1
            select case (ord)
                case(1)                                                         ! first derivative
                    ! first point
                    dvar(1:1) = inv_step*&
                        &(-11*var(1)+18*var(2)-9*var(3)+2*var(4))/6
                    ! second point
                    dvar(2:2) = inv_step*&
                        &(-2*var(1)-3*var(2)+6*var(3)-var(4))/6
                    ! middle points
                    ! (take 5 points instead of 4 because of symmetry)
                    dvar(3:max_n-2) = inv_step*&
                        &(var(1:max_n-4)-8*var(2:max_n-3)+&
                        &8*var(4:max_n-1)-var(5:max_n))/12
                    ! next-to-last point
                    dvar(max_n-1:max_n-1) = -inv_step*&                         ! - because odd number (1) of sign changes
                        &(-2*var(max_n)-3*var(max_n-1)+6*var(max_n-2)-&
                        &var(max_n-3))/6
                    ! last point
                    dvar(max_n:max_n) = -inv_step*&                             ! - because odd number (1) of sign changes
                        &(-11*var(max_n)+18*var(max_n-1)-9*var(max_n-2)+&
                        &2*var(max_n-3))/6
                case default
                    ! This you should never reach!
                    err_msg = 'Derivation of order '//trim(i2str(ord))//&
                        &' not possible in calc_deriv_equidistant'
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end subroutine
    end function calc_deriv_equidistant_real
    integer function calc_deriv_equidistant_complex(var,dvar,inv_step,ord,prec) &
        &result(ierr)                                                           ! equidistant complex version
        
        character(*), parameter :: rout_name = 'calc_deriv_equidistant_complex'
        
        ! input / output
        complex(dp), intent(in) :: var(:)
        real(dp), intent(in) :: inv_step
        complex(dp), intent(inout) :: dvar(:)
        integer, intent(in) :: ord, prec
        
        ! local variables
        real(dp), allocatable :: dvar_loc(:)
        
        ! set up local copy of component of dvar
        allocate(dvar_loc(size(var)))
        
        ! call real version for real part
        ierr = calc_deriv_equidistant_real(realpart(var),dvar_loc,&
            &inv_step,ord,prec)
        CHCKERR('')
        
        ! update dvar with local copy
        dvar = dvar_loc
        
        ! call real version for imaginary part
        ierr = calc_deriv_equidistant_real(imagpart(var),dvar_loc,&
            &inv_step,ord,prec)
        CHCKERR('')
        
        ! update dvar with local copy
        dvar = dvar + iu*dvar_loc
        
        ! deallocate local variables
        deallocate(dvar_loc)
    end function calc_deriv_equidistant_complex
    
    ! numerically derives  a function whose values  are given on a  regular, but
    ! not  equidistant, grid,  to  an order  specified by  ord  and a  precision
    ! specified by prec (which is the power of the step size to which the result
    ! is still correct. E.g.: for forward differences, prec = 0, and for central
    ! differences prec=1)
    integer function calc_deriv_regular_real(var,dvar,x,ord,prec) result(ierr)  ! regular, non-equidistant version
        character(*), parameter :: rout_name = 'calc_deriv_regular'
        
        ! input / output
        real(dp), intent(in) :: var(:)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: dvar(:)
        integer, intent(in) :: ord, prec
        
        ! local variables
        integer :: max_n
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: max_order(2) = [3,1]                                         ! maximum orders for precisions
        integer :: max_prec = 2                                                 ! maximum precision
        real(dp), allocatable, target :: delta(:)                               ! delta(i) = x(i+1)-x(i) (called delta_(i+1/2) in text)
        real(dp), pointer :: a(:), b(:), c(:), d(:)                             ! pointers to parts of delta
        
        ! initialize ierr
        ierr = 0
        
        max_n = size(var)
        
        ! tests
        if (size(dvar).ne.max_n) then
            err_msg = 'Derived vector has to be of the same length as input &
                &vector'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(x).ne.max_n) then
            err_msg = 'Abscissa vector has to be of the same length as &
                &ordinate vector'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (prec.lt.1 .or. prec.gt.max_prec) then
            err_msg = 'Precision '//trim(i2str(prec))//' not implemented'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (ord.lt.1 .or. ord.gt.max_order(prec)) then
            err_msg = 'For precision '//trim(i2str(prec))//&
                &', can only derive from order 1 up to order '//&
                &trim(i2str(max_order(prec)))//', not '//trim(i2str(ord))
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set up delta
        allocate(delta(max_n-1))
        delta(1:max_n-1) = x(2:max_n) - x(1:max_n-1)
        
        ! choose correct precision
        select case (prec)
            case (1)
                call prec1
            case (2)
                call prec2
            case default
                err_msg = 'Precision of order '//trim(i2str(prec))//&
                    &' not implemented'
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
        subroutine prec1
            ! test whether enough points are given
            if (max_n-2.lt.ord) then
                err_msg = 'For a derivative of order '//&
                    &trim(i2str(ord))//', with precision '//trim(i2str(prec))&
                    &//', at least '//trim(i2str(ord+2))//&
                    &' input values have to be passed'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! apply derivation rules precise up to order 1
            select case (ord)
                case(1)                                                         ! first derivative
                    ! first point
                    a => delta(1:1)
                    b => delta(2:2)
                    dvar(1:1) = (-(2*a+b)*b*var(1)+(a+b)**2*var(2)-&
                        &a**2*var(3)) / (a*(a+b)*b)
                    nullify(a,b)
                    ! middle points
                    a => delta(1:max_n-2)
                    b => delta(2:max_n-1)
                    dvar(2:max_n-1) = (-b**2*var(1:max_n-2)+&
                        &(b**2-a**2)*var(2:max_n-1)+ a**2*var(3:max_n)) / &
                        &(a*(a+b)*b)
                    nullify(a,b)
                    ! last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    dvar(max_n:max_n) = -(-(2*a+b)*b*var(max_n)+(a+b)**2&       ! - because odd number (1) of sign changes
                        &*var(max_n-1)-a**2*var(max_n-2)) / (a*(a+b)*b)
                    nullify(a,b)
                case(2)                                                         ! second derivative
                    ! first point
                    a => delta(1:1)
                    b => delta(2:2)
                    c => delta(3:3)
                    dvar(1:1) = 2*(b*(b+c)*c*(3*a+2*b+c)*var(1)-&
                        &(a+b)*(a+b+c)*c*(2*a+2*b+c)*var(2)+&
                        &(b+c)*(a+b+c)*a*(2*a+b+c)*var(3)-&
                        &a*(a+b)*b*(2*a+b)*var(4)) / &
                        &(a*b*c*(a+b)*(b+c)*(a+b+c))
                    nullify(a,b,c)
                    ! middle points
                    a => delta(1:max_n-2)
                    b => delta(2:max_n-1)
                    dvar(2:max_n-1) = 2*(b*var(1:max_n-2)-&
                        &(a+b)*var(2:max_n-1)+ a*var(3:max_n)) / &
                        &(a*(a+b)*b)
                    nullify(a,b)
                    ! last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    c => delta(max_n-3:max_n-3)
                    dvar(max_n:max_n) = 2*(b*(b+c)*c*(3*a+2*b+c)*var(max_n)-&
                        &(a+b)*(a+b+c)*c*(2*a+2*b+c)*var(max_n-1)+&
                        &(b+c)*(a+b+c)*a*(2*a+b+c)*var(max_n-2)-&
                        &a*(a+b)*b*(2*a+b)*var(max_n-3)) / &
                        &(a*b*c*(a+b)*(b+c)*(a+b+c))
                    nullify(a,b,c)
                case(3)                                                         ! third derivative
                    ! first point
                    a => delta(1:1)
                    b => delta(2:2)
                    c => delta(3:3)
                    d => delta(4:4)
                    dvar(1:1) = 6*(&
                        &-b*c*d*(b+c)*(c+d)*(b+c+d)*(4*a+3*b+2*c+d)*var(1:1)&
                        &+c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*(3*a+3*b+2*c+d)*&
                        &var(2:2)&
                        &-a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*(3*a+2*b+2*c+d)*&
                        &var(3:3)&
                        &+a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*(3*a+2*b+c+d)*&
                        &var(4:4)&
                        &-a*b*c*(b+c)*(a+b)*(a+b+c)*(3*a+2*b+c)*&
                        &var(5:5)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                    ! second point
                    a => delta(1:1)
                    b => delta(2:2)
                    c => delta(3:3)
                    d => delta(4:4)
                    dvar(2:2) = 6*(&
                        &-b*c*d*(b+c)*(c+d)*(b+c+d)*(3*b+2*c+d)*var(1:1)&
                        &-c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*(a-3*b-2*c-d)*&
                        &var(2:2)&
                        &+a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*(a-2*b-2*c-d)*&
                        &var(3:3)&
                        &-a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*(a-2*b-c-d)*&
                        &var(4:4)&
                        &+a*b*c*(b+c)*(a+b)*(a+b+c)*(a-2*b-c)*&
                        &var(5:5)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                    ! middle points
                    a => delta(1:max_n-4)
                    b => delta(2:max_n-3)
                    c => delta(3:max_n-2)
                    d => delta(4:max_n-1)
                    dvar(3:max_n-2) = 6*(&
                        &b*c*d*(b+c)*(c+d)*(b+c+d)*(b-2*c-d)*&
                        &var(1:max_n-4)&
                        &-c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*(a+b-2*c-d)*&
                        &var(2:max_n-3)&
                        &+a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*(a+2*b-2*c-d)*&
                        &var(3:max_n-2)&
                        &-a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*(a+2*b-c-d)*&
                        &var(4:max_n-1)&
                        &+a*b*c*(b+c)*(a+b)*(a+b+c)*(a+2*b-c)*&
                        &var(5:max_n)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                    ! next-to-last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    c => delta(max_n-3:max_n-3)
                    d => delta(max_n-4:max_n-4)
                    dvar(max_n-1:max_n-1) = -6*(&                               ! - because odd number (3) of sign changes
                        &-b*c*d*(b+c)*(c+d)*(b+c+d)*(3*b+2*c+d)*&
                        &var(max_n:max_n)&
                        &-c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*(a-3*b-2*c-d)*&
                        &var(max_n-1:max_n-1)&
                        &+a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*(a-2*b-2*c-d)*&
                        &var(max_n-2:max_n-2)&
                        &-a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*(a-2*b-c-d)*&
                        &var(max_n-3:max_n-3)&
                        &+a*b*c*(b+c)*(a+b)*(a+b+c)*(a-2*b-c)*&
                        &var(max_n-4:max_n-4)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                    ! next-to-last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    c => delta(max_n-3:max_n-3)
                    d => delta(max_n-4:max_n-4)
                    dvar(max_n-1:max_n-1) = -6*(&                               ! - because odd number (3) of sign changes
                        &-b*c*d*(b+c)*(c+d)*(b+c+d)*(3*b+2*c+d)*&
                        &var(max_n:max_n)&
                        &-c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*(a-3*b-2*c-d)*&
                        &var(max_n-1:max_n-1)&
                        &+a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*(a-2*b-2*c-d)*&
                        &var(max_n-2:max_n-2)&
                        &-a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*(a-2*b-c-d)*&
                        &var(max_n-3:max_n-3)&
                        &+a*b*c*(b+c)*(a+b)*(a+b+c)*(a-2*b-c)*&
                        &var(max_n-4:max_n-4)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                    ! last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    c => delta(max_n-3:max_n-3)
                    d => delta(max_n-4:max_n-4)
                    dvar(max_n:max_n) = -6*(&                                   ! -1 because odd number (3) of sign changes
                        &-b*c*d*(b+c)*(c+d)*(b+c+d)*(4*a+3*b+2*c+d)*&
                        &var(max_n:max_n)&
                        &+c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*(3*a+3*b+2*c+d)*&
                        &var(max_n-1:max_n-1)&
                        &-a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*(3*a+2*b+2*c+d)*&
                        &var(max_n-2:max_n-2)&
                        &+a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*(3*a+2*b+c+d)*&
                        &var(max_n-3:max_n-3)&
                        &-a*b*c*(b+c)*(a+b)*(a+b+c)*(3*a+2*b+c)*&
                        &var(max_n-4:max_n-4)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                case default
                    ! This you should never reach!
                    err_msg = 'Derivation of order '//trim(i2str(ord))//&
                        &' not possible in calc_deriv_regular'
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end subroutine
        
        subroutine prec2
            ! test whether enough points are given
            if (max_n-3.lt.ord) then
                err_msg = 'For a derivative of order '//&
                    &trim(i2str(ord))//', with precision '//trim(i2str(prec))&
                    &//', at least '//trim(i2str(ord+3))//&
                    &' input values have to be passed'
                ierr = 1
                CHCKERR(err_msg)
            end if
            
            ! apply derivation rules precise up to order 1
            select case (ord)
                case(1)                                                         ! first derivative
                    ! first point
                    a => delta(1:1)
                    b => delta(2:2)
                    c => delta(3:3)
                    dvar(1:1) = (&
                        &-b*(b+c)*c*(a*(3*a+2*b+c)+b*(2*a+b)+c*(a+b))*var(1)&
                        &+(a+b)*(a+b+c)*c*(a+b)*(a+b+c)*var(2)&
                        &-(b+c)*(a+b+c)*a*a*(a+b+c)*var(3)&
                        &+a*(a+b)*b*a*(a+b)*var(4)) / &
                        &(a*b*c*(a+b)*(b+c)*(a+b+c))
                    nullify(a,b,c)
                    ! second point
                    a => delta(1:1)
                    b => delta(2:2)
                    c => delta(3:3)
                    dvar(2:2) = (&
                        &-b*(b+c)*c*b*(b+c)*var(1)&
                        &+(a+b)*(a+b+c)*c*(b*c-a*c+b**2-2*a*b)*var(2)&
                        &+(b+c)*(a+b+c)*a*a*(b+c)*var(3)&
                        &-a*(a+b)*b*a*b*var(4)) / &
                        &(a*b*c*(a+b)*(b+c)*(a+b+c))
                    nullify(a,b,c)
                    ! middle points
                    ! (take 5 points instead of 4 because of symmetry)
                    a => delta(1:max_n-4)
                    b => delta(2:max_n-3)
                    c => delta(3:max_n-2)
                    d => delta(4:max_n-1)
                    dvar(3:max_n-2) = (&
                        &b*c*d*(b+c)*(c+d)*(b+c+d)*b*c*(c+d)*&
                        &var(1:max_n-4)&
                        &-c*d*(a+b)*(c+d)*(a+b+c)*(a+b+c+d)*c*(a+b)*(c+d)*&
                        &var(2:max_n-3)&
                        &+a*d*(b+c)*(a+b+c)*(b+c+d)*(a+b+c+d)*&
                        &(a*d*(c-b)+a*c**2-d*b**2-2*b*c*(a+b-c-d))*&
                        &var(3:max_n-2)&
                        &+a*b*(a+b)*(c+d)*(b+c+d)*(a+b+c+d)*b*(a+b)*(c+d)*&
                        &var(4:max_n-1)&
                        &-a*b*c*(b+c)*(a+b)*(a+b+c)*b*c*(a+b)*&
                        &var(5:max_n)) / &
                        &(a*b*c*d*(a+b)*(b+c)*(c+d)*(a+b+c)*(b+c+d)*(a+b+c+d))
                    nullify(a,b,c,d)
                    ! next-to-last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    c => delta(max_n-3:max_n-3)
                    dvar(max_n-1:max_n-1) = -(&                                 ! - because odd number (1) of sign changes
                        &-b*(b+c)*c*b*(b+c)*var(max_n)&
                        &+(a+b)*(a+b+c)*c*(b*c-a*c+b**2-2*a*b)*var(max_n-1)&
                        &+(b+c)*(a+b+c)*a*a*(b+c)*var(max_n-2)&
                        &-a*(a+b)*b*a*b*var(max_n-3)) / &
                        &(a*b*c*(a+b)*(b+c)*(a+b+c))
                    nullify(a,b,c)
                    ! last point
                    a => delta(max_n-1:max_n-1)
                    b => delta(max_n-2:max_n-2)
                    c => delta(max_n-3:max_n-3)
                    dvar(max_n:max_n) = -(&                                     ! - because odd number (1) of sign changes
                        &-b*(b+c)*c*(a*(3*a+2*b+c)+b*(2*a+b)+c*(a+b))*&
                        &var(max_n)&
                        &+(a+b)*(a+b+c)*c*(a+b)*(a+b+c)*var(max_n-1)&
                        &-(b+c)*(a+b+c)*a*a*(a+b+c)*var(max_n-2)&
                        &+a*(a+b)*b*a*(a+b)*var(max_n-3)) / &
                        &(a*b*c*(a+b)*(b+c)*(a+b+c))
                    nullify(a,b,c)
                case default
                    ! This you should never reach!
                    err_msg = 'Derivation of order '//trim(i2str(ord))//&
                        &' not possible in calc_deriv_regular'
                    ierr = 1
                    CHCKERR(err_msg)
            end select
        end subroutine
    end function calc_deriv_regular_real
    integer function calc_deriv_regular_complex(var,dvar,x,ord,prec) &
        &result(ierr)                                                           ! regular, non-equidistant, complex version
        
        character(*), parameter :: rout_name = 'calc_deriv_regular_complex'
        
        ! input / output
        complex(dp), intent(in) :: var(:)
        real(dp), intent(in) :: x(:)
        complex(dp), intent(inout) :: dvar(:)
        integer, intent(in) :: ord, prec
        
        ! local variables
        real(dp), allocatable :: dvar_loc(:)
        
        ! set up local copy of component of dvar
        allocate(dvar_loc(size(var)))
        
        ! call real version for real part
        ierr = calc_deriv_regular_real(realpart(var),dvar_loc,x,ord,prec)
        CHCKERR('')
        
        ! update dvar with local copy
        dvar = dvar_loc
        
        ! call real version for imaginary part
        ierr = calc_deriv_regular_real(imagpart(var),dvar_loc,x,ord,prec)
        CHCKERR('')
        
        ! update dvar with local copy
        dvar = dvar + iu*dvar_loc
        
        ! deallocate local variables
        deallocate(dvar_loc)
    end function calc_deriv_regular_complex
    
    ! Calculates the coefficients of a cubic  spline through a number of points,
    ! which  can  later  be  used  to  calculate  the  interpolating  values  or
    ! derivatives thereof.
    ! The information about the spline is saved in the spline_info object.
    ! Based on http://www.math.ntnu.no/emner/TMA4215/2008h/cubicsplines.pdf
    !!!! FOR NOW, JUST USE CUBIC SPLINE IN 1 D AND DON'T THINK ABOUT THE OTHER DIMENSIONS !!!!
    integer function calc_spline_3(x,y,spline_info) result(ierr)
        use num_vars, only: spline_type
        
        character(*), parameter :: rout_name = 'calc_spline_3'
        
        ! input / output
        real(dp), intent(in) :: x(:), y(:)                                      ! abscissa and ordenate
        type(spline_type), intent(inout) :: spline_info                         ! stores the spline information
        
        ! local variables
        integer :: n_pt                                                         ! nr. of points
        real(dp), allocatable :: h(:)                                           ! h(i) = x(i+1)-x(i)
        real(dp), allocatable :: b(:)                                           ! b(i) = (y(i+1)-y(i))/h(i)
        real(dp), allocatable :: v(:)                                           ! v(i) = 2(h(i-1)+h(i))
        real(dp), allocatable :: u(:)                                           ! u(i) = 6(b(i)-b(i-1))
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set n_pt
        n_pt = size(x)
        
        ! tests
        if (size(y).ne.n_pt) then
            err_msg = 'x and y need to have the same dimensions'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        if (allocated(spline_info%x) .or. allocated(spline_info%y) .or. &
            &allocated(spline_info%z)) then
            call writo('WARING: spline_info was already in use and will be &
                &overwritten')
            if (allocated(spline_info%x)) deallocate(spline_info%x)
            if (allocated(spline_info%y)) deallocate(spline_info%y)
            if (allocated(spline_info%z)) deallocate(spline_info%z)
        end if
        
        ! set up objects in spline_info
        allocate(spline_info%x(n_pt))
        allocate(spline_info%y(n_pt))
        allocate(spline_info%z(n_pt))
        spline_info%x = x
        spline_info%y = y
        
        ! precalculations in order to get z
        allocate(h(n_pt),b(n_pt),u(n_pt),v(n_pt))
        h(1:n_pt-1) = x(2:n_pt)-x(1:n_pt-1)                                     ! h(n_pt) has no meaning
        b(1:n_pt-1) = (y(2:n_pt)-y(1:n_pt-1))/h(1:n_pt-1)                       ! b(n_pt) has no meaning
        v(2:n_pt) = 2*(h(2:n_pt)+h(1:n_pt-1))                                   ! v(1) has no meaning
        u(2:n_pt) = 6*(b(2:n_pt)+b(1:n_pt-1))                                   ! u(1) has no meaning
        
        ierr = 1
        err_msg = 'SPLINES NOT YET IMPLEMENTED!!!'
        CHCKERR(err_msg)
    end function calc_spline_3

    ! numerically interpolates a function that is given on either FM to HM or to
    ! FM. If FM2HM is .true., the starting variable is FM, if .false., it is HM
    integer function conv_FHM(var,cvar,FM2HM) result(ierr)
        character(*), parameter :: rout_name = 'conv_FHM'
        
        ! input / output
        real(dp), intent(in) :: var(:)
        real(dp), intent(inout) :: cvar(:)
        logical, intent(in) :: FM2HM
        
        ! local variables
        integer :: max_n
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        max_n = size(var)
        
        ! tests
        if (size(cvar).ne.max_n) then
            err_msg = 'The converted vector has to be of the same length as &
                &the input vector'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        if (FM2HM) then                                                         ! FM to HM
            cvar(1) = 0.0_dp
            cvar(2:max_n) = (var(1:max_n-1)+var(2:max_n))/2.0_dp
        else                                                                    ! HM to FM
            cvar(1) = (3.0_dp*var(2)-var(3))/2.0_dp
            cvar(2:max_n-1) = (var(2:max_n-1)+var(3:max_n))/2.0_dp
            cvar(max_n) = (-var(max_n-1)+3.0_dp*var(max_n))/2.0_dp
        end if
    end function conv_FHM
    
    ! checks  whether the  derivatives requested  for a  certain subroutine  are
    ! valid
    integer function check_deriv(deriv,max_deriv,sr_name) result(ierr)
        character(*), parameter :: rout_name = 'check_deriv'
        
        ! input / output
        integer, intent(in) :: deriv(3)
        integer, intent(in) :: max_deriv
        character(len=*), intent(in) :: sr_name
        
        ! local variables
        integer :: id
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test the derivatives
        do id = 1, 3
            if (deriv(id).gt.max_deriv) then
                err_msg = 'Asking '//trim(sr_name)//' for a derivative&
                    & in the '//trim(i2str(id))//'th dimension of order '&
                    &//trim(i2str(deriv(id)))//', while the maximum order is '&
                    &//trim(i2str(max_deriv))
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do
    end function

    ! Integrates a function using the trapezoidal rule:
    !   int_1^n f(x) dx = sum_k=1^(n-1) {(f(k+1)+f(k))*(x(k+1)-x(k))/2},
    ! with n the number of points. So, n  points have to be specified as well as
    ! n values  for the function  to be interpolated. They  have to be  given in
    ! ascending order but the step size does not have to be constant
    ! this yields the following difference formula:
    !   int_1^n f(x) dx = int_1^(n-1) f(x) dx + (f(n)+f(n-1))*(x(n)-x(n-1))/2,
    ! which is used here
    integer function calc_int_regular(var,x,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_regular'
        
        ! input / output
        real(dp), intent(inout) :: var_int(:)                                   ! integrated variable
        real(dp), intent(in) :: var(:)                                          ! variable to be integrated
        real(dp), intent(in) :: x(:)                                            ! abscissa
        
        ! local variables
        integer :: n_max
        integer :: kd
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set n_max
        n_max = size(var)
        
        ! tests
        if (size(x).ne.n_max) then
            err_msg = 'The arrays of points and values are not of the same &
                &length'
            ierr = 1
            CHCKERR(err_msg)
        else if (n_max.lt.2) then
            err_msg = 'Asking to integrate with '//trim(i2str(n_max))&
                &//' points. Need at least 2'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! calculate integral for all points
        var_int = 0.0_dp
        
        do kd = 2, n_max
            var_int(kd) = var_int(kd-1) + &
                &(var(kd)+var(kd-1))/2 * (x(kd)-x(kd-1))
        end do
    end function calc_int_regular
    
    ! Integrates a function using the trapezoidal rule:
    !   int_1^n f(x) dx = sum_k=1^(n-1) {(f(k+1)+f(k))*delta_x/2},
    ! with n  the number of points,  which are assumed to be  equidistant with a
    ! given step size delta_x.
    integer function calc_int_equidistant(var,step_size,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_equidistant'
        
        ! input / output
        real(dp), intent(inout) :: var_int(:)                                   ! integrated variable
        real(dp), intent(in) :: var(:)                                          ! variable to be integrated
        real(dp), intent(in) :: step_size                                       ! step size of abscissa
        
        ! local variables
        integer :: n_max
        integer :: kd
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set n_max
        n_max = size(var)
        
        ! tests
        if (n_max.lt.2) then
            err_msg = 'Asking to integrate with '//trim(i2str(n_max))&
                &//' points. Need at least 2'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! calculate integral for all points
        var_int = 0.0_dp
        
        do kd = 2, n_max
            var_int(kd) = var_int(kd-1) + (var(kd)+var(kd-1))/2 * step_size
        end do
    end function calc_int_equidistant
    
    ! extrapolates  a   function,  using  linear  or   quadratic  interpolation,
    ! depending on  the number of  points and values  given. The data  should be
    ! given sorted, in ascending order, without duplicate points in var_points.
    ! It uses the following solution: A b = c
    ! | x_1^0 x_1^1 x_1^2 ... | | a_0 |   | var(1) |
    ! | x_2^0 x_2^1 x_2^2 ... | | a_1 | = | var(2) |
    ! | x_3^0 x_3^1 x_3^2 ... | | a_2 |   | var(3) |
    ! | ...   ...   ...   ... | | ... |   |  ...   |
    ! for the interpolating polynomial ext_var = a_0 + a_1 x + a_2 x^2 + ...
    ! There is an optional flag to look  for the k^th derivative of the function
    ! instead of the function itself, where k should be lower than the degree of
    ! the polynomial
    integer function calc_ext_var(ext_var,var,var_points,ext_point,deriv_in) &
        &result(ierr)
        character(*), parameter :: rout_name = 'ext_var'
        
        ! input / output
        real(dp), intent(inout) :: ext_var                                      ! output
        real(dp), intent(in) :: var(:), var_points(:)                           ! input function and points at which tabulated
        real(dp), intent(in) :: ext_point                                       ! point to which extrapolate
        integer, intent(in), optional :: deriv_in                               ! specifies an optional derivative
        
        ! local variables
        integer :: id, jd
        integer :: istat
        integer :: pol_deg
        integer :: deriv
        real(dp), allocatable :: A(:,:), LU(:,:)
        real(dp), allocatable :: a_i(:)
        real(dp) :: fact                                                        ! multiplicative factor, see below
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! determine the degree of the polynomial ext_var
        pol_deg = size(var)
        
        ! tests
        if (size(var).ne.size(var_points)) then
            err_msg = 'The size of the arrays containing the variable and &
                &the points do not match'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (present(deriv_in)) then 
            deriv = deriv_in
            if (deriv.ge.pol_deg) then                                          ! order of derivative too hight
                err_msg = 'Asking ext_var for '//trim(i2str(deriv_in))&
                    &//'th derivative, while using a polynomial of degree '&
                    &//trim(i2str(pol_deg))
                ierr = 1
                CHCKERR(err_msg)
            else if (deriv.lt.0) then
                err_msg = 'Asking ext_var for a derivative of order '&
                    &//trim(i2str(deriv_in))
                ierr = 1
                CHCKERR(err_msg)
            end if
        else
            deriv = 0
        end if
        
        ! fill matrix
        allocate(A(pol_deg,pol_deg))
        allocate(LU(pol_deg,pol_deg))
        A(:,1) = 1.0_dp
        col: do jd = 2, pol_deg
            row: do id = 1, pol_deg
                A(id,jd) = var_points(id)**(jd-1)
            end do row
        end do col
        
        ! solve the system for a_i
        allocate(a_i(pol_deg))
        a_i = var
        LU = A
        call dgesv( pol_deg, 1, LU, pol_deg, [(id,id=1,pol_deg)], &
            &a_i, pol_deg, istat)
        
        
        if (istat.ne.0) then
            err_msg = 'The solution could not be found. Are the extrapolating &
                &points independent?'
            ierr = 1
            if (ierr.ne.0) then
                call writo('The matrix equation to be solved was Ax=b, A = ')
                call print_ar_2(A)
                call writo('and b = ')
                call print_ar_1(var)
            end if
            CHCKERR(err_msg)
        end if
        
        ext_var = 0.0_dp
        ! for deriv of order i: 
        ! d^i/dx^i f = sum_j = i ^ pol_deg ( j*(j-1)*...*(j-i) a_j x^(j-i) )
        polynom: do id = 1+deriv, pol_deg
            ! construct the factor 
            fact = 1.0_dp
            d_dx: do jd = 1,deriv
                fact = fact * float((id-1)-(jd-1))
            end do d_dx
            ! update ext_var with current polynomial term a_i(id)
            ext_var = ext_var + fact*a_i(id)*ext_point**(id-1-deriv)
        end do polynom
    end function calc_ext_var
    
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    integer function add_arr_mult_3_3(arr_1,arr_2,arr_3,deriv) result(ierr)     ! Version with arr_1 and arr_2 in three coords.
        character(*), parameter :: rout_name = 'add_arr_mult_3_3'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: arr_2(1:,1:,1:,0:,0:,0:)
        real(dp), intent(out) :: arr_3(1:,1:,1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac(3)                                                  ! binomial factor for norm., pol. and tor. sum
        integer :: r, m, n                                                      ! current degree of norm., pol., tor. derivatives
        integer :: kd                                                           ! normal counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(arr_1,1).ne.size(arr_2,1) .or. size(arr_1,2).ne.size(arr_2,2) &
            &.or. size(arr_1,3).ne.size(arr_2,3)) then
            err_msg = 'Arrays 1 and 2 need to have the same size'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(arr_1,1).ne.size(arr_3,1) .or. size(arr_1,2).ne.size(arr_3,2) &
            &.or. size(arr_1,3).ne.size(arr_3,3)) then
            err_msg = 'Arrays 1 and 2 need to have the same size as the &
                &resulting array 3'
            ierr = 1
            CHCKERR(err_msg)
        end if
        do kd = 1,3
            if (size(arr_1,3+kd).lt.deriv(kd)+1 .or. &
                &size(arr_2,3+kd).lt.deriv(kd)+1) then
                err_msg = 'Arrays 1 and 2 do not provide the necessary &
                    &orders of derivatives to calculate a derivative of order &
                    &['//trim(i2str(deriv(1)))//','//trim(i2str(deriv(2)))//&
                    &','//trim(i2str(deriv(3)))//' in coordinate '&
                    &//trim(i2str(kd))
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do
        
        ! distribute normal and angular derivatives using binomial theorem
        ! normal derivatives
        do r = 0,deriv(1)
            if (r.eq.0) then                                                    ! first term in sum
                bin_fac(1) = 1.0_dp
            else
                bin_fac(1) = bin_fac(1)*(deriv(1)-(r-1))/r
            end if
            ! poloidal derivatives
            do m = 0,deriv(2)
                if (m.eq.0) then                                                ! first term in sum
                    bin_fac(2) = 1.0_dp
                else
                    bin_fac(2) = bin_fac(2)*(deriv(2)-(m-1))/m
                end if
                ! toroidal derivatives
                do n = 0,deriv(3)
                    if (n.eq.0) then                                            ! first term in sum
                        bin_fac(3) = 1.0_dp
                    else
                        bin_fac(3) = bin_fac(3)*(deriv(3)-(n-1))/n
                    end if
                    
                    ! current term in the tripple summation
                    arr_3 = arr_3 + &
                        &bin_fac(1)*bin_fac(2)*bin_fac(3) &
                        &* arr_1(:,:,:,r,m,n)&
                        &* arr_2(:,:,:,deriv(1)-r,deriv(2)-m,deriv(3)-n)
                end do
            end do
        end do
    end function add_arr_mult_3_3
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    integer function add_arr_mult_3_1(arr_1,arr_2,arr_3,deriv) result(ierr)     ! Version with arr_1 in three coords and arr_2 only in the flux coord.
        character(*), parameter :: rout_name = 'add_arr_mult_3_1'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,1:,0:,0:,0:)
        real(dp), intent(in) :: arr_2(1:,0:)
        real(dp), intent(out) :: arr_3(1:,1:,1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac                                                     ! binomial factor for norm. sum
        integer :: r                                                            ! current degree of norm. derivative
        integer :: kd                                                           ! normal counter
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(arr_1,3).ne.size(arr_2,1)) then
            err_msg = 'Arrays 1 and 2 need to have the same size in the flux &
                &coord.'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(arr_1,1).ne.size(arr_3,1) .or. size(arr_1,2).ne.size(arr_3,2) &
            &.or. size(arr_1,3).ne.size(arr_3,3)) then
            err_msg = 'Array 1 needs to have the same size as the resulting &
                &array 3'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! distribute   normal  derivatives   using  binomial   theorem,  angular
        ! derivatives only operate on the first array
        ! normal derivatives
        do r = 0,deriv(1)
            if (r.eq.0) then                                                    ! first term in sum
                bin_fac = 1.0_dp
            else
                bin_fac = bin_fac*(deriv(1)-(r-1))/r
            end if
            
            ! current term in the tripple summation
            do kd = 1, size(arr_1,3)
                arr_3(:,:,kd) = arr_3(:,:,kd) + bin_fac * &
                    &arr_1(:,:,kd,r,deriv(2),deriv(3))* arr_2(kd,deriv(1)-r)
            end do
        end do
    end function add_arr_mult_3_1
    ! Add to an  array (3) the product  of arrays (1) and  (2), with derivatives
    ! that are distributed between both acording to the binomial theorem
    integer function add_arr_mult_1_1(arr_1,arr_2,arr_3,deriv) result(ierr)     ! Version with arr_1 and arr_2 only in the flux coord.
        character(*), parameter :: rout_name = 'add_arr_mult_1_1'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,0:)
        real(dp), intent(in) :: arr_2(1:,0:)
        real(dp), intent(out) :: arr_3(1:)
        integer, intent(in) :: deriv(3)
        
        ! local variables
        real(dp) :: bin_fac                                                     ! binomial factor for norm. sum
        integer :: r                                                            ! current degree of norm. derivative
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(arr_1,1).ne.size(arr_2,1)) then
            err_msg = 'Arrays 1 and 2 need to have the same size in the &
                &flux coord.'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(arr_1,1).ne.size(arr_3,1)) then
            err_msg = 'Array 1 needs to have the same size as the resulting &
                &array 3'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        if (deriv(2).gt.0 .or. deriv(3).gt.0) then                              ! no addition to arr_3
            return
        else
            ! distribute  normal  derivatives  using binomial  theorem,  angular
            ! derivatives only operate on the first array
            ! normal derivatives
            do r = 0,deriv(1)
                if (r.eq.0) then                                                    ! first term in sum
                    bin_fac = 1.0_dp
                else
                    bin_fac = bin_fac*(deriv(1)-(r-1))/r
                end if
                
                ! current term in the tripple summation
                arr_3 = arr_3 + bin_fac * arr_1(:,r)* arr_2(:,deriv(1)-r)
            end do
        end if
    end function add_arr_mult_1_1

    ! Calculate determinant of a matrix which is defined on a 3D grid . The size
    ! of  the  matrix  (last  two  indices)  should  be  small,  as  the  direct
    ! formula employing cofactors  is used. The storage  convention described in
    ! metric_type is used.
    integer recursive function calc_det_2D(detA,A,n) result (ierr)
        character(*), parameter :: rout_name = 'calc_det_2D'
        
        ! input / output
        real(dp), intent(inout) :: detA(:,:,:)                                  ! output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! input
        integer, intent(in) :: n                                                ! size of matrix
        
        ! local variables
        logical :: sym                                                          ! .true. if matrix symmetric
        integer :: id, jd, kd                                                   ! counter
        integer :: nn                                                           ! nr. of elements in matrix
        real(dp) :: sgn                                                         ! holds either plus or minus one
        integer, allocatable :: slct(:)                                         ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: work(:,:,:)                                    ! work array
        integer, allocatable :: c_sub(:)                                        ! indices of submatrix in storage convention of metric_type
        integer, allocatable :: k_sub(:)                                        ! first indices of submatrix in 2D
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(detA,1).ne.size(A,1) .or. size(detA,2).ne.size(A,2) .or. &
            &size(detA,3).ne.size(A,3)) then
            err_msg = 'The output and input matrix have to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set total number of elements in matrix
        nn = size(A,4)
        
        ! set sym
        ierr = is_sym(n,nn,sym)
        CHCKERR('')
        
        ! intialize other quantities
        sgn = 1.0_dp
        allocate (idx(n))
        idx = [(id,id=1,n)]
        allocate (slct(n))
        slct = 1
        allocate (work(size(A,1),size(A,2),size(A,3)))
        
        detA = 0.0_dp
        
        if (n.eq.1) then                                                        ! shouldn't be used
            detA = A(:,:,:,c([1,1],sym,n))
        else if (n.eq.2) then
            if (plot_info) write(*,*) '   n = 2 sym = ', sym
            detA = A(:,:,:,c([1,1],sym,n))*A(:,:,:,c([2,2],sym,n))-&
                &A(:,:,:,c([1,2],sym,n))*A(:,:,:,c([2,1],sym,n))
            if (plot_info) then
                write(*,*) 'A(1,1):'
                call print_HDF5('X','X',A(:,:,:,c([1,1],sym,n)),[size(A,1),size(A,2),size(A,3)])
                read(*,*)
                write(*,*) 'A(2,2):'
                call print_HDF5('X','X',A(:,:,:,c([2,2],sym,n)),[size(A,1),size(A,2),size(A,3)])
                read(*,*)
                write(*,*) 'A(1,2):'
                call print_HDF5('X','X',A(:,:,:,c([1,2],sym,n)),[size(A,1),size(A,2),size(A,3)])
                read(*,*)
                write(*,*) 'A(2,1):'
                call print_HDF5('X','X',A(:,:,:,c([2,1],sym,n)),[size(A,1),size(A,2),size(A,3)])
                read(*,*)
            end if
        else
            if (plot_info) write(*,*) '   n = '//trim(i2str(n))//' sym = ', sym
            ! allocate coordinates of (n-1)x(n-1) submatrix
            allocate(c_sub((n-1)**2))
            allocate(k_sub(n**2)); k_sub = 0
            ! iterate over indices in second dimension
            do id = 1, n
                if (plot_info) then
                    write(*,*) 'id = ', id
                    write(*,*) 'previous det:'
                    call print_HDF5('X','X',detA,[size(A,1),size(A,2),size(A,3)])
                    read(*,*)
                end if
                slct(id) = 0
                ! set up coordinates of submatrix (2:n,pack(idx,slct.gt.0))
                ! (in general not symmetric, even if parent matrix is)
                k_sub = pack(idx,slct.gt.0)
                do kd = 1,n-1
                    do jd = 1,n-1
                        c_sub((kd-1)*(n-1)+jd) = c([jd+1,k_sub(kd)],sym,n)
                    end do
                end do
                ierr = calc_det_2D(work,A(:,:,:,c_sub),n-1)
                CHCKERR('')
                detA = detA + A(:,:,:,c([1,id],sym,n))*sgn*work
                sgn = -sgn
                slct(id) = 1
                if (plot_info) then
                    write(*,*) 'plotting A of ', c([1,id],sym,n)
                        call print_HDF5('X','X',A(:,:,:,c([1,id],sym,n)),&
                            &[size(A,1),size(A,2),size(A,3)])
                    read(*,*)
                    write(*,*) 'plotting subdet '
                        call print_HDF5('X','X',work,&
                            &[size(A,1),size(A,2),size(A,3)])
                    read(*,*)
                    write(*,*) 'sign = ', -sgn
                    write(*,*) 'updated det:'
                    call print_HDF5('X','X',detA,[size(A,1),size(A,2),size(A,3)])
                    read(*,*)
                end if
            end do
            ! deallocate c_sub and k_sub
            deallocate(c_sub,k_sub)
        end if
    end function calc_det_2D
    ! calculate determinant of a constant matrix
    ! (adapted from http://dualm.wordpress.com/2012/01/06/computing-determinant-
    ! in-fortran/)
    integer function calc_det_0D(det_0D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_det_0D'
        
        ! input / output
        real(dp), intent(inout) :: det_0D                                       ! output
        real(dp), intent(in) :: A(:,:)                                          ! input
        
        ! local variables
        integer :: n 
        real(dp) :: sgn
        integer :: id, info
        integer, allocatable :: ipiv(:)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: A_loc(:,:)                                     ! copy of A
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,1).ne.size(A,2)) then
            err_msg = 'The matrix A has to be square'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initialize
        n = size(A,1)
        allocate(ipiv(n))
        ipiv = 0
        
        ! set up local A
        allocate(A_loc(n,n))
        A_loc = A
        
        call dgetrf(n, n, A_loc, n, ipiv, info)
        
        det_0D = 1.0_dp
        do id = 1, n
            det_0D = det_0D*A_loc(id, id)
        end do
        
        sgn = 1.0_dp
        do id = 1, n
            if(ipiv(id).ne.id) then
                sgn = -sgn
            end if
        end do
        det_0D = sgn*det_0D
        
        ! deallocate local variables
        deallocate(A_loc)
    end function calc_det_0D
    
    ! calculate inverse of square matrix A  which has elements depending on a 3D
    ! grid. The storage convention described in metric_type is used.
    ! This method uses direct inversion using  Cramer's rule, since the matrix A
    ! is supposed to  be very small (i.e.  3x3) and since the inverse  has to be
    ! calculated at each of the points in the grid, this can be quite fast.
    integer function calc_inv_2D(inv_2D,A,n) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_2D'
        
        ! input / output
        real(dp), intent(inout) :: inv_2D(:,:,:,:)                              ! output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! input
        integer, intent(in) :: n                                                ! size of matrix
        
        ! local variables
        logical :: sym                                                          ! .true. if matrix symmetric
        real(dp), allocatable :: detA(:,:,:)                                    ! determinant of a submatrix of A
        integer :: id, jd, kd, ld                                               ! counters
        integer :: nn                                                           ! nr. of elements in matrix
        integer, allocatable :: slct(:,:)                                       ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: c_sub(:)                                        ! indices of submatrix in storage convention of metric_type
        integer, allocatable :: l_sub(:), k_sub(:)                              ! first and second indices of submatrix in 2D
        integer :: kd_min                                                       ! min. of kd
       
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(inv_2D,1).ne.size(A,1) .or. size(inv_2D,2).ne.size(A,2) .or. &
            &size(inv_2D,3).ne.size(A,3) .or. size(inv_2D,4).ne.size(A,4)) then
            err_msg = 'The output and input matrix have to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set total number of elements in matrix
        nn = size(A,4)
        
        ! set sym
        ierr = is_sym(n,nn,sym)
        CHCKERR('')
        
        ! initialize other quantitites
        allocate (idx(n))
        idx = [(id,id=1,n)]
        allocate (slct(n,2))
        slct = 1
        allocate(detA(size(A,1),size(A,2),size(A,3)))
        detA = 0.0_dp
        
        ! allocate coordinates of (n-1)x(n-1) submatrix
        allocate(c_sub((n-1)**2))
        allocate(l_sub(n**2)); l_sub = 0
        allocate(k_sub(n**2)); k_sub = 0
        
        ! set up kd_min
        kd_min = 1
        
        ! calculate cofactor matrix, replacing original elements
        ! use is made of detA
        do ld = 1,n
            ! adjust kd_min if symmetric
            if (sym) kd_min = ld
            do kd = kd_min,n
                if (plot_info) write(*,*) 'in calc_inv kd,ld = ', kd, ld
                ! set up coordinates of submatrix (strike out row i, column j)
                slct(ld,1) = 0                                                  ! take out column ld
                slct(kd,2) = 0                                                  ! take out row kd
                l_sub = pack(idx,slct(:,1).gt.0)
                k_sub = pack(idx,slct(:,2).gt.0)
                do jd = 1,n-1
                    do id = 1,n-1
                        c_sub((jd-1)*(n-1)+id) = c([l_sub(id),k_sub(jd)],sym,n) ! taking the transpose!
                    end do
                end do
                if (plot_info) write(*,*) 'c_sub = ', c_sub
                ierr = calc_det_2D(detA,A(:,:,:,c_sub),n-1)
                CHCKERR('')
                inv_2D(:,:,:,c([kd,ld],sym,n)) = (-1.0_dp)**(kd+ld)*detA
                slct(ld,1) = 1
                slct(kd,2) = 1
                if (plot_info) then
                    write(*,*) 'inv( = ',kd,ld,')'
                    call print_HDF5('X','X',inv_2D(:,:,:,c([kd,ld],sym,n)),&
                        &[size(inv_2D,1),size(inv_2D,2),size(inv_2D,3)])
                    write(*,*) 'last row: ', inv_2D(:,:,size(inv_2D,3),c([kd,ld],sym,n))
                    read(*,*)
                end if
            end do
        end do
        
        ! deallocate c_sub l_sub and k_sub
        deallocate(c_sub,l_sub,k_sub)
        
        ! calculate determinant in detA
        ierr = calc_det(detA,A,n)
        CHCKERR('')
        
        if (plot_info) then
            write(*,*) 'detA = '
            call print_HDF5('X','X',detA,&
                &[size(inv_2D,1),size(inv_2D,2),size(inv_2D,3)])
            read(*,*)
        end if
        ! divide by determinant
        do kd = 1,nn
            inv_2D(:,:,:,kd) = inv_2D(:,:,:,kd) / detA
            
            if (plot_info) then
                write(*,*) 'FINAL RESULT', kd
                call print_HDF5('X','X',inv_2D(:,:,:,kd),&
                    &[size(inv_2D,1),size(inv_2D,2),size(inv_2D,3)])
                read(*,*)
            end if
        end do
    end function calc_inv_2D
    ! calculate inverse of a constant square matrix
    integer function calc_inv_0D(inv_0D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_0D'
        
        ! input / output
        real(dp), intent(inout) :: inv_0D(:,:)                                  ! output
        real(dp), intent(in) :: A(:,:)                                          ! input
        
        ! local variables
        integer :: n 
        integer, allocatable :: ipiv(:)                                         ! pivot variable, used by Lapack
        real(dp), allocatable :: w(:)                                           ! work variable
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,1).ne.size(A,1)) then
            err_msg = 'The input matrix A has to be square'
            ierr = 1
            CHCKERR(err_msg)
        end if
        if (size(inv_0D,1).ne.size(A,1) .or. size(inv_0D,2).ne.size(A,2)) then
            err_msg = 'The output and input matrix have to have the same sizes'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! initialize
        n = size(A,1)
        allocate(ipiv(n))
        ipiv = 0
        allocate(w(n))
        
        inv_0D = A
        
        call dgetrf(n, n, inv_0D, n, ipiv, ierr)                                ! LU factorization
        err_msg = 'Lapack couldn''t find the LU factorization'
        CHCKERR(err_msg)
        
        call dgetri(n, inv_0D, n, ipiv, w, n, ierr)                             ! inverse of LU
        CHCKERR('Lapack couldn''t compute the inverse')
    end function calc_inv_0D
    
    ! Calculate matrix multiplication of two square matrices AB = A B which have
    ! elements  defined  on a  3D  grid.  The  storage convention  described  in
    ! metric_type is used.
    ! Optionally, A and/or B can be transposed.
    integer function calc_mult_2D_real(A,B,AB,n,transp) result(ierr)            ! real version with matrix defined on 3D grid
        character(*), parameter :: rout_name = 'calc_mult_2D_real'
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! input A
        real(dp), intent(in) :: B(:,:,:,:)                                      ! input B
        real(dp), intent(inout) :: AB(:,:,:,:)                                  ! output A B
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp(2)                              ! .true. if A and/or B transposed
        
        ! local variables
        logical :: sym(3)                                                       ! .true. if matrices symmetric
        integer :: nn(3)                                                        ! nr. of elements in matrices
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counters
        integer :: id_min                                                       ! min. of id
        integer :: c1, c2, c3                                                   ! coordinates
        integer :: ind(2,3)                                                     ! indices
        logical :: transp_loc(2)                                                ! local copy of transp
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,1).ne.size(B,1) .or. size(A,2).ne.size(B,2) .or. &
            &size(A,3).ne.size(B,3) .or. size(A,1).ne.size(AB,1) .or. &
            &size(A,2).ne.size(AB,2) .or. size(A,3).ne.size(AB,3)) then
            err_msg = 'The output and input matrices have to defined on the &
                &same grid'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set local transp
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        
        ! set total number of elements in matrices
        nn(1) = size(A,4)
        nn(2) = size(B,4)
        nn(3) = size(AB,4)
        
        ! set sym
        do id = 1,3
            ierr = is_sym(n,nn(id),sym(id))
            CHCKERR('')
        end do
        
        ! initialize AB
        AB = 0._dp
        
        ! set up id_min
        id_min = 1
        
        ! loop over all 2D rows and columns of AB
        do jd = 1,n
            ! adjust id_min if AB symmetric
            if (sym(3)) id_min = jd
            do id = id_min,n
                do kd = 1,n
                    ind(:,3) = [id,jd]
                    c3 = c(ind(:,3),sym(3),n)
                    if (transp_loc(1)) then
                        ind(:,1) = [kd,id]
                    else
                        ind(:,1) = [id,kd]
                    end if
                    c1 = c(ind(:,1),sym(1),n)
                    if (transp_loc(2)) then
                        ind(:,2) = [jd,kd]
                    else
                        ind(:,2) = [kd,jd]
                    end if
                    c2 = c(ind(:,2),sym(2),n)
                    AB(:,:,:,c3) = AB(:,:,:,c3) + &
                        &A(:,:,:,c1)*B(:,:,:,c2)
                end do
            end do
        end do
    end function calc_mult_2D_real
    integer function calc_mult_2D_complex(A,B,AB,n,transp) result(ierr)         ! complex version with matrix defined on 3D grid
        character(*), parameter :: rout_name = 'calc_mult_2D_complex'
        
        ! input / output
        complex(dp), intent(in) :: A(:,:,:,:)                                   ! input A
        complex(dp), intent(in) :: B(:,:,:,:)                                   ! input B
        complex(dp), intent(inout) :: AB(:,:,:,:)                               ! output A B
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp(2)                              ! .true. if A and/or B transposed
        
        ! local variables
        logical :: sym(3)                                                       ! .true. if matrices symmetric
        integer :: nn(3)                                                        ! nr. of elements in matrices
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counters
        integer :: id_min                                                       ! min. of id
        integer :: c1, c2, c3                                                   ! coordinates
        logical :: transp_loc(2)                                                ! local copy of transp
        integer :: ind(2,3)                                                     ! indices
        integer :: d(3)                                                         ! dimension of matrices A and B
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A,1).ne.size(B,1) .or. size(A,2).ne.size(B,2) .or. &
            &size(A,3).ne.size(B,3) .or. size(A,1).ne.size(AB,1) .or. &
            &size(A,2).ne.size(AB,2) .or. size(A,3).ne.size(AB,3)) then
            err_msg = 'The output and input matrices have to defined on the &
                &same grid'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! set local transp
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        
        ! set d
        d = [size(A,1),size(A,2),size(A,3)]
        
        ! set total number of elements in matrices
        nn(1) = size(A,4)
        nn(2) = size(B,4)
        nn(3) = size(AB,4)
        
        ! set sym
        do id = 1,3
            ierr = is_sym(n,nn(id),sym(id))
            CHCKERR('')
        end do
        
        ! initialize AB
        AB = 0._dp
        
        ! set up id_min
        id_min = 1
        
        ! loop over all 2D rows and columns of AB
        do jd = 1,n
            ! adjust id_min if AB symmetric
            if (sym(3)) id_min = jd
            do id = id_min,n
                do kd = 1,n
                    ind(:,3) = [id,jd]
                    c3 = c(ind(:,3),sym(3),n)
                    if (transp_loc(1)) then
                        ind(:,1) = [kd,id]
                    else
                        ind(:,1) = [id,kd]
                    end if
                    c1 = c(ind(:,1),sym(1),n)
                    if (transp_loc(2)) then
                        ind(:,2) = [jd,kd]
                    else
                        ind(:,2) = [kd,jd]
                    end if
                    c2 = c(ind(:,2),sym(2),n)
                    AB(:,:,:,c3) = AB(:,:,:,c3) + &
                        &con(A(:,:,:,c1),ind(:,1),sym(1),d)*&
                        &con(B(:,:,:,c2),ind(:,2),sym(2),d)
                end do
            end do
        end do
    end function calc_mult_2D_complex
    integer function calc_mult_0D_real(A,B,AB,n) result(ierr)                   ! real version with constant matrix
        character(*), parameter :: rout_name = 'calc_mult_0D_real'
        
        ! input / output
        real(dp), intent(in) :: A(:)                                            ! input A
        real(dp), intent(in) :: B(:)                                            ! input B
        real(dp), intent(inout) :: AB(:)                                        ! output A B
        integer, intent(in) :: n                                                ! size of matrix
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: loc_AB(:,:,:,:)                                ! local copy of C
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(A).ne.size(B) .or. size(A).ne.size(AB)) then
            err_msg = 'The matrices A, B and AB need to have the same shape'
            ierr = 1
            CHCKERR(err_msg)
        end if
        
        ! allocate local AB
        allocate(loc_AB(1,1,1,size(AB)))
        
        ! calculate using 2D version
        ierr = calc_mult_2D_real(reshape(A,[1,1,1,size(A)]),&
            &reshape(B,[1,1,1,size(B)]),loc_AB,n)
        
        ! update AB with local value
        AB = loc_AB(1,1,1,:)
        
        ! deallocate local variables
        deallocate(loc_AB)
    end function calc_mult_0D_real
    
    ! Converts a  (symmetric) matrix  A defined  on a 3D  grid with  the storage
    ! convention  described in  metric_type. If  the matrix  is stored  with n^2
    ! numbers, only  the lower diagonal elements  are kept in matrix  B. If only
    ! the  lower diagonal  elements are  stored, they  are copied  to the  upper
    ! diagonal ones of the matrix B as well.
    ! Note  that  the  routine  does  not check  whether  the  matrix is  indeed
    ! symmetric and that information may thus be lost after conversion.
    integer function conv_sym(A,B,n) result(ierr)
        character(*), parameter :: rout_name = 'conv_sym'
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! matrix to be converted
        real(dp), intent(inout) :: B(:,:,:,:)                                   ! converted matrix
        integer, intent(in) :: n                                                ! size of matrix
        
        ! local variables
        logical :: sym(2)                                                       ! .true. if matrices symmetric
        integer :: nn(2)                                                        ! nr. of elements in matrices
        integer :: id, jd                                                       ! counters
        integer :: id_min                                                       ! minimum value of id
        
        ! initialize ierr
        ierr = 0
        
        ! set total number of elements in matrix
        nn(1) = size(A,4)
        nn(2) = size(B,4)
        
        ! set sym
        do id = 1,2
            ierr = is_sym(n,nn(id),sym(id))
            CHCKERR('')
        end do
        
        ! copy A into B
        if ((sym(1).neqv.sym(2))) then                                          ! matrices A and B have different storage type
            ! set id_min
            id_min = 1
            ! convert according to sym of A
            do jd = 1,n
                if (sym(2)) id_min = jd
                do id = id_min,n
                    B(:,:,:,c([id,jd],sym(2),n)) = A(:,:,:,c([id,jd],sym(1),n))
                end do
            end do
        else
            B = A
        end if
    end function conv_sym
    
    ! finds the zero of a function using Newton-Rhapson iteration
    integer function calc_zero_NR(zero_NR,fun,dfun,guess) result(ierr)
        use num_vars, only: max_it_NR, tol_NR
        
        character(*), parameter :: rout_name = 'calc_zero_NR'
        
        ! input / output
        real(dp), intent(inout) :: zero_NR                                      ! output
        real(dp), intent(in) :: guess                                           ! first guess
        interface
            function fun(x)                                                     ! the function
                use num_vars, only: dp
                real(dp) :: fun
                real(dp), intent(in) :: x
            end function fun
            function dfun(x)                                                    ! derivative of the function
                use num_vars, only: dp
                real(dp) :: dfun
                real(dp), intent(in) :: x
            end function dfun
        end interface
        
        ! local variables
        integer :: jd
        real(dp) :: corr
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), parameter :: relax_fac = 0.5_dp                               ! factor for relaxation
        
        ! initialize ierr
        ierr = 0
        
        zero_NR = guess
        
        NR: do jd = 1,max_it_NR
            ! correction to theta_NR
            corr = -fun(zero_NR)/dfun(zero_NR)
            !write(*,*) 'jd, zero_NR = ', jd, zero_NR
            !write(*,*) 'fun, dfun = ', fun(zero_NR), dfun(zero_NR)
            zero_NR = zero_NR + relax_fac*corr
            
            ! check for convergence
            if (abs(corr).lt.tol_NR) then
                return
            else if (jd .eq. max_it_NR) then
                err_msg = 'Not converged after '//trim(i2str(jd))//&
                    &' iterations, with residual '//trim(r2str(corr))//&
                    &' and final value '//trim(r2str(zero_NR))
                zero_NR = 0.0_dp
                ierr = 1
                CHCKERR(err_msg)
            end if
        end do NR
    end function calc_zero_NR

    ! Convert  between  points  from  a  continuous grid  to  a  discrete  grid,
    ! providing either  the the limits on  the grid (lim_c and  lim_d), in which
    ! case the grid is assumed to be equidistant, or the grid values themselves,
    ! in which case the grid just has to be regular.
    ! The output is a  real value where the floored integer is  the index in the
    ! discrete grid  and the remainder  corresponds to the fraction  towards the
    ! next index.  If no solution  is found, a  negative value is  outputted, as
    ! well as a message.
    subroutine con2dis_equidistant(pt_c,pt_d,lim_c,lim_d)                       ! equidistant version
        ! input / output
        real(dp), intent(in) :: pt_c                                            ! point on continous grid
        real(dp), intent(inout) :: pt_d                                         ! point on discrete grid
        real(dp), intent(in) :: lim_c(2)                                        ! [min_c,max_c]
        integer, intent(in) :: lim_d(2)                                         ! [min_d,max_d]
        
        ! Calculate, using the formula:
        !    pt_c - min_c     pt_d - min_d
        !   ------------- =  ------------- ,
        !   max_c - min_c    max_d - min_d
        pt_d = lim_d(1) + (lim_d(2)-lim_d(1)) * (pt_c-lim_c(1)) / &
            &(lim_c(2)-lim_c(1))
    end subroutine con2dis_equidistant
    subroutine con2dis_regular(pt_c,pt_d,var_c)                                 ! regular grid version
        ! input / output
        real(dp), intent(in) :: pt_c                                            ! point on continous grid
        real(dp), intent(inout) :: pt_d                                         ! point on discrete grid
        real(dp), intent(in) :: var_c(:)                                        ! continous grid values
        
        ! local variables
        real(dp), allocatable :: var_c_inv(:)                                   ! inverted var_c
        integer :: size_c                                                       ! size of var_c
        integer :: ind_lo, ind_hi                                               ! lower and upper index comprising pt_c
        integer :: id                                                           ! counter
        
        ! set up var_c_inv
        size_c = size(var_c)
        allocate(var_c_inv(size_c))
        var_c_inv = var_c(size_c:1:-1)
        
        ! find the lower and upper index comprising pt_c
        ind_lo = 0
        ind_hi = size_c+1
        do id = 1,size_c
            if (var_c(id).le.pt_c) ind_lo = id
            if (var_c_inv(id).ge.pt_c) ind_hi = size_c+1-id
        end do
        
        ! tests
        if (ind_lo.eq.0 .or. ind_hi.eq.size_c+1) then                           ! not within range
            call writo('WARNING: pt_c not within range',persistent=.true.)
            pt_d = -1._dp
            return
        end if
        
        ! set output
        if (ind_lo.lt.ind_hi) then                                              ! valid output that does not correspond to a point on grid
            pt_d = ind_lo + (pt_c-var_c(ind_lo))/(var_c(ind_hi)-var_c(ind_lo))
        else if (ind_lo.eq.ind_hi) then                                         ! valid output that does correspond to a point on grid
            pt_d = ind_lo
        else                                                                    ! invalid output
            call writo('WARNING: ind_lo cannot be higher than ind_hi',&
                &persistent=.true.)
            pt_d = -1._dp
            return
        end if
        
        ! deallocate
        deallocate(var_c_inv)
    end subroutine con2dis_regular
    
    ! Convert  between  points  from  a  discrete grid  to  a  continuous  grid,
    ! providing either  the the limits on  the grid (lim_c and  lim_d), in which
    ! case the grid is assumed to be equidistant, or the grid values themselves,
    ! in which case the grid just has to be regular.
    ! The output is a real value. If  the discrete value lies outside the range,
    ! in the case of a regular grid, a negative value is outputted, as well as a
    ! message.
    subroutine dis2con_equidistant(pt_d,pt_c,lim_d,lim_c)                       ! equidistant version
        ! input / output
        integer, intent(in) :: pt_d                                             ! point on discrete grid
        real(dp), intent(inout) :: pt_c                                         ! point on continous grid
        integer, intent(in) :: lim_d(2)                                         ! [min_d,max_d]
        real(dp), intent(in) :: lim_c(2)                                        ! [min_c,max_c]
        
        ! Calculate, using the formula:
        !    pt_c - min_c     pt_d - min_d
        !   ------------- =  ------------- ,
        !   max_c - min_c    max_d - min_d
        pt_c = lim_c(1) + (lim_c(2)-lim_c(1)) * (pt_d-lim_d(1)) / &
            &(lim_d(2)-lim_d(1))
    end subroutine dis2con_equidistant
    subroutine dis2con_regular(pt_d,pt_c,var_c)                                 ! regular grid version
        ! input / output
        integer, intent(in) :: pt_d                                             ! point on discrete grid
        real(dp), intent(inout) :: pt_c                                         ! point on continous grid
        real(dp), intent(in) :: var_c(:)                                        ! continous grid values
        
        ! Check whether the discrete value lies inside the range
        if (pt_d.lt.1 .or. pt_d.gt.size(var_c)) then
            call writo('WARNING: pt_c not within range',persistent=.true.)
            pt_c = -1._dp
            return
        end if
        
        ! Return the continuous variable
        pt_c = var_c(pt_d)
    end subroutine dis2con_regular
    
    ! rounds  an  arry of  values  to limits,  with  a tolerance  1E-5 that  can
    ! optionally be modified
    integer function round_with_tol_arr(vals,lim_lo,lim_hi,tol) result(ierr)    ! array version
        character(*), parameter :: rout_name = 'round_with_tol_arr'
        
        ! input / output
        real(dp), intent(inout) :: vals(:)                                      ! values to be rounded
        real(dp), intent(in) :: lim_lo, lim_hi                                  ! limits on val
        real(dp), intent(in), optional :: tol                                   ! tolerance
        
        ! local variables
        real(dp) :: loc_tol                                                     ! local value of tol
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set loc_tol
        if (present(tol)) then
            loc_tol = tol
        else
            loc_tol = 1E-5_dp
        end if
        
        ! set possible error message
        err_msg = 'value has to lie within 0..1'
        
        if (minval(vals-lim_lo+loc_tol).lt.0) then
            ierr = 1
            CHCKERR(err_msg)
        else if (maxval(vals-lim_hi-loc_tol).gt.0) then
            ierr = 1
            CHCKERR(err_msg)
        else
            where (vals.lt.lim_lo) vals = lim_lo
            where (vals.gt.lim_hi) vals = lim_hi
        end if
    end function round_with_tol_arr
    integer function round_with_tol_ind(val,lim_lo,lim_hi,tol) result(ierr)     ! individual version
        character(*), parameter :: rout_name = 'round_with_tol_ind'
        
        ! input / output
        real(dp), intent(inout) :: val                                          ! value to be rounded
        real(dp), intent(in) :: lim_lo, lim_hi                                  ! limits on val
        real(dp), intent(in), optional :: tol                                   ! tolerance
        
        ! local variables
        real(dp) :: vals(1)                                                     ! array copy of val
        
        ! initialize ierr
        ierr = 0
        
        vals(1) = val
        
        ierr = round_with_tol_arr(vals,lim_lo,lim_hi,tol)
        CHCKERR('')
        
        val = vals(1)
    end function round_with_tol_ind
    
    ! Interpolate a  1D function  y(x) by  providing the array  y and  x_in. The
    ! result is stored in  y_out. The array x can be  optionally passed. If not,
    ! it is assumed to be the (equidistant) linear space between 0 and 1.
    ! Note: This function is assumed to be monotomous. If not, an error results.
    integer function interp_fun_2D(y_out,y,x_in,x) result(ierr)                 ! 2D version
        character(*), parameter :: rout_name = 'interp_fun_2D'
        
        ! input / output
        real(dp), intent(inout) :: y_out(:,:)                                   ! output y_out
        real(dp), intent(in) :: y(:,:,:)                                        ! y(x)
        real(dp), intent(in) :: x_in                                            ! input x_in
        real(dp), intent(in), optional :: x(:)                                  ! x(x)
        
        ! local variables
        integer :: n_pt                                                         ! nr. points in y
        real(dp) :: ind                                                         ! unrounded x-index
        integer :: ind_lo, ind_hi                                               ! lower and higher index of x_in in x(x)
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_pt
        n_pt = size(y,3)
        
        ! tests
        if (present(x)) then
            if (size(x).ne.n_pt) then
                err_msg = 'x and y need to have the same size'
                ierr = 1
                CHCKERR(err_msg)
            end if
        end if
        
        ! find the lower range of the index in the x array
        if (present(x)) then
            call con2dis(x_in,ind,x)
        else
            call con2dis(x_in,ind,[0._dp,1._dp],[1,n_pt])
        end if
        
        ! test if result has been found
        if (ind.le.0) then
            ierr = 1
            CHCKERR('see WARNING from con2dis above')
        end if
        
        ! set ind_lo and ind_hi
        ind_lo = floor(ind)
        ind_hi = ceiling(ind)
        
        ! calculate y_out
        y_out = y(:,:,ind_lo) + (y(:,:,ind_hi)-y(:,:,ind_lo)) * &
            &(ind-ind_lo)
    end function interp_fun_2D
    integer function interp_fun_1D(y_out,y,x_in,x) result(ierr)                 ! 1D version
        character(*), parameter :: rout_name = 'interp_fun_1D'
        
        ! input / output
        real(dp), intent(inout) :: y_out(:)                                     ! output y_out
        real(dp), intent(in) :: y(:,:)                                          ! y(x)
        real(dp), intent(in) :: x_in                                            ! input x_in
        real(dp), intent(in), optional :: x(:)                                  ! x(x)
        
        ! local variables
        real(dp), allocatable :: y_out_loc(:,:)
        
        ! allocate y_out_loc
        allocate(y_out_loc(1,size(y_out)))
        
        ! call 2D version
        ierr = interp_fun_2D(y_out_loc,reshape(y,[1,size(y_out),size(y)]),&
            &x_in,x)
        CHCKERR('')
        
        ! copy to y_out
        y_out = y_out_loc(1,1)
        
        ! clean up
        deallocate(y_out_loc)
    end function interp_fun_1D
    integer function interp_fun_0D(y_out,y,x_in,x) result(ierr)                 ! 0D version
        character(*), parameter :: rout_name = 'interp_fun_0D'
        
        ! input / output
        real(dp), intent(inout) :: y_out                                        ! output y_out
        real(dp), intent(in) :: y(:)                                            ! y(x)
        real(dp), intent(in) :: x_in                                            ! input x_in
        real(dp), intent(in), optional :: x(:)                                  ! x(x)
        
        ! local variables
        real(dp) :: y_out_loc(1,1)
        
        ! call 2D version
        ierr = interp_fun_2D(y_out_loc,reshape(y,[1,1,size(y)]),x_in,x)
        CHCKERR('')
        
        ! copy to y_out
        y_out = y_out_loc(1,1)
    end function interp_fun_0D
    
    ! convert 2D  coordinates (i,j) to  the storage convention used  in (square)
    ! metric matrices:
    !   (1 4 7)      (1    )
    !   (2 5 8)  or  (2 4  ) for symmetric matrices.
    !   (3 6 9)      (3 5 6)
    ! Optionally,  the  size of  the  (square) matrix  can be  changed from  its
    ! default value of 3x3 using n.
    ! The value of c is given by
    !   c = (j-1)*n + i
    ! for non-symmetric matrices, and by
    !   c = (j-1)*n + i - (j-1)*j/2   if i.ge.j
    !   c = (i-1)*n + j - (i-1)*i/2   if j.ge.i
    ! since sum_k=1^j-1 = (j-1)*j/2
    integer function c(ij,sym,n)
        ! input / output
        integer, intent(in) :: ij(2)                                            ! 2D coords. (i,j)
        logical, intent(in) :: sym                                              ! .true. if symmetric
        integer, intent(in), optional :: n                                      ! max of 2D coords.
        
        ! local variables
        integer :: n_loc
        
        ! set n_loc
        n_loc = 3
        if (present(n)) n_loc = n
        
        ! depending on symmetry
        if (sym) then                                                           ! symmetric
            if (ij(1).ge.ij(2)) then
                c = (ij(2)-1)*n_loc + ij(1) - (ij(2)-1)*ij(2)/2
            else
                c = (ij(1)-1)*n_loc + ij(2) - (ij(1)-1)*ij(1)/2
            end if
        else                                                                    ! asymmetric
            c = (ij(2)-1)*n_loc + ij(1)
        end if
    end function c
    
    ! determines whether a metric matrix is symmetric or not
    integer function is_sym(n,nn,sym) result(ierr)
        character(*), parameter :: rout_name = 'is_sym'
        
        ! input / output
        integer, intent(in) :: n                                                ! size of matrix
        integer, intent(in) :: nn                                               ! number of elements in matrix
        logical, intent(inout) :: sym                                           ! output
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set sym
        if (nn.eq.n**2) then
            sym = .false.
        else if (nn.eq.n*(n+1)/2) then
            sym = .true.
        else
            ierr = 1
            err_msg = 'The matrix corresponds neither to a symmetric nor to a &
                &normal matrix'
            CHCKERR('')
        end if
    end function is_sym

    ! Either takes the  complex conjugate of a square matrix  element A, defined
    ! on a  3D grid,  or not,  depending on  whether the  indices of  the matrix
    ! element  correspond  to the  upper  (conjugate)  or lower  (no  conjugate)
    ! triangular part.
    ! Note: Only for symmetric matrices does this have to be applied.
    function con_3D(A,c,sym,d) result(B)                                        ! 3D version
        complex(dp), intent(in) :: A(:,:,:)                                     ! input A
        integer, intent(in) :: c(2)                                             ! indices in square matrix
        logical, intent(in) :: sym                                              ! if A is symmetric
        integer, intent(in) :: d(3)                                             ! dimensions of matrix A
        complex(dp) :: B(d(1),d(2),d(3))                                        ! output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_3D
    function con_2D(A,c,sym,d) result(B)                                        ! 2D version
        complex(dp), intent(in) :: A(:,:)                                       ! input A
        integer, intent(in) :: c(2)                                             ! indices in square matrix
        logical, intent(in) :: sym                                              ! if A is symmetric
        integer, intent(in) :: d(2)                                             ! dimensions of matrix A
        complex(dp) :: B(d(1),d(2))                                             ! output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_2D
    function con_1D(A,c,sym,d) result(B)                                        ! 1D version
        complex(dp), intent(in) :: A(:)                                         ! input A
        integer, intent(in) :: c(2)                                             ! indices in square matrix
        logical, intent(in) :: sym                                              ! if A is symmetric
        integer, intent(in) :: d(1)                                             ! dimensions of matrix A
        complex(dp) :: B(d(1))                                                  ! output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_1D
    function con_0D(A,c,sym) result(B)                                          ! 0D version
        complex(dp), intent(in) :: A                                            ! input A
        integer, intent(in) :: c(2)                                             ! indices in square matrix
        logical, intent(in) :: sym                                              ! if A is symmetric
        complex(dp) :: B                                                        ! output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_0D
end module utilities
