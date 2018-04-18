!------------------------------------------------------------------------------!
!> Numerical utilities.
!------------------------------------------------------------------------------!
module num_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    
    implicit none
    private
    public calc_ext_var, calc_det, calc_int, add_arr_mult, c, &
        &check_deriv, calc_inv, calc_mult, calc_aux_utilities, derivs, &
        &con2dis, dis2con, round_with_tol, conv_mat, is_sym, con, &
        &calc_coeff_fin_diff, fac, d, m, f, bubble_sort, GCD, LCM, &
        &order_per_fun, shift_F, spline, solve_vand
#if ldebug
    public debug_con2dis_reg, debug_calc_coeff_fin_diff
#endif
    
    ! global variables
    integer, allocatable :: d(:,:,:)                                            !< 1-D array indices of derivatives
    integer, allocatable :: m(:,:)                                              !< 1-D array indices of metric indices
    integer, allocatable :: f(:,:)                                              !< 1-D array indices of Fourier mode combination indices
#if ldebug
    logical :: debug_con2dis_reg = .false.                                      !< plot debug information for con2dis_reg() \ldebug
    logical :: debug_calc_coeff_fin_diff = .false.                              !< plot debug information for calc_coeff_fin_diff() \ldebug
#endif
    
    ! interfaces
    
    !> \public Add to an  array (3) the product  of arrays (1) and  (2).
    !!
    !! The derivatives  are distributed  between both  acording to  the binomial
    !! theorem. Both arrays are given on a 3-D or 1-D grid.
    !!
    !! \return ierr
    interface add_arr_mult
        !> \public
        module procedure add_arr_mult_3_3
        !> \public
        module procedure add_arr_mult_3_1
        !> \public
        module procedure add_arr_mult_1_1
    end interface
    
    !> \public Calculate determinant of a matrix
    !!
    !! This  matrix can  be  defined on  a  3-D grid  or  constant. The  storage
    !! convention described in eq_vars.eq_2_type is used.
    !!
    !! In the former  case the size of  the matrix (last two  indices) should be
    !! small,  as the  direct  formula  employing cofactors  is  used through  a
    !! recursive formulation. 
    !!
    !! In the latter case, lapack routines are used.
    !!
    !! \see Adapted from
    !! <http://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/>
    !!
    !! \return ierr
    interface calc_det
        !> \public
        module procedure calc_det_0D
        !> \public
        module procedure calc_det_3D
    end interface
    
    !> \public  Calculate inverse of square matrix \c A.
    !! 
    !! This  matrix can  be  defined on  a  3-D grid  or  constant. The  storage
    !! convention described in eq_vars.eq_2_type is used.
    !!
    !! In the former  case the size of  the matrix (last two  indices) should be
    !! small, as direct inversion is performed using Cramer's rule.
    !!
    !! In the latter case, lapack routines are used.
    !!
    !! \return ierr
    interface calc_inv
        !> \public
        module procedure calc_inv_0D
        !> \public
        module procedure calc_inv_3D
    end interface
    
    !> \public Calculate matrix multiplication of two square matrices
    !! \f$\overline{\text{AB}} = \overline{\text{A}} \ \overline{\text{B}}\f$.
    !!
    !! This  matrix can  be  defined on  a  3-D grid  or  constant. The  storage
    !! convention described in eq_vars.eq_2_type is used.
    !!
    !! \return ierr
    interface calc_mult
        !> \public
        module procedure calc_mult_0D_real
        !> \public
        module procedure calc_mult_3D_real
        !> \public
        module procedure calc_mult_3D_complex
    end interface
    
    !> \public  Converts a  (symmetric)  matrix A  with  the storage  convention
    !! described in eq_vars.eq_2_type.
    !! 
    !! This matric can have elements depending on a 3-D grid or be constant.
    !!
    !!  - If  the matrix  is stored  with n^2 numbers,  only the  lower diagonal
    !!  elements are kept in matrix B.
    !!  - If only the lower diagonal elements are stored, they are copied to the
    !!  upper diagonal ones of the matrix B as well.
    !!
    !! Optionally, the transpose can be calculated.
    !!
    !! \note 
    !!  -# The routine does not check whether the matrix is indeed symmetric and
    !!  that information may thus be lost after conversion.
    !!  -#  This routine  makes a  copy of  A so  that by  providing A  as input
    !!  argument for both A and B overwrites A.
    !!
    !! \return ierr
    interface conv_mat
        !> \public
        module procedure conv_mat_3D
        !> \public
        module procedure conv_mat_3D_complex
        !> \public
        module procedure conv_mat_0D
        !> \public
        module procedure conv_mat_0D_complex
    end interface
    
    !> \public Integrates a function using the trapezoidal rule.
    !!
    !! This function can be defined on an equidistant grid or a regular one:
    !!  - For an equidistant grid,
    !!      \f[\int_1^n f(x) dx = \sum_{k=1}^{n-1}
    !!          {\left(f(k+1)+f(k)\right) \frac{\Delta_x}{2}},\f]
    !!  is used,  with \f$n\f$  the number  of points, which  are assumed  to be
    !!  equidistant with a given step size \f$Delta_x\f$.
    !!  - For a regular grid,
    !!      \f[\int_1^n f(x) dx = \sum_{k=1}^{n-1}
    !!          {\left(f(k+1)+f(k)\right) \frac{x(k+1)-x(k)}{2}}, \f]
    !!  is used, with \f$n\f$ the number of points.
    !!      - Therefore, \f$n\f$  points have to be specified as well as
    !!      \f$n\f$ values  for the function  to be integrated.
    !!      - They  have to be  given in
    !!      ascending order but the step size does not have to be constant.
    !!      - This yields the following difference formula:
    !!      \f[\int_1^n f(x) dx = \int_1^{n-1} f(x) dx +
    !!          \left(f(n)+f(n-1)\right) \frac{x(n)-x(n-1)}{2},\f]
    !!      which is used here
    !!
    !! \note For periodic function, the trapezoidal  rule works well only if the
    !! last point of the grid is included,  i.e. the point where the function is
    !! equal to the first point.
    !!
    !! \return ierr
    interface calc_int
        !> \public
        module procedure calc_int_eqd
        !> \public
        module procedure calc_int_reg
    end interface
    
    !>  \public  Rounds  an   arry  of  values  to  limits,   with  a  tolerance
    !! \f$10^{-5}\f$ that can optionally be modified.
    interface round_with_tol
        !> \public
        module procedure round_with_tol_ind
        !> \public
        module procedure round_with_tol_arr
    end interface
    
    !> \public Convert between points from a continuous grid to a discrete grid.
    !!
    !! This is done by providing either the the limits on the grid (\c lim_c and
    !! \c lim_d), in  which case the grid  is assumed to be  equidistant, or the
    !! grid values themselves, in which case the grid just has to be regular.
    !!
    !! The output is a real value where  the floored integer is the index in the
    !! discrete grid and  the remainder corresponds to the  fraction towards the
    !! next index.  If no solution is  found, a negative value  is outputted, as
    !! well as a message.
    !!
    !! \return ierr
    interface con2dis
        !> \public
        module procedure con2dis_eqd
        !> \public
        module procedure con2dis_reg
    end interface
    
    !> \public Convert between points from a discrete grid to a continuous grid.
    !!
    !! This is done by providing either the the limits on the grid (\c lim_c and
    !! \c lim_d), in  which case the grid  is assumed to be  equidistant, or the
    !! grid values themselves, in which case the grid just has to be regular.
    !!
    !! The output is a real value. If the discrete value lies outside the range,
    !! in the case of a regular grid,  a negative value is outputted, as well as
    !! a message.
    !!
    !! \return ierr
    interface dis2con
        !> \public
        module procedure dis2con_eqd
        !> \public
        module procedure dis2con_reg
    end interface
    
    !> \public Either takes the complex conjugate  of a square matrix element A,
    !! defined on a 3-D grid, or not.
    !!
    !! This  is done  depending on  whether the  indices of  the matrix  element
    !! correspond to  the upper (conjugate)  or lower (no  conjugate) triangular
    !! part.
    !!
    !! \note Only for symmetric matrices does this have to be applied.
    interface con
        !> \public
        module procedure con_3D
        !> \public
        module procedure con_2D
        !> \public
        module procedure con_1D
        !> \public
        module procedure con_0D
    end interface
    
    !> \public Sorting with the bubble sort routine.
    !!
    !! Optionally, the pivots can be given back.
    !!
    !! \note Adapted from <http://rosettacode.org/wiki/Category:Fortran>
    interface bubble_sort
        !> \public
        module procedure bubble_sort_int
        !> \public
        module procedure bubble_sort_real
    end interface
    
    !> \public Order  a periodic function  to include \f$0\ldots 2\pi\f$  and an
    !! overlap.
    !!
    !! \return ierr
    interface order_per_fun
        !> \public
        module procedure order_per_fun_1
        !> \public
        module procedure order_per_fun_2
    end interface
    
    !> \public  Wrapper to  the pspline  library, making  it easier  to use  for
    !! 1-D  applications  where  speed  is  not the  main  priority.  If  spline
    !! representations are to be reused, manually use the library.
    !!
    !! Order  1 (linear),  2 (akima  hermite)  or 3  (cubic) possible.  Boundary
    !! conditions are possible:
    !!  - -1: periodic
    !!  -  0: not-a-knot
    !!  -  1: prescribe first derivative
    !!  -  2: prescribe second derivative
    !!
    !! However, for order 2 boundary condition  2 is not available and for order
    !! 1 none of them.
    !!
    !! Furthermore, derivatives can be specified:
    !!  - up to 1 for order 1 and 2
    !!  - up to 2 for order 3
    !!
    !! Finally, extrapolation can be performed as well.
    !!
    !! \return ierr
    interface spline
        !> \public
        module procedure spline_real
        !> \public
        module procedure spline_complex
    end interface
    
contains
    !> \private regular version
    integer function calc_int_reg(var,x,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_reg'
        
        ! input / output
        real(dp), intent(inout) :: var_int(:)                                   !< integrated variable
        real(dp), intent(in) :: var(:)                                          !< variable to be integrated
        real(dp), intent(in) :: x(:)                                            !< abscissa
        
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
    end function calc_int_reg
    !> \private equidistant version
    integer function calc_int_eqd(var,step_size,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_eqd'
        
        ! input / output
        real(dp), intent(inout) :: var_int(:)                                   !< integrated variable
        real(dp), intent(in) :: var(:)                                          !< variable to be integrated
        real(dp), intent(in) :: step_size                                       !< step size of abscissa
        
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
    end function calc_int_eqd
    
    !> \private version with \c arr_1 and \c arr_2 in 3 coords.
    integer function add_arr_mult_3_3(arr_1,arr_2,arr_3,deriv) result(ierr)
        character(*), parameter :: rout_name = 'add_arr_mult_3_3'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,1:,0:,0:,0:)                        !< \c arr_1
        real(dp), intent(in) :: arr_2(1:,1:,1:,0:,0:,0:)                        !< \c arr_2
        real(dp), intent(out) :: arr_3(1:,1:,1:)                                !< \c arr_3
        integer, intent(in) :: deriv(3)                                         !< derivatives
        
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
                err_msg = 'Arrays 1 and 2 do not provide deriv. orders for a &
                    &derivative of order ['//trim(i2str(deriv(1)))//','//&
                    &trim(i2str(deriv(2)))//','//trim(i2str(deriv(3)))//'&
                    &] in coordinate '//trim(i2str(kd))
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
    !> \private version with \c arr_1 in 3  coords and \c arr_2 only in the flux
    !! coord.
    integer function add_arr_mult_3_1(arr_1,arr_2,arr_3,deriv) result(ierr)
        character(*), parameter :: rout_name = 'add_arr_mult_3_1'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,1:,1:,0:,0:,0:)                        !< \c arr_1
        real(dp), intent(in) :: arr_2(1:,0:)                                    !< \c arr_2
        real(dp), intent(out) :: arr_3(1:,1:,1:)                                !< \c arr_3
        integer, intent(in) :: deriv(3)                                         !< derivatives
        
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
    !> \private Version with \c arr_1 and \c arr_2 only in the flux coord.
    integer function add_arr_mult_1_1(arr_1,arr_2,arr_3,deriv) result(ierr)
        character(*), parameter :: rout_name = 'add_arr_mult_1_1'
        
        ! input / output
        real(dp), intent(in) :: arr_1(1:,0:)                                    !< \c arr_1
        real(dp), intent(in) :: arr_2(1:,0:)                                    !< \c arr_2
        real(dp), intent(out) :: arr_3(1:)                                      !< \c arr_3
        integer, intent(in) :: deriv(3)                                         !< derivatives
        
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

    !> private array version
    integer recursive function calc_det_3D(detA,A,n) result (ierr)
        character(*), parameter :: rout_name = 'calc_det_3D'
        
        ! input / output
        real(dp), intent(inout) :: detA(:,:,:)                                  !< output
        real(dp), intent(in) :: A(:,:,:,:)                                      !< input
        integer, intent(in) :: n                                                !< size of matrix
        
        ! local variables
        logical :: sym                                                          ! .true. if matrix symmetric
        integer :: id, jd, kd                                                   ! counter
        integer :: nn                                                           ! nr. of elements in matrix
        real(dp) :: sgn                                                         ! holds either plus or minus one
        integer, allocatable :: slct(:)                                         ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: work(:,:,:)                                    ! work array
        integer, allocatable :: c_sub(:)                                        ! indices of submatrix in storage convention of eq_vars.eq_2_type
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
            detA = A(:,:,:,c([1,1],sym,n))*A(:,:,:,c([2,2],sym,n))-&
                &A(:,:,:,c([1,2],sym,n))*A(:,:,:,c([2,1],sym,n))
        else
            ! allocate coordinates of (n-1)x(n-1) submatrix
            allocate(c_sub((n-1)**2))
            allocate(k_sub(n**2)); k_sub = 0
            ! iterate over indices in second dimension
            do id = 1, n
                slct(id) = 0
                ! set up coordinates of submatrix (2:n,pack(idx,slct.gt.0))
                ! (in general not symmetric, even if parent matrix is)
                k_sub = pack(idx,slct.gt.0)
                do kd = 1,n-1
                    do jd = 1,n-1
                        c_sub((kd-1)*(n-1)+jd) = c([jd+1,k_sub(kd)],sym,n)
                    end do
                end do
                ierr = calc_det_3D(work,A(:,:,:,c_sub),n-1)
                CHCKERR('')
                detA = detA + A(:,:,:,c([1,id],sym,n))*sgn*work
                sgn = -sgn
                slct(id) = 1
            end do
            ! deallocate c_sub and k_sub
            deallocate(c_sub,k_sub)
        end if
    end function calc_det_3D
    !> private constant version
    integer function calc_det_0D(det_0D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_det_0D'
        
        ! input / output
        real(dp), intent(inout) :: det_0D                                       !< output
        real(dp), intent(in) :: A(:,:)                                          !< input
        
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
        
        call dgetrf(n,n,A_loc,n,ipiv,info)
        
        det_0D = 1.0_dp
        do id = 1, n
            det_0D = det_0D*A_loc(id,id)
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
    
    !> \private array version
    integer function calc_inv_3D(inv_3D,A,n) result(ierr)
        use num_vars, only: tol_zero
        
        character(*), parameter :: rout_name = 'calc_inv_3D'
        
        ! input / output
        real(dp), intent(inout) :: inv_3D(:,:,:,:)                              !< output
        real(dp), intent(in) :: A(:,:,:,:)                                      !< input
        integer, intent(in) :: n                                                !< size of matrix
        
        ! local variables
        logical :: sym                                                          ! .true. if matrix symmetric
        real(dp), allocatable :: detA(:,:,:)                                    ! determinant of a submatrix of A
        integer :: id, jd, kd, ld                                               ! counters
        integer :: nn                                                           ! nr. of elements in matrix
        integer, allocatable :: slct(:,:)                                       ! 0 in between 1's selects which column to delete
        integer, allocatable :: idx(:)                                          ! counts from 1 to size(A)
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: c_sub(:)                                        ! indices of submatrix in storage convention of eq_vars.eq_2_type
        integer, allocatable :: l_sub(:), k_sub(:)                              ! first and second indices of submatrix in 2D
        integer :: kd_min                                                       ! min. of kd
       
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(inv_3D,1).ne.size(A,1) .or. size(inv_3D,2).ne.size(A,2) .or. &
            &size(inv_3D,3).ne.size(A,3) .or. size(inv_3D,4).ne.size(A,4)) then
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
                ierr = calc_det_3D(detA,A(:,:,:,c_sub),n-1)
                CHCKERR('')
                inv_3D(:,:,:,c([kd,ld],sym,n)) = (-1.0_dp)**(kd+ld)*detA
                slct(ld,1) = 1
                slct(kd,2) = 1
            end do
        end do
        
        ! deallocate c_sub l_sub and k_sub
        deallocate(c_sub,l_sub,k_sub)
        
        ! calculate determinant in detA
        ierr = calc_det(detA,A,n)
        CHCKERR('')
        
        ! limit determinant from being zero
        detA = sign(max(abs(detA),tol_zero),detA)
        
        ! divide by determinant
        do kd = 1,nn
            inv_3D(:,:,:,kd) = inv_3D(:,:,:,kd) / detA
        end do
    end function calc_inv_3D
    !> private constant version
    integer function calc_inv_0D(inv_0D,A) result(ierr)
        character(*), parameter :: rout_name = 'calc_inv_0D'
        
        ! input / output
        real(dp), intent(inout) :: inv_0D(:,:)                                  !< output
        real(dp), intent(in) :: A(:,:)                                          !< input
        
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
        
        call dgetrf(n,n,inv_0D,n,ipiv,ierr)                                     ! LU factorization
        err_msg = 'Lapack couldn''t find the LU factorization'
        CHCKERR(err_msg)
        
        call dgetri(n,inv_0D,n,ipiv,w,n,ierr)                                   ! inverse of LU
        CHCKERR('Lapack couldn''t compute the inverse')
    end function calc_inv_0D
    
    !> \private real array version
    integer function calc_mult_3D_real(A,B,AB,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'calc_mult_3D_real'
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)                                      !< input A
        real(dp), intent(in) :: B(:,:,:,:)                                      !< input B
        real(dp), intent(inout) :: AB(:,:,:,:)                                  !< output A B
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp(2)                              !< whether A and/or B transposed
        
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
    end function calc_mult_3D_real
    !> \private complex array version
    integer function calc_mult_3D_complex(A,B,AB,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'calc_mult_3D_complex'
        
        ! input / output
        complex(dp), intent(in) :: A(:,:,:,:)                                   !< input A
        complex(dp), intent(in) :: B(:,:,:,:)                                   !< input B
        complex(dp), intent(inout) :: AB(:,:,:,:)                               !< output A B
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp(2)                              !< whether A and/or B transposed
        
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
    end function calc_mult_3D_complex
    !> \private real constant version
    integer function calc_mult_0D_real(A,B,AB,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'calc_mult_0D_real'
        
        ! input / output
        real(dp), intent(in) :: A(:)                                            !< input A
        real(dp), intent(in) :: B(:)                                            !< input B
        real(dp), intent(inout) :: AB(:)                                        !< output A B
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp(2)                              !< whether A and/or B transposed
        
        ! local variables
        real(dp), allocatable :: loc_AB(:,:,:,:)                                ! local copy of C
        
        ! initialize ierr
        ierr = 0
        
        ! allocate local AB
        allocate(loc_AB(1,1,1,size(AB)))
        
        ! calculate using 2D version
        ierr = calc_mult_3D_real(reshape(A,[1,1,1,size(A)]),&
            &reshape(B,[1,1,1,size(B)]),loc_AB,n,transp)
        
        ! update AB with local value
        AB = loc_AB(1,1,1,:)
        
        ! deallocate local variables
        deallocate(loc_AB)
    end function calc_mult_0D_real
    
    !> \private version defined on 3-D grid
    integer function conv_mat_3D(A,B,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'conv_mat_3D'
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)                                      !< matrix to be converted
        real(dp), intent(inout) :: B(:,:,:,:)                                   !< converted matrix
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp                                 !< transpose
        
        ! local variables
        logical :: sym(2)                                                       ! .true. if matrices symmetric
        logical :: transp_loc                                                   ! local transp
        integer :: nn(2)                                                        ! nr. of elements in matrices
        integer :: id, jd                                                       ! counters
        integer :: id_min                                                       ! minimum value of id
        integer :: ind(2,2)                                                     ! indices
        real(dp), allocatable :: A_loc(:,:,:,:)                                 ! local copy of A
        
        ! initialize ierr
        ierr = 0
        
        ! set total number of elements in matrix
        nn(1) = size(A,4)
        nn(2) = size(B,4)
        
        ! set local transp
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        
        ! set local A
        allocate(A_loc(size(A,1),size(A,2),size(A,3),size(A,4)))
        A_loc = A
        
        ! set sym
        do id = 1,2
            ierr = is_sym(n,nn(id),sym(id))
            CHCKERR('')
        end do
        
        ! copy A into B
        if (sym(1).and.sym(2) .or. &
            &.not.sym(1).and..not.sym(2).and..not.transp_loc) then              ! both symmetric or both unsymmetric and no transpose needed
            B = A                                                               ! symmetric matrix equal to transpose or no tranpose needed
        else                                                                    ! at least one is not symmetric
            ! loop over row
            do jd = 1,n
                ! set id_min
                if (sym(2)) then
                    id_min = jd
                else
                    id_min = 1
                end if
                ! loop over column
                do id = id_min,n
                    ! set ind
                    if (transp_loc) then
                        ind(:,1) = [jd,id]
                    else
                        ind(:,1) = [id,jd]
                    end if
                    ind(:,2) = [id,jd]
                    ! copy elements of A into B
                    B(:,:,:,c(ind(:,2),sym(2),n)) = &
                        &A_loc(:,:,:,c(ind(:,1),sym(1),n))
                end do
            end do
        end if
    end function conv_mat_3D
    !> \private scalar version
    integer function conv_mat_0D(A,B,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'conv_mat_0D'
        
        ! input / output
        real(dp), intent(in) :: A(:)                                            !< matrix to be converted
        real(dp), intent(inout) :: B(:)                                         !< converted matrix
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp                                 !< transpose
        
        ! local variables
        real(dp), allocatable :: B_loc(:,:,:,:)                                 ! local B
        
        ! initialize ierr
        ierr = 0
        
        ! set up local B
        allocate(B_loc(1,1,1,size(B)))
        
        ! call 3D version
        ierr = conv_mat_3D(reshape(A,[1,1,1,size(A)]),B_loc,n,transp)
        CHCKERR('')
        
        ! copy to B
        B = B_loc(1,1,1,:)
    end function conv_mat_0D
    !> \private complex version defined on 3D grid
    integer function conv_mat_3D_complex(A,B,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'conv_mat_3D_complex'
        
        ! input / output
        complex(dp), intent(in) :: A(:,:,:,:)                                   !< matrix to be converted
        complex(dp), intent(inout) :: B(:,:,:,:)                                !< converted matrix
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp                                 !< transpose
        
        ! local variables
        logical :: sym(2)                                                       ! .true. if matrices symmetric
        logical :: transp_loc                                                   ! local transp
        integer :: nn(2)                                                        ! nr. of elements in matrices
        integer :: id, jd                                                       ! counters
        integer :: id_min                                                       ! minimum value of id
        integer :: ind(2,2)                                                     ! indices
        complex(dp), allocatable :: A_loc(:,:,:,:)                              ! local copy of A
        
        ! initialize ierr
        ierr = 0
        
        ! set total number of elements in matrix
        nn(1) = size(A,4)
        nn(2) = size(B,4)
        
        ! set local transp
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        
        ! set local A
        allocate(A_loc(size(A,1),size(A,2),size(A,3),size(A,4)))
        A_loc = A
        
        ! set sym
        do id = 1,2
            ierr = is_sym(n,nn(id),sym(id))
            CHCKERR('')
        end do
        
        ! copy A into B
        if (sym(1).and.sym(2) .or. &
            &.not.sym(1).and..not.sym(2).and..not.transp_loc) then              ! both symmetric or both unsymmetric and no transpose needed
            B = A                                                               ! symmetric matrix equal to transpose or no tranpose needed
        else                                                                    ! at least one is not symmetric
            ! loop over row
            do jd = 1,n
                ! set id_min
                if (sym(2)) then
                    id_min = jd
                else
                    id_min = 1
                end if
                ! loop over column
                do id = id_min,n
                    ! set ind
                    if (transp_loc) then
                        ind(:,1) = [jd,id]
                    else
                        ind(:,1) = [id,jd]
                    end if
                    ind(:,2) = [id,jd]
                    ! copy elements of A into B
                    B(:,:,:,c(ind(:,2),sym(2),n)) = &
                        &A_loc(:,:,:,c(ind(:,1),sym(1),n))
                end do
            end do
        end if
    end function conv_mat_3D_complex
    !> \private complex scalar version
    integer function conv_mat_0D_complex(A,B,n,transp) result(ierr)
        character(*), parameter :: rout_name = 'conv_mat_0D_complex'
        
        ! input / output
        complex(dp), intent(in) :: A(:)                                         !< matrix to be converted
        complex(dp), intent(inout) :: B(:)                                      !< converted matrix
        integer, intent(in) :: n                                                !< size of matrix
        logical, intent(in), optional :: transp                                 !< transpose
        
        ! local variables
        complex(dp), allocatable :: B_loc(:,:,:,:)                              ! local B
        
        ! initialize ierr
        ierr = 0
        
        ! set up local B
        allocate(B_loc(1,1,1,size(B)))
        
        ! call 3D version
        ierr = conv_mat_3D_complex(reshape(A,[1,1,1,size(A)]),B_loc,n,transp)
        CHCKERR('')
        
        ! copy to B
        B = B_loc(1,1,1,:)
    end function conv_mat_0D_complex

    !> \private equidistant version
    integer function con2dis_eqd(pt_c,pt_d,lim_c,lim_d) result(ierr)
        !character(*), parameter :: rout_name = 'con2dis_eqd'
        
        ! input / output
        real(dp), intent(in) :: pt_c                                            !< point on continous grid
        real(dp), intent(inout) :: pt_d                                         !< point on discrete grid
        real(dp), intent(in) :: lim_c(2)                                        !< <tt>[min_c,max_c]</tt>
        integer, intent(in) :: lim_d(2)                                         !< <tt>[min_d,max_d]</tt>
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate, using the formula:
        !    pt_c - min_c     pt_d - min_d
        !   ------------- =  ------------- ,
        !   max_c - min_c    max_d - min_d
        pt_d = lim_d(1) + (lim_d(2)-lim_d(1)) * (pt_c-lim_c(1)) / &
            &(lim_c(2)-lim_c(1))
    end function con2dis_eqd
    !> \private regular version
    integer function con2dis_reg(pt_c,pt_d,var_c) result(ierr)
        character(*), parameter :: rout_name = 'con2dis_reg'
        
        ! input / output
        real(dp), intent(in) :: pt_c                                            !< point on continous grid
        real(dp), intent(inout) :: pt_d                                         !< point on discrete grid
        real(dp), intent(in) :: var_c(:)                                        !< continous grid values
        
        ! local variables
        real(dp), allocatable :: var_c_loc(:)                                   ! local copy of var
        real(dp), allocatable :: var_c_inv(:)                                   ! inverted var_c
        integer :: size_c                                                       ! size of var_c
        integer :: ind_lo, ind_hi                                               ! lower and upper index comprising pt_c
        integer :: id                                                           ! counter
        real(dp) :: pt_c_loc                                                    ! local pt_c
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), parameter :: tol = 1.E-8_dp                                   ! tolerance for comparisons
        
        ! initialize ierr
        ierr = 0
        
        ! set up size_c
        size_c = size(var_c)
        
#if ldebug
        if (debug_con2dis_reg) &
            &call writo('finding the discrete index of '//trim(r2str(pt_c)))
#endif
        
        ! set up local var_c and pt_c_loc
        allocate(var_c_loc(size_c))
        if (var_c(1).lt.var_c(size_c)) then
            var_c_loc = var_c
            pt_c_loc = pt_c
        else
#if ldebug
            if (debug_con2dis_reg) &
                &call writo('setting var_c_loc = - var_c and pt_c_loc = -pt_c')
#endif
            var_c_loc = - var_c                                                 ! local var should be rising
            pt_c_loc = - pt_c                                                   ! invert pt_c as well
        end if
        
        ! set up inverse var_c
        allocate(var_c_inv(size_c))
        var_c_inv = var_c_loc(size_c:1:-1)
        
        ! find the lower and upper index comprising local pt_c
        ind_lo = 0
        ind_hi = size_c+1
        do id = 1,size_c
            if (pt_c_loc+tol*abs(pt_c_loc) .ge. var_c_loc(id)) then
                ind_lo = id
#if ldebug 
                if (debug_con2dis_reg) call writo('for iteration '//&
                    &trim(i2str(id))//'/'//trim(i2str(size_c))//', '//&
                    &trim(r2strt(pt_c_loc))//' + '//trim(r2strt(tol))//&
                    &'*abs('//trim(r2strt(pt_c_loc))//') >= '//&
                    &trim(r2strt(var_c_loc(id)))//' => ['//&
                    &trim(i2str(ind_lo))//','//trim(i2str(ind_hi))//']')
#endif
            end if
            if (pt_c_loc-tol*abs(pt_c_loc) .le. var_c_inv(id)) then
                ind_hi = size_c+1-id
#if ldebug 
                if (debug_con2dis_reg) call writo('for iteration '//&
                    &trim(i2str(id))//'/'//trim(i2str(size_c))//', '//&
                    &trim(r2strt(pt_c_loc))//' - '//trim(r2strt(tol))//&
                    &'*abs('//trim(r2strt(pt_c_loc))//') < '//&
                    &trim(r2strt(var_c_inv(id)))//' => ['//&
                    &trim(i2str(ind_lo))//','//trim(i2str(ind_hi))//']')
#endif
            end if
        end do
#if ldebug
        if (debug_con2dis_reg) then
            call writo('final ind_lo = '//trim(i2str(ind_lo))//', ind_hi = '//&
                &trim(i2str(ind_hi)))
            !!call print_ex_2D('var_c_loc, var_c_inv','',&
                !!&reshape([var_c_loc,var_c_inv],[size_c,2]))
        end if
#endif
        
        ! tests
        if (ind_lo.eq.0 .or. ind_hi.eq.size_c+1) then                           ! not within range
            pt_d = -1._dp
            ierr = 1
            err_msg = 'pt_c '//trim(r2str(pt_c))//' not within range '//&
                &trim(r2str(minval(var_c)))//'..'//trim(r2str(maxval(var_c)))
            CHCKERR(err_msg)
        end if
        
        ! set output
        if (ind_lo.lt.ind_hi) then                                              ! valid output that does not correspond to a point on grid
            pt_d = ind_lo + (pt_c_loc-var_c_loc(ind_lo))/&
                &(var_c_loc(ind_hi)-var_c_loc(ind_lo))
        else if (ind_lo.eq.ind_hi) then                                         ! valid output that does correspond to a point on grid
            pt_d = ind_lo
        else                                                                    ! invalid output
            pt_d = -1._dp
            ierr = 1
            err_msg = 'ind_lo cannot be higher than ind_hi'
            CHCKERR(err_msg)
        end if
        
        ! deallocate
        deallocate(var_c_loc,var_c_inv)
    end function con2dis_reg
    
    !> \private equidistant version
    integer function dis2con_eqd(pt_d,pt_c,lim_d,lim_c) result(ierr)
        !character(*), parameter :: rout_name = 'dis2con_eqd'
        
        ! input / output
        integer, intent(in) :: pt_d                                             !< point on discrete grid
        real(dp), intent(inout) :: pt_c                                         !< point on continous grid
        integer, intent(in) :: lim_d(2)                                         !< <tt>[min_d,max_d]</tt>
        real(dp), intent(in) :: lim_c(2)                                        !< <tt>[min_c,max_c]</tt>
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate, using the formula:
        !    pt_c - min_c     pt_d - min_d
        !   ------------- =  ------------- ,
        !   max_c - min_c    max_d - min_d
        pt_c = lim_c(1) + (lim_c(2)-lim_c(1)) * (pt_d-lim_d(1)) / &
            &(lim_d(2)-lim_d(1))
    end function dis2con_eqd
    !> \private regular version
    integer function dis2con_reg(pt_d,pt_c,var_c) result(ierr)
        character(*), parameter :: rout_name = 'dis2con_reg'
        
        ! input / output
        integer, intent(in) :: pt_d                                             !< point on discrete grid
        real(dp), intent(inout) :: pt_c                                         !< point on continous grid
        real(dp), intent(in) :: var_c(:)                                        !< continous grid values
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! Check whether the discrete value lies inside the range
        if (pt_d.lt.1 .or. pt_d.gt.size(var_c)) then
            pt_c = -1._dp
            ierr = 1
            err_msg = 'pt_c not within range'
            CHCKERR(err_msg)
        end if
        
        ! Return the continuous variable
        pt_c = var_c(pt_d)
    end function dis2con_reg
    
    !> \private array version
    integer function round_with_tol_arr(vals,lim_lo,lim_hi,tol) result(ierr)
        character(*), parameter :: rout_name = 'round_with_tol_arr'
        
        ! input / output
        real(dp), intent(inout) :: vals(:)                                      !< values to be rounded
        real(dp), intent(in) :: lim_lo                                          !< lower limit on val
        real(dp), intent(in) :: lim_hi                                          !< upper limit on val
        real(dp), intent(in), optional :: tol                                   !< tolerance
        
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
    !> \private individual version
    integer function round_with_tol_ind(val,lim_lo,lim_hi,tol) result(ierr)
        character(*), parameter :: rout_name = 'round_with_tol_ind'
        
        ! input / output
        real(dp), intent(inout) :: val                                          !< value to be rounded
        real(dp), intent(in) :: lim_lo                                          !< lower limit on val
        real(dp), intent(in) :: lim_hi                                          !< upper limit on val
        real(dp), intent(in), optional :: tol                                   !< tolerance
        
        ! local variables
        real(dp) :: vals(1)                                                     ! array copy of val
        
        ! initialize ierr
        ierr = 0
        
        vals(1) = val
        
        ierr = round_with_tol_arr(vals,lim_lo,lim_hi,tol)
        CHCKERR('')
        
        val = vals(1)
    end function round_with_tol_ind
    
    !> \private 3-D version
    function con_3D(A,c,sym,d) result(B)
        complex(dp), intent(in) :: A(:,:,:)                                     !< input A
        integer, intent(in) :: c(2)                                             !< indices in square matrix
        logical, intent(in) :: sym                                              !< if A is symmetric
        integer, intent(in) :: d(3)                                             !< dimensions of matrix A
        complex(dp) :: B(d(1),d(2),d(3))                                        !< output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_3D
    !> \private 2-D version
    function con_2D(A,c,sym,d) result(B)
        complex(dp), intent(in) :: A(:,:)                                       !< input A
        integer, intent(in) :: c(2)                                             !< indices in square matrix
        logical, intent(in) :: sym                                              !< if A is symmetric
        integer, intent(in) :: d(2)                                             !< dimensions of matrix A
        complex(dp) :: B(d(1),d(2))                                             !< output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_2D
    !> \private 1-D version
    function con_1D(A,c,sym,d) result(B)
        complex(dp), intent(in) :: A(:)                                         !< input A
        integer, intent(in) :: c(2)                                             !< indices in square matrix
        logical, intent(in) :: sym                                              !< if A is symmetric
        integer, intent(in) :: d(1)                                             !< dimensions of matrix A
        complex(dp) :: B(d(1))                                                  !< output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_1D
    !> \private 0-D version
    function con_0D(A,c,sym) result(B)
        complex(dp), intent(in) :: A                                            !< input A
        integer, intent(in) :: c(2)                                             !< indices in square matrix
        logical, intent(in) :: sym                                              !< if A is symmetric
        complex(dp) :: B                                                        !< output B
        
        ! determine whether to take complex conjugate or not
        if (c(2).gt.c(1) .and. sym) then                                        ! upper triangular matrix
            B = conjg(A)
        else
            B = A
        end if
    end function con_0D
    
    !> \private real version
    subroutine bubble_sort_real(a,piv)
        ! input / output
        real(dp), intent(inout) :: a(:)                                         !< vector to sort
        integer, intent(inout), optional :: piv(:)                              !< pivots
        
        ! local variables
        real(dp) :: temp                                                        ! temporary value
        integer :: id, jd                                                       ! counter
        logical :: swapped                                                      ! value swapped with next
        logical :: save_piv                                                     ! save pivots
        
        save_piv = .false.
        if (present(piv)) then
            if (size(piv).eq.size(a)) then
                save_piv = .true.
                piv = [(jd,jd=1,size(a))]
            end if
        end if
        
        do jd = size(a)-1,1,-1
            swapped = .false.
            do id = 1, jd
                if (a(id) .gt. a(id+1)) then
                    temp = a(id)
                    a(id) = a(id+1)
                    a(id+1) = temp
                    piv(id:id+1) = piv(id+1:id:-1)
                    swapped = .true.
                end if
            end do
            if (.not. swapped) exit
        end do
    end subroutine bubble_sort_real
    !> \private integer version
    subroutine bubble_sort_int(a,piv)
        ! input / output
        integer, intent(inout) :: a(:)                                          !< vector to sort
        integer, intent(inout), optional :: piv(:)                              !< pivots
        
        ! local variables
        integer :: temp                                                         ! temporary value
        integer :: id, jd                                                       ! counter
        logical :: swapped                                                      ! value swapped with next
        logical :: save_piv                                                     ! save pivots
        
        save_piv = .false.
        if (present(piv)) then
            if (size(piv).eq.size(a)) then
                save_piv = .true.
                piv = [(jd,jd=1,size(a))]
            end if
        end if
        
        do jd = size(a)-1,1,-1
            swapped = .false.
            do id = 1, jd
                if (a(id) .gt. a(id+1)) then
                    temp = a(id)
                    a(id) = a(id+1)
                    a(id+1) = temp
                    if (present(piv)) piv(id:id+1) = piv(id+1:id:-1)
                    swapped = .true.
                end if
            end do
            if (.not. swapped) exit
        end do
    end subroutine bubble_sort_int
    
    !> \private version with \c x and \c y separate
    integer function order_per_fun_1(x,y,x_out,y_out,overlap,tol) result(ierr)
        use num_vars, only: tol_zero
        
        character(*), parameter :: rout_name = 'order_per_fun_1'
        
        ! input / output
        real(dp), intent(in) :: x(:)                                            !< abscissa
        real(dp), intent(in) :: y(:)                                            !< ordinate
        real(dp), intent(inout), allocatable :: x_out(:)                        !< ordered x
        real(dp), intent(inout), allocatable :: y_out(:)                        !< ordered y
        integer, intent(in) :: overlap                                          !< overlap to include
        real(dp), intent(in), optional :: tol                                   !< tolerance for error
        
        ! local variables
        real(dp) :: tol_loc                                                     ! local tolerance
        real(dp), allocatable :: x_fund(:)                                      ! x on fundamental interval 0..2pi
        integer :: ml_x                                                         ! location of maximum of x
        integer :: lim(2)                                                       ! limits in x
        integer :: lim_loc(2)                                                   ! limits in loc_x
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set local tol
        tol_loc = tol_zero
        if (present(tol)) tol_loc = tol
        
        ! tests
        if (abs(cos(x(size(x)))-cos(x(1))).lt.tol_loc .and. &
            &abs(sin(x(size(x)))-sin(x(1))).lt.tol_loc) then
            ierr = 1
            err_msg = 'First and last point cannot be identical. Remove one.'
            CHCKERR(err_msg)
        end if
        
        ! initialize
        allocate(x_fund(size(x)))
        x_fund = mod(x,2*pi)
        where(x_fund.lt.0._dp) x_fund = x_fund+2*pi
        ml_x = maxloc(x_fund,1)
        
        allocate(x_out(size(x)+2*overlap))
        allocate(y_out(size(x)+2*overlap))
        
        ! left overlap
        lim = [max(ml_x-overlap+1,1),ml_x]                                      ! all possible points to the left of ml_x, u
        lim_loc = [overlap-lim(2)+lim(1),overlap]
        x_out(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))-2*pi
        y_out(lim_loc(1):lim_loc(2)) = y(lim(1):lim(2))
        lim_loc = [1,lim_loc(1)-1]
        lim = [size(x)-(lim_loc(2)-lim_loc(1)),size(x)]
        x_out(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))-2*pi
        y_out(lim_loc(1):lim_loc(2)) = y(lim(1):lim(2))
        
        ! bulk
        lim = [ml_x+1,size(x)]                                                  ! points after max, until last
        lim_loc = [overlap+1,overlap+lim(2)-lim(1)+1]
        x_out(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))
        y_out(lim_loc(1):lim_loc(2)) = y(lim(1):lim(2))
        lim = [1,ml_x]                                                          ! points from start, until max
        lim_loc = lim_loc(2) + [1,lim(2)-lim(1)+1]
        x_out(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))
        y_out(lim_loc(1):lim_loc(2)) = y(lim(1):lim(2))
        
        ! right overlap
        lim = [ml_x+1,min(ml_X+overlap,size(x))]                                ! all possible poinst to the right of ml_X
        lim_loc = lim_loc(2) + [1,lim(2)-lim(1)+1]
        x_out(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))+2*pi
        y_out(lim_loc(1):lim_loc(2)) = y(lim(1):lim(2))
        lim_loc = [lim_loc(2)+1,size(x)+2*overlap]                              ! remaining point loop around
        lim = [1,lim_loc(2)-lim_loc(1)+1]
        x_out(lim_loc(1):lim_loc(2)) = x_fund(lim(1):lim(2))+2*pi
        y_out(lim_loc(1):lim_loc(2)) = y(lim(1):lim(2))
    end function order_per_fun_1
    !> \private version with \c x and \c y together
    integer function order_per_fun_2(xy,xy_out,overlap,tol) result(ierr)
        character(*), parameter :: rout_name = 'order_per_fun_2'
        
        ! input / output
        real(dp), intent(in) :: xy(:,:)                                         !< abscissa and ordinate
        real(dp), intent(inout), allocatable :: xy_out(:,:)                     !< ordered xy
        integer, intent(in) :: overlap                                          !< overlap to include
        real(dp), intent(in), optional :: tol                                   !< tolerance for error
        
        ! local variables
        real(dp), allocatable :: x_out(:)                                       ! local x_out
        real(dp), allocatable :: y_out(:)                                       ! local y_out
        
        ! initialize ierr
        ierr = 0
        
        ierr = order_per_fun_1(xy(:,1),xy(:,2),x_out,y_out,overlap,tol)
        CHCKERR('')
        
        allocate(xy_out(size(x_out),2))
        xy_out(:,1) = x_out
        xy_out(:,2) = y_out
    end function order_per_fun_2
    
    !> \private real version
    integer function spline_real(x,y,xnew,ynew,ord,deriv,bcs,bcs_val,extrap) &
        &result(ierr)
        
        use EZspline_obj
        use EZspline
        
        character(*), parameter :: rout_name = 'spline_real'
        
        ! input / output
        real(dp), intent(in) :: x(:)                                            !< coordinates
        real(dp), intent(in) :: y(:)                                            !< function value
        real(dp), intent(in) :: xnew(:)                                         !< new coordinates
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
        real(dp) :: bcs_val_loc(2)                                              ! local bcs_val
        real(dp) :: lim_vals(0:3)                                               ! limit values at boundaries for extrapolation
        character(len=max_str_ln) :: err_msg                                    ! error message
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
        
#if ldebug
        ! check monotony
        do kd = 2,n
            if (x(kd).le.x(kd-1)) then
                ierr = 1
                err_msg = 'x is not monotonously increasing'
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
            if (xnew(kd).ge.x(1)) exit
        end do
        il(1) = kd
        do kd = nnew,1,-1
            if (xnew(kd).le.x(n)) exit
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
        f_spl%x1 = x
        
        ! set boundary condition
        if (present(bcs_val)) then
            f_spl%bcval1min = bcs_val_loc(1)
            f_spl%bcval1max = bcs_val_loc(2)
        end if
        
        ! set up 
        call EZspline_setup(f_spl,y,ierr,exact_dim=.true.)                      ! match exact dimensions
        call EZspline_error(ierr)
        CHCKERR('')
        
        ! interpolated part
        if (il(1).le.il(2)) then
            if (deriv_loc.eq.0) then
                ! interpolate
                call EZspline_interp(f_spl,nnew_interp,xnew(il(1):il(2)),&
                    &ynew(il(1):il(2)),ierr)
                call EZspline_error(ierr)
                CHCKERR('')
            else
                call EZspline_derivative(f_spl,deriv_loc,nnew_interp,&
                    &xnew(il(1):il(2)),ynew(il(1):il(2)),ierr) 
                call EZspline_error(ierr)
                CHCKERR('')
            end if
        end if
        
        ! extrapolated part
        if (extrap_loc) then
            ! extrapolate left
            if (il(1).gt.1) then
                ! set up limit values at first point
                ierr = setup_lim_vals(x(1),lim_vals)
                CHCKERR('')
                
                ! calculate extrapolation
                call calc_extrap(xnew(1:il(1)-1),x(1),lim_vals,ynew(1:il(1)-1))
            end if
            
            ! extrapolate right
            if (il(2).lt.nnew) then
                ! set up limit values at last point
                ierr = setup_lim_vals(x(n),lim_vals)
                CHCKERR('')
                
                ! calculate extrapolation
                call calc_extrap(xnew(il(2)+1:nnew),x(n),lim_vals,&
                    &ynew(il(2)+1:nnew))
            end if
        end if
        
        ! free
        call EZspline_free(f_spl,ierr)
        call EZspline_error(ierr)
        CHCKERR('')
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
    
    !> \public  Initialize utilities  for  fast future  reference, depending  on
    !! program style.
    !!
    !! Utilities initialized:
    !!  - derivatives
    !!  - metrics
    !!  - Fourier modes (optionally)
    !!
    !! If Fourier  modes are also initialized,  the quantity \c n_mod  has to be
    !! provided as well.
    subroutine calc_aux_utilities(n_mod)
        use num_vars, only: max_deriv
        
        ! input / output
        integer, intent(in), optional :: n_mod                                  !< n_mod for Fourier modes
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        
        ! common for all program stles
        
        ! derivatives d from 0 to max_deriv+1
        allocate(d(0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        do kd = 0,max_deriv+1
            do jd = 0,max_deriv+1
                do id = 0,max_deriv+1
                    if (id+jd+kd.le.max_deriv+1) then                           ! valid derivatives
                        d(id,jd,kd) = calc_derivs_1D_id([id,jd,kd],3)
                    else                                                        ! derivatives too high
                        d(id,jd,kd) = 0
                    end if
                end do
            end do
        end do
        
        ! metrics m from 1 to 3
        allocate(m(1:3,1:3))
        do jd = 1,3
            do id = 1,3
                m(id,jd) = calc_derivs_1D_id([id,jd],2)
            end do
        end do
        
        ! Fourier modes from 1 to n_mod
        if (present(n_mod)) then
            allocate(f(1:n_mod,1:n_mod))
            do jd = 1,n_mod
                do id = 1,n_mod
                    f(id,jd) = calc_derivs_1D_id([id,jd],2)
                end do
            end do
        end if
    end subroutine calc_aux_utilities
    
    !> Returns the  1D indices  for derivatives  of a  certain order  in certain
    !! dimensions.
    !!
    !! The  algorithm  works by  considering a  structure such  as the  following
    !! example in three dimensions:
    !!  1. (0,0,0)
    !!  2. (1,0,0)
    !!  3. (0,1,0)
    !!  4. (0,0,1)
    !!  5. (2,0,0)
    !!  6. (1,1,0)
    !!  7. (1,0,1)
    !!  8. (0,2,0)
    !!  9. (0,1,1)
    !!  10. (0,0,2)
    !!  11. (3,0,0)
    !!  12. (2,1,0)
    !!
    !! etc...
    !!
    !! By then  defining \f$d\f$ as  the vector  of derivatives, \f$n\f$  as the
    !! size  of  \f$d\f$, \f$D  =  \sum_i^n  d(i)\f$  the  total degree  of  the
    !! derivative,  and \f$I\f$  as the  index of  the last  nonzero element  in
    !! \f$d\f$, and extending \f$d\f$ to the left by considering \f$d(0)\f$ = 0,
    !! the following formula can be can be deduced for the displacements in this
    !! table with respect to index 1:
    !! 
    !! displacements in this table with respect to index 1:
    !!  \f[\sum_{j=0}^{I-1} \sum_{i=0}^
    !!      {D-(d(0)+...+d(j))-1}
    !!      \left(\begin{array}{c}n-j+i-1\\i\end{array}\right) ,
    !!  \f]
    !! making use of the binomial coefficients
    !!  \f[\left(\begin{array}{c}a\\b\end{array}\right) = \frac{a!}{b!(a-b)!}\f]
    !!
    !! It can be
    !! seen that  each of the terms  in the summation in  \f$j\f$ corresponds to
    !! the displacement in dimension \f$j\f$  and the binomial coefficient comes
    !! from the classic stars and bars problem.
    integer function calc_derivs_1D_id(deriv,dims) result(res)
        ! input / output
        integer, intent(in) :: deriv(:)                                         !< derivatives
        integer, intent(in) :: dims                                             !< nr. of dimensions
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: id_max                                                       ! index of last nonzero element in deriv
        integer :: tot_deriv                                                    ! total derivative
        integer, allocatable :: deriv_loc(:)                                    ! local extended deriv
        
        ! set up local deriv
        allocate(deriv_loc(0:size(deriv)))
        deriv_loc(0) = 0
        deriv_loc(1:size(deriv)) = deriv
        
        ! set up total derivative
        tot_deriv = sum(deriv)
        
        ! initialize res
        res = 1
        
        ! calculate  displacement  if total  deriv not  equal to  zero (assuming
        ! deriv > 0)
        if (tot_deriv.gt.0) then
            ! get index of last nonzero element in deriv
            id_max = 0
            do id = 1,size(deriv)
                if (deriv(id).gt.0) id_max = id
            end do
            
            ! iterate over dimensions
            do jd = 0,id_max-1
                ! iterate over derivatives of dimension jd (where deriv(jd) = 0)
                do id = 0,tot_deriv-sum(deriv_loc(0:jd))-1
                    res = res + fac(dims-jd+id-1)/(fac(id)*fac(dims-jd-1))
                end do
            end do
        end if
    end function calc_derivs_1D_id
    
    !> Calculate  all of  the  unique  derivatives in  certain  dimensions of  a
    !! certain total order.
    !!
    !! For example, for order 2 in 3 dimensions this would be:
    !!  - (2,0,0)
    !!  - (1,1,0)
    !!  - (1,0,1)
    !!  - (0,2,0)
    !!  - (0,1,1) and
    !!  - (0,0,2)
    !! 
    !! The total  number of  derivatives  is found  through the  stars and  bars
    !! problem to be the binomial coefficient 
    !!   \f[\left(\begin{array}{c}order+dims-1\\order\end{array}\right) =
    !!      \frac{(order+dims-1)!}{order!(dims-1)!} . \f]
    !! The number  of dimensions is taken to  be 3 by default but  can be changed
    !! optionally.
    recursive subroutine calc_derivs_of_order(res,order,dims)
        ! input / output
        integer, intent(inout), allocatable :: res(:,:)                         !< array of all unique derivatives
        integer, intent(in) :: order                                            !< order of derivative
        integer, intent(in), optional :: dims                                   !< nr. of dimensions
        
        ! local variables
        integer :: id                                                           ! counter
        integer, allocatable :: derivs_loc(:,:)                                 ! derivs of local subblock
        integer :: id_tot                                                       ! counter for total dimension of local derivs.
        integer :: dims_loc                                                     ! local version of dims
        
        ! set up local dims
        dims_loc = 3
        if (present(dims)) dims_loc = dims
        
        ! allocate resulting derivs
        if (allocated(res)) deallocate(res)
        allocate(res(&
            &dims_loc,fac(order+dims_loc-1)/(fac(order)*fac(dims_loc-1))))
        
        ! iterate over the first index if order greater than 0
        if (order.gt.0) then
            id_tot = 0
            do id = order,0,-1
                ! calculate  the  derivs  of  local  subblock  if  more  than  1
                ! dimension
                if (dims_loc.gt.1) then
                    call calc_derivs_of_order(derivs_loc,order-id,dims_loc-1)   ! calculate iteratively subblock
                    res(1,id_tot+1:id_tot+size(derivs_loc,2)) = id              ! set first dimension
                    res(2:dims_loc,id_tot+1:id_tot+size(derivs_loc,2)) &
                        &= derivs_loc                                           ! set next dimensions
                    id_tot = id_tot+size(derivs_loc,2)
                else 
                    res = order
                end if
            end do
        else
            res = 0
        end if
    end subroutine calc_derivs_of_order
    !> Returns derivatives of certain order.
    function derivs(order,dims) result(res)
        ! input / output
        integer, intent(in) :: order                                            !< order of derivative
        integer, intent(in), optional :: dims                                   !< nr. of dimensions
        integer, allocatable :: res(:,:)                                        !< array of all unique derivatives
        
        ! call the subroutine
        call calc_derivs_of_order(res,order,dims)
    end function derivs
    
    !> checks whether  the derivatives  requested for  a certain  subroutine are
    !! valid
    !!
    !! \return ierr
    integer function check_deriv(deriv,max_deriv,sr_name) result(ierr)
        character(*), parameter :: rout_name = 'check_deriv'
        
        ! input / output
        integer, intent(in) :: deriv(3)                                         !< derivative to be checked
        integer, intent(in) :: max_deriv                                        !< maximum derivative allowed
        character(len=*), intent(in) :: sr_name                                 !< name of subroutine
        
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

    !> Extrapolates a function.
    !!
    !! This is  done using linear  or quadratic interpolation, depending  on the
    !! number of points and values given.
    !!
    !! The data  should be given  sorted, in ascending order,  without duplicate
    !! points in var_points.
    !!
    !! It uses the following solution:
    !!  \f[\overline{\text{A}} \ \vec{b} = \vec{c}\f]
    !! which is written out to
    !!  \f[\left(\begin{array}{cccc}
    !!          x_1^0  & x_1^1  & x_1^2  & \ldots \\
    !!          x_2^0  & x_2^1  & x_2^2  & \ldots \\
    !!          x_3^0  & x_3^1  & x_3^2  & \ldots \\
    !!          \vdots & \vdots & \vdots & \ddots \\
    !!      \end{array}\right)
    !!      \left(\begin{array}{c}
    !!          a_0 \\ a_1 \\ a_2 \\ \vdots
    !!      \end{array}\right)
    !!      =
    !!      \left(\begin{array}{c}
    !!          v_0 \\ v_1 \\ v_2 \\ \vdots
    !!      \end{array}\right)
    !!  \f]
    !! where \f$v_{i-1}\f$ = \c var_i.
    !!
    !! This is solved for the coefficients of the interpolating polynomial
    !!
    !!  <tt>ext_var = a_0 + a_1 x + a_2 x^2 + ...</tt>
    !!
    !! There is an optional flag to look for the \f$k^\text{th}\f$ derivative of
    !! the  function instead  of the  function itself,  where \f$k\f$  should be
    !! lower than the degree of the polynomial.
    !!
    !! \return ierr
    integer function calc_ext_var(ext_var,var,var_points,ext_point,deriv_in) &
        &result(ierr)
        character(*), parameter :: rout_name = 'ext_var'
        
        ! input / output
        real(dp), intent(inout) :: ext_var                                      !< output
        real(dp), intent(in) :: var(:)                                          !< abscissa
        real(dp), intent(in) :: var_points(:)                                   !< ordinate
        real(dp), intent(in) :: ext_point                                       !< point to which extrapolate
        integer, intent(in), optional :: deriv_in                               !< specifies an optional derivative
        
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
    
    !> Calculate the coefficients for finite differences.
    !!
    !! These represent a  derivative of degree \c  deriv at an index  \c ind, by
    !! the weighted sum of a number \c nr of points with <tt>1<=ind<=nr</tt>.
    !!
    !! They are found by solving a Vandermonde system, multiplied by the faculty
    !! of the derivative.
    !!
    !! The result  returned in <tt>coeff(1:nr)</tt> are used in
    !!  \f[\frac{d f}{d x} =
    !!      \sum_{j=1}^n c(j) f(i+j) . \f]
    !!
    !! where \f$i\f$ = \c ind, \f$n\f$ = \c nr and \f$c\f$ = \c coeff
    !!
    !! \note They need to be divided by the step size before usage.
    !!
    !! Examples:
    !!  - symmetric finite differences for derivative \c deriv of order \c ord:
    !!      - \f$m = \text{ceiling}(\frac{p+d}{2})\f$
    !!          to guarantee the order
    !!      - \f$n   = 1+2m\f$
    !!      - \f$i  = 1+m\f$
    !!  - left finite differences for derivative \c deriv of order \c ord:
    !!      - \f$n   = 1+d+p \f$
    !!      - \f$i  = m\f$
    !! where \f$p\f$ = \c ord, \f$d\f$ = \c deriv and \f$m\f$ = \c nr_2.
    !!
    !! \return ierr
    integer function calc_coeff_fin_diff(deriv,nr,ind,coeff) result(ierr)
        character(*), parameter :: rout_name = 'calc_coeff_fin_diff'
        
        ! input / output
        integer, intent(in) :: deriv                                            !< degree of derivative
        integer, intent(in) :: nr                                               !< number of points
        integer, intent(in) :: ind                                              !< position of derivative
        real(dp), intent(inout), allocatable :: coeff(:)                        !< output coefficients
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: mat_loc(:)                                     ! elements of local Vandermonde matrix
        real(dp), allocatable :: rhs_loc(:)                                     ! local right-hand side
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! test  whether  number of  points  nr  and  position ind  positive  and
        ! consistent
        if (deriv.le.0) then
            ierr = 1
            err_msg = 'The degree of the derivative has to be positive'
            CHCKERR(err_msg)
        end if
        if (nr.le.0) then
            ierr = 1
            err_msg = 'The number of points in the finite differences has to &
                &be positive'
            CHCKERR(err_msg)
        end if
        if (ind.lt.1 .or. ind.gt.nr) then
            ierr = 1
            err_msg = 'The position of the derivative is not in range.'
            CHCKERR(err_msg)
        end if
        
        ! test whether number of points enough for derivative
        if (deriv.gt.nr-1) then
            ierr = 1
            err_msg = 'For derivatives of degree '//trim(i2str(deriv))//&
                &' minimally '//trim(i2str(deriv+1))//' points are necessary'
            CHCKERR(err_msg)
        end if
        
        ! set up output
        if (allocated(coeff)) deallocate(coeff)
        allocate(coeff(nr))
        allocate(mat_loc(nr))
        allocate(rhs_loc(nr))
        
        ! set up matrix and rhs
        do id = 1,nr
            mat_loc(id) = (id-ind)*1._dp
        end do
        rhs_loc = 0._dp
        rhs_loc(deriv+1) = 1._dp                                                ! looking for this derivative
        call solve_vand(nr,mat_loc,rhs_loc,coeff,transp=.true.)
        coeff = coeff*fac(deriv)
    end function calc_coeff_fin_diff
    
    !> Convert 2-D coordinates (i,j) to the storage convention used in matrices.
    !!
    !! Their size is by default taken to be 3:
    !!  \f[
    !!      \left(\begin{array}{ccc}
    !!          1 & 4 & 7 \\ 2 & 5 & 8 \\ 3 & 6 & 9
    !!      \end{array}\right)
    !!      \ \text{or} \
    !!      \left(\begin{array}{ccc}
    !!          1 &   &   \\ 2 & 4 &   \\ 3 & 5 & 6
    !!      \end{array}\right) \
    !!      \text{for symmetic matrices} \ .
    !!  \f]
    !! Optionally, this can be changed using \c n.
    !!
    !! The value of \f$c\f$ is then given by
    !!   \f[c = (j-1) n + i\f]
    !! for non-symmetric matrices, and by
    !!  \f[\begin{aligned}
    !!      c &= (j-1)n + i - (j-1)\frac{j}{2} &  \text{if} \ i > j \\
    !!      c &= (i-1)n + j - (i-1)\frac{i}{2} &  \text{if} \ j > i
    !!  \end{aligned}\f]
    !! since \f$\sum_{k=1}^{j-1} k = \frac{(j-1)j}{2}\f$.
    !!
    !! For local indices  in submatrices, the limits in both  dimensions have to
    !! be passed. The results for the  full matrix are then subtracted by amount
    !! corresponding to  the left,  above and  below parts  with respect  to the
    !! submatrix:
    !!  - left:
    !!      \f[\sum_{i=1}^{m(2)-1} n-i+1 =
    !!          (m(2)-1) \left(n+1-\frac{m(2)}{2}\right) \f]
    !!  - above (if positive):
    !!      \f[\sum_{i=1}^j \left(m(1)-m(2)+1-i\right) =
    !!          j \left(m(1)-m(2)+\frac{1}{2} -
    !!          \frac{j^*}{2}\right) \ , \ \text{with} \ j^* =
    !!          \text{min}\left(0,j,m(1)-m(2)+1\right) \f]
    !!  - below:
    !!      \f[\sum_{i=1}^{j-1} n-M(1) =
    !!          \left(n-M(1)\right) \left(j-1\right) \f]
    !!
    !! where \f$m\f$ = \c min and \f$M\f$ = \c max.
    !!
    !! \note
    !!  -# The  submatrix version is  not fast, so  results should be  saved and
    !!  reused.
    !!  -# No checks are done whether the indices make sense.
    integer function c(ij,sym,n,lim_n)
        ! input / output
        integer, intent(in) :: ij(2)                                            !< 2D coords. (i,j)
        logical, intent(in) :: sym                                              !< whether symmetric
        integer, intent(in), optional :: n                                      !< max of 2-D coords.
        integer, intent(in), optional :: lim_n(2,2)                             !< min. and max. of 2-D coords.
        
        ! local variables
        integer :: n_loc                                                        ! local n
        integer :: js                                                           ! j* = min(0,j,min(1)-min(2)+1)
        integer :: ij_loc(2)                                                    ! local ij
        integer :: kd                                                           ! dummy integer
        
        ! set local n
        n_loc = 3
        if (present(n)) n_loc = n
        
        ! set c depending on symmetry
        if (sym) then                                                           ! symmetric
            ! set local ij, refering to the total indices
            ij_loc = ij
            if (present(lim_n)) ij_loc = ij + lim_n(1,:) - 1                    ! displace with minimum indices
            
            ! swap indices if upper diagonal
            if (ij_loc(2).gt.ij_loc(1)) then
                kd = ij_loc(2)
                ij_loc(2) = ij_loc(1)
                ij_loc(1) = kd
            end if
            
            ! get c for whole matrix
            c = nint((ij_loc(2)-1)*n_loc + ij_loc(1) - &
                &(ij_loc(2)-1)*ij_loc(2)*0.5)
            
            if (present(lim_n)) then
                ! subtract if submatrix
                js = max(0,min(ij_loc(2)-lim_n(1,2)+1,lim_n(1,1)-lim_n(1,2)+1))
                c = c - nint((lim_n(1,2)-1)*(n_loc+1-lim_n(1,2)*0.5) &
                    &+js*(lim_n(1,1)-lim_n(1,2)+0.5-js*0.5) &
                    &+(n_loc-lim_n(2,1))*(ij_loc(2)-lim_n(1,2)))
            end if
        else                                                                    ! asymmetric
            if (present(lim_n)) then
                ! set c with local size
                c = (ij(2)-1)*(lim_n(2,1)-lim_n(1,1)+1) + ij(1)
            else
                ! set c with total size
                c = (ij(2)-1)*n_loc + ij(1)
            end if
        end if
    end function c
    
    !> Calculate factorial.
    recursive function fac(n)  result(fact)
        ! input / output
        integer, intent(in) :: n                                                !< input
        integer :: fact                                                         !< output
        
        ! calculate recursively
        if (n.eq.0) then
            fact = 1
        else
            fact = n * fac(n-1)
        end if
    end function fac
    
    !> Determines  whether a  matrix making  use  of the  storage convention  in
    !! eq_vars.eq_2_type is symmetric or not.
    !!
    !! \return ierr
    integer function is_sym(n,nn,sym) result(ierr)
        character(*), parameter :: rout_name = 'is_sym'
        
        ! input / output
        integer, intent(in) :: n                                                !< size of matrix
        integer, intent(in) :: nn                                               !< number of elements in matrix
        logical, intent(inout) :: sym                                           !< output
        
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
    
    !> Returns common multiple using the Euclid's algorithm.
    !!
    !! \see From
    !!  <https://rosettacode.org/wiki/Greatest_common_divisor#Recursive_Euclid_algorithm_3>
    recursive function LCM(u, v) result(res)
        integer, intent(in) :: u                                                !< input
        integer, intent(in) :: v                                                !< input
        integer :: res                                                          !< result
        
        res = u*v / GCD(u,v)
    end function LCM
    
    !> Returns least denominator using the GCD.
    !!
    !! \see From
    !!  <https://rosettacode.org/wiki/Least_common_multiple#Fortran>
    recursive function GCD(u, v) result(res)
        integer, intent(in) :: u                                                !< input
        integer, intent(in) :: v                                                !< input
        integer :: res                                                          !< result
     
        if (mod(u, v) /= 0) then
            res = GCD(v, mod(u, v))
        else
            res = v
        end if
    end function GCD
    
    !> Calculate multiplication  through shifting of fourier modes A  and B into
    !! C.
    !!
    !! \note
    !!  -# These all  are assumed to  have nonzero mode  numers starting from  0. Any
    !!  negative modes are converted into positive ones.
    !!  -# Modes that are larger than what C can hold are thrown away.
    subroutine shift_F(A,B,C)
        ! input / output
        real(dp), intent(in) :: A(0:,:)                                         !< input
        real(dp), intent(in) :: B(0:,:)                                         !< input
        real(dp), intent(inout) :: C(0:,:)                                      !< result
        
        ! local variables
        integer :: i_A, i_B, i_C                                                ! indices in A, B and C
        integer :: n_A, n_B, n_C                                                ! number of modes for A, B and C
        
        ! set number of modes
        n_A = size(A,1)-1
        n_B = size(B,1)-1
        n_C = size(C,1)-1
        
        ! initialize C
        C = 0._dp
        
        ! loop over A
        do i_A = 0,n_A
            do i_B = 0,n_B
                ! contribution to i_A + i_B
                i_C = i_A + i_B
                
                if (i_C.le.n_C) then
                    C(i_C,1) = C(i_C,1) + 0.5* A(i_A,1)*B(i_B,1)
                    C(i_C,2) = C(i_C,2) + 0.5* A(i_A,1)*B(i_B,2)
                    C(i_C,2) = C(i_C,2) + 0.5* A(i_A,2)*B(i_B,1)
                    C(i_C,1) = C(i_C,1) - 0.5* A(i_A,2)*B(i_B,2)
                end if
                
                ! contribution to i_A + i_B
                i_C = i_A - i_B
                
                if (i_C.ge.0) then
                    C(i_C,1) = C(i_C,1) + 0.5* A(i_A,1)*B(i_B,1)
                    C(i_C,2) = C(i_C,2) - 0.5* A(i_A,1)*B(i_B,2)
                    C(i_C,2) = C(i_C,2) + 0.5* A(i_A,2)*B(i_B,1)
                    C(i_C,1) = C(i_C,1) + 0.5* A(i_A,2)*B(i_B,2)
                else
                    C(-i_C,1) = C(-i_C,1) + 0.5* A(i_A,1)*B(i_B,1)
                    C(-i_C,2) = C(-i_C,2) + 0.5* A(i_A,1)*B(i_B,2)
                    C(-i_C,2) = C(-i_C,2) - 0.5* A(i_A,2)*B(i_B,1)
                    C(-i_C,1) = C(-i_C,1) + 0.5* A(i_A,2)*B(i_B,2)
                end if
            end do
        end do
    end subroutine shift_F
    
    !> Solve a Vandermonde system \f$\overline{\text{A}} \ \vec{X} = \vec{B}\f$.
    !!
    !! The Vandermonde matrix has the form
    !!  \f[A = \left(\begin{array}{ccccc}
    !!      1      & a_1    & a_1^2  & \cdots & a_1^n \\
    !!      1      & a_2    & a_2^2  & \cdots & a_2^n \\
    !!      1      & a_3    & a_3^2  & \cdots & a_3^n \\
    !!      \vdots & \vdots & \vdots & \ddots & \vdots \\
    !!      1      & a_n    & a_n^2  & \cdots & a_n^n
    !!  \end{array}\right)\f]
    !!
    !! \see Adapted from two routines \c dvand and \c pvand in
    !! <https://people.sc.fsu.edu/~jburkardt/f_src/vandermonde/vandermonde.html>
    !! based on the Björk-Pereyra Algorithm \cite Bjorck1970.
    subroutine solve_vand(n,a,b,x,transp)
        ! input / output
        integer, intent(in) :: n                                                !< matrix size
        real(dp), intent(in) :: a(n)                                            !< parameters of Vandermonde matrix
        real(dp), intent(in) :: b(n)                                            !< right-hand side
        real(dp), intent(inout) :: x(n)                                         !< solution
        logical, intent(in), optional :: transp                                 !< transposed Vandermonde matrix
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        logical :: transp_loc                                                   ! local transp
        
        ! initialize
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        x(1:n) = b(1:n)
        
        if (transp_loc) then
            do kd = 1, n - 1
                do jd = n, kd + 1, -1
                    x(jd) = x(jd) - a(kd) * x(jd-1)
                end do
            end do
            
            do kd = n - 1, 1, -1
                do jd = kd + 1, n
                    x(jd) = x(jd) / ( a(jd) - a(jd-kd) )
                end do
                do jd = kd, n - 1
                    x(jd) = x(jd) - x(jd+1)
                end do
            end do
        else
            do kd = 1, n - 1
                do jd = n, kd + 1, -1
                    x(jd) = ( x(jd) - x(jd-1) ) / ( a(jd) - a(jd-kd) )
                end do
            end do
            
            do kd = n - 1, 1, -1
                do jd = kd, n - 1
                    x(jd) = x(jd) - a(kd) * x(jd+1)
                end do
            end do
        end if
    end subroutine solve_vand
end module num_utilities
