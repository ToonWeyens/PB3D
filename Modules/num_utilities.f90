!------------------------------------------------------------------------------!
!   Numerical utilities                                                        !
!------------------------------------------------------------------------------!
module num_utilities
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, qp, iu, max_str_ln, pi
    
    implicit none
    private
    public calc_zero_HH, calc_ext_var, calc_det, calc_int, add_arr_mult, c, &
        &conv_FHM, check_deriv, calc_inv, calc_mult, calc_aux_utilities, &
        &derivs, con2dis, dis2con, round_with_tol, conv_mat, is_sym, &
        &con, calc_coeff_fin_diff, fac, d, m, f, calc_zero_Zhang
#if ldebug
    public debug_calc_zero, debug_con2dis_reg, debug_calc_coeff_fin_diff
#endif
    
    ! global variables
    integer, allocatable :: d(:,:,:)                                            ! 1D array indices of derivatives
    integer, allocatable :: m(:,:)                                              ! 1D array indices of metric indices
    integer, allocatable :: f(:,:)                                              ! 1D array indices of Fourier mode combination indices
#if ldebug
    logical :: debug_calc_zero = .false.                                        ! plot debug information for calc_zero
    logical :: debug_con2dis_reg = .false.                                      ! plot debug information for con2dis_reg
    logical :: debug_calc_coeff_fin_diff = .false.                              ! plot debug information for calc_coeff_fin_diff
#endif

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
        module procedure &
            &calc_mult_0D_real, calc_mult_2D_real, calc_mult_2D_complex
    end interface
    interface conv_mat
        module procedure &
            &conv_mat_3D, conv_mat_3D_complex, &
            &conv_mat_0D, conv_mat_0D_complex
    end interface
    interface calc_int
        module procedure calc_int_eqd, calc_int_reg
    end interface
    interface calc_zero_HH
        module procedure calc_zero_HH_0D, calc_zero_HH_3D
    end interface
    interface round_with_tol
        module procedure round_with_tol_arr, round_with_tol_ind
    end interface
    interface con2dis
        module procedure con2dis_eqd, con2dis_reg
    end interface
    interface dis2con
        module procedure dis2con_eqd, dis2con_reg
    end interface
    interface con
        module procedure con_3D, con_2D, con_1D, con_0D
    end interface
    
contains
    ! initialize utilities for fast future reference, depending on program style
    !   - derivatives
    !   - metrics
    !   - Fourier modes (optionally)
    ! If  Fourier modes  are also  initialized, the quantity  "n_mod" has  to be
    ! provided as well.
    subroutine calc_aux_utilities(n_mod)
        use num_vars, only: max_deriv
        
        ! input / output
        integer, intent(in), optional :: n_mod                                  ! n_mod for Fourier modes
        
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
    
    ! Calculate the  1D indices for  derivatives of  a certain order  in certain
    ! dimensions.
    ! The  algorithm  works by  considering a  structure such  as the  following
    ! example in three dimensions:
    !   1: (0,0,0)
    !   2: (1,0,0)
    !   3: (0,1,0)
    !   4: (0,0,1)
    !   5: (2,0,0)
    !   6: (1,1,0)
    !   7: (1,0,1)
    !   8: (0,2,0)
    !   9: (0,1,1)
    !   10: (0,0,2)
    !   11: (3,0,0)
    !   12: (2,1,0)
    ! etc...
    ! By then defining  n_dims = size(deriv), tot_deriv =  sum(deriv) and id_max
    ! as the  index of  the last  nonzero element in  deriv and  extending deriv
    ! towards the left  by considering deriv(0)=0, the following  formula can be
    ! deduced for the displacements in this table with respect to index 1:
    !   sum_jd=0^id_max-1 sum_id=0^(tot_deriv-(deriv(0)+...+deriv(j))-1) 
    !       (n_dims-jd+id-1,id) ,
    ! making use of  the binomial coefficients (a,b) = a!/(b!(a-b)!).  It can be
    ! seen that  each of  the terms in  the summation in  jd corresponds  to the
    ! displacement in dimension  jd and the binomial coefficient  comes from the
    ! classic stars and bars problem.
    integer function calc_derivs_1D_id(deriv,dims) result(res)
        ! input / output
        integer, intent(in) :: deriv(:)                                         ! derivatives
        integer, intent(in) :: dims                                             ! nr. of dimensions
        
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
    
    ! Calculate all of the unique derivatives in certain dimensions of a certain
    ! total order, e.g. for order 2 in 3 dimensions this would be:
    !   (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1) and (0,0,2)
    !  The total  number of  derivatives  is found  through the  stars and  bars
    ! problem to be the binomial coefficient 
    !   (order+dims-1,order) = (order+dims-1)!/(order!(dims-1)!)
    ! The number  of dimensions is taken to  be 3 by default but  can be changed
    ! optionally.
    recursive function derivs(order,dims) result(derivs_res)
        ! input / output
        integer, intent(in) :: order                                            ! order of derivative
        integer, intent(in), optional :: dims                                   ! nr. of dimensions
        integer, allocatable :: derivs_res(:,:)                                 ! array of all unique derivatives
        
        ! local variables
        integer :: id                                                           ! counter
        integer, allocatable :: derivs_loc(:,:)                                 ! derivs of local subblock
        integer :: id_tot                                                       ! counter for total dimension of local derivs.
        integer :: dims_loc                                                     ! local version of dims
        
        ! set up local dims
        dims_loc = 3
        if (present(dims)) dims_loc = dims
        
        ! allocate resulting derivs
        allocate(derivs_res(&
            &dims_loc,fac(order+dims_loc-1)/(fac(order)*fac(dims_loc-1))))
        
        ! iterate over the first index if order greater than 0
        if (order.gt.0) then
            id_tot = 0
            do id = order,0,-1
                ! calculate  the  derivs  of  local  subblock  if  more  than  1
                ! dimension
                if (dims_loc.gt.1) then
                    derivs_loc = derivs(order-id,dims_loc-1)                    ! calculate iteratively subblock
                    derivs_res(1,id_tot+1:id_tot+size(derivs_loc,2)) = id       ! set first dimension
                    derivs_res(2:dims_loc,id_tot+1:id_tot+size(derivs_loc,2)) &
                        &= derivs_loc                                           ! set next dimensions
                    id_tot = id_tot+size(derivs_loc,2)
                else 
                    derivs_res = order
                end if
            end do
        else
            derivs_res = 0
        end if
    end function derivs

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
    ! Note: For periodic  function, the trapezoidal rule works well  only if the
    ! last point of the  grid is included, i.e. the point  where the function is
    ! equal to the first point.
    integer function calc_int_reg(var,x,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_reg'
        
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
    end function calc_int_reg
    
    ! Integrates a function using the trapezoidal rule:
    !   int_1^n f(x) dx = sum_k=1^(n-1) {(f(k+1)+f(k))*delta_x/2},
    ! with n  the number of points,  which are assumed to be  equidistant with a
    ! given step size delta_x.
    ! Note: For periodic  function, the trapezoidal rule works well  only if the
    ! last point of the  grid is included, i.e. the point  where the function is
    ! equal to the first point.
    integer function calc_int_eqd(var,step_size,var_int) result(ierr)
        character(*), parameter :: rout_name = 'calc_int_eqd'
        
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
    end function calc_int_eqd
    
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
    ! that are distributed  between both acording to the  binomial theorem. Both
    ! arrays are given on a 3D or 1D grid.
    integer function add_arr_mult_3_3(arr_1,arr_2,arr_3,deriv) result(ierr)     ! Version with arr_1 and arr_2 in 3 coords.
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
    integer function add_arr_mult_3_1(arr_1,arr_2,arr_3,deriv) result(ierr)     ! Version with arr_1 in 3 coords and arr_2 only in the flux coord.
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
    ! met_type is used.
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
        integer, allocatable :: c_sub(:)                                        ! indices of submatrix in storage convention of met_type
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
                ierr = calc_det_2D(work,A(:,:,:,c_sub),n-1)
                CHCKERR('')
                detA = detA + A(:,:,:,c([1,id],sym,n))*sgn*work
                sgn = -sgn
                slct(id) = 1
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
    
    ! calculate inverse of square matrix A  which has elements depending on a 3D
    ! grid. The storage convention described in met_type is used.
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
        integer, allocatable :: c_sub(:)                                        ! indices of submatrix in storage convention of met_type
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
                ierr = calc_det_2D(detA,A(:,:,:,c_sub),n-1)
                CHCKERR('')
                inv_2D(:,:,:,c([kd,ld],sym,n)) = (-1.0_dp)**(kd+ld)*detA
                slct(ld,1) = 1
                slct(kd,2) = 1
            end do
        end do
        
        ! deallocate c_sub l_sub and k_sub
        deallocate(c_sub,l_sub,k_sub)
        
        ! calculate determinant in detA
        ierr = calc_det(detA,A,n)
        CHCKERR('')
        
        ! divide by determinant
        do kd = 1,nn
            inv_2D(:,:,:,kd) = inv_2D(:,:,:,kd) / detA
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
        
        call dgetrf(n,n,inv_0D,n,ipiv,ierr)                                     ! LU factorization
        err_msg = 'Lapack couldn''t find the LU factorization'
        CHCKERR(err_msg)
        
        call dgetri(n,inv_0D,n,ipiv,w,n,ierr)                                   ! inverse of LU
        CHCKERR('Lapack couldn''t compute the inverse')
    end function calc_inv_0D
    
    ! Calculate matrix multiplication of two square matrices AB = A B which have
    ! elements  defined  on a  3D  grid.  The  storage convention  described  in
    ! met_type is used for the fourth dimension.
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
    integer function calc_mult_0D_real(A,B,AB,n,transp) result(ierr)            ! real version with constant matrix
        character(*), parameter :: rout_name = 'calc_mult_0D_real'
        
        ! input / output
        real(dp), intent(in) :: A(:)                                            ! input A
        real(dp), intent(in) :: B(:)                                            ! input B
        real(dp), intent(inout) :: AB(:)                                        ! output A B
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp(2)                              ! .true. if A and/or B transposed
        
        ! local variables
        real(dp), allocatable :: loc_AB(:,:,:,:)                                ! local copy of C
        
        ! initialize ierr
        ierr = 0
        
        ! allocate local AB
        allocate(loc_AB(1,1,1,size(AB)))
        
        ! calculate using 2D version
        ierr = calc_mult_2D_real(reshape(A,[1,1,1,size(A)]),&
            &reshape(B,[1,1,1,size(B)]),loc_AB,n,transp)
        
        ! update AB with local value
        AB = loc_AB(1,1,1,:)
        
        ! deallocate local variables
        deallocate(loc_AB)
    end function calc_mult_0D_real
    
    ! Converts a  (symmetric) matrix  A defined  on a 3D  grid with  the storage
    ! convention  described  in met_type.  If  the  matrix  is stored  with  n^2
    ! numbers, only  the lower diagonal elements  are kept in matrix  B. If only
    ! the  lower diagonal  elements are  stored, they  are copied  to the  upper
    ! diagonal ones of the matrix B as well.
    ! Optionally, the transpose can be calculated.
    ! Note  that  the  routine  does  not check  whether  the  matrix is  indeed
    ! symmetric and that information may thus be lost after conversion.
    ! Note also that  this routine makes a copy  of A so that by  providing A as
    ! input argument for both A and B overwrites A.
    integer function conv_mat_3D(A,B,n,transp) result(ierr)                     ! version defined on 3D grid
        character(*), parameter :: rout_name = 'conv_mat_3D'
        
        ! input / output
        real(dp), intent(in) :: A(:,:,:,:)                                      ! matrix to be converted
        real(dp), intent(inout) :: B(:,:,:,:)                                   ! converted matrix
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp                                 ! transpose
        
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
    integer function conv_mat_0D(A,B,n,transp) result(ierr)                     ! scalar version
        character(*), parameter :: rout_name = 'conv_mat_0D'
        
        ! input / output
        real(dp), intent(in) :: A(:)                                            ! matrix to be converted
        real(dp), intent(inout) :: B(:)                                         ! converted matrix
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp                                 ! transpose
        
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
    integer function conv_mat_3D_complex(A,B,n,transp) result(ierr)             ! complex version defined on 3D grid
        character(*), parameter :: rout_name = 'conv_mat_3D_complex'
        
        ! input / output
        complex(dp), intent(in) :: A(:,:,:,:)                                   ! matrix to be converted
        complex(dp), intent(inout) :: B(:,:,:,:)                                ! converted matrix
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp                                 ! transpose
        
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
    integer function conv_mat_0D_complex(A,B,n,transp) result(ierr)             ! complex scalar version
        character(*), parameter :: rout_name = 'conv_mat_0D_complex'
        
        ! input / output
        complex(dp), intent(in) :: A(:)                                         ! matrix to be converted
        complex(dp), intent(inout) :: B(:)                                      ! converted matrix
        integer, intent(in) :: n                                                ! size of matrix
        logical, intent(in), optional :: transp                                 ! transpose
        
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
    
    ! Finds the  zero of  a function using  Householder iteration.  If something
    ! goes wrong,  by default multiple tries  can be attempted, by  lowering the
    ! relaxation  factor. If  still nothing  is  achieved, an  error message  is
    ! returned, that is empty otherwise.
    function calc_zero_HH_0D(zero,fun,ord,guess,relax_fac,max_nr_tries) &
        &result(err_msg)
        use num_vars, only: max_it_zero, tol_zero, relax_fac_HH, max_nr_tries_HH
        
        ! input / output
        real(dp), intent(inout) :: zero                                         ! output
        interface
            function fun(x,ord)                                                 ! the function
                use num_vars, only: dp
                real(dp), intent(in) :: x 
                integer, intent(in) :: ord
                real(dp) :: fun 
            end function fun
        end interface
        integer, intent(in) :: ord                                              ! order of solution
        real(dp), intent(in) :: guess                                           ! first guess
        real(dp), intent(in), optional :: relax_fac                             ! relaxation factor
        integer, intent(in), optional :: max_nr_tries                           ! max nr. of tries with different relaxation factors
        character(len=max_str_ln) :: err_msg                                    ! possible error message
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: max_nr_tries_loc                                             ! local max_nr_tries
        real(dp) :: corr                                                        ! correction
        real(dp), allocatable :: fun_vals(:)                                    ! function values
        real(dp) :: relax_fac_loc                                               ! local relaxation factor
#if ldebug
        real(dp), allocatable :: corrs(:)                                       ! corrections for all steps
        real(dp), allocatable :: values(:)                                      ! values for all steps
#endif
        
        ! initialize error message
        err_msg = ''
        
        ! set up zero
        zero = guess
        
        ! set up local max_nr_tries
        max_nr_tries_loc = max_nr_tries_HH
        if (present(max_nr_tries)) max_nr_tries_loc = max_nr_tries
        
        ! set up local relaxation factor
        relax_fac_loc = relax_fac_HH
        if (present(relax_fac)) relax_fac_loc = relax_fac
        
        ! initialize function values
        allocate(fun_vals(0:ord))
        
        ! possibly multiple tries with different relaxation factors
        do id = 1,max_nr_tries_loc
#if ldebug
            if (debug_calc_zero) then
                ! set up corrs
                allocate(corrs(max_it_zero))
                allocate(values(max_it_zero))
            end if
#endif
            
            HH: do jd = 1,max_it_zero
                ! calculate function values
                do kd = 0,ord
                    fun_vals(kd) = fun(zero,kd)
                end do
                
                ! correction to theta_HH
                select case (ord)
                    case (1)                                                    ! Newton-Rhapson
                        corr = -fun_vals(0)/fun_vals(1)
                    case (2)                                                    ! Halley
                        corr = -fun_vals(0)*fun_vals(1)/&
                            &(fun_vals(1)**2 - 0.5_dp*&
                            &fun_vals(0)*fun_vals(2))
                    case (3)                                                    ! 3rd order
                        corr = -(6*fun_vals(0)*fun_vals(1)**2-&
                            &3*fun_vals(0)**2*fun_vals(2))/&
                            &(6*fun_vals(1)**3 - 6*fun_vals(0)*&
                            &fun_vals(1)*fun_vals(2)+&
                            &fun_vals(0)**2*fun_vals(3))
                end select
#if ldebug
                if (debug_calc_zero) then
                    corrs(jd) = corr
                    values(jd) = zero
                end if
#endif
                zero = zero + relax_fac_loc/id*corr                             ! relaxation factor scaled by id
                
                ! check for convergence
                if (abs(corr).lt.tol_zero) then
#if ldebug
                    if (debug_calc_zero) &
                        &call plot_evolution(corrs(1:jd),values(1:jd))
#endif
                    return
                else if (jd .eq. max_it_zero) then
#if ldebug
                    if (debug_calc_zero) then
                        call plot_evolution(corrs,values)
                        deallocate(corrs)
                        deallocate(values)
                        if (id.lt.max_nr_tries_loc) &
                            &call writo('Trying again with relaxation factor &
                            &divided by '//trim(i2str(id)))
                    end if
#endif
                    if (id.eq.max_nr_tries_loc) then
                        err_msg = 'Not converged after '//trim(i2str(jd))//&
                            &' iterations, with residual '//&
                            &trim(r2strt(corr))//' and final value '//&
                            &trim(r2strt(zero))
                        zero = 0.0_dp
                    end if
                end if
            end do HH
        end do
#if ldebug
    contains
        ! plots corrections
        subroutine plot_evolution(corrs,values)
            ! input / output
            real(dp), intent(in) :: corrs(:)                                    ! corrections
            real(dp), intent(in) :: values(:)                                   ! values
            
            ! local variables
            real(dp) :: min_x, max_x                                            ! min. and max. x for which to plot the function
            real(dp) :: extra_x = 1._dp                                         ! how much bigger to take the plot interval than just min and max
            real(dp) :: delta_x                                                 ! interval between zero and guess
            integer :: n_x = 200                                                ! how many x values to plot
            real(dp), allocatable :: x_out(:,:)                                 ! output x of plot
            real(dp), allocatable :: y_out(:,:)                                 ! output y of plot
            
            ! set up min and max x
            min_x = min(zero,guess)
            max_x = max(zero,guess)
            delta_x = max_x-min_x
            min_x = min_x - delta_x*extra_x
            max_x = max_x + delta_x*extra_x
            
            ! set up output of plot
            allocate(x_out(n_x,2))
            allocate(y_out(n_x,2))
            x_out(:,1) = [(min_x+(max_x-min_x)*(jd-1._dp)/(n_x-1),jd=1,n_x)]
            x_out(:,2) = [(min_x+(max_x-min_x)*(jd-1._dp)/(n_x-1),jd=1,n_x)]
            y_out(:,1) = [(fun(min_x+(max_x-min_x)*(jd-1._dp)/(n_x-1),0),&
                &jd=1,n_x)]
            y_out(:,2) = [(fun(min_x+(max_x-min_x)*(jd-1._dp)/(n_x-1),1),&
                &jd=1,n_x)]
            
            ! plot function
            call writo('Initial guess: '//trim(r2str(guess)))
            call writo('Resulting zero: '//trim(r2str(zero)))
            call print_GP_2D('function and derivative','',y_out,x=x_out)
            
            ! user output
            call writo('Last correction '//trim(r2strt(corrs(size(corrs))))&
                &//' with tolerance '//trim(r2strt(tol_zero)))
            
            ! plot corrections and values
            call print_GP_2D('corrs','',corrs)
            call print_GP_2D('values','',values)
        end subroutine plot_evolution
#endif
    end function calc_zero_HH_0D
    function calc_zero_HH_3D(dims,zero,fun,ord,guess,relax_fac,&
        &max_nr_tries) result(err_msg)
        use num_vars, only: max_it_zero, tol_zero, relax_fac_HH, max_nr_tries_HH
        
        ! input / output
        integer, intent(in) :: dims(3)                                          ! dimensions of the problem
        real(dp), intent(inout) :: zero(dims(1),dims(2),dims(3))                ! output
        interface
            function fun(dims,x,ord)                                            ! the function
                use num_vars, only: dp
                integer, intent(in) :: dims(3)
                real(dp), intent(in) :: x(dims(1),dims(2),dims(3)) 
                integer, intent(in) :: ord
                real(dp) :: fun(dims(1),dims(2),dims(3)) 
            end function fun
        end interface
        integer, intent(in) :: ord                                              ! order of solution
        real(dp), intent(in) :: guess(dims(1),dims(2),dims(3))                  ! first guess
        real(dp), intent(in), optional :: relax_fac                             ! relaxation factor
        integer, intent(in), optional :: max_nr_tries                           ! max nr. of tries with different relaxation factors
        character(len=max_str_ln) :: err_msg                                    ! possible error message
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: max_nr_tries_loc                                             ! local max_nr_tries
        integer :: mc_ind(3)                                                    ! index of maximum correction
        real(dp) :: corr(dims(1),dims(2),dims(3))                               ! correction
        real(dp) :: relax_fac_loc                                               ! local relaxation factor
        real(dp), allocatable :: fun_vals(:,:,:,:)                              ! function values
#if ldebug
        real(dp), allocatable :: corrs(:,:,:,:)                                 ! corrections for all steps
        real(dp), allocatable :: values(:,:,:,:)                                ! values for all steps
#endif
        
        ! initialize error message
        err_msg = ''
        
        ! test
        if (ord.lt.1 .or. ord.gt.3) then
            zero = 0.0_dp
            err_msg = 'only orders 1 (Newton-Rhapson), 2 (Halley) and 3 &
                &implemented'
            return
        end if
        
        ! set up zero
        zero = guess
        
        ! set up local max_nr_tries
        max_nr_tries_loc = max_nr_tries_HH
        if (present(max_nr_tries)) max_nr_tries_loc = max_nr_tries
        
        ! set up local relaxation factor
        relax_fac_loc = relax_fac_HH
        if (present(relax_fac)) relax_fac_loc = relax_fac
        
        ! initialize function values
        allocate(fun_vals(dims(1),dims(2),dims(3),0:ord))
        
        ! possibly multiple tries with different relaxation factors
        do id = 1,max_nr_tries_loc
#if ldebug
            if (debug_calc_zero) then
                ! set up corrs
                allocate(corrs(dims(1),dims(2),dims(3),max_it_zero))
                allocate(values(dims(1),dims(2),dims(3),max_it_zero))
            end if
#endif
            
            HH: do jd = 1,max_it_zero
                ! calculate function values
                do kd = 0,ord
                    fun_vals(:,:,:,kd) = fun(dims,zero,kd)
                end do
                
                ! correction to theta_HH
                select case (ord)
                    case (1)                                                    ! Newton-Rhapson
                        corr = -fun_vals(:,:,:,0)/fun_vals(:,:,:,1)
                    case (2)                                                    ! Halley
                        corr = -fun_vals(:,:,:,0)*fun_vals(:,:,:,1)/&
                            &(fun_vals(:,:,:,1)**2 - 0.5_dp*&
                            &fun_vals(:,:,:,0)*fun_vals(:,:,:,2))
                    case (3)                                                    ! 3rd order
                        corr = -(6*fun_vals(:,:,:,0)*fun_vals(:,:,:,1)**2-&
                            &3*fun_vals(:,:,:,0)**2*fun_vals(:,:,:,2))/&
                            &(6*fun_vals(:,:,:,1)**3 - 6*fun_vals(:,:,:,0)*&
                            &fun_vals(:,:,:,1)*fun_vals(:,:,:,2)+&
                            &fun_vals(:,:,:,0)**2*fun_vals(:,:,:,3))
                end select
#if ldebug
                if (debug_calc_zero) then
                    corrs(:,:,:,jd) = corr
                    values(:,:,:,jd) = zero
                end if
#endif
                zero = zero + relax_fac_loc/id*corr                             ! relaxation factor scaled by id
                
                ! check for convergence
                if (maxval(abs(corr)).lt.tol_zero) then
#if ldebug
                    if (debug_calc_zero) call &
                        &plot_evolution(corrs(:,:,:,1:jd),values(:,:,:,1:jd))
#endif
                    return
                else if (jd .eq. max_it_zero) then
#if ldebug
                    if (debug_calc_zero) then
                        call plot_evolution(corrs,values)
                        deallocate(corrs)
                        deallocate(values)
                        if (id.lt.max_nr_tries_loc) &
                            &call writo('Trying again with relaxation factor &
                            &divided by '//trim(i2str(id)))
                    end if
#endif
                    if (id.eq.max_nr_tries_loc) then
                        mc_ind = maxloc(abs(corr))
                        err_msg = 'Not converged after '//trim(i2str(jd))//&
                            &' iterations, with maximum residual '//&
                            &trim(r2strt(corr(mc_ind(1),mc_ind(2),mc_ind(3))))&
                            &//' and final value '//trim(r2strt(&
                            &zero(mc_ind(1),mc_ind(2),mc_ind(3))))
                        zero = 0.0_dp
                    end if
                end if
            end do HH
        end do
#if ldebug
    contains
        ! plots corrections
        subroutine plot_evolution(corrs,values)
            ! input / output
            real(dp), intent(in) :: corrs(:,:,:,:)                              ! corrections
            real(dp), intent(in) :: values(:,:,:,:)                             ! values
            
            ! local variables
            integer :: n_corrs                                                  ! number of corrections
            character(len=max_str_ln) :: var_name(1)                            ! names of variable
            character(len=max_str_ln) :: file_name                              ! name of file
            
            ! initialize
            n_corrs = size(corrs,4)
            mc_ind = maxloc(abs(corrs(:,:,:,n_corrs)))
            
            ! user output
            call writo('Last maximum correction '//&
                &trim(r2strt(corrs(mc_ind(1),mc_ind(2),mc_ind(3),n_corrs)))&
                &//' and final value '//trim(r2strt(&
                &values(mc_ind(1),mc_ind(2),mc_ind(3),n_corrs)))//&
                &' with tolerance '//trim(r2strt(tol_zero)))
            
            ! plot corrs
            file_name = 'corrs'
            var_name = 'var'
            
            call plot_HDF5(var_name,file_name,corrs,col=2,&
                &description='corrections')
            
            ! plot values
            file_name = 'values'
            var_name = 'var'
            
            call plot_HDF5(var_name,file_name,values,col=2,&
                &description='values')
        end subroutine plot_evolution
#endif
    end function calc_zero_HH_3D
    
    ! Finds the zero of a function  using Zhang's method, which is simpler than 
    ! Brent's  method. Taken  from from  Steven Stage's  correction of  Zhang's 
    ! paper [http://www.cscjournals.org/manuscript/Journals/IJEA/Volume4/Issue1
    ! /IJEA-33.pdf].
    ! Unlike Householder,  Zhang's method needs  an interval "x_int_in"  to work
    ! in, not a guess. Note that this  interval needs to be so that the function
    ! values at either end are of different value.
    function calc_zero_Zhang(zero,fun,x_int_in) result(err_msg)
        use num_vars, only: max_it_zero, tol_zero
        
        ! input / output
        real(dp), intent(inout) :: zero                                         ! output
        interface
            function fun(x)
                use num_vars, only: dp
                real(dp), intent(in) :: x 
                real(dp) :: fun 
            end function fun
        end interface
        real(dp), intent(in) :: x_int_in(2)                                     ! interval
        character(len=max_str_ln) :: err_msg                                    ! possible error message
        
        ! local variables
        logical :: converged                                                    ! whether converged
        integer :: jd                                                           ! counter
        integer :: id_loc                                                       ! local id (either 1 or 2)
        real(dp) :: x_int(2)                                                    ! points of current interval [a,b]
        real(dp) :: x_mid                                                       ! x at midpoint [c]
        real(dp) :: x_rule                                                      ! x found with rule [s]
        real(dp) :: x_temp                                                      ! temporary x
        real(dp) :: fun_int(2)                                                  ! function at x_int
        real(dp) :: fun_mid                                                     ! function at x_mid
        real(dp) :: fun_rule                                                    ! function at x_rule
#if ldebug
        real(dp), allocatable :: x_ints(:,:)                                    ! bounds for all steps
#endif
        
        ! initialize error message and zero
        err_msg = ''
        zero = 0.0_dp
        
        ! initialize from input
        x_int(1) = minval(x_int_in)
        x_int(2) = maxval(x_int_in)
        fun_int = [fun(x_int(1)),fun(x_int(2))]
        
        ! test whether already zero
        if (abs(fun_int(1)).lt.tol_zero) then
            zero = x_int(1)
            return
        end if
        if (abs(fun_int(2)).lt.tol_zero) then
            zero = x_int(2)
            return
        end if
        
        ! test whether different sign
        if (product(fun_int).gt.0) then
            err_msg = 'Function values on starting interval have same sign'
            return
        end if
        
        ! intialize converged
        converged = .false.
        
#if ldebug
        if (debug_calc_zero) then
            ! set up x_ints
            allocate(x_ints(max_it_zero,2))
        end if
#endif
        
        ! iterations
        do jd = 1,max_it_zero
            ! calculate midpoint and function value
            x_mid = 0.5_dp*sum(x_int)
            fun_mid = fun(x_mid)
            
            if (abs(fun_int(1)-fun_mid).gt.tol_zero .and. &
                &abs(fun_int(2)-fun_mid).gt.tol_zero) then
                ! three  different  function   values:   use  inverse  quadratic
                ! interpolation
                x_rule = fun_mid*fun_int(2)*x_int(1)/&
                    &((fun_int(1)-fun_mid)*(fun_int(1)-fun_int(2))) + &
                    &fun_int(1)*fun_int(2)*x_mid/&
                    &((fun_mid-fun_int(1))*(fun_mid-fun_int(2))) + &
                    &fun_int(1)*fun_mid*x_int(2)/&
                    &((fun_int(2)-fun_mid)*(fun_int(2)-fun_mid))
                
                if (x_int(1).lt.x_rule .and. x_rule.lt.x_int(2)) then
                    ! found x within interval
                else
                    ! use bisection instead
                    if (fun_int(1)*fun_mid.lt.0) then                           ! [a,c] contains sign change
                        id_loc = 1
                    else if (fun_int(2)*fun_mid.lt.0) then                      ! [b,c] contains sign change
                        id_loc = 2
                    else
                        err_msg = 'Something went wrong...'
                        converged = .true.
                    end if
                    x_rule = 0.5_dp*(x_int(id_loc)+x_mid)
                end if
            else
                ! two function values overlap so use secant rule
                if (fun_int(1)*fun_mid.lt.0) then                               ! [a,c] contains sign change
                    id_loc = 1
                else if (fun_int(2)*fun_mid.lt.0) then                          ! [b,c] contains sign change
                    id_loc = 2
                else
                    err_msg = 'Something went wrong...'
                    converged = .true.
                end if
                x_rule = x_int(id_loc) - &
                    &fun_int(id_loc)*(x_int(id_loc)-x_mid)/&
                    &(fun_int(id_loc)-fun_mid)
            end if
            
            ! calculate fun(s)
            fun_rule = fun(x_rule)
            
            ! make sure x_mid < x_rule
            if (x_mid.gt.x_rule) then
                x_temp = x_rule
                x_rule = x_mid
                x_mid = x_temp
            end if
            
            ! check whether [c,s] contains root
            if (fun_mid*fun_rule.lt.0) then
                x_int(1) = x_mid
                x_int(2) = x_rule
            else
                if (fun_int(1)*fun_mid.lt.0) then
                    x_int(2) = x_mid
                else
                    x_int(1) = x_rule
                end if
            end if
            
#if ldebug
            if (debug_calc_zero) then
                x_ints(jd,:) = x_int
            end if
#endif
            
            ! check for convergence
            fun_int = [fun(x_int(1)),fun(x_int(2))]
            if (abs(fun_int(1)).lt.tol_zero) then
                zero = x_int(1)
                converged = .true.
            end if
            if (abs(fun_int(2)).lt.tol_zero) then
                zero = x_int(2)
                converged = .true.
            end if
            if ((x_int(2)-x_int(1)).lt.tol_zero) then
                zero = 0.5_dp*sum(x_int)
                converged = .true.
            end if
            
            if (converged) then
#if ldebug
                if (debug_calc_zero) call plot_evolution(x_ints(1:jd,:))
#endif
                return
            end if
        end do
#if ldebug
    contains
        ! plots corrections
        subroutine plot_evolution(x_ints)
            ! input / output
            real(dp), intent(in) :: x_ints(:,:)                                 ! x_intervals
            
            ! local variables
            real(dp) :: min_x, max_x                                            ! min. and max. x for which to plot the function
            integer :: n_x = 200                                                ! how many x values to plot
            integer :: n_ints                                                   ! nr. of intervals
            real(dp), allocatable :: x_out(:)                                   ! output x of plot
            real(dp), allocatable :: y_out(:)                                   ! output y of plot
            
            ! set up nr. of intervals
            n_ints = size(x_ints,1)
            
            ! set up min and max x
            min_x = x_int_in(1)
            max_x = x_int_in(2)
            
            ! set up output of plot
            allocate(x_out(n_x))
            allocate(y_out(n_x))
            x_out = [(min_x+(max_x-min_x)*(jd-1._dp)/(n_x-1),jd=1,n_x)]
            y_out = [(fun(min_x+(max_x-min_x)*(jd-1._dp)/(n_x-1)),jd=1,n_x)]
            
            ! plot function
            call writo('Initial interval: ['//trim(r2str(x_int_in(1)))//&
                &'..'//trim(r2str(x_int_in(2)))//']')
            call writo('Resulting zero: '//trim(r2str(zero)))
            call print_GP_2D('function','',y_out,x=x_out)
            
            ! user output
            call writo('Final interval ['//trim(r2str(x_ints(n_ints,1)))//&
                &'..'//trim(r2str(x_ints(n_ints,2)))//'] with tolerance '//&
                &trim(r2strt(tol_zero)))
            
            ! plot intervals
            call print_GP_2D('x_intervals','',x_ints)
        end subroutine plot_evolution
#endif
    end function calc_zero_Zhang

    ! Convert  between  points  from  a  continuous grid  to  a  discrete  grid,
    ! providing either  the the limits on  the grid (lim_c and  lim_d), in which
    ! case the grid is assumed to be equidistant, or the grid values themselves,
    ! in which case the grid just has to be regular.
    ! The output is a  real value where the floored integer is  the index in the
    ! discrete grid  and the remainder  corresponds to the fraction  towards the
    ! next index.  If no solution  is found, a  negative value is  outputted, as
    ! well as a message.
    integer function con2dis_eqd(pt_c,pt_d,lim_c,lim_d) result(ierr)            ! equidistant version
        !character(*), parameter :: rout_name = 'con2dis_eqd'
        
        ! input / output
        real(dp), intent(in) :: pt_c                                            ! point on continous grid
        real(dp), intent(inout) :: pt_d                                         ! point on discrete grid
        real(dp), intent(in) :: lim_c(2)                                        ! [min_c,max_c]
        integer, intent(in) :: lim_d(2)                                         ! [min_d,max_d]
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate, using the formula:
        !    pt_c - min_c     pt_d - min_d
        !   ------------- =  ------------- ,
        !   max_c - min_c    max_d - min_d
        pt_d = lim_d(1) + (lim_d(2)-lim_d(1)) * (pt_c-lim_c(1)) / &
            &(lim_c(2)-lim_c(1))
    end function con2dis_eqd
    integer function con2dis_reg(pt_c,pt_d,var_c) result(ierr)                  ! regular grid version
        character(*), parameter :: rout_name = 'con2dis_reg'
        
        ! input / output
        real(dp), intent(in) :: pt_c                                            ! point on continous grid
        real(dp), intent(inout) :: pt_d                                         ! point on discrete grid
        real(dp), intent(in) :: var_c(:)                                        ! continous grid values
        
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
            if (var_c_loc(id).le.pt_c_loc*(1+tol)) then
                ind_lo = id
#if ldebug 
                if (debug_con2dis_reg) call writo('for iteration '//&
                    &trim(i2str(id))//'/'//trim(i2str(size_c))//', ind_lo = '//&
                    &trim(i2str(ind_lo)))
#endif
            end if
            if (var_c_inv(id).ge.pt_c_loc*(1-tol)) then
                ind_hi = size_c+1-id
#if ldebug 
                if (debug_con2dis_reg) call writo('for iteration '//&
                    &trim(i2str(id))//'/'//trim(i2str(size_c))//', ind_hi = '//&
                    &trim(i2str(ind_hi)))
#endif
            end if
        end do
#if ldebug
        if (debug_con2dis_reg) then
            call writo('final ind_lo = '//trim(i2str(ind_lo))//', ind_hi = '//&
                &trim(i2str(ind_hi)))
            call print_GP_2D('var_c_loc, var_c_inv','',&
                &reshape([var_c_loc,var_c_inv],[size_c,2]))
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
            CHCKERR('')
        end if
        
        ! deallocate
        deallocate(var_c_loc,var_c_inv)
    end function con2dis_reg
    
    ! Convert  between  points  from  a  discrete grid  to  a  continuous  grid,
    ! providing either  the the limits on  the grid (lim_c and  lim_d), in which
    ! case the grid is assumed to be equidistant, or the grid values themselves,
    ! in which case the grid just has to be regular.
    ! The output is a real value. If  the discrete value lies outside the range,
    ! in the case of a regular grid, a negative value is outputted, as well as a
    ! message.
    integer function dis2con_eqd(pt_d,pt_c,lim_d,lim_c) result(ierr)            ! equidistant version
        !character(*), parameter :: rout_name = 'dis2con_eqd'
        
        ! input / output
        integer, intent(in) :: pt_d                                             ! point on discrete grid
        real(dp), intent(inout) :: pt_c                                         ! point on continous grid
        integer, intent(in) :: lim_d(2)                                         ! [min_d,max_d]
        real(dp), intent(in) :: lim_c(2)                                        ! [min_c,max_c]
        
        ! initialize ierr
        ierr = 0
        
        ! Calculate, using the formula:
        !    pt_c - min_c     pt_d - min_d
        !   ------------- =  ------------- ,
        !   max_c - min_c    max_d - min_d
        pt_c = lim_c(1) + (lim_c(2)-lim_c(1)) * (pt_d-lim_d(1)) / &
            &(lim_d(2)-lim_d(1))
    end function dis2con_eqd
    integer function dis2con_reg(pt_d,pt_c,var_c) result(ierr)                  ! regular grid version 
        character(*), parameter :: rout_name = 'dis2con_reg'
        
        ! input / output
        integer, intent(in) :: pt_d                                             ! point on discrete grid
        real(dp), intent(inout) :: pt_c                                         ! point on continous grid
        real(dp), intent(in) :: var_c(:)                                        ! continous grid values
        
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
    
    ! Calculate the  coefficients for central finite  differences representing a
    ! derivative of degree deriv and order  ord, referring to (half) the stencil
    ! width.
    ! As stated in [ADD REF], these are found by solving the matrix equation
    !   sum_(j=1)^(2ord+1) (j+1-ord)^(i-1)/(i-1)! coeff(j) = 2 for i = deriv ,
    !                                                        0 else .
    ! with i = 1..2ord+1.
    ! The result  returned in coeff(j-1-ord), j=1..2ord+1 are used in
    !   df/dx = sum_(j=-ord)^ord coeff(j-1-ord) f(j) .
    integer function calc_coeff_fin_diff(deriv,ord,coeff) result(ierr)
        character(*), parameter :: rout_name = 'calc_coeff_fin_diff'
        
        ! input / output
        integer, intent(in) :: deriv                                            ! degree of derivative
        integer, intent(in) :: ord                                              ! order
        real(dp), intent(inout), allocatable :: coeff(:)                        ! output coefficients
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: A(:,:)                                         ! matrix of problem
        real(dp), allocatable :: b(:)                                           ! rhs of problem
        integer, allocatable :: ipiv(:)                                         ! pivot indices
        integer :: id, jd                                                       ! counter
        integer :: fac                                                          ! (i-1)!
        
        ! initialize ierr
        ierr = 0
        
        ! test whether order and degree of derivative positive
        if (deriv.le.0) then
            ierr = 1
            err_msg = 'The degree of the derivative has to be positive'
            CHCKERR(err_msg)
        end if
        if (ord.le.0) then
            ierr = 1
            err_msg = 'The order of the finite differences has to be positive'
            CHCKERR(err_msg)
        end if
        ! test whether order high enough for degree of derivative
        if (deriv.gt.2*ord+1) then
            ierr = 1
            err_msg = 'For derivatives of degree '//trim(i2str(deriv))//&
                &' a minimal order of '//trim(i2str(ceiling((deriv-1.)/2)))&
                &//' is necsesary'
            CHCKERR(err_msg)
        end if
        
        ! set up output
        if (allocated(coeff)) deallocate(coeff)
        allocate(coeff(1:2*ord+1))
        
        ! set up matrix
        allocate(A(1:2*ord+1,1:2*ord+1))
        A = 0._dp
        do jd = 1,2*ord+1
            fac = 1
            do id = 1,size(A,1)
                A(id,jd) = (jd-1._dp-ord)**(id-1._dp)/fac
                fac = fac*id
            end do
        end do
        
        ! set up rhs
        allocate(b(1:2*ord+1))
        b = 0._dp
        b(1+deriv) = 1._dp
        
#if ldebug
        if (debug_calc_coeff_fin_diff) then
            write(*,*) 'solving system of linear equations A x = b'
            write(*,*) 'with A = '
            call print_ar_2(A)
            write(*,*) 'and b = '
            call print_ar_1(b)
        end if
#endif
        
        ! solve system
        allocate(ipiv(1:2*ord+1))
        call dgesv(size(b),1,A,size(b),ipiv,b,size(b),ierr)
        CHCKERR('Could not find solution')
        
        ! set solution
        coeff = b
        
#if ldebug
        if (debug_calc_coeff_fin_diff) then
            write(*,*) 'solution x = '
            call print_ar_1(coeff)
        end if
#endif
        
        ! test whether numerical coefficients matrix has right size
        if (size(coeff).ne.2*ord+1) then
            ierr = 1
            err_msg = 'Array of discretization coefficients needs to have &
                &the correct size 2*ord+1'
            CHCKERR(err_msg)
        end if
    end function calc_coeff_fin_diff
    
    ! Convert 2D coordinates (i,j) to the storage convention used in matrices.
    ! Their size is by default taken to be 3:
    !   (1 4 7)      (1    )
    !   (2 5 8)  or  (2 4  ) for symmetric matrices.
    !   (3 6 9)      (3 5 6)
    ! Optionally, this can be changed using 'n'.
    ! The value of c is then given by
    !   c = (j-1)*n + i
    ! for non-symmetric matrices, and by
    !   c = (j-1)*n + i - (j-1)*j/2   if i.ge.j
    !   c = (i-1)*n + j - (i-1)*i/2   if j.ge.i
    ! since sum_k=1^j-1 = (j-1)*j/2
    ! For submatrices, the limits in both dimensions have to be passed.
    ! The  results   for  the  full   matrix  are  then  subtracted   by  amount
    ! corresponding  to the  left, above  and below  parts with  respect to  the
    ! submatrix:
    !   - left: sum_i=1^(min(2)-1) (n-i+1)
    !       = (min(2)-1) (n+1-min(2)/2)
    !   - above: sum_i=1^j (min(1)-min(2)+1-i) if positive
    !       = j* (min(1)-min(2)+1/2 - j*/2), with j* = min(0,j,min(1)-min(2)+1)
    !   - below: sum_i=1^(j-1) (n-max(1))
    !       = (n-max(1)) (j-1)
    ! Note: the  submatrix version is not  fast, so results should  be saved and
    ! reused.
    ! Note: No checks are done whether the indices make sense.
    integer function c(ij,sym,n,lim_n)
        ! input / output
        integer, intent(in) :: ij(2)                                            ! 2D coords. (i,j)
        logical, intent(in) :: sym                                              ! .true. if symmetric
        integer, intent(in), optional :: n                                      ! max of 2D coords.
        integer, intent(in), optional :: lim_n(2,2)                             ! min. and max. of 2D coords.
        
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
    
    ! calculate factorial
    recursive function fac(n)  result(fact)
        ! input / output
        integer, intent(in) :: n
        integer :: fact
        
        ! calculate recursively
        if (n.eq.0) then
            fact = 1
        else
            fact = n * fac(n-1)
        end if
    end function fac
end module num_utilities