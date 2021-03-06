!------------------------------------------------------------------------------!
!> Numerical operations.
!------------------------------------------------------------------------------!
module num_ops
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    use num_utilities
    use output_ops
    
    implicit none
    private
    public calc_zero_HH, calc_zero_Zhang
#if ldebug
    public debug_calc_zero
#endif
    
    ! global variables
#if ldebug
    !> \ldebug
    logical :: debug_calc_zero = .false.                                        !< plot debug information for calc_zero
#endif

    ! interfaces
    
    !> \public Finds the zero of a function using Householder iteration.
    !!
    !! If something goes  wrong, by default multiple tries can  be attempted, by
    !! backtracking on the correction by  multiplying it by a relaxation factor.
    !! This can be done \c max_nr_backtracks times.
    !!
    !! If still nothing is achieved, an error message is returned,
    !! that is empty otherwise.
    interface calc_zero_HH
        !> \public
        module procedure calc_zero_HH_0D
        !> \public
        module procedure calc_zero_HH_3D
    end interface
    
contains
    !> \private 0-D version
    !! \param[inout] fun <tt>fun(x,ord)</tt> with
    !!  - <tt>x</tt> abscissa
    !!  - <tt>ord</tt> order of derivative
    !!  - <tt>fun</tt> ordinate
    function calc_zero_HH_0D(zero,fun,ord,guess,max_nr_backtracks,output) &
        &result(err_msg)
        use num_vars, only: max_it_zero, tol_zero, max_nr_backtracks_HH
        
        ! input / output
        real(dp), intent(inout) :: zero                                         !< output
        interface
            function fun(x,ord)                                                 ! the function
                use num_vars, only: dp
                real(dp), intent(in) :: x
                integer, intent(in) :: ord
                real(dp) :: fun
            end function fun
        end interface
        integer, intent(in) :: ord                                              !< order of solution
        real(dp), intent(in) :: guess                                           !< first guess
        integer, intent(in), optional :: max_nr_backtracks                      !< max nr. of tries with different relaxation factors
        logical, intent(in), optional :: output                                 !< give output on convergence
        character(len=max_str_ln) :: err_msg                                    !< possible error message
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: max_nr_backtracks_loc                                        ! local max_nr_backtracks
        logical :: output_loc                                                   ! local output
        logical :: relaxed_enough                                               ! whether step was relaxed enough to get a smaller new value
        real(dp) :: zero_new                                                    ! new zero calculated
        real(dp) :: fun_new                                                     ! function value at new zero
        real(dp) :: corr                                                        ! correction
        real(dp), allocatable :: fun_vals(:)                                    ! function values
#if ldebug
        real(dp), allocatable :: corrs(:)                                       ! corrections for all steps
        real(dp), allocatable :: values(:)                                      ! values for all steps
#endif
        
        ! initialize error message
        err_msg = ''
        
        ! set up local output
        output_loc = .false.
        if (present(output)) output_loc = output
        
        ! set up zero
        zero = guess
        
        ! set up local max_nr_backtracks
        max_nr_backtracks_loc = max_nr_backtracks_HH
        if (present(max_nr_backtracks)) &
            &max_nr_backtracks_loc = max_nr_backtracks
        
        ! initialize function values
        allocate(fun_vals(0:ord))
        
#if ldebug
        if (debug_calc_zero) then
            ! set up corrs
            allocate(corrs(max_it_zero))
            allocate(values(max_it_zero))
        end if
#endif
        
        ! loop over iterations
        HH: do jd = 1,max_it_zero
            ! calculate function values
            do kd = 0,ord
                fun_vals(kd) = fun(zero,kd)
            end do
            
            ! correction to theta_HH
            select case (ord)
                case (1)                                                        ! Newton-Rhapson
                    corr = -fun_vals(0)/fun_vals(1)
                case (2)                                                        ! Halley
                    corr = -fun_vals(0)*fun_vals(1)/&
                        &(fun_vals(1)**2 - 0.5_dp*&
                        &fun_vals(0)*fun_vals(2))
                case (3)                                                        ! 3rd order
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
            relaxed_enough = .false.
            do id = 1,max_nr_backtracks_loc
                zero_new = zero + corr                                          ! propose a new zero
                
                fun_new = fun(zero_new,0)                                       ! calculate function value there
                
                if (abs(fun_new)-abs(fun_vals(0)).le.0._dp) then                ! error has decreased
                    relaxed_enough = .true.
                    zero = zero_new
                    exit
                end if
                corr = corr*0.5_dp
            end do
            
            ! check for relaxation
            if (.not.relaxed_enough) then
                err_msg = trim(i2str(max_nr_backtracks_HH))//' backtracks were &
                    &not sufficient to get a better guess for the zero'
                zero = 0.0_dp
            else
                ! output
                if (output_loc) call writo(trim(i2str(jd))//' / '&
                    &//trim(i2str(max_it_zero))//&
                    &' - relative error: '//&
                    &trim(r2str(abs(corr)))//' (tol: '//&
                    &trim(r2str(tol_zero))//')')
                
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
                    end if
#endif
                    err_msg = 'Not converged after '//trim(i2str(jd))//&
                        &' iterations, with residual '//&
                        &trim(r2strt(corr))//' and final value '//&
                        &trim(r2strt(zero))
                    zero = 0.0_dp
                end if
            end if
        end do HH
#if ldebug
    contains
        ! plots corrections
        !> \private
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
            call print_ex_2D(['function  ','derivative'],'',y_out,x=x_out)
            
            ! user output
            call writo('Last correction '//trim(r2strt(corrs(size(corrs))))&
                &//' with tolerance '//trim(r2strt(tol_zero)))
            
            ! plot corrections and values
            call print_ex_2D('corrs','',corrs)
            call print_ex_2D('values','',values)
        end subroutine plot_evolution
#endif
    end function calc_zero_HH_0D
    !> \private 3-D version
    !! \param[inout] fun <tt>fun(dims,x,ord)</tt> with
    !!  - <tt>dims(3)</tt> dimension of abscissa and ordinate
    !!  - <tt>x(dims(1),dims(2),dims(3))</tt> abscissa
    !!  - <tt>ord</tt> order of derivative
    !!  - <tt>fun(dims(1),dims(2),dims(3))</tt> ordinate
    function calc_zero_HH_3D(dims,zero,fun,ord,guess,max_nr_backtracks,output) &
        &result(err_msg)
        use num_vars, only: max_it_zero, tol_zero, max_nr_backtracks_HH
        
        ! input / output
        integer, intent(in) :: dims(3)                                          !< dimensions of the problem
        real(dp), intent(inout) :: zero(dims(1),dims(2),dims(3))                !< output
        interface
            function fun(dims,x,ord)                                            ! the function
                use num_vars, only: dp
                integer, intent(in) :: dims(3)
                real(dp), intent(in) :: x(dims(1),dims(2),dims(3))
                integer, intent(in) :: ord
                real(dp) :: fun(dims(1),dims(2),dims(3))
            end function fun
        end interface
        integer, intent(in) :: ord                                              !< order of solution
        real(dp), intent(in) :: guess(dims(1),dims(2),dims(3))                  !< first guess
        integer, intent(in), optional :: max_nr_backtracks                      !< max nr. of backtracks
        logical, intent(in), optional :: output                                 !< give output on convergence
        character(len=max_str_ln) :: err_msg                                    !< possible error message
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: max_nr_backtracks_loc                                        ! local max_nr_backtracks
        integer :: mc_ind(3)                                                    ! index of maximum correction
        logical :: output_loc                                                   ! local output
        logical :: relaxed_enough                                               ! relaxed enough
        real(dp) :: zero_new(dims(1),dims(2),dims(3))                           ! new zeros calculated
        real(dp) :: fun_new(dims(1),dims(2),dims(3))                            ! function values at new zero
        real(dp) :: corr(dims(1),dims(2),dims(3))                               ! correction
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
                &implemented, not '//trim(i2str(ord))
            return
        end if
        
        ! set up local output
        output_loc = .false.
        if (present(output)) output_loc = output
        
        ! set up zero
        zero = guess
        
        ! set up local max_nr_backtracks
        max_nr_backtracks_loc = max_nr_backtracks_HH
        if (present(max_nr_backtracks)) &
            &max_nr_backtracks_loc = max_nr_backtracks
        
        ! initialize function values
        allocate(fun_vals(dims(1),dims(2),dims(3),0:ord))
        
#if ldebug
        if (debug_calc_zero) then
            ! set up corrs
            allocate(corrs(dims(1),dims(2),dims(3),max_it_zero))
            allocate(values(dims(1),dims(2),dims(3),max_it_zero))
        end if
#endif
        
        ! loop over iterations
        HH: do jd = 1,max_it_zero
            ! calculate function values
            do kd = 0,ord
                fun_vals(:,:,:,kd) = fun(dims,zero,kd)
            end do
            
            ! correction to theta_HH
            select case (ord)
                case (1)                                                        ! Newton-Rhapson
                    corr = -fun_vals(:,:,:,0)/fun_vals(:,:,:,1)
                case (2)                                                        ! Halley
                    corr = -fun_vals(:,:,:,0)*fun_vals(:,:,:,1)/&
                        &(fun_vals(:,:,:,1)**2 - 0.5_dp*&
                        &fun_vals(:,:,:,0)*fun_vals(:,:,:,2))
                case (3)                                                        ! 3rd order
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
            relaxed_enough = .false.
            do id = 1,max_nr_backtracks_loc
                zero_new = zero                                                 ! otherwise some values might stay unitialized
                where (abs(corr).gt.tol_zero) zero_new = zero + corr            ! propose a new zero
                
                fun_new = fun(dims,zero_new,0)                                  ! calculate function value there
                
#if ldebug
                if (debug_calc_zero) &
                    &call writo('backtrack '//trim(i2str(id))//' of max. '//&
                    &trim(i2str(max_nr_backtracks_loc))//' - criterion: '//&
                    &trim(r2str(maxval(abs(fun_new)-abs(fun_vals(:,:,:,0)))))//&
                    &' <= 0?')
#endif
                if (maxval(abs(fun_new)-abs(fun_vals(:,:,:,0))).le.0._dp) then  ! error has decreased
                    relaxed_enough = .true.
                    zero = zero_new
                    exit
                end if
                corr = corr*0.5_dp
            end do
            
            ! check for relaxation
            if (.not.relaxed_enough) then
                err_msg = trim(i2str(max_nr_backtracks_HH))//' backtracks were &
                    &not sufficient to get a better guess for the zero'
                zero = 0.0_dp
                exit
            else
                ! output
                if (output_loc) call writo(trim(i2str(jd))//' / '&
                    &//trim(i2str(max_it_zero))//&
                    &' - maximum relative error: '//&
                    &trim(r2str(maxval(abs(corr))))//' (tol: '//&
                    &trim(r2str(tol_zero))//')')
                
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
                    end if
#endif
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
#if ldebug
    contains
        ! plots corrections
        !> \private
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
                &descr='corrections')
            
            ! plot values
            file_name = 'values'
            var_name = 'var'
            
            call plot_HDF5(var_name,file_name,values,col=2,&
                &descr='values')
        end subroutine plot_evolution
#endif
    end function calc_zero_HH_3D
    
    !> Finds the zero of a function  using Zhang's method, which is simpler than
    !! Brent's method.
    !!
    !! Taken  from  from  Steven  Stage's  correction  of  Zhang's  paper 
    !! \cite zhang2011improvement.
    !!
    !! Unlike Householder, Zhang's method needs  an interval \c x_int_in to work
    !! in,  not  a guess.  Also,  it  does not  require  the  derivative of  the
    !! function.
    !!
    !! The routine returns  an error message if  no zero is found,  and which is
    !! empty otherwise.
    !!
    !! \note The interval \c x_int_in needs to be so that the function values at
    !! either end are of different value.
    !! \param[inout] fun <tt>fun(x)</tt> with
    !!  - <tt>x</tt> abscissa
    !!  - <tt>fun</tt> ordinate
    function calc_zero_Zhang(zero,fun,x_int_in) result(err_msg)
        use num_vars, only: max_it_zero, tol_zero
        
        ! input / output
        real(dp), intent(inout) :: zero                                         !< output
        interface
            function fun(x)                                                     ! the function
                use num_vars, only: dp
                real(dp), intent(in) :: x 
                real(dp) :: fun 
            end function fun
        end interface
        real(dp), intent(in) :: x_int_in(2)                                     !< interval
        character(len=max_str_ln) :: err_msg                                    !< possible error message
        
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
                    if (.not.converged) x_rule = 0.5_dp*(x_int(id_loc)+x_mid)
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
                if (.not.converged) x_rule = x_int(id_loc) - &
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
        !> \private
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
            call print_ex_2D('function','',y_out,x=x_out)
            
            ! user output
            call writo('Final interval ['//trim(r2str(x_ints(n_ints,1)))//&
                &'..'//trim(r2str(x_ints(n_ints,2)))//'] with tolerance '//&
                &trim(r2strt(tol_zero)))
            
            ! plot intervals
            call print_ex_2D(['x_intervals'],'',x_ints)
        end subroutine plot_evolution
#endif
    end function calc_zero_Zhang
end module num_ops
