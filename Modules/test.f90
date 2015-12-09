!------------------------------------------------------------------------------!
!   Generic tests                                                              !
!------------------------------------------------------------------------------!
module test
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln, iu
    
    implicit none
    private
#if ldebug
    public generic_tests
#endif

#if ldebug
contains
    ! performs generic tests
    integer function generic_tests() result(ierr)
        use input_ops, only: pause_prog, get_log
        use num_vars, only: prog_style
        
        character(*), parameter :: rout_name = 'generic_tests'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! select according to program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                call writo('Test calculation of derivatives?')
                if (get_log(.false.)) then
                    ierr = test_calc_deriv()
                    CHCKERR('')
                    call pause_prog
                end if
                
                call writo('Test conversion between full and half mesh?')
                if (get_log(.false.)) then
                    ierr = test_conv_FHM()
                    CHCKERR('')
                    call pause_prog
                end if
            case(2)                                                             ! POST
                call writo('Test calculation of volume integral?')
                if (get_log(.false.)) then
                    ierr = test_calc_int_vol()
                    CHCKERR('')
                    call pause_prog
                end if
            case default
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function generic_tests
    
    ! tests the calculation of derivatives
    integer function test_calc_deriv() result(ierr)
        use num_vars, only: rank
        use input_ops, only: get_log, get_int, get_real
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'test_calc_deriv'
        
        ! local variables
        integer :: loc_max, id, kd
        integer :: n_steps, grid_type
        real(dp) :: start_step
        real(dp), allocatable :: varin(:)
        real(dp), allocatable :: var_an(:,:), var_nm(:,:)
        logical :: ind_plots, eqd_plots
        real(dp), allocatable :: step_size(:)
        real(dp), allocatable :: max_error(:,:)
        real(dp), allocatable :: mean_error(:,:)
        real(dp), allocatable :: plot_y(:,:)
        real(dp), allocatable :: plot_x(:,:)
        integer :: input_type
        integer :: prec
        real(dp), allocatable :: x(:)                                           ! abscissa
        integer :: max_deriv
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Going to test the numerical calculation of derivatives')
        call lvl_ud(1)
        
        ! read user inputs
        call writo('how many steps?')
        n_steps = get_int(lim_lo=10)
        
        call writo('precision?')
        call lvl_ud(1)
        call writo('1: truncation error ~ Delta^2')
        call writo('2: truncation error ~ Delta^3')
        call lvl_ud(-1)
        prec = get_int(lim_lo=1,lim_hi=2)
        
        call writo('starting (max) step size?')
        start_step = get_real(lim_hi=0.1_dp)
        
        call writo('equidistant plot?')
        eqd_plots = get_log(.true.)
        
        if (.not.eqd_plots) then
            call writo('which type of x-axis?')
            call lvl_ud(1)
            call writo('1: linear x = (x-1)/(N-1)')
            call writo('2: quadratic x = ((x-1)/(N-1))^2')
            call lvl_ud(-1)
            grid_type = get_int(lim_lo=1,lim_hi=2)
        end if
        
        call writo('individual plots?')
        ind_plots = get_log(.false.)
        
        call writo('which function?')
        call lvl_ud(1)
        call writo('1: sin(pi*x) + 0.25*cos(2*pi*x)')
        call writo('2: 0.5+x-2*x^3+x^4')
        call writo('3: x^7')
        call writo('4: x^6')
        call writo('5: x^5')
        call writo('6: x^4')
        call lvl_ud(-1)
        input_type = get_int(lim_lo=1,lim_hi=6)
        
        ! only do the tests by the master
        if (rank.eq.0) then
            ! set up max_deriv
            select case (prec)
                case(1)
                    if (eqd_plots) then
                        max_deriv = 5
                    else
                        max_deriv = 3
                    end if
                case(2)
                    if (eqd_plots) then
                        max_deriv = 1
                    else
                        max_deriv = 1
                    end if
                case default
                    ierr = 1
                    err_msg = 'precision has to be 1 or 2'
                    CHCKERR(err_msg)
            end select
            
            ! initialize variables
            allocate(step_size(n_steps)); step_size = 0.0_dp
            allocate(max_error(n_steps,max_deriv)); max_error = 0.0_dp
            allocate(mean_error(n_steps,max_deriv)); mean_error = 0.0_dp
            allocate(plot_x(n_steps,2)); plot_x = 0.0_dp
            allocate(plot_y(n_steps,2)); plot_y = 0.0_dp
            
            ! n_steps calculations up to min_step
            do id = 1,n_steps
                loc_max = nint(id/start_step)                                   ! how many points for this step id
                
                if (allocated(x)) deallocate(x)
                step_size(id) = 1._dp/(loc_max-1._dp)                           ! equidistant step size for loc_max points
                
                if (eqd_plots) then
                    call writo('setting up equidistant grid')
                    x = [((kd-1._dp)*step_size(id),kd=1,loc_max)]
                else
                    call writo('setting up regular grid')
                    allocate(x(loc_max))
                    select case (grid_type)
                        case(1)
                            x = [((kd-1._dp)*step_size(id),kd=1,loc_max)]
                        case(2)
                            x = [(((kd-1._dp)*step_size(id))**2,kd=1,loc_max)]
                        case default
                            ierr = 1
                            err_msg = 'grid_type has to be 1 or 2'
                            CHCKERR(err_msg)
                    end select
                end if
                
                call writo(trim(i2str(id))//'/'//trim(i2str(n_steps))// &
                    &': calculating using '//trim(i2str(loc_max))//' points')
                
                ! set up input function
                if (allocated(varin)) deallocate(varin)
                if (allocated(var_an)) deallocate(var_an)
                if (allocated(var_nm)) deallocate(var_nm)
                
                allocate(varin(loc_max))
                allocate(var_an(loc_max,5))
                allocate(var_nm(loc_max,5))
                
                select case (input_type)
                    case (1)
                        varin = sin(pi*x)+0.25*cos(2*pi*x)
                        var_an(:,1) = &
                            &(pi)**1*cos(pi*x)-0.25*(2*pi)**1*sin(2*pi*x)
                        var_an(:,2) = &
                            &-(pi)**2*sin(pi*x)-0.25*(2*pi)**2*cos(2*pi*x)
                        var_an(:,3) = &
                            &-(pi)**3*cos(pi*x)+0.25*(2*pi)**3*sin(2*pi*x)
                        var_an(:,4) = &
                            &(pi)**4*sin(pi*x)+0.25*(2*pi)**4*cos(2*pi*x)
                        var_an(:,5) = &
                            &(pi)**5*cos(pi*x)-0.25*(2*pi)**5*sin(2*pi*x)
                    case (2)
                        varin = 0.5+x-2*x**3+x**4
                        var_an(:,1) = 1-6*x**2+4*x**3
                        var_an(:,2) = -12*x+12*x**2
                        var_an(:,3) = -12+24*x
                        var_an(:,4) = 24._dp
                        var_an(:,5) = 0.0_dp
                    case(3)
                        varin = x**7
                        var_an(:,1) = 7*x**6
                        var_an(:,2) = 6*7*x**5
                        var_an(:,3) = 5*6*7*x**4
                        var_an(:,4) = 4*5*6*7*x**3
                        var_an(:,5) = 3*4*5*6*7*x**2
                    case(4)
                        varin = x**6
                        var_an(:,1) = 6*x**5
                        var_an(:,2) = 5*6*x**4
                        var_an(:,3) = 4*5*6*x**3
                        var_an(:,4) = 3*4*5*6*x**2
                        var_an(:,5) = 2*3*4*5*6*x**1
                    case(5)
                        varin = x**5
                        var_an(:,1) = 5*x**4
                        var_an(:,2) = 4*5*x**3
                        var_an(:,3) = 3*4*5*x**2
                        var_an(:,4) = 2*3*4*5*x**1
                        var_an(:,5) = 1*2*3*4*5._dp
                    case(6)
                        varin = x**4
                        var_an(:,1) = 4*x**3
                        var_an(:,2) = 3*4*x**2
                        var_an(:,3) = 2*3*4*x**1
                        var_an(:,4) = 1*2*3*4._dp
                        var_an(:,5) = 0._dp
                    case default
                        ierr = 1
                        err_msg = 'Input type has to be between 1 and 6'
                        CHCKERR(err_msg)
                end select
                
                ! numerical derivatives
                if (eqd_plots) then
                    do kd = 1,max_deriv
                        ierr = calc_deriv(varin,var_nm(:,kd),&
                            &1.0_dp/step_size(id),kd,prec)
                        CHCKERR('')
                    end do
                else
                    do kd = 1,max_deriv
                        ierr = calc_deriv(varin,var_nm(:,kd),x,kd,prec)
                        CHCKERR('')
                    end do
                end if
                
                ! max and mean errors
                do kd = 1,max_deriv
                    max_error(id,kd) = maxval(abs(var_an(:,kd)-var_nm(:,kd)))/ &
                    &maxval(abs(var_an(:,kd)))
                    mean_error(id,kd) = sum(abs(var_an(:,kd)-var_nm(:,kd)))/&
                        &loc_max/ maxval(abs(var_an(:,kd)))
                end do
                
                ! plots
                if (ind_plots) then
                    call print_GP_2D('input variable','',varin,x=x)
                    
                    do kd = 1,max_deriv
                        call print_GP_2D('analytical deriv. ord. '//&
                            &trim(i2str(kd)),'',var_an(:,kd),x=x)
                        call print_GP_2D('numerical deriv. ord. '//&
                            &trim(i2str(kd)),'',var_nm(:,kd),x=x)
                        call print_GP_2D('diff deriv. ord. '//&
                            &trim(i2str(kd)),'',(var_an(:,kd)-var_nm(:,kd))&
                            &/maxval(abs(var_an(:,kd))),x=x)
                    end do
                end if
            end do
            
            ! plot errors
            do id = 1,max_deriv
                plot_x(:,1) = step_size
                plot_x(:,2) = step_size
                plot_y(:,1) = max_error(:,id)
                plot_y(:,2) = mean_error(:,id)
                call print_GP_2D('max. and mean error for deriv. of order '&
                    &//trim(i2str(id)),'',plot_y,x=plot_x)
                plot_x(:,1) = log10(step_size)
                plot_x(:,2) = log10(step_size)
                plot_y(:,1) = log10(max_error(:,id))
                plot_y(:,2) = log10(mean_error(:,id))
                call print_GP_2D('max. and mean error [log-log], for deriv. of &
                    &order '//trim(i2str(id)),'',plot_y,x=plot_x)
            end do
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_calc_deriv
    
    ! tests the conversion between full and half mesh
    integer function test_conv_FHM() result(ierr)
        use num_vars, only: rank
        use utilities, only: conv_FHM
        use input_ops, only: get_log, get_int
        
        character(*), parameter :: rout_name = 'test_conv_FHM'
        
        ! local variables
        integer :: id, kd
        integer :: n_max, step, n_r
        integer, parameter :: length = 50
        real(dp), allocatable :: varin(:), varout(:), varoutout(:), vardiff(:)
        real(dp) :: maxerr(length), averr(length), num_points(length)
        real(dp) :: plotvar(2,length)
        logical :: ind_plot, log_plot
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Going to test the numerical conversion between full and &
            &half mesh variables')
        call writo('The relative difference between an original FM variable&
            & var and h2f*f2h*var, should be decreasing with increasing &
            &number of radial points')
        call lvl_ud(1)
        
        ! get user inputs
        call writo('Do you want the individual plots?')
        ind_plot = get_log(.false.)
        
        call writo('n_max = ?')
        n_max = get_int(lim_lo=4*length)
        
        call writo('logarithmic plot?')
        log_plot = get_log(.false.)
        
        ! only do the tests by the master
        if (rank.eq.0) then
            ! set step
            step = n_max/length
            
            do id = 1,length
                n_r = id*step
                
                num_points(id) = n_r
                
                if (allocated(varin)) deallocate(varin)
                if (allocated(varout)) deallocate(varout)
                if (allocated(varoutout)) deallocate(varoutout)
                if (allocated(vardiff)) deallocate(vardiff)
                allocate(varin(n_r))
                allocate(varout(n_r))
                allocate(varoutout(n_r))
                allocate(vardiff(n_r))
                
                do kd = 1,n_r
                    ! some continuous curve:
                    varin(kd) = sin(float(kd)/n_r) + 0.5*cos(3*float(kd)/n_r)
                end do
                
                ierr = conv_FHM(varin,varout,.true.)
                CHCKERR('')
                ierr = conv_FHM(varout,varoutout,.false.)
                CHCKERR('')
                do kd = 1,n_r
                    vardiff(kd) = 2*(varin(kd)-varoutout(kd))/&
                        &(varin(kd)+varoutout(kd))
                end do
                maxerr(id) = maxval(abs(vardiff))
                averr(id) = sum(abs(vardiff))/size(vardiff)
                
                if (ind_plot) then
                    call print_GP_2D('input with '//trim(i2str(n_r))&
                        &//' radial points ('//trim(i2str(id))//'/'//&
                        &trim(i2str(length))//')','',varin)
                    call print_GP_2D('first output with '//trim(i2str(n_r))&
                        &//' radial points ('//trim(i2str(id))//'/'//&
                        &trim(i2str(length))//')','',varout)
                    call print_GP_2D('output with '//trim(i2str(n_r))&
                        &//' radial points ('//trim(i2str(id))//'/'//&
                        &trim(i2str(length))//')','',varoutout)
                    call print_GP_2D('rel. difference with '//trim(i2str(n_r))&
                        &//' radial points ('//trim(i2str(id))//'/'//&
                        &trim(i2str(length))//')','',vardiff)
                    call writo('max = '//trim(r2strt(maxerr(id)))&
                        &//', average: '//trim(r2strt(averr(id))))
                end if
            end do
            
            plotvar(1,:) = num_points
            if (log_plot) then
                plotvar(2,:) = log10(maxerr)
            else
                plotvar(2,:) = maxerr
            end if
            call print_GP_2D('maximum error as a function of numer of points',&
                &'',plotvar(2,:),x=plotvar(1,:))
            if (log_plot) then
                plotvar(2,:) = log10(averr)
            else
                plotvar(2,:) = averr
            end if
            call print_GP_2D('average error as a function of numer of points',&
                &'',plotvar(2,:),x=plotvar(1,:))
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_conv_FHM
    
    ! tests calculation of volume integral
    integer function test_calc_int_vol() result(ierr)
        use num_vars, only: rank
        use input_ops, only: get_int, get_log
        use grid_vars, only: create_grid
        use grid_ops, only: calc_eqd_grid, calc_int_vol, debug_calc_int_vol
        
        character(*), parameter :: rout_name = 'test_calc_int_vol'
        
        ! local variables
        real(dp), allocatable :: theta(:,:,:), zeta(:,:,:)                      ! angles
        complex(dp), allocatable :: fun(:,:,:,:)                                ! function to integrate
        real(dp), allocatable :: J(:,:,:)                                       ! Jacobian
        real(dp), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:)                   ! X,Y and Z for plot
        real(dp), allocatable :: norm(:)                                        ! normal variable
        real(dp), allocatable :: r(:)                                           ! radius
        real(dp) :: R_0                                                         ! geometric axis of torus
        complex(dp), allocatable :: fun_int(:,:,:)                              ! integrated function
        integer, allocatable :: dims(:,:)                                       ! dimensions of grid
        integer :: n_steps                                                      ! nr. of steps
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: eqd_grid_dum(:)                                ! equidistant grid dummy
        character(len=max_str_ln) :: var_name(3), file_name(3)                  ! name of variable and of file
        complex(dp), allocatable :: fun_int_plot(:,:)                           ! results for plotting
        
        ! initialize ierr
        ierr = 0
        
        if (rank.eq.0) then
            call writo('Checking integral of')
            call lvl_ud(1)
            call writo('f = 1 - r^2 + i cos(theta)')
            call lvl_ud(-1)
            call writo('on a torus')
            
            call writo('Debug mode?')
            if (get_log(.false.)) then
                debug_calc_int_vol = .true.
            else
                debug_calc_int_vol = .false.
            end if
            
            call writo('How many different grid sizes?')
            n_steps = get_int(lim_lo=1)
            
            allocate(dims(3,n_steps))
            allocate(eqd_grid_dum(n_steps))
            
            ! get dimensions
            if (n_steps.gt.1) then
                do kd = 1,size(dims,1)
                    call writo('Starting dimension '//trim(i2str(kd)))
                    dims(kd,1) = get_int(lim_lo=4)
                    call writo('Final dimension '//trim(i2str(kd)))
                    dims(kd,n_steps) = get_int(lim_lo=dims(kd,1)+n_steps)
                    ierr = calc_eqd_grid(eqd_grid_dum,1._dp*dims(kd,1),&
                        &1._dp*dims(kd,n_steps))
                    CHCKERR('')
                    dims(kd,:) = nint(eqd_grid_dum)
                end do
                call writo('Performing calculations for '//&
                    &trim(i2str(n_steps))//' grid sizes')
            else
                do kd = 1,size(dims)
                    call writo('Dimension '//trim(i2str(kd)))
                    dims(kd,1) = get_int(lim_lo=1)
                    if (dims(kd,1).eq.1) then
                        call lvl_ud(1)
                        call writo('Make sure the integrand does not depend on &
                            &dimension '//trim(i2str(kd)))
                        call writo('A factor such as 2pi can be missing from &
                            &resulting integral')
                        call lvl_ud(-1)
                    end if
                end do
            end if
            
            allocate(fun_int(2,1,n_steps))
            
            do jd = 1,n_steps
                ! user output
                call writo('grid size = ['//trim(i2str(dims(1,jd)))//' ,'//&
                    &trim(i2str(dims(2,jd)))//', '//trim(i2str(dims(3,jd)))//&
                    &']')
                call lvl_ud(1)
                
                ! set up variables
                allocate(theta(dims(1,jd),dims(2,jd),dims(3,jd)))
                allocate(zeta(dims(1,jd),dims(2,jd),dims(3,jd)))
                allocate(J(dims(1,jd),dims(2,jd),dims(3,jd)))
                allocate(fun(dims(1,jd),dims(2,jd),dims(3,jd),1))
                allocate(X(dims(1,jd),dims(2,jd),dims(3,jd)))
                allocate(Y(dims(1,jd),dims(2,jd),dims(3,jd)))
                allocate(Z(dims(1,jd),dims(2,jd),dims(3,jd)))
                allocate(norm(dims(3,jd)))
                allocate(r(dims(3,jd)))
                
                R_0 = 2._dp
                ierr = calc_eqd_grid(theta,0*pi,2*pi,1)
                CHCKERR('')
                ierr = calc_eqd_grid(zeta,0*pi,2*pi,2)
                CHCKERR('')
                ierr = calc_eqd_grid(norm,0._dp,1._dp)
                CHCKERR('')
                do kd = 1,dims(3,jd)
                    r(kd) = (kd-1._dp)/(dims(3,jd)-1)
                    fun(:,:,kd,1) = 1._dp - r(kd)**2 + iu*cos(theta(:,:,kd))
                    do id = 1,dims(1,jd)
                        J(id,:,kd) = r(kd)*(R_0+r(kd)*cos(theta(id,:,kd)))
                    end do
                    X(:,:,kd) = cos(zeta(:,:,kd))*(R_0+r(kd)*cos(theta(:,:,kd)))
                    Y(:,:,kd) = sin(zeta(:,:,kd))*(R_0+r(kd)*cos(theta(:,:,kd)))
                    Z(:,:,kd) = r(kd)*sin(theta(:,:,kd))
                end do
                
                var_name(1) = 'TEST_J_torus'
                var_name(2) = 'TEST_RE_fun_torus'
                var_name(3) = 'TEST_IM_fun_torus'
                file_name(1) = 'TEST_J_torus'
                file_name(2) = 'TEST_RE_fun_torus'
                file_name(3) = 'TEST_IM_fun_torus'
                if (n_steps.gt.1) then
                    do kd = 1,3
                        file_name(kd) = trim(file_name(kd))//'_'//trim(i2str(jd))
                    end do
                end if
                call plot_HDF5(var_name(1),file_name(1),J,X=X,Y=Y,Z=Z)
                call plot_HDF5(var_name(2),file_name(2),realpart(fun(:,:,:,1)),&
                    &X=X,Y=Y,Z=Z)
                call plot_HDF5(var_name(3),file_name(3),imagpart(fun(:,:,:,1)),&
                    &X=X,Y=Y,Z=Z)
                
                ! user output
                call writo('Analytically calculating integral')
                
                ! calculate integral analytically
                fun_int(1,1,jd) = R_0*pi**2 + iu*2*pi**2/3
                
                ! user output
                call writo('Numerically calculating integral')
                
                ! calculate integral numerically
                ierr = calc_int_vol(theta,zeta,norm,J,fun,fun_int(2,:,jd))
                CHCKERR('')
                
                ! output
                call writo('Analytical value: '//&
                    &trim(c2str(fun_int(1,1,jd))))
                call writo('Numerical value: '//&
                    &trim(c2str(fun_int(2,1,jd))))
                
                ! clean up
                deallocate(theta,zeta,J,fun)
                deallocate(X,Y,Z)
                deallocate(norm,r)
                
                call lvl_ud(-1)
            end do
            
            if (n_steps.gt.1) then
                allocate(fun_int_plot(n_steps,2))
                do jd = 1,n_steps
                    fun_int_plot(jd,:) = fun_int(:,1,jd)/fun_int(1,1,jd)
                end do
                call print_GP_2D('analytical vs. numerical integral, real','',&
                    &realpart(fun_int_plot))
                call print_GP_2D('analytical vs. numerical integral, imag.','',&
                    &imagpart(fun_int_plot))
            end if
        end if 
    end function test_calc_int_vol
#endif
end module test
