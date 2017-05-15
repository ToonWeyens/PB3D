!------------------------------------------------------------------------------!
!   Generic tests                                                              !
!------------------------------------------------------------------------------!
module test
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: pi, dp, max_str_ln, iu
    use input_utilities
    
    implicit none
    private
#if ldebug
    public generic_tests
#endif

#if ldebug
contains
    ! performs generic tests
    integer function generic_tests() result(ierr)
        use num_vars, only: prog_style
        
        character(*), parameter :: rout_name = 'generic_tests'
        
        ! initialize ierr
        ierr = 0
        
        call writo('Test interpolation?')
        if (get_log(.false.)) then
            ierr = test_interp()
            CHCKERR('')
            call pause_prog
        end if
        
        call writo('Test lock system?')
        if (get_log(.false.)) then
            ierr = test_lock()
            CHCKERR('')
            call pause_prog
        end if
        
        call writo('Test calculation of RZL?')
        if (get_log(.false.)) then
            ierr = test_calc_RZL()
            CHCKERR('')
            call pause_prog
        end if
        
        ! select according to program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                call writo('Test calculation of derivatives?')
                if (get_log(.false.)) then
                    ierr = test_calc_deriv()
                    CHCKERR('')
                    call pause_prog
                end if
            case(2)                                                             ! POST
                call writo('Test reading of subsets?')
                if (get_log(.false.)) then
                    ierr = test_read_HDF5_subset()
                    CHCKERR('')
                    call pause_prog
                end if
                
                call writo('Test calculation of volume integral?')
                if (get_log(.false.)) then
                    ierr = test_calc_int_vol()
                    CHCKERR('')
                    call pause_prog
                end if
        end select
    end function generic_tests
    
    ! test interpolation
    integer function test_interp() result(ierr)
        use MPI_utilities, only: wait_MPI
        use grid_vars, only: disc_type
        use grid_utilities, only: setup_interp_data, apply_disc
        
        character(*), parameter :: rout_name = 'test_interp'
        
        ! local variables
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_interp(:)
        real(dp), allocatable :: f(:)
        real(dp), allocatable :: f_interp(:)
        integer :: interp_style
        integer :: fun_style
        integer :: n
        integer :: n_interp
        integer :: ord
        real(dp) :: x_lim(2)
        real(dp) :: x_lim_interp(2)
        logical :: change_lims
        type(disc_type) :: interp_data
        logical :: is_trigon
        real(dp), allocatable :: x_plot(:,:)
        real(dp), allocatable :: f_plot(:,:)
        real(dp) :: gam
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
        call lvl_ud(1)
        
        ! user input
        call writo('which kind of interpolation?')
        call lvl_ud(1)
        call writo('1. polynomical interpolation')
        call writo('2. trigonometric interpolation')
        interp_style = get_int(lim_lo=1,lim_hi=2)
        call lvl_ud(-1)
        is_trigon = .false.
        select case (interp_style)
            case(1)
                x_lim = [0._dp,1._dp]
            case(2)
                x_lim = [0._dp,2*pi]
                is_trigon = .true.
        end select
        
        ! set up the function to be interpolated
        call writo('function to interpolate?')
        call lvl_ud(1)
        call writo('1: sin(x)')
        call writo('2: sin(x) + 0.25*cos(2*x)')
        call writo('3: sin(x) + 0.25*cos(2*x) + 2*sin(6*x) + 3*cos(7*x)')
        call writo('4: 0.5+x-2*x^3+x^4')
        call writo('5: x^7')
        call writo('6: x^6')
        call writo('7: x^5')
        call writo('8: x^4')
        call writo('9: square wave')
        call writo('10: sawtooth')
        call writo('11: triangle wave')
        call writo('12: sine composed by squares')
        fun_style = get_int(lim_lo=1,lim_hi=12)
        call lvl_ud(-1)
        
        call writo('order of interpolation?')
        ord = get_int(lim_lo=1)
        
        ! set up the fixed points and function
        call setup_x(x,'at which to tabulate',x_lim,ord+1)
        n = size(x)
        allocate(f(n))
        select case (fun_style)
            case (1)
                f = sin(x)
            case (2)
                f = sin(x) + 0.25*cos(2*x)
            case (3)
                f = sin(x) + 0.25*cos(2*x) + 2*sin(6*x) + 3*cos(7*x)
            case (4)
                f = 0.5 + x-2*x**3 + x**4
            case (5)
                f = x**7
            case (6)
                f = x**6
            case (7)
                f = x**5
            case (8)
                f = x**4
            case (9)
                do id = 1,n
                    if (x(id)-x_lim(1) .le. (x_lim(2)-x_lim(1))/2) then
                        f(id) = 0._dp
                    else
                        f(id) = 1._dp
                    end if
                end do
            case (10)
                f = x
            case (11)
                do id = 1,n
                    if (x(id)-x_lim(1) .le. (x_lim(2)-x_lim(1))/2) then
                        f(id) = x(id)
                    else
                        f(id) = x_lim(2)-x(id)
                    end if
                end do
            case (12)
                do id = 1,n
                    gam = (x_lim(2)-x_lim(1))/2
                    if (x(id)-x_lim(1) .le. (x_lim(2)-x_lim(1))/2) then
                        f(id) = -(x(id)-x_lim(1))/(x_lim(2)+x_lim(1)) * &
                            &(x(id)-gam)/(3*x_lim(1)+x_lim(2))*16
                    else
                        f(id) = (x(id)-x_lim(2))/(x_lim(2)-x_lim(1)) * &
                            &(x(id)-gam)/(x_lim(2)-x_lim(1))*16
                    end if
                end do
        end select
        
        ! set up interpolation points
        x_lim_interp = [x(1),x(n)]
        call setup_x(x_interp,'at which to interpolate',x_lim_interp,n_min=2)
        n_interp = size(x_interp)
        allocate(f_interp(n_interp))
        
        ! interpolate
        ierr = setup_interp_data(x,x_interp,interp_data,ord,is_trigon)
        CHCKERR('')
        ierr = apply_disc(f,interp_data,f_interp)
        CHCKERR('')
        
        ! visualize
        allocate(x_plot(max(n,n_interp),2))
        allocate(f_plot(max(n,n_interp),2))
        x_plot(:,1) = x(n)
        x_plot(:,2) = x_interp(n_interp)
        f_plot(:,1) = f(n)
        f_plot(:,2) = f_interp(n_interp)
        x_plot(1:n,1) = x
        f_plot(1:n,1) = f
        x_plot(1:n_interp,2) = x_interp
        f_plot(1:n_interp,2) = f_interp
        call print_ex_2D(['f       ','f_interp'],'',f_plot,x=x_plot)
        
        if (interp_style.eq.2) then
            call writo('Spectrum of f?')
            if (get_log(.true.)) call plot_fft(x,f)
            
            call writo('Spectrum of interpolated f?')
            call lvl_ud(1)
            call writo('Note: The end points are missing: not reliable &
                &spectrum', alert=.true.)
            call lvl_ud(-1)
            if (get_log(.false.)) call plot_fft(x_interp,f_interp)
        end if
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        call lvl_ud(-1)
        
        call lvl_ud(-1)
    contains
        subroutine setup_x(x,description,x_lim,n_min)
            ! input / output
            real(dp), intent(inout), allocatable :: x(:)
            character(len=*), intent(in) :: description
            real(dp), intent(inout) :: x_lim(2)
            integer, intent(in) :: n_min
            
            ! local variables
            integer :: n
            integer :: style
            integer :: id
            
            ! set up the number of points
            call writo('nr. of points '//description//'?')
            n = get_int(lim_lo=n_min)
            call writo('which sampling distribution?')
            call lvl_ud(1)
            call writo('1. equidistant grid')
            call writo('2. more points in center')
            call writo('3. more points at edges')
            style = get_int(lim_lo=1,lim_hi=3)
            call lvl_ud(-1)
            call writo('changed proposed limits ['//&
                &trim(r2strt(x_lim(1)))//','//trim(r2strt(x_lim(2)))//')?')
            change_lims = get_log(.false.)
            if (change_lims) then
                call writo('lower limit?')
                x_lim(1) = get_real()
                call writo('upper limit?')
                x_lim(2) = get_real(lim_lo=x_lim(1))
            end if
            
            allocate(x(n))
            x = [((id-1._dp)/n*2*pi-pi,id=1,n)]                                 ! equidistant base -pi..pi
            select case (style)
                case (1)
                    ! do nothing
                case (2)
                    x = x - 0.75*sin(x)
                case (3)
                    x = x + 0.75*sin(x)
            end select
            x = x_lim(1) + (x+pi)/(2*pi)*(x_lim(2)-x_lim(1))
        end subroutine setup_x
        
        subroutine plot_fft(x,f)
            ! input / output
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: f(:)
            
            ! local variables
            integer :: n
            integer :: m_F
            real(dp), allocatable :: f_loc(:)
            real(dp), allocatable :: work(:)
            real(dp), allocatable :: x_plot(:,:)
            real(dp), allocatable :: f_plot(:,:)
            
            ! fft
            call writo('!!! The following does not work for &
                &non-equidistant grids !!!',alert=.true.)
            n = size(x,1)
            allocate(work(3*n+15))
            allocate(f_loc(n))
            f_loc = f
            call dffti(n,work)
            call dfftf(n,f_loc,work)
            ! rescale
            f_loc = f_loc*2/n
            f_loc(1) = f_loc(1)/2
            ! separate cos and sine
            m_F = (n-1)/2
            ! plot
            allocate(x_plot(m_F,2))
            allocate(f_plot(m_F,2))
            x_plot = 0._dp
            f_plot = 0._dp
            f_plot(:,1) = f_loc(2:2*m_F:2)
            f_plot(:,2) = -f_loc(3:2*m_F+1:2)
            call print_ex_2D(['cos','sin'],'',f_plot(1:m_f,:))
        end subroutine
    end function test_interp
    
    ! test lock system
    integer function test_lock() result(ierr)
        use MPI_vars, only: lock_type
        use MPI_utilities, only: lock_req_acc, lock_return_acc, wait_MPI, &
            &lock_header, lock_wl_change, &
            &debug_lock
        use num_vars, only: rank
        
        character(*), parameter :: rout_name = 'test_lock'
        
        ! local variables
        type(lock_type) :: lock
        integer :: wait_time
        integer(kind=8) :: sleep_time(2)
        integer :: testcase
        integer :: n_cycles
        integer :: id
        logical :: blocking
        logical :: red_lat
        integer, allocatable :: wl(:)
        character(len=max_str_ln) :: title
        
        ! initialize ierr
        ierr = 0
        
        ! debug messages
        debug_lock = .true.
        
        call writo('which test?')
        call lvl_ud(1)
        call writo('1. blocking')
        call writo('2. non-blocking')
        call writo('3. blocking and non-blocking')
        testcase = get_int(lim_lo=1,lim_hi=3)
        call lvl_ud(-1)
        
        call writo('How many cycles?')
        n_cycles = get_int(lim_lo=1)
        
        if (testcase.eq.2 .or. testcase.eq.3) then
            call writo('Do you want to check whether locking and unlocking the &
                &window multiple times by the master helps reduce latency?')
            red_lat = get_log(.false.)
        end if
        
        ! set name
        select case (testcase)
            case (1)
                title = 'blocking'
            case (2)
                title = 'non-blocking'
            case (3)
                title = 'blocking and non-blocking'
        end select
        
        ! create lock
        ierr = lock%init(100)
        CHCKERR('')
        
        ! test blocking access
        call writo('Testing '//trim(title)//' access')
        call lvl_ud(1)
        call writo('Every process immediately requests access and then &
            &waits some seconds between 0 and 3')
        
        ! set name
        select case (testcase)
            case (1)
                call writo('Every process is blocking')
            case (2)
                call writo('Every process is non-blocking')
            case (3)
                call writo('In the first cycle, every process is non-blocking')
                call writo('Afterwards, each process becomes blocking &
                    &sequentially')
        end select
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        do id = 1,n_cycles
            select case (testcase)
                case (1)                                                        ! blocking
                    blocking = .true.
                case (2)                                                        ! non-blocking
                    blocking = .false.
                case (3)                                                        ! both
                    if (id.eq.rank+2) then
                        blocking = .true.
                    else
                        blocking = .false.
                    end if
            end select
            ierr = lock_req_acc(lock,blocking)
            CHCKERR('')
            
            wait_time = random_int([0,3])
            write(*,*) trim(lock_header(lock)),'sleeps ',&
                &wait_time,' seconds'
            if (red_lat .and. rank.eq.0 .and. wait_time.gt.0) then
                call system_clock(sleep_time(1))
                call system_clock(sleep_time(2))
                do while (1.E-9_dp*(sleep_time(2)-sleep_time(1)).lt.wait_time)
                    debug_lock = .false.
                    ierr = lock_wl_change(2,lock%blocking,lock,wl)
                    CHCKERR('')
                    call system_clock(sleep_time(2))
                end do
                debug_lock = .true.
            else
                call sleep(wait_time)
            endif
            write(*,*) trim(lock_header(lock)),'has slept ',&
                &wait_time,' seconds, returning lock'
            
            ierr = lock_return_acc(lock)
            CHCKERR('')
        end do
        
        ! clean up
        ierr = lock%dealloc()
        CHCKERR('')
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        call lvl_ud(-1)
    contains
        ! not very high quality random number between the limits
        integer function random_int(lim) result(res)
            ! input / output
            integer, intent(in) :: lim(2)
            
            ! local variables
            integer :: seed(12)
            real(dp) :: r_nr
            integer :: n_steps = 5
            integer :: id
            
            do id = 1,n_steps
                ! random number for seed
                call random_number(r_nr)
                
                ! set seed
                seed = int(1000_dp*r_nr*(rank+1))
                call random_seed(put=seed)
            end do
            
            ! random numbeer
            call random_number(r_nr)
            res = lim(1) + floor((lim(2)-lim(1)+1)*r_nr)
        end function random_int
    end function test_lock
    
    ! tests the calculation of derivatives
    integer function test_calc_deriv() result(ierr)
        use num_vars, only: rank
        use grid_vars, only: disc_type
        use grid_utilities, only: setup_deriv_data, apply_disc
        
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
        type(disc_type) :: deriv_data
        integer :: input_type
        integer :: prec
        real(dp), allocatable :: x(:)                                           ! abscissa
        integer :: max_deriv = 4
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Going to test the numerical calculation of derivatives')
        call lvl_ud(1)
        
        ! read user inputs
        call writo('precision?')
        prec = get_int(lim_lo=1)
        
        call writo('how many steps?')
        n_steps = get_int(lim_lo=10)
        
        call writo('starting (max) step size?')
        start_step = get_real(lim_hi=0.05_dp)
        
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
                        ierr = setup_deriv_data(step_size(id),loc_max,&
                            &deriv_data,kd,prec)
                        CHCKERR('')
                        ierr = apply_disc(varin,deriv_data,var_nm(:,kd))
                        CHCKERR('')
                    end do
                else
                    do kd = 1,max_deriv
                        ierr = setup_deriv_data(x,deriv_data,kd,prec)
                        CHCKERR('')
                        ierr = apply_disc(varin,deriv_data,var_nm(:,kd))
                        CHCKERR('')
                    end do
                end if
                call deriv_data%dealloc()
                
                ! max and mean errors
                do kd = 1,max_deriv
                    max_error(id,kd) = maxval(abs(var_an(:,kd)-var_nm(:,kd)))/ &
                    &maxval(abs(var_an(:,kd)))
                    mean_error(id,kd) = sum(abs(var_an(:,kd)-var_nm(:,kd)))/&
                        &loc_max/ maxval(abs(var_an(:,kd)))
                end do
                
                ! plots
                if (ind_plots) then
                    call print_ex_2D('input variable','',varin,x=x)
                    
                    do kd = 1,max_deriv
                        call print_ex_2D('analytical deriv. ord. '//&
                            &trim(i2str(kd)),'',var_an(:,kd),x=x)
                        call print_ex_2D('numerical deriv. ord. '//&
                            &trim(i2str(kd)),'',var_nm(:,kd),x=x)
                        call print_ex_2D('diff deriv. ord. '//&
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
                call print_ex_2D(['max. error','mean error']//&
                    &' for deriv. of order '//trim(i2str(id)),'',&
                    &plot_y,x=plot_x)
                plot_x(:,1) = log10(step_size)
                plot_x(:,2) = log10(step_size)
                plot_y(:,1) = log10(max_error(:,id))
                plot_y(:,2) = log10(mean_error(:,id))
                call print_ex_2D(['max. error','mean error']//&
                    &' for deriv. of order '//trim(i2str(id))//' [log-log]','',&
                    &plot_y,x=plot_x)
            end do
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_calc_deriv
    
    ! tests the calculation of R, Z, and L
    integer function test_calc_RZL() result(ierr)
        use PB3D_ops, only: reconstruct_PB3D_in
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'test_calc_RZL'
        
        ! initialize ierr
        ierr = 0
        
        ! preliminary things
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        
        ! select according to equilibrium style
        select case (eq_style)
            case(1)                                                             ! VMEC
                ierr = test_calc_RZL_VMEC()
                CHCKERR('')
            case(2)                                                             ! HELENA
                ! do nothing
        end select
        
        ! clean up
        call dealloc_in()
    contains
        integer function test_calc_RZL_VMEC() result(ierr)
            use num_vars, only: rank
            use grid_vars, only: n_r_eq, grid_type
            use grid_utilities, only: calc_eqd_grid
            use VMEC_utilities, only: calc_trigon_factors, fourier2real
            use VMEC_vars, only: R_V_c, R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, &
                &is_asym_V
            
            character(*), parameter :: rout_name = 'test_calc_RZL_VMEC'
            
            ! local variables
            integer :: id, jd                                                   ! counters
            integer :: dims(3)                                                  ! dimensions of grid
            integer :: norm_id                                                  ! normal point for which to do the analysis
            real(dp) :: grid_lims(2,2)                                          ! limits on grid
            real(dp), allocatable :: norm(:,:,:)                                ! copy of normal grid variable
            real(dp), allocatable :: R(:,:,:), Z(:,:,:), L(:,:,:)               ! local R, Z and L
            type(grid_type) :: grid                                             ! grid
            integer :: deriv(3)                                                 ! local derivative
            integer :: max_deriv(3)                                             ! maximum derivative
            
            ! user output
            call writo('Going to test the calculation of R, Z and L')
            call lvl_ud(1)
            
            ! initialize ierr
            ierr = 0
            
            if (rank.eq.0) then
                ! set normal point
                call writo('For which normal point do you want the analysis?')
                norm_id = get_int(lim_lo=1,lim_hi=n_r_eq)
                
                ! set dimensions [theta,zeta,r]
                dims = [1,1,1]
                do id = 2,3
                    call writo('grid size in dimension '//trim(i2str(id))//'?')
                    dims(id-1) = get_int(lim_lo=1)
                    
                    call writo('limits for angle '//trim(i2str(id))//' [pi]?')
                    grid_lims(1,id-1) = get_real()*pi
                    grid_lims(2,id-1) = get_real(lim_lo=grid_lims(1,id-1)/pi)*pi
                end do
                
                ierr = grid%init(dims)
                CHCKERR('')
                
                ! create grid
                allocate(norm(dims(1),dims(2),dims(3)))
                ierr = calc_eqd_grid(grid%theta_E,grid_lims(1,1),&
                    &grid_lims(2,1),1)
                CHCKERR('')
                ierr = calc_eqd_grid(grid%zeta_E,grid_lims(1,2),&
                    &grid_lims(2,2),2)
                CHCKERR('')
                norm = (norm_id-1._dp)/(n_r_eq-1)
                grid%r_E = norm(1,1,:)
                
                ! set maximum derivatives [r,theta,zeta]
                max_deriv(1) = 0
                do id = 2,3
                    call writo('maximum derivative in dimension '//&
                        &trim(i2str(id))//'?')
                    max_deriv(id) = get_int(lim_lo=0,lim_hi=3)
                end do
                
                ! calculate R, Z and L
                allocate(R(dims(1),dims(2),dims(3)))
                allocate(Z(dims(1),dims(2),dims(3)))
                allocate(L(dims(1),dims(2),dims(3)))
                ierr = calc_trigon_factors(grid%theta_E,grid%zeta_E,&
                    &grid%trigon_factors)
                CHCKERR('')
                do jd = 0,max_deriv(3)
                    do id = 0,max_deriv(2)
                        deriv = [0,id,jd]
                        call writo('Calculating for deriv = ['//&
                            &trim(i2str(deriv(1)))//','//&
                            &trim(i2str(deriv(2)))//','//&
                            &trim(i2str(deriv(3)))//']')
                        call lvl_ud(1)
                        ierr = fourier2real(R_V_c(:,norm_id:norm_id,deriv(1)),&
                            &R_V_s(:,norm_id:norm_id,deriv(1)),&
                            !&grid%theta_E,grid%zeta_E,&
                            &grid%trigon_factors,&
                            &R,sym=[.true.,is_asym_V],&
                            &deriv=[deriv(2),deriv(3)])
                        CHCKERR('')
                        ierr = fourier2real(Z_V_c(:,norm_id:norm_id,deriv(1)),&
                            &Z_V_s(:,norm_id:norm_id,deriv(1)),&
                            !&grid%theta_E,grid%zeta_E,&
                            &grid%trigon_factors,&
                            &Z,sym=[is_asym_V,.true.],&
                            &deriv=[deriv(2),deriv(3)])
                        CHCKERR('')
                        ierr = fourier2real(L_V_c(:,norm_id:norm_id,deriv(1)),&
                            &L_V_s(:,norm_id:norm_id,deriv(1)),&
                            !&grid%theta_E,grid%zeta_E,&
                            &grid%trigon_factors,&
                            &L,sym=[is_asym_V,.true.],&
                            &deriv=[deriv(2),deriv(3)])
                        CHCKERR('')
                        call plot_HDF5(['R','Z','L'],'TEST_RZL_0_'//&
                            &trim(i2str(id))//'_'//trim(i2str(jd)),&
                            &reshape([R,Z,L],[dims,3]),&
                            &X=reshape([grid%theta_E],[dims,1]),&
                            &Y=reshape([grid%zeta_E],[dims,1]),&
                            &Z=reshape([norm],[dims,1]))
                        call lvl_ud(-1)
                    end do
                end do
            end if
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end function test_calc_RZL_VMEC
    end function test_calc_RZL
    
    ! tests reading of HDF5 subset
    integer function test_read_HDF5_subset() result(ierr)
        use num_vars, only: rank, eq_style, n_procs
        use MPI_utilities, only: wait_MPI
        use grid_vars, only: grid_type
        use grid_utilities, only: copy_grid, calc_XYZ_grid, trim_grid
        use rich_vars, only: rich_lvl
        use input_utilities, only: dealloc_in
        use PB3D_ops, only: reconstruct_PB3D_in, reconstruct_PB3D_grid, &
            &reconstruct_PB3D_eq_1, reconstruct_PB3D_eq_2, &
            &reconstruct_PB3D_X_1, reconstruct_PB3D_sol
        use eq_vars, only: eq_1_type, eq_2_type
        use X_vars, only: X_1_type, &
            &m_X
        use sol_vars, only: sol_type
        use X_ops, only: setup_nm_X
        use VMEC_utilities, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'test_read_HDF5_subset'
        
        ! local variables
        type(grid_type) :: grid                                                 ! equilibrium or perturbation grid
        type(grid_type) :: grid_sub                                             ! grid for subset
        type(grid_type) :: grid_trim                                            ! trimmed grid for subset
        type(grid_type) :: grid_eq                                              ! equilibrium grid for XYZ calculation
        type(eq_1_type) :: eq_1                                                 ! flux equilibrium variables
        type(eq_2_type) :: eq_2                                                 ! metric equilibrium variables
        type(X_1_type) :: X_1                                                   ! vectorial perturbation variables
        type(sol_type) :: sol                                                   ! solution variables
        integer :: rich_lvl_name                                                ! either the Richardson level or zero, to append to names
        integer :: rich_lvl_name_eq                                             ! either the Richardson level or zero, to append to names
        integer :: id, jd, kd                                                   ! counters
        integer :: i_sub                                                        ! counter
        integer :: n_tot(3)                                                     ! n of equilibrium grid
        integer :: n_sub(3)                                                     ! nr. of subsets
        integer :: n_sub_norm                                                   ! nr. of subsets
        integer :: norm_ol = 4                                                  ! normal overlap
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        integer :: i_lim(2)                                                     ! normal limit for this process
        integer :: plot_dim(3)                                                  ! dimensions of plot
        integer :: plot_offset(3)                                               ! local offset of plot
        integer :: test_case                                                    ! which case to test (1: eq_2, 2: X_1)
        integer :: i_sec                                                        ! which mode for X_1 or sol
        integer, allocatable :: n_sub_loc(:)                                    ! number of subsets per process
        integer, allocatable :: n_sub_norm_loc(:)                               ! number of normal poins per process per subset
        integer, allocatable :: lims(:,:,:)                                     ! limits (3,2,product(n_sub))
        integer, allocatable :: lims_loc(:,:,:)                                 ! limits for current process
        real(dp), allocatable :: XYZ(:,:,:,:)                                   ! X, Y and Z on output grid
        real(dp), allocatable :: XYZ_i(:,:,:,:)                                 ! X, Y and Z for slab plot
        character(len=max_str_ln) :: description                                ! description for plot
        character(len=8) :: name_suff(2)                                        ! suffix of name, without and possibly with rank appended
        character(len=3) :: grid_name                                           ! either eq or X
        character :: XYZ_name(3)                                                ! X, Y and Z
        logical :: direct_grid_read                                             ! whether grid is directly read with subsets or copied
        logical :: divide_over_procs                                            ! whether a division over process should be made
        logical :: pers_out                                                     ! whether output is persistent (all procs) or not (only master)
        logical :: tot_rich                                                     ! whether variable is spread out over multiple Richardson levels or not
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Checking reading of HDF5 subsets')
        call lvl_ud(1)
        
        ! some preliminary things
        ierr = reconstruct_PB3D_in('in')                                        ! reconstruct miscellaneous PB3D output variables
        CHCKERR('')
        XYZ_name(1) = 'X'
        XYZ_name(2) = 'Y'
        XYZ_name(3) = 'Z'
        
        ! get information about which variables to reconstruct from user
        call writo('Which variable do you want to test?')
        call lvl_ud(1)
        call writo('1. eq_2: Jac_FD and kappa_n')
        call writo('2. X_1:  IM(U_0), RE(DU_1)')
        call writo('3. sol:  RE(sol_vec)')
        call lvl_ud(-1)
        test_case = get_int(lim_lo=1,lim_hi=3)
        
        ! set up  whether Richardson level  has to be  appended to the  name and
        ! whether tot_rich is used or not
        select case (eq_style) 
            case (1)                                                            ! VMEC
                select case (test_case)
                    case (1)                                                    ! eq_2
                        rich_lvl_name = rich_lvl                                ! append richardson level
                        tot_rich = .true.
                    case (2)                                                    ! X_1
                        rich_lvl_name = rich_lvl                                ! append richardson level
                        tot_rich = .true.
                    case (3)                                                    ! sol
                        rich_lvl_name = 0                                       ! do not append name
                        tot_rich = .false.
                end select
                rich_lvl_name_eq = rich_lvl                                     ! append richardson level
            case (2)                                                            ! HELENA
                rich_lvl_name = 0                                               ! do not append name
                rich_lvl_name_eq = 0                                            ! do not append name
                tot_rich = .false.
        end select
        select case (test_case)
            case (1)                                                            ! eq_2
                grid_name = 'eq'
            case (2)                                                            ! X_1
                grid_name = 'X'
            case (3)                                                            ! sol
                grid_name = 'sol'
        end select
        
        ! reconstruct full grid
        ierr = reconstruct_PB3D_grid(grid,trim(grid_name),&
            &rich_lvl=rich_lvl_name,tot_rich=.true.)
        CHCKERR('')
        n_tot = grid%n
        
        ! get direct read flag from user
        call writo('Read the grids directly using subsets?')
        direct_grid_read = get_log(.true.)
        
        ! get flag concerning division over MPI procs
        if (n_procs.gt.1) then
            call writo('Divide over MPI processes?')
            divide_over_procs = get_log(.true.)
        else
            divide_over_procs = .false.
        end if
        
        ! get limits from user
        call writo('Total '//trim(grid_name)//' grid size: ['//&
            &trim(i2str(n_tot(1)))//','//trim(i2str(n_tot(2)))//','//&
            &trim(i2str(n_tot(3)))//']')
        do id = 1,3
            if (n_tot(id).gt.0) then
                call writo('how many divisions in dimension '//trim(i2str(id)))
                n_sub(id) = get_int(lim_lo=0,lim_hi=n_tot(id)-1) + 1
            else
                n_sub(id) = 1
            end if
        end do
        allocate(lims(3,2,product(n_sub)))
        
        ! also construct equilibrium grid for XYZ calculation, using lims
        lims(1,:,1) = [0,0]
        lims(2,:,1) = [0,0]
        lims(3,:,1) = [-1,-1]
        ierr = reconstruct_PB3D_grid(grid_eq,'eq',rich_lvl=rich_lvl_name_eq,&
            &tot_rich=.true.,lim_pos=lims(:,:,1))
        CHCKERR('')
        call print_ex_2D('r_E','r_E_'//trim(i2str(rank)),grid_eq%r_E,draw=.false.)
        call draw_ex(['r_E'],'r_E_'//trim(i2str(rank)),1,1,.false.)
        
        ! set limits
        i_sub = 1
        do kd = 1,n_sub(3)
            do jd = 1,n_sub(2)
                do id = 1,n_sub(1)
                    lims(1,1,i_sub) = (id-1)*(n_tot(1)/n_sub(1)) + &
                        &min(mod(n_tot(1),n_sub(1)),id-1) + 1
                    lims(1,2,i_sub) = id*(n_tot(1)/n_sub(1)) + &
                        &min(mod(n_tot(1),n_sub(1)),id)
                    lims(2,1,i_sub) = (jd-1)*(n_tot(2)/n_sub(2)) + &
                        &min(mod(n_tot(2),n_sub(2)),jd-1) + 1
                    lims(2,2,i_sub) = jd*(n_tot(2)/n_sub(2)) + &
                        &min(mod(n_tot(2),n_sub(2)),jd)
                    lims(3,1,i_sub) = (kd-1)*(n_tot(3)/n_sub(3)) + &
                        &min(mod(n_tot(3),n_sub(3)),kd-1) + 1
                    lims(3,2,i_sub) = kd*(n_tot(3)/n_sub(3)) + &
                        &min(mod(n_tot(3),n_sub(3)),kd)
                    i_sub = i_sub+1
                end do
            end do
        end do
        
        ! divide under processes:  either just dividing over  all processes with
        ! full normal  subranges, or by letting  every process work on  the same
        ! subranges, dividing these normally.
        allocate(n_sub_loc(n_procs))
        allocate(n_sub_norm_loc(n_procs))
        if (divide_over_procs) then
            pers_out = rank.eq.0
            n_sub_loc = product(n_sub)
            allocate(lims_loc(3,2,product(n_sub)))
            lims_loc = lims
            call writo('All process perform every job together, dividing each &
                &one normally with overlap '//trim(i2str(norm_ol)))
        else
            pers_out = .true.
            n_sub_loc = product(n_sub)/n_procs
            n_sub_loc(1:mod(product(n_sub),n_procs)) = &
                &n_sub_loc(1:mod(product(n_sub),n_procs)) + 1
            allocate(lims_loc(3,2,n_sub_loc(rank+1)))
            lims_loc = &
                &lims(:,:,sum(n_sub_loc(1:rank))+1:sum(n_sub_loc(1:rank+1)))
            call writo('Process '//trim(i2str(rank))//' has subsets:',&
                &persistent=.true.)
        end if
        call lvl_ud(1)
        do i_sub = 1,size(lims_loc,3)
            call writo('subset '//trim(i2str(i_sub)),persistent=pers_out)
            call lvl_ud(1)
            do id = 1,3
                call writo('lim('//trim(i2str(id))//'&
                    &) = ['//trim(i2str(lims_loc(id,1,i_sub)))//'..'//&
                    &trim(i2str(lims_loc(id,2,i_sub)))//']',&
                    &persistent=pers_out)
            end do
            call lvl_ud(-1)
        end do
        call lvl_ud(-1)
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! for  X_1 and sol, set  up nm and  get information about which  mode to
        ! plot
        if (test_case.eq.2 .or. test_case.eq.3) then
            ! set up nm in full grids
            ierr = reconstruct_PB3D_eq_1(grid_eq,eq_1,'eq_1')
            CHCKERR('')
            ierr = setup_nm_X(grid_eq,grid,eq_1,plot_nm=.false.)                ! is necessary for X variables
            CHCKERR('')
            call eq_1%dealloc()
            ! get user input
            call writo('Which perturbation mode do you want to plot?')
            i_sec = get_int(lim_lo=1,lim_hi=size(m_X,2))
        end if
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! loop over all josb
        do i_sub = 1,size(lims_loc,3)
            ! user output
            call writo('Process '//trim(i2str(rank))//', subset '//&
                &trim(i2str(i_sub))//':',persistent=.true.)
            call lvl_ud(1)
            
            ! initialize description
            description = ''
            do id = 1,3
                description = trim(description)//' lim('//trim(i2str(id))//&
                    &') = ['//trim(i2str(lims_loc(id,1,i_sub)))//'..'//&
                    &trim(i2str(lims_loc(id,2,i_sub)))//']'
                if (id.lt.3) description = trim(description)// ','
            end do
            
            ! set up i_lim
            n_sub_norm = lims_loc(3,2,i_sub)-lims_loc(3,1,i_sub)+1
            if (divide_over_procs) then
                n_sub_norm_loc = n_sub_norm/n_procs
                n_sub_norm_loc(1:mod(n_sub_norm,n_procs)) = &
                    &n_sub_norm_loc(1:mod(n_sub_norm,n_procs)) + 1
                i_lim = [sum(n_sub_norm_loc(1:rank))+1,&
                    &sum(n_sub_norm_loc(1:rank+1))]
                if (rank+1.lt.n_procs) i_lim(2) = i_lim(2)+norm_ol
                description = trim(description)//', normal range = '//&
                    &trim(i2str(i_lim(1)))//'..'//trim(i2str(i_lim(2)))
            else
                n_sub_norm_loc = n_sub_norm
                i_lim = lims_loc(3,:,i_sub) - lims_loc(3,1,i_sub) + 1
            end if
            
            ! output description
            description = adjustl(description)
            call writo(trim(description),persistent=.true.)
            
            ! set up subset grid
            if (direct_grid_read) then
                ! check whether it coincides with a direct read using subsets
                ierr = reconstruct_PB3D_grid(grid_sub,trim(grid_name),&
                    &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                    &lim_pos=lims_loc(:,:,i_sub),grid_limits=i_lim)
                CHCKERR('')
            else
                ierr = copy_grid(grid,grid_sub,&
                    &lims_B=lims_loc(:,:,i_sub),i_lim=i_lim)
                CHCKERR('')
            end if
            
            ! trim grid
            ierr = trim_grid(grid_sub,grid_trim,norm_id)
            CHCKERR('')
            
            ! set up plot dimensions and local dimensions
            plot_dim = grid_trim%n
            plot_offset = [0,0,grid_trim%i_min-1]
            
            ! reconstruct metric equilibrium or vectorial perturbation variables
            ! of subset
            select case (test_case)
                case (1)                                                        ! eq_2
                    ierr = reconstruct_PB3D_eq_2(grid_sub,eq_2,'eq_2',&
                        &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                        &lim_pos=lims_loc(:,:,i_sub))
                    CHCKERR('')
                case (2)                                                        ! X_1
                    ierr = reconstruct_PB3D_X_1(grid_sub,X_1,'X_1',&
                        &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                        &lim_pos=lims_loc(:,:,i_sub),&
                        &lim_sec_X=[i_sec,i_sec])
                    CHCKERR('')
                case (3)                                                        ! sol
                    ierr = reconstruct_PB3D_sol(grid_sub,sol,'sol',&
                        &rich_lvl=rich_lvl,&
                        &lim_pos=lims_loc(3:3,:,i_sub),&
                        &lim_sec_sol=[i_sec,i_sec])
                    CHCKERR('')
            end select
            
            ! set up name suffix
            if (divide_over_procs) then
                name_suff(1) = '_'//trim(i2str(i_sub))
            else
                name_suff(1) = '_'//trim(i2str(sum(n_sub_loc(1:rank))+i_sub))
            end if
            if (divide_over_procs) then
                name_suff(2) = trim(name_suff(1))//'_'//trim(i2str(rank))
            else
                name_suff(2) = trim(name_suff(1))
            end if
            
            ! set up XYZ and XYZ_i
            if (test_case.eq.1 .or. test_case.eq.2) then
                if (eq_style.eq.1) then                                         ! if VMEC, calculate trigonometric factors of output grid
                    ierr = calc_trigon_factors(grid_trim%theta_E,&
                        &grid_trim%zeta_E,grid_trim%trigon_factors)
                    CHCKERR('')
                end if
                allocate(XYZ(grid_trim%n(1),grid_trim%n(2),&
                    &grid_trim%loc_n_r,3))
                ierr = calc_XYZ_grid(grid_eq,grid_trim,XYZ(:,:,:,1),&
                    &XYZ(:,:,:,2),XYZ(:,:,:,3))
                CHCKERR('')
                allocate(&
                    &XYZ_i(grid_trim%n(1),grid_trim%n(2),grid_trim%loc_n_r,3))
                do id = 1,grid_trim%n(1)
                    XYZ_i(id,:,:,1) = lims_loc(1,1,i_sub) + id - 2
                end do
                do jd = 1,grid_trim%n(2)
                    XYZ_i(:,jd,:,2) = lims_loc(2,1,i_sub) + jd - 2
                end do
                do kd = 1,grid_trim%loc_n_r
                    XYZ_i(:,:,kd,3) = lims_loc(3,1,i_sub) + kd - 2 + &
                        &grid_trim%i_min - 1
                end do
            else
                allocate(XYZ_i(1,grid_trim%loc_n_r,1,1))
                do kd = 1,grid_trim%loc_n_r
                    XYZ_i(1,kd,1,1) = lims_loc(3,1,i_sub) + kd - 2 + &
                        &grid_trim%i_min - 1
                end do
            end if
            
            ! plot
            select case (test_case)
                case (1)                                                        ! eq_2
                    call plot_HDF5('jac_FD',&
                        &'jac_FD'//trim(name_suff(1)),&
                        &eq_2%jac_FD(:,:,norm_id(1):norm_id(2),0,0,0),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &description=description)
                    call plot_HDF5('kappa_n',&
                        &'kappa_n'//trim(name_suff(1)),&
                        &eq_2%kappa_n(:,:,norm_id(1):norm_id(2)),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &description=description)
                case (2)                                                        ! X_1
                    call plot_HDF5('U_0_IM',&
                        &'U_0_IM'//trim(name_suff(1)),&
                        &ip(X_1%U_0(:,:,norm_id(1):norm_id(2),1)),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &description=description)
                    call plot_HDF5('DU_1_RE',&
                        &'DU_1_RE'//trim(name_suff(1)),&
                        &rp(X_1%DU_1(:,:,norm_id(1):norm_id(2),1)),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &description=description)
                case (3)                                                        ! sol
                    call print_ex_2D(['sol_RE'],'sol_RE'//trim(name_suff(2)),&
                        &rp(sol%vec(1,norm_id(1):norm_id(2),:)),&
                        &x=XYZ_i(1,:,:,1),draw=.false.)
                    call draw_ex(['sol_RE'],'sol_RE'//trim(name_suff(2)),&
                        &size(sol%vec,3),1,.false.)
            end select
            if (test_case.eq.1 .or. test_case.eq.2) then
                do id = 1,3
                    call plot_HDF5(XYZ_name(id),XYZ_name(id)//&
                        &trim(name_suff(2)),XYZ(:,:,:,id),&
                        &X=XYZ_i(:,:,:,1),Y=XYZ_i(:,:,:,2),Z=XYZ_i(:,:,:,3))
                end do
            end if
            
            ! clean up
            call grid_sub%dealloc()
            call grid_trim%dealloc()
            select case (test_case)
                case (1)                                                        ! eq_2
                    call eq_2%dealloc()
                case (2)                                                        ! X_1
                    call X_1%dealloc()
                case (3)                                                        ! sol
                    call sol%dealloc()
            end select
            if (test_case.eq.1 .or. test_case.eq.2) deallocate(XYZ)
            deallocate(XYZ_i)
            call lvl_ud(-1)
        end do
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! Now for entire region if wanted
        call writo('Do you want to plot the entire region at once?')
        if (get_log(.false.) .and. rank.eq.0) then
            ! reconstruct metric equilibrium or vectorial perturbation variables
            select case (test_case)
                case (1)                                                        ! eq_2
                    ierr = reconstruct_PB3D_eq_2(grid,eq_2,'eq_2',&
                        &rich_lvl=rich_lvl_name,tot_rich=.true.)
                    CHCKERR('')
                case (2)                                                        ! X_1
                    ierr = reconstruct_PB3D_X_1(grid,X_1,'X_1',&
                        &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                        &lim_sec_X=[i_sec,i_sec])
                    CHCKERR('')
                case (3)                                                        ! sol
                    ierr = reconstruct_PB3D_sol(grid,sol,'sol',&
                        &rich_lvl=rich_lvl,&
                        &lim_sec_sol=[i_sec,i_sec])
                    CHCKERR('')
            end select
            
            ! set up XYZ and XYZ_i
            if (test_case.eq.1 .or. test_case.eq.2) then
                if (eq_style.eq.1) then                                         ! if VMEC, calculate trigonometric factors of output grid
                    ierr = calc_trigon_factors(grid%theta_E,&
                        &grid%zeta_E,grid%trigon_factors)
                    CHCKERR('')
                end if
                allocate(XYZ(grid%n(1),grid%n(2),grid%loc_n_r,3))
                ierr = calc_XYZ_grid(grid_eq,grid,XYZ(:,:,:,1),&
                    &XYZ(:,:,:,2),XYZ(:,:,:,3))
                CHCKERR('')
                allocate(XYZ_i(grid%n(1),grid%n(2),grid%loc_n_r,3))
                do id = 1,grid%n(1)
                    XYZ_i(id,:,:,1) = id - 1
                end do
                do jd = 1,grid%n(2)
                    XYZ_i(:,jd,:,2) = jd - 1
                end do
                do kd = 1,grid%loc_n_r
                    XYZ_i(:,:,kd,3) = kd - 1
                end do
            else
                allocate(XYZ_i(1,grid%loc_n_r,1,1))
                do kd = 1,grid%loc_n_r
                    XYZ_i(1,kd,1,1) = kd - 1
                end do
            end if
            
            ! plot
            description = 'total range'
            select case (test_case)
                case (1)                                                        ! eq_2
                    call plot_HDF5('jac_FD','jac_FD_tot',&
                        &eq_2%jac_FD(:,:,:,0,0,0),&
                        &description=description,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3))
                    call plot_HDF5('kappa_n','kappa_n_tot',eq_2%kappa_n,&
                        &description=description,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3))
                case (2)                                                        ! X_1
                    call plot_HDF5('U_0_IM','U_0_IM_tot',&
                        &ip(X_1%U_0(:,:,:,1)),&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &description=description)
                    call plot_HDF5('DU_1_RE','DU_1_RE_tot',&
                        &rp(X_1%DU_1(:,:,:,1)),&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &description=description)
                case (3)                                                        ! sol
                    call print_ex_2D(['sol_RE'],'sol_RE_tot',&
                        &rp(sol%vec(1,:,:)),x=XYZ_i(1,:,:,1),draw=.false.)
                    call draw_ex(['sol_RE'],'sol_RE_tot',&
                        &size(sol%vec,3),1,.false.)
            end select
            if (test_case.eq.1 .or. test_case.eq.2) then
                do id = 1,3
                    call plot_HDF5(XYZ_name(id),XYZ_name(id)//'_tot',&
                        &XYZ(:,:,:,id),&
                        &X=XYZ_i(:,:,:,1),Y=XYZ_i(:,:,:,2),Z=XYZ_i(:,:,:,3))
                end do
            end if
            
            ! clean up
            call eq_2%dealloc()
            if (test_case.eq.1 .or. test_case.eq.2) deallocate(XYZ)
            deallocate(XYZ_i)
            call lvl_ud(-1)
        end if
        
        ! clean up
        call grid%dealloc()
        if (test_case.eq.2 .or. test_case.eq.3) call grid_eq%dealloc()
        call dealloc_in()
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Test complete')
    end function test_read_HDF5_subset
    
    ! tests calculation of volume integral
    integer function test_calc_int_vol() result(ierr)
        use num_vars, only: rank
        use grid_utilities, only: calc_eqd_grid, calc_int_vol, &
            &debug_calc_int_vol
        
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
                call plot_HDF5(var_name(2),file_name(2),rp(fun(:,:,:,1)),&
                    &X=X,Y=Y,Z=Z)
                call plot_HDF5(var_name(3),file_name(3),ip(fun(:,:,:,1)),&
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
                call print_ex_2D(['analytical integral, RE',&
                    &'numerical integral, RE '],'',rp(fun_int_plot))
                call print_ex_2D(['analytical integral, IM',&
                    &'numerical integral, IM '],'',ip(fun_int_plot))
            end if
        end if 
    end function test_calc_int_vol
#endif
end module test
