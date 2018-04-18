!------------------------------------------------------------------------------!
!> Generic tests.
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
    !> Performs generic tests.
    !!
    !! These tests  are run in interactive  fashion by running the  program with
    !! the <tt>--test</tt> command-line option.
    !! \see init_files()
    !!
    !! \note
    !!  -#    It   is    probably   also    recommended   to    run   it    with
    !!  <tt>--do_execute_command_line</tt>  in order  to generate  the plots  in
    !!  real time.
    !!  -#  It has  not  been  tested whether  testing  routines "pollutes"  the
    !!  further simulations.  It is  therefore best  to not  use the  tests when
    !!  doing a real run.
    !!
    !! \ldebug
    !!
    !! \return ierr
    integer function generic_tests() result(ierr)
        use num_vars, only: prog_style
        
        character(*), parameter :: rout_name = 'generic_tests'
        
        ! initialize ierr
        ierr = 0
        
        call writo('Test splines?')
        if (get_log(.false.)) then
            ierr = test_splines()
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
        
        call writo('Test calculation of toroidal functions?')
        if (get_log(.false.)) then
            ierr = test_tor_fun()
            CHCKERR('')
            call pause_prog
        end if
        
        ! select according to program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                ! do nothing
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
    
    !> Test splines.
    !!
    !! \return ierr
    integer function test_splines() result(ierr)
        use num_utilities, only: spline
        use num_vars, only: tol_zero
        
        character(*), parameter :: rout_name = 'test_splines'
        
        ! local variables
        integer :: kd, id
        integer :: nx
        integer :: nx_int
        integer :: ord
        integer :: max_deriv
        integer :: bcs(2)                                                       ! boundary conditions
        real(dp) :: bcs_val(2)                                                  ! value of derivative at boundary conditions
        real(dp) :: extrap_fraction
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: x_int(:)
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: y_int(:)
        logical :: ready
        logical :: extrap
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=5) :: side_str(2)
        
        ! initialize ierr
        ierr = 0
        
        ! user input
        call writo('testing splines')
        call lvl_ud(1)
        
        ! set up function
        nx = 10
        allocate(x(nx))
        allocate(y(nx,0:3))
        do kd = 1,nx
            x(kd) = (kd-1._dp)/(nx-1._dp)                                       ! 0..1
            y(kd,0) = sin(2*pi*x(kd))
            y(kd,1) = 2*pi*cos(2*pi*x(kd))
            y(kd,2) = -(2*pi)**2*sin(2*pi*x(kd))
            y(kd,3) = -(2*pi)**3*cos(2*pi*x(kd))
        end do
        
        ! set up interpolated function
        nx_int = 400
        allocate(x_int(nx_int))
        allocate(y_int(nx_int))
        
        ready = .false.
        do while (.not.ready)
            ! get order
            call writo('order?')
            call lvl_ud(1)
            call writo('1: linear')
            call writo('2: akima hermite')
            call writo('3: cubic')
            call lvl_ud(-1)
            ord = get_int(lim_lo=1,lim_hi=3)
            
            ! get extrapolation fraction
            call writo('test extrapolation?')
            extrap = get_log(.false.)
            if (extrap) then
                call writo('Extrapolation fraction?')
                extrap_fraction = get_real(lim_lo=0._dp,lim_hi=10._dp)
            else
                extrap_fraction = -tol_zero
            end if
            
            side_str(1) = 'left'
            side_str(2) = 'right'
            if (ord.gt.1) then
                do kd = 1,2
                    call writo('boundary condition '//trim(side_str(kd))//&
                        &' = ?')
                    call lvl_ud(1)
                    call writo('-1: periodic')
                    call writo(' 0: not-a-knot')
                    call writo(' 1: first derivative prescribed')
                    call writo(' 2: second derivative prescribed')
                    call lvl_ud(-1)
                    bcs(kd) = get_int(lim_lo=-1,lim_hi=2)
                    if (bcs(kd).eq.1 .or. bcs(kd).eq.2) then
                        call writo('value of derivative of order '//&
                            &trim(i2str(bcs(kd))))
                        bcs_val(kd) = get_real()
                    else
                        bcs_val(kd) = 0._dp
                    end if
                end do
            else
                bcs = 0
                bcs_val = 0._dp
            end if
            
            select case (ord)
                case (1:2)
                    max_deriv = 1
                case (3)
                    max_deriv = 3
            end select
            
            ! set up x_int
            do kd = 1,nx_int
                x_int(kd) = -extrap_fraction + &
                    &(1+2*extrap_fraction)*(kd-1._dp)/(nx_int-1._dp)            ! 0..1 with possible extrapolation fraction
            end do
            
            ! calculate spline and possibly plot
            do id = 0,max_deriv
                ierr = spline(x,y(:,0),x_int,y_int,ord=ord,deriv=id,bcs=bcs,&
                    &bcs_val=bcs_val,extrap=extrap)                             ! use y(0)
                CHCKERR('')
                call plot_y(x,x_int,y(:,id),y_int,id)                           ! use y(id)
            end do
            
            call writo('repeat?')
            ready = .not.get_log(.true.)
        end do
        
        call lvl_ud(-1)
    contains
        subroutine plot_y(x,x_int,y,y_int,deriv)
            ! input / output
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: x_int(:)
            real(dp), intent(in) :: y(:)
            real(dp), intent(in) :: y_int(:)
            integer, intent(in) :: deriv
            
            ! local variables
            integer :: nx, nx_int
            real(dp) :: minmax_y(2)
            real(dp), allocatable :: x_plot(:,:)
            real(dp), allocatable :: y_plot(:,:)
            character(len=max_str_ln), allocatable :: plot_names(:)
            
            ! set up plot variables
            nx = size(x)
            nx_int = size(x_int)
            allocate(x_plot(max(nx,nx_int),2))                                  ! (x_int,x)
            allocate(y_plot(max(nx,nx_int),2))                                  ! (y_int,y)
            
            ! x,y
            x_plot(1:nx_int,1) = x_int
            x_plot(nx_int+1:max(nx,nx_int),1) = x_int(nx_int)
            y_plot(1:nx_int,1) = y_int
            y_plot(nx_int+1:max(nx,nx_int),1) = y_int(nx_int)
            
            ! int x,y
            x_plot(1:nx,2) = x
            x_plot(nx+1:max(nx,nx_int),2) = x(nx)
            y_plot(1:nx,2) = y
            y_plot(nx+1:max(nx,nx_int),2) = y(nx)
            
            ! knots
            minmax_y(1) = minval(y_plot(:,1:2))
            minmax_y(2) = maxval(y_plot(:,1:2))
            
            ! plot
            allocate(plot_names(2))
            plot_names(1) = 'D'//trim(i2str(deriv))//'y_int'
            plot_names(2) = 'D'//trim(i2str(deriv))//'y'
            call print_ex_2D(plot_names,plot_names(1),y_plot,x=x_plot)
        end subroutine plot_y
    end function test_splines
    
    !> Test lock system.
    !!
    !! \return ierr
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
        !> \private
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
    
    !> Tests the calculation of R, Z, and L.
    !!
    !! \return ierr
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
        !> \private
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
    
    !> Test calculation of toroidal functions.
    !!
    !! \return ierr
    integer function test_tor_fun() result(ierr)
        use num_vars, only: rank
        use grid_utilities, only: calc_eqd_grid
        use dtorh, only: dtorh1
        
        character(*), parameter :: rout_name = 'test_tor_fun'
        
        ! local variables
        integer :: jd, kd                                                       ! counters
        integer :: npt                                                          ! number of points
        real(dp) :: lim(2)                                                      ! limits
        real(dp), allocatable :: x(:)                                           ! abscissa
        real(dp), allocatable :: f(:,:,:)                                       ! ordinate P and Q
        integer :: n                                                            ! degree
        integer :: m                                                            ! order
        integer :: max_n                                                        ! max. degree that can be calculated
        integer :: min_n                                                        ! max. degree that can be plot
        logical :: done                                                         ! whether done
        logical :: test_approx                                                  ! test approximation
        logical :: test_rel_diff                                                ! relative difference, not absolute
        character(len=max_str_ln), allocatable :: plot_name(:)                  ! names of plot
        character(len=max_str_ln), allocatable :: fun_names(:,:)                ! names of functions that are plot
        
        ! initialize ierr
        ierr = 0
        
        done = .false.
        test_approx = .false.
        
        do while (.not.done .and. rank.eq.0)
            call lvl_ud(1)
            
            ! user input
            call writo('Test approximation at singular point?')
            test_approx = get_log(.false.)
            test_rel_diff = .false.
            if (test_approx) then
                call writo('Relative difference in stead of absolute?')
                test_rel_diff = get_log(.false.)
            end if
            
            ! user output
            call writo('Lower bound to plot?')
            lim(1) = get_real(lim_lo=1._dp)
            call writo('Upper bound to plot?')
            lim(2) = get_real(lim_lo=lim(1))
            call writo('How many points?')
            npt = get_int(lim_lo=1)
            call writo('Max. degree? (subscript n in P_(n-0.5)^m)')
            n = get_int(lim_lo=0)
            call writo('Min. degree? (subscript n in P_(n-0.5)^m)')
            min_n = get_int(lim_lo=0,lim_hi=n)
            if (.not.test_approx) then
                call writo('Order? (superscript m in P_(n-0.5)^m)')
                m = get_int(lim_lo=0)
            else
                call writo('Order is 0 when testing approximation')
                m = 0
            end if
            
            ! set up
            allocate(x(npt))
            allocate(f(npt,0:n,2))
            ierr = calc_eqd_grid(x,lim(1),lim(2))
            CHCKERR('')
            if (test_approx) then
                allocate(plot_name(min_n:n))
                allocate(fun_names(min_n:n,2))
                if (test_rel_diff) then
                    do kd = min_n,n
                        plot_name(kd) = 'tor_Q_rel_diff_'//trim(i2str(kd))
                        fun_names(kd,1) = 'ReldiffQ_{'//trim(i2str(kd))//&
                            &'-0.5}^'//trim(i2str(m))
                        fun_names(kd,2) = 'Q_{'//trim(i2str(kd))//&
                            &'-0.5}^'//trim(i2str(m))
                    end do
                else
                    do kd = min_n,n
                        plot_name(kd) = 'tor_Q_comp_'//trim(i2str(kd))
                        fun_names(kd,1) = 'ApproxQ_{'//trim(i2str(kd))//&
                            &'-0.5}^'//trim(i2str(m))
                        fun_names(kd,2) = 'Q_{'//trim(i2str(kd))//&
                            &'-0.5}^'//trim(i2str(m))
                    end do
                end if
            else
                allocate(plot_name(2))
                allocate(fun_names(2,min_n:n))
                plot_name(1) = 'tor_P'
                plot_name(2) = 'tor_Q'
                do kd = min_n,n
                    fun_names(1,kd) = 'P_{'//trim(i2str(kd))//'-0.5}^'//&
                        &trim(i2str(m))
                    fun_names(2,kd) = 'Q_{'//trim(i2str(kd))//'-0.5}^'//&
                        &trim(i2str(m))
                end do
            end if
            
            ! calculate
            do kd = 1,npt
                ierr = dtorh1(x(kd),m,n,f(kd,:,1),f(kd,:,2),max_n)
                CHCKERR('')
                if (max_n.lt.n) call writo('For point x = '//&
                    &trim(r2str(x(kd)))//', the solution did not converge for &
                    &degree > '//trim(i2str(max_n)),warning=.true.)
                if (test_approx) then
                    f(kd,:,1) = -0.5*log((x(kd)-1._dp)/32)
                    do jd = 1,n
                        f(kd,jd:n,1) = f(kd,jd:n,1) - 2_dp/(2._dp*jd-1._dp)
                    end do
                    if (test_rel_diff) f(kd,:,1) = &
                        &2._dp*(f(kd,:,1)-f(kd,:,2))/(f(kd,:,1)+f(kd,:,2))
                end if
            end do
            
            ! plot
            do kd = lbound(plot_name,1),ubound(plot_name,1)
                if (test_approx) then
                    call print_ex_2D(fun_names(kd,:),trim(plot_name(kd)),&
                        &f(:,kd,:),x=reshape(x,[npt,1]))
                else
                    call print_ex_2D(fun_names(kd,:),trim(plot_name(kd)),&
                        &f(:,min_n:n,kd),x=reshape(x,[npt,1]))
                end if
                call draw_ex(fun_names(kd,:),trim(plot_name(kd)),&
                    &size(fun_names,2),1,.false.)
            end do
            
            ! clean up
            deallocate(x,f,fun_names,plot_name)
            
            ! user input
            call writo('Test again?')
            done = .not.get_log(.true.)
            
            call lvl_ud(-1)
        end do
    end function test_tor_fun
    
    !> Tests reading of HDF5 subset.
    !!
    !! \return ierr
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
        use X_vars, only: X_1_type, modes_type
        use sol_vars, only: sol_type
        use X_ops, only: init_modes, setup_modes
        use VMEC_utilities, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'test_read_HDF5_subset'
        
        ! local variables
        type(modes_type) :: mds                                                 ! general modes variables
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
        call print_ex_2D('r_E','r_E_'//trim(i2str(rank)),grid_eq%r_E,&
            &draw=.false.)
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
            ierr = init_modes(grid_eq,eq_1)
            CHCKERR('')
            ierr = setup_modes(mds,grid_eq,grid,plot_nm=.false.)
            CHCKERR('')
            call eq_1%dealloc()
            ! get user input
            call writo('Which perturbation mode do you want to plot?')
            i_sec = get_int(lim_lo=1,lim_hi=size(mds%m,2))
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
                    ierr = reconstruct_PB3D_X_1(mds,grid_sub,X_1,'X_1',&
                        &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                        &lim_pos=lims_loc(:,:,i_sub),&
                        &lim_sec_X=[i_sec,i_sec])
                    CHCKERR('')
                case (3)                                                        ! sol
                    ierr = reconstruct_PB3D_sol(mds,grid_sub,sol,'sol',&
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
                        &descr=description)
                    call plot_HDF5('kappa_n',&
                        &'kappa_n'//trim(name_suff(1)),&
                        &eq_2%kappa_n(:,:,norm_id(1):norm_id(2)),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &descr=description)
                case (2)                                                        ! X_1
                    call plot_HDF5('U_0_IM',&
                        &'U_0_IM'//trim(name_suff(1)),&
                        &ip(X_1%U_0(:,:,norm_id(1):norm_id(2),1)),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &descr=description)
                    call plot_HDF5('DU_1_RE',&
                        &'DU_1_RE'//trim(name_suff(1)),&
                        &rp(X_1%DU_1(:,:,norm_id(1):norm_id(2),1)),&
                        &tot_dim=plot_dim,loc_offset=plot_offset,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &descr=description)
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
                    ierr = reconstruct_PB3D_X_1(mds,grid,X_1,'X_1',&
                        &rich_lvl=rich_lvl_name,tot_rich=.true.,&
                        &lim_sec_X=[i_sec,i_sec])
                    CHCKERR('')
                case (3)                                                        ! sol
                    ierr = reconstruct_PB3D_sol(mds,grid,sol,'sol',&
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
                        &descr=description,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3))
                    call plot_HDF5('kappa_n','kappa_n_tot',eq_2%kappa_n,&
                        &descr=description,&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3))
                case (2)                                                        ! X_1
                    call plot_HDF5('U_0_IM','U_0_IM_tot',&
                        &ip(X_1%U_0(:,:,:,1)),&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &descr=description)
                    call plot_HDF5('DU_1_RE','DU_1_RE_tot',&
                        &rp(X_1%DU_1(:,:,:,1)),&
                        &X=XYZ(:,:,:,1),Y=XYZ(:,:,:,2),Z=XYZ(:,:,:,3),&
                        &descr=description)
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
    
    !> Tests calculation of volume integral.
    !!
    !! \return ierr
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
