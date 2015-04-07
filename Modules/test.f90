!------------------------------------------------------------------------------!
!   Test some routines and functions                                           !
!------------------------------------------------------------------------------!
module test
#include <PB3D_macros.h>
    use num_vars, only: dp, max_str_ln, pi, mu_0, output_i, iu
    use output_ops, only: print_GP_2D, print_GP_3D
    use messages, only: print_ar_1, print_ar_2, start_time, &
        &stop_time
    use input_ops, only: yes_no
    use str_ops, only: i2str, r2str, r2strt

    implicit none
    private
    public test_calc_deriv, test_conv_FHM
    
contains
    integer function test_calc_deriv() result(ierr)
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
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'test calc_deriv?'
        if(yes_no(.false.)) then
            output_i = 0
            
            ! read user input
            do
                write(*,*) 'how many steps ?'
                read(*,*) n_steps
                if (n_steps.lt.10) then
                    write(*,*) 'choose a value larger than 10'
                    cycle
                else
                    exit
                end if
            end do
            do
                write(*,*) 'precision ?'
                    write(*,*) '1: truncation error ~ Delta^2'
                    write(*,*) '2: truncation error ~ Delta^3'
                read(*,*) prec
                if (prec.lt.1 .or. prec.gt.2) then
                    write(*,*) 'choose a value from the correct range'
                    cycle
                else
                    exit
                end if
            end do
            do
                write(*,*) 'starting (max) step size ?'
                read(*,*) start_step
                if (start_step.gt.1E-1_dp) then
                    write(*,*) 'Choose a value smaller than 0.1'
                    cycle
                else
                    exit
                end if
            end do
            write(*,*) 'equidistant plot? [yes]'
            if(yes_no(.true.)) then
                eqd_plots = .true.
            else
                eqd_plots = .false.
            end if
            if (.not.eqd_plots) then
                do
                    write(*,*) 'which type of x-axis?'
                        write(*,*) '1: linear x = (x-1)/(N-1)'
                        write(*,*) '2: quadratic x = ((x-1)/(N-1))^2'
                    read(*,*) grid_type
                    if (grid_type.lt.1 .or. grid_type.gt.2) then
                        write(*,*) 'choose a value from the range 1-2'
                        cycle
                    else
                        exit
                    end if
                end do
            end if
            
            write(*,*) 'Individual plots?'
            if(yes_no(.false.)) then
                ind_plots = .true.
            else
                ind_plots = .false.
            end if
            
            do
                write(*,*) 'which function?'
                    write(*,*) '1: sin(pi*x) + 0.25*cos(2*pi*x)'
                    write(*,*) '2: 0.5+x-2*x^3+x^4'
                    write(*,*) '3: x^7'
                    write(*,*) '4: x^6'
                    write(*,*) '5: x^5'
                    write(*,*) '6: x^4'
                read(*,*) input_type
                if (input_type.lt.1 .or. input_type.gt.6) then
                    write(*,*) 'Choose a value from the range 1-6'
                    cycle
                else
                    exit
                end if
            end do
            
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
                    write(*,*) 'ERROR: in test_calc_deriv, how did &
                        &you get here?!??!'
            end select
            
            ! initialize
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
                    write(*,*) 'setting up equidistant grid'
                    x = [((kd-1._dp)*step_size(id),kd=1,loc_max)]
                else
                    write(*,*) 'setting up regular grid'
                    allocate(x(loc_max))
                    select case (grid_type)
                        case(1)
                            x = [((kd-1._dp)*step_size(id),kd=1,loc_max)]
                        case(2)
                            x = [(((kd-1._dp)*step_size(id))**2,kd=1,loc_max)]
                        case default
                            write(*,*) 'ERROR: you cannot have regular a grid &
                                &of type '//trim(i2str(grid_type))
                            exit
                    end select
                end if
                
                write(*,*) trim(i2str(id))//'/'//trim(i2str(n_steps))// &
                    &': calculating using '//trim(i2str(loc_max))//' points'
                
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
                        var_an(:,1) =  (pi)**1*cos(pi*x)-0.25*(2*pi)**1*sin(2*pi*x)
                        var_an(:,2) = -(pi)**2*sin(pi*x)-0.25*(2*pi)**2*cos(2*pi*x)
                        var_an(:,3) = -(pi)**3*cos(pi*x)+0.25*(2*pi)**3*sin(2*pi*x)
                        var_an(:,4) =  (pi)**4*sin(pi*x)+0.25*(2*pi)**4*cos(2*pi*x)
                        var_an(:,5) =  (pi)**5*cos(pi*x)-0.25*(2*pi)**5*sin(2*pi*x)
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
                        write(*,*) 'ERROR: in test_calc_deriv, how did &
                            &you get here?!??!'
                        stop
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
            
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_calc_deriv
    
    integer function test_conv_FHM() result(ierr)
        use utilities, only: conv_FHM
        
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
        
        write(*,*) 'test conv_FHM?'
        if(yes_no(.false.)) then
            output_i = 0
            write(*,*) 'Testing whether h2f*f2h = 1'
            write(*,*) 'The relative difference between an original FM variable&
                & var and h2f*f2h*var, should be decreasing with increasing &
                &number of radial points'
            write(*,*) ''
            write(*,*) 'Do you want the individual plots?'
            ind_plot = yes_no(.false.)
            write(*,*) ''
            do
                write(*,*) 'n_max = ?'
                read(*,*) n_max
                if (n_max.lt.4*length) then
                    write(*,*) 'n_max has to be larger than or equal to '&
                        &//trim(i2str(4*length))//'...'
                else 
                    exit
                end if
            end do
            write(*,*) 'logarithmic plot?'
            log_plot = yes_no(.false.)
            write(*,*) ''
            
            
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
                
                if(ind_plot) then
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
                    write(*,*) 'max = '//trim(r2strt(maxerr(id)))&
                        &//', average: '//trim(r2strt(averr(id)))
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
                
            write(*,*) 'Paused... press enter'
            read(*,*)
            
            write(*,*) 'Stopping'
            stop
        end if
    end function test_conv_FHM
end module test
